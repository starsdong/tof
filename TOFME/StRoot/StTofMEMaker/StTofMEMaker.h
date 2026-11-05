#ifndef StTofMEMaker_h
#define StTofMEMaker_h

// StTofMEMaker
//   TOFME step 2: TPC-track -> BTOF cell projection from the picoDst helix and
//   same-event (SE) / mixed-event (ME) cell-level matching.
//
//   * Geometry: data-driven per-module map (centre C, local axes U,V, cell step W) fitted from
//     production-matched tracks: hitPos = C + yLocal*U + zLocal*V + (cell-3.5)*W.
//     The maker always accumulates the normal-equation sums (hGeomSum); macros/makeGeom.C
//     solves them and writes the geometry text file read here (setGeomFile).
//   * Projection: global helix (StPicoTrack::helix(B)) -> cylinder at <R_module> ->
//     candidate modules -> plane intersection -> (module, cell, yLocal, zLocal).
//   * SE match: a BTOF hit in the projected cell of the same event.
//     ME match: a BTOF hit in the projected cell of another event of the same
//     run / Vz class / refMult class (FIFO buffer of nMix events).
//   * m2 (v2): every BTOF hit gets the time of flight of the picoDst BTofPidTraits matched to it
//     (hits not matched to any picoDst track have no time). For the hit picked for a track
//     (projected cell, else the nearer +-1 neighbour; timed hits first) beta = L/(c*tof) with
//     L = tofPathLength(primary vertex, projected point, curvature) as in StBTofMatchMaker,
//     m2 = p^2 (1/beta^2 - 1). Filled vs (centrality9, p, m2) for SE, ME and production beta.
//     Centrality: StRefMultCorr (Au+Au 19.6 GeV 2019), bad runs and pile-up excluded from m2 histos.
//   * Overlay ME (v3): the partner's hits are added on top of this event's hits; the production-like
//     choice (projected cell, then nearer neighbour; earlier time if the same cell) is redone and,
//     when a partner hit wins, its m2 is filled (hM2MEo). This includes the shadowing by the track's
//     own hit, i.e. it estimates the random-match contribution to the SE m2 without free parameters.
//   * Background diagnostics (v4): SE and overlay-ME m2 (0-80%) split by kNK category variables of the
//     track, event and chosen hit (see kCatName in the .cxx), with the counts needed to normalise ME per
//     category; time residual t - t_pi for SE and ME.

#include "StMaker.h"
#include "TString.h"
#include <vector>

class StPicoDstMaker;
class StPicoEvent;
class StPicoPhysicalHelix;
class TFile;
class TH1D;
class TH2D;
class TH3F;
class TH2F;
class StRefMultCorr;
class TVector3;

class StTofMEMaker : public StMaker {
 public:
  StTofMEMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);
  virtual ~StTofMEMaker();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();

  void setGeomFile(const char* f)   { mGeomFile = f; }
  void setTimePdfFile(const char* f) { mPdfFile = f; }   // v10: MuDst unmatched-hit time pdf (hTU3)
  void setNMix(Int_t n)             { mNMix = n; }
  void setVzCut(Float_t v)          { mVzCut = v; }
  void setVrCut(Float_t v)          { mVrCut = v; }
  void setZLocalCut(Float_t v)      { mZLocalCut = v; }
  void addTriggerId(UInt_t id)      { mTrigIds.push_back(id); }

  enum { kNMod = 3840, kNCell = 23040, kNPar = 12 };
  enum { kNK = 12 };   // number of v4 background category variables

 protected:   // protected so that StTofMEMuMaker (MuDst version) can reuse geometry, projection and mixing classes
  Bool_t passEvent(StPicoEvent* ev);
  Bool_t loadGeom();
  // returns module index (0..3839) or -1; status: 0 ok, 1 no cylinder crossing,
  // 2 no candidate module, 3 outside cells (y), 4 outside |zLocal| window
  Int_t  project(const StPicoPhysicalHelix& h, Int_t& cell, Double_t& yl, Double_t& zl,
                 TVector3& pos, Int_t& status);
  Int_t  multClass(Int_t refMult) const;
  // pick the hit for a projected track in a cell->tof map (0 no hit, <0 hit without time, >0 tof):
  // returns 0 no hit, 1 timed hit (tof, gcellOut set), 2 hit(s) but none with time
  Int_t  candOrder(Int_t gcell, Int_t cell, Double_t yl, Double_t pitch, Bool_t nb, Int_t* cand) const;
  // overlay of partner hits on this event; seSel = cell picked in this event (timed) or -1.
  // returns 0 this-event hit kept, 1 partner timed hit wins (tof), 2 partner hit wins but has no time,
  // 3 same cell in both with unknown time order, 4 no hit at all
  Int_t  overlay(const std::vector<float>& se, const std::vector<float>& me, const Int_t* cand, Int_t n,
                 Int_t seSel, Float_t& tof, Int_t& cellOut) const;
  Int_t  pickHit(const std::vector<float>& m, Int_t gcell, Int_t cell, Double_t yl, Double_t pitch,
                 Bool_t neighbours, Float_t& tof, Int_t& gcellOut) const;

  StPicoDstMaker* mPicoDstMaker;
  TString  mOutName, mGeomFile, mPdfFile;
  TFile*   mFile;
  Bool_t   mHaveGeom;
  Int_t    mNMix;
  Float_t  mVzCut, mVrCut, mZLocalCut;
  std::vector<UInt_t> mTrigIds;

  // geometry (transient)
  Double_t mGeo[kNMod][kNPar];  //! C(3) U(3) V(3) W(3)
  Int_t    mGeoFlag[kNMod];     //! 0 none, 1 fitted, 2 rotated copy
  Double_t mModPhi[kNMod];      //!
  Double_t mModZ[kNMod];        //!
  Double_t mTrayPhi[120];       //!
  Double_t mRmean;              //!
  Int_t    mCurrentRun;         //!
  StRefMultCorr* mRefMultCorr;  //!
  Bool_t   mBadRun;             //!

  // ---- histograms
  TH1D* hEvtCut;
  TH1D* hNMixed;          // number of ME partners per event
  TH2D* hGeomSum;         // normal-equation sums per module (hadd-able)
  TH1D* hProjStatus;      // projection status for QA tracks
  TH1D* hProjModFlag;     // geometry flag of projected module
  // production comparison (tracks with prod match flag>0)
  TH1D* hProdCmp;         // 0 prod matched, 1 my proj valid, 2 same module, 3 same cell
  TH2D* hDyVsP;           // my yLocal - prod yLocal (same cell)
  TH2D* hDzVsP;           // my zLocal - prod zLocal
  TH2D* hDposVsP;         // |my point - prod hitPos|
  TH2D* hProdCmpVsP;      // p vs category
  // SE / ME counts, vs pT, eta, phi, refMult, nHits
  TH1D *hPrjPt, *hSEPt, *hSENbPt, *hMETrPt, *hMEPt, *hMENbPt, *hProdPt;
  TH1D *hPrjEta, *hSEEta, *hMETrEta, *hMEEta, *hProdEta;
  TH1D *hPrjPhi, *hSEPhi, *hMETrPhi, *hMEPhi;
  TH1D *hPrjRef, *hSERef, *hMETrRef, *hMERef;
  TH1D *hSETrNh, *hSENh, *hMETrNh, *hMENh;     // vs number of BTOF hits (of the event providing hits)
  TH1D *hPrjYl, *hSEYl, *hMEYl, *hPrjZl, *hSEZl, *hMEZl;  // projected local coords
  TH1D *hProdYl;                               // production yLocal (for comparison)
  TH2D *hPrjTrayMod, *hSEtrayMod;
  // ---- v2: centrality and m2 (cent9 x p x m2); cent9 from StRefMultCorr: 0 = 70-80% ... 8 = 0-5%
  TH1D *hCentFlag, *hCent9, *hRefMultCorrH, *hTimedFrac, *hHitTof;
  TH2D *hNPrjM2;              // projected tracks used for m2 (cent9, p)
  TH2D *hNTrM2;               // ME trials used for m2 (partner with TOF time)
  TH3F *hSEStat, *hMEStat;    // cent9, p, status (0 none, 1 timed hit, 2 hit without time), cell+-1 rule
  TH3F *hM2SE, *hM2ME;        // cell+-1 rule (production-like)
  TH3F *hM2SEc, *hM2MEc;      // projected cell only
  TH3F *hM2Prod;              // production beta, projected QA tracks
  TH2D *hDLVsP, *hDm2VsP;     // own-hit validation vs production (L and m2)
  TH3F *hMEOvStat;            // overlay status (0..4, see overlay())
  TH3F *hM2MEo;               // overlay ME: m2 of partner hits that win the choice
  // ---- v4 background diagnostics (0-80%); index k = category variable
  TH3F *hSEcat[kNK], *hMEcat[kNK];   // (category, p, m2)
  TH2D *hNPcat[kNK], *hNTcat[kNK];   // projected tracks / trials per (category, p) (track & event categories)
  TH3F *hOvStcat[kNK];               // overlay status per (category, p, status)
  TH2F *hDtSE, *hDtSEown, *hDtME;    // (p, t - t_pi)
  // ---- v5 per projected cell (0.6 < p < 3, 0-80%): HI = 1.6 < m2 < 2.8, PR = 0.6 < m2 < 1.2
  TH1D *hCellPrj, *hCellOcc, *hCellSEnb, *hCellSEnbHI, *hCellSEnbPR, *hCellSEprHI, *hCellSEpr;
  TH1D *hCellTr, *hCellMEnb, *hCellMEnbHI, *hCellMEprHI, *hCellMEnbU;
  TH2F *hDtSEnb, *hDtMEnb;           // (p, t - t_pi) for neighbour-cell choices
  // ---- v6: neighbour-cell matches in detail (0-80%)
  TH3F *hNbYproj;     // SE: (|yLocal| of my projection, p, m2), split below by chosen cell
  TH3F *hPrYproj;
  TH3F *hNbYprod;     // SE own neighbour hit: (production yLocal signed towards projected cell, p, m2)
  TH3F *hNbState;     // SE neighbour choice: (projected-cell state 0 empty / 1 hit without time, p, m2)
  TH3F *hMeNbYproj, *hMePrYproj;   // ME overlay: (|yLocal| proj, p, m2)
  // ---- v7: overlay ME with the production choice among several candidate hits (largest ToT, StBTofMatchMaker).
  // picoDst has no ToT: ToT is assumed independent of the hit origin, so each of the nA + nB hits in the
  // candidate cells (projected, +-1) of this event (A) and the partner (B) wins with probability 1/(nA+nB).
  TH3F *hM2MEtot;    // (cent9, p, m2), weight 1/(nA+nB) per timed partner hit
  TH3F *hMEtotSt;    // (cent9, p, 0 timed / 1 untimed) summed weights of partner hits
  TH3F *hSEnA;       // (cent9, p, nA) number of hits of this event in the candidate cells
  // ---- v8: overlay ME with "earliest time wins" among the candidate cells (projected, +-1): the hit chosen
  // for the track in this event (A, timed) competes with the timed partner hits (B) in the candidate cells.
  TH3F *hM2MEearly;  // (cent9, p, m2) partner hit wins
  TH3F *hMEearlySt;  // (cent9, p, status) 0 A wins, 1 B timed wins, 2 B has only untimed hits (no decision), 3 no B hit
  TH3F *hDtSEc;      // (cent9, p, t - t_pi) same event
  // ---- v9: hits without a picoDst time get a time sampled from timed hits of earlier events in the same
  // (Vz bin, BTOF half, module) class (module position + Vz ~ path length). Reservoir of recent times per class.
  TH3F *hM2MEsamp;   // (cent9, p, m2) earliest time wins among ALL hits of A and B in the candidate cells, B wins
  TH3F *hMEsampSt;   // (cent9, p, status) 0 A wins, 1 B wins, 2 only unsampled hits, 3 no B hit
  TH3F *hM2MEclos;   // closure: as hM2MEearly (v8) but B's timed hits get resampled times
  TH3F *hTofRes;     // (Vz bin, half*32+module, tof) timed hits feeding the reservoirs
  TH1D *hSampFlag;   // untimed hits: 0 sampled, 1 reservoir too small
  TH2F *hNTrSamp;    // (cent9, p) ME trials where both events are fully sampled (normalisation of v9 histograms)
  // ---- v10: untimed hits get a time sampled from the MuDst time pdf of unmatched hits (Vz bin, half*32+module);
  // timed hits keep their picoDst time; earliest time wins among all hits of A and B in the candidate cells.
  TH3F *hM2MEmu;     // (cent9, p, m2) B wins
  TH3F *hMEmuSt;     // (cent9, p, status) 0 A wins, 1 B wins, 2 -, 3 no B hit
  TH1D *hMuSampT;    // sampled times (QA)

  ClassDef(StTofMEMaker, 10)
};

#endif
