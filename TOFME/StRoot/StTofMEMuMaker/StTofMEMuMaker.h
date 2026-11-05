#ifndef StTofMEMuMaker_h
#define StTofMEMuMaker_h

// StTofMEMuMaker
//   TOFME on MuDst: the conventional mixed-event estimate of random TPC-track <-> BTOF-hit matches, using
//   all BTOF hits with their own time and ToT (StMuBTofHit), which the picoDst lacks.
//   * Geometry, projection, candidate cells (projected cell, +-1 in the same module) and the mixing classes
//     (run x 14 Vz x 9 raw-refMult bins, FIFO depth nMix) are inherited from StTofMEMaker.
//   * Hit time: t = t_LE - off_cell - T0, off_cell = per-run cell offset (median of t_LE - tof - T0 of matched
//     hits, from macros/probeMuBTof.C pass 1), T0 = event median of (t_LE - off - tof) over matched hits.
//     Cells without an offset use the module (else tray) average offset.
//   * SE: production beta of the primary track (hM2SE). The production choice is emulated on this event's hits
//     (largest ToT among the candidate cells) and compared with the hit production used (hChoice).
//   * ME overlay: the partner's hits are added to this event's hits in the candidate cells; production rule
//     (largest ToT wins) -> hM2MEtot if a partner hit wins; earliest time wins -> hM2MEtime; partner hits
//     alone (no overlay, largest ToT) -> hM2MEplain. m2 with this track's p and path length and the partner
//     hit's time. Normalisation per projected track: ME * NPrj / NTr per (cent9, p) bin.

#include "StTofMEMaker/StTofMEMaker.h"

class StMuDstMaker;
class TH1D;
class TH2D;
class TH3F;

class StTofMEMuMaker : public StTofMEMaker {
 public:
  StTofMEMuMaker(const char* name, StMuDstMaker* mu, const char* outName);
  virtual ~StTofMEMuMaker();
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
  void setOffsetPattern(const char* p) { mOffPattern = p; }   // e.g. "output/mudst/off_%d.root" (hRes1 per run)

 private:
  Bool_t loadOffsets(Int_t run);
  StMuDstMaker* mMu;
  TString mOffPattern;
  Bool_t mHaveOff;

  TH1D *hEvtFlow, *hT0nMatch, *hCalRes, *hNHitEv, *hTotAll;
  TH2D *hNPrj, *hNTr;              // (cent9, p)
  TH3F *hM2SE;                     // production beta
  TH3F *hM2SEemu;                  // emulated production choice on this event's hits (largest ToT)
  TH3F *hM2MEtot, *hM2MEtime, *hM2MEplain;
  TH3F *hM2MEmerge;                // overlay, largest ToT wins, but A and B hits in the same cell merge into one hit
                                   // (earlier leading edge, owner = earlier hit, ToT = larger of the two)
  TH3F *hMEmergeSt;                // (cent9, p, 0 A wins / 1 B wins / 2 no hit)
  // ---- v2: full emulation of StBTofMatchMaker on the (overlaid) hit map: all production-valid global tracks are
  // projected; a hit cell claimed (cell within +-1 of a projection) by more than one track/hit entry is discarded;
  // among the remaining candidates of a track: one -> taken, several -> largest ToT below 40 ns.
  TH3F *hM2SEprodEmu, *hM2MEprod, *hMEprodSt;
  TH2D *hChoiceP;                  // as hChoice, for the full emulation
  TH1D *hNGlobProj;                // production-valid global tracks projected per event
  // ---- v3: the partner's own track projections also claim the partner's hits (a partner hit made by a
  // reconstructed particle is claimed by its own track as in the real event)
  TH3F *hM2MEprodB, *hMEprodBSt;
  TH3F *hQM2SE, *hQM2MEprodB;   // charge-signed q*m2
  TH3F *hM2SEneg, *hM2MEprodBneg; // m2 of negative tracks only
  TH3F *hMEtotSt;                  // (cent9, p, 0 this event wins / 1 partner wins / 2 no hit)
  TH2D *hChoice;                   // (p, code) 0 same hit as production, 1 other cell, 2 we pick none, 3 prod none
  TH2D *hDtSEprod, *hDtMEtot;      // (p, t - t_pi)
  TH1D *hCent9;

  ClassDef(StTofMEMuMaker, 4)
};
#endif
