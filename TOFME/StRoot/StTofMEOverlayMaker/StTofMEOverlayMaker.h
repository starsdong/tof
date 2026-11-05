#ifndef StTofMEOverlayMaker_h
#define StTofMEOverlayMaker_h
// StTofMEOverlayMaker: mixed-event estimate of random TPC-track <-> BTOF-hit matches with the PRODUCTION
// matching and calibration re-run on MuDst (StBTofMatchMaker, StVpdCalibMaker, StBTofCalibMaker in MuDst mode).
//   Chain: StMuDstMaker, St_db_Maker, btofMatch, vpdCalib, btofCalib (MuDst mode), this maker.
//   * SE: after the chain re-run on event A, every QA primary track whose global track has a match gets
//     m2 = p^2 (1/beta^2 - 1), beta = L / (c tof), tof = calibrated time of flight of the global BTofPidTraits,
//     L = path length from the primary vertex to the matched hit position (as with mUseEventVertex).
//   * ME: for each partner B (same run, Vz bin, refMult class; FIFO nMix) B's BTOF hits are added to A's hit
//     array with leading/trailing edges shifted by tStart_A - tStart_B; a B hit in a cell that already has an A hit
//     overlapping in time (|dLE| < ToT of the earlier hit) is merged into it (earlier LE, later TE).
//     The three production makers are re-run on the overlaid event; a QA track matched to an added hit, or to a
//     merged hit whose leading edge came from B, is a random match; its m2 is filled as for SE.
//     A's hit array is restored afterwards.
//   v2 modes (setMode): 0 overlay of partner hits only (v1); 1 overlay of partner hits AND partner global tracks
//   (nHitsFit>=15), so that partner hits of reconstructed particles are claimed by their own tracks; 2 pool: this
//   event's hits removed, partner hits only; 3 pool + partner tracks. Only this event's QA tracks are counted.
//   Normalisation: ME per QA track = ME * NTrk / NTrial in each (cent9, p) bin.
#include "StMaker.h"
#include "TString.h"
#include <vector>
class StMuDstMaker;
class StRefMultCorr;
class TFile;
class TH1D;
class TH2D;
class TH3F;

class StTofMEOverlayMaker : public StMaker {
 public:
  StTofMEOverlayMaker(const char* name, StMuDstMaker* mu, const char* outName);
  virtual ~StTofMEOverlayMaker() {}
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
  void setProdMakers(StMaker* match, StMaker* vpd, StMaker* calib) { mMatch = match; mVpd = vpd; mCalib = calib; }
  void setNMix(Int_t n) { mNMix = n; }
  void setMode(Int_t m) { mMode = m; }
  void addTriggerId(UInt_t id) { mTrig.push_back(id); }
 private:
  void rerun();
  StMuDstMaker* mMu;
  StMaker *mMatch, *mVpd, *mCalib;
  TString mOutName;
  TFile* mFile;
  Int_t mNMix, mRun, mMode;
  Bool_t mBadRun;
  StRefMultCorr* mRmc;
  std::vector<UInt_t> mTrig;
  TH1D *hEvtFlow, *hDtStart, *hNAdd, *hNMerge, *hCent9, *hNPartner;
  TH2D *hNTrk, *hNTrial;
  TH3F *hM2SE, *hM2ME, *hMESt, *hM2MEmerged;
  TH3F *hQM2SE, *hQM2ME;   // charge-signed q*m2
  TH3F *hM2SEneg, *hM2MEneg; // m2 of negative tracks only
  TH2D *hDtSE, *hDtME;
  TH3F *hMEStSE;      // status in the overlay for QA tracks matched in the same event
  TH1D *hNAddTrk;
  ClassDef(StTofMEOverlayMaker, 3)
};
#endif
