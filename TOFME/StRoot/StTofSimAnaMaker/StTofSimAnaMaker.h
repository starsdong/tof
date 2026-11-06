#ifndef StTofSimAnaMaker_h
#define StTofSimAnaMaker_h
// StTofSimAnaMaker : truth decomposition of BTOF matches in a simulation MuDst (HIJING + GEANT + bfc).
// For every TOF-matched QA primary track (same QA cuts as StTofMEOverlayMaker) the matched hit is classified with the
// MC record (StMuMcTrack / StMuMcVertex, hit and track idTruth):
//   0 own        hit made by the track's own MC particle
//   1 decay      hit made by a descendant via a decay (weak decayer -> mu, or decay vertex at R < 190 cm)
//   2 secondary  hit made by a descendant, first step an interaction (hadronic, conversion, delta ray, ...)
//   3 relative   hit made by an ancestor or a sibling-branch particle (common non-primary ancestor)
//   4 unrelated  no common ancestry (combinatorial, the part a mixed event models)
//   5 no truth   hit or track without idTruth
#include "StMaker.h"
#include "TString.h"
class StMuDstMaker;
class StMuDst;
class TFile;
class TH1D;
class TH2D;
class TH3F;

class StTofSimAnaMaker : public StMaker {
public:
  StTofSimAnaMaker(const char* name, StMuDstMaker* mu, const char* outName);
  Int_t Init();
  Int_t Make();
  Int_t Finish();
  enum { kNCat = 6 };

private:
  Int_t classify(StMuDst* d, Int_t trkId, Int_t hitId, Float_t& rVx, Int_t& gePid);
  StMuDstMaker* mMu;
  TString mOutName;
  TFile* mFile;
  TH1D* hEvt;            // 0 all, 1 |Vz|<70, 2 Vr<2, 3 with QA tracks
  TH1D* hRefMult;
  TH2D* hNTrk;           // refMult, p
  TH3F* hQM2[kNCat];     // refMult, p, q*m2 (true start time)
  TH3F* hQM2cal;         // all categories, calibrated (startless) timeOfFlight
  TH2D* hDt[kNCat];      // p, t - t_pi (true start)
  TH2D* hRVx[kNCat];     // p, radius of the vertex where the hit particle's branch leaves the track's family
  TH2D* hPid[kNCat];     // p, GEANT id of the particle that made the hit
  TH2D* hTStart;         // refMult, tStart(calib) - McClock
  TH2D* hTot[kNCat];     // p, ToT of the matched hit (ns)
  TH2D* hTotHi[kNCat];   // p, ToT of the matched hit for 1.6 < |q m2| < 2.8
  ClassDef(StTofSimAnaMaker, 2)
};
#endif
