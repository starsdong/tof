#ifndef StTofClosureMaker_h
#define StTofClosureMaker_h
// StTofClosureMaker: closure test of re-running StBTofMatchMaker + StVpdCalibMaker + StBTofCalibMaker on MuDst.
//   mode 0 (place BEFORE btofMatch): snapshot of the production BTOF info (per global track: matchFlag, cell,
//            tof, beta; event tStart).
//   mode 1 (place AFTER btofCalib): compares the re-run with the snapshot.
#include "StMaker.h"
#include "TString.h"
#include <vector>
class StMuDstMaker;
class TFile;
class TH1D;
class TH2D;

class StTofClosureMaker : public StMaker {
 public:
  StTofClosureMaker(const char* name, StMuDstMaker* mu, Int_t mode, const char* outName = "");
  virtual ~StTofClosureMaker() {}
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
 private:
  StMuDstMaker* mMu;
  Int_t mMode;
  TString mOutName;
  TFile* mFile;
  TH1D *hCat, *hDtof, *hDbeta, *hDtStart, *hNHit;
  TH2D *hDtofVsP, *hCatVsP;
  ClassDef(StTofClosureMaker, 1)
};
#endif
