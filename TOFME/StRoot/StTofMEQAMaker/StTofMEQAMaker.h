#ifndef StTofMEQAMaker_h
#define StTofMEQAMaker_h

// StTofMEQAMaker
//   Step 1 of the BTOF mixed-event matching study (TOFME):
//   read picoDst and fill basic QA histograms for events, TPC tracks,
//   BTOF hits and TPC-TOF matched (BTofPidTraits) information.

#include "StMaker.h"
#include "TString.h"
#include <vector>

class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class TFile;
class TH1F;
class TH2F;
class TProfile;

class StTofMEQAMaker : public StMaker {
 public:
  StTofMEQAMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);
  virtual ~StTofMEQAMaker();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();

  // configuration (call before Init)
  void setVzCut(Float_t v)          { mVzCut = v; }
  void setVrCut(Float_t v)          { mVrCut = v; }
  void setVzDiffCut(Float_t v)      { mVzDiffCut = v; }
  void addTriggerId(UInt_t id)      { mTrigIds.push_back(id); }
  void setNHitsFitMin(Int_t n)      { mNHitsFitMin = n; }
  void setDcaMax(Float_t v)         { mDcaMax = v; }
  void setEtaMax(Float_t v)         { mEtaMax = v; }

 private:
  Bool_t passEvent(StPicoEvent* ev);

  StPicoDstMaker* mPicoDstMaker;
  StPicoDst*      mPicoDst;
  TString         mOutName;
  TFile*          mFile;

  // cuts
  Float_t mVzCut, mVrCut, mVzDiffCut;
  std::vector<UInt_t> mTrigIds;
  Int_t   mNHitsFitMin;
  Float_t mNHitsRatioMin, mDcaMax, mEtaMax;

  // ---- event level
  TH1F* hEvtCut;        // cut flow
  TH1F* hTrigId;        // trigger ids (all events)
  TH1F* hVzAll;
  TH1F* hVz;
  TH2F* hVxVy;
  TH1F* hVzVpd;
  TH1F* hVzDiff;        // vz(TPC) - vz(VPD)
  TH1F* hRefMult;
  TH1F* hGRefMult;
  TH1F* hTofMult;       // btofTrayMultiplicity
  TH1F* hNBTofMatch;
  TH1F* hNBTofHits;     // # entries in BTofHit array
  TH2F* hRefMultVsTofMult;
  TH2F* hRefMultVsNBTofMatch;
  TH2F* hNBTofHitsVsNBTofMatch;
  TH2F* hNTofTrkVsNBTofHits;  // # TOF-matched tracks vs # BTof hits

  // ---- track level (primary, quality cuts)
  TH1F* hPt;
  TH1F* hEta;
  TH1F* hPhi;
  TH2F* hEtaPhi;
  TH1F* hNHitsFit;
  TH1F* hNHitsDedx;
  TH1F* hNHitsRatio;
  TH1F* hDca;
  TH2F* hDedxVsPq;
  TH2F* hNSigPiVsP;

  // ---- TOF matched tracks
  TH1F* hMatchFlag;
  TH1F* hBeta;
  TH2F* hInvBetaVsP;
  TH2F* hInvBetaVsPq;
  TH2F* hM2VsP;
  TH2F* hM2VsPq;
  TH1F* hM2;
  TH1F* hYLocal;
  TH1F* hZLocal;
  TH2F* hYLocalVsZLocal;
  TH2F* hYLocalVsP;
  TH2F* hZLocalVsP;
  TH2F* hTrayModMatched;   // matched cell: tray vs module
  TH2F* hHitPosPhiZ;       // matched hit global position phi vs z
  TH2F* hNSigTofPiVsP;
  TH2F* hNSigTofKVsP;
  TH2F* hNSigTofPVsP;
  TH2F* hDInvBetaPiVsP;    // 1/beta - 1/beta_pi(expected)

  // TOF matching efficiency (num/den) for |eta|<etaMax
  TH1F* hPtAll;   TH1F* hPtTof;
  TH1F* hEtaAll;  TH1F* hEtaTof;
  TH1F* hPhiAll;  TH1F* hPhiTof;
  TH2F* hEtaPhiTof;

  // ---- BTOF hits (all hits in event)
  TH2F* hTrayModHits;      // tray vs module occupancy (all hits)
  TH1F* hHitCell;          // cell 1-6

  ClassDef(StTofMEQAMaker, 1)
};

#endif
