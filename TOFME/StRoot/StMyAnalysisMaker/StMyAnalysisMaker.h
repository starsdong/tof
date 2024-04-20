#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

#include "StMaker.h"
class StPicoDst;
class StPicoDstMaker;
class TString;
class TH1F;
class TH2F;
class TRandom3;

class StMyAnalysisMaker : public StMaker {
  public:
    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StMyAnalysisMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    void    DeclareHistograms();
    void    WriteHistograms();
    
  private:
    static const Int_t NB = 10;
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    
    TString    mOutName;
    TRandom3   *gRandom;
    
    TH1F *hEventCount;
    TH1F *hVz;
    TH2F *hVxVy;
    TH1F *hVzGoodVxy;
    TH1F *hTriggerIds;
    TH1F *hRefMult;
    
    TH2F *hdEdxvsp;
    TH2F *hInvBetavsp;
    TH2F *hMass2vsp;
    
    TH2F *hMass2vspME;
    
    vector<Int_t> mMixPoolId[NB];
    vector<Float_t> mMixPoolTof[NB];
                    
    ClassDef(StMyAnalysisMaker, 1)
};

#endif
