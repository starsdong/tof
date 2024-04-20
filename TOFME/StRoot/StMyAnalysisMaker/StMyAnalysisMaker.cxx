#include "StMyAnalysisMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TRandom3.h"

ClassImp(StMyAnalysisMaker)

//-----------------------------------------------------------------------------
StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mOutName = outName;
  for(int i=0;i<NB;i++) {
    mMixPoolId[i].clear();
    mMixPoolTof[i].clear();
  }
  gRandom = new TRandom3();
}

//----------------------------------------------------------------------------- 
StMyAnalysisMaker::~StMyAnalysisMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Init() {
  DeclareHistograms();
  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Finish() {
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),"RECREATE");
    fout->cd();
    WriteHistograms();
    fout->Close();
  }
  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::DeclareHistograms() {
    hEventCount = new TH1F("eventCount","",20,0,20);
    hVz = new TH1F("vz","",3000,-300.,300.);
    hVzGoodVxy = new TH1F("vzGoodVxy","",3000,-300.,300.);
    hVxVy = new TH2F("vxvy","",250,-5.,5.,250,-5.,5.);
    hTriggerIds = new TH1F("triggerIds","",10,0,10);
    hRefMult = new TH1F("refMult","",1000,0,1000);
    hdEdxvsp = new TH2F("dEdxvsp","",1000,-5.,5.,500, 0., 10.);
    hInvBetavsp = new TH2F("InvBetavsp","",1000,-5.,5.,500, 0., 5.);    
    hMass2vsp = new TH2F("Mass2vsp","",1000,-5.,5.,500, -0.5, 2.0);
    hMass2vspME = new TH2F("Mass2vspME","",1000,-5.,5.,500, -0.5, 2.0);
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::WriteHistograms() {
  hEventCount->Write();
  hVz->Write();
  hVzGoodVxy->Write();
  hVxVy->Write();
  hTriggerIds->Write();
  hRefMult->Write();
  hdEdxvsp->Write();
  hInvBetavsp->Write();
  hMass2vsp->Write();
  hMass2vspME->Write();
}

//----------------------------------------------------------------------------- 
void StMyAnalysisMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Make() {
  hEventCount->Fill(0);
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  hEventCount->Fill(1);

  StPicoEvent *mPicoEvent = mPicoDst->event();
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }
  hEventCount->Fill(2);

/*  
  const Int_t nId = 4;
  Bool_t isMB = kFALSE;
  static const Int_t mTriggerIds[nId] = {810010, 810020, 810030, 810040};
  for(int i=0;i<nId;i++) {
    if(mPicoEvent->isTrigger(mTriggerIds[i])) {
      isMB = kTRUE;
      hTriggerIds->Fill(i);
    }    
  }
  
  if(!isMB) return kStOK;
  */
  hEventCount->Fill(3);
 
  TVector3 pVtx = mPicoEvent->primaryVertex();
  float vz = pVtx.Z();
  int refMult = mPicoEvent->refMult();
  hVz->Fill(pVtx[2]);
  hVxVy->Fill(pVtx[0], pVtx[1]);

// beam spot offset (-0.27, -0.22)
  float dVr = TMath::Sqrt((pVtx[0] + 0.27)*(pVtx[0] + 0.27) + (pVtx[1] + 0.22)*(pVtx[1] + 0.22));  
  if(dVr>1.0) return kStOK;
  hEventCount->Fill(4);
  hVzGoodVxy->Fill(pVtx[2]);
  if(fabs(vz)>200.) return kStOK;
  hEventCount->Fill(5);
  if(fabs(vz)>150.) return kStOK;
  hEventCount->Fill(6);
  if(fabs(vz)>70.) return kStOK;
  hEventCount->Fill(7);
  if(fabs(vz)>50.) return kStOK;
  hEventCount->Fill(8);
  if(fabs(vz)>20.) return kStOK;
  hEventCount->Fill(9);

  hRefMult->Fill(refMult);


  for(int i=0;i<NB;i++) {
//    LOG_INFO << " buffer " << i << " size=" << mMixPoolId[i].size() << " " << mMixPoolTof[i].size() << endm;
  }
   
  int i_Buffer = (Int_t)(gRandom->Rndm()/0.1);
//  LOG_INFO << " buffer index " << i_Buffer << endm;
  
  for(size_t i=0;i<mPicoDst->numberOfTracks();i++) {
    StPicoTrack *t = (StPicoTrack *)mPicoDst->track(i);
    if(!t) continue;
    TVector3 pMom = t->pMom();
    float p = pMom.Mag();
    if(p<0.1) continue;
    int q = t->charge();
    float dEdx = t->dEdx();
    hdEdxvsp->Fill(p*q, dEdx);
    
    int btofIndex = t->bTofPidTraitsIndex();
    if(btofIndex>=0) {
      StPicoBTofPidTraits *trait = mPicoDst->btofPidTraits(btofIndex);
      if(!trait) continue;
      float beta = trait->btofBeta();
      if(beta<=0) continue;
      hInvBetavsp->Fill(p*q, 1./beta);
      
      float m2 = p*p*(1./beta/beta - 1);
      hMass2vsp->Fill(p*q, m2);
      float btof = trait->btof();
      
      // Mixed-event matching
      int thisCellId = trait->btofCellId();
      for(int j=0;j<NB;j++) {
        for(size_t k=0;k<mMixPoolId[j].size();k++) {
          if(thisCellId != mMixPoolId[j][k]) continue;
//          LOG_INFO << " CellId Match!!!  TOF  = " << btof << " \t " << mMixPoolTof[j][k] << endm;
          float beta_p = beta * btof / mMixPoolTof[j][k];   // scale beta
          float m2_p = p*p*(1./beta_p/beta_p - 1);
          hMass2vspME->Fill(p*q, m2_p);
        }
      }
    }
  }  // end track loop


  // push this event into the buffer region
  if(mMixPoolId[i_Buffer].size()!=0) {
    mMixPoolId[i_Buffer].clear();
    mMixPoolTof[i_Buffer].clear();
  }

  for(size_t j=0;j<mPicoDst->numberOfBTofPidTraits();j++) {
    StPicoBTofPidTraits *t = mPicoDst->btofPidTraits(j);
    if(!t) continue;

    mMixPoolId[i_Buffer].push_back(t->btofCellId());
    mMixPoolTof[i_Buffer].push_back(t->btof());
  }
  
  return kStOK;
}

