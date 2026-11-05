#include "StTofMEQAMaker.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofHit.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"

#include <map>
#include <cmath>
#include <iostream>

ClassImp(StTofMEQAMaker)

namespace {
  const Float_t kMPi = 0.13957;
  const Float_t kMK  = 0.49368;
  const Float_t kMP  = 0.93827;
  const Float_t kTofRes = 0.013;   // nominal 1/beta resolution for nSigma_TOF
  std::map<UInt_t, Long64_t> gTrigCount;
  enum { kCutAll = 0, kCutTrig, kCutVz, kCutVr, kCutVzDiff, kNCut };
  const char* kCutName[kNCut] = {"all", "trigger", "|V_{z}|", "V_{r}", "|V_{z}-V_{z}^{VPD}|"};
}

//_____________________________________________________________________________
StTofMEQAMaker::StTofMEQAMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName)
  : StMaker(name), mPicoDstMaker(picoMaker), mPicoDst(0), mOutName(outName), mFile(0),
    mVzCut(70.), mVrCut(2.), mVzDiffCut(-1.),
    mNHitsFitMin(15), mNHitsRatioMin(0.52), mDcaMax(3.), mEtaMax(1.0)
{
}

StTofMEQAMaker::~StTofMEQAMaker() {}

//_____________________________________________________________________________
Int_t StTofMEQAMaker::Init()
{
  mFile = new TFile(mOutName.Data(), "RECREATE");
  mFile->cd();

  const Double_t pi = TMath::Pi();

  // event
  hEvtCut = new TH1F("hEvtCut", "event cut flow;;N_{events}", kNCut, 0, kNCut);
  for (int i = 0; i < kNCut; i++) hEvtCut->GetXaxis()->SetBinLabel(i + 1, kCutName[i]);
  hTrigId  = new TH1F("hTrigId", "trigger id;trigger id;N_{events}", 1, 0, 1);
  hTrigId->SetBit(TH1::kCanRebin);
  hVzAll   = new TH1F("hVzAll", "all events;V_{z} (cm);N_{events}", 400, -200, 200);
  hVz      = new TH1F("hVz", "selected;V_{z} (cm);N_{events}", 400, -200, 200);
  hVxVy    = new TH2F("hVxVy", ";V_{x} (cm);V_{y} (cm)", 200, -5, 5, 200, -5, 5);
  hVzVpd   = new TH1F("hVzVpd", ";V_{z}^{VPD} (cm);N_{events}", 400, -200, 200);
  hVzDiff  = new TH1F("hVzDiff", ";V_{z}^{TPC}-V_{z}^{VPD} (cm);N_{events}", 400, -20, 20);
  hRefMult = new TH1F("hRefMult", ";refMult;N_{events}", 600, 0, 600);
  hGRefMult= new TH1F("hGRefMult", ";gRefMult;N_{events}", 600, 0, 600);
  hTofMult = new TH1F("hTofMult", ";btofTrayMultiplicity;N_{events}", 1500, 0, 1500);
  hNBTofMatch = new TH1F("hNBTofMatch", ";nBTOFMatch;N_{events}", 500, 0, 500);
  hNBTofHits  = new TH1F("hNBTofHits", ";N BTofHit (pico);N_{events}", 1500, 0, 1500);
  hRefMultVsTofMult      = new TH2F("hRefMultVsTofMult", ";btofTrayMultiplicity;refMult", 300, 0, 1500, 300, 0, 600);
  hRefMultVsNBTofMatch   = new TH2F("hRefMultVsNBTofMatch", ";nBTOFMatch;refMult", 250, 0, 500, 300, 0, 600);
  hNBTofHitsVsNBTofMatch = new TH2F("hNBTofHitsVsNBTofMatch", ";nBTOFMatch;N BTofHit (pico)", 250, 0, 500, 300, 0, 1500);
  hNTofTrkVsNBTofHits    = new TH2F("hNTofTrkVsNBTofHits", ";N BTofHit (pico);N TOF-matched primary tracks", 300, 0, 1500, 250, 0, 500);

  // tracks
  hPt        = new TH1F("hPt", ";p_{T} (GeV/c);N_{tracks}", 200, 0, 5);
  hEta       = new TH1F("hEta", ";#eta;N_{tracks}", 200, -2, 2);
  hPhi       = new TH1F("hPhi", ";#phi;N_{tracks}", 360, -pi, pi);
  hEtaPhi    = new TH2F("hEtaPhi", ";#eta;#phi", 100, -2, 2, 180, -pi, pi);
  hNHitsFit  = new TH1F("hNHitsFit", ";nHitsFit;N_{tracks}", 80, 0, 80);
  hNHitsDedx = new TH1F("hNHitsDedx", ";nHitsDedx;N_{tracks}", 80, 0, 80);
  hNHitsRatio= new TH1F("hNHitsRatio", ";nHitsFit/nHitsMax;N_{tracks}", 110, 0, 1.1);
  hDca       = new TH1F("hDca", ";gDCA (cm);N_{tracks}", 300, 0, 10);
  hDedxVsPq  = new TH2F("hDedxVsPq", ";p#timesq (GeV/c);dE/dx (keV/cm)", 400, -4, 4, 400, 0, 20);
  hNSigPiVsP = new TH2F("hNSigPiVsP", ";p (GeV/c);n#sigma_{#pi}^{TPC}", 200, 0, 4, 200, -10, 10);

  // TOF matched
  hMatchFlag   = new TH1F("hMatchFlag", ";btofMatchFlag;N_{tracks}", 5, -0.5, 4.5);
  hBeta        = new TH1F("hBeta", ";#beta;N_{tracks}", 300, 0, 1.5);
  hInvBetaVsP  = new TH2F("hInvBetaVsP", ";p (GeV/c);1/#beta", 400, 0, 4, 400, 0.8, 3.0);
  hInvBetaVsPq = new TH2F("hInvBetaVsPq", ";p#timesq (GeV/c);1/#beta", 400, -4, 4, 400, 0.8, 3.0);
  hM2VsP       = new TH2F("hM2VsP", ";p (GeV/c);m^{2} (GeV^{2}/c^{4})", 400, 0, 4, 500, -0.5, 4.5);
  hM2VsPq      = new TH2F("hM2VsPq", ";p#timesq (GeV/c);m^{2} (GeV^{2}/c^{4})", 400, -4, 4, 500, -0.5, 4.5);
  hM2          = new TH1F("hM2", ";m^{2} (GeV^{2}/c^{4});N_{tracks}", 500, -0.5, 4.5);
  hYLocal      = new TH1F("hYLocal", ";y_{local} (cm);N_{tracks}", 200, -5, 5);
  hZLocal      = new TH1F("hZLocal", ";z_{local} (cm);N_{tracks}", 200, -5, 5);
  hYLocalVsZLocal = new TH2F("hYLocalVsZLocal", ";z_{local} (cm);y_{local} (cm)", 100, -5, 5, 100, -5, 5);
  hYLocalVsP   = new TH2F("hYLocalVsP", ";p (GeV/c);y_{local} (cm)", 100, 0, 4, 100, -5, 5);
  hZLocalVsP   = new TH2F("hZLocalVsP", ";p (GeV/c);z_{local} (cm)", 100, 0, 4, 100, -5, 5);
  hTrayModMatched = new TH2F("hTrayModMatched", "matched cells;tray;module", 120, 0.5, 120.5, 32, 0.5, 32.5);
  hHitPosPhiZ  = new TH2F("hHitPosPhiZ", "matched hit position;z (cm);#phi", 260, -260, 260, 360, -pi, pi);
  hNSigTofPiVsP = new TH2F("hNSigTofPiVsP", ";p (GeV/c);n#sigma_{#pi}^{TOF}", 200, 0, 4, 200, -10, 10);
  hNSigTofKVsP  = new TH2F("hNSigTofKVsP",  ";p (GeV/c);n#sigma_{K}^{TOF}",   200, 0, 4, 200, -10, 10);
  hNSigTofPVsP  = new TH2F("hNSigTofPVsP",  ";p (GeV/c);n#sigma_{p}^{TOF}",   200, 0, 4, 200, -10, 10);
  hDInvBetaPiVsP= new TH2F("hDInvBetaPiVsP", ";p (GeV/c);1/#beta-1/#beta_{#pi}^{exp}", 200, 0, 4, 300, -0.15, 0.45);

  hPtAll  = new TH1F("hPtAll",  "TPC tracks;p_{T} (GeV/c);N", 100, 0, 5);
  hPtTof  = new TH1F("hPtTof",  "TOF matched;p_{T} (GeV/c);N", 100, 0, 5);
  hEtaAll = new TH1F("hEtaAll", "TPC tracks;#eta;N", 100, -1.5, 1.5);
  hEtaTof = new TH1F("hEtaTof", "TOF matched;#eta;N", 100, -1.5, 1.5);
  hPhiAll = new TH1F("hPhiAll", "TPC tracks;#phi;N", 360, -pi, pi);
  hPhiTof = new TH1F("hPhiTof", "TOF matched;#phi;N", 360, -pi, pi);
  hEtaPhiTof = new TH2F("hEtaPhiTof", "TOF matched;#eta;#phi", 100, -1.5, 1.5, 180, -pi, pi);

  // BTOF hits
  hTrayModHits = new TH2F("hTrayModHits", "all BTOF hits;tray;module", 120, 0.5, 120.5, 32, 0.5, 32.5);
  hHitCell     = new TH1F("hHitCell", ";cell;N_{hits}", 6, 0.5, 6.5);

  LOG_INFO << "StTofMEQAMaker::Init  |Vz|<" << mVzCut << " Vr<" << mVrCut
           << " |Vz-VzVpd|<" << mVzDiffCut << " (neg=off)  nTrigIds=" << mTrigIds.size()
           << "  nHitsFit>=" << mNHitsFitMin << " DCA<" << mDcaMax << " |eta|<" << mEtaMax << endm;
  return kStOK;
}

//_____________________________________________________________________________
Bool_t StTofMEQAMaker::passEvent(StPicoEvent* ev)
{
  hEvtCut->Fill(Double_t(kCutAll) + 0.5);

  std::vector<unsigned int> trg = ev->triggerIds();
  for (size_t i = 0; i < trg.size(); i++) {
    gTrigCount[trg[i]]++;
    hTrigId->Fill(Form("%u", trg[i]), 1.);
  }
  TVector3 v = ev->primaryVertex();
  hVzAll->Fill(v.Z());

  if (!mTrigIds.empty()) {
    Bool_t ok = kFALSE;
    for (size_t i = 0; i < mTrigIds.size(); i++) if (ev->isTrigger(mTrigIds[i])) { ok = kTRUE; break; }
    if (!ok) return kFALSE;
  }
  hEvtCut->Fill(Double_t(kCutTrig) + 0.5);
  if (std::fabs(v.Z()) > mVzCut) return kFALSE;
  hEvtCut->Fill(Double_t(kCutVz) + 0.5);
  if (std::sqrt(v.X() * v.X() + v.Y() * v.Y()) > mVrCut) return kFALSE;
  hEvtCut->Fill(Double_t(kCutVr) + 0.5);
  hVzDiff->Fill(v.Z() - ev->vzVpd());
  if (mVzDiffCut > 0 && std::fabs(v.Z() - ev->vzVpd()) > mVzDiffCut) return kFALSE;
  hEvtCut->Fill(Double_t(kCutVzDiff) + 0.5);
  return kTRUE;
}

//_____________________________________________________________________________
Int_t StTofMEQAMaker::Make()
{
  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) { LOG_WARN << "no picoDst" << endm; return kStWarn; }
  StPicoEvent* ev = mPicoDst->event();
  if (!ev) { LOG_WARN << "no picoEvent" << endm; return kStWarn; }

  if (!passEvent(ev)) return kStOK;

  TVector3 pv = ev->primaryVertex();
  hVz->Fill(pv.Z());
  hVxVy->Fill(pv.X(), pv.Y());
  hVzVpd->Fill(ev->vzVpd());
  hRefMult->Fill(ev->refMult());
  hGRefMult->Fill(ev->grefMult());
  hTofMult->Fill(ev->btofTrayMultiplicity());
  hNBTofMatch->Fill(ev->nBTOFMatch());
  hRefMultVsTofMult->Fill(ev->btofTrayMultiplicity(), ev->refMult());
  hRefMultVsNBTofMatch->Fill(ev->nBTOFMatch(), ev->refMult());

  // ---- BTOF hits
  const Int_t nHits = mPicoDst->numberOfBTofHits();
  hNBTofHits->Fill(nHits);
  hNBTofHitsVsNBTofMatch->Fill(ev->nBTOFMatch(), nHits);
  for (Int_t i = 0; i < nHits; i++) {
    StPicoBTofHit* h = mPicoDst->btofHit(i);
    if (!h) continue;
    const Int_t id = h->id();               // (tray-1)*192+(module-1)*6+(cell-1)
    const Int_t tray = id / 192 + 1;
    const Int_t mod  = (id % 192) / 6 + 1;
    const Int_t cell = id % 6 + 1;          // StPicoBTofHit::cell() is not id%6+1
    hTrayModHits->Fill(tray, mod);
    hHitCell->Fill(cell);
  }

  // ---- tracks
  Int_t nTofTrk = 0;
  const Int_t nTrk = mPicoDst->numberOfTracks();
  for (Int_t i = 0; i < nTrk; i++) {
    StPicoTrack* t = mPicoDst->track(i);
    if (!t || !t->isPrimary()) continue;
    TVector3 mom = t->pMom();
    const Float_t pt  = mom.Perp();
    if (pt < 0.15) continue;
    const Float_t p   = mom.Mag();
    const Float_t eta = mom.PseudoRapidity();
    const Float_t phi = mom.Phi();
    const Int_t   q   = t->charge();
    const Int_t   nFit = t->nHitsFit();
    const Float_t ratio = (t->nHitsMax() > 0) ? Float_t(nFit) / t->nHitsMax() : 0.;
    const Float_t dca = t->gDCA(pv).Mag();

    hNHitsFit->Fill(nFit);
    hNHitsRatio->Fill(ratio);
    hDca->Fill(dca);
    if (nFit < mNHitsFitMin || ratio < mNHitsRatioMin || dca > mDcaMax) continue;
    hEta->Fill(eta);
    if (std::fabs(eta) > mEtaMax) continue;

    hPt->Fill(pt);
    hPhi->Fill(phi);
    hEtaPhi->Fill(eta, phi);
    hNHitsDedx->Fill(t->nHitsDedx());
    hDedxVsPq->Fill(p * q, t->dEdx());
    hNSigPiVsP->Fill(p, t->nSigmaPion());

    hPtAll->Fill(pt); hEtaAll->Fill(eta); hPhiAll->Fill(phi);

    if (!t->isTofTrack()) continue;
    StPicoBTofPidTraits* tof = mPicoDst->btofPidTraits(t->bTofPidTraitsIndex());
    if (!tof) continue;
    hMatchFlag->Fill(tof->btofMatchFlag());
    if (tof->btofMatchFlag() <= 0) continue;
    nTofTrk++;

    const Int_t cid = tof->btofCellId();
    hTrayModMatched->Fill(cid / 192 + 1, (cid % 192) / 6 + 1);
    hYLocal->Fill(tof->btofYLocal());
    hZLocal->Fill(tof->btofZLocal());
    hYLocalVsZLocal->Fill(tof->btofZLocal(), tof->btofYLocal());
    hYLocalVsP->Fill(p, tof->btofYLocal());
    hZLocalVsP->Fill(p, tof->btofZLocal());
    TVector3 hp = tof->btofHitPos();
    if (hp.Mag() > 0) hHitPosPhiZ->Fill(hp.Z(), hp.Phi());

    const Float_t beta = tof->btofBeta();
    hBeta->Fill(beta);
    if (beta <= 0) continue;
    // efficiency numerator: matched with valid beta
    hPtTof->Fill(pt); hEtaTof->Fill(eta); hPhiTof->Fill(phi); hEtaPhiTof->Fill(eta, phi);

    const Float_t ib = 1. / beta;
    const Float_t m2 = p * p * (ib * ib - 1.);
    hInvBetaVsP->Fill(p, ib);
    hInvBetaVsPq->Fill(p * q, ib);
    hM2VsP->Fill(p, m2);
    hM2VsPq->Fill(p * q, m2);
    hM2->Fill(m2);
    const Float_t ibPi = std::sqrt(p * p + kMPi * kMPi) / p;
    const Float_t ibK  = std::sqrt(p * p + kMK * kMK) / p;
    const Float_t ibP  = std::sqrt(p * p + kMP * kMP) / p;
    hDInvBetaPiVsP->Fill(p, ib - ibPi);
    hNSigTofPiVsP->Fill(p, (ib - ibPi) / kTofRes);
    hNSigTofKVsP->Fill(p, (ib - ibK) / kTofRes);
    hNSigTofPVsP->Fill(p, (ib - ibP) / kTofRes);
  }
  hNTofTrkVsNBTofHits->Fill(nHits, nTofTrk);

  return kStOK;
}

//_____________________________________________________________________________
Int_t StTofMEQAMaker::Finish()
{
  LOG_INFO << "StTofMEQAMaker::Finish  trigger id counts (all events):" << endm;
  for (std::map<UInt_t, Long64_t>::const_iterator it = gTrigCount.begin(); it != gTrigCount.end(); ++it)
    LOG_INFO << "   trgId " << it->first << " : " << it->second << endm;
  for (int i = 1; i <= hEvtCut->GetNbinsX(); i++)
    LOG_INFO << "   cut " << hEvtCut->GetXaxis()->GetBinLabel(i) << " : " << hEvtCut->GetBinContent(i) << endm;

  mFile->cd();
  hTrigId->LabelsDeflate("X");
  mFile->Write();
  mFile->Close();
  LOG_INFO << "StTofMEQAMaker: wrote " << mOutName << endm;
  return kStOK;
}
