#include "StTofClosureMaker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
#include "StEvent/StBTofHeader.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <cmath>

ClassImp(StTofClosureMaker)

namespace {
  struct Snap { Int_t flag, cell; Float_t tof, beta; };
  std::vector<Snap> gSnap;      // per global track
  Double_t gTStart = 0;
  Int_t gNHit = 0;
  void cellOfGlobal(const StMuDst* d, std::vector<Int_t>& v)
  {
    v.assign(d->numberOfGlobalTracks(), -1);
    for (UInt_t i = 0; i < d->numberOfBTofHit(); i++) {
      const StMuBTofHit* h = d->btofHit(i);
      const Int_t ig = h->index2Global();
      if (ig >= 0 && ig < (Int_t)v.size()) v[ig] = ((h->tray() - 1) * 32 + (h->module() - 1)) * 6 + (h->cell() - 1);
    }
  }
}

StTofClosureMaker::StTofClosureMaker(const char* name, StMuDstMaker* mu, Int_t mode, const char* outName)
  : StMaker(name), mMu(mu), mMode(mode), mOutName(outName), mFile(0) {}

Int_t StTofClosureMaker::Init()
{
  if (mMode != 1) return kStOK;
  mFile = new TFile(mOutName.Data(), "RECREATE");
  hCat = new TH1D("hCat", "0 same cell, 1 other cell, 2 only production, 3 only rerun;;global tracks", 4, -0.5, 3.5);
  hCatVsP = new TH2D("hCatVsP", ";p (GeV/c);category", 30, 0, 3, 4, -0.5, 3.5);
  hDtof = new TH1D("hDtof", "same cell: tof_{rerun} - tof_{prod};ns;tracks", 2000, -1, 1);
  hDtofVsP = new TH2D("hDtofVsP", "same cell;p (GeV/c);tof_{rerun} - tof_{prod} (ns)", 30, 0, 3, 400, -1, 1);
  hDbeta = new TH1D("hDbeta", "same cell: 1/#beta_{rerun} - 1/#beta_{prod};;tracks", 2000, -0.05, 0.05);
  hDtStart = new TH1D("hDtStart", "tStart_{rerun} - tStart_{prod};ns;events", 2000, -5, 5);
  hNHit = new TH1D("hNHit", "BTOF hits rerun - production;;events", 21, -10.5, 10.5);
  return kStOK;
}

Int_t StTofClosureMaker::Make()
{
  StMuDst* d = mMu->muDst();
  if (!d || !d->event()) return kStOK;
  std::vector<Int_t> cell;
  cellOfGlobal(d, cell);
  const Int_t nG = d->numberOfGlobalTracks();
  const StBTofHeader* hd = d->btofHeader();
  if (mMode == 0) {
    gSnap.assign(nG, Snap());
    for (Int_t i = 0; i < nG; i++) {
      const StMuTrack* t = d->globalTracks(i);
      Snap s; s.flag = 0; s.cell = -1; s.tof = 0; s.beta = 0;
      if (t) { const StMuBTofPidTraits& p = t->btofPidTraits(); s.flag = p.matchFlag(); s.tof = p.timeOfFlight(); s.beta = p.beta(); s.cell = cell[i]; }
      gSnap[i] = s;
    }
    gTStart = hd ? hd->tStart() : -999;
    gNHit = d->numberOfBTofHit();
    return kStOK;
  }
  if ((Int_t)gSnap.size() != nG) return kStOK;
  if (hd) hDtStart->Fill(hd->tStart() - gTStart);
  hNHit->Fill((Int_t)d->numberOfBTofHit() - gNHit);
  for (Int_t i = 0; i < nG; i++) {
    const StMuTrack* t = d->globalTracks(i);
    if (!t) continue;
    const StMuBTofPidTraits& p = t->btofPidTraits();
    const Bool_t a = gSnap[i].flag > 0, b = p.matchFlag() > 0;
    if (!a && !b) continue;
    const Double_t mom = t->p().mag();
    Int_t c = 3;
    if (a && b) c = (gSnap[i].cell == cell[i]) ? 0 : 1;
    else if (a) c = 2;
    hCat->Fill(c); hCatVsP->Fill(mom, c);
    if (c == 0 && gSnap[i].tof > 0 && p.timeOfFlight() > 0) {
      hDtof->Fill(p.timeOfFlight() - gSnap[i].tof); hDtofVsP->Fill(mom, p.timeOfFlight() - gSnap[i].tof);
      if (gSnap[i].beta > 0 && p.beta() > 0) hDbeta->Fill(1. / p.beta() - 1. / gSnap[i].beta);
    }
  }
  return kStOK;
}

Int_t StTofClosureMaker::Finish()
{
  if (mMode != 1) return kStOK;
  LOG_INFO << "StTofClosureMaker: tracks same cell " << hCat->GetBinContent(1) << " other cell " << hCat->GetBinContent(2)
           << " only prod " << hCat->GetBinContent(3) << " only rerun " << hCat->GetBinContent(4)
           << " | dtof mean " << hDtof->GetMean() << " rms " << hDtof->GetRMS() << " | dtStart mean " << hDtStart->GetMean()
           << " rms " << hDtStart->GetRMS() << endm;
  mFile->cd(); mFile->Write(); mFile->Close();
  return kStOK;
}
