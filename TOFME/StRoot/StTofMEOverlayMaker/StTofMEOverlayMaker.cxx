#include "StTofMEOverlayMaker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StEvent/StTriggerId.h"
#include "StEvent/StRunInfo.h"
#include "StEvent/StBTofHeader.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include <deque>
#include <vector>
#include <map>
#include <cmath>

ClassImp(StTofMEOverlayMaker)

namespace {
  const Int_t kNVz = 14, kNMult = 9;
  const Int_t kMultEdge[kNMult + 1] = {0, 10, 25, 50, 80, 120, 170, 230, 300, 100000};
  const Double_t kTMax = 51200., kC = 29.9792458, kMPi = 0.13957;
  const Int_t kNP = 64; const Double_t kPMax = 3.2;
  const Int_t kNM2 = 600; const Double_t kM2Lo = -1.0, kM2Hi = 5.0;
  struct HitR { UChar_t tray, module, cell; Double_t le, te; };   // LE/TE relative to the event tStart
  struct Mix { std::vector<HitR> hits; std::vector<StMuTrack> trk; };
  std::deque<Mix> gBuf[kNVz][kNMult];
  Int_t multClass(Int_t r) { for (int i = 0; i < kNMult; i++) if (r < kMultEdge[i + 1]) return i; return kNMult - 1; }
  Double_t wrap(Double_t t) { while (t < 0) t += kTMax; while (t >= kTMax) t -= kTMax; return t; }
  Int_t cellId(Int_t tray, Int_t module, Int_t cell) { return ((tray - 1) * 32 + (module - 1)) * 6 + (cell - 1); }
  Double_t pathL(const StThreeVectorF& b, const StThreeVectorF& e, Double_t curv)
  {
    const Double_t xd = e.x() - b.x(), yd = e.y() - b.y();
    const Double_t C = std::sqrt(xd * xd + yd * yd);
    Double_t sPerp = C;
    if (curv != 0) { const Double_t R = 1. / curv; sPerp = 2 * R * std::asin(C / (2 * R)); }
    const Double_t sZ = std::fabs(e.z() - b.z());
    return std::sqrt(sPerp * sPerp + sZ * sZ);
  }
  // m2 of a QA primary track from its global track's (re-run) BTOF pid traits; returns false if not matched
  Bool_t trackM2(const StMuTrack* t, const StThreeVectorF& pv, Double_t& m2, Double_t& dt)
  {
    const StMuTrack* g = t->globalTrack();
    if (!g) return kFALSE;
    const StMuBTofPidTraits& pid = g->btofPidTraits();
    if (pid.matchFlag() <= 0 || pid.timeOfFlight() <= 0) return kFALSE;
    const Double_t p = t->p().mag();
    const Double_t L = pathL(pv, pid.position(), t->helix().curvature());
    const Double_t beta = L / (kC * pid.timeOfFlight());
    m2 = p * p * (1. / (beta * beta) - 1.);
    dt = pid.timeOfFlight() - L / kC * std::sqrt(1. + kMPi * kMPi / (p * p));
    return kTRUE;
  }
  Bool_t qaTrack(const StMuTrack* t)
  {
    if (!t) return kFALSE;
    const StThreeVectorF pm = t->p();
    if (pm.perp() < 0.15 || std::fabs(pm.pseudoRapidity()) > 1.0) return kFALSE;
    const Int_t nF = t->nHitsFit(), nP = t->nHitsPoss();
    const Float_t r = nP > 0 ? Float_t(nF) / nP : 0.;
    if (nF < 15 || r < 0.52f || t->dcaGlobal().mag() > 3.) return kFALSE;
    return t->globalTrack() != 0;
  }
}

StTofMEOverlayMaker::StTofMEOverlayMaker(const char* name, StMuDstMaker* mu, const char* outName)
  : StMaker(name), mMu(mu), mMatch(0), mVpd(0), mCalib(0), mOutName(outName), mFile(0), mNMix(5), mRun(-1), mMode(0),
    mBadRun(kFALSE), mRmc(0) {}

Int_t StTofMEOverlayMaker::Init()
{
  mFile = new TFile(mOutName.Data(), "RECREATE");
  hEvtFlow = new TH1D("hEvtFlow", "0 all,1 trig,2 |Vz|<70,3 Vr<2,4 tStart ok,5 m2 (cent ok, no pile-up);;events", 6, -0.5, 5.5);
  hDtStart = new TH1D("hDtStart", "tStart_{A} - tStart_{B} (raw);ns;pairs", 400, -200, 200);
  hNAdd = new TH1D("hNAdd", "partner hits added per overlay;hits;overlays", 300, 0, 1500);
  hNMerge = new TH1D("hNMerge", "partner hits merged into this-event hits;hits;overlays", 100, 0, 100);
  hNPartner = new TH1D("hNPartner", ";partners;events", 11, -0.5, 10.5);
  hCent9 = new TH1D("hCent9", ";cent9;events", 9, -0.5, 8.5);
  hNTrk = new TH2D("hNTrk", "QA tracks;cent9;p (GeV/c)", 9, -0.5, 8.5, kNP, 0, kPMax);
  hNTrial = new TH2D("hNTrial", "QA tracks x partners;cent9;p (GeV/c)", 9, -0.5, 8.5, kNP, 0, kPMax);
  const char* ax = ";cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})";
  hM2SE = new TH3F("hM2SE", Form("SE (re-run production)%s", ax), 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hM2ME = new TH3F("hM2ME", Form("ME: matched to a partner hit%s", ax), 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hM2MEmerged = new TH3F("hM2MEmerged", Form("ME: matched to a merged hit with partner leading edge%s", ax), 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hMESt = new TH3F("hMESt", "overlay: 0 own-event hit, 1 partner hit, 2 no match;cent9;p (GeV/c);status", 9, -0.5, 8.5, kNP, 0, kPMax, 3, -0.5, 2.5);
  hDtSE = new TH2D("hDtSE", "SE;p (GeV/c);t - t_{#pi} (ns)", 32, 0, 3.2, 600, -3, 27);
  hDtME = new TH2D("hDtME", "ME;p (GeV/c);t - t_{#pi} (ns)", 32, 0, 3.2, 600, -3, 27);
  hMEStSE = new TH3F("hMEStSE", "QA tracks matched in SE, overlay status: 0 own-event hit, 1 partner hit, 2 no match;cent9;p (GeV/c);status", 9, -0.5, 8.5, kNP, 0, kPMax, 3, -0.5, 2.5);
  hNAddTrk = new TH1D("hNAddTrk", "partner global tracks added per overlay;tracks;overlays", 300, 0, 3000);
  mRmc = CentralityMaker::instance()->getRefMultCorr();
  LOG_INFO << "StTofMEOverlayMaker v2: mode " << mMode << " nMix " << mNMix << " sizeof(StMuTrack) " << sizeof(StMuTrack) << endm;
  return kStOK;
}

void StTofMEOverlayMaker::rerun()
{
  mMatch->Make();
  mVpd->Make();
  mCalib->Make();
}

Int_t StTofMEOverlayMaker::Make()
{
  StMuDst* d = mMu->muDst();
  if (!d || !d->event()) return kStOK;
  StMuEvent* ev = d->event();
  hEvtFlow->Fill(0);
  Bool_t trig = mTrig.empty();
  for (size_t i = 0; i < mTrig.size() && !trig; i++) if (ev->triggerIdCollection().nominal().isTrigger(mTrig[i])) trig = kTRUE;
  if (!trig) return kStOK;
  hEvtFlow->Fill(1);
  const StThreeVectorF pv = ev->primaryVertexPosition();
  if (std::fabs(pv.z()) >= 70.) return kStOK;
  hEvtFlow->Fill(2);
  if (pv.perp() >= 2.) return kStOK;
  hEvtFlow->Fill(3);
  const StBTofHeader* hd = d->btofHeader();
  if (!hd || hd->tStart() <= 0 || hd->tStart() > kTMax) return kStOK;    // no VPD start
  const Double_t tsA = hd->tStart();
  hEvtFlow->Fill(4);

  const Int_t run = ev->runId();
  if (run != mRun) {
    for (int i = 0; i < kNVz; i++) for (int j = 0; j < kNMult; j++) gBuf[i][j].clear();
    mRun = run; mRmc->init(run); mBadRun = mRmc->isBadRun(run);
  }
  const Int_t refMult = ev->refMult();
  Int_t nBTOFMatch = 0;
  const Int_t nPrim = d->numberOfPrimaryTracks();
  for (Int_t i = 0; i < nPrim; i++) { const StMuTrack* t = d->primaryTracks(i); if (t && t->btofPidTraits().matchFlag() > 0) nBTOFMatch++; }
  mRmc->initEvent(refMult, pv.z(), ev->runInfo().zdcCoincidenceRate());
  const Int_t cent9 = mRmc->getCentralityBin9();
  const Bool_t useM2 = !mBadRun && !mRmc->isPileUpEvent(refMult, nBTOFMatch, pv.z()) && cent9 >= 0 && cent9 <= 8;
  Int_t iVz = (Int_t)((pv.z() + 70.) / 140. * kNVz); if (iVz < 0) iVz = 0; if (iVz >= kNVz) iVz = kNVz - 1;
  std::deque<Mix>& buf = gBuf[iVz][multClass(refMult)];

  // snapshot of this event's hits (state after the chain's re-run on A alone)
  TClonesArray* arr = StMuDst::btofArray(muBTofHit);
  const Int_t n0 = arr->GetEntriesFast();
  std::vector<StMuBTofHit> orig; orig.reserve(n0);
  Mix cur; cur.hits.reserve(n0);
  std::map<Int_t, Int_t> cellToHit;
  for (Int_t i = 0; i < n0; i++) {
    const StMuBTofHit* h = (StMuBTofHit*)arr->UncheckedAt(i);
    orig.push_back(*h);
    HitR r; r.tray = h->tray(); r.module = h->module(); r.cell = h->cell();
    r.le = h->leadingEdgeTime() - tsA; r.te = h->trailingEdgeTime() - tsA;
    cur.hits.push_back(r);
    cellToHit[cellId(h->tray(), h->module(), h->cell())] = i;
  }
  TClonesArray* garr = StMuDst::array(muGlobal);
  const Int_t nG0 = garr->GetEntriesFast();
  if (mMode == 1 || mMode == 3) {
    for (Int_t i = 0; i < nG0; i++) {
      const StMuTrack* g = (const StMuTrack*)garr->UncheckedAt(i);
      if (!g || g->flag() <= 0 || g->flag() >= 1000 || g->nHitsFit(kTpcId) < 15) continue;
      cur.trk.push_back(*g);
    }
  }

  if (useM2) {
    hEvtFlow->Fill(5); hCent9->Fill(cent9); hNPartner->Fill(buf.size());
    std::vector<Int_t> qa; std::vector<Bool_t> qaSE;
    for (Int_t i = 0; i < nPrim; i++) {
      const StMuTrack* t = d->primaryTracks(i);
      if (!qaTrack(t)) continue;
      qa.push_back(i);
      const Double_t p = t->p().mag();
      hNTrk->Fill(cent9, p);
      Double_t m2, dt;
      const Bool_t mSE = trackM2(t, pv, m2, dt);
      qaSE.push_back(mSE);
      if (mSE) { hM2SE->Fill(cent9, p, m2); hDtSE->Fill(p, dt); }
    }
    for (size_t k = 0; k < buf.size(); k++) {
      const Mix& B = buf[k];
      // overlay
      std::vector<Int_t> merged;     // A hit indices whose leading edge now comes from B
      Int_t nAdd = 0;
      const Bool_t pool = (mMode == 2 || mMode == 3);
      if (pool) { arr->Clear(); cellToHit.clear(); }
      const Int_t nA = pool ? 0 : n0;  // hits with index >= nA come from the partner
      for (size_t t = 0; t < B.trk.size(); t++) new ((*garr)[garr->GetEntriesFast()]) StMuTrack(B.trk[t]);
      hNAddTrk->Fill(B.trk.size());
      for (size_t j = 0; j < B.hits.size(); j++) {
        const HitR& b = B.hits[j];
        const Double_t le = wrap(b.le + tsA), te = wrap(b.te + tsA);
        std::map<Int_t, Int_t>::const_iterator it = cellToHit.find(cellId(b.tray, b.module, b.cell));
        if (it != cellToHit.end()) {
          StMuBTofHit* a = (StMuBTofHit*)arr->UncheckedAt(it->second);
          const Double_t leA = a->leadingEdgeTime(), totA = a->tot(), totB = te - le;
          const Bool_t bFirst = le < leA;
          if (std::fabs(le - leA) < (bFirst ? totB : totA)) {     // overlapping pulses -> one hit
            if (bFirst) { a->setLeadingEdgeTime(le); merged.push_back(it->second); }
            if (te > a->trailingEdgeTime()) a->setTrailingEdgeTime(te);
            continue;
          }
        }
        StMuBTofHit* n = new ((*arr)[arr->GetEntriesFast()]) StMuBTofHit();
        n->setTrayIndex(b.tray); n->setModuleIndex(b.module); n->setCellIndex(b.cell);
        n->setLeadingEdgeTime(le); n->setTrailingEdgeTime(te);
        n->setAssociatedTrackId(-1); n->setIndex2Primary(-1); n->setIndex2Global(-1);
        nAdd++;
      }
      hNAdd->Fill(nAdd); hNMerge->Fill(merged.size());
      rerun();
      // which hit did each global track get?
      const Int_t nh = arr->GetEntriesFast();
      std::map<Int_t, Int_t> hitOfGlobal;
      for (Int_t i = 0; i < nh; i++) { const StMuBTofHit* h = (StMuBTofHit*)arr->UncheckedAt(i); if (h->index2Global() >= 0) hitOfGlobal[h->index2Global()] = i; }
      for (size_t q = 0; q < qa.size(); q++) {
        const StMuTrack* t = d->primaryTracks(qa[q]);
        const Double_t p = t->p().mag();
        hNTrial->Fill(cent9, p);
        Double_t m2, dt;
        if (!trackM2(t, pv, m2, dt)) { hMESt->Fill(cent9, p, 2); if (qaSE[q]) hMEStSE->Fill(cent9, p, 2); continue; }
        std::map<Int_t, Int_t>::const_iterator it = hitOfGlobal.find(t->index2Global());
        const Int_t hi = it != hitOfGlobal.end() ? it->second : -1;
        Bool_t fromB = hi >= nA;
        Bool_t fromMerge = kFALSE;
        for (size_t m = 0; m < merged.size() && !fromB; m++) if (merged[m] == hi) { fromB = kTRUE; fromMerge = kTRUE; }
        hMESt->Fill(cent9, p, fromB ? 1 : 0);
        if (qaSE[q]) hMEStSE->Fill(cent9, p, fromB ? 1 : 0);
        if (fromB) { hM2ME->Fill(cent9, p, m2); hDtME->Fill(p, dt); if (fromMerge) hM2MEmerged->Fill(cent9, p, m2); }
      }
      // restore A's hit array
      for (Int_t i = garr->GetEntriesFast() - 1; i >= nG0; i--) garr->RemoveAt(i);
      if (pool) {
        arr->Clear();
        for (Int_t i = 0; i < n0; i++) new ((*arr)[i]) StMuBTofHit(orig[i]);
        for (Int_t i = 0; i < n0; i++) { const StMuBTofHit* h = (StMuBTofHit*)arr->UncheckedAt(i); cellToHit[cellId(h->tray(), h->module(), h->cell())] = i; }
      } else {
        for (Int_t i = arr->GetEntriesFast() - 1; i >= n0; i--) arr->RemoveAt(i);
        for (Int_t i = 0; i < n0; i++) *((StMuBTofHit*)arr->UncheckedAt(i)) = orig[i];
      }
    }
  }
  buf.push_back(cur);
  while ((Int_t)buf.size() > mNMix) buf.pop_front();
  return kStOK;
}

Int_t StTofMEOverlayMaker::Finish()
{
  LOG_INFO << "StTofMEOverlayMaker::Finish mode " << mMode << " events " << hEvtFlow->GetBinContent(1) << " tStart ok " << hEvtFlow->GetBinContent(5)
           << " m2 " << hEvtFlow->GetBinContent(6) << " QA tracks " << hNTrk->GetEntries() << " trials " << hNTrial->GetEntries()
           << " SE " << hM2SE->GetEntries() << " ME " << hM2ME->GetEntries() << " (merged " << hM2MEmerged->GetEntries() << ")" << endm;
  mFile->cd(); mFile->Write(); mFile->Close();
  return kStOK;
}
