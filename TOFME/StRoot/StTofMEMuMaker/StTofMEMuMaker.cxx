#include "StTofMEMuMaker.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StEvent/StTriggerId.h"
#include "StEvent/StRunInfo.h"
#include "StEvent/StEnumerations.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TVector3.h"

#include <deque>
#include <vector>
#include <cmath>
#include <algorithm>

ClassImp(StTofMEMuMaker)

namespace {
  const Int_t kNVzM = 14, kNMultM = 9;
  const Int_t kNCellM = StTofMEMaker::kNCell;
  const Double_t kCm = 29.9792458, kMPiM = 0.13957;
  const Int_t kNPm = 64; const Double_t kPMaxM = 3.2;
  const Int_t kNM2m = 600; const Double_t kM2LoM = -1.0, kM2HiM = 5.0;
  struct MuMix {
    std::vector<float> tTot, totMax, tEarly;   // per cell: time of the largest-ToT hit, its ToT, earliest time (0 = no hit)
    std::vector<unsigned char> mult;           // hits in the cell
    std::vector<unsigned char> own;            // projections of this event's production-valid tracks covering the cell
    Int_t nHits;
  };
  std::deque<MuMix> gBufM[kNVzM][kNMultM];
  MuMix gCur;
  Double_t gOff[kNCellM];
  Double_t pathL(const TVector3& b, const TVector3& e, Double_t curv)
  {
    const Double_t xd = e.X() - b.X(), yd = e.Y() - b.Y();
    const Double_t C = std::sqrt(xd * xd + yd * yd);
    Double_t sPerp = C;
    if (curv != 0) { const Double_t R = 1. / curv; sPerp = 2 * R * std::asin(C / (2 * R)); }
    const Double_t sZ = std::fabs(e.Z() - b.Z());
    return std::sqrt(sPerp * sPerp + sZ * sZ);
  }
  struct PrjT { Int_t gi; Int_t n; Int_t cand[3]; };
  struct HitC { Int_t mult; Float_t tot, t; Int_t owner; };
  // hit in cell c of this event A overlaid with partner B (B = 0: A alone); same-cell hits merge
  // (earlier leading edge and its owner, larger ToT, multiplicity of A)
  HitC combHit(const MuMix& A, const MuMix* B, Int_t c)
  {
    HitC h; h.owner = -1; h.mult = 0; h.tot = 0; h.t = 0;
    const Bool_t a = A.totMax[c] > 0, b = B && B->totMax[c] > 0;
    if (a && b) { h.mult = A.mult[c]; h.tot = std::max(A.totMax[c], B->totMax[c]);
                  const Bool_t bF = B->tEarly[c] < A.tEarly[c]; h.owner = bF ? 1 : 0; h.t = bF ? B->tEarly[c] : A.tEarly[c]; }
    else if (a) { h.mult = A.mult[c]; h.tot = A.totMax[c]; h.t = A.tTot[c]; h.owner = 0; }
    else if (b) { h.mult = B->mult[c]; h.tot = B->totMax[c]; h.t = B->tTot[c]; h.owner = 1; }
    return h;
  }
  // StBTofMatchMaker steps C-E on the hit map for all projected tracks
  void prodChoice(const std::vector<PrjT>& P, const MuMix& A, const MuMix* B, std::vector<Int_t>& claims,
                  std::vector<Int_t>& resCell, std::vector<HitC>& resHit, Bool_t bTracks = kFALSE)
  {
    resCell.assign(P.size(), -1); resHit.resize(P.size());
    for (size_t i = 0; i < P.size(); i++)
      for (Int_t q = 0; q < P[i].n; q++) { const HitC h = combHit(A, B, P[i].cand[q]); if (h.mult > 0) claims[P[i].cand[q]] += h.mult; }
    if (bTracks && B) {   // claims of the partner's own tracks on cells with a partner hit (added once per cell)
      for (size_t i = 0; i < P.size(); i++)
        for (Int_t q = 0; q < P[i].n; q++) {
          const Int_t c = P[i].cand[q];
          if (claims[c] > 0 && claims[c] < 100000 && B->totMax[c] > 0 && B->own[c] > 0) {
            const HitC h = combHit(A, B, c);
            claims[c] = 100000 + claims[c] + B->own[c] * h.mult;   // mark as done
          }
        }
      for (size_t i = 0; i < P.size(); i++)
        for (Int_t q = 0; q < P[i].n; q++) { const Int_t c = P[i].cand[q]; if (claims[c] >= 100000) claims[c] -= 100000; }
    }
    for (size_t i = 0; i < P.size(); i++) {
      Int_t nc = 0, single = -1, best = -1; Float_t bestTot = 0; HitC hs, hb; hs.owner = hb.owner = -1;
      for (Int_t q = 0; q < P[i].n; q++) {
        const Int_t c = P[i].cand[q];
        const HitC h = combHit(A, B, c);
        if (h.mult <= 0 || claims[c] != 1) continue;
        nc++; single = c; hs = h;
        if (h.tot < 40. && h.tot > bestTot) { bestTot = h.tot; best = c; hb = h; }
      }
      if (nc == 1) { resCell[i] = single; resHit[i] = hs; }
      else if (nc > 1 && best >= 0) { resCell[i] = best; resHit[i] = hb; }
    }
    for (size_t i = 0; i < P.size(); i++) for (Int_t q = 0; q < P[i].n; q++) claims[P[i].cand[q]] = 0;
  }
  Double_t m2of(Double_t p, Double_t L, Double_t t) { const Double_t b = L / (kCm * t); return p * p * (1. / (b * b) - 1.); }
  Int_t cellOf(const StMuBTofHit* h) { return ((h->tray() - 1) * 32 + (h->module() - 1)) * 6 + (h->cell() - 1); }
  Bool_t prodTof(const StMuDst* d, const StMuBTofHit* h, Double_t& tof)
  {
    const Int_t ig = h->index2Global();
    if (ig < 0) return kFALSE;
    const StMuTrack* t = d->globalTracks(ig);
    if (!t) return kFALSE;
    const StMuBTofPidTraits& pt = t->btofPidTraits();
    if (pt.matchFlag() <= 0 || pt.timeOfFlight() <= 0) return kFALSE;
    tof = pt.timeOfFlight();
    return kTRUE;
  }
}

StTofMEMuMaker::StTofMEMuMaker(const char* name, StMuDstMaker* mu, const char* outName)
  : StTofMEMaker(name, 0, outName), mMu(mu), mOffPattern("output/mudst/off_%d.root"), mHaveOff(kFALSE) {}
StTofMEMuMaker::~StTofMEMuMaker() {}

Bool_t StTofMEMuMaker::loadOffsets(Int_t run)
{
  for (Int_t c = 0; c < kNCellM; c++) gOff[c] = 0;
  TFile* f = TFile::Open(Form(mOffPattern.Data(), run));
  TH2* r = f ? (TH2*)f->Get("hRes1") : 0;
  if (!r) { LOG_ERROR << "StTofMEMuMaker: no offsets for run " << run << endm; if (f) f->Close(); return kFALSE; }
  std::vector<Int_t> ok(kNCellM, 0);
  for (Int_t c = 0; c < kNCellM; c++) {
    TH1D* py = r->ProjectionY("py_off", c + 1, c + 1);
    if (py->GetEntries() >= 10) { Double_t q = 0.5, med = 0; py->GetQuantiles(1, &med, &q); gOff[c] = med; ok[c] = 1; }
    delete py;
  }
  f->Close();
  Int_t nMod = 0, nTray = 0, nOk = 0;
  for (Int_t c = 0; c < kNCellM; c++) nOk += ok[c];
  for (Int_t c = 0; c < kNCellM; c++) {
    if (ok[c]) continue;
    Double_t s = 0; Int_t n = 0;
    const Int_t m0 = (c / 6) * 6;
    for (Int_t k = m0; k < m0 + 6; k++) if (ok[k]) { s += gOff[k]; n++; }
    if (n) { gOff[c] = s / n; nMod++; continue; }
    const Int_t t0 = (c / 192) * 192;
    for (Int_t k = t0; k < t0 + 192; k++) if (ok[k]) { s += gOff[k]; n++; }
    if (n) { gOff[c] = s / n; nTray++; }
  }
  LOG_INFO << "StTofMEMuMaker: run " << run << " cell offsets: " << nOk << " measured, " << nMod << " module avg, "
           << nTray << " tray avg" << endm;
  return kTRUE;
}

Int_t StTofMEMuMaker::Init()
{
  mHaveGeom = loadGeom();
  mFile = new TFile(mOutName.Data(), "RECREATE");
  hEvtFlow = new TH1D("hEvtFlow", "0 all,1 trigger,2 |Vz|,3 Vr,4 offsets,5 T0 (>=3 matched),6 cent ok;;events", 7, -0.5, 6.5);
  hT0nMatch = new TH1D("hT0nMatch", ";matched hits used for T0;events", 200, 0, 400);
  hCalRes = new TH1D("hCalRes", "matched hits: t - tof_{prod};ns;hits", 400, -5, 5);
  hNHitEv = new TH1D("hNHitEv", ";BTOF hits per event;events", 300, 0, 1500);
  hTotAll = new TH1D("hTotAll", ";ToT (ns);hits", 300, 0, 60);
  hCent9 = new TH1D("hCent9", ";cent9;events", 9, -0.5, 8.5);
  hNPrj = new TH2D("hNPrj", "projected tracks;cent9;p (GeV/c)", 9, -0.5, 8.5, kNPm, 0, kPMaxM);
  hNTr = new TH2D("hNTr", "ME trials;cent9;p (GeV/c)", 9, -0.5, 8.5, kNPm, 0, kPMaxM);
  const char* ax = ";cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})";
  hM2SE = new TH3F("hM2SE", Form("SE production #beta%s", ax), 9, -0.5, 8.5, kNPm, 0, kPMaxM, kNM2m, kM2LoM, kM2HiM);
  hM2SEemu = new TH3F("hM2SEemu", Form("SE emulated choice (largest ToT)%s", ax), 9, -0.5, 8.5, kNPm, 0, kPMaxM, kNM2m, kM2LoM, kM2HiM);
  hM2MEtot = new TH3F("hM2MEtot", Form("ME overlay, largest ToT wins%s", ax), 9, -0.5, 8.5, kNPm, 0, kPMaxM, kNM2m, kM2LoM, kM2HiM);
  hM2MEtime = new TH3F("hM2MEtime", Form("ME overlay, earliest time wins%s", ax), 9, -0.5, 8.5, kNPm, 0, kPMaxM, kNM2m, kM2LoM, kM2HiM);
  hM2MEplain = new TH3F("hM2MEplain", Form("ME partner hits only, largest ToT%s", ax), 9, -0.5, 8.5, kNPm, 0, kPMaxM, kNM2m, kM2LoM, kM2HiM);
  hM2MEmerge = new TH3F("hM2MEmerge", Form("ME overlay, same-cell merge, largest ToT wins%s", ax), 9, -0.5, 8.5, kNPm, 0, kPMaxM, kNM2m, kM2LoM, kM2HiM);
  hMEmergeSt = new TH3F("hMEmergeSt", "overlay (merge): 0 this event, 1 partner, 2 no hit;cent9;p;status", 9, -0.5, 8.5, kNPm, 0, kPMaxM, 3, -0.5, 2.5);
  hM2SEprodEmu = new TH3F("hM2SEprodEmu", Form("SE full production emulation%s", ax), 9, -0.5, 8.5, kNPm, 0, kPMaxM, kNM2m, kM2LoM, kM2HiM);
  hM2MEprod = new TH3F("hM2MEprod", Form("ME overlay, full production emulation, partner hit chosen%s", ax), 9, -0.5, 8.5, kNPm, 0, kPMaxM, kNM2m, kM2LoM, kM2HiM);
  hMEprodSt = new TH3F("hMEprodSt", "full emulation overlay: 0 this-event hit, 1 partner hit, 2 none;cent9;p;status", 9, -0.5, 8.5, kNPm, 0, kPMaxM, 3, -0.5, 2.5);
  hChoiceP = new TH2D("hChoiceP", "full emulation: 0 same hit as prod, 1 other cell, 2 we none, 3 prod none;p (GeV/c);code", kNPm, 0, kPMaxM, 4, -0.5, 3.5);
  hM2MEprodB = new TH3F("hM2MEprodB", Form("ME overlay, production emulation incl. partner track claims%s", ax), 9, -0.5, 8.5, kNPm, 0, kPMaxM, kNM2m, kM2LoM, kM2HiM);
  hMEprodBSt = new TH3F("hMEprodBSt", "incl. partner track claims: 0 this-event hit, 1 partner hit, 2 none;cent9;p;status", 9, -0.5, 8.5, kNPm, 0, kPMaxM, 3, -0.5, 2.5);
  hNGlobProj = new TH1D("hNGlobProj", ";projected production-valid global tracks;events", 300, 0, 1500);
  hMEtotSt = new TH3F("hMEtotSt", "overlay (ToT): 0 this event, 1 partner, 2 no hit;cent9;p;status", 9, -0.5, 8.5, kNPm, 0, kPMaxM, 3, -0.5, 2.5);
  hChoice = new TH2D("hChoice", "0 same hit as prod, 1 other cell, 2 we none, 3 prod none;p (GeV/c);code", kNPm, 0, kPMaxM, 4, -0.5, 3.5);
  hDtSEprod = new TH2D("hDtSEprod", "SE production;p (GeV/c);t - t_{#pi} (ns)", 32, 0, 3.2, 600, -3, 27);
  hDtMEtot = new TH2D("hDtMEtot", "ME overlay ToT;p (GeV/c);t - t_{#pi} (ns)", 32, 0, 3.2, 600, -3, 27);
  mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  gCur.tTot.assign(kNCellM, 0.f); gCur.totMax.assign(kNCellM, 0.f); gCur.tEarly.assign(kNCellM, 0.f); gCur.mult.assign(kNCellM, 0); gCur.own.assign(kNCellM, 0);
  LOG_INFO << "StTofMEMuMaker::Init geometry=" << (mHaveGeom ? "yes" : "NO") << " nMix=" << mNMix << endm;
  return mHaveGeom ? kStOK : kStFatal;
}

Int_t StTofMEMuMaker::Make()
{
  StMuDst* d = mMu->muDst();
  if (!d) return kStWarn;
  StMuEvent* ev = d->event();
  if (!ev) return kStWarn;
  hEvtFlow->Fill(0);
  Bool_t trig = mTrigIds.empty();
  for (size_t i = 0; i < mTrigIds.size() && !trig; i++) if (ev->triggerIdCollection().nominal().isTrigger(mTrigIds[i])) trig = kTRUE;
  if (!trig) return kStOK;
  hEvtFlow->Fill(1);
  const StThreeVectorF v = ev->primaryVertexPosition();
  const TVector3 pv(v.x(), v.y(), v.z());
  if (std::fabs(pv.Z()) >= mVzCut) return kStOK;
  hEvtFlow->Fill(2);
  if (pv.Perp() >= mVrCut) return kStOK;
  hEvtFlow->Fill(3);

  const Int_t run = ev->runId();
  if (run != mCurrentRun) {
    for (int i = 0; i < kNVzM; i++) for (int j = 0; j < kNMultM; j++) gBufM[i][j].clear();
    mCurrentRun = run;
    mHaveOff = loadOffsets(run);
    mRefMultCorr->init(run);
    mBadRun = mRefMultCorr->isBadRun(run);
  }
  if (!mHaveOff) return kStOK;
  hEvtFlow->Fill(4);

  // ---- hits: calibrated times, event T0
  const Int_t nh = d->numberOfBTofHit();
  std::vector<Double_t> buf; buf.reserve(nh);
  for (Int_t i = 0; i < nh; i++) {
    StMuBTofHit* h = d->btofHit(i);
    const Int_t c = cellOf(h);
    Double_t tof;
    if (c < 0 || c >= kNCellM || !prodTof(d, h, tof)) continue;
    buf.push_back(h->leadingEdgeTime() - gOff[c] - tof);
  }
  if (buf.size() < 3) return kStOK;
  hEvtFlow->Fill(5);
  hT0nMatch->Fill(buf.size());
  const Double_t T0 = TMath::Median(buf.size(), &buf[0]);
  std::fill(gCur.tTot.begin(), gCur.tTot.end(), 0.f);
  std::fill(gCur.totMax.begin(), gCur.totMax.end(), 0.f);
  std::fill(gCur.tEarly.begin(), gCur.tEarly.end(), 0.f);
  std::fill(gCur.mult.begin(), gCur.mult.end(), 0);
  std::fill(gCur.own.begin(), gCur.own.end(), 0);
  for (Int_t i = 0; i < nh; i++) {
    StMuBTofHit* h = d->btofHit(i);
    const Int_t c = cellOf(h);
    if (c < 0 || c >= kNCellM) continue;
    Float_t t = h->leadingEdgeTime() - gOff[c] - T0;
    if (t < 0.01f) t = 0.01f;                       // keep the hit (0 means empty)
    const Float_t tot = h->tot();
    hTotAll->Fill(tot);
    Double_t tof;
    if (prodTof(d, h, tof)) hCalRes->Fill(t - tof);
    if (tot > gCur.totMax[c]) { gCur.totMax[c] = tot; gCur.tTot[c] = t; }
    if (gCur.tEarly[c] == 0.f || t < gCur.tEarly[c]) gCur.tEarly[c] = t;
    if (gCur.mult[c] < 255) gCur.mult[c]++;
  }
  gCur.nHits = nh;
  hNHitEv->Fill(nh);

  // ---- centrality
  const Int_t refMult = ev->refMult();
  Int_t nBTOFMatch = 0;
  const Int_t nPrim = d->numberOfPrimaryTracks();
  for (Int_t i = 0; i < nPrim; i++) { const StMuTrack* t = d->primaryTracks(i); if (t && t->btofPidTraits().matchFlag() > 0) nBTOFMatch++; }
  mRefMultCorr->initEvent(refMult, pv.Z(), ev->runInfo().zdcCoincidenceRate());
  const Int_t cent9 = mRefMultCorr->getCentralityBin9();
  const Bool_t useM2 = !mBadRun && !mRefMultCorr->isPileUpEvent(refMult, nBTOFMatch, pv.Z()) && cent9 >= 0 && cent9 <= 8;
  Int_t iVz = (Int_t)((pv.Z() + mVzCut) / (2 * mVzCut) * kNVzM);
  iVz = std::min(kNVzM - 1, std::max(0, iVz));
  std::deque<MuMix>& bufM = gBufM[iVz][multClass(refMult)];

  // projections of all production-valid global tracks (every event: also needed when it is a partner)
  const Double_t Bf = ev->magneticField();
  // ---- full production emulation: project all production-valid global tracks
  const Int_t nG = d->numberOfGlobalTracks();
  std::vector<PrjT> P; P.reserve(nG);
  std::vector<Int_t> prjOfGlobal(nG, -1);
  for (Int_t gi = 0; gi < nG; gi++) {
    const StMuTrack* g = d->globalTracks(gi);
    if (!g || g->flag() <= 0 || g->flag() >= 1000) continue;
    const Int_t nF = g->nHitsFit(kTpcId);
    if (nF < 15) continue;
    const float rr = (float)nF / (1.0 * g->nHitsPoss(kTpcId));
    if (rr < 0.52) continue;
    const StPhysicalHelixD oh = g->outerHelix();
    const StThreeVectorD gm = oh.momentum(Bf * kilogauss), go = oh.origin();
    StPicoPhysicalHelix hx(TVector3(gm.x(), gm.y(), gm.z()), TVector3(go.x(), go.y(), go.z()), Bf * kilogauss,
                           (Double_t)oh.charge(Bf * kilogauss));
    Int_t cell = 0, status = 0; Double_t yl = 0, zl = 0; TVector3 pos;
    const Int_t mi = project(hx, cell, yl, zl, pos, status);
    if (mi < 0) continue;
    const TVector3 U(mGeo[mi][3], mGeo[mi][4], mGeo[mi][5]), W(mGeo[mi][9], mGeo[mi][10], mGeo[mi][11]);
    PrjT pj; pj.gi = gi;
    pj.n = candOrder(mi * 6 + (cell - 1), cell, yl, W.Dot(U) / U.Dot(U), kTRUE, pj.cand);
    prjOfGlobal[gi] = P.size();
    P.push_back(pj);
  }
  hNGlobProj->Fill(P.size());
  for (size_t i = 0; i < P.size(); i++) for (Int_t q = 0; q < P[i].n; q++) if (gCur.own[P[i].cand[q]] < 255) gCur.own[P[i].cand[q]]++;
  if (useM2) {
    hEvtFlow->Fill(6); hCent9->Fill(cent9);
    // production hit cell of each global track
    std::vector<Int_t> hitOfGlobal(d->numberOfGlobalTracks(), -1);
    for (Int_t i = 0; i < nh; i++) {
      StMuBTofHit* h = d->btofHit(i);
      const Int_t ig = h->index2Global();
      if (ig >= 0 && ig < (Int_t)hitOfGlobal.size()) hitOfGlobal[ig] = cellOf(h);
    }
    static std::vector<Int_t> claims(kNCellM, 0);
    std::vector<Int_t> seCell; std::vector<HitC> seHit;
    prodChoice(P, gCur, 0, claims, seCell, seHit);
    std::vector<std::vector<Int_t> > meCell(bufM.size()); std::vector<std::vector<HitC> > meHit(bufM.size());
    for (size_t k = 0; k < bufM.size(); k++) prodChoice(P, gCur, &bufM[k], claims, meCell[k], meHit[k]);
    std::vector<std::vector<Int_t> > meCellB(bufM.size()); std::vector<std::vector<HitC> > meHitB(bufM.size());
    for (size_t k = 0; k < bufM.size(); k++) prodChoice(P, gCur, &bufM[k], claims, meCellB[k], meHitB[k], kTRUE);

    for (Int_t i = 0; i < nPrim; i++) {
      const StMuTrack* t = d->primaryTracks(i);
      if (!t) continue;
      const StThreeVectorF pm = t->p();
      const Double_t pt = pm.perp(), p = pm.mag();
      if (pt < 0.15 || p <= 0) continue;
      const Int_t nFit = t->nHitsFit(), nPoss = t->nHitsPoss();
      const Float_t ratio = nPoss > 0 ? Float_t(nFit) / nPoss : 0.;
      if (nFit < 15 || ratio < 0.52f || t->dcaGlobal().mag() > 3.) continue;
      if (std::fabs(pm.pseudoRapidity()) > 1.0) continue;
      const StMuTrack* g = t->globalTrack();
      if (!g) continue;
      // projection with the outer helix of the global track (as StBTofMatchMaker)
      const StPhysicalHelixD oh = g->outerHelix();
      const StThreeVectorD gm = oh.momentum(Bf * kilogauss), go = oh.origin();
      StPicoPhysicalHelix helix(TVector3(gm.x(), gm.y(), gm.z()), TVector3(go.x(), go.y(), go.z()), Bf * kilogauss,
                                (Double_t)oh.charge(Bf * kilogauss));
      Int_t cell = 0, status = 0; Double_t yl = 0, zl = 0; TVector3 pos;
      const Int_t mi = project(helix, cell, yl, zl, pos, status);
      if (mi < 0) continue;
      const Int_t gcell = mi * 6 + (cell - 1);
      const Double_t L = pathL(pv, pos, helix.curvature());
      const TVector3 U(mGeo[mi][3], mGeo[mi][4], mGeo[mi][5]), W(mGeo[mi][9], mGeo[mi][10], mGeo[mi][11]);
      const Double_t pitch = W.Dot(U) / U.Dot(U);
      Int_t cand[3] = {-1, -1, -1};
      const Int_t nCand = candOrder(gcell, cell, yl, pitch, kTRUE, cand);
      hNPrj->Fill(cent9, p);
      const Double_t tPi = L / kCm * std::sqrt(1. + kMPiM * kMPiM / (p * p));

      // SE: production beta
      const StMuBTofPidTraits& ptr = t->btofPidTraits();
      const Bool_t prodOk = ptr.matchFlag() > 0 && ptr.beta() > 0;
      if (prodOk) {
        const Double_t b = ptr.beta();
        hM2SE->Fill(cent9, p, p * p * (1. / (b * b) - 1.));
        hDtSEprod->Fill(p, ptr.timeOfFlight() - tPi);
      }
      // SE: emulated production choice on this event's hits
      Int_t cA = -1; Float_t totA = 0;
      for (Int_t q = 0; q < nCand; q++) if (gCur.totMax[cand[q]] > totA) { totA = gCur.totMax[cand[q]]; cA = cand[q]; }
      if (cA >= 0) hM2SEemu->Fill(cent9, p, m2of(p, L, gCur.tTot[cA]));
      { const Int_t gi = t->index2Global();
        const Int_t pc = (prodOk && gi >= 0 && gi < (Int_t)hitOfGlobal.size()) ? hitOfGlobal[gi] : -1;
        Int_t code = 3;
        if (prodOk && pc >= 0) code = (cA < 0) ? 2 : (cA == pc ? 0 : 1);
        else if (cA < 0) code = -1;
        if (code >= 0) hChoice->Fill(p, code); }
      const Int_t ip = (t->index2Global() >= 0 && t->index2Global() < nG) ? prjOfGlobal[t->index2Global()] : -1;
      if (ip >= 0) {
        if (seCell[ip] >= 0) hM2SEprodEmu->Fill(cent9, p, m2of(p, L, seHit[ip].t));
        const Int_t gi = t->index2Global();
        const Int_t pc = (prodOk && gi < (Int_t)hitOfGlobal.size()) ? hitOfGlobal[gi] : -1;
        Int_t code = -1;
        if (prodOk && pc >= 0) code = (seCell[ip] < 0) ? 2 : (seCell[ip] == pc ? 0 : 1);
        else if (seCell[ip] >= 0) code = 3;
        if (code >= 0) hChoiceP->Fill(p, code);
      }
      Float_t tEA = 1e9f;
      for (Int_t q = 0; q < nCand; q++) { const Float_t te = gCur.tEarly[cand[q]]; if (te > 0 && te < tEA) tEA = te; }

      // ME
      for (size_t k = 0; k < bufM.size(); k++) {
        const MuMix& B = bufM[k];
        hNTr->Fill(cent9, p);
        if (ip >= 0) {
          const Int_t o = meCell[k][ip] >= 0 ? meHit[k][ip].owner : -1;
          hMEprodSt->Fill(cent9, p, o < 0 ? 2 : o);
          if (o == 1) hM2MEprod->Fill(cent9, p, m2of(p, L, meHit[k][ip].t));
          const Int_t oB = meCellB[k][ip] >= 0 ? meHitB[k][ip].owner : -1;
          hMEprodBSt->Fill(cent9, p, oB < 0 ? 2 : oB);
          if (oB == 1) hM2MEprodB->Fill(cent9, p, m2of(p, L, meHitB[k][ip].t));
        }
        Int_t cB = -1; Float_t totB = 0, tEB = 1e9f;
        for (Int_t q = 0; q < nCand; q++) {
          const Int_t c = cand[q];
          if (B.totMax[c] > totB) { totB = B.totMax[c]; cB = c; }
          if (B.tEarly[c] > 0 && B.tEarly[c] < tEB) tEB = B.tEarly[c];
        }
        if (cB >= 0) hM2MEplain->Fill(cent9, p, m2of(p, L, B.tTot[cB]));
        Int_t st = 2;
        if (cB >= 0 && totB > totA) {
          st = 1;
          hM2MEtot->Fill(cent9, p, m2of(p, L, B.tTot[cB]));
          hDtMEtot->Fill(p, B.tTot[cB] - tPi);
        } else if (cA >= 0) st = 0;
        hMEtotSt->Fill(cent9, p, st);
        if (tEB < tEA) hM2MEtime->Fill(cent9, p, m2of(p, L, tEB));
        { // same-cell merge
          Float_t best = 0, bt = 0; Int_t owner = -1;
          for (Int_t q = 0; q < nCand; q++) {
            const Int_t c = cand[q];
            const Float_t a = gCur.totMax[c], b = B.totMax[c];
            if (a <= 0 && b <= 0) continue;
            Float_t tot, t; Int_t o;
            if (a > 0 && b > 0) { tot = std::max(a, b); const Bool_t bFirst = B.tEarly[c] < gCur.tEarly[c];
                                  o = bFirst ? 1 : 0; t = bFirst ? B.tEarly[c] : gCur.tEarly[c]; }
            else if (a > 0) { tot = a; o = 0; t = gCur.tTot[c]; }
            else { tot = b; o = 1; t = B.tTot[c]; }
            if (tot > best) { best = tot; owner = o; bt = t; }
          }
          hMEmergeSt->Fill(cent9, p, owner < 0 ? 2 : owner);
          if (owner == 1) hM2MEmerge->Fill(cent9, p, m2of(p, L, bt));
        }
      }
    }
  }
  bufM.push_back(gCur);
  while ((Int_t)bufM.size() > mNMix) bufM.pop_front();
  return kStOK;
}

Int_t StTofMEMuMaker::Finish()
{
  LOG_INFO << "StTofMEMuMaker::Finish events " << hEvtFlow->GetBinContent(1) << " with T0 " << hEvtFlow->GetBinContent(6)
           << " m2 events " << hEvtFlow->GetBinContent(7) << " projected " << hNPrj->GetEntries() << " trials " << hNTr->GetEntries()
           << " SE " << hM2SE->GetEntries() << " ME(ToT) " << hM2MEtot->GetEntries() << endm;
  mFile->cd(); mFile->Write(); mFile->Close();
  return kStOK;
}
