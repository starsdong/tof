#include "StTofSimAnaMaker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StEvent/StBTofHeader.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TClonesArray.h"
#include <cmath>
#include <map>
#include <vector>

ClassImp(StTofSimAnaMaker)

namespace {
  const Double_t kC = 29.9792458, kMPi = 0.13957;
  const Double_t kMcClock = 570. / kC;   // StVpdSimConfig::kMcClock (ns), added to every simulated BTOF time
  const Int_t kNR = 25; const Double_t kRMax = 500.;
  const Int_t kNP = 32; const Double_t kPMax = 3.2;
  const char* kCatName[6] = {"own", "decay", "secondary", "relative", "unrelated", "notruth"};
  std::map<Int_t, const StMuMcTrack*> gTrk;   // MC track id -> track
  std::map<Int_t, const StMuMcVertex*> gVx;   // MC vertex id -> vertex
  Double_t pathL(const StThreeVectorF& b, const StThreeVectorF& e, Double_t curv)
  {
    const Double_t xd = e.x() - b.x(), yd = e.y() - b.y();
    const Double_t C = std::sqrt(xd * xd + yd * yd);
    Double_t sPerp = C;
    if (curv != 0) { const Double_t R = 1. / curv; sPerp = 2 * R * std::asin(C / (2 * R)); }
    const Double_t sZ = std::fabs(e.z() - b.z());
    return std::sqrt(sPerp * sPerp + sZ * sZ);
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
  // parent MC track id of an MC track (0 if it starts at a primary/event vertex)
  Int_t parentOf(Int_t id)
  {
    std::map<Int_t, const StMuMcTrack*>::const_iterator it = gTrk.find(id);
    if (it == gTrk.end()) return 0;
    std::map<Int_t, const StMuMcVertex*>::const_iterator iv = gVx.find(it->second->IdVx());
    if (iv == gVx.end()) return 0;
    return iv->second->IdParTrk();
  }
  std::vector<Int_t> chain(Int_t id)   // id, parent, grandparent, ... (max 30)
  {
    std::vector<Int_t> c;
    for (Int_t k = 0; k < 30 && id > 0; k++) { c.push_back(id); id = parentOf(id); }
    return c;
  }
  Float_t vxR(Int_t trkId)   // radius of the start vertex of an MC track
  {
    std::map<Int_t, const StMuMcTrack*>::const_iterator it = gTrk.find(trkId);
    if (it == gTrk.end()) return -1;
    std::map<Int_t, const StMuMcVertex*>::const_iterator iv = gVx.find(it->second->IdVx());
    if (iv == gVx.end()) return -1;
    return iv->second->XyzV().perp();
  }
  Bool_t decays(Int_t gePid) { return gePid == 5 || gePid == 6 || gePid == 8 || gePid == 9 || gePid == 11 || gePid == 12 || gePid == 10 || gePid == 16 || gePid == 19 || gePid == 21 || gePid == 23 || gePid == 24; }
}

StTofSimAnaMaker::StTofSimAnaMaker(const char* name, StMuDstMaker* mu, const char* outName)
  : StMaker(name), mMu(mu), mOutName(outName), mFile(0) {}

Int_t StTofSimAnaMaker::Init()
{
  mFile = new TFile(mOutName.Data(), "RECREATE");
  hEvt = new TH1D("hEvt", "0 all, 1 |Vz|<70, 2 Vr<2, 3 analysed;;events", 4, -0.5, 3.5);
  hRefMult = new TH1D("hRefMult", ";refMult;events", 500, 0, 500);
  hNTrk = new TH2D("hNTrk", "QA tracks;refMult;p (GeV/c)", kNR, 0, kRMax, kNP, 0, kPMax);
  for (int c = 0; c < kNCat; c++) {
    hQM2[c] = new TH3F(Form("hQM2_%s", kCatName[c]), Form("matched QA tracks, %s;refMult;p (GeV/c);q #times m^{2} (GeV^{2}/c^{4})", kCatName[c]),
                       kNR, 0, kRMax, kNP, 0, kPMax, 500, -5, 5);
    hDt[c] = new TH2D(Form("hDt_%s", kCatName[c]), Form("%s;p (GeV/c);t - t_{#pi} (ns)", kCatName[c]), kNP, 0, kPMax, 600, -3, 27);
    hRVx[c] = new TH2D(Form("hRVx_%s", kCatName[c]), Form("%s;p (GeV/c);R of branch vertex (cm)", kCatName[c]), kNP, 0, kPMax, 300, 0, 600);
    hTot[c] = new TH2D(Form("hTot_%s", kCatName[c]), Form("%s;p (GeV/c);ToT (ns)", kCatName[c]), kNP, 0, kPMax, 200, 0, 50);
    hTotHi[c] = new TH2D(Form("hTotHi_%s", kCatName[c]), Form("%s, 1.6<|q m^{2}|<2.8;p (GeV/c);ToT (ns)", kCatName[c]), kNP, 0, kPMax, 200, 0, 50);
    hPid[c] = new TH2D(Form("hPid_%s", kCatName[c]), Form("%s;p (GeV/c);GEANT id of hit particle", kCatName[c]), kNP, 0, kPMax, 50, -0.5, 49.5);
  }
  hQM2cal = new TH3F("hQM2cal", "matched QA tracks, calibrated tof;refMult;p (GeV/c);q #times m^{2} (GeV^{2}/c^{4})", kNR, 0, kRMax, kNP, 0, kPMax, 500, -5, 5);
  hTStart = new TH2D("hTStart", ";refMult;tStart - McClock (ns)", kNR, 0, kRMax, 1000, -50, 50);
  return kStOK;
}

Int_t StTofSimAnaMaker::classify(StMuDst* d, Int_t trkId, Int_t hitId, Float_t& rVx, Int_t& gePid)
{
  rVx = -1; gePid = 0;
  if (trkId <= 0 || hitId <= 0) return 5;
  std::map<Int_t, const StMuMcTrack*>::const_iterator ih = gTrk.find(hitId);
  if (ih != gTrk.end()) gePid = ih->second->GePid();
  if (hitId == trkId) return 0;
  const std::vector<Int_t> ch = chain(hitId), ct = chain(trkId);
  // descendant: the track id appears in the hit particle's ancestry
  for (size_t k = 1; k < ch.size(); k++) if (ch[k] == trkId) {
    const Int_t child = ch[k - 1];   // first particle of the branch below the track
    rVx = vxR(child);
    std::map<Int_t, const StMuMcTrack*>::const_iterator it = gTrk.find(trkId), ic = gTrk.find(child);
    Bool_t isDecay = kFALSE;
    // GEANT3 ends a track at a hadronic interaction as well, so the end vertex alone does not identify a decay:
    // decay = weakly decaying parent AND (muon daughter OR vertex inside the TPC gas volume, R < 190 cm)
    if (it != gTrk.end() && ic != gTrk.end()) {
      const Int_t gd = ic->second->GePid();
      isDecay = decays(it->second->GePid()) && it->second->IdVxEnd() == ic->second->IdVx() && (gd == 5 || gd == 6 || (rVx >= 0 && rVx < 190.));
    }
    return isDecay ? 1 : 2;
  }
  // relative: hit particle is an ancestor of the track, or both share a non-primary ancestor
  for (size_t a = 1; a < ct.size(); a++) for (size_t b = 0; b < ch.size(); b++) if (ct[a] == ch[b]) { rVx = vxR(ch[b > 0 ? b - 1 : 0]); return 3; }
  return 4;
}

Int_t StTofSimAnaMaker::Make()
{
  StMuDst* d = mMu->muDst();
  if (!d || !d->event()) return kStOK;
  hEvt->Fill(0);
  StMuEvent* ev = d->event();
  const StThreeVectorF pv = ev->primaryVertexPosition();
  if (std::fabs(pv.z()) > 70.) return kStOK;
  hEvt->Fill(1);
  if (pv.perp() > 2.) return kStOK;
  hEvt->Fill(2);
  gTrk.clear(); gVx.clear();
  TClonesArray* mcT = d->mcArray(1); TClonesArray* mcV = d->mcArray(0);
  if (mcT) for (int i = 0; i < mcT->GetEntriesFast(); i++) { const StMuMcTrack* t = (const StMuMcTrack*)mcT->UncheckedAt(i); gTrk[t->Id()] = t; }
  if (mcV) for (int i = 0; i < mcV->GetEntriesFast(); i++) { const StMuMcVertex* v = (const StMuMcVertex*)mcV->UncheckedAt(i); gVx[v->Id()] = v; }
  const Int_t refMult = ev->refMult();
  hRefMult->Fill(refMult);
  if (d->btofHeader()) hTStart->Fill(refMult, d->btofHeader()->tStart() - kMcClock);
  hEvt->Fill(3);
  for (Int_t i = 0; i < (Int_t)d->numberOfPrimaryTracks(); i++) {
    const StMuTrack* t = d->primaryTracks(i);
    if (!qaTrack(t)) continue;
    const Double_t p = t->p().mag();
    hNTrk->Fill(refMult, p);
    const StMuTrack* g = t->globalTrack();
    const StMuBTofPidTraits& pid = g->btofPidTraits();
    if (pid.matchFlag() <= 0 || g->index2BTofHit() < 0) continue;
    const StMuBTofHit* h = d->btofHit(g->index2BTofHit());
    if (!h) continue;
    const Double_t L = pathL(pv, pid.position(), t->helix().curvature());
    const Double_t tTrue = h->leadingEdgeTime() - kMcClock;   // simulated hit time relative to the collision
    if (tTrue <= 0 || L <= 0) continue;
    const Double_t beta = L / (kC * tTrue);
    const Double_t m2 = p * p * (1. / (beta * beta) - 1.);
    const Double_t dt = tTrue - L / kC * std::sqrt(1. + kMPi * kMPi / (p * p));
    Float_t rVx; Int_t gePid;
    const Int_t cat = classify(d, g->idTruth(), h->idTruth(), rVx, gePid);
    hQM2[cat]->Fill(refMult, p, t->charge() * m2);
    hDt[cat]->Fill(p, dt);
    hRVx[cat]->Fill(p, rVx);
    hPid[cat]->Fill(p, gePid);
    hTot[cat]->Fill(p, h->tot());
    if (std::fabs(m2) > 1.6 && std::fabs(m2) < 2.8) hTotHi[cat]->Fill(p, h->tot());
    if (pid.timeOfFlight() > 0) {
      const Double_t bc = L / (kC * pid.timeOfFlight());
      hQM2cal->Fill(refMult, p, t->charge() * p * p * (1. / (bc * bc) - 1.));
    }
  }
  return kStOK;
}

Int_t StTofSimAnaMaker::Finish()
{
  LOG_INFO << "StTofSimAnaMaker::Finish events " << hEvt->GetBinContent(1) << " analysed " << hEvt->GetBinContent(4) << " QA tracks " << hNTrk->GetEntries() << endm;
  mFile->Write();
  mFile->Close();
  return kStOK;
}
