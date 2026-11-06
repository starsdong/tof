// totData.C : BTOF ToT of matched hits in data MuDst, for the ToT model of the simulation.
// hTotBg: ToT vs log10(beta*gamma) for identified pi, K, p (QA primary tracks, |m2 - m2_PDG| window, 0.2<p<3)
// hTotM2Lo: ToT vs m2 for 0.2<p<0.4 GeV/c (all matched QA tracks), hTotM2Mid: same for 0.4<p<1
// usage: root4star -b -q 'macros/totData.C(N, "list", "out.root")'
void totData(Int_t nEv = 5000, const char* in = "x.list", const char* out = "totdata.root")
{
  gROOT->Macro("loadMuDst.C");
  StChain* chain = new StChain();
  StMuDstMaker* mu = new StMuDstMaker(0, 0, "", in, "", 100);
  chain->Init();
  TFile* f = new TFile(out, "RECREATE");
  TH2D* hTotBg = new TH2D("hTotBg", "identified #pi,K,p;log_{10}(#beta#gamma);ToT (ns)", 60, -0.6, 1.4, 200, 0, 50);
  TH2D* hTotBgPi = new TH2D("hTotBgPi", "#pi;log_{10}(#beta#gamma);ToT (ns)", 60, -0.6, 1.4, 200, 0, 50);
  TH2D* hTotBgP = new TH2D("hTotBgP", "p;log_{10}(#beta#gamma);ToT (ns)", 60, -0.6, 1.4, 200, 0, 50);
  TH2D* hTotM2Lo = new TH2D("hTotM2Lo", "0.2<p<0.4;m^{2};ToT (ns)", 120, -1, 5, 200, 0, 50);
  TH2D* hTotM2Mid = new TH2D("hTotM2Mid", "0.4<p<1;m^{2};ToT (ns)", 120, -1, 5, 200, 0, 50);
  TH1D* hTotAll = new TH1D("hTotAll", "all BTOF hits;ToT (ns)", 200, 0, 50);
  const double mPDG[3] = {0.13957, 0.49368, 0.93827}, m2PDG[3] = {0.01948, 0.24372, 0.88035}, w[3] = {0.05, 0.06, 0.15};
  Int_t n = std::min((Long64_t)nEv, mu->chain()->GetEntries()), nUsed = 0;
  for (Int_t i = 0; i < n; i++) {
    chain->Clear(); if (chain->Make(i) >= kStEOF) break;
    StMuDst* d = mu->muDst(); if (!d || !d->event()) continue;
    const StThreeVectorF pv = d->event()->primaryVertexPosition();
    if (fabs(pv.z()) > 70 || pv.perp() > 2) continue;
    nUsed++;
    for (Int_t j = 0; j < (Int_t)d->numberOfBTofHit(); j++) hTotAll->Fill(d->btofHit(j)->tot());
    for (Int_t j = 0; j < (Int_t)d->numberOfPrimaryTracks(); j++) {
      StMuTrack* t = d->primaryTracks(j);
      if (!t || t->nHitsFit() < 15 || fabs(t->eta()) > 1 || t->pt() < 0.15 || t->dcaGlobal().mag() > 3) continue;
      const StMuBTofPidTraits& pid = t->btofPidTraits();
      if (pid.matchFlag() <= 0 || pid.beta() <= 0) continue;
      const StMuTrack* g = t->globalTrack(); if (!g || g->index2BTofHit() < 0) continue;
      const StMuBTofHit* h = d->btofHit(g->index2BTofHit()); if (!h) continue;
      const double p = t->p().mag(), b = pid.beta(), m2 = p * p * (1. / (b * b) - 1.), tot = h->tot();
      if (p > 0.2 && p < 0.4) hTotM2Lo->Fill(m2, tot);
      if (p > 0.4 && p < 1.0) hTotM2Mid->Fill(m2, tot);
      if (p < 0.2 || p > 3) continue;
      for (int k = 0; k < 3; k++) if (fabs(m2 - m2PDG[k]) < w[k]) {
        if (k == 1 && p > 1.2) continue;
        if (k == 2 && p > 2.5) continue;
        const double bg = log10(p / mPDG[k]);
        hTotBg->Fill(bg, tot); if (k == 0) hTotBgPi->Fill(bg, tot); if (k == 2) hTotBgP->Fill(bg, tot);
      }
    }
  }
  f->Write(); f->Close();
  printf("totData: done %d events used %d\n", n, nUsed);
}
