// probeMuBTof.C : BTOF hits in MuDst (all hits have leading/trailing edge times and ToT).
// Raw leading-edge times carry per-channel offsets (several ns), so the macro is run twice on the same events:
//  pass 1: per event T0 = median(t_LE - tof) over matched hits; writes hRes1(cell, t_LE - tof - T0) to <out>_pass1.
//  pass 2: reads the per-cell offsets (median of hRes1, cells with >= 10 entries); with c = t_LE - off_cell,
//          T0' = median(c - tof) over matched hits and every hit gets t = c - T0' (a time of flight; ToT slewing
//          not corrected). Matched = hit of a global track with matchFlag>0 and tof>0.
// root4star -b -q 'macros/probeMuBTof.C(-1,"list","output/mudst_probe.root",1)'  then the same with 2.
// Event selection: |Vz| < 70 cm, at least 3 matched hits.
void probeMuBTof(Int_t nEv = -1, const char* list = "filelists/mudst/test1.list", const char* out = "output/mudst_probe.root",
                 Int_t pass = 1)
{
  gROOT->Macro("loadMuDst.C");
  StChain* chain = new StChain();             // must exist before the makers so they attach to it
  StMuDstMaker* mu = new StMuDstMaker(0, 0, "", list, "MuDst", 10);
  mu->SetStatus("*", 0);
  mu->SetStatus("MuEvent*", 1); mu->SetStatus("GlobalTracks*", 1); mu->SetStatus("BTof*", 1); mu->SetStatus("PrimaryVertices*", 1);
  const Int_t kNC = 23040;
  TString p1name = TString(out).ReplaceAll(".root", "_pass1.root");
  static Double_t off[23040];
  static Int_t nOffEnt[23040];
  TH1D* hOff = 0;
  if (pass == 2) {
    TFile* fp = TFile::Open(p1name);
    TH2F* r = (TH2F*)fp->Get("hRes1");
    hOff = new TH1D("hOff", "cell offsets;cell id;offset (ns)", kNC, -0.5, kNC - 0.5);
    hOff->SetDirectory(0);
    for (Int_t c = 0; c < kNC; c++) {
      TH1D* py = r->ProjectionY("py", c + 1, c + 1);
      nOffEnt[c] = (Int_t)py->GetEntries(); off[c] = 0;
      if (nOffEnt[c] >= 10) { Double_t q = 0.5, med = 0; py->GetQuantiles(1, &med, &q); off[c] = med; }
      hOff->SetBinContent(c + 1, off[c]);
      delete py;
    }
    fp->Close();
  } else for (Int_t c = 0; c < kNC; c++) { off[c] = 0; nOffEnt[c] = 0; }

  TFile* f = new TFile(pass == 1 ? p1name.Data() : out, "RECREATE");
  if (hOff) hOff->SetDirectory(f);
  TH2F* hRes1   = new TH2F("hRes1", "pass 1;cell id;t_{LE} - tof - T_{0} (ns)", kNC, -0.5, kNC - 0.5, 300, -15, 15);
  TH1D* hNHit   = new TH1D("hNHit", ";BTOF hits per event;events", 300, 0, 1500);
  TH1D* hFracM  = new TH1D("hFracM", ";fraction of hits matched to a track;events", 100, 0, 1);
  TH1D* hDtM    = new TH1D("hDtM", "matched;t (ns);hits", 3000, -100, 200);
  TH1D* hDtU    = new TH1D("hDtU", "unmatched;t (ns);hits", 3000, -100, 200);
  TH1D* hDtUn   = new TH1D("hDtUn", "unmatched, cell without offset;t (ns);hits", 3000, -100, 200);
  TH1D* hTotM   = new TH1D("hTotM", "matched;ToT (ns);hits", 300, 0, 60);
  TH1D* hTotU   = new TH1D("hTotU", "unmatched;ToT (ns);hits", 300, 0, 60);
  TH2D* hTotDtU = new TH2D("hTotDtU", "unmatched;t (ns);ToT (ns)", 300, -50, 100, 120, 0, 60);
  TH2D* hTotDtM = new TH2D("hTotDtM", "matched;t (ns);ToT (ns)", 300, -50, 100, 120, 0, 60);
  TH1D* hRes2   = new TH1D("hRes2", "pass 2 matched: t - tof;ns;hits", 2000, -10, 10);
  TH2D* hDtUmult = new TH2D("hDtUmult", "unmatched;BTOF hits per event;t (ns)", 30, 0, 1500, 300, -50, 100);
  TH2D* hDtMmult = new TH2D("hDtMmult", "matched;BTOF hits per event;t (ns)", 30, 0, 1500, 300, -50, 100);
  TH2D* hDtUmod = new TH2D("hDtUmod", "unmatched;half*32+module;t (ns)", 64, -0.5, 63.5, 300, -50, 100);
  TH3F* hTU3 = new TH3F("hTU3", "unmatched;Vz bin;half*32+module;t (ns)", 14, -0.5, 13.5, 64, -0.5, 63.5, 1400, -20, 120);
  TH3F* hTM3 = new TH3F("hTM3", "matched;Vz bin;half*32+module;t (ns)", 14, -0.5, 13.5, 64, -0.5, 63.5, 1400, -20, 120);
  TH2D* hDtMmod = new TH2D("hDtMmod", "matched;half*32+module;t (ns)", 64, -0.5, 63.5, 300, -50, 100);

  chain->Init();
  Int_t nTot = mu->chain()->GetEntries();
  Int_t nRun = (nEv > 0 && nEv < nTot) ? nEv : nTot;
  cout << "probeMuBTof: " << nTot << " events in chain, processing " << nRun << ", pass " << pass << endl;
  static Double_t buf[20000];
  for (Int_t iev = 0; iev < nRun; iev++) {
    chain->Clear();
    if (chain->Make(iev) != kStOK) continue;
    StMuDst* d = mu->muDst();
    if (!d || !d->event()) continue;
    const Double_t vz = d->event()->primaryVertexPosition().z();
    if (std::fabs(vz) > 70) continue;
    Int_t iVz = (Int_t)((vz + 70.) / 140. * 14);
    if (iVz < 0) iVz = 0;
    if (iVz > 13) iVz = 13;
    const Int_t nh = d->numberOfBTofHit();
    if (nh <= 0) continue;
    Int_t nm = 0;
    for (Int_t i = 0; i < nh; i++) {
      StMuBTofHit* h = d->btofHit(i);
      const Int_t ig = h->index2Global();
      if (ig < 0) continue;
      StMuTrack* t = d->globalTracks(ig);
      if (!t) continue;
      const StMuBTofPidTraits& pt = t->btofPidTraits();
      if (pt.matchFlag() <= 0 || pt.timeOfFlight() <= 0) continue;
      const Int_t cid = ((h->tray() - 1) * 32 + (h->module() - 1)) * 6 + (h->cell() - 1);
      const Double_t o = (cid >= 0 && cid < kNC) ? off[cid] : 0;
      if (nm < 20000) buf[nm++] = h->leadingEdgeTime() - o - pt.timeOfFlight();
    }
    if (nm < 3) continue;
    const Double_t T = TMath::Median(nm, buf);
    hNHit->Fill(nh); hFracM->Fill(Double_t(nm) / nh);
    for (Int_t i = 0; i < nh; i++) {
      StMuBTofHit* h = d->btofHit(i);
      const Int_t cid = ((h->tray() - 1) * 32 + (h->module() - 1)) * 6 + (h->cell() - 1);
      if (cid < 0 || cid >= kNC) continue;
      Bool_t matched = kFALSE; Double_t tof = 0;
      const Int_t ig = h->index2Global();
      if (ig >= 0) {
        StMuTrack* t = d->globalTracks(ig);
        if (t) { const StMuBTofPidTraits& pt = t->btofPidTraits();
          if (pt.matchFlag() > 0 && pt.timeOfFlight() > 0) { matched = kTRUE; tof = pt.timeOfFlight(); } }
      }
      if (pass == 1) { if (matched) hRes1->Fill(cid, h->leadingEdgeTime() - tof - T); continue; }
      const Double_t tt = h->leadingEdgeTime() - off[cid] - T, tot = h->tot();
      const Int_t key = ((h->tray() <= 60) ? 0 : 32) + (h->module() - 1);
      if (matched) { hDtM->Fill(tt); hTotM->Fill(tot); hTotDtM->Fill(tt, tot); hRes2->Fill(tt - tof);
                     hDtMmult->Fill(nh, tt); hDtMmod->Fill(key, tt); hTM3->Fill(iVz, key, tt); }
      else if (nOffEnt[cid] >= 10) { hDtU->Fill(tt); hTotU->Fill(tot); hTotDtU->Fill(tt, tot); hDtUmult->Fill(nh, tt); hDtUmod->Fill(key, tt); hTU3->Fill(iVz, key, tt); }
      else hDtUn->Fill(tt);
    }
  }
  chain->Finish();
  f->cd(); f->Write(); f->Close();
  cout << "probeMuBTof: done pass " << pass << ", hits matched " << hDtM->GetEntries() << " unmatched " << hDtU->GetEntries()
       << " (cells without offset: " << hDtUn->GetEntries() << ", pass-1 entries " << hRes1->GetEntries() << ")" << endl;
}
