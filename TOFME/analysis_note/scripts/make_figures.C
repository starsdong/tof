// make_figures.C : note figures for the TOFME QA (uses the user's style.C/draw.C via ~/.rootlogon.C)
// usage (from analysis_note/): root -b -q -l 'scripts/make_figures.C("data/qa_test_10k.root")'
// writes figures/<name>.pdf and figures/<name>.png (overwritten in place)

static TFile* gF = 0;
static TString gTag;

static void label(double x, double y, const char* t, double size = 0.045, int color = 1)
{
  TLatex* l = new TLatex(x, y, t);
  l->SetNDC(); l->SetTextFont(42); l->SetTextSize(size); l->SetTextColor(color);
  l->Draw();
}
static TCanvas* canvas(const char* name, bool is2D = false)
{
  TCanvas* c = new TCanvas(name, "", 800, 600);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(is2D ? 0.15 : 0.04);
  c->SetTopMargin(0.09);
  c->SetBottomMargin(0.15);
  return c;
}
static void save(TCanvas* c, const char* name)
{
  c->SaveAs(Form("figures/%s.pdf", name));
  c->SaveAs(Form("figures/%s.png", name));
}
static TH1* H1(const char* n) { TH1* h = (TH1*)gF->Get(n); h->SetTitle(""); h->SetStats(0); return h; }
static TH2* H2(const char* n) { TH2* h = (TH2*)gF->Get(n); h->SetTitle(""); h->SetStats(0); return h; }
static void tag(double x = 0.15, double y = 0.935) { label(x, y, gTag, 0.042); }

static void draw1(const char* hn, const char* fn, bool logy, double cutLo = -1e9, double cutHi = -1e9,
                  double xmin = 1, double xmax = 0)
{
  TCanvas* c = canvas(fn);
  c->SetLogy(logy);
  TH1* h = H1(hn);
  if (xmin < xmax) h->GetXaxis()->SetRangeUser(xmin, xmax);
  h->SetLineWidth(2); h->SetLineColor(kBlue + 2);
  h->SetMaximum(h->GetMaximum() * (logy ? 30 : 1.35));
  if (logy) h->SetMinimum(0.5);
  h->Draw("hist");
  double y1 = logy ? 0.5 : 0, y2 = h->GetMaximum() / (logy ? 3 : 1.15);
  if (cutLo > -1e8) drawLine(cutLo, y1, cutLo, y2, 2, 2, kRed);
  if (cutHi > -1e8) drawLine(cutHi, y1, cutHi, y2, 2, 2, kRed);
  tag();
  save(c, fn);
}
static TCanvas* draw2(const char* hn, const char* fn, bool logz = true)
{
  TCanvas* c = canvas(fn, true);
  c->SetLogz(logz);
  TH2* h = H2(hn);
  h->GetZaxis()->SetLabelSize(0.04);
  h->Draw("colz");
  tag();
  return c;
}

void make_figures(const char* fn = "data/qa_test_10k.root")
{
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kBird);
  gStyle->SetTitleSize(0.06, "XYZ");
  gStyle->SetTitleOffset(1.05, "X");
  gStyle->SetTitleOffset(1.15, "Y");
  gStyle->SetLabelSize(0.05, "XY");
  gF = TFile::Open(fn);
  const double nSel = ((TH1*)gF->Get("hEvtCut"))->GetBinContent(4);
  gTag = Form("Au+Au 19.6 GeV, Run19 P24iy, %.0f events", nSel);

  // ---------------- event level
  {
    TCanvas* c = canvas("evt_vz"); c->SetLogy(0);
    TH1* a = H1("hVzAll"); TH1* s = H1("hVz");
    a->SetLineColor(kGray + 2); a->SetLineWidth(2); s->SetLineColor(kBlue + 2); s->SetLineWidth(2);
    s->SetFillColorAlpha(kBlue - 9, 0.5);
    a->SetMaximum(a->GetMaximum() * 1.4); a->SetMinimum(0);
    a->GetYaxis()->SetTitle("N_{events} / 1 cm");
    a->Draw("hist"); s->Draw("hist same");
    drawLine(-70, 0, -70, a->GetMaximum() / 1.4, 2, 2, kRed);
    drawLine(70, 0, 70, a->GetMaximum() / 1.4, 2, 2, kRed);
    TLegend* l = new TLegend(0.18, 0.72, 0.6, 0.88); l->SetBorderSize(0); l->SetTextSize(0.04);
    l->AddEntry(a, "all events (before cuts)", "l"); l->AddEntry(s, "selected (trigger, V_{z}, V_{r})", "f");
    l->Draw();
    label(0.15, 0.935, "Au+Au 19.6 GeV, Run19 P24iy, 10000 events read", 0.042);
    save(c, "evt_vz");
  }
  {
    TCanvas* c = draw2("hVxVy", "evt_vxvy");
    H2("hVxVy")->GetXaxis()->SetRangeUser(-3, 3); H2("hVxVy")->GetYaxis()->SetRangeUser(-3, 3);
    TEllipse* e = new TEllipse(0, 0, 2, 2); e->SetFillStyle(0); e->SetLineColor(kRed); e->SetLineStyle(2); e->SetLineWidth(2);
    e->Draw();
    label(0.18, 0.2, "V_{r} < 2 cm (dashed); |V_{z}|<70 cm applied", 0.04, kRed);
    save(c, "evt_vxvy");
  }
  draw1("hVzDiff", "evt_vzdiff", true);
  {
    TCanvas* c = canvas("evt_refmult"); c->SetLogy(1);
    TH1* r = H1("hRefMult"); TH1* g = H1("hGRefMult");
    r->SetLineColor(kBlue + 2); r->SetLineWidth(2); g->SetLineColor(kRed + 1); g->SetLineWidth(2); g->SetLineStyle(2);
    r->GetXaxis()->SetTitle("multiplicity"); r->GetXaxis()->SetRangeUser(0, 450);
    r->SetMaximum(r->GetMaximum() * 30); r->SetMinimum(0.5);
    r->Draw("hist"); g->Draw("hist same");
    TLegend* l = new TLegend(0.6, 0.62, 0.92, 0.8); l->SetBorderSize(0); l->SetTextSize(0.045);
    l->AddEntry(r, "refMult", "l"); l->AddEntry(g, "gRefMult", "l"); l->Draw();
    tag(); save(c, "evt_refmult");
  }
  { TCanvas* c = draw2("hRefMultVsTofMult", "evt_refmult_vs_tofmult"); save(c, "evt_refmult_vs_tofmult"); }
  { TCanvas* c = draw2("hRefMultVsNBTofMatch", "evt_refmult_vs_nmatch"); save(c, "evt_refmult_vs_nmatch"); }
  { TCanvas* c = draw2("hNBTofHitsVsNBTofMatch", "evt_nhits_vs_nmatch"); save(c, "evt_nhits_vs_nmatch"); }
  { TCanvas* c = draw2("hNTofTrkVsNBTofHits", "evt_ntoftrk_vs_nhits"); save(c, "evt_ntoftrk_vs_nhits"); }
  draw1("hNBTofHits", "evt_nbtofhits", true, -1e9, -1e9, 0, 1300);

  // ---------------- track level
  draw1("hNHitsFit", "trk_nhitsfit", false, 15);
  draw1("hNHitsRatio", "trk_nhitsratio", false, 0.52);
  draw1("hDca", "trk_dca", true, 3);
  draw1("hEta", "trk_eta", false, -1, 1, -1.8, 1.8);
  draw1("hPt", "trk_pt", true);
  draw1("hPhi", "trk_phi", false);
  draw1("hNHitsDedx", "trk_nhitsdedx", false);
  { TCanvas* c = draw2("hEtaPhi", "trk_etaphi", false); H2("hEtaPhi")->GetXaxis()->SetRangeUser(-1.2, 1.2); save(c, "trk_etaphi"); }
  { TCanvas* c = draw2("hDedxVsPq", "trk_dedx"); H2("hDedxVsPq")->GetXaxis()->SetRangeUser(-3, 3); save(c, "trk_dedx"); }
  { TCanvas* c = draw2("hNSigPiVsP", "trk_nsigpi"); save(c, "trk_nsigpi"); }

  // ---------------- TOF matched
  {
    TCanvas* c = canvas("tof_matchflag"); TH1* h = H1("hMatchFlag");
    h->SetLineWidth(2); h->SetLineColor(kBlue + 2); h->SetMinimum(0); h->SetMaximum(h->GetMaximum() * 1.35);
    h->GetXaxis()->SetNdivisions(5); h->Draw("hist"); tag(); save(c, "tof_matchflag");
  }
  draw1("hBeta", "tof_beta", true, 1.0, -1e9, 0, 1.4);
  {
    TCanvas* c = draw2("hInvBetaVsP", "tof_invbeta");
    const char* nm[] = {"#pi", "K", "p"}; double m[] = {0.13957, 0.49368, 0.93827};
    for (int i = 0; i < 3; i++) {
      TF1* f = new TF1(Form("ib%d", i), Form("sqrt(x*x+%f)/x", m[i] * m[i]), 0.1, 4);
      f->SetLineColor(kRed); f->SetLineStyle(2); f->SetLineWidth(3); f->SetNpx(400); f->Draw("same");
    }
    drawText(1.6, 1.08, "#pi", 42, 0.055, 0, kRed); drawText(0.36, 1.62, "K", 42, 0.055, 0, kRed); drawText(0.62, 1.95, "p", 42, 0.055, 0, kRed);
    save(c, "tof_invbeta");
  }
  { TCanvas* c = draw2("hInvBetaVsPq", "tof_invbeta_pq"); save(c, "tof_invbeta_pq"); }
  { TCanvas* c = draw2("hM2VsPq", "tof_m2_pq"); H2("hM2VsPq")->GetXaxis()->SetRangeUser(-3, 3); H2("hM2VsPq")->GetYaxis()->SetRangeUser(-0.3, 4.2); save(c, "tof_m2_pq"); }
  { TCanvas* c = draw2("hM2VsP", "tof_m2_p"); H2("hM2VsP")->GetXaxis()->SetRangeUser(0, 3); H2("hM2VsP")->GetYaxis()->SetRangeUser(-0.3, 4.2); save(c, "tof_m2_p"); }
  {
    TCanvas* c = canvas("tof_m2_proj"); c->SetLogy(1);
    TH2* hm = H2("hM2VsP");
    hm->GetXaxis()->SetRange(0, 0); hm->GetYaxis()->SetRange(0, 0);
    TH1* pm = hm->ProjectionY("pm_fig", hm->GetXaxis()->FindBin(0.4 + 1e-4), hm->GetXaxis()->FindBin(1.0 - 1e-4));
    pm->SetTitle(""); pm->GetYaxis()->SetTitle("N_{tracks} / (0.01 GeV^{2}/c^{4})"); pm->GetYaxis()->SetTitleSize(0.05); pm->GetYaxis()->SetTitleOffset(1.4);
    pm->SetMarkerStyle(20); pm->SetMarkerSize(0.6); pm->SetLineColor(1);
    pm->SetMaximum(pm->GetMaximum() * 30); pm->SetMinimum(0.5);
    pm->Draw("pe");
    double lo[] = {-0.06, 0.17, 0.75, 3.2}, hi[] = {0.10, 0.32, 1.0, 3.8};
    for (int i = 0; i < 4; i++) {
      TF1* g = new TF1(Form("gfig%d", i), "gaus", lo[i], hi[i]);
      g->SetLineColor(kRed); g->SetLineWidth(2);
      pm->Fit(g, "RQ+");
    }
    tag();
    label(0.45, 0.83, "0.4 < p < 1.0 GeV/c, Gaussian fits (red)", 0.04);
    save(c, "tof_m2_proj");
  }
  draw1("hYLocal", "tof_ylocal", true);
  draw1("hZLocal", "tof_zlocal", true);
  { TCanvas* c = draw2("hYLocalVsZLocal", "tof_ylz"); save(c, "tof_ylz"); }
  { TCanvas* c = draw2("hYLocalVsP", "tof_ylocal_p"); save(c, "tof_ylocal_p"); }
  { TCanvas* c = draw2("hZLocalVsP", "tof_zlocal_p"); save(c, "tof_zlocal_p"); }
  { TCanvas* c = draw2("hTrayModHits", "tof_occ_hits", false); save(c, "tof_occ_hits"); }
  { TCanvas* c = draw2("hTrayModMatched", "tof_occ_matched", false); save(c, "tof_occ_matched"); }
  { TCanvas* c = draw2("hHitPosPhiZ", "tof_hitpos"); save(c, "tof_hitpos"); }
  { TCanvas* c = draw2("hNSigTofPiVsP", "tof_nsig_pi"); save(c, "tof_nsig_pi"); }
  { TCanvas* c = draw2("hNSigTofKVsP", "tof_nsig_k"); save(c, "tof_nsig_k"); }
  { TCanvas* c = draw2("hNSigTofPVsP", "tof_nsig_p"); save(c, "tof_nsig_p"); }
  { TCanvas* c = draw2("hDInvBetaPiVsP", "tof_dinvbeta_pi"); save(c, "tof_dinvbeta_pi"); }

  // efficiency
  const char* num[] = {"hPtTof", "hEtaTof", "hPhiTof"};
  const char* den[] = {"hPtAll", "hEtaAll", "hPhiAll"};
  const char* out[] = {"tof_eff_pt", "tof_eff_eta", "tof_eff_phi"};
  for (int i = 0; i < 3; i++) {
    TH1* n = H1(num[i]); TH1* d = H1(den[i]);
    TGraphAsymmErrors* g = new TGraphAsymmErrors(n, d, "cl=0.683 b(1,1) mode");
    TCanvas* c = canvas(out[i]);
    TH1* fr = (TH1*)d->Clone(Form("fr_%d", i)); fr->Reset();
    fr->SetMinimum(0); fr->SetMaximum(1.0);
    if (i == 0) fr->GetXaxis()->SetRangeUser(0, 3.5);
    if (i == 1) fr->GetXaxis()->SetRangeUser(-1.1, 1.1);
    fr->GetYaxis()->SetTitle("#varepsilon_{match} = N_{TOF}/N_{TPC}");
    fr->Draw("axis");
    g->SetMarkerStyle(20); g->SetMarkerSize(i == 2 ? 0.5 : 1.0); g->SetMarkerColor(kBlue + 2); g->SetLineColor(kBlue + 2);
    g->Draw("p same");
    tag();
    label(0.18, 0.84, "primary, |#eta|<1, matchFlag>0, #beta>0", 0.04);
    save(c, out[i]);
  }
  printf("figures done\n");
}
