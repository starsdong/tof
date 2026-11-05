// make_m2.C : same-event vs mixed-event m^2 per centrality and momentum bin.
// usage (from analysis_note/): root -b -q -l 'scripts/make_m2.C("data/me_m2.root","m2")'
// Inputs from StTofMEMaker v2: hM2SE, hM2ME (cell+-1), hM2SEc, hM2MEc (cell), hM2Prod,
//   hNPrjM2 (projected tracks), hNTrM2 (ME trials), hMEStat/hSEStat (0 none, 1 timed hit, 2 hit w/o time).
// ME normalisation in each fine p bin j of centrality c (no free parameter):
//   overlay ME (hM2MEo): partner hits added on top of this event, production-like choice redone; m2 filled
//     when a partner hit wins.  ME_oa = MEo(c,j) * Nproj(c,j)/Ntrials(c,j)  (partner hits with time only)
//                               ME_ob = ME_oa * (n1 + n2)/n1   (n1/n2 = plain-ME trials with timed/untimed hit:
//                                       untimed random hits assumed to have the timed-hit time distribution)
//   plain ME (hM2ME, no shadowing by the track's own hit): ME_b = ME * Nproj/Ntrials * (n1+n2)/n1
// Outputs: figures/m2_<tag>_cent<c>.pdf/png, figures/m2_<tag>_allcent.pdf/png, logs/m2_numbers_<tag>.txt,
//          tables/m2_<tag>.tex (macros), tables/m2_contam_<tag>.tex
#include <fstream>

namespace m2 {
const int kNPb = 9;
const double kPb[kNPb + 1] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0};
const char* kCentName[9] = {"70-80%", "60-70%", "50-60%", "40-50%", "30-40%", "20-30%", "10-20%", "5-10%", "0-5%"};
const int kRebin = 4;
TFile* gF = 0;
TString gTag, gMac;

// projection over centrality bins [c1,c2] (0..8) and p range [p1,p2): SE, ME_a, ME_b, Prod, and counts
struct Proj { TH1D *se, *mea, *meb, *prod, *moa, *mob, *mtb, *mte, *mts, *mcl; double nProj, nTr, n1, n2; };
int gMode = 0;          // 0: overlay projected-first, 1: ToT-blind production choice, 2: earliest time wins,
                        // 3: earliest time wins with sampled times for untimed hits (v9)
bool gUseTot = false;   // true: main ME = ToT-blind production choice (hM2MEtot), else overlay projected-first
Proj project(const char* sname, const char* mname, int c1, int c2, double p1, double p2, const char* id)
{
  TH3* hs = (TH3*)gF->Get(sname); TH3* hm = (TH3*)gF->Get(mname); TH3* hp = (TH3*)gF->Get("hM2Prod");
  TH2* np = (TH2*)gF->Get("hNPrjM2"); TH2* nt = (TH2*)gF->Get("hNTrM2"); TH3* st = (TH3*)gF->Get("hMEStat");
  const int j1 = hs->GetYaxis()->FindBin(p1 + 1e-6), j2 = hs->GetYaxis()->FindBin(p2 - 1e-6);
  Proj r;
  r.se = hs->ProjectionZ(Form("se_%s", id), c1 + 1, c2 + 1, j1, j2, "e");
  r.prod = hp->ProjectionZ(Form("pr_%s", id), c1 + 1, c2 + 1, j1, j2, "e");
  r.mea = (TH1D*)r.se->Clone(Form("mea_%s", id)); r.mea->Reset();
  r.meb = (TH1D*)r.se->Clone(Form("meb_%s", id)); r.meb->Reset();
  r.moa = (TH1D*)r.se->Clone(Form("moa_%s", id)); r.moa->Reset();
  r.mob = (TH1D*)r.se->Clone(Form("mob_%s", id)); r.mob->Reset();
  r.mtb = (TH1D*)r.se->Clone(Form("mtb_%s", id)); r.mtb->Reset();
  r.mte = (TH1D*)r.se->Clone(Form("mte_%s", id)); r.mte->Reset();
  TH3* hte = (TH3*)gF->Get("hM2MEearly");
  r.mts = (TH1D*)r.se->Clone(Form("mts_%s", id)); r.mts->Reset();
  r.mcl = (TH1D*)r.se->Clone(Form("mcl_%s", id)); r.mcl->Reset();
  TH3* hsa = (TH3*)gF->Get("hM2MEsamp"); TH3* hcl = (TH3*)gF->Get("hM2MEclos"); TH2* nts = (TH2*)gF->Get("hNTrSamp");
  TH3* ht = (TH3*)gF->Get("hM2MEtot"); TH3* hts = (TH3*)gF->Get("hMEtotSt");
  TH3* ho = (TH3*)gF->Get("hM2MEo");
  r.nProj = r.nTr = r.n1 = r.n2 = 0;
  for (int c = c1; c <= c2; c++)
    for (int j = j1; j <= j2; j++) {
      const double nP = np->GetBinContent(c + 1, j), nT = nt->GetBinContent(c + 1, j);
      const double s1 = st->GetBinContent(c + 1, j, 2), s2 = st->GetBinContent(c + 1, j, 3);
      r.nProj += nP; r.nTr += nT; r.n1 += s1; r.n2 += s2;
      if (nT <= 0 || s1 <= 0) continue;
      TH1D* t = hm->ProjectionZ("tmp_me", c + 1, c + 1, j, j, "e");
      r.mea->Add(t, nP / nT);
      r.meb->Add(t, nP / nT * (s1 + s2) / s1);
      delete t;
      if (hte) {
        TH1D* o = hte->ProjectionZ("tmp_me8", c + 1, c + 1, j, j, "e");
        r.mte->Add(o, nP / nT);
        delete o;
      }
      if (hsa && nts && nts->GetBinContent(c + 1, j) > 0) {
        TH1D* o = hsa->ProjectionZ("tmp_ms", c + 1, c + 1, j, j, "e");
        r.mts->Add(o, nP / nts->GetBinContent(c + 1, j));
        delete o;
      }
      if (hcl && nts && nts->GetBinContent(c + 1, j) > 0) {
        TH1D* o = hcl->ProjectionZ("tmp_mc", c + 1, c + 1, j, j, "e");
        r.mcl->Add(o, nP / nts->GetBinContent(c + 1, j));
        delete o;
      }
      if (ht && hts) {
        const double t1 = hts->GetBinContent(c + 1, j, 1), t2 = hts->GetBinContent(c + 1, j, 2);
        if (t1 > 0) {
          TH1D* o = ht->ProjectionZ("tmp_mt", c + 1, c + 1, j, j, "e");
          r.mtb->Add(o, nP / nT * (t1 + t2) / t1);
          delete o;
        }
      }
      if (ho) {
        TH1D* o = ho->ProjectionZ("tmp_mo", c + 1, c + 1, j, j, "e");
        r.moa->Add(o, nP / nT);
        r.mob->Add(o, nP / nT * (s1 + s2) / s1);
        delete o;
      }
    }
  if (gMode == 1) { r.moa = r.mob; r.mob = r.mtb; }   // moa = previous overlay (all hits), drawn dashed
  if (gMode == 2) { r.moa = r.mtb; r.mob = r.mte; }   // main = earliest time wins, dashed = ToT-blind
  if (gMode == 3) { r.moa = r.mte; r.mob = r.mts; }   // main = sampled times, dashed = timed hits only (v8)
  return r;
}
double integ(TH1* h, double a, double b) { return h->Integral(h->FindBin(a + 1e-6), h->FindBin(b - 1e-6)); }

std::ofstream gL, gT;
void num(const char* name, double v, const char* fmt)
{
  TString s = Form(fmt, v);
  gL << name << " = " << s << "\n";
  gT << "\\newcommand{\\" << name << gMac << "}{" << s << "}\n";
}
}  // namespace m2

void drawPanel(m2::Proj& r, double p1, double p2, bool first, const char* centLab)
{
  using namespace m2;
  TH1D* se = (TH1D*)r.se->Clone(); se->Rebin(kRebin);
  TH1D* ma = (TH1D*)r.moa->Clone(); ma->Rebin(kRebin);
  TH1D* mb = (TH1D*)r.mob->Clone(); mb->Rebin(kRebin);
  TH1D* mp = (TH1D*)r.meb->Clone(); mp->Rebin(kRebin);
  if (gUseTot) { ma->SetLineColor(kBlue + 1); ma->SetLineWidth(2); ma->SetLineStyle(2); }
  se->SetTitle(""); se->GetXaxis()->SetTitle("m^{2} (GeV^{2}/c^{4})"); se->GetYaxis()->SetTitle("counts / 0.04 GeV^{2}/c^{4}");
  se->GetXaxis()->SetRangeUser(-0.6, 4.2);
  se->SetMarkerStyle(20); se->SetMarkerSize(0.5); se->SetMarkerColor(kBlack); se->SetLineColor(kBlack);
  double mx = se->GetMaximum(), mn = 1e30;
  for (int b = se->FindBin(-0.6); b <= se->FindBin(4.2); b++) {
    if (se->GetBinContent(b) > 0) mn = std::min(mn, se->GetBinContent(b));
    if (mb->GetBinContent(b) > 0) mn = std::min(mn, mb->GetBinContent(b));
  }
  se->SetMaximum(mx * 12); se->SetMinimum(std::max(0.3, mn * 0.5));
  se->GetXaxis()->SetLabelSize(0.06); se->GetYaxis()->SetLabelSize(0.06);
  se->GetXaxis()->SetTitleSize(0.065); se->GetYaxis()->SetTitleSize(0.06);
  se->GetXaxis()->SetTitleOffset(1.0); se->GetYaxis()->SetTitleOffset(1.25);
  se->Draw("pe");
  mb->SetLineColor(kRed + 1); mb->SetLineWidth(2); mb->SetFillColorAlpha(kRed - 9, 0.5); mb->SetFillStyle(1001);
  mb->Draw("hist same");

  if (gUseTot) ma->Draw("hist same");
  se->Draw("pe same");
  TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.065);
  l.DrawLatex(0.45, 0.82, Form("%.1f < p < %.1f GeV/c", p1, p2));
  if (first) {
    l.SetTextSize(0.06); l.DrawLatex(0.45, 0.74, centLab);
    TLegend* lg = new TLegend(0.36, 0.48, 0.97, 0.71); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.05);
    lg->AddEntry(se, "same event", "p");
    lg->AddEntry(mb, gMode == 3 ? "ME, sampled times" : gMode == 2 ? "ME, earliest time wins" : (gMode == 1 ? "ME, production choice" : "mixed event"), "f");
    if (gMode == 1) lg->AddEntry(ma, "ME, projected cell first", "l");
    if (gMode == 2) lg->AddEntry(ma, "ME, equal-chance choice", "l");
    if (gMode == 3) lg->AddEntry(ma, "ME, timed hits only", "l");
    lg->Draw();
  }
}

void make_m2(const char* fn = "data/me_m2.root", const char* tag = "m2", const char* lab = "", bool useTot = false)
{
  using namespace m2;
  gUseTot = useTot; gMode = useTot ? 1 : 0;
  if (TString(tag) == "v8") { gMode = 2; gUseTot = true; }
  if (TString(tag) == "v9") { gMode = 3; gUseTot = true; }
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
  gF = TFile::Open(fn);
  gTag = tag; gMac = tag; gMac.ReplaceAll("_", ""); gMac.ReplaceAll("test", "Test"); gMac.ReplaceAll("m2", "Mtwo");
  if (gTag == "v3") gMac = "Vthree";
  if (gTag == "v7") gMac = "Vseven";
  if (gTag == "v8") gMac = "Veight";
  if (gTag == "v9") gMac = "Vnine";
  gMac.ReplaceAll("0", "Z"); gMac.ReplaceAll("1", "O"); gMac.ReplaceAll("2", "T"); gMac.ReplaceAll("3", "Th");
  gL.open(Form("logs/m2_numbers_%s.txt", tag)); gT.open(Form("tables/m2_%s.tex", tag));
  gL << "# generated by scripts/make_m2.C from " << fn << "\n";
  gT << "% generated by scripts/make_m2.C from " << fn << " -- do not edit\n";
  TH1* cf = (TH1*)gF->Get("hCentFlag");
  num("mtEvtSel", cf->GetBinContent(1), "%.0f"); num("mtEvtBad", cf->GetBinContent(2), "%.0f");
  num("mtEvtPile", cf->GetBinContent(3), "%.0f"); num("mtEvtPeriph", cf->GetBinContent(4), "%.0f");
  num("mtEvtNoTime", cf->GetBinContent(5), "%.0f"); num("mtEvtUsed", cf->GetBinContent(6), "%.0f");
  num("mtTimedFrac", 100 * ((TH1*)gF->Get("hTimedFrac"))->GetMean(), "%.1f");
  TH2* dl = (TH2*)gF->Get("hDLVsP"); TH2* dm = (TH2*)gF->Get("hDm2VsP");
  num("mtDLMean", dl->GetMean(2), "%.3f"); num("mtDLRms", dl->GetRMS(2), "%.2f");
  num("mtDmRms", dm->GetRMS(2), "%.4f");
  { TH3* st = (TH3*)gF->Get("hMEStat"); TH3* ss = (TH3*)gF->Get("hSEStat");
    const double a1 = st->Integral(1, 9, 1, st->GetNbinsY(), 2, 2), a2 = st->Integral(1, 9, 1, st->GetNbinsY(), 3, 3);
    const double b1 = ss->Integral(1, 9, 1, ss->GetNbinsY(), 2, 2), b2 = ss->Integral(1, 9, 1, ss->GetNbinsY(), 3, 3);
    num("mtMETimedFrac", 100 * a1 / (a1 + a2), "%.1f"); num("mtSETimedFrac", 100 * b1 / (b1 + b2), "%.1f"); }

  // per-centrality canvases, 3x3 momentum bins
  for (int c = 8; c >= 0; c--) {
    TCanvas* cv = new TCanvas(Form("c%d", c), "", 1500, 1300);
    cv->Divide(3, 3, 0.001, 0.001);
    for (int ib = 0; ib < kNPb; ib++) {
      cv->cd(ib + 1); gPad->SetLogy(1);
      gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.15);
      Proj r = project("hM2SE", "hM2ME", c, c, kPb[ib], kPb[ib + 1], Form("c%d_%d", c, ib));
      drawPanel(r, kPb[ib], kPb[ib + 1], ib == 0, Form("Au+Au 19.6 GeV, %s", kCentName[c]));
    }
    cv->SaveAs(Form("figures/m2_%s_cent%d.pdf", tag, c)); cv->SaveAs(Form("figures/m2_%s_cent%d.png", tag, c));
  }
  { TCanvas* cv = new TCanvas("call", "", 1500, 1300);
    cv->Divide(3, 3, 0.001, 0.001);
    for (int ib = 0; ib < kNPb; ib++) {
      cv->cd(ib + 1); gPad->SetLogy(1);
      gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.15);
      Proj r = project("hM2SE", "hM2ME", 0, 8, kPb[ib], kPb[ib + 1], Form("all_%d", ib));
      drawPanel(r, kPb[ib], kPb[ib + 1], ib == 0, "Au+Au 19.6 GeV, 0-80%");
    }
    cv->SaveAs(Form("figures/m2_%s_allcent.pdf", tag)); cv->SaveAs(Form("figures/m2_%s_allcent.png", tag)); }

  // contamination table: ME_b / SE in PID windows; sideband closure SE/ME_b
  // windows: pi -0.1..0.1 ; K 0.17..0.33 ; p 0.6..1.2 ; sideband 1.6..2.8 (between p and d)
  const double win[4][2] = {{-0.1, 0.1}, {0.17, 0.33}, {0.6, 1.2}, {1.6, 2.8}};
  std::ofstream t(Form("tables/m2_contam_%s.tex", tag));
  t << "% generated by scripts/make_m2.C\n\\begin{tabular}{llccccc}\n\\toprule\n"
       " & & \\multicolumn{3}{c}{random fraction under peak (\\%)} & SE/ME & max\\\\\n"
       "\\cmidrule(lr){3-5}\n"
       "centrality & $p$ (GeV/$c$) & $\\pi$ & K & p & $1.6<m^2<2.8$ & ME/SE\\\\\n\\midrule\n";
  const int cSel[3] = {8, 4, 0};
  for (int ic = 0; ic < 3; ic++) {
    const int c = cSel[ic];
    for (int ib = 0; ib < kNPb; ib++) {
      Proj r = project("hM2SE", "hM2ME", c, c, kPb[ib], kPb[ib + 1], Form("t%d_%d", c, ib));
      double fr[4], frt[4];
      for (int w = 0; w < 3; w++) {
        const double s = integ(r.se, win[w][0], win[w][1]);
        fr[w] = s > 0 ? 100 * integ(r.mob, win[w][0], win[w][1]) / s : 0;
        frt[w] = s > 0 ? 100 * integ(r.moa, win[w][0], win[w][1]) / s : 0;
      }
      // largest ME/SE over m2 bins (0.04 wide, SE >= 50 counts, -0.3 < m2 < 3): the random part cannot exceed SE
      double mxA = 0, mxT = 0, mxP = 0;
      { TH1D* s4 = (TH1D*)r.se->Clone("s4"); s4->Rebin(kRebin); TH1D* a4 = (TH1D*)r.mob->Clone("a4"); a4->Rebin(kRebin);
        TH1D* t4 = (TH1D*)r.moa->Clone("t4"); t4->Rebin(kRebin); TH1D* p4 = (TH1D*)r.meb->Clone("p4"); p4->Rebin(kRebin);
        for (int k = s4->FindBin(-0.3); k <= s4->FindBin(3.0); k++) {
          const double sv = s4->GetBinContent(k); if (sv < 50) continue;
          mxA = std::max(mxA, a4->GetBinContent(k) / sv); mxT = std::max(mxT, t4->GetBinContent(k) / sv);
          mxP = std::max(mxP, p4->GetBinContent(k) / sv);
        }
        delete s4; delete a4; delete t4; delete p4; }
      const double sb = integ(r.mob, win[3][0], win[3][1]);
      const double sbr = sb > 0 ? integ(r.se, win[3][0], win[3][1]) / sb : 0;
      const double ng = integ(r.moa, win[3][0], win[3][1]);
      const double ngr = ng > 0 ? integ(r.se, win[3][0], win[3][1]) / ng : 0;
      TString cn = kCentName[c]; cn.ReplaceAll("-", "--"); cn.ReplaceAll("%", "\\%");
      t << (ib == 0 ? cn.Data() : "") << " & " << Form("%.1f--%.1f", kPb[ib], kPb[ib + 1])
        << Form(" & %.2f & %.2f & %.2f & %.1f & %.2f", fr[0], fr[1], fr[2], sbr, mxA) << "\\\\\n";
      gL << "cent " << kCentName[c] << " p " << kPb[ib] << "-" << kPb[ib + 1] << " ME_b/SE % pi " << fr[0] << " K " << fr[1]
         << " p " << fr[2] << " | SE/MEo_all sideband(1.6-2.8) " << sbr << " SE/MEo_timed " << ngr
         << " | MEo_all/ME_plain " << (r.meb->Integral() > 0 ? r.mob->Integral() / r.meb->Integral() : 0)
         << " | max ME/SE: overlay-all " << mxA << " overlay-timed " << mxT << " plain " << mxP
         << " | timed-only frac % pi " << frt[0] << " K " << frt[1] << " p " << frt[2]
         << " | ME_b/ME_a " << (r.n1 > 0 ? (r.n1 + r.n2) / r.n1 : 0) << "\n";
    }
    if (ic < 2) t << "\\midrule\n";
  }
  t << "\\bottomrule\n\\end{tabular}\n"; t.close();

  // summary: random fraction under the pi/K/p windows vs p for selected centralities
  { const int cs[5] = {8, 6, 4, 2, 0}; const int col[5] = {kRed + 1, kOrange + 7, kGreen + 2, kAzure + 2, kViolet + 1};
    const int mk[5] = {20, 21, 33, 34, 29};
    const char* sp[3] = {"#pi", "K", "p"};
    TCanvas* cv = new TCanvas("csum", "", 1800, 600); cv->Divide(3, 1, 0.001, 0.001);
    TGraph* g[3][5];
    for (int w = 0; w < 3; w++) for (int ic = 0; ic < 5; ic++) { g[w][ic] = new TGraph(); }
    for (int ic = 0; ic < 5; ic++)
      for (int ib = 0; ib < kNPb; ib++) {
        Proj r = project("hM2SE", "hM2ME", cs[ic], cs[ic], kPb[ib], kPb[ib + 1], Form("sum%d_%d", ic, ib));
        for (int w = 0; w < 3; w++) {
          const double sv = integ(r.se, win[w][0], win[w][1]);
          if (sv > 100) g[w][ic]->SetPoint(g[w][ic]->GetN(), 0.5 * (kPb[ib] + kPb[ib + 1]), 100 * integ(r.mob, win[w][0], win[w][1]) / sv);
        }
      }
    for (int w = 0; w < 3; w++) {
      cv->cd(w + 1); gPad->SetLogy(1); gPad->SetLeftMargin(0.21); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.05);
      TH1F* fr = new TH1F(Form("frs%d", w), Form(";p (GeV/c);random fraction under %s peak (%%)", sp[w]), 10, 0, 3.1);
      fr->SetMinimum(0.01); fr->SetMaximum(2000); fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05);
      fr->GetXaxis()->SetLabelSize(0.05); fr->GetYaxis()->SetLabelSize(0.05); fr->GetYaxis()->SetTitleOffset(1.75); fr->Draw("axis");
      TLegend* lg = new TLegend(0.22, 0.60, 0.95, 0.84); lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
      for (int ic = 0; ic < 5; ic++) {
        g[w][ic]->SetMarkerStyle(mk[ic]); g[w][ic]->SetMarkerColor(col[ic]); g[w][ic]->SetLineColor(col[ic]); g[w][ic]->SetMarkerSize(1.3);
        g[w][ic]->Draw("pl same"); lg->AddEntry(g[w][ic], kCentName[cs[ic]], "pl");
      }
      if (w == 0) lg->Draw();
      TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.05);
      l.DrawLatex(0.22, 0.88, Form("%s window: %.2f < m^{2} < %.2f GeV^{2}/c^{4}", sp[w], win[w][0], win[w][1]));
    }
    cv->SaveAs(Form("figures/m2_%s_frand.pdf", tag)); cv->SaveAs(Form("figures/m2_%s_frand.png", tag));
    // numbers for the text: 0-5% and 30-40%, p 0.2-0.4 and 0.8-1.0
    const char* wn[3] = {"Pi", "K", "P"};
    for (int ic = 0; ic < 5; ic += 2)
      for (int w = 0; w < 3; w++) {
        double xx, yy, lo = -1, mid = -1;
        for (int k = 0; k < g[w][ic]->GetN(); k++) { g[w][ic]->GetPoint(k, xx, yy); if (fabs(xx - 0.3) < 0.01) lo = yy; if (fabs(xx - 0.9) < 0.01) mid = yy; }
        num(Form("mtFr%sLoC%c", wn[w], 'a' + cs[ic]), lo, "%.1f"); num(Form("mtFr%sMidC%c", wn[w], 'a' + cs[ic]), mid, "%.1f");
      }
  }
  // valley test summary for 0-5%, 0.4 < p < 1.5 (6 momentum bins)
  { double mA = 0, mT = 0, mP = 0, mPlo = 1e9;
    for (int ib = 1; ib <= 5; ib++) {
      Proj r = project("hM2SE", "hM2ME", 8, 8, kPb[ib], kPb[ib + 1], Form("v_%d", ib));
      TH1D* s4 = (TH1D*)r.se->Clone("vs4"); s4->Rebin(kRebin); TH1D* a4 = (TH1D*)r.mob->Clone("va4"); a4->Rebin(kRebin);
      TH1D* t4 = (TH1D*)r.moa->Clone("vt4"); t4->Rebin(kRebin); TH1D* p4 = (TH1D*)r.meb->Clone("vp4"); p4->Rebin(kRebin);
      double pA = 0, pT = 0, pP = 0;
      for (int k = s4->FindBin(-0.3); k <= s4->FindBin(3.0); k++) {
        const double sv = s4->GetBinContent(k); if (sv < 50) continue;
        pA = std::max(pA, a4->GetBinContent(k) / sv); pT = std::max(pT, t4->GetBinContent(k) / sv); pP = std::max(pP, p4->GetBinContent(k) / sv);
      }
      mA = std::max(mA, pA); mT = std::max(mT, pT); mP = std::max(mP, pP); mPlo = std::min(mPlo, pP);
      delete s4; delete a4; delete t4; delete p4;
    }
    num("mtValleyAll", mA, "%.2f"); num("mtValleyTimed", mT, "%.2f"); num("mtValleyPlainMin", mPlo, "%.1f"); num("mtValleyPlainMax", mP, "%.1f"); }
  // integrated numbers (0-80%, 0.2<p<3): sideband closure and timed/all factor per centrality
  Proj all = project("hM2SE", "hM2ME", 0, 8, 0.2, 3.0, "allint");
  num("mtSbRatioAll", integ(all.se, 1.6, 2.8) / integ(all.mob, 1.6, 2.8), "%.2f");
  num("mtSbRatioTimed", integ(all.se, 1.6, 2.8) / integ(all.moa, 1.6, 2.8), "%.2f");
  num("mtShadowFrac", 100 * all.mob->Integral() / all.meb->Integral(), "%.0f");
  num("mtTimedScaleAll", (all.n1 + all.n2) / all.n1, "%.2f");
  num("mtSENorm", all.se->Integral(0, all.se->GetNbinsX() + 1), "%.0f");
  num("mtMENormFrac", 100 * all.mob->Integral(0, all.mob->GetNbinsX() + 1) / all.se->Integral(0, all.se->GetNbinsX() + 1), "%.1f");
  num("mtMENormFracTimed", 100 * all.moa->Integral(0, all.moa->GetNbinsX() + 1) / all.se->Integral(0, all.se->GetNbinsX() + 1), "%.1f");
  for (int c = 0; c < 9; c += 4) {
    Proj r = project("hM2SE", "hM2ME", c, c, 0.2, 3.0, Form("ci%d", c));
    TString nm = Form("mtTimedScaleC%c", 'a' + c);
    num(nm.Data(), (r.n1 + r.n2) / r.n1, "%.2f");
    num(Form("mtMENormFracC%c", 'a' + c), 100 * r.mob->Integral(0, r.mob->GetNbinsX() + 1) / r.se->Integral(0, r.se->GetNbinsX() + 1), "%.1f");
    num(Form("mtMENormFracTimedC%c", 'a' + c), 100 * r.moa->Integral(0, r.moa->GetNbinsX() + 1) / r.se->Integral(0, r.se->GetNbinsX() + 1), "%.1f");
  }
  gL.close(); gT.close();
  printf("make_m2 done (%s)\n", tag);
}
