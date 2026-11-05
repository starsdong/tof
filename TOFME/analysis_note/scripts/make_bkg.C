// make_bkg.C : background diagnostics of SE vs overlay-ME m2 (StTofMEMaker v4, 0-80%).
// usage (from analysis_note/): root -b -q -l 'scripts/make_bkg.C("data/me_v4.root")'
// For category variable k (see hSEcat<k> title) and category c, in fine p bins j:
//   ME_norm(c,j,m2) = ME(c,j,m2) * NP(c,j)/NT(c,j) * (n1+n2)/n1     (n1, n2: overlay status 1 / 2 in (c,j))
// hit categories (k = 7, 8, 9) use NP, NT, n1, n2 summed over categories (from k = 0).
// Windows (GeV^2/c^4): HI 1.6 < m2 < 2.8 (between p and d), LO m2 < -0.15, PR 0.6 < m2 < 1.2 (protons),
//   PI -0.1 < m2 < 0.1.
// Outputs: logs/bkg_v4.txt, tables/bkg_v4.tex (macros), tables/bkg_cat_v4.tex, figures/bkg_v4_*.pdf/png
#include <fstream>

namespace bk {
TFile* gF;
std::ofstream gL, gT;
const bool kHit[12] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0};
const char* kVar[12] = {"gDCA", "|V_{z}-V_{z}^{VPD}|", "|z_{local}|", "|y_{local}|", "nHitsFit", "|#eta|", "centrality",
                        "hit cell", "hit owner", "hits in cell", "n#sigma_{p}", "charge"};
const char* kVarTex[12] = {"gDCA", "$|V_z-V_z^{\\rm VPD}|$", "$|z_{\\rm local}|$", "$|y_{\\rm local}|$", "nHitsFit", "$|\\eta|$",
                           "centrality", "hit cell", "hit owner", "hits in cell", "$n\\sigma_{\\rm p}$", "charge"};
const char* kCat[12][4] = {{"<0.5 cm", "0.5-1 cm", "1-2 cm", "2-3 cm"}, {"<3 cm", "3-6 cm", ">6 cm", "no VPD"},
                           {"<1 cm", "1-2 cm", "2-2.6 cm", ">2.6 cm"}, {"<0.6 cm", "0.6-1.2 cm", ">1.2 cm", ""},
                           {"15-24", "25-34", ">=35", ""}, {"<0.3", "0.3-0.6", ">0.6", ""},
                           {"0-10%", "10-40%", "40-80%", ""}, {"projected", "neighbour", "", ""},
                           {"own", "own, shared", "other track", "other, shared"}, {"1", ">=2", "", ""},
                           {"<-3", "-3..3", "3..8", ">8"}, {"-", "+"}};
double integ(TH1* h, double a, double b) { return h->Integral(h->FindBin(a + 1e-6), h->FindBin(b - 1e-6)); }

// SE and normalised ME m2 for category c of variable k, p in [p1,p2); np returns projected tracks
void get(int k, int c, double p1, double p2, TH1D*& se, TH1D*& me, double& np, const char* id)
{
  TH3* hs = (TH3*)gF->Get(Form("hSEcat%d", k)); TH3* hm = (TH3*)gF->Get(Form("hMEcat%d", k));
  const int kn = kHit[k] ? 0 : k;
  TH2* NP = (TH2*)gF->Get(Form("hNPcat%d", kn)); TH2* NT = (TH2*)gF->Get(Form("hNTcat%d", kn));
  TH3* OS = (TH3*)gF->Get(Form("hOvStcat%d", kn));
  const int j1 = hs->GetYaxis()->FindBin(p1 + 1e-6), j2 = hs->GetYaxis()->FindBin(p2 - 1e-6);
  se = hs->ProjectionZ(Form("se_%s", id), c + 1, c + 1, j1, j2, "e");
  me = (TH1D*)se->Clone(Form("me_%s", id)); me->Reset();
  np = 0;
  const int nc = NP->GetNbinsX();
  for (int j = j1; j <= j2; j++) {
    double nP = 0, nT = 0, s1 = 0, s2 = 0;
    const int c1 = kHit[k] ? 1 : c + 1, c2 = kHit[k] ? nc : c + 1;
    for (int cc = c1; cc <= c2; cc++) {
      nP += NP->GetBinContent(cc, j); nT += NT->GetBinContent(cc, j);
      s1 += OS->GetBinContent(cc, j, 2); s2 += OS->GetBinContent(cc, j, 3);
    }
    if (!kHit[k]) np += nP; else np = -1;
    if (nT <= 0 || s1 <= 0) continue;
    TH1D* t = hm->ProjectionZ("tmpme", c + 1, c + 1, j, j, "e");
    me->Add(t, nP / nT * (s1 + s2) / s1);
    delete t;
  }
}
void num(const char* n, double v, const char* fmt)
{
  TString s = Form(fmt, v);
  gL << n << " = " << s << "\n"; gT << "\\newcommand{\\" << n << "}{" << s << "}\n";
}
}  // namespace bk

void drawCats(int k, double p1, double p2, const char* tag)
{
  using namespace bk;
  TH3* hs = (TH3*)gF->Get(Form("hSEcat%d", k));
  const int nb = hs->GetNbinsX();
  TCanvas* cv = new TCanvas(Form("cv%d_%s", k, tag), "", 400 * nb, 420);
  cv->Divide(nb, 1, 0.001, 0.001);
  for (int c = 0; c < nb; c++) {
    cv->cd(c + 1); gPad->SetLogy(1); gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.15);
    TH1D *se, *me; double np;
    get(k, c, p1, p2, se, me, np, Form("d%d_%d_%s", k, c, tag));
    se->Rebin(5); me->Rebin(5);
    se->SetTitle(""); se->GetXaxis()->SetTitle("m^{2} (GeV^{2}/c^{4})"); se->GetYaxis()->SetTitle("counts / 0.1 GeV^{2}/c^{4}");
    se->GetXaxis()->SetRangeUser(-1, 4.5); se->SetMarkerStyle(20); se->SetMarkerSize(0.6);
    se->GetXaxis()->SetTitleSize(0.06); se->GetYaxis()->SetTitleSize(0.055); se->GetXaxis()->SetLabelSize(0.055); se->GetYaxis()->SetLabelSize(0.055);
    se->GetYaxis()->SetTitleOffset(1.4);
    se->SetMinimum(0.5); se->SetMaximum(se->GetMaximum() * 20);
    se->Draw("pe");
    me->SetLineColor(kRed + 1); me->SetFillColorAlpha(kRed - 9, 0.5); me->SetLineWidth(2); me->Draw("hist same"); se->Draw("pe same");
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.055);
    l.DrawLatex(0.22, 0.88, Form("%s: %s", kVar[k], kCat[k][c]));
    if (c == 0) { l.SetTextSize(0.048); l.DrawLatex(0.22, 0.81, Form("%.1f < p < %.1f GeV/c, 0-80%%", p1, p2)); }
    if (c == nb - 1) {
      TLegend* lg = new TLegend(0.55, 0.62, 0.97, 0.78); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.05);
      lg->AddEntry(se, "same event", "p"); lg->AddEntry(me, "mixed event", "f"); lg->Draw();
    }
  }
  cv->SaveAs(Form("figures/bkg_v4_k%d_%s.pdf", k, tag)); cv->SaveAs(Form("figures/bkg_v4_k%d_%s.png", k, tag));
}

void make_bkg(const char* fn = "data/me_v4.root")
{
  using namespace bk;
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
  gF = TFile::Open(fn);
  gL.open("logs/bkg_v4.txt"); gT.open("tables/bkg_v4.tex");
  gL << "# generated by scripts/make_bkg.C from " << fn << "\n";
  gT << "% generated by scripts/make_bkg.C from " << fn << " -- do not edit\n";
  // ---- table: per category, two p ranges: HI and LO per 10^4 projected tracks for SE and ME, SE/ME, SE_HI/SE_PR
  const double pr[2][2] = {{0.6, 1.5}, {1.5, 3.0}};
  std::ofstream t("tables/bkg_cat_v4.tex");
  t << "% generated by scripts/make_bkg.C\n\\begin{tabular}{llrrrrrrr}\n\\toprule\n"
       " & & \\multicolumn{3}{c}{$0.6<p<1.5$\\,GeV/$c$} & \\multicolumn{3}{c}{$1.5<p<3$\\,GeV/$c$} & \\\\\n"
       "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\n"
       "variable & category & SE & ME & SE/ME & SE & ME & SE/ME & share\\\\\n\\midrule\n";
  for (int k = 0; k < 12; k++) {
    TH3* hs = (TH3*)gF->Get(Form("hSEcat%d", k));
    const int nb = hs->GetNbinsX();
    double tot = 0; for (int c = 0; c < nb; c++) tot += hs->Integral(c + 1, c + 1, 1, 32, 1, 300);
    for (int c = 0; c < nb; c++) {
      double v[2][3];
      for (int ip = 0; ip < 2; ip++) {
        TH1D *se, *me; double np;
        get(k, c, pr[ip][0], pr[ip][1], se, me, np, Form("t%d_%d_%d", k, c, ip));
        const double norm = kHit[k] ? 1. : (np > 0 ? 1e4 / np : 0);
        v[ip][0] = integ(se, 1.6, 2.8) * norm; v[ip][1] = integ(me, 1.6, 2.8) * norm;
        v[ip][2] = v[ip][1] > 0 ? v[ip][0] / v[ip][1] : 0;
        gL << "k" << k << " " << kVar[k] << " [" << kCat[k][c] << "] p " << pr[ip][0] << "-" << pr[ip][1]
           << (kHit[k] ? " counts" : " per1e4") << " HI SE " << v[ip][0] << " ME " << v[ip][1] << " SE/ME " << v[ip][2]
           << " | LO SE " << integ(se, -1, -0.15) * norm << " ME " << integ(me, -1, -0.15) * norm
           << " | PR SE " << integ(se, 0.6, 1.2) * norm << " | PI SE " << integ(se, -0.1, 0.1) * norm
           << " | HI/PR SE " << integ(se, 1.6, 2.8) / std::max(1., integ(se, 0.6, 1.2))
           << " | np " << np << "\n";
      }
      const double share = hs->Integral(c + 1, c + 1, 1, 32, 1, 300) / tot * 100;
      TString cn = kCat[k][c]; cn.ReplaceAll("%", "\\%"); cn.ReplaceAll(">=", "$\\ge$"); cn.ReplaceAll("<", "$<$"); cn.ReplaceAll(">", "$>$");
      t << (c == 0 ? kVarTex[k] : "") << " & " << cn.Data()
        << Form(" & %.2f & %.2f & %.1f & %.2f & %.2f & %.1f & %.1f\\%%", v[0][0], v[0][1], v[0][2], v[1][0], v[1][1], v[1][2], share)
        << (k == 7 || k == 8 || k == 9 ? "$^{*}$" : "") << "\\\\\n";
    }
    if (k < 11) t << "\\midrule\n";
  }
  t << "\\bottomrule\n\\end{tabular}\n"; t.close();

  // ---- figures: categories at 0.8 < p < 1.5
  const int kd[7] = {0, 1, 2, 10, 11, 6, 7};
  for (int i = 0; i < 7; i++) drawCats(kd[i], 0.8, 1.5, "p08");
  drawCats(11, 0.4, 0.8, "p04");
  drawCats(10, 0.4, 0.8, "p04");

  // ---- time residual t - t_pi: SE (all, own hit) vs ME (overlay, timed only, per projected track ~ shape)
  {
    TH2* ds = (TH2*)gF->Get("hDtSE"); TH2* dm = (TH2*)gF->Get("hDtME");
    TCanvas* cv = new TCanvas("cdt", "", 1500, 450); cv->Divide(3, 1, 0.001, 0.001);
    const double pb[4] = {0.4, 0.8, 1.5, 3.0};
    for (int i = 0; i < 3; i++) {
      cv->cd(i + 1); gPad->SetLogy(1); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.14);
      const int j1 = ds->GetXaxis()->FindBin(pb[i] + 1e-6), j2 = ds->GetXaxis()->FindBin(pb[i + 1] - 1e-6);
      TH1D* s = ds->ProjectionY(Form("dts%d", i), j1, j2); TH1D* m = dm->ProjectionY(Form("dtm%d", i), j1, j2);
      s->Rebin(4); m->Rebin(4);
      // scale ME to SE in the late region 8-25 ns (shape comparison only)
      const double a = integ(s, 8, 25), b = integ(m, 8, 25);
      if (b > 0) m->Scale(a / b);
      s->SetTitle(""); s->GetXaxis()->SetTitle("t - t_{#pi} (ns)"); s->GetYaxis()->SetTitle("counts / 0.2 ns");
      s->SetMarkerStyle(20); s->SetMarkerSize(0.5); s->SetMinimum(0.5); s->Draw("pe");
      s->GetXaxis()->SetTitleSize(0.055); s->GetYaxis()->SetTitleSize(0.055); s->GetXaxis()->SetLabelSize(0.05); s->GetYaxis()->SetLabelSize(0.05);
      m->SetLineColor(kRed + 1); m->SetLineWidth(2); m->Draw("hist same");
      TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.05);
      l.DrawLatex(0.45, 0.85, Form("%.1f < p < %.1f GeV/c", pb[i], pb[i + 1]));
      if (i == 0) { TLegend* lg = new TLegend(0.45, 0.64, 0.97, 0.80); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
        lg->AddEntry(s, "same event", "p"); lg->AddEntry(m, "mixed event (scaled, 8-25 ns)", "l"); lg->Draw(); }
      gL << "dt p " << pb[i] << "-" << pb[i + 1] << " SE frac t-tpi>8ns " << integ(s, 8, 27) / s->Integral()
         << " ME(scaled) frac " << integ(m, 8, 27) / s->Integral() << "\n";
    }
    cv->SaveAs("figures/bkg_v4_dt.pdf"); cv->SaveAs("figures/bkg_v4_dt.png");
  }
  gL.close(); gT.close();
  printf("make_bkg done\n");
}
