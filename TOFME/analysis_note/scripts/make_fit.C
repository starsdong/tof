// make_fit.C : random-match subtraction S = SE - ME (ME normalized per track) and PID fits of the m2 distributions.
// charge (5th argument): 0 all tracks, -1 negative, +1 positive (needs the <name>neg histograms of the v3/v4 makers).
// usage (analysis_note/):
//   root -b -q -l 'scripts/make_fit.C("data/ov1.root",0,"ovtfour","Overlay")'     method 2 (StTofMEOverlayMaker)
//   root -b -q -l 'scripts/make_fit.C("data/me_mu.root",1,"mufour","Emul")'       method 1 (StTofMEMuMaker)
// method 0: SE hM2SE, ME hM2ME, normalization hNTrk/hNTrial; method 1: SE hM2SE, ME hM2MEprodB, hNPrj/hNTr.
// Native 0.01 GeV^2/c^4 bins. The SE histogram is fitted (Poisson likelihood, function integrated over each bin) with
// signal + const + ME, the normalized ME fixed; this is equivalent to fitting S = SE - ME with correct statistics.
// Signal = four peaks (pi, K, p, d) in -0.6 < m2 < 4.8 GeV^2/c^4, each with one of four line shapes:
//   G  : Gaussian
//   G2 : double Gaussian (common mean; core width s1, wide width s2, core fraction f)
//   T  : Student-t, A/s * t_nu((m2 - mu)/s)
//   T2 : split Student-t (scale sL below and sR above the mean, common nu; normalized to A)
// plus a constant (free sign). The model with the smallest Pearson chi2/ndf (SE vs model) is "best"; the Gaussian
// family is {G, G2} and the Student-t family {T, T2}. Outputs: logs/fit_<tag>.txt, tables/fit_<tag>.tex (macro prefix
// fit<mp>), tables/fit_<tag>_tab.tex, figures/fit_<tag>_{cent,all}.pdf/png, data/sig_<tag>.root (SE, ME, S per bin).
#include <fstream>
namespace ft {
TFile* gF; int gMethod, gCharge = 0;
// charge 0: all tracks; -1: <name>neg; +1: <name> - <name>neg
TH3* get3(const char* name)
{
  if (gCharge == 0) return (TH3*)gF->Get(name);
  TH3* n = (TH3*)gF->Get(Form("%sneg", name));
  if (gCharge < 0) return n;
  TH3* a = (TH3*)((TH3*)gF->Get(name))->Clone(Form("%spos_%d", name, gRandom->Integer(1000000)));
  a->Add(n, -1); return a;
}
TH1D* proj(const char* name, int c1, int c2, double p1, double p2, bool me, const char* id)
{
  TH3* h = get3(name);
  const int j1 = h->GetYaxis()->FindBin(p1 + 1e-6), j2 = h->GetYaxis()->FindBin(p2 - 1e-6);
  if (!me) return h->ProjectionZ(Form("%s_%s", name, id), c1 + 1, c2 + 1, j1, j2, "e");
  TH2* nk = (TH2*)gF->Get(gMethod == 0 ? "hNTrk" : "hNPrj"); TH2* nt = (TH2*)gF->Get(gMethod == 0 ? "hNTrial" : "hNTr");
  TH1D* r = h->ProjectionZ(Form("%s_%s", name, id), 1, 1, 1, 1, "e"); r->Reset();
  for (int c = c1; c <= c2; c++) for (int j = j1; j <= j2; j++) {
    const double K = nk->GetBinContent(c + 1, j), T = nt->GetBinContent(c + 1, j);
    if (T <= 0) continue;
    TH1D* t = h->ProjectionZ("tmpft", c + 1, c + 1, j, j, "e"); r->Add(t, K / T); delete t;
  }
  return r;
}
const double kMu[4] = {0.01948, 0.24372, 0.88035, 3.5196};   // PDG m^2
const char* kSp[4] = {"pi", "K", "p", "d"};
const char* kModel[4] = {"G", "G2", "T", "T2"};
const double kLo = -0.6, kHi = 4.8;
const int kNPar = 21;   // per species k: [5k] A, [5k+1] mu, [5k+2] s1, [5k+3] s2, [5k+4] f or nu; [20] const
TH1* gME = 0;           // fixed random-match background (normalized ME) added to the fit model of SE
bool gAddME = true;
double meAt(double x) { return (gME && gAddME) ? gME->GetBinContent(gME->FindBin(x)) : 0.; }
double shape(int model, double x, const double* q)   // unit-area line shape of one species
{
  const double m = q[1], s1 = q[2], s2 = q[3], v = q[4], z = x - m;
  const double c = 1. / std::sqrt(2 * TMath::Pi());
  switch (model) {
    case 0: return c / s1 * std::exp(-0.5 * z * z / (s1 * s1));
    case 1: return v * c / s1 * std::exp(-0.5 * z * z / (s1 * s1)) + (1 - v) * c / s2 * std::exp(-0.5 * z * z / (s2 * s2));
    case 2: return TMath::Student(z / s1, v) / s1;
    default: { const double s = z < 0 ? s1 : s2; return 2. / (s1 + s2) * TMath::Student(z / s, v); }
  }
}
template <int M> double fModel(double* x, double* p)
{
  double s = 0;
  for (int k = 0; k < 4; k++) s += p[5 * k] * shape(M, x[0], p + 5 * k);
  return s * 0.01 + p[20] + meAt(x[0]);   // per 0.01 bin
}
template <int M> double fSig(double* x, double* p) { gAddME = false; const double v = fModel<M>(x, p); gAddME = true; return v; }
typedef double (*FPtr)(double*, double*);
FPtr kF[4] = {fModel<0>, fModel<1>, fModel<2>, fModel<3>};
FPtr kS[4] = {fSig<0>, fSig<1>, fSig<2>, fSig<3>};
// Pearson chi2 of the SE histogram against model f in [kLo, kHi]
double pearson(TH1* h, TF1* f, int& ndf)
{
  double c = 0; int n = 0;
  for (int b = h->FindBin(kLo + 1e-6); b <= h->FindBin(kHi - 1e-6); b++) {
    const double m = f->Integral(h->GetBinLowEdge(b), h->GetBinLowEdge(b + 1)) / h->GetBinWidth(b);
    if (m <= 0) continue;
    c += std::pow(h->GetBinContent(b) - m, 2) / m; n++;
  }
  ndf = n - f->GetNumberFreeParameters(); return c;
}
}

void make_fit(const char* fn = "data/ov1.root", int method = 0, const char* tag = "ovtfour", const char* mp = "Overlay", int charge = 0)
{
  using namespace ft;
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
  gF = TFile::Open(fn); gMethod = method; gCharge = charge;
  const char* qlab = charge < 0 ? ", negative" : (charge > 0 ? ", positive" : "");
  const char* meName = method == 0 ? "hM2ME" : "hM2MEprodB";
  std::ofstream L(Form("logs/fit_%s.txt", tag)); std::ofstream T(Form("tables/fit_%s.tex", tag));
  L << "# generated by scripts/make_fit.C from " << fn << " (method " << method << ", charge " << charge << ")\n";
  T << "% generated by scripts/make_fit.C -- do not edit\n";
  TFile* fo = new TFile(Form("data/sig_%s.root", tag), "RECREATE");
  const double pb[10] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0};
  const int cg[2][2] = {{8, 8}, {0, 8}};
  const char* cname[2] = {"0-5%", "0-80%"}; const char* ctag[2] = {"cent", "all"}; const char* cm[2] = {"Ci", "All"};
  std::vector<TString> tab;
  for (int g = 0; g < 2; g++) {
    TCanvas* cv = new TCanvas(Form("cf%d", g), "", 1500, 1300); cv->Divide(3, 3, 0.001, 0.001);
    int nBest[4] = {0, 0, 0, 0}, nTfam = 0, nFit = 0; double chiSum[4] = {0, 0, 0, 0};
    for (int i = 0; i < 9; i++) {
      const char* id = Form("%s_%d_%d", tag, g, i);
      const double pc = 0.5 * (pb[i] + pb[i + 1]);
      TH1D* se = proj("hM2SE", cg[g][0], cg[g][1], pb[i], pb[i + 1], false, id);
      TH1D* me = proj(meName, cg[g][0], cg[g][1], pb[i], pb[i + 1], true, id);
      TH1D* s = (TH1D*)se->Clone(Form("sig_%s", id)); s->Add(me, -1);
      fo->cd(); se->Write(); me->Write(); s->Write();
      // seeds
      double seedW[4], seedA[4];
      for (int k = 0; k < 4; k++) {
        seedW[k] = 0.012 + 0.035 * pc * pc * (k == 0 ? 0.6 : 1.0) + (k == 3 ? 0.05 : 0);
        seedA[k] = std::max(1., s->Integral(s->FindBin(kMu[k] - 1.5 * seedW[k]), s->FindBin(kMu[k] + 1.5 * seedW[k])) / 0.87);
      }
      const double dWmax = 0.08 + 0.12 * pc * pc;   // deuteron width bound (keeps it from absorbing a continuum)
      TF1* f[4]; double chi[4]; int ndf[4];
      for (int M = 0; M < 4; M++) {
        f[M] = new TF1(Form("f%s_%s", kModel[M], id), kF[M], kLo, kHi, kNPar);
        for (int k = 0; k < 4; k++) {
          const double wmax = k == 3 ? dWmax : 1.0;
          f[M]->SetParameter(5 * k, seedA[k]); f[M]->SetParLimits(5 * k, 0, 10 * s->Integral() + 10);
          f[M]->SetParameter(5 * k + 1, kMu[k]);
          f[M]->SetParLimits(5 * k + 1, kMu[k] - (k == 3 ? 0.18 : 0.12), kMu[k] + (k == 3 ? 0.18 : 0.12));
          f[M]->SetParameter(5 * k + 2, std::min(seedW[k], 0.9 * wmax)); f[M]->SetParLimits(5 * k + 2, 0.001, wmax);
          if (M == 0) { f[M]->FixParameter(5 * k + 3, 1); f[M]->FixParameter(5 * k + 4, 1); }
          if (M == 1) { f[M]->SetParameter(5 * k + 3, 2.5 * seedW[k]); f[M]->SetParLimits(5 * k + 3, 0.001, 1.5); f[M]->SetParameter(5 * k + 4, 0.85); f[M]->SetParLimits(5 * k + 4, 0.3, 1.0); }
          if (M == 2) { f[M]->FixParameter(5 * k + 3, 1); f[M]->SetParameter(5 * k + 4, 5); f[M]->SetParLimits(5 * k + 4, 1.0, 300); }
          if (M == 3) { f[M]->SetParameter(5 * k + 3, std::min(seedW[k], 0.9 * wmax)); f[M]->SetParLimits(5 * k + 3, 0.001, wmax); f[M]->SetParameter(5 * k + 4, 5); f[M]->SetParLimits(5 * k + 4, 1.0, 300); }
        }
        f[M]->SetParameter(20, 0); f[M]->SetParLimits(20, -100, 1e5);
      }
      gME = me;
      // G, then G2 and T seeded from G, T2 seeded from T
      se->Fit(f[0], "QRN0IL"); se->Fit(f[0], "QRN0IL");
      for (int M : {1, 2}) {
        for (int k = 0; k < 4; k++) { f[M]->SetParameter(5 * k, f[0]->GetParameter(5 * k)); f[M]->SetParameter(5 * k + 1, f[0]->GetParameter(5 * k + 1));
          f[M]->SetParameter(5 * k + 2, (M == 1 ? 0.9 : 1.0) * f[0]->GetParameter(5 * k + 2)); }
        f[M]->SetParameter(20, f[0]->GetParameter(20));
        se->Fit(f[M], "QRN0IL"); se->Fit(f[M], "QRN0IL");
      }
      for (int q = 0; q < kNPar; q++) if (q % 5 != 3 || q == 20) f[3]->SetParameter(q, f[2]->GetParameter(q));
      for (int k = 0; k < 4; k++) f[3]->SetParameter(5 * k + 3, f[2]->GetParameter(5 * k + 2));
      se->Fit(f[3], "QRN0IL"); se->Fit(f[3], "QRN0IL");
      int best = 0;
      for (int M = 0; M < 4; M++) { const double x = pearson(se, f[M], ndf[M]); chi[M] = ndf[M] > 0 ? x / ndf[M] : 1e9; chiSum[M] += chi[M]; if (chi[M] < chi[best]) best = M; }
      nBest[best]++; nFit++;
      const int bG = chi[1] < chi[0] ? 1 : 0, bT = chi[3] < chi[2] ? 3 : 2;
      if (chi[bT] < chi[bG]) nTfam++;
      TF1* fb = f[best];
      L << cname[g] << " p " << pb[i] << "-" << pb[i + 1] << " | SE " << se->Integral() << " ME " << me->Integral() << " S " << s->Integral()
        << " | chi2/ndf G " << chi[0] << " G2 " << chi[1] << " T " << chi[2] << " T2 " << chi[3] << " best " << kModel[best];
      double sumY = 0;
      for (int k = 0; k < 4; k++) {
        L << " | " << kSp[k] << " Y " << fb->GetParameter(5 * k) << " +- " << fb->GetParError(5 * k) << " mu " << fb->GetParameter(5 * k + 1)
          << " s1 " << fb->GetParameter(5 * k + 2) << " s2 " << fb->GetParameter(5 * k + 3) << " x " << fb->GetParameter(5 * k + 4);
        sumY += fb->GetParameter(5 * k);
      }
      L << " | const " << fb->GetParameter(20) << " per bin (" << 100 * fb->GetParameter(20) * (kHi - kLo) / 0.01 / std::max(1., s->Integral()) << " % of S)"
        << " | sumY/S " << sumY / std::max(1., s->Integral()) << "\n";
      tab.push_back(Form("%s & %.1f--%.1f & %.1f & %.1f & %.1f & %.1f & %s & %.0f & %.0f & %.0f & %.0f \\\\", TString(cname[g]).ReplaceAll("%", "\\%").Data(), pb[i], pb[i + 1],
                         chi[0], chi[1], chi[2], chi[3], kModel[best], fb->GetParameter(0), fb->GetParameter(5), fb->GetParameter(10), fb->GetParameter(15)));
      // draw: S with the best model of each family (signal part only)
      cv->cd(i + 1); gPad->SetLogy(1); gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.13);
      s->SetTitle(""); s->GetXaxis()->SetTitle("m^{2} (GeV^{2}/c^{4})"); s->GetYaxis()->SetTitle("SE - ME / 0.01 GeV^{2}/c^{4}");
      s->GetXaxis()->SetRangeUser(kLo, 4.5); s->SetMarkerStyle(20); s->SetMarkerSize(0.45); s->SetMinimum(0.5); s->SetMaximum(s->GetMaximum() * 30);
      s->GetXaxis()->SetTitleSize(0.055); s->GetYaxis()->SetTitleSize(0.05); s->GetXaxis()->SetLabelSize(0.05); s->GetYaxis()->SetLabelSize(0.05);
      s->Draw("pe");
      TF1* gs = new TF1(Form("gs_%s", id), kS[bG], kLo, kHi, kNPar); gs->SetParameters(f[bG]->GetParameters());
      TF1* ts = new TF1(Form("ts_%s", id), kS[bT], kLo, kHi, kNPar); ts->SetParameters(f[bT]->GetParameters());
      gs->SetNpx(1500); ts->SetNpx(1500);
      gs->SetLineColor(kBlue + 1); gs->SetLineWidth(2); gs->SetLineStyle(2); gs->Draw("same");
      ts->SetLineColor(kRed + 1); ts->SetLineWidth(2); ts->Draw("same");
      TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.052);
      l.DrawLatex(0.45, 0.88, Form("%.1f < p < %.1f GeV/c", pb[i], pb[i + 1]));
      l.SetTextSize(0.042);
      l.DrawLatex(0.45, 0.81, Form("#chi^{2}/ndf: %s %.1f, %s %.1f", kModel[bG], chi[bG], kModel[bT], chi[bT]));
      if (i == 0) {
        l.DrawLatex(0.45, 0.74, Form("%s, %s%s", mp, cname[g], qlab));
        TLegend* lg = new TLegend(0.40, 0.50, 0.97, 0.70); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
        lg->AddEntry(s, "SE - ME", "p"); lg->AddEntry(gs, "best Gaussian model", "l"); lg->AddEntry(ts, "best Student-t model", "l"); lg->Draw();
      }
    }
    cv->SaveAs(Form("figures/fit_%s_%s.pdf", tag, ctag[g])); cv->SaveAs(Form("figures/fit_%s_%s.png", tag, ctag[g]));
    L << cname[g] << " summary: best model G " << nBest[0] << " G2 " << nBest[1] << " T " << nBest[2] << " T2 " << nBest[3] << " of " << nFit
      << " | Student-t family better in " << nTfam << " | mean chi2/ndf G " << chiSum[0] / nFit << " G2 " << chiSum[1] / nFit << " T " << chiSum[2] / nFit << " T2 " << chiSum[3] / nFit << "\n";
    T << Form("\\newcommand{\\fit%sTfam%s}{%d}\n\\newcommand{\\fit%sNbins%s}{%d}\n", mp, cm[g], nTfam, mp, cm[g], nFit);
    for (int M = 0; M < 4; M++) T << Form("\\newcommand{\\fit%sChi%s%s}{%.1f}\n\\newcommand{\\fit%sBest%s%s}{%d}\n", mp, TString(kModel[M]).ReplaceAll("2", "two").Data(), cm[g], chiSum[M] / nFit, mp, TString(kModel[M]).ReplaceAll("2", "two").Data(), cm[g], nBest[M]);
  }
  fo->Close();
  std::ofstream TT(Form("tables/fit_%s_tab.tex", tag));
  TT << "% generated by scripts/make_fit.C -- do not edit\n";
  for (auto& r : tab) TT << r << "\n";
  L.close(); T.close(); TT.close(); printf("make_fit done\n");
}
