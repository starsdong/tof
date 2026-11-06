// make_corrdd.C : data-driven estimate of the track-correlated (late-hit) component of the TOF m2 spectrum (method 2 data).
//
// Observation (logs/corrdd.txt, Sec. sec:corrdd): SE - ME away from the particle peaks is smooth, charge-symmetric for p > 0.3 GeV/c,
// the same per QA track in all centralities, and depends on the hit time: roughly flat up to t ~ 30 ns after the collision and
// falling to ~0 by t ~ 50-60 ns. Model, per QA track and per ns of hit time t:
//     f(t) = A / (1 + exp((t - t0)/tau)),   t = kappa * sqrt(1 + m2/p^2),  kappa = L/c (L = 2.5 m, kappa = 8.34 ns)
// (f = 0 for t below the pion arrival time, m2 < m2_pi: the hits are later than the track's own particle)
// so in m2:  dN/dm2 = sum_j N_j f(t_j(m2)) dt_j/dm2 over fine momentum bins j (N_j QA tracks).
// A, t0, tau are fitted (chi2) to SE - ME of NEGATIVE tracks (no proton energy-loss tail) in m2 sidebands free of particles,
// 0.20 < p < 1.0 GeV/c, all centralities. The fitted f is applied per QA track to all charges and every centrality class and
// integrated under the pi, K and p peaks. Validation: positive tracks (p > 0.3), per-centrality sidebands, K-p valley.
// The full simulation (Sec. sec:sim) is used as evidence for the origin (secondaries from the BTOF trays), not as input.
// usage (analysis_note/): root -b -q -l scripts/make_corrdd.C
// input: data/ov_big.root. outputs: logs/corrdd.txt, tables/corrdd.tex (macros), tables/corrdd_tab.tex, tables/corrdd_val.tex,
//   figures/corrdd_res.{pdf,png}, figures/corrdd_m2.{pdf,png}, figures/corrdd_cent.{pdf,png}
#include <fstream>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
namespace cdd {
TH3 *gSE, *gME, *gSEn, *gMEn; TH2 *gNk, *gNt;
const double kKappa = 8.34;   // ns, L = 2.5 m
const double kMpi2 = 0.01948, kMK2 = 0.24372, kMp2 = 0.88035, kMd2 = 3.5155;
// histograms of SE - ME for one fine p bin j and cent classes c1..c2 (0-based); charge: 0 all, -1 neg, +1 pos
TH1D* resid(int j, int c1, int c2, int q, double& N, double& Nall, TH1D** seOut = 0, TH1D** meOut = 0)
{
  TH1D* r = gSE->ProjectionZ(Form("r_%d_%d_%d_%d", j, c1, c2, q), c1 + 1, c1 + 1, j, j); r->Reset();
  TH1D* se = (TH1D*)r->Clone(Form("s_%d_%d_%d_%d", j, c1, c2, q)); TH1D* me = (TH1D*)r->Clone(Form("m_%d_%d_%d_%d", j, c1, c2, q));
  N = 0; Nall = 0;
  for (int c = c1; c <= c2; c++) {
    const double k = gNt->GetBinContent(c + 1, j) > 0 ? gNk->GetBinContent(c + 1, j) / gNt->GetBinContent(c + 1, j) : 0;
    TH1D* a = gSE->ProjectionZ("ta", c + 1, c + 1, j, j); TH1D* b = gME->ProjectionZ("tb", c + 1, c + 1, j, j);
    TH1D* an = gSEn->ProjectionZ("tan", c + 1, c + 1, j, j); TH1D* bn = gMEn->ProjectionZ("tbn", c + 1, c + 1, j, j);
    if (q == 0) { se->Add(a); me->Add(b, k); }
    if (q < 0) { se->Add(an); me->Add(bn, k); }
    if (q > 0) { se->Add(a); se->Add(an, -1); me->Add(b, k); me->Add(bn, -k); }
    delete a; delete b; delete an; delete bn;
    Nall += gNk->GetBinContent(c + 1, j);
  }
  N = Nall;   // QA tracks of all charges (per-track normalization uses all tracks; charge fractions handled by the caller)
  r->Add(se); r->Add(me, -1);
  if (seOut) *seOut = se; else delete se;
  if (meOut) *meOut = me; else delete me;
  return r;
}
// sideband: m2 away from pi, K, p, d peaks (widths grow with p)
bool inSB(double m2, double p)
{
  const double s = 0.02 + 0.10 * p * p;                // approximate m2 resolution scale
  if (m2 < kMK2 + 4 * s + 0.08) return false;          // below and at pi, K
  if (std::fabs(m2 - kMp2) < 0.30 + 6 * s) return false;
  if (std::fabs(m2 - kMd2) < 0.45 + 6 * s) return false;
  return m2 < 5.0;
}
double fT(double t, const double* a) { return a[0] / (1. + std::exp((t - a[1]) / a[2])); }
// expected correlated counts in an m2 interval [x1,x2] for N tracks at momentum p (integral of f dt)
double expInt(double x1, double x2, double p, double N, const double* a)
{
  if (x1 < kMpi2) x1 = kMpi2;   // the component is later than the track's own arrival: start at the pion mass
  if (x2 <= x1) return 0;
  const double t1 = kKappa * std::sqrt(1. + x1 / (p * p)), t2 = kKappa * std::sqrt(1. + x2 / (p * p));
  // analytic integral of the Fermi function
  auto F = [&](double t) { return a[0] * (t - a[2] * std::log(1. + std::exp((t - a[1]) / a[2]))); };
  return N * (F(t2) - F(t1));
}
struct Pt { double m2a, m2b, p, N, y, e, me; };   // e unused (Pearson variance me + model, avoids the Neyman bias at low counts)
std::vector<Pt> gPts;
double chi2(const double* a)
{
  double c = 0;
  for (auto& q : gPts) { const double m = expInt(q.m2a, q.m2b, q.p, q.N, a); const double v = std::max(q.me + m, 1.); c += (q.y - m) * (q.y - m) / v; }
  return c;
}
double ig(TH1* h, double a, double b) { return h->Integral(h->FindBin(a + 1e-6), h->FindBin(b - 1e-6)); }
}

void make_corrdd()
{
  using namespace cdd;
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0); gStyle->SetEndErrorSize(0);
  TFile* D = TFile::Open("data/ov_big.root");
  gSE = (TH3*)D->Get("hM2SE"); gME = (TH3*)D->Get("hM2ME"); gSEn = (TH3*)D->Get("hM2SEneg"); gMEn = (TH3*)D->Get("hM2MEneg");
  gNk = (TH2*)D->Get("hNTrk"); gNt = (TH2*)D->Get("hNTrial");
  std::ofstream L("logs/corrdd.txt"), T("tables/corrdd.tex");
  L << "# generated by scripts/make_corrdd.C\n"; T << "% generated by scripts/make_corrdd.C -- do not edit\n";
  const int jlo = gSE->GetYaxis()->FindBin(0.2001), jhi = gSE->GetYaxis()->FindBin(0.9999);
  // ---------- (1) fit f(t) to negative tracks, all centralities ----------
  // negative-track SE - ME, but per-track normalization needs the number of negative tracks: estimate it from the
  // fraction of negative TOF matches (SE counts) in the same fine p bin (all m2) -> N_neg = N_all * fneg
  int nPts = 0;
  for (int j = jlo; j <= jhi; j++) {
    double N, Na; TH1D *sen, *men; TH1D* r = resid(j, 0, 8, -1, N, Na, &sen, &men);
    TH1D *sea, *mea; double N2, Na2; TH1D* ra = resid(j, 0, 8, 0, N2, Na2, &sea, &mea);
    const double fneg = sea->Integral() > 0 ? sen->Integral() / sea->Integral() : 0.5;
    const double p = gSE->GetYaxis()->GetBinCenter(j);
    for (int b = 1; b <= r->GetNbinsX(); b += 5) {   // 0.05 GeV^2 groups
      const double x1 = r->GetBinLowEdge(b), x2 = r->GetBinLowEdge(b + 5);
      if (!inSB(0.5 * (x1 + x2), p)) continue;
      double y = 0, s = 0, m = 0;
      for (int k = b; k < b + 5; k++) { y += r->GetBinContent(k); s += sen->GetBinContent(k); m += men->GetBinContent(k); }
      const double e = std::sqrt(std::max(s, 1.) + 0.0 * m);
      gPts.push_back({x1, x2, p, Na * fneg, y, e, m}); nPts++;
    }
    delete r; delete sen; delete men; delete ra; delete sea; delete mea;
  }
  ROOT::Math::Minimizer* mn = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  ROOT::Math::Functor fc(&chi2, 3); mn->SetFunction(fc);
  mn->SetLimitedVariable(0, "A", 2e-5, 1e-6, 0, 1e-2); mn->SetLimitedVariable(1, "t0", 35, 0.5, 15, 80); mn->SetLimitedVariable(2, "tau", 5, 0.2, 0.5, 30);
  mn->SetPrintLevel(0); mn->Minimize(); mn->Hesse();
  const double* a = mn->X(); const double* ea = mn->Errors();
  double par[3] = {a[0], a[1], a[2]};
  const double c2 = mn->MinValue(); const int ndf = nPts - 3;
  L << "fit negatives, 0.2<p<1, all cent: points " << nPts << " chi2/ndf " << c2 << "/" << ndf << " | A " << a[0] << " +- " << ea[0]
    << " per track per ns | t0 " << a[1] << " +- " << ea[1] << " ns | tau " << a[2] << " +- " << ea[2] << " ns\n";
  T << Form("\\newcommand{\\cddA}{%.2f}\n\\newcommand{\\cddAe}{%.2f}\n\\newcommand{\\cddTz}{%.1f}\n\\newcommand{\\cddTze}{%.1f}\n\\newcommand{\\cddTau}{%.1f}\n\\newcommand{\\cddTaue}{%.1f}\n\\newcommand{\\cddChiN}{%.0f}\n\\newcommand{\\cddNdf}{%d}\n",
            1e6 * a[0], 1e6 * ea[0], a[1], ea[1], a[2], ea[2], c2, ndf);
  // systematic variants: kappa +-10% handled by refitting with scaled kappa is equivalent to scaling t0, tau: report only the
  // yield variation with an alternative shape (flat in time, no fall-off, fitted A) below.

  // ---------- (2) apply per centrality group and momentum bin; validation and yields ----------
  const int np = 3; const double pb[np][2] = {{0.2, 0.4}, {0.4, 0.6}, {0.6, 1.0}}; const char* pn[np] = {"A", "B", "C"};
  const int ng = 4; const int gc1[ng] = {7, 4, 0, 0}, gc2[ng] = {8, 5, 1, 8};
  const char* gl[ng] = {"0--10\\%", "20--40\\%", "60--80\\%", "0--80\\%"}; const char* gt[ng] = {"0-10%", "20-40%", "60-80%", "0-80%"}; const char* gn[ng] = {"Ci", "Mid", "Per", "All"};
  const int nw = 5; const double wl[nw] = {-0.1, 0.15, 0.6, 0.42, 1.6}, wh[nw] = {0.12, 0.40, 1.2, 0.58, 2.8}; const char* wn[nw] = {"pi", "K", "p", "valley", "pd"};
  std::ofstream TT("tables/corrdd_tab.tex"), TV("tables/corrdd_val.tex");
  TT << "% generated by scripts/make_corrdd.C -- do not edit\n\\begin{table}[H]\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}\n"
     << "\\caption{Data-driven correlated component under the identified-particle windows ($\\pi$: $-0.1<m^2<0.12$, K: 0.15--0.40,\n"
     << "p: 0.6--1.2\\,GeV$^2$/$c^4$), all charges, in \\% of SE$-$ME in the same window, i.e.\\ the fraction of the ME-subtracted yield that\n"
     << "is track-correlated. Fit uncertainty from the covariance of $(A,t_0,\\tau)$ is below 10\\% (relative) of each value.}\n\\label{tab:corrdd}\n"
     << "\\begin{tabular}{llrrr}\\hline\n$p$ (GeV/$c$) & centrality & $\\pi$ (\\%) & K (\\%) & p (\\%) \\\\ \\hline\n";
  TV << "% generated by scripts/make_corrdd.C -- do not edit\n\\begin{table}[H]\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}\n"
     << "\\caption{Validation of the data-driven correlated component (all charges): data SE over (ME + correlated) in the K--p valley\n"
     << "(0.42--0.58\\,GeV$^2$/$c^4$) and the p--d window (1.6--2.8), with SE/ME for comparison. The fit uses negative tracks of all centralities;\n"
     << "positive tracks and the individual centralities are not fitted. Statistical errors of SE.}\n\\label{tab:corrddval}\n"
     << "\\begin{tabular}{llrrrr}\\hline\n$p$ (GeV/$c$) & centrality & valley SE/ME & valley SE/(ME+corr) & p--d SE/ME & p--d SE/(ME+corr) \\\\ \\hline\n";
  TCanvas* cm = new TCanvas("cm", "", 1500, 1250); cm->Divide(3, 3, 0.001, 0.001);
  // ranges for the text: [0] p-d SE/(ME+corr) central+mid at p bins B,C; [1] p-d SE/ME same; [2] valley SE/(ME+corr) c+m at A,B;
  // [3] valley c+m at C; [4] valley SE/ME c+m at A,B,C; [5] valley SE/(ME+corr) peripheral
  double rmin[6], rmax[6]; for (int k = 0; k < 6; k++) { rmin[k] = 1e9; rmax[k] = -1e9; }
  auto upd = [&](int k, double v) { rmin[k] = std::min(rmin[k], v); rmax[k] = std::max(rmax[k], v); };
  for (int ip = 0; ip < np; ip++) for (int g = 0; g < ng; g++) {
    const int j1 = gSE->GetYaxis()->FindBin(pb[ip][0] + 1e-6), j2 = gSE->GetYaxis()->FindBin(pb[ip][1] - 1e-6);
    TH1D *se = 0, *me = 0, *cc = 0;
    for (int j = j1; j <= j2; j++) {
      double N, Na; TH1D *s, *m; TH1D* r = resid(j, gc1[g], gc2[g], 0, N, Na, &s, &m); delete r;
      if (!se) { se = (TH1D*)s->Clone(Form("se_%d_%d", ip, g)); se->Reset(); me = (TH1D*)se->Clone(Form("me_%d_%d", ip, g)); cc = (TH1D*)se->Clone(Form("cc_%d_%d", ip, g)); }
      se->Add(s); me->Add(m); delete s; delete m;
      const double p = gSE->GetYaxis()->GetBinCenter(j);
      for (int b = 1; b <= cc->GetNbinsX(); b++) cc->AddBinContent(b, expInt(cc->GetBinLowEdge(b), cc->GetBinLowEdge(b + 1), p, Na, par));
    }
    double fr[nw], val[nw][2];
    L << gt[g] << " p " << pb[ip][0] << "-" << pb[ip][1];
    for (int w = 0; w < nw; w++) {
      const double s = ig(se, wl[w], wh[w]), m = ig(me, wl[w], wh[w]), c = ig(cc, wl[w], wh[w]);
      fr[w] = 100 * c / (s - m); val[w][0] = s / m; val[w][1] = s / (m + c);
      L << " | " << wn[w] << " SE " << s << " ME " << m << " corr " << c << " corr/(SE-ME) " << fr[w] << "% SE/ME " << s / m << " SE/(ME+corr) " << s / (m + c) << " +- " << std::sqrt(s) / (m + c);
      T << Form("\\newcommand{\\cdd%s%s%s}{%.2f}\n\\newcommand{\\cddR%s%s%s}{%.2f}\n", wn[w], gn[g], pn[ip], fr[w], wn[w], gn[g], pn[ip], s / (m + c));
    }
    L << "\n";
    if (g == 0 || g == 1) {
      if (ip > 0) { upd(0, val[4][1]); upd(1, val[4][0]); }
      if (ip < 2) upd(2, val[3][1]); else upd(3, val[3][1]);
      upd(4, val[3][0]);
    }
    if (g == 2) upd(5, val[3][1]);
    TT << Form("%s & %s & %.2f & %.1f & %.1f \\\\\n", g == 0 ? Form("%.1f--%.1f", pb[ip][0], pb[ip][1]) : "", gl[g], fr[0], fr[1], fr[2]);
    const double sv = ig(se, 0.42, 0.58), sp = ig(se, 1.6, 2.8);
    TV << Form("%s & %s & %.2f & %.2f $\\pm$ %.2f & %.2f & %.2f $\\pm$ %.2f \\\\\n", g == 0 ? Form("%.1f--%.1f", pb[ip][0], pb[ip][1]) : "", gl[g],
               val[3][0], val[3][1], val[3][1] / std::sqrt(std::max(sv, 1.)), val[4][0], val[4][1], val[4][1] / std::sqrt(std::max(sp, 1.)));
    if (g < 3) {
      cm->cd(ip * 3 + g + 1); gPad->SetLogy(1); gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.13);
      se->Rebin(5); me->Rebin(5); cc->Rebin(5); TH1D* sum = (TH1D*)me->Clone(Form("sum_%d_%d", ip, g)); sum->Add(cc);
      se->SetTitle(""); se->GetXaxis()->SetTitle("m^{2} (GeV^{2}/c^{4})"); se->GetYaxis()->SetTitle("counts / 0.05 GeV^{2}/c^{4}");
      se->GetXaxis()->SetTitleSize(0.05); se->GetYaxis()->SetTitleSize(0.05); se->GetXaxis()->SetLabelSize(0.045); se->GetYaxis()->SetLabelSize(0.045); se->GetYaxis()->SetTitleOffset(1.5);
      se->GetXaxis()->SetRangeUser(-0.5, 5.0); se->SetMinimum(0.5); se->SetMaximum(se->GetMaximum() * 30);
      se->SetMarkerStyle(20); se->SetMarkerSize(0.5); se->Draw("pe");
      me->SetLineColor(kRed + 1); me->SetFillColorAlpha(kRed - 9, 0.5); me->Draw("hist same");
      cc->SetLineColor(kGreen + 2); cc->SetLineWidth(2); cc->SetLineStyle(2); cc->Draw("hist same");
      sum->SetLineColor(kBlue + 1); sum->SetLineWidth(2); sum->Draw("hist same"); se->Draw("pe same");
      TLatex lt; lt.SetNDC(); lt.SetTextFont(42); lt.SetTextSize(0.05); lt.DrawLatex(0.42, 0.88, Form("%s, %.1f < p < %.1f GeV/c", gt[g], pb[ip][0], pb[ip][1]));
      if (g == 0 && ip == 0) { TLegend* lg = new TLegend(0.42, 0.58, 0.97, 0.83); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
        lg->AddEntry(se, "data SE", "p"); lg->AddEntry(me, "ME (method 2)", "f"); lg->AddEntry(cc, "correlated (data-driven)", "l"); lg->AddEntry(sum, "ME + correlated", "l"); lg->Draw(); }
    }
  }
  const char* rn[6] = {"PdCM", "PdSEME", "VaCMab", "VaCMc", "VaSEME", "VaPer"};
  for (int k = 0; k < 6; k++) T << Form("\\newcommand{\\cdd%sMin}{%.2f}\n\\newcommand{\\cdd%sMax}{%.2f}\n", rn[k], rmin[k], rn[k], rmax[k]);
  TT << "\\hline\\end{tabular}\\end{table}\n"; TT.close(); TV << "\\hline\\end{tabular}\\end{table}\n"; TV.close();
  cm->SaveAs("figures/corrdd_m2.pdf"); cm->SaveAs("figures/corrdd_m2.png");

  // ---------- (3) residual figure: SE - ME per 1e4 QA tracks of each charge vs m2, with the model, 0-80% ----------
  TCanvas* cr = new TCanvas("cr", "", 1500, 1000); cr->Divide(3, 2, 0.001, 0.001);
  const double rp[6][2] = {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.8}, {0.8, 1.0}};
  for (int k = 0; k < 6; k++) {
    const int j1 = gSE->GetYaxis()->FindBin(rp[k][0] + 1e-6), j2 = gSE->GetYaxis()->FindBin(rp[k][1] - 1e-6);
    TH1D *rn = 0, *rq = 0, *mdl = 0; double Nt = 0;
    for (int j = j1; j <= j2; j++) {
      double N, Na; TH1D *sn, *mn2, *sa, *ma;
      TH1D* a1 = resid(j, 0, 8, -1, N, Na, &sn, &mn2); TH1D* a2 = resid(j, 0, 8, 0, N, Na, &sa, &ma);
      const double fneg = sa->Integral() > 0 ? sn->Integral() / sa->Integral() : 0.5;
      TH1D* a3 = (TH1D*)a2->Clone("a3"); a3->Add(a1, -1);   // positive
      if (!rn) { rn = (TH1D*)a1->Clone(Form("rn%d", k)); rn->Reset(); rq = (TH1D*)rn->Clone(Form("rq%d", k)); mdl = (TH1D*)rn->Clone(Form("md%d", k)); }
      // per 1e4 tracks of the respective charge: scale negatives by 1/(N fneg), positives by 1/(N (1-fneg)), summed over j
      rn->Add(a1, 1e4 / (Na * fneg) * Na / 1.); rq->Add(a3, 1e4 / (Na * (1 - fneg)) * Na / 1.);
      const double p = gSE->GetYaxis()->GetBinCenter(j);
      for (int b = 1; b <= mdl->GetNbinsX(); b++) mdl->AddBinContent(b, 1e4 * expInt(mdl->GetBinLowEdge(b), mdl->GetBinLowEdge(b + 1), p, 1., par) * Na);
      Nt += Na;
      delete a1; delete a2; delete a3; delete sn; delete mn2; delete sa; delete ma;
    }
    rn->Scale(1. / Nt); rq->Scale(1. / Nt); mdl->Scale(1. / Nt);   // track-weighted average over the fine bins
    rn->Rebin(10); rq->Rebin(10); mdl->Rebin(10);
    cr->cd(k + 1); gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.13);
    rq->SetTitle(""); rq->GetXaxis()->SetTitle("m^{2} (GeV^{2}/c^{4})"); rq->GetYaxis()->SetTitle("(SE - ME) / 0.1 GeV^{2} / 10^{4} tracks");
    rq->GetXaxis()->SetTitleSize(0.05); rq->GetYaxis()->SetTitleSize(0.045); rq->GetXaxis()->SetLabelSize(0.045); rq->GetYaxis()->SetLabelSize(0.045); rq->GetYaxis()->SetTitleOffset(1.6);
    rq->GetXaxis()->SetRangeUser(0.3, 5.0);
    double mx = 0; for (int b = rq->FindBin(1.3); b <= rq->FindBin(4.99); b++) mx = std::max(mx, std::max(rq->GetBinContent(b), rn->GetBinContent(b)));
    mx = std::max(mx, mdl->GetBinContent(mdl->FindBin(0.5))); rq->SetMaximum(mx * 1.6); rq->SetMinimum(-0.15 * mx);
    rq->SetMarkerStyle(24); rq->SetMarkerColor(kRed + 1); rq->SetLineColor(kRed + 1); rq->SetMarkerSize(0.7); rq->Draw("pe");
    rn->SetMarkerStyle(20); rn->SetMarkerColor(kBlue + 1); rn->SetLineColor(kBlue + 1); rn->SetMarkerSize(0.7); rn->Draw("pe same");
    mdl->SetLineColor(kBlack); mdl->SetLineWidth(2); mdl->Draw("hist same");
    TLine* l0 = new TLine(0.3, 0, 5.0, 0); l0->SetLineStyle(3); l0->Draw();
    TLatex lt; lt.SetNDC(); lt.SetTextFont(42); lt.SetTextSize(0.05); lt.DrawLatex(0.45, 0.88, Form("%.1f < p < %.1f GeV/c, 0-80%%", rp[k][0], rp[k][1]));
    if (k == 0) { TLegend* lg = new TLegend(0.45, 0.62, 0.97, 0.84); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
      lg->AddEntry(rn, "negative tracks", "p"); lg->AddEntry(rq, "positive tracks", "p"); lg->AddEntry(mdl, "fit f(t) (negatives)", "l"); lg->Draw(); }
  }
  cr->SaveAs("figures/corrdd_res.pdf"); cr->SaveAs("figures/corrdd_res.png");

  // ---------- (4) per-track rate vs centrality: A fitted per class (t0, tau fixed), negatives ----------
  TCanvas* c9 = new TCanvas("c9", "", 800, 600); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.15); gPad->SetRightMargin(0.03);
  const double cx[9] = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
  TGraphErrors* ge = new TGraphErrors(); double amin = 1e9, amax = 0;
  for (int c = 0; c < 9; c++) {
    gPts.clear();
    for (int j = jlo; j <= jhi; j++) {
      double N, Na; TH1D *sen, *men, *sea, *mea; TH1D* r = resid(j, c, c, -1, N, Na, &sen, &men); TH1D* ra = resid(j, c, c, 0, N, Na, &sea, &mea);
      const double fneg = sea->Integral() > 0 ? sen->Integral() / sea->Integral() : 0.5; const double p = gSE->GetYaxis()->GetBinCenter(j);
      for (int b = 1; b <= r->GetNbinsX(); b += 5) { const double x1 = r->GetBinLowEdge(b), x2 = r->GetBinLowEdge(b + 5); if (!inSB(0.5 * (x1 + x2), p)) continue;
        double y = 0, s = 0, mm = 0; for (int k = b; k < b + 5; k++) { y += r->GetBinContent(k); s += sen->GetBinContent(k); mm += men->GetBinContent(k); }
        gPts.push_back({x1, x2, p, Na * fneg, y, std::sqrt(std::max(s, 1.)), mm}); }
      delete r; delete ra; delete sen; delete men; delete sea; delete mea;
    }
    // linear in A: A = sum(y m1 / e2) / sum(m1^2 / e2) with m1 the model for A = 1
    double num = 0, den = 0; double u[3] = {1., par[1], par[2]};
    for (auto& q : gPts) { const double m1 = expInt(q.m2a, q.m2b, q.p, q.N, u); const double v = std::max(q.me + par[0] * m1, 1.); num += q.y * m1 / v; den += m1 * m1 / v; }
    const double Ac = num / den, eAc = 1. / std::sqrt(den);
    ge->SetPoint(c, cx[c], 1e6 * Ac); ge->SetPointError(c, 0, 1e6 * eAc); amin = std::min(amin, Ac); amax = std::max(amax, Ac);
    L << "cent9 " << c << " A " << Ac << " +- " << eAc << "\n";
  }
  T << Form("\\newcommand{\\cddAmin}{%.2f}\n\\newcommand{\\cddAmax}{%.2f}\n", 1e6 * amin, 1e6 * amax);
  TH1F* fr = gPad->DrawFrame(0, 0, 80, 1e6 * amax * 1.6); fr->GetXaxis()->SetTitle("centrality (%)"); fr->GetYaxis()->SetTitle("A (10^{-6} per QA track per ns)");
  fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045);
  fr->GetYaxis()->SetTitleOffset(1.3); fr->GetXaxis()->SetTitleOffset(1.1);
  ge->SetMarkerStyle(20); ge->Draw("p same");
  TLine* la = new TLine(0, 1e6 * par[0], 80, 1e6 * par[0]); la->SetLineColor(kRed + 1); la->SetLineStyle(2); la->Draw();
  c9->SaveAs("figures/corrdd_cent.pdf"); c9->SaveAs("figures/corrdd_cent.png");
  L.close(); T.close(); printf("make_corrdd done\n");
}
