// root -b -q -l 'macros/drawTofMEQA.C("output/qa_test_10k.root","plots/qa_test_10k")'
// Uses the user's style.C / draw.C (copied into macros/).
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

void drawTofMEQA(const char* inName = "output/qa_test_10k.root", const char* outDir = "plots/qa_test_10k")
{
  gROOT->ProcessLine(".L macros/style.C");
  gROOT->ProcessLine("style();");
  gROOT->ProcessLine(".L macros/draw.C");
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetOptTitle(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile* f = TFile::Open(inName);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", inName); return; }
  const Double_t nEvt = ((TH1*)f->Get("hVz"))->GetEntries();
  TString tag = Form("Au+Au 19.6 GeV (Run19, P24iy), %.0f events", nEvt);

  TLatex lat; lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.045);

  // 1D histograms
  const char* h1[] = {"hEvtCut", "hTrigId", "hVzAll", "hVz", "hVzVpd", "hVzDiff", "hRefMult", "hGRefMult",
                      "hTofMult", "hNBTofMatch", "hNBTofHits", "hPt", "hEta", "hPhi", "hNHitsFit",
                      "hNHitsDedx", "hNHitsRatio", "hDca", "hMatchFlag", "hBeta", "hM2", "hYLocal", "hZLocal",
                      "hHitCell"};
  const int n1 = sizeof(h1) / sizeof(h1[0]);
  for (int i = 0; i < n1; i++) {
    TH1* h = (TH1*)f->Get(h1[i]);
    if (!h) { printf("missing %s\n", h1[i]); continue; }
    TCanvas* c = new TCanvas(Form("c_%s", h1[i]), "", 800, 600);
    Bool_t logy = TString(h1[i]).Contains("Mult") || TString(h1[i]).Contains("NBTof") || TString(h1[i]) == "hPt"
                  || TString(h1[i]) == "hDca" || TString(h1[i]) == "hM2" || TString(h1[i]) == "hBeta"
                  || TString(h1[i]) == "hVzDiff";
    c->SetLogy(logy);
    if (TString(h1[i]) == "hEvtCut" || TString(h1[i]) == "hTrigId") { c->SetBottomMargin(0.25); h->LabelsOption("v", "X"); }
    h->SetLineWidth(2);
    h->SetMinimum(logy ? 0.5 : 0);
    h->Draw("hist");
    lat.DrawLatex(0.22, 0.91, tag);
    c->SaveAs(Form("%s/%s.png", outDir, h1[i]));
    c->SaveAs(Form("%s/%s.pdf", outDir, h1[i]));
  }

  // 2D histograms
  const char* h2[] = {"hVxVy", "hRefMultVsTofMult", "hRefMultVsNBTofMatch", "hNBTofHitsVsNBTofMatch",
                      "hNTofTrkVsNBTofHits", "hEtaPhi", "hDedxVsPq", "hNSigPiVsP", "hInvBetaVsP", "hInvBetaVsPq",
                      "hM2VsP", "hM2VsPq", "hYLocalVsZLocal", "hYLocalVsP", "hZLocalVsP", "hTrayModMatched",
                      "hTrayModHits", "hHitPosPhiZ", "hNSigTofPiVsP", "hNSigTofKVsP", "hNSigTofPVsP",
                      "hDInvBetaPiVsP", "hEtaPhiTof"};
  const int n2 = sizeof(h2) / sizeof(h2[0]);
  for (int i = 0; i < n2; i++) {
    TH2* h = (TH2*)f->Get(h2[i]);
    if (!h) { printf("missing %s\n", h2[i]); continue; }
    TCanvas* c = new TCanvas(Form("c_%s", h2[i]), "", 800, 600);
    c->SetRightMargin(0.14);
    c->SetLogz(!TString(h2[i]).BeginsWith("hTrayMod"));
    h->Draw("colz");
    lat.DrawLatex(0.22, 0.91, tag);
    c->SaveAs(Form("%s/%s.png", outDir, h2[i]));
    c->SaveAs(Form("%s/%s.pdf", outDir, h2[i]));
  }

  // TOF matching efficiency (matched w/ beta>0 / TPC primary tracks)
  const char* eff[][2] = {{"hPtTof", "hPtAll"}, {"hEtaTof", "hEtaAll"}, {"hPhiTof", "hPhiAll"}};
  const char* effName[] = {"hEffPt", "hEffEta", "hEffPhi"};
  for (int i = 0; i < 3; i++) {
    TH1* num = (TH1*)f->Get(eff[i][0]);
    TH1* den = (TH1*)f->Get(eff[i][1]);
    if (!num || !den) continue;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(num, den, "cl=0.683 b(1,1) mode");
    TCanvas* c = new TCanvas(Form("c_%s", effName[i]), "", 800, 600);
    TH1* fr = (TH1*)den->Clone(Form("fr%d", i)); fr->Reset();
    fr->SetMinimum(0); fr->SetMaximum(1.0);
    fr->GetYaxis()->SetTitle("TOF matching efficiency");
    fr->Draw();
    g->SetMarkerStyle(20); g->SetMarkerSize(i == 2 ? 0.6 : 1.0);
    g->SetLineColor(1); g->SetMarkerColor(1);
    if (i == 0) fr->GetXaxis()->SetRangeUser(0, 3.);
    g->Draw("p same");
    lat.DrawLatex(0.22, 0.91, tag);
    c->SaveAs(Form("%s/%s.png", outDir, effName[i]));
    c->SaveAs(Form("%s/%s.pdf", outDir, effName[i]));
  }

  // summary numbers
  TH1* hc = (TH1*)f->Get("hEvtCut");
  for (int i = 1; i <= hc->GetNbinsX(); i++) printf("cut %-24s %10.0f\n", hc->GetXaxis()->GetBinLabel(i), hc->GetBinContent(i));
  TH1* pa = (TH1*)f->Get("hPtAll");
  TH1* pt = (TH1*)f->Get("hPtTof");
  printf("primary tracks (QA cuts): %.0f   TOF matched (beta>0): %.0f   ratio %.3f\n",
         pa->GetEntries(), pt->GetEntries(), pt->GetEntries() / pa->GetEntries());
  int b1 = pa->FindBin(0.5), b2 = pa->FindBin(1.5);
  printf("TOF match eff, 0.5<pT<1.5: %.3f\n", pt->Integral(b1, b2) / pa->Integral(b1, b2));
  printf("mean refMult %.1f, mean nBTOFMatch %.1f, mean N BTofHit %.1f\n",
         ((TH1*)f->Get("hRefMult"))->GetMean(), ((TH1*)f->Get("hNBTofMatch"))->GetMean(),
         ((TH1*)f->Get("hNBTofHits"))->GetMean());
}
