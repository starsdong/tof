void dumpEvtCut(const char* fn = "output/qa_test_10k.root")
{
  TFile* f = TFile::Open(fn);
  TH1* h = (TH1*)f->Get("hEvtCut");
  printf("entries %.0f nbins %d\n", h->GetEntries(), h->GetNbinsX());
  for (int i = 0; i <= h->GetNbinsX() + 1; i++)
    printf("bin %d [%s] %.0f\n", i, h->GetXaxis()->GetBinLabel(i), h->GetBinContent(i));
  TH1* v = (TH1*)f->Get("hVzAll");
  printf("hVzAll entries %.0f\n", v->GetEntries());
}
