void chkq() {
  for (const char* f : {"output/mudst/ovq_test.root", "output/mudst/muq_test.root"}) {
    TFile* F = TFile::Open(f); if (!F) continue;
    for (const char* n : {"hQM2SE", "hQM2ME", "hQM2MEprodB", "hM2SE"}) { TH3* h = (TH3*)F->Get(n); if (h) { TH1* z = h->ProjectionZ(); printf("%s %s entries %.0f neg %.0f pos %.0f\n", f, n, h->GetEntries(), z->Integral(1, z->FindBin(-0.0001)), z->Integral(z->FindBin(0.0001), z->GetNbinsX())); } }
  }
}
