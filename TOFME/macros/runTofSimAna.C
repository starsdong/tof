// runTofSimAna.C : truth decomposition of BTOF matches in simulation MuDsts (StTofSimAnaMaker)
// usage: root4star -b -q 'macros/runTofSimAna.C(N, "list or file", "out.root")'
void runTofSimAna(Int_t nEv = 1000000, const char* in = "test.list", const char* out = "simana.root")
{
  gROOT->Macro("loadMuDst.C");
  gSystem->Load("StTofSimAnaMaker");
  StChain* chain = new StChain();
  StMuDstMaker* mu = new StMuDstMaker(0, 0, "", in, "", 100000);
  StTofSimAnaMaker* ana = new StTofSimAnaMaker("simAna", mu, out);
  chain->Init();
  const Int_t n = std::min((Long64_t)nEv, mu->chain()->GetEntries());
  for (Int_t i = 0; i < n; i++) { chain->Clear(); if (chain->Make(i) >= kStEOF) break; }
  chain->Finish();
  printf("runTofSimAna: done %d events\n", n);
}
