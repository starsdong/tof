// reco_sim.C : full STAR reconstruction of a starsim .fzd file (TpcRS + Sti + BTOF/VPD simulators + BTOF matching/calibration)
// usage: root4star -b -q 'reco_sim.C(nevents, "in.fzd", "chainopts")'
void reco_sim(int nevents, const char* fzd, const char* opts)
{
  gROOT->LoadMacro("bfc.C");
  bfc(-1, opts, fzd);
  chain->Init();
  for (int i = 0; i < nevents; i++) {
    chain->Clear();
    int st = chain->Make(i + 1);
    if (st >= kStEOF) { printf("reco_sim: stop at event %d status %d\n", i, st); break; }
    printf("reco_sim: event %d done\n", i);
  }
  chain->Finish();
  printf("reco_sim: finished\n");
}
