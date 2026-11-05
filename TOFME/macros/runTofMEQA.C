// root4star -b -q 'macros/runTofMEQA.C(10000,"filelists/test.list","output/qa_test.root")'
void runTofMEQA(Long64_t nEvents = 10000,
                const char* inputList = "filelists/test.list",
                const char* outName = "output/qa_test.root",
                Bool_t useMBTrig = kTRUE)
{
  gROOT->Macro("loadMuDst.C");   // STAR base libs
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StTofMEQAMaker");

  StChain* chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, inputList, "picoDst");
  // read only the branches needed (much faster over XRootD)
  picoMaker->SetStatus("*", 0);
  const char* br[] = {"^Event$", "^Track$", "^BTofHit$", "^BTofPidTraits$"};
  for (int ib = 0; ib < 4; ib++) picoMaker->SetStatus(br[ib], 1);
  StTofMEQAMaker* qa = new StTofMEQAMaker("TofMEQA", picoMaker, outName);
  qa->setVzCut(70.);
  qa->setVrCut(2.);
  if (useMBTrig) {
    // Run19 19.6 GeV minimum-bias trigger ids
    const UInt_t mb[] = {640001, 640011, 640021, 640031, 640041, 640051};
    for (int it = 0; it < 6; it++) qa->addTriggerId(mb[it]);
  }

  if (chain->Init() == kStERR) { cout << "chain Init failed" << endl; return; }
  Long64_t nAvail = picoMaker->chain()->GetEntries();
  Long64_t nRun = (nEvents > 0 && nEvents < nAvail) ? nEvents : nAvail;
  cout << "runTofMEQA: " << nAvail << " events in chain, processing " << nRun << endl;

  for (Long64_t iev = 0; iev < nRun; iev++) {
    if (iev % 1000 == 0) cout << "  event " << iev << endl;
    chain->Clear();
    Int_t ret = chain->Make(iev);
    if (ret == kStEOF || ret == kStFatal) { cout << "Make returned " << ret << " at " << iev << endl; break; }
  }
  chain->Finish();
  delete chain;
}
