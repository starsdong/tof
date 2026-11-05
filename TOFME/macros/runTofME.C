// root4star -b -q 'macros/runTofME.C(10000,"filelists/test.list","output/me_test.root","geom/btof_geom_test.txt",5)'
// geomFile "" -> only accumulate the geometry sums (hGeomSum), no projection/mixing.
void runTofME(Long64_t nEvents = 10000,
              const char* inputList = "filelists/test.list",
              const char* outName = "output/me_test.root",
              const char* geomFile = "geom/btof_geom_test.txt",
              Int_t nMix = 5,
              const char* pdfFile = "")
{
  gROOT->Macro("loadMuDst.C");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StTofMEMaker");

  StChain* chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, inputList, "picoDst");
  picoMaker->SetStatus("*", 0);
  const char* br[] = {"^Event$", "^Track$", "^BTofHit$", "^BTofPidTraits$"};
  for (int ib = 0; ib < 4; ib++) picoMaker->SetStatus(br[ib], 1);

  StTofMEMaker* me = new StTofMEMaker("TofME", picoMaker, outName);
  me->setGeomFile(geomFile);
  me->setNMix(nMix);
  me->setTimePdfFile(pdfFile);
  me->setVzCut(70.);
  me->setVrCut(2.);
  me->setZLocalCut(3.05);
  const UInt_t mb[] = {640001, 640011, 640021, 640031, 640041, 640051};
  for (int it = 0; it < 6; it++) me->addTriggerId(mb[it]);

  if (chain->Init() == kStERR) { cout << "chain Init failed" << endl; return; }
  Long64_t nAvail = picoMaker->chain()->GetEntries();
  Long64_t nRun = (nEvents > 0 && nEvents < nAvail) ? nEvents : nAvail;
  cout << "runTofME: " << nAvail << " events in chain, processing " << nRun << endl;
  for (Long64_t iev = 0; iev < nRun; iev++) {
    if (iev % 2000 == 0) cout << "  event " << iev << endl;
    chain->Clear();
    Int_t ret = chain->Make(iev);
    if (ret == kStEOF || ret == kStFatal) { cout << "Make returned " << ret << " at " << iev << endl; break; }
  }
  chain->Finish();
  delete chain;
}
