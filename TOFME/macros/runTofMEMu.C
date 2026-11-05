// root4star -b -q 'macros/runTofMEMu.C(-1,"filelists/mudst/run_20060065.list","output/mudst/me_20060065.root","geom/btof_geom_v3.txt",5)'
void runTofMEMu(Long64_t nEvents = -1, const char* inputList = "filelists/mudst/test1.list",
                const char* outName = "output/mudst/me_test.root", const char* geomFile = "geom/btof_geom_v3.txt",
                Int_t nMix = 5, const char* offPattern = "output/mudst/off_%d.root")
{
  gROOT->Macro("loadMuDst.C");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StTofMEMaker");
  gSystem->Load("StTofMEMuMaker");
  StChain* chain = new StChain();
  StMuDstMaker* mu = new StMuDstMaker(0, 0, "", inputList, "MuDst", 100);
  mu->SetStatus("*", 0);
  const char* br[] = {"MuEvent*", "PrimaryVertices*", "PrimaryTracks*", "GlobalTracks*", "BTof*"};
  for (int i = 0; i < 5; i++) mu->SetStatus(br[i], 1);
  StTofMEMuMaker* me = new StTofMEMuMaker("TofMEMu", mu, outName);
  me->setGeomFile(geomFile);
  me->setNMix(nMix);
  me->setOffsetPattern(offPattern);
  me->setVzCut(70.); me->setVrCut(2.); me->setZLocalCut(3.05);
  const UInt_t mb[] = {640001, 640011, 640021, 640031, 640041, 640051};
  for (int it = 0; it < 6; it++) me->addTriggerId(mb[it]);
  if (chain->Init() == kStERR) { cout << "chain Init failed" << endl; return; }
  Long64_t nAvail = mu->chain()->GetEntries();
  Long64_t nRun = (nEvents > 0 && nEvents < nAvail) ? nEvents : nAvail;
  cout << "runTofMEMu: " << nAvail << " events in chain, processing " << nRun << endl;
  for (Long64_t iev = 0; iev < nRun; iev++) {
    if (iev % 2000 == 0) cout << "  event " << iev << endl;
    chain->Clear();
    Int_t ret = chain->Make(iev);
    if (ret == kStEOF || ret == kStFatal) { cout << "Make returned " << ret << " at " << iev << endl; break; }
  }
  chain->Finish();
  delete chain;
}
