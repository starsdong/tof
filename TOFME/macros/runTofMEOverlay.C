// Mixed event with the production BTOF matching + calibration re-run on overlaid events (StTofMEOverlayMaker).
// root4star -b -q 'macros/runTofMEOverlay.C(1000,"filelists/mudst/test1.list","output/mudst/ov_test.root","y2019")'
void loadLibsAgMLc()
{
  gSystem->Load("libStarAgmlUtil");
  gSystem->Load("libStarAgmlLib");
  gSystem->Load("libGeometry");
  gSystem->Load("libStarGeometry");
}
void runTofMEOverlay(Int_t nEvents = 1000, const char* list = "filelists/mudst/test1.list",
                   const char* out = "output/mudst/ov_test.root", const char* geoTag = "y2019", Bool_t startless = kFALSE, Int_t mode = 0)
{
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");
  loadLibsAgMLc();
  gSystem->Load("StBTofUtil");
  gSystem->Load("StVpdCalibMaker");
  gSystem->Load("StBTofCalibMaker");
  gSystem->Load("StBTofMatchMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StTofMEOverlayMaker");
  AgModule::SetStacker(new StarTGeoStacker);
  if (TClass::GetClass("StarGeometry")) StarGeometry::Construct(geoTag);
  else (new Geometry())->ConstructGeometry(geoTag);

  StChain* chain = new StChain();
  StMuDstMaker* mu = new StMuDstMaker(0, 0, "", list, "MuDst", 10);
  St_db_Maker* db = new St_db_Maker("db", "StarDb", "MySQL:StarDb", "$STAR/StarDb");
  StBTofMatchMaker* btofMatch = new StBTofMatchMaker();
  btofMatch->setMuDstIn();
  StVpdCalibMaker* vpdCalib = new StVpdCalibMaker();
  vpdCalib->setMuDstIn();
  if (startless) vpdCalib->setUseVpdStart(kFALSE);
  StBTofCalibMaker* btofCalib = new StBTofCalibMaker();
  btofCalib->setMuDstIn();
  StTofMEOverlayMaker* ov = new StTofMEOverlayMaker("tofMEov", mu, out);
  ov->setProdMakers(btofMatch, vpdCalib, btofCalib);
  ov->setNMix(5);
  ov->setMode(mode);
  const UInt_t mb[] = {640001, 640011, 640021, 640031, 640041, 640051};
  for (int it = 0; it < 6; it++) ov->addTriggerId(mb[it]);
  if (chain->Init() == kStERR) { cout << "chain Init failed" << endl; return; }
  Int_t nAvail = mu->chain()->GetEntries();
  Int_t nRun = (nEvents > 0 && nEvents < nAvail) ? nEvents : nAvail;
  cout << "runTofMEOverlay: processing " << nRun << " of " << nAvail << endl;
  for (Int_t i = 0; i < nRun; i++) {
    if (i % 500 == 0) cout << "  event " << i << endl;
    chain->Clear();
    Int_t r = chain->Make(i);
    if (r == kStEOF || r == kStFatal) break;
  }
  chain->Finish();
}
