// Closure test: re-run BTOF matching + calibration (MuDst mode) and compare with production.
// root4star -b -q 'macros/runTofClosure.C(1000,"filelists/mudst/test1.list","output/mudst/closure.root","y2019")'
void loadLibsAgMLc()
{
  gSystem->Load("libStarAgmlUtil");
  gSystem->Load("libStarAgmlLib");
  gSystem->Load("libGeometry");
  gSystem->Load("libStarGeometry");
}
void runTofClosure(Int_t nEvents = 1000, const char* list = "filelists/mudst/test1.list",
                   const char* out = "output/mudst/closure.root", const char* geoTag = "y2019", Bool_t startless = kFALSE)
{
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");
  loadLibsAgMLc();
  gSystem->Load("StBTofUtil");
  gSystem->Load("StVpdCalibMaker");
  gSystem->Load("StBTofCalibMaker");
  gSystem->Load("StBTofMatchMaker");
  gSystem->Load("StTofClosureMaker");
  AgModule::SetStacker(new StarTGeoStacker);
  if (TClass::GetClass("StarGeometry")) StarGeometry::Construct(geoTag);
  else (new Geometry())->ConstructGeometry(geoTag);

  StChain* chain = new StChain();
  StMuDstMaker* mu = new StMuDstMaker(0, 0, "", list, "MuDst", 10);
  St_db_Maker* db = new St_db_Maker("db", "StarDb", "MySQL:StarDb", "$STAR/StarDb");
  StTofClosureMaker* snap = new StTofClosureMaker("snap", mu, 0);
  StBTofMatchMaker* btofMatch = new StBTofMatchMaker();
  btofMatch->setMuDstIn();
  StVpdCalibMaker* vpdCalib = new StVpdCalibMaker();
  vpdCalib->setMuDstIn();
  if (startless) vpdCalib->setUseVpdStart(kFALSE);
  StBTofCalibMaker* btofCalib = new StBTofCalibMaker();
  btofCalib->setMuDstIn();
  StTofClosureMaker* cmp = new StTofClosureMaker("cmp", mu, 1, out);
  if (chain->Init() == kStERR) { cout << "chain Init failed" << endl; return; }
  Int_t nAvail = mu->chain()->GetEntries();
  Int_t nRun = (nEvents > 0 && nEvents < nAvail) ? nEvents : nAvail;
  cout << "runTofClosure: processing " << nRun << " of " << nAvail << endl;
  for (Int_t i = 0; i < nRun; i++) {
    chain->Clear();
    Int_t r = chain->Make(i);
    if (r == kStEOF || r == kStFatal) break;
  }
  chain->Finish();
}
