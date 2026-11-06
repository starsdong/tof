// starsim_hijing.C : HIJING 1.383 Au+Au at sqrt(s_NN) = 19.6 GeV through the full STAR GEANT3 geometry (AgML).
// Output: <out>.fzd (GEANT hits for bfc) and <out>.gener.root (generator record).
// usage: root4star -b -q 'starsim_hijing.C(nevents, seed, "out", "y2019a", bmin, bmax, sigmaZ)'
// Based on $STAR/StRoot/StarGenerator/macros/starsim.hijing.C (SL24y).
class St_geant_Maker;
St_geant_Maker* geant_maker = 0;
class StarPrimaryMaker;
StarPrimaryMaker* _primary = 0;
class StarHijing;
StarHijing* hijing = 0;

void command(TString cmd)
{
  if (!geant_maker) geant_maker = (St_geant_Maker*)chain->GetMaker("geant");
  geant_maker->Do(cmd);
}

void starsim_hijing(int nevents = 5, int seed = 1234, const char* out = "hijing_test", const char* tag = "y2019a",
                    double bmin = 0.0, double bmax = 14.5, double sigmaZ = 50.0)
{
  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = Form("%s geant gstar usexgeom agml ", tag);
    bfc(0, simple);
  }
  gSystem->Load("libVMC.so");
  gSystem->Load("StarGeneratorUtil.so");
  gSystem->Load("StarGeneratorEvent.so");
  gSystem->Load("StarGeneratorBase.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libHijing1_383.so");
  gSystem->Load("xgeometry.so");

  StarRandom::seed(seed);
  StarRandom::capture();

  _primary = new StarPrimaryMaker();
  _primary->SetFileName(Form("%s.gener.root", out));
  chain->AddBefore("geant", _primary);

  hijing = new StarHijing();
  hijing->SetTitle("Hijing 1.383");
  hijing->SetFrame("CMS", 19.6);
  hijing->SetBlue("Au");
  hijing->SetYell("Au");
  hijing->SetImpact(bmin, bmax);
  _primary->AddGenerator(hijing);
  // all generated particles go to GEANT (no eta cut: VPD at |eta| 4.2-5.1, and forward particles can make secondaries)
  _primary->SetCuts(1.0E-6, -1., -10., +10.);
  // vertex: Run19 19.6 GeV beam spot; the analysis cuts |Vz| < 70 cm and Vr < 2 cm
  _primary->SetVertex(0., 0., 0.);
  _primary->SetSigma(0.1, 0.1, sigmaZ);

  _primary->Init();
  HiParnt_t& hiparnt = hijing->hiparnt();
  {
    hiparnt.ihpr2(4) = 0;    // jet quenching off
    hiparnt.ihpr2(3) = 0;    // no triggered hard scattering
    hiparnt.ihpr2(8) = 10;   // max jets per nucleon
    hiparnt.ihpr2(11) = 1;   // baryon production
    hiparnt.ihpr2(12) = 1;   // HIJING does not decay pi0, K0s, Lambda, Sigma, ... (STAR default): GEANT decays them in the detector
  }

  command("gkine -4 0");
  command(Form("gfile o %s.fzd", out));
  for (int i = 0; i < nevents; i++) {
    chain->Clear();
    chain->Make();
  }
  command("call agexit");
}
