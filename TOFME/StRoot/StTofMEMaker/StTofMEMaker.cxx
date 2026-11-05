#include "StTofMEMaker.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofHit.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TRandom3.h"

#include <deque>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>

ClassImp(StTofMEMaker)

namespace {
  enum { kCutAll = 0, kCutTrig, kCutVz, kCutVr, kNCut };
  const char* kCutName[kNCut] = {"all", "trigger", "|V_{z}|", "V_{r}"};
  const Int_t kNVz = 14;                       // 10 cm bins in |Vz|<70
  const Int_t kNMult = 9;
  const Int_t kMultEdge[kNMult + 1] = {0, 10, 25, 50, 80, 120, 170, 230, 300, 100000};
  const Double_t kNoSol = 1e30;

  const Double_t kC = 29.9792458;               // cm/ns
  // m2 histogram binning
  const Int_t kNP = 64;   const Double_t kPMax = 3.2;
  const Int_t kNM2 = 600; const Double_t kM2Lo = -1.0, kM2Hi = 5.0;

  struct MixEvt {
    std::vector<float> map;  // per cell: 0 no hit, -1 hit without time, >0 time of flight (ns)
    std::vector<float> mapS; // v9: untimed hits replaced by sampled times (stays -1 if no sample)
    std::vector<float> mapR; // v9 closure: timed hits replaced by sampled times, untimed stay -1
    std::vector<float> mapU; // v10: untimed hits with times sampled from the MuDst unmatched-hit pdf
    Bool_t sampOK;           // v9: every hit of the event got a sampled time (reservoirs warmed up)
    std::vector<unsigned char> nMatch, nDup;
    Int_t nHits;
    Bool_t hasTime;          // event has a valid TOF start (at least one timed hit)
  };
  std::deque<MixEvt> gBuf[kNVz][kNMult];
  std::vector<float> gHitMap(StTofMEMaker::kNCell, 0.f);
  std::vector<unsigned char> gNMatch(StTofMEMaker::kNCell, 0);   // picoDst tracks matched to the cell
  std::vector<unsigned char> gNDup(StTofMEMaker::kNCell, 0);     // BTofHit entries in the cell
  const Double_t kMPi = 0.13957;
  // v9 time reservoirs: class = Vz bin x (half*32 + module); ring buffer of recent timed-hit tof values
  const Int_t kNKey = 64, kResSize = 3000, kResMin = 200;
  std::vector<float> gRes[kNVz][kNKey];
  Int_t gResPos[kNVz][kNKey];
  std::vector<float> gMapS(StTofMEMaker::kNCell, 0.f), gMapR(StTofMEMaker::kNCell, 0.f);
  std::vector<float> gMapU(StTofMEMaker::kNCell, 0.f);
  TH1D* gPdf[kNVz][kNKey];      // v10 MuDst unmatched-hit time pdfs (0 if not loaded)
  TRandom3 gRnd(4357);
  Int_t cellKey(Int_t cid) { const Int_t mi = cid / 6; return ((mi / 32) < 60 ? 0 : 32) + mi % 32; }
  Float_t sampleTof(Int_t iVz, Int_t cid)
  {
    const std::vector<float>& r = gRes[iVz][cellKey(cid)];
    if ((Int_t)r.size() < kResMin) return -1.f;
    return r[gRnd.Integer(r.size())];
  }
  void addTof(Int_t iVz, Int_t cid, Float_t t)
  {
    const Int_t k = cellKey(cid);
    std::vector<float>& r = gRes[iVz][k];
    if ((Int_t)r.size() < kResSize) r.push_back(t);
    else { r[gResPos[iVz][k]] = t; gResPos[iVz][k] = (gResPos[iVz][k] + 1) % kResSize; }
  }
  // v4 category variables: 0-6, 10, 11 depend on track/event; 7-9 on the chosen hit
  const Int_t kNCatBin[StTofMEMaker::kNK] = {4, 4, 4, 3, 3, 3, 3, 2, 4, 2, 4, 2};
  const char* kCatName[StTofMEMaker::kNK] = {
    "gDCA: <0.5, 0.5-1, 1-2, 2-3 cm", "|Vz-VzVpd|: <3, 3-6, >6 cm, no VPD", "|zLocal|: <1, 1-2, 2-2.6, >2.6 cm",
    "|yLocal|: <0.6, 0.6-1.2, >1.2 cm", "nHitsFit: 15-24, 25-34, >=35", "|eta|: <0.3, 0.3-0.6, >0.6",
    "centrality: 0-10, 10-40, 40-80%", "hit cell: projected, neighbour",
    "hit owner: own (1 track), own (shared), other (1 track), other (shared)", "cell hits: 1, >=2",
    "nSigmaProton: <-3, -3..3, 3..8, >8", "charge: -, +"};
  const Bool_t kHitCat[StTofMEMaker::kNK] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0};

  Double_t dPhi(Double_t a, Double_t b) { return TVector2::Phi_mpi_pi(a - b); }
  // same formula as StBTofUtil/tofPathLength.cc
  Double_t tofPathLength(const TVector3& b, const TVector3& e, Double_t curv)
  {
    const Double_t xd = e.X() - b.X(), yd = e.Y() - b.Y();
    const Double_t C = std::sqrt(xd * xd + yd * yd);
    Double_t sPerp = C;
    if (curv != 0) { const Double_t R = 1. / curv; sPerp = 2 * R * std::asin(C / (2 * R)); }
    const Double_t sZ = std::fabs(e.Z() - b.Z());
    return std::sqrt(sPerp * sPerp + sZ * sZ);
  }
  Double_t mass2(Double_t p, Double_t L, Double_t tof)
  {
    const Double_t beta = L / (kC * tof);
    return p * p * (1. / (beta * beta) - 1.);
  }
}

//_____________________________________________________________________________
StTofMEMaker::StTofMEMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName)
  : StMaker(name), mPicoDstMaker(picoMaker), mOutName(outName), mGeomFile(""), mFile(0),
    mHaveGeom(kFALSE), mNMix(5), mVzCut(70.), mVrCut(2.), mZLocalCut(3.05), mRmean(0), mCurrentRun(-1),
    mRefMultCorr(0), mBadRun(kFALSE)
{
  for (int i = 0; i < kNMod; i++) { mGeoFlag[i] = 0; for (int k = 0; k < kNPar; k++) mGeo[i][k] = 0; }
}

StTofMEMaker::~StTofMEMaker() {}

//_____________________________________________________________________________
Bool_t StTofMEMaker::loadGeom()
{
  if (mGeomFile.Length() == 0) return kFALSE;
  std::ifstream in(mGeomFile.Data());
  if (!in.good()) { LOG_ERROR << "StTofMEMaker: cannot open geometry " << mGeomFile << endm; return kFALSE; }
  std::string line;
  Int_t nOk = 0;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    Int_t idx, tray, mod, flag, n;
    ss >> idx >> tray >> mod >> flag >> n;
    if (idx < 0 || idx >= kNMod) continue;
    for (int k = 0; k < kNPar; k++) ss >> mGeo[idx][k];
    mGeoFlag[idx] = flag;
    if (flag > 0) nOk++;
  }
  // derived quantities
  Double_t rSum = 0; Int_t rN = 0;
  Double_t sx[120] = {0}, sy[120] = {0};
  for (int i = 0; i < kNMod; i++) {
    if (mGeoFlag[i] <= 0) continue;
    TVector3 c(mGeo[i][0], mGeo[i][1], mGeo[i][2]);
    mModPhi[i] = c.Phi(); mModZ[i] = c.Z();
    rSum += c.Perp(); rN++;
    sx[i / 32] += c.X(); sy[i / 32] += c.Y();
  }
  mRmean = rN ? rSum / rN : 215.;
  for (int t = 0; t < 120; t++) mTrayPhi[t] = (sx[t] != 0 || sy[t] != 0) ? atan2(sy[t], sx[t]) : -999.;
  LOG_INFO << "StTofMEMaker: geometry " << mGeomFile << " modules with geometry " << nOk
           << "  <R> = " << mRmean << " cm" << endm;
  return nOk > 0;
}

//_____________________________________________________________________________
Int_t StTofMEMaker::Init()
{
  mHaveGeom = loadGeom();
  for (int i = 0; i < kNVz; i++) for (int k = 0; k < kNKey; k++) gPdf[i][k] = 0;
  if (mPdfFile.Length()) {
    TFile* fp = TFile::Open(mPdfFile.Data());
    TH3* u = fp ? (TH3*)fp->Get("hTU3") : 0;
    if (!u) LOG_ERROR << "StTofMEMaker: no hTU3 in " << mPdfFile << endm;
    else {
      Int_t nFb = 0;
      for (int i = 0; i < kNVz; i++) for (int k = 0; k < kNKey; k++) {
        TH1D* h = u->ProjectionZ(Form("pdf_%d_%d", i, k), i + 1, i + 1, k + 1, k + 1);
        if (h->GetEntries() < 200) { delete h; h = u->ProjectionZ(Form("pdf_%d_%d", i, k), 1, kNVz, k + 1, k + 1); nFb++; }
        h->SetDirectory(0); gPdf[i][k] = h;
      }
      LOG_INFO << "StTofMEMaker: time pdfs loaded from " << mPdfFile << ", " << nFb << " classes use all Vz" << endm;
    }
    if (fp) fp->Close();
  }
  mFile = new TFile(mOutName.Data(), "RECREATE");
  const Double_t pi = TMath::Pi();

  hEvtCut = new TH1D("hEvtCut", "event cut flow;;N_{events}", kNCut, 0, kNCut);
  for (int i = 0; i < kNCut; i++) hEvtCut->GetXaxis()->SetBinLabel(i + 1, kCutName[i]);
  hNMixed = new TH1D("hNMixed", ";N mixed partners;N_{events}", 21, -0.5, 20.5);
  hGeomSum = new TH2D("hGeomSum", "geometry normal-equation sums;module index;slot", kNMod, 0, kNMod, 32, 0, 32);
  hProjStatus = new TH1D("hProjStatus", ";projection status;N_{tracks}", 6, -0.5, 5.5);
  const char* sl[] = {"ok", "no R crossing", "no module", "outside cells", "|z_{loc}|>cut", "no geometry"};
  for (int i = 0; i < 6; i++) hProjStatus->GetXaxis()->SetBinLabel(i + 1, sl[i]);
  hProjModFlag = new TH1D("hProjModFlag", ";geometry flag of projected module;N_{tracks}", 3, -0.5, 2.5);

  hProdCmp = new TH1D("hProdCmp", ";;N_{tracks}", 4, 0, 4);
  const char* pl[] = {"prod matched", "my proj ok", "same module", "same cell"};
  for (int i = 0; i < 4; i++) hProdCmp->GetXaxis()->SetBinLabel(i + 1, pl[i]);
  hProdCmpVsP = new TH2D("hProdCmpVsP", ";p (GeV/c);category", 40, 0, 4, 4, 0, 4);
  hDyVsP   = new TH2D("hDyVsP", ";p (GeV/c);y_{loc}^{proj}-y_{loc}^{prod} (cm)", 40, 0, 4, 200, -2, 2);
  hDzVsP   = new TH2D("hDzVsP", ";p (GeV/c);z_{loc}^{proj}-z_{loc}^{prod} (cm)", 40, 0, 4, 200, -2, 2);
  hDposVsP = new TH2D("hDposVsP", ";p (GeV/c);|#Deltar| proj - prod (cm)", 40, 0, 4, 200, 0, 4);

  const int nPt = 50; const double ptMax = 5.;
  hPrjPt  = new TH1D("hPrjPt",  "projected;p_{T} (GeV/c);N", nPt, 0, ptMax);
  hSEPt   = new TH1D("hSEPt",   "SE cell match;p_{T} (GeV/c);N", nPt, 0, ptMax);
  hSENbPt = new TH1D("hSENbPt", "SE cell or neighbour;p_{T} (GeV/c);N", nPt, 0, ptMax);
  hMETrPt = new TH1D("hMETrPt", "ME trials;p_{T} (GeV/c);N", nPt, 0, ptMax);
  hMEPt   = new TH1D("hMEPt",   "ME cell match;p_{T} (GeV/c);N", nPt, 0, ptMax);
  hMENbPt = new TH1D("hMENbPt", "ME cell or neighbour;p_{T} (GeV/c);N", nPt, 0, ptMax);
  hProdPt = new TH1D("hProdPt", "prod matched (beta>0), projected;p_{T} (GeV/c);N", nPt, 0, ptMax);
  hPrjEta = new TH1D("hPrjEta", ";#eta;N", 40, -1, 1);
  hSEEta  = new TH1D("hSEEta",  ";#eta;N", 40, -1, 1);
  hMETrEta= new TH1D("hMETrEta",";#eta;N", 40, -1, 1);
  hMEEta  = new TH1D("hMEEta",  ";#eta;N", 40, -1, 1);
  hProdEta= new TH1D("hProdEta",";#eta;N", 40, -1, 1);
  hPrjPhi = new TH1D("hPrjPhi", ";#phi;N", 120, -pi, pi);
  hSEPhi  = new TH1D("hSEPhi",  ";#phi;N", 120, -pi, pi);
  hMETrPhi= new TH1D("hMETrPhi",";#phi;N", 120, -pi, pi);
  hMEPhi  = new TH1D("hMEPhi",  ";#phi;N", 120, -pi, pi);
  hPrjRef = new TH1D("hPrjRef", ";refMult;N", 50, 0, 500);
  hSERef  = new TH1D("hSERef",  ";refMult;N", 50, 0, 500);
  hMETrRef= new TH1D("hMETrRef",";refMult;N", 50, 0, 500);
  hMERef  = new TH1D("hMERef",  ";refMult;N", 50, 0, 500);
  hSETrNh = new TH1D("hSETrNh", ";N BTOF hits;N", 60, 0, 1200);
  hSENh   = new TH1D("hSENh",   ";N BTOF hits;N", 60, 0, 1200);
  hMETrNh = new TH1D("hMETrNh", ";N BTOF hits (mixed event);N", 60, 0, 1200);
  hMENh   = new TH1D("hMENh",   ";N BTOF hits (mixed event);N", 60, 0, 1200);
  hPrjYl  = new TH1D("hPrjYl", ";y_{loc}^{proj} (cm);N", 100, -2.5, 2.5);
  hSEYl   = new TH1D("hSEYl",  ";y_{loc}^{proj} (cm);N", 100, -2.5, 2.5);
  hMEYl   = new TH1D("hMEYl",  ";y_{loc}^{proj} (cm);N", 100, -2.5, 2.5);
  hPrjZl  = new TH1D("hPrjZl", ";z_{loc}^{proj} (cm);N", 80, -4, 4);
  hSEZl   = new TH1D("hSEZl",  ";z_{loc}^{proj} (cm);N", 80, -4, 4);
  hMEZl   = new TH1D("hMEZl",  ";z_{loc}^{proj} (cm);N", 80, -4, 4);
  hProdYl = new TH1D("hProdYl", ";y_{loc}^{prod} (cm);N", 200, -5, 5);
  hPrjTrayMod = new TH2D("hPrjTrayMod", "projected;tray;module", 120, 0.5, 120.5, 32, 0.5, 32.5);
  hSEtrayMod  = new TH2D("hSEtrayMod", "SE matched;tray;module", 120, 0.5, 120.5, 32, 0.5, 32.5);

  // ---- v2
  hCentFlag = new TH1D("hCentFlag", "selected events;;N_{events}", 6, 0, 6);
  const char* cf[] = {"selected", "bad run", "pile-up", "cent. > 80%", "no TOF time", "used for m^{2}"};
  for (int i = 0; i < 6; i++) hCentFlag->GetXaxis()->SetBinLabel(i + 1, cf[i]);
  hCent9 = new TH1D("hCent9", ";centrality bin (0 = 70-80%, 8 = 0-5%);N_{events}", 9, -0.5, 8.5);
  hRefMultCorrH = new TH1D("hRefMultCorrH", ";refMult_{corr};N_{events}", 500, 0, 500);
  hTimedFrac = new TH1D("hTimedFrac", ";fraction of BTOF hits with time;N_{events}", 100, 0, 1);
  hHitTof = new TH1D("hHitTof", ";time of flight of timed hits (ns);N_{hits}", 400, -10, 70);
  hNPrjM2 = new TH2D("hNPrjM2", "projected tracks (m^{2} sample);cent9;p (GeV/c)", 9, -0.5, 8.5, kNP, 0, kPMax);
  hNTrM2  = new TH2D("hNTrM2", "ME trials (m^{2} sample);cent9;p (GeV/c)", 9, -0.5, 8.5, kNP, 0, kPMax);
  hSEStat = new TH3F("hSEStat", "SE status (cell#pm1);cent9;p (GeV/c);status", 9, -0.5, 8.5, kNP, 0, kPMax, 3, -0.5, 2.5);
  hMEStat = new TH3F("hMEStat", "ME status (cell#pm1);cent9;p (GeV/c);status", 9, -0.5, 8.5, kNP, 0, kPMax, 3, -0.5, 2.5);
  hM2SE  = new TH3F("hM2SE",  "SE, cell#pm1;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hM2ME  = new TH3F("hM2ME",  "ME, cell#pm1;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hM2SEc = new TH3F("hM2SEc", "SE, cell;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hM2MEc = new TH3F("hM2MEc", "ME, cell;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hM2Prod= new TH3F("hM2Prod","production #beta;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hDLVsP  = new TH2D("hDLVsP", ";p (GeV/c);L_{proj} - #beta_{prod} c t_{prod} (cm)", kNP, 0, kPMax, 200, -5, 5);
  hMEOvStat = new TH3F("hMEOvStat", "overlay ME status;cent9;p (GeV/c);status", 9, -0.5, 8.5, kNP, 0, kPMax, 5, -0.5, 4.5);
  hM2MEo = new TH3F("hM2MEo", "overlay ME, cell#pm1;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  for (int k = 0; k < kNK; k++) {
    const int nb = kNCatBin[k];
    hSEcat[k] = new TH3F(Form("hSEcat%d", k), Form("SE: %s;category;p (GeV/c);m^{2}", kCatName[k]), nb, -0.5, nb - 0.5, 32, 0, 3.2, 300, -1, 5);
    hMEcat[k] = new TH3F(Form("hMEcat%d", k), Form("ME overlay: %s;category;p (GeV/c);m^{2}", kCatName[k]), nb, -0.5, nb - 0.5, 32, 0, 3.2, 300, -1, 5);
    hNPcat[k] = new TH2D(Form("hNPcat%d", k), Form("projected: %s;category;p (GeV/c)", kCatName[k]), nb, -0.5, nb - 0.5, 32, 0, 3.2);
    hNTcat[k] = new TH2D(Form("hNTcat%d", k), Form("trials: %s;category;p (GeV/c)", kCatName[k]), nb, -0.5, nb - 0.5, 32, 0, 3.2);
    hOvStcat[k] = new TH3F(Form("hOvStcat%d", k), Form("overlay status: %s;category;p (GeV/c);status", kCatName[k]), nb, -0.5, nb - 0.5, 32, 0, 3.2, 5, -0.5, 4.5);
  }
  hDtSE    = new TH2F("hDtSE", "SE;p (GeV/c);t - t_{#pi} (ns)", 32, 0, 3.2, 600, -3, 27);
  hDtSEown = new TH2F("hDtSEown", "SE own hit, 1 track;p (GeV/c);t - t_{#pi} (ns)", 32, 0, 3.2, 600, -3, 27);
  hDtME    = new TH2F("hDtME", "ME overlay;p (GeV/c);t - t_{#pi} (ns)", 32, 0, 3.2, 600, -3, 27);
  hCellPrj    = new TH1D("hCellPrj", "projected (0.6<p<3);projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellOcc    = new TH1D("hCellOcc", "projected cell has a hit;projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellSEpr   = new TH1D("hCellSEpr", "SE hit in projected cell (timed);projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellSEprHI = new TH1D("hCellSEprHI", "SE projected cell, HI;projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellSEnb   = new TH1D("hCellSEnb", "SE neighbour cell (timed);projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellSEnbHI = new TH1D("hCellSEnbHI", "SE neighbour cell, HI;projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellSEnbPR = new TH1D("hCellSEnbPR", "SE neighbour cell, PR;projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellTr     = new TH1D("hCellTr", "ME trials;projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellMEnb   = new TH1D("hCellMEnb", "ME overlay neighbour (timed);projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellMEnbU  = new TH1D("hCellMEnbU", "ME overlay partner hit wins, untimed;projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellMEnbHI = new TH1D("hCellMEnbHI", "ME overlay neighbour, HI;projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hCellMEprHI = new TH1D("hCellMEprHI", "ME overlay projected, HI;projected cell id;N", kNCell, -0.5, kNCell - 0.5);
  hNbYproj = new TH3F("hNbYproj", "SE nb;|y_{loc}^{proj}| (cm);p;m^{2}", 7, 0, 1.75, 32, 0, 3.2, 300, -1, 5);
  hPrYproj = new TH3F("hPrYproj", "SE proj;|y_{loc}^{proj}| (cm);p;m^{2}", 7, 0, 1.75, 32, 0, 3.2, 300, -1, 5);
  hNbYprod = new TH3F("hNbYprod", "SE own nb;y_{loc}^{prod} towards proj. cell (cm);p;m^{2}", 16, -2, 2, 32, 0, 3.2, 300, -1, 5);
  hNbState = new TH3F("hNbState", "SE nb;projected cell state;p;m^{2}", 2, -0.5, 1.5, 32, 0, 3.2, 300, -1, 5);
  hMeNbYproj = new TH3F("hMeNbYproj", "ME nb;|y_{loc}^{proj}| (cm);p;m^{2}", 7, 0, 1.75, 32, 0, 3.2, 300, -1, 5);
  hMePrYproj = new TH3F("hMePrYproj", "ME proj;|y_{loc}^{proj}| (cm);p;m^{2}", 7, 0, 1.75, 32, 0, 3.2, 300, -1, 5);
  hM2MEtot = new TH3F("hM2MEtot", "overlay ME, ToT-blind production choice;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hMEtotSt = new TH3F("hMEtotSt", "overlay ME (ToT-blind) partner wins;cent9;p (GeV/c);0 timed, 1 untimed", 9, -0.5, 8.5, kNP, 0, kPMax, 2, -0.5, 1.5);
  hSEnA    = new TH3F("hSEnA", "hits of this event in candidate cells;cent9;p (GeV/c);n_{A}", 9, -0.5, 8.5, kNP, 0, kPMax, 4, -0.5, 3.5);
  hM2MEtot->Sumw2();
  hM2MEearly = new TH3F("hM2MEearly", "overlay ME, earliest time wins;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hMEearlySt = new TH3F("hMEearlySt", "overlay ME earliest: 0 A, 1 B timed, 2 B untimed only, 3 no B;cent9;p (GeV/c);status", 9, -0.5, 8.5, kNP, 0, kPMax, 4, -0.5, 3.5);
  hM2MEsamp = new TH3F("hM2MEsamp", "overlay ME, sampled times for untimed hits, earliest wins;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hMEsampSt = new TH3F("hMEsampSt", "v9: 0 A wins, 1 B wins, 2 only unsampled, 3 no B hit;cent9;p (GeV/c);status", 9, -0.5, 8.5, kNP, 0, kPMax, 4, -0.5, 3.5);
  hM2MEclos = new TH3F("hM2MEclos", "closure: v8 earliest wins with resampled B times;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hTofRes = new TH3F("hTofRes", "timed hits;Vz bin;half*32+module;tof (ns)", kNVz, -0.5, kNVz - 0.5, kNKey, -0.5, kNKey - 0.5, 400, 0, 80);
  hM2MEmu = new TH3F("hM2MEmu", "overlay ME, MuDst pdf for untimed hits, earliest wins;cent9;p (GeV/c);m^{2} (GeV^{2}/c^{4})", 9, -0.5, 8.5, kNP, 0, kPMax, kNM2, kM2Lo, kM2Hi);
  hMEmuSt = new TH3F("hMEmuSt", "v10: 0 A wins, 1 B wins, 3 no B hit;cent9;p (GeV/c);status", 9, -0.5, 8.5, kNP, 0, kPMax, 4, -0.5, 3.5);
  hMuSampT = new TH1D("hMuSampT", "v10 sampled times;t (ns);hits", 1400, -20, 120);
  hNTrSamp = new TH2F("hNTrSamp", "ME trials, both events sampled;cent9;p (GeV/c)", 9, -0.5, 8.5, kNP, 0, kPMax);
  hSampFlag = new TH1D("hSampFlag", "untimed hits: 0 sampled, 1 reservoir too small;flag;hits", 2, -0.5, 1.5);
  for (int i = 0; i < kNVz; i++) for (int k = 0; k < kNKey; k++) { gRes[i][k].clear(); gResPos[i][k] = 0; }
  hDtSEc = new TH3F("hDtSEc", "SE;cent9;p (GeV/c);t - t_{#pi} (ns)", 9, -0.5, 8.5, 32, 0, 3.2, 300, -3, 27);
  hDtSEnb = new TH2F("hDtSEnb", "SE neighbour;p (GeV/c);t - t_{#pi} (ns)", 32, 0, 3.2, 600, -3, 27);
  hDtMEnb = new TH2F("hDtMEnb", "ME neighbour;p (GeV/c);t - t_{#pi} (ns)", 32, 0, 3.2, 600, -3, 27);
  hDm2VsP = new TH2D("hDm2VsP", ";p (GeV/c);m^{2}_{SE} - m^{2}_{prod} (GeV^{2}/c^{4})", kNP, 0, kPMax, 200, -0.2, 0.2);
  mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();

  LOG_INFO << "StTofMEMaker::Init geometry=" << (mHaveGeom ? "yes" : "NO (geometry-sum mode)")
           << " nMix=" << mNMix << " |zLocal|<" << mZLocalCut << endm;
  return kStOK;
}

//_____________________________________________________________________________
Bool_t StTofMEMaker::passEvent(StPicoEvent* ev)
{
  hEvtCut->Fill(kCutAll + 0.5);
  if (!mTrigIds.empty()) {
    Bool_t ok = kFALSE;
    for (size_t i = 0; i < mTrigIds.size(); i++) if (ev->isTrigger(mTrigIds[i])) { ok = kTRUE; break; }
    if (!ok) return kFALSE;
  }
  hEvtCut->Fill(kCutTrig + 0.5);
  TVector3 v = ev->primaryVertex();
  if (std::fabs(v.Z()) >= mVzCut) return kFALSE;
  hEvtCut->Fill(kCutVz + 0.5);
  if (v.Perp() >= mVrCut) return kFALSE;
  hEvtCut->Fill(kCutVr + 0.5);
  return kTRUE;
}

Int_t StTofMEMaker::candOrder(Int_t gcell, Int_t cell, Double_t yl, Double_t pitch, Bool_t nb, Int_t* cand) const
{
  Int_t n = 0;
  cand[n++] = gcell;
  if (nb) {
    const Int_t up = (cell < 6) ? gcell + 1 : -1, dn = (cell > 1) ? gcell - 1 : -1;
    const Bool_t upFirst = std::fabs(yl - pitch) <= std::fabs(yl + pitch);   // nearer neighbour first
    const Int_t first = upFirst ? up : dn, second = upFirst ? dn : up;
    if (first >= 0) cand[n++] = first;
    if (second >= 0) cand[n++] = second;
  }
  return n;
}

Int_t StTofMEMaker::overlay(const std::vector<float>& se, const std::vector<float>& me, const Int_t* cand, Int_t n,
                            Int_t seSel, Float_t& tof, Int_t& cellOut) const
{
  for (Int_t k = 0; k < n; k++) {
    const Int_t c = cand[k];
    cellOut = c;
    const Bool_t sePresent = (seSel >= 0) ? (c == seSel) : (se[c] != 0.f);
    const Float_t mv = me[c];
    if (!sePresent && mv == 0.f) continue;
    if (mv == 0.f) return 0;                                      // this event's hit comes first
    if (!sePresent) { if (mv > 0) { tof = mv; return 1; } return 2; }
    if (seSel >= 0 && mv > 0) { if (mv < se[c]) { tof = mv; return 1; } return 0; }   // same cell: earlier wins
    return 3;
  }
  return 4;
}

Int_t StTofMEMaker::pickHit(const std::vector<float>& m, Int_t gcell, Int_t cell, Double_t yl, Double_t pitch,
                            Bool_t nb, Float_t& tof, Int_t& gcellOut) const
{
  Int_t cand[3];
  const Int_t n = candOrder(gcell, cell, yl, pitch, nb, cand);
  Bool_t any = kFALSE;
  for (Int_t k = 0; k < n; k++) {
    const Float_t v = m[cand[k]];
    if (v > 0) { tof = v; gcellOut = cand[k]; return 1; }
    if (v < 0) any = kTRUE;
  }
  return any ? 2 : 0;
}

Int_t StTofMEMaker::multClass(Int_t r) const
{
  for (int i = 0; i < kNMult; i++) if (r < kMultEdge[i + 1]) return i;
  return kNMult - 1;
}

//_____________________________________________________________________________
Int_t StTofMEMaker::project(const StPicoPhysicalHelix& h, Int_t& cell, Double_t& yl, Double_t& zl,
                            TVector3& pos, Int_t& status)
{
  std::pair<Double_t, Double_t> ss = h.pathLength(mRmean);
  Double_t s = kNoSol;
  if (ss.first > 0 && ss.first < 1e4) s = ss.first;
  if (ss.second > 0 && ss.second < 1e4 && ss.second < s) s = ss.second;
  if (s >= kNoSol) { status = 1; return -1; }
  TVector3 p0 = h.at(s);
  const Double_t phi0 = p0.Phi(), z0 = p0.Z();

  Int_t best = -1; Double_t bestAbsZ = 1e9;
  Bool_t anyPlane = kFALSE, anyInY = kFALSE;
  for (int t = 0; t < 120; t++) {
    if (mTrayPhi[t] < -900 || std::fabs(dPhi(phi0, mTrayPhi[t])) > 0.09) continue;
    for (int m = 0; m < 32; m++) {
      const int i = t * 32 + m;
      if (mGeoFlag[i] <= 0 || std::fabs(z0 - mModZ[i]) > 12.) continue;
      TVector3 C(mGeo[i][0], mGeo[i][1], mGeo[i][2]);
      TVector3 U(mGeo[i][3], mGeo[i][4], mGeo[i][5]);
      TVector3 V(mGeo[i][6], mGeo[i][7], mGeo[i][8]);
      TVector3 W(mGeo[i][9], mGeo[i][10], mGeo[i][11]);
      TVector3 N = U.Cross(V);
      if (N.Mag() <= 0) continue;
      N = N.Unit();
      Double_t sp = h.pathLength(C, N);
      if (!(sp > 0 && sp < 1e4)) continue;
      anyPlane = kTRUE;
      TVector3 x = h.at(sp);
      TVector3 d = x - C;
      // solve d = a U + b V (least squares in the plane)
      const Double_t uu = U.Dot(U), vv = V.Dot(V), uv = U.Dot(V);
      const Double_t du = d.Dot(U), dv = d.Dot(V), det = uu * vv - uv * uv;
      if (det <= 0) continue;
      const Double_t a = (du * vv - dv * uv) / det;   // in yLocal units (U ~ unit vector)
      const Double_t b = (dv * uu - du * uv) / det;   // in zLocal units
      const Double_t pitch = W.Dot(U) / uu;          // signed cell step in yLocal units
      if (pitch == 0) continue;
      const Int_t c = (Int_t)std::floor(a / pitch + 3.0) + 1;
      if (c < 1 || c > 6) continue;
      anyInY = kTRUE;
      if (std::fabs(b) > mZLocalCut) continue;
      if (std::fabs(b) < bestAbsZ) {
        bestAbsZ = std::fabs(b); best = i; cell = c;
        yl = a - (c - 3.5) * pitch; zl = b; pos = x;
      }
    }
  }
  if (best >= 0) { status = 0; return best; }
  status = !anyPlane ? 2 : (!anyInY ? 3 : 4);
  return -1;
}

//_____________________________________________________________________________
Int_t StTofMEMaker::Make()
{
  StPicoDst* dst = mPicoDstMaker->picoDst();
  if (!dst) return kStWarn;
  StPicoEvent* ev = dst->event();
  if (!ev) return kStWarn;
  if (!passEvent(ev)) return kStOK;

  if (ev->runId() != mCurrentRun) {         // never mix across runs (detector status)
    for (int i = 0; i < kNVz; i++) for (int j = 0; j < kNMult; j++) gBuf[i][j].clear();
    mCurrentRun = ev->runId();
    mRefMultCorr->init(mCurrentRun);
    mBadRun = mRefMultCorr->isBadRun(mCurrentRun);
    LOG_INFO << "StTofMEMaker: run " << mCurrentRun << (mBadRun ? " is a StRefMultCorr bad run" : "") << endm;
  }
  const TVector3 pv = ev->primaryVertex();
  const Float_t B = ev->bField();
  const Int_t refMult = ev->refMult();
  const Int_t iVz = std::min(kNVz - 1, std::max(0, (Int_t)((pv.Z() + mVzCut) / (2 * mVzCut) * kNVz)));
  const Int_t iMult = multClass(refMult);
  std::deque<MixEvt>& buf = gBuf[iVz][iMult];
  hNMixed->Fill(buf.size());

  // hit map of this event: -1 = hit, then the time of flight from any matched picoDst track
  std::fill(gHitMap.begin(), gHitMap.end(), 0.f);
  const Int_t nHits = dst->numberOfBTofHits();
  for (Int_t i = 0; i < nHits; i++) {
    StPicoBTofHit* hh = dst->btofHit(i);
    if (hh && hh->id() >= 0 && hh->id() < kNCell) gHitMap[hh->id()] = -1.f;
  }
  std::fill(gNMatch.begin(), gNMatch.end(), 0);
  std::fill(gNDup.begin(), gNDup.end(), 0);
  for (Int_t i = 0; i < nHits; i++) {
    StPicoBTofHit* hh = dst->btofHit(i);
    if (hh && hh->id() >= 0 && hh->id() < kNCell && gNDup[hh->id()] < 255) gNDup[hh->id()]++;
  }
  Int_t nTimed = 0;
  const Int_t nPid = dst->numberOfBTofPidTraits();
  for (Int_t i = 0; i < nPid; i++) {
    StPicoBTofPidTraits* pt = dst->btofPidTraits(i);
    if (!pt || pt->btofMatchFlag() <= 0) continue;
    { const Int_t cm = pt->btofCellId(); if (cm >= 0 && cm < kNCell && gNMatch[cm] < 255) gNMatch[cm]++; }
    if (pt->btofBeta() <= 0 || pt->btof() <= 0) continue;
    const Int_t cid = pt->btofCellId();
    if (cid < 0 || cid >= kNCell || gHitMap[cid] == 0.f) continue;   // only cells present in the hit list
    if (gHitMap[cid] < 0) { nTimed++; hHitTof->Fill(pt->btof()); }
    gHitMap[cid] = pt->btof();
  }
  const Bool_t hasTime = nTimed > 0;
  // v9: sampled maps from the reservoirs of earlier events, then feed this event's timed hits
  Bool_t sampOK = kTRUE;
  for (Int_t cid = 0; cid < kNCell; cid++) {
    const Float_t v = gHitMap[cid];
    gMapS[cid] = v; gMapR[cid] = v;
    if (v < 0) { const Float_t s = sampleTof(iVz, cid); gMapS[cid] = s; hSampFlag->Fill(s > 0 ? 0 : 1); if (s <= 0) sampOK = kFALSE; }
    else if (v > 0) { gMapR[cid] = sampleTof(iVz, cid); if (gMapR[cid] <= 0) sampOK = kFALSE; }
  }
  for (Int_t cid = 0; cid < kNCell; cid++) {      // v10 map
    const Float_t v = gHitMap[cid];
    gMapU[cid] = v;
    if (v < 0) {
      TH1D* pd = gPdf[iVz][cellKey(cid)];
      if (pd) { const Double_t t = pd->GetRandom(); hMuSampT->Fill(t); gMapU[cid] = std::max(0.01, t); }
    }
  }
  for (Int_t cid = 0; cid < kNCell; cid++)
    if (gHitMap[cid] > 0) { addTof(iVz, cid, gHitMap[cid]); hTofRes->Fill(iVz, cellKey(cid), gHitMap[cid]); }
  if (nHits > 0) hTimedFrac->Fill(Double_t(nTimed) / nHits);

  // centrality (StRefMultCorr); m2 histograms only for good, non-pile-up events with TOF time
  mRefMultCorr->initEvent(refMult, pv.Z(), ev->ZDCx());
  const Bool_t pileUp = mRefMultCorr->isPileUpEvent(refMult, ev->nBTOFMatch(), pv.Z());
  const Int_t cent9 = mRefMultCorr->getCentralityBin9();
  hCentFlag->Fill(0.5);
  Bool_t useM2 = kFALSE;
  if (mBadRun) hCentFlag->Fill(1.5);
  else if (pileUp) hCentFlag->Fill(2.5);
  else if (cent9 < 0 || cent9 > 8) hCentFlag->Fill(3.5);
  else if (!hasTime) hCentFlag->Fill(4.5);
  else {
    useM2 = kTRUE; hCentFlag->Fill(5.5);
    hCent9->Fill(cent9); hRefMultCorrH->Fill(mRefMultCorr->getRefMultCorr());
  }

  const Int_t nTrk = dst->numberOfTracks();
  for (Int_t i = 0; i < nTrk; i++) {
    StPicoTrack* t = dst->track(i);
    if (!t || !t->isPrimary()) continue;
    const TVector3 mom = t->pMom();
    const Float_t pt = mom.Perp();
    if (pt < 0.15) continue;
    const Int_t nFit = t->nHitsFit();
    const Float_t ratio = (t->nHitsMax() > 0) ? Float_t(nFit) / t->nHitsMax() : 0.;
    if (nFit < 15 || ratio < 0.52f || t->gDCA(pv).Mag() > 3.) continue;  // 0.52f: same float comparison as StTofMEQAMaker
    const Float_t eta = mom.PseudoRapidity();
    if (std::fabs(eta) > 1.0) continue;
    const Float_t p = mom.Mag(), phi = mom.Phi();

    // production match info
    StPicoBTofPidTraits* tof = t->isTofTrack() ? dst->btofPidTraits(t->bTofPidTraitsIndex()) : 0;
    const Bool_t prodMatch = tof && tof->btofMatchFlag() > 0;

    // geometry sums from production matches
    if (prodMatch) {
      const Int_t cid = tof->btofCellId();
      if (cid >= 0 && cid < kNCell) {
        const Int_t mi = cid / 6;
        const Double_t xr[4] = {1., tof->btofYLocal(), tof->btofZLocal(), (cid % 6 + 1) - 3.5};
        const TVector3 hp = tof->btofHitPos();
        const Double_t yv[3] = {hp.X(), hp.Y(), hp.Z()};
        for (int a = 0; a < 4; a++) {
          for (int b = 0; b < 4; b++) hGeomSum->AddBinContent(hGeomSum->GetBin(mi + 1, a * 4 + b + 1), xr[a] * xr[b]);
          for (int c = 0; c < 3; c++) hGeomSum->AddBinContent(hGeomSum->GetBin(mi + 1, 16 + a * 3 + c + 1), xr[a] * yv[c]);
        }
        for (int c = 0; c < 3; c++) hGeomSum->AddBinContent(hGeomSum->GetBin(mi + 1, 28 + c + 1), yv[c] * yv[c]);
        hGeomSum->AddBinContent(hGeomSum->GetBin(mi + 1, 31 + 1), 1.);
        hProdYl->Fill(tof->btofYLocal());
      }
    }
    if (!mHaveGeom) continue;

    // ---- projection
    StPicoPhysicalHelix helix = t->helix(B);
    Int_t cell = 0, status = 0; Double_t yl = 0, zl = 0; TVector3 pos;
    const Int_t mi = project(helix, cell, yl, zl, pos, status);
    hProjStatus->Fill(status);

    if (prodMatch) {
      const Int_t cid = tof->btofCellId();
      hProdCmp->Fill(0.5); hProdCmpVsP->Fill(p, 0.5);
      if (mi >= 0) {
        hProdCmp->Fill(1.5); hProdCmpVsP->Fill(p, 1.5);
        if (mi == cid / 6) {
          hProdCmp->Fill(2.5); hProdCmpVsP->Fill(p, 2.5);
          if (cell == cid % 6 + 1) {
            hProdCmp->Fill(3.5); hProdCmpVsP->Fill(p, 3.5);
            hDyVsP->Fill(p, yl - tof->btofYLocal());
            hDzVsP->Fill(p, zl - tof->btofZLocal());
            hDposVsP->Fill(p, (pos - tof->btofHitPos()).Mag());
          }
        }
      }
    }
    if (mi < 0) continue;

    hProjModFlag->Fill(mGeoFlag[mi]);
    const Int_t gcell = mi * 6 + (cell - 1);
    hPrjPt->Fill(pt); hPrjEta->Fill(eta); hPrjPhi->Fill(phi); hPrjRef->Fill(refMult);
    hPrjYl->Fill(yl); hPrjZl->Fill(zl);
    hPrjTrayMod->Fill(mi / 32 + 1, mi % 32 + 1);
    if (prodMatch && tof->btofBeta() > 0) { hProdPt->Fill(pt); hProdEta->Fill(eta); }

    // ---- same event
    hSETrNh->Fill(nHits);
    const Bool_t se = gHitMap[gcell] != 0.f;
    const Bool_t seNb = se || (cell > 1 && gHitMap[gcell - 1] != 0.f) || (cell < 6 && gHitMap[gcell + 1] != 0.f);
    if (se) {
      hSEPt->Fill(pt); hSEEta->Fill(eta); hSEPhi->Fill(phi); hSERef->Fill(refMult); hSENh->Fill(nHits);
      hSEYl->Fill(yl); hSEZl->Fill(zl); hSEtrayMod->Fill(mi / 32 + 1, mi % 32 + 1);
    }
    if (seNb) hSENbPt->Fill(pt);

    // ---- m2: same event, production, validation
    Double_t L = 0, pitch = 0;
    Int_t seSel = -1, nCand = 0, cand[3] = {-1, -1, -1};
    Int_t catv[kNK] = {0};
    Int_t nAcand = 0;
    if (useM2) {
      L = tofPathLength(pv, pos, helix.curvature());
      const TVector3 U(mGeo[mi][3], mGeo[mi][4], mGeo[mi][5]), W(mGeo[mi][9], mGeo[mi][10], mGeo[mi][11]);
      pitch = W.Dot(U) / U.Dot(U);
      hNPrjM2->Fill(cent9, p);
      // v4 track/event categories
      { const Double_t dca = t->gDCA(pv).Mag();
        catv[0] = dca < 0.5 ? 0 : (dca < 1 ? 1 : (dca < 2 ? 2 : 3));
        const Double_t vv = ev->vzVpd(), dvz = std::fabs(pv.Z() - vv);
        catv[1] = (std::fabs(vv) > 200 || vv == 0) ? 3 : (dvz < 3 ? 0 : (dvz < 6 ? 1 : 2));
        const Double_t az = std::fabs(zl);
        catv[2] = az < 1 ? 0 : (az < 2 ? 1 : (az < 2.6 ? 2 : 3));
        const Double_t ay = std::fabs(yl);
        catv[3] = ay < 0.6 ? 0 : (ay < 1.2 ? 1 : 2);
        catv[4] = nFit < 25 ? 0 : (nFit < 35 ? 1 : 2);
        catv[5] = std::fabs(eta) < 0.3 ? 0 : (std::fabs(eta) < 0.6 ? 1 : 2);
        catv[6] = cent9 >= 7 ? 0 : (cent9 >= 4 ? 1 : 2);
        const Double_t ns = t->nSigmaProton();
        catv[10] = ns < -3 ? 0 : (ns < 3 ? 1 : (ns < 8 ? 2 : 3));
        catv[11] = t->charge() > 0 ? 1 : 0;
        for (int k = 0; k < kNK; k++) if (!kHitCat[k]) hNPcat[k]->Fill(catv[k], p); }
      if (p > 0.6 && p < 3) { hCellPrj->Fill(gcell); if (gHitMap[gcell] != 0.f) hCellOcc->Fill(gcell); }
      Float_t tofH = 0; Int_t gh = -1;
      const Int_t st = pickHit(gHitMap, gcell, cell, yl, pitch, kTRUE, tofH, gh);
      hSEStat->Fill(cent9, p, st);
      if (st == 1) seSel = gh;
      nCand = candOrder(gcell, cell, yl, pitch, kTRUE, cand);
      nAcand = 0;
      for (Int_t q = 0; q < nCand; q++) if (gHitMap[cand[q]] != 0.f) nAcand++;
      hSEnA->Fill(cent9, p, nAcand);
      if (st == 1) {
        const Double_t m2 = mass2(p, L, tofH);
        hM2SE->Fill(cent9, p, m2);
        { const Bool_t own = prodMatch && tof->btofCellId() == gh;
          const Int_t nm = gNMatch[gh];
          catv[7] = (gh == gcell) ? 0 : 1;
          catv[8] = own ? (nm >= 2 ? 1 : 0) : (nm >= 2 ? 3 : 2);
          catv[9] = gNDup[gh] >= 2 ? 1 : 0;
          for (int k = 0; k < kNK; k++) hSEcat[k]->Fill(catv[k], p, m2);
          const Double_t dt = tofH - L / kC * std::sqrt(1. + kMPi * kMPi / (p * p));
          hDtSE->Fill(p, dt);
          hDtSEc->Fill(cent9, p, dt);
          if (catv[8] == 0) hDtSEown->Fill(p, dt);
          if (catv[7] == 1) {
            hDtSEnb->Fill(p, dt);
            hNbYproj->Fill(std::fabs(yl), p, m2);
            hNbState->Fill(gHitMap[gcell] < 0 ? 1 : 0, p, m2);
            if (catv[8] == 0) {   // own production hit: sign yLocal so that + points towards the projected cell
              const Double_t yp = tof->btofYLocal();
              const Int_t dirToProj = (gcell > gh) ? 1 : -1;      // projected cell is at higher/lower cell index
              hNbYprod->Fill(yp * dirToProj * (pitch > 0 ? 1 : -1), p, m2);
            }
          } else hPrYproj->Fill(std::fabs(yl), p, m2);
          if (p > 0.6 && p < 3) {
            const Bool_t hi = m2 > 1.6 && m2 < 2.8, prr = m2 > 0.6 && m2 < 1.2;
            if (catv[7] == 0) { hCellSEpr->Fill(gcell); if (hi) hCellSEprHI->Fill(gcell); }
            else { hCellSEnb->Fill(gcell); if (hi) hCellSEnbHI->Fill(gcell); if (prr) hCellSEnbPR->Fill(gcell); }
          } }
        if (prodMatch && tof->btofBeta() > 0 && gh == tof->btofCellId()) {
          const Double_t bp = tof->btofBeta();
          hDLVsP->Fill(p, L - bp * kC * tof->btof());
          hDm2VsP->Fill(p, m2 - p * p * (1. / (bp * bp) - 1.));
        }
      }
      Float_t tofC = 0; Int_t gc = -1;
      if (pickHit(gHitMap, gcell, cell, yl, pitch, kFALSE, tofC, gc) == 1) hM2SEc->Fill(cent9, p, mass2(p, L, tofC));
      if (prodMatch && tof->btofBeta() > 0) {
        const Double_t bp = tof->btofBeta();
        hM2Prod->Fill(cent9, p, p * p * (1. / (bp * bp) - 1.));
      }
    }

    // ---- mixed events
    for (size_t k = 0; k < buf.size(); k++) {
      const MixEvt& me = buf[k];
      hMETrPt->Fill(pt); hMETrEta->Fill(eta); hMETrPhi->Fill(phi); hMETrRef->Fill(refMult); hMETrNh->Fill(me.nHits);
      const Bool_t m0 = me.map[gcell] != 0.f;
      const Bool_t mNb = m0 || (cell > 1 && me.map[gcell - 1] != 0.f) || (cell < 6 && me.map[gcell + 1] != 0.f);
      if (m0) { hMEPt->Fill(pt); hMEEta->Fill(eta); hMEPhi->Fill(phi); hMERef->Fill(refMult); hMENh->Fill(me.nHits);
                hMEYl->Fill(yl); hMEZl->Fill(zl); }
      if (mNb) hMENbPt->Fill(pt);
      // m2 with the partner's hit time and this track's path length and momentum
      if (useM2 && me.hasTime) {
        hNTrM2->Fill(cent9, p);
        Float_t tofM = 0; Int_t gm = -1;
        const Int_t st = pickHit(me.map, gcell, cell, yl, pitch, kTRUE, tofM, gm);
        hMEStat->Fill(cent9, p, st);
        if (st == 1) hM2ME->Fill(cent9, p, mass2(p, L, tofM));
        if (pickHit(me.map, gcell, cell, yl, pitch, kFALSE, tofM, gm) == 1) hM2MEc->Fill(cent9, p, mass2(p, L, tofM));
        { // v8: earliest time wins (A's chosen timed hit, if any, vs B's timed hits in the candidate cells)
          Float_t best = (seSel >= 0) ? gHitMap[seSel] : 1e9f;
          Bool_t bWins = kFALSE, bAny = kFALSE, bUnt = kFALSE;
          for (Int_t q = 0; q < nCand; q++) {
            const Float_t mv = me.map[cand[q]];
            if (mv != 0.f) bAny = kTRUE;
            if (mv < 0) bUnt = kTRUE;
            if (mv > 0 && mv < best) { best = mv; bWins = kTRUE; }
          }
          Int_t es = 0;
          if (bWins) { es = 1; hM2MEearly->Fill(cent9, p, mass2(p, L, best)); }
          else if (!bAny) es = 3;
          else if (seSel < 0 && bUnt) es = 2;
          hMEearlySt->Fill(cent9, p, es);
        }
        if (gPdf[0][0]) { // v10: earliest time wins, untimed hits with MuDst-sampled times
          Float_t bestA = 1e9f, bestB = 1e9f; Bool_t bAny = kFALSE;
          for (Int_t q = 0; q < nCand; q++) {
            const Float_t a = gMapU[cand[q]], b = me.mapU[cand[q]];
            if (a > 0 && a < bestA) bestA = a;
            if (b > 0 && b < bestB) bestB = b;
            if (me.map[cand[q]] != 0.f) bAny = kTRUE;
          }
          Int_t us = 0;
          if (bestB < bestA) { us = 1; hM2MEmu->Fill(cent9, p, mass2(p, L, bestB)); }
          else if (!bAny) us = 3;
          hMEmuSt->Fill(cent9, p, us);
        }
        if (sampOK && me.sampOK) { // v9: earliest time wins among all hits of A (real or sampled) and B (real or sampled)
          hNTrSamp->Fill(cent9, p);
          Float_t bestA = 1e9f, bestB = 1e9f;
          Bool_t bAny = kFALSE;
          for (Int_t q = 0; q < nCand; q++) {
            const Float_t a = gMapS[cand[q]], b = me.mapS[cand[q]];
            if (a > 0 && a < bestA) bestA = a;
            if (me.map[cand[q]] != 0.f) bAny = kTRUE;
            if (b > 0 && b < bestB) bestB = b;
          }
          Int_t vs = 0;
          if (bestB < bestA) { vs = 1; hM2MEsamp->Fill(cent9, p, mass2(p, L, bestB)); }
          else if (!bAny) vs = 3;
          else if (bestB > 1e8f && bestA > 1e8f) vs = 2;
          hMEsampSt->Fill(cent9, p, vs);
          // closure: v8 rule (A's chosen timed hit) vs B's timed hits with resampled times
          Float_t best = (seSel >= 0) ? gHitMap[seSel] : 1e9f; Bool_t w = kFALSE;
          for (Int_t q = 0; q < nCand; q++) { const Float_t r = me.mapR[cand[q]]; if (r > 0 && r < best) { best = r; w = kTRUE; } }
          if (w) hM2MEclos->Fill(cent9, p, mass2(p, L, best));
        }
        { // v7: ToT-blind production choice among all hits of A and B in the candidate cells
          Int_t nB = 0;
          for (Int_t q = 0; q < nCand; q++) if (me.map[cand[q]] != 0.f) nB++;
          if (nB > 0) {
            const Double_t wq = 1. / (nAcand + nB);
            for (Int_t q = 0; q < nCand; q++) {
              const Float_t mv = me.map[cand[q]];
              if (mv > 0) { hM2MEtot->Fill(cent9, p, mass2(p, L, mv), wq); hMEtotSt->Fill(cent9, p, 0., wq); }
              else if (mv < 0) hMEtotSt->Fill(cent9, p, 1., wq);
            }
          }
        }
        Float_t tofO = 0; Int_t co = -1;
        const Int_t so = overlay(gHitMap, me.map, cand, nCand, seSel, tofO, co);
        hMEOvStat->Fill(cent9, p, so);
        for (int k = 0; k < kNK; k++) if (!kHitCat[k]) { hNTcat[k]->Fill(catv[k], p); hOvStcat[k]->Fill(catv[k], p, so); }
        if (p > 0.6 && p < 3) { hCellTr->Fill(gcell); if (so == 2) hCellMEnbU->Fill(gcell); }
        if (so == 1) {
          const Double_t m2o = mass2(p, L, tofO);
          hM2MEo->Fill(cent9, p, m2o);
          catv[7] = (co == gcell) ? 0 : 1;
          catv[8] = me.nMatch[co] >= 2 ? 3 : 2;
          catv[9] = me.nDup[co] >= 2 ? 1 : 0;
          for (int k = 0; k < kNK; k++) hMEcat[k]->Fill(catv[k], p, m2o);
          const Double_t dto = tofO - L / kC * std::sqrt(1. + kMPi * kMPi / (p * p));
          hDtME->Fill(p, dto);
          if (catv[7] == 1) { hDtMEnb->Fill(p, dto); hMeNbYproj->Fill(std::fabs(yl), p, m2o); }
          else hMePrYproj->Fill(std::fabs(yl), p, m2o);
          if (p > 0.6 && p < 3) {
            const Bool_t hi = m2o > 1.6 && m2o < 2.8;
            if (catv[7] == 1) { hCellMEnb->Fill(gcell); if (hi) hCellMEnbHI->Fill(gcell); }
            else if (hi) hCellMEprHI->Fill(gcell);
          }
        }
      }
    }
  }

  // push this event into the mixing buffer
  MixEvt cur; cur.map = gHitMap; cur.mapS = gMapS; cur.mapR = gMapR; cur.mapU = gMapU; cur.sampOK = sampOK; cur.nHits = nHits; cur.hasTime = hasTime; cur.nMatch = gNMatch; cur.nDup = gNDup;
  buf.push_back(cur);
  while ((Int_t)buf.size() > mNMix) buf.pop_front();
  return kStOK;
}

//_____________________________________________________________________________
Int_t StTofMEMaker::Finish()
{
  LOG_INFO << "StTofMEMaker::Finish events: all " << hEvtCut->GetBinContent(1) << " selected "
           << hEvtCut->GetBinContent(kNCut) << "  projected tracks " << hPrjPt->GetEntries()
           << "  ME trials " << hMETrPt->GetEntries() << "  m2 events " << hCentFlag->GetBinContent(6)
           << "  m2 SE " << hM2SE->GetEntries() << "  m2 ME " << hM2ME->GetEntries() << endm;
  mFile->cd();
  mFile->Write();
  mFile->Close();
  return kStOK;
}
