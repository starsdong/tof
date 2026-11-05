// makeGeom.C : solve the per-module BTOF geometry from hGeomSum (StTofMEMaker output)
//   hitPos = C + yLocal*U + zLocal*V + (cell-3.5)*W   (least squares, per module, per coordinate)
// Modules with < nMin entries are filled by rotating the same module of the nearest fitted tray
// of the same half (east/west) about the z axis (flag 2). A closure test predicts every fitted
// module from its neighbour tray the same way.
// usage: root -b -q -l 'macros/makeGeom.C+("output/me_geom_test.root","geom/btof_geom_test.txt")'
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TVector2.h"
#include <cstdio>
#include <cmath>

static const int NM = 3840;
static double G[NM][12];
static int flag[NM], nEnt[NM], src[NM];
static double rms[NM][3];
static double trayPhi[120]; static bool trayOk[120]; static double step[2];

static void rotCopy(int srcI, double alpha, double* out)
{
  for (int v = 0; v < 4; v++) {
    TVector3 x(G[srcI][v * 3], G[srcI][v * 3 + 1], G[srcI][v * 3 + 2]);
    x.RotateZ(alpha);
    out[v * 3] = x.X(); out[v * 3 + 1] = x.Y(); out[v * 3 + 2] = x.Z();
  }
}
static double alphaFor(int t, int ts)
{
  if (trayOk[t] && trayOk[ts]) return TVector2::Phi_mpi_pi(trayPhi[t] - trayPhi[ts]);
  return (t - ts) * step[t / 60];
}
// nearest tray (same half, != t) whose module m is fitted
static int findSrc(int t, int m)
{
  const int h = t / 60;
  for (int d = 1; d < 60; d++)
    for (int sgn = -1; sgn <= 1; sgn += 2) {
      int ts = t + sgn * d;
      if (ts < h * 60 || ts >= h * 60 + 60) continue;
      if (flag[ts * 32 + m] == 1) return ts;
    }
  return -1;
}

void makeGeom(const char* inName = "output/me_geom_test.root", const char* outName = "geom/btof_geom_test.txt",
              int nMin = 20)
{
  TFile* f = TFile::Open(inName);
  TH2D* hs = (TH2D*)f->Get("hGeomSum");
  if (!hs) { printf("no hGeomSum\n"); return; }
  int nFit = 0;
  for (int i = 0; i < NM; i++) {
    flag[i] = 0; src[i] = -1; nEnt[i] = (int)hs->GetBinContent(i + 1, 32);
    for (int k = 0; k < 12; k++) G[i][k] = 0;
    if (nEnt[i] < nMin) continue;
    TMatrixDSym A(4);
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) A(a, b) = hs->GetBinContent(i + 1, a * 4 + b + 1);
    double det = 0;
    TMatrixDSym Ai(A); Ai.Invert(&det);
    if (det == 0) continue;
    bool ok = true;
    for (int c = 0; c < 3; c++) {
      TVectorD r(4);
      for (int a = 0; a < 4; a++) r(a) = hs->GetBinContent(i + 1, 16 + a * 3 + c + 1);
      TVectorD beta = Ai * r;
      for (int a = 0; a < 4; a++) G[i][a * 3 + c] = beta(a);
      // residual sum of squares: yty - 2 beta.r + beta.A.beta
      double yty = hs->GetBinContent(i + 1, 28 + c + 1);
      double rss = yty - 2 * (beta * r) + beta * (A * beta);
      rms[i][c] = sqrt(std::max(0., rss / nEnt[i]));
      if (!std::isfinite(beta(0))) ok = false;
    }
    if (ok) { flag[i] = 1; nFit++; }
  }
  // tray phi from fitted module centres
  for (int t = 0; t < 120; t++) {
    double sx = 0, sy = 0; int n = 0;
    for (int m = 0; m < 32; m++) if (flag[t * 32 + m] == 1) { sx += G[t * 32 + m][0]; sy += G[t * 32 + m][1]; n++; }
    trayOk[t] = n > 0; trayPhi[t] = n ? atan2(sy, sx) : 0;
  }
  // phi step per tray number in each half (median of consecutive differences)
  for (int h = 0; h < 2; h++) {
    double sum = 0; int n = 0;
    for (int t = h * 60; t < h * 60 + 59; t++)
      if (trayOk[t] && trayOk[t + 1]) { sum += TVector2::Phi_mpi_pi(trayPhi[t + 1] - trayPhi[t]); n++; }
    step[h] = n ? sum / n : 0;
  }
  printf("fitted modules %d / %d ; phi step per tray: west %.5f east %.5f rad (6 deg = %.5f)\n",
         nFit, NM, step[0], step[1], TMath::Pi() / 30);

  // closure test on fitted modules
  TFile* fo = new TFile(TString(outName).ReplaceAll(".txt", "_qa.root"), "RECREATE");
  TH1D* hClC = new TH1D("hClosureC", "closure: |#DeltaC| rotated vs fitted;|#DeltaC| (cm);N_{modules}", 200, 0, 2);
  TH1D* hClDir = new TH1D("hClosureDir", "closure: angle(U_{rot},U_{fit});angle (mrad);N_{modules}", 100, 0, 20);
  TH1D* hRms = new TH1D("hFitRms", "fit residual RMS (max of x,y,z);RMS (cm);N_{modules}", 100, 0, 0.1);
  TH1D* hN = new TH1D("hFitN", "entries per fitted module;N;N_{modules}", 100, 0, 1000);
  TH1D* hPitch = new TH1D("hPitch", "cell pitch W.U/|U|^{2};pitch (cm);N_{modules}", 200, -4, 4);
  TH1D* hUmag = new TH1D("hUmag", "|U|;|U|;N_{modules}", 100, 0.9, 1.1);
  TH1D* hR = new TH1D("hR", "module centre radius;R (cm);N_{modules}", 200, 200, 230);
  for (int i = 0; i < NM; i++) {
    if (flag[i] != 1) continue;
    hRms->Fill(std::max(rms[i][0], std::max(rms[i][1], rms[i][2])));
    hN->Fill(nEnt[i]);
    TVector3 U(G[i][3], G[i][4], G[i][5]), W(G[i][9], G[i][10], G[i][11]), C(G[i][0], G[i][1], G[i][2]);
    hPitch->Fill(W.Dot(U) / U.Mag2()); hUmag->Fill(U.Mag()); hR->Fill(C.Perp());
    const int t = i / 32, m = i % 32;
    int best = findSrc(t, m);
    if (best < 0) continue;
    double out[12];
    rotCopy(best * 32 + m, alphaFor(t, best), out);
    TVector3 Cp(out[0], out[1], out[2]), Up(out[3], out[4], out[5]);
    hClC->Fill((Cp - C).Mag());
    hClDir->Fill(1000. * Up.Angle(U));
  }
  // fill missing modules
  int nRot = 0, nMiss = 0;
  for (int i = 0; i < NM; i++) {
    if (flag[i] == 1) continue;
    const int t = i / 32, m = i % 32;
    int best = findSrc(t, m);
    if (best < 0) { nMiss++; continue; }
    rotCopy(best * 32 + m, alphaFor(t, best), G[i]);
    flag[i] = 2; src[i] = best + 1; nRot++;
  }
  FILE* fp = fopen(outName, "w");
  fprintf(fp, "# BTOF module geometry from %s (StTofMEMaker hGeomSum, makeGeom.C)\n", inName);
  fprintf(fp, "# hitPos = C + yLocal*U + zLocal*V + (cell-3.5)*W ; flag 1=fitted 2=rotated copy of tray 'src'\n");
  fprintf(fp, "# idx tray module flag n Cx Cy Cz Ux Uy Uz Vx Vy Vz Wx Wy Wz src rmsx rmsy rmsz\n");
  for (int i = 0; i < NM; i++) {
    fprintf(fp, "%d %d %d %d %d", i, i / 32 + 1, i % 32 + 1, flag[i], nEnt[i]);
    for (int k = 0; k < 12; k++) fprintf(fp, " %.5f", G[i][k]);
    fprintf(fp, " %d %.4f %.4f %.4f\n", src[i], flag[i] == 1 ? rms[i][0] : -1., flag[i] == 1 ? rms[i][1] : -1., flag[i] == 1 ? rms[i][2] : -1.);
  }
  fclose(fp);
  printf("modules: fitted %d, rotated %d, missing %d -> %s\n", nFit, nRot, nMiss, outName);
  printf("fit RMS mean %.4f cm ; closure |dC| mean %.3f RMS %.3f cm, max bin edge %.3f; dir mean %.2f mrad\n",
         hRms->GetMean(), hClC->GetMean(), hClC->GetRMS(), hClC->GetXaxis()->GetBinUpEdge(hClC->FindLastBinAbove(0)),
         hClDir->GetMean());
  printf("pitch mean %.4f RMS %.4f ; |U| mean %.5f ; R mean %.3f RMS %.3f\n", hPitch->GetMean(), hPitch->GetRMS(),
         hUmag->GetMean(), hR->GetMean(), hR->GetRMS());
  fo->Write(); fo->Close();
}
