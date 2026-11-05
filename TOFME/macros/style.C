#include "TStyle.h"

void style()
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);
  gStyle->SetTickLength(0.04,"X");
  gStyle->SetTickLength(0.04,"Y");
  gStyle->SetGridWidth(1);
  gStyle->SetGridStyle(2);
  gStyle->SetGridColor(kGray+1);

  gStyle->SetFillColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(2);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.18);

  gStyle->SetNdivisions(108,"X");
  gStyle->SetNdivisions(108,"Y");
  gStyle->SetLabelOffset(0.005,"X");
  gStyle->SetLabelOffset(0.01,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetTitleSize(0.07,"X");
  gStyle->SetTitleSize(0.07,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");

  gStyle->SetPalette(1);

  gStyle->SetMarkerSize(1.8);
  gStyle->SetMarkerStyle(20);

  gStyle->SetLegendFillColor(10);
}
