
// #include <home/prottay/alice/AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C>
#include "YieldMean.C"
#include "../src/style.h"
#include "../src/common.h"
using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par)
{

  Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
  return (p);
}

TLegend *DrawLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{

  TLegend *legend = new TLegend(x1, y1, x2, y2);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  return legend;
}

void meanpT_dNdy_Levy()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat(0);

  int setnum = 0;

  // TFile *f1 = new TFile("spectra.root", "READ");
  // TFile *f2 = new TFile("spectra.root", "READ");
  string pathtofile = "/home/sawan/k892_postprocessing/output/LHC24_ab_INEL/common/spectra.root"; // only TPC
  // string pathtofile = "/home/sawan/k892_postprocessing/output/LHC24_ab_INEL/common/spectra_1.root"; // TPC+TOF
  TFile *f1 = new TFile(pathtofile.c_str(), "READ");
  TFile *f2 = new TFile(pathtofile.c_str(), "READ");

  if (f1->IsZombie() || f2->IsZombie())
  {
    cout << "Error: files not found" << endl;
    return;
  }

  TH1D *h1 = (TH1D *)f1->Get("lf-k892analysis/K892/0/hCorrectedYields");
  TH1D *h2 = (TH1D *)f2->Get("lf-k892analysis/K892/0/hCorrectedYields");

  cout << "stat. error: " << h2->GetBinError(4) << endl;
  for (int i = 1; i <= h2->GetNbinsX(); i++)
  {
    double systemerr = (0.09 * h2->GetBinContent(i));
    h2->SetBinError(i, systemerr);
  }
  cout << "Sytematic error: " << h2->GetBinError(4) << endl;

  TCanvas *c = new TCanvas("", "", 720, 720);
  c->SetLogy();
  SetCanvasStyle(c, 0.14, 0.01, 0.01, 0.12);
  c->cd();

  h1->Draw();

  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetXaxis()->SetTitleSize(0.06);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetNdivisions(511);
  h1->GetYaxis()->SetTitleOffset(1.35);
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetTitleSize(0.06);
  h1->GetYaxis()->SetLabelSize(0.06);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetNdivisions(509);
  h1->GetXaxis()->SetRangeUser(0.0, 30.0);

  h1->SetMarkerColor(4);
  h1->SetMarkerStyle(20);
  SetHistoQA(h1);
  h1->SetMarkerSize(1.3);
  h1->GetYaxis()->SetTitle("(1/N)d^{2}N/(dydp_{T}) (GeV/c)^{-1}");
  h1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitleOffset(1.35);

  // check_SDRY
  TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 10.0, 4);
  fitFcn->SetParameter(0, 5.0);
  fitFcn->SetParameter(1, 0.0007);
  fitFcn->FixParameter(2, 1.285);
  fitFcn->SetParameter(3, 0.3);
  fitFcn->SetParNames("n", "dn/dy", "mass", "T");

  /*************meanpT*****************byresonance*******************package*************************/
  Double_t min = 0;
  Double_t max = 5;
  Double_t loprecision = 0.01;
  Double_t hiprecision = 0.1;
  Option_t *opt = "RI+";
  TString logfilename = "log.root";
  Double_t minfit = 0;
  Double_t maxfit = 5;
  // Double_t maxfit=8.0;

  TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
  hout->Draw();

  // hout->GetXaxis()->SetRangeUser(0,20);

  for (Int_t ip = 0; ip < 9; ip++)
  {
    cout << hout->GetBinContent(ip + 1) << endl;
  }

  gPad->Update();
  c->cd();
  hout->SetLineWidth(2);
  TLegend *lp2 = new TLegend(0.68, 0.48, 0.88, 0.6);
  lp2->SetTextFont(42);
  lp2->SetTextSize(0.045);
  lp2->SetFillColor(0);
  lp2->SetLineColor(1);
  lp2->SetBorderSize(0);
  lp2->AddEntry(h1, "K^{*0} + #bar{K^{*0}}", "lep");
  hout->SetLineColor(kRed);
  lp2->AddEntry(hout, "Levy Fit", "l");
  lp2->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.045);
  lat.SetTextFont(42);
  lat.DrawLatex(0.65, 0.62, "pp #sqrt{s} = 13.6 TeV");
  // c->SaveAs(("check_all/levyfit" + to_string(setnum) + ".png").c_str());
  // c->SaveAs("levyfit_latest.png");
  // TFile *fyield = new TFile("yieldfit.root", "RECREATE");
  // c->Write("yieldfit");
}
