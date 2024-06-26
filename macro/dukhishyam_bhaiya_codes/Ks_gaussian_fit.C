#include <iostream>
#include "../src/common_glue.h"
using namespace std;
Int_t fReBinSpectra = 1;
// forward declarations
Double_t BreitWignerExpol(Double_t *x, Double_t *par);
Double_t Expol(Double_t *x, Double_t *par);
Double_t VoigtPoly2(Double_t *x, Double_t *par);
Double_t Voigt(Double_t *x, Double_t *par);
Double_t Poly2(Double_t *x, Double_t *par);
Double_t BWPoly2(Double_t *x, Double_t *par);
Double_t BW(Double_t *x, Double_t *par);
Double_t Guass(Double_t *x, Double_t *par);
TCanvas *DrawCanvas(TString opt = "c");
TPad *MyPad(Double_t x1 = 0.01, Double_t y1 = 0.01, Double_t x2 = 0.49, Double_t y2 = 0.49, Float_t lm = 0.2, Float_t rm = 0.05, Float_t tm = 0.05, Float_t bm = 0.2, TString name = "pad");
TGraphErrors *PlotGraph(Int_t NdataPoint, Int_t MarkerColor, Int_t MarkerStyle, Double_t *X, Double_t *ErX, Double_t *Y, Double_t *ErY);
TLegend *DrawLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2);
TLatex *DrawText(Double_t x = 0, Double_t y = 0, Int_t tColor = 2, TString name = "");
void MyStyle();
// void MakeInvQMPlot1_K0s_preliminiary();
void DrawFrame(TH1 *h, Int_t MCol, Int_t MSty);

void Ks_gaussian_fit()
{
  // MyStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  // TFile *fInputFile = new TFile("../AnalysisResultsKshortKshortpp13TeV_fulldataset_aod.root","Read");
  // TKey* key = (TKey*)fInputFile->GetListOfKeys()->At(1);
  // fInputList = (TList*)fInputFile->Get(key->GetName());
  // TH1F *hfsigbkg =(TH1F*)fInputList->FindObject("setK0s.K0s_Mass_k0s_mass");

  const string outputQAfolder_str = kSignalOutput + "/" + kchannel + "/" + kfoldername + "/QA";
  TFile *fInputFile = new TFile(("../" + outputQAfolder_str + "/KsInvMass.root").c_str(), "read");
  TH1F *hfsigbkg = (TH1F *)fInputFile->Get("ksmass");
  if (hfsigbkg == nullptr)
  {
    cout << "Histogram not found" << endl;
    return;
  }

  Int_t STCol = 1;
  Int_t STSty = 20;

  Int_t CBCol = 1;
  Int_t CBSty = 24;

  Float_t lowInv = 1.1;
  Float_t highInv = 2.0;

  Float_t maxst = 820E3;
  Float_t minst = .05E6;

  Float_t maxs = 2.9E+6;
  Float_t mins = 0.4E+6;

  // Double_t fitlow = 0.66;
  // Double_t fithi = 1.1;

  Double_t fitlow = 0.4915;
  Double_t fithi = 0.5043;
  Int_t pTBin1 = 12;

  TString ptbn = "0 #leq #it{p}_{T} < 20 GeV/#it{c}";
  TLegend *l2 = DrawLegend(0.25, 0.05, 0.4, 0.15);
  // TGaxis::SetMaxDigits(3);
  TGaxis::SetMaxDigits(2);

  TF1 *fSigBg1 = new TF1("fSigBg1", Guass, fitlow, fithi, 3);
  fSigBg1->SetParameter(0, 100000);
  fSigBg1->SetParameter(1, 0.47);
  fSigBg1->SetParameter(2, 0.005);

  fSigBg1->SetLineColor(1);
  fSigBg1->SetLineStyle(1);
  fSigBg1->SetLineWidth(2);
  fSigBg1->SetParNames("Norm", "Mass", "width");
  // fSigBg1->Draw();

  // return 0;
  TH1D *hSigFit1 = (TH1D *)hfsigbkg->Clone("hSigFit1");
  TH1D *hSigFit1_clone = (TH1D *)hfsigbkg->Clone("hSigFit1_clone");
  TH1D *hSigFit1_clone_dummy = (TH1D *)hfsigbkg->Clone("hSigFit1_clone_dummy");

  // TCanvas *c2 = new TCanvas("c2","c2",600,600);
  //  c2->cd();
  //  hSigFit1_clone->SetMarkerColor(kRed);
  // hSigFit1_clone->SetFillColor(kRed);

  /*
  Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
    hSigFit1_clone->SetLineColor(ci);
    hSigFit1_clone->SetMarkerStyle(20);
    hSigFit1_clone->GetXaxis()->SetRange(68,129);
    hSigFit1_clone->GetXaxis()->SetLabelFont(42);
    hSigFit1_clone->GetXaxis()->SetLabelSize(0.035);
    hSigFit1_clone->GetXaxis()->SetTitleSize(0.035);
    hSigFit1_clone->GetXaxis()->SetTitleFont(42);
    hSigFit1_clone->GetYaxis()->SetLabelFont(42);
    hSigFit1_clone->GetYaxis()->SetLabelSize(0.035);
    hSigFit1_clone->GetYaxis()->SetTitleSize(0.035);
    hSigFit1_clone->GetYaxis()->SetTitleFont(42);
    hSigFit1_clone->GetZaxis()->SetLabelFont(42);
    hSigFit1_clone->GetZaxis()->SetLabelSize(0.035);
    hSigFit1_clone->GetZaxis()->SetTitleSize(0.035);
    hSigFit1_clone->GetZaxis()->SetTitleFont(42);
  */
  //  hSigFit1_clone->Draw("E3HIST");

  // TH1F *h2 = new TH1F("h2","test hstack",100,-4,4);
  //  h2->FillRandom("gaus",15000);
  //  h2->SetFillColor(kBlue);
  //  h2->Draw();

  // return 0;

  cout << "checkpoint 1" << endl;

  hSigFit1->Fit(fSigBg1, "REI+");

  hSigFit1_clone_dummy->Fit(fSigBg1, "REI+");

  cout << "checkpoint 2" << endl;

  // hSigFit1->Fit(f2);
  /*
   Int_t fitres1;
   Int_t trials1 = 0;
   trials1 = 0;
   do{
    // fitres1 = hSigFit1->Fit(fSigBg1, "REBMS");
    Double_t *Partemp   = fSigBg1->GetParameters();
    fSigBg1->SetParameters(&Partemp[0]);
    Printf("Trial: %d %d", trials1++,fitres1);
    if(trials1 > 1) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres1 != 0);

  Double_t *Par1   = fSigBg1->GetParameters();
  Double_t *Erpar1 = fSigBg1->GetParErrors();

  Double_t yield, mass, width, yielderr, masserr,widtherr;
  yield = fSigBg1->GetParameter(0);
  yielderr = fSigBg1->GetParError(0);

  mass = fSigBg1->GetParameter(1);
  masserr = fSigBg1->GetParError(1);

  width = fSigBg1->GetParameter(2);
  widtherr = fSigBg1->GetParError(2);

  cout << mass << "\t-----width" << width << endl;
   */

  // fSigOnly1->SetParameters(&Par1[0]);
  //  fSigOnly1->SetParErrors(&Erpar1[0]);
  //  fBgOnly1->SetParameters(&Par1[3]);
  //  fBgOnly1->SetParErrors(&Erpar1[3]);

  TCanvas *c11 = DrawCanvas("c11");
  // c11->Range(0,0,1,1);
  c11->SetLeftMargin(0.15);
  c11->SetRightMargin(0.1);
  c11->SetTopMargin(0.08);
  c11->SetBottomMargin(0.1);
  c11->cd();
  // c11->SetLeftMargin(0.5);

  TPad *p3 = MyPad(0.01, 0.01, 0.99, 0.99, 0.14, 0.02, 0.06, 0.12, "pad3");
  p3->SetLeftMargin(0.16);
  p3->SetRightMargin(0.1);
  p3->SetBottomMargin(0.14);
  p3->SetTopMargin(0.1);

  cout << "checkpoint 2.1" << endl;

  DrawFrame(hSigFit1, STCol, STSty);
  cout << "checkpoint 2.2" << endl;
  DrawFrame(hSigFit1_clone_dummy, STCol, STSty);
  hSigFit1->SetMarkerSize(1.2);

  hSigFit1->GetXaxis()->SetRangeUser(0.469, 0.526);
  hSigFit1->GetYaxis()->SetRangeUser(0.001, 23.1e+7);

  hSigFit1_clone_dummy->GetXaxis()->SetRangeUser(0.469, 0.526);
  hSigFit1_clone_dummy->GetYaxis()->SetRangeUser(0.001, 23.1e+7);
  hSigFit1_clone->SetLineColor(0);
  hSigFit1_clone->SetLineWidth(0);

  //  hSigFit1->Draw("ep");
  hSigFit1_clone->SetFillColor(5);

  // hSigFit1->SetFillColor(42);

  //  hSigFit1_clone->GetXaxis()->SetRangeUser(0.48,0.51);
  //  hSigFit1_clone->GetXaxis()->SetRangeUser(0.488,0.507);
  hSigFit1_clone->GetXaxis()->SetRangeUser(0.482, 0.512);
  // hSigFit1_clone->GetXaxis()->SetRangeUser(0.482,0.512);
  hSigFit1_clone_dummy->GetXaxis()->SetRangeUser(0.469, 0.526);

  hSigFit1_clone_dummy->Draw("p");
  //   hSigFit1_clone->Draw("E3hist same");
  hSigFit1_clone->Draw("E3 hist same");
  hSigFit1->Draw("psame");

  cout << "checkpoint 3" << endl;

  // hSigFit1_clone->Draw("p same");

  // fSigOnly1->Draw("same l");
  //  fSigBg1->Draw("same");
  //  fBgOnly1->Draw("same l");

  TLegend *lp2 = DrawLegend(0.14, 0.64, 0.4, 0.88);
  lp2->SetTextSize(.04);
  lp2->SetFillStyle(0);
  lp2->AddEntry((TObject *)0, "#bf{ALICE Performance}", "");
  lp2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  lp2->AddEntry((TObject *)0, "V0M, 0-100%", "");
  lp2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
  lp2->Draw("same");

  // rotational
  TLegend *lp3 = DrawLegend(0.50, 0.68, 0.80, 0.82);
  lp3->SetFillStyle(0);
  //  lp3->AddEntry((TObject*)0,Form("%s",ptbn.Data()),"");
  //  lp3->AddEntry(fSigBg1,"Fit Function","l");
  lp3->AddEntry(hSigFit1, "Data (stat. uncert.)", "p");
  lp3->AddEntry(fSigBg1, "Gaussian Peak Fit", "l");
  lp3->AddEntry(hSigFit1_clone, "Signal", "f");
  lp3->SetTextSize(.04);
  lp3->Draw("same");

  p3->Modified();

  cout << "checkpoint 4" << endl;

  // c11->cd();
  // c11->SaveAs("K0s_selection_performance.png");
  // c11->Print("K0s_selection_performance.eps");

  //  c11->Print("Invmass_Preliminary_signal_k0sk0s_pp13TeV_0to20.root");
  // c11->Print("Invmass_Preliminary_signal_k0sk0s_pp13TeV_0to20.pdf");
}

TCanvas *DrawCanvas(TString opt = "c")
{
  TCanvas *c1 = new TCanvas(opt.Data(), opt.Data(), 10, 10, 600, 600);
  c1->cd(1);
  c1->SetLeftMargin(0.2);
  c1->SetRightMargin(0.15);
  c1->SetTopMargin(0.08);
  c1->SetBottomMargin(0.15);
  return c1;
}
TPad *MyPad(Double_t x1 = 0.01, Double_t y1 = 0.01, Double_t x2 = 0.49, Double_t y2 = 0.49, Float_t lm = 0.2, Float_t rm = 0.05, Float_t tm = 0.05, Float_t bm = 0.2, TString name = "pad")
{
  TPad *c1_1 = new TPad(name.Data(), name.Data(), x1, y1, x2, y2);
  c1_1->Draw();
  c1_1->cd();
  c1_1->Range(0, 0, 1, 1);
  c1_1->SetBorderSize(2);
  c1_1->SetBorderMode(0);
  c1_1->SetFillColor(10);
  c1_1->SetFrameFillColor(10);
  c1_1->SetFrameLineWidth(2);
  c1_1->SetLeftMargin(lm);
  c1_1->SetRightMargin(rm);
  c1_1->SetTopMargin(tm);
  c1_1->SetBottomMargin(bm);
  c1_1->SetTicks(1, 1);
  return c1_1;
}
TGraphErrors *PlotGraph(Int_t NdataPoint, Int_t MarkerColor, Int_t MarkerStyle, Double_t *X, Double_t *ErX, Double_t *Y, Double_t *ErY)
{
  MarkerColor = 1;
  MarkerStyle = 20;
  TGraphErrors *gr = new TGraphErrors(NdataPoint, X, Y, ErX, ErY);
  gr->SetTitle("");
  gr->SetMarkerStyle(MarkerStyle);
  gr->SetMarkerColor(MarkerColor);
  gr->SetMarkerSize(1.5);
  // gr->GetXaxis()->SetTitle("m_{T} - m_{0} (GeV/c^{2})");
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleOffset(1.35);
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->CenterTitle(true);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetNdivisions(509);
  // gr->GetYaxis()->SetTitle("#frac{1}{2#pi m_{T}N_{Ev}}#frac{dN^{2}}{dydm_{T}} (GeV/c)^{-1}");
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->CenterTitle(true);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleOffset(1.75);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetNdivisions(509);
  return gr;
}
//==================================
TLegend *DrawLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  // TLegend *legend = new TLegend(0.5,0.65,0.88,0.85);
  TLegend *legend = new TLegend(x1, y1, x2, y2);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  // legend->AddEntry(gr1,"(0 - 100) %","p");
  // legend->AddEntry(func1,"p_{0}[ 1 + 2 v_{2}^{obs} cos(2(#Phi - #Psi))]","l");
  return legend;
}
//=====================================
TLatex *DrawText(Double_t x = 0, Double_t y = 0, Int_t tColor = 2, TString name = "")
{
  TLatex *tex = new TLatex(x, y, name.Data());
  tex->SetTextSize(0.04);
  tex->SetTextColor(tColor);
  tex->SetTextFont(42);
  // tex->Draw();
  return tex;
}
//==================================
void MyStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFillColor(0);
  gStyle->SetLineColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleColor(1);
  gStyle->SetLineWidth(2);
}

Double_t BreitWignerExpol(Double_t *x, Double_t *par)
{
  return par[0] * par[2] / (2 * 3.14159) / (pow((x[0] - par[1]), 2) + pow(par[2], 2) / 4.) + (pow((x[0] - 0.63718), par[6])) * exp(par[3] + x[0] * par[4] + x[0] * x[0] * par[5]);
}

Double_t Expol(Double_t *x, Double_t *par)
{
  return (pow((x[0] - 0.63718), par[3])) * exp(par[0] + x[0] * par[1] + x[0] * x[0] * par[2]);
}

Double_t VoigtPoly2(Double_t *x, Double_t *par)
{
  return 0.001 * fReBinSpectra * par[0] * TMath::Voigt(x[0] - par[2], par[3], par[1], 4) + (par[4] + par[5] * x[0] + par[6] * x[0] * x[0]);
}
Double_t Voigt(Double_t *x, Double_t *par)
{
  return 0.001 * fReBinSpectra * par[0] * TMath::Voigt(x[0] - par[2], par[3], par[1], 4);
}
Double_t Poly2(Double_t *x, Double_t *par)
{
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}
Double_t BWPoly2(Double_t *x, Double_t *par)
{
  return (par[0] * par[1] * 0.01 * fReBinSpectra) / (2 * 3.14159) / ((x[0] - par[2]) * (x[0] - par[2]) + (par[1] * par[1]) / 4.) + par[3] + par[4] * x[0] + par[5] * x[0] * x[0];
}

Double_t BW(Double_t *x, Double_t *par)
{
  return (1.0 * 0.01 * fReBinSpectra / (2 * 3.14159)) * ((par[0] * par[1]) / ((x[0] - par[2]) * (x[0] - par[2]) + (par[1] / 2) * (par[1] / 2)));
}

// Double_t Guass(Double_t *x, Double_t *par)
//{
//  return (par[0]*(TMath::Exp(-0.5*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]))));
// }

Double_t Guass(Double_t *x, Double_t *par)
{
  return (par[0] * TMath::Exp(-0.5 * pow(((x[0] - par[1]) / par[2]), 2)));
}

void DrawFrame(TH1 *h, Int_t MCol, Int_t MSty)
{ // h->GetXaxis()->SetTitle("#it{M}_{K^{0}_{S}#pi^{#pm}} (GeV/#it{c^{2}})");
  // h->GetXaxis()->SetTitle("#it{M}_{K^{0}_{S}K^{0}_{S}} (GeV/#it{c^{2}})");
  h->GetXaxis()->SetTitle("#it{M}_{#pi^{+}#pi^{-}} (GeV/#it{c^{2}})");
  h->SetMarkerColor(MCol);
  h->SetMarkerStyle(MSty);
  h->SetMarkerSize(0.5);
  h->SetLineColor(MCol);
  h->SetLineWidth(2);
  // h->GetXaxis()->CenterTitle(true);
  h->GetXaxis()->SetNdivisions(509);
  h->GetYaxis()->SetNdivisions(509);
  h->GetXaxis()->SetLabelOffset(0.015);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleFont(42);
  // h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTickLength(0.04);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleOffset(2.2);
  h->GetYaxis()->CenterTitle(true);
  h->GetXaxis()->CenterTitle(true);
  // h->GetYaxis()->SetDecimals(true);
  // h->GetYaxis()->SetNdivisions(310);
  h->GetYaxis()->SetLabelOffset(0.035);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTickLength(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitle("Counts /(1 MeV/#it{c^{2}})");
  h->GetYaxis()->SetTitleFont(42);
  // h->Draw("");
  // return h;
}
