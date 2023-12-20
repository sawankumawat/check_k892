

/////////////////////////////ALL IN ONE MACRO USING MIXED EVENT COMBINATORIAL BKG////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "TArrow.h"

using namespace std;

void SetHistoStyle(TH1 *h, Int_t mcolor, Int_t mstyle, Float_t msize, Float_t Tsizex, Float_t Tsizey, Float_t Lsizex, Float_t Lsizey, Float_t Offsetx, Float_t Offsety)
{
  h->SetMarkerColor(mcolor);
  h->SetMarkerStyle(mstyle);
  h->SetMarkerSize(msize);
  h->GetXaxis()->SetTitleSize(Tsizex);
  h->GetYaxis()->SetTitleSize(Tsizey);
  h->GetXaxis()->SetLabelSize(Lsizex);
  h->GetYaxis()->SetLabelSize(Lsizey);
  h->GetXaxis()->SetTitleOffset(Offsetx);
  h->GetYaxis()->SetTitleOffset(Offsety);
}

//******Different fitting functions*****************************************************************************************

Double_t BreitWignerpoly2(Double_t *x, Double_t *par)
{
  double BW = (0.5 * par[2] * par[1] / TMath::Pi()) / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]);
  // double poly = par[5] + par[4]*x[0] + par[3]*x[0]*x[0];
  double poly3 = par[6] + par[5] * x[0] + par[4] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
  return (BW + poly3);
}

Double_t BW(Double_t *x, Double_t *par)
{
  return (0.5 * par[2] * par[1] / TMath::Pi() / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]));
}

Double_t BreitWignerpoly1(Double_t *x, Double_t *par)
{
  double BW = (0.5 * par[2] * par[1] / 3.14159) / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]);
  double poly1 = par[4] + par[3] * x[0];
  // double poly3 = par[0] + par[1]*x[0] + par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0] ;
  return (BW + poly1);
}

Double_t polynomial2(Double_t *x, Double_t *par)
{
  double poly1 = par[2] + par[1] * x[0] + par[0] * x[0] * x[0];
  return (poly1);
}

Double_t polynomial1(Double_t *x, Double_t *par)
{
  return (par[0] + par[1] * x[0]);
}

Double_t polynomial3(Double_t *x, Double_t *par)
{
  double poly3 = par[3] + par[2] * x[0] + par[1] * x[0] * x[0] + par[0] * x[0] * x[0] * x[0];
  return (poly3);
}

Double_t BWExpo(Double_t *x, Double_t *par)
{
  // return (1.0/(2*3.14159))*((par[2]*par[1])/((x[0]-par[0])*(x[0]-par[0]) +(par[1]/2)*(par[1]/2))) + TMath::Power((x[0]-(0.13957+0.49367)),par[6])*exp(par[5] + par[4]*x[0]+par[3]*x[0]*x[0]);

  // return par[0]*par[2]/(2*3.14159)/((x[0]-par[1])**2+par[2]**2/4.)+ ((x[0]- 0.63718)**par[6])*exp(par[3]+x[0]*par[4]+ x[0]*x[0]*par[5]);
  // return par[2]*par[1]/(2*3.14159)/((x[0]-par[0])**2+par[1]**2/4.)+ ((x[0]- 0.63718)**par[6])*exp(x[0]*x[0]*par[3]+x[0]*par[4]+par[5]);
  // return par[2]*par[1]/(2*3.14159)/((x[0]-par[0])**2+par[1]**2/4.)+ ((x[0]- 0.63718)**par[6])*exp(par[3]+x[0]*par[4]+ x[0]*x[0]*par[5]);

  return par[2] * par[1] / (2 * 3.14159) / (pow((x[0] - par[0]), 2) + pow(par[1], 2) / 4.) + (pow((x[0] - 0.63718), par[6])) * exp(par[3] + x[0] * par[4] + x[0] * x[0] * par[5]);
}

Double_t Expo(Double_t *x, Double_t *par)
{
  // return ((x[0]- 0.63718)**par[3])*exp(par[0] + x[0]*par[1] + x[0]*x[0]*par[2]);
  return (pow((x[0] - 0.63718), par[3])) * exp(par[0] + x[0] * par[1] + x[0] * x[0] * par[2]);
  // return TMath::Power((x[0]-(0.13957+0.49367)),par[3])*exp(par[2] + par[1]*x[0]+par[0]*x[0]*x[0]);
}

//******************************************************************************************************************************************

//*****TLegend Class******************************************************************

TLegend *DrawLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{

  TLegend *legend = new TLegend(x1, y1, x2, y2);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  return legend;
}

//*************************************************************************************

void allinonemacromix()

{

  TLatex *t2 = new TLatex();
  t2->SetNDC(); // to self adjust the text so that it remains in the box
  t2->SetTextSize(0.04);

  // const Int_t Npt=14;  //no. of pt bins(max)
  const Int_t Npt = 32; // no. of pt bins(max)

  //** For projection of required signals in different pt bins used as array elements**************************

  double Low_pt[Npt] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.5, 6.0, 7.0, 8.0, 10, 12};
  double High_pt[Npt] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.5, 6.0, 7.0, 8.0, 10, 12, 15};

  double x[Npt] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 4.7, 4.9, 5.25, 5.75, 6.5, 7.5, 9, 11, 13.5};
  //*****************************************************************************************************************

  //**Total no. of pt bins as an array and defining different histograms for extracting fitting parameters********************************************************************************************************************
  // double NPT[Npt+1] ={0.4,0.8,1.2,1.6,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0,16.0};
  double NPT[Npt + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.5, 6.0, 7.0, 8.0, 10, 12, 15};

  TH1F *hYieldpar = new TH1F("hYieldpar", "hYieldpar", Npt, NPT);                   // for ptspectra from fitting parameter directly
  TH1F *hintegral_yield = new TH1F("hintegral_yield", "hintegral_yield", Npt, NPT); // pt spectra from function integration
  TH1F *hChiSquare = new TH1F("hChiSquare", "hChiSquare", Npt, NPT);                // for chisquare
  TH1F *hsgnfcance = new TH1F("hsgnfcance", "hsgnfcance", Npt, NPT);                // for significance of signal
  TH1F *hmass = new TH1F("hmass", "hmass", Npt, NPT);                               // for mass from fit
  TH1F *hwidth = new TH1F("hwidth", "hwidth", Npt, NPT);                            // for width from fit
  TH1F *hYbincount = new TH1F("hYbincount", "hYbincount", Npt, NPT);                // Yield calculation using bin counting method
  TH1F *hFrac_stat_error = new TH1F("hFrac_stat_error", "hFrac_stat_error", Npt, NPT);
  TH1F *herrormass = new TH1F("herrormass", "", Npt, 0.5, 16);   // for error band in mass
  TH1F *herrorwidth = new TH1F("herrorwidth", "", Npt, 0.5, 16); // for error band in width

  ////Setting bin content for mass and width error/////////////////////////////////////////////////
  for (int i = 1; i <= Npt; i++)
  {
    herrormass->SetBinContent(i, 0.895);
    herrormass->SetBinError(i, 0.0015);
    herrorwidth->SetBinContent(i, 0.047);
    herrorwidth->SetBinError(i, 0.005);
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////

  //*************************************************************************************************************

  //**For normalisation **************************************************************************************

  // Right normalization
  double lownorm[Npt + 5] = {0.65, 0.65, 0.65, 0.65, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
  double highnorm[Npt + 5] = {0.70, 0.70, 0.70, 0.70, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

  // double lownorm[Npt + 5] = {1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1};
  // double highnorm[Npt + 5] = {1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

  // double lownorm[Npt + 5] = {1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05};
  // double highnorm[Npt + 5] = {1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15};

  // Left-normalization
  //  double lownorm[Npt + 5] = {0.68, 0.68,0.68, 0.68,0.68, 0.68,0.68, 0.68,0.68, 0.68,0.68, 0.68,0.68, 0.68,0.68, 0.68,0.68, 0.68,0.68, 0.68,0.68, 0.68,0.68, 0.68, 0.68};
  //  double highnorm[Npt + 5] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};

  //**********************************************************************************************************

  //**For fitting ranges in different pt bins******************************************************************************

  // For ME free width
  // double lowfitrange[Npt + 1] = {0.66, 0.69, 0.65, 0.72, 0.74, 0.745, 0.738, 0.75, 0.81, 0.72, 0.72, 0.68, 0.70, 0.74, 0.74, 0.71, 0.75, 0.76, 0.73, 0.74, 0.7, 0.69, 0.7, 0.67, 0.67, 0.67, 0.67, 0.67, 0.67, 0.67, 0.67, 0.67};

  // double highfitrange[Npt + 1] = {1.18, 1.21, 1.27, 1.15, 1.28, 1.03, 1.082, 1.20, 1.2, 1.08, 1.1, 1.1, 1.1, 1.1, 1.15, 1.12, 1.1, 1.1, 1.14, 1.1, 1.15, 1.08, 1.2, 1.1, 1.11, 1.11, 1.08, 1.08, 1.1, 1.15, 1.15, 1.15};

  // For ME fixed width
  //  double lowfitrange[Npt + 1] = {0.74, 0.72, 0.72, 0.72, 0.73, 0.745, 0.78, 0.79, 0.795, 0.785, 0.72, 0.68, 0.66, 0.74, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7};

  // double highfitrange[Npt + 1] = {1.15, 1.1, 1.1, 1.05, 1.2, 1.03, 1.04, 1.08, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1};

  // For LS free width
  double lowfitrange[Npt + 1] = {0.735, 0.70, 0.71, 0.725, 0.729, 0.74, 0.785, 0.79, 0.73, 0.785, 0.8, 0.75, 0.75, 0.745, 0.75, 0.76, 0.755, 0.75, 0.75, 0.735, 0.75, 0.74, 0.732, 0.745, 0.743, 0.74, 0.74, 0.74, 0.74, 0.74, 0.74, 0.74};

  double highfitrange[Npt + 1] = {1.2, 1.2, 1.2, 1.25, 1.25, 1.19, 1.02, 1.08, 1.1, 1.09, 1.1, 1.14, 1.14, 1.12, 1.12, 1.1, 1.12, 1.15, 1.12, 1.15, 1.2, 1.2, 1.18, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25};

  // For LS fixed width
  // double lowfitrange[Npt + 1]  = {0.735, 0.70, 0.71, 0.725, 0.729, 0.74, 0.785, 0.79, 0.73, 0.785, 0.8, 0.75, 0.75, 0.745, 0.75, 0.76, 0.755, 0.75, 0.75, 0.735, 0.75, 0.74, 0.732, 0.755, 0.743};

  // double highfitrange[Npt + 1] = {1.05 , 1.2 , 1.2,  1.25,  1.25,   1.19, 1.02, 1.08,  1.1,  1.09,  1.1, 1.14, 1.14, 1.12,  1.12, 1.1,  1.12, 1.15,  1.12, 1.15,  1.2,  1.2,  1.18,  1.25,  1.25};

  //************************************************************************************************************************

  //**For different centrality in different pt bins *********************************************************************

  double Low_Centr[Npt + 1] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double High_Centr[Npt + 1] = {100, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90};

  //**********************************************************************************************************************

  //**PDG Mass and Width of required signal for significance calculation ******************************************

  float masspdg = 0.895;
  float widthpdg = 0.047; // in 1 sigma

  //**************************************************************************

  //**Initialisation values for yield of BW and coefficients of POL2 and POL3***************************************************************

  Double_t Para2[Npt + 5] = {1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04};

  // 0-10

  Double_t RBPar3L[Npt + 5] = {-100, -100, -100, -100, -100, -100, -20, -100, -100, -100, -20, -20, -20, -20};
  Double_t RBPar3H[Npt + 5] = {100, 100, 100, 100, 100, 100, 200, 100, 100, 100, 200, 200, 200, 200};

  Double_t RBPar4L[Npt + 5] = {-100, -100, -100, -100, -100, -100, -20, -100, -100, -100, -20, -20, -20, -20};
  Double_t RBPar4H[Npt + 5] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100};

  Double_t RBPar5L[Npt + 5] = {-200, -200, -200, -200, -200, -200, -50, -200, -200, -200, -50, -50, -50, -50};
  Double_t RBPar5H[Npt + 5] = {200, 200, 200, 200, 200, 200, 30, 200, 200, 200, 30, 30, 30, 30};

  //****************************************************************************************************************************************

  //**Variables and arrays needed *******************************************************************************************

  Double_t integralsignalfunc[Npt]; // for calculation of area under signal only (after fitting and extraction BW signal)
  Double_t ptcenter[Npt];           // stores mean pt for each pt bin
  Double_t ptbinwidth[Npt];         // stores pt binwidth
  // TGraphErrors *g=new TGraphErrors();
  Double_t interror[100];                                   // stores error while area calculation of signal for yield
  const int N = 1;                                          // for centrality loop below used
  Double_t dy = 1.0;                                        // for rapidity difference
  TFitResultPtr r;                                          // for fitting using TFitter
  TH1D *fHistTotal[Npt];                                    // for sig+bg
  TH1D *fHistBkg[Npt];                                      // for mixedbg
  TH1D *fHistbkgLS[Npt];                                    // for like sign
  TH1D *fHistlike[Npt];                                     // for resultant like bg
  TH1D *histPP[Npt];                                        // for like pp bg
  TH1D *histMM[Npt];                                        // for like mm bg
  TH1D *hfbkg;                                              // for normalised mixed bkg
  TH1D *hfbkg_like;                                         // for resultant like bkg
  char name[100];                                           // for giving name to canvases and histograms
  Double_t lowpt, highpt;                                   // for extracting pt from the above defined pt array
  int lbin, hbin;                                           // corresponding bin for a given pt
  int rebin = 2;                                            // for rebinning  //=1 for no rebin
  int lownormbin1, lownormbin2, highnormbin1, highnormbin2; // normalisation bins when using mixed background
  Double_t sigbkg_integral, bkg_integral, normfactor;       // integrals and normalisation factor(scaling factor) using above normalisation bins

  // details of input root file////////////////////
  int Nbins_file = 900;
  float lowrange_file = 0.6;
  float highrange_file = 1.5;
  float binwidth_file = 0.02;
  ///////////////////////////////////////////////

  TString ResBkg = "LIKE"; // residual bkg to be used  //LIKE or MIX

  // fiiting parameters array///////////////////////
  Double_t Yield[Npt];
  Double_t Mass[Npt];
  Double_t Width[Npt];
  Double_t poly0[Npt];
  Double_t poly1[Npt];
  Double_t poly2[Npt];
  Double_t poly3[Npt];
  Double_t ErrorMass[Npt];
  Double_t ErrorWidth[Npt];
  Double_t ErrorYield[Npt];
  Double_t Chi2Ndf[Npt];
  /////////////////////////////////////////////////

  Double_t yieldcalc, yielderror;            // calculation of raw yield from function integration
   Double_t yieldcalc1[Npt], yielderror1[Npt]; 
  Double_t Yield_value_par, Yield_error_par; // calculation of raw yield from fitting parameter
  Double_t BR = 0.33;                        // branching ratio

  Double_t significance_den, significance_num, ratio, ratio2; // calculate signal+bkg and only signal integral respectively and then significance(ratio)
  int bmin, bmax;
  Double_t hBCError_1, bkgvalue, Integral_BW_withsigma, fYield_BinCount, YieldIntegral_BW, Yfraction_cBW, sum_tail_correction, Total_Ybincounting, Tail_correction_plusm, Tail_correction_minusm, Error_2, Final_pro_error;
  Double_t nlow, nhigh;
  Double_t Yield_bincount_hist;

  // float al=0.65;
  // float bh=1.2;

  float al = 0.65;
  // float bh=1.05;
  float bh = 1.3;
  //***************************************************************************************************************************

  //**Canvas definitions and initialisations********************************************************************************

  TCanvas *cgrid1 = new TCanvas();
  TCanvas *cgrid2 = new TCanvas();
  TCanvas *cgrid3 = new TCanvas();
  TCanvas *cgrid4 = new TCanvas();

  TCanvas *cinv[Npt]; // for output canvases on screen containing fitted signal after subtraction

  for (Int_t ip = 0; ip < Npt; ip++)
  {
    TString cName = TString::Format("cinv_pt_%2.1f-%2.1f", Low_pt[ip], High_pt[ip]);
    cinv[ip] = new TCanvas(Form("cinv%d", ip), cName.Data(), 10, 10, 1000, 800);
    cinv[ip]->SetLeftMargin(0.15);
    cinv[ip]->SetRightMargin(0.05);
    // cinv[ip]->SetBottomMargin(0.2);
    cinv[ip]->SetTopMargin(0.05);
    cinv[ip]->SetBottomMargin(0.15);
    cinv[ip]->SetBorderMode(0);
    cinv[ip]->SetBorderSize(2);
    cinv[ip]->SetFrameLineWidth(2);
    cinv[ip]->SetFillColor(0);
    cinv[ip]->SetTicks(1, 1);
  }

  TCanvas *c2 = new TCanvas("SigAfterSub", "SigAfterSub", 10, 10, 1000, 800); // single canvas containing fitted signal after subtraction
  c2->Divide(5, 5);
  TCanvas *c22 = new TCanvas("c22", "", 10, 10, 600, 600);
  c22->Divide(5, 5);

  TCanvas *c1 = new TCanvas("SigBeforeSub", "SigBeforeSub", 10, 10, 1000, 800); // single canvas containing signal with bkg(after norm in case of mix)
  c1->Divide(5, 5);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  // c1->SetBottomMargin(0.2);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.1);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  c1->SetFillColor(0);
  c1->SetTicks(1, 1);

  TCanvas *c11 = new TCanvas("c11", "", 10, 10, 1000, 800);
  c11->Divide(5, 5);
  c11->SetLeftMargin(0.15);
  c11->SetRightMargin(0.05);
  // c11->SetBottomMargin(0.2);
  c11->SetTopMargin(0.05);
  c11->SetBottomMargin(0.1);
  c11->SetBorderMode(0);
  c11->SetBorderSize(2);
  c11->SetFrameBorderMode(0);
  c11->SetFillColor(0);
  c11->SetTicks(1, 1);

  TCanvas *cSigbkg[Npt]; // for output canvases on screen containing signal with bkg(after norm in case of mix)

  for (Int_t ip = 0; ip < Npt; ip++)
  {
    TString cNam = TString::Format("cSigbkg_pt_%2.1f-%2.1f", Low_pt[ip], High_pt[ip]);
    cSigbkg[ip] = new TCanvas(Form("cSigbkg%d", ip), cNam.Data(), 10, 10, 1000, 800);
    cSigbkg[ip]->SetLeftMargin(0.15);
    cSigbkg[ip]->SetRightMargin(0.05);
    // cSigbkg[ip]->SetBottomMargin(0.2);
    cSigbkg[ip]->SetTopMargin(0.08);
    cSigbkg[ip]->SetBottomMargin(0.15);
    cSigbkg[ip]->SetBorderMode(0);
    cSigbkg[ip]->SetBorderSize(2);
    cSigbkg[ip]->SetFrameBorderMode(0);
    cSigbkg[ip]->SetFillColor(0);
    cSigbkg[ip]->SetTicks(1, 1);
  }

  TCanvas *cfit[Npt]; // loosely defined

  for (Int_t ip = 0; ip < Npt; ip++)
  {
    TString cName = TString::Format("cfit_pt_%2.1f-%2.1f", Low_pt[ip], High_pt[ip]);
    cfit[ip] = new TCanvas(Form("cfit%d", ip), cName.Data(), 10, 10, 1000, 800);
    cfit[ip]->SetLeftMargin(0.15);
    cfit[ip]->SetRightMargin(0.05);
    // cfit[ip]->SetBottomMargin(0.2);
    cfit[ip]->SetBottomMargin(0.1);
    cfit[ip]->SetTopMargin(0.05);
    cfit[ip]->SetBorderMode(0);
    cfit[ip]->SetBorderSize(2);
    cfit[ip]->SetFrameBorderMode(0);
    cfit[ip]->SetFillColor(0);
    cfit[ip]->SetTicks(1, 1);
  }

  //**********************************************************************************************************************************

  //**Histogram arrays initialisation*************************************************************************************************

  for (Int_t m = 0; m < Npt; m++)
  {
    sprintf(name, "fHistNum%d", m);
    fHistTotal[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
    sprintf(name, "fHistbkg%d", m);
    fHistBkg[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
    sprintf(name, "histPP_pt%d", m);
    histPP[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
    sprintf(name, "histMM_pt%d", m);
    histMM[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
    sprintf(name, "histlike_pt%d", m);
    fHistlike[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
  }

  //***************************************************************************************************************************************

  //**For accessing the histograms stored within root file*****************************************************************************

  TFile *fInputFile = new TFile("../data/LHC22m_pass4_temp/AnalysisResults.root", "Read"); // Reading the File
  TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(3);
  TList *fInputList = (TList *)fInputFile->Get(key->GetName());

  TFile *file2 = new TFile("../HEPData-ins1797443-v1-root.root", "READ");
  TKey *key2 = (TKey *)file2->GetListOfKeys()->At(3);
  TGraph *grsourav = (TGraph *)file2->Get("Table 4/Graph1D_y1");
  grsourav->SetTitle(0);
  grsourav->GetXaxis()->SetTitle(0);
  grsourav->GetYaxis()->SetTitle(0);

  const string kOutputName = "lf-k892analysis";
  const string kOutputName2 = "QAbefore";
  const string kOutputName3 = "QAafter";

  //**To calculate total number of events for which histograms were filled*************************************************************
  // TH1F* hEVent = (TH1F *) fInputList->FindObject("hAEventsVsMulti");
  // Double_t Event=hEVent->Integral(lc,hc);
  Double_t Event = 4.2276e8;
  //*************************************************************************************************************************

  cout << "*****************number of events********************:" << Event << endl;

  //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************
  TH3F *fHistNum = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassDS", kOutputName.c_str()));
  TH3F *fHistDen = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassME", kOutputName.c_str()));
  TH3F *fHistLS = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassLS", kOutputName.c_str()));

  // // QA before

  TH1F *histtrkpT_pi = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_pi", kOutputName.c_str(), kOutputName2.c_str()));
  TH1F *histtrkpT_ka = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_ka", kOutputName.c_str(), kOutputName2.c_str()));
  TH2F *histTofTpc = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_TPC_Map_pi_all", kOutputName.c_str(), kOutputName2.c_str()));
  TH2F *histTofpi = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_pi_all", kOutputName.c_str(), kOutputName2.c_str()));
  TH2F *histTpcpi = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigma_pi_all", kOutputName.c_str(), kOutputName2.c_str()));
  TH2F *histTofka = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_ka_all", kOutputName.c_str(), kOutputName2.c_str()));
  TH2F *histTpcka = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigmaka_all", kOutputName.c_str(), kOutputName2.c_str()));

  // // QA after

  TH1F *histtrkpT_pi2 = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_pi", kOutputName.c_str(), kOutputName3.c_str()));
  TH1F *histtrkpT_ka2 = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_ka", kOutputName.c_str(), kOutputName3.c_str()));
  TH2F *histTofTpc2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_TPC_Map_pi_all", kOutputName.c_str(), kOutputName3.c_str()));
  TH2F *histTofpi2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_pi_all", kOutputName.c_str(), kOutputName3.c_str()));
  TH2F *histTpcpi2 = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigma_pi_all", kOutputName.c_str(), kOutputName3.c_str()));
  TH2F *histTofka2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_ka_all", kOutputName.c_str(), kOutputName3.c_str()));
  TH2F *histTpcka2 = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigmaka_all", kOutputName.c_str(), kOutputName3.c_str()));

  //**********************************************************************************************************************************


  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1110);
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
  // gStyle->SetErrorX(0);

  //*********************************


  for (Int_t ip = 0; ip < Npt; ip++) // start pt bin loop
  // for (Int_t ip = 1; ip < 4; ip++)
  {
    lowpt = Low_pt[ip];
    highpt = High_pt[ip];

    fHistTotal[ip] = fHistNum->ProjectionZ(Form("hSig_%d", ip), -1, -1, fHistNum->GetYaxis()->FindBin(lowpt), fHistNum->GetYaxis()->FindBin(highpt), "E");
    fHistBkg[ip] = fHistDen->ProjectionZ(Form("hbkg_%d", ip), -1, -1, fHistDen->GetYaxis()->FindBin(lowpt), fHistDen->GetYaxis()->FindBin(highpt), "E");
    fHistbkgLS[ip] = fHistLS->ProjectionZ(Form("hbkgLS_%d", ip), -1, -1, fHistLS->GetYaxis()->FindBin(lowpt), fHistLS->GetYaxis()->FindBin(highpt), "E");

    fHistTotal[ip]->Rebin(rebin);
    fHistBkg[ip]->Rebin(rebin);
    fHistbkgLS[ip]->Rebin(rebin);

    //**Cloning sig+bkg histogram for like sign or mixed event subtraction *********************************************************
    TH1D *hfsig = (TH1D *)fHistTotal[ip]->Clone();
    //*****************************************************************************************************************************

    //**If mixed bkg is used *******************************************************************************

    if (ResBkg.CompareTo("MIX") == 0)
    {

      Double_t init = lownorm[ip];
      Double_t fin = highnorm[ip];

      TAxis *axis1 = fHistTotal[ip]->GetXaxis();
      lownormbin1 = axis1->FindBin(init);
      highnormbin1 = axis1->FindBin(fin);
      sigbkg_integral = (fHistTotal[ip]->Integral(lownormbin1, highnormbin1));

      TAxis *axis2 = fHistBkg[ip]->GetXaxis();
      lownormbin2 = axis2->FindBin(init);
      highnormbin2 = axis2->FindBin(fin);
      bkg_integral = (fHistBkg[ip]->Integral(lownormbin2, highnormbin2));

      normfactor = sigbkg_integral / bkg_integral; // scaling factor for mixed bkg
      hfbkg = (TH1D *)fHistBkg[ip]->Clone();

      hfbkg->Scale(normfactor);

      hfsig->Add(hfbkg, -1);
    }
    else if (ResBkg.CompareTo("MIX") != 0)
    {
      hfbkg_like = (TH1D *)fHistbkgLS[ip]->Clone();
      hfsig->Add(hfbkg_like, -1);
    }

    //****pt bincenter and pt binwidth**********************************************************************

    ptcenter[ip] = (Low_pt[ip] + High_pt[ip]) / 2;
    ptbinwidth[ip] = High_pt[ip] - Low_pt[ip];

    //******************************************************************************************************

    TF1 *fitFcn = new TF1("fitfunc", BreitWignerpoly2, lowfitrange[ip], highfitrange[ip], 7); // sig+bkg fit function
    TF1 *fitFcn1 = new TF1("fitfunc1", polynomial3, lowfitrange[ip], highfitrange[ip], 4);    // only residualbkg
    TF1 *fitFcn2 = new TF1("fitFcn2", BW, lowfitrange[ip], highfitrange[ip], 3);              // only signal

    fitFcn->SetParLimits(0, 0.83, 0.92); // Mass
    fitFcn->SetParLimits(2, 0, 100000);  // Yield
    fitFcn->SetParameter(1, 0.047);      // width
    // fitFcn->FixParameter(1, 0.047);      // width
    fitFcn->SetParameter(2, 1000); // yield

    fitFcn->SetParameter(3, 100);
    fitFcn->SetParameter(4, 100);
    fitFcn->SetParameter(5, 100);
    fitFcn->SetParameter(6, 100);

    /*
    fitFcn->SetParLimits(3,-200,200);
    fitFcn->SetParLimits(4,-200,200);
    fitFcn->SetParLimits(5,-200,200);
    */

    // fitFcn->SetParNames("Mass","Width","Yield","A","B","C","D");
    fitFcn->SetParNames("Mass", "Width", "Yield", "C", "B", "A");
    // fitFcn->SetParNames("Yield","Mass","Width","C","B","A","D");
    r = hfsig->Fit(fitFcn, "RS+"); // signal after bkg subtraction

    //****************************************************************************************************************

    //**Extraction of fitting parameters******************************************************************************

    Double_t *par = fitFcn->GetParameters();

    Mass[ip] = fitFcn->GetParameter(0);
    Width[ip] = fitFcn->GetParameter(1);
    Yield[ip] = fitFcn->GetParameter(2);
    /*poly0[ip]=fitFcn->GetParameter(3);
    poly1[ip]=fitFcn->GetParameter(4);
    poly2[ip]=fitFcn->GetParameter(5);*/
    poly2[ip] = fitFcn->GetParameter(3);
    poly1[ip] = fitFcn->GetParameter(4);
    poly0[ip] = fitFcn->GetParameter(5);
    poly3[ip] = fitFcn->GetParameter(6);

    fitFcn2->SetParameters(&par[0]);
    fitFcn1->SetParameters(&par[3]);

    ErrorMass[ip] = fitFcn->GetParError(0);
    ErrorWidth[ip] = fitFcn->GetParError(1);
    ErrorYield[ip] = fitFcn->GetParError(2);
    Chi2Ndf[ip] = (fitFcn->GetChisquare()) / (fitFcn->GetNDF());

    //******************************************************************************************************************

    //**ERROR BIN COUNTING METHOD CALCULATION*****************************************************************************

    TF1 *fitFcn2_plusm = new TF1("fitFcn2_plusm", BW, lowfitrange[ip], highfitrange[ip], 3);
    TF1 *fitFcn2_minusm = new TF1("fitFcn2_minusm", BW, lowfitrange[ip], highfitrange[ip], 3);
    fitFcn2_plusm->FixParameter(0, Mass[ip] + ErrorMass[ip]);
    fitFcn2_plusm->FixParameter(1, 0.047);
    fitFcn2_plusm->FixParameter(2, Yield[ip]);

    fitFcn2_minusm->FixParameter(0, Mass[ip] - ErrorMass[ip]);
    fitFcn2_minusm->FixParameter(1, 0.047);
    fitFcn2_minusm->FixParameter(2, Yield[ip]);

    //*********************************************************************************************************************

    //**Calculation of significance and storing chi2 and sig in respective histograms*****************************************************

    bmin = hfsig->GetXaxis()->FindBin(masspdg - 2 * widthpdg);
    bmax = hfsig->GetXaxis()->FindBin(masspdg + 2 * widthpdg);

    significance_den = TMath::Sqrt(fHistTotal[ip]->Integral(bmin, bmax));
    significance_num = (fitFcn2->Integral(masspdg - 2 * widthpdg, masspdg + 2 * widthpdg)) / (binwidth_file);

    ratio = significance_num / significance_den; // significance of signal

    hsgnfcance->SetBinContent(ip + 1, ratio);
    hChiSquare->SetBinContent(ip + 1, Chi2Ndf[ip]); // storing both significance and chi2 in histogram

    //*****************************************************************************************************************************************

    //**Calculation of Yield using bin counting method and storing it in histogram***********************************************************

    Yield_bincount_hist = hfsig->IntegralAndError(bmin, bmax, hBCError_1);
    bkgvalue = fitFcn1->Integral(hfsig->GetBinLowEdge(bmin), hfsig->GetBinLowEdge(bmax + 1));
    Integral_BW_withsigma = fitFcn2->Integral(hfsig->GetBinLowEdge(bmin), hfsig->GetBinLowEdge(bmax + 1));
    fYield_BinCount = Yield_bincount_hist - (bkgvalue / binwidth_file);
    YieldIntegral_BW = fitFcn2->Integral(0.635, 5) / binwidth_file;
    Yfraction_cBW = (Integral_BW_withsigma / YieldIntegral_BW);

    sum_tail_correction = (fitFcn2->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + fitFcn2->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / binwidth_file;

    nlow = (fitFcn2->Integral(0.635, hfsig->GetBinLowEdge(bmin))) / binwidth_file;
    nhigh = (fitFcn2->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / binwidth_file;
    nlow = nlow / (Event * ptbinwidth[ip] * dy * BR * 2 * ptcenter[ip]);
    nhigh = nhigh / (Event * ptbinwidth[ip] * dy * BR * 2 * ptcenter[ip]);

    Total_Ybincounting = (sum_tail_correction + fYield_BinCount) / (Event * ptbinwidth[ip] * dy * BR * 2 * ptcenter[ip]);

    cout << "***************************************************************" << endl;
    cout << "****fraction of nlow for bin***********:"
         << " " << ip << " " << nlow / Total_Ybincounting << endl;
    cout << "****fraction of nhigh for bin***********:"
         << " " << ip << " " << nhigh / Total_Ybincounting << endl;
    cout << "***************************************************************" << endl;
    Tail_correction_plusm = (fitFcn2_plusm->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + (fitFcn2_plusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5))) / binwidth_file;
    Tail_correction_minusm = ((fitFcn2_minusm->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + fitFcn2_minusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / binwidth_file);
    Error_2 = sum_tail_correction - Tail_correction_plusm;
    Final_pro_error = TMath::Sqrt(Error_2 * Error_2 + hBCError_1 * hBCError_1) / (Event * ptbinwidth[ip] * dy * BR * 2 * ptcenter[ip]);

    ////Uncorrected Yield/////////////////////////////////////////////////////////////////////////////////

    hYbincount->SetBinContent(ip + 1, Total_Ybincounting);
    hYbincount->SetBinError(ip + 1, Final_pro_error);
    cout << "--------Total Value from bin counting----------" << (sum_tail_correction + fYield_BinCount) << endl;
    cout << "--------Value from bin counting----------" << Total_Ybincounting << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    // Fractional stat error///////////////////////////////////////////////////////////////////////////////

    hFrac_stat_error->SetBinContent(ip + 1, Final_pro_error / Total_Ybincounting);
    cout << "--------Frac error from bin counting----------" << (Final_pro_error / Total_Ybincounting) << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    //****************************************************************************************************************************************

    //**Calculation for raw pt spectra using function integration and filling it in histogram*********************************************

    integralsignalfunc[ip] = (fitFcn2->Integral((masspdg - 5 * widthpdg), (masspdg + 5 * widthpdg)));
    TMatrixDSym cov = r->GetCovarianceMatrix();
    TMatrixDSym cov1;
    TMatrixDSym cov2;
    cov.GetSub(0, 2, 0, 2, cov1);
    cov.GetSub(3, 5, 3, 5, cov2);
    Double_t *b = cov1.GetMatrixArray();
    Double_t *a = cov2.GetMatrixArray();
    Double_t *para = fitFcn->GetParameters();
    interror[ip] = fitFcn2->IntegralError((masspdg - 5 * widthpdg), (masspdg + 5 * widthpdg), &para[0], b);

    yieldcalc = integralsignalfunc[ip] / (Event * ptbinwidth[ip] * dy * BR * 2 * binwidth_file * ptcenter[ip]); // raw yield calculation
    yieldcalc1[ip] = integralsignalfunc[ip] / (Event * ptbinwidth[ip] * dy * BR * 2 * binwidth_file * ptcenter[ip]); // raw yield calculation
    yielderror = interror[ip] / (Event * ptbinwidth[ip] * dy * BR * 2 * binwidth_file * ptcenter[ip]);          // raw yield error
    yielderror1 [ip]= interror[ip] / (Event * ptbinwidth[ip] * dy * BR * 2 * binwidth_file * ptcenter[ip]);          // raw yield error

    /*if (ip==10 || ip==11 || ip==12)
      {
      yieldcalc=Total_Ybincounting;
      yielderror=Final_pro_error;
      }
    */
    hintegral_yield->SetBinContent(ip + 1, yieldcalc);
    hintegral_yield->SetBinError(ip + 1, yielderror); // filling histogram including error
    

    cout << "Total Yield from function integration:" << integralsignalfunc[ip] << " error : " << interror[ip] << " " << ip << endl;
    cout << "Yield from function integration:" << yieldcalc << " error : " << interror[ip] << " " << ip << endl;
    cout << "--------Frac error from integral method----------" << (yielderror / yieldcalc) << endl;

    //*********************************************************************************************************************

    //**Calculation for raw pt spectra using fitting parameter and filling it in histogram*********************************************

    Yield_value_par = Yield[ip] / (Event * ptbinwidth[ip] * dy * BR * 2 * binwidth_file * ptcenter[ip]);
    Yield_error_par = ErrorYield[ip] / (Event * ptbinwidth[ip] * dy * BR * 2 * binwidth_file * ptcenter[ip]);

    hYieldpar->SetBinContent(ip + 1, Yield_value_par);
    hYieldpar->SetBinError(ip + 1, Yield_error_par);

    cout << "Total Yield from fitting:" << Yield[ip] << endl;
    cout << "Yield from fitting:" << Yield_value_par << endl;
    //***********************************************************************************************************************************

    //**Filling mass and width fitting parameter in histogram*******************************************************************************

    hmass->SetBinContent(ip + 1, Mass[ip]);
    hmass->SetBinError(ip + 1, ErrorMass[ip]);

    hwidth->SetBinContent(ip + 1, Width[ip]);
    hwidth->SetBinError(ip + 1, ErrorWidth[ip]);

    //*****************************************************************************************************************************

    //**Setting plot parameters style*************************************************************************************************

    // SetHistoStyle(hfsig,1,20,0.6,0.04,0.04,0.03,0.03,1.0,1.4);
    // SetHistoStyle(fHistTotal[ip],1,24,0.5,0.04,0.04,0.03,0.03,1.0,1.4);
    SetHistoStyle(hfsig, 1, 20, 1.0, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
    SetHistoStyle(fHistTotal[ip], 1, 24, 1.0, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

    hfsig->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
    hfsig->GetYaxis()->SetMaxDigits(2);
    hfsig->GetYaxis()->SetTitle(Form("Counts/(%0.2f GeV/#it{c}^{2})", binwidth_file));

    if (ResBkg.CompareTo("MIX") == 0)
    {

      SetHistoStyle(hfbkg, 2, 24, 0.5, 0.05, 0.05, 0.03, 0.03, 1.0, 1.4);
      SetHistoStyle(hfbkg, kRed, 24, 1.0, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

      hfbkg->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
      hfsig->GetYaxis()->SetMaxDigits(2);
      hfbkg->GetYaxis()->SetTitle(Form("Counts/(%0.2f GeV/#it{c}^{2})", binwidth_file));
    }

    else if (ResBkg.CompareTo("LIKE") == 0)
    {
      // SetHistoStyle(hfbkg_like,2,24,0.5,0.05,0.05,1.0,1.1);
      SetHistoStyle(hfbkg_like, 2, 24, 0.5, 0.05, 0.05, 0.03, 0.03, 1.0, 1.4);
      hfbkg_like->GetXaxis()->SetTitle("M_{K#pi} (GeV/#it{c}^{2})");
      hfsig->GetYaxis()->SetMaxDigits(2);
      hfbkg_like->GetYaxis()->SetTitle(Form("Counts/(%0.2f GeV/#it{c}^{2})", binwidth_file));
    }

    fitFcn1->SetLineColor(4);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);
    fitFcn2->SetLineColor(6);
    fitFcn2->SetLineStyle(2);
    fitFcn2->SetLineWidth(2);

    //*******************************************************************************************************************************

    //**Plot of histograms and graphs*********************************************************************************************

    cinv[ip]->cd();
    hfsig->Draw("e");
    fitFcn->Draw("same");
    fitFcn1->Draw("same");
    fitFcn2->Draw("same");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    t2->DrawLatex(0.65, 0.97,
                  Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
                       High_pt[ip]));
    hfsig->GetXaxis()->SetRangeUser(al, bh);
    // TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
    TLegend *leg = new TLegend(0.280554, 0.4312735, 0.492902, 0.3338954, NULL, "brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->SetTextSize(0.04);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->AddEntry(hfsig, "Data (stat. uncert)", "pe");
    leg->AddEntry(fitFcn, "Breit-Wigner Peak Fit", "l");
    // leg->AddEntry(fitFcn2,"BW","l");
    leg->AddEntry(fitFcn1, "Res.bkg", "l");
    // leg->Draw();

    TLegend *leg2ap = new TLegend(0.8, 0.1, 0.9, 0.2, NULL, "brNDC");
    leg2ap->SetBorderSize(0);
    leg2ap->SetTextFont(22);
    leg2ap->SetTextSize(0.04);
    leg2ap->SetLineColor(1);
    leg2ap->SetLineStyle(1);
    leg2ap->SetLineWidth(1);
    leg2ap->SetFillColor(0);
    leg2ap->AddEntry((TObject *)0, "ALICE", "");
    // leg2ap->Draw();

    TLegend *leg2 = new TLegend(0.6, 0.7, 0.8, 0.9, NULL, "brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.04);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
    // leg2->AddEntry((TObject*)0,"ALICE Preliminary","");
    leg2->AddEntry((TObject *)0, "pp #sqrt{#it{s}} = 5.02 TeV (0-1%) ", "");
    leg2->AddEntry((TObject *)0, "K*^{0} #rightarrow K#pi", "");
    sprintf(name, "%0.1f < #it{p}_{T}(GeV/#it{c}) < %0.1f", Low_pt[ip], High_pt[ip]);
    leg2->AddEntry((TObject *)0, name, "");
    leg2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
    // leg2->Draw();

    cinv[ip]->SaveAs(Form("../output/LHC22q/figures2/hfitsig_pt%d.png", ip + 1));

    // if (ip<12)
    if (ip < 25)
    {
      c2->cd(ip + 1);
      hfsig->Draw("e");

      hfsig->GetXaxis()->SetRangeUser(al, bh);
      TF1 *ftotal = (TF1 *)fitFcn->Clone();
      TF1 *fBW = (TF1 *)fitFcn2->Clone();
      TF1 *fpoly = (TF1 *)fitFcn1->Clone();
      ftotal->Draw("same");
      fBW->Draw("same");
      fpoly->Draw("same");
      sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
      TLatex *ltx = new TLatex(0.35, 0.93, name);
      ltx->SetNDC();
      ltx->SetTextFont(42);
      ltx->SetTextSize(0.07);
      ltx->Draw();

      TLegend *leg1 = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
      // TLegend *leg = new TLegend(0.280554,0.4312735,0.492902,0.3338954,NULL,"brNDC");
      leg1->SetBorderSize(0);
      leg1->SetTextFont(42);
      leg1->SetTextSize(0.05);
      leg1->SetLineColor(1);
      leg1->SetLineStyle(1);
      leg1->SetLineWidth(1);
      leg1->SetFillColor(0);
      leg1->AddEntry(fitFcn, "Breit-Wigner+Res.bkg", "l");
      leg1->AddEntry(fBW, "Breit-Wigner", "l");
      leg1->AddEntry(fpoly, "Res.bkg", "l");
      leg1->Draw();
    }

    else
    {
      // c22->cd(ip-11);
      c22->cd(ip - 14);
      hfsig->Draw("e");
      hfsig->GetXaxis()->SetRangeUser(al, bh);
      TF1 *ftotal = (TF1 *)fitFcn->Clone();
      TF1 *fBW = (TF1 *)fitFcn2->Clone();
      TF1 *fpoly = (TF1 *)fitFcn1->Clone();
      ftotal->Draw("same");
      fBW->Draw("same");
      fpoly->Draw("same");
      sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
      TLatex *ltx = new TLatex(0.35, 0.93, name);
      ltx->SetNDC();
      ltx->SetTextFont(42);
      ltx->SetTextSize(0.07);
      ltx->Draw();

      TLegend *leg3 = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
      leg3->SetBorderSize(0);
      leg3->SetTextFont(22);
      leg3->SetTextSize(0.05);
      leg3->SetLineColor(1);
      leg3->SetLineStyle(1);
      leg3->SetLineWidth(1);
      leg3->SetFillColor(0);
      leg3->AddEntry(fitFcn, "Breit-Wigner+Res.bkg", "l");
      leg3->AddEntry(fBW, "Breit-Wigner", "l");
      leg3->AddEntry(fpoly, "Res.bkg", "l");
      leg3->Draw();
    }

    ////////////////////////////////////////////////////////////////////////////

    // signal extraction before bkg subtraction/////////////////////////////////////////////////////////////////////////

    // if (ip<12)
    if (ip < 25)
    {

      // For mixed event bkg/////////////////////////////////
      if (ResBkg.CompareTo("MIX") == 0)
      {

        c1->cd(ip + 1);

        hfbkg->Draw("e");
        hfbkg->GetXaxis()->SetRangeUser(al, bh);
        fHistTotal[ip]->Draw("same");
        sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
        TLatex *ltx = new TLatex(0.35, 0.93, name);
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.07);
        // ltx->Draw();
        TLegend *leg4 = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
        leg4->SetBorderSize(0);
        leg4->SetTextFont(42);
        leg4->SetTextSize(0.03);
        leg4->SetLineColor(1);
        leg4->SetLineStyle(1);
        leg4->SetLineWidth(1);
        leg4->SetFillColor(0);
        leg4->AddEntry((TObject *)0, "ALICE", "");
        leg4->AddEntry((TObject *)0, "pp #sqrt{#it{s}} = 5.02 TeV ", "");
        leg4->AddEntry((TObject *)0, "K*^{0} , |#it{y}| < 0.5 ", "");
        sprintf(name, "%0.1f < #it{p}_{T}(GeV/#it{c}) < %0.1f", Low_pt[ip], High_pt[ip]);
        leg4->AddEntry((TObject *)0, name, "");
        leg4->AddEntry((TObject *)0, "0-1 %", "");
        // leg4->Draw();
        TLegend *legsb = new TLegend(0.720554, 0.6312735, 0.892902, 0.7238954, NULL, "brNDC");
        legsb->SetBorderSize(0);
        legsb->SetTextFont(42);
        legsb->SetTextSize(0.03);
        legsb->SetLineColor(1);
        legsb->SetLineStyle(1);
        legsb->SetLineWidth(1);
        legsb->SetFillColor(0);
        legsb->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
        legsb->AddEntry(hfbkg, "Normalised mixed bkg", "p");
        legsb->Draw();
        Double_t max = hfbkg->GetMaximum();
        Double_t min = fHistTotal[ip]->GetMaximum();
        /*
        if (max>min)
        hfbkg->GetYaxis()->SetRangeUser(0,1.2*max);
        else
        hfbkg->GetYaxis()->SetRangeUser(0,1.2*min);
        */
      }
      //////////////////////////////////////////////////////

      // For like sign bkg//////////////////////////////////
      else if (ResBkg.CompareTo("LIKE") == 0)
      {

        c1->cd(ip + 1);
        // hfbkg_like->GetXaxis()->SetRangeUser(al, bh);
        fHistTotal[ip]->Draw("e");
        hfbkg_like->Draw("same");
        sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
        TLatex *ltx = new TLatex(0.35, 0.93, name);
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.045);
        ltx->Draw();
        TLegend *leg5 = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
        leg5->SetBorderSize(0);
        leg5->SetTextFont(22);
        leg5->SetTextSize(0.05);
        leg5->SetLineColor(1);
        leg5->SetLineStyle(1);
        leg5->SetLineWidth(1);
        leg5->SetFillColor(0);
        leg5->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
        leg5->AddEntry(hfbkg_like, "LikeSignPairs", "p");
        leg5->Draw();
      }
      /////////////////////////////////////////////////////
    }

    else
    {

      // For mixed event bkg/////////////////////////////////
      if (ResBkg.CompareTo("MIX") == 0)
      {

        // c11->cd(ip-11);
        c11->cd(ip - 14);
        hfbkg->Draw("e");
        hfbkg->GetXaxis()->SetRangeUser(al, bh);
        fHistTotal[ip]->Draw("same");
        sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
        TLatex *ltx = new TLatex(0.35, 0.93, name);
        ltx->SetNDC();
        ltx->SetTextFont(42);
        ltx->SetTextSize(0.07);
        // ltx->Draw();
        TLegend *leg6 = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
        leg6->SetBorderSize(0);
        leg6->SetTextFont(42);
        leg6->SetTextSize(0.03);
        leg6->SetLineColor(1);
        leg6->SetLineStyle(1);
        leg6->SetLineWidth(1);
        leg6->SetFillColor(0);
        leg6->AddEntry((TObject *)0, "ALICE", "");
        leg6->AddEntry((TObject *)0, "pp #sqrt{#it{s}} = 5.02 TeV ", "");
        leg6->AddEntry((TObject *)0, "K*^{0} , |#it{y}| < 0.5", "");
        sprintf(name, "%0.1f < #it{p}_{T}(GeV/#it{c}) < %0.1f", Low_pt[ip], High_pt[ip]);
        leg6->AddEntry((TObject *)0, name, "");
        leg6->AddEntry((TObject *)0, "0-1 %", "");
        TLegend *legsb7 = new TLegend(0.720554, 0.6312735, 0.892902, 0.7238954, NULL, "brNDC");
        legsb7->SetBorderSize(0);
        legsb7->SetTextFont(42);
        legsb7->SetTextSize(0.03);
        legsb7->SetLineColor(1);
        legsb7->SetLineStyle(1);
        legsb7->SetLineWidth(1);
        legsb7->SetFillColor(0);
        legsb7->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
        legsb7->AddEntry(hfbkg, "Normalised mixed bkg", "p");
        legsb7->Draw();
        // leg6->Draw();

        Double_t max = hfbkg->GetMaximum();
        Double_t min = fHistTotal[ip]->GetMaximum();
        /*if (max>min)
        hfbkg->GetYaxis()->SetRangeUser(0,2*max);
        else
        hfbkg->GetYaxis()->SetRangeUser(0,2*min);
        */
      }

      //////////////////////////////////////////////////////

      // For like sign bkg//////////////////////////////////
      else if (ResBkg.CompareTo("LIKE") == 0)
      {

        // c11->cd(ip-11);
        c11->cd(ip - 14);
        // hfbkg_like->GetXaxis()->SetRangeUser(al, bh);
        fHistTotal[ip]->Draw("e");
        hfbkg_like->Draw("same");
        sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
        TLatex *ltx = new TLatex(0.35, 0.93, name);
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.05);
        ltx->Draw();
        TLegend *leg8 = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
        leg8->SetBorderSize(0);
        leg8->SetTextFont(22);
        leg8->SetTextSize(0.05);
        leg8->SetLineColor(1);
        leg8->SetLineStyle(1);
        leg8->SetLineWidth(1);
        leg8->SetFillColor(0);
        leg8->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
        leg8->AddEntry(hfbkg_like, "LikeSignPairs", "p");
        leg8->Draw();
      }

      /////////////////////////////////////////////////////
    }

    cSigbkg[ip]->cd();

    // For mixed event bkg///////////////////////////////////////////////////
    if (ResBkg.CompareTo("MIX") == 0)
    {
      hfbkg->GetXaxis()->SetRangeUser(al, bh);
      fHistTotal[ip]->Draw("e");
      fHistTotal[ip]->GetYaxis()->SetTitle("Counts");
      hfbkg->Draw("same");
      sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
      TLatex *ltx2 = new TLatex(0.35, 0.93, name);
      ltx2->SetNDC();
      ltx2->SetTextFont(42);
      ltx2->SetTextSize(0.05);
      // ltx2->Draw();

      TLegend *leg2P = new TLegend(0.50554, 0.512735, 0.55, 0.738954, NULL, "brNDC");

      leg2P->SetBorderSize(0);
      leg2P->SetTextFont(22);
      leg2P->SetTextSize(0.04);
      leg2P->SetLineColor(1);
      leg2P->SetLineStyle(1);
      leg2P->SetLineWidth(1);
      leg2P->SetFillColor(0);
      leg2P->AddEntry((TObject *)0, "ALICE", "");
      // leg2P->Draw();

      TLegend *leg29 = new TLegend(0.280554, 0.3312735, 0.492902, 0.5338954, NULL, "brNDC");

      leg29->SetBorderSize(0);
      leg29->SetTextFont(42);
      leg29->SetTextSize(0.04);
      leg29->SetLineColor(1);
      leg29->SetLineStyle(1);
      leg29->SetLineWidth(1);
      leg29->SetFillColor(0);
      // leg2->AddEntry((TObject*)0,"ALICE Preliminary","");
      leg29->AddEntry((TObject *)0, "pp #sqrt{#it{s}} = 5.02 TeV (0-1%)", "");
      leg29->AddEntry((TObject *)0, "K*^{0} #rightarrow K#pi", "");
      sprintf(name, "%0.1f < #it{p}_{T}(GeV/#it{c}) < %0.1f", Low_pt[ip], High_pt[ip]);
      leg29->AddEntry((TObject *)0, name, "");
      leg29->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
      // leg29->Draw();

      TLegend *leg112 = new TLegend(0.680554, 0.8112735, 0.892902, 0.9138954, NULL, "brNDC");
      leg112->SetBorderSize(0);
      leg112->SetTextFont(2);
      leg112->SetTextSize(0.035);
      leg112->SetLineColor(1);
      leg112->SetLineStyle(1);
      leg112->SetLineWidth(1);
      leg112->SetFillColor(0);
      leg112->SetFillStyle(0);
      leg112->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
      leg112->AddEntry(hfbkg, "Mixed-event bkg", "p");
      leg112->Draw();

      sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
      TLatex *ltx = new TLatex(0.35, 0.93, name);
      ltx->SetNDC();
      ltx->SetTextFont(22);
      ltx->SetTextSize(0.05);
      ltx->Draw();

      Double_t max = hfbkg->GetMaximum();
      Double_t min = fHistTotal[ip]->GetMaximum();
      /*if (max>min)
      hfbkg->GetYaxis()->SetRangeUser(0,max);
      else
      hfbkg->GetYaxis()->SetRangeUser(0,min);
      */
      cSigbkg[ip]->SaveAs(Form("../output/LHC22q/figures2/hsigbkg_pt%d.png", ip + 1));
    }

    ////////////////////////////////////////////////////////////////////////

    // For like sign bkg///////////////////////////////////////////////////
    else if (ResBkg.CompareTo("LIKE") == 0)
    {
      // hfbkg_like->GetXaxis()->SetRangeUser(al, bh);
      fHistTotal[ip]->Draw("e");
      fHistTotal[ip]->GetYaxis()->SetTitle("Counts");
      hfbkg_like->Draw("same");
      TLegend *leg11 = new TLegend(0.680554, 0.8112735, 0.892902, 0.9138954, NULL, "brNDC");
      leg11->SetBorderSize(0);
      leg11->SetTextFont(2);
      leg11->SetTextSize(0.04);
      leg11->SetLineColor(1);
      leg11->SetLineStyle(1);
      leg11->SetLineWidth(1);
      leg11->SetFillColor(0);
      leg11->SetFillStyle(0);
      leg11->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
      leg11->AddEntry(hfbkg_like, "LikeSignPairs", "p");
      leg11->Draw();
      sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
      TLatex *ltx = new TLatex(0.35, 0.93, name);
      ltx->SetNDC();
      ltx->SetTextFont(22);
      ltx->SetTextSize(0.05);
      ltx->Draw();

      cSigbkg[ip]->SaveAs(Form("../output/LHC22q/figures2/hsigbkg_pt%d.png", ip + 1));
    }

    ////////////////////////////////////////////////////////////////////////

    cfit[ip]->cd();

  } // pt loop ends

  TCanvas *cgr = new TCanvas("cgr", "graphs", 900, 800);
  cgr->Range(0, 0, 1, 1);
  cgr->SetBorderSize(2);
  cgr->SetBorderMode(0);
  cgr->SetFillColor(10);
  cgr->SetFrameFillColor(10);
  cgr->SetFrameLineWidth(2);
  cgr->SetLeftMargin(0.2);
  cgr->SetRightMargin(0.05);
  cgr->SetTopMargin(0.08);
  cgr->SetBottomMargin(0.2);
  cgr->SetTicks(1, 1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
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
  // gStyle->SetErrorX(0);

  // chisquare_NDF vs pt
  TGraph *grchi = new TGraph(Npt, x, Chi2Ndf);
  grchi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  // grchi->GetXaxis()->SetLabelFont(43);
  // grchi->GetXaxis()->SetLabelSize(21);
  // grchi->GetXaxis()->SetTitleFont(43);
  // grchi->GetXaxis()->SetTitleSize(28);
  grchi->GetYaxis()->SetTitle("#chi^{2}/NDF ");
  // grchi->GetYaxis()->SetLabelFont(43);
  // grchi->GetYaxis()->SetLabelSize(21);
  // grchi->GetYaxis()->SetTitleFont(43);
  // grchi->GetYaxis()->SetTitleSize(28);
  // grchi->SetMarkerStyle(20);
  // grchi->SetMarkerSize(1.5);
  grchi->SetMarkerStyle(20);
  grchi->SetMarkerColor(1);
  grchi->SetLineColor(1);
  grchi->SetMarkerSize(1.4);
  grchi->SetLineWidth(2);
  // grchi->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  grchi->GetXaxis()->CenterTitle(false);
  grchi->GetXaxis()->SetNdivisions(506);
  grchi->GetYaxis()->SetNdivisions(505);
  grchi->GetXaxis()->SetLabelOffset(0.015);
  grchi->GetXaxis()->SetLabelFont(42);
  grchi->GetXaxis()->SetTitleFont(42);
  grchi->GetXaxis()->SetLabelSize(0.04);
  grchi->GetXaxis()->SetTitleSize(0.04);
  grchi->GetXaxis()->SetTickLength(0.04);
  grchi->GetXaxis()->SetTitleOffset(1.2);
  grchi->GetYaxis()->SetTitleOffset(1.7);
  grchi->GetYaxis()->CenterTitle(true);
  grchi->GetYaxis()->SetDecimals(false);
  // grchi->GetYaxis()->SetNdivisions(310);
  grchi->GetYaxis()->SetLabelOffset(0.015);
  grchi->GetYaxis()->SetLabelFont(42);
  grchi->GetYaxis()->SetLabelSize(0.04);
  grchi->GetYaxis()->SetTickLength(0.04);
  grchi->GetYaxis()->SetTitleSize(0.04);
  grchi->GetYaxis()->SetTitleFont(42);
  grchi->Draw("AP");
  t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
  // cgr->SetGrid();
  cgr->SaveAs("../output/LHC22q/figures2/chindfvspt.png");

  // mass vs pt
  TGraphErrors *grmass = new TGraphErrors(Npt, x, Mass, 0, ErrorMass);
  grmass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  // grmass->GetXaxis()->SetLabelFont(43);
  // grmass->GetXaxis()->SetLabelSize(21);
  // grmass->GetXaxis()->SetTitleFont(43);
  // grmass->GetXaxis()->SetTitleSize(28);
  grmass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
  // grmass->GetYaxis()->SetLabelFont(43);
  // grmass->GetYaxis()->SetLabelSize(21);
  // grmass->GetYaxis()->SetTitleFont(43);
  // grmass->GetYaxis()->SetTitleSize(28);
  // grmass->SetMarkerStyle(20);
  // grmass->SetMarkerSize(1.5);
  grmass->SetMarkerStyle(20);
  grmass->SetMarkerColor(1);
  grmass->SetLineColor(1);
  grmass->SetMarkerSize(1.4);
  grmass->SetLineWidth(2);
  grmass->GetXaxis()->CenterTitle(false);
  grmass->GetXaxis()->SetNdivisions(506);
  grmass->GetYaxis()->SetNdivisions(505);
  grmass->GetXaxis()->SetLabelOffset(0.015);
  grmass->GetXaxis()->SetLabelFont(42);
  grmass->GetXaxis()->SetTitleFont(42);
  grmass->GetXaxis()->SetLabelSize(0.04);
  grmass->GetXaxis()->SetTitleSize(0.04);
  grmass->GetXaxis()->SetTickLength(0.04);
  grmass->GetXaxis()->SetTitleOffset(1.2);
  grmass->GetYaxis()->SetTitleOffset(1.9);
  grmass->GetYaxis()->CenterTitle(true);
  grmass->GetYaxis()->SetDecimals(false);
  grmass->GetYaxis()->SetLabelOffset(0.015);
  grmass->GetYaxis()->SetLabelFont(42);
  grmass->GetYaxis()->SetLabelSize(0.04);
  grmass->GetYaxis()->SetTickLength(0.04);
  grmass->GetYaxis()->SetTitleSize(0.04);
  grmass->GetYaxis()->SetTitleFont(42);
  grmass->Draw("ap");
  TLegend *massleg = new TLegend(0.65, 0.45, 0.9, 0.65);
  massleg->SetFillColor(0);
  t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
  TLine *line = new TLine(grmass->GetXaxis()->GetXmin(), 0.895, grmass->GetXaxis()->GetXmax(), 0.895);
  line->SetLineStyle(2);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->Draw();
  massleg->AddEntry(line, "PDG Mass", "l");
  massleg->Draw("l");
  cgr->SaveAs("../output/LHC22q/figures2/mass_pt.png");

  // Width vs pT
  TGraphErrors *grwidth = new TGraphErrors(Npt, x, Width, 0, ErrorWidth);
  grwidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  grwidth->GetYaxis()->SetTitle("Width (GeV)");
  grwidth->SetMarkerStyle(20);
  grwidth->SetMarkerColor(1);
  grwidth->SetLineColor(1);
  grwidth->SetMarkerSize(1.4);
  grwidth->SetLineWidth(2);
  // grwidth->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  grwidth->GetXaxis()->CenterTitle(false);
  grwidth->GetXaxis()->SetNdivisions(506);
  grwidth->GetYaxis()->SetNdivisions(505);
  grwidth->GetXaxis()->SetLabelOffset(0.015);
  grwidth->GetXaxis()->SetLabelFont(42);
  grwidth->GetXaxis()->SetTitleFont(42);
  grwidth->GetXaxis()->SetLabelSize(0.04);
  grwidth->GetXaxis()->SetTitleSize(0.04);
  grwidth->GetXaxis()->SetTickLength(0.04);
  grwidth->GetXaxis()->SetTitleOffset(1.2);
  grwidth->GetYaxis()->SetTitleOffset(1.7);
  grwidth->GetYaxis()->CenterTitle(true);
  grwidth->GetYaxis()->SetDecimals(false);
  // grwidth->GetYaxis()->SetNdivisions(310);
  grwidth->GetYaxis()->SetLabelOffset(0.015);
  grwidth->GetYaxis()->SetLabelFont(42);
  grwidth->GetYaxis()->SetLabelSize(0.04);
  grwidth->GetYaxis()->SetTickLength(0.04);
  grwidth->GetYaxis()->SetTitleSize(0.04);
  grwidth->GetYaxis()->SetTitleFont(42);
  grwidth->Draw("ap");
  TLegend *widthleg = new TLegend(0.65, 0.75, 0.9, 0.85);
  widthleg->SetFillColor(0);
  t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
  TLine *line2 = new TLine(grmass->GetXaxis()->GetXmin(), 0.047, grmass->GetXaxis()->GetXmax(), 0.047);
  line2->SetLineStyle(2);
  line2->SetLineColor(2);
  line2->SetLineWidth(3);
  line2->Draw();
  widthleg->AddEntry(line, "PDG Width", "l");
  widthleg->Draw("l");
  cgr->SaveAs("../output/LHC22q/figures2/width_pt.png");

  // yield vs pt

  double xerr[Npt] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5, 1, 1, 2, 2, 3};


  TCanvas *cyield1 = new TCanvas("cyield1", "yeild", 1000, 800);
    cyield1->Range(0, 0, 1, 1);
    cyield1->SetBorderSize(2);
    cyield1->SetBorderMode(0);
    cyield1->SetFillColor(10);
    cyield1->SetFrameFillColor(10);
    cyield1->SetFrameLineWidth(2);
    cyield1->SetLeftMargin(0.2);
    cyield1->SetRightMargin(0.05);
    cyield1->SetTopMargin(0.08);
    cyield1->SetBottomMargin(0.2);
    cyield1->SetTicks(1, 1);
    cyield1->SetLogy();

    TGraphErrors *gryield1 = new TGraphErrors(Npt, x, yieldcalc1, xerr, yielderror1 );
    TLegend *legyield = new TLegend(0.5,0.7,0.9,0.9);

    gryield1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gryield1->GetYaxis()->SetTitle("Yield (GeV/c)^{-1}");
    gryield1->SetMarkerStyle(20);
    gryield1->SetMarkerColor(1);
    gryield1->SetLineColor(1);
    gryield1->SetMarkerSize(1.4);
    gryield1->SetLineWidth(2);
    gryield1->GetXaxis()->CenterTitle(false);
    gryield1->GetXaxis()->SetNdivisions(506);
    gryield1->GetYaxis()->SetNdivisions(505);
    gryield1->GetXaxis()->SetLabelOffset(0.015);
    gryield1->GetXaxis()->SetLabelFont(42);
    gryield1->GetXaxis()->SetTitleFont(42);
    gryield1->GetXaxis()->SetLabelSize(0.04);
    gryield1->GetXaxis()->SetTitleSize(0.04);
    gryield1->GetXaxis()->SetTickLength(0.04);
    gryield1->GetXaxis()->SetTitleOffset(1.2);
    gryield1->GetYaxis()->SetTitleOffset(1.7);
    gryield1->GetYaxis()->CenterTitle(true);
    gryield1->GetYaxis()->SetDecimals(false);
    gryield1->GetYaxis()->SetLabelOffset(0.015);
    gryield1->GetYaxis()->SetLabelFont(42);
    gryield1->GetYaxis()->SetLabelSize(0.04);
    gryield1->GetYaxis()->SetTickLength(0.04);
    gryield1->GetYaxis()->SetTitleSize(0.04);
    gryield1->GetYaxis()->SetTitleFont(42);

    // gryield1->Draw("ap");
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");

    grsourav->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    grsourav->GetYaxis()->SetTitle("(1/N_{ev})* d^{2}N/(dp_{T} dy) (Gev/c)^{-1}");
    grsourav->GetYaxis()->SetRangeUser(10e-8, 1);
    grsourav->SetMarkerStyle(20);
    grsourav->SetMarkerColor(2);
    grsourav->SetLineColor(2);
    grsourav->SetMarkerSize(1.4);
    grsourav->SetLineWidth(2);
    grsourav->GetXaxis()->CenterTitle(false);
    grsourav->GetXaxis()->SetNdivisions(506);
    grsourav->GetYaxis()->SetNdivisions(505);
    grsourav->GetXaxis()->SetLabelOffset(0.015);
    grsourav->GetXaxis()->SetLabelFont(42);
    grsourav->GetXaxis()->SetTitleFont(42);
    grsourav->GetXaxis()->SetLabelSize(0.04);
    grsourav->GetXaxis()->SetTitleSize(0.04);
    grsourav->GetXaxis()->SetTickLength(0.04);
    grsourav->GetXaxis()->SetTitleOffset(1.2);
    grsourav->GetYaxis()->SetTitleOffset(1.7);
    grsourav->GetYaxis()->CenterTitle(true);
    grsourav->GetYaxis()->SetDecimals(false);
    grsourav->GetYaxis()->SetLabelOffset(0.015);
    grsourav->GetYaxis()->SetLabelFont(42);
    grsourav->GetYaxis()->SetLabelSize(0.04);
    grsourav->GetYaxis()->SetTickLength(0.04);
    grsourav->GetYaxis()->SetTitleSize(0.04);
    grsourav->GetYaxis()->SetTitleFont(42);
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    legyield->AddEntry(grsourav, "pp 13TeV (Published)");
    legyield->AddEntry(gryield1, "pp 13.6TeV (This Analysis)");
    
    grsourav->Draw("ap");
    gryield1->Draw("psame");
    legyield->SetTextSize(0.03); // Adjust the value to your desired size
    legyield->SetTextFont(2); // Set the legend text font to bold    
    legyield->Draw();
    cyield1->SaveAs("../output/LHC22q/figures2/yield1.png");













}
