Double_t rBreitWigner(Double_t *x, Double_t *par)
{

  /// double rBW1 = par[0]*par[2]/(2*3.14159)/((x[0]-par[1])**2+par[2]**2/4.);
  // double rBW2 = par[3]*par[5]/(2*3.14159)/((x[0]-par[4])**2+par[5]**2/4.);
  //  double rBW3 = par[6]*par[8]/(2*3.14159)/((x[0]-par[7])**2+par[9]**2/4.);
  //  double rBW4 = par[9]*par[11]/(2*3.14159)/((x[0]-par[10])**2+par[11]**2/4.);

  double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
  double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
  double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
  double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
  double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

  double Gamma = par[2];
  double Gamma1 = par[5];
  double Gamma2 = par[8];
  double Gamma3 = par[11];

  double rBW1 = par[0] * x[0] * par[1] * Gamma / (TMath::Power((x[0] * x[0] - par[1] * par[1]), 2.0) + par[1] * par[1] * Gamma * Gamma);
  double rBW2 = par[3] * x[0] * par[4] * Gamma1 / (TMath::Power((x[0] * x[0] - par[4] * par[4]), 2.0) + par[4] * par[4] * Gamma1 * Gamma1);
  double rBW3 = par[6] * x[0] * par[7] * Gamma2 / (TMath::Power((x[0] * x[0] - par[7] * par[7]), 2.0) + par[7] * par[7] * Gamma2 * Gamma2);
  double rBW4 = par[9] * x[0] * par[10] * Gamma2 / (TMath::Power((x[0] * x[0] - par[10] * par[10]), 2.0) + par[10] * par[10] * Gamma2 * Gamma2);

  double sigfun = par[12] * (TMath::Power(((5 * rBW1) - (3 * rBW2) + 2 * rBW3), 2)) + par[13] * (TMath::Power(rBW4, 2.0));
  double poly2 = par[14] * TMath::Power((x[0] - 2.0 * 0.4976), par[15]) * TMath::Exp(-(x[0] - 2.0 * 0.4976) * par[16]);

  return (sigfun + poly2);
}

Double_t rBW(Double_t *x, Double_t *par)
{
  //  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  // double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);

  // double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  // double rBW = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  //  double BF = TMath::Sqrt((x[0]*x[0]+par[3]*par[3]));
  //  double PS =  x[0]*TMath::Exp(-BF/par[4])/(BF);

  // double rBW = ;

  double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
  double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
  double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
  double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
  double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

  //  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  // double Gamma1 = par[5]*(TMath::Power(par[4]/x[0],1.0))*TMath::Power((npart1)/(dpart2),1.5);
  // double Gamma2 = par[8]*(TMath::Power(par[7]/x[0],1.0))*TMath::Power((npart1)/(dpart3),1.5);
  // double Gamma3 = par[11]*(TMath::Power(par[10]/x[0],1.0))*TMath::Power((npart1)/(dpart4),1.5);

  // double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  //   double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  //  double rBW3 = par[6]*x[0]*par[7]*Gamma2/(TMath::Power((x[0]*x[0]-par[7]*par[7]),2.0)+par[7]*par[7]*Gamma2*Gamma2);
  //  double rBW4 = par[9]*x[0]*par[10]*Gamma3/(TMath::Power((x[0]*x[0]-par[10]*par[10]),2.0)+par[10]*par[10]*Gamma3*Gamma3);
  //  double rBW1 = par[0]*par[1]*par[1]*Gamma/(TMath::Power(( -x[0]*x[0]+par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double Gamma = par[2];
  double Gamma1 = par[5];
  double Gamma2 = par[8];
  double Gamma3 = par[11];

  /*
  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);
  double  dpart2 = par[4]*par[4] - 4*(0.4976*0.4976);

  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);

  double Gamma1 = par[5]*(TMath::Power(par[4]/x[0],1.0))*TMath::Power((npart1)/(dpart2),1.5);

  double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  */

  /* double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);

  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);

  double Gamma1 = par[5]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);

  double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double rBW2 = par[3]*x[0]*par[4]*Gamma/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  */

  //  return(rBW1+rBW2+rBW3+rBW4);
  //  return(rBW1);

  double rBW1 = par[0] * x[0] * par[1] * Gamma / (TMath::Power((x[0] * x[0] - par[1] * par[1]), 2.0) + par[1] * par[1] * Gamma * Gamma);
  double rBW2 = par[3] * x[0] * par[4] * Gamma1 / (TMath::Power((x[0] * x[0] - par[4] * par[4]), 2.0) + par[4] * par[4] * Gamma1 * Gamma1);
  double rBW3 = par[6] * x[0] * par[7] * Gamma2 / (TMath::Power((x[0] * x[0] - par[7] * par[7]), 2.0) + par[7] * par[7] * Gamma2 * Gamma2);
  double rBW4 = par[9] * x[0] * par[10] * Gamma2 / (TMath::Power((x[0] * x[0] - par[10] * par[10]), 2.0) + par[10] * par[10] * Gamma2 * Gamma2);

  double sigfun = par[12] * (TMath::Power((5 * rBW1 - 3 * rBW2 + 2 * rBW3), 2)) + par[13] * (TMath::Power(rBW4, 2.0));

  //      return(rBW1+rBW2+rBW3+rBW4);
  return (sigfun);
}

Double_t BreitWigner(Double_t *x, Double_t *par)
{
  return (par[0] * par[2] / (2 * 3.14159) / ((x[0] - par[1]) * *2 + par[2] * *2 / 4.) + par[3] + x[0] * par[4] + x[0] * x[0] * par[5]);
}

Double_t polynomial2(Double_t *x, Double_t *par)
{
  // double poly2 = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];

  //  double poly2 = par[0]*TMath::Power(x[0],par[1])*TMath::Exp(-x[0]*par[2]);

  //  double poly2 = par[3]*TMath::Power(x[0],par[4])*TMath::Exp(-x[0]*par[5]);

  //  double poly2 = par[0]*TMath::Sqrt(TMath::Power(x[0] -(2*0.497),par[1]))*TMath::Power(par[2],1.5)*TMath::Exp(-par[2]*(TMath::Power(x[0]-(2*0.497)),par[1]));
  // double poly2 = par[3]*exp(x[0]*par[4]);

  // double poly2 = par[0]*exp(x[0]*par[1])+par[2]*exp(x[0]*par[3]*x[0])+par[4]*exp(x[0]*x[0]*x[0]*par[5]);

  double poly2 = par[0] * TMath::Power((x[0] - 2.0 * 0.4976), par[1]) * TMath::Exp(-(x[0] - 2.0 * 0.4976) * par[2]);
  return (poly2);
  // return (bg);

  //  return ((x[0]- 0.63718)**par[3])*exp(par[0] + x[0]*par[1] + x[0]*x[0]*par[2]);
}

Double_t BW(Double_t *x, Double_t *par)
{
  return par[0] * par[2] / (2 * 3.14159) / ((x[0] - par[1]) * *2 + par[2] * *2 / 4.);
  // return (0.5*par[0]*par[1]/TMath::Pi() /((x[0]-par[2])*(x[0]-par[2]) + 0.25*par[1]*par[1]));
}

Double_t dy = 1.0;
// const Int_t Npt=7;
const Int_t Npt = 1;
const Int_t N = 1;

// Double_t lownorm = 1.1;
//  Double_t lownorm = 1.2;
// Double_t hinorm = 1.15;

Double_t lownorm = 2.2;
// Double_t lownorm = 1.2;
Double_t hinorm = 2.3;

// Double_t lownorm1= 1.1;
//  Double_t lownorm = 1.2;
// Double_t hinorm1= 1.1;

Double_t lownorm1 = 2.2;
// Double_t lownorm =1.2;
Double_t hinorm1 = 2.3;

// Double_t hinorm = 1.3;
// Double_t hinorm = 1.3; // n0
// Double_t hinorm = 1.4; // n1
// Double_t hinorm = 1.3; // n2
Double_t lowsigmarange = 1.6;
Double_t highsigmarange = 2.0;
void Signal_4rBW_1to10_Default_free()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
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

  Int_t STCol = 1;
  Int_t STSty = 21;

  Int_t CBCol = 1;
  Int_t CBSty = 21;

  TString ptbn = "1 #leq #it{p}_{T} < 10 GeV/#it{c}";
  TLegend *l2 = DrawLegend(0.25, 0.05, 0.4, 0.15);

  Double_t fRawYield[Npt];
  Double_t fErRawYield[Npt];
  Double_t wRawYield[Npt];
  Double_t wErRawYield[Npt];

  Double_t d, c, dpt, mass, error_mass, width, error_width;
  Double_t chi;
  TFitResultPtr r;
  TH1D *fHistTotal[Npt];
  TH1D *fHistBkg[Npt];
  TH1D *fHistlike[Npt];
  TH1D *histPP[Npt];
  TH1D *histMM[Npt];
  char name[100];

  Double_t yield, Yield, mass[Npt], fyield[Npt], Width[Npt], yerror[Npt], xerror[Npt], xerr[Npt], yerr[Npt], dn[Npt];
  Double_t poly0[Npt], poly1[Npt], poly2[Npt], poly3[Npt];
  Double_t Par[Npt], pt, hBCError_1;

  Double_t yield_error_par, bmin, bmax, Yield_bincount_hist, Integral_BW_withsigma, fYield_BinCount, YieldIntegral_BW, Yfraction_cBW, Yield_value_par, sum_tail_correction, sum_tail_correction, Total_Ybincounting, Tail_correction_plusm, Tail_correction_minusm, Error_2, Final_pro_error, bkgvalue, Yield_value_par1710, yield_error_par1710;

  Double_t lc, hc, low_value, high_value, lbin, hbin, sig_integral, sig_integral1, bkg_integral1, bkg_integral, mT, mT_error, significance_num, significance_den, ratio, BinWidth;

  //   Double_t Low_pt[Npt+1] ={0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,10.0,12.0,15.0,20.0};

  //   Double_t Low_pt[Npt+1] ={0.8,1.2,1.6,1.8,2.0,2.4,2.8,3.2,4.0,5.0};

  //   Double_t Low_pt[Npt+1] ={0.8,1.2,1.6,1.8,2.0,2.4,2.8,3.2,4.0,5.0};

  //   Double_t Low_pt[Npt+1] ={0.5,0.7,1.2,1.6,2.0,2.5,3.0,3.5,4.0,5.0,7.0,10.0,15.0,20.0};

  //   Double_t Low_pt[Npt+1] ={0.8,1.2,1.6,2.0,2.5,3.0,4.0,5.0};

  //   Double_t Low_pt[Npt+1] ={1.0,2.0,3.0,4.0,6.0,8.0,10.0,15.0,20.0};

  //   Double_t Low_pt[Npt+1] ={0.0,1.0,1.2,1.6,1.8,2.0,2.5,3.0,4.0};

  // Double_t Low_pt[Npt+1] ={1.0,2.0,3.0,4.0,5.0,6.0,8.0,12.0};

  // Double_t Low_pt[Npt+1] ={1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0};

  Double_t Low_pt[Npt + 1] = {1.0, 10.0};

  // Double_t Low_pt[Npt+1] ={0.0,10.0};
  Double_t Low_Centr[N] = {1};
  Double_t High_Centr[N] = {100};

  // 0-100 %
  //  Double_t fitlow[Npt] = {.74,.72,.72,.66,.66,.66,.66,.66,.66};
  //  Double_t fithigh[Npt]= {1.07,1.1,1.1,1.1,1.09,1.09,1.09,1.09,1.1};
  // 10-30 %
  // Double_t fitlow[Npt] = {.74,.72,.72,.66,.66,.7,.66,.66,.66};
  // Double_t fithigh[Npt]= {1.07,1.1,1.1,1.1,1.09,1.09,1.09,1.09,1.1};

  // 1-10 %
  // Double_t fitlow[Npt] = {.74,.72,.72,.66,.66,.7,.66,.66,.66};
  // Double_t fithigh[Npt]= {1.07,1.1,1.1,1.1,1.09,1.09,1.09,1.09,1.1};

  // Double_t fitlow[Npt] = {.75,.72,.75,.66,.66,.72,.66,.66,.66};
  // Double_t fithigh[Npt]= {1.09,1.1,1.1,1.12,1.09,1.09,1.09,1.09,1.1};

  // Double_t fitlow[Npt] = {.75,.72,.75,.7,.66,.72,.66,.66,.66};
  // Double_t fithigh[Npt]= {1.09,1.1,1.1,1.09,1.05,1.09,1.09,1.09,1.1};

  // 0 to 30
  //  Double_t fitlow[Npt] = {.75,.72,.75,.67,.67,.72,.66,.66,.66};
  //  Double_t fithigh[Npt]= {1.09,1.1,1.1,1.09,1.05,1.09,1.09,1.09,1.1};

  // 70 to 90
  // Double_t fitlow[Npt] = {.74,.72,.71,.67,.67,.72,.66,.66,.66};

  // 0 to 90 kstar binning
  // Double_t fitlow[Npt] = {.67,.67,.68,.66,.67,.72,.66,.66,.66};
  // Double_t fithigh[Npt]= {1.1,1.1,1.1,1.09,1.05,1.09,1.09,1.09,1.1};

  /*   Double_t fitlow[Npt] = {.71,.67,.68,.66,.71,.71,.72};
       Double_t fithigh[Npt]= {1.1,1.14,1.09,1.05,1.055,1.17,1.07};
  */

  // 20 -40 %
  // Double_t fitlow[Npt] = {.71,.7,.66,.69,.66,.66,.66};
  // Double_t fithigh[Npt]= {1.1,1.1,1.09,1.11,1.1,1.1,1.13};
  // 40 -60 %
  //   Double_t fitlow[Npt] = {.67,.66,.67,.69,.66,.66,.69};
  // Double_t fithigh[Npt]= {1.1,1.19,1.13,1.11,1.1,1.1,1.1};

  //   Double_t fitlow[Npt] = {.72,.66,.67,.69,.66,.66,.69};
  // Double_t fithigh[Npt]= {1.06,1.19,1.13,1.11,1.1,1.1,1.1};

  // Double_t fitlow[Npt] = {.72,.72,.72,.72,.72,.72,.72};
  // Double_t fithigh[Npt]= {1.06,1.06,1.06,1.06,1.06,1.06,1.06};

  //   Double_t fitlow[Npt] = {1.4,1.42,1.42,1.42,1.42,1.42,1.42,1.42};
  // Double_t fithigh[Npt]= {1.6,1.7,1.7,1.7,1.7,1.7,1.7,1.7};

  //   Double_t fitlow[Npt] = {1.48,1.45,1.4,1.37,1.4,1.45,1.45};
  // Double_t fithigh[Npt]= {1.6,1.65,1.65,1.7,1.7,1.75,1.75};

  //   Double_t fitlow[Npt] = {1.46,1.45,1.4,1.42,1.38,1.4,1.42};
  // Double_t fithigh[Npt]= {1.6,1.68,1.68,1.7,1.7,1.72,1.72};

  //   Double_t fitlow[Npt] = {1.46,1.45,1.4,1.42,1.38,1.4,1.42};
  // Double_t fithigh[Npt]= {1.6,1.68,1.68,1.7,1.7,1.72,1.72};

  //   Double_t fitlow[Npt] = {1.45,1.45,1.4,1.42,1.42,1.46,1.42};
  //   Double_t fitlow[Npt] = {1.45,1.42,1.4,1.42,1.42,1.46,1.42};
  // Double_t fithigh[Npt]= {1.6,1.71,1.68,1.7,1.7,1.74,1.72};

  //  Double_t fitlow[Npt] = {1.45,1.45,1.4,1.42,1.42,1.46,1.42};
  // Double_t fithigh[Npt]= {1.61,1.71,1.68,1.7,1.7,1.74,1.72};

  //   Double_t fitlow[Npt] = {1.45,1.45,1.4,1.42,1.42,1.46,1.42};
  // Double_t fithigh[Npt]= {1.61,1.71,1.68,1.7,1.7,1.74,1.72};

  // Double_t fitlow[Npt] = {1.38};
  // Double_t fithigh[Npt]= {1.86};

  //   Double_t fitlow[Npt] = {1.211};

  //   Double_t fitlow[Npt] = {1.15};
  // Double_t fithigh[Npt]= {2.0};

  //   Double_t fitlow[Npt] = {1.15, 1.15,1.15,1.15,1.15,1.15,1.15};
  // Double_t fithigh[Npt]= {2.0,2.0,2.0,2.0,2.0,2.0,2.0};

  //   Double_t fitlow[Npt] = {1.24, 1.24,1.24,1.4,1.4,1.4,1.4};
  // Double_t fithigh[Npt]= {2.0,1.95,1.9,1.7,1.7,1.7,1.7};

  //   Double_t fitlow[Npt] = {1.24,1.24,1.24,1.26,1.3,1.26,1.25};
  // Double_t fithigh[Npt]= {2.0,2.0,2.0,1.95,1.98,1.95,1.95};

  //   Double_t fitlow[Npt] = {1.15,1.15,1.15,1.16};
  // Double_t fithigh[Npt]= {2.0,2.0,2.0,1.95};

  //   Double_t fitlow[Npt] = {1.15,1.3,1.22,1.28};
  // Double_t fithigh[Npt]= {1.9,1.87,2.0,1.95};

  //   Double_t fitlow[Npt] = {1.2,1.3,1.29,1.3,1.32};
  // Double_t fithigh[Npt]= {1.9,1.85,1.97,1.97,1.99};

  // fitting range: Default fiitting range
  //  Double_t fitlow[Npt] = {1.2,1.3,1.29,1.26,1.26};
  //  Double_t fithigh[Npt]= {1.95,1.87,1.97,1.95,1.95};

  //  Double_t fitlow[Npt] = {1.2,1.3,1.29,1.26,1.28};
  // Double_t fithigh[Npt]= {1.96,1.87,1.97,1.95,1.95};

  //  Double_t fitlow[Npt] = {1.2,1.3,1.3,1.26,1.28};
  // Double_t fithigh[Npt]= {1.96,1.87,1.95,1.95,1.95};

  // Double_t fitlow[Npt] = {1.25};
  //  Double_t fithigh[Npt]= {2.0};

  Double_t fitlow[Npt] = {1.05};
  Double_t fithigh[Npt] = {2.2};

  // Double_t fitlow[Npt] = {1.2,1.3,1.3,1.26,1.26};
  // Double_t fithigh[Npt]= {1.95,1.87,1.87,1.95,1.95};

  //  Double_t fitlow[Npt] = {1.2,1.3,1.3,1.26,1.26};
  // Double_t fithigh[Npt]= {1.96,1.87,1.88,1.95,1.95};

  // fitting range new
  // Double_t fitlow[Npt] = {1.2,1.3,1.25,1.26,1.26};
  // Double_t fithigh[Npt]= {1.95,1.85,1.9,1.95,1.95};
  // Double_t fitlow[Npt] = {1.2,1.3,1.25,1.26,1.26};
  // Double_t fithigh[Npt]= {1.95,1.85,1.85,1.95,1.95};
  // Double_t fitlow[Npt] = {1.2,1.3,1.2,1.23,1.23};
  // Double_t fithigh[Npt]= {1.95,1.85,1.8,1.9,1.9};
  // Double_t fitlow[Npt] = {1.2,1.3,1.26,1.26,1.26};
  // Double_t fithigh[Npt]= {1.9,1.85,1.97,1.97,1.99};
  // Double_t Para1[Npt]=  {100, 1.0e+02, 1.0e+02, 1.0e+02,100};

  Double_t Para1[Npt] = {100};

  // Double_t fitlow[Npt] = {1.42};
  // Double_t fithigh[Npt]= {1.7};

  // Double_t fitlow[Npt] = {.75};
  //  Double_t fithigh[Npt]= {1.09};

  /*
  Double_t Para1[Npt]=  {100, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02};

  Double_t RBPar1[Npt]= {10,1,1,1,1,1,1};

  Double_t RBPar2[Npt]= {10,10 ,10 ,10 ,10,10,10};

  Double_t RBPar3[Npt]= {4,4,4,4,4,4,4};

  Double_t RBPar4[Npt]= {10,10,10,10,10,10,10};
  */
  // Double_t Para1[Npt]=  {1000};
  // Double_t RBPar1[Npt]= {1000};
  // Double_t RBPar2[Npt]= {-100};
  // Double_t RBPar3[Npt]= {20};
  //  Double_t RBPar4[Npt]= {10};

  // Double_t RBPar5[Npt]= {20};

  // Double_t RBPar6[Npt]= {50};

  /*
  Double_t Para1[Npt]=  {100000,1000,1000,1000,1000,1000,1000};

  Double_t RBPar1[Npt]= {-10,1000,1000,1000,1000,1000,1000};

  Double_t RBPar2[Npt]= {-10,-100,-100,-100,-100,-100,-100};

  Double_t RBPar3[Npt]= {20,20,20,20,20,20,20};

  Double_t RBPar4[Npt]= {10,10,10,10,10,10,10};

  Double_t RBPar5[Npt]= {20,20,20,20,20,20,20};

  Double_t RBPar6[Npt]= {50,50,50,50,50,50,50};

  */

  /*

  Double_t Para1[Npt]=  {100, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02};

  Double_t RBPar1[Npt]= {10,1,1,1,1,1,1,1};

  Double_t RBPar2[Npt]= {10,10 ,10 ,10 ,10,10,10,10};

  Double_t RBPar3[Npt]= {4,4,4,4,4,4,4,4};

  Double_t RBPar4[Npt]= {10,10,10,10,10,10,10,10};
  */

  /*

  Double_t Para1[Npt]=  {100};

  Double_t RBPar1[Npt]= {10};

  Double_t RBPar2[Npt]= {10};

  Double_t RBPar3[Npt]= {4};

  Double_t RBPar4[Npt]= {10};
  */

  // fitFcn->SetParameters(4000000,0.8916,0.0508, 10, 4, 10, 1);

  // canvas
  TCanvas *cSigbkg[Npt];
  TCanvas *cinv[Npt];

  TCanvas *c1 = new TCanvas("c1", "", 10, 10, 600, 600);
  c1->Divide(4, 3);
  TCanvas *c11 = new TCanvas("c11", "", 10, 10, 600, 600);
  c11->Divide(4, 3);
  TCanvas *c2 = new TCanvas("c2", "", 10, 10, 600, 600);
  c2->Divide(4, 3);
  TCanvas *c22 = new TCanvas("c22", "", 10, 10, 600, 600);
  c22->Divide(4, 3);
  for (Int_t ip = 0; ip < Npt; ip++)

  // for(Int_t ip = 0; ip <2; ip++)
  {
    TString cName = TString::Format("cSigbkg_pt_%2.1f-%2.1f", Low_pt[ip], Low_pt[ip + 1]);
    cSigbkg[ip] = new TCanvas(Form("cSigbkg%d", ip), cName.Data(), 10, 10, 600, 600);
    cSigbkg[ip]->SetLeftMargin(0.2);
    cSigbkg[ip]->SetRightMargin(0.05);
    cSigbkg[ip]->SetBottomMargin(0.2);
  }

  for (Int_t ip = 0; ip < Npt; ip++)
  //  for(Int_t ip = 0; ip <2; ip++)
  {
    TString cName = TString::Format("cinv_pt_%2.1f-%2.1f", Low_pt[ip], Low_pt[ip + 1]);
    cinv[ip] = new TCanvas(Form("cinv%d", ip), cName.Data(), 10, 10, 600, 600);

    cinv[ip]->SetRightMargin(0.05);
    cinv[ip]->SetBottomMargin(0.2);
  }

  TH1F *hYbincount = new TH1F("hYbincount", "", Npt, Low_pt);
  TH1F *hYieldpar = new TH1F("hYieldpar", "", Npt, Low_pt);
  TH1F *hChiSquare = new TH1F("hChiSquare", "", Npt, Low_pt);
  TH1F *hsgnfcance = new TH1F("hsgnfcance", "", Npt, Low_pt);
  TH1F *hmass = new TH1F("hmass", "", Npt, Low_pt);
  TH1F *hwidth = new TH1F("hwidth", "", Npt, Low_pt);

  TH1F *hmass1525 = new TH1F("hmass1525", "hmass1525", Npt, Low_pt);

  TH1F *hwidth1525 = new TH1F("hwidth1525", "hwidth1525", Npt, Low_pt);

  TH1F *hmass1710 = new TH1F("hmass1710", "hmass1710", Npt, Low_pt);

  TH1F *hwidth1710 = new TH1F("hwidth1710", "hmass1710", Npt, Low_pt);

  TH1F *htransmass = new TH1F("htransmass", "", Npt, Low_pt);
  TH1F *hFrac_stat_error = new TH1F("hFrac_stat_error", "", Npt, Low_pt);

  TH1F *hYieldpar1525 = new TH1F("hYieldpar1525", "", Npt, Low_pt);

  TH1F *hYieldpar1710 = new TH1F("hYieldpar1710", "", Npt, Low_pt);

  TFile *fInputFile = new TFile("../../../AnalysisResultsKshortKshortpp13TeV_fulldataset_aod.root", "Read");

  TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(1);
  fInputList = (TList *)fInputFile->Get(key->GetName());
  // fInputList->ls();

  // return 0;
  for (Int_t j = 0; j < N; j++) // centrality loop
  {
    Double_t lc = Low_Centr[j];
    Double_t hc = High_Centr[j];
    TH1F *hEVent = (TH1F *)fInputList->FindObject("hAEventsVsMulti"); // Event calculation

    Double_t Event = hEVent->Integral(lc, hc);

    cout << Event << "\t" << "\t" << lc << "\t" << hc << endl;

    THnSparseF *fHistNum = (THnSparseF *)fInputList->FindObject("KStarPlusMinusppData_ChargeKstar_KStarPlusMinus");
    //	THnSparseF * fHistNum1 = (THnSparseF *) fInputList->FindObject("KStarPlusMinusppData_ChargeKstar_AKStarPlusMinus");
    //	fHistNum->Add(fHistNum1);

    THnSparseF *fHistDen = (THnSparseF *)fInputList->FindObject("KStarPlusMinusppData_ChargeKstar_KStarPlusMinusmix");
    //	THnSparseF * fHistDen1 = (THnSparseF *) fInputList->FindObject("KStarPlusMinusppData_ChargeKstar_AKStarPlusMinusmix");
    //	fHistDen->Add(fHistDen1);

    for (Int_t m = 0; m < Npt; m++)
    // for(Int_t m=0;m<2 ;m++)
    {
      sprintf(name, "fHistNum%d", m);
      fHistTotal[m] = new TH1D(name, "inv_mass", 200, 0.5, 2.5);
      sprintf(name, "fHistbkg%d", m);
      fHistBkg[m] = new TH1D(name, "inv_mass", 200, 0.5, 2.5);
      sprintf(name, "histPP_pt%d", m);
    }
    for (Int_t ip = 0; ip < Npt; ip++) // start pt bin loop
    // for(Int_t ip = 0; ip <3; ip++)//start pt bin loop
    //	for(Int_t ip =2; ip <3; ip++)//start pt bin loop
    {
      low_value = Low_pt[ip];
      high_value = Low_pt[ip + 1];
      ;

      sprintf(name, "fHistbkg%d", ip);
      sprintf(name, "fHistNum%d", ip);

      fHistNum->GetAxis(2)->SetRange(lc, hc);
      fHistNum->GetAxis(2)->SetBit(TAxis::kAxisRange);

      lbin = fHistNum->GetAxis(1)->FindBin(low_value + 1.e-7); // projection pt bin //histogram sig+bkg
      hbin = fHistNum->GetAxis(1)->FindBin(high_value - 1.e-7);

      fHistNum->GetAxis(1)->SetRange(lbin, hbin);
      fHistNum->GetAxis(1)->SetBit(TAxis::kAxisRange);
      fHistTotal[ip] = fHistNum->Projection(0, "E"); // projection x//histogram sig+bkg

      fHistDen->GetAxis(2)->SetRange(lc, hc);
      fHistDen->GetAxis(2)->SetBit(TAxis::kAxisRange);

      fHistDen->GetAxis(1)->SetRange(lbin, hbin);
      fHistDen->GetAxis(1)->SetBit(TAxis::kAxisRange);
      fHistBkg[ip] = fHistDen->Projection(0, "E"); // projectionx combine background histogram

      cout << fHistTotal[ip]->GetBinWidth(12) << endl;

      fHistTotal[ip]->Rebin(2);
      fHistBkg[ip]->Rebin(2);

      // return 0;
      fHistBkg[ip]->SetName(name);
      fHistTotal[ip]->SetName(name);

      Int_t lownormbin = fHistTotal[ip]->GetXaxis()->FindBin(lownorm); // change it to bin no.//findbin for bkg normalisation
      Int_t hinormbin = fHistTotal[ip]->GetXaxis()->FindBin(hinorm);

      Int_t lownormbin1 = fHistTotal[ip]->GetXaxis()->FindBin(lownorm1); // change it to bin no.//findbin for bkg normalisation
      Int_t hinormbin1 = fHistTotal[ip]->GetXaxis()->FindBin(hinorm1);

      sig_integral = fHistTotal[ip]->Integral(lownormbin, hinormbin); // using intrgral method find counts for normalisation
      bkg_integral = fHistBkg[ip]->Integral(lownormbin, hinormbin);

      sig_integral1 = fHistTotal[ip]->Integral(lownormbin1, hinormbin1); // using intrgral method find counts for normalisation
      bkg_integral1 = fHistBkg[ip]->Integral(lownormbin1, hinormbin1);

      // error in normalisation factor
      if (sig_integral <= 0 || bkg_integral <= 0)
        continue;

      Double_t normfactor = sig_integral / bkg_integral;
      Double_t normfactorerror = sqrt(normfactor / bkg_integral + normfactor * normfactor / bkg_integral);

      TH1D *hfsig = (TH1D *)fHistTotal[ip]->Clone();
      // TH1D *hfsigplusbkg=(TH1D*)fHistTotal[ip]->Clone();
      TH1D *hfbkg = (TH1D *)fHistBkg[ip]->Clone(); // mixed event

      hfbkg->Scale(sig_integral1 / bkg_integral1); // normalised bkg for mixed event

      /*   TFile *fout_bkg = new TFile("K0sK0s_signalplusbkg_0to100_1to10GeV.root","recreate");
      fout_bkg->cd();
      hfbkg->Write();
      hfsig->Write();
      */

      /*if(ip ==0)
        {
        hfbkg->Scale(sig_integral1/bkg_integral1);//normalised bkg for mixed event

        }

        else
              {

        hfbkg->Scale(sig_integral/bkg_integral);//normalised bkg for mixed event
        }
      */

      hfsig->Add(hfbkg, -1); // signal after subtraction of bkg

      TH1D *hfsig_new = (TH1D *)hfsig->Clone();

      TFile *fout_sig = new TFile("K0sK0s_signal_0to100_1to10GeV.root", "recreate");
      fout_sig->cd();
      hfsig->Write();

      pt = (Low_pt[ip] + Low_pt[ip + 1]) / 2;
      dpt = Low_pt[ip + 1] - Low_pt[ip];

      // cout<<pt[ip]<<endl;
      BinWidth = hfsig->GetBinWidth(15);
      // cout<< BinWidth <<endl;

      // fitting BW and poly2
      /*TF1 *fitFcn = new TF1("fitfunc",BreitWigner,fitlow[ip],fithigh[ip],6);
        TF1 *fitFcn1 = new TF1("fitfunc1",polynomial2,fitlow[ip],fithigh[ip],3);
        TF1 *fitFcn2 = new TF1("fitFcn2",BW,fitlow[ip],fithigh[ip],3);
      */

      TF1 *fitFcn = new TF1("fitfunc", rBreitWigner, fitlow[ip], fithigh[ip], 17);
      TF1 *fitFcn1 = new TF1("fitfunc1", polynomial2, fitlow[ip], fithigh[ip], 3);
      TF1 *fitFcn2 = new TF1("fitFcn2", rBW, fitlow[ip], fithigh[ip], 14);

      //    TF1 *fitFcn = new TF1("fitfunc",rBreitWigner,fitlow[ip],fithigh[ip],9);
      // TF1 *fitFcn1 = new TF1("fitfunc1",polynomial2,fitlow[ip],fithigh[ip],3);
      // TF1 *fitFcn2 = new TF1("fitFcn2",rBW,fitlow[ip],fithigh[ip],6);

      //	     fitFcn->SetParameter(0,10000); //yield
      // fitFcn->SetParLimits(1,1.2,1.4); //inv. mass peak range
      // fitFcn->SetParamet(1,1.270); //inv. mass peak range
      // fitFcn->SetParLimits(2,0.03,0.1); //width-fixed
      //	     fitFcn->FixParameter(2,0.10); //width-fixed

      // fitFcn->SetParameter(0,10000); //yield
      // fitFcn->FixParameter(1,1.270); //inv. mass peak range
      //	     fitFcn->SetParLimits(2,0.03,0.15); //width-fixed
      // fitFcn->FixParameter(2,0.10); //width-fixed

      fitFcn->SetParameter(0, 1000); // yield
      fitFcn->SetParLimits(1, 1.01, 1.3);
      fitFcn->SetParLimits(2, 0.1, 0.19);

      //  fitFcn->SetParLimits(1,1.15,1.3); //inv. mass peak range
      // fitFcn->SetParLimits(2,0.03,0.1); //width-fixed
      //	    fitFcn->SetParameter(1,1.270); //inv. mass peak range
      //  fitFcn->SetParameter(2,0.186); //width-fixed

      fitFcn->SetParameter(3, 1001);      // yield
      fitFcn->SetParLimits(4, 1.15, 1.4); // inv. mass peak range
      fitFcn->SetParLimits(5, 0.07, 0.1); // width-fixed

      /// fitFcn->SetParameter(4,1.320); //inv. mass peak range
      // fitFcn->SetParameter(5,0.11); //width-fixed

      fitFcn->SetParameter(6, 10000); // yield

      fitFcn->SetParLimits(7, 1.46, 1.56); // inv. mass peak range

      fitFcn->SetParLimits(8, 0.03, 0.26); // width-fixed

      //  fitFcn->SetParameter(7,1.525); //inv. mass peak range
      // fitFcn->SetParameter(8,0.073); //width-fixed

      fitFcn->SetParameter(9, 10000); // yield

      // fitFcn->SetParameter(10,1.710); //inv. mass peak range
      fitFcn->SetParLimits(10, 1.56, 1.76); // inv. mass peak range

      //  fitFcn->SetParameter(11,0.086); //width-fixed
      fitFcn->SetParLimits(11, 0.06, 0.16); // width-fixed

      /*
      fitFcn->SetParameter(12,0.00001); //width-fixed

      fitFcn->SetParameter(13,0.0001); //width-fixed

      fitFcn->SetParameter(14,60000); //width-fixed
      fitFcn->SetParameter(15,-0.56933); //width-fixed
      fitFcn->SetParameter(16,3.7); //width-fixed
      */

      fitFcn->SetParameter(12, 0.00001); // width-fixed

      fitFcn->SetParameter(13, 0.0001); // width-fixed

      fitFcn->SetParameter(14, 60000);  // width-fixed
      fitFcn->SetParameter(15, -0.126); // width-fixed
      fitFcn->SetParameter(16, 3.7);    // width-fixed

      r = hfsig->Fit(fitFcn, "REBMS+"); // signal after like subtraction likesign  bkg
      Double_t *par = fitFcn->GetParameters();
      Yield = fitFcn->GetParameter(0);
      // poly0[ip]=fitFcn->GetParameter(6);

      // poly1[ip]=fitFcn->GetParameter(12);
      // poly2[ip]=fitFcn->GetParameter(13);
      // poly3[ip]=fitFcn->GetParameter(14)
      // cout << "BKGFunc" <<  poly1[ip] << "\t" << poly2[ip] << "\t" << poly3[ip] << endl;
      // fitFcn1->SetParameter(12,poly1[ip]);
      // fitFcn1->SetParameter(13,poly2[ip]);
      // fitFcn1->SetParameter(14,poly3[ip]);

      fitFcn1->SetParameters(&par[14]);
      fitFcn2->SetParameters(&par[0]);

      /*************GetChiSquare********************/
      d = fitFcn->GetChisquare();
      c = fitFcn->GetNDF();
      chi = d / c;
      Double_t Yield1320, Yield1320error, Yield1525, Yield1525error, Yield1710, Yield1710error, mass1320, error_mass1320, width1320, error_width1320, mass1525, mass1710, error_mass1525, error_mass1710, width1525, width1710, error_width1525, error_width1710;

      /**************mass spectra************/
      Yield1320 = fitFcn->GetParameter(3);
      Yield1320error = fitFcn->GetParError(3);

      mass1320 = fitFcn->GetParameter(4);
      error_mass1320 = fitFcn->GetParError(4);

      width1320 = fitFcn->GetParameter(5);
      error_width1320 = fitFcn->GetParError(5);

      Yield1525 = fitFcn->GetParameter(6);
      Yield1525error = fitFcn->GetParError(6);

      mass1525 = fitFcn->GetParameter(7);
      error_mass1525 = fitFcn->GetParError(7);

      width1525 = fitFcn->GetParameter(8);
      error_width1525 = fitFcn->GetParError(8);

      cout << "mass upper " << mass1525 + width1525 << "\t--mass lower----" << mass1525 - width1525 << endl;

      Yield1710 = fitFcn->GetParameter(9);
      Yield1710error = fitFcn->GetParError(9);

      mass1710 = fitFcn->GetParameter(10);
      error_mass1710 = fitFcn->GetParError(10);

      width1710 = fitFcn->GetParameter(11);
      error_width1710 = fitFcn->GetParError(11);

      hmass1525->SetBinContent(ip + 1, mass1525);
      hmass1525->SetBinError(ip + 1, error_mass1525);

      hwidth1525->SetBinContent(ip + 1, width1525);
      hwidth1525->SetBinError(ip + 1, error_width1525);

      hmass1710->SetBinContent(ip + 1, mass1710);
      hmass1710->SetBinError(ip + 1, error_mass1710);

      hwidth1710->SetBinContent(ip + 1, width1710);
      hwidth1710->SetBinError(ip + 1, error_width1710);

      // hmass1525->Draw();
      // return 0;
      cout << "Mass1525" << hmass1525->GetBinContent(1) << endl;

      /***********transversemass************/
      mT = TMath::Sqrt(pt * pt + mass * mass);
      mT_error = TMath::Sqrt(pt * pt * dpt * dpt + mass * mass * error_mass * error_mass) / mT;
      htransmass->SetBinContent(ip + 1, mT);
      htransmass->SetBinError(ip + 1, mT_error);

      /*************ERROR BIN COUNTING METHOD CALCULATION************/
      TF1 *fitFcn2_plusm = new TF1("fitFcn2_plusm", BW, fitlow[ip], fithigh[ip], 3);
      TF1 *fitFcn2_minusm = new TF1("fitFcn2_plusm", BW, fitlow[ip], fithigh[ip], 3);

      fitFcn2_plusm->FixParameter(0, Yield);
      fitFcn2_plusm->FixParameter(1, width + error_width);
      fitFcn2_plusm->FixParameter(2, mass + error_mass);
      fitFcn2_minusm->FixParameter(0, Yield);
      fitFcn2_minusm->FixParameter(1, width - error_width);
      fitFcn2_minusm->FixParameter(2, mass - error_mass);

      /********different method for Yield coounting check******************/

      bmin = hfsig->GetXaxis()->FindBin(lowsigmarange);
      bmax = hfsig->GetXaxis()->FindBin(highsigmarange);
      Yield_bincount_hist = hfsig->IntegralAndError(bmin, bmax, hBCError_1);
      bkgvalue = fitFcn1->Integral(hfsig->GetBinLowEdge(bmin), hfsig->GetBinLowEdge(bmax + 1));
      Integral_BW_withsigma = fitFcn2->Integral(hfsig->GetBinLowEdge(bmin), hfsig->GetBinLowEdge(bmax + 1));
      fYield_BinCount = Yield_bincount_hist - (bkgvalue / BinWidth);
      YieldIntegral_BW = fitFcn2->Integral(1.25, 2) / BinWidth;
      Yfraction_cBW = (Integral_BW_withsigma / YieldIntegral_BW);

      sum_tail_correction = (fitFcn2->Integral(1.25, hfsig->GetBinLowEdge(bmin)) + fitFcn2->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / BinWidth;
      Total_Ybincounting = (sum_tail_correction + fYield_BinCount) / (Event * dpt * dy * 0.33 * 2);
      Tail_correction_plusm = (fitFcn2_plusm->Integral(1.25, hfsig->GetBinLowEdge(bmin)) + (fitFcn2_plusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5))) / BinWidth;
      Tail_correction_minusm = ((fitFcn2_minusm->Integral(1.25, hfsig->GetBinLowEdge(bmin)) + fitFcn2_minusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / BinWidth);
      Error_2 = sum_tail_correction - Tail_correction_plusm;
      Final_pro_error = TMath::Sqrt(Error_2 * Error_2 + hBCError_1 * hBCError_1) / (Event * dpt * dy * 0.33 * 2);

      Yield_value_par = fitFcn->GetParameter(0) / (Event * dpt * dy);
      yield_error_par = fitFcn->GetParError(0) / (Event * dpt * dy);

      Yield_value_par1710 = fitFcn->GetParameter(3) / (Event * dpt * dy);
      yield_error_par1710 = fitFcn->GetParError(3) / (Event * dpt * dy);

      // cout<< sum_tail_correction << "\t" <<Tail_correction_plusm  <<"\t" <<  Tail_correction_minusm << "\t" << hBCError_1 << "\t" <<  Error_2 << endl;

      /*****************uncorrected Yield**********************/

      hYbincount->SetBinContent(ip + 1, Total_Ybincounting);
      hYbincount->SetBinError(ip + 1, Final_pro_error);

      hYieldpar1525->SetBinContent(ip + 1, Yield_value_par);
      hYieldpar1525->SetBinError(ip + 1, yield_error_par);

      hYieldpar1710->SetBinContent(ip + 1, Yield_value_par1710);
      hYieldpar1710->SetBinError(ip + 1, yield_error_par1710);

      // fractional stat error
      hFrac_stat_error->SetBinContent(ip + 1, Final_pro_error / Total_Ybincounting);

      {
        significance_num = (fitFcn2->Integral(1.25, 2.0)) / (BinWidth);
        significance_den = TMath::Sqrt(fHistTotal[ip]->Integral(bmin, bmax));
        ratio = significance_num / significance_den;

        cout << ip << "\t" << BinWidth << "\t" << significance_num << "\t" << significance_den << "\t" << ratio << (Final_pro_error / Total_Ybincounting) << endl;

        // out<<ip<<"\t" <<  BinWidth << "\t" << significance_num << "\t" <<  significance_den <<"\t"<< ratio << (Final_pro_error/Total_Ybincounting) <<  endl;
        hsgnfcance->SetBinContent(ip + 1, ratio);
      }

      // significance for all three resonances

      // one set of resonances : a2(1320)/f2(1320) , keep accorrding mass peak 1320. Mass = 1318 +- 0.6 , width = 107 +- 5, 1 sigma = 1318 +- 107, 1318 +- 2*107, 1318 +- 3.0*107

      // range for sigma calculation 1 sigma
      //  Double_t lr1 = 1.104;
      // Double_t hr1 = 1.532;

      // Double_t lr1 = 1.104;
      // Double_t hr1 = 1.532;

      Double_t lr1 = mass1320 - width1320;
      Double_t hr1 = mass1320 + width1320;

      // bmin1=hfsig->GetXaxis()->FindBin(1.1);
      // bmax1=hfsig->GetXaxis()->FindBin(1.53);
      Double_t bmin1, bmax1;
      bmin1 = hfsig->GetXaxis()->FindBin(lr1);
      bmax1 = hfsig->GetXaxis()->FindBin(hr1);

      Double_t splusb_r11 = fHistTotal[ip]->Integral(bmin1, bmax1);
      Double_t splusb_r1 = TMath::Sqrt(fHistTotal[ip]->Integral(bmin1, bmax1));

      Double_t signal_r1 = fitFcn2->Integral(lr1, hr1) / (BinWidth);
      // Double_t signal =  fitFcn2->Integral(1.2,1.4)/(BinWidth);
      Double_t signi_r1 = (signal_r1 / splusb_r1);

      cout << "bin--low" << bmin1 << "bin--high" << bmax1 << "\t" << "signal+bkg" << splusb_r11 << "\t----sqrt(s+b)" << splusb_r1 << "\t--signal" << signal_r1 << "\t" << signi_r1 << endl;

      // one set of resonances : 1525 , keep accorrding mass peak 1320. Mass = 1525 +- 2.5 , width = 86 +- 5, 1 sigma = 1525 +-86, 1525 +- 2*86

      Double_t lr2 = mass1525 - width1525;
      Double_t hr2 = mass1525 + width1525;

      // Double_t lr2 = 1.439;
      // Double_t hr2 = 1.611;

      // bmin1=hfsig->GetXaxis()->FindBin(1.1);
      // bmax1=hfsig->GetXaxis()->FindBin(1.53);
      Double_t bmin2, bmax2;
      bmin2 = hfsig->GetXaxis()->FindBin(lr2);
      bmax2 = hfsig->GetXaxis()->FindBin(hr2);

      Double_t splusb_r22 = fHistTotal[ip]->Integral(bmin2, bmax2);
      Double_t splusb_r2 = TMath::Sqrt(fHistTotal[ip]->Integral(bmin2, bmax2));

      Double_t signal_r2 = fitFcn2->Integral(lr2, hr2) / (BinWidth);
      // Double_t signal =  fitFcn2->Integral(1.2,1.4)/(BinWidth);
      Double_t signi_r2 = (signal_r2 / splusb_r2);

      cout << "bin--low" << bmin2 << "bin--high" << bmax2 << "\t" << "signal+bkg" << splusb_r22 << "\t----sqrt(s+b)" << splusb_r2 << "\t--signal" << signal_r2 << "\t" << signi_r2 << endl;

      // one set of resonances : 1710 , keep accorrding mass peak 1320. Mass = 1704 +- 12 , width = 123 +-18, 1 sigma = 1710 +-86, 1525 +- 2*86

      //  Double_t lr3 = 1.587;
      // Double_t hr3 = 1.833;

      Double_t lr3 = mass1710 - width1710;
      Double_t hr3 = mass1710 + width1710;

      // bmin1=hfsig->GetXaxis()->FindBin(1.1);
      // bmax1=hfsig->GetXaxis()->FindBin(1.53);
      Double_t bmin3, bmax3;
      bmin3 = hfsig->GetXaxis()->FindBin(lr3);
      bmax3 = hfsig->GetXaxis()->FindBin(hr3);

      Double_t splusb_r33 = fHistTotal[ip]->Integral(bmin3, bmax3);
      Double_t splusb_r3 = TMath::Sqrt(fHistTotal[ip]->Integral(bmin3, bmax3));

      Double_t signal_r3 = fitFcn2->Integral(lr3, hr3) / (BinWidth);
      // Double_t signal =  fitFcn2->Integral(1.2,1.4)/(BinWidth);
      Double_t signi_r3 = (signal_r3 / splusb_r3);

      cout << "bin--low" << bmin3 << "bin--high" << bmax3 << "\t" << "signal+bkg" << splusb_r33 << "\t----sqrt(s+b)" << splusb_r3 << "\t--signal" << signal_r3 << "\t" << signi_r3 << endl;
      // one set of resonances : a2(1320)/f2(1320) , keep accorrding mass peak 1320. Mass = 1525 +- 2.5 , width = 86 +- 5

      // one set of resonances : a2(1320)/f2(1320) , keep accorrding mass peak 1320. Mass = 1704 +- 12 , width = 123 +- 18

      hChiSquare->SetBinContent(ip + 1, chi);
      /*****setting plot parametres ***********/

      SetHistoStyle(hfsig, 1, 20, 0.6, 0.05, 0.05, 1.0, 1.1);
      SetHistoStyle(hfbkg, 2, 24, 0.5, 0.05, 0.05, 1.0, 1.1);
      SetHistoStyle(fHistTotal[ip], 1, 24, 0.5, 0.05, 0.05, 1.0, 1.1);
      hfsig->GetXaxis()->SetTitle("m_{K^{0}_{s}K^{0}_{s}} (GeV/c^{2})");
      hfsig->GetYaxis()->SetTitle(Form("Counts/(%2.0f MeV/c^{2})", 1000 * hfsig->GetBinWidth(10)));

      // hfsig->GetXaxis()->SetRangeUser(1.1,2.0);
      hfbkg->GetXaxis()->SetRangeUser(0.9, 2.5);
      hfbkg->GetXaxis()->SetTitle("m_{inv} (GeV/c^{2})");
      hfbkg->GetYaxis()->SetTitle(Form("Counts/(%2.0f MeV/c^{2})", 1000 * hfbkg->GetBinWidth(10)));

      fitFcn1->SetLineColor(4);
      fitFcn1->SetLineStyle(2);
      fitFcn1->SetLineWidth(2);
      fitFcn2->SetLineColor(6);
      fitFcn2->SetLineStyle(2);
      fitFcn2->SetLineWidth(2);

      // fitFcn3->SetLineColor(1);
      // fitFcn3->SetLineStyle(2);
      // fitFcn3->SetLineWidth(2);

      // fitFcn4->SetLineColor(3);
      // fitFcn4->SetLineStyle(2);
      // fitFcn4->SetLineWidth(2);

      // fitFcn5->SetLineColor(7);
      // fitFcn5->SetLineStyle(2);
      // fitFcn5->SetLineWidth(2);

      // fitFcn6->SetLineColor(41);
      // fitFcn6->SetLineStyle(2);
      //  fitFcn6->SetLineWidth(2);

      /***************************************/
      // plot
      cinv[ip]->cd();
      hfsig->GetXaxis()->SetRangeUser(1.1, 2.0);
      hfsig->Draw("e");
      fitFcn->Draw("same");
      fitFcn1->Draw("same");
      fitFcn2->SetFillColor(2);
      fitFcn2->Draw("same");

      //	    fitFcn3->Draw("same");
      // fitFcn4->Draw("same");
      // fitFcn5->Draw("same");
      // fitFcn6->Draw("same");
      sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], Low_pt[ip + 1]);
      if (ip < 1)
      {
        c2->cd(ip + 1);
        hfsig->Draw("e");
        hfsig->GetXaxis()->SetRangeUser(1.1, 2.0);
        TF1 *ftotal = (TF1 *)fitFcn->Clone();
        TF1 *fBW = (TF1 *)fitFcn2->Clone();
        TF1 *fpoly = (TF1 *)fitFcn1->Clone();
        ftotal->Draw("same");
        fBW->Draw("same");
        fpoly->Draw("same");

        TLatex *ltx = new TLatex(0.35, 0.93, name);
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.07187006);
        ltx->Draw();

        TLegend *leg = new TLegend(0.68, 0.7, 0.89, 0.83, NULL, "brNDC");
        leg->SetBorderSize(0.3);
        leg->SetTextFont(22);
        leg->SetTextSize(0.05);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->AddEntry(fitFcn, "4rBW + Res.Bkg", "l");
        leg->AddEntry(fBW, "4rBW", "l");
        leg->AddEntry(fpoly, "Res.Bkg", "l");
        leg->Draw();
      }
      else
      {
        c22->cd(ip - 11);
        hfsig->GetXaxis()->SetRangeUser(1.1, 2.0);
        hfsig->Draw("e");

        TF1 *ftotal = (TF1 *)fitFcn->Clone();
        TF1 *fBW = (TF1 *)fitFcn2->Clone();
        TF1 *fpoly = (TF1 *)fitFcn1->Clone();
        ftotal->Draw("same");
        fBW->Draw("same");
        fpoly->Draw("same");
        TLatex *ltx = new TLatex(0.35, 0.93, name);
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.07187006);
        ltx->Draw();

        TLegend *leg = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
        leg->SetBorderSize(0.3);
        leg->SetTextFont(22);
        leg->SetTextSize(0.05);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->AddEntry(fitFcn, "4rBW+Res.Bkg", "l");
        leg->AddEntry(fBW, "rBW", "l");
        leg->AddEntry(fpoly, "Res.Bkg", "l");
        leg->Draw();
      }

      // signal extraction before bkg subtraction :

      if (ip < 12)
      {
        c1->cd(ip + 1);
        hfbkg->Draw("e");
        fHistTotal[ip]->Draw("samee");
        TLatex *ltx = new TLatex(0.35, 0.93, name);
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.07187006);
        ltx->Draw();
        TLegend *leg = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
        leg->SetBorderSize(0.3);
        leg->SetTextFont(22);
        leg->SetTextSize(0.05);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
        leg->AddEntry(hfbkg, "Normalised mixed bkg", "p");
        leg->Draw();
      }
      else
      {
        c11->cd(ip - 11);
        fHistTotal[ip]->GetXaxis()->SetRangeUser(0.9, 2.5);
        fHistTotal[ip]->Draw("pe");
        hfbkg->Draw("pe same");
        TLegend *leg = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
        leg->SetBorderSize(0.3);
        leg->SetTextFont(22);
        leg->SetTextSize(0.05);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
        leg->AddEntry(hfbkg, "Normalised mixed bkg", "p");
        leg->Draw();
        TLatex *ltx = new TLatex(0.35, 0.93, name);
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.07187006);
        ltx->Draw();
      }
      cSigbkg[ip]->cd();
      fHistTotal[ip]->GetXaxis()->SetRangeUser(0.9, 2.5);
      fHistTotal[ip]->Draw("e");
      hfbkg->Draw("same pe");
      // covariant matrix for yielderror  calculation
      TMatrixDSym mat = r->GetCovarianceMatrix();
      // mat.Print();
      TMatrixDSym mat1;
      TMatrixDSym mat2;
      mat.GetSub(0, 2, 0, 2, mat1);
      mat.GetSub(3, 5, 3, 5, mat2);
      Double_t *b = mat1.GetMatrixArray();
      Double_t *a = mat2.GetMatrixArray();
      // wRawYield[ip] = fitFcn2->Integral(0.66,5.0)/(BinWidth*2*Event*dpt*dy*0.66);
      wRawYield[ip] = fitFcn2->Integral(1.25, 5.0) / (BinWidth * 2 * Event * dpt * dy * 0.33);
      wErRawYield[ip] = fitFcn2->IntegralError(1.25, 5.0, &par[0], b) / (BinWidth * Event * dpt * dy * 0.33 * 2);
      // cout <<  wRawYield[ip] << "\t" << Total_Ybincounting <<"\t" << Yield_value_par << "\t" << wErRawYield[ip]  << "\t" << Final_pro_error <<"\t" << yield_error_par<<endl;
      //  out <<  wRawYield[ip] << "\t" << Total_Ybincounting <<"\t" << Yield_value_par << "\t" << wErRawYield[ip]  << "\t" << Final_pro_error <<"\t" << yield_error_par<<endl;
    }
  }
  // return 0;
  /*************hsignificance plot************/
  TCanvas *cs = new TCanvas("cs", "", 10, 10, 600, 600);
  cs->cd();
  SetHistoStyle(hsgnfcance, 2, 8, 1.5, 0.05, 0.05, 1.0, 1.1);
  hsgnfcance->GetYaxis()->SetTitle("S/#sqrt{S+B}");
  hsgnfcance->GetXaxis()->SetTitle("p_{T} GeV/c");
  hsgnfcance->Draw("p");
  // hsgnfcance->Print();
  TLegend *leg = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hsgnfcance, "Significance", "p");
  leg->Draw("");

  // hsgnfcance-->>GetYaxis()->SetTitle("Significance(S/#Sqrt(S+B))");

  TCanvas *cmass1525 = new TCanvas("cmass1525", "cmass1525", 10, 10, 600, 600);
  cmass1525->cd();
  SetHistoStyle(hmass1525, 2, 20, 1.5, 0.05, 0.05, 1.0, 1.1);
  TLine *line = new TLine(1.0, 1.525, 10, 1.525);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  hmass1525->GetYaxis()->SetTitleFont(42);
  hmass1525->GetYaxis()->SetTitle("Mass(GeV/c^{2})");
  hmass1525->GetXaxis()->SetTitle("p_{T} GeV/c");
  hmass1525->Draw("p");
  line->Draw("p same");
  TCanvas *cmass1710 = new TCanvas("cmass1710", "cmass1710", 10, 10, 600, 600);
  cmass1710->cd();
  SetHistoStyle(hmass1710, 2, 20, 1.5, 0.05, 0.05, 1.0, 1.1);
  TLine *line = new TLine(1.0, 1.710, 10, 1.710);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  hmass1710->GetYaxis()->SetTitleFont(42);
  hmass1710->GetYaxis()->SetTitle("Mass(GeV/c^{2})");
  hmass1710->GetXaxis()->SetTitle("p_{T} GeV/c");
  hmass1710->GetYaxis()->SetRangeUser(1.6, 1.8);
  hmass1710->Draw("p");
  line->Draw("p same");
  TLegend *leg = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hmass, "Mass Peak", "p");
  leg->AddEntry(line, "PDG", "l");
  leg->Draw("");

  TCanvas *cwidth1525 = new TCanvas("cwidth1525", "", 10, 10, 600, 600);
  cwidth1525->cd();
  SetHistoStyle(hwidth1525, 2, 20, 1.5, 0.05, 0.05, 1.0, 1.1);
  TLine *line = new TLine(1, 0.073, 10, 0.073);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);

  hwidth1525->GetYaxis()->SetTitleFont(42);
  hwidth1525->GetYaxis()->SetTitle("Width(GeV/c^{2})");
  hwidth1525->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hwidth1525->GetYaxis()->SetRangeUser(0.0, 0.3);
  hwidth1525->Draw("p");
  line->Draw("psame");

  TLegend *leg = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hwidth1525, "f2(1525)", "p");
  leg->AddEntry(line, "PDG", "l");
  leg->Draw("");

  TCanvas *cwidth1710 = new TCanvas("cwidth1710", "", 10, 10, 600, 600);
  cwidth1710->cd();
  SetHistoStyle(hwidth1710, 2, 20, 1.5, 0.05, 0.05, 1.0, 1.1);
  TLine *line = new TLine(1, 0.139, 10, 0.139);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);

  hwidth1710->GetYaxis()->SetTitleFont(42);
  hwidth1710->GetYaxis()->SetTitle("Width(GeV/c^{2})");
  hwidth1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hwidth1710->GetYaxis()->SetRangeUser(0.0, 0.3);

  TLegend *leg = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hwidth1710, "f0(1710)", "p");
  leg->AddEntry(line, "PDG", "l");
  leg->Draw("");

  hwidth1710->Draw("p");
  line->Draw("psame");

  TCanvas *cmT = new TCanvas("cmT", "", 10, 10, 600, 600);
  SetHistoStyle(htransmass, 2, 20, 1.5, 0.05, 0.05, 1.0, 1.1);
  htransmass->GetYaxis()->SetTitleFont(42);
  htransmass->GetYaxis()->SetTitle("m_{T} GeV/c)");
  htransmass->GetXaxis()->SetTitle("p_{T} GeV/c");
  htransmass->Draw();

  TCanvas *cptspectra1525 = new TCanvas("cptspectra1525", "", 10, 10, 600, 600);
  cptspectra1525->cd();
  SetHistoStyle(hYieldpar1525, 2, 20, 1.5, 0.05, 0.05, 1.0, 1.1);
  hYieldpar1525->GetYaxis()->SetTitleFont(42);
  hYieldpar1525->GetYaxis()->SetTitle("Yield");
  hYieldpar1525->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hYieldpar1525->Draw("p");

  TLegend *leg = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hYieldpar1525, "f2(1525)", "p");
  leg->Draw("");

  TCanvas *cptspectra1710 = new TCanvas("cptspectra1710", "", 10, 10, 600, 600);
  cptspectra1710->cd();
  SetHistoStyle(hYieldpar1710, 2, 20, 1.5, 0.05, 0.05, 1.0, 1.1);
  hYieldpar1710->GetYaxis()->SetTitleFont(42);
  hYieldpar1710->GetYaxis()->SetTitle("Yield");
  hYieldpar1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hYieldpar1710->Draw("p");

  TLegend *leg = new TLegend(0.680554, 0.7312735, 0.892902, 0.8338954, NULL, "brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hYieldpar1710, "f0(1710)", "p");
  leg->Draw("");

  // residual background sutraction method :

  TH1F *hSigFit1_SResBkg = (TH1F *)hfsig_new->Clone("hSigFit1new");
  hSigFit1_SResBkg->Reset();

  TH1F *hSig = (TH1F *)hfsig_new->Clone("hSigFit1_new");

  for (Int_t ip = 0; ip < hSigFit1_SResBkg->GetNbinsX(); ip++)
  {
    float x = hSigFit1_SResBkg->GetBinCenter(ip + 1);
    // float x=hSigFit1_SResBkg->GetBinLowEdge(ip+1);
    float y = fitFcn1->Eval(x); // value = hSigTemp1->GetBinContent(ip+1) - (fBgOnly1->Integral(ip+1,ip+2)/hSigFit1_SResBkg->GetBinWidth(2));
    //  cout << "\----ip " <<  ip <<  "---position__x" << x  <<  "bkg-----" <<  y << "\t signal------" << hSig->GetBinContent(ip+1) << endl;
    hSigFit1_SResBkg->SetBinContent(ip + 1, y);
  }

  //  TF1 *fSigOnly2 = new TF1("fitFcn2",rBW,fitlow,fithi,14);

  //  TF1 *fSigOnly2 = new TF1("fitFcn2",rBW,1.24,1.88,14);

  //  TF1 *fSigOnly2 = new TF1("fitFcn2",rBW,1.21,1.96,14);

  //  TF1 *fSigOnly2 = new TF1("fitFcn2",rBW,1.23,2.5,14);
  //  TF1 *fSigOnly2 = new TF1("fitFcn2",rBW,1.21,1.99,14);
  TF1 *fSigOnly2 = new TF1("fitFcn2", rBW, 1.12, 2.3, 14);

  fSigOnly2->SetLineColor(kRed);
  fSigOnly2->SetParameter(0, 1000);      // yield
  fSigOnly2->SetParLimits(1, 1.16, 1.4); // inv. mass peak range
  fSigOnly2->SetParLimits(2, 0.1, 0.26); // width-fixed

  fSigOnly2->SetParameter(3, 1000);       // yield
  fSigOnly2->SetParLimits(4, 1.26, 1.5);  // inv. mass peak range
  fSigOnly2->SetParLimits(5, 0.02, 0.26); // width-fixed

  fSigOnly2->SetParameter(6, 1000);      // yield
  fSigOnly2->SetParLimits(7, 1.46, 1.6); // inv. mass peak range
  fSigOnly2->SetParLimits(8, 0.01, 0.2); // width-fixed

  fSigOnly2->SetParameter(9, 100);        // yield
  fSigOnly2->SetParLimits(10, 1.6, 1.9);  // inv. mass peak range
  fSigOnly2->SetParameter(11, 0.123);     // width-fixed
  fSigOnly2->SetParLimits(11, 0.01, 0.2); // width-fixedfSigOnly2->SetParLimits(11,0.01,0.2); //width-fixed

  hSig->Add(hSigFit1_SResBkg, -1);
  fSigOnly2->SetNpx(10000);

  r = hSig->Fit(fSigOnly2, "REBMS+");
  hSig = DrawFrame(hSig, STCol, STSty, 1);
  hSig->GetXaxis()->SetRangeUser(1.12, 2.3);
  hSig->SetMarkerSize(0.8);
  hSig->GetYaxis()->SetRangeUser(-1000.0, 17000.0);
  // hSigFit1_SResBkg->Draw();
  hSig->Draw("p");
  hSig->SetTitle("");
  hSig->SetMarkerSize(1.0);
  // fSig->SetNpx(10000);
  //   hSig->GetXaxis()->SetRangeUser(1.22,1.92);
  //   hSig->GetXaxis()->SetRangeUser(1.24,2.1);

  //   fSigOnly1->Draw("p same");

  TLegend *lp2 = DrawLegend(0.6, 0.66, 0.88, 0.9);
  lp2->SetTextSize(.05);
  // lp2->AddEntry((TObject*)0,"#bf{ALICE Performance}","");
  lp2->AddEntry((TObject *)0, "pp #sqrt{#it{s}} = 13 TeV (0-100%)", "");
  // lp2->AddEntry((TObject*)0,"V0M (0-100%)","");
  lp2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
  lp2->AddEntry((TObject *)0, Form("%s", ptbn.Data()), "");
  // lp2->Draw("same");

  TLegend *lp3 = DrawLegend(0.8, 0.8, 0.9, 0.9);
  // lp3->AddEntry(hSig,"Data","p");
  lp3->AddEntry(fSigOnly2, "4rBW", "l");
  lp3->SetTextSize(.05);
  lp3->Draw("same");

  // Draw the function only
  TCanvas *csignal = new TCanvas("csignal", "", 10, 10, 600, 600);
  csignal->cd();

  TF1 *fSigOnly23 = new TF1("fitFcn2", rBW3, 1.12, 2.2, 9);

  fSigOnly23->FixParameter(0, 733.3);  // yield
  fSigOnly23->FixParameter(1, 1.251);  // inv. mass peak range
  fSigOnly23->FixParameter(2, 0.1542); // width-fixed

  fSigOnly23->FixParameter(3, 1072);    // yield
  fSigOnly23->FixParameter(4, 1.518);   // inv. mass peak range
  fSigOnly23->FixParameter(5, 0.09986); // width-fixed

  fSigOnly23->FixParameter(6, 414);    // yield
  fSigOnly23->FixParameter(7, 1.692);  // inv. mass peak range
  fSigOnly23->FixParameter(8, 0.1341); // width-fixed

  fSigOnly23->SetLineColor(4);

  fSigOnly23->Draw("l");
  fSigOnly2->Draw(" same l");

  TLegend *lp3 = DrawLegend(0.6, 0.6, 0.8, 0.8);
  // lp3->AddEntry(hSig,"Data","p");
  lp3->AddEntry(fSigOnly2, "4rBW", "l");
  lp3->AddEntry(fSigOnly23, "3rBW", "l");
  lp3->SetTextSize(.05);
  lp3->Draw("same");

  TCanvas *cbkg = new TCanvas("cbkg", "", 10, 10, 600, 600);
  cbkg->cd();

  TF1 *fBkgOnly23 = new TF1("fitFcn2", polynomial2, 1.05, 2.2, 9);

  fBkgOnly23->FixParameter(0, 191400.3); // yield
  fBkgOnly23->FixParameter(1, 0.054);    // inv. mass peak range
  fBkgOnly23->FixParameter(2, 4.872);    // width-fixed

  fBkgOnly23->SetLineColor(4);

  fitFcn1->Draw("l");
  fBkgOnly23->Draw(" same l");

  TLegend *lp3 = DrawLegend(0.4, 0.6, 0.6, 0.8);
  // lp3->AddEntry(hSig,"Data","p");
  lp3->AddEntry(fitFcn1, " Res.Bkg: 4rBW", "l");
  lp3->AddEntry(fBkgOnly23, "Res.Bkg: 3rBW", "l");
  lp3->SetTextSize(.05);
  lp3->Draw("same");
}

TCanvas *DrawCanvas(TString opt = "c")
{
  TCanvas *c1 = new TCanvas(opt.Data(), opt.Data(), 10, 10, 600, 600);
  c1->cd(1);
  c1->SetLeftMargin(0.25);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
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
TGraphErrors *PlotGraph(Int_t NdataPoint, Int_t MarkerColor = 1, Int_t MarkerStyle = 20, Double_t *X, Double_t *ErX, Double_t *Y, Double_t *ErY)
{
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
  gr->GetYaxis()->SetTitleOffset(1.7);
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
TLatex *DrawText(Double_t x = 0, Double_t y = 0, Int_t tColor = 2, TString name)
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

void SetHistoStyle(TH1 *h, Int_t mcolor, Int_t mstyle, Float_t msize, Float_t Tsizex, Float_t Tsizey, Float_t Offsetx, Float_t Offsety)
{
  h->SetMarkerColor(mcolor);
  h->SetMarkerStyle(mstyle);
  h->SetMarkerSize(msize);
  h->GetXaxis()->SetTitleSize(Tsizex);
  h->GetYaxis()->SetTitleSize(Tsizey);
  h->GetXaxis()->SetTitleOffset(Offsetx);
  h->GetYaxis()->SetTitleOffset(Offsety);
}

void *DrawFrame(TH1 *h, Int_t MCol, Int_t MSty, Bool_t mrk = 0)
{
  //  h->GetXaxis()->SetTitle("#it{M}_{K^{0}_{S}#pi^{#pm}} (GeV/#it{c^{2}})");

  h->GetXaxis()->SetTitle("#it{M}_{K^{0}_{S}K^{0}_{S}} (GeV/#it{c^{2}})");
  if (mrk)
  {
    h->SetMarkerColor(MCol);
    h->SetMarkerStyle(MSty);
    h->SetMarkerSize(0.5);
  }
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
  h->GetXaxis()->SetTickLength(0.01);
  h->GetYaxis()->SetTickLength(0.01);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleOffset(1.65);
  // h->GetYaxis()->CenterTitle(true);
  h->GetYaxis()->SetDecimals(true);
  // h->GetYaxis()->SetNdivisions(310);
  h->GetYaxis()->SetLabelOffset(0.015);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTickLength(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitle("Counts / (0.02 GeV/#it{c^{2}})");
  h->GetYaxis()->SetTitleFont(42);
  h->Draw("");
  return h;
}

void SetHistoStyle(TH1 *h, Int_t mcolor, Int_t mstyle, Float_t msize, Float_t Tsizex, Float_t Tsizey, Float_t Offsetx, Float_t Offsety)
{
  h->SetMarkerColor(mcolor);
  h->SetMarkerStyle(mstyle);
  h->SetMarkerSize(msize);
  h->GetXaxis()->SetTitleSize(Tsizex);
  h->GetYaxis()->SetTitleSize(Tsizey);
  h->GetXaxis()->SetTitleOffset(Offsetx);
  h->GetYaxis()->SetTitleOffset(Offsety);
}

Double_t rBW3(Double_t *x, Double_t *par)
{

  double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
  double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
  double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
  double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
  double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

  double Gamma = par[2];
  double Gamma1 = par[5];
  double Gamma2 = par[8];

  double rBW1 = par[0] * x[0] * par[1] * Gamma / (TMath::Power((x[0] * x[0] - par[1] * par[1]), 2.0) + par[1] * par[1] * Gamma * Gamma);
  double rBW2 = par[3] * x[0] * par[4] * Gamma1 / (TMath::Power((x[0] * x[0] - par[4] * par[4]), 2.0) + par[4] * par[4] * Gamma1 * Gamma1);
  double rBW3 = par[6] * x[0] * par[7] * Gamma2 / (TMath::Power((x[0] * x[0] - par[7] * par[7]), 2.0) + par[7] * par[7] * Gamma2 * Gamma2);
  double sigfun = rBW1 + rBW2 + rBW3;

  return (sigfun);
}
