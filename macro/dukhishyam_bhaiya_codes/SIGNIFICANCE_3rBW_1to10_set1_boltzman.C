Double_t rBreitWigner(Double_t *x, Double_t *par)
{

  ///double rBW1 = par[0]*par[2]/(2*3.14159)/((x[0]-par[1])**2+par[2]**2/4.);
  //double rBW2 = par[3]*par[5]/(2*3.14159)/((x[0]-par[4])**2+par[5]**2/4.);
  // double rBW3 = par[6]*par[8]/(2*3.14159)/((x[0]-par[7])**2+par[9]**2/4.);
  // double rBW4 = par[9]*par[11]/(2*3.14159)/((x[0]-par[10])**2+par[11]**2/4.);

  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);
  double  dpart2 = par[4]*par[4] - 4*(0.4976*0.4976);
  double  dpart3 = par[7]*par[7] - 4*(0.4976*0.4976);
  double  dpart4 = par[10]*par[10] - 4*(0.4976*0.4976);
 
  double Gamma = par[2];
  double Gamma1 = par[5];
  double Gamma2 = par[8];
  // double Gamma3 = par[11];

  double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  double rBW3 = par[6]*x[0]*par[7]*Gamma2/(TMath::Power((x[0]*x[0]-par[7]*par[7]),2.0)+par[7]*par[7]*Gamma2*Gamma2);
  //  double rBW4 = par[9]*x[0]*par[10]*Gamma2/(TMath::Power((x[0]*x[0]-par[10]*par[10]),2.0)+par[10]*par[10]*Gamma2*Gamma2);

  double sigfun = rBW1+rBW2+rBW3;
  //  double poly2 = par[9]*TMath::Power(x[0],par[10])*TMath::Exp(-x[0]*par[11]);


  double poly2 = par[9]*TMath::Power(x[0]-0.994,par[10]/2.0)*TMath::Power(par[11],3.0/2.0)*TMath::Exp(-par[11]*(TMath::Power(x[0]-0.994,par[10])));


  return(sigfun+poly2);
    
}

Double_t rBW(Double_t *x, Double_t *par)
{
  //  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  // double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);

  //double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  //double rBW = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  // double BF = TMath::Sqrt((x[0]*x[0]+par[3]*par[3]));
  // double PS =  x[0]*TMath::Exp(-BF/par[4])/(BF);

  // double rBW = ;

  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);
  double  dpart2 = par[4]*par[4] - 4*(0.4976*0.4976);
  double  dpart3 = par[7]*par[7] - 4*(0.4976*0.4976);
  double  dpart4 = par[10]*par[10] - 4*(0.4976*0.4976);
  
  //  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  // double Gamma1 = par[5]*(TMath::Power(par[4]/x[0],1.0))*TMath::Power((npart1)/(dpart2),1.5);
  // double Gamma2 = par[8]*(TMath::Power(par[7]/x[0],1.0))*TMath::Power((npart1)/(dpart3),1.5);
  // double Gamma3 = par[11]*(TMath::Power(par[10]/x[0],1.0))*TMath::Power((npart1)/(dpart4),1.5);


  
  //double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  //  double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  // double rBW3 = par[6]*x[0]*par[7]*Gamma2/(TMath::Power((x[0]*x[0]-par[7]*par[7]),2.0)+par[7]*par[7]*Gamma2*Gamma2);
  // double rBW4 = par[9]*x[0]*par[10]*Gamma3/(TMath::Power((x[0]*x[0]-par[10]*par[10]),2.0)+par[10]*par[10]*Gamma3*Gamma3);
  // double rBW1 = par[0]*par[1]*par[1]*Gamma/(TMath::Power(( -x[0]*x[0]+par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double Gamma = par[2];
  double Gamma1 = par[5];
  double Gamma2 = par[8];
  //double Gamma3 = par[11];
  
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

   
    double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
    double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
    double rBW3 = par[6]*x[0]*par[7]*Gamma2/(TMath::Power((x[0]*x[0]-par[7]*par[7]),2.0)+par[7]*par[7]*Gamma2*Gamma2);
    //  double rBW4 = par[9]*x[0]*par[10]*Gamma2/(TMath::Power((x[0]*x[0]-par[10]*par[10]),2.0)+par[10]*par[10]*Gamma2*Gamma2);

    double sigfun = rBW1+rBW2+rBW3;
      
      //      return(rBW1+rBW2+rBW3+rBW4);
      return(sigfun);


}



Double_t BreitWigner(Double_t *x, Double_t *par) 
{  
  return (par[0]*par[2]/(2*3.14159)/((x[0]-par[1])**2+par[2]**2/4.)+ par[3]+x[0]*par[4]+ x[0]*x[0]*par[5]);
  
}

Double_t polynomial2(Double_t *x, Double_t *par) 
{ 
  // double poly2 = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  
  //  double poly2 = par[0]*TMath::Power(x[0],par[1])*TMath::Exp(-x[0]*par[2]);
  
  //  double poly2 = par[3]*TMath::Power(x[0],par[4])*TMath::Exp(-x[0]*par[5]);
  
  //  double poly2 = par[0]*TMath::Sqrt(TMath::Power(x[0] -(2*0.497),par[1]))*TMath::Power(par[2],1.5)*TMath::Exp(-par[2]*(TMath::Power(x[0]-(2*0.497)),par[1]));
  // double poly2 = par[3]*exp(x[0]*par[4]);

  // double poly2 = par[0]*exp(x[0]*par[1])+par[2]*exp(x[0]*par[3]*x[0])+par[4]*exp(x[0]*x[0]*x[0]*par[5]);

  double poly2 = par[0]*TMath::Power(x[0]-0.994,par[1]/2.0)*TMath::Power(par[2],3.0/2.0)*TMath::Exp(-(TMath::Power(x[0]-0.994),par[1])*par[2]);
   return (poly2);
  // return (bg);

  //  return ((x[0]- 0.63718)**par[3])*exp(par[0] + x[0]*par[1] + x[0]*x[0]*par[2]);
}

Double_t BW(Double_t *x, Double_t *par) 
{
  return  par[0]*par[2]/(2*3.14159)/((x[0]-par[1])**2+par[2]**2/4.);
  //return (0.5*par[0]*par[1]/TMath::Pi() /((x[0]-par[2])*(x[0]-par[2]) + 0.25*par[1]*par[1]));
}

Double_t dy=1.0;
//const Int_t Npt=7;
const Int_t Npt=1;
const Int_t N=1;

//Double_t lownorm = 1.1;
// Double_t lownorm = 1.2;
//Double_t hinorm = 1.15;

Double_t lownorm = 2.2;
// Double_t lownorm = 1.2;
Double_t hinorm = 2.3;

//Double_t lownorm1= 1.1;
// Double_t lownorm = 1.2;
//Double_t hinorm1= 1.1;

Double_t lownorm1= 2.2;
// Double_t lownorm =1.2;
Double_t hinorm1= 2.3;


// Double_t hinorm = 1.3;
// Double_t hinorm = 1.3; // n0
// Double_t hinorm = 1.4; // n1
//Double_t hinorm = 1.3; // n2
Double_t lowsigmarange=1.6;
Double_t highsigmarange = 2.0;
void SIGNIFICANCE_3rBW_1to10_set1_boltzman()
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
  
  Double_t  fRawYield[Npt];  
  Double_t  fErRawYield[Npt];
  Double_t  wRawYield[Npt];  
  Double_t  wErRawYield[Npt];
  
  Double_t d,c,dpt,mass,error_mass, width, error_width;
  Double_t chi;
  TFitResultPtr r;
   TH1D *fHistTotal[Npt];
   TH1D *fHistBkg[Npt];
   TH1D *fHistlike[Npt];
   TH1D *histPP[Npt];
   TH1D *histMM[Npt];
   char name[100];
   
   Double_t yield,Yield, mass[Npt], fyield[Npt], Width[Npt], yerror[Npt], xerror[Npt], xerr[Npt], yerr[Npt], dn[Npt];
   Double_t  poly0[Npt], poly1[Npt],poly2[Npt],poly3[Npt];
   Double_t Par[Npt],pt,hBCError_1;

   Double_t yield_error_par,bmin,bmax,Yield_bincount_hist,Integral_BW_withsigma,fYield_BinCount,YieldIntegral_BW,Yfraction_cBW,Yield_value_par,sum_tail_correction,sum_tail_correction,Total_Ybincounting ,Tail_correction_plusm,Tail_correction_minusm ,Error_2,Final_pro_error,bkgvalue,Yield_value_par1710,yield_error_par1710;

 Double_t lc,hc,low_value,high_value,lbin,hbin,sig_integral,sig_integral1,bkg_integral1, bkg_integral,mT,mT_error, significance_num ,significance_den,ratio,BinWidth;

   //   Double_t Low_pt[Npt+1] ={0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,10.0,12.0,15.0,20.0};

   //   Double_t Low_pt[Npt+1] ={0.8,1.2,1.6,1.8,2.0,2.4,2.8,3.2,4.0,5.0};

   //   Double_t Low_pt[Npt+1] ={0.8,1.2,1.6,1.8,2.0,2.4,2.8,3.2,4.0,5.0};

   //   Double_t Low_pt[Npt+1] ={0.5,0.7,1.2,1.6,2.0,2.5,3.0,3.5,4.0,5.0,7.0,10.0,15.0,20.0};

   //   Double_t Low_pt[Npt+1] ={0.8,1.2,1.6,2.0,2.5,3.0,4.0,5.0};
   
   //   Double_t Low_pt[Npt+1] ={1.0,2.0,3.0,4.0,6.0,8.0,10.0,15.0,20.0};
   
   //   Double_t Low_pt[Npt+1] ={0.0,1.0,1.2,1.6,1.8,2.0,2.5,3.0,4.0};
   
   // Double_t Low_pt[Npt+1] ={1.0,2.0,3.0,4.0,5.0,6.0,8.0,12.0};
 
   // Double_t Low_pt[Npt+1] ={1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0};

 Double_t Low_pt[Npt+1] ={1.0,10.0};

  // Double_t Low_pt[Npt+1] ={0.0,10.0};
  Double_t Low_Centr[N]={1};
  Double_t High_Centr[N]={100};

  //0-100 % 
  // Double_t fitlow[Npt] = {.74,.72,.72,.66,.66,.66,.66,.66,.66};
  // Double_t fithigh[Npt]= {1.07,1.1,1.1,1.1,1.09,1.09,1.09,1.09,1.1};
  //10-30 %
  //Double_t fitlow[Npt] = {.74,.72,.72,.66,.66,.7,.66,.66,.66};
  //Double_t fithigh[Npt]= {1.07,1.1,1.1,1.1,1.09,1.09,1.09,1.09,1.1};
  
   //1-10 %
   //Double_t fitlow[Npt] = {.74,.72,.72,.66,.66,.7,.66,.66,.66};
   //Double_t fithigh[Npt]= {1.07,1.1,1.1,1.1,1.09,1.09,1.09,1.09,1.1};
  
   //Double_t fitlow[Npt] = {.75,.72,.75,.66,.66,.72,.66,.66,.66};
   //Double_t fithigh[Npt]= {1.09,1.1,1.1,1.12,1.09,1.09,1.09,1.09,1.1};


   // Double_t fitlow[Npt] = {.75,.72,.75,.7,.66,.72,.66,.66,.66};
   // Double_t fithigh[Npt]= {1.09,1.1,1.1,1.09,1.05,1.09,1.09,1.09,1.1};
 
   //0 to 30 
   // Double_t fitlow[Npt] = {.75,.72,.75,.67,.67,.72,.66,.66,.66};
   // Double_t fithigh[Npt]= {1.09,1.1,1.1,1.09,1.05,1.09,1.09,1.09,1.1};
 
  //70 to 90 
  //Double_t fitlow[Npt] = {.74,.72,.71,.67,.67,.72,.66,.66,.66};

  // 0 to 90 kstar binning
  // Double_t fitlow[Npt] = {.67,.67,.68,.66,.67,.72,.66,.66,.66};
  // Double_t fithigh[Npt]= {1.1,1.1,1.1,1.09,1.05,1.09,1.09,1.09,1.1};
  
 /*   Double_t fitlow[Npt] = {.71,.67,.68,.66,.71,.71,.72};
      Double_t fithigh[Npt]= {1.1,1.14,1.09,1.05,1.055,1.17,1.07};
 */
  
  
 // 20 -40 %
 // Double_t fitlow[Npt] = {.71,.7,.66,.69,.66,.66,.66};
  //Double_t fithigh[Npt]= {1.1,1.1,1.09,1.11,1.1,1.1,1.13};
  // 40 -60 %
  //   Double_t fitlow[Npt] = {.67,.66,.67,.69,.66,.66,.69};
  // Double_t fithigh[Npt]= {1.1,1.19,1.13,1.11,1.1,1.1,1.1};
  
  
  //   Double_t fitlow[Npt] = {.72,.66,.67,.69,.66,.66,.69};
   // Double_t fithigh[Npt]= {1.06,1.19,1.13,1.11,1.1,1.1,1.1};
   
   //Double_t fitlow[Npt] = {.72,.72,.72,.72,.72,.72,.72};
   //Double_t fithigh[Npt]= {1.06,1.06,1.06,1.06,1.06,1.06,1.06};
      
   //   Double_t fitlow[Npt] = {1.4,1.42,1.42,1.42,1.42,1.42,1.42,1.42};
   // Double_t fithigh[Npt]= {1.6,1.7,1.7,1.7,1.7,1.7,1.7,1.7};
   
   
   //   Double_t fitlow[Npt] = {1.48,1.45,1.4,1.37,1.4,1.45,1.45};
   // Double_t fithigh[Npt]= {1.6,1.65,1.65,1.7,1.7,1.75,1.75};
   

   //   Double_t fitlow[Npt] = {1.46,1.45,1.4,1.42,1.38,1.4,1.42};
   // Double_t fithigh[Npt]= {1.6,1.68,1.68,1.7,1.7,1.72,1.72};
   
   //   Double_t fitlow[Npt] = {1.46,1.45,1.4,1.42,1.38,1.4,1.42};
   //Double_t fithigh[Npt]= {1.6,1.68,1.68,1.7,1.7,1.72,1.72};
   
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
   //Double_t fithigh[Npt]= {2.0};
   
   
   //   Double_t fitlow[Npt] = {1.15, 1.15,1.15,1.15,1.15,1.15,1.15};
   // Double_t fithigh[Npt]= {2.0,2.0,2.0,2.0,2.0,2.0,2.0};
   
   //   Double_t fitlow[Npt] = {1.24, 1.24,1.24,1.4,1.4,1.4,1.4};
   //Double_t fithigh[Npt]= {2.0,1.95,1.9,1.7,1.7,1.7,1.7};


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
  //Double_t fithigh[Npt]= {1.96,1.87,1.95,1.95,1.95};
  
  //Double_t fitlow[Npt] = {1.25};
  // Double_t fithigh[Npt]= {2.0};


  Double_t fitlow[Npt] = {1.12};
  Double_t fithigh[Npt]= {2.3};

  // Double_t fitlow[Npt] = {1.2,1.3,1.3,1.26,1.26};
  // Double_t fithigh[Npt]= {1.95,1.87,1.87,1.95,1.95};

  //  Double_t fitlow[Npt] = {1.2,1.3,1.3,1.26,1.26};
  // Double_t fithigh[Npt]= {1.96,1.87,1.88,1.95,1.95};

  // fitting range new
  // Double_t fitlow[Npt] = {1.2,1.3,1.25,1.26,1.26};
   //Double_t fithigh[Npt]= {1.95,1.85,1.9,1.95,1.95};
  //Double_t fitlow[Npt] = {1.2,1.3,1.25,1.26,1.26};
  //Double_t fithigh[Npt]= {1.95,1.85,1.85,1.95,1.95};
  //Double_t fitlow[Npt] = {1.2,1.3,1.2,1.23,1.23};
  //Double_t fithigh[Npt]= {1.95,1.85,1.8,1.9,1.9};
  //Double_t fitlow[Npt] = {1.2,1.3,1.26,1.26,1.26};
  //Double_t fithigh[Npt]= {1.9,1.85,1.97,1.97,1.99};
  //Double_t Para1[Npt]=  {100, 1.0e+02, 1.0e+02, 1.0e+02,100};
  
  Double_t Para1[Npt]=  {100};

  // Double_t fitlow[Npt] = {1.42};
   //Double_t fithigh[Npt]= {1.7};

   //Double_t fitlow[Npt] = {.75};
   // Double_t fithigh[Npt]= {1.09};
   
   /*
   Double_t Para1[Npt]=  {100, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02};
   
   Double_t RBPar1[Npt]= {10,1,1,1,1,1,1};

   Double_t RBPar2[Npt]= {10,10 ,10 ,10 ,10,10,10};
   
   Double_t RBPar3[Npt]= {4,4,4,4,4,4,4};

   Double_t RBPar4[Npt]= {10,10,10,10,10,10,10};
   */  
   // Double_t Para1[Npt]=  {1000};
  //Double_t RBPar1[Npt]= {1000};
  //Double_t RBPar2[Npt]= {-100};
  //Double_t RBPar3[Npt]= {20};
  // Double_t RBPar4[Npt]= {10};

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
   
   //fitFcn->SetParameters(4000000,0.8916,0.0508, 10, 4, 10, 1);
   
   //canvas
    TCanvas *cSigbkg[Npt];
    TCanvas *cinv[Npt];
    
    TCanvas*c1=new TCanvas("c1","",10,10,600,600);
    c1->Divide(1,1);
    TCanvas*c11=new TCanvas("c11","",10,10,600,600);
    c11->Divide(1,1);
    TCanvas*c2=new TCanvas("c2","",10,10,600,600);
    c2->Divide(1,1);
    TCanvas*c22=new TCanvas("c22","",10,10,600,600);
    c22->Divide(1,1);
    for(Int_t ip = 0; ip <Npt; ip++)

      // for(Int_t ip = 0; ip <2; ip++)
      {
	TString cName = TString::Format("cSigbkg_pt_%2.1f-%2.1f",Low_pt[ip],Low_pt[ip+1]);
	cSigbkg[ip] = new TCanvas(Form("cSigbkg%d",ip),cName.Data(),10,10,600,600); 
	cSigbkg[ip]->SetLeftMargin(0.2);
	cSigbkg[ip]->SetRightMargin(0.05);
	cSigbkg[ip]->SetBottomMargin(0.2); 
      }
    
    for(Int_t ip = 0; ip <Npt; ip++)
      //  for(Int_t ip = 0; ip <2; ip++)
      {  
	TString cName = TString::Format("cinv_pt_%2.1f-%2.1f",Low_pt[ip],Low_pt[ip+1]);
	cinv[ip] = new TCanvas(Form("cinv%d",ip),cName.Data(),10,10,600,600); 
	
	cinv[ip]->SetRightMargin(0.05);
	cinv[ip]->SetBottomMargin(0.2);
     }
    
    TH1F*hYbincount = new TH1F("hYbincount","",Npt ,Low_pt);
    TH1F*hYieldpar = new TH1F("hYieldpar","",Npt ,Low_pt);
    TH1F*hChiSquare = new TH1F("hChiSquare","",Npt ,Low_pt);
    TH1F*hsgnfcance = new TH1F("hsgnfcance","",Npt ,Low_pt);
    TH1F*hmass = new TH1F("hmass","",Npt ,Low_pt);
    TH1F*hwidth = new TH1F("hwidth","",Npt ,Low_pt);

    TH1F*hmass1525 = new TH1F("hmass1525","hmass1525",Npt ,Low_pt);

    TH1F*hwidth1525 = new TH1F("hwidth1525","hwidth1525",Npt ,Low_pt);

    TH1F*hmass1710 = new TH1F("hmass1710","hmass1710",Npt ,Low_pt);

    TH1F*hwidth1710 = new TH1F("hwidth1710","hmass1710",Npt ,Low_pt);

    TH1F*htransmass = new TH1F("htransmass","",Npt ,Low_pt);
    TH1F*hFrac_stat_error = new TH1F("hFrac_stat_error","",Npt,Low_pt);
    
    TH1F*hYieldpar1525 = new TH1F("hYieldpar1525","",Npt ,Low_pt);

    TH1F*hYieldpar1710 = new TH1F("hYieldpar1710","",Npt ,Low_pt);
    
    TFile *fInputFile = new TFile("../../../AnalysisResultsKshortKshortpp13TeV_fulldataset_aod.root","Read"); 
    
    TKey* key = (TKey*)fInputFile->GetListOfKeys()->At(1);
    fInputList = (TList*)fInputFile->Get(key->GetName());
    //fInputList->ls();
    
    //return 0;
   for(Int_t j=0;j<N;j++) //centrality loop
     {
	Double_t lc=Low_Centr[j];
	Double_t hc=High_Centr[j];
	TH1F* hEVent = (TH1F *) fInputList->FindObject("hAEventsVsMulti");//Event calculation

	Double_t Event=hEVent->Integral(lc,hc);

	cout<<Event<<"\t"<<"\t"<<lc<<"\t"<< hc<<endl;
	
	//get invariantmass histogram from sparse histogram//MIXED EVENT
	//  THnSparseF * fHistNum = (THnSparseF *) fInputList->FindObject("KStarPlusMinusPbPbData_ChargeKstar_KStarPlusMinus");
	//	THnSparseF * fHistNum1 = (THnSparseF *) fInputList->FindObject("KStarPlusMinusPbPbData_ChargeKstar_AKStarPlusMinus");

	//	fHistNum->Add(fHistNum1);
	
	//	THnSparseF * fHistDen = (THnSparseF *) fInputList->FindObject("KStarPlusMinusPbPbData_ChargeKstar_KStarPlusMinusmix");
	//	THnSparseF * fHistDen1 = (THnSparseF *) fInputList->FindObject("KStarPlusMinusPbPbData_ChargeKstar_AKStarPlusMinusmix");
	//	fHistDen->Add(fHistDen1);

	THnSparseF * fHistNum = (THnSparseF *) fInputList->FindObject("KStarPlusMinusppData_ChargeKstar_KStarPlusMinus");
	//	THnSparseF * fHistNum1 = (THnSparseF *) fInputList->FindObject("KStarPlusMinusppData_ChargeKstar_AKStarPlusMinus");
	//	fHistNum->Add(fHistNum1);

	THnSparseF * fHistDen = (THnSparseF *) fInputList->FindObject("KStarPlusMinusppData_ChargeKstar_KStarPlusMinusmix");
	//	THnSparseF * fHistDen1 = (THnSparseF *) fInputList->FindObject("KStarPlusMinusppData_ChargeKstar_AKStarPlusMinusmix");
	//	fHistDen->Add(fHistDen1);
	
	
	for(Int_t m=0;m<Npt ;m++)
	//for(Int_t m=0;m<2 ;m++)
	  {
	    sprintf(name ,"fHistNum%d",m);
	    fHistTotal[m] = new TH1D(name,"inv_mass",200,0.5,2.5);
	    sprintf(name ,"fHistbkg%d",m);
	    fHistBkg[m]= new TH1D(name,"inv_mass",200,0.5,2.5);
	    sprintf(name,"histPP_pt%d",m);
	    
	  }
	for(Int_t ip = 0; ip <Npt; ip++)//start pt bin loop
	  // for(Int_t ip = 0; ip <3; ip++)//start pt bin loop
	  //	for(Int_t ip =2; ip <3; ip++)//start pt bin loop
	  {
	    low_value = Low_pt[ip];
	    high_value =Low_pt[ip+1]; ;
	    
      	    sprintf(name ,"fHistbkg%d",ip);
	    sprintf(name ,"fHistNum%d",ip);
	    
	    fHistNum->GetAxis(2)->SetRange(lc,hc);
	    fHistNum->GetAxis(2)->SetBit(TAxis::kAxisRange);
	    
	    lbin = fHistNum->GetAxis(1)->FindBin(low_value+1.e-7);//projection pt bin //histogram sig+bkg
	    hbin = fHistNum->GetAxis(1)->FindBin(high_value-1.e-7);
	    
	    fHistNum->GetAxis(1)->SetRange(lbin,hbin);
	    fHistNum->GetAxis(1)->SetBit(TAxis::kAxisRange);
	    fHistTotal[ip] = fHistNum->Projection(0,"E");//projection x//histogram sig+bkg
	    
	    fHistDen->GetAxis(2)->SetRange(lc,hc);
	    fHistDen->GetAxis(2)->SetBit(TAxis::kAxisRange);

	    fHistDen->GetAxis(1)->SetRange(lbin,hbin);
	    fHistDen->GetAxis(1)->SetBit(TAxis::kAxisRange);
	    fHistBkg[ip] = fHistDen->Projection(0,"E");//projectionx combine background histogram
	     
	     
	     cout <<  fHistTotal[ip]->GetBinWidth(12) << endl;

	     fHistTotal[ip]->Rebin(2);
	     fHistBkg[ip]->Rebin(2);
		 
	     
	     //sprintf(name,"histPP_pt%d",ip); //projection invaraint mass  like sign methode

	     /*
	     if(ip ==0 )
	       {
		 fHistTotal[ip]->Rebin(2);
		 fHistBkg[ip]->Rebin(2);
		 // fHistlike[ip]->Rebin(3);
	       }
		 
	     if(ip >=1 && ip<=5)
	       {
		 fHistTotal[ip]->Rebin(1);
		fHistBkg[ip]->Rebin(1);
		//fHistlike[ip]->Rebin(2);
		}
	     if(ip>5)
	       {
	       fHistTotal[ip]->Rebin(1);
	       fHistBkg[ip]->Rebin(1);
	       //fHistlike[ip]->Rebin(4);
	       }
	     */

	     // return 0;
	     fHistBkg[ip]->SetName(name);
	     fHistTotal[ip]->SetName(name);

	     Int_t lownormbin = fHistTotal[ip]->GetXaxis()->FindBin(lownorm); // change it to bin no.//findbin for bkg normalisation
	     Int_t hinormbin = fHistTotal[ip]->GetXaxis()->FindBin(hinorm);


	     Int_t lownormbin1 = fHistTotal[ip]->GetXaxis()->FindBin(lownorm1); // change it to bin no.//findbin for bkg normalisation
	     Int_t hinormbin1 = fHistTotal[ip]->GetXaxis()->FindBin(hinorm1);
	     
	     sig_integral= fHistTotal[ip]->Integral(lownormbin,hinormbin); //using intrgral method find counts for normalisation
	     bkg_integral= fHistBkg[ip]->Integral(lownormbin,hinormbin);
	     
	     
	     sig_integral1= fHistTotal[ip]->Integral(lownormbin1,hinormbin1); //using intrgral method find counts for normalisation
	     bkg_integral1= fHistBkg[ip]->Integral(lownormbin1,hinormbin1);
	     
	     // error in normalisation factor 
	     if(sig_integral <= 0 || bkg_integral <=0) continue;

	     Double_t normfactor  = sig_integral/bkg_integral;
	     Double_t normfactorerror = sqrt(normfactor/bkg_integral+normfactor*normfactor/bkg_integral);

	     TH1D *hfsig=(TH1D*)fHistTotal[ip]->Clone();
	     TH1D *hfbkg=(TH1D*)fHistBkg[ip]->Clone();   //mixed event 

	     hfbkg->Scale(sig_integral1/bkg_integral1);//normalised bkg for mixed event

	     TFile *fout_bkg = new TFile("K0sK0s_signalplusbkg_0to100_1to10GeV.root","recreate");
	     fout_bkg->cd();
	     hfbkg->Write();
	     hfsig->Write();

	     
	     /*if(ip ==0)
	       {
	       hfbkg->Scale(sig_integral1/bkg_integral1);//normalised bkg for mixed event 
	       
	       }
	       
	       else
               {
	       
	       hfbkg->Scale(sig_integral/bkg_integral);//normalised bkg for mixed event
	       }
	     */
	     hfsig->Add(hfbkg,-1);//signal after subtraction of bkg

	     TFile *fout_sig = new TFile("K0sK0s_signal_0to100_1to10GeV.root","recreate");
	     fout_sig->cd();
	     hfsig->Write();

	     pt=(Low_pt[ip]+Low_pt[ip+1])/2;
	     dpt=Low_pt[ip+1]-Low_pt[ip];
	     
	     //cout<<pt[ip]<<endl;
	     BinWidth=hfsig->GetBinWidth(15);
	     //cout<< BinWidth <<endl;

	     // fitting BW and poly2 
	     /*TF1 *fitFcn = new TF1("fitfunc",BreitWigner,fitlow[ip],fithigh[ip],6);
	       TF1 *fitFcn1 = new TF1("fitfunc1",polynomial2,fitlow[ip],fithigh[ip],3);
	       TF1 *fitFcn2 = new TF1("fitFcn2",BW,fitlow[ip],fithigh[ip],3);
	     */

	     
	     TF1 *fitFcn = new TF1("fitfunc",rBreitWigner,fitlow[ip],fithigh[ip],12);
             TF1 *fitFcn1 = new TF1("fitfunc1",polynomial2,fitlow[ip],fithigh[ip],3);
             TF1 *fitFcn2 = new TF1("fitFcn2",rBW,fitlow[ip],fithigh[ip],9);
	    
	     
	     //    TF1 *fitFcn = new TF1("fitfunc",rBreitWigner,fitlow[ip],fithigh[ip],9);
	     //TF1 *fitFcn1 = new TF1("fitfunc1",polynomial2,fitlow[ip],fithigh[ip],3);
	     // TF1 *fitFcn2 = new TF1("fitFcn2",rBW,fitlow[ip],fithigh[ip],6);


	     //	     fitFcn->SetParameter(0,10000); //yield
	     // fitFcn->SetParLimits(1,1.2,1.4); //inv. mass peak range
	     // fitFcn->SetParamet(1,1.270); //inv. mass peak range
	     //fitFcn->SetParLimits(2,0.03,0.1); //width-fixed
	     //	     fitFcn->FixParameter(2,0.10); //width-fixed


	      // fitFcn->SetParameter(0,10000); //yield
	      //fitFcn->FixParameter(1,1.270); //inv. mass peak range
	     //	     fitFcn->SetParLimits(2,0.03,0.15); //width-fixed
	     // fitFcn->FixParameter(2,0.10); //width-fixed


	    fitFcn->SetParameter(0,1000); //yield
	    fitFcn->SetParLimits(1,1.0,1.4);
	    //	    fitFcn->SetParLimits(2,0.06,0.18);
	    fitFcn->FixParameter(2,0.185);

	    //  fitFcn->SetParLimits(1,1.15,1.3); //inv. mass peak range
	     //fitFcn->SetParLimits(2,0.03,0.1); //width-fixed
	    //	    fitFcn->SetParameter(1,1.270); //inv. mass peak range
	    //  fitFcn->SetParameter(2,0.186); //width-fixed
	     
	   
	     
	     
	    fitFcn->SetParameter(3,10000); //yield

	    fitFcn->SetParLimits(4,1.46,1.56); //inv. mass peak range

	    // fitFcn->SetParLimits(5,0.03,0.16); //width-fixed/

	      fitFcn->SetParameter(5,0.086); //width-fixed
	     
	     //  fitFcn->SetParameter(7,1.525); //inv. mass peak range
	     //fitFcn->SetParameter(8,0.073); //width-fixed
	     
	    fitFcn->SetParameter(6,1000); //yield

	    fitFcn->SetParLimits(7,1.6,1.9); //inv. mass peak range

	    fitFcn->SetParLimits(8,0.08,0.16); //width-fixed

	    //	    fitFcn->SetParameter(8,0.123); //width-fixed

	    
	    fitFcn->SetParameter(9,10000000); //width-fixed
	      fitFcn->SetParameter(10,1.0); //width-fixed
	     fitFcn->SetParameter(11,0.0001); //width-fixed
	     
	    
	     /*
	     fitFcn->SetParameter(0,Para1[ip]); //yield
	     //	     fitFcn->SetParLimits(1,1.45,1.7); //inv. mass peak range
	     fitFcn->SetParLimits(1,1.45,1.66); //inv. mass peak range
	     //fitFcn->SetParameter(4,1.525); //inv. mass peak range
	     fitFcn->SetParLimits(2,0.03,0.15); //width-fixed
	     //    fitFcn->SetParameter(2,0.73); //width-fixed
	     
	     
	     fitFcn->SetParameter(3,1000); //yield
	     fitFcn->SetParLimits(4,1.6,1.8); //inv. mass peak range
	     //fitFcn->SetParameter(5,0.140); //width-fixed
	     //fitFcn->SetParLimits(5,0.08,0.16); //width-fixed
	     fitFcn->SetParLimits(5,0.08,0.171); //width-fixed
	     //fitFcn->SetParameter(5,0.16); //width-fixed

       	     
	     // fitFcn->SetParameter(6,1); //yield
	     // fitFcn->SetParameter(7,10); //width-fixed
	     // fitFcn->SetParameter(8,10); //width-fixed

	   //   fitFcn->SetParameter(9,1000); //yield
	    /// fitFcn->SetParameter(10,10); //width-fixed
	     // fitFcn->SetParameter(11,10); //width-fixed
	     */
	  
	     //	     fitFcn->SetParNames("Yield","Mass","width","B","n","C");
	     
	     // fitFcn->SetParNames("Yield","Mass","width","yieldf1710","mass1710","width1710""B","n","C");
	     //	     fitFcn->SetParNames("Y1","M1","#Gamma1","Y2","M2","Gamma2""Y3","M3","#Gamma3","Y4");
		     
	     
	     //	     fitFcn->SetParameter(3,RBPar1[ip]);
	     // fitFcn->SetParameter(4,RBPar2[ip]);
	     // fitFcn->SetParameter(5,RBPar3[ip]);


	    // fitFcn->SetParameter(3,1100);
	       // fitFcn->SetParameter(4,1.525);
	       // fitFcn->SetParameter(5,20);
	     
	   r=hfsig->Fit(fitFcn,"REBMS+"); //signal after like subtraction likesign  bkg 
	   Double_t *par = fitFcn->GetParameters();
	   Yield=fitFcn->GetParameter(0);
	   // poly0[ip]=fitFcn->GetParameter(6);        
	   
	   // poly1[ip]=fitFcn->GetParameter(12);
	   // poly2[ip]=fitFcn->GetParameter(13);
	   // poly3[ip]=fitFcn->GetParameter(14)
	   // cout << "BKGFunc" <<  poly1[ip] << "\t" << poly2[ip] << "\t" << poly3[ip] << endl;
	   // fitFcn1->SetParameter(12,poly1[ip]);
	   //fitFcn1->SetParameter(13,poly2[ip]);
	   // fitFcn1->SetParameter(14,poly3[ip]);
	   
	   fitFcn1->SetParameters(&par[9]);
	   fitFcn2->SetParameters(&par[0]);
	   
	   /*************GetChiSquare********************/
	   d=fitFcn->GetChisquare();
	   c=fitFcn->GetNDF();
	   chi=d/c;
	   Double_t Yield1320, Yield1320error, Yield1525, Yield1525error,Yield1710,Yield1710error, mass1320, error_mass1320, width1320,error_width1320,mass1525, mass1710, error_mass1525, error_mass1710, width1525, width1710,error_width1525,error_width1710;

	   
	   /**************mass spectra************/
           Yield1320= fitFcn->GetParameter(1);
	   Yield1320error= fitFcn->GetParError(1);
	   
	   mass1320=fitFcn->GetParameter(2); 
	   error_mass1320= fitFcn->GetParError(2);
	   
	   width1320=fitFcn->GetParameter(3); 
	   error_width1320= fitFcn->GetParError(3);


	   Yield1525= fitFcn->GetParameter(4);
	   Yield1525error= fitFcn->GetParError(4);
	   
	   mass1525=fitFcn->GetParameter(5); 
	   error_mass1525= fitFcn->GetParError(5);
	   
	   width1525=fitFcn->GetParameter(6); 
	   error_width1525= fitFcn->GetParError(6);

	   cout << "mass upper " << mass1525+ width1525 << "\t--mass lower----"<<  mass1525- width1525 << endl;
	   
	   Yield1710= fitFcn->GetParameter(7);
	   Yield1710error= fitFcn->GetParError(7);
	   
	   mass1710=fitFcn->GetParameter(8); 
	   error_mass1710= fitFcn->GetParError(8);
	   
	   width1710=fitFcn->GetParameter(9); 
	   error_width1710= fitFcn->GetParError(9);
	   
	   hmass1525->SetBinContent(ip+1,mass1525);
	   hmass1525->SetBinError(ip+1,error_mass1525);
	   
	   hwidth1525->SetBinContent(ip+1,width1525);
	   hwidth1525->SetBinError(ip+1,error_width1525);
	   
	   
	   hmass1710->SetBinContent(ip+1,mass1710);
	   hmass1710->SetBinError(ip+1,error_mass1710);
	   
	   hwidth1710->SetBinContent(ip+1,width1710);
	   hwidth1710->SetBinError(ip+1,error_width1710);

	   //hmass1525->Draw();
	   //return 0;
	   cout << "Mass1525" << hmass1525->GetBinContent(1) << endl;
	    
	    
	    /***********transversemass************/
	    mT = TMath::Sqrt(pt*pt + mass*mass);
	    mT_error = TMath::Sqrt(pt*pt*dpt*dpt+mass*mass*error_mass*error_mass)/mT;
	    htransmass->SetBinContent(ip+1,mT);
	    htransmass->SetBinError(ip+1,mT_error);
	    
	    /*************ERROR BIN COUNTING METHOD CALCULATION************/
	    TF1 *fitFcn2_plusm = new TF1("fitFcn2_plusm",BW,fitlow[ip],fithigh[ip],3);
	    TF1 *fitFcn2_minusm = new TF1("fitFcn2_plusm",BW,fitlow[ip],fithigh[ip],3); 
	    
	    fitFcn2_plusm->FixParameter(0, Yield);
	    fitFcn2_plusm->FixParameter(1, width +error_width);
	    fitFcn2_plusm->FixParameter(2, mass + error_mass);
	    fitFcn2_minusm->FixParameter(0, Yield);
	    fitFcn2_minusm->FixParameter(1,width -error_width);
	    fitFcn2_minusm->FixParameter(2, mass - error_mass);
	    
	    /********different method for Yield coounting check******************/
	   
	    bmin=hfsig->GetXaxis()->FindBin(lowsigmarange);
	    bmax=hfsig->GetXaxis()->FindBin(highsigmarange);
	    Yield_bincount_hist=hfsig->IntegralAndError (bmin,bmax,hBCError_1);
	    bkgvalue = fitFcn1->Integral(hfsig->GetBinLowEdge(bmin),hfsig->GetBinLowEdge(bmax+1));
	    Integral_BW_withsigma =fitFcn2->Integral(hfsig->GetBinLowEdge(bmin),hfsig->GetBinLowEdge(bmax+1));
	    fYield_BinCount = Yield_bincount_hist - (bkgvalue/BinWidth) ;
	    YieldIntegral_BW = fitFcn2->Integral(1.25,2)/BinWidth;
	    Yfraction_cBW =(Integral_BW_withsigma/YieldIntegral_BW);

	    sum_tail_correction=(fitFcn2->Integral(1.25,hfsig->GetBinLowEdge(bmin)) + fitFcn2->Integral(hfsig->GetBinLowEdge(bmax+1),5))/BinWidth;
	    Total_Ybincounting = (sum_tail_correction + fYield_BinCount)/(Event*dpt*dy*0.33*2);
	    Tail_correction_plusm =(fitFcn2_plusm->Integral(1.25,hfsig->GetBinLowEdge(bmin)) + (fitFcn2_plusm->Integral(hfsig->GetBinLowEdge(bmax+1),5)))/BinWidth; 
	    Tail_correction_minusm =((fitFcn2_minusm->Integral(1.25,hfsig->GetBinLowEdge(bmin))+ fitFcn2_minusm->Integral(hfsig->GetBinLowEdge(bmax+1),5))/BinWidth);
	    Error_2 =sum_tail_correction -Tail_correction_plusm;
	    Final_pro_error = TMath::Sqrt(Error_2*Error_2 +hBCError_1*hBCError_1)/(Event*dpt*dy*0.33*2);  

	    Yield_value_par=fitFcn->GetParameter(0)/(Event*dpt*dy);
	    yield_error_par = fitFcn->GetParError(0)/(Event*dpt*dy);

	    Yield_value_par1710=fitFcn->GetParameter(3)/(Event*dpt*dy);
	    yield_error_par1710 = fitFcn->GetParError(3)/(Event*dpt*dy); 
	    
	    // cout<< sum_tail_correction << "\t" <<Tail_correction_plusm  <<"\t" <<  Tail_correction_minusm << "\t" << hBCError_1 << "\t" <<  Error_2 << endl;
	    
	    /*****************uncorrected Yield**********************/
	    
	    hYbincount->SetBinContent(ip+1,Total_Ybincounting);
	    hYbincount->SetBinError(ip+1,Final_pro_error);

	    hYieldpar1525->SetBinContent(ip+1,Yield_value_par);
	    hYieldpar1525->SetBinError(ip+1, yield_error_par);

	    hYieldpar1710->SetBinContent(ip+1,Yield_value_par1710);
	    hYieldpar1710->SetBinError(ip+1, yield_error_par1710);

	    
	    //fractional stat error
	    hFrac_stat_error->SetBinContent(ip+1,Final_pro_error/Total_Ybincounting);
	    
	    {
	      significance_num = (fitFcn2->Integral(1.25,2.0))/(BinWidth);
	      significance_den = TMath::Sqrt(fHistTotal[ip]->Integral(bmin,bmax));
	      // ratio = significance_num/significance_den;

	      //  cout<<ip<<"\t" <<  BinWidth << "\t" << significance_num << "\t" <<  significance_den <<"\t"<< ratio << (Final_pro_error/Total_Ybincounting) <<  endl;

    //out<<ip<<"\t" <<  BinWidth << "\t" << significance_num << "\t" <<  significance_den <<"\t"<< ratio << (Final_pro_error/Total_Ybincounting) <<  endl;
	      hsgnfcance->SetBinContent(ip+1,ratio);
	    }
	   
  
	    // significance for all three resonances

	    // one set of resonances : a2(1320)/f2(1320) , keep accorrding mass peak 1320. Mass = 1318 +- 0.6 , width = 107 +- 5, 1 sigma = 1318 +- 107, 1318 +- 2*107, 1318 +- 3.0*107

	    // range for sigma calculation 1 sigma 
	    //  Double_t lr1 = 1.104;
	    //Double_t hr1 = 1.532;

	     Double_t lr1 = 1.104;
	     Double_t hr1 = 1.532;

	    //  Double_t lr1 = mass1320-width1320;
	    // Double_t hr1 = mass1320+width1320;
	    
	    // bmin1=hfsig->GetXaxis()->FindBin(1.1);
	    // bmax1=hfsig->GetXaxis()->FindBin(1.53);
	    Double_t bmin1, bmax1;
	    bmin1=hfsig->GetXaxis()->FindBin(lr1);
	    bmax1=hfsig->GetXaxis()->FindBin(hr1);
	    
	    Double_t splusb_r11 = fHistTotal[ip]->Integral(bmin1,bmax1);
	    Double_t splusb_r1 = TMath::Sqrt(fHistTotal[ip]->Integral(bmin1,bmax1));

	   Double_t signal_r1 =  fitFcn2->Integral(lr1,hr1)/(BinWidth);
	    // Double_t signal =  fitFcn2->Integral(1.2,1.4)/(BinWidth);
	   Double_t signi_r1 = (signal_r1/splusb_r1);

	   cout << "bin--low" << bmin1 << "bin--high" << bmax1 << "\t" << "signal+bkg" << splusb_r11 << "\t----sqrt(s+b)" << splusb_r1  << "\t--signal" <<  signal_r1 << "\t"  << signi_r1 << endl;

	 

	   // one set of resonances : 1525 , keep accorrding mass peak 1320. Mass = 1525 +- 2.5 , width = 86 +- 5, 1 sigma = 1525 +-86, 1525 +- 2*86

	   // Double_t lr2 = mass1525-width1525;
	   // Double_t hr2 = mass1525+width1525;

	   Double_t lr2 = 1.439;
	   Double_t hr2 = 1.611;

	   // bmin1=hfsig->GetXaxis()->FindBin(1.1);
	    //bmax1=hfsig->GetXaxis()->FindBin(1.53);
	    Double_t bmin2, bmax2;
	    bmin2=hfsig->GetXaxis()->FindBin(lr2);
	    bmax2=hfsig->GetXaxis()->FindBin(hr2);
	    
	    Double_t splusb_r22 = fHistTotal[ip]->Integral(bmin2,bmax2);
	    Double_t splusb_r2 = TMath::Sqrt(fHistTotal[ip]->Integral(bmin2,bmax2));


	    Double_t signal_r2 =  fitFcn2->Integral(lr2,hr2)/(BinWidth);
	   // Double_t signal =  fitFcn2->Integral(1.2,1.4)/(BinWidth);
	    Double_t signi_r2 = (signal_r2/splusb_r2);

	   cout << "bin--low" << bmin2 << "bin--high" << bmax2 << "\t" << "signal+bkg" << splusb_r22 << "\t----sqrt(s+b)" << splusb_r2  << "\t--signal" <<  signal_r2 << "\t"  << signi_r2 << endl;

	   // one set of resonances : 1710 , keep accorrding mass peak 1320. Mass = 1704 +- 12 , width = 123 +-18, 1 sigma = 1710 +-86, 1525 +- 2*86
	   

	    Double_t lr3 = 1.587;
	    Double_t hr3 = 1.833;
	   
	   /// Double_t lr3 = mass1710-width1710;
	   // Double_t hr3 = mass1710+width1710;

	   // bmin1=hfsig->GetXaxis()->FindBin(1.1);
	    //bmax1=hfsig->GetXaxis()->FindBin(1.53);
	    Double_t bmin3, bmax3;
	    bmin3=hfsig->GetXaxis()->FindBin(lr3);
	    bmax3=hfsig->GetXaxis()->FindBin(hr3);
	    
	    Double_t splusb_r33 = fHistTotal[ip]->Integral(bmin3,bmax3);
	    Double_t splusb_r3 = TMath::Sqrt(fHistTotal[ip]->Integral(bmin3,bmax3));


	    Double_t signal_r3 =  fitFcn2->Integral(lr3,hr3)/(BinWidth);
	    //Double_t signal =  fitFcn2->Integral(1.2,1.4)/(BinWidth);
	    Double_t signi_r3 = (signal_r3/splusb_r3);

	   cout << "bin--low" << bmin3 << "bin--high" << bmax3 << "\t" << "signal+bkg" << splusb_r33 << "\t----sqrt(s+b)" << splusb_r3  << "\t--signal" <<  signal_r3 << "\t"  << signi_r3 << endl;
	       // one set of resonances : a2(1320)/f2(1320) , keep accorrding mass peak 1320. Mass = 1525 +- 2.5 , width = 86 +- 5


	     

	       // one set of resonances : a2(1320)/f2(1320) , keep accorrding mass peak 1320. Mass = 1704 +- 12 , width = 123 +- 18
	    
	       
	    hChiSquare->SetBinContent(ip+1,chi);
	    /*****setting plot parametres ***********/
	    
	    SetHistoStyle(hfsig,1,20,0.6,0.05,0.05,1.0,1.1);
	    SetHistoStyle(hfbkg,2,24,0.5,0.05,0.05,1.0,1.1);
	    SetHistoStyle(fHistTotal[ip],1,24,0.5,0.05,0.05,1.0,1.1);
	    hfsig->GetXaxis()->SetTitle("m_{K^{0}_{s}K^{0}_{s}} (GeV/c^{2})");
	    hfsig->GetYaxis()->SetTitle(Form("Counts/(%2.0f MeV/c^{2})",1000*hfsig->GetBinWidth(10)));

	    //hfsig->GetXaxis()->SetRangeUser(1.1,2.0);
	    hfbkg->GetXaxis()->SetRangeUser(0.9,2.5);
	    hfbkg->GetXaxis()->SetTitle("m_{inv} (GeV/c^{2})");
	    hfbkg->GetYaxis()->SetTitle(Form("Counts/(%2.0f MeV/c^{2})",1000*hfbkg->GetBinWidth(10)));
	    
	    fitFcn1->SetLineColor(4);
	    fitFcn1->SetLineStyle(2);
	    fitFcn1->SetLineWidth(2);
	    fitFcn2->SetLineColor(6);
	    fitFcn2->SetLineStyle(2);
	    fitFcn2->SetLineWidth(2);


	    //fitFcn3->SetLineColor(1);
	    //fitFcn3->SetLineStyle(2);
	    //fitFcn3->SetLineWidth(2);


	    //fitFcn4->SetLineColor(3);
	    //fitFcn4->SetLineStyle(2);
	    //fitFcn4->SetLineWidth(2);
	    
	    
	    // fitFcn5->SetLineColor(7);
	    // fitFcn5->SetLineStyle(2);
	    // fitFcn5->SetLineWidth(2);
	    
	    //fitFcn6->SetLineColor(41);
	    //fitFcn6->SetLineStyle(2);
	    // fitFcn6->SetLineWidth(2);

	    
	    /***************************************/
	    //plot
	    cinv[ip]->cd();
	    hfsig->GetXaxis()->SetRangeUser(1.1,2.0);
	    hfsig->Draw("e");
	    fitFcn->Draw("same");
	    fitFcn1->Draw("same");
	    fitFcn2->SetFillColor(2);
	    fitFcn2->Draw("same");
           
	    
	    //	    fitFcn3->Draw("same");
	    // fitFcn4->Draw("same");
	    // fitFcn5->Draw("same");
	    // fitFcn6->Draw("same");
	   sprintf(name,"%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip],Low_pt[ip+1]);
	    if(ip<1)
	      {
		c2->cd(ip+1);
		hfsig->Draw("e");
		hfsig->GetXaxis()->SetRangeUser(1.1,2.0);
		TF1 *ftotal = (TF1*)fitFcn->Clone();
		TF1 *fBW = (TF1*)fitFcn2->Clone();
		TF1 *fpoly = (TF1*)fitFcn1->Clone();
		ftotal->Draw("same");
		fBW->Draw("same");
		fpoly->Draw("same");
		
		TLatex *ltx = new TLatex(0.35,0.93,name);
		ltx->SetNDC();
		ltx->SetTextFont(22);
		ltx->SetTextSize(0.07187006);
		ltx->Draw();
		
		TLegend *leg = new TLegend(0.68,0.7,0.89,0.83,NULL,"brNDC");
		leg->SetBorderSize(0.3);
		leg->SetTextFont(22);
		leg->SetTextSize(0.05);
		leg->SetLineColor(1);
		leg->SetLineStyle(1);
		leg->SetLineWidth(1);
		leg->SetFillColor(0);
		leg->AddEntry(fitFcn,"3rBW + Res.Bkg","l");
		leg->AddEntry(fBW,"3rBW","l");
		leg->AddEntry(fpoly,"Res.Bkg","l");
	  	leg->Draw();
		
	      }
          else
	    {
	      c22->cd(ip-11);
 	     hfsig->GetXaxis()->SetRangeUser(1.1,2.0);
	      hfsig->Draw("e");

	      TF1 *ftotal = (TF1*)fitFcn->Clone();
	      TF1 *fBW = (TF1*)fitFcn2->Clone();
	      TF1 *fpoly = (TF1*)fitFcn1->Clone();
	      ftotal->Draw("same");
	      fBW->Draw("same");
	      fpoly->Draw("same");
	      TLatex *ltx = new TLatex(0.35,0.93,name);
	      ltx->SetNDC();
	      ltx->SetTextFont(22);
	      ltx->SetTextSize(0.07187006);
	      ltx->Draw();
	      
	      TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
	      leg->SetBorderSize(0.3);
	      leg->SetTextFont(22);
	      leg->SetTextSize(0.05);
	      leg->SetLineColor(1);
	      leg->SetLineStyle(1);
	      leg->SetLineWidth(1);
	      leg->SetFillColor(0);
	      leg->AddEntry(fitFcn,"4rBW+Res.Bkg","l");
	      leg->AddEntry(fBW,"rBW","l");
	      leg->AddEntry(fpoly,"Res.Bkg","l");
	      leg->Draw();
	    }
       
       //signal extraction before bkg subtraction :
      
       
       if(ip<12)
	 {
	   c1->cd(ip+1);
	   hfbkg->Draw("e");
	   fHistTotal[ip]->Draw("samee");
	   TLatex *ltx = new TLatex(0.35,0.93,name);
	   ltx->SetNDC();
	   ltx->SetTextFont(22);
	   ltx->SetTextSize(0.07187006);
	   ltx->Draw();
	   TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
	   leg->SetBorderSize(0.3);
	   leg->SetTextFont(22);
	   leg->SetTextSize(0.05);
	   leg->SetLineColor(1);
	   leg->SetLineStyle(1);
	   leg->SetLineWidth(1);
	   leg->SetFillColor(0);
	   leg->AddEntry(fHistTotal[ip],"Sig+bkg","p");
	   leg->AddEntry(hfbkg,"Normalised mixed bkg","p");
	   leg->Draw();
	  }
       else
	 {
	   c11->cd(ip-11);
	   fHistTotal[ip]->GetXaxis()->SetRangeUser(0.9,2.5);
	   fHistTotal[ip]->Draw("pe");
	   hfbkg->Draw("pe same");
	   TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
	   leg->SetBorderSize(0.3);
	   leg->SetTextFont(22);
	   leg->SetTextSize(0.05);
	   leg->SetLineColor(1);
	   leg->SetLineStyle(1);
	   leg->SetLineWidth(1);
	   leg->SetFillColor(0);
	   leg->AddEntry(fHistTotal[ip],"Sig+bkg","p");
	   leg->AddEntry(hfbkg,"Normalised mixed bkg","p");
	   leg->Draw();
	   TLatex *ltx = new TLatex(0.35,0.93,name);
	   ltx->SetNDC();
	   ltx->SetTextFont(22);
	   ltx->SetTextSize(0.07187006);
	   ltx->Draw();
	 }
       cSigbkg[ip]->cd();
       fHistTotal[ip]->GetXaxis()->SetRangeUser(0.9,2.5);
       fHistTotal[ip]->Draw("e");
       hfbkg->Draw("same pe");              
       //covariant matrix for yielderror  calculation
       TMatrixDSym mat = r->GetCovarianceMatrix();
       //mat.Print();
       TMatrixDSym mat1;
       TMatrixDSym mat2;
       mat.GetSub(0,2,0,2,mat1);
       mat.GetSub(3,5,3,5,mat2);
       Double_t * b   = mat1.GetMatrixArray();
       Double_t * a   = mat2.GetMatrixArray();  	   
       // wRawYield[ip] = fitFcn2->Integral(0.66,5.0)/(BinWidth*2*Event*dpt*dy*0.66);
       wRawYield[ip] = fitFcn2->Integral(1.25,5.0)/(BinWidth*2*Event*dpt*dy*0.33); 
       wErRawYield[ip] = fitFcn2->IntegralError(1.25,5.0,&par[0],b)/(BinWidth*Event*dpt*dy*0.33*2);
       // cout <<  wRawYield[ip] << "\t" << Total_Ybincounting <<"\t" << Yield_value_par << "\t" << wErRawYield[ip]  << "\t" << Final_pro_error <<"\t" << yield_error_par<<endl;
       //  out <<  wRawYield[ip] << "\t" << Total_Ybincounting <<"\t" << Yield_value_par << "\t" << wErRawYield[ip]  << "\t" << Final_pro_error <<"\t" << yield_error_par<<endl;  
	 }
       
       }
   // return 0;
  /*************hsignificance plot************/
   TCanvas*cs = new TCanvas("cs","",10,10,600,600);
   cs->cd();
   SetHistoStyle(hsgnfcance,2,8,1.5,0.05,0.05,1.0,1.1);
   hsgnfcance->GetYaxis()->SetTitle("S/#sqrt{S+B}");
   hsgnfcance->GetXaxis()->SetTitle("p_{T} GeV/c");
   hsgnfcance->Draw("p");
   // hsgnfcance->Print();
   TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
   leg->SetBorderSize(0.3);
   leg->SetTextFont(22);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->AddEntry(hsgnfcance,"Significance","p");
   leg->Draw("");
   
   //hsgnfcance-->>GetYaxis()->SetTitle("Significance(S/#Sqrt(S+B))");

 
  TCanvas *cmass1525 = new TCanvas("cmass1525","cmass1525", 10,10,600,600);
  cmass1525->cd();
  SetHistoStyle(hmass1525,2,20,1.5,0.05,0.05,1.0,1.1);
  TLine *line = new TLine(1.0,1.525,10,1.525);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  hmass1525->GetYaxis()->SetTitleFont(42);   
  hmass1525->GetYaxis()->SetTitle("Mass(GeV/c^{2})");
  hmass1525->GetXaxis()->SetTitle("p_{T} GeV/c");
  hmass1525->Draw("p");
  line->Draw("p same");
  TCanvas *cmass1710 = new TCanvas("cmass1710","cmass1710", 10,10,600,600);
  cmass1710->cd();
  SetHistoStyle(hmass1710,2,20,1.5,0.05,0.05,1.0,1.1);
  TLine *line = new TLine(1.0,1.710,10,1.710);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  hmass1710->GetYaxis()->SetTitleFont(42);   
  hmass1710->GetYaxis()->SetTitle("Mass(GeV/c^{2})");
  hmass1710->GetXaxis()->SetTitle("p_{T} GeV/c");
  hmass1710->GetYaxis()->SetRangeUser(1.6,1.8);
  hmass1710->Draw("p");
  line->Draw("p same");
  TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hmass,"Mass Peak","p");
  leg->AddEntry(line,"PDG","l");
  leg->Draw("");

  TCanvas *cwidth1525 = new TCanvas("cwidth1525","", 10,10,600,600);
  cwidth1525->cd();
  SetHistoStyle(hwidth1525,2,20,1.5,0.05,0.05,1.0,1.1);
  TLine *line = new TLine(1,0.073,10,0.073);
  line->SetLineColor(kRed);
  line->SetLineStyle(2); 
  
  hwidth1525->GetYaxis()->SetTitleFont(42);   
  hwidth1525->GetYaxis()->SetTitle("Width(GeV/c^{2})");
  hwidth1525->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hwidth1525->GetYaxis()->SetRangeUser(0.0,0.3);
  hwidth1525->Draw("p");
  line->Draw("psame");

  TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hwidth1525,"f2(1525)","p");
  leg->AddEntry(line,"PDG","l");
  leg->Draw("");
  
  TCanvas *cwidth1710 = new TCanvas("cwidth1710","", 10,10,600,600);
  cwidth1710->cd();
  SetHistoStyle(hwidth1710,2,20,1.5,0.05,0.05,1.0,1.1);
  TLine *line = new TLine(1,0.139,10,0.139);
  line->SetLineColor(kRed);
  line->SetLineStyle(2); 
  
  hwidth1710->GetYaxis()->SetTitleFont(42);   
  hwidth1710->GetYaxis()->SetTitle("Width(GeV/c^{2})");
  hwidth1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hwidth1710->GetYaxis()->SetRangeUser(0.0,0.3);


  TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hwidth1710,"f0(1710)","p");
  leg->AddEntry(line,"PDG","l");
  leg->Draw("");

  hwidth1710->Draw("p");
  line->Draw("psame");
  
  
  TCanvas *cmT = new TCanvas("cmT","", 10,10,600,600);
    SetHistoStyle(htransmass,2,20,1.5,0.05,0.05,1.0,1.1);
    htransmass->GetYaxis()->SetTitleFont(42); 
    htransmass->GetYaxis()->SetTitle("m_{T} GeV/c)");
    htransmass->GetXaxis()->SetTitle("p_{T} GeV/c");
    htransmass->Draw();

 TCanvas *cptspectra1525 = new TCanvas("cptspectra1525","", 10,10,600,600);
  cptspectra1525->cd();
  SetHistoStyle(hYieldpar1525,2,20,1.5,0.05,0.05,1.0,1.1);
  hYieldpar1525->GetYaxis()->SetTitleFont(42);   
  hYieldpar1525->GetYaxis()->SetTitle("Yield");
  hYieldpar1525->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hYieldpar1525->Draw("p");
  
  TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hYieldpar1525,"f2(1525)","p");
  leg->Draw("");

 TCanvas *cptspectra1710 = new TCanvas("cptspectra1710","", 10,10,600,600);
  cptspectra1710->cd();
  SetHistoStyle(hYieldpar1710,2,20,1.5,0.05,0.05,1.0,1.1);
  hYieldpar1710->GetYaxis()->SetTitleFont(42);   
  hYieldpar1710->GetYaxis()->SetTitle("Yield");
  hYieldpar1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hYieldpar1710->Draw("p");
  
  TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.05);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(hYieldpar1710,"f0(1710)","p");
  leg->Draw("");
  
   /*
      TCanvas *ctranM = new TCanvas("ctranM","", 10,10,600,600);
      ctranM->cd();
      SetHistoStyle(htransmass,2,20,1.5,0.05,0.05,1.0,1.1);
      hm->GetYaxis()->SetTitle("m GeV/c)");
      hm->GetXaxis()->SetTitle("m_{T} GeV/c");
   hm->Draw();
    */

   /*TCanvas *comparision = new TCanvas("comparision","", 10,10,600,600);
     comparision->cd();
     hmass->Draw("p");
     // hm->Draw("same");
     // SetHistoStyle(hm,4,20,1.5,0.05,0.05,1.0,1.1);
     //hm->GetYaxis()->SetTitleFont(42);
     //hm->GetYaxis()->SetTitle("m GeV/c)");
     //hm->GetXaxis()->SetTitle("m_{T}  or p_{T} GeV/c");
     
     TLegend *leg = new TLegend(0.580554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
     leg->SetBorderSize(0.3);
     leg->SetTextFont(22);
     leg->SetTextSize(0.05);
     leg->SetLineColor(1);
     leg->SetLineStyle(1);
     leg->SetLineWidth(1);
     leg->SetFillColor(0);
     leg->AddEntry(hmass,"m vs p_{T}","p");
     // leg->AddEntry(hm,"m vs m_{T}","p");
     leg->Draw("");
     TCanvas *chisquare = new TCanvas("chisquare","", 10,10,600,600);
     chisquare->cd();
     SetHistoStyle(hChiSquare,1,20,1.5,0.05,0.05,1.0,1.1);
     hChiSquare->Draw("p");
     hChiSquare->GetYaxis()->SetTitleFont(42);
     hChiSquare->GetYaxis()->SetTitle("#chi^{2}/NDf");
     hChiSquare->GetXaxis()->SetTitle("p_{T} GeV/c");
     TLegend *leg = new TLegend(0.580554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
     leg->SetBorderSize(0);
     leg->AddEntry( hChiSquare,"#chi^{2}/NDf vs p_{T}","p");
     leg->Draw("");
     TH1F *hintegral_yield = new TH1F("hintegral_yield","",Npt,Low_pt);
     for(Int_t i=0;i<Npt;i++)
       {
	 hintegral_yield->SetBinContent(i+1,wRawYield[i]);
	 hintegral_yield->SetBinError(i+1,wErRawYield[i]);
       }
     
     TCanvas *g = new TCanvas("g","",10,10,600,600);
     g->cd();
     g->SetLogy();
     SetHistoStyle(hintegral_yield,6,20,1.0,0.05,0.05,1.0,1.1);
     // hintegral_yield->SetLogy(1);
     hintegral_yield->Draw("e");
     hintegral_yield->GetXaxis()->SetTitle("p_{T}in GeV/c");
     hintegral_yield->GetXaxis()->SetLabelFont(42);
     hintegral_yield->GetYaxis()->SetTitle("(1/Nevt) (d^{2}N/dp_{T}dy)");
     TCanvas *frac_stat = new TCanvas("frac_stat","",600,600);
     frac_stat ->cd();
     hFrac_stat_error->Draw();
     hFrac_stat_error->GetXaxis()->SetTitle("p_{T} GeV/c");
     hFrac_stat_error->GetYaxis()->SetTitle("Fractional Stat. Uncertainty");
     
     
    
   TCanvas *yield_comparision = new TCanvas("yield_comparision","", 10,10,700,800);
   yield_comparision->cd();
   TPad *pad1 = new TPad("pad1","pad1",0,0.2999,0.9999,0.99999);
   pad1->SetBottomMargin(0);
   pad1->SetLeftMargin(0.2);
   pad1->Draw();
   pad1->cd();
   hintegral_yield->GetYaxis()->SetTitle("(1/Nevt) (d^{2}N/dp_{T}dy)");
   hintegral_yield->GetXaxis()->SetTitle("p_{T}in GeV/c");
   
   hYbincount->SetMarkerStyle(20);
   hYieldpar->SetMarkerStyle(24);
   hintegral_yield->SetMarkerStyle(25);
   hYieldpar->SetMarkerColor(2);
   hYbincount->SetMarkerColor(4);
   hintegral_yield->SetMarkerColor(6);
   hYieldpar->SetMarkerSize(1.2);
   hintegral_yield->SetMarkerSize(1.2);
   hintegral_yield->GetYaxis()->SetTitleFont(42);
   
  hintegral_yield->GetYaxis()->SetTitleOffset(1.3);
  hintegral_yield->GetYaxis()->SetTitleSize(0.04);
  hintegral_yield->GetXaxis()->SetTitleSize(0.04);
  
  hintegral_yield->Draw("e");
  hYbincount->Draw("same");
  //hYieldpar->Draw("same");
  pad1->SetLogy();
  
  TLegend *leg = new TLegend(0.680554,0.7312735,0.892902,0.8338954,NULL,"brNDC");
  leg->SetBorderSize(0.3);
  leg->SetTextFont(22);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry( hintegral_yield,"Integral method","plf");
  leg->AddEntry( hYbincount,"BinCount method ","plf");
  //leg->AddEntry( hYieldpar,"Parmeter method","plf");
  leg->Draw("");
  yield_comparision->cd();
  TPad *bottom_pad = new TPad("bottom_pad","bottom_pad",0.0,0.0,0.9999,0.2999);
  bottom_pad->SetBottomMargin(0.1);
  bottom_pad->SetLeftMargin(0.2);
  bottom_pad->Draw();
  bottom_pad->cd();
 TH1F *h3 = (TH1F*)hintegral_yield ->Clone();
 TH1F *h4 = (TH1F*)hYbincount->Clone();
 h4->Divide(hintegral_yield);
 // h4->Divide( hYieldpar);
  h4->SetMinimum(0.8999);
  h4->SetMaximum(1.199);
  h4->SetTitle("");
  h4->SetMarkerStyle(hintegral_yield ->GetMarkerStyle());
  h4->SetMarkerSize(hintegral_yield ->GetMarkerSize());
  h4->SetMarkerColor(hintegral_yield ->GetMarkerColor());
  h4->GetXaxis()->SetLabelFont(44); //font in pixels                                                                                          
  h4->GetXaxis()->SetLabelSize(24); //in pixels                                                                                               
  h4->GetYaxis()->SetLabelFont(44); //font in pixels
  h4->GetYaxis()->SetLabelSize(24); //in pixels                                                                                                   
  h4->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h4->GetYaxis()->SetTitle("Yield Bincount/Yield integral");
  h4->GetXaxis()->CenterTitle(true);
  h4->GetYaxis()->CenterTitle(true);
  h4->GetXaxis()->SetTitleFont(42);
  h4->GetYaxis()->SetTitleFont(42);
  h4->GetXaxis()->SetTitleSize(0.15);
  h4->GetYaxis()->SetTitleSize(0.15);
  h4->GetXaxis()->SetTitleOffset(1.0);
  h4->GetYaxis()->SetTitleOffset(0.517);
  h4->GetXaxis()->SetTickLength(0.05);
  h4->Draw("same");
   */
   
  /*
  TFile *fout=new TFile("Rawyield_f1525_pp13TeV_pol2.root","RECREATE");
  fout->cd();
  c2->Write();
  c22->Write();
  c1->Write();
  c11->Write();
  g->Write();
  yield_comparision->Write();
  hintegral_yield ->Write();
  hYbincount->Write();
  hYieldpar->Write();
  hChiSquare->Write();
  hmass->Write();
  hwidth->Write();
  hsgnfcance->Write();
  */
  
 }


void SetHistoStyle(TH1 *h, Int_t mcolor, Int_t mstyle, Float_t msize,Float_t Tsizex,Float_t Tsizey,Float_t Offsetx, Float_t Offsety)
{
  h->SetMarkerColor(mcolor);
  h->SetMarkerStyle(mstyle);
  h->SetMarkerSize(msize);
  h->GetXaxis()->SetTitleSize(Tsizex);
  h->GetYaxis()->SetTitleSize(Tsizey);
  h->GetXaxis()->SetTitleOffset(Offsetx);
  h->GetYaxis()->SetTitleOffset(Offsety);
  
}


