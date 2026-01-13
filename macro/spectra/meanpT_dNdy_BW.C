
// #include </home/prottay/alice/AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C>
#include "YieldMean.C"

Double_t FuncLavy(Double_t *x, Double_t *par)
{

  Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
  return (p);

  // par[0]=n,par[1]=dndy,par[2]=mass and par[3]=Temp;
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

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;
Double_t
BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{
  /*                                                                                                                         x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
  */

  Double_t r = x[0];
  Double_t mt = p[0];
  Double_t pt = p[1];
  Double_t beta_max = p[2];
  Double_t temp_1 = 1. / p[3];
  Double_t n = p[4];

  Double_t beta = beta_max * TMath::Power(r, n);
  if (beta > 0.9999999999999999)
    beta = 0.9999999999999999;
  Double_t rho = TMath::ATanH(beta);
  Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
  if (argI0 > 700.)
    argI0 = 700.;
  Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
  return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
}

Double_t Boltzmann_Func(const Double_t *x, const Double_t *p)

{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm * mt * TMath::Exp(-mt / T) * TMath::Exp(mass / T) * (1.0 / (2 * T * T * T + 2 * T * T * mass + T * mass * mass));
}

Double_t Bose_Fuct(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm * (1 / (TMath::Exp(mt / T) - 1)) * (TMath::Exp(mass / T) - 1);
}

Double_t Exp_Fuct(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t m = p[0];
  //  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm * (TMath::Exp(-pt / T)) * TMath::Exp(m / T) * (1.0 / (T * T + T * m));
}

Double_t BGBlastWave_Func(const Double_t *x, const Double_t *p)
{

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t norm = p[4];

  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0.0, 5.0, 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  // Double_t integral = fBGBlastWave_Integrand->Integral(0.0, 5.0, (Double_t *)0, 1.e4);
  Double_t integral = fBGBlastWave_Integrand->Integral(0.0, 5.0, 1.e4);
  return norm * pt * integral;
}

void meanpT_dNdy_BW()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat(1111);

  // const Int_t Npt=17;
  const Int_t Npt = 21;
  Double_t pT, dpT, meanpT_m1_num, meanpT_m1_den, mean_pT_value, mean_pT_value1;
  Double_t sum_dndy = 0.0;
  Double_t value, dpt = 0.2, dy = 1.0;
  // Double_t Low_pt[Npt+1] ={0.4,0.8,1.2,1.6,2.2,2.6,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,10.0,12.0,16.0,20.0};
  Double_t Low_pt[Npt + 1] = {0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7, 8.0};

  TFile *f1 = new TFile("corrected_oo_pass2_final_cor1.root", "READ");
  TFile *f2 = new TFile("corrected_oo_pass2_final_cor1.root", "READ");
  TH1D *h1 = (TH1D *)f1->Get("hYbincount");
  TH1D *h2 = (TH1D *)f2->Get("hYbincount");

  TCanvas *c = new TCanvas("c", "c", 10, 10, 600, 600);
  c->cd();
  c->SetLeftMargin(0.2);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetFrameBorderMode(0);
  c->SetFillColor(0);

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

  h1->SetMarkerColor(4);
  h1->SetMarkerStyle(20);

  h1->GetYaxis()->SetTitle("(1/N)d^{2}N/(dydp_{T}) (GeV/c)^{-1}");
  h1->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  TLegend *lp2 = DrawLegend(0.6, 0.53, 0.8, 0.65);
  lp2->SetTextFont(42);
  lp2->SetTextSize(0.03);
  lp2->SetShadowColor(1);
  // lp2->AddEntry(h1,"K^{*#pm}","");
  lp2->Draw();

  TF1 *fitFcn = new TF1("fitFcn", BGBlastWave_Func, 0.0, 6, 5);

  fitFcn->FixParameter(0, 1.0195); // mass
  // fitFcn->SetParameter(1,0.78);//Beta max
  fitFcn->SetParLimits(1, 0.1, 1);    // Beta max   //0-3.5 arb par
  fitFcn->SetParameter(2, 0.15);      // T
  fitFcn->SetParLimits(2, 0.08, 0.9); // T
  fitFcn->FixParameter(3, 0.95);      // n
  fitFcn->SetParLimits(3, 0.1, 5);    // n
  fitFcn->SetParameter(4, 1.e6);      // norm

  fitFcn->SetParNames("mass", "beta_max", "T", "n", "norm");

  /*************meanpT*****************byresonance*******************package*************************/
  Double_t min = 0.4;
  Double_t max = 10;
  // Double_t max = 5.0;
  Double_t loprecision = 0.01;
  Double_t hiprecision = 0.1;
  // Option_t *opt = "REBMS+";
  Option_t *opt = "RI+";
  TString logfilename = "log.root";
  Double_t minfit = 0.0;
  Double_t maxfit = 6;

  TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
  hout->Draw();

  // hout->GetXaxis()->SetRangeUser(0,20);

  for (Int_t ip = 0; ip < 9; ip++)
  {
    cout << hout->GetBinContent(ip + 1) << endl;
  }
}
