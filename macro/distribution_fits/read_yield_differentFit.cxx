#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"
#include "TF1.h"
#include "TMath.h"
#include "../spectra/YieldMean.C"

using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par);
Double_t BGBlastWave_Integrand(const Double_t *x, const Double_t *p);
Double_t BGBlastWave_Func(const Double_t *x, const Double_t *p);
Double_t FuncBoltzmanndNdptTimesPt(Double_t *x, Double_t *par);
Double_t FuncExpdNdptTimesPt(Double_t *x, Double_t *par);
Double_t FuncMTExpdNdptTimesPt(Double_t *x, Double_t *par);
Double_t FuncBoseEinsteindNdptTimesPt(Double_t *x, Double_t *par);

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;

void read_yield_differentFit()
{
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra";
    TFile *file = new TFile((path + "/spectra.root").c_str(), "read");
    if (file->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hYield1525Corrected = (TH1F *)file->Get("hYield1525Corrected");
    TH1F *hYield1710Corrected = (TH1F *)file->Get("hYield1710Corrected");
    if (hYield1525Corrected == nullptr || hYield1710Corrected == nullptr)
    {
        cout << "Histograms not found" << endl;
        return;
    }

    // first fit with levy-Tsallis fit
    TH1F *h1f2 = (TH1F *)hYield1525Corrected->Clone("h1f2");
    TH1F *h2f2 = (TH1F *)hYield1525Corrected->Clone("h2f2");
    TH1F *h1f0 = (TH1F *)hYield1710Corrected->Clone("h1f0");
    TH1F *h2f0 = (TH1F *)hYield1710Corrected->Clone("h2f0");

    for (int i = 1; i <= h1f2->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerrf2 = (0.1 * h2f2->GetBinContent(i)); // 10% systematic error
        h2f2->SetBinError(i, systemerrf2);
        double systemerrf0 = (0.1 * h2f0->GetBinContent(i)); // 10% systematic error
        h2f0->SetBinError(i, systemerrf0);
    }

    //=====================================Fit functions==========================================

    Double_t min = 0.0;
    Double_t max = 5.0;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    // Option_t *opt = "REBMS0+";
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0.0;
    Double_t maxfit = 5.0;

    // //====================================Levy-Tsallis Fit==========================================

    // TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
    // fitFcn->SetParameter(0, 5.0);
    // // fitFcn->SetParameter(1, 0.05);
    // fitFcn->SetParameter(1, 0.5);
    // fitFcn->FixParameter(2, 1.525);
    // fitFcn->SetParameter(3, 0.35);
    // fitFcn->SetParNames("n", "dn/dy", "mass", "T");

    // TH1 *hout = YieldMean(h1f2, h2f2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // cout << "Yield dN/dy = " << hout->GetBinContent(1) << " +- " << hout->GetBinContent(2) << endl;
    // cout << "Mean pT = " << hout->GetBinContent(5) << " +- " << hout->GetBinContent(6) << endl;

    // TF1 *fitFcn2 = new TF1("fitfunc2", FuncLavy, 0.0, 15.0, 4);
    // fitFcn2->SetParameter(0, 5.0);
    // // fitFcn2->SetParameter(1, 0.05);
    // fitFcn2->SetParameter(1, 0.5);
    // fitFcn2->FixParameter(2, 1.710);
    // fitFcn2->SetParameter(3, 0.35);
    // fitFcn2->SetParNames("n", "dn/dy", "mass", "T");

    // TH1 *hout2 = YieldMean(h1f0, h2f0, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // cout << "Yield dN/dy of f0(1710) = " << hout2->GetBinContent(1) << " +- " << hout2->GetBinContent(2) << endl;
    // cout << "Mean pT of f0(1710) = " << hout2->GetBinContent(5) << " +- " << hout2->GetBinContent(6) << endl;

    // //====================================Blast-Wave Fit==========================================
    // TF1 *fitFcn = new TF1("fitFcn", BGBlastWave_Func, 0.0, 5, 5);

    // fitFcn->FixParameter(0, 1.525); // mass
    // // fitFcn->SetParameter(1,0.78);//Beta max
    // fitFcn->SetParLimits(1, 0.1, 1);    // Beta max   //0-3.5 arb par
    // fitFcn->SetParameter(2, 0.15);      // T
    // fitFcn->SetParLimits(2, 0.08, 0.9); // T
    // fitFcn->FixParameter(3, 0.95);      // n
    // fitFcn->SetParLimits(3, 0.1, 5);    // n
    // fitFcn->SetParameter(4, 1.e6);      // norm
    // fitFcn->SetParNames("mass", "beta_max", "T", "n", "norm");

    // TH1 *hout = YieldMean(h1f2, h2f2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // cout << "Yield dN/dy of f2(1525) = " << hout->GetBinContent(1) << " +- " << hout->GetBinContent(2) << endl;
    // cout << "Mean pT of f2(1525) = " << hout->GetBinContent(5) << " +- " << hout->GetBinContent(6) << endl;

    // TF1 *fitFcn2 = new TF1("fitFcn2", BGBlastWave_Func, 0.0, 5.0, 5);
    // fitFcn2->FixParameter(0, 1.710); // mass
    // // fitFcn2->SetParameter(1,0.78);//Beta max
    // fitFcn2->SetParLimits(1, 0.1, 1);    // Beta max   //0-3.5 arb par
    // fitFcn2->SetParameter(2, 0.15);      // T
    // fitFcn2->SetParLimits(2, 0.08, 0.9); // T
    // fitFcn2->FixParameter(3, 0.95);      // n
    // fitFcn2->SetParLimits(3, 0.1, 5);    // n
    // fitFcn2->SetParameter(4, 1.e6);      // norm
    // fitFcn2->SetParNames("mass", "beta_max", "T", "n", "norm");

    // TH1 *hout2 = YieldMean(h1f0, h2f0, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // cout << "Yield dN/dy of f0(1710) = " << hout2->GetBinContent(1) << " +- " << hout2->GetBinContent(2) << endl;
    // cout << "Mean pT of f0(1710) = " << hout2->GetBinContent(5) << " +- " << hout2->GetBinContent(6) << endl;

    // //====================================Boltzmann Fit==========================================
    // TF1 *fitFcn = new TF1("fitBoltzmann", FuncBoltzmanndNdptTimesPt, 0.0, 5.0, 3);
    // fitFcn->SetParameter(0, 1.e6);  // norm
    // fitFcn->SetParameter(1, 0.15);  // T
    // fitFcn->FixParameter(2, 1.525); // mass
    // fitFcn->SetParNames("norm", "T", "mass");

    // TH1 *hout = YieldMean(h1f2, h2f2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // cout << "Yield dN/dy of f2(1525) = " << hout->GetBinContent(1) << " +- " << hout->GetBinContent(2) << endl;
    // cout << "Mean pT of f2(1525) = " << hout->GetBinContent(5) << " +- " << hout->GetBinContent(6) << endl;

    // TF1 *fitFcn2 = new TF1("fitBoltzmann2", FuncBoltzmanndNdptTimesPt, 0.0, 5.0, 3);
    // fitFcn2->SetParameter(0, 1.e6);  // norm
    // fitFcn2->SetParameter(1, 0.15);  // T
    // fitFcn2->FixParameter(2, 1.710); // mass
    // fitFcn2->SetParNames("norm", "T", "mass");
    // TH1 *hout2 = YieldMean(h1f0, h2f0, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // cout << "Yield dN/dy of f0(1710) = " << hout2->GetBinContent(1) << " +- " << hout2->GetBinContent(2) << endl;
    // cout << "Mean pT of f0(1710) = " << hout2->GetBinContent(5) << " +- " << hout2->GetBinContent(6) << endl;

    // //================================Exponential Fit==========================================
    // TF1 *fitFcn = new TF1("fitExponential", FuncExpdNdptTimesPt, 0.0, 6.0, 2);
    // fitFcn->SetParameter(0, 1.e6);  // norm
    // fitFcn->SetParameter(1, 0.15);  // T
    // fitFcn->FixParameter(2, 1.525); // mass
    // fitFcn->SetParNames("norm", "T", "mass");

    // TH1 *hout = YieldMean(h1f2, h2f2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // cout << "Yield dN/dy of f2(1525) = " << hout->GetBinContent(1) << " +- " << hout->GetBinContent(2) << endl;
    // cout << "Mean pT of f2(1525) = " << hout->GetBinContent(5) << " +- " << hout->GetBinContent(6) << endl;

    // TF1 *fitFcn2 = new TF1("fitExponential2", FuncExpdNdptTimesPt, 0.0, 6.0, 2);
    // fitFcn2->SetParameter(0, 1.e6);  // norm
    // fitFcn2->SetParameter(1, 0.15);  // T
    // fitFcn2->FixParameter(2, 1.710); // mass
    // fitFcn2->SetParNames("norm", "T", "mass");
    // TH1 *hout2 = YieldMean(h1f0, h2f0, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // cout << "Yield dN/dy of f0(1710) = " << hout2->GetBinContent(1) << " +- " << hout2->GetBinContent(2) << endl;
    // cout << "Mean pT of f0(1710) = " << hout2->GetBinContent(5) << " +- " << hout2->GetBinContent(6) << endl;

    TCanvas *cFitf2 = new TCanvas("cFitf2", "Levy-Tsallis Fit for f2'(1525)", 720, 720);
    SetCanvasStyle(cFitf2, 0.14, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    h1f2->SetMarkerStyle(20);
    h1f2->SetMarkerColor(kBlue);
    h1f2->SetLineColor(kBlue);
    h1f2->Draw("pe");
    fitFcn->SetLineColor(kRed);
    fitFcn->Draw("same");

    TCanvas *cFitf0 = new TCanvas("cFitf0", "Levy-Tsallis Fit for f0(1710)", 720, 720);
    SetCanvasStyle(cFitf0, 0.14, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    h1f0->SetMarkerStyle(20);
    h1f0->SetMarkerColor(kBlue);
    h1f0->SetLineColor(kBlue);
    h1f0->Draw("pe");
    fitFcn2->SetLineColor(kRed);
    fitFcn2->Draw("same");
}

//=======================Fit functions=========================
Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}
Double_t BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
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

// Boltzmann in dN/dpT times pT (fix the mass while using the function)
Double_t FuncBoltzmanndNdptTimesPt(Double_t *x, Double_t *par)
{

    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    // set mass explicitly
    // const Double_t mass = 0.497611; // example: K0s

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    // return norm * pT * mT * TMath::Exp(-mT / T); // Not normalized
    return norm * pT * mT * TMath::Exp(-mT / T) * TMath::Exp(mass / T) * (1.0 / (2 * T * T * T + 2 * T * T * mass + T * mass * mass)); // normalized
}

// Exponential in dN/dpT times pT
Double_t FuncExpdNdptTimesPt(Double_t *x, Double_t *par)
{
    // par[0] = norm
    // par[1] = T
    // x[0]   = pT

    Double_t pT = x[0];
    Double_t norm = par[0];
    Double_t T = par[1];

    // Double_t p = norm * pT * TMath::Exp(-pT / T); // not normalized
    Double_t p = norm * pT * TMath::Exp(-pT / T) * (1.0 / (T * T)); // normalized
    return p;
}

// Exponential in dN/dmT times pT (fix the mass while using the function)
Double_t FuncMTExpdNdptTimesPt(Double_t *x, Double_t *par)
{

    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    // set mass explicitly (example: K0s)
    // const Double_t mass = 0.497611;

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    // return norm * pT * TMath::Exp(-mT / T); // not normalized
    return norm * pT * TMath::Exp(-mT / T) * TMath::Exp(mass / T) * (1.0 / (T * T + T * mass)); // normalized
}

// Bose-Einstein in dN/dpT times pT (fix the mass while using the function)
Double_t FuncBoseEinsteindNdptTimesPt(Double_t *x, Double_t *par)
{

    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    // set mass explicitly
    // const Double_t mass = 0.497611; // example: K0s

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    return norm * pT / (TMath::Exp(mT / T) - 1.0); // not normalized
    // return norm * pT / (TMath::Exp(mT / T) - 1.0) * (TMath::Exp(mass / T) - 1.0); // normalized (may not be correct)
}