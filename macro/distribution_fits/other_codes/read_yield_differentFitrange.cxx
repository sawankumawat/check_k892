#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"
#include "TF1.h"
#include "TMath.h"
#include "../spectra/YieldMean2.C"

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

void read_yield_differentFitrange()
{
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra";
    TFile *file = new TFile((path + "/ReweightedSpectra.root").c_str(), "read");
    if (file->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hYield1525Corrected = (TH1F *)file->Get("f21525_Reweighted_Yield");
    TH1F *hYield1710Corrected = (TH1F *)file->Get("f01710_Reweighted_Yield");
    if (hYield1525Corrected == nullptr || hYield1710Corrected == nullptr)
    {
        cout << "Histograms not found" << endl;
        return;
    }
    hYield1525Corrected->GetYaxis()->SetTitle("1/N_{evt} #times d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hYield1525Corrected->SetMaximum(4e-2);
    hYield1525Corrected->SetMinimum(5e-7);
    // hYield1525Corrected->SetMinimum(2e-5); // for zoomed view
    hYield1525Corrected->SetMarkerSize(1.5);
    hYield1525Corrected->GetYaxis()->SetTitleOffset(1.6);
    hYield1710Corrected->GetYaxis()->SetTitle("1/N_{evt} #times d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hYield1710Corrected->SetMaximum(4e-2);
    hYield1710Corrected->SetMinimum(9e-7);
    // hYield1710Corrected->SetMinimum(2e-5); // for zoomed view
    hYield1710Corrected->SetMarkerSize(1.5);
    hYield1710Corrected->GetYaxis()->SetTitleOffset(1.6);
    // hYield1525Corrected->GetXaxis()->SetRangeUser(0.0, 7.0); // for zoomed view
    // hYield1710Corrected->GetXaxis()->SetRangeUser(0.0, 7.0); // for zoomed view

    // Clone histograms for fitting
    TH1F *h1f2 = (TH1F *)hYield1525Corrected->Clone("h1f2");
    TH1F *h2f2 = (TH1F *)hYield1525Corrected->Clone("h2f2");
    TH1F *h1f0 = (TH1F *)hYield1710Corrected->Clone("h1f0");
    TH1F *h2f0 = (TH1F *)hYield1710Corrected->Clone("h2f0");

    // Add systematic errors
    for (int i = 1; i <= h1f2->GetNbinsX(); i++)
    {
        double systemerrf2 = (0.1 * h2f2->GetBinContent(i)); // 10% systematic error
        h2f2->SetBinError(i, systemerrf2);
        double systemerrf0 = (0.1 * h2f0->GetBinContent(i)); // 10% systematic error
        h2f0->SetBinError(i, systemerrf0);
    }

    // Fit parameters
    Double_t min = 0.0;
    Double_t max = 5.0;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    // Option_t *opt = "RI0+";
    Option_t *opt = "REBMS0+";
    TString logfilename = "log.root";
    Double_t minfit = 1.0;
    Double_t maxfit = 10.0;

    // Fit parameters for Levy
    Double_t min2 = 0.0;
    Double_t max2 = 10.0;
    // Option_t *opt2 = "RI0+";
    Option_t *opt2 = "REBMS0+";
    Double_t minfit2 = 1.0;
    Double_t maxfit2 = 10.0;

    // Fit parameters for Levy2
    Double_t min3 = 0.0;
    Double_t max3 = 15.0;
    Option_t *opt3 = "REBMS0+";
    Double_t minfit3 = 1.0;
    Double_t maxfit3 = 15.0;

    // Fit parameters for Exponential
    Double_t minExp = 0.0;
    Double_t maxExp = 10.0;
    Double_t minfitExp = 1.0;
    Double_t maxfitExp = 10.0;

    //====================================Levy-Tsallis Fit==========================================
    TF1 *fitLevyf2 = new TF1("fitLevyf2", FuncLavy, 0.0, 15.0, 4);
    fitLevyf2->SetParameter(0, 5.0);
    fitLevyf2->SetParameter(1, 0.5);
    fitLevyf2->FixParameter(2, 1.525);
    fitLevyf2->SetParameter(3, 0.35);
    fitLevyf2->SetParNames("n", "dn/dy", "mass", "T");
    // h1f2->Fit(fitLevyf2, opt, "", minfit, maxfit);
    TH1 *houtLevyf2 = YieldMean(h1f2, h2f2, fitLevyf2, min2, max2, loprecision, hiprecision, opt2, logfilename, minfit2, maxfit2);

    TF1 *fitLevyf0 = new TF1("fitLevyf0", FuncLavy, 0.0, 15.0, 4);
    fitLevyf0->SetParameter(0, 5.0);
    fitLevyf0->SetParameter(1, 0.5);
    fitLevyf0->FixParameter(2, 1.710);
    fitLevyf0->SetParameter(3, 0.35);
    fitLevyf0->SetParNames("n", "dn/dy", "mass", "T");
    // h1f0->Fit(fitLevyf0, opt, "", minfit, maxfit);
    TH1 *houtLevyf0 = YieldMean(h1f0, h2f0, fitLevyf0, min2, max2, loprecision, hiprecision, opt2, logfilename, minfit2, maxfit2);

    TF1 *fitLevyf2_v2 = new TF1("fitLevyf2_v2", FuncLavy, 0.0, 15.0, 4);
    fitLevyf2_v2->SetParameter(0, 5.0);
    fitLevyf2_v2->SetParameter(1, 0.5);
    fitLevyf2_v2->FixParameter(2, 1.525);
    fitLevyf2_v2->SetParameter(3, 0.35);
    fitLevyf2_v2->SetParNames("n", "dn/dy", "mass", "T");
    // h1f2->Fit(fitLevyf2_v2, opt, "", minfit, maxfit);
    TH1 *houtLevyf2_v2 = YieldMean(h1f2, h2f2, fitLevyf2_v2, min3, max3, loprecision, hiprecision, opt3, logfilename, minfit3, maxfit3);

    TF1 *fitLevyf0_v2 = new TF1("fitLevyf0_v2", FuncLavy, 0.0, 15.0, 4);
    fitLevyf0_v2->SetParameter(0, 5.0);
    fitLevyf0_v2->SetParameter(1, 0.5);
    fitLevyf0_v2->FixParameter(2, 1.710);
    fitLevyf0_v2->SetParameter(3, 0.35);
    fitLevyf0_v2->SetParNames("n", "dn/dy", "mass", "T");
    // h1f0->Fit(fitLevyf0_v2, opt, "", minfit, maxfit);
    TH1 *houtLevyf0_v2 = YieldMean(h1f0, h2f0, fitLevyf0_v2, min3, max3, loprecision, hiprecision, opt3, logfilename, minfit3, maxfit3);

    // //====================================Boltzmann Fit==========================================
    // TF1 *fitBoltzmannf2 = new TF1("fitBoltzmannf2", FuncBoltzmanndNdptTimesPt, 0.0, 10.0, 3);
    // fitBoltzmannf2->SetParameter(0, 0.01);
    // fitBoltzmannf2->SetParameter(1, 0.5);
    // fitBoltzmannf2->FixParameter(2, 1.525);
    // fitBoltzmannf2->SetParNames("norm", "T", "mass");
    // // h1f2->Fit(fitBoltzmannf2, opt, "", minfit, maxfit);
    // TH1 *houtBoltf2 = YieldMean(h1f2, h2f2, fitBoltzmannf2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    // TF1 *fitBoltzmannf0 = new TF1("fitBoltzmannf0", FuncBoltzmanndNdptTimesPt, 0.0, 10.0, 3);
    // fitBoltzmannf0->SetParameter(0, 0.1);
    // fitBoltzmannf0->SetParameter(1, 0.5);
    // fitBoltzmannf0->FixParameter(2, 1.710);
    // fitBoltzmannf0->SetParNames("norm", "T", "mass");
    // // h1f0->Fit(fitBoltzmannf0, opt, "", minfit, maxfit);
    // TH1 *houtBoltf0 = YieldMean(h1f0, h2f0, fitBoltzmannf0, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    // //====================================Exponential Fit==========================================
    // TF1 *fitExpf2 = new TF1("fitExpf2", FuncExpdNdptTimesPt, 0.0, 10.0, 2);
    // fitExpf2->SetParameter(0, 0.1);
    // fitExpf2->SetParameter(1, 0.5);
    // fitExpf2->SetParNames("norm", "T");
    // // h1f2->Fit(fitExpf2, opt, "", minfit, maxfit);
    // TH1 *houtExpo2 = YieldMean(h1f2, h2f2, fitExpf2, minExp, maxExp, loprecision, hiprecision, opt, logfilename, minfitExp, maxfitExp);

    // TF1 *fitExpf0 = new TF1("fitExpf0", FuncExpdNdptTimesPt, 0.0, 10.0, 2);
    // fitExpf0->SetParameter(0, 0.1);
    // fitExpf0->SetParameter(1, 0.5);
    // fitExpf0->SetParNames("norm", "T");
    // // h1f0->Fit(fitExpf0, opt, "", minfitExp, maxfitExp);
    // TH1 *houtExpof0 = YieldMean(h1f0, h2f0, fitExpf0, minExp, maxExp, loprecision, hiprecision, opt, logfilename, minfitExp, maxfitExp);

    // // // //====================================mT Exponential Fit==========================================
    // // // TF1 *fitMTExpf2 = new TF1("fitMTExpf2", FuncMTExpdNdptTimesPt, 0.0, 10.0, 3);
    // // // fitMTExpf2->SetParameter(0, 0.1);
    // // // fitMTExpf2->SetParameter(1, 0.5);
    // // // fitMTExpf2->FixParameter(2, 1.525);
    // // // fitMTExpf2->SetParNames("norm", "T", "mass");
    // // // h1f2->Fit(fitMTExpf2, opt, "", minfit, maxfit);
    // // // TF1 *fitMTExpf0 = new TF1("fitMTExpf0", FuncMTExpdNdptTimesPt, 0.0, 10.0, 3);
    // // // fitMTExpf0->SetParameter(0, 0.1);
    // // // fitMTExpf0->SetParameter(1, 0.5);
    // // // fitMTExpf0->FixParameter(2, 1.710);
    // // // fitMTExpf0->SetParNames("norm", "T", "mass");
    // // // h1f0->Fit(fitMTExpf0, opt, "", minfit, maxfit);

    // //====================================Bose-Einstein Fit==========================================
    // TF1 *fitBoseEinsteinf2 = new TF1("fitBoseEinsteinf2", FuncBoseEinsteindNdptTimesPt, 0.0, 10.0, 3);
    // fitBoseEinsteinf2->SetParameter(0, 0.5);
    // fitBoseEinsteinf2->SetParameter(1, 0.5);
    // fitBoseEinsteinf2->FixParameter(2, 1.525);
    // fitBoseEinsteinf2->SetParNames("norm", "T", "mass");
    // // h1f2->Fit(fitBoseEinsteinf2, opt, "", minfit, maxfit);
    // TH1 *houtBEf2 = YieldMean(h1f2, h2f2, fitBoseEinsteinf2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    // TF1 *fitBoseEinsteinf0 = new TF1("fitBoseEinsteinf0", FuncBoseEinsteindNdptTimesPt, 0.0, 10.0, 3);
    // fitBoseEinsteinf0->SetParameter(0, 0.5);
    // fitBoseEinsteinf0->SetParameter(1, 0.5);
    // fitBoseEinsteinf0->FixParameter(2, 1.710);
    // fitBoseEinsteinf0->SetParNames("norm", "T", "mass");
    // // h1f0->Fit(fitBoseEinsteinf0, opt, "", minfit, maxfit);
    // TH1 *houtBEf0 = YieldMean(h1f0, h2f0, fitBoseEinsteinf0, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    // //====================================Blast-Wave Fit==========================================
    // TF1 *fitBWf2 = new TF1("fitBWf2", BGBlastWave_Func, 0.0, 7.0, 5);
    // fitBWf2->FixParameter(0, 1.525);
    // fitBWf2->FixParameter(1, 0.9);
    // // fitBWf2->SetParLimits(1, 0.1, 1);
    // fitBWf2->SetParameter(2, 0.15);
    // fitBWf2->SetParLimits(2, 0.05, 0.9);
    // fitBWf2->SetParameter(3, 0.95);
    // fitBWf2->SetParLimits(3, 0.1, 10);
    // fitBWf2->SetParameter(4, 0.1);
    // fitBWf2->SetParNames("mass", "beta_max", "T", "n", "norm");
    // // h1f2->Fit(fitBWf2, opt, "", minfit, maxfit);
    // TH1 *houtBWf2 = YieldMean(h1f2, h2f2, fitBWf2, minExp, maxExp, loprecision, hiprecision, opt, logfilename, minfitExp, maxfitExp);

    // TF1 *fitBWf0 = new TF1("fitBWf0", BGBlastWave_Func, 0.0, 7.0, 5);
    // fitBWf0->FixParameter(0, 1.710);
    // fitBWf0->FixParameter(1, 0.9);
    // // fitBWf0->SetParLimits(1, 0.1, 1);
    // fitBWf0->SetParameter(2, 0.15);
    // fitBWf0->SetParLimits(2, 0.05, 0.9);
    // fitBWf0->SetParameter(3, 0.95);
    // fitBWf0->SetParLimits(3, 0.1, 10);
    // fitBWf0->SetParameter(4, 0.1);
    // fitBWf0->SetParNames("mass", "beta_max", "T", "n", "norm");
    // // h1f0->Fit(fitBWf0, opt, "", minfit, maxfit);
    // TH1 *houtBWf0 = YieldMean(h1f0, h2f0, fitBWf0, minExp, maxExp, loprecision, hiprecision, opt, logfilename, minfitExp, maxfitExp);

    //====================================Plotting All Fits on Same Canvas for f2(1525)==========================================
    TCanvas *cFitf2All = new TCanvas("cFitf2All", "All Fits for f2'(1525)", 720, 720);
    SetCanvasStyle(cFitf2All, 0.17, 0.03, 0.05, 0.14);
    gPad->SetLogy();

    h1f2->SetMarkerStyle(20);
    h1f2->SetMarkerColor(kBlack);
    h1f2->SetLineColor(kBlack);
    h1f2->SetMinimum(1e-7);
    h1f2->Draw("pe");

    // Set different colors for eachloprecision fit function
    fitLevyf2->SetLineColor(kRed);
    fitLevyf2->SetLineWidth(2);
    fitLevyf2->Draw("same");

    fitLevyf2_v2->SetLineColor(kBlue);
    fitLevyf2_v2->SetLineWidth(2);
    fitLevyf2_v2->Draw("same");

    // fitBoltzmannf2->SetLineColor(kBlue);
    // fitBoltzmannf2->SetLineWidth(2);
    // fitBoltzmannf2->Draw("same");

    // fitExpf2->SetLineColor(kGreen + 2);
    // fitExpf2->SetLineWidth(2);
    // fitExpf2->Draw("same");

    // // fitMTExpf2->SetLineColor(kMagenta);
    // // fitMTExpf2->SetLineWidth(2);
    // // fitMTExpf2->Draw("same");

    // fitBoseEinsteinf2->SetLineColor(kCyan);
    // fitBoseEinsteinf2->SetLineWidth(2);
    // fitBoseEinsteinf2->Draw("same");

    // // fitBWf2->SetLineColor(kMagenta);
    // // fitBWf2->SetLineWidth(2);
    // // fitBWf2->Draw("same");

    // Create legend
    TLegend *legf2 = new TLegend(0.55, 0.55, 0.95, 0.92);
    legf2->AddEntry((TObject *)0, "ALICE", "");
    legf2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legf2->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    legf2->AddEntry(h1f2, "f_{2}'(1525) spectra", "pe");
    legf2->AddEntry((TObject *)0, "Levy-Tsallis fit", "");
    legf2->AddEntry(fitLevyf2, "Fit: 0-10 GeV/c", "l");
    legf2->AddEntry(fitLevyf2_v2, "Fit: 0-15 GeV/c", "l");

    // legf2->AddEntry(fitBoltzmannf2, "Boltzmann", "l");
    // legf2->AddEntry(fitExpf2, "Exponential", "l");
    // // legf2->AddEntry(fitMTExpf2, "mT Exponential", "l");
    // legf2->AddEntry(fitBoseEinsteinf2, "Bose-Einstein", "l");
    // // legf2->AddEntry(fitBWf2, "Blast-Wave", "l");
    legf2->SetBorderSize(0);
    legf2->SetFillStyle(0);
    legf2->SetTextSize(0.03);
    legf2->Draw();
    cFitf2All->SaveAs((path + "/DifferentFitFunc/FitRangeCompare_f2.png").c_str());

    //====================================Plotting All Fits on Same Canvas for f0(1710)==========================================
    TCanvas *cFitf0All = new TCanvas("cFitf0All", "All Fits for f0(1710)", 720, 720);
    SetCanvasStyle(cFitf0All, 0.17, 0.03, 0.05, 0.14);
    gPad->SetLogy();

    h1f0->SetMarkerStyle(20);
    h1f0->SetMarkerColor(kBlack);
    h1f0->SetLineColor(kBlack);
    h1f0->SetMinimum(5e-7);
    h1f0->Draw("pe");

    // Set different colors for each fit function
    fitLevyf0->SetLineColor(kRed);
    fitLevyf0->SetLineWidth(2);
    fitLevyf0->Draw("same");

    fitLevyf0_v2->SetLineColor(kBlue);
    fitLevyf0_v2->SetLineWidth(2);
    fitLevyf0_v2->Draw("same");

    // fitBoltzmannf0->SetLineColor(kBlue);
    // fitBoltzmannf0->SetLineWidth(2);
    // fitBoltzmannf0->Draw("same");

    // fitExpf0->SetLineColor(kGreen + 2);
    // fitExpf0->SetLineWidth(2);
    // fitExpf0->Draw("same");

    // // fitMTExpf0->SetLineColor(kMagenta);
    // // fitMTExpf0->SetLineWidth(2);
    // // fitMTExpf0->Draw("same");

    // fitBoseEinsteinf0->SetLineColor(kCyan);
    // fitBoseEinsteinf0->SetLineWidth(2);
    // fitBoseEinsteinf0->Draw("same");

    // // fitBWf0->SetLineColor(kMagenta);
    // // fitBWf0->SetLineWidth(2);
    // // fitBWf0->Draw("same");

    // Create legend
    TLegend *legf0 = new TLegend(0.55, 0.6, 0.95, 0.92);
    legf0->AddEntry((TObject *)0, "ALICE", "");
    legf0->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legf0->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    legf0->AddEntry(h1f0, "f_{0}(1710) spectra", "pe");
    legf0->AddEntry((TObject *)0, "Levy-Tsallis fit", "");
    legf0->AddEntry(fitLevyf0, "Fit: 0-10 GeV/c", "l");
    legf0->AddEntry(fitLevyf0_v2, "Fit: 0-15 GeV/c", "l");
    // legf0->AddEntry(fitLevyf0, "Levy-Tsallis", "l");
    // legf0->AddEntry(fitBoltzmannf0, "Boltzmann", "l");
    // legf0->AddEntry(fitExpf0, "Exponential", "l");
    // // legf0->AddEntry(fitMTExpf0, "mT Exponential", "l");
    // legf0->AddEntry(fitBoseEinsteinf0, "Bose-Einstein", "l");
    // // legf0->AddEntry(fitBWf0, "Blast-Wave", "l");
    legf0->SetBorderSize(0);
    legf0->SetFillStyle(0);
    legf0->SetTextSize(0.03);
    legf0->Draw();
    cFitf0All->SaveAs((path + "/DifferentFitFunc/FitRangeCompare_f0.png").c_str());

    // Print Chi2/NDF, dN/dy and <pT> for each fit function
    // Lambda function to avoid code repetition
    auto printFitResults = [](TF1 *fitFunc, TH1 *hout, const string &fitName, const string &mesonName)
    {
        cout << fitName << " " << mesonName << " , Chi2 " << fitFunc->GetChisquare() << " NDF " << fitFunc->GetNDF() << endl;
        cout << "dN/dy : " << hout->GetBinContent(1) << " +/- " << hout->GetBinContent(2) << endl;
        cout << "<pT> : " << hout->GetBinContent(5) << " +/- " << hout->GetBinContent(6) << endl;
        cout << "Low pT extrapolation contribution (%): " << hout->GetBinContent(10) * 100 << endl;
        cout << endl;
    };
    cout << "Fit results for f2(1525):" << endl;
    printFitResults(fitLevyf2, houtLevyf2, "Levy-Tsallis", "f2(1525)");
    printFitResults(fitLevyf2_v2, houtLevyf2_v2, "Levy-Tsallis (0-15 GeV/c)", "f2(1525)");
    // printFitResults(fitBoltzmannf2, houtBoltf2, "Boltzmann", "f2(1525)");
    // printFitResults(fitExpf2, houtExpo2, "Exponential", "f2(1525)");
    // printFitResults(fitBoseEinsteinf2, houtBEf2, "Bose-Einstein", "f2(1525)");
    // printFitResults(fitBWf2, houtBWf2, "Blast-Wave", "f2(1525)");

    cout << "Fit results for f0(1710):" << endl;
    printFitResults(fitLevyf0, houtLevyf0, "Levy-Tsallis", "f0(1710)");
    printFitResults(fitLevyf0_v2, houtLevyf0_v2, "Levy-Tsallis (0-15 GeV/c)", "f0(1710)");
    // printFitResults(fitBoltzmannf0, houtBoltf0, "Boltzmann", "f0(1710)");
    // printFitResults(fitExpf0, houtExpof0, "Exponential", "f0(1710)");
    // printFitResults(fitBoseEinsteinf0, houtBEf0, "Bose-Einstein", "f0(1710)");
    // printFitResults(fitBWf0, houtBWf0, "Blast-Wave", "f0(1710)");
}

//=======================Fit functions=========================
Double_t FuncLavy(Double_t *x, Double_t *par)
{
    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

Double_t BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{
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
        fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0.0, 7.0, 5);
    fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
    Double_t integral = fBGBlastWave_Integrand->Integral(0.0, 7.0, 1.e4);
    return norm * pt * integral;
}

Double_t FuncBoltzmanndNdptTimesPt(Double_t *x, Double_t *par)
{
    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    return norm * pT * mT * TMath::Exp(-mT / T) * TMath::Exp(mass / T) * (1.0 / (2 * T * T * T + 2 * T * T * mass + T * mass * mass));
}

Double_t FuncExpdNdptTimesPt(Double_t *x, Double_t *par)
{
    Double_t pT = x[0];
    Double_t norm = par[0];
    Double_t T = par[1];

    Double_t p = norm * pT * TMath::Exp(-pT / T) * (1.0 / (T * T));
    return p;
}

Double_t FuncMTExpdNdptTimesPt(Double_t *x, Double_t *par)
{
    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    return norm * pT * TMath::Exp(-mT / T) * TMath::Exp(mass / T) * (1.0 / (T * T + T * mass));
}

Double_t FuncBoseEinsteindNdptTimesPt(Double_t *x, Double_t *par)
{
    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    return norm * pT / (TMath::Exp(mT / T) - 1.0) * (TMath::Exp(mass / T) - 1.0);
}
