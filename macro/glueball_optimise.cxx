#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/common_glue.h"
#include "src/fitting_range_glue.h"

using namespace std;

Double_t BWsumMassDepWidth(double *x, double *par);
Double_t exponential_bkg_3(double *x, double *par);
Double_t BWsumMassDepWidth_exponential(double *x, double *par);
Double_t single_BW_mass_dep_spin0(double *x, double *par);
Double_t single_BW_mass_dep_spin2(double *x, double *par);
int colors[] = {kGreen + 4, 28, kMagenta, kBlue, kBlack};

void glueball_optimise()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    bool QAplots = false;
    TString inputFile = "582256";
    // TString inputFile = "586332";
    TString outputPath = "../output/glueball/LHC22o_pass7_small/" + inputFile + "/";
    gSystem->Exec("mkdir -p " + outputPath);
    TFile *file = new TFile(("../data/glueball/LHC22o_pass7_small/" + inputFile + ".root").Data(), "READ");
    if (file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << "\n";
        return;
    }
    THnSparseF *hSparseGlue = (THnSparseF *)file->Get("higher-mass-resonances/hglueball/h3glueInvMassDS");
    THnSparseF *hSparseRot = (THnSparseF *)file->Get("higher-mass-resonances/hglueball/h3glueInvMassRot"); // Rotated background
    // THnSparseF *hSparseRot = (THnSparseF *)file->Get("higher-mass-resonances/hglueball/h3glueInvMassME"); // Mixed-evet background
    if (hSparseGlue == nullptr || hSparseRot == nullptr)
    {
        std::cerr << "Error: Could not find required histograms in file\n";
        return;
    }
    THnSparseF *hSparseGlueClone = (THnSparseF *)hSparseGlue->Clone("hSparseGlueClone");
    THnSparseF *hSparseRotClone = (THnSparseF *)hSparseRot->Clone("hSparseRotClone");
    TH1D *ptCorr = hSparseGlue->Projection(5);
    TH1D *AngSep = hSparseGlue->Projection(4);
    TH1D *deltaM = hSparseGlue->Projection(3);

    vector<vector<float>> ptbins = {
        {1.0, 2.0},
        {2.0, 3.0},
        {3.0, 5.0},
        {5.0, 7.0},
        {7.0, 10.0}};

    float deltaRcut[] = {0.2, 0.2, 0.15, 0.15, 0.1};
    int totalbins = ptbins.size();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.SetTextFont(42);

    for (int ipt = 0; ipt < totalbins; ipt++)
    {
        int multLow = hSparseGlue->GetAxis(0)->FindBin(0.0 + 0.0001);
        int multHigh = hSparseGlue->GetAxis(0)->FindBin(100.0 - 0.0001);

        int ptLow = hSparseGlue->GetAxis(1)->FindBin(ptbins[ipt][0] + 0.01);
        int ptHigh = hSparseGlue->GetAxis(1)->FindBin(ptbins[ipt][1] - 0.01);

        int deltaMlow = hSparseGlue->GetAxis(3)->FindBin(0.0 + 0.00001);
        int deltaMhigh = hSparseGlue->GetAxis(3)->FindBin(0.01 - 0.00001);

        int angSepLow = hSparseGlue->GetAxis(4)->FindBin(deltaRcut[ipt] + 0.0001);
        int angSepHigh = hSparseGlue->GetAxis(4)->FindBin(50.0 - 0.0001);

        int ptCorrLow = hSparseGlue->GetAxis(5)->FindBin(0.0 + 0.0001);
        int ptCorrHigh = hSparseGlue->GetAxis(5)->FindBin(50.0 - 0.0001);

        // hSparseGlue->GetAxis(0)->SetRange(multLow, multHigh);
        hSparseGlue->GetAxis(1)->SetRange(ptLow, ptHigh);
        hSparseGlue->GetAxis(3)->SetRange(deltaMlow, deltaMhigh);
        hSparseGlue->GetAxis(4)->SetRange(angSepLow, angSepHigh);
        // hSparseGlue->GetAxis(5)->SetRange(ptCorrLow, ptCorrHigh);
        // hSparseRot->GetAxis(0)->SetRange(multLow, multHigh);
        hSparseRot->GetAxis(1)->SetRange(ptLow, ptHigh);
        hSparseRot->GetAxis(3)->SetRange(deltaMlow, deltaMhigh);
        hSparseRot->GetAxis(4)->SetRange(angSepLow, angSepHigh);
        // hSparseRot->GetAxis(5)->SetRange(ptCorrLow, ptCorrHigh);

        // hSparseGlueClone->GetAxis(0)->SetRange(multLow, multHigh);
        hSparseGlueClone->GetAxis(1)->SetRange(ptLow, ptHigh);
        // hSparseGlueClone->GetAxis(3)->SetRange(deltaMlow, deltaMhigh);
        // hSparseGlueClone->GetAxis(4)->SetRange(angSepLow, angSepHigh);
        // hSparseGlueClone->GetAxis(5)->SetRange(ptCorrLow, ptCorrHigh);
        // hSparseRotClone->GetAxis(0)->SetRange(multLow, multHigh);
        hSparseRotClone->GetAxis(1)->SetRange(ptLow, ptHigh);
        // hSparseRotClone->GetAxis(3)->SetRange(deltaMlow, deltaMhigh);
        // hSparseRotClone->GetAxis(4)->SetRange(angSepLow, angSepHigh);
        // hSparseRotClone->GetAxis(5)->SetRange(ptCorrLow, ptCorrHigh);

        TH1D *hGlueMass = hSparseGlue->Projection(2, "E");
        hGlueMass->SetName(Form("hGlueMass_pt%d", ipt));
        TH1D *hRotMass = hSparseRot->Projection(2, "E");
        hRotMass->SetName(Form("hRotMass_pt%d", ipt));
        TH1D *hGlueMassClone = (TH1D *)hSparseGlueClone->Projection(2, "E");
        hGlueMassClone->SetName(Form("hGlueMassClone_pt%d", ipt));
        TH1D *hRotMassClone = (TH1D *)hSparseRotClone->Projection(2, "E");
        hRotMassClone->SetName(Form("hRotMassClone_pt%d", ipt));

        float normLow = 2.8;
        float normHigh = 2.9;
        float normLowClone = 2.8;
        float normHighClone = 2.9;

        auto signalCounts = hGlueMass->Integral(hGlueMass->GetXaxis()->FindBin(normLow + 0.0001), hGlueMass->GetXaxis()->FindBin(normHigh - 0.0001));
        auto bkgCounts = hRotMass->Integral(hRotMass->GetXaxis()->FindBin(normLow + 0.0001), hRotMass->GetXaxis()->FindBin(normHigh - 0.0001));
        auto normFactor = signalCounts / bkgCounts;
        hRotMass->Scale(normFactor);

        auto signalCountsClone = hGlueMassClone->Integral(hGlueMassClone->GetXaxis()->FindBin(normLowClone + 0.0001), hGlueMassClone->GetXaxis()->FindBin(normHighClone - 0.0001));
        auto bkgCountsClone = hRotMassClone->Integral(hRotMassClone->GetXaxis()->FindBin(normLowClone + 0.0001), hRotMassClone->GetXaxis()->FindBin(normHighClone - 0.0001));
        auto normFactorClone = signalCountsClone / bkgCountsClone;
        hRotMassClone->Scale(normFactorClone);
        TH1F *hSignal = (TH1F *)hGlueMass->Clone(Form("hSignal_pt%d", ipt));
        hSignal->Add(hRotMass, -1);
        TH1F *hSignalClone = (TH1F *)hGlueMassClone->Clone(Form("hSignalClone_pt%d", ipt));
        hSignalClone->Add(hRotMassClone, -1);

        int rebinFactor = 2;
        // TCanvas *cGlueballMass = new TCanvas("", "Glueball Invariant Mass", 720, 720);
        // SetCanvasStyle(cGlueballMass, 0.14, 0.03, 0.05, 0.13);
        // SetHistoQA(hGlueMassClone);
        // hGlueMassClone->GetXaxis()->SetTitle("Invariant Mass (GeV/#it{c}^{2})");
        // hGlueMassClone->SetMarkerStyle(20);
        // hGlueMassClone->GetYaxis()->SetMaxDigits(3);
        // hGlueMassClone->GetYaxis()->SetTitleOffset(1.4);
        // hGlueMassClone->SetMarkerSize(1.0);
        // hGlueMassClone->Rebin(rebinFactor);
        // hGlueMassClone->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hGlueMassClone->GetBinWidth(1) * 1000));
        // hGlueMassClone->Draw("pe");
        // SetHistoQA(hGlueMass);
        // hGlueMass->Rebin(rebinFactor);
        // hGlueMass->SetLineColor(kBlue);
        // hGlueMass->SetMarkerColor(kBlue);
        // hGlueMass->SetMarkerStyle(20);
        // hGlueMass->SetMarkerSize(1.0);
        // hGlueMass->Draw("pE same");
        // hRotMass->SetMarkerStyle(21);
        // hRotMass->SetMarkerColor(kRed);
        // hRotMass->SetLineColor(kRed);
        // hRotMass->SetMarkerSize(1.0);
        // hRotMass->Draw("pE same");
        // auto statsLost = 100 * (hGlueMassClone->GetEntries() - hGlueMass->GetEntries()) / hGlueMassClone->GetEntries();
        // TLegend *legend = new TLegend(0.45, 0.75, 0.85, 0.9);
        // legend->SetTextSize(0.035);
        // legend->SetTextFont(42);
        // legend->SetFillStyle(0);
        // legend->SetBorderSize(0);
        // // legend->AddEntry(hGlueMass, "Same event K^{0}_{s}K^{0}_{s} pair", "lp");
        // // legend->AddEntry(hRotMass, "Rotated event K^{0}_{s}K^{0}_{s} pair", "lp");
        // legend->AddEntry(hGlueMassClone, "No cut", "lp");
        // legend->AddEntry(hGlueMass, "With cut", "lp");
        // legend->AddEntry((TObject *)0, Form("Statistics lost: %.1f %%", statsLost), "");
        // legend->Draw();

        TCanvas *cSignal = new TCanvas("", "Glueball Signal", 720, 720);
        SetCanvasStyle(cSignal, 0.14, 0.03, 0.05, 0.13);
        SetHistoQA(hSignalClone);
        hSignalClone->GetXaxis()->SetTitle("Invariant Mass (GeV/#it{c}^{2})");
        hSignalClone->SetMarkerStyle(20);
        hSignalClone->GetYaxis()->SetMaxDigits(3);
        hSignalClone->GetYaxis()->SetTitleOffset(1.4);
        hSignalClone->SetMarkerSize(1.0);
        hSignalClone->Rebin(rebinFactor);
        hSignalClone->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hSignalClone->GetBinWidth(1) * 1000));
        hSignalClone->GetXaxis()->SetRangeUser(1.0, 2.2);
        hSignalClone->Draw("pe");
        SetHistoQA(hSignal);
        hSignal->SetLineColor(kBlue);
        hSignal->SetMarkerColor(kBlue);
        hSignal->SetMarkerStyle(20);
        hSignal->SetMarkerSize(1.0);
        hSignal->Rebin(rebinFactor);
        hSignal->GetXaxis()->SetRangeUser(1.0, 2.2);
        hSignal->Draw("pE same");
        // legend->Draw();
        // TLine *lineaty0 = new TLine(1.0, 0, 2.2, 0);
        // lineaty0->SetLineColor(kBlack);
        // lineaty0->SetLineStyle(2);
        // lineaty0->Draw("same");
        lat.DrawLatex(0.5, 0.8, Form("p_{T} bin: %.1f - %.1f GeV/#it{c}", ptbins[ipt][0], ptbins[ipt][1]));
        auto statsLostSignal = 100 * (hSignalClone->GetEntries() - hSignal->GetEntries()) / hSignalClone->GetEntries();
        lat.DrawLatex(0.5, 0.74, Form("Statistics lost: %.1f %%", statsLostSignal));
    }

    /*
    // Fitting the distribution
    TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_exponential, 1.1, 2.20, 16);
    string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
    for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
    {
        BEexpol->SetParName(i, parnames[i].c_str());
    }
    double parameters[] = {214, f1270Mass, f1270Width, 50, a1320Mass, a1320Width, 450, f1525Mass, f1525Width, 130, f1710Mass, f1710Width};
    int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

    for (int i = 0; i < size_fitparams; i++)
    {
        BEexpol->SetParameter(i, parameters[i]);
    }

    double initial_param_bkg[] = {1.37518e8, 2.0, 10.67, 0.8}; // rebin twice (2022 dataset)
    BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]);
    BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]);
    BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]);
    BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]);

    // Limits on the amplitudes
    // BEexpol->SetParLimits(0, 0, 1e6);
    // BEexpol->SetParLimits(3, 0, 1e6);
    // BEexpol->SetParLimits(6, 0, 1e6);
    // BEexpol->SetParLimits(9, 0, 1e6);

    // Fixed parameters
    BEexpol->FixParameter(2, f1270Width);
    BEexpol->FixParameter(5, a1320Width);
    BEexpol->FixParameter(8, f1525Width);
    BEexpol->FixParameter(11, f1710Width);

    // Free parameters
    BEexpol->FixParameter(1, f1270Mass);
    BEexpol->FixParameter(4, a1320Mass);
    BEexpol->FixParameter(7, f1525Mass);
    BEexpol->FixParameter(10, f1710Mass);

    TFitResultPtr fitResultptr = hSignal->Fit("BEexpol", "REMBS");
    double *obtained_parameters = BEexpol->GetParameters();

    TF1 *singlefits[4];
    for (int i = 0; i < 4; i++)
    {
        singlefits[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, 1.00, 3.0, 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, 1.00, 3.0, 3);
        singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
        singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
        singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
        singlefits[i]->SetLineColor(colors[i]);
        singlefits[i]->SetLineStyle(2);
        singlefits[i]->Draw("same");
    }

    TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
    for (int i = 0; i < 4; i++)
    {
        expol->SetParameter(i, obtained_parameters[size_fitparams + i]);
    }
    expol->SetLineColor(3);
    expol->SetLineStyle(2);
    expol->Draw("same");

    TLegend *ltemp = new TLegend(0.25, 0.52, 0.55, 0.87);
    ltemp->SetFillStyle(0);
    ltemp->SetBorderSize(0);
    ltemp->SetTextFont(42);
    ltemp->SetTextSize(0.03);
    ltemp->AddEntry(hSignal, "Data (stat. uncert.)", "lpe");
    ltemp->AddEntry(BEexpol, "4rBW + Residual BG", "l");
    ltemp->AddEntry(expol, "Residual BG", "l");
    ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
    ltemp->AddEntry(singlefits[1], "a_{2}(1320)^{0}", "l");
    ltemp->AddEntry(singlefits[2], "f'_{2}(1525)", "l");
    ltemp->AddEntry(singlefits[3], "f_{0}(1710)", "l");
    ltemp->Draw();
    */

    if (QAplots)
    {
        TCanvas *cptCorr = new TCanvas("cptCorr", "p_{T} correlation", 720, 720);
        SetCanvasStyle(cptCorr, 0.14, 0.03, 0.05, 0.13);
        SetHistoQA(ptCorr);
        ptCorr->GetXaxis()->SetTitle("pt correlation");
        ptCorr->GetYaxis()->SetTitleOffset(1.4);
        ptCorr->GetYaxis()->SetTitle("Counts");
        ptCorr->Draw("pe");

        TCanvas *cAngSep = new TCanvas("cAngSep", "Angular Separation", 720, 720);
        SetCanvasStyle(cAngSep, 0.14, 0.03, 0.05, 0.13);
        SetHistoQA(AngSep);
        AngSep->GetXaxis()->SetTitle("Angular Separation (rad)");
        AngSep->GetYaxis()->SetTitleOffset(1.4);
        AngSep->GetYaxis()->SetTitle("Counts");
        AngSep->Draw("pe");

        TCanvas *cdeltaM = new TCanvas("cdeltaM", "Delta M", 720, 720);
        SetCanvasStyle(cdeltaM, 0.14, 0.03, 0.05, 0.13);
        SetHistoQA(deltaM);
        deltaM->GetXaxis()->SetTitle("#DeltaM (GeV/#it{c}^{2})");
        deltaM->GetYaxis()->SetTitleOffset(1.4);
        deltaM->GetYaxis()->SetTitle("Counts");
        deltaM->Draw("pe");
    }
} //******** end of main function ***********

//==============================================================//
//======================== Fit functions =======================//
//==============================================================//

Double_t BWsumMassDepWidth(double *x, double *par)
{
    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
    double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    // double phase_space = (x[0] / TMath::Sqrt(x[0] * x[0] + 15 * 15)) * (TMath::Exp(-TMath::Sqrt(x[0] * x[0] + 15 * 15) / 0.16)); // 160 MeV is the kinetic freeze-out temperature

    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2] * (TMath::Power(par[1] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double yield1320 = par[3];
    double mass1320 = par[4];
    double width1320 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double yield1525 = par[6];
    double mass1525 = par[7];
    double width1525 = par[8] * (TMath::Power(par[7] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n1);
    double yield1710 = par[9];
    double mass1710 = par[10];
    double width1710 = par[11] * (TMath::Power(par[10] / x[0], 1.0)) * TMath::Power((npart1) / (dpart4), n2);

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1320 = yield1320 * mass1320 * width1320 * x[0] / (pow((x[0] * x[0] - mass1320 * mass1320), 2) + pow(mass1320 * width1320, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    double fit = (fit1270 + fit1320 + fit1525 + fit1710);
    return fit;
}

Double_t exponential_bkg_3(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[3])));
}

Double_t BWsumMassDepWidth_exponential(double *x, double *par)
{
    // return (BWsumMassDepWidth(x, par) + expol_chkstar(x, &par[12]));
    return (BWsumMassDepWidth(x, par) + exponential_bkg_3(x, &par[12]));
}

Double_t single_BW_mass_dep_spin0(double *x, double *par)
{
    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    int spin = 0;
    double n1 = (2.0 * spin + 1.0) / 2.0;

    double yield = par[0];
    double mass = par[1];
    double width = par[2] * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));

    return fit;
}

Double_t single_BW_mass_dep_spin2(double *x, double *par)
{
    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    int spin = 2;
    double n1 = (2.0 * spin + 1.0) / 2.0;

    double yield = par[0];
    double mass = par[1];
    double width = par[2] * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));

    return fit;
}