#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/initializations.h"
#include "spectra/YieldMean.C"

TFile *OpenFile(const string &path);
TH1D *GetHisto(TFile *f, const string &name);
Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

using namespace std;

void plot_spectra_INELgt0()
{
    bool isRaw = false;

    string filePath = "../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/";
    TFile *fINELgt0 = (isRaw) ? OpenFile(filePath + "yield_0_100.root") : OpenFile(filePath + "corrected_spectra_0_100.root");
    TH1D *hSpectraINELgt0 = (isRaw) ? GetHisto(fINELgt0, "mult_0-100/yield_integral") : GetHisto(fINELgt0, "mult_0-100/corrected_spectra_Integral_final");

    TFile *fSysUncertINELTemp = OpenFile("../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots/SysUncert.root");
    TH1D *hSysINELTemp = GetHisto(fSysUncertINELTemp, "hTotalSysSmoothed_0_100");

    TH1F *h1 = (TH1F *)hSpectraINELgt0->Clone("h1");
    TH1F *h2 = (TH1F *)hSpectraINELgt0->Clone("h2");

    for (int i = 1; i <= h2->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr = (h2->GetBinContent(i) * hSysINELTemp->GetBinContent(i));
        h2->SetBinError(i, systemerr);
    }

    Double_t min = 0.0;
    Double_t max = 20.0;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.5;
    Option_t *opt = "RI0+";
    TString logfilename = "log_fit.root";
    Double_t minfit = 0.0;
    Double_t maxfit = 20.0;

    TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 20.0, 4);
    fitFcn->SetParameter(0, 5.0);
    fitFcn->SetParameter(1, 0.05);
    // fitFcn->SetParameter(1, 0.5);
    fitFcn->FixParameter(2, 0.895);
    fitFcn->SetParameter(3, 0.15);
    fitFcn->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn->SetLineColor(kRed + 1);
    // h2->Fit(fitFcn, "RI0+");

    double meanpT, meanpT_errStat, meanpT_errSys;
    double yield, yield_errStat, yield_errSys;

    if (!isRaw)
    {
        TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
        fitFcn->SetLineColor(kRed + 1);
        fitFcn->SetLineWidth(3);
        fitFcn->SetLineStyle(2);

        meanpT = hout->GetBinContent(5);
        meanpT_errStat = hout->GetBinContent(6);
        meanpT_errSys = hout->GetBinContent(7);
        yield = hout->GetBinContent(1);
        yield_errStat = hout->GetBinContent(2);
        yield_errSys = hout->GetBinContent(3);

        cout << "meanpT: " << meanpT << " +/- " << meanpT_errStat << " (stat) +/- " << meanpT_errSys << " (sys)" << endl;
        cout << "yield: " << yield << " +/- " << yield_errStat << " (stat) +/- " << yield_errSys << " (sys)" << endl;
    }

    TCanvas *cSpectraINEL = new TCanvas("cSpectraINEL", "", 720, 720);
    SetCanvasStyle(cSpectraINEL, 0.17, 0.06, 0.01, 0.14);
    gPad->SetLogy();
    SetHistoQA(hSpectraINELgt0);
    hSpectraINELgt0->SetMaximum(0.2);
    hSpectraINELgt0->SetMinimum(8e-9);
    hSpectraINELgt0->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSpectraINELgt0->GetYaxis()->SetTitle("1/#it{N}_{Ev} d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    hSpectraINELgt0->GetYaxis()->SetTitleOffset(1.6);
    hSpectraINELgt0->SetStats(0);
    hSpectraINELgt0->SetMarkerStyle(20);
    hSpectraINELgt0->SetMarkerSize(1.2);
    hSpectraINELgt0->SetMarkerColor(kBlue - 1);
    hSpectraINELgt0->SetLineColor(kBlue - 1);
    hSpectraINELgt0->Draw("pe");
    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(1.2);
    h2->SetMarkerColor(kBlack);
    h2->SetLineColor(kBlack);
    h2->SetFillStyle(0);
    h2->SetLineWidth(2);
    if (!isRaw)
    {
        h2->Draw("e2 same");
        fitFcn->Draw("l same");
    }
    TLegend *leg = new TLegend(0.68, 0.8, 0.9, 0.96);
    leg->SetTextSize(0.032);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hSpectraINELgt0, "INEL Spectra", "p");
    if (!isRaw)
    {
        leg->AddEntry(fitFcn, "L#acute{e}vy-Tsallis", "l");
    }

    // leg->AddEntry((TObject *)0, "|y| < 0.5", "");
    // leg->AddEntry((TObject *)0, "INEL", "");
    // leg->AddEntry((TObject *)0, "K*(892)^{0}", "");
    leg->Draw();

    TLatex *lat = new TLatex();
    lat->SetTextSize(0.032);
    lat->SetTextFont(42);
    lat->SetNDC();
    lat->DrawLatex(0.37, 0.91, "ALICE");
    lat->DrawLatex(0.37, 0.85, "pp, #sqrt{s} = 13.6 TeV");
    lat->DrawLatex(0.37, 0.79, "|y| < 0.5");
    lat->DrawLatex(0.37, 0.73, "K*(892)^{0}");
    lat->SetTextSize(0.03);
    if (!isRaw)
        lat->DrawLatex(0.2, 0.2, "Uncertainties: Stat. (bars), Sys. (boxes)");

    // (isRaw) ? cSpectraINEL->SaveAs((filePath + "/Rawspectra_INEL_0-100.png").c_str()) : cSpectraINEL->SaveAs((filePath + "/spectraFit_INEL_0-100.png").c_str());
}
TH1D *GetHisto(TFile *f, const string &name)
{
    TH1D *histo = (TH1D *)f->Get(name.c_str());

    if (!histo || histo == nullptr)
    {
        cout << "Error: histo " << name << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetHistoQA(histo);
    histo->SetTitle(0);
    return histo;
}

TFile *OpenFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (f->IsZombie())
    {
        cout << "Error: File not found: " << path << endl;
        return nullptr;
    }
    return f;
}
