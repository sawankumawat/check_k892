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

void plot_spectra_INEL()
{
    bool isRaw = false;

    string filePath = "../output/kstar/LHC22o_pass7/788987/kstarqa/hInvMass/ROTATED/";
    TFile *fINEL = (isRaw) ? OpenFile(filePath + "yield_0_120.root") : OpenFile(filePath + "corrected_spectra_0_120.root");
    TH1D *hSpectraINEL = (isRaw) ? GetHisto(fINEL, "mult_0-120/yield_integral") : GetHisto(fINEL, "mult_0-120/corrected_spectra_Integral_final");

    if (!isRaw)
    {
        hSpectraINEL->Scale(1.0 / 0.514912); // Event loss removed
        hSpectraINEL->Scale(1.0 / 0.677);    // Visible cross section ratio
        double ratioSel_TVX = 9235908442.0 / 12478172244.0;
        hSpectraINEL->Scale(ratioSel_TVX); // selected events ==> TVX events
    }

    TFile *fSysUncertINELTemp = OpenFile("../output/kstar/LHC22o_pass7/756343/kstarqa/hInvMass/SystematicsPlots/SysUncert.root");
    TH1D *hSysINELTemp = GetHisto(fSysUncertINELTemp, "hTotalSysSmoothed_0_120");

    TH1F *h1 = (TH1F *)hSpectraINEL->Clone("h1");
    TH1F *h2 = (TH1F *)hSpectraINEL->Clone("h2");

    for (int i = 1; i <= h2->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr = (h2->GetBinContent(i) * hSysINELTemp->GetBinContent(i));
        h2->SetBinError(i, systemerr);
    }

    Double_t min = 0.0;
    Double_t max = 30.0;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.5;
    Option_t *opt = "RI0+";
    TString logfilename = "log_fit.root";
    Double_t minfit = 0.0;
    Double_t maxfit = 30.0;

    TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 30.0, 4);
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
    SetHistoQA(hSpectraINEL);
    hSpectraINEL->SetMaximum(0.9);
    hSpectraINEL->SetMinimum(8e-9);
    hSpectraINEL->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSpectraINEL->GetYaxis()->SetTitle("1/#it{N}_{Ev} d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    hSpectraINEL->GetYaxis()->SetTitleOffset(1.6);
    hSpectraINEL->SetStats(0);
    hSpectraINEL->SetMarkerStyle(20);
    hSpectraINEL->SetMarkerSize(1.2);
    hSpectraINEL->SetMarkerColor(kBlue - 1);
    hSpectraINEL->SetLineColor(kBlue - 1);
    hSpectraINEL->Draw("pe");
    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(1.2);
    h2->SetMarkerColor(kBlack);
    h2->SetLineColor(kBlack);
    h2->SetFillStyle(0);
    h2->SetLineWidth(2);
    // if (!isRaw)
    // {
    //     h2->Draw("e2 same");
    //     fitFcn->Draw("l same");
    // }
    TLegend *leg = new TLegend(0.68, 0.8, 0.9, 0.96);
    leg->SetTextSize(0.032);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hSpectraINEL, "INEL Spectra", "p");
    // if (!isRaw)
    // {
    //     leg->AddEntry(fitFcn, "L#acute{e}vy-Tsallis", "l");
    // }

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
    // if (!isRaw)
    //     lat->DrawLatex(0.2, 0.2, "Uncertainties: Stat. (bars), Sys. (boxes)");

    (isRaw) ? cSpectraINEL->SaveAs((filePath + "/Rawspectra_INEL_0-120.png").c_str()) : cSpectraINEL->SaveAs((filePath + "/spectra_INEL_0-120.png").c_str());
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
