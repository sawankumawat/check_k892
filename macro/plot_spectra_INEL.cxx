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
    string filePath = "../output/kstar/LHC22o_pass7/708297/kstarqa/hInvMass/";
    TFile *fINEL = OpenFile(filePath + "corrected_spectra_0_120.root");
    TH1D *hSpectraINEL = GetHisto(fINEL, "mult_0-120/corrected_spectra_Integral_final");

    TH1F *h1 = (TH1F *)hSpectraINEL->Clone("h1");
    TH1F *h2 = (TH1F *)hSpectraINEL->Clone("h2");

    for (int i = 1; i <= h2->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr = (h2->GetBinContent(i) * 0.08); // Assuming 8% systematic uncertainty
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

    TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    fitFcn->SetLineColor(kRed + 1);
    fitFcn->SetLineWidth(3);
    fitFcn->SetLineStyle(2);

    TCanvas *cSpectraINEL = new TCanvas("cSpectraINEL", "", 720, 720);
    SetCanvasStyle(cSpectraINEL, 0.16, 0.06, 0.01, 0.14);
    gPad->SetLogy();
    SetHistoQA(hSpectraINEL);
    hSpectraINEL->SetMaximum(0.2);
    hSpectraINEL->SetMinimum(8e-9);
    hSpectraINEL->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSpectraINEL->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    hSpectraINEL->GetYaxis()->SetTitleOffset(1.6);
    hSpectraINEL->SetStats(0);
    hSpectraINEL->SetMarkerStyle(20);
    hSpectraINEL->SetMarkerSize(1.2);
    hSpectraINEL->Draw("pe");
    fitFcn->Draw("l same");
    TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.88);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry((TObject *)0, "K*(892)^{0}", "");
    leg->AddEntry(hSpectraINEL, "pp, #sqrt{s} = 13.6 TeV", "p");
    leg->AddEntry(fitFcn, "L#acute{e}vy-Tsallis", "l");
    leg->Draw();
    cSpectraINEL->SaveAs((filePath + "/spectraFit_INEL_0-120.png").c_str());
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
