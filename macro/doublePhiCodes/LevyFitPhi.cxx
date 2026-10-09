using namespace std;
#include "../src/style.h"
#include "../spectra/YieldMean.C"

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void LevyFitPhi()
{
    TFile *file = new TFile("Phi13TeV_INEL.root");
    TH1D *hPhi = (TH1D *)file->Get("Table 5/Hist1D_y1");
    TH1D *hPhiStat = (TH1D *)file->Get("Table 5/Hist1D_y1_e1");
    TH1D *hPhiSys = (TH1D *)file->Get("Table 5/Hist1D_y1_e2");

    int totalPoints = hPhi->GetNbinsX();
    for (int i = 1; i <= totalPoints; i++)
    {
        double x = hPhi->GetBinCenter(i);
        double y = hPhi->GetBinContent(i);
        double yErr = hPhiSys->GetBinContent(i);
        double yStatErr = hPhiStat->GetBinContent(i);

        hPhiSys->SetBinContent(i, y);
        hPhiSys->SetBinError(i, yErr);

        hPhiStat->SetBinContent(i, y);
        hPhiStat->SetBinError(i, yStatErr);
    }

    TCanvas *c1 = new TCanvas("c1", "Levy Fit for Phi Meson", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.05, 0.05, 0.15);
    gPad->SetLogy();
    SetHistoQA(hPhiSys);
    hPhiSys->Draw("E1");

    TF1 *fitFcn1 = new TF1("fitfunc", FuncLavy, 0.0, 11.0, 4);
    fitFcn1->SetParameter(0, 5.0);
    fitFcn1->SetParameter(1, 0.07);
    fitFcn1->FixParameter(2, 1.01946);
    fitFcn1->SetParameter(3, 0.3);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn1->SetLineColor(kBlack);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);

    TH1D *h1 = (TH1D *)hPhiStat->Clone("h1");
    TH1D *h2 = (TH1D *)hPhiSys->Clone("h2");

    // hPhiSys->Fit(fitFcn1, "R");

    /*************meanpT*****************byresonance*******************package*************************/
    Double_t min = 0.0;
    Double_t max = 11;
    Double_t loprecision = 0.001;
    Double_t hiprecision = 0.1;
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0.0;
    Double_t maxfit = 11;

    TH1 *hout = YieldMean(h1, h2, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    c1->cd();
    fitFcn1->Draw("same");
    cout << "Levy Fit Parameters for Phi Meson: " << endl;
    cout << "n: " << fitFcn1->GetParameter(0) << endl;
    cout << "dn/dy: " << fitFcn1->GetParameter(1) << endl;
    cout << "mass: " << fitFcn1->GetParameter(2) << endl;
    cout << "T: " << fitFcn1->GetParameter(3) << endl;
}