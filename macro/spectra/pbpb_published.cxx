#include <iostream>
#include "style.h"
#include "fitfunc.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.12);
    pad2->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.04);
}

void pbpb_published()
{
    TFile *fspec = new TFile("spectra.root", "read");
    TH1D *h1 = (TH1D *)fspec->Get("lf-k892analysis/K892/0/hCorrectedYields");

    float ptbin[] = {0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12, 16, 20};
    float values[] = {3.06447, 2.71771, 2.25728, 1.63196, 0.895651, 0.487635, 0.19833, 0.107873, 0.0401064, 0.0112694, 0.00382887, 0.0011622, 0.000398064, 0.000156181, 3.47542e-05};
    float value_err[] = {0.321725, 0.266959, 0.2196, 0.186548, 0.117899, 0.0674138, 0.0222473, 0.0115173, 0.0042446, 0.00118458, 0.000384392, 0.000117726, 4.57373e-05, 1.95185e-05, 4.83167e-06};

    cout << "The size of ptbin is: " << sizeof(ptbin) / sizeof(ptbin[0]) << endl;
    cout << "The size of values is: " << sizeof(values) / sizeof(values[0]) << endl;

    TH1F *hist = new TH1F("hist", "hist", sizeof(values) / sizeof(values[0]), ptbin);
    for (int i = 0; i < sizeof(values) / sizeof(values[0]); i++)
    {
        hist->SetBinContent(i + 1, values[i]);
        hist->SetBinError(i + 1, value_err[i]);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    SetHistoQA(hist);

    SetHistoStyle(hist, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    hist->GetYaxis()->SetTitleSize(0.04 / pad1Size);

    c1->cd(1);
    gPad->SetLogy(1);
    hist->SetMarkerSize(1.5);
    h1->SetMarkerSize(1.5);
    hist->SetMarkerStyle(22);
    hist->Draw("ep");
    h1->SetMarkerColor(kRed);
    h1->SetLineColor(kRed);
    h1->Draw("same");

    TF1 *fLevyTsallis = new TF1("levyTsallis", levyTsallis, 0, 10, 3); // Adjust range and parameter count as needed
    gStyle->SetOptFit(1);
    // Set initial parameter values
    fLevyTsallis->SetParameters(0.055025, 6.91579, 0.219293); // Example values: q, T, N, mass
    fLevyTsallis->SetParNames("dN/dy", "T", "n");
    fLevyTsallis->SetLineStyle(2);

    // Fit the function to the histogram
    h1->Fit("levyTsallis", "RM"); // "R" option to use the function range

    TH1F *hratio = new TH1F("", "", sizeof(values) / sizeof(values[0]), ptbin);
    TH1F *hratio2 = new TH1F("", "", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());

    for (int i = 0; i < sizeof(values) / sizeof(values[0]) - 3; i++)
    {
        double valpublished = hist->GetBinContent(i + 1);
        double valanalysis = fLevyTsallis->Eval((ptbin[i] + ptbin[i + 1]) / 2.0);
        double ratioval = valpublished / valanalysis;
        // double ratioerr = sqrt(pow(value_err[i] / valpublished, 2) + pow(fLevyTsallis->GetParError(0) / valanalysis, 2));
        double ratioerr = value_err[i] / valpublished;
        hratio->SetBinContent(i + 1, ratioval);
        hratio->SetBinError(i + 1, ratioerr);
    }

    for(int i = 0; i< h1->GetNbinsX(); i++){
        double funcvalue = fLevyTsallis->Eval(h1->GetBinCenter(i+1));
        double ratioval = h1->GetBinContent(i+1) / funcvalue;
        double ratioerr = h1->GetBinError(i+1) / funcvalue;
        hratio2->SetBinContent(i+1, ratioval);
        hratio2->SetBinError(i+1, ratioerr);
    }

    c1->cd(2);
    SetHistoStyle(hratio, 1, 53, 1, 0.05, 0.05, 0.04 / pad2Size, 0.04 / pad2Size, 1.13, 1.8);
    hratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);

    hratio->Draw("ep");
    hratio2->SetMarkerSize(1.5);
    hratio->SetMarkerStyle(22);
    hratio2->SetMarkerStyle(20);
    hratio2->SetMarkerColor(kRed);
    hratio2->SetLineColor(kRed);
    hratio2->Draw("ep same");

}