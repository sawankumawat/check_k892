#include <iostream>
#include "../src/style.h"
#include "../src/initializations.h"
#include "YieldMean.C"

using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
}

void compare_yield()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TFile *fpublished = new TFile("pp13TeV.root", "READ");
    TGraph *grpublished = (TGraph *)fpublished->Get("Table 4/Graph1D_y1");
    grpublished->Scale(1.0 / 0.7448);
    if (!grpublished)
    {
        cout << "Graph not found" << endl;
        return;
    }

    string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/447406/kstarqa_id21631/hInvMass/";
    TFile *spectra1 = new TFile((path + "corrected_spectra.root").c_str(), "READ");
    TFile *spectra2 = new TFile((path + "corrected_spectra.root").c_str(), "READ");

    if (spectra1->IsZombie() || spectra2->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    TH1D *h1 = (TH1D *)spectra1->Get("corrected_spectra_Integral");
    TH1D *h2 = (TH1D *)spectra2->Get("corrected_spectra_Integral");

    if (h1 == nullptr || h2 == nullptr)
    {
        cout << "Histogram not found" << endl;
        return;
    }

    // // // h1->Scale(0.5);
    // TF1 *fLevyTsallis = new TF1("levyTsallis", levyTsallis, 0, 15.5, 3); // Adjust range and parameter count as needed
    // // Set initial parameter values
    // fLevyTsallis->SetParameters(0.155025, 6.91579, 0.219293); // Example values: q, T, N, mass
    // fLevyTsallis->SetParNames("dN/dy", "T", "n");
    // fLevyTsallis->SetLineWidth(2);
    // fLevyTsallis->SetLineColor(kBlue);

    // // //Fit the function to the histogram
    // h1->Fit("levyTsallis", "REBMS"); // "R" option to use the function range

    for (int i = 1; i <= h1->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr = (0.1 * h1->GetBinContent(i));
        h1->SetBinError(i, systemerr);
    }

    TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 20.0, 4);
    fitFcn->SetParameter(0, 5.0);
    fitFcn->SetParameter(1, 0.0007);
    fitFcn->FixParameter(2, 1.285);
    fitFcn->SetParameter(3, 0.3);
    fitFcn->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn->SetLineColor(kBlue);

    /*************meanpT*****************byresonance*******************package*************************/
    Double_t min = 0;
    Double_t max = 20;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "RI+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 20;
    // Double_t maxfit=8.0;

    TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // hout->Draw();

    // TCanvas *c = new TCanvas();
    // c->SetLogy(1);
    // h1->Draw();

    cout << "run 2 bins from graph size is " << grpublished->GetN() << endl;
    cout << "run 3 bins are " << h1->GetNbinsX() << endl; // one extra bin from 15 GeV/c2 to 20 GeV/c2

    TGraphErrors *gratio = new TGraphErrors();
    for (int i = 0; i < grpublished->GetN(); i++)
    {
        double x_run2, yield_run2, x_error, y_error_run2;
        grpublished->GetPoint(i, x_run2, yield_run2);
        x_error = grpublished->GetErrorX(i);
        y_error_run2 = grpublished->GetErrorY(i);

        double bincontent = h1->GetBinContent(i + 1);
        double binerror = h1->GetBinError(i + 1);
        cout << "Bin " << i + 1 << " content: " << bincontent << " error: " << binerror << ", run 2 content: " << yield_run2 << " error: " << y_error_run2 << endl;
        double ratio = bincontent / yield_run2;
        double error = sqrt(pow(binerror / yield_run2, 2) + pow(bincontent * y_error_run2 / (yield_run2 * yield_run2), 2));
        cout<<"Ratio error is "<<error<<endl;
        gratio->SetPoint(i, x_run2, ratio);
        gratio->SetPointError(i, x_error, error);

        // double thisanalysis = fitFcn->Eval(x);
        // cout<<"x is "<<x<<" y is "<<published<<endl;
        // gratio->SetPoint(i, x, thisanalysis / published);
        // double error = (thisanalysis / (published * published)) * sqrt(pow(gr->GetErrorY(i), 2));
        // gratio->SetPointError(i, 0, error);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    SetHistoQA(h1);
    h1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    h1->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    h1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    h1->SetMaximum(h1->GetMaximum() * 1.5);
    // h1->SetMinimum(h1->GetMinimum() * 50);
    h1->GetYaxis()->SetTitleOffset(1.30);
    h1->GetXaxis()->SetTitleOffset(1.02);
    h1->SetMarkerStyle(22);
    h1->SetMarkerSize(1.5);
    h1->GetXaxis()->SetRangeUser(0, 15.5);
    h1->Draw("pe");
    gPad->SetLogy(1);
    grpublished->SetMarkerStyle(29);
    grpublished->SetMarkerSize(1.5);
    grpublished->SetMarkerColor(kRed);
    grpublished->SetLineColor(kRed);
    grpublished->Draw("pe same");

    TLegend *leg = new TLegend(0.48, 0.5, 0.8, 0.8);
    SetLegendStyle(leg);
    // leg->AddEntry(h1, "pp 13.6 TeV (TPC + TOF)", "lpe");
    leg->AddEntry(h1, "pp 13.6 TeV", "lpe");
    leg->AddEntry(fitFcn, "Levy-Tsallis fit (pp 13.6 TeV)", "l");
    leg->AddEntry(grpublished, "pp 13 TeV (Published)", "lpe");
    leg->SetTextSize(0.05);
    leg->Draw();

    c1->cd(2);
    TH1F *hdummy = (TH1F *)h1->Clone();
    for (int i = 0; i < hdummy->GetNbinsX(); i++)
    {
        hdummy->SetBinContent(i + 1, 0);
        hdummy->SetBinError(i + 1, 0);
    }

    SetGrapherrorStyle(gratio);
    gratio->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    gratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    gratio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    gratio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    gratio->SetMarkerStyle(23);
    gratio->SetMarkerSize(1.5);
    gratio->SetMarkerColor(kRed);
    gratio->SetLineColor(kRed);
    gratio->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
    gratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gratio->GetXaxis()->CenterTitle(1);
    gratio->GetYaxis()->SetTitleOffset(0.6);
    gratio->GetXaxis()->SetTitleOffset(1.1);
    gratio->GetYaxis()->SetNdivisions(506);
    gratio->GetXaxis()->SetRangeUser(0, 15.5);
    gratio->Draw("ap");

    TLine *line = new TLine(0, 1, 15, 1);
    line->SetLineStyle(2);
    line->Draw();

    c1->SaveAs((path + "Ratio_publishedKstar.png").c_str());
}