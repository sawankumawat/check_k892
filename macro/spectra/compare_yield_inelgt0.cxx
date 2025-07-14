#include <iostream>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../levy_fits/YieldMean.C"

using namespace std;

string names[] = {
    "", // Dump
    "",
    "_DCAxy_loose",
    "_all_manual_cuts",
    "_bins_cent",
    "_manual_trk_sel",
    "_mult_calib",
    "_no_pvcontributor",
    "_tof_high_pt",
    "_tpc_cluster120"};

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
    pad1->SetLeftMargin(0.14);
    pad2->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.02);
}

void compare_yield_inelgt0()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TFile *fpub = new TFile("pp13TeV_INELgt0.root", "READ");
    if (fpub->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }

    // choose any table except Table 1 for the INELgt0 case, since for Table 1, there are less bins.

    TGraphErrors *g_vom = (TGraphErrors *)fpub->Get(Form("Table %d/Graph1D_y1", 4));                   // vom1
    TGraphErrors *g_ratio_vom_inelgt0 = (TGraphErrors *)fpub->Get(Form("Table %d/Graph1D_y1", 4 + 9)); // ratio: vom1 / INEL>0

    if (g_vom == nullptr || g_ratio_vom_inelgt0 == nullptr)
    {
        cout << "Graph not found" << endl;
        return;
    }
    TGraphErrors *gr_inelgt0 = new TGraphErrors();
    for (int i = 0; i < g_vom->GetN(); i++)
    {
        double x, y, xerr, yerr;
        g_vom->GetPoint(i, x, y);
        xerr = g_vom->GetErrorX(i);
        yerr = g_vom->GetErrorY(i);

        double x1, y1, xerr1, yerr1;
        g_ratio_vom_inelgt0->GetPoint(i, x1, y1);
        xerr1 = g_ratio_vom_inelgt0->GetErrorX(i);
        yerr1 = g_ratio_vom_inelgt0->GetErrorY(i);

        if (x != x1)
        {
            cout << "x values do not match" << endl;
            return;
        }

        double inelgt0 = y / y1;
        double errory = sqrt(pow(yerr / y1, 2) + pow(y * yerr1 / (y1 * y1), 2));

        gr_inelgt0->SetPoint(i, x, inelgt0);
        gr_inelgt0->SetPointError(i, xerr, errory);
    }

    int spectrano = 0;
    // TFile *spectra = new TFile(("spectra_" + to_string(spectrano) + ".root").c_str(), "READ");
    TFile *spectra1 = new TFile("/home/sawan/k892_postprocessing/output/LHC22o_pass7_INEL/common/spectra.root", "READ");
    TFile *spectra2 = new TFile("/home/sawan/k892_postprocessing/output/LHC22o_pass7_INEL/common/spectra_1.root", "READ");

    if (spectra1->IsZombie() || spectra2->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    TH1D *h1 = (TH1D *)spectra2->Get(("lf-k892analysis_TOF" + names[spectrano] + "/K892/0/hCorrectedYields").c_str());
    TH1D *h2 = (TH1D *)spectra1->Get(("lf-k892analysis" + names[spectrano] + "/K892/0/hCorrectedYields").c_str());
    if (h1 == nullptr || h2 == nullptr)
    {
        cout << "Histogram not found" << endl;
        return;
    }

    h1->Scale(0.5);
    h2->Scale(0.5);
    TF1 *fLevyTsallis = new TF1("levyTsallis", levyTsallis, 0, 13.5, 3); // Adjust range and parameter count as needed
    // Set initial parameter values
    fLevyTsallis->SetParameters(0.155025, 6.91579, 0.219293); // Example values: q, T, N, mass
    fLevyTsallis->SetParNames("dN/dy", "T", "n");
    fLevyTsallis->SetLineWidth(2);
    fLevyTsallis->SetLineColor(kBlue);

    // // //Fit the function to the histogram
    // h1->Fit("levyTsallis", "REBMS"); // "R" option to use the function range

    for (int i = 1; i <= h1->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr = (0.1 * h1->GetBinContent(i));
        h1->SetBinError(i, systemerr);
    }

    TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
    fitFcn->SetParameter(0, 5.0);
    fitFcn->SetParameter(1, 0.0007);
    fitFcn->FixParameter(2, 1.285);
    fitFcn->SetParameter(3, 0.3);
    fitFcn->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn->SetLineColor(kBlue);

    /*************meanpT*****************byresonance*******************package*************************/
    Double_t min = 0;
    Double_t max = 15;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "RI+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 20;
    // Double_t maxfit=8.0;

    TH1 *hout = YieldMean(h1, h1, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    // TCanvas *c = new TCanvas();
    // c->SetLogy(1);
    // h1->Draw();

    cout << "run 2 bins from graph size is " << gr_inelgt0->GetN() << endl;
    cout << "run 3 bins are " << h1->GetNbinsX() << endl; // one extra bin from 15 GeV/c2 to 20 GeV/c2

    TGraphErrors *gratio = new TGraphErrors();
    for (int i = 0; i < gr_inelgt0->GetN(); i++)
    {
        double x_run2, yield_run2, x_error, y_error_run2;
        gr_inelgt0->GetPoint(i, x_run2, yield_run2);
        x_error = gr_inelgt0->GetErrorX(i);
        y_error_run2 = gr_inelgt0->GetErrorY(i);

        // double bincontent = h1->GetBinContent(i + 1);
        // double binerror = h1->GetBinError(i + 1);
        // double ratio = bincontent / yield_run2;
        // double error = sqrt(pow(binerror / yield_run2, 2) + pow(bincontent * y_error_run2 / (yield_run2 * yield_run2), 2));
        // gratio->SetPoint(i, x_run2, ratio);
        // gratio->SetPointError(i, x_error, error);

        double thisanalysis = fitFcn->Eval(x_run2);
        cout << "run 2 x is " << x_run2 << " y run2 is " << yield_run2 << "run 3 " << thisanalysis << endl;
        gratio->SetPoint(i, x_run2, thisanalysis / yield_run2);
        double error = sqrt(pow(thisanalysis * y_error_run2 / (yield_run2 * yield_run2), 2));
        gratio->SetPointError(i, x_error, error);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    SetHistoStyle(h1, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    h1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    h1->SetMaximum(h1->GetMaximum() * 1.5);
    h1->SetMinimum(h1->GetMinimum() * 10);
    h1->GetYaxis()->SetTitleOffset(1.15);
    h1->GetXaxis()->SetTitleOffset(1.02);
    h1->SetMarkerStyle(22);
    h1->SetMarkerSize(1);
    h1->GetXaxis()->SetRangeUser(0, 10);
    h1->Draw("pe");
    gPad->SetLogy(1);
    gr_inelgt0->SetMarkerStyle(29);
    gr_inelgt0->SetMarkerSize(1);
    gr_inelgt0->SetMarkerColor(kRed);
    gr_inelgt0->SetLineColor(kRed);
    gr_inelgt0->Draw("pe same");

    TLegend *leg = new TLegend(0.46, 0.7, 0.9, 0.91);
    SetLegendStyle(leg);
    // leg->AddEntry(h1, "pp 13.6 TeV (TPC + TOF)", "lpe");
    leg->AddEntry(h1, "pp 13.6 TeV (TPC only)", "lpe");
    leg->AddEntry(fitFcn, "Levy-Tsallis fit (pp 13.6 TeV)", "l");
    leg->AddEntry(gr_inelgt0, "pp 13 TeV (Published)", "lpe");
    leg->SetTextSize(0.05);
    leg->Draw();

    c1->cd(2);
    TH1F *hdummy = (TH1F *)h1->Clone();
    for (int i = 0; i < hdummy->GetNbinsX(); i++)
    {
        hdummy->SetBinContent(i + 1, 0);
        hdummy->SetBinError(i + 1, 0);
    }

    SetgrgrStyle(gratio, 1, 53, 1, 0.05, 0.05, 0.04 / pad2Size, 0.04 / pad2Size, 1.13, 1.8);
    gratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    gratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    gratio->SetMarkerStyle(23);
    gratio->SetMarkerSize(1);
    gratio->SetMarkerColor(kRed);
    gratio->SetLineColor(kRed);
    gratio->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
    gratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gratio->GetXaxis()->CenterTitle(1);
    gratio->GetYaxis()->SetTitleOffset(0.45);
    gratio->GetYaxis()->SetNdivisions(506);
    gratio->GetXaxis()->SetRangeUser(0, 10);
    gratio->Draw("ap");

    TLine *line = new TLine(0, 1, 15, 1);
    line->SetLineStyle(2);
    line->Draw();

    // TCanvas *ctemp = new TCanvas("ctemp", "ctemp", 850, 900);
    // h1->Draw();
    // fitFcn->SetLineColor(8);
    // fitFcn->Draw("same");

    // c1->SaveAs("/home/sawan/TPC_only.png");
}