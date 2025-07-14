#include <iostream>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/common.h"
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
    pad1->SetLeftMargin(0.15);
    pad2->SetLeftMargin(0.15);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.04);
}

void yield_compare_inelgt0()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TGraphErrors *gr = new TGraphErrors();
    int multi_case = 2; // 0: 0-1, 1: 1-5, 2: 5-10% multiplicity
    if (multi_case == 0)
    {
        TFile *fpub = new TFile("HEPData-ins1762348-v1-Table_1.root", "READ");
        TH1* y = (TH1*)fpub->Get("Table 1/Hist1D_y1");
        TH1* ye1 = (TH1*)fpub->Get("Table 1/Hist1D_y1_e1");
        TH1* ye2 = (TH1*)fpub->Get("Table 1/Hist1D_y1_e2");
        // TGraphErrors *gr = new TGraphErrors(y->GetNbinsX());
        for (int i = 0; i < y->GetNbinsX(); i++)
        {
            double x = y->GetBinCenter(i + 1);
            double yval = y->GetBinContent(i + 1);
            // double yerr = sqrt(pow(ye1->GetBinContent(i + 1), 2) + pow(ye2->GetBinContent(i + 1), 2));
            double yerr = ye2->GetBinContent(i + 1); // consider only the systematic error at the moment.
            gr->SetPoint(i, x, yval);
            gr->SetPointError(i, y->GetBinWidth(i + 1) / 2, yerr);
        }
        TFile *fpub2 = new TFile("HEPData-ins1762348-v1-Table_10.root", "READ"); // ratio
        TH1* y2 = (TH1*)fpub2->Get("Table 10/Hist1D_y1");
        TH1* ye1_2 = (TH1*)fpub2->Get("Table 10/Hist1D_y1_e1");
        TH1* ye2_2 = (TH1*)fpub2->Get("Table 10/Hist1D_y1_e2");
        // Divide first graph by the second
        for (int i = 0; i < gr->GetN(); i++)
        {
            double x, y;
            gr->GetPoint(i, x, y);
            double yerr = gr->GetErrorY(i);
            double y2val = y2->GetBinContent(i + 1);
            double y2err = ye2_2->GetBinContent(i + 1);
            gr->SetPoint(i, x, y / y2val);
            gr->SetPointError(i, gr->GetErrorX(i), y / y2val * sqrt(pow(yerr / y, 2) + pow(y2err / y2val, 2)));
        }
    }
    if (multi_case == 1)
    {
        TFile *fpub = new TFile("HEPData-ins1762348-v1-Table_2.root", "READ");
        TH1* y = (TH1*)fpub->Get("Table 2/Hist1D_y1");
        TH1* ye1 = (TH1*)fpub->Get("Table 2/Hist1D_y1_e1");
        TH1* ye2 = (TH1*)fpub->Get("Table 2/Hist1D_y1_e2");
        // TGraphErrors *gr = new TGraphErrors(y->GetNbinsX());
        for (int i = 0; i < y->GetNbinsX(); i++)
        {
            double x = y->GetBinCenter(i + 1);
            double yval = y->GetBinContent(i + 1);
            // double yerr = sqrt(pow(ye1->GetBinContent(i + 1), 2) + pow(ye2->GetBinContent(i + 1), 2));
            double yerr = ye2->GetBinContent(i + 1); // consider only the systematic error at the moment.
            gr->SetPoint(i, x, yval);
            gr->SetPointError(i, y->GetBinWidth(i + 1) / 2, yerr);
        }
        TFile *fpub2 = new TFile("HEPData-ins1762348-v1-Table_11.root", "READ"); // ratio
        TH1* y2 = (TH1*)fpub2->Get("Table 11/Hist1D_y1");
        TH1* ye1_2 = (TH1*)fpub2->Get("Table 11/Hist1D_y1_e1");
        TH1* ye2_2 = (TH1*)fpub2->Get("Table 11/Hist1D_y1_e2");
        // Divide first graph by the second
        for (int i = 0; i < gr->GetN(); i++)
        {
            double x, y;
            gr->GetPoint(i, x, y);
            double yerr = gr->GetErrorY(i);
            double y2val = y2->GetBinContent(i + 1);
            double y2err = ye2_2->GetBinContent(i + 1);
            gr->SetPoint(i, x, y / y2val);
            gr->SetPointError(i, gr->GetErrorX(i), y / y2val * sqrt(pow(yerr / y, 2) + pow(y2err / y2val, 2)));
        }
    }
    if (multi_case == 2)
    {
        TFile *fpub = new TFile("HEPData-ins1762348-v1-Table_3.root", "READ");
        TH1* y = (TH1*)fpub->Get("Table 3/Hist1D_y1");
        TH1* ye1 = (TH1*)fpub->Get("Table 3/Hist1D_y1_e1");
        TH1* ye2 = (TH1*)fpub->Get("Table 3/Hist1D_y1_e2");
        // TGraphErrors *gr = new TGraphErrors(y->GetNbinsX());
        for (int i = 0; i < y->GetNbinsX(); i++)
        {
            double x = y->GetBinCenter(i + 1);
            double yval = y->GetBinContent(i + 1);
            // double yerr = sqrt(pow(ye1->GetBinContent(i + 1), 2) + pow(ye2->GetBinContent(i + 1), 2));
            double yerr = ye2->GetBinContent(i + 1); // consider only the systematic error at the moment.
            gr->SetPoint(i, x, yval);
            gr->SetPointError(i, y->GetBinWidth(i + 1) / 2, yerr);
        }
        TFile *fpub2 = new TFile("HEPData-ins1762348-v1-Table_12.root", "READ"); // ratio
        TH1* y2 = (TH1*)fpub2->Get("Table 12/Hist1D_y1");
        TH1* ye1_2 = (TH1*)fpub2->Get("Table 12/Hist1D_y1_e1");
        TH1* ye2_2 = (TH1*)fpub2->Get("Table 12/Hist1D_y1_e2");
        // Divide first graph by the second
        for (int i = 0; i < gr->GetN(); i++)
        {
            double x, y;
            gr->GetPoint(i, x, y);
            double yerr = gr->GetErrorY(i);
            double y2val = y2->GetBinContent(i + 1);
            double y2err = ye2_2->GetBinContent(i + 1);
            gr->SetPoint(i, x, y / y2val);
            gr->SetPointError(i, gr->GetErrorX(i), y / y2val * sqrt(pow(yerr / y, 2) + pow(y2err / y2val, 2)));
        }
    }

    TFile *spectra = new TFile("/Users/blim/cernbox/SWAN_projects/resonances/K892/output/LHC22o_pass6_small_TPC2TOF3Cut_tune/common/spectra_1_0.root", "READ");
    TH1D *h1 = (TH1D *)spectra->Get("lf-k892analysis/K892/0/hCorrectedYields");
    if (!h1)
    {
        cout << "Histogram not found" << endl;
        return;
    }

    TF1 *fLevyTsallis = new TF1("levyTsallis", levyTsallis, 0, 16, 3); // Adjust range and parameter count as needed
    // Set initial parameter values
    fLevyTsallis->SetParameters(0.155025, 6.91579, 0.219293); // Example values: q, T, N, mass
    fLevyTsallis->SetParNames("dN/dy", "T", "n");
    fLevyTsallis->SetLineWidth(2);
    fLevyTsallis->SetLineColor(kBlue);

    // Fit the function to the histogram
    h1->Fit("levyTsallis", "RM"); // "R" option to use the function range
    TCanvas *c = new TCanvas("c", "c", 850, 900); // 850, 900
    c->SetLogy(1);
    h1->GetXaxis()->SetRangeUser(0, 15);
    h1->Draw();

    TGraphErrors *gratio = new TGraphErrors();
    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        double published = y;
        double thisanalysis = fLevyTsallis->Eval(x);
        gratio->SetPoint(i, x, thisanalysis / published);
        double error = thisanalysis / published * sqrt(2*pow(gr->GetErrorY(i) / y, 2));
        gratio->SetPointError(i, gr->GetErrorX(i), error);
    }
    gratio->SetPoint(gr->GetN(), h1->GetBinCenter(h1->GetNbinsX()), 999);

    // TCanvas *c2 = new TCanvas();
    // gratio->SetMarkerStyle(20);
    // gratio->SetMarkerSize(1.5);
    // gratio->Draw("ap");

    TCanvas *c1 = new TCanvas("c1", "c1", 850, 900);
    SetCanvasStyle(c1, 0.16, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    SetHistoStyle(h1, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    h1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    h1->SetMaximum(5e-1);
    h1->SetMinimum(7e-7);
    h1->GetYaxis()->SetTitleOffset(1.02);
    h1->GetXaxis()->SetTitleOffset(1.02);
    h1->SetMarkerStyle(22);
    h1->SetMarkerSize(1.5);
    h1->GetYaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetTitleOffset(1.7);
    h1->GetYaxis()->SetLabelSize(0.045);
    // h1->Scale(1.5);
    h1->Draw();
    c1->SetLogy(1);
    gPad->SetLogy(1);
    gr->SetMarkerStyle(29);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed);
    gr->Draw("psame");

    TLegend *leg = new TLegend(0.4, 0.75, 0.98, 0.98);
    SetLegendStyle(leg);
    leg->AddEntry(h1, "pp 13.6 TeV (This Analysis)", "lpe");
    leg->AddEntry(fLevyTsallis, "Levy-Tsallis Fit (pp 13.6 TeV)", "l");
    leg->AddEntry(gr, "pp 13 TeV (Published)", "lpe");
    leg->SetTextSize(0.05);
    leg->Draw();

    c1->cd(2);
    SetgrgrStyle(gratio, 1, 53, 1, 0.05, 0.05, 0.04 / pad2Size, 0.04 / pad2Size, 1.13, 1.8);
    gratio->SetMinimum(0.4);
    gratio->GetYaxis()->SetRangeUser(0.4, 1.5);
    gratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    gratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    gratio->SetMarkerStyle(23);
    gratio->SetMarkerSize(1.5);
    gratio->SetMarkerColor(kRed);
    gratio->SetLineColor(kRed);
    gratio->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
    gratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gratio->GetXaxis()->CenterTitle(1);
    // gratio->GetYaxis()->CenterTitle(1);
    gratio->GetYaxis()->SetTitleOffset(0.4);
    gratio->GetYaxis()->SetNdivisions(505);
    gratio->GetYaxis()->SetTitleSize(0.1);
    gratio->GetYaxis()->SetTitleOffset(0.7);
    gratio->GetYaxis()->SetLabelSize(0.1);
    gratio->Draw("ap");

    // Add a line at y=1
    TLine *line = new TLine(0, 1, 16.5, 1);
    line->SetLineStyle(2);
    line->Draw("same");

    c1->SaveAs("yield_compare_INELgt0_TOF_TPC_tune.png");
    c1->SaveAs("yield_compare_INELgt0_TOF_TPC_tune.pdf");
}