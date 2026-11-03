#include <iostream>
#include "src/style.h"

using namespace std;
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void plotEfficiency()
{
    TString pathINELgt0 = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED";
    TString pathINEL = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/756343/kstarqa/hInvMass/ROTATED";

    TFile *fileINEL = new TFile((pathINEL + "/corrected_spectra_0_120.root").Data(), "read");
    TH1D *hEffINEL = (TH1D *)fileINEL->Get("mult_0-120/heff");

    TH1D *hEffINELgt0[11];

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

    const int kColorsMult[12] = {
        kRed + 2,
        kRed - 4,
        kOrange - 3,
        kOrange - 2,
        kSpring - 1,
        kGreen + 2,
        kCyan + 1,
        kAzure + 1,
        kAzure + 2,
        kBlue + 2,
        kViolet + 1,
        kGreen + 4};

    int markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 32, 47};

    for (int imult = 0; imult < 11; imult++)
    {

        float multlow, multhigh;
        if (imult == 0)
        {
            multlow = 0;
            multhigh = 100; // for all multiplicity
        }
        else
        {
            multlow = mult_classes[imult - 1];
            multhigh = mult_classes[imult];
        }

        TFile *fileINELgt0 = new TFile((pathINELgt0 + Form("/corrected_spectra_%.0f_%.0f.root", multlow, multhigh)).Data(), "read");
        hEffINELgt0[imult] = (TH1D *)fileINELgt0->Get(Form("mult_%.0f-%.0f/heff", multlow, multhigh));
    }

    TCanvas *cefficiency = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cefficiency, 0.16, 0.06, 0.01, 0.14);
    double pad1Size, pad2Size;
    canvas_style(cefficiency, pad1Size, pad2Size);
    cefficiency->cd(1);
    for (int imult = 0; imult < 11; imult++)
    {
        SetHistoQA(hEffINELgt0[imult]);
        hEffINELgt0[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hEffINELgt0[imult]->GetYaxis()->SetTitle("Acceptance x Efficiency");
        hEffINELgt0[imult]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
        hEffINELgt0[imult]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        hEffINELgt0[imult]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
        hEffINELgt0[imult]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
        hEffINELgt0[imult]->GetYaxis()->SetTitleOffset(1.7 * pad1Size);
        hEffINELgt0[imult]->SetMaximum(0.95);
        hEffINELgt0[imult]->SetMarkerStyle(markers[imult]);
        hEffINELgt0[imult]->SetMarkerSize(1.2);
        hEffINELgt0[imult]->SetLineColor(kColorsMult[imult]);
        hEffINELgt0[imult]->SetMarkerColor(kColorsMult[imult]);
        hEffINELgt0[imult]->Draw("pe same");
    }
    SetHistoQA(hEffINEL);
    hEffINEL->SetMarkerStyle(49);
    hEffINEL->SetMarkerSize(1.2);
    hEffINEL->SetLineColor(kColorsMult[11]);
    hEffINEL->SetMarkerColor(kColorsMult[11]);
    hEffINEL->Draw("pe same");

    TLegend *legall = new TLegend(0.5, 0.65, 0.92, 0.92);
    legall->SetTextSize(0.035);
    legall->SetNColumns(3);
    legall->SetFillStyle(0);
    legall->SetBorderSize(0);
    legall->AddEntry(hEffINEL, "INEL", "p");
    legall->AddEntry(hEffINELgt0[0], "0-100%", "p");
    for (int imult = 1; imult < 11; imult++)
    {
        legall->AddEntry(hEffINELgt0[imult], Form("%.0f-%.0f%%", mult_classes[imult - 1], mult_classes[imult]), "p");
    }
    legall->Draw();

    TLatex *latexMC = new TLatex();
    latexMC->SetTextSize(0.04);
    latexMC->SetTextFont(42);
    latexMC->SetNDC();
    latexMC->DrawLatex(0.23, 0.89, "ALICE");
    latexMC->DrawLatex(0.23, 0.83, "pp #sqrt{#it{s}} = 13.6 TeV");
    latexMC->DrawLatex(0.23, 0.77, "|y| < 0.5");
    latexMC->DrawLatex(0.23, 0.70, "K*(892)^{0}");

    cefficiency->cd(2);
    gPad->SetGridy(1);
    for (int imult = 1; imult < 11; imult++)
    {
        TH1F *hRatio1 = (TH1F *)hEffINELgt0[imult]->Clone(Form("hRatio1_%d", imult));
        SetHistoQA(hRatio1);
        hRatio1->Divide(hEffINELgt0[0]);
        hRatio1->GetYaxis()->SetTitleSize(0.025 / pad2Size);
        hRatio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        hRatio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        hRatio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        hRatio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hRatio1->GetYaxis()->SetTitle("Ratio to min. bias");
        hRatio1->GetYaxis()->SetTitleOffset(2.5 * pad2Size);
        hRatio1->SetMarkerStyle(markers[imult]);
        hRatio1->SetMarkerSize(1.0);
        hRatio1->GetYaxis()->SetNdivisions(505);
        hRatio1->GetYaxis()->SetRangeUser(0.69, 1.29);
        hRatio1->SetLineColor(kColorsMult[imult]);
        hRatio1->SetMarkerColor(kColorsMult[imult]);
        hRatio1->Draw("ep same");
    }

    cefficiency->SaveAs(pathINELgt0 + "/efficiencyWithINEL.png");
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.08, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}