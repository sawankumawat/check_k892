#include <iostream>
#include "src/style.h"

using namespace std;
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void compare_YieldMethod()
{
    TFile *fYield = new TFile("../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/yield.root", "read");
    if (fYield->IsZombie())
    {
        cout << "File not found" << endl;
        return;
    }
    TH1F *hYieldFunctionIntegral = (TH1F *)fYield->Get("mult_0-100/yield_integral");
    TH1F *hYieldBinCount = (TH1F *)fYield->Get("mult_0-100/yield_bincount");
    if (hYieldFunctionIntegral == nullptr || hYieldBinCount == nullptr)
    {
        cout << "Histogram not found" << endl;
        return;
    }

    TCanvas *cPlot = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cPlot, 0.15, 0.005, 0.08, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cPlot, pad1Size, pad2Size);
    cPlot->cd(1);
    gPad->SetLogy(1);
    SetHistoQA(hYieldFunctionIntegral);
    hYieldFunctionIntegral->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldFunctionIntegral->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldFunctionIntegral->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hYieldFunctionIntegral->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hYieldFunctionIntegral->GetYaxis()->SetTitleOffset(1.7 * pad1Size);
    hYieldFunctionIntegral->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hYieldFunctionIntegral->GetYaxis()->SetTitle("Raw Yield");
    hYieldFunctionIntegral->SetMaximum(hYieldFunctionIntegral->GetMaximum() * 1.6);
    hYieldFunctionIntegral->Draw("ep");
    hYieldBinCount->SetMarkerStyle(25);
    hYieldBinCount->SetMarkerColor(kBlue);
    hYieldBinCount->SetLineColor(kBlue);
    hYieldBinCount->Draw("ep same");
    TLegend *leg = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hYieldFunctionIntegral, "Function integration", "p");
    leg->AddEntry(hYieldBinCount, "Bin counting", "p");
    leg->Draw();

    cPlot->cd(2);
    gPad->SetGridy(1);
    TH1F *hRatio1 = (TH1F *)hYieldFunctionIntegral->Clone("hRatio1");
    SetHistoQA(hRatio1);
    hRatio1->Divide(hYieldBinCount);
    hRatio1->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hRatio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hRatio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hRatio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio1->GetYaxis()->SetTitle("Ratio");
    hRatio1->GetYaxis()->SetTitleOffset(1.7 * pad2Size);
    hRatio1->SetMarkerStyle(20);
    hRatio1->SetMarkerSize(1.0);
    hRatio1->GetYaxis()->SetNdivisions(505);
    hRatio1->GetYaxis()->SetRangeUser(0.89, 1.09);
    for (int i = 1; i <= hRatio1->GetNbinsX(); i++)
    {
        hRatio1->SetBinError(i, 0.0);
    }
    hRatio1->Draw("p");
    TLine *line = new TLine(0, 1, 10, 1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw();
    cPlot->SaveAs("../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/compare_yield_methods.pdf");
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