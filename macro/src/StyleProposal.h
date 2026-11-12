// ========================================================================
// Style
// ========================================================================

#ifndef STYLEPROPOSAL_H
#define STYLEPROPOSAL_H

#include <TColor.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TPad.h>
#include <TLegend.h>

std::vector<float> centralitybinning = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
std::vector<std::string> multiplicityClasses = {"0-1%", "1-5%", "5-10%", "10-15%", "15-20%", "20-30%", "30-40%", "40-50%", "50-70%", "70-100%"};

gStyle->SetOptStat(0);
gStyle->SetOptFit(0);

// Colors used for the multiplicity classes
const int kColorsMult[10] = {
    kRed + 2,
    kRed - 4,
    kOrange - 3,
    kOrange - 2,
    kSpring - 1,
    kGreen + 2,
    kCyan + 1,
    kAzure + 1,
    kAzure + 2,
    kBlue + 2};

// Colors commonly used for various plots
const int colors[] = {kMagenta - 3,
                      kBlue,
                      kCyan + 1,
                      kAzure + 1,
                      kGreen + 2,
                      kRed,
                      kOrange + 7,
                      kViolet + 1};

// ========================================================================
// Legend
// ========================================================================

void SetLegend(TLegend *legend)
{
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    legend->SetTextAlign(12);
}

// ========================================================================
// Plots with one pad (square canvas)
// ------------------------------------------------------------------------
// TCanvas *c = new TCanvas("c", "c", 700, 700);
// ========================================================================

void SetPad(TPad *pad)
{
    pad->SetBottomMargin(0.17);
    pad->SetTopMargin(0.05);
    pad->SetLeftMargin(0.17);
    pad->SetRightMargin(0.05);
    pad->SetTickx();
    pad->SetTicky();
    pad->Draw();
}

void SetPad(TVirtualPad *pad)
{
    pad->SetBottomMargin(0.17);
    pad->SetTopMargin(0.05);
    pad->SetLeftMargin(0.17);
    pad->SetRightMargin(0.05);
    pad->SetTickx();
    pad->SetTicky();
    pad->Draw();
}

// ========================================================================
// Histogram style
// ========================================================================

void SetHistStyle(TH1 *hist,
                  int color = kBlack,
                  double alpha = 1.0,
                  int markerStyle = 20,
                  double markerSize = 1.0,
                  int lineWidth = 1,
                  int lineStyle = 1)
{
    hist->SetLineColorAlpha(color, alpha);
    hist->SetMarkerColorAlpha(color, alpha);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(markerSize);
    hist->SetLineWidth(lineWidth);
    hist->SetLineStyle(lineStyle);
    hist->SetTitle("");

    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetXaxis()->SetTitleOffset(1.05);
    hist->GetXaxis()->SetLabelSize(0.045);
    hist->GetXaxis()->SetLabelOffset(0.01);

    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetTitleOffset(1.5);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetLabelOffset(0.01);

    hist->GetXaxis()->SetTickLength(0.025);
    hist->GetYaxis()->SetTickLength(0.025);

    hist->GetYaxis()->SetNdivisions(505);
}

// ========================================================================
// Function style
// ========================================================================

void SetFunctionStyle(TF1 *func,
                      int color = kBlack,
                      double alpha = 1.0,
                      int lineWidth = 2,
                      int lineStyle = 1)
{
    func->SetLineColorAlpha(color, alpha);
    func->SetLineWidth(lineWidth);
    func->SetLineStyle(lineStyle);
    func->SetFillStyle(0);
    func->SetFillColor(0);
    func->SetNpx(500);
}

// ========================================================================
// Style for ratio plots
// ------------------------------------------------------------------------
// Made for Canvas with two pads, top pad for the main plot
// and bottom pad for the ratio plot.
// TCanvas *c = new TCanvas("c", "c", 700, 850);
// TPad *top = new TPad("top", "top", 0, .30, 1, 1);
// SetTopPad(top);
// TPad *bot = new TPad("bot", "bot", 0, 0, 1, .30);
// SetBottomPad(bot);
// ========================================================================

void SetTopPad(TPad *padTop, bool logY = false)
{
    padTop->SetBottomMargin(0.00);
    padTop->SetTopMargin(0.06);
    padTop->SetLeftMargin(0.16);
    padTop->SetRightMargin(0.05);
    padTop->SetTickx();
    padTop->SetTicky();
    padTop->SetLogy(logY);
    padTop->Draw();
}

void SetBottomPad(TPad *padBottom, bool logY = false)
{
    padBottom->SetTopMargin(0.00);
    padBottom->SetBottomMargin(0.35);
    padBottom->SetLeftMargin(0.16);
    padBottom->SetRightMargin(0.05);
    padBottom->SetTickx();
    padBottom->SetTicky();
    padBottom->SetLogy(logY);
    padBottom->Draw();
}

void SetTopPadHistStyle(TH1 *hist,
                        int color = kBlack,
                        double alpha = 1.0,
                        int markerStyle = 20,
                        double markerSize = 1.0,
                        int lineWidth = 1,
                        int lineStyle = 1)
{
    hist->SetLineColorAlpha(color, alpha);
    hist->SetMarkerColorAlpha(color, alpha);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(markerSize);
    hist->SetLineWidth(lineWidth);
    hist->SetLineStyle(lineStyle);
    hist->SetTitle("");

    hist->GetYaxis()->SetTitleSize(0.065);
    hist->GetYaxis()->SetTitleOffset(1.0);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelOffset(0.007);

    hist->GetYaxis()->CenterTitle(false);
    hist->GetXaxis()->CenterTitle(false);

    hist->GetXaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetTickLength(0.022);

    hist->GetYaxis()->SetNdivisions(505);
}

void SetRatioHistStyle(TH1 *hist,
                       int color = kBlack,
                       double alpha = 1.0,
                       int markerStyle = 20,
                       double markerSize = 1.0,
                       int lineWidth = 1,
                       int lineStyle = 1)
{
    hist->SetLineColorAlpha(color, alpha);
    hist->SetMarkerColorAlpha(color, alpha);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(markerSize);
    hist->SetLineWidth(lineWidth);
    hist->SetLineStyle(lineStyle);
    hist->SetTitle("");

    hist->GetXaxis()->SetTitleSize(0.12);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.12);
    hist->GetXaxis()->SetLabelOffset(0.01);

    hist->GetYaxis()->SetTitleSize(0.12);
    hist->GetYaxis()->SetTitleOffset(0.55);
    hist->GetYaxis()->SetLabelSize(0.11);
    hist->GetYaxis()->SetLabelOffset(0.007);

    hist->GetXaxis()->SetTickLength(0.04);
    hist->GetYaxis()->SetTickLength(0.03);

    hist->GetYaxis()->SetNdivisions(505);
}

void DrawUnityLine(double xMin, double xMax, int lineStyle = 2, int lineWidth = 1, int lineColor = kBlack)
{
    TLine *line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineStyle(lineStyle);
    line->SetLineWidth(lineWidth);
    line->SetLineColor(lineColor);
    line->Draw("SAME");
}

// ========================================================================
// Graph style
// ========================================================================

void SetGraphStyle(TGraph *graph,
                   int color = kBlack,
                   double alpha = 1.0,
                   int markerStyle = 20,
                   double markerSize = 1.0,
                   int lineWidth = 1,
                   int lineStyle = 1)
{
    graph->SetLineColorAlpha(color, alpha);
    graph->SetMarkerColorAlpha(color, alpha);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetLineWidth(lineWidth);
    graph->SetLineStyle(lineStyle);
    graph->SetTitle("");

    graph->GetXaxis()->SetTitleSize(0.055);
    graph->GetXaxis()->SetTitleOffset(1.05);
    graph->GetXaxis()->SetLabelSize(0.045);
    graph->GetXaxis()->SetLabelOffset(0.01);

    graph->GetYaxis()->SetTitleSize(0.055);
    graph->GetYaxis()->SetTitleOffset(1.5);
    graph->GetYaxis()->SetLabelSize(0.045);
    graph->GetYaxis()->SetLabelOffset(0.01);

    graph->GetXaxis()->SetTickLength(0.025);
    graph->GetYaxis()->SetTickLength(0.025);

    graph->GetYaxis()->SetNdivisions(505);
}

void SetRatioGraphStyle(TGraph *graph,
                        int color = kBlack,
                        double alpha = 1.0,
                        int markerStyle = 20,
                        double markerSize = 1.0,
                        int lineWidth = 2,
                        int lineStyle = 1)
{
    graph->SetLineColorAlpha(color, alpha);
    graph->SetMarkerColorAlpha(color, alpha);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetLineWidth(lineWidth);
    graph->SetLineStyle(lineStyle);

    graph->SetTitle("");

    graph->GetXaxis()->SetTitleSize(0.12);
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetLabelSize(0.12);
    graph->GetXaxis()->SetLabelOffset(0.01);

    graph->GetYaxis()->SetTitleSize(0.12);
    graph->GetYaxis()->SetTitleOffset(0.55);
    graph->GetYaxis()->SetLabelSize(0.11);
    graph->GetYaxis()->SetLabelOffset(0.007);

    graph->GetXaxis()->SetTickLength(0.04);
    graph->GetYaxis()->SetTickLength(0.03);

    graph->GetYaxis()->SetNdivisions(505);
}

#endif