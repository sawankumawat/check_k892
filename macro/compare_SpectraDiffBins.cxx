#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"
#include "spectra/YieldMean.C"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size, double &pad3Size);

Double_t FuncLavy(Double_t *x, Double_t *par)
{
    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void compare_SpectraDiffBins()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    bool isSameBins = false;

    TString legend1 = "INEL > 0";
    TString legend2 = "INEL";

    string path1 = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED"; // INEL>0
    string path2 = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/756343/kstarqa/hInvMass/ROTATED"; // INEL
    TString outputPath = path1 + "/spectra_compare";
    gSystem->mkdir(outputPath, kTRUE);

    TFile *fspectra1 = new TFile((path1 + "/corrected_spectra_0_100.root").c_str(), "read");
    TFile *fspectra2 = new TFile((path2 + "/corrected_spectra_0_120.root").c_str(), "read");

    TH1F *hSpectra1 = (TH1F *)fspectra1->Get("mult_0-100/corrected_spectra_Integral_final");
    TH1F *hSpectra2 = (TH1F *)fspectra2->Get("mult_0-120/corrected_spectra_Integral_final");

    // Fit the spectra with Levy function
    TF1 *fitFcn1 = new TF1("fitfunc", FuncLavy, 0.0, 20.0, 4);
    fitFcn1->SetParameter(0, 5.0);
    fitFcn1->SetParameter(1, 0.07);
    fitFcn1->FixParameter(2, 0.89556);
    fitFcn1->SetParameter(3, 0.3);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn1->SetLineColor(kRed + 1);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);

    Double_t min = 0;
    Double_t max = 20;
    Double_t loprecision = 0.001;
    Double_t hiprecision = 0.5;
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 20;


    double pT_bins[29 + 1] = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0}; // INEL with even higher pT bins

    TH1F *hClone1 = (TH1F *)hSpectra1->Clone("hClone1");
    TH1F *hClone2 = (TH1F *)hSpectra1->Clone("hClone2");

    for (int i = 1; i <= hClone1->GetNbinsX(); i++)
    {
        double systemerr1 = (hClone2->GetBinContent(i) * 0.08);
        hClone2->SetBinError(i, systemerr1);
    }

    TH1 *hout = YieldMean(hClone1, hClone2, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    // ============================================================
    // Ratio histograms -- use INEL>0 binning
    // ============================================================

    TH1F *hratioINEL = (TH1F *)hSpectra1->Clone("hratioINEL");
    TH1F *hratioDataFit = (TH1F *)hSpectra1->Clone("hratioDataFit");

    hratioINEL->Reset();
    hratioDataFit->Reset();

    // ============================================================
    // Loop over INEL>0 bins
    // ============================================================

    for (int ibin = 1; ibin <= hSpectra1->GetNbinsX(); ibin++)
    {
        double xLow = hSpectra1->GetBinLowEdge(ibin);
        double xHigh = xLow + hSpectra1->GetBinWidth(ibin);
        double xCenter = hSpectra1->GetBinCenter(ibin);

        double yData = hSpectra1->GetBinContent(ibin);
        double eData = hSpectra1->GetBinError(ibin);

        // ========================================================
        // 1. DATA / FIT
        // ========================================================

        double yFit = fitFcn1->Eval(xCenter);

        if (yFit > 0.0 && yData != 0.0)
        {
            double ratio = yData / yFit;
            double error = eData / std::abs(yFit);

            hratioDataFit->SetBinContent(ibin, ratio);
            hratioDataFit->SetBinError(ibin, error);
        }

        // ========================================================
        // 2. INEL > 0 / INEL
        //
        // Exact matching bin:
        //      use measured INEL point
        //
        // No matching bin:
        //      use INEL fit
        // ========================================================

        double yINEL = 0.0;
        double eINEL = 0.0;

        int matchingBin = -1;

        for (int jbin = 1; jbin <= hSpectra2->GetNbinsX(); jbin++)
        {
            double inelLow =
                hSpectra2->GetBinLowEdge(jbin);

            double inelHigh =
                inelLow + hSpectra2->GetBinWidth(jbin);

            if (std::abs(xLow - inelLow) < 1e-5 &&
                std::abs(xHigh - inelHigh) < 1e-5)
            {
                matchingBin = jbin;
                break;
            }
        }

        // --------------------------------------------------------
        // Use data if bins match
        // --------------------------------------------------------

        if (matchingBin >= 0)
        {
            yINEL = hSpectra2->GetBinContent(matchingBin);
            eINEL = hSpectra2->GetBinError(matchingBin);
        }

        // --------------------------------------------------------
        // Otherwise use INEL fit
        // --------------------------------------------------------

        else
        {
            yINEL = fitFcn1->Eval(xCenter);
            eINEL = 0.0;
        }

        // --------------------------------------------------------
        // Calculate INEL>0 / INEL
        // --------------------------------------------------------

        if (yINEL > 0.0 && yData != 0.0)
        {
            double ratio = yData / yINEL;

            double error;

            if (matchingBin >= 0 && eINEL > 0.0)
            {
                // Data / data
                error = ratio *
                        std::sqrt(
                            std::pow(eData / yData, 2) +
                            std::pow(eINEL / yINEL, 2));
            }
            else
            {
                // Data / fit
                error = eData / std::abs(yINEL);
            }

            hratioINEL->SetBinContent(ibin, ratio);
            hratioINEL->SetBinError(ibin, error);
        }
    }

    // Canvas & Layout
    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    double pad1Size, pad2Size, pad3Size;
    canvas_style(c1, pad1Size, pad2Size, pad3Size);

    // --- Panel 1: Spectra ---
    c1->cd(1);
    gPad->SetLogy(1);

    SetHistoQA(hSpectra1);

    hSpectra1->GetYaxis()->SetTitleSize(0.03 / pad1Size);
    hSpectra1->GetXaxis()->SetTitleSize(0.03 / pad1Size);
    hSpectra1->GetYaxis()->SetLabelSize(0.03 / pad1Size);
    hSpectra1->GetXaxis()->SetLabelSize(0.03 / pad1Size);

    hSpectra1->SetMaximum(hSpectra1->GetMaximum() * 3);
    hSpectra1->SetMinimum(hSpectra1->GetMinimum() * 0.1);

    hSpectra1->GetYaxis()->SetTitleOffset(1.20);

    hSpectra1->SetMarkerStyle(20);
    hSpectra1->SetLineColor(kRed + 1);
    hSpectra1->SetMarkerColor(kRed + 1);
    hSpectra1->Draw("pe");

    hSpectra2->SetMarkerColor(kBlue + 1);
    hSpectra2->SetLineColor(kBlue + 1);
    hSpectra2->SetMarkerStyle(21);
    hSpectra2->Draw("pe same");

    // INEL > 0 fit
    fitFcn1->SetLineColor(kRed + 1);
    fitFcn1->Draw("same");

    TLegend *leg = new TLegend(0.64, 0.73, 0.9, 0.94);
    SetLegendStyle(leg);
    leg->SetTextSize(0.045);

    leg->AddEntry(hSpectra1, Form("%s", legend1.Data()), "p");
    leg->AddEntry(hSpectra2, Form("%s", legend2.Data()), "p");
    leg->AddEntry(fitFcn1, "Levy-Tsallis (INEL > 0)", "l");

    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.SetTextFont(42);

    lat.DrawLatex(0.45, 0.90, "ALICE");
    lat.DrawLatex(0.45, 0.83, "|y| < 0.5");
    lat.DrawLatex(0.45, 0.76, "K*(892)^{0}");

    // ============================================================
    // Panel 2: Data / Fit
    // ============================================================

    c1->cd(2);

    SetHistoQA(hratioDataFit);

    gPad->SetGrid(1, 1);

    hratioDataFit->GetYaxis()->SetTitleSize(0.025 / pad2Size);
    hratioDataFit->GetXaxis()->SetTitleSize(0.03 / pad2Size);

    hratioDataFit->GetYaxis()->SetLabelSize(0.03 / pad2Size);
    hratioDataFit->GetXaxis()->SetLabelSize(0.00);

    hratioDataFit->SetMarkerStyle(20);
    hratioDataFit->SetMarkerSize(1.0);

    hratioDataFit->SetMarkerColor(kRed + 1);
    hratioDataFit->SetLineColor(kRed + 1);

    hratioDataFit->GetYaxis()->SetTitle("Data / Fit");
    hratioDataFit->GetYaxis()->SetTitleOffset(0.6);
    hratioDataFit->GetYaxis()->SetNdivisions(506);

    hratioDataFit->SetMaximum(1.49);
    hratioDataFit->SetMinimum(0.51);

    hratioDataFit->Draw("pe");

    TLine *line1 = new TLine(0, 1, 20, 1);
    line1->SetLineStyle(2);
    line1->SetLineWidth(2);
    line1->SetLineColor(kBlack);
    line1->Draw();

    // ============================================================
    // Panel 3: INEL > 0 / INEL
    // ============================================================

    c1->cd(3);

    SetHistoQA(hratioINEL);

    gPad->SetGrid(1, 1);

    hratioINEL->GetYaxis()->SetTitleSize(0.025 / pad3Size);
    hratioINEL->GetXaxis()->SetTitleSize(0.03 / pad3Size);

    hratioINEL->GetYaxis()->SetLabelSize(0.03 / pad3Size);
    hratioINEL->GetXaxis()->SetLabelSize(0.03 / pad3Size);

    hratioINEL->SetMarkerStyle(21);
    hratioINEL->SetMarkerSize(1.0);

    hratioINEL->SetMarkerColor(kBlue + 1);
    hratioINEL->SetLineColor(kBlue + 1);

    hratioINEL->GetYaxis()->SetTitle("INEL > 0 / INEL");
    hratioINEL->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

    hratioINEL->GetYaxis()->SetTitleOffset(0.6);
    hratioINEL->GetXaxis()->SetTitleOffset(1.1);

    hratioINEL->GetYaxis()->SetNdivisions(506);

    hratioINEL->SetMaximum(1.49);
    hratioINEL->SetMinimum(0.51);

    hratioINEL->Draw("pe");

    TLine *line2 = new TLine(0, 1, 20, 1);
    line2->SetLineStyle(2);
    line2->SetLineWidth(2);
    line2->SetLineColor(kBlack);
    line2->Draw();
    c1->SaveAs(outputPath + "/SpectraRatio.png");
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size, double &pad3Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 3, 0, 0);

    pad1Size = 0.50; // Top panel (Spectra)
    pad2Size = 0.25; // Middle panel (Data1 / Fit - Red)
    pad3Size = 0.25; // Bottom panel (Data2 / Fit - Blue)

    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    TPad *pad3 = (TPad *)c->GetPad(3);

    pad1->SetPad(0, 0.50, 1, 1.00);
    pad2->SetPad(0, 0.25, 1, 0.50);
    pad3->SetPad(0, 0.00, 1, 0.25);

    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad3->SetRightMargin(0.06);

    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad3->SetLeftMargin(0.16);

    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);

    pad2->SetTopMargin(0.001);
    pad2->SetBottomMargin(0.001);

    pad3->SetTopMargin(0.001);
    pad3->SetBottomMargin(0.33);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad3->SetTickx(1);
    pad3->SetTicky(1);
}