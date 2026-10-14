#include <iostream>
#include <cmath>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../spectra/YieldMean.C"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

Double_t FuncLavy(Double_t *x, Double_t *par)
{
    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
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

void SpectraRatioKshortRun2_run3()
{
    gStyle->SetOptStat(0);
    TFile *fKshortRun2 = OpenFile("K0s_Run2_Results/K0s_Run2_Spectra.root");
    TFile *fKshort = OpenFile("K0s_Run3_Results/Spectra_K0Short.root");

    if (fKshortRun2->IsZombie() || fKshort->IsZombie())
    {
        cout << "Error: One or both files could not be opened." << endl;
        return;
    }

    TH1F *hK0sRun2 = (TH1F *)fKshortRun2->Get("Table 1/Hist1D_y10");
    TH1F *hK0sRun2_stat = (TH1F *)fKshortRun2->Get("Table 1/Hist1D_y10_e1");
    TH1F *hK0sRun2_sys = (TH1F *)fKshortRun2->Get("Table 1/Hist1D_y10_e2");

    for (int ibin = 1; ibin <= hK0sRun2->GetNbinsX(); ++ibin)
    {
        double binContent = hK0sRun2->GetBinContent(ibin);
        double binErrorStat = hK0sRun2_stat->GetBinContent(ibin);
        double binErrorSys = hK0sRun2_sys->GetBinContent(ibin);

        hK0sRun2_stat->SetBinContent(ibin, binContent);
        hK0sRun2_stat->SetBinError(ibin, binErrorStat);

        hK0sRun2_sys->SetBinContent(ibin, binContent);
        hK0sRun2_sys->SetBinError(ibin, binErrorSys);
    }

    // TH1F *hKshort_stat = (TH1F *)fKshort->Get("hK0sMB");
    // TH1F *hKshort_sys = (TH1F *)fKshort->Get("hK0sMBSyst");

    TH1F *hKshort_stat = (TH1F *)fKshort->Get("hK0s_70_100");
    TH1F *hKshort_sys = (TH1F *)fKshort->Get("hK0sSyst_70_100");

    TH1F *hKstarK0sRatio_stat = (TH1F *)hK0sRun2_stat->Clone("hKstarK0sRatio_stat");
    TH1F *hKstarK0sRatio_sys = (TH1F *)hK0sRun2_sys->Clone("hKstarK0sRatio_sys");

    if (hK0sRun2_stat == nullptr || hK0sRun2_sys == nullptr || hKshort_stat == nullptr || hKshort_sys == nullptr)
    {
        cout << "Error: One or more histograms could not be retrieved." << endl;
        return;
    }

    TF1 *fitFcn1 = new TF1("fitfunc", FuncLavy, 0.0, 10.0, 4);
    fitFcn1->SetParameter(0, 6.6);
    fitFcn1->SetParameter(1, 0.2);
    fitFcn1->SetParameter(2, 0.493677);
    fitFcn1->SetParameter(3, 0.18);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn1->SetLineColor(kRed + 1);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);

    Double_t min = 0.0;
    Double_t max = 10;
    Double_t loprecision = 0.001;
    Double_t hiprecision = 0.5;
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 10;

    TH1 *hout = YieldMean(hKshort_stat, hKshort_sys, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    auto getKaonErr = [&](TH1F *hErr, double xLow, double xHigh) -> double
    {
        const double eps = 1e-6;
        const int nBins = hErr->GetNbinsX();

        // Kaon histogram range
        const double kaonMin = hErr->GetBinLowEdge(1);
        const double kaonMax = hErr->GetBinLowEdge(nBins) + hErr->GetBinWidth(nBins);

        // ---------------------------------------------------------
        // Rule 3:
        // If the K* bin is outside the kaon histogram range,
        // use the error from the last kaon bin.
        // ---------------------------------------------------------
        if (xLow < kaonMin - eps || xHigh > kaonMax + eps)
        {
            return hErr->GetBinError(nBins);
        }

        // ---------------------------------------------------------
        // Rule 1:
        // Check whether K* bin edges exactly match ONE kaon bin.
        // ---------------------------------------------------------
        for (int ibin = 1; ibin <= nBins; ++ibin)
        {
            double kaonLow = hErr->GetBinLowEdge(ibin);
            double kaonHigh = kaonLow + hErr->GetBinWidth(ibin);

            if (std::abs(xLow - kaonLow) < eps &&
                std::abs(xHigh - kaonHigh) < eps)
            {
                return hErr->GetBinError(ibin);
            }
        }

        // ---------------------------------------------------------
        // Rule 2:
        // K* bin edges do not exactly match.
        // Find all kaon bins that overlap the K* bin and
        // average their uncertainties.
        // ---------------------------------------------------------
        double sumErr = 0.0;
        int nOverlap = 0;

        for (int ibin = 1; ibin <= nBins; ++ibin)
        {
            double kaonLow = hErr->GetBinLowEdge(ibin);
            double kaonHigh = kaonLow + hErr->GetBinWidth(ibin);

            // Check for overlap
            if (kaonHigh > xLow + eps && kaonLow < xHigh - eps)
            {
                sumErr += hErr->GetBinError(ibin);
                ++nOverlap;
            }
        }

        if (nOverlap > 0)
        {
            return sumErr / nOverlap;
        }

        // ---------------------------------------------------------
        // Safety fallback
        // ---------------------------------------------------------
        return hErr->GetBinError(nBins);
    };

    auto getKaonYield = [&](TH1F *hKaon, TF1 *fitFcn,
                            double xLow, double xHigh) -> double
    {
        const double eps = 1e-6;
        const int nBins = hKaon->GetNbinsX();

        // ---------------------------------------------------------
        // Check whether K* bin edges exactly match ONE K0S bin
        // ---------------------------------------------------------
        for (int ibin = 1; ibin <= nBins; ++ibin)
        {
            double kaonLow = hKaon->GetBinLowEdge(ibin);
            double kaonHigh = kaonLow + hKaon->GetBinWidth(ibin);

            if (std::abs(xLow - kaonLow) < eps &&
                std::abs(xHigh - kaonHigh) < eps)
            {
                // Exact bin match:
                // use the measured K0S value directly
                return hKaon->GetBinContent(ibin);
            }
        }

        // ---------------------------------------------------------
        // No exact bin match:
        // use the fitted value at the K* bin center
        // ---------------------------------------------------------
        double xCenter = 0.5 * (xLow + xHigh);

        return fitFcn->Eval(xCenter);
    };

    int TotalBinsKstar = hK0sRun2_stat->GetNbinsX();

    for (int ibin = 1; ibin <= TotalBinsKstar; ++ibin)
    {
        // ---------------------------------------------------------
        // K* value and uncertainties
        // ---------------------------------------------------------
        double kstarYield = hK0sRun2_stat->GetBinContent(ibin);
        double kstar_stat_err = hK0sRun2_stat->GetBinError(ibin);
        double kstar_sys_err = hK0sRun2_sys->GetBinError(ibin); // total systematic error

        // ---------------------------------------------------------
        // K* bin boundaries
        // ---------------------------------------------------------
        double xLow = hK0sRun2_stat->GetBinLowEdge(ibin);
        double xHigh = xLow + hK0sRun2_stat->GetBinWidth(ibin);
        double xCenter = hK0sRun2_stat->GetBinCenter(ibin);

        // ---------------------------------------------------------
        // Kaon value:
        // If K* bin exactly matches a K0S bin, use the measured
        // K0S bin content directly. Otherwise, use the fit value.
        // ---------------------------------------------------------
        double kaonYield = getKaonYield(hKshort_stat, fitFcn1, xLow, xHigh);

        // ---------------------------------------------------------
        // Protect against invalid values
        // ---------------------------------------------------------
        if (kstarYield <= 0 || kaonYield <= 0)
        {
            hKstarK0sRatio_stat->SetBinContent(ibin, 0);
            hKstarK0sRatio_stat->SetBinError(ibin, 0);

            hKstarK0sRatio_sys->SetBinContent(ibin, 0);
            hKstarK0sRatio_sys->SetBinError(ibin, 0);

            continue;
        }

        // ---------------------------------------------------------
        // Ratio
        // ---------------------------------------------------------
        double ratio = kstarYield / kaonYield;
        // cout << "Bin " << ibin << ", bin value " << xCenter << ", Kstar yield: " << kstarYield << ", Kshort yield : " << kaonYield << endl;

        // =========================================================
        // 1. STATISTICAL UNCERTAINTY
        // =========================================================

        double kaon_stat_err = getKaonErr(hKshort_stat, xLow, xHigh);
        double rel_kstar_stat = kstar_stat_err / kstarYield;
        double rel_kaon_stat = kaon_stat_err / kaonYield;

        double ratio_stat_err = ratio * std::sqrt(rel_kstar_stat * rel_kstar_stat + rel_kaon_stat * rel_kaon_stat);

        // =========================================================
        // 2. SYSTEMATIC UNCERTAINTY
        //
        // K*   -> total systematic
        // KaonYield -> total systematic
        // =========================================================

        double kaon_sys_err = getKaonErr(hKshort_sys, xLow, xHigh);
        double rel_kstar_sys = kstar_sys_err / kstarYield;
        double rel_kaon_sys = kaon_sys_err / kaonYield;

        double ratio_sys_err = ratio * std::sqrt(rel_kstar_sys * rel_kstar_sys + rel_kaon_sys * rel_kaon_sys);

        hKstarK0sRatio_stat->SetBinContent(ibin, ratio);
        hKstarK0sRatio_stat->SetBinError(ibin, ratio_stat_err);

        hKstarK0sRatio_sys->SetBinContent(ibin, ratio);
        hKstarK0sRatio_sys->SetBinError(ibin, ratio_sys_err);
    }

    TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 720, 720);
    SetCanvasStyle(cRatio, 0.16, 0.06, 0.06, 0.13);
    double pad1Size, pad2Size;
    canvas_style(cRatio, pad1Size, pad2Size);
    cRatio->cd(1);
    gPad->SetLogy();
    SetHistoQA(hK0sRun2_stat);
    SetHistoQA(hKshort_stat);
    hK0sRun2_stat->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hK0sRun2_stat->GetYaxis()->SetTitle("1/N_{ev} d^{2}N/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}");
    hK0sRun2_stat->GetYaxis()->SetTitleOffset(1.60 * pad1Size);
    hK0sRun2_stat->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hK0sRun2_stat->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hK0sRun2_stat->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hK0sRun2_stat->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    // hK0sRun2_stat->GetYaxis()->SetRangeUser(0.08, 1.25);
    hK0sRun2_stat->GetXaxis()->SetRangeUser(0.0, 10.0);
    hK0sRun2_stat->SetMarkerStyle(20);
    hK0sRun2_stat->SetMarkerSize(1.25);
    hK0sRun2_stat->SetMarkerColor(kRed + 1);
    hK0sRun2_stat->SetLineColor(kRed + 1);
    // hK0sRun2_stat->Draw("pe");
    hK0sRun2_sys->SetFillStyle(0);
    hK0sRun2_sys->SetLineColor(kRed + 1);
    hK0sRun2_sys->SetMarkerColor(kRed + 1);
    // hK0sRun2_sys->Draw("e2 same");
    hKshort_stat->SetMarkerStyle(22);
    hKshort_stat->SetMarkerSize(1.25);
    hKshort_stat->SetMarkerColor(kBlue + 1);
    hKshort_stat->SetLineColor(kBlue + 1);
    // hKshort_stat->Draw("pe");
    hKshort_sys->SetFillStyle(0);
    hKshort_sys->SetLineColor(kBlue + 1);
    hKshort_sys->SetMarkerColor(kBlue + 1);
    // hKshort_sys->Draw("e2 same");

    hKshort_stat->Draw("pe");
    hKshort_sys->Draw("e2 same");

    hK0sRun2_sys->Draw("e2 same");
    hK0sRun2_stat->Draw("pe same");

    TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetHeader("Minimum Bias (70 - 100%)");
    leg->AddEntry(hK0sRun2_stat, "K^{*0} Run 2", "pe");
    leg->AddEntry(hKshort_stat, "K^{0}_{S} Run 3", "pe");
    leg->Draw("same");

    cRatio->cd(2);
    gPad->SetGridy();
    SetHistoQA(hKstarK0sRatio_stat);
    hKstarK0sRatio_stat->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hKstarK0sRatio_stat->GetYaxis()->SetTitle("Run 2  / Run 3");
    hKstarK0sRatio_stat->GetYaxis()->SetTitleOffset(1.60 * pad2Size);
    hKstarK0sRatio_stat->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hKstarK0sRatio_stat->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hKstarK0sRatio_stat->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hKstarK0sRatio_stat->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hKstarK0sRatio_stat->GetYaxis()->SetRangeUser(0.2, 1.08);
    hKstarK0sRatio_stat->GetXaxis()->SetRangeUser(0.0, 10.0);
    hKstarK0sRatio_stat->SetMarkerStyle(20);
    hKstarK0sRatio_stat->SetMarkerSize(1.25);
    hKstarK0sRatio_stat->SetMarkerColor(kRed + 1);
    hKstarK0sRatio_stat->SetLineColor(kRed + 1);
    hKstarK0sRatio_stat->GetYaxis()->SetNdivisions(507);
    hKstarK0sRatio_stat->Draw("pe");
    hKstarK0sRatio_sys->SetFillStyle(0);
    hKstarK0sRatio_sys->SetLineColor(kRed + 1);
    hKstarK0sRatio_sys->SetMarkerColor(kRed + 1);
    hKstarK0sRatio_sys->Draw("e2 same");

    TLine *line = new TLine(0, 1, 10, 1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("same");

    cRatio->SaveAs("Plots/pTDiff_Ratios/KshortRun2_vsRun3_Class70_100.png");
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