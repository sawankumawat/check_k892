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

void SpectraRatio_pTDiffKstarKshort()
{
    gStyle->SetOptStat(0);
    TFile *fKstar = OpenFile("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/KstarSpectra.root");
    TFile *fKshort = OpenFile("K0s_Run3_Results/Spectra_K0Short.root");

    TH1F *hKstar_statHM = (TH1F *)fKstar->Get("SpectraStat_0_1");
    TH1F *hKstar_sysHM = (TH1F *)fKstar->Get("SpectraSys_0_1");
    TH1F *hKstar_sysUncorrHM = (TH1F *)fKstar->Get("SpectraUncorrSys_0_1");

    TH1F *hKstar_statLM = (TH1F *)fKstar->Get("SpectraStat_70_100");
    TH1F *hKstar_sysLM = (TH1F *)fKstar->Get("SpectraSys_70_100");
    TH1F *hKstar_sysUncorrLM = (TH1F *)fKstar->Get("SpectraUncorrSys_70_100");

    TH1F *hKshort_statHM = (TH1F *)fKshort->Get("hK0s_0_1");
    TH1F *hKshort_sysHM = (TH1F *)fKshort->Get("hK0sSyst_0_1");
    TH1F *hKshort_sysUncorrHM = (TH1F *)fKshort->Get("hK0sSyst_0_1"); // Does not have uncorrelated systematic, so using total systematic

    TH1F *hKshort_statLM = (TH1F *)fKshort->Get("hK0s_70_100");
    TH1F *hKshort_sysLM = (TH1F *)fKshort->Get("hK0sSyst_70_100");
    TH1F *hKshort_sysUncorrLM = (TH1F *)fKshort->Get("hK0sSyst_70_100"); // Does not have uncorrelated systematic, so using total systematic

    cout << "Number of kaon bins in HM: " << hKshort_statHM->GetNbinsX() << endl;
    cout << "Number of kaon bins in LM: " << hKshort_statLM->GetNbinsX() << endl;

    int TotalBinsKaon = hKshort_statHM->GetNbinsX();

    TF1 *fitFcn1 = new TF1("fitfunc", FuncLavy, 0.0, 10.0, 4);
    fitFcn1->SetParameter(0, 6.6);
    fitFcn1->SetParameter(1, 0.2);
    fitFcn1->SetParameter(2, 0.493677);
    fitFcn1->SetParameter(3, 0.18);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn1->SetLineColor(kRed + 1);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);

    TF1 *fitFcn2 = new TF1("fitfunc2", FuncLavy, 0.0, 10.0, 4);
    fitFcn2->SetParameter(0, 6.6);
    fitFcn2->SetParameter(1, 0.2);
    fitFcn2->SetParameter(2, 0.493677);
    fitFcn2->SetParameter(3, 0.18);
    fitFcn2->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn2->SetLineColor(kBlue + 1);
    fitFcn2->SetLineStyle(2);
    fitFcn2->SetLineWidth(2);

    Double_t min = 0.0;
    Double_t max = 10;
    Double_t loprecision = 0.001;
    Double_t hiprecision = 0.5;
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 10;

    TH1 *hout = YieldMean(hKshort_statHM, hKshort_sysHM, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *hout2 = YieldMean(hKshort_statLM, hKshort_sysLM, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    fitFcn1->SetNpx(500);
    fitFcn2->SetNpx(500);

    TCanvas *c1 = new TCanvas("c1", "HM histograms", 720, 720);
    SetCanvasStyle(c1, 0.18, 0.06, 0.06, 0.13);
    hKstar_statHM->SetMarkerStyle(21);
    hKstar_statHM->SetMarkerSize(1.3);
    hKstar_statHM->SetMaximum(5.0);
    hKstar_statHM->SetMinimum(9e-5);
    hKstar_statHM->GetXaxis()->SetRangeUser(0, 10.0);
    hKstar_statHM->Draw("pe");
    SetHistoQA(hKshort_statHM);
    gPad->SetLogy(1);
    hKshort_statHM->SetMarkerColor(kBlue + 1);
    hKshort_statHM->SetLineColor(kBlue + 1);
    hKshort_statHM->SetMarkerSize(1.3);
    hKshort_statHM->Draw("pe same");
    fitFcn1->Draw("same");

    TLegend *legend = new TLegend(0.7, 0.7, 0.80, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetHeader("FT0M: 0-1%");
    legend->AddEntry(hKstar_statHM, "K*(892)^{0}", "p");
    legend->AddEntry(hKshort_statHM, "K^{0}_{S}", "p");
    legend->AddEntry(fitFcn1, "Levy-Tsallis", "l");
    legend->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.SetTextFont(42);
    lat.DrawLatex(0.45, 0.89, "ALICE");
    lat.DrawLatex(0.45, 0.84, "|y| < 0.5");
    lat.DrawLatex(0.45, 0.79, "pp #sqrt{#it{s}} = 13.6 TeV");
    c1->SaveAs("Plots/KstarK0sSpectra_HM.png");

    TCanvas *c2 = new TCanvas("c2", "LM histograms", 720, 720);
    SetCanvasStyle(c2, 0.18, 0.06, 0.06, 0.13);
    hKstar_statLM->SetMarkerStyle(21);
    hKstar_statLM->SetMarkerSize(1.3);
    hKstar_statLM->SetMaximum(0.7);
    hKstar_statLM->SetMinimum(2e-6);
    hKstar_statLM->GetXaxis()->SetRangeUser(0, 10.0);
    hKstar_statLM->Draw("pe");
    SetHistoQA(hKshort_statLM);
    gPad->SetLogy(1);
    hKshort_statLM->SetMarkerColor(kBlue + 1);
    hKshort_statLM->SetLineColor(kBlue + 1);
    hKshort_statLM->SetMarkerSize(1.3);
    hKshort_statLM->Draw("pe same");
    fitFcn2->Draw("same");

    TLegend *legend2 = new TLegend(0.7, 0.7, 0.80, 0.92);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextFont(42);
    legend2->SetTextSize(0.03);
    legend2->SetHeader("FT0M: 70-100%");
    legend2->AddEntry(hKstar_statLM, "K*(892)^{0}", "p");
    legend2->AddEntry(hKshort_statLM, "K^{0}_{S}", "p");
    legend2->AddEntry(fitFcn2, "Levy-Tsallis", "l");
    legend2->Draw();

    lat.DrawLatex(0.45, 0.89, "ALICE");
    lat.DrawLatex(0.45, 0.84, "|y| < 0.5");
    lat.DrawLatex(0.45, 0.79, "pp #sqrt{#it{s}} = 13.6 TeV");
    c2->SaveAs("Plots/KstarK0sSpectra_LM.png");

    TH1F *hKstarK0sRatioHM_stat = (TH1F *)hKstar_statHM->Clone("hKstarK0sRatioHM_stat");
    TH1F *hKstarK0sRatioHM_sys = (TH1F *)hKstar_sysHM->Clone("hKstarK0sRatioHM_sys");
    TH1F *hKstarK0sRatioHM_Uncorr = (TH1F *)hKstar_sysUncorrHM->Clone("hKstarK0sRatioHM_Uncorr");

    TH1F *hKstarK0sRatioLM_stat = (TH1F *)hKstar_statLM->Clone("hKstarK0sRatioLM_stat");
    TH1F *hKstarK0sRatioLM_sys = (TH1F *)hKstar_sysLM->Clone("hKstarK0sRatioLM_sys");
    TH1F *hKstarK0sRatioLM_Uncorr = (TH1F *)hKstar_sysUncorrLM->Clone("hKstarK0sRatioLM_Uncorr");

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

    int TotalBinsKstar = hKstar_statHM->GetNbinsX();

    for (int ibin = 1; ibin <= TotalBinsKstar; ++ibin)
    {
        // ---------------------------------------------------------
        // K* value and uncertainties
        // ---------------------------------------------------------
        double kstarYieldHM = hKstar_statHM->GetBinContent(ibin);
        double kstar_stat_errHM = hKstar_statHM->GetBinError(ibin);
        double kstar_sys_errHM = hKstar_sysUncorrHM->GetBinError(ibin); // uncorrelated systematic error

        double kstarYieldLM = hKstar_statLM->GetBinContent(ibin);
        double kstar_stat_errLM = hKstar_statLM->GetBinError(ibin);
        double kstar_sys_errLM = hKstar_sysUncorrLM->GetBinError(ibin); // uncorrelated systematic error

        // ---------------------------------------------------------
        // K* bin boundaries
        // ---------------------------------------------------------
        double xLow = hKstar_statHM->GetBinLowEdge(ibin);
        double xHigh = xLow + hKstar_statHM->GetBinWidth(ibin);
        double xCenter = hKstar_statHM->GetBinCenter(ibin);

        // // ---------------------------------------------------------
        // // Kaon value from fit
        // // ---------------------------------------------------------
        // double kaonYieldHM = fitFcn1->Eval(xCenter);
        // double kaonYieldLM = fitFcn2->Eval(xCenter);

        // ---------------------------------------------------------
        // Kaon value:
        // If K* bin exactly matches a K0S bin, use the measured
        // K0S bin content directly. Otherwise, use the fit value.
        // ---------------------------------------------------------
        double kaonYieldHM = getKaonYield(hKshort_statHM, fitFcn1, xLow, xHigh);
        double kaonYieldLM = getKaonYield(hKshort_statLM, fitFcn2, xLow, xHigh);

        // ---------------------------------------------------------
        // Protect against invalid values
        // ---------------------------------------------------------
        if (kstarYieldHM <= 0 || kaonYieldHM <= 0)
        {
            hKstarK0sRatioHM_stat->SetBinContent(ibin, 0);
            hKstarK0sRatioHM_stat->SetBinError(ibin, 0);

            hKstarK0sRatioHM_sys->SetBinContent(ibin, 0);
            hKstarK0sRatioHM_sys->SetBinError(ibin, 0);

            hKstarK0sRatioHM_Uncorr->SetBinContent(ibin, 0);
            hKstarK0sRatioHM_Uncorr->SetBinError(ibin, 0);

            continue;
        }

        if (kstarYieldLM <= 0 || kaonYieldLM <= 0)
        {
            hKstarK0sRatioLM_stat->SetBinContent(ibin, 0);
            hKstarK0sRatioLM_stat->SetBinError(ibin, 0);

            hKstarK0sRatioLM_sys->SetBinContent(ibin, 0);
            hKstarK0sRatioLM_sys->SetBinError(ibin, 0);

            hKstarK0sRatioLM_Uncorr->SetBinContent(ibin, 0);
            hKstarK0sRatioLM_Uncorr->SetBinError(ibin, 0);

            continue;
        }

        // ---------------------------------------------------------
        // Ratio
        // ---------------------------------------------------------
        double ratioHM = kstarYieldHM / kaonYieldHM;
        double ratioLM = kstarYieldLM / kaonYieldLM;

        // =========================================================
        // 1. STATISTICAL UNCERTAINTY
        // =========================================================

        double kaon_stat_errHM = getKaonErr(hKshort_statHM, xLow, xHigh);
        double rel_kstar_statHM = kstar_stat_errHM / kstarYieldHM;
        double rel_kaon_statHM = kaon_stat_errHM / kaonYieldHM;

        double kaon_stat_errLM = getKaonErr(hKshort_statLM, xLow, xHigh);
        double rel_kstar_statLM = kstar_stat_errLM / kstarYieldLM;
        double rel_kaon_statLM = kaon_stat_errLM / kaonYieldLM;

        double ratio_stat_errHM = ratioHM * std::sqrt(rel_kstar_statHM * rel_kstar_statHM + rel_kaon_statHM * rel_kaon_statHM);
        double ratio_stat_errLM = ratioLM * std::sqrt(rel_kstar_statLM * rel_kstar_statLM + rel_kaon_statLM * rel_kaon_statLM);

        // =========================================================
        // 2. SYSTEMATIC UNCERTAINTY
        //
        // K*   -> uncorrelated systematic
        // KaonYieldHM -> total systematic
        // =========================================================

        double kaon_sys_errHM = getKaonErr(hKshort_sysHM, xLow, xHigh);
        double rel_kstar_sysHM = kstar_sys_errHM / kstarYieldHM;
        double rel_kaon_sysHM = kaon_sys_errHM / kaonYieldHM;

        double kaon_sys_errLM = getKaonErr(hKshort_sysLM, xLow, xHigh);
        double rel_kstar_sysLM = kstar_sys_errLM / kstarYieldLM;
        double rel_kaon_sysLM = kaon_sys_errLM / kaonYieldLM;

        double ratio_sys_errHM = ratioHM * std::sqrt(rel_kstar_sysHM * rel_kstar_sysHM + rel_kaon_sysHM * rel_kaon_sysHM);
        double ratio_sys_errLM = ratioLM * std::sqrt(rel_kstar_sysLM * rel_kstar_sysLM + rel_kaon_sysLM * rel_kaon_sysLM);

        // =========================================================
        // 3. SYSTEMATIC UNCERTAINTY
        //
        // K*   -> uncorrelated systematic
        // KaonYieldHM -> uncorrelated systematic
        // =========================================================
        double kaon_sysUncorr_errHM = getKaonErr(hKshort_sysUncorrHM, xLow, xHigh);
        double rel_kaon_sysUncorrHM = kaon_sysUncorr_errHM / kaonYieldHM;

        double kaon_sysUncorr_errLM = getKaonErr(hKshort_sysUncorrLM, xLow, xHigh);
        double rel_kaon_sysUncorrLM = kaon_sysUncorr_errLM / kaonYieldLM;

        double ratio_sysUncorr_errHM = ratioHM * std::sqrt(rel_kstar_sysHM * rel_kstar_sysHM + rel_kaon_sysUncorrHM * rel_kaon_sysUncorrHM);
        double ratio_sysUncorr_errLM = ratioLM * std::sqrt(rel_kstar_sysLM * rel_kstar_sysLM + rel_kaon_sysUncorrLM * rel_kaon_sysUncorrLM);

        // =========================================================
        // Store
        // =========================================================

        hKstarK0sRatioHM_stat->SetBinContent(ibin, ratioHM);
        hKstarK0sRatioHM_stat->SetBinError(ibin, ratio_stat_errHM);

        hKstarK0sRatioHM_sys->SetBinContent(ibin, ratioHM);
        hKstarK0sRatioHM_sys->SetBinError(ibin, ratio_sys_errHM);

        hKstarK0sRatioHM_Uncorr->SetBinContent(ibin, ratioHM);
        hKstarK0sRatioHM_Uncorr->SetBinError(ibin, ratio_sysUncorr_errHM);

        hKstarK0sRatioLM_stat->SetBinContent(ibin, ratioLM);
        hKstarK0sRatioLM_stat->SetBinError(ibin, ratio_stat_errLM);

        hKstarK0sRatioLM_sys->SetBinContent(ibin, ratioLM);
        hKstarK0sRatioLM_sys->SetBinError(ibin, ratio_sys_errLM);

        hKstarK0sRatioLM_Uncorr->SetBinContent(ibin, ratioLM);
        hKstarK0sRatioLM_Uncorr->SetBinError(ibin, ratio_sysUncorr_errLM);
    }

    TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 720, 720);
    SetCanvasStyle(cRatio, 0.16, 0.06, 0.06, 0.13);
    double pad1Size, pad2Size;
    canvas_style(cRatio, pad1Size, pad2Size);
    cRatio->cd(1);

    SetHistoQA(hKstarK0sRatioHM_stat);
    hKstarK0sRatioHM_stat->GetYaxis()->SetTitle("K^{*0} / K^{0}_{S}");
    hKstarK0sRatioHM_stat->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hKstarK0sRatioHM_stat->SetMarkerStyle(20);
    hKstarK0sRatioHM_stat->SetMarkerSize(1.45);
    hKstarK0sRatioHM_stat->GetYaxis()->SetTitleOffset(1.60 * pad1Size);
    hKstarK0sRatioHM_stat->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetYaxis()->SetRangeUser(0.08, 1.25);
    hKstarK0sRatioHM_stat->GetXaxis()->SetRangeUser(0.0, 10.0);
    hKstarK0sRatioHM_stat->SetMarkerColor(kRed + 1);
    hKstarK0sRatioHM_stat->SetLineColor(kRed + 1);
    hKstarK0sRatioHM_stat->Draw("pe");
    SetHistoQA(hKstarK0sRatioHM_sys);
    hKstarK0sRatioHM_sys->SetFillStyle(0);
    hKstarK0sRatioHM_sys->SetLineColor(kRed + 1);
    hKstarK0sRatioHM_sys->SetMarkerColor(kRed + 1);
    hKstarK0sRatioHM_sys->Draw("e2 same");

    SetHistoQA(hKstarK0sRatioLM_stat);
    SetHistoQA(hKstarK0sRatioLM_sys);
    hKstarK0sRatioLM_stat->SetMarkerStyle(21);
    hKstarK0sRatioLM_stat->SetMarkerSize(1.45);
    hKstarK0sRatioLM_stat->SetMarkerColor(kBlue + 1);
    hKstarK0sRatioLM_stat->SetLineColor(kBlue + 1);
    hKstarK0sRatioLM_stat->Draw("pe same");
    hKstarK0sRatioLM_sys->SetFillStyle(0);
    hKstarK0sRatioLM_sys->SetLineColor(kBlue + 1);
    hKstarK0sRatioLM_sys->SetMarkerColor(kBlue + 1);
    hKstarK0sRatioLM_sys->Draw("e2 same");

    TLegend *legendRatio = new TLegend(0.64, 0.78, 0.80, 0.92);
    legendRatio->SetBorderSize(0);
    legendRatio->SetFillStyle(0);
    legendRatio->SetTextFont(42);
    legendRatio->SetTextSize(0.05);
    legendRatio->AddEntry(hKstarK0sRatioHM_stat, "FT0M: 0-1%", "p");
    legendRatio->AddEntry(hKstarK0sRatioLM_stat, "FT0M: 70-100%", "p");
    legendRatio->Draw();

    lat.SetTextSize(0.05);
    lat.DrawLatex(0.22, 0.89, "ALICE");
    lat.DrawLatex(0.22, 0.83, "|y| < 0.5");
    lat.DrawLatex(0.22, 0.76, "pp #sqrt{#it{s}} = 13.6 TeV");

    cRatio->cd(2);
    TH1F *hDoubleRatioStat = (TH1F *)hKstarK0sRatioHM_stat->Clone("hDoubleRatioStat");
    hDoubleRatioStat->Divide(hKstarK0sRatioLM_stat);
    TH1F *hDoubleRatioSys = (TH1F *)hKstarK0sRatioHM_Uncorr->Clone("hDoubleRatioSys");
    hDoubleRatioSys->Divide(hKstarK0sRatioLM_Uncorr);
    SetHistoQA(hDoubleRatioStat);
    hDoubleRatioStat->GetYaxis()->SetTitle("#frac{High mult.}{Low mult.}");
    hDoubleRatioStat->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hDoubleRatioStat->SetMarkerStyle(20);
    hDoubleRatioStat->SetMarkerSize(1.45);
    hDoubleRatioStat->GetYaxis()->SetTitleOffset(1.90 * pad2Size);
    hDoubleRatioStat->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    hDoubleRatioStat->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hDoubleRatioStat->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hDoubleRatioStat->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hDoubleRatioStat->GetYaxis()->SetRangeUser(0.45, 1.75);
    hDoubleRatioStat->GetXaxis()->SetRangeUser(0.0, 10.0);
    hDoubleRatioStat->SetMarkerColor(kRed + 1);
    hDoubleRatioStat->SetLineColor(kRed + 1);
    hDoubleRatioStat->GetYaxis()->SetNdivisions(505);
    hDoubleRatioStat->Draw("pe");
    SetHistoQA(hDoubleRatioSys);
    hDoubleRatioSys->SetFillStyle(0);
    hDoubleRatioSys->SetLineColor(kRed + 1);
    hDoubleRatioSys->SetMarkerColor(kRed + 1);
    hDoubleRatioSys->Draw("e2 same");

    TLine *line = new TLine(0.0, 1.0, 10.0, 1.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    cRatio->SaveAs("Plots/KstarK0sRatio.png");

    // ---------------------------------------------------------
    // Calculate significance with respect to unity
    // significance = (value - 1) / systematic error
    // ---------------------------------------------------------
    TH1F *hDoubleRatioSignificance = (TH1F *)hDoubleRatioSys->Clone("hDoubleRatioSignificance");

    hDoubleRatioSignificance->Reset();

    int nBins = hDoubleRatioSys->GetNbinsX();

    for (int ibin = 1; ibin <= nBins; ++ibin)
    {
        double value = hDoubleRatioSys->GetBinContent(ibin);
        double sysErr = hDoubleRatioSys->GetBinError(ibin);
        double statErr = hDoubleRatioStat->GetBinError(ibin);

        double significance = 0.0;

        if (sysErr > 0)
        {
            significance = abs(value - 1.0) / sysErr;
        }

        hDoubleRatioSignificance->SetBinContent(ibin, significance);
        hDoubleRatioSignificance->SetBinError(ibin, 0.0);
    }

    for (int ibin = 1; ibin <= nBins; ++ibin)
    {
        double pt = hDoubleRatioSignificance->GetBinCenter(ibin);
        double significance = hDoubleRatioSignificance->GetBinContent(ibin);

        cout << "pT = " << pt
             << " GeV/c, significance = "
             << significance << endl;
    }
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