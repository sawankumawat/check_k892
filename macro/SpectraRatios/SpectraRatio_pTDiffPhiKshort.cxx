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

void SpectraRatio_pTDiffPhiKshort()
{
    bool plotEPOS = true;
    gStyle->SetOptStat(0);
    TFile *fPhi = OpenFile("PhiResults/PhiAnalysisResults_pprun3.root");
    TFile *fKshort = OpenFile("K0s_Run3_Results/Spectra_K0Short.root"); // Run 3
    // TFile *fKshort = OpenFile("K0s_Run2_Results/K0s_Run2_Spectra.root"); // Run 2

    TH1F *hKstar_statHM = (TH1F *)fPhi->Get("Spectra/01/hStat");
    TH1F *hKstar_sysHM = (TH1F *)fPhi->Get("Spectra/01/hSys");
    TH1F *hKstar_sysUncorrHM = (TH1F *)fPhi->Get("Spectra/01/hSys");

    TH1F *hKstar_statLM = (TH1F *)fPhi->Get("Spectra/10/hStat");
    TH1F *hKstar_sysLM = (TH1F *)fPhi->Get("Spectra/10/hSys");
    TH1F *hKstar_sysUncorrLM = (TH1F *)fPhi->Get("Spectra/10/hSys");

    // Run 3 K0s results
    TH1F *hKshort_statHM = (TH1F *)fKshort->Get("hK0s_0_1");
    TH1F *hKshort_sysHM = (TH1F *)fKshort->Get("hK0sSyst_0_1");
    TH1F *hKshort_sysUncorrHM = (TH1F *)fKshort->Get("hK0sSyst_0_1"); // Does not have uncorrelated systematic, so using total systematic

    TH1F *hKshort_statLM = (TH1F *)fKshort->Get("hK0s_70_100");
    TH1F *hKshort_sysLM = (TH1F *)fKshort->Get("hK0sSyst_70_100");
    TH1F *hKshort_sysUncorrLM = (TH1F *)fKshort->Get("hK0sSyst_70_100"); // Does not have uncorrelated systematic, so using total systematic

    // ////Run 2 K0s results
    // TH1F *hKshortHM = (TH1F *)fKshort->Get("Table 1/Hist1D_y1");
    // TH1F *hKshort_statHM = (TH1F *)fKshort->Get("Table 1/Hist1D_y1_e1");
    // TH1F *hKshort_sysHM = (TH1F *)fKshort->Get("Table 1/Hist1D_y1_e2");
    // TH1F *hKshort_sysUncorrHM = (TH1F *)fKshort->Get("Table 1/Hist1D_y1_e3");

    // TH1F *hKshortLM = (TH1F *)fKshort->Get("Table 1/Hist1D_y10");
    // TH1F *hKshort_statLM = (TH1F *)fKshort->Get("Table 1/Hist1D_y10_e1");
    // TH1F *hKshort_sysLM = (TH1F *)fKshort->Get("Table 1/Hist1D_y10_e2");
    // TH1F *hKshort_sysUncorrLM = (TH1F *)fKshort->Get("Table 1/Hist1D_y10_e3");

    // for (int ibin = 1; ibin <= hKshortLM->GetNbinsX() - 1; ++ibin)
    // {
    //     double binContent = hKshortLM->GetBinContent(ibin);
    //     double binErrorStat = hKshort_statLM->GetBinContent(ibin);
    //     double binErrorSys = hKshort_sysLM->GetBinContent(ibin);
    //     double binErrorSysUncorr = hKshort_sysUncorrLM->GetBinContent(ibin);

    //     hKshort_statLM->SetBinContent(ibin, binContent);
    //     hKshort_statLM->SetBinError(ibin, binErrorStat);

    //     hKshort_sysLM->SetBinContent(ibin, binContent);
    //     hKshort_sysLM->SetBinError(ibin, binErrorSys);

    //     hKshort_sysUncorrLM->SetBinContent(ibin, binContent);
    //     hKshort_sysUncorrLM->SetBinError(ibin, binErrorSysUncorr);

    //     double binContentHM = hKshortHM->GetBinContent(ibin);
    //     double binErrorStatHM = hKshort_statHM->GetBinContent(ibin);
    //     double binErrorSysHM = hKshort_sysHM->GetBinContent(ibin);
    //     double binErrorSysUncorrHM = hKshort_sysUncorrHM->GetBinContent(ibin);

    //     hKshort_statHM->SetBinContent(ibin, binContentHM);
    //     hKshort_statHM->SetBinError(ibin, binErrorStatHM);

    //     hKshort_sysHM->SetBinContent(ibin, binContentHM);
    //     hKshort_sysHM->SetBinError(ibin, binErrorSysHM);

    //     hKshort_sysUncorrHM->SetBinContent(ibin, binContentHM);
    //     hKshort_sysUncorrHM->SetBinError(ibin, binErrorSysUncorrHM);

    //     cout << "K0s bin " << ibin << ", yield LM: " << binContent << ", yield HM: " << binContentHM << endl;
    // }

    cout << "Number of kaon bins in HM: " << hKshort_statHM->GetNbinsX() << endl;
    cout << "Number of kaon bins in LM: " << hKshort_statLM->GetNbinsX() << endl;

    int TotalBinsKaon = hKshort_statHM->GetNbinsX();

    TF1 *fitFcn1 = new TF1("fitfunc", FuncLavy, 0.0, 10.0, 4);
    fitFcn1->SetParameter(0, 4.6);
    fitFcn1->SetParameter(1, 0.2);
    fitFcn1->SetParameter(2, 0.493677);
    fitFcn1->SetParameter(3, 0.18);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn1->SetLineColor(kRed + 1);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);

    TF1 *fitFcn2 = new TF1("fitfunc2", FuncLavy, 0.0, 10.0, 4);
    fitFcn2->SetParameter(0, 4.6);
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
    Double_t minfit = 0.1;
    Double_t maxfit = 10;

    TH1 *hout = YieldMean(hKshort_statHM, hKshort_sysHM, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *hout2 = YieldMean(hKshort_statLM, hKshort_sysLM, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    cout << "Fit parameters for HM: " << endl;
    for (int i = 0; i < fitFcn1->GetNpar(); ++i)
    {
        cout << "Parameter " << i << ": " << fitFcn1->GetParameter(i) << endl;
    }

    cout << "Fit parameters for LM: " << endl;
    for (int i = 0; i < fitFcn2->GetNpar(); ++i)
    {
        cout << "Parameter " << i << ": " << fitFcn2->GetParameter(i) << endl;
    }

    fitFcn1->SetNpx(500);
    fitFcn2->SetNpx(500);

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.SetTextFont(42);

    TCanvas *c1 = new TCanvas("c1", "HM histograms", 720, 720);
    SetCanvasStyle(c1, 0.18, 0.06, 0.06, 0.13);
    hKstar_statHM->SetMarkerStyle(21);
    hKstar_statHM->SetMarkerSize(1.3);
    hKstar_statHM->SetMaximum(5.0);
    hKstar_statHM->SetMinimum(9e-5);
    hKstar_statHM->GetXaxis()->SetRangeUser(0.5, 6.0);
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
    legend->AddEntry(hKstar_statHM, "#Phi(1020)", "p");
    legend->AddEntry(hKshort_statHM, "K^{0}_{S}", "p");
    legend->AddEntry(fitFcn1, "Levy-Tsallis", "l");
    legend->Draw();

    lat.DrawLatex(0.45, 0.89, "ALICE");
    lat.DrawLatex(0.45, 0.84, "|y| < 0.5");
    lat.DrawLatex(0.45, 0.79, "pp #sqrt{#it{s}} = 13.6 TeV");
    // c1->SaveAs("Plots/KstarK0sSpectra_HM.png");

    TCanvas *c2 = new TCanvas("c2", "LM histograms", 720, 720);
    SetCanvasStyle(c2, 0.18, 0.06, 0.06, 0.13);
    hKstar_statLM->SetMarkerStyle(21);
    hKstar_statLM->SetMarkerSize(1.3);
    hKstar_statLM->SetMaximum(0.7);
    hKstar_statLM->SetMinimum(2e-6);
    hKstar_statLM->GetXaxis()->SetRangeUser(0.5, 6.0);
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
    legend2->AddEntry(hKstar_statLM, "#Phi(1020)", "p");
    legend2->AddEntry(hKshort_statLM, "K^{0}_{S}", "p");
    legend2->AddEntry(fitFcn2, "Levy-Tsallis", "l");
    legend2->Draw();

    lat.DrawLatex(0.45, 0.89, "ALICE");
    lat.DrawLatex(0.45, 0.84, "|y| < 0.5");
    lat.DrawLatex(0.45, 0.79, "pp #sqrt{#it{s}} = 13.6 TeV");
    // c2->SaveAs("Plots/KstarK0sSpectra_LM.png");

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
        double kstar_sys_UncorrerrHM = hKstar_sysUncorrHM->GetBinError(ibin); // uncorrelated systematic error
        double kstar_sys_errHM = hKstar_sysHM->GetBinError(ibin);             // total systematic error

        double kstarYieldLM = hKstar_statLM->GetBinContent(ibin);
        double kstar_stat_errLM = hKstar_statLM->GetBinError(ibin);
        double kstar_sys_UncorrerrLM = hKstar_sysUncorrLM->GetBinError(ibin); // uncorrelated systematic error
        double kstar_sys_errLM = hKstar_sysLM->GetBinError(ibin);             // total systematic error

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
        cout << "Bin " << ibin << ", bin value " << xCenter << ", Kstar yield: " << kstarYieldHM << ", Kshort yield : " << kaonYieldHM << endl;

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
        // K*   -> total systematic
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

        double kstar_sysUncorr_errHM = hKstar_sysUncorrHM->GetBinError(ibin); // uncorrelated systematic error
        double rel_kstar_sysUncorrHM = kstar_sysUncorr_errHM / kstarYieldHM;

        double kstar_sysUncorr_errLM = hKstar_sysUncorrLM->GetBinError(ibin); // uncorrelated systematic error
        double rel_kstar_sysUncorrLM = kstar_sysUncorr_errLM / kstarYieldLM;

        double kaon_sysUncorr_errHM = getKaonErr(hKshort_sysUncorrHM, xLow, xHigh);
        double rel_kaon_sysUncorrHM = kaon_sysUncorr_errHM / kaonYieldHM;

        double kaon_sysUncorr_errLM = getKaonErr(hKshort_sysUncorrLM, xLow, xHigh);
        double rel_kaon_sysUncorrLM = kaon_sysUncorr_errLM / kaonYieldLM;

        double ratio_sysUncorr_errHM = ratioHM * std::sqrt(rel_kstar_sysUncorrHM * rel_kstar_sysUncorrHM + rel_kaon_sysUncorrHM * rel_kaon_sysUncorrHM);
        double ratio_sysUncorr_errLM = ratioLM * std::sqrt(rel_kstar_sysUncorrLM * rel_kstar_sysUncorrLM + rel_kaon_sysUncorrLM * rel_kaon_sysUncorrLM);

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
    hKstarK0sRatioHM_stat->GetYaxis()->SetTitle("#Phi(1020) / K^{0}_{S}");
    hKstarK0sRatioHM_stat->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hKstarK0sRatioHM_stat->SetMarkerStyle(20);
    hKstarK0sRatioHM_stat->SetMarkerSize(1.45);
    hKstarK0sRatioHM_stat->GetYaxis()->SetTitleOffset(1.90 * pad1Size);
    hKstarK0sRatioHM_stat->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetYaxis()->SetRangeUser(0.02, 0.49);
    hKstarK0sRatioHM_stat->GetXaxis()->SetRangeUser(0.5, 6.0);
    hKstarK0sRatioHM_stat->SetMarkerColor(kRed + 1);
    hKstarK0sRatioHM_stat->SetLineColor(kRed + 1);
    hKstarK0sRatioHM_stat->Draw("pe");
    SetHistoQA(hKstarK0sRatioHM_sys);
    hKstarK0sRatioHM_sys->SetFillStyle(0);
    hKstarK0sRatioHM_sys->SetLineColor(kRed + 1);
    hKstarK0sRatioHM_sys->SetMarkerColor(kRed + 1);
    hKstarK0sRatioHM_sys->Draw("e2 same");
    hKstarK0sRatioHM_Uncorr->SetFillStyle(1001);
    hKstarK0sRatioHM_Uncorr->SetMarkerStyle(20);
    hKstarK0sRatioHM_Uncorr->SetMarkerColor(kRed + 1);
    hKstarK0sRatioHM_Uncorr->SetFillColorAlpha(kRed + 1, 0.35);
    // hKstarK0sRatioHM_Uncorr->Draw("e2 same");

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
    hKstarK0sRatioLM_Uncorr->SetFillStyle(1001);
    hKstarK0sRatioLM_Uncorr->SetMarkerStyle(21);
    hKstarK0sRatioLM_Uncorr->SetMarkerColor(kBlue + 1);
    hKstarK0sRatioLM_Uncorr->SetFillColorAlpha(kBlue + 1, 0.35);
    // hKstarK0sRatioLM_Uncorr->Draw("e2 same");

    // Now open the spectra file in EPOS
    TFile *fEPOS_spectra = OpenFile("/home/sawan/Storage/EPOS4_Resonance/yield_outfile_light_oct2026_fixed.root");

    // TH1F *hEPOS_UrQMDOFF_PhiHM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent9");
    // TH1F *hEPOS_UrQMDOFF_PhiLM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent0");
    // TH1F *hEPOS_UrQMDON_PhiHM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent9");
    // TH1F *hEPOS_UrQMDON_PhiLM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent0");

    TH1F *hEPOS_UrQMDOFF_PhiHM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent9");
    TH1F *hEPOS_UrQMDOFF_PhiLM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent0");
    TH1F *hEPOS_UrQMDON_PhiHM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent9");
    TH1F *hEPOS_UrQMDON_PhiLM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent0");

    TH1F *hEPOS_UrQMDOFF_K0sHM = (TH1F *)fEPOS_spectra->Get("K0s_Mult_OFF/hK0s_Pt_OFF_Cent9");
    TH1F *hEPOS_UrQMDOFF_K0sLM = (TH1F *)fEPOS_spectra->Get("K0s_Mult_OFF/hK0s_Pt_OFF_Cent0");
    TH1F *hEPOS_UrQMDON_K0sHM = (TH1F *)fEPOS_spectra->Get("K0s_Mult_OFF/hK0s_Pt_OFF_Cent9");
    TH1F *hEPOS_UrQMDON_K0sLM = (TH1F *)fEPOS_spectra->Get("K0s_Mult_OFF/hK0s_Pt_OFF_Cent0");

    //     TH1F *hEPOS_UrQMDOFF_K0sHM = (TH1F *)fEPOS_spectra->Get("K0s_Mult_status0/hK0s_Pt_status0_Cent9");
    // TH1F *hEPOS_UrQMDOFF_K0sLM = (TH1F *)fEPOS_spectra->Get("K0s_Mult_status0/hK0s_Pt_status0_Cent0");
    // TH1F *hEPOS_UrQMDON_K0sHM = (TH1F *)fEPOS_spectra->Get("K0s_Mult_status0/hK0s_Pt_status0_Cent9");
    // TH1F *hEPOS_UrQMDON_K0sLM = (TH1F *)fEPOS_spectra->Get("K0s_Mult_status0/hK0s_Pt_status0_Cent0");

    double phiBR = 0.489; // Branching ratio for Phi(1020) -> K+ K-
    hEPOS_UrQMDON_PhiHM->Scale(1 / phiBR);
    hEPOS_UrQMDON_PhiLM->Scale(1 / phiBR);

    if (!hEPOS_UrQMDOFF_PhiHM || !hEPOS_UrQMDOFF_PhiLM || !hEPOS_UrQMDON_PhiHM || !hEPOS_UrQMDON_PhiLM)
    {
        std::cerr << "Error: One or more Phi histograms not found in the EPOS spectra file." << std::endl;
        return;
    }
    if (!hEPOS_UrQMDOFF_K0sHM || !hEPOS_UrQMDOFF_K0sLM || !hEPOS_UrQMDON_K0sHM || !hEPOS_UrQMDON_K0sLM)
    {
        std::cerr << "Error: One or more K0s histograms not found in the EPOS spectra file." << std::endl;
        return;
    }

    // TH1F *hEPOS_UrQMDOFF_KstarK0sRatioHM = (TH1F *)hEPOS_UrQMDOFF_PhiHM->Clone("hEPOS_UrQMDOFF_KstarK0sRatioHM");
    // hEPOS_UrQMDOFF_KstarK0sRatioHM->Divide(hEPOS_UrQMDOFF_K0sHM);
    // hEPOS_UrQMDOFF_KstarK0sRatioHM->Scale(0.5); // Since kstar was average but not K0s
    // TH1F *hEPOS_UrQMDOFF_KstarK0sRatioLM = (TH1F *)hEPOS_UrQMDOFF_PhiLM->Clone("hEPOS_UrQMDOFF_KstarK0sRatioLM");
    // hEPOS_UrQMDOFF_KstarK0sRatioLM->Divide(hEPOS_UrQMDOFF_K0sLM);
    // hEPOS_UrQMDOFF_KstarK0sRatioLM->Scale(0.5); // Since kstar was average but not K0s

    // ============================================================
    // EPOS K* / K0s ratio using pT-integrated bins
    // ============================================================

    // double mergedBins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0}; // Note: since the graph is used for draw, it will draw upto bin center. so if last bin is 8-10, then it will draw upto 9 only.
    double mergedBins[] = {0.3, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 4.0, 5.0, 6.0}; // Note: since the graph is used for draw, it will draw upto bin center. so if last bin is 8-10, then it will draw upto 9 only.

    const int nMergedBins = sizeof(mergedBins) / sizeof(double) - 1;

    // ------------------------------------------------------------
    // Function to calculate integrated yield and its statistical
    // uncertainty in a pT interval
    // ------------------------------------------------------------
    auto getIntegral = [](TH1F *h, double xLow, double xHigh, double &yield, double &error)
    {
        int binLow = h->GetXaxis()->FindBin(xLow + 1e-6);
        int binHigh = h->GetXaxis()->FindBin(xHigh - 1e-6);

        yield = h->IntegralAndError(binLow, binHigh, error, "width");
    };

    // ------------------------------------------------------------
    // Create empty ratio histograms
    // ------------------------------------------------------------
    TH1F *hEPOS_UrQMDOFF_KstarK0sRatioHM = new TH1F("hEPOS_UrQMDOFF_KstarK0sRatioHM", ";#it{p}_{T} (GeV/#it{c});K^{*0}/K^{0}_{S}", nMergedBins, mergedBins);

    TH1F *hEPOS_UrQMDOFF_KstarK0sRatioLM = new TH1F("hEPOS_UrQMDOFF_KstarK0sRatioLM", ";#it{p}_{T} (GeV/#it{c});K^{*0}/K^{0}_{S}", nMergedBins, mergedBins);

    // ------------------------------------------------------------
    // Calculate integrated ratios
    // ------------------------------------------------------------
    for (int i = 0; i < nMergedBins; ++i)
    {
        double xLow = mergedBins[i];
        double xHigh = mergedBins[i + 1];

        double kstarHM, kstarErrHM;
        double k0sHM, k0sErrHM;

        double kstarLM, kstarErrLM;
        double k0sLM, k0sErrLM;

        // HM
        getIntegral(hEPOS_UrQMDOFF_PhiHM, xLow, xHigh, kstarHM, kstarErrHM);
        getIntegral(hEPOS_UrQMDOFF_K0sHM, xLow, xHigh, k0sHM, k0sErrHM);
        // LM
        getIntegral(hEPOS_UrQMDOFF_PhiLM, xLow, xHigh, kstarLM, kstarErrLM);
        getIntegral(hEPOS_UrQMDOFF_K0sLM, xLow, xHigh, k0sLM, k0sErrLM);
        // --------------------------------------------------------
        // HM ratio
        // --------------------------------------------------------
        if (kstarHM > 0 && k0sHM > 0)
        {
            double ratioHM = kstarHM / k0sHM;
            double relKstarHM = kstarErrHM / kstarHM;
            double relK0sHM = k0sErrHM / k0sHM;

            double ratioErrHM = ratioHM * sqrt(relKstarHM * relKstarHM + relK0sHM * relK0sHM) / 2.0;

            hEPOS_UrQMDOFF_KstarK0sRatioHM->SetBinContent(i + 1, ratioHM);
            hEPOS_UrQMDOFF_KstarK0sRatioHM->SetBinError(i + 1, ratioErrHM);
        }

        // --------------------------------------------------------
        // LM ratio
        // --------------------------------------------------------
        if (kstarLM > 0 && k0sLM > 0)
        {
            double ratioLM = kstarLM / k0sLM;

            double relKstarLM = kstarErrLM / kstarLM;
            double relK0sLM = k0sErrLM / k0sLM;

            double ratioErrLM = ratioLM * sqrt(relKstarLM * relKstarLM + relK0sLM * relK0sLM) / 2.0;
            hEPOS_UrQMDOFF_KstarK0sRatioLM->SetBinContent(i + 1, ratioLM);
            hEPOS_UrQMDOFF_KstarK0sRatioLM->SetBinError(i + 1, ratioErrLM);
        }
    }

    hEPOS_UrQMDOFF_KstarK0sRatioHM->Smooth(2);
    hEPOS_UrQMDOFF_KstarK0sRatioLM->Smooth(2);

    // Build a line + band that extend to the outer bin edges
    auto makeCurve = [](TH1 *h, Color_t col, Style_t lstyle, TGraph *&gLine, TGraphErrors *&gBand)
    {
        gLine = new TGraph();
        gBand = new TGraphErrors();
        const int n = h->GetNbinsX();

        // first/last filled bins
        int first = -1, last = -1;
        for (int i = 1; i <= n; ++i)
            if (h->GetBinContent(i) > 0)
            {
                if (first < 0)
                    first = i;
                last = i;
            }
        if (first < 0)
            return;

        auto add = [&](double x, int ibin)
        {
            gLine->SetPoint(gLine->GetN(), x, h->GetBinContent(ibin));
            int k = gBand->GetN();
            gBand->SetPoint(k, x, h->GetBinContent(ibin));
            gBand->SetPointError(k, 0, h->GetBinError(ibin));
        };

        add(h->GetXaxis()->GetXmin(), first); // left edge
        for (int i = first; i <= last; ++i)
            if (h->GetBinContent(i) > 0)
                add(h->GetBinCenter(i), i);  // bin centers
        add(h->GetXaxis()->GetXmax(), last); // right edge (10 GeV/c)

        gLine->SetLineColor(col);
        gLine->SetLineWidth(2);
        gLine->SetLineStyle(lstyle);
        gBand->SetFillColorAlpha(col, 0.3);
        gBand->SetLineColor(col);
    };

    TGraph *gHM, *gLM;
    TGraphErrors *bHM, *bLM;
    makeCurve(hEPOS_UrQMDOFF_KstarK0sRatioHM, kRed + 1, 1, gHM, bHM);
    makeCurve(hEPOS_UrQMDOFF_KstarK0sRatioLM, kBlue + 1, 1, gLM, bLM);
    if (plotEPOS)
    {
        bHM->Draw("3 same");
        bLM->Draw("3 same"); // bands
        gHM->Draw("L same");
        gLM->Draw("L same"); // lines
    }

    // ============================================================
    // UrQMD ON K* / K0s ratio using pT-integrated bins
    // ============================================================

    // ------------------------------------------------------------
    // Create UrQMD ON ratio histograms
    // ------------------------------------------------------------
    TH1F *hEPOS_UrQMDON_KstarK0sRatioHM = new TH1F("hEPOS_UrQMDON_KstarK0sRatioHM", ";#it{p}_{T} (GeV/#it{c});K^{*0}/K^{0}_{S}", nMergedBins, mergedBins);
    TH1F *hEPOS_UrQMDON_KstarK0sRatioLM = new TH1F("hEPOS_UrQMDON_KstarK0sRatioLM", ";#it{p}_{T} (GeV/#it{c});K^{*0}/K^{0}_{S}", nMergedBins, mergedBins);

    // ------------------------------------------------------------
    // Calculate integrated ratios
    // ------------------------------------------------------------
    for (int i = 0; i < nMergedBins; ++i)
    {
        double xLow = mergedBins[i];
        double xHigh = mergedBins[i + 1];

        double kstarHM, kstarErrHM;
        double k0sHM, k0sErrHM;

        double kstarLM, kstarErrLM;
        double k0sLM, k0sErrLM;

        // --------------------------------------------------------
        // HM
        // --------------------------------------------------------
        getIntegral(hEPOS_UrQMDON_PhiHM, xLow, xHigh, kstarHM, kstarErrHM);
        getIntegral(hEPOS_UrQMDON_K0sHM, xLow, xHigh, k0sHM, k0sErrHM);

        // --------------------------------------------------------
        // LM
        // --------------------------------------------------------
        getIntegral(hEPOS_UrQMDON_PhiLM, xLow, xHigh, kstarLM, kstarErrLM);
        getIntegral(hEPOS_UrQMDON_K0sLM, xLow, xHigh, k0sLM, k0sErrLM);

        // ========================================================
        // HM ratio
        // ========================================================
        if (kstarHM > 0 && k0sHM > 0)
        {
            double ratioHM = kstarHM / k0sHM;
            double relKstarHM = kstarErrHM / kstarHM;

            double relK0sHM = k0sErrHM / k0sHM;

            double ratioErrHM = ratioHM * sqrt(relKstarHM * relKstarHM + relK0sHM * relK0sHM) / 2.0;

            hEPOS_UrQMDON_KstarK0sRatioHM->SetBinContent(i + 1, ratioHM);
            hEPOS_UrQMDON_KstarK0sRatioHM->SetBinError(i + 1, ratioErrHM); // Scale down the error for visualization
        }

        // ========================================================
        // LM ratio
        // ========================================================
        if (kstarLM > 0 && k0sLM > 0)
        {
            double ratioLM = kstarLM / k0sLM;
            double relKstarLM = kstarErrLM / kstarLM;

            double relK0sLM = k0sErrLM / k0sLM;

            double ratioErrLM = ratioLM * sqrt(relKstarLM * relKstarLM + relK0sLM * relK0sLM) / 2.0;

            hEPOS_UrQMDON_KstarK0sRatioLM->SetBinContent(i + 1, ratioLM);
            hEPOS_UrQMDON_KstarK0sRatioLM->SetBinError(i + 1, ratioErrLM);
        }
    }

    hEPOS_UrQMDON_KstarK0sRatioHM->Smooth(2);
    hEPOS_UrQMDON_KstarK0sRatioLM->Smooth(2);

    TGraph *gHM_UrQMDON, *gLM_UrQMDON;
    TGraphErrors *bHM_UrQMDON, *bLM_UrQMDON;
    makeCurve(hEPOS_UrQMDON_KstarK0sRatioHM, kRed + 1, 7, gHM_UrQMDON, bHM_UrQMDON);
    makeCurve(hEPOS_UrQMDON_KstarK0sRatioLM, kBlue + 1, 7, gLM_UrQMDON, bLM_UrQMDON);
    if (plotEPOS)
    {
        bHM_UrQMDON->Draw("3 same");
        bLM_UrQMDON->Draw("3 same"); // bands
        gHM_UrQMDON->Draw("L same");
        gLM_UrQMDON->Draw("L same"); // lines
    }

    TLegend *legendRatio = new TLegend(0.5, 0.04, 0.80, 0.36);
    legendRatio->SetBorderSize(0);
    legendRatio->SetFillStyle(0);
    legendRatio->SetTextFont(42);
    legendRatio->SetTextSize(0.04);
    legendRatio->AddEntry(hKstarK0sRatioHM_stat, "FT0M: 0-1%", "p");
    legendRatio->AddEntry(hKstarK0sRatioLM_stat, "FT0M: 70-100%", "p");
    if (plotEPOS)
    {
        legendRatio->AddEntry(gHM, "EPOS UrQMD OFF: 0-1%", "l");
        legendRatio->AddEntry(gLM, "EPOS UrQMD OFF: 70-100%", "l");
        legendRatio->AddEntry(gHM_UrQMDON, "EPOS UrQMD ON: 0-1%", "l");
        legendRatio->AddEntry(gLM_UrQMDON, "EPOS UrQMD ON: 70-100%", "l");
    }
    legendRatio->Draw();

    lat.SetTextSize(0.05);
    lat.DrawLatex(0.22, 0.89, "ALICE");
    lat.DrawLatex(0.22, 0.83, "|y| < 0.5");
    lat.DrawLatex(0.22, 0.76, "pp #sqrt{#it{s}} = 13.6 TeV");

    cRatio->cd(2);
    // gPad->SetGridy();
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
    hDoubleRatioStat->GetYaxis()->SetRangeUser(0.32, 1.75);
    hDoubleRatioStat->GetXaxis()->SetRangeUser(0.5, 6.0);
    hDoubleRatioStat->SetMarkerColor(kRed + 1);
    hDoubleRatioStat->SetLineColor(kRed + 1);
    hDoubleRatioStat->GetYaxis()->SetNdivisions(505);
    hDoubleRatioStat->Draw("pe");
    SetHistoQA(hDoubleRatioSys);
    hDoubleRatioSys->SetFillStyle(0);
    hDoubleRatioSys->SetLineColor(kRed + 1);
    hDoubleRatioSys->SetMarkerColor(kRed + 1);
    hDoubleRatioSys->Draw("e2 same");

    TH1F *hRatio_HM_LM = (TH1F *)hEPOS_UrQMDOFF_KstarK0sRatioHM->Clone("hRatio_HM_LM");
    hRatio_HM_LM->Divide(hEPOS_UrQMDOFF_KstarK0sRatioLM);
    TH1F *hRatio_HM_LM_ON = (TH1F *)hEPOS_UrQMDON_KstarK0sRatioHM->Clone("hRatio_HM_LM_ON");
    hRatio_HM_LM_ON->Divide(hEPOS_UrQMDON_KstarK0sRatioLM);

    TGraph *gRatio_HM_LM_On, *gRatio_HM_LM_Off;
    TGraphErrors *bRatio_HM_LM_On, *bRatio_HM_LM_Off;
    makeCurve(hRatio_HM_LM, kRed + 1, 1, gRatio_HM_LM_Off, bRatio_HM_LM_Off);
    makeCurve(hRatio_HM_LM_ON, kBlue + 1, 7, gRatio_HM_LM_On, bRatio_HM_LM_On);

    if (plotEPOS)
    {
        bRatio_HM_LM_Off->Draw("3 same");
        bRatio_HM_LM_On->Draw("3 same"); // bands
        gRatio_HM_LM_Off->Draw("L same");
        gRatio_HM_LM_On->Draw("L same"); // lines
    }

    TLine *line = new TLine(0.5, 1.0, 6.0, 1.0);
    line->SetLineColor(kGray + 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    // ============================================
    // Pad 3: Significance of deviation from unity
    // ============================================
    cRatio->cd(3);
    // gPad->SetGridy();

    TH1F *hSignificance = (TH1F *)hDoubleRatioStat->Clone("hSignificance");
    // hSignificance->Reset("ICES");
    double maximumSignificance = 0.0;

    for (int i = 1; i <= hDoubleRatioStat->GetNbinsX(); i++)
    {
        double ratio = hDoubleRatioStat->GetBinContent(i);

        double statErr = hDoubleRatioStat->GetBinError(i);
        double sysErr = hDoubleRatioSys->GetBinError(i);

        double totalErr = sqrt(statErr * statErr + sysErr * sysErr);

        double significance = 0.0;

        if (totalErr > 0.0)
            significance = abs(ratio - 1.0) / totalErr;

        hSignificance->SetBinContent(i, significance);
        hSignificance->SetBinError(i, 0.0);

        if (significance > maximumSignificance)
            maximumSignificance = significance;
    }
    cout << "Maximum significance: " << maximumSignificance << endl;

    SetHistoQA(hSignificance);
    hSignificance->SetTitle("");
    hSignificance->SetLineColor(kRed + 1);

    hSignificance->GetYaxis()->SetTitle("n#sigma");
    hSignificance->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");

    hSignificance->GetXaxis()->SetRangeUser(0.0, 10.0);
    hSignificance->GetYaxis()->SetTitleSize(0.04 / 0.20);
    hSignificance->GetYaxis()->SetLabelSize(0.035 / 0.20);
    hSignificance->GetXaxis()->SetTitleSize(0.04 / 0.20);
    hSignificance->GetXaxis()->SetLabelSize(0.035 / 0.20);

    hSignificance->GetXaxis()->SetRangeUser(0.5, 6.0);
    hSignificance->GetYaxis()->SetTitleOffset(1.60 * 0.20);
    hSignificance->GetYaxis()->SetNdivisions(505);

    hSignificance->SetMinimum(0.0);
    hSignificance->SetMaximum(13.0);
    hSignificance->Draw("HIST");

    TLine *phiLine = new TLine(0.5, 3.58291, 6.0, 3.58291);
    phiLine->SetLineColor(kGray + 1);
    phiLine->SetLineStyle(2);
    phiLine->SetLineWidth(3);
    phiLine->Draw("same");

    cRatio->SaveAs("Plots/pTDiff_Ratios/PhiK0sRatio_significance.png");
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);

    c->Divide(1, 3, 0, 0);

    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    TPad *pad3 = (TPad *)c->GetPad(3);

    // Relative pad heights
    pad1Size = 0.60;
    pad2Size = 0.20;
    double pad3Size = 0.20;

    // Set pad positions: (x1, y1, x2, y2)
    pad1->SetPad(0, pad2Size + pad3Size, 1, 1);
    pad2->SetPad(0, pad3Size, 1, pad2Size + pad3Size);
    pad3->SetPad(0, 0, 1, pad3Size);

    // Right margins
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad3->SetRightMargin(0.06);

    // Left margins
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad3->SetLeftMargin(0.16);

    // Top margins
    pad1->SetTopMargin(0.02);
    pad2->SetTopMargin(0.001);
    pad3->SetTopMargin(0.00);

    // Bottom margins
    pad1->SetBottomMargin(0.001);
    pad2->SetBottomMargin(0.00);
    pad3->SetBottomMargin(0.20);

    // Ticks on all pads
    pad1->SetTickx(1);
    pad1->SetTicky(1);

    pad2->SetTickx(1);
    pad2->SetTicky(1);

    pad3->SetTickx(1);
    pad3->SetTicky(1);
}
