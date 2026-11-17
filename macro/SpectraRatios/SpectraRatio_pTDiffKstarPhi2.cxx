#include <iostream>
#include <cmath>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../spectra/YieldMean.C"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size, double &pad3Size);

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

void SpectraRatio_pTDiffKstarPhi2()
{
    bool plotEPOS = false;
    gStyle->SetOptStat(0);

    TFile *fKstar = OpenFile("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/KstarSpectra.root");
    TFile *fPhi = OpenFile("PhiResults/PhiAnalysisResults_pprun3.root");

    TH1F *hKstar_statHM = (TH1F *)fKstar->Get("SpectraStat_0_1");
    TH1F *hKstar_sysHM = (TH1F *)fKstar->Get("SpectraSys_0_1");
    TH1F *hKstar_sysUncorrHM = (TH1F *)fKstar->Get("SpectraUncorrSys_0_1");

    TH1F *hKstar_statLM = (TH1F *)fKstar->Get("SpectraStat_70_100");
    TH1F *hKstar_sysLM = (TH1F *)fKstar->Get("SpectraSys_70_100");
    TH1F *hKstar_sysUncorrLM = (TH1F *)fKstar->Get("SpectraUncorrSys_70_100");

    TH1F *hPhi_statHM = (TH1F *)fPhi->Get("Spectra/01/hStat");
    TH1F *hPhi_sysHM = (TH1F *)fPhi->Get("Spectra/01/hSys");
    TH1F *hPhi_sysUncorrHM = (TH1F *)fPhi->Get("Spectra/01/hSys"); // Does not have uncorrelated systematic, so using total systematic

    TH1F *hPhi_statLM = (TH1F *)fPhi->Get("Spectra/10/hStat");
    TH1F *hPhi_sysLM = (TH1F *)fPhi->Get("Spectra/10/hSys");
    TH1F *hPhi_sysUncorrLM = (TH1F *)fPhi->Get("Spectra/10/hSys"); // Does not have uncorrelated systematic, so using total systematic

    cout << "Number of kaon bins in HM: " << hPhi_statHM->GetNbinsX() << endl;
    cout << "Number of kaon bins in LM: " << hPhi_statLM->GetNbinsX() << endl;

    int TotalBinsKaon = hPhi_statHM->GetNbinsX();

    TF1 *fitFcn1 = new TF1("fitfunc", FuncLavy, 0.0, 10.0, 4);
    fitFcn1->SetParameter(0, 6.6);
    fitFcn1->SetParameter(1, 0.2);
    fitFcn1->SetParameter(2, 1.019);
    fitFcn1->SetParameter(3, 0.18);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn1->SetLineColor(kRed + 1);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);

    TF1 *fitFcn2 = new TF1("fitfunc2", FuncLavy, 0.0, 10.0, 4);
    fitFcn2->SetParameter(0, 6.6);
    fitFcn2->SetParameter(1, 0.2);
    fitFcn2->SetParameter(2, 1.019);
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

    TH1 *hout = YieldMean(hPhi_statHM, hPhi_sysHM, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *hout2 = YieldMean(hPhi_statLM, hPhi_sysLM, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

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
    SetHistoQA(hPhi_statHM);
    gPad->SetLogy(1);
    hPhi_statHM->SetMarkerColor(kBlue + 1);
    hPhi_statHM->SetLineColor(kBlue + 1);
    hPhi_statHM->SetMarkerSize(1.3);
    hPhi_statHM->Draw("pe same");
    fitFcn1->Draw("same");

    TLegend *legend = new TLegend(0.7, 0.7, 0.80, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetHeader("FT0M: 0-1%");
    legend->AddEntry(hKstar_statHM, "K*(892)^{0}", "p");
    legend->AddEntry(hPhi_statHM, "K^{0}_{S}", "p");
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
    SetHistoQA(hPhi_statLM);
    gPad->SetLogy(1);
    hPhi_statLM->SetMarkerColor(kBlue + 1);
    hPhi_statLM->SetLineColor(kBlue + 1);
    hPhi_statLM->SetMarkerSize(1.3);
    hPhi_statLM->Draw("pe same");
    fitFcn2->Draw("same");

    TLegend *legend2 = new TLegend(0.7, 0.7, 0.80, 0.92);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextFont(42);
    legend2->SetTextSize(0.03);
    legend2->SetHeader("FT0M: 70-100%");
    legend2->AddEntry(hKstar_statLM, "K*(892)^{0}", "p");
    legend2->AddEntry(hPhi_statLM, "K^{0}_{S}", "p");
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
        double kaonYieldHM = getKaonYield(hPhi_statHM, fitFcn1, xLow, xHigh);
        double kaonYieldLM = getKaonYield(hPhi_statLM, fitFcn2, xLow, xHigh);

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

        double kaon_stat_errHM = getKaonErr(hPhi_statHM, xLow, xHigh);
        double rel_kstar_statHM = kstar_stat_errHM / kstarYieldHM;
        double rel_kaon_statHM = kaon_stat_errHM / kaonYieldHM;

        double kaon_stat_errLM = getKaonErr(hPhi_statLM, xLow, xHigh);
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

        double kaon_sys_errHM = getKaonErr(hPhi_sysHM, xLow, xHigh);
        double rel_kstar_sysHM = kstar_sys_errHM / kstarYieldHM;
        double rel_kaon_sysHM = kaon_sys_errHM / kaonYieldHM;

        double kaon_sys_errLM = getKaonErr(hPhi_sysLM, xLow, xHigh);
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

        double kaon_sysUncorr_errHM = getKaonErr(hPhi_sysUncorrHM, xLow, xHigh);
        double rel_kaon_sysUncorrHM = kaon_sysUncorr_errHM / kaonYieldHM;

        double kaon_sysUncorr_errLM = getKaonErr(hPhi_sysUncorrLM, xLow, xHigh);
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
    double pad1Size, pad2Size, pad3Size;
    canvas_style(cRatio, pad1Size, pad2Size, pad3Size);
    cRatio->cd(1);

    SetHistoQA(hKstarK0sRatioHM_stat);
    hKstarK0sRatioHM_stat->GetYaxis()->SetTitle("K^{*0} / #Phi");
    hKstarK0sRatioHM_stat->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hKstarK0sRatioHM_stat->SetMarkerStyle(20);
    hKstarK0sRatioHM_stat->SetMarkerSize(1.45);
    hKstarK0sRatioHM_stat->GetYaxis()->SetTitleOffset(1.60 * pad1Size);
    hKstarK0sRatioHM_stat->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    // hKstarK0sRatioHM_stat->GetYaxis()->SetRangeUser(-2.8, 9.25);
    hKstarK0sRatioHM_stat->GetYaxis()->SetRangeUser(0.8, 5.25);
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
    // TFile *fEPOS_spectra = OpenFile("/home/sawan/Storage/EPOS4_Resonance/yield_outfile_latest2.root");
    // TFile *fEPOS_spectra = OpenFile("/home/sawan/Storage/EPOS4_Resonance/yield_outfile_light_oct2026_fixed.root");
    TFile *fEPOS_spectra = OpenFile("/home/sawan/Storage/EPOS4_Resonance/yield_outfile_light_Bongfix4.root");

    TH1F *hEPOS_UrQMDOFF_KstarHM = (TH1F *)fEPOS_spectra->Get("Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent9");
    TH1F *hEPOS_UrQMDOFF_KstarLM = (TH1F *)fEPOS_spectra->Get("Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent0");
    TH1F *hEPOS_UrQMDON_KstarHM = (TH1F *)fEPOS_spectra->Get("Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent9");
    TH1F *hEPOS_UrQMDON_KstarLM = (TH1F *)fEPOS_spectra->Get("Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent0");

    TH1F *hEPOS_UrQMDOFF_PhiHM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent9");
    TH1F *hEPOS_UrQMDOFF_PhiLM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent0");
    TH1F *hEPOS_UrQMDON_PhiHM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent9");
    TH1F *hEPOS_UrQMDON_PhiLM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent0");

    // Correct for branching ratio of phi for UrQMD ON case only
    double phiBR = 0.415; // PDG 2022
    hEPOS_UrQMDON_PhiHM->Scale(1.0 / phiBR);
    hEPOS_UrQMDON_PhiLM->Scale(1.0 / phiBR);

    hEPOS_UrQMDON_KstarHM->Scale(1.0 / 0.667); // Correct for branching ratio of K* -> K pi
    hEPOS_UrQMDON_KstarLM->Scale(1.0 / 0.667); // Correct for branching ratio of K* -> K pi

    // ============================================================
    // EPOS K*/Phi ratios in the ORIGINAL pT bins (no rebinning)
    // ============================================================

    auto CalculateRatio = [](TH1F *hKstar, TH1F *hPhi, const char *name) -> TH1F *
    {
        if (!hKstar || !hPhi)
        {
            cout << "Error: Missing EPOS input histogram for " << name << endl;
            return nullptr;
        }

        if (hKstar->GetNbinsX() != hPhi->GetNbinsX())
        {
            cout << "Error: Kstar and Phi have different numbers of bins for " << name << endl;
            return nullptr;
        }

        const int nBins = hKstar->GetNbinsX();
        for (int ibin = 1; ibin <= nBins; ++ibin)
        {
            const double kLow = hKstar->GetXaxis()->GetBinLowEdge(ibin);
            const double kHigh = hKstar->GetXaxis()->GetBinUpEdge(ibin);
            const double pLow = hPhi->GetXaxis()->GetBinLowEdge(ibin);
            const double pHigh = hPhi->GetXaxis()->GetBinUpEdge(ibin);

            if (fabs(kLow - pLow) > 1e-6 || fabs(kHigh - pHigh) > 1e-6)
            {
                cout << "Error: Kstar and Phi pT bin edges differ for " << name << ", bin " << ibin << endl;
                return nullptr;
            }
        }

        TH1F *hRatio = (TH1F *)hKstar->Clone(name);
        hRatio->SetDirectory(nullptr);
        hRatio->Reset("ICES");

        for (int ibin = 1; ibin <= nBins; ++ibin)
        {
            const double kstar = hKstar->GetBinContent(ibin);
            const double kstarErr = hKstar->GetBinError(ibin);
            const double phi = hPhi->GetBinContent(ibin);
            const double phiErr = hPhi->GetBinError(ibin);

            if (kstar <= 0.0 || phi <= 0.0)
                continue;

            // Keep the factor 0.5 as in the original analysis.
            const double ratio = 0.5 * kstar / phi;
            const double relKstar = kstarErr / kstar;
            const double relPhi = phiErr / phi;
            const double ratioErr = ratio * sqrt(relKstar * relKstar + relPhi * relPhi);

            hRatio->SetBinContent(ibin, ratio);
            hRatio->SetBinError(ibin, ratioErr);
        }

        return hRatio;
    };

    // Build a line through the original bin centres and a shaded statistical band.
    auto makeCurve = [](TH1 *h, Color_t col, Style_t lstyle, TGraph *&gLine, TGraphErrors *&gBand)
    {
        gLine = new TGraph();
        gBand = new TGraphErrors();

        if (!h)
            return;

        for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin)
        {
            const double y = h->GetBinContent(ibin);
            if (y <= 0.0)
                continue;

            const double x = h->GetXaxis()->GetBinCenter(ibin);
            const double ex = 0.5 * h->GetXaxis()->GetBinWidth(ibin);
            const double ey = h->GetBinError(ibin);

            const int ip = gLine->GetN();
            gLine->SetPoint(ip, x, y);

            const int ib = gBand->GetN();
            gBand->SetPoint(ib, x, y);
            gBand->SetPointError(ib, ex, ey);
        }

        gLine->SetLineColor(col);
        gLine->SetLineWidth(3);
        gLine->SetLineStyle(lstyle);

        gBand->SetFillColorAlpha(col, 0.20);
        gBand->SetLineColor(col);
        gBand->SetLineWidth(1);
    };

    // -------------------------
    // UrQMD OFF
    // -------------------------
    TH1F *hEPOS_UrQMDOFF_KstarK0sRatioHM = CalculateRatio(hEPOS_UrQMDOFF_KstarHM, hEPOS_UrQMDOFF_PhiHM, "hEPOS_UrQMDOFF_KstarPhiRatioHM");
    TH1F *hEPOS_UrQMDOFF_KstarK0sRatioLM = CalculateRatio(hEPOS_UrQMDOFF_KstarLM, hEPOS_UrQMDOFF_PhiLM, "hEPOS_UrQMDOFF_KstarPhiRatioLM");

    // -------------------------
    // UrQMD ON
    // -------------------------
    TH1F *hEPOS_UrQMDON_KstarK0sRatioHM = CalculateRatio(hEPOS_UrQMDON_KstarHM, hEPOS_UrQMDON_PhiHM, "hEPOS_UrQMDON_KstarPhiRatioHM");
    TH1F *hEPOS_UrQMDON_KstarK0sRatioLM = CalculateRatio(hEPOS_UrQMDON_KstarLM, hEPOS_UrQMDON_PhiLM, "hEPOS_UrQMDON_KstarPhiRatioLM");

    if (!hEPOS_UrQMDOFF_KstarK0sRatioHM || !hEPOS_UrQMDOFF_KstarK0sRatioLM || !hEPOS_UrQMDON_KstarK0sRatioHM || !hEPOS_UrQMDON_KstarK0sRatioLM)
    {
        cout << "Error: Could not create all EPOS Kstar/Phi ratio histograms. Check the input histograms and their pT binning." << endl;
        return;
    }

    // Keep the histograms' bin contents/errors untouched; graphs are used for drawing.
    TGraph *gHM = nullptr, *gLM = nullptr, *gHM_UrQMDON = nullptr, *gLM_UrQMDON = nullptr;
    TGraphErrors *bHM = nullptr, *bLM = nullptr, *bHM_UrQMDON = nullptr, *bLM_UrQMDON = nullptr;

    makeCurve(hEPOS_UrQMDOFF_KstarK0sRatioHM, kRed + 1, 1, gHM, bHM);
    makeCurve(hEPOS_UrQMDOFF_KstarK0sRatioLM, kBlue + 1, 1, gLM, bLM);
    makeCurve(hEPOS_UrQMDON_KstarK0sRatioHM, kRed + 1, 7, gHM_UrQMDON, bHM_UrQMDON);
    makeCurve(hEPOS_UrQMDON_KstarK0sRatioLM, kBlue + 1, 7, gLM_UrQMDON, bLM_UrQMDON);

    if (plotEPOS)
    {
        // Draw uncertainty bands first, then the corresponding central-value lines.
        bHM->Draw("3 same");
        bLM->Draw("3 same");
        bHM_UrQMDON->Draw("3 same");
        bLM_UrQMDON->Draw("3 same");

        gHM->Draw("L same");
        gLM->Draw("L same");
        gHM_UrQMDON->Draw("L same");
        gLM_UrQMDON->Draw("L same");
    }

    TLegend *legendRatio = new TLegend(0.5, 0.02, 0.80, 0.22);
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
    hDoubleRatioStat->GetYaxis()->SetRangeUser(0.53, 1.25);
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
    makeCurve(hRatio_HM_LM, kBlue + 1, 1, gRatio_HM_LM_Off, bRatio_HM_LM_Off);
    makeCurve(hRatio_HM_LM_ON, kRed + 1, 7, gRatio_HM_LM_On, bRatio_HM_LM_On);
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

    hSignificance->GetXaxis()->SetRangeUser(0.0, 6.0);
    hSignificance->GetYaxis()->SetTitleSize(0.04 / pad3Size);
    hSignificance->GetYaxis()->SetLabelSize(0.035 / pad3Size);
    hSignificance->GetXaxis()->SetTitleSize(0.04 / pad3Size);
    hSignificance->GetXaxis()->SetLabelSize(0.035 / pad3Size);

    hSignificance->GetXaxis()->SetRangeUser(0.0, 6.0);
    hSignificance->GetYaxis()->SetTitleOffset(1.60 * pad3Size);
    hSignificance->GetYaxis()->SetNdivisions(505);

    hSignificance->SetMinimum(0.0);
    hSignificance->SetMaximum(4.8);

    hSignificance->Draw("HIST");

    // print the values of significance for each bin
    for (int i = 1; i <= hSignificance->GetNbinsX(); i++)
    {
        cout << "Bin " << i << ", pT = " << hSignificance->GetXaxis()->GetBinCenter(i) << ", Significance = " << hSignificance->GetBinContent(i) << endl;
    }

    cRatio->SaveAs("Plots/pTDiff_Ratios/KstarPhiRatio.png");
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size, double &pad3Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);

    c->Divide(1, 3, 0, 0);

    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    TPad *pad3 = (TPad *)c->GetPad(3);

    // Relative pad heights
    pad1Size = 0.60;
    pad2Size = 0.15;
    pad3Size = 0.25;

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
    pad3->SetBottomMargin(0.43);

    // Ticks on all pads
    pad1->SetTickx(1);
    pad1->SetTicky(1);

    pad2->SetTickx(1);
    pad2->SetTicky(1);

    pad3->SetTickx(1);
    pad3->SetTicky(1);
}
