#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <cmath>

void compare_systematics()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    // Open both ROOT files
    TFile *f1 = TFile::Open("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots/SysUncert.root", "READ");
    TFile *f2 = TFile::Open("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots/SysUncert_newWidth.root", "READ");

    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie())
    {
        std::cerr << "Error: Could not open input files!" << std::endl;
        return;
    }

    // Define multiplicity boundaries (10 classes + 1 Minimum Bias)
    std::vector<std::pair<int, int>> multBins = {
        {0, 1}, {1, 5}, {5, 10}, {10, 15}, {15, 20}, {20, 30}, {30, 40}, {40, 50}, {50, 70}, {70, 100}, {0, 100} // Minimum Bias
    };

    // 4 Separated Canvases (4x3 grid each)
    TCanvas *cSigExtRaw = new TCanvas("cSigExtRaw", "Signal Extraction (Raw): WidthVar1 vs WidthVarNew", 1600, 1200);
    cSigExtRaw->Divide(4, 3);

    TCanvas *cSigExtSmooth = new TCanvas("cSigExtSmooth", "Signal Extraction (Smoothed): WidthVar1 vs WidthVarNew", 1600, 1200);
    cSigExtSmooth->Divide(4, 3);

    TCanvas *cTotalSysRaw = new TCanvas("cTotalSysRaw", "Total Systematic (Raw): WidthVar1 vs WidthVarNew", 1600, 1200);
    cTotalSysRaw->Divide(4, 3);

    TCanvas *cTotalSysSmooth = new TCanvas("cTotalSysSmooth", "Total Systematic (Smoothed): WidthVar1 vs WidthVarNew", 1600, 1200);
    cTotalSysSmooth->Divide(4, 3);

    // Lambda function to compute quadrature total systematic
    auto calcTotalSys = [](TH1D *hTot, TH1D *hSig, TH1D *hTrk, TH1D *hPID, TH1D *hMat, TH1D *hHad)
    {
        for (int b = 1; b <= hTot->GetNbinsX(); ++b)
        {
            double e1 = hSig ? hSig->GetBinContent(b) : 0;
            double e2 = hTrk ? hTrk->GetBinContent(b) : 0;
            double e3 = hPID ? hPID->GetBinContent(b) : 0;
            double e4 = hMat ? hMat->GetBinContent(b) : 0;
            double e5 = hHad ? hHad->GetBinContent(b) : 0;
            double total = std::sqrt(e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4 + e5 * e5);
            hTot->SetBinContent(b, total);
        }
    };

    for (size_t i = 0; i < multBins.size(); ++i)
    {
        int multLow = multBins[i].first;
        int multHigh = multBins[i].second;

        // Construct histogram names
        TString nameRaw = Form("hSignalExtTotalSys_%d_%d", multLow, multHigh);
        TString nameSmooth = Form("hSignalExtTotalSysSmoothed_%d_%d", multLow, multHigh);

        // Fetch Signal Extraction Histograms
        TH1D *hSigExt_WV1_Raw = (TH1D *)f1->Get(nameRaw);
        TH1D *hSigExt_WV1_Smooth = (TH1D *)f1->Get(nameSmooth);
        TH1D *hSigExt_WVNew_Raw = (TH1D *)f2->Get(nameRaw);
        TH1D *hSigExt_WVNew_Smooth = (TH1D *)f2->Get(nameSmooth);

        // Fetch other systematics for WidthVar1
        TH1D *hTrack_WV1_Raw = (TH1D *)f1->Get(Form("hTrackSelTotalSys_%d_%d", multLow, multHigh));
        TH1D *hPID_WV1_Raw = (TH1D *)f1->Get(Form("hPIDTotalSys_%d_%d", multLow, multHigh));
        TH1D *hMat_WV1_Raw = (TH1D *)f1->Get(Form("hMaterialBudgetTotalSys_%d_%d", multLow, multHigh));
        TH1D *hHad_WV1_Raw = (TH1D *)f1->Get(Form("hHadronicInteractionTotalSys_%d_%d", multLow, multHigh));

        TH1D *hTrack_WV1_Smooth = (TH1D *)f1->Get(Form("hTrackSelTotalSysSmoothed_%d_%d", multLow, multHigh));
        TH1D *hPID_WV1_Smooth = (TH1D *)f1->Get(Form("hPIDTotalSysSmoothed_%d_%d", multLow, multHigh));
        TH1D *hMat_WV1_Smooth = (TH1D *)f1->Get(Form("hMaterialBudgetTotalSysSmoothed_%d_%d", multLow, multHigh));
        TH1D *hHad_WV1_Smooth = (TH1D *)f1->Get(Form("hHadronicInteractionTotalSysSmoothed_%d_%d", multLow, multHigh));

        // Fetch other systematics for WidthVarNew
        TH1D *hTrack_WVNew_Raw = (TH1D *)f2->Get(Form("hTrackSelTotalSys_%d_%d", multLow, multHigh));
        TH1D *hPID_WVNew_Raw = (TH1D *)f2->Get(Form("hPIDTotalSys_%d_%d", multLow, multHigh));
        TH1D *hMat_WVNew_Raw = (TH1D *)f2->Get(Form("hMaterialBudgetTotalSys_%d_%d", multLow, multHigh));
        TH1D *hHad_WVNew_Raw = (TH1D *)f2->Get(Form("hHadronicInteractionTotalSys_%d_%d", multLow, multHigh));

        TH1D *hTrack_WVNew_Smooth = (TH1D *)f2->Get(Form("hTrackSelTotalSysSmoothed_%d_%d", multLow, multHigh));
        TH1D *hPID_WVNew_Smooth = (TH1D *)f2->Get(Form("hPIDTotalSysSmoothed_%d_%d", multLow, multHigh));
        TH1D *hMat_WVNew_Smooth = (TH1D *)f2->Get(Form("hMaterialBudgetTotalSysSmoothed_%d_%d", multLow, multHigh));
        TH1D *hHad_WVNew_Smooth = (TH1D *)f2->Get(Form("hHadronicInteractionTotalSysSmoothed_%d_%d", multLow, multHigh));

        if (!hSigExt_WV1_Raw || !hSigExt_WV1_Smooth || !hSigExt_WVNew_Raw || !hSigExt_WVNew_Smooth)
        {
            std::cerr << "Warning: Missing histograms for class " << multLow << "-" << multHigh << "%" << std::endl;
            continue;
        }

        // ==========================================
        // 1. RAW SIGNAL EXTRACTION CANVAS
        // ==========================================
        cSigExtRaw->cd(i + 1);
        hSigExt_WV1_Raw->SetLineColor(kRed);
        hSigExt_WV1_Raw->SetLineWidth(2);
        hSigExt_WV1_Raw->SetTitle(Form("Raw SigExt: Mult %d-%d%%;p_{T} (GeV/c);Uncertainty", multLow, multHigh));

        hSigExt_WVNew_Raw->SetLineColor(kBlue);
        hSigExt_WVNew_Raw->SetLineWidth(2);

        hSigExt_WV1_Raw->Draw("HIST");
        hSigExt_WVNew_Raw->Draw("HIST SAME");

        if (i == 0)
        {
            TLegend *leg1 = new TLegend(0.45, 0.70, 0.88, 0.88);
            leg1->SetBorderSize(0);
            leg1->AddEntry(hSigExt_WV1_Raw, "WidthVar1 (Raw)", "l");
            leg1->AddEntry(hSigExt_WVNew_Raw, "WidthVarNew (Raw)", "l");
            leg1->Draw();
        }

        // ==========================================
        // 2. SMOOTHED SIGNAL EXTRACTION CANVAS
        // ==========================================
        cSigExtSmooth->cd(i + 1);
        hSigExt_WV1_Smooth->SetLineColor(kRed);
        hSigExt_WV1_Smooth->SetLineWidth(2);
        hSigExt_WV1_Smooth->SetTitle(Form("Smoothed SigExt: Mult %d-%d%%;p_{T} (GeV/c);Uncertainty", multLow, multHigh));

        hSigExt_WVNew_Smooth->SetLineColor(kBlue);
        hSigExt_WVNew_Smooth->SetLineWidth(2);

        hSigExt_WV1_Smooth->Draw("HIST");
        hSigExt_WVNew_Smooth->Draw("HIST SAME");

        if (i == 0)
        {
            TLegend *leg2 = new TLegend(0.25, 0.70, 0.88, 0.88);
            leg2->SetBorderSize(0);
            leg2->SetTextSize(0.06);
            leg2->AddEntry(hSigExt_WV1_Smooth, "WidthVar Old", "l");
            leg2->AddEntry(hSigExt_WVNew_Smooth, "WidthVar New", "l");
            leg2->Draw();
        }

        // --- Calculate Total Systematics ---
        TH1D *hTot_WV1_Raw = (TH1D *)hSigExt_WV1_Raw->Clone(Form("hTot_WV1_Raw_%d_%d", multLow, multHigh));
        TH1D *hTot_WV1_Smooth = (TH1D *)hSigExt_WV1_Smooth->Clone(Form("hTot_WV1_Smooth_%d_%d", multLow, multHigh));
        TH1D *hTot_WVNew_Raw = (TH1D *)hSigExt_WVNew_Raw->Clone(Form("hTot_WVNew_Raw_%d_%d", multLow, multHigh));
        TH1D *hTot_WVNew_Smooth = (TH1D *)hSigExt_WVNew_Smooth->Clone(Form("hTot_WVNew_Smooth_%d_%d", multLow, multHigh));

        calcTotalSys(hTot_WV1_Raw, hSigExt_WV1_Raw, hTrack_WV1_Raw, hPID_WV1_Raw, hMat_WV1_Raw, hHad_WV1_Raw);
        calcTotalSys(hTot_WV1_Smooth, hSigExt_WV1_Smooth, hTrack_WV1_Smooth, hPID_WV1_Smooth, hMat_WV1_Smooth, hHad_WV1_Smooth);
        calcTotalSys(hTot_WVNew_Raw, hSigExt_WVNew_Raw, hTrack_WVNew_Raw, hPID_WVNew_Raw, hMat_WVNew_Raw, hHad_WVNew_Raw);
        calcTotalSys(hTot_WVNew_Smooth, hSigExt_WVNew_Smooth, hTrack_WVNew_Smooth, hPID_WVNew_Smooth, hMat_WVNew_Smooth, hHad_WVNew_Smooth);

        // ==========================================
        // 3. RAW TOTAL SYSTEMATIC CANVAS
        // ==========================================
        cTotalSysRaw->cd(i + 1);
        hTot_WV1_Raw->SetLineColor(kRed);
        hTot_WV1_Raw->SetLineWidth(2);
        hTot_WV1_Raw->SetTitle(Form("Raw Total Sys: Mult %d-%d%%;p_{T} (GeV/c);Total Uncertainty", multLow, multHigh));

        hTot_WVNew_Raw->SetLineColor(kBlue);
        hTot_WVNew_Raw->SetLineWidth(2);

        hTot_WV1_Raw->Draw("HIST");
        hTot_WVNew_Raw->Draw("HIST SAME");

        if (i == 0)
        {
            TLegend *leg3 = new TLegend(0.45, 0.70, 0.88, 0.88);
            leg3->SetBorderSize(0);
            leg3->AddEntry(hTot_WV1_Raw, "WidthVar1 (Raw)", "l");
            leg3->AddEntry(hTot_WVNew_Raw, "WidthVarNew (Raw)", "l");
            leg3->Draw();
        }

        // ==========================================
        // 4. SMOOTHED TOTAL SYSTEMATIC CANVAS
        // ==========================================
        cTotalSysSmooth->cd(i + 1);
        hTot_WV1_Smooth->SetLineColor(kRed);
        hTot_WV1_Smooth->SetLineWidth(2);
        hTot_WV1_Smooth->SetTitle(Form("Smoothed Total Sys: Mult %d-%d%%;p_{T} (GeV/c);Total Uncertainty", multLow, multHigh));

        hTot_WVNew_Smooth->SetLineColor(kBlue);
        hTot_WVNew_Smooth->SetLineWidth(2);

        hTot_WV1_Smooth->Draw("HIST");
        hTot_WVNew_Smooth->Draw("HIST SAME");

        if (i == 0)
        {
            TLegend *leg4 = new TLegend(0.45, 0.70, 0.88, 0.88);
            leg4->SetBorderSize(0);
            leg4->AddEntry(hTot_WV1_Smooth, "WidthVar1 (Smoothed)", "l");
            leg4->AddEntry(hTot_WVNew_Smooth, "WidthVarNew (Smoothed)", "l");
            leg4->Draw();
        }
    }

    // Optional: Save files
    // cSigExtRaw->SaveAs("SignalExtraction_Raw_Comparison.png");
    cSigExtSmooth->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots/CompareSystematics/SignalExtraction_Smoothed_Comparison.png");
    // cTotalSysRaw->SaveAs("TotalSystematic_Raw_Comparison.png");
    cTotalSysSmooth->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots/CompareSystematics/TotalSystematic_Smoothed_Comparison.png");
}