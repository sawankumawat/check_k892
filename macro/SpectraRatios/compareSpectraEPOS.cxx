#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"

// Include custom style headers if present in your environment
#include "../src/style.h"

using namespace std;

TFile *OpenFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (!f || f->IsZombie())
    {
        cout << "Error: File not found or corrupted: " << path << endl;
        return nullptr;
    }
    return f;
}

// Convert TH1F to TGraphErrors for smooth EPOS curve drawing
void makeEPOSCurve(TH1F *h, Color_t col, Style_t lstyle, TGraph *&gLine, TGraphErrors *&gBand)
{
    gLine = new TGraph();
    gBand = new TGraphErrors();
    if (!h)
        return;

    int nBins = h->GetNbinsX();
    int first = -1, last = -1;
    for (int i = 1; i <= nBins; ++i)
    {
        if (h->GetBinContent(i) > 0)
        {
            if (first < 0)
                first = i;
            last = i;
        }
    }
    if (first < 0)
        return;

    auto addPoint = [&](double x, int ibin)
    {
        int idxLine = gLine->GetN();
        gLine->SetPoint(idxLine, x, h->GetBinContent(ibin));
        int idxBand = gBand->GetN();
        gBand->SetPoint(idxBand, x, h->GetBinContent(ibin));
        gBand->SetPointError(idxBand, 0, h->GetBinError(ibin));
    };

    addPoint(h->GetXaxis()->GetXmin(), first);
    for (int i = first; i <= last; ++i)
    {
        if (h->GetBinContent(i) > 0)
            addPoint(h->GetBinCenter(i), i);
    }
    addPoint(h->GetXaxis()->GetXmax(), last);

    gLine->SetLineColor(col);
    gLine->SetLineWidth(2);
    gLine->SetLineStyle(lstyle);

    gBand->SetFillColorAlpha(col, 0.25);
    gBand->SetLineColor(col);
}

void compareSpectraEPOS()
{
    gStyle->SetOptStat(0);

    // =========================================================================
    // 1. INPUT FILES
    // =========================================================================
    TFile *fKstar = OpenFile("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/KstarSpectra.root");
    TFile *fPhi = OpenFile("PhiResults/PhiAnalysisResults_pprun3.root");
    TFile *fEPOS = OpenFile("/home/sawan/Storage/EPOS4_Resonance/yield_outfile_light_oct2026_fixed.root");

    if (!fKstar || !fPhi || !fEPOS)
        return;

    // Multiplicity classes setup
    // 10 classes + Minimum Bias (0-100%)
    vector<string> multLabels = {
        "0-1%", "1-5%", "5-10%", "10-15%", "15-20%", "20-30%",
        "30-40%", "40-50%", "50-70%", "70-100%", "MB (0-100%)"};

    // Paths inside ROOT files
    vector<string> kstarHistNames = {
        "SpectraStat_0_1", "SpectraStat_1_5", "SpectraStat_5_10", "SpectraStat_10_15", "SpectraStat_15_20", "SpectraStat_20_30",
        "SpectraStat_30_40", "SpectraStat_40_50", "SpectraStat_50_70", "SpectraStat_70_100", "SpectraStat_0_100"};
    vector<string> kstarSysNames = {
        "SpectraSys_0_1", "SpectraSys_1_5", "SpectraSys_5_10", "SpectraSys_10_15", "SpectraSys_15_20", "SpectraSys_20_30",
        "SpectraSys_30_40", "SpectraSys_40_50", "SpectraSys_50_70", "SpectraSys_70_100", "SpectraSys_0_100"};

    vector<string> phiHistNames = {
        "Spectra/01/hStat", "Spectra/02/hStat", "Spectra/03/hStat", "Spectra/04/hStat", "Spectra/05/hStat",
        "Spectra/06/hStat", "Spectra/07/hStat", "Spectra/08/hStat", "Spectra/09/hStat", "Spectra/10/hStat", "Spectra/MB/hStat"};
    vector<string> phiSysNames = {
        "Spectra/01/hSys", "Spectra/02/hSys", "Spectra/03/hSys", "Spectra/04/hSys", "Spectra/05/hSys",
        "Spectra/06/hSys", "Spectra/07/hSys", "Spectra/08/hSys", "Spectra/09/hSys", "Spectra/10/hSys", "Spectra/MB/hSys"};

    // EPOS Cent indices (Cent9 = 0-1%, Cent0 = 70-100%)
    vector<string> eposOffKstar = {
        "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent9", "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent8",
        "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent7", "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent6",
        "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent5", "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent4",
        "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent3", "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent2",
        "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent1", "Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent0", "Kstar_MB_Spectra/hKstar_Pt_MB_OFF"};
    vector<string> eposOnKstar = {
        "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent9", "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent8",
        "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent7", "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent6",
        "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent5", "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent4",
        "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent3", "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent2",
        "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent1", "Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent0", "Kstar_MB_Spectra/hKstar_Pt_MB_ON_reg_res"};

    vector<string> eposOffPhi = {
        "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent9", "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent8",
        "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent7", "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent6",
        "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent5", "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent4",
        "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent3", "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent2",
        "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent1", "Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent0", "Phi_MB_Spectra/hPhi_Pt_MB_OFF"};
    vector<string> eposOnPhi = {
        "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent9", "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent8",
        "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent7", "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent6",
        "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent5", "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent4",
        "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent3", "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent2",
        "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent1", "Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent0", "Phi_MB_Spectra/hPhi_Pt_MB_ON_reg_res"};

    // Color palette for distinction across multiplicity bins
    Color_t colors[] = {kRed + 1, kOrange + 7, kYellow + 2, kGreen + 2, kCyan + 2,
                        kBlue + 1, kViolet + 1, kMagenta + 2, kGray + 2, kBlack, kBrown + 2};

    // =========================================================================
    // 2. CANVAS 1: K* SPECTRA COMPARISON
    // =========================================================================
    TCanvas *cKstar = new TCanvas("cKstar", "Kstar Spectra Multiplicity Comparison", 800, 900);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);

    TLegend *legKstar = new TLegend(0.55, 0.65, 0.92, 0.94);
    legKstar->SetBorderSize(0);
    legKstar->SetFillStyle(0);
    legKstar->SetNColumns(2);
    legKstar->SetTextSize(0.02);

    bool firstDraw = true;

    for (size_t i = 0; i < multLabels.size(); ++i)
    {
        TH1F *hStat_raw = (TH1F *)fKstar->Get(kstarHistNames[i].c_str());
        TH1F *hSys_raw  = (TH1F *)fKstar->Get(kstarSysNames[i].c_str());

        TH1F *hEPOS_OFF_raw = (TH1F *)fEPOS->Get(eposOffKstar[i].c_str());
        TH1F *hEPOS_ON_raw  = (TH1F *)fEPOS->Get(eposOnKstar[i].c_str());

        if (!hStat_raw) continue;

        TH1F *hStat     = (TH1F *)hStat_raw->Clone(Form("hStat_Kstar_%zu", i));
        TH1F *hSys      = hSys_raw ? (TH1F *)hSys_raw->Clone(Form("hSys_Kstar_%zu", i)) : nullptr;
        TH1F *hEPOS_OFF = hEPOS_OFF_raw ? (TH1F *)hEPOS_OFF_raw->Clone(Form("hEPOS_OFF_Kstar_%zu", i)) : nullptr;
        TH1F *hEPOS_ON  = hEPOS_ON_raw ? (TH1F *)hEPOS_ON_raw->Clone(Form("hEPOS_ON_Kstar_%zu", i)) : nullptr;

        // Apply visual scaling to separate spectra visually
        double scaleFactor = pow(2.0, (int)multLabels.size() - 1 - (int)i);
        hStat->Scale(scaleFactor);
        if (hSys)
            hSys->Scale(scaleFactor);

        // Scale EPOS K* by 0.5 along with the visual separation scaleFactor
        if (hEPOS_OFF)
            hEPOS_OFF->Scale(0.5 * scaleFactor);
        else
            cout << "Warning: Couldn't find " << eposOffKstar[i] << endl;

        if (hEPOS_ON)
            hEPOS_ON->Scale(0.5 * scaleFactor);
        else
            cout << "Warning: Couldn't find " << eposOnKstar[i] << endl;

        hStat->SetMarkerStyle(20);
        hStat->SetMarkerSize(1.1);
        hStat->SetMarkerColor(colors[i]);
        hStat->SetLineColor(colors[i]);

        if (firstDraw)
        {
            hStat->SetTitle(";#it{p}_{T} (GeV/#it{c}); d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1} [Scaled]");
            hStat->GetXaxis()->SetRangeUser(0.0, 8.0);
            hStat->GetYaxis()->SetRangeUser(4e-6, 2e5);
            hStat->GetXaxis()->SetTitleSize(0.04);
            hStat->GetYaxis()->SetTitleSize(0.04);
            hStat->GetYaxis()->SetTitleOffset(1.5);
            hStat->Draw("pe");
            firstDraw = false;
        }
        else
        {
            hStat->Draw("pe same");
        }

        if (hSys)
        {
            hSys->SetFillStyle(0);
            hSys->SetLineColor(colors[i]);
            hSys->SetMarkerColor(colors[i]);
            hSys->SetMarkerStyle(21);
            hSys->SetMarkerSize(1.1);
            hSys->Draw("e2 same");
        }

        // Draw EPOS curves if available
        if (hEPOS_OFF)
        {
            TGraph *gLine = nullptr;
            TGraphErrors *gBand = nullptr;
            makeEPOSCurve(hEPOS_OFF, colors[i], 1, gLine, gBand);
            if (gLine)
                gLine->Draw("L same");
        }
        if (hEPOS_ON)
        {
            TGraph *gLine = nullptr;
            TGraphErrors *gBand = nullptr;
            makeEPOSCurve(hEPOS_ON, colors[i], 7, gLine, gBand);
            if (gLine)
                gLine->Draw("L same");
        }

        legKstar->AddEntry(hStat, Form("%s (#times 2^{%d})", multLabels[i].c_str(), (int)(multLabels.size() - 1 - i)), "pe");
    }

    legKstar->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.02);
    lat.DrawLatex(0.20, 0.90, "ALICE Run 3 pp #sqrt{#it{s}} = 13.6 TeV");
    lat.DrawLatex(0.20, 0.86, "K*(892)^{0} Spectrum comparison");
    lat.DrawLatex(0.20, 0.82, "Solid: EPOS OFF | Dashed: EPOS ON");

    cKstar->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/SpectraCompareWithModel/Kstar_Multiplicity_Spectra_EPOS.png");

    // =========================================================================
    // 3. CANVAS 2: PHI SPECTRA COMPARISON
    // =========================================================================
    TCanvas *cPhi = new TCanvas("cPhi", "Phi Spectra Multiplicity Comparison", 800, 900);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);

    TLegend *legPhi = new TLegend(0.55, 0.65, 0.92, 0.94);
    legPhi->SetBorderSize(0);
    legPhi->SetFillStyle(0);
    legPhi->SetTextSize(0.02);
    legPhi->SetNColumns(2);

    firstDraw = true;
    double phiBR = 0.489; // Branching ratio scaling for UrQMD ON case if required

    for (size_t i = 0; i < multLabels.size(); ++i)
    {
        TH1F *hStat_raw = (TH1F *)fPhi->Get(phiHistNames[i].c_str());
        TH1F *hSys_raw  = (TH1F *)fPhi->Get(phiSysNames[i].c_str());

        TH1F *hEPOS_OFF_raw = (TH1F *)fEPOS->Get(eposOffPhi[i].c_str());
        TH1F *hEPOS_ON_raw  = (TH1F *)fEPOS->Get(eposOnPhi[i].c_str());

        if (!hStat_raw) continue;

        TH1F *hStat     = (TH1F *)hStat_raw->Clone(Form("hStat_Phi_%zu", i));
        TH1F *hSys      = hSys_raw ? (TH1F *)hSys_raw->Clone(Form("hSys_Phi_%zu", i)) : nullptr;
        TH1F *hEPOS_OFF = hEPOS_OFF_raw ? (TH1F *)hEPOS_OFF_raw->Clone(Form("hEPOS_OFF_Phi_%zu", i)) : nullptr;
        TH1F *hEPOS_ON  = hEPOS_ON_raw ? (TH1F *)hEPOS_ON_raw->Clone(Form("hEPOS_ON_Phi_%zu", i)) : nullptr;

        double scaleFactor = pow(2.0, (int)multLabels.size() - 1 - (int)i);
        hStat->Scale(scaleFactor);
        if (hSys)
            hSys->Scale(scaleFactor);
        if (hEPOS_OFF)
            hEPOS_OFF->Scale(scaleFactor);
        else
            cout << "Warning: Couldn't find " << eposOffPhi[i] << endl;

        if (hEPOS_ON)
        {
            hEPOS_ON->Scale(1.0 / phiBR); // Scale BR for UrQMD ON
            hEPOS_ON->Scale(scaleFactor);
        }
        else
            cout << "Warning: Couldn't find " << eposOnPhi[i] << endl;

        hStat->SetMarkerStyle(21);
        hStat->SetMarkerSize(1.1);
        hStat->SetMarkerColor(colors[i]);
        hStat->SetLineColor(colors[i]);

        if (firstDraw)
        {
            hStat->SetTitle(";#it{p}_{T} (GeV/#it{c}); d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1} [Scaled]");
            hStat->GetXaxis()->SetRangeUser(0.0, 8.0);
            hStat->GetYaxis()->SetRangeUser(2e-6, 1e5);
            hStat->GetXaxis()->SetTitleSize(0.04);
            hStat->GetYaxis()->SetTitleSize(0.04);
            hStat->GetYaxis()->SetTitleOffset(1.5);
            hStat->Draw("pe");
            firstDraw = false;
        }
        else
        {
            hStat->Draw("pe same");
        }

        if (hSys)
        {
            hSys->SetFillStyle(0);
            hSys->SetLineColor(colors[i]);
            hSys->SetMarkerColor(colors[i]);
            hSys->SetMarkerStyle(21);
            hSys->SetMarkerSize(1.1);
            hSys->Draw("e2 same");
        }

        if (hEPOS_OFF)
        {
            TGraph *gLine = nullptr;
            TGraphErrors *gBand = nullptr;
            makeEPOSCurve(hEPOS_OFF, colors[i], 1, gLine, gBand);
            if (gLine)
                gLine->Draw("L same");
        }
        if (hEPOS_ON)
        {
            TGraph *gLine = nullptr;
            TGraphErrors *gBand = nullptr;
            makeEPOSCurve(hEPOS_ON, colors[i], 7, gLine, gBand);
            if (gLine)
                gLine->Draw("L same");
        }

        legPhi->AddEntry(hStat, Form("%s (#times 2^{%d})", multLabels[i].c_str(), (int)(multLabels.size() - 1 - i)), "pe");
    }

    legPhi->Draw();

    lat.DrawLatex(0.20, 0.90, "ALICE Run 3 pp #sqrt{#it{s}} = 13.6 TeV");
    lat.DrawLatex(0.20, 0.86, "#Phi(1020) Spectrum comparison");
    lat.DrawLatex(0.20, 0.82, "Solid: EPOS OFF | Dashed: EPOS ON");

    cPhi->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/SpectraCompareWithModel/Phi_Multiplicity_Spectra_EPOS.png");
}