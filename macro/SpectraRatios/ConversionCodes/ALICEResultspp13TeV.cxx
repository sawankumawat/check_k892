#include <iostream>
#include <iomanip>
#include "../../src/style.h"

void ALICEResultspp13TeV()
{
    TFile *fALICE = new TFile("../HEP_data/pp13TeV_INELgt0.root", "read");
    if (fALICE->IsZombie())
    {
        cout << "Error: ALICE file not found" << endl;
        return;
    }

    // Phi has 10 multiplicity bins, while Kstar has 9 bins. So we need to read them separately
    int TableMeanpt = 40;    // for Phi
    int TableMeanYield = 42; // for Phi
    int TablePhi_KshortRatio = 44;
    int TablePhi_PiRatio = 46;
    int Table2Phi_KaRatio = 57;
    int TablePhi_PrRatio = 58;

    int TableKstar_KaRatio = 55;
    int TableKstar_PiRatio = 56;
    int TableKstar_KshortRatio = 43;

    TGraphErrors *gMeanYieldRun2 = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TableMeanYield));
    TH1D *hMeanYieldRun2_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TableMeanYield)); // stat error
    TH1D *hMeanYieldRun2_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TableMeanYield));  // sys error
    // TH1D *hMeanYieldRun2_sysUncorr = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e3", TableMeanYield)); // uncorrelated sys error

    TGraphErrors *gMeanpTRun2 = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TableMeanpt));
    TH1D *hMeanpTRun2_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TableMeanpt)); // stat error
    TH1D *hMeanpTRun2_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TableMeanpt));  // sys error
    // TH1D *hMeanpTRun2_sysUncorr = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e3", TableMeanpt)); // uncorrelated sys error

    TGraphErrors *gPhi_KshortRatio = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TablePhi_KshortRatio));
    TH1D *hPhi_KshortRatio_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TablePhi_KshortRatio)); // stat error
    TH1D *hPhi_KshortRatio_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TablePhi_KshortRatio));  // sys error

    TGraphErrors *gPhi_PiRatio = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TablePhi_PiRatio));
    TH1D *hPhi_PiRatio_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TablePhi_PiRatio)); // stat error
    TH1D *hPhi_PiRatio_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TablePhi_PiRatio));  // sys error

    TGraphErrors *gPhi_KaRatio = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", Table2Phi_KaRatio));
    TH1D *hPhi_KaRatio_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", Table2Phi_KaRatio)); // stat error
    TH1D *hPhi_KaRatio_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", Table2Phi_KaRatio));  // sys error

    TGraphErrors *gPhi_PrRatio = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TablePhi_PrRatio));
    TH1D *hPhi_PrRatio_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TablePhi_PrRatio)); // stat error
    TH1D *hPhi_PrRatio_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TablePhi_PrRatio));  // sys error

    TGraphErrors *gMeanYieldRun2_stat = new TGraphErrors(gMeanYieldRun2->GetN());
    TGraphErrors *gMeanYieldRun2_sys = new TGraphErrors(gMeanYieldRun2->GetN());

    TGraphErrors *gMeanpTRun2_stat = new TGraphErrors(gMeanpTRun2->GetN());
    TGraphErrors *gMeanpTRun2_sys = new TGraphErrors(gMeanpTRun2->GetN());

    TGraphErrors *gPhi_KshortRatio_stat = new TGraphErrors(gPhi_KshortRatio->GetN());
    TGraphErrors *gPhi_KshortRatio_sys = new TGraphErrors(gPhi_KshortRatio->GetN());

    TGraphErrors *gPhi_PiRatio_stat = new TGraphErrors(gPhi_PiRatio->GetN());
    TGraphErrors *gPhi_PiRatio_sys = new TGraphErrors(gPhi_PiRatio->GetN());

    TGraphErrors *gPhi_KaRatio_stat = new TGraphErrors(gPhi_KaRatio->GetN());
    TGraphErrors *gPhi_KaRatio_sys = new TGraphErrors(gPhi_KaRatio->GetN());

    TGraphErrors *gPhi_PrRatio_stat = new TGraphErrors(gPhi_PrRatio->GetN());
    TGraphErrors *gPhi_PrRatio_sys = new TGraphErrors(gPhi_PrRatio->GetN());

    for (int i = 0; i < gMeanYieldRun2->GetN(); ++i)
    {
        double x = 0.0, y = 0.0;
        gMeanYieldRun2->GetPoint(i, x, y);
        int bin = hMeanYieldRun2_stat->FindBin(x);
        double ex = gMeanYieldRun2->GetErrorX(i);
        double ey_stat = hMeanYieldRun2_stat->GetBinContent(bin);
        double ey_sys = hMeanYieldRun2_sys->GetBinContent(bin);
        // double ey_sysUncorr = hMeanYieldRun2_sysUncorr->GetBinContent(bin);

        gMeanYieldRun2_stat->SetPoint(i, x, y);
        gMeanYieldRun2_stat->SetPointError(i, ex, ey_stat);

        gMeanYieldRun2_sys->SetPoint(i, x, y);
        gMeanYieldRun2_sys->SetPointError(i, ex, ey_sys);

        // mean pT
        gMeanpTRun2->GetPoint(i, x, y);
        double pt_stat = hMeanpTRun2_stat->GetBinContent(bin);
        double pt_sys = hMeanpTRun2_sys->GetBinContent(bin);
        // double pt_sysUncorr = hMeanpTRun2_sysUncorr->GetBinContent(bin);
        gMeanpTRun2_stat->SetPoint(i, x, y);
        gMeanpTRun2_stat->SetPointError(i, ex, pt_stat);

        gMeanpTRun2_sys->SetPoint(i, x, y);
        gMeanpTRun2_sys->SetPointError(i, ex, pt_sys);

        // Phi/Kshort ratio
        gPhi_KshortRatio->GetPoint(i, x, y);
        double phiKshort_stat = hPhi_KshortRatio_stat->GetBinContent(bin);
        double phiKshort_sys = hPhi_KshortRatio_sys->GetBinContent(bin);
        gPhi_KshortRatio_stat->SetPoint(i, x, y);
        gPhi_KshortRatio_stat->SetPointError(i, ex, phiKshort_stat);
        gPhi_KshortRatio_sys->SetPoint(i, x, y);
        gPhi_KshortRatio_sys->SetPointError(i, ex, phiKshort_sys);

        // Phi/Pi ratio
        gPhi_PiRatio->GetPoint(i, x, y);
        double phiPi_stat = hPhi_PiRatio_stat->GetBinContent(bin);
        double phiPi_sys = hPhi_PiRatio_sys->GetBinContent(bin);
        gPhi_PiRatio_stat->SetPoint(i, x, y);
        gPhi_PiRatio_stat->SetPointError(i, ex, phiPi_stat);
        gPhi_PiRatio_sys->SetPoint(i, x, y);
        gPhi_PiRatio_sys->SetPointError(i, ex, phiPi_sys);

        // 2*Phi/Ka ratio
        gPhi_KaRatio->GetPoint(i, x, y);
        double phiKa_stat = hPhi_KaRatio_stat->GetBinContent(bin);
        double phiKa_sys = hPhi_KaRatio_sys->GetBinContent(bin);
        gPhi_KaRatio_stat->SetPoint(i, x, y);
        gPhi_KaRatio_stat->SetPointError(i, ex, phiKa_stat);
        gPhi_KaRatio_sys->SetPoint(i, x, y);
        gPhi_KaRatio_sys->SetPointError(i, ex, phiKa_sys);

        // Phi/Pr ratio
        gPhi_PrRatio->GetPoint(i, x, y);
        double phiPr_stat = hPhi_PrRatio_stat->GetBinContent(bin);
        double phiPr_sys = hPhi_PrRatio_sys->GetBinContent(bin);
        gPhi_PrRatio_stat->SetPoint(i, x, y);
        gPhi_PrRatio_stat->SetPointError(i, ex, phiPr_stat);
        gPhi_PrRatio_sys->SetPoint(i, x, y);
        gPhi_PrRatio_sys->SetPointError(i, ex, phiPr_sys);
    }

    TGraphErrors *gKstar_KaRatio = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TableKstar_KaRatio));
    TH1D *hKstar_KaRatio_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TableKstar_KaRatio)); // stat error
    TH1D *hKstar_KaRatio_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TableKstar_KaRatio));  // sys error

    TGraphErrors *gKstar_PiRatio = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TableKstar_PiRatio));
    TH1D *hKstar_PiRatio_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TableKstar_PiRatio)); // stat error
    TH1D *hKstar_PiRatio_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TableKstar_PiRatio));  // sys error

    TGraphErrors *gKstar_KshortRatio = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TableKstar_KshortRatio));
    TH1D *hKstar_KshortRatio_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TableKstar_KshortRatio)); // stat error
    TH1D *hKstar_KshortRatio_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TableKstar_KshortRatio));  // sys error

    TGraphErrors *gKstar_KaRatio_stat = new TGraphErrors(gKstar_KaRatio->GetN());
    TGraphErrors *gKstar_KaRatio_sys = new TGraphErrors(gKstar_KaRatio->GetN());

    TGraphErrors *gKstar_PiRatio_stat = new TGraphErrors(gKstar_PiRatio->GetN());
    TGraphErrors *gKstar_PiRatio_sys = new TGraphErrors(gKstar_PiRatio->GetN());

    TGraphErrors *gKstar_KshortRatio_stat = new TGraphErrors(gKstar_KshortRatio->GetN());
    TGraphErrors *gKstar_KshortRatio_sys = new TGraphErrors(gKstar_KshortRatio->GetN());

    for (int i = 0; i < gKstar_KaRatio->GetN(); ++i)
    {
        double x = 0.0, y = 0.0;
        gKstar_KaRatio->GetPoint(i, x, y);
        int bin = hKstar_KaRatio_stat->FindBin(x);
        double ex = gKstar_KaRatio->GetErrorX(i);
        double ey_stat = hKstar_KaRatio_stat->GetBinContent(bin);
        double ey_sys = hKstar_KaRatio_sys->GetBinContent(bin);

        gKstar_KaRatio_stat->SetPoint(i, x, y);
        gKstar_KaRatio_stat->SetPointError(i, ex, ey_stat);

        gKstar_KaRatio_sys->SetPoint(i, x, y);
        gKstar_KaRatio_sys->SetPointError(i, ex, ey_sys);

        // Kstar/Pi ratio
        gKstar_PiRatio->GetPoint(i, x, y);
        double kstarPi_stat = hKstar_PiRatio_stat->GetBinContent(bin);
        double kstarPi_sys = hKstar_PiRatio_sys->GetBinContent(bin);
        gKstar_PiRatio_stat->SetPoint(i, x, y);
        gKstar_PiRatio_stat->SetPointError(i, ex, kstarPi_stat);
        gKstar_PiRatio_sys->SetPoint(i, x, y);
        gKstar_PiRatio_sys->SetPointError(i, ex, kstarPi_sys);

        // Kstar/Kshort ratio
        gKstar_KshortRatio->GetPoint(i, x, y);
        double kstarKshort_stat = hKstar_KshortRatio_stat->GetBinContent(bin);
        double kstarKshort_sys = hKstar_KshortRatio_sys->GetBinContent(bin);
        gKstar_KshortRatio_stat->SetPoint(i, x, y);
        gKstar_KshortRatio_stat->SetPointError(i, ex, kstarKshort_stat);
        gKstar_KshortRatio_sys->SetPoint(i, x, y);
        gKstar_KshortRatio_sys->SetPointError(i, ex, kstarKshort_sys);
    }

    TFile *fOutput = new TFile("pp13TeVALICE.root", "recreate");
    gMeanYieldRun2_stat->Write("gMeanYieldRun2_stat");
    gMeanYieldRun2_sys->Write("gMeanYieldRun2_sys");

    gMeanpTRun2_stat->Write("gMeanpTRun2_stat");
    gMeanpTRun2_sys->Write("gMeanpTRun2_sys");

    gPhi_KshortRatio_stat->Write("gPhi_KshortRatio_stat");
    gPhi_KshortRatio_sys->Write("gPhi_KshortRatio_sys");

    gPhi_PiRatio_stat->Write("gPhi_PiRatio_stat");
    gPhi_PiRatio_sys->Write("gPhi_PiRatio_sys");

    gPhi_KaRatio_stat->Write("gPhi_KaRatio_stat");
    gPhi_KaRatio_sys->Write("gPhi_KaRatio_sys");

    gPhi_PrRatio_stat->Write("gPhi_PrRatio_stat");
    gPhi_PrRatio_sys->Write("gPhi_PrRatio_sys");

    gKstar_KaRatio_stat->Write("gKstar_KaRatio_stat");
    gKstar_KaRatio_sys->Write("gKstar_KaRatio_sys");

    gKstar_PiRatio_stat->Write("gKstar_PiRatio_stat");
    gKstar_PiRatio_sys->Write("gKstar_PiRatio_sys");

    gKstar_KshortRatio_stat->Write("gKstar_KshortRatio_stat");
    gKstar_KshortRatio_sys->Write("gKstar_KshortRatio_sys");

    fOutput->Close();
}