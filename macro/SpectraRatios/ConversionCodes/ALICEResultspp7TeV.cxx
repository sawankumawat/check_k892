#include <iostream>
#include <iomanip>
#include "../../src/style.h"

void ALICEResultspp7TeV()
{
    TFile *fALICE = new TFile("../HEP_data/HEPData_7TeV_KstarPhi_INELgt0.root", "read");
    if (fALICE->IsZombie())
    {
        cout << "Error: ALICE file not found" << endl;
        return;
    }

    // Phi has 10 multiplicity bins, while Kstar has 9 bins. So we need to read them separately
    int TableKstar_PiRatio = 96;
    int TablePhi_PiRatio = 97;

    TGraphErrors *gKstar_PiRatio = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TableKstar_PiRatio));
    TH1D *hKstarPiRatio_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TableKstar_PiRatio)); // stat error
    TH1D *hKstarPiRatio_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TableKstar_PiRatio));  // sys error

    TGraphErrors *gPhi_PiRatio = (TGraphErrors *)fALICE->Get(Form("Table %d/Graph1D_y1", TablePhi_PiRatio));
    TH1D *hPhi_PiRatio_stat = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e1", TablePhi_PiRatio)); // stat error
    TH1D *hPhi_PiRatio_sys = (TH1D *)fALICE->Get(Form("Table %d/Hist1D_y1_e2", TablePhi_PiRatio));  // sys error

    TGraphErrors *gKstar_PiRatio_stat = new TGraphErrors(gKstar_PiRatio->GetN());
    TGraphErrors *gKstar_PiRatio_sys = new TGraphErrors(gKstar_PiRatio->GetN());

    TGraphErrors *gPhi_PiRatio_stat = new TGraphErrors(gPhi_PiRatio->GetN());
    TGraphErrors *gPhi_PiRatio_sys = new TGraphErrors(gPhi_PiRatio->GetN());

    for (int i = 0; i < gKstar_PiRatio->GetN(); ++i)
    {
        double x = 0.0, y = 0.0;
        gKstar_PiRatio->GetPoint(i, x, y);
        int bin = hKstarPiRatio_stat->FindBin(x);
        double ex = gKstar_PiRatio->GetErrorX(i);
        double ey_stat = hKstarPiRatio_stat->GetBinContent(bin);
        double ey_sys = hKstarPiRatio_sys->GetBinContent(bin);
        // double ey_sysUncorr = hKstarPiRatio_sysUncorr->GetBinContent(bin);

        gKstar_PiRatio_stat->SetPoint(i, x, y);
        gKstar_PiRatio_stat->SetPointError(i, ex, ey_stat);

        gKstar_PiRatio_sys->SetPoint(i, x, y);
        gKstar_PiRatio_sys->SetPointError(i, ex, ey_sys);

        // Phi/Pi ratio
        gPhi_PiRatio->GetPoint(i, x, y);
        double phiPi_stat = hPhi_PiRatio_stat->GetBinContent(bin);
        double phiPi_sys = hPhi_PiRatio_sys->GetBinContent(bin);
        gPhi_PiRatio_stat->SetPoint(i, x, y);
        gPhi_PiRatio_stat->SetPointError(i, ex, phiPi_stat);
        gPhi_PiRatio_sys->SetPoint(i, x, y);
        gPhi_PiRatio_sys->SetPointError(i, ex, phiPi_sys);
    }

    TFile *fOutput = new TFile("pp7TeVALICE.root", "recreate");
    gKstar_PiRatio_stat->Write("gKstar_PiRatio_stat");
    gKstar_PiRatio_sys->Write("gKstar_PiRatio_sys");

    gPhi_PiRatio_stat->Write("gPhi_PiRatio_stat");
    gPhi_PiRatio_sys->Write("gPhi_PiRatio_sys");

    fOutput->Close();
}