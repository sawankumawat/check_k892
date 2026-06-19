#include <iostream>
#include <iomanip>
#include "../../src/style.h"

void plotResultsRho()
{
    TFile *fResults = new TFile("../Rho_Run3_Results/DrawdNdy_0315_13p6TeV_TGraphErrors.root", "read");
    if (fResults->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }
    TGraphErrors *hYieldStat = (TGraphErrors *)fResults->Get("g13TeV_dNdy_stat");
    TGraphErrors *hMeanpTStat = (TGraphErrors *)fResults->Get("g13TeV_meanpT_stat");
    TGraphErrors *hYieldSys = (TGraphErrors *)fResults->Get("g13TeV_dNdy_sys");
    TGraphErrors *hMeanpTSys = (TGraphErrors *)fResults->Get("g13TeV_meanpT_sys");

    if (hYieldStat == nullptr || hMeanpTStat == nullptr || hYieldSys == nullptr || hMeanpTSys == nullptr)
    {
        cout << "Error: graphs not found" << endl;
        return;
    }

    cout << "Total points in graphs: " << hMeanpTStat->GetN() << endl;
    for (int i = 0; i < hMeanpTStat->GetN(); i++)
    {
        double x, y;
        hMeanpTStat->GetPoint(i, x, y);
        cout << "Point " << i << ": x = " << x << ", y = " << y << endl;
    }

    TFile *fOutput = new TFile("../Rho_Run3_Results/Sawan/ResultsRho.root", "recreate");
    hMeanpTStat->Write("gMeanpTRun3_stat");
    hYieldStat->Write("gMeanYieldRun3_stat");
    hMeanpTSys->Write("gMeanpTRun3_sys");
    hYieldSys->Write("gMeanYieldRun3_sys");
    fOutput->Close();
}