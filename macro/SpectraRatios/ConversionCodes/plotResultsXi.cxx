#include <iostream>
#include <iomanip>
#include "../../src/style.h"

void plotResultsXi()
{
    TFile *fResults = new TFile("../XiOmegaResults/Xi.root", "read");
    if (fResults->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }

    TH1F *hYieldStat = (TH1F *)fResults->Get("hIntegratedYield");
    TH1F *hYieldSys = (TH1F *)fResults->Get("hIntegratedYieldSyst");
    TH1F *hYieldSysUncorr = (TH1F *)fResults->Get("hIntegratedYieldSystUncorr");

    // First bin is 0-100% and then in decreasing order i.e highest multiplicity to lowest multiplicity
    TGraphErrors *gYieldStat = new TGraphErrors(10);
    TGraphErrors *gYieldSys = new TGraphErrors(10);
    TGraphErrors *gYieldSysUncorr = new TGraphErrors(10);

    double dnch_detaRun3[] = {21.212, 18.069, 15.483, 13.673, 12.317, 10.707, 8.964, 7.496, 5.699, 3.296}; // updated uncertainties (matching correction)
    double dnch_detaRun3_err[] = {0.38, 0.25, 0.22, 0.19, 0.17, 0.15, 0.13, 0.11, 0.09, 0.06};             // using same as previous one

    for (int i = 1; i < 11; i++)
    {
        double x = dnch_detaRun3[i - 1];
        double ex = dnch_detaRun3_err[i - 1];
        double y = hYieldStat->GetBinContent(i + 1);
        double ey_stat = hYieldStat->GetBinError(i + 1);
        double ey_sys = hYieldSys->GetBinError(i + 1);
        double ey_sysUncorr = hYieldSysUncorr->GetBinError(i + 1);

        gYieldStat->SetPoint(i - 1, x, y);
        gYieldStat->SetPointError(i - 1, ex, ey_stat);

        gYieldSys->SetPoint(i - 1, x, y);
        gYieldSys->SetPointError(i - 1, ex, ey_sys);

        gYieldSysUncorr->SetPoint(i - 1, x, y);
        gYieldSysUncorr->SetPointError(i - 1, ex, ey_sysUncorr);
    }

    TFile *fOutput = new TFile("../XiOmegaResults/Sawan/ResultsXi.root", "recreate");
    gYieldStat->Write("gMeanYieldRun3_stat");
    gYieldSys->Write("gMeanYieldRun3_sys");
    gYieldSysUncorr->Write("gMeanYieldRun3_sysuncorr");
}