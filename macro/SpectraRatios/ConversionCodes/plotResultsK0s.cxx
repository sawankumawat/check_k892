#include <iostream>
#include <iomanip>
#include "../../src/style.h"

void plotResultsK0s()
{
    TFile *fSpectra = new TFile("../K0s_Run3_Results/Results_K0Short_LHC24apass1_TVX_TF_ITSnoMC_MBEff_FT0M0100.root", "read");
    TFile *fResults = new TFile("../K0s_Run3_Results/YieldsIntegratedK0Short_LHC24.root", "read");
    if (fSpectra->IsZombie() || fResults->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }
    TH1F *hCorrectedMinBiasSpectra = (TH1F *)fSpectra->Get("hPtCorrected");
    TGraphErrors *hMeanpT = (TGraphErrors *)fResults->Get("MeanptNchStat");
    TGraphErrors *hMeanpTsys = (TGraphErrors *)fResults->Get("MeanptNchSys");
    TGraphErrors *hMeanpTsysUncorr = (TGraphErrors *)fResults->Get("MeanptNchSysUncorr");
    TGraphErrors *hYield = (TGraphErrors *)fResults->Get("YieldsNchStat");
    TGraphErrors *hYieldsys = (TGraphErrors *)fResults->Get("YieldsNchSys");
    TGraphErrors *hYieldsysUncorr = (TGraphErrors *)fResults->Get("YieldsNchSysUncorr");
    if (hCorrectedMinBiasSpectra == nullptr || hMeanpT == nullptr || hYield == nullptr || hMeanpTsys == nullptr || hYieldsys == nullptr || hMeanpTsysUncorr == nullptr || hYieldsysUncorr == nullptr)
    {
        cout << "Error: histograms not found" << endl;
        return;
    }

    cout << "Total points in graphs: " << hMeanpT->GetN() << endl;
    for (int i = 0; i < hMeanpT->GetN(); i++)
    {
        double x, y;
        hMeanpT->GetPoint(i, x, y);
        cout << "Point " << i << ": x = " << x << ", y = " << y << endl;
    }

    TFile *fOutput = new TFile("../K0s_Run3_Results/Sawan/ResultsK0s.root", "recreate");
    hCorrectedMinBiasSpectra->Write("hCorrectedMinBiasSpectra");
    hMeanpT->Write("gMeanpTRun3_stat");
    hMeanpTsys->Write("gMeanpTRun3_sys");
    hMeanpTsysUncorr->Write("gMeanpTRun3_sysuncorr");
    hYield->Write("gMeanYieldRun3_stat");
    hYieldsys->Write("gMeanYieldRun3_sys");
    hYieldsysUncorr->Write("gMeanYieldRun3_sysuncorr");
    fOutput->Close();
}