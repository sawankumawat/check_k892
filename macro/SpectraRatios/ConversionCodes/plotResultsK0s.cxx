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
    TGraphErrors *hYield = (TGraphErrors *)fResults->Get("YieldsNchStat");
    if (hCorrectedMinBiasSpectra == nullptr || hMeanpT == nullptr || hYield == nullptr)
    {
        cout << "Error: histograms not found" << endl;
        return;
    }

    TFile *fOutput = new TFile("../K0s_Run3_Results/Sawan/ResultsK0s.root", "recreate");
    hCorrectedMinBiasSpectra->Write("hCorrectedMinBiasSpectra");
    hMeanpT->Write("gMeanpTRun3");
    hYield->Write("gMeanYieldRun3");
    hMeanpT->Write("gMeanpTRun3_sys");
    hYield->Write("gMeanYieldRun3_sys");
    fOutput->Close();
}