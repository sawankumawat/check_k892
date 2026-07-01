#include <iostream>
#include <vector>
#include "../../src/style.h"

using namespace std;

void temp(){
    TFile *fPiKp = new TFile("ppRun3_PiKpHEP.root", "UPDATE"); // commented it for protection
    TFile *fPi = new TFile("../PiKp_Run3_Results/Sawan/Pi_results.root", "READ");
    TFile *fKa = new TFile("../PiKp_Run3_Results/Sawan/Ka_results.root", "READ");
    TFile *fPr = new TFile("../PiKp_Run3_Results/Sawan/Pr_results.root", "READ");

    if (fPiKp->IsZombie() || fPi->IsZombie() || fKa->IsZombie() || fPr->IsZombie())
    {
        cout << "Error: PiKp or Pi/Ka/Pr files not found" << endl;
        return;
    }

    TH1D *hPiSysUncorr = (TH1D *)fPi->Get("gMeanYieldRun3_sysuncorr");
    TH1D *hKaSysUncorr = (TH1D *)fKa->Get("gMeanYieldRun3_sysuncorr");
    TH1D *hPrSysUncorr = (TH1D *)fPr->Get("gMeanYieldRun3_sysuncorr");

    if (!hPiSysUncorr || !hKaSysUncorr || !hPrSysUncorr)
    {
        cout << "Error: SysUncorr histograms not found in Pi/Ka/Pr files" << endl;
        return;
    }

    fPiKp->cd();

    hPiSysUncorr->Write("gPion_MeanYield_sysuncorr", TObject::kOverwrite);
    hKaSysUncorr->Write("gKaon_MeanYield_sysuncorr", TObject::kOverwrite);
    hPrSysUncorr->Write("gProton_MeanYield_sysuncorr", TObject::kOverwrite);
}