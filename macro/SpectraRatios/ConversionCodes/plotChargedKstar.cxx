#include <iostream>
#include <iomanip>
#include "../../src/style.h"

void plotChargedKstar()
{
    TFile *fSpectra = new TFile("../ChargedKstar_Results/outfile_hist_dNdy_default.root", "read");
    if (fSpectra->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }
    TH1D *hSpectra = (TH1D *)fSpectra->Get("h1d_pTspectra_MB_file00");
    if (hSpectra == nullptr)
    {
        cout << "Error: histogram not found" << endl;
        return;
    }
    TGraphErrors *gYield = (TGraphErrors *)fSpectra->Get("gdNdyStat");
    TGraphErrors *gMeanpT = (TGraphErrors *)fSpectra->Get("gMpTStat");
    TGraphErrors *gYieldSys = (TGraphErrors *)fSpectra->Get("gdNdySyst");
    TGraphErrors *gMeanpTSys = (TGraphErrors *)fSpectra->Get("gMpTSyst");
    if (gYield == nullptr || gMeanpT == nullptr || gYieldSys == nullptr || gMeanpTSys == nullptr)
    {
        cout << "Error: histograms not found" << endl;
        return;
    }

    cout << "Total points in graphs: " << gMeanpT->GetN() << endl;
    for (int i = 0; i < gMeanpT->GetN(); i++)
    {
        double x, y;
        gMeanpT->GetPoint(i, x, y);
        cout << "Point " << i + 1 << ": x = " << x << ", y = " << y << endl;
    }

    TFile *fOutput = new TFile("../ChargedKstar_Results/Sawan/ResultsChargedKstar.root", "recreate");
    gMeanpTSys->Write("gMeanpTRun3_sys");
    gYieldSys->Write("gMeanYieldRun3_sys");
    gMeanpT->Write("gMeanpTRun3");
    gYield->Write("gMeanYieldRun3");
    fOutput->Close();
}