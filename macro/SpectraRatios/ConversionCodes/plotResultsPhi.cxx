#include <iostream>
#include <vector>
#include "../../src/style.h"

using namespace std;

void plotResultsPhi()
{
    TFile *fResults = new TFile("../PhiResults/PhiAnalysisResults_pprun3.root", "read");
    if (fResults->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }

    TGraphErrors *hYieldStat = (TGraphErrors *)fResults->Get("YieldMean/dNdy/gYieldStat");
    TGraphErrors *hMeanpTStat = (TGraphErrors *)fResults->Get("YieldMean/MeanPt/gMeanStat");
    TGraphErrors *hYieldSys = (TGraphErrors *)fResults->Get("YieldMean/dNdy/gYieldSys");
    TGraphErrors *hMeanpTSys = (TGraphErrors *)fResults->Get("YieldMean/MeanPt/gMeanSys");

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

    TFile *fOutput = new TFile("../PhiResults/Sawan/ResultsPhi.root", "recreate");
    hMeanpTStat->Write("gMeanpTRun3_stat");
    hYieldStat->Write("gMeanYieldRun3_stat");
    hMeanpTSys->Write("gMeanpTRun3_sys");
    hYieldSys->Write("gMeanYieldRun3_sys");
    fOutput->Close();
}