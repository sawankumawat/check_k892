#include <iostream>
#include <iomanip>
#include "../../src/style.h"

void plotResultsLambdaStar()
{
    TFile *fResuls = new TFile("../LambdaRun3/Lambda1520_Yield_MeanPt_13_6TeV.root", "read");
    if (fResuls->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }
    TGraphErrors *gYield = (TGraphErrors *)fResuls->Get("Yield_Stat");
    TGraphErrors *gMeanPt = (TGraphErrors *)fResuls->Get("MeanPt_Stat");
    TGraphErrors *gYieldSys = (TGraphErrors *)fResuls->Get("Yield_SysTot");
    TGraphErrors *gMeanPtSys = (TGraphErrors *)fResuls->Get("MeanPt_SysTot");
    if (gYield == nullptr || gMeanPt == nullptr || gYieldSys == nullptr || gMeanPtSys == nullptr)
    {
        cout << "Error: graphs not found" << endl;
        return;
    }

    cout << "Total points in graphs: " << gYield->GetN() << endl;
    for (int i = 0; i < gYield->GetN(); i++)
    {
        double x, y;
        gYield->GetPoint(i, x, y);
        cout << "Point " << i << ": x = " << x << ", y = " << y << endl;
    }

    TFile *fOutput = new TFile("../LambdaRun3/Sawan/ResultsLambda1520.root", "recreate");
    gYield->Write("gMeanYieldRun3_stat");
    gMeanPt->Write("gMeanpTRun3_stat");
    gYieldSys->Write("gMeanYieldRun3_sys");
    gMeanPtSys->Write("gMeanpTRun3_sys");
    fOutput->Close();
}