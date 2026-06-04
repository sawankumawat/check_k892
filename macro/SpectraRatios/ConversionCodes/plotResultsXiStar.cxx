#include <iostream>
#include <iomanip>
#include "../../src/style.h"

void plotResultsXiStar()
{
    TFile *fYield = new TFile("../XiStar_Results/XiStar_pp13p6_Draw02_Yield_Oct02.root", "read");
    TFile *fPt = new TFile("../XiStar_Results/XiStar_pp13p6_Draw02_pT_Oct02.root", "read");
    if (fYield->IsZombie() || fPt->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }
    TGraphErrors *gYield = (TGraphErrors *)fYield->Get("gdNdy_Stat");
    TGraphErrors *gMeanpT = (TGraphErrors *)fPt->Get("gpT_Stat");
    TGraphErrors *gYieldSys = (TGraphErrors *)fYield->Get("gdNdy_TotalSys");
    TGraphErrors *gMeanpTSys = (TGraphErrors *)fPt->Get("gpT_TotalSys");
    if (gYield == nullptr || gMeanpT == nullptr || gYieldSys == nullptr || gMeanpTSys == nullptr)
    {
        cout << "Error: histograms not found" << endl;
        return;
    }
    gYield->RemovePoint(0);
    gYield->RemovePoint(0);
    gMeanpT->RemovePoint(0);
    gMeanpT->RemovePoint(0);
    gYieldSys->RemovePoint(0);
    gYieldSys->RemovePoint(0);
    gMeanpTSys->RemovePoint(0);
    gMeanpTSys->RemovePoint(0);
    cout << "Total points in graphs: " << gMeanpT->GetN() << endl;
    for (int i = 0; i < gMeanpT->GetN(); i++)
    {
        double x, y;
        gMeanpT->GetPoint(i, x, y);
        cout << "Point " << i << ": x = " << x << ", y = " << y << endl;
    }

    TFile *fOutput = new TFile("../XiStar_Results/Sawan/ResultsXiStar.root", "recreate");
    gMeanpTSys->Write("gMeanpTRun3_sys");
    gYieldSys->Write("gMeanYieldRun3_sys");
    gMeanpT->Write("gMeanpTRun3");
    gYield->Write("gMeanYieldRun3");
    fOutput->Close();
}