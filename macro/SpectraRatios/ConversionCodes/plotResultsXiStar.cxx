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
    cout << "Total points in graphs: " << gYield->GetN() << endl;
    for (int i = 0; i < gYield->GetN(); i++)
    {
        double x, y;
        gYield->GetPoint(i, x, y);
        cout << "Point " << i << ": x = " << x << ", y = " << y << endl;
    }

    TGraphErrors *gYieldReverse = new TGraphErrors(gYield->GetN());
    TGraphErrors *gMeanpTReverse = new TGraphErrors(gMeanpT->GetN());
    TGraphErrors *gYieldSysReverse = new TGraphErrors(gYieldSys->GetN());
    TGraphErrors *gMeanpTSysReverse = new TGraphErrors(gMeanpTSys->GetN());
    for (int i = 0; i < gYield->GetN(); i++)
    {
        double x, y;
        gYield->GetPoint(i, x, y);
        double xErr = gYield->GetErrorX(i);
        double yErr = gYield->GetErrorY(i);
        gYieldReverse->SetPoint(gYield->GetN() - 1 - i, x, y);
        gYieldReverse->SetPointError(gYield->GetN() - 1 - i, xErr, yErr);

        double xsys, ysys;
        gYieldSys->GetPoint(i, xsys, ysys);
        double xSysErr = gYieldSys->GetErrorX(i);
        double ySysErr = gYieldSys->GetErrorY(i);
        gYieldSysReverse->SetPoint(gYieldSys->GetN() - 1 - i, xsys, ysys);
        gYieldSysReverse->SetPointError(gYieldSys->GetN() - 1 - i, xSysErr, ySysErr);

        double xpt, ypt;
        gMeanpT->GetPoint(i, xpt, ypt);
        double xPtErr = gMeanpT->GetErrorX(i);
        double yPtErr = gMeanpT->GetErrorY(i);
        gMeanpTReverse->SetPoint(gMeanpT->GetN() - 1 - i, xpt, ypt);
        gMeanpTReverse->SetPointError(gMeanpT->GetN() - 1 - i, xPtErr, yPtErr);

        double xptSys, yptSys;
        gMeanpTSys->GetPoint(i, xptSys, yptSys);
        double xPtSysErr = gMeanpTSys->GetErrorX(i);
        double yPtSysErr = gMeanpTSys->GetErrorY(i);
        gMeanpTSysReverse->SetPoint(gMeanpTSys->GetN() - 1 - i, xptSys, yptSys);
        gMeanpTSysReverse->SetPointError(gMeanpTSys->GetN() - 1 - i, xPtSysErr, yPtSysErr);
    }
    cout << endl;
    cout << "Reversed points in gYieldReverse: " << gYieldReverse->GetN() << endl;
    for (int i = 0; i < gYieldReverse->GetN(); i++)
    {
        double x, y;
        gYieldReverse->GetPoint(i, x, y);
        cout << "Point " << i << ": x = " << x << ", y = " << y << endl;
    }

    TFile *fOutput = new TFile("../XiStar_Results/Sawan/ResultsXiStar.root", "recreate");
    // gMeanpTSys->Write("gMeanpTRun3_sys");
    // gYieldSys->Write("gMeanYieldRun3_sys");
    // gMeanpT->Write("gMeanpTRun3");
    // gYield->Write("gMeanYieldRun3");
    gMeanpTSysReverse->Write("gMeanpTRun3_sys");
    gYieldSysReverse->Write("gMeanYieldRun3_sys");
    gMeanpTReverse->Write("gMeanpTRun3");
    gYieldReverse->Write("gMeanYieldRun3");
    fOutput->Close();
}