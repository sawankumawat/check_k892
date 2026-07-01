#include <iostream>
#include <iomanip>
#include "../../src/style.h"

using namespace std;

void LambdaStarRun2()
{
    std::vector<std::vector<double>> data = {
        {18.567, 0.28, 0.0499482, 0.000398665, 0.002845, 0.0023779, 1.67, 0.00868841, 0.0479393, 0.0453262},
        {11.51, 0.17, 0.0289709, 0.000202132, 0.00193851, 0.00170737, 1.48091, 0.00640869, 0.0483154, 0.0464396},
        {7.275, 0.11, 0.0185552, 0.000156064, 0.00141009, 0.00126209, 1.28134, 0.00620787, 0.0440264, 0.0419529},
        {4.64, 0.07, 0.0113894, 0.000110996, 0.000879422, 0.000762244, 1.12545, 0.00610979, 0.0405509, 0.0386648},
        {2.55, 0.04, 0.00532433, 6.43593e-05, 0.000526658, 0.000487455, 0.978536, 0.00605513, 0.0421375, 0.0409937},
        {6.94, 0.10, 0.0179722, 9.0899e-05, 0.00103682, 0.000836964, 1.4128, 0.00425586, 0.0461922, 0.0435410}};

    // dN/dy
    TGraphErrors *gYieldStat = new TGraphErrors();
    TGraphErrors *gYieldSys = new TGraphErrors();
    TGraphErrors *gYieldSysUncorr = new TGraphErrors();

    // <pT>
    TGraphErrors *gMeanPtStat = new TGraphErrors();
    TGraphErrors *gMeanPtSys = new TGraphErrors();
    TGraphErrors *gMeanPtSysUncorr = new TGraphErrors();

    for (size_t i = 0; i < data.size() - 1; ++i)
    {
        double mult = data[i][0];
        double multErr = data[i][1];

        double dNdy = data[i][2];
        double dNdyStat = data[i][3];
        double dNdySys = data[i][4];
        double dNdySysUnc = data[i][5];

        double meanPt = data[i][6];
        double meanPtStat = data[i][7];
        double meanPtSys = data[i][8];
        double meanPtSysUnc = data[i][9];

        // dN/dy
        gYieldStat->SetPoint(i, mult, dNdy);
        gYieldStat->SetPointError(i, multErr, dNdyStat);

        gYieldSys->SetPoint(i, mult, dNdy);
        gYieldSys->SetPointError(i, multErr, dNdySys);

        gYieldSysUncorr->SetPoint(i, mult, dNdy);
        gYieldSysUncorr->SetPointError(i, multErr, dNdySysUnc);

        // <pT>
        gMeanPtStat->SetPoint(i, mult, meanPt);
        gMeanPtStat->SetPointError(i, multErr, meanPtStat);

        gMeanPtSys->SetPoint(i, mult, meanPt);
        gMeanPtSys->SetPointError(i, multErr, meanPtSys);

        gMeanPtSysUncorr->SetPoint(i, mult, meanPt);
        gMeanPtSysUncorr->SetPointError(i, multErr, meanPtSysUnc);
    }

    TFile fout("pp13TeV_LambdaStar.root", "RECREATE");

    gYieldStat->Write("gLambda_MeanYield_stat");
    gYieldSys->Write("gLambda_MeanYield_sys");
    gYieldSysUncorr->Write("gLambda_MeanYield_sysuncorr");

    gMeanPtStat->Write("gLambda_MeanpT_stat");
    gMeanPtSys->Write("gLambda_MeanpT_sys");
    gMeanPtSysUncorr->Write("gLambda_MeanpT_sysuncorr");

    // LambdaStar/Lambda ratio (from XY scan) (with systematic errors) (stat errors give by hand temporarily)
    //  x,     dx,       y,         -dy,    +dy
    // 18.567  0.28	   0.9407	0.07692	0.07692
    // 11.51   0.17    0.8813	0.08132	0.08132
    // 7.275   0.11    0.9231	0.09011	0.09011
    // 4.64	   0.07    0.9473	0.09451	0.09451
    // 2.55    0.04     1.0	     0.1341	0.1341

    // Lambda*/Lambda ratio

    double mult[] = {18.567, 11.51, 7.275, 4.64, 2.55};
    double multErr[] = {0.28, 0.17, 0.11, 0.07, 0.04};

    double ratio[] = {0.9407, 0.8813, 0.9231, 0.9473, 1.0000};

    // Systematic uncertainties
    double sysErr[] = {0.07692, 0.08132, 0.09011, 0.09451, 0.13410};

    // Statistical uncertainties (2%)
    double statErr[5];
    for (int i = 0; i < 5; ++i)
        statErr[i] = 0.02 * ratio[i];

    // Create graphs
    TGraphErrors *gLambdaStarLambdaStat = new TGraphErrors(5);
    TGraphErrors *gLambdaStarLambdaSys = new TGraphErrors(5);

    for (int i = 0; i < 5; ++i)
    {
        gLambdaStarLambdaStat->SetPoint(i, mult[i], ratio[i]);
        gLambdaStarLambdaStat->SetPointError(i, multErr[i], statErr[i]);

        gLambdaStarLambdaSys->SetPoint(i, mult[i], ratio[i]);
        gLambdaStarLambdaSys->SetPointError(i, multErr[i], sysErr[i]);
    }

    // Optional names
    gLambdaStarLambdaStat->SetName("gLambdaStar_LambdaRatio_stat");
    gLambdaStarLambdaSys->SetName("gLambdaStar_LambdaRatio_sys");

    gLambdaStarLambdaStat->Write();
    gLambdaStarLambdaSys->Write();

    fout.Close();
}