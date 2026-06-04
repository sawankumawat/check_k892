#include <iostream>
#include <iomanip>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/initializations.h"

int colors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

void compareRatio()
{
    // TFile *fData1 = new TFile("levy_fit.root", "read");
    TFile *fData1 = new TFile("PiKp_Run3_Results/Sawan/Pr_results.root", "read");
    TFile *fData2 = new TFile("PiKp_Run3_Results/Sawan/Pi_results.root", "read");
    TFile *fModel = new TFile("ModelResults.root", "read");
    if (fData1->IsZombie() || fModel->IsZombie() || fData2->IsZombie())
    {
        cout << "Error: Data or Model file not found" << endl;
        return;
    }

    // const vector<string> modelNames = {"EPOS", "Pythia_Monash", "Pythia_CR", "Pythia_Ropes", "Pythia_Shoving"};
    const vector<string> modelNames = {"EPOS", "Pythia_Monash", "Pythia_Ropes", "Pythia_Shoving"};
    TGraphErrors *gMeanYieldModels1[modelNames.size()];
    TGraphErrors *gMeanYieldModels2[modelNames.size()];
    TGraphErrors *gMeanYieldModelsRatio[modelNames.size()];
    for (size_t i = 0; i < modelNames.size(); ++i)
    {
        TDirectory *modelDir = (TDirectory *)fModel->Get(modelNames[i].c_str());
        if (modelDir == nullptr)
        {
            cout << "Error: Model directory not found for " << modelNames[i] << endl;
            return;
        }
        gMeanYieldModels1[i] = (TGraphErrors *)modelDir->Get("Proton/gMeanYield_Proton");
        gMeanYieldModels2[i] = (TGraphErrors *)modelDir->Get("Pion/gMeanYield_Pion");
        if (gMeanYieldModels1[i] == nullptr || gMeanYieldModels2[i] == nullptr)
        {
            cout << "Error: Model graphs not found for " << modelNames[i] << endl;
            return;
        }
        gMeanYieldModelsRatio[i] = (TGraphErrors *)gMeanYieldModels1[i]->Clone(Form("%s_Ratio", modelNames[i].c_str()));
        gMeanYieldModelsRatio[i]->SetTitle(Form("%s Proton/Pion Ratio", modelNames[i].c_str()));
        for (int j = 0; j < gMeanYieldModels1[i]->GetN(); ++j)
        {
            double x1, y1, x2, y2;
            gMeanYieldModels1[i]->GetPoint(j, x1, y1);
            gMeanYieldModels2[i]->GetPoint(j, x2, y2);
            if (fabs(x1 - x2) > 1e-5)
            {
                cout << "Error: Mismatched x values for model " << modelNames[i] << " at point " << j << ": " << x1 << " vs " << x2 << endl;
                return;
            }
            double ratio = (y2 != 0) ? y1 / y2 : 0.0;
            gMeanYieldModelsRatio[i]->SetPoint(j, x1, ratio);
            double err1 = gMeanYieldModels1[i]->GetErrorY(j);
            double err2 = gMeanYieldModels2[i]->GetErrorY(j);
            double ratioErr = (y2 != 0) ? ratio * sqrt(pow(err1 / y1, 2) + pow(err2 / y2, 2)) : 0.0;
            gMeanYieldModelsRatio[i]->SetPointError(j, 0.0, ratioErr);
        }
    }

    TGraphErrors *gMeanYieldRun3Data1 = (TGraphErrors *)fData1->Get("gMeanYieldRun3");
    TGraphErrors *gMeanYieldRun3Data2 = (TGraphErrors *)fData2->Get("gMeanYieldRun3");
    if (gMeanYieldRun3Data1 == nullptr || gMeanYieldRun3Data2 == nullptr)
    {
        cout << "Error: Data <p_T> graphs not found" << endl;
        return;
    }
    TGraphErrors *gMeanYieldRun3DataRatio = (TGraphErrors *)gMeanYieldRun3Data1->Clone("Data_Ratio");
    gMeanYieldRun3DataRatio->SetTitle("Run 3 Proton/Pion Ratio");
    for (int j = 0; j < gMeanYieldRun3Data1->GetN(); ++j)
    {
        double x1, y1, x2, y2;
        gMeanYieldRun3Data1->GetPoint(j, x1, y1);
        gMeanYieldRun3Data2->GetPoint(j, x2, y2);
        if (fabs(x1 - x2) > 1e-5)
        {
            cout << "Error: Mismatched x values for data at point " << j << ": " << x1 << " vs " << x2 << endl;
            return;
        }
        double ratio = (y2 != 0) ? y1 / y2 : 0.0;
        gMeanYieldRun3DataRatio->SetPoint(j, x1, ratio);
        double err1 = gMeanYieldRun3Data1->GetErrorY(j);
        double err2 = gMeanYieldRun3Data2->GetErrorY(j);
        double ratioErr = (y2 != 0) ? ratio * sqrt(pow(err1 / y1, 2) + pow(err2 / y2, 2)) : 0.0;
        gMeanYieldRun3DataRatio->SetPointError(j, 0.0, ratioErr);
    }

    gStyle->SetCanvasPreferGL(kTRUE);

    TCanvas *cMeanYieldvsNch = new TCanvas("cMeanYieldvsNch", "cMeanYieldvsNch", 720, 720);
    SetCanvasStyle(cMeanYieldvsNch, 0.15, 0.03, 0.03, 0.15);
    SetGraphStyle(gMeanYieldRun3DataRatio);
    gPad->SetLogx();
    gMeanYieldRun3DataRatio->SetTitle(0);
    // gMeanYieldRun3DataRatio->GetYaxis()->SetTitle("(K^{+} + K^{-}) / (#pi^{+} + #pi^{-})");
    gMeanYieldRun3DataRatio->GetYaxis()->SetTitle("(p + #bar{p}) / (#pi^{+} + #pi^{-})");
    gMeanYieldRun3DataRatio->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMeanYieldRun3DataRatio->SetMarkerStyle(20);
    gMeanYieldRun3DataRatio->GetYaxis()->SetTitleOffset(1.7);
    gMeanYieldRun3DataRatio->SetMarkerSize(1.2);
    gMeanYieldRun3DataRatio->SetMarkerColor(kRed);
    gMeanYieldRun3DataRatio->SetLineColor(kRed);
    // gMeanYieldRun3DataRatio->GetYaxis()->SetRangeUser(0, 0.79);
    gMeanYieldRun3DataRatio->SetMinimum(0.01);
    gMeanYieldRun3DataRatio->SetMaximum(0.18);
    gMeanYieldRun3DataRatio->GetXaxis()->SetLimits(0, 28.9);
    gMeanYieldRun3DataRatio->Draw("AP");

    TLegend *legMeanYield = new TLegend(0.2, 0.18, 0.45, 0.45);
    legMeanYield->SetTextSize(0.04);
    legMeanYield->SetBorderSize(0);
    legMeanYield->SetFillStyle(0);
    legMeanYield->AddEntry(gMeanYieldRun3DataRatio, "Data", "p");
    for (size_t i = 0; i < modelNames.size(); ++i)
    {
        SetGrapherrorStyle(gMeanYieldModelsRatio[i]);
        gMeanYieldModelsRatio[i]->SetMarkerStyle(20 + i);
        gMeanYieldModelsRatio[i]->SetMarkerSize(1.2);
        gMeanYieldModelsRatio[i]->SetMarkerColor(colors[i]);
        gMeanYieldModelsRatio[i]->SetLineColor(colors[i]);
        gMeanYieldModelsRatio[i]->SetFillColorAlpha(colors[i], 0.5);
        gMeanYieldModelsRatio[i]->Draw("E3 same");
        gMeanYieldModelsRatio[i]->Draw("l same");
        legMeanYield->AddEntry(gMeanYieldModelsRatio[i], modelNames[i].c_str(), "l");
    }
    legMeanYield->Draw();
}