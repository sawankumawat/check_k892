#include <iostream>
#include <iomanip>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/initializations.h"

int colors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

void compareDataModel()
{
    // TFile *fData = new TFile("levy_fit.root", "read");
    TFile *fData = new TFile("PiKp_Run3_Results/Sawan/Ka_results.root", "read");
    TFile *fModel = new TFile("ModelResults.root", "read");
    if (fData->IsZombie() || fModel->IsZombie())
    {
        cout << "Error: Data or Model file not found" << endl;
        return;
    }

    // const vector<string> modelNames = {"EPOS", "Pythia_Monash", "Pythia_CR", "Pythia_Ropes", "Pythia_Shoving"};
    const vector<string> modelNames = {"Pythia_Monash", "Pythia_CR", "Pythia_Ropes", "Pythia_Shoving"};
    TGraphErrors *gMeanpTModels[modelNames.size()];
    TGraphErrors *gMeanYieldModels[modelNames.size()];
    for (size_t i = 0; i < modelNames.size(); ++i)
    {
        TDirectory *modelDir = (TDirectory *)fModel->Get(modelNames[i].c_str());
        if (modelDir == nullptr)
        {
            cout << "Error: Model directory not found for " << modelNames[i] << endl;
            return;
        }
        gMeanpTModels[i] = (TGraphErrors *)modelDir->Get("Kstar/gMeanpT_Kstar");
        gMeanYieldModels[i] = (TGraphErrors *)modelDir->Get("Kstar/gMeanYield_Kstar");
        if (gMeanpTModels[i] == nullptr || gMeanYieldModels[i] == nullptr)
        {
            cout << "Error: Model graphs not found for " << modelNames[i] << endl;
            return;
        }
    }

    TGraphErrors *gMeanpTRun3Data = (TGraphErrors *)fData->Get("gMeanpTRun3");
    TGraphErrors *gMeanYieldRun3Data = (TGraphErrors *)fData->Get("gMeanYieldRun3");

    // // divide the value of gMeanYieldRun3Data by 2
    // for (int i = 0; i < gMeanYieldRun3Data->GetN(); ++i)
    // {
    //     double x, y;
    //     gMeanYieldRun3Data->GetPoint(i, x, y);
    //     gMeanYieldRun3Data->SetPoint(i, x, y / 2.0);
    //     double ex = gMeanYieldRun3Data->GetErrorX(i);
    //     double ey = gMeanYieldRun3Data->GetErrorY(i);
    //     gMeanYieldRun3Data->SetPointError(i, ex, ey / 2.0);
    // }

    // TGraphErrors *gMeanpTRun2Data = (TGraphErrors *)fData->Get("gMeanpTRun2");
    // TGraphErrors *gMeanYieldRun2Data = (TGraphErrors *)fData->Get("gMeanYieldRun2");

    if (gMeanpTRun3Data == nullptr)
    {
        cout << "Error: Data <p_T> graphs not found" << endl;
        return;
    }
    gStyle->SetCanvasPreferGL(kTRUE);

    TCanvas *cMeanpTvsNch = new TCanvas("cMeanpTvsNch", "cMeanpTvsNch", 720, 720);
    SetCanvasStyle(cMeanpTvsNch, 0.15, 0.03, 0.03, 0.15);
    SetGraphStyle(gMeanpTRun3Data);
    gMeanpTRun3Data->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/c)");
    gMeanpTRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanpTRun3Data->SetMarkerStyle(20);
    gMeanpTRun3Data->SetMarkerSize(1.2);
    gMeanpTRun3Data->SetMarkerColor(kRed);
    gMeanpTRun3Data->SetLineColor(kRed);
    gMeanpTRun3Data->GetYaxis()->SetRangeUser(0.41, 1.89);
    gMeanpTRun3Data->GetXaxis()->SetLimits(0, 28.9);
    gMeanpTRun3Data->Draw("AP");
    // gMeanpTRun2Data->SetMarkerStyle(25);
    // gMeanpTRun2Data->SetMarkerSize(1.2);
    // gMeanpTRun2Data->SetMarkerColor(kBlue);
    // gMeanpTRun2Data->SetLineColor(kBlue);
    // gMeanpTRun2Data->Draw("P same");
    TLegend *legMeanpT = new TLegend(0.2, 0.70, 0.45, 0.92);
    legMeanpT->SetTextSize(0.035);
    legMeanpT->SetBorderSize(0);
    legMeanpT->SetFillStyle(0);
    legMeanpT->AddEntry(gMeanpTRun3Data, "Run 3 (13.6 TeV)", "p");
    // legMeanpT->Add Entry(gMeanpTRun2Data, "Run 2 (13 TeV)", "p");
    for (size_t i = 0; i < modelNames.size(); ++i)
    {
        SetGrapherrorStyle(gMeanpTModels[i]);
        gMeanpTModels[i]->SetMarkerStyle(20 + i);
        gMeanpTModels[i]->SetMarkerSize(1.2);
        gMeanpTModels[i]->SetMarkerColor(colors[i]);
        gMeanpTModels[i]->SetLineColor(colors[i]);
        gMeanpTModels[i]->SetFillColorAlpha(colors[i], 0.5);
        gMeanpTModels[i]->Draw("E3 same");
        gMeanpTModels[i]->Draw("l");
        legMeanpT->AddEntry(gMeanpTModels[i], modelNames[i].c_str(), "l");
    }
    legMeanpT->Draw();

    TCanvas *cMeanYieldvsNch = new TCanvas("cMeanYieldvsNch", "cMeanYieldvsNch", 720, 720);
    SetCanvasStyle(cMeanYieldvsNch, 0.15, 0.03, 0.03, 0.15);
    SetGraphStyle(gMeanYieldRun3Data);
    gMeanYieldRun3Data->GetYaxis()->SetTitle("dN/dy");
    gMeanYieldRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanYieldRun3Data->SetMarkerStyle(20);
    gMeanYieldRun3Data->SetMarkerSize(1.2);
    gMeanYieldRun3Data->SetMarkerColor(kRed);
    gMeanYieldRun3Data->SetLineColor(kRed);
    gMeanYieldRun3Data->GetYaxis()->SetRangeUser(0, 0.69);
    gMeanYieldRun3Data->GetXaxis()->SetLimits(0, 28.9);
    gMeanYieldRun3Data->Draw("AP");
    // gMeanYieldRun2Data->SetMarkerStyle(25);
    // gMeanYieldRun2Data->SetMarkerSize(1.2);
    // gMeanYieldRun2Data->SetMarkerColor(kBlue);
    // gMeanYieldRun2Data->SetLineColor(kBlue);
    // gMeanYieldRun2Data->Draw("P same");
    TLegend *legMeanYield = new TLegend(0.2, 0.60, 0.45, 0.92);
    legMeanYield->SetTextSize(0.04);
    legMeanYield->SetBorderSize(0);
    legMeanYield->SetFillStyle(0);
    legMeanYield->AddEntry(gMeanYieldRun3Data, "Run 3 (13.6 TeV)", "p");
    // legMeanYield->AddEntry(gMeanYieldRun2Data, "Run 2 (13 TeV)", "p");
    for (size_t i = 0; i < modelNames.size(); ++i)
    {
        SetGrapherrorStyle(gMeanYieldModels[i]);
        gMeanYieldModels[i]->SetMarkerStyle(20 + i);
        gMeanYieldModels[i]->SetMarkerSize(1.2);
        gMeanYieldModels[i]->SetMarkerColor(colors[i]);
        gMeanYieldModels[i]->SetLineColor(colors[i]);
        gMeanYieldModels[i]->SetFillColorAlpha(colors[i], 0.5);
        gMeanYieldModels[i]->Draw("E3 same");
        gMeanYieldModels[i]->Draw("l same");
        legMeanYield->AddEntry(gMeanYieldModels[i], modelNames[i].c_str(), "l");
    }
    legMeanYield->Draw();
}