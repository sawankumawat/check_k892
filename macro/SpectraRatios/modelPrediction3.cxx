#include <iostream>
#include <iomanip>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/initializations.h"

void modelPrediction3()
{
    //_id53739 for Pythia Monash with rescattering
    // TFile *fmodel = new TFile("ModelRootFiles/Pythia_Monash.root", "read");
    TFile *fmodel = new TFile("ModelRootFiles/EPOS.root", "read");
    // TFile *fmodel = new TFile("/home/sawan/Downloads/epos.root", "read");

    if (fmodel->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }

    vector<string> particleTypes = {"Kstar", "Phi", "Pion", "Kaon", "Proton"};
    int selectParticle = 1;

    TH2D *hNch05VsFT0M = (TH2D *)fmodel->Get("mc-particle-prediction/multiplicity/vsETA05/FT0AC");                                             // ETA0p5 (y) vs FT0M (x)
    TH2D *hParticlevsFT0M = (TH2D *)fmodel->Get(Form("mc-particle-prediction/prediction/pt/FT0AC/%s", particleTypes[selectParticle].c_str())); // FT0M (y) vs pT (x) for the selected particle
    if (hNch05VsFT0M == nullptr || hParticlevsFT0M == nullptr)
    {
        cout << "Error: histograms not found" << endl;
        return;
    }

    // float percentilesRun3[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    float percentilesRun3[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    const int nPercentiles = sizeof(percentilesRun3) / sizeof(percentilesRun3[0]);
    TGraphErrors *gMeanpTvsNch = new TGraphErrors(nPercentiles - 1);
    TGraphErrors *gMeanYield = new TGraphErrors(nPercentiles - 1);

    TH1D *hFT0M = (TH1D *)hNch05VsFT0M->ProjectionX("hFT0M");
    double totalArea = hFT0M->Integral();
    vector<int> boundaryBins(nPercentiles, -1);
    // 0% boundary = highest FT0M bin
    boundaryBins[0] = hFT0M->GetNbinsX();

    double cumulativePercent = 0.0;
    int percentileIndex = 1;

    for (int ibin = hFT0M->GetNbinsX(); ibin >= 1 && percentileIndex < nPercentiles; --ibin)
    {
        cumulativePercent += hFT0M->GetBinContent(ibin) * 100.0 / totalArea;

        while (percentileIndex < nPercentiles && cumulativePercent >= percentilesRun3[percentileIndex])
        {
            boundaryBins[percentileIndex] = ibin;
            percentileIndex++;
        }
    }

    // 100% boundary = lowest FT0M bin
    boundaryBins[nPercentiles - 1] = 1;

    // for (int i = 0; i < nPercentiles - 1; i++)
    // {
    //     cout << percentilesRun3[i] << "-" << percentilesRun3[i + 1] << "% : bins " << boundaryBins[i] << " -> " << boundaryBins[i + 1] << endl;
    // }
    // cout << endl;

    for (int icent = 0; icent < nPercentiles - 1; icent++)
    {
        int binHigh = boundaryBins[icent];    // larger FT0M
        int binLow = boundaryBins[icent + 1]; // smaller FT0M

        // Project Nch distribution for this centrality class
        TH1D *hNch = hNch05VsFT0M->ProjectionY(Form("hNch_%d", icent), binLow, binHigh);

        double meanNch = hNch->GetMean();

        // number of events in this FT0M interval
        double nEvents = hFT0M->Integral(binLow, binHigh);

        // pT spectrum in this FT0M interval
        TH1D *hPt = hParticlevsFT0M->ProjectionX(Form("hPt_%d", icent), binLow, binHigh);
        hPt->Scale(1.0 / nEvents);

        double frac = nEvents / totalArea;
        double BR = 0.666;                   // for kstar, for example
        double yield = hPt->Integral() / BR; // yield per event, corrected for branching ratio
        double meanPt = hPt->GetMean();

        gMeanpTvsNch->SetPoint(icent, meanNch, meanPt);
        gMeanpTvsNch->SetPointError(icent, 0, 0);
        gMeanYield->SetPoint(icent, meanNch, yield);
        gMeanYield->SetPointError(icent, 0, 0);

        cout << percentilesRun3[icent] << "-" << percentilesRun3[icent + 1] << "%   " << "Mean Nch = " << meanNch << "   Events = " << nEvents << "  <pT> = " << meanPt << "  Yield = " << yield << endl;
    }

    string KstarPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/";
    // TFile *fData = new TFile((KstarPath + "Results.root").c_str(), "read");
    TFile *fData = new TFile("PiKp_Run3_Results/Sawan/Pi_results.root", "read");

    if (fData->IsZombie())
    {
        cout << "Error: Data file not found" << endl;
        return;
    }
    TGraphErrors *gMeanpTRun3Data = (TGraphErrors *)fData->Get("gMeanpTRun3");
    TGraphErrors *gMeanYieldRun3Data = (TGraphErrors *)fData->Get("gMeanYieldRun3");

    if (gMeanpTRun3Data == nullptr)
    {
        cout << "Error: Data <p_T> graphs not found" << endl;
        return;
    }
    gStyle->SetCanvasPreferGL(kTRUE);

    // TCanvas *cMeanpTvsNch = new TCanvas("cMeanpTvsNch", "cMeanpTvsNch", 720, 720);
    // SetCanvasStyle(cMeanpTvsNch, 0.15, 0.03, 0.03, 0.15);
    // SetGraphStyle(gMeanpTvsNch);
    // gMeanpTvsNch->SetMarkerStyle(20);
    // gMeanpTvsNch->SetFillColorAlpha(kGreen + 2, 0.6);
    // gMeanpTvsNch->SetLineColor(kGreen + 2);
    // gMeanpTvsNch->SetLineWidth(2);
    // gMeanpTRun3Data->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/c)");
    // gMeanpTRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    // gMeanpTRun3Data->SetMarkerStyle(24);
    // gMeanpTRun3Data->SetMarkerSize(1.2);
    // gMeanpTRun3Data->SetMarkerColor(kRed);
    // gMeanpTRun3Data->SetLineColor(kRed);
    // gMeanpTRun3Data->GetYaxis()->SetRangeUser(0.31, 1.79);
    // gMeanpTRun3Data->GetXaxis()->SetLimits(0, 28.9);
    // gMeanpTRun3Data->Draw("AP");
    // gMeanpTvsNch->Draw("E3 same");
    // gMeanpTvsNch->Draw("l same");

    // TCanvas *cMeanYieldvsNch = new TCanvas("cMeanYieldvsNch", "cMeanYieldvsNch", 720, 720);
    // SetCanvasStyle(cMeanYieldvsNch, 0.15, 0.03, 0.03, 0.15);
    // SetGraphStyle(gMeanYield);
    // gMeanYield->SetMarkerStyle(20);
    // gMeanYield->SetFillColorAlpha(kGreen + 2, 0.6);
    // gMeanYield->SetLineColor(kGreen + 2);
    // gMeanYield->SetLineWidth(2);
    // gMeanYieldRun3Data->GetYaxis()->SetTitle("dN/dy");
    // gMeanYieldRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    // gMeanYieldRun3Data->SetMarkerStyle(24);
    // gMeanYieldRun3Data->SetMarkerSize(1.2);
    // gMeanYieldRun3Data->SetMarkerColor(kRed);
    // gMeanYieldRun3Data->SetLineColor(kRed);
    // gMeanYieldRun3Data->GetYaxis()->SetRangeUser(0, 1.89);
    // gMeanYieldRun3Data->GetXaxis()->SetLimits(0, 28.9);
    // gMeanYieldRun3Data->Draw("AP");
    // gMeanYield->Draw("E3 same");
    // gMeanYield->Draw("l same");

    // cout << "\nParticle type: " << particleTypes[selectParticle] << endl;
    // cout << "Data file used: " << fData->GetName() << endl;
}