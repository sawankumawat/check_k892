#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TMath.h"

using namespace std;

void SaveModelPrediction2()
{
    const vector<string> modelNames = {"EPOS", "Pythia_CR", "Pythia_Monash", "Pythia_Ropes", "Pythia_Shoving", "Pythia_Monash_Rescattering"};
    const vector<string> particleTypes = {"Kstar", "Phi", "Pion", "Kaon", "Proton", "PionMinus", "KaonMinus", "AntiProton", "Xi1530", "Kshort", "KstarPM"};
    const map<string, double> branchingRatios = {
        {"Kstar", 0.666},
        {"Phi", 0.492},
        {"Pion", 1.0},
        {"Kaon", 1.0},
        {"Proton", 1.0},
        {"PionMinus", 1.0},
        {"KaonMinus", 1.0},
        {"AntiProton", 1.0},
        {"Xi1530", 0.666},
        {"Kshort", 1.0},
        {"KstarPM", 0.333}};

    float percentilesRun3[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    const int nPercentiles = sizeof(percentilesRun3) / sizeof(percentilesRun3[0]);

    TFile *fOutput = new TFile("ModelResults.root", "recreate");
    if (fOutput->IsZombie())
    {
        cout << "Error: cannot create output file ModelResults2.root" << endl;
        return;
    }

    for (const auto &modelName : modelNames)
    {
        TFile *fModel = new TFile(Form("ModelRootFiles/%s.root", modelName.c_str()), "read");
        if (fModel->IsZombie())
        {
            cout << "Error: file not found for model " << modelName << endl;
            delete fModel;
            continue;
        }

        TH2D *hNch05VsFT0M = (TH2D *)fModel->Get("mc-particle-prediction/multiplicity/vsETA05/FT0AC");
        if (hNch05VsFT0M == nullptr)
        {
            cout << "Error: multiplicity histogram not found for model " << modelName << endl;
            fModel->Close();
            delete fModel;
            continue;
        }

        TH1D *hFT0M = (TH1D *)hNch05VsFT0M->ProjectionX(Form("hFT0M_%s", modelName.c_str()));
        hFT0M->SetDirectory(nullptr);
        const double totalArea = hFT0M->Integral();
        if (totalArea <= 0)
        {
            cout << "Warning: zero FT0M integral for model " << modelName << endl;
            delete hFT0M;
            fModel->Close();
            delete fModel;
            continue;
        }

        vector<int> boundaryBins(nPercentiles, -1);
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

        boundaryBins[nPercentiles - 1] = 1;

        // Fill any unresolved boundaries with nearest valid bin to avoid invalid projections.
        for (int i = 1; i < nPercentiles - 1; ++i)
        {
            if (boundaryBins[i] < 1)
            {
                boundaryBins[i] = boundaryBins[i - 1];
            }
        }

        fOutput->cd();
        TDirectory *modelDir = fOutput->mkdir(modelName.c_str());
        if (modelDir == nullptr)
        {
            cout << "Error: cannot create output directory for model " << modelName << endl;
            delete hFT0M;
            fModel->Close();
            delete fModel;
            continue;
        }

        cout << "\nModel: " << modelName << endl;
        for (int i = 0; i < nPercentiles - 1; ++i)
        {
            cout << percentilesRun3[i] << "-" << percentilesRun3[i + 1] << "% : bins " << boundaryBins[i] << " -> " << boundaryBins[i + 1] << endl;
        }

        for (const auto &particleType : particleTypes)
        {
            TH2D *hParticlevsFT0M = (TH2D *)fModel->Get(Form("mc-particle-prediction/prediction/pt/FT0AC/%s", particleType.c_str()));
            if (hParticlevsFT0M == nullptr)
            {
                cout << "Warning: particle histogram not found for " << modelName << "/" << particleType << endl;
                continue;
            }

            modelDir->cd();
            TDirectory *particleDir = modelDir->mkdir(particleType.c_str());
            if (particleDir == nullptr)
            {
                cout << "Warning: cannot create particle directory for " << modelName << "/" << particleType << endl;
                continue;
            }
            particleDir->cd();

            TGraphErrors *gMeanpTvsNch = new TGraphErrors(nPercentiles - 1);
            TGraphErrors *gMeanYieldvsNch = new TGraphErrors(nPercentiles - 1);

            const double br = branchingRatios.count(particleType) ? branchingRatios.at(particleType) : 1.0;

            TH1D *hPtMinimumBias = hParticlevsFT0M->ProjectionX(Form("hPt_%s_MinBias_%s", modelName.c_str(), particleType.c_str()), 1, hParticlevsFT0M->GetYaxis()->GetNbins());
            hPtMinimumBias->SetDirectory(nullptr);
            const double nEventsMinimumBias = hFT0M->Integral();
            if (nEventsMinimumBias > 0)
            {
                hPtMinimumBias->Scale(1.0 / nEventsMinimumBias);
            }
            hPtMinimumBias->Write(Form("hPt_%s_MinBias", particleType.c_str()));

            for (int icent = 0; icent < nPercentiles - 1; ++icent)
            {
                int binHigh = boundaryBins[icent];
                int binLow = boundaryBins[icent + 1];
                if (binLow > binHigh)
                {
                    std::swap(binLow, binHigh);
                }

                TH1D *hNch = hNch05VsFT0M->ProjectionY(Form("hNch_%s_%s_%d", modelName.c_str(), particleType.c_str(), icent), binLow, binHigh);
                hNch->SetDirectory(nullptr);
                const double meanNch = hNch->GetMean();

                const double nEvents = hFT0M->Integral(binLow, binHigh);
                if (nEvents <= 0)
                {
                    gMeanpTvsNch->SetPoint(icent, meanNch, 0.0);
                    gMeanpTvsNch->SetPointError(icent, 0.0, 0.0);
                    gMeanYieldvsNch->SetPoint(icent, meanNch, 0.0);
                    gMeanYieldvsNch->SetPointError(icent, 0.0, 0.0);
                    delete hNch;
                    continue;
                }

                TH1D *hPt = hParticlevsFT0M->ProjectionX(Form("hPt_%s_%s_%d", modelName.c_str(), particleType.c_str(), icent), binLow, binHigh);
                hPt->SetDirectory(nullptr);
                hPt->Scale(1.0 / nEvents);
                hPt->Write(Form("hPt_%s_%d_%0.0f_%0.0f", particleType.c_str(), icent, percentilesRun3[icent], percentilesRun3[icent + 1]));

                const double meanPt = hPt->GetMean();
                const double yield = hPt->Integral() / br;

                gMeanpTvsNch->SetPoint(icent, meanNch, meanPt);
                gMeanpTvsNch->SetPointError(icent, 0.0, 0.0);
                gMeanYieldvsNch->SetPoint(icent, meanNch, yield);
                gMeanYieldvsNch->SetPointError(icent, 0.0, 0.0);

                delete hPt;
                delete hNch;
            }

            gMeanpTvsNch->Write(Form("gMeanpT_%s", particleType.c_str()));
            gMeanYieldvsNch->Write(Form("gMeanYield_%s", particleType.c_str()));

            delete hPtMinimumBias;
            delete gMeanpTvsNch;
            delete gMeanYieldvsNch;
        }

        delete hFT0M;
        fModel->Close();
        delete fModel;
    }

    fOutput->Close();
    delete fOutput;
}