#include <iostream>
#include <iomanip>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/initializations.h"
#include "TDirectory.h"

void SaveModelPrediction()
{
    vector<string> modelNames = {"EPOS", "Pythia_CR", "Pythia_Monash", "Pythia_Ropes", "Pythia_Shoving"};
    vector<string> particleTypes = {"Kstar", "Phi", "Pion", "Kaon", "Proton"};
    TFile *fOutput = new TFile("ModelResults.root", "recreate");

    for (const auto &modelName : modelNames)
    {

        TFile *fmodel = new TFile(Form("ModelRootFiles/%s.root", modelName.c_str()), "read");
        if (fmodel->IsZombie())
        {
            cout << "Error: file not found" << endl;
            return;
        }

        TH2D *hNch05VsFT0M = (TH2D *)fmodel->Get("mc-particle-prediction/multiplicity/vsETA05/FT0AC");
        if (hNch05VsFT0M == nullptr)
        {
            cout << "Error: histograms not found" << endl;
            return;
        }
        double dnch_detaRun3[] = {21.78, 18.48, 15.76, 13.89, 12.50, 10.86, 9.09, 7.63, 5.87, 3.69};
        double dnch_detaRun3_err[] = {0.38, 0.25, 0.22, 0.19, 0.17, 0.15, 0.13, 0.11, 0.09, 0.06};
        int nBins = sizeof(dnch_detaRun3) / sizeof(dnch_detaRun3[0]);
        vector<double> xValuesFT0Mcorrespondingdnch_detaRun3;
        vector<double> ft0mBinCenter;
        vector<double> meanNchPerFt0mBin;
        // const int nXBins = hNch05VsFT0M->GetNbinsX();
        const int nXBins = 157; // After 157th bin, all entries are mostly 0.
        vector<int> nEventsFT0M = {0};

        // Build per-bin (not cumulative) mean Nch map
        for (int xbin = 1; xbin <= nXBins; ++xbin)
        {
            TH1D *projY = hNch05VsFT0M->ProjectionY(Form("py_%d", xbin), xbin, xbin); // Single bin projection
            ft0mBinCenter.push_back(hNch05VsFT0M->GetXaxis()->GetBinCenter(xbin));
            meanNchPerFt0mBin.push_back(projY->GetMean());
            // cout << "FT0M bin " << xbin << " (center=" << hNch05VsFT0M->GetXaxis()->GetBinCenter(xbin) << ") has mean Nch in |eta|<0.5 as " << projY->GetMean() << endl;
        }

        // For each experimental dNch/deta, find the FT0M multiplicity by interpolation
        cout << left
             << setw(15) << "dNch/deta"
             << setw(12) << "y1"
             << setw(10) << "y2"
             << setw(15) << "FT0M value"
             << setw(18) << "bin range"
             << setw(12) << "entries" << endl;
        cout << string(90, '-') << endl;
        TGraphErrors *gMeanpT = new TGraphErrors(nBins);
        TGraphErrors *gMeanYield = new TGraphErrors(nBins);

        fOutput->cd();
        TDirectory *modelDir = fOutput->mkdir(modelName.c_str());
        modelDir->cd();

        vector<double> matchedFt0m(nBins, -1.0);
        vector<double> matchedY1(nBins, -1.0);
        vector<double> matchedY2(nBins, -1.0);
        vector<int> matchedBin(nBins, -1);
        vector<double> targetNchVec(nBins, -1.0);

        for (int ich = nBins - 1; ich >= 0; --ich)
        {
            const double targetNch = dnch_detaRun3[ich];
            targetNchVec[ich] = targetNch;
            // int matchedBin = -1;
            // double matchedY1 = -1.0;
            // double matchedY2 = -1.0;
            // Find the bin where mean Nch is closest to target, or bracket for interpolation
            for (size_t ibin = 0; ibin + 1 < meanNchPerFt0mBin.size(); ++ibin)
            {
                const double y1 = meanNchPerFt0mBin[ibin];
                const double y2 = meanNchPerFt0mBin[ibin + 1];

                // Check if targetNch is between y1 and y2 (or equal to one of them)
                if ((targetNch >= y1 && targetNch <= y2) || (targetNch >= y2 && targetNch <= y1))
                {
                    const double x1 = ft0mBinCenter[ibin];
                    const double x2 = ft0mBinCenter[ibin + 1];
                    matchedY1[ich] = y1;
                    matchedY2[ich] = y2;

                    if (fabs(y2 - y1) > 1e-10) // Avoid division by zero
                    {
                        matchedFt0m[ich] = x1 + (targetNch - y1) * (x2 - x1) / (y2 - y1);
                    }
                    else
                    {
                        matchedFt0m[ich] = x1;
                    }
                    matchedBin[ich] = ibin;
                    break;
                }
            }
        }

        for (const auto &particleType : particleTypes)
        {
            TH2D *hParticlevsFT0M = (TH2D *)fmodel->Get(Form("mc-particle-prediction/prediction/pt/FT0AC/%s", particleType.c_str()));

            for (int ich = nBins - 1; ich >= 0; --ich)
            {
                TH1D *hPtDistribution = nullptr;

                if (matchedFt0m[ich] >= 0.0)
                {

                    const int matchedFt0mBin = hParticlevsFT0M->GetYaxis()->FindBin(matchedFt0m[ich]);
                    hPtDistribution = hParticlevsFT0M->ProjectionX(Form("hPtDist_dNch_%.2f", targetNchVec[ich]), matchedFt0mBin, matchedFt0mBin);
                    TH1F *hFT0M1Ddist = (TH1F *)hNch05VsFT0M->ProjectionY(Form("hFT0M1D_dNch_%.2f", targetNchVec[ich]), matchedBin[ich], matchedBin[ich]);
                    nEventsFT0M.push_back(hFT0M1Ddist->GetEntries());
                    int entriesTemp = hFT0M1Ddist->GetEntries();

                    xValuesFT0Mcorrespondingdnch_detaRun3.push_back(matchedFt0m[ich]);
                    // cout << left
                    //      << setw(12) << targetNch
                    //      << setw(12) << matchedY1[ich]
                    //      << setw(15) << matchedY2[ich]
                    //      << setw(15) << matchedFt0m[ich]
                    //      << setw(2) << matchedBin[ich] << "-" << setw(12) << matchedBin[ich] + 1
                    //      << setw(12) << entriesTemp << endl;

                    if (hPtDistribution != nullptr)
                    {
                        // cout << "    pT projection: entries = " << hPtDistribution->GetEntries()
                        //      << ", mean pT = " << hPtDistribution->GetMean()
                        //      << ", RMS = " << hPtDistribution->GetRMS() << endl;
                        gMeanpT->SetPoint(nBins - 1 - ich, targetNchVec[ich], hPtDistribution->GetMean());
                        gMeanpT->SetPointError(nBins - 1 - ich, 0, hPtDistribution->GetRMS() / sqrt(hPtDistribution->GetEntries()));
                        double meanYield = hPtDistribution->Integral() / nEventsFT0M.back(); // Normalize by number of events in that FT0M bin
                        gMeanYield->SetPoint(nBins - 1 - ich, targetNchVec[ich], meanYield);
                        gMeanYield->SetPointError(nBins - 1 - ich, 0, sqrt(meanYield * (1 + meanYield) / nEventsFT0M.back())); // Poisson error on yield per event
                    }
                }
            }
            gMeanYield->Write(Form("gMeanYield_%s", particleType.c_str()));
            gMeanpT->Write(Form("gMeanpT_%s", particleType.c_str()));
        }
    }
}