#include <iostream>
#include <iomanip>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/initializations.h"

void modelPrediction()
{
    //_id53739 for Pythia Monash with rescattering
    TFile *fmodel = new TFile("ModelRootFiles/Pythia_Monash.root", "read");
    // TFile *fmodel = new TFile("ModelRootFiles/EPOS.root", "read");
    if (fmodel->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }
    TFile *fEPOSLocal = new TFile("EPOS_output_NoUrQMD.root", "read");
    TGraphErrors *gKstarYield = (TGraphErrors *)fEPOSLocal->Get("kstar_vs_mult_urqmdOFF");
    if (gKstarYield == nullptr)
    {
        cout << "Error: histogram not found in EPOS local file" << endl;
    }

    vector<string> particleTypes = {"Kstar", "Phi", "Pion", "Kaon", "Proton"};
    int selectParticle = 0;

    TH2D *hNch05VsFT0M = (TH2D *)fmodel->Get("mc-particle-prediction/multiplicity/vsETA05/FT0AC");
    TH2D *hParticlevsFT0M = (TH2D *)fmodel->Get(Form("mc-particle-prediction/prediction/pt/FT0AC/%s", particleTypes[selectParticle].c_str()));
    if (hNch05VsFT0M == nullptr || hParticlevsFT0M == nullptr)
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
    TGraphErrors *gMeanpTvsNch = new TGraphErrors(nBins);
    TGraphErrors *gMeanYield = new TGraphErrors(nBins);

    for (int ich = nBins - 1; ich >= 0; --ich)
    {
        const double targetNch = dnch_detaRun3[ich];
        double matchedFt0m = -1.0;
        int matchedBin = -1;
        double matchedY1 = -1.0;
        double matchedY2 = -1.0;
        TH1D *hPtDistribution = nullptr;
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
                matchedY1 = y1;
                matchedY2 = y2;

                if (fabs(y2 - y1) > 1e-10) // Avoid division by zero
                {
                    matchedFt0m = x1 + (targetNch - y1) * (x2 - x1) / (y2 - y1);
                }
                else
                {
                    matchedFt0m = x1;
                }
                matchedBin = ibin;
                break;
            }
        }

        if (matchedFt0m >= 0.0)
        {
            const int matchedFt0mBin = hParticlevsFT0M->GetYaxis()->FindBin(matchedFt0m);
            hPtDistribution = hParticlevsFT0M->ProjectionX(Form("hPtDist_dNch_%.2f", targetNch), matchedFt0mBin, matchedFt0mBin);
            TH1F *hFT0M1Ddist = (TH1F *)hNch05VsFT0M->ProjectionY(Form("hFT0M1D_dNch_%.2f", targetNch), matchedBin, matchedBin);
            nEventsFT0M.push_back(hFT0M1Ddist->GetEntries());
            int entriesTemp = hFT0M1Ddist->GetEntries();

            xValuesFT0Mcorrespondingdnch_detaRun3.push_back(matchedFt0m);
            cout << left
                 << setw(12) << targetNch
                 << setw(12) << matchedY1
                 << setw(15) << matchedY2
                 << setw(15) << matchedFt0m
                 << setw(2) << matchedBin << "-" << setw(12) << matchedBin + 1
                 << setw(12) << entriesTemp << endl;

            if (hPtDistribution != nullptr)
            {
                // cout << "    pT projection: entries = " << hPtDistribution->GetEntries()
                //      << ", mean pT = " << hPtDistribution->GetMean()
                //      << ", RMS = " << hPtDistribution->GetRMS() << endl;
                gMeanpTvsNch->SetPoint(nBins - 1 - ich, targetNch, hPtDistribution->GetMean());
                gMeanpTvsNch->SetPointError(nBins - 1 - ich, 0, hPtDistribution->GetRMS() / sqrt(hPtDistribution->GetEntries()));
                double meanYield = hPtDistribution->Integral() / nEventsFT0M.back(); // Normalize by number of events in that FT0M bin
                gMeanYield->SetPoint(nBins - 1 - ich, targetNch, meanYield);
                gMeanYield->SetPointError(nBins - 1 - ich, 0, sqrt(meanYield * (1 + meanYield) / nEventsFT0M.back())); // Poisson error on yield per event
            }
        }
        else
        {
            cout << left
                 << setw(12) << targetNch
                 << setw(12) << "n/a"
                 << setw(12) << "n/a"
                 << setw(12) << "n/a"
                 << setw(12) << "n/a"
                 << setw(12) << "n/a"
                 << setw(12) << "not found" << endl;
        }
    }

    // TFile *fData = new TFile("../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/Results.root", "read");
    TFile *fData = new TFile("PiKp_Run3_Results/Sawan/Ka_results.root", "read");
    if (fData->IsZombie())
    {
        cout << "Error: Data file not found" << endl;
        return;
    }

    TGraphErrors *gMeanpTRun3Data = (TGraphErrors *)fData->Get("gMeanpTRun3");
    TGraphErrors *gMeanYieldRun3Data = (TGraphErrors *)fData->Get("gMeanYieldRun3");

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
    SetGraphStyle(gMeanpTvsNch);
    gMeanpTvsNch->SetMarkerStyle(20);
    // gMeanpTvsNch->SetMarkerSize(1.2);
    // gMeanpTvsNch->SetFillStyle(3002);
    // gMeanpTvsNch->SetMarkerColor(kGreen + 2);
    gMeanpTvsNch->SetFillColorAlpha(kGreen + 2, 0.6);
    gMeanpTvsNch->SetLineColor(kGreen + 2);
    gMeanpTvsNch->SetLineWidth(2);
    // gMeanpTvsNch->Draw("AP");
    gMeanpTRun3Data->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/c)");
    gMeanpTRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanpTRun3Data->SetMarkerStyle(24);
    gMeanpTRun3Data->SetMarkerSize(1.2);
    gMeanpTRun3Data->SetMarkerColor(kRed);
    gMeanpTRun3Data->SetLineColor(kRed);
    gMeanpTRun3Data->GetYaxis()->SetRangeUser(0.31, 1.79);
    gMeanpTRun3Data->GetXaxis()->SetLimits(0, 28.9);
    gMeanpTRun3Data->Draw("AP");
    gMeanpTvsNch->Draw("E3 same");
    gMeanpTvsNch->Draw("l same");

    // gMeanpTRun2Data->SetMarkerStyle(25);
    // gMeanpTRun2Data->SetMarkerSize(1.2);
    // gMeanpTRun2Data->SetMarkerColor(kBlue);
    // gMeanpTRun2Data->SetLineColor(kBlue);
    // gMeanpTRun2Data->Draw("P same");

    TCanvas *cMeanYieldvsNch = new TCanvas("cMeanYieldvsNch", "cMeanYieldvsNch", 720, 720);
    SetCanvasStyle(cMeanYieldvsNch, 0.15, 0.03, 0.03, 0.15);
    SetGraphStyle(gMeanYield);
    gMeanYield->SetMarkerStyle(20);
    // gMeanYield->SetFillStyle(3002);
    // gMeanYield->SetMarkerSize(1.2);
    // gMeanYield->SetMarkerColor(kGreen + 2);
    gMeanYield->SetFillColorAlpha(kGreen + 2, 0.6);
    gMeanYield->SetLineColor(kGreen + 2);
    gMeanYield->SetLineWidth(2);
    // gMeanYield->Draw("AP");
    gMeanYieldRun3Data->GetYaxis()->SetTitle("dN/dy");
    gMeanYieldRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanYieldRun3Data->SetMarkerStyle(24);
    gMeanYieldRun3Data->SetMarkerSize(1.2);
    gMeanYieldRun3Data->SetMarkerColor(kRed);
    gMeanYieldRun3Data->SetLineColor(kRed);
    gMeanYieldRun3Data->GetYaxis()->SetRangeUser(0, 1.89);
    gMeanYieldRun3Data->GetXaxis()->SetLimits(0, 28.9);
    gMeanYieldRun3Data->Draw("AP");
    gMeanYield->Draw("E3 same");
    gMeanYield->Draw("l same");
    // gKstarYield->Draw("l same");

    // gMeanYieldRun2Data->SetMarkerStyle(25);
    // gMeanYieldRun2Data->SetMarkerSize(1.2);
    // gMeanYieldRun2Data->SetMarkerColor(kBlue);
    // gMeanYieldRun2Data->SetLineColor(kBlue);
    // gMeanYieldRun2Data->Draw("P same");

    cout << "\nParticle type: " << particleTypes[selectParticle] << endl;
    cout << "Data file used: " << fData->GetName() << endl;
}