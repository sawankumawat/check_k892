#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void systematicsTable()
{
    TFile *file = TFile::Open("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots/SysUncert.root", "READ");

    if (!file || file->IsZombie())
    {
        cout << "Error: Cannot open file " << endl;
        return;
    }

    // ============================================================
    // Multiplicity classes
    // ============================================================

    vector<int> multLow = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    vector<int> multHigh = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};

    const int nMult = multLow.size();

    // ============================================================
    // Get first histogram to determine pT binning
    // ============================================================

    TH1D *hFirst =
        (TH1D *)file->Get(
            Form("hSignalExtTotalSysSmoothed_%d_%d",
                 multLow[0], multHigh[0]));

    if (!hFirst)
    {
        cout << "Error: Could not find first histogram." << endl;
        file->Close();
        return;
    }

    int nPtBins = hFirst->GetNbinsX();

    // ============================================================
    // STEP 1
    //
    // Calculate multiplicity-averaged uncertainty for every pT bin
    // ============================================================

    vector<double> avgSignal(nPtBins, 0.);
    vector<double> avgTrack(nPtBins, 0.);
    vector<double> avgPID(nPtBins, 0.);
    vector<double> avgMaterial(nPtBins, 0.);
    vector<double> avgHadronic(nPtBins, 0.);
    vector<double> avgTotal(nPtBins, 0.);

    vector<int> nValid(nPtBins, 0);

    for (int ipt = 1; ipt <= nPtBins; ipt++)
    {
        double sumSignal = 0.;
        double sumTrack = 0.;
        double sumPID = 0.;
        double sumMaterial = 0.;
        double sumHadronic = 0.;
        double sumTotal = 0.;

        int nValidBin = 0;

        // --------------------------------------------------------
        // Loop over multiplicity classes
        // --------------------------------------------------------

        for (int im = 0; im < nMult; im++)
        {
            int low = multLow[im];
            int high = multHigh[im];

            TH1D *hSignal =
                (TH1D *)file->Get(
                    Form("hSignalExtTotalSysSmoothed_%d_%d",
                         low, high));

            TH1D *hTrack =
                (TH1D *)file->Get(
                    Form("hTrackSelTotalSysSmoothed_%d_%d",
                         low, high));

            TH1D *hPID =
                (TH1D *)file->Get(
                    Form("hPIDTotalSysSmoothed_%d_%d",
                         low, high));

            TH1D *hMaterial =
                (TH1D *)file->Get(
                    Form("hMaterialBudgetTotalSysSmoothed_%d_%d",
                         low, high));

            TH1D *hHadronic =
                (TH1D *)file->Get(
                    Form("hHadronicInteractionTotalSysSmoothed_%d_%d",
                         low, high));

            TH1D *hTotal =
                (TH1D *)file->Get(
                    Form("hTotalSysSmoothed_%d_%d",
                         low, high));

            if (!hSignal || !hTrack || !hPID ||
                !hMaterial || !hHadronic || !hTotal)
            {
                cout << "Warning: Missing histogram for "
                     << low << "-" << high
                     << "% multiplicity at pT bin "
                     << ipt << endl;

                continue;
            }

            sumSignal += hSignal->GetBinContent(ipt);
            sumTrack += hTrack->GetBinContent(ipt);
            sumPID += hPID->GetBinContent(ipt);
            sumMaterial += hMaterial->GetBinContent(ipt);
            sumHadronic += hHadronic->GetBinContent(ipt);
            sumTotal += hTotal->GetBinContent(ipt);

            nValidBin++;
        }

        // --------------------------------------------------------
        // Multiplicity average
        // --------------------------------------------------------

        if (nValidBin > 0)
        {
            avgSignal[ipt - 1] = sumSignal / nValidBin;
            avgTrack[ipt - 1] = sumTrack / nValidBin;
            avgPID[ipt - 1] = sumPID / nValidBin;
            avgMaterial[ipt - 1] = sumMaterial / nValidBin;
            avgHadronic[ipt - 1] = sumHadronic / nValidBin;
            avgTotal[ipt - 1] = sumTotal / nValidBin;

            nValid[ipt - 1] = nValidBin;
        }
    }

    // ============================================================
    // Print multiplicity-averaged uncertainty for every pT bin
    // ============================================================

    cout << endl;
    cout << "=========================================================================="
         << endl;

    cout << "Multiplicity-averaged systematic uncertainty"
         << endl;

    cout << "=========================================================================="
         << endl;

    cout << "pT range"
         << "\tSignal"
         << "\tTrack"
         << "\tPID"
         << "\tMaterial"
         << "\tHadronic"
         << "\tTotal"
         << endl;

    for (int ipt = 1; ipt <= nPtBins; ipt++)
    {
        double lowPt =
            hFirst->GetBinLowEdge(ipt);

        double highPt =
            hFirst->GetBinLowEdge(ipt + 1);

        cout << Form("%.2f-%.2f", lowPt, highPt)
             << "\t"
             << Form(" %.2f", avgSignal[ipt - 1] * 100.)
             << "\t"
             << Form("%.2f", avgTrack[ipt - 1] * 100.)
             << "\t"
             << Form("%.2f", avgPID[ipt - 1] * 100.)
             << "\t"
             << Form("%.2f", avgMaterial[ipt - 1] * 100.)
             << "\t"
             << Form("%.2f", avgHadronic[ipt - 1] * 100.)
             << "\t"
             << Form("%.2f", avgTotal[ipt - 1] * 100.)
             << endl;
    }

    // ============================================================
    // STEP 2
    //
    // Average pT bins into:
    //
    // 0 - 1.2
    // 1.2 - 3.0
    // 3.0 - 20.0
    //
    // ============================================================

    // ============================================================
    // pT ranges
    // ============================================================

    const int nPtRanges = 3;

    double ptLow[nPtRanges] = {0.0, 1.0, 4.0};
    double ptHigh[nPtRanges] = {1.0, 4.0, 20.0};

    TString ptRangeName[nPtRanges] =
        {
            "0-1.0",
            "1.0-4.0",
            "4.0-20"};

    cout << endl;
    cout << "=========================================================================="
         << endl;

    cout << "Multiplicity-averaged AND pT-averaged systematic uncertainty"
         << endl;

    cout << "=========================================================================="
         << endl;

    cout << "pT range"
         << "\tSignal"
         << "\tTrack"
         << "\tPID"
         << "\tMaterial"
         << "\tHadronic"
         << "\tTotal"
         << endl;

    for (int irange = 0; irange < nPtRanges; irange++)
    {
        double sumSignal = 0.;
        double sumTrack = 0.;
        double sumPID = 0.;
        double sumMaterial = 0.;
        double sumHadronic = 0.;
        double sumTotal = 0.;

        int nBins = 0;

        for (int ipt = 1; ipt <= nPtBins; ipt++)
        {
            double binCenter =
                hFirst->GetBinCenter(ipt);

            if (binCenter >= ptLow[irange] &&
                binCenter < ptHigh[irange])
            {
                sumSignal += avgSignal[ipt - 1];
                sumTrack += avgTrack[ipt - 1];
                sumPID += avgPID[ipt - 1];
                sumMaterial += avgMaterial[ipt - 1];
                sumHadronic += avgHadronic[ipt - 1];
                sumTotal += avgTotal[ipt - 1];

                nBins++;
            }
        }

        // --------------------------------------------------------
        // pT average
        // --------------------------------------------------------

        double finalSignal = (nBins > 0) ? sumSignal / nBins : 0.;
        double finalTrack = (nBins > 0) ? sumTrack / nBins : 0.;
        double finalPID = (nBins > 0) ? sumPID / nBins : 0.;
        double finalMaterial = (nBins > 0) ? sumMaterial / nBins : 0.;
        double finalHadronic = (nBins > 0) ? sumHadronic / nBins : 0.;
        double finalTotal = (nBins > 0) ? sumTotal / nBins : 0.;
        // --------------------------------------------------------
        // Print final result
        // --------------------------------------------------------

        cout << ptRangeName[irange]
             << "\t"
             << Form("          %.2f", finalSignal * 100.)
             << "\t"
             << Form("%.2f", finalTrack * 100.)
             << "\t"
             << Form("%.2f", finalPID * 100.)
             << "\t"
             << Form("%.2f", finalMaterial * 100.)
             << "\t"
             << Form("         %.2f", finalHadronic * 100.)
             << "\t"
             << Form("       %.2f", finalTotal * 100.)
             << endl;
    }

    file->Close();
}