#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TH1.h"
#include "TDirectory.h"

void systematics_rebin()
{
    string path = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots";
    TString inputFile = path + "/SysUncert.root";
    TString outputFile = path + "/SysUncert_coarse.root";

    double ptBins[] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 10.0, 12.0, 15.0, 20.0};
    const int nBins = sizeof(ptBins) / sizeof(ptBins[0]) - 1;

    const char *histNames[] = {
        "Fit_range",
        "Norm_range",
        "Combinatorial_bkg",
        "Residual_bkg",
        "Width_fix_free",
        "PID",
        "TrkSel",
        "Hadronic_interaction",
        "Material_budget"};

    const int nSources = sizeof(histNames) / sizeof(histNames[0]);

    // TString shortNames[] = {"Fit", "Norm", "CombBkg", "ResBkg", "Width", "PID", "DCA", "PVCont", "HI", "MatBgt"};
    TString shortNames[] = {"Fit", "Norm", "CombBkg", "ResBkg", "Width", "PID", "TrkSel", "HI", "MB"};

    TFile *fin = TFile::Open(inputFile, "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Error opening " << inputFile << std::endl;
        return;
    }

    TDirectory *dir = fin->GetDirectory("SystematicSources");
    if (!dir)
    {
        std::cerr << "SystematicSources folder not found!" << std::endl;
        fin->Close();
        return;
    }

    TFile *fout = TFile::Open(outputFile, "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "Error creating " << outputFile << std::endl;
        fin->Close();
        return;
    }

    TDirectory *outDir = fout->mkdir("SystematicSources");
    outDir->cd();

    std::vector<TH1 *> hCoarse(nSources, nullptr);

    for (int i = 0; i < nSources; i++)
    {
        TH1 *h = dynamic_cast<TH1 *>(dir->Get(histNames[i]));

        if (!h)
        {
            std::cerr << "Histogram not found: " << histNames[i] << std::endl;
            continue;
        }

        hCoarse[i] = new TH1D(histNames[i], h->GetTitle(), nBins, ptBins);

        for (int newBin = 1; newBin <= nBins; newBin++)
        {
            double low = ptBins[newBin - 1];
            double high = ptBins[newBin];

            double weightedSum = 0.0;
            double totalWidth = 0.0;

            for (int oldBin = 1; oldBin <= h->GetNbinsX(); oldBin++)
            {
                double oldLow = h->GetXaxis()->GetBinLowEdge(oldBin);
                double oldHigh = h->GetXaxis()->GetBinUpEdge(oldBin);

                if (oldLow >= low && oldHigh <= high)
                {
                    double value = h->GetBinContent(oldBin);
                    double width = oldHigh - oldLow;

                    weightedSum += value * width;
                    totalWidth += width;
                }
            }

            double average = totalWidth > 0.0 ? weightedSum / totalWidth : 0.0;
            hCoarse[i]->SetBinContent(newBin, average);
        }

        hCoarse[i]->SetDirectory(outDir);
        hCoarse[i]->Write(histNames[i], TObject::kOverwrite);

        std::cout << "Rebinned (width-weighted average): " << histNames[i] << std::endl;
    }

    // ============================================================
    // Print systematic uncertainty table
    // ============================================================

    std::cout << "\n\n========================================================================================\n";
    std::cout << "               SYSTEMATIC UNCERTAINTIES (0-100%) Run 3 data (pp 13.6 TeV) \n";
    std::cout << "========================================================================================\n";

    std::cout << std::left << std::setw(8) << "pT";
    for (int i = 0; i < nSources; i++)
        std::cout << std::right << std::setw(8) << shortNames[i];
    std::cout << std::right << std::setw(8) << "Total";
    std::cout << "\n";

    std::cout << std::string(8 + 8 * (nSources + 1), '-') << "\n";

    // Values for each pT bin
    for (int bin = 1; bin <= nBins; bin++)
    {
        double ptLow = ptBins[bin - 1];
        double ptHigh = ptBins[bin];
        std::string ptRange = Form("%.1f-%.1f", ptLow, ptHigh);

        std::cout << std::left << std::setw(8) << ptRange;

        double totalSquared = 0.0;

        for (int i = 0; i < nSources; i++)
        {
            double value = hCoarse[i] ? hCoarse[i]->GetBinContent(bin) * 100.0 : 0.0;
            totalSquared += value * value;
            std::cout << std::right << std::setw(8) << std::fixed << std::setprecision(1) << value;
        }

        double total = std::sqrt(totalSquared);
        std::cout << std::right << std::setw(8) << std::fixed << std::setprecision(1) << total;
        std::cout << "\n";
    }

    // ============================================================
    // Average uncertainty over all pT bins
    // ============================================================

    std::cout << std::string(8 + 8 * (nSources + 1), '-') << "\n";
    std::cout << std::left << std::setw(8) << "Average";

    for (int i = 0; i < nSources; i++)
    {
        double sum = 0.0;
        int count = 0;

        if (hCoarse[i])
        {
            for (int bin = 1; bin <= nBins; bin++)
            {
                sum += hCoarse[i]->GetBinContent(bin);
                count++;
            }
        }

        double average = (count > 0) ? sum / count : 0.0;
        std::cout << std::right << std::setw(8) << std::fixed << std::setprecision(1) << average * 100.0;
    }

    // Average total uncertainty
    double totalSum = 0.0;

    for (int bin = 1; bin <= nBins; bin++)
    {
        double totalSquared = 0.0;

        for (int i = 0; i < nSources; i++)
        {
            if (hCoarse[i])
            {
                double value = hCoarse[i]->GetBinContent(bin) * 100.0;
                totalSquared += value * value;
            }
        }

        totalSum += std::sqrt(totalSquared);
    }

    double averageTotal = totalSum / nBins;
    std::cout << std::right << std::setw(8) << std::fixed << std::setprecision(1) << averageTotal;
    std::cout << "\n";

    std::cout << "========================================================================================\n";

    fout->Close();
    fin->Close();

    std::cout << "\nCreated: " << outputFile << std::endl;
}