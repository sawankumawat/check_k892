#include <iostream>
#include <vector>
#include <cmath>
#include "src/style.h"
#include "src/initializations.h"
#include "SystematicHelper/helper.cxx"

void openTFile(TFile *&file, const string &path);
void openTH1D(TH1D *&hist, TFile *file, const string &histPath);
TH1D *RebinHistogram(TH1 *hOrig, const double *newBins, int nNewBins, const char *newName = "");

// Helper function to calculate R-value, Barlow's Sigma, and relative difference |R - 1|
double CalculateRValueAndBarlowSigma(double Y_var_mult, double err_var_mult,
                                     double Y_def_mult, double err_def_mult,
                                     double Y_var_MB, double err_var_MB,
                                     double Y_def_MB, double err_def_MB,
                                     bool &passBarlow)
{
    passBarlow = false;
    if (Y_def_mult <= 0 || Y_def_MB <= 0 || Y_var_MB <= 0 || Y_var_mult <= 0)
    {
        return 0.0;
    }

    double A_mult = Y_var_mult / Y_def_mult;
    double A_MB = Y_var_MB / Y_def_MB;
    double R = A_mult / A_MB;

    // Relative error of A using Barlow's subtraction prescription
    double relA_mult_sq = std::abs(std::pow(err_var_mult / Y_var_mult, 2) - std::pow(err_def_mult / Y_def_mult, 2));
    double relA_MB_sq = std::abs(std::pow(err_var_MB / Y_var_MB, 2) - std::pow(err_def_MB / Y_def_MB, 2));

    if ((relA_mult_sq + relA_MB_sq) < 0)
    {
        cout << "Error: Negative value under square root for Barlow's sigma calculation. Setting sigma_R to 0." << std::endl;
    }

    double sigma_R = R * std::sqrt(relA_mult_sq + relA_MB_sq);

    // Test Barlow's condition: |R - 1| > sigma_R
    double diff = std::abs(R - 1.0);
    if (diff > sigma_R)
    {
        passBarlow = true;
        return diff; // Relative uncertainty ratio |R - 1|
    }

    return 0.0;
}

// Helper function to compute RMS of variations passing Barlow's check within a category
double GetRMSOfVariations(const std::vector<double> &uncertList)
{
    if (uncertList.empty())
        return 0.0;
    double sumSq = 0.0;
    for (double val : uncertList)
    {
        sumSq += val * val;
    }
    return std::sqrt(sumSq / uncertList.size());
}

void uncorrelatedSystematics()
{
    int lineColors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

    string basePathSigExt = "../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/";
    string basePathSigExtpol2 = "../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/";
    string basePathCommon = "../output/kstar/LHC22o_pass7/";
    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    int nmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    TFile *SysUncertainties = new TFile((basePathSigExt + "SystematicsPlots/UnCorrSystematics.root").c_str(), "RECREATE");
    double newpTbins[] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0};
    int nNewpTbins = sizeof(newpTbins) / sizeof(newpTbins[0]) - 1;

    // Load 0-100% Minimum Bias (MB) Default spectrum
    TFile *fMB_Default = nullptr;
    TH1D *hMB_DefaultSpectra = nullptr;
    openTFile(fMB_Default, basePathSigExt + "ROTATED/corrected_spectra_0_100.root");
    openTH1D(hMB_DefaultSpectra, fMB_Default, "mult_0-100/corrected_spectra_Integral_final");
    TH1D *hDefMB_rebinned = RebinHistogram(hMB_DefaultSpectra, newpTbins, nNewpTbins, "hDefMB_rebinned");

    vector<string> normVars = {"Norm1", "Norm2"};
    vector<string> fitRangeVars = {"FitRange1", "FitRange2"};
    vector<string> CombinatorialBkgVars = {"LIKE"};
    vector<string> ResidualBkgVars = {"pol2"};
    vector<string> BinCounting = {"BinCounting"};
    vector<string> widthVars = {"WidthFree"};

    // Load MB variation histograms
    std::vector<TH1D *> hNormMB(normVars.size()), hFitMB(fitRangeVars.size()),
        hRotMB(CombinatorialBkgVars.size()), hPolMB(ResidualBkgVars.size()),
        hBCMB(BinCounting.size()), hWidthMB(widthVars.size());

    for (size_t i = 0; i < normVars.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExt + normVars[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
        hNormMB[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hNormMB_%zu", i));
    }
    for (size_t i = 0; i < fitRangeVars.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExt + fitRangeVars[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
        hFitMB[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hFitMB_%zu", i));
    }
    for (size_t i = 0; i < CombinatorialBkgVars.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExt + CombinatorialBkgVars[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
        hRotMB[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hRotMB_%zu", i));
    }
    for (size_t i = 0; i < ResidualBkgVars.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExtpol2 + ResidualBkgVars[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
        hPolMB[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hPolMB_%zu", i));
    }
    for (size_t i = 0; i < BinCounting.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExt + BinCounting[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_BinCount_final");
        hBCMB[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hBCMB_%zu", i));
    }
    for (size_t i = 0; i < widthVars.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExt + widthVars[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
        hWidthMB[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hWidthMB_%zu", i));
    }

    for (int imult = 0; imult < nmultbins; imult++)
    {
        int multLow = mult_classes[imult];
        int multHigh = mult_classes[imult + 1];

        TString savePath = (basePathSigExt + Form("SystematicsPlots/Uncorrelated/mult_%d-%d/", multLow, multHigh)).c_str();
        if (gSystem->mkdir(savePath, kTRUE))
        {
            std::cout << "Folder " << savePath << " created successfully." << std::endl;
        }

        string multRangeUnderscore = to_string(multLow) + "_" + to_string(multHigh);
        string multRangeDash = to_string(multLow) + "-" + to_string(multHigh);
        string correctedFileName = string("corrected_spectra_") + multRangeUnderscore + ".root";
        string correctedFileNameDefault = string("ROTATED/corrected_spectra_") + multRangeUnderscore + ".root";
        string multDir = string("mult_") + multRangeDash + "/";

        // Load multiplicity-specific variations
        std::vector<TH1D *> hNormMult(normVars.size()), hFitMult(fitRangeVars.size()),
            hCombBkgMult(CombinatorialBkgVars.size()), hPolMult(ResidualBkgVars.size()),
            hBCMult(BinCounting.size()), hWidthMult(widthVars.size());

        TFile *fDefault;
        TH1D *hSpectraDefault;
        openTFile(fDefault, basePathSigExt + correctedFileNameDefault);
        openTH1D(hSpectraDefault, fDefault, (multDir + "corrected_spectra_Integral_final").c_str());
        TH1D *hDefMult_rebinned = RebinHistogram(hSpectraDefault, newpTbins, nNewpTbins, "hDefMult_rebinned");

        for (size_t i = 0; i < normVars.size(); i++)
        {
            TFile *f;
            openTFile(f, basePathSigExt + normVars[i] + "/" + correctedFileName);
            TH1D *h;
            openTH1D(h, f, (multDir + "corrected_spectra_Integral_final").c_str());
            hNormMult[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hNormMult_%zu", i));
        }
        for (size_t i = 0; i < fitRangeVars.size(); i++)
        {
            TFile *f;
            openTFile(f, basePathSigExt + fitRangeVars[i] + "/" + correctedFileName);
            TH1D *h;
            openTH1D(h, f, (multDir + "corrected_spectra_Integral_final").c_str());
            hFitMult[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hFitMult_%zu", i));
        }
        for (size_t i = 0; i < CombinatorialBkgVars.size(); i++)
        {
            TFile *f;
            openTFile(f, basePathSigExt + CombinatorialBkgVars[i] + "/" + correctedFileName);
            TH1D *h;
            openTH1D(h, f, (multDir + "corrected_spectra_Integral_final").c_str());
            hCombBkgMult[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hCombBkgMult_%zu", i));
        }
        for (size_t i = 0; i < ResidualBkgVars.size(); i++)
        {
            TFile *f;
            openTFile(f, basePathSigExtpol2 + ResidualBkgVars[i] + "/" + correctedFileName);
            TH1D *h;
            openTH1D(h, f, (multDir + "corrected_spectra_Integral_final").c_str());
            hPolMult[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hPolMult_%zu", i));
        }
        for (size_t i = 0; i < BinCounting.size(); i++)
        {
            TFile *f;
            openTFile(f, basePathSigExt + BinCounting[i] + "/" + correctedFileName);
            TH1D *h;
            openTH1D(h, f, (multDir + "corrected_spectra_BinCount_final").c_str());
            hBCMult[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hBCMult_%zu", i));
        }
        for (size_t i = 0; i < widthVars.size(); i++)
        {
            TFile *f;
            openTFile(f, basePathSigExt + widthVars[i] + "/" + correctedFileName);
            TH1D *h;
            openTH1D(h, f, (multDir + "corrected_spectra_Integral_final").c_str());
            hWidthMult[i] = RebinHistogram(h, newpTbins, nNewpTbins, Form("hWidthMult_%zu", i));
        }

        // Store relative uncorrelated uncertainty calculated per merged bin
        std::vector<double> relUncorrInMergedBin(nNewpTbins + 1, 0.0);

        // 1. Calculate relative uncorrelated uncertainties in merged pT bins
        for (int ipt = 1; ipt <= nNewpTbins; ipt++)
        {
            double Y_def_mult = hDefMult_rebinned->GetBinContent(ipt);
            double err_def_mult = hDefMult_rebinned->GetBinError(ipt);
            double Y_def_MB = hDefMB_rebinned->GetBinContent(ipt);
            double err_def_MB = hDefMB_rebinned->GetBinError(ipt);

            auto evaluateCategory = [&](const std::vector<TH1D *> &hVarsMult, const std::vector<TH1D *> &hVarsMB)
            {
                std::vector<double> passedUncert;
                for (size_t i = 0; i < hVarsMult.size(); i++)
                {
                    bool pass = false;
                    double uncert = CalculateRValueAndBarlowSigma(
                        hVarsMult[i]->GetBinContent(ipt), hVarsMult[i]->GetBinError(ipt),
                        Y_def_mult, err_def_mult,
                        hVarsMB[i]->GetBinContent(ipt), hVarsMB[i]->GetBinError(ipt),
                        Y_def_MB, err_def_MB, pass);
                    if (pass)
                        passedUncert.push_back(uncert);
                }
                return GetRMSOfVariations(passedUncert);
            };

            // Calculate grouped RMS for each category
            double rmsNorm = evaluateCategory(hNormMult, hNormMB);
            double rmsFit = evaluateCategory(hFitMult, hFitMB);
            // double rmsCombBkg = evaluateCategory(hCombBkgMult, hRotMB);
            double rmsPol = evaluateCategory(hPolMult, hPolMB);
            double rmsBC = evaluateCategory(hBCMult, hBCMB);
            double rmsWidth = evaluateCategory(hWidthMult, hWidthMB);

            // Total relative uncorrelated uncertainty in merged pT bin ipt
            // relUncorrInMergedBin[ipt] = std::sqrt(rmsNorm * rmsNorm + rmsFit * rmsFit +
            //                                       rmsCombBkg * rmsCombBkg + rmsPol * rmsPol +
            //                                       rmsBC * rmsBC + rmsWidth * rmsWidth);

            relUncorrInMergedBin[ipt] = std::sqrt(rmsNorm * rmsNorm + rmsFit * rmsFit + rmsPol * rmsPol + rmsBC * rmsBC + rmsWidth * rmsWidth);
        }

        // 2. Propagate back to original fine pT bins
        TH1D *hUncorrOriginal = (TH1D *)hSpectraDefault->Clone(Form("hUncorrelatedUncertainty_mult_%d_%d", multLow, multHigh));
        hUncorrOriginal->Reset();
        hUncorrOriginal->SetTitle("Uncorrelated Systematic Uncertainty; p_{T} (GeV/c); Uncorrelated Error");

        for (int j = 1; j <= hSpectraDefault->GetNbinsX(); j++)
        {
            double ptCenter = hSpectraDefault->GetBinCenter(j);
            int mergedBinIdx = hDefMult_rebinned->FindBin(ptCenter);

            if (mergedBinIdx >= 1 && mergedBinIdx <= nNewpTbins)
            {
                double relRatio = relUncorrInMergedBin[mergedBinIdx];
                double origYield = hSpectraDefault->GetBinContent(j);

                // Multiply relative ratio by original pT bin yield to obtain absolute uncorrelated error
                double absUncorrelatedError = relRatio * origYield;

                hUncorrOriginal->SetBinContent(j, absUncorrelatedError);
                hUncorrOriginal->SetBinError(j, 0.0);
            }
        }

        SysUncertainties->cd();
        hUncorrOriginal->Write();
    }
}

void openTFile(TFile *&file, const string &path)
{
    file = new TFile(path.c_str(), "READ");
    if (file->IsZombie())
    {
        cout << "Error opening file: " << path << endl;
        file = nullptr;
    }
}

void openTH1D(TH1D *&hist, TFile *file, const string &histPath)
{
    if (file != nullptr)
    {
        hist = (TH1D *)file->Get(histPath.c_str());
        if (hist == nullptr)
        {
            cout << "Error reading histogram: " << histPath << " from file: " << file->GetName() << endl;
        }
    }
    else
    {
        cout << "File is not open. Cannot read histogram: " << histPath << endl;
    }
}

TH1D *RebinHistogram(TH1 *hOrig, const double *newBins, int nNewBins, const char *newName = "")
{
    if (!hOrig)
    {
        std::cerr << "Error: Null histogram pointer passed to RebinHistogram!" << std::endl;
        return nullptr;
    }

    TString name = newName;
    if (name.IsNull())
    {
        name = Form("%s_rebinned", hOrig->GetName());
    }

    // For nNewBins, give the number of new bins, which is one less than the number of edges in newBins array.
    TH1D *hRebinned = (TH1D *)hOrig->Rebin(nNewBins, name.Data(), newBins);

    return hRebinned;
}