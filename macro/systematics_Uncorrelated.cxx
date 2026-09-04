#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include "src/style.h"
#include "src/initializations.h"
#include "SystematicHelper/helper.cxx"

void openTFile(TFile *&file, const string &path);
void openTH1D(TH1D *&hist, TFile *file, const string &histPath);

TH1D *smooth(TH1D *hist1, int n = 2)
{
    TH1D *hsmooth = (TH1D *)hist1->Clone();
    for (int i = 0; i < n; i++)
    {
        for (int j = 1; j < hsmooth->GetNbinsX() - 1; j++)
        {
            double bin1 = hsmooth->GetBinContent(j);
            double bin2 = hsmooth->GetBinContent(j + 1);
            double bin3 = hsmooth->GetBinContent(j + 2);
            double avg = (bin1 + bin2 + bin3) / 3;
            hsmooth->SetBinContent(j + 1, avg);
        }
    }
    return hsmooth;
}

// ---------------------------------------------------------------------
// Manual Rebinning matching uncorrSys.cpp logic (density <-> yield)
// ---------------------------------------------------------------------
TH1D *RebinHistogram(TH1 *h, const TString &name, const std::vector<double> &wideBins)
{
    if (!h)
        return nullptr;

    TH1D *hOut = new TH1D(name, h->GetTitle(), wideBins.size() - 1, wideBins.data());
    hOut->SetDirectory(0);

    for (int i = 1; i <= hOut->GetNbinsX(); ++i)
    {
        double ptLow = hOut->GetXaxis()->GetBinLowEdge(i);
        double ptHigh = hOut->GetXaxis()->GetBinUpEdge(i);
        double wideWidth = ptHigh - ptLow;

        double sumYield = 0.0;
        double sumErrSq = 0.0;

        for (int j = 1; j <= h->GetNbinsX(); ++j)
        {
            double fineCenter = h->GetXaxis()->GetBinCenter(j);
            if (fineCenter >= ptLow && fineCenter < ptHigh)
            {
                double fineWidth = h->GetXaxis()->GetBinWidth(j);

                // Convert yield density back to absolute yield before summing
                double yield = h->GetBinContent(j) * fineWidth;
                double err = h->GetBinError(j) * fineWidth;

                sumYield += yield;
                sumErrSq += err * err;
            }
        }

        // Convert back to density for the wide bin
        hOut->SetBinContent(i, sumYield / wideWidth);
        hOut->SetBinError(i, std::sqrt(sumErrSq) / wideWidth);
    }
    return hOut;
}

// ---------------------------------------------------------------------
// R-factor & Barlow check matching uncorrSys.cpp logic
// ---------------------------------------------------------------------
double CalculateRValueAndBarlowSigma(double yVarMult, double eVarMult, double yDefMult, double eDefMult,
                                     double yVarMB, double eVarMB, double yDefMB, double eDefMB,
                                     double &R, double &sigmaRB)
{
    R = 1.0;
    sigmaRB = 0.0;
    if (yVarMult <= 0 || yDefMult <= 0 || yVarMB <= 0 || yDefMB <= 0)
        return 0.0;

    double relVarMult = eVarMult / yVarMult;
    double relDefMult = eDefMult / yDefMult;
    double relVarMB = eVarMB / yVarMB;
    double relDefMB = eDefMB / yDefMB;

    double relA_mult = std::sqrt(std::abs(relVarMult * relVarMult - relDefMult * relDefMult));
    double relA_MB = std::sqrt(std::abs(relVarMB * relVarMB - relDefMB * relDefMB));
    double relR = std::sqrt(relA_mult * relA_mult + relA_MB * relA_MB);

    double Amult = yVarMult / yDefMult;
    double AMB = yVarMB / yDefMB;
    R = Amult / AMB;
    sigmaRB = std::abs(R) * relR;

    if (std::isnan(R) || std::isnan(sigmaRB) || std::isinf(sigmaRB))
    {
        R = 1.0;
        sigmaRB = 0.0;
        return -999.0;
    }

    double delta = pow(std::abs(R - 1.0), 1);

    return (delta > sigmaRB) ? delta : -999.0;
}

void systematics_Uncorrelated()
{
    int lineColors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

    TFile *fTotalSys = new TFile("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots/SysUncert.root", "READ");
    if (fTotalSys->IsZombie())
    {
        cout << "Error: SysUncert.root file not found" << endl;
        return;
    }

    string basePathSigExt = "../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/";
    string basePathSigExtpol2 = "../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/";
    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    int nmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    TFile *SysUncertainties = new TFile((basePathSigExt + "SystematicsPlots/UnCorrSystematics.root").c_str(), "RECREATE");
    // std::vector<double> widePtBins = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0};
    std::vector<double> widePtBins = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0};
    int nWideBins = widePtBins.size() - 1;

    // Load 0-100% Minimum Bias (MB) Default spectrum
    TFile *fMB_Default = nullptr;
    TH1D *hMB_DefaultSpectra = nullptr;
    openTFile(fMB_Default, basePathSigExt + "ROTATED/corrected_spectra_0_100.root");
    openTH1D(hMB_DefaultSpectra, fMB_Default, "mult_0-100/corrected_spectra_Integral_final");
    TH1D *hDefMB_rebinned = RebinHistogram(hMB_DefaultSpectra, "hDefMB_rebinned", widePtBins);

    vector<string> normVars = {"Norm1", "Norm2"};
    vector<string> fitRangeVars = {"FitRange1", "FitRange2"};
    // vector<string> CombinatorialBkgVars = {"LIKE"};
    vector<string> ResidualBkgVars = {"pol2"};
    vector<string> BinCounting = {"BinCounting"};
    vector<string> widthVars = {"WidthFree"};

    // // Load MB variation histograms
    // std::vector<TH1D *> hNormMB(normVars.size()), hFitMB(fitRangeVars.size()),
    //     hLikeMB(CombinatorialBkgVars.size()), hPolMB(ResidualBkgVars.size()),
    //     hBCMB(BinCounting.size()), hWidthMB(widthVars.size());

    std::vector<TH1D *> hNormMB(normVars.size()), hFitMB(fitRangeVars.size()), hPolMB(ResidualBkgVars.size()), hBCMB(BinCounting.size()), hWidthMB(widthVars.size());

    for (size_t i = 0; i < normVars.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExt + normVars[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
        hNormMB[i] = RebinHistogram(h, Form("hNormMB_%zu", i), widePtBins);
    }
    for (size_t i = 0; i < fitRangeVars.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExt + fitRangeVars[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
        hFitMB[i] = RebinHistogram(h, Form("hFitMB_%zu", i), widePtBins);
    }
    // for (size_t i = 0; i < CombinatorialBkgVars.size(); i++)
    // {
    //     TFile *f;
    //     openTFile(f, basePathSigExt + CombinatorialBkgVars[i] + "/corrected_spectra_0_100.root");
    //     TH1D *h;
    //     openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
    //     hLikeMB[i] = RebinHistogram(h, Form("hLikeMB_%zu", i), widePtBins);
    // }
    for (size_t i = 0; i < ResidualBkgVars.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExtpol2 + ResidualBkgVars[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
        hPolMB[i] = RebinHistogram(h, Form("hPolMB_%zu", i), widePtBins);
    }
    for (size_t i = 0; i < BinCounting.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExt + BinCounting[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_BinCount_final");
        hBCMB[i] = RebinHistogram(h, Form("hBCMB_%zu", i), widePtBins);
    }
    for (size_t i = 0; i < widthVars.size(); i++)
    {
        TFile *f;
        openTFile(f, basePathSigExt + widthVars[i] + "/corrected_spectra_0_100.root");
        TH1D *h;
        openTH1D(h, f, "mult_0-100/corrected_spectra_Integral_final");
        hWidthMB[i] = RebinHistogram(h, Form("hWidthMB_%zu", i), widePtBins);
    }

    // Dynamic structures to process sources following uncorrSys.cpp logic
    struct VarGroup
    {
        string name;
        std::vector<TH1D *> hMB;
        std::vector<string> names;
        bool isPol2;
        bool isBC;
    };
    std::vector<VarGroup> sourceGroups = {
        {"Norm", hNormMB, normVars, false, false},
        {"FitRange", hFitMB, fitRangeVars, false, false},
        // {"CombinatorialBkg", hLikeMB, CombinatorialBkgVars, false, false},
        {"ResidualBkg", hPolMB, ResidualBkgVars, true, false},
        {"BinCounting", hBCMB, BinCounting, false, true},
        {"Width", hWidthMB, widthVars, false, false}};

    std::map<string, std::vector<std::vector<double>>> uncorrBySource;

    // Process each systematic source independently across all multiplicity classes
    for (auto &grp : sourceGroups)
    {
        std::vector<std::vector<std::vector<double>>> uncorrByVar;

        for (size_t ivar = 0; ivar < grp.names.size(); ivar++)
        {
            std::vector<std::vector<double>> uncorr_thisVar(nWideBins, std::vector<double>(nmultbins, 0.0));

            for (int imult = 0; imult < nmultbins; imult++)
            {
                int multLow = mult_classes[imult];
                int multHigh = mult_classes[imult + 1];
                string multRangeUnderscore = to_string(multLow) + "_" + to_string(multHigh);
                string multRangeDash = to_string(multLow) + "-" + to_string(multHigh);
                string correctedFileName = string("corrected_spectra_") + multRangeUnderscore + ".root";
                string correctedFileNameDefault = string("ROTATED/corrected_spectra_") + multRangeUnderscore + ".root";
                string multDir = string("mult_") + multRangeDash + "/";

                TFile *fDefault = nullptr;
                TH1D *hSpectraDefault = nullptr;
                openTFile(fDefault, basePathSigExt + correctedFileNameDefault);
                openTH1D(hSpectraDefault, fDefault, (multDir + "corrected_spectra_Integral_final").c_str());
                TH1D *hDefMult = RebinHistogram(hSpectraDefault, "hDefMult", widePtBins);

                TFile *fVar = nullptr;
                TH1D *hVarRaw = nullptr;
                string varPath = (grp.isPol2 ? basePathSigExtpol2 : basePathSigExt) + grp.names[ivar] + "/" + correctedFileName;
                string histNameMult = grp.isBC ? (multDir + "corrected_spectra_BinCount_final") : (multDir + "corrected_spectra_Integral_final");

                openTFile(fVar, varPath);
                openTH1D(hVarRaw, fVar, histNameMult.c_str());
                TH1D *hVarMult = RebinHistogram(hVarRaw, "hVarMult", widePtBins);

                for (int p = 0; p < nWideBins; ++p)
                {
                    double yDef = hDefMult->GetBinContent(p + 1), eDef = hDefMult->GetBinError(p + 1);
                    double yVar = hVarMult->GetBinContent(p + 1), eVar = hVarMult->GetBinError(p + 1);
                    double yDefMB = hDefMB_rebinned->GetBinContent(p + 1), eDefMB = hDefMB_rebinned->GetBinError(p + 1);
                    double yVarMB = grp.hMB[ivar]->GetBinContent(p + 1), eVarMB = grp.hMB[ivar]->GetBinError(p + 1);

                    double R, sigmaRB;
                    uncorr_thisVar[p][imult] = CalculateRValueAndBarlowSigma(yVar, eVar, yDef, eDef, yVarMB, eVarMB, yDefMB, eDefMB, R, sigmaRB);
                }

                if (fDefault)
                {
                    fDefault->Close();
                    delete fDefault;
                }
                if (fVar)
                {
                    fVar->Close();
                    delete fVar;
                }
                delete hDefMult;
                delete hVarMult;
            }
            uncorrByVar.push_back(uncorr_thisVar);
        }

        // Combine variations inside source via RMS
        std::vector<std::vector<double>> uncorr_source(nWideBins, std::vector<double>(nmultbins, 0.0));
        for (int p = 0; p < nWideBins; ++p)
        {
            for (int c = 0; c < nmultbins; ++c)
            {
                double sumSq = 0.0;
                for (auto &uv : uncorrByVar)
                {
                    if (uv[p][c] > 0)
                        sumSq += uv[p][c] * uv[p][c];
                }

                uncorr_source[p][c] = std::sqrt(sumSq / uncorrByVar.size());
            }
        }
        uncorrBySource[grp.name] = uncorr_source;
    }

    // Quadrature-sum across distinct sources per bin
    std::vector<std::vector<double>> combinedFrac(nWideBins, std::vector<double>(nmultbins, 0.0));
    for (int p = 0; p < nWideBins; ++p)
    {
        for (int c = 0; c < nmultbins; ++c)
        {
            double sumSq = 0.0;
            for (auto &kv : uncorrBySource)
                sumSq += kv.second[p][c] * kv.second[p][c];
            combinedFrac[p][c] = std::sqrt(sumSq);
        }
    }

    // Take MAX across multiplicity classes per wide bin (ALICE Note Prescription)
    std::vector<double> finalFracWide(nWideBins, 0.0);
    for (int p = 0; p < nWideBins; ++p)
    {
        double maxVal = 0.0;
        for (int c = 0; c < nmultbins; ++c)
            maxVal = std::max(maxVal, combinedFrac[p][c]);
        finalFracWide[p] = maxVal;
    }

    // Save and plot output per multiplicity class
    for (int imult = 0; imult < nmultbins; imult++)
    {
        int multLow = mult_classes[imult];
        int multHigh = mult_classes[imult + 1];

        TString savePath = (basePathSigExt + Form("SystematicsPlots/Uncorrelated/mult_%d-%d/", multLow, multHigh)).c_str();
        if (gSystem->mkdir(savePath, kTRUE))
        {
            std::cout << "Folder " << savePath << " created successfully." << std::endl;
        }

        TH1D *hTotalUncert = (TH1D *)fTotalSys->Get(Form("hTotalSysSmoothed_%d_%d", multLow, multHigh));
        TH1D *hSigExtUncert = (TH1D *)fTotalSys->Get(Form("hSignalExtTotalSysSmoothed_%d_%d", multLow, multHigh));
        if (!hTotalUncert || !hSigExtUncert)
        {
            std::cout << "Error: Could not retrieve hTotalSys_" << multLow << "_" << multHigh << " from fTotalSys" << std::endl;
            return;
        }

        string multRangeUnderscore = to_string(multLow) + "_" + to_string(multHigh);
        string multRangeDash = to_string(multLow) + "-" + to_string(multHigh);
        string correctedFileNameDefault = string("ROTATED/corrected_spectra_") + multRangeUnderscore + ".root";
        string multDir = string("mult_") + multRangeDash + "/";

        TFile *fDefault;
        TH1D *hSpectraDefault;
        openTFile(fDefault, basePathSigExt + correctedFileNameDefault);
        openTH1D(hSpectraDefault, fDefault, (multDir + "corrected_spectra_Integral_final").c_str());

        // Rebin fine total systematic to wide bins to determine scaling fraction
        TH1D *hTotalSysWide = RebinHistogram(hTotalUncert, Form("hTotWide_%d_%d", multLow, multHigh), widePtBins);

        TH1D *hUncorrOriginal = (TH1D *)hTotalUncert->Clone(Form("hUncorrelatedUncertainty_mult_%d_%d", multLow, multHigh));
        hUncorrOriginal->Reset();
        hUncorrOriginal->SetTitle("Uncorrelated Systematic Uncertainty; p_{T} (GeV/c); Uncorrelated Error");

        for (int j = 1; j <= hSpectraDefault->GetNbinsX(); j++)
        {
            double ptCenter = hSpectraDefault->GetBinCenter(j);
            int wideBin = -1;
            for (int p = 0; p < nWideBins; ++p)
            {
                if (ptCenter >= widePtBins[p] && ptCenter < widePtBins[p + 1])
                {
                    wideBin = p;
                    break;
                }
            }

            if (wideBin >= 0)
            {
                // double relUncorrWide = finalFracWide[wideBin];
                double relUncorrWide = combinedFrac[wideBin][imult];
                double relTotWide = hTotalSysWide->GetBinContent(wideBin + 1);

                // 1. Calculate true fraction
                double trueFraction = (relTotWide > 0) ? (relUncorrWide / relTotWide) : 0.0;

                // 2. Cap fraction at 1.0 so uncorrelated never exceeds total systematic
                if (trueFraction > 1.0)
                {
                    trueFraction = 1.0;
                }

                // 3. Scale fine-binned relative total uncertainty
                double fineTotalSys = hTotalUncert->GetBinContent(j);
                double finalRelUncorrFine = trueFraction * fineTotalSys;

                // Safety guardrail
                if (finalRelUncorrFine > fineTotalSys)
                {
                    finalRelUncorrFine = fineTotalSys;
                }

                // hUncorrOriginal->SetBinContent(j, relUncorrWide);
                hUncorrOriginal->SetBinContent(j, finalRelUncorrFine);
                hUncorrOriginal->SetBinError(j, 0.0);
            }
        }

        SysUncertainties->cd();
        hUncorrOriginal->Write();

        // Smooth the uncorrelated histogram for better visualization
        int Iterations = 2;
        TH1D *hUncorrSmoothed = smooth(hUncorrOriginal, Iterations);
        hUncorrSmoothed->SetName(Form("hUncorrelatedUncertaintySmoothed_mult_%d_%d", multLow, multHigh));
        hUncorrSmoothed->SetTitle("Smoothed Uncorrelated Systematic Uncertainty; p_{T} (GeV/c); Uncorrelated Error");
        SysUncertainties->cd();
        hUncorrSmoothed->Write();

        // -------------------------------------------------------------
        // Plotting logic (Preserved original visualization)
        // -------------------------------------------------------------
        TCanvas *cUncorr = new TCanvas(Form("cUncorr_mult_%d_%d", multLow, multHigh), "Uncorrelated Systematic Uncertainty", 720, 720);
        SetCanvasStyle(cUncorr, 0.15, 0.03, 0.06, 0.15);
        SetHistoQA(hTotalUncert);
        hTotalUncert->SetLineColor(kBlue);
        hTotalUncert->SetMinimum(0.0);
        hTotalUncert->SetMaximum(0.2);
        hTotalUncert->Draw("HIST");

        SetHistoQA(hUncorrOriginal);
        hUncorrOriginal->SetLineColor(kRed);
        hUncorrOriginal->Draw("HIST SAME");

        hUncorrSmoothed->SetLineColor(kGreen + 2);
        hUncorrSmoothed->Draw("HIST SAME");

        TLegend *legUncorr = new TLegend(0.2, 0.7, 0.85, 0.85);
        legUncorr->SetBorderSize(0);
        legUncorr->SetFillStyle(0);
        legUncorr->SetTextSize(0.027);
        legUncorr->AddEntry(hUncorrOriginal, "Uncorrelated Systematic Uncertainty", "l");
        legUncorr->AddEntry(hTotalUncert, "Total Systematic Uncertainty", "l");
        legUncorr->AddEntry(hUncorrSmoothed, "Smoothed Uncorrelated Systematic Uncertainty", "l");
        legUncorr->Draw();

        cUncorr->SaveAs((basePathSigExt + Form("SystematicsPlots/Uncorrelated/hUncorrelatedUncertainty_mult_%d_%d.png", multLow, multHigh)).c_str());

        if (fDefault)
        {
            fDefault->Close();
            delete fDefault;
        }
        delete hTotalSysWide;
    }

    SysUncertainties->Close();
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