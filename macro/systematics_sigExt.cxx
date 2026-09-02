#include <iostream>
#include "src/style.h"
#include "src/initializations.h"
#include "SystematicHelper/helper.cxx"

//// Variations
// Train: 679906 (Default, FT0C, FV0A, TPC1p5_combined2, TPC2p5_combined3p5)
// Train2: 682963 (Default, DCAvar1, DCAvar2, NoPVContributor)

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

void systematics_sigExt()
{
    // 1.         Default                         (1)
    // 2.         Norm variation                  (2)
    // 3.         Fit variation                   (2)
    // 4.         Like-sign bkg                   (1)
    // 5.         pol2                            (1)
    // 6.         Bin counting                    (1)
    // 7.         WidthFree                       (1) //Upto here only we will calculate for all multiplicities
    // 8.         DCA variations                  (2)
    // 9.         PV contributor                  (1)
    // 10.        PID variations                  (2) // Currently using 1
    // 11.        Material Budget                 (2)
    // 12.        Hadronic cross section          (1)

    int lineColors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

    string basePathSigExt = "../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/";
    string basePathSigExtpol2 = "../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/";
    string basePathCommon = "../output/kstar/LHC22o_pass7/";
    string pathPIDAndMultEst = "679906/";
    string pathTrackSel = "679906/";
    // string pathTrackSel = "682963/";
    string basePathTrackSel = basePathCommon + pathTrackSel + "kstarqa_";
    string basePathPIDAndMultEst = basePathCommon + pathPIDAndMultEst + "kstarqa_";

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    // float mult_classes[] = {0, 100.0};
    int nmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1; // number of multiplicity bins

    TFile *fDefault0100;
    TH1D *hSpectraDefault0100;
    TFile *SysUncertainties = new TFile((basePathSigExt + "SystematicsPlots/SysUncert.root").c_str(), "RECREATE");

    ////For signal extraction variations
    vector<string> normVars = {"Norm1", "Norm2"};
    vector<string> fitRangeVars = {"FitRange1", "FitRange2"};
    vector<string> CombinatorialBkgVars = {"LIKE"};
    vector<string> ResidualBkgVars = {"pol2"};
    vector<string> BinCounting = {"BinCounting"};
    vector<string> widthVars = {"WidthFree"};
    int totalSizeSigExt = normVars.size() + fitRangeVars.size() + CombinatorialBkgVars.size() + ResidualBkgVars.size() + BinCounting.size() + widthVars.size();

    // -------------------------------------------------------------
    // Pre-calculate Combinatorial Bkg Uncertainty for Minimum Bias (0-100%)
    // -------------------------------------------------------------
    TFile *fMB_Default = nullptr;
    TH1D *hMB_DefaultSpectra = nullptr;
    openTFile(fMB_Default, basePathSigExt + "ROTATED/corrected_spectra_0_100.root");
    openTH1D(hMB_DefaultSpectra, fMB_Default, "mult_0-100/corrected_spectra_Integral_final");

    vector<TH1 *> MB_CombinatorialBkgVariationHists;
    for (size_t i = 0; i < CombinatorialBkgVars.size(); i++)
    {
        TFile *fMB_Comb = nullptr;
        TH1D *hMB_CombSpectra = nullptr;
        openTFile(fMB_Comb, basePathSigExt + CombinatorialBkgVars[i] + "/corrected_spectra_0_100.root");
        openTH1D(hMB_CombSpectra, fMB_Comb, "mult_0-100/corrected_spectra_Integral_final");
        MB_CombinatorialBkgVariationHists.push_back(hMB_CombSpectra);
    }

    HistogramOperations operationsInitial;

    // Calculate the fixed Minimum Bias relative uncertainty for Combinatorial Background
    TH1D *hRelUncertCombinatorialBkgVars_MB = operationsInitial.RelativeUncertainty(hMB_DefaultSpectra, MB_CombinatorialBkgVariationHists);
    // -------------------------------------------------------------
    // -------------------------------------------------------------

    for (int imult = 0; imult < nmultbins + 1; imult++)
    {
        int multLow, multHigh;

        if (imult == 0)
        {
            multLow = 0;
            multHigh = 100;
        }
        else
        {
            multLow = mult_classes[imult - 1];
            multHigh = mult_classes[imult];
        }

        TString savePath = (basePathSigExt + Form("SystematicsPlots/mult_%d-%d/", multLow, multHigh)).c_str();
        if (gSystem->mkdir(savePath, kTRUE))
        {
            std::cout << "Folder " << savePath << " created successfully." << std::endl;
        }

        // Build mult-range strings used in filenames and histogram paths
        string multRangeUnderscore = to_string(multLow) + "_" + to_string(multHigh);
        string multRangeDash = to_string(multLow) + "-" + to_string(multHigh);
        string correctedFileNameDefault = string("ROTATED/corrected_spectra_") + multRangeUnderscore + ".root";
        string correctedFileName = string("corrected_spectra_") + multRangeUnderscore + ".root";
        string multDir = string("mult_") + multRangeDash + "/";
        cout << "Multiplicity directory: " << multDir << endl;

        TFile *fDefault;
        TFile *fnormVars[normVars.size()];
        TFile *ffitRangeVars[fitRangeVars.size()];
        TFile *fCombinatorialBkgVars[CombinatorialBkgVars.size()];
        TFile *fResidualBkgVars[ResidualBkgVars.size()];
        TFile *fBinCounting[BinCounting.size()];
        TFile *fwidthVars[widthVars.size()];

        TH1D *hSpectraDefault, *hEfficiency_Default;
        TH1D *hSpectraNormVars[normVars.size()];
        TH1D *hSpectraFitRangeVars[fitRangeVars.size()];
        TH1D *hSpectraCombinatorialBkgVars[CombinatorialBkgVars.size()];
        TH1D *hSpectraResidualBkgVars[ResidualBkgVars.size()];
        TH1D *hSpectraBinCounting[BinCounting.size()];
        TH1D *hSpectraWidthVars[widthVars.size()];

        openTFile(fDefault, basePathSigExt + correctedFileNameDefault);
        openTH1D(hSpectraDefault, fDefault, (multDir + "corrected_spectra_Integral_final").c_str());
        openTH1D(hEfficiency_Default, fDefault, (multDir + "heff").c_str());
        openTFile(fDefault0100, basePathSigExt + "ROTATED/corrected_spectra_0_100.root");
        openTH1D(hSpectraDefault0100, fDefault0100, "mult_0-100/corrected_spectra_Integral_final");

        // Signal extraction
        for (int i = 0; i < normVars.size(); i++)
        {
            openTFile(fnormVars[i], basePathSigExt + normVars[i] + "/" + correctedFileName);
            openTH1D(hSpectraNormVars[i], fnormVars[i], (multDir + "corrected_spectra_Integral_final").c_str());
        }
        for (int i = 0; i < fitRangeVars.size(); i++)
        {
            openTFile(ffitRangeVars[i], basePathSigExt + fitRangeVars[i] + "/" + correctedFileName);
            openTH1D(hSpectraFitRangeVars[i], ffitRangeVars[i], (multDir + "corrected_spectra_Integral_final").c_str());
        }
        for (int i = 0; i < CombinatorialBkgVars.size(); i++)
        {
            openTFile(fCombinatorialBkgVars[i], basePathSigExt + CombinatorialBkgVars[i] + "/" + correctedFileName);
            openTH1D(hSpectraCombinatorialBkgVars[i], fCombinatorialBkgVars[i], (multDir + "corrected_spectra_Integral_final").c_str());
        }
        for (int i = 0; i < ResidualBkgVars.size(); i++)
        {
            openTFile(fResidualBkgVars[i], basePathSigExtpol2 + ResidualBkgVars[i] + "/" + correctedFileName);
            openTH1D(hSpectraResidualBkgVars[i], fResidualBkgVars[i], (multDir + "corrected_spectra_Integral_final").c_str());
        }
        for (int i = 0; i < BinCounting.size(); i++)
        {
            openTFile(fBinCounting[i], basePathSigExt + BinCounting[i] + "/" + correctedFileName);
            openTH1D(hSpectraBinCounting[i], fBinCounting[i], (multDir + "corrected_spectra_BinCount_final").c_str());
        }
        for (int i = 0; i < widthVars.size(); i++)
        {
            openTFile(fwidthVars[i], basePathSigExt + widthVars[i] + "/" + correctedFileName);
            openTH1D(hSpectraWidthVars[i], fwidthVars[i], (multDir + "corrected_spectra_Integral_final").c_str());
        }

        vector<TH1 *> NormVariationHists;
        for (int i = 0; i < normVars.size(); i++)
        {
            NormVariationHists.push_back(hSpectraNormVars[i]);
        }
        vector<TH1 *> FitRangeVariationHists;
        for (int i = 0; i < fitRangeVars.size(); i++)
        {
            FitRangeVariationHists.push_back(hSpectraFitRangeVars[i]);
        }
        vector<TH1 *> CombinatorialBkgVariationHists;
        for (int i = 0; i < CombinatorialBkgVars.size(); i++)
        {
            CombinatorialBkgVariationHists.push_back(hSpectraCombinatorialBkgVars[i]);
        }
        vector<TH1 *> ResidualBkgVariationHists;
        for (int i = 0; i < ResidualBkgVars.size(); i++)
        {
            ResidualBkgVariationHists.push_back(hSpectraResidualBkgVars[i]);
        }
        vector<TH1 *> BinCountingVariationHists;
        for (int i = 0; i < BinCounting.size(); i++)
        {
            BinCountingVariationHists.push_back(hSpectraBinCounting[i]);
        }
        vector<TH1 *> WidthVariationHists;
        for (int i = 0; i < widthVars.size(); i++)
        {
            WidthVariationHists.push_back(hSpectraWidthVars[i]);
        }

        vector<size_t> AllVarSizes = {
            NormVariationHists.size(),
            FitRangeVariationHists.size(),
            CombinatorialBkgVariationHists.size(),
            ResidualBkgVariationHists.size(),
            BinCountingVariationHists.size(),
            WidthVariationHists.size()};
        int totalCombinations = AllVarSizes.size();
        cout << "Total number of systematic sources are " << totalCombinations << endl;

        int AllVariationsSize = 0;
        for (size_t i = 0; i < AllVarSizes.size(); i++)
        {
            AllVariationsSize += AllVarSizes[i];
        }
        cout << "Total number of variations are " << AllVariationsSize << endl;

        vector<vector<TH1 *>> AllVariationHists = {
            NormVariationHists,
            FitRangeVariationHists,
            CombinatorialBkgVariationHists,
            ResidualBkgVariationHists,
            BinCountingVariationHists,
            WidthVariationHists};

        vector<vector<string>> AllVariationNames = {
            normVars,
            fitRangeVars,
            CombinatorialBkgVars,
            ResidualBkgVars,
            BinCounting,
            widthVars};

        vector<TH1 *> VarRatios[AllVariationsSize];
        vector<TH1 *> histBarlow[AllVariationsSize];

        HistogramOperations operations;
        bool checkbar = true;
        vector<bool> barlowPassed;
        int counter = 0;

        for (size_t i = 0; i < AllVariationHists.size(); i++)
        {
            for (size_t j = 0; j < AllVariationHists[i].size(); j++)
            {
                // if (i == AllVariationHists.size() - 1)
                // {
                //     VarRatios[i].push_back(operations.CalculateRatio(hSpectraDefault0100, AllVariationHists[i][j]));
                //     histBarlow[i].push_back(operations.barlowcheck(hSpectraDefault0100, AllVariationHists[i][j], checkbar));
                // }
                // else
                {
                    VarRatios[i].push_back(operations.CalculateRatio(hSpectraDefault, AllVariationHists[i][j]));
                    histBarlow[i].push_back(operations.barlowcheck(hSpectraDefault, AllVariationHists[i][j], checkbar));
                }
                cout << "Systematic source " << i + 1 << ", Variation " << j + 1 << ", Barlow check " << (checkbar ? "Passed" : "Failed") << endl;
                barlowPassed.push_back(checkbar);
                counter++;
            }
        }

        TCanvas *cPlotBarlowAll = new TCanvas("", "Barlow Checks", 1080, 720);
        SetCanvasStyle(cPlotBarlowAll, 0.15, 0.03, 0.06, 0.15);
        cPlotBarlowAll->Divide(3, 3);
        TCanvas *cRatioAll = new TCanvas("", "Spectra ratio (Variation / Default)", 1080, 720);
        SetCanvasStyle(cRatioAll, 0.15, 0.03, 0.06, 0.15);
        cRatioAll->Divide(3, 3);
        counter = 0;
        TLatex latBarlow;
        latBarlow.SetNDC();
        latBarlow.SetTextSize(0.075);
        latBarlow.SetTextColor(kBlue);
        latBarlow.SetTextFont(20);
        gStyle->SetTitleFontSize(0.08);
        gStyle->SetTitleFont(62, "t");
        TLine *line1 = new TLine();
        line1->SetLineColor(kRed);
        line1->SetLineStyle(2);
        line1->SetLineWidth(2);
        for (size_t i = 0; i < AllVariationHists.size(); i++)
        {
            for (size_t j = 0; j < AllVariationHists[i].size(); j++)
            {
                cPlotBarlowAll->cd(counter + 1);
                gPad->SetLeftMargin(0.06);
                gPad->SetBottomMargin(0.14);
                gPad->SetRightMargin(0.03);
                gPad->SetTopMargin(0.1);
                SetHistoQA(histBarlow[i][j]);
                histBarlow[i][j]->SetTitle(Form("%s", AllVariationNames[i][j].c_str()));
                histBarlow[i][j]->GetXaxis()->SetTitle("#Delta/#sigma");
                histBarlow[i][j]->GetXaxis()->SetTitleSize(0.06);
                histBarlow[i][j]->SetStats(0);
                histBarlow[i][j]->SetMaximum(histBarlow[i][j]->GetMaximum() * 1.5);
                histBarlow[i][j]->Draw("HIST");
                float mean = histBarlow[i][j]->GetMean();
                float rms = histBarlow[i][j]->GetRMS();
                int countstotal = histBarlow[i][j]->Integral();
                float deltaR1 = histBarlow[i][j]->Integral(histBarlow[i][j]->FindBin(-1), histBarlow[i][j]->FindBin(1)) / countstotal;
                float deltaR2 = histBarlow[i][j]->Integral(histBarlow[i][j]->FindBin(-2), histBarlow[i][j]->FindBin(2)) / countstotal;
                if (barlowPassed[counter])
                    latBarlow.SetTextColor(kRed);
                else
                    latBarlow.SetTextColor(kBlue);
                latBarlow.DrawLatex(0.1, 0.82, Form("%s", barlowPassed[counter] ? "Passed" : "Failed"));
                latBarlow.DrawLatex(0.1, 0.72, Form("Mean: %.3f", mean));
                latBarlow.DrawLatex(0.1, 0.62, Form("RMS: %.3f", rms));
                latBarlow.DrawLatex(0.1, 0.52, Form("|n|<1: %.2f%%", deltaR1 * 100));
                latBarlow.DrawLatex(0.1, 0.42, Form("|n|<2: %.2f%%", deltaR2 * 100));

                cRatioAll->cd(counter + 1);
                gPad->SetGridx();
                gPad->SetLeftMargin(0.15);
                gPad->SetBottomMargin(0.14);
                gPad->SetRightMargin(0.03);
                gPad->SetTopMargin(0.1);
                SetHistoQA(VarRatios[i][j]);
                VarRatios[i][j]->GetXaxis()->SetLabelSize(0.06);
                VarRatios[i][j]->GetYaxis()->SetLabelSize(0.06);
                VarRatios[i][j]->GetYaxis()->SetTitle("Variation / Default");
                VarRatios[i][j]->GetYaxis()->SetTitleOffset(1.6);
                VarRatios[i][j]->SetTitle(Form("%s", AllVariationNames[i][j].c_str()));
                VarRatios[i][j]->SetStats(0);
                VarRatios[i][j]->SetMaximum(VarRatios[i][j]->GetMaximum() * 1.02);
                VarRatios[i][j]->Draw("HIST");

                line1->SetX1(0);
                line1->SetX2(VarRatios[i][j]->GetXaxis()->GetXmax());
                line1->SetY1(1);
                line1->SetY2(1);
                line1->Draw("SAME");
                counter++;
            }
        }
        TH1D *hRelUncertNormVars = operations.RelativeUncertainty(hSpectraDefault, NormVariationHists);
        TH1D *hRelUncertFitRangeVars = operations.RelativeUncertainty(hSpectraDefault, FitRangeVariationHists);
        // TH1D *hRelUncertCombinatorialBkgVars = operations.RelativeUncertainty(hSpectraDefault, CombinatorialBkgVariationHists);
        TH1D *hRelUncertCombinatorialBkgVars = (TH1D *)hRelUncertCombinatorialBkgVars_MB->Clone(Form("hRelUncertCombinatorialBkgVars_mult_%d_%d", multLow, multHigh));
        TH1D *hRelUncertResidualBkgVars = operations.RelativeUncertainty(hSpectraDefault, ResidualBkgVariationHists);
        TH1D *hRelUncertBinCounting = operations.RelativeUncertainty(hSpectraDefault, BinCountingVariationHists);
        TH1D *hRelUncertWidthVars = operations.RelativeUncertainty(hSpectraDefault, WidthVariationHists);

        TCanvas *cRelUncert = new TCanvas("", "Relative Uncertainties", 1080, 720);
        SetCanvasStyle(cRelUncert, 0.15, 0.03, 0.06, 0.15);
        cRelUncert->Divide(3, 2);
        vector<TH1D *> relUncertHists = {
            (TH1D *)hRelUncertNormVars->Clone("hRelUncertNormVars_clone"),
            (TH1D *)hRelUncertFitRangeVars->Clone("hRelUncertFitRangeVars_clone"),
            (TH1D *)hRelUncertCombinatorialBkgVars->Clone("hRelUncertCombinatorialBkgVars_clone"),
            (TH1D *)hRelUncertResidualBkgVars->Clone("hRelUncertResidualBkgVars_clone"),
            (TH1D *)hRelUncertBinCounting->Clone("hRelUncertBinCounting_clone"),
            (TH1D *)hRelUncertWidthVars->Clone("hRelUncertWidthVars_clone")};

        vector<string> relUncertNames = {
            "Norm. Range",
            "Fit Range",
            "Combinatorial Bkg",
            "Residual Bkg",
            "Yield Extraction",
            "Width Variation"};

        for (size_t i = 0; i < relUncertHists.size(); i++)
        {
            cRelUncert->cd(i + 1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.03);
            gPad->SetTopMargin(0.15);
            SetHistoQA(relUncertHists[i]);
            relUncertHists[i]->GetXaxis()->SetLabelSize(0.06);
            relUncertHists[i]->GetYaxis()->SetLabelSize(0.06);
            // relUncertHists[i]->GetYaxis()->SetTitle("Relative Uncertainty");
            relUncertHists[i]->GetYaxis()->SetTitle("");
            relUncertHists[i]->SetTitle(Form("%s", relUncertNames[i].c_str()));
            relUncertHists[i]->SetStats(0);
            // relUncertHists[i]->GetYaxis()->SetMaxDigits(3);
            relUncertHists[i]->SetMaximum(relUncertHists[i]->GetMaximum() * 1.05);
            relUncertHists[i]->Draw("HIST");
        }

        vector<TH1D *> vecSignalExt = {hRelUncertNormVars, hRelUncertFitRangeVars, hRelUncertCombinatorialBkgVars, hRelUncertResidualBkgVars, hRelUncertWidthVars}; // Bin counting excluded as it passed Barlow

        TH1D *hSignalExtTotalSys = operations.sigma(vecSignalExt);
        TH1D *hSignalExtTotalSysClone = (TH1D *)hSignalExtTotalSys->Clone();

        string SigExtNames[] = {"Norm. range", "Fit Range", "Combinatorial Bkg", "Residual Bkg", "Width fix/free"};

        TCanvas *cSigExtAll = new TCanvas("", "Systematic Uncertainties from all sources", 1080, 720);
        SetCanvasStyle(cSigExtAll, 0.14, 0.03, 0.06, 0.13);
        TLegend *legSigExt = new TLegend(0.17, 0.6, 0.5, 0.88);
        legSigExt->SetBorderSize(0);
        legSigExt->SetFillStyle(0);
        legSigExt->SetTextSize(0.03);
        legSigExt->SetTextFont(42);
        legSigExt->SetHeader("Signal Extraction");
        legSigExt->AddEntry((TObject *)0, Form("Multiplicity: %d-%d", multLow, multHigh), "");
        for (int i = 0; i < vecSignalExt.size(); i++)
        {
            SetHistoQA(vecSignalExt[i]);
            vecSignalExt[i]->GetYaxis()->SetTitle("Relative Uncertainty");
            vecSignalExt[i]->SetStats(0);
            vecSignalExt[i]->SetMaximum(0.21);
            vecSignalExt[i]->SetMinimum(0);
            vecSignalExt[i]->SetLineColor(lineColors[i]);
            vecSignalExt[i]->Draw("HIST SAME");
            legSigExt->AddEntry(vecSignalExt[i], Form("%s", SigExtNames[i].c_str()), "l");
        }
        hSignalExtTotalSysClone->SetLineColor(1);
        hSignalExtTotalSysClone->SetLineWidth(3);
        hSignalExtTotalSysClone->Draw("HIST SAME");
        legSigExt->AddEntry(hSignalExtTotalSysClone, "Total", "l");
        legSigExt->Draw();

        // smoothing procedure on separate sources
        int Iterations = 2;
        TH1D *hSignalExtTotalSysSmoothed = operations.smooth(hSignalExtTotalSys, Iterations);

        // Save all the plots
        cPlotBarlowAll->SaveAs(savePath + "BarlowChecks_AllVariations.png");
        cRatioAll->SaveAs(savePath + "Ratio_AllVariations.png");
        cRelUncert->SaveAs(savePath + "RelativeUncertainties_AllSources.png");
        cSigExtAll->SaveAs(savePath + "SignalExtractionSystematics.png");
    }
}