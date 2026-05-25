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

void systematics()
{
    // 1.         Default                         (1)
    // 2.         Norm variation                  (2)
    // 3.         Fit variation                   (2)
    // 4.         Rotatinal bkg                   (1)
    // 5.         pol2                            (1)
    // 6.         Bin counting                    (1)
    // 7.         WidthFree                       (1)
    // 8.         DCA variations                  (2)
    // 9.         PV contributor                  (1)
    // 10.        PID variations                  (2)// Currently using 1
    // 11.        Multiplicity estimator          (2)
    // 12.        Material Budget                 (2)
    // 13.        ITS-TPS matching                (1)
    // 14.        Hadronic cross section          (1)

    int lineColors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

    string basePathSigExt = "../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/";
    string basePathCommon = "../output/kstar/LHC22o_pass7/";
    string pathPIDAndMultEst = "679906/";
    string pathTrackSel = "682963/";
    string basePathTrackSel = basePathCommon + pathTrackSel + "kstarqa_";
    string basePathPIDAndMultEst = basePathCommon + pathPIDAndMultEst + "kstarqa_";
    TString savePath = basePathSigExt + "SystematicsPlots/";
    int multLow = 0;
    int multHigh = 100;

    // Build mult-range strings used in filenames and histogram paths
    string multRangeUnderscore = to_string(multLow) + "_" + to_string(multHigh);
    string multRangeDash = to_string(multLow) + "-" + to_string(multHigh);
    string correctedFileName = string("corrected_spectra_") + multRangeUnderscore + ".root";
    string multDir = string("mult_") + multRangeDash + "/";

    ////For signal extraction variations
    vector<string> normVars = {"Norm1", "Norm2"};
    vector<string> fitRangeVars = {"FitRange1", "FitRange2"};
    vector<string> CombinatorialBkgVars = {"ROTATED"};
    vector<string> ResidualBkgVars = {"pol2"};
    vector<string> BinCounting = {"BinCounting"};
    vector<string> widthVars = {"WidthFree"};
    int totalSizeSigExt = normVars.size() + fitRangeVars.size() + CombinatorialBkgVars.size() + ResidualBkgVars.size() + BinCounting.size() + widthVars.size();

    //// For track selections variations
    vector<string> DCAvars = {"DCAvar1", "DCAvar2"};
    vector<string> PVcontributorVars = {"NoPVContributor"};
    int totalSizeTrackSel = DCAvars.size() + PVcontributorVars.size();

    //// For PID variations
    // vector<string> PIDVars = {"TPC1p5_combined2", "TPC2p5_combined3p5"};
    vector<string> PIDVars = {"TPC2p5_combined3p5"};

    //// Multiplicity estimator variations
    vector<string> multEstVars = {"FT0C", "FV0A"};

    //// Material budget variations
    vector<string> materialBudgetVars = {"MaterialBudgetMinus10", "MaterialBudgetPlus10"};

    TFile *fDefault;
    TFile *fnormVars[normVars.size()];
    TFile *ffitRangeVars[fitRangeVars.size()];
    TFile *fCombinatorialBkgVars[CombinatorialBkgVars.size()];
    TFile *fResidualBkgVars[ResidualBkgVars.size()];
    TFile *fBinCounting[BinCounting.size()];
    TFile *fwidthVars[widthVars.size()];
    TFile *fDCAvars[DCAvars.size()];
    TFile *fPVcontributorVars[PVcontributorVars.size()];
    TFile *fPIDVars[PIDVars.size()];
    TFile *fmultEstVars[multEstVars.size()];
    TFile *fmaterialBudgetVars[materialBudgetVars.size()];

    TH1D *hSpectraDefault, *hEfficiency_Default;
    TH1D *hSpectraNormVars[normVars.size()];
    TH1D *hSpectraFitRangeVars[fitRangeVars.size()];
    TH1D *hSpectraCombinatorialBkgVars[CombinatorialBkgVars.size()];
    TH1D *hSpectraResidualBkgVars[ResidualBkgVars.size()];
    TH1D *hSpectraBinCounting[BinCounting.size()];
    TH1D *hSpectraWidthVars[widthVars.size()];
    TH1D *hSpectraDCAvars[DCAvars.size()];
    TH1D *hSpectraPVcontributorVars[PVcontributorVars.size()];
    TH1D *hSpectraPIDVars[PIDVars.size()];
    TH1D *hSpectraMultEstVars[multEstVars.size()];
    TH1D *hSpectraMaterialBudgetVars[materialBudgetVars.size()];

    openTFile(fDefault, basePathSigExt + correctedFileName);
    openTH1D(hSpectraDefault, fDefault, (multDir + "corrected_spectra_Integral_final").c_str());
    openTH1D(hEfficiency_Default, fDefault, (multDir + "heff").c_str());

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
        openTFile(fResidualBkgVars[i], basePathSigExt + ResidualBkgVars[i] + "/" + correctedFileName);
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

    // Track selection
    for (int i = 0; i < DCAvars.size(); i++)
    {
        openTFile(fDCAvars[i], basePathTrackSel + DCAvars[i] + "/hInvMass/" + correctedFileName);
        openTH1D(hSpectraDCAvars[i], fDCAvars[i], (multDir + "corrected_spectra_Integral_final").c_str());
    }
    for (int i = 0; i < PVcontributorVars.size(); i++)
    {
        openTFile(fPVcontributorVars[i], basePathTrackSel + PVcontributorVars[i] + "/hInvMass/" + correctedFileName);
        openTH1D(hSpectraPVcontributorVars[i], fPVcontributorVars[i], (multDir + "corrected_spectra_Integral_final").c_str());
    }

    // PID
    for (int i = 0; i < PIDVars.size(); i++)
    {
        openTFile(fPIDVars[i], basePathPIDAndMultEst + PIDVars[i] + "/hInvMass/" + correctedFileName);
        openTH1D(hSpectraPIDVars[i], fPIDVars[i], (multDir + "corrected_spectra_Integral_final").c_str());
    }

    // Multiplicity estimator
    for (int i = 0; i < multEstVars.size(); i++)
    {
        openTFile(fmultEstVars[i], basePathPIDAndMultEst + multEstVars[i] + "/hInvMass/" + correctedFileName);
        openTH1D(hSpectraMultEstVars[i], fmultEstVars[i], (multDir + "corrected_spectra_Integral_final").c_str());
    }

    // Material budget
    for (int i = 0; i < materialBudgetVars.size(); i++)
    {
        openTFile(fmaterialBudgetVars[i], basePathSigExt + materialBudgetVars[i] + "/" + correctedFileName);
        openTH1D(hSpectraMaterialBudgetVars[i], fmaterialBudgetVars[i], (multDir + "heff").c_str());
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
    vector<TH1 *> DCAVariationHists;
    for (int i = 0; i < DCAvars.size(); i++)
    {
        DCAVariationHists.push_back(hSpectraDCAvars[i]);
    }
    vector<TH1 *> PVContributorVariationHists;
    for (int i = 0; i < PVcontributorVars.size(); i++)
    {
        PVContributorVariationHists.push_back(hSpectraPVcontributorVars[i]);
    }
    vector<TH1 *> PIDVariationHists;
    for (int i = 0; i < PIDVars.size(); i++)
    {
        PIDVariationHists.push_back(hSpectraPIDVars[i]);
    }
    vector<TH1 *> MultEstVariationHists;
    for (int i = 0; i < multEstVars.size(); i++)
    {
        MultEstVariationHists.push_back(hSpectraMultEstVars[i]);
    }
    vector<TH1 *> MaterialBudgetVariationHists;
    for (int i = 0; i < materialBudgetVars.size(); i++)
    {
        MaterialBudgetVariationHists.push_back(hSpectraMaterialBudgetVars[i]);
    }

    vector<size_t> AllVarSizes = {
        NormVariationHists.size(),
        FitRangeVariationHists.size(),
        CombinatorialBkgVariationHists.size(),
        ResidualBkgVariationHists.size(),
        BinCountingVariationHists.size(),
        WidthVariationHists.size(),
        DCAVariationHists.size(),
        PVContributorVariationHists.size(),
        PIDVariationHists.size(),
        MultEstVariationHists.size()};
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
        WidthVariationHists,
        DCAVariationHists,
        PVContributorVariationHists,
        PIDVariationHists,
        MultEstVariationHists};

    vector<vector<string>> AllVariationNames = {
        normVars,
        fitRangeVars,
        CombinatorialBkgVars,
        ResidualBkgVars,
        BinCounting,
        widthVars,
        DCAvars,
        PVcontributorVars,
        PIDVars,
        multEstVars,
        materialBudgetVars};

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
            VarRatios[i].push_back(operations.CalculateRatio(hSpectraDefault, AllVariationHists[i][j]));
            histBarlow[i].push_back(operations.barlowcheck(hSpectraDefault, AllVariationHists[i][j], checkbar));
            cout << "Systematic source " << i + 1 << ", Variation " << j + 1 << ", Barlow check " << (checkbar ? "Passed" : "Failed") << endl;
            barlowPassed.push_back(checkbar);
            counter++;
        }
    }

    TCanvas *cPlotBarlowAll = new TCanvas("", "Barlow Checks", 1080, 720);
    SetCanvasStyle(cPlotBarlowAll, 0.15, 0.03, 0.06, 0.15);
    cPlotBarlowAll->Divide(4, 4);
    TCanvas *cRatioAll = new TCanvas("", "Spectra ratio (Variation / Default)", 1080, 720);
    SetCanvasStyle(cRatioAll, 0.15, 0.03, 0.06, 0.15);
    cRatioAll->Divide(4, 4);
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
            gPad->SetLeftMargin(0.09);
            gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.03);
            gPad->SetTopMargin(0.1);
            SetHistoQA(VarRatios[i][j]);
            VarRatios[i][j]->GetXaxis()->SetLabelSize(0.06);
            VarRatios[i][j]->GetYaxis()->SetLabelSize(0.06);
            VarRatios[i][j]->GetYaxis()->SetTitle("Variation / Default");
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
    TH1D *hRelUncertCombinatorialBkgVars = operations.RelativeUncertainty(hSpectraDefault, CombinatorialBkgVariationHists);
    TH1D *hRelUncertResidualBkgVars = operations.RelativeUncertainty(hSpectraDefault, ResidualBkgVariationHists);
    TH1D *hRelUncertBinCounting = operations.RelativeUncertainty(hSpectraDefault, BinCountingVariationHists);
    TH1D *hRelUncertWidthVars = operations.RelativeUncertainty(hSpectraDefault, WidthVariationHists);
    TH1D *hRelUncertDCAvars = operations.RelativeUncertainty(hSpectraDefault, DCAVariationHists);
    TH1D *hRelUncertPVcontributorVars = operations.RelativeUncertainty(hSpectraDefault, PVContributorVariationHists);
    TH1D *hRelUncertPIDVars = operations.RelativeUncertainty(hSpectraDefault, PIDVariationHists);
    TH1D *hRelUncertMultEstVars = operations.RelativeUncertainty(hSpectraDefault, MultEstVariationHists);
    TH1D *hRelUncertMaterialBudgetVars = operations.RelativeUncertainty(hEfficiency_Default, MaterialBudgetVariationHists);

    TCanvas *cRelUncert = new TCanvas("", "Relative Uncertainties", 1080, 720);
    SetCanvasStyle(cRelUncert, 0.15, 0.03, 0.06, 0.15);
    cRelUncert->Divide(4, 3);
    vector<TH1D *> relUncertHists = {
        (TH1D *)hRelUncertNormVars->Clone("hRelUncertNormVars_clone"),
        (TH1D *)hRelUncertFitRangeVars->Clone("hRelUncertFitRangeVars_clone"),
        (TH1D *)hRelUncertCombinatorialBkgVars->Clone("hRelUncertCombinatorialBkgVars_clone"),
        (TH1D *)hRelUncertResidualBkgVars->Clone("hRelUncertResidualBkgVars_clone"),
        (TH1D *)hRelUncertBinCounting->Clone("hRelUncertBinCounting_clone"),
        (TH1D *)hRelUncertWidthVars->Clone("hRelUncertWidthVars_clone"),
        (TH1D *)hRelUncertDCAvars->Clone("hRelUncertDCAvars_clone"),
        (TH1D *)hRelUncertPVcontributorVars->Clone("hRelUncertPVcontributorVars_clone"),
        (TH1D *)hRelUncertPIDVars->Clone("hRelUncertPIDVars_clone"),
        (TH1D *)hRelUncertMultEstVars->Clone("hRelUncertMultEstVars_clone"),
        (TH1D *)hRelUncertMaterialBudgetVars->Clone("hRelUncertMaterialBudgetVars_clone")};

    vector<string> relUncertNames = {
        "Norm. Range",
        "Fit Range",
        "Combinatorial Bkg",
        "Residual Bkg",
        "Yield Extraction",
        "Width Variation",
        "DCA Selection",
        "PV Contributor",
        "PID Selection",
        "Mult. Estimator",
        "Material Budget"};

    for (size_t i = 0; i < relUncertHists.size(); i++)
    {
        cRelUncert->cd(i + 1);
        gPad->SetLeftMargin(0.09);
        gPad->SetBottomMargin(0.14);
        gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.1);
        SetHistoQA(relUncertHists[i]);
        relUncertHists[i]->GetXaxis()->SetLabelSize(0.06);
        relUncertHists[i]->GetYaxis()->SetLabelSize(0.06);
        relUncertHists[i]->GetYaxis()->SetTitle("Relative Uncertainty");
        relUncertHists[i]->SetTitle(Form("%s", relUncertNames[i].c_str()));
        relUncertHists[i]->SetStats(0);
        relUncertHists[i]->SetMaximum(relUncertHists[i]->GetMaximum() * 1.05);
        relUncertHists[i]->Draw("HIST");
    }

    vector<TH1D *> vecSignalExt = {hRelUncertNormVars, hRelUncertFitRangeVars, hRelUncertCombinatorialBkgVars, hRelUncertResidualBkgVars, hRelUncertWidthVars}; // Bin counting excluded as it passed Barlow
    vector<TH1D *> vecTrackSel = {hRelUncertDCAvars, hRelUncertPVcontributorVars};
    vector<TH1D *> vecPID = {hRelUncertPIDVars};
    vector<TH1D *> vecMultEst = {hRelUncertMultEstVars};
    vector<TH1D *> vecMaterialBudget = {hRelUncertMaterialBudgetVars};
    TH1D *hITSTPCMatchingRelUncert = (TH1D *)hRelUncertMaterialBudgetVars->Clone("hITSTPCMatchingRelUncert");
    for (int ibin = 1; ibin <= hITSTPCMatchingRelUncert->GetNbinsX(); ibin++)
    {
        hITSTPCMatchingRelUncert->SetBinContent(ibin, 0.02); // Taking 2% relative uncertainty.
    }
    vector<TH1D *> vecITSTPCMatching = {hITSTPCMatchingRelUncert};

    TH1D *hHadronicInteractionRelUncert = (TH1D *)hRelUncertMaterialBudgetVars->Clone("hHadronicInteractionRelUncert");

    float pTValues[] = {0.4046, 0.7965, 1.201, 1.606, 1.997, 4.008};
    float relUncertHI[] = {0.02713, 0.02373, 0.01916, 0.01979, 0.02042, 0.01182};

    int nValues = sizeof(pTValues) / sizeof(float);

    // First set everything to zero
    for (int ibin = 1; ibin <= hHadronicInteractionRelUncert->GetNbinsX(); ibin++)
    {
        hHadronicInteractionRelUncert->SetBinContent(ibin, 0.);
    }

    // Fill only matching pT bins
    for (int ibin = 1; ibin <= hHadronicInteractionRelUncert->GetNbinsX(); ibin++)
    {
        double binCenter = hHadronicInteractionRelUncert->GetBinCenter(ibin);

        for (int i = 0; i < nValues; i++)
        {
            if (fabs(binCenter - pTValues[i]) < 0.2) // tolerance
            {
                hHadronicInteractionRelUncert->SetBinContent(ibin, relUncertHI[i]);
                break;
            }
        }
    }
    vector<TH1D *> vecHadronicInteraction = {hHadronicInteractionRelUncert};

    TH1D *hSignalExtTotalSys = operations.sigma(vecSignalExt);
    TH1D *hTrackSelTotalSys = operations.sigma(vecTrackSel);
    TH1D *hPIDTotalSys = operations.sigma(vecPID);
    TH1D *hMultEstTotalSys = operations.sigma(vecMultEst);
    TH1D *hMaterialBudgetTotalSys = operations.sigma(vecMaterialBudget);
    TH1D *hITSTPCMatchingTotalSys = operations.sigma(vecITSTPCMatching);
    TH1D *hHadronicInteractionTotalSys = operations.sigma(vecHadronicInteraction);

    vector<TH1D *> vecTotal = {hSignalExtTotalSys, hTrackSelTotalSys, hPIDTotalSys, hMultEstTotalSys, hMaterialBudgetTotalSys, hITSTPCMatchingTotalSys, hHadronicInteractionTotalSys};
    TH1D *hTotalSys = operations.sigma(vecTotal);

    TH1D *hSignalExtTotalSysClone = (TH1D *)hSignalExtTotalSys->Clone();
    TH1D *hTrackSelTotalSysClone = (TH1D *)hTrackSelTotalSys->Clone();

    string SigExtNames[] = {"Norm. range", "Fit Range", "Combinatorial Bkg", "Residual Bkg", "Width fix/free"};

    TCanvas *cSigExtAll = new TCanvas("", "Systematic Uncertainties from all sources", 720, 720);
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
        vecSignalExt[i]->SetMaximum(0.305);
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

    TCanvas *cTrackSelAll = new TCanvas("", "Systematic Uncertainties from Track Selection", 720, 720);
    SetCanvasStyle(cTrackSelAll, 0.14, 0.03, 0.06, 0.13);
    TLegend *legTrackSel = new TLegend(0.17, 0.6, 0.5, 0.88);
    legTrackSel->SetBorderSize(0);
    legTrackSel->SetFillStyle(0);
    legTrackSel->SetTextSize(0.03);
    legTrackSel->SetTextFont(42);
    legTrackSel->SetHeader("Track Selection");
    legTrackSel->AddEntry((TObject *)0, Form("Multiplicity: %d-%d", multLow, multHigh), "");
    for (int i = 0; i < vecTrackSel.size(); i++)
    {
        SetHistoQA(vecTrackSel[i]);
        vecTrackSel[i]->GetYaxis()->SetTitle("Relative Uncertainty");
        vecTrackSel[i]->SetStats(0);
        vecTrackSel[i]->SetMaximum(0.105);
        vecTrackSel[i]->SetMinimum(0);
        vecTrackSel[i]->SetLineColor(lineColors[i]);
        vecTrackSel[i]->Draw("HIST SAME");
        legTrackSel->AddEntry(vecTrackSel[i], Form("%s", AllVariationNames[i + 6][0].c_str()), "l");
    }
    SetHistoQA(hTrackSelTotalSysClone);
    hTrackSelTotalSysClone->SetLineColor(kBlack);
    hTrackSelTotalSysClone->SetLineWidth(3);
    hTrackSelTotalSysClone->Draw("HIST SAME");
    legTrackSel->AddEntry(hTrackSelTotalSysClone, "Total", "l");
    legTrackSel->Draw();

    TCanvas *cTotalSys = new TCanvas("", "Total Systematic Uncertainties", 720, 720);
    SetCanvasStyle(cTotalSys, 0.14, 0.03, 0.06, 0.13);
    TLegend *legTotal = new TLegend(0.17, 0.6, 0.5, 0.88);
    legTotal->SetBorderSize(0);
    legTotal->SetFillStyle(0);
    legTotal->SetTextSize(0.03);
    legTotal->SetTextFont(42);
    legTotal->SetHeader("Total Systematic");
    legTotal->AddEntry((TObject *)0, Form("Multiplicity: %d-%d", multLow, multHigh), "");
    for (int i = 0; i < vecTotal.size(); i++)
    {
        SetHistoQA(vecTotal[i]);
        vecTotal[i]->GetYaxis()->SetTitle("Relative Uncertainty");
        vecTotal[i]->SetStats(0);
        vecTotal[i]->SetMaximum(0.305);
        vecTotal[i]->SetMinimum(0);
        vecTotal[i]->SetLineColor(lineColors[i]);
        vecTotal[i]->Draw("HIST SAME");
        if (i == 0)
            legTotal->AddEntry(vecTotal[i], "Signal Extraction", "l");
        else if (i == 1)
            legTotal->AddEntry(vecTotal[i], "Track Selection", "l");
        else if (i == 2)
            legTotal->AddEntry(vecTotal[i], "PID", "l");
        else if (i == 3)
            legTotal->AddEntry(vecTotal[i], "Multiplicity Estimator", "l");
        else if (i == 4)
            legTotal->AddEntry(vecTotal[i], "Material Budget", "l");
        else if (i == 5)
            legTotal->AddEntry(vecTotal[i], "ITS-TPC Matching", "l");
        else if (i == 6)
            legTotal->AddEntry(vecTotal[i], "Hadronic Interaction", "l");
    }
    hTotalSys->SetLineColor(kBlack);
    hTotalSys->SetLineWidth(3);
    hTotalSys->Draw("HIST SAME");
    legTotal->AddEntry(hTotalSys, "Total", "l");
    legTotal->Draw();

    // smoothing procedure on separate sources
    int Iterations = 2;
    TH1D *hSignalExtTotalSysSmoothed = operations.smooth(hSignalExtTotalSys, Iterations);
    TH1D *hTrackSelTotalSysSmoothed = operations.smooth(hTrackSelTotalSys, Iterations);
    TH1D *hPIDTotalSysSmoothed = operations.smooth(hPIDTotalSys, Iterations);
    TH1D *hMultEstTotalSysSmoothed = operations.smooth(hMultEstTotalSys, Iterations);
    TH1D *hMaterialBudgetTotalSysSmoothed = operations.smooth(hMaterialBudgetTotalSys, Iterations);
    TH1D *hITSTPCMatchingTotalSysSmoothed = operations.smooth(hITSTPCMatchingTotalSys, Iterations);
    TH1D *hHadronicInteractionTotalSysSmoothed = operations.smooth(hHadronicInteractionTotalSys, Iterations);

    vector<TH1D *> smoothedTotalVec = {hSignalExtTotalSysSmoothed, hTrackSelTotalSysSmoothed, hPIDTotalSysSmoothed, hMultEstTotalSysSmoothed, hMaterialBudgetTotalSysSmoothed, hITSTPCMatchingTotalSysSmoothed, hHadronicInteractionTotalSysSmoothed};
    TH1D *hTotalSysSmoothed = operations.smooth(hTotalSys);

    TCanvas *cSmoothedTotalSys = new TCanvas("", "Smoothed Total Systematic Uncertainties", 720, 720);
    SetCanvasStyle(cSmoothedTotalSys, 0.14, 0.03, 0.06, 0.13);
    for (int i = 0; i < smoothedTotalVec.size(); i++)
    {
        SetHistoQA(smoothedTotalVec[i]);
        smoothedTotalVec[i]->GetYaxis()->SetTitle("Relative Uncertainty");
        smoothedTotalVec[i]->SetStats(0);
        smoothedTotalVec[i]->SetMaximum(0.305);
        smoothedTotalVec[i]->SetMinimum(0);
        smoothedTotalVec[i]->SetLineColor(lineColors[i]);
        smoothedTotalVec[i]->SetLineWidth(2);
        smoothedTotalVec[i]->Draw("HIST SAME");
    }
    hTotalSysSmoothed->SetLineColor(kBlack);
    hTotalSysSmoothed->SetLineWidth(3);
    hTotalSysSmoothed->Draw("HIST SAME");
    legTotal->Draw();

    // Save all the plots
    cPlotBarlowAll->SaveAs(savePath + "BarlowChecks_AllVariations.png");
    cRatioAll->SaveAs(savePath + "Ratio_AllVariations.png");
    cRelUncert->SaveAs(savePath + "RelativeUncertainties_AllSources.png");
    cSigExtAll->SaveAs(savePath + "SignalExtractionSystematics.png");
    cTrackSelAll->SaveAs(savePath + "TrackSelectionSystematics.png");
    cTotalSys->SaveAs(savePath + "TotalSystematics.png");
    cSmoothedTotalSys->SaveAs(savePath + "SmoothedTotalSystematics.png");
}