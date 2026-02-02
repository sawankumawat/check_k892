#include "../src/style.h"

using namespace std;

void checkVectorSizeMatch(const vector<string> &varNames, const vector<int> &nVariations, const string &categoryName);
void calculateRelativeUncertainty(TH1F *hDefault, TH1F *hVariation, TH1F *hRelUncertainty, string variationName, string type);

void systematics_pt()
{
    // //========= Systematic information =================
    // 1. Signal extraction
    // 1a. Fit range: 3 variations (1 in left and 2 in right)
    // 1b. Normalization range: 2 variations (left and right)
    // 1c. Fit function: 4 variations (2 in background and 2 in signal)
    // 1d. All parameters free/fixed: 2 variations

    // 2. Track selection
    // 2a. DCA to primary vertex: 1 Variation
    // 2b. TPC PID: 2 Variations
    // 2c. TPC min clusters: 2 Variations

    // 3. Topological selections
    // 3a. CosPA: 2 Variations
    // 3b. Decay radius: 1 Variation
    // 3c. DCA between daughters: 2 Variations
    // 3d. Mass tolerance: 2 Variations
    // 3e. Lambda rejection: 2 Variations
    // 3f. Lifetime cut: 2 Variations

    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/systematic2022_new/KsKs_Channel/higher-mass-resonances";

    vector<string> SignalExtractionVars = {"_fitLow1p07", "_fitHigh2p17", "_fitHigh2p25", "_normLeft", "_normRight", "_FitChKstar", "_FitExpoHERA", "_ConstWidth", "_Fit3rBW", "_AllParametersFree", "_AllParametersFixed"};
    vector<string> TrackSelectionVars = {"_DCA0p1", "_TPCPID2", "_TPCPID5", "_TPCMinCls100", "_TPCMinCls60"};
    vector<string> TopologicalSelectionVars = {"_cospa0p95", "_cospa0p992", "_decay_rad1p0", "_DCAv0dau0p3", "_DCAv0dau1p0", "_Ks_selection2p5", "_Ks_selection5", "_lambda_rej4", "_lambda_rej6", "_lifetime15", "_lifetime25"};

    vector<int> nVariationsSignalExtraction = {3, 2, 4, 2};           // total 11
    vector<int> nVariationsTrackSelection = {1, 2, 2};                // total 5
    vector<int> nVariationsTopologicalSelection = {2, 1, 2, 2, 2, 2}; // total 11

    // Check size matches
    checkVectorSizeMatch(SignalExtractionVars, nVariationsSignalExtraction, "Signal Extraction");
    checkVectorSizeMatch(TrackSelectionVars, nVariationsTrackSelection, "Track Selection");
    checkVectorSizeMatch(TopologicalSelectionVars, nVariationsTopologicalSelection, "Topological Selection");

    // We have to do variations in the mass and raw yield of f2(1525) and f0(1710) resonances
    TFile *fDefault = new TFile((path + "/fits/FitParam.root").c_str(), "READ");
    if (fDefault->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hMass_f1525_Default = (TH1F *)fDefault->Get("Mult_0_100/hMass_1525");
    TH1F *hMass_f1710_Default = (TH1F *)fDefault->Get("Mult_0_100/hMass_1710");
    TH1F *hYield_f1525_Default = (TH1F *)fDefault->Get("Mult_0_100/hYield_1525");
    TH1F *hYield_f1710_Default = (TH1F *)fDefault->Get("Mult_0_100/hYield_1710");
    if (hMass_f1525_Default == nullptr || hMass_f1710_Default == nullptr || hYield_f1525_Default == nullptr || hYield_f1710_Default == nullptr)
    {
        cout << "Error: One of the default histograms not found!" << endl;
        return;
    }

    int totalVariationsSignalExtraction = SignalExtractionVars.size();
    int totalVariationsTrackSelection = TrackSelectionVars.size();
    int totalVariationsTopologicalSelection = TopologicalSelectionVars.size();
    TFile *fVarSignalExtraction[totalVariationsSignalExtraction];
    TFile *fVarTrackSelection[totalVariationsTrackSelection];
    TFile *fVarTopologicalSelection[totalVariationsTopologicalSelection];
    for (int i = 0; i < totalVariationsSignalExtraction; i++)
    {
        fVarSignalExtraction[i] = new TFile((path + "/fits/FitParam" + SignalExtractionVars[i] + ".root").c_str(), "READ");
        if (fVarSignalExtraction[i]->IsZombie())
        {
            cout << "Error opening file for variation: " << SignalExtractionVars[i] << endl;
            return;
        }
    }
    for (int i = 0; i < totalVariationsTrackSelection; i++)
    {
        fVarTrackSelection[i] = new TFile((path + TrackSelectionVars[i] + "/fits/FitParam" + TrackSelectionVars[i] + ".root").c_str(), "READ");
        if (fVarTrackSelection[i]->IsZombie())
        {
            cout << "Error opening file for variation: " << TrackSelectionVars[i] << endl;
            return;
        }
    }
    for (int i = 0; i < totalVariationsTopologicalSelection; i++)
    {
        fVarTopologicalSelection[i] = new TFile((path + TopologicalSelectionVars[i] + "/fits/FitParam" + TopologicalSelectionVars[i] + ".root").c_str(), "READ");
        if (fVarTopologicalSelection[i]->IsZombie())
        {
            cout << "Error opening file for variation: " << TopologicalSelectionVars[i] << endl;
            return;
        }
    }
    // Relative uncertainty histograms storage
    vector<TH1F *> hMassVar_f1525_relUncert_sigExt, hMassVar_f1710_relUncert_sigExt, hYieldVar_f1525_relUncert_sigExt, hYieldVar_f1710_relUncert_sigExt;
    vector<TH1F *> hMassVar_f1525_relUncert_trackSel, hMassVar_f1710_relUncert_trackSel, hYieldVar_f1525_relUncert_trackSel, hYieldVar_f1710_relUncert_trackSel;
    vector<TH1F *> hMassVar_f1525_relUncert_topSel, hMassVar_f1710_relUncert_topSel, hYieldVar_f1525_relUncert_topSel, hYieldVar_f1710_relUncert_topSel;
    TCanvas *cSigExtract[4], *cTrackSelect[4], *cTopSelect[4];
    string CanvasTypes[4] = {"Mass1525", "Mass1710", "Yield1525", "Yield1710"};
    for (int i = 0; i < 4; i++)
    {
        cSigExtract[i] = new TCanvas(Form("cSigExtract_%d", i), Form("Signal Extraction Systematics %d", i), 2880, 1440);
        cSigExtract[i]->Divide(4, 3);
        cTrackSelect[i] = new TCanvas(Form("cTrackSelect_%d", i), Form("Track Selection Systematics %d", i), 2880, 1440);
        cTrackSelect[i]->Divide(3, 2);
        cTopSelect[i] = new TCanvas(Form("cTopSelect_%d", i), Form("Topological Selection Systematics %d", i), 2880, 1440);
        cTopSelect[i]->Divide(4, 3);
    }

    // Load histograms for Signal Extraction variations
    for (int i = 0; i < totalVariationsSignalExtraction; i++)
    {
        TH1F *hMass_f1525_Var = (TH1F *)fVarSignalExtraction[i]->Get("Mult_0_100/hMass_1525");
        TH1F *hMass_f1710_Var = (TH1F *)fVarSignalExtraction[i]->Get("Mult_0_100/hMass_1710");
        TH1F *hYield_f1525_Var = (TH1F *)fVarSignalExtraction[i]->Get("Mult_0_100/hYield_1525");
        TH1F *hYield_f1710_Var = (TH1F *)fVarSignalExtraction[i]->Get("Mult_0_100/hYield_1710");
        if (hMass_f1525_Var == nullptr || hMass_f1710_Var == nullptr || hYield_f1525_Var == nullptr || hYield_f1710_Var == nullptr)
        {
            cout << "Error: One of the histograms not found for variation: " << SignalExtractionVars[i] << endl;
            return;
        }

        TH1F *hMassRelUncert_1525 = (TH1F *)hMass_f1525_Default->Clone(Form("hMassRelUncert_1525_%s", SignalExtractionVars[i].c_str()));
        hMassRelUncert_1525->SetTitle(Form("Variation: %s", SignalExtractionVars[i].c_str()));
        TH1F *hMassRelUncert_1710 = (TH1F *)hMass_f1710_Default->Clone(Form("hMassRelUncert_1710_%s", SignalExtractionVars[i].c_str()));
        hMassRelUncert_1710->SetTitle(Form("Variation: %s", SignalExtractionVars[i].c_str()));
        TH1F *hYieldRelUncert_1525 = (TH1F *)hYield_f1525_Default->Clone(Form("hYieldRelUncert_1525_%s", SignalExtractionVars[i].c_str()));
        hYieldRelUncert_1525->SetTitle(Form("Variation: %s", SignalExtractionVars[i].c_str()));
        TH1F *hYieldRelUncert_1710 = (TH1F *)hYield_f1710_Default->Clone(Form("hYieldRelUncert_1710_%s", SignalExtractionVars[i].c_str()));
        hYieldRelUncert_1710->SetTitle(Form("Variation: %s", SignalExtractionVars[i].c_str()));

        calculateRelativeUncertainty(hMass_f1525_Default, hMass_f1525_Var, hMassRelUncert_1525, SignalExtractionVars[i], "_Mass1525");
        calculateRelativeUncertainty(hMass_f1710_Default, hMass_f1710_Var, hMassRelUncert_1710, SignalExtractionVars[i], "_Mass1710");
        calculateRelativeUncertainty(hYield_f1525_Default, hYield_f1525_Var, hYieldRelUncert_1525, SignalExtractionVars[i], "_Yield1525");
        calculateRelativeUncertainty(hYield_f1710_Default, hYield_f1710_Var, hYieldRelUncert_1710, SignalExtractionVars[i], "_Yield1710");

        hMassVar_f1525_relUncert_sigExt.push_back(hMassRelUncert_1525);
        hMassVar_f1710_relUncert_sigExt.push_back(hMassRelUncert_1710);
        hYieldVar_f1525_relUncert_sigExt.push_back(hYieldRelUncert_1525);
        hYieldVar_f1710_relUncert_sigExt.push_back(hYieldRelUncert_1710);

        for (int itype = 0; itype < 4; itype++)
        {
            cSigExtract[itype]->cd(i + 1);
            if (itype == 0)
                hMassRelUncert_1525->Draw("HIST");
            else if (itype == 1)
                hMassRelUncert_1710->Draw("HIST");
            else if (itype == 2)
                hYieldRelUncert_1525->Draw("HIST");
            else if (itype == 3)
                hYieldRelUncert_1710->Draw("HIST");
        }
    }
    for (int i = 0; i < 4; i++)
    {
        cSigExtract[i]->SaveAs(Form("%s/fits/SystematicPlots/SignalExtrct_Sys_%s.png", path.c_str(), CanvasTypes[i].c_str()));
    }
    // Similarly, load histograms for Track Selection and Topological Selection variations
    for (int i = 0; i < totalVariationsTrackSelection; i++)
    {
        TH1F *hMass_f1525_Var = (TH1F *)fVarTrackSelection[i]->Get("Mult_0_100/hMass_1525");
        TH1F *hMass_f1710_Var = (TH1F *)fVarTrackSelection[i]->Get("Mult_0_100/hMass_1710");
        TH1F *hYield_f1525_Var = (TH1F *)fVarTrackSelection[i]->Get("Mult_0_100/hYield_1525");
        TH1F *hYield_f1710_Var = (TH1F *)fVarTrackSelection[i]->Get("Mult_0_100/hYield_1710");
        if (hMass_f1525_Var == nullptr || hMass_f1710_Var == nullptr || hYield_f1525_Var == nullptr || hYield_f1710_Var == nullptr)
        {
            cout << "Error: One of the histograms not found for variation: " << TrackSelectionVars[i] << endl;
            return;
        }

        TH1F *hMassRelUncert_1525 = (TH1F *)hMass_f1525_Default->Clone(Form("hMassRelUncert_1525_%s", TrackSelectionVars[i].c_str()));
        hMassRelUncert_1525->SetTitle(Form("Variation: %s", TrackSelectionVars[i].c_str()));
        TH1F *hMassRelUncert_1710 = (TH1F *)hMass_f1710_Default->Clone(Form("hMassRelUncert_1710_%s", TrackSelectionVars[i].c_str()));
        hMassRelUncert_1710->SetTitle(Form("Variation: %s", TrackSelectionVars[i].c_str()));
        TH1F *hYieldRelUncert_1525 = (TH1F *)hYield_f1525_Default->Clone(Form("hYieldRelUncert_1525_%s", TrackSelectionVars[i].c_str()));
        hYieldRelUncert_1525->SetTitle(Form("Variation: %s", TrackSelectionVars[i].c_str()));
        TH1F *hYieldRelUncert_1710 = (TH1F *)hYield_f1710_Default->Clone(Form("hYieldRelUncert_1710_%s", TrackSelectionVars[i].c_str()));
        hYieldRelUncert_1710->SetTitle(Form("Variation: %s", TrackSelectionVars[i].c_str()));

        calculateRelativeUncertainty(hMass_f1525_Default, hMass_f1525_Var, hMassRelUncert_1525, TrackSelectionVars[i], "_Mass1525");
        calculateRelativeUncertainty(hMass_f1710_Default, hMass_f1710_Var, hMassRelUncert_1710, TrackSelectionVars[i], "_Mass1710");
        calculateRelativeUncertainty(hYield_f1525_Default, hYield_f1525_Var, hYieldRelUncert_1525, TrackSelectionVars[i], "_Yield1525");
        calculateRelativeUncertainty(hYield_f1710_Default, hYield_f1710_Var, hYieldRelUncert_1710, TrackSelectionVars[i], "_Yield1710");

        hMassVar_f1525_relUncert_trackSel.push_back(hMassRelUncert_1525);
        hMassVar_f1710_relUncert_trackSel.push_back(hMassRelUncert_1710);
        hYieldVar_f1525_relUncert_trackSel.push_back(hYieldRelUncert_1525);
        hYieldVar_f1710_relUncert_trackSel.push_back(hYieldRelUncert_1710);

        for (int itype = 0; itype < 4; itype++)
        {
            cTrackSelect[itype]->cd(i + 1);
            if (itype == 0)
                hMassRelUncert_1525->Draw("HIST");
            else if (itype == 1)
                hMassRelUncert_1710->Draw("HIST");
            else if (itype == 2)
                hYieldRelUncert_1525->Draw("HIST");
            else if (itype == 3)
                hYieldRelUncert_1710->Draw("HIST");
        }
    }
    for (int i = 0; i < 4; i++)
    {
        cTrackSelect[i]->SaveAs(Form("%s/fits/SystematicPlots/TrackSelect_Sys_%s.png", path.c_str(), CanvasTypes[i].c_str()));
    }
    for (int i = 0; i < totalVariationsTopologicalSelection; i++)
    {
        TH1F *hMass_f1525_Var = (TH1F *)fVarTopologicalSelection[i]->Get("Mult_0_100/hMass_1525");
        TH1F *hMass_f1710_Var = (TH1F *)fVarTopologicalSelection[i]->Get("Mult_0_100/hMass_1710");
        TH1F *hYield_f1525_Var = (TH1F *)fVarTopologicalSelection[i]->Get("Mult_0_100/hYield_1525");
        TH1F *hYield_f1710_Var = (TH1F *)fVarTopologicalSelection[i]->Get("Mult_0_100/hYield_1710");
        if (hMass_f1525_Var == nullptr || hMass_f1710_Var == nullptr || hYield_f1525_Var == nullptr || hYield_f1710_Var == nullptr)
        {
            cout << "Error: One of the histograms not found for variation: " << TopologicalSelectionVars[i] << endl;
            return;
        }

        TH1F *hMassRelUncert_1525 = (TH1F *)hMass_f1525_Default->Clone(Form("hMassRelUncert_1525_%s", TopologicalSelectionVars[i].c_str()));
        hMassRelUncert_1525->SetTitle(Form("Variation: %s", TopologicalSelectionVars[i].c_str()));
        TH1F *hMassRelUncert_1710 = (TH1F *)hMass_f1710_Default->Clone(Form("hMassRelUncert_1710_%s", TopologicalSelectionVars[i].c_str()));
        hMassRelUncert_1710->SetTitle(Form("Variation: %s", TopologicalSelectionVars[i].c_str()));
        TH1F *hYieldRelUncert_1525 = (TH1F *)hYield_f1525_Default->Clone(Form("hYieldRelUncert_1525_%s", TopologicalSelectionVars[i].c_str()));
        hYieldRelUncert_1525->SetTitle(Form("Variation: %s", TopologicalSelectionVars[i].c_str()));
        TH1F *hYieldRelUncert_1710 = (TH1F *)hYield_f1710_Default->Clone(Form("hYieldRelUncert_1710_%s", TopologicalSelectionVars[i].c_str()));
        hYieldRelUncert_1710->SetTitle(Form("Variation: %s", TopologicalSelectionVars[i].c_str()));

        calculateRelativeUncertainty(hMass_f1525_Default, hMass_f1525_Var, hMassRelUncert_1525, TopologicalSelectionVars[i], "_Mass1525");
        calculateRelativeUncertainty(hMass_f1710_Default, hMass_f1710_Var, hMassRelUncert_1710, TopologicalSelectionVars[i], "_Mass1710");
        calculateRelativeUncertainty(hYield_f1525_Default, hYield_f1525_Var, hYieldRelUncert_1525, TopologicalSelectionVars[i], "_Yield1525");
        calculateRelativeUncertainty(hYield_f1710_Default, hYield_f1710_Var, hYieldRelUncert_1710, TopologicalSelectionVars[i], "_Yield1710");

        hMassVar_f1525_relUncert_topSel.push_back(hMassRelUncert_1525);
        hMassVar_f1710_relUncert_topSel.push_back(hMassRelUncert_1710);
        hYieldVar_f1525_relUncert_topSel.push_back(hYieldRelUncert_1525);
        hYieldVar_f1710_relUncert_topSel.push_back(hYieldRelUncert_1710);

        for (int itype = 0; itype < 4; itype++)
        {
            cTopSelect[itype]->cd(i + 1);
            if (itype == 0)
                hMassRelUncert_1525->Draw("HIST");
            else if (itype == 1)
                hMassRelUncert_1710->Draw("HIST");
            else if (itype == 2)
                hYieldRelUncert_1525->Draw("HIST");
            else if (itype == 3)
                hYieldRelUncert_1710->Draw("HIST");
        }
    }
    for (int i = 0; i < 4; i++)
    {
        cTopSelect[i]->SaveAs(Form("%s/fits/SystematicPlots/TopologicalSelect_Sys_%s.png", path.c_str(), CanvasTypes[i].c_str()));
    }
}

void checkVectorSizeMatch(const vector<string> &varNames, const vector<int> &nVariations, const string &categoryName)
{
    int totalVariations = 0;
    for (const auto &nVar : nVariations)
    {
        totalVariations += nVar;
    }

    if (varNames.size() != totalVariations)
    {
        cerr << "Error: Mismatch in number of variations for " << categoryName << " category." << endl;
        cerr << "Expected " << totalVariations << " variations, but got " << varNames.size() << " variation names." << endl;
        exit(EXIT_FAILURE);
    }
}

void calculateRelativeUncertainty(TH1F *hDefault, TH1F *hVariation, TH1F *hRelUncertainty, string variationName, string type)
{
    int nBins = hDefault->GetNbinsX();
    for (int bin = 1; bin <= nBins; bin++)
    {
        float defaultValue = hDefault->GetBinContent(bin);
        float variationValue = hVariation->GetBinContent(bin);
        if (defaultValue != 0)
        {
            float relUncertainty = fabs(variationValue - defaultValue) / defaultValue;
            hRelUncertainty->SetBinContent(bin, relUncertainty);
        }
    }
}