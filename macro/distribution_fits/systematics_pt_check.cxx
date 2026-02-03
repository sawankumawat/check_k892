#include "../src/style.h"

using namespace std;

void checkVectorSizeMatch(const vector<string> &varNames, const vector<int> &nVariations, const string &categoryName);
void calculateRelativeUncertainty(TH1F *hDefault, TH1F *hVariation, TH1F *hRelUncertainty);
void checkBarlowSignificance(TH1F *hDefault, TH1F *hVariation, TH1F *hBarlow, string variationName, string type, bool &barlowPassed);
void processVariations(TFile **fVariations, const vector<string> &varNames, TH1F *defaultHist, TCanvas *canvas, TCanvas *canvasBarlow, vector<TH1F *> *resultVector, const string &histPath, const string &suffix, const string &path, const string &categoryName, const string &categoryType);
void averageQuadratureSum(const vector<TH1F *> &hVariationHistos, const int nVariations, TH1F *&hAvgQuadratureSum, const char *name);
void QuadratureSum(const vector<TH1F *> &hVariationHistos, TH1F *&hQuadratureSum);

void systematics_pt()
{
    bool saveRelUncertHisto = true;
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
    Color_t lineColors[6] = {kRed, kBlue, kGreen + 2, kMagenta, kCyan + 2, kOrange + 7};

    vector<string> SignalExtractionVars = {"_fitLow1p07", "_fitHigh2p17", "_fitHigh2p25", "_normLeft", "_normRight", "_FitChKstar", "_FitExpoHERA", "_ConstWidth", "_Fit3rBW", "_AllParametersFree", "_AllParametersFixed"};
    vector<string> TrackSelectionVars = {"_DCA0p1", "_TPCPID2", "_TPCPID5", "_TPCMinCls100", "_TPCMinCls60"};
    vector<string> TopologicalSelectionVars = {"_cospa0p95", "_cospa0p992", "_decay_rad1p0", "_DCAv0dau0p3", "_DCAv0dau1p0", "_Ks_selection2p5", "_Ks_selection5", "_lambda_rej4", "_lambda_rej6", "_lifetime15", "_lifetime25"};

    vector<int> nVarSigExt = {3, 2, 4, 2};    // total 11
    vector<int> nVarTrk = {1, 2, 2};          // total 5
    vector<int> nVarTop = {2, 1, 2, 2, 2, 2}; // total 11

    vector<string> SourcesSignalExtraction = {"Fit Range", "Normalization Range", "Fit Function", "Fit Parameters Free/Fixed"};
    vector<string> SourcesTrackSelection = {"DCA to PV", "TPC PID", "TPC Min Clusters"};
    vector<string> SourcesTopologicalSelection = {"CosPA", "Decay Radius", "DCA between Daughters", "Mass Tolerance", "Lambda Rejection", "Lifetime Cut"};

    // Check size matches
    checkVectorSizeMatch(SignalExtractionVars, nVarSigExt, "Signal Extraction");
    checkVectorSizeMatch(TrackSelectionVars, nVarTrk, "Track Selection");
    checkVectorSizeMatch(TopologicalSelectionVars, nVarTop, "Topological Selection");

    // We have to do variations in the mass and raw yield of f2(1525) and f0(1710) resonances
    TFile *fDefault = new TFile((path + "/fits/FitParam.root").c_str(), "READ");
    // TFile *fCorrYield = new TFile((path + "/mult_0-100/Spectra/ReweightedSpectra_.root").c_str(), "READ");
    TFile *fCorrYield = new TFile((path + "/mult_0-100/Spectra/spectra_.root").c_str(), "READ");
    if (fDefault->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hMass_f1525_Default = (TH1F *)fDefault->Get("Mult_0_100/hMass_1525");
    TH1F *hMass_f1710_Default = (TH1F *)fDefault->Get("Mult_0_100/hMass_1710");
    TH1F *hYield_f1525_Default = (TH1F *)fDefault->Get("Mult_0_100/hYield_1525");
    TH1F *hYield_f1710_Default = (TH1F *)fDefault->Get("Mult_0_100/hYield_1710");
    // TH1F *hYield_f1525_Default = (TH1F *)fCorrYield->Get("f21525_Reweighted_Yield");
    // TH1F *hYield_f1710_Default = (TH1F *)fCorrYield->Get("f01710_Reweighted_Yield");
    // TH1F *hYield_f1525_Default = (TH1F *)fCorrYield->Get("hYield1525Corrected");
    // TH1F *hYield_f1710_Default = (TH1F *)fCorrYield->Get("hYield1710Corrected");
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
    TFile *fVarSignalCorrYield[totalVariationsSignalExtraction];
    TFile *fVarTrackCorrYield[totalVariationsTrackSelection];
    TFile *fVarTopologicalCorrYield[totalVariationsTopologicalSelection];
    for (int i = 0; i < totalVariationsSignalExtraction; i++)
    {
        fVarSignalExtraction[i] = new TFile((path + "/fits/FitParam" + SignalExtractionVars[i] + ".root").c_str(), "READ");
        // fVarSignalCorrYield[i] = new TFile((path + "/mult_0-100/Spectra/ReweightedSpectra_" + SignalExtractionVars[i] + ".root").c_str(), "READ");
        fVarSignalCorrYield[i] = new TFile((path + "/mult_0-100/Spectra/spectra_" + SignalExtractionVars[i] + ".root").c_str(), "READ");
        if (fVarSignalExtraction[i]->IsZombie() || fVarSignalCorrYield[i]->IsZombie())
        {
            cout << "Error opening file for variation: " << SignalExtractionVars[i] << endl;
            return;
        }
    }
    for (int i = 0; i < totalVariationsTrackSelection; i++)
    {
        fVarTrackSelection[i] = new TFile((path + TrackSelectionVars[i] + "/fits/FitParam" + TrackSelectionVars[i] + ".root").c_str(), "READ");
        // fVarTrackCorrYield[i] = new TFile((path + TrackSelectionVars[i] + "/mult_0-100/Spectra/ReweightedSpectra_" + TrackSelectionVars[i] + ".root").c_str(), "READ");
        fVarTrackCorrYield[i] = new TFile((path + TrackSelectionVars[i] + "/mult_0-100/Spectra/spectra_" + TrackSelectionVars[i] + ".root").c_str(), "READ");
        if (fVarTrackSelection[i]->IsZombie() || fVarTrackCorrYield[i]->IsZombie())
        {
            cout << "Error opening file for variation: " << TrackSelectionVars[i] << endl;
            return;
        }
    }
    for (int i = 0; i < totalVariationsTopologicalSelection; i++)
    {
        fVarTopologicalSelection[i] = new TFile((path + TopologicalSelectionVars[i] + "/fits/FitParam" + TopologicalSelectionVars[i] + ".root").c_str(), "READ");
        // fVarTopologicalCorrYield[i] = new TFile((path + TopologicalSelectionVars[i] + "/mult_0-100/Spectra/ReweightedSpectra_" + TopologicalSelectionVars[i] + ".root").c_str(), "READ");
        fVarTopologicalCorrYield[i] = new TFile((path + TopologicalSelectionVars[i] + "/mult_0-100/Spectra/spectra_" + TopologicalSelectionVars[i] + ".root").c_str(), "READ");
        if (fVarTopologicalSelection[i]->IsZombie() || fVarTopologicalCorrYield[i]->IsZombie())
        {
            cout << "Error opening file for variation: " << TopologicalSelectionVars[i] << endl;
            return;
        }
    }

    // Arrays for each category type
    string CanvasTypes[4] = {"Mass1525", "Mass1710", "Yield1525", "Yield1710"};
    // string histPaths[4] = {"Mult_0_100/hMass_1525", "Mult_0_100/hMass_1710", "f21525_Reweighted_Yield", "f01710_Reweighted_Yield"};
    string histPaths[4] = {"Mult_0_100/hMass_1525", "Mult_0_100/hMass_1710", "hYield1525Corrected", "hYield1710Corrected"};
    TH1F *defaultHists[4] = {hMass_f1525_Default, hMass_f1710_Default, hYield_f1525_Default, hYield_f1710_Default};
    string suffixes[4] = {"_Mass1525", "_Mass1710", "_Yield1525", "_Yield1710"};
    vector<TFile **> fVarSignalExtractionVec = {fVarSignalExtraction, fVarSignalExtraction, fVarSignalExtraction, fVarSignalExtraction};
    vector<TFile **> fVarTrackSelectionVec = {fVarTrackSelection, fVarTrackSelection, fVarTrackSelection, fVarTrackSelection};
    vector<TFile **> fVarTopologicalSelectionVec = {fVarTopologicalSelection, fVarTopologicalSelection, fVarTopologicalSelection, fVarTopologicalSelection};
    // vector<TFile **> fVarSignalExtractionVec = {fVarSignalExtraction, fVarSignalExtraction, fVarSignalCorrYield, fVarSignalCorrYield};
    // vector<TFile **> fVarTrackSelectionVec = {fVarTrackSelection, fVarTrackSelection, fVarTrackCorrYield, fVarTrackCorrYield};
    // vector<TFile **> fVarTopologicalSelectionVec = {fVarTopologicalSelection, fVarTopologicalSelection, fVarTopologicalCorrYield, fVarTopologicalCorrYield};

    // Loop over 4 category types
    for (int typeIdx = 0; typeIdx < 4; typeIdx++)
    {
        TH1F *defaultHist = defaultHists[typeIdx];
        string histPath = histPaths[typeIdx];
        string suffix = suffixes[typeIdx];
        string canvasType = CanvasTypes[typeIdx];

        // Relative uncertainty histograms storage for this type
        vector<TH1F *> hMassVar_relUncert_sigExt, hMassVar_relUncert_trackSel, hMassVar_relUncert_topSel;

        // Create canvases for this type
        TCanvas *cSigExtract = new TCanvas(Form("cSigExtract_%s", canvasType.c_str()), Form("Signal Extraction Systematics %s", canvasType.c_str()), 2880, 1440);
        cSigExtract->Divide(4, 3);
        TCanvas *cTrackSelect = new TCanvas(Form("cTrackSelect_%s", canvasType.c_str()), Form("Track Selection Systematics %s", canvasType.c_str()), 2880, 1440);
        cTrackSelect->Divide(3, 2);
        TCanvas *cTopSelect = new TCanvas(Form("cTopSelect_%s", canvasType.c_str()), Form("Topological Selection Systematics %s", canvasType.c_str()), 2880, 1440);
        cTopSelect->Divide(4, 3);

        TCanvas *cSigExtract_Barlow = new TCanvas(Form("cSigExtract_Barlow_%s", canvasType.c_str()), Form("Signal Extraction Barlow Significance %s", canvasType.c_str()), 2880, 1440);
        cSigExtract_Barlow->Divide(4, 3);
        TCanvas *cTrackSelect_Barlow = new TCanvas(Form("cTrackSelect_Barlow_%s", canvasType.c_str()), Form("Track Selection Barlow Significance %s", canvasType.c_str()), 2880, 1440);
        cTrackSelect_Barlow->Divide(3, 2);
        TCanvas *cTopSelect_Barlow = new TCanvas(Form("cTopSelect_Barlow_%s", canvasType.c_str()), Form("Topological Selection Barlow Significance %s", canvasType.c_str()), 2880, 1440);
        cTopSelect_Barlow->Divide(4, 3);

        // Process variations for Signal Extraction
        processVariations(fVarSignalExtractionVec[typeIdx], SignalExtractionVars, defaultHist, cSigExtract, cSigExtract_Barlow, &hMassVar_relUncert_sigExt, histPath, suffix, path, "SignalExtrct", "Signal Extraction");
        cout << "Size of relative uncertainty vector for signal extraction: " << hMassVar_relUncert_sigExt.size() << endl;

        if (saveRelUncertHisto)
        {
            cSigExtract->SaveAs(Form("%s/fits/SystematicPlots/SignalExtrct_Sys_%s.png", path.c_str(), canvasType.c_str()));
            cSigExtract_Barlow->SaveAs(Form("%s/fits/SystematicPlots/SignalExtrct_Barlow_%s.png", path.c_str(), canvasType.c_str()));
        }

        // Process variations for Track Selection
        processVariations(fVarTrackSelectionVec[typeIdx], TrackSelectionVars, defaultHist, cTrackSelect, cTrackSelect_Barlow, &hMassVar_relUncert_trackSel, histPath, suffix, path, "TrackSelect", "Track Selection");

        if (saveRelUncertHisto)
        {
            cTrackSelect->SaveAs(Form("%s/fits/SystematicPlots/TrackSelect_Sys_%s.png", path.c_str(), canvasType.c_str()));
            cTrackSelect_Barlow->SaveAs(Form("%s/fits/SystematicPlots/TrackSelect_Barlow_%s.png", path.c_str(), canvasType.c_str()));
        }

        // Process variations for Topological Selection
        processVariations(fVarTopologicalSelectionVec[typeIdx], TopologicalSelectionVars, defaultHist, cTopSelect, cTopSelect_Barlow, &hMassVar_relUncert_topSel, histPath, suffix, path, "TopologicalSelect", "Topological Selection");

        if (saveRelUncertHisto)
        {
            cTopSelect->SaveAs(Form("%s/fits/SystematicPlots/TopologicalSelect_Sys_%s.png", path.c_str(), canvasType.c_str()));
            cTopSelect_Barlow->SaveAs(Form("%s/fits/SystematicPlots/TopologicalSelect_Barlow_%s.png", path.c_str(), canvasType.c_str()));
        }

        // Now lets calculate the average quadrature sum for each systematic source in signal extraction
        TH1F *hSigExt_source1;
        TH1F *hSigExt_source2;
        TH1F *hSigExt_source3;
        TH1F *hSigExt_source4;
        TH1F *hSigExt_total;
        vector<TH1F *> vec_source1(hMassVar_relUncert_sigExt.begin(), hMassVar_relUncert_sigExt.begin() + nVarSigExt[0]);
        vector<TH1F *> vec_source2(hMassVar_relUncert_sigExt.begin() + nVarSigExt[0],
                                   hMassVar_relUncert_sigExt.begin() + nVarSigExt[0] + nVarSigExt[1]);
        vector<TH1F *> vec_source3(hMassVar_relUncert_sigExt.begin() + nVarSigExt[0] + nVarSigExt[1],
                                   hMassVar_relUncert_sigExt.begin() + nVarSigExt[0] + nVarSigExt[1] + nVarSigExt[2]);
        vector<TH1F *> vec_source4(hMassVar_relUncert_sigExt.begin() + nVarSigExt[0] + nVarSigExt[1] + nVarSigExt[2],
                                   hMassVar_relUncert_sigExt.end());

        averageQuadratureSum(vec_source1, nVarSigExt[0], hSigExt_source1, Form("hSigExt_AvgQuadSum_Source1_%s", suffix.c_str()));
        averageQuadratureSum(vec_source2, nVarSigExt[1], hSigExt_source2, Form("hSigExt_AvgQuadSum_Source2_%s", suffix.c_str()));
        averageQuadratureSum(vec_source3, nVarSigExt[2], hSigExt_source3, Form("hSigExt_AvgQuadSum_Source3_%s", suffix.c_str()));
        averageQuadratureSum(vec_source4, nVarSigExt[3], hSigExt_source4, Form("hSigExt_AvgQuadSum_Source4_%s", suffix.c_str()));
        vector<TH1F *> vec_total = {hSigExt_source1, hSigExt_source2, hSigExt_source3, hSigExt_source4};
        QuadratureSum(vec_total, hSigExt_total);

        TCanvas *cSigExt_sys = new TCanvas("", "", 720, 720);
        SetHistoQA(hSigExt_source1);
        hSigExt_source1->SetTitle(Form("%s", canvasType.c_str()));
        hSigExt_source1->SetLineColor(lineColors[0]);
        hSigExt_source1->SetMaximum(0.8);
        hSigExt_source1->SetMinimum(0);
        hSigExt_source1->Draw("HIST");
        hSigExt_source2->SetLineColor(lineColors[1]);
        hSigExt_source2->Draw("HIST SAME");
        hSigExt_source3->SetLineColor(lineColors[2]);
        hSigExt_source3->Draw("HIST SAME");
        hSigExt_source4->SetLineColor(lineColors[3]);
        hSigExt_source4->Draw("HIST SAME");
        hSigExt_total->SetLineColor(kBlack);
        hSigExt_total->SetLineWidth(3);
        hSigExt_total->Draw("HIST SAME");
        TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.85);
        leg->AddEntry(hSigExt_source1, SourcesSignalExtraction[0].c_str(), "l");
        leg->AddEntry(hSigExt_source2, SourcesSignalExtraction[1].c_str(), "l");
        leg->AddEntry(hSigExt_source3, SourcesSignalExtraction[2].c_str(), "l");
        leg->AddEntry(hSigExt_source4, SourcesSignalExtraction[3].c_str(), "l");
        leg->AddEntry(hSigExt_total, "Total Systematic", "l");
        leg->Draw();
        cSigExt_sys->SaveAs(Form("%s/fits/SystematicPlots/SignalExtrct_TotalSys_%s.png", path.c_str(), canvasType.c_str()));

        // Now lets calculate the average quadrature sum for each systematic source in track selection
        TH1F *hTrkSel_source1;
        TH1F *hTrkSel_source2;
        TH1F *hTrkSel_source3;
        TH1F *hTrkSel_total;
        vector<TH1F *> vec_trk_source1(hMassVar_relUncert_trackSel.begin(), hMassVar_relUncert_trackSel.begin() + nVarTrk[0]);
        vector<TH1F *> vec_trk_source2(hMassVar_relUncert_trackSel.begin() + nVarTrk[0],
                                       hMassVar_relUncert_trackSel.begin() + nVarTrk[0] + nVarTrk[1]);
        vector<TH1F *> vec_trk_source3(hMassVar_relUncert_trackSel.begin() + nVarTrk[0] + nVarTrk[1],
                                       hMassVar_relUncert_trackSel.end());

        averageQuadratureSum(vec_trk_source1, nVarTrk[0], hTrkSel_source1, Form("hTrkSel_AvgQuadSum_Source1_%s", suffix.c_str()));
        averageQuadratureSum(vec_trk_source2, nVarTrk[1], hTrkSel_source2, Form("hTrkSel_AvgQuadSum_Source2_%s", suffix.c_str()));
        averageQuadratureSum(vec_trk_source3, nVarTrk[2], hTrkSel_source3, Form("hTrkSel_AvgQuadSum_Source3_%s", suffix.c_str()));
        vector<TH1F *> vec_trk_total = {hTrkSel_source1, hTrkSel_source2, hTrkSel_source3};
        QuadratureSum(vec_trk_total, hTrkSel_total);

        TCanvas *cTrkSel_sys = new TCanvas("", "", 720, 720);
        SetHistoQA(hTrkSel_source1);
        hTrkSel_source1->SetTitle(Form("%s", canvasType.c_str()));
        hTrkSel_source1->SetLineColor(lineColors[0]);
        hTrkSel_source1->SetMaximum(0.8);
        hTrkSel_source1->SetMinimum(0);
        hTrkSel_source1->Draw("HIST");
        hTrkSel_source2->SetLineColor(lineColors[1]);
        hTrkSel_source2->Draw("HIST SAME");
        hTrkSel_source3->SetLineColor(lineColors[2]);
        hTrkSel_source3->Draw("HIST SAME");
        hTrkSel_total->SetLineColor(kBlack);
        hTrkSel_total->SetLineWidth(3);
        hTrkSel_total->Draw("HIST SAME");
        TLegend *legTrk = new TLegend(0.55, 0.65, 0.85, 0.85);
        legTrk->AddEntry(hTrkSel_source1, SourcesTrackSelection[0].c_str(), "l");
        legTrk->AddEntry(hTrkSel_source2, SourcesTrackSelection[1].c_str(), "l");
        legTrk->AddEntry(hTrkSel_source3, SourcesTrackSelection[2].c_str(), "l");
        legTrk->AddEntry(hTrkSel_total, "Total Systematic", "l");
        legTrk->Draw();
        cTrkSel_sys->SaveAs(Form("%s/fits/SystematicPlots/TrackSelect_TotalSys_%s.png", path.c_str(), canvasType.c_str()));

        // Now lets calculate the average quadrature sum for each systematic source in topological selection
        TH1F *hTopSel_source1;
        TH1F *hTopSel_source2;
        TH1F *hTopSel_source3;
        TH1F *hTopSel_source4;
        TH1F *hTopSel_source5;
        TH1F *hTopSel_source6;
        TH1F *hTopSel_total;
        vector<TH1F *> vec_top_source1(hMassVar_relUncert_topSel.begin(), hMassVar_relUncert_topSel.begin() + nVarTop[0]);
        vector<TH1F *> vec_top_source2(hMassVar_relUncert_topSel.begin() + nVarTop[0],
                                       hMassVar_relUncert_topSel.begin() + nVarTop[0] + nVarTop[1]);
        vector<TH1F *> vec_top_source3(hMassVar_relUncert_topSel.begin() + nVarTop[0] + nVarTop[1],
                                       hMassVar_relUncert_topSel.begin() + nVarTop[0] + nVarTop[1] + nVarTop[2]);
        vector<TH1F *> vec_top_source4(hMassVar_relUncert_topSel.begin() + nVarTop[0] + nVarTop[1] + nVarTop[2],
                                       hMassVar_relUncert_topSel.begin() + nVarTop[0] + nVarTop[1] + nVarTop[2] + nVarTop[3]);
        vector<TH1F *> vec_top_source5(hMassVar_relUncert_topSel.begin() + nVarTop[0] + nVarTop[1] + nVarTop[2] + nVarTop[3],
                                       hMassVar_relUncert_topSel.begin() + nVarTop[0] + nVarTop[1] + nVarTop[2] + nVarTop[3] + nVarTop[4]);
        vector<TH1F *> vec_top_source6(hMassVar_relUncert_topSel.begin() + nVarTop[0] + nVarTop[1] + nVarTop[2] + nVarTop[3] + nVarTop[4],
                                       hMassVar_relUncert_topSel.end());

        averageQuadratureSum(vec_top_source1, nVarTop[0], hTopSel_source1, Form("hTopSel_AvgQuadSum_Source1_%s", suffix.c_str()));
        averageQuadratureSum(vec_top_source2, nVarTop[1], hTopSel_source2, Form("hTopSel_AvgQuadSum_Source2_%s", suffix.c_str()));
        averageQuadratureSum(vec_top_source3, nVarTop[2], hTopSel_source3, Form("hTopSel_AvgQuadSum_Source3_%s", suffix.c_str()));
        averageQuadratureSum(vec_top_source4, nVarTop[3], hTopSel_source4, Form("hTopSel_AvgQuadSum_Source4_%s", suffix.c_str()));
        averageQuadratureSum(vec_top_source5, nVarTop[4], hTopSel_source5, Form("hTopSel_AvgQuadSum_Source5_%s", suffix.c_str()));
        averageQuadratureSum(vec_top_source6, nVarTop[5], hTopSel_source6, Form("hTopSel_AvgQuadSum_Source6_%s", suffix.c_str()));
        vector<TH1F *> vec_top_total = {hTopSel_source1, hTopSel_source2, hTopSel_source3, hTopSel_source4, hTopSel_source5, hTopSel_source6};
        QuadratureSum(vec_top_total, hTopSel_total);

        TCanvas *cTopSel_sys = new TCanvas("", "", 720, 720);
        SetHistoQA(hTopSel_source1);
        hTopSel_source1->SetTitle(Form("%s", canvasType.c_str()));
        hTopSel_source1->SetLineColor(lineColors[0]);
        hTopSel_source1->SetMaximum(0.8);
        hTopSel_source1->SetMinimum(0);
        hTopSel_source1->Draw("HIST");
        hTopSel_source2->SetLineColor(lineColors[1]);
        hTopSel_source2->Draw("HIST SAME");
        hTopSel_source3->SetLineColor(lineColors[2]);
        hTopSel_source3->Draw("HIST SAME");
        hTopSel_source4->SetLineColor(lineColors[3]);
        hTopSel_source4->Draw("HIST SAME");
        hTopSel_source5->SetLineColor(lineColors[4]);
        hTopSel_source5->Draw("HIST SAME");
        hTopSel_source6->SetLineColor(lineColors[5]);
        hTopSel_source6->Draw("HIST SAME");
        hTopSel_total->SetLineColor(kBlack);
        hTopSel_total->SetLineWidth(3);
        hTopSel_total->Draw("HIST SAME");
        TLegend *legTop = new TLegend(0.55, 0.45, 0.85, 0.85);
        legTop->AddEntry(hTopSel_source1, SourcesTopologicalSelection[0].c_str(), "l");
        legTop->AddEntry(hTopSel_source2, SourcesTopologicalSelection[1].c_str(), "l");
        legTop->AddEntry(hTopSel_source3, SourcesTopologicalSelection[2].c_str(), "l");
        legTop->AddEntry(hTopSel_source4, SourcesTopologicalSelection[3].c_str(), "l");
        legTop->AddEntry(hTopSel_source5, SourcesTopologicalSelection[4].c_str(), "l");
        legTop->AddEntry(hTopSel_source6, SourcesTopologicalSelection[5].c_str(), "l");
        legTop->AddEntry(hTopSel_total, "Total Systematic", "l");
        legTop->Draw();
        cTopSel_sys->SaveAs(Form("%s/fits/SystematicPlots/TopologicalSelect_TotalSys_%s.png", path.c_str(), canvasType.c_str()));

        // Combining all three systematic categories to get total systematic uncertainty
        TH1F *hTotalSys;
        vector<TH1F *> vec_allSys = {hSigExt_total, hTrkSel_total, hTopSel_total};
        QuadratureSum(vec_allSys, hTotalSys);
        TCanvas *cTotalSys = new TCanvas("", "", 720, 720);
        SetHistoQA(hSigExt_total);
        hSigExt_total->SetTitle(Form("%s", canvasType.c_str()));
        hSigExt_total->SetLineColor(lineColors[0]);
        hSigExt_total->SetMaximum(0.8);
        hSigExt_total->SetMinimum(0);
        hSigExt_total->Draw("HIST");
        hTrkSel_total->SetLineColor(lineColors[1]);
        hTrkSel_total->Draw("HIST SAME");
        hTopSel_total->SetLineColor(lineColors[2]);
        hTopSel_total->Draw("HIST SAME");
        hTotalSys->SetLineColor(kBlack);
        hTotalSys->SetLineWidth(3);
        hTotalSys->Draw("HIST SAME");
        TLegend *legTotal = new TLegend(0.55, 0.55, 0.85, 0.85);
        legTotal->AddEntry(hSigExt_total, "Signal Extraction", "l");
        legTotal->AddEntry(hTrkSel_total, "Track Selection", "l");
        legTotal->AddEntry(hTopSel_total, "Topological Selection", "l");
        legTotal->AddEntry(hTotalSys, "Total Systematic", "l");
        legTotal->Draw();
        cTotalSys->SaveAs(Form("%s/fits/SystematicPlots/Total_Systematics_%s.png", path.c_str(), canvasType.c_str()));
    }
}
//==========================================================================
//===========================End of main function===========================
//==========================================================================

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

void calculateRelativeUncertainty(TH1F *hDefault, TH1F *hVariation, TH1F *hRelUncertainty)
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

void checkBarlowSignificance(TH1F *hDefault, TH1F *hVariation, TH1F *hBarlow, string variationName, string type, bool &barlowPassed)
{
    int nBins = hDefault->GetNbinsX();
    for (int bin = 1; bin <= nBins; bin++)
    {
        float defaultValue = hDefault->GetBinContent(bin);
        float variationValue = hVariation->GetBinContent(bin);
        float defaultError = hDefault->GetBinError(bin);
        float variationError = hVariation->GetBinError(bin);

        float delta = variationValue - defaultValue;
        float sigma = std::sqrt(
            std::fabs(defaultError * defaultError - variationError * variationError));

        if (sigma > 0)
        {
            float barlowSignificance = delta / sigma;
            hBarlow->Fill(barlowSignificance);
        }
    }

    // ---- Barlow criteria ----
    double mean = hBarlow->GetMean();
    double rms = hBarlow->GetRMS();
    double total = hBarlow->Integral();

    double frac1 = (total > 0)
                       ? hBarlow->Integral(hBarlow->FindBin(-1),
                                           hBarlow->FindBin(1)) /
                             total
                       : 0;

    double frac2 = (total > 0)
                       ? hBarlow->Integral(hBarlow->FindBin(-2),
                                           hBarlow->FindBin(2)) /
                             total
                       : 0;

    bool cond1 = (std::abs(mean) < 0.1);
    bool cond2 = (rms < 1.1);
    bool cond3 = (frac1 > 0.60);
    bool cond4 = (frac2 > 0.90);

    int passed = cond1 + cond2 + cond3 + cond4;

    // ---- Print result ----
    if (passed >= 3)
    {
        // std::cout << "✅ Barlow check PASSED for set " << iset << std::endl;
        barlowPassed = true;
    }
    else
    {
        // std::cout << "❌ Barlow check FAILED for set " << iset << std::endl;
        barlowPassed = false;
    }

    // std::cout << "Mean = " << mean
    //           << ", RMS = " << rms << std::endl;
    // std::cout << "|n| < 1 fraction = " << frac1 * 100 << "%" << std::endl;
    // std::cout << "|n| < 2 fraction = " << frac2 * 100 << "%" << std::endl;
}

void processVariations(TFile **fVariations, const vector<string> &varNames, TH1F *defaultHist, TCanvas *canvas, TCanvas *canvasBarlow, vector<TH1F *> *resultVector, const string &histPath, const string &suffix, const string &path, const string &categoryName, const string &categoryType)
{

    for (int i = 0; i < varNames.size(); i++)
    {
        // Load variation histogram for this category
        TH1F *hist = (TH1F *)fVariations[i]->Get(histPath.c_str());
        if (hist == nullptr)
        {
            cout << "Error: Histogram " << histPath << " not found for variation: " << varNames[i] << endl;
            return;
        }

        // ==========Now check Barlow significance==============
        TH1F *hBarlow = new TH1F(Form("hBarlow_%s_%s", suffix.c_str(), varNames[i].c_str()), Form("Barlow Significance %s %s; n; Counts", suffix.c_str(), varNames[i].c_str()), 400, -100, 100);
        bool barlowPassed = true;
        checkBarlowSignificance(defaultHist, hist, hBarlow, varNames[i], suffix, barlowPassed);

        // Barlow criteria
        double total = hBarlow->Integral();
        double frac1 = (total > 0)
                           ? hBarlow->Integral(hBarlow->FindBin(-1), hBarlow->FindBin(1)) / total
                           : 0;

        double frac2 = (total > 0)
                           ? hBarlow->Integral(hBarlow->FindBin(-2), hBarlow->FindBin(2)) / total
                           : 0;

        canvasBarlow->cd(i + 1);
        hBarlow->Draw("HIST");
        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.06);
        if (barlowPassed)
        {
            lat.SetTextColor(kRed);
            lat.DrawLatex(0.17, 0.84, "Barlow PASSED");
        }
        else
        {
            lat.SetTextColor(kGreen + 2);
            lat.DrawLatex(0.17, 0.84, "Barlow Not PASSED");
        }
        lat.DrawLatex(0.17, 0.76, Form("Mean = %.3f", hBarlow->GetMean()));
        lat.DrawLatex(0.17, 0.68, Form("RMS = %.3f", hBarlow->GetRMS()));
        lat.DrawLatex(0.17, 0.60, Form("|n|<1 frac = %.2f%%", frac1 * 100));
        lat.DrawLatex(0.17, 0.52, Form("|n|<2 frac = %.2f%%", frac2 * 100));

        // Create relative uncertainty histogram and calculate
        TH1F *hRelUncert = (TH1F *)defaultHist->Clone(Form("hRelUncert_%s_%s", suffix.c_str(), varNames[i].c_str()));
        hRelUncert->SetTitle(Form("Variation: %s", varNames[i].c_str()));

        calculateRelativeUncertainty(defaultHist, hist, hRelUncert);
        resultVector->push_back(hRelUncert);

        // Draw on canvas
        canvas->cd(i + 1);
        hRelUncert->Draw("HIST");
    }
}

void averageQuadratureSum(const vector<TH1F *> &hVariationHistos, const int nVariations, TH1F *&hAvgQuadratureSum, const char *name)
{
    if (hVariationHistos.empty())
    {
        cerr << "Error: No histograms provided for averaging quadrature sum." << endl;
        return;
    }
    cout << "Category: " << name << endl;
    cout << "Total variations for average quadrature sum: " << nVariations << endl;
    cout << "Size of vector provided: " << hVariationHistos.size() << endl;

    int nBins = hVariationHistos[0]->GetNbinsX();
    hAvgQuadratureSum = (TH1F *)hVariationHistos[0]->Clone(name);
    hAvgQuadratureSum->Reset();
    // hAvgQuadratureSum->SetTitle(0);
    hAvgQuadratureSum->SetStats(0);
    hAvgQuadratureSum->SetLineWidth(2);

    for (int bin = 1; bin <= nBins; bin++)
    {
        double sumSquares = 0.0;
        for (const auto &hVar : hVariationHistos)
        {
            double value = hVar->GetBinContent(bin);
            sumSquares += value * value;
        }
        double avgSquare = sumSquares / nVariations;
        hAvgQuadratureSum->SetBinContent(bin, sqrt(avgSquare));
    }
}
void QuadratureSum(const vector<TH1F *> &hVariationHistos, TH1F *&hQuadratureSum)
{
    if (hVariationHistos.empty())
    {
        cerr << "Error: No histograms provided for quadrature sum." << endl;
        return;
    }
    cout << "Size of vector provided for total quadrature sum: " << hVariationHistos.size() << endl;

    int nBins = hVariationHistos[0]->GetNbinsX();
    hQuadratureSum = (TH1F *)hVariationHistos[0]->Clone("hQuadratureSum");
    hQuadratureSum->Reset();
    hQuadratureSum->SetTitle(0);
    hQuadratureSum->SetStats(0);
    hQuadratureSum->SetLineWidth(2);
    hQuadratureSum->SetLineColor(kBlack);

    for (int bin = 1; bin <= nBins; bin++)
    {
        double sumSquares = 0.0;
        for (const auto &hVar : hVariationHistos)
        {
            double value = hVar->GetBinContent(bin);
            sumSquares += value * value;
        }
        hQuadratureSum->SetBinContent(bin, sqrt(sumSquares));
    }
}
