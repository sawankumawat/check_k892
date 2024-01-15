#include "src/common.h"
#include "src/logger.h"
#include "src/DrawingHelper.C"
#include "src/SystematicHelper.cxx"
TLatex *t = new TLatex();
TH1 *GetSpectra(int sys, int iAnti, int centbin, TString sysExt = "", int sysvar = 0);
TH1 *ErrorSmoothing(TH1 *hinput, double smoothingCriteron = 1.5, int binstart = 2, int binend = kNPtBins-1);
TH1 *PlotSystematics_multi(int iAnti, int centbin, int correl_error_skip, bool DrawQAplots);
void Systematics()
{
    gErrorIgnoreLevel = kError; // Suppressing warning outputs
    TFile *outputFile = (TFile *)new TFile(kSystematicsOutput.data(), "recreate");
    LogCustom log(TLogLevel::INFO);
    log.Get(TLogLevel::INFO) << "Output file: " << kSystematicsOutput.data();
    t->SetNDC();
    t->SetTextSize(0.05);
    for (int iAnti = 0; iAnti < kAntiBinSize; iAnti++)
    {
        for (int iCent = 0; iCent < kNMultiplicityBins+1; iCent++)
        {
            auto test = (TH1 *)PlotSystematics_multi(iAnti,  // iAnti
                                                     iCent,  // centbin
                                                     0,      // correl_error_skip
                                                     kTRUE); // DrawQAplots
            outputFile->cd();
            test->Write(Form("hSystematic_%d_%d", iAnti, iCent));
            auto test_nocor =
                (TH1 *)PlotSystematics_multi(iAnti,   // iAnti
                                             iCent,   // centbin
                                             2,       // correl_error_skip
                                             kFALSE); // DrawQAplots
            outputFile->cd();
            test_nocor->Write(
                Form("hSystematic_NoCor_%d_%d", iAnti, iCent));
            auto test_cor =
                (TH1 *)PlotSystematics_multi(iAnti,   // iAnti
                                             iCent,   // centbin
                                             1,       // correl_error_skip
                                             kFALSE); // DrawQAplots
            outputFile->cd();
            test_cor->Write(
                Form("hSystematic_Cor_%d_%d", iAnti, iCent));
            auto test_MultiCor =
                (TH1 *)PlotSystematics_multi(iAnti,   // iAnti
                                             iCent,   // centbin
                                             3,       // correl_error_skip
                                             kFALSE); // DrawQAplots
            outputFile->cd();
            test_MultiCor->Write(
                Form("hSystematic_MultiCor_%d_%d", iAnti, iCent));
            auto test_MultiUnCor =
                (TH1 *)PlotSystematics_multi(iAnti,   // iAnti
                                             iCent,   // centbin
                                             4,       // correl_error_skip
                                             kFALSE); // DrawQAplots
            outputFile->cd();
            test_MultiUnCor->Write(
                Form("hSystematic_MultiUnCor_%d_%d", iAnti, iCent));
            if (kIsINEL) // INEL
                break;
        }
    }
}
TH1 *PlotSystematics_multi(int iAnti, int centbin, int correl_error_skip, bool DrawQAplots)
{
    // correl_error_skip = 0: total error (cor+uncor)
    // correl_error_skip = 1: cor error (cor)
    // correl_error_skip = 2: uncor error (uncor)
    // correl_error_skip = 3: Multicor error (cor)
    // correl_error_skip = 4: Multiuncor error (uncor)
    LogCustom log(TLogLevel::INFO);
    log.Get(TLogLevel::INFO) << "Checking systematics for: " << kParticleType[iAnti].data() << " in " << centbin << " multiplicity bin" << " with correl_error_skip = " << correl_error_skip;

    bool onlyCorrelatedError = (correl_error_skip == 1) || (correl_error_skip == 3);
    bool skipCorrelatedError = (correl_error_skip == 2) || (correl_error_skip == 4);
    int doMultCorCorrection = (correl_error_skip > 2);
    int isMultUnCor = (correl_error_skip == 4);
    bool smoothing = true;
    bool isINEL = false;
    bool use0100TopolError = true;
    bool setBarlow = false;
    bool useMBSigExtStudy = false;

    std::vector<TH1 *> totalsystematic;
    std::vector<TString> totalsystematicname;

    TString fvarfile;
    // Deafult Spectra
    // Given multiplicity
    auto hBase = (TH1 *)GetSpectra(1, iAnti, centbin);
    auto hBaseMB = (TH1 *)GetSpectra(1, iAnti, 0);
    
    TFile *fout_systematic_details = (correl_error_skip == 1)
                                         ? new TFile(Form("%s_cor_details_%d_%d.root", kSystematicsDetailOutput.data(), iAnti, centbin), "RECREATE")
                                     : (correl_error_skip == 2)
                                         ? new TFile(Form("%s_uncor_details_%d_%d.root", kSystematicsDetailOutput.data(), iAnti, centbin), "RECREATE")
                                     : (correl_error_skip == 3)
                                         ? new TFile(Form("%s_Multicor_details_%d_%d.root", kSystematicsDetailOutput.data(), iAnti, centbin), "RECREATE")
                                     : (correl_error_skip == 4)
                                         ? new TFile(Form("%s_Multiuncor_details_%d_%d.root", kSystematicsDetailOutput.data(), iAnti, centbin), "RECREATE")
                                         : new TFile(Form("%s_details_%d_%d.root", kSystematicsDetailOutput.data(), iAnti, centbin), "RECREATE");
    const TArrayD *ptbinarray = hBase->GetXaxis()->GetXbins();
    vector<double> ptbin;
    for (int i = 0; i < ptbinarray->GetSize(); i++)
        ptbin.push_back(ptbinarray->GetAt(i));
    
    // Signal Extraction variations
    // -----------------------------------------------------------
    std::vector<TH1 *> SigExtErrorHistos;
    // Type1 variation: No detailed study needed. ----------------
    std::vector<TH1 *> hSigExtSys_variations_type1; // spectra
    std::vector<TString> SigExtSys_type1bins = {"LSBkg", "pol3"};
    
    for (int i = 0; i < SigExtSys_type1bins.size(); i++)
    {
        auto hBkgSyst = (TH1 *)GetSpectra(1, iAnti, centbin, SigExtSys_type1bins[i], 0);
        hSigExtSys_variations_type1.push_back(hBkgSyst);
    }
    auto Systematic_SigExt = (useMBSigExtStudy)
                                 ? SystematicHelper(hBaseMB)
                                 : SystematicHelper(hBase);
    Systematic_SigExt.AddVariationtHistograms(hSigExtSys_variations_type1);
    Systematic_SigExt.SetBarlowCheck(0);
    Systematic_SigExt.SetDefaultXRange({kpTbin[0], kLastpTbin});
    // Systematic_SigExt.SetVerbose(); // for debuging
    Systematic_SigExt.SetDefaultYRange({0.5, 1.5});
    Systematic_SigExt.SetVariationNames(SigExtSys_type1bins);
    Systematic_SigExt.SetYaxisTitle("Ratio");
    Systematic_SigExt.InitAbsDiffRatioColors();

    auto hSigExtErrors = Systematic_SigExt.GetAbsDiffRatio();
    auto hSigExtErrorMax = Systematic_SigExt.GetMaxAbsDiffRatio();
    auto hSigExtErrorSum = Systematic_SigExt.GetSumAbsDiffRatio();
    auto hSigExtErrorStd = Systematic_SigExt.GetStdevRatio();

    std::vector<double> SigExtErrorSum_variation1;
    for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
    {
        SigExtErrorSum_variation1.push_back(0);
        for (auto const &hInput : hSigExtErrors)
        {
            if (abs(hInput->GetBinContent(bin + 1)) > abs(SigExtErrorSum_variation1[bin]))
                SigExtErrorSum_variation1[bin] = abs(hInput->GetBinContent(bin + 1));
        }
    }
    auto hSigExtTotalError_Type1 =
        (TH1 *)MakeHistfromArray("SigExtErrorSum1", SigExtErrorSum_variation1, kpTbin);
    
    // Fit systematics
    std::vector<TH1 *> hSignalExtsys_variation;
    std::vector<TString> SignalExtSystematic_bins = {"fitrangeL", "fitrangeR", "fitrangeLR", "normrangeL", "normrangeR", "normrangeLR"};
    std::vector<TString> SignalExtSystematic_bins_used;
    int sysvarrange = 2;
    for (int signalExt_cutvar = 0;
         signalExt_cutvar < SignalExtSystematic_bins.size(); signalExt_cutvar++)
    {
        for (int signalExt_cutvarrange = -sysvarrange;
             signalExt_cutvarrange < sysvarrange + 1; signalExt_cutvarrange++)
        {
            if (signalExt_cutvarrange == 0)
                continue;
            auto hFitSyst = (TH1 *)GetSpectra(1, iAnti, centbin, SignalExtSystematic_bins[signalExt_cutvar], signalExt_cutvarrange);
            hSignalExtsys_variation.push_back(hFitSyst);
            SignalExtSystematic_bins_used.push_back(Form("%s_%d", SignalExtSystematic_bins[signalExt_cutvar].Data(), signalExt_cutvarrange));
        }
    }

    auto Systematic_SysExt = SystematicHelper(hBase);
    Systematic_SysExt.AddVariationtHistograms(hSignalExtsys_variation);
    // Systematic_SysExt.SetVerbose(); // for debuging
    Systematic_SysExt.SetBarlowCheck(0); // balrow cut is not applied in signal extraction part
    Systematic_SysExt.SetDefaultXRange({kpTbin[0], kLastpTbin});
    Systematic_SysExt.SetDefaultYRange({0.1, 1.2});
    Systematic_SysExt.SetYaxisTitle("Ratio");
    Systematic_SysExt.InitAbsDiffRatioColors(false);
    Systematic_SysExt.SetVariationNames(SignalExtSystematic_bins_used);
    Systematic_SysExt.InitVariationColors();

    auto hSysExtErrors = Systematic_SysExt.GetAbsDiffRatio(); // use sum
    auto hSysExtErrorsMax = Systematic_SysExt.GetMaxAbsDiffRatio();
    auto hSysExtErrorSum = Systematic_SysExt.GetSumAbsDiffRatio();
    auto hSysExtErrorStd = Systematic_SysExt.GetStdevRatio();

    // Error sum
    std::vector<double> SigExtErrorSum;
    for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
    {
        SigExtErrorSum.push_back(0);
        for (auto const &hInput : hSysExtErrors)
        {
            if (abs(hInput->GetBinContent(bin + 1)) > abs(SigExtErrorSum[bin]))
                SigExtErrorSum[bin] = abs(hInput->GetBinContent(bin + 1));
        }
    }
    auto hSigExtTotalError_Type2 =
        (TH1 *)MakeHistfromArray("SigExtErrorSum2", SigExtErrorSum, kpTbin);
    hSigExtErrorStd->SetName("hSigExtErrorStd");
    hSysExtErrorStd->SetName("hSysExtErrorStd");

    // TrackSelection Cut variations
    // -----------------------------------------------------------
    std::vector<TH1 *> hTrksys_variation;
    std::vector<TH1 *> hTrksys_variation_selected;
    std::vector<TString> TrackSelectionCutSystematic_bins = {"_PVTrk", "_Trk_loose"};
    std::vector<TString> TrackSelectionCutSystematic_bins_selected = {"_PVTrk", "_Trk_loose"};
    
    for (int TrkSel_cutvar = 0;
         TrkSel_cutvar < TrackSelectionCutSystematic_bins.size();
         TrkSel_cutvar++)
    {
        auto it = std::find(SystematicBinsName.begin(), SystematicBinsName.end(),
                            TrackSelectionCutSystematic_bins[TrkSel_cutvar]);
        auto temp = (TH1 *)GetSpectra(std::distance(SystematicBinsName.begin(), it), iAnti, centbin);
        hTrksys_variation.push_back(temp);
    }
    for (int TrkSel_cutvar = 0;
         TrkSel_cutvar < TrackSelectionCutSystematic_bins_selected.size();
         TrkSel_cutvar++)
    {
        auto it = std::find(SystematicBinsName.begin(), SystematicBinsName.end(),
                            TrackSelectionCutSystematic_bins_selected[TrkSel_cutvar]);
        auto temp = (TH1 *)GetSpectra(std::distance(SystematicBinsName.begin(), it), iAnti, centbin);
        hTrksys_variation_selected.push_back(temp);
    }
    auto Systematic_TrkSelCuts = SystematicHelper(hBase);
    Systematic_TrkSelCuts.AddVariationtHistograms(hTrksys_variation);
    // Systematic_TrkSelCuts.SetVerbose(); // for debuging
    Systematic_TrkSelCuts.SetBarlowCheck(setBarlow); // balrow cut
    Systematic_TrkSelCuts.SetDefaultXRange({kpTbin[0], kLastpTbin});
    Systematic_TrkSelCuts.SetDefaultYRange({0.2, 0.2});
    Systematic_TrkSelCuts.SetYaxisTitle("Ratio");
    Systematic_TrkSelCuts.InitAbsDiffRatioColors(false);
    // Error sum for the proper error selection
    std::vector<double> TrkSelErrorSum;
    for (auto const &hInput : hTrksys_variation)
    {
        auto tempvalue = 0.0;
        for (int bin = 0; bin < hInput->GetNbinsX(); bin++)
        {
            tempvalue += abs(1 - hInput->GetBinContent(bin + 1) / hBase->GetBinContent(bin + 1));
        }
        TrkSelErrorSum.push_back(tempvalue);
    }
    std::vector<TString> TrackSelectionCutSystematic_bins_with_errorsum;
    for (int index = 0; index < TrackSelectionCutSystematic_bins.size(); index++)
    {
        TString temp("");
        temp.Form("%s - Error sum: %.2f%%", TrackSelectionCutSystematic_bins[index].Data(), TrkSelErrorSum[index] * 100);
        TrackSelectionCutSystematic_bins_with_errorsum.push_back(temp);
    }
    Systematic_TrkSelCuts.SetVariationNames(TrackSelectionCutSystematic_bins_with_errorsum);
    Systematic_TrkSelCuts.InitVariationColors();

    auto hTrkSelErrors = Systematic_TrkSelCuts.GetAbsDiffRatio(); // use sum
    auto hTrkSelErrorsMax = Systematic_TrkSelCuts.GetMaxAbsDiffRatio();
    auto hTrkSelErrorSum =
        Systematic_TrkSelCuts.GetSumAbsDiffRatio();

    // Selected TrkSel cuts
    auto Systematic_TrkSelCuts_selected = SystematicHelper(hBase);
    Systematic_TrkSelCuts_selected.AddVariationtHistograms(hTrksys_variation_selected);
    // Systematic_TrkSelCuts_selected.SetVerbose(); // for debuging
    Systematic_TrkSelCuts_selected.SetBarlowCheck(setBarlow); // balrow cut
    Systematic_TrkSelCuts_selected.SetDefaultXRange({kpTbin[0], kLastpTbin});
    Systematic_TrkSelCuts_selected.SetDefaultYRange({0.2, 0.2});
    Systematic_TrkSelCuts_selected.SetYaxisTitle("Ratio");
    Systematic_TrkSelCuts_selected.InitAbsDiffRatioColors(false);
    Systematic_TrkSelCuts_selected.SetVariationNames(TrackSelectionCutSystematic_bins_selected);
    Systematic_TrkSelCuts_selected.InitVariationColors();

    auto hTrkSelErrors_selected = Systematic_TrkSelCuts_selected.GetStdevRatio(); // use sum
    auto hTrkSelErrorsMax_selected = Systematic_TrkSelCuts_selected.GetMaxAbsDiffRatio();
    auto hTrkSelErrorSum_selected =
        Systematic_TrkSelCuts_selected.GetSumAbsDiffRatio();
    hTrkSelErrors_selected->SetName("hTrkSelErrorSum");

    // Topological Cut variations
    // -----------------------------------------------------------
    std::vector<TH1 *> hTopolsys_variation;
    std::vector<TH1 *> hTopolsys_variation_selected;
    std::vector<TString> TopologicalCutSystematic_bins = {
        "_DCAxy_loose", "_DCAz_loose", "_DCA_loose"};
    std::vector<TString> TopologicalCutSystematic_bins_selected = {
        "_DCAxy_loose", "_DCAz_loose", "_DCA_loose"};
    
    for (int topological_cutvar = 0;
         topological_cutvar < TopologicalCutSystematic_bins.size();
         topological_cutvar++)
    {
        auto it = std::find(SystematicBinsName.begin(), SystematicBinsName.end(),
                            TopologicalCutSystematic_bins[topological_cutvar]);
        auto temp = (TH1 *)GetSpectra(std::distance(SystematicBinsName.begin(), it), iAnti, centbin);
        hTopolsys_variation.push_back(temp);
    }
    for (int topological_cutvar = 0;
         topological_cutvar < TopologicalCutSystematic_bins_selected.size();
         topological_cutvar++)
    {
        auto it = std::find(SystematicBinsName.begin(), SystematicBinsName.end(),
                            TopologicalCutSystematic_bins_selected[topological_cutvar]);
        auto temp = (TH1 *)GetSpectra(std::distance(SystematicBinsName.begin(), it), iAnti, centbin);
        hTopolsys_variation_selected.push_back(temp);
    }
    auto Systematic_TopolCuts = SystematicHelper(hBase);
    Systematic_TopolCuts.AddVariationtHistograms(hTopolsys_variation);
    // Systematic_TopolCuts.SetVerbose(); // for debuging
    Systematic_TopolCuts.SetBarlowCheck(setBarlow); // balrow cut
    Systematic_TopolCuts.SetDefaultXRange({kpTbin[0], kLastpTbin});
    Systematic_TopolCuts.SetDefaultYRange({0.2, 0.2});
    Systematic_TopolCuts.SetYaxisTitle("Ratio");
    Systematic_TopolCuts.InitAbsDiffRatioColors(false);
    // Error sum for the proper error selection
    std::vector<double> TopolErrorSum;
    for (auto const &hInput : hTopolsys_variation)
    {
        auto tempvalue = 0.0;
        for (int bin = 0; bin < hInput->GetNbinsX(); bin++)
        {
            tempvalue += abs(1 - hInput->GetBinContent(bin + 1) / hBase->GetBinContent(bin + 1));
        }
        TopolErrorSum.push_back(tempvalue);
    }
    std::vector<TString> TopologicalCutSystematic_bins_with_errorsum;
    for (int index = 0; index < TopologicalCutSystematic_bins.size(); index++)
    {
        TString temp("");
        temp.Form("%s - Error sum: %.2f%%", TopologicalCutSystematic_bins[index].Data(), TopolErrorSum[index] * 100);
        TopologicalCutSystematic_bins_with_errorsum.push_back(temp);
    }
    Systematic_TopolCuts.SetVariationNames(TopologicalCutSystematic_bins_with_errorsum);
    Systematic_TopolCuts.InitVariationColors();

    auto hTopolCutErrors = Systematic_TopolCuts.GetSumAbsDiffRatio(); // use sum
    auto hTopolCutErrorsMax = Systematic_TopolCuts.GetMaxAbsDiffRatio();
    auto hTopolCutErrorSum =
        Systematic_TopolCuts.GetSumAbsDiffRatio();

    // Selected Topol cuts
    auto Systematic_TopolCuts_selected = SystematicHelper(hBase);
    Systematic_TopolCuts_selected.AddVariationtHistograms(hTopolsys_variation_selected);
    // Systematic_TopolCuts_selected.SetVerbose(); // for debuging
    Systematic_TopolCuts_selected.SetBarlowCheck(setBarlow); // balrow cut
    Systematic_TopolCuts_selected.SetDefaultXRange({kpTbin[0], kLastpTbin});
    Systematic_TopolCuts_selected.SetDefaultYRange({0.2, 0.2});
    Systematic_TopolCuts_selected.SetYaxisTitle("Ratio");
    Systematic_TopolCuts_selected.InitAbsDiffRatioColors(false);
    Systematic_TopolCuts_selected.SetVariationNames(TopologicalCutSystematic_bins_selected);
    Systematic_TopolCuts_selected.InitVariationColors();

    auto hTopolCutErrors_selected = Systematic_TopolCuts_selected.GetSumAbsDiffRatio(); // use sum
    auto hTopolCutErrorsMax_selected = Systematic_TopolCuts_selected.GetMaxAbsDiffRatio();
    auto hTopolCutErrorSum_selected =
        Systematic_TopolCuts_selected.GetSumAbsDiffRatio();
    hTopolCutErrors_selected->SetName("hTopolCutErrorSum");

    // PID Cut variations
    // -----------------------------------------------------------
    std::vector<TH1 *> hPIDCutsys_variation;
    std::vector<TString> PIDCutSystematic_bins = {"_PID20", "_PID25", "_PID35", "_PID40"};
    
    for (int topological_cutvar = 0;
         topological_cutvar < PIDCutSystematic_bins.size();
         topological_cutvar++)
    {
        auto it = std::find(SystematicBinsName.begin(), SystematicBinsName.end(), PIDCutSystematic_bins[topological_cutvar]);
        auto temp = (TH1 *)GetSpectra(std::distance(SystematicBinsName.begin(), it), iAnti, centbin);
        hPIDCutsys_variation.push_back(temp);
    }

    auto Systematic_PIDCuts = SystematicHelper(hBase);
    Systematic_PIDCuts.AddVariationtHistograms(hPIDCutsys_variation);
    // Systematic_PIDCuts.SetVerbose(); // for debuging
    Systematic_PIDCuts.SetBarlowCheck(setBarlow); // balrow cut
    Systematic_PIDCuts.SetDefaultXRange({kpTbin[0], kLastpTbin});
    Systematic_PIDCuts.SetDefaultYRange({0.2, 0.2});
    Systematic_PIDCuts.SetYaxisTitle("Ratio");
    Systematic_PIDCuts.InitAbsDiffRatioColors(false);
    Systematic_PIDCuts.SetVariationNames(PIDCutSystematic_bins);
    Systematic_PIDCuts.InitVariationColors();

    auto hPIDCutErrors = Systematic_PIDCuts.GetAbsDiffRatio();
    // auto hPIDCutErrorMax = Systematic_PIDCuts.GetMaxAbsDiffRatio();
    // auto hPIDCutErrorSum = Systematic_PIDCuts.GetSumAbsDiffRatio();

    // Error sum
    std::vector<double> PIDCutErrorSum;
    for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
    {
        PIDCutErrorSum.push_back(0);
        for (auto const &hInput : hPIDCutsys_variation)
        {
            if (abs(hInput->GetBinContent(bin + 1)) > abs(PIDCutErrorSum[bin]))
                PIDCutErrorSum[bin] = abs(hInput->GetBinContent(bin + 1));
        }
    }
    auto hPIDCutErrorMax =
        (TH1 *)MakeHistfromArray("PIDCutErrorSum", PIDCutErrorSum, kpTbin);
    hPIDCutErrorMax->SetName("hPIDCutErrorMax");

    // Other Systematics
    // -----------------------------------------------------------
    // 1. Material Budget (Kaon + pion)
    // 
    // Further study is needed
    // FIXME: need to be updated
    std::vector<double> materialbudget;
    for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
    {
        materialbudget.push_back(kMaterialBudgetError);
    }
    auto hSysMaterialBudget =
        (TH1 *)MakeHistfromArray("hSysMaterialBudget", materialbudget, kpTbin);
    // 2. Mult. indipendent efficiencies
    
    // Further study is needed
    // Currently multiplicity-dependent efficiency has not been studied
    vector<double> multiindepeffi;
    std::vector<double> values = {0.047497416423538, 0.038637794134421734, 0.039450601658603326, 0.02417357921212393, 0.047592353968654785, 0.039459592736600664, 0.04602622819499959, 0.02838256212691799, 0.034094916124434195, 0.014598486425621598, 0.028644755489325713, 0.02129979811198308, 0.03542126946191168};
    for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
        multiindepeffi.push_back(values[bin]);
    auto hSysMultiIndepEffi =
        (TH1 *)MakeHistfromArray("hSysMultiIndepEffi", multiindepeffi, kpTbin);
    // 3. Tracking efficiencies
    // Unokwn
    // vector<double> trackingeffi;
    // for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
    //     trackingeffi.push_back(0.03);
    // auto hSysTrackingEffi =
    //     (TH1 *)MakeHistfromArray("hSysTrackingEffi", trackingeffi, kpTbin);
    
    // R values for multiplicity correlated uncertatainty
    // -----------------------------------------------------------
    if (fout_systematic_details->IsOpen() && doMultCorCorrection && centbin > 0)
    {
        // MB spectra and errors
        TFile *fout_systematic0100 = new TFile(Form("%s_Multiuncor_details_%d_0.root", kSystematicsDetailOutput.data(), iAnti), "READ");
        TH1 *hSigExtErrorStd0100 = (TH1 *)fout_systematic0100->Get("hSigExtErrorStd");           //
        TH1 *hSysExtErrorStd0100 = (TH1 *)fout_systematic0100->Get("hSysExtErrorStd");           //
        TH1 *hTrkSelErrors_selected0100 = (TH1 *)fout_systematic0100->Get("hTrkSelErrorSum"); //
        TH1 *hTopolCutErrorSum0100 = (TH1 *)fout_systematic0100->Get("hTopolCutErrorSum");       //
        TH1 *hPIDCutErrorMax0100 = (TH1 *)fout_systematic0100->Get("hPIDCutErrorMax");           //
        vector<double> listRvaluesForMutiCorUnc_SigExt;
        vector<double> listRvaluesForMutiCorUnc_SysExt;
        vector<double> listRvaluesForMutiCorUnc_TrkSel;
        vector<double> listRvaluesForMutiCorUnc_Topol;
        vector<double> listRvaluesForMutiCorUnc_PID;

        for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
        {
            listRvaluesForMutiCorUnc_SigExt.push_back((1 + hSigExtErrorStd->GetBinContent(bin + 1)) / (1 + hSigExtErrorStd0100->GetBinContent(bin + 1)));
            listRvaluesForMutiCorUnc_SysExt.push_back((1 + hSysExtErrorStd->GetBinContent(bin + 1)) / (1 + hSysExtErrorStd0100->GetBinContent(bin + 1)));
            listRvaluesForMutiCorUnc_TrkSel.push_back((1 + hTrkSelErrors_selected->GetBinContent(bin + 1)) / (1 + hTrkSelErrors_selected0100->GetBinContent(bin + 1)));
            listRvaluesForMutiCorUnc_Topol.push_back((1 + hTopolCutErrors_selected->GetBinContent(bin + 1)) / (1 + hTopolCutErrorSum0100->GetBinContent(bin + 1)));
            listRvaluesForMutiCorUnc_PID.push_back((1 + hPIDCutErrorMax->GetBinContent(bin + 1)) / (1 + hPIDCutErrorMax0100->GetBinContent(bin + 1)));
        }
        // for (auto i : listRvaluesForMutiCorUnc_SigExt)
        //     std::cout << i << std::endl;
        auto hRvaluesForMutiCorUnc_SigExt =
            (TH1 *)MakeHistfromArray("RvaluesForMutiCorUnc_SigExt", listRvaluesForMutiCorUnc_SigExt, ptbin);
        auto hRvaluesForMutiCorUnc_SysExt =
            (TH1 *)MakeHistfromArray("RvaluesForMutiCorUnc_SysExt", listRvaluesForMutiCorUnc_SysExt, ptbin);
        auto hRvaluesForMutiCorUnc_TrkSel =
            (TH1 *)MakeHistfromArray("RvaluesForMutiCorUnc_TrkSel", listRvaluesForMutiCorUnc_TrkSel, ptbin);
        auto hRvaluesForMutiCorUnc_Topol =
            (TH1 *)MakeHistfromArray("RvaluesForMutiCorUnc_Topol", listRvaluesForMutiCorUnc_Topol, ptbin);
        auto hRvaluesForMutiCorUnc_PID =
            (TH1 *)MakeHistfromArray("RvaluesForMutiCorUnc_PID", listRvaluesForMutiCorUnc_PID, ptbin);
        // Multiplicity correlated uncertainties
        for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
        {
            if (isMultUnCor)
            {
                hSigExtErrorStd->SetBinContent(bin + 1, hSigExtErrorStd->GetBinContent(bin + 1) * (TMath::Abs(1 - listRvaluesForMutiCorUnc_SigExt[bin])));
                hSysExtErrorStd->SetBinContent(bin + 1, hSysExtErrorStd->GetBinContent(bin + 1) * (TMath::Abs(1 - listRvaluesForMutiCorUnc_SysExt[bin])));
                hTrkSelErrors_selected->SetBinContent(bin + 1, hTrkSelErrors_selected->GetBinContent(bin + 1) * (TMath::Abs(1 - listRvaluesForMutiCorUnc_TrkSel[bin])));
                hTopolCutErrors_selected->SetBinContent(bin + 1, hTopolCutErrors_selected->GetBinContent(bin + 1) * (TMath::Abs(1 - listRvaluesForMutiCorUnc_Topol[bin])));
                hPIDCutErrorMax->SetBinContent(bin + 1, hPIDCutErrorMax->GetBinContent(bin + 1) * (TMath::Abs(1 - listRvaluesForMutiCorUnc_PID[bin])));
            }
            else
            {
                hSigExtErrorStd->SetBinContent(bin + 1, hSigExtErrorStd->GetBinContent(bin + 1) * (1 - TMath::Abs(1 - listRvaluesForMutiCorUnc_SigExt[bin])));
                hSysExtErrorStd->SetBinContent(bin + 1, hSysExtErrorStd->GetBinContent(bin + 1) * (1 - TMath::Abs(1 - listRvaluesForMutiCorUnc_SysExt[bin])));
                hTrkSelErrors_selected->SetBinContent(bin + 1, hTrkSelErrors_selected->GetBinContent(bin + 1) * (1 - TMath::Abs(1 - listRvaluesForMutiCorUnc_TrkSel[bin])));
                hTopolCutErrors_selected->SetBinContent(bin + 1, hTopolCutErrors_selected->GetBinContent(bin + 1) * (1 - TMath::Abs(1 - listRvaluesForMutiCorUnc_Topol[bin])));
                hPIDCutErrorMax->SetBinContent(bin + 1, hPIDCutErrorMax->GetBinContent(bin + 1) * (1 - TMath::Abs(1 - listRvaluesForMutiCorUnc_PID[bin])));
            }
        }
        // if (fout_systematic_details->IsOpen())
        // {
        //     fout_systematic_details->cd();
        //     hRvaluesForMutiCorUnc_BkgFunc->Write();
        //     hRvaluesForMutiCorUnc_SigExt->Write();
        //     hRvaluesForMutiCorUnc_SysExt->Write();
        //     hRvaluesForMutiCorUnc_Topol->Write();
        //     hRvaluesForMutiCorUnc_PID->Write();
        // }
    }
    // -----------------------------------------------------------


    // Total systematic uncertainty
    // -----------------------------------------------------------//
    if (!onlyCorrelatedError)
    {
        // totalsystematic.push_back(hSigExtTotalError_Type1); // max
        totalsystematic.push_back(hSigExtErrorStd); // std
        totalsystematicname.push_back("Signal Extraction method variation");
        // totalsystematic.push_back(hSigExtTotalError_Type2); // max
        totalsystematic.push_back(hSysExtErrorStd); // std
        totalsystematicname.push_back("Signal Extraction range variation");
        totalsystematic.push_back(hTrkSelErrors_selected); // sum
        totalsystematicname.push_back("Track Selection variation");
        totalsystematic.push_back(hTopolCutErrors_selected);
        totalsystematicname.push_back("Topological Cut variation");
        totalsystematic.push_back(hPIDCutErrorMax);
        totalsystematicname.push_back("TPC PID Cut variation");
        if (!skipCorrelatedError)
        {
            totalsystematic.push_back(hSysMaterialBudget);
            totalsystematicname.push_back("Material budget");
            totalsystematic.push_back(hSysMultiIndepEffi);
            totalsystematicname.push_back("MC correction");
            // totalsystematic.push_back(hSysTrackingEffi);
            // totalsystematicname.push_back("Tracking Efficiency");
            // totalsystematic.push_back(hSysbranchingRatio);
            // totalsystematicname.push_back("Branching Ratio");
        }
    }
    else
    {
        totalsystematic.push_back(hSysMaterialBudget);
            totalsystematicname.push_back("Material budget");
        totalsystematic.push_back(hSysMultiIndepEffi);
        totalsystematicname.push_back("MC correction");
        // totalsystematic.push_back(hSysTrackingEffi);
        // totalsystematicname.push_back("Tracking Efficiency");
        // totalsystematic.push_back(hSysbranchingRatio);
        // totalsystematicname.push_back("Branching Ratio");
    }

    // Smoothing
    std::vector<TH1 *> totalsystematic_beforesmoothing;
    if (smoothing)
    {
        for (int sysbin = 0; sysbin < totalsystematic.size(); sysbin++)
        {
            // keep old syserror
            auto tempOldSyst = (TH1 *)totalsystematic[sysbin]->Clone();
            totalsystematic_beforesmoothing.push_back(tempOldSyst);
            // smoothed error
            totalsystematic[sysbin] = ErrorSmoothing(totalsystematic[sysbin], 1.5);
            // std::cout << "[Smoothing] before: " << totalsystematic[sysbin]->GetBinContent(9) << ", after: " << totalsystematic[sysbin]->GetBinContent(9) << std::endl;
        }
    }

    std::vector<double> totalSysError; // Total sys error sum array
    for (int bin = 0; bin < totalsystematic[0]->GetNbinsX(); bin++)
    {
        double TotError = 0.;
        for (int sysbin = 0; sysbin < totalsystematic.size(); sysbin++)
            TotError += pow(totalsystematic[sysbin]->GetBinContent(bin + 1), 2);
        totalSysError.push_back(sqrt(TotError));
    }
    auto hSysSpectra = (TH1 *)hBase->Clone(); // Spectrum with systematic error
    for (int bin = 0; bin < hSysSpectra->GetNbinsX(); bin++)
        hSysSpectra->SetBinError(bin + 1, totalSysError[bin] * hSysSpectra->GetBinContent(bin + 1));
    // SAVE
    if (fout_systematic_details->IsOpen())
    {
        fout_systematic_details->cd();
        if (!onlyCorrelatedError)
        {
            hSigExtErrorStd->Write();
            hSysExtErrorStd->Write();
            hTrkSelErrors_selected->Write();
            hTopolCutErrors_selected->Write();
            hPIDCutErrorMax->Write();
            if (!skipCorrelatedError)
            {
                hSysMaterialBudget->Write();
                hSysMultiIndepEffi->Write();
                // hSysTrackingEffi->Write();
                // hSysbranchingRatio->Write();
            }
        }
        else
        {
            hSysMaterialBudget->Write();
            hSysMultiIndepEffi->Write();
            // hSysTrackingEffi->Write();
            // hSysbranchingRatio->Write();
        }
    }

    // QA Plots
    // -----------------------------------------------------------
    if (DrawQAplots)
    {
        auto hSysTotalSum = // systematic error fraction
            (TH1 *)MakeHistfromArray(Form("hSysTotalSum_%d", centbin), totalSysError, kpTbin);
        hSysTotalSum->SetLineColor(kBlack);
        hSysTotalSum->SetLineStyle(2);
        hSysTotalSum->SetLineWidth(2);
        hSysTotalSum->SetName("hSysTotalSum");

        std::vector<double> baseStatError = GetBinErrorFrations(hBase);
        auto hStatError = // statistical error fraction
            (TH1 *)MakeHistfromArray(Form("hSysTotalSum_%d", centbin), baseStatError, kpTbin);
        hStatError->SetFillColorAlpha(kBlack, 0.2);
        hStatError->SetLineColorAlpha(kBlack, 0.2);
        hStatError->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hStatError->SetMaximum(0.3);
        hStatError->SetMinimum(0);
        hStatError->GetXaxis()->SetRangeUser(kpTbin[0], kLastpTbin);

        // Initialize the total systematic errors
        // color, axis range
        int sysColorPallet = GetSerialColors(totalsystematic.size());
        for (int sysbin = 0; sysbin < totalsystematic.size(); sysbin++)
        {
            totalsystematic[sysbin]->SetMarkerColor(
                sysColorPallet + totalsystematic.size() - sysbin - 1);
            totalsystematic[sysbin]->SetLineColor(
                sysColorPallet + totalsystematic.size() - sysbin - 1);
            totalsystematic[sysbin]->GetXaxis()->SetRangeUser(kpTbin[0], kLastpTbin);
            totalsystematic[sysbin]->SetTitle("");
            totalsystematic[sysbin]->SetLineWidth(2);
        }
        TCanvas *cQA = new TCanvas(Form("cFitResults_centbin_%d", centbin), Form("cFitResults_centbin_%d", centbin), 960, 720);
        TGaxis::SetMaxDigits(3);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetLegendBorderSize(0);
        cQA->SetTickx();
        cQA->SetTicky();
        cQA->SetTopMargin(0.05);
        cQA->SetLeftMargin(0.10);
        // cQA->SetBottomMargin(0.01);
        cQA->SetRightMargin(0.01);
        cQA->SetFillStyle(0);
        cQA->Draw();
        cQA->cd();
        auto legend_TotalSys = new TLegend(.6, .63, .9, .88);
        legend_TotalSys->SetBorderSize(0);
        legend_TotalSys->SetFillStyle(0);
        legend_TotalSys->AddEntry(hStatError, "Statistical Error", "F");
        hStatError->Draw("BAR");
        for (int i = 0; i < totalsystematic.size(); i++)
        {
            legend_TotalSys->AddEntry(totalsystematic[i], totalsystematicname[i], "L");
            totalsystematic[i]->Draw("hist same");
        }
        legend_TotalSys->AddEntry(hSysTotalSum, "Total systematic uncertainty", "L");
        hSysTotalSum->Draw("hist same");
        legend_TotalSys->Draw();
        if (!centbin)
            t->DrawLatex(0.13, 0.88, Form("%s MB", kParticleType[iAnti].c_str()));
        else
            t->DrawLatex(0.13, 0.88, Form("%s %.1f-%.1f", kParticleType[iAnti].c_str(), kMultiplicityBins[centbin - 1], kMultiplicityBins[centbin]));
        SaveCanvas(cQA, Form("cQA_Fullsystematic_%d_%d", iAnti, centbin), Form("%ssys/Full/", kFiguresFolder.data()), kOutputType.data());
        // Smoothing Check
        if (smoothing)
        {
            vector<double> totalSysError_noSmoothing;
            for (int bin = 0; bin < totalsystematic[0]->GetNbinsX(); bin++)
            {
                double TotError = 0.;
                for (int sysbin = 0; sysbin < totalsystematic.size(); sysbin++)
                {
                    TotError += pow(totalsystematic_beforesmoothing[sysbin]->GetBinContent(bin + 1),2);
                }
                totalSysError_noSmoothing.push_back(sqrt(TotError));
            }

            // Initialize the total systematic errors
            // color, axis range
            for (int sysbin = 0; sysbin < totalsystematic.size(); sysbin++)
            {
                totalsystematic_beforesmoothing[sysbin]->SetMarkerColor(
                    sysColorPallet + totalsystematic.size() - sysbin - 1);
                totalsystematic_beforesmoothing[sysbin]->SetLineColor(
                    sysColorPallet + totalsystematic.size() - sysbin - 1);
                totalsystematic_beforesmoothing[sysbin]
                    ->GetXaxis()
                    ->SetRangeUser(kpTbin[0], kLastpTbin);
                totalsystematic_beforesmoothing[sysbin]->SetTitle("");
                totalsystematic_beforesmoothing[sysbin]->SetLineWidth(2);
            }

            auto hSysTotalSum_noSmoothing = // systematic error fraction
                (TH1 *)MakeHistfromArray(
                    Form("hSysTotalSum_%d", centbin),
                    totalSysError_noSmoothing, kpTbin);
            hSysTotalSum_noSmoothing->SetLineColor(kBlack);
            hSysTotalSum_noSmoothing->SetMaximum(0.3);
            hSysTotalSum_noSmoothing->SetLineStyle(2);
            hSysTotalSum_noSmoothing->SetLineWidth(2);
            cQA->cd();
            hStatError->Draw("BAR");
            for (int i = 0; i < totalsystematic_beforesmoothing.size(); i++)
            {
                totalsystematic_beforesmoothing[i]->Draw("hist same");
            }
            hSysTotalSum_noSmoothing->Draw("hist same");
            legend_TotalSys->Draw();
            TString NoSmoothtext;
            if (kIsINEL)
                NoSmoothtext = "cQA_TotalError_INEL_noSmoothing";
            else
                NoSmoothtext = Form("cQA_TotalError_%d_noSmoothing", centbin);
            if (setBarlow)
                NoSmoothtext += "_Balrow";
            SaveCanvas(cQA, NoSmoothtext.Data(), Form("%ssys/NoSmoothing/", kFiguresFolder.data()), kOutputType.data());
        }

        // QA Hist for Signal Extraion(Ratio)
        Systematic_SysExt.SetDefaultYmin(0.8);
        auto hSysExtErrors_Ratio = Systematic_SysExt.GetRatio();
        Systematic_SysExt.InitRatioColors();
        //
        Systematic_SigExt.SetDefaultYmin(0.8);
        auto hSigExtErrors_Ratio = Systematic_SigExt.GetRatio();
        Systematic_SigExt.InitRatioColors();
        Systematic_SigExt.InitVariationColors();
        //
        Systematic_TrkSelCuts_selected.SetDefaultYmax(0.05);
        auto hTrkSelCutsErrors_Ratio = Systematic_TrkSelCuts_selected.GetRatio();
        Systematic_TrkSelCuts_selected.InitRatioColors();

        SaveCanvas(Systematic_SigExt.GetQAPlot("StdevRatio"),
                   Form("cQA_SigExt_%d_%d", iAnti, centbin),
                   Form("%ssys/SigExt/", kFiguresFolder.data()), kOutputType.data());
        SaveCanvas(Systematic_SysExt.GetQAPlot("StdevRatio"),
                   Form("cQA_SysExt_%d_%d", iAnti, centbin),
                   Form("%ssys/SysExt/", kFiguresFolder.data()), kOutputType.data());
        SaveCanvas(Systematic_TrkSelCuts_selected.GetQAPlot("AbsDiffRatio"),
                   Form("cQA_TrkSelCuts_%d_%d", iAnti, centbin),
                   Form("%ssys/TrkSelCuts/", kFiguresFolder.data()), kOutputType.data());

        SaveCanvas(Systematic_TopolCuts.GetQAPlot("AbsDiffRatio"),
                   Form("cQA_ToplogicalCuts_%d_%d", iAnti, centbin),
                   Form("%ssys/Topol/", kFiguresFolder.data()), kOutputType.data());
        SaveCanvas(Systematic_PIDCuts.GetQAPlot("AbsDiffRatio"),
                   Form("cQA_PIDCuts_%d_%d", iAnti, centbin),
                   Form("%ssys/Topol/", kFiguresFolder.data()), kOutputType.data());
        if (fout_systematic_details->IsOpen())
        {
            fout_systematic_details->cd();
            hSysTotalSum->Write();
        }
    }
    if (fout_systematic_details->IsOpen())
        fout_systematic_details->Close();
    return hSysSpectra;

    return hBase;
}
TH1 *GetSpectra(int sys, int iAnti, int centbin, TString sysExt, int sysvar)
{
    LogCustom log(TLogLevel::INFO);
    TString SystematicDirectory = Form("%s%s", kFilterListNames.data(), SystematicBinsName[sys].c_str());
    log.Get(TLogLevel::INFO) << "sys: " << sys << ", caption: " << SystematicBinsName[sys] << ", iAnti: " << iAnti << ", centbin: " << centbin << ", sysExt: " << sysExt << ", sysvar: " << sysvar << ", SystematicDirectory: " << SystematicDirectory.Data();
    bool RunFitSys = false;
    if ((sysvar >= 0) && (sysExt.Length() > 0))
        RunFitSys = true;
    auto isAntiSum = (kAntiBinSize < 2) ? true : false;
    TString spectraPath = kSpectraOutputPart.data();
    spectraPath += Form("_%d", sys);
    if (!isAntiSum)
        spectraPath += Form("_%d", iAnti);
    if (RunFitSys)
        spectraPath += Form("_%s_%d", sysExt.Data(), sysvar);
    TFile *inputfile = new TFile(Form("%s.root", spectraPath.Data()));
    TString hPath = Form("%s/%s/%i/hCorrectedYields", SystematicDirectory.Data(), kParticleType[iAnti].c_str(), centbin);
    log.Get(TLogLevel::INFO) << "hPath: " << hPath.Data();
    auto temp = (TH1 *)inputfile->Get(hPath.Data());
    gROOT->cd();
    TH1D *hReturn = (TH1D *)temp->Clone();
    inputfile->Close();
    return hReturn;
}
TH1 *ErrorSmoothing(TH1 *hinput,
                    double smoothingCriteron, // 1.5
                    int binstart,             // 2
                    int binend)
{ // 8
    bool vebose = false;
    auto htemp = (TH1 *)hinput->Clone();
    for (int bin = 1; bin < htemp->GetNbinsX(); bin++)
    {
        if (bin < binstart)
            continue; // like 0 - 0.8 GeV/c bin
        if (binend < bin)
            continue;                                       // like 8.8 - 15.0 GeV/c bin
        double Error1 = abs(htemp->GetBinContent(bin + 1)); // current bin
        double Error0 = abs(htemp->GetBinContent(bin));     // post bin
        double Error2 = abs(htemp->GetBinContent(bin + 2)); // next bin
        double smoothedValue = abs(htemp->GetBinContent(bin + 1));

        double averageOfNeighbor = (Error1 + Error0 + Error2) / 3;

        if (vebose)
            std::cout << "[Smoothing] bin: " << bin << ", error: " << Error1 << ", average of neighbor: " << averageOfNeighbor;

        if (Error1 > smoothingCriteron * averageOfNeighbor)
        {
            // only if the value is 'significantly' larger than neighbor.
            smoothedValue = averageOfNeighbor;
            if (vebose)
                std::cout << ", smoothed! -> " << smoothedValue;
        }
        if (vebose)
            std::cout << std::endl;
        htemp->SetBinContent(bin + 1, smoothedValue);
    }

    return htemp;
}