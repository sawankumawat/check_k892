#include <AliPWGFunc.h>
#include "src/common.h"
#include "src/logger.h"
#include "src/DrawingHelper.C"
void Normalise(TH1 *h);
void Spectra(int chosenDir = 0, int chosenAnti = -1, int chosenCent = -1, string sys = "", int sysvar = 0)
{
    gErrorIgnoreLevel = kError; // Suppressing warning outputs
    LogCustom log(TLogLevel::INFO);
    bool RunSpecificFolder(false), RunSpecificAnti(false), RunSpecificCent(false), RunFitSys(false);
    if (chosenDir > 0)
        RunSpecificFolder = true;
    if (chosenAnti >= 0)
        RunSpecificAnti = true;
    if (chosenCent >= 0)
        RunSpecificCent = true;
    if (sys.length() > 0)
        RunFitSys = true;
    TString sysOption = sys;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    // Output file naming and figure folder
    TString outputFileName = kSpectraOutputPart.data();
    TString inputSignalFileName = kSignalOutputPart.data();
    TString inputEfficiencyFileName = kEfficiencyOutputPart.data();
    TString figurePath = kFiguresFolder.data();
    if ((chosenDir == 0) && (chosenAnti < 0) && (chosenCent < 0) && !RunFitSys)
    {
        outputFileName = kSpectraOutput.data();
        inputSignalFileName = kSignalOutput.data();
        inputEfficiencyFileName = kEfficiencyOutput.data();
    }
    else
    {
        outputFileName += Form("_%d", chosenDir);
        inputSignalFileName += Form("_%d", chosenDir);
        inputEfficiencyFileName += Form("_%d", chosenDir);
        figurePath += Form("/sys%d", chosenDir);
        if (chosenAnti != -1)
        {
            outputFileName += Form("_%d", chosenAnti);
            inputSignalFileName += Form("_%d", chosenAnti);
            inputEfficiencyFileName += Form("_%d", chosenAnti);
            figurePath += Form("_part%d", chosenAnti);
        }
        if (chosenCent != -1)
        {
            outputFileName += Form("_%d", chosenCent);
            inputSignalFileName += Form("_%d", chosenCent);
            // inputEfficiencyFileName += Form("_%d", chosenCent);
        }
        if (RunFitSys)
        {
            outputFileName += Form("_%s_%d", sys.c_str(), sysvar);
            inputSignalFileName += Form("_%s_%d", sys.c_str(), sysvar);
            figurePath += Form("_%s_%d", sys.c_str(), sysvar);
        }
        inputSignalFileName += ".root";
        outputFileName += ".root";
        inputEfficiencyFileName += ".root";
    }
    figurePath += "/";
    TFile *outputFile = (TFile *)new TFile(outputFileName, "recreate");
    TFile InitFile(kInitOutput.data());
    TFile SignalFile(inputSignalFileName);
    TFile EfficiencyFile(inputEfficiencyFileName);

    log.Get(TLogLevel::INFO) << "Input signal file: " << inputSignalFileName;
    log.Get(TLogLevel::INFO) << "Input efficiency file: " << inputEfficiencyFileName;
    log.Get(TLogLevel::INFO) << "Output file: " << outputFileName;
    log.Get(TLogLevel::INFO) << "Figure path: " << figurePath;

    TCanvas *cFitCheck = GetCanvas("cFitCheck");
    AliPWGFunc *fm = new AliPWGFunc;
    fm->SetVarType(AliPWGFunc::VarType_t(AliPWGFunc::kdNdpt));
    TF1 *func = fm->GetLevi(masspdg, 0.4, 750, 3);
    func->SetParLimits(1, 0.0001, 2e7);
    // Systematic directory name
    const TString SystematicDirectory = Form("%s%s", kFilterListNames.data(), SystematicBins[chosenDir].c_str());
    const TString SystematicDirectoryName = Form("%s%s", kFilterListNames.data(), SystematicBinsName[chosenDir].c_str());
    for (auto list_key : *SignalFile.GetListOfKeys())
    {
        /// Preliminary operation to read the list and create an output dir
        log.Get(TLogLevel::INFO) << "Checking " << list_key->GetName();
        if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos)
        {
            log.Get(TLogLevel::INFO) << "  Analysis folder not found, stop processing";
            continue;
        }
        log.Get(TLogLevel::INFO) << "  Analysis folder found";
        if (RunSpecificFolder && string(list_key->GetName()) != SystematicDirectory) // For multi-core run
        {
            log.Get(TLogLevel::INFO) << "  Not a target folder, skip";
            continue;
        }
        log.Get(TLogLevel::INFO) << "  Start processing " << list_key->GetName();

        TDirectory *baseDir = outputFile->mkdir(SystematicDirectoryName);
        outputFile->cd(list_key->GetName());

        // Read the number of events from the FT0M distribution
        auto hFT0M = (TH1F *)InitFile.Get(Form("%s/hCentFT0M", kFilterListNames.data()));
        log.Get(TLogLevel::INFO) << "  Total Number of events: " << hFT0M->GetEntries();

        auto Efficiencylist = (TDirectory *)EfficiencyFile.Get(SystematicDirectoryName);
        auto Signallist = (TDirectory *)SignalFile.Get(list_key->GetName());

        for (int iAnti = 0; iAnti < kAntiBinSize; iAnti++)
        {
            if (RunSpecificAnti && (chosenAnti != iAnti))
                continue;
            auto isAntiSum = (kAntiBinSize < 2) ? true : false;
            TDirectory *partDir = baseDir->mkdir(kParticleType[iAnti].c_str());
            partDir->cd();
            for (int iCent = 0; iCent < kNMultiplicityBins + 1; iCent++)
            {
                double multiStart{kMultiplicityBins[iCent - 1]}, multiEnd{kMultiplicityBins[iCent]};
                if (RunSpecificCent && (chosenCent != iCent))
                    continue;
                log.Get(TLogLevel::INFO) << "Cent bin: " << iCent;
                if (!kIsINEL && (iCent < 1))
                {
                    multiStart = 0.0;
                    multiEnd = 100.0;
                    log.Get(TLogLevel::INFO) << "MB MODE";
                }
                if (kIsINEL)
                {
                    log.Get(TLogLevel::INFO) << "INEL";
                    multiStart = 0.0;
                    multiEnd = 110.0;
                    if (iCent > 0)
                    {
                        log.Get(TLogLevel::INFO) << "Skip (INEL)";
                        continue; // skip INEL>0
                    }
                }
                TDirectory *centDir = partDir->mkdir(Form("%d", iCent));
                centDir->cd();

                // Number of event in this centrality bin
                auto numberOfEvents = hFT0M->Integral(hFT0M->FindBin(multiStart), hFT0M->FindBin(multiEnd));
                log.Get(TLogLevel::INFO) << "Number of events in this centrality bin: " << numberOfEvents;
                if (kIsINEL)
                {
                    numberOfEvents /= kINELnorm;
                }

                // Read the efficiency
                TString efficiencyName = "hEfficiency";
                if (!isAntiSum)
                    efficiencyName += Form("_%d", iAnti);
                efficiencyName += Form("_cen%d", iCent);
                log.Get(TLogLevel::INFO) << "Efficiency name: " << efficiencyName;
                auto hEfficiency = (TH1F *)Efficiencylist->Get(efficiencyName);
                // auto hEfficiency = (TH1F *)Efficiencylist->Get("hEfficiency"); // for the moment, we use MB efficiency
                if (!hEfficiency)
                {
                    log.Get(TLogLevel::INFO) << "Efficiency not found, stop processing";
                    return;
                }

                // Read the raw yield
                TString rawYieldName = Form("%s/%i/hRawYields", kParticleType[iAnti].data(), iCent);
                auto hRawYields = (TH1F *)Signallist->Get(rawYieldName);
                if (!hRawYields)
                {
                    log.Get(TLogLevel::INFO) << "Raw yield not found!! stop processing";
                    return;
                }

                auto hCorrectedYields = (TH1F *)hRawYields->Clone("hCorrectedYields");
                hCorrectedYields->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                hCorrectedYields->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
                // Normalise(hCorrectedYields);
                for (int i = 1; i <= hCorrectedYields->GetNbinsX(); ++i)
                {
                    if (i < 3) {
                        log.Get(TLogLevel::INFO) << "pT bin: " << i << ", " << hCorrectedYields->GetBinLowEdge(i)
                                                  << " - " << hCorrectedYields->GetBinLowEdge(i + 1);
                        log.Get(TLogLevel::INFO) << "Number of events: " << numberOfEvents << " Raw yields: " << hCorrectedYields->GetBinContent(i) << " efficiency: " << hEfficiency->GetBinContent(i) << " bin widht: " << hCorrectedYields->GetBinWidth(i) << " branching ratio: " << branchingratio;
                        log.Get(TLogLevel::INFO) << "-> corrected value: " << hCorrectedYields->GetBinContent(i) / (numberOfEvents*hEfficiency->GetBinContent(i)* hCorrectedYields->GetBinWidth(i)*2*branchingratio);
                    }   
                    hCorrectedYields->SetBinContent(i, hCorrectedYields->GetBinContent(i) / (hEfficiency->GetBinContent(i)* hCorrectedYields->GetBinWidth(i)));
                    auto originalErr = hCorrectedYields->GetBinError(i) / hCorrectedYields->GetBinWidth(i);
                    auto effiErr = hEfficiency->GetBinError(i) / hEfficiency->GetBinContent(i);
                    auto errSum = sqrt(originalErr * originalErr + effiErr * effiErr);
                    hCorrectedYields->SetBinError(i, errSum);
                }
                // hCorrectedYields->Divide(hEfficiency);
                hCorrectedYields->Scale(1.0 / numberOfEvents);
                hCorrectedYields->Scale(1.0 / branchingratio);
                // if (isAntiSum)
                //     hCorrectedYields->Scale(0.5);

                // hCorrectedYields->Fit(func, "RISqM");
                // auto fitYields = func->Integral(0.0, 10.0);
                // log.Get(TLogLevel::INFO) << "Fit yields: " << fitYields;

                cFitCheck->SetLogy(0);
                cFitCheck->cd();
                hCorrectedYields->Draw();
                TString canvCorSpectraName = "hCorrectedYields";
                if (!isAntiSum)
                    canvCorSpectraName += Form("_%s", kParticleType[iAnti].data());
                canvCorSpectraName += Form("_cen%i", iCent);
                SaveCanvas(cFitCheck, canvCorSpectraName, figurePath, kOutputType.data());
                cFitCheck->SetLogy(1);
                hCorrectedYields->Draw();
                canvCorSpectraName += "_log";
                SaveCanvas(cFitCheck, canvCorSpectraName, figurePath, kOutputType.data());
                hCorrectedYields->Write();

            } // end of loop over centrality
        }     // end of loop over particle type
    }         // end of loop over the list
}
void Normalise(TH1 *h)
{
    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
        h->SetBinContent(i, h->GetBinContent(i) / h->GetBinWidth(i));
        h->SetBinError(i, h->GetBinError(i) / h->GetBinWidth(i));
    }
}