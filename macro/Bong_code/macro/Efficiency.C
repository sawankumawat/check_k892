#include "src/common.h"
#include "src/logger.h"
#include "src/DrawingHelper.C"

void Efficiency(int chosenDir = 0, int chosenAnti = -1, int chosenCent = -1)
{
    gErrorIgnoreLevel = kError; // Suppressing warning outputs
    LogCustom log(TLogLevel::INFO);
    bool RunSpecificFolder(false), RunSpecificAnti(false), RunSpecificCent(false);
    if (chosenDir > 0)
        RunSpecificFolder = true;
    if (chosenAnti >= 0)
        RunSpecificAnti = true;
    if (chosenCent >= 0)
        RunSpecificCent = true;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    TFile inputFile(kMCfilename.data());
    TFile *outputFile;
    TString outputFileName = kEfficiencyOutputPart.data();
    TString figurePath = kFiguresFolder.data();
    if ((chosenDir == 0) && (chosenAnti < 0))
    {
        outputFile = (TFile *)new TFile(kEfficiencyOutput.data(), "recreate");
    }
    else
    {
        outputFileName += Form("_%d", chosenDir);
        figurePath += Form("/sys%d", chosenDir);
        if (chosenAnti != -1)
        {
            outputFileName += Form("_%d", chosenAnti);
            figurePath += Form("_part%d", chosenAnti);
        }
        if (chosenCent != -1)
            outputFileName += Form("_%d", chosenCent);
        outputFileName += ".root";
        outputFile = (TFile *)new TFile(outputFileName, "recreate");
    }
    log.Get(TLogLevel::INFO) << "Output file: " << outputFileName;
    log.Get(TLogLevel::INFO) << "Figure path: " << figurePath;

    TCanvas *cEffiCheck = GetCanvas("cEffiCheck");
    TH1F *hGen_temp = nullptr;
    TH1F *hRec_temp = nullptr;
    TString SystematicDirectory = Form("%s%s", kFilterListNamesMC.data(),
                                       SystematicBinsName[chosenDir].c_str());
    for (auto list_key : *inputFile.GetListOfKeys())
    {
        /// Preliminary operation to read the list and create an output dir
        log.Get(TLogLevel::INFO) << "Checking " << list_key->GetName();
        if (string(list_key->GetName()).find(kFilterListNamesMC.data()) ==
            string::npos)
            continue;
        log.Get(TLogLevel::INFO) << "  Analysis folder found";
        if (RunSpecificFolder && string(list_key->GetName()) != SystematicDirectory){ // For multi-core run
            continue;
        }
        TDirectory *baseDir = outputFile->mkdir(list_key->GetName());
        outputFile->cd(list_key->GetName());
        log.Get(TLogLevel::INFO) << Form("%s/%s/k892Gen", list_key->GetName(), kOutputName.c_str());
        auto th1Gen = (TH3F *)inputFile.Get(Form("%s/k892Gen", list_key->GetName()));
        auto th1Gen_Anti = (TH3F *)inputFile.Get(Form("%s/k892GenAnti", list_key->GetName()));
        auto th1Rec = (TH3F *)inputFile.Get(Form("%s/k892Rec", list_key->GetName()));
        auto th1Rec_Anti = (TH3F *)inputFile.Get(Form("%s/k892RecAnti", list_key->GetName()));

        std::vector<TH3F *> hGen = {th1Gen, th1Gen_Anti};
        std::vector<TH3F *> hRec = {th1Rec, th1Rec_Anti};

        std::vector<Double_t> receffi = {};
        std::vector<Double_t> receffi_err = {};

        for (int iAnti = 0; iAnti < kAntiBinSize; iAnti++)
        {
            if (RunSpecificAnti && (chosenAnti != iAnti))
                continue;
            auto isAntiSum = (kAntiBinSize < 2) ? true : false;

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
                hGen_temp = (TH1F *)hGen[iAnti]->Clone(Form("hGen_%d_%d", iAnti, iCent));
                hRec_temp = (TH1F *)hRec[iAnti]->Clone(Form("hRec_%d_%d", iAnti, iCent));
                if (isAntiSum)
                {
                    if (iAnti > 0)
                        continue; // skip anti particle if requested
                    log.Get(TLogLevel::INFO) << "Normal + Anti";
                    hGen_temp->Add(hGen[abs(1 - iAnti)]); // 0: 1, 1: 0
                    hRec_temp->Add(hRec[abs(1 - iAnti)]);
                }
                receffi.clear();
                receffi_err.clear();
                for (Int_t iPt = 0, nPt = kpTbin.size() - 1; iPt < nPt; ++iPt)
                {
                    log.Get(TLogLevel::INFO) << "pT bin: " << iPt << ", " << kpTbin[iPt] << " - " << kpTbin[iPt + 1];
                    auto Gen = hGen_temp->Integral(hGen_temp->GetXaxis()->FindBin(kpTbin[iPt] + 0.001),
                                                   hGen_temp->GetXaxis()->FindBin(kpTbin[iPt + 1] - 0.001));
                    auto Rec = hRec_temp->Integral(hRec_temp->GetXaxis()->FindBin(kpTbin[iPt] + 0.001),
                                                   hRec_temp->GetXaxis()->FindBin(kpTbin[iPt + 1] - 0.001));

                    log.Get(TLogLevel::INFO) << "Gen: " << Gen << ", Rec: " << Rec << " Ratio: " << Rec / Gen;
                    receffi.push_back(Rec / Gen);
                    receffi_err.push_back(sqrt(receffi[iPt] * (1 - receffi[iPt]) / Gen));

                } // pT loop
                // Efficiency
                auto hEfficiecny = MakeHistfromArray("Effi", receffi, kpTbin, receffi_err);
                hEfficiecny->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                hEfficiecny->GetYaxis()->SetTitle("Acceptance x Efficiency");
                hEfficiecny->SetMaximum(0.6);
                hEfficiecny->SetMinimum(1e-5);
                TString outputFigName = "hEfficiency";
                if (!isAntiSum)
                    outputFigName += Form("_%d", iAnti);
                outputFigName += Form("_cen%i", iCent);
                hEfficiecny->Write(outputFigName);
                cEffiCheck->cd();
                hEfficiecny->Draw();
                SaveCanvas(cEffiCheck, outputFigName, figurePath, kOutputType.data());

                if (kIsINEL) // INEL
                    break;
                if (RunSpecificCent)
                    break;
            } // Anti
            baseDir->Close();
        }
    }
    outputFile->Close();
}