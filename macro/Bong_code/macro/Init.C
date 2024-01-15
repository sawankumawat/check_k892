#include "src/common.h"
#include "src/logger.h"
#include "src/DrawingHelper.C"
void Init(int chosenDir = 0)
{
    gErrorIgnoreLevel = kError; // Suppressing warning outputs
    LogCustom log(TLogLevel::INFO);
    // Make a folder structure
    TString path = Form("%s/common/", kBaseOutputDir.c_str());
    MakeFolder(path);
    path = Form("%s/systematic_details/", kBaseOutputDir.c_str());
    MakeFolder(path);
    path = Form("%s/images/", kBaseOutputDir.c_str());
    MakeFolder(path);
    path = Form("%s/figures/", kBaseOutputDir.c_str());
    MakeFolder(path);
        
    TFile inputDataFile(kDataFilename.data());
    TFile inputMCFile(kMCfilename.data());
    TString outputFileName = kInitOutputPart.data();
    if (chosenDir != 0)
        outputFileName = outputFileName + Form("_%d", chosenDir);
    outputFileName += ".root";
    TFile *outputFile = (TFile *)new TFile(outputFileName, "recreate");
    TString figurePath = kFiguresFolder.data();
    gStyle->SetOptStat(1);
    gStyle->SetOptTitle(1);
    gStyle->SetLegendBorderSize(0);
    // Output file naming and figure folder

    log.Get(TLogLevel::INFO) << "Input data file: " << kDataFilename;
    log.Get(TLogLevel::INFO) << "Output file: " << outputFileName;
    log.Get(TLogLevel::INFO) << "Figure folder: " << figurePath;

    TCanvas *cEventQA = GetCanvas("cEventQA");

    int nDir = 0;
    for (auto list_key : *inputDataFile.GetListOfKeys())
    {
        /// Preliminary operation to read the list and create an output dir
        log.Get(TLogLevel::INFO) << "Checking " << list_key->GetName();
        if (string(list_key->GetName()).find(kFilterListResoTaskNames.data()) == string::npos)
        {
            continue;
        }
        log.Get(TLogLevel::INFO) << "  Analysis folder found";
        nDir++;
        if ((chosenDir > 0) && (nDir != chosenDir))
            continue;
        log.Get(TLogLevel::INFO) << "  Start processing " << list_key->GetName();

        TDirectory *baseDir = outputFile->mkdir(kFilterListNames.data());
        baseDir->cd();

        auto hPosZ = (TH1F *)inputDataFile.Get(Form("%s/Event/posZ", list_key->GetName()));
        auto hCentFT0M = (TH1F *)inputDataFile.Get(Form("%s/Event/CentFT0M", list_key->GetName()));

        // Basic QA Plots
        // ------------------------------------------------------------------------------------------
        // Vertex Z distribution
        cEventQA->cd();
        hPosZ->Draw();
        SaveCanvas(cEventQA, "hVertexPosZ", figurePath, kOutputType.data());
        hPosZ->Write("hVertexPosZ");

        // Multiplicity distribution
        cEventQA->cd();
        hCentFT0M->Draw();
        SaveCanvas(cEventQA, "hCentFT0M", figurePath, kOutputType.data());
        hCentFT0M->Write("hCentFT0M");

        if (kDoMCQA)
        {
            // MC Plots
            auto hPosZ_MC = (TH1F *)inputMCFile.Get(Form("%s/Event/posZ", list_key->GetName()));
            auto hCentFT0M_MC = (TH1F *)inputMCFile.Get(Form("%s/Event/CentFT0M", list_key->GetName()));

            // Vertex Z distribution
            cEventQA->cd();
            hPosZ_MC->Draw();
            SaveCanvas(cEventQA, "hVertexPosZ_MC", figurePath, kOutputType.data());
            hPosZ_MC->Write("hVertexPosZ_MC");

            // Multiplicity distribution
            cEventQA->cd();
            hCentFT0M_MC->Draw();
            SaveCanvas(cEventQA, "hCentFT0M_MC", figurePath, kOutputType.data());
            hCentFT0M_MC->Write("hCentFT0M_MC");
        }
    }
    outputFile->Close();
}