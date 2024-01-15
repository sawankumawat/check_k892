#include "src/common.h"
#include "src/logger.h"
#include "src/DrawingHelper.C"

std::vector<TString> runNumbers = {
    "519906",
    "519905",
    "519904",
    "519903",
    "519507",
    "519506",
    "519504",
    "519503",
    "519502",
    "519499",
    "519498",
    "519497",
    "519045",
    "519043",
    "519041",
    "518547",
    "518546",
    "518543",
    "518541",
    "517753",
    "517751",
    "517748",
    "517737",
    "517690",
    "517685",
    "517679",
    "517678",
    "517677",
    "517623",
    "517620",
    "517619"
};

void EffiCheck()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gErrorIgnoreLevel = kError; // Suppressing warning outputs
    LogCustom log(TLogLevel::INFO);
    log.Get(TLogLevel::INFO) << "Efficiecny check";
    // Output file under the kBaseOutputDir
    TString figurePath = kFiguresFolder.data();
    TFile outputFile(Form("%s/common/efficheck.root", kBaseOutputDir.c_str()), "recreate");

    // loop over runNumbers and print out the run number
    TCanvas *cEffiCheck = GetCanvas("cEffiCheck");
    TCanvas *cEffiCheckAll = GetCanvas("cEffiCheckAll");
    TH1F *hGen_temp = nullptr;
    TH1F *hRec_temp = nullptr;
    int sysColorPallet = GetSerialColors(runNumbers.size());
    int runBin = 0;
    TLegend *leg = new TLegend(0.15, 0.5, 0.5, 0.9);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
    for (auto runNumber : runNumbers)
    {
        log.Get(TLogLevel::INFO) << "Input Run number: " << runNumber;
        // input file path: kBaseInputDir + "mc/LHC22h1c1/" + runNumber + ".root"
        auto inputFile = TFile::Open(Form("%s/mc/LHC22h1d1/%s.root",
                                          kBaseInputDir.c_str(), runNumber.Data()));
        log.Get(TLogLevel::INFO) << "Input file: " << inputFile->GetName();

        outputFile.cd();
        log.Get(TLogLevel::INFO) << Form("%s/k892Gen", kOutputName.c_str());
        auto th1Gen =
            (TH3F *)inputFile->Get(Form("%s/k892Gen", kOutputName.c_str()));
        auto th1Gen_Anti =
            (TH3F *)inputFile->Get(Form("%s/k892GenAnti", kOutputName.c_str()));
        auto th1Rec =
            (TH3F *)inputFile->Get(Form("%s/k892Rec", kOutputName.c_str()));
        auto th1Rec_Anti =
            (TH3F *)inputFile->Get(Form("%s/k892RecAnti", kOutputName.c_str()));

        std::vector<TH3F *> hGen = {th1Gen, th1Gen_Anti};
        std::vector<TH3F *> hRec = {th1Rec, th1Rec_Anti};

        std::vector<Double_t> receffi = {};
        std::vector<Double_t> receffi_err = {};
        auto iCent = 0;
        for (int iAnti = 0; iAnti < kAntiBinSize; iAnti++)
        {
            hGen_temp = (TH1F *)hGen[iAnti]->Clone(Form("hGen_%d_%d", iAnti, iCent));
            hRec_temp = (TH1F *)hRec[iAnti]->Clone(Form("hRec_%d_%d", iAnti, iCent));
            auto isAntiSum = (kAntiBinSize < 2) ? true : false;
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
                log.Get(TLogLevel::INFO) << "pT bin: " << iPt << ", " << kpTbin[iPt]
                                         << " - " << kpTbin[iPt + 1];
                auto Gen = hGen_temp->Integral(
                    hGen_temp->GetXaxis()->FindBin(kpTbin[iPt] + 0.001),
                    hGen_temp->GetXaxis()->FindBin(kpTbin[iPt + 1] - 0.001));
                auto Rec = hRec_temp->Integral(
                    hRec_temp->GetXaxis()->FindBin(kpTbin[iPt] + 0.001),
                    hRec_temp->GetXaxis()->FindBin(kpTbin[iPt + 1] - 0.001));

                log.Get(TLogLevel::INFO)
                    << "Gen: " << Gen << ", Rec: " << Rec << " Ratio: " << Rec / Gen;
                receffi.push_back(Rec / Gen);
                receffi_err.push_back(sqrt(receffi[iPt] * (1 - receffi[iPt]) / Gen));
            } // pT loop
            // Efficiency
            auto hEfficiency = MakeHistfromArray("Effi", receffi, kpTbin, receffi_err);
            hEfficiency->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hEfficiency->GetYaxis()->SetTitle("Acceptance x Efficiency x BR");
            hEfficiency->SetMaximum(0.6);
            hEfficiency->SetMinimum(1e-5);
            hEfficiency->SetLineColor(sysColorPallet + runNumbers.size() - runBin - 1);
            hEfficiency->SetMarkerColor(sysColorPallet + runNumbers.size() - runBin - 1);
            hEfficiency->SetMarkerStyle(20);
            hEfficiency->SetMarkerSize(0.5);
            TString outputFigName = "hEfficiency";
            outputFigName += Form("_LHC22h1d1_%s", runNumber.Data());
            if (!isAntiSum)
                outputFigName += Form("_%d", iAnti);
            outputFigName += Form("_cen%i", iCent);
            hEfficiency->Write(outputFigName);
            cEffiCheck->cd();

            TLegend *leg_run = new TLegend(0.15, 0.5, 0.5, 0.9);
            leg_run->SetBorderSize(0);
            leg_run->SetFillStyle(0);
            leg_run->AddEntry(hEfficiency, Form("%s", runNumber.Data()), "lpe");
            leg->AddEntry(hEfficiency, Form("%s", runNumber.Data()), "lpe");

            hEfficiency->Draw("PZE");
            leg_run->Draw();
            SaveCanvas(cEffiCheck, outputFigName, figurePath, kOutputType.data());
            cEffiCheckAll->cd();
            if (runBin == 0)
                hEfficiency->Draw("PZE");
            else
                hEfficiency->Draw("PZE same");
            runBin++;
        }
    }
    cEffiCheckAll->cd();
    leg->Draw();
    SaveCanvas(cEffiCheckAll, "hEfficiencyAllRuns_LHC22h1d1", figurePath, kOutputType.data());    
}
