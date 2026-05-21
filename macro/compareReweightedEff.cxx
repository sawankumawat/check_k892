#include "src/style.h"

int FindHighestIndex(TFile *file, const string &baseHistoName)
{
    int maxIndex = -1;
    for (int i = 10; i >= 0; i--) // Check from i10 down to i0
    {
        string histoName = baseHistoName + "i" + to_string(i);
        TObject *obj = file->Get(histoName.c_str());
        if (obj != nullptr)
        {
            maxIndex = i;
            cout << "Highest index found: " << maxIndex << endl;
            break; // Found the highest index
        }
    }
    return maxIndex;
}

void compareReweightedEff()
{
    string path = "../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass"; // 2024 data
    TString reweightEffPath = path + "/ReweightEfficiency";

    TFile *fspectra = new TFile((path + "/corrected_spectra.root").c_str(), "read");
    if (fspectra->IsZombie())
    {
        cout << "File not found" << endl;
        return;
    }
    TH1F *hMultiplicity = (TH1F *)fspectra->Get("Multiplicity");
    int entries = hMultiplicity->Integral(hMultiplicity->FindBin(0 + 0.001), hMultiplicity->FindBin(100 - 0.001));
    entries = entries * hMultiplicity->GetBinWidth(1) * 1.3;

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    int multLow, multHigh;

    // for (int i = 0; i < numofmultbins + 1; i++)
    for (int i = 0; i < 1; i++)
    {
        if (i == 0)
        {
            multLow = 0;
            multHigh = 100;
        }
        else
        {
            multLow = (int)mult_classes[i - 1];
            multHigh = (int)mult_classes[i];
        }
        TFile *fReweightedSpectra = new TFile(Form("%s/ReweightFactor_mult_%d-%d.root", reweightEffPath.Data(), multLow, multHigh), "read");
        if (fReweightedSpectra->IsZombie())
        {
            cout << "File not found: " << Form("%s/ReweightFactor_mult_%d-%d.root", reweightEffPath.Data(), multLow, multHigh) << endl;
            return;
        }

        int maxIndexReweighted = FindHighestIndex(fReweightedSpectra, "hlossCorrected_integral_");
        string indexStr = "i" + to_string(maxIndexReweighted);

        TH1F *hReweightedSpectra = (TH1F *)fReweightedSpectra->Get(Form("hlossCorrected_integral_%s", indexStr.c_str()));
        TH1F *hGenReweighted = (TH1F *)fReweightedSpectra->Get(Form("hk892GenpTCalib1_proj_0_%s", indexStr.c_str()));
        TH1F *hRecReweighted = (TH1F *)fReweightedSpectra->Get(Form("h2KstarRecptCalib1_proj_0_%s", indexStr.c_str()));

        TH1F *hOriginalSpectra = (TH1F *)fReweightedSpectra->Get(Form("hlossCorrected_integral_i%d", 0));
        TH1F *hGen = (TH1F *)fReweightedSpectra->Get(Form("hk892GenpTCalib1_proj_0_i%d", 0));
        TH1F *hRec = (TH1F *)fReweightedSpectra->Get(Form("h2KstarRecptCalib1_proj_0_i%d", 0));
        TH1F *hCorrFactor = (TH1F *)fReweightedSpectra->Get(Form("hlossCorrected_integral_correction_%s", indexStr.c_str()));

        if (hReweightedSpectra == nullptr || hGenReweighted == nullptr || hRecReweighted == nullptr || hOriginalSpectra == nullptr || hGen == nullptr || hRec == nullptr)
        {
            cout << "Histogram not found in reweighting file: " << endl;
            return;
        }
        hGen->Scale(1.0 / entries);
        hRec->Scale(1.0 / entries);

        TCanvas *cCompareBefore = new TCanvas("", "", 720, 720);
        SetCanvasStyle(cCompareBefore, 0.18, 0.06, 0.01, 0.14);
        gPad->SetLogy();
        SetHistoQA(hOriginalSpectra);
        hOriginalSpectra->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hOriginalSpectra->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
        hOriginalSpectra->GetYaxis()->SetTitleOffset(1.7);
        hOriginalSpectra->GetXaxis()->SetRangeUser(0.0, 5.0);
        hOriginalSpectra->SetMinimum(1e-4);
        hOriginalSpectra->SetMaximum(hOriginalSpectra->GetMaximum() * 5);
        hOriginalSpectra->SetMarkerSize(1.5);
        hOriginalSpectra->Draw("ep");
        SetHistoQA(hGen);
        hGen->SetMarkerStyle(24);
        hGen->SetMarkerColor(kGreen + 2);
        hGen->SetLineColor(kGreen + 2);
        hGen->Draw("ep same");
        hRec->SetMarkerStyle(24);
        hRec->SetMarkerColor(kRed + 2);
        hRec->SetLineColor(kRed + 2);
        hRec->Draw("ep same");

        TCanvas *cCompareAfter = new TCanvas("", "", 720, 720);
        SetCanvasStyle(cCompareAfter, 0.18, 0.06, 0.01, 0.14);
        gPad->SetLogy();
        SetHistoQA(hReweightedSpectra);
        hReweightedSpectra->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hReweightedSpectra->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
        hReweightedSpectra->GetYaxis()->SetTitleOffset(1.7);
        hReweightedSpectra->GetXaxis()->SetRangeUser(0.0, 5.0);
        hReweightedSpectra->SetMinimum(1e-4);
        hReweightedSpectra->SetMaximum(hReweightedSpectra->GetMaximum() * 5);
        hReweightedSpectra->SetMarkerSize(1.5);
        hReweightedSpectra->Draw("ep");
        SetHistoQA(hGenReweighted);
        hGenReweighted->SetLineColor(kGreen + 2);
        hGenReweighted->SetMarkerStyle(24);
        hGenReweighted->SetMarkerColor(kGreen + 2);
        hGenReweighted->Draw("p same");
        SetHistoQA(hRecReweighted);
        hRecReweighted->SetLineColor(kRed);
        hRecReweighted->SetMarkerStyle(24);
        hRecReweighted->SetMarkerColor(kRed);
        hRecReweighted->Draw("p same");

        TLegend *legAfter = new TLegend(0.5, 0.7, 0.9, 0.9);
        legAfter->SetBorderSize(0);
        legAfter->SetTextFont(42);
        legAfter->SetTextSize(0.036);
        legAfter->SetHeader(Form("pp #sqrt{#it{s}} = 13.6 TeV, %d-%d%%", multLow, multHigh));
        legAfter->AddEntry(hReweightedSpectra, "Data", "ep");
        legAfter->AddEntry(hGenReweighted, "Generated MC", "p");
        legAfter->AddEntry(hRecReweighted, "Reconstructed MC", "p");
        legAfter->Draw();
        cCompareAfter->SaveAs(reweightEffPath + Form("/CompareReweight_mult_%d-%d.pdf", multLow, multHigh));
        cCompareBefore->SaveAs(reweightEffPath + Form("/CompareWoReweight_mult_%d-%d.pdf", multLow, multHigh));

        TCanvas *cCorrFactor = new TCanvas("", "", 720, 720);
        SetCanvasStyle(cCorrFactor, 0.18, 0.06, 0.01, 0.14);
        SetHistoQA(hCorrFactor);
        hCorrFactor->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hCorrFactor->GetYaxis()->SetTitle("Re-weighting factor");
        hCorrFactor->GetYaxis()->SetTitleOffset(1.7);
        hCorrFactor->GetXaxis()->SetRangeUser(0.0, 5.0);
        hCorrFactor->SetMinimum(0.91);
        hCorrFactor->SetMaximum(1.09);
        hCorrFactor->SetMarkerSize(1.5);
        hCorrFactor->Draw("HIST");
        TLegend *legCorr = new TLegend(0.5, 0.7, 0.9, 0.9);
        legCorr->SetBorderSize(0);
        legCorr->SetTextFont(42);
        legCorr->SetTextSize(0.036);
        legCorr->SetHeader(Form("pp #sqrt{#it{s}} = 13.6 TeV, %d-%d%%", multLow, multHigh));
        // legCorr->AddEntry(hCorrFactor, "Re-weighting factor", "l");
        legCorr->Draw();
        cCorrFactor->SaveAs(reweightEffPath + Form("/ReweightingFactor_mult_%d-%d.pdf", multLow, multHigh));
    }

    // TCanvas *cCorrectionFactor = new TCanvas("", "", 720, 720);
    // SetCanvasStyle(cCorrectionFactor, 0.18, 0.06, 0.01, 0.14);

    // for (int i = 0; i < numofmultbins + 1; i++)
    // {
    //     if (i == 0)
    //     {
    //         multLow = 0;
    //         multHigh = 100;
    //     }
    //     else
    //     {
    //         multLow = (int)mult_classes[i - 1];
    //         multHigh = (int)mult_classes[i];
    //     }
    //     TFile *fReweightedSpectra = new TFile(Form("%s/ReweightFactor_mult_%d-%d.root", reweightEffPath.Data(), multLow, multHigh), "read");
    //     if (fReweightedSpectra->IsZombie())
    //     {
    //         cout << "File not found: " << Form("%s/ReweightFactor_mult_%d-%d.root", reweightEffPath.Data(), multLow, multHigh) << endl;
    //         return;
    //     }

    //     int maxIndexReweighted = FindHighestIndex(fReweightedSpectra, "hlossCorrected_integral_");
    //     string indexStr = "i" + to_string(maxIndexReweighted);

    //     TH1F *hCorrFactor = (TH1F *)fReweightedSpectra->Get(Form("hlossCorrected_integral_correction_%s", indexStr.c_str()));

    //     SetHistoQA(hCorrFactor);
    //     hCorrFactor->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    //     hCorrFactor->GetYaxis()->SetTitle("Re-weighting factor");
    //     hCorrFactor->GetYaxis()->SetTitleOffset(1.7);
    //     hCorrFactor->GetXaxis()->SetRangeUser(0.0, 5.0);
    //     hCorrFactor->SetMinimum(0.91);
    //     hCorrFactor->SetMaximum(1.09);
    //     hCorrFactor->SetMarkerSize(1.5);
    //     hCorrFactor->Draw("HIST same PLC PMC");
    // }
}