#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

TH1F *RebinFirst10Bins(TH1F *h, const char *newName = "hRebinned")
{
    if (!h)
    {
        std::cerr << "Error: input histogram is null!" << std::endl;
        return nullptr;
    }

    const int nBins = h->GetNbinsX();

    if (nBins < 10)
    {
        std::cerr << "Error: histogram has fewer than 10 bins!" << std::endl;
        return nullptr;
    }

    // First 10 bins are merged pairwise
    const int nNewBins = 5 + (nBins - 10);

    std::vector<double> newEdges;
    newEdges.reserve(nNewBins + 1);

    // ------------------------------------------------------------------
    // First 10 bins -> 5 bins
    // ------------------------------------------------------------------

    newEdges.push_back(h->GetXaxis()->GetBinLowEdge(1));

    newEdges.push_back(h->GetXaxis()->GetBinUpEdge(2));
    newEdges.push_back(h->GetXaxis()->GetBinUpEdge(4));
    newEdges.push_back(h->GetXaxis()->GetBinUpEdge(6));
    newEdges.push_back(h->GetXaxis()->GetBinUpEdge(8));
    newEdges.push_back(h->GetXaxis()->GetBinUpEdge(10));

    // ------------------------------------------------------------------
    // Keep all remaining bins unchanged
    // ------------------------------------------------------------------

    for (int i = 11; i <= nBins; ++i)
    {
        newEdges.push_back(h->GetXaxis()->GetBinUpEdge(i));
    }

    // ------------------------------------------------------------------
    // Create new histogram
    // ------------------------------------------------------------------

    TH1F *hNew = new TH1F(
        newName,
        h->GetTitle(),
        nNewBins,
        newEdges.data());

    hNew->Sumw2();

    // Copy axis titles
    hNew->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    hNew->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());

    // ------------------------------------------------------------------
    // Merge bins 1+2, 3+4, ..., 9+10
    // ------------------------------------------------------------------

    int newBin = 1;

    for (int oldBin = 1; oldBin <= 10; oldBin += 2)
    {

        double content =
            (h->GetBinContent(oldBin) +
             h->GetBinContent(oldBin + 1)) /
            2.0;

        double error = std::sqrt(
                           std::pow(h->GetBinError(oldBin), 2) +
                           std::pow(h->GetBinError(oldBin + 1), 2)) /
                       2.0;

        hNew->SetBinContent(newBin, content);
        hNew->SetBinError(newBin, error);

        ++newBin;
    }

    // ------------------------------------------------------------------
    // Copy bins 11 onwards unchanged
    // ------------------------------------------------------------------

    for (int oldBin = 11; oldBin <= nBins; ++oldBin)
    {

        hNew->SetBinContent(
            newBin,
            h->GetBinContent(oldBin));

        hNew->SetBinError(
            newBin,
            h->GetBinError(oldBin));

        ++newBin;
    }

    // Copy underflow
    hNew->SetBinContent(0, h->GetBinContent(0));
    hNew->SetBinError(0, h->GetBinError(0));

    // Copy overflow
    hNew->SetBinContent(nNewBins + 1, h->GetBinContent(nBins + 1));

    hNew->SetBinError(nNewBins + 1, h->GetBinError(nBins + 1));

    return hNew;
}

void compare_rawCorrecYield()
{
    bool isCorrectedYield = true;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    bool isINEL = false;

    string defaultName = "BW + pol3";
    string variationName = "Voigt + Template";

    // string path1 = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED"; // Default1
    // string path2 = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/751768/kstarqa/hInvMass/WidthFree";   // Variation

    string path1 = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/ROTATED"; // Default
    string path2 = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED"; // Variation

    TString outputPath = path2;

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    TH1F *hmult1[numofmultbins + 1];
    TH1F *hmult2[numofmultbins + 1];

    TH1F *hefficiency1[numofmultbins + 1];
    TH1F *hefficiency2[numofmultbins + 1];

    // for (int imult = 0; imult < numofmultbins + 1; imult++)
    for (int imult = 0; imult < 1; imult++)
    {
        double multlow = (imult == 0) ? 0 : mult_classes[imult - 1];
        double multhigh = (imult == 0) ? (isINEL) ? 120 : 100 : mult_classes[imult];

        TFile *fspectra1 = (isCorrectedYield) ? new TFile((path1 + Form("/corrected_spectra_%d_%d.root", (int)multlow, (int)multhigh)).c_str(), "read") : new TFile((path1 + Form("/yield_%d_%d.root", (int)multlow, (int)multhigh)).c_str(), "read");
        TFile *fspectra2 = (isCorrectedYield) ? new TFile((path2 + Form("/corrected_spectra_%d_%d.root", (int)multlow, (int)multhigh)).c_str(), "read") : new TFile((path2 + Form("/yield_%d_%d.root", (int)multlow, (int)multhigh)).c_str(), "read");

        if (fspectra1->IsZombie() || fspectra2->IsZombie())
        {
            cout << "Error: files not found" << endl;
            return;
        }

        hmult1[imult] = (isCorrectedYield) ? (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral_final", multlow, multhigh)) : (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));
        // hmult2[imult] = (isCorrectedYield) ? (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral_final", multlow, multhigh)) : (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));

        TH1F *hmult2_temp = (isCorrectedYield) ? (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral_final", multlow, multhigh)) : (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));

        // Lets rebin the variation histogram
        hmult2[imult] = RebinFirst10Bins(hmult2_temp, Form("hmult2_rebinned_mult_%.0f-%.0f", multlow, multhigh));

        if (hmult1[imult] == nullptr)
        {
            cout << "Histogram others not found" << endl;
            return;
        }

        TH1F *hratio1 = (TH1F *)hmult2[imult]->Clone(Form("ratio_mult_%.0f-%.0f", multlow, multhigh));
        hratio1->Divide(hmult1[imult]);

        TCanvas *c1 = new TCanvas(Form("c1_mult_%.0f-%.0f", multlow, multhigh), "c1", 720, 720);
        SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
        double pad1Size, pad2Size;
        canvas_style(c1, pad1Size, pad2Size);
        c1->cd(1);
        gPad->SetLogy(1);
        SetHistoStyle(hmult1[imult], 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
        hmult1[imult]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        hmult1[imult]->SetMaximum(hmult1[imult]->GetMaximum() * 10);
        hmult1[imult]->SetMinimum(hmult1[imult]->GetMinimum() * 0.5);
        hmult1[imult]->GetYaxis()->SetTitleOffset(1.30);
        hmult1[imult]->GetXaxis()->SetTitleOffset(1.02);
        hmult1[imult]->SetMarkerStyle(20);
        hmult1[imult]->SetMarkerSize(1);
        hmult1[imult]->Draw("pe");
        hmult2[imult]->SetMarkerStyle(21);
        hmult2[imult]->SetMarkerSize(1);
        hmult2[imult]->SetMarkerColor(kBlue);
        hmult2[imult]->SetLineColor(kBlue);
        hmult2[imult]->SetLineWidth(2);
        hmult2[imult]->Draw("pe same");

        TLegend *leg = new TLegend(0.65, 0.75, 0.89, 0.93);
        SetLegendStyle(leg);
        leg->SetTextSize(0.04);
        leg->SetHeader(Form("FT0M: %.0f-%.0f%%", multlow, multhigh));
        leg->AddEntry(hmult1[imult], Form("%s", defaultName.c_str()), "p");
        leg->AddEntry(hmult2[imult], Form("%s", variationName.c_str()), "p");
        leg->Draw();

        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextFont(42);
        latex->SetTextSize(0.04);
        latex->DrawLatex(0.35, 0.90, "ALICE");
        latex->DrawLatex(0.35, 0.84, "pp #sqrt{#it{s}} = 13.6 TeV");
        latex->DrawLatex(0.35, 0.78, "|y| < 0.5");

        c1->cd(2);
        TH1F *hdummy = (TH1F *)hmult1[0]->Clone();
        for (int i = 0; i < hdummy->GetNbinsX(); i++)
        {
            hdummy->SetBinContent(i + 1, 0);
            hdummy->SetBinError(i + 1, 0);
        }

        SetHistoQA(hratio1);
        gPad->SetGrid(1, 1);
        hratio1->GetYaxis()->SetTitleSize(0.035 / pad2Size);
        hratio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        hratio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        hratio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        hratio1->SetMarkerStyle(20);
        hratio1->SetMarkerSize(1.0);
        hratio1->SetMarkerColor(kBlue);
        hratio1->SetLineColor(kBlue);
        // hratio1->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
        // hratio1->GetYaxis()->SetTitle("#frac{2023 data}{2022 data}");
        hratio1->GetYaxis()->SetTitle("Var/Def");
        hratio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hratio1->GetXaxis()->CenterTitle(1);
        hratio1->GetYaxis()->SetTitleOffset(0.6);
        hratio1->GetXaxis()->SetTitleOffset(1.1);
        hratio1->GetYaxis()->SetNdivisions(510);
        // hratio1->SetMaximum(hratio1->GetMaximum() * 1.3);
        // hratio1->SetMinimum(hratio1->GetMinimum() * 0.7);
        // hratio1->GetXaxis()->SetRangeUser(0, 20);
        // hratio1->GetYaxis()->SetRangeUser(0.98, 1.11);
        hratio1->GetYaxis()->SetRangeUser(0.34, 1.69);
        hratio1->GetYaxis()->SetNdivisions(510);
        hratio1->Draw("p");

        TLine *line = (isINEL) ? new TLine(0, 1, 30, 1) : new TLine(0, 1, 20, 1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(1);
        line->Draw();
        if (isCorrectedYield)
        {
            c1->SaveAs(outputPath + Form("/CorrectedYield_%.0f-%.0f.png", multlow, multhigh));
        }
        else
        {
            c1->SaveAs(outputPath + Form("/RawYield_%.0f-%.0f.png", multlow, multhigh));
        }
    }
}

// string path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/459845/kstarqa/hInvMass"; // 2022 data
// string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/466180/kstarqa_id33593/hInvMass"; // 2024 data
// string path3 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/459908/kstarqa/hInvMass";         // 2023 (135 kHz) data
// string path3 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/LHC23z/kstarqa/hInvMass";         // 2023 (450 kHz) data
// string path4 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/LHC23ls/kstarqa/hInvMass"; // 2023 (650 kHz) data

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}