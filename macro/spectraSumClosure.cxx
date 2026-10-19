#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include "src/style.h"

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void spectraSumClosure()
{
    int mult_classes[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    int totalMultBins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;
    TH1D *hSpectra[totalMultBins];
    double eventsperMultBin[totalMultBins];
    double minimumBiasEvents = 0;
    // double eventLoss[] = {0.707802, 0.680640, 0.674949, 0.667825, 0.659018, 0.641315, 0.620520, 0.603644, 0.559658, 0.394608}; // Default
    double eventLoss[] = {0.723197, 0.718152, 0.714390, 0.711143, 0.708311, 0.703969, 0.697424, 0.688330, 0.662817, 0.492254}; // lf-mc-centrality

    // double eventLossMinBias = 0.594455; // Default
    double eventLossMinBias = 0.592687; // lf-mc-centrality

    double eventLossSum = 0.0;
    double rawEventsSum = 0.0;

    for (int i = 0; i < 10; i++)
    {
        eventLossSum += eventLoss[i];
    }

    TFile *fInputFile = new TFile("/home/sawan/Storage/check_k892/data/kstar/LHC22o_pass7/749276.root", "read");
    string multpath = "kstarqa";
    TH1F *hmult = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hMultiplicity", multpath.c_str()));
    minimumBiasEvents = hmult->Integral(hmult->GetXaxis()->FindBin(0.0 + 1e-5), hmult->GetXaxis()->FindBin(100.0 - 1e-5));

    double totalBefore30 = 0.0;

    for (int imult = 0; imult < totalMultBins; imult++)
    {
        int multLow = mult_classes[imult];
        int multHigh = mult_classes[imult + 1];

        eventsperMultBin[imult] = hmult->Integral(hmult->GetXaxis()->FindBin(multLow + 1e-5), hmult->GetXaxis()->FindBin(multHigh - 1e-5));
        // cout << "Total events in multiplicity bin " << multLow << "-" << multHigh << " = " << eventsperMultBin[imult] << ", corrected events = " << eventsperMultBin[imult] * eventLoss[imult] << endl;
        rawEventsSum += eventsperMultBin[imult];

        // if (imult < 6)
        //     totalBefore30 += eventsperMultBin[imult] * eventLoss[imult];

        // if (imult == totalMultBins - 1)
        //     cout << "Total events before 30 = " << totalBefore30 << endl;

        TFile *fSpectra = new TFile(Form("../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/corrected_spectra_%d_%d.root", multLow, multHigh), "read");
        hSpectra[imult] = (TH1D *)fSpectra->Get(Form("mult_%d-%d/corrected_spectra_Integral_final", multLow, multHigh));

        // TFile *fSpectra = new TFile(Form("../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/yield_%d_%d.root", multLow, multHigh), "read");
        // hSpectra[imult] = (TH1D *)fSpectra->Get(Form("mult_%d-%d/yield_integral", multLow, multHigh));
    }

    // compare with minimum bias spectra
    TFile *fMinBias = new TFile("../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/corrected_spectra_0_100.root", "read");
    TH1D *hMinBias = (TH1D *)fMinBias->Get("mult_0-100/corrected_spectra_Integral_final");

    // TFile *fMinBias = new TFile("../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/yield_0_100.root", "read");
    // TH1D *hMinBias = (TH1D *)fMinBias->Get("mult_0-100/yield_integral");

    double pad1Size, pad2Size;

    TCanvas *cSpectra = new TCanvas("cSpectra", "Spectra Sum Closure", 720, 720);
    canvas_style(cSpectra, pad1Size, pad2Size);
    cSpectra->cd(1);
    gPad->SetLogy();

    TH1D *hSum = (TH1D *)hSpectra[0]->Clone("hSum");

    hSum->Scale((eventsperMultBin[0] / eventLoss[0]) / (minimumBiasEvents / eventLossMinBias));
    // hSum->Scale(eventsperMultBin[0] / eventLoss[0]);
    // hSum->Scale(eventLoss[0]);
    // hSum->Scale((mult_classes[1] - mult_classes[0]) / eventLoss[0]);
    // hSum->Scale(mult_classes[1] - mult_classes[0]);

    for (int imult = 1; imult < totalMultBins; imult++)
    {
        // double weight = (mult_classes[imult + 1] - mult_classes[imult]);
        // double weight = eventLoss[imult];
        double weight = (eventsperMultBin[imult] / eventLoss[imult]) / (minimumBiasEvents / eventLossMinBias);
        // double weight = (eventsperMultBin[imult] / eventLoss[imult]);
        // double weight = (mult_classes[imult + 1] - mult_classes[imult]) / eventLoss[imult];

        hSum->Add(hSpectra[imult], weight);
    }

    // hSum->Scale(1.0 / 100);
    // hSum->Scale(eventLossSum / rawEventsSum);

    hSum->SetLineColor(kRed);
    hSum->SetLineWidth(2);
    hSum->SetMarkerColor(kRed);
    hSum->SetMarkerStyle(20);
    hSum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hSum->GetYaxis()->SetTitle("d^{2}N/dydp_{T} (GeV/c)^{-1}");
    hSum->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hSum->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hSum->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hSum->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hSum->GetYaxis()->SetTitleOffset(1.3);
    hSum->Draw("pe");
    hMinBias->SetLineColor(kBlue);
    hMinBias->SetLineWidth(2);
    hMinBias->SetMarkerColor(kBlue);
    hMinBias->SetMarkerStyle(21);
    hMinBias->Draw("pe same");

    TLegend *legend = new TLegend(0.55, 0.75, 0.9, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(hSum, "Sum over mult. classes", "lp");
    legend->AddEntry(hMinBias, "Minimum Bias", "lp");
    legend->Draw();

    cSpectra->cd(2);
    gPad->SetGridy(1);
    TH1D *hRatio = (TH1D *)hSum->Clone("hRatio");
    hRatio->Divide(hMinBias);
    hRatio->SetLineColor(kBlack);
    hRatio->SetLineWidth(2);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hRatio->GetYaxis()->SetTitle("Sum / MinBias");
    hRatio->GetYaxis()->SetTitleSize(0.03 / pad2Size);
    hRatio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatio->GetYaxis()->SetLabelSize(0.03 / pad2Size);
    hRatio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hRatio->GetYaxis()->SetTitleOffset(0.55);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->SetMinimum(0.55);
    hRatio->SetMaximum(1.45);
    hRatio->Draw("pe");

    cSpectra->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/SpectraClosure/SpectraClosure.png");
}

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

    pad1->SetTicks(1, 1);
    pad2->SetTicks(1, 1);
}