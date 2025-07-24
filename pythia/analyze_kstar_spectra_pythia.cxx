#include <iostream>
#include "../macro/src/style.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void analyze_kstar_spectra_pythia()
{
    gStyle->SetOptStat(0);

    TFile *fpp13000GeV = new TFile("100M_events/kstar_pp_13TeV.root");
    TFile *fpp136000GeV = new TFile("100M_events/kstar_pp_13.6TeV.root");

    if (fpp13000GeV->IsZombie() || fpp136000GeV->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    TH1F *hmultpp13 = (TH1F *)fpp13000GeV->Get("hmult");
    TH1F *hmultpp136 = (TH1F *)fpp136000GeV->Get("hmult");
    TH1F *spectrapp13 = (TH1F *)fpp13000GeV->Get("hkstarpp3");
    TH1F *spectrapp136 = (TH1F *)fpp136000GeV->Get("hkstarpp3");

    if (hmultpp13 == nullptr || hmultpp136 == nullptr || spectrapp13 == nullptr || spectrapp136 == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    // Now lets the use the binning as in the data histogra
    double ptbins[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    int nptbins = sizeof(ptbins) / sizeof(double) - 1;

    spectrapp13->Scale(1. / spectrapp13->GetEntries());
    spectrapp136->Scale(1. / spectrapp136->GetEntries());

    TH1F *ratio_rebin = new TH1F("", "", nptbins, ptbins);
    TH1F *spectrapp136_rebin = (TH1F *)spectrapp136->Rebin(nptbins, "spectrapp136_rebin", ptbins);
    TH1F *spectrapp13_rebin = (TH1F *)spectrapp13->Rebin(nptbins, "spectrapp13_rebin", ptbins);
    double deviation_below10GeV = 0.0;
    double deviation_above10GeV = 0.0;
    for (int i = 0; i < nptbins; i++)
    {
        int binmin = spectrapp13->FindBin(ptbins[i] + 0.0001);
        int binmax = spectrapp13->FindBin(ptbins[i + 1] - 0.0001);

        double content13 = spectrapp13->Integral(binmin, binmax);
        double content136 = spectrapp136->Integral(binmin, binmax);
        spectrapp13_rebin->SetBinContent(i + 1, content13 / (ptbins[i + 1] - ptbins[i]));
        spectrapp136_rebin->SetBinContent(i + 1, content136 / (ptbins[i + 1] - ptbins[i]));

        // Get errors for the integrals
        double error13 = 0.0;
        double error136 = 0.0;
        spectrapp13->IntegralAndError(binmin, binmax, error13);
        spectrapp136->IntegralAndError(binmin, binmax, error136);

        ratio_rebin->SetBinContent(i + 1, content136 / content13);
        double ratioError = sqrt(pow(error136 / content13, 2) + pow((error13 * content136) / (content13 * content13), 2));
        ratio_rebin->SetBinError(i + 1, ratioError);
        if (ptbins[i] >= 10)
            deviation_above10GeV += abs(content136 / content13 - 1);
        else
            deviation_below10GeV += abs(content136 / content13 - 1);
    }
    // deviation /= nptbins;
    deviation_above10GeV /= 3;
    deviation_below10GeV /= 25;
    cout << "Deviation above 10 GeV/c %: " << deviation_above10GeV * 100 << endl;
    cout << "Deviation below 10 GeV/c %: " << deviation_below10GeV * 100 << endl;

    TCanvas *cspectra = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cspectra, 0.14, 0.02, 0.02, 0.14);
    double pad1Size, pad2Size;
    canvas_style(cspectra, pad1Size, pad2Size);

    cspectra->cd(1);
    gPad->SetLogy();
    SetHistoQA(spectrapp13_rebin);
    SetHistoQA(spectrapp136_rebin);
    spectrapp13_rebin->SetLineColor(kRed);
    spectrapp13_rebin->SetMarkerColor(kRed);
    spectrapp13_rebin->SetMarkerStyle(20);
    spectrapp136_rebin->SetMarkerStyle(25);
    spectrapp136_rebin->SetMarkerColor(kBlue);
    spectrapp136_rebin->SetLineColor(kBlue);

    spectrapp13_rebin->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    spectrapp13_rebin->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/#it{c})^{-1}]");
    spectrapp13_rebin->GetYaxis()->SetTitleOffset(1.4);
    spectrapp13_rebin->GetXaxis()->SetTitleSize(0.042 / pad1Size);
    spectrapp13_rebin->GetYaxis()->SetTitleSize(0.044 / pad1Size);
    spectrapp13_rebin->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    spectrapp13_rebin->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    spectrapp13_rebin->GetXaxis()->SetRangeUser(0, 10);
    spectrapp13_rebin->Draw("pe");
    spectrapp136_rebin->Draw("pe same");

    TLegend *leg = new TLegend(0.52, 0.55, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04 / pad1Size);
    leg->AddEntry((TObject *)0, "Pythia simulations", "");
    leg->AddEntry((TObject *)0, "K*^{0} (892) spectra", "");
    leg->AddEntry(spectrapp136_rebin, "pp 13.6 TeV", "p");
    leg->AddEntry(spectrapp13_rebin, "pp 13 TeV", "p");
    leg->Draw();

    cspectra->cd(2);
    // TH1F *histratio = (TH1F *)spectrapp136->Clone();
    // histratio->Divide(spectrapp13); // ratio
    SetHistoQA(ratio_rebin);
    ratio_rebin->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    ratio_rebin->GetYaxis()->SetTitle("#frac{pp 13.6 TeV}{pp 13 TeV}");
    ratio_rebin->GetYaxis()->SetTitleOffset(0.66);
    ratio_rebin->GetXaxis()->SetTitleOffset(1.1);
    ratio_rebin->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    ratio_rebin->GetYaxis()->SetTitleSize(0.037 / pad2Size);
    ratio_rebin->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    ratio_rebin->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    ratio_rebin->GetXaxis()->SetRangeUser(0, 10);
    ratio_rebin->GetYaxis()->SetRangeUser(0.945, 1.07);
    ratio_rebin->GetXaxis()->SetNdivisions(510);
    ratio_rebin->GetYaxis()->SetNdivisions(505);
    ratio_rebin->Draw("pe");
    TLine *line = new TLine(0, 1, 10, 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    cspectra->SaveAs("spectra_kstar_ratio_pythia.png");
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.05);
    pad2->SetRightMargin(0.05);
    pad2->SetBottomMargin(0.35);
    pad1->SetLeftMargin(0.175);
    pad2->SetLeftMargin(0.175);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.003);
    pad2->SetTopMargin(0.04);
    
    // Set ticks on individual pads
    pad1->SetTicks(1, 1);
    pad2->SetTicks(1, 1);
}
