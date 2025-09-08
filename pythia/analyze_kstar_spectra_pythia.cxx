#include <iostream>
#include "../macro/src/style.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void analyze_kstar_spectra_pythia()
{
    gStyle->SetOptStat(0);

    // TFile *fpp13000GeV = new TFile("100M_events/kstar_pp_13TeV.root");
    // TFile *fpp136000GeV = new TFile("100M_events/kstar_pp_13.6TeV.root");

    TFile *fpp13000GeV = new TFile("500M_events_mult/13TeV/pythiaYield13.root");
    TFile *fpp136000GeV = new TFile("500M_events_mult/13p6TeV/pythiaYield13p6.root");

    if (fpp13000GeV->IsZombie() || fpp136000GeV->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    // TH1F *hmultpp13 = (TH1F *)fpp13000GeV->Get("hmult");
    // TH1F *hmultpp136 = (TH1F *)fpp136000GeV->Get("hmult");
    // TH1F *spectrapp13 = (TH1F *)fpp13000GeV->Get("hkstarpp3");
    // TH1F *spectrapp136 = (TH1F *)fpp136000GeV->Get("hkstarpp3");

    // TH1F *hmultpp13 = (TH1F *)fpp13000GeV->Get("mult_dist");
    // TH1F *hmultpp136 = (TH1F *)fpp136000GeV->Get("mult_dist");
    // TH1F *spectrapp13 = (TH1F *)fpp13000GeV->Get("kstar_ptspectra");
    // TH1F *spectrapp136 = (TH1F *)fpp136000GeV->Get("kstar_ptspectra");
    // TH1F *spectraycut13 = (TH1F *)fpp13000GeV->Get("kstar_dist_ycut");
    // TH1F *spectraycut136 = (TH1F *)fpp136000GeV->Get("kstar_dist_ycut");
    // TH1F *spectraPhiycut13 = (TH1F *)fpp13000GeV->Get("phi_ptspectra");
    // TH1F *spectraPhiycut136 = (TH1F *)fpp136000GeV->Get("phi_ptspectra");

    TH1F *hmultpp13 = (TH1F *)fpp13000GeV->Get("hmultV0");
    TH1F *hmultpp136 = (TH1F *)fpp136000GeV->Get("hmultFT0");
    TH1F *spectrapp13 = (TH1F *)fpp13000GeV->Get("hKstarV0");
    TH1F *spectrapp136 = (TH1F *)fpp136000GeV->Get("hKstarFT0");
    TH1F *spectraycut13 = (TH1F *)fpp13000GeV->Get("hKstarV0ycut");
    TH1F *spectraycut136 = (TH1F *)fpp136000GeV->Get("hKstarFT0ycut");
    TH1F *spectraPhiycut13 = (TH1F *)fpp13000GeV->Get("hPhiV0ycut");
    TH1F *spectraPhiycut136 = (TH1F *)fpp136000GeV->Get("hPhiFT0ycut");

    if (hmultpp13 == nullptr || hmultpp136 == nullptr || spectrapp13 == nullptr || spectrapp136 == nullptr || spectraycut13 == nullptr || spectraycut136 == nullptr || spectraPhiycut13 == nullptr || spectraPhiycut136 == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    // Now lets the use the binning as in the data histogra
    double ptbins[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    int nptbins = sizeof(ptbins) / sizeof(double) - 1;

    spectrapp13->Scale(1. / spectrapp13->GetEntries());
    spectrapp136->Scale(1. / spectrapp136->GetEntries());
    spectraycut13->Scale(1. / spectraycut13->GetEntries());
    spectraycut136->Scale(1. / spectraycut136->GetEntries());
    spectraPhiycut13->Scale(1. / spectraPhiycut13->GetEntries());
    spectraPhiycut136->Scale(1. / spectraPhiycut136->GetEntries());

    TH1F *ratioSpectraKstar = new TH1F("", "", nptbins, ptbins);
    TH1F *ratioSpectraKstarYcut = new TH1F("", "", nptbins, ptbins);
    TH1F *ratioSpectraPhiYcut = new TH1F("", "", nptbins, ptbins);

    TH1F *hkstarSpectra13p6 = (TH1F *)spectrapp136->Rebin(nptbins, "hkstarSpectra13p6", ptbins);
    TH1F *hkstarSpectra13 = (TH1F *)spectrapp13->Rebin(nptbins, "hkstarSpectra13", ptbins);
    TH1F *hphiSpectra13p6 = (TH1F *)spectraPhiycut136->Rebin(nptbins, "hphiSpectra13p6", ptbins);
    TH1F *hphiSpectra13 = (TH1F *)spectraPhiycut13->Rebin(nptbins, "hphiSpectra13", ptbins);
    TH1F *hkstarSpectraYcut13p6 = (TH1F *)spectraycut136->Rebin(nptbins, "hkstarSpectraYcut13p6", ptbins);
    TH1F *hkstarSpectraYcut13 = (TH1F *)spectraycut13->Rebin(nptbins, "hkstarSpectraYcut13", ptbins);
    double deviation_below10GeV = 0.0;
    double deviation_above10GeV = 0.0;
    for (int i = 0; i < nptbins; i++)
    {
        int binmin = spectrapp13->FindBin(ptbins[i] + 0.0001);
        int binmax = spectrapp13->FindBin(ptbins[i + 1] - 0.0001);

        double yieldKstar13 = spectrapp13->Integral(binmin, binmax);
        double yieldKstar13p6 = spectrapp136->Integral(binmin, binmax);
        double yieldKstarYcut13 = spectraycut13->Integral(binmin, binmax);
        double yieldKstarYcut13p6 = spectraycut136->Integral(binmin, binmax);
        double yieldPhi13 = spectraPhiycut13->Integral(binmin, binmax);
        double yieldPhi13p6 = spectraPhiycut136->Integral(binmin, binmax);

        hkstarSpectra13->SetBinContent(i + 1, yieldKstar13 / (ptbins[i + 1] - ptbins[i]));
        hkstarSpectra13p6->SetBinContent(i + 1, yieldKstar13p6 / (ptbins[i + 1] - ptbins[i]));
        hkstarSpectraYcut13->SetBinContent(i + 1, yieldKstarYcut13 / (ptbins[i + 1] - ptbins[i]));
        hkstarSpectraYcut13p6->SetBinContent(i + 1, yieldKstarYcut13p6 / (ptbins[i + 1] - ptbins[i]));
        hphiSpectra13->SetBinContent(i + 1, yieldPhi13 / (ptbins[i + 1] - ptbins[i]));
        hphiSpectra13p6->SetBinContent(i + 1, yieldPhi13p6 / (ptbins[i + 1] - ptbins[i]));

        // Get errors for the integrals
        double yieldErrKstar13 = 0.0;
        double yieldErrKstar13p6 = 0.0;
        double yieldErrKstarYcut13 = 0.0;
        double yieldErrKstarYcut13p6 = 0.0;
        double yieldErrPhi13 = 0.0;
        double yieldErrPhi13p6 = 0.0;

        spectrapp13->IntegralAndError(binmin, binmax, yieldErrKstar13);
        spectrapp136->IntegralAndError(binmin, binmax, yieldErrKstar13p6);
        spectraycut13->IntegralAndError(binmin, binmax, yieldErrKstarYcut13);
        spectraycut136->IntegralAndError(binmin, binmax, yieldErrKstarYcut13p6);
        spectraPhiycut13->IntegralAndError(binmin, binmax, yieldErrPhi13);
        spectraPhiycut136->IntegralAndError(binmin, binmax, yieldErrPhi13p6);

        ratioSpectraKstar->SetBinContent(i + 1, yieldKstar13p6 / yieldKstar13);
        ratioSpectraKstarYcut->SetBinContent(i + 1, yieldKstarYcut13p6 / yieldKstarYcut13);
        ratioSpectraPhiYcut->SetBinContent(i + 1, yieldPhi13p6 / yieldPhi13);

        // Error propagation for ratio
        double yieldKstarRatioErr = sqrt(pow(yieldErrKstar13p6 / yieldKstar13, 2) + pow((yieldErrKstar13 * yieldKstar13p6) / (yieldKstar13 * yieldKstar13), 2));
        double yieldKstarYcutRatioErr = sqrt(pow(yieldErrKstarYcut13p6 / yieldKstarYcut13, 2) + pow((yieldErrKstarYcut13 * yieldKstarYcut13p6) / (yieldKstarYcut13 * yieldKstarYcut13), 2));
        double yieldPhiRatioErr = sqrt(pow(yieldErrPhi13p6 / yieldPhi13, 2) + pow((yieldErrPhi13 * yieldPhi13p6) / (yieldPhi13 * yieldPhi13), 2));

        ratioSpectraKstar->SetBinError(i + 1, yieldKstarRatioErr);
        ratioSpectraKstarYcut->SetBinError(i + 1, yieldKstarYcutRatioErr);
        ratioSpectraPhiYcut->SetBinError(i + 1, yieldPhiRatioErr);

        if (ptbins[i] >= 10)
            deviation_above10GeV += abs(yieldKstar13p6 / yieldKstar13 - 1);
        else
            deviation_below10GeV += abs(yieldKstar13p6 / yieldKstar13 - 1);
    }
    // deviation /= nptbins;
    deviation_above10GeV /= 3;
    deviation_below10GeV /= 25;
    cout << "Deviation above 10 GeV/c %: " << deviation_above10GeV * 100 << endl;
    cout << "Deviation below 10 GeV/c %: " << deviation_below10GeV * 100 << endl;

    TCanvas *cspectraKstar = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cspectraKstar, 0.14, 0.02, 0.02, 0.14);
    double pad1Size, pad2Size;
    canvas_style(cspectraKstar, pad1Size, pad2Size);

    cspectraKstar->cd(1);
    gPad->SetLogy();
    SetHistoQA(hkstarSpectra13);
    SetHistoQA(hkstarSpectra13p6);
    hkstarSpectra13->SetLineColor(kRed);
    hkstarSpectra13->SetMarkerColor(kRed);
    hkstarSpectra13->SetMarkerStyle(20);
    hkstarSpectra13p6->SetMarkerStyle(25);
    hkstarSpectra13p6->SetMarkerColor(kBlue);
    hkstarSpectra13p6->SetLineColor(kBlue);

    hkstarSpectra13->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    hkstarSpectra13->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/#it{c})^{-1}]");
    hkstarSpectra13->GetYaxis()->SetTitleOffset(1.4);
    hkstarSpectra13->GetXaxis()->SetTitleSize(0.042 / pad1Size);
    hkstarSpectra13->GetYaxis()->SetTitleSize(0.044 / pad1Size);
    hkstarSpectra13->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hkstarSpectra13->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hkstarSpectra13->GetXaxis()->SetRangeUser(0, 15);
    hkstarSpectra13->Draw("pe");
    hkstarSpectra13p6->Draw("pe same");

    TLegend *leg = new TLegend(0.52, 0.55, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03 / pad1Size);
    leg->AddEntry((TObject *)0, "Pythia simulations", "");
    leg->AddEntry((TObject *)0, "K*^{0}(892) spectra", "");
    leg->AddEntry(hkstarSpectra13p6, "pp 13.6 TeV", "p");
    leg->AddEntry(hkstarSpectra13, "pp 13 TeV", "p");
    leg->Draw();

    cspectraKstar->cd(2);
    // TH1F *histratio = (TH1F *)spectrapp136->Clone();
    // histratio->Divide(spectrapp13); // ratio
    gPad->SetGridy();
    SetHistoQA(ratioSpectraKstar);
    ratioSpectraKstar->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    ratioSpectraKstar->GetYaxis()->SetTitle("#frac{pp 13.6 TeV}{pp 13 TeV}");
    ratioSpectraKstar->GetYaxis()->SetTitleOffset(0.66);
    ratioSpectraKstar->GetXaxis()->SetTitleOffset(1.1);
    ratioSpectraKstar->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    ratioSpectraKstar->GetYaxis()->SetTitleSize(0.037 / pad2Size);
    ratioSpectraKstar->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    ratioSpectraKstar->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    ratioSpectraKstar->GetXaxis()->SetRangeUser(0, 15);
    ratioSpectraKstar->GetYaxis()->SetRangeUser(0.945, 1.07);
    ratioSpectraKstar->GetXaxis()->SetNdivisions(510);
    ratioSpectraKstar->GetYaxis()->SetNdivisions(505);
    ratioSpectraKstar->Draw("pe");
    TLine *line = new TLine(0, 1, 15, 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
    cspectraKstar->SaveAs("pythiaKstarSpectra.png");

    TCanvas *cspectraKstarYcut = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cspectraKstarYcut, 0.14, 0.02, 0.02, 0.14);
    canvas_style(cspectraKstarYcut, pad1Size, pad2Size);
    cspectraKstarYcut->cd(1);
    gPad->SetLogy();
    SetHistoQA(hkstarSpectraYcut13);
    SetHistoQA(hkstarSpectraYcut13p6);
    hkstarSpectraYcut13->SetLineColor(kRed);
    hkstarSpectraYcut13->SetMarkerColor(kRed);
    hkstarSpectraYcut13->SetMarkerStyle(20);
    hkstarSpectraYcut13p6->SetMarkerStyle(25);
    hkstarSpectraYcut13p6->SetMarkerColor(kBlue);
    hkstarSpectraYcut13p6->SetLineColor(kBlue);
    hkstarSpectraYcut13->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    hkstarSpectraYcut13->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/#it{c})^{-1}]");
    hkstarSpectraYcut13->GetYaxis()->SetTitleOffset(1.4);
    hkstarSpectraYcut13->GetXaxis()->SetTitleSize(0.042 / pad1Size);
    hkstarSpectraYcut13->GetYaxis()->SetTitleSize(0.044 / pad1Size);
    hkstarSpectraYcut13->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hkstarSpectraYcut13->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hkstarSpectraYcut13->GetXaxis()->SetRangeUser(0, 15);
    hkstarSpectraYcut13->Draw("pe");
    hkstarSpectraYcut13p6->Draw("pe same");
    leg->Clear();
    leg->AddEntry((TObject *)0, "Pythia simulations", "");
    leg->AddEntry((TObject *)0, "K*^{0}(892) spectra |y|<0.5", "");
    leg->AddEntry(hkstarSpectraYcut13p6, "pp 13.6 TeV", "p");
    leg->AddEntry(hkstarSpectraYcut13, "pp 13 TeV", "p");
    leg->Draw();

    cspectraKstarYcut->cd(2);
    gPad->SetGridy();
    SetHistoQA(ratioSpectraKstarYcut);
    ratioSpectraKstarYcut->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    ratioSpectraKstarYcut->GetYaxis()->SetTitle("#frac{pp 13.6 TeV}{pp 13 TeV}");
    ratioSpectraKstarYcut->GetYaxis()->SetTitleOffset(0.66);
    ratioSpectraKstarYcut->GetXaxis()->SetTitleOffset(1.1);
    ratioSpectraKstarYcut->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    ratioSpectraKstarYcut->GetYaxis()->SetTitleSize(0.037 / pad2Size);
    ratioSpectraKstarYcut->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    ratioSpectraKstarYcut->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    ratioSpectraKstarYcut->GetXaxis()->SetRangeUser(0, 15);
    ratioSpectraKstarYcut->GetYaxis()->SetRangeUser(0.945, 1.07);
    ratioSpectraKstarYcut->GetXaxis()->SetNdivisions(510);
    ratioSpectraKstarYcut->GetYaxis()->SetNdivisions(505);
    ratioSpectraKstarYcut->Draw("pe");
    line->Draw();
    cspectraKstarYcut->SaveAs("pythiaKstarSpectraYcut.png");

    TCanvas *cspectraPhiYcut = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cspectraPhiYcut, 0.14, 0.02, 0.02, 0.14);
    canvas_style(cspectraPhiYcut, pad1Size, pad2Size);
    cspectraPhiYcut->cd(1);
    gPad->SetLogy();
    SetHistoQA(hphiSpectra13);
    SetHistoQA(hphiSpectra13p6);
    hphiSpectra13->SetLineColor(kRed);
    hphiSpectra13->SetMarkerColor(kRed);
    hphiSpectra13->SetMarkerStyle(20);
    hphiSpectra13p6->SetMarkerStyle(25);
    hphiSpectra13p6->SetMarkerColor(kBlue);
    hphiSpectra13p6->SetLineColor(kBlue);
    hphiSpectra13->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    hphiSpectra13->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/#it{c})^{-1}]");
    hphiSpectra13->GetYaxis()->SetTitleOffset(1.4);
    hphiSpectra13->GetXaxis()->SetTitleSize(0.042 / pad1Size);
    hphiSpectra13->GetYaxis()->SetTitleSize(0.044 / pad1Size);
    hphiSpectra13->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hphiSpectra13->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hphiSpectra13->GetXaxis()->SetRangeUser(0, 15);
    hphiSpectra13->Draw("pe");
    hphiSpectra13p6->Draw("pe same");
    leg->Clear();
    leg->AddEntry((TObject *)0, "Pythia simulations", "");
    leg->AddEntry((TObject *)0, "#phi(1020) spectra |y|<0.5", "");
    leg->AddEntry(hphiSpectra13p6, "pp 13.6 TeV", "p");
    leg->AddEntry(hphiSpectra13, "pp 13 TeV", "p");
    leg->Draw();

    cspectraPhiYcut->cd(2);
    gPad->SetGridy();
    SetHistoQA(ratioSpectraPhiYcut);
    ratioSpectraPhiYcut->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    ratioSpectraPhiYcut->GetYaxis()->SetTitle("#frac{pp 13.6 TeV}{pp 13 TeV}");
    ratioSpectraPhiYcut->GetYaxis()->SetTitleOffset(0.66);
    ratioSpectraPhiYcut->GetXaxis()->SetTitleOffset(1.1);
    ratioSpectraPhiYcut->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    ratioSpectraPhiYcut->GetYaxis()->SetTitleSize(0.037 / pad2Size);
    ratioSpectraPhiYcut->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    ratioSpectraPhiYcut->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    ratioSpectraPhiYcut->GetXaxis()->SetRangeUser(0, 15);
    ratioSpectraPhiYcut->GetYaxis()->SetRangeUser(0.945, 1.07);
    ratioSpectraPhiYcut->GetXaxis()->SetNdivisions(510);
    ratioSpectraPhiYcut->GetYaxis()->SetNdivisions(505);
    ratioSpectraPhiYcut->Draw("pe");
    line->Draw();
    cspectraPhiYcut->SaveAs("pythiaPhiSpectraYcut.png");

    TCanvas *cmult = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cmult, 0.14, 0.02, 0.05, 0.14);
    gPad->SetLogy();
    hmultpp13->SetLineColor(kRed);
    hmultpp13->SetMarkerColor(kRed);
    hmultpp13->SetMarkerStyle(20);
    hmultpp136->SetMarkerStyle(25);
    hmultpp136->SetMarkerColor(kBlue);
    hmultpp136->SetLineColor(kBlue);
    hmultpp13->GetXaxis()->SetTitle("Multiplicity");
    hmultpp13->GetYaxis()->SetTitle("Entries");
    hmultpp13->GetYaxis()->SetTitleOffset(1.4);
    hmultpp13->GetXaxis()->SetTitleSize(0.042);
    hmultpp13->GetYaxis()->SetTitleSize(0.044);
    hmultpp13->GetXaxis()->SetLabelSize(0.04);
    hmultpp13->GetYaxis()->SetLabelSize(0.04);
    hmultpp13->GetXaxis()->SetRangeUser(0, 320);
    hmultpp13->SetMaximum(hmultpp13->GetMaximum() * 3);
    hmultpp13->Draw("pe");
    hmultpp136->Draw("pe same");
    TLegend *leg2 = new TLegend(0.52, 0.7, 0.92, 0.92);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->AddEntry((TObject *)0, "Pythia simulations", "");
    leg2->AddEntry((TObject *)0, "Multiplicity distribution", "");
    leg2->AddEntry(hmultpp136, "pp 13.6 TeV (FT0M)", "p");
    leg2->AddEntry(hmultpp13, "pp 13 TeV (V0M)", "p");
    leg2->Draw();
    cmult->SaveAs("pythiaMultiplicity.png");
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
