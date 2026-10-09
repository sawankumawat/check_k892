#include <iostream>
#include "../src/style.h"

TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const std::string &name);
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void RatioMinBiasTrigger()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    TFile *fInputTrigger = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults_MorePhiBins.root"); // Trigger
    TFile *fInputMinBias = OpenFile("/home/sawan/Storage/check_k892/data/doublePhi/758805.root");

    THnSparseF *hUnlikeTrigger = GetHisto<THnSparseF>(fInputTrigger, "doublephimeson/SEMassDoublePhi");
    THnSparseF *hUnlikeMinBias = GetHisto<THnSparseF>(fInputMinBias, "doublephimeson/SEMassDoublePhi");

    int lowpTTrigger = hUnlikeTrigger->GetAxis(1)->FindBin(6.0 + 0.001);
    int highpTTrigger = hUnlikeTrigger->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowpTMinBias = hUnlikeMinBias->GetAxis(1)->FindBin(6.0 + 0.001);
    int highpTMinBias = hUnlikeMinBias->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowDeltaMTrigger = hUnlikeTrigger->GetAxis(6)->FindBin(0.0 + 0.00001);
    int highDeltaMTrigger = hUnlikeTrigger->GetAxis(6)->FindBin(0.005 - 0.00001);

    int lowDeltaMMinBias = hUnlikeMinBias->GetAxis(6)->FindBin(0.0 + 0.00001);
    int highDeltaMMinBias = hUnlikeMinBias->GetAxis(6)->FindBin(0.005 - 0.00001);

    hUnlikeTrigger->GetAxis(1)->SetRange(lowpTTrigger, highpTTrigger);
    // hUnlikeTrigger->GetAxis(6)->SetRange(lowDeltaMTrigger, highDeltaMTrigger);

    hUnlikeMinBias->GetAxis(1)->SetRange(lowpTMinBias, highpTMinBias);
    // hUnlikeMinBias->GetAxis(6)->SetRange(lowDeltaMMinBias, highDeltaMMinBias);

    TH1D *hInvMassTrigger = hUnlikeTrigger->Projection(0, "E");
    SetHistoQA(hInvMassTrigger);
    hInvMassTrigger->Rebin(8);
    hInvMassTrigger->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hInvMassTrigger->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hInvMassTrigger->GetBinWidth(1) * 1000));

    TH1D *hInvMassMinBias = hUnlikeMinBias->Projection(0, "E");
    SetHistoQA(hInvMassMinBias);
    hInvMassMinBias->Rebin(8);
    hInvMassMinBias->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hInvMassMinBias->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hInvMassMinBias->GetBinWidth(1) * 1000));

    int binIntegralLow = hInvMassTrigger->FindBin(2.8 + 0.0001);
    int binIntegralHigh = hInvMassTrigger->FindBin(2.9 - 0.0001);

    int binIntegralLow2 = hInvMassTrigger->FindBin(2.45 + 0.0001);
    int binIntegralHigh2 = hInvMassTrigger->FindBin(2.55 - 0.0001);

    double integralTrigger = hInvMassTrigger->Integral(binIntegralLow, binIntegralHigh);
    double integralMinBias = hInvMassMinBias->Integral(binIntegralLow, binIntegralHigh);
    double integralTrigger2 = hInvMassTrigger->Integral(binIntegralLow2, binIntegralHigh2);
    double integralMinBias2 = hInvMassMinBias->Integral(binIntegralLow2, binIntegralHigh2);
    double ratio = integralTrigger / integralMinBias;
    double ratio2 = integralTrigger2 / integralMinBias2;
    double averageRatio = (ratio + ratio2) / 2.0;
    // hInvMassTrigger->Scale(1. / averageRatio);
    cout << "Average ratio: " << averageRatio << endl;

    TH1D *hRatio = (TH1D *)hInvMassTrigger->Clone("hRatio");
    hRatio->Divide(hInvMassMinBias);

    TCanvas *cInvMass = new TCanvas("cInvMass", "Invariant Mass", 720, 720);
    SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cInvMass, pad1Size, pad2Size);
    cInvMass->cd(1);
    hInvMassTrigger->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hInvMassTrigger->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hInvMassTrigger->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hInvMassTrigger->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hInvMassTrigger->GetYaxis()->SetTitleOffset(1.90 * pad1Size);
    hInvMassTrigger->SetLineColor(kRed);
    hInvMassTrigger->SetMarkerColor(kRed);
    hInvMassTrigger->SetMarkerStyle(20);
    hInvMassTrigger->SetMinimum(0.09e6);
    hInvMassTrigger->SetMaximum(0.218e6);
    hInvMassTrigger->Draw("pe");
    hInvMassMinBias->SetLineColor(kBlue);
    hInvMassMinBias->SetMarkerColor(kBlue);
    hInvMassMinBias->SetMarkerStyle(21);
    hInvMassMinBias->Draw("pe same");

    TLegend *legend = new TLegend(0.5, 0.68, 0.9, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.045);
    legend->SetHeader("#phi#phi raw invariant mass (2026 data)");
    legend->AddEntry(hInvMassTrigger, "Triggered data", "p");
    legend->AddEntry(hInvMassMinBias, "MinBias data", "p");
    legend->Draw();

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.045);
    latex->DrawLatex(0.57, 0.61, "#it{p}_{T}^{#phi#phi} > 6 GeV/#it{c}");

    cInvMass->cd(2);
    SetHistoQA(hRatio);
    hRatio->GetYaxis()->SetTitleSize(0.025 / pad2Size);
    hRatio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hRatio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hRatio->GetYaxis()->SetTitleOffset(2.8 * pad2Size);
    hRatio->GetXaxis()->SetTitleOffset(1.1);
    hRatio->GetYaxis()->SetTitle("Trigger / MinBias");
    hRatio->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMinimum(0.85);
    hRatio->SetMaximum(1.15);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->Draw("pe");

    TLine *line = new TLine(2.4, 0.9, 3.0, 0.9);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw();

    // cInvMass->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Ratio_TriggerToMinBias.png");
}

//==============End of the main code==================

TFile *OpenFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (f->IsZombie())
    {
        cout << "Error: File not found: " << path << endl;
        return nullptr;
    }
    return f;
}

template <typename T>
T *GetHisto(TFile *f, const std::string &name)
{
    T *histo = dynamic_cast<T *>(f->Get(name.c_str()));

    if (!histo)
    {
        std::cout << "Error: histo " << name
                  << " not found in file " << f->GetName() << std::endl;
        return nullptr;
    }

    return histo;
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.08, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.03);
    pad2->SetRightMargin(0.03);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}