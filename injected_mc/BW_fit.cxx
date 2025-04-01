#include <iostream>
using namespace std;
#include "style.h"
// #include "../macro/src/fitfunc.h"
#include "../macro/src/common_glue.h"

Double_t RelativisticBW(double *x, double *par)
{
    double norm = par[0];
    double mass = par[1];
    double width = par[2];

    Int_t j1 = 0; // spin
    double n1 = (2.0 * j1 + 1.0) / 2.0;

    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = mass * mass - 4 * (0.4976 * 0.4976);

    double mass_dep_width = width * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = norm * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));

    return fit;
}

Double_t RBW_massDepWidth(double *x, double *par)
{
    double norm = par[0];
    double mass = par[1];
    double width = par[2];

    Int_t j1 = 0; // spin
    double n1 = (2.0 * j1 + 1.0) / 2.0;

    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = mass * mass - 4 * (0.4976 * 0.4976);

    double mass_dep_width = width * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = norm * mass * mass_dep_width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * mass_dep_width, 2));

    return fit;
}

Double_t fitGaussian(Double_t *x, Double_t *par)
{
    double mass = par[1];
    double width = par[2];
    double norm = par[0];
    double fit = norm * TMath::Exp(-0.5 * TMath::Power((x[0] - mass) / width, 2));
    return fit;
}

void BW_fit()
{
    gStyle->SetOptStat(0);
    // gStyle->SetOptStat(1110);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1111);

    TFile *f = new TFile("/home/sawan/check_k892/mc/LHC24l1/356039.root", "read");

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    TH1F *recMass = (TH1F *)f->Get("higher-mass-resonances/hMChists/Recf1710_mass");
    TH1F *genMass = (TH1F *)f->Get("higher-mass-resonances/hMChists/Genf1710_mass");
    if (recMass == nullptr || genMass == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.15, 0.05, 0.05, 0.15);
    SetHistostyle2(recMass);
    recMass->SetTitle("Rec f_{0}(1710) Mass");
    recMass->GetYaxis()->SetTitle("Counts");
    recMass->GetXaxis()->SetRangeUser(1.49, 2.01);
    recMass->Draw("pe");

    TLatex lat2;
    lat2.SetNDC();
    lat2.SetTextSize(0.03);
    lat2.SetTextFont(42);
    lat2.DrawLatex(0.58, 0.4, "Reconstructed f_{0}(1710) Mass");

    cout << "bin width is " << recMass->GetBinWidth(1) << endl;

    // TF1 *fit = new TF1("fit", RelativisticBW, 1.5, 1.9, 3);
    // fit->SetParNames("Norm", "Mass", "Width");
    // fit->SetParameter(0, 100);
    // fit->SetParameter(1, f1710Mass);
    // fit->SetParameter(2, f1710Width);
    // fit->SetParLimits(0, 0.0, 200.0);
    // fit->SetParLimits(2, 0.0, 2.0);
    // recMass->Fit("fit", "REBMS0");
    // fit->Draw("same");

    TF1 *fitgaus = new TF1("fitgaus", fitGaussian, 1.709, 1.76, 3);
    fitgaus->SetParNames("Norm", "Mass", "Width");
    fitgaus->SetParameter(0, 1000);
    fitgaus->SetParameter(1, f1710Mass);
    fitgaus->SetParameter(2, f1710Width);
    // fitgaus->SetParLimits(0, 0.0, 200.0);
    fitgaus->SetParLimits(2, 0.0, 2.0);
    recMass->Fit("fitgaus", "REBMS");

    // TPaveStats *st = (TPaveStats *)fit->FindObject("stats");
    // st->SetX1NDC(0.6);
    // st->SetX2NDC(0.99);
    // st->SetY1NDC(0.4);
    // st->SetY2NDC(0.9);
    // st->Draw();

    c->SaveAs("/home/sawan/Music/reconstructed_mass.png");
    // c->Close();

    TCanvas *cgen = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cgen, 0.15, 0.05, 0.05, 0.15);
    SetHistostyle2(genMass);
    genMass->GetYaxis()->SetTitle("Counts");
    genMass->GetXaxis()->SetRangeUser(1.49, 2.01);
    // genMass->GetXaxis()->SetRangeUser(1.19, 1.89);
    // genMass->GetXaxis()->SetRangeUser(1.09, 1.69);
    genMass->Draw();
    float generatedmass = genMass->GetMean();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.SetTextFont(42);
    lat.DrawLatex(0.58, 0.6, "Generated f_{0}(1710) Mass");

    cout << "bin width of generated mass is " << genMass->GetBinWidth(1) << endl;

    // cgen->SaveAs("/home/sawan/Music/generated_mass.png");
    cgen->Close();

    // lets check with differential pT
    THnSparseF *RecMasspt = (THnSparseF *)f->Get("higher-mass-resonances/hMChists/Recf1710_pt1");
    if (RecMasspt == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    // axis: multiplicity, pt, mass
    double ptbins[] = {0, 1, 2, 3, 5, 8, 12, 20};
    // int ptbins[] = {0, 20};
    int size = sizeof(ptbins) / sizeof(ptbins[0]);
    TH1F *hmass[size];
    TCanvas *cdiff = new TCanvas("", "", 1080, 720);
    cdiff->Divide(3, 2);
    float fitrangelow[] = {1.712, 1.712, 1.712, 1.712, 1.712, 1.712, 1.701};
    float fitrangehigh[] = {1.755, 1.755, 1.755, 1.755, 1.755, 1.755, 1.765};
    float f1710masses[6];
    float f1710mass_errors[6];
    float f1710widths[6];
    float f1710width_errors[6];
    float f1710mass_difference[6];

    for (int i = 0; i < size - 1; i++)
    {
        cdiff->cd(i + 1);
        int lowptbin = RecMasspt->GetAxis(1)->FindBin(ptbins[i]);
        int highptbin = RecMasspt->GetAxis(1)->FindBin(ptbins[i + 1]);
        RecMasspt->GetAxis(1)->SetRange(lowptbin, highptbin);
        hmass[i] = (TH1F *)RecMasspt->Projection(2);
        hmass[i]->SetName(Form("hmass_%d", i));
        hmass[i]->GetXaxis()->SetRangeUser(1.49, 2.01);
        hmass[i]->Draw();
        TF1 *fitgaus1 = new TF1("fitgaus1", fitGaussian, fitrangelow[i], fitrangehigh[i], 3);
        fitgaus1->SetParNames("Norm", "Mass", "Width");
        fitgaus1->SetParameter(0, (i + 1) * 50);
        fitgaus1->SetParameter(1, f1710Mass);
        fitgaus1->SetParameter(2, f1710Width);
        fitgaus1->SetParLimits(2, 0.001, 0.075);
        hmass[i]->Fit("fitgaus1", "REBMS");
        f1710masses[i] = fitgaus1->GetParameter(1);
        f1710mass_errors[i] = fitgaus1->GetParError(1);
        f1710widths[i] = fitgaus1->GetParameter(2);
        f1710width_errors[i] = fitgaus1->GetParError(2);
        f1710mass_difference[i] = f1710masses[i] - generatedmass;
        TLatex lat3;
        lat3.SetNDC();
        lat3.SetTextSize(0.045);
        lat3.SetTextFont(42);
        lat3.DrawLatex(0.58, 0.3, Form("p_{T} %.0f - %.0f GeV/c", ptbins[i], ptbins[i + 1]));
    }

    cdiff->SaveAs("/home/sawan/Music/mass_plot_allpt.png");

    TCanvas *cmassfit = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cmassfit, 0.20, 0.05, 0.05, 0.15);
    TH1F *hmassfit = new TH1F("hmassfit", "f_{0}(1710) Mass vs p_{T}", size-1, ptbins);
    SetHistoQA(hmassfit);
    hmassfit->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    hmassfit->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hmassfit->GetYaxis()->SetTitleOffset(2.0);
    hmassfit->SetMarkerStyle(20);
    hmassfit->SetMarkerSize(1.5);
    hmassfit->SetMarkerColor(kBlue);
    hmassfit->SetLineColor(kBlue);
    hmassfit->SetLineWidth(2);
    for (int i = 0; i < size - 1; i++)
    {
        hmassfit->SetBinContent(i + 1, f1710masses[i]);
        hmassfit->SetBinError(i + 1, f1710mass_errors[i]);
    }
    hmassfit->GetYaxis()->SetRangeUser(1.725, 1.749);
    hmassfit->Draw("pe");
    TLine *lgenmass = new TLine(0, generatedmass, 20, generatedmass);
    lgenmass->SetLineColor(kRed);
    lgenmass->SetLineStyle(2);
    lgenmass->SetLineWidth(2);
    lgenmass->Draw("same");

    TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry(hmassfit, "Reconstructed Mass", "lpe");
    leg->AddEntry(lgenmass, "Generated Mass", "l");
    leg->Draw();

    cmassfit->SaveAs("/home/sawan/Music/mass_pt.png");

    // TCanvas *cwidthfit = new TCanvas("", "", 720, 720);
    // SetCanvasStyle(cwidthfit, 0.15, 0.05, 0.05, 0.15);
    // TH1F *hwidthfit = new TH1F("hwidthfit", "f_{0}(1710) Width vs p_{T}", size - 1, 0, size - 1);
    // hwidthfit->GetYaxis()->SetTitle("Width (GeV/c^{2})");
    // hwidthfit->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // hwidthfit->SetMarkerStyle(20);
    // hwidthfit->SetMarkerSize(1.5);
    // hwidthfit->SetMarkerColor(kRed);
    // hwidthfit->SetLineColor(kRed);
    // hwidthfit->SetLineWidth(2);
    // for (int i = 0; i < size - 1; i++)
    // {
    //     hwidthfit->SetBinContent(i + 1, f1710widths[i]);
    //     hwidthfit->SetBinError(i + 1, f1710width_errors[i]);
    // }
    // // hwidthfit->GetYaxis()->SetRangeUser(0.02, 0.06);
    // hwidthfit->Draw("pe");

    TCanvas *cmassdiff = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cmassdiff, 0.15, 0.05, 0.05, 0.15);
    TH1F *hmassdiff = new TH1F("hmassdiff", "f_{0}(1710) Mass Difference vs p_{T}", size - 1, ptbins);
    SetHistoQA(hmassdiff);
    hmassdiff->GetYaxis()->SetTitle("Mass Difference (MeV/c^{2})");
    hmassdiff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hmassdiff->SetMarkerStyle(20);
    hmassdiff->SetMarkerSize(1.8);
    hmassdiff->SetMarkerColor(2);
    hmassdiff->SetLineColor(2);
    hmassdiff->SetLineWidth(2);
    for (int i = 0; i < size - 1; i++)
    {
        hmassdiff->SetBinContent(i + 1, f1710mass_difference[i] * 1000);
        hmassdiff->SetBinError(i + 1, 0);
    }
    hmassdiff->GetYaxis()->SetRangeUser(0.00, 7.99);
    hmassdiff->Draw("pe");

    TLatex lat4;
    lat4.SetNDC();
    lat4.SetTextSize(0.04);
    lat4.SetTextFont(42);
    lat4.DrawLatex(0.5, 0.8, "f_{0}(1710) resonance");
    lat4.DrawLatex(0.5, 0.7, "Gen. Mass - Rec. Mass");

    cmassdiff->SaveAs("/home/sawan/Music/mass_difference_pt.png");
}

// // TF1 *fit2 = new TF1("fit2", RBW_massDepWidth, 1.5, 1.9, 3);
// TF1 *fit2 = new TF1("fit2", RBW_massDepWidth, 1.2, 1.8, 3);
// fit2->SetParameter(0, 1);
// fit2->SetParameter(1, f1525Mass);
// fit2->SetParameter(2, f1525Width);
// recMass->Fit("fit2", "REBMS0");
// // fit2->Draw("same");