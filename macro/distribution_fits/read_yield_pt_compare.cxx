#include "../src/style.h"
#include "../src/common_glue.h"
#include "TF1.h"
#include "TMath.h"
#include "../spectra/YieldMean.C"
using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void read_yield_pt_compare()
{
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/435450/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/";
    string path2 = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/";
    TFile *file1 = new TFile((path + "WidthFree/mult_0-100/Spectra/spectra.root").c_str(), "read");
    TFile *file2 = new TFile((path + "mult_0-100/Spectra/spectra.root").c_str(), "read");
    TFile *file3 = new TFile((path2 + "WidthFree/mult_0-100/Spectra/spectra.root").c_str(), "read");
    TFile *file4 = new TFile((path2 + "mult_0-100/Spectra/spectra.root").c_str(), "read");

    // TFile *file1 = new TFile((path + "mult_0-100/Spectra/spectra_pt0.root").c_str(), "read");
    // TFile *file2 = new TFile((path + "mult_0-100/Spectra/spectra_pt1.root").c_str(), "read");
    // TFile *file3 = new TFile((path + "mult_0-100/Spectra/spectra_pt2.root").c_str(), "read");
    // TFile *file4 = new TFile((path + "mult_0-100/Spectra/spectra.root").c_str(), "read");
    if (file1->IsZombie() || file2->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TString outputPath = path + "/SpectraCompare";
    if (gSystem->mkdir(outputPath, kTRUE))
    {
        std::cout << "Folder " << outputPath << " created successfully." << std::endl;
    }

    TH1F *hMass1710_1 = (TH1F *)file1->Get("hMass1710");
    TH1F *hWidth1710_1 = (TH1F *)file1->Get("hWidth1710");
    TH1F *hMass1525_1 = (TH1F *)file1->Get("hMass1525");
    TH1F *hMass1270_1 = (TH1F *)file1->Get("hMass1270");
    TH1F *hMass1320_1 = (TH1F *)file1->Get("hMass1320");
    TH1F *hYield1525Corrected_1 = (TH1F *)file1->Get("hYield1525Corrected");
    TH1F *hYield1710Corrected_1 = (TH1F *)file1->Get("hYield1710Corrected");
    TH1F *hYieldRatio_1 = (TH1F *)file1->Get("hYieldRatio");
    TH1F *hYield1525Raw_1 = (TH1F *)file1->Get("hYield1525Raw");
    TH1F *hYield1710Raw_1 = (TH1F *)file1->Get("hYield1710Raw");
    TGraphErrors *gMeanPtvsMass_1 = (TGraphErrors *)file1->Get("gMeanPt_f0f2");

    TH1F *hMass1710_2 = (TH1F *)file2->Get("hMass1710");
    TH1F *hWidth1710_2 = (TH1F *)file2->Get("hWidth1710");
    TH1F *hMass1525_2 = (TH1F *)file2->Get("hMass1525");
    TH1F *hMass1270_2 = (TH1F *)file2->Get("hMass1270");
    TH1F *hMass1320_2 = (TH1F *)file2->Get("hMass1320");
    TH1F *hYield1525Corrected_2 = (TH1F *)file2->Get("hYield1525Corrected");
    TH1F *hYield1710Corrected_2 = (TH1F *)file2->Get("hYield1710Corrected");
    TH1F *hYieldRatio_2 = (TH1F *)file2->Get("hYieldRatio");
    TH1F *hYield1525Raw_2 = (TH1F *)file2->Get("hYield1525Raw");
    TH1F *hYield1710Raw_2 = (TH1F *)file2->Get("hYield1710Raw");
    TGraphErrors *gMeanPtvsMass_2 = (TGraphErrors *)file2->Get("gMeanPt_f0f2");

    TH1F *hMass1710_3 = (TH1F *)file3->Get("hMass1710");
    TH1F *hWidth1710_3 = (TH1F *)file3->Get("hWidth1710");
    TH1F *hMass1525_3 = (TH1F *)file3->Get("hMass1525");
    TH1F *hMass1270_3 = (TH1F *)file3->Get("hMass1270");
    TH1F *hMass1320_3 = (TH1F *)file3->Get("hMass1320");
    TH1F *hYield1525Corrected_3 = (TH1F *)file3->Get("hYield1525Corrected");
    TH1F *hYield1710Corrected_3 = (TH1F *)file3->Get("hYield1710Corrected");
    TH1F *hYieldRatio_3 = (TH1F *)file3->Get("hYieldRatio");
    TH1F *hYield1525Raw_3 = (TH1F *)file3->Get("hYield1525Raw");
    TH1F *hYield1710Raw_3 = (TH1F *)file3->Get("hYield1710Raw");
    TGraphErrors *gMeanPtvsMass_3 = (TGraphErrors *)file3->Get("gMeanPt_f0f2");

    TH1F *hMass1710_4 = (TH1F *)file4->Get("hMass1710");
    TH1F *hWidth1710_4 = (TH1F *)file4->Get("hWidth1710");
    TH1F *hMass1525_4 = (TH1F *)file4->Get("hMass1525");
    TH1F *hMass1270_4 = (TH1F *)file4->Get("hMass1270");
    TH1F *hMass1320_4 = (TH1F *)file4->Get("hMass1320");
    TH1F *hYield1525Corrected_4 = (TH1F *)file4->Get("hYield1525Corrected");
    TH1F *hYield1710Corrected_4 = (TH1F *)file4->Get("hYield1710Corrected");
    TH1F *hYieldRatio_4 = (TH1F *)file4->Get("hYieldRatio");
    TH1F *hYield1525Raw_4 = (TH1F *)file4->Get("hYield1525Raw");
    TH1F *hYield1710Raw_4 = (TH1F *)file4->Get("hYield1710Raw");
    TGraphErrors *gMeanPtvsMass_4 = (TGraphErrors *)file4->Get("gMeanPt_f0f2");

    TCanvas *cMass1710 = new TCanvas("cMass1710", "Mass vs #it{p}_{T} for f_{0}(1710)", 720, 720);
    SetCanvasStyle(cMass1710, 0.18, 0.03, 0.05, 0.14);
    hMass1710_1->GetYaxis()->SetRangeUser(1.64, 1.84);
    hMass1710_1->SetMarkerStyle(20);
    hMass1710_1->SetMarkerColor(kRed);
    hMass1710_1->SetLineColor(kRed);
    hMass1710_1->Draw("pe");
    hMass1710_2->SetMarkerStyle(21);
    hMass1710_2->SetMarkerColor(kBlue);
    hMass1710_2->SetLineColor(kBlue);
    hMass1710_2->Draw("pe same");
    hMass1710_3->SetMarkerStyle(22);
    hMass1710_3->SetMarkerColor(kGreen + 2);
    hMass1710_3->SetLineColor(kGreen + 2);
    hMass1710_3->Draw("pe same");
    hMass1710_4->SetMarkerStyle(23);
    hMass1710_4->SetMarkerColor(kMagenta);
    hMass1710_4->SetLineColor(kMagenta);
    hMass1710_4->Draw("pe same");
    TLine *line1710Mass = new TLine(0, f1710Mass, 12, f1710Mass);
    line1710Mass->SetLineStyle(2);
    line1710Mass->SetLineColor(2);
    line1710Mass->Draw();
    // TLegend *leg1710Mass = new TLegend(0.55, 0.7, 0.85, 0.93);
    TLegend *leg1710Mass = new TLegend(0.5, 0.6, 0.85, 0.93);
    leg1710Mass->SetBorderSize(0);
    leg1710Mass->SetFillStyle(0);
    leg1710Mass->SetTextSize(0.03);
    leg1710Mass->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    leg1710Mass->AddEntry(hMass1710_1, "f_{0}(1710) Width Free (2023)", "pe");
    leg1710Mass->AddEntry(hMass1710_2, "f_{0}(1710) Width Fixed (2023)", "pe");
    leg1710Mass->AddEntry(hMass1710_3, "f_{0}(1710) Width Free (2022)", "pe");
    leg1710Mass->AddEntry(hMass1710_4, "f_{0}(1710) Width Fixed (2022)", "pe");
    leg1710Mass->AddEntry(line1710Mass, "PDG value", "l");
    leg1710Mass->Draw();
    cMass1710->SaveAs(outputPath + "/Mass1710.png");

    TCanvas *cWidth1710 = new TCanvas("cWidth1710", "Width vs #it{p}_{T} for f_{0}(1710)", 720, 720);
    SetCanvasStyle(cWidth1710, 0.18, 0.03, 0.05, 0.14);
    hWidth1710_1->GetYaxis()->SetRangeUser(0.05, 0.33);
    hWidth1710_1->SetMarkerStyle(20);
    hWidth1710_1->SetMarkerColor(kRed);
    hWidth1710_1->SetLineColor(kRed);
    hWidth1710_1->Draw("pe");
    hWidth1710_3->SetMarkerStyle(21);
    hWidth1710_3->SetMarkerColor(kBlue);
    hWidth1710_3->SetLineColor(kBlue);
    hWidth1710_3->Draw("pe same");
    TLine *line1710Width = new TLine(0, f1710Width, 12, f1710Width);
    line1710Width->SetLineStyle(2);
    line1710Width->SetLineColor(2);
    line1710Width->Draw();
    TLegend *leg1710Width = new TLegend(0.55, 0.75, 0.85, 0.93);
    leg1710Width->SetBorderSize(0);
    leg1710Width->SetFillStyle(0);
    leg1710Width->SetTextSize(0.03);
    leg1710Width->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    leg1710Width->AddEntry(hWidth1710_1, "f_{0}(1710) Width (2023)", "pe");
    leg1710Width->AddEntry(hWidth1710_3, "f_{0}(1710) Width (2022)", "pe");
    leg1710Width->AddEntry(line1710Width, "PDG value", "l");
    leg1710Width->Draw();
    cWidth1710->SaveAs(outputPath + "/Width1710.png");

    TCanvas *cMass1525 = new TCanvas("cMass1525", "Mass vs #it{p}_{T} for f_{2}(1525)", 720, 720);
    SetCanvasStyle(cMass1525, 0.18, 0.03, 0.05, 0.14);
    hMass1525_1->GetYaxis()->SetRangeUser(1.49, 1.59);
    hMass1525_1->SetMarkerStyle(20);
    hMass1525_1->SetMarkerColor(kRed);
    hMass1525_1->SetLineColor(kRed);
    hMass1525_1->Draw("pe");
    hMass1525_2->SetMarkerStyle(21);
    hMass1525_2->SetMarkerColor(kBlue);
    hMass1525_2->SetLineColor(kBlue);
    hMass1525_2->Draw("pe same");
    hMass1525_3->SetMarkerStyle(22);
    hMass1525_3->SetMarkerColor(kGreen + 2);
    hMass1525_3->SetLineColor(kGreen + 2);
    hMass1525_3->Draw("pe same");
    hMass1525_4->SetMarkerStyle(23);
    hMass1525_4->SetMarkerColor(kMagenta);
    hMass1525_4->SetLineColor(kMagenta);
    hMass1525_4->Draw("pe same");
    TLine *line1525Mass = new TLine(0, f1525Mass, 12, f1525Mass);
    line1525Mass->SetLineStyle(2);
    line1525Mass->SetLineColor(2);
    line1525Mass->Draw();
    leg1710Mass->Draw();
    cMass1525->SaveAs(outputPath + "/Mass1525.png");

    TCanvas *cMass1270 = new TCanvas("cMass1270", "Mass vs #it{p}_{T} for f_{2}(1270)", 720, 720);
    SetCanvasStyle(cMass1270, 0.18, 0.03, 0.05, 0.14);
    hMass1270_1->GetYaxis()->SetRangeUser(1.15, 1.46);
    hMass1270_1->SetMarkerStyle(20);
    hMass1270_1->SetMarkerColor(kRed);
    hMass1270_1->SetLineColor(kRed);
    hMass1270_1->Draw("pe");
    hMass1270_2->SetMarkerStyle(21);
    hMass1270_2->SetMarkerColor(kBlue);
    hMass1270_2->SetLineColor(kBlue);
    hMass1270_2->Draw("pe same");
    hMass1270_3->SetMarkerStyle(22);
    hMass1270_3->SetMarkerColor(kGreen + 2);
    hMass1270_3->SetLineColor(kGreen + 2);
    hMass1270_3->Draw("pe same");
    hMass1270_4->SetMarkerStyle(23);
    hMass1270_4->SetMarkerColor(kMagenta);
    hMass1270_4->SetLineColor(kMagenta);
    hMass1270_4->Draw("pe same");
    TLine *line1270Mass = new TLine(0, f1270Mass, 12, f1270Mass);
    line1270Mass->SetLineStyle(2);
    line1270Mass->SetLineColor(2);
    line1270Mass->Draw();
    leg1710Mass->Draw();
    cMass1270->SaveAs(outputPath + "/Mass1270.png");

    TCanvas *cMass1320 = new TCanvas("cMass1320", "Mass vs #it{p}_{T} for a_{2}(1320)", 720, 720);
    SetCanvasStyle(cMass1320, 0.18, 0.03, 0.05, 0.14);
    hMass1320_1->GetYaxis()->SetRangeUser(1.26, 1.42);
    hMass1320_1->SetMarkerStyle(20);
    hMass1320_1->SetMarkerColor(kRed);
    hMass1320_1->SetLineColor(kRed);
    hMass1320_1->Draw("pe");
    hMass1320_2->SetMarkerStyle(21);
    hMass1320_2->SetMarkerColor(kBlue);
    hMass1320_2->SetLineColor(kBlue);
    hMass1320_2->Draw("pe same");
    hMass1320_3->SetMarkerStyle(22);
    hMass1320_3->SetMarkerColor(kGreen + 2);
    hMass1320_3->SetLineColor(kGreen + 2);
    hMass1320_3->Draw("pe same");
    hMass1320_4->SetMarkerStyle(23);
    hMass1320_4->SetMarkerColor(kMagenta);
    hMass1320_4->SetLineColor(kMagenta);
    hMass1320_4->Draw("pe same");
    TLine *line1320Mass = new TLine(0, a1320Mass, 12, a1320Mass);
    line1320Mass->SetLineStyle(2);
    line1320Mass->SetLineColor(2);
    line1320Mass->Draw();
    leg1710Mass->Draw();
    cMass1320->SaveAs(outputPath + "/Mass1320.png");

    TCanvas *cYieldCorrectedf1525 = new TCanvas("cYieldCorrectedf1525", "Yield vs #it{p}_{T} for f_{2}(1525)", 720, 720);
    SetCanvasStyle(cYieldCorrectedf1525, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    hYield1525Corrected_1->SetMarkerStyle(20);
    hYield1525Corrected_1->SetMaximum(0.01);
    hYield1525Corrected_1->SetMinimum(1e-6);
    hYield1525Corrected_1->SetMarkerColor(kRed);
    hYield1525Corrected_1->SetLineColor(kRed);
    hYield1525Corrected_1->Draw("pe");
    hYield1525Corrected_2->SetMarkerStyle(21);
    hYield1525Corrected_2->SetMarkerColor(kBlue);
    hYield1525Corrected_2->SetLineColor(kBlue);
    hYield1525Corrected_2->Draw("pe same");
    hYield1525Corrected_3->SetMarkerStyle(22);
    hYield1525Corrected_3->SetMarkerColor(kGreen + 2);
    hYield1525Corrected_3->SetLineColor(kGreen + 2);
    hYield1525Corrected_3->Draw("pe same");
    hYield1525Corrected_4->SetMarkerStyle(23);
    hYield1525Corrected_4->SetMarkerColor(kMagenta);
    hYield1525Corrected_4->SetLineColor(kMagenta);
    hYield1525Corrected_4->Draw("pe same");
    // TLegend *legYield = new TLegend(0.55, 0.75, 0.85, 0.93);
    TLegend *legYield = new TLegend(0.5, 0.6, 0.85, 0.93);
    legYield->SetBorderSize(0);
    legYield->SetFillStyle(0);
    legYield->SetTextSize(0.03);
    legYield->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    legYield->AddEntry(hYield1525Corrected_4, "Combined fit (p_{T} > 2 GeV/c)", "pe");
    legYield->AddEntry(hYield1525Corrected_1, "Lowest p_{T} bin: 0-1 GeV/c", "pe");
    legYield->AddEntry(hYield1525Corrected_2, "Lowest p_{T} bin: 1-2 GeV/c", "pe");
    legYield->AddEntry(hYield1525Corrected_3, "Lowest p_{T} bin: 2-3 GeV/c", "pe");
    // legYield->AddEntry(hYield1525Corrected_1, "f_{0}(1710) Width Free (2023)", "pe");
    // legYield->AddEntry(hYield1525Corrected_2, "f_{0}(1710) Width Fixed (2023)", "pe");
    // legYield->AddEntry(hYield1525Corrected_3, "f_{0}(1710) Width Free (2022)", "pe");
    // legYield->AddEntry(hYield1525Corrected_4, "f_{0}(1710) Width Fixed (2022)", "pe");

    // /*************meanpT*****************byresonance*******************package*************************/
    Double_t min = 0.0;
    Double_t max = 10.0;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "REBMS0+";
    TString logfilename = "log.root";
    Double_t minfit = 2.0;
    Double_t maxfit = 10.0;

    TH1F *h1 = (TH1F *)hYield1525Corrected_1->Clone("h1");
    TH1F *h2 = (TH1F *)hYield1525Corrected_2->Clone("h2");
    TH1F *h3 = (TH1F *)hYield1710Corrected_1->Clone("h3");
    TH1F *h4 = (TH1F *)hYield1710Corrected_2->Clone("h4");
    TH1F *h5 = (TH1F *)hYield1525Corrected_3->Clone("h5");
    TH1F *h6 = (TH1F *)hYield1525Corrected_4->Clone("h6");
    TH1F *h7 = (TH1F *)hYield1710Corrected_3->Clone("h7");
    TH1F *h8 = (TH1F *)hYield1710Corrected_4->Clone("h8");

    for (int i = 1; i <= h1->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr = (0.1 * h2->GetBinContent(i));
        h1->SetBinError(i, systemerr);
        h2->SetBinError(i, systemerr);
        h3->SetBinError(i, systemerr);
        h4->SetBinError(i, systemerr);
        h5->SetBinError(i, systemerr);
        h6->SetBinError(i, systemerr);
        h7->SetBinError(i, systemerr);
        h8->SetBinError(i, systemerr);
    }

    TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 10.0, 4);
    fitFcn->SetParameter(0, 5.0);
    fitFcn->SetParameter(1, 0.5);
    fitFcn->FixParameter(2, 1.525);
    fitFcn->SetParameter(3, 0.35);
    fitFcn->SetParNames("n", "dn/dy", "mass", "T");

    TH1 *hout = YieldMean(h1, h1, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // cout << "Yield dN/dy = " << hout->GetBinContent(1) << " +- " << hout->GetBinContent(2) << endl;
    // cout << "Mean pT = " << hout->GetBinContent(5) << " +- " << hout->GetBinContent(6) << endl;
    fitFcn->SetLineColor(2);
    fitFcn->SetLineStyle(2);
    cYieldCorrectedf1525->cd();
    fitFcn->Draw("l same");

    TF1 *fitFcn1 = new TF1("fitfunc", FuncLavy, 0.0, 10.0, 4);
    fitFcn1->SetParameter(0, 5.0);
    fitFcn1->SetParameter(1, 0.5);
    fitFcn1->FixParameter(2, 1.525);
    fitFcn1->SetParameter(3, 0.35);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");

    TH1 *hout2 = YieldMean(h2, h2, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    fitFcn1->SetLineColor(kBlue);
    fitFcn1->SetLineStyle(2);
    cYieldCorrectedf1525->cd();
    fitFcn1->Draw("l same");
    legYield->AddEntry(fitFcn, "Levy-Tsallis Fit", "l");
    legYield->Draw();

    TF1 *fitFcn_temp1 = new TF1("fitfunc_temp1", FuncLavy, 0.0, 10.0, 4);
    fitFcn_temp1->SetParameter(0, 5.0);
    fitFcn_temp1->SetParameter(1, 0.5);
    fitFcn_temp1->FixParameter(2, 1.525);
    fitFcn_temp1->SetParameter(3, 0.35);
    fitFcn_temp1->SetParNames("n", "dn/dy", "mass", "T");
    TH1F *hout_temp1 = (TH1F *)YieldMean(h5, h5, fitFcn_temp1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    fitFcn_temp1->SetLineColor(kGreen + 2);
    fitFcn_temp1->SetLineStyle(2);
    cYieldCorrectedf1525->cd();
    fitFcn_temp1->Draw("l same");

    TF1 *fitFcn_temp2 = new TF1("fitfunc_temp2", FuncLavy, 0.0, 10.0, 4);
    fitFcn_temp2->SetParameter(0, 5.0);
    fitFcn_temp2->SetParameter(1, 0.5);
    fitFcn_temp2->FixParameter(2, 1.525);
    fitFcn_temp2->SetParameter(3, 0.35);
    fitFcn_temp2->SetParNames("n", "dn/dy", "mass", "T");
    TH1F *hout_temp2 = (TH1F *)YieldMean(h6, h6, fitFcn_temp2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    fitFcn_temp2->SetLineColor(kMagenta);
    fitFcn_temp2->SetLineStyle(2);
    cYieldCorrectedf1525->cd();
    fitFcn_temp2->Draw("l same");
    cYieldCorrectedf1525->SaveAs(outputPath + "/CorrectedYieldf2.png");

    TCanvas *cYieldCorrected1710 = new TCanvas("cYieldCorrected1710", "Yield vs #it{p}_{T} for f_{0}(1710)", 720, 720);
    SetCanvasStyle(cYieldCorrected1710, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    hYield1710Corrected_1->SetMarkerStyle(20);
    hYield1710Corrected_1->SetMarkerColor(kRed);
    hYield1710Corrected_1->SetLineColor(kRed);
    hYield1710Corrected_1->SetMarkerStyle(20);
    hYield1710Corrected_1->SetMaximum(0.008);
    hYield1710Corrected_1->SetMinimum(1e-6);
    hYield1710Corrected_1->Draw("pe");
    hYield1710Corrected_2->SetMarkerStyle(21);
    hYield1710Corrected_2->SetMarkerColor(kBlue);
    hYield1710Corrected_2->SetLineColor(kBlue);
    hYield1710Corrected_2->Draw("pe same");
    hYield1710Corrected_3->SetMarkerStyle(22);
    hYield1710Corrected_3->SetMarkerColor(kGreen + 2);
    hYield1710Corrected_3->SetLineColor(kGreen + 2);
    hYield1710Corrected_3->Draw("pe same");
    hYield1710Corrected_4->SetMarkerStyle(23);
    hYield1710Corrected_4->SetMarkerColor(kMagenta);
    hYield1710Corrected_4->SetLineColor(kMagenta);
    hYield1710Corrected_4->Draw("pe same");
    legYield->Draw();

    fitFcn->FixParameter(2, 1.710);

    TH1F *hout3 = (TH1F *)YieldMean(h3, h3, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    fitFcn->SetLineColor(2);
    fitFcn->SetLineStyle(2);
    cYieldCorrected1710->cd();
    fitFcn->Draw("l same");

    TF1 *fitFcn2 = new TF1("fitfunc2", FuncLavy, 0.0, 10.0, 4);
    fitFcn2->SetParameter(0, 5.0);
    // fitFcn2->SetParameter(1, 0.05);
    fitFcn2->SetParameter(1, 0.5);
    fitFcn2->FixParameter(2, 1.710);
    fitFcn2->SetParameter(3, 0.35);
    fitFcn2->SetParNames("n", "dn/dy", "mass", "T");

    TH1F *hout4 = (TH1F *)YieldMean(h4, h4, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    fitFcn2->SetLineColor(kBlue);
    fitFcn2->SetLineStyle(2);
    cYieldCorrected1710->cd();
    fitFcn2->Draw("l same");

    TF1 *fitFcn_temp3 = new TF1("fitfunc_temp3", FuncLavy, 0.0, 10.0, 4);
    fitFcn_temp3->SetParameter(0, 5.0);
    fitFcn_temp3->SetParameter(1, 0.5);
    fitFcn_temp3->FixParameter(2, 1.710);
    fitFcn_temp3->SetParameter(3, 0.35);
    fitFcn_temp3->SetParNames("n", "dn/dy", "mass", "T");
    TH1F *hout_temp3 = (TH1F *)YieldMean(h7, h7, fitFcn_temp3, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    fitFcn_temp3->SetLineColor(kGreen + 2);
    fitFcn_temp3->SetLineStyle(2);
    cYieldCorrected1710->cd();
    fitFcn_temp3->Draw("l same");
    TF1 *fitFcn_temp4 = new TF1("fitfunc_temp4", FuncLavy, 0.0, 10.0, 4);
    fitFcn_temp4->SetParameter(0, 5.0);
    fitFcn_temp4->SetParameter(1, 0.5);
    fitFcn_temp4->FixParameter(2, 1.710);
    fitFcn_temp4->SetParameter(3, 0.35);
    fitFcn_temp4->SetParNames("n", "dn/dy", "mass", "T");
    TH1F *hout_temp4 = (TH1F *)YieldMean(h8, h8, fitFcn_temp4, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    fitFcn_temp4->SetLineColor(kMagenta);
    fitFcn_temp4->SetLineStyle(2);
    cYieldCorrected1710->cd();
    fitFcn_temp4->Draw("l same");
    cYieldCorrected1710->SaveAs(outputPath + "/CorrectedYieldf0.png");

    TCanvas *cYieldRatio = new TCanvas("cYieldRatio", "Yield ratio vs #it{p}_{T} for f_{0}(1710)/f_{2}'(1525)", 720, 720);
    SetCanvasStyle(cYieldRatio, 0.18, 0.03, 0.05, 0.14);
    hYieldRatio_1->SetMarkerStyle(20);
    hYieldRatio_1->SetMaximum(2.5);
    hYieldRatio_1->SetLineColor(kRed);
    hYieldRatio_1->SetMarkerColor(kRed);
    hYieldRatio_1->Draw("pe");
    hYieldRatio_2->SetMarkerStyle(21);
    hYieldRatio_2->SetMarkerColor(kBlue);
    hYieldRatio_2->SetLineColor(kBlue);
    hYieldRatio_2->Draw("pe same");
    hYieldRatio_3->SetMarkerStyle(22);
    hYieldRatio_3->SetMarkerColor(kGreen + 2);
    hYieldRatio_3->SetLineColor(kGreen + 2);
    hYieldRatio_3->Draw("pe same");
    hYieldRatio_4->SetMarkerStyle(23);
    hYieldRatio_4->SetMarkerColor(kMagenta);
    hYieldRatio_4->SetLineColor(kMagenta);
    hYieldRatio_4->Draw("pe same");
    // TLegend *legYieldRatio = new TLegend(0.55, 0.75, 0.85, 0.93);
    TLegend *legYieldRatio = new TLegend(0.55, 0.65, 0.87, 0.93);
    legYieldRatio->SetBorderSize(0);
    legYieldRatio->SetFillStyle(0);
    legYieldRatio->SetTextSize(0.03);
    legYieldRatio->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    legYieldRatio->AddEntry(hYieldRatio_1, "f_{0}(1710) Width Free (2023)", "pe");
    legYieldRatio->AddEntry(hYieldRatio_2, "f_{0}(1710) Width Fixed (2023)", "pe");
    legYieldRatio->AddEntry(hYieldRatio_3, "f_{0}(1710) Width Free (2022)", "pe");
    legYieldRatio->AddEntry(hYieldRatio_4, "f_{0}(1710) Width Fixed (2022)", "pe");
    legYieldRatio->Draw();
    cYieldRatio->SaveAs(outputPath + "/YieldRatiof0f2.png");

    TCanvas *cRawYieldf2 = new TCanvas("cRawYieldf2", "Raw Yield vs #it{p}_{T} for f_{2}(1525)", 720, 720);
    SetCanvasStyle(cRawYieldf2, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    hYield1525Raw_3->SetMarkerStyle(20);
    hYield1525Raw_3->SetMarkerColor(kRed);
    hYield1525Raw_3->SetLineColor(kRed);
    hYield1525Raw_3->Draw("pe");
    hYield1525Raw_4->SetMarkerStyle(21);
    hYield1525Raw_4->SetMarkerColor(kMagenta);
    hYield1525Raw_4->SetLineColor(kMagenta);
    hYield1525Raw_4->Draw("pe same");
    // hYield1525Raw_3->SetMarkerStyle(22);
    // hYield1525Raw_3->SetMarkerColor(kGreen + 2);
    // hYield1525Raw_3->SetLineColor(kGreen + 2);
    // hYield1525Raw_3->Draw("pe same");
    // hYield1525Raw_4->SetMarkerStyle(23);
    // hYield1525Raw_4->SetMarkerColor(kBlue);
    // hYield1525Raw_4->SetLineColor(kBlue);
    // hYield1525Raw_4->Draw("pe same");
    TLegend *legRawYieldf2 = new TLegend(0.63, 0.75, 0.89, 0.93);
    legRawYieldf2->SetBorderSize(0);
    legRawYieldf2->SetFillStyle(0);
    legRawYieldf2->SetTextSize(0.03);
    legRawYieldf2->SetHeader("f_2(1525) Raw p_{T} spectra");
    // legRawYieldf2->AddEntry(hYield1525Raw_1, "Single rBW fit", "pe");
    // legRawYieldf2->AddEntry(hYield1525Raw_4, "4 rBW fit", "pe");
    legRawYieldf2->AddEntry(hYield1525Raw_3, "f_{0} width free", "pe");
    legRawYieldf2->AddEntry(hYield1525Raw_4, "f_{0} width fixed", "pe");
    legRawYieldf2->Draw();
    cRawYieldf2->SaveAs(outputPath + "/RawYieldf2.png");

    TCanvas *cRawYieldRatio = new TCanvas("cRawYieldRatio", "Raw Yield ratio vs #it{p}_{T}", 720, 720);
    SetCanvasStyle(cRawYieldRatio, 0.18, 0.03, 0.05, 0.14);
    TH1F *hRawYieldRatio = (TH1F *)hYield1710Raw_3->Clone("hRawYieldRatio");
    hRawYieldRatio->Divide(hYield1710Raw_4);
    hRawYieldRatio->SetMarkerStyle(20);
    hRawYieldRatio->SetMaximum(1.5);
    hRawYieldRatio->SetMinimum(0.0);
    hRawYieldRatio->SetLineColor(kRed);
    hRawYieldRatio->SetMarkerColor(kRed);
    hRawYieldRatio->GetYaxis()->SetTitle("Raw yield f_{0}(1710) (Width Free / Width Fixed)");
    hRawYieldRatio->Draw("pe");
    cRawYieldRatio->SaveAs(outputPath + "/RawYieldRatioWidthfreeFix.png");

    TCanvas *cRawYieldf0 = new TCanvas("cRawYieldf0", "Raw Yield vs #it{p}_{T} for f_{0}(1710)", 720, 720);
    SetCanvasStyle(cRawYieldf0, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    hYield1710Raw_3->SetMarkerStyle(20);
    hYield1710Raw_3->SetMarkerColor(kRed);
    hYield1710Raw_3->SetLineColor(kRed);
    hYield1710Raw_3->Draw("pe");
    hYield1710Raw_4->SetMarkerStyle(21);
    hYield1710Raw_4->SetMarkerColor(kBlue);
    hYield1710Raw_4->SetLineColor(kBlue);
    hYield1710Raw_4->Draw("pe same");
    // hYield1710Raw_3->SetMarkerStyle(22);
    // hYield1710Raw_3->SetMarkerColor(kGreen + 2);
    // hYield1710Raw_3->SetLineColor(kGreen + 2);
    // hYield1710Raw_3->Draw("pe same");
    // hYield1710Raw_4->SetMarkerStyle(23);
    // hYield1710Raw_4->SetMarkerColor(kMagenta);
    // hYield1710Raw_4->SetLineColor(kMagenta);
    // hYield1710Raw_4->Draw("pe same");
    legRawYieldf2->Draw();
    cRawYieldf0->SaveAs(outputPath + "/RawYieldf0.png");

    TFile *flightFlavourHadrons = new TFile("../spectra/LightFlavourHadronsProduction.root", "read");
    if (flightFlavourHadrons->IsZombie())
    {
        cout << "Error opening light flavour hadrons production file" << endl;
        return;
    }
    // There are graphs of <pT> in the file in table 26 to 34 for different particles (pion, kaon, K0s, K*(892), phi, proton, Lambda, sigma, omega)
    //(π+/π−,K+/K−,KS0​,K∗(892),ϕ(1020),p/pˉ​,Λ/Λˉ,Σ+/Σ−,Ω−/Ωˉ+)

    int totalParticles = 9; // 5 mesons and 4 baryons
    string particlesLatex[9] = {"#pi", "K", "K^{0}_{S}", "K^{*0}", "#phi", "p", "#Lambda", "#Xi^{-}", "#Omega"};
    double particleMass[9] = {0.13957, 0.49367, 0.49761, 0.89166, 1.01946, 0.93827, 1.11568, 1.3217, 1.67245}; // in GeV/c2
    int colors[9] = {kBlack, kBlue, kGreen + 2, kOrange + 7, kViolet + 7, kCyan + 2, kMagenta + 2, kGray + 2, kPink + 2};
    int markers[9] = {21, 22, 23, 33, 34, 43, 45, 29, 39};
    vector<vector<float>> dNdyvalues_13TeV = {
        {4.775, 0.001, 0.243, 1.0},  // 4.775 ± 0.001 ± 0.243
        {6.205, 0.004, 0.303, 1e-1}, // (6.205 ± 0.004 ± 0.303) × 10⁻¹
        {3.192, 0.004, 0.111, 1e-1}, // (3.192 ± 0.004 ± 0.111) × 10⁻¹
        {2.098, 0.016, 0.200, 1e-1}, // (2.098 ± 0.016 ± 0.200) × 10⁻¹
        {2.750, 0.002, 0.188, 1e-1}, // (2.750 ± 0.002 ± 0.188) × 10⁻¹
        {3.734, 0.040, 0.213, 1e-2}, // (3.734 ± 0.040 ± 0.213) × 10⁻²
        {1.807, 0.005, 0.102, 1e-1}, // (1.807 ± 0.005 ± 0.102) × 10⁻¹
        {1.980, 0.012, 0.082, 1e-2}, // (1.980 ± 0.012 ± 0.082) × 10⁻²
        {1.846, 0.046, 0.122, 1e-3}  // (1.846 ± 0.046 ± 0.122) × 10⁻³
    };

    double meanPtAt13TeV[9];
    double meanPtAt13TeV_err[9];
    TGraphErrors *gMeanPt[9];
    TGraphErrors *gMeanPtvsMass_mesons = new TGraphErrors();
    TGraphErrors *gMeanPtvsMass_baryons = new TGraphErrors();
    for (int i = 0; i < totalParticles; i++)
    {
        gMeanPt[i] = (TGraphErrors *)flightFlavourHadrons->Get(Form("Table %d/Graph1D_y1", 26 + i));
        if (gMeanPt[i] == nullptr)
        {
            cout << "Error reading graph for " << particlesLatex[i] << endl;
            return;
        }
        gMeanPt[i]->SetMarkerColor(colors[i]);
        gMeanPt[i]->SetLineColor(colors[i]);
        gMeanPt[i]->SetMarkerStyle(markers[i]);
        gMeanPt[i]->SetMarkerSize(1);

        // Now from graph store the value of mean pT at 13 TeV
        int nPoints = gMeanPt[i]->GetN();
        double *x = gMeanPt[i]->GetX();
        double *y = gMeanPt[i]->GetY();
        for (int j = 0; j < nPoints; j++)
        {
            if (x[j] == 13.0)
            {
                meanPtAt13TeV[i] = y[j];
                meanPtAt13TeV_err[i] = gMeanPt[i]->GetErrorY(j);

                cout << "Mean pT of " << particlesLatex[i] << " at 13 TeV is " << meanPtAt13TeV[i] << " +- " << meanPtAt13TeV_err[i] << endl;
                break;
            }
        }
        if (i < 5)
        {
            int idx = gMeanPtvsMass_mesons->GetN();
            gMeanPtvsMass_mesons->SetPoint(idx, particleMass[i], meanPtAt13TeV[i]);
            gMeanPtvsMass_mesons->SetPointError(idx, 0, meanPtAt13TeV_err[i]);
        }
        else
        {
            int idx = gMeanPtvsMass_baryons->GetN();
            gMeanPtvsMass_baryons->SetPoint(idx, particleMass[i], meanPtAt13TeV[i]);
            gMeanPtvsMass_baryons->SetPointError(idx, 0, meanPtAt13TeV_err[i]);
        }
    }
    SetGrapherrorStyle(gMeanPtvsMass_mesons);
    gMeanPtvsMass_mesons->SetMarkerStyle(22);
    gMeanPtvsMass_mesons->SetMarkerColor(kBlue);
    gMeanPtvsMass_mesons->SetLineColor(kBlue);
    gMeanPtvsMass_mesons->SetMarkerSize(1.5);
    SetGrapherrorStyle(gMeanPtvsMass_baryons);
    gMeanPtvsMass_baryons->SetMarkerStyle(22);
    gMeanPtvsMass_baryons->SetMarkerColor(kRed);
    gMeanPtvsMass_baryons->SetLineColor(kRed);
    gMeanPtvsMass_baryons->SetMarkerSize(1.5);
    TCanvas *cMeanPt = new TCanvas("cMeanPt", "Mean pT vs mass", 720, 720);
    SetCanvasStyle(cMeanPt, 0.14, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    gMeanPtvsMass_mesons->GetXaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    gMeanPtvsMass_mesons->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/#it{c})");
    gMeanPtvsMass_mesons->GetYaxis()->SetTitleOffset(1.3);
    gMeanPtvsMass_mesons->SetMinimum(0.0);
    gMeanPtvsMass_mesons->SetMaximum(3.15);
    gMeanPtvsMass_mesons->GetXaxis()->SetLimits(0, 1.8);
    gMeanPtvsMass_mesons->Draw("AP");
    gMeanPtvsMass_baryons->Draw("P SAME");

    TF1 *pol1_meson = new TF1("pol1_meson", "pol1", 0.1, 1.8);
    pol1_meson->SetLineColor(kYellow + 2);
    pol1_meson->SetLineStyle(2);
    gMeanPtvsMass_mesons->Fit(pol1_meson, "R");

    TF1 *pol1_baryon = new TF1("pol1_baryon", "pol1", 0.1, 1.8);
    pol1_baryon->SetLineColor(kYellow + 2);
    pol1_baryon->SetLineStyle(2);
    gMeanPtvsMass_baryons->Fit(pol1_baryon, "R");

    // Draw the last marker (f2(1525)) and its error bar
    double f2_mass, f2_meanpt_1, f2_meanpt_2, f2_meanpt_3, f2_meanpt_4;
    gMeanPtvsMass_1->GetPoint(0, f2_mass, f2_meanpt_1);
    double f2_meanpt_err_1 = gMeanPtvsMass_1->GetErrorY(0);
    gMeanPtvsMass_2->GetPoint(0, f2_mass, f2_meanpt_2);
    double f2_meanpt_err_2 = gMeanPtvsMass_2->GetErrorY(0);
    gMeanPtvsMass_3->GetPoint(0, f2_mass, f2_meanpt_3);
    double f2_meanpt_err_3 = gMeanPtvsMass_3->GetErrorY(0);
    gMeanPtvsMass_4->GetPoint(0, f2_mass, f2_meanpt_4);
    double f2_meanpt_err_4 = gMeanPtvsMass_4->GetErrorY(0);
    int f2_marker = 25;      // choose a unique marker style for f2(1525)
    int f2_color = kMagenta; // choose a unique color for f2(1525)
    int f2_marker2 = 24;     // hollow circle marker
    int f2_marker3 = 29;     // star marker
    cout << "f2 mean pT single BW " << f2_meanpt_3 << " +/- " << f2_meanpt_err_3 << endl;
    cout << "f2 mean pT 4 BW " << f2_meanpt_4 << " +/- " << f2_meanpt_err_4 << endl;

    TMarker *marker_f2_1 = new TMarker(f2_mass, f2_meanpt_1, f2_marker);
    marker_f2_1->SetMarkerColor(f2_color);
    marker_f2_1->SetMarkerSize(1.5);
    marker_f2_1->Draw("SAME");
    TLine *errBar_f2 = new TLine(f2_mass, f2_meanpt_1 - f2_meanpt_err_1, f2_mass, f2_meanpt_1 + f2_meanpt_err_1);
    errBar_f2->SetLineColor(f2_color);
    errBar_f2->SetLineWidth(2);
    errBar_f2->Draw("SAME");
    double cap_f2 = 0.01;
    TLine *capLow_f2 = new TLine(f2_mass - cap_f2, f2_meanpt_1 - f2_meanpt_err_1, f2_mass + cap_f2, f2_meanpt_1 - f2_meanpt_err_1);
    TLine *capHigh_f2 = new TLine(f2_mass - cap_f2, f2_meanpt_1 + f2_meanpt_err_1, f2_mass + cap_f2, f2_meanpt_1 + f2_meanpt_err_1);
    capLow_f2->SetLineColor(f2_color);
    capHigh_f2->SetLineColor(f2_color);
    capLow_f2->SetLineWidth(2);
    capHigh_f2->SetLineWidth(2);
    capLow_f2->Draw("SAME");
    capHigh_f2->Draw("SAME");

    TMarker *marker_f2_2 = new TMarker(f2_mass, f2_meanpt_2, f2_marker2);
    marker_f2_2->SetMarkerColor(kRed + 2);
    marker_f2_2->SetMarkerSize(1.5);
    marker_f2_2->Draw("SAME");
    TLine *errBar_f2_2 = new TLine(f2_mass, f2_meanpt_2 - f2_meanpt_err_2, f2_mass, f2_meanpt_2 + f2_meanpt_err_2);
    errBar_f2_2->SetLineColor(kRed + 2);
    errBar_f2_2->SetLineWidth(2);
    errBar_f2_2->Draw("SAME");
    double cap_f2_2 = 0.01;
    TLine *capLow_f2_2 = new TLine(f2_mass - cap_f2_2, f2_meanpt_2 - f2_meanpt_err_2, f2_mass + cap_f2_2, f2_meanpt_2 - f2_meanpt_err_2);
    TLine *capHigh_f2_2 = new TLine(f2_mass - cap_f2_2, f2_meanpt_2 + f2_meanpt_err_2, f2_mass + cap_f2_2, f2_meanpt_2 + f2_meanpt_err_2);
    capLow_f2_2->SetLineColor(kRed + 2);
    capHigh_f2_2->SetLineColor(kRed + 2);
    capLow_f2_2->SetLineWidth(2);
    capHigh_f2_2->SetLineWidth(2);
    capLow_f2_2->Draw("SAME");
    capHigh_f2_2->Draw("SAME");

    TMarker *marker_f2_3 = new TMarker(f2_mass, f2_meanpt_3, f2_marker3);
    marker_f2_3->SetMarkerColor(kGreen + 2);
    marker_f2_3->SetMarkerSize(1.7);
    marker_f2_3->Draw("SAME");
    TLine *errBar_f2_3 = new TLine(f2_mass, f2_meanpt_3 - f2_meanpt_err_3, f2_mass, f2_meanpt_3 + f2_meanpt_err_3);
    errBar_f2_3->SetLineColor(kGreen + 2);
    errBar_f2_3->SetLineWidth(2);
    errBar_f2_3->Draw("SAME");
    double cap_f2_3 = 0.01;
    TLine *capLow_f2_3 = new TLine(f2_mass - cap_f2_3, f2_meanpt_3 - f2_meanpt_err_3, f2_mass + cap_f2_3, f2_meanpt_3 - f2_meanpt_err_3);
    TLine *capHigh_f2_3 = new TLine(f2_mass - cap_f2_3, f2_meanpt_3 + f2_meanpt_err_3, f2_mass + cap_f2_3, f2_meanpt_3 + f2_meanpt_err_3);
    capLow_f2_3->SetLineColor(kGreen + 2);
    capHigh_f2_3->SetLineColor(kGreen + 2);
    capLow_f2_3->SetLineWidth(2);
    capHigh_f2_3->SetLineWidth(2);
    capLow_f2_3->Draw("SAME");
    capHigh_f2_3->Draw("SAME");

    TMarker *marker_f2_4 = new TMarker(f2_mass, f2_meanpt_4, f2_marker2);
    marker_f2_4->SetMarkerColor(kRed + 2);
    marker_f2_4->SetMarkerSize(1.5);
    marker_f2_4->Draw("SAME");
    TLine *errBar_f2_4 = new TLine(f2_mass, f2_meanpt_4 - f2_meanpt_err_4, f2_mass, f2_meanpt_4 + f2_meanpt_err_4);
    errBar_f2_4->SetLineColor(kRed + 2);
    errBar_f2_4->SetLineWidth(2);
    errBar_f2_4->Draw("SAME");
    double cap_f2_4 = 0.01;
    TLine *capLow_f2_4 = new TLine(f2_mass - cap_f2_4, f2_meanpt_4 - f2_meanpt_err_4, f2_mass + cap_f2_4, f2_meanpt_4 - f2_meanpt_err_4);
    TLine *capHigh_f2_4 = new TLine(f2_mass - cap_f2_4, f2_meanpt_4 + f2_meanpt_err_4, f2_mass + cap_f2_4, f2_meanpt_4 + f2_meanpt_err_4);
    capLow_f2_4->SetLineColor(kRed + 2);
    capHigh_f2_4->SetLineColor(kRed + 2);
    capLow_f2_4->SetLineWidth(2);
    capHigh_f2_4->SetLineWidth(2);
    capLow_f2_4->Draw("SAME");
    capHigh_f2_4->Draw("SAME");

    // Draw the last marker (f0(1710)) and its error bar
    double f0_mass, f0_meanpt_1, f0_meanpt_2;
    gMeanPtvsMass_1->GetPoint(1, f0_mass, f0_meanpt_1);
    double f0_meanpt_err_1 = gMeanPtvsMass_1->GetErrorY(1);
    gMeanPtvsMass_2->GetPoint(1, f0_mass, f0_meanpt_2);
    double f0_meanpt_err_2 = gMeanPtvsMass_2->GetErrorY(1);
    int f0_marker = 21;       // choose a unique marker style for f0(1710)
    int f0_color = kBlue + 2; // choose a unique color for f0(1710)
    TMarker *marker_f0_1 = new TMarker(f0_mass, f0_meanpt_1, f0_marker);
    marker_f0_1->SetMarkerColor(f0_color);
    marker_f0_1->SetMarkerSize(1.5);
    marker_f0_1->Draw("SAME");
    TLine *errBar_f0 = new TLine(f0_mass, f0_meanpt_1 - f0_meanpt_err_1, f0_mass, f0_meanpt_1 + f0_meanpt_err_1);
    errBar_f0->SetLineColor(f0_color);
    errBar_f0->SetLineWidth(2);
    errBar_f0->Draw("SAME");
    double cap_f0 = 0.01;
    TLine *capLow_f0 = new TLine(f0_mass - cap_f0, f0_meanpt_1 - f0_meanpt_err_1, f0_mass + cap_f0, f0_meanpt_1 - f0_meanpt_err_1);
    TLine *capHigh_f0 = new TLine(f0_mass - cap_f0, f0_meanpt_1 + f0_meanpt_err_1, f0_mass + cap_f0, f0_meanpt_1 + f0_meanpt_err_1);
    capLow_f0->SetLineColor(f0_color);
    capHigh_f0->SetLineColor(f0_color);
    capLow_f0->SetLineWidth(2);
    capHigh_f0->SetLineWidth(2);
    capLow_f0->Draw("SAME");
    capHigh_f0->Draw("SAME");

    TMarker *marker_f0_2 = new TMarker(f0_mass, f0_meanpt_2, f0_marker);
    marker_f0_2->SetMarkerColor(kGreen + 1);
    marker_f0_2->SetMarkerSize(1.5);
    marker_f0_2->Draw("SAME");
    TLine *errBar_f0_2 = new TLine(f0_mass, f0_meanpt_2 - f0_meanpt_err_2, f0_mass, f0_meanpt_2 + f0_meanpt_err_2);
    errBar_f0_2->SetLineColor(kGreen + 1);
    errBar_f0_2->SetLineWidth(2);
    errBar_f0_2->Draw("SAME");
    double cap_f0_2 = 0.01;
    TLine *capLow_f0_2 = new TLine(f0_mass - cap_f0_2, f0_meanpt_2 - f0_meanpt_err_2, f0_mass + cap_f0_2, f0_meanpt_2 - f0_meanpt_err_2);
    TLine *capHigh_f0_2 = new TLine(f0_mass - cap_f0_2, f0_meanpt_2 + f0_meanpt_err_2, f0_mass + cap_f0_2, f0_meanpt_2 + f0_meanpt_err_2);
    capLow_f0_2->SetLineColor(kGreen + 1);
    capHigh_f0_2->SetLineColor(kGreen + 1);
    capLow_f0_2->SetLineWidth(2);
    capHigh_f0_2->SetLineWidth(2);
    capLow_f0_2->Draw("SAME");
    capHigh_f0_2->Draw("SAME");

    // Add particle names below each point
    TLatex latex;
    latex.SetTextAlign(22);
    latex.SetTextSize(0.035);
    for (int i = 0; i < totalParticles; i++)
    {
        double x = particleMass[i];
        double y = meanPtAt13TeV[i];
        // latex.SetTextColor(colors[i]);
        latex.SetTextColor(kBlack);
        if (i == 2)
            latex.DrawLatex(x + 0.08, y + 0.22, particlesLatex[i].c_str());
        else if (i <= 4)
            latex.DrawLatex(x, y + 0.22, particlesLatex[i].c_str());
        else
            latex.DrawLatex(x, y - 0.20, particlesLatex[i].c_str());
    }
    // latex.SetTextColor(kRed); // match f2(1525) marker color
    // latex.DrawLatex(1.5173, f2_meanpt_1 - 0.22, "f'_{2}(1525)");
    latex.DrawLatex(1.5173, f2_meanpt_1 - 0.5, "f'_{2}(1525)");
    // latex.SetTextColor(kBlue); // match f0(1710) marker color
    latex.DrawLatex(1.68, f0_meanpt_1 + 0.33, "f_{0}(1710)");

    // TLegend *legend4 = new TLegend(0.18, 0.75, 0.85, 0.92);
    // legend4->SetBorderSize(0);
    // legend4->SetFillStyle(0);
    // legend4->SetTextSize(0.033);
    // legend4->SetNColumns(2);
    // legend4->AddEntry(gMeanPtvsMass_mesons, "Light flavour hadrons (13 TeV)", "p");
    // legend4->AddEntry((TObject *)0, "", "");
    // legend4->AddEntry((TObject *)0, "f_{0}(1710) Width free", "");
    // legend4->AddEntry((TObject *)0, "f_{0}(1710) Width fixed", "");
    // legend4->AddEntry(marker_f2_1, "f'_{2}(1525)", "p");
    // legend4->AddEntry(marker_f2_2, "f'_{2}(1525)", "p");
    // legend4->AddEntry(marker_f0_1, "f_{0}(1710)", "p");
    // legend4->AddEntry(marker_f0_2, "f_{0}(1710)", "p");
    // legend4->Draw();

    TLegend *legendTemp1 = new TLegend(0.18, 0.75, 0.85, 0.92);
    legendTemp1->SetBorderSize(0);
    legendTemp1->SetFillStyle(0);
    legendTemp1->SetTextSize(0.033);
    legendTemp1->AddEntry(gMeanPtvsMass_mesons, "Light flavour hadrons (13 TeV)", "p");
    legendTemp1->AddEntry(marker_f2_1, "Lowest p_{T} bin: 0-1 GeV/c", "p");
    legendTemp1->AddEntry(marker_f2_2, "Lowest p_{T} bin: 1-2 GeV/c", "p");
    legendTemp1->AddEntry(marker_f2_3, "Lowest p_{T} bin: 2-3 GeV/c", "p");
    legendTemp1->Draw();

    cMeanPt->SaveAs(outputPath + "/MeanPt_vs_Mass.png");
}