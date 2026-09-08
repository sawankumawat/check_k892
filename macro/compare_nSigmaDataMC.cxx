#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

void compare_nSigmaDataMC()
{
    gStyle->SetOptStat(0);
    TString dataPath = "/home/sawan/Storage/check_k892/data/kstar/LHC22o_pass7";
    TString mcPath = "/home/sawan/Storage/check_k892/mc/LHC24f3c";

    TFile *fData = new TFile(dataPath + "/751768.root", "read");
    // TFile *fMC = new TFile(mcPath + "/750013.root", "read");
    TFile *fMC = new TFile(mcPath + "/734009.root", "read");

    // TFile *fData = new TFile(dataPath + "/679906.root", "read");
    // TFile *fMC = new TFile(mcPath + "/679945.root", "read");

    string multpath = "kstarqa";

    TH3F *hNsigmaTPCKaon = (TH3F *)fData->Get(Form("%s/hPID/Before/hTPCnsigKa_mult_pt", multpath.c_str()));
    TH3F *hNsigmaTPCPion = (TH3F *)fData->Get(Form("%s/hPID/Before/hTPCnsigPi_mult_pt", multpath.c_str()));
    TH3F *hNsigmaTOFKaon = (TH3F *)fData->Get(Form("%s/hPID/Before/hTOFnsigKa_mult_pt", multpath.c_str()));
    TH3F *hNsigmaTOFPion = (TH3F *)fData->Get(Form("%s/hPID/Before/hTOFnsigPi_mult_pt", multpath.c_str()));
    if (hNsigmaTPCKaon == nullptr || hNsigmaTPCPion == nullptr || hNsigmaTOFKaon == nullptr || hNsigmaTOFPion == nullptr)
    {
        cerr << "PID histograms before selection not found!!!!!!!!!!!!" << endl;
        return;
    }

    TH3F *hNsigmaTPCKaonMC = (TH3F *)fMC->Get(Form("%s/hPID/Before/hTPCnsigKa_mult_pt", multpath.c_str()));
    TH3F *hNsigmaTPCPionMC = (TH3F *)fMC->Get(Form("%s/hPID/Before/hTPCnsigPi_mult_pt", multpath.c_str()));
    TH3F *hNsigmaTOFKaonMC = (TH3F *)fMC->Get(Form("%s/hPID/Before/hTOFnsigKa_mult_pt", multpath.c_str()));
    TH3F *hNsigmaTOFPionMC = (TH3F *)fMC->Get(Form("%s/hPID/Before/hTOFnsigPi_mult_pt", multpath.c_str()));
    if (hNsigmaTPCKaonMC == nullptr || hNsigmaTPCPionMC == nullptr || hNsigmaTOFKaonMC == nullptr || hNsigmaTOFPionMC == nullptr)
    {
        cerr << "PID histograms before selection not found!!!!!!!!!!!!" << endl;
        return;
    }

    // double lowpT = 0.8;
    // double highpT = 1.0;

    double lowpT = 0.0;
    double highpT = 100.0;

    double multLow = 70.0;
    double multHigh = 100.0;

    int lowBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multLow + 0.001);
    int highBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multHigh - 0.001);
    int lowBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(lowpT + 0.001);
    int highBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(highpT - 0.001);

    TH1F *h1DNsigmaTPCKaon_pt = (TH1F *)hNsigmaTPCKaon->ProjectionX(Form("h1D_TPC_Ka_pt_%d", 0), lowBinMult, highBinMult, lowBinpT, highBinpT);
    TH1F *h1DNsigmaTPCPion_pt = (TH1F *)hNsigmaTPCPion->ProjectionX(Form("h1D_TPC_Pi_pt_%d", 0), lowBinMult, highBinMult, lowBinpT, highBinpT);
    TH1F *h1DNsigmaTOFKaon_pt = (TH1F *)hNsigmaTOFKaon->ProjectionX(Form("h1D_TOF_Ka_pt_%d", 0), lowBinMult, highBinMult, lowBinpT, highBinpT);
    TH1F *h1DNsigmaTOFPion_pt = (TH1F *)hNsigmaTOFPion->ProjectionX(Form("h1D_TOF_Pi_pt_%d", 0), lowBinMult, highBinMult, lowBinpT, highBinpT);

    TH1F *h1DNsigmaTPCKaonMC_pt = (TH1F *)hNsigmaTPCKaonMC->ProjectionX(Form("h1D_TPC_Ka_MC_pt_%d", 0), lowBinMult, highBinMult, lowBinpT, highBinpT);
    TH1F *h1DNsigmaTPCPionMC_pt = (TH1F *)hNsigmaTPCPionMC->ProjectionX(Form("h1D_TPC_Pi_MC_pt_%d", 0), lowBinMult, highBinMult, lowBinpT, highBinpT);
    TH1F *h1DNsigmaTOFKaonMC_pt = (TH1F *)hNsigmaTOFKaonMC->ProjectionX(Form("h1D_TOF_Ka_MC_pt_%d", 0), lowBinMult, highBinMult, lowBinpT, highBinpT);
    TH1F *h1DNsigmaTOFPionMC_pt = (TH1F *)hNsigmaTOFPionMC->ProjectionX(Form("h1D_TOF_Pi_MC_pt_%d", 0), lowBinMult, highBinMult, lowBinpT, highBinpT);

    int integralBinLow = h1DNsigmaTPCKaon_pt->FindBin(-1.0 + 0.001);
    int integralBinHigh = h1DNsigmaTPCKaon_pt->FindBin(1.0 - 0.001);

    TCanvas *cNsigmaTPCPion = new TCanvas("cNsigmaTPCPion", "TPC nSigma Pion", 720, 720);
    SetCanvasStyle(cNsigmaTPCPion, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(h1DNsigmaTPCPion_pt);
    h1DNsigmaTPCPion_pt->GetXaxis()->SetTitle("#it{n#sigma}_{TPC}^{#pi}");
    h1DNsigmaTPCPion_pt->GetYaxis()->SetTitle("Counts");
    h1DNsigmaTPCPion_pt->GetXaxis()->SetRangeUser(-2, 2);
    h1DNsigmaTPCPion_pt->SetMaximum(h1DNsigmaTPCPion_pt->GetMaximum() * 1.2);
    h1DNsigmaTPCPion_pt->Draw("pe");
    SetHistoQA(h1DNsigmaTPCPionMC_pt);
    h1DNsigmaTPCPionMC_pt->SetMarkerColor(kRed);
    h1DNsigmaTPCPionMC_pt->SetLineColor(kRed);
    h1DNsigmaTPCPionMC_pt->Scale(h1DNsigmaTPCPion_pt->Integral(integralBinLow, integralBinHigh) / h1DNsigmaTPCPionMC_pt->Integral(integralBinLow, integralBinHigh));
    h1DNsigmaTPCPionMC_pt->Draw("pe same");

    TLegend *legend = new TLegend(0.75, 0.75, 0.9, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(h1DNsigmaTPCPion_pt, "Data", "p");
    legend->AddEntry(h1DNsigmaTPCPionMC_pt, "MC", "p");
    legend->Draw();
    // cNsigmaTPCPion->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/PID_Plots/nSigmaTPCPi_DataMC.png");

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.03);
    latex->DrawLatex(0.2, 0.9, "ALICE");
    latex->DrawLatex(0.2, 0.85, "pp, #sqrt{#it{s}} = 13.6 TeV");
    latex->DrawLatex(0.2, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", lowpT, highpT));
    latex->DrawLatex(0.2, 0.75, "FT0M: 0-100%");

    TCanvas *cNsigmaTPCKaon = new TCanvas("cNsigmaTPCKaon", "TPC nSigma Kaon", 720, 720);
    SetCanvasStyle(cNsigmaTPCKaon, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(h1DNsigmaTPCKaon_pt);
    h1DNsigmaTPCKaon_pt->GetXaxis()->SetTitle("#it{n#sigma}_{TPC}^{K}");
    h1DNsigmaTPCKaon_pt->GetYaxis()->SetTitle("Counts");
    h1DNsigmaTPCKaon_pt->GetXaxis()->SetRangeUser(-2, 2);
    h1DNsigmaTPCKaon_pt->SetMaximum(h1DNsigmaTPCKaon_pt->GetMaximum() * 1.2);
    h1DNsigmaTPCKaon_pt->Draw("pe");
    SetHistoQA(h1DNsigmaTPCKaonMC_pt);
    h1DNsigmaTPCKaonMC_pt->SetMarkerColor(kRed);
    h1DNsigmaTPCKaonMC_pt->SetLineColor(kRed);
    h1DNsigmaTPCKaonMC_pt->Scale(h1DNsigmaTPCKaon_pt->Integral(integralBinLow, integralBinHigh) / h1DNsigmaTPCKaonMC_pt->Integral(integralBinLow, integralBinHigh));
    h1DNsigmaTPCKaonMC_pt->Draw("pe same");

    legend->Draw();
    latex->DrawLatex(0.17, 0.9, "ALICE");
    latex->DrawLatex(0.17, 0.85, "pp, #sqrt{#it{s}} = 13.6 TeV");
    latex->DrawLatex(0.17, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", lowpT, highpT));
    latex->DrawLatex(0.17, 0.75, "FT0M: 0-100%");
    // cNsigmaTPCKaon->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/PID_Plots/nSigmaTPCKa_DataMC.png");

    TCanvas *cNsigmaTOFPion = new TCanvas("cNsigmaTOFPion", "TOF nSigma Pion", 720, 720);
    SetCanvasStyle(cNsigmaTOFPion, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(h1DNsigmaTOFPion_pt);
    h1DNsigmaTOFPion_pt->GetXaxis()->SetTitle("#it{n#sigma}_{TOF}^{#pi}");
    h1DNsigmaTOFPion_pt->GetYaxis()->SetTitle("Counts");
    h1DNsigmaTOFPion_pt->SetMaximum(h1DNsigmaTOFPion_pt->GetMaximum() * 1.2);
    h1DNsigmaTOFPion_pt->GetXaxis()->SetRangeUser(-2, 2);
    h1DNsigmaTOFPion_pt->Draw("pe");
    SetHistoQA(h1DNsigmaTOFPionMC_pt);
    h1DNsigmaTOFPionMC_pt->SetMarkerColor(kRed);
    h1DNsigmaTOFPionMC_pt->SetLineColor(kRed);
    h1DNsigmaTOFPionMC_pt->Scale(h1DNsigmaTOFPion_pt->Integral(integralBinLow, integralBinHigh) / h1DNsigmaTOFPionMC_pt->Integral(integralBinLow, integralBinHigh));
    h1DNsigmaTOFPionMC_pt->Draw("pe same");

    legend->Draw();
    latex->DrawLatex(0.2, 0.9, "ALICE");
    latex->DrawLatex(0.2, 0.85, "pp, #sqrt{#it{s}} = 13.6 TeV");
    latex->DrawLatex(0.2, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", lowpT, highpT));
    latex->DrawLatex(0.2, 0.75, "FT0M: 0-100%");
    // cNsigmaTOFPion->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/PID_Plots/nSigmaTOFPi_DataMC.png");

    TCanvas *cNsigmaTOFKaon = new TCanvas("cNsigmaTOFKaon", "TOF nSigma Kaon", 720, 720);
    SetCanvasStyle(cNsigmaTOFKaon, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(h1DNsigmaTOFKaon_pt);
    h1DNsigmaTOFKaon_pt->GetXaxis()->SetTitle("#it{n#sigma}_{TOF}^{K}");
    h1DNsigmaTOFKaon_pt->GetYaxis()->SetTitle("Counts");
    h1DNsigmaTOFKaon_pt->GetXaxis()->SetRangeUser(-2, 2);
    h1DNsigmaTOFKaon_pt->SetMaximum(h1DNsigmaTOFKaon_pt->GetMaximum() * 1.2);
    h1DNsigmaTOFKaon_pt->Draw("pe");
    SetHistoQA(h1DNsigmaTOFKaonMC_pt);
    h1DNsigmaTOFKaonMC_pt->SetMarkerColor(kRed);
    h1DNsigmaTOFKaonMC_pt->SetLineColor(kRed);
    h1DNsigmaTOFKaonMC_pt->Scale(h1DNsigmaTOFKaon_pt->Integral(integralBinLow, integralBinHigh) / h1DNsigmaTOFKaonMC_pt->Integral(integralBinLow, integralBinHigh));
    h1DNsigmaTOFKaonMC_pt->Draw("pe same");

    legend->Draw();
    latex->DrawLatex(0.17, 0.9, "ALICE");
    latex->DrawLatex(0.17, 0.85, "pp, #sqrt{#it{s}} = 13.6 TeV");
    latex->DrawLatex(0.17, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", lowpT, highpT));
    latex->DrawLatex(0.17, 0.75, "FT0M: 0-100%");
    // cNsigmaTOFKaon->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/PID_Plots/nSigmaTOFKa_DataMC.png");
}