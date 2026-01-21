#include <iostream>
using namespace std;
#include "src/style.h"

void YieldRatioThermalFist()
{
    gStyle->SetOptStat(0);
    TString path = "/home/sawan/Thermal-FIST/src/examples/cpc";
    TFile *f1710 = new TFile(path + "/Thermalf0.root", "read");
    TFile *f1525 = new TFile(path + "/Thermalf2.root", "read");
    TFile *fK0s = new TFile(path + "/ThermalK0s.root", "read");
    TFile *fPion = new TFile(path + "/ThermalPionDecay.root", "read");
    if (f1710->IsZombie() || f1525->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TTree *t1710 = (TTree *)f1710->Get("f01710");
    TTree *t1525 = (TTree *)f1525->Get("f21525");
    TTree *tK0s = (TTree *)fK0s->Get("k0s");
    TTree *tPion = (TTree *)fPion->Get("pionDecayed");

    double pt1710, pt1525, ptK0s, ptPion;
    t1710->SetBranchAddress("pT", &pt1710);
    t1525->SetBranchAddress("pT", &pt1525);
    tK0s->SetBranchAddress("pT", &ptK0s);
    tPion->SetBranchAddress("pT", &ptPion);
    int totalEvents1710 = t1710->GetEntries();
    int totalEvents1525 = t1525->GetEntries();
    int totalEventsK0s = tK0s->GetEntries();
    int totalEventsPion = tPion->GetEntries();
    cout << "Total f1710 events: " << totalEvents1710 << endl;
    cout << "Total f1525 events: " << totalEvents1525 << endl;
    cout << "Total K0s events: " << totalEventsK0s << endl;
    cout << "Total Pion events: " << totalEventsPion << endl;

    float ptBins[] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0}; // 2022 dataset
    int nBins = sizeof(ptBins) / sizeof(ptBins[0]) - 1;
    TH1F *hPt1710 = new TH1F("hPt1710", "f0(1710) pT distribution", nBins, ptBins);
    TH1F *hPt1525 = new TH1F("hPt1525", "f2'(1525) pT distribution", nBins, ptBins);
    TH1F *hPtK0s = new TH1F("hPtK0s", "K0s pT distribution", nBins, ptBins);
    TH1F *hPtPion = new TH1F("hPtPion", "Pion pT distribution", nBins, ptBins);

    for (int i = 0; i < totalEvents1710; i++)
    {
        t1710->GetEntry(i);
        hPt1710->Fill(pt1710);
    }
    for (int i = 0; i < totalEvents1525; i++)
    {
        t1525->GetEntry(i);
        hPt1525->Fill(pt1525);
    }
    for (int i = 0; i < totalEventsK0s; i++)
    {
        tK0s->GetEntry(i);
        hPtK0s->Fill(ptK0s);
    }
    for (int i = 0; i < totalEventsPion; i++)
    {
        tPion->GetEntry(i);
        hPtPion->Fill(ptPion);
    }
    hPt1710->Scale(1.0 / totalEvents1710);
    hPt1525->Scale(1.0 / totalEvents1525);
    hPtK0s->Scale(1.0 / totalEventsK0s);
    hPtPion->Scale(1.0 / totalEventsPion);
    TH1F *hRatio = (TH1F *)hPt1710->Clone("hRatio");
    hRatio->Divide(hPt1525);
    TH1F *hRatiof0_k0s = (TH1F *)hPt1710->Clone("hRatiof0_k0s");
    hRatiof0_k0s->Divide(hPtK0s);
    TH1F *hRatiof2pion = (TH1F *)hPt1525->Clone("hRatiof2pion");
    hRatiof2pion->Divide(hPtPion);
    TH1F *hRatiof0_pion = (TH1F *)hPt1710->Clone("hRatiof0_pion");
    hRatiof0_pion->Divide(hPtPion);
    TH1F *hRatiof2_k0s = (TH1F *)hPt1525->Clone("hRatiof2_k0s");
    hRatiof2_k0s->Divide(hPtK0s);

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.06);
    TCanvas *cpT = new TCanvas("cpT", "pT distributions", 1080, 720);
    SetCanvasStyle(cpT, 0.18, 0.03, 0.05, 0.14);
    cpT->Divide(2, 2);
    cpT->cd(1);
    gPad->SetLogy();
    SetHistoQA(hPt1710);
    hPt1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPt1710->GetYaxis()->SetTitle("1/N_{ev} * dN/(d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hPt1710->GetYaxis()->SetTitleOffset(1.6);
    hPt1710->SetMaximum(hPt1710->GetMaximum() * 1.5);
    hPt1710->SetMarkerStyle(20);
    hPt1710->Draw("pe");
    lat.DrawLatex(0.2, 0.85, "f_{0}(1710)");
    cpT->cd(2);
    gPad->SetLogy();
    SetHistoQA(hPt1525);
    hPt1525->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPt1525->GetYaxis()->SetTitle("1/N_{ev} * dN/(d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hPt1525->GetYaxis()->SetTitleOffset(1.6);
    hPt1525->SetMaximum(hPt1525->GetMaximum() * 1.5);
    hPt1525->SetMarkerStyle(21);
    hPt1525->Draw("pe");
    lat.DrawLatex(0.2, 0.85, "f_{2}'(1525)");
    cpT->cd(3);
    gPad->SetLogy();
    SetHistoQA(hPtK0s);
    hPtK0s->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPtK0s->GetYaxis()->SetTitle("1/N_{ev} * dN/(d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hPtK0s->GetYaxis()->SetTitleOffset(1.6);
    hPtK0s->SetMaximum(hPtK0s->GetMaximum() * 1.5);
    hPtK0s->SetMarkerStyle(22);
    hPtK0s->Draw("pe");
    lat.DrawLatex(0.2, 0.85, "K^{0}_{S}");
    cpT->cd(4);
    gPad->SetLogy();
    SetHistoQA(hPtPion);
    hPtPion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPtPion->GetYaxis()->SetTitle("1/N_{ev} * dN/(d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hPtPion->GetYaxis()->SetTitleOffset(1.6);
    hPtPion->SetMaximum(hPtPion->GetMaximum() * 1.5);
    hPtPion->SetMarkerStyle(23);
    hPtPion->Draw("pe");
    lat.DrawLatex(0.2, 0.85, "Pion (decayed)");

    TCanvas *cRatio = new TCanvas("cRatio", "Yield ratio f0(1710)/f2'(1525)", 720, 720);
    SetCanvasStyle(cRatio, 0.18, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    SetHistoQA(hRatio);
    hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio->GetYaxis()->SetTitle("f_{0}(1710)/f_{2}'(1525)");
    hRatio->GetYaxis()->SetTitleOffset(1.6);
    hRatio->SetMaximum(hRatio->GetMaximum() * 1.5);
    hRatio->SetMarkerStyle(22);
    hRatio->Draw("pe");
    cRatio->SaveAs(path + "/YieldRatio_f0_f2p.png");

    TCanvas *cRatio_k0s = new TCanvas("cRatio_k0s", "Yield ratio f0(1710)/K0s", 720, 720);
    SetCanvasStyle(cRatio_k0s, 0.18, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    SetHistoQA(hRatiof0_k0s);
    hRatiof0_k0s->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatiof0_k0s->GetYaxis()->SetTitle("f_{0}(1710)/K^{0}_{S}");
    hRatiof0_k0s->GetYaxis()->SetTitleOffset(1.6);
    hRatiof0_k0s->SetMaximum(hRatiof0_k0s->GetMaximum() * 1.5);
    hRatiof0_k0s->SetMarkerStyle(22);
    hRatiof0_k0s->Draw("pe");
    SetHistoQA(hRatiof2_k0s);
    hRatiof2_k0s->SetMarkerStyle(23);
    hRatiof2_k0s->SetLineColor(kRed);
    hRatiof2_k0s->SetMarkerColor(kRed);
    hRatiof2_k0s->Draw("pe same");
    TLegend *leg_k0s = new TLegend(0.7, 0.7, 0.9, 0.85);
    leg_k0s->AddEntry(hRatiof0_k0s, "f_{0}(1710)/K^{0}_{S}", "pe");
    leg_k0s->AddEntry(hRatiof2_k0s, "f_{2}'(1525)/K^{0}_{S}", "pe");
    leg_k0s->SetBorderSize(0);
    leg_k0s->SetTextSize(0.035);
    leg_k0s->SetFillStyle(0);
    leg_k0s->Draw();

    TCanvas *cRatio_pion = new TCanvas("cRatio_pion", "Yield ratio f0(1710)/Pion", 720, 720);
    SetCanvasStyle(cRatio_pion, 0.18, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    SetHistoQA(hRatiof0_pion);
    hRatiof0_pion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatiof0_pion->GetYaxis()->SetTitle("f_{0}(1710)/Pion");
    hRatiof0_pion->GetYaxis()->SetTitleOffset(1.6);
    hRatiof0_pion->SetMaximum(hRatiof0_pion->GetMaximum() * 1.5);
    hRatiof0_pion->SetMarkerStyle(22);
    hRatiof0_pion->Draw("pe");
    SetHistoQA(hRatiof2pion);
    hRatiof2pion->SetMarkerStyle(23);
    hRatiof2pion->SetLineColor(kRed);
    hRatiof2pion->SetMarkerColor(kRed);
    hRatiof2pion->Draw("pe same");
    leg_k0s->Draw();
}