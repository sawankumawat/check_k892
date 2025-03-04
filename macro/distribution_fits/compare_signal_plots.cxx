#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/common_glue.h"
#include "../src/style.h"
using namespace std;

void compare_signal_plots() {
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/351471/KsKs_Channel/";
    gStyle->SetOptStat(0);
    int colors[] = {2, 4, 6, 28, 46};
    int markers[] = {20, 21, 22, 23, 24};
    string path1 = path + "higher-mass-resonances_2sigmaKs/";
    string path2 = path + "higher-mass-resonances_3sigmaKs/";
    string path3 = path + "higher-mass-resonances/";
    string path4 = path + "higher-mass-resonances_5SigmaKs/";

    TFile *file1 = new TFile((path1 + "hglue_ROTATED_norm_2.20_2.30_pt_1.00_10.00.root").c_str(), "READ");
    TFile *file2 = new TFile((path2 + "hglue_ROTATED_norm_2.20_2.30_pt_1.00_10.00.root").c_str(), "READ");
    TFile *file3 = new TFile((path3 + "hglue_ROTATED_norm_2.20_2.30_pt_1.00_10.00.root").c_str(), "READ");
    TFile *file4 = new TFile((path4 + "hglue_ROTATED_norm_2.20_2.30_pt_1.00_10.00.root").c_str(), "READ");
    // TFile *file5 = new TFile((path + "hglue_ROTATED_norm_2.70_2.80.root").c_str(), "READ");

    if(file1->IsZombie() || file2->IsZombie() || file3->IsZombie() || file4->IsZombie()) {
        cout << "Error opening file" << endl;
        return;
    }

    TH1F *h1 = (TH1F *)file1->Get("ksks_subtracted_invmass_pt_1.0_10.0");
    TH1F *h2 = (TH1F *)file2->Get("ksks_subtracted_invmass_pt_1.0_10.0");
    TH1F *h3 = (TH1F *)file3->Get("ksks_subtracted_invmass_pt_1.0_10.0");
    TH1F *h4 = (TH1F *)file4->Get("ksks_subtracted_invmass_pt_1.0_10.0");
    // TH1F *h5 = (TH1F *)file5->Get("ksks_subtracted_invmass_pt_1.0_10.0");

    if(h1 == nullptr || h2 == nullptr || h3 == nullptr || h4 == nullptr) {
        cout << "Error opening histogram" << endl;
        return;
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
    h1->GetXaxis()->SetRangeUser(1.00, 2.50);
    SetHistoQA(h1);
    h1->SetMarkerSize(0.8);
    h1->SetMarkerColor(colors[0]);
    h1->SetLineColor(colors[0]);
    h1->SetMarkerStyle(markers[0]);
    h1->Draw();
    SetHistoQA(h2);
    h2->GetXaxis()->SetRangeUser(1.00, 2.50);
    h2->SetMarkerSize(0.8);
    h2->SetMarkerColor(colors[1]);
    h2->SetLineColor(colors[1]);
    h2->SetMarkerStyle(markers[1]);
    h2->Draw("same");
    SetHistoQA(h3);
    h3->GetXaxis()->SetRangeUser(1.00, 2.50);
    h3->SetMarkerSize(0.8);
    h3->SetMarkerColor(colors[2]);
    h3->SetLineColor(colors[2]);
    h3->SetMarkerStyle(markers[2]);
    h3->Draw("same");
    SetHistoQA(h4);
    h4->GetXaxis()->SetRangeUser(1.00, 2.50);
    h4->SetMarkerSize(0.8);
    h4->SetMarkerColor(colors[3]);
    h4->SetLineColor(colors[3]);
    h4->SetMarkerStyle(markers[3]);
    // h4->Draw("same");
    // SetHistoQA(h5);
    // h5->GetXaxis()->SetRangeUser(1.00, 2.50);
    // h5->SetMarkerColor(colors[3]);
    // h5->SetLineColor(colors[3]);
    // h5->SetMarkerStyle(markers[3]);
    // // h5->Draw("same");

    TLegend *leg = new TLegend(0.50, 0.62, 0.92, 0.92);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    // leg->SetHeader("Normalization ranges");
    leg->SetHeader("K^{0}_{s} selection mass window");
    leg->AddEntry(h1, "2#sigma", "lpe");
    leg->AddEntry(h2, "3#sigma", "lpe");
    leg->AddEntry(h3, "4#sigma", "lpe");
    // leg->AddEntry(h4, "5#sigma", "lpe");
    // leg->AddEntry(h5, "2.70 - 2.80", "lpe");
    leg->Draw("same");

    c->SaveAs((path1 + "compare_norm_ranges.png").c_str());

}