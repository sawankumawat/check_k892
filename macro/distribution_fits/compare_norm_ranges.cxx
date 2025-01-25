#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/common_glue.h"
#include "../src/style.h"
using namespace std;

void compare_norm_ranges() {
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/";
    int colors[] = {2, 4, 6, 28, 46};
    int markers[] = {20, 21, 22, 23, 24};
    TFile *file1 = new TFile((path + "hglue_ROTATED_norm_2.30_2.40.root").c_str(), "READ");
    TFile *file2 = new TFile((path + "hglue_ROTATED_norm_2.40_2.50.root").c_str(), "READ");
    TFile *file3 = new TFile((path + "hglue_ROTATED_norm_2.50_2.60.root").c_str(), "READ");
    TFile *file4 = new TFile((path + "hglue_ROTATED_norm_2.60_2.70.root").c_str(), "READ");
    TFile *file5 = new TFile((path + "hglue_ROTATED_norm_2.70_2.80.root").c_str(), "READ");
    if(file1->IsZombie() || file2->IsZombie() || file3->IsZombie() || file4->IsZombie()) {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *h1 = (TH1F *)file1->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    TH1F *h2 = (TH1F *)file2->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    TH1F *h3 = (TH1F *)file3->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    TH1F *h4 = (TH1F *)file4->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    TH1F *h5 = (TH1F *)file5->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    if(h1 == nullptr || h2 == nullptr || h3 == nullptr || h4 == nullptr) {
        cout << "Error opening histogram" << endl;
        return;
    }
    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
    h1->GetXaxis()->SetRangeUser(1.00, 2.80);
    SetHistoQA(h1);
    h1->SetMarkerColor(colors[0]);
    h1->SetLineColor(colors[0]);
    h1->SetMarkerStyle(markers[0]);
    // h1->Draw();
    SetHistoQA(h2);
    h2->GetXaxis()->SetRangeUser(1.00, 2.80);
    h2->SetMarkerColor(colors[1]);
    h2->SetLineColor(colors[1]);
    h2->SetMarkerStyle(markers[1]);
    h2->Draw();
    SetHistoQA(h3);
    h3->GetXaxis()->SetRangeUser(1.00, 2.80);
    h3->SetMarkerColor(colors[2]);
    h3->SetLineColor(colors[2]);
    h3->SetMarkerStyle(markers[2]);
    h3->Draw("same");
    SetHistoQA(h4);
    h4->GetXaxis()->SetRangeUser(1.00, 2.80);
    h4->SetMarkerColor(colors[3]);
    h4->SetLineColor(colors[3]);
    h4->SetMarkerStyle(markers[3]);
    h4->Draw("same");
    SetHistoQA(h5);
    h5->GetXaxis()->SetRangeUser(1.00, 2.80);
    h5->SetMarkerColor(colors[3]);
    h5->SetLineColor(colors[3]);
    h5->SetMarkerStyle(markers[3]);
    // h5->Draw("same");

    TLegend *leg = new TLegend(0.60, 0.67, 0.92, 0.92);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->SetHeader("Normalization ranges");
    // leg->AddEntry(h1, "2.30 - 2.40", "lpe");
    leg->AddEntry(h2, "2.40 - 2.50", "lpe");
    leg->AddEntry(h3, "2.50 - 2.60", "lpe");
    leg->AddEntry(h4, "2.60 - 2.70", "lpe");
    // leg->AddEntry(h5, "2.70 - 2.80", "lpe");
    leg->Draw("same");

    c->SaveAs((path + "check/compare_norm_ranges.pdf").c_str());

}