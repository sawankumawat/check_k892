#include<iostream>
#include "../../src/style.h"
using namespace std;

void temp(){
    gStyle->SetOptStat(0);

    TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/341913/KsKs_Channel/higher-mass-resonances/hglue_ROTATED_norm_2.50_2.60.root", "READ");
    TFile *f2ks = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/341913/KsKs_Channel/higher-mass-resonances_twoKsOnly4sigma/hglue_ROTATED_norm_2.50_2.60.root", "READ");

    TH1F *h = (TH1F *)f->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    TH1F *h2ks = (TH1F *)f2ks->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    if(h == nullptr || h2ks == nullptr){
        cout << "Error opening histogram" << endl;
        return;
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.17, 0.03, 0.05, 0.15);
    SetHistoQA(h);
    SetHistoQA(h2ks);
    h->GetYaxis()->SetTitleOffset(1.7);
    h->GetXaxis()->SetRangeUser(0.99, 2.69);
    h2ks->GetXaxis()->SetRangeUser(0.99, 2.69);
    h->Draw("pe");
    h2ks->SetLineColor(2);
    h2ks->SetMarkerStyle(22);
    h2ks->SetMarkerColor(2);
    h2ks->Draw("pe same");

    TLegend *leg = new TLegend(0.4, 0.75, 0.9, 0.9);
    leg->SetTextFont(42);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(h, "All selected events", "lpe");
    leg->AddEntry(h2ks, "Selected events with two K^{0}_{s}", "lpe");
    leg->Draw();

    c->SaveAs("temp.png");

    
}