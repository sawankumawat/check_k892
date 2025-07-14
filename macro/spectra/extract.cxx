#include<iostream>
#include "style.h"
using namespace std;

void extract() {
    gStyle->SetOptStat(0);
    TFile *f = new TFile("NKS_ME.root", "read");
    TH1F *hmass = (TH1F*)f->Get("hmass");
    TH1F *hwidth = (TH1F*)f->Get("hwidth");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    SetCanvasStyle(c1, 0.17, 0.03, 0.03, 0.13);

    SetHistoQA(hmass);
    SetHistoQA(hwidth);
    hmass->Draw("ep");
    hmass->SetMarkerSize(2);
    hmass->GetYaxis()->SetRangeUser(0.885, 0.905);
    hmass->GetXaxis()->SetTitleOffset(1.1);
    hmass->GetYaxis()->SetTitleOffset(1.7);

    TLegend *massleg = DrawLegend(0.6, 0.7, 0.9, 0.9);
    massleg->SetTextFont(42);
    massleg->SetTextSize(0.04);
    massleg->AddEntry(hmass, "pp 13.6 TeV", "l");
    TLine *line = new TLine(hmass->GetXaxis()->GetXmin(), 0.895, hmass->GetXaxis()->GetXmax(), 0.895);
    line->SetLineStyle(2);
    line->SetLineColor(4);
    line->SetLineWidth(3);
    line->Draw();
    massleg->AddEntry(line, "PDG Mass", "l");
    massleg->Draw("l");
    c1->SaveAs("mass.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
    SetCanvasStyle(c2, 0.13, 0.03, 0.03, 0.13);
    SetHistoQA(hwidth);
    hwidth->Draw("ep");
    hwidth->SetMarkerSize(2);
    // hwidth->GetYaxis()->SetRangeUser(0.001, 0.01);
    hwidth->GetXaxis()->SetTitleOffset(1.1);
    hwidth->GetYaxis()->SetTitleOffset(1.3);


}