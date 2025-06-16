#include <iostream>
using namespace std;
#include "style.h"

void check_mc_local()
{
    // gStyle->SetOptStat(1110);
    gStyle->SetOptStat(0);
    TFile *f1 = new TFile("/home/sawan/alice/practice/AnalysisResults_f0_1710.root", "read");
    TFile *f2 = new TFile("/home/sawan/alice/practice/AnalysisResults_Lambda.root", "read");
    if (f1->IsZombie() || f2->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }
    TH2F *hptvsy_f1710 = (TH2F *)f1->Get("alice3-single-particle/particle/Ptvzy");
    TH2F *hptvsy_Lambda = (TH2F *)f2->Get("alice3-single-particle/particle/Ptvzy");
    TH2F *hptvsEta_f1710 = (TH2F *)f1->Get("alice3-single-particle/particle/PtvsEta");
    TH2F *hptvsEta_Lambda = (TH2F *)f2->Get("alice3-single-particle/particle/PtvsEta");
    if (!hptvsy_f1710 || !hptvsy_Lambda || !hptvsEta_f1710 || !hptvsEta_Lambda)
    {
        cout << "Error reading histograms" << endl;
        return;
    }

    SetHistoQA2(hptvsy_f1710);
    SetHistoQA2(hptvsy_Lambda);
    SetHistoQA2(hptvsEta_f1710);
    SetHistoQA2(hptvsEta_Lambda);

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.04);

    TCanvas *c1 = new TCanvas("c1", "Ptvsy f_{0}(1710)", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.10, 0.05, 0.14);
    hptvsy_f1710->Draw("colz");
    lat.DrawLatex(0.20, 0.80, "f_{0}(1710) candidate");
    c1->SaveAs("Ptvsy_f1710.png");
    

    TCanvas *c2 = new TCanvas("c2", "Ptvsy Lambda", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.10, 0.05, 0.14);
    hptvsy_Lambda->Draw("colz");
    lat.DrawLatex(0.20, 0.80, "#Lambda candidate");
    c2->SaveAs("Ptvsy_Lambda.png");

    TCanvas *c3 = new TCanvas("c3", "PtvsEta f_{0}(1710)", 720, 720);
    SetCanvasStyle(c3, 0.15, 0.10, 0.05, 0.14);
    hptvsEta_f1710->Draw("colz");
    lat.DrawLatex(0.20, 0.80, "f_{0}(1710) candidate");
    c3->SaveAs("PtvsEta_f1710.png");

    TCanvas *c4 = new TCanvas("c4", "PtvsEta Lambda", 720, 720);
    SetCanvasStyle(c4, 0.15, 0.10, 0.05, 0.14);
    hptvsEta_Lambda->Draw("colz");
    lat.DrawLatex(0.20, 0.80, "#Lambda candidate");
    c4->SaveAs("PtvsEta_Lambda.png");
}