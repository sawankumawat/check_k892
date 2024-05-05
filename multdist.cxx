#include <iostream>
using namespace std;
#include "macro/src/style.h"

void multdist()
{
    gStyle->SetOptStat(0);
    TFile *f = new TFile("AnalysisResults.root");
    TH1F *multdist = (TH1F *)f->Get("kstarqa/histos/multdist");
    TH1F *multint = (TH1F *)f->Get("kstarqa/eventSelection/hmult");
    if (!multdist || !multint)
    {
        cout << "Error: histogram not found" << endl;
        return;
    }
    SetHistoQA(multdist);
    SetHistoQA(multint);
    multdist->GetXaxis()->SetTitle("Multiplicity");
    multdist->GetYaxis()->SetTitle("No. of events");
    multdist->GetXaxis()->SetNdivisions(506);
    TCanvas *c = new TCanvas("c", "c", 720, 720);
    SetCanvasStyle(c, 0.15, 0.07, 0.05, 0.13);
    c->SetLogy();
    multdist->Draw();
    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.05, 0.05, 0.13);
    multint->GetXaxis()->SetTitle("Multiplicity percentile");
    multint->GetYaxis()->SetTitle("No. of events");
    for (int i = 0; i < multint->GetNbinsX(); i++)
    {
        multint->SetBinContent(i, multint->GetBinContent(i));
        if (i > multint->GetNbinsX() - 10)
            multint->SetBinContent(i, 0);
    }

    multint->Draw();

    c->SaveAs("multdist.png");
    c1->SaveAs("multint.png");
}