#include <iostream>
using namespace std;
#include "src/style.h"

void temp_code()
{
    gStyle->SetOptStat(0);
    // gStyle->SetOptStat(1110);
    TFile *f = new TFile("../data/glueball/LHC22o_pass7_small/248267.root", "READ");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string path = "strangeness_tutorial/kzeroShort/NksProduced";
    TH1F *h = (TH1F *)f->Get(path.c_str());
    if(h == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    TCanvas *c = new TCanvas("c", "c", 720, 720);
    SetCanvasStyle(c, 0.15, 0.03, 0.03, 0.14);
    SetHistoQA(h);
    gPad->SetLogy();
    h->GetXaxis()->SetTitle("N_{K_{S}/Event}");
    h->GetYaxis()->SetTitle("Counts");
    h->Draw("hist");
    c->SaveAs("distribution_fits/saved/NksProduced.png");


}