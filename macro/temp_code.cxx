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
    if (h == nullptr)
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
    h->SetMaximum(1e11);
    h->Draw("hist text");
    c->SaveAs("distribution_fits/saved/NksProduced.png");

    // cout<<"Fraction of events in which more than 1 Ks are produced is "<<h->Integral(h->FindBin(2), h->GetNbinsX())/h->Integral()<<endl;

    string path_tpc_energyloss = "strangeness_tutorial/kzeroShort/dE_by_dx_TPC";
    TH2F *h_tpc_energyloss = (TH2F *)f->Get(path_tpc_energyloss.c_str());
    if (h_tpc_energyloss == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    TCanvas *c_tpc_energyloss = new TCanvas("c_tpc_energyloss", "c_tpc_energyloss", 720, 720);
    SetCanvasStyle(c_tpc_energyloss, 0.15, 0.13, 0.03, 0.14);
    SetHistoQA(h_tpc_energyloss);
    gPad->SetLogz();
    gPad->SetLogx();
    h_tpc_energyloss->GetXaxis()->SetRangeUser(0.1, 100);
    h_tpc_energyloss->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_tpc_energyloss->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
    h_tpc_energyloss->GetXaxis()->SetTitleOffset(1.3);
    h_tpc_energyloss->Draw("colz");
    c_tpc_energyloss->SaveAs("distribution_fits/saved/dE_by_dx_TPC.png");
}