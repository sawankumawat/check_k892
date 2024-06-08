#include <iostream>
#include "../src/style.h"
#include "../src/common.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);

    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.009);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.002);

    pad2->SetLeftMargin(0.1392758);
    pad2->SetRightMargin(0.009749304);
    pad2->SetTopMargin(0.1714559);
    pad2->SetBottomMargin(0.3304598);
}

void compare_rawyield()
{
    // published spectra
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TFile *fpub = new TFile("pp13TeV.root", "READ");
    TGraph *gr = (TGraph *)fpub->Get("Table 4/Graph1D_y1");
    // gr->Scale(1.0/0.7448);
    if (!gr)
    {
        cout << "Graph not found" << endl;
        return;
    }
    int spectrano = 0;
    //*********Root file*************
    // TFile *spectra = new TFile(("spectra_" + to_string(spectrano) + ".root").c_str(), "READ");
    TFile *fspectra = new TFile("/home/sawan/k892_postprocessing/output/LHC22o_pass6_small_INEL/common/spectra_1_LSBkg_0.root", "READ");
    TFile *fspectraqa = new TFile(("../" + kSignalOutput + "/yield.root").c_str(), "READ");
    if (fspectra->IsZombie() || fspectraqa->IsZombie())
    {
        cout << "File not found" << endl;
        return;
    }


    TH1D *h1 = (TH1D *)fspectra->Get("lf-k892analysis/K892/0/hraw1Yields");
    TH1D *hqa = (TH1D *)fspectraqa->Get("yield_integral");

    if (h1 == nullptr || hqa == nullptr)
    {
        cout << "Histogram not found" << endl;
        return;
    }
    TCanvas *csimple = new TCanvas();
    h1->Draw();
    TCanvas *csimple1 = new TCanvas();
    hqa->Draw();

    TH1F *hratio = (TH1F *)hqa->Clone();
    hratio->Divide(h1);
  

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    c1->SetLogy(1);
    gPad->SetLogy(1);
    h1->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
    h1->GetYaxis()->CenterTitle(1);
    h1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1->SetTitle(0);
    SetHistoStyle(h1, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    h1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    h1->GetXaxis()->SetTitleOffset(1.02);
    h1->GetYaxis()->SetTitleOffset(1.15);
    h1->SetMinimum(1e-6);
    h1->SetMarkerStyle(29);
    h1->SetMarkerSize(1);
    h1->SetMarkerColor(kRed);
    h1->SetLineColor(kRed);
    h1->GetXaxis()->SetRangeUser(0, 15);
    h1->Draw("pe");
    SetHistoStyle(hqa, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    hqa->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hqa->SetMaximum(hqa->GetMaximum() * 2);
    hqa->GetYaxis()->SetTitleOffset(1.15);
    hqa->GetXaxis()->SetTitleOffset(1.02);
    hqa->SetMarkerStyle(22);
    hqa->SetMarkerSize(1);
    hqa->SetLineColor(1);
    hqa->SetMarkerColor(1);
    hqa->SetMinimum(7e-6);
    hqa->Draw("pe same");
  

    TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.8);
    SetLegendStyle(leg);
    leg->AddEntry(gr, "raw yield (kstarqa code)", "lpe");
    leg->AddEntry(h1, "raw yield (k892analysis code)", "lpe");
    leg->SetTextSize(0.05);
    leg->Draw();

    c1->cd(2);

    SetHistoStyle(hratio, 1, 53, 1, 0.05, 0.05, 0.04 / pad2Size, 0.04 / pad2Size, 1.13, 1.8);
    hratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hratio->SetMarkerStyle(23);
    hratio->SetMarkerSize(1.0);
    hratio->SetMarkerColor(kRed);
    hratio->SetLineColor(kRed);
    hratio->GetYaxis()->SetTitle("#frac{kstarqa}{k892analysis}");
    hratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hratio->GetXaxis()->CenterTitle(1);
    hratio->GetYaxis()->SetTitleOffset(0.45);
    hratio->GetYaxis()->SetNdivisions(505);
    hratio->GetXaxis()->SetRangeUser(0, 15);
    hratio->Draw("ep");
}