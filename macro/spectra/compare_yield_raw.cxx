#include <iostream>
#include "../src/style.h"

using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void compare_yield_raw()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    string path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/658306/kstarqa_INELgt0/hInvMass";
    string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/661905/kstarqa/hInvMass";
    TFile *spectra1 = new TFile((path1 + "/yield_INEL.root").c_str(), "READ");
    TFile *spectra2 = new TFile((path2 + "/yield_INEL.root").c_str(), "READ");

    if (spectra1->IsZombie() || spectra2->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    TH1D *h1 = (TH1D *)spectra1->Get("mult_0-120/yield_integral");
    TH1D *h2 = (TH1D *)spectra2->Get("mult_0-120/yield_integral");

    if (h1 == nullptr || h2 == nullptr)
    {
        cout << "Histogram not found" << endl;
        return;
    }

    cout << "No. of bins in hist 1 is " << h1->GetNbinsX() << ", in hist 2 is " << h2->GetNbinsX() << endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    gPad->SetLogy();
    SetHistoQA(h1);
    h1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    h1->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    h1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    h1->SetMaximum(h1->GetMaximum() * 1.5);
    h1->GetYaxis()->SetTitleOffset(1.30);
    h1->GetXaxis()->SetTitleOffset(1.02);
    h1->SetMarkerStyle(22);
    h1->SetMarkerSize(1.5);
    h1->GetXaxis()->SetRangeUser(0, 15.5);
    h1->Draw("pe");
    SetHistoQA(h2);
    h2->SetMarkerStyle(24);
    h2->SetMarkerSize(1.5);
    h2->SetMarkerColor(kBlue);
    h2->SetLineColor(kBlue);
    h2->Draw("pe same");
    TLegend *leg = new TLegend(0.48, 0.5, 0.8, 0.8);
    SetLegendStyle(leg);
    leg->AddEntry(h1, "Hist 1", "lpe");
    leg->AddEntry(h2, "Hist 2", "lpe");
    leg->SetTextSize(0.05);
    leg->Draw();

    c1->cd(2);
    TH1D *hratio = (TH1D *)h2->Clone("hratio");
    hratio->Divide(h1);
    SetHistoQA(hratio);
    hratio->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    hratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hratio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hratio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hratio->GetYaxis()->SetTitle("h2 / h1");
    hratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hratio->GetYaxis()->SetTitleOffset(0.6);
    hratio->GetXaxis()->SetTitleOffset(1.1);
    hratio->GetYaxis()->SetNdivisions(506);
    hratio->GetYaxis()->SetRangeUser(0.9, 1.1);
    hratio->Draw("pe");

    TLine *line = new TLine(0, 1, 15, 1);
    line->SetLineStyle(2);
    line->Draw();
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}