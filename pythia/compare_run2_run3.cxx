#include <iostream>
#include "../macro/src/style.h"
#include "THnSparse.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void compare_run2_run3()
{
    TFile *f = new TFile("kstar_all.root", "read");
    THnSparseF *hsparse_run3 = (THnSparseF *)f->Get("hkstar_spraseFT0");
    THnSparseF *hsparse_run2 = (THnSparseF *)f->Get("hkstar_spraseFv0");
    if (hsparse_run3 == nullptr || hsparse_run2 == nullptr)
    {
        cout << "Error reading sparse histograms" << endl;
        return;
    }
    TH1D *hrun2pT = hsparse_run2->Projection(1);
    TH1D *hrun3pT = hsparse_run3->Projection(1);
    if (hrun2pT == nullptr || hrun3pT == nullptr)
    {
        cout << "Error projecting histograms" << endl;
        return;
    }
    TH1D *hratio = (TH1D *)hrun3pT->Clone();
    hratio->Divide(hrun2pT); // ratio

    TCanvas *cspectra = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cspectra, 0.14, 0.02, 0.02, 0.14);
    double pad1Size, pad2Size;
    canvas_style(cspectra, pad1Size, pad2Size);

    cspectra->cd(1);
    gPad->SetLogy();
    SetHistoQA(hrun3pT);
    SetHistoQA(hrun2pT);
    hrun2pT->SetLineColor(kBlue);
    hrun2pT->SetMarkerColor(kBlue);
    hrun3pT->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hrun3pT->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
    hrun3pT->GetYaxis()->SetTitleOffset(1.4);
    hrun3pT->GetXaxis()->SetTitleSize(0.042 / pad1Size);
    hrun3pT->GetYaxis()->SetTitleSize(0.042 / pad1Size);
    hrun3pT->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hrun3pT->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hrun3pT->GetXaxis()->SetRangeUser(0, 10);
    hrun3pT->Draw();
    hrun2pT->Draw("same");

    TLegend *leg = new TLegend(0.52, 0.55, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04 / pad1Size);
    leg->AddEntry(hrun3pT, "Run 3 (FT0)", "l");
    leg->AddEntry(hrun2pT, "Run 2 (FV0)", "l");
    leg->Draw();

    cspectra->cd(2);
    SetHistoQA(hratio);
    hratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hratio->GetYaxis()->SetTitle("#frac{Run 3}{Run 2}");
    hratio->GetYaxis()->SetTitleOffset(0.6);
    hratio->GetXaxis()->SetTitleOffset(1.1);
    hratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hratio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hratio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hratio->GetXaxis()->SetRangeUser(0, 10);
    hratio->GetYaxis()->SetRangeUser(0.945, 1.07);
    hratio->GetXaxis()->SetNdivisions(510);
    hratio->GetYaxis()->SetNdivisions(505);
    hratio->Draw("pe");
    TLine *line = new TLine(0, 1, 10, 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.05);
    pad2->SetRightMargin(0.05);
    pad2->SetBottomMargin(0.35);
    pad1->SetLeftMargin(0.175);
    pad2->SetLeftMargin(0.175);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.003);
    pad2->SetTopMargin(0.04);
}