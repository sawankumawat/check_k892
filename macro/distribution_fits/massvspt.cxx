#include "../src/style.h"
using namespace std;

void massvspt()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    double ptbins[] = {2.0, 3.0, 5.0, 8.0, 12.0};
    const int size = sizeof(ptbins) / sizeof(ptbins[0]);
    // double masses[] = {1.80074, 1.69425, 1.69193, 1.6828, 1.69281};
    // double masses_error[] = {0.0176463, 0.00896966, 0.00375625, 0.00322409, 0.00384157};
    double masses[] = {1.69425, 1.69193, 1.6828, 1.69281};
    double masses_error[] = {0.00896966, 0.00375625, 0.00322409, 0.00384157};
    double widths[] = {0.169626, 0.146556, 0.184768, 0.161962};
    double widths_error[] = {0.0534454, 0.0135399, 0.0110638, 0.0121257};

    TH1F *hmassvspt = new TH1F("hmassvspt", "Mass vs pT", size - 1, ptbins);
    TH1F *hwidthvspt = new TH1F("hwidthvspt", "Width vs pT", size - 1, ptbins);
    for (int i = 0; i < size - 1; i++)
    {
        hmassvspt->SetBinContent(i + 1, masses[i]);
        hmassvspt->SetBinError(i + 1, masses_error[i]);
        hwidthvspt->SetBinContent(i + 1, widths[i]);
        hwidthvspt->SetBinError(i + 1, widths_error[i]);
    }

    TH1F *hdummy = new TH1F("", "", 12, 1.5, 12.5);

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.16, 0.03, 0.05, 0.15);
    SetHistoQA(hdummy);
    SetHistoQA(hmassvspt);
    hdummy->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hdummy->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hdummy->GetYaxis()->SetRangeUser(1.655, 1.749);
    hdummy->GetYaxis()->SetTitleOffset(1.5);
    hmassvspt->SetMarkerStyle(20);
    hmassvspt->SetMarkerSize(1.0);
    hmassvspt->SetLineWidth(2);
    hmassvspt->SetMarkerColor(kBlack);
    hmassvspt->SetLineColor(kBlack);
    hdummy->Draw();
    hmassvspt->Draw("same E1");

    TLegend *leg = new TLegend(0.5, 0.65, 0.85, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg->AddEntry((TObject *)0, "FT0M, 0-100%", "");
    leg->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
    leg->AddEntry(hmassvspt, "f_{0}(1710) Mass", "p");
    leg->Draw();

    c->SaveAs("/home/sawan/Music/massvspt_data.png");

    TH1F *hdummy2 = new TH1F("", "", 12, 1.5, 12.5);
    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.16, 0.03, 0.05, 0.15);
    SetHistoQA(hdummy2);
    SetHistoQA(hwidthvspt);
    hdummy2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hdummy2->GetYaxis()->SetTitle("Width (GeV/#it{c}^{2})");
    hdummy2->GetYaxis()->SetRangeUser(0.051, 0.339);
    hdummy2->GetYaxis()->SetTitleOffset(1.5);
    hwidthvspt->SetMarkerStyle(20);
    hwidthvspt->SetMarkerSize(1.0);
    hwidthvspt->SetLineWidth(2);
    hwidthvspt->SetMarkerColor(kBlack);
    hwidthvspt->SetLineColor(kBlack);
    hdummy2->Draw();
    hwidthvspt->Draw("same E1");

    TLegend *leg2 = new TLegend(0.5, 0.65, 0.85, 0.9);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.035);
    leg2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg2->AddEntry((TObject *)0, "FT0M, 0-100%", "");
    leg2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
    leg2->AddEntry(hwidthvspt, "f_{0}(1710) Width", "p");
    leg2->Draw();

    c2->SaveAs("/home/sawan/Music/widthvspt_data.png");
}