#include <iostream>
#include <cmath>
#include "src/style.h"
#include "src/fitfunc.h"

using namespace std;
void SetHistogramProperties(TH1F *histogram, int color, int marker)
{
    histogram->SetMarkerColor(color);
    histogram->SetLineColor(color);
    histogram->SetLineWidth(2);
    histogram->SetMarkerSize(2);
    histogram->SetMarkerStyle(marker);
    histogram->SetTitle(0);
}

void compare_pbpb()
{
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    const int nof = 3;

    TFile *file[nof];
    TH1F *hmass[nof];
    TH1F *hwidth[nof];
    TH1F *hyield[nof];

    file[0] = new TFile("/home/sawan/check_k892/output/pbpb/kstar/LHC23zzf_pass2_QC/LHC23zzf_pass2_QC.root", "READ");
    file[1] = new TFile("/home/sawan/check_k892/output/pbpb/kstar/LHC23zzh_pass2_small/LHC23zzh_pass2_small.root", "READ");
    file[2] = new TFile("/home/sawan/check_k892/output/pbpb/kstar/LHC23zzg_apass2/LHC23zzg_apass2.root", "READ");

    string datanames[nof] = {"LHC23zzf_pass2_QC", "LHC23zzh_pass2_small", "LHC23zzg_apass2"};

    if (file[0]->IsZombie() || file[1]->IsZombie() || file[2]->IsZombie())
    {
        cerr << "File not found" << endl;
        return;
    }

    TCanvas *cmass = new TCanvas("", "", 1200, 1000);
    SetCanvasStyle2(cmass, 0.15, 0.05, 0.05, 0.15);
    int color[nof + 10] = {1, 2, 4, 6, 8};
    // int marker[nof + 10] = {53, 55, 54, 56, 57};
    // int marker[nof + 10] = {8, 71};
    int marker[3] = {8,8,8};

    // mass vs pT

    for (int i = 0; i < nof; i++)
    {
        hmass[i] = (TH1F *)file[i]->Get("mass");
        hwidth[i] = (TH1F *)file[i]->Get("width");
        hyield[i] = (TH1F *)file[i]->Get("yield");
        SetHistogramProperties(hmass[i], color[i], marker[i]);
        SetHistogramProperties(hwidth[i], color[i], marker[i]);
        SetHistogramProperties(hyield[i], color[i], marker[i]);
    }

    hmass[2]->SetBinContent(1, -999); // to remove the first bin from the LHC23zzh_pass2_small dataset
    hyield[2]->SetBinContent(1, -999);
     hmass[2]->SetBinContent(2, -999); // to remove the first bin from the LHC23zzh_pass2_small dataset
    hyield[2]->SetBinContent(2, -999);


    TH1F *dummy = new TH1F("", "", 16 * 5, -1, 16);
    SetHistoStyle(dummy, 1, 20, 1, 0.05, 0.05, 0.045, 0.045, 1.13, 1.8);
    dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    dummy->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    dummy->GetYaxis()->SetTitleOffset(1.4);

    dummy->SetMaximum(0.945);
    dummy->SetMinimum(0.855);
    dummy->GetXaxis()->SetRangeUser(-0.01, 11);
    dummy->Draw();
    for (int i = 0; i < nof; i++)
    {
        hmass[i]->Draw("same");
        /* code */
    }

    TLegend *massleg = new TLegend(0.45, 0.6, 0.9, 0.91);
    massleg->SetFillStyle(0);
    massleg->SetBorderSize(0);
    massleg->SetTextFont(42);
    massleg->SetTextSize(0.04);
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    // TLine *line = new TLine(hmass[1]->GetXaxis()->GetXmin(), 0.895, hmass[1]->GetXaxis()->GetXmax(), 0.895);
    TLine *line = new TLine(-0.01, 0.895, 16.1, 0.895);
    SetLineStyle(line, 28);
    line->SetLineWidth(3);
    massleg->AddEntry(line, "PDG Mass (0.895 GeV/c^{2})", "l");
    for (int i = 0; i < nof; i++)
    {
        massleg->AddEntry(hmass[i], datanames[i].c_str(), "lep");
    }
    // massleg->AddEntry(hmass[0], "pp 13.6 TeV low IR", "lep");
    // massleg->AddEntry(hmass[1], "pp 13.6 TeV high IR", "lep");
    // massleg->AddEntry(hmass[2], "Pb-Pb 5.36 TeV", "lep");
    massleg->Draw("l");
    cmass->SaveAs("/home/sawan/check_k892/output/compare/mass.png");

    // yield plots
    TCanvas *cyield = new TCanvas("", "", 1200, 1000);
    cyield->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)cyield->GetPad(1);
    TPad *pad2 = (TPad *)cyield->GetPad(2);
    double pad2Size = 0.3; // Size of the first pad
    double pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.05);
    pad2->SetRightMargin(0.05);
    pad2->SetBottomMargin(0.3);
    pad1->SetLeftMargin(0.12);
    pad2->SetLeftMargin(0.12);

    cyield->cd(1);
    SetCanvasStyle2(cyield, 0.15, 0.05, 0.05, 0.15);
    gPad->SetLogy();
    TH1F *dummy3 = new TH1F("", "", 17 * 5, -1, 16);
    SetHistoStyle(dummy3, 1, 20, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    dummy3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    dummy3->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    dummy3->GetYaxis()->SetTitleOffset(1);
    dummy3->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    dummy3->SetMaximum(5);
    dummy3->SetMinimum(15e-6);
    dummy3->GetXaxis()->SetRangeUser(-0.01, 11);
    dummy3->Draw();
    for (int i = 0; i < nof; i++)
    {
        hyield[i]->Draw("epsame");
    }

    // hyield[0]->Draw("same");
    // hyield[2]->Draw("same");
    TLegend *yieldleg = new TLegend(0.48, 0.6, 0.9, 0.91);
    yieldleg->SetFillStyle(0);
    yieldleg->SetBorderSize(0);
    yieldleg->SetTextFont(42);
    yieldleg->SetTextSize(0.04/pad1Size);
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");

    for (int i = 0; i < nof; i++)
    {
        yieldleg->AddEntry(hmass[i], datanames[i].c_str(), "lep");
    }

    // yieldleg->AddEntry(hmass[0], "pp 13.6 TeV low IR", "lep");
    // yieldleg->AddEntry(hmass[1], "pp 13.6 TeV high IR", "lep");
    // yieldleg->AddEntry(hmass[2], "Pb-Pb 5.36 TeV", "lep");
    yieldleg->Draw("l");

    // ratio of yields
    // TCanvas *cyield_ratio = new TCanvas("", "", 1200, 1000);
    // SetCanvasStyle2(cyield_ratio, 0.15, 0.05, 0.05, 0.15);
    cyield->cd(2);
    TH1F *dummy4 = new TH1F("", "", 17 * 5, -1, 16);
    TH1F *hratio = new TH1F("", "", hyield[0]->GetNbinsX(), hyield[0]->GetXaxis()->GetXmin(), hyield[0]->GetXaxis()->GetXmax());
    SetHistoStyle(dummy4, 1, 20, 1, 0.05, 0.05, 0.04 / pad2Size, 0.04 / pad2Size, 1.13, 1.8);
    dummy4->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    dummy4->GetYaxis()->SetTitle("Ratio");
    dummy4->GetYaxis()->SetTitleOffset(0.45);
    dummy4->GetXaxis()->SetTitleOffset(1);
    dummy4->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    dummy4->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    dummy4->SetMaximum(1.3);
    dummy4->SetMinimum(0.5);
    dummy4->GetXaxis()->SetRangeUser(-0.01, 11);
    dummy4->GetYaxis()->SetNdivisions(505); // Reduce the number of divisions in the y-axis

    dummy4->Draw();
    for (int i = 0; i < nof; i++)
    {
        TH1F *hclone = (TH1F *)hyield[i]->Clone();
        hclone->Divide(hyield[0]);
        hclone->SetMarkerStyle(8);
        hclone->SetMarkerColor(color[i]);
        hclone->SetMarkerSize(2);
        hclone->SetLineColor(color[i]);
        hclone->Draw("epsame");
    }
    // yieldleg->Draw("l");
    cyield->SaveAs("/home/sawan/check_k892/output/compare/yield.png");

    //  for (int i = 0; i < hyield[0]->GetNbinsX(); i++)
    // {
    //     float ratio = hyield[1]->GetBinContent(i + 1) / hyield[0]->GetBinContent(i + 1);
    //     float ratio_err = ratio * sqrt(pow(hyield[0]->GetBinError(i + 1) / hyield[0]->GetBinContent(i + 1), 2) + pow(hyield[1]->GetBinError(i + 1) / hyield[1]->GetBinContent(i + 1), 2));
    //     cout << "ratio " << ratio << endl;
    //     hratio->SetBinContent(i + 1, ratio);
    //     hratio->SetBinError(i + 1, ratio_err);
    // }

    // hratio->SetMarkerStyle(53);
    // hratio->SetMarkerColor(1);
    // hratio->SetLineColor(1);
    // hratio->SetLineWidth(2);
    // hratio->SetMarkerSize(2);
    // hratio->Draw("same");
    // TLatex lat;
    // lat.SetNDC();
    // lat.SetTextSize(0.06);
    // lat.SetTextFont(42);
    // lat.DrawLatex(0.5, 0.88, "pp 13.6 TeV");
    // cyield_ratio->SaveAs("/home/sawan/check_k892/output/compare/yield_ratio.png");
}