#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"

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
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
}
void compare_dndy_meanpt()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // string path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/448490/kstarqa/hInvMass"; // 2022 data
    // string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/449695/kstarqa/hInvMass"; // 2023 data
    // string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/451003/kstarqa/hInvMass"; // 2024 data

    string path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/480447/kstarqa/hInvMass"; // 2023 data
    string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/480657/kstarqa/hInvMass"; // 2024 data
    string path3 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/480657/kstarqa/hInvMass"; // dummy data
    TString outputPath = path2;

    TFile *fspectra1 = new TFile((path1 + "/levy_fit.root").c_str(), "read");
    TFile *fspectra2 = new TFile((path2 + "/levy_fit.root").c_str(), "read");
    TFile *fspectra3 = new TFile((path3 + "/levy_fit.root").c_str(), "read");

    if (fspectra1->IsZombie() || fspectra2->IsZombie() || fspectra3->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    TGraphErrors *gMeanYieldRun2 = (TGraphErrors *)fspectra1->Get("gMeanYieldRun2");
    TGraphErrors *gMeanYieldRun3_file1 = (TGraphErrors *)fspectra1->Get("gMeanYieldRun3");
    TGraphErrors *gMeanYieldRun3_file2 = (TGraphErrors *)fspectra2->Get("gMeanYieldRun3");
    TGraphErrors *gMeanYieldRun3_file3 = (TGraphErrors *)fspectra3->Get("gMeanYieldRun3");
    TGraphErrors *gMeanpTRun2 = (TGraphErrors *)fspectra1->Get("gMeanpTRun2");
    TGraphErrors *gMeanpTRun3_file1 = (TGraphErrors *)fspectra1->Get("gMeanpTRun3");
    TGraphErrors *gMeanpTRun3_file2 = (TGraphErrors *)fspectra2->Get("gMeanpTRun3");
    TGraphErrors *gMeanpTRun3_file3 = (TGraphErrors *)fspectra3->Get("gMeanpTRun3");

    TCanvas *cMeanYield = new TCanvas("cMeanYield", "cMeanYield", 720, 720);
    SetCanvasStyle(cMeanYield, 0.15, 0.03, 0.03, 0.15);
    SetGrapherrorStyle(gMeanYieldRun2);
    gMeanYieldRun2->SetMarkerStyle(20);
    gMeanYieldRun2->SetMarkerSize(1.2);
    gMeanYieldRun2->SetMarkerColor(kRed);
    gMeanYieldRun2->SetLineColor(kRed);
    gMeanYieldRun2->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanYieldRun2->GetYaxis()->SetTitle("dN/dy");
    gMeanYieldRun2->Draw("AP");
    SetGrapherrorStyle(gMeanYieldRun3_file1);
    gMeanYieldRun3_file1->SetMarkerStyle(21);
    gMeanYieldRun3_file1->SetMarkerSize(1.2);
    gMeanYieldRun3_file1->SetMarkerColor(kBlue);
    gMeanYieldRun3_file1->SetLineColor(kBlue);
    gMeanYieldRun3_file1->Draw("P same");
    SetGrapherrorStyle(gMeanYieldRun3_file2);
    gMeanYieldRun3_file2->SetMarkerStyle(22);
    gMeanYieldRun3_file2->SetMarkerSize(1.2);
    gMeanYieldRun3_file2->SetMarkerColor(kGreen + 2);
    gMeanYieldRun3_file2->SetLineColor(kGreen + 2);
    gMeanYieldRun3_file2->Draw("P same");
    SetGrapherrorStyle(gMeanYieldRun3_file3);
    gMeanYieldRun3_file3->SetMarkerStyle(23);
    gMeanYieldRun3_file3->SetMarkerSize(1.2);
    gMeanYieldRun3_file3->SetMarkerColor(kMagenta);
    gMeanYieldRun3_file3->SetLineColor(kMagenta);
    // gMeanYieldRun3_file3->Draw("P same");
    TLegend *legMeanYield = new TLegend(0.2, 0.72, 0.45, 0.93);
    legMeanYield->SetTextSize(0.04);
    legMeanYield->SetBorderSize(0);
    legMeanYield->SetFillStyle(0);
    // legMeanYield->AddEntry(gMeanYieldRun3_file3, "2024 data", "p");
    legMeanYield->AddEntry(gMeanYieldRun3_file2, "2023 data", "p");
    legMeanYield->AddEntry(gMeanYieldRun3_file1, "2022 data", "p");
    legMeanYield->AddEntry(gMeanYieldRun2, "Run 2", "p");
    legMeanYield->Draw();
    cMeanYield->SaveAs(outputPath + "/mean_yield_compare.png");

    TCanvas *cMeanpT = new TCanvas("cMeanpT", "cMeanpT", 720, 720);
    SetCanvasStyle(cMeanpT, 0.15, 0.03, 0.03, 0.15);
    SetGrapherrorStyle(gMeanpTRun2);
    gMeanpTRun2->SetMarkerStyle(20);
    gMeanpTRun2->SetMarkerSize(1.2);
    gMeanpTRun2->SetMarkerColor(kRed);
    gMeanpTRun2->SetLineColor(kRed);
    gMeanpTRun2->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanpTRun2->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    gMeanpTRun2->SetMaximum(2.19);
    gMeanpTRun2->Draw("AP");
    SetGrapherrorStyle(gMeanpTRun3_file1);
    gMeanpTRun3_file1->SetMarkerStyle(21);
    gMeanpTRun3_file1->SetMarkerSize(1.2);
    gMeanpTRun3_file1->SetMarkerColor(kBlue);
    gMeanpTRun3_file1->SetLineColor(kBlue);
    gMeanpTRun3_file1->Draw("P same");
    SetGrapherrorStyle(gMeanpTRun3_file2);
    gMeanpTRun3_file2->SetMarkerStyle(22);
    gMeanpTRun3_file2->SetMarkerSize(1.2);
    gMeanpTRun3_file2->SetMarkerColor(kGreen + 2);
    gMeanpTRun3_file2->SetLineColor(kGreen + 2);
    gMeanpTRun3_file2->Draw("P same");
    SetGrapherrorStyle(gMeanpTRun3_file3);
    gMeanpTRun3_file3->SetMarkerStyle(23);
    gMeanpTRun3_file3->SetMarkerSize(1.2);
    gMeanpTRun3_file3->SetMarkerColor(kMagenta);
    gMeanpTRun3_file3->SetLineColor(kMagenta);
    // gMeanpTRun3_file3->Draw("P same");
    legMeanYield->Draw();
    cMeanpT->SaveAs(outputPath + "/mean_pT_compare.png");
}