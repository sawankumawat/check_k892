#include <iostream>
#include "TDatabasePDG.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"

void compare_channels()
{
    string channelKsKs = "KsKs_Channel";
    string channelKK = "KK_Channel";
    const string outputfolder_str1 = "../" + kSignalOutput + "/" + channelKsKs + "/strangeness_tutorial";
    const string outputfolder_str2 = "../" + kSignalOutput + "/" + channelKK + "/kaonkaonAnalysisRun3_Deep_angle";

    // pT loop ***************************************************
    const string resbkg1 = "MIX";
    const string resbkg2 = "LIKE";
    for (Int_t ip = pt_start; ip < pt_end; ip++)
    {
        gStyle->SetOptStat(0);
        TFile *fKsKs = new TFile((outputfolder_str1 + "/hglue_" + resbkg1 + Form("_%.1f_%.1f.root", pT_bins[ip], pT_bins[ip + 1])).c_str(), "READ");
        TFile *fKK = new TFile((outputfolder_str2 + "/hglue_" + resbkg2 + Form("_%.1f_%.1f.root", pT_bins[ip], pT_bins[ip + 1])).c_str(), "READ");
        if (fKsKs->IsZombie() || fKK->IsZombie())
        {
            cout << "File not found" << endl;
            return;
        }

        TH1F *hinvMassKsKs = (TH1F *)fKsKs->Get("ksks_invmass");
        TH1F *hinvMassKK = (TH1F *)fKK->Get("ksks_invmass");
        if (hinvMassKsKs == nullptr || hinvMassKK == nullptr)
        {
            cout << "Histogram not found" << endl;
            return;
        }

        TCanvas *c1 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c1, 0.12, 0.03, 0.05, 0.14);
        SetHistoQA(hinvMassKsKs);
        SetHistoQA(hinvMassKK);
        // hinvMassKsKs->Rebin(2);
        hinvMassKsKs->GetYaxis()->SetTitleOffset(1.1);
        hinvMassKsKs->GetXaxis()->SetTitleOffset(1.3);
        hinvMassKsKs->GetXaxis()->SetRangeUser(1.1, 2.8);
        hinvMassKsKs->SetMarkerSize(0.5);
        hinvMassKsKs->SetMarkerStyle(20);
        hinvMassKsKs->SetMarkerColor(kRed);
        hinvMassKsKs->SetLineColor(kRed);
        hinvMassKsKs->SetLineWidth(2);
        hinvMassKsKs->Scale(700);
        hinvMassKsKs->Draw("E");
        hinvMassKK->SetMarkerSize(0.5);
        hinvMassKK->SetMarkerStyle(20);
        hinvMassKK->SetMarkerColor(kBlue);
        hinvMassKK->SetLineColor(kBlue);
        hinvMassKK->SetLineWidth(2);
        hinvMassKK->Draw("E SAME");

        TLegend *leg = new TLegend(0.5, 0.8, 0.7, 0.9);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.045);
        leg->AddEntry(hinvMassKsKs, "KsKs #times 700", "lep");
        leg->AddEntry(hinvMassKK, "KK", "lep");
        leg->Draw("SAME");

        c1->SaveAs("compare_channels.png");
        
    }
}