#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <TArrow.h>
#include "../../src/common_glue.h"
#include "../../src/fitting_range_glue.h"
#include "../../src/style.h"
using namespace std;

void plot_ressub()
{

    // TFile *file1 = new TFile("3rBW_plots_boltzmann.root", "READ");
    // TFile *file2 = new TFile("3rBW_plots_expol.root", "READ");
    // TFile *file3 = new TFile("3rBW_plots_exponential.root", "READ");

    // TFile *file1 = new TFile("4rBW_plots_boltzmann.root", "READ");
    TFile *file1 = new TFile("3rBW_plots_expol.root", "READ");
    TFile *file2 = new TFile("4rBW_plots_expol.root", "READ");
    TFile *file3 = new TFile("4rBW_plots_exponential.root", "READ");

    if (file1->IsZombie() || file2->IsZombie() || file3->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hsubtracted1 = (TH1F *)file1->Get("3BW");
    TH1F *hsubtracted2 = (TH1F *)file2->Get("3BW");
    TH1F *hsubtracted3 = (TH1F *)file3->Get("3BW");

    if (hsubtracted1 == nullptr || hsubtracted2 == nullptr || hsubtracted3 == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }
    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.145, 0.03, 0.05, 0.14);
    hsubtracted1->GetXaxis()->SetRangeUser(1.0, 2.5);
    hsubtracted1->SetMaximum(hsubtracted1->GetMaximum() * 1.5);
    hsubtracted1->SetLineColor(kRed);
    hsubtracted1->SetMarkerColor(kRed);
    hsubtracted1->SetStats(0);
    hsubtracted1->SetMaximum(hsubtracted1->GetMaximum() * 0.7);
    hsubtracted1->SetMarkerStyle(22);
    hsubtracted1->Draw("pe");
    hsubtracted2->SetLineColor(kBlue);
    hsubtracted2->SetMarkerColor(kBlue);
    hsubtracted2->SetMarkerStyle(23);
    hsubtracted2->SetStats(0);
    hsubtracted2->Draw("pe SAME");
    hsubtracted3->SetLineColor(kGreen);
    hsubtracted3->SetMarkerColor(kGreen);
    hsubtracted3->SetStats(0);
    // hsubtracted3->Draw("pe SAME");

    TLegend *leg = new TLegend(0.5, 0.67, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->AddEntry(hsubtracted1, "3rBW + expol1", "lpe");
    leg->AddEntry(hsubtracted2, "4rBW + expol1", "lpe");
    // leg->AddEntry(hsubtracted3, "4rBW + expol1", "lpe");
    leg->Draw("same");

    c->SaveAs("plot_ressub_compare.png");
}