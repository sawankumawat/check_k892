#include <iostream>
#include "../src/style.h"
using namespace std;

void compare_rebins()
{
    TFile *frebin1 = new TFile("saved/output_rebin1.root", "READ");
    TFile *frebin2 = new TFile("saved/output_rebin2.root", "READ");
    TFile *frebin3 = new TFile("saved/output_rebin3.root", "READ");

    if (frebin1->IsZombie() || frebin2->IsZombie() || frebin3->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string histnames[] = {"ks_mass_fit", "ks_width_fit", "ks_yield_int", "ks_yield_bin"};
    TH1F *hrebin1[4], *hrebin2[4], *hrebin3[4];

    for (int ih = 0; ih < 4; ih++)
    {
        hrebin1[ih] = (TH1F *)frebin1->Get(histnames[ih].c_str());
        hrebin2[ih] = (TH1F *)frebin2->Get(histnames[ih].c_str());
        hrebin3[ih] = (TH1F *)frebin3->Get(histnames[ih].c_str());
        if (hrebin1[ih] == nullptr || hrebin2[ih] == nullptr || hrebin3[ih] == nullptr)
        {
            cout << "Error reading histogram" << endl;
            return;
        }
        SetHistoQA(hrebin1[ih]);
        SetHistoQA(hrebin2[ih]);
        SetHistoQA(hrebin3[ih]);
        hrebin1[ih]->SetLineColor(1);
        hrebin1[ih]->SetMarkerColor(1);
        hrebin1[ih]->SetMarkerStyle(20);
        hrebin2[ih]->SetLineColor(2);
        hrebin2[ih]->SetMarkerColor(2);
        hrebin3[ih]->SetMarkerStyle(21);
        hrebin3[ih]->SetLineColor(4);
        hrebin3[ih]->SetMarkerColor(4);
        hrebin3[ih]->SetMarkerStyle(22);
    }

    TCanvas *cmass = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cmass, 0.18, 0.05, 0.03, 0.125);
    hrebin1[0]->GetYaxis()->SetRangeUser(0.494, 0.50);
    hrebin1[0]->GetYaxis()->SetTitleOffset(1.8);
    hrebin1[0]->GetYaxis()->SetMaxDigits(3);
    hrebin1[0]->Draw("E1");
    hrebin2[0]->Draw("E1 SAME");
    hrebin3[0]->Draw("E1 SAME");

    TLegend *leg = new TLegend(0.2, 0.6, 0.6, 0.9);
    SetLegendStyle(leg);
    leg->SetTextSize(0.045);
    leg->AddEntry(hrebin1[0], "Rebin 1", "ple");
    leg->AddEntry(hrebin2[0], "Rebin 2", "ple");
    leg->AddEntry(hrebin3[0], "Rebin 3", "ple");
    leg->Draw();
    cmass->SaveAs("saved/mass_rebins_comp.png");

    TCanvas *cwidth = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cwidth, 0.18, 0.05, 0.07, 0.125);
    // hrebin1[1]->GetYaxis()->SetRangeUser(0.012, 0.018);
    hrebin1[1]->GetYaxis()->SetTitleOffset(1.8);
    hrebin1[1]->GetYaxis()->SetMaxDigits(3);
    hrebin1[1]->Draw("E1");
    hrebin2[1]->Draw("E1 SAME");
    hrebin3[1]->Draw("E1 SAME");
    leg->Draw();
    cwidth->SaveAs("saved/width_rebins_comp.png");

    TCanvas *cyield = new TCanvas("", "", 720, 720);
    gPad->SetLogy();
    SetCanvasStyle(cyield, 0.18, 0.05, 0.07, 0.125);
    hrebin1[2]->GetYaxis()->SetTitleOffset(1.8);
    hrebin1[2]->GetYaxis()->SetMaxDigits(3);
    hrebin1[2]->Draw("E1");
    hrebin2[2]->Draw("E1 SAME");
    hrebin3[2]->Draw("E1 SAME");
    TLegend *leg2 = new TLegend(0.6, 0.6, 0.9, 0.9);
    SetLegendStyle(leg2);
    leg2->SetTextSize(0.045);
    leg2->AddEntry(hrebin1[2], "Rebin 1", "ple");
    leg2->AddEntry(hrebin2[2], "Rebin 2", "ple");
    leg2->AddEntry(hrebin3[2], "Rebin 3", "ple");
    leg2->Draw();
    cyield->SaveAs("saved/yield_rebins_comp.png");
}