#include <iostream>
#include "../src/style.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void compare_rebins()
{
    // Open files and check for errors
    TFile *files[] = {
        new TFile("saved/pass6/output_rebin1.root", "READ"),
        new TFile("saved/pass6/output_rebin2.root", "READ"),
        new TFile("saved/pass6/output_rebin3.root", "READ")};

    for (auto &file : files)
    {
        if (file->IsZombie())
        {
            cout << "Error opening file" << endl;
            return;
        }
    }

    // Define histogram names and containers
    string histnames[] = {"ks_mass_fit", "ks_width_fit", "ks_yield_int", "ks_yield_bin"};
    TH1F *hist_rebin[3][4], *hratio[3][4];
    int colors[] = {kBlack, kRed, kBlue};

    // Load histograms and set styles
    for (int itype = 0; itype < 4; itype++)
    {
        for (int irebin = 0; irebin < 3; irebin++)
        {
            hist_rebin[irebin][itype] = (TH1F *)files[irebin]->Get(histnames[itype].c_str());
            if (!hist_rebin[irebin][itype])
            {
                cout << "Error reading histogram" << endl;
                return;
            }
            SetHistoQA(hist_rebin[irebin][itype]);
            hist_rebin[irebin][itype]->SetLineColor(colors[irebin]);
            hist_rebin[irebin][itype]->SetMarkerColor(colors[irebin]);
            hist_rebin[irebin][itype]->SetMarkerStyle(20 + irebin);
        }
    }

    // Create ratio histograms
    for (int itype = 0; itype < 4; itype++)
    {
        hratio[0][itype] = (TH1F *)hist_rebin[0][itype]->Clone();
        hratio[1][itype] = (TH1F *)hist_rebin[1][itype]->Clone();
        hratio[2][itype] = (TH1F *)hist_rebin[2][itype]->Clone();

        hratio[0][itype]->Divide(hist_rebin[0][itype]);
        hratio[1][itype]->Divide(hist_rebin[0][itype]);
        hratio[2][itype]->Divide(hist_rebin[0][itype]);
    }

    // Draw histograms and save to files
    string canvas_titles[] = {"mass_rebins_comp.png", "width_rebins_comp.png", "yield_rebins_comp.png"};
    double yaxis_ranges[] = {0.494, 0, 0}; // Adjust as necessary
    bool log_scale[] = {false, false, true};
    float yaxis_range_ratio[3][2] = {{0.998, 1.0019}, {0.9, 1.09}, {0.6, 1.39}};

    for (int itype = 0; itype < 3; itype++)
    {
        TCanvas *c = new TCanvas("", "", 720, 720);
        if (itype == 0)
            SetCanvasStyle(c, 0.29, 0.05, 0.07, 0.125);
        else
            SetCanvasStyle(c, 0.18, 0.05, 0.07, 0.125);
        double pad1Size, pad2Size;
        canvas_style(c, pad1Size, pad2Size);
        c->cd(1);
        if (log_scale[itype])
            gPad->SetLogy();

        (itype == 0) ? hist_rebin[0][itype]->GetYaxis()->SetTitleOffset(1.4): hist_rebin[0][itype]->GetYaxis()->SetTitleOffset(1.2);
        if (yaxis_ranges[itype] > 0)
            hist_rebin[0][itype]->GetYaxis()->SetRangeUser(yaxis_ranges[itype], 0.50);
        hist_rebin[0][itype]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
        hist_rebin[0][itype]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        hist_rebin[0][itype]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
        hist_rebin[0][itype]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
        hist_rebin[0][itype]->GetYaxis()->SetMaxDigits(3);
        hist_rebin[0][itype]->Draw("E1");
        hist_rebin[1][itype]->Draw("E1 SAME");
        hist_rebin[2][itype]->Draw("E1 SAME");

        TLegend *leg;
        if (itype == 2) // Third iteration
        {
            leg = new TLegend(0.5, 0.6, 0.9, 0.9);
        }
        else
        {
            leg = new TLegend(0.2, 0.6, 0.6, 0.9);
        }
        SetLegendStyle(leg);
        leg->SetTextSize(0.045);
        leg->AddEntry(hist_rebin[0][itype], "Rebin 1", "ple");
        leg->AddEntry(hist_rebin[1][itype], "Rebin 2", "ple");
        leg->AddEntry(hist_rebin[2][itype], "Rebin 3", "ple");
        leg->Draw();

        c->cd(2);
        for (int irebin = 0; irebin < 3; irebin++)
        {

            (itype == 0) ? hratio[irebin][itype]->GetYaxis()->SetTitleOffset(0.58) : hratio[irebin][itype]->GetYaxis()->SetTitleOffset(0.5);
            hratio[irebin][itype]->GetYaxis()->SetRangeUser(yaxis_range_ratio[itype][0], yaxis_range_ratio[itype][1]);
            hratio[irebin][itype]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
            hratio[irebin][itype]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
            hratio[irebin][itype]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
            hratio[irebin][itype]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
            hratio[irebin][itype]->GetYaxis()->SetNdivisions(505);
            hratio[irebin][itype]->GetYaxis()->SetMaxDigits(3);
            hratio[irebin][itype]->GetYaxis()->SetTitle("Rebin1/ RebinX");
            hratio[irebin][itype]->GetYaxis()->CenterTitle(0);
            hratio[irebin][itype]->Draw("E1 SAME");
        }

        c->SaveAs(("saved/pass6/" + canvas_titles[itype]).c_str());
    }
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.17, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.02);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.15);
    pad2->SetLeftMargin(0.15);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.0001);
    pad2->SetTopMargin(0.001);
}
