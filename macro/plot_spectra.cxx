#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/initializations.h"

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
    pad1->SetRightMargin(0.02);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.31);
    pad1->SetLeftMargin(0.12);
    pad2->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.04);
}

void plot_spectra()
{
    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);
    string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/447406/kstarqa_id21631/hInvMass";

    TFile *fspectra = new TFile((path + "/corrected_spectra.root").c_str(), "read");

    if (fspectra->IsZombie())
    {
        cout << "File not found" << endl;
        return;
    }
    int markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 32, 47};
    float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;
    TH1F *hmult[numofmultbins + 1];
    // hmult[0] = (TH1F *)fspectra->Get("mult_0-100/yield_integral");
    hmult[0] = (TH1F *)fspectra->Get("mult_0-100/corrected_spectra_Integral");

    if (hmult[0] == nullptr)
    {
        cout << "Histogram 1 not found" << endl;
        return;
    }

    for (int i = 1; i < numofmultbins + 1; i++)
    {
        hmult[i] = (TH1F *)fspectra->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", mult_classes[i - 1], mult_classes[i]));
        if (hmult[i] == nullptr)
        {
            cout << "Histogram others not found" << endl;
            return;
        }
    }

    TH1F *hratios[numofmultbins];
    for (int i = 1; i < numofmultbins + 1; i++)
    {
        hratios[i] = (TH1F *)hmult[i]->Clone(Form("hratio%d", i));
        hratios[i]->Divide(hmult[0]);
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c, pad1Size, pad2Size);
    c->cd(1);
    c->SetLogy();
    gPad->SetLogy();
    SetHistoQA(hmult[0]);
    hmult[1]->GetYaxis()->SetTitle(hmult[1]->GetYaxis()->GetTitle());
    hmult[1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    int multiplicationFactors[] = {9, 8, 7, 5, 4, 3, 2, 1, 0}; // 6 is skipped because the multiplicity bin 10-15 and 15-20 are merged (in run 2)

    for (int i = 1; i < numofmultbins + 1; i++)
    {
        SetHistoQA(hmult[i]);
        hmult[i]->SetMarkerSize(1.2);
        // hmult[i]->Scale(pow(2, numofmultbins + 1 - i));
        // hmult[i]->Scale(pow(2, numofmultbins - i));
        hmult[i]->Scale(pow(2, multiplicationFactors[i - 1]));

        hmult[i]->GetXaxis()->SetTitleSize(0.045);
        hmult[i]->GetYaxis()->SetTitleSize(0.045);
        hmult[i]->GetYaxis()->SetTitleOffset(1.3);
        hmult[i]->SetMaximum(hmult[1]->GetMaximum() * 15);
        hmult[i]->SetMinimum(2e-6);
        hmult[i]->SetMarkerStyle(markers[i - 1]);
        hmult[i]->Draw("pe same PLC PMC");
    }
    hmult[0]->SetMarkerStyle(markers[numofmultbins]);
    hmult[0]->SetMarkerSize(1.2);
    hmult[0]->Draw("pe same PLC PMC");

    TLegend *leg = new TLegend(0.50, 0.65, 0.98, 0.98);
    SetLegendStyle(leg);
    leg->SetNColumns(3);
    leg->AddEntry(hmult[0], "0-100%", "lpe");
    for (int i = 1; i < numofmultbins + 1; i++)
    {
        // leg->AddEntry(hmult[i], Form("%.0f-%.0f%%(#times2^{%d})", mult_classes[i - 1], mult_classes[i], numofmultbins + 1 - i), "lpe");
        // leg->AddEntry(hmult[i], Form("%.0f-%.0f%%(#times2^{%d})", mult_classes[i - 1], mult_classes[i], numofmultbins - i), "lpe");
        leg->AddEntry(hmult[i], Form("%.0f-%.0f%%(#times2^{%d})", mult_classes[i - 1], mult_classes[i], multiplicationFactors[i - 1]), "lpe");
    }
    leg->SetTextSize(0.03);
    leg->Draw();

    c->cd(2);
    gPad->SetLogy(1);
    hratios[1]->SetMinimum(0.1);
    hratios[1]->SetMaximum(10);
    hratios[1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hratios[1]->GetYaxis()->SetTitle("Ratio to 0-100%");
    hratios[1]->GetYaxis()->SetNdivisions(505);
    // hratios[0]->GetXaxis()->SetRangeUser(0, 25);

    for (int i = 1; i < numofmultbins; i++)
    {
        SetHistoQA(hratios[i]);
        hratios[i]->GetXaxis()->SetTitleSize(0.045 * pad1Size / pad2Size);
        hratios[i]->GetYaxis()->SetTitleSize(0.045 * pad1Size / pad2Size);
        hratios[i]->GetXaxis()->SetLabelSize(0.04 * pad1Size / pad2Size);
        hratios[i]->GetYaxis()->SetLabelSize(0.04 * pad1Size / pad2Size);
        hratios[i]->GetYaxis()->SetTitleOffset(0.55);
        hratios[i]->SetMarkerStyle(markers[i]);
        hratios[i]->Draw("pe same PLC PMC");
    }
    TLine *line = new TLine(0, 1, 20, 1);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(3);
    line->Draw();

    TString outputfolder = kSignalOutput + "/" + kfoldername;
    c->SaveAs(outputfolder + "/spectra.png");
}