#include <iostream>
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

void mc_closure()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    string path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/MC_closure/476106/kstarqa/hInvMass"; // path for yield.root file (from rec MC)
    string path2 = "/home/sawan/check_k892/data/kstar/LHC22o_pass7/MC_closure";                           // MC file path

    TFile *fspectra1 = new TFile((path1 + "/yield.root").c_str(), "read");
    TFile *fspectra2 = new TFile((path2 + "/476106.root").c_str(), "read");

    if (fspectra1->IsZombie() || fspectra2->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    // float mult_classes[] = {0};
    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    TH1F *hmult1[numofmultbins + 1];
    TH1F *hmult2[numofmultbins + 1];

    THnSparseF *hSparseRec = (THnSparseF *)fspectra2->Get("kstarqa/hInvMass/h2KstarRecpt2");
    // THnSparseF *hSparseRec = (THnSparseF *)fspectra2->Get("kstarqa/hInvMass/hk892GenpT");
    if (hSparseRec == nullptr)
    {
        cout << "Error reading efficiency histogram MC" << endl;
        return;
    }
    TH1D *h1recmult = (TH1D *)fspectra2->Get("kstarqa/hInvMass/h1RecMult");

    for (int imult = 0; imult < numofmultbins + 1; imult++)
    {
        double multlow = (imult == 0) ? 0 : mult_classes[imult - 1];
        double multhigh = (imult == 0) ? 100 : mult_classes[imult];

        // hmult1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));
        hmult1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/yield_bincount", multlow, multhigh));
        if (hmult1[imult] == nullptr)
        {
            cout << "Histogram hmult1 not found" << endl;
            return;
        }
        hmult2[imult] = (TH1F *)hmult1[imult]->Clone();

        int lowbinMultRec = hSparseRec->GetAxis(1)->FindBin(multlow + 1e-3);
        int highbinMultRec = hSparseRec->GetAxis(1)->FindBin(multhigh - 1e-3);
        hSparseRec->GetAxis(1)->SetRange(lowbinMultRec, highbinMultRec);

        TH1D *h1rec = hSparseRec->Projection(0, "E");
        int entries = h1recmult->Integral(h1recmult->GetXaxis()->FindBin(multlow + 1e-3), h1recmult->GetXaxis()->FindBin(multhigh - 1e-3));
        cout<<"multiplicity class "<<multlow<<" - "<<multhigh<<" : "<<entries<<endl;

        for (int i = 0; i < hmult1[imult]->GetNbinsX(); i++)
        {
            float lowpt = pT_bins[i];
            float highpt = pT_bins[i + 1];
            float ptbinwidth = highpt - lowpt;
            float BR = 0.76;

            double nrec = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt + 1e-3), h1rec->GetXaxis()->FindBin(highpt - 1e-3)) / (ptbinwidth * BR * entries);
            hmult2[imult]->SetBinContent(i + 1, nrec);
            hmult2[imult]->SetBinError(i + 1, 0);
        }

        TH1F *hratio1 = (TH1F *)hmult1[imult]->Clone(Form("ratio_mult_%.0f-%.0f", multlow, multhigh));
        hratio1->Divide(hmult2[imult]);

        TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
        SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
        double pad1Size, pad2Size;
        canvas_style(c1, pad1Size, pad2Size);
        c1->cd(1);
        gPad->SetLogy(1);
        SetHistoStyle(hmult1[imult], 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
        hmult1[imult]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        hmult1[imult]->SetMaximum(hmult1[imult]->GetMaximum() * 10);
        hmult1[imult]->SetMinimum(hmult1[imult]->GetMinimum() * 0.8);
        hmult1[imult]->GetYaxis()->SetTitleOffset(1.30);
        hmult1[imult]->GetXaxis()->SetTitleOffset(1.02);
        hmult1[imult]->SetMarkerStyle(20);
        hmult1[imult]->SetMarkerSize(1);
        hmult1[imult]->GetXaxis()->SetRangeUser(0, 10);
        hmult1[imult]->Draw("pe");
        hmult2[imult]->SetMarkerStyle(21);
        hmult2[imult]->SetMarkerSize(1);
        hmult2[imult]->SetMarkerColor(kBlue);
        hmult2[imult]->SetLineColor(kBlue);
        hmult2[imult]->SetLineWidth(2);
        hmult2[imult]->Draw("pe same");

        TLegend *leg = new TLegend(0.46, 0.64, 0.9, 0.91);
        SetLegendStyle(leg);
        leg->SetHeader(Form("Multiplicity: %.0f-%.0f%%", multlow, multhigh));
        leg->AddEntry(hmult1[imult], "Reconstructed Yield (Fit)", "lpe");
        leg->AddEntry(hmult2[imult], "True Yield (reconstructed K*)", "lpe");
        leg->SetTextSize(0.04);
        leg->Draw();

        c1->cd(2);
        TH1F *hdummy = (TH1F *)hmult1[0]->Clone();
        for (int i = 0; i < hdummy->GetNbinsX(); i++)
        {
            hdummy->SetBinContent(i + 1, 0);
            hdummy->SetBinError(i + 1, 0);
        }

        SetHistoQA(hratio1);
        hratio1->GetYaxis()->SetTitleSize(0.035 / pad2Size);
        hratio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        hratio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        hratio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        hratio1->SetMarkerStyle(20);
        hratio1->SetMarkerSize(1.0);
        // hratio1->SetMarkerColor(kBlue);
        // hratio1->SetLineColor(kBlue);
        hratio1->GetYaxis()->SetTitle("#frac{Rec Yield (Fit)}{True Yield}");
        hratio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hratio1->GetXaxis()->CenterTitle(1);
        hratio1->GetYaxis()->SetTitleOffset(0.6);
        hratio1->GetXaxis()->SetTitleOffset(1.1);
        hratio1->GetYaxis()->SetNdivisions(505);
        hratio1->SetMaximum(1.5);
        hratio1->SetMinimum(0.5);
        hratio1->GetXaxis()->SetRangeUser(0, 10);
        hratio1->Draw("p");

        TLine *line = new TLine(0, 1, 10, 1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(1);
        line->Draw();

        c1->SaveAs((path1 + Form("/MCclosure_%.0f-%.0f.png", multlow, multhigh)).c_str());
    }
}