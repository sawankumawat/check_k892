#include <iostream>
#include "TDatabasePDG.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"

TDatabasePDG *pdg = new TDatabasePDG();
double parameter0(double mass, double width);
TF1 *BW(TH1 *h, double mass, double width, double lowrange, double highrange);
TF1 *BW3expo(TH1 *h, double *parameters, double lowfitrange, double highfitrange);
TF1 *BW3pol3(TH1 *h, double *parameters, double lowfitrange, double highfitrange);
TF1 *BW3boltzman(TH1 *h, double *parameters, double lowfitrange, double highfitrange);
TF1 *BW3fit(TH1 *h, double *parameters, double lowfitrange, double highfitrange);

void rBW_fits()
{
    // constant parameters ****************************************
    const string kResBkg = "MIX";
    // const string kResBkg = "LIKE";
    const int rebin = 1;
    const string outputfolder_str = "../" + kSignalOutput + "/" + kchannel + "/" + kfoldername;
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    // double f1710Mass = pdg->GetParticle(10331)->Mass();
    // double f1710Width = pdg->GetParticle(10331)->Width();

    // pT loop ***************************************************
    for (Int_t ip = pt_start; ip < pt_end; ip++)
    {
        TFile *f = new TFile((outputfolder_str + "/hglue_" + kResBkg + Form("_%.1f_%.1f.root", pT_bins[ip], pT_bins[ip + 1])).c_str(), "READ");

        if (f->IsZombie())
        {
            cout << "Error opening file" << endl;
            return;
        }

        TH1F *hinvMass = (TH1F *)f->Get("ksks_invmass");
        if (hinvMass == nullptr)
        {
            cout << "Error opening histogram" << endl;
            return;
        }
        hinvMass->Rebin(rebin);

        TCanvas *c1 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c1, 0.12, 0.03, 0.05, 0.14);
        hinvMass->GetYaxis()->SetTitleOffset(1.1);
        hinvMass->GetXaxis()->SetTitleOffset(1.3);
        hinvMass->Draw();

        // Fitting *********************************************
        double parameters1[9] = {50, f1270Mass, f1270Width, 25, f1525Mass, f1525Width, 25, f1710Mass, f1710Width};
        TF1 *f3pol3 = BW3pol3(hinvMass, parameters1, 1.1, 2.2);
        f3pol3->Draw("same");

        //making the size of textbox optimal
        gPad->Update();
        TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
        ptstats->SetX1NDC(0.6);
        ptstats->SetX2NDC(0.99);
        ptstats->SetY1NDC(0.3);
        ptstats->SetY2NDC(0.95);
        ptstats->Draw("same");

        // individual functions drawing
        // BW1
        TF1 *frBW1 = new TF1("frBW1", RelativisticBW, 1.1, 2.2, 3);
        TF1 *frBW2 = new TF1("frBW2", RelativisticBW, 1.1, 2.2, 3);
        TF1 *frBW3 = new TF1("frBW3", RelativisticBW, 1.1, 2.2, 3);
        frBW1->SetParameter(0, f3pol3->GetParameter(0));
        frBW1->SetParameter(1, f3pol3->GetParameter(1));
        frBW1->SetParameter(2, f3pol3->GetParameter(2));
        // frBW1->SetLineWidth(2);
        frBW1->SetLineStyle(2);
        frBW1->SetLineColor(3);
        frBW1->Draw("same");

        // BW2
        frBW2->SetParameter(0, f3pol3->GetParameter(3));
        frBW2->SetParameter(1, f3pol3->GetParameter(4));
        frBW2->SetParameter(2, f3pol3->GetParameter(5));
        // frBW2->SetLineWidth(2);
        frBW2->SetLineStyle(2);
        frBW2->SetLineColor(4);
        frBW2->Draw("same");

        // BW3
        frBW3->SetParameter(0, f3pol3->GetParameter(6));
        frBW3->SetParameter(1, f3pol3->GetParameter(7));
        frBW3->SetParameter(2, f3pol3->GetParameter(8));
        // frBW3->SetLineWidth(2);
        frBW3->SetLineStyle(2);
        frBW3->SetLineColor(6);
        frBW3->Draw("same");

        //exponential background
        TF1 *fpol3 = new TF1("fpol3", polynomial3, 1.1, 2.2, 4);
        fpol3->SetParameter(0, f3pol3->GetParameter(9));
        fpol3->SetParameter(1, f3pol3->GetParameter(10));
        fpol3->SetParameter(2, f3pol3->GetParameter(11));
        fpol3->SetParameter(3, f3pol3->GetParameter(12));
        // fpol3->SetLineWidth(2);
        fpol3->SetLineStyle(2);
        fpol3->SetLineColor(8);
        fpol3->Draw("same");

        // // save plot ***********************************************
        TLegend *lfit = new TLegend(0.3, 0.65, 0.55, 0.94);
        lfit->SetFillColor(0);
        // lfit->SetBorderSize(0);
        lfit->SetFillStyle(0);
        lfit->SetTextFont(42);
        lfit->SetTextSize(0.04);
        lfit->AddEntry(hinvMass, "Data", "lpe");
        lfit->AddEntry(f3pol3, "3rBW + pol3", "l");
        lfit->AddEntry(frBW1, "rBW(f_{2}(1270))", "l");
        lfit->AddEntry(frBW2, "rBW(f_{2}(1525))", "l");
        lfit->AddEntry(frBW3, "rBW(f_{0}(1710))", "l");
        lfit->AddEntry(fpol3, "pol3", "l");
        t2->SetNDC();
        t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", 0.0, 30.0));
    }
}

double parameter0(double mass, double width)
{
    double gamma = mass * TMath::Sqrt((mass * mass + width * width));
    double norm = 2 * sqrt(2) * width * gamma / (TMath::Pi() * TMath::Sqrt(mass * mass + gamma));
    return norm;
}

TF1 *BW(TH1 *h, double mass, double width, double lowrange, double highrange)
{
    TF1 *f1bw = new TF1("f1bw", RelativisticBW, lowrange, highrange, 3);
    f1bw->SetParNames("norm", "mass", "width");
    std::cout << "mass = " << mass << ", width = " << width << std::endl;
    f1bw->SetParameter(0, parameter0(mass, width));
    f1bw->SetParameter(1, mass);
    f1bw->SetParLimits(1, mass - 1 * width, mass + 1 * width);
    f1bw->SetParameter(2, width);
    f1bw->SetParLimits(2, 0.0, 10);
    f1bw->SetLineStyle(1);
    f1bw->SetLineWidth(2);
    h->Fit("f1bw", "REBMS0");
    return f1bw;
}

TF1 *BW3expo(TH1 *h, double *parameters, double lowfitrange, double highfitrange)
{
    TF1 *f3bwexpo = new TF1("f3bwexpo", expo3bW, lowfitrange, highfitrange, 12);

   // fit names here
    f3bwexpo->SetParName(0, "norm1");
    f3bwexpo->SetParName(1, "mass1");
    f3bwexpo->SetParName(2, "width1");
    f3bwexpo->SetParName(3, "norm2");
    f3bwexpo->SetParName(4, "mass2");
    f3bwexpo->SetParName(5, "width2");
    f3bwexpo->SetParName(6, "norm3");
    f3bwexpo->SetParName(7, "mass3");
    f3bwexpo->SetParName(8, "width3");
    f3bwexpo->SetParName(9, "A");
    f3bwexpo->SetParName(10, "B");
    f3bwexpo->SetParName(11, "C");
    f3bwexpo->SetParName(12, "D");

    // fit parameters here
    f3bwexpo->SetParameter(0, parameters[0]);
    f3bwexpo->SetParameter(1, parameters[1]);
    f3bwexpo->SetParLimits(1, parameters[1] - 0.08, parameters[1] + 0.08);
    f3bwexpo->SetParameter(2, parameters[2]);
    f3bwexpo->SetParLimits(2, parameters[2] - 0.01, parameters[2] + 0.01);
    f3bwexpo->SetParameter(3, parameters[3]);
    f3bwexpo->SetParameter(4, parameters[4]);
    f3bwexpo->SetParameter(5, parameters[5]);
    f3bwexpo->SetParameter(6, parameters[6]);
    f3bwexpo->SetParameter(7, parameters[7]);
    f3bwexpo->SetParLimits(7, parameters[7] - 0.08, parameters[7] + 0.08);
    f3bwexpo->FixParameter(8, parameters[8]);
    f3bwexpo->SetLineStyle(1);
    f3bwexpo->SetLineWidth(2);

    h->Fit("f3bwexpo", "REBMS0");
    return f3bwexpo;
}

TF1 *BW3pol3(TH1 *h, double *parameters, double lowfitrange, double highfitrange)
{
    TF1 *f3pol3 = new TF1("f3pol3", pol3bW3, lowfitrange, highfitrange, 13);

    // fit names here
    f3pol3->SetParName(0, "norm1");
    f3pol3->SetParName(1, "mass1");
    f3pol3->SetParName(2, "width1");
    f3pol3->SetParName(3, "norm2");
    f3pol3->SetParName(4, "mass2");
    f3pol3->SetParName(5, "width2");
    f3pol3->SetParName(6, "norm3");
    f3pol3->SetParName(7, "mass3");
    f3pol3->SetParName(8, "width3");
    f3pol3->SetParName(9, "A");
    f3pol3->SetParName(10, "B");
    f3pol3->SetParName(11, "C");
    f3pol3->SetParName(12, "D");

    // fit parameters here
    f3pol3->SetParameter(0, parameters[0]);
    f3pol3->SetParameter(1, parameters[1]);
    f3pol3->SetParLimits(1, parameters[1] - 0.08, parameters[1] + 0.08);
    f3pol3->SetParameter(2, parameters[2]);
    f3pol3->SetParLimits(2, parameters[2] - 0.01, parameters[2] + 0.01);
    f3pol3->SetParameter(3, parameters[3]);
    f3pol3->SetParameter(4, parameters[4]);
    f3pol3->SetParameter(5, parameters[5]);
    f3pol3->SetParameter(6, parameters[6]);
    f3pol3->SetParameter(7, parameters[7]);
    f3pol3->SetParLimits(7, parameters[7] - 0.08, parameters[7] + 0.08);
    f3pol3->FixParameter(8, parameters[8]);
    f3pol3->SetLineStyle(1);
    f3pol3->SetLineWidth(2);

    h->Fit("f3pol3", "REBMS0");
    return f3pol3;
}

TF1 *BW3fit(TH1 *h, double *parameters, double lowfitrange, double highfitrange)
{
    TF1 *f3bw3 = new TF1("f3bw3", BW3, lowfitrange, highfitrange, 9);
    f3bw3->FixParameter(0, 50);
    f3bw3->FixParameter(1, parameters[1]);
    f3bw3->FixParameter(2, parameters[2]);
    f3bw3->FixParameter(3, 25);
    f3bw3->FixParameter(4, parameters[4]);
    f3bw3->FixParameter(5, parameters[5]);
    f3bw3->FixParameter(6, 25);
    f3bw3->FixParameter(7, parameters[7]);
    f3bw3->FixParameter(8, parameters[8]);
    return f3bw3;
}

TF1 *BW3boltzman(TH1 *h, double *parameters, double lowfitrange, double highfitrange)
{
    TF1 *f3bw = new TF1("f3bw", BW3boltzman, lowfitrange, highfitrange, 12);

    // fit names here
    f3bw->SetParName(0, "norm1");
    f3bw->SetParName(1, "mass1");
    f3bw->SetParName(2, "width1");
    f3bw->SetParName(3, "norm2");
    f3bw->SetParName(4, "mass2");
    f3bw->SetParName(5, "width2");
    f3bw->SetParName(6, "norm3");
    f3bw->SetParName(7, "mass3");
    f3bw->SetParName(8, "width3");
    f3bw->SetParName(9, "A");
    f3bw->SetParName(10, "B");
    f3bw->SetParName(11, "C");

    // fit parameters here
    f3bw->SetParameter(0, parameters[0]);
    f3bw->SetParameter(1, parameters[1]);
    f3bw->SetParameter(2, parameters[2]);
    f3bw->SetParameter(3, parameters[3]);
    f3bw->SetParameter(4, parameters[4]);
    f3bw->SetParameter(5, parameters[5]);
    f3bw->SetParameter(6, parameters[6]);
    f3bw->SetParameter(7, parameters[7]);
    f3bw->SetParameter(8, parameters[8]);
    f3bw->SetLineStyle(1);
    f3bw->SetLineWidth(2);

    h->Fit("f3bw", "REBMS0");
    return f3bw;
}
