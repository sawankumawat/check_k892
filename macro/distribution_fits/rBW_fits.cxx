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
TF1 *draw_individual_functions(double *parameters, TLegend *lfit);

void rBW_fits()
{
    // // *********************** constant parameters *****************************
    const string kResBkg = "MIX";
    // const string kResBkg = "ROTATED";
    // const string kResBkg = "LIKE";
    const string kbgfitfunction = "pol3";
    // const string kbgfitfunction = "expol";
    // const string kbgfitfunction = "Boltzman";

    const int rebin = 1;
    bool testing = false;
    bool saveplots = false;
    const string outputfolder_str = "../" + kSignalOutput + "/" + kchannel + "/" + kfoldername;
    // double f1710Mass = pdg->GetParticle(10331)->Mass();
    // double f1710Width = pdg->GetParticle(10331)->Width();
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    // // *********************** constant parameters *****************************

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
        float binwidth = hinvMass->GetBinWidth(1);
        hinvMass->GetYaxis()->SetTitle(Form("Counts / %.2f MeV/c^{2}", binwidth));

        TCanvas *c1 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c1, 0.12, 0.03, 0.05, 0.14);
        hinvMass->GetYaxis()->SetTitleOffset(1.1);
        hinvMass->GetXaxis()->SetTitleOffset(1.3);
        hinvMass->GetXaxis()->SetRangeUser(1.0, 2.4);
        TH1F *hinvMassResSub = (TH1F *)hinvMass->Clone();
        hinvMass->Draw();
        double parameters1[9] = {50, f1270Mass, f1270Width, 25, f1525Mass, f1525Width, 25, f1710Mass, f1710Width};

        if (testing)
        {
            // TF1 *f3bw = new TF1("f3bw", BW3, 1.1, 2.18, 9);
            // f3bw->SetParameter(0, 120000);
            // f3bw->SetParameter(1, f1270Mass - 0.06);
            // f3bw->SetParameter(2, f1270Width);
            // f3bw->SetParameter(3, 35000);
            // f3bw->SetParameter(4, a1320Mass + 0.03);
            // f3bw->SetParameter(5, a1320Width);
            // f3bw->SetParameter(6, 25000);
            // f3bw->SetParameter(7, f1525Mass);
            // f3bw->SetParameter(8, f1525Width);
            // f3bw->Draw("same");

            // TF1 *f1bw = new TF1("f1bw", BreitWignerpoly3, 1.3, 1.42, 7);
            // f1bw->SetParameter(0, 5000);
            // f1bw->SetParLimits(0, 0, 1e7);
            // f1bw->SetParameter(1, 1.355);
            // f1bw->SetParLimits(1, 1.3, 1.4);
            // f1bw->SetParameter(2, 0.2);
            // f1bw->SetParLimits(2, 0.1, 0.3);
            // hinvMass->Fit("f1bw", "REBMS");
            // f1bw->Draw("same");

            TF1 *f1bw = new TF1("f1bw", RelativisticBW, 1.34, 1.38, 3);
            f1bw->SetParameter(0, 2e4);
            f1bw->SetParameter(1, 1.355);
            f1bw->SetParameter(2, 0.1);
            hinvMass->Fit("f1bw", "REBMS");
            f1bw->Draw("same");
        }

        if (!testing)
        {
            // Fitting *********************************************
            TF1 *f3pol3;
            if (kchannel == "KsKs_Channel" && kResBkg == "MIX" && kbgfitfunction == "pol3")
            {
                f3pol3 = BW3pol3(hinvMass, parameters1, 1.1, 2.18);
                f3pol3->SetParameter(0, parameters1[0]);
                f3pol3->SetParameter(1, parameters1[1]);
                f3pol3->SetParLimits(1, parameters1[1] - 0.08, parameters1[1] + 0.08);
                f3pol3->SetParameter(2, parameters1[2]);
                f3pol3->SetParLimits(2, parameters1[2] - 0.01, parameters1[2] + 0.01);
                f3pol3->SetParameter(3, parameters1[3]);
                f3pol3->SetParameter(4, parameters1[4]);
                f3pol3->SetParameter(5, parameters1[5]);
                f3pol3->SetParameter(6, parameters1[6]);
                f3pol3->SetParameter(7, parameters1[7]);
                f3pol3->SetParLimits(7, parameters1[7] - 0.08, parameters1[7] + 0.08);
                f3pol3->FixParameter(8, parameters1[8]);
            }

            if (kchannel == "KsKs_Channel" && kResBkg == "MIX" && kbgfitfunction == "expol")
            {
                // FIXME: need to modify the parameters here for the successfull fit
                f3pol3 = BW3expo(hinvMass, parameters1, 1.1, 2.18);
                f3pol3->SetParameter(0, parameters1[0]);
                f3pol3->SetParameter(1, parameters1[1]);
                f3pol3->SetParLimits(1, parameters1[1] - 0.08, parameters1[1] + 0.08);
                f3pol3->SetParameter(2, parameters1[2]);
                f3pol3->SetParLimits(2, parameters1[2] - 0.01, parameters1[2] + 0.01);
                f3pol3->SetParameter(3, parameters1[3]);
                f3pol3->SetParameter(4, parameters1[4]);
                f3pol3->SetParameter(5, parameters1[5]);
                f3pol3->SetParameter(6, parameters1[6]);
                f3pol3->SetParameter(7, parameters1[7]);
                f3pol3->SetParLimits(7, parameters1[7] - 0.08, parameters1[7] + 0.08);
                f3pol3->FixParameter(8, parameters1[8]);
            }

            if (kchannel == "KsKs_Channel" && kResBkg == "ROTATED")
            {
                f3pol3 = BW3pol3(hinvMass, parameters1, 1.1, 2.18);
                f3pol3->SetParameter(0, parameters1[0]);
                f3pol3->SetParLimits(0, 0, 1e6);
                f3pol3->SetParameter(1, parameters1[1]);
                f3pol3->SetParLimits(1, parameters1[1] - 0.09, parameters1[1] + 0.09);
                f3pol3->SetParameter(2, parameters1[2]);
                f3pol3->SetParLimits(2, parameters1[2] - 0.05, parameters1[2] + 0.05);
                f3pol3->SetParameter(3, parameters1[3]);
                f3pol3->SetParameter(4, parameters1[4]);
                f3pol3->SetParLimits(4, parameters1[4] - 0.1, parameters1[4] + 0.1);
                f3pol3->SetParameter(5, parameters1[5]);
                f3pol3->SetParLimits(5, parameters1[5] - 0.01, parameters1[5] + 0.01);
                f3pol3->SetParameter(6, parameters1[6]);
                f3pol3->SetParameter(7, parameters1[7]);
                f3pol3->SetParLimits(7, parameters1[7] - 0.08, parameters1[7] + 0.08);
                f3pol3->FixParameter(8, parameters1[8]);
            }
            if (kchannel == "KK_Channel")
            {
                // f3pol3 = BW3pol3(hinvMass, parameters1, 1.1, 1.95);
                f3pol3 = new TF1("f3pol3", pol3bW3, 1.15, 1.7, 13);
                // hinvMass->SetMaximum(8e6);
                f3pol3->SetParameter(0, 3.666e4);
                f3pol3->SetParLimits(0, 0, 1e7);
                f3pol3->SetParameter(1, f1270Mass);
                // f3pol3->SetParLimits(1, 1.15, 1.27);
                f3pol3->SetParameter(2, f1270Width);
                f3pol3->SetParLimits(2, 0.05, 0.15);
                f3pol3->SetParameter(3, 2e4);
                f3pol3->SetParLimits(3, 0, 1e7);
                f3pol3->SetParameter(4, a1320Mass);
                // f3pol3->SetParLimits(4, 1.3, 1.4);
                f3pol3->SetParameter(5, a1320Width - 0.05);
                f3pol3->SetParameter(6, 35000);
                f3pol3->SetParLimits(6, 0, 1e7);
                f3pol3->FixParameter(7, f1525Mass);
                f3pol3->SetParLimits(7, 1.5, 1.58);
                f3pol3->SetParameter(8, f1525Width);
            }
            f3pol3->SetLineStyle(1);
            f3pol3->SetLineWidth(2);
            hinvMass->Fit("f3pol3", "REBMS0");
            f3pol3->Draw("same");

            // making the size of textbox optimal
            gPad->Update();
            TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
            ptstats->SetX1NDC(0.6);
            ptstats->SetX2NDC(0.99);
            ptstats->SetY1NDC(0.3);
            ptstats->SetY2NDC(0.95);
            ptstats->Draw("same");

            // // //Drawing the individual rBW functions and the legend
            TLegend *lfit = new TLegend(0.25, 0.6, 0.55, 0.9);
            lfit->SetFillColor(0);
            // lfit->SetBorderSize(0);
            lfit->SetFillStyle(0);
            lfit->SetTextFont(42);
            lfit->SetTextSize(0.035);
            lfit->AddEntry(hinvMass, "Data", "lpe");
            lfit->AddEntry(f3pol3, "3rBW + pol3", "l");
            double *parameters2 = f3pol3->GetParameters();
            TF1 *residualpol3 = draw_individual_functions(parameters2, lfit);
            lfit->Draw("same");
            t2->SetNDC();
            t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", 0.0, 30.0));
            if (saveplots)
            {
                string channel_temp = (kchannel == "KsKs_Channel") ? "KsKs_fit" : "KK_fit";
                c1->SaveAs(("saved/" + channel_temp + "/BW3pol3_" + kResBkg + Form("_%.1f_%.1f.png", pT_bins[ip], pT_bins[ip + 1])).c_str());
            }

            // ************* subtracting residual background *************
            TCanvas *c2 = new TCanvas("", "", 720, 720);
            SetCanvasStyle(c2, 0.12, 0.03, 0.05, 0.14);
            hinvMassResSub->Add(residualpol3, -1);
            hinvMassResSub->GetXaxis()->SetRangeUser(1.1, 2.18);
            hinvMassResSub->Draw();
            TF1 *f3bw3 = new TF1("f3bw3", BW3, 1.1, 2.18, 9);
            f3bw3->SetParameter(0, f3pol3->GetParameter(0));
            f3bw3->SetParameter(1, f3pol3->GetParameter(1));
            f3bw3->SetParameter(2, f3pol3->GetParameter(2));
            f3bw3->SetParameter(3, f3pol3->GetParameter(3));
            f3bw3->SetParameter(4, f3pol3->GetParameter(4));
            f3bw3->SetParameter(5, f3pol3->GetParameter(5));
            f3bw3->SetParameter(6, f3pol3->GetParameter(6));
            f3bw3->SetParameter(7, f3pol3->GetParameter(7));
            f3bw3->SetParameter(8, f3pol3->GetParameter(8));
            hinvMassResSub->Fit("f3bw3", "REBMS");
            hinvMassResSub->SetMarkerSize(0.8);
            hinvMassResSub->SetMaximum(hinvMassResSub->GetMaximum() * 1.9);

            TLegend *lfit2 = new TLegend(0.15, 0.72, 0.55, 0.89);
            lfit2->SetFillColor(0);
            lfit2->SetFillStyle(0);
            lfit2->SetTextFont(42);
            lfit2->SetBorderSize(0);
            lfit2->SetTextSize(0.035);
            lfit2->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
            lfit2->AddEntry(hinvMassResSub, "KsKs invariant mass", "lpe");
            lfit2->AddEntry(f3bw3, "3rBW fit", "l");
            lfit2->Draw("same");
            if (saveplots)
            {
                string channel_temp = (kchannel == "KsKs_Channel") ? "KsKs_fit" : "KK_fit";
                c2->SaveAs(("saved/" + channel_temp + "/BW3pol3_residual_sub_" + kResBkg + Form("_%.1f_%.1f.png", pT_bins[ip], pT_bins[ip + 1])).c_str());
            }
        }
    }
} //*******************************end of main function ***************************************

//**********************************all functions are defined here **********************************
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

    // // fit parameters here
    // f3pol3->SetParameter(0, parameters[0]);
    // f3pol3->SetParameter(1, parameters[1]);
    // f3pol3->SetParLimits(1, parameters[1] - 0.08, parameters[1] + 0.08);
    // f3pol3->SetParameter(2, parameters[2]);
    // f3pol3->SetParLimits(2, parameters[2] - 0.01, parameters[2] + 0.01);
    // f3pol3->SetParameter(3, parameters[3]);
    // f3pol3->SetParameter(4, parameters[4]);
    // f3pol3->SetParameter(5, parameters[5]);
    // f3pol3->SetParameter(6, parameters[6]);
    // f3pol3->SetParameter(7, parameters[7]);
    // f3pol3->SetParLimits(7, parameters[7] - 0.08, parameters[7] + 0.08);
    // f3pol3->FixParameter(8, parameters[8]);
    // f3pol3->SetLineStyle(1);
    // f3pol3->SetLineWidth(2);
    // h->Fit("f3pol3", "REBMS0");

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

TF1 *draw_individual_functions(double *parameters, TLegend *lfit)
{
    // BW1
    TF1 *frBW1 = new TF1("frBW1", RelativisticBW, 1.1, 2.2, 3);
    TF1 *frBW2 = new TF1("frBW2", RelativisticBW, 1.1, 2.2, 3);
    TF1 *frBW3 = new TF1("frBW3", RelativisticBW, 1.1, 2.2, 3);
    frBW1->SetParameter(0, parameters[0]);
    frBW1->SetParameter(1, parameters[1]);
    frBW1->SetParameter(2, parameters[2]);
    frBW1->SetLineStyle(2);
    frBW1->SetLineColor(3);
    frBW1->Draw("same");

    // BW2
    frBW2->SetParameter(0, parameters[3]);
    frBW2->SetParameter(1, parameters[4]);
    frBW2->SetParameter(2, parameters[5]);
    frBW2->SetLineStyle(2);
    frBW2->SetLineColor(4);
    frBW2->Draw("same");

    // BW3
    frBW3->SetParameter(0, parameters[6]);
    frBW3->SetParameter(1, parameters[7]);
    frBW3->SetParameter(2, parameters[8]);
    frBW3->SetLineStyle(2);
    frBW3->SetLineColor(6);
    frBW3->Draw("same");

    // pol3 background
    TF1 *fpol3 = new TF1("fpol3", polynomial3, 1.1, 2.2, 4);
    fpol3->SetParameter(0, parameters[9]);
    fpol3->SetParameter(1, parameters[10]);
    fpol3->SetParameter(2, parameters[11]);
    fpol3->SetParameter(3, parameters[12]);
    fpol3->SetLineStyle(2);
    fpol3->SetLineColor(28);
    fpol3->Draw("same");

    lfit->AddEntry(frBW1, "rBW(f_{2}(1270))", "l");
    lfit->AddEntry(frBW2, "rBW(f_{2}(1525))", "l");
    lfit->AddEntry(frBW3, "rBW(f_{0}(1710))", "l");
    lfit->AddEntry(fpol3, "pol3", "l");

    return fpol3;
}
