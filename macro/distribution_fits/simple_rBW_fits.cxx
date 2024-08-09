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
TF1 *draw_individual_functions(TF1 *fit, double *parameters, TLegend *lfit, bool drawpol2 = true);

void simple_rBW_fits()
{
    // // *********************** constant parameters *****************************
    const string kResBkg = "MIX";
    // const string kResBkg = "ROTATED";
    // const string kResBkg = "LIKE";
    const string kbgfitfunction = "pol3";
    // const string kbgfitfunction = "expol";
    // const string kbgfitfunction = "Boltzman";

    const int rebin = 1;
    bool saveplots = true;
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);

    const string outputfolder_str = "../" + kSignalOutput + "/" + kchannel + "/" + kfoldername;
    TString fits_folder = outputfolder_str + "/fits";
    const string fits_folder_str = outputfolder_str + "/fits/";
    float me_var_low[] = {1.10, 1.90, 2.00, 2.10, 2.20, 2.30};
    float me_var_high[] = {1.15, 2.00, 2.10, 2.20, 2.30, 2.50};
    TCanvas *cmevar = new TCanvas("", "", 1200, 720);
    SetCanvasStyle(cmevar, 0.12, 0.03, 0.05, 0.14);
    cmevar->Divide(3, 2);

    for (int ime = 0; ime < 6; ime++)
    {

        for (int ip = 0; ip < Npt; ip++)
        {

            TFile *f = new TFile((outputfolder_str + "/hglue_" + kResBkg + Form("_%.1f_%.1f_norm_%.2f_%.2f_.root", pT_bins[ip], pT_bins[ip + 1], me_var_low[ime], me_var_high[ime])).c_str(), "READ");

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

            // TCanvas *c1 = new TCanvas("", "", 720, 720);
            // SetCanvasStyle(c1, 0.12, 0.03, 0.05, 0.14);
            // cmevar->Close();
            cmevar->cd(ime + 1);
            hinvMass->GetYaxis()->SetTitleOffset(1.1);
            hinvMass->GetXaxis()->SetTitleOffset(1.3);
            hinvMass->GetXaxis()->SetRangeUser(1.0, 2.4);
            TH1F *hinvMassResSub = (TH1F *)hinvMass->Clone();
            // hinvMass->Draw();
            double parameters1[9] = {50, f1270Mass, f1270Width, 25, f1525Mass, f1525Width, 25, f1710Mass, f1710Width};

            // Fitting *********************************************
            TF1 *f3pol3;
            struct FitParams
            {
                double low;
                double high;
                double param0_low_limit;
                double param1_limit;
                double param2_limit;
                double param3_limit;
                double param6_limit;
                double param7_limit;
            };

            std::vector<FitParams> bwfit_params_me = {
                // {1.1, 2.15, -1, 0.08, 0.01, -1, -1, 0.08}, // for testing
                {1.09, 2.18, 0, 0.08, 0.01, 0, 0, 0.06},    // norm range 1
                {1.09, 2.13, -1, 0.08, 0.01, -1, -1, 0.08}, // norm range 2
                {1.09, 2.15, -1, 0.08, 0.01, -1, -1, 0.08}, // norm range 3
                {1.09, 2.12, 1, 0.08, 0.01, -1, -1, 0.08},  // norm range 4
                {1.1, 2.15, 1, 0.08, 0.01, -1, -1, 0.08},   // norm range 5
                {1.11, 2.15, -1, 0.08, 0.01, -1, -1, 0.08}  // norm range 6
            };

            const auto &iter_bin = bwfit_params_me[ime]; // variation for the norm range

            f3pol3 = BW3pol3(hinvMass, parameters1, iter_bin.low, iter_bin.high);
            f3pol3->SetParameter(0, parameters1[0]);
            if (iter_bin.param0_low_limit != -1)
            {
                f3pol3->SetParLimits(0, iter_bin.param0_low_limit, 1e9);
            }
            f3pol3->SetParameter(1, parameters1[1]);
            if (iter_bin.param1_limit != -1)
            {
                f3pol3->SetParLimits(1, parameters1[1] - iter_bin.param1_limit, parameters1[1] + iter_bin.param1_limit);
            }
            f3pol3->SetParameter(2, parameters1[2]);
            if (iter_bin.param2_limit != -1)
            {
                f3pol3->SetParLimits(2, parameters1[2] - iter_bin.param2_limit, parameters1[2] + iter_bin.param2_limit);
            }
            f3pol3->SetParameter(3, parameters1[3]);
            if (iter_bin.param3_limit != -1)
            {
                f3pol3->SetParLimits(3, parameters1[3] - iter_bin.param3_limit, 1e9);
            }
            f3pol3->SetParameter(4, parameters1[4]);
            f3pol3->SetParameter(5, parameters1[5]);
            f3pol3->SetParameter(6, parameters1[6]);
            if (iter_bin.param6_limit != -1)
            {
                f3pol3->SetParLimits(6, parameters1[6] - iter_bin.param6_limit, 1e9);
            }
            f3pol3->SetParameter(7, parameters1[7]);
            if (iter_bin.param7_limit != -1)
            {
                f3pol3->SetParLimits(7, parameters1[7] - iter_bin.param7_limit, parameters1[7] + iter_bin.param7_limit);
            }
            f3pol3->FixParameter(8, parameters1[8]);

            f3pol3->SetLineStyle(1);
            f3pol3->SetLineWidth(2);
            hinvMass->Fit("f3pol3", "REBMS0");
            hinvMass->Draw();
            hinvMass->GetXaxis()->SetRangeUser(f3pol3->GetXmin(), f3pol3->GetXmax());
            hinvMass->SetMaximum(hinvMass->GetMaximum() * 2.0);
            f3pol3->Draw("same");

            // making the size of textbox optimal
            gPad->Update();
            TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
            ptstats->SetX1NDC(0.6);
            ptstats->SetX2NDC(0.99);
            ptstats->SetY1NDC(0.5);
            ptstats->SetY2NDC(0.95);
            ptstats->Draw("same");

            // // //Drawing the individual rBW functions and the legend
            TLegend *lfit = new TLegend(0.15, 0.65, 0.55, 0.9);
            lfit->SetFillColor(0);
            lfit->SetBorderSize(0);
            lfit->SetFillStyle(0);
            lfit->SetTextFont(42);
            lfit->SetTextSize(0.025);
            lfit->AddEntry((TObject *)0, Form("Norm range: %.2f - %.2f GeV/c^{2}", me_var_low[ime], me_var_low[ime]), "");
            lfit->AddEntry(hinvMass, "Data", "lpe");
            lfit->AddEntry(f3pol3, "3rBW + pol3", "l");
            double *parameters2 = f3pol3->GetParameters();
            TF1 *residualpol3 = draw_individual_functions(f3pol3, parameters2, lfit);
            lfit->Draw("same");
            t2->SetNDC();
            t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip], pT_bins[ip + 1]));
            // if (saveplots)
            // {
            //     c1->SaveAs((fits_folder_str + "_" + kResBkg + Form("_%.1f_%.1f_me_range_compare_%.2f_%.2f_.png", pT_bins[ip], pT_bins[ip + 1], me_var_low[ime], me_var_high[ime])).c_str());
            // }

            // ************* subtracting residual background *************
            TCanvas *c2 = new TCanvas("", "", 2432, 85, 720, 720);
            SetCanvasStyle(c2, 0.12, 0.03, 0.05, 0.14);
            hinvMassResSub->Add(residualpol3, -1);
            float rangelow = f3pol3->GetXmin();
            float rangehigh = f3pol3->GetXmax();
            hinvMassResSub->GetXaxis()->SetRangeUser(rangelow + 0.01, 2.18);
            hinvMassResSub->Draw();
            TF1 *f3bw3 = new TF1("f3bw3", BW3, rangelow, rangehigh, 9);
            f3bw3->SetParNames("Norm", "Mass_{f1270}", "#Gamma_{f1270}", "Norm_{f1525}", "Mass_{f1525}", "#Gamma_{f1525}", "Norm_{f1710}", "Mass_{f1710}", "#Gamma_{f1710}");
            f3bw3->SetParameter(0, f3pol3->GetParameter(0));
            f3bw3->SetParameter(1, f3pol3->GetParameter(1));
            f3bw3->FixParameter(2, f3pol3->GetParameter(2));
            // f3bw3->SetParLimits(2, f3pol3->GetParameter(2) - 0.02, f3pol3->GetParameter(2) + 0.02);
            f3bw3->SetParameter(3, f3pol3->GetParameter(3));
            f3bw3->SetParameter(4, f3pol3->GetParameter(4));
            f3bw3->SetParameter(5, f3pol3->GetParameter(5));
            f3bw3->SetParameter(6, f3pol3->GetParameter(6));
            f3bw3->SetParameter(7, f3pol3->GetParameter(7));
            f3bw3->SetParameter(8, f3pol3->GetParameter(8));
            hinvMassResSub->Fit("f3bw3", "REBMS");
            hinvMassResSub->SetMarkerSize(0.8);
            hinvMassResSub->SetMaximum(hinvMassResSub->GetMaximum() * 1.9);

            gPad->Update();
            TPaveStats *ptstats2 = (TPaveStats *)hinvMassResSub->FindObject("stats");
            ptstats2->SetX1NDC(0.6);
            ptstats2->SetX2NDC(0.99);
            ptstats2->SetY1NDC(0.55);
            ptstats2->SetY2NDC(0.95);
            ptstats2->Draw("same");

            TLegend *lfit2 = new TLegend(0.17, 0.67, 0.57, 0.92);
            lfit2->SetFillColor(0);
            lfit2->SetFillStyle(0);
            lfit2->SetTextFont(42);
            lfit2->SetBorderSize(0);
            lfit2->SetTextSize(0.03);
            lfit2->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
            lfit2->AddEntry(hinvMassResSub, "KsKs invariant mass", "lpe");
            lfit2->AddEntry(f3bw3, "3rBW fit", "l");
            lfit2->Draw("same");
            draw_individual_functions(f3pol3, parameters2, lfit2, false); // draw the individual functions

            t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip], pT_bins[ip + 1]));
            // if (saveplots)
            // {
            //     c2->SaveAs((fits_folder_str + "_" + kResBkg + "res_sub_" + Form("_%.1f_%.1f.png", pT_bins[ip], pT_bins[ip + 1])).c_str());
            // }
        } // end of loop over pT bins
    } // end of me variations
    cmevar->SaveAs((fits_folder_str + "_" + kResBkg + "_me_range_compare.png").c_str());

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

TF1 *draw_individual_functions(TF1 *fit, double *parameters, TLegend *lfit, bool drawpol2 = true)
{
    // BW1
    float rangelow = fit->GetXmin();
    float rangehigh = fit->GetXmax();
    TF1 *frBW1 = new TF1("frBW1", RelativisticBW, rangelow, rangehigh, 3);
    TF1 *frBW2 = new TF1("frBW2", RelativisticBW, rangelow, rangehigh, 3);
    TF1 *frBW3 = new TF1("frBW3", RelativisticBW, rangelow, rangehigh, 3);
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
    TF1 *fpol3 = new TF1("fpol3", polynomial3, rangelow, rangehigh, 4);
    fpol3->SetParameter(0, parameters[9]);
    fpol3->SetParameter(1, parameters[10]);
    fpol3->SetParameter(2, parameters[11]);
    fpol3->SetParameter(3, parameters[12]);
    fpol3->SetLineStyle(2);
    fpol3->SetLineColor(28);
    if (drawpol2)
    {
        fpol3->Draw("same");
        lfit->AddEntry(fpol3, "pol3", "l");
    }
    lfit->AddEntry(frBW1, "rBW(f_{2}(1270))", "l");
    lfit->AddEntry(frBW2, "rBW(f_{2}(1525))", "l");
    lfit->AddEntry(frBW3, "rBW(f_{0}(1710))", "l");

    return fpol3;
}
