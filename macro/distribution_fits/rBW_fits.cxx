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
    bool saveplots = true;
    // double f1710Mass = pdg->GetParticle(10331)->Mass();
    // double f1710Width = pdg->GetParticle(10331)->Width();
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);

    const string outputfolder_str = "../" + kSignalOutput + "/" + kchannel + "/" + kfoldername;
    TString fits_folder = outputfolder_str + "/fits";
    const string fits_folder_str = outputfolder_str + "/fits/";
    // Create the folder to save the fitted distributions
    if (gSystem->mkdir(fits_folder, kTRUE))
    {
        std::cout << "Folder " << fits_folder << " created successfully." << std::endl;
    }
    else
    {
        std::cout << "Creating folder " << fits_folder << std::endl;
    }

    // // *********************** constant parameters *****************************
    TGraphErrors *massf1270 = new TGraphErrors(Npt);
    TGraphErrors *widthf1270 = new TGraphErrors(Npt);
    TGraphErrors *massf1525 = new TGraphErrors(Npt);
    TGraphErrors *widthf1525 = new TGraphErrors(Npt);
    TGraphErrors *massf1710 = new TGraphErrors(Npt);
    TGraphErrors *widthf1710 = new TGraphErrors(Npt);

    // pT loop ***************************************************
    for (Int_t ip = pt_start; ip < pt_end; ip++)
    {
        TFile *f = new TFile((outputfolder_str + "/hglue_" + kResBkg + Form("_%.1f_%.1f.root", pT_bins[ip], pT_bins[ip + 1])).c_str(), "READ");
        // TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC220_pass6_small/230281/KsKs_Channel/strangeness_tutorial/hglue_MIX_0.0_30.0_norm_1.9_2.0_.root", "READ");

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
                struct FitParams
                {
                    double low;
                    double high;
                    double param0_low_limit;
                    double param1_limit;
                };

                // Define the fit parameters for each pT bin
                std::vector<FitParams> bwfit_params_me = {
                    {1.1, 2.15, -1, 0.08}, // for testing purpose for single 
                    // {1.1, 2.15, -1, 0.08}, // for full pT range
                    {1.09, 2.18, 3, 0.09}, // pT 1 to 2
                    {1.11, 2.16, 0, 0.09}, // pT 2 to 3
                    {1.12, 1.95, 0, 0.09}, // pT 3 to 4
                    {1.1, 2.15, -1, 0.08}, // pT 4 to 6
                    {1.1, 2.15, -1, 0.08}  // pT 6 to 12
                };

                const auto &iter_bin = bwfit_params_me[ip];

                // // for all pT bins
                f3pol3 = BW3pol3(hinvMass, parameters1, iter_bin.low, iter_bin.high);
                f3pol3->SetParameter(0, parameters1[0]);
                if (iter_bin.param0_low_limit != -1)
                {
                    f3pol3->SetParLimits(0, iter_bin.param0_low_limit, 1e9);
                }
                f3pol3->SetParameter(1, parameters1[1]);
                f3pol3->SetParLimits(1, parameters1[1] - iter_bin.param1_limit, parameters1[1] + iter_bin.param1_limit);
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

            if (kchannel == "KsKs_Channel" && kResBkg == "ROTATED" && kbgfitfunction == "pol3")
            {
                struct FitParams2
                {
                    double low;
                    double high;
                    double param0_low_limit; // norm for f1270
                    double param1_limit;     // mass for f1270
                    double param2_limit;     // width for f1270
                    double param4_limit;     // mass for f1525
                    double param5_limit;     // width for f1525
                    double param7_limit;     // mass for f1710
                };
                // Define the fit parameters for each pT bin
                std::vector<FitParams2> bwfit_params_rot = {
                    // {1.05, 2.15, 0, 0.1, 0.02, 0.1, 0.01, 0.08}, // for testing purpose for single pT bin
                    {1.11, 2.15, 1, 0.07, -1.0, 0.1, 0.01, 0.08}, // pT 1 to 2
                    {1.11, 2.15, 1, 0.07, -1.0, 0.1, 0.01, 0.08}, // pT 2 to 3
                    {1.11, 2.15, 1, 0.07, -1.0, 0.1, 0.01, 0.08}, // pT 3 to 4
                    {1.09, 2.15, 0, 0.07, 0.05, 0.1, 0.01, 0.08}, // pT 4 to 6
                    {1.05, 2.15, 0, 0.1, 0.02, 0.1, 0.01, 0.08}   // pT 6 to 12 and min bias
                };

                const auto &iter_bin = bwfit_params_rot[ip];

                f3pol3 = BW3pol3(hinvMass, parameters1, iter_bin.low, iter_bin.high);
                f3pol3->SetParameter(0, parameters1[0]);
                f3pol3->SetParLimits(0, iter_bin.param0_low_limit, 1e9);
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
                f3pol3->SetParameter(4, parameters1[4]);
                if (iter_bin.param4_limit != -1)
                {
                    f3pol3->SetParLimits(4, parameters1[4] - iter_bin.param4_limit, parameters1[4] + iter_bin.param4_limit);
                }
                f3pol3->SetParameter(5, parameters1[5]);
                if (iter_bin.param5_limit != -1)
                {
                    f3pol3->SetParLimits(5, parameters1[5] - iter_bin.param5_limit, parameters1[5] + iter_bin.param5_limit);
                }
                f3pol3->SetParameter(6, parameters1[6]);
                f3pol3->SetParameter(7, parameters1[7]);
                if (iter_bin.param7_limit != -1)
                {
                    f3pol3->SetParLimits(7, parameters1[7] - iter_bin.param7_limit, parameters1[7] + iter_bin.param7_limit);
                }
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
            hinvMass->SetMaximum(hinvMass->GetMaximum() * 1.4);
            hinvMass->Fit("f3pol3", "REBMS0");
            f3pol3->Draw("same");

            // Setting the values of the mass and width of the resonances into Tgraphs
            massf1270->SetPoint(ip, (pT_bins[ip] + pT_bins[ip + 1]) / 2, f3pol3->GetParameter(1));
            massf1270->SetPointError(ip, (pT_bins[ip + 1] - pT_bins[ip]) / 2.0, f3pol3->GetParError(1));
            widthf1270->SetPoint(ip, (pT_bins[ip] + pT_bins[ip + 1]) / 2, f3pol3->GetParameter(2));
            widthf1270->SetPointError(ip, (pT_bins[ip + 1] - pT_bins[ip]) / 2.0, f3pol3->GetParError(2));
            massf1525->SetPoint(ip, (pT_bins[ip] + pT_bins[ip + 1]) / 2, f3pol3->GetParameter(4));
            massf1525->SetPointError(ip, (pT_bins[ip + 1] - pT_bins[ip]) / 2.0, f3pol3->GetParError(4));
            widthf1525->SetPoint(ip, (pT_bins[ip] + pT_bins[ip + 1]) / 2, f3pol3->GetParameter(5));
            widthf1525->SetPointError(ip, (pT_bins[ip + 1] - pT_bins[ip]) / 2.0, f3pol3->GetParError(5));
            massf1710->SetPoint(ip, (pT_bins[ip] + pT_bins[ip + 1]) / 2, f3pol3->GetParameter(7));
            massf1710->SetPointError(ip, (pT_bins[ip + 1] - pT_bins[ip]) / 2.0, f3pol3->GetParError(7));

            // making the size of textbox optimal
            gPad->Update();
            TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
            ptstats->SetX1NDC(0.6);
            ptstats->SetX2NDC(0.99);
            ptstats->SetY1NDC(0.5);
            ptstats->SetY2NDC(0.95);
            ptstats->Draw("same");

            // // //Drawing the individual rBW functions and the legend
            TLegend *lfit = new TLegend(0.25, 0.72, 0.55, 0.92);
            lfit->SetFillColor(0);
            // lfit->SetBorderSize(0);
            lfit->SetFillStyle(0);
            lfit->SetTextFont(42);
            lfit->SetTextSize(0.025);
            lfit->AddEntry(hinvMass, "Data", "lpe");
            lfit->AddEntry(f3pol3, "3rBW + pol3", "l");
            double *parameters2 = f3pol3->GetParameters();
            TF1 *residualpol3 = draw_individual_functions(f3pol3, parameters2, lfit);
            lfit->Draw("same");
            t2->SetNDC();
            t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip], pT_bins[ip + 1]));
            if (saveplots)
            {
                c1->SaveAs((fits_folder_str + "_" + kResBkg + Form("_%.1f_%.1f.png", pT_bins[ip], pT_bins[ip + 1])).c_str());
            }

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
            if (saveplots)
            {
                c2->SaveAs((fits_folder_str + "_" + kResBkg + "res_sub_" + Form("_%.1f_%.1f.png", pT_bins[ip], pT_bins[ip + 1])).c_str());
            }
        }

    } //************************ end of pT loop **************************************

    if (Npt > 1)
    {
        TFile *file_fitparams = new TFile((fits_folder_str + "fitparams_" + kResBkg + ".root").c_str(), "RECREATE");
        // ************* Drawing the mass and width graphs *************
        TCanvas *c3 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.14);
        SetGrapherrorStyle(massf1270);
        massf1270->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
        massf1270->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        massf1270->SetMaximum(f1270Mass + 0.15);
        massf1270->SetMinimum(f1270Mass - 0.1);
        massf1270->Draw("ape");
        massf1270->Write("massf1270");
        TLine *linepdg = new TLine(0, f1270Mass, 12, f1270Mass);
        linepdg->SetLineColor(kRed);
        linepdg->SetLineStyle(2);
        linepdg->SetLineWidth(2);
        linepdg->Draw("same");
        TLegend *lfit3 = new TLegend(0.17, 0.75, 0.47, 0.92);
        SetLegendStyle(lfit3);
        lfit3->SetTextSize(0.035);
        lfit3->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
        lfit3->AddEntry(massf1270, "BW fit K_{s}K_{s} inv. mass", "lpe");
        lfit3->AddEntry(linepdg, "f1270 PDG mass", "l");
        lfit3->Draw("same");
        if (saveplots)
        {
            c3->SaveAs((fits_folder_str + "_" + kResBkg + "massf1270.png").c_str());
        }

        TCanvas *c4 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c4, 0.15, 0.03, 0.05, 0.14);
        SetGrapherrorStyle(widthf1270);
        widthf1270->GetYaxis()->SetTitle("Width (GeV/c^{2})");
        widthf1270->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        widthf1270->SetMaximum(f1270Width + 0.3);
        widthf1270->SetMinimum(f1270Width - 0.3);
        widthf1270->Draw("ape");
        widthf1270->Write("widthf1270");
        linepdg->SetY1(f1270Width);
        linepdg->SetY2(f1270Width);
        linepdg->Draw("same");
        lfit3->Clear();
        lfit3->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
        lfit3->AddEntry(widthf1270, "BW fit K_{s}K_{s} inv. mass", "lpe");
        lfit3->AddEntry(linepdg, "f1270 PDG width", "l");
        lfit3->Draw("same");
        if (saveplots)
            c4->SaveAs((fits_folder_str + "_" + kResBkg + "widthf1270.png").c_str());

        TCanvas *c5 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c5, 0.15, 0.03, 0.05, 0.14);
        SetGrapherrorStyle(massf1525);
        massf1525->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
        massf1525->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        massf1525->SetMaximum(f1525Mass + 0.15);
        massf1525->SetMinimum(f1525Mass - 0.1);
        massf1525->Draw("ape");
        massf1525->Write("massf1525");
        linepdg->SetY1(f1525Mass);
        linepdg->SetY2(f1525Mass);
        linepdg->Draw("same");
        lfit3->Clear();
        lfit3->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
        lfit3->AddEntry(massf1525, "BW fit K_{s}K_{s} inv. mass", "lpe");
        lfit3->AddEntry(linepdg, "f1525 PDG mass", "l");
        lfit3->Draw("same");
        if (saveplots)
            c5->SaveAs((fits_folder_str + "_" + kResBkg + "massf1525.png").c_str());

        TCanvas *c6 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c6, 0.15, 0.03, 0.05, 0.14);
        SetGrapherrorStyle(widthf1525);
        widthf1525->GetYaxis()->SetTitle("Width (GeV/c^{2})");
        widthf1525->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        widthf1525->SetMaximum(f1525Width + 0.15);
        widthf1525->SetMinimum(f1525Width - 0.1);
        widthf1525->Draw("ape");
        widthf1525->Write("widthf1525");
        linepdg->SetY1(f1525Width);
        linepdg->SetY2(f1525Width);
        linepdg->Draw("same");
        lfit3->Clear();
        lfit3->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
        lfit3->AddEntry(widthf1525, "BW fit K_{s}K_{s} inv. mass", "lpe");
        lfit3->AddEntry(linepdg, "f1525 PDG width", "l");
        lfit3->Draw("same");
        if (saveplots)
            c6->SaveAs((fits_folder_str + "_" + kResBkg + "widthf1525.png").c_str());

        TCanvas *c7 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c7, 0.15, 0.03, 0.05, 0.14);
        SetGrapherrorStyle(massf1710);
        massf1710->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
        massf1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        massf1710->SetMaximum(f1710Mass + 0.15);
        massf1710->SetMinimum(f1710Mass - 0.1);
        massf1710->Draw("ape");
        massf1710->Write("massf1710");
        linepdg->SetY1(f1710Mass);
        linepdg->SetY2(f1710Mass);
        linepdg->Draw("same");
        lfit3->Clear();
        lfit3->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");

        lfit3->AddEntry(massf1710, "BW fit K_{s}K_{s} inv. mass", "lpe");
        lfit3->AddEntry(linepdg, "f1710 PDG mass", "l");
        lfit3->Draw("same");
        if (saveplots)
            c7->SaveAs((fits_folder_str + "_" + kResBkg + "massf1710.png").c_str());

        // Since the width of f1710 is fixed, so it the graph will be empty for it
        TCanvas *c8 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c8, 0.15, 0.03, 0.05, 0.14);
        SetGrapherrorStyle(widthf1710);
        widthf1710->GetYaxis()->SetTitle("Width (GeV/c^{2})");
        widthf1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        widthf1710->SetMaximum(f1710Width + 0.1);
        widthf1710->SetMinimum(f1710Width - 0.1);
        widthf1710->Draw("ape");
        widthf1710->Write("widthf1710");
        linepdg->SetY1(f1710Width);
        linepdg->SetY2(f1710Width);
        linepdg->Draw("same");
        lfit3->Clear();
        lfit3->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
        lfit3->AddEntry(widthf1710, "BW fit K_{s}K_{s} inv. mass", "lpe");
        lfit3->AddEntry(linepdg, "f1710 PDG width", "l");
        lfit3->Draw("same");
        if (saveplots)
            c8->SaveAs((fits_folder_str + "_" + kResBkg + "widthf1710.png").c_str());
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
