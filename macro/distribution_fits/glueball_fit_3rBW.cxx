#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

Double_t single_BW_hera(double *x, double *par);
Double_t single_BW(double *x, double *par);
Double_t BWsum_hera(double *x, double *par);
Double_t BWsum(double *x, double *par);

Double_t exponential_bkg_1(double *x, double *par); // 3 parameters
Double_t exponential_bkg_2(double *x, double *par); // 3 parameters
Double_t exponential_bkg_3(double *x, double *par); // 4 parameters
Double_t exponential_bkg_4(double *x, double *par); // 5 parameters
Double_t exponential_bkg_5(double *x, double *par); // 3 parameters
Double_t exponential_bkg_6(double *x, double *par); // 4 parameters
Double_t expol_chkstar(double *x, double *par);     // 4 parameters

Double_t Boltzmann_bkg_1(double *x, double *par); // 3 parameters
Double_t Boltzmann_bkg_2(double *x, double *par); // 4 parameters

Double_t single_BW_expol3(double *x, double *par);
Double_t single_BW_expol3_hera(double *x, double *par);
Double_t BWsum_expol3(double *x, double *par);
Double_t BWsum_expol3_hera(double *x, double *par);
Double_t BWsum_expol_chkstar(double *x, double *par);

Double_t single_BW_boltzman_1(double *x, double *par);
Double_t single_BW_boltzman_2(double *x, double *par);
Double_t BWsum_boltzman_1(double *x, double *par);
Double_t BWsum_boltzman_2(double *x, double *par);

void glueball_fit_3rBW()
{
    ofstream file;
    file.open("fit_params.txt");
    string savepath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/fits/3rBW_fits";
    // string savepath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/294059/KK_Channel/kaonkaonAnalysisRun3/fits/3rBW_fits";

    TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60_fullpt.root", "READ"); // full pT range
    // TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60.root", "READ");        // pT differential range

    // TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/294059/KK_Channel/kaonkaonAnalysisRun3/hglue_LIKE_norm_2.50_2.60_pt_0.0_30.0.root", "READ"); // KK channel

    int colors[] = {4, 6, 28, 46};
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);
    t2->SetTextFont(42);

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    // #define b_expol
// #define b_boltzman
#define b_expol1
    // #define residual_subtracted

    for (int ipt = 0; ipt < Npt; ipt++)
    {

        // TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", 0.0, 30.0));
        TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", pT_bins[ipt], pT_bins[ipt + 1]));

        if (hinvMass == nullptr)
        {
            cout << "Error opening histogram" << endl;
            return;
        }
        TCanvas *c = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
        hinvMass->Rebin(2);
        hinvMass->GetXaxis()->SetRangeUser(1.00, 2.50);
        hinvMass->GetXaxis()->SetTitle("M_{K^{0}_{s}K^{0}_{s}} (GeV/c^{2})");
        hinvMass->Draw();
        // t2->DrawLatex(0.29, 0.96, Form("%.1f < #it{p}_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));
        TH1F *hsubtracted_res = (TH1F *)hinvMass->Clone("hsubtracted_res");
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1111);
        vector<tuple<float, int, float, float>> fit_parameters;

// // //************************************************************************ */
// // // **************** For BW sum with expol HERA ****************************

// // Default fitting range is 1.02 to 2.20. Four types of fitting range variations: extend left (1.0), extend right (2.50), large range (1.0 to 2.50), small range (1.05 to 2.15)
#ifdef b_expol
        TF1 *BEexpol = new TF1("BEexpol", BWsum_expol3, 1.05, 2.30, 13); // expol 3
        string parnames[] = {"norm1270", "mass1270", "width1270", "norm1525", "mass1525", "width1525", "norm1710", "mass1710", "width1710", "expol1", "expol2", "expol3", "expol4"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        // double parameters[] = {10, f1270Mass, f1270Width, 10, f1525Mass, f1525Width, 5, f1710Mass, f1710Width};
        double parameters[] = {8384, f1270Mass, f1270Width, 7858, f1525Mass, f1525Width, 3218, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {4, 5 * f1525Width}, {7, 10 * f1710WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        // for rotational bkg with pt range 0-30 GeV/c
        BEexpol->SetParameter(size_fitparams + 0, 5.562e5);  // 5.562e5   // Fix
        BEexpol->SetParameter(size_fitparams + 1, -0.09379); // -0.09379  //Free
        BEexpol->SetParameter(size_fitparams + 2, 2.569);    // 2.569     // Fix
        BEexpol->SetParameter(size_fitparams + 3, 1.098);    // 1.098     // Free

        // // for rotational bkg with pt range 1-30 GeV/c
        // BEexpol->FixParameter(size_fitparams + 0, 5.927e5);  // 5.562e5   // Fix
        // BEexpol->SetParameter(size_fitparams + 1, -0.05466); // -0.09379  //Free
        // BEexpol->FixParameter(size_fitparams + 2, 3.26);     // 2.569     // Fix
        // BEexpol->SetParameter(size_fitparams + 3, 0.9221);   // 1.098     // Free

        // for ME bkg with pt range 1-30 GeV/c
        // Till now the problem with ME data is that, it goes negative. Now the expol fit does not go to neagtive values. Need to fix this.
        // BEexpol->SetParameter(size_fitparams + 0, -0); // 5.562e5   // Fix
        // BEexpol->SetParameter(size_fitparams + 1, 0.4);      // -0.09379  //Free
        // BEexpol->SetParameter(size_fitparams + 2, 15.5);     // 2.569     // Fix
        // BEexpol->SetParameter(size_fitparams + 3, 2.8);      // 1.098     // Free
        // // BEexpol->SetParLimits(6, 1500, 10000);
        // // BEexpol->SetParLimits(0, 0, 10000);
        // // BEexpol->SetParLimits(3, 0, 10000);

        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(2, f1270Width);
        // BEexpol->FixParameter(4, f1525Mass);
        // BEexpol->FixParameter(5, f1525Width);
        // BEexpol->FixParameter(7, f1710Mass);
        // BEexpol->FixParameter(8, f1710Width);

        // BEexpol->FixParameter(1, a1320Mass);
        // BEexpol->FixParameter(2, a1320Width);

        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        // status codes: 4000 successful, 4 call limit, 4910 failed
        string fitstatus = "Successfull";
        if (fitResultptr->Status() != 4000)
        {
            cout << "Fit failed or call limit" << endl;
            fitstatus = "Failed";
        }
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;

        double *obtained_parameters = BEexpol->GetParameters();
        TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
        TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
        for (int i = 0; i < 4; i++)
        {
            expol->SetParameter(i, obtained_parameters[size_fitparams + i]);
            expol_clone->SetParameter(i, obtained_parameters[size_fitparams + i]);
        }
        expol->SetLineColor(3);
        expol->SetLineStyle(2);
        expol_clone->SetLineColor(3);
        expol_clone->SetLineStyle(2);
        expol->Draw("same");

        TF1 *onlyBW = new TF1("onlyBW", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        onlyBW_clone->SetParNames("norm1270", "mass1270", "width1270", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710");
        for (int i = 0; i < 9; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        onlyBW->Draw("same");

        TLegend *ltemp = new TLegend(0.20, 0.67, 0.52, 0.92);
        ltemp->SetFillStyle(0);
        ltemp->SetTextFont(42);
        ltemp->SetTextSize(0.035);
        ltemp->AddEntry(hinvMass, "Data", "lpe");
        ltemp->AddEntry(BEexpol, "3rBW + expol", "l");
        ltemp->AddEntry(onlyBW, "3rBW", "l");
        ltemp->AddEntry(expol, "expol", "l");
        ltemp->Draw("same");
#endif

        // // // ************************************************************************************
        // // // **************** For BW sum with Boltzmann ****************************
#ifdef b_boltzman

        // // int iteration = 0;
        // // for (int ipar1 = 700000; ipar1 < 800000; ipar1 += 10000) // loop for expol parameter 1
        // // {
        // //     for (double ipar2 = 0.4; ipar2 <= 0.8; ipar2 += 0.02) // loop for expol parameter 2
        // //     {
        // //         for (double ipar3 = 3.8; ipar3 < 4.6; ipar3 += 0.1) // loop for expol parameter 3
        // //         {

        // Default fitting range is 1.02 to 2.20. Four types of fitting range variations: extend left (1.0), extend right (2.50), large range (1.0 to 2.50), small range (1.05 to 2.15)

        TF1 *BEexpol = new TF1("BEexpol", BWsum_boltzman_1, 1.05, 2.30, 12); // Boltzmann
        string parnames[] = {"norm1270", "mass1270", "width1270", "norm1525", "mass1525", "width1525", "norm1710", "mass1710", "width1710", "Boltzmann1", "Boltzmann2", "Boltzmann3", "Boltzmann4"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        // double parameters[] = {100, f1270Mass, f1270Width, 100, f1525Mass, f1525Width, 50, f1710Mass, f1710Width};
        double parameters[] = {6000, f1270Mass, f1270Width, 8000, f1525Mass, f1525Width, 4000, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {2, 10 * f1270WidthErr}, {4, 5 * f1525Width}, {7, 20 * f1710WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        BEexpol->SetParameter(size_fitparams + 0, 7.618e5); // expol 3  // 7.618e5  // 5.562e5
        BEexpol->SetParameter(size_fitparams + 1, 0.6456);  // expol 3  //  0.6456  // -0.09379
        BEexpol->SetParameter(size_fitparams + 2, 4.238);   // expol 3  //4.238  // 2.569

        // BEexpol->SetParameter(size_fitparams + 0, ipar1); // expol 3  // 710000  // 5.562e5
        // BEexpol->FixParameter(size_fitparams + 1, ipar2); // expol 3  // -0.03  // -0.09379
        // BEexpol->FixParameter(size_fitparams + 2, ipar3); // expol 3  // 2.78  // 2.569

        // BEexpol->FixParameter(0, 6998);
        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(2, f1270Width);
        // BEexpol->FixParameter(3, 7128);
        // BEexpol->FixParameter(4, f1525Mass);
        // BEexpol->FixParameter(5, f1525Width);
        // BEexpol->FixParameter(6, 4058);
        // BEexpol->FixParameter(7, f1710Mass);
        // BEexpol->FixParameter(8, f1710Width);

        // BEexpol->FixParameter(1, a1320Mass);
        // BEexpol->FixParameter(2, a1320Width);

        hinvMass->Fit("BEexpol", "REMBS");
        // //             float chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
        // //             fit_parameters.push_back(make_tuple(chi2ndf, ipar1, ipar2, ipar3));
        // //             iteration++;
        // //             cout << "Iteration: " << iteration << endl;
        // //         }
        // //     }
        // // }

        // // // sort in asceding order w.r.t to the third array i.e. chi2/NDF
        // // sort(fit_parameters.begin(), fit_parameters.end(),
        // //      [](const auto &a, const auto &b)
        // //      {
        // //          return get<0>(a) < get<0>(b);
        // //      });
        //
        // // for (int i = 0; i < 20; i++)
        // // {
        // //     float best_ipar1 = std::get<1>(fit_parameters[i]);
        // //     float best_ipar2 = std::get<2>(fit_parameters[i]);
        // //     float best_ipar3 = std::get<3>(fit_parameters[i]);
        // //     float best_chi2ndf = std::get<0>(fit_parameters[i]);
        // //     cout << "ipar1: " << best_ipar1 << ",  ipar2: " << best_ipar2 << ",  ipar3: " << best_ipar3 << ", chi2/NDF: " << best_chi2ndf << endl;
        // // }

        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        // status codes: 4000 successful, 4 call limit, 4910 failed
        string fitstatus = "Successfull";
        // cout<<"status code is "<<fitResultptr->Status()<<endl;
        if (fitResultptr->Status() != 4070)
        {
            cout << "Fit failed or call limit" << endl;
            fitstatus = "Failed";
        }
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;

        double *obtained_parameters = BEexpol->GetParameters();
        TF1 *expol = new TF1("expol", Boltzmann_bkg_1, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
        TF1 *expol_clone = new TF1("expol_clone", Boltzmann_bkg_1, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
        for (int i = 0; i < 4; i++)
        {
            expol->SetParameter(i, obtained_parameters[size_fitparams + i]);
            expol_clone->SetParameter(i, obtained_parameters[size_fitparams + i]);
        }
        expol->SetLineColor(3);
        expol->SetLineStyle(2);
        expol_clone->SetLineColor(3);
        expol_clone->SetLineStyle(2);
        expol->Draw("same");

        TF1 *onlyBW = new TF1("onlyBW", BWsum_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        onlyBW_clone->SetParNames("norm1270", "mass1270", "width1270", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710");
        for (int i = 0; i < 9; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        onlyBW->Draw("same");

        TLegend *ltemp = new TLegend(0.20, 0.67, 0.52, 0.92);
        ltemp->SetFillStyle(0);
        ltemp->SetTextFont(42);
        ltemp->SetTextSize(0.035);
        ltemp->AddEntry(hinvMass, "Data", "lpe");
        ltemp->AddEntry(BEexpol, "3rBW + Boltzmann", "l");
        ltemp->AddEntry(onlyBW, "3rBW", "l");
        ltemp->AddEntry(expol, "Boltzmann", "l");
        ltemp->Draw("same");
#endif

        // // // //************************************************************************ */
        // // // // **************** For BW sum with exp + pol2 as used in Charged kstar **************************
#ifdef b_expol1
        TF1 *BEexpol = new TF1("BEexpol", BWsum_expol_chkstar, 1.03, 2.20, 13); // expol 1 (charked star)
        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "A", "n", "b", "c"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        double parameters[] = {8384, f1270Mass, f1270Width, 7858, f1525Mass, f1525Width, 3218, f1710Mass, f1710Width};
        // double parameters[] = {8384, a1320Mass, a1320Width, 7858, f1525Mass, f1525Width, 3218, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }

        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {4, 5 * f1525Width}, {7, 20 * f1710WidthErr}};
        // vector<vector<float>> par_limits = {{1, 3 * a1320Mass}, {4, 5 * f1525Width}, {7, 20 * f1710WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        // for rotational bkg with pt range 0-30 GeV/c (KsKs channel)
        BEexpol->SetParameter(size_fitparams + 0, 3.86556e+02);  // 5.562e5   // Fix
        BEexpol->SetParameter(size_fitparams + 1, -5.10155e-02); // -0.09379  //Free
        BEexpol->SetParameter(size_fitparams + 2, 1.01604e+01);  // 2.569     // Fix
        BEexpol->SetParameter(size_fitparams + 3, -2.73381e+00); // 1.098     // Free

        // // for like sign bkg with pt range 0-30 GeV/c (KK channel)
        // BEexpol->SetParameter(size_fitparams + 0, 355.443); //
        // BEexpol->SetParameter(size_fitparams + 1, -0.388);  //
        // BEexpol->SetParameter(size_fitparams + 2, 10.107);  //
        // BEexpol->SetParameter(size_fitparams + 3, -1.697);  //

        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(2, f1270Width);
        // BEexpol->FixParameter(4, f1525Mass);
        // BEexpol->FixParameter(5, f1525Width);
        // BEexpol->FixParameter(7, f1710Mass);
        // BEexpol->FixParameter(8, f1710Width);

        // BEexpol->FixParameter(1, a1320Mass);
        // BEexpol->FixParameter(2, a1320Width);

        hinvMass->Fit("BEexpol", "REBMS");
        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        // status codes: 4000 successful, 4 call limit, 4910 failed
        string fitstatus = "Successfull";
        // cout<<"fit status code "<<fitResultptr->Status()<<endl;
        if (fitResultptr->Status() != 4140)
        {
            cout << "Fit failed or call limit !!!!!!!" << endl;
            fitstatus = "Failed";
        }
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;

        double *obtained_parameters = BEexpol->GetParameters();
        TF1 *expol = new TF1("expol", expol_chkstar, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
        TF1 *expol_clone = new TF1("expol_clone", expol_chkstar, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
        for (int i = 0; i < 4; i++)
        {
            expol->SetParameter(i, obtained_parameters[size_fitparams + i]);
            expol_clone->SetParameter(i, obtained_parameters[size_fitparams + i]);
        }
        expol->SetLineColor(3);
        expol->SetLineStyle(2);
        expol_clone->SetLineColor(3);
        expol_clone->SetLineStyle(2);
        expol->Draw("same");

        TF1 *onlyBW = new TF1("onlyBW", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        for (int i = 0; i < 9; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        // onlyBW->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits[4];
        for (int i = 0; i < 4; i++)
        {
            singlefits[i] = new TF1(Form("singlef%d", i), single_BW, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
            singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
            singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
            singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
            singlefits[i]->SetLineColor(colors[i]);
            singlefits[i]->SetLineStyle(2);
            singlefits[i]->SetLineWidth(2);
            singlefits[i]->Draw("same");
        }

        TLegend *ltemp = new TLegend(0.25, 0.55, 0.55, 0.87);
        ltemp->SetFillStyle(0);
        ltemp->SetBorderSize(0);
        ltemp->SetTextFont(42);
        ltemp->SetTextSize(0.03);
        ltemp->AddEntry((TObject *)0, "", "");
        ltemp->AddEntry((TObject *)0, "", "");
        ltemp->AddEntry(hinvMass, "Data (stat. uncert.)", "lpe");
        ltemp->AddEntry(BEexpol, "3rBW + Residual BG", "l");
        // ltemp->AddEntry(onlyBW, "4rBW", "l");
        ltemp->AddEntry(expol, "Residual BG", "l");
        ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
        ltemp->AddEntry(singlefits[2], "f'_{2}(1525)", "l");
        ltemp->AddEntry(singlefits[3], "f_{0}(1710)", "l");
        ltemp->Draw("same");

        TLatex lat1;
        lat1.SetNDC();
        lat1.SetTextSize(0.03);
        lat1.SetTextFont(42);
        lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
        lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
        lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));
#endif

        // // //******************************************************************************************
        // // //********************************* common for all fits ***************************************
        // // //****************************common for all fits************************************
        // // //************************** printing fit parameters in the file *******************************
        gPad->Update();
        TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
        ptstats->SetX1NDC(0.6);
        ptstats->SetX2NDC(0.99);
        ptstats->SetY1NDC(0.4);
        ptstats->SetY2NDC(0.92);
        ptstats->Draw("same");
        c->SaveAs((savepath + Form("/rBWfit_pt_%.2f_%.2f.png", pT_bins[ipt], pT_bins[ipt + 1])).c_str());

        double chi2_ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
        double fitnorm1525 = BEexpol->GetParameter(3);
        double fitnorm1525_err = BEexpol->GetParError(3);
        double fitnorm1710 = BEexpol->GetParameter(6);
        double fitnorm1710_err = BEexpol->GetParError(6);
        double fitmass1525 = BEexpol->GetParameter(4);
        double fitmass1525_err = BEexpol->GetParError(4);
        double fitmass1710 = BEexpol->GetParameter(7);
        double fitmass1710_err = BEexpol->GetParError(7);
        double fitwidth1525 = BEexpol->GetParameter(5);
        double fitwidth1525_err = BEexpol->GetParError(5);
        double fitwidth1710 = BEexpol->GetParameter(8);
        double fitwidth1710_err = BEexpol->GetParError(8);
        double fitnorm1270 = BEexpol->GetParameter(0);
        double fitnorm1270_err = BEexpol->GetParError(0);
        double fitmass1270 = BEexpol->GetParameter(1);
        double fitmass1270_err = BEexpol->GetParError(1);
        double fitwidth1270 = BEexpol->GetParameter(2);
        double fitwidth1270_err = BEexpol->GetParError(2);
        double fitrangelow = BEexpol->GetXmin();
        double fitrangehigh = BEexpol->GetXmax();
        double expol1 = BEexpol->GetParameter(9);
        double expol2 = BEexpol->GetParameter(10);
        double expol3 = BEexpol->GetParameter(11);
        double expol4 = BEexpol->GetParameter(12);

        file << fitstatus << endl;
        file << std::fixed << std::setprecision(2);
        file << fitrangelow << " - " << fitrangehigh << endl;
        file << std::fixed << std::setprecision(1);
        file << chi2_ndf << endl;
        file << std::fixed << std::setprecision(0);
        file << fitnorm1525 << " ± " << fitnorm1525_err << endl;
        file << std::fixed << std::setprecision(3);
        file << fitmass1525 << " ± " << fitmass1525_err << endl;
        file << fitwidth1525 << " ± " << fitwidth1525_err << endl;
        file << std::fixed << std::setprecision(0);
        file << fitnorm1710 << " ± " << fitnorm1710_err << endl;
        file << std::fixed << std::setprecision(3);
        file << fitmass1710 << " ± " << fitmass1710_err << endl;
        file << fitwidth1710 << " ± " << fitwidth1710_err << endl;
        file << std::fixed << std::setprecision(0);
        file << fitnorm1270 << " ± " << fitnorm1270_err << endl;
        file << std::fixed << std::setprecision(3);
        file << fitmass1270 << " ± " << fitmass1270_err << endl;
        file << fitwidth1270 << " ± " << fitwidth1270_err << endl;
        file << endl;
#if defined(b_expol) || defined(b_expol1)
        file << expol1 << ", " << expol2 << ", " << expol3 << ", " << expol4 << endl;
#endif
#ifdef b_boltzman
        file << expol1 << ", " << expol2 << ", " << expol3 << endl;
#endif
        file << endl;
        file << endl;
        file << std::fixed << std::setprecision(4);
        file << "fit masses" << endl;
        file << fitmass1270 << ", " << fitmass1525 << ", " << fitmass1710 << endl;
        file << "fit masses errors" << endl;
        file << fitmass1270_err << ", " << fitmass1525_err << ", " << fitmass1710_err << endl;
        file << "fit widths" << endl;
        file << fitwidth1270 << ", " << fitwidth1525 << ", " << fitwidth1710 << endl;
        file << "fit widths errors" << endl;
        file << fitwidth1270_err << ", " << fitwidth1525_err << ", " << fitwidth1710_err << endl;

#ifdef residual_subtracted
        // Now subtract the residual background and plot
        TCanvas *c2 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c2, 0.14, 0.03, 0.05, 0.14);
        TH1F *hsubtracted = (TH1F *)hinvMass->Clone("hsubtracted");
        expol_clone->SetRange(0.99, 2.99);
        hsubtracted->Add(expol_clone, -1);
        hsubtracted->GetXaxis()->SetRangeUser(1.0, 2.5);
        hsubtracted->SetMaximum(hsubtracted->GetMaximum() * 1.5);
        hsubtracted->Draw();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }
        hsubtracted->Fit("onlyBW_clone", "REBMS");
        double *obtained_parameters2 = onlyBW_clone->GetParameters();
        TLine *line = new TLine(0.99, 0, 2.5, 0);
        line->SetLineColor(28);
        line->SetLineStyle(2);
        line->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits1[3];
        for (int i = 0; i < 3; i++)
        {
            singlefits1[i] = new TF1(Form("singlef%d", i), single_BW, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
            singlefits1[i]->SetParameter(0, obtained_parameters2[3 * i]);
            singlefits1[i]->SetParameter(1, obtained_parameters2[3 * i + 1]);
            singlefits1[i]->SetParameter(2, obtained_parameters2[3 * i + 2]);
            singlefits1[i]->SetLineColor(colors[i]);
            singlefits1[i]->SetLineStyle(2);
            singlefits1[i]->Draw("same");
        }

        TLegend *ltemp2 = new TLegend(0.20, 0.67, 0.42, 0.92);
        ltemp2->SetFillStyle(0);
        ltemp2->SetTextFont(42);
        ltemp2->SetTextSize(0.035);
        ltemp2->AddEntry(hsubtracted, "Data", "lpe");
        ltemp2->AddEntry(onlyBW_clone, "3rBW", "l");
        ltemp2->AddEntry(singlefits1[0], "f1270", "l");
        ltemp2->AddEntry(singlefits1[1], "f1525", "l");
        ltemp2->AddEntry(singlefits1[2], "f1710", "l");
        ltemp2->Draw("same");
        c2->SaveAs((savepath + "/boltzmann/rBWfit_residual.png").c_str());

#endif

        // **********************************************************************************************
        // *******************subtract the resonance peaks and fit the residual background*****************

        // // TCanvas *c3 = new TCanvas("", "", 720, 720);
        // // // Here we will subtract the resonances peaks and plot the residual background and the fit it
        // // SetCanvasStyle(c3, 0.14, 0.03, 0.05, 0.14);
        // // onlyBW->SetRange(f1270Mass - 3 * f1270Width, f1710Mass + 3 * f1710Width);
        // // hsubtracted_res->Add(onlyBW, -1);
        // // hsubtracted_res->Draw();
        // // hsubtracted_res->Fit("expol", "REBMSI");

        // // // Yield calculation
        // // double yield1270 = onlyBW_clone->Integral(f1270Mass - 2 * f1270Width, f1270Mass + 2 * f1270Width);
        // // double yield1320 = onlyBW_clone->Integral(a1320Mass - 2 * a1320Width, a1320Mass + 2 * a1320Width);
        // // double yield1525 = onlyBW_clone->Integral(f1525Mass - 2 * f1525Width, f1525Mass + 2 * f1525Width);
        // // double yield1710 = onlyBW_clone->Integral(f1710Mass - 2 * f1710Width, f1710Mass + 2 * f1710Width);

        // // double yield1270_err = onlyBW_clone->IntegralError((f1270Mass - 3 * f1270Width), (f1270Mass + 3 * f1270Width));
        // // double yield1320_err = onlyBW_clone->IntegralError((a1320Mass - 3 * a1320Width), (a1320Mass + 3 * a1320Width));
        // // double yield1525_err = onlyBW_clone->IntegralError((f1525Mass - 3 * f1525Width), (f1525Mass + 3 * f1525Width));
        // // double yield1710_err = onlyBW_clone->IntegralError((f1710Mass - 3 * f1710Width), (f1710Mass + 3 * f1710Width));

        // // cout << "Yield 1270: " << yield1270 << " +- " << yield1270_err << endl;
        // // cout << "Yield 1320: " << yield1320 << " +- " << yield1320_err << endl;
        // // cout << "Yield 1525: " << yield1525 << " +- " << yield1525_err << endl;
        // // cout << "Yield 1710: " << yield1710 << " +- " << yield1710_err << endl;
    }
}
// end of main program

// *****************************************************************************************************
//************************************Fit functions************************************************* */

// We will define the single BW and the sum of 3 BWs
Double_t single_BW_hera(double *x, double *par)
{
    // normalization factor is missing and how to add it I am not sure
    double amplitude = par[0];
    double mass = par[1];
    double width = par[2];

    double den = (x[0] * x[0] - mass * mass) * (x[0] * x[0] - mass * mass) + mass * mass * width * width;
    double realnum = (mass * mass - x[0] * x[0]) * mass * TMath::Sqrt(width);
    double imagnum = mass * mass * width * TMath::Sqrt(width);

    double real3BW = realnum / den;
    double imag3BW = imagnum / den;
    double sig1 = amplitude * (real3BW * real3BW + imag3BW * imag3BW);

    return sig1;
}

Double_t BWsum_hera(double *x, double *par) // taken 4 resonances here
{

    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2];
    double yield1525 = par[3];
    double mass1525 = par[4];
    double width1525 = par[5];
    double yield1710 = par[6];
    double mass1710 = par[7];
    double width1710 = par[8];

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270);
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525);
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710);

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270);
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525);
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710);

    double fit = yield1270 * (realnum1270 * realnum1270 + imagnum1270 * imagnum1270) / (den1270 * den1270) + yield1525 * (realnum1525 * realnum1525 + imagnum1525 * imagnum1525) / (den1525 * den1525) + yield1710 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / (den1710 * den1710);

    return fit;
}

Double_t single_BW(double *x, double *par)
{
    double yield = par[0];
    double mass = par[1];
    double width = par[2];

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));
    return fit;
}

Double_t BWsum(double *x, double *par)
{
    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2];
    double yield1525 = par[3];
    double mass1525 = par[4];
    double width1525 = par[5];
    double yield1710 = par[6];
    double mass1710 = par[7];
    double width1710 = par[8];

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    double fit = fit1270 + fit1525 + fit1710;
    return fit;
}

// Now we will define the functions for the exponential background

Double_t exponential_bkg_1(double *x, double *par) // 3 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * (x[0] - 2.0 * 0.497)));
}

Double_t exponential_bkg_2(double *x, double *par) // 3 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[1])));
}

Double_t exponential_bkg_3(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[3])));
}

Double_t exponential_bkg_4(double *x, double *par) // 5 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * ((pow((x[0] - 2.0 * 0.497), par[3])) + pow((x[0] - 2.0 * 0.497), par[4]))));
}

Double_t exponential_bkg_5(double *x, double *par) // 3 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * x[0]));
}

Double_t exponential_bkg_6(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] + par[1] * x[0] + par[3] * x[0] * x[0]));
}

Double_t expol_chkstar(double *x, double *par)
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] + x[0] * par[3]));
}

// Now we will define the functions for the Boltzmann background
Double_t Boltzmann_bkg_1(double *x, double *par) // 3 parameters
{
    double norm = par[0];
    double massks = 0.497;
    double n = par[1];
    double C = par[2];
    double y = norm * pow((x[0] - 2 * massks), n / 2) * pow(C, 3 / 2) * exp(-C * pow((x[0] - 2 * massks), n));
    return y;
}

Double_t Boltzmann_bkg_2(double *x, double *par) // 4 parameters
{
    double norm = par[0];
    double massks = 0.497;
    double n = par[1];
    double C = par[2];
    double n1 = par[3];
    double y = norm * pow((x[0] - 2 * massks), n / 2) * pow(C, 3 / 2) * exp(-C * pow((x[0] - 2 * massks), n1));
    return y;
}

// Now we will define the functions for the sum of the Breit-Wigner and the exponential background

Double_t single_BW_expol3(double *x, double *par)
{
    return (single_BW(x, par) + exponential_bkg_3(x, &par[3]));
}

Double_t single_BW_expol3_hera(double *x, double *par)
{
    return (single_BW_hera(x, par) + exponential_bkg_3(x, &par[3]));
}

Double_t BWsum_expol3(double *x, double *par)
{
    return (BWsum(x, par) + exponential_bkg_3(x, &par[9]));
}

Double_t BWsum_expol3_hera(double *x, double *par)
{
    return (BWsum_hera(x, par) + exponential_bkg_3(x, &par[9]));
}

// Now we will define the functions for the sum of the Breit-Wigner and the Boltzmann background

Double_t single_BW_boltzman_1(double *x, double *par)
{
    return (single_BW(x, par) + Boltzmann_bkg_1(x, &par[3]));
}

Double_t single_BW_boltzman_2(double *x, double *par)
{
    return (single_BW(x, par) + Boltzmann_bkg_2(x, &par[3]));
}

Double_t BWsum_boltzman_1(double *x, double *par)
{
    return (BWsum(x, par) + Boltzmann_bkg_1(x, &par[9]));
}

Double_t BWsum_boltzman_2(double *x, double *par)
{
    return (BWsum(x, par) + Boltzmann_bkg_2(x, &par[9]));
}

Double_t BWsum_expol_chkstar(double *x, double *par)
{
    return (BWsum(x, par) + expol_chkstar(x, &par[9]));
}
