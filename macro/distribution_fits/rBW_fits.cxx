#include <iostream>
#include "TDatabasePDG.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
#include <tuple>

TDatabasePDG *pdg = new TDatabasePDG();
double parameter0(double mass, double width);
TF1 *BW(TH1 *h, double mass, double width, double lowrange, double highrange);
TF1 *BW3expo(TH1 *h, double *parameters, double lowfitrange, double highfitrange);
TF1 *BW3pol3(TH1 *h, double *parameters, double lowfitrange, double highfitrange);
TF1 *BW3boltzman(TH1 *h, double *parameters, double lowfitrange, double highfitrange);
TF1 *BW3fit(TH1 *h, double *parameters, double lowfitrange, double highfitrange);
TF1 *CoherentBWexpol(TH1 *h, double *parameters, double lowfitrange, double highfitrange);  // this one is onld
TF1 *CoherentBWexpol2(TH1 *h, double *parameters, double lowfitrange, double highfitrange); // this one is from HERA thesis which includes both real and imaginary parts
TF1 *draw_individual_functions(TF1 *fit, double *parameters, TLegend *lfit, bool drawpol2 = true, string kbgfitfunction = "pol3");
void setParameterNames(TF1 *f3bw3);

void rBW_fits()
{
    // // *********************** constant parameters *****************************
    // const string kResBkg = "MIX";
    const string kResBkg = "ROTATED";
    // const string kResBkg = "LIKE";
    // const string kbgfitfunction = "pol3";
    // const string kbgfitfunction = "expol";
    // const string kbgfitfunction = "Boltzman";
    // const string kbgfitfunction = "CoherentBWsum";
    const string kbgfitfunction = "CoherentBWsum2";

    const int rebin = 2;
    bool testing = false;
    bool saveplots = true;
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);

    // path to the output folder
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

    // // *********************** Initializing graphs for pt different studies *****************************
    TGraphErrors *massf1270 = new TGraphErrors(Npt);
    TGraphErrors *widthf1270 = new TGraphErrors(Npt);
    TGraphErrors *massf1525 = new TGraphErrors(Npt);
    TGraphErrors *widthf1525 = new TGraphErrors(Npt);
    TGraphErrors *massf1710 = new TGraphErrors(Npt);
    TGraphErrors *widthf1710 = new TGraphErrors(Npt);
    TGraphErrors *yield1270 = new TGraphErrors(Npt);
    TGraphErrors *yield1525 = new TGraphErrors(Npt);
    TGraphErrors *yield1710 = new TGraphErrors(Npt);
    TGraphErrors *yield_bc[3] = {new TGraphErrors(Npt), new TGraphErrors(Npt), new TGraphErrors(Npt)};

    // TFile *f = new TFile((outputfolder_str + "/hglue_" + kResBkg + ".root").c_str(), "READ");
    // TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/253148/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60..root", "READ");
    // TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/253148/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60_full_ptrange_0.01MeV..root", "READ");
    TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60..root", "READ");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    // pT loop ***************************************************
    TFile *file_fitparams = new TFile((fits_folder_str + "fitparams_" + kResBkg + "_" + kbgfitfunction + ".root").c_str(), "RECREATE"); // file to save plots and parameters
    for (Int_t ip = pt_start; ip < pt_end; ip++)
    {

        TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", pT_bins[ip], pT_bins[ip + 1]));
        if (hinvMass == nullptr)
        {
            cout << "Error opening histogram" << endl;
            return;
        }
        hinvMass->Rebin(rebin);
        float binwidth = hinvMass->GetBinWidth(1);
        int noofevents = hinvMass->Integral();
        hinvMass->GetYaxis()->SetTitle(Form("Counts / %.2f MeV/c^{2}", binwidth));

        TCanvas *c1 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c1, 0.12, 0.03, 0.05, 0.14);
        hinvMass->GetYaxis()->SetTitleOffset(1.1);
        hinvMass->GetXaxis()->SetTitleOffset(1.3);
        hinvMass->GetXaxis()->SetRangeUser(1.0, 2.4);
        TH1F *hinvMassResSub = (TH1F *)hinvMass->Clone();
        hinvMass->Draw();
        // c1->SaveAs((fits_folder_str + Form("hinvMass_withoutfit_pt_%.1f_%.1f.png", pT_bins[ip], pT_bins[ip + 1])).c_str());
        double parameters1[9] = {50, f1270Mass, f1270Width, 25, f1525Mass, f1525Width, 25, f1710Mass, f1710Width};
        double parameters2[9] = {50, a1320Mass, a1320Width, 25, f1525Mass, f1525Width, 25, f1710Mass, f1710Width};
        double parameter_coherent[12] = {5000, f1270Mass, f1270Width, 500, a1320Mass, a1320Width, 250, f1525Mass, f1525Width, 25, f1710Mass, f1710Width};
        TGraphErrors *mass_all = new TGraphErrors();
        TGraphErrors *width_all = new TGraphErrors();
        TGraphErrors *yield_all = new TGraphErrors();
        mass_all->SetMarkerStyle(20);
        width_all->SetMarkerStyle(20);
        yield_all->SetMarkerStyle(20);

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
            // if (kchannel == "KsKs_Channel" && kbgfitfunction == "CoherentBWsum")
            // {
            //     // f3pol3->SetParameter(12, parameter_coherent[12]);  // overall amplitude
            //     // f3pol3->SetParameter(13, parameter_coherent[13]);  // norm of f1710

            //     // expol parameters
            //     f3pol3->SetParameter(14, 1e5);
            //     f3pol3->SetParameter(15, 0.2);
            //     f3pol3->SetParameter(16, 5);
            // }

            if (kchannel == "KsKs_Channel" && kbgfitfunction == "CoherentBWsum2")
            {

                double parameter_temp[] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, f1710Mass, f1710Width, 1, 5000, 500};

                struct FitParams
                {
                    double low;
                    double high;
                    double param0_limit;  // mass for f1270
                    double param1_limit;  // width for f1270
                    double param2_limit;  // mass for a1320
                    double param3_limit;  // width for a1320
                    double param4_limit;  // mass for f1525
                    double param5_limit;  // width for f1525
                    double param6_limit;  // mass for f1710
                    double param7_limit;  // width for f1710
                    double param8_limit;  // a0
                    double param9_limit;  // a1
                    double param10_limit; // a3
                };

                // Define the fit parameters for each pT bin // -1 free, -2 fixed
                std::vector<FitParams> bwfit_params_me = {
                    {1.03, 2.2, f1270Width, -1, 3 * a1320Width, -1, 3 * f1525Width, -1, 3 * f1710Width, -1, -2, -1, -1}, //
                };
                const auto &iter_bin = bwfit_params_me[ip];

                // Vector to store (low, high, chi2/NDF)
                std::vector<std::tuple<double, double, double>> chi2_results;
                std::vector<std::tuple<float, float, float, float>> chi2_results2;

                // int iteration = 0;
                // for (int ipar1 = 4e5; ipar1 < 9e5; ipar1 += 4e4) // loop for expol parameter 1
                // {
                //     for (double ipar2 = -0.1; ipar2 <= -0.01; ipar2 += 0.01) // loop for expol parameter 2
                //     {
                //         for (double ipar3 = 2; ipar3 < 4; ipar3 += 0.02) // loop for expol parameter 3
                //         {

                f3pol3 = CoherentBWexpol2(hinvMass, parameter_temp, iter_bin.low, iter_bin.high); //

                std::vector<std::pair<int, double>> limits = {
                    {0, iter_bin.param0_limit},
                    {1, iter_bin.param1_limit},
                    {2, iter_bin.param2_limit},
                    {3, iter_bin.param3_limit},
                    {4, iter_bin.param4_limit},
                    {5, iter_bin.param5_limit},
                    {6, iter_bin.param6_limit},
                    {7, iter_bin.param7_limit},
                    {8, iter_bin.param8_limit},
                    {9, iter_bin.param9_limit},
                    {10, iter_bin.param10_limit}};

                for (const auto &[param_idx, limit] : limits)
                {
                    if (parameter_temp[param_idx] != -1)
                        f3pol3->SetParameter(param_idx, parameter_temp[param_idx]); // Set parameter

                    if (limit != -1 && limit != 0 && limit != -2)
                    {
                        f3pol3->SetParLimits(param_idx, parameter_temp[param_idx] - limit, parameter_temp[param_idx] + limit); // Set limit
                    }
                    else if (limit == 0)
                    {
                        f3pol3->SetParLimits(param_idx, 0, 1e9); // Set limit for normalization
                    }
                    else if (limit == -2)
                    {
                        f3pol3->FixParameter(param_idx, parameter_temp[param_idx]); // Fix the parameter
                    }
                }

                //             // // Loop in expol parameters
                //             f3pol3->SetParameter(11, ipar1);
                //             f3pol3->SetParameter(12, ipar2);
                //             f3pol3->SetParameter(13, ipar3);
                //             hinvMass->Fit("f3pol3", "REBMS"); // use for loop
                //             double chi2ndf = f3pol3->GetChisquare() / f3pol3->GetNDF();
                //             cout << "chi2ndf: " << chi2ndf << endl;
                //             chi2_results2.push_back(std::make_tuple(ipar1, ipar2, ipar3, chi2ndf));
                //             iteration++;
                //             cout << "Iteration: " << iteration << endl;
                //             f3pol3->Clear();
                //         } // ipar3 loop end
                //     } // ipar2 loop end
                // } // ipar1 loop end

                // cout << "Total iterations: " << iteration << endl;

                // // sort in asceding order w.r.t to the third array i.e. chi2/NDF
                // sort(chi2_results2.begin(), chi2_results2.end(),
                //      [](const auto &a, const auto &b)
                //      {
                //          return get<3>(a) < get<3>(b);
                //      });

                // for (int i = 0; i < 30; i++)
                // {
                //     float best_ipar1 = std::get<0>(chi2_results2[i]);
                //     float best_ipar2 = std::get<1>(chi2_results2[i]);
                //     float best_ipar3 = std::get<2>(chi2_results2[i]);
                //     float best_chi2ndf = std::get<3>(chi2_results2[i]);
                //     cout << "ipar1: " << best_ipar1 << ",  ipar2: " << best_ipar2 << ",  ipar3: " << best_ipar3 << ",  chi2/NDF: " << best_chi2ndf << endl;
                // }

                // temporary fixing
                // f3pol3->FixParameter(6, 1.689);
                // f3pol3->FixParameter(7, 0.224);

                // //for expol parameters
                f3pol3->SetParameter(11, 760000); // default
                f3pol3->SetParameter(12, -0.09);  // This parameter shifts the left peak up or down. A lower value shifts up and vice versa
                f3pol3->SetParameter(13, 3.7);    // This scales the distribution up or down. A lower value scans up and vice versa
            }

            if (kchannel == "KsKs_Channel" && (kbgfitfunction == "Boltzman"))
            {

                // double parameter_temp[sizeof(parameters1) / sizeof(parameters1[0])]; // for f1270
                // std::copy(std::begin(parameters1), std::end(parameters1), parameter_temp);
                // cout<<"used parameter is f1270<<endl;

                double parameter_temp[sizeof(parameters2) / sizeof(parameters2[0])]; // for a1320
                std::copy(std::begin(parameters2), std::end(parameters2), parameter_temp);
                cout << "used parameter is a1320" << endl;

                struct FitParams
                {
                    double low;
                    double high;
                    double param0_limit; // norm for f1270
                    double param1_limit; // mass for f1270
                    double param2_limit; // width for f1270
                    double param3_limit; // norm for f1525
                    double param4_limit; // mass for f1525
                    double param5_limit; // width for f1525
                    double param6_limit; // norm for f1710
                    double param7_limit; // mass for f1710
                    double param8_limit; // width for f1710
                };

                // Define the fit parameters for each pT bin // -1 free, -2 fixed
                std::vector<FitParams> bwfit_params_me = {
                    {1.09, 2.17, 0, -2, 0.08, 0, 10 * f1525Width, -1, 0, -1, -1}, // free mass and width of a1320
                    // {1.09, 2.17, 0, 10 * f1270Width, 0.08, 0, 10 * f1525Width, -1, 0, -1, -1}, // free mass and width of f1270
                    // {1.09, 2.17, 0, -2, -2, 0, 10 * f1525Width, -1, 0, -1, -1},                // fix mass and width of f1270
                };
                const auto &iter_bin = bwfit_params_me[ip];

                // Vector to store (low, high, chi2/NDF)
                std::vector<std::tuple<double, double, double>> chi2_results;
                std::vector<std::tuple<float, float, float, float>> chi2_results2;

                // // // Loop through low and high ranges
                // for (double low = 0.98; low <= 1.10; low += 0.01)
                // {
                //     for (double high = 2.1; high <= 2.8; high += 0.01)
                //     {
                // Update the fitting range
                // f3pol3 = BW3boltzman(hinvMass, parameter_temp, low, high); // use for fitting range loop
                f3pol3 = BW3boltzman(hinvMass, parameter_temp, iter_bin.low, iter_bin.high); // without fitting range loop

                // int iteration = 0;
                // for (int ipar1 = 1e4; ipar1 < 1e7; ipar1 *= 2) // loop for expol parameter 1
                // {
                //     for (double ipar2 = -1; ipar2 <= 1; ipar2 += 0.1) // loop for expol parameter 2
                //     {
                //         for (double ipar3 = 2; ipar3 < 8; ipar3 += 0.2) // loop for expol parameter 3
                //         {

                std::vector<std::pair<int, double>> limits = {
                    {0, iter_bin.param0_limit},
                    {1, iter_bin.param1_limit},
                    {2, iter_bin.param2_limit},
                    {3, iter_bin.param3_limit},
                    {4, iter_bin.param4_limit},
                    {5, iter_bin.param5_limit},
                    {6, iter_bin.param6_limit},
                    {7, iter_bin.param7_limit},
                    {8, iter_bin.param8_limit},
                };

                for (const auto &[param_idx, limit] : limits)
                {
                    f3pol3->SetParameter(param_idx, parameter_temp[param_idx]); // Set parameter
                    if (limit != -1 && limit != 0 && limit != -2)
                    {
                        f3pol3->SetParLimits(param_idx, parameter_temp[param_idx] - limit, parameter_temp[param_idx] + limit); // Set limit
                    }
                    else if (limit == 0)
                    {
                        f3pol3->SetParLimits(param_idx, 0, 1e9); // Set limit for normalization
                    }
                    else if (limit == -2)
                    {
                        f3pol3->FixParameter(param_idx, parameter_temp[param_idx]); // Fix the parameter
                    }
                }

                //             // // Loop in expol parameters
                //             f3pol3->SetParameter(9, ipar1);
                //             f3pol3->SetParameter(10, ipar2);
                //             f3pol3->SetParameter(11, ipar3);
                //             hinvMass->Fit("f3pol3", "REBMS"); // use for loop
                //             double chi2ndf = f3pol3->GetChisquare() / f3pol3->GetNDF();
                //             cout << "chi2ndf: " << chi2ndf << endl;
                //             chi2_results2.push_back(std::make_tuple(ipar1, ipar2, ipar3, chi2ndf));
                //             iteration++;
                //             cout << "Iteration: " << iteration << endl;
                //         } // ipar3 loop end
                //     } // ipar2 loop end
                // } // ipar1 loop end

                // cout << "Total iterations: " << iteration << endl;
                // // Sort the vector to get the 10 lowest chi2/NDF values
                // size_t n = std::min(chi2_results2.size(), static_cast<size_t>(10));

                // // Partially sort to get the 10 smallest elements based on chi2/NDF
                // std::partial_sort(chi2_results2.begin(), chi2_results2.begin() + n, chi2_results2.end(),
                //                   [](const auto &a, const auto &b)
                //                   {
                //                       return std::get<3>(a) < std::get<3>(b);
                //                   });

                // // Output the 10 best (low, high) combinations with their chi2/NDF values
                // std::cout << "Top " << n << " lowest chi2/NDF values:\n";
                // for (size_t i = 0; i < n; ++i)
                // {
                //     float best_ipar1 = std::get<0>(chi2_results2[i]);
                //     float best_ipar2 = std::get<1>(chi2_results2[i]);
                //     float best_ipar3 = std::get<2>(chi2_results2[i]);
                //     float best_chi2ndf = std::get<3>(chi2_results2[i]);
                //     cout << " ipar1: " << best_ipar1 << ",  ipar2: " << best_ipar2 << ",  ipar3: " << best_ipar3 << ",  chi2/NDF: " << best_chi2ndf << endl;
                // } // here ends the loop for expol parameters

                // //for Boltzman parameters
                f3pol3->SetParameter(9, 520000);
                f3pol3->SetParameter(10, 0.5); // n
                f3pol3->SetParameter(11, 5.5); // c

                //         // Fit the function
                //         hinvMass->Fit("f3pol3", "REBMS"); // use for loop

                //         // Get chi2 and NDF (number of degrees of freedom)
                //         double chi2 = f3pol3->GetChisquare();
                //         double ndf = f3pol3->GetNDF();
                //         double chi2_ndf = chi2 / ndf;
                //         cout << "chi2/NDF is: " << chi2 / ndf << endl;

                //         // Store (low, high, chi2/NDF)
                //         chi2_results.push_back(std::make_tuple(low, high, chi2_ndf));
                //     } // high range loop
                // } // low range loop

                // // Sort the vector to get the 10 lowest chi2/NDF values
                // size_t n = std::min(chi2_results.size(), static_cast<size_t>(10)); // Take only the top 10 or less if size is smaller

                // // Partially sort to get the 10 smallest elements based on chi2/NDF
                // std::partial_sort(chi2_results.begin(), chi2_results.begin() + n, chi2_results.end(),
                //                   [](const auto &a, const auto &b)
                //                   {
                //                       return std::get<2>(a) < std::get<2>(b);
                //                   });

                // for (size_t i = 0; i < n; ++i)
                // {
                //     double low = std::get<0>(chi2_results[i]);
                //     double high = std::get<1>(chi2_results[i]);
                //     double chi2_ndf = std::get<2>(chi2_results[i]);
                //     std::cout << "Low: " << low << ", High: " << high << ", chi2/NDF: " << chi2_ndf << std::endl;
                // }
            }

            if (kchannel == "KsKs_Channel" && (kbgfitfunction == "expol"))
            {
                double parameter_temp[sizeof(parameters1) / sizeof(parameters1[0])]; // for f1270
                std::copy(std::begin(parameters1), std::end(parameters1), parameter_temp);
                cout << "used parameter is f1270" << endl;

                // double parameter_temp[sizeof(parameters2) / sizeof(parameters2[0])]; // for a1320
                // std::copy(std::begin(parameters2), std::end(parameters2), parameter_temp);
                // cout<<"used parameter is a1320"<<endl;

                struct FitParams
                {
                    double low;
                    double high;
                    double param0_limit; // norm for f1270
                    double param1_limit; // mass for f1270
                    double param2_limit; // width for f1270
                    double param3_limit; // norm for f1525
                    double param4_limit; // mass for f1525
                    double param5_limit; // width for f1525
                    double param6_limit; // norm for f1710
                    double param7_limit; // mass for f1710
                    double param8_limit; // width for f1710
                };

                // Define the fit parameters for each pT bin
                std::vector<FitParams> bwfit_params_me = {
                    {1.01, 2.2, -1, 5 * f1270Mass, 0.1, -1, 5 * f1525Width, -1, -1, -1, -1}, // free mass and width of f1270

                    // {1.01, 2.2, 0, 10 * a1320Width, 0.05, 0, 10 * f1525Width, -1, 0, -1, -1}, // free mass and width of f1270

                    // {1.01, 2.2, 0, -2, -2, 0, 10 * f1525Width, -1, 0, -1, -1}, // fix mass and width of f1270
                };
                const auto &iter_bin = bwfit_params_me[ip];

                // Vector to store (low, high, chi2/NDF)
                vector<tuple<double, double, double>> chi2_results;
                vector<tuple<float, float, float, float>> chi2_results2;

                // // // Loop through low and high ranges
                // for (double low = 0.98; low <= 1.10; low += 0.01)
                // {
                //     for (double high = 2.1; high <= 2.8; high += 0.1)
                //     {
                //         // Update the fitting range
                //         f3pol3 = BW3expo(hinvMass, parameter_temp, low, high); // use for fitting range loop

                // int iteration = 0;
                // for (int ipar1 = 1e4; ipar1 < 1e7; ipar1 *= 10) // loop for expol parameter 1
                // {
                //     for (double ipar2 = -1; ipar2 <= 1; ipar2 += 1) // loop for expol parameter 2
                //     {
                //         for (double ipar3 = 2; ipar3 < 8; ipar3 += 2) // loop for expol parameter 3
                //         {
                f3pol3 = BW3expo(hinvMass, parameter_temp, iter_bin.low, iter_bin.high); // without fitting range loop

                std::vector<std::pair<int, double>> limits = {
                    {0, iter_bin.param0_limit},
                    {1, iter_bin.param1_limit},
                    {2, iter_bin.param2_limit},
                    {3, iter_bin.param3_limit},
                    {4, iter_bin.param4_limit},
                    {5, iter_bin.param5_limit},
                    {6, iter_bin.param6_limit},
                    {7, iter_bin.param7_limit},
                    {8, iter_bin.param8_limit},
                };

                // Loop through the parameter index and par limits
                for (const auto &[param_idx, limit] : limits)
                {
                    f3pol3->SetParameter(param_idx, parameter_temp[param_idx]); // Set parameter
                    if (limit != -1 && limit != 0 && limit != -2)
                    {
                        f3pol3->SetParLimits(param_idx, parameter_temp[param_idx] - limit, parameter_temp[param_idx] + limit); // Set limit
                    }
                    else if (limit == 0)
                    {
                        f3pol3->SetParLimits(param_idx, 0, 1e9); // Set limit for normalization
                    }
                    else if (limit == -2)
                    {
                        f3pol3->FixParameter(param_idx, parameter_temp[param_idx]); // Fix the parameter
                    }
                }

                //             // // Loop in expol parameters
                //             f3pol3->SetParameter(9, ipar1);
                //             f3pol3->SetParameter(10, ipar2);
                //             f3pol3->SetParameter(11, ipar3);
                //             hinvMass->Fit("f3pol3", "REBMS"); // use for loop
                //             double chi2ndf = f3pol3->GetChisquare() / f3pol3->GetNDF();
                //             cout << "chi2ndf: " << chi2ndf << endl;
                //             chi2_results2.push_back(std::make_tuple(ipar1, ipar2, ipar3, chi2ndf));
                //             iteration++;
                //             cout << "Iteration: " << iteration << endl;
                //             f3pol3->Clear();
                //         } // ipar3 loop end
                //     } // ipar2 loop end
                // } // ipar1 loop end

                // cout << "Total iterations: " << iteration << endl;

                // // sort in asceding order w.r.t to the third array i.e. chi2/NDF
                // sort(chi2_results2.begin(), chi2_results2.end(),
                //      [](const auto &a, const auto &b)
                //      {
                //          return get<3>(a) < get<3>(b);
                //      });

                // for (int i = 0; i < 20; i++)
                // {
                //     float best_ipar1 = std::get<0>(chi2_results2[i]);
                //     float best_ipar2 = std::get<1>(chi2_results2[i]);
                //     float best_ipar3 = std::get<2>(chi2_results2[i]);
                //     float best_chi2ndf = std::get<3>(chi2_results2[i]);
                //     cout << "ipar1: " << best_ipar1 << ",  ipar2: " << best_ipar2 << ",  ipar3: " << best_ipar3 << ",  chi2/NDF: " << best_chi2ndf << endl;
                // }

                // for expol parameters
                f3pol3->SetParameter(9, 1e6); // default 1e6
                f3pol3->SetParameter(10, 0);
                f3pol3->SetParameter(11, 6);

                //         // Fit the function
                //         hinvMass->Fit("f3pol3", "REBMS"); // use for loop

                //         // Get chi2 and NDF (number of degrees of freedom)
                //         double chi2 = f3pol3->GetChisquare();
                //         double ndf = f3pol3->GetNDF();
                //         double chi2_ndf = chi2 / ndf;
                //         cout << "chi2/NDF is: " << chi2_ndf << endl;

                //         // Store (low, high, chi2/NDF)
                //         chi2_results.push_back(std::make_tuple(low, high, chi2_ndf));
                //         f3pol3->Clear();
                //     } // high range loop
                // } // low range loop

                // // Sort in ascending order based on chi2/NDF
                // sort(chi2_results.begin(), chi2_results.end(),
                //                   [](const auto &a, const auto &b)
                //                   {
                //                       return std::get<2>(a) < std::get<2>(b);
                //                   });

                // for (int i = 0; i < 20; i++)
                // {
                //     double low = std::get<0>(chi2_results[i]);
                //     double high = std::get<1>(chi2_results[i]);
                //     double chi2_ndf = std::get<2>(chi2_results[i]);
                //     std::cout << "Low: " << low << ", High: " << high << ", chi2/NDF: " << chi2_ndf << std::endl;
                // }
                // cout << endl;
            }

            // Drawing the fit
            f3pol3->SetLineStyle(1);
            f3pol3->SetLineWidth(2);
            // hinvMass->SetMaximum(hinvMass->GetMaximum() * 1.4);
            // hinvMass->SetMinimum(-1000);
            hinvMass->Fit("f3pol3", "REBMS");
            cout << "chi2/ndf: " << f3pol3->GetChisquare() / f3pol3->GetNDF() << endl;
            f3pol3->Draw("same");
            int mass_index[3] = {1, 4, 7};
            int width_index[3] = {2, 5, 8};

            for (int i = 0; i < 3; i++)
            {
                mass_all->SetPoint(i, i + 1, f3pol3->GetParameter(mass_index[i]));
                mass_all->SetPointError(i, 0, f3pol3->GetParError(mass_index[i]));
                width_all->SetPoint(i, i + 1, f3pol3->GetParameter(width_index[i]));
                width_all->SetPointError(i, 0, f3pol3->GetParError(width_index[i]));
            }

            // making the size of textbox optimal
            gPad->Update();
            TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
            ptstats->SetX1NDC(0.6);
            ptstats->SetX2NDC(0.99);
            ptstats->SetY1NDC(0.5);
            ptstats->SetY2NDC(0.95);
            ptstats->Draw("same");

            // Initializting the individual rBW functions
            TF1 *fitFuncs[3] = {
                new TF1("fit1270", RelativisticBW, f3pol3->GetXmin(), f3pol3->GetXmax(), 3),
                new TF1("fit1525", RelativisticBW, f3pol3->GetXmin(), f3pol3->GetXmax(), 3),
                new TF1("fit1710", RelativisticBW, f3pol3->GetXmin(), f3pol3->GetXmax(), 3)};

            // Initializing the background functions
            TF1 *fitpol3 = new TF1("fitpol3", polynomial3, f3pol3->GetXmin(), f3pol3->GetXmax(), 4);
            TF1 *fitexpol = new TF1("fitexpol", exponential_bkg, f3pol3->GetXmin(), f3pol3->GetXmax(), 3);
            TF1 *fitboltzman = new TF1("fitboltzman", Boltzman, f3pol3->GetXmin(), f3pol3->GetXmax(), 3);

            // Setting the parameters of the individual rBW functions
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    fitFuncs[i]->SetParameter(j, f3pol3->GetParameter(i * 3 + j));
            }

            // Setting the parameters of the background functions
            for (int i = 0; i < 4; i++)
            {
                if (kbgfitfunction == "pol3")
                    fitpol3->SetParameter(i, f3pol3->GetParameter(i + 9));
                else if (kbgfitfunction == "expol" && i < 3)
                    fitexpol->SetParameter(i, f3pol3->GetParameter(i + 9));
                else if (kbgfitfunction == "Boltzman" && i < 3)
                    fitboltzman->SetParameter(i, f3pol3->GetParameter(i + 9));
            }

            double *parameters3 = f3pol3->GetParameters();
            double parameter_expo[3] = {parameters3[11], parameters3[12], parameters3[13]};
            double parameter_BW[11] = {parameters3[0], parameters3[1], parameters3[2], parameters3[3], parameters3[4], parameters3[5], parameters3[6], parameters3[7], parameters3[8], parameters3[9], parameters3[10]};
            TF1 *residualbkg = new TF1("residualbkg", exponential_bkg, f3pol3->GetXmin(), f3pol3->GetXmax(), 3);
            residualbkg->SetParameter(0, parameters3[11]);
            residualbkg->SetParameter(1, parameters3[12]);
            residualbkg->SetParameter(2, parameters3[13]);
            residualbkg->SetLineColor(kBlue);
            residualbkg->SetLineStyle(2);
            residualbkg->SetLineWidth(2);
            residualbkg->Draw("same");
            TF1 *BWfits = new TF1("BWfits", choerentBW_fitfunction_wo_bkg, f3pol3->GetXmin(), f3pol3->GetXmax(), 11);
            for (int i = 0; i < 11; i++)
            {
                BWfits->SetParameter(i, parameter_BW[i]);
            }
            BWfits->SetParLimits(0, parameter_BW[0] - 2 * f1270Width, parameter_BW[0] + 2 * f1270Width);
            BWfits->SetLineColor(kMagenta);
            BWfits->SetLineStyle(2);
            BWfits->SetLineWidth(2);
            BWfits->Draw("same");

            TLegend *ltemp = new TLegend(0.20, 0.67, 0.52, 0.92);
            ltemp->SetFillStyle(0);
            ltemp->SetTextFont(42);
            // ltemp->SetBorderSize(0);
            ltemp->SetTextSize(0.035);
            ltemp->AddEntry(hinvMass, "Data", "lpe");
            ltemp->AddEntry(f3pol3, Form("4rBW + %s", "expol"), "l");
            ltemp->AddEntry(BWfits, "4rBW", "l");
            ltemp->AddEntry(residualbkg, Form("%s", "expol"), "l");
            ltemp->Draw("same");

            // // // //Drawing the individual rBW functions and the legend
            // TLegend *lfit = new TLegend(0.25, 0.72, 0.55, 0.92);
            // lfit->SetFillStyle(0);
            // lfit->SetTextFont(42);
            // lfit->SetTextSize(0.025);
            // lfit->AddEntry(hinvMass, "Data", "lpe");
            // lfit->AddEntry(f3pol3, Form("3rBW + %s", kbgfitfunction.c_str()), "l");
            // double *parameters3 = f3pol3->GetParameters();
            // TF1 *residualbkg = draw_individual_functions(f3pol3, parameters3, lfit, true, kbgfitfunction);
            // lfit->Draw("same");
            t2->SetNDC();
            // t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip], pT_bins[ip + 1]));
            if (saveplots)
            {
                c1->SaveAs((fits_folder_str + "_" + kResBkg + "_" + kbgfitfunction + Form("_%.1f_%.1f.png", pT_bins[ip], pT_bins[ip + 1])).c_str());
            }

            // ************* subtract residual background and the plot *************
            TCanvas *c2 = new TCanvas("", "", 2432, 85, 720, 720);
            SetCanvasStyle(c2, 0.12, 0.03, 0.05, 0.14);
            residualbkg->SetRange(0.99, 2.99);    // extrpolating the fitting range
            hinvMassResSub->Add(residualbkg, -1); // subtract the residual background
            float rangelow = f3pol3->GetXmin();
            float rangehigh = f3pol3->GetXmax();
            hinvMassResSub->GetXaxis()->SetRangeUser(0.95, 2.6);
            hinvMassResSub->SetMinimum(-100);
            hinvMassResSub->Draw();

            TF1 *f3bw3 = new TF1("f3bw3", choerentBW_fitfunction_wo_bkg, f3pol3->GetXmin(), f3pol3->GetXmax(), 11);
            f3bw3->SetParNames("mass f1270", "width f1270", "mass a1320", "width a1320", "mass f1525", "width f1525", "mass f1710", "width f1710", "a0", "a1", "a3");
            for (int i = 0; i < 11; i++)
            {
                f3bw3->SetParameter(i, parameter_BW[i]);
            }

            f3bw3->SetLineColor(kRed);
            f3bw3->SetLineStyle(1);
            hinvMassResSub->Fit("f3bw3", "REBMS");
            hinvMassResSub->SetMarkerSize(0.8);
            hinvMassResSub->SetMaximum(hinvMassResSub->GetMaximum() * 1.3);
            cout << "Subtracted plot Chi2/ndf: " << f3bw3->GetChisquare() / f3bw3->GetNDF() << endl;

            // for (int i = 0; i < 3; i++)
            // {
            //     mass_all->SetPoint(i + 3, i + 4, f3bw3->GetParameter(mass_index[i]));
            //     mass_all->SetPointError(i + 3, 0, f3bw3->GetParError(mass_index[i]));
            //     width_all->SetPoint(i + 3, i + 4, f3bw3->GetParameter(width_index[i]));
            //     width_all->SetPointError(i + 3, 0, f3bw3->GetParError(width_index[i]));
            // }
            // mass_all->Write("mass_all");
            // width_all->Write("width_all");

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
            // draw_individual_functions(f3pol3, parameters3, lfit2, false, kbgfitfunction); // draw the individual functions

            // save the plot
            t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip], pT_bins[ip + 1]));
            if (saveplots)
            {
                c2->SaveAs((fits_folder_str + "_" + kResBkg + "_" + kbgfitfunction + "res_sub_" + Form("_%.1f_%.1f.png", pT_bins[ip], pT_bins[ip + 1])).c_str());
            }

            // **************************** Compute the yield and plot *********************************

            auto setPoint = [&](TGraphErrors *graph, int paramIndex)
            {
                double pT_center = (pT_bins[ip] + pT_bins[ip + 1]) / 2;
                double pT_width = (pT_bins[ip + 1] - pT_bins[ip]) / 2.0;
                graph->SetPoint(ip, pT_center, f3pol3->GetParameter(paramIndex));
                graph->SetPointError(ip, pT_width, f3pol3->GetParError(paramIndex));
            };

            setPoint(massf1270, 1);
            setPoint(widthf1270, 2);
            setPoint(massf1525, 4);
            setPoint(widthf1525, 5);
            setPoint(massf1710, 7);
            setPoint(widthf1710, 8);

            float ptbinwidth = pT_bins[ip + 1] - pT_bins[ip];
            double normFactor = binwidth * ptbinwidth * noofevents;

            TVirtualFitter::SetFitter(0); // Reset the fitter to avoid conflicts
            auto computeYield = [&](TF1 *func, TGraphErrors *graph, int index)
            {
                TVirtualFitter::Fitter(func); // Set the fitter for this TF1
                double yield = func->Integral(f3pol3->GetXmin(), f3pol3->GetXmax()) / normFactor;
                // double yieldErr = func->IntegralError(f3pol3->GetXmin(), f3pol3->GetXmax()) / normFactor;
                double pT_center = (pT_bins[ip + 1] + pT_bins[ip]) / 2;
                double pT_width = (pT_bins[ip + 1] - pT_bins[ip]) / 2.0;
                graph->SetPoint(ip, pT_center, yield);
                // graph->SetPointError(ip, pT_width, yieldErr);
                // std::cout << "yield" << index << ": " << yield * normFactor << std::endl;
            };

            computeYield(fitFuncs[0], yield1270, 1270);
            computeYield(fitFuncs[1], yield1525, 1525);
            computeYield(fitFuncs[2], yield1710, 1710);

            int resonanceno[] = {1270, 1525, 1710};
            double mass[] = {f1270Mass, f1525Mass, f1710Mass};
            double width[] = {f1270Width, f1525Width, f1710Width};
            int bin_min[3], bin_max[3];
            double bc_errors[3], bkgValues[3], intBW[3], sumTailCorr[3], totalY[3], yields[3];

            for (int i = 0; i < 3; i++)
            {
                bin_min[i] = hinvMass->FindBin(mass[i] - 2 * width[i]);
                bin_max[i] = hinvMass->FindBin(mass[i] + 2 * width[i]);
                yields[i] = hinvMass->IntegralAndError(bin_min[i], bin_max[i], bc_errors[i]);
            }

            for (int i = 0; i < 3; ++i)
            {
                // cout << "bc_errors: " << bc_errors[i] << endl;
                bkgValues[i] = fitpol3->Integral(hinvMass->GetBinLowEdge(bin_min[i]), hinvMass->GetBinLowEdge(bin_max[i] + 1));
                intBW[i] = fitFuncs[i]->Integral(hinvMass->GetBinLowEdge(bin_min[i]), hinvMass->GetBinLowEdge(bin_max[i] + 1));
                sumTailCorr[i] = (fitFuncs[i]->Integral(f3pol3->GetXmin(), hinvMass->GetBinLowEdge(bin_min[i])) + fitFuncs[i]->Integral(hinvMass->GetBinLowEdge(bin_max[i] + 1), f3pol3->GetXmax())) / binwidth;
                totalY[i] = (sumTailCorr[i] + yields[i] - (bkgValues[i] / binwidth)) / (ptbinwidth * noofevents);
                // std::cout << "Total_Ybincounting" << resonanceno[i] << " " << totalY[i] << std::endl;
                yield_bc[i]->SetPoint(ip, (pT_bins[ip] + pT_bins[ip + 1]) / 2, totalY[i]);
                yield_bc[i]->SetPointError(ip, (pT_bins[ip + 1] - pT_bins[ip]) / 2, bc_errors[i] / (ptbinwidth * noofevents));
            }

            // c2->Close();
        }

    } //************************ end of pT loop **************************************

    // ************* Plot the graphs of mass, width and yield for pT different studies ***************
    if (Npt > 1)
    {
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

        TCanvas *c9 = new TCanvas("", "", 720, 720);
        gPad->SetLogy();
        SetCanvasStyle(c9, 0.15, 0.03, 0.05, 0.14);
        SetGrapherrorStyle(yield1270);
        SetGrapherrorStyle(yield_bc[0]);
        yield_bc[0]->SetMarkerStyle(21);
        yield_bc[0]->SetMarkerColor(kRed);
        yield_bc[0]->SetLineColor(kRed);
        yield1270->GetYaxis()->SetTitle("Yield");
        yield1270->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        yield1270->SetMaximum(0.5);
        // yield1270->SetMinimum(0);
        yield1270->Draw("ape");
        // yield_bc[0]->Draw("pe same");
        yield1270->Write("yield1270");

        TLegend *lfit4 = new TLegend(0.65, 0.65, 0.9, 0.9);
        lfit4->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
        lfit4->AddEntry(yield1270, "Function integration", "lpe");
        // lfit4->AddEntry(yield_bc[0], "Bin counting", "lpe");
        lfit4->Draw("same");
        if (saveplots)
            c9->SaveAs((fits_folder_str + "_" + kResBkg + "yield1270.png").c_str());

        TCanvas *c10 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c10, 0.15, 0.03, 0.05, 0.14);
        gPad->SetLogy();
        SetGrapherrorStyle(yield1525);
        SetGrapherrorStyle(yield_bc[1]);
        yield_bc[1]->SetMarkerStyle(21);
        yield_bc[1]->SetMarkerColor(kRed);
        yield_bc[1]->SetLineColor(kRed);
        yield1525->GetYaxis()->SetTitle("Yield");
        yield1525->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        yield1525->SetMaximum(0.5);
        // yield1525->SetMinimum(0);
        yield1525->Draw("ape");
        // yield_bc[1]->Draw("pe same");
        yield1525->Write("yield1525");
        lfit4->Clear();
        lfit4->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
        lfit4->AddEntry(yield1525, "Function integration", "lpe");
        // lfit4->AddEntry(yield_bc[1], "Bin counting", "lpe");
        lfit4->Draw("same");
        if (saveplots)
            c10->SaveAs((fits_folder_str + "_" + kResBkg + "yield1525.png").c_str());

        TCanvas *c11 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c11, 0.15, 0.03, 0.05, 0.14);
        gPad->SetLogy();
        SetGrapherrorStyle(yield1710);
        SetGrapherrorStyle(yield_bc[2]);
        yield_bc[2]->SetMarkerStyle(21);
        yield_bc[2]->SetMarkerColor(kRed);
        yield_bc[2]->SetLineColor(kRed);
        yield1710->GetYaxis()->SetTitle("Yield");
        yield1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        yield1710->SetMaximum(0.5);
        // yield1710->SetMinimum(0);
        yield1710->Draw("ape");
        // yield_bc[2]->Draw("pe same");
        yield1710->Write("yield1710");
        lfit4->Clear();
        lfit4->AddEntry((TObject *)0, "Run 3 pp #sqrt{s} = 13 TeV", "");
        lfit4->AddEntry(yield1710, "Function integration", "lpe");
        // lfit4->AddEntry(yield_bc[2], "Bin counting", "lpe");
        lfit4->Draw("same");
        if (saveplots)
            c11->SaveAs((fits_folder_str + "_" + kResBkg + "yield1710.png").c_str());
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
    TF1 *f3pol3 = new TF1("f3pol3", expo3bW, lowfitrange, highfitrange, 12);

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

    // h->Fit("f3pol3", "REBMS0");
    return f3pol3;
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

TF1 *CoherentBWexpol(TH1 *h, double *parameters, double lowfitrange, double highfitrange)
{
    TF1 *f3pol3 = new TF1("f3pol3", coherentBWsum_expol, lowfitrange, highfitrange, 17);
    f3pol3->SetParName(0, "norm1");
    f3pol3->SetParName(1, "mass1");
    f3pol3->SetParName(2, "width1");
    f3pol3->SetParName(3, "norm2");
    f3pol3->SetParName(4, "mass2");
    f3pol3->SetParName(5, "width2");
    f3pol3->SetParName(6, "norm3");
    f3pol3->SetParName(7, "mass3");
    f3pol3->SetParName(8, "width3");
    f3pol3->SetParName(9, "norm4");
    f3pol3->SetParName(10, "mass4");
    f3pol3->SetParName(11, "width4");
    f3pol3->SetParName(12, "Amplitude"); // overall normalization
    f3pol3->SetParName(13, "Amp f1710"); // overall Normalization of the f1710 during the coherence
    f3pol3->SetParName(14, "A");         // expol first parameter
    f3pol3->SetParName(15, "B");
    f3pol3->SetParName(16, "C");

    return f3pol3;
}

TF1 *BW3boltzman(TH1 *h, double *parameters, double lowfitrange, double highfitrange)
{
    TF1 *f3pol3 = new TF1("f3pol3", BW3boltzman, lowfitrange, highfitrange, 12);

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

    // h->Fit("f3pol3", "REBMS0");
    return f3pol3;
}

TF1 *CoherentBWexpol2(TH1 *h, double *parameters, double lowfitrange, double highfitrange)
{
    TF1 *f3pol3 = new TF1("f3pol3", choerentBW_fitfunction, lowfitrange, highfitrange, 14);
    f3pol3->SetParName(0, "mass1");
    f3pol3->SetParName(1, "width1");
    f3pol3->SetParName(2, "mass2");
    f3pol3->SetParName(3, "width2");
    f3pol3->SetParName(4, "mass3");
    f3pol3->SetParName(5, "width3");
    f3pol3->SetParName(6, "mass4");
    f3pol3->SetParName(7, "width4");
    f3pol3->SetParName(8, "a0");
    f3pol3->SetParName(9, "a1");
    f3pol3->SetParName(10, "a3");
    f3pol3->SetParName(11, "expol 1");
    f3pol3->SetParName(12, "expol 2");
    f3pol3->SetParName(13, "expol 3");

    return f3pol3;
}

TF1 *draw_individual_functions(TF1 *fit, double *parameters, TLegend *lfit, bool drawpol2 = true, string kbgfitfunction = "pol3")
{
    // BW1
    float rangelow = fit->GetXmin();
    float rangehigh = fit->GetXmax();
    TF1 *frBW1 = new TF1("frBW1", RelativisticBW, rangelow, rangehigh, 3);
    TF1 *frBW2 = new TF1("frBW2", RelativisticBW, rangelow, rangehigh, 3);
    TF1 *frBW3 = new TF1("frBW3", RelativisticBW, rangelow, rangehigh, 3);
    TF1 *frBW4 = new TF1("frBW4", RelativisticBW, rangelow, rangehigh, 3);

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

    // BW4
    if (kbgfitfunction == "CoherentBWsum")
    {
        frBW4->SetParameter(0, parameters[9]);
        frBW4->SetParameter(1, parameters[10]);
        frBW4->SetParameter(2, parameters[11]);
        frBW4->SetLineStyle(2);
        frBW4->SetLineColor(9);
        frBW4->Draw("same");
    }

    // pol3 background
    TF1 *fpol3 = new TF1("fpol3", polynomial3, rangelow, rangehigh, 4);
    fpol3->SetParameter(0, parameters[9]);
    fpol3->SetParameter(1, parameters[10]);
    fpol3->SetParameter(2, parameters[11]);
    fpol3->SetParameter(3, parameters[12]);
    fpol3->SetLineStyle(2);
    fpol3->SetLineColor(28);

    // expol background
    TF1 *fexpol = new TF1("fexpol", exponential_bkg, rangelow, rangehigh, 3);
    if (kbgfitfunction != "CoherentBWsum")
    {
        fexpol->SetParameter(0, parameters[9]);
        fexpol->SetParameter(1, parameters[10]);
        fexpol->SetParameter(2, parameters[11]);
    }
    else
    {
        fexpol->SetParameter(0, parameters[14]);
        fexpol->SetParameter(1, parameters[15]);
        fexpol->SetParameter(2, parameters[16]);
    }
    fexpol->SetLineStyle(2);
    fexpol->SetLineColor(28);

    // Boltzman background
    TF1 *fBoltzman = new TF1("fBoltzman", Boltzman, rangelow, rangehigh, 3);
    fBoltzman->SetParameter(0, parameters[9]);
    fBoltzman->SetParameter(1, parameters[10]);
    fBoltzman->SetParameter(2, parameters[11]);
    fBoltzman->SetLineStyle(2);
    fBoltzman->SetLineColor(28);

    if (drawpol2 && kbgfitfunction == "pol3")
    {
        fpol3->Draw("same");
        lfit->AddEntry(fpol3, "pol3", "l");
    }
    else if (drawpol2 && (kbgfitfunction == "expol" || kbgfitfunction == "CoherentBWsum"))
    {
        fexpol->Draw("same");
        lfit->AddEntry(fexpol, "expol", "l");
    }
    else if (drawpol2 && kbgfitfunction == "Boltzman")
    {
        fBoltzman->Draw("same");
        lfit->AddEntry(fBoltzman, "Boltzman", "l");
    }
    lfit->AddEntry(frBW1, "rBW(f_{2}(1270))", "l");
    lfit->AddEntry(frBW2, "rBW(f_{2}(1525))", "l");
    lfit->AddEntry(frBW3, "rBW(f_{0}(1710))", "l");

    if (kbgfitfunction == "pol3")
    {
        return fpol3;
    }
    else if (kbgfitfunction == "expol")
    {
        return fexpol;
    }
    else if (kbgfitfunction == "Boltzman")
    {
        return fBoltzman;
    }
    else
    {
        return fexpol;
    }
}

void setParameterNames(TF1 *f3bw3)
{
    f3bw3->SetParNames("Norm", "Mass_{f1270}", "#Gamma_{f1270}", "Norm_{f1525}", "Mass_{f1525}", "#Gamma_{f1525}", "Norm_{f1710}", "Mass_{f1710}", "#Gamma_{f1710}");
}
