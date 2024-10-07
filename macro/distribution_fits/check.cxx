#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
using namespace std;


Double_t singleBW(double *x, double *par)
{
    // total 11 + 3 parameters, 2 for each resonance (total 4 resonance), 3 for normalization and 3 for background
    double mass1270 = par[0];
    double width1270 = par[1];

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270);

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270);

    // real part 
    double real = realnum1270 / den1270;

    // imaginary part 
    double imag = imagnum1270 / den1270;

    // cross section of first 3 BW
    double sig1 = (real * real + imag * imag);

    return sig1;
}

int check()
{
    // // //checking the usage of tuple
    // vector<tuple<int, int, int, int>> check;

    // int arr1[5] = {1, 2, 3, 4, 5};      // First element of the tuple
    // int arr2[5] = {10, 20, 30, 40, 50}; // Second element of the tuple
    // int arr3[5] = {9, 8, 8, 2, 6};      // Third element of the tuple
    // int arr4[5] = {5, 10, 50, 20, 25};  // Fourth element of the tuple

    // // Create tuples using elements at corresponding positions in the arrays
    // for (int i = 0; i < 5; i++)
    // {
    //     check.push_back(make_tuple(arr1[i], arr2[i], arr3[i], arr4[i]));
    // }

    // // sort in asceding order w.r.t to the third array i.e. chi2/NDF
    // sort(check.begin(), check.end(),
    //      [](const auto &a, const auto &b)
    //      {
    //          return std::get<3>(a) < std::get<3>(b);
    //      });

    // for (int i = 0; i < 5; i++)
    // {
    //     float best_ipar1 = std::get<0>(check[i]);
    //     float best_ipar2 = std::get<1>(check[i]);
    //     float best_ipar3 = std::get<2>(check[i]);
    //     float best_chi2ndf = std::get<3>(check[i]);
    //     cout << "ipar1: " << best_ipar1 << ",  ipar2: " << best_ipar2 << ",  ipar3: " << best_ipar3 << ",  chi2/NDF: " << best_chi2ndf << endl;
    // }

    // checking the fitting function expol
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    // TF1 *expol = new TF1("expol", exponential_bkg, 0, 10, 3);
    // expol->SetParameter(0, -1000); // This tells the curvature. Positive and negative values has opposite curvatures. Magnitude of this parameters scales the overall curve
    // expol->SetParameter(1, 0);
    // expol->SetParameter(2, -0.5); // This parameter tells the shape. Positive value means exponential decay and negative value means exponential rise. Larger the value, faster the decay or rise. 0 means flat line.
    // expol->Draw();

    // TF1 *coherentBW = new TF1("coherentBW", choerentBW_fitfunction_wo_bkg, 1, 2.2, 11);
    // double parameter_temp[11] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, f1710Mass, f1710Width, 1, 1, 1000};
    // // double parameter_temp[11] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, 1, 10};
    // for (int i = 0; i < 11; i++)
    // {
    //     coherentBW->FixParameter(i, parameter_temp[i]);
    // }
    
    // coherentBW->Draw();

    // TF1 *coherentBW = new TF1("coherentBW", coherentBWsum, -3, 3, 14);
    // double parameter_temp[14] = {100, f1270Mass, f1270Width, 10, a1320Mass, a1320Width, 1, f1525Mass, f1525Width, 1, f1710Mass, f1710Width, 1, 1};
    // for (int i = 0; i < 14; i++)
    // {
    //     coherentBW->SetParameter(i, parameter_temp[i]);
    // }
    // coherentBW->Draw();

    TF1 *sBW1 = new TF1("sBW1", singleBW, 0, 3, 2);
    sBW1->SetParameter(0, f1270Mass);
    sBW1->SetParameter(1, f1270Width);
    sBW1->GetYaxis()->SetRangeUser(0, 14);
    sBW1->Draw();

    TF1 *sBW2 = new TF1("sBW2", singleBW, 0, 3, 2);
    sBW2->SetParameter(0, a1320Mass);
    sBW2->SetParameter(1, a1320Width);
    sBW2->SetLineColor(1);
    sBW2->Draw("same");

    TF1 *sBW3 = new TF1("sBW3", singleBW, 0, 3, 2);
    sBW3->SetParameter(0, f1525Mass);
    sBW3->SetParameter(1, f1525Width);
    sBW3->SetLineColor(3);
    sBW3->Draw("same");

    TF1 *sBW4 = new TF1("sBW4", singleBW, 0, 3, 2);
    sBW4->SetParameter(0, f1710Mass);
    sBW4->SetParameter(1, f1710Width);
    sBW4->SetLineColor(4);
    sBW4->Draw("same");

    TLegend *leg = new TLegend(0.1, 0.6, 0.4, 0.9);
    leg->AddEntry(sBW1, "1270", "l");
    leg->AddEntry(sBW2, "1320", "l");
    leg->AddEntry(sBW3, "1525", "l");
    leg->AddEntry(sBW4, "1710", "l");
    leg->Draw("same");



    return 0;
}