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

Double_t BWchecks(double *x, double *par)
{
    // total 11 + 3 parameters, 2 for each resonance (total 4 resonance), 3 for normalization and 3 for background
    double mass1270 = par[0];
    double width1270 = par[1];
    double mass1320 = par[2];
    double width1320 = par[3];
    double mass1525 = par[4];
    double width1525 = par[5];
    double mass1710 = par[6];
    double width1710 = par[7];
    double a0 = par[8];
    double a1 = par[9];
    double a3 = par[10]; // we took a3 so as to not confuse with a2(1320)

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1320 = (x[0] * x[0] - mass1320 * mass1320) * (x[0] * x[0] - mass1320 * mass1320) + mass1320 * mass1320 * width1320 * width1320;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270);
    double realnum1320 = (mass1320 * mass1320 - x[0] * x[0]) * mass1320 * TMath::Sqrt(width1320);
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525);
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710);

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270);
    double imagnum1320 = mass1320 * mass1320 * width1320 * TMath::Sqrt(width1320);
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525);
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710);

    // double real3BW = 5 * realnum1270 / den1270 - 3 * realnum1320 / den1320 + 2 * a0 * realnum1525 / den1525;
    // double imag3BW = 5 * imagnum1270 / den1270 - 3 * imagnum1320 / den1320 + 2 * a0 * imagnum1525 / den1525;
    // double sig1 = (real3BW * real3BW + imag3BW * imag3BW);
    // double fit = a1 * sig1 + a3 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / den1710 * den1710;

    double fit1 = realnum1270 / den1270 - realnum1320 / den1320;
    double fit2 = imagnum1270 / den1270 - imagnum1320 / den1320;

    // double fit1 = realnum1270 / den1270 + realnum1320 / den1320 + realnum1525 / den1525;
    // double fit2 = imagnum1270 / den1270 + imagnum1320 / den1320 + imagnum1525 / den1525;

    // double fit1 = realnum1270 / den1270 + realnum1320 / den1320 + realnum1525 / den1525 + realnum1710 / den1710;
    // double fit2 = imagnum1270 / den1270 + imagnum1320 / den1320 + imagnum1525 / den1525 + imagnum1710 / den1710;

    double fit = fit1 * fit1 + fit2 * fit2;

    return fit;
}

Double_t RelativisticBW1(double *x, double *par)
{
    return (x[0] * par[0] * par[1] / (pow((x[0] * x[0] - par[0] * par[0]), 2) + pow((par[0] * par[1]), 2)));
}

Double_t RelativisticBW2(double *x, double *par)
{
    return (x[0] * par[0] * par[1] * par[2] / (pow((x[0] * x[0] - par[0] * par[0]), 2) + pow((par[0] * par[1]), 2)));
}

Double_t combineBW(double *x, double *par)
{
    // return (RelativisticBW1(x, &par[0]) + RelativisticBW1(x, &par[2]) + RelativisticBW1(x, &par[4]));
    // return (RelativisticBW1(x, &par[0]) + RelativisticBW1(x, &par[2]) + RelativisticBW1(x, &par[4]) + RelativisticBW1(x, &par[6]));
    return (5 * RelativisticBW1(x, &par[0]) - 3 * RelativisticBW1(x, &par[2]));
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

    // Here we are checking with the expol function
    //  TF1 *expol = new TF1("expol", exponential_bkg, 0, 10, 3);
    //  expol->SetParameter(0, -1000); // This tells the curvature. Positive and negative values has opposite curvatures. Magnitude of this parameters scales the overall curve
    //  expol->SetParameter(1, 0);
    //  expol->SetParameter(2, -0.5); // This parameter tells the shape. Positive value means exponential decay and negative value means exponential rise. Larger the value, faster the decay or rise. 0 means flat line.
    //  expol->Draw();

    // Here we are checking with BW which has both imag and real parts taken from HERA
    //  TF1 *coherentBW = new TF1("coherentBW", choerentBW_fitfunction_wo_bkg, 1, 2.2, 11);
    //  double parameter_temp[11] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, f1710Mass, f1710Width, 1, 1, 1000};
    //  // double parameter_temp[11] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, 1, 10};
    //  for (int i = 0; i < 11; i++)
    //  {
    //      coherentBW->FixParameter(i, parameter_temp[i]);
    //  }

    // coherentBW->Draw();

    // Here we are checking with the BW function with only the magnitude TF1 *coherentBW = new TF1("coherentBW", coherentBWsum, -3, 3, 14);
    // double parameter_temp[14] = {100, f1270Mass, f1270Width, 10, a1320Mass, a1320Width, 1, f1525Mass, f1525Width, 1, f1710Mass, f1710Width, 1, 1};
    // for (int i = 0; i < 14; i++)
    // {
    //     coherentBW->SetParameter(i, parameter_temp[i]);
    // }
    // coherentBW->Draw();

    // Here we are checking with the combinations of the BW functions using both real and imag parts (this is unsuccesful)
    TF1 *BWf1270_and_f1320 = new TF1("BWf1270_and_f1320", BWchecks, 1, 2.5, 8);
    double parameter_temp[11] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, f1710Mass, f1710Width};
    // double parameter_temp[11] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width};
    for (int i = 0; i < 11; i++)
    {
        BWf1270_and_f1320->SetParameter(i, parameter_temp[i]);
    }
    BWf1270_and_f1320->Draw();

    // // //Here we are checking with the single BW functions
    TF1 *sBW1 = new TF1("sBW1", singleBW, 0, 3, 2);
    sBW1->SetParameter(0, f1270Mass);
    sBW1->SetParameter(1, f1270Width);
    sBW1->SetLineColor(1);
    sBW1->SetLineStyle(2);
    sBW1->Draw("same");

    TF1 *sBW2 = new TF1("sBW2", singleBW, 0, 3, 2);
    sBW2->SetParameter(0, a1320Mass);
    sBW2->SetParameter(1, a1320Width);
    sBW2->SetLineColor(3);
    sBW2->SetLineStyle(2);
    sBW2->Draw("same");

    // TF1 *sBW3 = new TF1("sBW3", RelativisticBW1, 0, 3, 2);
    // sBW3->SetParameter(0, f1525Mass);
    // sBW3->SetParameter(1, f1525Width);
    // sBW3->SetLineColor(4);
    // sBW3->SetLineStyle(2);
    // sBW3->Draw("same");

    // TF1 *sBW4 = new TF1("sBW4", RelativisticBW1, 0, 3, 2);
    // sBW4->SetParameter(0, f1710Mass);
    // sBW4->SetParameter(1, f1710Width);
    // sBW4->SetLineColor(6);
    // sBW4->SetLineStyle(2);
    // sBW4->Draw("same");

    TLegend *leg = new TLegend(0.65, 0.6, 0.92, 0.9);
    leg->AddEntry(sBW1, "f1270", "l");
    leg->AddEntry(sBW2, "a1320", "l");
    // leg->AddEntry(sBW3, "1525", "l");
    // leg->AddEntry(sBW4, "1710", "l");
    leg->AddEntry(BWf1270_and_f1320, "f1270 + a1320 + f1525 + f1710", "l");
    leg->Draw("same");

    // // Here we are checking with the combinations of the BW functions using only the magnitude
    // TF1 *BWf1270_and_f1320 = new TF1("BWf1270_and_f1320", combineBW, 1, 2.5, 8);
    // // double parameter_temp[4] = {f1525Mass, f1525Width, f1710Mass, f1710Width};
    // double parameter_temp[8] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, f1710Mass, f1710Width};
    // for (int i = 0; i < 8; i++)
    // {
    //     BWf1270_and_f1320->SetParameter(i, parameter_temp[i]);
    // }
    // BWf1270_and_f1320->SetLineWidth(2);
    // BWf1270_and_f1320->GetYaxis()->SetRangeUser(-30, 30);
    // BWf1270_and_f1320->Draw();

    // TF1 *BWf1270 = new TF1("BWf1270", RelativisticBW2, 1, 2.5, 3);
    // BWf1270->SetParameter(0, f1270Mass);
    // BWf1270->SetParameter(1, f1270Width);
    // BWf1270->SetParameter(2, 5);
    // BWf1270->SetLineColor(1);
    // BWf1270->SetLineStyle(2);
    // BWf1270->Draw("same");

    // TF1 *BWf1320 = new TF1("BWf1320", RelativisticBW2, 1, 2.5, 3);
    // BWf1320->SetParameter(0, a1320Mass);
    // BWf1320->SetParameter(1, a1320Width);
    // BWf1320->SetParameter(2, -3);
    // BWf1320->SetLineColor(3);
    // BWf1320->SetLineStyle(2);
    // BWf1320->Draw("same");

    // // TF1 *BWf1525 = new TF1("BWf1525",RelativisticBW1 , 1, 2.5, 2);
    // // BWf1525->SetParameter(0, f1525Mass);
    // // BWf1525->SetParameter(1, f1525Width);
    // // BWf1525->SetLineColor(4);
    // // BWf1525->SetLineStyle(2);
    // // BWf1525->Draw("same");

    // // TF1 *BWf1710 = new TF1("BWf1710",RelativisticBW1 , 1, 2.5, 2);
    // // BWf1710->SetParameter(0, f1710Mass);
    // // BWf1710->SetParameter(1, f1710Width);
    // // BWf1710->SetLineColor(6);
    // // BWf1710->SetLineStyle(2);
    // // BWf1710->Draw("same");

    // TLegend *leg = new TLegend(0.55, 0.6, 0.92, 0.9);
    // leg->AddEntry(BWf1270, "f1270", "l");
    // leg->AddEntry(BWf1320, "f1320", "l");
    // // leg->AddEntry(BWf1525, "f1525", "l");
    // // leg->AddEntry(BWf1710, "f1710", "l");
    // leg->AddEntry(BWf1270_and_f1320, "f1270 - a1320", "l");
    // leg->Draw("same");

    c1->SaveAs("/home/sawan/Pictures/BWchecks.png");

    return 0;
}