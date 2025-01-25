#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <TArrow.h>
#include "../../src/common_glue.h"
#include "../../src/fitting_range_glue.h"
#include "../../src/style.h"
using namespace std;

Double_t single_BW_hera(double *x, double *par);
Double_t single_BW(double *x, double *par);
Double_t BWsum_hera(double *x, double *par);
Double_t BWsum(double *x, double *par);
Double_t BWsum2(double *x, double *par);

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

void plot_ressub2()
{

    ofstream outfile("fit_results.txt");

    // TFile *file1 = new TFile("3rBW_plots_boltzmann.root", "READ");
    // TFile *file2 = new TFile("3rBW_plots_expol.root", "READ");
    // TFile *file3 = new TFile("3rBW_plots_exponential.root", "READ");

    // // TFile *file1 = new TFile("4rBW_plots_boltzmann.root", "READ");
    TFile *file1 = new TFile("3rBW_plots_expol.root", "READ");
    TFile *file2 = new TFile("4rBW_plots_expol.root", "READ");
    TFile *file3 = new TFile("4rBW_plots_exponential.root", "READ");

    if (file1->IsZombie() || file2->IsZombie() || file3->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hsubtracted1 = (TH1F *)file1->Get("3BW");
    TH1F *hsubtracted2 = (TH1F *)file2->Get("3BW");
    TH1F *hsubtracted3 = (TH1F *)file3->Get("3BW");

    if (hsubtracted1 == nullptr || hsubtracted2 == nullptr || hsubtracted3 == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }

    TCanvas *c2 = new TCanvas("", "", 720, 720);

    // Now lets fit the distribution
    TF1 *fit3BW = new TF1("fit3BW", BWsum, 1.05, 2.20, 9);
    fit3BW->SetParameters(6646, 1.281, 0.156, 7547, 1.514, 0.080, 4553, 1.692, 0.211); // expol 1
    // fit3BW->SetParameters(14233, 1.283, 0.211, 9206, 1.512, 0.090, 8208, 1.683, 0.266); // Boltzmann
    fit3BW->SetLineStyle(2);
    hsubtracted1->Fit("fit3BW", "REBMS");

    double chi2ndf1 = fit3BW->GetChisquare() / fit3BW->GetNDF();
    cout << "Chi2/NDF1 is " << chi2ndf1 << endl;
    outfile << "Boltzmann fitting parameters" << endl;
    for (int i = 0; i < 9; i++)
    {
        if (i == 0 || i == 3 || i == 6)
        {
            outfile << std::fixed << std::setprecision(0);
        }
        else
        {
            outfile << std::fixed << std::setprecision(3);
        }
        outfile << fit3BW->GetParameter(i) << " ± " << fit3BW->GetParError(i) << endl;
    }

    // fit3BW->SetLineColor(kBlue);
    // hsubtracted2->Fit("fit3BW", "REBMS");
    outfile << "Expol fitting parameters" << endl;
    TF1 *fit4BW = new TF1("fit4BW", BWsum2, 1.05, 2.20, 12);
    // fit4BW->SetParameters(3312, 1.218, 0.141, 3959, 1.305, 0.104, 8233, 1.513, 0.085, 4431, 1.694, 0.203); // expol 1
    vector<float> parameters = {3312, 1.218, 0.141, 3959, 1.305, 0.104, 8233, 1.513, 0.085, 4431, 1.694, 0.203};
    for (int i = 0; i < 12; i++)
    {
        fit4BW->SetParameter(i, parameters[i]);
    }

    fit4BW->SetLineStyle(2);
    fit4BW->SetLineColor(kBlue);
    hsubtracted2->Fit("fit4BW", "REBMS");

    double chi2ndf2 = fit4BW->GetChisquare() / fit4BW->GetNDF();
    cout << "Chi2/NDF2 is " << chi2ndf2 << endl;

    for (int i = 0; i < 12; i++)
    {
        if (i == 0 || i == 3 || i == 6 || i == 9)
        {
            outfile << std::fixed << std::setprecision(0);
        }
        else
        {
            outfile << std::fixed << std::setprecision(3);
        }
        outfile << fit4BW->GetParameter(i) << " ± " << fit4BW->GetParError(i) << endl;
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.145, 0.03, 0.05, 0.14);
    hsubtracted1->GetXaxis()->SetRangeUser(1.0, 2.5);
    // hsubtracted1->SetMaximum(hsubtracted1->GetMaximum() * 0.7);
    hsubtracted1->SetLineColor(kRed);
    hsubtracted1->SetMarkerColor(kRed);
    hsubtracted1->SetStats(0);
    hsubtracted1->SetMarkerStyle(22);
    hsubtracted1->SetMarkerSize(1.0);
    hsubtracted1->Draw("pe");
    hsubtracted2->SetLineColor(kBlue);
    hsubtracted2->SetMarkerColor(kBlue);
    hsubtracted2->SetMarkerStyle(23);
    hsubtracted2->SetMarkerSize(1.0);
    hsubtracted2->SetStats(0);
    hsubtracted2->Draw("pe SAME");
    // hsubtracted3->SetLineColor(kGreen);
    // hsubtracted3->SetMarkerColor(kGreen);
    // hsubtracted3->SetStats(0);
    // hsubtracted3->Draw("pe SAME");

    TLegend *leg = new TLegend(0.5, 0.55, 0.9, 0.75);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    // leg->AddEntry(hsubtracted1, "3rBW + Boltzmann", "lpe");
    // leg->AddEntry(hsubtracted2, "3rBW + Expol1", "lpe");
    leg->AddEntry(hsubtracted1, "3rBW + Expol1", "lpe");
    leg->AddEntry(hsubtracted2, "4rBW + Expol1", "lpe");
    leg->AddEntry(fit3BW, Form("#chi^{2}/ndf = %.2f", chi2ndf1), "l");
    leg->AddEntry(fit4BW, Form("#chi^{2}/ndf = %.2f", chi2ndf2), "l");
    leg->Draw("same");

    TLegend *leg2 = new TLegend(0.25, 0.75, 0.9, 0.90);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.04);
    leg2->AddEntry((TObject *)0, "M_{inv} after residual", "");
    leg2->AddEntry((TObject *)0, "background subtraction", "");
    leg2->AddEntry((TObject *)0, "Original distribution was fit with:", "");
    leg2->Draw("same");

    c->SaveAs("plot_ressub_compare.png");
}

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

Double_t BWsum2(double *x, double *par)
{
    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2];
    double yield1320 = par[3];
    double mass1320 = par[4];
    double width1320 = par[5];
    double yield1525 = par[6];
    double mass1525 = par[7];
    double width1525 = par[8];
    double yield1710 = par[9];
    double mass1710 = par[10];
    double width1710 = par[11];

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1320 = yield1320 * mass1320 * width1320 * x[0] / (pow((x[0] * x[0] - mass1320 * mass1320), 2) + pow(mass1320 * width1320, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    double fit = fit1270 + fit1320 + fit1525 + fit1710;
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