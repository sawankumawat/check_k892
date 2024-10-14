#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

Double_t choerentBW_fitfunction_wo_bkg1(double *x, double *par) // taken 4 resonances here
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
    double real3BW = -realnum1270 / den1270 + realnum1320 / den1320 + a0 * realnum1525 / den1525;

    // double imag3BW = 5 * imagnum1270 / den1270 - 3 * imagnum1320 / den1320 + 2 * a0 * imagnum1525 / den1525;
    double imag3BW = -imagnum1270 / den1270 + imagnum1320 / den1320 + a0 * imagnum1525 / den1525;

    double sig1 = (real3BW * real3BW + imag3BW * imag3BW);

    double fit = a1 * sig1 + a3 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / (den1710 * den1710);
    // double fit_with_bkg = fit + exponential_bkg_1(x, &par[11]);

    return fit;
}

Double_t choerentBW_fitfunction_wo_bkg2(double *x, double *par) // taken 3 resonances here
{
    double mass1270 = par[0];
    double width1270 = par[1];
    double mass1525 = par[2];
    double width1525 = par[3];
    double mass1710 = par[4];
    double width1710 = par[5];
    double a0 = par[6];
    double a1 = par[7];
    double a3 = par[8];

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270);
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525);
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710);

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270);
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525);
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710);

    double real3BW = a0 * realnum1270 / den1270 + a1 * realnum1525 / den1525 + a3 * realnum1710 / den1710;
    double imag3BW = a0 * imagnum1270 / den1270 + a1 * imagnum1525 / den1525 + a3 * imagnum1710 / den1710;

    double fit = (real3BW * real3BW + imag3BW * imag3BW);

    return fit;
}

Double_t BWsum(double *x, double *par)
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
    double norm_interference = par[12];
    double norm_f1710 = par[13];

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1320 = yield1320 * mass1320 * width1320 * x[0] / (pow((x[0] * x[0] - mass1320 * mass1320), 2) + pow(mass1320 * width1320, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    // double fit = fit1270 + fit1320 + fit1525 + fit1710;
    double fit = norm_interference * pow((fit1320 - fit1270 + fit1525), 2) + norm_f1710 * fit1710 * fit1710;
    return fit;
}

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

Double_t BWandexpol1(double *x, double *par)
{
    return choerentBW_fitfunction_wo_bkg1(x, par) + exponential_bkg_1(x, &par[11]);
}

Double_t BWandexpol2(double *x, double *par)
{
    return choerentBW_fitfunction_wo_bkg1(x, par) + exponential_bkg_2(x, &par[11]);
}

Double_t BWandexpol3(double *x, double *par)
{
    return choerentBW_fitfunction_wo_bkg1(x, par) + exponential_bkg_3(x, &par[11]);
}

Double_t BWandexpol4(double *x, double *par)
{
    return choerentBW_fitfunction_wo_bkg1(x, par) + exponential_bkg_4(x, &par[11]);
}

Double_t BWandexpol5(double *x, double *par)
{
    return choerentBW_fitfunction_wo_bkg1(x, par) + exponential_bkg_5(x, &par[11]);
}

Double_t BWandexpol6(double *x, double *par)
{
    return choerentBW_fitfunction_wo_bkg1(x, par) + exponential_bkg_6(x, &par[11]);
}

Double_t BWsumexpol1(double *x, double *par)
{
    return BWsum(x, par) + exponential_bkg_1(x, &par[14]);
}

Double_t BWsumexpol3(double *x, double *par)
{
    return BWsum(x, par) + exponential_bkg_3(x, &par[14]);
}

void check2()
{
    string savepath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/fits/fit_func1_check2_code";
    TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60..root", "READ");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", pT_bins[0], pT_bins[1]));
    if (hinvMass == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }
    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
    hinvMass->Rebin(2);
    hinvMass->Draw();
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);

    // // // // **************** For Coherent BW sum ****************************8
    // // TF1 *BEexpol = new TF1("BEexpol", BWandexpol1, 1.03, 2.30, 15); // expol 1
    // TF1 *BEexpol = new TF1("BEexpol", BWandexpol3, 1.03, 2.30, 15); // expol 3
    // // TF1 *BEexpol = new TF1("BEexpol", BWandexpol4, 1.03, 2.30, 16); // expol 4
    // // TF1 *BEexpol = new TF1("BEexpol", BWandexpol6, 1.03, 2.3, 16); // expol 6
    // BEexpol->SetParNames("mass1270", "width1270", "mass1320", "width1320", "mass1525", "width1525", "mass1710", "width1710", "a0", "a1", "a3");
    // double parameters[] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, f1710Mass, f1710Width, 0.8, 2300, 3000};
    // for (int i = 0; i < 11; i++)
    // {
    //     BEexpol->SetParameter(i, parameters[i]);
    // }
    // vector<vector<float>> par_limits = {{0, 5 * f1270Width}, {2, 5 * a1320Width}, {4, 5 * f1525Width}, {6, 5 * f1710Width}};
    // int limits_size = par_limits.size();
    // for (int i = 0; i < limits_size; i++)
    // {
    //     int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
    //     BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
    // }
    // // BEexpol->SetParameter(11, 710000); // expol 1
    // // BEexpol->SetParameter(12, -0.03); // expol 1
    // // BEexpol->SetParameter(13, 2.78); // expol 1

    // BEexpol->SetParameter(11, 710000); // expol 3
    // BEexpol->SetParameter(12, -0.03);  // expol 3
    // BEexpol->SetParameter(13, 2.78);   // expol 3
    // BEexpol->SetParameter(14, 1.15);   // expol 3

    // // BEexpol->SetParameter(11, 710000); // expol 4
    // // BEexpol->SetParameter(12, -0.03);  // expol 4
    // // BEexpol->SetParameter(13, 2.78);   // expol 4
    // // BEexpol->SetParameter(14, 1.13);   // expol 4
    // // BEexpol->SetParameter(15, 2); // expol 4

    // // BEexpol->SetParameter(11, 710000); // expol 6
    // // BEexpol->SetParameter(12, -0.03);  // expol 6
    // // BEexpol->SetParameter(13, 2.7);    // expol 6
    // // BEexpol->SetParameter(14, 1.16);   // expol 6
    // // BEexpol->SetParameter(15, 2.05);   // expol 6

    // hinvMass->Fit("BEexpol", "REBMS");

    // double *obtained_parameters = BEexpol->GetParameters();
    // // TF1 *expol = new TF1("expol", exponential_bkg_1, 1.03, 2.2, 3); // expol 1
    // TF1 *expol = new TF1("expol", exponential_bkg_3, 1.03, 2.2, 4); // expol 3
    // // TF1 *expol = new TF1("expol", exponential_bkg_4, 1.03, 2.2, 5); // expol 4
    // // TF1 *expol = new TF1("expol", exponential_bkg_6, 1.03, 2.2, 5); // expol 6
    // expol->SetParameter(0, obtained_parameters[11]);
    // expol->SetParameter(1, obtained_parameters[12]);
    // expol->SetParameter(2, obtained_parameters[13]);
    // expol->SetParameter(3, obtained_parameters[14]); // for expol 3
    // // expol->SetParameter(4, obtained_parameters[15]); // for expol 4 and 6
    // expol->SetLineColor(3);
    // expol->SetLineStyle(2);
    // expol->Draw("same");

    // TF1 *onlyBW = new TF1("onlyBW", choerentBW_fitfunction_wo_bkg1, 1.03, 2.2, 11);
    // TF1 *onlyBW_clone = new TF1("onlyBW_clone", choerentBW_fitfunction_wo_bkg1, 1.03, 2.2, 11);
    // onlyBW_clone->SetParNames("mass1270", "width1270", "mass1320", "width1320", "mass1525", "width1525", "mass1710", "width1710", "a0", "a1", "a3");
    // for (int i = 0; i < 11; i++)
    // {
    //     // onlyBW->SetParameter(i, parameters[i]);
    //     onlyBW->SetParameter(i, obtained_parameters[i]);
    //     onlyBW_clone->SetParameter(i, obtained_parameters[i]);
    // }
    // onlyBW->SetLineColor(4);
    // onlyBW->SetLineStyle(2);
    // onlyBW->Draw("same");

    // TLegend *ltemp = new TLegend(0.20, 0.67, 0.52, 0.92);
    // ltemp->SetFillStyle(0);
    // ltemp->SetTextFont(42);
    // ltemp->SetTextSize(0.035);
    // ltemp->AddEntry(hinvMass, "Data", "lpe");
    // ltemp->AddEntry(BEexpol, "4rBW + expol", "l");
    // ltemp->AddEntry(onlyBW, "4rBW", "l");
    // ltemp->AddEntry(expol, "expol", "l");
    // ltemp->Draw("same");

    ////// ***************** For BW sum with exponential background ******************
    // TF1 *BEexpol = new TF1("BEexpol", BWsumexpol1, 1.03, 2.30, 17);
    TF1 *BEexpol = new TF1("BEexpol", BWsumexpol3, 1.03, 2.30, 18);
    BEexpol->SetParNames("yield1270", "mass1270", "width1270", "yield1320", "mass1320", "width1320", "yield1525", "mass1525", "width1525");
    double parameters[] = {10, f1270Mass, f1270Width, 10, a1320Mass, a1320Width, 10, f1525Mass, f1525Width, 10, f1710Mass, f1710Width, 7, 50};
    for (int i = 0; i < 14; i++)
    {
        BEexpol->SetParameter(i, parameters[i]);
    }
    vector<vector<float>> par_limits = {{1, 5 * f1270Width}, {4, 5 * a1320Width}, {7, 5 * f1525Width}, {10, 5 * f1710Width}};
    int limits_size = par_limits.size();
    for (int i = 0; i < limits_size; i++)
    {
        int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
        BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
    }
    BEexpol->SetParameter(14, 610000); // expol 3
    BEexpol->SetParameter(15, -0.071); // expol 3
    BEexpol->SetParameter(16, 2.597);  // expol 3
    BEexpol->SetParameter(17, 1.121);  // expol 3
    // BEexpol->FixParameter(1, f1270Mass);
    // BEexpol->FixParameter(11, f1710Width);
    hinvMass->Fit("BEexpol", "REBMS");

    double *obtained_parameters = BEexpol->GetParameters();
    // TF1 *expol = new TF1("expol", exponential_bkg_1, 1.03, 2.2, 3);
    TF1 *expol = new TF1("expol", exponential_bkg_3, 1.03, 2.2, 4);
    expol->SetParameter(0, obtained_parameters[14]);
    expol->SetParameter(1, obtained_parameters[15]);
    expol->SetParameter(2, obtained_parameters[16]);
    expol->SetParameter(3, obtained_parameters[17]);
    expol->SetLineColor(3);
    expol->SetLineStyle(2);
    expol->Draw("same");

    TF1 *onlyBW = new TF1("onlyBW", BWsum, 1.03, 2.2, 14);
    TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum, 1.03, 2.2, 14);
    onlyBW_clone->SetParNames("yield1270", "mass1270", "width1270", "yield1320", "mass1320", "width1320", "yield1525", "mass1525", "width1525", "yield1710", "mass1710");
    for (int i = 0; i < 14; i++)
    {
        // onlyBW->SetParameter(i, parameters[i]);
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
    ltemp->AddEntry(BEexpol, "4rBW + expol", "l");
    ltemp->AddEntry(onlyBW, "4rBW", "l");
    ltemp->AddEntry(expol, "expol", "l");
    ltemp->Draw("same");

    // ************ common for all fits ******************

    gPad->Update();
    TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
    ptstats->SetX1NDC(0.5);
    ptstats->SetX2NDC(0.99);
    ptstats->SetY1NDC(0.4);
    ptstats->SetY2NDC(0.92);
    ptstats->Draw("same");
    c->SaveAs((savepath + "/rBWfit.png").c_str());

    // Now subtract the residual background and plot
    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.14, 0.03, 0.05, 0.14);
    TH1F *hsubtracted = (TH1F *)hinvMass->Clone("hsubtracted");
    expol->SetRange(0.99, 2.99);
    hsubtracted->Add(expol, -1);
    hsubtracted->GetXaxis()->SetRangeUser(0.98, 2.5);
    hsubtracted->Draw();
    // onlyBW_clone->SetParLimits(1, f1270Mass - 0.5 * f1270Width, f1270Mass + 0.5 * f1270Width);
    // onlyBW_clone->SetParLimits(4, a1320Mass - 0.5 * a1320Width, a1320Mass + 0.5 * a1320Width);
    hsubtracted->Fit("onlyBW_clone", "REBMS");
    TLine *line = new TLine(0.99, 0, 2.5, 0);
    line->SetLineColor(28);
    line->SetLineStyle(2);
    line->Draw("same");
    c2->SaveAs((savepath + "/rBWfit_residual.png").c_str());
}
