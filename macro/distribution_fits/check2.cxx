#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

Double_t choerentBW_fitfunction_wo_bkg1(double *x, double *par)
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
    double real3BW = realnum1270 / den1270 + realnum1320 / den1320 + a0 * realnum1525 / den1525;

    // double imag3BW = 5 * imagnum1270 / den1270 - 3 * imagnum1320 / den1320 + 2 * a0 * imagnum1525 / den1525;
    double imag3BW = imagnum1270 / den1270 + imagnum1320 / den1320 + a0 * imagnum1525 / den1525;

    double sig1 = (real3BW * real3BW + imag3BW * imag3BW);

    double fit = a1 * sig1 + a3 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / (den1710 * den1710);
    // double fit_with_bkg = fit + exponential_bkg(x, &par[11]);

    return fit;
}

Double_t choerentBW_fitfunction_wo_bkg2(double *x, double *par)
{
    // total 11 + 3 parameters, 2 for each resonance (total 4 resonance), 3 for normalization and 3 for background
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

Double_t BWandexpol1(double *x, double *par)
{
    return choerentBW_fitfunction_wo_bkg1(x, par) + exponential_bkg(x, &par[11]);
}

Double_t BWandexpol2(double *x, double *par)
{
    return choerentBW_fitfunction_wo_bkg2(x, par) + exponential_bkg(x, &par[9]);
}

void check2()
{
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
    hinvMass->Rebin(2);
    hinvMass->Draw();
    gStyle->SetOptFit(1111);

    TF1 *BEexpol = new TF1("BEexpol", BWandexpol2, 1.02, 2.3, 12);
    BEexpol->SetParNames("mass1270", "width1270", "mass1525", "width1525", "mass1710", "width1710", "a0", "a1", "a3", "expol 1", "expol 2");
    double parameters[] = {f1270Mass, f1270Width, f1525Mass, f1525Width, f1710Mass, f1710Width, 10, 10, 10};
    for (int i = 0; i < 9; i++)
    {
        BEexpol->SetParameter(i, parameters[i]);
    }
    BEexpol->SetParameter(9, 760000);
    BEexpol->SetParameter(10, -0.09);
    BEexpol->SetParameter(11, 3.7);
    hinvMass->Fit("BEexpol", "REBMS");

    double *obtained_parameters = BEexpol->GetParameters();
    TF1 *expol = new TF1("expol", exponential_bkg, 1, 2.2, 3);
    expol->SetParameter(0, obtained_parameters[9]);
    expol->SetParameter(1, obtained_parameters[10]);
    expol->SetParameter(2, obtained_parameters[11]);
    expol->SetLineColor(3);
    expol->SetLineStyle(2);
    expol->Draw("same");

    TF1 *onlyBW = new TF1("onlyBW", choerentBW_fitfunction_wo_bkg2, 1.02, 2.2, 9);
    for (int i = 0; i < 9; i++)
    {
        // onlyBW->SetParameter(i, parameters[i]);
        onlyBW->SetParameter(i, BEexpol->GetParameter(i));
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
}
