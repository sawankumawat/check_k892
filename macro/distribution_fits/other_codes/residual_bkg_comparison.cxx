#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/style.h"
using namespace std;

Double_t exponential_bkg_3(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[3])));
}

Double_t expol_chkstar(double *x, double *par)
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] + x[0] * par[3]));
}

Double_t Boltzmann_bkg_1(double *x, double *par) // 3 parameters
{
    double norm = par[0];
    double massks = 0.497;
    double n = par[1];
    double C = par[2];
    double y = norm * pow((x[0] - 2 * massks), n / 2) * pow(C, 3 / 2) * exp(-C * pow((x[0] - 2 * massks), n));
    return y;
}

void residual_bkg_comparison()
{
    gStyle->SetOptStat(0);
    TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60_fullpt.root", "READ"); // full pT range
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", 0.0, 30.0));
    // if (pT_bins[0] != 0.0 || pT_bins[1] != 30.0)
    // {
    //     cout << "Not full pT range" << endl;
    //     return;
    // }
    if (hinvMass == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
    hinvMass->Rebin(2);
    SetHistoQA(hinvMass);
    hinvMass->SetMarkerSize(0.62);
    hinvMass->GetXaxis()->SetRangeUser(1.00, 2.50);
    hinvMass->Draw();
    float fitrange_low = 1.02;
    float fitrange_high = 2.20;
    // TF1 *expol = new TF1("expol", exponential_bkg_3, fitrange_low, fitrange_high, 4);
    // // expol->SetParameters(617042.324, -0.069, 2.709, 1.027); // 1.02 - 2.20
    // // expol->SetParameters(1375426.548, 0.147, 3.544, 0.806); // 1.05 - 2.20
    // expol->SetParameters(904279.119, 0.047, 3.135, 0.941); // 1.05 - 2.30
    // // expol->SetParameters(5.62894e+05, -9.27599e-02, 2.60563e+00, 1.07645e+00); // 1.00 - 2.20
    // // expol->SetParameters(1.34241e+06, 1.58099e-01, 3.49088e+00, 8.40530e-01);  // 1.10 - 2.20
    // // expol->SetParameters(5.72804e+05, -8.62894e-02, 2.64582e+00, 1.09212e+00); // 1.02 - 2.30
    // // expol->SetParameters(5.38386e+05, -1.01490e-01, 2.60553e+00, 1.15508e+00); // 1.02 - 2.40
    // // expol->SetParameters(5.13488e+05, -1.12616e-01, 2.68930e+00, 1.22087e+00); // 1.02 - 2.50
    // expol->SetLineWidth(2);
    // // expol->SetLineStyle(2);
    // expol->Draw("same");

    TF1 *boltzmann = new TF1("boltzmann", Boltzmann_bkg_1, fitrange_low, fitrange_high, 3);
    boltzmann->SetParameters(770718.278, 0.592, 4.494); // 1.02 - 2.20
    // boltzmann->SetParameters(708544.710, 0.665, 4.406); // 1.05 - 2.20
    // boltzmann->SetParameters(690831.678, 0.685, 4.536); // 1.05 - 2.30
    // boltzmann->SetParameters(8.43199e+05, 5.21042e-01, 4.55467e+00); // 1.00 - 2.20
    // boltzmann->SetParameters(6.78135e+05, 7.17350e-01, 4.27840e+00); // 1.10 - 2.20
    // boltzmann->SetParameters(7.57609e+05, 6.05557e-01, 4.59917e+00); // 1.02 - 2.30
    // boltzmann->SetParameters(7.39720e+05, 6.23074e-01, 4.75354e+00); // 1.02 - 2.40
    // boltzmann->SetParameters(7.20030e+05, 6.42877e-01, 4.93774e+00); // 1.02 - 2.50
    boltzmann->SetLineColor(4);
    boltzmann->SetLineWidth(2);
    // boltzmann->SetLineStyle(3);
    boltzmann->Draw("same");

    TF1 *expol_pol1 = new TF1("expol_pol1", expol_chkstar, fitrange_low, fitrange_high, 4);
    expol_pol1->SetParameters(385.461, -0.059, 10.158, -2.758); // 1.02 - 2.20
    // expol_pol1->SetParameters(440.224, -0.007, 10.301, -2.882); // 1.05 - 2.20
    // expol_pol1->SetParameters(460.899, 0.004, 10.356, -2.952); // 1.05 - 2.30
    // expol_pol1->SetParameters(3.65735e+02, -7.73484e-02, 1.01066e+01, -2.71538e+00); // 1.00 - 2.20
    // expol_pol1->SetParameters(4.09960e+02, -3.82645e-02, 1.02214e+01, -2.81221e+00); // 1.10 - 2.20
    // expol_pol1->SetParameters(4.01343e+02, -5.58932e-02, 1.01986e+01, -2.82765e+00); // 1.02 - 2.30
    // expol_pol1->SetParameters(4.40388e+02, -4.88387e-02, 1.03002e+01, -2.99961e+00); // 1.02 - 2.40
    // expol_pol1->SetParameters(5.01609e+02, -3.62839e-02, 1.04606e+01, -3.24896e+00); // 1.02 - 2.50
    expol_pol1->SetLineColor(2);
    expol_pol1->SetLineWidth(2);
    // expol_pol1->SetLineStyle(4);
    expol_pol1->Draw("same");

    TLegend *leg = new TLegend(0.45, 0.67, 0.9, 0.92);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->AddEntry((TObject *)0, Form("Fit range: %.2f - %.2f GeV/c^{2}", fitrange_low, fitrange_high), "");
    leg->AddEntry(hinvMass, "Data", "lpe");
    // leg->AddEntry(expol, "Exponential function", "l");
    leg->AddEntry(boltzmann, "Boltzmann function", "l");
    leg->AddEntry(expol_pol1, "Expol 1 function", "l");
    leg->Draw("same");
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/fits/3rBW_fits/";
    c->SaveAs(Form("%s/res_bkg_comp_%.2f_%.2f.png", path.c_str(), fitrange_low, fitrange_high));
}