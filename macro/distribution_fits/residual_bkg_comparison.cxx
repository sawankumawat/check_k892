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
    hinvMass->SetMarkerSize(0.6);
    hinvMass->GetXaxis()->SetRangeUser(1.00, 2.50);
    hinvMass->Draw();
    TF1 *expol = new TF1("expol", exponential_bkg_3, 1.02, 2.30, 4);
    expol->SetParameters(5.56200e+05, -9.45594e-02, 2.56900e+00, 1.10242e+00);
    expol->SetLineWidth(2);
    // expol->SetLineStyle(2);
    expol->Draw("same");
    TF1 *boltzmann = new TF1("boltzmann", Boltzmann_bkg_1, 1.02, 2.30, 3);
    // boltzmann->SetParameters(7.54894e+05, 6.45600e-01, 4.23800e+00);
    boltzmann->SetParameters(7.70718e+05, 5.92481e-01, 4.49443e+00);
    boltzmann->SetLineColor(4);
    boltzmann->SetLineWidth(2);
    // boltzmann->SetLineStyle(3);
    boltzmann->Draw("same");

    TF1 *expol_pol1 = new TF1("expol_pol1", expol_chkstar, 1.02, 2.30, 4);
    expol_pol1->SetParameters(3.85461e+02, -5.94977e-02, 1.01575e+01, -2.75798e+00);
    expol_pol1->SetLineColor(8);
    expol_pol1->SetLineWidth(2);
    // expol_pol1->SetLineStyle(4);
    expol_pol1->Draw("same");

    TLegend *leg = new TLegend(0.55, 0.67, 0.9, 0.92);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->AddEntry(hinvMass, "Data", "lpe");
    leg->AddEntry(expol, "Exponential function", "l");
    leg->AddEntry(boltzmann, "Boltzmann function", "l");
    leg->AddEntry(expol_pol1, "Expol 1 function", "l");
    leg->Draw("same");

    c->SaveAs("/home/sawan/Music/residual_bkg_comparison.png");
}