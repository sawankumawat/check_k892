#include <iostream>
using namespace std;
#include "style.h"
// #include "../macro/src/fitfunc.h"
#include "../macro/src/common_glue.h"

Double_t RelativisticBW(double *x, double *par)
{
    double norm = par[0];
    double mass = par[1];
    double width = par[2];

    Int_t j1 = 0; // spin
    double n1 = (2.0 * j1 + 1.0) / 2.0;

    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = mass * mass - 4 * (0.4976 * 0.4976);

    double mass_dep_width = width * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = norm * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));

    return fit;
}

Double_t RBW_massDepWidth(double *x, double *par)
{
    double norm = par[0];
    double mass = par[1];
    double width = par[2];

    Int_t j1 = 0; // spin
    double n1 = (2.0 * j1 + 1.0) / 2.0;

    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = mass * mass - 4 * (0.4976 * 0.4976);

    double mass_dep_width = width * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = norm * mass * mass_dep_width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * mass_dep_width, 2));

    return fit;
}

void BW_fit()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1111);

    TFile *f = new TFile("/home/sawan/check_k892/mc/LHC24l1/337948.root", "read");

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    TH1F *recMass = (TH1F *)f->Get("higher-mass-resonances_a01320/hMChists/Recf1710_mass");
    if (recMass == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.15, 0.05, 0.05, 0.15);
    SetHistostyle2(recMass);
    recMass->SetTitle("Rec f_{0}(1710) Mass");
    recMass->GetYaxis()->SetTitle("Counts");
    // recMass->GetXaxis()->SetRangeUser(1.49, 2.01);
    // recMass->GetXaxis()->SetRangeUser(1.19, 1.89);
    recMass->GetXaxis()->SetRangeUser(1.09, 1.69);
    recMass->Draw("pe");

    // TF1 *fit = new TF1("fit", RelativisticBW, 1.5, 1.9, 3);
    // TF1 *fit = new TF1("fit", RelativisticBW, 1.2, 1.8, 3);
    TF1 *fit = new TF1("fit", RelativisticBW, 1.0, 1.6, 3);
    fit->SetParNames("Norm", "Mass", "Width");
    fit->SetParameter(0, 100);
    fit->SetParameter(1, a1320Mass);
    fit->SetParameter(2, a1320Width);
    fit->SetParLimits(0, 0.0, 200.0);
    fit->SetParLimits(2, 0.0, 2.0);
    recMass->Fit("fit", "REBMS0");
    // fit->SetLineColor(4);
    fit->Draw("same"); 

    // TPaveStats *st = (TPaveStats *)fit->FindObject("stats");
    // st->SetX1NDC(0.6);
    // st->SetX2NDC(0.99);
    // st->SetY1NDC(0.4);
    // st->SetY2NDC(0.9);
    // st->Draw();

    // // TF1 *fit2 = new TF1("fit2", RBW_massDepWidth, 1.5, 1.9, 3);
    // TF1 *fit2 = new TF1("fit2", RBW_massDepWidth, 1.2, 1.8, 3);
    // fit2->SetParameter(0, 1);
    // fit2->SetParameter(1, f1525Mass);
    // fit2->SetParameter(2, f1525Width);
    // recMass->Fit("fit2", "REBMS0");
    // // fit2->Draw("same");

    c->SaveAs("plots/BW_fit_a21320.png");
}