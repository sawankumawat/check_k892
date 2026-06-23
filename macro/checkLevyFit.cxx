#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>

Double_t FuncLavy(Double_t *x, Double_t *par)
{
    Double_t pT = x[0];

    Double_t n = par[0];
    Double_t dNdy = par[1];
    Double_t mass = par[2];
    Double_t T = par[3];

    Double_t mT = sqrt(mass * mass + pT * pT);

    Double_t p =
        (n - 1.) * (n - 2.) * dNdy * pT /
        (pow(1. + (mT - mass) / (n * T), n) *
         (n * T * (n * T + mass * (n - 2.))));

    return p;
}

void checkLevyFit()
{
    gStyle->SetOptStat(0);

    // Default parameters
    double n0 = 7.0;
    double dNdy0 = 0.06;
    double mass0 = 0.895;
    double T0 = 0.2;

    TCanvas *c1 = new TCanvas("c1", "Levy-Tsallis comparison", 1440, 720);
    c1->Divide(2, 1);

    // =====================================================
    // Vary n
    // =====================================================
    c1->cd(1);
    gPad->SetLogy(1);
    TF1 *fDef = new TF1("fDef", FuncLavy, 0., 10., 4);
    fDef->SetParameters(n0, dNdy0, mass0, T0);
    fDef->SetLineWidth(3);
    fDef->SetLineColor(kBlack);
    fDef->SetTitle("Effect of n;p_{T} (GeV/c);dN/(dp_{T}dy)");

    TF1 *fNlow = new TF1("fNlow", FuncLavy, 0., 10., 4);
    fNlow->SetParameters(4.0, dNdy0, mass0, T0);
    fNlow->SetLineColor(kRed);

    TF1 *fNhigh = new TF1("fNhigh", FuncLavy, 0., 10., 4);
    fNhigh->SetParameters(10.0, dNdy0, mass0, T0);
    fNhigh->SetLineColor(kBlue);

    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.02);

    fDef->Draw();
    fNlow->Draw("same");
    fNhigh->Draw("same");

    // TLegend *leg1 = new TLegend(0.55, 0.70, 0.88, 0.88);
    // leg1->AddEntry(fDef, Form("Default n = %.2f", n0), "l");
    // leg1->AddEntry(fNlow, "n = 4.0", "l");
    // leg1->AddEntry(fNhigh, "n = 10.0", "l");
    // leg1->Draw();

    TLegend *leg1 = new TLegend(0.35, 0.65, 0.90, 0.88);
    leg1->SetTextSize(0.04);
    leg1->AddEntry(fDef, Form("n=%.0f, <p_{T}>=%.3f GeV/c", n0, fDef->Mean(0, 10)), "l");
    leg1->AddEntry(fNlow, Form("n=4, <p_{T}>=%.3f GeV/c", fNlow->Mean(0, 10)), "l");
    leg1->AddEntry(fNhigh, Form("n=10, <p_{T}>=%.3f GeV/c", fNhigh->Mean(0, 10)), "l");
    leg1->Draw();

    // =====================================================
    // Vary T
    // =====================================================
    c1->cd(2);
    gPad->SetLogy(1);
    TF1 *fTdef = new TF1("fTdef", FuncLavy, 0., 10., 4);
    fTdef->SetParameters(n0, dNdy0, mass0, T0);
    fTdef->SetLineWidth(3);
    fTdef->SetLineColor(kBlack);
    fTdef->SetTitle("Effect of T;p_{T} (GeV/c);dN/(dp_{T}dy)");

    TF1 *fTlow = new TF1("fTlow", FuncLavy, 0., 10., 4);
    fTlow->SetParameters(n0, dNdy0, mass0, 0.15);
    fTlow->SetLineColor(kRed);

    TF1 *fThigh = new TF1("fThigh", FuncLavy, 0., 10., 4);
    fThigh->SetParameters(n0, dNdy0, mass0, 0.30);
    fThigh->SetLineColor(kBlue);

    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.02);

    fTdef->Draw();
    fTlow->Draw("same");
    fThigh->Draw("same");

    TLegend *leg2 = new TLegend(0.4, 0.70, 0.88, 0.88);
    // leg2->AddEntry(fTdef, Form("Default T = %.3f", T0), "l");
    // leg2->AddEntry(fTlow, "T = 0.15 GeV", "l");
    // leg2->AddEntry(fThigh, "T = 0.30 GeV", "l");
    // leg2->Draw();
    leg2->SetTextSize(0.04);

    leg2->AddEntry(fTdef, Form("T=%.3f, <p_{T}>=%.3f", T0, fTdef->Mean(0, 10)), "l");
    leg2->AddEntry(fTlow, Form("T=0.15, <p_{T}>=%.3f", fTlow->Mean(0, 10)), "l");
    leg2->AddEntry(fThigh, Form("T=0.30, <p_{T}>=%.3f", fThigh->Mean(0, 10)), "l");
    leg2->Draw();

    c1->Update();

    c1->SaveAs("LevyFitComparison.png");
}