#include <iostream>
#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "TF2.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"

// Global pointer to access data histogram inside Minuit FCN
TH2D *hData = nullptr;

// Real Spherical Harmonics Intensity Function
double Intensity(double *x, double *par)
{
    double costheta = x[0];
    double phi = x[1];

    double t00 = par[0];
    double t20 = par[1];
    double t22 = par[2];
    double t21 = par[3];

    double sintheta = sqrt(1.0 - costheta * costheta);

    // Real Spherical Harmonics Re(Y_L^M)
    double y00 = 1.0 / sqrt(4.0 * TMath::Pi());
    double y20 = sqrt(5.0 / (16.0 * TMath::Pi())) * (3.0 * costheta * costheta - 1.0);
    double y22 = sqrt(15.0 / (32.0 * TMath::Pi())) * sintheta * sintheta * cos(2.0 * phi);
    double y21 = -sqrt(15.0 / (32.0 * TMath::Pi())) * cos(phi) * (2.0 * costheta * sintheta);

    double I0 = 1.0 * t00 * y00 +
                1.0 * t20 * y20 +
                2.0 * t22 * y22 +
                2.0 * t21 * y21;

    return (I0 > 0.0) ? I0 : 1e-12;
}

// Binned NLL FCN for Minuit
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
    double nll = 0.0;
    int nBinsX = hData->GetNbinsX();
    int nBinsY = hData->GetNbinsY();

    double dCosTheta = hData->GetXaxis()->GetBinWidth(1);
    double dPhi = hData->GetYaxis()->GetBinWidth(1);
    double dOmega = dCosTheta * dPhi;

    for (int ix = 1; ix <= nBinsX; ++ix)
    {
        for (int iy = 1; iy <= nBinsY; ++iy)
        {
            double N_data = hData->GetBinContent(ix, iy);
            double x[2] = {hData->GetXaxis()->GetBinCenter(ix),
                           hData->GetYaxis()->GetBinCenter(iy)};

            double mu = Intensity(x, par) * dOmega;

            if (mu > 0.0)
            {
                if (N_data > 0.0)
                {
                    nll += mu - N_data * log(mu);
                }
                else
                {
                    nll += mu;
                }
            }
            else if (N_data > 0.0)
            {
                nll += 1e9;
            }
        }
    }
    f = nll;
}

// Mass-dependent parameter shapes for true input generation
double GetTrueT00(double mass) { return 100000.0 * exp(-0.5 * std::pow((mass - 1.5) / 0.3, 2)); }
double GetTrueT20(double mass) { return 20000.0 * sin(TMath::Pi() * (mass - 1.0)); }
double GetTrueT22(double mass) { return 15000.0 * (mass - 1.0); }
double GetTrueT21(double mass) { return 7000.0 * cos(TMath::Pi() * (mass - 1.0)); }

void PWAToyFullMassRange()
{
    gStyle->SetOptStat(0);

    // Binning setup: 40 MeV/c^2 mass bins from 1.0 to 2.0 GeV/c^2
    const double mMin = 1.0;
    const double mMax = 2.0;
    const double binWidth = 0.04;
    const int nMassBins = std::round((mMax - mMin) / binWidth);

    // Spectra Histograms
    TH1D *hGen = new TH1D("hGen", "Generated Yield;M(K^{0}_{s}K^{0}_{s}) [GeV/c^{2}];Events / 40 MeV/c^{2}", nMassBins, mMin, mMax);
    TH1D *hPW = new TH1D("hPW", "Reconstructed PW Yield;M(K^{0}_{s}K^{0}_{s}) [GeV/c^{2}];Events / 40 MeV/c^{2}", nMassBins, mMin, mMax);

    TRandom3 randGen(42);

    std::cout << "Starting ALICE-style PWA self-consistency check over " << nMassBins << " mass bins..." << std::endl;

    for (int im = 1; im <= nMassBins; ++im)
    {
        double mass = hGen->GetXaxis()->GetBinCenter(im);

        // 1. True input parameters at mass bin center
        double true_t00 = GetTrueT00(mass);
        double true_t20 = GetTrueT20(mass);
        double true_t22 = GetTrueT22(mass);
        double true_t21 = GetTrueT21(mass);
        double truePars[4] = {true_t00, true_t20, true_t22, true_t21};

        // 2. Generate Toy Data histogram for this mass bin
        if (hData)
            delete hData;
        hData = new TH2D(Form("hData_bin_%d", im), "Toy Data;cos(#theta);#phi",
                         50, -1.0, 1.0, 50, -TMath::Pi(), TMath::Pi());

        double dOmega = hData->GetXaxis()->GetBinWidth(1) * hData->GetYaxis()->GetBinWidth(1);
        int totalGenerated = 0;

        for (int ix = 1; ix <= hData->GetNbinsX(); ++ix)
        {
            for (int iy = 1; iy <= hData->GetNbinsY(); ++iy)
            {
                double x[2] = {hData->GetXaxis()->GetBinCenter(ix),
                               hData->GetYaxis()->GetBinCenter(iy)};
                double mu = Intensity(x, truePars) * dOmega;
                int n_obs = (mu > 0.0) ? randGen.Poisson(mu) : 0;
                hData->SetBinContent(ix, iy, n_obs);
                totalGenerated += n_obs;
            }
        }

        // Store true generated yield in histogram
        hGen->SetBinContent(im, totalGenerated);
        hGen->SetBinError(im, sqrt(totalGenerated));

        // 3. Minuit Fit to extract tLM moments
        TMinuit minuit(4);
        minuit.SetPrintLevel(-1);
        minuit.SetFCN(fcn);

        int ierflg = 0;
        minuit.mnparm(0, "t00", true_t00 * 0.8, 10.0, 0, 0, ierflg);
        minuit.mnparm(1, "t20", true_t20 * 0.8, 10.0, 0, 0, ierflg);
        minuit.mnparm(2, "t22", true_t22 * 0.8, 10.0, 0, 0, ierflg);
        minuit.mnparm(3, "t21", true_t21 * 0.8, 10.0, 0, 0, ierflg);

        double arglist[10];
        arglist[0] = 0.5;
        minuit.mnexcm("SET ERR", arglist, 1, ierflg);

        arglist[0] = 5000;
        arglist[1] = 1.0;
        minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

        // Get fitted parameters
        double fitPars[4];
        double fitErr[4];

        for (int i = 0; i < 4; i++)
        {
            minuit.GetParameter(i, fitPars[i], fitErr[i]);
        }

        // Reconstruct PW yield
        double PWyield = 0.0;

        for (int ix = 1; ix <= hData->GetNbinsX(); ++ix)
        {
            for (int iy = 1; iy <= hData->GetNbinsY(); ++iy)
            {
                double x[2] = {
                    hData->GetXaxis()->GetBinCenter(ix),
                    hData->GetYaxis()->GetBinCenter(iy)};

                PWyield += Intensity(x, fitPars) * dOmega;
            }
        }

        hPW->SetBinContent(im, PWyield);
        hPW->SetBinError(im, 0.0);
    }

    // 5. Compute PW / Gen Ratio
    TH1D *hRatio = (TH1D *)hPW->Clone("hRatio");
    hRatio->SetTitle("Self-Consistency Check;M(K^{0}_{s}K^{0}_{s}) [GeV/c^{2}];PW / Gen.");
    hRatio->SetTitle(0);
    hRatio->Divide(hGen);

    // 6. Draw Yields and Ratio Canvas (ALICE Style)
    TCanvas *cCheck = new TCanvas("cCheck", "PWA Mass Consistency Check", 1000, 800);
    cCheck->Divide(1, 2);

    // Top Pad: Yield Spectra
    cCheck->cd(1);
    gPad->SetPad(0.0, 0.35, 1.0, 1.0);
    gPad->SetBottomMargin(0.00);
    gPad->SetGrid();

    hGen->SetLineColor(kBlack);
    hGen->SetLineWidth(2);
    hGen->Draw("HIST E");

    hPW->SetMarkerStyle(20);
    hPW->SetMarkerSize(0.9);
    hPW->SetMarkerColor(kRed + 1);
    hPW->SetLineColor(kRed + 1);
    hPW->Draw("E1 SAME");

    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.85);
    leg->AddEntry(hGen, "Generated (Gen)", "l");
    leg->AddEntry(hPW, "Reconstructed PW", "p");
    leg->Draw();

    // Bottom Pad: Ratio Plot (PW / Gen)
    cCheck->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.35);
    gPad->SetTopMargin(0.00);
    gPad->SetBottomMargin(0.25);
    gPad->SetGrid();

    hRatio->SetMinimum(0.8);
    hRatio->SetMaximum(1.2);
    hRatio->GetYaxis()->SetTitleSize(0.08);
    hRatio->GetYaxis()->SetTitleOffset(0.5);
    hRatio->GetYaxis()->SetLabelSize(0.07);
    hRatio->GetXaxis()->SetTitleSize(0.09);
    hRatio->GetXaxis()->SetLabelSize(0.07);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.9);
    hRatio->SetMarkerColor(kBlue + 1);
    hRatio->SetLineColor(kBlue + 1);
    hRatio->Draw("E1");

    TLine *unity = new TLine(mMin, 1.0, mMax, 1.0);
    unity->SetLineStyle(2);
    unity->SetLineWidth(2);
    unity->SetLineColor(kBlack);
    unity->Draw();

    cCheck->Update();
}