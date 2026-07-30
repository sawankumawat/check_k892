#include <iostream>
#include <cmath>
#include "TH2D.h"
#include "TF2.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TStyle.h"

//==============================================================
// In this code we additionally include the efficiency term
//==============================================================

// Global pointers for data and efficiency histograms
TH2D *hData = nullptr;
TH2D *hEff = nullptr; // Efficiency map: epsilon_j = N_rec / N_gen

// Toy efficiency function for simulation: e.g., lower efficiency at extreme cos(theta)
double EfficiencyFunc(double costheta, double phi)
{
    return 0.8 * (1.0 - 0.5 * costheta * costheta);
}

// Differential intensity I_0(Omega)
double Intensity(double *x, double *par)
{
    double costheta = x[0];
    double phi = x[1];

    double t00 = par[0];
    double t20 = par[1];
    double t22 = par[2];

    double sintheta = sqrt(1.0 - costheta * costheta);

    double y00 = 1.0 / sqrt(4.0 * TMath::Pi());
    double y20 = sqrt(5.0 / (16.0 * TMath::Pi())) * (3.0 * costheta * costheta - 1.0);
    double y22 = sqrt(15.0 / (32.0 * TMath::Pi())) * sintheta * sintheta * cos(2.0 * phi);

    double I0 = 1.0 * t00 * y00 +
                1.0 * t20 * y20 +
                2.0 * t22 * y22;

    return (I0 > 0.0) ? I0 : 1e-12;
}

// Binned Negative Log-Likelihood Function with Efficiency Included
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
            double eff = hEff->GetBinContent(ix, iy); // Retrieve epsilon_j

            double x[2] = {hData->GetXaxis()->GetBinCenter(ix),
                           hData->GetYaxis()->GetBinCenter(iy)};

            // Expected count in bin j: mu_j = I_0(Omega_j) * epsilon_j * Delta(Omega)
            double I0 = Intensity(x, par);
            double mu = I0 * eff * dOmega;

            // Binned Poisson Likelihood
            if (mu > 0)
            {
                if (N_data > 0)
                {
                    nll += mu - N_data * log(mu);
                }
                else
                {
                    nll += mu;
                }
            }
            else
            {
                if (N_data > 0)
                    nll += 1e9; // Penalty
            }
        }
    }
    f = nll;
}

void PWAToy2()
{
    gStyle->SetOptStat(0);

    // 1. True parameters
    double true_t00 = 10000.0;
    double true_t20 = 2000.0;
    double true_t22 = 1500.0;

    TF2 *fInt = new TF2("fInt", Intensity, -1.0, 1.0, -TMath::Pi(), TMath::Pi(), 3);
    fInt->SetParameters(true_t00, true_t20, true_t22);

    int nBinsX = 50;
    int nBinsY = 50;

    // 2. Build 2D Efficiency Map from Toy MC
    TH2D *hMCGen = new TH2D("hMCGen", "Generated MC", nBinsX, -1.0, 1.0, nBinsY, -TMath::Pi(), TMath::Pi());
    TH2D *hMCRec = new TH2D("hMCRec", "Reconstructed MC", nBinsX, -1.0, 1.0, nBinsY, -TMath::Pi(), TMath::Pi());
    hEff = new TH2D("hEff", "Efficiency Map", nBinsX, -1.0, 1.0, nBinsY, -TMath::Pi(), TMath::Pi());

    // Populate MC histograms assuming flat generated phase space (like DRgen in the analysis note)
    int nMCGenEventsPerBin = 10000; // High statistics for precise efficiency
    for (int ix = 1; ix <= nBinsX; ++ix)
    {
        for (int iy = 1; iy <= nBinsY; ++iy)
        {
            double costheta = hMCGen->GetXaxis()->GetBinCenter(ix);
            double phi = hMCGen->GetYaxis()->GetBinCenter(iy);

            double eff_val = EfficiencyFunc(costheta, phi);

            hMCGen->SetBinContent(ix, iy, nMCGenEventsPerBin);
            hMCRec->SetBinContent(ix, iy, nMCGenEventsPerBin * eff_val);
        }
    }

    // Divide reconstructed by generated to get efficiency epsilon_j
    hEff->Divide(hMCRec, hMCGen);

    // 3. Generate Toy Data taking efficiency into account
    hData = new TH2D("hData", "Accepted Toy Data;cos(#theta);#phi",
                     nBinsX, -1.0, 1.0, nBinsY, -TMath::Pi(), TMath::Pi());

    double dCosTheta = hData->GetXaxis()->GetBinWidth(1);
    double dPhi = hData->GetYaxis()->GetBinWidth(1);
    double dOmega = dCosTheta * dPhi;

    TRandom3 randGen(0);
    int totalEvents = 0;

    for (int ix = 1; ix <= nBinsX; ++ix)
    {
        for (int iy = 1; iy <= nBinsY; ++iy)
        {
            double x[2] = {hData->GetXaxis()->GetBinCenter(ix),
                           hData->GetYaxis()->GetBinCenter(iy)};

            double eff_j = hEff->GetBinContent(ix, iy);

            // Expected count WITH efficiency folding
            double mu = Intensity(x, fInt->GetParameters()) * eff_j * dOmega;

            int n_obs = (mu > 0) ? randGen.Poisson(mu) : 0;
            hData->SetBinContent(ix, iy, n_obs);
            totalEvents += n_obs;
        }
    }
    std::cout << "Generated total observed data events (after efficiency): " << totalEvents << std::endl;

    // 4. Setup Minuit Fit
    TMinuit *minuit = new TMinuit(3);
    minuit->SetFCN(fcn);

    int ierflg = 0;
    minuit->mnparm(0, "t00", 8000.0, 10.0, 0, 0, ierflg);
    minuit->mnparm(1, "t20", 500.0, 10.0, 0, 0, ierflg);
    minuit->mnparm(2, "t22", 500.0, 10.0, 0, 0, ierflg);

    double arglist[10];
    arglist[0] = 0.5; // Error definition for -lnL
    minuit->mnexcm("SET ERR", arglist, 1, ierflg);

    arglist[0] = 5000;
    arglist[1] = 1.0;

    std::cout << "\nStarting Minuit Fit with Efficiency Correction..." << std::endl;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // 5. Results
    double fit_t00, err_t00, fit_t20, err_t20, fit_t22, err_t22;
    minuit->GetParameter(0, fit_t00, err_t00);
    minuit->GetParameter(1, fit_t20, err_t20);
    minuit->GetParameter(2, fit_t22, err_t22);

    std::cout << "\n============================================\n";
    std::cout << "         EFFICIENCY-CORRECTED FIT RESULTS     \n";
    std::cout << "============================================\n";
    std::cout << "Parameter | Input Value | Fit Result +/- Error\n";
    std::cout << "--------------------------------------------\n";
    std::cout << Form("t00       | %10.1f | %8.1f +/- %5.1f\n", true_t00, fit_t00, err_t00);
    std::cout << Form("t20       | %10.1f | %8.1f +/- %5.1f\n", true_t20, fit_t20, err_t20);
    std::cout << Form("t22       | %10.1f | %8.1f +/- %5.1f\n", true_t22, fit_t22, err_t22);
    std::cout << "============================================\n";

    // 6. Draw plots
    TCanvas *c1 = new TCanvas("c1", "PWA Efficiency Fit", 1200, 400);
    c1->Divide(3, 1);
    c1->cd(1);
    hData->Draw("COLZ");
    c1->cd(2);
    hEff->Draw("COLZ");
    c1->cd(3);
    TH2D *hFitModel = (TH2D *)hData->Clone("hFitModel");
    hFitModel->SetTitle("Fit Model Expected (#mu_{j});cos(#theta);#phi");
    for (int ix = 1; ix <= nBinsX; ++ix)
    {
        for (int iy = 1; iy <= nBinsY; ++iy)
        {
            double x[2] = {hFitModel->GetXaxis()->GetBinCenter(ix), hFitModel->GetYaxis()->GetBinCenter(iy)};
            double pars[3] = {fit_t00, fit_t20, fit_t22};
            double mu = Intensity(x, pars) * hEff->GetBinContent(ix, iy) * dOmega;
            hFitModel->SetBinContent(ix, iy, mu);
        }
    }
    hFitModel->Draw("COLZ");

    // Additional fit-overlay plots similar to PWAToy1.cxx
    // double dOmega2 = dOmega; // expected bin solid angle
    double dOmega2 = hData->GetXaxis()->GetBinWidth(1) * hData->GetYaxis()->GetBinWidth(1);

    // 2D intensity contour (fitted parameters)
    TF2 *fFit2D = (TF2 *)fInt->Clone("fFit2D");
    fFit2D->SetParameters(fit_t00, fit_t20, fit_t22);
    fFit2D->SetLineColor(kRed + 1);
    fFit2D->SetLineWidth(2);

    // 3D fitted surface including efficiency (uses toy EfficiencyFunc)
    TF2 *fFit3D = new TF2("fFit3D", [fit_t00, fit_t20, fit_t22, dOmega2](double *x, double *p)
                          {
                              double par[3] = {fit_t00, fit_t20, fit_t22};
                              double costheta = x[0];
                              double phi = x[1];
                              double mu = Intensity(x, par) * EfficiencyFunc(costheta, phi) * dOmega2;
                              return mu; }, -1.0, 1.0, -TMath::Pi(), TMath::Pi(), 0);

    fFit3D->SetLineColor(kRed);
    fFit3D->SetFillColorAlpha(kRed, 0.35);

    TCanvas *c2 = new TCanvas("c2", "PWA Efficiency Fit Overlays", 1200, 500);
    c2->Divide(2, 1);
    c2->cd(1);
    gPad->SetRightMargin(0.15);
    hData->Draw("COLZ");
    fFit2D->Draw("CONT3 SAME");

    c2->cd(2);
    gPad->SetRightMargin(0.15);
    hData->Draw("LEGO2");
    fFit3D->Draw("SURF SAME");

    c2->Update();
}