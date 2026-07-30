#include <iostream>
#include <cmath>
#include "TH2D.h"
#include "TF2.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"

// Global pointer to access the data histogram inside the Minuit FCN
TH2D *hData = nullptr;

// Definition of the Real Spherical Harmonics used for the waves
// x[0] = cos(theta), x[1] = phi
// par[0] = t00, par[1] = t20, par[2] = t22
double Intensity(double *x, double *par)
{
    double costheta = x[0];
    double phi = x[1];

    double t00 = par[0];
    double t20 = par[1];
    double t22 = par[2];

    double sintheta = sqrt(1.0 - costheta * costheta);

    // Calculate Real Spherical Harmonics Re(Y_L^M)
    // L=0, M=0
    double y00 = 1.0 / sqrt(4.0 * TMath::Pi());
    // L=2, M=0
    double y20 = sqrt(5.0 / (16.0 * TMath::Pi())) * (3.0 * costheta * costheta - 1.0);
    // L=2, M=2
    double y22 = sqrt(15.0 / (32.0 * TMath::Pi())) * sintheta * sintheta * cos(2.0 * phi);

    // n_LM factors: n00=1, n20=1, n22=2
    double I0 = 1.0 * t00 * y00 +
                1.0 * t20 * y20 +
                2.0 * t22 * y22;

    // Intensity must be physically positive
    return (I0 > 0.0) ? I0 : 1e-12;
}

// Binned Negative Log-Likelihood Function for Minuit
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
    double nll = 0.0;

    int nBinsX = hData->GetNbinsX();
    int nBinsY = hData->GetNbinsY();

    double dCosTheta = hData->GetXaxis()->GetBinWidth(1);
    double dPhi = hData->GetYaxis()->GetBinWidth(1);
    double dOmega = dCosTheta * dPhi; // Solid angle per bin

    for (int ix = 1; ix <= nBinsX; ++ix)
    {
        for (int iy = 1; iy <= nBinsY; ++iy)
        {

            double N_data = hData->GetBinContent(ix, iy);

            // Get bin centers
            double x[2] = {hData->GetXaxis()->GetBinCenter(ix),
                           hData->GetYaxis()->GetBinCenter(iy)};

            // Calculate expected counts (mu) in this bin
            double I0 = Intensity(x, par);
            double mu = I0 * dOmega; // epsilon = 1

            // Accumulate Poisson negative log-likelihood: -ln(P) = mu - N * ln(mu)
            if (mu > 0)
            {
                if (N_data > 0)
                {
                    nll += mu - N_data * log(mu);
                }
                else
                {
                    nll += mu; // N_data = 0
                }
            }
            else
            {
                if (N_data > 0)
                    nll += 1e9; // Heavy penalty for non-physical parameters
            }
        }
    }
    f = nll; // Pass the negative log-likelihood back to Minuit
}

// Main Macro
void Angular2Dfit3()
{
    gStyle->SetOptStat(0);

    // 1. Define True Parameters
    double true_t00 = 100000.0; // Needs to be large to keep baseline intensity positive
    double true_t20 = 20000.0;
    double true_t22 = 15000.0;

    // 2. Setup the theoretical intensity function
    TF2 *fInt = new TF2("fInt", Intensity, -1.0, 1.0, -TMath::Pi(), TMath::Pi(), 3);
    fInt->SetParameters(true_t00, true_t20, true_t22);
    fInt->SetNpx(100);
    fInt->SetNpy(100);

    // // 3. Generate Toy Data (Binned)
    // int nEvents = 500000;
    // hData = new TH2D("hData", "Toy Data (#cos#theta vs #phi);cos(#theta);#phi",
    //                  50, -1.0, 1.0, 50, -TMath::Pi(), TMath::Pi());

    // std::cout << "Generating " << nEvents << " toy events..." << std::endl;
    // hData->FillRandom("fInt", nEvents);

    // 3. Generate Toy Data (Binned) using Poisson statistics based on exact parameters
    hData = new TH2D("hData", "Toy Data (#cos#theta vs #phi);cos(#theta);#phi",
                     50, -1.0, 1.0, 50, -TMath::Pi(), TMath::Pi());

    std::cout << "Generating toy events from exact parameter scale..." << std::endl;

    double dCosTheta = hData->GetXaxis()->GetBinWidth(1);
    double dPhi = hData->GetYaxis()->GetBinWidth(1);
    double dOmega = dCosTheta * dPhi; // Solid angle per bin

    TRandom3 randGen(0); // Random number generator for Poisson sampling
    int totalGeneratedEvents = 0;

    for (int ix = 1; ix <= hData->GetNbinsX(); ++ix)
    {
        for (int iy = 1; iy <= hData->GetNbinsY(); ++iy)
        {

            // Get bin center
            double x[2] = {hData->GetXaxis()->GetBinCenter(ix),
                           hData->GetYaxis()->GetBinCenter(iy)};

            // Calculate exact expected mu for this bin
            double mu = Intensity(x, fInt->GetParameters()) * dOmega;

            // Draw an observed integer number of events from a Poisson distribution
            int n_obs = 0;
            if (mu > 0)
            {
                n_obs = randGen.Poisson(mu);
            }

            hData->SetBinContent(ix, iy, n_obs);
            totalGeneratedEvents += n_obs;
        }
    }

    std::cout << "Total generated events: " << totalGeneratedEvents << " (Expected ~35449)" << std::endl;

    // 4. Setup TMinuit for fitting
    TMinuit *minuit = new TMinuit(3); // 3 parameters
    minuit->SetFCN(fcn);

    // 5. Initialize fit parameters (purposely offset from true values)
    int ierflg = 0;
    minuit->mnparm(0, "t00", 80000.0, 10.0, 0, 0, ierflg);
    minuit->mnparm(1, "t20", 20000.0, 10.0, 0, 0, ierflg);
    minuit->mnparm(2, "t22", 5000.0, 10.0, 0, 0, ierflg);

    // 6. Run the Minimizer
    double arglist[10];
    arglist[0] = 0.5; // Error definition for log-likelihood (0.5 for -lnL)
    minuit->mnexcm("SET ERR", arglist, 1, ierflg);

    arglist[0] = 5000; // max calls
    arglist[1] = 1.0;  // tolerance

    std::cout << "\nStarting Minuit Fit..." << std::endl;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // 7. Extract and Print Results
    double fit_t00, err_t00, fit_t20, err_t20, fit_t22, err_t22;
    minuit->GetParameter(0, fit_t00, err_t00);
    minuit->GetParameter(1, fit_t20, err_t20);
    minuit->GetParameter(2, fit_t22, err_t22);

    std::cout << "\n============================================\n";
    std::cout << "                 FIT RESULTS                  \n";
    std::cout << "============================================\n";
    std::cout << "Parameter | Input Value | Fit Result +/- Error\n";
    std::cout << "--------------------------------------------\n";
    std::cout << Form("t00       | %10.1f | %8.1f +/- %5.1f\n", true_t00, fit_t00, err_t00);
    std::cout << Form("t20       | %10.1f | %8.1f +/- %5.1f\n", true_t20, fit_t20, err_t20);
    std::cout << Form("t22       | %10.1f | %8.1f +/- %5.1f\n", true_t22, fit_t22, err_t22);
    std::cout << "============================================\n";

    // 8. Draw the toy data histogram with the fitted intensity overlaid
    // Create a wrapper/scaled TF2 to match histogram bin counts on the Z axis
    double dOmega2 = hData->GetXaxis()->GetBinWidth(1) * hData->GetYaxis()->GetBinWidth(1);

    // Standard intensity TF2 for 2D contour (cd 1)
    TF2 *fFit2D = (TF2 *)fInt->Clone("fFit2D");
    fFit2D->SetParameters(fit_t00, fit_t20, fit_t22);
    fFit2D->SetLineColor(kRed + 1);
    fFit2D->SetLineWidth(2);

    // Scaled intensity TF2 for 3D surface matching histogram counts (cd 2)
    // Intensity * dOmega yields expected bin contents (mu)
    TF2 *fFit3D = new TF2("fFit3D", [fit_t00, fit_t20, fit_t22, dOmega2](double *x, double *p)
                          {
    double par[3] = {fit_t00, fit_t20, fit_t22};
    return Intensity(x, par) * dOmega2; }, -1.0, 1.0, -TMath::Pi(), TMath::Pi(), 0);

    fFit3D->SetLineColor(kRed);
    fFit3D->SetFillColorAlpha(kRed, 0.35); // Semi-transparent red surface

    TCanvas *c1 = new TCanvas("c1", "PWA Binned Toy Fits", 1400, 600);
    c1->Divide(2, 1);

    // Left Pad: 2D Color plot with contour overlay
    c1->cd(1);
    gPad->SetRightMargin(0.15);
    hData->Draw("COLZ");
    fFit2D->Draw("CONT3 SAME");

    // Right Pad: 3D LEGO plot with 3D Fitted Surface overlay
    c1->cd(2);
    gPad->SetRightMargin(0.15);
    hData->Draw("LEGO2");
    fFit3D->Draw("SURF SAME"); // Overlays fit surface directly in 3D view

    c1->Update();
}