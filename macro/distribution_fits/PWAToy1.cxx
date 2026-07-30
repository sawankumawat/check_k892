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
    double t21 = par[3];

    double sintheta = sqrt(1.0 - costheta * costheta);

    // Calculate Real Spherical Harmonics Re(Y_L^M)
    // L=0, M=0
    double y00 = 1.0 / sqrt(4.0 * TMath::Pi());
    // L=2, M=0
    double y20 = sqrt(5.0 / (16.0 * TMath::Pi())) * (3.0 * costheta * costheta - 1.0);
    // L=2, M=2
    double y22 = sqrt(15.0 / (32.0 * TMath::Pi())) * sintheta * sintheta * cos(2.0 * phi);
    // L=2, M=1
    double y21 = -sqrt(15.0 / (32.0 * TMath::Pi())) * cos(phi) * (2 * costheta * sqrt(1 - costheta * costheta));

    // n_LM factors: n00=1, n20=1, n22=2
    double I0 = 1.0 * t00 * y00 +
                1.0 * t20 * y20 +
                2.0 * t22 * y22 +
                2.0 * t21 * y21; // Include t21 term if needed

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
                    nll += mu - N_data * log(mu); // can also use log(I0) since the log(ϵj*dΩj) term is constant and can be ignored in minimization
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
void PWAToy1()
{
    gStyle->SetOptStat(0);

    // 1. Define True Parameters
    double true_t00 = 100000.0; // Needs to be large to keep baseline intensity positive
    double true_t20 = 20000.0;
    double true_t22 = 15000.0;
    double true_t21 = 7000.0; // Define the fourth parameter

    // 2. Setup the theoretical intensity function
    TF2 *fInt = new TF2("fInt", Intensity, -1.0, 1.0, -TMath::Pi(), TMath::Pi(), 4);
    fInt->SetParameters(true_t00, true_t20, true_t22, true_t21);
    fInt->SetNpx(100);
    fInt->SetNpy(100);

    // // 3. Generate Toy Data (Binned)
    // int nEvents = 500000;
    // hData = new TH2D("hData", "Toy Data (#cos#theta vs #phi);cos(#theta);#phi",
    //                  50, -1.0, 1.0, 50, -TMath::Pi(), TMath::Pi());

    // std::cout << "Generating " << nEvents << " toy events..." << std::endl;
    // hData->FillRandom("fInt", nEvents);

    // 3. Generate Toy Data (Binned) using Poisson statistics based on exact parameters
    hData = new TH2D("hData", "Toy Data (#cos#theta vs #phi);cos(#theta);#phi", 50, -1.0, 1.0, 50, -TMath::Pi(), TMath::Pi());

    std::cout << "Generating toy events from exact parameter scale..." << std::endl;

    double dCosTheta = hData->GetXaxis()->GetBinWidth(1);
    double dPhi = hData->GetYaxis()->GetBinWidth(1);
    double dOmega = dCosTheta * dPhi; // Solid angle per bin

    TRandom3 randGen(0); // Random number generator for Poisson sampling
    int totalGeneratedEvents = 0;

    for (int ix = 1; ix <= hData->GetNbinsX(); ++ix) // Loop over cos(theta) bins
    {
        for (int iy = 1; iy <= hData->GetNbinsY(); ++iy) // Loop over phi bins
        {

            // Get bin center
            double x[2] = {hData->GetXaxis()->GetBinCenter(ix),
                           hData->GetYaxis()->GetBinCenter(iy)};

            // Calculate exact expected mu (average counts) in this bin
            double mu = Intensity(x, fInt->GetParameters()) * dOmega; //efficiency = 1

            // Draw an observed integer number of events from a Poisson distribution
            int n_obs = 0;
            if (mu > 0)
            {
                n_obs = randGen.Poisson(mu);
            }

            hData->SetBinContent(ix, iy, n_obs);
            totalGeneratedEvents += n_obs;

            // Here the events are generated based on the expected intensity (function of partial amptudes tLM) using Poisson statistics and then binned into the 2D histogram hData.
        }
    }

    std::cout << "Total generated events: " << totalGeneratedEvents << " (Expected ~35449)" << std::endl;

    // 4. Setup TMinuit for fitting
    TMinuit *minuit = new TMinuit(4); // 4 parameters
    minuit->SetFCN(fcn);              // Set the function to minimize (using negative log-likelihood)

    // 5. Initialize fit parameters (purposely offset from true values)
    int ierflg = 0; // Integer error flag for Minuit. If ierflg != 0, an error occurred during the parameter initialization or minimization.
    // minuit->mnparm(parameter_index, parameter_name, initial_value, step_size, lower_limit, upper_limit, ierflg);
    minuit->mnparm(0, "t00", 80000.0, 10.0, 0, 0, ierflg); // Upper and lower limit 0 i.e. there is no limit on parameter
    minuit->mnparm(1, "t20", 20000.0, 10.0, 0, 0, ierflg);
    minuit->mnparm(2, "t22", 5000.0, 10.0, 0, 0, ierflg);
    minuit->mnparm(3, "t21", 9000.0, 10.0, 0, 0, ierflg);

    if (ierflg != 0)
    {
        std::cerr << "Error creating parameter. ierflg = " << ierflg << std::endl;
    }

    // 6. Run the Minimizer
    double arglist[10];
    arglist[0] = 0.5;                              // In NLL, 0.5 corresponds to 1-sigma error, just like value 1 in chi-square minimization.
    minuit->mnexcm("SET ERR", arglist, 1, ierflg); // 1 is telling that 1 argument is passed to the command

    arglist[0] = 5000; // max calls for function evaluation
    arglist[1] = 1.0;  // tolerance

    std::cout << "\nStarting Minuit Fit..." << std::endl;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg); // MIGRAD with 2 arguments: max calls and tolerance. It is a minimization algorithm that uses the gradient of the function to find the minimum.

    // 7. Extract and Print Results
    double fit_t00, err_t00, fit_t20, err_t20, fit_t22, err_t22, fit_t21, err_t21;
    minuit->GetParameter(0, fit_t00, err_t00);
    minuit->GetParameter(1, fit_t20, err_t20);
    minuit->GetParameter(2, fit_t22, err_t22);
    minuit->GetParameter(3, fit_t21, err_t21);

    // --- Calculate Chi2 and NDF ---
    double chi2 = 0.0;
    int nBinsUsed = 0;
    double fitPars[4] = {fit_t00, fit_t20, fit_t22, fit_t21};

    int nBinsX = hData->GetNbinsX();
    int nBinsY = hData->GetNbinsY();
    double dOmega_val = hData->GetXaxis()->GetBinWidth(1) * hData->GetYaxis()->GetBinWidth(1);

    for (int ix = 1; ix <= nBinsX; ++ix)
    {
        for (int iy = 1; iy <= nBinsY; ++iy)
        {
            double N_data = hData->GetBinContent(ix, iy);

            double x[2] = {hData->GetXaxis()->GetBinCenter(ix),
                           hData->GetYaxis()->GetBinCenter(iy)};

            // Expected model count in this bin
            double mu = Intensity(x, fitPars) * dOmega_val; // Multiply by efficiency here if efficiency is used

            // Using Pearson Chi2 (variance = mu)
            if (mu > 0)
            {
                double diff = N_data - mu;
                chi2 += (diff * diff) / mu;
                nBinsUsed++;
            }
        }
    }

    int nFitParams = 4; // t00, t20, t22, t21
    int ndf = nBinsUsed - nFitParams;
    double chi2_ndf = chi2 / ndf;

    std::cout << "\n============================================\n";
    std::cout << "                 FIT RESULTS                  \n";
    std::cout << "============================================\n";
    std::cout << "Parameter | Input Value | Fit Result +/- Error\n";
    std::cout << "--------------------------------------------\n";
    std::cout << Form("t00       | %10.1f | %8.1f +/- %5.1f\n", true_t00, fit_t00, err_t00);
    std::cout << Form("t20       | %10.1f | %8.1f +/- %5.1f\n", true_t20, fit_t20, err_t20);
    std::cout << Form("t22       | %10.1f | %8.1f +/- %5.1f\n", true_t22, fit_t22, err_t22);
    std::cout << Form("t21       | %10.1f | %8.1f +/- %5.1f\n", true_t21, fit_t21, err_t21);
    std::cout << "--------------------------------------------\n";
    std::cout << Form("Chi2 / NDF : %.2f / %d = %.3f\n", chi2, ndf, chi2_ndf);
    std::cout << Form("p-value    : %.4e\n", TMath::Prob(chi2, ndf));
    std::cout << "============================================\n";

    // 8. Visualizations
    double dOmega2 = hData->GetXaxis()->GetBinWidth(1) * hData->GetYaxis()->GetBinWidth(1);

    // --- 2D and 3D Fit Functions ---
    TF2 *fFit2D = (TF2 *)fInt->Clone("fFit2D");
    fFit2D->SetParameters(fit_t00, fit_t20, fit_t22, fit_t21);
    fFit2D->SetLineColor(kRed + 1);
    fFit2D->SetLineWidth(2);

    // Lambda function which captures the fitted parameters and dOmega for 3D surface plot (its longer version which we usually write is writte in the end)
    TF2 *fFit3D = new TF2("fFit3D", [fit_t00, fit_t20, fit_t22, fit_t21, dOmega2](double *x, double *p)
                          {
        double par[4] = {fit_t00, fit_t20, fit_t22, fit_t21};
        return Intensity(x, par) * dOmega2; }, -1.0, 1.0, -TMath::Pi(), TMath::Pi(), 0);
    fFit3D->SetLineColor(kRed);
    fFit3D->SetFillColorAlpha(kRed, 0.35);

    // --- 1D Projections and Fit Models ---
    // Extract 1D data projections
    TH1D *hData_cosTheta = hData->ProjectionX("hData_cosTheta");
    hData_cosTheta->SetTitle("cos(#theta) Projection;cos(#theta);Events");
    hData_cosTheta->SetMarkerStyle(20);
    hData_cosTheta->SetMarkerSize(0.8);

    TH1D *hData_phi = hData->ProjectionY("hData_phi");
    hData_phi->SetTitle("#phi Projection;#phi;Events");
    hData_phi->SetMarkerStyle(20);
    hData_phi->SetMarkerSize(0.8);

    // Create 1D histograms for fit projections
    TH1D *hFit_cosTheta = (TH1D *)hData_cosTheta->Clone("hFit_cosTheta");
    hFit_cosTheta->Reset();
    hFit_cosTheta->SetLineColor(kRed + 1);
    hFit_cosTheta->SetLineWidth(2);

    TH1D *hFit_phi = (TH1D *)hData_phi->Clone("hFit_phi");
    hFit_phi->Reset();
    hFit_phi->SetLineColor(kRed + 1);
    hFit_phi->SetLineWidth(2);

    // Populate 1D fit histograms by integrating/summing predicted counts across 2D bins
    // int nBinsX = hData->GetNbinsX();
    // int nBinsY = hData->GetNbinsY();
    // double fitPars[4] = {fit_t00, fit_t20, fit_t22, fit_t21};

    for (int ix = 1; ix <= nBinsX; ++ix)
    {
        for (int iy = 1; iy <= nBinsY; ++iy)
        {
            double x[2] = {hData->GetXaxis()->GetBinCenter(ix),
                           hData->GetYaxis()->GetBinCenter(iy)};

            // Expected count in 2D bin (ix, iy)
            double mu = Intensity(x, fitPars) * dOmega2;

            // Accumulate expected counts into 1D projections
            hFit_cosTheta->AddBinContent(ix, mu);
            hFit_phi->AddBinContent(iy, mu);
        }
    }

    // --- Draw Everything on a 2x2 Canvas ---
    TCanvas *c1 = new TCanvas("c1", "PWA Binned Fits and 1D Projections", 1200, 1000);
    c1->Divide(2, 2);

    // Pad 1: 2D Color Plot with Contour
    c1->cd(1);
    gPad->SetRightMargin(0.15);
    hData->Draw("COLZ");
    fFit2D->Draw("CONT3 SAME");

    // Pad 2: 3D Surface Overlay
    c1->cd(2);
    gPad->SetRightMargin(0.15);
    hData->Draw("LEGO2");
    fFit3D->Draw("SURF SAME");

    // Pad 3: 1D Projection in cos(theta)
    c1->cd(3);
    gPad->SetLeftMargin(0.12);
    hData_cosTheta->Draw("E");        // Data points with Poisson errors
    hFit_cosTheta->Draw("HIST SAME"); // Fit model smooth curve

    // Pad 4: 1D Projection in phi
    c1->cd(4);
    gPad->SetLeftMargin(0.12);
    hData_phi->Draw("E");        // Data points with Poisson errors
    hFit_phi->Draw("HIST SAME"); // Fit model smooth curve

    c1->Update();
}

// // The lambda function is equivalent to the below:
// double g_t00;
// double g_t20;
// double g_t22;
// double g_t21;
// double g_dOmega2;

// double MyFunction(double *x, double *p)
// {
//     double par[4];

//     par[0] = g_t00;
//     par[1] = g_t20;
//     par[2] = g_t22;
//     par[3] = g_t21;

//     return Intensity(x, par) * g_dOmega2;
// }
// TF2 *fFit3D = new TF2("fFit3D", MyFunction, -1, 1, -TMath::Pi(), TMath::Pi(), 0);