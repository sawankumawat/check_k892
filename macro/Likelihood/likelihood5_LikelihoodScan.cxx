#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"

//// We have calculated the best parameters values the minimizes -ln(L), now we want to calculate the uncertainties on these parameters. We can do this by looking at the change in the log-likelihood function as we vary each parameter while keeping the others fixed at their best fit values. The 1-sigma uncertainty corresponds to the parameter values for which the log-likelihood increases by 0.5 from its minimum value.
//// For fast calculation we will only use two free paramters L(Ns, Nb)

using namespace std;

//--------------------------------------------------
// Signal shape
//--------------------------------------------------

double gaussian(double x, double mean, double sigma)
{
    return TMath::Gaus(x, mean, sigma, true);
}

//--------------------------------------------------
// Background shape
//--------------------------------------------------

double background(double x)
{
    return (0.1 + 0.02 * x) / 2.0; // Factor 2.0 is for normalization over the range 0-10
}

//--------------------------------------------------
// Total model
//--------------------------------------------------

double totalFunction(double *x, double *par)
{
    double xx = x[0];

    double Ns = par[0];
    double Nb = par[1];

    double mean = 5.0;
    double sigma = 0.3;

    double s = gaussian(xx, mean, sigma);
    double b = background(xx);

    // expected counts per bin
    double bw = 0.1;

    return (Ns * s + Nb * b) * bw;
}

//--------------------------------------------------
// Main function
//--------------------------------------------------

void likelihood5_LikelihoodScan()
{
    //--------------------------------------------------
    // Create toy histogram
    //--------------------------------------------------

    TH1D *h = new TH1D("h", "Toy data", 100, 0, 10);

    TRandom3 r(0);

    //--------------------------------------------------
    // Generate signal
    //--------------------------------------------------

    for (int i = 0; i < 1000; i++)
    {
        double x = r.Gaus(5.0, 0.3);
        h->Fill(x);
    }

    //--------------------------------------------------
    // Generate background
    //--------------------------------------------------

    TF1 *fbkg = new TF1("fbkg", "0.1 + 0.02*x", 0, 10);

    for (int i = 0; i < 9000; i++)
    {
        double x = fbkg->GetRandom();
        h->Fill(x);
    }

    //--------------------------------------------------
    // Fixed shape parameters
    //--------------------------------------------------

    double mean = 5.0;
    double sigma = 0.3;

    //--------------------------------------------------
    // Variables to store best fit
    //--------------------------------------------------

    double bestLogL = -1e30;

    double bestNs = 0;
    double bestNb = 0;

    //--------------------------------------------------
    // Scan parameter space
    //--------------------------------------------------

    for (double Ns = 950; Ns <= 1030; Ns += 1)
    {
        for (double Nb = 8850; Nb <= 10000; Nb += 1)
        {
            double logL = 0;

            //--------------------------------------------------
            // Loop over histogram bins
            //--------------------------------------------------

            for (int i = 1; i <= h->GetNbinsX(); i++)
            {
                double x = h->GetBinCenter(i);
                double bw = h->GetBinWidth(i);

                // observed counts
                double n = h->GetBinContent(i);

                //--------------------------------------------------
                // Expected counts
                //--------------------------------------------------

                double s = gaussian(x, mean, sigma);

                double b = background(x);

                double mu = (Ns * s + Nb * b) * bw;

                //--------------------------------------------------
                // Poisson log likelihood
                //--------------------------------------------------

                if (mu > 0)
                {
                    logL += n * log(mu) - mu;
                }
            }

            //--------------------------------------------------
            // Check if likelihood improved
            //--------------------------------------------------

            if (logL > bestLogL)
            {
                bestLogL = logL;

                bestNs = Ns;
                bestNb = Nb;
            }
        }
    }

    //--------------------------------------------------
    // Print best fit
    //--------------------------------------------------

    cout << "Best fit results" << endl;
    cout << "Signal yield     = " << bestNs << endl;
    cout << "Background yield = " << bestNb << endl;
    cout << "Maximum logL     = " << bestLogL << endl;

    //--------------------------------------------------
    // Draw histogram and the best fit
    //--------------------------------------------------
    TCanvas *cLikelihoodFit = new TCanvas("cLikelihoodFit", "cLikelihoodFit", 720, 720);
    h->SetLineColor(kBlack);
    h->Draw();
    TF1 *f = new TF1("f", totalFunction, 0, 10, 2);
    f->SetParameter(0, bestNs);
    f->SetParameter(1, bestNb);
    f->SetLineStyle(2);
    f->Draw("same");

    //--------------------------------------------------
    // To compute the significance we compare the signal+background model with background only model.
    //--------------------------------------------------

    //--------------------------------------------------
    // BACKGROUND ONLY FIT
    //--------------------------------------------------

    double bestLogL_B = -1e30;
    double bestNb_B = 0;

    for (double Nb = 8850; Nb <= 10000; Nb += 1)
    {
        double logL = 0;

        for (int i = 1; i <= h->GetNbinsX(); i++)
        {
            double x = h->GetBinCenter(i);
            double bw = h->GetBinWidth(i);

            double n = h->GetBinContent(i);

            //--------------------------------------
            // signal fixed to ZERO
            //--------------------------------------

            double b = background(x);

            double mu = (Nb * b) * bw;

            if (mu > 0)
            {
                logL += n * log(mu) - mu;
            }
        }

        //--------------------------------------
        // keep best background-only likelihood
        //--------------------------------------

        if (logL > bestLogL_B)
        {
            bestLogL_B = logL;
            bestNb_B = Nb;
        }
    }

    //--------------------------------------------------
    // LIKELIHOOD RATIO (only significance estimation)
    //--------------------------------------------------

    double q = 2.0 * (bestLogL - bestLogL_B);
    double significance = sqrt(q);

    cout << endl;
    cout << "===== SIGNIFICANCE =====" << endl;
    cout << "Best logL (S+B) = " << bestLogL << endl;
    cout << "Best logL (B)   = " << bestLogL_B << endl;
    cout << "q = " << q << endl;
    cout << "Significance = " << significance << " sigma" << endl;

    //------------------------------------------------
    //===============================================
    // Likelihood scan for error estimation
    //================================================
    //-------------------------------------------------

    TH1D *hScan = new TH1D("hScan", "-2#Delta logL scan;N_{s};-2#Delta logL", 260, 900, 1160);

    //--------------------------------------------------
    // Variables for uncertainty extraction
    //--------------------------------------------------

    double prevNs = 0;
    double prevDelta = 0;

    bool foundLeft = false;
    bool foundRight = false;

    double leftCrossing = 0;
    double rightCrossing = 0;

    //--------------------------------------------------
    // Scan likelihood for parameter Ns
    //--------------------------------------------------

    for (double Ns = 900; Ns <= 1160; Ns += 1)
    {
        double logL = 0;

        //--------------------------------------------------
        // Compute likelihood for this Ns
        //--------------------------------------------------

        for (int i = 1; i <= h->GetNbinsX(); i++)
        {
            double x = h->GetBinCenter(i);

            double bw = h->GetBinWidth(i);

            double n = h->GetBinContent(i);

            //------------------------------------------
            // Fixed background
            //------------------------------------------

            double s = gaussian(x, mean, sigma);

            double b = background(x);

            double mu = (Ns * s + bestNb * b) * bw; // Ns free, Nb fixed to best fit value

            if (mu > 0)
            {
                logL += n * log(mu) - mu;
            }
        }

        //--------------------------------------------------
        // Delta(-2logL)
        //--------------------------------------------------

        double deltaNLL = -2.0 * (logL - bestLogL);

        hScan->Fill(Ns, deltaNLL);

        //--------------------------------------------------
        // LEFT crossing
        //--------------------------------------------------

        if (Ns < bestNs)
        {
            if (prevDelta > 1.0 && deltaNLL <= 1.0 && !foundLeft) // Left side which goes from higher to lower LogL values (since bestLogL is in middle and we start from lower Ns values). so we check crossing from gt 1 to le 1 using the previous Delta and current Delta values.
            {
                double fraction = (1.0 - prevDelta) / (deltaNLL - prevDelta); // linear interpolation (x-x1)/ (x2-x1) = (y-y1)/(y2-y1)

                leftCrossing = prevNs + fraction * (Ns - prevNs);

                foundLeft = true;
            }
        }

        //--------------------------------------------------
        // RIGHT crossing
        //--------------------------------------------------

        if (Ns > bestNs)
        {
            if (prevDelta < 1.0 && deltaNLL >= 1.0 && !foundRight) // Right side which goes from best LogL value i.e low value to the higher value as we increase Ns. So we check crossing from lt 1 to ge 1
            {
                double fraction = (1.0 - prevDelta) / (deltaNLL - prevDelta);

                rightCrossing = prevNs + fraction * (Ns - prevNs);

                foundRight = true;
            }
        }

        //--------------------------------------------------
        // Store previous point
        //--------------------------------------------------

        prevNs = Ns;
        prevDelta = deltaNLL;
    }

    TCanvas *cLikelihoodScanNs = new TCanvas("cLikelihoodScanNs", "cLikelihoodScanNs", 720, 720);
    hScan->SetLineColor(kBlue);
    hScan->SetMarkerStyle(20);
    hScan->SetMinimum(-0.1);
    hScan->Draw("HIST");
    TLine *line = new TLine(900, 1.0, 1160, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw("same");
    TLine *lineBest = new TLine(bestNs, 0, bestNs, hScan->GetMaximum());
    lineBest->SetLineStyle(2);
    lineBest->SetLineColor(kGreen + 2);
    lineBest->Draw("same");

    //--------------------------------------------------
    // Print uncertainties
    //--------------------------------------------------

    double errLow = bestNs - leftCrossing;
    double errHigh = rightCrossing - bestNs;

    cout << endl;
    cout << "===== LIKELIHOOD SCAN UNCERTAINTY =====" << endl;
    cout << "Best Ns = " << bestNs << endl;
    cout << "Left crossing  = " << leftCrossing << endl;
    cout << "Right crossing = " << rightCrossing << endl;
    cout << endl;
    cout << "Ns = " << bestNs << " + " << errHigh << " - " << errLow << endl;

    //// Compare these error with Root fit with option "L"
    TF1 *fFit = new TF1("fFit", totalFunction, 0, 10, 2);
    fFit->SetParameter(0, bestNs);
    fFit->SetParameter(1, bestNb);
    h->Fit(fFit, "L0");
    double errNsFit = fFit->GetParError(0);
    cout << "Ns from fit = " << fFit->GetParameter(0) << " +/- " << errNsFit << endl;
    
    //// Fit with additional option "E" to get the error estimation from profile likelihood from root instead of assuming parabolic shape at minimum and using hessian matrix (manual calculation for profile likelihood is shown in next file likelihood6_profileLikelihoodRatio.cxx)
    
    TF1 *fFitE = new TF1("fFitE", totalFunction, 0, 10, 2);
    fFitE->SetParameter(0, bestNs);
    fFitE->SetParameter(1, bestNb);
    TFitResultPtr fitResultE =h->Fit(fFitE, "LES0");
    double errNsFitLow = fitResultE->LowerError(0);
    double errNsFitHigh = fitResultE->UpperError(0);
    cout << "Ns from fit with option E = " << fFitE->GetParameter(0) << " + " << errNsFitHigh << " - " << -errNsFitLow << endl;
}