#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"

//// In this code we are doing the same thing as in likelihood1_oneParam.cxx but now we are scanning the parameter guess Ns and Nb to find the values that maximize the likelihood. Then we are fitting the model with the best parameters we got. We are also comparing our results with the ROOT fit option "L" which does the same thing and also provides uncertainties on the parameters.

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

void likelihood2_minimize()
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

    // for (int i = 0; i < 9000; i++)
    // {
    //     double x = r.Uniform(0, 10);

    //     double y = r.Uniform(0, 0.3);

    //     if (y < background(x))
    //         h->Fill(x);
    // }

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

    for (double Ns = 900; Ns <= 1100; Ns += 1)
    {
        for (double Nb = 8000; Nb <= 9500; Nb += 1)
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

    cout << "Minimum -logL     = " << bestLogL << endl;

    //--------------------------------------------------
    // Draw histogram and the best fit
    //--------------------------------------------------
    h->SetLineColor(kBlack);
    h->Draw();
    TF1 *f = new TF1("f", totalFunction, 0, 10, 2);
    f->SetParameter(0, bestNs);
    f->SetParameter(1, bestNb);
    f->SetLineStyle(2);
    f->Draw("same");

    //--------------------------------------------------
    // Compare the results with ROOT L fit option which does the same thing and also provides uncertainties on the parameters.
    //--------------------------------------------------

    TF1 *fROOT = new TF1("fROOT", totalFunction, 0, 10, 2);
    fROOT->SetParameter(0, 1000);
    fROOT->SetParameter(1, 9000);
    TFitResultPtr fitResult = h->Fit(fROOT, "LS0");

    cout << endl;
    cout << "===== ROOT OPTION L FIT =====" << endl;

    cout << "Signal     = " << fROOT->GetParameter(0) << endl;
    cout << "Background = " << fROOT->GetParameter(1) << endl;
    cout << "Minimum -logL = " << fitResult->MinFcnValue() << endl; // ROOT returns -logL (with additional factorial factor that we ingnored as it is constant. That factor does not matter for minimization).
}