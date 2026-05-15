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

void likelihood2_minimize_nonExtended()
{
    //--------------------------------------------------
    // Create toy histogram
    //--------------------------------------------------

    TH1D *h = new TH1D("h", "Toy data", 100, 0, 10);

    TRandom3 r(0);

    //--------------------------------------------------
    // Generate signal
    //--------------------------------------------------

    for (int i = 0; i < 2000; i++)
    {
        double x = r.Gaus(5.0, 0.3);
        h->Fill(x);
    }

    //--------------------------------------------------
    // Generate background
    //--------------------------------------------------

    TF1 *fbkg = new TF1("fbkg", "0.1 + 0.02*x", 0, 10);

    for (int i = 0; i < 8000; i++)
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

    double bestfs = 0;

    //--------------------------------------------------
    // Scan parameter space
    //--------------------------------------------------

    for (double fs = 0.00; fs <= 0.5; fs += 0.01)
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

            double mu = (fs * s + (1 - fs) * b) * bw;

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

            bestfs = fs;
        }
    }

    //--------------------------------------------------
    // Print best fit
    //--------------------------------------------------

    cout << "Best fit results" << endl;

    cout << "Signal fraction  = " << bestfs << endl;

    cout << "Minimum -logL     = " << bestLogL << endl;

    //--------------------------------------------------
    // Draw histogram and the best fit
    //--------------------------------------------------
    TCanvas *c = new TCanvas("c", "c", 720, 720);
    h->SetLineColor(kBlack);
    h->Draw();
    double entries = h->Integral();
    TF1 *f = new TF1("f", totalFunction, 0, 10, 2);
    f->SetParameter(0, bestfs * entries);
    f->SetParameter(1, (1 - bestfs) * entries);
    f->SetLineStyle(2);
    f->Draw("same");

}