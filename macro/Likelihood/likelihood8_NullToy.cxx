#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"

//// In this code we will calculate the p-value by generating the toy events using null model and checking which toy events have a test statistic (q) equal to or more than the observed value. This is a more direct way to compute the p-value without relying on the asymptotic approximation of the test statistic distribution.

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

void likelihood8_NullToy()
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
    // LIKELIHOOD RATIO SIGNIFICANCE
    //--------------------------------------------------

    double q = 2.0 * (bestLogL - bestLogL_B);

    double significance = sqrt(q);

    //--------------------------------------------------
    // p-value from Gaussian tail
    //--------------------------------------------------

    double pvalue = 0.5 * TMath::Erfc(significance / sqrt(2.0)); // Two-sided p-value for Gaussian distribution

    cout << endl;
    cout << "===== SIGNIFICANCE =====" << endl;
    cout << "Best logL (S+B) = " << bestLogL << endl;
    cout << "Best logL (B)   = " << bestLogL_B << endl;
    cout << "q = " << q << endl;
    cout << "Significance = " << significance << " sigma" << endl;
    cout << "p-value = " << pvalue << endl;

    //--------------------------------------------------
    // TOY MONTE CARLO
    //--------------------------------------------------

    int ntoys = 1000; // For 5σ significance, we should atleast generate 10^7 events. Generating low events just for testing purpose.

    TH1D *hq = new TH1D("hq", "q distribution from toys", 100, 0, 50);

    int nExtreme = 0;

    //--------------------------------------------------
    // Loop over toys
    //--------------------------------------------------

    for (int itoy = 0; itoy < ntoys; itoy++)
    {
        //--------------------------------------------------
        // Create toy histogram
        //--------------------------------------------------

        TH1D *htoy = new TH1D(Form("htoy_%d", itoy), "", 100, 0, 10);

        //--------------------------------------------------
        // Generate BACKGROUND ONLY toy
        //--------------------------------------------------

        for (int i = 0; i < 9000; i++)
        {
            double x = fbkg->GetRandom();
            htoy->Fill(x); // instead of Gaussian+pol2, just generating for pol2.
        }

        //--------------------------------------------------
        // SIGNAL + BACKGROUND FIT
        //--------------------------------------------------

        double toyBestLogL = -1e30;

        for (double Ns = 0; Ns <= 1000; Ns += 10)
        {
            for (double Nb = 8500; Nb <= 9500; Nb += 10)
            {
                double logL = 0;

                for (int i = 1; i <= htoy->GetNbinsX(); i++)
                {
                    double x = htoy->GetBinCenter(i);

                    double bw = htoy->GetBinWidth(i);

                    double n = htoy->GetBinContent(i);

                    double s = gaussian(x, mean, sigma);

                    double b = background(x);

                    double mu = (Ns * s + Nb * b) * bw;

                    if (mu > 0)
                    {
                        logL += n * log(mu) - mu;
                    }
                }

                if (logL > toyBestLogL)
                {
                    toyBestLogL = logL;
                }
            }
        }

        //--------------------------------------------------
        // BACKGROUND ONLY FIT
        //--------------------------------------------------

        double toyBestLogL_B = -1e30;

        for (double Nb = 8500; Nb <= 9500; Nb += 10)
        {
            double logL = 0;

            for (int i = 1; i <= htoy->GetNbinsX(); i++)
            {
                double x = htoy->GetBinCenter(i);

                double bw = htoy->GetBinWidth(i);

                double n = htoy->GetBinContent(i);

                double b = background(x);

                double mu = (Nb * b) * bw;

                if (mu > 0)
                {
                    logL += n * log(mu) - mu;
                }
            }

            if (logL > toyBestLogL_B)
            {
                toyBestLogL_B = logL;
            }
        }

        //--------------------------------------------------
        // Compute toy q
        //--------------------------------------------------

        double qtoy = 2.0 * (toyBestLogL - toyBestLogL_B);

        hq->Fill(qtoy);

        //--------------------------------------------------
        // Count extreme toys
        //--------------------------------------------------

        if (qtoy >= q)
        {
            nExtreme++;
        }

        delete htoy;
    }

    //--------------------------------------------------
    // Toy MC p-value
    //--------------------------------------------------

    double pToy = (double)nExtreme / ntoys;

    cout << endl;

    cout << "===== TOY MC =====" << endl;

    cout << "Observed q = "
         << q << endl;

    cout << "Extreme toys = "
         << nExtreme << endl;

    cout << "Toy p-value = "
         << pToy << endl;

    //--------------------------------------------------
    // Draw q distribution
    //--------------------------------------------------

    TCanvas *c2 = new TCanvas("c2", "Toy q distribution", 720, 720);
    hq->Draw();
}