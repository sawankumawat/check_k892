#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"

//// Now we calculate the profile likelihood ratio.
//// Earlier in likelihood scan we calculate the likelihood each time by varying the paramter of interest, but in profile likelihood ratio we do not keep the other parameters fixed at their best fit values, instead we minimize the likelihood with respect to the other parameters for each value of the parameter of interest. This is more accurate way to calculate the uncertainties as it takes into account the correlations between the parameters.

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

void likelihood6_profileLikelihoodRatio()
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

    //------------------------------------------------
    //===============================================
    // Profile likelihood ratios
    //================================================
    //------------------------------------------------

    TGraph *gProfile = new TGraph();

    int p = 0;

    // scan fixed Ns values
    for (double Ns_fixed = 900; Ns_fixed <= 1150; Ns_fixed += 1)
    {
        //------------------------------------------
        // For this fixed Ns:
        // maximize likelihood wrt Nb
        //------------------------------------------

        double bestProfileLogL = -1e30;
        double bestProfileNb = 0;

        for (double Nb = 8500; Nb <= 10000; Nb += 1)
        {
            double logL = 0;

            //--------------------------------------
            // Loop over bins
            //--------------------------------------

            for (int i = 1; i <= h->GetNbinsX(); i++)
            {
                double x = h->GetBinCenter(i);
                double bw = h->GetBinWidth(i);

                double n = h->GetBinContent(i);

                double s = gaussian(x, mean, sigma);
                double b = background(x);

                //--------------------------------------
                // Ns fixed
                //--------------------------------------

                double mu = (Ns_fixed * s + Nb * b) * bw;

                if (mu > 0)
                {
                    logL += n * log(mu) - mu;
                }
            }

            //--------------------------------------
            // keep best Nb for this Ns
            //--------------------------------------

            if (logL > bestProfileLogL)
            {
                bestProfileLogL = logL;
                bestProfileNb = Nb;
            }
        }

        //------------------------------------------
        // profile likelihood ratio
        //------------------------------------------

        double q = 2.0 * (bestLogL - bestProfileLogL);

        gProfile->SetPoint(p, Ns_fixed, q);

        p++;

        // cout << "Ns = " << Ns_fixed << "  profiled Nb = " << bestProfileNb << "  q = " << q << endl;
    }

    TCanvas *cProfile = new TCanvas("cProfile", "Profile Likelihood", 700, 600);

    gProfile->SetTitle("Profile Likelihood Scan;N_{s};-2ln#lambda");
    gProfile->SetLineWidth(2);
    gProfile->Draw("AL");
    TLine *line1Sigma = new TLine(900, 1, 1150, 1);
    line1Sigma->SetLineColor(kRed);
    line1Sigma->SetLineStyle(2);
    line1Sigma->Draw("same");
    TLine *lineBestFit = new TLine(bestNs, 0, bestNs, 5);
    lineBestFit->SetLineColor(kBlue);
    lineBestFit->SetLineStyle(2);
    lineBestFit->Draw("same");
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry(line1Sigma, "1#sigma", "l");
    leg->AddEntry(lineBestFit, "Best N_{s} value", "l");
    leg->Draw();
}