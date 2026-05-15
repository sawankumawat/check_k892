#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"

using namespace std;

// In this code we create a histogram with Gaussian + pol1 background. After then we calculate the likelihood function for the signal+background model. To just start with, we are not doing any fitting and then minimizing the parameters (Ns and Nb) to maximize the likelihood. We know their exact values and are just calculating the likelihood for those values. This is just to understand how to calculate the likelihood function and then we will move to the next step of fitting and maximizing the likelihood in the likelihood2_minimize.cxx code.

//// This is extended likelihood function in which we use the total number of events in each bin using the poisson distribution instead of considering fixed number of events in each bin. This is more realistic case from which we can also get the number of signal and background events along with the shape.

double gaussian(double x, double mean, double sigma)
{
    return TMath::Gaus(x, mean, sigma, true); // Already normalized Gaussian
}

double background(double x)
{
    // simple linear background
    return (0.1 + 0.02 * x) / 2.0; // Factor 2.0 is for normalization over the range 0-10
}

void likelihood1_oneParam()
{
    //------------------------------------------
    // Create toy histogram
    //------------------------------------------

    TH1D *h = new TH1D("h", "toy data", 100, 0, 10);

    TRandom3 r(0);

    //------------------------------------------
    // Generate signal events
    //------------------------------------------

    for (int i = 0; i < 1000; i++)
    {
        double x = r.Gaus(5.0, 0.3);
        h->Fill(x);
    }

    //------------------------------------------
    // Generate background events
    //------------------------------------------

    // for (int i = 0; i < 9000; i++)
    // {
    //     double x = r.Uniform(0, 10);

    //     // accept/reject for linear background
    //     double y = r.Uniform(0, 0.3);

    //     if (y < background(x))
    //         h->Fill(x);
    // }

    // // We can also fill the background using the TF1 random number generator, which is more efficient and also ensures that we get the correct number of background events (9000 in this case) without fluctuations, because the above logic can reject some events leading to reduced background statistics.

    TF1 *fbkg = new TF1("fbkg", "0.1 + 0.02*x", 0, 10);

    for (int i = 0; i < 9000; i++)
    {
        double x = fbkg->GetRandom();
        h->Fill(x);
    }

    //------------------------------------------
    // Model parameters
    //------------------------------------------

    // // The highest value of likelihood will be obtained for the true values of Ns and Nb, which are 1000 and 9000 respectively in this case. We can vary these parameters and see how the likelihood changes. It should decrease as we move away from the true values.

    double Ns = 1000; // expected signal yield
    double Nb = 9000; // expected background yield

    double mean = 5.0;
    double sigma = 0.3;

    //------------------------------------------
    // Calculate log-likelihood
    //------------------------------------------

    double logL = 0;

    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        double x = h->GetBinCenter(i);
        double bw = h->GetBinWidth(i);

        // observed counts
        double n = h->GetBinContent(i);

        //--------------------------------------
        // expected counts in this bin
        //--------------------------------------

        double s = gaussian(x, mean, sigma);

        double b = background(x);

        double mu = (Ns * s + Nb * b) * bw;

        //--------------------------------------
        // Poisson log likelihood
        //--------------------------------------

        if (mu > 0)
        {
            logL += n * log(mu) - mu; // Ignoring the n! term from poisson function which is constant (does not depend on the parameters). n is the observed counts (probability of observing n counts when the mean (expected counts) is mu).
        }
    }

    //------------------------------------------
    // Print result
    //------------------------------------------
    cout << "Ns = " << Ns << ", Nb = " << Nb << endl;
    cout << "Log Likelihood = " << logL << endl;

    // h->Draw();

    //// The above was extended likelihood function. Now let us see the non-extended likelihood function method.
    //// In non-extended likelihood function, we do not fit total normalization. Total number of observed events is fixed and only relative fraction of signal/background matters. So instead of L(Ns, Nb) we will have L(fs) where fs is signal fraction and background will be (1-fs).
    //// Extended Likelihood: L = Poisson(Nobs | Nexp) * Product_over_bins(PDF)
    //// Non-extended Likelihood: L = Product_over_bins(PDF)

    double fs = 0.1; // 10% signal fraction

    //------------------------------------------
    // Calculate shape-only log likelihood
    //------------------------------------------

    double logL_nonExtended = 0;

    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        double x = h->GetBinCenter(i);

        // observed counts
        double n = h->GetBinContent(i);

        //--------------------------------------
        // normalized PDFs
        //--------------------------------------

        double s = gaussian(x, mean, sigma);

        double b = background(x);

        //--------------------------------------
        // total normalized model
        //--------------------------------------

        double p = fs * s + (1.0 - fs) * b;

        //--------------------------------------
        // shape-only likelihood
        //--------------------------------------

        if (p > 0)
        {
            logL_nonExtended += n * log(p); // Here n is fixed which is total number of events in a bin.
        }
    }
    cout << "\nSignal fraction = " << fs << endl;
    cout << "Non-extended Log Likelihood = " << logL_nonExtended << endl; // This value does not matter. Only the change in it matters. (dont compare with the extended likelihood since we excluded factorial term in it)
}