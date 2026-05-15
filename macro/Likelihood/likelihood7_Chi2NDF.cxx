#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"

//// Now we calculate the profile likelihood ratio.
//// Earlier in likelihood scan we calculate the likelihood each time by varying the paramter of interest, but in profile likelihood ratio we do not keep the other parameters fixed at their best fit values, instead we minimize the likelihood with respect to the other parameters for each value of the parameter of interest. This is more accurate way to calculate the uncertainties as it takes into account the correlations between the parameters.

///// To understand profile likelihood ratio better for error estimation, we will keep the mean and sigma of the BW free

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
    double mean = par[2];
    double sigma = par[3];

    // double mean = 5.0;
    // double sigma = 0.3;

    double s = gaussian(xx, mean, sigma);
    double b = background(xx);

    // expected counts per bin
    double bw = 0.1;

    return (Ns * s + Nb * b) * bw;
}

//--------------------------------------------------
// Main function
//--------------------------------------------------

void likelihood7_Chi2NDF()
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

    // double mean = 5.0;
    // double sigma = 0.3;

    //--------------------------------------------------
    // Variables to store best fit
    //--------------------------------------------------

    double bestLogL = -1e30;

    double bestNs = 0;
    double bestNb = 0;
    double bestMean = 0;
    double bestSigma = 0;

    //--------------------------------------------------
    // Scan parameter space
    //--------------------------------------------------

    for (double Ns = 980; Ns <= 1030; Ns += 2)
    {
        for (double Nb = 8700; Nb <= 9200; Nb += 2)
        {
            for (double mean = 4.85; mean <= 5.15; mean += 0.01)
            {
                for (double sigma = 0.27; sigma <= 0.33; sigma += 0.01)
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
                        bestMean = mean;
                        bestSigma = sigma;
                    }
                }
            }
        }
    }

    //--------------------------------------------------
    // Print best fit
    //--------------------------------------------------

    cout << "Best fit results" << endl;
    cout << "Best signal yield     = " << bestNs << endl;
    cout << "Best background yield = " << bestNb << endl;
    cout << "Best mean        = " << bestMean << endl;
    cout << "Best sigma       = " << bestSigma << endl;
    cout << "Maximum logL     = " << bestLogL << endl;

    //--------------------------------------------------
    // Draw histogram and the best fit
    //--------------------------------------------------
    TCanvas *cLikelihoodFit = new TCanvas("cLikelihoodFit", "cLikelihoodFit", 720, 720);
    h->SetLineColor(kBlack);
    h->Draw();
    TF1 *f = new TF1("f", totalFunction, 0, 10, 4);
    f->SetParameter(0, bestNs);
    f->SetParameter(1, bestNb);
    f->SetParameter(2, bestMean);
    f->SetParameter(3, bestSigma);
    f->SetLineStyle(2);
    f->Draw("same");

    //// Now we will calculate the Chi2/NDF manually.
    //// The Chi2 is defined as sum over bins of (observed - expected)^2 / expected, also called Pearson's Chi2.
    //// In likelihood we calculate it as follow:
    //// We know from poisson ln(L) = summation(ni * ln(mu_i) - mu_i), where ni is the observed counts in bin i and mu_i is the mean (expected counts in bin i).
    //// In saturate model i.e. model which perfectly describes the data, "mu_i = n_i", i.e. expected counts are real counts. So ln(L_saturated) = summation(ni * ln(ni) - ni).
    //// So the likelihood ratio test statistic is defined as q = -2 * (ln(L_model) - ln(L_saturated)) = -2 * (summation(ni * ln(mu_i) - mu_i) - summation(ni * ln(ni) - ni)) = -2 * summation(ni * ln(mu_i/ni) - (mu_i - ni)) = 2 * summation(ni * ln(ni/mu_i) + (mu_i - ni))
    //// For large statistics n_i is close to mu_i. so q = 2 * summation(ni * ln(ni/mu_i) + (mu_i - ni)) ≈ 2 * summation((ni - mu_i)^2 / mu_i) = Chi2, which is the Pearson's Chi2.
    //// Detailed calculation is given here https://chatgpt.com/s/t_6a05df772b908191b3725caca15b2414

    //--------------------------------------------------
    // Goodness of fit (Poisson deviance)
    //--------------------------------------------------

    double deviance = 0.0;

    int nBinsUsed = 0;

    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        double x = h->GetBinCenter(i);
        double bw = h->GetBinWidth(i);

        double n = h->GetBinContent(i);

        //--------------------------------------------------
        // Expected counts from best-fit model
        //--------------------------------------------------

        double s = gaussian(x, bestMean, bestSigma);

        double b = background(x);

        double mu = (bestNs * s + bestNb * b) * bw;

        //--------------------------------------------------
        // Poisson deviance contribution
        //--------------------------------------------------

        if (n > 0 && mu > 0)
        {
            deviance += 2.0 * (mu - n + n * log(n / mu));
        }
        else if (mu > 0)
        {
            deviance += 2.0 * mu;
        }

        nBinsUsed++;
    }

    //--------------------------------------------------
    // Degrees of freedom
    //--------------------------------------------------

    int nParameters = 4;

    int ndf = nBinsUsed - nParameters;

    //--------------------------------------------------
    // Reduced deviance
    //--------------------------------------------------

    double chi2byNDF = deviance / ndf;

    //--------------------------------------------------
    // p-value from chi-square distribution
    //--------------------------------------------------

    double pvalue = TMath::Prob(deviance, ndf);

    //--------------------------------------------------
    // Print results
    //--------------------------------------------------

    cout << endl;
    cout << "========== Goodness of Fit ==========" << endl;
    cout << "Deviance            = " << deviance << endl;
    cout << "NDF                  = " << ndf << endl;
    cout << "Reduced Deviance (≈ χ^2/NDF)    = " << chi2byNDF << endl;
    cout << "p-value              = " << pvalue << endl;
    cout << "=====================================" << endl;

    //// Now lets calculate from root fit
    TF1 *fRoot = new TF1("fRoot", totalFunction, 0, 10, 4);
    fRoot->SetParameter(0, 1000);
    fRoot->SetParameter(1, 9000);
    fRoot->SetParameter(2, 5);
    fRoot->SetParameter(3, 0.3);
    TFitResultPtr fitResult = h->Fit(fRoot, "LES0");
    double chi2Root = fitResult->Chi2();
    int ndfRoot = fitResult->Ndf();
    double pvalueRoot = fitResult->Prob();
    cout << endl;
    cout << "========== Goodness of Fit from ROOT fit====================" << endl;
    cout << "Chi2                 = " << chi2Root << endl;
    cout << "NDF                  = " << ndfRoot << endl;
    cout << "Reduced Chi2 (χ^2/NDF) = " << chi2Root / ndfRoot << endl;
    cout << "p-value              = " << pvalueRoot << endl;
}
