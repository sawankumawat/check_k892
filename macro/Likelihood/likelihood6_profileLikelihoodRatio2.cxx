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

void likelihood6_profileLikelihoodRatio2()
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

    //------------------------------------------------
    //===============================================
    // Profile likelihood ratios
    //================================================
    //------------------------------------------------

    TGraph *gProfileNs = new TGraph();

    int p = 0;
    double prevNs = 0;
    double prevDelta = 0;

    bool foundLeft = false;
    bool foundRight = false;

    double leftCrossing = 0;
    double rightCrossing = 0;

    // scan fixed Ns values for profile likelihood ratio
    for (double Ns_fixed = 970; Ns_fixed <= 1100; Ns_fixed += 2)
    {
        //------------------------------------------
        // For this fixed Ns:
        // maximize likelihood wrt Nb
        //------------------------------------------

        double bestProfileLogL = -1e30;
        double bestProfileNb = 0;
        double bestProfileMean = 0;
        double bestProfileSigma = 0;

        for (double Nb = 8500; Nb <= 9300; Nb += 2)
        {
            for (double mean = 4.85; mean <= 5.15; mean += 0.02)
            {
                for (double sigma = 0.25; sigma <= 0.35; sigma += 0.02)
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
                        bestProfileMean = mean;
                        bestProfileSigma = sigma;
                    }
                }
            }
        }

        //------------------------------------------
        // profile likelihood ratio
        //------------------------------------------

        double q = -2.0 * (bestProfileLogL - bestLogL);
        double deltaNLL = q;

        gProfileNs->SetPoint(p, Ns_fixed, q);

        p++;

        // cout << "Ns = " << Ns_fixed << "  profiled Nb = " << bestProfileNb << "  q = " << q << endl;

        //--------------------------------------------------
        // LEFT crossing
        //--------------------------------------------------

        if (Ns_fixed < bestNs)
        {
            if (prevDelta > 1.0 && deltaNLL <= 1.0 && !foundLeft)
            {
                double fraction = (1.0 - prevDelta) / (deltaNLL - prevDelta);

                leftCrossing = prevNs + fraction * (Ns_fixed - prevNs);

                foundLeft = true;
            }
        }

        //--------------------------------------------------
        // RIGHT crossing
        //--------------------------------------------------

        if (Ns_fixed > bestNs)
        {
            if (prevDelta < 1.0 && deltaNLL >= 1.0 && !foundRight)
            {
                double fraction = (1.0 - prevDelta) / (deltaNLL - prevDelta);

                rightCrossing = prevNs + fraction * (Ns_fixed - prevNs);

                foundRight = true;
            }
        }

        //--------------------------------------------------
        // Store previous point
        //--------------------------------------------------

        prevNs = Ns_fixed;
        prevDelta = deltaNLL;
    }

    //--------------------------------------------------
    // Print uncertainties
    //--------------------------------------------------

    double errLow = bestNs - leftCrossing;
    double errHigh = rightCrossing - bestNs;
    cout << endl;
    cout << "===== PROFILE LIKELIHOOD RATIO UNCERTAINTY =====" << endl;
    cout << "Best Ns = " << bestNs << endl;
    cout << "Left crossing  = " << leftCrossing << endl;
    cout << "Right crossing = " << rightCrossing << endl;
    cout << endl;
    cout << "Ns = " << bestNs << " + " << errHigh << " - " << errLow << endl;

    TCanvas *cProfile = new TCanvas("cProfile", "Profile Likelihood", 700, 600);
    gProfileNs->SetTitle("Profile Likelihood Scan;N_{s};-2ln#lambda");
    gProfileNs->SetLineWidth(2);
    gProfileNs->Draw("AL");
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

    // similary calculate the profile likelihood ratio for mean
    double prevMean = 0;
    p = 0, prevDelta = 0;
    foundLeft = false, foundRight = false, leftCrossing = 0, rightCrossing = 0;
    TGraph *gProfileMean = new TGraph();

    // scan fixed Ns values for profile likelihood ratio
    for (double mean_fixed = 4.85; mean_fixed <= 5.15; mean_fixed += 0.01)
    {
        //------------------------------------------
        // For this fixed Ns:
        // maximize likelihood wrt Nb
        //------------------------------------------

        double bestProfileLogL = -1e30;
        double bestProfileNb = 0;
        double bestProfileMean = 0;
        double bestProfileSigma = 0;

        for (double Nb = 8700; Nb <= 9200; Nb += 4)
        {
            for (double Ns = 950; Ns <= 1070; Ns += 4)
            {
                for (double sigma = 0.25; sigma <= 0.35; sigma += 0.02)
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

                        double s = gaussian(x, mean_fixed, sigma);
                        double b = background(x);

                        //--------------------------------------
                        // mean fixed
                        //--------------------------------------

                        double mu = (Ns * s + Nb * b) * bw;

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
                        bestProfileMean = mean_fixed;
                        bestProfileSigma = sigma;
                    }
                }
            }
        }

        //------------------------------------------
        // profile likelihood ratio
        //------------------------------------------

        double q = -2.0 * (bestProfileLogL - bestLogL);
        double deltaNLL = q;

        gProfileMean->SetPoint(p, mean_fixed, q);

        p++;

        // cout << "Ns = " << Ns << "  profiled Nb = " << bestProfileNb << "  q = " << q << endl;

        //--------------------------------------------------
        // LEFT crossing
        //--------------------------------------------------

        if (mean_fixed < bestMean)
        {
            if (prevDelta > 1.0 && deltaNLL <= 1.0 && !foundLeft)
            {
                double fraction = (1.0 - prevDelta) / (deltaNLL - prevDelta);

                leftCrossing = prevMean + fraction * (mean_fixed - prevMean);

                foundLeft = true;
            }
        }

        //--------------------------------------------------
        // RIGHT crossing
        //--------------------------------------------------

        if (mean_fixed > bestMean)
        {
            if (prevDelta < 1.0 && deltaNLL >= 1.0 && !foundRight)
            {
                double fraction = (1.0 - prevDelta) / (deltaNLL - prevDelta);

                rightCrossing = prevMean + fraction * (mean_fixed - prevMean);

                foundRight = true;
            }
        }

        //--------------------------------------------------
        // Store previous point
        //--------------------------------------------------

        prevMean = mean_fixed;
        prevDelta = deltaNLL;
    }

    double errLowMean = bestMean - leftCrossing;
    double errHighMean = rightCrossing - bestMean;
    cout << endl;
    cout << "===== PROFILE LIKELIHOOD RATIO UNCERTAINTY FOR MEAN =====" << endl;
    cout << "Best mean = " << bestMean << endl;
    cout << "Left crossing  = " << leftCrossing << endl;
    cout << "Right crossing = " << rightCrossing << endl;
    cout << endl;
    cout << "Mean = " << bestMean << " + " << errHighMean << " - " << errLowMean << endl;

    TCanvas *cProfileMean = new TCanvas("cProfileMean", "Profile Likelihood for Mean", 720, 720);
    gProfileMean->SetTitle("Profile Likelihood Scan for Mean;Mean;-2ln#lambda");
    gProfileMean->SetLineWidth(2);
    gProfileMean->Draw("AL");
    TLine *line1SigmaMean = new TLine(4.85, 1, 5.15, 1);
    line1SigmaMean->SetLineColor(kRed);
    line1SigmaMean->SetLineStyle(2);
    line1SigmaMean->Draw("same");
    TLine *lineBestFitMean = new TLine(bestMean, 0, bestMean, 5);
    lineBestFitMean->SetLineColor(kBlue);
    lineBestFitMean->SetLineStyle(2);
    lineBestFitMean->Draw("same");
    TLegend *legMean = new TLegend(0.6, 0.7, 0.9, 0.9);
    legMean->SetFillStyle(0);
    legMean->SetBorderSize(0);
    legMean->SetTextFont(42);
    legMean->SetTextSize(0.035);
    legMean->AddEntry(line1SigmaMean, "1#sigma", "l");
    legMean->AddEntry(lineBestFitMean, "Best mean value", "l");
    legMean->Draw();

    // similary calculate the profile likelihood ratio for sigma
    double prevSigma = 0;
    p = 0, prevDelta = 0;
    foundLeft = false, foundRight = false, leftCrossing = 0, rightCrossing = 0;
    TGraph *gProfileSigma = new TGraph();

    // scan fixed Ns values for profile likelihood ratio
    for (double sigma_fixed = 0.26; sigma_fixed <= 0.34; sigma_fixed += 0.005)
    {
        //------------------------------------------
        // For this fixed sigma:
        // maximize likelihood wrt Nb
        //------------------------------------------
        double bestProfileLogL = -1e30;
        double bestProfileNb = 0;
        double bestProfileMean = 0;
        double bestProfileSigma = 0;

        for (double Nb = 8700; Nb <= 9200; Nb += 4)
        {
            for (double Ns = 950; Ns <= 1070; Ns += 4)
            {
                for (double mean = 4.85; mean <= 5.15; mean += 0.02)
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

                        double s = gaussian(x, mean, sigma_fixed);
                        double b = background(x);

                        //--------------------------------------
                        // sigma fixed
                        //--------------------------------------

                        double mu = (Ns * s + Nb * b) * bw;

                        if (mu > 0)
                        {
                            logL += n * log(mu) - mu;
                        }
                    }
                    //--------------------------------------
                    // keep best Nb for this sigma
                    //--------------------------------------

                    if (logL > bestProfileLogL)
                    {
                        bestProfileLogL = logL;
                        bestProfileNb = Nb;
                        bestProfileMean = mean;
                        bestProfileSigma = sigma_fixed;
                    }
                }
            }
        }
        //------------------------------------------
        // profile likelihood ratio
        //------------------------------------------
        double q = -2.0 * (bestProfileLogL - bestLogL);
        double deltaNLL = q;
        gProfileSigma->SetPoint(p, sigma_fixed, q);
        p++;

        //--------------------------------------------------
        // LEFT crossing
        //--------------------------------------------------
        if (sigma_fixed < bestSigma)
        {
            if (prevDelta > 1.0 && deltaNLL <= 1.0 && !foundLeft)
            {
                double fraction = (1.0 - prevDelta) / (deltaNLL - prevDelta);

                leftCrossing = prevSigma + fraction * (sigma_fixed - prevSigma);

                foundLeft = true;
            }
        }

        //--------------------------------------------------
        // Right crossing
        //--------------------------------------------------
        if (sigma_fixed > bestSigma)
        {
            if (prevDelta < 1.0 && deltaNLL >= 1.0 && !foundRight)
            {
                double fraction = (1.0 - prevDelta) / (deltaNLL - prevDelta);

                rightCrossing = prevSigma + fraction * (sigma_fixed - prevSigma);

                foundRight = true;
            }
        }
        //-------------------------------------------------- Store previous point
        //--------------------------------------------------
        prevSigma = sigma_fixed;
        prevDelta = deltaNLL;
    }

    double errLowSigma = bestSigma - leftCrossing;
    double errHighSigma = rightCrossing - bestSigma;
    cout << endl;
    cout << "===== PROFILE LIKELIHOOD RATIO UNCERTAINTY FOR SIGMA =====" << endl;
    cout << "Best sigma = " << bestSigma << endl;
    cout << "Left crossing  = " << leftCrossing << endl;
    cout << "Right crossing = " << rightCrossing << endl;
    cout << endl;
    cout << "Sigma = " << bestSigma << " + " << errHighSigma << " - " << errLowSigma << endl;

    TCanvas *cProfileSigma = new TCanvas("cProfileSigma", "Profile Likelihood for Sigma", 720, 720);
    gProfileSigma->SetTitle("Profile Likelihood Scan for Sigma;Sigma;-2ln#lambda");
    gProfileSigma->SetLineWidth(2);
    gProfileSigma->Draw("AL");
    TLine *line1SigmaSigma = new TLine(0.25, 1, 0.35, 1);
    line1SigmaSigma->SetLineColor(kRed);
    line1SigmaSigma->SetLineStyle(2);
    line1SigmaSigma->Draw("same");
    TLine *lineBestFitSigma = new TLine(bestSigma, 0, bestSigma, 5);
    lineBestFitSigma->SetLineColor(kBlue);
    lineBestFitSigma->SetLineStyle(2);
    lineBestFitSigma->Draw("same");
    TLegend *legSigma = new TLegend(0.6, 0.7, 0.9, 0.9);
    legSigma->SetFillStyle(0);
    legSigma->SetBorderSize(0);
    legSigma->SetTextFont(42);
    legSigma->SetTextSize(0.035);
    legSigma->AddEntry(line1SigmaSigma, "1#sigma", "l");
    legSigma->AddEntry(lineBestFitSigma, "Best sigma value", "l");
    legSigma->Draw();

    // Complete summary of the results obtained manually
    cout << endl;
    cout << "========== Summary of results from manual method==============" << endl;
    cout << "Best signal yield     = " << bestNs << " + " << errHigh << " - " << errLow << endl;
    cout << "Best background yield = " << bestNb << endl; // profile likelihood ratio not performed for it
    cout << "Best mean        = " << bestMean << " + " << errHighMean << " - " << errLowMean << endl;
    cout << "Best sigma       = " << bestSigma << " + " << errHighSigma << " - " << errLowSigma << endl;

    // Now claculate the same thing using root method.
    TF1 *fRoot = new TF1("fRoot", totalFunction, 0, 10, 4);
    fRoot->SetParameter(0, 1000);
    fRoot->SetParameter(1, 9000);
    fRoot->SetParameter(2, 5.00);
    fRoot->SetParameter(3, 0.30);
    TFitResultPtr fitResult = h->Fit(fRoot, "LES0");
    cout << endl;
    cout << "========== Summary of results from ROOT method==============" << endl;
    cout << "Best signal yield     = " << fRoot->GetParameter(0) << " + " << fitResult->LowerError(0) << " - " << fitResult->UpperError(0) << endl;
    cout << "Best background yield = " << fRoot->GetParameter(1) << " + " << fitResult->LowerError(1) << " - " << fitResult->UpperError(1) << endl;
    cout << "Best mean        = " << fRoot->GetParameter(2) << " + " << fitResult->LowerError(2) << " - " << fitResult->UpperError(2) << endl;
    cout << "Best sigma       = " << fRoot->GetParameter(3) << " + " << fitResult->LowerError(3) << " - " << fitResult->UpperError(3) << endl;

    TCanvas *cFitRoot = new TCanvas("cFitRoot", "Fit with ROOT", 720, 720);
    h->SetLineColor(kBlack);
    h->Draw();
    fRoot->SetLineStyle(2);
    fRoot->Draw("same");
}

//// Results
// Processing likelihood6_profileLikelihoodRatio2.cxx...
// Best fit results
// Best signal yield     = 1020
// Best background yield = 8980
// Best mean        = 4.99
// Best sigma       = 0.29
// Maximum logL     = 36809.2

// ===== PROFILE LIKELIHOOD RATIO UNCERTAINTY =====
// Best Ns = 1020
// Left crossing  = 973.197
// Right crossing = 1066.52

// Ns = 1020 + 46.5165 - 46.8025

// ===== PROFILE LIKELIHOOD RATIO UNCERTAINTY FOR MEAN =====
// Best mean = 4.99
// Left crossing  = 4.9743
// Right crossing = 5.00511

// Mean = 4.99 + 0.0151103 - 0.0156955

// ===== PROFILE LIKELIHOOD RATIO UNCERTAINTY FOR SIGMA =====
// Best sigma = 0.29
// Left crossing  = 0.276028
// Right crossing = 0.301192

// Sigma = 0.29 + 0.0111919 - 0.013972


/////================comparing uncertainties between manual and root method=================
//========== Best fit results==============
// Best signal yield     = 1030 + 36.2068 - 48.4212
// Best background yield = 8968
// Best mean        = 5.02 + 0.0166951 - 0.0128802
// Best sigma       = 0.3 + 0.0112386 - 0.0169526
// ****************************************
// Minimizer is Minuit2 / Migrad
// MinFCN                    =      45.6811
// Chi2                      =      91.3621
// NDf                       =           96
// Edm                       =  1.63017e-08
// NCalls                    =           81
// p0                        =       1028.6   +/-   52.1845       -51.6593     +52.7214      (Minos) 
// p1                        =      8971.39   +/-   103.275       -102.935     +103.621      (Minos) 
// p2                        =      5.02209   +/-   0.0165449     -0.0165433   +0.0165688    (Minos) 
// p3                        =     0.296443   +/-   0.0154775     -0.0150476   +0.0159441    (Minos) 

// Results from ROOT fit with option LES0 (profile likelihood ratio method in ROOT)
// Best signal yield     = 1028.6 + -51.6593 - 52.7214
// Best background yield = 8971.39 + -102.935 - 103.621
// Best mean        = 5.02209 + -0.0165433 - 0.0165688
// Best sigma       = 0.296443 + -0.0150476 - 0.0159441

