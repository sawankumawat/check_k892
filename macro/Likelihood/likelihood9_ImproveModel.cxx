#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"

//// In this code we will calculate the p-value by generating the toy events using null model and checking which toy events have a test statistic (q) equal to or more than the observed value. This is a more direct way to compute the p-value without relying on the asymptotic approximation of the test statistic distribution.
//// Results are pasted in the end of the code.

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

double background_pol1(double x)
{
    return (0.1 + 0.02 * x) / 2.0; // Factor 2.0 is for normalization over the range 0-10
}

double background_pol2(double x, double p0, double p1, double p2)
{
    double poly = p0 + p1 * x + p2 * x * x;
    double norm = (p0 * 10.0 + p1 * pow(10.0, 2) / 2.0 + p2 * pow(10.0, 3) / 3.0); // normalized for range 0-10
    return poly / norm;
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
    double b = background_pol1(xx);

    // expected counts per bin
    double bw = 0.1;

    return (Ns * s + Nb * b) * bw;
}

double totalFunction_pol2(double *x, double *par)
{
    double xx = x[0];

    double Ns = par[0];
    double Nb = par[1];

    double mean = 5.0;
    double sigma = 0.3;

    double s = gaussian(xx, mean, sigma);
    double b = background_pol2(xx, par[2], par[3], par[4]);

    // expected counts per bin
    double bw = 0.1;

    return (Ns * s + Nb * b) * bw;
}

//--------------------------------------------------
// Main function
//--------------------------------------------------

void likelihood9_ImproveModel()
{
    //--------------------------------------------------
    // Create toy histogram
    //--------------------------------------------------

    TH1D *h = new TH1D("h", "Toy data", 100, 0, 10); // Note that background is normalized for this range. If you change the range here, make sure to update the normalization in the background function.

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

    // TF1 *fbkg = new TF1("fbkg", "0.1 + 0.02*x", 0, 10); // pol1 background
    TF1 *fbkg = new TF1("fbkg", "0.1 + 0.02*x + 0.01*x*x", 0, 10); // pol2 background

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

    for (double Ns = 920; Ns <= 1050; Ns += 5)
    {
        for (double Nb = 8850; Nb <= 10000; Nb += 5)
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

                double b = background_pol1(x);

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
    // POL2 BACKGROUND MODEL to check improvement in likelihood
    //--------------------------------------------------

    double bestLogL_pol2 = -1e30;

    double bestNs_pol2 = 0;
    double bestNb_pol2 = 0;

    double bestP2 = 0;

    //--------------------------------------------------
    // scan additional cubic parameter
    //--------------------------------------------------

    for (double Ns = 920; Ns <= 1050; Ns += 5)
    {
        for (double Nb = 8850; Nb <= 10000; Nb += 5)
        {
            //--------------------------------------------------
            // NEW parameter
            //--------------------------------------------------

            for (double p2 = 0.008; p2 <= 0.012; p2 += 0.0005)
            {
                double logL = 0;

                for (int i = 1; i <= h->GetNbinsX(); i++)
                {
                    double x = h->GetBinCenter(i);

                    double bw = h->GetBinWidth(i);

                    double n = h->GetBinContent(i);

                    double s = gaussian(x, mean, sigma);

                    //--------------------------------------------------
                    // POL2 background
                    //--------------------------------------------------

                    double b = background_pol2(x, 0.1, 0.02, p2);

                    double mu = (Ns * s + Nb * b) * bw;

                    if (mu > 0)
                    {
                        logL += n * log(mu) - mu;
                    }
                }

                //--------------------------------------------------
                // check improvement
                //--------------------------------------------------

                if (logL > bestLogL_pol2)
                {
                    bestLogL_pol2 = logL;

                    bestNs_pol2 = Ns;

                    bestNb_pol2 = Nb;

                    bestP2 = p2;
                }
            }
        }
    }

    //--------------------------------------------------
    // PROFILE LIKELIHOOD MODEL TEST
    //--------------------------------------------------

    double q_model = 2.0 * (bestLogL_pol2 - bestLogL);

    //--------------------------------------------------
    // significance of extra parameter
    //--------------------------------------------------

    double Z_model = sqrt(q_model);

    //--------------------------------------------------
    // p-value
    //--------------------------------------------------

    double p_model = TMath::Prob(q_model, 1);

    cout << endl;

    cout << "===== POL1 vs POL2 =====" << endl;
    cout << "Best logL (pol1) = " << bestLogL << endl;
    cout << "Best logL (pol2) = " << bestLogL_pol2 << endl;
    cout << "Best square parameter = " << bestP2 << endl;
    cout << "q_model = " << q_model << endl;
    cout << "significance = " << Z_model << " sigma" << endl;
    cout << "p-value = " << p_model << endl;

    TF1 *f_pol2 = new TF1("f_pol2", totalFunction_pol2, 0, 10, 5);
    f_pol2->SetParameter(0, bestNs_pol2);
    f_pol2->SetParameter(1, bestNb_pol2);
    f_pol2->SetParameter(2, 0.1);
    f_pol2->SetParameter(3, 0.02);
    f_pol2->SetParameter(4, bestP2);
    f_pol2->SetLineStyle(3);
    f_pol2->SetLineColor(kBlue);
    f_pol2->Draw("same");
}

//===============Results===============
// 1. For data generated with pol1 background

// Best fit results
// Signal yield     = 1030
// Background yield = 8970
// Maximum logL     = 36806.7
// Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1

// ===== POL1 vs POL2 =====
// Best logL (pol1) = 36806.7
// Best logL (pol2) = 36806.9
// Best square parameter = 6e-05
// q_model = 0.410779
// significance = 0.64092 sigma
// p-value = 0.521574
// == The results are as expected. The improvement in likelihood is very small and the p-value is large, indicating that the extra parameter in the pol2 model does not significantly improve the fit compared to the pol1 model, which is consistent with the fact that the data was generated with a pol1 background. This demonstrates how we can use the likelihood ratio test to compare nested models and assess whether the more complex model provides a significantly better fit to the data.

// 2. For data generated with pol2 background
// Best fit results
// Signal yield     = 950
// Background yield = 9112
// Maximum logL     = 37617
// Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1

// ===== POL1 vs POL2 =====
// Best logL (pol1) = 37617
// Best logL (pol2) = 38309.4
// Best square parameter = 0.0105
// q_model = 1384.75
// significance = 37.2122 sigma
// p-value = 4.32945e-303
// == Again, the results are as expected. The improvement in likelihood is very large and the p-value is extremely small, indicating that the extra parameter in the pol2 model significantly improves the fit compared to the pol1 model, which is consistent with the fact that the data was generated with a pol2 background.
