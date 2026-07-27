#include "TH2D.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TMath.h"
#include <complex>

using namespace std;

//------------------------------------------------------------
// Real spherical harmonics
//------------------------------------------------------------

// Y00
double Y00(double c)
{
    return sqrt(1. / (4. * TMath::Pi()));
}

// Real Y22 = sqrt(15/32pi) sin²(theta) cos(2phi)
double Y22(double c, double phi)
{
    double s2 = 1.0 - c * c;

    return sqrt(15. / (32. * TMath::Pi())) * s2 * cos(2. * phi);
}

//------------------------------------------------------------
// Fit function
//------------------------------------------------------------
Double_t Intensity(Double_t *x, Double_t *par)
{
    double c = x[0];
    double phi = x[1];

    double r0 = par[0];
    double r2 = par[1];
    double delta = par[2];

    double A0 = r0 * Y00(c);

    double A2 = r2 * Y22(c, phi);

    return A0 * A0 + A2 * A2 + 2. * A0 * A2 * cos(delta);
}

//------------------------------------------------------------
void Angular2Dfit()
{
    TRandom3 rnd(0);

    TH2D *h = new TH2D("h",
                       ";cos#theta;#phi",
                       50, -1, 1,
                       72, -TMath::Pi(), TMath::Pi());

    //--------------------------------------------------------
    // TRUE parameters
    //--------------------------------------------------------

    double r0True = 0.6;
    double r2True = 0.8;
    double phaseTrue = 0.7;

    //--------------------------------------------------------
    // Maximum intensity
    //--------------------------------------------------------

    double maxI = 0;

    for (int i = 0; i < 300; i++)
    {
        double c = -1. + 2. * i / 299.;

        for (int j = 0; j < 300; j++)
        {
            double phi = -TMath::Pi() + 2. * TMath::Pi() * j / 299.;

            double xx[2] = {c, phi};
            double pp[3] = {r0True, r2True, phaseTrue};

            double val = Intensity(xx, pp);

            if (val > maxI)
                maxI = val;
        }
    }

    //--------------------------------------------------------
    // Generate events
    //--------------------------------------------------------

    int Nev = 2000000;
    int n = 0;

    while (n < Nev)
    {
        double c = rnd.Uniform(-1, 1);
        double phi = rnd.Uniform(-TMath::Pi(), TMath::Pi());

        double xx[2] = {c, phi};
        double pp[3] = {r0True, r2True, phaseTrue};

        double I = Intensity(xx, pp);

        if (rnd.Uniform(0, maxI) < I)
        {
            h->Fill(c, phi);
            n++;
        }
    }

    //--------------------------------------------------------
    // Fit
    //--------------------------------------------------------

    TF2 *fit = new TF2("fit",
                       Intensity,
                       -1, 1,
                       -TMath::Pi(), TMath::Pi(),
                       3);

    fit->SetParNames("r0", "r2", "phase");

    fit->SetParameters(0.5, 0.9, 0.7);

    fit->SetParLimits(0, 0, 1000);
    fit->SetParLimits(1, 0, 1000);
    fit->SetParLimits(2, -TMath::Pi(), TMath::Pi());

    h->Fit(fit, "REBMS");

    printf("\nTrue values\n");
    printf("r0 = %.2f\n", r0True);
    printf("r2 = %.2f\n", r2True);
    printf("Ratio r2/r0 = %.2f\n", r2True / r0True);
    printf("phase = %.2f\n\n", phaseTrue);

    printf("Fit values\n");
    printf("r0 = %.2f\n", fit->GetParameter(0));
    printf("r2 = %.2f\n", fit->GetParameter(1));
    printf("Ratio r2/r0 = %.2f\n", fit->GetParameter(1) / fit->GetParameter(0));
    printf("phase = %.2f\n", fit->GetParameter(2));

    //--------------------------------------------------------
    // Draw
    //--------------------------------------------------------

    TCanvas *c1 = new TCanvas("c1", "", 800, 700);

    h->Draw("COLZ");

    TCanvas *c2 = new TCanvas("c2", "", 800, 700);

    fit->Draw("SURF1");
}