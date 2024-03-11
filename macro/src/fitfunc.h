//******Different fitting functions*****************************************************************************************

#include <cmath>

Double_t BreitWignerpoly3(Double_t *x, Double_t *par)
{
    double BW = (0.5 * par[2] * par[1] / TMath::Pi()) / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]); // parameter is BW mass, 1 is width, 2 is the yield. x[0] is the invariant mass
    double poly3 = par[6] + par[5] * x[0] + par[4] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
    return (BW + poly3);
}

Double_t voigtpol2(Double_t *x, Double_t *par)
{
    double vgt = par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3], 4);
    double poly2 = par[6] + par[5] * x[0] + par[4] * x[0] * x[0];
    return (vgt + poly2);
}

Double_t voigt(Double_t *x, Double_t *par)
{
    return (par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]));
}

Double_t voigt_phi(Double_t *x, Double_t *par)
{
    double vgt = (par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]));
    double srv = (par[6] + par[5] * x[0] + par[4] * x[0] * x[0]);
    return (vgt + srv);
}

Double_t phi_bkg(Double_t *x, Double_t *par)
{
    double srv = (par[6] + par[5] * x[0] + par[4] * x[0] * x[0]);
    return (srv);
}

Double_t BreitWignerpoly2(Double_t *x, Double_t *par)
{
    double BW = (0.5 * par[2] * par[1] / TMath::Pi()) / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]); // parameter 0 is invariant mass, 1 is width, 2 is the yield
    double poly2 = par[5] + par[4] * x[0] + par[3] * x[0] * x[0];
    return (BW + poly2);
}

Double_t BW(Double_t *x, Double_t *par)
{
    return (0.5 * par[2] * par[1] / TMath::Pi() / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]));
}

Double_t BreitWignerpoly1(Double_t *x, Double_t *par)
{
    double BW = (0.5 * par[2] * par[1] / 3.14159) / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]);
    double poly1 = par[4] + par[3] * x[0];
    // double poly3 = par[0] + par[1]*x[0] + par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0] ;
    return (BW + poly1);
}

Double_t polynomial2(Double_t *x, Double_t *par)
{
    double poly2 = par[2] + par[1] * x[0] + par[0] * x[0] * x[0];
    return (poly2);
}

Double_t polynomial1(Double_t *x, Double_t *par)
{
    return (par[0] + par[1] * x[0]);
}

Double_t polynomial3(Double_t *x, Double_t *par)
{
    double poly3 = par[3] + par[2] * x[0] + par[1] * x[0] * x[0] + par[0] * x[0] * x[0] * x[0];
    return (poly3);
}

Double_t BWExpo(Double_t *x, Double_t *par)
{
    double BW = (0.5 * par[2] * par[1] / TMath::Pi()) / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]);
    double expo = (pow((x[0] - 0.63718), par[3])) * exp(-par[4] - x[0] * par[5] - x[0] * x[0] * par[6]);
    return (BW + expo);
}

Double_t Expo(Double_t *x, Double_t *par)
{
    // return ((x[0]- 0.63718)**par[3])*exp(par[0] + x[0]*par[1] + x[0]*x[0]*par[2]);
    double expo = (pow((x[0] - 0.63267), par[0])) * exp(-par[3] - x[0] * par[2] - x[0] * x[0] * par[1]);
    return (expo);
}

//*****************************glueball fit functions****************

Double_t RelativisticBW(double *x, double *par)
{
    return (par[0] / (pow((x[0] * x[0] - par[1] * par[1]), 2) + pow((par[1] * par[2]), 2))); // par[0] is the normalization constant, par[1] is the mass, par[2] is the width
}

double exponential_bkg(double *x, double *par)
{
    return (par[0] * pow((x[0] - 2 * 0.497), par[1]) * TMath::Exp(-par[2] * (x[0] - 2 * 0.497)));
}

double gluefit2bW(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + exponential_bkg(x, &par[6]));
}

double gluefit3bW(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + RelativisticBW(x, &par[6]) + exponential_bkg(x, &par[9]));
}