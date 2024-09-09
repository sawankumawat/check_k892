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
    return (x[0] * par[0] * par[1] / (pow((x[0] * x[0] - par[1] * par[1]), 2) + pow((par[1] * par[2]), 2))); // par[0] is the normalization constant, par[1] is the mass, par[2] is the width
}

Double_t exponential_bkg(double *x, double *par)
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * (x[0] - 2.0 * 0.497)));
}

Double_t BW3(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + RelativisticBW(x, &par[6]));
}

Double_t expo2bW(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + exponential_bkg(x, &par[6]));
}

Double_t expo3bW(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + RelativisticBW(x, &par[6]) + exponential_bkg(x, &par[9]));
}

Double_t pol2bW3(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + RelativisticBW(x, &par[6]) + polynomial2(x, &par[9]));
}

Double_t pol3bW3(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + RelativisticBW(x, &par[6]) + polynomial3(x, &par[9]));
}

Double_t pol2bW2(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + polynomial2(x, &par[6]));
}

Double_t pol3bW2(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + polynomial3(x, &par[6]));
}

Double_t gauspol2(double *x, double *par)
{
    double gaus = par[0] * TMath::Gaus(x[0], par[1], par[2]);
    double pol2 = par[5] + par[4] * x[0] + par[3] * x[0] * x[0];
    return (gaus + pol2);
}

Double_t CrystalBall(double *x, double *par)
{
    // par[0] normalization
    // par[1] mean of gaussian
    // par[2] sigma of gaussian
    // par[3] alpha
    // par[4] n

    double t = (x[0] - par[1]) / par[2];
    double absAlpha_L = fabs(par[3]);
    double n = par[4];
    double y1 = 0;

    if (t < -absAlpha_L)
    {
        double a = exp(-0.5 * absAlpha_L * absAlpha_L) * TMath::Power(n / absAlpha_L, n);
        double b = (n / absAlpha_L) - absAlpha_L;
        y1 = par[0] * (a / TMath::Power(b - t, n));
    }
    else if (t >= -absAlpha_L)
    {
        y1 = par[0] * exp(-0.5 * t * t);
    }
    return y1;
}

Double_t DoubleCrystalBall(Double_t *x, Double_t *par)
{
    // par[0] normalization
    // par[1] mean of gaussian
    // par[2] sigma of gaussian
    // par[3] alpha left
    // par[4] n1
    // par[5] alpha right
    // par[6] n2

    double t = (x[0] - par[1]) / par[2];
    double absAlpha_L = fabs(par[3]);
    double n1 = par[4];
    double absAlpha_R = fabs(par[5]);
    double n2 = par[6];
    double y1 = 0;

    if ((t >= -absAlpha_L) && (t < absAlpha_R))
    {
        y1 = par[0] * exp(-0.5 * t * t);
    }
    else if (t < -absAlpha_L)
    {
        double a = exp(-0.5 * absAlpha_L * absAlpha_L) * TMath::Power(n1 / absAlpha_L, n1);
        double b = n1 / absAlpha_L - absAlpha_L;
        y1 = par[0] * (a / TMath::Power(b - t, n1));
    }
    else if (t >= absAlpha_R)
    {
        double a = exp(-0.5 * absAlpha_R * absAlpha_R) * TMath::Power(n2 / absAlpha_R, n2);
        double b = n2 / absAlpha_R - absAlpha_R;
        y1 = par[0] * (a / TMath::Power(b + t, n2));
    }

    return y1;
}

Double_t CrystalBallpol2(double *x, double *par)
{
    double CB = CrystalBall(x, &par[0]);
    double pol2 = par[7] + par[6] * x[0] + par[5] * x[0] * x[0];
    return (CB + pol2);
}

Double_t DoubleCrystalBallpol1(double *x, double *par)
{
    double DCB = DoubleCrystalBall(x, &par[0]);
    double pol1 = par[8] + par[7] * x[0];
    return (DCB + pol1);
}

Double_t DoubleCrystalBallpol2(double *x, double *par)
{
    double DCB = DoubleCrystalBall(x, &par[0]);
    double pol2 = par[9] + par[8] * x[0] + par[7] * x[0] * x[0];
    return (DCB + pol2);
}

Double_t DoubleCrystalBallpol3(double *x, double *par)
{
    double DCB = DoubleCrystalBall(x, &par[0]);
    double pol3 = par[10] + par[9] * x[0] + par[8] * x[0] * x[0] + par[7] * x[0] * x[0] * x[0];
    return (DCB + pol3);
}

Double_t Boltzman(double *x, double *par)
{
    double norm = par[0];
    double massks = 0.497;
    double n = par[1];
    double C = par[2];
    double y = norm * pow((x[0] - 2 * massks), n / 2) * pow(C, 3 / 2) * exp(-C * pow((x[0] - 2 * massks), n));
    return y;
}

Double_t BW3boltzman(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + RelativisticBW(x, &par[6]) + Boltzman(x, &par[9]));
}

Double_t BW2boltzman(double *x, double *par)
{
    return (RelativisticBW(x, &par[0]) + RelativisticBW(x, &par[3]) + Boltzman(x, &par[6]));
}

Double_t coherentBWsum(double *x, double *par)
{
    return par[12] * TMath::Power(5 * RelativisticBW(x, &par[0]) - 3 * RelativisticBW(x, &par[3]) + 2 * RelativisticBW(x, &par[6]), 2) + par[13] * (TMath::Power(RelativisticBW(x, &par[9]), 2));
}

Double_t coherentBWsum_expol(double *x, double *par)
{
    return coherentBWsum(x, &par[0]) + exponential_bkg(x, &par[14]);
}