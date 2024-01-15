//******Different fitting functions*****************************************************************************************

Double_t BreitWignerpoly3(Double_t *x, Double_t *par)
{
    double BW = (0.5 * par[2] * par[1] / TMath::Pi()) / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]); // parameter 0 is invariant mass, 1 is width, 2 is the yield
    double poly3 = par[6] + par[5] * x[0] + par[4] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
    return (BW + poly3);
}

Double_t voigtpol2(Double_t *x, Double_t *par)
{
    double vgt = par[0]*TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    double poly2 = par[6] + par[5]*x[0] + par[4]*x[0]*x[0];
    return (vgt+poly2);
}

Double_t voigt(Double_t *x, Double_t *par)
{
    return(par[0]*TMath::Voigt(x[0]-par[1],par[2],par[3]));

}

Double_t voigt_phi(Double_t *x, Double_t *par)
{
    double vgt = (par[0]*TMath::Voigt(x[0]-par[1],par[2],par[3]));
    double srv = (par[6] + par[5]*x[0] + par[4] * x[0]*x[0]);
    return (vgt+srv);
}

Double_t phi_bkg(Double_t *x, Double_t *par)
{
    double srv = (par[6] + par[5]*x[0] + par[4] * x[0]*x[0]);
    return (srv);
}

Double_t BreitWignerpoly2(Double_t *x, Double_t *par)
{
    double BW = (0.5 * par[2] * par[1] / TMath::Pi()) / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]); // parameter 0 is invariant mass, 1 is width, 2 is the yield
    double poly2 = par[5] + par[4]*x[0] + par[3]*x[0]*x[0];
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
    double expo = (pow((x[0] - 0.63718), par[3])) * exp(par[4] + x[0] * par[5] + x[0] * x[0] * par[6]);
    return (BW+expo);
}

Double_t Expo(Double_t *x, Double_t *par)
{
    // return ((x[0]- 0.63718)**par[3])*exp(par[0] + x[0]*par[1] + x[0]*x[0]*par[2]);
    return (pow((x[0] - 0.63718), par[0])) * exp(par[1] + x[0] * par[2] + x[0] * x[0] * par[3]);
    // return TMath::Power((x[0]-(0.13957+0.49367)),par[3])*exp(par[2] + par[1]*x[0]+par[0]*x[0]*x[0]);
}