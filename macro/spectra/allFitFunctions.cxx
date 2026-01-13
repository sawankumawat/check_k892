
// Levy-Tsallis function
Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

// Exponential in dN/dpT times pT
Double_t FuncExpdNdptTimesPt(Double_t *x, Double_t *par)
{
    // par[0] = norm
    // par[1] = T
    // x[0]   = pT

    Double_t pT = x[0];
    Double_t norm = par[0];
    Double_t T = par[1];

    // Double_t p = norm * pT * TMath::Exp(-pT / T); // not normalized
    Double_t p = norm * pT * TMath::Exp(-pT / T) * (1.0 / (T * T)); // normalized
    return p;
}

// Exponential in dN/dmT times pT (fix the mass while using the function)
Double_t FuncMTExpdNdptTimesPt(Double_t *x, Double_t *par)
{

    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    // set mass explicitly (example: K0s)
    // const Double_t mass = 0.497611;

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    // return norm * pT * TMath::Exp(-mT / T); // not normalized
    return norm * pT * TMath::Exp(-mT / T) * TMath::Exp(mass / T) * (1.0 / (T * T + T * mass)); // normalized
}
// Boltzmann in dN/dpT times pT (fix the mass while using the function)
Double_t FuncBoltzmanndNdptTimesPt(Double_t *x, Double_t *par)
{

    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    // set mass explicitly
    // const Double_t mass = 0.497611; // example: K0s

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    // return norm * pT * mT * TMath::Exp(-mT / T); // Not normalized
    return norm * pT * mT * TMath::Exp(-mT / T) * TMath::Exp(mass / T) * (1.0 / (2 * T * T * T + 2 * T * T * mass + T * mass * mass)); // normalized
}

// Bose-Einstein in dN/dpT times pT (fix the mass while using the function)
Double_t FuncBoseEinsteindNdptTimesPt(Double_t *x, Double_t *par)
{

    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    // set mass explicitly
    // const Double_t mass = 0.497611; // example: K0s

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    return norm * pT / (TMath::Exp(mT / T) - 1.0); // not normalized
    // return norm * pT / (TMath::Exp(mT / T) - 1.0) * (TMath::Exp(mass / T) - 1.0); // normalized (may not be correct)
}
// Fermi-Dirac in dN/dpT times pT (fix the mass while using the function)
Double_t FuncFermiDiracdNdptTimesPt(Double_t *x, Double_t *par)
{

    const Double_t pT = x[0];
    const Double_t norm = par[0];
    const Double_t T = par[1];
    const Double_t mass = par[2];

    // set mass explicitly
    // const Double_t mass = 0.497611; // example: K0s

    const Double_t mT = TMath::Sqrt(pT * pT + mass * mass);

    return norm * pT / (TMath::Exp(mT / T) + 1.0); // not normalized
}
