#include <iostream>
#include "../src/style.h"

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;
Double_t BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{
    //  x[0] -> r (radius)
    //  p[0] -> mT (transverse mass)
    //  p[1] -> pT (transverse momentum)
    //  p[2] -> beta_max (surface velocity)
    //  p[3] -> T (freezout temperature)
    //  p[4] -> n (velocity profile)

    Double_t r = x[0];
    Double_t mt = p[0];
    Double_t pt = p[1];
    Double_t beta_max = p[2];
    Double_t temp_1 = 1. / p[3];
    Double_t n = p[4];

    Double_t beta = beta_max * TMath::Power(r, n);
    if (beta > 0.9999999999999999)
        beta = 0.9999999999999999;
    Double_t rho = TMath::ATanH(beta);
    Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
    if (argI0 > 700.)
        argI0 = 700.;
    Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
    return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
}

Double_t BGBlastWave_Func(const Double_t *x, const Double_t *p)
{

    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t beta_max = p[1];
    Double_t temp = p[2];
    Double_t n = p[3];
    Double_t norm = p[4];

    if (!fBGBlastWave_Integrand)
        fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0.0, 5.0, 5);
    fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
    // Double_t integral = fBGBlastWave_Integrand->Integral(0.0, 5.0, (Double_t *)0, 1.e4);
    Double_t integral = fBGBlastWave_Integrand->Integral(0.0, 1.0, 1.e4);
    return norm * pt * integral;
}

double BlastWave(double pt, double mass, double beta, double T, double n, double norm)
{
    double mt = sqrt(pt * pt + mass * mass);

    if (!fBGBlastWave_Integrand)
        fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0.0, 1.0, 5);

    fBGBlastWave_Integrand->SetParameters(mt, pt, beta, T, n);

    double integral = fBGBlastWave_Integrand->Integral(0., 1.);

    return norm * pt * integral;
}

double Chi2Hist(TH1D *h, double mass, double beta, double T, double n, double norm)
{
    double chi2 = 0.0;

    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        double pt = h->GetBinCenter(i);
        double y = h->GetBinContent(i);
        double err = h->GetBinError(i);

        if (err <= 0)
            continue;

        double model = BlastWave(pt, mass, beta, T, n, norm);

        chi2 += pow((y - model) / err, 2);
    }

    return chi2;
}

void GlobalChi2(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    double beta = par[0];
    double T = par[1];
    double n = par[2];

    double normPi = par[3];
    double normK = par[4];
    double normP = par[5];

    fval = 0;
    fval += Chi2Hist(gPi, 0.13957, beta, T, n, normPi);
    fval += Chi2Hist(gK, 0.49367, beta, T, n, normK);
    fval += Chi2Hist(gP, 0.93827, beta, T, n, normP);
}

TH1D *GetHisto(TFile *f, const string &name)
{
    TH1D *histo = (TH1D *)f->Get(name.c_str());

    if (!histo || histo == nullptr)
    {
        cout << "Error: histo " << name
             << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetHistoQA(histo);
    histo->SetTitle(0);
    return histo;
}

TFile *OpenFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (f->IsZombie())
    {
        cout << "Error: File not found: " << path << endl;
        return nullptr;
    }
    return f;
}

void simultaneousBlastWaveFit()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    string KstarPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/";
    TFile *fPion = OpenFile("PiKp_Run3_Results/Sawan/Pi_results.root");
    TFile *fProton = OpenFile("PiKp_Run3_Results/Sawan/Pr_results.root");
    TFile *fKaon = OpenFile("PiKp_Run3_Results/Sawan/Ka_results.root");

    if (!fPion || !fProton || !fKaon)
    {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    TH1D *hSpectraKstar[10], *hSpectraPion[10], *hSpectraProton[10], *hSpectraKaon[10];
    int multClases[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};

    for (int i = 0; i < 10; i++)
    {
        TFile *fKstarMult = OpenFile(KstarPath + Form("corrected_spectra_%d_%d.root", multClases[i], multClases[i + 1]));
        hSpectraKstar[i] = GetHisto(fKstarMult, Form("mult_%d-%d/corrected_spectra_Integral_final", multClases[i], multClases[i + 1]));
        hSpectraPion[i] = GetHisto(fPion, Form("hCorrectedSpectra_%d_%d", multClases[i], multClases[i + 1]));
        hSpectraProton[i] = GetHisto(fProton, Form("hCorrectedSpectra_%d_%d", multClases[i], multClases[i + 1]));
        hSpectraKaon[i] = GetHisto(fKaon, Form("hCorrectedSpectra_%d_%d", multClases[i], multClases[i + 1]));
    }

    TCanvas *cPiKp = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cPiKp, 0.12, 0.03, 0.03, 0.14);
    gPad->SetLogy();
    hSpectraPion[0]->SetMarkerColor(kRed);
    hSpectraPion[0]->SetLineColor(kRed);
    hSpectraPion[0]->SetMarkerStyle(20);
    hSpectraPion[0]->GetYaxis()->SetTitle("dN/dy");
    hSpectraPion[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    hSpectraPion[0]->GetYaxis()->SetTitleOffset(1.0);
    hSpectraPion[0]->GetXaxis()->SetTitleOffset(1.2);
    hSpectraPion[0]->SetMinimum(1e-2);
    hSpectraPion[0]->Draw("PE");
    hSpectraProton[0]->SetMarkerColor(kBlue);
    hSpectraProton[0]->SetLineColor(kBlue);
    hSpectraProton[0]->SetMarkerStyle(21);
    hSpectraProton[0]->Draw("PE same");
    hSpectraKaon[0]->SetMarkerColor(kGreen + 2);
    hSpectraKaon[0]->SetLineColor(kGreen + 2);
    hSpectraKaon[0]->SetMarkerStyle(22);
    hSpectraKaon[0]->Draw("PE same");
    hSpectraKstar[0]->SetMarkerColor(kMagenta);
    hSpectraKstar[0]->SetLineColor(kMagenta);
    hSpectraKstar[0]->SetMarkerStyle(23);
    hSpectraKstar[0]->Draw("PE same");

    TVirtualFitter *minuit = TVirtualFitter::Fitter(0, 6);

    minuit->SetFCN(GlobalChi2);

    double arglist[10];
    int ierflg = 0;

    minuit->SetParameter(0, "beta_s", 0.8, 0.01, 0.1, 0.99);
    minuit->SetParameter(1, "T", 0.10, 0.001, 0.05, 0.20);
    minuit->SetParameter(2, "n", 1.0, 0.01, 0.1, 3.0);
    minuit->SetParameter(3, "NormPi", 100, 1, 0, 1e10);
    minuit->SetParameter(4, "NormK", 10, 1, 0, 1e10);
    minuit->SetParameter(5, "NormP", 5, 1, 0, 1e10);

    arglist[0] = 10000;
    arglist[1] = 1;

    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    double beta, betaErr;
    double T, TErr;
    double n, nErr;

    minuit->GetParameter(0, beta, betaErr);
    minuit->GetParameter(1, T, TErr);
    minuit->GetParameter(2, n, nErr);

    cout << "beta_s = " << beta << " +/- " << betaErr << endl;
    cout << "T = " << T << " +/- " << TErr << endl;
    cout << "n = " << n << " +/- " << nErr << endl;

    double betaAvg = 2.0 * beta / (n + 2.0);
}