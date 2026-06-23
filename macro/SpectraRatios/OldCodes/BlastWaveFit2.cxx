// ============================================================
//  BlastWaveFit.C
//  Simultaneous Blast Wave fit to pi, K, p pT spectra
//  Extracts: T_kin (kinetic freeze-out temperature)
//            <beta> (mean transverse flow velocity)
//
//  Usage in ROOT:
//    root -l BlastWaveFit.C
//  or:
//    root -l -b -q 'BlastWaveFit.C("yourfile.root")'
//
//  The macro expects three TH1 histograms named:
//    hPi, hK, hP   (1/N d²N/dpT dy vs pT in GeV/c)
//  You can adapt LoadSpectra() to match your file structure.
// ============================================================

#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TMatrixDSym.h>
#include <TVirtualFitter.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <iostream>
#include <cmath>
#include "../src/style.h"
using namespace std;

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

// ============================================================
//  GLOBAL state for TMinuit FCN
// ============================================================
TH1 *gHPi = nullptr, *gHK = nullptr, *gHP = nullptr;

// Particle masses (GeV/c^2)
const double kMpi = 0.13957;
const double kMK = 0.49368;
const double kMp = 0.93827;

// ============================================================
//  Blast Wave integrand
//
//  BW(pT; T, beta_s, n, mass) =
//    pT * Int_0^R r dr * mT * I0(pT*sinh(rho)/T) * K1(mT*cosh(rho)/T)
//
//  rho(r) = tanh^{-1}[beta_s * (r/R)^n]   (flow rapidity profile)
//  <beta>  = 2/(2+n) * beta_s              (mean surface velocity)
//
//  The radial integral is evaluated numerically (Simpson, 50 steps).
// ============================================================
struct BWParams
{
    double pT;
    double T;
    double betaS;
    double n;
    double mass;
};

// integrand over r in [0,1]   (normalised radius xi = r/R)
double BWintegrand(double xi, void *p)
{
    BWParams *par = (BWParams *)p;
    double pT = par->pT;
    double T = par->T;
    double betaS = par->betaS;
    double n = par->n;
    double mass = par->mass;

    double mT = TMath::Sqrt(pT * pT + mass * mass);
    double beta = betaS * TMath::Power(xi, n); // local flow velocity
    if (beta >= 1.0)
        beta = 0.9999;
    double rho = TMath::ATanH(beta); // flow rapidity

    double argI0 = pT * TMath::SinH(rho) / T;
    double argK1 = mT * TMath::CosH(rho) / T;

    // protect against overflow
    double i0 = (argI0 < 500) ? TMath::BesselI0(argI0) : TMath::Exp(argI0) / TMath::Sqrt(2 * TMath::Pi() * argI0);
    double k1;
    if (argK1 < 500)
    {
        k1 = TMath::BesselK1(argK1);
    }
    else
    {
        k1 = TMath::Sqrt(TMath::Pi() / (2 * argK1)) * TMath::Exp(-argK1);
    }

    return xi * mT * i0 * k1;
}

// Simpson integration over xi in [0, 1]
double BWintegrate(double pT, double T, double betaS, double n, double mass,
                   int nsteps = 50)
{
    BWParams par = {pT, T, betaS, n, mass};
    double h = 1.0 / nsteps;
    double sum = 0;
    for (int i = 0; i <= nsteps; ++i)
    {
        double xi = i * h;
        double fxi = BWintegrand(xi, &par);
        double w = (i == 0 || i == nsteps) ? 1 : (i % 2 == 0 ? 2 : 4);
        sum += w * fxi;
    }
    return sum * h / 3.0;
}

// Full BW spectrum: dN/dpT = pT * Norm * Int
// params: [0]=Norm, [1]=T, [2]=betaS, [3]=n, [4]=mass
double BlastWaveFunc(double *x, double *par)
{
    double pT = x[0];
    double norm = par[0];
    double T = par[1];
    double betaS = par[2];
    double n = par[3];
    double mass = par[4];

    if (T <= 0 || betaS <= 0 || betaS >= 1)
        return 0;

    double integral = BWintegrate(pT, T, betaS, n, mass);
    return norm * pT * integral;
}

// ============================================================
//  Chi-square contribution from one histogram given parameters
// ============================================================
double Chi2Species(TH1 *h, double norm, double T, double betaS, double n, double mass,
                   double ptMin, double ptMax)
{
    double chi2 = 0;
    int ndof = 0;
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
    {
        double pt = h->GetBinCenter(ib);
        if (pt < ptMin || pt > ptMax)
            continue;
        double dat = h->GetBinContent(ib);
        double err = h->GetBinError(ib);
        if (err <= 0 || dat <= 0)
            continue;

        double dpar[5] = {norm, T, betaS, n, mass};
        double theory = BlastWaveFunc(&pt, dpar);
        chi2 += (dat - theory) * (dat - theory) / (err * err);
        ++ndof;
    }
    return chi2;
}

// ============================================================
//  TMinuit FCN
//  Parameters:
//    0 = T_kin     (shared)
//    1 = beta_s    (shared surface velocity)
//    2 = n         (shared flow profile exponent)
//    3 = Norm_pi   (free per species)
//    4 = Norm_K    (free per species)
//    5 = Norm_p    (free per species)
// ============================================================
// pT fit ranges for each species (GeV/c)  -- adjust to your data
double gPtMinPi = 0.50, gPtMaxPi = 2.0;
double gPtMinK = 0.30, gPtMaxK = 1.5;
double gPtMinP = 0.30, gPtMaxP = 3.0;

void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    double T = par[0];
    double betaS = par[1];
    double n = par[2];
    double normPi = par[3];
    double normK = par[4];
    double normP = par[5];

    if (T <= 0 || betaS <= 0 || betaS >= 1 || n <= 0)
    {
        f = 1e30;
        return;
    }

    double chi2 = 0;
    chi2 += Chi2Species(gHPi, normPi, T, betaS, n, kMpi, gPtMinPi, gPtMaxPi);
    chi2 += Chi2Species(gHK, normK, T, betaS, n, kMK, gPtMinK, gPtMaxK);
    chi2 += Chi2Species(gHP, normP, T, betaS, n, kMp, gPtMinP, gPtMaxP);
    f = chi2;
}

// ============================================================
//  Convenience: estimate normalisations from histogram integral
// ============================================================
double EstimateNorm(TH1 *h, double T, double betaS, double n, double mass,
                    double ptMin, double ptMax)
{
    // rough integral of BW shape over fit range
    double bwInt = 0;
    int nb = 200;
    double dp = (ptMax - ptMin) / nb;
    for (int i = 0; i < nb; ++i)
    {
        double pt = ptMin + (i + 0.5) * dp;
        double dpar[5] = {1.0, T, betaS, n, mass};
        bwInt += BlastWaveFunc(&pt, dpar) * dp;
    }

    // data integral
    double datInt = 0;
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
    {
        double pt = h->GetBinCenter(ib);
        if (pt < ptMin || pt > ptMax)
            continue;
        if (h->GetBinContent(ib) <= 0)
            continue;
        datInt += h->GetBinContent(ib) * h->GetBinWidth(ib);
    }
    return (bwInt > 0) ? datInt / bwInt : 1e-3;
}

// ============================================================
//  Create dummy spectra for self-contained demo
//  Replace this with your real histograms
// ============================================================
void BlastWaveFitw(TH1D *&hPi, TH1D *&hK, TH1D *&hP)
{
    // "True" parameters we'll try to recover
    double T_true = 0.10;     // GeV
    double betaS_true = 0.80; // surface velocity
    double n_true = 0.90;

    const int NB = 20;
    double ptBins[NB + 1];
    for (int i = 0; i <= NB; ++i)
        ptBins[i] = 0.05 + i * 0.15;

    hPi = new TH1D("hPi", "#pi^{#pm} spectrum;p_{T} (GeV/c);1/N d^{2}N/dp_{T}dy", NB, ptBins);
    hK = new TH1D("hK", "K^{#pm} spectrum;p_{T} (GeV/c);1/N d^{2}N/dp_{T}dy", NB, ptBins);
    hP = new TH1D("hP", "p+#bar{p} spectrum;p_{T} (GeV/c);1/N d^{2}N/dp_{T}dy", NB, ptBins);

    // Fill histograms from "true" BW + Gaussian noise
    auto fillHist = [&](TH1D *h, double mass, double norm, double ptMin, double ptMax)
    {
        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
        {
            double pt = h->GetBinCenter(ib);
            if (pt < ptMin || pt > ptMax)
            {
                h->SetBinContent(ib, 0);
                h->SetBinError(ib, 0);
                continue;
            }
            double dpar[5] = {norm, T_true, betaS_true, n_true, mass};
            double val = BlastWaveFunc(&pt, dpar);
            // 5% statistical noise
            double noise = val * 0.05 * gRandom->Gaus(0, 1);
            double err = val * 0.05;
            h->SetBinContent(ib, TMath::Max(val + noise, 0.0));
            h->SetBinError(ib, err);
        }
    };

    fillHist(hPi, kMpi, 15.0, gPtMinPi, gPtMaxPi);
    fillHist(hK, kMK, 2.5, gPtMinK, gPtMaxK);
    fillHist(hP, kMp, 1.2, gPtMinP, gPtMaxP);
}

// ============================================================
//  MAIN function
// ============================================================
void BlastWaveFit2(const char *filename = "")
{

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // --------------------------------------------------------
    //  1. Load spectra
    // --------------------------------------------------------
    TH1D *hPi = nullptr, *hK = nullptr, *hP = nullptr;

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
    hSpectraPion[0]->GetYaxis()->SetTitle("d^{2}N/(2#pi #it{p}_{T}dp_{T}dy) [(GeV/c)^{-2}]");
    hSpectraPion[0]->GetXaxis()->SetTitle("#it{p}_{T} [GeV/c]");
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

    hPi = hSpectraPion[0];
    hK = hSpectraKaon[0];
    hP = hSpectraProton[0];

    // Set globals for FCN
    gHPi = hPi;
    gHK = hK;
    gHP = hP;

    // --------------------------------------------------------
    //  2. Initial parameter guesses
    // --------------------------------------------------------
    double T0 = 0.100;     // T_kin start (GeV)
    double betaS0 = 0.700; // surface beta start
    double n0 = 1.000;     // profile exponent start

    double normPi0 = EstimateNorm(hPi, T0, betaS0, n0, kMpi, gPtMinPi, gPtMaxPi);
    double normK0 = EstimateNorm(hK, T0, betaS0, n0, kMK, gPtMinK, gPtMaxK);
    double normP0 = EstimateNorm(hP, T0, betaS0, n0, kMp, gPtMinP, gPtMaxP);

    std::cout << "\nInitial norms: Npi=" << normPi0 << "  NK=" << normK0 << "  Np=" << normP0 << std::endl;

    // --------------------------------------------------------
    //  3. TMinuit setup
    // --------------------------------------------------------
    TMinuit *minuit = new TMinuit(6);
    minuit->SetFCN(FCN);
    minuit->SetPrintLevel(1); // set to -1 to suppress printout

    double arglist[10];
    Int_t ierflg = 0;

    // Error definition: 1 = chi^2
    arglist[0] = 1.0;
    minuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Define parameters: (index, name, start, step, min, max)
    minuit->mnparm(0, "T_kin", T0, 0.005, 0.01, 0.30, ierflg);
    minuit->mnparm(1, "beta_s", betaS0, 0.01, 0.10, 0.999, ierflg);
    minuit->mnparm(2, "n", n0, 0.05, 0.1, 3.0, ierflg);
    minuit->mnparm(3, "Norm_pi", normPi0, normPi0 * 0.1, 0, 0, ierflg);
    minuit->mnparm(4, "Norm_K", normK0, normK0 * 0.1, 0, 0, ierflg);
    minuit->mnparm(5, "Norm_p", normP0, normP0 * 0.1, 0, 0, ierflg);

    // --------------------------------------------------------
    //  4. Minimise
    // --------------------------------------------------------
    // MIGRAD for minimum
    arglist[0] = 5000; // max calls
    arglist[1] = 0.1;  // tolerance
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // MINOS for asymmetric errors
    minuit->mnexcm("MINOS", arglist, 0, ierflg);

    // --------------------------------------------------------
    //  5. Extract results
    // --------------------------------------------------------
    double T_fit, T_err, betaS_fit, betaS_err, n_fit, n_err;
    double normPi_fit, normPi_err, normK_fit, normK_err, normP_fit, normP_err;

    minuit->GetParameter(0, T_fit, T_err);
    minuit->GetParameter(1, betaS_fit, betaS_err);
    minuit->GetParameter(2, n_fit, n_err);
    minuit->GetParameter(3, normPi_fit, normPi_err);
    minuit->GetParameter(4, normK_fit, normK_err);
    minuit->GetParameter(5, normP_fit, normP_err);

    // Mean beta: <beta> = 2/(2+n) * beta_s   (for power-law profile)
    double betaMean = 2.0 / (2.0 + n_fit) * betaS_fit;
    // Propagate error (simplified, ignoring correlation)
    double dBdBs = 2.0 / (2.0 + n_fit);
    double dBdn = -2.0 * betaS_fit / ((2.0 + n_fit) * (2.0 + n_fit));
    double betaMean_err = TMath::Sqrt(dBdBs * dBdBs * betaS_err * betaS_err +
                                      dBdn * dBdn * n_err * n_err);

    // chi2 and NDF
    double chi2_final;
    int ndf_total;
    double edm, errdef;
    int nvpar, nparx;
    minuit->mnstat(chi2_final, edm, errdef, nvpar, nparx, ierflg);

    // Count NDF
    auto countBins = [](TH1 *h, double ptMin, double ptMax)
    {
        int n = 0;
        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
        {
            double pt = h->GetBinCenter(ib);
            if (pt >= ptMin && pt <= ptMax && h->GetBinContent(ib) > 0 && h->GetBinError(ib) > 0)
                ++n;
        }
        return n;
    };
    int nbins_total = countBins(hPi, gPtMinPi, gPtMaxPi) + countBins(hK, gPtMinK, gPtMaxK) + countBins(hP, gPtMinP, gPtMaxP);
    ndf_total = nbins_total - 6; // 6 free parameters

    std::cout << "\n========================================" << std::endl;
    std::cout << "  BLAST WAVE FIT RESULTS" << std::endl;
    std::cout << "========================================" << std::endl;
    std::printf("  T_kin      = %.4f +/- %.4f  GeV\n", T_fit, T_err);
    std::printf("  beta_s     = %.4f +/- %.4f\n", betaS_fit, betaS_err);
    std::printf("  n          = %.4f +/- %.4f\n", n_fit, n_err);
    std::printf("  <beta>     = %.4f +/- %.4f\n", betaMean, betaMean_err);
    std::printf("  chi2/NDF   = %.2f / %d = %.3f\n", chi2_final, ndf_total,
                (ndf_total > 0) ? chi2_final / ndf_total : 0);
    std::printf("  Norm_pi    = %.4e +/- %.4e\n", normPi_fit, normPi_err);
    std::printf("  Norm_K     = %.4e +/- %.4e\n", normK_fit, normK_err);
    std::printf("  Norm_p     = %.4e +/- %.4e\n", normP_fit, normP_err);
    std::cout << "========================================\n"
              << std::endl;

    // --------------------------------------------------------
    //  6. Build TF1 curves for plotting
    // --------------------------------------------------------
    auto makeFitCurve = [&](double norm, double mass, const char *name,
                            double ptMin, double ptMax) -> TF1 *
    {
        TF1 *f1 = new TF1(name, BlastWaveFunc, ptMin, ptMax, 5);
        f1->SetParameters(norm, T_fit, betaS_fit, n_fit, mass);
        f1->SetLineWidth(2);
        return f1;
    };

    TF1 *fPi = makeFitCurve(normPi_fit, kMpi, "fPi", gPtMinPi, gPtMaxPi);
    TF1 *fK = makeFitCurve(normK_fit, kMK, "fK", gPtMinK, gPtMaxK);
    TF1 *fP = makeFitCurve(normP_fit, kMp, "fP", gPtMinP, gPtMaxP);

    fPi->SetLineColor(kBlue + 1);
    fK->SetLineColor(kGreen + 2);
    fP->SetLineColor(kRed + 1);

    hPi->SetMarkerColor(kBlue + 1);
    hPi->SetLineColor(kBlue + 1);
    hPi->SetMarkerStyle(20);
    hK->SetMarkerColor(kGreen + 2);
    hK->SetLineColor(kGreen + 2);
    hK->SetMarkerStyle(21);
    hP->SetMarkerColor(kRed + 1);
    hP->SetLineColor(kRed + 1);
    hP->SetMarkerStyle(22);

    // --------------------------------------------------------
    //  7. Draw
    // --------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "Blast Wave Fit", 720, 720);
    c1->SetLogy();
    SetCanvasStyle(c1, 0.12, 0.03, 0.03, 0.14);

    // Find global y-range
    double ymax = TMath::Max(hPi->GetMaximum(), TMath::Max(hK->GetMaximum(), hP->GetMaximum()));
    double ymin = 1e-5 * ymax;

    hPi->GetYaxis()->SetRangeUser(ymin, ymax * 5);
    hPi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPi->GetYaxis()->SetTitle("1/N d^{2}N/dp_{T}dy (GeV/c)^{-1}");
    hPi->SetTitle("Simultaneous Blast Wave Fit");

    hPi->Draw("E1");
    hK->Draw("E1 same");
    hP->Draw("E1 same");
    fPi->Draw("same");
    fK->Draw("same");
    fP->Draw("same");

    TLegend *leg = new TLegend(0.60, 0.62, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hPi, "#pi^{#pm}", "ep");
    leg->AddEntry(hK, "K^{#pm}", "ep");
    leg->AddEntry(hP, "p+#bar{p}", "ep");
    leg->AddEntry(fPi, "BW fit", "l");
    leg->Draw();

    // Results box
    TPaveText *pt = new TPaveText(0.15, 0.12, 0.52, 0.38, "NDC");
    pt->SetBorderSize(1);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.033);
    pt->AddText(Form("T_{kin} = %.0f #pm %.0f MeV", T_fit * 1000, T_err * 1000));
    pt->AddText(Form("#LT#beta#GT = %.3f #pm %.3f", betaMean, betaMean_err));
    pt->AddText(Form("#beta_{s} = %.3f #pm %.3f", betaS_fit, betaS_err));
    pt->AddText(Form("n = %.2f #pm %.2f", n_fit, n_err));
    pt->AddText(Form("#chi^{2}/NDF = %.1f/%d = %.2f",
                     chi2_final, ndf_total, (ndf_total > 0) ? chi2_final / ndf_total : 0));
    pt->Draw();

    c1->Update();
    // c1->SaveAs("BlastWaveFit.pdf");

    // --------------------------------------------------------
    //  8. Correlation matrix (optional)
    // --------------------------------------------------------
    TCanvas *c2 = new TCanvas("c2", "Correlation matrix", 600, 500);
    TMatrixDSym corr(6);
    minuit->mnemat(corr.GetMatrixArray(), 6);

    // Print correlations for shared parameters
    std::cout << "\nCorrelation matrix (T, betaS, n):" << std::endl;
    double cov[36];
    minuit->mnemat(cov, 6);
    double sigT = T_err, sigBs = betaS_err, sigN = n_err;
    if (sigT > 0 && sigBs > 0)
        std::printf("  rho(T, betaS) = %.3f\n", cov[1] / (sigT * sigBs));
    if (sigT > 0 && sigN > 0)
        std::printf("  rho(T, n)     = %.3f\n", cov[2] / (sigT * sigN));
    if (sigBs > 0 && sigN > 0)
        std::printf("  rho(betaS, n) = %.3f\n", cov[7] / (sigBs * sigN));

    std::cout << "\nDone. Output saved to BlastWaveFit.pdf" << std::endl;
}