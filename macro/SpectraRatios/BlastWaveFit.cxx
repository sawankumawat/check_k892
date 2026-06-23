#include <iostream>
#include "TVirtualFitter.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include <cmath>
#include "../src/style.h"

using namespace std;

TH1D *gPi = nullptr;
TH1D *gK = nullptr;
TH1D *gP = nullptr;
TH1D *gKstar = nullptr;

double PiMin = 0.5;
double PiMax = 1.0;

double KMin = 0.3;
double KMax = 1.2;

double PMin = 0.3;
double PMax = 2.5;

double KstarMin = 0.2;
double KstarMax = 1.0;

static TF1 *fBGBlastWave_Integrand = NULL;
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

    // Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
    // if (argI0 > 700.)
    //     argI0 = 700.;
    // Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
    // return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);

    Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
    Double_t argK1 = mt * TMath::CosH(rho) * temp_1;

    Double_t i0, k1;

    // I0
    if (argI0 < 500.)
    {
        i0 = TMath::BesselI0(argI0);
    }
    else
    {
        i0 = TMath::Exp(argI0) / TMath::Sqrt(2.0 * TMath::Pi() * argI0);
    }

    // K1
    if (argK1 < 500.)
    {
        k1 = TMath::BesselK1(argK1);
    }
    else
    {
        k1 = TMath::Sqrt(TMath::Pi() / (2.0 * argK1)) * TMath::Exp(-argK1);
    }

    return r * mt * i0 * k1;
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

double Chi2Hist(TH1D *h, double mass, double beta, double T, double n, double norm, double ptMin, double ptMax)
{
    double chi2 = 0.0;

    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        double pt = h->GetBinCenter(i);

        if (pt < ptMin || pt > ptMax)
            continue;

        double y = h->GetBinContent(i);
        double err = h->GetBinError(i);

        if (err <= 0)
            continue;

        double model = BlastWave(pt, mass, beta, T, n, norm);
        chi2 += pow((y - model) / err, 2);
    }

    return chi2;
}

void GlobalChi2(Int_t &npar, Double_t *gin, Double_t &fval, Double_t *par, Int_t iflag)
{
    double beta = par[0];
    double T = par[1];
    double n = par[2];

    double normPi = par[3];
    double normK = par[4];
    double normP = par[5];

    fval = 0;

    fval += Chi2Hist(gPi, 0.13957, beta, T, n, normPi, PiMin, PiMax);
    fval += Chi2Hist(gK, 0.49367, beta, T, n, normK, KMin, KMax);
    fval += Chi2Hist(gP, 0.93827, beta, T, n, normP, PMin, PMax);
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

// Helper to count points within a specific pT range for a single histogram
int CountPointsInRange(const TH1 *hist, double ptMin, double ptMax);
// Main helper function to calculate total NDF
int CalculateNDF(const TH1 *gPi, const TH1 *gK, const TH1 *gP, double piMin, double piMax, double kMin, double kMax, double pMin, double pMax);

void BlastWaveFit()
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
    TH1D *hSpectraPionSys[10], *hSpectraProtonSys[10], *hSpectraKaonSys[10];
    int multClases[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};

    for (int i = 0; i < 10; i++)
    {
        TFile *fKstarMult = OpenFile(KstarPath + Form("corrected_spectra_%d_%d.root", multClases[i], multClases[i + 1]));
        hSpectraKstar[i] = GetHisto(fKstarMult, Form("mult_%d-%d/corrected_spectra_Integral_final", multClases[i], multClases[i + 1]));
        hSpectraPion[i] = GetHisto(fPion, Form("hCorrectedSpectraStat_%d_%d", multClases[i], multClases[i + 1]));
        hSpectraProton[i] = GetHisto(fProton, Form("hCorrectedSpectraStat_%d_%d", multClases[i], multClases[i + 1]));
        hSpectraKaon[i] = GetHisto(fKaon, Form("hCorrectedSpectraStat_%d_%d", multClases[i], multClases[i + 1]));

        hSpectraPionSys[i] = GetHisto(fPion, Form("hCorrectedSpectraSys_%d_%d", multClases[i], multClases[i + 1]));
        hSpectraProtonSys[i] = GetHisto(fProton, Form("hCorrectedSpectraSys_%d_%d", multClases[i], multClases[i + 1]));
        hSpectraKaonSys[i] = GetHisto(fKaon, Form("hCorrectedSpectraSys_%d_%d", multClases[i], multClases[i + 1]));
    }

    double betaArr[10], betaErrArr[10];
    double betaAvgArr[10], betaAvgErrArr[10];
    double TArr[10], TErrArr[10];
    double nArr[10], nErrArr[10];
    double chi2ndfArr[10];

    TFile *fOutput = new TFile("BlastWaveFits/BlastWaveFitResults.root", "recreate");

    for (int im = 0; im < 10; im++)
    {
        TCanvas *cPiKp = new TCanvas(Form("cPiKp_%d-%d", multClases[im], multClases[im + 1]), "", 720, 720);
        SetCanvasStyle(cPiKp, 0.15, 0.03, 0.03, 0.14);
        gPad->SetLogy();
        hSpectraPion[im]->SetMarkerColor(kBlue + 1);
        hSpectraPion[im]->SetLineColor(kBlue + 1);
        hSpectraPion[im]->SetMarkerStyle(20);
        hSpectraPion[im]->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}");
        hSpectraPion[im]->GetXaxis()->SetTitle("#it{p}_{T} [GeV/c]");
        hSpectraPion[im]->GetXaxis()->SetTitleOffset(1.2);
        hSpectraPion[im]->SetMinimum(9e-4);
        hSpectraPion[im]->SetMaximum(4e2);
        hSpectraPion[im]->SetMarkerSize(1.3);
        hSpectraPion[im]->Draw("PE");
        hSpectraKaon[im]->SetMarkerColor(kRed + 1);
        hSpectraKaon[im]->SetLineColor(kRed + 1);
        hSpectraKaon[im]->SetMarkerStyle(21);
        hSpectraKaon[im]->SetMarkerSize(1.2);
        hSpectraKaon[im]->Draw("PE same");
        hSpectraProton[im]->SetMarkerColor(kGreen + 2);
        hSpectraProton[im]->SetLineColor(kGreen + 2);
        hSpectraProton[im]->SetMarkerStyle(22);
        hSpectraProton[im]->SetMarkerSize(1.4);
        hSpectraProton[im]->Draw("PE same");
        hSpectraKstar[im]->SetMarkerColor(kBrown);
        hSpectraKstar[im]->SetLineColor(kBrown);
        hSpectraKstar[im]->SetMarkerStyle(24);
        hSpectraKstar[im]->SetMarkerSize(1.1);
        // hSpectraKstar[im]->Draw("PE same");

        gPi = hSpectraPion[im];
        gK = hSpectraKaon[im];
        gP = hSpectraProton[im];

        TVirtualFitter *minuit = TVirtualFitter::Fitter(0, 6);

        minuit->SetFCN(GlobalChi2);

        double arglist[10];
        int ierflg = 0;

        minuit->SetParameter(0, "beta_s", 0.8, 0.01, 0.1, 0.99); // parameter index, name, initial value, step size, lower limit, upper limit
        minuit->SetParameter(1, "T", 0.16, 0.001, 0.01, 0.20);
        minuit->SetParameter(2, "n", 1.5, 0.01, 0.01, 15.0);
        minuit->SetParameter(3, "NormPi", 4100, 1, 1, 1e5);
        minuit->SetParameter(4, "NormK", 2500, 1, 1, 1e5);
        minuit->SetParameter(5, "NormP", 9500, 1, 1, 1e5);

        arglist[0] = 10000;
        arglist[1] = 1;

        minuit->ExecuteCommand("MIGRAD", arglist, 2);

        betaArr[im] = minuit->GetParameter(0);
        TArr[im] = minuit->GetParameter(1);
        nArr[im] = minuit->GetParameter(2);
        betaErrArr[im] = minuit->GetParError(0);
        TErrArr[im] = minuit->GetParError(1);
        nErrArr[im] = minuit->GetParError(2);

        betaAvgArr[im] = 2.0 * betaArr[im] / (nArr[im] + 2.0);
        double dBdbeta = 2.0 / (nArr[im] + 2.0);
        double dBdn = -2.0 * betaArr[im] / pow(nArr[im] + 2.0, 2);
        betaAvgErrArr[im] = sqrt(pow(dBdbeta * betaErrArr[im], 2) + pow(dBdn * nErrArr[im], 2));

        // cout << "Blast Wave Fit Results for Mult Class " << multClases[im] << "-" << multClases[im + 1] << ":" << endl;
        // cout << "Norm Pion = " << normPi << endl;
        // cout << "Norm Kaon = " << normK << endl;
        // cout << "Norm Proton = " << normP << endl;
        // cout << "beta_s = " << beta << " +/- " << betaErr << endl;
        // cout << "T = " << T << " +/- " << TErr << endl;
        // cout << "n = " << n << " +/- " << nErr << endl;
        // cout << "beta_avg = " << betaAvg << " +/- " << betaAvgErr << endl;

        double normPi = minuit->GetParameter(3);
        double normK = minuit->GetParameter(4);
        double normP = minuit->GetParameter(5);

        TF1 *fPi = new TF1(Form("fPi_%d-%d", multClases[im], multClases[im + 1]), BGBlastWave_Func, PiMin, 1.2, 5);
        fPi->SetParameters(0.13957, betaArr[im], TArr[im], nArr[im], normPi);

        TF1 *fK = new TF1(Form("fK_%d-%d", multClases[im], multClases[im + 1]), BGBlastWave_Func, 0.3, KMax, 5);
        fK->SetParameters(0.49367, betaArr[im], TArr[im], nArr[im], normK);

        TF1 *fP = new TF1(Form("fP_%d-%d", multClases[im], multClases[im + 1]), BGBlastWave_Func, PMin, PMax, 5);
        fP->SetParameters(0.93827, betaArr[im], TArr[im], nArr[im], normP);

        fPi->SetLineColor(kBlue + 1);
        fK->SetLineColor(kRed + 1);
        fP->SetLineColor(kGreen + 2);

        fPi->SetLineWidth(3);
        fK->SetLineWidth(3);
        fP->SetLineWidth(3);

        fPi->Draw("same");
        fK->Draw("same");
        fP->Draw("same");

        double amin, edm, errdef;
        int nvpar, nparx;
        minuit->GetStats(amin, edm, errdef, nvpar, nparx);

        int ndf = CalculateNDF(gPi, gK, gP, PiMin, PiMax, KMin, KMax, PMin, PMax);
        chi2ndfArr[im] = amin / ndf;

        cout << "chi2/NDF = " << amin << "/" << ndf << " = " << amin / ndf << endl;

        TPaveText *pt = new TPaveText(0.55, 0.65, 0.88, 0.88, "NDC");
        pt->AddText(Form("#chi^{2}/NDF = %.1f/%d = %.2f", amin, ndf, amin / ndf));
        pt->AddText(Form("#beta_{s} = %.3f #pm %.3f", betaArr[im], betaErrArr[im]));
        pt->AddText(Form("T = %.3f #pm %.3f GeV", TArr[im], TErrArr[im]));
        pt->AddText(Form("n = %.3f #pm %.3f", nArr[im], nErrArr[im]));
        pt->AddText(Form("<#beta> = %.3f #pm %.3f", betaAvgArr[im], betaAvgErrArr[im]));
        pt->Draw();
        cPiKp->SaveAs(Form("BlastWaveFits/WithoutKstar/BlastWaveFit_%d-%d.png", multClases[im], multClases[im + 1]));
        cPiKp->Write(Form("BlastWaveFit_%d-%d", multClases[im], multClases[im + 1]));
    }

    double dnch_detaRun3[] = {21.78, 18.48, 15.76, 13.89, 12.50, 10.86, 9.09, 7.63, 5.87, 3.69};
    double dnch_detaRun3_err[] = {0.38, 0.25, 0.22, 0.19, 0.17, 0.15, 0.13, 0.11, 0.09, 0.06};
    double chi2err[10] = {0};
    TGraphErrors *grBetaAvg = new TGraphErrors(10, dnch_detaRun3, betaAvgArr, dnch_detaRun3_err, betaAvgErrArr);
    TGraphErrors *grT = new TGraphErrors(10, dnch_detaRun3, TArr, dnch_detaRun3_err, TErrArr);
    TGraphErrors *grN = new TGraphErrors(10, dnch_detaRun3, nArr, dnch_detaRun3_err, nErrArr);
    TGraphErrors *grChi2 = new TGraphErrors(10, dnch_detaRun3, chi2ndfArr, dnch_detaRun3_err, chi2err);

    TCanvas *cAvgBeta = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cAvgBeta, 0.15, 0.03, 0.03, 0.14);
    SetGraphErrorStyle(grBetaAvg);
    grBetaAvg->SetMarkerColor(kBlue + 1);
    grBetaAvg->SetLineColor(kBlue + 1);
    grBetaAvg->SetMarkerStyle(20);
    grBetaAvg->SetMarkerSize(1.3);
    grBetaAvg->GetXaxis()->SetTitle("#it{d}#it{N}_{ch}/d#it{#eta}");
    grBetaAvg->GetYaxis()->SetTitle("<#beta> (c)");
    grBetaAvg->GetXaxis()->SetTitleOffset(1.2);
    grBetaAvg->SetMinimum(0.1);
    grBetaAvg->SetMaximum(0.6);
    grBetaAvg->Draw("AP");
    grBetaAvg->Write("grBetaAvg");
    cAvgBeta->SaveAs("BlastWaveFits/WithoutKstar/BlastWaveFit_AvgBeta.png");

    TCanvas *cT = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cT, 0.15, 0.03, 0.03, 0.14);
    SetGraphErrorStyle(grT);
    grT->SetMarkerColor(kBlue + 1);
    grT->SetLineColor(kBlue + 1);
    grT->SetMarkerStyle(21);
    grT->SetMarkerSize(1.2);
    grT->GetXaxis()->SetTitle("#it{d}#it{N}_{ch}/d#it{#eta}");
    grT->GetYaxis()->SetTitle("T (GeV)");
    grT->GetXaxis()->SetTitleOffset(1.2);
    grT->SetMinimum(0.12);
    grT->SetMaximum(0.18);
    grT->Draw("AP");
    grT->Write("grT");
    cT->SaveAs("BlastWaveFits/WithoutKstar/BlastWaveFit_T.png");

    TCanvas *cN = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cN, 0.15, 0.03, 0.03, 0.14);
    SetGraphErrorStyle(grN);
    grN->SetMarkerColor(kBlue + 1);
    grN->SetLineColor(kBlue + 1);
    grN->SetMarkerStyle(22);
    grN->SetMarkerSize(1.4);
    grN->GetXaxis()->SetTitle("#it{d}#it{N}_{ch}/d#it{#eta}");
    grN->GetYaxis()->SetTitle("n");
    grN->GetXaxis()->SetTitleOffset(1.2);
    grN->SetMinimum(1.0);
    grN->SetMaximum(6.0);
    grN->Draw("AP");
    grN->Write("grN");
    cN->SaveAs("BlastWaveFits/WithoutKstar/BlastWaveFit_n.png");

    TCanvas *cChi2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cChi2, 0.15, 0.03, 0.03, 0.14);
    SetGraphErrorStyle(grChi2);
    gPad->SetLogy();
    grChi2->SetMarkerColor(kBlue + 1);
    grChi2->SetLineColor(kBlue + 1);
    grChi2->SetMarkerStyle(23);
    grChi2->SetMarkerSize(1.4);
    grChi2->GetXaxis()->SetTitle("#it{d}#it{N}_{ch}/d#it{#eta}");
    grChi2->GetYaxis()->SetTitle("#chi^{2}/NDF");
    grChi2->GetXaxis()->SetTitleOffset(1.2);
    grChi2->SetMinimum(1000.0);
    grChi2->SetMaximum(100000);
    grChi2->Draw("AP");
    grChi2->Write("grChi2");
    cChi2->SaveAs("BlastWaveFits/WithoutKstar/BlastWaveFit_Chi2.png");
}

// Helper to count points within a specific pT range for a single histogram
int CountPointsInRange(const TH1 *hist, double ptMin, double ptMax)
{
    if (!hist)
        return 0; // Safety check for null pointers

    int points = 0;
    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        double pt = hist->GetBinCenter(i);
        if (pt >= ptMin && pt <= ptMax)
        {
            points++;
        }
    }
    return points;
}

// Main helper function to calculate total NDF
int CalculateNDF(const TH1 *gPi, const TH1 *gK, const TH1 *gP, double piMin, double piMax, double kMin, double kMax, double pMin, double pMax)
{
    int totalPoints = 0;

    // Sum up points from pions, kaons, and protons with their respective limits
    totalPoints += CountPointsInRange(gPi, piMin, piMax);
    totalPoints += CountPointsInRange(gK, kMin, kMax);
    totalPoints += CountPointsInRange(gP, pMin, pMax);

    // NDF = Total data points - Number of free parameters (6)
    int ndf = totalPoints - 6;

    return ndf;
}