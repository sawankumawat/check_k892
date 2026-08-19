#include <iostream>
#include <cmath>
#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"

#include "src/style.h"
#include "src/fitfunc.h"

TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const string &name);

double BWShape(double m, double m0, double gamma)
{
    double denominator = (m - m0) * (m - m0) + gamma * gamma / 4.0;
    return gamma / (TMath::Pi() * 2 * denominator);
}

double FitFunc2DBW(double *x, double *p)
{
    double m1 = x[0];
    double m2 = x[1];
    double mPDG = 1.0198; // PDG mass of phi meson in GeV/c^2

    double n_SS = p[0];
    double n_nonSS = p[1];

    double m0 = p[2];
    double gamma = p[3];

    double sig1 = BWShape(m1, m0, gamma);
    double sig2 = BWShape(m2, m0, gamma);

    // //pol3 background function
    // auto Bkg = [&](double m)
    // {
    //     double z = m - mPDG;
    //     return p[4] + p[5] * z + p[6] * z * z;
    // };

    // exp(pol3) background function
    auto Bkg = [&](double m)
    {
        double z = m - mPDG;
        return pow(m, p[4]) * TMath::Exp(p[5] * m + p[6] * m * m + p[7] * m * m * m);
        // return pow(z, p[4]) * TMath::Exp(p[5] * z + p[6] * z * z + p[7] * z * z * z);
    };

    double bkg1 = Bkg(m1);
    double bkg2 = Bkg(m2);

    double shape_SS = sig1 * sig2;

    // Merge non-SS background shape into a single component
    double shape_nonSS = sig1 * bkg2 + bkg1 * sig2 + bkg1 * bkg2;
    return n_SS * shape_SS + n_nonSS * shape_nonSS;
}

void doublePhiTemplateBWpol3()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template";

    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults_MorePhiBins.root");
    if (!fInput)
        return;

    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassDoublePhi");
    if (!hUnlike)
        return;

    int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowDeltaM = hUnlike->GetAxis(6)->FindBin(0.0 + 0.00001);
    int highDeltaM = hUnlike->GetAxis(6)->FindBin(0.005 - 0.00001);

    hUnlike->GetAxis(1)->SetRange(lowpT, highpT);

    TH3D *h3D = hUnlike->Projection(0, 4, 5, "E");
    int rebin = 12;

    int totalBins = h3D->GetNbinsX() / rebin;
    // totalBins = 1;                             // For testing, only process the first bin
    double interval = (2.9 - 2.5) / totalBins; // Calculate the interval for each bin

    double last_pars[8] = {0.0};
    bool has_valid_seed = false;

    // ================================================
    // 1D Histograms for SS and Non-SS Yields
    // ================================================
    TH1D *h_N_SS = new TH1D("h_N_SS", "SS Template (True #phi#phi Yield);#it{M}_{#phi#phi} (GeV/#it{c}^{2});selected N_{SS}", totalBins, 2.5, 2.9);
    TH1D *h_N_nonSS = new TH1D("h_N_nonSS", "Non-SS Template (Background Yield);#it{M}_{#phi#phi} (GeV/#it{c}^{2});selected N_{nonSS}", totalBins, 2.5, 2.9);

    for (int ibin = 0; ibin < totalBins; ibin++)
    {
        h3D->GetXaxis()->SetRange(0, -1); // Reset the range for each iteration

        double massLow = 2.5 + ibin * interval + 0.00001;
        double massHigh = 2.5 + (ibin + 1) * interval - 0.00001;

        // Exclude region 2.65 to 2.73 (signal region)
        if (massLow > 2.65 && massHigh < 2.73)
        {
            h_N_SS->SetBinContent(ibin + 1, 0);
            h_N_SS->SetBinError(ibin + 1, 0);
            h_N_nonSS->SetBinContent(ibin + 1, 0);
            h_N_nonSS->SetBinError(ibin + 1, 0);
            continue;
        }

        int lowInvMassBin = h3D->GetXaxis()->FindBin(massLow);
        int highInvMassBin = h3D->GetXaxis()->FindBin(massHigh);
        h3D->GetXaxis()->SetRange(lowInvMassBin, highInvMassBin);

        TH2D *h2D_Mass = (TH2D *)h3D->Project3D("yz");
        h2D_Mass->SetName(Form("h2D_Mass_bin%d", ibin));
        if (!h2D_Mass)
        {
            cout << "Error: Failed to project 2D histogram from 3D histogram." << endl;
            return;
        }

        double binWidthX = h2D_Mass->GetXaxis()->GetBinWidth(1);
        double binWidthY = h2D_Mass->GetYaxis()->GetBinWidth(1);
        double binArea = binWidthX * binWidthY;
        double totalIntegral = h2D_Mass->Integral();

        // ================================================
        // 2D FIT using BW + pol2 background model
        // ================================================
        TF2 *f2D = new TF2(Form("f2D_bin%d", ibin), FitFunc2DBW, 1.0, 1.04, 1.0, 1.04, 8);

        // // Seed Yields (scaled by bin area)
        // f2D->SetParameter(0, (0.5 * totalIntegral) * binArea); // N_SS seed
        // f2D->SetParLimits(0, 0.0, totalIntegral * binArea * 10.0);

        // f2D->SetParameter(1, (0.5 * totalIntegral) * binArea); // N_nonSS seed
        // f2D->SetParLimits(1, 0.0, totalIntegral * binArea * 10.0);

        // ================================================
        // SEQUENTIAL SEEDING LOGIC
        // ================================================
        // if (!has_valid_seed)
        {
            // Default initial guesses (for bin 0 or after a fit failure)
            f2D->SetParameter(0, (0.5 * totalIntegral) * binArea); // N_SS
            f2D->SetParameter(1, (0.5 * totalIntegral) * binArea); // N_nonSS

            f2D->SetParameter(2, 1.019); // Mass peak
            f2D->SetParLimits(2, 1.016, 1.025);
            f2D->SetParameter(3, 0.00425); // Width
            f2D->SetParLimits(3, 0.003, 0.009);

            // f2D->SetParameter(4, 1.0); // Pol2 p0
            // f2D->SetParameter(5, 1.0); // Pol2 p1
            // f2D->SetParameter(6, 1.0); // Pol2 p2
            // f2D->SetParameter(7, 1.0); // Pol2 p3

            f2D->FixParameter(4, 166.3); // Pol2 p0
            f2D->FixParameter(5, 68.4); // Pol2 p1
            f2D->FixParameter(6, 12.4); // Pol2 p2
            f2D->FixParameter(7, -77.2); // Pol2 p3
        }
        // else
        // {
        //     // Seed parameters from the previous valid converged bin
        //     for (int p = 0; p < 7; p++)
        //     {
        //         f2D->SetParameter(p, last_pars[p]);
        //     }
        //     // Preserve parameter limits for signal mass & width
        //     f2D->SetParLimits(2, 1.016, 1.025);
        //     f2D->SetParLimits(3, 0.003, 0.009);
        // }

        f2D->SetNpx(1000);
        f2D->SetNpy(1000);

        TVirtualFitter::SetDefaultFitter("Minuit2");
        TVirtualFitter::SetMaxIterations(20000);
        ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");

        // Perform the 2D fit and store result
        h2D_Mass->Fit(f2D, "Q0ERN");
        TFitResultPtr fitResult = h2D_Mass->Fit(f2D, "Q0ERSN");
        // TFitResultPtr fitResult = h2D_Mass->Fit(f2D, "0REBMS");

        // Convert fitted density parameters back to integrated bin counts
        double n_SS = f2D->GetParameter(0) / binArea;
        double err_SS = f2D->GetParError(0) / binArea;

        double n_nonSS = f2D->GetParameter(1) / binArea;
        double err_nonSS = f2D->GetParError(1) / binArea;

        // Fill the 1D pair-mass histograms
        h_N_SS->SetBinContent(ibin + 1, n_SS);
        h_N_SS->SetBinError(ibin + 1, err_SS);

        h_N_nonSS->SetBinContent(ibin + 1, n_nonSS);
        h_N_nonSS->SetBinError(ibin + 1, err_nonSS);

        int fitStatus = static_cast<int>(fitResult);
        int covQual = fitResult.Get() ? fitResult->CovMatrixStatus() : -1;
        bool fitValid = fitResult.Get() && fitStatus == 0 && covQual >= 2;

        // ================================================
        // UPDATE SEED STORE ONLY ON VALID FITS
        // ================================================
        if (fitValid)
        {
            for (int p = 0; p < 8; p++)
            {
                last_pars[p] = f2D->GetParameter(p);
            }
            has_valid_seed = true;
        }

        // Print fit quality
        cout << "========== 2D FIT RESULTS ==========" << endl;
        cout << "Bin range " << Form("%.2f - %.2f", 2.5 + ibin * interval, 2.5 + (ibin + 1) * interval) << endl;
        cout << "Fit status        : " << fitStatus << endl;
        cout << "Covariance quality: " << covQual << endl;
        cout << "Fit valid         : " << boolalpha << fitValid << endl;
        cout << "Mass peak   : " << f2D->GetParameter(2) << " +/- " << f2D->GetParError(2) << endl;
        cout << "Width       : " << f2D->GetParameter(3) << " +/- " << f2D->GetParError(3) << endl;
        cout << "N_SS        : " << f2D->GetParameter(0) << " +/- " << f2D->GetParError(0) << endl;
        cout << "N_nonSS     : " << f2D->GetParameter(1) << " +/- " << f2D->GetParError(1) << endl;
        cout << "Chi2        : " << f2D->GetChisquare() << endl;
        cout << "NDF         : " << f2D->GetNDF() << endl;
        cout << "Chi2 / NDF  : " << f2D->GetChisquare() / f2D->GetNDF() << endl;
        cout << "Bkg par0   : " << f2D->GetParameter(4) << " +/- " << f2D->GetParError(4) << endl;
        cout << "Bkg par1   : " << f2D->GetParameter(5) << " +/- " << f2D->GetParError(5) << endl;
        cout << "Bkg par2   : " << f2D->GetParameter(6) << " +/- " << f2D->GetParError(6) << endl;
        cout << "Bkg par3   : " << f2D->GetParameter(7) << " +/- " << f2D->GetParError(7) << endl;

        cout << "====================================" << endl;
        cout << endl;

        // Draw 2D mass distribution with 2D fit overlay
        TCanvas *c2D = new TCanvas(Form("c2D_bin%d", ibin), "2D Mass Distribution", 720, 720);
        SetCanvasStyle(c2D, 0.18, 0.2, 0.07, 0.15);
        SetHistoQA2D(h2D_Mass);
        h2D_Mass->GetXaxis()->SetTitle("#it{M}_{#phi2} (GeV/#it{c}^{2})");
        h2D_Mass->GetYaxis()->SetTitle("#it{M}_{#phi1} (GeV/#it{c}^{2})");
        h2D_Mass->GetYaxis()->SetTitleOffset(1.9);
        h2D_Mass->GetXaxis()->SetNdivisions(505);
        h2D_Mass->GetXaxis()->SetRangeUser(1.0, 1.04);
        h2D_Mass->GetYaxis()->SetRangeUser(1.0, 1.04);
        h2D_Mass->Draw("colz");
        f2D->Draw("SAME");
        c2D->SaveAs(savepath + Form("/2D_Mass_bin%d.png", ibin));
    }

    // ================================================
    // Drawing & Saving 1D SS and Non-SS Templates
    // ================================================
    TCanvas *cSS = new TCanvas("cSS", "SS Template", 800, 600);
    SetCanvasStyle(cSS, 0.15, 0.05, 0.08, 0.12);
    SetHistoQA(h_N_SS);
    h_N_SS->SetMarkerStyle(20);
    h_N_SS->SetMarkerSize(0.8);
    // h_N_SS->SetMinimum(h_N_SS->GetMaximum() * 0.1);
    // h_N_SS->SetMaximum(h_N_SS->GetMaximum() * 1.5);
    h_N_SS->Draw("E1");
    cSS->SaveAs(savepath + "/SS_Template.png");

    TCanvas *cNonSS = new TCanvas("cNonSS", "Non-SS Template", 800, 600);
    SetCanvasStyle(cNonSS, 0.15, 0.05, 0.08, 0.12);
    SetHistoQA(h_N_nonSS);
    h_N_nonSS->SetMarkerStyle(20);
    h_N_nonSS->SetMarkerSize(0.8);
    // h_N_nonSS->SetMinimum(h_N_nonSS->GetMaximum() * 0.1);
    // h_N_nonSS->SetMaximum(h_N_nonSS->GetMaximum() * 1.5);
    h_N_nonSS->Draw("E1");
    cNonSS->SaveAs(savepath + "/NonSS_Template.png");

    TFile *fOutput = new TFile(savepath + "/DoublePhiBackgroundTemplates.root", "RECREATE");
    h_N_SS->Write("h_N_SS");
    h_N_nonSS->Write("h_N_nonSS");

    //========================================================================
    //=======To visualize 1D phi projections with the fits=================
    //========================================================================

    // TH1D *hPhi1 = h2D_Mass->ProjectionY("hPhi1", 1, h2D_Mass->GetNbinsX(), "E");
    // TH1D *hPhi2 = h2D_Mass->ProjectionX("hPhi2", 1, h2D_Mass->GetNbinsY(), "E");

    // TCanvas *c1D = new TCanvas("c1D", "c1D", 720, 720);
    // SetCanvasStyle(c1D, 0.18, 0.13, 0.07, 0.15);
    // SetHistoQA(hPhi1);
    // hPhi1->GetXaxis()->SetTitle("#it{M}_{#phi1} (GeV/#it{c}^{2})");
    // hPhi1->GetYaxis()->SetTitle("Counts");
    // hPhi1->Draw("pe");
    // TF1 *fPhi1 = new TF1("fPhi1", FitVoightianpol2, 1.0, 1.04, 7);
    // fPhi1->FixParameter(1, f2D->GetParameter(3));
    // fPhi1->FixParameter(2, f2D->GetParameter(4));
    // fPhi1->FixParameter(3, f2D->GetParameter(5));
    // fPhi1->SetParameter(0, 100);
    // fPhi1->SetParameter(4, 1.0);
    // fPhi1->SetParameter(5, 1.0);
    // fPhi1->SetParameter(6, 1.0);
    // hPhi1->Fit(fPhi1, "REBMS");

    // TCanvas *c1D_2 = new TCanvas("c1D_2", "c1D_2", 720, 720);
    // SetCanvasStyle(c1D_2, 0.18, 0.13, 0.07, 0.15);
    // SetHistoQA(hPhi2);
    // hPhi2->GetXaxis()->SetTitle("#it{M}_{#phi2} (GeV/#it{c}^{2})");
    // hPhi2->GetYaxis()->SetTitle("Counts");
    // hPhi2->Draw("pe");
}

TFile *OpenFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (!f || f->IsZombie())
    {
        cout << "Error: File not found: " << path << endl;
        return nullptr;
    }
    return f;
}

template <typename T>
T *GetHisto(TFile *f, const string &name)
{
    T *histo = dynamic_cast<T *>(f->Get(name.c_str()));
    if (!histo)
    {
        cout << "Error: histo " << name
             << " not found in file " << f->GetName() << endl;
        return nullptr;
    }
    return histo;
}