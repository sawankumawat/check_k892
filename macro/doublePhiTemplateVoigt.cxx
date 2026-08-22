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

double VoigtianShape(double m, double m0, double sigma, double gamma)
{
    return TMath::Voigt(m - m0, sigma, gamma);
}

double FitVoightianpol2(double *x, double *par)
{
    double m = x[0];
    double amp = par[0];
    double m0 = par[1];
    double sigma = par[2];
    double gamma = par[3];

    double pol2 = par[4] + par[5] * m + par[6] * m * m;

    return amp * VoigtianShape(m, m0, sigma, gamma) + pol2;
}

// 2D Fit Function Model
double FitFunc2D(double *x, double *p)
{
    double m1 = x[0];
    double m2 = x[1];
    double mPDG = 1.0198; // PDG mass of phi meson in GeV/c^2

    double n_SS = p[0];
    double n_SB = p[1];
    double n_BS = p[2];
    double n_BB = p[3];

    double m0 = p[4];
    double sigma = p[5]; // Resolution
    double gamma = p[6]; // Width

    double sig1 = VoigtianShape(m1, m0, sigma, gamma);
    double sig2 = VoigtianShape(m2, m0, sigma, gamma);

    // exp(pol3) background function
    auto Bkg = [&](double m)
    {
        return pow(m, p[7]) * TMath::Exp(p[8] * m + p[9] * m * m + p[10] * m * m * m);
    };

    double bkg1 = Bkg(m1);
    double bkg2 = Bkg(m2);

    double shape_SS = sig1 * sig2;
    double shape_SB = sig1 * bkg2;
    double shape_BS = bkg1 * sig2;
    double shape_BB = bkg1 * bkg2;

    return n_SS * shape_SS + n_SB * shape_SB + n_BS * shape_BS + n_BB * shape_BB;
}

void doublePhiTemplateVoigt()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template";

    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults_MorePhiBins.root");

    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassDoublePhi");

    int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowDeltaM = hUnlike->GetAxis(6)->FindBin(0.0 + 0.00001);
    int highDeltaM = hUnlike->GetAxis(6)->FindBin(0.005 - 0.00001);

    hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    TH3D *h3D_full = hUnlike->Projection(0, 4, 5, "E");

    hUnlike->GetAxis(6)->SetRange(lowDeltaM, highDeltaM); // Uncommented: essential for DeltaM selection
    TH3D *h3D_cut = hUnlike->Projection(0, 4, 5, "E");

    TH1D *hInvMass = hUnlike->Projection(0, "E");
    TCanvas *cInvMass = new TCanvas("cInvMass", "Invariant Mass Distribution", 720, 720);
    SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
    hInvMass->Rebin(8);
    hInvMass->Draw("pe");
    cInvMass->SaveAs(savepath + "/PhiInvMass.png");

    hUnlike->GetAxis(1)->SetRange(0, -1); // Reset pT range for further analysis

    // 2. Fix the bin-beating effect by iterating strictly by bin index
    int rebin = 10;
    int startBin = h3D_full->GetXaxis()->FindBin(2.5 + 0.0001);
    int endBin = h3D_full->GetXaxis()->FindBin(2.9 - 0.0001);
    int nBinsInRange = endBin - startBin + 1;
    int totalBins = nBinsInRange / rebin;
    // totalBins = 1; // For testing, set to 1. Remove this line for full analysis.

    // ================================================
    // 1D Histograms for SS and Non-SS Yields
    // ================================================
    TH1D *h_N_SS = new TH1D("h_N_SS", "SS Template (True #phi#phi Yield);#it{M}_{#phi#phi} (GeV/#it{c}^{2});selected N_{SS}", totalBins, 2.5, 2.9);
    TH1D *h_N_nonSS = new TH1D("h_N_nonSS", "Non-SS Template (Background Yield);#it{M}_{#phi#phi} (GeV/#it{c}^{2});selected N_{nonSS}", totalBins, 2.5, 2.9);
    TH1D *h_N_Total = new TH1D("h_N_Total", "Total Yield;#it{M}_{#phi#phi} (GeV/#it{c}^{2});selected N_{Total}", totalBins, 2.5, 2.9);

    for (int ibin = 0; ibin < totalBins; ibin++)
    {
        int lowInvMassBin = startBin + ibin * rebin;
        int highInvMassBin = startBin + (ibin + 1) * rebin - 1;
        h3D_full->GetXaxis()->SetRange(lowInvMassBin, highInvMassBin);

        TH2D *h2D_Mass = (TH2D *)h3D_full->Project3D("yz");
        h2D_Mass->SetName(Form("h2D_Mass_bin%d", ibin));
        if (!h2D_Mass)
        {
            cout << "Error: Failed to project 2D histogram." << endl;
            return;
        }

        // ================================================
        // 2D FIT using BW + Expol background model
        // ================================================
        // FIX: Increased parameter count from 9 to 16
        TF2 *f2D = new TF2(Form("f2D_bin%d", ibin), FitFunc2D, 1.0, 1.04, 1.0, 1.04, 11);

        f2D->SetParameter(0, 0.006); // N_SS
        f2D->SetParameter(1, 0.001); // N_BS
        f2D->SetParameter(2, 0.001); // N_SB
        f2D->SetParameter(3, 0.002); // N_BB

        f2D->FixParameter(4, 1.01983); // Mass peak
        f2D->SetParameter(5, 0.0012);  // Resolution
        f2D->FixParameter(6, 0.0043);  // Width

        // f2D->SetParameter(7,  166.3); //  p0
        // f2D->SetParameter(8, 68.4); //  p1
        // f2D->SetParameter(9, 12.4); //  p2
        // f2D->SetParameter(10, -77.2); //  p3

        f2D->FixParameter(7, 166.3);  // p0
        f2D->FixParameter(8, 68.4);   // p1
        f2D->FixParameter(9, 12.4);   // p2
        f2D->FixParameter(10, -77.2); // p3

        f2D->SetNpx(1000); // Reduced for speed, increase if fit drawing looks jagged
        f2D->SetNpy(1000);

        TVirtualFitter::SetDefaultFitter("Minuit2");
        TVirtualFitter::SetMaxIterations(20000);
        ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");

        h2D_Mass->Fit(f2D, "QERN");                            // Fit without drawing
        TFitResultPtr fitResult = h2D_Mass->Fit(f2D, "QERSN"); // Fit without drawing
        Double_t *par = f2D->GetParameters();                  // Full fit parameters

        int fitStatus = static_cast<int>(fitResult);
        int covQual = fitResult.Get() ? fitResult->CovMatrixStatus() : -1;
        bool fitValid = fitResult.Get() && fitStatus == 0 && covQual >= 2;

        // Print fit quality
        cout << "========== 2D FIT RESULTS ==========" << endl;
        cout << "Bin range " << h3D_full->GetXaxis()->GetBinLowEdge(lowInvMassBin) << " - " << h3D_full->GetXaxis()->GetBinUpEdge(highInvMassBin) << endl;
        cout << "Fit status        : " << fitStatus << endl;
        cout << "Covariance quality: " << covQual << endl;
        cout << "Fit valid         : " << boolalpha << fitValid << endl;
        cout << "Resolution   : " << f2D->GetParameter(5) << " +/- " << f2D->GetParError(5) << endl;
        cout << "N_SS        : " << f2D->GetParameter(0) << " +/- " << f2D->GetParError(0) << endl;
        cout << "N_SB     : " << f2D->GetParameter(1) << " +/- " << f2D->GetParError(1) << endl;
        cout << "N_BS     : " << f2D->GetParameter(2) << " +/- " << f2D->GetParError(2) << endl;
        cout << "N_BB     : " << f2D->GetParameter(3) << " +/- " << f2D->GetParError(3) << endl;
        cout << "Chi2 / NDF  : " << f2D->GetChisquare() / f2D->GetNDF() << endl;
        cout << "Bkg par0   : " << f2D->GetParameter(7) << " +/- " << f2D->GetParError(7) << endl;
        cout << "Bkg par1   : " << f2D->GetParameter(8) << " +/- " << f2D->GetParError(8) << endl;
        cout << "Bkg par2   : " << f2D->GetParameter(9) << " +/- " << f2D->GetParError(9) << endl;
        cout << "Bkg par3   : " << f2D->GetParameter(10) << " +/- " << f2D->GetParError(10) << endl;

        cout << "====================================" << endl;
        cout << endl;

        // =========================================================================
        // STEP 2: Project dataset WITH DeltaM < 0.005 selection (Integrated n_k entries)
        // =========================================================================
        h3D_cut->GetXaxis()->SetRange(lowInvMassBin, highInvMassBin);

        TH2D *h2D_cut = (TH2D *)h3D_cut->Project3D("yz");
        h2D_cut->SetName(Form("h2D_cut_bin%d", ibin));

        // =========================================================================
        // STEP 3: Sum local posterior weights over selected entries (Eq. 5 and 6)
        // =========================================================================
        double yieldSS_cut = 0.0;
        double yieldNonSS_cut = 0.0;

        for (int ix = 1; ix <= h2D_cut->GetNbinsX(); ++ix)
        {
            for (int iy = 1; iy <= h2D_cut->GetNbinsY(); ++iy)
            {
                double n_k = h2D_cut->GetBinContent(ix, iy); // Measured entries in cut region
                if (n_k <= 0)
                    continue;

                double m2 = h2D_cut->GetXaxis()->GetBinCenter(ix);
                double m1 = h2D_cut->GetYaxis()->GetBinCenter(iy);

                // // Compute shapes using full-fit parameters
                // double S1 = BWShape(m1, par[4], par[5]);
                // double S2 = BWShape(m2, par[6], par[7]);
                // double B1 = std::max(0.0, ExpolBkg(m1, par[8], par[9], par[10], par[11]));
                // double B2 = std::max(0.0, ExpolBkg(m2, par[12], par[13], par[14], par[15]));

                double S1 = VoigtianShape(m1, f2D->GetParameter(4), f2D->GetParameter(5), f2D->GetParameter(6));
                double S2 = VoigtianShape(m2, f2D->GetParameter(4), f2D->GetParameter(5), f2D->GetParameter(6));
                double B1 = std::max(0.0, pow(m1, f2D->GetParameter(7)) * TMath::Exp(f2D->GetParameter(8) * m1 + f2D->GetParameter(9) * m1 * m1 + f2D->GetParameter(10) * m1 * m1 * m1));
                double B2 = std::max(0.0, pow(m2, f2D->GetParameter(7)) * TMath::Exp(f2D->GetParameter(8) * m2 + f2D->GetParameter(9) * m2 * m2 + f2D->GetParameter(10) * m2 * m2 * m2));

                double valSS = f2D->GetParameter(0) * S1 * S2;
                double valSB = f2D->GetParameter(1) * S1 * B2;
                double valBS = f2D->GetParameter(2) * B1 * S2;
                double valBB = f2D->GetParameter(3) * B1 * B2;

                double totalF = valSS + valSB + valBS + valBB;
                double P_SS = (totalF > 0) ? (valSS / totalF) : 0.0;

                yieldSS_cut += n_k * P_SS;
                yieldNonSS_cut += n_k * (1.0 - P_SS);
            }
        }

        // Fill the 1D yield histograms
        h_N_SS->SetBinContent(ibin + 1, yieldSS_cut);
        h_N_nonSS->SetBinContent(ibin + 1, yieldNonSS_cut);
        h_N_Total->SetBinContent(ibin + 1, yieldSS_cut + yieldNonSS_cut);
        cout << "Signal component (SS) yield in bin " << ibin << ": " << yieldSS_cut << endl;
        cout << "Background component (Non-SS) yield in bin " << ibin << ": " << yieldNonSS_cut << endl;

        // Draw and save
        TCanvas *c2D = new TCanvas(Form("c2D_bin%d", ibin), "2D Mass Distribution", 720, 720);
        // SetCanvasStyle(c2D, 0.18, 0.2, 0.07, 0.15); // Assuming this is defined in style.h
        // SetHistoQA2D(h2D_Mass);                     // Assuming this is defined in style.h
        h2D_Mass->GetXaxis()->SetTitle("#it{M}_{#phi2} (GeV/#it{c}^{2})");
        h2D_Mass->GetYaxis()->SetTitle("#it{M}_{#phi1} (GeV/#it{c}^{2})");
        h2D_Mass->GetYaxis()->SetTitleOffset(1.9);
        h2D_Mass->GetXaxis()->SetNdivisions(505);
        h2D_Mass->GetXaxis()->SetRangeUser(1.0, 1.04);
        h2D_Mass->GetYaxis()->SetRangeUser(1.0, 1.04);
        h2D_Mass->Draw("colz");
        f2D->Draw("SAME");
        c2D->SaveAs(savepath + Form("/2DFits/Voigt/2D_Mass_bin%d.png", ibin));

        delete c2D;
        delete f2D;
    }

    // ================================================
    // Drawing & Saving 1D SS and Non-SS Templates
    // ================================================
    TCanvas *cSS = new TCanvas("cSS", "SS Template", 720, 720);
    SetCanvasStyle(cSS, 0.15, 0.05, 0.08, 0.12);
    SetHistoQA(h_N_SS);
    h_N_SS->SetMarkerStyle(20);
    h_N_SS->SetMarkerSize(0.8);
    h_N_SS->Draw("pe");
    cSS->SaveAs(savepath + "/SS_Template_Voigt.png");

    TCanvas *cNonSS = new TCanvas("cNonSS", "Non-SS Template", 720, 720);
    SetCanvasStyle(cNonSS, 0.15, 0.05, 0.08, 0.12);
    SetHistoQA(h_N_nonSS);
    h_N_nonSS->SetMarkerStyle(20);
    h_N_nonSS->SetMarkerSize(0.8);
    h_N_nonSS->GetYaxis()->SetRangeUser(1000, 2300);
    h_N_nonSS->Draw("pe");
    cNonSS->SaveAs(savepath + "/NonSS_Template_Voigt.png");

    TCanvas *cTotal = new TCanvas("cTotal", "Total Yield", 720, 720);
    SetCanvasStyle(cTotal, 0.15, 0.05, 0.08, 0.12);
    SetHistoQA(h_N_Total);
    h_N_Total->SetMarkerStyle(20);
    h_N_Total->SetMarkerSize(0.8);
    h_N_Total->Draw("pe");
    cTotal->SaveAs(savepath + "/Total_Yield_Voigt.png");

    TFile *fOutput = new TFile(savepath + "/PhiPhiBkgTemplate_Voigt.root", "RECREATE");
    h_N_SS->Write("h_N_SS");
    h_N_nonSS->Write("h_N_nonSS");
    fOutput->Close();
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