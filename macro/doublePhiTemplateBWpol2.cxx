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
    double n_SB = p[1];
    double n_BB = p[2];

    double m0 = p[3];
    double gamma = p[4];

    double sig1 = BWShape(m1, m0, gamma);
    double sig2 = BWShape(m2, m0, gamma);

    auto Bkg = [&](double m)
    {
        double z = m - mPDG;
        return p[5] + p[6] * z + p[7] * z * z;
    };

    double bkg1 = Bkg(m1);
    double bkg2 = Bkg(m2);

    double shape_SS = sig1 * sig2;
    double shape_SB = sig1 * bkg2 + bkg1 * sig2;
    double shape_BB = bkg1 * bkg2;

    return n_SS * shape_SS + n_SB * shape_SB + n_BB * shape_BB;
}

void doublePhiTemplateBWpol2()
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
    // hUnlike->GetAxis(6)->SetRange(lowDeltaM, highDeltaM);

    TH3D *h3D = hUnlike->Projection(0, 4, 5, "E");
    int rebin = 10;

    int totalBins = h3D->GetNbinsX() / rebin;
    // totalBins = 1;                             // For testing, only process the first bin
    double interval = (2.9 - 2.5) / totalBins; // Calculate the interval for each bin

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

        // // Exclude region 2.65 to 2.73 (signal region)
        // if (massLow > 2.65 && massHigh < 2.73)
        // {
        //     h_N_SS->SetBinContent(ibin + 1, 0);
        //     h_N_SS->SetBinError(ibin + 1, 0);
        //     h_N_nonSS->SetBinContent(ibin + 1, 0);
        //     h_N_nonSS->SetBinError(ibin + 1, 0);
        //     continue;
        // }

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
        // Default initial guesses (for bin 0 or after a fit failure)
        f2D->SetParameter(0, 1.0); // N_SS
        f2D->SetParameter(1, 1.0); // N_BS
        f2D->SetParameter(2, 1.0); // N_BB

        f2D->FixParameter(3, 1.01985); // Mass peak
        f2D->FixParameter(4, 0.0072); // Width

        // f2D->SetParameter(3, 1.01983); // Mass peak
        // f2D->SetParLimits(3, 1.016, 1.025);
        // f2D->SetParameter(4, 0.007077); // Width
        // f2D->SetParLimits(4, 0.003, 0.009);

        // f2D->SetParameter(5, 137.3); // Pol2 p0
        // f2D->SetParameter(6, 2077.1);  // Pol2 p1
        // f2D->SetParameter(7, -10778.7);  // Pol2 p2

        f2D->FixParameter(5, 137.3); // Pol2 p0
        f2D->FixParameter(6, 2077.1);  // Pol2 p1
        f2D->FixParameter(7, -10778.7);  // Pol2 p2

        f2D->SetNpx(1000);
        f2D->SetNpy(1000);

        TVirtualFitter::SetDefaultFitter("Minuit2");
        TVirtualFitter::SetMaxIterations(20000);
        ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");

        // Perform the 2D fit and store result
        h2D_Mass->Fit(f2D, "Q0ERN");
        TFitResultPtr fitResult = h2D_Mass->Fit(f2D, "Q0ERSN");

        // Convert fitted density parameters back to integrated bin counts
        double n_SS = f2D->GetParameter(0) / binArea;
        double err_SS = f2D->GetParError(0) / binArea;

        double n_SB = f2D->GetParameter(1) / binArea;
        double err_SB = f2D->GetParError(1) / binArea;

        double n_BB = f2D->GetParameter(2) / binArea;
        double err_BB = f2D->GetParError(2) / binArea;

        double m0 = f2D->GetParameter(3);
        double gamma = f2D->GetParameter(4);
        double p5 = f2D->GetParameter(5), p6 = f2D->GetParameter(6);
        double p7 = f2D->GetParameter(7), p8 = f2D->GetParameter(8);

        // same 1D shapes as inside FitFunc2DBW
        TF1 fSig("fSig", [=](double *x, double *)
                 { return BWShape(x[0], m0, gamma); }, 1.0, 1.04, 0);
        TF1 fBkg("fBkg", [=](double *x, double *)
                 { return p5 + p6 * (x[0] - 1.0198) + p7 * (x[0] - 1.0198) * (x[0] - 1.0198); }, 1.0, 1.04, 0);
        // the signal box: |m - mPDG| < 0.005
        double winLo = 1.01983 - 0.005, winHi = 1.01983 + 0.005;

        double Isig_win = fSig.Integral(winLo, winHi);
        double Ibkg_win = fBkg.Integral(winLo, winHi);

        // yields WITHIN the deltaM<0.005 box, obtained by integrating the fitted 2D pdf
        double N_SS_win = n_SS * Isig_win * Isig_win;
        double N_SB_win = n_SB * 2.0 * Isig_win * Ibkg_win; // sig1*bkg2 + bkg1*sig2
        double N_BB_win = n_BB * Ibkg_win * Ibkg_win;

        h_N_SS->SetBinContent(ibin + 1, N_SS_win);
        h_N_SS->SetBinError(ibin + 1, err_SS * Isig_win * Isig_win);

        h_N_nonSS->SetBinContent(ibin + 1, N_SB_win + N_BB_win);
        h_N_nonSS->SetBinError(ibin + 1, sqrt(err_SB * err_SB * 4.0 * Isig_win * Isig_win * Ibkg_win * Ibkg_win + err_BB * err_BB * Ibkg_win * Ibkg_win));

        int fitStatus = static_cast<int>(fitResult);
        int covQual = fitResult.Get() ? fitResult->CovMatrixStatus() : -1;
        bool fitValid = fitResult.Get() && fitStatus == 0 && covQual >= 2;

        // Print fit quality
        cout << "========== 2D FIT RESULTS ==========" << endl;
        cout << "Bin range " << Form("%.2f - %.2f", 2.5 + ibin * interval, 2.5 + (ibin + 1) * interval) << endl;
        cout << "Fit status        : " << fitStatus << endl;
        cout << "Covariance quality: " << covQual << endl;
        cout << "Fit valid         : " << boolalpha << fitValid << endl;
        cout << "Mass peak   : " << f2D->GetParameter(3) << " +/- " << f2D->GetParError(3) << endl;
        cout << "Width       : " << f2D->GetParameter(4) << " +/- " << f2D->GetParError(4) << endl;
        cout << "N_SS        : " << f2D->GetParameter(0) << " +/- " << f2D->GetParError(0) << endl;
        cout << "N_SB     : " << f2D->GetParameter(1) << " +/- " << f2D->GetParError(1) << endl;
        cout << "N_BB     : " << f2D->GetParameter(2) << " +/- " << f2D->GetParError(2) << endl;
        cout << "Chi2        : " << f2D->GetChisquare() << endl;
        cout << "NDF         : " << f2D->GetNDF() << endl;
        cout << "Chi2 / NDF  : " << f2D->GetChisquare() / f2D->GetNDF() << endl;
        cout << "Bkg par0   : " << f2D->GetParameter(5) << " +/- " << f2D->GetParError(5) << endl;
        cout << "Bkg par1   : " << f2D->GetParameter(6) << " +/- " << f2D->GetParError(6) << endl;
        cout << "Bkg par2   : " << f2D->GetParameter(7) << " +/- " << f2D->GetParError(7) << endl;

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
        c2D->SaveAs(savepath + Form("/2DFits/2D_Mass_bin%d.png", ibin));

        //==============================================
        //======To visualize fit using 1D phi===========
        //==============================================

        // ================================================
        // 1D PROJECTIONS & COMPONENT OVERLAYS (CORRECTED)
        // ================================================
        TH1D *h1D_m1 = h2D_Mass->ProjectionY(Form("h1D_m1_bin%d", ibin));
        TH1D *h1D_m2 = h2D_Mass->ProjectionX(Form("h1D_m2_bin%d", ibin));

        // Full integrals over the complete fit range [1.0, 1.04] for projection
        double fitMin = 1.0, fitMax = 1.04;
        double Isig_full = fSig.Integral(fitMin, fitMax);
        double Ibkg_full = fBkg.Integral(fitMin, fitMax);

        // Reusable plotting lambda for 1D projections
        auto Draw1DProj = [&](TH1D *h1D, const char *titleX, const char *canvasName)
        {
            TCanvas *c1D = new TCanvas(canvasName, canvasName, 800, 600);
            SetCanvasStyle(c1D, 0.15, 0.05, 0.12, 0.12);
            SetHistoQA(h1D);
            h1D->SetMinimum(0);
            h1D->GetXaxis()->SetTitle(titleX);
            h1D->GetYaxis()->SetTitle("Counts / bin");
            h1D->SetMarkerStyle(20);
            h1D->SetMarkerSize(0.8);
            h1D->Draw("E");

            // Bin width of the projected 1D histogram axis
            double binWidth = h1D->GetXaxis()->GetBinWidth(1);

            // Total 1D Fit (Solid Line)
            TF1 *fTotal1D = new TF1(Form("fTotal1D_%s", canvasName), [=](double *x, double *)
                                    {
                double m = x[0];
                double sig = BWShape(m, m0, gamma);
                double bkg = p5 + p6 * (m - 1.0198) + p7 * (m - 1.0198) * (m - 1.0198);
                
                // Integrate out m2 over [1.0, 1.04] and scale by 1D bin width
                double term_SS = n_SS * sig * Isig_full;
                double term_SB = n_SB * (sig * Ibkg_full + bkg * Isig_full);
                double term_BB = n_BB * bkg * Ibkg_full;
                
                return (term_SS + term_SB + term_BB) * binWidth; }, fitMin, fitMax, 0);

            fTotal1D->SetLineColor(kRed);
            fTotal1D->SetLineWidth(3);
            fTotal1D->SetLineStyle(1); // Solid line

            // N_SS Component (Dotted Line)
            TF1 *fSS1D = new TF1(Form("fSS1D_%s", canvasName), [=](double *x, double *)
                                 {
                double m = x[0];
                double sig = BWShape(m, m0, gamma);
                return n_SS * sig * Isig_full * binWidth; }, fitMin, fitMax, 0);
            fSS1D->SetLineColor(kBlue);
            fSS1D->SetLineWidth(2);
            fSS1D->SetLineStyle(2); // Dotted line

            // N_SB Component (Dotted Line)
            TF1 *fSB1D = new TF1(Form("fSB1D_%s", canvasName), [=](double *x, double *)
                                 {
                double m = x[0];
                double sig = BWShape(m, m0, gamma);
                double bkg = p5 + p6 * (m - 1.0198) + p7 * (m - 1.0198) * (m - 1.0198);
                return n_SB * (sig * Ibkg_full + bkg * Isig_full) * binWidth; }, fitMin, fitMax, 0);
            fSB1D->SetLineColor(kGreen + 2);
            fSB1D->SetLineWidth(2);
            fSB1D->SetLineStyle(2); // Dotted line

            // N_BB Component (Dotted Line)
            TF1 *fBB1D = new TF1(Form("fBB1D_%s", canvasName), [=](double *x, double *)
                                 {
                double m = x[0];
                double bkg = p5 + p6 * (m - 1.0198) + p7 * (m - 1.0198) * (m - 1.0198);
                return n_BB * bkg * Ibkg_full * binWidth; }, fitMin, fitMax, 0);
            fBB1D->SetLineColor(kMagenta);
            fBB1D->SetLineWidth(2);
            fBB1D->SetLineStyle(2); // Dotted line

            // Draw functions on canvas
            fTotal1D->Draw("SAME");
            fSS1D->Draw("SAME");
            fSB1D->Draw("SAME");
            fBB1D->Draw("SAME");

            TLegend *leg = new TLegend(0.65, 0.60, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(h1D, "Data", "pe");
            leg->AddEntry(fTotal1D, "Total Fit", "l");
            leg->AddEntry(fSS1D, "N_{SS}", "l");
            leg->AddEntry(fSB1D, "N_{SB}", "l");
            leg->AddEntry(fBB1D, "N_{BB}", "l");
            leg->Draw();

            c1D->SaveAs(savepath + Form("/1DFits/%s.png", canvasName));
        };

        // Draw 1D distributions for phi1 and phi2
        Draw1DProj(h1D_m1, "#it{M}_{#phi1} (GeV/#it{c}^{2})", Form("1D_Mass_phi1_bin%d", ibin));
        Draw1DProj(h1D_m2, "#it{M}_{#phi2} (GeV/#it{c}^{2})", Form("1D_Mass_phi2_bin%d", ibin));
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
    h_N_SS->Draw("pe");
    cSS->SaveAs(savepath + "/SS_Template.png");

    TCanvas *cNonSS = new TCanvas("cNonSS", "Non-SS Template", 800, 600);
    SetCanvasStyle(cNonSS, 0.15, 0.05, 0.08, 0.12);
    SetHistoQA(h_N_nonSS);
    h_N_nonSS->SetMarkerStyle(20);
    h_N_nonSS->SetMarkerSize(0.8);
    // h_N_nonSS->SetMinimum(h_N_nonSS->GetMaximum() * 0.1);
    // h_N_nonSS->SetMaximum(h_N_nonSS->GetMaximum() * 1.5);
    h_N_nonSS->Draw("pe");
    cNonSS->SaveAs(savepath + "/NonSS_Template.png");

    TFile *fOutput = new TFile(savepath + "/DoublePhiBackgroundTemplates.root", "RECREATE");
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