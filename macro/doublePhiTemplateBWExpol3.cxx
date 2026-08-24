#include "THnSparse.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TF2.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "Math/MinimizerOptions.h"
#include <iostream>
#include <algorithm>

#include "src/style.h"
#include "src/fitfunc.h"

using namespace std;

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
    double n_BS = p[2];
    double n_BB = p[3];

    double m0 = p[4];
    double gamma = p[5];

    double sig1 = BWShape(m1, m0, gamma);
    double sig2 = BWShape(m2, m0, gamma);

    // exp(pol3) background function
    auto Bkg = [&](double m)
    {
        double z = m - mPDG;
        return pow(m, p[6]) * TMath::Exp(p[7] * m + p[8] * m * m + p[9] * m * m * m);
        // return pow(z, p[4]) * TMath::Exp(p[5] * z + p[6] * z * z + p[7] * z * z * z);
    };

    double bkg1 = Bkg(m1);
    double bkg2 = Bkg(m2);

    double shape_SS = sig1 * sig2;
    double shape_SB = sig1 * bkg2;
    double shape_BS = bkg1 * sig2;
    double shape_BB = bkg1 * bkg2;

    return n_SS * shape_SS + n_SB * shape_SB + n_BS * shape_BS + n_BB * shape_BB;
}

void doublePhiTemplateBWExpol3()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template";
    TString suffix = "_ExtendedFitRange";

    // //=========================================================
    // //============Using processOpti5 data======================
    // //=========================================================

    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults_MorePhiBins.root");
    // THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassDoublePhi");

    // int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    // int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    // int lowDeltaM = hUnlike->GetAxis(6)->FindBin(0.0 + 0.00001);
    // int highDeltaM = hUnlike->GetAxis(6)->FindBin(0.005 - 0.00001);

    // hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    // TH3D *h3D_full = hUnlike->Projection(0, 4, 5, "E");

    // hUnlike->GetAxis(6)->SetRange(lowDeltaM, highDeltaM); // Uncommented: essential for DeltaM selection
    // TH3D *h3D_cut = hUnlike->Projection(0, 4, 5, "E");

    // TH1D *hInvMass = hUnlike->Projection(0, "E");
    // TCanvas *cInvMass = new TCanvas("cInvMass", "Invariant Mass Distribution", 720, 720);
    // SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
    // hInvMass->Rebin(12);
    // // hInvMass->Draw("pe");
    // // cInvMass->SaveAs(savepath + "/PhiInvMass.png");

    //=========================================================
    //============Using processOpti8 data======================
    //=========================================================

    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResults_WithPhiMasses.root");
    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassPhiPhiRefitted");
    // Axes: InvMass, pT, deltaM, Chi2, FitProb, Phi1Mass, Phi2Mass

    int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowDeltaM = hUnlike->GetAxis(2)->FindBin(0.0 + 0.00001);
    int highDeltaM = hUnlike->GetAxis(2)->FindBin(0.005 - 0.00001);

    int lowChi2 = hUnlike->GetAxis(3)->FindBin(0.0 + 0.00001);
    int highChi2 = hUnlike->GetAxis(3)->FindBin(25.0 - 0.00001);

    int lowFitProb = hUnlike->GetAxis(4)->FindBin(0.3 + 0.00001);
    int highFitProb = hUnlike->GetAxis(4)->FindBin(2.0 - 0.00001);

    hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    TH3D *h3D_full = hUnlike->Projection(0, 5, 6, "E");
    
    hUnlike->GetAxis(2)->SetRange(lowDeltaM, highDeltaM);
    // hUnlike->GetAxis(3)->SetRange(lowChi2, highChi2);
    hUnlike->GetAxis(4)->SetRange(lowFitProb, highFitProb);
    TH3D *h3D_cut = hUnlike->Projection(0, 5, 6, "E");

    TH1D *hInvMass = hUnlike->Projection(0, "E");
    TCanvas *cInvMass = new TCanvas("cInvMass", "Invariant Mass Distribution", 720, 720);
    SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
    hInvMass->Rebin(10);
    // hInvMass->Draw("pe");
    // cInvMass->SaveAs(savepath + "/PhiInvMass.png");

    // 2. Fix the bin-beating effect by iterating strictly by bin index
    int rebin = 12;
    int startBin = h3D_full->GetXaxis()->FindBin(2.41 + 0.0001);
    int endBin = h3D_full->GetXaxis()->FindBin(2.95 - 0.0001);
    int nBinsInRange = endBin - startBin + 1;
    int totalBins = nBinsInRange / rebin;

    double startBinValue = h3D_full->GetXaxis()->GetBinLowEdge(startBin);
    double endBinValue = h3D_full->GetXaxis()->GetBinUpEdge(endBin);

    cout << "Start bin " << startBin << " with value " << startBinValue << endl;
    cout << "End bin " << endBin << " with value " << endBinValue << endl;
    cout << "Total bins in range: " << nBinsInRange << endl;
    cout << "Total bins after rebinning: " << totalBins << endl;

    // ================================================
    // 1D Histograms for SS and Non-SS Yields
    // ================================================
    TH1D *h_N_SS = new TH1D("h_N_SS", "SS Template (True #phi#phi Yield);#it{M}_{#phi#phi} (GeV/#it{c}^{2});selected N_{SS}", totalBins, startBinValue, endBinValue);
    TH1D *h_N_nonSS = new TH1D("h_N_nonSS", "Non-SS Template (Background Yield);#it{M}_{#phi#phi} (GeV/#it{c}^{2});selected N_{nonSS}", totalBins, startBinValue, endBinValue);
    TH1D *h_N_Total = new TH1D("h_N_Total", "Total Yield;#it{M}_{#phi#phi} (GeV/#it{c}^{2});selected N_{Total}", totalBins, startBinValue, endBinValue);

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
        TF2 *f2D = new TF2(Form("f2D_bin%d", ibin), FitFunc2DBW, 1.0, 1.04, 1.0, 1.04, 10);

        f2D->SetParameter(0, 0.006); // N_SS
        f2D->SetParameter(1, 0.001); // N_BS
        f2D->SetParameter(2, 0.001); // N_SB
        f2D->SetParameter(3, 0.002); // N_BB

        // f2D->SetParLimits(0, 0, 1e9);
        // f2D->SetParLimits(1, 0, 1e9);
        // f2D->SetParLimits(2, 0, 1e9);
        // f2D->SetParLimits(3, 0, 1e9);

        f2D->FixParameter(4, 1.01983);  // Mass peak
        f2D->FixParameter(5, 0.007077); // Width

        // f2D->SetParameter(6,  166.3); // Pol2 p0
        // f2D->SetParameter(7, 68.4); // Pol2 p1
        // f2D->SetParameter(8, 12.4); // Pol2 p2
        // f2D->SetParameter(9, -77.2); // Pol2 p3

        f2D->FixParameter(6, 166.3); // Pol2 p0
        f2D->FixParameter(7, 68.4);  // Pol2 p1
        f2D->FixParameter(8, 12.4);  // Pol2 p2
        f2D->FixParameter(9, -77.2); // Pol2 p3

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
        cout << "N_SS        : " << f2D->GetParameter(0) << " +/- " << f2D->GetParError(0) << endl;
        cout << "N_SB     : " << f2D->GetParameter(1) << " +/- " << f2D->GetParError(1) << endl;
        cout << "N_BS     : " << f2D->GetParameter(2) << " +/- " << f2D->GetParError(2) << endl;
        cout << "N_BB     : " << f2D->GetParameter(3) << " +/- " << f2D->GetParError(3) << endl;
        cout << "Chi2 / NDF  : " << f2D->GetChisquare() / f2D->GetNDF() << endl;
        cout << "Bkg par0   : " << f2D->GetParameter(6) << " +/- " << f2D->GetParError(6) << endl;
        cout << "Bkg par1   : " << f2D->GetParameter(7) << " +/- " << f2D->GetParError(7) << endl;
        cout << "Bkg par2   : " << f2D->GetParameter(8) << " +/- " << f2D->GetParError(8) << endl;
        cout << "Bkg par3   : " << f2D->GetParameter(9) << " +/- " << f2D->GetParError(9) << endl;

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

                double S1 = BWShape(m1, f2D->GetParameter(4), f2D->GetParameter(5));
                double S2 = BWShape(m2, f2D->GetParameter(4), f2D->GetParameter(5));
                double B1 = std::max(0.0, pow(m1, f2D->GetParameter(6)) * TMath::Exp(f2D->GetParameter(7) * m1 + f2D->GetParameter(8) * m1 * m1 + f2D->GetParameter(9) * m1 * m1 * m1));
                double B2 = std::max(0.0, pow(m2, f2D->GetParameter(6)) * TMath::Exp(f2D->GetParameter(7) * m2 + f2D->GetParameter(8) * m2 * m2 + f2D->GetParameter(9) * m2 * m2 * m2));

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
        c2D->SaveAs(savepath + Form("/2DFits/2D_Mass_bin%d.png", ibin));

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
    h_N_SS->GetYaxis()->SetRangeUser(1300, 3750);
    h_N_SS->Draw("pe");
    cSS->SaveAs(savepath + "/SS_Template" + suffix + ".png");

    TCanvas *cNonSS = new TCanvas("cNonSS", "Non-SS Template", 720, 720);
    SetCanvasStyle(cNonSS, 0.15, 0.05, 0.08, 0.12);
    SetHistoQA(h_N_nonSS);
    h_N_nonSS->SetMarkerStyle(20);
    h_N_nonSS->SetMarkerSize(0.8);
    h_N_nonSS->GetYaxis()->SetRangeUser(1300, 2950);
    h_N_nonSS->Draw("pe");
    cNonSS->SaveAs(savepath + "/NonSS_Template" + suffix + ".png");

    TCanvas *cTotal = new TCanvas("cTotal", "Total Yield", 720, 720);
    SetCanvasStyle(cTotal, 0.15, 0.05, 0.08, 0.12);
    SetHistoQA(h_N_Total);
    h_N_Total->SetMarkerStyle(20);
    h_N_Total->SetMarkerSize(0.8);
    h_N_Total->GetYaxis()->SetRangeUser(2800, 5750);
    h_N_Total->Draw("pe");
    hInvMass->SetMarkerStyle(25);
    hInvMass->SetMarkerSize(0.8);
    hInvMass->SetMarkerColor(kRed);
    hInvMass->SetLineColor(kRed);
    hInvMass->Draw("pe SAME");

    TLegend *legClosure = new TLegend(0.5, 0.75, 0.9, 0.9);
    legClosure->SetBorderSize(0);
    legClosure->SetFillStyle(0);
    legClosure->SetTextSize(0.035);
    legClosure->SetTextFont(42);
    legClosure->SetHeader("Closure Test");
    legClosure->AddEntry(h_N_Total, "SS + non-SS", "pe");
    legClosure->AddEntry(hInvMass, "Data", "pe");
    legClosure->Draw();
    // cTotal->SaveAs(savepath + "/Total_Yield" + suffix + ".png");

    TFile *fOutput = new TFile(savepath + Form("/PhiPhiBkgTemplate_BW%s.root", suffix.Data()), "RECREATE");
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