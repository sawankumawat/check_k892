#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "src/style.h"

using namespace std;

// Function declarations
TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const std::string &name);

Double_t expPol3(Double_t *x, Double_t *par)
{
    return (pow(x[0], par[0])) * TMath::Exp(
                                     par[1] * x[0] +
                                     par[2] * x[0] * x[0] +
                                     par[3] * x[0] * x[0] * x[0]);
}

Double_t breitWigner(Double_t *x, Double_t *par)
{
    double m = x[0];
    double amp = par[0];
    double mass = par[1];
    double width = par[2];

    double denominator = (m - mass) * (m - mass) + width * width / 4.0;
    return amp * width / (TMath::Pi() * 2 * denominator);
}

Double_t BWExpol(Double_t *x, Double_t *par)
{
    return breitWigner(x, par) + expPol3(x, &par[3]);
}

// Global fixed mass value from merged dataset
const double FIXED_MASS = 2.6902;

// Combined Multi-Histogram Global Objective Structure for Stage 3
struct GlobalChi2
{
    GlobalChi2(ROOT::Fit::Chi2FCN<ROOT::Math::IMultiGenFunction> &f1,
               ROOT::Fit::Chi2FCN<ROOT::Math::IMultiGenFunction> &f2,
               ROOT::Fit::Chi2FCN<ROOT::Math::IMultiGenFunction> &f3)
        : fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

    // Global parameters mapping (19 total parameters): (Shared mass and width free)
    // [0]: Shared Mass
    // [1..6]   Period 1: Yield_1, Width_1, p0_1, p1_1, p2_1, p3_1
    // [7..12]  Period 2: Yield_2, Width_2, p0_2, p1_2, p2_2, p3_2
    // [13..18] Period 3: Yield_3, Width_3, p0_3, p1_3, p2_3, p3_3
    double operator()(const double *par) const
    {
        // Function signature: BWExpol(x, par) where par = [Yield, Mass, Width, p0, p1, p2, p3]
        double p1[7] = {par[1], par[0], par[2], par[3], par[4], par[5], par[6]};
        double p2[7] = {par[7], par[0], par[8], par[9], par[10], par[11], par[12]};
        double p3[7] = {par[13], par[0], par[14], par[15], par[16], par[17], par[18]};

        return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
    }

    ROOT::Fit::Chi2FCN<ROOT::Math::IMultiGenFunction> *fChi2_1;
    ROOT::Fit::Chi2FCN<ROOT::Math::IMultiGenFunction> *fChi2_2;
    ROOT::Fit::Chi2FCN<ROOT::Math::IMultiGenFunction> *fChi2_3;
};

void doublePhiSimultaneousFit_SharedMass()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    vector<string> periods = {"26afai", "26ac", "26adaeag"};
    vector<string> periodLabels = {"LHC26(af+ai)", "LHC26ac", "LHC26(ad+ae+ag)"};
    int nPeriods = periods.size();

    int rebinFactor = 8;
    float fitRangeLow = 2.41;
    float fitRangeHigh = 2.90;

    vector<TH1D *> hInvMassVec(nPeriods);
    vector<TH1D *> hBkgVec(nPeriods);
    vector<TF1 *> fitBkgVec(nPeriods);
    vector<TF1 *> fitIndivSBVec(nPeriods);

    // Arrays to store background parameters and individual fit results per period
    double bkgPars[3][4];
    double indivYields[3];
    double indivWidths[3];

    cout << "\n=======================================================" << endl;
    cout << " STAGE 1 & STAGE 2: INDIVIDUAL PERIOD PRE-FITS" << endl;
    cout << "=======================================================" << endl;

    for (int i = 0; i < nPeriods; ++i)
    {
        string filePath = Form("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/All/AnalysisResults_LHC%s.root", periods[i].c_str());
        TFile *fInput = OpenFile(filePath);
        if (!fInput)
            return;

        THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassPhiPhiRefitted");

        int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
        int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);
        int lowDeltaM = hUnlike->GetAxis(2)->FindBin(0.0 + 0.00001);
        int highDeltaM = hUnlike->GetAxis(2)->FindBin(0.005 - 0.00001);
        int lowChi2 = hUnlike->GetAxis(3)->FindBin(0.0 + 0.00001);
        int highChi2 = hUnlike->GetAxis(3)->FindBin(25.0 - 0.00001);
        int lowFitProb = hUnlike->GetAxis(4)->FindBin(0.3 + 0.00001);
        int highFitProb = hUnlike->GetAxis(4)->FindBin(2.0 - 0.00001);

        hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
        hUnlike->GetAxis(2)->SetRange(lowDeltaM, highDeltaM);
        hUnlike->GetAxis(3)->SetRange(lowChi2, highChi2);
        hUnlike->GetAxis(4)->SetRange(lowFitProb, highFitProb);

        TH1D *hInvMass = hUnlike->Projection(0, "E");
        hInvMass->SetName(Form("hInvMass_%s", periods[i].c_str()));
        hInvMass->Rebin(rebinFactor);
        hInvMass->GetXaxis()->SetRangeUser(fitRangeLow, fitRangeHigh);
        hInvMassVec[i] = hInvMass;

        // -----------------------------------------------------------------
        // Stage 1: Fit Background-Only with Signal Region Excluded
        // -----------------------------------------------------------------
        TH1D *hBkg = (TH1D *)hInvMass->Clone(Form("hBkg_%s", periods[i].c_str()));
        hBkg->Reset();

        for (int b = 1; b <= hInvMass->GetNbinsX(); b++)
        {
            double x = hInvMass->GetBinCenter(b);
            if (x >= 2.62 && x <= 2.73)
                continue; // Exclude signal region
            hBkg->SetBinContent(b, hInvMass->GetBinContent(b));
            hBkg->SetBinError(b, hInvMass->GetBinError(b));
        }
        hBkgVec[i] = hBkg;

        TF1 *fitBkg = new TF1(Form("fitBkg_%s", periods[i].c_str()), expPol3, fitRangeLow, fitRangeHigh, 4);
        fitBkg->SetParNames("p0", "p1", "p2", "p3");
        fitBkg->SetParameters(-2.7e2, 64.0, 40.0, -6.5);

        TFitResultPtr resBkg = hBkg->Fit(fitBkg, "REBMSQ0");
        fitBkgVec[i] = fitBkg;

        for (int p = 0; p < 4; ++p)
        {
            bkgPars[i][p] = fitBkg->GetParameter(p);
        }

        cout << Form("\n--- Stage 1: Bkg-Only Fit (LHC%s) ---", periods[i].c_str()) << endl;
        cout << "Fit Status       : " << resBkg->Status() << " (" << (resBkg->IsValid() ? "Valid" : "Failed") << ")" << endl;
        cout << "Cov Matrix Status: " << resBkg->CovMatrixStatus() << endl;
        cout << Form("Chi2/NDF         : %.2f / %d = %.3f", fitBkg->GetChisquare(), fitBkg->GetNDF(), fitBkg->GetChisquare() / fitBkg->GetNDF()) << endl;

        // -----------------------------------------------------------------
        // Stage 2: Fit S+B Individually per Period
        // Fixed Mass = 2.6902, Fixed Bkg from Stage 1, Floating Width & Yield
        // -----------------------------------------------------------------
        TF1 *fitIndivSB = new TF1(Form("fitIndivSB_%s", periods[i].c_str()), BWExpol, fitRangeLow, fitRangeHigh, 7);
        fitIndivSB->SetParNames("SignalYield", "Mass", "Width", "p0", "p1", "p2", "p3");

        fitIndivSB->SetParameter(0, 100.0);
        fitIndivSB->FixParameter(1, FIXED_MASS);
        fitIndivSB->SetParameter(2, 0.03);
        fitIndivSB->SetParLimits(2, 0.005, 0.15);

        for (int p = 0; p < 4; ++p)
        {
            fitIndivSB->FixParameter(3 + p, bkgPars[i][p]);
        }

        TFitResultPtr resIndivSB = hInvMass->Fit(fitIndivSB, "RLSQ0S");
        fitIndivSBVec[i] = fitIndivSB;

        indivYields[i] = fitIndivSB->GetParameter(0);
        indivWidths[i] = fitIndivSB->GetParameter(2);

        cout << Form("--- Stage 2: Individual S+B Fit (LHC%s) ---", periods[i].c_str()) << endl;
        cout << "Fit Status       : " << resIndivSB->Status() << " (" << (resIndivSB->IsValid() ? "Valid" : "Failed") << ")" << endl;
        cout << "Cov Matrix Status: " << resIndivSB->CovMatrixStatus() << endl;
        cout << Form("Chi2/NDF         : %.2f / %d = %.3f", fitIndivSB->GetChisquare(), fitIndivSB->GetNDF(), fitIndivSB->GetChisquare() / fitIndivSB->GetNDF()) << endl;
        cout << Form("Results          : Yield = %.2f +/- %.2f | Width = %.4f +/- %.4f GeV/c2",
                     indivYields[i], fitIndivSB->GetParError(0),
                     indivWidths[i], fitIndivSB->GetParError(2))
             << endl;
    }

    // -----------------------------------------------------------------
    // Stage 3: Simultaneous Global Fit Across All Periods
    // Shared Mass, Free Floating Widths, Independent Yields & Backgrounds
    // -----------------------------------------------------------------
    cout << "\n=======================================================" << endl;
    cout << " STAGE 3: SIMULTANEOUS GLOBAL FIT ACROSS ALL PERIODS" << endl;
    cout << "=======================================================" << endl;

    vector<TF1 *> fitSimulVec(nPeriods);
    for (int i = 0; i < nPeriods; ++i)
    {
        fitSimulVec[i] = new TF1(Form("fitSimul_%s", periods[i].c_str()), BWExpol, fitRangeLow, fitRangeHigh, 7);
    }

    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange range;
    range.SetRange(fitRangeLow, fitRangeHigh);

    ROOT::Fit::BinData data1(opt, range);
    ROOT::Fit::BinData data2(opt, range);
    ROOT::Fit::BinData data3(opt, range);

    ROOT::Fit::FillData(data1, hInvMassVec[0]);
    ROOT::Fit::FillData(data2, hInvMassVec[1]);
    ROOT::Fit::FillData(data3, hInvMassVec[2]);

    ROOT::Math::WrappedMultiTF1 fitFunction1(*fitSimulVec[0], 1);
    ROOT::Math::WrappedMultiTF1 fitFunction2(*fitSimulVec[1], 1);
    ROOT::Math::WrappedMultiTF1 fitFunction3(*fitSimulVec[2], 1);

    ROOT::Fit::Chi2FCN<ROOT::Math::IMultiGenFunction> chi2_1(data1, fitFunction1);
    ROOT::Fit::Chi2FCN<ROOT::Math::IMultiGenFunction> chi2_2(data2, fitFunction2);
    ROOT::Fit::Chi2FCN<ROOT::Math::IMultiGenFunction> chi2_3(data3, fitFunction3);

    GlobalChi2 globalChi2(chi2_1, chi2_2, chi2_3);

    ROOT::Fit::Fitter fitter;
    // Expanded to 19 total parameters
    fitter.Config().SetParamsSettings(19, std::vector<double>(19, 0.0).data());

    // Parameter [0]: Shared floating Mass across all periods
    fitter.Config().ParSettings(0).SetName("SharedMass");
    fitter.Config().ParSettings(0).SetValue(FIXED_MASS); // Initialized to 2.6902
    fitter.Config().ParSettings(0).SetLimits(2.65, 2.73);

    // Seed per-period parameters using Stage 1 & Stage 2 results
    for (int i = 0; i < 3; ++i)
    {
        int baseIndex = 1 + i * 6;

        // Yield initialized from Stage 2
        fitter.Config().ParSettings(baseIndex + 0).SetName(Form("Yield_%s", periods[i].c_str()));
        fitter.Config().ParSettings(baseIndex + 0).SetValue(indivYields[i]);
        fitter.Config().ParSettings(baseIndex + 0).SetLimits(0.0, 1.0e6);

        // Independent Width initialized from Stage 2
        fitter.Config().ParSettings(baseIndex + 1).SetName(Form("Width_%s", periods[i].c_str()));
        fitter.Config().ParSettings(baseIndex + 1).SetValue(indivWidths[i]);
        fitter.Config().ParSettings(baseIndex + 1).SetLimits(0.005, 0.15);

        // Background shape initialized from Stage 1
        for (int p = 1; p <= 4; ++p)
        {
            fitter.Config().ParSettings(baseIndex + 2 + (p - 1)).SetName(Form("bkg_p%d_%s", p - 1, periods[i].c_str()));
            fitter.Config().ParSettings(baseIndex + 2 + (p - 1)).SetValue(bkgPars[i][p - 1]);
        }
    }

    int totalNPoints = data1.Size() + data2.Size() + data3.Size();
    fitter.FitFCN(19, globalChi2, 0, totalNPoints, true);
    ROOT::Fit::FitResult result = fitter.Result();

    double sharedMass = result.Value(0);
    double sharedMassErr = result.Error(0);
    double globalChi2Val = result.Chi2();
    int globalNDF = result.Ndf();

    cout << "\n--- Stage 3: Simultaneous Global Fit Summary ---" << endl;
    cout << "Fit Status       : " << result.Status() << " (" << (result.IsValid() ? "Valid" : "Failed") << ")" << endl;
    cout << "Cov Matrix Status: " << result.CovMatrixStatus() << endl;
    cout << Form("Global Chi2/NDF  : %.2f / %d = %.3f", globalChi2Val, globalNDF, globalChi2Val / globalNDF) << endl;
    cout << Form("Shared Mass      : %.4f +/- %.4f GeV/c2", sharedMass, sharedMassErr) << endl;

    for (int i = 0; i < nPeriods; ++i)
    {
        int baseIndex = 1 + i * 6;
        cout << Form("  -> Period LHC%s: Yield = %.1f +/- %.1f | Free Width = %.4f +/- %.4f GeV/c2",
                     periods[i].c_str(),
                     result.Value(baseIndex + 0), result.Error(baseIndex + 0),
                     result.Value(baseIndex + 1), result.Error(baseIndex + 1))
             << endl;
    }
    cout << "=======================================================\n"
         << endl;

    // Plotting per-period spectra with simultaneous fit results
    TCanvas *cMulti = new TCanvas("cMulti", "Period-wise Mass Spectra", 1350, 450);
    cMulti->Divide(3, 1);

    for (int i = 0; i < nPeriods; ++i)
    {
        cMulti->cd(i + 1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.02);
        gPad->SetTopMargin(0.06);
        gPad->SetBottomMargin(0.14);
        int baseIndex = 1 + i * 6;
        double simulYield = result.Value(baseIndex + 0);
        double simulYieldErr = result.Error(baseIndex + 0);
        double indivWidth = result.Value(baseIndex + 1);
        double indivWidthErr = result.Error(baseIndex + 1);

        fitSimulVec[i]->SetParameters(simulYield, sharedMass, indivWidth,
                                      result.Value(baseIndex + 2),
                                      result.Value(baseIndex + 3),
                                      result.Value(baseIndex + 4),
                                      result.Value(baseIndex + 5));

        fitSimulVec[i]->SetLineColor(kRed + 1);
        fitSimulVec[i]->SetLineWidth(2);
        SetHistoQA(hInvMassVec[i]);
        hInvMassVec[i]->SetMarkerSize(0.5);
        hInvMassVec[i]->SetLineWidth(1);
        hInvMassVec[i]->GetYaxis()->SetTitle(Form("Counts / %.2f MeV/#it{c}^{2}", hInvMassVec[i]->GetBinWidth(1) * 1000.0));
        hInvMassVec[i]->GetXaxis()->SetTitle("M_{#phi#phi} (GeV/#it{c}^{2})");
        hInvMassVec[i]->Draw("PE");
        fitSimulVec[i]->Draw("SAME");

        TF1 *fitBkgDraw = new TF1(Form("fitBkgDraw_%s", periods[i].c_str()), expPol3, fitRangeLow, fitRangeHigh, 4);
        fitBkgDraw->SetParameters(result.Value(baseIndex + 2), result.Value(baseIndex + 3),
                                  result.Value(baseIndex + 4), result.Value(baseIndex + 5));
        fitBkgDraw->SetLineColor(kBlue);
        fitBkgDraw->SetLineStyle(2);
        fitBkgDraw->Draw("SAME");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.33, 0.89, Form("Period: %s", periodLabels[i].c_str()));
        latex.SetTextSize(0.035);
        latex.SetTextFont(42);
        latex.DrawLatex(0.33, 0.82, Form("#Gamma = %.4f #pm %.4f GeV/c^{2}", indivWidth, indivWidthErr));
        if (i == 0)
        {
            latex.DrawLatex(0.33, 0.77, Form("Shared M = %.4f #pm %.4f GeV/c^{2}", sharedMass, sharedMassErr));
            latex.DrawLatex(0.33, 0.72, Form("#Chi^{2}/NDF  : %.2f / %d = %.3f", globalChi2Val, globalNDF, globalChi2Val / globalNDF));
        }
        // latex.DrawLatex(0.33, 0.67, Form("Yield = %.1f #pm %.1f", simulYield, simulYieldErr));
    }

    cMulti->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/TetraquarkFit/PeriodWise26Data/LHC26_PeriodStability_SharedMass.png");
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

template <typename T>
T *GetHisto(TFile *f, const std::string &name)
{
    T *histo = dynamic_cast<T *>(f->Get(name.c_str()));
    if (!histo)
    {
        std::cout << "Error: histo " << name << " not found in file " << f->GetName() << std::endl;
        return nullptr;
    }
    return histo;
}