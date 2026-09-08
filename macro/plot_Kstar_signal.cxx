#include <iostream>
#include <cmath>
#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include "src/style.h"

void plot_Kstar_signal()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TString savepath = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/InvMassPlot";
    // 1. Configurable Parameters (Single pT & Multiplicity Bin)
    const std::string dataFilePath = "/home/sawan/Storage/check_k892/data/kstar/LHC22o_pass7/749276.root";
    const std::string templFilePath = "template/buildTemplate/template/ROTATED/SignalMinusTrue.root";
    const std::string folderName = "kstarqa/hInvMass";

    double multLow = 15.0, multHigh = 20.0;
    double ptLow = 2.5, ptHigh = 3.0;
    double normLow = 1.4, normHigh = 1.5;
    double fitLow = 0.75, fitHigh = 1.06;
    int rebinVal = 1;

    double massPDG = 0.89555;
    double widthPDG = 0.0473;

    // 2. Load Data & Reflection Template
    TFile *fInput = TFile::Open(dataFilePath.c_str(), "READ");
    TFile *fTempl = TFile::Open(templFilePath.c_str(), "READ");
    if (!fInput || fInput->IsZombie() || !fTempl || fTempl->IsZombie())
    {
        std::cerr << "Error opening data or template file!" << std::endl;
        return;
    }

    THnSparseF *hSparseNum = (THnSparseF *)fInput->Get(Form("%s/h3KstarInvMassUnlikeSign", folderName.c_str()));
    THnSparseF *hSparseRot = (THnSparseF *)fInput->Get(Form("%s/h3KstarInvMassRotated", folderName.c_str()));

    // Template histogram name (matching path conventions in your code)
    TH1D *hReflRaw = (TH1D *)fTempl->Get(Form("%d-%d/hSigminusTrue_pt_%.1f_%.1f", (int)multLow, (int)multHigh, ptLow, ptHigh));
    if (!hSparseNum || !hSparseRot || !hReflRaw)
    {
        std::cerr << "Histograms or Template not found!" << std::endl;
        cout << Form("%d-%d/hSigminusTrue_pt_%.1f_%.1f", multLow, multHigh, ptLow, ptHigh) << endl;
        return;
    }

    // 3. Projections & Background Subtraction
    hSparseNum->GetAxis(0)->SetRangeUser(multLow, multHigh);
    hSparseNum->GetAxis(1)->SetRangeUser(ptLow, ptHigh);
    hSparseRot->GetAxis(0)->SetRangeUser(multLow, multHigh);
    hSparseRot->GetAxis(1)->SetRangeUser(ptLow, ptHigh);

    TH1D *hSigBkg = hSparseNum->Projection(2, "E");
    TH1D *hRotBkg = hSparseRot->Projection(2, "E");

    hSigBkg->Rebin(rebinVal);
    hRotBkg->Rebin(rebinVal);

    // Normalize Rotational Background
    int bNormLo = hSigBkg->FindBin(normLow + 1e-5);
    int bNormHi = hSigBkg->FindBin(normHigh - 1e-5);
    double normFactor = hSigBkg->Integral(bNormLo, bNormHi) / hRotBkg->Integral(bNormLo, bNormHi);

    TH1D *hBkgScaled = (TH1D *)hRotBkg->Clone("hBkgScaled");
    hBkgScaled->Scale(normFactor);

    // =========================================================================
    // PLOT 1: Signal + Rotational Background & Normalization Region
    // =========================================================================
    TCanvas *c1 = new TCanvas("c1_SigBkg", "Signal + Bkg", 700, 600);
    SetCanvasStyle(c1, 0.14, 0.05, 0.06, 0.13);
    SetHistoQA(hSigBkg);
    SetHistoQA(hBkgScaled);

    double binWidthFit = hSigBkg->GetBinWidth(1);

    hSigBkg->SetTitle(Form(Form(";M_{K^{#pm}#pi^{#mp}} (GeV/#it{c}^{2});Counts/(%.3f GeV/#it{c}^{2})", binWidthFit)));
    hSigBkg->SetMarkerStyle(20);
    hSigBkg->SetMarkerSize(1.0);
    hSigBkg->SetMaximum(hSigBkg->GetMaximum() * 1.18);
    hSigBkg->SetMinimum(0.13e6);
    hSigBkg->Draw("E");

    hBkgScaled->SetMarkerStyle(24);
    hBkgScaled->SetMarkerSize(1.0);
    hBkgScaled->SetMarkerColor(kRed);
    hBkgScaled->SetLineColor(kRed);
    hBkgScaled->Draw("E SAME");

    TH1D *hNormRegion = (TH1D *)hBkgScaled->Clone("hNormRegion");
    hNormRegion->SetFillColorAlpha(kRed, 0.35);
    hNormRegion->SetFillStyle(3001);
    for (int ibin = 1; ibin <= hNormRegion->GetNbinsX(); ibin++)
    {
        double bc = hNormRegion->GetBinCenter(ibin);
        if (bc < normLow || bc > normHigh)
            hNormRegion->SetBinContent(ibin, -999);
    }
    hNormRegion->SetLineWidth(0);
    hNormRegion->Draw("BAR SAME");

    TLegend *leg1 = new TLegend(0.36, 0.20, 0.70, 0.40);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.04);
    leg1->AddEntry(hSigBkg, "Unlike-sign pairs", "p");
    leg1->AddEntry(hBkgScaled, "Rotational background", "p");
    leg1->AddEntry(hNormRegion, "Normalization region", "f");
    leg1->Draw();

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextFont(42);
    ltx.SetTextSize(0.04);
    // ltx.DrawLatex(0.18, 0.82, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));
    ltx.DrawLatex(0.17, 0.85, "ALICE");
    ltx.DrawLatex(0.17, 0.79, "pp #sqrt{#it{s}} = 13.6 TeV");
    ltx.DrawLatex(0.17, 0.73, Form("|#it{y}| < 0.5"));

    ltx.DrawLatex(0.59, 0.85, "K*(892)^{0} #rightarrow K^{#pm} + #pi^{#mp}");
    ltx.DrawLatex(0.59, 0.79, Form("FT0M class: (%d - %d)%%", (int)multLow, (int)multHigh));
    ltx.DrawLatex(0.59, 0.73, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));

    c1->SaveAs(savepath + "/SigPlusBkg.png");   
    // =========================================================================
    // PLOT 2: Background-Subtracted Signal Fit (Voigtian + Template + Pol3)
    // =========================================================================
    TH1D *hSubSignal = (TH1D *)hSigBkg->Clone("hSubSignal");
    hSubSignal->Add(hBkgScaled, -1.0);

    TH1D *hReflection = (TH1D *)hReflRaw->Clone("hReflection");
    hReflection->Rebin(rebinVal);
    for (int ib = 1; ib <= hReflection->GetNbinsX(); ib++)
    {
        if (hReflection->GetBinContent(ib) < 0)
        {
            hReflection->SetBinContent(ib, 0.0);
            hReflection->SetBinError(ib, 0.0);
        }
    }

    const double reflBinWidth = hReflection->GetBinWidth(1);
    int bReflLo = hReflection->GetXaxis()->FindBin(fitLow + 1e-6);
    int bReflHi = hReflection->GetXaxis()->FindBin(fitHigh - 1e-6);
    const double reflNorm = std::max(hReflection->Integral(bReflLo, bReflHi) * reflBinWidth, 1e-12);
    const double totalInFit = std::max(1.0, hSubSignal->Integral(hSubSignal->FindBin(fitLow), hSubSignal->FindBin(fitHigh)));

    // Chebyshev Pol3 Residual Bkg Density Evaluator Lambda
    auto evalBkgDensity = [fitLow, fitHigh](double x, const double *par) -> double
    {
        double t = 2.0 * (x - fitLow) / (fitHigh - fitLow) - 1.0;
        double T1 = t, T2 = 2.0 * t * t - 1.0, T3 = 4.0 * t * t * t - 3.0 * t;
        double chebRaw = 1.0 + par[0] * T1 + par[1] * T2 + par[2] * T3;
        double chebInt = ((fitHigh - fitLow) / 2.0) * (2.0 - (2.0 / 3.0) * par[1]);
        return chebRaw / std::max(chebInt, 1e-12);
    };

    // Total Composite Fit Function
    TF1 *fTotal = new TF1("fTotal", [hReflection, reflNorm, fitLow, fitHigh, binWidthFit, evalBkgDensity](double *xx, double *pp) -> double
                          {
        double x = xx[0];
        double Nsig = pp[0], mass = pp[1], width = pp[2], sigma = pp[3];
        double Ncorr = pp[4], Nbkg = pp[5];

        // 1. Normalized Voigtian Signal Density
        double voigtRaw = TMath::Voigt(x - mass, sigma, width);
        int nSteps = 40;
        double h_step = (fitHigh - fitLow) / nSteps;
        double vInt = TMath::Voigt(fitLow - mass, sigma, width) + TMath::Voigt(fitHigh - mass, sigma, width);
        for (int j = 1; j < nSteps; j++) {
            double xj = fitLow + j * h_step;
            vInt += TMath::Voigt(xj - mass, sigma, width) * ((j % 2 == 0) ? 2.0 : 4.0);
        }
        vInt *= h_step / 3.0;
        double sigDensity = voigtRaw / std::max(vInt, 1e-12);

        // 2. MC Correlated Template Density
        double corrDensity = (hReflection && hReflection->GetEntries() > 0) ? hReflection->Interpolate(x) / reflNorm : 0.0;
        if (corrDensity < 0) corrDensity = 0.0;

        // 3. Residual Background Density
        double bkgDensity = evalBkgDensity(x, &pp[6]);

        return binWidthFit * (Nsig * sigDensity + Ncorr * corrDensity + Nbkg * bkgDensity); }, fitLow, fitHigh, 9);

    fTotal->SetParNames("Nsig", "mass", "width", "sigma", "Ncorr", "Nbkg", "c1", "c2", "c3");
    fTotal->SetParameters(totalInFit * 0.4, massPDG, widthPDG, 0.005, totalInFit * 0.2, totalInFit * 0.3, -0.5, 0.1, -0.05);
    fTotal->FixParameter(2, widthPDG); // Fix width to PDG as in standard analysis
    fTotal->SetParLimits(1, massPDG - 0.01, massPDG + 0.01);
    fTotal->SetParLimits(3, 0.0005, 0.02);

    TCanvas *c2 = new TCanvas("c2_SubtractedFit", "Template Fit", 700, 600);
    SetCanvasStyle(c2, 0.14, 0.05, 0.06, 0.13);
    SetHistoQA(hSubSignal);
    hSubSignal->Fit(fTotal, "R Q S");

    hSubSignal->SetTitle(Form(";M_{K^{#pm}#pi^{#mp}} (GeV/#it{c}^{2});Counts/(%.3f GeV/#it{c}^{2})", binWidthFit));
    hSubSignal->GetXaxis()->SetRangeUser(0.75, 1.05);
    hSubSignal->SetMarkerStyle(20);
    hSubSignal->SetMarkerSize(1.0);
    hSubSignal->SetMinimum(-0.08 * hSubSignal->GetMinimum());
    hSubSignal->SetMaximum(hSubSignal->GetMaximum() * 1.3);
    hSubSignal->Draw("E");

    // Extract Components for Visualization
    TF1 *fSig = new TF1("fSig", [fitLow, fitHigh, binWidthFit](double *xx, double *pp) -> double
                        {
        double x = xx[0];
        double voigtRaw = TMath::Voigt(x - pp[1], pp[3], pp[2]);
        int nSteps = 40; double h_step = (fitHigh - fitLow) / nSteps;
        double vInt = TMath::Voigt(fitLow - pp[1], pp[3], pp[2]) + TMath::Voigt(fitHigh - pp[1], pp[3], pp[2]);
        for (int j = 1; j < nSteps; j++) vInt += TMath::Voigt(fitLow + j * h_step - pp[1], pp[3], pp[2]) * ((j % 2 == 0) ? 2.0 : 4.0);
        return binWidthFit * pp[0] * voigtRaw / std::max(vInt * h_step / 3.0, 1e-12); }, fitLow, fitHigh, 4);
    fSig->SetParameters(fTotal->GetParameter(0), fTotal->GetParameter(1), fTotal->GetParameter(2), fTotal->GetParameter(3));
    fSig->SetLineColor(kMagenta + 1);
    // fSig->SetNpx(5000);
    fSig->Draw("SAME");

    TF1 *fCorr = new TF1("fCorr", [hReflection, reflNorm, binWidthFit](double *xx, double *pp) -> double
                         { return binWidthFit * pp[0] * (hReflection->Interpolate(xx[0]) / reflNorm); }, fitLow, fitHigh, 1);
    fCorr->SetParameter(0, fTotal->GetParameter(4));
    fCorr->SetLineColor(kGreen + 2);
    fCorr->SetNpx(5000);
    fCorr->Draw("SAME");

    TF1 *fBkg = new TF1("fBkg", [fitLow, fitHigh, binWidthFit, evalBkgDensity](double *xx, double *pp) -> double
                        { return binWidthFit * pp[0] * evalBkgDensity(xx[0], &pp[1]); }, fitLow, fitHigh, 4);
    fBkg->SetParameter(0, fTotal->GetParameter(5));
    fBkg->SetParameters(fTotal->GetParameter(5), fTotal->GetParameter(6), fTotal->GetParameter(7), fTotal->GetParameter(8));
    fBkg->SetLineColor(kBlue);
    // fBkg->SetNpx(5000);
    fBkg->Draw("SAME");

    fTotal->SetLineColor(kRed);
    fTotal->SetLineWidth(2);
    fTotal->SetNpx(5000);
    fTotal->Draw("SAME");

    // Combined Total Background (MC Reflection + Residual Background)
    TF1 *fTotalBkg = new TF1("fTotalBkg", [fCorr, fBkg](double *xx, double *pp) -> double
                             { return fCorr->Eval(xx[0]) + fBkg->Eval(xx[0]); }, fitLow, fitHigh, 0);

    fTotalBkg->SetLineColor(kCyan + 2);
    fTotalBkg->SetLineStyle(2); // Dashed line to distinguish from Total Fit
    fTotalBkg->SetLineWidth(2);
    fTotalBkg->SetNpx(5000);
    fTotalBkg->Draw("SAME");

    TLegend *leg2 = new TLegend(0.60, 0.63, 0.95, 0.92);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->AddEntry(hSubSignal, "Background subtracted", "p");
    leg2->AddEntry(fTotal, "Signal+Background fit", "l");
    leg2->AddEntry(fSig, "Voigtian (Signal)", "l");
    leg2->AddEntry(fCorr, "MC Reflection Template", "l");
    leg2->AddEntry(fBkg, "Chebyshev T_{3}", "l");
    leg2->AddEntry(fTotalBkg, "Template + Chebyshev", "l");
    leg2->Draw();
    ltx.SetTextSize(0.03);
    ltx.DrawLatex(0.690544, 0.607639, "(Background)");
    c2->SaveAs(savepath + "/SubtractedFit.png");
}