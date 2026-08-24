#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "src/style.h"

TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const string &name);

void fitPhiPhiTemplateVoigt()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    // Paths
    TString inFilePath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template/PhiPhiBkgTemplate_BW_ExtendedFitRange.root";
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template";
    TString suffix = "VoigtExpol3_extendedFitRange";
    // TString suffix = "Voigt";
    // TString suffix = "pol2";

    double fitLow = 2.41;
    double fitHigh = 2.95;

    // Open input ROOT file
    TFile *fInput = OpenFile(inFilePath.Data());

    // Retrieve the templates
    TH1D *h_N_SS = GetHisto<TH1D>(fInput, "h_N_SS");
    TH1D *h_N_nonSS = GetHisto<TH1D>(fInput, "h_N_nonSS");

    // Exclusion region bounds for the signal window
    double excLow = 2.650;
    double excHigh = 2.730;

    // =========================================================================
    // STEP 1: Fit Sideband Regions with Chebyshev-3 Polynomials
    // =========================================================================

    // Chebyshev-3 sideband fit function with window rejection
    auto chebyshev3_sideband = [excLow, excHigh, fitLow, fitHigh](double *x, double *p) -> double
    {
        double m = x[0];
        if (m > excLow && m < excHigh)
        {
            TF1::RejectPoint();
            return 0.0;
        }
        // Normalize x from [fitLow, fitHigh] to [-1, 1]
        double x_norm = (2.0 * (m - fitLow) / (fitHigh - fitLow)) - 1.0;
        double T0 = 1.0;
        double T1 = x_norm;
        double T2 = 2.0 * x_norm * x_norm - 1.0;
        double T3 = 4.0 * x_norm * x_norm * x_norm - 3.0 * x_norm;

        return p[0] * T0 + p[1] * T1 + p[2] * T2 + p[3] * T3;
    };

    // Full range evaluation for interpolation through signal region
    auto chebyshev3_full = [fitLow, fitHigh](double *x, double *p) -> double
    {
        double m = x[0];
        double x_norm = (2.0 * (m - fitLow) / (fitHigh - fitLow)) - 1.0;
        double T0 = 1.0;
        double T1 = x_norm;
        double T2 = 2.0 * x_norm * x_norm - 1.0;
        double T3 = 4.0 * x_norm * x_norm * x_norm - 3.0 * x_norm;

        return p[0] * T0 + p[1] * T1 + p[2] * T2 + p[3] * T3;
    };

    TF1 *f_SS_sideband = new TF1("f_SS_sideband", chebyshev3_sideband, fitLow, fitHigh, 4);
    TF1 *f_nonSS_sideband = new TF1("f_nonSS_sideband", chebyshev3_sideband, fitLow, fitHigh, 4);

    f_SS_sideband->SetParameters(1500, -300, 50, -10);
    f_nonSS_sideband->SetParameters(1500, -50, 10, -5);

    h_N_SS->Fit(f_SS_sideband, "R0Q");
    h_N_nonSS->Fit(f_nonSS_sideband, "R0Q");

    // Extrapolated templates across full mass range
    TF1 *f_SS_template = new TF1("f_SS_template", chebyshev3_full, fitLow, fitHigh, 4);
    TF1 *f_nonSS_template = new TF1("f_nonSS_template", chebyshev3_full, fitLow, fitHigh, 4);

    f_SS_template->SetParameters(f_SS_sideband->GetParameters());
    f_nonSS_template->SetParameters(f_nonSS_sideband->GetParameters());

    // Plot Step 1 Results
    TCanvas *cTemplates = new TCanvas("cTemplates", "Background Templates", 1280, 620);
    SetCanvasStyle(cTemplates, 0.15, 0.03, 0.05, 0.15);
    cTemplates->Divide(2, 1);

    // Left sideband SS
    TF1 *f_SS_left = new TF1("f_SS_left", chebyshev3_full, fitLow, excLow, 4);
    f_SS_left->SetParameters(f_SS_sideband->GetParameters());
    f_SS_left->SetLineColor(kRed);
    f_SS_left->SetLineWidth(3);

    // Right sideband SS
    TF1 *f_SS_right = new TF1("f_SS_right", chebyshev3_full, excHigh, fitHigh, 4);
    f_SS_right->SetParameters(f_SS_sideband->GetParameters());
    f_SS_right->SetLineColor(kRed);
    f_SS_right->SetLineWidth(3);

    // Left sideband non-SS
    TF1 *f_nonSS_left = new TF1("f_nonSS_left", chebyshev3_full, fitLow, excLow, 4);
    f_nonSS_left->SetParameters(f_nonSS_sideband->GetParameters());
    f_nonSS_left->SetLineColor(kRed);
    f_nonSS_left->SetLineWidth(3);

    // Right sideband non-SS
    TF1 *f_nonSS_right = new TF1("f_nonSS_right", chebyshev3_full, excHigh, fitHigh, 4);
    f_nonSS_right->SetParameters(f_nonSS_sideband->GetParameters());
    f_nonSS_right->SetLineColor(kRed);
    f_nonSS_right->SetLineWidth(3);

    cTemplates->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.10);

    h_N_SS->SetTitle("SS template (Chebyshev-3);M_{#phi#phi} (GeV/c^{2});selected N_{SS}");
    h_N_SS->SetMarkerStyle(20);
    SetHistoQA(h_N_SS);
    h_N_SS->SetTitle("SS template (Chebyshev-3)");
    h_N_SS->SetMaximum(h_N_SS->GetMaximum() * 1.5);
    h_N_SS->SetMinimum(h_N_SS->GetMinimum() * 0.7);
    h_N_SS->Draw("PE");
    f_SS_sideband->SetLineColor(kRed);
    f_SS_sideband->SetLineWidth(3);
    // f_SS_sideband->Draw("SAME");
    f_SS_left->Draw("SAME");
    f_SS_right->Draw("SAME");
    f_SS_template->SetLineStyle(2); // Dotted line interpolation
    f_SS_template->SetLineColor(kRed);
    f_SS_template->Draw("SAME");

    TLegend *leg = new TLegend(0.35, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.033);
    leg->AddEntry(h_N_SS, "Extracted N_{SS}", "p");
    leg->AddEntry(f_SS_sideband, "Chebyshev-3 fit (side-bands)", "l");
    leg->AddEntry(f_SS_template, "Interpolated in excluded window", "l");
    leg->Draw();

    cTemplates->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.10);

    h_N_nonSS->SetTitle("Non-SS template (Chebyshev-3);M_{#phi#phi} (GeV/c^{2});selected N_{non-SS}");
    h_N_nonSS->SetMarkerStyle(20);
    SetHistoQA(h_N_nonSS);
    h_N_nonSS->SetTitle("Non-SS template (Chebyshev-3)");
    h_N_nonSS->Draw("PE");
    f_nonSS_sideband->SetLineColor(kRed);
    f_nonSS_sideband->SetLineWidth(3);
    // f_nonSS_sideband->Draw("SAME");
    f_nonSS_left->Draw("SAME");
    f_nonSS_right->Draw("SAME");
    f_nonSS_template->SetLineStyle(2);
    f_nonSS_template->SetLineColor(kRed);
    f_nonSS_template->Draw("SAME");

    TLegend *leg3 = new TLegend(0.5, 0.75, 0.88, 0.88);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.033);
    leg3->AddEntry(h_N_nonSS, "Extracted N_{non-SS}", "p");
    leg3->Draw();

    cTemplates->SaveAs(savepath + "/Bkg_Templates=" + suffix + ".png");

    // =========================================================================
    // STEP 2: Combine Templates to Form Total Background Fit
    // =========================================================================

    // TFile *fInputData = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults_MorePhiBins.root");
    // THnSparseF *hUnlike = GetHisto<THnSparseF>(fInputData, "doublephimeson/SEMassDoublePhi");

    // int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    // int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    // int lowDeltaM = hUnlike->GetAxis(6)->FindBin(0.0 + 0.00001);
    // int highDeltaM = hUnlike->GetAxis(6)->FindBin(0.005 - 0.00001);

    // hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    // hUnlike->GetAxis(6)->SetRange(lowDeltaM, highDeltaM);
    // TH1D *hInvMass = hUnlike->Projection(0, "E");
    // hInvMass->Rebin(8);
    // hInvMass->GetXaxis()->SetRangeUser(2.5, 2.9);

    TFile *fInputData = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResults_WithPhiMasses.root");
    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInputData, "doublephimeson/SEMassPhiPhiRefitted");
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
    hUnlike->GetAxis(2)->SetRange(lowDeltaM, highDeltaM);
    hUnlike->GetAxis(3)->SetRange(lowChi2, highChi2);
    hUnlike->GetAxis(4)->SetRange(lowFitProb, highFitProb);

    TH1D *hInvMass = hUnlike->Projection(0, "E");
    hInvMass->Rebin(8);
    hInvMass->GetXaxis()->SetRangeUser(fitLow, fitHigh);

    TH1D *h_Data = (TH1D *)hInvMass->Clone("h_Data");
    TH1D *h_Data_NoPeak = (TH1D *)hInvMass->Clone("h_Data_NoPeak");
    int totalBins = h_Data->GetNbinsX();
    for (int i = 1; i <= totalBins; i++)
    {
        double x = h_Data->GetBinCenter(i);
        if (x >= excLow && x <= excHigh)
        {
            h_Data_NoPeak->SetBinContent(i, 0);
            h_Data_NoPeak->SetBinError(i, 0);
        }
        else
        {
            h_Data_NoPeak->SetBinContent(i, h_Data->GetBinContent(i));
            h_Data_NoPeak->SetBinError(i, h_Data->GetBinError(i));
        }
    }

    // Total Background function with window rejection
    auto total_bkg_sideband_func = [f_SS_template, f_nonSS_template, excLow, excHigh](double *x, double *p) -> double
    {
        double m = x[0];
        if (m > excLow && m < excHigh)
        {
            TF1::RejectPoint();
            return 0.0;
        }
        return p[0] * f_SS_template->Eval(m) + p[1] * f_nonSS_template->Eval(m);
    };

    // Full range total background evaluation
    auto total_bkg_full_func = [f_SS_template, f_nonSS_template](double *x, double *p) -> double
    {
        double m = x[0];
        return p[0] * f_SS_template->Eval(m) + p[1] * f_nonSS_template->Eval(m);
    };

    TF1 *f_TotalBkg_Sideband = new TF1("f_TotalBkg_Sideband", total_bkg_sideband_func, fitLow, fitHigh, 2);
    f_TotalBkg_Sideband->SetParameters(1.0, 1.0);
    f_TotalBkg_Sideband->SetParNames("scale_SS", "scale_nonSS");
    h_Data_NoPeak->Fit(f_TotalBkg_Sideband, "R0Q");
    int parameter1 = f_TotalBkg_Sideband->GetParameter(0);
    int parameter2 = f_TotalBkg_Sideband->GetParameter(1);
    h_Data->Fit(f_TotalBkg_Sideband, "R0Q");

    TF1 *f_TotalBkg_Full = new TF1("f_TotalBkg_Full", total_bkg_full_func, fitLow, fitHigh, 2);
    f_TotalBkg_Full->SetParameters(f_TotalBkg_Sideband->GetParameters());
    f_TotalBkg_Full->SetLineColor(kRed);
    f_TotalBkg_Full->SetLineWidth(3);

    // Individual component components for visualization
    TF1 *f_SS_Component = new TF1("f_SS_Component", [f_SS_template](double *x, double *p)
                                  { return p[0] * f_SS_template->Eval(x[0]); }, fitLow, fitHigh, 1);
    f_SS_Component->SetParameter(0, f_TotalBkg_Sideband->GetParameter(0));
    f_SS_Component->SetLineColor(kBlue);
    f_SS_Component->SetLineStyle(3);

    TF1 *f_nonSS_Component = new TF1("f_nonSS_Component", [f_nonSS_template](double *x, double *p)
                                     { return p[0] * f_nonSS_template->Eval(x[0]); }, fitLow, fitHigh, 1);
    f_nonSS_Component->SetParameter(0, f_TotalBkg_Sideband->GetParameter(1));
    f_nonSS_Component->SetLineColor(kGreen + 2);
    f_nonSS_Component->SetLineStyle(7);

    // Plot Step 2 Results
    TCanvas *cBkgFit = new TCanvas("cBkgFit", "Total Background Fit", 720, 720);
    SetCanvasStyle(cBkgFit, 0.15, 0.03, 0.05, 0.15);
    h_Data->SetMarkerStyle(20);
    SetHistoQA(h_Data);
    h_Data->SetMinimum(480);
    h_Data->SetMaximum(5300);
    h_Data->GetXaxis()->SetTitle("M_{#phi#phi} (GeV/#it{c}^{2})");
    h_Data->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", h_Data->GetBinWidth(1) * 1000));
    h_Data->Draw("PE");

    f_TotalBkg_Full->Draw("SAME");
    f_SS_Component->SetLineWidth(3);
    f_SS_Component->Draw("SAME");
    f_nonSS_Component->SetLineWidth(3);
    f_nonSS_Component->Draw("SAME");

    TLegend *leg2 = new TLegend(0.55, 0.73, 0.88, 0.92);
    leg2->AddEntry(h_Data, "Data", "pe");
    leg2->AddEntry(f_TotalBkg_Full, "Total background", "l");
    leg2->AddEntry(f_SS_Component, "SS shape", "l");
    leg2->AddEntry(f_nonSS_Component, "non-SS background", "l");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->Draw();

    cBkgFit->SaveAs(savepath + "/Bkg_Fit=" + suffix + ".png");

    // Save fitted background models for Step 3 signal extraction
    TFile *fOut = new TFile(savepath + "/Step1_Step2_Background_Results.root", "RECREATE");
    f_TotalBkg_Full->Write("f_TotalBkg_Full");
    f_SS_Component->Write("f_SS_Component");
    f_nonSS_Component->Write("f_nonSS_Component");
    h_Data->Write("h_Data_1D");
    fOut->Close();

    // =========================================================================
    // STEP 3: Background Subtraction & Signal Extraction
    // =========================================================================

    // 1. Bin-by-bin background subtraction (Data - Background)
    TH1D *h_Subtracted = (TH1D *)h_Data->Clone("h_Subtracted");
    h_Subtracted->SetTitle("Background-subtracted pair mass;M_{#phi#phi} (GeV/#it{c}^{2});data - background");

    for (int bin = 1; bin <= h_Data->GetNbinsX(); ++bin)
    {
        double binCenter = h_Data->GetBinCenter(bin);
        double dataVal = h_Data->GetBinContent(bin);
        double dataErr = h_Data->GetBinError(bin);

        // Evaluate background function at the bin center
        double bkgVal = f_TotalBkg_Full->Eval(binCenter);

        // Subtraction keeping the data statistical error
        h_Subtracted->SetBinContent(bin, dataVal - bkgVal);
        h_Subtracted->SetBinError(bin, dataErr);
    }

    // 2. Define Voigtian function here
    auto VoigtFunction = [](double *x, double *p) -> double
    {
        double m = x[0];
        double N = p[0];
        double m0 = p[1];
        double sigma = p[2]; // Gaussian width
        double gamma = p[3]; // Lorentzian width

        return N * TMath::Voigt(m - m0, sigma, gamma);
    };

    // Fit function range focused around the signal peak region
    TF1 *f_Signal = new TF1("f_Signal", VoigtFunction, fitLow, fitHigh, 4);
    // TF1 *f_Signal = new TF1("f_Signal", VoigtFunction, 2.6, 2.8, 4);
    f_Signal->SetParNames("N_X", "M_X", "Sigma_X", "Gamma_X");

    // Parameter initialization based on PDF slide
    f_Signal->SetParameter(0, 1100.0);       // Yield initial guess
    f_Signal->SetParameter(1, 2.6905);       // M_X initial guess
    f_Signal->SetParLimits(1, 2.650, 2.730); // Mass peak search range
    f_Signal->FixParameter(2, 0.012);        // Gaussian width
    f_Signal->SetParameter(3, 0.012);         // Resonance width
    f_Signal->SetParLimits(3, 0.005, 0.05);  // Reasonance width
    f_Signal->SetLineColor(kMagenta + 2);
    f_Signal->SetLineWidth(2);

    // Fit signal peak
    SetHistoQA(h_Subtracted);
    h_Subtracted->GetYaxis()->SetRangeUser(-290, 480);
    h_Subtracted->Fit(f_Signal, "REN");

    // 3. Extract parameters and statistical significance
    double yield = f_Signal->GetParameter(0) / h_Subtracted->GetBinWidth(1); // Normalize yield by bin width
    double yieldErr = f_Signal->GetParError(0) / h_Subtracted->GetBinWidth(1);
    double mass = f_Signal->GetParameter(1);
    double massErr = f_Signal->GetParError(1);
    double width = f_Signal->GetParameter(3);
    double widthErr = f_Signal->GetParError(3);
    double significance = yield / yieldErr; // Raw statistical significance

    // 4. Plot Step 3 Results
    TCanvas *cSignal = new TCanvas("cSignal", "Signal Extraction", 720, 720);
    SetCanvasStyle(cSignal, 0.15, 0.05, 0.08, 0.12);

    SetHistoQA(h_Subtracted);
    h_Subtracted->SetMarkerStyle(20);
    h_Subtracted->SetMarkerSize(0.8);
    h_Subtracted->GetXaxis()->SetRangeUser(fitLow, fitHigh);
    h_Subtracted->Draw("PE");

    f_Signal->Draw("SAME");

    // Add zero reference line
    TF1 *f_Zero = new TF1("f_Zero", "0", fitLow, fitHigh);
    f_Zero->SetLineStyle(2);
    f_Zero->SetLineColor(kGray + 2);
    f_Zero->Draw("SAME");

    // Display extracted fit parameter text overlay
    TLegend *legSig = new TLegend(0.11, 0.72, 0.4, 0.91);
    legSig->SetBorderSize(0);
    legSig->SetFillStyle(0);
    legSig->SetTextSize(0.032);
    legSig->AddEntry((TObject *)0, Form("N_{X} = %.1f #pm %.1f", yield, yieldErr), "");
    legSig->AddEntry((TObject *)0, Form("M_{X} = %.4f #pm %.4f", mass, massErr), "");
    legSig->AddEntry((TObject *)0, Form("#Gamma_{X} = %.4f #pm %.4f", width, widthErr), "");
    legSig->AddEntry((TObject *)0, Form("Stat. Significance = %.2f", significance), "");
    legSig->Draw();

    TLegend *legSig2 = new TLegend(0.55, 0.78, 0.88, 0.91);
    legSig2->SetBorderSize(0);
    legSig2->SetFillStyle(0);
    legSig2->SetTextSize(0.032);
    legSig2->AddEntry(h_Subtracted, "Data - Background", "pe");
    legSig2->AddEntry(f_Signal, "Voigtian Fit", "l");
    legSig2->Draw();

    cSignal->SaveAs(savepath + "/TetraquarkPeakFit=" + suffix + ".png");

    // Output values to terminal
    std::cout << "\n========== STEP 3 FIT RESULTS ==========" << std::endl;
    std::cout << "Yield (N_X)    : " << yield << " +/- " << yieldErr << std::endl;
    std::cout << "Mass (M_X)     : " << mass << " +/- " << massErr << " GeV/c^2" << std::endl;
    std::cout << "Width (Gamma_X): " << width << " +/- " << widthErr << " GeV/c^2" << std::endl;
    std::cout << "Significance   : " << significance << std::endl;
    std::cout << "========================================" << std::endl;

    std::cout << "Code completed successfully!" << std::endl;
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