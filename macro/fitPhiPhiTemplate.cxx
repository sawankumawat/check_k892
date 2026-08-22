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

void fitPhiPhiTemplate()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    // Paths
    TString inFilePath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template/DoublePhiBackgroundTemplates.root";
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template";

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
    auto chebyshev3_sideband = [excLow, excHigh](double *x, double *p) -> double
    {
        double m = x[0];
        if (m > excLow && m < excHigh)
        {
            TF1::RejectPoint();
            return 0.0;
        }
        // Normalize x from [2.5, 2.9] to [-1, 1]
        double x_norm = (2.0 * (m - 2.5) / (2.9 - 2.5)) - 1.0;
        double T0 = 1.0;
        double T1 = x_norm;
        double T2 = 2.0 * x_norm * x_norm - 1.0;
        double T3 = 4.0 * x_norm * x_norm * x_norm - 3.0 * x_norm;

        return p[0] * T0 + p[1] * T1 + p[2] * T2 + p[3] * T3;
    };

    // Full range evaluation for interpolation through signal region
    auto chebyshev3_full = [](double *x, double *p) -> double
    {
        double m = x[0];
        double x_norm = (2.0 * (m - 2.5) / (2.9 - 2.5)) - 1.0;
        double T0 = 1.0;
        double T1 = x_norm;
        double T2 = 2.0 * x_norm * x_norm - 1.0;
        double T3 = 4.0 * x_norm * x_norm * x_norm - 3.0 * x_norm;

        return p[0] * T0 + p[1] * T1 + p[2] * T2 + p[3] * T3;
    };

    TF1 *f_SS_sideband = new TF1("f_SS_sideband", chebyshev3_sideband, 2.5, 2.9, 4);
    TF1 *f_nonSS_sideband = new TF1("f_nonSS_sideband", chebyshev3_sideband, 2.5, 2.9, 4);

    f_SS_sideband->SetParameters(1500, -300, 50, -10);
    f_nonSS_sideband->SetParameters(1500, -50, 10, -5);

    h_N_SS->Fit(f_SS_sideband, "R0Q");
    h_N_nonSS->Fit(f_nonSS_sideband, "R0Q");

    // Extrapolated templates across full mass range
    TF1 *f_SS_template = new TF1("f_SS_template", chebyshev3_full, 2.5, 2.9, 4);
    TF1 *f_nonSS_template = new TF1("f_nonSS_template", chebyshev3_full, 2.5, 2.9, 4);

    f_SS_template->SetParameters(f_SS_sideband->GetParameters());
    f_nonSS_template->SetParameters(f_nonSS_sideband->GetParameters());

    // Plot Step 1 Results
    TCanvas *cTemplates = new TCanvas("cTemplates", "Background Templates", 1280, 480);
    SetCanvasStyle(cTemplates, 0.15, 0.03, 0.05, 0.15);
    cTemplates->Divide(2, 1);

    cTemplates->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.10);

    h_N_SS->SetTitle("SS template (Chebyshev-3);M_{#phi#phi} (GeV/c^{2});selected N_{SS}");
    h_N_SS->SetMarkerStyle(20);
    SetHistoQA(h_N_SS);
    h_N_SS->SetTitle("SS template (Chebyshev-3)");
    h_N_SS->Draw("PE");
    f_SS_sideband->SetLineColor(kRed);
    f_SS_sideband->SetLineWidth(2);
    f_SS_sideband->Draw("SAME");
    f_SS_template->SetLineStyle(2); // Dotted line interpolation
    f_SS_template->SetLineColor(kRed);
    f_SS_template->Draw("SAME");

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
    f_nonSS_sideband->SetLineWidth(2);
    f_nonSS_sideband->Draw("SAME");
    f_nonSS_template->SetLineStyle(2);
    f_nonSS_template->SetLineColor(kRed);
    f_nonSS_template->Draw("SAME");

    cTemplates->SaveAs(savepath + "/Background_Templates.png");

    // =========================================================================
    // STEP 2: Combine Templates to Form Total Background Fit
    // =========================================================================

    TFile *fInputData = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults_MorePhiBins.root");
    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInputData, "doublephimeson/SEMassDoublePhi");

    int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowDeltaM = hUnlike->GetAxis(6)->FindBin(0.0 + 0.00001);
    int highDeltaM = hUnlike->GetAxis(6)->FindBin(0.005 - 0.00001);

    hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    hUnlike->GetAxis(6)->SetRange(lowDeltaM, highDeltaM);
    TH1D *hInvMass = hUnlike->Projection(0, "E");
    hInvMass->Rebin(8);
    hInvMass->GetXaxis()->SetRangeUser(2.5, 2.9);

    // Construct Total Mass spectrum (Sum of extracted SS and non-SS yields as pseudo-data)
    // TH1D *h_Data = (TH1D *)h_N_SS->Clone("h_Data");
    // h_Data->Add(h_N_nonSS);
    // h_Data->SetTitle("Background Template Fit;M_{#phi#phi} (GeV/c^{2});Counts");

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

    TF1 *f_TotalBkg_Sideband = new TF1("f_TotalBkg_Sideband", total_bkg_sideband_func, 2.5, 2.9, 2);
    f_TotalBkg_Sideband->SetParameters(1.0, 1.0);
    f_TotalBkg_Sideband->SetParNames("scale_SS", "scale_nonSS");
    h_Data_NoPeak->Fit(f_TotalBkg_Sideband, "R0Q");
    int parameter1 = f_TotalBkg_Sideband->GetParameter(0);
    int parameter2 = f_TotalBkg_Sideband->GetParameter(1);
    h_Data->Fit(f_TotalBkg_Sideband, "R0Q");

    TF1 *f_TotalBkg_Full = new TF1("f_TotalBkg_Full", total_bkg_full_func, 2.5, 2.9, 2);
    f_TotalBkg_Full->SetParameters(f_TotalBkg_Sideband->GetParameters());
    f_TotalBkg_Full->SetLineColor(kRed);
    f_TotalBkg_Full->SetLineWidth(2);

    // Individual component components for visualization
    TF1 *f_SS_Component = new TF1("f_SS_Component", [f_SS_template](double *x, double *p)
                                  { return p[0] * f_SS_template->Eval(x[0]); }, 2.5, 2.9, 1);
    f_SS_Component->SetParameter(0, f_TotalBkg_Sideband->GetParameter(0));
    f_SS_Component->SetLineColor(kBlue);
    f_SS_Component->SetLineStyle(3);

    TF1 *f_nonSS_Component = new TF1("f_nonSS_Component", [f_nonSS_template](double *x, double *p)
                                     { return p[0] * f_nonSS_template->Eval(x[0]); }, 2.5, 2.9, 1);
    f_nonSS_Component->SetParameter(0, f_TotalBkg_Sideband->GetParameter(1));
    f_nonSS_Component->SetLineColor(kGreen + 2);
    f_nonSS_Component->SetLineStyle(7);

    // Plot Step 2 Results
    TCanvas *cBkgFit = new TCanvas("cBkgFit", "Total Background Fit", 720, 720);
    SetCanvasStyle(cBkgFit, 0.15, 0.03, 0.05, 0.15);
    h_Data->SetMarkerStyle(20);
    SetHistoQA(h_Data);
    h_Data->SetMinimum(980);
    h_Data->SetMaximum(4300);
    h_Data->GetXaxis()->SetTitle("M_{#phi#phi} (GeV/#it{c}^{2})");
    h_Data->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", h_Data->GetBinWidth(1) * 1000));
    h_Data->Draw("PE");

    f_TotalBkg_Full->Draw("SAME");
    f_SS_Component->SetLineWidth(3);
    f_SS_Component->Draw("SAME");
    f_nonSS_Component->SetLineWidth(3);
    f_nonSS_Component->Draw("SAME");

    TLegend *leg = new TLegend(0.55, 0.73, 0.88, 0.92);
    leg->AddEntry(h_Data, "Data", "pe");
    leg->AddEntry(f_TotalBkg_Full, "Total background", "l");
    leg->AddEntry(f_SS_Component, "SS shape", "l");
    leg->AddEntry(f_nonSS_Component, "non-SS background", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    cBkgFit->SaveAs(savepath + "/Total_Background_Fit.png");

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

    // 2. Define Breit-Wigner signal fit function
    // p[0] = Signal Yield (N_X), p[1] = Mass (M_X), p[2] = Width (Gamma_X)
    auto BW = [](double *x, double *p) -> double
    {
        double m = x[0];
        double N = p[0];
        double m0 = p[1];
        double gamma = p[2];

        double denominator = (m - m0) * (m - m0) + (gamma * gamma) / 4.0;
        return (N * (gamma / 2.0)) / denominator;
    };

    // Fit function range focused around the signal peak region
    TF1 *f_Signal = new TF1("f_Signal", BW, 2.50, 2.90, 3);
    f_Signal->SetParNames("N_X", "M_X", "Gamma_X");

    // Parameter initialization based on PDF slide
    f_Signal->SetParameter(0, 1100.0);       // Yield initial guess
    f_Signal->SetParameter(1, 2.6905);       // M_X initial guess
    f_Signal->SetParLimits(1, 2.650, 2.730); // Mass peak search range
    f_Signal->SetParameter(2, 0.0269);       // Gamma_X initial guess
    f_Signal->SetParLimits(2, 0.005, 0.080); // Reasonable width limits
    f_Signal->SetLineColor(kMagenta + 2);
    f_Signal->SetLineWidth(2);

    // Fit signal peak
    SetHistoQA(h_Subtracted);
    h_Subtracted->GetYaxis()->SetRangeUser(-230, 290);
    h_Subtracted->Fit(f_Signal, "R0Q");

    // 3. Extract parameters and statistical significance
    double yield = f_Signal->GetParameter(0) / h_Subtracted->GetBinWidth(1); // Normalize yield by bin width
    double yieldErr = f_Signal->GetParError(0) / h_Subtracted->GetBinWidth(1);
    double mass = f_Signal->GetParameter(1);
    double massErr = f_Signal->GetParError(1);
    double width = f_Signal->GetParameter(2);
    double widthErr = f_Signal->GetParError(2);
    double significance = yield / yieldErr; // Raw statistical significance

    // 4. Plot Step 3 Results
    TCanvas *cSignal = new TCanvas("cSignal", "Signal Extraction", 720, 720);
    SetCanvasStyle(cSignal, 0.15, 0.05, 0.08, 0.12);

    SetHistoQA(h_Subtracted);
    h_Subtracted->SetMarkerStyle(20);
    h_Subtracted->SetMarkerSize(0.8);
    h_Subtracted->GetXaxis()->SetRangeUser(2.5, 2.9);
    h_Subtracted->Draw("PE");

    f_Signal->Draw("SAME");

    // Add zero reference line
    TF1 *f_Zero = new TF1("f_Zero", "0", 2.5, 2.9);
    f_Zero->SetLineStyle(2);
    f_Zero->SetLineColor(kGray + 2);
    f_Zero->Draw("SAME");

    // Display extracted fit parameter text overlay
    TLegend *legSig = new TLegend(0.15, 0.76, 0.5, 0.91);
    legSig->SetBorderSize(0);
    legSig->SetFillStyle(0);
    legSig->SetTextSize(0.032);
    // legSig->AddEntry((TObject*)0, Form("N_{X} = %.1f #pm %.1f", yield, yieldErr), "");
    legSig->AddEntry((TObject *)0, Form("M_{X} = %.4f #pm %.4f", mass, massErr), "");
    legSig->AddEntry((TObject *)0, Form("#Gamma_{X} = %.4f #pm %.4f", width, widthErr), "");
    legSig->AddEntry((TObject*)0, Form("yield / error = %.2f", significance), "");
    legSig->Draw();

    cSignal->SaveAs(savepath + "/Signal_Extraction.png");

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