#include <iostream>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/initializations.h"
#include "YieldMean.C"

using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par)
{
    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void compare_yield_inel2()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    bool isSameBins = false;

    double f_norm_run3 = 0.677; //(52.8 / 77.904)

    string path1 = "../../output/kstar/LHC22o_pass7/788987/kstarqa/hInvMass/ROTATED";
    string pathSys = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/756343/kstarqa/hInvMass/SystematicsPlots";
    TString outputPath = path1 + "/spectra_compare";
    gSystem->mkdir(outputPath, kTRUE);

    TFile *fSys = new TFile((pathSys + "/SysUncert.root").c_str(), "READ");
    TH1D *hSys = (TH1D *)fSys->Get("hTotalSysSmoothed_0_120");

    //==============First lets calculate signal loss for INEL=======================
    TFile *fSigLoss = new TFile("../../mc/LHC24f3c/801014.root", "READ"); // 2024 MC INEL lf-mc-centrality

    double eventLoss;

    TFile *fspectra1 = new TFile((path1 + "/corrected_spectra_0_120.root").c_str(), "read");

    if (fspectra1->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    TFile *fpub = new TFile("pp13TeV_INEL.root", "READ");
    if (fpub->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }

    TH1F *hmultBase = (TH1F *)fspectra1->Get("mult_0-120/corrected_spectra_Integral_final");

    if (hmultBase == nullptr)
    {
        cout << "Error: histograms not found in the corrected spectra root file" << endl;
        return;
    }

    TH1F *hDen = (TH1F *)fSigLoss->Get("kstarqa/hInvMass/MCcorrections/MultiplicityGen");
    TH1F *hNum = (TH1F *)fSigLoss->Get("kstarqa/hInvMass/MCcorrections/MultiplicityRec");
    TH1F *hNumExcludingTVX = (TH1F *)fSigLoss->Get("kstarqa/hInvMass/MCcorrections/MultiplicityRecINEL");

    double eventLossNum = hNum->Integral();
    double eventLossDen = hDen->Integral();
    double eventLossNumExcludingTVX = hNumExcludingTVX->Integral();
    eventLoss = eventLossNum / eventLossDen;
    double eventLossExcludingTVX = eventLossNumExcludingTVX / eventLossDen;
    cout << "Event Loss for INEL: " << eventLoss << endl;
    cout << "Event Loss for INEL excluding TVX: " << eventLossExcludingTVX << endl;

    // Apply common scaling to base histogram
    hmultBase->Scale(2.0);             // for INEL it is K* + K*_bar and not their average
    hmultBase->Scale(f_norm_run3);     // visible cross section ratio
    hmultBase->Scale(1.0 / eventLoss); // removing event loss

    double ratioSel_TVX = 9235908442.0 / 12478172244.0;

    // --- 1. Data-based corrected spectrum ---
    TH1F *hDataCorr = (TH1F *)hmultBase->Clone("hDataCorr");
    hDataCorr->Scale(ratioSel_TVX); // selected events ==> TVX events

    // --- 2. MC-based corrected spectrum ---
    TH1F *hMCCorr = (TH1F *)hmultBase->Clone("hMCCorr");
    hMCCorr->Scale(eventLossExcludingTVX); // Event loss without TriggerTVX

    // ==============================================================================
    //   Save corrected INEL spectra (Data-based) with statistical and systematic uncertainties
    // ==============================================================================

    TH1F *hINEL_stat = (TH1F *)hDataCorr->Clone("hINEL_stat");
    hINEL_stat->SetTitle("INEL spectrum with statistical uncertainties");
    hINEL_stat->SetDirectory(0);

    TH1F *hINEL_sys = (TH1F *)hDataCorr->Clone("hINEL_sys");
    hINEL_sys->SetTitle("INEL spectrum with systematic uncertainties");
    hINEL_sys->SetDirectory(0);

    if (hSys->GetNbinsX() != hINEL_sys->GetNbinsX())
    {
        cout << "Error: Systematic histogram and spectrum have different bin counts!" << endl;
        return;
    }

    for (int i = 1; i <= hINEL_sys->GetNbinsX(); i++)
    {
        double sysError = hSys->GetBinContent(i) * hINEL_sys->GetBinContent(i);
        hINEL_sys->SetBinError(i, sysError);
    }

    // Prepare systematic error histograms for both options
    TH1F *hDataCorr_sys = (TH1F *)hDataCorr->Clone("hDataCorr_sys");
    TH1F *hMCCorr_sys = (TH1F *)hMCCorr->Clone("hMCCorr_sys");

    for (int i = 1; i <= hSys->GetNbinsX(); i++)
    {
        double err1 = hSys->GetBinContent(i) * hDataCorr->GetBinContent(i);
        hDataCorr_sys->SetBinError(i, err1);

        double err2 = hSys->GetBinContent(i) * hMCCorr->GetBinContent(i);
        hMCCorr_sys->SetBinError(i, err2);
    }

    TFile *fINEL = new TFile((path1 + "/INEL_spectra.root").c_str(), "RECREATE");
    hDataCorr->Write("hINEL_DataCorr_stat");
    hDataCorr_sys->Write("hINEL_DataCorr_sys");
    hMCCorr->Write("hINEL_MCCorr_stat");
    hMCCorr_sys->Write("hINEL_MCCorr_sys");
    fINEL->Close();

    cout << "Saved INEL spectra with statistical and systematic uncertainties: " << path1 << "/INEL_spectra.root" << endl;

    TGraphErrors *gRun2_minBias = (TGraphErrors *)fpub->Get("Table 4/Graph1D_y1");
    if (gRun2_minBias == nullptr)
    {
        cout << "Run2 minBias graph not found" << endl;
        return;
    }

    // --- Levy Fit Functions ---
    TF1 *fitFcn1 = new TF1("fitfunc1", FuncLavy, 0.0, 15.0, 4);
    fitFcn1->SetParameter(0, 5.0);
    fitFcn1->SetParameter(1, 0.07);
    fitFcn1->FixParameter(2, 0.895);
    fitFcn1->SetParameter(3, 0.3);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");

    TF1 *fitFcn2 = new TF1("fitfunc2", FuncLavy, 0.0, 15.0, 4);
    fitFcn2->SetParameter(0, 5.0);
    fitFcn2->SetParameter(1, 0.07);
    fitFcn2->FixParameter(2, 0.895);
    fitFcn2->SetParameter(3, 0.3);
    fitFcn2->SetParNames("n", "dn/dy", "mass", "T");

    Double_t min = 0;
    Double_t max = 15;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 15;

    // Perform fits using YieldMean procedure
    TH1 *hout1 = YieldMean(hDataCorr, hDataCorr_sys, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *hout2 = YieldMean(hMCCorr, hMCCorr_sys, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    // --- Graphs for Ratios vs Published 13 TeV Data ---
    TGraphErrors *gratio1 = new TGraphErrors();
    TGraphErrors *gratio1_sys = new TGraphErrors();

    TGraphErrors *gratio2 = new TGraphErrors();
    TGraphErrors *gratio2_sys = new TGraphErrors();

    if (!isSameBins)
    {
        for (int i = 0; i < gRun2_minBias->GetN(); i++)
        {
            double x_run2, yield_run2, x_error, y_error_run2;
            gRun2_minBias->GetPoint(i, x_run2, yield_run2);
            x_error = gRun2_minBias->GetErrorX(i);
            y_error_run2 = gRun2_minBias->GetErrorY(i);

            int sysBin = hSys->GetXaxis()->FindBin(x_run2);

            // Data-based ratio calculation
            double eval1 = fitFcn1->Eval(x_run2);
            double sysError_fit1 = hSys->GetBinContent(sysBin) * eval1;
            double error1_stat = sqrt(pow(eval1 * y_error_run2 / (yield_run2 * yield_run2), 2));
            double error1_sys = sqrt(pow(sysError_fit1 / yield_run2, 2) + pow(eval1 * y_error_run2 / (yield_run2 * yield_run2), 2));

            gratio1->SetPoint(i, x_run2, eval1 / yield_run2);
            gratio1->SetPointError(i, x_error, error1_stat);
            gratio1_sys->SetPoint(i, x_run2, eval1 / yield_run2);
            double xBand = (x_error > 0) ? x_error : 0.10;
            gratio1_sys->SetPointError(i, xBand, error1_sys);

            // MC-based ratio calculation
            double eval2 = fitFcn2->Eval(x_run2);
            double sysError_fit2 = hSys->GetBinContent(sysBin) * eval2;
            double error2_stat = sqrt(pow(eval2 * y_error_run2 / (yield_run2 * yield_run2), 2));
            double error2_sys = sqrt(pow(sysError_fit2 / yield_run2, 2) + pow(eval2 * y_error_run2 / (yield_run2 * yield_run2), 2));

            gratio2->SetPoint(i, x_run2, eval2 / yield_run2);
            gratio2->SetPointError(i, x_error, error2_stat);
            gratio2_sys->SetPoint(i, x_run2, eval2 / yield_run2);
            gratio2_sys->SetPointError(i, xBand, error2_sys);
        }
    }
    else
    {
        if (gRun2_minBias->GetN() != hDataCorr->GetNbinsX())
        {
            cout << "Error: Number of points in Run2 minBias graph does not match number of bins in histogram." << endl;
            return;
        }
        for (int i = 0; i < gRun2_minBias->GetN(); i++)
        {
            double x_run2, yield_run2, x_error, y_error_run2;
            gRun2_minBias->GetPoint(i, x_run2, yield_run2);
            x_error = gRun2_minBias->GetErrorX(i);
            y_error_run2 = gRun2_minBias->GetErrorY(i);

            // Data-based
            double yieldError_run3_1 = hDataCorr->GetBinError(i + 1);
            double sysError_run3_1 = hSys->GetBinError(i + 1) * hDataCorr->GetBinContent(i + 1);
            double yield_Run3_1 = hDataCorr->GetBinContent(i + 1);

            gratio1->SetPoint(i, x_run2, yield_Run3_1 / yield_run2);
            gratio1_sys->SetPoint(i, x_run2, yield_Run3_1 / yield_run2);
            double error1 = sqrt(pow(yieldError_run3_1 / yield_run2, 2) + pow(yield_Run3_1 * y_error_run2 / (yield_run2 * yield_run2), 2));
            double error1_sys = sqrt(pow(sysError_run3_1 / yield_run2, 2) + pow(yield_Run3_1 * y_error_run2 / (yield_run2 * yield_run2), 2));
            gratio1->SetPointError(i, x_error, error1);
            double xBand = 0.5 * hDataCorr->GetXaxis()->GetBinWidth(i + 1);
            gratio1_sys->SetPointError(i, xBand, error1_sys);

            // MC-based
            double yieldError_run3_2 = hMCCorr->GetBinError(i + 1);
            double sysError_run3_2 = hSys->GetBinError(i + 1) * hMCCorr->GetBinContent(i + 1);
            double yield_Run3_2 = hMCCorr->GetBinContent(i + 1);

            gratio2->SetPoint(i, x_run2, yield_Run3_2 / yield_run2);
            gratio2_sys->SetPoint(i, x_run2, yield_Run3_2 / yield_run2);
            double error2 = sqrt(pow(yieldError_run3_2 / yield_run2, 2) + pow(yield_Run3_2 * y_error_run2 / (yield_run2 * yield_run2), 2));
            double error2_sys = sqrt(pow(sysError_run3_2 / yield_run2, 2) + pow(yield_Run3_2 * y_error_run2 / (yield_run2 * yield_run2), 2));
            gratio2->SetPointError(i, x_error, error2);
            gratio2_sys->SetPointError(i, xBand, error2_sys);
        }
    }

    // --- Plotting Canvas 1: Spectra and Ratios vs Published 13 TeV ---
    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);

    // --- Pad 1: Spectra Comparison ---
    c1->cd(1);
    SetHistoStyle(hDataCorr, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    hDataCorr->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hDataCorr->SetMaximum(hDataCorr->GetMaximum() * 5);
    hDataCorr->SetMinimum(hDataCorr->GetMinimum() * 0.15);
    hDataCorr->GetYaxis()->SetTitleOffset(1.30);
    hDataCorr->GetXaxis()->SetTitleOffset(1.02);
    hDataCorr->GetXaxis()->SetRangeUser(0, 15);

    // Data-based style (Blue)
    hDataCorr->SetMarkerStyle(20);
    hDataCorr->SetMarkerSize(1);
    hDataCorr->SetLineColor(kBlue);
    hDataCorr->SetMarkerColor(kBlue);
    hDataCorr->Draw("pe");

    hDataCorr_sys->SetMarkerColor(kBlue);
    hDataCorr_sys->SetLineColor(kBlue);
    hDataCorr_sys->SetFillStyle(0);
    hDataCorr_sys->SetLineWidth(1);
    hDataCorr_sys->Draw("e2 same");

    fitFcn1->SetLineColor(kBlue);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);
    fitFcn1->Draw("same");

    // MC-based style (Red)
    hMCCorr->SetMarkerStyle(21);
    hMCCorr->SetMarkerSize(1);
    hMCCorr->SetLineColor(kRed);
    hMCCorr->SetMarkerColor(kRed);
    hMCCorr->Draw("pe same");

    hMCCorr_sys->SetMarkerColor(kRed);
    hMCCorr_sys->SetLineColor(kRed);
    hMCCorr_sys->SetFillStyle(0);
    hMCCorr_sys->SetLineWidth(1);
    hMCCorr_sys->Draw("e2 same");

    fitFcn2->SetLineColor(kRed);
    fitFcn2->SetLineStyle(2);
    fitFcn2->SetLineWidth(2);
    fitFcn2->Draw("same");

    gPad->SetLogy(1);
    gRun2_minBias->SetMarkerStyle(22);
    gRun2_minBias->SetMarkerSize(1);
    gRun2_minBias->SetMarkerColor(kBlack);
    gRun2_minBias->SetLineColor(kBlack);
    gRun2_minBias->SetLineWidth(2);
    gRun2_minBias->Draw("pe same");

    TLegend *leg = new TLegend(0.20, 0.04, 0.65, 0.30);
    SetLegendStyle(leg);
    leg->SetTextSize(0.04);
    leg->AddEntry(hDataCorr, "pp 13.6 TeV (Data-based corr)", "p");
    leg->AddEntry(hMCCorr, "pp 13.6 TeV (MC-based corr)", "p");
    leg->AddEntry(gRun2_minBias, "pp 13 TeV (Eur. Phys. J. C 81:256)", "p");
    leg->AddEntry(fitFcn1, "Levy-Tsallis Fit", "l");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.SetTextFont(42);
    lat.DrawLatex(0.65, 0.90, "ALICE");
    lat.DrawLatex(0.65, 0.83, "|y| < 0.5");
    lat.DrawLatex(0.65, 0.76, "INEL");
    lat.DrawLatex(0.65, 0.69, "K*(892)^{0}");

    // --- Pad 2: Ratios vs Published ---
    c1->cd(2);
    SetGraphErrorStyle(gratio1);
    gratio1->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    gratio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    gratio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    gratio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    gratio1->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
    gratio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gratio1->GetXaxis()->CenterTitle(1);
    gratio1->GetYaxis()->SetTitleOffset(0.6);
    gratio1->GetXaxis()->SetTitleOffset(1.1);
    gratio1->GetYaxis()->SetNdivisions(506);
    gratio1->GetXaxis()->SetRangeUser(0, 15);
    gratio1->SetMinimum(0.2);
    gratio1->SetMaximum(1.95);

    // Draw Data-based Ratio (Blue)
    gratio1->SetMarkerStyle(20);
    gratio1->SetMarkerSize(1.0);
    gratio1->SetMarkerColor(kBlue);
    gratio1->SetLineColor(kBlue);
    gratio1->Draw("ap");

    gratio1_sys->SetMarkerStyle(20);
    gratio1_sys->SetMarkerSize(1.0);
    gratio1_sys->SetMarkerColor(kBlue);
    gratio1_sys->SetLineColor(kBlue);
    gratio1_sys->SetFillColor(kBlue);
    gratio1_sys->SetLineWidth(2);
    gratio1_sys->SetFillStyle(0);
    gratio1_sys->Draw("e2 same");

    // Draw MC-based Ratio (Red)
    SetGraphErrorStyle(gratio2);
    gratio2->SetMarkerStyle(21);
    gratio2->SetMarkerSize(1.0);
    gratio2->SetMarkerColor(kRed);
    gratio2->SetLineColor(kRed);
    gratio2->Draw("p same");

    gratio2_sys->SetMarkerStyle(21);
    gratio2_sys->SetMarkerSize(1.0);
    gratio2_sys->SetMarkerColor(kRed);
    gratio2_sys->SetLineColor(kRed);
    gratio2_sys->SetFillColor(kRed);
    gratio2_sys->SetLineWidth(2);
    gratio2_sys->SetFillStyle(0);
    gratio2_sys->Draw("e2 same");

    TLine *line = new TLine(0, 1.01, 15, 1.01);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(1);
    line->Draw();

    // Grey band with 20% uncertainty around 1
    TBox *box = new TBox(0, 0.9, 15, 1.1);
    box->SetFillColor(kGray + 2);
    box->SetFillStyle(3003);
    box->Draw("same");

    TString saveName = "YieldRatioINEL_newMC.png";
    c1->SaveAs(outputPath + "/" + saveName);

    // ==============================================================================
    //   Canvas 2: Direct Ratio of Data-based / MC-based Histograms (No Levy Fit)
    // ==============================================================================

    TCanvas *c2 = new TCanvas("c2", "c2 Direct Ratio", 720, 600);
    SetCanvasStyle(c2, 0.15, 0.05, 0.05, 0.15);
    c2->cd();

    TH1F *hDirectRatio = (TH1F *)hDataCorr->Clone("hDirectRatio");
    hDirectRatio->Divide(hMCCorr);

    TH1F *hDirectRatio_sys = (TH1F *)hDataCorr_sys->Clone("hDirectRatio_sys");
    hDirectRatio_sys->Divide(hMCCorr_sys);

    SetHistoStyle(hDirectRatio, 1, 53, 1, 0.04, 0.04, 0.04, 0.04, 1.2, 1.2);
    hDirectRatio->SetMarkerStyle(20);
    hDirectRatio->SetMarkerSize(1.0);
    hDirectRatio->SetMarkerColor(kBlack);
    hDirectRatio->SetLineColor(kBlack);
    hDirectRatio->GetYaxis()->SetTitle("#frac{Data-based}{MC-based}");
    hDirectRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hDirectRatio->GetXaxis()->CenterTitle(1);
    hDirectRatio->GetXaxis()->SetRangeUser(0, 15);
    hDirectRatio->SetMinimum(0.5);
    hDirectRatio->SetMaximum(1.5);
    hDirectRatio->Draw("pe");

    hDirectRatio_sys->SetMarkerStyle(20);
    hDirectRatio_sys->SetMarkerSize(1.0);
    hDirectRatio_sys->SetMarkerColor(kBlack);
    hDirectRatio_sys->SetLineColor(kBlack);
    hDirectRatio_sys->SetFillColor(kBlack);
    hDirectRatio_sys->SetLineWidth(1);
    hDirectRatio_sys->SetFillStyle(0);
    hDirectRatio_sys->Draw("e2 same");

    TLine *line2 = new TLine(0, 1.0, 15, 1.0);
    line2->SetLineStyle(2);
    line2->SetLineWidth(2);
    line2->SetLineColor(kRed);
    line2->Draw("same");

    TLegend *leg2 = new TLegend(0.20, 0.75, 0.55, 0.88);
    SetLegendStyle(leg2);
    leg2->SetTextSize(0.04);
    leg2->AddEntry(hDirectRatio, "Data-based / MC-based", "pe");
    leg2->Draw();

    lat.DrawLatex(0.65, 0.85, "ALICE");
    lat.DrawLatex(0.65, 0.78, "|y| < 0.5");
    lat.DrawLatex(0.65, 0.71, "INEL pp 13.6 TeV");

    TString saveName2 = "DirectRatio_Data_vs_MC.png";
    c2->SaveAs(outputPath + "/" + saveName2);
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}