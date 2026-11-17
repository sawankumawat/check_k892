#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

TFile *OpenFile(const string &path);
TGraphErrors *GetGraph(TFile *f, const string &name);
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void fitdNdy_meanpT()
{
    string KstarPath = "../../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/";
    TFile *fKstar = OpenFile(KstarPath + "Results.root");
    TFile *fPhi = OpenFile("PhiResults/Sawan/ResultsPhi.root");
    TFile *fpp13TeV = OpenFile("ConversionCodes/pp13TeV_KstarPhi.root");

    TGraphErrors *gMPtKstar[3], *gMYieldKstar[3], *gKstar_MeanpT_13TeV[3], *gKstar_Yield_13TeV[3], *gPhi_MeanpT_13TeV[3], *gMPtPhi[3];

    for (int i = 0; i < 3; i++)
    {
        string suffix;
        if (i == 0)
            suffix = "_stat";
        else if (i == 1)
            suffix = "_sys";
        else
            suffix = "_sysuncorr";

        gMPtKstar[i] = GetGraph(fKstar, Form("gMeanpTRun3%s", suffix.c_str()));
        gMYieldKstar[i] = GetGraph(fKstar, Form("gMeanYieldRun3%s", suffix.c_str()));

        gKstar_MeanpT_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_MeanpT%s", suffix.c_str()));
        gKstar_Yield_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_MeanYield%s", suffix.c_str()));

        gMPtPhi[i] = GetGraph(fPhi, Form("gMeanpTRun3%s", suffix.c_str()));
        gPhi_MeanpT_13TeV[i] = GetGraph(fpp13TeV, Form("gPhi_MeanpT%s", suffix.c_str()));
    }

    TF1 *fpol1Run2 = new TF1("fpol1Run2", "[0] + [1]*x", 0, 25);
    fpol1Run2->SetLineColor(kBlue + 1);
    fpol1Run2->SetLineWidth(2);
    fpol1Run2->SetLineStyle(2);
    TF1 *fpol1Run3 = new TF1("fpol1Run3", "[0] + [1]*x", 0, 25);
    fpol1Run3->SetLineColor(kRed + 1);
    fpol1Run3->SetLineWidth(2);

    TCanvas *cdNdyKstar = new TCanvas("cdNdyKstar", "cdNdyKstar", 720, 720);
    SetCanvasStyle(cdNdyKstar, 0.15, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cdNdyKstar, pad1Size, pad2Size);
    cdNdyKstar->cd(1);

    gMYieldKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMYieldKstar[0]->GetYaxis()->SetTitle("dN/dy");
    // gMYieldKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldKstar[0]->GetXaxis()->SetLimits(0, 25.5);
    gMYieldKstar[0]->SetLineWidth(3);
    gMYieldKstar[0]->GetYaxis()->SetRangeUser(0, 0.71);
    gMYieldKstar[0]->SetMarkerColor(kRed + 1);
    gMYieldKstar[0]->SetLineColor(kRed + 1);
    gMYieldKstar[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    gMYieldKstar[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    gMYieldKstar[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    gMYieldKstar[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    gMYieldKstar[0]->Draw("APE");
    gMYieldKstar[1]->SetFillStyle(0);
    gMYieldKstar[1]->SetLineColor(kRed + 1);
    gMYieldKstar[1]->SetLineWidth(3);
    gMYieldKstar[1]->Draw("5 same");
    gMYieldKstar[2]->SetFillStyle(1001);
    gMYieldKstar[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gMYieldKstar[2]->Draw("E2 same");
    gMYieldKstar[1]->Fit(fpol1Run3, "R");

    gKstar_Yield_13TeV[0]->SetMarkerStyle(21);
    gKstar_Yield_13TeV[0]->SetMarkerColor(kBlue + 1);
    gKstar_Yield_13TeV[0]->SetLineColor(kBlue + 1);
    gKstar_Yield_13TeV[0]->SetLineWidth(3);
    gKstar_Yield_13TeV[0]->Draw("P same");
    gKstar_Yield_13TeV[1]->SetLineColor(kBlue + 1);
    gKstar_Yield_13TeV[1]->SetFillStyle(0);
    gKstar_Yield_13TeV[1]->SetLineWidth(3);
    gKstar_Yield_13TeV[1]->Draw("5 same");
    gKstar_Yield_13TeV[2]->SetFillStyle(1001);
    gKstar_Yield_13TeV[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gKstar_Yield_13TeV[2]->Draw("E2 same");
    gKstar_Yield_13TeV[1]->Fit(fpol1Run2, "R");

    double slopeRun3 = fpol1Run3->GetParameter(1);
    double slopeRun3Err = fpol1Run3->GetParError(1);
    double interceptRun3 = fpol1Run3->GetParameter(0);
    double interceptRun3Err = fpol1Run3->GetParError(0);

    double slopeRun2 = fpol1Run2->GetParameter(1);
    double slopeRun2Err = fpol1Run2->GetParError(1);
    double interceptRun2 = fpol1Run2->GetParameter(0);
    double interceptRun2Err = fpol1Run2->GetParError(0);

    TLegend *legdNdyKstar = new TLegend(0.6, 0.76, 0.9, 0.9);
    legdNdyKstar->SetTextSize(0.04);
    legdNdyKstar->SetBorderSize(0);
    // legdNdyKstar->SetHeader("Uncertainty: Total systematic");
    // legdNdyKstar->SetHeader("Uncertainty: Statistical");
    legdNdyKstar->SetHeader("Uncertainty: Uncorrelated");
    legdNdyKstar->AddEntry(gMYieldKstar[0], "Run 3 (13.6 TeV)", "lp");
    legdNdyKstar->AddEntry(gKstar_Yield_13TeV[0], "Run 2 (13 TeV)", "lp");
    legdNdyKstar->Draw();

    TLatex *latdNdyKstar = new TLatex();
    latdNdyKstar->SetTextSize(0.04);
    latdNdyKstar->SetTextFont(22);
    latdNdyKstar->DrawLatexNDC(0.20, 0.9, "Fit (pp 13.6 TeV):");
    latdNdyKstar->SetTextFont(42);
    latdNdyKstar->DrawLatexNDC(0.20, 0.85, Form("Slope = %.4f #pm %.4f", slopeRun3, slopeRun3Err));
    latdNdyKstar->DrawLatexNDC(0.20, 0.8, Form("Intercept = %.4f #pm %.4f", interceptRun3, interceptRun3Err));
    latdNdyKstar->SetTextFont(22);
    latdNdyKstar->DrawLatexNDC(0.20, 0.72, "Fit (pp 13 TeV):");
    latdNdyKstar->SetTextFont(42);
    latdNdyKstar->DrawLatexNDC(0.20, 0.67, Form("Slope = %.4f #pm %.4f", slopeRun2, slopeRun2Err));
    latdNdyKstar->DrawLatexNDC(0.20, 0.62, Form("Intercept = %.4f #pm %.4f", interceptRun2, interceptRun2Err));

    cdNdyKstar->cd(2);
    // plot the ratio of run3/run2 where run2 values are taken from the fit
    TGraphErrors *gRatio = new TGraphErrors();
    for (int i = 0; i < gMYieldKstar[0]->GetN(); i++)
    {
        double x = gMYieldKstar[0]->GetX()[i];
        double yRun3 = gMYieldKstar[0]->GetY()[i];
        double yRun2 = fpol1Run2->Eval(x);
        double ratio = yRun3 / yRun2;

        double errRun3 = gMYieldKstar[0]->GetEY()[i];
        double errRun2 = fpol1Run2->GetParError(0) + fpol1Run2->GetParError(1) * x; // simplistic error propagation
        double errRatio = ratio * sqrt(pow(errRun3 / yRun3, 2) + pow(errRun2 / yRun2, 2));

        gRatio->SetPoint(i, x, ratio);
        gRatio->SetPointError(i, 0, errRatio);
    }
    SetGraphErrorStyle(gRatio);
    gRatio->SetMarkerStyle(21);
    gRatio->SetMarkerColor(kRed);
    gRatio->SetLineColor(kRed);
    gRatio->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gRatio->GetYaxis()->SetTitle("Run 3 / Run 2");
    gRatio->GetXaxis()->SetLimits(0, 25.5);
    gRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    gRatio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    gRatio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    gRatio->GetYaxis()->SetTitleSize(0.032 / pad2Size);
    gRatio->GetXaxis()->SetTitleSize(0.035 / pad2Size);
    gRatio->GetYaxis()->SetTitleOffset(1.9 * pad2Size);
    gRatio->GetYaxis()->SetNdivisions(505);
    gRatio->Draw("APE");
    TLine *line = new TLine(0, 1, 25.5, 1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("same");
    cdNdyKstar->SaveAs("Plots/dNdyKstarFit_ratio.png");

    //============================================================
    // ================<pT> K*0 fit ======================
    //============================================================

    // Fit <pT> 13 TeV with sigmoid (same function, parameters, and fit option)
    TF1 *fsigmoidRun2 = new TF1("fSigmoid", "[0]+[1]/(1+TMath::Exp(-(x-[2])/[3]))", 2.3, 28);
    fsigmoidRun2->SetParameters(-40, 40, -32, 9.0);
    fsigmoidRun2->SetLineColor(kBlue + 1);
    fsigmoidRun2->SetLineWidth(2);
    fsigmoidRun2->SetLineStyle(2);
    TFitResultPtr fitResultRun2 = gKstar_MeanpT_13TeV[1]->Fit(fsigmoidRun2, "RS");

    // Create Canvas using the same helper as dN/dy
    TCanvas *cMeanpT = new TCanvas("cMeanpT", "cMeanpT", 720, 720);
    SetCanvasStyle(cMeanpT, 0.15, 0.03, 0.03, 0.15);
    double pad1Size_pT, pad2Size_pT;
    canvas_style(cMeanpT, pad1Size_pT, pad2Size_pT);

    // ==========================================
    // TOP PANEL: <pT> vs <dNch/deta> + 13 TeV Fit
    // ==========================================
    cMeanpT->cd(1);

    gMPtKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMPtKstar[0]->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    gMPtKstar[0]->GetXaxis()->SetLimits(0, 27.5);
    gMPtKstar[0]->GetYaxis()->SetRangeUser(0.65, 1.78);
    gMPtKstar[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size_pT);
    gMPtKstar[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size_pT);
    gMPtKstar[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size_pT);
    gMPtKstar[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size_pT);

    gMPtKstar[0]->SetLineWidth(3);
    gMPtKstar[0]->SetMarkerColor(kRed + 1);
    gMPtKstar[0]->SetLineColor(kRed + 1);
    gMPtKstar[0]->Draw("APE");

    gMPtKstar[1]->SetFillStyle(0);
    gMPtKstar[1]->SetLineColor(kRed + 1);
    gMPtKstar[1]->SetLineWidth(3);
    gMPtKstar[1]->Draw("5 same");

    gMPtKstar[2]->SetFillStyle(1001);
    gMPtKstar[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gMPtKstar[2]->Draw("E2 same");

    gKstar_MeanpT_13TeV[0]->SetMarkerStyle(21);
    gKstar_MeanpT_13TeV[0]->SetMarkerColor(kBlue + 1);
    gKstar_MeanpT_13TeV[0]->SetLineColor(kBlue + 1);
    gKstar_MeanpT_13TeV[0]->SetLineWidth(3);
    gKstar_MeanpT_13TeV[0]->Draw("P same");

    gKstar_MeanpT_13TeV[1]->SetLineColor(kBlue + 1);
    gKstar_MeanpT_13TeV[1]->SetFillStyle(0);
    gKstar_MeanpT_13TeV[1]->SetLineWidth(3);
    gKstar_MeanpT_13TeV[1]->Draw("5 same");

    gKstar_MeanpT_13TeV[2]->SetFillStyle(1001);
    gKstar_MeanpT_13TeV[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gKstar_MeanpT_13TeV[2]->Draw("E2 same");

    // Legend in Top Panel including the 13 TeV Fit Function
    TLegend *leg = new TLegend(0.18, 0.72, 0.60, 0.93);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gMPtKstar[0], "Run 3 (13.6 TeV)", "lp");
    leg->AddEntry(gKstar_MeanpT_13TeV[0], "Run 2 (13 TeV)", "lp");
    leg->AddEntry(fsigmoidRun2, "Sigmoid Fit: p_{0} + #frac{p_{1}}{1 + exp(-(x - p_{2})/p_{3})}", "l");
    leg->Draw();

    TLatex latexMeanpT;
    latexMeanpT.SetTextSize(0.03 / pad1Size_pT);
    latexMeanpT.SetTextFont(22);
    latexMeanpT.SetTextFont(42);
    latexMeanpT.DrawLatexNDC(0.78, 0.87, "K*(892)^{0}");

    // ==========================================
    // BOTTOM PANEL: Ratio (13.6 TeV / 13 TeV Fit)
    // ==========================================
    cMeanpT->cd(2);

    if (!gMPtKstar[0] || !gMPtKstar[1] || !gKstar_MeanpT_13TeV[1])
    {
        cout << "Error: One of the MeanpT graphs is null!" << endl;
        return;
    }

    int nPointsRun3 = gMPtKstar[0]->GetN();
    int nPointsRun2 = gKstar_MeanpT_13TeV[1]->GetN();
    TGraphErrors *gRatio_meanpT = new TGraphErrors(nPointsRun3);
    TGraphErrors *gRatio_meanpT_sys = new TGraphErrors(nPointsRun3);

    // Compute exact 1-sigma fit confidence intervals at the 13.6 TeV x-positions
    std::vector<double> xVals(nPointsRun3), fitErrs(nPointsRun3, 0.0);
    for (int i = 0; i < nPointsRun3; i++)
    {
        double x = 0, y = 0;
        gMPtKstar[0]->GetPoint(i, x, y);
        xVals[i] = x;
    }
    if (fitResultRun2.Get())
    {
        fitResultRun2->GetConfidenceIntervals(nPointsRun3, 1, 1, xVals.data(), fitErrs.data(), 0.683, false);
    }

    for (int i = 0; i < nPointsRun3; i++)
    {
        double x = 0, yRun3 = 0;
        gMPtKstar[0]->GetPoint(i, x, yRun3);

        double yRun2Fit = fsigmoidRun2->Eval(x);
        double ratio = (yRun2Fit != 0) ? (yRun3 / yRun2Fit) : 0.0;

        // Safe error getter (works for both TGraphErrors & TGraphAsymmErrors)
        double errRun3 = gMPtKstar[1]->GetErrorY(i);

        // Find nearest point in 13 TeV (gKstar_MeanpT_13TeV[1]) for total uncertainty
        int nearestIdx = 0;
        double minDist = 1e9;
        for (int j = 0; j < nPointsRun2; j++)
        {
            double x13 = 0, y13 = 0;
            gKstar_MeanpT_13TeV[1]->GetPoint(j, x13, y13);
            double dist = fabs(x - x13);
            if (dist < minDist)
            {
                minDist = dist;
                nearestIdx = j;
            }
        }

        double x13Nearest = 0, yRun2Nearest = 0;
        gKstar_MeanpT_13TeV[1]->GetPoint(nearestIdx, x13Nearest, yRun2Nearest);
        double errRun2Nearest = gKstar_MeanpT_13TeV[1]->GetErrorY(nearestIdx);

        // 3. Fit uncertainty at x from covariance matrix
        double errFit = fitErrs[i];

        // Propagate relative uncertainties in quadrature
        double relErrRun3 = (yRun3 != 0) ? (errRun3 / yRun3) : 0.0;
        double relErrRun2 = (yRun2Nearest != 0) ? (errRun2Nearest / yRun2Nearest) : 0.0;
        double relErrFit = (yRun2Fit != 0) ? (errFit / yRun2Fit) : 0.0;
        // double errRatio = ratio * sqrt(pow(relErrRun3, 2) + pow(relErrRun2, 2) + pow(relErrFit, 2));
        double errRatio = ratio * sqrt(pow(relErrRun3, 2) + pow(relErrRun2, 2));
        double errxValue = gMPtKstar[0]->GetErrorX(i);

        gRatio_meanpT->SetPoint(i, x, ratio);
        gRatio_meanpT->SetPointError(i, 0, 0);

        gRatio_meanpT_sys->SetPoint(i, x, ratio);
        gRatio_meanpT_sys->SetPointError(i, errxValue, errRatio);
    }

    // Apply identical ratio styling as dN/dy
    SetGraphErrorStyle(gRatio_meanpT);
    gRatio_meanpT->SetMarkerStyle(21);
    gRatio_meanpT->SetMarkerColor(kRed);
    gRatio_meanpT->SetLineColor(kRed);
    gRatio_meanpT->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gRatio_meanpT->GetYaxis()->SetTitle("Run 3 / Run 2");
    gRatio_meanpT->GetXaxis()->SetLimits(0, 27.5);
    gRatio_meanpT->GetYaxis()->SetRangeUser(0.82, 1.18);
    gRatio_meanpT->GetYaxis()->SetLabelSize(0.04 / pad2Size_pT);
    gRatio_meanpT->GetXaxis()->SetLabelSize(0.04 / pad2Size_pT);
    gRatio_meanpT->GetYaxis()->SetTitleSize(0.032 / pad2Size_pT);
    gRatio_meanpT->GetXaxis()->SetTitleSize(0.035 / pad2Size_pT);
    gRatio_meanpT->GetYaxis()->SetTitleOffset(1.9 * pad2Size_pT);
    gRatio_meanpT->GetYaxis()->SetNdivisions(505);
    gRatio_meanpT->Draw("AP");
    gRatio_meanpT_sys->SetLineColor(kRed);
    gRatio_meanpT_sys->SetFillStyle(0);
    gRatio_meanpT_sys->SetLineWidth(2);
    gRatio_meanpT_sys->Draw("5 same");

    TLine *lineUnity = new TLine(0, 1, 27.5, 1);
    lineUnity->SetLineColor(kBlack);
    lineUnity->SetLineStyle(2);
    lineUnity->SetLineWidth(2);
    lineUnity->Draw("same");

    cMeanpT->SaveAs("Plots/MeanpT_KstarFit_ratio.png");

    //============================================================
    // ================ <pT> Phi fit =============================
    //============================================================

    // Fit <pT> 13 TeV Phi with sigmoid (same function, parameters, and fit option)
    TF1 *fsigmoidPhiRun2 = new TF1("fSigmoidPhi", "[0]+[1]/(1+TMath::Exp(-(x-[2])/[3]))", 2.3, 28);
    fsigmoidPhiRun2->SetParameters(-40, 40, -32, 9.0);
    fsigmoidPhiRun2->SetLineColor(kBlue + 1);
    fsigmoidPhiRun2->SetLineWidth(2);
    fsigmoidPhiRun2->SetLineStyle(2);
    TFitResultPtr fitResultPhiRun2 = gPhi_MeanpT_13TeV[1]->Fit(fsigmoidPhiRun2, "RS");

    // Create Canvas using the same helper
    TCanvas *cMeanpTPhi = new TCanvas("cMeanpTPhi", "cMeanpTPhi", 720, 720);
    SetCanvasStyle(cMeanpTPhi, 0.15, 0.03, 0.03, 0.15);
    double pad1Size_pTPhi, pad2Size_pTPhi;
    canvas_style(cMeanpTPhi, pad1Size_pTPhi, pad2Size_pTPhi);

    // ==========================================
    // TOP PANEL: Phi <pT> vs <dNch/deta> + 13 TeV Fit
    // ==========================================
    cMeanpTPhi->cd(1);

    gMPtPhi[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMPtPhi[0]->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    gMPtPhi[0]->GetXaxis()->SetLimits(0, 27.5);
    gMPtPhi[0]->GetYaxis()->SetRangeUser(0.75, 1.85);
    gMPtPhi[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size_pTPhi);
    gMPtPhi[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size_pTPhi);
    gMPtPhi[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size_pTPhi);
    gMPtPhi[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size_pTPhi);

    gMPtPhi[0]->SetLineWidth(3);
    gMPtPhi[0]->SetMarkerColor(kRed + 1);
    gMPtPhi[0]->SetLineColor(kRed + 1);
    gMPtPhi[0]->Draw("APE");

    gMPtPhi[1]->SetFillStyle(0);
    gMPtPhi[1]->SetLineColor(kRed + 1);
    gMPtPhi[1]->SetLineWidth(3);
    gMPtPhi[1]->Draw("5 same");

    gMPtPhi[2]->SetFillStyle(1001);
    gMPtPhi[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gMPtPhi[2]->Draw("E2 same");

    gPhi_MeanpT_13TeV[0]->SetMarkerStyle(21);
    gPhi_MeanpT_13TeV[0]->SetMarkerColor(kBlue + 1);
    gPhi_MeanpT_13TeV[0]->SetLineColor(kBlue + 1);
    gPhi_MeanpT_13TeV[0]->SetLineWidth(3);
    gPhi_MeanpT_13TeV[0]->Draw("P same");

    gPhi_MeanpT_13TeV[1]->SetLineColor(kBlue + 1);
    gPhi_MeanpT_13TeV[1]->SetFillStyle(0);
    gPhi_MeanpT_13TeV[1]->SetLineWidth(3);
    gPhi_MeanpT_13TeV[1]->Draw("5 same");

    gPhi_MeanpT_13TeV[2]->SetFillStyle(1001);
    gPhi_MeanpT_13TeV[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gPhi_MeanpT_13TeV[2]->Draw("E2 same");

    // Legend in Top Panel including the 13 TeV Fit Function
    TLegend *legPhi = new TLegend(0.18, 0.72, 0.60, 0.93);
    legPhi->SetTextSize(0.03);
    legPhi->SetBorderSize(0);
    legPhi->SetFillStyle(0);
    legPhi->AddEntry(gMPtPhi[0], "#phi Run 3 (13.6 TeV)", "lp");
    legPhi->AddEntry(gPhi_MeanpT_13TeV[0], "#phi Run 2 (13 TeV)", "lp");
    legPhi->AddEntry(fsigmoidPhiRun2, "Sigmoid Fit: p_{0} + #frac{p_{1}}{1 + exp(-(x - p_{2})/p_{3})}", "l");
    legPhi->Draw();

    latexMeanpT.DrawLatexNDC(0.78, 0.87, "#phi(1020)");

    // ==========================================
    // BOTTOM PANEL: Ratio (13.6 TeV / 13 TeV Fit)
    // ==========================================
    cMeanpTPhi->cd(2);

    if (!gMPtPhi[0] || !gMPtPhi[1] || !gPhi_MeanpT_13TeV[1])
    {
        cout << "Error: One of the Phi MeanpT graphs is null!" << endl;
        return;
    }

    int nPointsPhiRun3 = gMPtPhi[0]->GetN();
    int nPointsPhiRun2 = gPhi_MeanpT_13TeV[1]->GetN();
    TGraphErrors *gRatio_meanpTPhi = new TGraphErrors(nPointsPhiRun3);
    TGraphErrors *gRatio_meanpTPhi_sys = new TGraphErrors(nPointsPhiRun3);

    // Compute exact 1-sigma fit confidence intervals at the 13.6 TeV x-positions
    std::vector<double> xValsPhi(nPointsPhiRun3), fitErrsPhi(nPointsPhiRun3, 0.0);
    for (int i = 0; i < nPointsPhiRun3; i++)
    {
        double x = 0, y = 0;
        gMPtPhi[0]->GetPoint(i, x, y);
        xValsPhi[i] = x;
    }
    if (fitResultPhiRun2.Get())
    {
        fitResultPhiRun2->GetConfidenceIntervals(nPointsPhiRun3, 1, 1, xValsPhi.data(), fitErrsPhi.data(), 0.683, false);
    }

    for (int i = 0; i < nPointsPhiRun3; i++)
    {
        double x = 0, yRun3 = 0;
        gMPtPhi[0]->GetPoint(i, x, yRun3);

        double yRun2Fit = fsigmoidPhiRun2->Eval(x);
        double ratio = (yRun2Fit != 0) ? (yRun3 / yRun2Fit) : 0.0;

        // Total systematic uncertainty of 13.6 TeV point (gMPtPhi[1])
        double errRun3 = gMPtPhi[1]->GetErrorY(i);

        // Find nearest point in 13 TeV (gPhi_MeanpT_13TeV[1]) for total uncertainty
        int nearestIdx = 0;
        double minDist = 1e9;
        for (int j = 0; j < nPointsPhiRun2; j++)
        {
            double x13 = 0, y13 = 0;
            gPhi_MeanpT_13TeV[1]->GetPoint(j, x13, y13);
            double dist = fabs(x - x13);
            if (dist < minDist)
            {
                minDist = dist;
                nearestIdx = j;
            }
        }

        double x13Nearest = 0, yRun2Nearest = 0;
        gPhi_MeanpT_13TeV[1]->GetPoint(nearestIdx, x13Nearest, yRun2Nearest);
        double errRun2Nearest = gPhi_MeanpT_13TeV[1]->GetErrorY(nearestIdx);

        // Fit uncertainty at x from covariance matrix
        double errFit = fitErrsPhi[i];

        // Propagate relative uncertainties in quadrature
        double relErrRun3 = (yRun3 != 0) ? (errRun3 / yRun3) : 0.0;
        double relErrRun2 = (yRun2Nearest != 0) ? (errRun2Nearest / yRun2Nearest) : 0.0;
        double relErrFit = (yRun2Fit != 0) ? (errFit / yRun2Fit) : 0.0;

        // Uncomment the line below if you want to include the fit uncertainty in the systematic box
        // double errRatio = ratio * sqrt(pow(relErrRun3, 2) + pow(relErrRun2, 2) + pow(relErrFit, 2));
        double errRatio = ratio * sqrt(pow(relErrRun3, 2) + pow(relErrRun2, 2));
        double errxValue = gMPtPhi[0]->GetErrorX(i);

        gRatio_meanpTPhi->SetPoint(i, x, ratio);
        gRatio_meanpTPhi->SetPointError(i, 0, 0);

        gRatio_meanpTPhi_sys->SetPoint(i, x, ratio);
        gRatio_meanpTPhi_sys->SetPointError(i, errxValue, errRatio);
    }

    // Apply identical ratio styling
    SetGraphErrorStyle(gRatio_meanpTPhi);
    gRatio_meanpTPhi->SetMarkerStyle(21);
    gRatio_meanpTPhi->SetMarkerColor(kRed);
    gRatio_meanpTPhi->SetLineColor(kRed);
    gRatio_meanpTPhi->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gRatio_meanpTPhi->GetYaxis()->SetTitle("Run 3 / Run 2");
    gRatio_meanpTPhi->GetXaxis()->SetLimits(0, 27.5);
    gRatio_meanpTPhi->GetYaxis()->SetRangeUser(0.82, 1.18);
    gRatio_meanpTPhi->GetYaxis()->SetLabelSize(0.04 / pad2Size_pTPhi);
    gRatio_meanpTPhi->GetXaxis()->SetLabelSize(0.04 / pad2Size_pTPhi);
    gRatio_meanpTPhi->GetYaxis()->SetTitleSize(0.032 / pad2Size_pTPhi);
    gRatio_meanpTPhi->GetXaxis()->SetTitleSize(0.035 / pad2Size_pTPhi);
    gRatio_meanpTPhi->GetYaxis()->SetTitleOffset(1.9 * pad2Size_pTPhi);
    gRatio_meanpTPhi->GetYaxis()->SetNdivisions(505);
    gRatio_meanpTPhi->Draw("AP");

    gRatio_meanpTPhi_sys->SetLineColor(kRed);
    gRatio_meanpTPhi_sys->SetFillStyle(0);
    gRatio_meanpTPhi_sys->SetLineWidth(2);
    gRatio_meanpTPhi_sys->Draw("5 same");

    TLine *lineUnityPhi = new TLine(0, 1, 27.5, 1);
    lineUnityPhi->SetLineColor(kBlack);
    lineUnityPhi->SetLineStyle(2);
    lineUnityPhi->SetLineWidth(2);
    lineUnityPhi->Draw("same");

    cMeanpTPhi->SaveAs("Plots/MeanpTPhiFit_ratio.png");
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

TGraphErrors *GetGraph(TFile *f, const string &name)
{
    TGraphErrors *graph = (TGraphErrors *)f->Get(name.c_str());

    if (!graph || graph == nullptr)
    {
        cout << "Error: graph " << name
             << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetGraphErrorStyle(graph);
    graph->SetTitle(0);
    return graph;
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
    pad2->SetBottomMargin(0.35);
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