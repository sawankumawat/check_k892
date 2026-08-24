#include <iostream>
#include "src/style.h"

TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const std::string &name);

Double_t voigt(Double_t *x, Double_t *par)
{
    return (par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]));
}

Double_t polynomial2(Double_t *x, Double_t *par)
{
    double z = x[0] - 1.0019; // Shift the x-axis to center around the phi meson mass
    // double poly2 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
    double poly2 = par[0] + par[1] * z + par[2] * z * z;
    return (poly2);
}
Double_t voigtpol2(Double_t *x, Double_t *par)
{
    double vgt = par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3], 4);
    double poly2 = polynomial2(x, &par[4]);
    return (vgt + poly2);
}

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

TGraph *smoothGraph(TGraph *g, int n = 3)
{
    int N = g->GetN();
    TGraph *gs = new TGraph();

    for (int i = 0; i < N; i++)
    {
        double x, y;
        g->GetPoint(i, x, y);

        double sum = 0;
        int count = 0;

        for (int j = TMath::Max(0, i - n);
             j <= TMath::Min(N - 1, i + n); j++)
        {

            double xj, yj;
            g->GetPoint(j, xj, yj);
            sum += yj;
            count++;
        }

        gs->SetPoint(i, x, sum / count);
    }

    return gs;
}

void SignificanceTimesPurity()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString suffix = "25_aiam";
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/PhiInvMass";

    ////=====New===========
    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/LHC25/AnalysisResults25_aiam.root"); // 2025 ai+am

    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassPhiPhiRefitted");

    double intervals[] = {0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023};
    double totalIntervals = sizeof(intervals) / sizeof(intervals[0]);

    int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowChi2 = hUnlike->GetAxis(3)->FindBin(0.0 + 0.00001);
    int highChi2 = hUnlike->GetAxis(3)->FindBin(10.0 - 0.00001);

    int lowProbability = hUnlike->GetAxis(4)->FindBin(0.0 + 0.00001);
    int highProbability = hUnlike->GetAxis(4)->FindBin(1.0 - 0.00001);

    TGraph *gSigTimesPurity = new TGraph();

    for (int intervalBin = 0; intervalBin < totalIntervals; ++intervalBin)
    {
        double deltaMValue = intervals[intervalBin];

        int lowDeltaM = hUnlike->GetAxis(2)->FindBin(0.0 + 0.00001);
        int highDeltaM = hUnlike->GetAxis(2)->FindBin(deltaMValue - 0.00001);

        hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
        hUnlike->GetAxis(2)->SetRange(lowDeltaM, highDeltaM);

        TH1D *hInvMass = hUnlike->Projection(0, "E");
        hInvMass->SetName(Form("hInvMass_deltaM_%d", intervalBin));
        SetHistoQA(hInvMass);
        hInvMass->Rebin(8);
        hInvMass->GetXaxis()->SetRangeUser(2.4, 2.95);
        hInvMass->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
        hInvMass->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hInvMass->GetBinWidth(1) * 1000));

        TH1D *hBkg = (TH1D *)hInvMass->Clone(Form("hBkg_deltaM_%d", intervalBin));
        hBkg->Reset();

        for (int i = 1; i <= hInvMass->GetNbinsX(); i++)
        {

            double x = hInvMass->GetBinCenter(i);

            // Exclude signal region
            if (x >= 2.65 && x <= 2.75)
                continue;

            hBkg->SetBinContent(i, hInvMass->GetBinContent(i));
            hBkg->SetBinError(i, hInvMass->GetBinError(i));
        }

        TCanvas *cOnlyBkg = new TCanvas(Form("cOnlyBkg_deltaM_%d", intervalBin), "Background Only", 720, 720);
        SetCanvasStyle(cOnlyBkg, 0.15, 0.03, 0.05, 0.15);
        hBkg->Draw("pe");

        TF1 *fitBkg = new TF1("fitBkg", expPol3, 2.42, 2.95, 4);
        fitBkg->SetParNames("p0", "p1", "p2", "p3");
        // fitBkg->SetParameters(-300.0, 70, 46, -9);
        fitBkg->SetParameters(1, 1, 1, 1);
        if (intervalBin == 4 || intervalBin == 5 || intervalBin == 6)
        {
            fitBkg->SetParameters(-200.0, 40, 23, -5);
        }
        else
        {
            fitBkg->SetParameters(1, 1, 1, 1);
        }
        hBkg->Fit(fitBkg, "R");
        cOnlyBkg->SaveAs(savepath + Form("/PhiPhivsDeltaMFits/Bkg_deltaM_%d_%s.png", intervalBin, suffix.Data()));

        TCanvas *cInvMass = new TCanvas(Form("cInvMass_deltaM_%d", intervalBin), "Invariant Mass", 720, 720);
        SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
        hInvMass->SetMaximum(hInvMass->GetMaximum() * 1.2);
        hInvMass->Draw("pe");

        // //======================Fit function (BW)=========================
        TF1 *fitFunc = new TF1(Form("fitFunc_deltaM_%d", intervalBin), BWExpol, 2.5, 2.9, 7);
        fitFunc->SetParNames("Yield", "Mass", "Width", "p0", "p1", "p2", "p3");

        fitFunc->SetParameter(0, 100);  // Yield
        fitFunc->SetParameter(1, 2.7);  // Mass
        fitFunc->SetParameter(2, 0.03); // Width

        fitFunc->SetParameter(3, fitBkg->GetParameter(0)); // p0
        fitFunc->SetParameter(4, fitBkg->GetParameter(1)); // p1
        fitFunc->SetParameter(5, fitBkg->GetParameter(2)); // p2
        fitFunc->SetParameter(6, fitBkg->GetParameter(3)); // p3

        fitFunc->SetParLimits(0, 0, 1e4);
        fitFunc->SetParLimits(1, 2.65, 2.75);
        fitFunc->SetParLimits(2, 0.01, 0.05);

        // hInvMass->Fit(fitFunc, "REBMS");

        hInvMass->Fit(fitFunc, "QERN");
        TFitResultPtr fitResult = hInvMass->Fit(fitFunc, "ERSN");

        TF1 *fitBkgFinal = new TF1(Form("fitBkgFinal_deltaM_%d", intervalBin), expPol3, 2.42, 2.95, 4);
        fitBkgFinal->SetParameters(fitFunc->GetParameter(3), fitFunc->GetParameter(4), fitFunc->GetParameter(5), fitFunc->GetParameter(6));
        fitBkgFinal->SetLineColor(kBlue);
        fitBkgFinal->SetLineStyle(2);
        fitBkgFinal->Draw("same");

        TF1 *fitBW = new TF1(Form("fitBW_deltaM_%d", intervalBin), breitWigner, 2.42, 2.95, 3);
        fitBW->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2));
        fitBW->SetLineColor(kRed);
        fitBW->SetLineStyle(2);
        fitBW->Draw("same");
        cInvMass->SaveAs(savepath + Form("/PhiPhivsDeltaMFits/DoublePhiFit_deltaM_%d_%s.png", intervalBin, suffix.Data()));

        double IntegralLow = fitFunc->GetParameter(1) - 1.0 * fitFunc->GetParameter(2);
        double IntegralHigh = fitFunc->GetParameter(1) + 1.0 * fitFunc->GetParameter(2);

        int binLow = hInvMass->GetXaxis()->FindBin(IntegralLow);
        int binHigh = hInvMass->GetXaxis()->FindBin(IntegralHigh);

        double signalCounts = fitBW->Integral(IntegralLow, IntegralHigh) / hInvMass->GetBinWidth(1);
        double BkgCounts = fitBkgFinal->Integral(IntegralLow, IntegralHigh) / hInvMass->GetBinWidth(1);
        double SigBkg = hInvMass->Integral(binLow, binHigh);

        // double purity = (SigBkg - BkgCounts) / SigBkg;
        // double significance = signalCounts / sqrt(SigBkg);

        // double purity = (signalCounts) / (signalCounts + BkgCounts);
        // double significance = signalCounts / sqrt(SigBkg);

        double purity = signalCounts / SigBkg;
        double significance = signalCounts / sqrt(SigBkg);

        double SignificanceTimesPurity = significance * purity;
        gSigTimesPurity->SetPoint(intervalBin, deltaMValue, SignificanceTimesPurity * 100);
    }

    TCanvas *cSigTimesPurity = new TCanvas("cSigTimesPurity", "Significance times Purity", 720, 720);
    SetCanvasStyle(cSigTimesPurity, 0.15, 0.03, 0.09, 0.15);
    TGraph *gSmooth = smoothGraph(gSigTimesPurity, 2);
    SetGraphStyle(gSigTimesPurity);
    gSigTimesPurity->SetMarkerStyle(20);
    gSigTimesPurity->SetMarkerSize(0.8);
    gSigTimesPurity->GetXaxis()->SetTitle("#Delta#it{M}_{#phi} (GeV/#it{c}^{2})");
    gSigTimesPurity->GetYaxis()->SetTitle("S/(#sqrt{S+B}) #times S/(S+B) (%)");
    gSigTimesPurity->Draw("APL");
    gSmooth->SetLineColor(kRed);
    gSmooth->SetLineWidth(2);
    gSmooth->Draw("L SAME");
    TLatex *latex3 = new TLatex();
    latex3->SetNDC();
    latex3->SetTextFont(22);
    latex3->SetTextSize(0.05);
    latex3->DrawLatex(0.19, 0.93, Form("Significance x purity vs #Delta#it{M_{#phi}}"));
    latex3->SetTextFont(42);
    latex3->SetTextSize(0.035);
    latex3->DrawLatex(0.5, 0.85, Form("#it{p}_{T}^{#phi#phi} > 9 GeV/#it{c}"));
    latex3->DrawLatex(0.5, 0.79, Form("|y^{#phi#phi}| < 0.8"));
    latex3->DrawLatex(0.5, 0.73, Form("2.5 < #it{M}_{#phi#phi} < 2.9 GeV/#it{c}^{2}"));
    // cSigTimesPurity->SaveAs(savepath + "/SignificanceTimesPurity_" + suffix + ".png");

    TFile *fOutput = new TFile(savepath + "/SignificanceTimesPurity_" + suffix + ".root", "recreate");
    gSigTimesPurity->Write("gSigTimesPurity");
    gSmooth->Write("gSmooth");
    fOutput->Close();
}

// //====================From histogram SEMassUnlike_AllVars===============================
// THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassUnlike_AllVars");

// // It has 10 axes:
// // 0: Invariant Mass
// // 1: Pt
// // 2: DeltaR Phi
// // 3: Min DeltaR Kaon
// // 4: Phi1 Mass
// // 5: Phi2 Mass
// // 6: DeltaM
// // 7: DeltaM Normalized
// // 8: Pt Correlation
// // 9: PhiPhi Vector Size

// TH1D *hPt = hUnlike->Projection(1, "E");
// TH1D *hDeltaR = hUnlike->Projection(2, "E");
// TH1D *hMinDeltaRKaon = hUnlike->Projection(3, "E");
// TH1D *hPhi1Mass = hUnlike->Projection(4, "E");
// TH1D *hPhi2Mass = hUnlike->Projection(5, "E");
// TH1D *hDeltaM = hUnlike->Projection(6, "E");
// TH1D *hDeltaMNormalized = hUnlike->Projection(7, "E");
// TH1D *hPtCorrelation = hUnlike->Projection(8, "E");
// TH1D *hPhiPhiVectorSize = hUnlike->Projection(9, "E");

// int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
// int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

// int lowDeltaR = hUnlike->GetAxis(2)->FindBin(0.0 + 0.0001);
// int highDeltaR = hUnlike->GetAxis(2)->FindBin(0.55 - 0.0001);

// int lowMinDeltaRKaon = hUnlike->GetAxis(3)->FindBin(0.00 + 0.0001);
// int highMinDeltaRKaon = hUnlike->GetAxis(3)->FindBin(0.15 - 0.0001);

// int lowPhi1Mass = hUnlike->GetAxis(4)->FindBin(1.005 + 0.0001);
// int highPhi1Mass = hUnlike->GetAxis(4)->FindBin(1.035 - 0.0001);

// int lowPhi2Mass = hUnlike->GetAxis(5)->FindBin(1.005 + 0.0001);
// int highPhi2Mass = hUnlike->GetAxis(5)->FindBin(1.035 - 0.0001);

// int lowDeltaM = hUnlike->GetAxis(6)->FindBin(0.0 + 0.00001);
// int highDeltaM = hUnlike->GetAxis(6)->FindBin(0.005 - 0.00001);

// int lowDeltaMNormalized = hUnlike->GetAxis(7)->FindBin(0.0 + 0.00001);
// int highDeltaMNormalized = hUnlike->GetAxis(7)->FindBin(2.5 - 0.00001);

// int lowPtCorrelation = hUnlike->GetAxis(8)->FindBin(0.5 + 0.0001);
// int highPtCorrelation = hUnlike->GetAxis(8)->FindBin(4.5 - 0.0001);

// int lowPhiPhiVectorSize = hUnlike->GetAxis(9)->FindBin(1.0 + 0.0001);
// int highPhiPhiVectorSize = hUnlike->GetAxis(9)->FindBin(3.0 - 0.0001);

// hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
// // hUnlike->GetAxis(2)->SetRange(lowDeltaR, highDeltaR);
// // hUnlike->GetAxis(3)->SetRange(lowMinDeltaRKaon, highMinDeltaRKaon);
// // hUnlike->GetAxis(4)->SetRange(lowPhi1Mass, highPhi1Mass);
// // hUnlike->GetAxis(5)->SetRange(lowPhi2Mass, highPhi2Mass);
// hUnlike->GetAxis(6)->SetRange(lowDeltaM, highDeltaM);
// // hUnlike->GetAxis(7)->SetRange(lowDeltaMNormalized, highDeltaMNormalized);
// // hUnlike->GetAxis(8)->SetRange(lowPtCorrelation, highPtCorrelation);
// // hUnlike->GetAxis(9)->SetRange(lowPhiPhiVectorSize, highPhiPhiVectorSize);

// //====================From histogram SEMassDoublePhi===============================
// //It has 9 axes:
// 0: Invariant Mass
// 1: Pt
// 2: Phi Pt asymmetry
// 3: Rapidity
// 4: Phi1 Mass
// 5: Phi2 Mass
// 6: DeltaM
// 7: nKaon TOF
// 8: Combined PID 4kaon

//==============End of the main code==================

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
        std::cout << "Error: histo " << name
                  << " not found in file " << f->GetName() << std::endl;
        return nullptr;
    }

    return histo;
}
