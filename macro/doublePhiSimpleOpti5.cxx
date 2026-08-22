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

void doublePhiSimpleOpti5()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/PhiInvMass";
    ////======Pair=========
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/Pair/processopti5/AnalysisResults.root");

    ////=====New===========
    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults_PhiMassInPair.root");
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/LHC25/AnalysisResults_LHC25_PID2003.root"); //2025 data

    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassDoublePhi");

    // // TH1D *hPt = hUnlike->Projection(1, "E");
    // // TH1D *hPhiPtAsymmetry = hUnlike->Projection(2, "E");
    // // TH1D *hRapidity = hUnlike->Projection(3, "E");
    // // TH1D *hPhi1Mass = hUnlike->Projection(4, "E");
    // // TH1D *hPhi2Mass = hUnlike->Projection(5, "E");
    // // TH1D *hDeltaM = hUnlike->Projection(6, "E");
    // // TH1D *hnKaonTOF = hUnlike->Projection(7, "E");
    // // TH1D *hCombinedPID4Kaon = hUnlike->Projection(8, "E");

    // // vector<TH1D *> histos = {hPt, hPhiPtAsymmetry, hRapidity, hPhi1Mass, hPhi2Mass, hDeltaM, hnKaonTOF, hCombinedPID4Kaon};
    // // //Draw all histograms using loop
    // // for (auto histo : histos)
    // // {
    // //     SetHistoQA(histo);
    // //     TCanvas *c = new TCanvas(Form("c%s", histo->GetName()), histo->GetTitle(), 720, 720);
    // //     SetCanvasStyle(c, 0.15, 0.03, 0.05, 0.15);
    // //     histo->Draw("pe");
    // // }

    // // double deltaMLowValue = 0.002;
    // // double deltaMHighValue = 0.03;
    // // int totalIntervals = (deltaMHighValue - deltaMLowValue) / 0.001; // 29 intervals of 0.001 width

    // double intervals[] = {0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023};
    // double totalIntervals = sizeof(intervals) / sizeof(intervals[0]);

    // int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    // int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    // int lowPhiPtAsymmetry = hUnlike->GetAxis(2)->FindBin(0.0 + 0.0001);
    // int highPhiPtAsymmetry = hUnlike->GetAxis(2)->FindBin(1.0 - 0.0001);

    // int lowRapidity = hUnlike->GetAxis(3)->FindBin(-0.8 + 0.0001);
    // int highRapidity = hUnlike->GetAxis(3)->FindBin(0.8 - 0.0001);

    // int lowPhi1Mass = hUnlike->GetAxis(4)->FindBin(1.005 + 0.0001);
    // int highPhi1Mass = hUnlike->GetAxis(4)->FindBin(1.035 - 0.0001);

    // int lowPhi2Mass = hUnlike->GetAxis(5)->FindBin(1.005 + 0.0001);
    // int highPhi2Mass = hUnlike->GetAxis(5)->FindBin(1.035 - 0.0001);

    // int lownKaonTOF = hUnlike->GetAxis(7)->FindBin(0.0 + 0.0001);
    // int highnKaonTOF = hUnlike->GetAxis(7)->FindBin(4.0 - 0.0001);

    // int lowCombinedPID4Kaon = hUnlike->GetAxis(8)->FindBin(0.0 + 0.0001);
    // int highCombinedPID4Kaon = hUnlike->GetAxis(8)->FindBin(4.0 - 0.0001);

    // TGraph *gSigTimesPurity = new TGraph();

    // for (int intervalBin = 0; intervalBin < totalIntervals; ++intervalBin)
    // {
    //     // double deltaMValue = deltaMLowValue + intervalBin * 0.001;
    //     double deltaMValue = intervals[intervalBin];

    //     int lowDeltaM = hUnlike->GetAxis(6)->FindBin(0.0 + 0.00001);
    //     int highDeltaM = hUnlike->GetAxis(6)->FindBin(deltaMValue - 0.00001);

    //     hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    //     hUnlike->GetAxis(6)->SetRange(lowDeltaM, highDeltaM);

    //     TH1D *hInvMass = hUnlike->Projection(0, "E");
    //     hInvMass->SetName(Form("hInvMass_deltaM_%d", intervalBin));
    //     SetHistoQA(hInvMass);
    //     hInvMass->Rebin(8);
    //     hInvMass->GetXaxis()->SetRangeUser(2.5, 2.9);
    //     hInvMass->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    //     hInvMass->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hInvMass->GetBinWidth(1) * 1000));

    //     TH1D *hBkg = (TH1D *)hInvMass->Clone(Form("hBkg_deltaM_%d", intervalBin));
    //     hBkg->Reset();

    //     for (int i = 1; i <= hInvMass->GetNbinsX(); i++)
    //     {

    //         double x = hInvMass->GetBinCenter(i);

    //         // Exclude signal region
    //         if (x >= 2.65 && x <= 2.75)
    //             continue;

    //         hBkg->SetBinContent(i, hInvMass->GetBinContent(i));
    //         hBkg->SetBinError(i, hInvMass->GetBinError(i));
    //     }

    //     TCanvas *cOnlyBkg = new TCanvas(Form("cOnlyBkg_deltaM_%d", intervalBin), "Background Only", 720, 720);
    //     SetCanvasStyle(cOnlyBkg, 0.15, 0.03, 0.05, 0.15);
    //     hBkg->Draw("pe");

    //     TF1 *fitBkg = new TF1("fitBkg", expPol3, 2.5, 2.9, 4);
    //     fitBkg->SetParNames("p0", "p1", "p2", "p3");
    //     // fitBkg->SetParameters(-300.0, 70, 46, -9);
    //     fitBkg->SetParameters(1, 1, 1, 1);
    //     if (intervalBin == 4 || intervalBin == 5 || intervalBin == 6)
    //     {
    //         fitBkg->SetParameters(-200.0, 40, 23, -5);
    //     }
    //     else
    //     {
    //         fitBkg->SetParameters(1, 1, 1, 1);
    //     }
    //     hBkg->Fit(fitBkg, "R");
    //     cOnlyBkg->SaveAs(savepath + Form("/PhiPhivsDeltaMFits/BackgroundOnly_deltaM_%d.png", intervalBin));

    //     TCanvas *cInvMass = new TCanvas(Form("cInvMass_deltaM_%d", intervalBin), "Invariant Mass", 720, 720);
    //     SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
    //     hInvMass->SetMaximum(hInvMass->GetMaximum() * 1.2);
    //     hInvMass->Draw("pe");

    //     // //======================Fit function (BW)=========================
    //     TF1 *fitFunc = new TF1(Form("fitFunc_deltaM_%d", intervalBin), BWExpol, 2.5, 2.9, 7);
    //     fitFunc->SetParNames("Yield", "Mass", "Width", "p0", "p1", "p2", "p3");

    //     fitFunc->SetParameter(0, 100);  // Yield
    //     fitFunc->SetParameter(1, 2.7);  // Mass
    //     fitFunc->SetParameter(2, 0.03); // Width

    //     fitFunc->SetParameter(3, fitBkg->GetParameter(0)); // p0
    //     fitFunc->SetParameter(4, fitBkg->GetParameter(1)); // p1
    //     fitFunc->SetParameter(5, fitBkg->GetParameter(2)); // p2
    //     fitFunc->SetParameter(6, fitBkg->GetParameter(3)); // p3

    //     fitFunc->SetParLimits(0, 0, 1e4);
    //     fitFunc->SetParLimits(1, 2.65, 2.75);
    //     fitFunc->SetParLimits(2, 0.01, 0.05);

    //     // hInvMass->Fit(fitFunc, "REBMS");

    //     hInvMass->Fit(fitFunc, "Q0ERN");
    //     TFitResultPtr fitResult = hInvMass->Fit(fitFunc, "ERSN");

    //     TF1 *fitBkgFinal = new TF1(Form("fitBkgFinal_deltaM_%d", intervalBin), expPol3, 2.5, 2.9, 4);
    //     fitBkgFinal->SetParameters(fitFunc->GetParameter(3), fitFunc->GetParameter(4), fitFunc->GetParameter(5), fitFunc->GetParameter(6));
    //     fitBkgFinal->SetLineColor(kBlue);
    //     fitBkgFinal->SetLineStyle(2);
    //     fitBkgFinal->Draw("same");

    //     TF1 *fitBW = new TF1(Form("fitBW_deltaM_%d", intervalBin), breitWigner, 2.5, 2.9, 3);
    //     fitBW->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2));
    //     fitBW->SetLineColor(kRed);
    //     fitBW->SetLineStyle(2);
    //     fitBW->Draw("same");
    //     cInvMass->SaveAs(savepath + Form("/PhiPhivsDeltaMFits/DoublePhiFit_deltaM_%d.png", intervalBin));

    //     double IntegralLow = fitFunc->GetParameter(1) - 1.0 * fitFunc->GetParameter(2);
    //     double IntegralHigh = fitFunc->GetParameter(1) + 1.0 * fitFunc->GetParameter(2);

    //     int binLow = hInvMass->GetXaxis()->FindBin(IntegralLow);
    //     int binHigh = hInvMass->GetXaxis()->FindBin(IntegralHigh);

    //     double signalCounts = fitBW->Integral(IntegralLow, IntegralHigh) / hInvMass->GetBinWidth(1);
    //     double BkgCounts = fitBkgFinal->Integral(IntegralLow, IntegralHigh) / hInvMass->GetBinWidth(1);
    //     double SigBkg = hInvMass->Integral(binLow, binHigh);

    //     // double purity = (SigBkg - BkgCounts) / SigBkg;
    //     // double significance = signalCounts / sqrt(SigBkg);

    //     // double purity = (signalCounts) / (signalCounts + BkgCounts);
    //     // double significance = signalCounts / sqrt(SigBkg);

    //     double purity = signalCounts / SigBkg;
    //     double significance = signalCounts / sqrt(SigBkg);

    //     double SignificanceTimesPurity = significance * purity;
    //     gSigTimesPurity->SetPoint(intervalBin, deltaMValue, SignificanceTimesPurity * 100);
    // }

    // TCanvas *cSigTimesPurity = new TCanvas("cSigTimesPurity", "Significance times Purity", 720, 720);
    // SetCanvasStyle(cSigTimesPurity, 0.15, 0.03, 0.09, 0.15);
    // TGraph *gSmooth = smoothGraph(gSigTimesPurity, 2);
    // SetGraphStyle(gSigTimesPurity);
    // gSigTimesPurity->SetMarkerStyle(20);
    // gSigTimesPurity->SetMarkerSize(0.8);
    // gSigTimesPurity->GetXaxis()->SetTitle("#Delta#it{M}_{#phi} (GeV/#it{c}^{2})");
    // gSigTimesPurity->GetYaxis()->SetTitle("S/(#sqrt{S+B}) #times S/(S+B) (%)");
    // gSigTimesPurity->Draw("APL");
    // gSmooth->SetLineColor(kRed);
    // gSmooth->SetLineWidth(2);
    // gSmooth->Draw("L SAME");
    // TLatex *latex3 = new TLatex();
    // latex3->SetNDC();
    // latex3->SetTextFont(22);
    // latex3->SetTextSize(0.05);
    // latex3->DrawLatex(0.19, 0.93, Form("Significance x purity vs #Delta#it{M_{#phi}}"));
    // latex3->SetTextFont(42);
    // latex3->SetTextSize(0.035);
    // latex3->DrawLatex(0.5, 0.85, Form("#it{p}_{T}^{#phi#phi} > 9 GeV/#it{c}"));
    // latex3->DrawLatex(0.5, 0.79, Form("|y^{#phi#phi}| < 0.8"));
    // latex3->DrawLatex(0.5, 0.73, Form("2.5 < #it{M}_{#phi#phi} < 2.9 GeV/#it{c}^{2}"));
    // cSigTimesPurity->SaveAs(savepath + "/SignificanceTimesPurity.png");

    // //=========================================
    // //=======Phi Inv mass fit (1D)=============
    // //=========================================

    TH3F *hPhiMassVsPt = GetHisto<TH3F>(fInput, "doublephimeson/hPhiMass");
    TH2F *hPhiPhiMass = (TH2F *)hPhiMassVsPt->Project3D("yx");
    TH2F *hPhiMassVsPt2D = GetHisto<TH2F>(fInput, "doublephimeson/hPhiMassVsPt");

    TCanvas *cPhiVsPt = new TCanvas("cPhiVsPt", "Phi Mass vs Pt", 1080, 720);
    cPhiVsPt->Divide(4, 3);
    double pTBins[14] = {0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(22);
    latex->SetTextSize(0.06);
    vector<pair<double, double>> fitMass, fitResolution, fitWidth;
    TGraphErrors *gPurity = new TGraphErrors();
    for (int ibin = 0; ibin < 13; ibin++)
    {
        // int lowpTBin = hPhiMassVsPt->GetZaxis()->FindBin(pTBins[ibin] + 0.0001);
        // int highpTBin = hPhiMassVsPt->GetZaxis()->FindBin(pTBins[ibin + 1] - 0.0001);
        // TH1D *hPhi1 = hPhiMassVsPt->ProjectionY(Form("hPhi1_pT_%.1f_%.1f", pTBins[ibin], pTBins[ibin + 1]), lowpTBin, highpTBin, -1, -1, "E");

        int lowpTBin = hPhiMassVsPt2D->GetYaxis()->FindBin(pTBins[ibin] + 0.0001);
        int highpTBin = hPhiMassVsPt2D->GetYaxis()->FindBin(pTBins[ibin + 1] - 0.0001);
        TH1D *hPhi1 = hPhiMassVsPt2D->ProjectionX(Form("hPhi1_pT_%.1f_%.1f", pTBins[ibin], pTBins[ibin + 1]), lowpTBin, highpTBin, "E");

        SetHistoQA(hPhi1);
        cPhiVsPt->cd(ibin + 1);
        hPhi1->SetMarkerSize(0.8);
        hPhi1->Draw("pe");

        float fitRangeLow = 1.001;
        float fitRangeHigh = 1.039;

        //******Fitting for Phi*********************
        TF1 *fitFcn = new TF1("fitfunc", voigtpol2, fitRangeLow, fitRangeHigh, 7);       // sig+bkg fit function
        TF1 *fitFcnBkg = new TF1("fitfunc1", polynomial2, fitRangeLow, fitRangeHigh, 3); // only residualbkg
        TF1 *fitFcnSig = new TF1("fitFcnSig", voigt, fitRangeLow, fitRangeHigh, 4);      // only signal

        // for voigtian distribution
        fitFcn->SetParameter(0, 5000);          // yield
        fitFcn->SetParLimits(0, 0, 1e6);        // yield
        fitFcn->SetParameter(1, 1.0198);         // mass peak
        fitFcn->SetParLimits(1, 1.017, 1.026); // mass peak
        fitFcn->SetParameter(2, 0.0012);        //  Gaussian width (Detector resolution)
        fitFcn->SetParLimits(2, 0.0008, 0.008);  // Gaussian width.
        // fitFcn->SetParameter(3, 0.0042);   //lorentzian width (Resonance width)
        fitFcn->FixParameter(3, 0.0042); // lorentzian width

        fitFcn->SetParameter(4, 4e4);      // Pol2 p0
        fitFcn->SetParLimits(4, 1e2, 3e6);   // Pol2 p0
        fitFcn->SetParameter(5, 4.1e5);      // Pol2 p1
        fitFcn->SetParLimits(5, 1e2, 3e7);   // Pol2 p1
        fitFcn->SetParameter(6, -5.0e6);     // Pol2 p2
        fitFcn->SetParLimits(6, -1e9, -1e2); // Pol2 p2

        hPhi1->Fit("fitfunc", "REBMS");
        TVirtualFitter::SetDefaultFitter("Minuit2");
        TVirtualFitter::SetMaxIterations(20000);
        ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");

        // Perform the 2D fit and store result
        hPhi1->Fit(fitFcn, "Q0ERN");
        TFitResultPtr r = hPhi1->Fit(fitFcn, "E0RSN");
        int fitStatus = static_cast<int>(r);
        int covQual = r.Get() ? r->CovMatrixStatus() : -1;

        fitFcnBkg->SetParameters(fitFcn->GetParameter(4), fitFcn->GetParameter(5), fitFcn->GetParameter(6));
        fitFcnSig->SetParameters(fitFcn->GetParameter(0), fitFcn->GetParameter(1), fitFcn->GetParameter(2), fitFcn->GetParameter(3));
        fitFcnBkg->SetLineColor(kBlue);
        fitFcnSig->SetLineColor(kGreen + 2);
        fitFcnBkg->SetLineStyle(2);
        fitFcnSig->SetLineStyle(2);
        // fitFcnSig->SetNpx(10000);
        fitFcnBkg->Draw("same");
        fitFcnSig->Draw("same");
        latex->DrawLatex(0.25, 0.92, Form("%.1f < #it{p}_{T} < %.1f GeV/c", pTBins[ibin], pTBins[ibin + 1]));
        // latex->DrawLatex(0.12, 0.82, Form("#Chi^{2}/NDF = %d", static_cast<int>(r->Chi2() / r->Ndf())));
        latex->DrawLatex(0.12, 0.82, Form("#Gamma = %.4f", fitFcnSig->GetParameter(2)));
        latex->DrawLatex(0.12, 0.74, Form("M_{#Phi} = %.4f", fitFcnSig->GetParameter(1)));
        latex->DrawLatex(0.12, 0.66, Form("Fit Status = %d", fitStatus));

        fitMass.push_back({fitFcnSig->GetParameter(1), fitFcnSig->GetParError(1)});
        fitResolution.push_back({fitFcnSig->GetParameter(2) * 1000, fitFcnSig->GetParError(2) * 1000});
        fitWidth.push_back({fitFcnSig->GetParameter(3), fitFcnSig->GetParError(3)});

        double intLow = fitFcnSig->GetParameter(1) - 1 * 0.005;
        double intHigh = fitFcnSig->GetParameter(1) + 1 * 0.005;

        double sigBkg = hPhi1->Integral(hPhi1->GetXaxis()->FindBin(intLow), hPhi1->GetXaxis()->FindBin(intHigh));
        double signal = fitFcnSig->Integral(intLow, intHigh) / hPhi1->GetBinWidth(1);
        double bkgOnly = fitFcnBkg->Integral(intLow, intHigh) / hPhi1->GetBinWidth(1);

        double signalOnly = sigBkg - bkgOnly;
        double purity = signalOnly / sigBkg;
        double purity2 = signal / (signal + bkgOnly);

        gPurity->SetPoint(ibin, (pTBins[ibin] + pTBins[ibin + 1]) / 2.0, purity);
        gPurity->SetPointError(ibin, (pTBins[ibin + 1] - pTBins[ibin]) / 2.0, 0);
    }
    cPhiVsPt->SaveAs(savepath + "/PhiInvMassVsPt.png");

    TCanvas *cPurityVsPt = new TCanvas("cPurityVsPt", "Purity vs Pt", 720, 720);
    SetCanvasStyle(cPurityVsPt, 0.15, 0.03, 0.05, 0.15);
    SetGraphErrorStyle(gPurity);
    gPurity->SetMarkerStyle(20);
    gPurity->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gPurity->GetYaxis()->SetTitle("Purity (S/(S+B))");
    gPurity->GetYaxis()->SetRangeUser(0.0, 1.0);
    gPurity->Draw("AP");
    latex->SetTextSize(0.035);
    latex->DrawLatex(0.25, 0.88, "Fit function: Voigtian + Pol2");
    latex->DrawLatex(0.25, 0.81, "Purity window: #it{M}_{#Phi} #pm 0.005 GeV/#it{c}^{2}");
    // cPurityVsPt->SaveAs(savepath + "/PhiPurityVsPt.png");

    TCanvas *cMassVsPt = new TCanvas("cMassVsPt", "Mass vs Pt", 720, 720);
    SetCanvasStyle(cMassVsPt, 0.20, 0.03, 0.05, 0.15);
    TGraphErrors *gMassVsPt = new TGraphErrors(fitMass.size());
    for (size_t i = 0; i < fitMass.size(); i++)
    {
        gMassVsPt->SetPoint(i, (pTBins[i] + pTBins[i + 1]) / 2.0, fitMass[i].first);
        gMassVsPt->SetPointError(i, (pTBins[i + 1] - pTBins[i]) / 2.0, fitMass[i].second);
    }
    gMassVsPt->SetMarkerStyle(20);
    gMassVsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gMassVsPt->GetYaxis()->SetTitle("M_{#Phi} (GeV/#it{c}^{2})");
    gMassVsPt->GetYaxis()->SetRangeUser(1.0157, 1.0223);
    SetGraphErrorStyle(gMassVsPt);
    gMassVsPt->GetYaxis()->SetTitleOffset(2.1);
    gMassVsPt->Draw("AP");
    TLine *linePDGMass = new TLine(pTBins[0], 1.019461, pTBins[13], 1.019461);
    linePDGMass->SetLineColor(kRed);
    linePDGMass->SetLineStyle(2);
    linePDGMass->SetLineWidth(2);
    linePDGMass->Draw("same");
    latex->SetTextSize(0.035);
    latex->DrawLatex(0.25, 0.88, "Fit function: Voigtian + Pol2");
    latex->DrawLatex(0.25, 0.81, "Purity window: #it{M}_{#Phi} #pm 0.005 GeV/#it{c}^{2}");
    // cMassVsPt->SaveAs(savepath + "/PhiMassVsPt.png");

    TCanvas *cResolutionVsPt = new TCanvas("cResolutionVsPt", "Width vs Pt", 720, 720);
    SetCanvasStyle(cResolutionVsPt, 0.15, 0.03, 0.05, 0.15);
    TGraphErrors *gResolutionvsPt = new TGraphErrors(fitResolution.size());
    for (size_t i = 0; i < fitResolution.size(); i++)
    {
        gResolutionvsPt->SetPoint(i, (pTBins[i] + pTBins[i + 1]) / 2.0, fitResolution[i].first);
        gResolutionvsPt->SetPointError(i, (pTBins[i + 1] - pTBins[i]) / 2.0, fitResolution[i].second);
    }
    gResolutionvsPt->SetMarkerStyle(20);
    gResolutionvsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gResolutionvsPt->GetYaxis()->SetTitle("Resolution (MeV/#it{c}^{2})");
    SetGraphErrorStyle(gResolutionvsPt);
    gResolutionvsPt->SetMinimum(0.6);
    gResolutionvsPt->SetMaximum(4.4);
    gResolutionvsPt->Draw("AP");
    latex->DrawLatex(0.25, 0.88, "Fit function: Voigtian + Pol2");
    latex->DrawLatex(0.25, 0.81, "Purity window: #it{M}_{#Phi} #pm 0.005 GeV/#it{c}^{2}");
    // cResolutionVsPt->SaveAs(savepath + "/PhiResolutionVsPt.png");

    // // Phi mass correlation plot
    // SetHistoQA(hPhiPhiMass);
    // TCanvas *cPhiPhiMassCorr = new TCanvas("cPhiPhiMassCorr", "Phi Mass vs Pt", 720, 720);
    // SetCanvasStyle(cPhiPhiMassCorr, 0.15, 0.15, 0.05, 0.15);
    // SetHistoQA(hPhiPhiMass);
    // hPhiPhiMass->GetXaxis()->SetTitle("M_{#phi1} (GeV/#it{c})");
    // hPhiPhiMass->GetYaxis()->SetTitle("M_{#phi2} (GeV/#it{c})");
    // hPhiPhiMass->GetYaxis()->SetTitleOffset(1.4);
    // hPhiPhiMass->GetYaxis()->SetNdivisions(505);
    // hPhiPhiMass->GetXaxis()->SetNdivisions(505);
    // hPhiPhiMass->GetZaxis()->SetMaxDigits(3);
    // hPhiPhiMass->Draw("colz");
    // TEllipse *circle = new TEllipse(1.0198, 1.0198, 0.005);
    // circle->SetFillStyle(0); // no fill
    // circle->SetLineColor(kRed);
    // circle->SetLineWidth(3);
    // circle->Draw("same");
    // latex->SetTextSize(0.035);
    // latex->SetTextColor(kWhite);
    // latex->DrawLatex(0.2, 0.91, "Red Circle: #DeltaM < 0.005 GeV/#it{c}^{2}");
    // latex->DrawLatex(0.2, 0.85, "#it{p}_{T}^{#phi#phi} > 9 GeV/#it{c}");
    // cPhiPhiMassCorr->SaveAs(savepath + "/PhiMassCorrelation.png");

    TFile *fPhiParams = new TFile(savepath + "/PhiParams26.root", "recreate");
    gMassVsPt->Write("gMassVsPt");
    gResolutionvsPt->Write("gResolutionVsPt");
    gPurity->Write("gPurity");
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
