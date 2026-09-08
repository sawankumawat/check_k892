#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"

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

Double_t voigtFunction(Double_t *x, Double_t *par)
{
    return (par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]));
}

Double_t BWExpol(Double_t *x, Double_t *par)
{
    return breitWigner(x, par) + expPol3(x, &par[3]);
}

Double_t VoigtExpol(Double_t *x, Double_t *par)
{
    double vgt = voigtFunction(x, par);
    double poly3 = expPol3(x, &par[4]);
    return (vgt + poly3);
}

void doublePhiOpti8()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    // TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/TetraquarkFit/pTCutVariation26";
    // TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/TetraquarkFit/pTCutVariation25";
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/TetraquarkFit/PeriodWise25Data";
    // string suffix = "25_aiam";
    string suffix = "25";
    // string suffix = "26";
    // string suffix = "26afai";
    // string suffix = "26ac";
    // string suffix = "26adaeag";
    int rebinFactor = 8;
    float fitRangeLow = 2.41;
    float fitRangeHigh = 2.96;

    ////=====New===========
    ////====2026 data========
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResultsRefit2.root");
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResults_LHC26_PID2003.root");
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResults26_latest.root");
    // TFile *fInput = OpenFile(Form("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/All/AnalysisResults_LHC%s.root", suffix.c_str()));

    ////======2025 data==========
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/LHC25/AnalysisResultsLHC25.root");
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/LHC25/AnalysisResultsKaShift2.root");
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/LHC25/AnalysisResults_LHC25_PID2003.root");
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti9/LHC25/AnalysisResults_LHC25_KaShifted_BhaiyaCode2.root"); // opti9, shifted
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/LHC25/AnalysisResults25_aiam.root"); //opti8, 25 ai+am
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti9/LHC25/AnalysisResults25_aiamShifted.root"); // opti9, shifted
    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti10/LHC25/AnalysisResults25aiam_morepTbins.root"); // opti10, shifted
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti9/LHC25/AnalysisResults25ai.root"); // opti9, shifted

    // THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassPhiPhiRefitted");
    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassPhiPhiShifted");
    // THnSparseF *hRot = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassPhiPhiRotational");

    TH1D *hPt = hUnlike->Projection(1, "E");
    TH1D *hDeltaM = hUnlike->Projection(2, "E");

    double ptCut = 9.0;

    int lowpT = hUnlike->GetAxis(1)->FindBin(ptCut + 0.001);
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

    // hRot->GetAxis(1)->SetRange(lowpT, highpT);

    TH1D *hInvMass = hUnlike->Projection(0, "E");
    SetHistoQA(hInvMass);
    hInvMass->Rebin(rebinFactor);
    hInvMass->GetXaxis()->SetRangeUser(fitRangeLow, fitRangeHigh);
    hInvMass->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hInvMass->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hInvMass->GetBinWidth(1) * 1000));

    int lowDeltaMBkg = hUnlike->GetAxis(2)->FindBin(0.01 + 0.00001);
    int highDeltaMBkg = hUnlike->GetAxis(2)->FindBin(100 - 0.00001);

    int lowChi2Bkg = hUnlike->GetAxis(3)->FindBin(0.0 + 0.00001);
    int highChi2Bkg = hUnlike->GetAxis(3)->FindBin(30.0 - 0.00001);

    int lowFitProbBkg = hUnlike->GetAxis(4)->FindBin(0.0 + 0.00001);
    int highFitProbBkg = hUnlike->GetAxis(4)->FindBin(2.0 - 0.00001);

    hUnlike->GetAxis(2)->SetRange(lowDeltaMBkg, highDeltaMBkg);
    hUnlike->GetAxis(3)->SetRange(lowChi2Bkg, highChi2Bkg);
    hUnlike->GetAxis(4)->SetRange(lowFitProbBkg, highFitProbBkg);

    TH1D *hInvMassBkg = hUnlike->Projection(0, "E");
    SetHistoQA(hInvMassBkg);
    hInvMassBkg->Rebin(rebinFactor);
    hInvMassBkg->GetXaxis()->SetRangeUser(fitRangeLow, fitRangeHigh);

    // TH1D *h1DRot = hRot->Projection(0, "E");
    // h1DRot->Rebin(rebinFactor);
    // SetHistoQA(h1DRot);

    TH1D *hBkg = (TH1D *)hInvMass->Clone("hBkg");
    hBkg->Reset();

    for (int i = 1; i <= hInvMass->GetNbinsX(); i++)
    {
        double x = hInvMass->GetBinCenter(i);

        // Exclude signal region
        if (x >= 2.62 && x <= 2.73)
            continue;

        hBkg->SetBinContent(i, hInvMass->GetBinContent(i));
        hBkg->SetBinError(i, hInvMass->GetBinError(i));
    }

    TCanvas *cBkg = new TCanvas("cBkg", "Background", 720, 720);
    SetCanvasStyle(cBkg, 0.15, 0.03, 0.05, 0.15);
    hBkg->Draw("pe");

    const double FIT_MIN = fitRangeLow;
    const double FIT_MAX = fitRangeHigh;

    TF1 *fitBkg = new TF1("fitBkg", expPol3, FIT_MIN, FIT_MAX, 4);
    fitBkg->SetParNames("p0", "p1", "p2", "p3");
    fitBkg->SetParameters(-2.7e2, 64.0, 40.0, -6.5);
    hBkg->Fit(fitBkg, "REBM");

    TCanvas *cInvMass = new TCanvas("cInvMass", "Invariant Mass", 720, 720);
    SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
    // hInvMass->GetYaxis()->SetRangeUser(1.4e3, 3.3e3);
    hInvMass->SetMinimum(hInvMass->GetMinimum() * 0.7);
    hInvMass->Draw("pe");

    // TLine *lineat2p65 = new TLine(2.65, 971, 2.65, 2058);
    // lineat2p65->SetLineColor(kRed);
    // lineat2p65->SetLineStyle(2);
    // lineat2p65->SetLineWidth(2);
    // lineat2p65->Draw();

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(22);
    latex->SetTextSize(0.03);

    // TLatex *latex2 = new TLatex();
    // latex2->SetNDC();
    // latex2->SetTextFont(42);
    // latex2->SetTextSize(0.03);
    // latex2->DrawLatex(0.55, 0.9, "LHC25_skimmed");
    // latex2->DrawLatex(0.55, 0.85, "Before Kaon Momentum shift");

    // // Normalize the rotational background
    // int lowNormBin = h1DRot->GetXaxis()->FindBin(2.85 + 0.00001);
    // int highNormBin = h1DRot->GetXaxis()->FindBin(2.9 - 0.00001);
    // double SigCounts = hInvMass->Integral(lowNormBin, highNormBin);
    // double RotCounts = h1DRot->Integral(lowNormBin, highNormBin);
    // double scaleFactor = SigCounts / RotCounts;
    // h1DRot->Scale(scaleFactor);
    // h1DRot->SetLineColor(kRed);
    // h1DRot->SetMarkerColor(kRed);
    // h1DRot->Draw("pe same");

    // // Normalize the background from deltaM>0.005
    int lowNormBin = hInvMassBkg->GetXaxis()->FindBin(2.85 + 0.00001);
    int highNormBin = hInvMassBkg->GetXaxis()->FindBin(2.90 - 0.00001);
    double SigCounts = hInvMass->Integral(lowNormBin, highNormBin);
    double BkgCounts = hInvMassBkg->Integral(lowNormBin, highNormBin);
    double scaleFactor = SigCounts / BkgCounts;
    hInvMassBkg->Scale(scaleFactor);
    hInvMassBkg->SetLineColor(kRed);
    hInvMassBkg->SetMarkerColor(kRed);
    // hInvMassBkg->Draw("HISTe same");

    // TLegend *legTemp = new TLegend(0.35, 0.75, 0.9, 0.92);
    // legTemp->SetBorderSize(0);
    // legTemp->SetFillStyle(0);
    // legTemp->SetTextFont(42);
    // legTemp->SetTextSize(0.03);
    // legTemp->AddEntry((TObject *)0, "#Delta#it{M}_{#phi} < 0.005", "");
    // legTemp->AddEntry((TObject *)0, "#it{p}_{T}^{#phi#phi} > 9 GeV/#it{c}", "");
    // legTemp->AddEntry(hInvMassBkg, "#Delta#it{M}_{#phi} #geq 0.005", "l");
    // legTemp->AddEntry((TObject *)0, "Normalized in 2.85 < #it{M}_{#phi#phi} < 2.90 GeV/#it{c}^{2}", "");
    // legTemp->Draw();
    // cInvMass->SaveAs(savepath + "/InvariantMassWithoutFit.png");

    TF1 *fInitialCombinedFit = new TF1("fInitialCombinedFit", BWExpol, FIT_MIN, FIT_MAX, 7);
    fInitialCombinedFit->SetParNames("SignalYield", "Mass", "Width", "p0", "p1", "p2", "p3");
    fInitialCombinedFit->SetParameter(0, 100.0); // Signal yield
    fInitialCombinedFit->SetParameter(1, 2.70);  // Mass
    fInitialCombinedFit->SetParameter(2, 0.03);  // Width

    fInitialCombinedFit->SetParameter(3, fitBkg->GetParameter(0));
    fInitialCombinedFit->SetParameter(4, fitBkg->GetParameter(1));
    fInitialCombinedFit->SetParameter(5, fitBkg->GetParameter(2));
    fInitialCombinedFit->SetParameter(6, fitBkg->GetParameter(3));

    fInitialCombinedFit->SetParLimits(0, 0.0, 1.0e6);
    fInitialCombinedFit->SetParLimits(1, 2.65, 2.75);
    fInitialCombinedFit->SetParLimits(2, 0.01, 0.05);
    fInitialCombinedFit->SetLineStyle(2);
    fInitialCombinedFit->SetLineColor(kMagenta);
    hInvMass->Fit(fInitialCombinedFit, "REBMS");

    //============================================================
    // Likelihood fit (S + B)
    //============================================================

    TF1 *fitFunc = new TF1("fitFunc", BWExpol, FIT_MIN, FIT_MAX, 7);
    fitFunc->SetParNames("SignalYield", "Mass", "Width", "p0", "p1", "p2", "p3");

    fitFunc->SetParameter(0, fInitialCombinedFit->GetParameter(0)); // Signal yield
    fitFunc->SetParameter(1, fInitialCombinedFit->GetParameter(1)); // Mass
    fitFunc->SetParameter(2, fInitialCombinedFit->GetParameter(2)); // Width

    fitFunc->FixParameter(3, fInitialCombinedFit->GetParameter(3));
    fitFunc->FixParameter(4, fInitialCombinedFit->GetParameter(4));
    fitFunc->FixParameter(5, fInitialCombinedFit->GetParameter(5));
    fitFunc->FixParameter(6, fInitialCombinedFit->GetParameter(6));

    // fitFunc->FixParameter(3, -299.779);
    // fitFunc->FixParameter(4, 63.05);
    // fitFunc->FixParameter(5, 37.87);
    // fitFunc->FixParameter(6, -7.178);

    fitFunc->SetParLimits(0, 0.0, 1.0e6);
    fitFunc->SetParLimits(1, 2.65, 2.75);
    fitFunc->SetParLimits(2, 0.01, 0.15);

    TFitResultPtr fitResultSB = hInvMass->Fit(fitFunc, "RLS");

    if (fitResultSB.Get() == nullptr)
    {
        cout << "ERROR: S+B likelihood fit failed!" << endl;
        return;
    }

    double signalYield = fitFunc->GetParameter(0) / hInvMass->GetBinWidth(1);
    double signalYieldErr = fitFunc->GetParError(0) / hInvMass->GetBinWidth(1);
    double massFit = fitFunc->GetParameter(1);
    double massErr = fitFunc->GetParError(1);
    double widthFit = fitFunc->GetParameter(2);
    double widthErr = fitFunc->GetParError(2);

    cout << endl;
    cout << "==========================================" << endl;
    cout << "       S + B LIKELIHOOD FIT" << endl;
    cout << "==========================================" << endl;

    cout << "Mass  = " << massFit << " +/- " << massErr << " GeV/c2" << endl;
    cout << "Width = " << widthFit << " +/- " << widthErr << " GeV/c2" << endl;
    cout << "Signal yield S = " << signalYield << " +/- " << signalYieldErr << endl;

    double logL_SB = fitResultSB->MinFcnValue();

    //============================================================
    //      Background-only likelihood fit
    //============================================================
    TF1 *fitBOnly = new TF1("fitBOnly", expPol3, FIT_MIN, FIT_MAX, 4);
    fitBOnly->SetParNames("p0", "p1", "p2", "p3");
    fitBOnly->SetParameters(fInitialCombinedFit->GetParameter(3), fInitialCombinedFit->GetParameter(4), fInitialCombinedFit->GetParameter(5), fInitialCombinedFit->GetParameter(6));
    TFitResultPtr fitResultB = hInvMass->Fit(fitBOnly, "RLS0");

    if (fitResultB.Get() == nullptr)
    {
        cout << "ERROR: Background-only likelihood fit failed!" << endl;
        return;
    }
    double logL_B = fitResultB->MinFcnValue();

    //============================================================
    // Likelihood-ratio test statistic and significance
    //============================================================
    double q0 = 2.0 * (logL_B - logL_SB);

    if (q0 < 0)
        q0 = 0.0;

    double significance = TMath::Sqrt(q0);

    //============================================================
    // p-value
    //
    // For a one-sided discovery test:
    // p0 = 0.5 * erfc(Z/sqrt(2))
    //
    // The factor 0.5 accounts for the physical boundary
    // S >= 0 in the asymptotic Cowan et al. treatment.
    //============================================================

    double pValue = 0.5 * TMath::Erfc(significance / TMath::Sqrt(2.0));

    cout << endl;
    cout << "==========================================" << endl;
    cout << "       LIKELIHOOD-RATIO RESULTS" << endl;
    cout << "==========================================" << endl;

    cout << "S+B MinFcnValue = " << logL_SB << endl;
    cout << "B-only MinFcnValue = " << logL_B << endl;
    cout << "q0 = " << q0 << endl;
    cout << "Signal yield S = " << signalYield << " +/- " << signalYieldErr << endl;
    cout << "Significance Z = " << significance << " sigma" << endl;
    cout << "p0 = " << pValue << endl;
    cout << "==========================================" << endl;

    //============================================================
    // Draw S+B fit
    //============================================================

    // S+B total fit
    fitFunc->SetLineColor(kRed + 1);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("same");

    // Background from S+B fit
    TF1 *fitBkgFinal = new TF1("fitBkgFinal", expPol3, FIT_MIN, FIT_MAX, 4);

    fitBkgFinal->SetParameters(fitFunc->GetParameter(3), fitFunc->GetParameter(4), fitFunc->GetParameter(5), fitFunc->GetParameter(6));

    fitBkgFinal->SetLineColor(kBlue);
    fitBkgFinal->SetLineStyle(2);
    fitBkgFinal->Draw("same");

    // Signal component
    TF1 *fitSignal = new TF1("fitSignal", breitWigner, FIT_MIN, FIT_MAX, 3);
    fitSignal->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2));

    fitSignal->SetLineColor(kRed);
    fitSignal->SetLineStyle(2);
    fitSignal->Draw("same");
    cout << "Signal in +-3sigma is " << fitSignal->Integral(massFit - 10 * widthFit, massFit + 10 * widthFit) << endl;
    cout << "bkg p0 " << fitBkgFinal->GetParameter(0) << " p1 " << fitBkgFinal->GetParameter(1) << " p2 " << fitBkgFinal->GetParameter(2) << " p3 " << fitBkgFinal->GetParameter(3) << endl;

    TLegend *legend = new TLegend(0.51, 0.72, 0.9, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(hInvMass, Form("#Delta M < 0.005, #it{p}_{T}^{#phi#phi} > %.1f GeV/#it{c}", ptCut), "pe");
    legend->AddEntry(fitFunc, "BW + expol3", "l");
    legend->AddEntry(fitBkgFinal, "expol3 (bkg)", "l");
    legend->AddEntry(fitSignal, "Breit-Wigner", "l");
    legend->Draw();

    // latex->DrawLatex(0.25, 0.85, "LHC25(am+ai)_skimmed");
    latex->DrawLatex(0.2, 0.42, Form("M = %.4f #pm %.4f GeV/c^{2}", massFit, massErr));
    latex->DrawLatex(0.2, 0.37, Form("#Gamma = %.4f #pm %.4f GeV/c^{2}", widthFit, widthErr));
    // latex->DrawLatex(0.2, 0.35, Form("N_{sig} = %.1f #pm %.1f", signalYield, signalYieldErr));
    latex->DrawLatex(0.2, 0.32, Form("#Chi^{2}/NDF = %.2f / %d", fitFunc->GetChisquare(), fitFunc->GetNDF()));
    latex->DrawLatex(0.2, 0.27, Form("p-value = %.3e", pValue));
    latex->DrawLatex(0.2, 0.22, Form("Significance (Z) = %.2f #sigma", significance));
    // cInvMass->SaveAs(savepath + "/InvariantMassFit" + suffix + "_" + Form("%.1f", ptCut) + ".png");
}

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
