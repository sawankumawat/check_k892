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

Double_t BWExpol(Double_t *x, Double_t *par)
{
    return breitWigner(x, par) + expPol3(x, &par[3]);
}

void doublePhiSimpleOpti7()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests";
    ////======Pair=========
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/Pair/processopti5/AnalysisResults.root");

    ////=====New===========
    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti7/AnalysisResults.root");

    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassPhiPhi");

    TH1D *hPt = hUnlike->Projection(1, "E");
    TH1D *hDeltaM = hUnlike->Projection(2, "E");

    int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowDeltaM = hUnlike->GetAxis(2)->FindBin(0.0 + 0.00001);
    int highDeltaM = hUnlike->GetAxis(2)->FindBin(0.005 - 0.00001);

    hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    hUnlike->GetAxis(2)->SetRange(lowDeltaM, highDeltaM);

    TH1D *hInvMass = hUnlike->Projection(0, "E");
    SetHistoQA(hInvMass);
    hInvMass->Rebin(8);
    hInvMass->GetXaxis()->SetRangeUser(2.5, 2.9);
    hInvMass->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hInvMass->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hInvMass->GetBinWidth(1) * 1000));

    TH1D *hBkg = (TH1D *)hInvMass->Clone("hBkg");
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

    TCanvas *cBkg = new TCanvas("cBkg", "Background", 720, 720);
    SetCanvasStyle(cBkg, 0.15, 0.03, 0.05, 0.15);
    hBkg->Draw("pe");

    const double FIT_MIN = 2.5;
    const double FIT_MAX = 2.9;

    TF1 *fitBkg = new TF1("fitBkg", expPol3, FIT_MIN, FIT_MAX, 4);
    fitBkg->SetParNames("p0", "p1", "p2", "p3");
    fitBkg->SetParameters(1.0, 1.0, 1.0, 1.0);
    hBkg->Fit(fitBkg, "REBM");

    TCanvas *cInvMass = new TCanvas("cInvMass", "Invariant Mass", 720, 720);
    SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
    hInvMass->Draw("pe");

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

    fitFunc->SetParameter(3, fInitialCombinedFit->GetParameter(3));
    fitFunc->SetParameter(4, fInitialCombinedFit->GetParameter(4));
    fitFunc->SetParameter(5, fInitialCombinedFit->GetParameter(5));
    fitFunc->SetParameter(6, fInitialCombinedFit->GetParameter(6));

    fitFunc->SetParLimits(0, 0.0, 1.0e6);
    fitFunc->SetParLimits(1, 2.65, 2.75);
    fitFunc->SetParLimits(2, 0.01, 0.05);

    TFitResultPtr fitResultSB = hInvMass->Fit(fitFunc, "RLS");

    if (fitResultSB.Get() == nullptr)
    {
        cout << "ERROR: S+B likelihood fit failed!" << endl;
        return;
    }

    double signalYield = fitFunc->GetParameter(0);
    double signalYieldErr = fitFunc->GetParError(0);
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
    fitFunc->SetLineColor(kBlack);
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

    TLegend *legend = new TLegend(0.65, 0.75, 0.9, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(hInvMass, "Data", "pe");
    legend->AddEntry(fitFunc, "BW + expol3", "l");
    legend->AddEntry(fitBkgFinal, "expol3 (bkg)", "l");
    legend->AddEntry(fitSignal, "Breit-Wigner", "l");
    legend->Draw();

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(22);
    latex->SetTextSize(0.03);
    latex->DrawLatex(0.2, 0.45, Form("M = %.4f #pm %.4f GeV/c^{2}", massFit, massErr));
    latex->DrawLatex(0.2, 0.40, Form("#Gamma = %.4f #pm %.4f GeV/c^{2}", widthFit, widthErr));
    // latex->DrawLatex(0.2, 0.35, Form("N_{sig} = %.1f #pm %.1f", signalYield, signalYieldErr));
    latex->DrawLatex(0.2, 0.35, Form("#Chi^{2}/NDF = %.2f / %d", fitFunc->GetChisquare(), fitFunc->GetNDF()));
    latex->DrawLatex(0.2, 0.30, Form("p-value = %.3e", pValue));
    latex->DrawLatex(0.2, 0.25, Form("Significance (Z) = %.2f #sigma", significance));
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
