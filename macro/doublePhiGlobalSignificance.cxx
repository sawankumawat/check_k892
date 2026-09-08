#include <iostream>
#include "TRandom3.h" // Added for toy pseudo-experiments
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

void doublePhiGlobalSignificance()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/TetraquarkFit";
    string suffix = "25_aiamShifted";
    int rebinFactor = 6;
    float fitRangeLow = 2.41;
    float fitRangeHigh = 2.96;

    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResults26_latest.root");

    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassPhiPhiRefitted");

    TH1D *hPt = hUnlike->Projection(1, "E");
    TH1D *hDeltaM = hUnlike->Projection(2, "E");

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
    SetHistoQA(hInvMass);
    hInvMass->Rebin(rebinFactor);
    hInvMass->GetXaxis()->SetRangeUser(fitRangeLow, fitRangeHigh);
    hInvMass->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hInvMass->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hInvMass->GetBinWidth(1) * 1000));

    int lowDeltaMBkg = hUnlike->GetAxis(2)->FindBin(0.005 + 0.00001);
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

    TH1D *hBkg = (TH1D *)hInvMass->Clone("hBkg");
    hBkg->Reset();

    for (int i = 1; i <= hInvMass->GetNbinsX(); i++)
    {
        double x = hInvMass->GetBinCenter(i);
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
    hInvMass->SetMinimum(hInvMass->GetMinimum() * 0.7);
    hInvMass->Draw("pe");

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(22);
    latex->SetTextSize(0.03);

    int lowNormBin = hInvMassBkg->GetXaxis()->FindBin(2.85 + 0.00001);
    int highNormBin = hInvMassBkg->GetXaxis()->FindBin(2.90 - 0.00001);
    double SigCounts = hInvMass->Integral(lowNormBin, highNormBin);
    double BkgCounts = hInvMassBkg->Integral(lowNormBin, highNormBin);
    double scaleFactor = SigCounts / BkgCounts;
    hInvMassBkg->Scale(scaleFactor);
    hInvMassBkg->SetLineColor(kRed);
    hInvMassBkg->SetMarkerColor(kRed);

    TF1 *fInitialCombinedFit = new TF1("fInitialCombinedFit", BWExpol, FIT_MIN, FIT_MAX, 7);
    fInitialCombinedFit->SetParNames("SignalYield", "Mass", "Width", "p0", "p1", "p2", "p3");
    fInitialCombinedFit->SetParameter(0, 100.0);
    fInitialCombinedFit->SetParameter(1, 2.70);
    fInitialCombinedFit->SetParameter(2, 0.03);

    fInitialCombinedFit->SetParameter(3, fitBkg->GetParameter(0));
    fInitialCombinedFit->SetParameter(4, fitBkg->GetParameter(1));
    fInitialCombinedFit->SetParameter(5, fitBkg->GetParameter(2));
    fInitialCombinedFit->SetParameter(6, fitBkg->GetParameter(3));

    fInitialCombinedFit->SetParLimits(0, 0.0, 1.0e6);
    fInitialCombinedFit->SetParLimits(1, 2.65, 2.75);
    fInitialCombinedFit->SetParLimits(2, 0.01, 0.05);
    fInitialCombinedFit->SetLineStyle(2);
    fInitialCombinedFit->SetLineColor(kMagenta);
    TFitResultPtr fitResultInitial = hInvMass->Fit(fInitialCombinedFit, "REBMS");

    // Print fit status and Covariance matrix status
    double fitStatus = fitResultInitial->Status();
    double covMatrixStatus = fitResultInitial->CovMatrixStatus();
    cout << "Initial Combined Fit Status: " << fitStatus << endl;
    cout << "Initial Combined Fit Covariance Matrix Status: " << covMatrixStatus << endl;

    //============================================================
    // Likelihood fit (S + B)
    //============================================================
    TF1 *fitFunc = new TF1("fitFunc", BWExpol, FIT_MIN, FIT_MAX, 7);
    fitFunc->SetParNames("SignalYield", "Mass", "Width", "p0", "p1", "p2", "p3");

    fitFunc->SetParameter(0, fInitialCombinedFit->GetParameter(0));
    fitFunc->SetParameter(1, fInitialCombinedFit->GetParameter(1));
    fitFunc->SetParameter(2, fInitialCombinedFit->GetParameter(2));

    fitFunc->SetParameter(3, fInitialCombinedFit->GetParameter(3));
    fitFunc->SetParameter(4, fInitialCombinedFit->GetParameter(4));
    fitFunc->SetParameter(5, fInitialCombinedFit->GetParameter(5));
    fitFunc->SetParameter(6, fInitialCombinedFit->GetParameter(6));

    fitFunc->SetParLimits(0, 0.0, 1.0e6);
    fitFunc->SetParLimits(1, 2.65, 2.75);
    fitFunc->SetParLimits(2, 0.01, 0.15);

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

    double logL_SB = fitResultSB->MinFcnValue();

    //============================================================
    // Background-only likelihood fit
    //============================================================
    TCanvas *cBkgOnlyFitForToys = new TCanvas("cBkgOnlyFitForToys", "Background-only Fit for Toys", 720, 720);
    SetCanvasStyle(cBkgOnlyFitForToys, 0.15, 0.03, 0.05, 0.15);
    hInvMass->Draw("pe");

    TF1 *fitBOnly = new TF1("fitBOnly", expPol3, FIT_MIN, FIT_MAX, 4);
    fitBOnly->SetParNames("p0", "p1", "p2", "p3");
    fitBOnly->SetParameters(fInitialCombinedFit->GetParameter(3), fInitialCombinedFit->GetParameter(4), fInitialCombinedFit->GetParameter(5), fInitialCombinedFit->GetParameter(6));
    TFitResultPtr fitResultB = hInvMass->Fit(fitBOnly, "RLS");
    fitBOnly->SetLineColor(kBlue + 1);
    fitBOnly->SetLineWidth(2);
    fitBOnly->Draw("same");

    if (fitResultB.Get() == nullptr)
    {
        cout << "ERROR: Background-only likelihood fit failed!" << endl;
        return;
    }
    double logL_B = fitResultB->MinFcnValue();

    // Local Test Statistic
    double q0_data = 2.0 * (logL_B - logL_SB);
    if (q0_data < 0)
        q0_data = 0.0;
    double localSignificance = TMath::Sqrt(q0_data);
    double localPValue = 0.5 * TMath::Erfc(localSignificance / TMath::Sqrt(2.0));

    //============================================================
    // GLOBAL SIGNIFICANCE (Background-only Toys)
    //============================================================
    const int nToys = 1e7; // Increase to 10000+ for higher precision
    int countExceedingData = 0;

    TRandom3 randGen(42);
    TH1D *hToy = (TH1D *)hInvMass->Clone("hToy");

    TF1 *fitToySB = new TF1("fitToySB", BWExpol, FIT_MIN, FIT_MAX, 7);
    fitToySB->SetParNames("SignalYield", "Mass", "Width", "p0", "p1", "p2", "p3");

    TF1 *fitToyBOnly = new TF1("fitToyBOnly", expPol3, FIT_MIN, FIT_MAX, 4);
    fitToyBOnly->SetParNames("p0", "p1", "p2", "p3");

    cout << "\nRunning " << nToys << " background-only pseudo-experiments for global significance..." << endl;

    for (int itoy = 0; itoy < nToys; ++itoy)
    {
        // 1. Generate toy dataset from B-only PDF using Poisson fluctuation[cite: 1]
        hToy->Reset();
        for (int ibin = 1; ibin <= hToy->GetNbinsX(); ++ibin)
        {
            double binCenter = hToy->GetBinCenter(ibin);
            double binWidth = hToy->GetBinWidth(ibin);
            double expectedBkg = fitBOnly->Eval(binCenter) * binWidth;
            double toyCounts = randGen.Poisson(expectedBkg);

            hToy->SetBinContent(ibin, toyCounts);
            hToy->SetBinError(ibin, TMath::Sqrt(toyCounts));
        }

        // 2. Fit toy with B-only model
        fitToyBOnly->SetParameters(fitBOnly->GetParameters());
        TFitResultPtr resToyB = hToy->Fit(fitToyBOnly, "RLS0Q");

        // 3. Fit toy with S+B model (floating Mass & Width) (for LEE, need to float Mass and Width in full range)
        fitToySB->SetParameters(100.0, 2.70, 0.03, fitBOnly->GetParameter(0), fitBOnly->GetParameter(1), fitBOnly->GetParameter(2), fitBOnly->GetParameter(3));

        fitToySB->SetParLimits(0, 0.0, 1.0e6);
        fitToySB->SetParLimits(1, 2.65, 2.75); // Scan the mass in the +- 5sigma range from theory
        fitToySB->SetParLimits(2, 0.01, 0.15); // Allow width to vary in a reasonable range

        TFitResultPtr resToySB = hToy->Fit(fitToySB, "RLS0Q");

        if (resToyB.Get() != nullptr && resToySB.Get() != nullptr && resToyB->IsValid() && resToySB->IsValid())
        {
            double q0_toy = 2.0 * (resToyB->MinFcnValue() - resToySB->MinFcnValue());
            if (q0_toy < 0)
                q0_toy = 0.0;

            if (q0_toy >= q0_data)
                countExceedingData++;
        }
    }

    cout << "Number of toys exceeding data q0: " << countExceedingData << " out of " << nToys << endl;

    double globalPValue = static_cast<double>(countExceedingData) / nToys;
    double globalSignificance = 0.0;
    if (globalPValue > 0)
        globalSignificance = TMath::ErfInverse(1.0 - 2.0 * globalPValue) * TMath::Sqrt(2.0);

    cout << endl;
    cout << "==========================================" << endl;
    cout << "       LOCAL vs GLOBAL RESULTS" << endl;
    cout << "==========================================" << endl;
    cout << "Data q0 = " << q0_data << endl;
    cout << "Local p-value = " << localPValue << " | Local Z = " << localSignificance << " sigma" << endl;
    cout << "Global p-value = " << globalPValue << " | Global Z = " << globalSignificance << " sigma" << endl;
    cout << "==========================================" << endl;

    //============================================================
    // Draw S+B fit
    //============================================================
    fitFunc->SetLineColor(kRed + 1);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("same");

    TF1 *fitBkgFinal = new TF1("fitBkgFinal", expPol3, FIT_MIN, FIT_MAX, 4);
    fitBkgFinal->SetParameters(fitFunc->GetParameter(3), fitFunc->GetParameter(4), fitFunc->GetParameter(5), fitFunc->GetParameter(6));
    fitBkgFinal->SetLineColor(kBlue);
    fitBkgFinal->SetLineStyle(2);
    fitBkgFinal->Draw("same");

    TF1 *fitSignal = new TF1("fitSignal", breitWigner, FIT_MIN, FIT_MAX, 3);
    fitSignal->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2));
    fitSignal->SetLineColor(kRed);
    fitSignal->SetLineStyle(2);
    fitSignal->Draw("same");

    TLegend *legend = new TLegend(0.53, 0.72, 0.9, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(hInvMass, "#Delta M < 0.005, #it{p}_{T}^{#phi#phi} > 9 GeV/#it{c}", "pe");
    legend->AddEntry(fitFunc, "BW + expol3", "l");
    legend->AddEntry(fitBkgFinal, "expol3 (bkg)", "l");
    legend->AddEntry(fitSignal, "Breit-Wigner", "l");
    legend->Draw();

    latex->DrawLatex(0.25, 0.85, "LHC25(am+ai)_skimmed");
    latex->DrawLatex(0.2, 0.42, Form("M = %.4f #pm %.4f GeV/c^{2}", massFit, massErr));
    latex->DrawLatex(0.2, 0.37, Form("#Gamma = %.4f #pm %.4f GeV/c^{2}", widthFit, widthErr));
    latex->DrawLatex(0.2, 0.32, Form("#Chi^{2}/NDF = %.2f / %d", fitFunc->GetChisquare(), fitFunc->GetNDF()));
    latex->DrawLatex(0.2, 0.27, Form("Local p = %.3e (%.2f #sigma)", localPValue, localSignificance));
    latex->DrawLatex(0.2, 0.22, Form("Global p = %.3e (%.2f #sigma)", globalPValue, globalSignificance));

    delete hToy;
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
        std::cout << "Error: histo " << name
                  << " not found in file " << f->GetName() << std::endl;
        return nullptr;
    }

    return histo;
}