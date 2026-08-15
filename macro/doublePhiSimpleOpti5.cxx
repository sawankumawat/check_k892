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

void doublePhiSimpleOpti5()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests";
    ////======Pair=========
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/Pair/processopti5/AnalysisResults.root");

    ////=====New===========
    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults.root");

    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassDoublePhi");

    TH1D *hPt = hUnlike->Projection(1, "E");
    TH1D *hPhiPtAsymmetry = hUnlike->Projection(2, "E");
    TH1D *hRapidity = hUnlike->Projection(3, "E");
    TH1D *hPhi1Mass = hUnlike->Projection(4, "E");
    TH1D *hPhi2Mass = hUnlike->Projection(5, "E");
    TH1D *hDeltaM = hUnlike->Projection(6, "E");
    TH1D *hnKaonTOF = hUnlike->Projection(7, "E");
    TH1D *hCombinedPID4Kaon = hUnlike->Projection(8, "E");

    // vector<TH1D *> histos = {hPt, hPhiPtAsymmetry, hRapidity, hPhi1Mass, hPhi2Mass, hDeltaM, hnKaonTOF, hCombinedPID4Kaon};
    // //Draw all histograms using loop
    // for (auto histo : histos)
    // {
    //     SetHistoQA(histo);
    //     TCanvas *c = new TCanvas(Form("c%s", histo->GetName()), histo->GetTitle(), 720, 720);
    //     SetCanvasStyle(c, 0.15, 0.03, 0.05, 0.15);
    //     histo->Draw("pe");
    // }

    TCanvas *cPt = new TCanvas("cPt", "Pt", 720, 720);
    SetCanvasStyle(cPt, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hPt);
    hPt->Draw("pe");

    int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowPhiPtAsymmetry = hUnlike->GetAxis(2)->FindBin(0.0 + 0.0001);
    int highPhiPtAsymmetry = hUnlike->GetAxis(2)->FindBin(1.0 - 0.0001);

    int lowRapidity = hUnlike->GetAxis(3)->FindBin(-0.8 + 0.0001);
    int highRapidity = hUnlike->GetAxis(3)->FindBin(0.8 - 0.0001);

    int lowPhi1Mass = hUnlike->GetAxis(4)->FindBin(1.005 + 0.0001);
    int highPhi1Mass = hUnlike->GetAxis(4)->FindBin(1.035 - 0.0001);

    int lowPhi2Mass = hUnlike->GetAxis(5)->FindBin(1.005 + 0.0001);
    int highPhi2Mass = hUnlike->GetAxis(5)->FindBin(1.035 - 0.0001);

    int lowDeltaM = hUnlike->GetAxis(6)->FindBin(0.0 + 0.00001);
    int highDeltaM = hUnlike->GetAxis(6)->FindBin(0.005 - 0.00001);

    int lownKaonTOF = hUnlike->GetAxis(7)->FindBin(0.0 + 0.0001);
    int highnKaonTOF = hUnlike->GetAxis(7)->FindBin(4.0 - 0.0001);

    int lowCombinedPID4Kaon = hUnlike->GetAxis(8)->FindBin(0.0 + 0.0001);
    int highCombinedPID4Kaon = hUnlike->GetAxis(8)->FindBin(4.0 - 0.0001);

    hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    // hUnlike->GetAxis(2)->SetRange(lowPhiPtAsymmetry, highPhiPtAsymmetry);
    // hUnlike->GetAxis(3)->SetRange(lowRapidity, highRapidity);
    // hUnlike->GetAxis(4)->SetRange(lowPhi1Mass, highPhi1Mass);
    // hUnlike->GetAxis(5)->SetRange(lowPhi2Mass, highPhi2Mass);
    hUnlike->GetAxis(6)->SetRange(lowDeltaM, highDeltaM);
    // hUnlike->GetAxis(7)->SetRange(lownKaonTOF, highnKaonTOF);
    // hUnlike->GetAxis(8)->SetRange(lowCombinedPID4Kaon, highCombinedPID4Kaon);

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

    TCanvas *cOnlyBkg = new TCanvas("cOnlyBkg", "Only Background", 720, 720);
    SetCanvasStyle(cOnlyBkg, 0.15, 0.03, 0.05, 0.15);
    hBkg->Draw("pe");

    TF1 *fitBkg = new TF1("fitBkg", expPol3, 2.5, 2.9, 4);
    fitBkg->SetParNames("p0", "p1", "p2", "p3");
    fitBkg->SetParameters(1, 1, 1, 1);
    hBkg->Fit(fitBkg, "R");

    TCanvas *cInvMass = new TCanvas("cInvMass", "Invariant Mass", 720, 720);
    SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
    hInvMass->SetMaximum(4100);
    hInvMass->Draw("pe");

    // //======================Fit function (BW)=========================
    TF1 *fitFunc = new TF1("fitFunc", BWExpol, 2.5, 2.9, 7);
    fitFunc->SetParNames("Mass", "Width", "Yield", "p0", "p1", "p2", "p3");

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

    hInvMass->Fit(fitFunc, "REBMS");

    TF1 *fitBkgFinal = new TF1("fitBkgFinal", expPol3, 2.5, 2.9, 4);
    fitBkgFinal->SetParameters(fitFunc->GetParameter(3), fitFunc->GetParameter(4), fitFunc->GetParameter(5), fitFunc->GetParameter(6));
    fitBkgFinal->SetLineColor(kBlue);
    fitBkgFinal->SetLineStyle(2);
    fitBkgFinal->Draw("same");

    TF1 *fitBW = new TF1("fitBW", breitWigner, 2.5, 2.9, 3);
    fitBW->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2));
    fitBW->SetLineColor(kRed);
    fitBW->SetLineStyle(2);
    fitBW->Draw("same");
    cInvMass->SaveAs(savepath + "/DoublePhiFit.png");

    // Phi mass correlation plot
    TH3F *hPhiMassVsPt = GetHisto<TH3F>(fInput, "doublephimeson/hPhiMass");
    TH2F *hPhiPhiMass = (TH2F *)hPhiMassVsPt->Project3D("yx");
    SetHistoQA(hPhiPhiMass);
    TCanvas *cPhiPhiMassCorr = new TCanvas("cPhiPhiMassCorr", "Phi Mass vs Pt", 720, 720);
    SetCanvasStyle(cPhiPhiMassCorr, 0.18, 0.03, 0.05, 0.15);
    hPhiPhiMass->GetXaxis()->SetTitle("M_{#phi1} (GeV/#it{c})");
    hPhiPhiMass->GetYaxis()->SetTitle("M_{#phi2} (GeV/#it{c})");
    hPhiPhiMass->Draw("colz");

    TCanvas *cPhiVsPt = new TCanvas("cPhiVsPt", "Phi Mass vs Pt", 1080, 720);
    cPhiVsPt->Divide(4, 3);
    double pTBins[13] = {0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 20.0};
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(22);
    latex->SetTextSize(0.06);
    vector<pair<double, double>> fitMass, fitResolution, fitWidth;
    for (int ibin = 0; ibin < 12; ibin++)
    {
        int lowpTBin = hPhiMassVsPt->GetZaxis()->FindBin(pTBins[ibin] + 0.0001);
        int highpTBin = hPhiMassVsPt->GetZaxis()->FindBin(pTBins[ibin + 1] - 0.0001);
        TH1D *hPhi1 = hPhiMassVsPt->ProjectionY(Form("hPhi1_pT_%.1f_%.1f", pTBins[ibin], pTBins[ibin + 1]), lowpTBin, highpTBin, -1, -1, "E");
        SetHistoQA(hPhi1);
        cPhiVsPt->cd(ibin + 1);
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
        fitFcn->SetParameter(1, 1.019);         // mass peak
        fitFcn->SetParLimits(1, 1.006, 1.0316); // mass peak
        fitFcn->SetParameter(2, 0.0012);        //  Gaussian width (Detector resolution)
        fitFcn->SetParLimits(2, 0.001, 0.009);  // Gaussian width.
        // fitFcn->SetParameter(3, 0.0042);   //lorentzian width (Resonance width)
        fitFcn->FixParameter(3, 0.0042); // lorentzian width

        TFitResultPtr r = hPhi1->Fit("fitfunc", "REBMS");
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

        fitMass.push_back({fitFcnSig->GetParameter(1), fitFcnSig->GetParError(1)});
        fitResolution.push_back({fitFcnSig->GetParameter(2) * 1000, fitFcnSig->GetParError(2) * 1000});
        fitWidth.push_back({fitFcnSig->GetParameter(3), fitFcnSig->GetParError(3)});
    }
    cPhiVsPt->SaveAs(savepath + "/PhiMassVsPt.pdf");

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
    gMassVsPt->GetYaxis()->SetRangeUser(1.0194, 1.0202);
    SetGraphErrorStyle(gMassVsPt);
    gMassVsPt->GetYaxis()->SetTitleOffset(2.1);
    gMassVsPt->Draw("AP");
    cMassVsPt->SaveAs(savepath + "/PhiMassVsPt.png");

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
    gResolutionvsPt->SetMinimum(1.1);
    gResolutionvsPt->SetMaximum(2.9);
    gResolutionvsPt->Draw("AP");
    cResolutionVsPt->SaveAs(savepath + "/PhiResolutionVsPt.png");
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
