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

void phiFitParameters()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    // TString suffix = "25_aiam";
    // TString suffix = "25_aiamShifted_opti10";
    // TString suffix = "26";
    TString suffix = "26ai";
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/PhiInvMass";

    ////=====New===========
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResults26_latest.root"); // 2026
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/LHC25/AnalysisResults25_aiam.root"); // 2025 ai+am
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti9/LHC25/AnalysisResults25_aiamShifted.root"); // 2025 ai+am Opti9 (Ka momentum shifted)
    // TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti10/LHC25/AnalysisResults25_aiamShifted.root"); // 2025 ai+am Opti10 (Ka momentum shifted)

    ////=============Separate files====================
    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/All/AnalysisResults_LHC26ai.root");

    // //=========================================
    // //=======Phi Inv mass fit (1D)=============
    // //=========================================

    TH3F *hPhiMassVsPt = GetHisto<TH3F>(fInput, "doublephimeson/hPhiMass");
    // TH3F *hPhiMassVsPt = GetHisto<TH3F>(fInput, "doublephimeson/hPhiMassShifted");
    TH2F *hPhiPhiMass = (TH2F *)hPhiMassVsPt->Project3D("yx");
    TH2F *hPhiMassVsPt2D = GetHisto<TH2F>(fInput, "doublephimeson/hPhiMassVsPt");
    // TH2F *hPhiMassVsPt2D = GetHisto<TH2F>(fInput, "doublephimeson/hPhiMassVsPtShifted");

    TCanvas *cPhiVsPt = new TCanvas("cPhiVsPt", "Phi Mass vs Pt", 1280, 720);
    cPhiVsPt->Divide(4, 4);
    double pTBins[15] = {0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(22);
    latex->SetTextSize(0.06);
    vector<pair<double, double>> fitMass, fitResolution, fitWidth;
    TGraphErrors *gPurity = new TGraphErrors();
    for (int ibin = 0; ibin < 14; ibin++)
    {

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

        // For voigtian distribution
        fitFcn->SetParameter(0, 5000);          // yield
        fitFcn->SetParLimits(0, 0, 1e6);        // yield
        fitFcn->SetParameter(1, 1.0198);        // mass peak
        fitFcn->SetParLimits(1, 1.0175, 1.022);  // mass peak
        fitFcn->SetParameter(2, 0.0012);        //  Gaussian width (Detector resolution)
        fitFcn->SetParLimits(2, 0.0008, 0.006); // Gaussian width.
        // fitFcn->SetParameter(3, 0.0042);   //lorentzian width (Resonance width)
        fitFcn->FixParameter(3, 0.0042); // lorentzian width

        // if (ibin == 3)
        //     fitFcn->SetParLimits(0, 0, 1e3);


        fitFcn->SetParameter(4, 4e4);        // Pol2 p0
        fitFcn->SetParLimits(4, 1e2, 3e6);   // Pol2 p0
        fitFcn->SetParameter(5, 4.1e5);      // Pol2 p1
        fitFcn->SetParLimits(5, 1e2, 3e7);   // Pol2 p1
        fitFcn->SetParameter(6, -5.0e6);     // Pol2 p2
        fitFcn->SetParLimits(6, -1e9, -1e2); // Pol2 p2

        // if (ibin == 12 || ibin == 13)
        // {
        //     hPhi1->Rebin(2);
        // }


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
    // cPhiVsPt->SaveAs(savepath + "/PhiInvMassVsPt_" + suffix + ".png");

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
    // cPurityVsPt->SaveAs(savepath + "/PhiPurityVsPt_" + suffix + ".png");

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
    TLine *linePDGMass = new TLine(pTBins[0], 1.019461, pTBins[14], 1.019461);
    linePDGMass->SetLineColor(kRed);
    linePDGMass->SetLineStyle(2);
    linePDGMass->SetLineWidth(2);
    linePDGMass->Draw("same");
    latex->SetTextSize(0.035);
    latex->DrawLatex(0.25, 0.88, "Fit function: Voigtian + Pol2");
    latex->DrawLatex(0.25, 0.81, "Purity window: #it{M}_{#Phi} #pm 0.005 GeV/#it{c}^{2}");
    // cMassVsPt->SaveAs(savepath + "/PhiMassVsPt_" + suffix + ".png");

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
    // cResolutionVsPt->SaveAs(savepath + "/PhiResolutionVsPt_" + suffix + ".png");

    // Phi mass correlation plot
    SetHistoQA(hPhiPhiMass);
    TCanvas *cPhiPhiMassCorr = new TCanvas("cPhiPhiMassCorr", "Phi Mass vs Pt", 720, 720);
    SetCanvasStyle(cPhiPhiMassCorr, 0.15, 0.15, 0.05, 0.15);
    SetHistoQA(hPhiPhiMass);
    hPhiPhiMass->GetXaxis()->SetTitle("M_{#phi1} (GeV/#it{c})");
    hPhiPhiMass->GetYaxis()->SetTitle("M_{#phi2} (GeV/#it{c})");
    hPhiPhiMass->GetYaxis()->SetTitleOffset(1.4);
    hPhiPhiMass->GetYaxis()->SetNdivisions(505);
    hPhiPhiMass->GetXaxis()->SetNdivisions(505);
    hPhiPhiMass->GetZaxis()->SetMaxDigits(3);
    hPhiPhiMass->Draw("colz");
    TEllipse *circle = new TEllipse(1.0198, 1.0198, 0.005);
    circle->SetFillStyle(0); // no fill
    circle->SetLineColor(kRed);
    circle->SetLineWidth(3);
    circle->Draw("same");
    latex->SetTextSize(0.035);
    latex->SetTextColor(kWhite);
    latex->DrawLatex(0.2, 0.91, "Red Circle: #DeltaM < 0.005 GeV/#it{c}^{2}");
    latex->DrawLatex(0.2, 0.85, "#it{p}_{T}^{#phi#phi} > 9 GeV/#it{c}");
    // cPhiPhiMassCorr->SaveAs(savepath + "/PhiMassCorrelation_" + suffix + ".png");

    TFile *fPhiParams = new TFile(savepath + "/PhiParams" + suffix + ".root", "recreate");
    gMassVsPt->Write("gMassVsPt");
    gResolutionvsPt->Write("gResolutionVsPt");
    gPurity->Write("gPurity");
    hPhiPhiMass->Write("hPhiPhiMassCorrelation");
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
