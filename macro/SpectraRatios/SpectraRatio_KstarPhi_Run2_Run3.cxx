#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"

// Custom style headers
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../spectra/YieldMean.C"

using namespace std;

// Levy-Tsallis fit function for yield interpolation/extraction
Double_t FuncLavy(Double_t *x, Double_t *par)
{
    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] /
                 (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) *
                    (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

TFile *OpenFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (!f || f->IsZombie())
    {
        cout << "Error: File not found or corrupted: " << path << endl;
        return nullptr;
    }
    return f;
}

// Helper to get particle uncertainty across bin matching/overlapping
double getParticleErr(TH1F *hErr, double xLow, double xHigh)
{
    const double eps = 1e-6;
    const int nBins = hErr->GetNbinsX();
    const double pMin = hErr->GetBinLowEdge(1);
    const double pMax = hErr->GetBinLowEdge(nBins) + hErr->GetBinWidth(nBins);

    if (xLow < pMin - eps || xHigh > pMax + eps)
    {
        return hErr->GetBinError(nBins);
    }

    for (int ibin = 1; ibin <= nBins; ++ibin)
    {
        double kaonLow = hErr->GetBinLowEdge(ibin);
        double kaonHigh = kaonLow + hErr->GetBinWidth(ibin);

        if (std::abs(xLow - kaonLow) < eps && std::abs(xHigh - kaonHigh) < eps)
        {
            return hErr->GetBinError(ibin);
        }
    }

    double sumErr = 0.0;
    int nOverlap = 0;
    for (int ibin = 1; ibin <= nBins; ++ibin)
    {
        double kaonLow = hErr->GetBinLowEdge(ibin);
        double kaonHigh = kaonLow + hErr->GetBinWidth(ibin);

        if (kaonHigh > xLow + eps && kaonLow < xHigh - eps)
        {
            sumErr += hErr->GetBinError(ibin);
            ++nOverlap;
        }
    }

    if (nOverlap > 0)
        return sumErr / nOverlap;
    return hErr->GetBinError(nBins);
}

// Helper to get particle yield (either exact bin content or Levy fit evaluation)
double getParticleYield(TH1F *hParticle, TF1 *fitFcn, double xLow, double xHigh)
{
    const double eps = 1e-6;
    const int nBins = hParticle->GetNbinsX();

    for (int ibin = 1; ibin <= nBins; ++ibin)
    {
        double pLow = hParticle->GetBinLowEdge(ibin);
        double pHigh = pLow + hParticle->GetBinWidth(ibin);

        if (std::abs(xLow - pLow) < eps && std::abs(xHigh - pHigh) < eps)
        {
            return hParticle->GetBinContent(ibin);
        }
    }
    double xCenter = 0.5 * (xLow + xHigh);
    return fitFcn->Eval(xCenter);
}

// Function to compute ratio (K* / Phi) and uncertainties
void CalculateRatio(TH1F *hKstarStat, TH1F *hKstarSys, TH1F *hPhiStat, TH1F *hPhiSys,
                    TF1 *fitPhi, TH1F *hRatioStat, TH1F *hRatioSys)
{
    int nBins = hKstarStat->GetNbinsX();
    for (int ibin = 1; ibin <= nBins; ++ibin)
    {
        double kstarYield = hKstarStat->GetBinContent(ibin);
        double kstarStatErr = hKstarStat->GetBinError(ibin);
        double kstarSysErr = hKstarSys->GetBinError(ibin);

        double xLow = hKstarStat->GetBinLowEdge(ibin);
        double xHigh = xLow + hKstarStat->GetBinWidth(ibin);

        double phiYield = getParticleYield(hPhiStat, fitPhi, xLow, xHigh);

        if (kstarYield <= 0 || phiYield <= 0)
        {
            hRatioStat->SetBinContent(ibin, 0);
            hRatioStat->SetBinError(ibin, 0);
            hRatioSys->SetBinContent(ibin, 0);
            hRatioSys->SetBinError(ibin, 0);
            continue;
        }

        double ratio = kstarYield / phiYield;
        double phiStatErr = getParticleErr(hPhiStat, xLow, xHigh);
        double phiSysErr = getParticleErr(hPhiSys, xLow, xHigh);

        double relKstarStat = kstarStatErr / kstarYield;
        double relPhiStat = phiStatErr / phiYield;
        double ratioStatErr = ratio * std::sqrt(relKstarStat * relKstarStat + relPhiStat * relPhiStat);

        double relKstarSys = kstarSysErr / kstarYield;
        double relPhiSys = phiSysErr / phiYield;
        double ratioSysErr = ratio * std::sqrt(relKstarSys * relKstarSys + relPhiSys * relPhiSys);

        hRatioStat->SetBinContent(ibin, ratio);
        hRatioStat->SetBinError(ibin, ratioStatErr);
        hRatioSys->SetBinContent(ibin, ratio);
        hRatioSys->SetBinError(ibin, ratioSysErr);
    }
}

void SpectraRatio_KstarPhi_Run2_Run3()
{
    // =========================================================================
    // SWITCH FOR EPOS MODEL PREDICTIONS
    // =========================================================================
    bool plotEPOS = true; // Set to true to draw EPOS lines/bands, false to hide

    gStyle->SetOptStat(0);

    // =========================================================================
    // 1. RUN 3 DATA LOADING & RATIO CALCULATION
    // =========================================================================
    TFile *fKstarRun3 = OpenFile("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/KstarSpectra.root");
    TFile *fPhiRun3 = OpenFile("PhiResults/PhiAnalysisResults_pprun3.root");

    TH1F *hKstarRun3_statHM = (TH1F *)fKstarRun3->Get("SpectraStat_0_1");
    TH1F *hKstarRun3_sysHM = (TH1F *)fKstarRun3->Get("SpectraSys_0_1");
    TH1F *hKstarRun3_statLM = (TH1F *)fKstarRun3->Get("SpectraStat_70_100");
    TH1F *hKstarRun3_sysLM = (TH1F *)fKstarRun3->Get("SpectraSys_70_100");

    TH1F *hPhiRun3_statHM = (TH1F *)fPhiRun3->Get("Spectra/01/hStat");
    TH1F *hPhiRun3_sysHM = (TH1F *)fPhiRun3->Get("Spectra/01/hSys");
    TH1F *hPhiRun3_statLM = (TH1F *)fPhiRun3->Get("Spectra/10/hStat");
    TH1F *hPhiRun3_sysLM = (TH1F *)fPhiRun3->Get("Spectra/10/hSys");

    TF1 *fitPhiRun3_HM = new TF1("fitPhiRun3_HM", FuncLavy, 0.0, 10.0, 4);
    fitPhiRun3_HM->SetParameters(6.6, 0.2, 1.019, 0.18);
    YieldMean(hPhiRun3_statHM, hPhiRun3_sysHM, fitPhiRun3_HM, 0.0, 10.0, 0.001, 0.5, "RI0+", "log.root", 0.0, 10.0);

    TF1 *fitPhiRun3_LM = new TF1("fitPhiRun3_LM", FuncLavy, 0.0, 10.0, 4);
    fitPhiRun3_LM->SetParameters(6.6, 0.2, 1.019, 0.18);
    YieldMean(hPhiRun3_statLM, hPhiRun3_sysLM, fitPhiRun3_LM, 0.0, 10.0, 0.001, 0.5, "RI0+", "log.root", 0.0, 10.0);

    TH1F *hRatioRun3_HM_stat = (TH1F *)hKstarRun3_statHM->Clone("hRatioRun3_HM_stat");
    TH1F *hRatioRun3_HM_sys = (TH1F *)hKstarRun3_sysHM->Clone("hRatioRun3_HM_sys");
    TH1F *hRatioRun3_LM_stat = (TH1F *)hKstarRun3_statLM->Clone("hRatioRun3_LM_stat");
    TH1F *hRatioRun3_LM_sys = (TH1F *)hKstarRun3_sysLM->Clone("hRatioRun3_LM_sys");

    CalculateRatio(hKstarRun3_statHM, hKstarRun3_sysHM, hPhiRun3_statHM, hPhiRun3_sysHM, fitPhiRun3_HM, hRatioRun3_HM_stat, hRatioRun3_HM_sys);
    CalculateRatio(hKstarRun3_statLM, hKstarRun3_sysLM, hPhiRun3_statLM, hPhiRun3_sysLM, fitPhiRun3_LM, hRatioRun3_LM_stat, hRatioRun3_LM_sys);

    // =========================================================================
    // 2. RUN 2 DATA LOADING & RATIO CALCULATION
    // =========================================================================
    TFile *fRun2 = OpenFile("/home/sawan/Storage/check_k892/macro/SpectraRatios/HEP_data/pp13TeV_INELgt0.root");

    TH1F *hKstarRun2_HM = (TH1F *)fRun2->Get("Table 2/Hist1D_y1");
    TH1F *hKstarRun2_statHM = (TH1F *)fRun2->Get("Table 2/Hist1D_y1_e1");
    TH1F *hKstarRun2_sysHM = (TH1F *)fRun2->Get("Table 2/Hist1D_y1_e2");

    TH1F *hKstarRun2_LM = (TH1F *)fRun2->Get("Table 9/Hist1D_y1");
    TH1F *hKstarRun2_statLM = (TH1F *)fRun2->Get("Table 9/Hist1D_y1_e1");
    TH1F *hKstarRun2_sysLM = (TH1F *)fRun2->Get("Table 9/Hist1D_y1_e2");

    TH1F *hPhiRun2_HM = (TH1F *)fRun2->Get("Table 20/Hist1D_y1");
    TH1F *hPhiRun2_statHM = (TH1F *)fRun2->Get("Table 20/Hist1D_y1_e1");
    TH1F *hPhiRun2_sysHM = (TH1F *)fRun2->Get("Table 20/Hist1D_y1_e2");

    TH1F *hPhiRun2_LM = (TH1F *)fRun2->Get("Table 28/Hist1D_y1");
    TH1F *hPhiRun2_statLM = (TH1F *)fRun2->Get("Table 28/Hist1D_y1_e1");
    TH1F *hPhiRun2_sysLM = (TH1F *)fRun2->Get("Table 28/Hist1D_y1_e2");

    for (int i = 1; i <= hKstarRun2_HM->GetNbinsX(); i++)
    {
        hKstarRun2_statHM->SetBinContent(i, hKstarRun2_HM->GetBinContent(i));
        hKstarRun2_sysHM->SetBinContent(i, hKstarRun2_HM->GetBinContent(i));
    }
    for (int i = 1; i <= hKstarRun2_LM->GetNbinsX(); i++)
    {
        hKstarRun2_statLM->SetBinContent(i, hKstarRun2_LM->GetBinContent(i));
        hKstarRun2_sysLM->SetBinContent(i, hKstarRun2_LM->GetBinContent(i));
    }
    for (int i = 1; i <= hPhiRun2_HM->GetNbinsX(); i++)
    {
        hPhiRun2_statHM->SetBinContent(i, hPhiRun2_HM->GetBinContent(i));
        hPhiRun2_sysHM->SetBinContent(i, hPhiRun2_HM->GetBinContent(i));
    }
    for (int i = 1; i <= hPhiRun2_LM->GetNbinsX(); i++)
    {
        hPhiRun2_statLM->SetBinContent(i, hPhiRun2_LM->GetBinContent(i));
        hPhiRun2_sysLM->SetBinContent(i, hPhiRun2_LM->GetBinContent(i));
    }

    TF1 *fitPhiRun2_HM = new TF1("fitPhiRun2_HM", FuncLavy, 0.5, 6.0, 4);
    fitPhiRun2_HM->SetParameters(6.6, 0.2, 1.019, 0.18);
    YieldMean(hPhiRun2_statHM, hPhiRun2_sysHM, fitPhiRun2_HM, 0.0, 6.0, 0.001, 0.5, "RI0+", "log.root", 0.0, 6.0);

    TF1 *fitPhiRun2_LM = new TF1("fitPhiRun2_LM", FuncLavy, 0.0, 6.0, 4);
    fitPhiRun2_LM->SetParameters(6.6, 0.2, 1.019, 0.18);
    YieldMean(hPhiRun2_statLM, hPhiRun2_sysLM, fitPhiRun2_LM, 0.0, 6.0, 0.001, 0.5, "RI0+", "log.root", 0.0, 6.0);

    TH1F *hRatioRun2_HM_stat = (TH1F *)hKstarRun2_statHM->Clone("hRatioRun2_HM_stat");
    TH1F *hRatioRun2_HM_sys = (TH1F *)hKstarRun2_sysHM->Clone("hRatioRun2_HM_sys");
    TH1F *hRatioRun2_LM_stat = (TH1F *)hKstarRun2_statLM->Clone("hRatioRun2_LM_stat");
    TH1F *hRatioRun2_LM_sys = (TH1F *)hKstarRun2_sysLM->Clone("hRatioRun2_LM_sys");

    SetHistoQA(hRatioRun2_HM_stat);
    SetHistoQA(hRatioRun2_HM_sys);
    SetHistoQA(hRatioRun2_LM_stat);
    SetHistoQA(hRatioRun2_LM_sys);

    CalculateRatio(hKstarRun2_statHM, hKstarRun2_sysHM, hPhiRun2_statHM, hPhiRun2_sysHM, fitPhiRun2_HM, hRatioRun2_HM_stat, hRatioRun2_HM_sys);
    CalculateRatio(hKstarRun2_statLM, hKstarRun2_sysLM, hPhiRun2_statLM, hPhiRun2_sysLM, fitPhiRun2_LM, hRatioRun2_LM_stat, hRatioRun2_LM_sys);

    // =========================================================================
    // 3. EPOS PREPARATION (UrQMD ON vs OFF)
    // =========================================================================
    TGraph *gHM_EPOS_OFF = nullptr, *gLM_EPOS_OFF = nullptr;
    TGraphErrors *bHM_EPOS_OFF = nullptr, *bLM_EPOS_OFF = nullptr;
    TGraph *gHM_EPOS_ON = nullptr, *gLM_EPOS_ON = nullptr;
    TGraphErrors *bHM_EPOS_ON = nullptr, *bLM_EPOS_ON = nullptr;

    if (plotEPOS)
    {
        TFile *fEPOS_spectra = OpenFile("/home/sawan/Storage/EPOS4_Resonance/yield_outfile_light_oct2026_fixed.root");

        TH1F *hEPOS_OFF_KstarHM = (TH1F *)fEPOS_spectra->Get("Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent9");
        TH1F *hEPOS_OFF_KstarLM = (TH1F *)fEPOS_spectra->Get("Kstar_Mult_Spectra_OFF/hKstar_Pt_OFF_Cent0");
        TH1F *hEPOS_ON_KstarHM = (TH1F *)fEPOS_spectra->Get("Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent9");
        TH1F *hEPOS_ON_KstarLM = (TH1F *)fEPOS_spectra->Get("Kstar_Mult_Spectra_ON_reg_res/hKstar_Pt_ON_reg_res_Cent0");

        TH1F *hEPOS_OFF_PhiHM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent9");
        TH1F *hEPOS_OFF_PhiLM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_OFF/hPhi_Pt_OFF_Cent0");
        TH1F *hEPOS_ON_PhiHM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent9");
        TH1F *hEPOS_ON_PhiLM = (TH1F *)fEPOS_spectra->Get("Phi_Mult_Spectra_ON_reg_res/hPhi_Pt_ON_reg_res_Cent0");

        // Correct branching ratio for UrQMD ON phi spectra
        double phiBR = 0.489;
        hEPOS_ON_PhiHM->Scale(1.0 / phiBR);
        hEPOS_ON_PhiLM->Scale(1.0 / phiBR);

        double mergedBins[] = {0.3, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 4.0, 5.0, 10.0};
        const int nMergedBins = sizeof(mergedBins) / sizeof(double) - 1;

        auto getIntegral = [](TH1F *h, double xLow, double xHigh, double &yield, double &error)
        {
            int binLow = h->GetXaxis()->FindBin(xLow + 1e-6);
            int binHigh = h->GetXaxis()->FindBin(xHigh - 1e-6);
            yield = h->IntegralAndError(binLow, binHigh, error, "width");
        };

        TH1F *hEPOS_RatioHM_OFF = new TH1F("hEPOS_RatioHM_OFF", ";#it{p}_{T};Ratio", nMergedBins, mergedBins);
        TH1F *hEPOS_RatioLM_OFF = new TH1F("hEPOS_RatioLM_OFF", ";#it{p}_{T};Ratio", nMergedBins, mergedBins);
        TH1F *hEPOS_RatioHM_ON = new TH1F("hEPOS_RatioHM_ON", ";#it{p}_{T};Ratio", nMergedBins, mergedBins);
        TH1F *hEPOS_RatioLM_ON = new TH1F("hEPOS_RatioLM_ON", ";#it{p}_{T};Ratio", nMergedBins, mergedBins);

        for (int i = 0; i < nMergedBins; ++i)
        {
            double xLow = mergedBins[i], xHigh = mergedBins[i + 1];
            double kHM, kErrHM, phiHM, phiErrHM;
            double kLM, kErrLM, phiLM, phiErrLM;

            // OFF
            getIntegral(hEPOS_OFF_KstarHM, xLow, xHigh, kHM, kErrHM);
            getIntegral(hEPOS_OFF_PhiHM, xLow, xHigh, phiHM, phiErrHM);
            if (kHM > 0 && phiHM > 0)
            {
                double r = 0.5 * kHM / phiHM;
                double err = r * sqrt(pow(kErrHM / kHM, 2) + pow(phiErrHM / phiHM, 2)) / 2.0;
                hEPOS_RatioHM_OFF->SetBinContent(i + 1, r);
                hEPOS_RatioHM_OFF->SetBinError(i + 1, err);
            }

            getIntegral(hEPOS_OFF_KstarLM, xLow, xHigh, kLM, kErrLM);
            getIntegral(hEPOS_OFF_PhiLM, xLow, xHigh, phiLM, phiErrLM);
            if (kLM > 0 && phiLM > 0)
            {
                double r = 0.5 * kLM / phiLM;
                double err = r * sqrt(pow(kErrLM / kLM, 2) + pow(phiErrLM / phiLM, 2)) / 2.0;
                hEPOS_RatioLM_OFF->SetBinContent(i + 1, r);
                hEPOS_RatioLM_OFF->SetBinError(i + 1, err);
            }

            // ON
            getIntegral(hEPOS_ON_KstarHM, xLow, xHigh, kHM, kErrHM);
            getIntegral(hEPOS_ON_PhiHM, xLow, xHigh, phiHM, phiErrHM);
            if (kHM > 0 && phiHM > 0)
            {
                double r = 0.5 * kHM / phiHM;
                double err = r * sqrt(pow(0.5 * kErrHM / kHM, 2) + pow(phiErrHM / phiHM, 2)) / 2.0;
                hEPOS_RatioHM_ON->SetBinContent(i + 1, r);
                hEPOS_RatioHM_ON->SetBinError(i + 1, err);
            }

            getIntegral(hEPOS_ON_KstarLM, xLow, xHigh, kLM, kErrLM);
            getIntegral(hEPOS_ON_PhiLM, xLow, xHigh, phiLM, phiErrLM);
            if (kLM > 0 && phiLM > 0)
            {
                double r = 0.5 * kLM / phiLM;
                double err = r * sqrt(pow(0.5 * kErrLM / kLM, 2) + pow(phiErrLM / phiLM, 2)) / 2.0;
                hEPOS_RatioLM_ON->SetBinContent(i + 1, r);
                hEPOS_RatioLM_ON->SetBinError(i + 1, err);
            }
        }

        hEPOS_RatioHM_OFF->Smooth(2);
        hEPOS_RatioLM_OFF->Smooth(2);
        hEPOS_RatioHM_ON->Smooth(2);
        hEPOS_RatioLM_ON->Smooth(2);

        auto makeCurve = [](TH1 *h, Color_t col, Style_t lstyle, TGraph *&gLine, TGraphErrors *&gBand)
        {
            gLine = new TGraph();
            gBand = new TGraphErrors();
            int n = h->GetNbinsX(), first = -1, last = -1;
            for (int i = 1; i <= n; ++i)
                if (h->GetBinContent(i) > 0)
                {
                    if (first < 0)
                        first = i;
                    last = i;
                }
            if (first < 0)
                return;

            auto add = [&](double x, int ibin)
            {
                gLine->SetPoint(gLine->GetN(), x, h->GetBinContent(ibin));
                int k = gBand->GetN();
                gBand->SetPoint(k, x, h->GetBinContent(ibin));
                gBand->SetPointError(k, 0, h->GetBinError(ibin));
            };

            add(h->GetXaxis()->GetXmin(), first);
            for (int i = first; i <= last; ++i)
                if (h->GetBinContent(i) > 0)
                    add(h->GetBinCenter(i), i);
            add(h->GetXaxis()->GetXmax(), last);

            gLine->SetLineColor(col);
            gLine->SetLineWidth(2);
            gLine->SetLineStyle(lstyle);
            gBand->SetFillColorAlpha(col, 0.2);
            gBand->SetLineColor(col);
        };

        makeCurve(hEPOS_RatioHM_OFF, kGray + 2, 1, gHM_EPOS_OFF, bHM_EPOS_OFF);
        makeCurve(hEPOS_RatioLM_OFF, kGray + 2, 1, gLM_EPOS_OFF, bLM_EPOS_OFF);
        makeCurve(hEPOS_RatioHM_ON, kGreen + 2, 7, gHM_EPOS_ON, bHM_EPOS_ON);
        makeCurve(hEPOS_RatioLM_ON, kGreen + 2, 7, gLM_EPOS_ON, bLM_EPOS_ON);
    }

    // =========================================================================
    // 4. CANVAS & TWO-PANEL CANVAS PLOTTING
    // =========================================================================
    TCanvas *cCompare = new TCanvas("cCompare", "Kstar/Phi Ratio Comparison Run2 vs Run3", 1200, 600);
    cCompare->Divide(2, 1, 0.001, 0.001);

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);

    // -------------------------------------------------------------------------
    // Left Panel: 0-1% Multiplicity Comparison
    // -------------------------------------------------------------------------
    cCompare->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.06);
    gPad->SetBottomMargin(0.13);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    hRatioRun3_HM_stat->SetTitle(";#it{p}_{T} (GeV/#it{c}); K^{*0} / #Phi");
    hRatioRun3_HM_stat->GetXaxis()->SetRangeUser(0.5, 6.0);
    hRatioRun3_HM_stat->GetYaxis()->SetRangeUser(0.0, 5.0);
    hRatioRun3_HM_stat->GetXaxis()->SetTitleSize(0.045);
    hRatioRun3_HM_stat->GetYaxis()->SetTitleSize(0.045);
    hRatioRun3_HM_stat->GetXaxis()->SetLabelSize(0.04);
    hRatioRun3_HM_stat->GetYaxis()->SetLabelSize(0.04);
    hRatioRun3_HM_stat->GetYaxis()->SetTitleOffset(1.3);

    hRatioRun3_HM_stat->SetMarkerStyle(20);
    hRatioRun3_HM_stat->SetMarkerSize(1.3);
    hRatioRun3_HM_stat->SetMarkerColor(kRed + 1);
    hRatioRun3_HM_stat->SetLineColor(kRed + 1);
    hRatioRun3_HM_stat->Draw("pe");

    hRatioRun3_HM_sys->SetFillStyle(0);
    hRatioRun3_HM_sys->SetMarkerColor(kRed + 1);
    hRatioRun3_HM_sys->SetMarkerStyle(20);
    hRatioRun3_HM_sys->SetMarkerSize(1.3);
    hRatioRun3_HM_sys->SetLineColor(kRed + 1);
    hRatioRun3_HM_sys->Draw("e2 same");

    // Run 2 setup (Blue Marker / Box)
    hRatioRun2_HM_stat->SetMarkerStyle(24);
    hRatioRun2_HM_stat->SetMarkerSize(1.3);
    hRatioRun2_HM_stat->SetMarkerColor(kBlue + 1);
    hRatioRun2_HM_stat->SetLineColor(kBlue + 1);
    hRatioRun2_HM_stat->Draw("pe same");

    hRatioRun2_HM_sys->SetFillStyle(0);
    hRatioRun2_HM_sys->SetMarkerColor(kBlue + 1);
    hRatioRun2_HM_sys->SetMarkerStyle(24);
    hRatioRun2_HM_sys->SetMarkerSize(1.3);
    hRatioRun2_HM_sys->SetLineColor(kBlue + 1);
    hRatioRun2_HM_sys->Draw("e2 same");

    if (plotEPOS)
    {
        bHM_EPOS_OFF->Draw("3 same");
        gHM_EPOS_OFF->Draw("L same");
        bHM_EPOS_ON->Draw("3 same");
        gHM_EPOS_ON->Draw("L same");
    }

    TLegend *legHM = new TLegend(0.48, 0.72, 0.90, 0.90);
    legHM->SetBorderSize(0);
    legHM->SetFillStyle(0);
    legHM->SetTextSize(0.038);
    legHM->SetHeader("FT0M / V0M: 0-1%");
    legHM->AddEntry(hRatioRun3_HM_stat, "Run 3 (pp #sqrt{#it{s}} = 13.6 TeV)", "pe");
    legHM->AddEntry(hRatioRun2_HM_stat, "Run 2 (pp #sqrt{#it{s}} = 13 TeV)", "pe");
    if (plotEPOS)
    {
        legHM->AddEntry(gHM_EPOS_OFF, "EPOS4 UrQMD OFF", "l");
        legHM->AddEntry(gHM_EPOS_ON, "EPOS4 UrQMD ON", "l");
    }
    legHM->Draw();

    lat.SetTextSize(0.045);
    lat.DrawLatex(0.20, 0.88, "ALICE");
    lat.SetTextSize(0.038);
    lat.DrawLatex(0.20, 0.82, "|y| < 0.5");

    // -------------------------------------------------------------------------
    // Right Panel: 70-100% Multiplicity Comparison
    // -------------------------------------------------------------------------
    cCompare->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.06);
    gPad->SetBottomMargin(0.13);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    hRatioRun3_LM_stat->SetTitle(";#it{p}_{T} (GeV/#it{c}); K^{*0} / #Phi");
    hRatioRun3_LM_stat->GetXaxis()->SetRangeUser(0.5, 6.0);
    hRatioRun3_LM_stat->GetYaxis()->SetRangeUser(0.0, 5.0);
    hRatioRun3_LM_stat->GetXaxis()->SetTitleSize(0.045);
    hRatioRun3_LM_stat->GetYaxis()->SetTitleSize(0.045);
    hRatioRun3_LM_stat->GetXaxis()->SetLabelSize(0.04);
    hRatioRun3_LM_stat->GetYaxis()->SetLabelSize(0.04);
    hRatioRun3_LM_stat->GetYaxis()->SetTitleOffset(1.3);

    hRatioRun3_LM_stat->SetMarkerStyle(21);
    hRatioRun3_LM_stat->SetMarkerSize(1.3);
    hRatioRun3_LM_stat->SetMarkerColor(kRed + 1);
    hRatioRun3_LM_stat->SetLineColor(kRed + 1);
    hRatioRun3_LM_stat->Draw("pe");

    hRatioRun3_LM_sys->SetFillStyle(0);
    hRatioRun3_LM_sys->SetLineColor(kRed + 1);
    hRatioRun3_LM_sys->SetMarkerColor(kRed + 1);
    hRatioRun3_LM_sys->SetMarkerStyle(21);
    hRatioRun3_LM_sys->SetMarkerSize(1.3);
    hRatioRun3_LM_sys->Draw("e2 same");

    // Run 2 setup (Blue Marker / Box)
    hRatioRun2_LM_stat->SetMarkerStyle(25);
    hRatioRun2_LM_stat->SetMarkerSize(1.3);
    hRatioRun2_LM_stat->SetMarkerColor(kBlue + 1);
    hRatioRun2_LM_stat->SetLineColor(kBlue + 1);
    hRatioRun2_LM_stat->Draw("pe same");

    hRatioRun2_LM_sys->SetFillStyle(0);
    hRatioRun2_LM_sys->SetLineColor(kBlue + 1);
    hRatioRun2_LM_sys->SetMarkerColor(kBlue + 1);
    hRatioRun2_LM_sys->SetMarkerStyle(25);
    hRatioRun2_LM_sys->SetMarkerSize(1.3);
    hRatioRun2_LM_sys->Draw("e2 same");

    if (plotEPOS)
    {
        bLM_EPOS_OFF->Draw("3 same");
        gLM_EPOS_OFF->Draw("L same");
        bLM_EPOS_ON->Draw("3 same");
        gLM_EPOS_ON->Draw("L same");
    }

    TLegend *legLM = new TLegend(0.48, 0.72, 0.90, 0.90);
    legLM->SetBorderSize(0);
    legLM->SetFillStyle(0);
    legLM->SetTextSize(0.038);
    legLM->SetHeader("FT0M / V0M: 70-100%");
    legLM->AddEntry(hRatioRun3_LM_stat, "Run 3 (pp #sqrt{#it{s}} = 13.6 TeV)", "pe");
    legLM->AddEntry(hRatioRun2_LM_stat, "Run 2 (pp #sqrt{#it{s}} = 13 TeV)", "pe");
    if (plotEPOS)
    {
        legLM->AddEntry(gLM_EPOS_OFF, "EPOS4 UrQMD OFF", "l");
        legLM->AddEntry(gLM_EPOS_ON, "EPOS4 UrQMD ON", "l");
    }
    legLM->Draw();

    lat.SetTextSize(0.045);
    lat.DrawLatex(0.20, 0.88, "ALICE");
    lat.SetTextSize(0.038);
    lat.DrawLatex(0.20, 0.82, "|y| < 0.5");

    cCompare->SaveAs("Plots/pTDiff_Ratios/QA/KstarPhiRatio_Run2_vs_Run3_ComparisonWithModel.png");
}