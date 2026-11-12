#include <iostream>
#include <iomanip>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "../src/style.h"
using namespace std;

//--------Important info---------------------//
//=====In data:
// 1. K*0, K*+- is averaged
// 2. Pi,K,p is averaged
//=====In EPOS:
// 1. Only K*0 is averaged.
//===== In Pythia (Local):
// 1. Nothing is averaged. All as summed.

// Important: Do not scale any graph without cloning it first, otherwise the original graph will be scaled and it will affect other plots.
TGraphErrors *SmoothGraph(const TGraphErrors *gr);
TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false, float CorrFactorDen = 1);
TGraphErrors **MakeRatioUncorr(TGraphErrors **numerator, TGraphErrors **denominator, float CorrFactorDen = 1, int nGraphs = 3);
TGraphErrors *GetGraph(TFile *f, const string &name);
TFile *OpenFile(const string &path);
TGraphErrors *MakeDataModelRatio(TGraphErrors *data, TGraphErrors *model);
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);
void CalculateKstarPhiRatioRun2(TGraphErrors *kstar[3], TGraphErrors *phi[3], TGraphErrors *ratio[3]);

TGraphErrors *NormalizeGraphByX(TGraphErrors *gr, const char *name)
{
    if (!gr)
        return nullptr;

    TGraphErrors *grNorm = new TGraphErrors();

    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);

        double ex = gr->GetErrorX(i);
        double ey = gr->GetErrorY(i);

        if (x <= 0)
            continue;

        double yNorm = y / x;
        double eyNorm = ey / x;

        grNorm->SetPoint(grNorm->GetN(), x, yNorm);
        grNorm->SetPointError(grNorm->GetN() - 1, ex, eyNorm);

        cout << name << " Point " << i
             << ": y = " << yNorm
             << ", ey = " << eyNorm << endl;
    }

    grNorm->SetName(name);
    grNorm->SetLineColor(gr->GetLineColor());
    grNorm->SetLineStyle(gr->GetLineStyle());
    grNorm->SetLineWidth(gr->GetLineWidth());
    grNorm->SetMarkerColor(gr->GetMarkerColor());
    grNorm->SetMarkerStyle(gr->GetMarkerStyle());
    grNorm->SetMarkerSize(gr->GetMarkerSize());
    grNorm->SetFillColor(gr->GetFillColor());
    grNorm->SetFillStyle(gr->GetFillStyle());
    return grNorm;
}

TGraphErrors *NormalizeYieldByMB(TGraphErrors *gr, double yieldMB, double yieldMBErr, double multMB, double multMBErr, const char *name)
{
    if (!gr || yieldMB <= 0 || multMB <= 0)
        return nullptr;

    TGraphErrors *grNorm = new TGraphErrors();
    grNorm->SetName(name);

    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        double ex = gr->GetErrorX(i);
        double ey = gr->GetErrorY(i);

        double xNorm = x / multMB;
        double yNorm = y / yieldMB;
        // double exNorm = sqrt(pow(ex / multMB, 2) + pow(x * multMBErr / (multMB * multMB), 2));
        double exNorm = ex / multMB;
        double eyNorm = sqrt(pow(ey / yieldMB, 2) + pow(y * yieldMBErr / (yieldMB * yieldMB), 2));

        int ip = grNorm->GetN();
        grNorm->SetPoint(ip, xNorm, yNorm);
        grNorm->SetPointError(ip, exNorm, eyNorm);
    }

    grNorm->SetLineColor(gr->GetLineColor());
    grNorm->SetLineStyle(gr->GetLineStyle());
    grNorm->SetLineWidth(gr->GetLineWidth());
    grNorm->SetMarkerColor(gr->GetMarkerColor());
    grNorm->SetMarkerStyle(gr->GetMarkerStyle());
    grNorm->SetMarkerSize(gr->GetMarkerSize());
    return grNorm;
}

TGraphAsymmErrors *NormalizeYieldSysByMB(TGraphErrors *gr, TGraphAsymmErrors *grMB, double yieldMB, double multMB, double multMBErr, const char *name)
{
    if (!gr || !grMB || yieldMB <= 0 || multMB <= 0 || grMB->GetN() < 1)
        return nullptr;

    double xMB, yMB;
    grMB->GetPoint(0, xMB, yMB);
    TGraphAsymmErrors *grNorm = new TGraphAsymmErrors();
    grNorm->SetName(name);

    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        double ex = gr->GetErrorX(i);
        double ey = gr->GetErrorY(i);

        double xNorm = x / multMB;
        double yNorm = y / yieldMB;
        // double exNorm = sqrt(pow(ex / multMB, 2) + pow(x * multMBErr / (multMB * multMB), 2));
        double exNorm = ex / multMB;
        double eyLow = sqrt(pow(ey / yieldMB, 2) + pow(y * grMB->GetErrorYlow(0) / (yieldMB * yieldMB), 2));
        double eyHigh = sqrt(pow(ey / yieldMB, 2) + pow(y * grMB->GetErrorYhigh(0) / (yieldMB * yieldMB), 2));

        int ip = grNorm->GetN();
        grNorm->SetPoint(ip, xNorm, yNorm);
        grNorm->SetPointError(ip, exNorm, exNorm, eyLow, eyHigh);
    }

    return grNorm;
}

void setStyle(TGraphErrors *gr, int color, int style)
{
    gr->SetLineColor(color);
    gr->SetLineStyle(style);
}

void ScaleGraph(TGraph *gr, double scale);
TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance = 0.5, int noOfGraphs = 3);
TGraphErrors *DivideByMultModel(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5);
void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax);

void ParticleRatio_QA()
{
    bool isSavePlots = true;
    string KstarPath = "../../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/";
    TFile *fKstar = OpenFile(KstarPath + "Results.root");

    TFile *fPion = OpenFile("PiKp_Run3_Results/Sawan/Pi_results.root");
    TFile *fProton = OpenFile("PiKp_Run3_Results/Sawan/Pr_results.root");
    TFile *fKaon = OpenFile("PiKp_Run3_Results/Sawan/Ka_results.root");

    TFile *fPiKp = OpenFile("ConversionCodes/ppRun3_PiKpHEP.root");
    TFile *fKstarRun2 = OpenFile("ConversionCodes/pp13TeV_KstarPhi.root");
    TFile *fPhiRun2 = OpenFile("ConversionCodes/pp13TeV_KstarPhi.root");
    TFile *fChKstar = OpenFile("ConversionCodes/pp13TeV_ChKstar.root");
    TFile *fMultStrange = OpenFile("ConversionCodes/pp13TeV_MultStrange.root");
    TFile *fXiStarSigma = OpenFile("ConversionCodes/pp13TeV_XiStarSigma.root");
    TFile *fLambda1520Run2 = OpenFile("ConversionCodes/pp13TeV_LambdaStar.root");

    // Available in run3
    TFile *fLambda1520 = OpenFile("LambdaRun3/Sawan/ResultsLambda1520.root");
    TFile *fRho = OpenFile("Rho_Run3_Results/Sawan/ResultsRho.root");
    TFile *fKshortRun3 = OpenFile("K0s_Run3_Results/Sawan/ResultsK0s.root");
    TFile *fPhiRun3 = OpenFile("PhiResults/Sawan/ResultsPhi.root");

    // ALICE Run2 results
    TFile *fpp13TeV = OpenFile("ConversionCodes/pp13TeV_KstarPhi.root");
    TFile *fpp7TeV = OpenFile("ConversionCodes/pp7TeVALICE.root");

    TGraphErrors *gMPtKstar[3], *gMYieldKstar[3], *gMPtPion[3], *gMPtProton[3], *gMPtKaon[3], *gMPtKstarRun2[3], *gMPtPhiRun2[3], *gMPtChKstarRun2[3], *gMPtKshortRun2[3], *gMPtLambdaRun2[3], *gMPtXiRun2[3], *gMPtOmegaRun2[3], *gMPtXiStarRun2[3], *gMPtSigmaRun2[3], *gMPtLambda1520[3], *gMPtRho[3], *gMPtPhiRun3[3];
    TGraphErrors *gMYieldPion[3], *gMYieldProton[3], *gMYieldKaon[3], *gMYieldKstarRun2[3], *gMYieldPhiRun2[3], *gMYieldChKstarRun2[3], *gMYieldKshortRun2[3], *gMYieldLambdaRun2[3], *gMYieldXiRun2[3], *gMYieldOmegaRun2[3], *gMYieldXiStarRun2[3], *gMYieldSigmaRun2[3], *gMYieldLambda1520[3], *gMYieldRho[3], *gKstar_MeanpT_13TeV[3], *gKstar_Yield_13TeV[3], *gMYieldPhiRun3[3];

    TGraphErrors *gMPtKshortRun3[3], *gMYieldKshortRun3[3];

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

        // Minimum-bias reference graphs saved in the Kstar Results.root file.

        gMPtPhiRun3[i] = GetGraph(fPhiRun3, Form("gMeanpTRun3%s", suffix.c_str()));
        gMYieldPhiRun3[i] = GetGraph(fPhiRun3, Form("gMeanYieldRun3%s", suffix.c_str()));

        gKstar_MeanpT_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_MeanpT%s", suffix.c_str()));
        gKstar_Yield_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_MeanYield%s", suffix.c_str()));

        gMPtPion[i] = GetGraph(fPion, Form("gMeanpTRun3%s", suffix.c_str()));
        gMPtKaon[i] = GetGraph(fKaon, Form("gMeanpTRun3%s", suffix.c_str()));
        gMPtProton[i] = GetGraph(fProton, Form("gMeanpTRun3%s", suffix.c_str()));

        // gMPtPion[i] = GetGraph(fPiKp, Form("gPion_MeanpT%s", suffix.c_str()));
        // gMPtKaon[i] = GetGraph(fPiKp, Form("gKaon_MeanpT%s", suffix.c_str()));
        // gMPtProton[i] = GetGraph(fPiKp, Form("gProton_MeanpT%s", suffix.c_str()));

        gMPtKstarRun2[i] = GetGraph(fKstarRun2, Form("gKstar_MeanpT%s", suffix.c_str()));
        gMPtPhiRun2[i] = GetGraph(fPhiRun2, Form("gPhi_MeanpT%s", suffix.c_str()));
        gMPtChKstarRun2[i] = GetGraph(fChKstar, Form("gChKstar_MeanpT%s", suffix.c_str()));
        gMPtKshortRun2[i] = GetGraph(fMultStrange, Form("gKshort_MeanpT%s", suffix.c_str()));
        gMPtLambdaRun2[i] = GetGraph(fMultStrange, Form("gLambda_MeanpT%s", suffix.c_str()));
        gMPtXiRun2[i] = GetGraph(fMultStrange, Form("gXi_MeanpT%s", suffix.c_str()));
        gMPtOmegaRun2[i] = GetGraph(fMultStrange, Form("gOmega_MeanpT%s", suffix.c_str()));
        gMPtXiStarRun2[i] = GetGraph(fXiStarSigma, Form("gXi1530_MeanpT%s", suffix.c_str()));
        gMPtSigmaRun2[i] = GetGraph(fXiStarSigma, Form("gSigma1385_MeanpT%s", suffix.c_str()));
        gMPtLambda1520[i] = GetGraph(fLambda1520Run2, Form("gLambda_MeanpT%s", suffix.c_str()));

        gMYieldPion[i] = GetGraph(fPion, Form("gMeanYieldRun3%s", suffix.c_str()));
        gMYieldKaon[i] = GetGraph(fKaon, Form("gMeanYieldRun3%s", suffix.c_str()));
        gMYieldProton[i] = GetGraph(fProton, Form("gMeanYieldRun3%s", suffix.c_str()));

        // gMYieldPion[i] = GetGraph(fPiKp, Form("gPion_MeanYield%s", suffix.c_str()));
        // gMYieldKaon[i] = GetGraph(fPiKp, Form("gKaon_MeanYield%s", suffix.c_str()));
        // gMYieldProton[i] = GetGraph(fPiKp, Form("gProton_MeanYield%s", suffix.c_str()));

        ////Results from run2 (ALICE)
        gMYieldKstarRun2[i] = GetGraph(fKstarRun2, Form("gKstar_MeanYield%s", suffix.c_str()));
        gMYieldPhiRun2[i] = GetGraph(fPhiRun2, Form("gPhi_MeanYield%s", suffix.c_str()));
        gMYieldChKstarRun2[i] = GetGraph(fChKstar, Form("gChKstar_MeanYield%s", suffix.c_str()));
        gMYieldKshortRun2[i] = GetGraph(fMultStrange, Form("gKshort_MeanYield%s", suffix.c_str()));
        gMYieldLambdaRun2[i] = GetGraph(fMultStrange, Form("gLambda_MeanYield%s", suffix.c_str()));
        gMYieldXiRun2[i] = GetGraph(fMultStrange, Form("gXi_MeanYield%s", suffix.c_str()));
        gMYieldOmegaRun2[i] = GetGraph(fMultStrange, Form("gOmega_MeanYield%s", suffix.c_str()));
        gMYieldXiStarRun2[i] = GetGraph(fXiStarSigma, Form("gXi1530_MeanYield%s", suffix.c_str()));
        gMYieldSigmaRun2[i] = GetGraph(fXiStarSigma, Form("gSigma1385_MeanYield%s", suffix.c_str()));
        gMYieldLambda1520[i] = GetGraph(fLambda1520Run2, Form("gLambda_MeanYield%s", suffix.c_str()));

        // Present for run3 (In run2 the Lambda paper is in IRC review)
        gMPtRho[i] = GetGraph(fRho, Form("gMeanpTRun3%s", suffix.c_str()));
        gMYieldRho[i] = GetGraph(fRho, Form("gMeanYieldRun3%s", suffix.c_str()));

        gMPtKshortRun3[i] = GetGraph(fKshortRun3, Form("gMeanpTRun3%s", suffix.c_str()));
        gMYieldKshortRun3[i] = GetGraph(fKshortRun3, Form("gMeanYieldRun3%s", suffix.c_str()));
    }

    TGraphErrors *gMeanYieldMB_stat = GetGraph(fKstar, "gMeanYieldMB_stat");
    TGraphAsymmErrors *gMeanYieldMB_sys = (TGraphAsymmErrors *)fKstar->Get("gMeanYieldMB_sys");
    if (!gMeanYieldMB_stat || !gMeanYieldMB_sys || gMeanYieldMB_stat->GetN() < 1 || gMeanYieldMB_sys->GetN() < 1)
    {
        cout << "Error: minimum-bias yield graphs not found in " << fKstar->GetName() << endl;
        return;
    }

    TGraphErrors *gKstarKshortRatio_13TeV[2], *gPhiKshortRatio_13TeV[2], *gPhiPiRatio_13TeV[2], *gKstarKaRatio_13TeV[2], *gKstarPiRatio_13TeV[2], *gPhiKaRatio_13TeV[2], *gPhiPrRatio_13TeV[2], *gPhiLambdaRatio_13TeV[2];
    TGraphErrors *gKstarPiRatio_7TeV[2], *gPhiPiRatio_7TeV[2];

    for (int i = 0; i < 2; i++)
    {
        gKstarKshortRatio_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_KshortRatio%s", (i == 0) ? "_stat" : "_sys"));
        gPhiKshortRatio_13TeV[i] = GetGraph(fpp13TeV, Form("gPhi_KshortRatio%s", (i == 0) ? "_stat" : "_sys"));
        gPhiPiRatio_13TeV[i] = GetGraph(fpp13TeV, Form("gPhi_PiRatio%s", (i == 0) ? "_stat" : "_sys"));
        gKstarKaRatio_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_KaRatio%s", (i == 0) ? "_stat" : "_sys"));
        gKstarPiRatio_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_PiRatio%s", (i == 0) ? "_stat" : "_sys"));
        gPhiKaRatio_13TeV[i] = GetGraph(fpp13TeV, Form("gPhi_KaRatio%s", (i == 0) ? "_stat" : "_sys"));
        gPhiPrRatio_13TeV[i] = GetGraph(fpp13TeV, Form("gPhi_PrRatio%s", (i == 0) ? "_stat" : "_sys"));
        gPhiLambdaRatio_13TeV[i] = GetGraph(fpp13TeV, Form("gPhi_LambdaRatio%s", (i == 0) ? "_stat" : "_sys"));

        gKstarPiRatio_7TeV[i] = GetGraph(fpp7TeV, Form("gKstar_PiRatio%s", (i == 0) ? "_stat" : "_sys"));
        gPhiPiRatio_7TeV[i] = GetGraph(fpp7TeV, Form("gPhi_PiRatio%s", (i == 0) ? "_stat" : "_sys"));
    }

    //================================================
    //==============Models from hyperloop outputs======
    //=================================================

    enum Model
    {
        kPythiaMonash,
        kPythiaRopes,
        kPythiaShoving,
        kNModels
    };
    enum Particle
    {
        kKstar,
        kPhi,
        kPion,
        kKaon,
        kProton,
        kPionMinus,
        kKaonMinus,
        kAntiProton,
        kXi1530,
        kKshort,
        kKstarPM,
        kNParticles
    };

    struct ModelStyle
    {
        Color_t color;
        int style;
    };

    ModelStyle modelStyle[kNModels] = {
        {kRed + 1, 4},  // Pythia Monash
        {kCyan + 1, 7}, // Pythia Ropes
        {kMagenta, 3}   // Pythia Shoving
    };

    const char *modelLabel[kNModels] = {
        "Pythia Monash",
        "Pythia Ropes",
        "Pythia Shoving"};

    TFile *ModelsHyperloop = new TFile("ModelRootFiles/ModelResults2.root", "read");
    if (ModelsHyperloop->IsZombie())
    {
        cout << "Error: Hyperloop model file not found" << endl;
        return;
    }

    string hyperloopModels[kNModels] = {"Pythia_Monash", "Pythia_Ropes", "Pythia_Shoving"};
    string particles[kNParticles] = {"Kstar", "Phi", "Pion", "Kaon", "Proton", "PionMinus", "KaonMinus", "AntiProton", "Xi1530", "Kshort", "KstarPM"};
    int colorsPythia[] = {kCyan + 1, kYellow + 1, kMagenta, kRed + 1, kGreen + 2, kBlue + 1};
    int lineStylesPythia[] = {1, 1, 1, 2, 2, 2};

    vector<int> modelsToPlot = {
        kPythiaMonash,
        kPythiaShoving,
        kPythiaRopes};

    TGraphErrors *gMYield[kNModels][kNParticles];
    TGraphErrors *gMeanPt[kNModels][kNParticles];
    TH1D *hSpectra[kNModels][kNParticles];

    for (int iModel = 0; iModel < kNModels; iModel++)
    {
        for (int iParticle = 0; iParticle < kNParticles; iParticle++)
        {
            gMYield[iModel][iParticle] = GetGraph(ModelsHyperloop, hyperloopModels[iModel] + "/" + particles[iParticle] + "/gMeanYield_" + particles[iParticle]);
            gMeanPt[iModel][iParticle] = GetGraph(ModelsHyperloop, hyperloopModels[iModel] + "/" + particles[iParticle] + "/gMeanpT_" + particles[iParticle]);
            hSpectra[iModel][iParticle] = (TH1D *)ModelsHyperloop->Get((hyperloopModels[iModel] + "/" + particles[iParticle] + "/hPt_" + particles[iParticle] + "_MinBias").c_str());

            if (gMYield[iModel][iParticle] == nullptr || gMeanPt[iModel][iParticle] == nullptr || hSpectra[iModel][iParticle] == nullptr)
            {
                cout << "Error: Graphs or histogram for model " << hyperloopModels[iModel] << " and particle " << particles[iParticle] << " not found" << endl;
                return;
            }

            SetGraphErrorStyle(gMYield[iModel][iParticle]);
            SetGraphErrorStyle(gMeanPt[iModel][iParticle]);
            SetHistoQA(hSpectra[iModel][iParticle]);

            // Apply model style once
            gMYield[iModel][iParticle]->SetLineColor(modelStyle[iModel].color);
            gMYield[iModel][iParticle]->SetLineStyle(modelStyle[iModel].style);

            gMeanPt[iModel][iParticle]->SetLineColor(modelStyle[iModel].color);
            gMeanPt[iModel][iParticle]->SetLineStyle(modelStyle[iModel].style);

            int nPoints = gMYield[iModel][iParticle]->GetN();

            // RestrictModelXaxis(gMYield[iModel][iParticle], 3.3, 22.5);
            // RestrictModelXaxis(gMeanPt[iModel][iParticle], 3.3, 22.5);
        }

        // Add pion with pionMinus, kaon with kaonMinus, proton with anti-proton
        for (int i = 0; i < gMYield[iModel][kPion]->GetN(); i++)
        {
            double x, yPion, yPionMinus;
            gMYield[iModel][kPion]->GetPoint(i, x, yPion);
            gMYield[iModel][kPionMinus]->GetPoint(i, x, yPionMinus);
            double yPionTotal = yPion + yPionMinus;
            gMYield[iModel][kPion]->SetPoint(i, x, yPionTotal);

            // Kaon
            double yKaon, yKaonMinus;
            gMYield[iModel][kKaon]->GetPoint(i, x, yKaon);
            gMYield[iModel][kKaonMinus]->GetPoint(i, x, yKaonMinus);
            double yKaonTotal = yKaon + yKaonMinus;
            gMYield[iModel][kKaon]->SetPoint(i, x, yKaonTotal);

            // Proton
            double yProton, yAntiProton;
            gMYield[iModel][kProton]->GetPoint(i, x, yProton);
            gMYield[iModel][kAntiProton]->GetPoint(i, x, yAntiProton);
            double yProtonTotal = yProton + yAntiProton;
            gMYield[iModel][kProton]->SetPoint(i, x, yProtonTotal);
        }
    }

    //================================================
    //==============EPOS Model========================
    //================================================

    TFile *fEPOS_New = OpenFile("/home/sawan/Storage/EPOS4_Resonance/FinalEPOSResults_BongFix3.root");
    TFile *fEPOS_New2 = OpenFile("/home/sawan/Storage/EPOS4_Resonance/FinalEPOSResults_BongFix3.root");

    TGraphErrors *gEPOS_KstardNdy_UrQMDOff = GetGraph(fEPOS_New, "g_dNdy_opt0");
    TGraphErrors *gEPOS_KstardNdy_UrQMDOn = GetGraph(fEPOS_New, "g_dNdy_opt1");
    TGraphErrors *gEPOS_KstarMeanpT_UrQMDOff = GetGraph(fEPOS_New, "g_pT_opt0");
    TGraphErrors *gEPOS_KstarMeanpT_UrQMDOn = GetGraph(fEPOS_New, "g_pT_opt1");

    TGraphErrors *gEPOS_KstardNdy_UrQMDOffMB = GetGraph(fEPOS_New2, "g_dNdy_MB_opt0");
    TGraphErrors *gEPOS_KstardNdy_UrQMDOnMB = GetGraph(fEPOS_New2, "g_dNdy_MB_opt1");
    TGraphErrors *gEPOS_KstarMeanpT_UrQMDOffMB = GetGraph(fEPOS_New2, "g_pT_MB_opt0");
    TGraphErrors *gEPOS_KstarMeanpT_UrQMDOnMB = GetGraph(fEPOS_New2, "g_pT_MB_opt1");

    TGraphErrors *gEPOS_PhidNdy_UrQMDOff = GetGraph(fEPOS_New, "g_phi_dNdy_opt0");
    TGraphErrors *gEPOS_PhidNdy_UrQMDOn = GetGraph(fEPOS_New, "g_phi_dNdy_opt1");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOff = GetGraph(fEPOS_New, "g_phi_pT_opt0");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOn = GetGraph(fEPOS_New, "g_phi_pT_opt1");

    TGraphErrors *gEPOS_PhidNdy_UrQMDOffMB = GetGraph(fEPOS_New2, "g_phi_dNdy_MB_opt0");
    TGraphErrors *gEPOS_PhidNdy_UrQMDOnMB = GetGraph(fEPOS_New2, "g_phi_dNdy_MB_opt1");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOffMB = GetGraph(fEPOS_New2, "g_phi_pT_MB_opt0");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOnMB = GetGraph(fEPOS_New2, "g_phi_pT_MB_opt1");

    for (int i = gEPOS_PhidNdy_UrQMDOff->GetN() - 1; i >= 0; --i)
    {
        gEPOS_PhidNdy_UrQMDOff->SetPointError(i, 0, 0);   // Remove error
        gEPOS_PhidNdy_UrQMDOn->SetPointError(i, 0, 0);    // Remove error
        gEPOS_PhiMeanpT_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_PhiMeanpT_UrQMDOn->SetPointError(i, 0, 0);  // Remove error

        gEPOS_PhidNdy_UrQMDOffMB->SetPointError(i, 0, 0);   // Remove error
        gEPOS_PhidNdy_UrQMDOnMB->SetPointError(i, 0, 0);    // Remove error
        gEPOS_PhiMeanpT_UrQMDOffMB->SetPointError(i, 0, 0); // Remove error
        gEPOS_PhiMeanpT_UrQMDOnMB->SetPointError(i, 0, 0);  // Remove error
    }

    int counter = 0;
    for (auto model : modelsToPlot)
    {

        for (auto particle : {kKstar, kPhi, kPion, kKaon, kProton, kKshort})
        {
            gMYield[model][particle]->SetLineColor(colorsPythia[counter]);
            gMYield[model][particle]->SetLineStyle(lineStylesPythia[counter]);
            gMYield[model][particle]->SetLineWidth(3);

            gMeanPt[model][particle]->SetLineColor(colorsPythia[counter]);
            gMeanPt[model][particle]->SetLineStyle(lineStylesPythia[counter]);
            gMeanPt[model][particle]->SetLineWidth(3);
        }
        counter++;
    }

    ////======================================================
    //    ================Plots=====================
    //======================================================

    TCanvas *cdNdyKstar = new TCanvas("cdNdyKstar", "cdNdyKstar", 720, 720);
    SetCanvasStyle(cdNdyKstar, 0.15, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cdNdyKstar, pad1Size, pad2Size);
    cdNdyKstar->cd(1);

    gMYieldKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMYieldKstar[0]->GetYaxis()->SetTitle("dN/dy");
    gMYieldKstar[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    gMYieldKstar[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    gMYieldKstar[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    gMYieldKstar[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);

    gMYieldKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldKstar[0]->SetLineWidth(3);
    gMYieldKstar[0]->GetYaxis()->SetRangeUser(0, 0.71);
    gMYieldKstar[0]->SetMarkerColor(kRed + 1);
    gMYieldKstar[0]->SetLineColor(kRed + 1);
    gMYieldKstar[0]->Draw("APE");
    gMYieldKstar[1]->SetFillStyle(0);
    gMYieldKstar[1]->SetLineColor(kRed + 1);
    gMYieldKstar[1]->SetLineWidth(3);
    gMYieldKstar[1]->Draw("5 same");
    gMYieldKstar[2]->SetFillStyle(1001);
    gMYieldKstar[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gMYieldKstar[2]->Draw("E2 same");
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

    int TotalPointsRun3 = gMYieldKstar[0]->GetN();
    int TotalPointsRun2 = gKstar_Yield_13TeV[0]->GetN();

    for (auto model : modelsToPlot)
    {
        ScaleGraph(gMYield[model][kKstar], 0.5);
        gMYield[model][kKstar]->Draw("l same");
    }

    int totalPointsEPOS = gEPOS_KstardNdy_UrQMDOff->GetN();
    for (int i = 0; i < totalPointsEPOS; i++)
    {
        double x, y;
        gEPOS_KstardNdy_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_KstardNdy_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
    }

    gEPOS_KstardNdy_UrQMDOff->SetLineColor(kGreen + 2);
    gEPOS_KstardNdy_UrQMDOff->SetLineStyle(2);
    gEPOS_KstardNdy_UrQMDOff->SetLineWidth(3);
    gEPOS_KstardNdy_UrQMDOff->Draw("l same");
    gEPOS_KstardNdy_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_KstardNdy_UrQMDOn->SetLineStyle(9);
    gEPOS_KstardNdy_UrQMDOn->SetLineWidth(3);
    gEPOS_KstardNdy_UrQMDOn->Draw("l same");

    TLegend *legend = new TLegend(0.4, 0.82, 0.8, 0.91);
    SetLegendStyle(legend);
    legend->SetTextSize(0.035);
    legend->AddEntry(gMYieldKstar[0], "pp, #sqrt{s} = 13.6 TeV", "P");
    legend->AddEntry(gKstar_Yield_13TeV[0], "pp, #sqrt{s} = 13 TeV", "P");
    legend->Draw();

    // TLegend *legend2 = new TLegend(0.2, 0.66, 0.5, 0.92);
    TLegend *legend2 = new TLegend(0.45, 0.66, 0.9, 0.92);
    SetLegendStyle(legend2);
    legend2->SetTextSize(0.035);
    legend2->SetHeader("Models: pp #sqrt{s} = 13.6 TeV");

    legend2->AddEntry(gEPOS_KstardNdy_UrQMDOff, "EPOS4", "L");
    legend2->AddEntry(gEPOS_KstardNdy_UrQMDOn, "EPOS4 + UrQMD", "L");

    for (auto model : modelsToPlot)
    {
        legend2->AddEntry(gMYield[model][kKstar], modelLabel[model], "L");
    }
    legend->Draw();
    // legend2->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.22, 0.9, "ALICE");
    latex.DrawLatex(0.22, 0.78, "|y| < 0.5");
    latex.SetTextFont(22);
    latex.DrawLatex(0.22, 0.84, "K* (892)^{0}");
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.52, 0.2, "Uncertainties:");
    latex.SetTextFont(42);
    latex.DrawLatex(0.52, 0.15, "stat. (bars), total sys. (open box)");
    latex.DrawLatex(0.52, 0.1, "uncorr. sys. (shaded box)");

    cdNdyKstar->cd(2);
    cdNdyKstar->cd(2)->SetGridy();

    TGraphErrors *gRatioDataModelKstardNdy[kNModels];
    counter = 0;
    for (auto model : modelsToPlot)
    {
        gRatioDataModelKstardNdy[model] = MakeDataModelRatio(gMYieldKstar[0], gMYield[model][kKstar]);
        gRatioDataModelKstardNdy[model]->SetLineColor(colorsPythia[counter]);
        gRatioDataModelKstardNdy[model]->SetMarkerColor(colorsPythia[counter]);
        gRatioDataModelKstardNdy[model]->SetMarkerStyle(20);
        gRatioDataModelKstardNdy[model]->SetLineWidth(2);
        gRatioDataModelKstardNdy[model]->GetYaxis()->SetTitle("Data / Model");
        gRatioDataModelKstardNdy[model]->GetYaxis()->SetTitleOffset(0.5);
        gRatioDataModelKstardNdy[model]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
        gRatioDataModelKstardNdy[model]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
        gRatioDataModelKstardNdy[model]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        gRatioDataModelKstardNdy[model]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        gRatioDataModelKstardNdy[model]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        gRatioDataModelKstardNdy[model]->GetXaxis()->SetLimits(0, 27);
        gRatioDataModelKstardNdy[model]->GetYaxis()->SetRangeUser(0.5, 1.3);
        gRatioDataModelKstardNdy[model]->GetYaxis()->SetNdivisions(505);
        if (counter == 0)
        {
            gRatioDataModelKstardNdy[model]->Draw("APE");
        }
        else
        {
            gRatioDataModelKstardNdy[model]->Draw("PE same");
        }
        counter++;
    }
    TGraphErrors *gRatioDataModelKstardNdyEPOSON = MakeDataModelRatio(gMYieldKstar[0], gEPOS_KstardNdy_UrQMDOn);
    TGraphErrors *gRatioDataModelKstardNdyEPOSOFF = MakeDataModelRatio(gMYieldKstar[0], gEPOS_KstardNdy_UrQMDOff);

    gRatioDataModelKstardNdyEPOSON->SetLineColor(kBlue - 6);
    gRatioDataModelKstardNdyEPOSON->SetMarkerColor(kBlue - 6);
    gRatioDataModelKstardNdyEPOSON->SetMarkerStyle(20);
    gRatioDataModelKstardNdyEPOSON->SetLineWidth(2);
    gRatioDataModelKstardNdyEPOSON->Draw("PE same");
    gRatioDataModelKstardNdyEPOSOFF->SetLineColor(kGreen + 2);
    gRatioDataModelKstardNdyEPOSOFF->SetMarkerColor(kGreen + 2);
    gRatioDataModelKstardNdyEPOSOFF->SetMarkerStyle(20);
    gRatioDataModelKstardNdyEPOSOFF->SetLineWidth(2);
    gRatioDataModelKstardNdyEPOSOFF->Draw("PE same");

    TLine *line = new TLine(0, 1, 27, 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    if (isSavePlots)
    {
        cdNdyKstar->SaveAs("Plots/QA/MeanYield_Kstar.png");
    }

    //===================================================
    //  ================dN/dy Phi======================
    //===================================================

    TCanvas *cdNdyPhi = new TCanvas("cdNdyPhi", "cdNdyPhi", 720, 720);
    SetCanvasStyle(cdNdyPhi, 0.15, 0.03, 0.03, 0.15);
    canvas_style(cdNdyPhi, pad1Size, pad2Size);
    cdNdyPhi->cd(1);

    gMYieldPhiRun3[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMYieldPhiRun3[0]->GetYaxis()->SetTitle("dN/dy");
    gMYieldPhiRun3[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldPhiRun3[0]->SetLineWidth(3);
    gMYieldPhiRun3[0]->GetYaxis()->SetRangeUser(0, 0.26);
    gMYieldPhiRun3[0]->SetMarkerColor(kRed + 1);
    gMYieldPhiRun3[0]->SetLineColor(kRed + 1);
    gMYieldPhiRun3[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    gMYieldPhiRun3[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    gMYieldPhiRun3[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    gMYieldPhiRun3[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    gMYieldPhiRun3[0]->Draw("APE");
    gMYieldPhiRun3[1]->SetFillStyle(0);
    gMYieldPhiRun3[1]->SetLineColor(kRed + 1);
    gMYieldPhiRun3[1]->SetLineWidth(3);
    gMYieldPhiRun3[1]->Draw("5 same");
    gMYieldPhiRun3[2]->SetFillStyle(1001);
    gMYieldPhiRun3[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gMYieldPhiRun3[2]->Draw("E2 same");
    gMYieldPhiRun2[0]->SetMarkerStyle(21);
    gMYieldPhiRun2[0]->SetMarkerColor(kBlue + 1);
    gMYieldPhiRun2[0]->SetLineColor(kBlue + 1);
    gMYieldPhiRun2[0]->SetLineWidth(3);
    gMYieldPhiRun2[0]->Draw("P same");
    gMYieldPhiRun2[1]->SetLineColor(kBlue + 1);
    gMYieldPhiRun2[1]->SetFillStyle(0);
    gMYieldPhiRun2[1]->SetLineWidth(3);
    gMYieldPhiRun2[1]->Draw("5 same");
    gMYieldPhiRun2[2]->SetFillStyle(1001);
    gMYieldPhiRun2[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gMYieldPhiRun2[2]->Draw("E2 same");

    int totalPointsEPOS_Phi = gEPOS_PhidNdy_UrQMDOff->GetN();
    for (int i = 0; i < totalPointsEPOS_Phi; i++)
    {
        gEPOS_PhidNdy_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_PhidNdy_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
    }

    gEPOS_PhidNdy_UrQMDOff->SetLineColor(kGreen + 2);
    gEPOS_PhidNdy_UrQMDOff->SetLineStyle(2);
    gEPOS_PhidNdy_UrQMDOff->SetLineWidth(3);
    gEPOS_PhidNdy_UrQMDOff->Draw("l same");
    gEPOS_PhidNdy_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_PhidNdy_UrQMDOn->SetLineStyle(9);
    gEPOS_PhidNdy_UrQMDOn->SetLineWidth(3);
    gEPOS_PhidNdy_UrQMDOn->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMYield[model][kPhi]->Draw("l same");
    }

    legend2->Draw();
    latex.SetTextFont(22);
    latex.DrawLatex(0.76, 0.88, "#Phi (1020)");

    cdNdyPhi->cd(2);
    cdNdyPhi->cd(2)->SetGridy();

    TGraphErrors *gRatioDataModelPhidNdy[kNModels];
    counter = 0;
    for (auto model : modelsToPlot)
    {
        gRatioDataModelPhidNdy[model] = MakeDataModelRatio(gMYieldPhiRun3[0], gMYield[model][kPhi]);
        gRatioDataModelPhidNdy[model]->SetLineColor(colorsPythia[counter]);
        gRatioDataModelPhidNdy[model]->SetMarkerColor(colorsPythia[counter]);
        gRatioDataModelPhidNdy[model]->SetMarkerStyle(20);
        gRatioDataModelPhidNdy[model]->SetLineWidth(2);
        gRatioDataModelPhidNdy[model]->GetYaxis()->SetTitle("Data / Model");
        gRatioDataModelPhidNdy[model]->GetYaxis()->SetTitleOffset(0.5);
        gRatioDataModelPhidNdy[model]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
        gRatioDataModelPhidNdy[model]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
        gRatioDataModelPhidNdy[model]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        gRatioDataModelPhidNdy[model]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        gRatioDataModelPhidNdy[model]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        gRatioDataModelPhidNdy[model]->GetXaxis()->SetLimits(0, 27);
        gRatioDataModelPhidNdy[model]->GetYaxis()->SetRangeUser(0.6, 1.7);
        gRatioDataModelPhidNdy[model]->GetYaxis()->SetNdivisions(508);
        if (counter == 0)
        {
            gRatioDataModelPhidNdy[model]->Draw("APE");
        }
        else
        {
            gRatioDataModelPhidNdy[model]->Draw("PE same");
        }
        counter++;
    }

    TGraphErrors *gRatioDataModelPhidNdyEPOSON = MakeDataModelRatio(gMYieldPhiRun3[0], gEPOS_PhidNdy_UrQMDOn);
    TGraphErrors *gRatioDataModelPhidNdyEPOSOFF = MakeDataModelRatio(gMYieldPhiRun3[0], gEPOS_PhidNdy_UrQMDOff);

    gRatioDataModelPhidNdyEPOSON->SetLineColor(kBlue - 6);
    gRatioDataModelPhidNdyEPOSON->SetMarkerColor(kBlue - 6);
    gRatioDataModelPhidNdyEPOSON->SetMarkerStyle(20);
    gRatioDataModelPhidNdyEPOSON->SetLineWidth(2);
    gRatioDataModelPhidNdyEPOSON->Draw("PE same");
    gRatioDataModelPhidNdyEPOSOFF->SetLineColor(kGreen + 2);
    gRatioDataModelPhidNdyEPOSOFF->SetMarkerColor(kGreen + 2);
    gRatioDataModelPhidNdyEPOSOFF->SetMarkerStyle(20);
    gRatioDataModelPhidNdyEPOSOFF->SetLineWidth(2);
    gRatioDataModelPhidNdyEPOSOFF->Draw("PE same");

    TLine *linePhi = new TLine(0, 1, 27, 1);
    linePhi->SetLineStyle(2);
    linePhi->SetLineWidth(2);
    linePhi->Draw("same");

    if (isSavePlots)
    {
        cdNdyPhi->SaveAs("Plots/QA/MeanYield_Phi.png");
    }

    //===================================================
    //  ================<pT> K*0======================
    //===================================================
    TCanvas *cMeanPtKstar = new TCanvas("cMeanPtKstar", "cMeanPtKstar", 720, 720);
    SetCanvasStyle(cMeanPtKstar, 0.15, 0.03, 0.03, 0.15);
    canvas_style(cMeanPtKstar, pad1Size, pad2Size);
    cMeanPtKstar->cd(1);

    SetGraphErrorStyle(gMPtKstar[0]);
    gMPtKstar[0]->SetTitle(0);
    gMPtKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMPtKstar[0]->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/#it{c})");
    gMPtKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMPtKstar[0]->GetYaxis()->SetRangeUser(0.58, 1.76);
    gMPtKstar[0]->SetMarkerColor(kRed);
    gMPtKstar[0]->SetLineColor(kRed);
    gMPtKstar[0]->SetLineWidth(3);
    gMPtKstar[0]->GetYaxis()->SetTitleOffset(1.3);
    gMPtKstar[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    gMPtKstar[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    gMPtKstar[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    gMPtKstar[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    gMPtKstar[0]->Draw("APE");
    gMPtKstar[1]->SetFillStyle(0);
    gMPtKstar[1]->SetLineColor(kRed);
    gMPtKstar[1]->SetLineWidth(3);
    gMPtKstar[1]->Draw("5 same");
    gMPtKstar[2]->SetFillStyle(1001);
    gMPtKstar[2]->SetFillColorAlpha(kRed, 0.3);
    gMPtKstar[2]->Draw("E2 same");

    gKstar_MeanpT_13TeV[0]->SetMarkerStyle(21);
    gKstar_MeanpT_13TeV[0]->SetMarkerColor(kBlue);
    gKstar_MeanpT_13TeV[0]->SetLineColor(kBlue);
    gKstar_MeanpT_13TeV[0]->SetLineWidth(3);
    gKstar_MeanpT_13TeV[0]->Draw("P same");
    gKstar_MeanpT_13TeV[1]->SetLineColor(kBlue);
    gKstar_MeanpT_13TeV[1]->SetFillStyle(0);
    gKstar_MeanpT_13TeV[1]->SetLineWidth(3);
    gKstar_MeanpT_13TeV[1]->Draw("5 same");
    gKstar_MeanpT_13TeV[2]->SetFillStyle(1001);
    gKstar_MeanpT_13TeV[2]->SetFillColorAlpha(kBlue, 0.3);
    gKstar_MeanpT_13TeV[2]->Draw("E2 same");

    int totalPointsEPOS_MeanPt = gEPOS_KstarMeanpT_UrQMDOff->GetN();
    for (int i = 0; i < totalPointsEPOS_MeanPt; i++)
    {
        gEPOS_KstarMeanpT_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_KstarMeanpT_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
    }

    gEPOS_KstarMeanpT_UrQMDOff->SetLineColor(kGreen + 2);
    gEPOS_KstarMeanpT_UrQMDOff->SetLineStyle(2);
    gEPOS_KstarMeanpT_UrQMDOff->SetLineWidth(3);
    gEPOS_KstarMeanpT_UrQMDOff->Draw("l same");
    gEPOS_KstarMeanpT_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_KstarMeanpT_UrQMDOn->SetLineStyle(9);
    gEPOS_KstarMeanpT_UrQMDOn->SetLineWidth(3);
    gEPOS_KstarMeanpT_UrQMDOn->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMeanPt[model][kKstar]->Draw("l same");
    }

    legend->Draw();
    // legend2->Draw();
    latex.DrawLatex(0.22, 0.9, "ALICE");
    latex.DrawLatex(0.22, 0.78, "|y| < 0.5");
    latex.SetTextFont(22);
    latex.DrawLatex(0.22, 0.84, "K* (892)^{0}");
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.52, 0.2, "Uncertainties:");
    latex.SetTextFont(42);
    latex.DrawLatex(0.52, 0.15, "stat. (bars), total sys. (open box)");
    latex.DrawLatex(0.52, 0.1, "uncorr. sys. (shaded box)");

    cMeanPtKstar->cd(2);
    cMeanPtKstar->cd(2)->SetGridy();

    TGraphErrors *gRatioDataModelKstarMeanPt[kNModels];
    counter = 0;
    for (auto model : modelsToPlot)
    {
        gRatioDataModelKstarMeanPt[model] = MakeDataModelRatio(gMPtKstar[0], gMeanPt[model][kKstar]);
        gRatioDataModelKstarMeanPt[model]->SetLineColor(colorsPythia[counter]);
        gRatioDataModelKstarMeanPt[model]->SetMarkerColor(colorsPythia[counter]);
        gRatioDataModelKstarMeanPt[model]->SetMarkerStyle(20);
        gRatioDataModelKstarMeanPt[model]->SetLineWidth(2);
        gRatioDataModelKstarMeanPt[model]->GetYaxis()->SetTitle("Data / Model");
        gRatioDataModelKstarMeanPt[model]->GetYaxis()->SetTitleOffset(0.5);
        gRatioDataModelKstarMeanPt[model]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
        gRatioDataModelKstarMeanPt[model]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
        gRatioDataModelKstarMeanPt[model]->GetYaxis()->SetLabelSize(0.035 / pad2Size);
        gRatioDataModelKstarMeanPt[model]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        gRatioDataModelKstarMeanPt[model]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        gRatioDataModelKstarMeanPt[model]->GetXaxis()->SetLimits(0, 27);
        gRatioDataModelKstarMeanPt[model]->GetYaxis()->SetRangeUser(0.9, 1.5);
        gRatioDataModelKstarMeanPt[model]->GetYaxis()->SetNdivisions(508);
        if (counter == 0)
        {
            gRatioDataModelKstarMeanPt[model]->Draw("APE");
        }
        else
        {
            gRatioDataModelKstarMeanPt[model]->Draw("PE same");
        }
        counter++;
    }

    TGraphErrors *gRatioDataModelKstarMeanPtEPOSON = MakeDataModelRatio(gMPtKstar[0], gEPOS_KstarMeanpT_UrQMDOn);
    TGraphErrors *gRatioDataModelKstarMeanPtEPOSOFF = MakeDataModelRatio(gMPtKstar[0], gEPOS_KstarMeanpT_UrQMDOff);

    gRatioDataModelKstarMeanPtEPOSON->SetLineColor(kBlue - 6);
    gRatioDataModelKstarMeanPtEPOSON->SetMarkerColor(kBlue - 6);
    gRatioDataModelKstarMeanPtEPOSON->SetMarkerStyle(20);
    gRatioDataModelKstarMeanPtEPOSON->SetLineWidth(2);
    gRatioDataModelKstarMeanPtEPOSON->Draw("PE same");
    gRatioDataModelKstarMeanPtEPOSOFF->SetLineColor(kGreen + 2);
    gRatioDataModelKstarMeanPtEPOSOFF->SetMarkerColor(kGreen + 2);
    gRatioDataModelKstarMeanPtEPOSOFF->SetMarkerStyle(20);
    gRatioDataModelKstarMeanPtEPOSOFF->SetLineWidth(2);
    gRatioDataModelKstarMeanPtEPOSOFF->Draw("PE same");

    TLine *lineMeanPtKstar = new TLine(0, 1, 27, 1);
    lineMeanPtKstar->SetLineStyle(2);
    lineMeanPtKstar->SetLineWidth(2);
    lineMeanPtKstar->Draw("same");

    if (isSavePlots)
    {
        cMeanPtKstar->SaveAs("Plots/QA/MeanPt_Kstar.png");
    }

    //===================================================
    // ============== Normalized dN/dy Kstar ============
    //===================================================

    TGraphErrors *gNormYieldKstar[3];
    TGraphErrors *gNormYieldKstarRun2[3];

    for (int i = 0; i < 3; i++)
    {
        gNormYieldKstar[i] = NormalizeGraphByX(gMYieldKstar[i], Form("gNormYieldKstar_%d", i));
        gNormYieldKstarRun2[i] = NormalizeGraphByX(gKstar_Yield_13TeV[i], Form("gNormYieldKstarRun2_%d", i));
    }

    TCanvas *cNormYieldKstar = new TCanvas("cNormYieldKstar", "cNormYieldKstar", 720, 720);
    SetCanvasStyle(cNormYieldKstar, 0.17, 0.03, 0.03, 0.15);
    SetGraphStyle(gNormYieldKstar[0]);
    gNormYieldKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gNormYieldKstar[0]->GetYaxis()->SetTitle("(dN/dy) / #LT dN_{ch}/d#eta #GT");
    gNormYieldKstar[0]->GetXaxis()->SetLimits(0, 27);
    gNormYieldKstar[0]->SetMarkerColor(kRed + 1);
    gNormYieldKstar[0]->SetLineColor(kRed + 1);
    gNormYieldKstar[0]->SetLineWidth(3);
    gNormYieldKstar[0]->GetYaxis()->SetRangeUser(0.0155, 0.031);
    gNormYieldKstar[0]->GetYaxis()->SetTitleOffset(1.9);
    gNormYieldKstar[0]->Draw("APE");

    gNormYieldKstar[1]->SetFillStyle(0);
    gNormYieldKstar[1]->SetLineColor(kRed + 1);
    gNormYieldKstar[1]->SetLineWidth(3);
    gNormYieldKstar[1]->Draw("5 same");

    gNormYieldKstar[2]->SetFillStyle(1001);
    gNormYieldKstar[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gNormYieldKstar[2]->Draw("E2 same");

    // Run 2 data
    gNormYieldKstarRun2[0]->SetMarkerStyle(21);
    gNormYieldKstarRun2[0]->SetMarkerColor(kBlue + 1);
    gNormYieldKstarRun2[0]->SetLineColor(kBlue + 1);
    gNormYieldKstarRun2[0]->Draw("P same");

    gNormYieldKstarRun2[1]->SetFillStyle(0);
    gNormYieldKstarRun2[1]->SetLineColor(kBlue + 1);
    gNormYieldKstarRun2[1]->Draw("5 same");

    gNormYieldKstarRun2[2]->SetFillStyle(1001);
    gNormYieldKstarRun2[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gNormYieldKstarRun2[2]->Draw("E2 same");

    // Models
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gr = NormalizeGraphByX(gMYield[model][kKstar], Form("gNormModelYield_%d", model));
        gr->Draw("l same");
    }

    // Normalized EPOS graphs
    TGraphErrors *gNormEPOS_KstardNdy_UrQMDOff = NormalizeGraphByX(gEPOS_KstardNdy_UrQMDOff, "gNormEPOS_KstardNdy_UrQMDOff");
    TGraphErrors *gNormEPOS_KstardNdy_UrQMDOn = NormalizeGraphByX(gEPOS_KstardNdy_UrQMDOn, "gNormEPOS_KstardNdy_UrQMDOn");

    TGraphErrors *gNormEPOS_KstarMeanpT_UrQMDOff = NormalizeGraphByX(gEPOS_KstarMeanpT_UrQMDOff, "gNormEPOS_KstarMeanpT_UrQMDOff");
    TGraphErrors *gNormEPOS_KstarMeanpT_UrQMDOn = NormalizeGraphByX(gEPOS_KstarMeanpT_UrQMDOn, "gNormEPOS_KstarMeanpT_UrQMDOn");

    // EPOS4
    gNormEPOS_KstardNdy_UrQMDOff->SetLineColor(kGreen + 2);
    gNormEPOS_KstardNdy_UrQMDOff->SetLineStyle(2);
    gNormEPOS_KstardNdy_UrQMDOff->SetLineWidth(3);
    gNormEPOS_KstardNdy_UrQMDOff->Draw("l same");

    // EPOS4 + UrQMD
    gNormEPOS_KstardNdy_UrQMDOn->SetLineColor(kBlue - 6);
    gNormEPOS_KstardNdy_UrQMDOn->SetLineStyle(9);
    gNormEPOS_KstardNdy_UrQMDOn->SetLineWidth(3);
    gNormEPOS_KstardNdy_UrQMDOn->Draw("l same");

    TLegend *legNormYield = new TLegend(0.4, 0.82, 0.8, 0.91);
    SetLegendStyle(legNormYield);
    legNormYield->AddEntry(gNormYieldKstar[0], "pp, #sqrt{s} = 13.6 TeV", "P");
    legNormYield->AddEntry(gNormYieldKstarRun2[0], "pp, #sqrt{s} = 13 TeV", "P");
    legNormYield->Draw();

    latex.DrawLatex(0.22, 0.9, "ALICE");
    latex.SetTextFont(22);
    latex.DrawLatex(0.22, 0.84, "K* (892)^{0}");
    latex.SetTextFont(42);
    latex.DrawLatex(0.22, 0.78, "|y| < 0.5");

    if (isSavePlots)
        cNormYieldKstar->SaveAs("Plots/NormalizedYield_Kstar.png");

    //=========================================================
    // ===== dN/dy normalized to MB vs dNch/deta normalized MB
    //=========================================================

    double xMB, yieldMB;
    gMeanYieldMB_stat->GetPoint(0, xMB, yieldMB);
    const double dNchdetaMB = 7.087;
    const double dNchdetaMBErr = 0.18;
    const double dNchdetaMB_run2 = 6.89;
    double yieldMBStatErr = gMeanYieldMB_stat->GetErrorY(0);

    TGraphErrors *gNormYieldKstarToMB_stat = NormalizeYieldByMB(gMYieldKstar[0], yieldMB, yieldMBStatErr, dNchdetaMB, dNchdetaMBErr, "gNormYieldKstarToMB_stat");
    TGraphAsymmErrors *gNormYieldKstarToMB_sys = NormalizeYieldSysByMB(gMYieldKstar[1], gMeanYieldMB_sys, yieldMB, dNchdetaMB, dNchdetaMBErr, "gNormYieldKstarToMB_sys");
    TGraphAsymmErrors *gNormYieldKstarToMB_sysUncorr = new TGraphAsymmErrors();
    gNormYieldKstarToMB_sysUncorr->SetName("gNormYieldKstarToMB_sysuncorr");

    TGraphErrors *gNormYieldKstarToMBRun2_stat = NormalizeYieldByMB(gKstar_Yield_13TeV[0], yieldMB, yieldMBStatErr, dNchdetaMB_run2, dNchdetaMBErr, "gNormYieldKstarToMBRun2_stat");
    TGraphAsymmErrors *gNormYieldKstarToMBRun2_sys = NormalizeYieldSysByMB(gKstar_Yield_13TeV[1], gMeanYieldMB_sys, yieldMB, dNchdetaMB_run2, dNchdetaMBErr, "gNormYieldKstarToMBRun2_sys");
    TGraphAsymmErrors *gNormYieldKstarToMBRun2_sysUncorr = new TGraphAsymmErrors();
    gNormYieldKstarToMBRun2_sysUncorr->SetName("gNormYieldKstarToMBRun2_sysuncorr");

    TGraphErrors *gMBYieldPythiaKstar[3], *gMBYieldPythiaPhi[3];
    TGraphErrors *gNormalizeYieldKstarPythia[3], *gNormalizeYieldPhiPythia[3];

    string tempModelNames[3] = {"Pythia_Monash", "Pythia_Ropes", "Pythia_Shoving"};

    for (int i = 0; i < 3; i++)
    {
        gMBYieldPythiaKstar[i] = GetGraph(ModelsHyperloop, Form("%s/Kstar/gMeanYield_Kstar_MB", tempModelNames[i].c_str()));
        gMBYieldPythiaPhi[i] = GetGraph(ModelsHyperloop, Form("%s/Phi/gMeanYield_Phi_MB", tempModelNames[i].c_str()));

        if (!gMBYieldPythiaKstar[i] || !gMBYieldPythiaPhi[i])
        {
            cout << "Error: Could not retrieve MB yield graph for model " << tempModelNames[i] << endl;
            return;
        }

        // Kstar MB values
        double modelMultMB_Kstar, modelYieldMB_Kstar;
        gMBYieldPythiaKstar[i]->GetPoint(0, modelMultMB_Kstar, modelYieldMB_Kstar);

        double modelMultMBErr_Kstar = gMBYieldPythiaKstar[i]->GetErrorX(0);
        double modelYieldMBErr_Kstar = gMBYieldPythiaKstar[i]->GetErrorY(0);

        // Phi MB values
        double modelMultMB_Phi, modelYieldMB_Phi;
        gMBYieldPythiaPhi[i]->GetPoint(0, modelMultMB_Phi, modelYieldMB_Phi);

        double modelMultMBErr_Phi = gMBYieldPythiaPhi[i]->GetErrorX(0);
        double modelYieldMBErr_Phi = gMBYieldPythiaPhi[i]->GetErrorY(0);

        // Normalize multiplicity-dependent model graphs using their own MB values
        gNormalizeYieldKstarPythia[i] = NormalizeYieldByMB(gMYield[i][kKstar], modelYieldMB_Kstar, modelYieldMBErr_Kstar, modelMultMB_Kstar, modelMultMBErr_Kstar, Form("gNormYieldKstarPythia_%d", i));

        gNormalizeYieldPhiPythia[i] = NormalizeYieldByMB(gMYield[i][kPhi], modelYieldMB_Phi, modelYieldMBErr_Phi, modelMultMB_Phi, modelMultMBErr_Phi, Form("gNormYieldPhiPythia_%d", i));
    }

    double yieldMB_EPOSOff, yieldMBStatErr_EPOSOff, xMBOff, xMBOff_err, xMBOn, xMBOn_err, yieldMB_EPOSOn, yieldMBStatErr_EPOSOn;
    gEPOS_KstardNdy_UrQMDOffMB->GetPoint(0, xMBOff, yieldMB_EPOSOff);
    yieldMBStatErr_EPOSOff = gEPOS_KstardNdy_UrQMDOffMB->GetErrorY(0);
    gEPOS_KstardNdy_UrQMDOnMB->GetPoint(0, xMBOn, yieldMB_EPOSOn);
    yieldMBStatErr_EPOSOn = gEPOS_KstardNdy_UrQMDOnMB->GetErrorY(0);
    xMBOff_err = gEPOS_KstardNdy_UrQMDOffMB->GetErrorX(0);
    xMBOn_err = gEPOS_KstardNdy_UrQMDOnMB->GetErrorX(0);

    TGraphErrors *gNormYieldKstarToMB_EPOSOn = NormalizeYieldByMB(gEPOS_KstardNdy_UrQMDOn, yieldMB_EPOSOn, yieldMBStatErr_EPOSOn, xMBOn, xMBOn_err, "gNormYieldKstarToMB_EPOSOn");

    TGraphErrors *gNormYieldKstarToMB_EPOSOff = NormalizeYieldByMB(gEPOS_KstardNdy_UrQMDOff, yieldMB_EPOSOff, yieldMBStatErr_EPOSOff, xMBOff, xMBOff_err, "gNormYieldKstarToMB_EPOSOff");

    double yieldMB_EPOSOff_Phi, yieldMBStatErr_EPOSOff_Phi, xMBOff_Phi, xMBOff_err_Phi, xMBOn_Phi, xMBOn_err_Phi, yieldMB_EPOSOn_Phi, yieldMBStatErr_EPOSOn_Phi;
    gEPOS_PhidNdy_UrQMDOffMB->GetPoint(0, xMBOff_Phi, yieldMB_EPOSOff_Phi);
    yieldMBStatErr_EPOSOff_Phi = gEPOS_PhidNdy_UrQMDOffMB->GetErrorY(0);
    gEPOS_PhidNdy_UrQMDOnMB->GetPoint(0, xMBOn_Phi, yieldMB_EPOSOn_Phi);
    yieldMBStatErr_EPOSOn_Phi = gEPOS_PhidNdy_UrQMDOnMB->GetErrorY(0);
    xMBOff_err_Phi = gEPOS_PhidNdy_UrQMDOffMB->GetErrorX(0);
    xMBOn_err_Phi = gEPOS_PhidNdy_UrQMDOnMB->GetErrorX(0);

    TGraphErrors *gNormYieldPhiToMB_EPOSOn = NormalizeYieldByMB(gEPOS_PhidNdy_UrQMDOn, yieldMB_EPOSOn_Phi, yieldMBStatErr_EPOSOn_Phi, xMBOn_Phi, xMBOn_err_Phi, "gNormYieldPhiToMB_EPOSOn");
    TGraphErrors *gNormYieldPhiToMB_EPOSOff = NormalizeYieldByMB(gEPOS_PhidNdy_UrQMDOff, yieldMB_EPOSOff_Phi, yieldMBStatErr_EPOSOff_Phi, xMBOff_Phi, xMBOff_err_Phi, "gNormYieldPhiToMB_EPOSOff");

    for (int i = 0; i < gMYieldKstar[2]->GetN(); i++)
    {
        double x, y;
        gMYieldKstar[2]->GetPoint(i, x, y);
        double ex = gMYieldKstar[2]->GetErrorX(i);
        double ey = gMYieldKstar[2]->GetErrorY(i);
        double xNorm = x / dNchdetaMB;
        double yNorm = y / yieldMB;
        // double exNorm = sqrt(pow(ex / dNchdetaMB, 2) + pow(x * dNchdetaMBErr / (dNchdetaMB * dNchdetaMB), 2));
        double exNorm = ex / dNchdetaMB;
        double eyNorm = ey / yieldMB;
        int ip = gNormYieldKstarToMB_sysUncorr->GetN();
        gNormYieldKstarToMB_sysUncorr->SetPoint(ip, xNorm, yNorm);
        gNormYieldKstarToMB_sysUncorr->SetPointError(ip, exNorm, exNorm, eyNorm, eyNorm);
    }

    for (int i = 0; i < gKstar_Yield_13TeV[2]->GetN(); i++)
    {
        double x, y;
        gKstar_Yield_13TeV[2]->GetPoint(i, x, y);
        double ex = gKstar_Yield_13TeV[2]->GetErrorX(i);
        double ey = gKstar_Yield_13TeV[2]->GetErrorY(i);
        double xNorm = x / dNchdetaMB_run2;
        double yNorm = y / yieldMB;
        // double exNorm = sqrt(pow(ex / dNchdetaMB, 2) + pow(x * dNchdetaMBErr / (dNchdetaMB * dNchdetaMB), 2));
        double exNorm = ex / dNchdetaMB_run2;
        double eyNorm = ey / yieldMB;
        int ip = gNormYieldKstarToMBRun2_sysUncorr->GetN();
        gNormYieldKstarToMBRun2_sysUncorr->SetPoint(ip, xNorm, yNorm);
        gNormYieldKstarToMBRun2_sysUncorr->SetPointError(ip, exNorm, exNorm, eyNorm, eyNorm);
    }

    // Fit the Run-2 normalized yield and evaluate it at the Run-3 multiplicity points.
    double fitXmin = 1e9, fitXmax = -1e9;
    for (int i = 0; i < gNormYieldKstarToMBRun2_stat->GetN(); i++)
    {
        double x, y;
        gNormYieldKstarToMBRun2_stat->GetPoint(i, x, y);
        if (y > 0)
        {
            fitXmin = TMath::Min(fitXmin, x);
            fitXmax = TMath::Max(fitXmax, x);
        }
    }

    if (fitXmin >= fitXmax)
    {
        cout << "Error: insufficient positive Run-2 points for the normalized-yield fit." << endl;
        return;
    }

    TF1 *fitNormYieldRun2 = new TF1("fitNormYieldRun2", "pol1", fitXmin, fitXmax);
    TFitResultPtr fitResultRun2 = gNormYieldKstarToMBRun2_stat->Fit(fitNormYieldRun2, "QRS0");
    fitNormYieldRun2->SetLineColor(kBlue + 1);
    fitNormYieldRun2->SetLineStyle(2);
    fitNormYieldRun2->SetLineWidth(2);

    // The MB yield is the common normalization denominator for Run 2 and Run 3.
    // Its relative systematic uncertainty therefore cancels in Run-3 / Run-2.
    double mbRelSysLow = gMeanYieldMB_sys->GetErrorYlow(0) / yieldMB;
    double mbRelSysHigh = gMeanYieldMB_sys->GetErrorYhigh(0) / yieldMB;

    // Interpolate the numerator-only systematic uncertainty after removing the
    // common MB normalization component in quadrature from the normalized graph.
    auto GetInterpolatedSys = [mbRelSysLow, mbRelSysHigh](TGraphAsymmErrors *gr, double x, double &errLow, double &errHigh)
    {
        errLow = 0.0;
        errHigh = 0.0;
        if (!gr || gr->GetN() < 1)
            return;

        int iLow = -1, iHigh = -1;
        double xLow = -1e30, xHigh = 1e30;
        for (int ip = 0; ip < gr->GetN(); ip++)
        {
            double gx, gy;
            gr->GetPoint(ip, gx, gy);
            if (gx <= x && gx > xLow)
            {
                xLow = gx;
                iLow = ip;
            }
            if (gx >= x && gx < xHigh)
            {
                xHigh = gx;
                iHigh = ip;
            }
        }

        if (iLow < 0 || iHigh < 0)
            return;

        auto numeratorOnlyError = [&](int ip, bool isLow)
        {
            double gx, gy;
            gr->GetPoint(ip, gx, gy);
            double totalErr = isLow ? gr->GetErrorYlow(ip) : gr->GetErrorYhigh(ip);
            double commonMBErr = gy * (isLow ? mbRelSysLow : mbRelSysHigh);
            return sqrt(TMath::Max(0.0, totalErr * totalErr - commonMBErr * commonMBErr));
        };

        if (iLow == iHigh || TMath::AreEqualAbs(xLow, xHigh, 1e-12))
        {
            errLow = numeratorOnlyError(iLow, true);
            errHigh = numeratorOnlyError(iLow, false);
            return;
        }

        double w = (x - xLow) / (xHigh - xLow);
        errLow = (1.0 - w) * numeratorOnlyError(iLow, true) + w * numeratorOnlyError(iHigh, true);
        errHigh = (1.0 - w) * numeratorOnlyError(iLow, false) + w * numeratorOnlyError(iHigh, false);
    };

    TGraphErrors *gRatioNormYieldRun3Run2_stat = new TGraphErrors();
    gRatioNormYieldRun3Run2_stat->SetName("gRatioNormYieldRun3Run2_stat");
    TGraphAsymmErrors *gRatioNormYieldRun3Run2_sys = new TGraphAsymmErrors();
    gRatioNormYieldRun3Run2_sys->SetName("gRatioNormYieldRun3Run2_sys");

    for (int i = 0; i < gNormYieldKstarToMB_stat->GetN(); i++)
    {
        double x, yRun3;
        gNormYieldKstarToMB_stat->GetPoint(i, x, yRun3);
        if (x < fitXmin || x > fitXmax)
            continue;

        double yRun2 = fitNormYieldRun2->Eval(x);
        if (yRun3 <= 0 || yRun2 <= 0)
            continue;

        // Statistical uncertainty: Run-3 point and Run-2 fit covariance.
        double fitVariance = 0.0;
        if (fitResultRun2.Get() && fitResultRun2->NPar() >= 2)
            fitVariance = fitResultRun2->CovMatrix(0, 0) + x * x * fitResultRun2->CovMatrix(1, 1) + 2.0 * x * fitResultRun2->CovMatrix(0, 1);

        double eyRun3Stat = gNormYieldKstarToMB_stat->GetErrorY(i);
        double eyRun2Stat = sqrt(TMath::Max(0.0, fitVariance));
        double ratio = yRun3 / yRun2;
        double ratioStatErr = ratio * sqrt(pow(eyRun3Stat / yRun3, 2) + pow(eyRun2Stat / yRun2, 2));

        int ip = gRatioNormYieldRun3Run2_stat->GetN();
        gRatioNormYieldRun3Run2_stat->SetPoint(ip, x, ratio);
        gRatioNormYieldRun3Run2_stat->SetPointError(ip, gNormYieldKstarToMB_stat->GetErrorX(i), ratioStatErr);

        // Total systematic uncertainty: combine Run-3 and interpolated Run-2
        // relative systematic uncertainties in quadrature.
        double run3SysLowTotal = gNormYieldKstarToMB_sys->GetErrorYlow(i);
        double run3SysHighTotal = gNormYieldKstarToMB_sys->GetErrorYhigh(i);
        double run3SysLowMB = yRun3 * mbRelSysLow;
        double run3SysHighMB = yRun3 * mbRelSysHigh;
        double run3SysLow = sqrt(TMath::Max(0.0, run3SysLowTotal * run3SysLowTotal - run3SysLowMB * run3SysLowMB));
        double run3SysHigh = sqrt(TMath::Max(0.0, run3SysHighTotal * run3SysHighTotal - run3SysHighMB * run3SysHighMB));

        double run2SysLow = 0.0, run2SysHigh = 0.0;
        GetInterpolatedSys(gNormYieldKstarToMBRun2_sys, x, run2SysLow, run2SysHigh);

        double ratioSysLow = ratio * sqrt(pow(run3SysLow / yRun3, 2) + pow(run2SysLow / yRun2, 2));
        double ratioSysHigh = ratio * sqrt(pow(run3SysHigh / yRun3, 2) + pow(run2SysHigh / yRun2, 2));

        gRatioNormYieldRun3Run2_sys->SetPoint(ip, x, ratio);
        gRatioNormYieldRun3Run2_sys->SetPointError(ip, gNormYieldKstarToMB_stat->GetErrorX(i), gNormYieldKstarToMB_stat->GetErrorX(i), ratioSysLow, ratioSysHigh);
    }

    TCanvas *cNormYieldKstarToMB = new TCanvas("cNormYieldKstarToMB", "cNormYieldKstarToMB", 720, 800);
    TPad *padNormTop = new TPad("padNormTop", "padNormTop", 0, 0.30, 1, 1);
    TPad *padNormBottom = new TPad("padNormBottom", "padNormBottom", 0, 0, 1, 0.30);
    padNormTop->SetBottomMargin(0.02);
    padNormTop->SetLeftMargin(0.17);
    padNormTop->SetRightMargin(0.03);
    padNormTop->SetTopMargin(0.07);
    padNormBottom->SetTopMargin(0.03);
    padNormBottom->SetBottomMargin(0.30);
    padNormBottom->SetLeftMargin(0.17);
    padNormBottom->SetRightMargin(0.03);
    padNormTop->Draw();
    padNormBottom->Draw();

    padNormTop->cd();
    SetGraphErrorStyle(gNormYieldKstarToMB_stat);
    gNormYieldKstarToMB_stat->GetXaxis()->SetTitle("");
    gNormYieldKstarToMB_stat->GetYaxis()->SetTitle("(dN/dy) / (dN/dy)_{MB}");
    gNormYieldKstarToMB_stat->GetXaxis()->SetLimits(0, 3.9);
    gNormYieldKstarToMB_stat->GetYaxis()->SetRangeUser(0, 4.4);
    gNormYieldKstarToMB_stat->GetYaxis()->SetTitleOffset(1.55);
    gNormYieldKstarToMB_stat->GetXaxis()->SetLabelSize(0);
    gNormYieldKstarToMB_stat->SetMarkerColor(kRed + 1);
    gNormYieldKstarToMB_stat->SetLineColor(kRed + 1);
    gNormYieldKstarToMB_stat->SetMarkerStyle(20);
    gNormYieldKstarToMB_stat->SetMarkerSize(1.2);
    gNormYieldKstarToMB_stat->SetLineWidth(2);
    gNormYieldKstarToMB_stat->Draw("APE");

    gNormYieldKstarToMB_sys->SetFillStyle(0);
    gNormYieldKstarToMB_sys->SetLineColor(kRed + 1);
    gNormYieldKstarToMB_sys->SetLineWidth(2);
    gNormYieldKstarToMB_sys->Draw("5 same");
    gNormYieldKstarToMB_sysUncorr->SetFillColorAlpha(kRed + 1, 0.3);
    gNormYieldKstarToMB_sysUncorr->SetFillStyle(1001);
    gNormYieldKstarToMB_sysUncorr->SetLineColor(kRed + 1);
    gNormYieldKstarToMB_sysUncorr->Draw("E2 same");

    gNormYieldKstarToMBRun2_stat->SetMarkerStyle(21);
    gNormYieldKstarToMBRun2_stat->SetMarkerColor(kBlue + 1);
    gNormYieldKstarToMBRun2_stat->SetLineColor(kBlue + 1);
    gNormYieldKstarToMBRun2_stat->SetMarkerSize(1.2);
    gNormYieldKstarToMBRun2_stat->SetLineWidth(2);
    gNormYieldKstarToMBRun2_stat->Draw("P same");
    gNormYieldKstarToMBRun2_sys->SetFillStyle(0);
    gNormYieldKstarToMBRun2_sys->SetLineColor(kBlue + 1);
    gNormYieldKstarToMBRun2_sys->SetLineWidth(2);
    gNormYieldKstarToMBRun2_sys->Draw("5 same");
    gNormYieldKstarToMBRun2_sysUncorr->SetFillColorAlpha(kBlue + 1, 0.3);
    gNormYieldKstarToMBRun2_sysUncorr->SetFillStyle(1001);
    gNormYieldKstarToMBRun2_sysUncorr->SetLineColor(kBlue + 1);
    gNormYieldKstarToMBRun2_sysUncorr->Draw("E2 same");
    fitNormYieldRun2->Draw("same");

    for (int iModel = 0; iModel < 3; iModel++)
    {
        ScaleGraph(gNormalizeYieldKstarPythia[iModel], 2.0);
        gNormalizeYieldKstarPythia[iModel]->SetLineColor(colorsPythia[iModel]);
        gNormalizeYieldKstarPythia[iModel]->SetLineWidth(2);
        gNormalizeYieldKstarPythia[iModel]->Draw("l same");
    }

    gNormYieldKstarToMB_EPOSOn->SetLineColor(kBlue - 6);
    gNormYieldKstarToMB_EPOSOn->SetLineStyle(9);
    gNormYieldKstarToMB_EPOSOn->SetLineWidth(2);
    gNormYieldKstarToMB_EPOSOn->Draw("l same");
    gNormYieldKstarToMB_EPOSOff->SetLineColor(kGreen + 2);
    gNormYieldKstarToMB_EPOSOff->SetLineStyle(2);
    gNormYieldKstarToMB_EPOSOff->SetLineWidth(2);
    gNormYieldKstarToMB_EPOSOff->Draw("l same");

    TLegend *legNormYieldMB = new TLegend(0.42, 0.74, 0.82, 0.90);
    SetLegendStyle(legNormYieldMB);
    legNormYieldMB->AddEntry(gNormYieldKstarToMB_stat, "Run 3, #sqrt{s} = 13.6 TeV", "P");
    legNormYieldMB->AddEntry(gNormYieldKstarToMBRun2_stat, "Run 2, #sqrt{s} = 13 TeV", "P");
    legNormYieldMB->AddEntry(fitNormYieldRun2, "Run 2 linear fit", "L");
    legNormYieldMB->Draw();

    latex.DrawLatex(0.22, 0.9, "ALICE");
    latex.SetTextFont(22);
    latex.DrawLatex(0.22, 0.84, "K* (892)^{0}");
    latex.SetTextFont(42);
    latex.DrawLatex(0.22, 0.78, "|y| < 0.5");

    padNormBottom->cd();
    gRatioNormYieldRun3Run2_stat->SetTitle("");
    gRatioNormYieldRun3Run2_stat->SetMarkerStyle(20);
    gRatioNormYieldRun3Run2_stat->SetMarkerSize(1.1);
    gRatioNormYieldRun3Run2_stat->SetMarkerColor(kBlack);
    gRatioNormYieldRun3Run2_stat->SetLineColor(kBlack);
    gRatioNormYieldRun3Run2_stat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT / #LT dN_{ch}/d#eta #GT_{MB}");
    gRatioNormYieldRun3Run2_stat->GetYaxis()->SetTitle("Run 3 / Run 2");
    gRatioNormYieldRun3Run2_stat->GetXaxis()->SetLimits(0, 3.9);
    gRatioNormYieldRun3Run2_stat->GetYaxis()->SetRangeUser(0.5, 1.5);
    gRatioNormYieldRun3Run2_stat->GetXaxis()->SetTitleSize(0.11);
    gRatioNormYieldRun3Run2_stat->GetYaxis()->SetTitleSize(0.10);
    gRatioNormYieldRun3Run2_stat->GetXaxis()->SetLabelSize(0.09);
    gRatioNormYieldRun3Run2_stat->GetYaxis()->SetLabelSize(0.085);
    gRatioNormYieldRun3Run2_stat->GetXaxis()->SetTitleOffset(1.05);
    gRatioNormYieldRun3Run2_stat->GetYaxis()->SetTitleOffset(0.72);
    gRatioNormYieldRun3Run2_stat->GetYaxis()->SetNdivisions(505);
    gRatioNormYieldRun3Run2_stat->Draw("APE");

    gRatioNormYieldRun3Run2_sys->SetFillColorAlpha(kGray + 1, 0.35);
    gRatioNormYieldRun3Run2_sys->SetFillStyle(1001);
    gRatioNormYieldRun3Run2_sys->SetLineColor(kGray + 2);
    gRatioNormYieldRun3Run2_sys->SetLineWidth(1);
    gRatioNormYieldRun3Run2_sys->Draw("E2 same");
    gRatioNormYieldRun3Run2_stat->Draw("PE same");

    TLine *lineNormRun3Run2 = new TLine(0, 1, 3.5, 1);
    lineNormRun3Run2->SetLineStyle(2);
    lineNormRun3Run2->SetLineWidth(2);
    lineNormRun3Run2->Draw("same");

    cNormYieldKstarToMB->cd();
    if (isSavePlots)
        cNormYieldKstarToMB->SaveAs("Plots/NormalizedYieldToMB_Kstar_Run3Run2.png");

    //===================================================
    // Normalized dN/dy Kstar with MB, without ratio panel
    //===================================================
    TCanvas *cNormMBYieldKstar = new TCanvas("cNormMBYieldKstar", "cNormMBYieldKstar", 720, 720);
    SetCanvasStyle(cNormMBYieldKstar, 0.15, 0.03, 0.03, 0.15);
    SetGraphErrorStyle(gNormYieldKstarToMB_stat);
    gNormYieldKstarToMB_stat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT / #LT dN_{ch}/d#eta #GT_{MB}");
    gNormYieldKstarToMB_stat->GetYaxis()->SetTitle("(dN/dy) / (dN/dy)_{MB}");
    gNormYieldKstarToMB_stat->GetXaxis()->SetLimits(0, 3.9);
    gNormYieldKstarToMB_stat->GetYaxis()->SetRangeUser(0, 4.9);
    gNormYieldKstarToMB_stat->Draw("APE");

    gNormYieldKstarToMB_sys->SetFillStyle(0);
    gNormYieldKstarToMB_sys->SetLineColor(kRed + 1);
    gNormYieldKstarToMB_sys->SetLineWidth(2);
    gNormYieldKstarToMB_sys->Draw("5 same");
    gNormYieldKstarToMB_sysUncorr->SetFillColorAlpha(kRed + 1, 0.3);
    gNormYieldKstarToMB_sysUncorr->SetFillStyle(1001);
    gNormYieldKstarToMB_sysUncorr->SetLineColor(kRed + 1);
    gNormYieldKstarToMB_sysUncorr->Draw("E2 same");

    gNormYieldKstarToMBRun2_stat->SetMarkerStyle(21);
    gNormYieldKstarToMBRun2_stat->SetMarkerColor(kBlue + 1);
    gNormYieldKstarToMBRun2_stat->SetLineColor(kBlue + 1);
    gNormYieldKstarToMBRun2_stat->SetMarkerSize(1.2);
    gNormYieldKstarToMBRun2_stat->SetLineWidth(2);
    gNormYieldKstarToMBRun2_stat->Draw("P same");
    gNormYieldKstarToMBRun2_sys->SetFillStyle(0);
    gNormYieldKstarToMBRun2_sys->SetLineColor(kBlue + 1);
    gNormYieldKstarToMBRun2_sys->SetLineWidth(2);
    gNormYieldKstarToMBRun2_sys->Draw("5 same");
    gNormYieldKstarToMBRun2_sysUncorr->SetFillColorAlpha(kBlue + 1, 0.3);
    gNormYieldKstarToMBRun2_sysUncorr->SetFillStyle(1001);
    gNormYieldKstarToMBRun2_sysUncorr->SetLineColor(kBlue + 1);
    gNormYieldKstarToMBRun2_sysUncorr->Draw("E2 same");

    for (int iModel = 0; iModel < 3; iModel++)
    {
        // ScaleGraph(gNormalizeYieldKstarPythia[iModel], 2.0);
        gNormalizeYieldKstarPythia[iModel]->SetLineColor(colorsPythia[iModel]);
        gNormalizeYieldKstarPythia[iModel]->SetLineStyle(lineStylesPythia[iModel]);
        gNormalizeYieldKstarPythia[iModel]->SetLineWidth(3);
        gNormalizeYieldKstarPythia[iModel]->Draw("l same");
    }

    gNormYieldKstarToMB_EPOSOn->SetLineColor(kBlue - 6);
    gNormYieldKstarToMB_EPOSOn->SetLineStyle(9);
    gNormYieldKstarToMB_EPOSOn->SetLineWidth(3);
    gNormYieldKstarToMB_EPOSOn->Draw("l same");
    gNormYieldKstarToMB_EPOSOff->SetLineColor(kGreen + 2);
    gNormYieldKstarToMB_EPOSOff->SetLineStyle(2);
    gNormYieldKstarToMB_EPOSOff->SetLineWidth(3);
    gNormYieldKstarToMB_EPOSOff->Draw("l same");

    TLegend *legNormYieldMBOnly = new TLegend(0.42, 0.82, 0.82, 0.93);
    SetLegendStyle(legNormYieldMBOnly);
    legNormYieldMBOnly->SetTextSize(0.04);
    legNormYieldMBOnly->SetTextFont(42);
    legNormYieldMBOnly->AddEntry(gNormYieldKstarToMB_stat, "pp, #sqrt{s} = 13.6 TeV", "P");
    legNormYieldMBOnly->AddEntry(gNormYieldKstarToMBRun2_stat, "pp, #sqrt{s} = 13 TeV", "P");
    legNormYieldMBOnly->Draw();

    TLatex latexNormMB;
    latexNormMB.SetNDC();
    latexNormMB.SetTextSize(0.04);
    latexNormMB.SetTextFont(42);
    latexNormMB.DrawLatex(0.22, 0.9, "ALICE");
    latexNormMB.DrawLatex(0.22, 0.84, "|y| < 0.5");
    latexNormMB.SetTextFont(22);
    latexNormMB.DrawLatex(0.22, 0.78, "K* (892)^{0}");

    if (isSavePlots)
        cNormMBYieldKstar->SaveAs("Plots/NormalizedYieldToMB_Kstar.png");

    //===================================================
    // Normalized dN/dy Phi with MB, without ratio panel
    //===================================================
    TFile *filePhiMB = new TFile("PhiResults/PAR_plots.root", "READ");
    if (!filePhiMB || filePhiMB->IsZombie())
    {
        cout << "Error: Could not open PhiResults/PAR_plots.root" << endl;
        return;
    }
    TGraphErrors *gNormYieldPhiToMB_stat = (TGraphErrors *)filePhiMB->Get("SelfNormalizedYield/Run3_Stat");
    TGraphAsymmErrors *gNormYieldPhiToMB_sys = (TGraphAsymmErrors *)filePhiMB->Get("SelfNormalizedYield/Run3_TotalSys");
    TGraphAsymmErrors *gNormYieldPhiToMB_sysUncorr = (TGraphAsymmErrors *)filePhiMB->Get("SelfNormalizedYield/Run3_UncorrelatedSys");
    TGraphErrors *gNormYieldPhiToMBRun2_stat = (TGraphErrors *)filePhiMB->Get("SelfNormalizedYield/Run2_Stat");
    TGraphAsymmErrors *gNormYieldPhiToMBRun2_sys = (TGraphAsymmErrors *)filePhiMB->Get("SelfNormalizedYield/Run2_TotalSys");
    TGraphAsymmErrors *gNormYieldPhiToMBRun2_sysUncorr = (TGraphAsymmErrors *)filePhiMB->Get("SelfNormalizedYield/Run2_UncorrelatedSys");

    if (!gNormYieldPhiToMB_stat || !gNormYieldPhiToMB_sys || !gNormYieldPhiToMB_sysUncorr ||
        !gNormYieldPhiToMBRun2_stat || !gNormYieldPhiToMBRun2_sys || !gNormYieldPhiToMBRun2_sysUncorr)
    {
        cout << "Error: Could not retrieve normalized Phi yield graphs from file." << endl;
        return;
    }

    TCanvas *cNormMBYieldPhi = new TCanvas("cNormMBYieldPhi", "cNormMBYieldPhi", 720, 720);
    SetCanvasStyle(cNormMBYieldPhi, 0.15, 0.03, 0.03, 0.15);
    SetGraphErrorStyle(gNormYieldPhiToMB_stat);
    gNormYieldPhiToMB_stat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT / #LT dN_{ch}/d#eta #GT_{MB}");
    gNormYieldPhiToMB_stat->GetYaxis()->SetTitle("(dN/dy) / (dN/dy)_{MB}");
    gNormYieldPhiToMB_stat->GetXaxis()->SetLimits(0, 3.9);
    gNormYieldPhiToMB_stat->GetYaxis()->SetRangeUser(0, 4.9);
    gNormYieldPhiToMB_stat->Draw("APE");

    gNormYieldPhiToMB_sys->SetFillStyle(0);
    gNormYieldPhiToMB_sys->SetLineColor(kRed + 1);
    gNormYieldPhiToMB_sys->SetLineWidth(2);
    gNormYieldPhiToMB_sys->Draw("5 same");
    gNormYieldPhiToMB_sysUncorr->SetFillColorAlpha(kRed + 1, 0.3);
    gNormYieldPhiToMB_sysUncorr->SetFillStyle(1001);
    gNormYieldPhiToMB_sysUncorr->SetLineColor(kRed + 1);
    gNormYieldPhiToMB_sysUncorr->Draw("E2 same");

    gNormYieldPhiToMBRun2_stat->SetMarkerStyle(21);
    gNormYieldPhiToMBRun2_stat->SetMarkerColor(kBlue + 1);
    gNormYieldPhiToMBRun2_stat->SetLineColor(kBlue + 1);
    gNormYieldPhiToMBRun2_stat->SetMarkerSize(1.2);
    gNormYieldPhiToMBRun2_stat->SetLineWidth(2);
    gNormYieldPhiToMBRun2_stat->Draw("P same");
    gNormYieldPhiToMBRun2_sys->SetFillStyle(0);
    gNormYieldPhiToMBRun2_sys->SetLineColor(kBlue + 1);
    gNormYieldPhiToMBRun2_sys->SetLineWidth(2);
    gNormYieldPhiToMBRun2_sys->Draw("5 same");
    gNormYieldPhiToMBRun2_sysUncorr->SetFillColorAlpha(kBlue + 1, 0.3);
    gNormYieldPhiToMBRun2_sysUncorr->SetFillStyle(1001);
    gNormYieldPhiToMBRun2_sysUncorr->SetLineColor(kBlue + 1);
    gNormYieldPhiToMBRun2_sysUncorr->Draw("E2 same");

    TLegend *legNormYieldMBOnlyPhi = new TLegend(0.22, 0.68, 0.62, 0.93);
    SetLegendStyle(legNormYieldMBOnlyPhi);
    legNormYieldMBOnlyPhi->SetTextSize(0.035);
    legNormYieldMBOnlyPhi->SetTextFont(42);

    TString pythiaModelNames[3] = {"PYTHIA8 Monash", "PYTHIA8 Shoving", "PYTHIA8 Ropes"};

    for (int iModel = 0; iModel < 3; iModel++)
    {
        // ScaleGraph(gNormalizeYieldPhiPythia[iModel], 2.0);
        gNormalizeYieldPhiPythia[iModel]->SetLineColor(colorsPythia[iModel]);
        gNormalizeYieldPhiPythia[iModel]->SetLineStyle(lineStylesPythia[iModel]);
        gNormalizeYieldPhiPythia[iModel]->SetLineWidth(3);
        gNormalizeYieldPhiPythia[iModel]->Draw("l same");
        legNormYieldMBOnlyPhi->AddEntry(gNormalizeYieldPhiPythia[iModel], Form("%s", pythiaModelNames[iModel].Data()), "L");
    }

    gNormYieldPhiToMB_EPOSOn->SetLineColor(kBlue - 6);
    gNormYieldPhiToMB_EPOSOn->SetLineStyle(9);
    gNormYieldPhiToMB_EPOSOn->SetLineWidth(3);
    gNormYieldPhiToMB_EPOSOn->Draw("l same");
    gNormYieldPhiToMB_EPOSOff->SetLineColor(kGreen + 2);
    gNormYieldPhiToMB_EPOSOff->SetLineStyle(2);
    gNormYieldPhiToMB_EPOSOff->SetLineWidth(3);
    gNormYieldPhiToMB_EPOSOff->Draw("l same");

    legNormYieldMBOnlyPhi->AddEntry(gNormYieldPhiToMB_EPOSOn, "EPOS, UrQMD on", "L");
    legNormYieldMBOnlyPhi->AddEntry(gNormYieldPhiToMB_EPOSOff, "EPOS, UrQMD off", "L");
    legNormYieldMBOnlyPhi->Draw();

    latexNormMB.SetTextFont(22);
    latexNormMB.DrawLatex(0.65, 0.78, "#phi (1020)");

    cNormMBYieldPhi->SaveAs("Plots/NormalizedYieldToMB_Phi.png");

    //===================================================
    // ============== Normalized Mean pT Kstar ==========
    //===================================================

    TGraphErrors *gNormMeanPtKstar[3], *gNormMeanPtKstarRun2[3];

    for (int i = 0; i < 3; i++)
    {
        gNormMeanPtKstar[i] = NormalizeGraphByX(gMPtKstar[i], Form("gNormMeanPtKstar_%d", i));
        gNormMeanPtKstarRun2[i] = NormalizeGraphByX(gMPtKstarRun2[i], Form("gNormMeanPtKstarRun2_%d", i));
    }

    TCanvas *cNormMeanPtKstar = new TCanvas("cNormMeanPtKstar", "cNormMeanPtKstar", 720, 720);
    SetCanvasStyle(cNormMeanPtKstar, 0.15, 0.03, 0.03, 0.15);
    SetGraphErrorStyle(gNormMeanPtKstar[0]);
    gNormMeanPtKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gNormMeanPtKstar[0]->GetYaxis()->SetTitle("#LT p_{T} #GT / #LT dN_{ch}/d#eta #GT");
    gNormMeanPtKstar[0]->GetXaxis()->SetLimits(0, 27);
    gNormMeanPtKstar[0]->SetMarkerColor(kRed + 1);
    gNormMeanPtKstar[0]->SetLineColor(kRed + 1);
    gNormMeanPtKstar[0]->SetLineWidth(3);
    gNormMeanPtKstar[0]->GetYaxis()->SetRangeUser(0.035, 0.34);
    gNormMeanPtKstar[0]->Draw("APE");

    gNormMeanPtKstar[1]->SetFillStyle(0);
    gNormMeanPtKstar[1]->SetLineColor(kRed + 1);
    gNormMeanPtKstar[1]->SetLineWidth(3);
    gNormMeanPtKstar[1]->Draw("5 same");

    gNormMeanPtKstar[2]->SetFillStyle(1001);
    gNormMeanPtKstar[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gNormMeanPtKstar[2]->Draw("E2 same");

    gNormMeanPtKstarRun2[0]->SetMarkerStyle(21);
    gNormMeanPtKstarRun2[0]->SetMarkerColor(kBlue + 1);
    gNormMeanPtKstarRun2[0]->SetLineColor(kBlue + 1);
    gNormMeanPtKstarRun2[0]->Draw("P same");

    gNormMeanPtKstarRun2[1]->SetFillStyle(0);
    gNormMeanPtKstarRun2[1]->SetLineColor(kBlue + 1);
    gNormMeanPtKstarRun2[1]->SetLineWidth(3);
    gNormMeanPtKstarRun2[1]->Draw("5 same");

    gNormMeanPtKstarRun2[2]->SetFillStyle(1001);
    gNormMeanPtKstarRun2[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gNormMeanPtKstarRun2[2]->Draw("E2 same");

    // Models
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gr = NormalizeGraphByX(gMeanPt[model][kKstar], Form("gNormModelMeanPt_%d", model));
        gr->Draw("l same");
    }

    // EPOS4
    gNormEPOS_KstarMeanpT_UrQMDOff->SetLineColor(kGreen + 2);
    gNormEPOS_KstarMeanpT_UrQMDOff->SetLineStyle(2);
    gNormEPOS_KstarMeanpT_UrQMDOff->SetLineWidth(3);
    gNormEPOS_KstarMeanpT_UrQMDOff->Draw("l same");

    // EPOS4 + UrQMD
    gNormEPOS_KstarMeanpT_UrQMDOn->SetLineColor(kBlue - 6);
    gNormEPOS_KstarMeanpT_UrQMDOn->SetLineStyle(9);
    gNormEPOS_KstarMeanpT_UrQMDOn->SetLineWidth(3);
    gNormEPOS_KstarMeanpT_UrQMDOn->Draw("l same");

    // latex.DrawLatex(0.22, 0.9, "ALICE");
    // latex.SetTextFont(22);
    // latex.DrawLatex(0.22, 0.84, "K* (892)^{0}");
    // latex.SetTextFont(42);
    // latex.DrawLatex(0.22, 0.78, "|y| < 0.5");

    legend2->Draw();

    if (isSavePlots)
        cNormMeanPtKstar->SaveAs("Plots/NormalizedMeanPt_Kstar.png");

    //===================================================
    //  ================<pT> Phi======================
    //===================================================
    TCanvas *cMeanPtPhi = new TCanvas("cMeanPtPhi", "cMeanPtPhi", 720, 720);
    SetCanvasStyle(cMeanPtPhi, 0.15, 0.03, 0.03, 0.15);
    canvas_style(cMeanPtPhi, pad1Size, pad2Size);
    cMeanPtPhi->cd(1);
    SetGraphErrorStyle(gMPtPhiRun3[0]);
    gMPtPhiRun3[0]->SetTitle(0);
    gMPtPhiRun3[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMPtPhiRun3[0]->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/#it{c})");
    gMPtPhiRun3[0]->GetXaxis()->SetLimits(0, 27);
    gMPtPhiRun3[0]->GetYaxis()->SetRangeUser(0.58, 1.76);
    gMPtPhiRun3[0]->SetMarkerColor(kRed);
    gMPtPhiRun3[0]->SetLineColor(kRed);
    gMPtPhiRun3[0]->SetLineWidth(3);
    gMPtPhiRun3[0]->GetYaxis()->SetTitleOffset(1.3);
    gMPtPhiRun3[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    gMPtPhiRun3[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    gMPtPhiRun3[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    gMPtPhiRun3[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    gMPtPhiRun3[0]->Draw("APE");
    gMPtPhiRun3[1]->SetFillStyle(0);
    gMPtPhiRun3[1]->SetLineColor(kRed);
    gMPtPhiRun3[1]->SetLineWidth(3);
    gMPtPhiRun3[1]->Draw("5 same");
    gMPtPhiRun3[2]->SetFillStyle(1001);
    gMPtPhiRun3[2]->SetFillColorAlpha(kRed, 0.3);
    gMPtPhiRun3[2]->Draw("E2 same");
    gMPtPhiRun2[0]->SetMarkerStyle(21);
    gMPtPhiRun2[0]->SetMarkerColor(kBlue);
    gMPtPhiRun2[0]->SetLineColor(kBlue);
    gMPtPhiRun2[0]->SetLineWidth(3);
    gMPtPhiRun2[0]->Draw("P same");
    gMPtPhiRun2[1]->SetLineColor(kBlue);
    gMPtPhiRun2[1]->SetFillStyle(0);
    gMPtPhiRun2[1]->SetLineWidth(3);
    gMPtPhiRun2[1]->Draw("5 same");
    gMPtPhiRun2[2]->SetFillStyle(1001);
    gMPtPhiRun2[2]->SetFillColorAlpha(kBlue, 0.3);
    gMPtPhiRun2[2]->Draw("E2 same");

    for (int i = 0; i < totalPointsEPOS_Phi; i++)
    {
        gEPOS_PhiMeanpT_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_PhiMeanpT_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
    }

    gEPOS_PhiMeanpT_UrQMDOff->SetLineColor(kGreen + 2);
    gEPOS_PhiMeanpT_UrQMDOff->SetLineStyle(2);
    gEPOS_PhiMeanpT_UrQMDOff->SetLineWidth(3);
    gEPOS_PhiMeanpT_UrQMDOff->Draw("l same");
    gEPOS_PhiMeanpT_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_PhiMeanpT_UrQMDOn->SetLineStyle(9);
    gEPOS_PhiMeanpT_UrQMDOn->SetLineWidth(3);
    gEPOS_PhiMeanpT_UrQMDOn->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMeanPt[model][kKstar]->Draw("l same");
    }

    // legend->Draw();
    legend2->Draw();
    latex.SetTextSize(0.03);
    latex.SetTextFont(22);
    latex.DrawLatex(0.8, 0.90, "#Phi (1020)");

    cMeanPtPhi->cd(2);
    cMeanPtPhi->cd(2)->SetGridy();

    TGraphErrors *gRatioDataModelPhiMeanPt[kNModels];
    counter = 0;
    for (auto model : modelsToPlot)
    {
        gRatioDataModelPhiMeanPt[model] = MakeDataModelRatio(gMPtPhiRun3[0], gMeanPt[model][kPhi]);
        gRatioDataModelPhiMeanPt[model]->SetLineColor(colorsPythia[counter]);
        gRatioDataModelPhiMeanPt[model]->SetMarkerColor(colorsPythia[counter]);
        gRatioDataModelPhiMeanPt[model]->SetMarkerStyle(20);
        gRatioDataModelPhiMeanPt[model]->SetLineWidth(2);
        gRatioDataModelPhiMeanPt[model]->GetYaxis()->SetTitle("Data / Model");
        gRatioDataModelPhiMeanPt[model]->GetYaxis()->SetTitleOffset(0.5);
        gRatioDataModelPhiMeanPt[model]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
        gRatioDataModelPhiMeanPt[model]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
        gRatioDataModelPhiMeanPt[model]->GetYaxis()->SetLabelSize(0.035 / pad2Size);
        gRatioDataModelPhiMeanPt[model]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        gRatioDataModelPhiMeanPt[model]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        gRatioDataModelPhiMeanPt[model]->GetXaxis()->SetLimits(0, 27);
        gRatioDataModelPhiMeanPt[model]->GetYaxis()->SetRangeUser(0.9, 1.5);
        gRatioDataModelPhiMeanPt[model]->GetYaxis()->SetNdivisions(508);
        if (counter == 0)
        {
            gRatioDataModelPhiMeanPt[model]->Draw("APE");
        }
        else
        {
            gRatioDataModelPhiMeanPt[model]->Draw("PE same");
        }
        counter++;
    }

    TGraphErrors *gRatioDataModelPhiMeanPtEPOSON = MakeDataModelRatio(gMPtPhiRun3[0], gEPOS_PhiMeanpT_UrQMDOn);
    TGraphErrors *gRatioDataModelPhiMeanPtEPOSOFF = MakeDataModelRatio(gMPtPhiRun3[0], gEPOS_PhiMeanpT_UrQMDOff);

    gRatioDataModelPhiMeanPtEPOSON->SetLineColor(kBlue - 6);
    gRatioDataModelPhiMeanPtEPOSON->SetMarkerColor(kBlue - 6);
    gRatioDataModelPhiMeanPtEPOSON->SetMarkerStyle(20);
    gRatioDataModelPhiMeanPtEPOSON->SetLineWidth(2);
    gRatioDataModelPhiMeanPtEPOSON->Draw("PE same");
    gRatioDataModelPhiMeanPtEPOSOFF->SetLineColor(kGreen + 2);
    gRatioDataModelPhiMeanPtEPOSOFF->SetMarkerColor(kGreen + 2);
    gRatioDataModelPhiMeanPtEPOSOFF->SetMarkerStyle(20);
    gRatioDataModelPhiMeanPtEPOSOFF->SetLineWidth(2);
    gRatioDataModelPhiMeanPtEPOSOFF->Draw("PE same");

    TLine *lineMeanPtPhi = new TLine(0, 1, 27, 1);
    lineMeanPtPhi->SetLineStyle(2);
    lineMeanPtPhi->SetLineWidth(2);
    lineMeanPtPhi->Draw("same");

    cMeanPtPhi->SaveAs("Plots/QA/MeanPt_Phi.png");

    //====================================================
    //   ================Kstar/K Ratio ===================
    //====================================================

    TGraphErrors *gdNdy_Kstar_On = GetGraph(fEPOS_New2, "g_dNdy_opt1");
    TGraphErrors *gdNdy_Kstar_Off = GetGraph(fEPOS_New2, "g_dNdy_opt0");
    TGraphErrors *gdNdy_K0s_On = GetGraph(fEPOS_New2, "g_stable_dNdy_opt1_K0s");
    TGraphErrors *gdNdy_K0s_Off = GetGraph(fEPOS_New2, "g_stable_dNdy_opt0_K0s");
    TGraphErrors *gdNdy_Kaon_On = GetGraph(fEPOS_New2, "g_stable_dNdy_opt1_kaon");
    TGraphErrors *gdNdy_Kaon_Off = GetGraph(fEPOS_New2, "g_stable_dNdy_opt0_kaon");
    TGraphErrors *gdNdy_Pion_On = GetGraph(fEPOS_New2, "g_stable_dNdy_opt1_pi");
    TGraphErrors *gdNdy_Pion_Off = GetGraph(fEPOS_New2, "g_stable_dNdy_opt0_pi");
    TGraphErrors *gdNdy_Phi_On = GetGraph(fEPOS_New2, "g_phi_dNdy_opt1");
    TGraphErrors *gdNdy_Phi_Off = GetGraph(fEPOS_New2, "g_phi_dNdy_opt0");

    TGraphErrors *gKstarKaRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarKaRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gKstarPiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarPiRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gKstarK0sRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarK0sRatio_EPOSNew_UrQMDOn = new TGraphErrors();

    TGraphErrors *gPhiKaRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gPhiKaRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gPhiPiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gPhiPiRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gPhiK0sRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gPhiK0sRatio_EPOSNew_UrQMDOn = new TGraphErrors();

    TGraphErrors *gKstarPhiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarPhiRatio_EPOSNew_UrQMDOn = new TGraphErrors();

    int graphPointsEPOSNew = gdNdy_Kaon_On->GetN();
    for (int i = 0; i < graphPointsEPOSNew; i++)
    {
        double xK0s, yK0s, yKaon, yPion, yKstarOff, yKstarOn;
        double yK0sON, yKaonON, yPionON, yPhiOff, yPhiOn;

        gdNdy_Phi_Off->GetPoint(i, xK0s, yPhiOff);
        gdNdy_Phi_On->GetPoint(i, xK0s, yPhiOn);
        cout << "xPhi: " << xK0s << ", yPhiOff: " << yPhiOff << ", yPhiOn: " << yPhiOn << endl;

        gdNdy_Kstar_Off->GetPoint(i, xK0s, yKstarOff);
        gdNdy_Kstar_On->GetPoint(i, xK0s, yKstarOn);
        cout << "xKstar: " << xK0s << ", yKstarOff: " << yKstarOff << ", yKstarOn: " << yKstarOn << endl;

        gdNdy_K0s_Off->GetPoint(i, xK0s, yK0s);
        gdNdy_Kaon_Off->GetPoint(i, xK0s, yKaon);
        gdNdy_Pion_Off->GetPoint(i, xK0s, yPion);

        gdNdy_K0s_On->GetPoint(i, xK0s, yK0sON);
        gdNdy_Kaon_On->GetPoint(i, xK0s, yKaonON);
        gdNdy_Pion_On->GetPoint(i, xK0s, yPionON);

        gKstarKaRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yKaon);
        gKstarKaRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yKaon);
        gKstarPiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yPion);
        gKstarPiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yPion);
        gKstarK0sRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yK0s);
        gKstarK0sRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yK0s);

        gKstarPhiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yPhiOff);
        gKstarPhiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yPhiOn);
    }

    int graphPointsEPOSNewPhi = gdNdy_Phi_On->GetN();
    for (int i = 0; i < graphPointsEPOSNewPhi; i++)
    {
        double xK0s, yK0s, yKaon, yPion, yPhiOff, yPhiOn;
        double yK0sON, yKaonON, yPionON;
        gdNdy_Phi_Off->GetPoint(i, xK0s, yPhiOff);
        gdNdy_Phi_On->GetPoint(i, xK0s, yPhiOn);

        gdNdy_K0s_Off->GetPoint(i, xK0s, yK0s);
        gdNdy_Kaon_Off->GetPoint(i, xK0s, yKaon);
        gdNdy_Pion_Off->GetPoint(i, xK0s, yPion);

        gdNdy_K0s_On->GetPoint(i, xK0s, yK0sON);
        gdNdy_Kaon_On->GetPoint(i, xK0s, yKaonON);
        gdNdy_Pion_On->GetPoint(i, xK0s, yPionON);

        gPhiKaRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yKaon);
        gPhiKaRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yKaon);
        gPhiPiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yPion);
        gPhiPiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yPion);
        gPhiK0sRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yK0s);
        gPhiK0sRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yK0s);
    }

    TGraphErrors *gKstarKaRatio_UrQMDOffSmooth = SmoothGraph(gKstarKaRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarKaRatio_UrQMDOnSmooth = SmoothGraph(gKstarKaRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gKstarPiRatio_UrQMDOffSmooth = SmoothGraph(gKstarPiRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarPiRatio_UrQMDOnSmooth = SmoothGraph(gKstarPiRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gKstarK0sRatio_UrQMDOffSmooth = SmoothGraph(gKstarK0sRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarK0sRatio_UrQMDOnSmooth = SmoothGraph(gKstarK0sRatio_EPOSNew_UrQMDOn);

    TGraphErrors *gPhiKaRatio_UrQMDOffSmooth = SmoothGraph(gPhiKaRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gPhiKaRatio_UrQMDOnSmooth = SmoothGraph(gPhiKaRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gPhiPiRatio_UrQMDOffSmooth = SmoothGraph(gPhiPiRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gPhiPiRatio_UrQMDOnSmooth = SmoothGraph(gPhiPiRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gPhiK0sRatio_UrQMDOffSmooth = SmoothGraph(gPhiK0sRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gPhiK0sRatio_UrQMDOnSmooth = SmoothGraph(gPhiK0sRatio_EPOSNew_UrQMDOn);

    TGraphErrors *gKstarPhiRatio_UrQMDOffSmooth = SmoothGraph(gKstarPhiRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarPhiRatio_UrQMDOnSmooth = SmoothGraph(gKstarPhiRatio_EPOSNew_UrQMDOn);

    TCanvas *cRatioKstarKaon = new TCanvas("cRatioKstarKaon", "cRatioKstarKaon", 720, 720);
    SetCanvasStyle(cRatioKstarKaon, 0.15, 0.03, 0.03, 0.15);

    for (int i3 = 0; i3 < TotalPointsRun3; i3++)
    {
        double xkstar, ykstar, ysyskstar, yuncorrkstar;
        double xkaon, ykaon, ysyskaon, yuncorrkaon;
        gMYieldKstar[0]->GetPoint(i3, xkstar, ykstar);
        ysyskstar = gMYieldKstar[1]->GetErrorY(i3);
        yuncorrkstar = gMYieldKstar[2]->GetErrorY(i3);

        gMYieldKaon[0]->GetPoint(i3, xkaon, ykaon);
        ysyskaon = gMYieldKaon[1]->GetErrorY(i3);
        yuncorrkaon = gMYieldKaon[2]->GetErrorY(i3);

        cout << "Point " << i3 << ": K*0 yield = " << ykstar << " +/- " << ysyskstar << " (sys) +/- " << yuncorrkstar << " (uncorr)" << endl;
        cout << "Point " << i3 << ": Kaon yield = " << ykaon << " +/- " << ysyskaon << " (sys) +/- " << yuncorrkaon << " (uncorr)" << endl;
        double ratio = ykstar / ykaon;
        double ratioErrorsys = ratio * sqrt(pow(yuncorrkstar / ykstar, 2) + pow(ysyskaon / ykaon, 2));
        double ratioErroruncorr = ratio * sqrt(pow(yuncorrkstar / ykstar, 2) + pow(yuncorrkaon / ykaon, 2));
        cout << "Point " << i3 << ": K*0/Kaon ratio = " << ratio << " +/- " << ratioErrorsys << " (sys) +/- " << ratioErroruncorr << " (uncorr)" << endl;
        cout << endl;
    }

    TGraphErrors **gRatioKstarKaon = MakeRatioUncorr(gMYieldKstar, gMYieldKaon, 1.0, 3);
    gRatioKstarKaon[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioKstarKaon[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphErrorStyle(gRatioKstarKaon[0]);
    gRatioKstarKaon[0]->GetYaxis()->SetRangeUser(0.245, 0.44);
    gRatioKstarKaon[0]->GetXaxis()->SetLimits(1.8, 23.5);
    gRatioKstarKaon[0]->SetMarkerColor(kRed);
    gRatioKstarKaon[0]->SetLineColor(kRed);
    gRatioKstarKaon[0]->GetYaxis()->SetNdivisions(505);
    gRatioKstarKaon[0]->Draw("APE");
    gRatioKstarKaon[1]->SetFillStyle(0);
    gRatioKstarKaon[1]->SetLineColor(kRed);
    gRatioKstarKaon[1]->Draw("5 same");
    gRatioKstarKaon[2]->SetFillStyle(1001);
    gRatioKstarKaon[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioKstarKaon[2]->Draw("E2 same");

    // gKstarKaRatio_13TeV[0]->SetMarkerStyle(21);
    // gKstarKaRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gKstarKaRatio_13TeV[0]->SetLineColor(kBlue);
    // gKstarKaRatio_13TeV[0]->Draw("P same");
    // gKstarKaRatio_13TeV[1]->SetLineColor(kBlue);
    // gKstarKaRatio_13TeV[1]->SetFillStyle(0);
    // gKstarKaRatio_13TeV[1]->Draw("5 same");

    gKstarKaRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    gKstarKaRatio_UrQMDOffSmooth->SetLineStyle(2);
    gKstarKaRatio_UrQMDOffSmooth->SetLineWidth(3);
    gKstarKaRatio_UrQMDOffSmooth->Draw("LX same");
    gKstarKaRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    gKstarKaRatio_UrQMDOnSmooth->SetLineStyle(9);
    gKstarKaRatio_UrQMDOnSmooth->SetLineWidth(3);
    gKstarKaRatio_UrQMDOnSmooth->Draw("LX same");

    ////Calculation of suppression significance between HM and LM for K*0/Kaon ratio
    int totalPointsKstarKaonRatio = gRatioKstarKaon[0]->GetN();
    double xLM, xHM, yHM, yLM, yHMError, yLMError, yHMStatError, yLMStatError, yLMUncorr, yHMUncorr;
    gRatioKstarKaon[1]->GetPoint(totalPointsKstarKaonRatio - 1, xLM, yLM);
    gRatioKstarKaon[1]->GetPoint(0, xHM, yHM);
    yLMError = gRatioKstarKaon[1]->GetErrorY(totalPointsKstarKaonRatio - 1);
    yHMError = gRatioKstarKaon[1]->GetErrorY(0);
    yLMStatError = gRatioKstarKaon[0]->GetErrorY(totalPointsKstarKaonRatio - 1);
    yHMStatError = gRatioKstarKaon[0]->GetErrorY(0);
    yLMUncorr = gRatioKstarKaon[2]->GetErrorY(totalPointsKstarKaonRatio - 1);
    yHMUncorr = gRatioKstarKaon[2]->GetErrorY(0);

    double sigmaDiffHM_LM = abs(yHM - yLM) / sqrt(pow(yHMError, 2) + pow(yLMError, 2) + pow(yHMStatError, 2) + pow(yLMStatError, 2));
    double sigmaDiffHM_LM_uncorr = abs(yHM - yLM) / sqrt(pow(yHMUncorr, 2) + pow(yLMUncorr, 2));
    double sigmaDiffHM_LM_uncorrStat = abs(yHM - yLM) / sqrt(pow(yHMStatError, 2) + pow(yLMStatError, 2) + pow(yHMUncorr, 2) + pow(yLMUncorr, 2));
    cout << "K*0/Kaon ratio difference between HM and LM: " << sigmaDiffHM_LM << " sigma" << endl;
    cout << "K*0/Kaon ratio difference between HM and LM (uncorrelated): " << sigmaDiffHM_LM_uncorr << " sigma" << endl;
    cout << "K*0/Kaon ratio difference between HM and LM (uncorrelated, statistical): " << sigmaDiffHM_LM_uncorrStat << " sigma" << endl;

    TLegend *legendRatio = new TLegend(0.65, 0.85, 0.9, 0.92);
    SetLegendStyle(legendRatio);
    legendRatio->SetTextSize(0.04);
    legendRatio->AddEntry(gRatioKstarKaon[0], "K*^{0} / K", "p");

    TLegend *legendRatio2 = new TLegend(0.22, 0.66, 0.5, 0.92);
    SetLegendStyle(legendRatio2);
    legendRatio2->SetTextSize(0.03);
    legendRatio2->AddEntry(gKstarKaRatio_UrQMDOffSmooth, "EPOS4", "l");
    legendRatio2->AddEntry(gKstarKaRatio_UrQMDOnSmooth, "EPOS4 + UrQMD", "l");

    counter = 0;
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarKaModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKaon], true);
        ScaleGraph(gRatioKstarKaModel, 2.0); // because Kstar yield was already average during dN/dy plotting, but kaon was remaining
        setStyle(gRatioKstarKaModel, colorsPythia[counter], lineStylesPythia[counter]);
        gRatioKstarKaModel->SetLineWidth(3);
        gRatioKstarKaModel->Draw("l same");
        legendRatio2->AddEntry(gRatioKstarKaModel, modelLabel[model], "l");
        counter++;
    }

    legendRatio->Draw();
    // legendRatio2->Draw();
    latex.DrawLatex(0.31, 0.9, "ALICE");
    latex.DrawLatex(0.31, 0.85, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.31, 0.80, "|y| < 0.5");
    if (isSavePlots)
    {
        cRatioKstarKaon->SaveAs("Plots/QA/Ratio_KstarKaon_Run3.png");
    }

    //====================================================
    //   ================Kstar/Phi Ratio ===================
    //====================================================
    TCanvas *cRatioKstarPhi = new TCanvas("cRatioKstarPhi", "cRatioKstarPhi", 720, 720);
    SetCanvasStyle(cRatioKstarPhi, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors **gRatioKstarPhi = MakeRatioUncorr(gMYieldKstar, gMYieldPhiRun3, 1.0, 3);
    // cout << "Number of points in kstar run2 yield " << gKstar_Yield_13TeV[0]->GetN() << endl;
    // cout << "Number of points in phi run2 yield " << gMYieldPhiRun2[0]->GetN() << endl;
    // TGraphErrors **gRatioKstarPhiRun2 = MakeRatioUncorr(gKstar_Yield_13TeV, gMYieldPhiRun2, 1.0, 3);
    TGraphErrors *gRatioKstarPhiRun2[3];
    CalculateKstarPhiRatioRun2(gKstar_Yield_13TeV, gMYieldPhiRun2, gRatioKstarPhiRun2);

    gRatioKstarPhi[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioKstarPhi[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphErrorStyle(gRatioKstarPhi[0]);
    gRatioKstarPhi[0]->GetYaxis()->SetRangeUser(1.2, 6.15);
    gRatioKstarPhi[0]->GetXaxis()->SetLimits(1.2, 27.5);
    gRatioKstarPhi[0]->SetMarkerColor(kRed);
    gRatioKstarPhi[0]->SetLineColor(kRed);
    gRatioKstarPhi[0]->GetYaxis()->SetNdivisions(505);
    gRatioKstarPhi[0]->Draw("APE");
    gRatioKstarPhi[1]->SetFillStyle(0);
    gRatioKstarPhi[1]->SetLineColor(kRed);
    gRatioKstarPhi[1]->Draw("5 same");
    gRatioKstarPhi[2]->SetFillStyle(1001);
    gRatioKstarPhi[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioKstarPhi[2]->Draw("E2 same");
    gRatioKstarPhiRun2[0]->SetMarkerStyle(21);
    gRatioKstarPhiRun2[0]->SetMarkerColor(kBlue);
    gRatioKstarPhiRun2[0]->SetLineColor(kBlue);
    gRatioKstarPhiRun2[0]->SetLineWidth(3);
    gRatioKstarPhiRun2[0]->Draw("PE same");
    gRatioKstarPhiRun2[1]->SetLineColor(kBlue);
    gRatioKstarPhiRun2[1]->SetFillStyle(0);
    gRatioKstarPhiRun2[1]->SetLineWidth(3);
    gRatioKstarPhiRun2[1]->Draw("5 same");

    // calculate the significance of the difference between HM and LM for K*0/Phi ratio
    int totalPointsKstarPhiRatio = gRatioKstarPhi[0]->GetN();
    double xLM_KstarPhi, xHM_KstarPhi, yHM_KstarPhi, yLM_KstarPhi, yHMError_KstarPhi, yLMError_KstarPhi, yHMStatError_KstarPhi, yLMStatError_KstarPhi, yLMUncorr_KstarPhi, yHMUncorr_KstarPhi;
    gRatioKstarPhi[1]->GetPoint(totalPointsKstarPhiRatio - 1, xLM_KstarPhi, yLM_KstarPhi);
    gRatioKstarPhi[1]->GetPoint(0, xHM_KstarPhi, yHM_KstarPhi);
    yLMError_KstarPhi = gRatioKstarPhi[1]->GetErrorY(totalPointsKstarPhiRatio - 1);
    yHMError_KstarPhi = gRatioKstarPhi[1]->GetErrorY(0);
    yLMStatError_KstarPhi = gRatioKstarPhi[0]->GetErrorY(totalPointsKstarPhiRatio - 1);
    yHMStatError_KstarPhi = gRatioKstarPhi[0]->GetErrorY(0);
    yLMUncorr_KstarPhi = gRatioKstarPhi[2]->GetErrorY(totalPointsKstarPhiRatio - 1);
    yHMUncorr_KstarPhi = gRatioKstarPhi[2]->GetErrorY(0);
    double significance_KstarPhi = fabs(yHM_KstarPhi - yLM_KstarPhi) / sqrt(pow(yHMError_KstarPhi, 2) + pow(yLMError_KstarPhi, 2));
    cout << "Significance of K*0/Phi ratio difference: " << significance_KstarPhi << endl;
    double significance_KstarPhi_uncorr = fabs(yHM_KstarPhi - yLM_KstarPhi) / sqrt(pow(yHMUncorr_KstarPhi, 2) + pow(yLMUncorr_KstarPhi, 2));
    cout << "Significance of K*0/Phi ratio difference (uncorrelated): " << significance_KstarPhi_uncorr << endl;
    double significance_KstarPhi_uncorrStat = fabs(yHM_KstarPhi - yLM_KstarPhi) / sqrt(pow(yHMStatError_KstarPhi, 2) + pow(yLMStatError_KstarPhi, 2) + pow(yHMUncorr_KstarPhi, 2) + pow(yLMUncorr_KstarPhi, 2));
    cout << "Significance of K*0/Phi ratio difference (uncorrelated, statistical): " << significance_KstarPhi_uncorrStat << endl;

    // gKstarPhiRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    // gKstarPhiRatio_UrQMDOffSmooth->SetLineStyle(2);
    // gKstarPhiRatio_UrQMDOffSmooth->SetLineWidth(3);
    // gKstarPhiRatio_UrQMDOffSmooth->Draw("LX same");
    // gKstarPhiRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    // gKstarPhiRatio_UrQMDOnSmooth->SetLineStyle(9);
    // gKstarPhiRatio_UrQMDOnSmooth->SetLineWidth(3);
    // gKstarPhiRatio_UrQMDOnSmooth->Draw("LX same");

    TLegend *legendRatioKstarPhi = new TLegend(0.52, 0.73, 0.9, 0.94);
    SetLegendStyle(legendRatioKstarPhi);
    legendRatioKstarPhi->SetTextSize(0.03);
    legendRatioKstarPhi->AddEntry(gRatioKstarPhi[0], "K*^{0} / #Phi (Run 3)", "p");
    legendRatioKstarPhi->AddEntry(gRatioKstarPhiRun2[0], "K*^{0} / #Phi (Run 2)", "p");
    // legendRatioKstarPhi->AddEntry(gKstarPhiRatio_UrQMDOffSmooth, "EPOS4", "l");
    // legendRatioKstarPhi->AddEntry(gKstarPhiRatio_UrQMDOnSmooth, "EPOS4 + UrQMD", "l");

    // counter = 0;
    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioKstarPhiModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kPhi], true);
    //     setStyle(gRatioKstarPhiModel, colorsPythia[counter], lineStylesPythia[counter]);
    //     gRatioKstarPhiModel->SetLineWidth(3);
    //     gRatioKstarPhiModel->Draw("l same");
    //     counter++;
    //     legendRatioKstarPhi->AddEntry(gRatioKstarPhiModel, modelLabel[model], "l");
    // }

    // legendRatio->Draw();
    // legendRatio2->Draw();
    latex.DrawLatex(0.22, 0.9, "ALICE");
    latex.DrawLatex(0.22, 0.85, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.22, 0.80, "|y| < 0.5");
    // latex.DrawLatex(0.22, 0.75, "K*^{0} / #Phi");

    legendRatioKstarPhi->Draw();

    if (isSavePlots)
    {
        cRatioKstarPhi->SaveAs("Plots/QA/Ratio_KstarPhi_Run3.png");
    }

    //====================================================
    //   ================Phi/K Ratio ===================
    //====================================================
    TCanvas *cRatioPhiK = new TCanvas("cRatioPhiK", "cRatioPhiK", 720, 720);
    SetCanvasStyle(cRatioPhiK, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors **gRatioPhiK = MakeRatioUncorr(gMYieldPhiRun3, gMYieldKaon, 1.0, 3);
    gRatioPhiK[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioPhiK[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphErrorStyle(gRatioPhiK[0]);
    gRatioPhiK[0]->GetYaxis()->SetRangeUser(0.07, 0.17);
    gRatioPhiK[0]->GetXaxis()->SetLimits(1.8, 23.5);
    // gRatioPhiK[0]->GetXaxis()->SetLimits(1.8, 26.5);
    gRatioPhiK[0]->SetMarkerColor(kRed);
    gRatioPhiK[0]->SetLineColor(kRed);
    gRatioPhiK[0]->GetYaxis()->SetNdivisions(505);
    gRatioPhiK[0]->Draw("APE");
    gRatioPhiK[1]->SetFillStyle(0);
    gRatioPhiK[1]->SetLineColor(kRed);
    gRatioPhiK[1]->Draw("5 same");
    gRatioPhiK[2]->SetFillStyle(1001);
    gRatioPhiK[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioPhiK[2]->Draw("E2 same");

    // gPhiKaRatio_13TeV[0]->SetMarkerStyle(21);
    // gPhiKaRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gPhiKaRatio_13TeV[0]->SetLineColor(kBlue);
    // gPhiKaRatio_13TeV[0]->Draw("P same");
    // gPhiKaRatio_13TeV[1]->SetLineColor(kBlue);
    // gPhiKaRatio_13TeV[1]->SetFillStyle(0);
    // gPhiKaRatio_13TeV[1]->Draw("5 same");

    gPhiKaRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    gPhiKaRatio_UrQMDOffSmooth->SetLineStyle(2);
    gPhiKaRatio_UrQMDOffSmooth->SetLineWidth(3);
    gPhiKaRatio_UrQMDOffSmooth->Draw("LX same");
    gPhiKaRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    gPhiKaRatio_UrQMDOnSmooth->SetLineStyle(9);
    gPhiKaRatio_UrQMDOnSmooth->SetLineWidth(3);
    gPhiKaRatio_UrQMDOnSmooth->Draw("LX same");

    counter = 0;
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiKaModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kKaon], true);
        setStyle(gRatioPhiKaModel, colorsPythia[counter], lineStylesPythia[counter]);
        ScaleGraph(gRatioPhiKaModel, 2.0); // because Phi yield was already average during dN/dy plotting, but kaon was remaining
        gRatioPhiKaModel->SetLineWidth(3);
        gRatioPhiKaModel->Draw("l same");
        counter++;
    }

    legendRatio2->Draw();
    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.DrawLatex(0.68, 0.88, "ALICE");
    latex.DrawLatex(0.68, 0.83, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.68, 0.78, "|y| < 0.5");
    latex.DrawLatex(0.68, 0.73, "#Phi / K");

    if (isSavePlots)
    {
        cRatioPhiK->SaveAs("Plots/QA/Ratio_PhiKaon_Run3.png");
    }

    //==================================================
    // ===================Kstar / Pion Ratio ==================
    //==================================================
    TCanvas *cRatioKstarPion = new TCanvas("cRatioKstarPion", "cRatioKstarPion", 720, 720);
    SetCanvasStyle(cRatioKstarPion, 0.15, 0.03, 0.06, 0.15);
    TGraphErrors **gRatioKstarPion = MakeRatioUncorr(gMYieldKstar, gMYieldPion, 1.0, 3);
    gRatioKstarPion[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioKstarPion[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphErrorStyle(gRatioKstarPion[0]);
    // gRatioKstarPion[0]->GetYaxis()->SetRangeUser(0.029, 0.158);
    gRatioKstarPion[0]->GetYaxis()->SetRangeUser(0.031, 0.062);
    gRatioKstarPion[0]->GetXaxis()->SetLimits(1.8, 23.5);
    // gRatioKstarPion[0]->GetXaxis()->SetLimits(1.8, 26.5);
    gRatioKstarPion[0]->SetMarkerColor(kRed);
    gRatioKstarPion[0]->SetLineColor(kRed);
    gRatioKstarPion[0]->GetYaxis()->SetNdivisions(508);
    gRatioKstarPion[0]->GetYaxis()->SetMaxDigits(2);
    gRatioKstarPion[0]->Draw("APE");
    gRatioKstarPion[1]->SetFillStyle(0);
    gRatioKstarPion[1]->SetLineColor(kRed);
    gRatioKstarPion[1]->Draw("5 same");
    gRatioKstarPion[2]->SetFillStyle(1001);
    gRatioKstarPion[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioKstarPion[2]->Draw("E2 same");

    // gKstarPiRatio_13TeV[0]->SetMarkerStyle(21);
    // gKstarPiRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gKstarPiRatio_13TeV[0]->SetLineColor(kBlue);
    // gKstarPiRatio_13TeV[0]->Draw("P same");
    // gKstarPiRatio_13TeV[1]->SetLineColor(kBlue);
    // gKstarPiRatio_13TeV[1]->SetFillStyle(0);
    // gKstarPiRatio_13TeV[1]->Draw("5 same");

    gKstarPiRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    gKstarPiRatio_UrQMDOffSmooth->SetLineStyle(2);
    gKstarPiRatio_UrQMDOffSmooth->SetLineWidth(3);
    gKstarPiRatio_UrQMDOffSmooth->Draw("LX same");
    gKstarPiRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    gKstarPiRatio_UrQMDOnSmooth->SetLineStyle(9);
    gKstarPiRatio_UrQMDOnSmooth->SetLineWidth(3);
    gKstarPiRatio_UrQMDOnSmooth->Draw("LX same");

    counter = 0;
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarPiModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kPion], true);
        ScaleGraph(gRatioKstarPiModel, 2.0); // because Kstar yield was already average during dN/dy plotting, but pion was remaining
        setStyle(gRatioKstarPiModel, colorsPythia[counter], lineStylesPythia[counter]);
        gRatioKstarPiModel->SetLineWidth(3);
        gRatioKstarPiModel->Draw("l same");
        counter++;
    }

    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{#pi}");
    TLegend *legendRatio3 = new TLegend(0.65, 0.83, 0.75, 0.9);
    SetLegendStyle(legendRatio3);
    legendRatio3->SetTextSize(0.04);
    legendRatio3->Clear();
    legendRatio3->AddEntry(gRatioKstarPion[0], "K*^{0} / #pi", "p");
    legendRatio3->Draw();
    // legendRatio2->Draw();
    latex.SetTextSize(0.035);
    latex.SetTextFont(22);
    latex.DrawLatex(0.22, 0.88, "ALICE");
    latex.DrawLatex(0.22, 0.83, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.22, 0.78, "|y| < 0.5");

    if (isSavePlots)
    {
        cRatioKstarPion->SaveAs("Plots/QA/Ratio_KstarPion_Run3.png");
    }

    //==================================================
    // ===================Phi / Pion Ratio =============
    //==================================================
    TCanvas *cRatioPhiPion = new TCanvas("cRatioPhiPion", "cRatioPhiPion", 720, 720);
    SetCanvasStyle(cRatioPhiPion, 0.15, 0.03, 0.06, 0.15);
    TGraphErrors **gRatioPhiPion = MakeRatioUncorr(gMYieldPhiRun3, gMYieldPion, 1.0, 3);
    gRatioPhiPion[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioPhiPion[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphErrorStyle(gRatioPhiPion[0]);
    gRatioPhiPion[0]->GetYaxis()->SetRangeUser(7e-3, 28e-3);
    gRatioPhiPion[0]->GetXaxis()->SetLimits(1.8, 23.5);
    // gRatioPhiPion[0]->GetXaxis()->SetLimits(1.8, 26.5);
    gRatioPhiPion[0]->SetMarkerColor(kRed);
    gRatioPhiPion[0]->SetLineColor(kRed);
    gRatioPhiPion[0]->GetYaxis()->SetNdivisions(508);
    gRatioPhiPion[0]->GetYaxis()->SetMaxDigits(2);
    gRatioPhiPion[0]->Draw("APE");
    gRatioPhiPion[1]->SetFillStyle(0);
    gRatioPhiPion[1]->SetLineColor(kRed);
    gRatioPhiPion[1]->Draw("5 same");
    gRatioPhiPion[2]->SetFillStyle(1001);
    gRatioPhiPion[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioPhiPion[2]->Draw("E2 same");

    // ScaleGraph(gPhiPiRatio_13TeV[0], 2.0); // To average the pion yield
    // ScaleGraph(gPhiPiRatio_13TeV[1], 2.0);
    // gPhiPiRatio_13TeV[0]->SetMarkerStyle(21);
    // gPhiPiRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gPhiPiRatio_13TeV[0]->SetLineColor(kBlue);
    // gPhiPiRatio_13TeV[0]->Draw("P same");
    // gPhiPiRatio_13TeV[1]->SetLineColor(kBlue);
    // gPhiPiRatio_13TeV[1]->SetFillStyle(0);
    // gPhiPiRatio_13TeV[1]->Draw("5 same");

    gPhiPiRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    gPhiPiRatio_UrQMDOffSmooth->SetLineStyle(2);
    gPhiPiRatio_UrQMDOffSmooth->SetLineWidth(3);
    gPhiPiRatio_UrQMDOffSmooth->Draw("LX same");
    gPhiPiRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    gPhiPiRatio_UrQMDOnSmooth->SetLineStyle(9);
    gPhiPiRatio_UrQMDOnSmooth->SetLineWidth(3);
    gPhiPiRatio_UrQMDOnSmooth->Draw("LX same");

    counter = 0;
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiPiModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kPion], true);
        ScaleGraph(gRatioPhiPiModel, 2.0); // To average the pion yield
        // setStyle(gRatioPhiPiModel, modelStyle[model].color, modelStyle[model].style);
        setStyle(gRatioPhiPiModel, colorsPythia[counter], lineStylesPythia[counter]);
        gRatioPhiPiModel->SetLineWidth(3);
        gRatioPhiPiModel->Draw("l same");
        counter++;
    }

    legendRatio2->Draw();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.82, 0.85, "#Phi/#pi");

    cRatioPhiPion->SaveAs("Plots/QA/Ratio_PhiPion_Run3.png");

    //======================================================
    //================Kstar/Kshort Ratio====================
    //======================================================
    TCanvas *cRatioKstarKshort = new TCanvas("cRatioKstarKshort", "cRatioKstarKshort", 720, 720);
    SetCanvasStyle(cRatioKstarKshort, 0.15, 0.03, 0.06, 0.15);
    TGraphErrors **gRatioKstarKshort = MakeRatioUncorr(gMYieldKstar, gMYieldKshortRun3, 1.0, 3);
    gRatioKstarKshort[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioKstarKshort[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphErrorStyle(gRatioKstarKshort[0]);
    gRatioKstarKshort[0]->GetYaxis()->SetRangeUser(0.245, 0.54);
    gRatioKstarKshort[0]->GetXaxis()->SetLimits(1.8, 23.5);
    gRatioKstarKshort[0]->SetMarkerColor(kRed);
    gRatioKstarKshort[0]->SetLineColor(kRed);
    gRatioKstarKshort[0]->GetYaxis()->SetNdivisions(505);
    gRatioKstarKshort[0]->GetYaxis()->SetMaxDigits(3);
    gRatioKstarKshort[0]->Draw("APE");
    gRatioKstarKshort[1]->SetFillStyle(0);
    gRatioKstarKshort[1]->SetLineColor(kRed);
    gRatioKstarKshort[1]->Draw("5 same");
    gRatioKstarKshort[2]->SetFillStyle(1001);
    gRatioKstarKshort[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioKstarKshort[2]->Draw("E2 same");

    cout << "\n Information for K*0 and Kshort yields and their errors " << endl;
    for (int i3 = 0; i3 < TotalPointsRun3; i3++)
    {
        double xkstar, ykstar, ysyskstar, yuncorrkstar;
        double xkaon, ykaon, ysyskaon, yuncorrkaon;
        gMYieldKstar[0]->GetPoint(i3, xkstar, ykstar);
        ysyskstar = gMYieldKstar[1]->GetErrorY(i3);
        yuncorrkstar = gMYieldKstar[2]->GetErrorY(i3);

        gMYieldKshortRun3[0]->GetPoint(i3, xkaon, ykaon);
        ysyskaon = gMYieldKshortRun3[1]->GetErrorY(i3);
        yuncorrkaon = gMYieldKshortRun3[2]->GetErrorY(i3);

        cout << "Point " << i3 << ": K*0 yield = " << ykstar << " +/- " << ysyskstar << " (sys) +/- " << yuncorrkstar << " (uncorr)" << endl;
        cout << "Point " << i3 << ": Kshort yield = " << ykaon << " +/- " << ysyskaon << " (sys) +/- " << yuncorrkaon << " (uncorr)" << endl;
        double ratio = ykstar / ykaon;
        double ratioErrorsys = ratio * sqrt(pow(ysyskstar / ykstar, 2) + pow(ysyskaon / ykaon, 2));
        double ratioErroruncorr = ratio * sqrt(pow(yuncorrkstar / ykstar, 2) + pow(yuncorrkaon / ykaon, 2));
        cout << "Point " << i3 << ": K*0/Kshort ratio = " << ratio << " +/- " << ratioErrorsys << " (sys) +/- " << ratioErroruncorr << " (uncorr)" << endl;
        cout << endl;
    }

    ////Calculation of suppression significance between HM and LM for K*0/Kshort ratio
    int totalPointsKstarKshortRatio = gRatioKstarKshort[0]->GetN();
    double xLM_KstarKshort, xHM_KstarKshort, yHM_KstarKshort, yLM_KstarKshort, yHMError_KstarKshort, yLMError_KstarKshort, yHMStatError_KstarKshort, yLMStatError_KstarKshort, yLMUncorr_KstarKshort, yHMUncorr_KstarKshort;
    gRatioKstarKshort[1]->GetPoint(totalPointsKstarKshortRatio - 1, xLM_KstarKshort, yLM_KstarKshort);
    gRatioKstarKshort[1]->GetPoint(0, xHM_KstarKshort, yHM_KstarKshort);
    yLMError_KstarKshort = gRatioKstarKshort[1]->GetErrorY(totalPointsKstarKshortRatio - 1);
    yHMError_KstarKshort = gRatioKstarKshort[1]->GetErrorY(0);
    yLMStatError_KstarKshort = gRatioKstarKshort[0]->GetErrorY(totalPointsKstarKshortRatio - 1);
    yHMStatError_KstarKshort = gRatioKstarKshort[0]->GetErrorY(0);
    yLMUncorr_KstarKshort = gRatioKstarKshort[2]->GetErrorY(totalPointsKstarKshortRatio - 1);
    yHMUncorr_KstarKshort = gRatioKstarKshort[2]->GetErrorY(0);
    double sigmaDiffHM_LM_KstarKshort = abs(yHM_KstarKshort - yLM_KstarKshort) / sqrt(pow(yHMError_KstarKshort, 2) + pow(yLMError_KstarKshort, 2) + pow(yHMStatError_KstarKshort, 2) + pow(yLMStatError_KstarKshort, 2));
    double sigmaDiffHM_LM_uncorr_KstarKshort = abs(yHM_KstarKshort - yLM_KstarKshort) / sqrt(pow(yHMUncorr_KstarKshort, 2) + pow(yLMUncorr_KstarKshort, 2));
    cout << "K*0/Kshort ratio difference between HM and LM: " << sigmaDiffHM_LM_KstarKshort << " sigma" << endl;
    cout << "K*0/Kshort ratio difference between HM and LM (uncorrelated): " << sigmaDiffHM_LM_uncorr_KstarKshort << " sigma" << endl;

    gKstarK0sRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    gKstarK0sRatio_UrQMDOffSmooth->SetLineStyle(2);
    gKstarK0sRatio_UrQMDOffSmooth->SetLineWidth(3);
    gKstarK0sRatio_UrQMDOffSmooth->Draw("LX same");
    gKstarK0sRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    gKstarK0sRatio_UrQMDOnSmooth->SetLineStyle(9);
    gKstarK0sRatio_UrQMDOnSmooth->SetLineWidth(3);
    gKstarK0sRatio_UrQMDOnSmooth->Draw("LX same");

    counter = 0;
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarKshortModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKshort], true);
        setStyle(gRatioKstarKshortModel, colorsPythia[counter], lineStylesPythia[counter]);
        gRatioKstarKshortModel->SetLineWidth(3);
        gRatioKstarKshortModel->Draw("l same");
        counter++;
    }

    TLegend *legendRatio4 = new TLegend(0.75, 0.83, 0.9, 0.9);
    SetLegendStyle(legendRatio4);
    legendRatio4->SetTextSize(0.04);
    legendRatio4->AddEntry(gRatioKstarKshort[0], "K*^{0} / K_{S}^{0}", "p");
    legendRatio4->Draw();
    legendRatio2->Draw();

    if (isSavePlots)
    {
        cRatioKstarKshort->SaveAs("Plots/QA/Ratio_KstarKshort_Run3.png");
    }

    //=====================================================================
    // ===================Phi/K0s ratio ================================
    //=====================================================================
    TCanvas *cRatioPhiK0s = new TCanvas("cRatioPhiK0s", "cRatioPhiK0s", 720, 720);
    SetCanvasStyle(cRatioPhiK0s, 0.15, 0.03, 0.06, 0.15);
    TGraphErrors **gRatioPhiK0s = MakeRatioUncorr(gMYieldPhiRun3, gMYieldKshortRun3, 1.0, 3);
    gRatioPhiK0s[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioPhiK0s[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphErrorStyle(gRatioPhiK0s[0]);
    gRatioPhiK0s[0]->GetYaxis()->SetRangeUser(0.07, 0.18);
    gRatioPhiK0s[0]->GetXaxis()->SetLimits(1.8, 23.5);
    gRatioPhiK0s[0]->SetMarkerColor(kRed);
    gRatioPhiK0s[0]->SetLineColor(kRed);
    gRatioPhiK0s[0]->GetYaxis()->SetNdivisions(505);
    gRatioPhiK0s[0]->GetYaxis()->SetMaxDigits(3);
    gRatioPhiK0s[0]->Draw("APE");
    gRatioPhiK0s[1]->SetFillStyle(0);
    gRatioPhiK0s[1]->SetLineColor(kRed);
    gRatioPhiK0s[1]->Draw("5 same");
    gRatioPhiK0s[2]->SetFillStyle(1001);
    gRatioPhiK0s[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioPhiK0s[2]->Draw("E2 same");

    gPhiK0sRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    gPhiK0sRatio_UrQMDOffSmooth->SetLineStyle(2);
    gPhiK0sRatio_UrQMDOffSmooth->SetLineWidth(3);
    gPhiK0sRatio_UrQMDOffSmooth->Draw("LX same");
    gPhiK0sRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    gPhiK0sRatio_UrQMDOnSmooth->SetLineStyle(9);
    gPhiK0sRatio_UrQMDOnSmooth->SetLineWidth(3);
    gPhiK0sRatio_UrQMDOnSmooth->Draw("LX same");

    counter = 0;
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiKshortModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kKshort], true);
        // setStyle(gRatioPhiKshortModel, modelStyle[model].color, modelStyle[model].style);
        setStyle(gRatioPhiKshortModel, colorsPythia[counter], lineStylesPythia[counter]);
        gRatioPhiKshortModel->SetLineWidth(3);
        gRatioPhiKshortModel->Draw("l same");
        counter++;
    }

    TLegend *legendRatio5 = new TLegend(0.75, 0.83, 0.9, 0.9);
    SetLegendStyle(legendRatio5);
    legendRatio5->SetTextSize(0.04);
    legendRatio5->AddEntry(gRatioPhiK0s[0], "#Phi / K^{0}_{S}", "p");
    legendRatio5->Draw();
    legendRatio2->Draw();
    cRatioPhiK0s->SaveAs("Plots/QA/Ratio_PhiK0s_Run3.png");

    //======================================================================
    // ==========Double YieldRatio (Mult/LM): Kstar/K,Kstar/Pi ============
    //======================================================================
    TCanvas *cYieldLMRatio = new TCanvas("cYieldLMRatio", "cYieldLMRatio", 720, 720);
    SetCanvasStyle(cYieldLMRatio, 0.15, 0.03, 0.03, 0.15);

    TGraphErrors **gYieldKstarKaLMRatio = DivideByMult(gRatioKstarKaon, -1, 0.5, 3);
    TGraphErrors **gYieldKstarPiLMRatio = DivideByMult(gRatioKstarPion, -1, 0.5, 3);

    gYieldKstarKaLMRatio[0]->GetXaxis()->SetTickLength(0.02);
    gYieldKstarKaLMRatio[0]->GetYaxis()->SetTickLength(0.02);
    gYieldKstarKaLMRatio[0]->SetMaximum(1.33);
    gYieldKstarKaLMRatio[0]->SetMinimum(0.75);
    gYieldKstarKaLMRatio[0]->SetMarkerColor(kRed);
    gYieldKstarKaLMRatio[0]->SetLineColor(kRed);
    gYieldKstarKaLMRatio[0]->SetMarkerStyle(20);
    gYieldKstarKaLMRatio[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gYieldKstarKaLMRatio[0]->GetYaxis()->SetTitle("Y / Y_{LM}");
    gYieldKstarKaLMRatio[0]->SetFillStyle(0);
    gYieldKstarKaLMRatio[1]->SetFillStyle(0);
    gYieldKstarKaLMRatio[0]->Draw("APE");
    gYieldKstarKaLMRatio[1]->SetLineColor(kRed);
    gYieldKstarKaLMRatio[1]->Draw("5 same");

    gYieldKstarPiLMRatio[0]->SetMarkerColor(kGreen + 2);
    gYieldKstarPiLMRatio[0]->SetLineColor(kGreen + 2);
    gYieldKstarPiLMRatio[0]->SetMarkerStyle(21);
    gYieldKstarPiLMRatio[0]->SetFillStyle(0);
    gYieldKstarPiLMRatio[1]->SetFillStyle(0);
    gYieldKstarPiLMRatio[0]->Draw("PE same");
    gYieldKstarPiLMRatio[1]->SetLineColor(kGreen + 2);
    gYieldKstarPiLMRatio[1]->SetFillStyle(0);
    gYieldKstarPiLMRatio[1]->Draw("5 same");

    TGraphErrors *gEPOSNew_KstarKaLMRatio_UrQMDOff = DivideByMultModel(gKstarKaRatio_UrQMDOffSmooth, 2.5);
    TGraphErrors *gEPOSNew_KstarKaLMRatio_UrQMDOn = DivideByMultModel(gKstarKaRatio_UrQMDOnSmooth, 2.5);
    TGraphErrors *gEPOSNew_KstarPiLMRatio_UrQMDOff = DivideByMultModel(gKstarPiRatio_UrQMDOffSmooth, 2.5);
    TGraphErrors *gEPOSNew_KstarPiLMRatio_UrQMDOn = DivideByMultModel(gKstarPiRatio_UrQMDOnSmooth, 2.5);

    gEPOSNew_KstarKaLMRatio_UrQMDOff->SetLineColor(kGreen - 2);
    gEPOSNew_KstarKaLMRatio_UrQMDOff->SetLineStyle(9);
    gEPOSNew_KstarKaLMRatio_UrQMDOff->SetLineWidth(3);
    gEPOSNew_KstarKaLMRatio_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOSNew_KstarKaLMRatio_UrQMDOn->SetLineStyle(2);
    gEPOSNew_KstarKaLMRatio_UrQMDOn->SetLineWidth(3);
    gEPOSNew_KstarPiLMRatio_UrQMDOff->SetLineColor(kMagenta + 1);
    gEPOSNew_KstarPiLMRatio_UrQMDOff->SetLineStyle(9);
    gEPOSNew_KstarPiLMRatio_UrQMDOff->SetLineWidth(3);
    gEPOSNew_KstarPiLMRatio_UrQMDOn->SetLineColor(kRed + 2);
    gEPOSNew_KstarPiLMRatio_UrQMDOn->SetLineStyle(2);
    gEPOSNew_KstarPiLMRatio_UrQMDOn->SetLineWidth(3);

    gEPOSNew_KstarKaLMRatio_UrQMDOff->Draw("LX same");
    gEPOSNew_KstarKaLMRatio_UrQMDOn->Draw("LX same");
    gEPOSNew_KstarPiLMRatio_UrQMDOff->Draw("LX same");
    gEPOSNew_KstarPiLMRatio_UrQMDOn->Draw("LX same");

    TLegend *legendYieldLMRatio = new TLegend(0.2, 0.89, 0.95, 0.98);
    SetLegendStyle(legendYieldLMRatio);
    legendYieldLMRatio->SetTextSize(0.035);
    legendYieldLMRatio->SetNColumns(2);
    legendYieldLMRatio->AddEntry(gYieldKstarPiLMRatio[0], "K*^{0} / #pi", "P");
    legendYieldLMRatio->AddEntry(gYieldKstarKaLMRatio[0], "K*^{0} / K", "P");

    TLegend *legendYieldLMRatio2 = new TLegend(0.2, 0.78, 0.95, 0.92);
    SetLegendStyle(legendYieldLMRatio2);
    legendYieldLMRatio2->SetTextSize(0.028);
    legendYieldLMRatio2->SetNColumns(2);

    legendYieldLMRatio2->AddEntry(gEPOSNew_KstarPiLMRatio_UrQMDOff, "EPOS4", "l");
    legendYieldLMRatio2->AddEntry(gEPOSNew_KstarKaLMRatio_UrQMDOff, "EPOS4", "l");
    legendYieldLMRatio2->AddEntry(gEPOSNew_KstarPiLMRatio_UrQMDOn, "EPOS4 + UrQMD", "l");
    legendYieldLMRatio2->AddEntry(gEPOSNew_KstarKaLMRatio_UrQMDOn, "EPOS4 + UrQMD", "l");

    // counter = -1;
    // for (auto model : modelsToPlot)
    // {
    //     counter++;

    //     if (counter < 2) // plot only ropes
    //         continue;

    //     TGraphErrors *gYieldKstarKaPythiaModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKaon], true);
    //     ScaleGraph(gYieldKstarKaPythiaModel, 2.0); // because Kstar yield was already average during dN/dy plotting, but kaon was remaining
    //     setStyle(gYieldKstarKaPythiaModel, colorsPythia[counter], lineStylesPythia[counter]);

    //     TGraphErrors *gYieldKstarKaLMRatioPythia = (TGraphErrors *)gYieldKstarKaPythiaModel->Clone(Form("gYieldKstarKaLMRatioPythia_%d", counter));
    //     gYieldKstarKaLMRatioPythia = DivideByMultModel(gYieldKstarKaPythiaModel, 3.69);
    //     gYieldKstarKaLMRatioPythia->SetLineWidth(3);
    //     gYieldKstarKaLMRatioPythia->SetLineColor(colorsPythia[counter] + 1);
    //     gYieldKstarKaLMRatioPythia->SetLineStyle(1);
    //     gYieldKstarKaLMRatioPythia->Draw("l same");

    //     TGraphErrors *gYieldKstarPiPythiaModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kPion], true);
    //     ScaleGraph(gYieldKstarPiPythiaModel, 2.0); // because Kstar yield was already average during dN/dy plotting, but pion was remaining
    //     setStyle(gYieldKstarPiPythiaModel, colorsPythia[counter], lineStylesPythia[counter]);

    //     TGraphErrors *gYieldKstarPiLMRatioPythia = (TGraphErrors *)gYieldKstarPiPythiaModel->Clone(Form("gYieldKstarPiLMRatioPythia_%d", counter));
    //     gYieldKstarPiLMRatioPythia = DivideByMultModel(gYieldKstarPiPythiaModel, 3.69);
    //     gYieldKstarPiLMRatioPythia->SetLineWidth(3);
    //     gYieldKstarPiLMRatioPythia->SetLineColor(colorsPythia[counter] + 1);
    //     gYieldKstarPiLMRatioPythia->SetLineStyle(7);
    //     gYieldKstarPiLMRatioPythia->Draw("l same");

    //     legendYieldLMRatio2->AddEntry(gYieldKstarKaLMRatioPythia, Form("%s ", modelLabelLocal[counter]), "l");
    //     legendYieldLMRatio2->AddEntry(gYieldKstarPiLMRatioPythia, Form("%s ", modelLabelLocal[counter]), "l");
    // }

    legendYieldLMRatio2->Draw();
    legendYieldLMRatio->Draw();

    if (isSavePlots)
    {
        cYieldLMRatio->SaveAs("Plots/QA/Yield_LMRatio_KstarKaPi_Run3.png");
    }

    //======================================================================
    // ==========Double YieldRatio (Mult/LM): Phi/K, Phi/Pi ==========
    //======================================================================
    TCanvas *cYieldLMRatioPhi = new TCanvas("cYieldLMRatioPhi", "cYieldLMRatioPhi", 720, 720);
    SetCanvasStyle(cYieldLMRatioPhi, 0.15, 0.03, 0.03, 0.15);

    TGraphErrors **gYieldPhiKaLMRatio = DivideByMult(gRatioPhiK, -1, 0.5, 3);
    TGraphErrors **gYieldPhiPiLMRatio = DivideByMult(gRatioPhiPion, -1, 0.5, 3);

    gYieldPhiKaLMRatio[0]->GetXaxis()->SetTickLength(0.02);
    gYieldPhiKaLMRatio[0]->GetYaxis()->SetTickLength(0.02);
    gYieldPhiKaLMRatio[0]->SetMaximum(1.60);
    gYieldPhiKaLMRatio[0]->SetMinimum(0.95);
    gYieldPhiKaLMRatio[0]->SetMarkerColor(kRed);
    gYieldPhiKaLMRatio[0]->SetLineColor(kRed);
    gYieldPhiKaLMRatio[0]->SetMarkerStyle(20);
    gYieldPhiKaLMRatio[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gYieldPhiKaLMRatio[0]->GetYaxis()->SetTitle("Y / Y_{LM}");
    gYieldPhiKaLMRatio[0]->SetFillStyle(0);
    gYieldPhiKaLMRatio[1]->SetFillStyle(0);
    gYieldPhiKaLMRatio[0]->Draw("APE");
    gYieldPhiKaLMRatio[1]->SetLineColor(kRed);
    gYieldPhiKaLMRatio[1]->Draw("5 same");

    gYieldPhiPiLMRatio[0]->SetMarkerColor(kGreen + 2);
    gYieldPhiPiLMRatio[0]->SetLineColor(kGreen + 2);
    gYieldPhiPiLMRatio[0]->SetMarkerStyle(21);
    gYieldPhiPiLMRatio[0]->SetFillStyle(0);
    gYieldPhiPiLMRatio[1]->SetFillStyle(0);
    gYieldPhiPiLMRatio[0]->Draw("PE same");
    gYieldPhiPiLMRatio[1]->SetLineColor(kGreen + 2);
    gYieldPhiPiLMRatio[1]->SetFillStyle(0);
    gYieldPhiPiLMRatio[1]->Draw("5 same");

    TGraphErrors *gEPOS_PhiKaLMRatio_UrQMDOff = DivideByMultModel(gPhiKaRatio_UrQMDOffSmooth, 2.5);
    TGraphErrors *gEPOS_PhiKaLMRatio_UrQMDOn = DivideByMultModel(gPhiKaRatio_UrQMDOnSmooth, 2.5);
    TGraphErrors *gEPOS_PhiPiLMRatio_UrQMDOff = DivideByMultModel(gPhiPiRatio_UrQMDOffSmooth, 2.5);
    TGraphErrors *gEPOS_PhiPiLMRatio_UrQMDOn = DivideByMultModel(gPhiPiRatio_UrQMDOnSmooth, 2.5);

    gEPOS_PhiKaLMRatio_UrQMDOff->SetLineColor(kGreen - 2);
    gEPOS_PhiKaLMRatio_UrQMDOff->SetLineStyle(9);
    gEPOS_PhiKaLMRatio_UrQMDOff->SetLineWidth(3);
    gEPOS_PhiKaLMRatio_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_PhiKaLMRatio_UrQMDOn->SetLineStyle(2);
    gEPOS_PhiKaLMRatio_UrQMDOn->SetLineWidth(3);
    gEPOS_PhiPiLMRatio_UrQMDOff->SetLineColor(kMagenta + 1);
    gEPOS_PhiPiLMRatio_UrQMDOff->SetLineStyle(9);
    gEPOS_PhiPiLMRatio_UrQMDOff->SetLineWidth(3);
    gEPOS_PhiPiLMRatio_UrQMDOn->SetLineColor(kRed + 2);
    gEPOS_PhiPiLMRatio_UrQMDOn->SetLineStyle(2);
    gEPOS_PhiPiLMRatio_UrQMDOn->SetLineWidth(3);
    gEPOS_PhiKaLMRatio_UrQMDOff->Draw("LX same");
    gEPOS_PhiKaLMRatio_UrQMDOn->Draw("LX same");
    gEPOS_PhiPiLMRatio_UrQMDOff->Draw("LX same");
    gEPOS_PhiPiLMRatio_UrQMDOn->Draw("LX same");

    TLegend *legendYieldLMRatioPhi = new TLegend(0.2, 0.78, 0.95, 0.95);
    SetLegendStyle(legendYieldLMRatioPhi);
    legendYieldLMRatioPhi->SetTextSize(0.028);
    legendYieldLMRatioPhi->SetNColumns(2);
    legendYieldLMRatioPhi->AddEntry(gYieldPhiPiLMRatio[0], "#phi / #pi", "P");
    legendYieldLMRatioPhi->AddEntry(gYieldPhiKaLMRatio[0], "#phi / K", "P");
    legendYieldLMRatioPhi->AddEntry(gEPOS_PhiKaLMRatio_UrQMDOff, "#phi / K (UrQMD Off)", "L");
    legendYieldLMRatioPhi->AddEntry(gEPOS_PhiPiLMRatio_UrQMDOff, "#phi / #pi (UrQMD Off)", "L");
    legendYieldLMRatioPhi->AddEntry(gEPOS_PhiKaLMRatio_UrQMDOn, "#phi / K (UrQMD On)", "L");
    legendYieldLMRatioPhi->AddEntry(gEPOS_PhiPiLMRatio_UrQMDOn, "#phi / #pi (UrQMD On)", "L");

    // counter = -1;
    // for (auto model : modelsToPlot)
    // {
    //     counter++;

    //     if (counter < 2) // Plot only ropes
    //         continue;

    //     TGraphErrors *gYieldPhiKaPythiaModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kKaon], true);
    //     ScaleGraph(gYieldPhiKaPythiaModel, 2.0); // because Phi yield was already average during dN/dy plotting, but kaon was remaining
    //     setStyle(gYieldPhiKaPythiaModel, colorsPythia[counter], lineStylesPythia[counter]);

    //     TGraphErrors *gYieldPhiKaLMRatioPythia = (TGraphErrors *)gYieldPhiKaPythiaModel->Clone(Form("gYieldPhiKaLMRatioPythia_%d", counter));
    //     gYieldPhiKaLMRatioPythia = DivideByMultModel(gYieldPhiKaPythiaModel, 3.69);
    //     gYieldPhiKaLMRatioPythia->SetLineWidth(3);
    //     gYieldPhiKaLMRatioPythia->SetLineColor(colorsPythia[counter] + 1);
    //     gYieldPhiKaLMRatioPythia->SetLineStyle(1);
    //     gYieldPhiKaLMRatioPythia->Draw("l same");

    //     TGraphErrors *gYieldPhiPiPythiaModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kPion], true);
    //     ScaleGraph(gYieldPhiPiPythiaModel, 2.0); // because Phi yield was already average during dN/dy plotting, but pion was remaining
    //     setStyle(gYieldPhiPiPythiaModel, colorsPythia[counter], lineStylesPythia[counter]);

    //     TGraphErrors *gYieldPhiPiLMRatioPythia = (TGraphErrors *)gYieldPhiPiPythiaModel->Clone(Form("gYieldPhiPiLMRatioPythia_%d", counter));
    //     gYieldPhiPiLMRatioPythia = DivideByMultModel(gYieldPhiPiPythiaModel, 3.69);
    //     gYieldPhiPiLMRatioPythia->SetLineWidth(3);
    //     gYieldPhiPiLMRatioPythia->SetLineColor(colorsPythia[counter] + 1);
    //     gYieldPhiPiLMRatioPythia->SetLineStyle(7);
    //     gYieldPhiPiLMRatioPythia->Draw("l same");

    //     legendYieldLMRatioPhi->AddEntry(gYieldPhiKaLMRatioPythia, Form("%s ", modelLabelLocal[counter]), "l");
    //     legendYieldLMRatioPhi->AddEntry(gYieldPhiPiLMRatioPythia, Form("%s ", modelLabelLocal[counter]), "l");
    // }

    legendYieldLMRatioPhi->Draw("same");

    TLine *lineDoubleRatioPhi = new TLine(gYieldPhiKaLMRatio[0]->GetXaxis()->GetXmin(), 1.0, gYieldPhiKaLMRatio[0]->GetXaxis()->GetXmax(), 1.0);
    lineDoubleRatioPhi->SetLineStyle(2);
    lineDoubleRatioPhi->SetLineWidth(2);
    if (isSavePlots)
    {
        cYieldLMRatioPhi->SaveAs("Plots/QA/Yield_LMRatio_PhiKaPi_Run3.png");
    }
}

TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false, float CorrFactorDen = 1)
{
    if (numerator == nullptr || denominator == nullptr)
    {
        cout << "Error: null input graph while building ratio" << endl;
        return nullptr;
    }

    if (numerator->GetN() != denominator->GetN())
    {
        cout << "Error: ratio inputs have different number of points" << endl;
        return nullptr;
    }

    TGraphErrors *ratio = new TGraphErrors(numerator->GetN());
    for (int i = 0; i < numerator->GetN(); ++i)
    {
        double xNumerator = 0;
        double yNumerator = 0;
        double xDenominator = 0;
        double yDenominator = 0;
        numerator->GetPoint(i, xNumerator, yNumerator);
        denominator->GetPoint(i, xDenominator, yDenominator);

        // if (std::round(xNumerator * 10.0) != std::round(xDenominator * 10.0))
        // {
        //     cout << "Warning!!!!! x values do not match for ratio at point "
        //          << i << ": " << xNumerator << " vs " << xDenominator
        //          << ", for graphs " << numerator->GetName()
        //          << " and " << denominator->GetName() << endl;
        // }

        yDenominator *= CorrFactorDen; // Apply correction factor to denominator if needed

        double yRatio = (yDenominator != 0) ? yNumerator / yDenominator : 0;
        ratio->SetPoint(i, xNumerator, yRatio);

        double xErrorNumerator = numerator->GetErrorX(i);
        double numeratorError = numerator->GetErrorY(i);
        double denominatorError = denominator->GetErrorY(i);
        denominatorError *= CorrFactorDen; // Apply correction factor to denominator error if needed
        double yRatioError = 0;
        if (yDenominator != 0)
        {
            yRatioError = sqrt(pow(numeratorError / yDenominator, 2) + pow(yNumerator * denominatorError / (yDenominator * yDenominator), 2));
        }
        if (isModel)
            ratio->SetPointError(i, 0, 0);
        else
            ratio->SetPointError(i, xErrorNumerator, yRatioError);
    }

    SetGraphErrorStyle(ratio);

    return ratio;
}

TGraphErrors **MakeRatioUncorr(TGraphErrors **numerator, TGraphErrors **denominator, float CorrFactorDen = 1, int nGraphs = 3)
{
    if (!numerator || !denominator)
    {
        cout << "Error: null input graph array" << endl;
        return nullptr;
    }

    if (nGraphs < 3)
    {
        cout << "Error: nGraphs must be at least 3" << endl;
        return nullptr;
    }

    TGraphErrors **ratio = new TGraphErrors *[3];

    for (int ig = 0; ig < 3; ig++)
    {
        if (!numerator[ig] || !denominator[ig])
        {
            cout << "Error: null graph at index " << ig << endl;
            delete[] ratio;
            return nullptr;
        }

        if (numerator[ig]->GetN() != denominator[ig]->GetN())
        {
            cout << "Error: different number of points in numerator " << "and denominator graphs" << endl;

            cout << "ig = " << ig << ", numerator = " << numerator[ig]->GetN() << ", denominator = " << denominator[ig]->GetN() << endl;

            delete[] ratio;
            return nullptr;
        }

        ratio[ig] = new TGraphErrors(numerator[ig]->GetN());

        for (int i = 0; i < numerator[ig]->GetN(); i++)
        {
            double xNum, yNum;
            double xDen, yDen;

            numerator[ig]->GetPoint(i, xNum, yNum);
            denominator[ig]->GetPoint(i, xDen, yDen);

            // Apply correction factor to denominator
            yDen *= CorrFactorDen;

            double yRatio = (yDen != 0.) ? yNum / yDen : 0.;

            ratio[ig]->SetPoint(i, xNum, yRatio);

            // ---------------------------------------------------------
            // Select the errors according to the requested definition
            // ---------------------------------------------------------

            double numErr = 0.;
            double denErr = 0.;

            if (ig == 0)
            {
                // Ratio 0:
                // stat(num) / stat(den)
                numErr = numerator[0]->GetErrorY(i);
                denErr = denominator[0]->GetErrorY(i);
            }
            else if (ig == 1)
            {
                // Ratio 1:
                // uncorr(num) / total(den) (wrong)
                // total (num) / total (den) (correct)
                numErr = numerator[1]->GetErrorY(i);
                denErr = denominator[1]->GetErrorY(i);
            }
            else if (ig == 2)
            {
                // Ratio 2:
                // uncorr(num) / uncorr(den)
                numErr = numerator[2]->GetErrorY(i);
                denErr = denominator[2]->GetErrorY(i);
            }

            // Correction factor also scales the denominator uncertainty
            denErr *= CorrFactorDen;

            double yRatioErr = 0.;

            if (yDen != 0.)
            {
                yRatioErr = sqrt(
                    pow(numErr / yDen, 2) +
                    pow(yNum * denErr / (yDen * yDen), 2));
            }

            double xErr = numerator[ig]->GetErrorX(i);

            ratio[ig]->SetPointError(i, xErr, yRatioErr);
        }

        SetGraphErrorStyle(ratio[ig]);
    }

    return ratio;
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

void ScaleGraph(TGraph *gr, double scale)
{
    if (!gr)
        return;

    for (int i = 0; i < gr->GetN(); ++i)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        gr->SetPoint(i, x, y * scale);
    }

    if (auto *ge = dynamic_cast<TGraphErrors *>(gr))
    {
        for (int i = 0; i < ge->GetN(); ++i)
        {
            ge->SetPointError(i, ge->GetErrorX(i), ge->GetErrorY(i) * scale);
        }
    }
}

TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance = 0.5, int noOfGraphs = 3)
{
    if (!gr || !gr[0] || !gr[2])
    {
        cout << "Error: null graph array or required graphs" << endl;
        return nullptr;
    }

    if (noOfGraphs < 3)
    {
        cout << "Error: noOfGraphs must be at least 3" << endl;
        return nullptr;
    }

    // Output:
    // grCopy[0] = ratio with statistical uncertainty
    // grCopy[1] = ratio with uncorrelated uncertainty
    TGraphErrors **grCopy = new TGraphErrors *[2];

    // ------------------------------------------------------------
    // Find the reference multiplicity point
    // ------------------------------------------------------------

    int refIndex = -1;

    if (WhichMultPoint >= 0)
    {
        for (int i = 0; i < gr[0]->GetN(); i++)
        {
            double x, y;
            gr[0]->GetPoint(i, x, y);

            if (fabs(x - WhichMultPoint) < tolerance)
            {
                refIndex = i;
                break;
            }
        }

        if (refIndex < 0)
        {
            cout << "Error: Could not find multiplicity point "
                 << WhichMultPoint << endl;

            delete[] grCopy;
            return nullptr;
        }
    }
    else
    {
        // Use the lowest multiplicity point
        double xFirst, yFirst;
        double xLast, yLast;

        gr[0]->GetPoint(0, xFirst, yFirst);
        gr[0]->GetPoint(gr[0]->GetN() - 1, xLast, yLast);

        refIndex = (xFirst < xLast) ? 0 : gr[0]->GetN() - 1;
    }

    // ------------------------------------------------------------
    // Reference value
    // ------------------------------------------------------------

    double xRef, yRef;
    gr[0]->GetPoint(refIndex, xRef, yRef);

    if (yRef == 0.)
    {
        cout << "Error: Reference value is zero" << endl;
        delete[] grCopy;
        return nullptr;
    }

    // Reference uncertainties
    double yRefStatErr = gr[0]->GetErrorY(refIndex);
    double yRefUncorrErr = gr[2]->GetErrorY(refIndex);

    // ------------------------------------------------------------
    // Create the two output graphs
    // ------------------------------------------------------------

    grCopy[0] = (TGraphErrors *)gr[0]->Clone();
    grCopy[1] = (TGraphErrors *)gr[2]->Clone();

    // ------------------------------------------------------------
    // Calculate ratios
    // ------------------------------------------------------------

    for (int i = 0; i < gr[0]->GetN(); i++)
    {
        double x, y;
        gr[0]->GetPoint(i, x, y);

        double xerr = gr[0]->GetErrorX(i);

        // ========================================================
        // 1. STATISTICAL RATIO
        // ========================================================

        double yStatErr = gr[0]->GetErrorY(i);

        double ratio = y / yRef;

        double ratioStatErr = sqrt(
            pow(yStatErr / yRef, 2) +
            pow(y * yRefStatErr / (yRef * yRef), 2));

        grCopy[0]->SetPoint(i, x, ratio);
        grCopy[0]->SetPointError(i, xerr, ratioStatErr);

        // ========================================================
        // 2. UNCORRELATED SYSTEMATIC RATIO
        // ========================================================

        double yUncorrErr = gr[2]->GetErrorY(i);

        double ratioUncorrErr = sqrt(pow(yUncorrErr / yRef, 2) + pow(y * yRefUncorrErr / (yRef * yRef), 2));

        grCopy[1]->SetPoint(i, x, ratio);
        grCopy[1]->SetPointError(i, xerr, ratioUncorrErr);
    }

    SetGraphErrorStyle(grCopy[0]);
    SetGraphErrorStyle(grCopy[1]);

    return grCopy;
}

void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax)
{
    for (int i = gr->GetN() - 1; i >= 0; --i)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        if (x < xMin || x > xMax)
        {
            gr->RemovePoint(i);
        }
    }
}

TGraphErrors *DivideByMultModel(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5)
{
    TGraphErrors *grCopy = (TGraphErrors *)gr->Clone();

    double yGivenMult = 1.0;
    double yGivenMultErr = 0.0;

    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);

        if (fabs(x - WhichMultPoint) < tolerance)
        {
            yGivenMult = y;
            yGivenMultErr = gr->GetErrorY(i);
            break;
        }
    }

    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);

        double yerr = gr->GetErrorY(i);
        double xerr = gr->GetErrorX(i);

        // Normalization point
        if (fabs(x - WhichMultPoint) < tolerance)
        {
            grCopy->SetPoint(i, x, 1.0);
            grCopy->SetPointError(i, xerr, 0.0);
            continue;
        }

        double ratio = y / yGivenMult;

        double ratioErr = sqrt(pow(yerr / yGivenMult, 2) + pow(yGivenMultErr * y / (yGivenMult * yGivenMult), 2));

        grCopy->SetPoint(i, x, ratio);
        grCopy->SetPointError(i, xerr, ratioErr);
    }

    return grCopy;
}

TGraphErrors *SmoothGraph(const TGraphErrors *gr)
{
    int n = gr->GetN();

    TGraphErrors *grSmooth = new TGraphErrors(n);

    for (int i = 0; i < n; i++)
    {

        double x, y;
        gr->GetPoint(i, x, y);

        double ySmooth = 0.0;
        double errSmooth = 0.0;
        int nPoints = 0;

        // Average neighbouring points
        for (int j = i - 1; j <= i + 1; j++)
        {

            if (j < 0 || j >= n)
                continue;

            double xj, yj;
            gr->GetPoint(j, xj, yj);

            double eyj = gr->GetErrorY(j);

            ySmooth += yj;
            errSmooth += eyj * eyj;

            nPoints++;
        }

        ySmooth /= nPoints;

        // Error propagation for the average
        errSmooth = std::sqrt(errSmooth) / nPoints;

        grSmooth->SetPoint(i, x, ySmooth);
        grSmooth->SetPointError(
            i,
            gr->GetErrorX(i),
            errSmooth);
    }

    return grSmooth;
}

TGraphErrors *MakeDataModelRatio(TGraphErrors *data, TGraphErrors *model)
{
    if (!data || !model)
        return nullptr;

    TGraphErrors *ratio = new TGraphErrors();

    for (int i = 0; i < data->GetN(); i++)
    {
        double x, yData;
        data->GetPoint(i, x, yData);

        double yModel = model->Eval(x);

        if (yModel <= 0)
            continue;

        double ratioValue = yData / yModel;
        double ratioError = data->GetErrorY(i) / yModel;

        int n = ratio->GetN();
        ratio->SetPoint(n, x, ratioValue);
        ratio->SetPointError(n, data->GetErrorX(i), ratioError);
    }

    return ratio;
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.08, 0.15);
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

void CalculateKstarPhiRatioRun2(TGraphErrors *kstar[3], TGraphErrors *phi[3], TGraphErrors *ratio[3])
{
    // Mapping of phi multiplicity bins to Kstar multiplicity bins
    // Phi:   0-1, 1-5, 5-10, 10-15, 15-20, 20-30, 30-40, 40-50, 50-70, 70-100
    // Kstar: 0-1, 1-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-70, 70-100

    const int nPhi = phi[0]->GetN();

    int kstarBin[11] = {0, 1, 2, 3, 3, 4, 5, 6, 7, 8, 9};

    const char *graphNames[3] = {"Ratio_Stat", "Ratio_Sys", "Ratio_UncorrSys"};

    for (int j = 0; j < 3; j++)
    {
        ratio[j] = new TGraphErrors(nPhi);
        ratio[j]->SetName(graphNames[j]);
        ratio[j]->SetTitle(graphNames[j]);

        for (int i = 0; i < nPhi; i++)
        {

            double xK, yK, xP, yP;

            kstar[j]->GetPoint(kstarBin[i], xK, yK);
            phi[j]->GetPoint(i, xP, yP);

            double errK = kstar[j]->GetErrorY(kstarBin[i]);
            double errP = phi[j]->GetErrorY(i);

            double ex = phi[j]->GetErrorX(i);

            if (yP == 0 || yK == 0)
            {
                ratio[j]->SetPoint(i, xP, 0);
                ratio[j]->SetPointError(i, ex, 0);
                continue;
            }

            double r = yK / yP;

            double errR = r * TMath::Sqrt(TMath::Power(errK / yK, 2) + TMath::Power(errP / yP, 2));

            ratio[j]->SetPoint(i, xP, r);
            ratio[j]->SetPointError(i, ex, errR);
        }
    }
}
