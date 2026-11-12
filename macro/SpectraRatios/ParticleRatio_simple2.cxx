#include <iostream>
#include <iomanip>
#include "../src/style.h"
#include "../src/StyleProposal.h"
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

void setStyle(TGraphErrors *gr, int color, int style)
{
    gr->SetLineColor(color);
    gr->SetLineStyle(style);
}

void ScaleGraph(TGraph *gr, double scale);
TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance = 0.5, int noOfGraphs = 3);
TGraphErrors *DivideByMultModel(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5);
void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax);

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
        double exNorm = sqrt(pow(ex / multMB, 2) + pow(x * multMBErr / (multMB * multMB), 2));
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
        double exNorm = sqrt(pow(ex / multMB, 2) + pow(x * multMBErr / (multMB * multMB), 2));
        double eyLow = sqrt(pow(ey / yieldMB, 2) + pow(y * grMB->GetErrorYlow(0) / (yieldMB * yieldMB), 2));
        double eyHigh = sqrt(pow(ey / yieldMB, 2) + pow(y * grMB->GetErrorYhigh(0) / (yieldMB * yieldMB), 2));

        int ip = grNorm->GetN();
        grNorm->SetPoint(ip, xNorm, yNorm);
        grNorm->SetPointError(ip, exNorm, exNorm, eyLow, eyHigh);
        // grNorm->SetPointError(ip, ex, ex, eyLow, eyHigh);
    }

    return grNorm;
}

void ParticleRatio_simple2()
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
    TFile *fXiRun3 = OpenFile("XiOmegaResults/Sawan/ResultsXi.root");

    // ALICE Run2 results
    TFile *fpp13TeV = OpenFile("ConversionCodes/pp13TeV_KstarPhi.root");
    TFile *fpp7TeV = OpenFile("ConversionCodes/pp7TeVALICE.root");

    TGraphErrors *gMPtKstar[3], *gMYieldKstar[3], *gMPtPion[3], *gMPtProton[3], *gMPtKaon[3], *gMPtKstarRun2[3], *gMPtPhiRun2[3], *gMPtChKstarRun2[3], *gMPtKshortRun2[3], *gMPtLambdaRun2[3], *gMPtXiRun2[3], *gMPtOmegaRun2[3], *gMPtXiStarRun2[3], *gMPtSigmaRun2[3], *gMPtLambda1520[3], *gMPtRho[3], *gMPtPhiRun3[3];
    TGraphErrors *gMYieldPion[3], *gMYieldProton[3], *gMYieldKaon[3], *gMYieldKstarRun2[3], *gMYieldPhiRun2[3], *gMYieldChKstarRun2[3], *gMYieldKshortRun2[3], *gMYieldLambdaRun2[3], *gMYieldXiRun2[3], *gMYieldOmegaRun2[3], *gMYieldXiStarRun2[3], *gMYieldSigmaRun2[3], *gMYieldLambda1520[3], *gMYieldRho[3], *gKstar_MeanpT_13TeV[3], *gKstar_Yield_13TeV[3], *gMYieldPhiRun3[3], *gMYieldXiRun3[3];

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

        gMYieldXiRun3[i] = GetGraph(fXiRun3, Form("gMeanYieldRun3%s", suffix.c_str()));
    }

    TGraphErrors *gKstarKshortRatio_13TeV[2], *gPhiKshortRatio_13TeV[2], *gPhiPiRatio_13TeV[2], *gKstarKaRatio_13TeV[2], *gKstarPiRatio_13TeV[2], *gPhiKaRatio_13TeV[2], *gPhiPrRatio_13TeV[2], *gPhiLambdaRatio_13TeV[2];
    TGraphErrors *gKstarPiRatio_7TeV[2], *gPhiPiRatio_7TeV[2];

    TGraphErrors *gMeanYieldMB_stat = GetGraph(fKstar, "gMeanYieldMB_stat");
    TGraphAsymmErrors *gMeanYieldMB_sys = (TGraphAsymmErrors *)fKstar->Get("gMeanYieldMB_sys");
    if (!gMeanYieldMB_stat || !gMeanYieldMB_sys || gMeanYieldMB_stat->GetN() < 1 || gMeanYieldMB_sys->GetN() < 1)
    {
        cout << "Error: minimum-bias yield graphs not found in " << fKstar->GetName() << endl;
        return;
    }

    TFile *fPhiMinBias = OpenFile("PhiResults/PhiAnalysisResults_pprun3.root");
    TGraphErrors *gMeanYieldMB_Phi_stat = GetGraph(fPhiMinBias, "YieldMean/dNdy/gdNdyStatMB");
    TGraphErrors *gMeanYieldMB_Phi_sys = GetGraph(fPhiMinBias, "YieldMean/dNdy/gdNdySysTotalMB");
    if (!gMeanYieldMB_Phi_stat || !gMeanYieldMB_Phi_sys)
    {
        cout << "Error: minimum-bias yield graphs for Phi not found in " << fPhiMinBias->GetName() << endl;
        return;
    }

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
    // RestrictModelXaxis(gEPOS_KstardNdy_UrQMDOff, 3.3, 23.5);
    // RestrictModelXaxis(gEPOS_KstardNdy_UrQMDOn, 3.3, 24.5);
    // RestrictModelXaxis(gEPOS_KstarMeanpT_UrQMDOff, 3.3, 23.5);
    // RestrictModelXaxis(gEPOS_KstarMeanpT_UrQMDOn, 3.3, 24.5);

    TGraphErrors *gEPOS_PhidNdy_UrQMDOff = GetGraph(fEPOS_New, "g_phi_dNdy_opt0");
    TGraphErrors *gEPOS_PhidNdy_UrQMDOn = GetGraph(fEPOS_New, "g_phi_dNdy_opt1");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOff = GetGraph(fEPOS_New, "g_phi_pT_opt0");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOn = GetGraph(fEPOS_New, "g_phi_pT_opt1");
    // RestrictModelXaxis(gEPOS_PhidNdy_UrQMDOff, 3.3, 23.5);
    // RestrictModelXaxis(gEPOS_PhidNdy_UrQMDOn, 3.3, 24.5);
    // RestrictModelXaxis(gEPOS_PhiMeanpT_UrQMDOff, 3.3, 23.5);
    // RestrictModelXaxis(gEPOS_PhiMeanpT_UrQMDOn, 3.3, 24.5);

    for (int i = gEPOS_PhidNdy_UrQMDOff->GetN() - 1; i >= 0; --i)
    {
        gEPOS_PhidNdy_UrQMDOff->SetPointError(i, 0, 0);   // Remove error
        gEPOS_PhidNdy_UrQMDOn->SetPointError(i, 0, 0);    // Remove error
        gEPOS_PhiMeanpT_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_PhiMeanpT_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
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

    // ========================================================================
    // CHANGE 1: Side-by-side dN/dy plot of K* and Phi with zero gap
    // ========================================================================
    TCanvas *cdNdy_SideBySide = new TCanvas("cdNdy_SideBySide", "dN/dy K* and Phi", 1000, 600);

    // Left Pad: K* (892)0
    TPad *padLeft = new TPad("padLeft", "padLeft", 0.0, 0.0, 0.5, 1.0);
    padLeft->SetLeftMargin(0.18);
    padLeft->SetRightMargin(0.00); // Zero gap
    padLeft->SetTopMargin(0.06);
    padLeft->SetBottomMargin(0.15);
    padLeft->SetTickx();
    // padLeft->SetTicky();
    padLeft->Draw();

    // Right Pad: Phi (1020)
    TPad *padRight = new TPad("padRight", "padRight", 0.5, 0.0, 1.0, 1.0);
    padRight->SetLeftMargin(0.00); // Zero gap
    padRight->SetRightMargin(0.18);
    padRight->SetTopMargin(0.06);
    padRight->SetBottomMargin(0.15);
    padRight->SetTickx();
    // padRight->SetTicky();
    padRight->Draw();

    // --- Drawing Left Pad (K*) ---
    padLeft->cd();
    gMYieldKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMYieldKstar[0]->GetYaxis()->SetTitle("dN/dy");
    gMYieldKstar[0]->GetXaxis()->SetLimits(0, 29);
    gMYieldKstar[0]->GetYaxis()->SetRangeUser(0, 0.81);
    gMYieldKstar[0]->SetLineWidth(3);
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

    for (auto model : modelsToPlot)
    {
        ScaleGraph(gMYield[model][kKstar], 0.5);
        gMYield[model][kKstar]->Draw("l same");
    }

    int totalPointsEPOS = gEPOS_KstardNdy_UrQMDOff->GetN();
    for (int i = 0; i < totalPointsEPOS; i++)
    {
        gEPOS_KstardNdy_UrQMDOff->SetPointError(i, 0, 0);
        gEPOS_KstardNdy_UrQMDOn->SetPointError(i, 0, 0);
    }
    gEPOS_KstardNdy_UrQMDOff->SetLineColor(kGreen + 2);
    gEPOS_KstardNdy_UrQMDOff->SetLineStyle(2);
    gEPOS_KstardNdy_UrQMDOff->SetLineWidth(3);
    gEPOS_KstardNdy_UrQMDOff->Draw("l same");

    gEPOS_KstardNdy_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_KstardNdy_UrQMDOn->SetLineStyle(9);
    gEPOS_KstardNdy_UrQMDOn->SetLineWidth(3);
    gEPOS_KstardNdy_UrQMDOn->Draw("l same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.24, 0.88, "ALICE");
    latex.DrawLatex(0.24, 0.82, "|y| < 0.5");
    latex.SetTextFont(22);
    latex.DrawLatex(0.24, 0.76, "K* (892)^{0}");

    TLegend *legendData = new TLegend(0.5, 0.8, 0.70, 0.92);
    legendData->SetBorderSize(0);
    legendData->SetFillStyle(0);
    legendData->SetTextFont(42);
    legendData->SetTextSize(0.035);
    legendData->AddEntry(gMYieldKstar[0], "pp #sqrt{s} = 13.6 TeV", "p");
    legendData->AddEntry(gKstar_Yield_13TeV[0], "pp #sqrt{s} = 13 TeV", "p");
    legendData->Draw("same");

    // --- Drawing Right Pad (Phi) ---
    padRight->cd();
    gMYieldPhiRun3[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMYieldPhiRun3[0]->GetYaxis()->SetTitle("dN/dy");
    gMYieldPhiRun3[0]->GetYaxis()->SetTitleOffset(1.2);
    gMYieldPhiRun3[0]->GetYaxis()->SetAxisColor(kBlack);
    gMYieldPhiRun3[0]->GetYaxis()->SetLabelColor(kBlack);
    gMYieldPhiRun3[0]->GetXaxis()->SetLimits(-1, 29);
    gMYieldPhiRun3[0]->GetYaxis()->SetRangeUser(0, 0.28);
    gMYieldPhiRun3[0]->GetYaxis()->SetTitleOffset(1.6);
    gMYieldPhiRun3[0]->SetLineWidth(3);
    gMYieldPhiRun3[0]->SetMarkerColor(kRed + 1);
    gMYieldPhiRun3[0]->SetLineColor(kRed + 1);
    gMYieldPhiRun3[0]->Draw("APE Y+"); // Y+ places Y-axis on the right side

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
        gEPOS_PhidNdy_UrQMDOff->SetPointError(i, 0, 0);
        gEPOS_PhidNdy_UrQMDOn->SetPointError(i, 0, 0);
    }
    gEPOS_PhidNdy_UrQMDOff->SetLineColor(kGreen + 2);
    gEPOS_PhidNdy_UrQMDOff->SetLineStyle(2);
    gEPOS_PhidNdy_UrQMDOff->SetLineWidth(3);
    gEPOS_PhidNdy_UrQMDOff->Draw("l same");

    gEPOS_PhidNdy_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_PhidNdy_UrQMDOn->SetLineStyle(9);
    gEPOS_PhidNdy_UrQMDOn->SetLineWidth(3);
    gEPOS_PhidNdy_UrQMDOn->Draw("l same");

    TLegend *legendModels = new TLegend(0.05, 0.65, 0.58, 0.92);
    legendModels->SetBorderSize(0);
    legendModels->SetFillStyle(0);
    legendModels->SetTextFont(42);
    legendModels->SetTextSize(0.035);
    legendModels->AddEntry(gEPOS_KstardNdy_UrQMDOff, "EPOS4 UrQMD OFF", "l");
    legendModels->AddEntry(gEPOS_KstardNdy_UrQMDOn, "EPOS4 UrQMD ON", "l");

    for (auto model : modelsToPlot)
    {
        gMYield[model][kPhi]->Draw("l same");
        legendModels->AddEntry(gMYield[model][kPhi], modelLabel[model], "l");
    }
    legendModels->Draw("same");

    latex.SetTextFont(22);
    latex.DrawLatex(0.55, 0.82, "#phi (1020)");

    if (isSavePlots)
    {
        cdNdy_SideBySide->SaveAs("Plots/MeanYield_Kstar_Phi_SideBySide.png");
    }

    // ========================================================================
    // CHANGE 1: Side-by-side <pT> plot of K* and Phi with zero gap
    // ========================================================================
    TCanvas *cMeanPt_SideBySide = new TCanvas("cMeanPt_SideBySide", "<pT> K* and Phi", 1000, 600);
    cMeanPt_SideBySide->Divide(2, 1);
    padLeft->Draw();
    padRight->Draw();

    // --- Drawing Left Pad (K*) ---
    padLeft->cd();
    gMPtKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMPtKstar[0]->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    gMPtKstar[0]->GetXaxis()->SetLimits(0, 29);
    gMPtKstar[0]->GetYaxis()->SetRangeUser(0.62, 1.82);
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

    gMPtKstarRun2[0]->SetMarkerStyle(21);
    gMPtKstarRun2[0]->SetMarkerColor(kBlue + 1);
    gMPtKstarRun2[0]->SetLineColor(kBlue + 1);
    gMPtKstarRun2[0]->SetLineWidth(3);
    gMPtKstarRun2[0]->Draw("P same");

    gMPtKstarRun2[1]->SetLineColor(kBlue + 1);
    gMPtKstarRun2[1]->SetFillStyle(0);
    gMPtKstarRun2[1]->SetLineWidth(3);
    gMPtKstarRun2[1]->Draw("5 same");

    gMPtKstarRun2[2]->SetFillStyle(1001);
    gMPtKstarRun2[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gMPtKstarRun2[2]->Draw("E2 same");

    int totalPointsEPOS_MPt = gEPOS_KstarMeanpT_UrQMDOff->GetN();
    for (int i = 0; i < totalPointsEPOS_MPt; i++)
    {
        gEPOS_KstarMeanpT_UrQMDOff->SetPointError(i, 0, 0);
        gEPOS_KstarMeanpT_UrQMDOn->SetPointError(i, 0, 0);
    }
    gEPOS_KstarMeanpT_UrQMDOff->SetLineColor(kGreen + 2);
    gEPOS_KstarMeanpT_UrQMDOff->SetLineStyle(2);
    gEPOS_KstarMeanpT_UrQMDOff->SetLineWidth(3);
    gEPOS_KstarMeanpT_UrQMDOff->Draw("l same");

    gEPOS_KstarMeanpT_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_KstarMeanpT_UrQMDOn->SetLineStyle(9);
    gEPOS_KstarMeanpT_UrQMDOn->SetLineWidth(3);
    gEPOS_KstarMeanpT_UrQMDOn->Draw("l same");

    counter = 0;
    for (auto model : modelsToPlot)
    {
        gMeanPt[model][kKstar]->Draw("l same");
        counter++;
    }

    latex.SetTextFont(42);
    latex.DrawLatex(0.24, 0.88, "ALICE");
    latex.DrawLatex(0.24, 0.82, "|y| < 0.5");
    latex.SetTextFont(22);
    latex.DrawLatex(0.24, 0.76, "K* (892)^{0}");

    legendData->Draw("same");

    // --- Drawing Right Pad (Phi) ---
    padRight->cd();
    gMPtPhiRun3[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMPtPhiRun3[0]->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    gMPtPhiRun3[0]->GetXaxis()->SetLimits(-1, 29);
    gMPtPhiRun3[0]->GetYaxis()->SetRangeUser(0.62, 1.82);
    gMPtPhiRun3[0]->SetLineWidth(3);
    gMPtPhiRun3[0]->SetMarkerColor(kBlue + 1);
    gMPtPhiRun3[0]->SetLineColor(kBlue + 1);
    gMPtPhiRun3[0]->Draw("APE Y+"); // Y+ places Y-axis on the right side

    gMPtPhiRun3[1]->SetFillStyle(0);
    gMPtPhiRun3[1]->SetLineColor(kBlue + 1);
    gMPtPhiRun3[1]->SetLineWidth(3);
    gMPtPhiRun3[1]->Draw("5 same");

    gMPtPhiRun3[2]->SetFillStyle(1001);
    gMPtPhiRun3[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gMPtPhiRun3[2]->Draw("E2 same");

    gMPtPhiRun2[0]->SetMarkerStyle(21);
    gMPtPhiRun2[0]->SetMarkerColor(kRed + 1);
    gMPtPhiRun2[0]->SetLineColor(kRed + 1);
    gMPtPhiRun2[0]->SetLineWidth(3);
    gMPtPhiRun2[0]->Draw("P same");

    gMPtPhiRun2[1]->SetLineColor(kRed + 1);
    gMPtPhiRun2[1]->SetFillStyle(0);
    gMPtPhiRun2[1]->SetLineWidth(3);
    gMPtPhiRun2[1]->Draw("5 same");

    gMPtPhiRun2[2]->SetFillStyle(1001);
    gMPtPhiRun2[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gMPtPhiRun2[2]->Draw("E2 same");

    int totalPointsEPOS_MPt_Phi = gEPOS_PhiMeanpT_UrQMDOff->GetN();
    for (int i = 0; i < totalPointsEPOS_MPt_Phi; i++)
    {
        gEPOS_PhiMeanpT_UrQMDOff->SetPointError(i, 0, 0);
        gEPOS_PhiMeanpT_UrQMDOn->SetPointError(i, 0, 0);
    }
    gEPOS_PhiMeanpT_UrQMDOff->SetLineColor(kGreen + 2);
    gEPOS_PhiMeanpT_UrQMDOff->SetLineStyle(2);
    gEPOS_PhiMeanpT_UrQMDOff->SetLineWidth(3);
    gEPOS_PhiMeanpT_UrQMDOff->Draw("l same");

    gEPOS_PhiMeanpT_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_PhiMeanpT_UrQMDOn->SetLineStyle(9);
    gEPOS_PhiMeanpT_UrQMDOn->SetLineWidth(3);
    gEPOS_PhiMeanpT_UrQMDOn->Draw("l same");

    TLegend *legendModels_MPt = new TLegend(0.05, 0.65, 0.58, 0.92);
    legendModels_MPt->SetBorderSize(0);
    legendModels_MPt->SetFillStyle(0);
    legendModels_MPt->SetTextFont(42);
    legendModels_MPt->SetTextSize(0.035);
    legendModels_MPt->AddEntry(gEPOS_KstarMeanpT_UrQMDOff, "EPOS4 UrQMD OFF", "l");
    legendModels_MPt->AddEntry(gEPOS_KstarMeanpT_UrQMDOn, "EPOS4 UrQMD ON", "l");

    counter = 0;
    for (auto model : modelsToPlot)
    {
        gMeanPt[model][kPhi]->Draw("l same");
        legendModels_MPt->AddEntry(gMeanPt[model][kPhi], modelLabel[model], "l");
        counter++;
    }
    legendModels_MPt->Draw("same");

    latex.SetTextFont(22);
    latex.DrawLatex(0.55, 0.82, "#phi (1020)");

    if (isSavePlots)
    {
        cMeanPt_SideBySide->SaveAs("Plots/MeanPt_Kstar_Phi_SideBySide.png");
    }

    // //=========================================================
    // // ===== dN/dy normalized to MB vs dNch/deta normalized MB
    // //=========================================================

    // double xMB, yieldMB;
    // gMeanYieldMB_stat->GetPoint(0, xMB, yieldMB);
    // const double dNchdetaMB = 7.087;
    // const double dNchdetaMBErr = 0.18;
    // double yieldMBStatErr = gMeanYieldMB_stat->GetErrorY(0);

    // TGraphErrors *gNormYieldKstarToMB_stat = NormalizeYieldByMB(gMYieldKstar[0], yieldMB, yieldMBStatErr, dNchdetaMB, dNchdetaMBErr, "gNormYieldKstarToMB_stat");
    // TGraphAsymmErrors *gNormYieldKstarToMB_sys = NormalizeYieldSysByMB(gMYieldKstar[1], gMeanYieldMB_sys, yieldMB, dNchdetaMB, dNchdetaMBErr, "gNormYieldKstarToMB_sys");
    // TGraphAsymmErrors *gNormYieldKstarToMB_sysUncorr = new TGraphAsymmErrors();
    // gNormYieldKstarToMB_sysUncorr->SetName("gNormYieldKstarToMB_sysuncorr");

    // TGraphErrors *gNormYieldKstarToMBRun2_stat = NormalizeYieldByMB(gKstar_Yield_13TeV[0], yieldMB, yieldMBStatErr, dNchdetaMB, dNchdetaMBErr, "gNormYieldKstarToMBRun2_stat");
    // TGraphAsymmErrors *gNormYieldKstarToMBRun2_sys = NormalizeYieldSysByMB(gKstar_Yield_13TeV[1], gMeanYieldMB_sys, yieldMB, dNchdetaMB, dNchdetaMBErr, "gNormYieldKstarToMBRun2_sys");
    // TGraphAsymmErrors *gNormYieldKstarToMBRun2_sysUncorr = new TGraphAsymmErrors();
    // gNormYieldKstarToMBRun2_sysUncorr->SetName("gNormYieldKstarToMBRun2_sysuncorr");

    // for (int i = 0; i < gMYieldKstar[2]->GetN(); i++)
    // {
    //     double x, y;
    //     gMYieldKstar[2]->GetPoint(i, x, y);
    //     double ex = gMYieldKstar[2]->GetErrorX(i);
    //     double ey = gMYieldKstar[2]->GetErrorY(i);
    //     double xNorm = x / dNchdetaMB;
    //     double yNorm = y / yieldMB;
    //     double exNorm = sqrt(pow(ex / dNchdetaMB, 2) + pow(x * dNchdetaMBErr / (dNchdetaMB * dNchdetaMB), 2));
    //     double eyNorm = ey / yieldMB;
    //     int ip = gNormYieldKstarToMB_sysUncorr->GetN();
    //     gNormYieldKstarToMB_sysUncorr->SetPoint(ip, xNorm, yNorm);
    //     gNormYieldKstarToMB_sysUncorr->SetPointError(ip, exNorm, exNorm, eyNorm, eyNorm);
    // }

    // for (int i = 0; i < gKstar_Yield_13TeV[2]->GetN(); i++)
    // {
    //     double x, y;
    //     gKstar_Yield_13TeV[2]->GetPoint(i, x, y);
    //     double ex = gKstar_Yield_13TeV[2]->GetErrorX(i);
    //     double ey = gKstar_Yield_13TeV[2]->GetErrorY(i);
    //     double xNorm = x / dNchdetaMB;
    //     double yNorm = y / yieldMB;
    //     double exNorm = sqrt(pow(ex / dNchdetaMB, 2) + pow(x * dNchdetaMBErr / (dNchdetaMB * dNchdetaMB), 2));
    //     double eyNorm = ey / yieldMB;
    //     int ip = gNormYieldKstarToMBRun2_sysUncorr->GetN();
    //     gNormYieldKstarToMBRun2_sysUncorr->SetPoint(ip, xNorm, yNorm);
    //     gNormYieldKstarToMBRun2_sysUncorr->SetPointError(ip, exNorm, exNorm, eyNorm, eyNorm);
    // }

    // // Fit the Run-2 normalized yield and evaluate it at the Run-3 multiplicity points.
    // double fitXmin = 1e9, fitXmax = -1e9;
    // for (int i = 0; i < gNormYieldKstarToMBRun2_stat->GetN(); i++)
    // {
    //     double x, y;
    //     gNormYieldKstarToMBRun2_stat->GetPoint(i, x, y);
    //     if (y > 0)
    //     {
    //         fitXmin = TMath::Min(fitXmin, x);
    //         fitXmax = TMath::Max(fitXmax, x);
    //     }
    // }

    // if (fitXmin >= fitXmax)
    // {
    //     cout << "Error: insufficient positive Run-2 points for the normalized-yield fit." << endl;
    //     return;
    // }

    // TF1 *fitNormYieldRun2 = new TF1("fitNormYieldRun2", "pol1", fitXmin, fitXmax);
    // TFitResultPtr fitResultRun2 = gNormYieldKstarToMBRun2_stat->Fit(fitNormYieldRun2, "QRS0");
    // fitNormYieldRun2->SetLineColor(kBlue + 1);
    // fitNormYieldRun2->SetLineStyle(2);
    // fitNormYieldRun2->SetLineWidth(2);

    // // The MB yield is the common normalization denominator for Run 2 and Run 3.
    // // Its relative systematic uncertainty therefore cancels in Run-3 / Run-2.
    // double mbRelSysLow = gMeanYieldMB_sys->GetErrorYlow(0) / yieldMB;
    // double mbRelSysHigh = gMeanYieldMB_sys->GetErrorYhigh(0) / yieldMB;

    // // Interpolate the numerator-only systematic uncertainty after removing the
    // // common MB normalization component in quadrature from the normalized graph.
    // auto GetInterpolatedSys = [mbRelSysLow, mbRelSysHigh](TGraphAsymmErrors *gr, double x, double &errLow, double &errHigh)
    // {
    //     errLow = 0.0;
    //     errHigh = 0.0;
    //     if (!gr || gr->GetN() < 1)
    //         return;

    //     int iLow = -1, iHigh = -1;
    //     double xLow = -1e30, xHigh = 1e30;
    //     for (int ip = 0; ip < gr->GetN(); ip++)
    //     {
    //         double gx, gy;
    //         gr->GetPoint(ip, gx, gy);
    //         if (gx <= x && gx > xLow)
    //         {
    //             xLow = gx;
    //             iLow = ip;
    //         }
    //         if (gx >= x && gx < xHigh)
    //         {
    //             xHigh = gx;
    //             iHigh = ip;
    //         }
    //     }

    //     if (iLow < 0 || iHigh < 0)
    //         return;

    //     auto numeratorOnlyError = [&](int ip, bool isLow)
    //     {
    //         double gx, gy;
    //         gr->GetPoint(ip, gx, gy);
    //         double totalErr = isLow ? gr->GetErrorYlow(ip) : gr->GetErrorYhigh(ip);
    //         double commonMBErr = gy * (isLow ? mbRelSysLow : mbRelSysHigh);
    //         return sqrt(TMath::Max(0.0, totalErr * totalErr - commonMBErr * commonMBErr));
    //     };

    //     if (iLow == iHigh || TMath::AreEqualAbs(xLow, xHigh, 1e-12))
    //     {
    //         errLow = numeratorOnlyError(iLow, true);
    //         errHigh = numeratorOnlyError(iLow, false);
    //         return;
    //     }

    //     double w = (x - xLow) / (xHigh - xLow);
    //     errLow = (1.0 - w) * numeratorOnlyError(iLow, true) + w * numeratorOnlyError(iHigh, true);
    //     errHigh = (1.0 - w) * numeratorOnlyError(iLow, false) + w * numeratorOnlyError(iHigh, false);
    // };

    // TGraphErrors *gRatioNormYieldRun3Run2_stat = new TGraphErrors();
    // gRatioNormYieldRun3Run2_stat->SetName("gRatioNormYieldRun3Run2_stat");
    // TGraphAsymmErrors *gRatioNormYieldRun3Run2_sys = new TGraphAsymmErrors();
    // gRatioNormYieldRun3Run2_sys->SetName("gRatioNormYieldRun3Run2_sys");

    // for (int i = 0; i < gNormYieldKstarToMB_stat->GetN(); i++)
    // {
    //     double x, yRun3;
    //     gNormYieldKstarToMB_stat->GetPoint(i, x, yRun3);
    //     if (x < fitXmin || x > fitXmax)
    //         continue;

    //     double yRun2 = fitNormYieldRun2->Eval(x);
    //     if (yRun3 <= 0 || yRun2 <= 0)
    //         continue;

    //     // Statistical uncertainty: Run-3 point and Run-2 fit covariance.
    //     double fitVariance = 0.0;
    //     if (fitResultRun2.Get() && fitResultRun2->NPar() >= 2)
    //         fitVariance = fitResultRun2->CovMatrix(0, 0) + x * x * fitResultRun2->CovMatrix(1, 1) + 2.0 * x * fitResultRun2->CovMatrix(0, 1);

    //     double eyRun3Stat = gNormYieldKstarToMB_stat->GetErrorY(i);
    //     double eyRun2Stat = sqrt(TMath::Max(0.0, fitVariance));
    //     double ratio = yRun3 / yRun2;
    //     double ratioStatErr = ratio * sqrt(pow(eyRun3Stat / yRun3, 2) + pow(eyRun2Stat / yRun2, 2));

    //     int ip = gRatioNormYieldRun3Run2_stat->GetN();
    //     gRatioNormYieldRun3Run2_stat->SetPoint(ip, x, ratio);
    //     gRatioNormYieldRun3Run2_stat->SetPointError(ip, gNormYieldKstarToMB_stat->GetErrorX(i), ratioStatErr);

    //     // Total systematic uncertainty: combine Run-3 and interpolated Run-2
    //     // relative systematic uncertainties in quadrature.
    //     double run3SysLowTotal = gNormYieldKstarToMB_sys->GetErrorYlow(i);
    //     double run3SysHighTotal = gNormYieldKstarToMB_sys->GetErrorYhigh(i);
    //     double run3SysLowMB = yRun3 * mbRelSysLow;
    //     double run3SysHighMB = yRun3 * mbRelSysHigh;
    //     double run3SysLow = sqrt(TMath::Max(0.0, run3SysLowTotal * run3SysLowTotal - run3SysLowMB * run3SysLowMB));
    //     double run3SysHigh = sqrt(TMath::Max(0.0, run3SysHighTotal * run3SysHighTotal - run3SysHighMB * run3SysHighMB));

    //     double run2SysLow = 0.0, run2SysHigh = 0.0;
    //     GetInterpolatedSys(gNormYieldKstarToMBRun2_sys, x, run2SysLow, run2SysHigh);

    //     double ratioSysLow = ratio * sqrt(pow(run3SysLow / yRun3, 2) + pow(run2SysLow / yRun2, 2));
    //     double ratioSysHigh = ratio * sqrt(pow(run3SysHigh / yRun3, 2) + pow(run2SysHigh / yRun2, 2));

    //     gRatioNormYieldRun3Run2_sys->SetPoint(ip, x, ratio);
    //     gRatioNormYieldRun3Run2_sys->SetPointError(ip, gNormYieldKstarToMB_stat->GetErrorX(i), gNormYieldKstarToMB_stat->GetErrorX(i), ratioSysLow, ratioSysHigh);
    // }

    // TCanvas *cNormYieldKstarToMB = new TCanvas("cNormYieldKstarToMB", "cNormYieldKstarToMB", 720, 800);
    // SetCanvasStyle(cNormYieldKstarToMB, 0.15, 0.05, 0.05, 0.15);
    // SetGraphErrorStyle(gNormYieldKstarToMB_stat);
    // gNormYieldKstarToMB_stat->GetXaxis()->SetTitle("");
    // gNormYieldKstarToMB_stat->GetYaxis()->SetTitle("(dN/dy) / (dN/dy)_{MB}");
    // gNormYieldKstarToMB_stat->GetXaxis()->SetLimits(0, 3.9);
    // gNormYieldKstarToMB_stat->GetYaxis()->SetRangeUser(0, 4.4);
    // gNormYieldKstarToMB_stat->GetYaxis()->SetTitleOffset(1.55);
    // gNormYieldKstarToMB_stat->GetXaxis()->SetLabelSize(0);
    // gNormYieldKstarToMB_stat->SetMarkerColor(kRed + 1);
    // gNormYieldKstarToMB_stat->SetLineColor(kRed + 1);
    // gNormYieldKstarToMB_stat->SetMarkerStyle(20);
    // gNormYieldKstarToMB_stat->SetMarkerSize(1.2);
    // gNormYieldKstarToMB_stat->SetLineWidth(2);
    // gNormYieldKstarToMB_stat->Draw("APE");

    // gNormYieldKstarToMB_sys->SetFillStyle(0);
    // gNormYieldKstarToMB_sys->SetLineColor(kRed + 1);
    // gNormYieldKstarToMB_sys->SetLineWidth(2);
    // gNormYieldKstarToMB_sys->Draw("5 same");
    // gNormYieldKstarToMB_sysUncorr->SetFillColorAlpha(kRed + 1, 0.3);
    // gNormYieldKstarToMB_sysUncorr->SetFillStyle(1001);
    // gNormYieldKstarToMB_sysUncorr->SetLineColor(kRed + 1);
    // gNormYieldKstarToMB_sysUncorr->Draw("E2 same");

    // gNormYieldKstarToMBRun2_stat->SetMarkerStyle(21);
    // gNormYieldKstarToMBRun2_stat->SetMarkerColor(kBlue + 1);
    // gNormYieldKstarToMBRun2_stat->SetLineColor(kBlue + 1);
    // gNormYieldKstarToMBRun2_stat->SetMarkerSize(1.2);
    // gNormYieldKstarToMBRun2_stat->SetLineWidth(2);
    // gNormYieldKstarToMBRun2_stat->Draw("P same");
    // gNormYieldKstarToMBRun2_sys->SetFillStyle(0);
    // gNormYieldKstarToMBRun2_sys->SetLineColor(kBlue + 1);
    // gNormYieldKstarToMBRun2_sys->SetLineWidth(2);
    // gNormYieldKstarToMBRun2_sys->Draw("5 same");
    // gNormYieldKstarToMBRun2_sysUncorr->SetFillColorAlpha(kBlue + 1, 0.3);
    // gNormYieldKstarToMBRun2_sysUncorr->SetFillStyle(1001);
    // gNormYieldKstarToMBRun2_sysUncorr->SetLineColor(kBlue + 1);
    // gNormYieldKstarToMBRun2_sysUncorr->Draw("E2 same");
    // fitNormYieldRun2->Draw("same");

    // TLegend *legNormYieldMB = new TLegend(0.42, 0.74, 0.82, 0.90);
    // SetLegendStyle(legNormYieldMB);
    // legNormYieldMB->AddEntry(gNormYieldKstarToMB_stat, "Run 3, #sqrt{s} = 13.6 TeV", "P");
    // legNormYieldMB->AddEntry(gNormYieldKstarToMBRun2_stat, "Run 2, #sqrt{s} = 13 TeV", "P");
    // legNormYieldMB->AddEntry(fitNormYieldRun2, "Run 2 linear fit", "L");
    // legNormYieldMB->Draw();

    // latex.DrawLatex(0.22, 0.9, "ALICE");
    // latex.SetTextFont(22);
    // latex.DrawLatex(0.22, 0.84, "K* (892)^{0}");
    // latex.SetTextFont(42);
    // latex.DrawLatex(0.22, 0.78, "|y| < 0.5");

    // if (isSavePlots)
    //     cNormYieldKstarToMB->SaveAs("Plots/NormalizedYieldToMB_Kstar_Run3Run2.png");

    //=========================================================
    // Normalized dN/dy: Kstar and Phi side-by-side
    //=========================================================

    const double dNchdetaMB = 7.087;
    const double dNchdetaMBErr = 0.18;

    // ---------------- Kstar minimum-bias normalization ----------------
    double xMB_Kstar, yieldMB_Kstar;
    gMeanYieldMB_stat->GetPoint(0, xMB_Kstar, yieldMB_Kstar);
    double yieldMBStatErr_Kstar = gMeanYieldMB_stat->GetErrorY(0);

    TGraphErrors *gNormKstar_stat = NormalizeYieldByMB(gMYieldKstar[0], yieldMB_Kstar, yieldMBStatErr_Kstar, dNchdetaMB, dNchdetaMBErr, "gNormKstar_stat");

    TGraphAsymmErrors *gNormKstar_sys = NormalizeYieldSysByMB(gMYieldKstar[1], gMeanYieldMB_sys, yieldMB_Kstar, dNchdetaMB, dNchdetaMBErr, "gNormKstar_sys");

    // ---------------- Phi minimum-bias normalization ----------------
    double xMB_Phi, yieldMB_Phi;
    gMeanYieldMB_Phi_stat->GetPoint(0, xMB_Phi, yieldMB_Phi);
    double yieldMBStatErr_Phi = gMeanYieldMB_Phi_stat->GetErrorY(0);

    // Convert Phi MB systematic graph into TGraphAsymmErrors
    // because NormalizeYieldSysByMB() accepts TGraphAsymmErrors.
    TGraphAsymmErrors *gMeanYieldMB_Phi_sys_asymm = new TGraphAsymmErrors();

    for (int i = 0; i < gMeanYieldMB_Phi_sys->GetN(); i++)
    {
        double x, y;
        gMeanYieldMB_Phi_sys->GetPoint(i, x, y);
        double ex = gMeanYieldMB_Phi_sys->GetErrorX(i);
        double ey = gMeanYieldMB_Phi_sys->GetErrorY(i);

        gMeanYieldMB_Phi_sys_asymm->SetPoint(i, x, y);
        gMeanYieldMB_Phi_sys_asymm->SetPointError(i, ex, ex, ey, ey);
    }

    TGraphErrors *gNormPhi_stat = NormalizeYieldByMB(gMYieldPhiRun3[0], yieldMB_Phi, yieldMBStatErr_Phi, dNchdetaMB, dNchdetaMBErr, "gNormPhi_stat");

    TGraphAsymmErrors *gNormPhi_sys = NormalizeYieldSysByMB(gMYieldPhiRun3[1], gMeanYieldMB_Phi_sys_asymm, yieldMB_Phi, dNchdetaMB, dNchdetaMBErr, "gNormPhi_sys");

    if (!gNormKstar_stat || !gNormKstar_sys || !gNormPhi_stat || !gNormPhi_sys)
    {
        cout << "Error: Failed to create normalized yield graphs." << endl;
        return;
    }

    // ---------------- Canvas ----------------
    TCanvas *cNormYield_SideBySide = new TCanvas("cNormYield_SideBySide", "Normalized dN/dy", 1200, 550);
    cNormYield_SideBySide->Divide(2, 1, 0.01, 0.01);

    TLatex latexNorm;
    latexNorm.SetNDC();
    latexNorm.SetTextFont(42);
    latexNorm.SetTextSize(0.045);

    // ---------------- Left: Kstar ----------------
    cNormYield_SideBySide->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    gPad->SetTicks(1, 1);

    gNormKstar_stat->SetMarkerStyle(20);
    gNormKstar_stat->SetMarkerSize(1.1);
    gNormKstar_stat->SetMarkerColor(kRed + 1);
    gNormKstar_stat->SetLineColor(kRed + 1);
    gNormKstar_stat->SetLineWidth(2);

    gNormKstar_stat->SetTitle("");
    gNormKstar_stat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT / #LT dN_{ch}/d#eta #GT_{MB}");
    gNormKstar_stat->GetYaxis()->SetTitle("(dN/dy) / (dN/dy)_{MB}");
    gNormKstar_stat->GetXaxis()->SetLimits(0, 4.5);
    gNormKstar_stat->GetYaxis()->SetRangeUser(0, 2.5);
    gNormKstar_stat->GetXaxis()->SetTitleSize(0.045);
    gNormKstar_stat->GetYaxis()->SetTitleSize(0.045);
    gNormKstar_stat->GetXaxis()->SetLabelSize(0.04);
    gNormKstar_stat->GetYaxis()->SetLabelSize(0.04);
    gNormKstar_stat->GetYaxis()->SetTitleOffset(1.45);

    gNormKstar_stat->Draw("APE");

    gNormKstar_sys->SetFillColorAlpha(kRed + 1, 0.25);
    gNormKstar_sys->SetLineColor(kRed + 1);
    gNormKstar_sys->SetFillStyle(1001);
    gNormKstar_sys->Draw("5 same");

    gNormKstar_stat->Draw("P same");

    latexNorm.DrawLatex(0.22, 0.88, "ALICE");
    latexNorm.DrawLatex(0.22, 0.82, "pp, #sqrt{s} = 13.6 TeV");
    latexNorm.SetTextFont(22);
    latexNorm.DrawLatex(0.22, 0.75, "K* (892)^{0}");

    // ---------------- Right: Phi ----------------
    cNormYield_SideBySide->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    gPad->SetTicks(1, 1);

    gNormPhi_stat->SetMarkerStyle(20);
    gNormPhi_stat->SetMarkerSize(1.1);
    gNormPhi_stat->SetMarkerColor(kBlue + 1);
    gNormPhi_stat->SetLineColor(kBlue + 1);
    gNormPhi_stat->SetLineWidth(2);

    gNormPhi_stat->SetTitle("");
    gNormPhi_stat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT / #LT dN_{ch}/d#eta #GT_{MB}");
    gNormPhi_stat->GetYaxis()->SetTitle("(dN/dy) / (dN/dy)_{MB}");
    gNormPhi_stat->GetXaxis()->SetLimits(0, 4.5);
    gNormPhi_stat->GetYaxis()->SetRangeUser(0, 2.5);
    gNormPhi_stat->GetXaxis()->SetTitleSize(0.045);
    gNormPhi_stat->GetYaxis()->SetTitleSize(0.045);
    gNormPhi_stat->GetXaxis()->SetLabelSize(0.04);
    gNormPhi_stat->GetYaxis()->SetLabelSize(0.04);
    gNormPhi_stat->GetYaxis()->SetTitleOffset(1.45);

    gNormPhi_stat->Draw("APE");

    gNormPhi_sys->SetFillColorAlpha(kBlue + 1, 0.25);
    gNormPhi_sys->SetLineColor(kBlue + 1);
    gNormPhi_sys->SetFillStyle(1001);
    gNormPhi_sys->Draw("5 same");

    gNormPhi_stat->Draw("P same");

    latexNorm.SetTextFont(42);
    latexNorm.DrawLatex(0.22, 0.88, "ALICE");
    latexNorm.DrawLatex(0.22, 0.82, "pp, #sqrt{s} = 13.6 TeV");
    latexNorm.SetTextFont(22);
    latexNorm.DrawLatex(0.22, 0.75, "#phi (1020)");

    if (isSavePlots)
    {
        cNormYield_SideBySide->SaveAs("Plots/Normalized_dNdy_Kstar_Phi_SideBySide.png");
    }

    // ========================================================================
    // Prepare Ratios for 5-Panel Plot
    // ========================================================================
    TFile *fEPOS_New2_Ratio = OpenFile("/home/sawan/Storage/EPOS4_Resonance/FinalEPOSResults_BongFix3.root");
    TGraphErrors *gdNdy_Kstar_On = GetGraph(fEPOS_New2_Ratio, "g_dNdy_opt1");
    TGraphErrors *gdNdy_Kstar_Off = GetGraph(fEPOS_New2_Ratio, "g_dNdy_opt0");
    TGraphErrors *gdNdy_Phi_On = GetGraph(fEPOS_New2_Ratio, "g_phi_dNdy_opt1");
    TGraphErrors *gdNdy_Phi_Off = GetGraph(fEPOS_New2_Ratio, "g_phi_dNdy_opt0");
    TGraphErrors *gdNdy_K0s_On = GetGraph(fEPOS_New2_Ratio, "g_stable_dNdy_opt1_K0s");
    TGraphErrors *gdNdy_K0s_Off = GetGraph(fEPOS_New2_Ratio, "g_stable_dNdy_opt0_K0s");
    TGraphErrors *gdNdy_Pion_On = GetGraph(fEPOS_New2_Ratio, "g_stable_dNdy_opt1_pi");
    TGraphErrors *gdNdy_Pion_Off = GetGraph(fEPOS_New2_Ratio, "g_stable_dNdy_opt0_pi");
    TGraphErrors *gdNdy_Xi_On = GetGraph(fEPOS_New2_Ratio, "g_stable_dNdy_opt1_Xi");
    TGraphErrors *gdNdy_Xi_Off = GetGraph(fEPOS_New2_Ratio, "g_stable_dNdy_opt0_Xi");

    TGraphErrors *gKstarPiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarPiRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gKstarK0sRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarK0sRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gPhiPiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gPhiPiRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gPhiK0sRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gPhiK0sRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gPhiXiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gPhiXiRatio_EPOSNew_UrQMDOn = new TGraphErrors();

    int graphPointsEPOSNew = gdNdy_Kstar_On->GetN();
    for (int i = 0; i < graphPointsEPOSNew; i++)
    {
        double xK0s, yK0s, yPion, yKstarOff, yKstarOn;
        gdNdy_Kstar_Off->GetPoint(i, xK0s, yKstarOff);
        gdNdy_Kstar_On->GetPoint(i, xK0s, yKstarOn);
        gdNdy_K0s_Off->GetPoint(i, xK0s, yK0s);
        gdNdy_Pion_Off->GetPoint(i, xK0s, yPion);

        gKstarPiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yPion);
        gKstarPiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yPion);
        gKstarK0sRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yK0s);
        gKstarK0sRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yK0s);
    }

    int graphPointsEPOSNewPhi = gdNdy_Phi_On->GetN();
    for (int i = 0; i < graphPointsEPOSNewPhi; i++)
    {
        double xK0s, yK0s, yPion, yPhiOff, yPhiOn, yXiOff;
        gdNdy_Phi_Off->GetPoint(i, xK0s, yPhiOff);
        gdNdy_Phi_On->GetPoint(i, xK0s, yPhiOn);
        gdNdy_K0s_Off->GetPoint(i, xK0s, yK0s);
        gdNdy_Pion_Off->GetPoint(i, xK0s, yPion);
        gdNdy_Xi_Off->GetPoint(i, xK0s, yXiOff);

        gPhiPiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yPion);
        gPhiPiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yPion);
        gPhiK0sRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yK0s);
        gPhiK0sRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yK0s);
        gPhiXiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yXiOff);
        gPhiXiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yXiOff);
    }

    TGraphErrors *gKstarPiRatio_UrQMDOffSmooth = SmoothGraph(gKstarPiRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarPiRatio_UrQMDOnSmooth = SmoothGraph(gKstarPiRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gKstarK0sRatio_UrQMDOffSmooth = SmoothGraph(gKstarK0sRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarK0sRatio_UrQMDOnSmooth = SmoothGraph(gKstarK0sRatio_EPOSNew_UrQMDOn);

    TGraphErrors *gPhiPiRatio_UrQMDOffSmooth = SmoothGraph(gPhiPiRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gPhiPiRatio_UrQMDOnSmooth = SmoothGraph(gPhiPiRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gPhiK0sRatio_UrQMDOffSmooth = SmoothGraph(gPhiK0sRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gPhiK0sRatio_UrQMDOnSmooth = SmoothGraph(gPhiK0sRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gPhiXiRatio_UrQMDOffSmooth = SmoothGraph(gPhiXiRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gPhiXiRatio_UrQMDOnSmooth = SmoothGraph(gPhiXiRatio_EPOSNew_UrQMDOn);

    TGraphErrors **gRatioKstarKshort = MakeRatioUncorr(gMYieldKstar, gMYieldKshortRun3, 1.0, 3);
    TGraphErrors **gRatioPhiK0s = MakeRatioUncorr(gMYieldPhiRun3, gMYieldKshortRun3, 1.0, 3);
    TGraphErrors **gRatioKstarPion = MakeRatioUncorr(gMYieldKstar, gMYieldPion, 1.0, 3);
    TGraphErrors **gRatioPhiPion = MakeRatioUncorr(gMYieldPhiRun3, gMYieldPion, 1.0, 3);
    TGraphErrors **gRatioPhiXi = MakeRatioUncorr(gMYieldPhiRun3, gMYieldXiRun3, 1.0, 3);

    // ========================================================================
    // CHANGE 2: 4-Panel Plot (2x2 ALICE-style layout)
    // ========================================================================

    TCanvas *c5Panel = new TCanvas("c5Panel", "Particle Ratios 4-Panel Plot", 1200, 950);

    TPad *pad1 = new TPad("pad1", "Kstar/K0s", 0.0, 0.495, 0.5, 0.90);
    TPad *pad2 = new TPad("pad2", "Phi/K0s", 0.5, 0.495, 1.0, 0.90);
    TPad *pad3 = new TPad("pad3", "Kstar/Pi", 0.0, 0.05, 0.5, 0.495);
    TPad *pad4 = new TPad("pad4", "Phi/Pi", 0.5, 0.05, 1.0, 0.495);

    TPad *pads[4] = {pad1, pad2, pad3, pad4};

    for (int i = 0; i < 4; i++)
    {
        pads[i]->SetFillStyle(4000);
        pads[i]->SetTicks(1, 1);
        pads[i]->SetBorderMode(0);
        pads[i]->SetFrameBorderMode(0);
        pads[i]->Draw();
    }

    // Left column
    pad1->SetLeftMargin(0.16);
    pad1->SetRightMargin(0.0);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.0);

    pad3->SetLeftMargin(0.16);
    pad3->SetRightMargin(0.0);
    pad3->SetTopMargin(0.0);
    pad3->SetBottomMargin(0.15);

    // Right column
    pad2->SetLeftMargin(0.0);
    pad2->SetRightMargin(0.16);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.0);

    pad4->SetLeftMargin(0.0);
    pad4->SetRightMargin(0.16);
    pad4->SetTopMargin(0.0);
    pad4->SetBottomMargin(0.15);

    auto SetupRatioGraph = [](TGraphErrors **grArr, double yMin, double yMax, const char *ratioTitle, bool rightPanel, bool bottomPanel)
    {
        // grArr[0]->GetXaxis()->SetTitle(bottomPanel ? "#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}" : "");
        // grArr[0]->GetYaxis()->SetTitle(rightPanel ? "" : "dN/dy ratio");

        grArr[0]->GetXaxis()->SetTitle("");
        grArr[0]->GetYaxis()->SetTitle("");

        SetGraphErrorStyle(grArr[0]);

        grArr[0]->GetYaxis()->SetRangeUser(yMin, yMax);
        grArr[0]->GetXaxis()->SetLimits(1.8, 23.5);

        grArr[0]->SetMarkerColor(kRed);
        grArr[0]->SetLineColor(kRed);

        // Axis fonts and sizes
        grArr[0]->GetXaxis()->SetTitleFont(42);
        grArr[0]->GetYaxis()->SetTitleFont(42);
        grArr[0]->GetXaxis()->SetLabelFont(42);
        grArr[0]->GetYaxis()->SetLabelFont(42);

        grArr[0]->GetXaxis()->SetTitleSize(0.055);
        grArr[0]->GetYaxis()->SetTitleSize(0.060);
        grArr[0]->GetXaxis()->SetLabelSize(0.060);
        grArr[0]->GetYaxis()->SetLabelSize(0.060);

        grArr[0]->GetXaxis()->SetTitleOffset(1.25);
        grArr[0]->GetYaxis()->SetTitleOffset(1.25);

        grArr[0]->GetXaxis()->SetNdivisions(505);
        grArr[0]->GetYaxis()->SetNdivisions(505);

        // Hide X-axis labels and title in upper panels
        if (!bottomPanel)
        {
            grArr[0]->GetXaxis()->SetLabelSize(0);
            grArr[0]->GetXaxis()->SetTitleSize(0);
        }

        // Move Y-axis labels to the right for right-hand panels
        if (rightPanel)
        {
            grArr[0]->GetYaxis()->SetLabelSize(0);
            // grArr[0]->GetYaxis()->SetTitleSize(0);
            grArr[0]->GetYaxis()->SetTickLength(0);
        }

        grArr[0]->Draw("APE");

        grArr[1]->SetFillStyle(0);
        grArr[1]->SetLineColor(kRed);
        grArr[1]->Draw("5 same");

        grArr[2]->SetFillStyle(1001);
        grArr[2]->SetFillColorAlpha(kRed, 0.35);
        grArr[2]->Draw("E2 same");
    };

    auto DrawRightAxis = [](double yMin, double yMax)
    {
        TGaxis *axis = new TGaxis(23.5, yMin, 23.5, yMax, yMin, yMax, 505, "+L");

        axis->SetLabelFont(42);
        axis->SetLabelSize(0.045);
        axis->SetTickSize(0.025);
        axis->SetLineWidth(1);
        axis->Draw();
    };

    auto DrawRatioLabel = [](const char *label, double x, double y)
    {
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.055);
        latex.DrawLatex(x, y, label);
    };

    // Helper lambda for smooth EPOS curves
    auto DrawEPOSCurves = [](TGraphErrors *gOff, TGraphErrors *gOn)
    {
        gOff->SetLineColor(kGreen + 2);
        gOff->SetLineStyle(2);
        gOff->SetLineWidth(3);
        gOff->Draw("LX same");

        gOn->SetLineColor(kBlue - 6);
        gOn->SetLineStyle(9);
        gOn->SetLineWidth(3);
        gOn->Draw("LX same");
    };

    // Helper lambda for PYTHIA model curves
    auto DrawModelCurves = [&](int numerator, int denominator, double scale = 1.0)
    {
        counter = 0;

        for (auto model : modelsToPlot)
        {
            TGraphErrors *gRatioModel = MakeRatio(gMYield[model][numerator], gMYield[model][denominator], true);

            if (scale != 1.0)
                ScaleGraph(gRatioModel, scale);

            setStyle(gRatioModel, colorsPythia[counter], lineStylesPythia[counter]);
            gRatioModel->SetLineWidth(3);
            gRatioModel->Draw("L same");

            counter++;
        }
    };

    TLatex panelLatex;
    panelLatex.SetNDC();
    panelLatex.SetTextFont(42);
    panelLatex.SetTextSize(0.065);

    // ---------------- Panel (a): K*/K0s ----------------
    pad1->cd();

    SetupRatioGraph(gRatioKstarKshort, 0.245, 0.54, "K*^{0}/K_{S}^{0}", false, false);

    DrawEPOSCurves(gKstarK0sRatio_EPOSNew_UrQMDOff, gKstarK0sRatio_EPOSNew_UrQMDOn);
    DrawModelCurves(kKstar, kKshort);

    DrawRatioLabel("(a) K*^{0}/K_{S}^{0}", 0.65, 0.88);

    TLatex *latexRightAxis = new TLatex();
    latexRightAxis->SetNDC();
    latexRightAxis->SetTextFont(42);
    latexRightAxis->SetTextSize(0.045);
    latexRightAxis->DrawLatex(0.22, 0.88, "ALICE Preliminary");
    latexRightAxis->DrawLatex(0.22, 0.82, "pp, #sqrt{s} = 13.6 TeV");
    latexRightAxis->DrawLatex(0.22, 0.76, "|y| < 0.5");
    latexRightAxis->Draw();

    // ---------------- Panel (b): Phi/K0s ----------------
    pad2->cd();

    SetupRatioGraph(gRatioPhiK0s, 0.07, 0.18, "#phi/K_{S}^{0}", true, false);

    DrawEPOSCurves(gPhiK0sRatio_EPOSNew_UrQMDOff, gPhiK0sRatio_EPOSNew_UrQMDOn);
    DrawModelCurves(kPhi, kKshort);

    DrawRatioLabel("(b) #phi/K_{S}^{0}", 0.08, 0.88);
    DrawRightAxis(0.07, 0.18);

    // ---------------- Panel (c): K*/Pi ----------------
    pad3->cd();

    SetupRatioGraph(gRatioKstarPion, 0.031, 0.062, "K*^{0}/#pi", false, true);

    DrawEPOSCurves(gKstarPiRatio_EPOSNew_UrQMDOff, gKstarPiRatio_EPOSNew_UrQMDOn);
    DrawModelCurves(kKstar, kPion, 2.0);

    DrawRatioLabel("(c) K*^{0}/#pi", 0.65, 0.88);

    // ---------------- Panel (d): Phi/Pi ----------------
    pad4->cd();

    SetupRatioGraph(gRatioPhiPion, 0.007, 0.028, "#phi/#pi", true, true);

    DrawEPOSCurves(gPhiPiRatio_EPOSNew_UrQMDOff, gPhiPiRatio_EPOSNew_UrQMDOn);
    DrawModelCurves(kPhi, kPion, 2.0);

    DrawRatioLabel("(d) #phi/#pi", 0.08, 0.88);
    DrawRightAxis(0.007, 0.028);

    // ========================================================================
    // Common model legend at the top of the canvas
    // ========================================================================

    c5Panel->cd();

    TLegend *leg5Panel = new TLegend(0.08, 0.895, 0.94, 0.985);
    SetLegendStyle(leg5Panel);
    leg5Panel->SetNColumns(3);
    leg5Panel->SetTextSize(0.027);
    // leg5Panel->SetHeader("Models");

    leg5Panel->AddEntry(gKstarK0sRatio_UrQMDOffSmooth, "EPOS4 UrQMD OFF", "l");
    leg5Panel->AddEntry(gKstarK0sRatio_UrQMDOnSmooth, "EPOS4 UrQMD ON", "l");

    for (int iM = 0; iM < kNModels; iM++)
    {
        TGraphErrors *gDummyModel = new TGraphErrors();
        setStyle(gDummyModel, colorsPythia[iM], lineStylesPythia[iM]);
        gDummyModel->SetLineWidth(3);
        leg5Panel->AddEntry(gDummyModel, modelLabel[iM], "l");
    }

    leg5Panel->Draw();

    // Explicitly keep all axes linear
    for (int i = 0; i < 4; i++)
    {
        pads[i]->SetLogx(0);
        pads[i]->SetLogy(0);
        pads[i]->Modified();
    }

    // ========================================================================
    // Common axis titles
    // ========================================================================

    c5Panel->cd();

    TLatex axisTitle;
    axisTitle.SetNDC();
    axisTitle.SetTextFont(42);
    axisTitle.SetTextSize(0.035);

    // Left Y-axis title (common for panels a and c)
    axisTitle.SetTextAngle(90);
    axisTitle.SetTextAlign(22);
    axisTitle.DrawLatex(0.02, 0.475, "dN/dy ratio");

    // Right Y-axis title (common for panels b and d)
    axisTitle.SetTextAngle(90);
    axisTitle.DrawLatex(0.978, 0.475, "dN/dy ratio");

    // Common X-axis title (for all four panels)
    axisTitle.SetTextAngle(0);
    axisTitle.SetTextAlign(22);
    axisTitle.DrawLatex(0.50, 0.03, "#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");

    c5Panel->Modified();
    c5Panel->Update();

    if (isSavePlots)
    {
        c5Panel->SaveAs("Plots/Ratios_5Panel_3x2Canvas.png");
        // c5Panel->SaveAs("Plots/Ratios_5Panel_3x2Canvas.pdf");
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
