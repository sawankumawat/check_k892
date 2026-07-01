#include <iostream>
#include <iomanip>
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

TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false, float CorrFactorDen = 1);
TGraphErrors **MakeRatioUncorr(TGraphErrors **numerator, TGraphErrors **denominator, float CorrFactorDen = 1, int nGraphs = 3);
TGraphErrors *GetGraph(TFile *f, const string &name);
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
void canvas_style(TCanvas *c, double pad1Size = 0.7, double pad2Size = 0.3, double leftMargin = 0.16, double rightMargin = 0.06, double topMargin = 0.02, double bottomMargin = 0.35);

void setStyle(TGraphErrors *gr, int color, int style)
{
    gr->SetLineColor(color);
    gr->SetLineStyle(style);
}

void ScaleGraph(TGraph *gr, double scale);
void AddGraph(TGraph *gr1, TGraph *gr2);
TGraphErrors **AddConstantToGraph(TGraphErrors **gr, double constant, int noOfGraphs = 2);
TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance = 0.5, int noOfGraphs = 3);
TGraphErrors *DivideByMultModel(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5);
int FindGraphXPoint(TGraphErrors *gr, double xValue);
void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax);

void ParticleRatioWithRun2()
{
    bool isSavePlots = true;
    string KstarPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/";
    TFile *fKstar = OpenFile(KstarPath + "Results.root");
    if (fKstar->IsZombie())
    {
        cout << "Error: Kstar file not found" << endl;
        return;
    }

    TFile *fPion = OpenFile("PiKp_Run3_Results/Sawan/Pi_results.root");
    TFile *fProton = OpenFile("PiKp_Run3_Results/Sawan/Pr_results.root");
    TFile *fKaon = OpenFile("PiKp_Run3_Results/Sawan/Ka_results.root");

    TFile *fPiKp = OpenFile("ConversionCodes/ppRun3_PiKpHEP.root");

    TFile *fKstarRun2 = OpenFile("ConversionCodes/pp13TeVALICE.root");
    TFile *fPhi = OpenFile("ConversionCodes/pp13TeVALICE.root");
    TFile *fChKstar = OpenFile("ConversionCodes/pp13TeV_ChKstar.root");
    TFile *fMultStrange = OpenFile("ConversionCodes/pp13TeV_MultStrange.root");
    TFile *fXiStarSigma = OpenFile("ConversionCodes/pp13TeV_XiStarSigma.root");
    TFile *fLambda1520Run2 = OpenFile("ConversionCodes/pp13TeV_LambdaStar.root");

    // Only available in run3
    TFile *fLambda1520 = OpenFile("LambdaRun3/Sawan/ResultsLambda1520.root");
    TFile *fRho = OpenFile("Rho_Run3_Results/Sawan/ResultsRho.root");

    TGraphErrors *gMPtKstar[3], *gMYieldKstar[3], *gMPtPion[3], *gMPtProton[3], *gMPtKaon[3], *gMPtKstarRun2[3], *gMPtPhiRun2[3], *gMPtChKstarRun2[3], *gMPtKshortRun2[3], *gMPtLambdaRun2[3], *gMPtXiRun2[3], *gMPtOmegaRun2[3], *gMPtXiStarRun2[3], *gMPtSigmaRun2[3], *gMPtLambda1520[3], *gMPtRho[3];
    TGraphErrors *gMYieldPion[3], *gMYieldProton[3], *gMYieldKaon[3], *gMYieldKstarRun2[3], *gMYieldPhiRun2[3], *gMYieldChKstarRun2[3], *gMYieldKshortRun2[3], *gMYieldLambdaRun2[3], *gMYieldXiRun2[3], *gMYieldOmegaRun2[3], *gMYieldXiStarRun2[3], *gMYieldSigmaRun2[3], *gMYieldLambda1520[3], *gMYieldRho[3];

    TGraphErrors *gYieldPhiKaRatio[3], *gYieldChKstarKshortRatio[3], *gYieldSigmaLambdaRatio[3], *gYieldXiStarXiRatio[3], *gYieldLambda1520LambdaRatio[3];

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

        gMPtPion[i] = GetGraph(fPion, Form("gMeanpTRun3%s", suffix.c_str()));
        gMPtKaon[i] = GetGraph(fKaon, Form("gMeanpTRun3%s", suffix.c_str()));
        gMPtProton[i] = GetGraph(fProton, Form("gMeanpTRun3%s", suffix.c_str()));

        // gMPtPion[i] = GetGraph(fPiKp, Form("gPion_MeanpT%s", suffix.c_str()));
        // gMPtKaon[i] = GetGraph(fPiKp, Form("gKaon_MeanpT%s", suffix.c_str()));
        // gMPtProton[i] = GetGraph(fPiKp, Form("gProton_MeanpT%s", suffix.c_str()));

        gMPtKstarRun2[i] = GetGraph(fKstarRun2, Form("gKstar_MeanpT%s", suffix.c_str()));
        gMPtPhiRun2[i] = GetGraph(fPhi, Form("gPhi_MeanpT%s", suffix.c_str()));
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

        gMYieldKstarRun2[i] = GetGraph(fKstarRun2, Form("gKstar_MeanYield%s", suffix.c_str()));
        gMYieldPhiRun2[i] = GetGraph(fPhi, Form("gPhi_MeanYield%s", suffix.c_str()));
        gMYieldChKstarRun2[i] = GetGraph(fChKstar, Form("gChKstar_MeanYield%s", suffix.c_str()));
        gMYieldKshortRun2[i] = GetGraph(fMultStrange, Form("gKshort_MeanYield%s", suffix.c_str()));
        gMYieldLambdaRun2[i] = GetGraph(fMultStrange, Form("gLambda_MeanYield%s", suffix.c_str()));
        gMYieldXiRun2[i] = GetGraph(fMultStrange, Form("gXi_MeanYield%s", suffix.c_str()));
        gMYieldOmegaRun2[i] = GetGraph(fMultStrange, Form("gOmega_MeanYield%s", suffix.c_str()));
        gMYieldXiStarRun2[i] = GetGraph(fXiStarSigma, Form("gXi1530_MeanYield%s", suffix.c_str()));
        gMYieldSigmaRun2[i] = GetGraph(fXiStarSigma, Form("gSigma1385_MeanYield%s", suffix.c_str()));
        gMYieldLambda1520[i] = GetGraph(fLambda1520Run2, Form("gLambda_MeanYield%s", suffix.c_str()));

        // Particle ratios in run2
        gYieldPhiKaRatio[i] = GetGraph(fPhi, Form("gPhi_KaRatio%s", suffix.c_str()));
        gYieldChKstarKshortRatio[i] = GetGraph(fChKstar, Form("gChKstarKaRatio%s", suffix.c_str()));
        gYieldSigmaLambdaRatio[i] = GetGraph(fXiStarSigma, Form("gSigma1385_to_LambdaRatio%s", suffix.c_str()));
        gYieldXiStarXiRatio[i] = GetGraph(fXiStarSigma, Form("gXi1530_to_XiRatio%s", suffix.c_str()));
        gYieldLambda1520LambdaRatio[i] = GetGraph(fLambda1520Run2, Form("gLambdaStar_LambdaRatio%s", suffix.c_str()));

        // Only present for run3 (In run2 the Lambda paper is in IRC review)
        // gMPtLambda1520[i] = GetGraph(fLambda1520, Form("gMeanpTRun3%s", suffix.c_str()));
        // gMYieldLambda1520[i] = GetGraph(fLambda1520, Form("gMeanYieldRun3%s", suffix.c_str()));
        gMPtRho[i] = GetGraph(fRho, Form("gMeanpTRun3%s", suffix.c_str()));
        gMYieldRho[i] = GetGraph(fRho, Form("gMeanYieldRun3%s", suffix.c_str()));
    }

    // ALICE Run2 results
    TFile *fpp13TeV = new TFile("ConversionCodes/pp13TeVALICE.root", "read");
    TFile *fpp7TeV = new TFile("ConversionCodes/pp7TeVALICE.root", "read");
    if (fpp13TeV->IsZombie() || fpp7TeV->IsZombie())
    {
        cout << "Error: ALICE reference files not found" << endl;
        return;
    }

    TGraphErrors *gKstar_MeanpT_13TeV[2], *gKstar_Yield_13TeV[2], *gKstarKshortRatio_13TeV[2], *gPhiKshortRatio_13TeV[2], *gPhiPiRatio_13TeV[2], *gKstarKaRatio_13TeV[2], *gKstarPiRatio_13TeV[2], *gPhiKaRatio_13TeV[2], *gPhiPrRatio_13TeV[2], *gPhiLambdaRatio_13TeV[2];
    TGraphErrors *gKstarPiRatio_7TeV[2], *gPhiPiRatio_7TeV[2];

    for (int i = 0; i < 2; i++)
    {
        gKstar_MeanpT_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_MeanpT%s", (i == 0) ? "_stat" : "_sys"));
        gKstar_Yield_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_MeanYield%s", (i == 0) ? "_stat" : "_sys"));
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

    // /*
    //======================================================
    //    ===========EPOS local model files===========
    //======================================================
    // TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA.root", "read");
    // TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA_ptCut2.root", "read");
    // TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA_NopTCut.root", "read");
    TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA_INELgt0Correct.root", "read"); // latest one
    // TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA_ptCut_FinerBins.root", "read"); //Used for results
    // TFile *fEPOS = new TFile("EPOS_finalQA_CorrectpTCutPiKp.root", "read"); // temporary for checks
    if (fEPOS->IsZombie())
    {
        cout << "Error: EPOS file not found" << endl;
        return;
    }

    //======================
    // Long lived particles
    //======================
    TGraphErrors *gMYieldKshortEPOS_IST0 = GetGraph(fEPOS, "IST0/kshort_vs_mult");
    TGraphErrors *gMYieldProtonEPOS_IST0 = GetGraph(fEPOS, "IST0/proton_vs_mult");
    TGraphErrors *gMYieldPionEPOS_IST0 = GetGraph(fEPOS, "IST0/pion_vs_mult");
    TGraphErrors *gMYieldKaonEPOS_IST0 = GetGraph(fEPOS, "IST0/kaon_vs_mult");
    RestrictModelXaxis(gMYieldKshortEPOS_IST0, 3.1, 23.5);
    RestrictModelXaxis(gMYieldProtonEPOS_IST0, 3.1, 23.5);
    RestrictModelXaxis(gMYieldPionEPOS_IST0, 3.1, 23.5);
    RestrictModelXaxis(gMYieldKaonEPOS_IST0, 3.1, 23.5);

    enum Particle2
    {
        kPhi_epos,
        kKstar_epos,
        kKshort_epos,
        kProton_epos,
        kPion_epos,
        kKaon_epos,
        kKstarPM_epos,
        kXi1530_epos,
        kNParticles_epos
    };

    string ParticleNames[] = {"phi", "kstar", "kshort", "proton", "pion", "kaon", "chargedkstar", "xistar"};

    enum ITYvalues
    {
        kITY0,  // Without ITY condition
        kITY80, // With ITY = 80 condition (UrQMD ON)
        kITY81, // With ITY = 81 condition (Rescattering)
        kNITYvalues
    };
    int ITYvaluesArray[] = {0, 80, 81};
    const int kNIST = 10; // IST0-IST9

    TGraphErrors *gEPOS_Yield[kNIST][kNITYvalues][kNParticles_epos];
    TGraphErrors *gEPOS_MeanPt[kNIST][kNITYvalues][kNParticles_epos];

    for (int ist = 0; ist < kNIST; ist++)
    {
        for (int ity = 0; ity < kNITYvalues; ity++)
        {
            for (int ipart = 0; ipart < kNParticles_epos; ipart++)
            {
                if (ity == 0)
                {
                    gEPOS_Yield[ist][ity][ipart] = GetGraph(fEPOS, Form("IST%d/%s_vs_mult", ist, ParticleNames[ipart].c_str()));
                    gEPOS_MeanPt[ist][ity][ipart] = GetGraph(fEPOS, Form("IST%d/meanpt_%s_vs_mult", ist, ParticleNames[ipart].c_str()));
                    RestrictModelXaxis(gEPOS_Yield[ist][ity][ipart], 3.1, 23.5);
                    RestrictModelXaxis(gEPOS_MeanPt[ist][ity][ipart], 3.1, 23.5);
                }
                else if (ist == 9)
                {
                    gEPOS_Yield[ist][ity][ipart] = GetGraph(fEPOS, Form("IST%d_ITY%d/%s_vs_mult", ist, ITYvaluesArray[ity], ParticleNames[ipart].c_str()));
                    gEPOS_MeanPt[ist][ity][ipart] = GetGraph(fEPOS, Form("IST%d_ITY%d/meanpt_%s_vs_mult", ist, ITYvaluesArray[ity], ParticleNames[ipart].c_str()));
                    RestrictModelXaxis(gEPOS_Yield[ist][ity][ipart], 3.1, 23.5);
                    RestrictModelXaxis(gEPOS_MeanPt[ist][ity][ipart], 3.1, 23.5);
                }
            }
        }
    }

    //===============================================
    // ===========Pythia tested locally================
    //===============================================
    TFile *fPythiaMonash = new TFile("../../pythia/Pythia_MonashLocal.root", "read");
    TFile *fPythiaMonashNoCR = new TFile("../../pythia/Pythia_MonashWoCRLocal.root", "read");
    TFile *fPythiaShoving = new TFile("../../pythia/Pythia_ShovingLocal.root", "read");
    TFile *fPythiaMonashRescattering = new TFile("../../pythia/Pythia_MonashRescatteringLocal.root", "read");
    TFile *fPythiaRopes = new TFile("../../pythia/Pythia_RopesLocal.root", "read");
    TFile *fPythiaRescattering = new TFile("../../pythia/Pythia_RescatteringLocal.root", "read");
    if (fPythiaMonash->IsZombie() || fPythiaMonashNoCR->IsZombie() || fPythiaShoving->IsZombie() || fPythiaMonashRescattering->IsZombie() || fPythiaRopes->IsZombie() || fPythiaRescattering->IsZombie())
    {
        cout << "Error: Pythia local files not found" << endl;
        return;
    }
    // enum PythiaModel
    // {
    //     kPythiaMonashLocal,
    //     kPythiaMonashNoCRLocal,
    //     kPythiaShovingLocal,
    //     kPythiaRopesLocal,
    //     kPythiaMonashRescatteringLocal,
    //     kPythiaRescatteringLocal,
    //     kNPythiaModels
    // };
    enum PythiaModel
    {
        kPythiaMonashLocal,
        kPythiaShovingLocal,
        kNPythiaModels
    };

    // const char *modelLabelLocal[kNPythiaModels] = {
    //     "Pythia Monash",
    //     "Pythia Monash No CR",
    //     "Pythia Shoving",
    //     "Pythia Ropes",
    //     "Pythia Monash Rescattering",
    //     "Pythia Rescattering"};
    const char *modelLabelLocal[kNPythiaModels] = {
        "Pythia Monash",
        "Pythia Shoving"};
    // vector<TFile *> fPythiaModels = {fPythiaMonash, fPythiaMonashNoCR, fPythiaShoving, fPythiaRopes, fPythiaMonashRescattering, fPythiaRescattering};
    vector<TFile *> fPythiaModels = {fPythiaMonash, fPythiaShoving};
    int lineStylesPythia[kNPythiaModels + 10] = {1, 1, 2, 2, 2, 2};
    TGraphErrors *gPythiaYieldLocal[kNPythiaModels][kNParticles_epos];
    TGraphErrors *gPythiaMeanPtLocal[kNPythiaModels][kNParticles_epos];

    int colorsPythia[kNPythiaModels + 10] = {kCyan + 1, kYellow + 1, kMagenta, kGreen + 2, kBlue + 1, kRed + 1};

    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        for (int ipart = 0; ipart < kNParticles_epos; ipart++)
        {
            gPythiaYieldLocal[imodel][ipart] = GetGraph(fPythiaModels[imodel], Form("gMeanYield_%s", ParticleNames[ipart].c_str()));
            gPythiaMeanPtLocal[imodel][ipart] = GetGraph(fPythiaModels[imodel], Form("gMeanpT_%s", ParticleNames[ipart].c_str()));
            RestrictModelXaxis(gPythiaYieldLocal[imodel][ipart], 3.1, 23.5);
            RestrictModelXaxis(gPythiaMeanPtLocal[imodel][ipart], 3.1, 23.5);
            gPythiaYieldLocal[imodel][ipart]->SetLineColor(colorsPythia[imodel]);
            gPythiaMeanPtLocal[imodel][ipart]->SetLineColor(colorsPythia[imodel]);
            gPythiaYieldLocal[imodel][ipart]->SetLineStyle(lineStylesPythia[imodel]);
            gPythiaMeanPtLocal[imodel][ipart]->SetLineStyle(lineStylesPythia[imodel]);
            gPythiaYieldLocal[imodel][ipart]->SetLineWidth(3);
            gPythiaMeanPtLocal[imodel][ipart]->SetLineWidth(3);
        }
    }

    //================================================
    //==============Models from hyperloop outputs======
    //=================================================

    enum Model
    {
        kEPOS_Hydro,
        kPythiaCR,
        kPythiaMonash,
        kPythiaRopes,
        kPythiaShoving,
        kPythiaMonashRescattering,
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
        {kGray + 1, 2},  // EPOS
        {kBlue + 1, 2},  // Pythia CR
        {kRed + 1, 4},   // Pythia Monash
        {kCyan + 1, 7},  // Pythia Ropes
        {kMagenta, 3},   // Pythia Shoving
        {kOrange + 1, 2} // Pythia Monash Rescattering
    };

    const char *modelLabel[kNModels] = {
        "EPOS Hydro",
        "Pythia CR",
        "Pythia Monash",
        "Pythia Ropes",
        "Pythia Shoving",
        "Pythia Monash Rescattering"};

    TFile *ModelsHyperloop = new TFile("ModelRootFiles/ModelResults.root", "read");
    if (ModelsHyperloop->IsZombie())
    {
        cout << "Error: Hyperloop model file not found" << endl;
        return;
    }

    string hyperloopModels[kNModels] = {"EPOS_Hydro", "Pythia_CR", "Pythia_Monash", "Pythia_Ropes", "Pythia_Shoving", "Pythia_Monash_Rescattering"};
    string particles[kNParticles] = {"Kstar", "Phi", "Pion", "Kaon", "Proton", "PionMinus", "KaonMinus", "AntiProton", "Xi1530", "Kshort", "KstarPM"};

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

            RestrictModelXaxis(gMYield[iModel][iParticle], 3.3, 22.5);
            RestrictModelXaxis(gMeanPt[iModel][iParticle], 3.3, 22.5);
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

    ////======================================================
    //    ================Plots=====================
    //======================================================

    // vector<int> modelsToPlot = {
    //     kPythiaMonash,
    //     kPythiaShoving,
    //     kPythiaRopes};
    vector<int> modelsToPlot = {};

    TCanvas *cdNdyKstar = new TCanvas("cdNdyKstar", "cdNdyKstar", 720, 720);
    SetCanvasStyle(cdNdyKstar, 0.15, 0.03, 0.03, 0.15);
    gMYieldKstar[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldKstar[0]->GetYaxis()->SetTitle("dN/dy");
    gMYieldKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldKstar[0]->SetLineWidth(3);
    gMYieldKstar[0]->GetYaxis()->SetRangeUser(0, 0.68);
    gMYieldKstar[0]->SetMarkerColor(kRed + 1);
    gMYieldKstar[0]->SetLineColor(kRed + 1);
    gMYieldKstar[0]->Draw("APE");
    gMYieldKstar[1]->SetFillStyle(0);
    gMYieldKstar[1]->SetLineColor(kRed + 1);
    gMYieldKstar[1]->SetLineWidth(3);
    gMYieldKstar[1]->Draw("5 same");
    gKstar_Yield_13TeV[0]->SetMarkerStyle(21);
    gKstar_Yield_13TeV[0]->SetMarkerColor(kBlue + 1);
    gKstar_Yield_13TeV[0]->SetLineColor(kBlue + 1);
    gKstar_Yield_13TeV[0]->SetLineWidth(3);
    gKstar_Yield_13TeV[0]->Draw("P same");
    gKstar_Yield_13TeV[1]->SetLineColor(kBlue + 1);
    gKstar_Yield_13TeV[1]->SetFillStyle(0);
    gKstar_Yield_13TeV[1]->SetLineWidth(3);
    gKstar_Yield_13TeV[1]->Draw("5 same");

    // cout << "\n==============================================================\n";
    // cout << " Mult      Run3 syst(%)      Run2 syst(%)\n";
    // cout << "==============================================================\n";
    // double maxRun3 = -1;
    // double maxRun2 = -1;

    // int nRun3 = gMYieldKstar[1]->GetN();
    // int nRun2 = gKstar_Yield_13TeV[1]->GetN();

    // int nPoints = TMath::Min(nRun3, nRun2);

    // for (int i = 0; i < nPoints; i++)
    // {
    //     double x3, y3;
    //     double x2, y2;

    //     gMYieldKstar[1]->GetPoint(i, x3, y3);
    //     gKstar_Yield_13TeV[1]->GetPoint(i, x2, y2);

    //     double err3 = gMYieldKstar[1]->GetErrorY(i);
    //     double err2 = gKstar_Yield_13TeV[1]->GetErrorY(i);

    //     double perc3 = 100.0 * err3 / y3;
    //     double perc2 = 100.0 * err2 / y2;

    //     if (perc3 > maxRun3)
    //         maxRun3 = perc3;

    //     if (perc2 > maxRun2)
    //         maxRun2 = perc2;

    //     printf("%8.2f    %8.2f%%      %8.2f%%\n",
    //            x3, perc3, perc2);
    // }

    // cout << "==============================================================\n";
    // cout << Form("Maximum Run3 systematic uncertainty = %.2f %%\n", maxRun3);
    // cout << Form("Maximum Run2 systematic uncertainty = %.2f %%\n", maxRun2);
    // cout << "==============================================================\n";

    ScaleGraph(gEPOS_Yield[9][kITY0][kKstar_epos], 0.5); // In EPOS it is sum not average for all particles in latest root file.
    gEPOS_Yield[9][kITY0][kKstar_epos]->SetLineStyle(9);
    gEPOS_Yield[9][kITY0][kKstar_epos]->SetLineWidth(3);
    gEPOS_Yield[9][kITY0][kKstar_epos]->Draw("l same");
    ScaleGraph(gEPOS_Yield[9][kITY80][kKstar_epos], 0.5);
    gEPOS_Yield[9][kITY80][kKstar_epos]->SetLineColor(kBlue - 6);
    gEPOS_Yield[9][kITY80][kKstar_epos]->SetLineStyle(1);
    gEPOS_Yield[9][kITY80][kKstar_epos]->SetLineWidth(3);
    gEPOS_Yield[9][kITY80][kKstar_epos]->Draw("l same");
    // gEPOS_Yield[9][kITY81][kKstar_epos]->SetLineColor(kGreen + 2);
    // gEPOS_Yield[9][kITY81][kKstar_epos]->SetLineStyle(4);
    // gEPOS_Yield[9][kITY81][kKstar_epos]->Draw("l same");
    // gEPOS_Yield[6][kITY0][kKstar_epos]->SetLineColor(kMagenta);
    // gEPOS_Yield[6][kITY0][kKstar_epos]->SetLineStyle(2);
    // gEPOS_Yield[6][kITY0][kKstar_epos]->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMYield[model][kKstar]->Draw("l same");
    }

    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        ScaleGraph(gPythiaYieldLocal[imodel][kKstar_epos], 0.5);
        gPythiaYieldLocal[imodel][kKstar_epos]->Draw("l same");
    }

    TLegend *legend = new TLegend(0.18, 0.75, 0.48, 0.9);
    SetLegendStyle(legend);
    legend->SetTextSize(0.027);
    // legend->SetHeader("K* (892)^{0}");
    legend->AddEntry(gMYieldKstar[0], "pp, #sqrt{s} = 13.6 TeV", "P");
    legend->AddEntry(gKstar_Yield_13TeV[0], "pp, #sqrt{s} = 13 TeV", "P");
    legend->Draw();

    TLegend *legend2 = new TLegend(0.5, 0.72, 0.8, 0.92);
    SetLegendStyle(legend2);
    legend2->SetTextSize(0.027);
    legend2->AddEntry(gEPOS_Yield[9][kITY0][kKstar_epos], "EPOS UrQMD OFF", "L");
    legend2->AddEntry(gEPOS_Yield[9][kITY80][kKstar_epos], "EPOS UrQMD ON", "L");
    // legend2->AddEntry(gEPOS_Yield[9][kITY81][kKstar_epos], "EPOS Rescattering", "L");

    for (auto model : modelsToPlot)
    {
        legend2->AddEntry(gMYield[model][kKstar], modelLabel[model], "L");
    }
    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        legend2->AddEntry(gPythiaYieldLocal[imodel][kKstar], modelLabelLocal[imodel], "L");
    }
    legend->Draw();
    legend2->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.027);
    // latex.DrawLatex(0.7, 0.85, "#frac{K* (892)^{0} + #bar{K}* (892)^{0}}{2}");
    latex.DrawLatex(0.28, 0.9, "K* (892)^{0}");
    // if (isSavePlots)
    {
        cdNdyKstar->SaveAs("Plots/MeanYield_Kstar_EPOS_UrQMDON.png");
    }

    //===================================================
    //  ================<pT> K*0======================
    //===================================================
    TCanvas *cMeanPtKstar = new TCanvas("cMeanPtKstar", "cMeanPtKstar", 720, 720);
    SetCanvasStyle(cMeanPtKstar, 0.15, 0.03, 0.03, 0.15);
    SetGraphErrorStyle(gMPtKstar[0]);
    gMPtKstar[0]->SetTitle(0);
    gMPtKstar[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMPtKstar[0]->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/#it{c})");
    gMPtKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMPtKstar[0]->GetYaxis()->SetRangeUser(0.32, 1.96);
    gMPtKstar[0]->SetMarkerColor(kRed);
    gMPtKstar[0]->SetLineColor(kRed);
    gMPtKstar[0]->SetLineWidth(3);
    gMPtKstar[0]->Draw("APE");
    gMPtKstar[1]->SetFillStyle(0);
    gMPtKstar[1]->SetLineColor(kRed);
    gMPtKstar[1]->SetLineWidth(3);
    gMPtKstar[1]->Draw("5 same");
    gKstar_MeanpT_13TeV[0]->SetMarkerStyle(21);
    gKstar_MeanpT_13TeV[0]->SetMarkerColor(kBlue);
    gKstar_MeanpT_13TeV[0]->SetLineColor(kBlue);
    gKstar_MeanpT_13TeV[0]->SetLineWidth(3);
    gKstar_MeanpT_13TeV[0]->Draw("P same");
    gKstar_MeanpT_13TeV[1]->SetLineColor(kBlue);
    gKstar_MeanpT_13TeV[1]->SetFillStyle(0);
    gKstar_MeanpT_13TeV[1]->SetLineWidth(3);
    gKstar_MeanpT_13TeV[1]->Draw("5 same");

    gEPOS_MeanPt[9][kITY0][kKstar_epos]->SetLineStyle(6);
    gEPOS_MeanPt[9][kITY0][kKstar_epos]->Draw("C same");
    gEPOS_MeanPt[9][kITY80][kKstar_epos]->SetLineColor(kBlue);
    gEPOS_MeanPt[9][kITY80][kKstar_epos]->SetLineStyle(2);
    gEPOS_MeanPt[9][kITY80][kKstar_epos]->Draw("C same");
    // gEPOS_MeanPt[9][kITY81][kKstar_epos]->SetLineColor(kGreen + 2);
    // gEPOS_MeanPt[9][kITY81][kKstar_epos]->SetLineStyle(4);
    // gEPOS_MeanPt[9][kITY81][kKstar_epos]->Draw("l same");
    // gEPOS_MeanPt[6][kITY0][kKstar_epos]->SetLineColor(kMagenta);
    // gEPOS_MeanPt[6][kITY0][kKstar_epos]->SetLineStyle(2);
    // gEPOS_MeanPt[6][kITY0][kKstar_epos]->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMeanPt[model][kKstar]->Draw("l same");
    }
    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        gPythiaMeanPtLocal[imodel][kKstar_epos]->Draw("l same");
    }
    legend->Draw();
    legend2->Draw();
    // latex.DrawLatex(0.7, 0.85, "#frac{K* (892)^{0} + #bar{K}* (892)^{0}}{2}");
    latex.DrawLatex(0.28, 0.9, "K* (892)^{0}");
    // if (isSavePlots)
    {
        cMeanPtKstar->SaveAs("Plots/MeanPt_Kstar_EPOS_UrQMDON.png");
    }

    //===================================================
    //  ================<pT> all particles======================
    //===================================================
    TCanvas *cMeanPtAll = new TCanvas("cMeanPtAll", "cMeanPtAll", 720, 720);
    SetCanvasStyle(cMeanPtAll, 0.15, 0.03, 0.03, 0.15);
    gMPtKstar[0]->SetMaximum(2.08);
    gMPtKstar[0]->SetMinimum(0.27);
    gMPtKstar[0]->Draw("APE");
    gMPtKstar[1]->Draw("5 same");
    gMPtPion[0]->SetMarkerColor(kGreen + 2);
    gMPtPion[0]->SetLineColor(kGreen + 2);
    gMPtPion[0]->Draw("PE same");
    gMPtPion[1]->SetLineColor(kGreen + 2);
    gMPtPion[1]->SetFillStyle(0);
    gMPtPion[1]->Draw("5 same");
    gMPtKaon[0]->SetMarkerColor(kBlue);
    gMPtKaon[0]->SetLineColor(kBlue);
    gMPtKaon[0]->Draw("PE same");
    gMPtKaon[1]->SetLineColor(kBlue);
    gMPtKaon[1]->SetFillStyle(0);
    gMPtKaon[1]->Draw("5 same");
    gMPtProton[0]->SetMarkerColor(kMagenta);
    gMPtProton[0]->SetLineColor(kMagenta);
    gMPtProton[0]->Draw("PE same");
    gMPtProton[1]->SetLineColor(kMagenta);
    gMPtProton[1]->SetFillStyle(0);
    gMPtProton[1]->Draw("5 same");
    gMPtKshortRun2[0]->SetMarkerColor(kCyan + 1);
    gMPtKshortRun2[0]->SetLineColor(kCyan + 1);
    gMPtKshortRun2[0]->Draw("PE same");
    gMPtKshortRun2[1]->SetLineColor(kCyan + 1);
    gMPtKshortRun2[1]->SetFillStyle(0);
    gMPtKshortRun2[1]->Draw("5 same");
    gMPtChKstarRun2[0]->SetMarkerColor(kOrange + 1);
    gMPtChKstarRun2[0]->SetLineColor(kOrange + 1);
    gMPtChKstarRun2[0]->Draw("PE same");
    gMPtChKstarRun2[1]->SetLineColor(kOrange + 1);
    gMPtChKstarRun2[1]->SetFillStyle(0);
    gMPtChKstarRun2[1]->Draw("5 same");
    gMPtXiStarRun2[0]->SetMarkerColor(kViolet + 2);
    gMPtXiStarRun2[0]->SetLineColor(kViolet + 2);
    gMPtXiStarRun2[0]->Draw("PE same");
    gMPtXiStarRun2[1]->SetLineColor(kViolet + 2);
    gMPtXiStarRun2[1]->SetFillStyle(0);
    gMPtXiStarRun2[1]->Draw("5 same");
    gMPtPhiRun2[0]->SetMarkerColor(kGreen - 1);
    gMPtPhiRun2[0]->SetLineColor(kGreen - 1);
    gMPtPhiRun2[0]->Draw("PE same");
    gMPtPhiRun2[1]->SetLineColor(kGreen - 1);
    gMPtPhiRun2[1]->SetFillStyle(0);
    gMPtPhiRun2[1]->Draw("5 same");
    gMPtLambdaRun2[0]->SetMarkerColor(kBrown);
    gMPtLambdaRun2[0]->SetLineColor(kBrown);
    gMPtLambdaRun2[0]->Draw("PE same");
    gMPtLambdaRun2[1]->SetLineColor(kBrown);
    gMPtLambdaRun2[1]->SetFillStyle(0);
    gMPtLambdaRun2[1]->Draw("5 same");
    gMPtXiRun2[0]->SetMarkerColor(kGray + 2);
    gMPtXiRun2[0]->SetLineColor(kGray + 2);
    gMPtXiRun2[0]->Draw("PE same");
    gMPtXiRun2[1]->SetLineColor(kGray + 2);
    gMPtXiRun2[1]->SetFillStyle(0);
    gMPtXiRun2[1]->Draw("5 same");
    gMPtOmegaRun2[0]->SetMarkerColor(kPink + 1);
    gMPtOmegaRun2[0]->SetLineColor(kPink + 1);
    gMPtOmegaRun2[0]->Draw("PE same");
    gMPtOmegaRun2[1]->SetLineColor(kPink + 1);
    gMPtOmegaRun2[1]->SetFillStyle(0);
    gMPtOmegaRun2[1]->Draw("5 same");
    gMPtSigmaRun2[0]->SetMarkerColor(kAzure + 7);
    gMPtSigmaRun2[0]->SetLineColor(kAzure + 7);
    gMPtSigmaRun2[0]->Draw("PE same");
    gMPtSigmaRun2[1]->SetLineColor(kAzure + 7);
    gMPtSigmaRun2[1]->SetFillStyle(0);
    gMPtSigmaRun2[1]->Draw("5 same");
    gMPtLambda1520[0]->SetMarkerColor(43);
    gMPtLambda1520[0]->SetLineColor(43);
    gMPtLambda1520[0]->Draw("PE same");
    gMPtLambda1520[1]->SetLineColor(43);
    gMPtLambda1520[1]->SetFillStyle(0);
    gMPtLambda1520[1]->Draw("5 same");

    // // Run3 particles
    // gMPtRho[0]->SetMarkerColor(kGray + 2);
    // gMPtRho[0]->SetLineColor(kGray + 2);
    // gMPtRho[0]->Draw("PE same");
    // gMPtRho[1]->SetLineColor(kGray + 2);
    // gMPtRho[1]->SetFillStyle(0);
    // gMPtRho[1]->Draw("5 same");

    TLegend *legendMeanPt = new TLegend(0.2, 0.78, 0.7, 0.95);
    SetLegendStyle(legendMeanPt);
    legendMeanPt->SetTextSize(0.027);
    legendMeanPt->SetNColumns(4);
    legendMeanPt->AddEntry(gMPtPion[0], "Pion", "P");
    legendMeanPt->AddEntry(gMPtKaon[0], "Kaon", "P");
    legendMeanPt->AddEntry(gMPtKshortRun2[0], "K_{S}^{0}", "P");
    legendMeanPt->AddEntry(gMPtChKstarRun2[0], "K*^{#pm}", "P");
    legendMeanPt->AddEntry(gMPtKstar[0], "K* (892)^{0}", "P");
    legendMeanPt->AddEntry(gMPtPhiRun2[0], "#phi (1020)", "P");
    legendMeanPt->AddEntry(gMPtProton[0], "Proton", "P");
    legendMeanPt->AddEntry(gMPtSigmaRun2[0], "#Sigma", "P");
    legendMeanPt->AddEntry(gMPtLambdaRun2[0], "#Lambda", "P");
    legendMeanPt->AddEntry(gMPtXiRun2[0], "#Xi", "P");
    legendMeanPt->AddEntry(gMPtLambda1520[0], "#Lambda(1520)", "P");
    legendMeanPt->AddEntry(gMPtXiStarRun2[0], "#Xi(1530)", "P");
    legendMeanPt->AddEntry(gMPtOmegaRun2[0], "#Omega", "P");

    // legendMeanPt->AddEntry(gMPtRho[0], "#rho (770)", "P");
    legendMeanPt->Draw();
    if (isSavePlots)
    {
        cMeanPtAll->SaveAs("Plots/MeanPt_AllParticles_Run3.png");
    }

    //===================================================
    //  ================<pT> mesons and baryons======================
    //===================================================
    TCanvas *cMeanPtMesonsBaryons = new TCanvas("cMeanPtMesonsBaryons", "cMeanPtMesonsBaryons", 720, 720);
    SetCanvasStyle(cMeanPtMesonsBaryons, 0.15, 0.03, 0.03, 0.15);
    double pad1Size = 0.5, pad2Size = 0.5;
    canvas_style(cMeanPtMesonsBaryons, pad1Size, pad2Size, 0.13, 0.02, 0.02, 0.25);
    cMeanPtMesonsBaryons->cd(1);
    TGraphErrors *gMPtKstarClone = (TGraphErrors *)gMPtKstar[0]->Clone("gMPtKstarClone");
    TGraphErrors *gMPtKstarClone_sys = (TGraphErrors *)gMPtKstar[1]->Clone("gMPtKstarClone_sys");
    gMPtKstarClone->SetMaximum(1.79);
    gMPtKstarClone->SetMinimum(0.31);
    gMPtKstarClone->GetXaxis()->SetLimits(0, 27);
    gMPtKstarClone->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    gMPtKstarClone->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    gMPtKstarClone->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    gMPtKstarClone->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    gMPtKstarClone->GetYaxis()->SetTitleOffset(1.5 * pad1Size);
    gMPtKstarClone->Draw("APE");
    gMPtKstarClone_sys->Draw("5 same");
    gMPtPion[0]->Draw("PE same");
    gMPtPion[1]->Draw("5 same");
    gMPtKaon[0]->Draw("PE same");
    gMPtKaon[1]->Draw("5 same");
    gMPtKshortRun2[0]->Draw("PE same");
    gMPtKshortRun2[1]->Draw("5 same");
    gMPtChKstarRun2[0]->Draw("PE same");
    gMPtChKstarRun2[1]->Draw("5 same");
    gMPtPhiRun2[0]->Draw("PE same");
    gMPtPhiRun2[1]->Draw("5 same");

    TLegend *legendMeanPtMesons = new TLegend(0.17, 0.79, 0.65, 0.94);
    SetLegendStyle(legendMeanPtMesons);
    legendMeanPtMesons->SetNColumns(3);
    legendMeanPtMesons->SetTextSize(0.027 / pad1Size);
    legendMeanPtMesons->AddEntry(gMPtPion[0], "Pion", "P");
    legendMeanPtMesons->AddEntry(gMPtKaon[0], "Kaon", "P");
    legendMeanPtMesons->AddEntry(gMPtKshortRun2[0], "K_{S}^{0}", "P");
    legendMeanPtMesons->AddEntry(gMPtChKstarRun2[0], "K*^{#pm}", "P");
    legendMeanPtMesons->AddEntry(gMPtKstar[0], "K*^{0} (run3)", "P");
    legendMeanPtMesons->AddEntry(gMPtPhiRun2[0], "#phi(1020)", "P");
    legendMeanPtMesons->Draw();
    TLatex latexMesons;
    latexMesons.SetNDC();
    // latexMesons.SetTextFont(42);
    latexMesons.SetTextSize(0.06);
    latexMesons.DrawLatex(0.82, 0.07, "Mesons");

    cMeanPtMesonsBaryons->cd(2);
    gMPtProton[0]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    gMPtProton[0]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    gMPtProton[0]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    gMPtProton[0]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    gMPtProton[0]->GetYaxis()->SetTitleOffset(1.5 * pad2Size);
    gMPtProton[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMPtProton[0]->GetXaxis()->SetLimits(0, 27);
    gMPtProton[0]->SetMaximum(2.07);
    gMPtProton[0]->SetMinimum(0.63);
    gMPtProton[0]->Draw("APE");
    gMPtProton[1]->Draw("5 same");
    gMPtXiStarRun2[0]->Draw("PE same");
    gMPtXiStarRun2[1]->Draw("5 same");
    gMPtLambdaRun2[0]->Draw("PE same");
    gMPtLambdaRun2[1]->Draw("5 same");
    gMPtXiRun2[0]->Draw("PE same");
    gMPtXiRun2[1]->Draw("5 same");
    gMPtOmegaRun2[0]->Draw("PE same");
    gMPtOmegaRun2[1]->Draw("5 same");
    gMPtSigmaRun2[0]->Draw("PE same");
    gMPtSigmaRun2[1]->Draw("5 same");
    gMPtLambda1520[0]->Draw("PE same");
    gMPtLambda1520[1]->Draw("5 same");

    TLegend *legendMeanPtBaryons = new TLegend(0.17, 0.76, 0.55, 0.96);
    SetLegendStyle(legendMeanPtBaryons);
    legendMeanPtBaryons->SetNColumns(3);
    legendMeanPtBaryons->SetTextSize(0.027 / pad2Size);
    legendMeanPtBaryons->AddEntry(gMPtProton[0], "Proton", "P");
    legendMeanPtBaryons->AddEntry(gMPtLambdaRun2[0], "#Lambda", "P");
    legendMeanPtBaryons->AddEntry(gMPtSigmaRun2[0], "#Sigma(1385)^{#pm}", "P");
    legendMeanPtBaryons->AddEntry(gMPtXiRun2[0], "#Xi^{0}", "P");
    legendMeanPtBaryons->AddEntry(gMPtLambda1520[0], "#Lambda(1520)", "P");
    legendMeanPtBaryons->AddEntry(gMPtXiStarRun2[0], "#Xi*(1530)^{0}", "P");
    legendMeanPtBaryons->AddEntry(gMPtOmegaRun2[0], "#Omega", "P");
    legendMeanPtBaryons->Draw();
    latexMesons.DrawLatex(0.82, 0.3, "Baryons");
    if (isSavePlots)
        cMeanPtMesonsBaryons->SaveAs("Plots/MeanPt_Mesons_Baryons_Run3.png");

    //===================================================
    //  ========= <pT>/particle mass ======================
    //===================================================
    double PDGMasses[] = {
        0.8956,  // 0. K*0
        1.0195,  // 1. phi
        0.1396,  // 2. pion
        0.494,   // 3. Kaon
        0.938,   // 4. proton
        0.498,   // 5. Kshort
        0.8955,  // 6. charged K*
        1.5318,  // 7. Xi(1530)
        1.5192,  // 8. Lambda(1520)
        0.7753,  // 9. Rho(770)
        1.31486, // 10. Xi(Ground state)
        1.38283, // 11. Sigma(1385)
        1.67245, // 12. Omega
        1.115683 // 13. Lambda (ground state)
    };
    double PDGMassesErrors[] = {0.0002, 0.000016, 0.00000018, 0.000015, 29e-9, 0.000013, 0.0002, 0.00034, 0.00019, 0.0002, 0.0002, 0.00034, 0.00032, 0.000006};
    int nq[] = {2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3};
    string ParticleNames2[] = {"K*0", "phi", "pion", "Kaon", "proton", "Kshort", "charged K*", "Xi(1530)", "Lambda(1520)", "Rho(770)", "Xi(Ground state)", "Sigma(1385)", "Omega", "Lambda (ground state)"};

    vector<TGraphErrors **>
        gMeanPtAll = {gMPtKstar, gMPtPhiRun2, gMPtPion, gMPtKaon, gMPtProton, gMPtKshortRun2, gMPtChKstarRun2, gMPtXiStarRun2, gMPtLambda1520, gMPtRho, gMPtXiRun2, gMPtSigmaRun2, gMPtOmegaRun2, gMPtLambdaRun2};

    TCanvas *cMeanPtMass = new TCanvas("cMeanPtMass", "cMeanPtMass", 720, 720);
    SetCanvasStyle(cMeanPtMass, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors **gMPtKstarClone2 = new TGraphErrors *[2];
    TGraphErrors **gMPtPionClone = new TGraphErrors *[2];
    TGraphErrors **gMPtKaonClone = new TGraphErrors *[2];
    TGraphErrors **gMPtProtonClone = new TGraphErrors *[2];
    TGraphErrors **gMPtKshortRun2Clone = new TGraphErrors *[2];
    TGraphErrors **gMPtChKstarRun2Clone = new TGraphErrors *[2];
    TGraphErrors **gMPtPhiRun2Clone = new TGraphErrors *[2];
    TGraphErrors **gMPtLambdaRun2Clone = new TGraphErrors *[2];
    TGraphErrors **gMPtXiRun2Clone = new TGraphErrors *[2];
    TGraphErrors **gMPtXiStarRun2Clone = new TGraphErrors *[2];
    TGraphErrors **gMPtSigmaRun2Clone = new TGraphErrors *[2];
    TGraphErrors **gMPtOmegaRun2Clone = new TGraphErrors *[2];
    TGraphErrors **gMPtLambda1520Clone = new TGraphErrors *[2];
    for (int i = 0; i < 2; i++)
    {

        gMPtKstarClone2[i] = (TGraphErrors *)gMPtKstar[i]->Clone(Form("gMPtKstarClone_%d", i));
        gMPtPionClone[i] = (TGraphErrors *)gMPtPion[i]->Clone(Form("gMPtPionClone_%d", i));
        gMPtKaonClone[i] = (TGraphErrors *)gMPtKaon[i]->Clone(Form("gMPtKaonClone_%d", i));
        gMPtProtonClone[i] = (TGraphErrors *)gMPtProton[i]->Clone(Form("gMPtProtonClone_%d", i));
        gMPtKshortRun2Clone[i] = (TGraphErrors *)gMPtKshortRun2[i]->Clone(Form("gMPtKshortRun2Clone_%d", i));
        gMPtChKstarRun2Clone[i] = (TGraphErrors *)gMPtChKstarRun2[i]->Clone(Form("gMPtChKstarRun2Clone_%d", i));
        gMPtPhiRun2Clone[i] = (TGraphErrors *)gMPtPhiRun2[i]->Clone(Form("gMPtPhiRun2Clone_%d", i));
        gMPtLambdaRun2Clone[i] = (TGraphErrors *)gMPtLambdaRun2[i]->Clone(Form("gMPtLambdaRun2Clone_%d", i));
        gMPtXiRun2Clone[i] = (TGraphErrors *)gMPtXiRun2[i]->Clone(Form("gMPtXiRun2Clone_%d", i));
        gMPtXiStarRun2Clone[i] = (TGraphErrors *)gMPtXiStarRun2[i]->Clone(Form("gMPtXiStarRun2Clone_%d", i));
        gMPtSigmaRun2Clone[i] = (TGraphErrors *)gMPtSigmaRun2[i]->Clone(Form("gMPtSigmaRun2Clone_%d", i));
        gMPtOmegaRun2Clone[i] = (TGraphErrors *)gMPtOmegaRun2[i]->Clone(Form("gMPtOmegaRun2Clone_%d", i));
        gMPtLambda1520Clone[i] = (TGraphErrors *)gMPtLambda1520[i]->Clone(Form("gMPtLambda1520Clone_%d", i));

        ScaleGraph(gMPtKstarClone2[i], 1.0 / PDGMasses[0]);
        ScaleGraph(gMPtPionClone[i], 0.55 / PDGMasses[2]);
        ScaleGraph(gMPtKaonClone[i], 1.0 / PDGMasses[3]);
        ScaleGraph(gMPtProtonClone[i], 1.0 / PDGMasses[4]);
        ScaleGraph(gMPtKshortRun2Clone[i], 1.0 / PDGMasses[5]);
        ScaleGraph(gMPtChKstarRun2Clone[i], 1.0 / PDGMasses[6]);
        ScaleGraph(gMPtPhiRun2Clone[i], 1.0 / PDGMasses[1]);
        ScaleGraph(gMPtLambdaRun2Clone[i], 1.0 / PDGMasses[13]);
        ScaleGraph(gMPtXiRun2Clone[i], 1.0 / PDGMasses[10]);
        ScaleGraph(gMPtXiStarRun2Clone[i], 1.0 / PDGMasses[7]);
        ScaleGraph(gMPtSigmaRun2Clone[i], 1.0 / PDGMasses[11]);
        ScaleGraph(gMPtOmegaRun2Clone[i], 1.0 / PDGMasses[12]);
        ScaleGraph(gMPtLambda1520Clone[i], 1.0 / PDGMasses[8]);
    }
    gMPtKstarClone2[0]->SetMaximum(2.69);
    gMPtKstarClone2[0]->SetMinimum(0.52);
    gMPtKstarClone2[0]->GetYaxis()->SetTitle("<#it{p}_{T}>/m (GeV/#it{c}^{2})");
    gMPtKstarClone2[0]->Draw("APE");
    gMPtPionClone[0]->Draw("PE same");
    gMPtKaonClone[0]->Draw("PE same");
    gMPtProtonClone[0]->Draw("PE same");
    gMPtKshortRun2Clone[0]->Draw("PE same");
    gMPtChKstarRun2Clone[0]->Draw("PE same");
    gMPtPhiRun2Clone[0]->Draw("PE same");
    gMPtLambdaRun2Clone[0]->Draw("PE same");
    gMPtXiRun2Clone[0]->Draw("PE same");
    gMPtLambda1520Clone[0]->Draw("PE same");
    gMPtXiStarRun2Clone[0]->Draw("PE same");
    gMPtSigmaRun2Clone[0]->Draw("PE same");
    gMPtOmegaRun2Clone[0]->Draw("PE same");

    gMPtKstarClone2[1]->Draw("5 same");
    gMPtPionClone[1]->Draw("5 same");
    gMPtKaonClone[1]->Draw("5 same");
    gMPtProtonClone[1]->Draw("5 same");
    gMPtKshortRun2Clone[1]->Draw("5 same");
    gMPtChKstarRun2Clone[1]->Draw("5 same");
    gMPtPhiRun2Clone[1]->Draw("5 same");
    gMPtLambdaRun2Clone[1]->Draw("5 same");
    gMPtXiRun2Clone[1]->Draw("5 same");
    gMPtLambda1520Clone[1]->Draw("5 same");
    gMPtXiStarRun2Clone[1]->Draw("5 same");
    gMPtSigmaRun2Clone[1]->Draw("5 same");
    gMPtOmegaRun2Clone[1]->Draw("5 same");

    TLegend *legendMeanPtMass = new TLegend(0.18, 0.8, 0.69, 0.94);
    SetLegendStyle(legendMeanPtMass);
    legendMeanPtMass->SetNColumns(4);
    legendMeanPtMass->SetTextSize(0.027);
    legendMeanPtMass->AddEntry(gMPtPionClone[0], "Pion #times 0.55", "P");
    legendMeanPtMass->AddEntry(gMPtKaonClone[0], "Kaon", "P");
    legendMeanPtMass->AddEntry(gMPtKshortRun2Clone[0], "K_{S}^{0}", "P");
    legendMeanPtMass->AddEntry(gMPtChKstarRun2Clone[0], "K*^{#pm}", "P");
    legendMeanPtMass->AddEntry(gMPtKstarClone2[0], "K* (892)^{0}", "P");
    legendMeanPtMass->AddEntry(gMPtPhiRun2Clone[0], "#phi(1020)", "P");
    legendMeanPtMass->AddEntry(gMPtProtonClone[0], "Proton", "P");
    legendMeanPtMass->AddEntry(gMPtSigmaRun2Clone[0], "#Sigma", "P");
    legendMeanPtMass->AddEntry(gMPtLambdaRun2Clone[0], "#Lambda", "P");
    legendMeanPtMass->AddEntry(gMPtXiRun2Clone[0], "#Xi", "P");
    legendMeanPtMass->AddEntry(gMPtLambda1520Clone[0], "#Lambda(1520)", "P");
    legendMeanPtMass->AddEntry(gMPtXiStarRun2Clone[0], "#Xi(1530)", "P");
    legendMeanPtMass->AddEntry(gMPtOmegaRun2Clone[0], "#Omega", "P");
    legendMeanPtMass->Draw();

    if (isSavePlots)
    {
        cMeanPtMass->SaveAs("Plots/MeanPt_MassScaled_Run3.png");
    }

    // //=====================================================================
    // //  ====== <pT> vs mass/n_q for for different mult classes ===========
    // //=====================================================================
    TCanvas *cMeanPtMassNq = new TCanvas("cMeanPtMassNq", "cMeanPtMassNq", 1440, 720);
    SetCanvasStyle(cMeanPtMassNq, 0.15, 0.03, 0.03, 0.15);
    cMeanPtMassNq->Divide(4, 3);
    TGraphErrors *gMPtMultClasses[10][2]; // 10 multiplicity classes, 2 for statistical and systematic errors
    int multClasses[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    TLatex latexMPtMassNq;
    latexMPtMassNq.SetNDC();
    latexMPtMassNq.SetTextFont(42);
    latexMPtMassNq.SetTextSize(0.04);
    int sizePDGMasses = sizeof(PDGMasses) / sizeof(PDGMasses[0]);

    // Print the dNch/deta values for Xi(1530) and Charged K* for each multiplicity class
    int totalPointsKstar = gMeanPtAll[0][0]->GetN();
    int totalPointsXiStar = gMeanPtAll[7][0]->GetN();
    int totalPointsChKstar = gMeanPtAll[6][0]->GetN();
    int totalPointsPhi = gMeanPtAll[1][0]->GetN();

    for (int iKstar = 0; iKstar < totalPointsKstar; iKstar++)
    {
        double xKstar, yKstar;
        gMeanPtAll[0][0]->GetPoint(iKstar, xKstar, yKstar);
        cout << "Multiplicity class " << iKstar << ": dNch/deta for K*0 = " << xKstar << endl;
    }
    cout << endl;
    for (int iXi = 0; iXi < totalPointsXiStar; iXi++)
    {
        double xXi, yXi;
        gMeanPtAll[7][0]->GetPoint(iXi, xXi, yXi);
        cout << "Multiplicity class " << iXi << ": dNch/deta for Xi(1530) = " << xXi << endl;
    }
    cout << endl;
    for (int iChKstar = 0; iChKstar < totalPointsChKstar; iChKstar++)
    {
        double xChKstar, yChKstar;
        gMeanPtAll[6][0]->GetPoint(iChKstar, xChKstar, yChKstar);
        cout << "Multiplicity class " << iChKstar << ": dNch/deta for Charged K* = " << xChKstar << endl;
    }
    cout << endl;
    for (int iPhi = 0; iPhi < totalPointsPhi; iPhi++)
    {
        double xPhi, yPhi;
        gMeanPtAll[1][0]->GetPoint(iPhi, xPhi, yPhi);
        cout << "Multiplicity class " << iPhi << ": dNch/deta for Phi = " << xPhi << endl;
    }


    for (int iMult = 0; iMult < 10; iMult++)
    {
        gMPtMultClasses[iMult][0] = new TGraphErrors();
        gMPtMultClasses[iMult][1] = new TGraphErrors();
        gMPtMultClasses[iMult][0]->SetName(Form("gMPtMultClass_%d_stat", iMult));
        gMPtMultClasses[iMult][1]->SetName(Form("gMPtMultClass_%d_sys", iMult));

        for (int iMasses = 0; iMasses < sizePDGMasses; iMasses++)
        {
            if (iMult == 0)
                cout << "Particle " << ParticleNames2[iMasses] << ", No. of Mult Classes " << gMeanPtAll[iMasses][0]->GetN() << endl;

            double x, y;
            gMeanPtAll[iMasses][0]->GetPoint(iMult, x, y);
            double xerr = gMeanPtAll[iMasses][0]->GetErrorX(iMult);
            double yerr = gMeanPtAll[iMasses][0]->GetErrorY(iMult);
            double yerr_sys = gMeanPtAll[iMasses][1]->GetErrorY(iMult);

            gMPtMultClasses[iMult][0]->SetPoint(iMasses, PDGMasses[iMasses] / nq[iMasses], y);
            gMPtMultClasses[iMult][1]->SetPoint(iMasses, PDGMasses[iMasses] / nq[iMasses], y);
            gMPtMultClasses[iMult][0]->SetPointError(iMasses, 0, yerr);
            gMPtMultClasses[iMult][1]->SetPointError(iMasses, 0, yerr_sys);
        }
        cMeanPtMassNq->cd(iMult + 1);
        gMPtMultClasses[iMult][0]->SetMarkerStyle(20);
        gMPtMultClasses[iMult][0]->SetMarkerColor(kBlue);
        gMPtMultClasses[iMult][0]->SetLineColor(kBlue);
        gMPtMultClasses[iMult][0]->GetXaxis()->SetTitle("Mass/n_{q} (GeV/#it{c}^{2})");
        gMPtMultClasses[iMult][0]->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/#it{c})");
        gMPtMultClasses[iMult][0]->SetMaximum(2.5);
        gMPtMultClasses[iMult][0]->SetMinimum(0.0);
        gMPtMultClasses[iMult][0]->Draw("APE");
        gMPtMultClasses[iMult][1]->SetLineColor(kBlue);
        gMPtMultClasses[iMult][1]->SetFillStyle(0);
        gMPtMultClasses[iMult][1]->Draw("5 same");
        latexMPtMassNq.DrawLatex(0.6, 0.8, Form("Multiplicity: %d-%d", multClasses[iMult], multClasses[iMult + 1]));
    }
    // cout << "Total size of graph is " << gMPtMultClasses[0][0]->GetN() << endl;
    // for (int iMult = 0; iMult < 10; iMult++)
    // {
    //     cout << "Multiplicity: " << multClasses[iMult] << "-" << multClasses[iMult + 1] << endl;
    //     for (int iMasses = 0; iMasses < sizePDGMasses; iMasses++)
    //     {
    //         double x, y;
    //         gMPtMultClasses[iMult][0]->GetPoint(iMasses, x, y);
    //         cout << "Mass/n_q: " << x << ", <pT>: " << y << endl;
    //     }
    // }
    cMeanPtMassNq->SaveAs("Plots/MeanPt_MassNq_Run3.png");

    //=================================================================================
    // =====================(<pT> - m) vs <dNch/deta> (hadrons) =========================
    //=================================================================================
    TCanvas *cMeanPtMassScaled = new TCanvas("cMeanPtMassScaled", "cMeanPtMassScaled", 720, 720);
    SetCanvasStyle(cMeanPtMassScaled, 0.15, 0.03, 0.03, 0.15);

    TGraphErrors **gMPtKstarScaled = AddConstantToGraph(gMPtKstar, -PDGMasses[0]);
    TGraphErrors **gMPtPionScaled = AddConstantToGraph(gMPtPion, -PDGMasses[2]);
    TGraphErrors **gMPtKaonScaled = AddConstantToGraph(gMPtKaon, -PDGMasses[3]);
    TGraphErrors **gMPtProtonScaled = AddConstantToGraph(gMPtProton, -PDGMasses[4]);
    TGraphErrors **gMPtKshortRun2Scaled = AddConstantToGraph(gMPtKshortRun2, -PDGMasses[5]);
    TGraphErrors **gMPtChKstarRun2Scaled = AddConstantToGraph(gMPtChKstarRun2, -PDGMasses[6]);
    TGraphErrors **gMPtPhiRun2Scaled = AddConstantToGraph(gMPtPhiRun2, -PDGMasses[1]);
    TGraphErrors **gMPtLambdaRun2Scaled = AddConstantToGraph(gMPtLambdaRun2, -PDGMasses[13]);
    TGraphErrors **gMPtXiRun2Scaled = AddConstantToGraph(gMPtXiRun2, -PDGMasses[10]);
    TGraphErrors **gMPtLambda1520Scaled = AddConstantToGraph(gMPtLambda1520, -PDGMasses[8]);
    TGraphErrors **gMPtXiStarRun2Scaled = AddConstantToGraph(gMPtXiStarRun2, -PDGMasses[7]);
    TGraphErrors **gMPtSigmaRun2Scaled = AddConstantToGraph(gMPtSigmaRun2, -PDGMasses[11]);
    TGraphErrors **gMPtOmegaRun2Scaled = AddConstantToGraph(gMPtOmegaRun2, -PDGMasses[12]);

    gMPtKstarScaled[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMPtKstarScaled[0]->GetYaxis()->SetTitle("(<#it{p}_{T}> - m) (GeV/#it{c})");
    gMPtKstarScaled[0]->SetMaximum(0.78);
    gMPtKstarScaled[0]->SetMinimum(-0.65);
    gMPtKstarScaled[0]->Draw("APE");
    gMPtPionScaled[0]->Draw("PE same");
    gMPtKaonScaled[0]->Draw("PE same");
    gMPtProtonScaled[0]->Draw("PE same");
    gMPtKshortRun2Scaled[0]->Draw("PE same");
    gMPtChKstarRun2Scaled[0]->Draw("PE same");
    gMPtPhiRun2Scaled[0]->Draw("PE same");
    gMPtLambdaRun2Scaled[0]->Draw("PE same");
    gMPtXiRun2Scaled[0]->Draw("PE same");
    gMPtLambda1520Scaled[0]->Draw("PE same");
    gMPtXiStarRun2Scaled[0]->Draw("PE same");
    gMPtSigmaRun2Scaled[0]->Draw("PE same");
    gMPtOmegaRun2Scaled[0]->Draw("PE same");

    gMPtKstarScaled[1]->Draw("5 same");
    gMPtPionScaled[1]->Draw("5 same");
    gMPtKaonScaled[1]->Draw("5 same");
    gMPtProtonScaled[1]->Draw("5 same");
    gMPtKshortRun2Scaled[1]->Draw("5 same");
    gMPtChKstarRun2Scaled[1]->Draw("5 same");
    gMPtPhiRun2Scaled[1]->Draw("5 same");
    gMPtLambdaRun2Scaled[1]->Draw("5 same");
    gMPtXiRun2Scaled[1]->Draw("5 same");
    gMPtLambda1520Scaled[1]->Draw("5 same");
    gMPtXiStarRun2Scaled[1]->Draw("5 same");
    gMPtSigmaRun2Scaled[1]->Draw("5 same");
    gMPtOmegaRun2Scaled[1]->Draw("5 same");

    TLegend *legendMeanPtMassScaled = new TLegend(0.18, 0.77, 0.72, 0.95);
    SetLegendStyle(legendMeanPtMassScaled);
    legendMeanPtMassScaled->SetNColumns(4);
    legendMeanPtMassScaled->SetTextSize(0.027);
    legendMeanPtMassScaled->AddEntry(gMPtPionScaled[0], "Pion", "P");
    legendMeanPtMassScaled->AddEntry(gMPtKaonScaled[0], "Kaon", "P");
    legendMeanPtMassScaled->AddEntry(gMPtKshortRun2Scaled[0], "K_{S}^{0}", "P");
    legendMeanPtMassScaled->AddEntry(gMPtChKstarRun2Scaled[0], "K*^{#pm}", "P");
    legendMeanPtMassScaled->AddEntry(gMPtKstarScaled[0], "K* (892)^{0}", "P");
    legendMeanPtMassScaled->AddEntry(gMPtPhiRun2Scaled[0], "#phi(1020)", "P");
    legendMeanPtMassScaled->AddEntry(gMPtProtonScaled[0], "Proton", "P");
    legendMeanPtMassScaled->AddEntry(gMPtSigmaRun2Scaled[0], "#Sigma", "P");
    legendMeanPtMassScaled->AddEntry(gMPtLambdaRun2Scaled[0], "#Lambda", "P");
    legendMeanPtMassScaled->AddEntry(gMPtXiRun2Scaled[0], "#Xi", "P");
    legendMeanPtMassScaled->AddEntry(gMPtLambda1520Scaled[0], "#Lambda(1520)", "P");
    legendMeanPtMassScaled->AddEntry(gMPtXiStarRun2Scaled[0], "#Xi(1530)", "P");
    legendMeanPtMassScaled->AddEntry(gMPtOmegaRun2Scaled[0], "#Omega", "P");
    legendMeanPtMassScaled->Draw();
    // if (isSavePlots)
    {
        cMeanPtMassScaled->SaveAs("Plots/MeanPt_Mass_nq_scaled_Run3.png");
    }

    //=================================================================================
    // =====================(<pT> - m)/n_q vs <dNch/deta> (hadrons) =========================
    //=================================================================================
    TCanvas *cMeanPtMassnqScaled = new TCanvas("cMeanPtMassnqScaled", "cMeanPtMassnqScaled", 720, 720);
    SetCanvasStyle(cMeanPtMassnqScaled, 0.15, 0.03, 0.03, 0.15);

    for (int i = 0; i < 2; i++)
    {
        ScaleGraph(gMPtKstarScaled[i], 1.0 / 2.0);
        ScaleGraph(gMPtPionScaled[i], 1.0 / 2.0);
        ScaleGraph(gMPtKaonScaled[i], 1.0 / 2.0);
        ScaleGraph(gMPtProtonScaled[i], 1.0 / 3.0);
        ScaleGraph(gMPtKshortRun2Scaled[i], 1.0 / 2.0);
        ScaleGraph(gMPtChKstarRun2Scaled[i], 1.0 / 2.0);
        ScaleGraph(gMPtPhiRun2Scaled[i], 1.0 / 2.0);
        ScaleGraph(gMPtLambdaRun2Scaled[i], 1.0 / 3.0);
        ScaleGraph(gMPtXiRun2Scaled[i], 1.0 / 3.0);
        ScaleGraph(gMPtLambda1520Scaled[i], 1.0 / 3.0);
        ScaleGraph(gMPtXiStarRun2Scaled[i], 1.0 / 3.0);
        ScaleGraph(gMPtSigmaRun2Scaled[i], 1.0 / 3.0);
        ScaleGraph(gMPtOmegaRun2Scaled[i], 1.0 / 3.0);
    }

    gMPtKstarScaled[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMPtKstarScaled[0]->GetYaxis()->SetTitle("(<#it{p}_{T}> - m)/n_{q} (GeV/#it{c})");
    gMPtKstarScaled[0]->SetMaximum(0.37);
    gMPtKstarScaled[0]->SetMinimum(-0.23);
    gMPtKstarScaled[0]->Draw("APE");
    gMPtPionScaled[0]->Draw("PE same");
    gMPtKaonScaled[0]->Draw("PE same");
    gMPtProtonScaled[0]->Draw("PE same");
    gMPtKshortRun2Scaled[0]->Draw("PE same");
    gMPtChKstarRun2Scaled[0]->Draw("PE same");
    gMPtPhiRun2Scaled[0]->Draw("PE same");
    gMPtLambdaRun2Scaled[0]->Draw("PE same");
    gMPtXiRun2Scaled[0]->Draw("PE same");
    gMPtLambda1520Scaled[0]->Draw("PE same");
    gMPtXiStarRun2Scaled[0]->Draw("PE same");
    gMPtSigmaRun2Scaled[0]->Draw("PE same");
    gMPtOmegaRun2Scaled[0]->Draw("PE same");

    gMPtKstarScaled[1]->Draw("5 same");
    gMPtPionScaled[1]->Draw("5 same");
    gMPtKaonScaled[1]->Draw("5 same");
    gMPtProtonScaled[1]->Draw("5 same");
    gMPtKshortRun2Scaled[1]->Draw("5 same");
    gMPtChKstarRun2Scaled[1]->Draw("5 same");
    gMPtPhiRun2Scaled[1]->Draw("5 same");
    gMPtLambdaRun2Scaled[1]->Draw("5 same");
    gMPtXiRun2Scaled[1]->Draw("5 same");
    gMPtLambda1520Scaled[1]->Draw("5 same");
    gMPtXiStarRun2Scaled[1]->Draw("5 same");
    gMPtSigmaRun2Scaled[1]->Draw("5 same");
    gMPtOmegaRun2Scaled[1]->Draw("5 same");

    legendMeanPtMassScaled->Draw();
    // if (isSavePlots)
    {
        cMeanPtMassnqScaled->SaveAs("Plots/MeanPt_MassScaled_nq_Run3.png");
    }

    //==================================================================================
    // ============<pt>_mult / <pt>_lowestMult vs <dNch/deta> (mesons) ================
    //==================================================================================
    TCanvas *cMeanPtRatio = new TCanvas("cMeanPtRatio", "cMeanPtRatio", 720, 720);
    SetCanvasStyle(cMeanPtRatio, 0.15, 0.03, 0.03, 0.15);
    pad1Size = 0.5, pad2Size = 0.5;
    canvas_style(cMeanPtRatio, pad1Size, pad2Size, 0.13, 0.02, 0.02, 0.25);
    cMeanPtRatio->cd(1);

    TGraphErrors **gMPtKstarLMRatio = DivideByMult(gMPtKstar, 3.69, 0.5, 2);
    TGraphErrors **gMPtPionLMRatio = DivideByMult(gMPtPion, 3.69, 0.5, 3);
    TGraphErrors **gMPtKaonLMRatio = DivideByMult(gMPtKaon, 3.69, 0.5, 3);
    TGraphErrors **gMPtProtonLMRatio = DivideByMult(gMPtProton, 3.69, 0.5, 3);
    TGraphErrors **gMPtKshortLMRatio = DivideByMult(gMPtKshortRun2, -1.0, 0.5, 3);
    TGraphErrors **gMPtChKstarLMRatio = DivideByMult(gMPtChKstarRun2, -1.0, 0.5, 3);
    TGraphErrors **gMPtPhiLMRatio = DivideByMult(gMPtPhiRun2, -1.0, 0.5, 3);
    TGraphErrors **gMPtLambdaLMRatio = DivideByMult(gMPtLambdaRun2, -1.0, 0.5, 3);
    TGraphErrors **gMPtXiLMRatio = DivideByMult(gMPtXiRun2, -1.0, 0.5, 3);
    TGraphErrors **gMPtLambda1520LMRatio = DivideByMult(gMPtLambda1520, -1.0, 0.5, 3);
    TGraphErrors **gMPtXiStarLMRatio = DivideByMult(gMPtXiStarRun2, -1.0, 0.5, 3);
    TGraphErrors **gMPtSigmaLMRatio = DivideByMult(gMPtSigmaRun2, -1.0, 0.5, 3);
    TGraphErrors **gMPtOmegaLMRatio = DivideByMult(gMPtOmegaRun2, -1.0, 0.5, 3);

    gMPtKstarLMRatio[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    gMPtKstarLMRatio[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    gMPtKstarLMRatio[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    gMPtKstarLMRatio[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    gMPtKstarLMRatio[0]->GetYaxis()->SetTitleOffset(1.5 * pad1Size);
    gMPtKstarLMRatio[0]->SetMaximum(1.85);
    gMPtKstarLMRatio[0]->SetMinimum(0.88);
    gMPtKstarLMRatio[0]->GetYaxis()->SetTitle("<#it{p}_{T}>/<#it{p}_{T}>_{LM}");
    gMPtKstarLMRatio[0]->Draw("APE");
    gMPtPionLMRatio[0]->Draw("PE same");
    gMPtKaonLMRatio[0]->Draw("PE same");
    gMPtKshortLMRatio[0]->Draw("PE same");
    gMPtChKstarLMRatio[0]->Draw("PE same");
    gMPtPhiLMRatio[0]->Draw("PE same");

    gMPtKstarLMRatio[1]->Draw("5 same");
    gMPtPionLMRatio[1]->Draw("5 same");
    gMPtKaonLMRatio[1]->Draw("5 same");
    gMPtKshortLMRatio[1]->Draw("5 same");
    gMPtChKstarLMRatio[1]->Draw("5 same");
    gMPtPhiLMRatio[1]->Draw("5 same");

    TLegend *legendMeanPtRatioMesons = new TLegend(0.18, 0.8, 0.69, 0.94);
    SetLegendStyle(legendMeanPtRatioMesons);
    legendMeanPtRatioMesons->SetNColumns(4);
    legendMeanPtRatioMesons->SetTextSize(0.027 / pad1Size);
    legendMeanPtRatioMesons->AddEntry(gMPtPionLMRatio[0], "Pion", "P");
    legendMeanPtRatioMesons->AddEntry(gMPtKaonLMRatio[0], "Kaon", "P");
    legendMeanPtRatioMesons->AddEntry(gMPtKstarLMRatio[0], "K* (892)^{0}", "P");
    legendMeanPtRatioMesons->AddEntry(gMPtKshortLMRatio[0], "K_{S}^{0}", "P");
    legendMeanPtRatioMesons->AddEntry(gMPtChKstarLMRatio[0], "K*^{#pm}", "P");
    legendMeanPtRatioMesons->AddEntry(gMPtPhiLMRatio[0], "#phi (1020)", "P");
    legendMeanPtRatioMesons->Draw();

    cMeanPtRatio->cd(2);
    SetGraphErrorStyle(gMPtProtonLMRatio[0]);
    gMPtProtonLMRatio[0]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    gMPtProtonLMRatio[0]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    gMPtProtonLMRatio[0]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    gMPtProtonLMRatio[0]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    gMPtProtonLMRatio[0]->GetYaxis()->SetTitleOffset(1.5 * pad2Size);
    gMPtProtonLMRatio[0]->SetMaximum(2.05);
    gMPtProtonLMRatio[0]->SetMinimum(0.82);
    gMPtProtonLMRatio[0]->SetMarkerColor(kMagenta);
    gMPtProtonLMRatio[0]->SetLineColor(kMagenta);
    gMPtProtonLMRatio[0]->GetYaxis()->SetTitle("<#it{p}_{T}>/<#it{p}_{T}>_{LM}");
    gMPtProtonLMRatio[0]->Draw("APE");
    gMPtXiStarLMRatio[0]->Draw("PE same");
    gMPtLambdaLMRatio[0]->Draw("PE same");
    gMPtXiLMRatio[0]->Draw("PE same");
    gMPtLambda1520LMRatio[0]->Draw("PE same");
    gMPtOmegaLMRatio[0]->Draw("PE same");
    gMPtSigmaLMRatio[0]->Draw("PE same");

    gMPtProtonLMRatio[1]->Draw("5 same");
    gMPtXiStarLMRatio[1]->Draw("5 same");
    gMPtLambdaLMRatio[1]->Draw("5 same");
    gMPtXiLMRatio[1]->Draw("5 same");
    gMPtLambda1520LMRatio[1]->Draw("5 same");
    gMPtOmegaLMRatio[1]->Draw("5 same");
    gMPtSigmaLMRatio[1]->Draw("5 same");

    TLegend *legendMeanPtRatioBaryons = new TLegend(0.18, 0.76, 0.69, 0.95);
    SetLegendStyle(legendMeanPtRatioBaryons);
    legendMeanPtRatioBaryons->SetNColumns(3);
    legendMeanPtRatioBaryons->SetTextSize(0.027 / pad2Size);
    legendMeanPtRatioBaryons->AddEntry(gMPtProtonLMRatio[0], "Proton", "P");
    legendMeanPtRatioBaryons->AddEntry(gMPtLambdaLMRatio[0], "#Lambda(1520)", "P");
    legendMeanPtRatioBaryons->AddEntry(gMPtXiLMRatio[0], "#Xi", "P");
    legendMeanPtRatioBaryons->AddEntry(gMPtLambda1520LMRatio[0], "#Lambda(1520)", "P");
    legendMeanPtRatioBaryons->AddEntry(gMPtXiStarLMRatio[0], "#Xi(1530)", "P");
    legendMeanPtRatioBaryons->AddEntry(gMPtOmegaLMRatio[0], "#Omega", "P");
    legendMeanPtRatioBaryons->AddEntry(gMPtSigmaLMRatio[0], "#Sigma(1385)", "P");
    legendMeanPtRatioBaryons->Draw();

    if (isSavePlots)
    {
        cMeanPtRatio->SaveAs("Plots/MeanPt_LowestMultRatio_Run3.png");
    }

    // //===================================================
    // // =====<pT>_mult / <pT>_lowestMult vs Mass============
    // //===================================================
    // TCanvas *cMeanPtRatioVsMass = new TCanvas("cMeanPtRatioVsMass", "cMeanPtRatioVsMass", 720, 720);
    // SetCanvasStyle(cMeanPtRatioVsMass, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gMeanPtLMRatioVsMass[2];
    // for (int i = 0; i < 2; i++)
    // {
    //     gMeanPtLMRatioVsMass[i] = new TGraphErrors();
    // }
    // for (int i = 0; i < 2; i++)
    // {
    //     gMeanPtLMRatioVsMass[i]->SetPoint(0, PDGMasses[0], gMPtKstarLMRatio[i]->GetY()[FindGraphXPoint(gMPtKstarLMRatio[i], 21.78)]); // K*0
    //     gMeanPtLMRatioVsMass[i]->SetPointError(0, PDGMassesErrors[0], gMPtKstarLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtKstarLMRatio[i], 21.78)));
    //     gMeanPtLMRatioVsMass[i]->SetPoint(1, PDGMasses[2], gMPtPionLMRatio[i]->GetY()[FindGraphXPoint(gMPtPionLMRatio[i], 21.78)]); // Pion
    //     gMeanPtLMRatioVsMass[i]->SetPointError(1, PDGMassesErrors[2], gMPtPionLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtPionLMRatio[i], 21.78)));
    //     gMeanPtLMRatioVsMass[i]->SetPoint(2, PDGMasses[3], gMPtKaonLMRatio[i]->GetY()[FindGraphXPoint(gMPtKaonLMRatio[i], 21.78)]); // Kaon
    //     gMeanPtLMRatioVsMass[i]->SetPointError(2, PDGMassesErrors[3], gMPtKaonLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtKaonLMRatio[i], 21.78)));
    //     gMeanPtLMRatioVsMass[i]->SetPoint(3, PDGMasses[4], gMPtProtonLMRatio[i]->GetY()[FindGraphXPoint(gMPtProtonLMRatio[i], 21.78)]); // Proton
    //     gMeanPtLMRatioVsMass[i]->SetPointError(3, PDGMassesErrors[4], gMPtProtonLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtProtonLMRatio[i], 21.78)));
    //     gMeanPtLMRatioVsMass[i]->SetPoint(4, PDGMasses[5], gMPtKshortLMRatio[i]->GetY()[FindGraphXPoint(gMPtKshortLMRatio[i], 21.78)]); // Kshort
    //     gMeanPtLMRatioVsMass[i]->SetPointError(4, PDGMassesErrors[5], gMPtKshortLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtKshortLMRatio[i], 21.78)));
    //     gMeanPtLMRatioVsMass[i]->SetPoint(5, PDGMasses[6], gMPtChKstarLMRatio[i]->GetY()[FindGraphXPoint(gMPtChKstarLMRatio[i], 21.78)]); // Charged K*
    //     gMeanPtLMRatioVsMass[i]->SetPointError(5, PDGMassesErrors[6], gMPtChKstarLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtChKstarLMRatio[i], 21.78)));
    //     gMeanPtLMRatioVsMass[i]->SetPoint(6, PDGMasses[7], gMPtXiStarLMRatio[i]->GetY()[FindGraphXPoint(gMPtXiStarLMRatio[i], 21.78)]); // Xi(1530)
    //     gMeanPtLMRatioVsMass[i]->SetPointError(6, PDGMassesErrors[7], gMPtXiStarLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtXiStarLMRatio[i], 21.78)));
    //     gMeanPtLMRatioVsMass[i]->SetPoint(7, PDGMasses[8], gMPtLambda1520LMRatio[i]->GetY()[FindGraphXPoint(gMPtLambda1520LMRatio[i], 21.78)]); // Lambda(1520)
    //     gMeanPtLMRatioVsMass[i]->SetPointError(7, PDGMassesErrors[8], gMPtLambda1520LMRatio[i]->GetErrorY(FindGraphXPoint(gMPtLambda1520LMRatio[i], 21.78)));
    //     gMeanPtLMRatioVsMass[i]->SetPoint(8, PDGMasses[9], gMPtRhoLMRatio[i]->GetY()[FindGraphXPoint(gMPtRhoLMRatio[i], 21.78)]); // Rho(770)
    //     gMeanPtLMRatioVsMass[i]->SetPointError(8, PDGMassesErrors[9], gMPtRhoLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtRhoLMRatio[i], 21.78)));
    //     // gMeanPtLMRatioVsMass[i]->SetPoint(8, PDGMasses[1], gMPtPhiLMRatio[i]->GetY()[FindGraphXPoint(gMPtPhiLMRatio[i], 21.78)]); // Phi
    //     // gMeanPtLMRatioVsMass[i]->SetPointError(8, PDGMassesErrors[1], gMPtPhiLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtPhiLMRatio[i], 21.78)));
    // }
    // SetGraphErrorStyle(gMeanPtLMRatioVsMass[0]);
    // gMeanPtLMRatioVsMass[0]->GetXaxis()->SetTitle("Particle Mass (GeV/#it{c}^{2})");
    // gMeanPtLMRatioVsMass[0]->GetYaxis()->SetTitle("<#it{p}_{T}>_{HM} / <#it{p}_{T}>_{LM}");
    // gMeanPtLMRatioVsMass[0]->SetMaximum(1.83);
    // gMeanPtLMRatioVsMass[0]->SetMinimum(1.24);
    // gMeanPtLMRatioVsMass[0]->SetMarkerColor(kBlue);
    // gMeanPtLMRatioVsMass[0]->SetLineColor(kBlue);
    // gMeanPtLMRatioVsMass[0]->SetMarkerStyle(21);
    // gMeanPtLMRatioVsMass[0]->Draw("AP");
    // gMeanPtLMRatioVsMass[1]->SetMarkerColor(kBlue);
    // gMeanPtLMRatioVsMass[1]->SetLineColor(kBlue);
    // gMeanPtLMRatioVsMass[1]->Draw("5 same");
    // if (isSavePlots)
    // {
    //     cMeanPtRatioVsMass->SaveAs("Plots/MeanPt_LowestMultRatio_vs_Mass_Run3.png");
    // }

    //====================================================
    //   ================Kstar/K Ratio ===================
    //====================================================

    TCanvas *cRatioKstarKaon = new TCanvas("cRatioKstarKaon", "cRatioKstarKaon", 720, 720);
    SetCanvasStyle(cRatioKstarKaon, 0.15, 0.03, 0.03, 0.15);

    TGraphErrors **gRatioKstarKaon = MakeRatioUncorr(gMYieldKstar, gMYieldKaon, 1.0, 3);
    gRatioKstarKaon[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarKaon[0]->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarKaon[0]);
    // gRatioKstarKaon[0]->GetYaxis()->SetRangeUser(0.25, 0.43);
    gRatioKstarKaon[0]->GetYaxis()->SetRangeUser(0.19, 0.46);
    gRatioKstarKaon[0]->GetXaxis()->SetLimits(0, 27);
    gRatioKstarKaon[0]->SetMarkerColor(kRed);
    gRatioKstarKaon[0]->SetLineColor(kRed);
    gRatioKstarKaon[0]->Draw("APE");
    gRatioKstarKaon[1]->SetFillStyle(0);
    gRatioKstarKaon[1]->SetLineColor(kRed);
    gRatioKstarKaon[1]->Draw("5 same");

    gKstarKaRatio_13TeV[0]->SetMarkerStyle(21);
    gKstarKaRatio_13TeV[0]->SetMarkerColor(kBlue);
    gKstarKaRatio_13TeV[0]->SetLineColor(kBlue);
    gKstarKaRatio_13TeV[0]->Draw("P same");
    gKstarKaRatio_13TeV[1]->SetLineColor(kBlue);
    gKstarKaRatio_13TeV[1]->SetFillStyle(0);
    gKstarKaRatio_13TeV[1]->Draw("5 same");

    int totalPointsKstarKaonRatio = gRatioKstarKaon[0]->GetN();
    double xLM, xHM, yHM, yLM, yHMError, yLMError, yHMStatError, yLMStatError;
    gRatioKstarKaon[1]->GetPoint(totalPointsKstarKaonRatio - 1, xLM, yLM);
    gRatioKstarKaon[1]->GetPoint(0, xHM, yHM);
    yLMError = gRatioKstarKaon[1]->GetErrorY(totalPointsKstarKaonRatio - 1);
    yHMError = gRatioKstarKaon[1]->GetErrorY(0);
    yLMStatError = gRatioKstarKaon[0]->GetErrorY(totalPointsKstarKaonRatio - 1);
    yHMStatError = gRatioKstarKaon[0]->GetErrorY(0);
    double sigmaDiffHM_LM = abs(yHM - yLM) / sqrt(pow(yHMError, 2) + pow(yLMError, 2) + pow(yHMStatError, 2) + pow(yLMStatError, 2));
    cout << "K*0/Kaon ratio difference between HM and LM: " << sigmaDiffHM_LM << " sigma" << endl;

    // K* is (K* + anit_K*)/2 but K = (K^+ + K^-), so we need to divide the denoimator by 2 as well.
    TGraphErrors *gRatioKstarKa_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gMYieldKaonEPOS_IST0, true, 1.0);
    TGraphErrors *gRatioKstarKa_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gMYieldKaonEPOS_IST0, true, 1.0);

    gRatioKstarKa_IST9->SetLineStyle(2);
    gRatioKstarKa_IST9->Draw("l same");
    gRatioKstarKa_IST9_ITY80->SetLineColor(kBlue - 6);
    gRatioKstarKa_IST9_ITY80->SetLineStyle(2);
    gRatioKstarKa_IST9_ITY80->Draw("l same");

    TLegend *legendRatio = new TLegend(0.2, 0.75, 0.5, 0.85);
    SetLegendStyle(legendRatio);
    legendRatio->SetTextSize(0.027);
    legendRatio->AddEntry(gRatioKstarKaon[0], "pp, #sqrt{s} = 13.6 TeV", "p");
    legendRatio->AddEntry(gKstarKaRatio_13TeV[0], "pp, #sqrt{s} = 13 TeV", "p");

    TLegend *legendRatio2 = new TLegend(0.55, 0.72, 0.8, 0.92);
    SetLegendStyle(legendRatio2);
    legendRatio2->SetTextSize(0.027);
    legendRatio2->AddEntry(gRatioKstarKa_IST9, "EPOS UrQMD OFF", "l");
    legendRatio2->AddEntry(gRatioKstarKa_IST9_ITY80, "EPOS UrQMD ON", "l");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarKaModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKaon], true);
        setStyle(gRatioKstarKaModel, modelStyle[model].color, modelStyle[model].style);
        gRatioKstarKaModel->Draw("l same");
        legendRatio2->AddEntry(gRatioKstarKaModel, modelLabel[model], "l");
    }

    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        TGraphErrors *gRatioKstarKaPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kKaon_epos], true, 0.5);
        setStyle(gRatioKstarKaPythiaModel, colorsPythia[imodel], lineStylesPythia[imodel]);
        gRatioKstarKaPythiaModel->Draw("l same");
        legendRatio2->AddEntry(gRatioKstarKaPythiaModel, modelLabelLocal[imodel], "l");
    }

    legendRatio->Draw();
    legendRatio2->Draw();
    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0} + #bar{K}^{*0}}{K^{+} + K^{-}}");
    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{K}");
    if (isSavePlots)
    {
        cRatioKstarKaon->SaveAs("Plots/Ratio_KstarKaon_Run3.png");
    }

    //==================================================
    // ===================Kstar / Pion Ratio ==================
    //==================================================
    TCanvas *cRatioKstarPion = new TCanvas("cRatioKstarPion", "cRatioKstarPion", 720, 720);
    SetCanvasStyle(cRatioKstarPion, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors **gRatioKstarPion = MakeRatioUncorr(gMYieldKstar, gMYieldPion, 1.0, 3);
    gRatioKstarPion[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarPion[0]->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarPion[0]);
    // gRatioKstarPion[0]->GetYaxis()->SetRangeUser(0.026, 0.059);
    gRatioKstarPion[0]->GetYaxis()->SetRangeUser(0.032, 0.072);
    gRatioKstarPion[0]->GetXaxis()->SetLimits(0, 27);
    gRatioKstarPion[0]->SetMarkerColor(kRed);
    gRatioKstarPion[0]->SetLineColor(kRed);
    gRatioKstarPion[0]->Draw("APE");
    gRatioKstarPion[1]->SetFillStyle(0);
    gRatioKstarPion[1]->SetLineColor(kRed);
    gRatioKstarPion[1]->Draw("5 same");
    gKstarPiRatio_13TeV[0]->SetMarkerStyle(21);
    gKstarPiRatio_13TeV[0]->SetMarkerColor(kBlue);
    gKstarPiRatio_13TeV[0]->SetLineColor(kBlue);
    gKstarPiRatio_13TeV[0]->Draw("P same");
    gKstarPiRatio_13TeV[1]->SetLineColor(kBlue);
    gKstarPiRatio_13TeV[1]->SetFillStyle(0);
    gKstarPiRatio_13TeV[1]->Draw("5 same");
    gKstarPiRatio_7TeV[0]->SetMarkerStyle(22);
    gKstarPiRatio_7TeV[0]->SetMarkerColor(kGreen + 2);
    gKstarPiRatio_7TeV[0]->SetLineColor(kGreen + 2);
    // gKstarPiRatio_7TeV[0]->Draw("P same");
    gKstarPiRatio_7TeV[1]->SetLineColor(kGreen + 2);
    gKstarPiRatio_7TeV[1]->SetFillStyle(0);
    // gKstarPiRatio_7TeV[1]->Draw("5 same");

    // K* is (K* + anit_K*)/2, but for pion it is (Pi^+ + Pi^-), so no need a factor 2 correction.
    TGraphErrors *gRatioKstarPi_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gMYieldPionEPOS_IST0, true, 1.0);
    TGraphErrors *gRatioKstarPi_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gMYieldPionEPOS_IST0, true, 1.0);

    gRatioKstarPi_IST9->SetLineStyle(2);
    gRatioKstarPi_IST9->Draw("l same");
    gRatioKstarPi_IST9_ITY80->SetLineColor(kBlue);
    gRatioKstarPi_IST9_ITY80->SetLineStyle(2);
    gRatioKstarPi_IST9_ITY80->Draw("l same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarPiModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kPion], true);
        setStyle(gRatioKstarPiModel, modelStyle[model].color, modelStyle[model].style);
        gRatioKstarPiModel->Draw("l same");
    }

    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        TGraphErrors *gRatioKstarPiPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kPion_epos], true, 0.5);
        setStyle(gRatioKstarPiPythiaModel, colorsPythia[imodel], lineStylesPythia[imodel]);
        gRatioKstarPiPythiaModel->Draw("l same");
    }

    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{#pi}");
    TLegend *legendRatio3 = (TLegend *)legendRatio->Clone("legendRatio3");
    legendRatio3->AddEntry(gKstarPiRatio_7TeV[0], "pp, #sqrt{s} = 7 TeV", "P");
    legendRatio3->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
        cRatioKstarPion->SaveAs("Plots/Ratio_KstarPion_Run3.png");
    }

    //======================================================================
    // ==Double YieldRatio (Mult/LM): Kstar/K, ChKstar/Kaon, XiStar/Kaon^2 ==
    // ==Double YieldRatio (Mult/LM): Kstar/Pi, ChKstar/Pi, XiStar/Pi^2 ==
    //======================================================================
    TCanvas *cYieldLMRatio = new TCanvas("cYieldLMRatio", "cYieldLMRatio", 720, 720);
    SetCanvasStyle(cYieldLMRatio, 0.15, 0.03, 0.03, 0.15);

    TGraphErrors **gYieldKstarKaLMRatio = DivideByMult(gRatioKstarKaon, -1, 0.5, 2);
    TGraphErrors **gYieldKstarPiLMRatio = DivideByMult(gRatioKstarPion, -1, 0.5, 2);

    // int totalPoints = gYieldKstarKaLMRatio[1]->GetN();
    // for (int j = 0; j < totalPoints - 1; j++)
    // {
    //     double x, yKa, yPi;
    //     gYieldKstarKaLMRatio[1]->GetPoint(j, x, yKa);
    //     gYieldKstarPiLMRatio[1]->GetPoint(j, x, yPi);
    //     double percentageDifference = abs(yKa - yPi) / yKa * 100.0;
    //     double DifferenceSigma = abs(yKa - yPi) / sqrt(pow(gYieldKstarKaLMRatio[1]->GetErrorY(j), 2) + pow(gYieldKstarPiLMRatio[1]->GetErrorY(j), 2) + pow(gYieldKstarKaLMRatio[0]->GetErrorY(j), 2) + pow(gYieldKstarPiLMRatio[0]->GetErrorY(j), 2));
    //     cout << "Multiplicity: " << x << ", K* / K: " << yKa << ", K* / Pi: " << yPi << ", Percentage Difference: " << std::round(percentageDifference * 10.0) / 10.0 << "%, Difference in Sigma: " << std::round(DifferenceSigma * 10.0) / 10.0 << endl;
    // }

    gYieldKstarKaLMRatio[0]->SetMaximum(1.33);
    gYieldKstarKaLMRatio[0]->SetMinimum(0.75);
    gYieldKstarKaLMRatio[0]->SetMarkerColor(kRed);
    gYieldKstarKaLMRatio[0]->SetLineColor(kRed);
    gYieldKstarKaLMRatio[0]->SetMarkerStyle(20);
    gYieldKstarKaLMRatio[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gYieldKstarKaLMRatio[0]->GetYaxis()->SetTitle("Y/Y_{LM}");
    gYieldKstarKaLMRatio[0]->Draw("APE");
    gYieldKstarKaLMRatio[1]->Draw("5 same");

    gYieldKstarPiLMRatio[0]->SetMarkerColor(kGreen + 2);
    gYieldKstarPiLMRatio[0]->SetLineColor(kGreen + 2);
    gYieldKstarPiLMRatio[0]->SetMarkerStyle(21);
    gYieldKstarPiLMRatio[0]->Draw("PE same");
    gYieldKstarPiLMRatio[1]->SetLineColor(kGreen + 2);
    gYieldKstarPiLMRatio[1]->SetFillStyle(0);
    gYieldKstarPiLMRatio[1]->Draw("5 same");

    TGraphErrors *gEPOS_KstarKaLMRatio_IST9 = DivideByMultModel(gRatioKstarKa_IST9, 3.69);
    TGraphErrors *gEPOS_KstarKaLMRatio_IST9_ITY80 = DivideByMultModel(gRatioKstarKa_IST9_ITY80, 3.69);
    TGraphErrors *gEPOS_KstarPiLMRatio_IST9 = DivideByMultModel(gRatioKstarPi_IST9, 3.69);
    TGraphErrors *gEPOS_KstarPiLMRatio_IST9_ITY80 = DivideByMultModel(gRatioKstarPi_IST9_ITY80, 3.69);

    gEPOS_KstarKaLMRatio_IST9->SetLineStyle(2);
    gEPOS_KstarKaLMRatio_IST9_ITY80->SetLineStyle(1);
    gEPOS_KstarPiLMRatio_IST9->SetLineStyle(2);
    gEPOS_KstarPiLMRatio_IST9_ITY80->SetLineStyle(1);

    gEPOS_KstarKaLMRatio_IST9->SetLineColor(kRed + 2);
    gEPOS_KstarKaLMRatio_IST9_ITY80->SetLineColor(kBlue - 2);
    gEPOS_KstarPiLMRatio_IST9->SetLineColor(kGreen + 2);
    gEPOS_KstarPiLMRatio_IST9_ITY80->SetLineColor(kMagenta - 2);

    gEPOS_KstarKaLMRatio_IST9->Draw("l same");
    gEPOS_KstarKaLMRatio_IST9_ITY80->Draw("l same");
    gEPOS_KstarPiLMRatio_IST9->Draw("l same");
    gEPOS_KstarPiLMRatio_IST9_ITY80->Draw("l same");

    TLegend *legendYieldLMRatio = new TLegend(0.2, 0.86, 0.95, 0.95);
    SetLegendStyle(legendYieldLMRatio);
    legendYieldLMRatio->SetTextSize(0.03);
    legendYieldLMRatio->SetNColumns(2);
    legendYieldLMRatio->AddEntry(gYieldKstarKaLMRatio[0], "K*^{0}/K", "P");
    legendYieldLMRatio->AddEntry(gYieldKstarPiLMRatio[0], "K*^{0}/#pi", "P");

    TLegend *legendYieldLMRatio2 = new TLegend(0.2, 0.8, 0.95, 0.89);
    SetLegendStyle(legendYieldLMRatio2);
    legendYieldLMRatio2->SetTextSize(0.03);
    legendYieldLMRatio2->SetNColumns(2);
    legendYieldLMRatio2->AddEntry(gEPOS_KstarKaLMRatio_IST9, "EPOS UrQMD OFF", "l");
    legendYieldLMRatio2->AddEntry(gEPOS_KstarPiLMRatio_IST9, "EPOS UrQMD OFF", "l");
    legendYieldLMRatio2->AddEntry(gEPOS_KstarKaLMRatio_IST9_ITY80, "EPOS UrQMD ON", "l");
    legendYieldLMRatio2->AddEntry(gEPOS_KstarPiLMRatio_IST9_ITY80, "EPOS UrQMD ON", "l");

    // // for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    // for (int imodel = 3; imodel < 4; imodel++)
    // {
    //     TGraphErrors *gRatioKstarKaPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kKaon_epos], true, 0.5);
    //     setStyle(gRatioKstarKaPythiaModel, colorsPythia[imodel], lineStylesPythia[imodel]);
    //     // gRatioKstarKaPythiaModel->Draw("l same");
    //     // legendRatio2->AddEntry(gRatioKstarKaPythiaModel, modelLabelLocal[imodel], "l");
    //     TGraphErrors *gYieldKstarKaLMRatioPythia = (TGraphErrors *)gRatioKstarKaPythiaModel->Clone(Form("gYieldKstarKaLMRatioPythia_%d", imodel));
    //     gYieldKstarKaLMRatioPythia = DivideByMultModel(gRatioKstarKaPythiaModel, 3.69);
    //     gYieldKstarKaLMRatioPythia->SetLineWidth(3);
    //     gYieldKstarKaLMRatioPythia->SetLineColor(kMagenta + 1);
    //     gYieldKstarKaLMRatioPythia->SetLineStyle(2);
    //     gYieldKstarKaLMRatioPythia->Draw("l same");

    //     TGraphErrors *gRatioKstarPiPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kPion_epos], true, 0.5);
    //     TGraphErrors *gYieldKstarPiLMRatioPythia = (TGraphErrors *)gRatioKstarPiPythiaModel->Clone(Form("gYieldKstarPiLMRatioPythia_%d", imodel));
    //     gYieldKstarPiLMRatioPythia = DivideByMultModel(gRatioKstarPiPythiaModel, 3.69);
    //     gYieldKstarPiLMRatioPythia->SetLineWidth(3);
    //     gYieldKstarPiLMRatioPythia->SetLineColor(kCyan + 1);
    //     gYieldKstarPiLMRatioPythia->Draw("l same");
    //     legendYieldLMRatio2->AddEntry(gYieldKstarKaLMRatioPythia, Form("%s ", modelLabelLocal[imodel]), "l");
    //     legendYieldLMRatio2->AddEntry(gYieldKstarPiLMRatioPythia, Form("%s ", modelLabelLocal[imodel]), "l");
    // }
    legendYieldLMRatio2->Draw();
    legendYieldLMRatio->Draw();

    TLine *lineDoubleRatio = new TLine(gYieldKstarKaLMRatio[0]->GetXaxis()->GetXmin(), 1.0, gYieldKstarKaLMRatio[0]->GetXaxis()->GetXmax(), 1.0);
    lineDoubleRatio->SetLineStyle(2);
    lineDoubleRatio->SetLineWidth(2);
    // lineDoubleRatio->Draw("same");
    if (isSavePlots)
    {
        cYieldLMRatio->SaveAs("Plots/Yield_LMRatio_KstarKaPi_Run3.png");
    }

    //======================================================================
    // ==YieldRatio (Mult/LM): Kstar, ChKstar, Lambda1520, XiStar, Phi, Rho==
    //======================================================================
    TCanvas *cYieldLMRatio2 = new TCanvas("cYieldLMRatio2", "cYieldLMRatio2", 720, 720);
    SetCanvasStyle(cYieldLMRatio2, 0.15, 0.03, 0.03, 0.15);
    // canvas_style(cYieldLMRatio2);
    // cYieldLMRatio2->cd(1);

    TGraphErrors **gYieldKstarLMRatio = DivideByMult(gMYieldKstar, -1.0, 0.5, 2);
    // TGraphErrors **gYieldKstarRun2LMRatio = DivideByMult(gMYieldKstarRun2, -1.0, 0.5, 3);
    TGraphErrors **gYieldChKstarLMRatio = DivideByMult(gMYieldChKstarRun2, -1.0, 0.5, 3);
    TGraphErrors **gYieldXiStarLMRatio = DivideByMult(gMYieldXiStarRun2, -1.0, 0.5, 3);
    TGraphErrors **gYieldPhiLMRatio = DivideByMult(gMYieldPhiRun2, -1.0, 0.5, 3);
    TGraphErrors **gYieldSigmaLMRatio = DivideByMult(gMYieldSigmaRun2, -1.0, 0.5, 3);
    TGraphErrors **gYieldLambda1520LMRatio = DivideByMult(gMYieldLambda1520, -1.0, 0.5, 3); // Run2

    // These are from run3
    // TGraphErrors **gYieldLambda1520LMRatio = DivideByMult(gMYieldLambda1520, -1.0, 0.5, 2);
    TGraphErrors **gYieldRhoLMRatio = DivideByMult(gMYieldRho, -1.0, 0.5, 2);

    gYieldKstarLMRatio[0]->SetMaximum(21.5);
    gYieldKstarLMRatio[0]->SetMinimum(0.0);
    gYieldKstarLMRatio[0]->SetMarkerColor(kRed);
    gYieldKstarLMRatio[0]->SetLineColor(kRed);
    gYieldKstarLMRatio[0]->SetMarkerStyle(20);
    gYieldKstarLMRatio[0]->GetXaxis()->SetLimits(0, 27);
    gYieldKstarLMRatio[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gYieldKstarLMRatio[0]->GetYaxis()->SetTitle("Y/Y_{LM}");
    // pad1Size = 0.7, pad2Size = 0.3;
    // gYieldKstarLMRatio[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    // gYieldKstarLMRatio[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    // gYieldKstarLMRatio[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    // gYieldKstarLMRatio[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    // gYieldKstarLMRatio[0]->GetYaxis()->SetTitleOffset(1.7 * pad1Size);
    gYieldKstarLMRatio[0]->Draw("APE");
    gYieldKstarLMRatio[1]->SetLineColor(kRed);
    gYieldKstarLMRatio[1]->SetFillStyle(0);
    gYieldKstarLMRatio[1]->Draw("5 same");
    gYieldChKstarLMRatio[0]->SetMarkerColor(kBlue);
    gYieldChKstarLMRatio[0]->SetLineColor(kBlue);
    gYieldChKstarLMRatio[0]->SetMarkerStyle(21);
    gYieldChKstarLMRatio[0]->Draw("PE same");
    gYieldChKstarLMRatio[1]->SetLineColor(kBlue);
    gYieldChKstarLMRatio[1]->SetFillStyle(0);
    gYieldChKstarLMRatio[1]->Draw("5 same");
    gYieldXiStarLMRatio[0]->SetMarkerColor(kBrown);
    gYieldXiStarLMRatio[0]->SetLineColor(kBrown);
    gYieldXiStarLMRatio[0]->SetMarkerStyle(22);
    gYieldXiStarLMRatio[0]->Draw("PE same");
    gYieldXiStarLMRatio[1]->SetLineColor(kBrown);
    gYieldXiStarLMRatio[1]->SetFillStyle(0);
    gYieldXiStarLMRatio[1]->Draw("5 same");
    gYieldLambda1520LMRatio[0]->SetMarkerColor(kGreen + 2);
    gYieldLambda1520LMRatio[0]->SetLineColor(kGreen + 2);
    gYieldLambda1520LMRatio[0]->SetMarkerStyle(22);
    gYieldLambda1520LMRatio[0]->Draw("PE same");
    gYieldLambda1520LMRatio[1]->SetLineColor(kGreen + 2);
    gYieldLambda1520LMRatio[1]->SetFillStyle(0);
    gYieldLambda1520LMRatio[1]->Draw("5 same");
    gYieldPhiLMRatio[0]->SetMarkerColor(kMagenta);
    gYieldPhiLMRatio[0]->SetLineColor(kMagenta);
    gYieldPhiLMRatio[0]->SetMarkerStyle(47);
    gYieldPhiLMRatio[0]->Draw("PE same");
    gYieldPhiLMRatio[1]->SetLineColor(kMagenta);
    gYieldPhiLMRatio[1]->SetFillStyle(0);
    gYieldPhiLMRatio[1]->Draw("5 same");
    // gYieldRhoLMRatio[0]->SetMarkerColor(kCyan + 1);
    // gYieldRhoLMRatio[0]->SetLineColor(kCyan + 1);
    // gYieldRhoLMRatio[0]->SetMarkerStyle(25);
    // gYieldRhoLMRatio[0]->Draw("PE same");
    // gYieldRhoLMRatio[1]->SetLineColor(kCyan + 1);
    // gYieldRhoLMRatio[1]->SetFillStyle(0);
    // gYieldRhoLMRatio[1]->Draw("5 same");
    gYieldSigmaLMRatio[0]->SetMarkerColor(kAzure + 1);
    gYieldSigmaLMRatio[0]->SetLineColor(kAzure + 1);
    gYieldSigmaLMRatio[0]->SetMarkerStyle(26);
    gYieldSigmaLMRatio[0]->Draw("PE same");
    gYieldSigmaLMRatio[1]->SetLineColor(kAzure + 1);
    gYieldSigmaLMRatio[1]->SetFillStyle(0);
    gYieldSigmaLMRatio[1]->Draw("5 same");

    TLegend *legendYieldLMRatio3 = new TLegend(0.2, 0.78, 0.6, 0.9);
    SetLegendStyle(legendYieldLMRatio3);
    // legendYieldLMRatio3->SetTextSize(0.027 / pad1Size);
    legendYieldLMRatio3->SetTextSize(0.027);
    legendYieldLMRatio3->SetNColumns(3);
    // legendYieldLMRatio3->AddEntry(gYieldRhoLMRatio[0], "#rho", "P");
    legendYieldLMRatio3->AddEntry(gYieldKstarLMRatio[0], "K*^{0}", "P");
    legendYieldLMRatio3->AddEntry(gYieldChKstarLMRatio[0], "K*^{#pm}", "P");
    legendYieldLMRatio3->AddEntry(gYieldSigmaLMRatio[0], "#Sigma(1385)", "P");
    legendYieldLMRatio3->AddEntry(gYieldLambda1520LMRatio[0], "#Lambda(1520)", "P");
    legendYieldLMRatio3->AddEntry(gYieldXiStarLMRatio[0], "#Xi(1530)", "P");
    legendYieldLMRatio3->AddEntry(gYieldPhiLMRatio[0], "#phi", "P");
    legendYieldLMRatio3->Draw();

    // cYieldLMRatio2->cd(2);
    // // These are run2, so will use run2 kstar
    // TGraphErrors **gRatio2ChKstarKstar = MakeRatioUncorr(gYieldChKstarLMRatio, gYieldKstarRun2LMRatio, 1.0, 2);
    // // TGraphErrors **gRatio2XiStarKstar = MakeRatioUncorr(gYieldXiStarLMRatio, gYieldKstarRun2LMRatio, 1.0, 2);
    // // TGraphErrors **gRatio2PhiKstar = MakeRatioUncorr(gYieldPhiLMRatio, gYieldKstarRun2LMRatio, 1.0, 2);
    // // // These are run3, so will use run3 kstar
    // // TGraphErrors **gRatio2Lambda1520Kstar = MakeRatioUncorr(gYieldLambda1520LMRatio, gYieldKstarLMRatio, 1.0, 2);
    // // TGraphErrors **gRatio2RhoKstar = MakeRatioUncorr(gYieldRhoLMRatio, gYieldKstarLMRatio, 1.0, 2);

    // gRatio2ChKstarKstar[0]->SetMaximum(2.58);
    // gRatio2ChKstarKstar[0]->SetMinimum(0.65);
    // gRatio2ChKstarKstar[0]->GetYaxis()->SetNdivisions(505);
    // gRatio2ChKstarKstar[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatio2ChKstarKstar[0]->GetYaxis()->SetTitle("Ratio to K*^{0}");
    // gRatio2ChKstarKstar[0]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    // gRatio2ChKstarKstar[0]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    // gRatio2ChKstarKstar[0]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    // gRatio2ChKstarKstar[0]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    // gRatio2ChKstarKstar[0]->GetYaxis()->SetTitleOffset(1.7 * pad2Size);
    // gRatio2ChKstarKstar[0]->GetXaxis()->SetTitleOffset(1.15);
    // gRatio2ChKstarKstar[0]->SetMarkerColor(kBlue);
    // gRatio2ChKstarKstar[0]->SetLineColor(kBlue);
    // gRatio2ChKstarKstar[0]->SetMarkerStyle(21);
    // gRatio2ChKstarKstar[0]->Draw("APE");
    // gRatio2ChKstarKstar[1]->SetLineColor(kBlue);
    // gRatio2ChKstarKstar[1]->SetFillStyle(0);
    // gRatio2ChKstarKstar[1]->Draw("5 same");
    // // gRatio2XiStarKstar[0]->SetMarkerColor(kOrange - 2);
    // // gRatio2XiStarKstar[0]->SetLineColor(kOrange - 2);
    // // gRatio2XiStarKstar[0]->SetMarkerStyle(22);
    // // gRatio2XiStarKstar[0]->Draw("PE same");
    // // gRatio2XiStarKstar[1]->SetLineColor(kOrange - 2);
    // // gRatio2XiStarKstar[1]->SetFillStyle(0);
    // // gRatio2XiStarKstar[1]->Draw("5 same");
    // // gRatio2Lambda1520Kstar[0]->SetMarkerColor(kGreen + 2);
    // // gRatio2Lambda1520Kstar[0]->SetLineColor(kGreen + 2);
    // // gRatio2Lambda1520Kstar[0]->SetMarkerStyle(22);
    // // gRatio2Lambda1520Kstar[0]->Draw("PE same");
    // // gRatio2Lambda1520Kstar[1]->SetLineColor(kGreen + 2);
    // // gRatio2Lambda1520Kstar[1]->SetFillStyle(0);
    // // gRatio2Lambda1520Kstar[1]->Draw("5 same");
    // // gRatio2PhiKstar[0]->SetMarkerColor(kMagenta);
    // // gRatio2PhiKstar[0]->SetLineColor(kMagenta);
    // // gRatio2PhiKstar[0]->SetMarkerStyle(47);
    // // gRatio2PhiKstar[0]->Draw("PE same");
    // // gRatio2PhiKstar[1]->SetLineColor(kMagenta);
    // // gRatio2PhiKstar[1]->SetFillStyle(0);
    // // gRatio2PhiKstar[1]->Draw("5 same");
    // // gRatio2RhoKstar[0]->SetMarkerColor(kCyan + 1);
    // // gRatio2RhoKstar[0]->SetLineColor(kCyan + 1);
    // // gRatio2RhoKstar[0]->SetMarkerStyle(25);
    // // gRatio2RhoKstar[0]->Draw("PE same");
    // // gRatio2RhoKstar[1]->SetLineColor(kCyan + 1);
    // // gRatio2RhoKstar[1]->SetFillStyle(0);
    // // gRatio2RhoKstar[1]->Draw("5 same");

    // TLine *line = new TLine(gRatio2ChKstarKstar[0]->GetXaxis()->GetXmin(), 1.0, gRatio2ChKstarKstar[0]->GetXaxis()->GetXmax(), 1.0);
    // line->SetLineStyle(2);
    // line->SetLineColor(kRed);
    // line->SetLineWidth(2);
    // line->Draw("same");
    if (isSavePlots)
    {
        cYieldLMRatio2->SaveAs("Plots/Yield_LMRatio2.png");
    }

    //======================================================================
    // ==Yield Ratio (HM/LM) vs lifetime: Kstar, ChKstar, Sigma, XiStar, Phi
    //======================================================================
    TCanvas *cYieldLifetime = new TCanvas("cYieldLifetime", "cYieldLifetime", 720, 720);
    SetCanvasStyle(cYieldLifetime, 0.15, 0.03, 0.03, 0.15);

    // K*+-, K*0, Sigma(1385), Lambda(1520), XiStar, Phi
    float lifetime[6] = {3.9, 4.2, 5.25, 12.6, 22.0, 46.2};

    TGraphErrors *gYieldVsLifetime[2];
    for (int i = 0; i < 2; i++)
    {
        gYieldVsLifetime[i] = new TGraphErrors(6);
    }

    //------------------------------------------------------------------
    // Fill graphs
    //------------------------------------------------------------------
    for (int j = 0; j < 2; j++)
    {
        gYieldVsLifetime[j]->SetPoint(0, lifetime[0], gYieldChKstarLMRatio[j]->GetY()[gYieldChKstarLMRatio[j]->GetN() - 1]);
        gYieldVsLifetime[j]->SetPointError(0, 0.3, gYieldChKstarLMRatio[j]->GetErrorY(gYieldChKstarLMRatio[j]->GetN() - 1));

        gYieldVsLifetime[j]->SetPoint(1, lifetime[1], gYieldKstarLMRatio[j]->GetY()[0]);
        gYieldVsLifetime[j]->SetPointError(1, 0.3, gYieldKstarLMRatio[j]->GetErrorY(0));

        gYieldVsLifetime[j]->SetPoint(2, lifetime[2], gYieldSigmaLMRatio[j]->GetY()[0]);
        gYieldVsLifetime[j]->SetPointError(2, 0.3, gYieldSigmaLMRatio[j]->GetErrorY(0));

        gYieldVsLifetime[j]->SetPoint(3, lifetime[3], gYieldLambda1520LMRatio[j]->GetY()[0]);
        gYieldVsLifetime[j]->SetPointError(3, 0.3, gYieldLambda1520LMRatio[j]->GetErrorY(0));

        gYieldVsLifetime[j]->SetPoint(4, lifetime[4], gYieldXiStarLMRatio[j]->GetY()[0]);
        gYieldVsLifetime[j]->SetPointError(4, 0.3, gYieldXiStarLMRatio[j]->GetErrorY(0));

        gYieldVsLifetime[j]->SetPoint(5, lifetime[5], gYieldPhiLMRatio[j]->GetY()[0]);
        gYieldVsLifetime[j]->SetPointError(5, 0.3, gYieldPhiLMRatio[j]->GetErrorY(0));
    }

    //------------------------------------------------------------------
    // Draw only axes from the first graph
    //------------------------------------------------------------------
    SetGraphErrorStyle(gYieldVsLifetime[0]);

    gYieldVsLifetime[0]->SetMinimum(3.5);
    gYieldVsLifetime[0]->SetMaximum(21.5);

    gYieldVsLifetime[0]->SetMarkerSize(0);
    gYieldVsLifetime[0]->SetLineColor(0);

    gYieldVsLifetime[0]->GetXaxis()->SetTitle("Lifetime (fm/c)");
    gYieldVsLifetime[0]->GetYaxis()->SetTitle("Y_{HM}/Y_{LM}");

    gYieldVsLifetime[0]->Draw("AP");

    //------------------------------------------------------------------
    // Particle styles
    //------------------------------------------------------------------
    Color_t colors[6] = {
        kRed + 1,
        kBlue + 1,
        kAzure + 7,
        kGreen + 2,
        kBrown,
        kMagenta};

    Style_t markers[6] = {
        20, // K*±
        21, // K*0
        22, // Sigma
        25, // Lambda1520
        23, // Xi*
        33  // Phi
    };

    const char *names[6] = {
        "K^{*#pm}",
        "K^{*0}",
        "#Sigma(1385)",
        "#Lambda(1520)",
        "#Xi^{*}",
        "#phi"};

    //------------------------------------------------------------------
    // Create one graph per particle
    //------------------------------------------------------------------
    TGraphErrors *gParticleStat[6];
    TGraphErrors *gParticleSyst[6];

    for (int i = 0; i < 6; i++)
    {
        double xStat, yStat;
        double xSyst, ySyst;

        gYieldVsLifetime[0]->GetPoint(i, xStat, yStat);
        gYieldVsLifetime[1]->GetPoint(i, xSyst, ySyst);

        gParticleStat[i] = new TGraphErrors(1);
        gParticleSyst[i] = new TGraphErrors(1);

        // Statistical
        gParticleStat[i]->SetPoint(0, xStat, yStat);
        gParticleStat[i]->SetPointError(0, gYieldVsLifetime[0]->GetErrorX(i), gYieldVsLifetime[0]->GetErrorY(i));

        // Systematic
        gParticleSyst[i]->SetPoint(0, xSyst, ySyst);
        gParticleSyst[i]->SetPointError(0, gYieldVsLifetime[1]->GetErrorX(i), gYieldVsLifetime[1]->GetErrorY(i));

        //-----------------------------
        // Systematic band
        //-----------------------------
        gParticleSyst[i]->SetFillColorAlpha(colors[i], 0.25);
        gParticleSyst[i]->SetLineColor(colors[i]);
        gParticleSyst[i]->SetMarkerSize(0);

        //-----------------------------
        // Statistical point
        //-----------------------------
        gParticleStat[i]->SetMarkerStyle(markers[i]);
        gParticleStat[i]->SetMarkerSize(2.0);

        gParticleStat[i]->SetMarkerColor(colors[i]);
        gParticleStat[i]->SetLineColor(colors[i]);
        gParticleStat[i]->SetLineWidth(2);
    }

    //------------------------------------------------------------------
    // Draw systematic uncertainties first
    //------------------------------------------------------------------
    for (int i = 0; i < 6; i++)
    {
        gParticleSyst[i]->Draw("2 SAME");
    }

    //------------------------------------------------------------------
    // Draw statistical uncertainties + markers
    //------------------------------------------------------------------
    for (int i = 0; i < 6; i++)
    {
        gParticleStat[i]->Draw("PZ SAME");
    }

    //------------------------------------------------------------------
    // Labels
    //------------------------------------------------------------------
    TLatex latex2;
    latex2.SetTextSize(0.03);
    latex2.SetTextAlign(23);

    for (int i = 0; i < 6; i++)
    {
        double x, y;
        gParticleStat[i]->GetPoint(0, x, y);

        latex2.SetTextColor(colors[i]);

        if (i == 0)
            latex2.DrawLatex(x - 0.5, y - 0.8, names[i]);

        else if (i == 1)
            latex2.DrawLatex(x + 1.8, y - 0.4, names[i]);

        else if (i == 2)
            latex2.DrawLatex(x, y - 2.2, names[i]);

        else if (i == 3)
            latex2.DrawLatex(x, y - 1.0, names[i]);

        else if (i == 4)
            latex2.DrawLatex(x, y - 2.2, names[i]);

        else if (i == 5)
            latex2.DrawLatex(x, y - 1.0, names[i]);
    }
    if (isSavePlots)
    {
        cYieldLifetime->SaveAs("Plots/Yield_Lifetime_Ratio.png");
    }

    //======================================================================
    // ==Double YieldRatio (Mult/LM): K*0/K, K*+-/K0s, Sigma/Lambda, Xi*/Xi, Phi/K
    //======================================================================
    TCanvas *cYieldLMRatio3 = new TCanvas("cYieldLMRatio3", "cYieldLMRatio3", 720, 720);
    SetCanvasStyle(cYieldLMRatio3, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors **gYieldChKstarK0sLMRatio = DivideByMult(gYieldChKstarKshortRatio, -1, 0.5, 3);
    TGraphErrors **gYieldSigmaLambdaLMRatio = DivideByMult(gYieldSigmaLambdaRatio, -1, 0.5, 3);
    TGraphErrors **gYieldXiStarXiLMRatio = DivideByMult(gYieldXiStarXiRatio, -1, 0.5, 3);
    TGraphErrors **gYieldPhiKaLMRatio = DivideByMult(gYieldPhiKaRatio, -1, 0.5, 3);
    TGraphErrors **gYieldLambda1520LambdaLMRatio = DivideByMult(gYieldLambda1520LambdaRatio, -1, 0.5, 2);

    gYieldKstarKaLMRatio[0]->SetMaximum(1.89);
    gYieldKstarKaLMRatio[0]->SetMinimum(0.48);
    gYieldKstarKaLMRatio[0]->SetLineColor(kRed);
    gYieldKstarKaLMRatio[0]->SetMarkerColor(kRed);
    gYieldKstarKaLMRatio[0]->Draw("APE");
    gYieldKstarKaLMRatio[1]->SetLineColor(kRed);
    gYieldKstarKaLMRatio[1]->Draw("5 same");
    gYieldChKstarK0sLMRatio[0]->SetMarkerColor(kBlue + 1);
    gYieldChKstarK0sLMRatio[0]->SetLineColor(kBlue + 1);
    gYieldChKstarK0sLMRatio[0]->SetMarkerStyle(21);
    gYieldChKstarK0sLMRatio[0]->Draw("PE same");
    gYieldChKstarK0sLMRatio[1]->SetLineColor(kBlue + 1);
    gYieldChKstarK0sLMRatio[1]->SetFillStyle(0);
    gYieldChKstarK0sLMRatio[1]->Draw("5 same");
    gYieldSigmaLambdaLMRatio[0]->SetMarkerColor(kAzure + 7);
    gYieldSigmaLambdaLMRatio[0]->SetLineColor(kAzure + 7);
    gYieldSigmaLambdaLMRatio[0]->SetMarkerStyle(22);
    gYieldSigmaLambdaLMRatio[0]->Draw("PE same");
    gYieldSigmaLambdaLMRatio[1]->SetLineColor(kAzure + 7);
    gYieldSigmaLambdaLMRatio[1]->SetFillStyle(0);
    gYieldSigmaLambdaLMRatio[1]->Draw("5 same");
    gYieldLambda1520LambdaLMRatio[0]->SetMarkerColor(kGreen + 2);
    gYieldLambda1520LambdaLMRatio[0]->SetLineColor(kGreen + 2);
    gYieldLambda1520LambdaLMRatio[0]->SetMarkerStyle(22);
    gYieldLambda1520LambdaLMRatio[0]->Draw("PE same");
    gYieldLambda1520LambdaLMRatio[1]->SetLineColor(kGreen + 2);
    gYieldLambda1520LambdaLMRatio[1]->SetFillStyle(0);
    gYieldLambda1520LambdaLMRatio[1]->Draw("5 same");
    gYieldXiStarXiLMRatio[0]->SetMarkerColor(kBrown);
    gYieldXiStarXiLMRatio[0]->SetLineColor(kBrown);
    gYieldXiStarXiLMRatio[0]->SetMarkerStyle(23);
    gYieldXiStarXiLMRatio[0]->Draw("PE same");
    gYieldXiStarXiLMRatio[1]->SetLineColor(kBrown);
    gYieldXiStarXiLMRatio[1]->SetFillStyle(0);
    gYieldXiStarXiLMRatio[1]->Draw("5 same");
    gYieldPhiKaLMRatio[0]->SetMarkerColor(kMagenta);
    gYieldPhiKaLMRatio[0]->SetLineColor(kMagenta);
    gYieldPhiKaLMRatio[0]->SetMarkerStyle(47);
    gYieldPhiKaLMRatio[0]->Draw("PE same");
    gYieldPhiKaLMRatio[1]->SetLineColor(kMagenta);
    gYieldPhiKaLMRatio[1]->SetFillStyle(0);
    gYieldPhiKaLMRatio[1]->Draw("5 same");

    TLegend *legendYieldLMRatio4 = new TLegend(0.2, 0.78, 0.8, 0.9);
    SetLegendStyle(legendYieldLMRatio4);
    legendYieldLMRatio4->SetTextSize(0.027);
    legendYieldLMRatio4->SetNColumns(2);
    legendYieldLMRatio4->AddEntry(gYieldKstarKaLMRatio[0], "K*^{0}/K", "P");
    legendYieldLMRatio4->AddEntry(gYieldChKstarK0sLMRatio[0], "K*^{#pm}/K_{S}^{0}", "P");
    legendYieldLMRatio4->AddEntry(gYieldSigmaLambdaLMRatio[0], "#Sigma(1385)/#Lambda", "P");
    legendYieldLMRatio4->AddEntry(gYieldLambda1520LambdaLMRatio[0], "#Lambda(1520)/#Lambda", "P");
    legendYieldLMRatio4->AddEntry(gYieldXiStarXiLMRatio[0], "#Xi(1530)/#Xi", "P");
    legendYieldLMRatio4->AddEntry(gYieldPhiKaLMRatio[0], "#phi/K", "P");
    legendYieldLMRatio4->Draw();
    if (isSavePlots)
    {
        cYieldLMRatio3->SaveAs("Plots/Yield_LMRatio3.png");
    }

    // // /*
    // //====================================================
    // //   ================Kstar/ Kshort Ratio ===================
    // //====================================================
    // TCanvas *cRatioKstarKshort = new TCanvas("cRatioKstarKshort", "cRatioKstarKshort", 720, 720);
    // SetCanvasStyle(cRatioKstarKshort, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioKstarKshort = MakeRatio(gMYieldKstar[0], gMYieldKshort[0], false);
    // TGraphErrors *gRatioKstarKshort_sys = MakeRatio(gMYieldKstar[1], gMYieldKshort[1], false);
    // gRatioKstarKshort->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioKstarKshort->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioKstarKshort);
    // gRatioKstarKshort->GetYaxis()->SetRangeUser(0.18, 0.46);
    // gRatioKstarKshort->GetXaxis()->SetLimits(0, 27);
    // gRatioKstarKshort->SetMarkerColor(kRed);
    // gRatioKstarKshort->SetLineColor(kRed);
    // gRatioKstarKshort->Draw("APE");
    // gRatioKstarKshort_sys->SetFillStyle(0);
    // gRatioKstarKshort_sys->SetLineColor(kRed);
    // gRatioKstarKshort_sys->Draw("5 same");
    // gKstarKshortRatio_13TeV[0]->SetMarkerStyle(21);
    // gKstarKshortRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gKstarKshortRatio_13TeV[0]->SetLineColor(kBlue);
    // gKstarKshortRatio_13TeV[0]->Draw("P same");
    // gKstarKshortRatio_13TeV[1]->SetLineColor(kBlue);
    // gKstarKshortRatio_13TeV[1]->SetFillStyle(0);
    // gKstarKshortRatio_13TeV[1]->Draw("5 same");

    // // K* is (K* + anit_K*)/2 and Kshort does not have antiparticle, so no need for correction.
    // TGraphErrors *gRatioKstarKshort_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gMYieldKshortEPOS_IST0, true);
    // TGraphErrors *gRatioKstarKshort_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gMYieldKshortEPOS_IST0, true);

    // gRatioKstarKshort_IST9->SetLineStyle(2);
    // gRatioKstarKshort_IST9->Draw("l same");
    // gRatioKstarKshort_IST9_ITY80->SetLineColor(kBlue);
    // gRatioKstarKshort_IST9_ITY80->SetLineStyle(2);
    // gRatioKstarKshort_IST9_ITY80->Draw("l same");

    // // Here K*0 is not averaged.
    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioKstarKshortModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKshort], true, 2);
    //     setStyle(gRatioKstarKshortModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioKstarKshortModel->Draw("l same");
    //     // legendRatio->AddEntry(gRatioKstarKshortModel, modelLabel[model], "l");
    // }

    // for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    // {
    //     TGraphErrors *gRatioKstarKshortPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kKshort_epos], true);
    //     setStyle(gRatioKstarKshortPythiaModel, modelStyle[imodel].color, modelStyle[imodel].style);
    //     gRatioKstarKshortPythiaModel->Draw("l same");
    // }

    // // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0} + #bar{K}^{*0}}{2K_{S}^{0}}");
    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{K_{S}^{0}}");
    // legendRatio->Draw();
    // legendRatio2->Draw();
    // if (isSavePlots)
    // {
    //     cRatioKstarKshort->SaveAs("Plots/Ratio_KstarKshort_Run3.png");
    // }

    // //====================================================
    // // ====================Kstar/Phi Ratio ===================
    // //====================================================

    // TCanvas *cRatioKstarPhi = new TCanvas("cRatioKstarPhi", "cRatioKstarPhi", 720, 720);
    // SetCanvasStyle(cRatioKstarPhi, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioKstarPhi = MakeRatio(gMYieldKstar[0], gMYieldPhi[0], false);
    // TGraphErrors *gRatioKstarPhi_sys = MakeRatio(gMYieldKstar[1], gMYieldPhi[1], false);
    // gRatioKstarPhi->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioKstarPhi->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioKstarPhi);
    // gRatioKstarPhi->GetYaxis()->SetRangeUser(1.2, 6.7);
    // gRatioKstarPhi->GetXaxis()->SetLimits(0, 27);
    // gRatioKstarPhi->Draw("APE");
    // gRatioKstarPhi_sys->SetFillStyle(0);
    // gRatioKstarPhi_sys->Draw("5 same");

    // ////Here Kstar is average
    // TGraphErrors *gRatioKstarPhi_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gEPOS_Yield[9][kITY0][kPhi_epos], true);
    // TGraphErrors *gRatioKstarPhi_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gEPOS_Yield[9][kITY80][kPhi_epos], true);

    // gRatioKstarPhi_IST9->SetLineStyle(2);
    // gRatioKstarPhi_IST9->Draw("l same");
    // gRatioKstarPhi_IST9_ITY80->SetLineColor(kBlue);
    // gRatioKstarPhi_IST9_ITY80->SetLineStyle(2);
    // gRatioKstarPhi_IST9_ITY80->Draw("l same");

    // // Here K*0 is not averaged.
    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioKstarPhiModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kPhi], true, 2.0);
    //     setStyle(gRatioKstarPhiModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioKstarPhiModel->Draw("l same");
    // }
    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{#phi}");
    // legendRatio->Draw();
    // legendRatio2->Draw();
    // if (isSavePlots)
    // {
    //     cRatioKstarPhi->SaveAs("Plots/Ratio_KstarPhi_Run3.png");
    // }

    // //====================================================
    // //  ==================Kstar / Charged Kstar Ratio===================
    // //====================================================

    // TCanvas *cRatioKstarChargedKstar = new TCanvas("cRatioKstarChargedKstar", "cRatioKstarChargedKstar", 720, 720);
    // SetCanvasStyle(cRatioKstarChargedKstar, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioKstarChargedKstar = MakeRatio(gMYieldKstar[0], gMYieldChKstar[0], false);
    // TGraphErrors *gRatioKstarChargedKstar_sys = MakeRatio(gMYieldKstar[1], gMYieldChKstar[1], false);
    // gRatioKstarChargedKstar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioKstarChargedKstar->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioKstarChargedKstar);
    // gRatioKstarChargedKstar->GetYaxis()->SetRangeUser(0.5, 1.4);
    // gRatioKstarChargedKstar->GetXaxis()->SetLimits(0, 27);
    // gRatioKstarChargedKstar->Draw("APE");
    // gRatioKstarChargedKstar_sys->SetFillStyle(0);
    // gRatioKstarChargedKstar_sys->Draw("5 same");

    // // factor 2 is due to average in K*0 and another factor 2 is due wrong BR taken 0.333 instead of 0.666.
    // TGraphErrors *gRatioKstarChargedKstar_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gEPOS_Yield[9][kITY0][kKstarPM_epos], true, 0.25);
    // TGraphErrors *gRatioKstarChargedKstar_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gEPOS_Yield[9][kITY80][kKstarPM_epos], true, 0.25);

    // gRatioKstarChargedKstar_IST9->SetLineStyle(2);
    // gRatioKstarChargedKstar_IST9->Draw("l same");
    // gRatioKstarChargedKstar_IST9_ITY80->SetLineColor(kBlue);
    // gRatioKstarChargedKstar_IST9_ITY80->SetLineStyle(2);
    // gRatioKstarChargedKstar_IST9_ITY80->Draw("l same");

    // // Some issue in K*+- simulation. It has very less yield than K*0 somehow (a facotr of 700 difference.)
    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioKstarChargedKstarModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKstarPM], true);
    //     setStyle(gRatioKstarChargedKstarModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioKstarChargedKstarModel->Draw("l same");
    // }
    // for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    // {
    //     TGraphErrors *gRatioKstarChargedKstarPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kKstarPM_epos], true, 0.5);
    //     setStyle(gRatioKstarChargedKstarPythiaModel, modelStyle[imodel].color, modelStyle[imodel].style);
    //     gRatioKstarChargedKstarPythiaModel->Draw("l same");
    // }
    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{K^{* #pm}}");
    // legendRatio->Draw();
    // legendRatio2->Draw();
    // if (isSavePlots)
    // {
    //     cRatioKstarChargedKstar->SaveAs("Plots/Ratio_KstarChargedKstar_Run3.png");
    // }

    // //================================================
    // // ==================Kstar/ Xi* Ratio===================
    // //====================================================
    // TCanvas *cRatioKstarXiStar = new TCanvas("cRatioKstarXiStar", "cRatioKstarXiStar", 720, 720);
    // SetCanvasStyle(cRatioKstarXiStar, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioKstarXiStar = MakeRatio(gMYieldKstar[0], gMYieldXiStar[0], false);
    // TGraphErrors *gRatioKstarXiStar_sys = MakeRatio(gMYieldKstar[1], gMYieldXiStar[1], false);
    // gRatioKstarXiStar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioKstarXiStar->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioKstarXiStar);
    // gRatioKstarXiStar->GetYaxis()->SetRangeUser(15, 83);
    // gRatioKstarXiStar->GetXaxis()->SetLimits(0, 27);
    // gRatioKstarXiStar->Draw("APE");
    // gRatioKstarXiStar_sys->SetFillStyle(0);
    // gRatioKstarXiStar_sys->Draw("5 same");

    // TGraphErrors *gRatioKstarXiStar_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gEPOS_Yield[9][kITY0][kXi1530_epos], true, 0.5);
    // TGraphErrors *gRatioKstarXiStar_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gEPOS_Yield[9][kITY80][kXi1530_epos], true, 0.5);
    // TGraphErrors *gRatioKstarXiStar_IST9_ITY81 = MakeRatio(gEPOS_Yield[9][kITY81][kKstar_epos], gEPOS_Yield[9][kITY81][kXi1530_epos], true, 0.5);
    // TGraphErrors *gRatioKstarXiStar_IST6 = MakeRatio(gEPOS_Yield[6][kITY0][kKstar_epos], gEPOS_Yield[6][kITY0][kXi1530_epos], true, 0.5);

    // gRatioKstarXiStar_IST9->SetLineStyle(2);
    // gRatioKstarXiStar_IST9->Draw("l same");
    // gRatioKstarXiStar_IST9_ITY80->SetLineColor(kBlue);
    // gRatioKstarXiStar_IST9_ITY80->SetLineStyle(2);
    // gRatioKstarXiStar_IST9_ITY80->Draw("l same");
    // gRatioKstarXiStar_IST9_ITY81->SetLineColor(kGreen + 2);
    // gRatioKstarXiStar_IST9_ITY81->SetLineStyle(2);

    // // In this model, the Kstar is not divided by 2
    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioKstarXiStarModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kXi1530], true, 2);
    //     setStyle(gRatioKstarXiStarModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioKstarXiStarModel->Draw("l same");
    // }
    // for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    // {
    //     TGraphErrors *gRatioKstarXiStarPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kXi1530_epos], true);
    //     setStyle(gRatioKstarXiStarPythiaModel, modelStyle[imodel].color, modelStyle[imodel].style);
    //     gRatioKstarXiStarPythiaModel->Draw("l same");
    // }
    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{#Xi(1530)^{0}}");
    // legendRatio->Clear();
    // legendRatio->AddEntry(gRatioKstarXiStar, "pp, #sqrt{s} = 13.6 TeV", "p");
    // legendRatio->Draw();
    // legendRatio2->Draw();
    // if (isSavePlots)
    // {
    //     cRatioKstarXiStar->SaveAs("Plots/Ratio_KstarXiStar_Run3.png");
    // }

    // /*
    // //====================================================
    // // =================Kaon / Pion Ratio ==================
    // //====================================================
    // TCanvas *cRatioKaonPion = new TCanvas("cRatioKaonPion", "cRatioKaonPion", 720, 720);
    // SetCanvasStyle(cRatioKaonPion, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioKaonPion = MakeRatio(gMYieldKaon[0], gMYieldPion[0], false);
    // TGraphErrors *gRatioKaonPion_sys = MakeRatio(gMYieldKaon[1], gMYieldPion[1], false);
    // gRatioKaonPion->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioKaonPion->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioKaonPion);
    // gRatioKaonPion->GetYaxis()->SetRangeUser(0.05, 0.25);
    // gRatioKaonPion->GetXaxis()->SetLimits(0, 27);
    // gRatioKaonPion->SetMarkerColor(kRed);
    // gRatioKaonPion->SetLineColor(kRed);
    // gRatioKaonPion->Draw("APE");
    // gRatioKaonPion_sys->SetFillStyle(0);
    // gRatioKaonPion_sys->SetLineColor(kRed);
    // gRatioKaonPion_sys->Draw("5 same");

    // TGraphErrors *gRatioKaonPion_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKaon_epos], gEPOS_Yield[9][kITY0][kPion_epos], true);
    // TGraphErrors *gRatioKaonPion_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKaon_epos], gEPOS_Yield[9][kITY80][kPion_epos], true);

    // gRatioKaonPion_IST9->SetLineStyle(2);
    // gRatioKaonPion_IST9->Draw("l same");
    // gRatioKaonPion_IST9_ITY80->SetLineColor(kBlue);
    // gRatioKaonPion_IST9_ITY80->SetLineStyle(2);
    // gRatioKaonPion_IST9_ITY80->Draw("l same");

    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioKaonPionModel = MakeRatio(gMYield[model][kKaon], gMYield[model][kPion], true);
    //     setStyle(gRatioKaonPionModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioKaonPionModel->Draw("l same");
    // }
    // latex.DrawLatex(0.28, 0.88, "#frac{K}{#pi}");
    // legendRatio->Draw();
    // legendRatio2->Draw();
    // if (isSavePlots)
    //{
    // cRatioKaonPion->SaveAs("Plots/Ratio_KaonPion_Run3.png");
    // }

    //    //=================================================
    //     // ====================Proton / Pion Ratio ==================
    //     //=================================================
    //     TCanvas *cRatioProtonPion = new TCanvas("cRatioProtonPion", "cRatioProtonPion", 720, 720);
    //     SetCanvasStyle(cRatioProtonPion, 0.15, 0.03, 0.03, 0.15);
    //     TGraphErrors *gRatioProtonPion = MakeRatio(gMYieldProton[0], gMYieldPion[0], false);
    //     TGraphErrors *gRatioProtonPion_sys = MakeRatio(gMYieldProton[1], gMYieldPion[1], false);
    //     gRatioProtonPion->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    //     gRatioProtonPion->GetYaxis()->SetTitle("dN/dy");
    //     SetGraphErrorStyle(gRatioProtonPion);
    //     gRatioProtonPion->GetYaxis()->SetRangeUser(0.02, 0.12);
    //     gRatioProtonPion->GetXaxis()->SetLimits(0, 27);
    //     gRatioProtonPion->SetMarkerColor(kRed);
    //     gRatioProtonPion->SetLineColor(kRed);
    //     gRatioProtonPion->Draw("APE");
    //     gRatioProtonPion_sys->SetFillStyle(0);
    //     gRatioProtonPion_sys->SetLineColor(kRed);
    //     gRatioProtonPion_sys->Draw("5 same");

    //     TGraphErrors *gRatioProtonPion_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kProton_epos], gEPOS_Yield[9][kITY0][kPion_epos], true);
    //     TGraphErrors *gRatioProtonPion_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kProton_epos], gEPOS_Yield[9][kITY80][kPion_epos], true);

    //     gRatioProtonPion_IST9->SetLineStyle(2);
    //     gRatioProtonPion_IST9->Draw("l same");
    //     gRatioProtonPion_IST9_ITY80->SetLineColor(kBlue);
    //     gRatioProtonPion_IST9_ITY80->SetLineStyle(2);
    //     gRatioProtonPion_IST9_ITY80->Draw("l same");

    //     for (auto model : modelsToPlot)
    //     {
    //         TGraphErrors *gRatioProtonPionModel = MakeRatio(gMYield[model][kProton], gMYield[model][kPion], true);
    //         setStyle(gRatioProtonPionModel, modelStyle[model].color, modelStyle[model].style);
    //         gRatioProtonPionModel->Draw("l same");
    //     }
    //     latex.DrawLatex(0.28, 0.88, "#frac{p}{#pi}");
    //     legendRatio->Draw();
    //     legendRatio2->Draw();
    // if (isSavePlots)
    // {
    //     cRatioProtonPion->SaveAs("Plots/Ratio_ProtonPion_Run3.png");
    //}
    //

    //====================================================
    // ==================Pion yeild======================
    //====================================================
    TCanvas *cPionYield = new TCanvas("cPionYield", "cPionYield", 720, 720);
    SetCanvasStyle(cPionYield, 0.15, 0.03, 0.03, 0.15);
    gMYieldPion[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldPion[0]->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    SetGraphErrorStyle(gMYieldPion[0]);
    gMYieldPion[0]->GetYaxis()->SetRangeUser(0.0, 13);
    gMYieldPion[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldPion[0]->SetMarkerColor(kRed);
    gMYieldPion[0]->SetLineColor(kRed);
    gMYieldPion[0]->Draw("APE");
    gMYieldPion[1]->SetFillStyle(0);
    gMYieldPion[1]->SetLineColor(kRed);
    gMYieldPion[1]->Draw("5 same");

    gMYieldPionEPOS_IST0->SetLineStyle(2);
    gMYieldPionEPOS_IST0->SetLineWidth(3);
    ScaleGraph(gMYieldPionEPOS_IST0, 0.5); // In model, Pi,K,p are not averaged.
    gMYieldPionEPOS_IST0->SetLineColor(kRed + 1);
    gMYieldPionEPOS_IST0->Draw("l same");

    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        setStyle(gPythiaYieldLocal[imodel][kPion_epos], modelStyle[imodel].color, modelStyle[imodel].style);
        ScaleGraph(gPythiaYieldLocal[imodel][kPion_epos], 0.5); // In pythia nothing is averaged, but in data Pi,K,p is averaged.
        gPythiaYieldLocal[imodel][kPion_epos]->Draw("l same");
    }
    TLegend *legTemp = new TLegend(0.2, 0.65, 0.55, 0.85);
    SetLegendStyle(legTemp);
    legTemp->SetTextSize(0.03);
    legTemp->AddEntry(gMYieldPion[0], "Data", "p");
    legTemp->AddEntry(gMYieldPionEPOS_IST0, "EPOS", "l");
    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        legTemp->AddEntry(gPythiaYieldLocal[imodel][kPion_epos], modelLabelLocal[imodel], "l");
    }
    legTemp->Draw();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.28, 0.88, "#pi^{#pm}");
    // if (isSavePlots)
    {
        cPionYield->SaveAs("Plots/PionYield_Run3.png");
    }

    //====================================================
    // ==================Kaon yeild======================
    //====================================================
    TCanvas *cKaonYield = new TCanvas("cKaonYield", "cKaonYield", 720, 720);
    SetCanvasStyle(cKaonYield, 0.15, 0.03, 0.03, 0.15);
    gMYieldKaon[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldKaon[0]->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    SetGraphErrorStyle(gMYieldKaon[0]);
    gMYieldKaon[0]->GetYaxis()->SetRangeUser(0.0, 2.2);
    gMYieldKaon[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldKaon[0]->SetMarkerColor(kRed);
    gMYieldKaon[0]->SetLineColor(kRed);
    gMYieldKaon[0]->Draw("APE");
    gMYieldKaon[1]->SetFillStyle(0);
    gMYieldKaon[1]->SetLineColor(kRed);
    gMYieldKaon[1]->Draw("5 same");

    gMYieldKaonEPOS_IST0->SetLineStyle(2);
    gMYieldKaonEPOS_IST0->SetLineWidth(3);
    ScaleGraph(gMYieldKaonEPOS_IST0, 0.5);
    gMYieldKaonEPOS_IST0->SetLineColor(kRed + 1);
    gMYieldKaonEPOS_IST0->Draw("l same");

    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        setStyle(gPythiaYieldLocal[imodel][kKaon_epos], modelStyle[imodel].color, modelStyle[imodel].style);
        ScaleGraph(gPythiaYieldLocal[imodel][kKaon_epos], 0.5); // In pythia nothing is averaged, but in data Pi,K,p is averaged.
        gPythiaYieldLocal[imodel][kKaon_epos]->Draw("l same");
    }
    legTemp->Draw();
    latex.DrawLatex(0.28, 0.88, "K^{#pm}");
    // if (isSavePlots)
    {
        cKaonYield->SaveAs("Plots/KaonYield_Run3.png");
    }

    //====================================================
    // ==================Proton yeild======================
    //====================================================
    TCanvas *cProtonYield = new TCanvas("cProtonYield", "cProtonYield", 720, 720);
    SetCanvasStyle(cProtonYield, 0.15, 0.03, 0.03, 0.15);
    gMYieldProton[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldProton[0]->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    SetGraphErrorStyle(gMYieldProton[0]);
    gMYieldProton[0]->GetYaxis()->SetRangeUser(0.0, 1.05);
    gMYieldProton[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldProton[0]->SetMarkerColor(kRed);
    gMYieldProton[0]->SetLineColor(kRed);
    gMYieldProton[0]->Draw("APE");
    gMYieldProton[1]->SetFillStyle(0);
    gMYieldProton[1]->SetLineColor(kRed);
    gMYieldProton[1]->Draw("5 same");

    gMYieldProtonEPOS_IST0->SetLineStyle(2);
    gMYieldProtonEPOS_IST0->SetLineWidth(3);
    ScaleGraph(gMYieldProtonEPOS_IST0, 0.5);
    gMYieldProtonEPOS_IST0->SetLineColor(kRed + 1);
    gMYieldProtonEPOS_IST0->Draw("l same");

    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        setStyle(gPythiaYieldLocal[imodel][kProton_epos], modelStyle[imodel].color, modelStyle[imodel].style);
        ScaleGraph(gPythiaYieldLocal[imodel][kProton_epos], 0.5); // In pythia nothing is averaged, but in data Pi,K,p is averaged.
        gPythiaYieldLocal[imodel][kProton_epos]->Draw("l same");
    }
    legTemp->Draw();
    latex.DrawLatex(0.28, 0.88, "p");
    // if (isSavePlots)
    {
        cProtonYield->SaveAs("Plots/ProtonYield_Run3.png");
    }

    //====================================================
    // ==================Phi yeild======================
    //====================================================
    TCanvas *cPhiYield = new TCanvas("cPhiYield", "cPhiYield", 720, 720);
    SetCanvasStyle(cPhiYield, 0.15, 0.03, 0.03, 0.15);
    gMYieldPhiRun2[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldPhiRun2[0]->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gMYieldPhiRun2[0]);
    gMYieldPhiRun2[0]->GetYaxis()->SetRangeUser(0.0, 0.25);
    gMYieldPhiRun2[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldPhiRun2[0]->SetMarkerColor(kRed);
    gMYieldPhiRun2[0]->SetLineColor(kRed);
    gMYieldPhiRun2[0]->Draw("APE");
    gMYieldPhiRun2[1]->SetFillStyle(0);
    gMYieldPhiRun2[1]->SetLineColor(kRed);
    gMYieldPhiRun2[1]->Draw("5 same");

    gEPOS_Yield[9][kITY0][kPhi_epos]->SetLineStyle(9);
    gEPOS_Yield[9][kITY0][kPhi_epos]->Draw("l same");
    gEPOS_Yield[9][kITY80][kPhi_epos]->SetLineColor(kBlue - 6);
    gEPOS_Yield[9][kITY80][kPhi_epos]->SetLineStyle(1);
    gEPOS_Yield[9][kITY80][kPhi_epos]->Draw("l same");

    TLegend *legTempPhi = new TLegend(0.2, 0.65, 0.55, 0.85);
    SetLegendStyle(legTempPhi);
    legTempPhi->SetTextSize(0.03);
    legTempPhi->AddEntry(gMYieldPhiRun2[0], "Data", "p");
    legTempPhi->AddEntry(gEPOS_Yield[9][kITY0][kPhi_epos], "EPOS UrQMD OFF", "l");
    legTempPhi->AddEntry(gEPOS_Yield[9][kITY80][kPhi_epos], "EPOS UrQMD ON", "l");
    legTempPhi->Draw();
    latex.DrawLatex(0.28, 0.88, "#phi");
    // if (isSavePlots)
    {
        cPhiYield->SaveAs("Plots/PhiYield_Run3.png");
    }

    // */
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
    TGraphErrors **ratio = new TGraphErrors *[2];

    for (int ig = 0; ig < 2; ig++)
    {
        if (!numerator[ig] || !denominator[ig])
        {
            cout << "Error: null graph" << endl;
            return nullptr;
        }

        if (numerator[ig]->GetN() != denominator[ig]->GetN())
        {
            cout << "Error: different number of points in numerator and denominator graphs" << endl;
            cout << "In numerator graph points " << numerator[ig]->GetN() << " and in denominator graph points " << denominator[ig]->GetN() << endl;
            return nullptr;
        }

        ratio[ig] = new TGraphErrors(numerator[ig]->GetN());

        for (int i = 0; i < numerator[ig]->GetN(); i++)
        {
            double xNum, yNum;
            double xDen, yDen;

            numerator[ig]->GetPoint(i, xNum, yNum);
            denominator[ig]->GetPoint(i, xDen, yDen);

            yDen *= CorrFactorDen;

            double yRatio = (yDen != 0.) ? yNum / yDen : 0.;
            ratio[ig]->SetPoint(i, xNum, yRatio);

            double xErr = numerator[ig]->GetErrorX(i);
            double numErr = numerator[ig]->GetErrorY(i);
            double denErr;

            if (ig == 1 && nGraphs >= 3)
                denErr = denominator[2]->GetErrorY(i);
            else
                denErr = denominator[ig]->GetErrorY(i);

            denErr *= CorrFactorDen;

            double yRatioErr = 0.;

            if (yDen != 0.)
            {
                yRatioErr = sqrt(pow(numErr / yDen, 2) + pow(yNum * denErr / (yDen * yDen), 2));
            }

            ratio[ig]->SetPointError(i, xErr, yRatioErr);
        }

        SetGraphErrorStyle(ratio[ig]);
    }

    return ratio;
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

void canvas_style(TCanvas *c, double pad1Size = 0.7, double pad2Size = 0.3, double leftMargin = 0.16, double rightMargin = 0.06, double topMargin = 0.02, double bottomMargin = 0.35)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    // pad2Size = 0.3; // Size of the first pad
    // pad1Size = 1 - pad2Size;

    pad1->SetPad(0, pad2Size, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, pad2Size); // x1, y1, x2, y2 (bottom pad)
    pad1->SetRightMargin(rightMargin);
    pad2->SetRightMargin(rightMargin);
    pad1->SetLeftMargin(leftMargin);
    pad2->SetLeftMargin(leftMargin);
    pad1->SetTopMargin(topMargin);
    pad2->SetTopMargin(0.001);
    pad2->SetBottomMargin(bottomMargin);
    pad1->SetBottomMargin(0.001);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
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
            ge->SetPointError(i,
                              ge->GetErrorX(i),
                              ge->GetErrorY(i) * scale);
        }
    }
}

void AddGraph(TGraph *gr1, TGraph *gr2)
{
    if (!gr1 || !gr2 || gr1->GetN() != gr2->GetN())
        return;

    for (int i = 0; i < gr1->GetN(); ++i)
    {
        double x1, y1, x2, y2;
        gr1->GetPoint(i, x1, y1);
        gr2->GetPoint(i, x2, y2);
        if (x1 == x2)
            gr1->SetPoint(i, x1, y1 + y2);
    }

    if (auto *ge1 = dynamic_cast<TGraphErrors *>(gr1))
    {
        auto *ge2 = dynamic_cast<TGraphErrors *>(gr2);
        if (!ge2)
            return;

        for (int i = 0; i < ge1->GetN(); ++i)
        {
            double ex1 = ge1->GetErrorX(i);
            double ey1 = ge1->GetErrorY(i);
            double ex2 = ge2->GetErrorX(i);
            double ey2 = ge2->GetErrorY(i);
            ge1->SetPointError(i, ex1, sqrt(ey1 * ey1 + ey2 * ey2));
        }
    }
}

TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance = 0.5, int noOfGraphs = 3)
{
    TGraphErrors **grCopy = new TGraphErrors *[noOfGraphs];

    for (int ig = 0; ig < noOfGraphs; ig++)
    {
        grCopy[ig] = (TGraphErrors *)gr[ig]->Clone();

        double yGivenMult = 1.0;
        double yGivenMultErr = 0.0;

        if (WhichMultPoint >= 0)
        {
            for (int i = 0; i < gr[ig]->GetN(); i++)
            {
                double x, y;
                gr[ig]->GetPoint(i, x, y);

                if (fabs(x - WhichMultPoint) < tolerance)
                {
                    yGivenMult = y;

                    if (ig == 1 && noOfGraphs >= 3)
                        yGivenMultErr = gr[2]->GetErrorY(i);
                    else
                        yGivenMultErr = gr[ig]->GetErrorY(i);

                    break;
                }
            }
        }
        else
        {
            double xMult1, xMult2, yMultPoint1, yMultLastPoint, yMultErr1, yMultErrLastPoint;
            gr[ig]->GetPoint(0, xMult1, yMultPoint1);
            gr[ig]->GetPoint(gr[ig]->GetN() - 1, xMult2, yMultLastPoint);
            yGivenMult = (xMult1 < xMult2) ? yMultPoint1 : yMultLastPoint;

            if (ig == 1 && noOfGraphs >= 3)
                yGivenMultErr = (xMult1 < xMult2) ? gr[2]->GetErrorY(0) : gr[2]->GetErrorY(gr[2]->GetN() - 1);
            else
                yGivenMultErr = (xMult1 < xMult2) ? gr[ig]->GetErrorY(0) : gr[ig]->GetErrorY(gr[ig]->GetN() - 1);
        }

        for (int i = 0; i < gr[ig]->GetN(); i++)
        {
            double x, y;
            gr[ig]->GetPoint(i, x, y);

            double xerr = gr[ig]->GetErrorX(i);

            double yerr;

            if (ig == 1 && noOfGraphs >= 3)
                yerr = gr[1]->GetErrorY(i);
            else
                yerr = gr[ig]->GetErrorY(i);

            // if (fabs(x - WhichMultPoint) < tolerance)
            // {
            //     grCopy[ig]->SetPoint(i, x, 1.0);
            //     grCopy[ig]->SetPointError(i, xerr, 0.0);
            //     continue;
            // }

            double ratio = y / yGivenMult;

            double ratioErr = sqrt(pow(yerr / yGivenMult, 2) + pow(y * yGivenMultErr / (yGivenMult * yGivenMult), 2));

            grCopy[ig]->SetPoint(i, x, ratio);
            grCopy[ig]->SetPointError(i, xerr, ratioErr);
        }
    }

    return grCopy;
}

int FindGraphXPoint(TGraphErrors *gr, double xValue)
{
    int xPoint = 0;
    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        if (abs(x - xValue) < 0.5)
        {
            // cout << "Found point at x = " << xValue << ": y = " << y << " +/- " << gr->GetErrorY(i) << endl;
            xPoint = i;
            break;
        }
    }
    return xPoint;
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

TGraphErrors **AddConstantToGraph(TGraphErrors **gr, double constant, int noOfGraphs = 2)
{
    TGraphErrors **grCopy = new TGraphErrors *[noOfGraphs];

    for (int ig = 0; ig < noOfGraphs; ig++)
    {
        if (gr[ig] == nullptr)
        {
            cout << "Error: null graph" << endl;
            return nullptr;
        }
        grCopy[ig] = (TGraphErrors *)gr[ig]->Clone();

        for (int i = 0; i < gr[ig]->GetN(); i++)
        {
            double x, y;
            gr[ig]->GetPoint(i, x, y);

            double yerr = gr[ig]->GetErrorY(i);
            double xerr = gr[ig]->GetErrorX(i);

            double newY = y + constant;

            grCopy[ig]->SetPoint(i, x, newY);
            grCopy[ig]->SetPointError(i, xerr, yerr);
        }
    }

    return grCopy;
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

/*
  //====================================================
    //  ================Kstar / (Kaon + Pion) Ratio ===================
    //====================================================
    TCanvas *cRatioKstarKaonPion = new TCanvas("cRatioKstarKaonPion", "cRatioKstarKaonPion", 720, 720);
    SetCanvasStyle(cRatioKstarKaonPion, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarKaonPion1 = MakeRatio(gMYieldKstar[0], gMYieldKaon[0], false);
    TGraphErrors *gRatioKstarKaonPion = MakeRatio(gRatioKstarKaonPion1, gMYieldPion[0], false);
    TGraphErrors *gRatioKstarKaonPion_sys1 = MakeRatio(gMYieldKstar[1], gMYieldKaon[1], false);
    TGraphErrors *gRatioKstarKaonPion_sys = MakeRatio(gRatioKstarKaonPion_sys1, gMYieldPion[1], false);
    gRatioKstarKaonPion->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarKaonPion->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarKaonPion);
    gRatioKstarKaonPion->GetYaxis()->SetRangeUser(0.09, 0.39);
    gRatioKstarKaonPion->GetXaxis()->SetLimits(0, 27);
    gRatioKstarKaonPion->SetMarkerColor(kRed);
    gRatioKstarKaonPion->SetLineColor(kRed);
    gRatioKstarKaonPion->Draw("APE");
    gRatioKstarKaonPion_sys->SetFillStyle(0);
    gRatioKstarKaonPion_sys->SetLineColor(kRed);
    gRatioKstarKaonPion_sys->Draw("5 same");

    TGraphErrors *gRatioKstarKaonPion_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gMYieldKaonEPOS_IST0, true);
    gRatioKstarKaonPion_IST9 = MakeRatio(gRatioKstarKaonPion_IST9, gMYieldPionEPOS_IST0, true);
    TGraphErrors *gRatioKstarKaonPion_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gMYieldKaonEPOS_IST0, true);
    gRatioKstarKaonPion_IST9_ITY80 = MakeRatio(gRatioKstarKaonPion_IST9_ITY80, gMYieldPionEPOS_IST0, true);

    gRatioKstarKaonPion_IST9->SetLineStyle(2);
    gRatioKstarKaonPion_IST9->Draw("l same");
    gRatioKstarKaonPion_IST9_ITY80->SetLineColor(kBlue);
    gRatioKstarKaonPion_IST9_ITY80->SetLineStyle(2);
    gRatioKstarKaonPion_IST9_ITY80->Draw("l same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarKaonPionModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKaon], true);
        gRatioKstarKaonPionModel = MakeRatio(gRatioKstarKaonPionModel, gMYield[model][kPion], true);
        setStyle(gRatioKstarKaonPionModel, modelStyle[model].color, modelStyle[model].style);
        gRatioKstarKaonPionModel->Draw("l same");
    }

    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{(K * #pi)}");
    legendRatio->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
    cRatioKstarKaonPion->SaveAs("Plots/Ratio_KstarKaonPion_Run3.png");
    }

        //====================================================================
    //============Xistar dN/dy plot===========================
    //====================================================================
    TCanvas *cXiStarYield = new TCanvas("cXiStarYield", "cXiStarYield", 720, 720);
    SetCanvasStyle(cXiStarYield, 0.15, 0.03, 0.03, 0.15);
    gMYieldXiStar[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldXiStar[0]->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gMYieldXiStar[0]);
    // gMYieldXiStar[0]->GetYaxis()->SetRangeUser(0.0005, 0.01);
    gMYieldXiStar[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldXiStar[0]->SetMarkerColor(kRed);
    gMYieldXiStar[0]->SetLineColor(kRed);
    gMYieldXiStar[0]->Draw("APE");
    gMYieldXiStar[1]->SetFillStyle(0);
    gMYieldXiStar[1]->SetLineColor(kRed);
    gMYieldXiStar[1]->Draw("5 same");

       //===================================================
    // ======================Phi/ Pi Ratio===================
    //====================================================

    // ******** There is mistake in either 7 TeV or 13 TeV, because one of them is exactly off by a factor of 2 *********
    TCanvas *cRatioPhiPion = new TCanvas("cRatioPhiPion", "cRatioPhiPion", 720, 720);
    SetCanvasStyle(cRatioPhiPion, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioPhiPion = MakeRatio(gMYieldPhi[0], gMYieldPion[0], false);
    TGraphErrors *gRatioPhiPion_sys = MakeRatio(gMYieldPhi[1], gMYieldPion[1], false);
    gRatioPhiPion->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioPhiPion->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioPhiPion);
    gRatioPhiPion->GetYaxis()->SetRangeUser(0.002, 0.032);
    gRatioPhiPion->GetXaxis()->SetLimits(0, 27);
    gRatioPhiPion->SetMarkerColor(kRed);
    gRatioPhiPion->SetLineColor(kRed);
    gRatioPhiPion->Draw("APE");
    gRatioPhiPion_sys->SetFillStyle(0);
    gRatioPhiPion_sys->SetLineColor(kRed);
    gRatioPhiPion_sys->Draw("5 same");
    gPhiPiRatio_13TeV[0]->SetMarkerStyle(21);
    gPhiPiRatio_13TeV[0]->SetMarkerColor(kBlue);
    gPhiPiRatio_13TeV[0]->SetLineColor(kBlue);
    gPhiPiRatio_13TeV[0]->Draw("P same");
    gPhiPiRatio_13TeV[1]->SetLineColor(kBlue);
    gPhiPiRatio_13TeV[1]->SetFillStyle(0);
    gPhiPiRatio_13TeV[1]->Draw("5 same");
    // ScaleGraph(gPhiPiRatio_7TeV[0], 0.5); // In 7 TeV it was just K*0 + K*0bar and not its average
    // ScaleGraph(gPhiPiRatio_7TeV[1], 0.5);
    gPhiPiRatio_7TeV[0]->SetMarkerStyle(22);
    gPhiPiRatio_7TeV[0]->SetMarkerColor(kGreen + 2);
    gPhiPiRatio_7TeV[0]->SetLineColor(kGreen + 2);
    gPhiPiRatio_7TeV[0]->Draw("P same");
    gPhiPiRatio_7TeV[1]->SetLineColor(kGreen + 2);
    gPhiPiRatio_7TeV[1]->SetFillStyle(0);
    gPhiPiRatio_7TeV[1]->Draw("5 same");

    TGraphErrors *gRatioPhiPion_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gMYieldPionEPOS_IST0, true);
    TGraphErrors *gRatioPhiPion_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gMYieldPionEPOS_IST0, true);

    gRatioPhiPion_IST9->SetLineStyle(2);
    gRatioPhiPion_IST9->Draw("l same");
    gRatioPhiPion_IST9_ITY80->SetLineColor(kBlue);
    gRatioPhiPion_IST9_ITY80->SetLineStyle(2);
    gRatioPhiPion_IST9_ITY80->Draw("l same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiPionModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kPion], true);
        setStyle(gRatioPhiPionModel, modelStyle[model].color, modelStyle[model].style);
        gRatioPhiPionModel->Draw("l same");
    }
    latex.DrawLatex(0.28, 0.88, "#frac{#phi}{#pi}");
    legendRatio3->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
      cRatioPhiPion->SaveAs("Plots/Ratio_PhiPion_Run3.png");
    }

    //====================================================
    //   ================Phi/K Ratio===================
    //====================================================
    TCanvas *cRatioPhiKaon = new TCanvas("cRatioPhiKaon", "cRatioPhiKaon", 720, 720);
    SetCanvasStyle(cRatioPhiKaon, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioPhiKaon = MakeRatio(gMYieldPhi[0], gMYieldKaon[0], false);
    TGraphErrors *gRatioPhiKaon_sys = MakeRatio(gMYieldPhi[1], gMYieldKaon[1], false);
    gRatioPhiKaon->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioPhiKaon->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioPhiKaon);
    gRatioPhiKaon->GetYaxis()->SetRangeUser(0.0, 0.32);
    gRatioPhiKaon->GetXaxis()->SetLimits(0, 27);
    gRatioPhiKaon->SetMarkerColor(kRed);
    gRatioPhiKaon->SetLineColor(kRed);
    gRatioPhiKaon->Draw("APE");
    gRatioPhiKaon_sys->SetFillStyle(0);
    gRatioPhiKaon_sys->SetLineColor(kRed);
    gRatioPhiKaon_sys->Draw("5 same");
    gPhiKaRatio_13TeV[0]->SetMarkerStyle(21);
    gPhiKaRatio_13TeV[0]->SetMarkerColor(kBlue);
    gPhiKaRatio_13TeV[0]->SetLineColor(kBlue);
    gPhiKaRatio_13TeV[0]->Draw("P same");
    gPhiKaRatio_13TeV[1]->SetLineColor(kBlue);
    gPhiKaRatio_13TeV[1]->SetFillStyle(0);
    gPhiKaRatio_13TeV[1]->Draw("5 same");

    ////K = (K^+ + K^-) /2, so no need for correction.
    TGraphErrors *gRatioPhiKaon_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gMYieldKaonEPOS_IST0, true);
    TGraphErrors *gRatioPhiKaon_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gMYieldKaonEPOS_IST0, true);

    gRatioPhiKaon_IST9->SetLineStyle(2);
    gRatioPhiKaon_IST9->Draw("l same");
    gRatioPhiKaon_IST9_ITY80->SetLineColor(kBlue);
    gRatioPhiKaon_IST9_ITY80->SetLineStyle(2);
    gRatioPhiKaon_IST9_ITY80->Draw("l same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiKaonModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kKaon], true);
        setStyle(gRatioPhiKaonModel, modelStyle[model].color, modelStyle[model].style);
        gRatioPhiKaonModel->Draw("l same");
    }

    latex.DrawLatex(0.28, 0.88, "#frac{#phi}{K}");
    legendRatio->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
    cRatioPhiKaon->SaveAs("Plots/Ratio_PhiKaon_Run3.png");
    }

    //===================================================
    //   ================Phi/Kshort Ratio===================
    //====================================================

    TCanvas *cRatioPhiKshort = new TCanvas("cRatioPhiKshort", "cRatioPhiKshort", 720, 720);
    SetCanvasStyle(cRatioPhiKshort, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioPhiKshort = MakeRatio(gMYieldPhi[0], gMYieldKshort[0], false);
    TGraphErrors *gRatioPhiKshort_sys = MakeRatio(gMYieldPhi[1], gMYieldKshort[1], false);
    gRatioPhiKshort->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioPhiKshort->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioPhiKshort);
    gRatioPhiKshort->GetYaxis()->SetRangeUser(0.0, 0.32);
    gRatioPhiKshort->GetXaxis()->SetLimits(0, 27);
    gRatioPhiKshort->SetMarkerColor(kRed);
    gRatioPhiKshort->SetLineColor(kRed);
    gRatioPhiKshort->Draw("APE");
    gRatioPhiKshort_sys->SetFillStyle(0);
    gRatioPhiKshort_sys->SetLineColor(kRed);
    gRatioPhiKshort_sys->Draw("5 same");
    gPhiKshortRatio_13TeV[0]->SetMarkerStyle(21);
    gPhiKshortRatio_13TeV[0]->SetMarkerColor(kBlue);
    gPhiKshortRatio_13TeV[0]->SetLineColor(kBlue);
    gPhiKshortRatio_13TeV[0]->Draw("P same");
    gPhiKshortRatio_13TeV[1]->SetLineColor(kBlue);
    gPhiKshortRatio_13TeV[1]->SetFillStyle(0);
    gPhiKshortRatio_13TeV[1]->Draw("5 same");

    TGraphErrors *gRatioPhiKshort_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gMYieldKshortEPOS_IST0, true, 0.5);
    TGraphErrors *gRatioPhiKshort_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gMYieldKshortEPOS_IST0, true, 0.5);

    gRatioPhiKshort_IST9->SetLineStyle(2);
    gRatioPhiKshort_IST9->Draw("l same");
    gRatioPhiKshort_IST9_ITY80->SetLineColor(kBlue);
    gRatioPhiKshort_IST9_ITY80->SetLineStyle(2);
    gRatioPhiKshort_IST9_ITY80->Draw("l same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiKshortModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kKshort], true);
        setStyle(gRatioPhiKshortModel, modelStyle[model].color, modelStyle[model].style);
        gRatioPhiKshortModel->Draw("l same");
    }
    latex.DrawLatex(0.28, 0.88, "#frac{#phi}{K_{S}^{0}}");
    legendRatio->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
    cRatioPhiKshort->SaveAs("Plots/Ratio_PhiKshort_Run3.png");
    }

    //===================================================
    //  ==================Phi/ Proton Ratio===================
    //====================================================
    TCanvas *cRatioPhiProton = new TCanvas("cRatioPhiProton", "cRatioPhiProton", 720, 720);
    SetCanvasStyle(cRatioPhiProton, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioPhiProton = MakeRatio(gMYieldPhi[0], gMYieldProton[0], false, 2);
    TGraphErrors *gRatioPhiProton_sys = MakeRatio(gMYieldPhi[1], gMYieldProton[1], false, 2);
    gRatioPhiProton->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioPhiProton->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioPhiProton);
    gRatioPhiProton->GetYaxis()->SetRangeUser(0.0, 0.259);
    gRatioPhiProton->GetXaxis()->SetLimits(0, 27);
    gRatioPhiProton->SetMarkerColor(kRed);
    gRatioPhiProton->SetLineColor(kRed);
    gRatioPhiProton->Draw("APE");
    gRatioPhiProton_sys->SetFillStyle(0);
    gRatioPhiProton_sys->SetLineColor(kRed);
    gRatioPhiProton_sys->Draw("5 same");
    gPhiPrRatio_13TeV[0]->SetMarkerStyle(21);
    gPhiPrRatio_13TeV[0]->SetMarkerColor(kBlue);
    gPhiPrRatio_13TeV[0]->SetLineColor(kBlue);
    gPhiPrRatio_13TeV[0]->Draw("P same");
    gPhiPrRatio_13TeV[1]->SetLineColor(kBlue);
    gPhiPrRatio_13TeV[1]->SetFillStyle(0);
    gPhiPrRatio_13TeV[1]->Draw("5 same");

    TGraphErrors *gRatioPhiProton_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gMYieldProtonEPOS_IST0, true);
    TGraphErrors *gRatioPhiProton_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gMYieldProtonEPOS_IST0, true);

    gRatioPhiProton_IST9->SetLineStyle(2);
    gRatioPhiProton_IST9->Draw("l same");
    gRatioPhiProton_IST9_ITY80->SetLineColor(kBlue);
    gRatioPhiProton_IST9_ITY80->SetLineStyle(2);
    gRatioPhiProton_IST9_ITY80->Draw("l same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiProtonModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kProton], true);
        setStyle(gRatioPhiProtonModel, modelStyle[model].color, modelStyle[model].style);
        gRatioPhiProtonModel->Draw("l same");
    }
    latex.DrawLatex(0.28, 0.88, "#frac{#phi}{p}");
    legendRatio->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
        cRatioPhiProton->SaveAs("Plots/Ratio_PhiProton_Run3.png");
    }

    //================================================
    // ==================Phi/ Xi* Ratio===================
    //====================================================

    TCanvas *cRatioPhiXiStar = new TCanvas("cRatioPhiXiStar", "cRatioPhiXiStar", 720, 720);
    SetCanvasStyle(cRatioPhiXiStar, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioPhiXiStar = MakeRatio(gMYieldPhi[0], gMYieldXiStar[0], false);
    TGraphErrors *gRatioPhiXiStar_sys = MakeRatio(gMYieldPhi[1], gMYieldXiStar[1], false);
    gRatioPhiXiStar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioPhiXiStar->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioPhiXiStar);
    gRatioPhiXiStar->GetYaxis()->SetRangeUser(0.0, 209);
    gRatioPhiXiStar->GetXaxis()->SetLimits(0, 27);
    gRatioPhiXiStar->Draw("APE");
    gRatioPhiXiStar_sys->SetFillStyle(0);
    gRatioPhiXiStar_sys->Draw("5 same");

    TGraphErrors *gRatioPhiXiStar_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gEPOS_Yield[9][kITY0][kXi1530_epos], true);
    TGraphErrors *gRatioPhiXiStar_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gEPOS_Yield[9][kITY80][kXi1530_epos], true);
    TGraphErrors *gRatioPhiXiStar_IST9_ITY81 = MakeRatio(gEPOS_Yield[9][kITY81][kPhi_epos], gEPOS_Yield[9][kITY81][kXi1530_epos], true);
    TGraphErrors *gRatioPhiXiStar_IST6 = MakeRatio(gEPOS_Yield[6][kITY0][kPhi_epos], gEPOS_Yield[6][kITY0][kXi1530_epos], true);

    gRatioPhiXiStar_IST9->SetLineStyle(2);
    gRatioPhiXiStar_IST9->Draw("l same");
    gRatioPhiXiStar_IST9_ITY80->SetLineColor(kBlue);
    gRatioPhiXiStar_IST9_ITY80->SetLineStyle(2);
    gRatioPhiXiStar_IST9_ITY80->Draw("l same");
    gRatioPhiXiStar_IST9_ITY81->SetLineColor(kGreen + 2);
    gRatioPhiXiStar_IST9_ITY81->SetLineStyle(2);
    gRatioPhiXiStar_IST9_ITY81->Draw("l same");
    gRatioPhiXiStar_IST6->SetLineColor(kMagenta);
    gRatioPhiXiStar_IST6->SetLineStyle(5);
    // gRatioPhiXiStar_IST6->Draw("l same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiXiStarModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kXi1530], true);
        setStyle(gRatioPhiXiStarModel, modelStyle[model].color, modelStyle[model].style);
        gRatioPhiXiStarModel->Draw("l same");
    }
    latex.DrawLatex(0.28, 0.88, "#frac{#phi}{2#Xi^{*0}}");
    legendRatio->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
    cRatioPhiXiStar->SaveAs("Plots/Ratio_PhiXiStar_Run3.png");
    }
*/