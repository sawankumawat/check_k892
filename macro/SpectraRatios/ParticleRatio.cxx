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
TGraphErrors *GetGraph(TFile *f, const string &name);
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void setStyle(TGraphErrors *gr, int color, int style)
{
    gr->SetLineColor(color);
    gr->SetLineStyle(style);
}

void ScaleGraph(TGraph *gr, double scale);
void AddGraph(TGraph *gr1, TGraph *gr2);
TGraphErrors *DivideByMult(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5);
int FindGraphXPoint(TGraphErrors *gr, double xValue);
void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax);

void ParticleRatio()
{
    bool isSavePlots = true;
    string KstarPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/";
    TFile *fKstar = new TFile((KstarPath + "Results.root").c_str(), "read");
    if (fKstar->IsZombie())
    {
        cout << "Error: Kstar file not found" << endl;
        return;
    }

    TFile *fPion = new TFile("PiKp_Run3_Results/Sawan/Pi_results.root", "read");
    TFile *fProton = new TFile("PiKp_Run3_Results/Sawan/Pr_results.root", "read");
    TFile *fKaon = new TFile("PiKp_Run3_Results/Sawan/Ka_results.root", "read");

    TFile *fPhi = new TFile("ConversionCodes/pp13TeVALICE.root", "read");
    TFile *fChKstar = new TFile("ChargedKstar_Results/Sawan/ResultsChargedKstar.root", "read");
    TFile *fXiStar = new TFile("XiStar_Results/Sawan/ResultsXiStar.root", "read");
    TFile *fKshort = new TFile("K0s_Run3_Results/Sawan/ResultsK0s.root", "read");
    TFile *fLambda1520 = new TFile("LambdaRun3/Sawan/ResultsLambda1520.root", "read");
    TFile *fRho = new TFile("Rho_Run3_Results/Sawan/ResultsRho.root", "read");

    if (fPion->IsZombie() || fProton->IsZombie() || fKaon->IsZombie() || fPhi->IsZombie() || fChKstar->IsZombie() || fXiStar->IsZombie() || fKshort->IsZombie() || fLambda1520->IsZombie() || fRho->IsZombie())
    {
        cout << "Error: Particle files not found" << endl;
        return;
    }

    // FIXME: Add 3rd column for uncorrelated systematic error.
    TGraphErrors *gMPtKstar[2], *gMYieldKstar[2], *gMPtPhi[2], *gMPtPion[2], *gMPtProton[2], *gMPtKaon[2], *gMPtChKstar[2], *gMPtXiStar[2], *gMPtKshort[2], *gMPtLambda1520[2], *gMYieldPhi[2], *gMYieldPion[2], *gMYieldProton[2], *gMYieldKaon[2], *gMYieldChKstar[2], *gMYieldXiStar[2], *gMYieldKshort[2], *gMYieldLambda1520[2], *gMYieldRho[2], *gMPtRho[2];

    for (int i = 0; i < 2; i++)
    {
        gMPtKstar[i] = GetGraph(fKstar, Form("gMeanpTRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldKstar[i] = GetGraph(fKstar, Form("gMeanYieldRun3%s", (i == 0 ? "_stat" : "_sys")));

        gMPtPhi[i] = GetGraph(fPhi, Form("gPhi_MeanpT%s", (i == 0 ? "_stat" : "_sys")));
        gMPtPion[i] = GetGraph(fPion, Form("gMeanpTRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMPtProton[i] = GetGraph(fProton, Form("gMeanpTRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMPtKaon[i] = GetGraph(fKaon, Form("gMeanpTRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMPtChKstar[i] = GetGraph(fChKstar, Form("gMeanpTRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMPtXiStar[i] = GetGraph(fXiStar, Form("gMeanpTRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMPtKshort[i] = GetGraph(fKshort, Form("gMeanpTRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMPtLambda1520[i] = GetGraph(fLambda1520, Form("gMeanpTRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMPtRho[i] = GetGraph(fRho, Form("gMeanpTRun3%s", (i == 0 ? "_stat" : "_sys")));

        gMYieldPhi[i] = GetGraph(fPhi, Form("gPhi_MeanYield%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldPion[i] = GetGraph(fPion, Form("gMeanYieldRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldProton[i] = GetGraph(fProton, Form("gMeanYieldRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldKaon[i] = GetGraph(fKaon, Form("gMeanYieldRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldChKstar[i] = GetGraph(fChKstar, Form("gMeanYieldRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldXiStar[i] = GetGraph(fXiStar, Form("gMeanYieldRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldKshort[i] = GetGraph(fKshort, Form("gMeanYieldRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldLambda1520[i] = GetGraph(fLambda1520, Form("gMeanYieldRun3%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldRho[i] = GetGraph(fRho, Form("gMeanYieldRun3%s", (i == 0 ? "_stat" : "_sys")));
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
    TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA_INELgt0Correct.root", "read");
    // TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA_ptCut_FinerBins.root", "read"); //Used for results
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
    if (fPythiaMonash->IsZombie() || fPythiaMonashNoCR->IsZombie() || fPythiaShoving->IsZombie() || fPythiaMonashRescattering->IsZombie())
    {
        cout << "Error: Pythia local files not found" << endl;
        return;
    }
    enum PythiaModel
    {
        kPythiaMonashLocal,
        kPythiaMonashNoCRLocal,
        kPythiaShovingLocal,
        kNPythiaModels
    };

    const char *modelLabelLocal[kNPythiaModels] = {
        "Pythia Monash",
        "Pythia Monash No CR",
        "Pythia Shoving"};
    vector<TFile *> fPythiaModels = {fPythiaMonash, fPythiaMonashNoCR, fPythiaShoving};
    int lineStylesPythia[kNPythiaModels] = {1, 1, 2};
    TGraphErrors *gPythiaYieldLocal[kNPythiaModels][kNParticles_epos];
    TGraphErrors *gPythiaMeanPtLocal[kNPythiaModels][kNParticles_epos];

    int colorsPythia[kNPythiaModels] = {kCyan + 1, kYellow + 1, kMagenta};

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
        {kBlue, 2},      // Pythia CR
        {kRed, 4},       // Pythia Monash
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
    gMYieldKstar[0]->SetMarkerColor(kRed);
    gMYieldKstar[0]->SetLineColor(kRed);
    gMYieldKstar[0]->Draw("APE");
    gMYieldKstar[1]->SetFillStyle(0);
    gMYieldKstar[1]->SetLineColor(kRed);
    gMYieldKstar[1]->SetLineWidth(3);
    gMYieldKstar[1]->Draw("5 same");
    gKstar_Yield_13TeV[0]->SetMarkerStyle(21);
    gKstar_Yield_13TeV[0]->SetMarkerColor(kBlue);
    gKstar_Yield_13TeV[0]->SetLineColor(kBlue);
    gKstar_Yield_13TeV[0]->SetLineWidth(3);
    gKstar_Yield_13TeV[0]->Draw("P same");
    gKstar_Yield_13TeV[1]->SetLineColor(kBlue);
    gKstar_Yield_13TeV[1]->SetFillStyle(0);
    gKstar_Yield_13TeV[1]->SetLineWidth(3);
    gKstar_Yield_13TeV[1]->Draw("5 same");

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

    ScaleGraph(gPythiaYieldLocal[kPythiaMonashLocal][kKstar_epos], 0.5); // In pythia nothing is averaged. So we have to divide by 2 here.
    ScaleGraph(gPythiaYieldLocal[kPythiaMonashNoCRLocal][kKstar_epos], 0.5);
    ScaleGraph(gPythiaYieldLocal[kPythiaShovingLocal][kKstar_epos], 0.5);
    gPythiaYieldLocal[kPythiaMonashLocal][kKstar_epos]->Draw("l same");
    gPythiaYieldLocal[kPythiaMonashNoCRLocal][kKstar_epos]->Draw("l same");
    gPythiaYieldLocal[kPythiaShovingLocal][kKstar_epos]->Draw("l same");

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
    //  ================<pT> K*======================
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

    gPythiaMeanPtLocal[kPythiaMonashLocal][kKstar_epos]->Draw("l same");
    gPythiaMeanPtLocal[kPythiaMonashNoCRLocal][kKstar_epos]->Draw("l same");
    gPythiaMeanPtLocal[kPythiaShovingLocal][kKstar_epos]->Draw("l same");

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
    if (isSavePlots)
    {
        cMeanPtKstar->SaveAs("Plots/MeanPt_Kstar_EPOS_UrQMDON.png");
    }

    //===================================================
    //  ================<pT> all particles======================
    //===================================================
    TCanvas *cMeanPtAll = new TCanvas("cMeanPtAll", "cMeanPtAll", 720, 720);
    SetCanvasStyle(cMeanPtAll, 0.15, 0.03, 0.03, 0.15);
    gMPtKstar[0]->SetMaximum(2.09);
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
    gMPtKshort[0]->SetMarkerColor(kCyan + 1);
    gMPtKshort[0]->SetLineColor(kCyan + 1);
    gMPtKshort[0]->Draw("PE same");
    gMPtKshort[1]->SetLineColor(kCyan + 1);
    gMPtKshort[1]->SetFillStyle(0);
    gMPtKshort[1]->Draw("5 same");
    gMPtChKstar[0]->SetMarkerColor(kOrange + 1);
    gMPtChKstar[0]->SetLineColor(kOrange + 1);
    gMPtChKstar[0]->Draw("PE same");
    gMPtChKstar[1]->SetLineColor(kOrange + 1);
    gMPtChKstar[1]->SetFillStyle(0);
    gMPtChKstar[1]->Draw("5 same");
    gMPtXiStar[0]->SetMarkerColor(kViolet + 2);
    gMPtXiStar[0]->SetLineColor(kViolet + 2);
    gMPtXiStar[0]->Draw("PE same");
    gMPtXiStar[1]->SetLineColor(kViolet + 2);
    gMPtXiStar[1]->SetFillStyle(0);
    gMPtXiStar[1]->Draw("5 same");
    // gMPtPhi[0]->SetMarkerColor(kOrange + 1);
    // gMPtPhi[0]->SetLineColor(kOrange + 1);
    // gMPtPhi[0]->Draw("PE same");
    // gMPtPhi[1]->SetLineColor(kOrange + 1);
    // gMPtPhi[1]->SetFillStyle(0);
    // gMPtPhi[1]->Draw("5 same");
    gMPtLambda1520[0]->SetMarkerColor(kBrown);
    gMPtLambda1520[0]->SetLineColor(kBrown);
    gMPtLambda1520[0]->Draw("PE same");
    gMPtLambda1520[1]->SetLineColor(kBrown);
    gMPtLambda1520[1]->SetFillStyle(0);
    gMPtLambda1520[1]->Draw("5 same");
    gMPtRho[0]->SetMarkerColor(kGray + 2);
    gMPtRho[0]->SetLineColor(kGray + 2);
    gMPtRho[0]->Draw("PE same");
    gMPtRho[1]->SetLineColor(kGray + 2);
    gMPtRho[1]->SetFillStyle(0);
    gMPtRho[1]->Draw("5 same");

    TLegend *legendMeanPt = new TLegend(0.2, 0.77, 0.7, 0.92);
    SetLegendStyle(legendMeanPt);
    legendMeanPt->SetTextSize(0.027);
    legendMeanPt->SetNColumns(4);
    legendMeanPt->AddEntry(gMPtPion[0], "Pion", "P");
    legendMeanPt->AddEntry(gMPtKaon[0], "Kaon", "P");
    legendMeanPt->AddEntry(gMPtKshort[0], "K_{S}^{0}", "P");
    legendMeanPt->AddEntry(gMPtChKstar[0], "K*^{#pm}", "P");
    legendMeanPt->AddEntry(gMPtKstar[0], "K* (892)^{0}", "P");
    legendMeanPt->AddEntry(gMPtProton[0], "Proton", "P");
    legendMeanPt->AddEntry(gMPtXiStar[0], "#Xi(1530)", "P");
    // legendMeanPt->AddEntry(gMPtPhi[0], "#phi (1020)", "P");
    legendMeanPt->AddEntry(gMPtLambda1520[0], "#Lambda(1520)", "P");
    legendMeanPt->AddEntry(gMPtRho[0], "#rho (770)", "P");
    legendMeanPt->Draw();
    if (isSavePlots)
    {
        cMeanPtAll->SaveAs("Plots/MeanPt_AllParticles_Run3.png");
    }

    //===================================================
    //  =========<pT> all particles by their masses======================
    //===================================================
    double PDGMasses[] = {0.8956, 1.0195, 0.1396, 0.494, 0.938, 0.498, 0.8955, 1.5318, 1.5192, 0.7753}; // K*0, phi, pion, kaon, proton, Kshort, charged K*, Xi(1530), Lambda(1520), Rho(770)
    double PDGMassesErrors[] = {0.0002, 0.000016, 0.00000018, 0.000015, 29e-9, 0.000013, 0.0002, 0.00034, 0.00019, 0.0002};
    TCanvas *cMeanPtMass = new TCanvas("cMeanPtMass", "cMeanPtMass", 720, 720);
    SetCanvasStyle(cMeanPtMass, 0.15, 0.03, 0.03, 0.15);
    for (int i = 0; i < 2; i++)
    {
        ScaleGraph(gMPtKstar[i], 1.0 / PDGMasses[0]);
        ScaleGraph(gMPtPion[i], 1.0 / PDGMasses[2]);
        ScaleGraph(gMPtKaon[i], 1.0 / PDGMasses[3]);
        ScaleGraph(gMPtProton[i], 1.0 / PDGMasses[4]);
        ScaleGraph(gMPtKshort[i], 1.0 / PDGMasses[5]);
        ScaleGraph(gMPtChKstar[i], 1.0 / PDGMasses[6]);
        ScaleGraph(gMPtXiStar[i], 1.0 / PDGMasses[7]);
        ScaleGraph(gMPtLambda1520[i], 1.0 / PDGMasses[8]);
        ScaleGraph(gMPtRho[i], 1.0 / PDGMasses[9]);
        // ScaleGraph(gMPtPhi[i], 1.0 / PDGMasses[1]);
    }
    gMPtKstar[0]->SetMaximum(4.89);
    gMPtKstar[0]->SetMinimum(0.32 / PDGMasses[0]);
    gMPtKstar[0]->GetYaxis()->SetTitle("<#it{p}_{T}>/m (GeV/#it{c}^{2})");
    gMPtKstar[0]->Draw("APE");
    gMPtPion[0]->Draw("PE same");
    gMPtKaon[0]->Draw("PE same");
    gMPtProton[0]->Draw("PE same");
    gMPtKshort[0]->Draw("PE same");
    gMPtChKstar[0]->Draw("PE same");
    gMPtXiStar[0]->Draw("PE same");
    gMPtLambda1520[0]->Draw("PE same");
    // gMPtPhi[0]->Draw("PE same");
    gMPtRho[0]->Draw("PE same");

    gMPtKstar[1]->Draw("5 same");
    gMPtPion[1]->Draw("5 same");
    gMPtKaon[1]->Draw("5 same");
    gMPtProton[1]->Draw("5 same");
    gMPtKshort[1]->Draw("5 same");
    gMPtChKstar[1]->Draw("5 same");
    gMPtXiStar[1]->Draw("5 same");
    // gMPtPhi[1]->Draw("5 same");
    gMPtLambda1520[1]->Draw("5 same");
    gMPtRho[1]->Draw("5 same");

    legendMeanPt->Draw();

    if (isSavePlots)
    {
        cMeanPtMass->SaveAs("Plots/MeanPt_MassScaled_Run3.png");
    }

    cout << "Total bins for K*0 mean pT: " << gMPtXiStar[0]->GetN() << endl;

    for (int iBins = 0; iBins < gMPtXiStar[0]->GetN(); iBins++)
    {
        double x, y;
        gMPtXiStar[0]->GetPoint(iBins, x, y);
        cout << "Multiplicity: " << x << ", Mean pT/m for K*0: " << y << endl;
    }

    //==================================================================================
    // ===================<pt>_mult / <pt>_lowestMult vs <dNch/deta> =====================
    //==================================================================================
    TCanvas *cMeanPtRatio = new TCanvas("cMeanPtRatio", "cMeanPtRatio", 720, 720);
    SetCanvasStyle(cMeanPtRatio, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gMPtKstarLMRatio[2], *gMPtChKstarLMRatio[2], *gMPtPionLMRatio[2], *gMPtKaonLMRatio[2], *gMPtProtonLMRatio[2], *gMPtKshortLMRatio[2], *gMPtXiStarLMRatio[2], *gMPtLambda1520LMRatio[2], *gMPtRhoLMRatio[2];
    for (int i = 0; i < 2; i++)
    {
        gMPtKstarLMRatio[i] = (TGraphErrors *)gMPtKstar[i]->Clone(Form("gMPtKstarLMRatio_%d", i));
        gMPtPionLMRatio[i] = (TGraphErrors *)gMPtPion[i]->Clone(Form("gMPtPionLMRatio_%d", i));
        gMPtKaonLMRatio[i] = (TGraphErrors *)gMPtKaon[i]->Clone(Form("gMPtKaonLMRatio_%d", i));
        gMPtProtonLMRatio[i] = (TGraphErrors *)gMPtProton[i]->Clone(Form("gMPtProtonLMRatio_%d", i));
        gMPtKshortLMRatio[i] = (TGraphErrors *)gMPtKshort[i]->Clone(Form("gMPtKshortLMRatio_%d", i));
        gMPtChKstarLMRatio[i] = (TGraphErrors *)gMPtChKstar[i]->Clone(Form("gMPtChKstarLMRatio_%d", i));
        gMPtXiStarLMRatio[i] = (TGraphErrors *)gMPtXiStar[i]->Clone(Form("gMPtXiStarLMRatio_%d", i));
        gMPtLambda1520LMRatio[i] = (TGraphErrors *)gMPtLambda1520[i]->Clone(Form("gMPtLambda1520LMRatio_%d", i));
        gMPtRhoLMRatio[i] = (TGraphErrors *)gMPtRho[i]->Clone(Form("gMPtRhoLMRatio_%d", i));

        gMPtKstarLMRatio[i] = DivideByMult(gMPtKstar[i], 3.69);
        gMPtPionLMRatio[i] = DivideByMult(gMPtPion[i], 3.69);
        gMPtKaonLMRatio[i] = DivideByMult(gMPtKaon[i], 3.69);
        gMPtProtonLMRatio[i] = DivideByMult(gMPtProton[i], 3.69);
        gMPtKshortLMRatio[i] = DivideByMult(gMPtKshort[i], 3.69);
        gMPtChKstarLMRatio[i] = DivideByMult(gMPtChKstar[i], 3.69);
        gMPtXiStarLMRatio[i] = DivideByMult(gMPtXiStar[i], 3.69);
        // gMPtPhiLMRatio[i] = DivideByMult(gMPtPhi[i], 3.69);
        gMPtLambda1520LMRatio[i] = DivideByMult(gMPtLambda1520[i], 3.69);
        gMPtRhoLMRatio[i] = DivideByMult(gMPtRho[i], 3.69);
    }
    gMPtKstarLMRatio[0]->SetMaximum(1.75);
    gMPtKstarLMRatio[0]->SetMinimum(0.88);
    gMPtKstarLMRatio[0]->GetYaxis()->SetTitle("<#it{p}_{T}>/<#it{p}_{T}>_{LM}");
    gMPtKstarLMRatio[0]->Draw("APE");
    gMPtPionLMRatio[0]->Draw("PE same");
    gMPtKaonLMRatio[0]->Draw("PE same");
    gMPtProtonLMRatio[0]->Draw("PE same");
    gMPtKshortLMRatio[0]->Draw("PE same");
    gMPtChKstarLMRatio[0]->Draw("PE same");
    // gMPtXiStarLMRatio[0]->Draw("PE same");
    // gMPtPhiLMRatio[0]->Draw("PE same");
    // gMPtLambda1520LMRatio[0]->Draw("PE same");
    gMPtRhoLMRatio[0]->Draw("PE same");

    gMPtKstarLMRatio[1]->Draw("5 same");
    gMPtPionLMRatio[1]->Draw("5 same");
    gMPtKaonLMRatio[1]->Draw("5 same");
    gMPtProtonLMRatio[1]->Draw("5 same");
    gMPtKshortLMRatio[1]->Draw("5 same");
    gMPtChKstarLMRatio[1]->Draw("5 same");
    // gMPtXiStarLMRatio[1]->Draw("5 same");
    // gMPtPhiLMRatio[1]->Draw("5 same");
    // gMPtLambda1520LMRatio[1]->Draw("5 same");
    gMPtRhoLMRatio[1]->Draw("5 same");

    legendMeanPt->Draw();
    if (isSavePlots)
    {
        cMeanPtRatio->SaveAs("Plots/MeanPt_LowestMultRatio_Run3.png");
    }

    //===================================================
    // =====<pT>_mult / <pT>_lowestMult vs Mass============
    //===================================================

    TCanvas *cMeanPtRatioVsMass = new TCanvas("cMeanPtRatioVsMass", "cMeanPtRatioVsMass", 720, 720);
    SetCanvasStyle(cMeanPtRatioVsMass, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gMeanPtLMRatioVsMass[2];
    for (int i = 0; i < 2; i++)
    {
        gMeanPtLMRatioVsMass[i] = new TGraphErrors();
    }

    for (int i = 0; i < 2; i++)
    {
        gMeanPtLMRatioVsMass[i]->SetPoint(0, PDGMasses[0], gMPtKstarLMRatio[i]->GetY()[FindGraphXPoint(gMPtKstarLMRatio[i], 21.78)]); // K*0
        gMeanPtLMRatioVsMass[i]->SetPointError(0, PDGMassesErrors[0], gMPtKstarLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtKstarLMRatio[i], 21.78)));
        gMeanPtLMRatioVsMass[i]->SetPoint(1, PDGMasses[2], gMPtPionLMRatio[i]->GetY()[FindGraphXPoint(gMPtPionLMRatio[i], 21.78)]); // Pion
        gMeanPtLMRatioVsMass[i]->SetPointError(1, PDGMassesErrors[2], gMPtPionLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtPionLMRatio[i], 21.78)));
        gMeanPtLMRatioVsMass[i]->SetPoint(2, PDGMasses[3], gMPtKaonLMRatio[i]->GetY()[FindGraphXPoint(gMPtKaonLMRatio[i], 21.78)]); // Kaon
        gMeanPtLMRatioVsMass[i]->SetPointError(2, PDGMassesErrors[3], gMPtKaonLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtKaonLMRatio[i], 21.78)));
        gMeanPtLMRatioVsMass[i]->SetPoint(3, PDGMasses[4], gMPtProtonLMRatio[i]->GetY()[FindGraphXPoint(gMPtProtonLMRatio[i], 21.78)]); // Proton
        gMeanPtLMRatioVsMass[i]->SetPointError(3, PDGMassesErrors[4], gMPtProtonLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtProtonLMRatio[i], 21.78)));
        gMeanPtLMRatioVsMass[i]->SetPoint(4, PDGMasses[5], gMPtKshortLMRatio[i]->GetY()[FindGraphXPoint(gMPtKshortLMRatio[i], 21.78)]); // Kshort
        gMeanPtLMRatioVsMass[i]->SetPointError(4, PDGMassesErrors[5], gMPtKshortLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtKshortLMRatio[i], 21.78)));
        gMeanPtLMRatioVsMass[i]->SetPoint(5, PDGMasses[6], gMPtChKstarLMRatio[i]->GetY()[FindGraphXPoint(gMPtChKstarLMRatio[i], 21.78)]); // Charged K*
        gMeanPtLMRatioVsMass[i]->SetPointError(5, PDGMassesErrors[6], gMPtChKstarLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtChKstarLMRatio[i], 21.78)));
        gMeanPtLMRatioVsMass[i]->SetPoint(6, PDGMasses[7], gMPtXiStarLMRatio[i]->GetY()[FindGraphXPoint(gMPtXiStarLMRatio[i], 21.78)]); // Xi(1530)
        gMeanPtLMRatioVsMass[i]->SetPointError(6, PDGMassesErrors[7], gMPtXiStarLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtXiStarLMRatio[i], 21.78)));
        gMeanPtLMRatioVsMass[i]->SetPoint(7, PDGMasses[8], gMPtLambda1520LMRatio[i]->GetY()[FindGraphXPoint(gMPtLambda1520LMRatio[i], 21.78)]); // Lambda(1520)
        gMeanPtLMRatioVsMass[i]->SetPointError(7, PDGMassesErrors[8], gMPtLambda1520LMRatio[i]->GetErrorY(FindGraphXPoint(gMPtLambda1520LMRatio[i], 21.78)));
        gMeanPtLMRatioVsMass[i]->SetPoint(8, PDGMasses[9], gMPtRhoLMRatio[i]->GetY()[FindGraphXPoint(gMPtRhoLMRatio[i], 21.78)]); // Rho(770)
        gMeanPtLMRatioVsMass[i]->SetPointError(8, PDGMassesErrors[9], gMPtRhoLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtRhoLMRatio[i], 21.78)));
        // gMeanPtLMRatioVsMass[i]->SetPoint(8, PDGMasses[1], gMPtPhiLMRatio[i]->GetY()[FindGraphXPoint(gMPtPhiLMRatio[i], 21.78)]); // Phi
        // gMeanPtLMRatioVsMass[i]->SetPointError(8, PDGMassesErrors[1], gMPtPhiLMRatio[i]->GetErrorY(FindGraphXPoint(gMPtPhiLMRatio[i], 21.78)));
    }
    SetGraphErrorStyle(gMeanPtLMRatioVsMass[0]);
    gMeanPtLMRatioVsMass[0]->GetXaxis()->SetTitle("Particle Mass (GeV/#it{c}^{2})");
    gMeanPtLMRatioVsMass[0]->GetYaxis()->SetTitle("<#it{p}_{T}>_{HM} / <#it{p}_{T}>_{LM}");
    gMeanPtLMRatioVsMass[0]->SetMaximum(1.83);
    gMeanPtLMRatioVsMass[0]->SetMinimum(1.24);
    gMeanPtLMRatioVsMass[0]->SetMarkerColor(kBlue);
    gMeanPtLMRatioVsMass[0]->SetLineColor(kBlue);
    gMeanPtLMRatioVsMass[0]->SetMarkerStyle(21);
    gMeanPtLMRatioVsMass[0]->Draw("AP");
    gMeanPtLMRatioVsMass[1]->SetMarkerColor(kBlue);
    gMeanPtLMRatioVsMass[1]->SetLineColor(kBlue);
    gMeanPtLMRatioVsMass[1]->Draw("5 same");
    if (isSavePlots)
    {
        cMeanPtRatioVsMass->SaveAs("Plots/MeanPt_LowestMultRatio_vs_Mass_Run3.png");
    }

    //====================================================
    //   ================Kstar/K Ratio ===================
    //====================================================

    TCanvas *cRatioKstarKaon = new TCanvas("cRatioKstarKaon", "cRatioKstarKaon", 720, 720);
    SetCanvasStyle(cRatioKstarKaon, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarKaon[2];
    gRatioKstarKaon[0] = MakeRatio(gMYieldKstar[0], gMYieldKaon[0], false);
    gRatioKstarKaon[1] = MakeRatio(gMYieldKstar[1], gMYieldKaon[1], false);
    gRatioKstarKaon[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarKaon[0]->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarKaon[0]);
    gRatioKstarKaon[0]->GetYaxis()->SetRangeUser(0.25, 0.43);
    // gRatioKstarKaon[0]->GetYaxis()->SetRangeUser(0.0, 2.66);
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
    double xLM, xHM, yHM, yLM, yHMError, yLMError;
    gRatioKstarKaon[1]->GetPoint(totalPointsKstarKaonRatio - 1, xLM, yLM);
    gRatioKstarKaon[1]->GetPoint(0, xHM, yHM);
    yLMError = gRatioKstarKaon[1]->GetErrorY(totalPointsKstarKaonRatio - 1);
    yHMError = gRatioKstarKaon[1]->GetErrorY(0);
    double sigmaDiffHM_LM = abs(yHM - yLM) / sqrt(pow(yHMError, 2) + pow(yLMError, 2));
    cout << "K*0/Kaon ratio difference between HM and LM: " << sigmaDiffHM_LM << " sigma" << endl;

    // K* is (K* + anit_K*)/2 but K = (K^+ + K^-), so we need to divide the denoimator by 2 as well.
    TGraphErrors *gRatioKstarKa_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gMYieldKaonEPOS_IST0, true, 0.5);
    TGraphErrors *gRatioKstarKa_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gMYieldKaonEPOS_IST0, true, 0.5);

    // gRatioKstarKa_IST9->SetLineStyle(2);
    // gRatioKstarKa_IST9->Draw("l same");
    // gRatioKstarKa_IST9_ITY80->SetLineColor(kBlue - 6);
    // gRatioKstarKa_IST9_ITY80->SetLineStyle(2);
    // gRatioKstarKa_IST9_ITY80->Draw("l same");

    TLegend *legendRatio = new TLegend(0.2, 0.75, 0.5, 0.85);
    SetLegendStyle(legendRatio);
    legendRatio->SetTextSize(0.027);
    legendRatio->AddEntry(gRatioKstarKaon[0], "pp, #sqrt{s} = 13.6 TeV", "p");
    legendRatio->AddEntry(gKstarKaRatio_13TeV[0], "pp, #sqrt{s} = 13 TeV", "p");

    TLegend *legendRatio2 = new TLegend(0.55, 0.72, 0.8, 0.92);
    SetLegendStyle(legendRatio2);
    legendRatio2->SetTextSize(0.027);
    // legendRatio2->AddEntry(gRatioKstarKa_IST9, "EPOS UrQMD OFF", "l");
    // legendRatio2->AddEntry(gRatioKstarKa_IST9_ITY80, "EPOS UrQMD ON", "l");

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
    TGraphErrors *gRatioKstarPion[2];
    gRatioKstarPion[0] = MakeRatio(gMYieldKstar[0], gMYieldPion[0], false);
    gRatioKstarPion[1] = MakeRatio(gMYieldKstar[1], gMYieldPion[1], false);
    gRatioKstarPion[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarPion[0]->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarPion[0]);
    gRatioKstarPion[0]->GetYaxis()->SetRangeUser(0.026, 0.059);
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
    gKstarPiRatio_7TeV[0]->Draw("P same");
    gKstarPiRatio_7TeV[1]->SetLineColor(kGreen + 2);
    gKstarPiRatio_7TeV[1]->SetFillStyle(0);
    gKstarPiRatio_7TeV[1]->Draw("5 same");

    // K* is (K* + anit_K*)/2, but for pion it is (Pi^+ + Pi^-), so no need a factor 2 correction.
    TGraphErrors *gRatioKstarPi_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gMYieldPionEPOS_IST0, true, 2.0);
    TGraphErrors *gRatioKstarPi_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gMYieldPionEPOS_IST0, true, 2.0);

    // gRatioKstarPi_IST9->SetLineStyle(2);
    // gRatioKstarPi_IST9->Draw("l same");
    // gRatioKstarPi_IST9_ITY80->SetLineColor(kBlue);
    // gRatioKstarPi_IST9_ITY80->SetLineStyle(2);
    // gRatioKstarPi_IST9_ITY80->Draw("l same");

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

    TGraphErrors *gRatioChKstarKaon[2], *gRatioXiStarKaon[2], *gRatioChKstarPion[2], *gRatioXiStarPiPi[2], *gRatioLambda1520Kaon[2], *gRatioLambda1520Pion[2], *gRatioLambda1520Proton[2], *gRatioPhiKaon[2];

    TGraphErrors *gYieldKstarKaLMRatio[2], *gYieldChKstarKaLMRatio[2], *gYieldXiStarKaLMRatio[2], *gYieldKstarPiLMRatio[2], *gYieldChKstarPiLMRatio[2], *gYieldXiStarPiPiLMRatio[2], *gYieldLambda1520KaLMRatio[2], *gYieldLambda1520PiLMRatio[2], *gYieldLambda1520ProtonLMRatio[2], *gYieldPhiKaLMRatio[2];
    for (int i = 0; i < 2; i++)
    {
        gYieldKstarKaLMRatio[i] = (TGraphErrors *)gRatioKstarKaon[i]->Clone(Form("gYieldKstarKaLMRatio_%d", i));
        gYieldKstarKaLMRatio[i] = DivideByMult(gRatioKstarKaon[i], 3.69);
        gYieldKstarPiLMRatio[i] = (TGraphErrors *)gRatioKstarPion[i]->Clone(Form("gYieldKstarPiLMRatio_%d", i));
        gYieldKstarPiLMRatio[i] = DivideByMult(gRatioKstarPion[i], 3.69);

        int totalPoints = gYieldKstarKaLMRatio[i]->GetN();
        for (int j = 0; j < totalPoints - 1; j++)
        {
            double x, yKa, yPi;
            gYieldKstarKaLMRatio[i]->GetPoint(j, x, yKa);
            gYieldKstarPiLMRatio[i]->GetPoint(j, x, yPi);
            double percentageDifference = abs(yKa - yPi) / yKa * 100.0;
            double DifferenceSigma = abs(yKa - yPi) / sqrt(pow(gYieldKstarKaLMRatio[i]->GetErrorY(j), 2) + pow(gYieldKstarPiLMRatio[i]->GetErrorY(j), 2));
            cout << "Multiplicity: " << x << ", K*/K: " << yKa << ", K*/Pi: " << yPi << ", Percentage Difference: " << std::round(percentageDifference * 10.0) / 10.0 << "%, Difference in Sigma: " << std::round(DifferenceSigma * 10.0) / 10.0 << endl;
        }
        cout << endl;

        gRatioChKstarKaon[i] = MakeRatio(gMYieldChKstar[i], gMYieldKaon[i], false);
        gYieldChKstarKaLMRatio[i] = (TGraphErrors *)gRatioChKstarKaon[i]->Clone(Form("gYieldChKstarKaLMRatio_%d", i));
        gYieldChKstarKaLMRatio[i] = DivideByMult(gRatioChKstarKaon[i], 3.69);
        gRatioChKstarPion[i] = MakeRatio(gMYieldChKstar[i], gMYieldPion[i], false);
        gYieldChKstarPiLMRatio[i] = (TGraphErrors *)gRatioChKstarPion[i]->Clone(Form("gYieldChKstarPiLMRatio_%d", i));
        gYieldChKstarPiLMRatio[i] = DivideByMult(gRatioChKstarPion[i], 3.69);

        gRatioXiStarKaon[i] = MakeRatio(gMYieldXiStar[i], gMYieldKaon[i], false);
        // gRatioXiStarKaon[i] = MakeRatio(gRatioXiStarKaon[i], gMYieldKaon[i], false); // Divide by Kaon again to get XiStar/Kaon^2
        gYieldXiStarKaLMRatio[i] = (TGraphErrors *)gRatioXiStarKaon[i]->Clone(Form("gYieldXiStarKaLMRatio_%d", i));
        gYieldXiStarKaLMRatio[i] = DivideByMult(gRatioXiStarKaon[i], 3.69);
        gRatioXiStarPiPi[i] = MakeRatio(gMYieldXiStar[i], gMYieldPion[i], false);
        // gRatioXiStarPiPi[i] = MakeRatio(gRatioXiStarPiPi[i], gMYieldPion[i], false); // Divide by pion again
        gYieldXiStarPiPiLMRatio[i] = (TGraphErrors *)gRatioXiStarPiPi[i]->Clone(Form("gYieldXiStarPiPiLMRatio_%d", i));
        gYieldXiStarPiPiLMRatio[i] = DivideByMult(gRatioXiStarPiPi[i], 3.69);

        gRatioLambda1520Kaon[i] = MakeRatio(gMYieldLambda1520[i], gMYieldKaon[i], false);
        gYieldLambda1520KaLMRatio[i] = (TGraphErrors *)gRatioLambda1520Kaon[i]->Clone(Form("gYieldLambda1520KaLMRatio_%d", i));
        gYieldLambda1520KaLMRatio[i] = DivideByMult(gRatioLambda1520Kaon[i], 3.69);
        gRatioLambda1520Pion[i] = MakeRatio(gMYieldLambda1520[i], gMYieldPion[i], false);
        gYieldLambda1520PiLMRatio[i] = (TGraphErrors *)gRatioLambda1520Pion[i]->Clone(Form("gYieldLambda1520PiLMRatio_%d", i));
        gYieldLambda1520PiLMRatio[i] = DivideByMult(gRatioLambda1520Pion[i], 3.69);
        gRatioLambda1520Proton[i] = MakeRatio(gMYieldLambda1520[i], gMYieldProton[i], false);
        gYieldLambda1520ProtonLMRatio[i] = (TGraphErrors *)gRatioLambda1520Proton[i]->Clone(Form("gYieldLambda1520ProtonLMRatio_%d", i));
        gYieldLambda1520ProtonLMRatio[i] = DivideByMult(gRatioLambda1520Proton[i], 3.69);

        gRatioPhiKaon[i] = MakeRatio(gMYieldPhi[i], gMYieldKaon[i], false);
        gYieldPhiKaLMRatio[i] = (TGraphErrors *)gRatioPhiKaon[i]->Clone(Form("gYieldPhiKaLMRatio_%d", i));
        gYieldPhiKaLMRatio[i] = DivideByMult(gRatioPhiKaon[i], 3.69);
    }

    gYieldKstarKaLMRatio[0]->SetMaximum(1.39);
    gYieldKstarKaLMRatio[0]->SetMinimum(0.65);
    gYieldKstarKaLMRatio[0]->SetMarkerColor(kRed);
    gYieldKstarKaLMRatio[0]->SetLineColor(kRed);
    gYieldKstarKaLMRatio[0]->SetMarkerStyle(20);
    gYieldKstarKaLMRatio[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gYieldKstarKaLMRatio[0]->GetYaxis()->SetTitle("Y/Y_{LM}");
    gYieldKstarKaLMRatio[0]->Draw("APE");
    gYieldKstarKaLMRatio[1]->Draw("5 same");
    // gYieldChKstarKaLMRatio[0]->SetMarkerColor(kBlue);
    // gYieldChKstarKaLMRatio[0]->SetLineColor(kBlue);
    // gYieldChKstarKaLMRatio[0]->SetMarkerStyle(22);
    // gYieldChKstarKaLMRatio[0]->Draw("PE same");
    // gYieldChKstarKaLMRatio[1]->SetLineColor(kBlue);
    // gYieldChKstarKaLMRatio[1]->SetFillStyle(0);
    // gYieldChKstarKaLMRatio[1]->Draw("5 same");

    // gYieldXiStarKaLMRatio[0]->SetMarkerColor(kOrange - 2);
    // gYieldXiStarKaLMRatio[0]->SetLineColor(kOrange - 2);
    // gYieldXiStarKaLMRatio[0]->SetMarkerStyle(22);
    // // gYieldXiStarKaLMRatio[0]->Draw("PE same");
    // gYieldXiStarKaLMRatio[1]->SetLineColor(kOrange - 2);
    // gYieldXiStarKaLMRatio[1]->SetFillStyle(0);
    // // gYieldXiStarKaLMRatio[1]->Draw("5 same");

    gYieldKstarPiLMRatio[0]->SetMarkerColor(kGreen + 2);
    gYieldKstarPiLMRatio[0]->SetLineColor(kGreen + 2);
    gYieldKstarPiLMRatio[0]->SetMarkerStyle(21);
    gYieldKstarPiLMRatio[0]->Draw("PE same");
    gYieldKstarPiLMRatio[1]->SetLineColor(kGreen + 2);
    gYieldKstarPiLMRatio[1]->SetFillStyle(0);
    gYieldKstarPiLMRatio[1]->Draw("5 same");
    // gYieldChKstarPiLMRatio[0]->SetMarkerColor(kMagenta);
    // gYieldChKstarPiLMRatio[0]->SetLineColor(kMagenta);
    // gYieldChKstarPiLMRatio[0]->SetMarkerStyle(47);
    // gYieldChKstarPiLMRatio[0]->Draw("PE same");
    // gYieldChKstarPiLMRatio[1]->SetLineColor(kMagenta);
    // gYieldChKstarPiLMRatio[1]->SetFillStyle(0);
    // gYieldChKstarPiLMRatio[1]->Draw("5 same");
    // gYieldXiStarPiPiLMRatio[0]->SetMarkerColor(kOrange - 2);
    // gYieldXiStarPiPiLMRatio[0]->SetLineColor(kOrange - 2);
    // gYieldXiStarPiPiLMRatio[0]->SetMarkerStyle(26);
    // gYieldXiStarPiPiLMRatio[0]->Draw("PE same");
    // gYieldXiStarPiPiLMRatio[1]->SetLineColor(kOrange - 2);
    // gYieldXiStarPiPiLMRatio[1]->SetFillStyle(0);
    // gYieldXiStarPiPiLMRatio[1]->Draw("5 same");

    TLegend *legendYieldLMRatio = new TLegend(0.2, 0.86, 0.8, 0.95);
    SetLegendStyle(legendYieldLMRatio);
    legendYieldLMRatio->SetTextSize(0.03);
    legendYieldLMRatio->SetNColumns(2);
    legendYieldLMRatio->AddEntry(gYieldKstarKaLMRatio[0], "K*^{0}/K", "P");
    legendYieldLMRatio->AddEntry(gYieldKstarPiLMRatio[0], "K*^{0}/#pi", "P");
    // legendYieldLMRatio->AddEntry(gYieldXiStarKaLMRatio[0], "#Xi(1530)/K", "P");
    // legendYieldLMRatio->AddEntry(gYieldChKstarKaLMRatio[0], "K*^{#pm}/K", "P");
    // legendYieldLMRatio->AddEntry(gYieldChKstarPiLMRatio[0], "K*^{#pm}/#pi", "P");
    // legendYieldLMRatio->AddEntry(gYieldXiStarPiPiLMRatio[0], "#Xi(1530)/#pi", "P");

    TLegend *legendYieldLMRatio2 = new TLegend(0.2, 0.8, 0.8, 0.89);
    SetLegendStyle(legendYieldLMRatio2);
    legendYieldLMRatio2->SetTextSize(0.03);
    legendYieldLMRatio2->SetNColumns(2);

    // for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    for (int imodel = 1; imodel < 2; imodel++)
    {
        TGraphErrors *gRatioKstarKaPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kKaon_epos], true, 0.5);
        setStyle(gRatioKstarKaPythiaModel, colorsPythia[imodel], lineStylesPythia[imodel]);
        // gRatioKstarKaPythiaModel->Draw("l same");
        // legendRatio2->AddEntry(gRatioKstarKaPythiaModel, modelLabelLocal[imodel], "l");
        TGraphErrors *gYieldKstarKaLMRatioPythia = (TGraphErrors *)gRatioKstarKaPythiaModel->Clone(Form("gYieldKstarKaLMRatioPythia_%d", imodel));
        gYieldKstarKaLMRatioPythia = DivideByMult(gRatioKstarKaPythiaModel, 3.69);
        gYieldKstarKaLMRatioPythia->SetLineWidth(3);
        gYieldKstarKaLMRatioPythia->SetLineColor(kMagenta + 1);
        gYieldKstarKaLMRatioPythia->SetLineStyle(2);
        gYieldKstarKaLMRatioPythia->Draw("l same");

        TGraphErrors *gRatioKstarPiPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kPion_epos], true, 0.5);
        TGraphErrors *gYieldKstarPiLMRatioPythia = (TGraphErrors *)gRatioKstarPiPythiaModel->Clone(Form("gYieldKstarPiLMRatioPythia_%d", imodel));
        gYieldKstarPiLMRatioPythia = DivideByMult(gRatioKstarPiPythiaModel, 3.69);
        gYieldKstarPiLMRatioPythia->SetLineWidth(3);
        gYieldKstarPiLMRatioPythia->SetLineColor(kCyan + 1);
        gYieldKstarPiLMRatioPythia->Draw("l same");
        legendYieldLMRatio2->AddEntry(gYieldKstarKaLMRatioPythia, Form("%s ", modelLabelLocal[imodel]), "l");
        legendYieldLMRatio2->AddEntry(gYieldKstarPiLMRatioPythia, Form("%s ", modelLabelLocal[imodel]), "l");
    }
    legendYieldLMRatio->Draw();
    legendYieldLMRatio2->Draw();

    TLine *lineDoubleRatio = new TLine(gYieldKstarKaLMRatio[0]->GetXaxis()->GetXmin(), 1.0, gYieldKstarKaLMRatio[0]->GetXaxis()->GetXmax(), 1.0);
    lineDoubleRatio->SetLineStyle(2);
    lineDoubleRatio->SetLineWidth(2);
    // lineDoubleRatio->Draw("same");
    if (isSavePlots)
    {
        cYieldLMRatio->SaveAs("Plots/Yield_LMRatio_KstarKaPi_Run3.png");
    }

    //======================================================================
    // ==Double YieldRatio (Mult/LM): Kstar, ChKstar, Lambda1520, XiStar, Phi, Rho==
    //======================================================================
    TCanvas *cYieldLMRatio2 = new TCanvas("cYieldLMRatio2", "cYieldLMRatio2", 720, 720);
    SetCanvasStyle(cYieldLMRatio2, 0.15, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cYieldLMRatio2, pad1Size, pad2Size);
    cYieldLMRatio2->cd(1);

    TGraphErrors *gYieldKstarLMRatio[2], *gYieldChKstarLMRatio[2], *gYieldXiStarLMRatio[2], *gYieldLambda1520LMRatio[2], *gYieldPhiLMRatio[2], *gYieldRhoLMRatio[2];
    TGraphErrors *gRatio2KstarKstar[2], *gRatio2ChKstarKstar[2], *gRatio2XiStarKstar[2], *gRatio2Lambda1520Kstar[2], *gRatio2PhiKstar[2], *gRatio2RhoKstar[2];
    for (int i = 0; i < 2; i++)
    {
        gYieldKstarLMRatio[i] = (TGraphErrors *)gMYieldKstar[i]->Clone(Form("gYieldKstarLMRatio_%d", i));
        gYieldKstarLMRatio[i] = DivideByMult(gMYieldKstar[i], 3.69);
        gYieldChKstarLMRatio[i] = (TGraphErrors *)gMYieldChKstar[i]->Clone(Form("gYieldChKstarLMRatio_%d", i));
        gYieldChKstarLMRatio[i] = DivideByMult(gMYieldChKstar[i], 3.69);
        gYieldXiStarLMRatio[i] = (TGraphErrors *)gMYieldXiStar[i]->Clone(Form("gYieldXiStarLMRatio_%d", i));
        gYieldXiStarLMRatio[i] = DivideByMult(gMYieldXiStar[i], 3.69);
        gYieldLambda1520LMRatio[i] = (TGraphErrors *)gMYieldLambda1520[i]->Clone(Form("gYieldLambda1520LMRatio_%d", i));
        gYieldLambda1520LMRatio[i] = DivideByMult(gMYieldLambda1520[i], 3.69);
        gYieldPhiLMRatio[i] = (TGraphErrors *)gMYieldPhi[i]->Clone(Form("gYieldPhiLMRatio_%d", i));
        gYieldPhiLMRatio[i] = DivideByMult(gMYieldPhi[i], 2.5);
        gYieldRhoLMRatio[i] = (TGraphErrors *)gMYieldRho[i]->Clone(Form("gYieldRhoLMRatio_%d", i));
        gYieldRhoLMRatio[i] = DivideByMult(gMYieldRho[i], 3.69);

        gRatio2KstarKstar[i] = MakeRatio(gYieldKstarLMRatio[i], gYieldKstarLMRatio[i], false);
        gRatio2ChKstarKstar[i] = MakeRatio(gYieldChKstarLMRatio[i], gYieldKstarLMRatio[i], false);
        gRatio2XiStarKstar[i] = MakeRatio(gYieldXiStarLMRatio[i], gYieldKstarLMRatio[i], false);
        gRatio2Lambda1520Kstar[i] = MakeRatio(gYieldLambda1520LMRatio[i], gYieldKstarLMRatio[i], false);
        gRatio2PhiKstar[i] = MakeRatio(gYieldPhiLMRatio[i], gYieldKstarLMRatio[i], false);
        gRatio2RhoKstar[i] = MakeRatio(gYieldRhoLMRatio[i], gYieldKstarLMRatio[i], false);
    }
    gYieldKstarLMRatio[0]->SetMaximum(16.5);
    gYieldKstarLMRatio[0]->SetMinimum(0.0);
    gYieldKstarLMRatio[0]->SetMarkerColor(kRed);
    gYieldKstarLMRatio[0]->SetLineColor(kRed);
    gYieldKstarLMRatio[0]->SetMarkerStyle(20);
    gYieldKstarLMRatio[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gYieldKstarLMRatio[0]->GetYaxis()->SetTitle("Y/Y_{LM}");
    gYieldKstarLMRatio[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    gYieldKstarLMRatio[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    gYieldKstarLMRatio[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    gYieldKstarLMRatio[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    gYieldKstarLMRatio[0]->GetYaxis()->SetTitleOffset(1.7 * pad1Size);
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
    gYieldXiStarLMRatio[0]->SetMarkerColor(kOrange - 2);
    gYieldXiStarLMRatio[0]->SetLineColor(kOrange - 2);
    gYieldXiStarLMRatio[0]->SetMarkerStyle(22);
    gYieldXiStarLMRatio[0]->Draw("PE same");
    gYieldXiStarLMRatio[1]->SetLineColor(kOrange - 2);
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
    gYieldRhoLMRatio[0]->SetMarkerColor(kCyan + 1);
    gYieldRhoLMRatio[0]->SetLineColor(kCyan + 1);
    gYieldRhoLMRatio[0]->SetMarkerStyle(25);
    gYieldRhoLMRatio[0]->Draw("PE same");
    gYieldRhoLMRatio[1]->SetLineColor(kCyan + 1);
    gYieldRhoLMRatio[1]->SetFillStyle(0);
    gYieldRhoLMRatio[1]->Draw("5 same");

    TLegend *legendYieldLMRatio3 = new TLegend(0.2, 0.78, 0.8, 0.9);
    SetLegendStyle(legendYieldLMRatio3);
    legendYieldLMRatio3->SetTextSize(0.027 / pad1Size);
    legendYieldLMRatio3->SetNColumns(3);
    legendYieldLMRatio3->AddEntry(gYieldKstarLMRatio[0], "K*^{0}", "P");
    legendYieldLMRatio3->AddEntry(gYieldChKstarLMRatio[0], "K*^{#pm}", "P");
    legendYieldLMRatio3->AddEntry(gYieldXiStarLMRatio[0], "#Xi(1530)", "P");
    legendYieldLMRatio3->AddEntry(gYieldLambda1520LMRatio[0], "#Lambda(1520)", "P");
    legendYieldLMRatio3->AddEntry(gYieldPhiLMRatio[0], "#phi", "P");
    legendYieldLMRatio3->AddEntry(gYieldRhoLMRatio[0], "#rho", "P");
    legendYieldLMRatio3->Draw();

    cYieldLMRatio2->cd(2);

    gRatio2ChKstarKstar[0]->SetMaximum(2.58);
    gRatio2ChKstarKstar[0]->SetMinimum(0.65);
    gRatio2ChKstarKstar[0]->GetYaxis()->SetNdivisions(505);
    gRatio2ChKstarKstar[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatio2ChKstarKstar[0]->GetYaxis()->SetTitle("Ratio to K*^{0}");
    gRatio2ChKstarKstar[0]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    gRatio2ChKstarKstar[0]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    gRatio2ChKstarKstar[0]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    gRatio2ChKstarKstar[0]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    gRatio2ChKstarKstar[0]->GetYaxis()->SetTitleOffset(1.7 * pad2Size);
    gRatio2ChKstarKstar[0]->GetXaxis()->SetTitleOffset(1.15);
    gRatio2ChKstarKstar[0]->SetMarkerColor(kBlue);
    gRatio2ChKstarKstar[0]->SetLineColor(kBlue);
    gRatio2ChKstarKstar[0]->SetMarkerStyle(21);
    gRatio2ChKstarKstar[0]->Draw("APE");
    gRatio2ChKstarKstar[1]->SetLineColor(kBlue);
    gRatio2ChKstarKstar[1]->SetFillStyle(0);
    gRatio2ChKstarKstar[1]->Draw("5 same");
    gRatio2XiStarKstar[0]->SetMarkerColor(kOrange - 2);
    gRatio2XiStarKstar[0]->SetLineColor(kOrange - 2);
    gRatio2XiStarKstar[0]->SetMarkerStyle(22);
    gRatio2XiStarKstar[0]->Draw("PE same");
    gRatio2XiStarKstar[1]->SetLineColor(kOrange - 2);
    gRatio2XiStarKstar[1]->SetFillStyle(0);
    gRatio2XiStarKstar[1]->Draw("5 same");
    gRatio2Lambda1520Kstar[0]->SetMarkerColor(kGreen + 2);
    gRatio2Lambda1520Kstar[0]->SetLineColor(kGreen + 2);
    gRatio2Lambda1520Kstar[0]->SetMarkerStyle(22);
    gRatio2Lambda1520Kstar[0]->Draw("PE same");
    gRatio2Lambda1520Kstar[1]->SetLineColor(kGreen + 2);
    gRatio2Lambda1520Kstar[1]->SetFillStyle(0);
    gRatio2Lambda1520Kstar[1]->Draw("5 same");
    gRatio2PhiKstar[0]->SetMarkerColor(kMagenta);
    gRatio2PhiKstar[0]->SetLineColor(kMagenta);
    gRatio2PhiKstar[0]->SetMarkerStyle(47);
    gRatio2PhiKstar[0]->Draw("PE same");
    gRatio2PhiKstar[1]->SetLineColor(kMagenta);
    gRatio2PhiKstar[1]->SetFillStyle(0);
    gRatio2PhiKstar[1]->Draw("5 same");
    gRatio2RhoKstar[0]->SetMarkerColor(kCyan + 1);
    gRatio2RhoKstar[0]->SetLineColor(kCyan + 1);
    gRatio2RhoKstar[0]->SetMarkerStyle(25);
    gRatio2RhoKstar[0]->Draw("PE same");
    gRatio2RhoKstar[1]->SetLineColor(kCyan + 1);
    gRatio2RhoKstar[1]->SetFillStyle(0);
    gRatio2RhoKstar[1]->Draw("5 same");

    TLine *line = new TLine(gRatio2ChKstarKstar[0]->GetXaxis()->GetXmin(), 1.0, gRatio2ChKstarKstar[0]->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("same");
    if (isSavePlots)
    {
        cYieldLMRatio2->SaveAs("Plots/Yield_LMRatio2.png");
    }

    //======================================================================
    // ==Yield Ratio (HM/LM) vs lifetime: Kstar, ChKstar, Lambda1520, XiStar, Phi, Rho==
    //======================================================================
    TCanvas *cYieldLifetime = new TCanvas("cYieldLifetime", "cYieldLifetime", 720, 720);
    SetCanvasStyle(cYieldLifetime, 0.15, 0.03, 0.03, 0.15);
    float lifetime[6] = {3.9, 4.2, 12.6, 22, 46.2, 1.3}; // K*+-, K*0, Lambda1520, XiStar, Phi, Rho
    TGraphErrors *gYieldVsLifetime[2];
    for (int i = 0; i < 2; i++)
    {
        gYieldVsLifetime[i] = new TGraphErrors(6);
    }

    for (int j = 0; j < 2; j++)
    {
        gYieldVsLifetime[j]->SetPoint(0, lifetime[0], gYieldChKstarLMRatio[j]->GetY()[FindGraphXPoint(gYieldChKstarLMRatio[j], 21.78)]); // K*+-
        gYieldVsLifetime[j]->SetPointError(0, 0.3, gYieldChKstarLMRatio[j]->GetErrorY(FindGraphXPoint(gYieldChKstarLMRatio[j], 21.78)));
        gYieldVsLifetime[j]->SetPoint(1, lifetime[1], gYieldKstarLMRatio[j]->GetY()[FindGraphXPoint(gYieldKstarLMRatio[j], 21.78)]); // K*0
        gYieldVsLifetime[j]->SetPointError(1, 0.3, gYieldKstarLMRatio[j]->GetErrorY(FindGraphXPoint(gYieldKstarLMRatio[j], 21.78)));
        gYieldVsLifetime[j]->SetPoint(2, lifetime[2], gYieldLambda1520LMRatio[j]->GetY()[FindGraphXPoint(gYieldLambda1520LMRatio[j], 21.78)]); // Lambda1520
        gYieldVsLifetime[j]->SetPointError(2, 0.3, gYieldLambda1520LMRatio[j]->GetErrorY(FindGraphXPoint(gYieldLambda1520LMRatio[j], 21.78)));
        gYieldVsLifetime[j]->SetPoint(3, lifetime[3], gYieldXiStarLMRatio[j]->GetY()[FindGraphXPoint(gYieldXiStarLMRatio[j], 21.78)]); // XiStar
        gYieldVsLifetime[j]->SetPointError(3, 0.3, gYieldXiStarLMRatio[j]->GetErrorY(FindGraphXPoint(gYieldXiStarLMRatio[j], 21.78)));
        gYieldVsLifetime[j]->SetPoint(4, lifetime[4], gYieldPhiLMRatio[j]->GetY()[FindGraphXPoint(gYieldPhiLMRatio[j], 25.78)]); // Phi
        gYieldVsLifetime[j]->SetPointError(4, 0.3, gYieldPhiLMRatio[j]->GetErrorY(FindGraphXPoint(gYieldPhiLMRatio[j], 25.78)));
        gYieldVsLifetime[j]->SetPoint(5, lifetime[5], gYieldRhoLMRatio[j]->GetY()[FindGraphXPoint(gYieldRhoLMRatio[j], 21.78)]); // Rho
        gYieldVsLifetime[j]->SetPointError(5, 0.3, gYieldRhoLMRatio[j]->GetErrorY(FindGraphXPoint(gYieldRhoLMRatio[j], 21.78)));
    }
    SetGraphErrorStyle(gYieldVsLifetime[0]);
    gYieldVsLifetime[0]->SetMinimum(4.2);
    gYieldVsLifetime[0]->SetMaximum(17.5);
    gYieldVsLifetime[0]->SetMarkerColor(kRed);
    gYieldVsLifetime[0]->SetLineColor(kRed);
    gYieldVsLifetime[0]->SetMarkerStyle(20);
    gYieldVsLifetime[0]->GetXaxis()->SetTitle("Lifetime (fm/c)");
    gYieldVsLifetime[0]->GetYaxis()->SetTitle("Y_{HM}/Y_{LM}");
    gYieldVsLifetime[0]->Draw("APE");
    gYieldVsLifetime[1]->SetLineColor(kRed);
    gYieldVsLifetime[1]->SetFillStyle(0);
    gYieldVsLifetime[1]->Draw("5 same");

    TLatex latex2;
    latex2.SetTextSize(0.03);
    latex2.SetTextAlign(23); // center horizontally, top vertically

    const char *names[6] = {
        "K^{*#pm}",
        "K^{*0}",
        "#Lambda(1520)",
        "#Xi^{*}",
        "#phi",
        "#rho"};

    for (int i = 0; i < 6; i++)
    {
        double x, y;
        gYieldVsLifetime[0]->GetPoint(i, x, y);

        // place label slightly below point
        if (i == 0)
            latex2.DrawLatex(x - 0.5, y - 0.7, names[i]);
        if (i == 1)
            latex2.DrawLatex(x + 1.9, y - 0.6, names[i]);
        if (i == 2)
            latex2.DrawLatex(x, y - 0.9, names[i]);
        if (i == 3)
            latex2.DrawLatex(x, y - 3.8, names[i]);
        if (i == 4)
            latex2.DrawLatex(x, y - 2.0, names[i]);
        if (i == 5)
            latex2.DrawLatex(x, y - 0.7, names[i]);
    }
    if (isSavePlots)
    {
        cYieldLifetime->SaveAs("Plots/Yield_Lifetime_Ratio.png");
    }

    // /*
    //====================================================
    //   ================Kstar/ Kshort Ratio ===================
    //====================================================
    TCanvas *cRatioKstarKshort = new TCanvas("cRatioKstarKshort", "cRatioKstarKshort", 720, 720);
    SetCanvasStyle(cRatioKstarKshort, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarKshort = MakeRatio(gMYieldKstar[0], gMYieldKshort[0], false);
    TGraphErrors *gRatioKstarKshort_sys = MakeRatio(gMYieldKstar[1], gMYieldKshort[1], false);
    gRatioKstarKshort->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarKshort->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarKshort);
    gRatioKstarKshort->GetYaxis()->SetRangeUser(0.18, 0.46);
    gRatioKstarKshort->GetXaxis()->SetLimits(0, 27);
    gRatioKstarKshort->SetMarkerColor(kRed);
    gRatioKstarKshort->SetLineColor(kRed);
    gRatioKstarKshort->Draw("APE");
    gRatioKstarKshort_sys->SetFillStyle(0);
    gRatioKstarKshort_sys->SetLineColor(kRed);
    gRatioKstarKshort_sys->Draw("5 same");
    gKstarKshortRatio_13TeV[0]->SetMarkerStyle(21);
    gKstarKshortRatio_13TeV[0]->SetMarkerColor(kBlue);
    gKstarKshortRatio_13TeV[0]->SetLineColor(kBlue);
    gKstarKshortRatio_13TeV[0]->Draw("P same");
    gKstarKshortRatio_13TeV[1]->SetLineColor(kBlue);
    gKstarKshortRatio_13TeV[1]->SetFillStyle(0);
    gKstarKshortRatio_13TeV[1]->Draw("5 same");

    // K* is (K* + anit_K*)/2 and Kshort does not have antiparticle, so no need for correction.
    TGraphErrors *gRatioKstarKshort_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gMYieldKshortEPOS_IST0, true);
    TGraphErrors *gRatioKstarKshort_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gMYieldKshortEPOS_IST0, true);

    gRatioKstarKshort_IST9->SetLineStyle(2);
    gRatioKstarKshort_IST9->Draw("l same");
    gRatioKstarKshort_IST9_ITY80->SetLineColor(kBlue);
    gRatioKstarKshort_IST9_ITY80->SetLineStyle(2);
    gRatioKstarKshort_IST9_ITY80->Draw("l same");

    // Here K*0 is not averaged.
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarKshortModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKshort], true, 2);
        setStyle(gRatioKstarKshortModel, modelStyle[model].color, modelStyle[model].style);
        gRatioKstarKshortModel->Draw("l same");
        // legendRatio->AddEntry(gRatioKstarKshortModel, modelLabel[model], "l");
    }

    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        TGraphErrors *gRatioKstarKshortPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kKshort_epos], true);
        setStyle(gRatioKstarKshortPythiaModel, modelStyle[imodel].color, modelStyle[imodel].style);
        gRatioKstarKshortPythiaModel->Draw("l same");
    }

    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0} + #bar{K}^{*0}}{2K_{S}^{0}}");
    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{K_{S}^{0}}");
    legendRatio->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
        cRatioKstarKshort->SaveAs("Plots/Ratio_KstarKshort_Run3.png");
    }

    //====================================================
    // ====================Kstar/Phi Ratio ===================
    //====================================================

    TCanvas *cRatioKstarPhi = new TCanvas("cRatioKstarPhi", "cRatioKstarPhi", 720, 720);
    SetCanvasStyle(cRatioKstarPhi, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarPhi = MakeRatio(gMYieldKstar[0], gMYieldPhi[0], false);
    TGraphErrors *gRatioKstarPhi_sys = MakeRatio(gMYieldKstar[1], gMYieldPhi[1], false);
    gRatioKstarPhi->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarPhi->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarPhi);
    gRatioKstarPhi->GetYaxis()->SetRangeUser(1.2, 6.7);
    gRatioKstarPhi->GetXaxis()->SetLimits(0, 27);
    gRatioKstarPhi->Draw("APE");
    gRatioKstarPhi_sys->SetFillStyle(0);
    gRatioKstarPhi_sys->Draw("5 same");

    ////Here Kstar is average
    TGraphErrors *gRatioKstarPhi_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gEPOS_Yield[9][kITY0][kPhi_epos], true);
    TGraphErrors *gRatioKstarPhi_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gEPOS_Yield[9][kITY80][kPhi_epos], true);

    gRatioKstarPhi_IST9->SetLineStyle(2);
    gRatioKstarPhi_IST9->Draw("l same");
    gRatioKstarPhi_IST9_ITY80->SetLineColor(kBlue);
    gRatioKstarPhi_IST9_ITY80->SetLineStyle(2);
    gRatioKstarPhi_IST9_ITY80->Draw("l same");

    // Here K*0 is not averaged.
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarPhiModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kPhi], true, 2.0);
        setStyle(gRatioKstarPhiModel, modelStyle[model].color, modelStyle[model].style);
        gRatioKstarPhiModel->Draw("l same");
    }
    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{#phi}");
    legendRatio->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
        cRatioKstarPhi->SaveAs("Plots/Ratio_KstarPhi_Run3.png");
    }

    //====================================================
    //  ==================Kstar / Charged Kstar Ratio===================
    //====================================================

    TCanvas *cRatioKstarChargedKstar = new TCanvas("cRatioKstarChargedKstar", "cRatioKstarChargedKstar", 720, 720);
    SetCanvasStyle(cRatioKstarChargedKstar, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarChargedKstar = MakeRatio(gMYieldKstar[0], gMYieldChKstar[0], false);
    TGraphErrors *gRatioKstarChargedKstar_sys = MakeRatio(gMYieldKstar[1], gMYieldChKstar[1], false);
    gRatioKstarChargedKstar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarChargedKstar->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarChargedKstar);
    gRatioKstarChargedKstar->GetYaxis()->SetRangeUser(0.5, 1.4);
    gRatioKstarChargedKstar->GetXaxis()->SetLimits(0, 27);
    gRatioKstarChargedKstar->Draw("APE");
    gRatioKstarChargedKstar_sys->SetFillStyle(0);
    gRatioKstarChargedKstar_sys->Draw("5 same");

    // factor 2 is due to average in K*0 and another factor 2 is due wrong BR taken 0.333 instead of 0.666.
    TGraphErrors *gRatioKstarChargedKstar_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gEPOS_Yield[9][kITY0][kKstarPM_epos], true, 0.25);
    TGraphErrors *gRatioKstarChargedKstar_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gEPOS_Yield[9][kITY80][kKstarPM_epos], true, 0.25);

    gRatioKstarChargedKstar_IST9->SetLineStyle(2);
    gRatioKstarChargedKstar_IST9->Draw("l same");
    gRatioKstarChargedKstar_IST9_ITY80->SetLineColor(kBlue);
    gRatioKstarChargedKstar_IST9_ITY80->SetLineStyle(2);
    gRatioKstarChargedKstar_IST9_ITY80->Draw("l same");

    // Some issue in K*+- simulation. It has very less yield than K*0 somehow (a facotr of 700 difference.)
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarChargedKstarModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKstarPM], true);
        setStyle(gRatioKstarChargedKstarModel, modelStyle[model].color, modelStyle[model].style);
        gRatioKstarChargedKstarModel->Draw("l same");
    }
    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        TGraphErrors *gRatioKstarChargedKstarPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kKstarPM_epos], true, 0.5);
        setStyle(gRatioKstarChargedKstarPythiaModel, modelStyle[imodel].color, modelStyle[imodel].style);
        gRatioKstarChargedKstarPythiaModel->Draw("l same");
    }
    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{K^{* #pm}}");
    legendRatio->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
        cRatioKstarChargedKstar->SaveAs("Plots/Ratio_KstarChargedKstar_Run3.png");
    }

    //================================================
    // ==================Kstar/ Xi* Ratio===================
    //====================================================
    TCanvas *cRatioKstarXiStar = new TCanvas("cRatioKstarXiStar", "cRatioKstarXiStar", 720, 720);
    SetCanvasStyle(cRatioKstarXiStar, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarXiStar = MakeRatio(gMYieldKstar[0], gMYieldXiStar[0], false);
    TGraphErrors *gRatioKstarXiStar_sys = MakeRatio(gMYieldKstar[1], gMYieldXiStar[1], false);
    gRatioKstarXiStar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarXiStar->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarXiStar);
    gRatioKstarXiStar->GetYaxis()->SetRangeUser(15, 83);
    gRatioKstarXiStar->GetXaxis()->SetLimits(0, 27);
    gRatioKstarXiStar->Draw("APE");
    gRatioKstarXiStar_sys->SetFillStyle(0);
    gRatioKstarXiStar_sys->Draw("5 same");

    TGraphErrors *gRatioKstarXiStar_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gEPOS_Yield[9][kITY0][kXi1530_epos], true, 0.5);
    TGraphErrors *gRatioKstarXiStar_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gEPOS_Yield[9][kITY80][kXi1530_epos], true, 0.5);
    TGraphErrors *gRatioKstarXiStar_IST9_ITY81 = MakeRatio(gEPOS_Yield[9][kITY81][kKstar_epos], gEPOS_Yield[9][kITY81][kXi1530_epos], true, 0.5);
    TGraphErrors *gRatioKstarXiStar_IST6 = MakeRatio(gEPOS_Yield[6][kITY0][kKstar_epos], gEPOS_Yield[6][kITY0][kXi1530_epos], true, 0.5);

    gRatioKstarXiStar_IST9->SetLineStyle(2);
    gRatioKstarXiStar_IST9->Draw("l same");
    gRatioKstarXiStar_IST9_ITY80->SetLineColor(kBlue);
    gRatioKstarXiStar_IST9_ITY80->SetLineStyle(2);
    gRatioKstarXiStar_IST9_ITY80->Draw("l same");
    gRatioKstarXiStar_IST9_ITY81->SetLineColor(kGreen + 2);
    gRatioKstarXiStar_IST9_ITY81->SetLineStyle(2);

    // In this model, the Kstar is not divided by 2
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarXiStarModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kXi1530], true, 2);
        setStyle(gRatioKstarXiStarModel, modelStyle[model].color, modelStyle[model].style);
        gRatioKstarXiStarModel->Draw("l same");
    }
    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        TGraphErrors *gRatioKstarXiStarPythiaModel = MakeRatio(gPythiaYieldLocal[imodel][kKstar_epos], gPythiaYieldLocal[imodel][kXi1530_epos], true);
        setStyle(gRatioKstarXiStarPythiaModel, modelStyle[imodel].color, modelStyle[imodel].style);
        gRatioKstarXiStarPythiaModel->Draw("l same");
    }
    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{#Xi(1530)^{0}}");
    legendRatio->Clear();
    legendRatio->AddEntry(gRatioKstarXiStar, "pp, #sqrt{s} = 13.6 TeV", "p");
    legendRatio->Draw();
    legendRatio2->Draw();
    if (isSavePlots)
    {
        cRatioKstarXiStar->SaveAs("Plots/Ratio_KstarXiStar_Run3.png");
    }

    /*
    // //===================================================
    // // ======================Phi/ Pi Ratio===================
    // //====================================================

    // // ******** There is mistake in either 7 TeV or 13 TeV, because one of them is exactly off by a factor of 2 *********
    // TCanvas *cRatioPhiPion = new TCanvas("cRatioPhiPion", "cRatioPhiPion", 720, 720);
    // SetCanvasStyle(cRatioPhiPion, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioPhiPion = MakeRatio(gMYieldPhi[0], gMYieldPion[0], false);
    // TGraphErrors *gRatioPhiPion_sys = MakeRatio(gMYieldPhi[1], gMYieldPion[1], false);
    // gRatioPhiPion->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioPhiPion->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioPhiPion);
    // gRatioPhiPion->GetYaxis()->SetRangeUser(0.002, 0.032);
    // gRatioPhiPion->GetXaxis()->SetLimits(0, 27);
    // gRatioPhiPion->SetMarkerColor(kRed);
    // gRatioPhiPion->SetLineColor(kRed);
    // gRatioPhiPion->Draw("APE");
    // gRatioPhiPion_sys->SetFillStyle(0);
    // gRatioPhiPion_sys->SetLineColor(kRed);
    // gRatioPhiPion_sys->Draw("5 same");
    // gPhiPiRatio_13TeV[0]->SetMarkerStyle(21);
    // gPhiPiRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gPhiPiRatio_13TeV[0]->SetLineColor(kBlue);
    // gPhiPiRatio_13TeV[0]->Draw("P same");
    // gPhiPiRatio_13TeV[1]->SetLineColor(kBlue);
    // gPhiPiRatio_13TeV[1]->SetFillStyle(0);
    // gPhiPiRatio_13TeV[1]->Draw("5 same");
    // // ScaleGraph(gPhiPiRatio_7TeV[0], 0.5); // In 7 TeV it was just K*0 + K*0bar and not its average
    // // ScaleGraph(gPhiPiRatio_7TeV[1], 0.5);
    // gPhiPiRatio_7TeV[0]->SetMarkerStyle(22);
    // gPhiPiRatio_7TeV[0]->SetMarkerColor(kGreen + 2);
    // gPhiPiRatio_7TeV[0]->SetLineColor(kGreen + 2);
    // gPhiPiRatio_7TeV[0]->Draw("P same");
    // gPhiPiRatio_7TeV[1]->SetLineColor(kGreen + 2);
    // gPhiPiRatio_7TeV[1]->SetFillStyle(0);
    // gPhiPiRatio_7TeV[1]->Draw("5 same");

    // TGraphErrors *gRatioPhiPion_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gMYieldPionEPOS_IST0, true);
    // TGraphErrors *gRatioPhiPion_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gMYieldPionEPOS_IST0, true);

    // gRatioPhiPion_IST9->SetLineStyle(2);
    // gRatioPhiPion_IST9->Draw("l same");
    // gRatioPhiPion_IST9_ITY80->SetLineColor(kBlue);
    // gRatioPhiPion_IST9_ITY80->SetLineStyle(2);
    // gRatioPhiPion_IST9_ITY80->Draw("l same");

    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioPhiPionModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kPion], true);
    //     setStyle(gRatioPhiPionModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioPhiPionModel->Draw("l same");
    // }
    // latex.DrawLatex(0.28, 0.88, "#frac{#phi}{#pi}");
    // legendRatio3->Draw();
    // legendRatio2->Draw();
    // if (isSavePlots)
    // {
    //   cRatioPhiPion->SaveAs("Plots/Ratio_PhiPion_Run3.png");
    // }

    // //====================================================
    // //   ================Phi/K Ratio===================
    // //====================================================
    // TCanvas *cRatioPhiKaon = new TCanvas("cRatioPhiKaon", "cRatioPhiKaon", 720, 720);
    // SetCanvasStyle(cRatioPhiKaon, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioPhiKaon = MakeRatio(gMYieldPhi[0], gMYieldKaon[0], false);
    // TGraphErrors *gRatioPhiKaon_sys = MakeRatio(gMYieldPhi[1], gMYieldKaon[1], false);
    // gRatioPhiKaon->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioPhiKaon->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioPhiKaon);
    // gRatioPhiKaon->GetYaxis()->SetRangeUser(0.0, 0.32);
    // gRatioPhiKaon->GetXaxis()->SetLimits(0, 27);
    // gRatioPhiKaon->SetMarkerColor(kRed);
    // gRatioPhiKaon->SetLineColor(kRed);
    // gRatioPhiKaon->Draw("APE");
    // gRatioPhiKaon_sys->SetFillStyle(0);
    // gRatioPhiKaon_sys->SetLineColor(kRed);
    // gRatioPhiKaon_sys->Draw("5 same");
    // gPhiKaRatio_13TeV[0]->SetMarkerStyle(21);
    // gPhiKaRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gPhiKaRatio_13TeV[0]->SetLineColor(kBlue);
    // gPhiKaRatio_13TeV[0]->Draw("P same");
    // gPhiKaRatio_13TeV[1]->SetLineColor(kBlue);
    // gPhiKaRatio_13TeV[1]->SetFillStyle(0);
    // gPhiKaRatio_13TeV[1]->Draw("5 same");

    // ////K = (K^+ + K^-) /2, so no need for correction.
    // TGraphErrors *gRatioPhiKaon_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gMYieldKaonEPOS_IST0, true);
    // TGraphErrors *gRatioPhiKaon_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gMYieldKaonEPOS_IST0, true);

    // gRatioPhiKaon_IST9->SetLineStyle(2);
    // gRatioPhiKaon_IST9->Draw("l same");
    // gRatioPhiKaon_IST9_ITY80->SetLineColor(kBlue);
    // gRatioPhiKaon_IST9_ITY80->SetLineStyle(2);
    // gRatioPhiKaon_IST9_ITY80->Draw("l same");

    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioPhiKaonModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kKaon], true);
    //     setStyle(gRatioPhiKaonModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioPhiKaonModel->Draw("l same");
    // }

    // latex.DrawLatex(0.28, 0.88, "#frac{#phi}{K}");
    // legendRatio->Draw();
    // legendRatio2->Draw();
    //if (isSavePlots)
    //{
    // cRatioPhiKaon->SaveAs("Plots/Ratio_PhiKaon_Run3.png");
    //}

    // //===================================================
    // //   ================Phi/Kshort Ratio===================
    // //====================================================

    // TCanvas *cRatioPhiKshort = new TCanvas("cRatioPhiKshort", "cRatioPhiKshort", 720, 720);
    // SetCanvasStyle(cRatioPhiKshort, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioPhiKshort = MakeRatio(gMYieldPhi[0], gMYieldKshort[0], false);
    // TGraphErrors *gRatioPhiKshort_sys = MakeRatio(gMYieldPhi[1], gMYieldKshort[1], false);
    // gRatioPhiKshort->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioPhiKshort->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioPhiKshort);
    // gRatioPhiKshort->GetYaxis()->SetRangeUser(0.0, 0.32);
    // gRatioPhiKshort->GetXaxis()->SetLimits(0, 27);
    // gRatioPhiKshort->SetMarkerColor(kRed);
    // gRatioPhiKshort->SetLineColor(kRed);
    // gRatioPhiKshort->Draw("APE");
    // gRatioPhiKshort_sys->SetFillStyle(0);
    // gRatioPhiKshort_sys->SetLineColor(kRed);
    // gRatioPhiKshort_sys->Draw("5 same");
    // gPhiKshortRatio_13TeV[0]->SetMarkerStyle(21);
    // gPhiKshortRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gPhiKshortRatio_13TeV[0]->SetLineColor(kBlue);
    // gPhiKshortRatio_13TeV[0]->Draw("P same");
    // gPhiKshortRatio_13TeV[1]->SetLineColor(kBlue);
    // gPhiKshortRatio_13TeV[1]->SetFillStyle(0);
    // gPhiKshortRatio_13TeV[1]->Draw("5 same");

    // TGraphErrors *gRatioPhiKshort_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gMYieldKshortEPOS_IST0, true, 0.5);
    // TGraphErrors *gRatioPhiKshort_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gMYieldKshortEPOS_IST0, true, 0.5);

    // gRatioPhiKshort_IST9->SetLineStyle(2);
    // gRatioPhiKshort_IST9->Draw("l same");
    // gRatioPhiKshort_IST9_ITY80->SetLineColor(kBlue);
    // gRatioPhiKshort_IST9_ITY80->SetLineStyle(2);
    // gRatioPhiKshort_IST9_ITY80->Draw("l same");

    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioPhiKshortModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kKshort], true);
    //     setStyle(gRatioPhiKshortModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioPhiKshortModel->Draw("l same");
    // }
    // latex.DrawLatex(0.28, 0.88, "#frac{#phi}{K_{S}^{0}}");
    // legendRatio->Draw();
    // legendRatio2->Draw();
    //if (isSavePlots)
    //{
    // cRatioPhiKshort->SaveAs("Plots/Ratio_PhiKshort_Run3.png");
    //}

    // //===================================================
    // //  ==================Phi/ Proton Ratio===================
    // //====================================================
    // TCanvas *cRatioPhiProton = new TCanvas("cRatioPhiProton", "cRatioPhiProton", 720, 720);
    // SetCanvasStyle(cRatioPhiProton, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioPhiProton = MakeRatio(gMYieldPhi[0], gMYieldProton[0], false, 2);
    // TGraphErrors *gRatioPhiProton_sys = MakeRatio(gMYieldPhi[1], gMYieldProton[1], false, 2);
    // gRatioPhiProton->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioPhiProton->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioPhiProton);
    // gRatioPhiProton->GetYaxis()->SetRangeUser(0.0, 0.259);
    // gRatioPhiProton->GetXaxis()->SetLimits(0, 27);
    // gRatioPhiProton->SetMarkerColor(kRed);
    // gRatioPhiProton->SetLineColor(kRed);
    // gRatioPhiProton->Draw("APE");
    // gRatioPhiProton_sys->SetFillStyle(0);
    // gRatioPhiProton_sys->SetLineColor(kRed);
    // gRatioPhiProton_sys->Draw("5 same");
    // gPhiPrRatio_13TeV[0]->SetMarkerStyle(21);
    // gPhiPrRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gPhiPrRatio_13TeV[0]->SetLineColor(kBlue);
    // gPhiPrRatio_13TeV[0]->Draw("P same");
    // gPhiPrRatio_13TeV[1]->SetLineColor(kBlue);
    // gPhiPrRatio_13TeV[1]->SetFillStyle(0);
    // gPhiPrRatio_13TeV[1]->Draw("5 same");

    // TGraphErrors *gRatioPhiProton_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gMYieldProtonEPOS_IST0, true);
    // TGraphErrors *gRatioPhiProton_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gMYieldProtonEPOS_IST0, true);

    // gRatioPhiProton_IST9->SetLineStyle(2);
    // gRatioPhiProton_IST9->Draw("l same");
    // gRatioPhiProton_IST9_ITY80->SetLineColor(kBlue);
    // gRatioPhiProton_IST9_ITY80->SetLineStyle(2);
    // gRatioPhiProton_IST9_ITY80->Draw("l same");

    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioPhiProtonModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kProton], true);
    //     setStyle(gRatioPhiProtonModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioPhiProtonModel->Draw("l same");
    // }
    // latex.DrawLatex(0.28, 0.88, "#frac{#phi}{p}");
    // legendRatio->Draw();
    // legendRatio2->Draw();
    // if (isSavePlots)
    // {
    //     cRatioPhiProton->SaveAs("Plots/Ratio_PhiProton_Run3.png");
    // }

    // //================================================
    // // ==================Phi/ Xi* Ratio===================
    // //====================================================

    // TCanvas *cRatioPhiXiStar = new TCanvas("cRatioPhiXiStar", "cRatioPhiXiStar", 720, 720);
    // SetCanvasStyle(cRatioPhiXiStar, 0.15, 0.03, 0.03, 0.15);
    // TGraphErrors *gRatioPhiXiStar = MakeRatio(gMYieldPhi[0], gMYieldXiStar[0], false);
    // TGraphErrors *gRatioPhiXiStar_sys = MakeRatio(gMYieldPhi[1], gMYieldXiStar[1], false);
    // gRatioPhiXiStar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    // gRatioPhiXiStar->GetYaxis()->SetTitle("dN/dy");
    // SetGraphErrorStyle(gRatioPhiXiStar);
    // gRatioPhiXiStar->GetYaxis()->SetRangeUser(0.0, 209);
    // gRatioPhiXiStar->GetXaxis()->SetLimits(0, 27);
    // gRatioPhiXiStar->Draw("APE");
    // gRatioPhiXiStar_sys->SetFillStyle(0);
    // gRatioPhiXiStar_sys->Draw("5 same");

    // TGraphErrors *gRatioPhiXiStar_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kPhi_epos], gEPOS_Yield[9][kITY0][kXi1530_epos], true);
    // TGraphErrors *gRatioPhiXiStar_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kPhi_epos], gEPOS_Yield[9][kITY80][kXi1530_epos], true);
    // TGraphErrors *gRatioPhiXiStar_IST9_ITY81 = MakeRatio(gEPOS_Yield[9][kITY81][kPhi_epos], gEPOS_Yield[9][kITY81][kXi1530_epos], true);
    // TGraphErrors *gRatioPhiXiStar_IST6 = MakeRatio(gEPOS_Yield[6][kITY0][kPhi_epos], gEPOS_Yield[6][kITY0][kXi1530_epos], true);

    // gRatioPhiXiStar_IST9->SetLineStyle(2);
    // gRatioPhiXiStar_IST9->Draw("l same");
    // gRatioPhiXiStar_IST9_ITY80->SetLineColor(kBlue);
    // gRatioPhiXiStar_IST9_ITY80->SetLineStyle(2);
    // gRatioPhiXiStar_IST9_ITY80->Draw("l same");
    // gRatioPhiXiStar_IST9_ITY81->SetLineColor(kGreen + 2);
    // gRatioPhiXiStar_IST9_ITY81->SetLineStyle(2);
    // gRatioPhiXiStar_IST9_ITY81->Draw("l same");
    // gRatioPhiXiStar_IST6->SetLineColor(kMagenta);
    // gRatioPhiXiStar_IST6->SetLineStyle(5);
    // // gRatioPhiXiStar_IST6->Draw("l same");

    // for (auto model : modelsToPlot)
    // {
    //     TGraphErrors *gRatioPhiXiStarModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kXi1530], true);
    //     setStyle(gRatioPhiXiStarModel, modelStyle[model].color, modelStyle[model].style);
    //     gRatioPhiXiStarModel->Draw("l same");
    // }
    // latex.DrawLatex(0.28, 0.88, "#frac{#phi}{2#Xi^{*0}}");
    // legendRatio->Draw();
    // legendRatio2->Draw();
    //if (isSavePlots)
    //{
    // cRatioPhiXiStar->SaveAs("Plots/Ratio_PhiXiStar_Run3.png");
    // }

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
    //if (isSavePlots)
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
    //if (isSavePlots)
   // {
    //     cRatioProtonPion->SaveAs("Plots/Ratio_ProtonPion_Run3.png");
    //}
    */

    //====================================================
    // ==================Pion yeild======================
    //====================================================
    TCanvas *cPionYield = new TCanvas("cPionYield", "cPionYield", 720, 720);
    SetCanvasStyle(cPionYield, 0.15, 0.03, 0.03, 0.15);
    gMYieldPion[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldPion[0]->GetYaxis()->SetTitle("dN/dy");
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
    // ScaleGraph(gMYieldPionEPOS_IST0, 2.0);
    gMYieldPionEPOS_IST0->Draw("l same");

    // for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    // {
    //     setStyle(gPythiaYieldLocal[imodel][kPion_epos], modelStyle[imodel].color, modelStyle[imodel].style);
    //     ScaleGraph(gPythiaYieldLocal[imodel][kPion_epos], 0.5); // In pythia nothing is averaged, but in data Pi,K,p is averaged.
    //     gPythiaYieldLocal[imodel][kPion_epos]->Draw("l same");
    // }
    TLegend *legTemp = new TLegend(0.2, 0.65, 0.55, 0.85);
    SetLegendStyle(legTemp);
    legTemp->SetTextSize(0.03);
    legTemp->AddEntry(gMYieldPion[0], "Data", "p");
    legTemp->AddEntry(gMYieldPionEPOS_IST0, "EPOS", "l");
    // for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    // {
    //     legTemp->AddEntry(gPythiaYieldLocal[imodel][kPion_epos], modelLabelLocal[imodel], "l");
    // }
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
    gMYieldKaon[0]->GetYaxis()->SetTitle("dN/dy");
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
    // ScaleGraph(gMYieldKaonEPOS_IST0, 0.5);
    gMYieldKaonEPOS_IST0->Draw("l same");

    // for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    // {
    //     setStyle(gPythiaYieldLocal[imodel][kKaon_epos], modelStyle[imodel].color, modelStyle[imodel].style);
    //     ScaleGraph(gPythiaYieldLocal[imodel][kKaon_epos], 0.5); // In pythia nothing is averaged, but in data Pi,K,p is averaged.
    //     gPythiaYieldLocal[imodel][kKaon_epos]->Draw("l same");
    // }
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
    gMYieldProton[0]->GetYaxis()->SetTitle("dN/dy");
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
    // ScaleGraph(gMYieldProtonEPOS_IST0, 0.5);
    gMYieldProtonEPOS_IST0->Draw("l same");

    // for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    // {
    //     setStyle(gPythiaYieldLocal[imodel][kProton_epos], modelStyle[imodel].color, modelStyle[imodel].style);
    //     ScaleGraph(gPythiaYieldLocal[imodel][kProton_epos], 0.5); // In pythia nothing is averaged, but in data Pi,K,p is averaged.
    //     gPythiaYieldLocal[imodel][kProton_epos]->Draw("l same");
    // }
    legTemp->Draw();
    latex.DrawLatex(0.28, 0.88, "(p + #bar{p})/2");
    // if (isSavePlots)
    {
        cProtonYield->SaveAs("Plots/ProtonYield_Run3.png");
    }

    //====================================================
    // ==================Phi yeild======================
    //====================================================
    TCanvas *cPhiYield = new TCanvas("cPhiYield", "cPhiYield", 720, 720);
    SetCanvasStyle(cPhiYield, 0.15, 0.03, 0.03, 0.15);
    gMYieldPhi[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldPhi[0]->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gMYieldPhi[0]);
    gMYieldPhi[0]->GetYaxis()->SetRangeUser(0.0, 0.25);
    gMYieldPhi[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldPhi[0]->SetMarkerColor(kRed);
    gMYieldPhi[0]->SetLineColor(kRed);
    gMYieldPhi[0]->Draw("APE");
    gMYieldPhi[1]->SetFillStyle(0);
    gMYieldPhi[1]->SetLineColor(kRed);
    gMYieldPhi[1]->Draw("5 same");

    gEPOS_Yield[9][kITY0][kPhi_epos]->SetLineStyle(9);
    gEPOS_Yield[9][kITY0][kPhi_epos]->Draw("l same");
    gEPOS_Yield[9][kITY80][kPhi_epos]->SetLineColor(kBlue - 6);
    gEPOS_Yield[9][kITY80][kPhi_epos]->SetLineStyle(1);
    gEPOS_Yield[9][kITY80][kPhi_epos]->Draw("l same");

    TLegend *legTempPhi = new TLegend(0.2, 0.65, 0.55, 0.85);
    SetLegendStyle(legTempPhi);
    legTempPhi->SetTextSize(0.03);
    legTempPhi->AddEntry(gMYieldPhi[0], "Data", "p");
    legTempPhi->AddEntry(gEPOS_Yield[9][kITY0][kPhi_epos], "EPOS UrQMD OFF", "l");
    legTempPhi->AddEntry(gEPOS_Yield[9][kITY80][kPhi_epos], "EPOS UrQMD ON", "l");
    legTempPhi->Draw();
    latex.DrawLatex(0.28, 0.88, "#phi");
    // if (isSavePlots)
    {
        cPhiYield->SaveAs("Plots/PhiYield_Run3.png");
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

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.35);
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

TGraphErrors *DivideByMult(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5)
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


*/