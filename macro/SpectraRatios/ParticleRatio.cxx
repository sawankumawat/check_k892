#include <iostream>
#include <iomanip>
#include "../src/style.h"

TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false, float CorrFactorDen = 1);
TGraphErrors *GetGraph(TFile *f, const string &name);

void setStyle(TGraphErrors *gr, int color, int style)
{
    gr->SetLineColor(color);
    gr->SetLineStyle(style);
}

void ScaleGraph(TGraphErrors *gr, double scale)
{
    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        gr->SetPoint(i, x, y * scale);
        gr->SetPointError(i, gr->GetErrorX(i), gr->GetErrorY(i) * scale);
    }
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

void ParticleRatio()
{
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

    if (fPion->IsZombie() || fProton->IsZombie() || fKaon->IsZombie() || fPhi->IsZombie() || fChKstar->IsZombie() || fXiStar->IsZombie() || fKshort->IsZombie())
    {
        cout << "Error: Particle files not found" << endl;
        return;
    }

    // FIXME: Add 3rd column for uncorrelated systematic error.
    TGraphErrors *gMPtKstar[2], *gMYieldKstar[2], *gMPtPhi[2], *gMPtPion[2], *gMPtProton[2], *gMPtKaon[2], *gMPtChKstar[2], *gMPtXiStar[2], *gMPtKshort[2], *gMYieldPhi[2], *gMYieldPion[2], *gMYieldProton[2], *gMYieldKaon[2], *gMYieldChKstar[2], *gMYieldXiStar[2], *gMYieldKshort[2];

    for (int i = 0; i < 2; i++)
    {
        gMPtKstar[i] = GetGraph(fKstar, Form("gMeanpTRun3%s", (i == 0 ? "" : "_sys")));
        gMYieldKstar[i] = GetGraph(fKstar, Form("gMeanYieldRun3%s", (i == 0 ? "" : "_sys")));

        gMPtPhi[i] = GetGraph(fPhi, Form("gPhi_MeanpT%s", (i == 0 ? "_stat" : "_sys")));
        gMPtPion[i] = GetGraph(fPion, Form("gMeanpTRun3%s", (i == 0 ? "" : "_sys")));
        gMPtProton[i] = GetGraph(fProton, Form("gMeanpTRun3%s", (i == 0 ? "" : "_sys")));
        gMPtKaon[i] = GetGraph(fKaon, Form("gMeanpTRun3%s", (i == 0 ? "" : "_sys")));
        gMPtChKstar[i] = GetGraph(fChKstar, Form("gMeanpTRun3%s", (i == 0 ? "" : "_sys")));
        gMPtXiStar[i] = GetGraph(fXiStar, Form("gMeanpTRun3%s", (i == 0 ? "" : "_sys")));
        gMPtKshort[i] = GetGraph(fKshort, Form("gMeanpTRun3%s", (i == 0 ? "" : "_sys")));

        gMYieldPhi[i] = GetGraph(fPhi, Form("gPhi_MeanYield%s", (i == 0 ? "_stat" : "_sys")));
        gMYieldPion[i] = GetGraph(fPion, Form("gMeanYieldRun3%s", (i == 0 ? "" : "_sys")));
        gMYieldProton[i] = GetGraph(fProton, Form("gMeanYieldRun3%s", (i == 0 ? "" : "_sys")));
        gMYieldKaon[i] = GetGraph(fKaon, Form("gMeanYieldRun3%s", (i == 0 ? "" : "_sys")));
        gMYieldChKstar[i] = GetGraph(fChKstar, Form("gMeanYieldRun3%s", (i == 0 ? "" : "_sys")));
        gMYieldXiStar[i] = GetGraph(fXiStar, Form("gMeanYieldRun3%s", (i == 0 ? "" : "_sys")));
        gMYieldKshort[i] = GetGraph(fKshort, Form("gMeanYieldRun3%s", (i == 0 ? "" : "_sys")));
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
    // TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA_ptCut.root", "read");
    TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA_ptCut_FinerBins.root", "read");
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

    vector<int> modelsToPlot = {
        kPythiaMonash,
        kPythiaShoving,
        kPythiaRopes};
    // vector<int> modelsToPlot = {};

    TCanvas *cdNdyKstar = new TCanvas("cdNdyKstar", "cdNdyKstar", 720, 720);
    SetCanvasStyle(cdNdyKstar, 0.15, 0.03, 0.03, 0.15);
    gMYieldKstar[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldKstar[0]->GetYaxis()->SetTitle("dN/dy");
    gMYieldKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldKstar[0]->GetYaxis()->SetRangeUser(0, 0.8);
    gMYieldKstar[0]->SetMarkerColor(kRed);
    gMYieldKstar[0]->SetLineColor(kRed);
    gMYieldKstar[0]->Draw("APE");
    gMYieldKstar[1]->SetFillStyle(0);
    gMYieldKstar[1]->SetLineColor(kRed);
    gMYieldKstar[1]->Draw("5 same");
    gKstar_Yield_13TeV[0]->SetMarkerStyle(21);
    gKstar_Yield_13TeV[0]->SetMarkerColor(kBlue);
    gKstar_Yield_13TeV[0]->SetLineColor(kBlue);
    gKstar_Yield_13TeV[0]->Draw("P same");
    gKstar_Yield_13TeV[1]->SetLineColor(kBlue);
    gKstar_Yield_13TeV[1]->SetFillStyle(0);
    gKstar_Yield_13TeV[1]->Draw("5 same");

    gEPOS_Yield[9][kITY0][kKstar_epos]->SetLineStyle(6);
    gEPOS_Yield[9][kITY0][kKstar_epos]->Draw("l same");
    gEPOS_Yield[9][kITY80][kKstar_epos]->SetLineColor(kBlue);
    gEPOS_Yield[9][kITY80][kKstar_epos]->SetLineStyle(2);
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
    legend->Draw();
    legend2->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.027);
    // latex.DrawLatex(0.7, 0.85, "#frac{K* (892)^{0} + #bar{K}* (892)^{0}}{2}");
    latex.DrawLatex(0.28, 0.9, "K* (892)^{0}");
    cdNdyKstar->SaveAs("Plots/MeanYield_Kstar_EPOS_UrQMDON.png");

    // /*
    TCanvas *cMeanPtKstar = new TCanvas("cMeanPtKstar", "cMeanPtKstar", 720, 720);
    SetCanvasStyle(cMeanPtKstar, 0.15, 0.03, 0.03, 0.15);
    SetGraphErrorStyle(gMPtKstar[0]);
    gMPtKstar[0]->SetTitle(0);
    gMPtKstar[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMPtKstar[0]->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    gMPtKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMPtKstar[0]->GetYaxis()->SetRangeUser(0.32, 1.96);
    gMPtKstar[0]->SetMarkerColor(kRed);
    gMPtKstar[0]->SetLineColor(kRed);
    gMPtKstar[0]->Draw("APE");
    gMPtKstar[1]->SetFillStyle(0);
    gMPtKstar[1]->SetLineColor(kRed);
    gMPtKstar[1]->Draw("5 same");
    gKstar_MeanpT_13TeV[0]->SetMarkerStyle(21);
    gKstar_MeanpT_13TeV[0]->SetMarkerColor(kBlue);
    gKstar_MeanpT_13TeV[0]->SetLineColor(kBlue);
    gKstar_MeanpT_13TeV[0]->Draw("P same");
    gKstar_MeanpT_13TeV[1]->SetLineColor(kBlue);
    gKstar_MeanpT_13TeV[1]->SetFillStyle(0);
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
    legend->Draw();
    legend2->Draw();
    // latex.DrawLatex(0.7, 0.85, "#frac{K* (892)^{0} + #bar{K}* (892)^{0}}{2}");
    latex.DrawLatex(0.28, 0.9, "K* (892)^{0}");
    cMeanPtKstar->SaveAs("Plots/MeanPt_Kstar_EPOS_UrQMDON.png");

    //====================================================
    //   ================Kstar/K Ratio ===================
    //====================================================

    TCanvas *cRatioKstarKaon = new TCanvas("cRatioKstarKaon", "cRatioKstarKaon", 720, 720);
    SetCanvasStyle(cRatioKstarKaon, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarKaon = MakeRatio(gMYieldKstar[0], gMYieldKaon[0], false);
    TGraphErrors *gRatioKstarKaon_sys = MakeRatio(gMYieldKstar[1], gMYieldKaon[1], false);
    gRatioKstarKaon->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarKaon->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarKaon);
    gRatioKstarKaon->GetYaxis()->SetRangeUser(0.18, 0.86);
    gRatioKstarKaon->GetXaxis()->SetLimits(0, 27);
    gRatioKstarKaon->SetMarkerColor(kRed);
    gRatioKstarKaon->SetLineColor(kRed);
    gRatioKstarKaon->Draw("APE");
    gRatioKstarKaon_sys->SetFillStyle(0);
    gRatioKstarKaon_sys->SetLineColor(kRed);
    gRatioKstarKaon_sys->Draw("5 same");
    gKstarKaRatio_13TeV[0]->SetMarkerStyle(21);
    gKstarKaRatio_13TeV[0]->SetMarkerColor(kBlue);
    gKstarKaRatio_13TeV[0]->SetLineColor(kBlue);
    gKstarKaRatio_13TeV[0]->Draw("P same");
    gKstarKaRatio_13TeV[1]->SetLineColor(kBlue);
    gKstarKaRatio_13TeV[1]->SetFillStyle(0);
    gKstarKaRatio_13TeV[1]->Draw("5 same");

    // K* is (K* + anit_K*)/2 but K = (K^+ + K^-), so we need correction by factor 2.
    TGraphErrors *gRatioKstarKa_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gMYieldKaonEPOS_IST0, true, 0.5);
    TGraphErrors *gRatioKstarKa_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gMYieldKaonEPOS_IST0, true, 0.5);

    gRatioKstarKa_IST9->SetLineStyle(2);
    gRatioKstarKa_IST9->Draw("l same");
    gRatioKstarKa_IST9_ITY80->SetLineColor(kBlue);
    gRatioKstarKa_IST9_ITY80->SetLineStyle(2);
    gRatioKstarKa_IST9_ITY80->Draw("l same");

    TLegend *legendRatio = new TLegend(0.2, 0.75, 0.5, 0.85);
    SetLegendStyle(legendRatio);
    legendRatio->SetTextSize(0.027);
    legendRatio->AddEntry(gRatioKstarKaon, "pp, #sqrt{s} = 13.6 TeV", "p");
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

    legendRatio->Draw();
    legendRatio2->Draw();
    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0} + #bar{K}^{*0}}{K^{+} + K^{-}}");
    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{K}");
    cRatioKstarKaon->SaveAs("Plots/Ratio_KstarKaon_Run3.png");

    //==================================================
    // ===================Kstar / Pion Ratio ==================
    //==================================================
    TCanvas *cRatioKstarPion = new TCanvas("cRatioKstarPion", "cRatioKstarPion", 720, 720);
    SetCanvasStyle(cRatioKstarPion, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarPion = MakeRatio(gMYieldKstar[0], gMYieldPion[0], false);
    TGraphErrors *gRatioKstarPion_sys = MakeRatio(gMYieldKstar[1], gMYieldPion[1], false);
    gRatioKstarPion->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarPion->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarPion);
    gRatioKstarPion->GetYaxis()->SetRangeUser(0.016, 0.069);
    gRatioKstarPion->GetXaxis()->SetLimits(0, 27);
    gRatioKstarPion->SetMarkerColor(kRed);
    gRatioKstarPion->SetLineColor(kRed);
    gRatioKstarPion->Draw("APE");
    gRatioKstarPion_sys->SetFillStyle(0);
    gRatioKstarPion_sys->SetLineColor(kRed);
    gRatioKstarPion_sys->Draw("5 same");
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

    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{#pi}");
    TLegend *legendRatio3 = (TLegend *)legendRatio->Clone("legendRatio3");
    legendRatio3->AddEntry(gKstarPiRatio_7TeV[0], "pp, #sqrt{s} = 7 TeV", "P");
    legendRatio3->Draw();
    legendRatio2->Draw();
    cRatioKstarPion->SaveAs("Plots/Ratio_KstarPion_Run3.png");

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

    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0} + #bar{K}^{*0}}{2K_{S}^{0}}");
    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{K_{S}^{0}}");
    legendRatio->Draw();
    legendRatio2->Draw();
    cRatioKstarKshort->SaveAs("Plots/Ratio_KstarKshort_Run3.png");

    // //===================================================
    // // ======================Phi/ Pi Ratio===================
    // //====================================================

    // //******** There is mistake in either 7 TeV or 13 TeV, because one of them is exactly off by a factor of 2 *********
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
    // cRatioPhiPion->SaveAs("Plots/Ratio_PhiPion_Run3.png");

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
    // cRatioPhiKaon->SaveAs("Plots/Ratio_PhiKaon_Run3.png");

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
    cRatioPhiKshort->SaveAs("Plots/Ratio_PhiKshort_Run3.png");

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
    // cRatioPhiProton->SaveAs("Plots/Ratio_PhiProton_Run3.png");

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
    cRatioKstarPhi->SaveAs("Plots/Ratio_KstarPhi_Run3.png");

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
    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{K^{* #pm}}");
    legendRatio->Draw();
    legendRatio2->Draw();
    cRatioKstarChargedKstar->SaveAs("Plots/Ratio_KstarChargedKstar_Run3.png");

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
    gRatioKstarXiStar->GetYaxis()->SetRangeUser(0.0, 179);
    gRatioKstarXiStar->GetXaxis()->SetLimits(0, 27);
    gRatioKstarXiStar->Draw("APE");
    gRatioKstarXiStar_sys->SetFillStyle(0);
    gRatioKstarXiStar_sys->Draw("5 same");

    TGraphErrors *gRatioKstarXiStar_IST9 = MakeRatio(gEPOS_Yield[9][kITY0][kKstar_epos], gEPOS_Yield[9][kITY0][kXi1530_epos], true);
    TGraphErrors *gRatioKstarXiStar_IST9_ITY80 = MakeRatio(gEPOS_Yield[9][kITY80][kKstar_epos], gEPOS_Yield[9][kITY80][kXi1530_epos], true);
    TGraphErrors *gRatioKstarXiStar_IST9_ITY81 = MakeRatio(gEPOS_Yield[9][kITY81][kKstar_epos], gEPOS_Yield[9][kITY81][kXi1530_epos], true);
    TGraphErrors *gRatioKstarXiStar_IST6 = MakeRatio(gEPOS_Yield[6][kITY0][kKstar_epos], gEPOS_Yield[6][kITY0][kXi1530_epos], true);

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
    latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{#Xi(1530)^{0}}");
    legendRatio->Clear();
    legendRatio->AddEntry(gRatioKstarXiStar, "pp, #sqrt{s} = 13.6 TeV", "p");
    legendRatio->Draw();
    legendRatio2->Draw();
    cRatioKstarXiStar->SaveAs("Plots/Ratio_KstarXiStar_Run3.png");

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
    // cRatioPhiXiStar->SaveAs("Plots/Ratio_PhiXiStar_Run3.png");

    // */

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
    // cRatioKaonPion->SaveAs("Plots/Ratio_KaonPion_Run3.png");

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
    //     cRatioProtonPion->SaveAs("Plots/Ratio_ProtonPion_Run3.png");
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