#include <iostream>
#include <iomanip>
#include "../src/style.h"

TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false);
TGraphErrors *GetGraph(TFile *f, const string &name);

void setStyle(TGraphErrors *gr, int color, int style)
{
    gr->SetLineColor(color);
    gr->SetLineStyle(style);
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

    TGraphErrors *gMPtKstar = GetGraph(fKstar, "gMeanpTRun3_sys");
    TGraphErrors *gMYieldKstar = GetGraph(fKstar, "gMeanYieldRun3_sys");

    TGraphErrors *gMPtPhi = GetGraph(fPhi, "gPhi_MeanpT_sys");
    TGraphErrors *gMPtPion = GetGraph(fPion, "gMeanpTRun3");
    TGraphErrors *gMPtProton = GetGraph(fProton, "gMeanpTRun3");
    TGraphErrors *gMPtKaon = GetGraph(fKaon, "gMeanpTRun3");
    TGraphErrors *gMPtChKstar = GetGraph(fChKstar, "gMeanpTRun3");
    TGraphErrors *gMPtXiStar = GetGraph(fXiStar, "gMeanpTRun3");
    TGraphErrors *gMPtKshort = GetGraph(fKshort, "gMeanpTRun3");

    TGraphErrors *gMYieldPhi = GetGraph(fPhi, "gPhi_MeanYield_sys");
    TGraphErrors *gMYieldPion = GetGraph(fPion, "gMeanYieldRun3");
    TGraphErrors *gMYieldProton = GetGraph(fProton, "gMeanYieldRun3");
    TGraphErrors *gMYieldKaon = GetGraph(fKaon, "gMeanYieldRun3");
    TGraphErrors *gMYieldChKstar = GetGraph(fChKstar, "gMeanYieldRun3");
    TGraphErrors *gMYieldXiStar = GetGraph(fXiStar, "gMeanYieldRun3");
    TGraphErrors *gMYieldKshort = GetGraph(fKshort, "gMeanYieldRun3");

    //======================================================
    //    ===========EPOS local model files===========
    //======================================================
    TFile *fEPOS = new TFile("ModelRootFiles/EPOS_finalQA.root", "read");
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

    //======================
    // Resonances (IST9)
    //======================
    TGraphErrors *gMYieldKstarEPOS_IST9 = GetGraph(fEPOS, "IST9/kstar_vs_mult");
    TGraphErrors *gMYieldPhiEPOS_IST9 = GetGraph(fEPOS, "IST9/phi_vs_mult");
    TGraphErrors *gMYieldXiStarEPOS_IST9 = GetGraph(fEPOS, "IST9/xistar_vs_mult");
    TGraphErrors *gMYieldChargedKstarEPOS_IST9 = GetGraph(fEPOS, "IST9/chargedkstar_vs_mult");

    TGraphErrors *gMYieldKstarEPOS_IST9_ITY80 = GetGraph(fEPOS, "IST9_ITY80/kstar_vs_mult");
    TGraphErrors *gMYieldKstarEPOS_IST9_ITY81 = GetGraph(fEPOS, "IST9_ITY81/kstar_vs_mult");
    TGraphErrors *gMYieldKstarEPOS_IST1 = GetGraph(fEPOS, "IST1/kstar_vs_mult"); // NoUrQMD
    TGraphErrors *gMYieldKstarEPOS_IST6_orig = GetGraph(fEPOS, "IST6/kstar_vs_mult");
    TGraphErrors *gMYieldKstarEPOS_IST6 = new TGraphErrors(gMYieldKstarEPOS_IST6_orig->GetN());

    for (int ibin = 0; ibin < gMYieldKstarEPOS_IST6_orig->GetN(); ibin++)
    {
        double x, y1, y2;
        gMYieldKstarEPOS_IST6_orig->GetPoint(ibin, x, y1);
        gMYieldKstarEPOS_IST1->GetPoint(ibin, x, y2);
        double y3 = y1 + y2;
        y3 = y3 * 0.666; // BR correction in not needed for UrQMD Off case
        // double y2 = y;
        gMYieldKstarEPOS_IST6->SetPoint(ibin, x, y3);
        gMYieldKstarEPOS_IST6->SetPointError(ibin, 0, 0);
    }

    TGraphErrors *gMeanPtKstarEPOS_IST9 = GetGraph(fEPOS, "IST9/meanpt_kstar_vs_mult");
    TGraphErrors *gMeanPtKstarEPOS_IST9_ITY80 = GetGraph(fEPOS, "IST9_ITY80/meanpt_kstar_vs_mult");
    TGraphErrors *gMeanPtKstarEPOS_IST9_ITY81 = GetGraph(fEPOS, "IST9_ITY81/meanpt_kstar_vs_mult");
    TGraphErrors *gMeanPtKstarEPOS_IST6 = GetGraph(fEPOS, "IST6/meanpt_kstar_vs_mult"); // NoUrQMD

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
        {kBlack, 2},      // EPOS
        {kBlue, 2},       // Pythia CR
        {kRed, 2},        // Pythia Monash
        {kGreen + 2, 2},  // Pythia Ropes
        {kOrange + 2, 2}, // Pythia Shoving
        {kMagenta + 1, 2} // Pythia Monash Rescattering
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
        }
    }

    ////======================================================
    //    ================Plots=====================
    //======================================================

    vector<int> modelsToPlot = {
        kPythiaMonash,
        kPythiaShoving};

    TCanvas *cdNdyKstar = new TCanvas("cdNdyKstar", "cdNdyKstar", 720, 720);
    SetCanvasStyle(cdNdyKstar, 0.15, 0.03, 0.03, 0.15);
    gMYieldKstar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldKstar->GetYaxis()->SetTitle("dN/dy");
    gMYieldKstar->GetXaxis()->SetLimits(0, 27);
    gMYieldKstar->GetYaxis()->SetRangeUser(0, 0.8);
    gMYieldKstar->Draw("AP");
    gMYieldKstarEPOS_IST9->SetLineStyle(2);
    gMYieldKstarEPOS_IST9->Draw("l same");
    gMYieldKstarEPOS_IST9_ITY80->SetLineColor(kBlue);
    gMYieldKstarEPOS_IST9_ITY80->SetLineStyle(2);
    gMYieldKstarEPOS_IST9_ITY80->Draw("l same");
    gMYieldKstarEPOS_IST9_ITY81->SetLineColor(kGreen + 2);
    gMYieldKstarEPOS_IST9_ITY81->SetLineStyle(4);
    gMYieldKstarEPOS_IST9_ITY81->Draw("l same");
    gMYieldKstarEPOS_IST6->SetLineColor(kMagenta);
    gMYieldKstarEPOS_IST6->SetLineStyle(2);
    SetGraphErrorStyle(gMYieldKstarEPOS_IST6);
    // gMYieldKstarEPOS_IST6->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMYield[model][kKstar]->Draw("l same");
    }

    TLegend *legend = new TLegend(0.2, 0.72, 0.5, 0.92);
    SetLegendStyle(legend);
    legend->SetTextSize(0.027);
    legend->SetHeader("#frac{K^{*0} + #bar{K}^{*0}}{2}", "C");
    legend->AddEntry(gMYieldKstar, "Data", "P");

    // legend->AddEntry(gMYieldKstarEPOS_IST9, "Status = 9", "L");
    // legend->AddEntry(gMYieldKstarEPOS_IST9_ITY80, "Status = 9 && ity = 80", "L");
    // legend->AddEntry(gMYieldKstarEPOS_IST9_ITY81, "Status = 9 && ity = 81", "L");
    // legend->AddEntry(gMYieldKstarEPOS_IST6, "Status = 7", "L");

    legend->AddEntry(gMYieldKstarEPOS_IST9, "EPOS UrQMD OFF", "L");
    legend->AddEntry(gMYieldKstarEPOS_IST9_ITY80, "EPOS UrQMD ON", "L");
    legend->AddEntry(gMYieldKstarEPOS_IST9_ITY81, "EPOS Rescattering", "L");
    // legend->AddEntry(gMYieldKstarEPOS_IST6, "EPOS UrQMD OFF", "L");

    for (auto model : modelsToPlot)
    {
        legend->AddEntry(gMYield[model][kKstar], modelLabel[model], "L");
    }

    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    cdNdyKstar->SaveAs("Plots/MeanYield_Kstar_EPOS_UrQMDON.png");

    TCanvas *cMeanPtKstar = new TCanvas("cMeanPtKstar", "cMeanPtKstar", 720, 720);
    SetCanvasStyle(cMeanPtKstar, 0.15, 0.03, 0.03, 0.15);
    SetGraphErrorStyle(gMPtKstar);
    SetGraphErrorStyle(gMeanPtKstarEPOS_IST9);
    SetGraphErrorStyle(gMeanPtKstarEPOS_IST9_ITY80);
    SetGraphErrorStyle(gMeanPtKstarEPOS_IST9_ITY81);
    SetGraphErrorStyle(gMeanPtKstarEPOS_IST6);
    gMPtKstar->SetTitle(0);
    gMPtKstar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMPtKstar->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    gMPtKstar->GetXaxis()->SetLimits(0, 27);
    gMPtKstar->GetYaxis()->SetRangeUser(0.0, 1.89);
    gMPtKstar->Draw("AP");
    gMeanPtKstarEPOS_IST9->SetLineStyle(2);
    gMeanPtKstarEPOS_IST9->Draw("l same");
    gMeanPtKstarEPOS_IST9_ITY80->SetLineColor(kBlue);
    gMeanPtKstarEPOS_IST9_ITY80->SetLineStyle(2);
    gMeanPtKstarEPOS_IST9_ITY80->Draw("l same");
    gMeanPtKstarEPOS_IST9_ITY81->SetLineColor(kGreen + 2);
    gMeanPtKstarEPOS_IST9_ITY81->SetLineStyle(4);
    gMeanPtKstarEPOS_IST9_ITY81->Draw("l same");
    gMeanPtKstarEPOS_IST6->SetLineColor(kMagenta);
    gMeanPtKstarEPOS_IST6->SetLineStyle(3);
    // gMeanPtKstarEPOS_IST6->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMeanPt[model][kKstar]->Draw("l same");
    }
    legend->Draw();
    cMeanPtKstar->SaveAs("Plots/MeanPt_Kstar_EPOS_UrQMDON.png");

    //====================================================
    //   ================K*/K Ratio===================
    //====================================================

    TCanvas *cRatioKstarKaon = new TCanvas("cRatioKstarKaon", "cRatioKstarKaon", 720, 720);
    SetCanvasStyle(cRatioKstarKaon, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarKaon = MakeRatio(gMYieldKstar, gMYieldKaon, "false");
    gRatioKstarKaon->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarKaon->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarKaon);
    gRatioKstarKaon->GetYaxis()->SetRangeUser(0.0, 0.52);
    gRatioKstarKaon->GetXaxis()->SetLimits(0, 27);
    gRatioKstarKaon->Draw("APE");

    TGraphErrors *gRatioKstarKa_IST9 = MakeRatio(gMYieldKstarEPOS_IST9, gMYieldKaonEPOS_IST0, "true");
    TGraphErrors *gRatioKstarKa_IST9_ITY80 = MakeRatio(gMYieldKstarEPOS_IST9_ITY80, gMYieldKaonEPOS_IST0, "true");
    TGraphErrors *gRatioKstarKa_IST9_ITY81 = MakeRatio(gMYieldKstarEPOS_IST9_ITY81, gMYieldKaonEPOS_IST0, "true");
    TGraphErrors *gRatioKstarKa_IST6 = MakeRatio(gMYieldKstarEPOS_IST6, gMYieldKaonEPOS_IST0, "true");

    gRatioKstarKa_IST9->SetLineColor(kRed);
    gRatioKstarKa_IST9->SetLineStyle(2);
    gRatioKstarKa_IST9->Draw("l same");
    gRatioKstarKa_IST9_ITY80->SetLineColor(kBlue);
    gRatioKstarKa_IST9_ITY80->SetLineStyle(2);
    gRatioKstarKa_IST9_ITY80->Draw("l same");
    gRatioKstarKa_IST9_ITY81->SetLineColor(kGreen + 2);
    gRatioKstarKa_IST9_ITY81->SetLineStyle(4);
    gRatioKstarKa_IST9_ITY81->Draw("l same");
    gRatioKstarKa_IST6->SetLineColor(kMagenta);
    gRatioKstarKa_IST6->SetLineStyle(5);
    // gRatioKstarKa_IST6->Draw("l same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);

    TLegend *legendRatio = new TLegend(0.2, 0.72, 0.5, 0.92);
    SetLegendStyle(legendRatio);
    legendRatio->SetTextSize(0.027);
    legendRatio->AddEntry(gRatioKstarKaon, "Data", "p");
    // legendRatio->AddEntry(gRatioKstarKa_IST6, "EPOS UrQMD OFF", "l");
    legendRatio->AddEntry(gRatioKstarKa_IST9, "EPOS UrQMD OFF", "l");
    legendRatio->AddEntry(gRatioKstarKa_IST9_ITY80, "EPOS UrQMD ON", "l");
    legendRatio->AddEntry(gRatioKstarKa_IST9_ITY81, "EPOS Rescattering", "l");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarKaModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKaon], true);
        setStyle(gRatioKstarKaModel, modelStyle[model].color, modelStyle[model].style);
        gRatioKstarKaModel->Draw("l same");
        legendRatio->AddEntry(gRatioKstarKaModel, modelLabel[model], "l");
    }

    legendRatio->SetBorderSize(0);
    latex.DrawLatex(0.6, 0.8, "#frac{K^{*0} + #bar{K}^{*0}}{2K}");
    legendRatio->Draw();
    cRatioKstarKaon->SaveAs("Plots/Ratio_KstarKaon_Run3.png");

    //====================================================
    //   ================K*/Kshort Ratio===================
    //====================================================
    TCanvas *cRatioKstarKshort = new TCanvas("cRatioKstarKshort", "cRatioKstarKshort", 720, 720);
    SetCanvasStyle(cRatioKstarKshort, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarKshort = MakeRatio(gMYieldKstar, gMYieldKshort, "false");
    gRatioKstarKshort->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarKshort->GetYaxis()->SetTitle("dN/dy");
    SetGraphErrorStyle(gRatioKstarKshort);
    gRatioKstarKshort->GetYaxis()->SetRangeUser(0.0, 0.52);
    gRatioKstarKshort->GetXaxis()->SetLimits(0, 27);
    gRatioKstarKshort->Draw("APE");

    TGraphErrors *gRatioKstarKshort_IST9 = MakeRatio(gMYieldKstarEPOS_IST9, gMYieldKshortEPOS_IST0, "true");
    TGraphErrors *gRatioKstarKshort_IST9_ITY80 = MakeRatio(gMYieldKstarEPOS_IST9_ITY80, gMYieldKshortEPOS_IST0, "true");
    TGraphErrors *gRatioKstarKshort_IST9_ITY81 = MakeRatio(gMYieldKstarEPOS_IST9_ITY81, gMYieldKshortEPOS_IST0, "true");
    TGraphErrors *gRatioKstarKshort_IST6 = MakeRatio(gMYieldKstarEPOS_IST6, gMYieldKshortEPOS_IST0, "true");

    gRatioKstarKshort_IST9->SetLineColor(kRed);
    gRatioKstarKshort_IST9->SetLineStyle(2);
    gRatioKstarKshort_IST9->Draw("l same");
    gRatioKstarKshort_IST9_ITY80->SetLineColor(kBlue);
    gRatioKstarKshort_IST9_ITY80->SetLineStyle(2);
    gRatioKstarKshort_IST9_ITY80->Draw("l same");
    gRatioKstarKshort_IST9_ITY81->SetLineColor(kGreen + 2);
    gRatioKstarKshort_IST9_ITY81->SetLineStyle(2);
    gRatioKstarKshort_IST9_ITY81->Draw("l same");
    gRatioKstarKshort_IST6->SetLineColor(kMagenta);
    gRatioKstarKshort_IST6->SetLineStyle(5);
    // gRatioKstarKshort_IST6->Draw("l same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarKshortModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKshort], true);
        setStyle(gRatioKstarKshortModel, modelStyle[model].color, modelStyle[model].style);
        gRatioKstarKshortModel->Draw("l same");
        // legendRatio->AddEntry(gRatioKstarKshortModel, modelLabel[model], "l");
    }

    latex.DrawLatex(0.6, 0.8, "#frac{K^{*0} + #bar{K}^{*0}}{2K_{S}^{0}}");
    legendRatio->Draw();
    cRatioKstarKshort->SaveAs("Plots/Ratio_KstarKshort_Run3.png");
}

TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false)
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

        double yRatio = (yDenominator != 0) ? yNumerator / yDenominator : 0;
        ratio->SetPoint(i, xNumerator, yRatio);

        double xErrorNumerator = numerator->GetErrorX(i);
        double numeratorError = numerator->GetErrorY(i);
        double denominatorError = denominator->GetErrorY(i);
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