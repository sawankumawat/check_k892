#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

TFile *OpenFile(const string &path);
TH1D *GetHisto(TFile *f, const string &name);
void ScaleGraph(TGraph *gr, double scale);
TGraphErrors *GetGraph(TFile *f, const string &name);
void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax);

void compareYieldSimple()
{
    string KstarPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/";
    TFile *fKstar = OpenFile(KstarPath + "Results.root");

    TFile *fPion = OpenFile("PiKp_Run3_Results/Sawan/Pi_results.root");
    TFile *fProton = OpenFile("PiKp_Run3_Results/Sawan/Pr_results.root");
    TFile *fKaon = OpenFile("PiKp_Run3_Results/Sawan/Ka_results.root");
    TFile *fPiKp = OpenFile("ConversionCodes/ppRun3_PiKpHEP.root");

    TGraphErrors *gMPtKstar[3], *gMYieldKstar[3], *gMPtPion[3], *gMPtProton[3], *gMPtKaon[3];
    TGraphErrors *gMYieldPion[3], *gMYieldProton[3], *gMYieldKaon[3];

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

        gMYieldPion[i] = GetGraph(fPion, Form("gMeanYieldRun3%s", suffix.c_str()));
        gMYieldKaon[i] = GetGraph(fKaon, Form("gMeanYieldRun3%s", suffix.c_str()));
        gMYieldProton[i] = GetGraph(fProton, Form("gMeanYieldRun3%s", suffix.c_str()));

        // gMYieldPion[i] = GetGraph(fPiKp, Form("gPion_MeanYield%s", suffix.c_str()));
        // gMYieldKaon[i] = GetGraph(fPiKp, Form("gKaon_MeanYield%s", suffix.c_str()));
        // gMYieldProton[i] = GetGraph(fPiKp, Form("gProton_MeanYield%s", suffix.c_str()));
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
    // enum Particle
    // {
    //     kKstar,
    //     kPhi,
    //     kPion,
    //     kKaon,
    //     kProton,
    //     kPionMinus,
    //     kKaonMinus,
    //     kAntiProton,
    //     kXi1530,
    //     kKshort,
    //     kKstarPM,
    //     kNParticles
    // };

    enum Particle
    {
        kKstar,
        kPion,
        kKaon,
        kProton,
        kPionMinus,
        kKaonMinus,
        kAntiProton,
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
        {kGreen + 2, 2}, // EPOS
        {kBlue + 1, 2},  // Pythia CR
        {kMagenta, 4},   // Pythia Monash
        {kCyan + 1, 7},  // Pythia Ropes
        {kRed + 1, 3},   // Pythia Shoving
        {kBlue + 1, 2}   // Pythia Monash Rescattering
    };

    const char *modelLabel[kNModels] = {
        "EPOS Hydro",
        "Pythia CR",
        "Pythia Monash",
        "Pythia Ropes",
        "Pythia Shoving",
        "Pythia Monash Rescattering"};

    TFile *ModelsHyperloop = new TFile("ModelRootFiles/ModelResults2.root", "read");
    if (ModelsHyperloop->IsZombie())
    {
        cout << "Error: Hyperloop model file not found" << endl;
        return;
    }

    string hyperloopModels[kNModels] = {"EPOS_Hydro", "Pythia_CR", "Pythia_Monash2", "Pythia_Ropes2", "Pythia_Shoving2", "Pythia_Monash_Rescattering"};
    // string particles[kNParticles] = {"Kstar", "Phi", "Pion", "Kaon", "Proton", "PionMinus", "KaonMinus", "AntiProton", "Xi1530", "Kshort", "KstarPM"};
    string particles[kNParticles] = {"Kstar", "Pion", "Kaon", "Proton", "PionMinus", "KaonMinus", "AntiProton", "Kshort", "KstarPM"};

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

            gMYield[iModel][iParticle]->SetLineWidth(3);
            gMeanPt[iModel][iParticle]->SetLineWidth(3);

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
            double yPionTotal = (yPion + yPionMinus) / 2.0;
            gMYield[iModel][kPion]->SetPoint(i, x, yPionTotal);

            // Kaon
            double yKaon, yKaonMinus;
            gMYield[iModel][kKaon]->GetPoint(i, x, yKaon);
            gMYield[iModel][kKaonMinus]->GetPoint(i, x, yKaonMinus);
            double yKaonTotal = (yKaon + yKaonMinus) / 2.0;
            gMYield[iModel][kKaon]->SetPoint(i, x, yKaonTotal);

            // Proton
            double yProton, yAntiProton;
            gMYield[iModel][kProton]->GetPoint(i, x, yProton);
            gMYield[iModel][kAntiProton]->GetPoint(i, x, yAntiProton);
            double yProtonTotal = (yProton + yAntiProton) / 2.0;
            gMYield[iModel][kProton]->SetPoint(i, x, yProtonTotal);
        }
    }

    // vector<int> modelsToPlot = {
    //     kEPOS_Hydro,
    //     kPythiaMonash,
    //     kPythiaMonashRescattering};

    vector<int> modelsToPlot = {
        kPythiaMonash,
        kPythiaRopes,
        kPythiaShoving};

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

    for (auto model : modelsToPlot)
    {
        ScaleGraph(gMYield[model][kKstar], 1.0 / (2.0));
        gMYield[model][kKstar]->Draw("l same");
    }

    TLegend *legend = new TLegend(0.18, 0.75, 0.48, 0.9);
    SetLegendStyle(legend);
    legend->SetTextSize(0.03);
    legend->AddEntry(gMYieldKstar[0], "pp, #sqrt{s} = 13.6 TeV", "P");
    legend->Draw();

    TLegend *legend2 = new TLegend(0.5, 0.72, 0.8, 0.92);
    SetLegendStyle(legend2);
    legend2->SetTextSize(0.027);

    for (auto model : modelsToPlot)
    {
        legend2->AddEntry(gMYield[model][kKstar], modelLabel[model], "L");
    }
    legend->Draw();
    legend2->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(22);
    latex.SetTextSize(0.027);
    latex.DrawLatex(0.28, 0.9, "K* (892)^{0}");
    cdNdyKstar->SaveAs("Plots/YieldCompareHyperloopModels/Kstar_dNdy.png");

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
    gMPtKstar[0]->GetYaxis()->SetRangeUser(0.52, 2.13);
    gMPtKstar[0]->SetMarkerColor(kRed);
    gMPtKstar[0]->SetLineColor(kRed);
    gMPtKstar[0]->SetLineWidth(3);
    gMPtKstar[0]->Draw("APE");
    gMPtKstar[1]->SetFillStyle(0);
    gMPtKstar[1]->SetLineColor(kRed);
    gMPtKstar[1]->SetLineWidth(3);
    gMPtKstar[1]->Draw("5 same");

    for (auto model : modelsToPlot)
    {
        gMeanPt[model][kKstar]->Draw("l same");
    }

    legend->Draw();
    legend2->Draw();
    latex.DrawLatex(0.28, 0.9, "K* (892)^{0}");
    cMeanPtKstar->SaveAs("Plots/YieldCompareHyperloopModels/Kstar_MeanPt.png");

    // //=============EPOS model ran locally================
    // TFile *fEPOSLocal = OpenFile("dNdy_vs_Nch_PercentileSliced.root");
    // TGraphErrors *gMYieldPionEPOS = GetGraph(fEPOSLocal, "gYieldVsNch_#pi");
    // TGraphErrors *gMYieldKaonEPOS = GetGraph(fEPOSLocal, "gYieldVsNch_K");
    // TGraphErrors *gMYieldProtonEPOS = GetGraph(fEPOSLocal, "gYieldVsNch_p");

    // SetGraphErrorStyle(gMYieldPionEPOS);
    // SetGraphErrorStyle(gMYieldKaonEPOS);
    // SetGraphErrorStyle(gMYieldProtonEPOS);

    //====================================================
    // ==================Pion yeild======================
    //====================================================
    TCanvas *cPionYield = new TCanvas("cPionYield", "cPionYield", 720, 720);
    SetCanvasStyle(cPionYield, 0.15, 0.03, 0.03, 0.15);
    gMYieldPion[0]->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldPion[0]->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    SetGraphErrorStyle(gMYieldPion[0]);
    gMYieldPion[0]->GetYaxis()->SetRangeUser(0.0, 15);
    gMYieldPion[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldPion[0]->SetMarkerColor(kRed);
    gMYieldPion[0]->SetLineColor(kRed);
    gMYieldPion[0]->Draw("APE");
    gMYieldPion[1]->SetFillStyle(0);
    gMYieldPion[1]->SetLineColor(kRed);
    gMYieldPion[1]->Draw("5 same");
    // ScaleGraph(gMYieldPionEPOS, 0.5);
    // gMYieldPionEPOS->SetLineColor(kGray + 3);
    // gMYieldPionEPOS->SetLineWidth(3);
    // gMYieldPionEPOS->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMYield[model][kPion]->Draw("l same");
    }

    legend->Draw();
    // legend2->AddEntry(gMYieldPionEPOS, "EPOS (local)", "L");
    legend2->Draw();
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.28, 0.9, "#pi");
    cPionYield->SaveAs("Plots/YieldCompareHyperloopModels/Pion_dNdy.png");

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
    // ScaleGraph(gMYieldKaonEPOS, 0.5);
    // gMYieldKaonEPOS->SetLineColor(kGray + 3);
    // gMYieldKaonEPOS->SetLineWidth(3);
    // gMYieldKaonEPOS->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMYield[model][kKaon]->Draw("l same");
    }

    legend->Draw();
    legend2->Draw();
    latex.DrawLatex(0.28, 0.9, "K");
    cKaonYield->SaveAs("Plots/YieldCompareHyperloopModels/Kaon_dNdy.png");

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
    // ScaleGraph(gMYieldProtonEPOS, 0.5);
    // gMYieldProtonEPOS->SetLineColor(kGray + 3);
    // gMYieldProtonEPOS->SetLineWidth(3);
    // gMYieldProtonEPOS->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMYield[model][kProton]->Draw("l same");
    }
    legend->Draw();
    legend2->Draw();
    latex.DrawLatex(0.28, 0.9, "p");
    cProtonYield->SaveAs("Plots/YieldCompareHyperloopModels/Proton_dNdy.png");
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

TH1D *GetHisto(TFile *f, const string &name)
{
    TH1D *histo = (TH1D *)f->Get(name.c_str());

    if (!histo || histo == nullptr)
    {
        cout << "Error: histo " << name << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetHistoQA(histo);
    histo->SetTitle(0);
    return histo;
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

TGraphErrors *GetGraph(TFile *f, const string &name)
{
    TGraphErrors *graph = (TGraphErrors *)f->Get(name.c_str());

    if (!graph || graph == nullptr)
    {
        cout << "Error: graph " << name << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetGraphErrorStyle(graph);
    graph->SetTitle(0);
    return graph;
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