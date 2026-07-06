#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

TFile *OpenFile(const string &path);
TH1D *GetHisto(TFile *f, const string &name);
void ScaleGraph(TGraph *gr, double scale);

void makeGraphXaxisCube(TGraph *gr)
{
    if (!gr)
        return;

    for (int i = 0; i < gr->GetN(); ++i)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        gr->SetPoint(i, x * x * x, y);
    }

    if (auto *ge = dynamic_cast<TGraphErrors *>(gr))
    {
        for (int i = 0; i < ge->GetN(); ++i)
        {
            ge->SetPointError(i, ge->GetErrorX(i) * ge->GetErrorX(i) * ge->GetErrorX(i), ge->GetErrorY(i));
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

void compareEPOS()
{
    string KstarPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/";
    TFile *fKstar = OpenFile(KstarPath + "Results.root");
    TFile *fPikp = OpenFile("ConversionCodes/ppRun3_PiKpHEP.root");
    TFile *fEPOSHyunji = OpenFile("ModelRootFiles/EPOS_pp13TeV_rhoKstar_Hyunji2.root");
    TFile *fEPOS = OpenFile("ModelRootFiles/EPOS_finalQA_INELgt0Correct.root");
    TFile *fEPOSHyperloop = OpenFile("ModelRootFiles/ModelResults.root");

    TGraphErrors *gYieldPionData = GetGraph(fPikp, "gPion_MeanYield_stat");
    TGraphErrors *gYieldKaonData = GetGraph(fPikp, "gKaon_MeanYield_stat");
    TGraphErrors *gYieldProtonData = GetGraph(fPikp, "gProton_MeanYield_stat");

    TGraphErrors *gdNdyKstarHloop = GetGraph(fEPOSHyperloop, "EPOS_Hydro/Kstar/gMeanYield_Kstar");
    TGraphErrors *gMPtKstarHloop = GetGraph(fEPOSHyperloop, "EPOS_Hydro/Kstar/gMeanpT_Kstar");
    TGraphErrors *gdNdyPionHloop = GetGraph(fEPOSHyperloop, "EPOS_Hydro/Pion/gMeanYield_Pion");
    TGraphErrors *gdNdyKaonHloop = GetGraph(fEPOSHyperloop, "EPOS_Hydro/Kaon/gMeanYield_Kaon");
    TGraphErrors *gdNdyProtonHloop = GetGraph(fEPOSHyperloop, "EPOS_Hydro/Proton/gMeanYield_Proton");

    TGraphErrors *gMPtPionData = GetGraph(fPikp, "gPion_MeanpT_stat");
    TGraphErrors *gMPtKaonData = GetGraph(fPikp, "gKaon_MeanpT_stat");
    TGraphErrors *gMPtProtonData = GetGraph(fPikp, "gProton_MeanpT_stat");

    TGraphErrors *gdNdyKstarHyunji_NoUrQMD = GetGraph(fEPOSHyunji, "c1_dNdy/Kstar892/UrQMD_OFF");
    TGraphErrors *gdNdyKstarHyunji_ITY80 = GetGraph(fEPOSHyunji, "c1_dNdy/Kstar892/UrQMD_reg");
    TGraphErrors *gdNdyKstarHyunji_ITY81 = GetGraph(fEPOSHyunji, "c1_dNdy/Kstar892/UrQMD_reg_res");

    TGraphErrors *gdNdyKstarEPOS_NoUrQMD = GetGraph(fEPOS, "IST9/kstar_vs_mult");
    TGraphErrors *gdNdyKstarEPOS_ITY80 = GetGraph(fEPOS, "IST9_ITY80/kstar_vs_mult");
    TGraphErrors *gdNdyKstarEPOS_ITY81 = GetGraph(fEPOS, "IST9_ITY81/kstar_vs_mult");

    TGraphErrors *gMPtKstarHyunji_NoUrQMD = GetGraph(fEPOSHyunji, "c2_meanpT/Kstar892/UrQMD_OFF");
    TGraphErrors *gMPtKstarHyunji_ITY80 = GetGraph(fEPOSHyunji, "c2_meanpT/Kstar892/UrQMD_reg");
    TGraphErrors *gMPtKstarHyunji_ITY81 = GetGraph(fEPOSHyunji, "c2_meanpT/Kstar892/UrQMD_reg_res");

    TGraphErrors *gMPtKstarEPOS_NoUrQMD = GetGraph(fEPOS, "IST9/meanpt_kstar_vs_mult");
    TGraphErrors *gMPtKstarEPOS_ITY80 = GetGraph(fEPOS, "IST9_ITY80/meanpt_kstar_vs_mult");
    TGraphErrors *gMPtKstarEPOS_ITY81 = GetGraph(fEPOS, "IST9_ITY81/meanpt_kstar_vs_mult");

    TGraphErrors *gdNdyKaonHyunji_NoUrQMD = GetGraph(fEPOSHyunji, "c1_dNdy/kaon/UrQMD_OFF");
    TGraphErrors *gdNdyKaonHyunji = GetGraph(fEPOSHyunji, "c1_dNdy/kaon/UrQMD_ON");
    TGraphErrors *gdNdyPionHyunji_NoUrQMD = GetGraph(fEPOSHyunji, "c1_dNdy/pion/UrQMD_OFF");
    TGraphErrors *gdNdyPionHyunji = GetGraph(fEPOSHyunji, "c1_dNdy/pion/UrQMD_ON");
    TGraphErrors *gdNdyProtonHyunji_NoUrQMD = GetGraph(fEPOSHyunji, "c1_dNdy/proton/UrQMD_OFF");
    TGraphErrors *gdNdyProtonHyunji = GetGraph(fEPOSHyunji, "c1_dNdy/proton/UrQMD_ON");

    TGraphErrors *gdNdyKaon_IST0 = GetGraph(fEPOS, "IST0/kaon_vs_mult");
    TGraphErrors *gdNdyPion_IST0 = GetGraph(fEPOS, "IST0/pion_vs_mult");
    TGraphErrors *gdNdyProton_IST0 = GetGraph(fEPOS, "IST0/proton_vs_mult");
    ScaleGraph(gdNdyKaon_IST0, 0.5);
    ScaleGraph(gdNdyPion_IST0, 0.5);
    ScaleGraph(gdNdyProton_IST0, 0.5);

    TGraphErrors *gMYieldKstar[2], *gMPtKstar[2];
    for (int i = 0; i < 2; i++)
    {
        string suffix;
        if (i == 0)
            suffix = "_stat";
        else if (i == 1)
            suffix = "_sys";

        gMYieldKstar[i] = GetGraph(fKstar, Form("gMeanYieldRun3%s", suffix.c_str()));
        gMPtKstar[i] = GetGraph(fKstar, Form("gMeanpTRun3%s", suffix.c_str()));
    }

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

    ScaleGraph(gdNdyKstarEPOS_NoUrQMD, 0.5);
    ScaleGraph(gdNdyKstarEPOS_ITY80, 0.5);
    gdNdyKstarEPOS_NoUrQMD->SetLineColor(kBlue + 1);
    gdNdyKstarEPOS_NoUrQMD->SetLineWidth(2);
    gdNdyKstarEPOS_NoUrQMD->Draw("L same");
    gdNdyKstarEPOS_ITY80->SetLineColor(kBlue + 1);
    gdNdyKstarEPOS_ITY80->SetLineWidth(2);
    gdNdyKstarEPOS_ITY80->SetLineStyle(2);
    gdNdyKstarEPOS_ITY80->Draw("L same");

    makeGraphXaxisCube(gdNdyKstarHyunji_NoUrQMD);
    ScaleGraph(gdNdyKstarHyunji_NoUrQMD, 0.5);
    gdNdyKstarHyunji_NoUrQMD->SetLineColor(kGreen + 1);
    gdNdyKstarHyunji_NoUrQMD->SetLineWidth(2);
    gdNdyKstarHyunji_NoUrQMD->Draw("L same");
    makeGraphXaxisCube(gdNdyKstarHyunji_ITY80);
    ScaleGraph(gdNdyKstarHyunji_ITY80, 0.5);
    gdNdyKstarHyunji_ITY80->SetLineColor(kGreen + 1);
    gdNdyKstarHyunji_ITY80->SetLineWidth(2);
    gdNdyKstarHyunji_ITY80->SetLineStyle(2);
    gdNdyKstarHyunji_ITY80->Draw("L same");
    makeGraphXaxisCube(gdNdyKstarHyunji_ITY81);
    ScaleGraph(gdNdyKstarHyunji_ITY81, 0.5);
    gdNdyKstarHyunji_ITY81->SetLineColor(kRed - 5);
    gdNdyKstarHyunji_ITY81->SetLineWidth(3);
    gdNdyKstarHyunji_ITY81->SetLineStyle(6);
    // gdNdyKstarHyunji_ITY81->Draw("L same");

    gdNdyKstarHloop->SetLineColor(kMagenta);
    gdNdyKstarHloop->SetLineWidth(2);
    gdNdyKstarHloop->Draw("L same");

    TLegend *legdNdyKstar = new TLegend(0.17, 0.73, 0.52, 0.95);
    legdNdyKstar->SetBorderSize(0);
    legdNdyKstar->SetFillStyle(0);
    legdNdyKstar->SetTextSize(0.027);
    legdNdyKstar->AddEntry(gMYieldKstar[0], "ALICE Run 3", "lp");
    legdNdyKstar->AddEntry(gdNdyKstarEPOS_NoUrQMD, "EPOS (No UrQMD) Sawan", "l");
    legdNdyKstar->AddEntry(gdNdyKstarEPOS_ITY80, "EPOS (UrQMD) Sawan", "l");
    legdNdyKstar->AddEntry(gdNdyKstarHyunji_NoUrQMD, "EPOS (No UrQMD) Hyunji", "l");
    legdNdyKstar->AddEntry(gdNdyKstarHyunji_ITY80, "EPOS (UrQMD) Hyunji (reg)", "l");
    // legdNdyKstar->AddEntry(gdNdyKstarHyunji_ITY81, "EPOS (UrQMD) Hyunji (reg+res)", "l");
    legdNdyKstar->AddEntry(gdNdyKstarHloop, "EPOS Hyperloop", "l");
    legdNdyKstar->Draw();
    cdNdyKstar->SaveAs("Plots/dNdyKstar_Comparison.png");

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

    gMPtKstarEPOS_NoUrQMD->SetLineColor(kBlue + 1);
    gMPtKstarEPOS_NoUrQMD->SetLineWidth(2);
    gMPtKstarEPOS_NoUrQMD->Draw("L same");
    gMPtKstarEPOS_ITY80->SetLineColor(kBlue + 1);
    gMPtKstarEPOS_ITY80->SetLineWidth(2);
    gMPtKstarEPOS_ITY80->SetLineStyle(2);
    gMPtKstarEPOS_ITY80->Draw("L same");

    makeGraphXaxisCube(gMPtKstarHyunji_NoUrQMD);
    gMPtKstarHyunji_NoUrQMD->SetLineColor(kGreen + 1);
    gMPtKstarHyunji_NoUrQMD->SetLineWidth(2);
    gMPtKstarHyunji_NoUrQMD->Draw("L same");
    makeGraphXaxisCube(gMPtKstarHyunji_ITY80);
    gMPtKstarHyunji_ITY80->SetLineColor(kGreen + 1);
    gMPtKstarHyunji_ITY80->SetLineWidth(2);
    gMPtKstarHyunji_ITY80->SetLineStyle(2);
    gMPtKstarHyunji_ITY80->Draw("L same");
    makeGraphXaxisCube(gMPtKstarHyunji_ITY81);
    gMPtKstarHyunji_ITY81->SetLineColor(kRed - 5);
    gMPtKstarHyunji_ITY81->SetLineWidth(3);
    gMPtKstarHyunji_ITY81->SetLineStyle(6);
    gMPtKstarHyunji_ITY81->Draw("L same");
    legdNdyKstar->Draw();
    cMeanPtKstar->SaveAs("Plots/MeanPtKstar_Comparison.png");

    TCanvas *cYieldPion = new TCanvas("cYieldPion", "cYieldPion", 720, 720);
    SetCanvasStyle(cYieldPion, 0.15, 0.03, 0.03, 0.15);

    gYieldPionData->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gYieldPionData->GetYaxis()->SetTitle("dN/dy");
    gYieldPionData->GetXaxis()->SetLimits(0, 27);
    gYieldPionData->GetYaxis()->SetRangeUser(0, 20);
    gYieldPionData->SetLineWidth(3);
    gYieldPionData->SetMarkerColor(kBlack);
    gYieldPionData->SetLineColor(kBlack);
    gYieldPionData->Draw("APE");
    gdNdyPion_IST0->SetLineColor(kBlue + 1);
    gdNdyPion_IST0->SetLineWidth(2);
    gdNdyPion_IST0->Draw("L same");
    makeGraphXaxisCube(gdNdyPionHyunji_NoUrQMD);
    gdNdyPionHyunji_NoUrQMD->SetLineColor(kGreen + 1);
    gdNdyPionHyunji_NoUrQMD->SetLineWidth(2);
    gdNdyPionHyunji_NoUrQMD->Draw("L same");
    makeGraphXaxisCube(gdNdyPionHyunji);
    gdNdyPionHyunji->SetLineColor(kRed - 5);
    gdNdyPionHyunji->SetLineWidth(3);
    gdNdyPionHyunji->SetLineStyle(6);
    gdNdyPionHyunji->Draw("L same");

    TLegend *legYieldPion = new TLegend(0.17, 0.73, 0.52, 0.95);
    legYieldPion->SetBorderSize(0);
    legYieldPion->SetFillStyle(0);
    legYieldPion->SetTextSize(0.027);
    legYieldPion->AddEntry(gYieldPionData, "pp 13.6 TeV, Pion", "lp");
    legYieldPion->AddEntry(gdNdyPion_IST0, "EPOS (No UrQMD) Sawan", "l");
    legYieldPion->AddEntry(gdNdyPionHyunji_NoUrQMD, "EPOS (No UrQMD) Hyunji", "l");
    legYieldPion->AddEntry(gdNdyPionHyunji, "EPOS (UrQMD) Hyunji", "l");
    legYieldPion->Draw();
    cYieldPion->SaveAs("Plots/dNdyPion_Comparison.png");

    TCanvas *cYieldKaon = new TCanvas("cYieldKaon", "cYieldKaon", 720, 720);
    SetCanvasStyle(cYieldKaon, 0.15, 0.03, 0.03, 0.15);
    gYieldKaonData->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gYieldKaonData->GetYaxis()->SetTitle("dN/dy");
    gYieldKaonData->GetXaxis()->SetLimits(0, 27);
    gYieldKaonData->GetYaxis()->SetRangeUser(0, 2);
    gYieldKaonData->SetLineWidth(3);
    gYieldKaonData->SetMarkerColor(kBlack);
    gYieldKaonData->SetLineColor(kBlack);
    gYieldKaonData->Draw("APE");
    gdNdyKaon_IST0->SetLineColor(kBlue + 1);
    gdNdyKaon_IST0->SetLineWidth(2);
    gdNdyKaon_IST0->Draw("L same");
    makeGraphXaxisCube(gdNdyKaonHyunji_NoUrQMD);
    gdNdyKaonHyunji_NoUrQMD->SetLineColor(kGreen + 1);
    gdNdyKaonHyunji_NoUrQMD->SetLineWidth(2);
    gdNdyKaonHyunji_NoUrQMD->Draw("L same");
    makeGraphXaxisCube(gdNdyKaonHyunji);
    gdNdyKaonHyunji->SetLineColor(kRed - 5);
    gdNdyKaonHyunji->SetLineWidth(3);
    gdNdyKaonHyunji->SetLineStyle(6);
    gdNdyKaonHyunji->Draw("L same");

    TLegend *legYieldKaon = new TLegend(0.17, 0.73, 0.52, 0.95);
    legYieldKaon->SetBorderSize(0);
    legYieldKaon->SetFillStyle(0);
    legYieldKaon->SetTextSize(0.027);
    legYieldKaon->AddEntry(gYieldKaonData, "pp 13.6 TeV, Kaon", "lp");
    legYieldKaon->AddEntry(gdNdyKaon_IST0, "EPOS (No UrQMD) Sawan", "l");
    legYieldKaon->AddEntry(gdNdyKaonHyunji_NoUrQMD, "EPOS (No UrQMD) Hyunji", "l");
    legYieldKaon->AddEntry(gdNdyKaonHyunji, "EPOS (UrQMD) Hyunji", "l");
    legYieldKaon->Draw();
    cYieldKaon->SaveAs("Plots/dNdyKaon_Comparison.png");

    TCanvas *cYieldProton = new TCanvas("cYieldProton", "cYieldProton", 720, 720);
    SetCanvasStyle(cYieldProton, 0.15, 0.03, 0.03, 0.15);
    gYieldProtonData->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gYieldProtonData->GetYaxis()->SetTitle("dN/dy");
    gYieldProtonData->GetXaxis()->SetLimits(0, 27);
    gYieldProtonData->GetYaxis()->SetRangeUser(0, 1);
    gYieldProtonData->SetLineWidth(3);
    gYieldProtonData->SetMarkerColor(kBlack);
    gYieldProtonData->SetLineColor(kBlack);
    gYieldProtonData->Draw("APE");
    gdNdyProton_IST0->SetLineColor(kBlue + 1);
    gdNdyProton_IST0->SetLineWidth(2);
    gdNdyProton_IST0->Draw("L same");
    makeGraphXaxisCube(gdNdyProtonHyunji_NoUrQMD);
    gdNdyProtonHyunji_NoUrQMD->SetLineColor(kGreen + 1);
    gdNdyProtonHyunji_NoUrQMD->SetLineWidth(2);
    gdNdyProtonHyunji_NoUrQMD->Draw("L same");
    makeGraphXaxisCube(gdNdyProtonHyunji);
    gdNdyProtonHyunji->SetLineColor(kRed - 5);
    gdNdyProtonHyunji->SetLineWidth(3);
    gdNdyProtonHyunji->SetLineStyle(6);
    gdNdyProtonHyunji->Draw("L same");

    TLegend *legYieldProton = new TLegend(0.17, 0.73, 0.52, 0.95);
    legYieldProton->SetBorderSize(0);
    legYieldProton->SetFillStyle(0);
    legYieldProton->SetTextSize(0.027);
    legYieldProton->AddEntry(gYieldProtonData, "pp 13.6 TeV, Proton", "lp");
    legYieldProton->AddEntry(gdNdyProton_IST0, "EPOS (No UrQMD) Sawan", "l");
    legYieldProton->AddEntry(gdNdyProtonHyunji_NoUrQMD, "EPOS (No UrQMD) Hyunji", "l");
    legYieldProton->AddEntry(gdNdyProtonHyunji, "EPOS (UrQMD) Hyunji", "l");
    legYieldProton->Draw();
    cYieldProton->SaveAs("Plots/dNdyProton_Comparison.png");
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