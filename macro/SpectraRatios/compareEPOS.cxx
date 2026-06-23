#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

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
    TFile *fEPOSHyunji = OpenFile("ModelRootFiles/EPOS_pp13TeV_rhoKstar_Hyunji.root");
    TFile *fEPOS = OpenFile("ModelRootFiles/EPOS_finalQA_INELgt0Correct.root");

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
    // ScaleGraph(gdNdyKstarHyunji_NoUrQMD, 0.5);
    gdNdyKstarHyunji_NoUrQMD->SetLineColor(kGreen + 1);
    gdNdyKstarHyunji_NoUrQMD->SetLineWidth(2);
    gdNdyKstarHyunji_NoUrQMD->Draw("L same");
    makeGraphXaxisCube(gdNdyKstarHyunji_ITY80);
    // ScaleGraph(gdNdyKstarHyunji_ITY80, 0.5);
    gdNdyKstarHyunji_ITY80->SetLineColor(kGreen + 1);
    gdNdyKstarHyunji_ITY80->SetLineWidth(2);
    gdNdyKstarHyunji_ITY80->SetLineStyle(2);
    gdNdyKstarHyunji_ITY80->Draw("L same");
    makeGraphXaxisCube(gdNdyKstarHyunji_ITY81);
    // ScaleGraph(gdNdyKstarHyunji_ITY81, 0.5);
    gdNdyKstarHyunji_ITY81->SetLineColor(kRed - 5);
    gdNdyKstarHyunji_ITY81->SetLineWidth(3);
    gdNdyKstarHyunji_ITY81->SetLineStyle(6);
    gdNdyKstarHyunji_ITY81->Draw("L same");

    TLegend *legdNdyKstar = new TLegend(0.17, 0.73, 0.52, 0.95);
    legdNdyKstar->SetBorderSize(0);
    legdNdyKstar->SetFillStyle(0);
    legdNdyKstar->SetTextSize(0.027);
    legdNdyKstar->AddEntry(gMYieldKstar[0], "ALICE Run 3", "lp");
    legdNdyKstar->AddEntry(gdNdyKstarEPOS_NoUrQMD, "EPOS (No UrQMD) Sawan", "l");
    legdNdyKstar->AddEntry(gdNdyKstarEPOS_ITY80, "EPOS (UrQMD) Sawan", "l");
    legdNdyKstar->AddEntry(gdNdyKstarHyunji_NoUrQMD, "EPOS (No UrQMD) Hyunji", "l");
    legdNdyKstar->AddEntry(gdNdyKstarHyunji_ITY80, "EPOS (UrQMD) Hyunji (reg)", "l");
    legdNdyKstar->AddEntry(gdNdyKstarHyunji_ITY81, "EPOS (UrQMD) Hyunji (reg+res)", "l");
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
}