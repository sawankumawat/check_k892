#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

TFile *OpenFile(const string &path);
TGraphErrors *GetGraph(TFile *f, const string &name);

void plotdNdy()
{
    string KstarPath = "../../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/";
    TFile *fKstar = OpenFile(KstarPath + "Results.root");
    // TFile *fKstar_old = OpenFile(KstarPath + "Backup/Results.root");
    // TFile *fKstar_old = OpenFile(KstarPath + "Results_OlddNchdeta.root");
    TFile *fKstar_old = OpenFile(KstarPath + "Results_beforeSpectraOptimise.root");
    TFile *fpp13TeV = OpenFile("ConversionCodes/pp13TeV_KstarPhi.root");

    TGraphErrors *gMPtKstar[3], *gMYieldKstar[3], *gKstar_MeanpT_13TeV[3], *gKstar_Yield_13TeV[3], *gMPtKstarOld[3], *gMYieldKstarOld[3];

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

        gKstar_MeanpT_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_MeanpT%s", suffix.c_str()));
        gKstar_Yield_13TeV[i] = GetGraph(fpp13TeV, Form("gKstar_MeanYield%s", suffix.c_str()));

        gMPtKstarOld[i] = GetGraph(fKstar_old, Form("gMeanpTRun3%s", suffix.c_str()));
        gMYieldKstarOld[i] = GetGraph(fKstar_old, Form("gMeanYieldRun3%s", suffix.c_str()));
    }

    TGraphErrors *MYieldKstarMinBias = new TGraphErrors();
    MYieldKstarMinBias->SetPoint(0, 7.10, 0.131965);
    MYieldKstarMinBias->SetPointError(0, 0.18, 0.000139658);
    TGraphErrors *MYieldKstarMinBiasSys = new TGraphErrors();
    MYieldKstarMinBiasSys->SetPoint(0, 7.10, 0.131965);
    MYieldKstarMinBiasSys->SetPointError(0, 0.18, 0.00745802);

    TCanvas *cdNdyKstar = new TCanvas("cdNdyKstar", "cdNdyKstar", 720, 720);
    SetCanvasStyle(cdNdyKstar, 0.15, 0.03, 0.03, 0.15);
    gMYieldKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMYieldKstar[0]->GetYaxis()->SetTitle("dN/dy");
    // gMYieldKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldKstar[0]->GetXaxis()->SetLimits(0, 25.5);
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

    gMYieldKstarOld[0]->SetMarkerStyle(22);
    gMYieldKstarOld[0]->SetMarkerSize(1.5);
    gMYieldKstarOld[0]->SetMarkerColor(kBlue);
    gMYieldKstarOld[0]->SetLineColor(kBlue);
    gMYieldKstarOld[0]->SetLineWidth(3);
    gMYieldKstarOld[0]->Draw("PE same");
    gMYieldKstarOld[1]->SetLineColor(kBlue);
    gMYieldKstarOld[1]->SetFillStyle(0);
    gMYieldKstarOld[1]->SetLineWidth(3);
    gMYieldKstarOld[1]->Draw("5 same");
    gMYieldKstarOld[2]->SetFillStyle(1001);
    gMYieldKstarOld[2]->SetFillColorAlpha(kBlue, 0.3);
    gMYieldKstarOld[2]->Draw("E2 same");

    // SetGraphErrorStyle(MYieldKstarMinBias);
    // MYieldKstarMinBias->SetMarkerStyle(22);
    // MYieldKstarMinBias->SetMarkerSize(1.5);
    // MYieldKstarMinBias->SetMarkerColor(kBlue);
    // MYieldKstarMinBias->SetLineColor(kBlue);
    // MYieldKstarMinBias->SetLineWidth(3);
    // MYieldKstarMinBias->Draw("pe same");
    // SetGraphErrorStyle(MYieldKstarMinBiasSys);
    // MYieldKstarMinBiasSys->SetMarkerStyle(22);
    // MYieldKstarMinBiasSys->SetMarkerColor(kBlue);
    // MYieldKstarMinBiasSys->SetLineColor(kBlue);
    // MYieldKstarMinBiasSys->SetLineWidth(3);
    // MYieldKstarMinBiasSys->SetFillStyle(0);
    // MYieldKstarMinBiasSys->Draw("5 same");

    // gKstar_Yield_13TeV[0]->SetMarkerStyle(21);
    // gKstar_Yield_13TeV[0]->SetMarkerColor(kBlue + 1);
    // gKstar_Yield_13TeV[0]->SetLineColor(kBlue + 1);
    // gKstar_Yield_13TeV[0]->SetLineWidth(3);
    // gKstar_Yield_13TeV[0]->Draw("P same");
    // gKstar_Yield_13TeV[1]->SetLineColor(kBlue + 1);
    // gKstar_Yield_13TeV[1]->SetFillStyle(0);
    // gKstar_Yield_13TeV[1]->SetLineWidth(3);
    // gKstar_Yield_13TeV[1]->Draw("5 same");
    //  gKstar_Yield_13TeV[2]->SetFillStyle(1001);
    // gKstar_Yield_13TeV[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    // gKstar_Yield_13TeV[2]->Draw("E2 same");
    // gKstar_Yield_13TeV[2]->Fit(fpol1Run2, "R");

    TLegend *legdNdyKstar = new TLegend(0.6, 0.8, 0.9, 0.9);
    legdNdyKstar->SetTextSize(0.03);
    legdNdyKstar->SetBorderSize(0);
    legdNdyKstar->SetFillStyle(0);
    // legdNdyKstar->AddEntry(gMYieldKstar[0], "pp 13.6 TeV", "lp");
    // legdNdyKstar->AddEntry(gKstar_Yield_13TeV[0], "Run 2 (13 TeV)", "lp");
    // legdNdyKstar->AddEntry(MYieldKstarMinBias, "Min Bias (0-100%)", "lp");

    // legdNdyKstar->AddEntry(gMYieldKstar[0], "New dN_{ch}/d#eta", "lp");
    // legdNdyKstar->AddEntry(gMYieldKstarOld[0], "Old dN_{ch}/d#eta", "lp");

    legdNdyKstar->AddEntry(gMYieldKstar[0], "Optimized spectra", "lp");
    legdNdyKstar->AddEntry(gMYieldKstarOld[0], "Old spectra", "lp");

    legdNdyKstar->Draw();
    cdNdyKstar->SaveAs("Plots/dNdyKstar_compare_optmiseSpectra.png");

    TCanvas *cMeanpTKstar = new TCanvas("cMeanpTKstar", "cMeanpTKstar", 720, 720);
    SetCanvasStyle(cMeanpTKstar, 0.15, 0.03, 0.03, 0.15);
    gMPtKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMPtKstar[0]->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    gMPtKstar[0]->GetXaxis()->SetLimits(0, 25.5);
    gMPtKstar[0]->SetLineWidth(3);
    gMPtKstar[0]->GetYaxis()->SetRangeUser(0.7, 1.78);
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

    gMPtKstarOld[0]->SetMarkerStyle(22);
    gMPtKstarOld[0]->SetMarkerSize(1.5);
    gMPtKstarOld[0]->SetMarkerColor(kBlue);
    gMPtKstarOld[0]->SetLineColor(kBlue);
    gMPtKstarOld[0]->SetLineWidth(3);
    gMPtKstarOld[0]->Draw("PE same");
    gMPtKstarOld[1]->SetLineColor(kBlue);
    gMPtKstarOld[1]->SetFillStyle(0);
    gMPtKstarOld[1]->SetLineWidth(3);
    gMPtKstarOld[1]->Draw("5 same");
    gMPtKstarOld[2]->SetFillStyle(1001);
    gMPtKstarOld[2]->SetFillColorAlpha(kBlue, 0.3);
    gMPtKstarOld[2]->Draw("E2 same");

    legdNdyKstar->Draw();
    cMeanpTKstar->SaveAs("Plots/MeanpTKstar_compare_optmiseSpectra.png");
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