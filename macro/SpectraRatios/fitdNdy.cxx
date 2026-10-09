#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

TFile *OpenFile(const string &path);
TGraphErrors *GetGraph(TFile *f, const string &name);

void fitdNdy()
{
    string KstarPath = "../../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/";
    TFile *fKstar = OpenFile(KstarPath + "Results.root");
    TFile *fpp13TeV = OpenFile("ConversionCodes/pp13TeV_KstarPhi.root");

    TGraphErrors *gMPtKstar[3], *gMYieldKstar[3], *gKstar_MeanpT_13TeV[3], *gKstar_Yield_13TeV[3];

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
    }

    TF1 *fpol1Run2 = new TF1("fpol1Run2", "[0] + [1]*x", 0, 25);
    fpol1Run2->SetLineColor(kBlue + 1);
    fpol1Run2->SetLineWidth(2);
    fpol1Run2->SetLineStyle(2);
    TF1 *fpol1Run3 = new TF1("fpol1Run3", "[0] + [1]*x", 0, 25);
    fpol1Run3->SetLineColor(kRed + 1);
    fpol1Run3->SetLineWidth(2);

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
    // gMYieldKstar[1]->SetFillStyle(0);
    // gMYieldKstar[1]->SetLineColor(kRed + 1);
    // gMYieldKstar[1]->SetLineWidth(3);
    // gMYieldKstar[1]->Draw("5 same");
    gMYieldKstar[2]->SetFillStyle(1001);
    gMYieldKstar[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gMYieldKstar[2]->Draw("E2 same");
    gMYieldKstar[2]->Fit(fpol1Run3, "R");

    gKstar_Yield_13TeV[0]->SetMarkerStyle(21);
    gKstar_Yield_13TeV[0]->SetMarkerColor(kBlue + 1);
    gKstar_Yield_13TeV[0]->SetLineColor(kBlue + 1);
    gKstar_Yield_13TeV[0]->SetLineWidth(3);
    gKstar_Yield_13TeV[0]->Draw("P same");
    gKstar_Yield_13TeV[1]->SetLineColor(kBlue + 1);
    // gKstar_Yield_13TeV[1]->SetFillStyle(0);
    // gKstar_Yield_13TeV[1]->SetLineWidth(3);
    // gKstar_Yield_13TeV[1]->Draw("5 same");
     gKstar_Yield_13TeV[2]->SetFillStyle(1001);
    gKstar_Yield_13TeV[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gKstar_Yield_13TeV[2]->Draw("E2 same");
    gKstar_Yield_13TeV[2]->Fit(fpol1Run2, "R");

    double slopeRun3 = fpol1Run3->GetParameter(1);
    double slopeRun3Err = fpol1Run3->GetParError(1);
    double interceptRun3 = fpol1Run3->GetParameter(0);
    double interceptRun3Err = fpol1Run3->GetParError(0);

    double slopeRun2 = fpol1Run2->GetParameter(1);
    double slopeRun2Err = fpol1Run2->GetParError(1);
    double interceptRun2 = fpol1Run2->GetParameter(0);
    double interceptRun2Err = fpol1Run2->GetParError(0);

    TLegend *legdNdyKstar = new TLegend(0.6, 0.8, 0.9, 0.9);
    legdNdyKstar->SetTextSize(0.03);
    legdNdyKstar->SetBorderSize(0);
    // legdNdyKstar->SetHeader("Uncertainty: Total systematic");
    // legdNdyKstar->SetHeader("Uncertainty: Statistical");
    legdNdyKstar->SetHeader("Uncertainty: Uncorrelated");
    legdNdyKstar->AddEntry(gMYieldKstar[0], "Run 3 (13.6 TeV)", "lp");
    legdNdyKstar->AddEntry(gKstar_Yield_13TeV[0], "Run 2 (13 TeV)", "lp");
    legdNdyKstar->Draw();

    TLatex *latdNdyKstar = new TLatex();
    latdNdyKstar->SetTextSize(0.03);
    latdNdyKstar->SetTextFont(22);
    latdNdyKstar->DrawLatexNDC(0.19, 0.9, "Fit (pp 13.6 TeV):");
    latdNdyKstar->SetTextFont(42);
    latdNdyKstar->DrawLatexNDC(0.19, 0.85, Form("Slope = %.4f #pm %.4f", slopeRun3, slopeRun3Err));
    latdNdyKstar->DrawLatexNDC(0.19, 0.8, Form("Intercept = %.4f #pm %.4f", interceptRun3, interceptRun3Err));
    latdNdyKstar->SetTextFont(22);
    latdNdyKstar->DrawLatexNDC(0.19, 0.75, "Fit (pp 13 TeV):");
    latdNdyKstar->SetTextFont(42);
    latdNdyKstar->DrawLatexNDC(0.19, 0.7, Form("Slope = %.4f #pm %.4f", slopeRun2, slopeRun2Err));
    latdNdyKstar->DrawLatexNDC(0.19, 0.65, Form("Intercept = %.4f #pm %.4f", interceptRun2, interceptRun2Err));
    cdNdyKstar->SaveAs("dNdyKstarFitUncorr.png");
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