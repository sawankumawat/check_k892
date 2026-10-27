#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

// Helper function to format and draw Data vs MC comparison in a canvas pad
void DrawDataMCComparison(TH1F *hData, TH1F *hMC, TString xTitle, double lowpT, double highpT, bool useNsigmaRange = true)
{
    if (!hData || !hMC)
        return;

    SetHistoQA(hData);
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetYaxis()->SetTitle("Normalized counts");

    if (useNsigmaRange)
    {
        hData->GetXaxis()->SetRangeUser(-2.5, 2.5);
    }
    else
    {
        hData->GetXaxis()->SetRangeUser(-0.1, 0.1);
    }

    // Normalization range calculation
    int normBinLow = useNsigmaRange ? hData->FindBin(-1.0 + 0.001) : 1;
    int normBinHigh = useNsigmaRange ? hData->FindBin(1.0 - 0.001) : hData->GetNbinsX();

    double intData = hData->Integral(normBinLow, normBinHigh);
    double intMC = hMC->Integral(normBinLow, normBinHigh);

    if (intData > 0)
        hData->Scale(1.0 / intData);
    if (intMC > 0)
        hMC->Scale(1.0 / intMC);

    hData->SetMaximum(hData->GetMaximum() * 1.3);
    hData->Draw("pe");

    SetHistoQA(hMC);
    hMC->SetMarkerColor(kRed);
    hMC->SetLineColor(kRed);
    hMC->Draw("pe same");

    TLegend *legend = new TLegend(0.68, 0.72, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
    legend->AddEntry(hData, "Data", "p");
    legend->AddEntry(hMC, "MC", "p");
    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.20, 0.88, "ALICE");
    latex.DrawLatex(0.20, 0.83, "pp, #sqrt{#it{s}} = 13.6 TeV");
    latex.DrawLatex(0.20, 0.78, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", lowpT, highpT));
    latex.DrawLatex(0.20, 0.73, "FT0M: 0-100%");
}

void compare_PID_DataMC()
{
    gStyle->SetOptStat(0);
    TString dataPath = "/home/sawan/Storage/check_k892/data/kstar/LHC22o_pass7";
    TString mcPath = "/home/sawan/Storage/check_k892/mc/LHC24f3c";

    TFile *fData = new TFile(dataPath + "/751768.root", "read");
    TFile *fMC = new TFile(mcPath + "/CalibratedMCchecks/780946.root", "read");
    // TFile *fMC = new TFile(mcPath + "/750013.root", "read");
    // TFile *fMC = new TFile(mcPath + "/734009.root", "read");

    // TFile *fData = new TFile(dataPath + "/679906.root", "read");
    // TFile *fMC = new TFile(mcPath + "/679945.root", "read");

    // TFile *fData = new TFile(dataPath + "/679906.root", "read"); // LHC24f4d
    // TFile *fMC = new TFile(mcPath + "/762508.root", "read"); // LHC24f4e

    if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie())
    {
        cerr << "Error opening the root files!" << endl;
        return;
    }

    // string multpathMC = "kstarqa_WithoutTOFShift";
    string multpathMC = "kstarqa_ServiceWagon";
    string multpathData = "kstarqa";

    // Fetch DCA histograms (Data & MC)
    TH3F *hDCAxy3D = (TH3F *)fData->Get(Form("%s/eventSelection/hDcaxy_cent_pt", multpathData.c_str()));
    TH3F *hDCAz3D = (TH3F *)fData->Get(Form("%s/eventSelection/hDcaz_cent_pt", multpathData.c_str()));
    TH3F *hDCAxy3DMC = (TH3F *)fMC->Get(Form("%s/eventSelection/hDcaxy_cent_pt", multpathMC.c_str()));
    TH3F *hDCAz3DMC = (TH3F *)fMC->Get(Form("%s/eventSelection/hDcaz_cent_pt", multpathMC.c_str()));

    // Fetch PID histograms (Data)
    TH3F *hNsigmaTPCKaon = (TH3F *)fData->Get(Form("%s/hPID/Before/hTPCnsigKa_mult_pt", multpathData.c_str()));
    TH3F *hNsigmaTPCPion = (TH3F *)fData->Get(Form("%s/hPID/Before/hTPCnsigPi_mult_pt", multpathData.c_str()));
    TH3F *hNsigmaTOFKaon = (TH3F *)fData->Get(Form("%s/hPID/Before/hTOFnsigKa_mult_pt", multpathData.c_str()));
    TH3F *hNsigmaTOFPion = (TH3F *)fData->Get(Form("%s/hPID/Before/hTOFnsigPi_mult_pt", multpathData.c_str()));

    // Fetch PID histograms (MC)
    TH3F *hNsigmaTPCKaonMC = (TH3F *)fMC->Get(Form("%s/hPID/Before/hTPCnsigKa_mult_pt", multpathMC.c_str()));
    TH3F *hNsigmaTPCPionMC = (TH3F *)fMC->Get(Form("%s/hPID/Before/hTPCnsigPi_mult_pt", multpathMC.c_str()));
    TH3F *hNsigmaTOFKaonMC = (TH3F *)fMC->Get(Form("%s/hPID/Before/hTOFnsigKa_mult_pt", multpathMC.c_str()));
    TH3F *hNsigmaTOFPionMC = (TH3F *)fMC->Get(Form("%s/hPID/Before/hTOFnsigPi_mult_pt", multpathMC.c_str()));

    if (!hNsigmaTPCKaon || !hNsigmaTPCPion || !hNsigmaTOFKaon || !hNsigmaTOFPion ||
        !hNsigmaTPCKaonMC || !hNsigmaTPCPionMC || !hNsigmaTOFKaonMC || !hNsigmaTOFPionMC)
    {
        cerr << "PID histograms Before selection not found!" << endl;
        return;
    }

    // Binning definitions
    double pTrangesForTPCandTOF[] = {0.0, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0};
    const int nPtBins = 6;

    double multLow = 70.0;
    double multHigh = 100.0;
    int lowBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multLow + 0.001);
    int highBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multHigh - 0.001);

    // Create Divided Canvases (3 x 2)
    TCanvas *cDCAxy_ptBin = new TCanvas("cDCAxy_ptBin", "DCAxy vs pT Bins", 1500, 1000);
    TCanvas *cDCAz_ptBin = new TCanvas("cDCAz_ptBin", "DCAz vs pT Bins", 1500, 1000);
    TCanvas *cTPCPion_ptBin = new TCanvas("cTPCPion_ptBin", "TPC nSigma Pion vs pT Bins", 1500, 1000);
    TCanvas *cTPCKaon_ptBin = new TCanvas("cTPCKaon_ptBin", "TPC nSigma Kaon vs pT Bins", 1500, 1000);
    TCanvas *cTOFPion_ptBin = new TCanvas("cTOFPion_ptBin", "TOF nSigma Pion vs pT Bins", 1500, 1000);
    TCanvas *cTOFKaon_ptBin = new TCanvas("cTOFKaon_ptBin", "TOF nSigma Kaon vs pT Bins", 1500, 1000);

    cDCAxy_ptBin->Divide(3, 2);
    cDCAz_ptBin->Divide(3, 2);
    cTPCPion_ptBin->Divide(3, 2);
    cTPCKaon_ptBin->Divide(3, 2);
    cTOFPion_ptBin->Divide(3, 2);
    cTOFKaon_ptBin->Divide(3, 2);

    // Loop over 6 pT bins
    for (int i = 0; i < nPtBins; ++i)
    {
        double lowpT = pTrangesForTPCandTOF[i];
        double highpT = pTrangesForTPCandTOF[i + 1];

        int lowBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(lowpT + 0.001);
        int highBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(highpT - 0.001);

        // 1D Projections (Data)
        TH1F *h1DDCAxy = hDCAxy3D ? (TH1F *)hDCAxy3D->ProjectionX(Form("h1D_DCAxy_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT) : nullptr;
        TH1F *h1DDCAz = hDCAz3D ? (TH1F *)hDCAz3D->ProjectionX(Form("h1D_DCAz_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT) : nullptr;

        TH1F *h1DTPCPion = (TH1F *)hNsigmaTPCPion->ProjectionX(Form("h1D_TPC_Pi_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        TH1F *h1DTPCKaon = (TH1F *)hNsigmaTPCKaon->ProjectionX(Form("h1D_TPC_Ka_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        TH1F *h1DTOFPion = (TH1F *)hNsigmaTOFPion->ProjectionX(Form("h1D_TOF_Pi_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        TH1F *h1DTOFKaon = (TH1F *)hNsigmaTOFKaon->ProjectionX(Form("h1D_TOF_Ka_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);

        // 1D Projections (MC)
        TH1F *h1DDCAxyMC = hDCAxy3DMC ? (TH1F *)hDCAxy3DMC->ProjectionX(Form("h1D_DCAxy_MC_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT) : nullptr;
        TH1F *h1DDCAzMC = hDCAz3DMC ? (TH1F *)hDCAz3DMC->ProjectionX(Form("h1D_DCAz_MC_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT) : nullptr;

        TH1F *h1DTPCPionMC = (TH1F *)hNsigmaTPCPionMC->ProjectionX(Form("h1D_TPC_Pi_MC_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        TH1F *h1DTPCKaonMC = (TH1F *)hNsigmaTPCKaonMC->ProjectionX(Form("h1D_TPC_Ka_MC_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        TH1F *h1DTOFPionMC = (TH1F *)hNsigmaTOFPionMC->ProjectionX(Form("h1D_TOF_Pi_MC_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        TH1F *h1DTOFKaonMC = (TH1F *)hNsigmaTOFKaonMC->ProjectionX(Form("h1D_TOF_Ka_MC_bin_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);

        // 1. DCA xy Canvas
        cDCAxy_ptBin->cd(i + 1);
        DrawDataMCComparison(h1DDCAxy, h1DDCAxyMC, "DCA_{xy} (cm)", lowpT, highpT, false);

        // 2. DCA z Canvas
        cDCAz_ptBin->cd(i + 1);
        DrawDataMCComparison(h1DDCAz, h1DDCAzMC, "DCA_{z} (cm)", lowpT, highpT, false);

        // 3. TPC Pion Canvas
        cTPCPion_ptBin->cd(i + 1);
        DrawDataMCComparison(h1DTPCPion, h1DTPCPionMC, "#it{n#sigma}_{TPC}^{#pi}", lowpT, highpT, true);

        // 4. TPC Kaon Canvas
        cTPCKaon_ptBin->cd(i + 1);
        DrawDataMCComparison(h1DTPCKaon, h1DTPCKaonMC, "#it{n#sigma}_{TPC}^{K}", lowpT, highpT, true);

        // 5. TOF Pion Canvas
        cTOFPion_ptBin->cd(i + 1);
        gPad->SetGridx();
        DrawDataMCComparison(h1DTOFPion, h1DTOFPionMC, "#it{n#sigma}_{TOF}^{#pi}", lowpT, highpT, true);

        // 6. TOF Kaon Canvas
        cTOFKaon_ptBin->cd(i + 1);
        gPad->SetGridx();
        DrawDataMCComparison(h1DTOFKaon, h1DTOFKaonMC, "#it{n#sigma}_{TOF}^{K}", lowpT, highpT, true);
    }
}