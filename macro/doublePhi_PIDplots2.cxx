#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "src/style.h"

void doublePhi_PIDplots2()
{
    gStyle->SetOptStat(0);
    TString subWagon = "";

    TFile *file = new TFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults.root");
    if (file->IsZombie() || !file)
    {
        std::cerr << "File not found " << std::endl;
        return;
    }

    TH2D *hnSigmaTPCKaonPlus = (TH2D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTPCKaonPlusBefore");
    TH2D *hnSigmaTPCKaonMinus = (TH2D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTPCKaonMinusBefore");
    TH3D *hnSigmaTPCTOFKaonBefore = (TH3D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTPCTOFKaonBefore");
    if (!hnSigmaTPCKaonPlus || !hnSigmaTPCKaonMinus || !hnSigmaTPCTOFKaonBefore)
    {
        std::cerr << "Error: Could not find required histograms in file\n";
        return;
    }

    TCanvas *cPIDTPCKaPos = new TCanvas("cPIDTPCKaPos", "PIDTPCKa Positive", 1080, 720);
    SetCanvasStyle(cPIDTPCKaPos, 0.1, 0.05, 0.06, 0.17);
    cPIDTPCKaPos->cd();
    TPad *mainPadPos = new TPad("mainPadPos", "", 0.0, 0.0, 1.0, 0.96);
    mainPadPos->SetMargin(0.0, 0.0, 0.0, 0.0);
    mainPadPos->Draw();
    mainPadPos->cd();
    mainPadPos->Divide(4, 4, 0.001, 0.001);

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.09);

    float pT_bins[] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 5.0, 10.0};
    int total_pT_bins = sizeof(pT_bins) / sizeof(pT_bins[0]) - 1;
    float pion_contamination_peak_mean[] = {-7, -7, -6, -5, -3, -2, -1, 1, 1, 2, 3, 3, 3, 3, 3, 3};

    // Canvas to display FOM vs nSigma cut for each pT bin
    TCanvas *cFOM = new TCanvas("cFOM", "Figure of Merit vs nSigma Cut", 1080, 720);
    SetCanvasStyle(cFOM, 0.1, 0.05, 0.06, 0.17);
    cFOM->cd(); // <--- FIX: Switch focus to cFOM instead of cPIDTPCKaPos

    TPad *mainPadFOM = new TPad("mainPadFOM", "", 0.0, 0.0, 1.0, 0.96);
    mainPadFOM->SetMargin(0.0, 0.0, 0.0, 0.0);
    mainPadFOM->Draw();
    mainPadFOM->cd();
    mainPadFOM->Divide(4, 4, 0.001, 0.001);

    std::cout << "\n========================================================" << std::endl;
    std::cout << " OPTIMAL TPC nSigma CUTS (MAXIMIZING FOM = Purity x Eff) " << std::endl;
    std::cout << "========================================================" << std::endl;

    for (Int_t ip = 0; ip < total_pT_bins; ip++)
    {
        int binLow = hnSigmaTPCKaonPlus->GetYaxis()->FindBin(pT_bins[ip] + 0.001);
        int binHigh = hnSigmaTPCKaonPlus->GetYaxis()->FindBin(pT_bins[ip + 1] - 0.001);
        TH1D *hKaonTPC1Dpos = hnSigmaTPCKaonPlus->ProjectionX(Form("hKaonTPC1Dpos_%.1f_%.1f", pT_bins[ip], pT_bins[ip + 1]), binLow, binHigh);
        SetHistoQA(hKaonTPC1Dpos);

        mainPadPos->cd(ip + 1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.02);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.07);
        hKaonTPC1Dpos->GetXaxis()->SetTitle("TPC N_{#sigma} Kaon");
        hKaonTPC1Dpos->GetYaxis()->SetTitle("Counts");
        hKaonTPC1Dpos->GetXaxis()->SetRangeUser(-5.5, 5.5);
        hKaonTPC1Dpos->SetMaximum(hKaonTPC1Dpos->GetMaximum() * 7);
        hKaonTPC1Dpos->SetLineColor(kBlue + 1);
        hKaonTPC1Dpos->SetLineWidth(2);
        hKaonTPC1Dpos->Draw("HIST");

        lat.DrawLatexNDC(0.2, 0.83, Form("p_{T} = %.2f - %.2f GeV/c", pT_bins[ip], pT_bins[ip + 1]));

        // Direct Double Gaussian Fitting Setup
        TF1 *gausFitKaon = new TF1("gausFitKaon", "gaus", -3, 3);
        gausFitKaon->SetParameter(0, hKaonTPC1Dpos->GetMaximum());
        gausFitKaon->FixParameter(1, 0);
        gausFitKaon->SetParameter(2, 1);
        gausFitKaon->SetParLimits(2, 0.8, 1.8);
        hKaonTPC1Dpos->Fit(gausFitKaon, "REBMS0Q");

        TF1 *gausFitPion = new TF1("gausFitPion", "gaus", -5, 5);
        gausFitPion->SetParameter(0, hKaonTPC1Dpos->GetMaximum());
        gausFitPion->SetParameter(1, pion_contamination_peak_mean[ip]);
        gausFitPion->SetParLimits(1, pion_contamination_peak_mean[ip] - 0.6, pion_contamination_peak_mean[ip] + 0.6);
        gausFitPion->SetParameter(2, 1);
        gausFitPion->SetParLimits(2, 0.5, 1.5);
        hKaonTPC1Dpos->Fit(gausFitPion, "REBMS0Q");

        TF1 *doubleGausFit = new TF1("doubleGausFit", "gaus(0)+gaus(3)", -5, 5);
        doubleGausFit->SetParameter(0, gausFitKaon->GetParameter(0));
        doubleGausFit->FixParameter(1, gausFitKaon->GetParameter(1));
        doubleGausFit->SetParLimits(1, gausFitKaon->GetParameter(1) - 1, gausFitKaon->GetParameter(1) + 1);
        doubleGausFit->SetParameter(2, gausFitKaon->GetParameter(2));
        doubleGausFit->SetParameter(3, gausFitPion->GetParameter(0));
        doubleGausFit->SetParameter(4, gausFitPion->GetParameter(1));
        doubleGausFit->SetParLimits(4, gausFitPion->GetParameter(1) - 0.8, gausFitPion->GetParameter(1) + 0.8);
        doubleGausFit->SetParameter(5, gausFitPion->GetParameter(2));
        doubleGausFit->SetLineColor(kRed);
        hKaonTPC1Dpos->Fit(doubleGausFit, "REBMSQ");
        doubleGausFit->Draw("same");

        TF1 *fitFucnGausTempKaon = new TF1("fitFucnGausTempKaon", "gaus", -6, 6);
        fitFucnGausTempKaon->SetParameters(doubleGausFit->GetParameter(0), doubleGausFit->GetParameter(1), doubleGausFit->GetParameter(2));
        fitFucnGausTempKaon->SetLineColor(kBlue);
        fitFucnGausTempKaon->SetLineStyle(2);

        TF1 *fitFucnGausTempPion = new TF1("fitFucnGausTempPion", "gaus", -6, 6);
        fitFucnGausTempPion->SetParameters(doubleGausFit->GetParameter(3), doubleGausFit->GetParameter(4), doubleGausFit->GetParameter(5));
        fitFucnGausTempPion->SetLineColor(kGreen + 2);
        fitFucnGausTempPion->SetLineStyle(2);

        fitFucnGausTempPion->Draw("same");
        fitFucnGausTempKaon->Draw("same");

        // --- Figure of Merit Calculation ---
        double signalTotal = fitFucnGausTempKaon->Integral(-5.0, 5.0);
        std::vector<double> nSigmaCutVal, fomVal;
        double maxFOM = -1.0;
        double bestCutMin = -2.0, bestCutMax = 2.0;

        // Scan asymmetric cut ranges around Kaon Mean
        double kMean = doubleGausFit->GetParameter(1);
        for (double cutWindow = 0.5; cutWindow <= 3.5; cutWindow += 0.05)
        {
            double lowCut = kMean - cutWindow;
            double highCut = kMean + cutWindow;

            double S_cut = fitFucnGausTempKaon->Integral(lowCut, highCut);
            double B_cut = fitFucnGausTempPion->Integral(lowCut, highCut);

            if (S_cut + B_cut <= 0 || signalTotal <= 0)
                continue;

            double purity = S_cut / (S_cut + B_cut);
            double efficiency = S_cut / signalTotal;
            double fom = purity * efficiency;

            nSigmaCutVal.push_back(cutWindow);
            fomVal.push_back(fom);

            if (fom > maxFOM)
            {
                maxFOM = fom;
                bestCutMin = lowCut;
                bestCutMax = highCut;
            }
        }

        // Output optimal cut for current pT range
        std::cout << Form("pT [%.2f - %.2f GeV/c] -> Best Cut: [%.2f, %.2f] nSigma | Max FOM: %.4f",
                          pT_bins[ip], pT_bins[ip + 1], bestCutMin, bestCutMax, maxFOM)
                  << std::endl;

        // Draw FOM Graph
        mainPadFOM->cd(ip + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.02);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.07);

        TGraph *gFOM = new TGraph(nSigmaCutVal.size(), &nSigmaCutVal[0], &fomVal[0]);
        SetGraphStyle(gFOM);
        gFOM->SetTitle(0);
        gFOM->GetXaxis()->SetTitle("|#Delta N_{#sigma}| Window");
        gFOM->GetYaxis()->SetTitle("FOM (Purity #times Efficiency)");
        gFOM->SetMarkerStyle(20);
        gFOM->SetMarkerSize(0.6);
        gFOM->SetMarkerColor(kRed + 1);
        gFOM->SetLineColor(kRed + 1);
        gFOM->Draw("APL");

        lat.DrawLatexNDC(0.2, 0.83, Form("p_{T} = %.2f - %.2f", pT_bins[ip], pT_bins[ip + 1]));
        lat.DrawLatexNDC(0.2, 0.70, Form("Best Cut: [%.2f, %.2f]", bestCutMin, bestCutMax));
    }

    cPIDTPCKaPos->Modified();
    cPIDTPCKaPos->Update();

    cFOM->Modified();
    cFOM->Update();
}