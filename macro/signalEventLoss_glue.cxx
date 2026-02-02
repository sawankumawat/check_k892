#include <iostream>
using namespace std;
#include "src/style.h"

void signalEventLoss_glue()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Input file
    TFile *fInputFile = new TFile("/home/sawan/Downloads/AnalysisResults.root", "read");
    if (fInputFile->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string histpathf01710 = "higher-mass-resonances/hMChists";
    string histpathf21525 = "higher-mass-resonances_f21525/hMChists";

    double ptbins[] = {1, 2, 3, 5, 7, 10, 15};
    int sizePtBins = sizeof(ptbins) / sizeof(ptbins[0]) - 1;
    TH1F *hEventLossf0 = new TH1F("hEventLossf0", "Event Loss vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); Event Loss", sizePtBins, ptbins);
    TH1F *hEventLossf2 = new TH1F("hEventLossf2", "Event Loss vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); Event Loss", sizePtBins, ptbins);
    TH1F *hSignalLossf0 = new TH1F("hSignalLossf0", "Signal Loss vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); Signal Loss", sizePtBins, ptbins);
    TH1F *hSignalLossf2 = new TH1F("hSignalLossf2", "Signal Loss vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); Signal Loss", sizePtBins, ptbins);
    TH1F *hEventbySignalLoss = new TH1F("hEventbySignalLoss", "Event Loss / Signal Loss vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); Event Loss / Signal Loss", sizePtBins, ptbins);

    // //Event loss histograms
    // Method 1 (from eventsignalloss function)
    TH1F *hEventLossNumf0 = (TH1F *)fInputFile->Get(Form("%s/MCcorrections/MultiplicityRec", histpathf01710.c_str()));
    TH1F *hEventLossDenf0 = (TH1F *)fInputFile->Get(Form("%s/MCcorrections/MultiplicityGen", histpathf01710.c_str()));

    // // Method 2 (from generated process function)
    // TH1F *hEventLossNumf0 = (TH1F *)fInputFile->Get(Form("%s/MCcorrections/hGenNo", histpathf01710.c_str()));

    // // Signal loss histograms
    // Method 1 (from signalsignalloss function)
    TH2F *h2DSignalLossNumf0 = (TH2F *)fInputFile->Get(Form("%s/MCcorrections/hSignalLossNumerator", histpathf01710.c_str()));
    TH2F *h2DSignalLossDenf0 = (TH2F *)fInputFile->Get(Form("%s/MCcorrections/hSignalLossDenominator", histpathf01710.c_str()));
    int multBinlow = h2DSignalLossDenf0->GetYaxis()->FindBin(0 + 0.001);
    int multBinhigh = h2DSignalLossDenf0->GetYaxis()->FindBin(100 - 0.001);

    TH1F *hSignalLossNumf0 = (TH1F *)h2DSignalLossNumf0->ProjectionX("hSignalLossNumf0", multBinlow, multBinhigh);
    TH1F *hSignalLossDenf0 = (TH1F *)h2DSignalLossDenf0->ProjectionX("hSignalLossDenf0", multBinlow, multBinhigh);

    // // Method 2 (from generated process function)
    // THnSparseF *hSparseSignalLossNumf0 = (THnSparseF *)fInputFile->Get(Form("%s/Genf1710", histpathf01710.c_str()));
    // TH2F *h2DSignalLossDenf0 = (TH2F *)fInputFile->Get(Form("%s/MCcorrections/hSignalLossDenominator4", histpathf01710.c_str()));
    // int multBinlow = h2DSignalLossDenf0->GetYaxis()->FindBin(0 + 0.001);
    // int multBinhigh = h2DSignalLossDenf0->GetYaxis()->FindBin(100 - 0.001);

    // hSparseSignalLossNumf0->GetAxis(0)->SetRange(multBinlow, multBinhigh);  // setting multiplicity range
    // TH1F *hSignalLossNumf0 = (TH1F *)hSparseSignalLossNumf0->Projection(1); // projecting pt axis
    // TH1F *hSignalLossDenf0 = (TH1F *)h2DSignalLossDenf0->ProjectionX("hSignalLossDenf0", multBinlow, multBinhigh);

    for (int ipt = 0; ipt < sizePtBins; ipt++)
    {
        // // Event loss (method 1)
        double eventLossNum = hEventLossNumf0->Integral();
        double eventLossDen = hEventLossDenf0->Integral();
        double eventLoss = 0.0;
        if (eventLossDen > 0)
            eventLoss = eventLossNum / eventLossDen;
        hEventLossf0->SetBinContent(ipt + 1, eventLoss);

        // // // Event loss (method 2)
        // double eventLossDen = hEventLossNumf0->GetBinContent(1);
        // double eventLossNum = hEventLossNumf0->GetBinContent(2);
        // double eventLoss = 0.0;
        // if (eventLossDen > 0)
        //     eventLoss = eventLossNum / eventLossDen;
        // hEventLossf0->SetBinContent(ipt + 1, eventLoss);

        if (ipt == 0)
            cout << "Event loss is " << eventLoss << endl;

        // Signal loss (method 1)
        int lowpTbin = hSignalLossNumf0->GetXaxis()->FindBin(ptbins[ipt] + 0.01);
        int highpTbin = hSignalLossNumf0->GetXaxis()->FindBin(ptbins[ipt + 1] - 0.01);
        double signalLossNum = hSignalLossNumf0->Integral(lowpTbin, highpTbin);
        double signalLossDen = hSignalLossDenf0->Integral(lowpTbin, highpTbin);
        double signalLoss = 0.0;
        if (signalLossDen > 0)
            signalLoss = signalLossNum / signalLossDen;
        hSignalLossf0->SetBinContent(ipt + 1, signalLoss);

        cout << "pT bin " << ipt << " signal loss " << signalLoss << endl;

        // Event by Signal loss
        double eventbySignalLoss = 0.0;
        if (signalLoss > 0)
            eventbySignalLoss = eventLoss / signalLoss;
        hEventbySignalLoss->SetBinContent(ipt + 1, eventbySignalLoss);
    }

    TCanvas *cLoss = new TCanvas("cLoss", "Event Loss vs pT", 720, 720);
    SetCanvasStyle(cLoss, 0.14, 0.03, 0.05, 0.14);
    SetHistoQA(hEventLossf0);
    hEventLossf0->GetYaxis()->SetRangeUser(0, 1);
    hEventLossf0->Draw("HIST");
    hSignalLossf0->SetMarkerStyle(20);
    hSignalLossf0->SetMarkerColor(kRed);
    hSignalLossf0->SetLineColor(kRed);
    hSignalLossf0->Draw("p same");

    TCanvas *cRatioLoss = new TCanvas("cRatioLoss", "Event Loss / Signal Loss vs pT", 720, 720);
    SetCanvasStyle(cRatioLoss, 0.14, 0.03, 0.05, 0.14);
    SetHistoQA(hEventbySignalLoss);
    hEventbySignalLoss->GetYaxis()->SetRangeUser(0.6, 1.4);
    hEventbySignalLoss->Draw("p");
}