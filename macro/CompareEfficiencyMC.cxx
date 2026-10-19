
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TColor.h"

using namespace std;

void CompareEfficiencyMC()
{
    // ================= MC FILES =================
    vector<string> files = {
        "/home/sawan/Storage/check_k892/mc/LHC24f3c/750013.root",
        // "/home/sawan/Storage/check_k892/mc/LHC24f3c/CalibratedMCchecks/780945.root",
        "/home/sawan/Storage/check_k892/mc/LHC24f3c/CalibratedMCchecks/780946.root",
        "/home/sawan/Storage/check_k892/mc/LHC24f3c/CalibratedMCchecks/780947.root",
        "/home/sawan/Storage/check_k892/mc/LHC24f3c/CalibratedMCchecks/780948.root"};

    // vector<string> labels = {"MC 1", "MC 2", "MC 3", "MC 4", "MC 5"};
    vector<string> labels = {"MC 1", "MC 3", "MC 4", "MC 5"};

    const int nFiles = files.size();

    int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 2, kCyan};
    int markers[] = {20, 21, 22, 23, 24};

    // ================= pT BINNING =================
    const int Npt = 23;
    double pT_bins[Npt + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};

    // ================= HISTOGRAM PATHS =================
    // Histogram directory for each corresponding file
    vector<TString> genPaths = {
        "kstarqa/hInvMass",
        // "kstarqa_Ercoles/hInvMass",
        "kstarqa_ServiceWagon/hInvMass",
        "kstarqa_serviceWagon_recalibMode0/hInvMass",
        "kstarqa_serviceWagon_recalibMode1/hInvMass"};

    gStyle->SetOptStat(0);

    // ================= OUTPUT HISTOGRAMS =================
    vector<TH1D *> hEfficiency(nFiles, nullptr);
    vector<TH1D *> hEventLoss(nFiles, nullptr);
    vector<TH1D *> hSignalLoss(nFiles, nullptr);
    vector<TH1D *> hEventLossOverSignalLoss(nFiles, nullptr);

    // ======================================================
    // PROCESS EACH MC FILE
    // ======================================================
    for (int i = 0; i < nFiles; ++i)
    {

        cout << "\n====================================\n";
        cout << "Processing: " << files[i] << endl;
        cout << "====================================\n";

        TFile *f = TFile::Open(files[i].c_str(), "READ");

        if (!f || f->IsZombie())
        {
            cerr << "ERROR: Cannot open " << files[i] << endl;
            if (f)
            {
                f->Close();
                delete f;
            }
            continue;
        }

        TString genPath = genPaths[i];

        TString genPtPath = genPath + "/hk892GenpTCalib1";
        TString recPtPath = genPath + "/h2KstarRecptCalib1";

        TString eventDenPath = genPath + "/MCcorrections/MultiplicityGen";
        TString eventNumPath = genPath + "/MCcorrections/MultiplicityRec";

        TString signalDenPath = genPath + "/MCcorrections/hSignalLossDenominator";
        TString signalNumPath = genPath + "/MCcorrections/hSignalLossNumerator";

        // ================= READ HISTOGRAMS =================
        auto *hGenPt = dynamic_cast<THnSparse *>(f->Get(genPtPath));
        auto *hRecPt = dynamic_cast<THnSparse *>(f->Get(recPtPath));

        auto *hEventDen = dynamic_cast<TH1 *>(f->Get(eventDenPath));
        auto *hEventNum = dynamic_cast<TH1 *>(f->Get(eventNumPath));

        auto *hSignalDen = dynamic_cast<TH2 *>(f->Get(signalDenPath));
        auto *hSignalNum = dynamic_cast<TH2 *>(f->Get(signalNumPath));

        if (!hGenPt || !hRecPt || !hEventDen || !hEventNum ||
            !hSignalDen || !hSignalNum)
        {

            cerr << "ERROR: One or more required histograms are missing "
                 << "in " << files[i] << endl;

            f->Close();
            delete f;
            continue;
        }

        // ==================================================
        // MINIMUM BIAS: 0-100%
        // ==================================================
        const double multLow = 0.0;
        const double multHigh = 100.0;

        // ================= EFFICIENCY =====================
        // Axis 0: pT
        // Axis 1: multiplicity

        int genMultLow = hGenPt->GetAxis(1)->FindBin(multLow + 0.001);
        int genMultHigh = hGenPt->GetAxis(1)->FindBin(multHigh - 0.001);

        int recMultLow = hRecPt->GetAxis(1)->FindBin(multLow + 0.001);
        int recMultHigh = hRecPt->GetAxis(1)->FindBin(multHigh - 0.001);

        hGenPt->GetAxis(1)->SetRange(genMultLow, genMultHigh);
        hRecPt->GetAxis(1)->SetRange(recMultLow, recMultHigh);

        TH1D *hGenProjection = hGenPt->Projection(0, "E");
        TH1D *hRecProjection = hRecPt->Projection(0, "E");

        hEfficiency[i] = new TH1D(Form("hEfficiency_%d", i), ";#it{p}_{T} (GeV/#it{c});Acceptance #times efficiency", Npt, pT_bins);
        hEfficiency[i]->SetDirectory(nullptr);
        hEfficiency[i]->Sumw2();

        for (int j = 0; j < Npt; ++j)
        {

            double lowPt = pT_bins[j];
            double highPt = pT_bins[j + 1];

            // Integrate each pT interval without double-counting
            // bins at adjacent pT boundaries.
            int genLow = hGenProjection->GetXaxis()->FindBin(lowPt + 1e-6);
            int genHigh = hGenProjection->GetXaxis()->FindBin(highPt - 1e-6);

            int recLow = hRecProjection->GetXaxis()->FindBin(lowPt + 1e-6);
            int recHigh = hRecProjection->GetXaxis()->FindBin(highPt - 1e-6);

            double ngen = hGenProjection->Integral(genLow, genHigh);
            double nrec = hRecProjection->Integral(recLow, recHigh);

            double efficiency = 0.0;
            double efficiencyError = 0.0;

            if (ngen > 0.0)
            {

                efficiency = nrec / ngen;

                // Same binomial-style uncertainty prescription
                // as in the original efficiency.cxx.
                efficiencyError = sqrt(abs(
                    ((nrec + 1.0) / (ngen + 2.0)) *
                    ((nrec + 2.0) / (ngen + 3.0) -
                     (nrec + 1.0) / (ngen + 2.0))));
            }

            hEfficiency[i]->SetBinContent(j + 1, efficiency);
            hEfficiency[i]->SetBinError(j + 1, efficiencyError);
        }

        // Reset multiplicity ranges.
        hGenPt->GetAxis(1)->SetRange(0, 0);
        hRecPt->GetAxis(1)->SetRange(0, 0);

        // ================= EVENT LOSS ====================
        // Numerator: reconstructed events
        // Denominator: generated events

        int eventDenLow = hEventDen->GetXaxis()->FindBin(multLow + 0.001);
        int eventDenHigh = hEventDen->GetXaxis()->FindBin(multHigh - 0.001);

        int eventNumLow = hEventNum->GetXaxis()->FindBin(multLow + 0.001);
        int eventNumHigh = hEventNum->GetXaxis()->FindBin(multHigh - 0.001);

        double nGenEvents = hEventDen->Integral(eventDenLow, eventDenHigh);
        double nRecEvents = hEventNum->Integral(eventNumLow, eventNumHigh);

        double eventLoss = 0.0;

        if (nGenEvents > 0.0)
        {
            eventLoss = nRecEvents / nGenEvents;
        }
        else
        {
            cerr << "WARNING: Generated event count is zero in "
                 << files[i] << endl;
        }

        hEventLoss[i] = new TH1D(Form("hEventLoss_%d", i), ";#it{p}_{T} (GeV/#it{c});Event loss factor", Npt, pT_bins);
        hEventLoss[i]->SetDirectory(nullptr);

        // Event loss is a single minimum-bias factor,
        // displayed as a constant across pT.
        for (int j = 1; j <= Npt; ++j)
        {
            hEventLoss[i]->SetBinContent(j, eventLoss);
            hEventLoss[i]->SetBinError(j, 0.0);
        }

        // ================= SIGNAL LOSS ===================
        // Numerator: generated K*0 in reconstructed events
        // Denominator: all generated K*0

        int signalNumLow = hSignalNum->GetYaxis()->FindBin(multLow + 0.01);
        int signalNumHigh = hSignalNum->GetYaxis()->FindBin(multHigh - 0.01);

        int signalDenLow = hSignalDen->GetYaxis()->FindBin(multLow + 0.01);
        int signalDenHigh = hSignalDen->GetYaxis()->FindBin(multHigh - 0.01);

        TH1D *hSignalNumPt = hSignalNum->ProjectionX(Form("hSignalNumPt_%d", i), signalNumLow, signalNumHigh);

        TH1D *hSignalDenPt = hSignalDen->ProjectionX(Form("hSignalDenPt_%d", i), signalDenLow, signalDenHigh);

        hSignalLoss[i] = new TH1D(Form("hSignalLoss_%d", i), ";#it{p}_{T} (GeV/#it{c});Signal loss factor", Npt, pT_bins);
        hSignalLoss[i]->SetDirectory(nullptr);

        for (int j = 0; j < Npt; ++j)
        {

            double lowPt = pT_bins[j];
            double highPt = pT_bins[j + 1];

            int numLow = hSignalNumPt->GetXaxis()->FindBin(lowPt + 0.001);
            int numHigh = hSignalNumPt->GetXaxis()->FindBin(highPt - 0.001);

            int denLow = hSignalDenPt->GetXaxis()->FindBin(lowPt + 0.001);
            int denHigh = hSignalDenPt->GetXaxis()->FindBin(highPt - 0.001);

            double nSignalNum = hSignalNumPt->Integral(numLow, numHigh);
            double nSignalDen = hSignalDenPt->Integral(denLow, denHigh);

            double signalLoss = 0.0;

            if (nSignalDen > 0.0)
            {
                signalLoss = nSignalNum / nSignalDen;
            }

            hSignalLoss[i]->SetBinContent(j + 1, signalLoss);
            hSignalLoss[i]->SetBinError(j + 1, 0.0);

            // =============== EVENT LOSS / SIGNAL LOSS ===============

            hEventLossOverSignalLoss[i] = new TH1D(
                Form("hEventLossOverSignalLoss_%d", i),
                ";#it{p}_{T} (GeV/#it{c});Event loss / Signal loss",
                Npt, pT_bins);

            hEventLossOverSignalLoss[i]->SetDirectory(nullptr);

            for (int j = 1; j <= Npt; ++j)
            {
                double eventLossValue = hEventLoss[i]->GetBinContent(j);
                double signalLossValue = hSignalLoss[i]->GetBinContent(j);

                double ratio = 0.0;

                if (signalLossValue > 0.0)
                {
                    ratio = eventLossValue / signalLossValue;
                }

                hEventLossOverSignalLoss[i]->SetBinContent(j, ratio);
                hEventLossOverSignalLoss[i]->SetBinError(j, 0.0);
            }

            hEventLossOverSignalLoss[i]->SetLineColor(colors[i]);
            hEventLossOverSignalLoss[i]->SetMarkerColor(colors[i]);
            hEventLossOverSignalLoss[i]->SetMarkerStyle(markers[i]);
            hEventLossOverSignalLoss[i]->SetMarkerSize(0.9);
            hEventLossOverSignalLoss[i]->SetLineWidth(2);
        }

        // ================= STYLE ==========================
        hEfficiency[i]->SetLineColor(colors[i]);
        hEfficiency[i]->SetMarkerColor(colors[i]);
        hEfficiency[i]->SetMarkerStyle(markers[i]);
        hEfficiency[i]->SetMarkerSize(0.9);
        hEfficiency[i]->SetLineWidth(2);

        hEventLoss[i]->SetLineColor(colors[i]);
        hEventLoss[i]->SetMarkerColor(colors[i]);
        hEventLoss[i]->SetMarkerStyle(markers[i]);
        hEventLoss[i]->SetLineWidth(2);

        hSignalLoss[i]->SetLineColor(colors[i]);
        hSignalLoss[i]->SetMarkerColor(colors[i]);
        hSignalLoss[i]->SetMarkerStyle(markers[i]);
        hSignalLoss[i]->SetMarkerSize(0.9);
        hSignalLoss[i]->SetLineWidth(2);

        // ================= PRINT RESULTS =================
        cout << "\n"
             << labels[i] << endl;
        cout << "Generated events = " << nGenEvents << endl;
        cout << "Reconstructed events = " << nRecEvents << endl;
        cout << "Minimum-bias event loss = " << eventLoss << endl;

        cout << "\n pT range       Efficiency       Signal loss\n";

        for (int j = 0; j < Npt; ++j)
        {
            cout << pT_bins[j] << " - " << pT_bins[j + 1]
                 << "     "
                 << hEfficiency[i]->GetBinContent(j + 1)
                 << "     "
                 << hSignalLoss[i]->GetBinContent(j + 1)
                 << endl;
        }

        delete hGenProjection;
        delete hRecProjection;
        delete hSignalNumPt;
        delete hSignalDenPt;

        f->Close();
        delete f;
    }

    // ======================================================
    // DRAW EFFICIENCY COMPARISON
    // ======================================================
    TCanvas *cEff = new TCanvas("cEff", "Efficiency comparison", 850, 700);
    TLegend *legEff = new TLegend(0.58, 0.68, 0.88, 0.88);
    legEff->SetBorderSize(0);
    legEff->SetFillStyle(0);

    bool first = true;

    for (int i = 0; i < nFiles; ++i)
    {
        if (!hEfficiency[i])
            continue;

        hEfficiency[i]->SetMinimum(0.0);
        hEfficiency[i]->SetMaximum(1.15);

        hEfficiency[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hEfficiency[i]->GetYaxis()->SetTitle("Acceptance #times efficiency");

        hEfficiency[i]->Draw(first ? "E1" : "E1 SAME");
        legEff->AddEntry(hEfficiency[i], labels[i].c_str(), "lp");

        first = false;
    }

    legEff->Draw();
    cEff->SaveAs("Compare_Efficiency.png");

    // ======================================================
    // DRAW EVENT LOSS COMPARISON
    // ======================================================
    TCanvas *cEvent = new TCanvas("cEvent", "Event loss comparison", 850, 700);
    TLegend *legEvent = new TLegend(0.58, 0.68, 0.88, 0.88);
    legEvent->SetBorderSize(0);
    legEvent->SetFillStyle(0);

    first = true;

    for (int i = 0; i < nFiles; ++i)
    {
        if (!hEventLoss[i])
            continue;

        hEventLoss[i]->SetMinimum(0.0);
        hEventLoss[i]->SetMaximum(1.15);

        hEventLoss[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hEventLoss[i]->GetYaxis()->SetTitle("Event loss factor");

        hEventLoss[i]->Draw(first ? "HIST" : "HIST SAME");
        legEvent->AddEntry(hEventLoss[i],
                           Form("%s: %.5f", labels[i].c_str(),
                                hEventLoss[i]->GetBinContent(1)),
                           "l");

        first = false;
    }

    legEvent->Draw();
    cEvent->SaveAs("Compare_EventLoss.png");

    // ======================================================
    // DRAW SIGNAL LOSS COMPARISON
    // ======================================================
    TCanvas *cSignal = new TCanvas("cSignal", "Signal loss comparison", 850, 700);
    TLegend *legSignal = new TLegend(0.58, 0.68, 0.88, 0.88);
    legSignal->SetBorderSize(0);
    legSignal->SetFillStyle(0);

    first = true;

    for (int i = 0; i < nFiles; ++i)
    {
        if (!hSignalLoss[i])
            continue;

        hSignalLoss[i]->SetMinimum(0.0);
        hSignalLoss[i]->SetMaximum(1.15);

        hSignalLoss[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hSignalLoss[i]->GetYaxis()->SetTitle("Signal loss factor");

        hSignalLoss[i]->Draw(first ? "E1" : "E1 SAME");
        legSignal->AddEntry(hSignalLoss[i], labels[i].c_str(), "lp");

        first = false;
    }

    legSignal->Draw();
    cSignal->SaveAs("Compare_SignalLoss.png");

    // ======================================================
    // DRAW EVENT LOSS / SIGNAL LOSS COMPARISON
    // ======================================================

    TCanvas *cRatio = new TCanvas(
        "cRatio", "Event loss / Signal loss", 850, 700);

    TLegend *legRatio = new TLegend(0.58, 0.68, 0.88, 0.88);
    legRatio->SetBorderSize(0);
    legRatio->SetFillStyle(0);

    bool firstRatio = true;

    for (int i = 0; i < nFiles; ++i)
    {
        if (!hEventLossOverSignalLoss[i])
            continue;

        hEventLossOverSignalLoss[i]->SetMinimum(0.0);
        hEventLossOverSignalLoss[i]->SetMaximum(1.15);

        hEventLossOverSignalLoss[i]->GetXaxis()->SetTitle(
            "#it{p}_{T} (GeV/#it{c})");

        hEventLossOverSignalLoss[i]->GetYaxis()->SetTitle(
            "Event loss / Signal loss");

        hEventLossOverSignalLoss[i]->Draw(
            firstRatio ? "E1" : "E1 SAME");

        legRatio->AddEntry(
            hEventLossOverSignalLoss[i], labels[i].c_str(), "lp");

        firstRatio = false;
    }

    legRatio->Draw();
    cRatio->SaveAs("Compare_EventLossOverSignalLoss.png");

    // ======================================================
    // RATIO TO MC1: (EVENT LOSS / SIGNAL LOSS)_MCi
    //             / (EVENT LOSS / SIGNAL LOSS)_MC1
    // ======================================================

    vector<TH1D *> hRatioToMC1(nFiles, nullptr);

    if (!hEventLossOverSignalLoss[0])
    {
        cerr << "ERROR: MC1 event-loss/signal-loss histogram is missing!"
             << endl;
    }
    else
    {
        for (int i = 0; i < nFiles; ++i)
        {
            if (!hEventLossOverSignalLoss[i])
                continue;

            hRatioToMC1[i] = (TH1D *)hEventLossOverSignalLoss[i]->Clone(
                Form("hRatioToMC1_%d", i));

            hRatioToMC1[i]->SetDirectory(nullptr);
            hRatioToMC1[i]->Divide(hEventLossOverSignalLoss[0]);

            hRatioToMC1[i]->SetLineColor(colors[i]);
            hRatioToMC1[i]->SetMarkerColor(colors[i]);
            hRatioToMC1[i]->SetMarkerStyle(markers[i]);
            hRatioToMC1[i]->SetMarkerSize(0.9);
            hRatioToMC1[i]->SetLineWidth(2);
        }
    }

    // ================= DRAW RATIO TO MC1 =================

    TCanvas *cRatioMC1 = new TCanvas(
        "cRatioMC1", "Ratio to MC1", 850, 700);

    TLegend *legRatioMC1 = new TLegend(0.58, 0.68, 0.88, 0.88);
    legRatioMC1->SetBorderSize(0);
    legRatioMC1->SetFillStyle(0);

    TH1F *frameMC1 = cRatioMC1->DrawFrame(
        pT_bins[0], 0.5, pT_bins[Npt], 1.5);

    frameMC1->SetTitle("");
    frameMC1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    frameMC1->GetYaxis()->SetTitle(
        "(Event loss / Signal loss) / MC1");

    TLine *unityLine = new TLine(
        pT_bins[0], 1.0, pT_bins[Npt], 1.0);

    unityLine->SetLineStyle(2);
    unityLine->SetLineColor(kGray + 2);
    unityLine->Draw();

    for (int i = 0; i < nFiles; ++i)
    {
        if (!hRatioToMC1[i])
            continue;

        hRatioToMC1[i]->Draw("E1 SAME");
        legRatioMC1->AddEntry(hRatioToMC1[i], labels[i].c_str(), "lp");
    }

    legRatioMC1->Draw();
    cRatioMC1->SaveAs("Compare_EventLossOverSignalLoss_RatioToMC1.png");

    cout << "\nComparison completed.\n";
    cout << "Saved plots:\n";
    cout << "  Compare_Efficiency.png\n";
    cout << "  Compare_EventLoss.png\n";
    cout << "  Compare_SignalLoss.png\n";
}
