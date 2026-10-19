#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>

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
#include "TLine.h"

using namespace std;

void CompareEfficiencyMC_MultBins()
{
    // ================= MC FILES =================
    vector<string> files = {
        "/home/sawan/Storage/check_k892/mc/LHC24f3c/750013.root",
        "/home/sawan/Storage/check_k892/mc/LHC24f3c/CalibratedMCchecks/780946.root"
        // "/home/sawan/Storage/check_k892/mc/LHC24f3c/CalibratedMCchecks/780947.root",
        // "/home/sawan/Storage/check_k892/mc/LHC24f3c/CalibratedMCchecks/780948.root"
    };

    vector<string> labels = {
        "Old MC calibration", "lf-mc-centrality (new)"
        // "MC 3", "MC 4", "MC 5"
    };
    const int nFiles = files.size();

    if ((int)labels.size() != nFiles)
    {
        cerr << "ERROR: files and labels have different sizes." << endl;
        return;
    }

    int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 1, kOrange + 7};
    int markers[] = {20, 21, 22, 23, 24, 25};

    // ================= pT BINNING =================
    const int Npt = 23;
    double pT_bins[Npt + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};

    // ================= MULTIPLICITY CLASSES =================
    const int Nmult = 10;
    double multBins[Nmult + 1] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};

    // Histogram directory for each corresponding file.
    vector<TString> genPaths = {
        "kstarqa/hInvMass",
        "kstarqa_ServiceWagon/hInvMass"
        // "kstarqa_serviceWagon_recalibMode0/hInvMass",
        // "kstarqa_serviceWagon_recalibMode1/hInvMass"
    };

    if ((int)genPaths.size() != nFiles)
    {
        cerr << "ERROR: genPaths and files have different sizes." << endl;
        return;
    }

    gStyle->SetOptStat(0);

    // Histograms are indexed [multiplicity class][MC file].
    vector<vector<TH1D *>> hEfficiency(Nmult, vector<TH1D *>(nFiles, nullptr));
    vector<vector<TH1D *>> hEventLoss(Nmult, vector<TH1D *>(nFiles, nullptr));
    vector<vector<TH1D *>> hSignalLoss(Nmult, vector<TH1D *>(nFiles, nullptr));
    vector<vector<TH1D *>> hEventLossOverSignalLoss(Nmult, vector<TH1D *>(nFiles, nullptr));

    // ======================================================
    // PROCESS EACH MC FILE
    // ======================================================
    for (int i = 0; i < nFiles; ++i)
    {
        cout << "\n====================================\nProcessing: " << files[i] << "\n====================================\n";

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

        auto *hGenPt = dynamic_cast<THnSparse *>(f->Get(genPtPath));
        auto *hRecPt = dynamic_cast<THnSparse *>(f->Get(recPtPath));
        auto *hEventDen = dynamic_cast<TH1 *>(f->Get(eventDenPath));
        auto *hEventNum = dynamic_cast<TH1 *>(f->Get(eventNumPath));
        auto *hSignalDen = dynamic_cast<TH2 *>(f->Get(signalDenPath));
        auto *hSignalNum = dynamic_cast<TH2 *>(f->Get(signalNumPath));

        if (!hGenPt || !hRecPt || !hEventDen || !hEventNum || !hSignalDen || !hSignalNum)
        {
            cerr << "ERROR: One or more required histograms are missing in " << files[i] << endl;
            f->Close();
            delete f;
            continue;
        }

        for (int m = 0; m < Nmult; ++m)
        {
            const double multLow = multBins[m];
            const double multHigh = multBins[m + 1];
            const TString multLabel = Form("%.0f-%.0f%%", multLow, multHigh);

            cout << "\n"
                 << labels[i] << " multiplicity " << multLabel << endl;

            // ================= EFFICIENCY =====================
            int genMultLow = hGenPt->GetAxis(1)->FindBin(multLow + 1e-4);
            int genMultHigh = hGenPt->GetAxis(1)->FindBin(multHigh - 1e-4);
            int recMultLow = hRecPt->GetAxis(1)->FindBin(multLow + 1e-4);
            int recMultHigh = hRecPt->GetAxis(1)->FindBin(multHigh - 1e-4);

            hGenPt->GetAxis(1)->SetRange(genMultLow, genMultHigh);
            hRecPt->GetAxis(1)->SetRange(recMultLow, recMultHigh);

            TH1D *hGenProjection = hGenPt->Projection(0, "E");
            TH1D *hRecProjection = hRecPt->Projection(0, "E");

            hEfficiency[m][i] = new TH1D(Form("hEfficiency_m%d_mc%d", m, i),
                                         Form("%s;#it{p}_{T} (GeV/#it{c});Acceptance #times efficiency", multLabel.Data()), Npt, pT_bins);
            hEfficiency[m][i]->SetDirectory(nullptr);
            hEfficiency[m][i]->Sumw2();

            for (int j = 0; j < Npt; ++j)
            {
                double lowPt = pT_bins[j];
                double highPt = pT_bins[j + 1];
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
                    efficiencyError = sqrt(abs(((nrec + 1.0) / (ngen + 2.0)) *
                                               ((nrec + 2.0) / (ngen + 3.0) - (nrec + 1.0) / (ngen + 2.0))));
                }
                hEfficiency[m][i]->SetBinContent(j + 1, efficiency);
                hEfficiency[m][i]->SetBinError(j + 1, efficiencyError);
                hEfficiency[m][i]->SetMaximum(1.02);
                hEfficiency[m][i]->SetMinimum(0.0);
            }

            delete hGenProjection;
            delete hRecProjection;
            hGenPt->GetAxis(1)->SetRange(0, 0);
            hRecPt->GetAxis(1)->SetRange(0, 0);

            // ================= EVENT LOSS ====================
            int eventDenLow = hEventDen->GetXaxis()->FindBin(multLow + 1e-4);
            int eventDenHigh = hEventDen->GetXaxis()->FindBin(multHigh - 1e-4);
            int eventNumLow = hEventNum->GetXaxis()->FindBin(multLow + 1e-4);
            int eventNumHigh = hEventNum->GetXaxis()->FindBin(multHigh - 1e-4);
            double nGenEvents = hEventDen->Integral(eventDenLow, eventDenHigh);
            double nRecEvents = hEventNum->Integral(eventNumLow, eventNumHigh);
            double eventLoss = (nGenEvents > 0.0) ? nRecEvents / nGenEvents : 0.0;

            hEventLoss[m][i] = new TH1D(Form("hEventLoss_m%d_mc%d", m, i),
                                        Form("%s;#it{p}_{T} (GeV/#it{c});Event loss factor", multLabel.Data()), Npt, pT_bins);
            hEventLoss[m][i]->SetDirectory(nullptr);
            for (int j = 1; j <= Npt; ++j)
                hEventLoss[m][i]->SetBinContent(j, eventLoss);
            hEventLoss[m][i]->SetMaximum(1.22);
            hEventLoss[m][i]->SetMinimum(0.35);

            // ================= SIGNAL LOSS ===================
            int signalNumLow = hSignalNum->GetYaxis()->FindBin(multLow + 1e-4);
            int signalNumHigh = hSignalNum->GetYaxis()->FindBin(multHigh - 1e-4);
            int signalDenLow = hSignalDen->GetYaxis()->FindBin(multLow + 1e-4);
            int signalDenHigh = hSignalDen->GetYaxis()->FindBin(multHigh - 1e-4);

            TH1D *hSignalNumPt = hSignalNum->ProjectionX(Form("hSignalNumPt_m%d_mc%d", m, i), signalNumLow, signalNumHigh);
            TH1D *hSignalDenPt = hSignalDen->ProjectionX(Form("hSignalDenPt_m%d_mc%d", m, i), signalDenLow, signalDenHigh);

            hSignalLoss[m][i] = new TH1D(Form("hSignalLoss_m%d_mc%d", m, i),
                                         Form("%s;#it{p}_{T} (GeV/#it{c});Signal loss factor", multLabel.Data()), Npt, pT_bins);
            hSignalLoss[m][i]->SetDirectory(nullptr);

            for (int j = 0; j < Npt; ++j)
            {
                double lowPt = pT_bins[j];
                double highPt = pT_bins[j + 1];
                int numLow = hSignalNumPt->GetXaxis()->FindBin(lowPt + 1e-6);
                int numHigh = hSignalNumPt->GetXaxis()->FindBin(highPt - 1e-6);
                int denLow = hSignalDenPt->GetXaxis()->FindBin(lowPt + 1e-6);
                int denHigh = hSignalDenPt->GetXaxis()->FindBin(highPt - 1e-6);
                double nSignalNum = hSignalNumPt->Integral(numLow, numHigh);
                double nSignalDen = hSignalDenPt->Integral(denLow, denHigh);
                double signalLoss = (nSignalDen > 0.0) ? nSignalNum / nSignalDen : 0.0;
                hSignalLoss[m][i]->SetBinContent(j + 1, signalLoss);
                hSignalLoss[m][i]->SetBinError(j + 1, 0.0);
                hSignalLoss[m][i]->SetMaximum(0.94);
                hSignalLoss[m][i]->SetMinimum(0.47);
            }
            delete hSignalNumPt;
            delete hSignalDenPt;

            // ============== EVENT LOSS / SIGNAL LOSS ==============
            hEventLossOverSignalLoss[m][i] = new TH1D(Form("hEventLossOverSignalLoss_m%d_mc%d", m, i),
                                                      Form("%s;#it{p}_{T} (GeV/#it{c});Event loss / Signal loss", multLabel.Data()), Npt, pT_bins);
            hEventLossOverSignalLoss[m][i]->SetDirectory(nullptr);
            for (int j = 1; j <= Npt; ++j)
            {
                double eventValue = hEventLoss[m][i]->GetBinContent(j);
                double signalValue = hSignalLoss[m][i]->GetBinContent(j);
                hEventLossOverSignalLoss[m][i]->SetBinContent(j, signalValue > 0.0 ? eventValue / signalValue : 0.0);
                hEventLossOverSignalLoss[m][i]->SetBinError(j, 0.0);
                hEventLossOverSignalLoss[m][i]->SetMaximum(1.3);
                hEventLossOverSignalLoss[m][i]->SetMinimum(0.55);
            }

            // ================= STYLE ==========================
            TH1D *histograms[] = {hEfficiency[m][i], hEventLoss[m][i], hSignalLoss[m][i], hEventLossOverSignalLoss[m][i]};
            for (TH1D *h : histograms)
            {
                h->SetLineColor(colors[i]);
                h->SetMarkerColor(colors[i]);
                h->SetMarkerStyle(markers[i]);
                h->SetMarkerSize(0.8);
                h->SetLineWidth(2);
            }
        }

        f->Close();
        delete f;
    }

    // ================= DRAW MULTI-PAD CANVAS (5x2) =================
    struct PlotSpec
    {
        const char *name;
        const char *title;
        vector<vector<TH1D *>> *hists;
        bool isEventLoss;
    };

    PlotSpec plots[] = {
        {"Efficiency", "Acceptance #times efficiency", &hEfficiency, false},
        {"SignalLoss", "Signal loss factor", &hSignalLoss, false},
        {"EventLoss", "Event loss factor", &hEventLoss, true},
        {"EventLossOverSignalLoss", "Event loss / Signal loss", &hEventLossOverSignalLoss, false}};

    for (const auto &plot : plots)
    {
        TCanvas *c = new TCanvas(Form("c_%s_Divided", plot.name), plot.title, 2000, 900);
        c->Divide(5, 2);

        for (int m = 0; m < Nmult; ++m)
        {
            c->cd(m + 1);
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.04);
            gPad->SetTopMargin(0.08);
            gPad->SetBottomMargin(0.12);

            TLegend *leg = new TLegend(0.40, 0.65, 0.90, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.035);

            bool first = true;
            for (int i = 0; i < nFiles; ++i)
            {
                TH1D *h = (*(plot.hists))[m][i];
                if (!h)
                    continue;

                h->GetXaxis()->SetTitleSize(0.045);
                h->GetYaxis()->SetTitleSize(0.045);
                h->GetXaxis()->SetLabelSize(0.04);
                h->GetYaxis()->SetLabelSize(0.04);

                h->Draw(first ? (plot.isEventLoss ? "HIST" : "E1") : (plot.isEventLoss ? "HIST SAME" : "E1 SAME"));

                TString legendLabel = labels[i].c_str();
                if (plot.isEventLoss)
                    legendLabel += Form(" (%.4f)", h->GetBinContent(1));
                leg->AddEntry(h, legendLabel, plot.isEventLoss ? "l" : "lp");

                first = false;
            }

            leg->Draw();
        }

        TString output = Form("Compare_%s_MultBins_Divided.png", plot.name);
        c->SaveAs(output);
        delete c;
    }

    cout << "\nComparison completed for all multiplicity classes.\n";
    cout << "Generated 4 divided (5x2) canvas plots.\n";
}