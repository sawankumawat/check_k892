#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TAxis.h"

using namespace std;

void Plot_dNdy_vs_Nch()
{
    // Open the output file generated from the tree scan
    TFile *infile = TFile::Open("yield_outfile_pp13TeV2.root", "READ");
    if (!infile || infile->IsZombie())
    {
        cout << "Error: Cannot open yield_outfile_pp13TeV.root!" << endl;
        return;
    }

    // 1. Get the 2D correlation histogram (X-axis: FT0 forward, Y-axis: Mid-rapidity Nch)
    TH2D *hNch2D = (TH2D *)infile->Get("hNch_opt0");
    if (!hNch2D)
    {
        cout << "Error: Missing hNch_opt0 correlation histogram!" << endl;
        return;
    }

    // 2. Project FT0 forward distribution (X-axis) to set up percentile boundaries
    TH1D *hFT0 = hNch2D->ProjectionX("hFT0_proj");
    double totalEvents = hFT0->Integral();
    if (totalEvents <= 0)
    {
        cout << "Error: FT0 projection has zero events!" << endl;
        return;
    }

    // Define ALICE Run 3 style percentile classes
    vector<float> percentiles = {0, 1, 2, 3, 4, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    const int nPerc = percentiles.size();

    vector<int> boundaryBins(nPerc, -1);
    boundaryBins[0] = hFT0->GetNbinsX();

    double cumulativePercent = 0.0;
    int percIdx = 1;

    // Scan from highest FT0 multiplicity down to zero
    for (int ibin = hFT0->GetNbinsX(); ibin >= 1 && percIdx < nPerc; --ibin)
    {
        cumulativePercent += hFT0->GetBinContent(ibin) * 100.0 / totalEvents;
        while (percIdx < nPerc && cumulativePercent >= percentiles[percIdx])
        {
            boundaryBins[percIdx] = ibin;
            percIdx++;
        }
    }
    boundaryBins[nPerc - 1] = 1;

    // Correct unassigned boundaries
    for (int i = 1; i < nPerc - 1; ++i)
    {
        if (boundaryBins[i] < 1)
            boundaryBins[i] = boundaryBins[i - 1];
    }

    // Particle Channels: 0: Pion, 1: Kaon, 3: Proton
    const int nParticles = 3;
    int channelIndices[nParticles] = {0, 1, 3};
    string particleNames[nParticles] = {"#pi", "K", "p"};
    int colors[nParticles] = {kRed + 1, kBlue + 1, kGreen + 2};
    int markers[nParticles] = {20, 21, 22};

    TGraphErrors *gYieldVsNch[nParticles];
    for (int i = 0; i < nParticles; i++)
    {
        gYieldVsNch[i] = new TGraphErrors(nPerc - 1);
        gYieldVsNch[i]->SetMarkerStyle(markers[i]);
        gYieldVsNch[i]->SetMarkerSize(1.2);
        gYieldVsNch[i]->SetMarkerColor(colors[i]);
        gYieldVsNch[i]->SetLineColor(colors[i]);
        gYieldVsNch[i]->SetLineWidth(2);
    }

    double eta_window = 1.0; // |eta| < 0.5 -> delta_eta = 1.0
    double y_window = 1.0;   // |y| < 0.5   -> delta_y   = 1.0

    // 3. Loop over Percentile Classes
    for (int icent = 0; icent < nPerc - 1; ++icent)
    {
        int binHigh = boundaryBins[icent];
        int binLow = boundaryBins[icent + 1];
        if (binLow > binHigh)
            swap(binLow, binHigh);

        // Compute mean mid-rapidity Nch (<dNch/deta>) for this FT0 slice
        TH1D *hNchMidProj = hNch2D->ProjectionY(Form("hMid_%d", icent), binLow, binHigh);
        double meanNch = hNchMidProj->GetMean() / eta_window; // Correct mean calculation
        double nEvents = hFT0->Integral(binLow, binHigh);
        delete hNchMidProj;

        if (nEvents <= 0)
            continue;

        // Extract yields for each particle species inside this percentile slice
        for (int i = 0; i < nParticles; i++)
        {
            int ich = channelIndices[i];
            TH2D *hParticleVsFT0 = (TH2D *)infile->Get(Form("hDec_opt0_ch%d", ich));
            if (!hParticleVsFT0)
                continue;

            // Project pT distribution for events falling into FT0 bin slice [binLow, binHigh]
            TH1D *hPt = hParticleVsFT0->ProjectionY(Form("hPt_%d_%d", i, icent), binLow, binHigh);
            double rawYield = hPt->Integral();
            delete hPt;

            double dNdy = rawYield / (nEvents * y_window);
            double dNdy_err = sqrt(rawYield) / (nEvents * y_window);

            gYieldVsNch[i]->SetPoint(icent, meanNch, dNdy);
            gYieldVsNch[i]->SetPointError(icent, 0.0, dNdy_err);
        }
    }

    TFile *fOutput = new TFile("dNdy_vs_Nch_PercentileSliced.root", "RECREATE");
    for (int i = 0; i < nParticles; i++)
    {
        gYieldVsNch[i]->Write(Form("gYieldVsNch_%s", particleNames[i].c_str()));
    }

    // 4. Plotting Setup
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "Corrected dN/dy vs Multiplicity", 800, 600);
    c1->SetMargin(0.12, 0.05, 0.12, 0.05);
    c1->SetGrid();

    TH1F *hFrame = c1->DrawFrame(0.0, 0.001, 55.0, 25.0);
    hFrame->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT_{|#eta|<0.5}");
    hFrame->GetYaxis()->SetTitle("dN/dy");
    hFrame->GetXaxis()->SetTitleSize(0.045);
    hFrame->GetYaxis()->SetTitleSize(0.045);

    TLegend *leg = new TLegend(0.20, 0.70, 0.40, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    for (int i = 0; i < nParticles; i++)
    {
        gYieldVsNch[i]->Draw("PE1 SAME");
        leg->AddEntry(gYieldVsNch[i], particleNames[i].c_str(), "pe");
    }

    leg->Draw();
    // c1->SaveAs("dNdy_vs_Nch_PercentileSliced.png");
    cout << "Plot successfully generated with percentile slicing!" << endl;
}