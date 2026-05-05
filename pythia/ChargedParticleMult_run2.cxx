#include <iostream>
#include <vector>
#include "../macro/src/style.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void ChargedParticleMult_run2()
{
    gStyle->SetOptStat(0);

    int colors[] = {kRed, kOrange, kYellow, kGreen, kCyan, kBlue, kViolet, kMagenta, kPink - 1, kGray};

    TFile *fpp13TeV = new TFile("NewSimulationChargeMult/13TeV.root", "READ");
    if (fpp13TeV->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    TH2D *hmultFT0vskstar_pt = (TH2D *)fpp13TeV->Get("hmultFT0vskstar_pt");
    TH2D *hmultFT0vsphi_pt = (TH2D *)fpp13TeV->Get("hmultFT0vsphi_pt");
    TH2D *hmultV0vshmultV0Eta0p5 = (TH2D *)fpp13TeV->Get("hmultV0vshmultV0Eta0p5");
    TH1D *hmultpp13 = (TH1D *)hmultV0vshmultV0Eta0p5->ProjectionX("hmultpp13");
    // TH1D *hmultpp13 = (TH1D *)fpp13TeV->Get("hmultFT0"); // Less bins (do not use)

    if (hmultpp13 == nullptr || hmultV0vshmultV0Eta0p5 == nullptr)
    {
        cout << "Error reading histograms" << endl;
        return;
    }

    double totalEvents = hmultpp13->Integral(1, hmultpp13->GetNbinsX());

    // ALICE-style percentiles
    double percentiles[10] = {0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00};
    // double percentiles[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    vector<int> regionEdges;
    double cumulative = 0;
    int idx = 0;

    // Integrate from HIGH multiplicity → LOW
    for (int i = hmultpp13->GetNbinsX(); i >= 1; i--)
    {
        cumulative += hmultpp13->GetBinContent(i);

        while (idx < 10 && cumulative >= percentiles[idx] * totalEvents)
        {
            regionEdges.push_back(i);
            idx++;
        }
    }

    // Safety (in case of rounding)
    cout << "Size of regionEdges: " << regionEdges.size() << endl;
    while ((int)regionEdges.size() < 10)
        regionEdges.push_back(1);

    // Create region histograms
    TH1D *hRegions[10];
    for (int i = 0; i < 10; i++)
    {
        hRegions[i] = (TH1D *)hmultpp13->Clone(Form("hRegion_%d", i));
        hRegions[i]->Reset();
        hRegions[i]->SetFillColor(colors[i]);
        hRegions[i]->SetLineColor(kBlack);
        hRegions[i]->SetFillStyle(3001);
    }

    // Fill regions (high → low)
    int prevBin = hmultpp13->GetNbinsX();

    for (int r = 0; r < 10; r++)
    {
        int lowBin = regionEdges[r];

        for (int b = lowBin; b <= prevBin; b++)
        {
            hRegions[r]->SetBinContent(b, hmultpp13->GetBinContent(b));
        }

        prevBin = lowBin - 1;
    }

    // -------------------------------
    // Plot multiplicity distribution
    // -------------------------------
    TCanvas *cMult = new TCanvas("cMult", "Multiplicity", 720, 720);
    SetCanvasStyle(cMult, 0.14, 0.01, 0.01, 0.14);
    SetHistoQA(hmultpp13);
    gPad->SetLogy();
    gPad->SetLogx();
    hmultpp13->SetMaximum(hmultpp13->GetMaximum() * 20);
    hmultpp13->GetXaxis()->SetTitle("N_{ch}");
    hmultpp13->GetYaxis()->SetTitle("Events");
    hmultpp13->Draw();

    const char *labels[10] = {"0-1%", "1-5%", "5-10%", "10-15%", "15-20%", "20-30%", "30-40%", "40-50%", "50-70%", "70-100%"};

    TLegend *leg = new TLegend(0.3, 0.87, 0.95, 0.98);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->SetNColumns(5);

    for (int i = 0; i < 10; i++)
    {
        hRegions[i]->Draw("HIST SAME");
        leg->AddEntry(hRegions[i], labels[i], "f");

        const double eventsInClass = hRegions[i]->Integral(1, hRegions[i]->GetNbinsX());
        cout << labels[i] << " events = " << eventsInClass
             << " (" << 100.0 * eventsInClass / totalEvents << "%)" << endl;
    }
    leg->Draw();

    // // -------------------------------------------
    // // Plot multiplicity percentiles
    // // -------------------------------------------

    // TCanvas *cPercentiles = new TCanvas("cPercentiles", "Multiplicity Percentiles", 720, 720);
    // SetCanvasStyle(cPercentiles, 0.14, 0.03, 0.06, 0.14);
    // TH1D *hMultiplicityPercentile = new TH1D("hMultiplicityPercentile", "Multiplicity Percentiles", 10, 0, 10);
    // for (int i = 0; i < 10; i++)
    // {
    //     hMultiplicityPercentile->SetBinContent(i + 1, hRegions[i]->Integral(1, hRegions[i]->GetNbinsX()));
    //     hMultiplicityPercentile->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    // }
    // hMultiplicityPercentile->GetXaxis()->SetTitle("Percentile class");
    // hMultiplicityPercentile->GetYaxis()->SetTitle("Events");
    // SetHistoQA(hMultiplicityPercentile);
    // hMultiplicityPercentile->SetMinimum(0);
    // hMultiplicityPercentile->Draw("HIST");

    // -------------------------------------------
    // ⟨dN/dη⟩ in |η| < 0.5 for each percentile
    // -------------------------------------------
    cout << "\n===== <dN/deta>|eta|<0.5 =====" << endl;

    int highBin = hmultpp13->GetNbinsX();
    TH1D *hChMultEta0p5[10];

    for (int i = 0; i < 10; i++)
    {
        int lowBin = regionEdges[i];

        hmultV0vshmultV0Eta0p5->GetXaxis()->SetRange(lowBin, highBin);

        // TH1D *proj = (TH1D *)hmultV0vshmultV0Eta0p5->ProjectionY(Form("proj_%d", i));
        hChMultEta0p5[i] = (TH1D *)hmultV0vshmultV0Eta0p5->ProjectionY(Form("hChMultEta0p5_%d", i));

        cout << labels[i] << " : " << hChMultEta0p5[i]->GetMean() << endl;

        highBin = lowBin - 1;
    }

    TCanvas *cChMultEta0p5 = new TCanvas("cChMultEta0p5", "cChMultEta0p5", 720, 720);
    SetCanvasStyle(cChMultEta0p5, 0.14, 0.03, 0.01, 0.14);
    gPad->SetLogy();
    // gPad->SetLogx();
    for (int i = 0; i < 10; i++)
    {
        SetHistoQA(hChMultEta0p5[i]);
        hChMultEta0p5[i]->GetXaxis()->SetRangeUser(0, 75);
        hChMultEta0p5[i]->SetMaximum(20e7);
        hChMultEta0p5[i]->SetLineColor(colors[i]);
        hChMultEta0p5[i]->SetFillColor(colors[i]);
        hChMultEta0p5[i]->SetFillStyle(3001);
        hChMultEta0p5[i]->GetXaxis()->SetTitle("dN/d#eta in |#eta|<0.5");
        hChMultEta0p5[i]->GetYaxis()->SetTitle("Entries");
        hChMultEta0p5[i]->Draw("same");
    }
    leg->Draw();

    // Get the K* pT for the given multiplicity classes
    TCanvas *cKstarPt = new TCanvas("cKstarPt", "Kstar pT vs Multiplicity", 720, 720);
    SetCanvasStyle(cKstarPt, 0.14, 0.03, 0.01, 0.14);
    gPad->SetLogy();
    highBin = hmultpp13->GetNbinsX();
    TH1D *hKstarPt[10];
    const double deltaY = 1.0; // Acceptance |y| < 0.5

    for (int i = 0; i < 10; i++)
    {
        int lowBin = regionEdges[i];
        hmultFT0vskstar_pt->GetXaxis()->SetRange(lowBin, highBin);

        hKstarPt[i] = (TH1D *)hmultFT0vskstar_pt->ProjectionY(Form("hKstarPt_%d", i));
        const double classEvents = hRegions[i]->Integral(1, hRegions[i]->GetNbinsX());
        const double rawYield = hKstarPt[i]->Integral(1, hKstarPt[i]->GetNbinsX());
        const double dNdy = (classEvents > 0.0) ? rawYield / (classEvents * deltaY) : 0.0;

        cout << "dN_{ch}/d#eta " << hChMultEta0p5[i]->GetMean()
             << " : <pT> = " << hKstarPt[i]->GetMean()
             << ", raw yield = " << rawYield
             << ", N_events(class) = " << classEvents
             << ", dN/dy = " << dNdy << endl;
        highBin = lowBin - 1;
    }
    SetHistoQA(hKstarPt[0]);
    hKstarPt[0]->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    hKstarPt[0]->GetYaxis()->SetTitle("Entries");
    hKstarPt[0]->Draw();
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.05);
    pad2->SetRightMargin(0.05);
    pad2->SetBottomMargin(0.35);
    pad1->SetLeftMargin(0.175);
    pad2->SetLeftMargin(0.175);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.003);
    pad2->SetTopMargin(0.04);

    // Set ticks on individual pads
    pad1->SetTicks(1, 1);
    pad2->SetTicks(1, 1);
}