#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TLatex.h"
#include "src/style.h"

// Helper function to apply custom X-axis labels
void SetCustomBinLabels(TH1D *h)
{
    std::vector<TString> labels = {
        "All events",
        "TriggerTVX",
        "|V_{Z}|<10 cm",
        "NoTFBorder",
        "NoITSROF"};

    for (size_t i = 0; i < labels.size(); ++i)
    {
        h->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    }
}

// Helper function to draw percentage labels on top of each bin
void DrawBinPercentages(TH1D *hInput)
{
    double firstBinContent = hInput->GetBinContent(1);
    if (firstBinContent <= 0)
        return;

    TLatex text;
    text.SetTextSize(0.035);
    text.SetTextAlign(21); // Centered horizontally, bottom-aligned vertically
    text.SetTextColor(kBlack);

    for (int i = 1; i <= 5; ++i) // Only process the first 5 bins
    {
        double content = hInput->GetBinContent(i);
        double percentage = (content / firstBinContent) * 100.0;

        double xPos = hInput->GetBinCenter(i);
        double yPos = content + (hInput->GetMaximum() * 0.02); // Position slightly above the bar

        text.DrawLatex(xPos, yPos, Form("%.1f%%", percentage));
    }
}

void compareEventsDataMC()
{
    gStyle->SetOptStat(0);

    TString dataPath = "/home/sawan/Storage/check_k892/data/kstar/LHC22o_pass7";
    TString mcPath = "/home/sawan/Storage/check_k892/mc/LHC24f3c";

    TFile *fData = TFile::Open(dataPath + "/794581.root", "READ");
    TFile *fMC = TFile::Open(mcPath + "/793191.root", "READ");

    if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie())
    {
        std::cout << "Error opening files!" << std::endl;
        return;
    }

    TH1D *hEventsData = (TH1D *)fData->Get("kstarqa/eventSelection/hEventCut");
    TH1D *hEventsMC = (TH1D *)fMC->Get("kstarqa/eventSelection/hEventCut");

    if (!hEventsData || !hEventsMC)
    {
        std::cout << "Error retrieving histograms!" << std::endl;
        return;
    }

    // Set custom labels for the bins
    SetCustomBinLabels(hEventsData);
    SetCustomBinLabels(hEventsMC);

    // Create Canvas with 2 separate pads side by side
    TCanvas *cEvents = new TCanvas("cEvents", "Event Selection Comparison", 1400, 600);
    cEvents->Divide(2, 1);

    // --- Pad 1: Data ---
    cEvents->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.22); // Increased margin so vertical labels don't get cut off
    gPad->SetGridx(1);
    SetHistoQA(hEventsData);

    // Set view range to cover the 5 bins (bins 1 to 5)
    hEventsData->GetXaxis()->SetRange(1, 5);
    hEventsData->SetLineColor(kRed);
    hEventsData->SetLineWidth(2);

    // Set X-axis labels to vertical (90 degrees)
    hEventsData->GetXaxis()->LabelsOption("v");
    hEventsData->GetXaxis()->SetLabelSize(0.045);

    hEventsData->GetYaxis()->SetTitle("Events");
    hEventsData->GetYaxis()->SetRangeUser(0, hEventsData->GetMaximum() * 1.15);
    hEventsData->SetTitle("Data - Event Selection");
    hEventsData->Draw("HIST");
    DrawBinPercentages(hEventsData);

    // --- Pad 2: MC ---
    cEvents->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.22); // Increased margin so vertical labels don't get cut off
    gPad->SetGridx(1);
    SetHistoQA(hEventsMC);

    // Set view range to cover the 5 bins (bins 1 to 5)
    hEventsMC->GetXaxis()->SetRange(1, 5);
    hEventsMC->SetLineColor(kBlue);
    hEventsMC->SetLineWidth(2);

    // Set X-axis labels to vertical (90 degrees)
    hEventsMC->GetXaxis()->LabelsOption("v");
    hEventsMC->GetXaxis()->SetLabelSize(0.045);

    hEventsMC->GetYaxis()->SetTitle("Events");
    hEventsMC->GetYaxis()->SetRangeUser(0, hEventsMC->GetMaximum() * 1.15);
    hEventsMC->SetTitle("MC - Event Selection");
    hEventsMC->Draw("HIST");
    DrawBinPercentages(hEventsMC);

    cEvents->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/EventSelectionComparison.pdf");

    //Also make Vz distribution from MC

    TCanvas *cVz = new TCanvas("cVz", "Vz Distribution", 720, 720);
    SetCanvasStyle(cVz, 0.16, 0.05, 0.06, 0.15);
    TH1D *hVzMC = (TH1D *)fMC->Get("kstarqa/eventSelection/hVertexZRec");
    SetHistoQA(hVzMC);
    hVzMC->GetXaxis()->SetTitle("Vertex Z (cm)");
    hVzMC->GetYaxis()->SetTitle("Counts");
    hVzMC->Draw();
    cVz->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/VertexZ_MC.pdf");
}