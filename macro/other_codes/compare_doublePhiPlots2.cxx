#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"
using namespace std;

std::vector<int> markerStyles = {20, 21, 22, 23, 29, 33, 34, 47, 48, 49};
std::vector<int> vibrantColors = {
    kBlack,      // #656364
    kBlue,       // #578dff
    kGreen + 3,  // #86c8dd
    kMagenta,    // #adad7d
    kRed,        // #c849a9
    kOrange + 3, // #c91116
    kCyan + 3,   // #1845fb
    kSpring + 3  // #f5e002
};
std::vector<int> lineStyles = {1, 2, 5, 4, 9};

//==================================Notes========================================//
// No effect of deep angle cut is seen on the phi invariant mass and fit parameters
//==================================Notes========================================//

void compare_doublePhiPlots2()
{
    gStyle->SetOptStat(0);
    const int TotalFiles = 2;
    vector<TFile *> compareFiles(TotalFiles);
    TString dataFilePath[TotalFiles];
    dataFilePath[0] = "../output/doublePhi/508501/";            //
    dataFilePath[1] = "../output/doublePhi/LHC25ac_allWagons/"; //
    TString outputFilePath = dataFilePath[1] + "ComparePlots/";
    gSystem->Exec("mkdir -p " + outputFilePath);
    TString fileNames[TotalFiles] = {"PhiPhi.root", "PhiPhi.root"};
    TString LegendNames[TotalFiles] = {"2023 Data", "Triggered Data"};
    TH1F *hPhiMassFit[TotalFiles];
    TH1F *hPhiMassResolutionFit[TotalFiles];
    TH1F *hPhiYieldFit[TotalFiles];
    TH1F *hPurity[TotalFiles];
    TH1F *hNumPhi[TotalFiles];
    TH1F *hPhiPhiInvMass[TotalFiles];

    TCanvas *cCompareMass = new TCanvas("cCompareMass", "Compare Phi Mass Fit Parameters", 720, 720);
    SetCanvasStyle(cCompareMass, 0.21, 0.03, 0.05, 0.15);
    TCanvas *cCompareMassResolution = new TCanvas("cCompareMassResolution", "Compare Phi Mass Resolution Fit Parameters", 720, 720);
    SetCanvasStyle(cCompareMassResolution, 0.15, 0.03, 0.05, 0.15);
    TCanvas *cCompareYield = new TCanvas("cCompareYield", "Compare Phi Yield Fit Parameters", 720, 720);
    SetCanvasStyle(cCompareYield, 0.15, 0.03, 0.05, 0.15);
    TCanvas *cComparePurity = new TCanvas("cComparePurity", "Compare Phi Purity", 720, 720);
    SetCanvasStyle(cComparePurity, 0.15, 0.03, 0.05, 0.15);
    TCanvas *cCompareNumPhi = new TCanvas("cCompareNumPhi", "Compare Number of #Phi mesons in an event", 720, 720);
    SetCanvasStyle(cCompareNumPhi, 0.15, 0.05, 0.05, 0.15);
    TCanvas *cComparePhiPhiInvMass = new TCanvas("cComparePhiPhiInvMass", "Compare #Phi#Phi Invariant Mass", 720, 720);
    SetCanvasStyle(cComparePhiPhiInvMass, 0.15, 0.05, 0.05, 0.15);

    TLegend *legMass = new TLegend(0.5, 0.7, 0.9, 0.91);
    legMass->SetFillStyle(0);
    legMass->SetBorderSize(0);
    legMass->SetTextFont(42);
    legMass->SetTextSize(0.03);

    TLegend *legMassRes = new TLegend(0.5, 0.74, 0.9, 0.91);
    legMassRes->SetFillStyle(0);
    legMassRes->SetBorderSize(0);
    legMassRes->SetTextFont(42);
    legMassRes->SetTextSize(0.03);

    TLegend *legMassRes2 = new TLegend(0.5, 0.74, 0.9, 0.91);
    legMassRes2->SetFillStyle(0);
    legMassRes2->SetBorderSize(0);
    legMassRes2->SetTextFont(42);
    legMassRes2->SetTextSize(0.03);

    TLegend *legPhiPhiMass = new TLegend(0.5, 0.7, 0.9, 0.91);
    legPhiPhiMass->SetFillStyle(0);
    legPhiPhiMass->SetBorderSize(0);
    legPhiPhiMass->SetTextFont(42);
    legPhiPhiMass->SetTextSize(0.03);

    for (int i = 0; i < TotalFiles; i++)
    {
        compareFiles[i] = new TFile(dataFilePath[i] + fileNames[i]);
        if (!compareFiles[i] || compareFiles[i]->IsZombie())
        {
            cerr << "Error: Could not open file " << fileNames[i] << "\n";
            return;
        }
        hNumPhi[i] = (TH1F *)compareFiles[i]->Get("hNumPhi");
        hPhiMassFit[i] = (TH1F *)compareFiles[i]->Get("FittedPhiMass");
        hPhiMassResolutionFit[i] = (TH1F *)compareFiles[i]->Get("FittedPhiMassResolution");
        hPhiYieldFit[i] = (TH1F *)compareFiles[i]->Get("FittedPhiYield");
        hPurity[i] = (TH1F *)compareFiles[i]->Get("PhiPurity");
        hPhiPhiInvMass[i] = (TH1F *)compareFiles[i]->Get("rawInvMass_PhiPhi");
        if (hPhiPhiInvMass[i] == nullptr)
        {
            cerr << "Error: Could not find one or more histograms in file " << fileNames[i] << "\n";
            return;
        }

        // /*
        cCompareMass->cd();
        SetHistoQA(hPhiMassFit[i]);
        hPhiMassFit[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hPhiMassFit[i]->GetYaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
        hPhiMassFit[i]->SetMarkerStyle(markerStyles[i]);
        hPhiMassFit[i]->SetMarkerColor(vibrantColors[i]);
        hPhiMassFit[i]->SetLineColor(vibrantColors[i]);
        hPhiMassFit[i]->SetMarkerSize(1.3);
        hPhiMassFit[i]->SetLineWidth(3);
        hPhiMassFit[i]->GetYaxis()->SetTitleOffset(2.1);
        if (i == 0)
        {
            hPhiMassFit[i]->GetYaxis()->SetRangeUser(1.0175, 1.0215);
            // hPhiMassFit[i]->SetMaximum(hPhiMassFit[i]->GetMaximum() * 1.2);
            hPhiMassFit[i]->Draw("pe");
        }
        else
        {
            hPhiMassFit[i]->Draw("pe same");
        }
        TLine *linePDG = new TLine(0.4, 1.019461, 10.0, 1.019461);
        linePDG->SetLineColor(kRed);
        linePDG->SetLineStyle(2);
        linePDG->SetLineWidth(3);
        if (TotalFiles > 1)
            legMass->AddEntry(hPhiMassFit[i], LegendNames[i], "lpe");
        if (i == TotalFiles - 1)
        {
            linePDG->Draw("same");
            legMass->AddEntry(linePDG, "#phi PDG mass", "l");
            legMass->Draw();
        }

        cCompareMassResolution->cd();
        SetHistoQA(hPhiMassResolutionFit[i]);
        hPhiMassResolutionFit[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hPhiMassResolutionFit[i]->GetYaxis()->SetTitle("#sigma_{#phi} (GeV/#it{c}^{2})");
        hPhiMassResolutionFit[i]->SetMarkerStyle(markerStyles[i]);
        hPhiMassResolutionFit[i]->SetMarkerColor(vibrantColors[i]);
        hPhiMassResolutionFit[i]->SetLineColor(vibrantColors[i]);
        hPhiMassResolutionFit[i]->SetMarkerSize(1.3);
        hPhiMassResolutionFit[i]->SetLineWidth(3);
        if (i == 0)
        {
            hPhiMassResolutionFit[i]->GetYaxis()->SetRangeUser(0.0006, 0.0043);
            // hPhiMassResolutionFit[i]->SetMaximum(hPhiMassResolutionFit[i]->GetMaximum() * 1.2);
            hPhiMassResolutionFit[i]->Draw("pe");
        }
        else
        {
            hPhiMassResolutionFit[i]->Draw("pe same");
        }
        if (TotalFiles > 1)
            legMassRes->AddEntry(hPhiMassResolutionFit[i], LegendNames[i], "lpe");
        if (i == TotalFiles - 1)
        {
            legMassRes->Draw();
        }

        cCompareYield->cd();
        gPad->SetLogy();
        SetHistoQA(hPhiYieldFit[i]);
        hPhiYieldFit[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hPhiYieldFit[i]->GetYaxis()->SetTitle("Yield");
        hPhiYieldFit[i]->SetMarkerStyle(markerStyles[i]);
        hPhiYieldFit[i]->SetMarkerColor(vibrantColors[i]);
        hPhiYieldFit[i]->SetLineColor(vibrantColors[i]);
        hPhiYieldFit[i]->SetMarkerSize(1.3);
        hPhiYieldFit[i]->SetLineWidth(3);
        if (i == 0)
        {
            // hPhiYieldFit[i]->GetYaxis()->SetRangeUser(1, 10000);
            hPhiYieldFit[i]->SetMaximum(hPhiYieldFit[i]->GetMaximum() * 5);
            hPhiYieldFit[i]->Draw("pe");
        }
        else
        {
            hPhiYieldFit[i]->Draw("pe same");
        }
        if (i == TotalFiles - 1)
        {
            legMassRes->Draw();
        }

        cComparePurity->cd();
        SetHistoQA(hPurity[i]);
        hPurity[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hPurity[i]->GetYaxis()->SetTitle("Purity (%)");
        // hPurity[i]->SetMarkerStyle(markerStyles[i]);
        // hPurity[i]->SetMarkerColor(vibrantColors[i]);
        hPurity[i]->SetLineColor(vibrantColors[i]);
        hPurity[i]->SetLineStyle(lineStyles[i]);
        hPurity[i]->SetMarkerSize(1.3);
        hPurity[i]->SetLineWidth(3);
        if (i == 0)
        {
            hPurity[i]->GetYaxis()->SetRangeUser(0, 100);
            hPurity[i]->Draw("HIST");
        }
        else
        {
            hPurity[i]->Draw("HIST same");
        }
        if (TotalFiles > 1)
            legMassRes2->AddEntry(hPurity[i], LegendNames[i], "l");
        if (i == TotalFiles - 1)
        {
            legMassRes2->Draw();
        }

        cCompareNumPhi->cd();
        SetHistoQA(hNumPhi[i]);
        hNumPhi[i]->GetXaxis()->SetTitle("Number of #Phi mesons in an event");
        hNumPhi[i]->GetYaxis()->SetTitle("Counts");
        // hNumPhi[i]->SetMarkerStyle(markerStyles[i]);
        // hNumPhi[i]->SetMarkerColor(vibrantColors[i]);
        hNumPhi[i]->SetLineColor(vibrantColors[i]);
        hNumPhi[i]->SetLineStyle(lineStyles[i]);
        hNumPhi[i]->SetMarkerSize(1.3);
        hNumPhi[i]->SetLineWidth(3);
        if (i == 0)
        {
            hNumPhi[i]->GetXaxis()->SetRangeUser(0, 20.0);
            hNumPhi[i]->GetYaxis()->SetRangeUser(1, hNumPhi[i]->GetMaximum() * 15);
            gPad->SetLogy();
            hNumPhi[i]->Draw("HIST");
        }
        else
        {
            hNumPhi[i]->Draw("HIST same");
        }
        if (i == TotalFiles - 1)
        {
            legMassRes->Draw();
        }
        // */

        cComparePhiPhiInvMass->cd();
        SetHistoQA(hPhiPhiInvMass[i]);
        hPhiPhiInvMass[i]->GetXaxis()->SetTitle("#it{M}_{#Phi#Phi} (GeV/#it{c}^{2})");
        hPhiPhiInvMass[i]->GetYaxis()->SetTitle("Counts");
        hPhiPhiInvMass[i]->SetMarkerStyle(markerStyles[i]);
        hPhiPhiInvMass[i]->SetMarkerColor(vibrantColors[i]);
        hPhiPhiInvMass[i]->SetLineColor(vibrantColors[i]);
        hPhiPhiInvMass[i]->SetLineStyle(lineStyles[i]);
        hPhiPhiInvMass[i]->SetMarkerSize(1.3);
        hPhiPhiInvMass[i]->SetLineWidth(2);
        if (i == 0)
        {
            hPhiPhiInvMass[i]->GetXaxis()->SetRangeUser(2.68, 2.84);
            // hPhiPhiInvMass[i]->GetYaxis()->SetRangeUser(1, hPhiPhiInvMass[i]->GetMaximum() * 1.5);
            hPhiPhiInvMass[i]->Draw("pe");
        }
        else
        {
            hPhiPhiInvMass[i]->Draw("pe same");
        }
        if (TotalFiles > 1)
            legPhiPhiMass->AddEntry(hPhiPhiInvMass[i], LegendNames[i], "pe");
        if (i == TotalFiles - 1)
        {
            legPhiPhiMass->Draw();
        }
    }
    cCompareMass->SaveAs(outputFilePath + "Compare_Phi_Mass_Fit.png");
    cCompareMassResolution->SaveAs(outputFilePath + "Compare_Phi_Mass_Resolution_Fit.png");
    cCompareYield->SaveAs(outputFilePath + "Compare_Phi_Yield_Fit.png");
    cComparePurity->SaveAs(outputFilePath + "Compare_Phi_Purity.png");
    cCompareNumPhi->SaveAs(outputFilePath + "Compare_NumPhi.png");
    cComparePhiPhiInvMass->SaveAs(outputFilePath + "Compare_PhiPhi_InvMass.png");
}