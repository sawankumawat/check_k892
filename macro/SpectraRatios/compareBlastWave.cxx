#include "../src/style.h"

void compareBlastWave()
{
    TFile *fWithoutKstar = new TFile("BlastWaveFits/BlastWaveFitResults.root", "read");
    TFile *fWithKstar = new TFile("BlastWaveFits/BlastWaveFitResultsKstar.root", "read");
    if (!fWithoutKstar || !fWithKstar)
    {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }
    TGraphErrors *grBetaAvgWithoutKstar = (TGraphErrors *)fWithoutKstar->Get("grBetaAvg");
    TGraphErrors *grBetaAvgWithKstar = (TGraphErrors *)fWithKstar->Get("grBetaAvg");
    TGraphErrors *grTWithoutKstar = (TGraphErrors *)fWithoutKstar->Get("grT");
    TGraphErrors *grTWithKstar = (TGraphErrors *)fWithKstar->Get("grT");
    TGraphErrors *grNWithoutKstar = (TGraphErrors *)fWithoutKstar->Get("grN");
    TGraphErrors *grNWithKstar = (TGraphErrors *)fWithKstar->Get("grN");
    TGraphErrors *grChi2WithoutKstar = (TGraphErrors *)fWithoutKstar->Get("grChi2");
    TGraphErrors *grChi2WithKstar = (TGraphErrors *)fWithKstar->Get("grChi2");

    TCanvas *cAvgBeta = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cAvgBeta, 0.15, 0.03, 0.03, 0.14);
    SetGraphErrorStyle(grBetaAvgWithoutKstar);
    grBetaAvgWithoutKstar->SetMarkerColor(kBlue + 1);
    grBetaAvgWithoutKstar->SetLineColor(kBlue + 1);
    grBetaAvgWithoutKstar->SetMarkerStyle(20);
    grBetaAvgWithoutKstar->SetMarkerSize(1.2);
    grBetaAvgWithoutKstar->GetXaxis()->SetTitle("#it{d}#it{N}_{ch}/d#it{#eta}");
    grBetaAvgWithoutKstar->GetYaxis()->SetTitle("<#it{#beta_{T}}>");
    grBetaAvgWithoutKstar->GetXaxis()->SetTitleOffset(1.2);
    grBetaAvgWithoutKstar->SetMinimum(0.1);
    grBetaAvgWithoutKstar->SetMaximum(0.6);
    grBetaAvgWithoutKstar->Draw("AP");  
    grBetaAvgWithKstar->SetMarkerColor(kRed + 1);
    grBetaAvgWithKstar->SetLineColor(kRed + 1);
    grBetaAvgWithKstar->SetMarkerStyle(21);
    grBetaAvgWithKstar->SetMarkerSize(1.2);
    grBetaAvgWithKstar->Draw("P same");
    TLegend *legend = new TLegend(0.35, 0.78, 0.68, 0.93);
    SetLegendStyle(legend);
    legend->SetTextSize(0.035);
    legend->SetHeader("Simultaneous Blast-Wave Fit");
    legend->AddEntry(grBetaAvgWithoutKstar, "#pi,K,p", "lp");
    legend->AddEntry(grBetaAvgWithKstar, "#pi,K,p,K^{*0}", "lp");
    legend->Draw();
    cAvgBeta->SaveAs("BlastWaveFits/AvgBetaComparison.png");

    TCanvas *cT = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cT, 0.15, 0.03, 0.03, 0.14);
    SetGraphErrorStyle(grTWithoutKstar);
    grTWithoutKstar->SetMarkerColor(kBlue + 1);
    grTWithoutKstar->SetLineColor(kBlue + 1);
    grTWithoutKstar->SetMarkerStyle(20);
    grTWithoutKstar->SetMarkerSize(1.2);
    grTWithoutKstar->GetXaxis()->SetTitle("#it{d}#it{N}_{ch}/d#it{#eta}");
    grTWithoutKstar->GetYaxis()->SetTitle("T_{kin} (GeV)");
    grTWithoutKstar->GetXaxis()->SetTitleOffset(1.2);
    grTWithoutKstar->SetMinimum(0.12);
    grTWithoutKstar->SetMaximum(0.18);
    grTWithoutKstar->Draw("AP");
    grTWithKstar->SetMarkerColor(kRed + 1);
    grTWithKstar->SetLineColor(kRed + 1);
    grTWithKstar->SetMarkerStyle(21);
    grTWithKstar->SetMarkerSize(1.2);
    grTWithKstar->Draw("P same");
    legend->Draw();
    cT->SaveAs("BlastWaveFits/TComparison.png");

    TCanvas *cN = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cN, 0.15, 0.03, 0.03, 0.14);
    SetGraphErrorStyle(grNWithoutKstar);
    grNWithoutKstar->SetMarkerColor(kBlue + 1);
    grNWithoutKstar->SetLineColor(kBlue + 1);
    grNWithoutKstar->SetMarkerStyle(20);
    grNWithoutKstar->SetMarkerSize(1.2);
    grNWithoutKstar->GetXaxis()->SetTitle("#it{d}#it{N}_{ch}/d#it{#eta}");
    grNWithoutKstar->GetYaxis()->SetTitle("n");
    grNWithoutKstar->GetXaxis()->SetTitleOffset(1.2);
    grNWithoutKstar->SetMinimum(1.0);
    grNWithoutKstar->SetMaximum(6.0);
    grNWithoutKstar->Draw("AP");
    grNWithKstar->SetMarkerColor(kRed + 1);
    grNWithKstar->SetLineColor(kRed + 1);
    grNWithKstar->SetMarkerStyle(21);
    grNWithKstar->SetMarkerSize(1.2);
    grNWithKstar->Draw("P same");
    legend->Draw();
    cN->SaveAs("BlastWaveFits/NComparison.png");

    TCanvas *cChi2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cChi2, 0.15, 0.03, 0.03, 0.14);
    SetGraphErrorStyle(grChi2WithoutKstar);
    gPad->SetLogy();
    grChi2WithoutKstar->SetMarkerColor(kBlue + 1);
    grChi2WithoutKstar->SetLineColor(kBlue + 1);
    grChi2WithoutKstar->SetMarkerStyle(20);
    grChi2WithoutKstar->SetMarkerSize(1.2);
    grChi2WithoutKstar->GetXaxis()->SetTitle("#it{d}#it{N}_{ch}/d#it{#eta}");
    grChi2WithoutKstar->GetYaxis()->SetTitle("#chi^{2}/NDF");
    grChi2WithoutKstar->GetXaxis()->SetTitleOffset(1.2);
    grChi2WithoutKstar->SetMinimum(1000.0);
    grChi2WithoutKstar->SetMaximum(100000);
    grChi2WithoutKstar->Draw("AP");
    grChi2WithKstar->SetMarkerColor(kRed + 1);
    grChi2WithKstar->SetLineColor(kRed + 1);
    grChi2WithKstar->SetMarkerStyle(21);
    grChi2WithKstar->SetMarkerSize(1.2);
    grChi2WithKstar->Draw("P same");
    legend->Draw();
    cChi2->SaveAs("BlastWaveFits/Chi2Comparison.png");
}