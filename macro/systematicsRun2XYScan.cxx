#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <iostream>
#include <vector>

void systematicsRun2XYScan()
{
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    const int nBins = 14;
    double ptBins[nBins + 1] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 10.0};

    double fitting[nBins] = {1.9, 1.6, 1.1, 1.2, 1.8, 0.6, 0.3, 0.9, 1.1, 1.4, 2.9, 0.7, 1.6, 1.0};
    double normalization[nBins] = {1.1, 1.1, 0.1, 0.1, 0.1, 0.3, 0.5, 0.6, 1.3, 0.3, 1.6, 1.4, 0.3, 0.3};
    double resBkg[nBins] = {0.1, 0.5, 1.0, 1.7, 2.2, 3.1, 0.8, 1.2, 0.9, 0.6, 5.7, 5.7, 2.8, 0.2};
    double freeWidth[nBins] = {6.8, 4.6, 4.7, 5.1, 5.3, 4.9, 4.9, 4.9, 3.6, 5.3, 4.6, 3.7, 6.9, 5.3};
    double PID[nBins] = {7.1, 3.3, 3.2, 3.8, 3.6, 5.5, 2.9, 5.3, 1.9, 1.8, 1.6, 1.7, 1.8, 2.5};
    double trackCut[nBins] = {3.5, 2.0, 2.0, 1.5, 1.5, 1.8, 1.8, 1.8, 1.9, 2.0, 2.0, 2.4, 2.2, 2.0};
    double vz[nBins] = {1.8, 1.8, 1.2, 1.1, 1.0, 1.0, 1.0, 1.0, 1.9, 1.8, 0.7, 1.7, 1.8, 2.5};
    double tracking[nBins] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    double hadronic[nBins] = {2.7, 2.5, 2.0, 2.1, 2.1, 1.2, 1.0, 1.0, 1.3, 0.0, 0.0, 0.0, 0.0, 0.0};
    double material[nBins] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double total[nBins] = {11.9, 7.3, 6.8, 7.5, 7.5, 7.4, 6.3, 6.2, 5.2, 6.6, 8.6, 7.9, 8.5, 7.4};

    auto makeHist = [&](const char *name, const char *title, double *values, int color, int style, int width)
    {
        TH1D *h = new TH1D(name, title, nBins, ptBins);
        for (int i = 0; i < nBins; i++)
            h->SetBinContent(i + 1, values[i] / 100.0);
        h->SetLineColor(color);
        h->SetLineStyle(style);
        h->SetLineWidth(width);
        h->SetFillStyle(0);
        return h;
    };

    TH1D *hFitting = makeHist("hFitting", "Fitting range", fitting, kRed, 2, 3);
    TH1D *hNorm = makeHist("hNorm", "Normalization region", normalization, kBlue, 2, 3);
    TH1D *hResBkg = makeHist("hResBkg", "Res. Bkg.", resBkg, kViolet + 1, 2, 3);
    TH1D *hFreeWidth = makeHist("hFreeWidth", "Free Width", freeWidth, kOrange + 7, 3, 3);
    TH1D *hPID = makeHist("hPID", "PID", PID, kGreen + 2, 3, 4);
    TH1D *hTrack = makeHist("hTrack", "Track cut", trackCut, kBlack, 1, 1);
    TH1D *hVz = makeHist("hVz", "V_{z}", vz, kGreen + 3, 1, 1);
    TH1D *hTracking = makeHist("hTracking", "Tracking efficiency", tracking, kMagenta, 1, 1);
    TH1D *hHadronic = makeHist("hHadronic", "Hadronic interaction", hadronic, kBlue, 1, 1);
    TH1D *hMaterial = makeHist("hMaterial", "Material budget", material, kRed, 1, 1);
    TH1D *hTotal = makeHist("hTotal", "Total", total, kBlack, 3, 4);

    TCanvas *c = new TCanvas("c", "Systematic uncertainties", 1000, 700);
    c->SetLeftMargin(0.10);
    c->SetRightMargin(0.03);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.05);

    TH1D *hFrame = new TH1D("hFrame", "", nBins, ptBins);
    hFrame->SetMinimum(0);
    hFrame->SetMaximum(0.3);
    hFrame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hFrame->GetYaxis()->SetTitle("Relative uncertainty");
    hFrame->GetXaxis()->SetTitleSize(0.055);
    hFrame->GetYaxis()->SetTitleSize(0.055);
    hFrame->GetXaxis()->SetLabelSize(0.045);
    hFrame->GetYaxis()->SetLabelSize(0.045);
    hFrame->GetYaxis()->SetTitleOffset(0.85);
    hFrame->Draw();

    hFitting->Draw("HIST SAME");
    hNorm->Draw("HIST SAME");
    hResBkg->Draw("HIST SAME");
    hFreeWidth->Draw("HIST SAME");
    hPID->Draw("HIST SAME");
    hTrack->Draw("HIST SAME");
    hVz->Draw("HIST SAME");
    hTracking->Draw("HIST SAME");
    hHadronic->Draw("HIST SAME");
    hMaterial->Draw("HIST SAME");
    hTotal->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.73, 0.50, 0.98, 0.96);
    leg->SetTextSize(0.033);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hFitting, "Fitting range", "l");
    leg->AddEntry(hNorm, "Normalization region", "l");
    leg->AddEntry(hResBkg, "Res. Bkg.", "l");
    leg->AddEntry(hFreeWidth, "Free Width", "l");
    leg->AddEntry(hPID, "PID", "l");
    leg->AddEntry(hTrack, "Track cut", "l");
    leg->AddEntry(hVz, "V_{z}", "l");
    leg->AddEntry(hTracking, "Tracking efficiency", "l");
    leg->AddEntry(hHadronic, "Hadronic interaction", "l");
    leg->AddEntry(hMaterial, "Material budget", "l");
    leg->AddEntry(hTotal, "Total", "l");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.DrawLatex(0.77, 0.43, "multiplicity 0-100%");
    latex.DrawLatex(0.77, 0.39, "K^{*0}, pp 13 TeV");

    // ============================================================
    // Print systematic uncertainties vs pT
    // ============================================================

    string shortNames[] = {"Fit", "Norm", "ResBkg", "Width", "PID", "TrkSel", "Vz", "TrkEff", "HI", "MB", "Total"};

    double *values[] = {
        fitting,
        normalization,
        resBkg,
        freeWidth,
        PID,
        trackCut,
        vz,
        tracking,
        hadronic,
        material,
        total};

    const int nSources = 11;

    std::cout << "\n\n=================================================================================================\n";
    std::cout << "             Systematic uncertainties (0-100%) in Run 2 data (pp 13 TeV)\n";
    std::cout << "=================================================================================================\n";

    std::cout << std::left << std::setw(9) << "pT";

    for (int j = 0; j < nSources; j++)
        std::cout << std::right << std::setw(8) << shortNames[j];

    std::cout << "\n";

    std::cout << std::string(9 + 8 * nSources, '-') << "\n";

    // pT-bin values
    for (int i = 0; i < nBins; i++)
    {
        std::string ptRange = Form("%.1f-%.1f", ptBins[i], ptBins[i + 1]);

        std::cout << std::left << std::setw(9) << ptRange;

        for (int j = 0; j < nSources; j++)
        {
            std::cout << std::right << std::setw(8) << std::fixed << std::setprecision(1) << values[j][i];
        }

        std::cout << "\n";
    }

    // Average
    std::cout << std::string(9 + 8 * nSources, '-') << "\n";

    std::cout << std::left << std::setw(9) << "Avg.";

    for (int j = 0; j < nSources; j++)
    {
        double sum = 0.0;

        for (int i = 0; i < nBins; i++)
            sum += values[j][i];

        double average = sum / nBins;

        std::cout << std::right << std::setw(8) << std::fixed << std::setprecision(1) << average;
    }

    std::cout << "\n";
    std::cout << "=================================================================================================\n";
}