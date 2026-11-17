#include <iostream>
#include "../src/style.h"

double pTBins[] = {0.0, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};

void PlotMassShift(TH2 *h25ac, TH2 *h25ai, TH2 *h26, bool momentumOnX = true, TString particle = "Phi", TString whichShift = "Mass", TString xAxis = "#it{p}_{T} (GeV/#it{c})", double maxPlotRange = 15.0);

void momentumShiftMC()
{
    gStyle->SetOptStat(0);

    TFile *file25ac = TFile::Open("/home/sawan/Storage/check_k892/mc/PhiPhi/PhiMomentumShift/AnalysisResults_LHC25ac.root");
    TFile *file25ai = TFile::Open("/home/sawan/Storage/check_k892/mc/PhiPhi/PhiMomentumShift/AnalysisResults_LHC25ai.root");
    TFile *file26 = TFile::Open("/home/sawan/Storage/check_k892/mc/PhiPhi/PhiMomentumShift/AnalysisResults_LHC26.root");

    if (!file25ac || !file25ai || !file26 || file25ac->IsZombie() || file25ai->IsZombie() || file26->IsZombie())
    {
        std::cerr << "Error opening ROOT files!" << std::endl;
        return;
    }

    TH1D *hGenMass25ac = (TH1D *)file25ac->Get("kstarqa/hShift/h1D_MassPhiGen");
    TH1D *hGenMass25ai = (TH1D *)file25ai->Get("kstarqa/hShift/h1D_MassPhiGen");
    TH1D *hGenMass26 = (TH1D *)file26->Get("kstarqa/hShift/h1D_MassPhiGen");

    TH1D *hRecMass25ac = (TH1D *)file25ac->Get("kstarqa/hShift/h1D_MassPhiRec");
    TH1D *hRecMass25ai = (TH1D *)file25ai->Get("kstarqa/hShift/h1D_MassPhiRec");
    TH1D *hRecMass26 = (TH1D *)file26->Get("kstarqa/hShift/h1D_MassPhiRec");

    TCanvas *cGenMassPhi = new TCanvas("cGenMassPhi", "Generated Mass Distribution", 720, 720);
    SetCanvasStyle(cGenMassPhi, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hGenMass25ac);
    SetHistoQA(hGenMass25ai);
    SetHistoQA(hGenMass26);
    hGenMass25ac->Scale(1.0 / hGenMass25ac->Integral());
    hGenMass25ai->Scale(1.0 / hGenMass25ai->Integral());
    hGenMass26->Scale(1.0 / hGenMass26->Integral());
    hGenMass25ac->SetLineColor(kBlue);
    hGenMass25ai->SetLineColor(kRed);
    hGenMass26->SetLineColor(kGreen + 2);
    hGenMass25ac->SetMarkerColor(kBlue);
    hGenMass25ai->SetMarkerColor(kRed);
    hGenMass26->SetMarkerColor(kGreen + 2);
    hGenMass25ac->SetLineWidth(2);
    hGenMass25ai->SetLineWidth(2);
    hGenMass26->SetLineWidth(2);
    hGenMass25ac->GetXaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
    hGenMass25ac->GetYaxis()->SetTitle("Counts");
    hGenMass25ac->Draw("pe");
    hGenMass25ai->Draw("pe same");
    hGenMass26->Draw("pe same");

    TLegend *legendGenMass = new TLegend(0.5, 0.7, 0.88, 0.88);
    legendGenMass->SetBorderSize(0);
    legendGenMass->SetFillStyle(0);
    legendGenMass->SetTextFont(42);
    legendGenMass->SetTextSize(0.025);
    legendGenMass->AddEntry(hGenMass25ac, Form("LHC25ac: #it{M}_{#phi} = %.2f MeV/#it{c}^{2}", hGenMass25ac->GetMean() * 1000), "pe");
    legendGenMass->AddEntry(hGenMass25ai, Form("LHC25ai: #it{M}_{#phi} = %.2f MeV/#it{c}^{2}", hGenMass25ai->GetMean() * 1000), "pe");
    legendGenMass->AddEntry(hGenMass26, Form("LHC26: #it{M}_{#phi} = %.2f MeV/#it{c}^{2}", hGenMass26->GetMean() * 1000), "pe");
    legendGenMass->Draw();
    cGenMassPhi->SaveAs("MomentumShift/GeneratedMassPhi.png");

    TCanvas *cRecMassPhi = new TCanvas("cRecMassPhi", "Reconstructed Mass Distribution", 720, 720);
    SetCanvasStyle(cRecMassPhi, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hRecMass25ac);
    SetHistoQA(hRecMass25ai);
    SetHistoQA(hRecMass26);
    hRecMass25ac->Scale(1.0 / hRecMass25ac->Integral());
    hRecMass25ai->Scale(1.0 / hRecMass25ai->Integral());
    hRecMass26->Scale(1.0 / hRecMass26->Integral());
    hRecMass25ac->SetLineColor(kBlue);
    hRecMass25ai->SetLineColor(kRed);
    hRecMass26->SetLineColor(kGreen + 2);
    hRecMass25ac->SetMarkerColor(kBlue);
    hRecMass25ai->SetMarkerColor(kRed);
    hRecMass26->SetMarkerColor(kGreen + 2);
    hRecMass25ac->SetLineWidth(2);
    hRecMass25ai->SetLineWidth(2);
    hRecMass26->SetLineWidth(2);
    hRecMass25ac->GetXaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
    hRecMass25ac->GetYaxis()->SetTitle("Counts");
    hRecMass25ac->Draw("pe");
    hRecMass25ai->Draw("pe same");
    hRecMass26->Draw("pe same");

    TLegend *legendRecMass = new TLegend(0.5, 0.7, 0.88, 0.88);
    legendRecMass->SetBorderSize(0);
    legendRecMass->SetFillStyle(0);
    legendRecMass->SetTextFont(42);
    legendRecMass->SetTextSize(0.025);
    legendRecMass->AddEntry(hRecMass25ac, Form("LHC25ac: #it{M}_{#phi} = %.2f MeV/#it{c}^{2}", hRecMass25ac->GetMean() * 1000), "pe");
    legendRecMass->AddEntry(hRecMass25ai, Form("LHC25ai: #it{M}_{#phi} = %.2f MeV/#it{c}^{2}", hRecMass25ai->GetMean() * 1000), "pe");
    legendRecMass->AddEntry(hRecMass26, Form("LHC26: #it{M}_{#phi} = %.2f MeV/#it{c}^{2}", hRecMass26->GetMean() * 1000), "pe");
    legendRecMass->Draw();
    cRecMassPhi->SaveAs("MomentumShift/ReconstructedMassPhi.png");

    // Phi reconstructed mass vs pT (without fitting)
    TH2D *hRecMassVsPt25ac = (TH2D *)file25ac->Get("kstarqa/hShift/h2D_pt_massShiftPhi");
    TH2D *hRecMassVsPt25ai = (TH2D *)file25ai->Get("kstarqa/hShift/h2D_pt_massShiftPhi");
    TH2D *hRecMassVsPt26 = (TH2D *)file26->Get("kstarqa/hShift/h2D_pt_massShiftPhi");

    int nBins = sizeof(pTBins) / sizeof(pTBins[0]) - 1;
    TH1D *hRecMassVspT25ac = new TH1D("hRecMassVspT25ac", "Reconstructed Mass vs pT (LHC25ac)", nBins, pTBins);
    TH1D *hRecMassVspT25ai = new TH1D("hRecMassVspT25ai", "Reconstructed Mass vs pT (LHC25ai)", nBins, pTBins);
    TH1D *hRecMassVspT26 = new TH1D("hRecMassVspT26", "Reconstructed Mass vs pT (LHC26)", nBins, pTBins);

    for (int i = 0; i < nBins; ++i)
    {
        double lowpT = pTBins[i];
        double highpT = pTBins[i + 1];
        int lowpTbin = hRecMassVsPt25ac->GetXaxis()->FindBin(lowpT + 1e-6);   // Add a small value to avoid binning issues
        int highpTbin = hRecMassVsPt25ac->GetXaxis()->FindBin(highpT - 1e-6); // Subtract a small value to avoid

        TH1D *hMassSlice25ac = hRecMassVsPt25ac->ProjectionY(Form("hMassSlice25ac_%.1f_%.1f", lowpT, highpT), lowpTbin, highpTbin);
        TH1D *hMassSlice25ai = hRecMassVsPt25ai->ProjectionY(Form("hMassSlice25ai_%.1f_%.1f", lowpT, highpT), lowpTbin, highpTbin);
        TH1D *hMassSlice26 = hRecMassVsPt26->ProjectionY(Form("hMassSlice26_%.1f_%.1f", lowpT, highpT), lowpTbin, highpTbin);

        double meanShift25ac = hMassSlice25ac->GetMean();
        double meanShift25ai = hMassSlice25ai->GetMean();
        double meanShift26 = hMassSlice26->GetMean();

        double recMass25ac = (hGenMass25ac->GetMean() - meanShift25ac) * 1000;
        double recMass25ai = (hGenMass25ai->GetMean() - meanShift25ai) * 1000;
        double recMass26 = (hGenMass26->GetMean() - meanShift26) * 1000;

        cout << "pT bin: [" << lowpT << ", " << highpT << "] GeV/c" << std::endl;
        cout << "LHC25ac: Reconstructed Mass = " << recMass25ac << " MeV/c^2" << std::endl;
        cout << "LHC25ai: Reconstructed Mass = " << recMass25ai << " MeV/c^2" << std::endl;
        cout << "LHC26: Reconstructed Mass = " << recMass26 << " MeV/c^2" << std::endl;
        cout << "----------------------------------------" << std::endl;

        hRecMassVspT25ac->SetBinContent(i + 1, recMass25ac);
        hRecMassVspT25ai->SetBinContent(i + 1, recMass25ai);
        hRecMassVspT26->SetBinContent(i + 1, recMass26);

        hRecMassVspT25ac->SetBinError(i + 1, (hMassSlice25ac->GetMeanError() - meanShift25ac) * 1000);
        hRecMassVspT25ai->SetBinError(i + 1, (hMassSlice25ai->GetMeanError() - meanShift25ai) * 1000);
        hRecMassVspT26->SetBinError(i + 1, (hMassSlice26->GetMeanError() - meanShift26) * 1000);
    }

    TCanvas *cRecMassVspT = new TCanvas("cRecMassVspT", "Reconstructed Mass vs pT", 720, 720);
    SetCanvasStyle(cRecMassVspT, 0.19, 0.06, 0.06, 0.13);
    SetHistoQA(hRecMassVspT25ac);
    SetHistoQA(hRecMassVspT25ai);
    SetHistoQA(hRecMassVspT26);
    hRecMassVspT25ac->SetLineColor(kBlue);
    hRecMassVspT25ai->SetLineColor(kRed);
    hRecMassVspT26->SetLineColor(kGreen + 2);
    hRecMassVspT25ac->SetMarkerColor(kBlue);
    hRecMassVspT25ai->SetMarkerColor(kRed);
    hRecMassVspT26->SetMarkerColor(kGreen + 2);
    hRecMassVspT25ac->SetLineWidth(2);
    hRecMassVspT25ai->SetLineWidth(2);
    hRecMassVspT26->SetLineWidth(2);
    hRecMassVspT25ac->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRecMassVspT25ac->GetYaxis()->SetTitle("Reconstructed Mass (MeV/#it{c}^{2})");
    hRecMassVspT25ac->SetMaximum(hRecMassVspT25ac->GetMaximum() * 1.0015); // Add some space above the highest bin
    hRecMassVspT25ac->SetMinimum(hRecMassVspT25ac->GetMinimum() * 0.999);  // Add some space below the lowest bin
    hRecMassVspT25ac->GetYaxis()->SetMaxDigits(3);                         // Set maximum digits on Y-axis to 3
    hRecMassVspT25ac->GetYaxis()->SetTitleOffset(2.1);                     // Adjust Y-axis title offset
    hRecMassVspT25ac->Draw("pe");
    hRecMassVspT25ai->Draw("pe same");
    hRecMassVspT26->Draw("pe same");

    TLine *lineGen = new TLine(0, 1019.6, 20, 1019.6);
    lineGen->SetLineColor(kBlack);
    lineGen->SetLineStyle(2);
    lineGen->SetLineWidth(2);
    lineGen->Draw("same");

    legendRecMass->AddEntry(lineGen, "Gen. Mass #it{M}_{#phi} = 1019.6 MeV/#it{c}^{2}", "l");
    legendRecMass->Draw();
    cRecMassVspT->SaveAs("MomentumShift/ReconstructedMassVsPt.png");

    // Phi mass shift: pT on X-axis
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_pt_massShiftPhi"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_pt_massShiftPhi"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_pt_massShiftPhi"),
        true, "Phi", "Mass", "#it{p}_{T} (GeV/#it{c})", 15.0);

    // Phi pT shift
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpt_vs_pt_Phi"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpt_vs_pt_Phi"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpt_vs_pt_Phi"),
        false, "Phi", "p_{T}", "p_{T} (GeV/#it{c})", 10.0);

    // Phi px shift
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpx_vs_px_Phi"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpx_vs_px_Phi"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpx_vs_px_Phi"),
        false, "Phi", "p_{x}", "p_{x} (GeV/#it{c})", 10.0);

    // Phi py shift
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpy_vs_py_Phi"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpy_vs_py_Phi"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpy_vs_py_Phi"),
        false, "Phi", "p_{y}", "p_{y} (GeV/#it{c})", 10.0);

    // Phi pz shift
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpz_vs_pz_Phi"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpz_vs_pz_Phi"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpz_vs_pz_Phi"),
        false, "Phi", "p_{z}", "p_{z} (GeV/#it{c})", 10.0);

    // Kaon pT shift (All kaons)
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpt_vs_pt_AllKaons"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpt_vs_pt_AllKaons"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpt_vs_pt_AllKaons"),
        false, "Kaon", "p_{T}", "p_{T} (GeV/#it{c})", 10.0);

    // Kaon px shift (All kaons)
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpx_vs_px_AllKaons"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpx_vs_px_AllKaons"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpx_vs_px_AllKaons"),
        false, "Kaon", "p_{x}", "p_{x} (GeV/#it{c})", 10.0);

    // Kaon py shift (All kaons)
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpy_vs_py_AllKaons"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpy_vs_py_AllKaons"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpy_vs_py_AllKaons"),
        false, "Kaon", "p_{y}", "p_{y} (GeV/#it{c})", 10.0);

    // Kaon pz shift (All kaons)
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpz_vs_pz_AllKaons"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpz_vs_pz_AllKaons"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpz_vs_pz_AllKaons"),
        false, "Kaon", "p_{z}", "p_{z} (GeV/#it{c})", 10.0);

    // kaon pT shift (Kaons from Phi)
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpt_vs_pt_KaonFromPhi"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpt_vs_pt_KaonFromPhi"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpt_vs_pt_KaonFromPhi"),
        false, "Kaon (#phi daughters)", "p_{T}", "p_{T} (GeV/#it{c})", 10.0);

    // kaon px shift (Kaons from Phi)
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpx_vs_px_KaonFromPhi"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpx_vs_px_KaonFromPhi"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpx_vs_px_KaonFromPhi"),
        false, "Kaon (#phi daughters)", "p_{x}", "p_{x} (GeV/#it{c})", 10.0);

    // kaon py shift (Kaons from Phi)
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpy_vs_py_KaonFromPhi"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpy_vs_py_KaonFromPhi"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpy_vs_py_KaonFromPhi"),
        false, "Kaon (#phi daughters)", "p_{y}", "p_{y} (GeV/#it{c})", 10.0);

    // kaon pz shift (Kaons from Phi)
    PlotMassShift(
        (TH2 *)file25ac->Get("kstarqa/hShift/h2D_dpz_vs_pz_KaonFromPhi"),
        (TH2 *)file25ai->Get("kstarqa/hShift/h2D_dpz_vs_pz_KaonFromPhi"),
        (TH2 *)file26->Get("kstarqa/hShift/h2D_dpz_vs_pz_KaonFromPhi"),
        false, "Kaon (#phi daughters)", "p_{z}", "p_{z} (GeV/#it{c})", 10.0);
}

void PlotMassShift(TH2 *h25ac, TH2 *h25ai, TH2 *h26, bool momentumOnX = true, TString particle = "Phi", TString whichShift = "Mass", TString xLabel = "#it{p}_{T} (GeV/#it{c})", double maxPlotRange = 15.0)
{
    if (!h25ac || !h25ai || !h26 || h25ac == nullptr || h25ai == nullptr || h26 == nullptr)
    {
        std::cerr << "Error: One or more histograms are null!" << std::endl;
        return;
    }

    int totalBins = sizeof(pTBins) / sizeof(pTBins[0]) - 1;

    TH2 *h2D[3] = {h25ac, h25ai, h26};
    TString labels[3] = {"LHC25ac", "LHC25ai", "LHC26"};
    int colors[3] = {kBlue, kRed, kGreen + 2};
    int markers[3] = {20, 21, 22};

    TH1D *hMean[3];

    TCanvas *cMassShift = new TCanvas(Form("cMassShift_%s_%s", particle.Data(), whichShift.Data()), Form("%s %s Shift", particle.Data(), whichShift.Data()), 1280, 720);
    cMassShift->Divide(4, 4);

    for (int j = 0; j < 3; j++)
    {
        hMean[j] = new TH1D(Form("hMeanMassShift_%s_%s", particle.Data(), labels[j].Data()), "", totalBins, pTBins);
        hMean[j]->SetDirectory(nullptr);
        SetHistoQA(hMean[j]);
        hMean[j]->SetLineColor(colors[j]);
        hMean[j]->SetMarkerColor(colors[j]);
        hMean[j]->SetMarkerStyle(markers[j]);
        hMean[j]->SetMarkerSize(1.3);
    }

    for (int i = 0; i < totalBins; i++)
    {
        double lowPt = pTBins[i];
        double highPt = pTBins[i + 1];

        TH1D *hShift[3] = {nullptr, nullptr, nullptr};

        for (int j = 0; j < 3; j++)
        {
            TAxis *pTAxis = momentumOnX ? h2D[j]->GetXaxis() : h2D[j]->GetYaxis();

            int lowBin = pTAxis->FindFixBin(lowPt + 1e-6);
            int highBin = pTAxis->FindFixBin(highPt - 1e-6);

            lowBin = std::max(1, lowBin);
            highBin = std::min(pTAxis->GetNbins(), highBin);

            if (lowBin > highBin)
                continue;

            TString histName = Form("hShift_%s_%s_bin%d", particle.Data(), labels[j].Data(), i);

            if (momentumOnX)
                hShift[j] = h2D[j]->ProjectionY(histName, lowBin, highBin, "e");
            else
                hShift[j] = h2D[j]->ProjectionX(histName, lowBin, highBin, "e");

            hShift[j]->SetDirectory(nullptr);

            // Calculate mean and statistical error before normalization
            double mean = hShift[j]->GetMean() * 1000;           // Convert to MeV/c^2
            double meanError = hShift[j]->GetMeanError() * 1000; // Convert to MeV/c^2

            hMean[j]->SetBinContent(i + 1, mean);
            hMean[j]->SetBinError(i + 1, meanError);

            double integral = hShift[j]->Integral();
            if (integral > 0)
                hShift[j]->Scale(1.0 / integral);

            SetHistoQA(hShift[j]);
            hShift[j]->SetLineColor(colors[j]);
            hShift[j]->SetMarkerColor(colors[j]);
            hShift[j]->SetLineWidth(2);
        }

        cMassShift->cd(i + 1);
        gPad->SetTicks(1, 1);

        double maxY = 0;
        for (int j = 0; j < 3; j++)
        {
            if (hShift[j])
                maxY = std::max(maxY, hShift[j]->GetMaximum()); // for plotting all three histograms on the same canvas
        }

        bool firstDraw = true;

        for (int j = 0; j < 3; j++)
        {
            if (!hShift[j])
                continue;

            hShift[j]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", lowPt, highPt));
            hShift[j]->GetXaxis()->SetTitle(xLabel);
            hShift[j]->GetYaxis()->SetTitle("Normalized counts");
            hShift[j]->GetXaxis()->SetRangeUser(-0.02, 0.02);
            hShift[j]->SetMaximum(maxY * 1.5);
            hShift[j]->SetMinimum(hShift[j]->GetMinimum() * 0.8);

            if (firstDraw)
            {
                hShift[j]->Draw("hist");
                firstDraw = false;
            }
            else
                hShift[j]->Draw("hist same");
        }
    }

    TCanvas *cMean = new TCanvas(Form("cMeanMassShift_%s_%s", particle.Data(), whichShift.Data()), Form("Mean %s %s Shift", particle.Data(), whichShift.Data()), 720, 720);
    SetCanvasStyle(cMean, 0.16, 0.06, 0.06, 0.13);

    hMean[2]->GetXaxis()->SetTitle(xLabel);
    hMean[2]->GetYaxis()->SetTitle(Form("Mean %s shift (MeV/#it{c}^{2})", whichShift.Data()));
    hMean[2]->GetXaxis()->SetRangeUser(0, maxPlotRange);

    double maxMeanY = 0;
    for (int j = 0; j < 3; j++)
        maxMeanY = std::max(maxMeanY, hMean[j]->GetMaximum());

    // hMean[2]->SetMaximum(maxMeanY * 0.6);
    hMean[2]->SetMaximum(hMean[2]->GetMaximum() * 1.9); // Set maximum based on the last histogram's maximum
    hMean[2]->Draw("pe");
    hMean[0]->Draw("pe same");
    hMean[1]->Draw("pe same");

    TLegend *legend = new TLegend(0.55, 0.7, 0.92, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.022);
    legend->SetHeader(Form("%s %s shift (Gen - Rec)", particle.Data(), whichShift.Data()));

    for (int j = 0; j < 3; j++)
        legend->AddEntry(hMean[j], labels[j], "lep");

    legend->Draw();

    cMassShift->Update();
    cMean->Update();

    cMassShift->SaveAs(Form("MomentumShift/GaussianDistributions/MomentumShift_%s_%s.png", particle.Data(), whichShift.Data()));
    cMean->SaveAs(Form("MomentumShift/MeanShifts/MeanMomentumShift_%s_%s.png", particle.Data(), whichShift.Data()));
}
