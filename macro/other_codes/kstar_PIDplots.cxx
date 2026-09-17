
#include <iostream>
#include "src/style.h"
// #include "src/fitfunc.h"
// #include "src/initializations.h"

void kstar_PIDplots()
{
    gStyle->SetOptStat(0);
    string filepath = "/home/sawan/check_k892/data/kstar/LHC22o_pass7";
    // string filepath = "/home/sawan/alice/practice";
    // TFile *file = new TFile((filepath + "/NewPID.root").c_str());
    TFile *file = new TFile((filepath + "/586469.root").c_str());
    if (file->IsZombie())
    {
        cerr << "File not found " << endl;
        return;
    }
    TH3F *h3PIDTPCPi = (TH3F *)file->Get("kstarqa/hPID/Before/hTPCnsigPi_mult_pt");
    TH3F *h3PIDTPCKa = (TH3F *)file->Get("kstarqa/hPID/Before/hTPCnsigKa_mult_pt");
    TH3F *h3PIDTPCTOFKa = (TH3F *)file->Get("kstarqa/hPID/Before/hNsigma_TPC_TOF_Ka_before");
    TH3F *h3PIDTPCTOFPi = (TH3F *)file->Get("kstarqa/hPID/Before/hNsigma_TPC_TOF_Pi_before");
    if (h3PIDTPCPi == nullptr || h3PIDTPCKa == nullptr)
    {
        cerr << "PID histograms before selection not found!!!!!!!!!!!!" << endl;
        return;
    }

    // TH2F *h2PIDTPCPiCut = (TH2F *)file->Get("kstarqa/hPID/AdditionalChecks/hTPCnsigPi_pt");
    // TH2F *h2PIDTPCKaCut = (TH2F *)file->Get("kstarqa/hPID/AdditionalChecks/hTPCnsigKa_pt");
    // TH2F *h2PIDTPCPiCutPr = (TH2F *)file->Get("kstarqa/hPID/AdditionalChecks/hTOFnsigKa_ptPrCut");
    // TH2F *h2PIDTPCKaCutPr = (TH2F *)file->Get("kstarqa/hPID/AdditionalChecks/hTOFnsigPi_ptPrCut");
    // if (h2PIDTPCPiCut == nullptr || h2PIDTPCKaCut == nullptr)
    // {
    //     cerr << "PID histograms after selection not found!!!!!!!!!!!!" << endl;
    //     return;
    // }

    TCanvas *cPIDTPCPi = new TCanvas("cPIDTPCPi", "PIDTPCPi", 1440, 720);
    SetCanvasStyle(cPIDTPCPi, 0.1, 0.05, 0.06, 0.17);
    cPIDTPCPi->Divide(5, 4, 0, 0);

    TCanvas *cPIDTPCKa = new TCanvas("cPIDTPCKa", "PIDTPCKa", 1440, 720);
    SetCanvasStyle(cPIDTPCKa, 0.1, 0.05, 0.06, 0.17);
    cPIDTPCKa->Divide(5, 4, 0, 0);

    TCanvas *cPIDTPCTOFKa = new TCanvas("cPIDTPCTOFKa", "PIDTPCTOFKa", 1440, 720);
    SetCanvasStyle(cPIDTPCTOFKa, 0.1, 0.05, 0.06, 0.17);
    cPIDTPCTOFKa->Divide(5, 4, 0, 0);

    TCanvas *cPIDTPCTOFPi = new TCanvas("cPIDTPCTOFPi", "PIDTPCTOFPi", 1440, 720);
    SetCanvasStyle(cPIDTPCTOFPi, 0.1, 0.05, 0.06, 0.17);
    cPIDTPCTOFPi->Divide(5, 4, 0, 0);

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.09);

    float pT_bins[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0};
    int total_pT_bins = sizeof(pT_bins) / sizeof(pT_bins[0]) - 1;
    float pion_contamination_peak_mean[] = {-7, -7, -6, -5, -3, -2, -1, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    float kaon_contamination_peak_mean[] = {6, 6, 6, 6, 6, 6, 6, 6, 5, -2, -2.5, -2.5, -2.5, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8};
    vector<float> purityPionVec, purityKaonVec, purityPionVecSelected, purityKaonVecSelected;
    vector<vector<float>> pionPIDrange = {{-2, 2}, {-2, 2}, {-2, 2}, {-2, 2}, {-2, 2}, {-2, 2}, {-2, 2}, {-2, 2}, {-1.5, 2}, {-1.5, 2}, {-1, 2}, {-1, 2}, {-1, 2}, {-1, 2}, {-1, 2}, {-0.5, 2}, {-0.5, 2}, {-0.5, 2}, {-0.5, 2}, {-0.5, 2}};
    vector<vector<float>> kaonPIDrange = {{-2, 2}, {-2, 2}, {-2, 2}, {-1.5, 2}, {-1.0, 2}, {-0.5, 2}, {0, 2}, {-2, 0}, {-2, 0}, {-2, 0}, {-2, 1}, {-2, 1}, {-2, 1}, {-2, 1}, {-2, 1.0}, {-2, 1.0}, {-2, 1.0}, {-2, 1.0}, {-2, 1.0}, {-2, 1.0}};

    for (Int_t ip = 0; ip < total_pT_bins; ip++) // start pt bin loop
    {
        int binpTlow = h3PIDTPCPi->GetZaxis()->FindBin(pT_bins[ip] + 1e-3);
        int binpThigh = h3PIDTPCPi->GetZaxis()->FindBin(pT_bins[ip + 1] - 1e-3);
        TH1F *hNsigmaTPCKaon = (TH1F *)h3PIDTPCKa->ProjectionX(Form("hNsigmaTPCKaon_ptbin%d", ip), -1, -1, binpTlow, binpThigh);
        TH1F *hNsigmaTPCPion = (TH1F *)h3PIDTPCPi->ProjectionX(Form("hNsigmaTPCPion_ptbin%d", ip), -1, -1, binpTlow, binpThigh);
        // TH1F *hNsigmaTPCKaonCut = (TH1F *)h2PIDTPCKaCut->ProjectionX(Form("hNsigmaTPCKaonCut_ptbin%d", ip), binpTlow, binpThigh);
        // TH1F *hNsigmaTPCPionCut = (TH1F *)h2PIDTPCPiCut->ProjectionX(Form("hNsigmaTPCPionCut_ptbin%d", ip), binpTlow, binpThigh);
        // TH1F *hNsigmaTPCKaonCutPr = (TH1F *)h2PIDTPCKaCutPr->ProjectionX(Form("hNsigmaTPCKaonCutPr_ptbin%d", ip), binpTlow, binpThigh);
        // TH1F *hNsigmaTPCPionCutPr = (TH1F *)h2PIDTPCPiCutPr->ProjectionX(Form("hNsigmaTPCPionCutPr_ptbin%d", ip), binpTlow, binpThigh);
        int binptLowNew = h3PIDTPCTOFKa->GetZaxis()->FindBin(pT_bins[ip] + 1e-3);
        int binptHighNew = h3PIDTPCTOFKa->GetZaxis()->FindBin(pT_bins[ip + 1] - 1e-3);
        h3PIDTPCTOFKa->GetZaxis()->SetRange(binptLowNew, binptHighNew);
        h3PIDTPCTOFPi->GetZaxis()->SetRange(binptLowNew, binptHighNew);
        TH2F *hNsigmaTPCTOFKaon = (TH2F *)h3PIDTPCTOFKa->Project3D("xy");
        TH2F *hNsigmaTPCTOFPion = (TH2F *)h3PIDTPCTOFPi->Project3D("xy");
        SetHistoQA(hNsigmaTPCKaon);
        SetHistoQA(hNsigmaTPCPion);
        // SetHistoQA(hNsigmaTPCKaonCut);
        // SetHistoQA(hNsigmaTPCPionCut);
        // SetHistoQA(hNsigmaTPCKaonCutPr);
        // SetHistoQA(hNsigmaTPCPionCutPr);
        cPIDTPCPi->cd(ip + 1);
        // gPad->SetGridx();
        gPad->SetLogy();
        hNsigmaTPCPion->GetXaxis()->SetTitle("TPC N_{#sigma} Pion");
        hNsigmaTPCPion->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCPion->GetXaxis()->SetRangeUser(-5.5, 5.5);
        // hNsigmaTPCPion->Rebin(2);
        // hNsigmaTPCKaonCut->Rebin(2);
        // hNsigmaTPCKaonCutPr->Rebin(2);
        hNsigmaTPCPion->SetMaximum(hNsigmaTPCPion->GetMaximum() * 10);
        hNsigmaTPCPion->Draw("pe");
        // hNsigmaTPCKaonCut->SetLineColor(kRed);
        // hNsigmaTPCKaonCut->Draw("same");
        // hNsigmaTPCKaonCutPr->SetLineColor(kBlue);
        // hNsigmaTPCKaonCutPr->Draw("same");
        lat.DrawLatexNDC(0.15, 0.88, Form("p_{T} = %.2f - %.2f GeV/c", pT_bins[ip], pT_bins[ip + 1]));
        TLine *lineat0 = new TLine(0, 0, 0, hNsigmaTPCPion->GetMaximum() * 1.1);
        lineat0->SetLineColor(kGreen + 2);
        lineat0->SetLineStyle(2);
        lineat0->Draw("same");

        /*
        TF1 *gausFitPion = new TF1("gausFitPion", "gaus", -3, 3);
        gausFitPion->SetParameter(0, hNsigmaTPCPion->GetMaximum());
        gausFitPion->FixParameter(1, 0);
        gausFitPion->SetParameter(2, 1);
        gausFitPion->SetParLimits(2, 0.8, 1.8);
        gausFitPion->SetLineStyle(2);
        gausFitPion->SetLineColor(kBlue);
        hNsigmaTPCPion->Fit(gausFitPion, "REBMSQ");
        TF1 *gausFitKaon = new TF1("gausFitKaon", "gaus", -5, 5);
        gausFitKaon->SetParameter(0, hNsigmaTPCPion->GetMaximum());
        gausFitKaon->SetParameter(1, kaon_contamination_peak_mean[ip]);
        gausFitKaon->SetParLimits(1, kaon_contamination_peak_mean[ip] - 0.7, kaon_contamination_peak_mean[ip] + 0.7);
        gausFitKaon->SetParameter(2, 1);
        if (ip < 9)
            gausFitKaon->SetParLimits(2, 0.5, 1.0);
        else
            gausFitKaon->SetParLimits(2, 0.5, 1.5);
        gausFitKaon->SetLineColor(kGreen + 2);
        gausFitKaon->SetLineStyle(2);
        hNsigmaTPCPion->Fit(gausFitKaon, "REBMSQ");
        TF1 *doubleGausFit = new TF1("doubleGausFit", "gaus(0)+gaus(3)", -5, 5);
        doubleGausFit->SetParameter(0, gausFitPion->GetParameter(0));
        doubleGausFit->SetParameter(1, gausFitPion->GetParameter(1));
        doubleGausFit->SetParLimits(1, gausFitPion->GetParameter(1) - 1, gausFitPion->GetParameter(1) + 1);
        doubleGausFit->SetParameter(2, gausFitPion->GetParameter(2));
        doubleGausFit->SetParameter(3, gausFitKaon->GetParameter(0));
        doubleGausFit->SetParameter(4, gausFitKaon->GetParameter(1));
        doubleGausFit->SetParLimits(4, gausFitKaon->GetParameter(1) - 1, gausFitKaon->GetParameter(1) + 1);
        doubleGausFit->SetParameter(5, gausFitKaon->GetParameter(2));
        doubleGausFit->SetParLimits(5, 0.5, 1.3);
        doubleGausFit->SetLineColor(kRed);
        hNsigmaTPCPion->Fit(doubleGausFit, "REBMSQ");
        doubleGausFit->Draw("same");

        TF1 *fitFucnGausTempPion = new TF1("fitFucnGausTempPion", "gaus", -6, 6);
        fitFucnGausTempPion->SetParameter(0, doubleGausFit->GetParameter(0));
        fitFucnGausTempPion->SetParameter(1, doubleGausFit->GetParameter(1));
        fitFucnGausTempPion->SetParameter(2, doubleGausFit->GetParameter(2));
        fitFucnGausTempPion->SetLineColor(kBlue);
        fitFucnGausTempPion->SetLineStyle(2);
        float purityPion = fitFucnGausTempPion->Integral(-2, 2) * 100 / doubleGausFit->Integral(-2, 2);
        // cout << "Purity Pion for pT bin " << ip << " (" << pT_bins[ip] << " - " << pT_bins[ip + 1] << " GeV/c) : " << purityPion << " %" << endl;
        purityPionVec.push_back(purityPion);
        float purityPionSelected = fitFucnGausTempPion->Integral(pionPIDrange[ip][0], pionPIDrange[ip][1]) * 100 / doubleGausFit->Integral(pionPIDrange[ip][0], pionPIDrange[ip][1]);
        purityPionVecSelected.push_back(purityPionSelected);

        TF1 *fitFucnGausTempKaon = new TF1("fitFucnGausTempKaon", "gaus", -6, 6);
        fitFucnGausTempKaon->SetParameter(0, doubleGausFit->GetParameter(3));
        fitFucnGausTempKaon->SetParameter(1, doubleGausFit->GetParameter(4));
        fitFucnGausTempKaon->SetParameter(2, doubleGausFit->GetParameter(5));
        fitFucnGausTempKaon->SetLineColor(kGreen + 2);
        fitFucnGausTempKaon->SetLineStyle(2);
        fitFucnGausTempKaon->Draw("same");
        fitFucnGausTempPion->Draw("same");

        TLegend *legPID = new TLegend(0.4, 0.25, 0.6, 0.5);
        legPID->SetBorderSize(0);
        legPID->SetTextFont(42);
        legPID->SetTextSize(0.08);
        // legPID->SetFillStyle(0);
        legPID->AddEntry(gausFitPion, "Gaussian 1", "l");
        legPID->AddEntry(gausFitKaon, "Gaussian 2", "l");
        legPID->AddEntry(doubleGausFit, "Double Gaussian", "l");
        if (ip == 0)
            legPID->Draw("same");

        */

        cPIDTPCKa->cd(ip + 1);
        // gPad->SetGridx();
        gPad->SetLogy();
        hNsigmaTPCKaon->GetXaxis()->SetTitle("TPC N_{#sigma} Kaon");
        hNsigmaTPCKaon->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCKaon->GetXaxis()->SetRangeUser(-5.5, 5.5);
        // hNsigmaTPCKaon->Rebin(2);
        // hNsigmaTPCKaonCut->Rebin(2);
        // hNsigmaTPCKaonCutPr->Rebin(2);
        hNsigmaTPCKaon->SetMaximum(hNsigmaTPCKaon->GetMaximum() * 10);
        hNsigmaTPCKaon->Draw("pe");
        // hNsigmaTPCKaonCut->SetLineColor(kRed);
        // hNsigmaTPCKaonCut->Draw("same");
        // hNsigmaTPCKaonCutPr->SetLineColor(kBlue);
        // hNsigmaTPCKaonCutPr->Draw("same");
        lat.DrawLatexNDC(0.15, 0.88, Form("p_{T} = %.2f - %.2f GeV/c", pT_bins[ip], pT_bins[ip + 1]));
        TLine *lineat0ka = new TLine(0, 0, 0, hNsigmaTPCKaon->GetMaximum() * 1.1);
        lineat0ka->SetLineColor(kGreen + 2);
        lineat0ka->SetLineStyle(2);
        lineat0ka->Draw("same");

        /*
        TF1 *gausFitKaon2 = new TF1("gausFitKaon2", "gaus", -3, 3);
        gausFitKaon2->SetParameter(0, hNsigmaTPCKaon->GetMaximum());
        gausFitKaon2->FixParameter(1, 0);
        gausFitKaon2->SetParameter(2, 1);
        gausFitKaon2->SetParLimits(2, 0.8, 1.8);
        gausFitKaon2->SetLineStyle(2);
        gausFitKaon2->SetLineColor(kBlue);
        hNsigmaTPCKaon->Fit(gausFitKaon2, "REBMS0Q");
        TF1 *gausFitPion2 = new TF1("gausFitPion2", "gaus", -5, 5);
        gausFitPion2->SetParameter(0, hNsigmaTPCKaon->GetMaximum());
        gausFitPion2->SetParameter(1, pion_contamination_peak_mean[ip]);
        gausFitPion2->SetParLimits(1, pion_contamination_peak_mean[ip] - 0.6, pion_contamination_peak_mean[ip] + 0.6);
        gausFitPion2->SetParameter(2, 1);
        gausFitPion2->SetParLimits(2, 0.5, 1.5);
        gausFitPion2->SetLineColor(kGreen + 2);
        gausFitPion2->SetLineStyle(2);
        hNsigmaTPCKaon->Fit(gausFitPion2, "REBMS0Q");

        TF1 *doubleGausFit2 = new TF1("doubleGausFit2", "gaus(0)+gaus(3)", -5, 5);
        doubleGausFit2->SetParameter(0, gausFitKaon2->GetParameter(0));
        doubleGausFit2->FixParameter(1, gausFitKaon2->GetParameter(1));
        doubleGausFit2->SetParLimits(1, gausFitKaon2->GetParameter(1) - 1, gausFitKaon2->GetParameter(1) + 1);
        doubleGausFit2->SetParameter(2, gausFitKaon2->GetParameter(2));
        doubleGausFit2->SetParameter(3, gausFitPion2->GetParameter(0));
        doubleGausFit2->SetParameter(4, gausFitPion2->GetParameter(1));
        doubleGausFit2->SetParLimits(4, gausFitPion2->GetParameter(1) - 0.8, gausFitPion2->GetParameter(1) + 0.8);
        doubleGausFit2->SetParameter(5, gausFitPion2->GetParameter(2));
        doubleGausFit2->SetLineColor(kRed);
        hNsigmaTPCKaon->Fit(doubleGausFit2, "REBMSQ");
        doubleGausFit2->Draw("same");
        if (ip == 0)
            legPID->Draw("same");

        TF1 *fitFucnGausTempKaon2 = new TF1("fitFucnGausTempKaon2", "gaus", -6, 6);
        fitFucnGausTempKaon2->SetParameter(0, doubleGausFit2->GetParameter(0));
        fitFucnGausTempKaon2->SetParameter(1, doubleGausFit2->GetParameter(1));
        fitFucnGausTempKaon2->SetParameter(2, doubleGausFit2->GetParameter(2));
        fitFucnGausTempKaon2->SetLineColor(kBlue);
        fitFucnGausTempKaon2->SetLineStyle(2);
        float purityKaon = fitFucnGausTempKaon2->Integral(-2, 2) * 100 / doubleGausFit2->Integral(-2, 2);
        cout << "Purity Kaon for pT bin " << ip << " (" << pT_bins[ip] << "-" << pT_bins[ip + 1] << " GeV/c) : " << purityKaon << " %" << endl;
        purityKaonVec.push_back(purityKaon);
        float purityKaonSelected = fitFucnGausTempKaon2->Integral(kaonPIDrange[ip][0], kaonPIDrange[ip][1]) * 100 / doubleGausFit2->Integral(kaonPIDrange[ip][0], kaonPIDrange[ip][1]);
        purityKaonVecSelected.push_back(purityKaonSelected);

        if (purityKaon < 10)
        {
            cout << "Mean of Kaon Gauss " << gausFitKaon2->GetParameter(1) << " , mean of pion Gauss " << gausFitPion2->GetParameter(1) << endl;
        }

        TF1 *fitFucnGausTempPion2 = new TF1("fitFucnGausTempPion2", "gaus", -6, 6);
        fitFucnGausTempPion2->SetParameter(0, doubleGausFit2->GetParameter(3));
        fitFucnGausTempPion2->SetParameter(1, doubleGausFit2->GetParameter(4));
        fitFucnGausTempPion2->SetParameter(2, doubleGausFit2->GetParameter(5));
        fitFucnGausTempPion2->SetLineColor(kGreen + 2);
        fitFucnGausTempPion2->SetLineStyle(2);
        fitFucnGausTempPion2->Draw("same");
        fitFucnGausTempKaon2->Draw("same");
        */

        cPIDTPCTOFKa->cd(ip + 1);
        hNsigmaTPCTOFKaon->GetXaxis()->SetTitle("TPC N_{#sigma} Kaon");
        hNsigmaTPCTOFKaon->GetYaxis()->SetTitle("TOF N_{#sigma} Kaon");
        // hNsigmaTPCTOFKaon->RebinX(5);
        // hNsigmaTPCTOFKaon->RebinY(5);
        hNsigmaTPCTOFKaon->GetXaxis()->SetRangeUser(-5, 5);
        hNsigmaTPCTOFKaon->GetYaxis()->SetRangeUser(-5, 5);
        hNsigmaTPCTOFKaon->Draw("COLZ");
        lat.DrawLatexNDC(0.15, 0.85, Form("p_{T} = %.2f - %.2f GeV/c", pT_bins[ip], pT_bins[ip + 1]));
        TLine *line0x = new TLine(0, -5, 0, 5);
        line0x->SetLineColor(kGreen + 2);
        line0x->SetLineStyle(2);
        line0x->Draw("same");
        TLine *line0y = new TLine(-5, 0, 5, 0);
        line0y->SetLineColor(kGreen + 2);
        line0y->SetLineStyle(2);
        line0y->Draw("same");

        cPIDTPCTOFPi->cd(ip + 1);
        hNsigmaTPCTOFPion->GetXaxis()->SetTitle("TPC N_{#sigma} Pion");
        hNsigmaTPCTOFPion->GetYaxis()->SetTitle("TOF N_{#sigma} Pion");
        // hNsigmaTPCTOFPion->RebinX(5);
        // hNsigmaTPCTOFPion->RebinY(5);
        hNsigmaTPCTOFPion->GetXaxis()->SetRangeUser(-5, 5);
        hNsigmaTPCTOFPion->GetYaxis()->SetRangeUser(-5, 5);
        hNsigmaTPCTOFPion->Draw("COLZ");
        lat.DrawLatexNDC(0.15, 0.85, Form("p_{T} = %.2f - %.2f GeV/c", pT_bins[ip], pT_bins[ip + 1]));
        TLine *line0xpi = new TLine(0, -5, 0, 5);
        line0xpi->SetLineColor(kGreen + 2);
        line0xpi->SetLineStyle(2);
        line0xpi->Draw("same");
        TLine *line0ypi = new TLine(-5, 0, 5, 0);
        line0ypi->SetLineColor(kGreen + 2);
        line0ypi->SetLineStyle(2);
        line0ypi->Draw("same");
    }
    // cPIDTPCPi->SaveAs("PIDchecks/kstar_TPCnsigma_Pion_beforeFit.png");
    // cPIDTPCKa->SaveAs("PIDchecks/kstar_TPCnsigma_Kaon_beforeFit.png");

    // TCanvas *cPurityPion = new TCanvas("cPurityPion", "PurityPion", 720, 720);
    // SetCanvasStyle(cPurityPion, 0.13, 0.03, 0.03, 0.13);
    // TGraph *gPurityPion = new TGraph(total_pT_bins, pT_bins, &purityPionVec[0]);
    // TGraph *gPurityPionSelected = new TGraph(total_pT_bins, pT_bins, &purityPionVecSelected[0]);
    // SetGraphStyle(gPurityPion);
    // SetGraphStyle(gPurityPionSelected);
    // gPurityPion->SetTitle(0);
    // gPurityPion->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // gPurityPion->GetYaxis()->SetTitle("Pion Purity (%)");
    // gPurityPion->GetYaxis()->SetRangeUser(30, 105);
    // gPurityPion->SetMarkerStyle(20);
    // gPurityPion->Draw("APL");
    // gPurityPionSelected->SetMarkerStyle(21);
    // gPurityPionSelected->SetMarkerColor(kRed);
    // gPurityPionSelected->SetLineColor(kRed);
    // gPurityPionSelected->Draw("PLsame");
    // TLegend *legPurityPion = new TLegend(0.4, 0.2, 0.7, 0.4);
    // legPurityPion->SetBorderSize(0);
    // legPurityPion->SetTextFont(42);
    // legPurityPion->SetTextSize(0.04);
    // legPurityPion->AddEntry(gPurityPion, "|N_{#sigma_{TPC}}| < #pm2", "p");
    // legPurityPion->AddEntry(gPurityPionSelected, "#it{p_{T}} dependent selection", "p");
    // legPurityPion->Draw("same");
    // TLine *lineat90Pion = new TLine(0, 90, 11.0, 90);
    // lineat90Pion->SetLineColor(kGray + 2);
    // lineat90Pion->SetLineStyle(2);
    // lineat90Pion->Draw("same");
    // cPurityPion->SaveAs("PIDchecks/kstar_PionPurity_vs_pT.png");

    // TCanvas *cPurityKaon = new TCanvas("cPurityKaon", "PurityKaon", 720, 720);
    // SetCanvasStyle(cPurityKaon, 0.13, 0.03, 0.03, 0.13);
    // TGraph *gPurityKaon = new TGraph(total_pT_bins, pT_bins, &purityKaonVec[0]);
    // TGraph *gPurityKaonSelected = new TGraph(total_pT_bins, pT_bins, &purityKaonVecSelected[0]);
    // SetGraphStyle(gPurityKaon);
    // SetGraphStyle(gPurityKaonSelected);
    // gPurityKaon->SetTitle(0);
    // gPurityKaon->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // gPurityKaon->GetYaxis()->SetTitle("Kaon Purity (%)");
    // gPurityKaon->GetYaxis()->SetRangeUser(30, 105);
    // gPurityKaon->SetMarkerStyle(20);
    // gPurityKaon->Draw("APL");
    // gPurityKaonSelected->SetMarkerStyle(21);
    // gPurityKaonSelected->SetMarkerColor(kRed);
    // gPurityKaonSelected->SetLineColor(kRed);
    // gPurityKaonSelected->Draw("PLsame");
    // legPurityPion->Draw("same");
    // lineat90Pion->Draw("same");
    // cPurityKaon->SaveAs("PIDchecks/kstar_KaonPurity_vs_pT.png");
}