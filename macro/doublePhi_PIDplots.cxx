
#include <iostream>
#include "src/style.h"
// #include "src/fitfunc.h"
// #include "src/initializations.h"

void doublePhi_PIDplots()
{
    gStyle->SetOptStat(0);
    // TString subWagon = "";
    // TString subWagon = "_PID1000";
    TString subWagon = "_pid2_opti5";

    TString filepath = "../data/doublePhi/";
    // TString filepath = "/home/sawan/alice/practice/DoublePhiAnalyzedFiles/";

    // TFile *file = new TFile(filepath + "/Analysis_ac_ah_MergedPID4ME.root");
    // TFile *file = new TFile(filepath + "/AnalysisResults_ac_ah_new.root");
    TFile *file = new TFile(filepath + "/680311.root");
    if (file->IsZombie())
    {
        cerr << "File not found " << endl;
        return;
    }

    TH2D *hnSigmaTPCKaonPlus = (TH2D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTPCKaonPlusBefore");
    TH2D *hnSigmaTPCKaonMinus = (TH2D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTPCKaonMinusBefore");
    TH3D *hnSigmaTPCTOFKaonBefore = (TH3D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTPCTOFKaonBefore");
    if (hnSigmaTPCKaonPlus == nullptr || hnSigmaTPCKaonMinus == nullptr || hnSigmaTPCTOFKaonBefore == nullptr)
    {
        std::cerr << "Error: Could not find histogram 'doublephimeson" + subWagon + "/hnsigmaTPCKaonPlus' in file\n";
        return;
    }

    TCanvas *cPIDTPCKaPos = new TCanvas("cPIDTPCKaPos", "PIDTPCKa", 1440, 800);
    SetCanvasStyle(cPIDTPCKaPos, 0.1, 0.05, 0.06, 0.17);
    cPIDTPCKaPos->cd();
    TPad *mainPadPos = new TPad("mainPadPos", "", 0.0, 0.0, 1.0, 0.96);
    mainPadPos->SetMargin(0.0, 0.0, 0.0, 0.0);
    mainPadPos->Draw();
    mainPadPos->cd();
    mainPadPos->Divide(4, 4, 0.001, 0.001);

    TCanvas *cPIDTPCKaNeg = new TCanvas("cPIDTPCKaNeg", "PIDTPCKa", 1440, 800);
    SetCanvasStyle(cPIDTPCKaNeg, 0.1, 0.05, 0.06, 0.17);
    cPIDTPCKaNeg->cd();
    TPad *mainPadNeg = new TPad("mainPadNeg", "", 0.0, 0.0, 1.0, 0.96);
    mainPadNeg->SetMargin(0.0, 0.0, 0.0, 0.0);
    mainPadNeg->Draw();
    mainPadNeg->cd();
    mainPadNeg->Divide(4, 4, 0.001, 0.001);

    TCanvas *cPIDTPCTOFKa = new TCanvas("cPIDTPCTOFKa", "PIDTPCTOFKa", 1440, 800);
    SetCanvasStyle(cPIDTPCTOFKa, 0.1, 0.05, 0.06, 0.17);
    cPIDTPCTOFKa->cd();
    TPad *mainPadTPCTOF = new TPad("mainPadTPCTOF", "", 0.0, 0.0, 1.0, 0.96);
    mainPadTPCTOF->SetMargin(0.0, 0.0, 0.0, 0.0);
    mainPadTPCTOF->Draw();
    mainPadTPCTOF->cd();
    mainPadTPCTOF->Divide(3, 2, 0.001, 0.001);

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.09);

    TLatex lat2;
    lat2.SetNDC();
    lat2.SetTextFont(42);
    lat2.SetTextSize(0.05);

    float pT_bins[] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 5.0, 10.0};
    int total_pT_bins = sizeof(pT_bins) / sizeof(pT_bins[0]) - 1;
    float pion_contamination_peak_mean[] = {-7, -7, -6, -5, -3, -2, -1, 1, 1, 2, 3, 3, 3, 3, 3, 3};
    vector<vector<float>> kaonPIDrange = {{-2.7, 2.3}, {-2.5, 2.5}, {-1.5, 2.5}, {-0.75, 2.5}, {-0.3, 2.5}, {-0.0, 2.5}, {-0.0, 2.5}, {-2, 0.0}, {-1.5, 0.0}, {-1.5, 0.0}, {-1, 0.0}, {-2.5, 1}, {-2.5, 1}, {-2.5, 1}, {-2.5, 1}, {-2.5, 1.5}};

    float pTbins_forCombinedTPCTOF[] = {0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0};
    int nPtBins_forCombinedTPCTOF = sizeof(pTbins_forCombinedTPCTOF) / sizeof(float) - 1;

    vector<float> purityPionVec, purityKaonVec, purityPionVecSelected, purityKaonVecSelected;

    for (Int_t ip = 0; ip < total_pT_bins; ip++) // start pt bin loop
    {
        int binLow = hnSigmaTPCKaonPlus->GetYaxis()->FindBin(pT_bins[ip] + 0.001);
        int binHigh = hnSigmaTPCKaonPlus->GetYaxis()->FindBin(pT_bins[ip + 1] - 0.001);
        TH1D *hKaonTPC1Dpos = hnSigmaTPCKaonPlus->ProjectionX(Form("hKaonTPC1Dpos_%.1f_%.1f", pT_bins[ip], pT_bins[ip + 1]), binLow, binHigh);
        TH1D *hKaonTPC1Dneg = hnSigmaTPCKaonMinus->ProjectionX(Form("hKaonTPC1Dneg_%.1f_%.1f", pT_bins[ip], pT_bins[ip + 1]), binLow, binHigh);
        SetHistoQA(hKaonTPC1Dpos);
        SetHistoQA(hKaonTPC1Dneg);

        mainPadPos->cd(ip + 1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.02);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.07);
        hKaonTPC1Dpos->GetXaxis()->SetTitle("TPC N_{#sigma} Kaon");
        hKaonTPC1Dpos->GetYaxis()->SetTitle("Counts");
        hKaonTPC1Dpos->GetXaxis()->SetRangeUser(-5.5, 5.5);
        hKaonTPC1Dpos->SetMaximum(hKaonTPC1Dpos->GetMaximum() * 7);
        hKaonTPC1Dpos->SetLineColor(kBlue + 1);
        hKaonTPC1Dpos->SetLineWidth(2);
        hKaonTPC1Dpos->Draw("HIST");

        lat.DrawLatexNDC(0.2, 0.83, Form("p_{T} = %.2f - %.2f GeV/c", pT_bins[ip], pT_bins[ip + 1]));
        // lat.DrawLatexNDC(0.2, 0.69, Form("%.1f < n#sigma_{TPC} < %.1f", kaonPIDrange[ip][0], kaonPIDrange[ip][1]));
        TLine *lineat0ka = new TLine(0, 0, 0, hKaonTPC1Dpos->GetMaximum() * 1.1);
        lineat0ka->SetLineColor(kGreen + 2);
        lineat0ka->SetLineStyle(2);
        lineat0ka->Draw("same");

        TF1 *gausFitKaon = new TF1("gausFitKaon", "gaus", -3, 3);
        gausFitKaon->SetParameter(0, hKaonTPC1Dpos->GetMaximum());
        gausFitKaon->FixParameter(1, 0);
        gausFitKaon->SetParameter(2, 1);
        gausFitKaon->SetParLimits(2, 0.8, 1.8);
        gausFitKaon->SetLineStyle(2);
        gausFitKaon->SetLineColor(kBlue);
        hKaonTPC1Dpos->Fit(gausFitKaon, "REBMS0Q");

        TF1 *gausFitPion = new TF1("gausFitPion", "gaus", -5, 5);
        gausFitPion->SetParameter(0, hKaonTPC1Dpos->GetMaximum());
        gausFitPion->SetParameter(1, pion_contamination_peak_mean[ip]);
        gausFitPion->SetParLimits(1, pion_contamination_peak_mean[ip] - 0.6, pion_contamination_peak_mean[ip] + 0.6);
        gausFitPion->SetParameter(2, 1);
        gausFitPion->SetParLimits(2, 0.5, 1.5);
        gausFitPion->SetLineColor(kGreen + 2);
        gausFitPion->SetLineStyle(2);
        hKaonTPC1Dpos->Fit(gausFitPion, "REBMS0Q");

        TF1 *doubleGausFit = new TF1("doubleGausFit", "gaus(0)+gaus(3)", -5, 5);
        doubleGausFit->SetParameter(0, gausFitKaon->GetParameter(0));
        doubleGausFit->FixParameter(1, gausFitKaon->GetParameter(1));
        doubleGausFit->SetParLimits(1, gausFitKaon->GetParameter(1) - 1, gausFitKaon->GetParameter(1) + 1);
        doubleGausFit->SetParameter(2, gausFitKaon->GetParameter(2));
        doubleGausFit->SetParameter(3, gausFitPion->GetParameter(0));
        doubleGausFit->SetParameter(4, gausFitPion->GetParameter(1));
        doubleGausFit->SetParLimits(4, gausFitPion->GetParameter(1) - 0.8, gausFitPion->GetParameter(1) + 0.8);
        doubleGausFit->SetParameter(5, gausFitPion->GetParameter(2));
        doubleGausFit->SetLineColor(kRed);
        hKaonTPC1Dpos->Fit(doubleGausFit, "REBMSQ");
        doubleGausFit->Draw("same");

        if (ip == 0)
        {
            TLegend *legPID = new TLegend(0.4, 0.25, 0.6, 0.5);
            legPID->SetBorderSize(0);
            legPID->SetTextFont(42);
            legPID->SetTextSize(0.08);
            // legPID->SetFillStyle(0);
            legPID->AddEntry(gausFitPion, "Gaussian 1", "l");
            legPID->AddEntry(gausFitKaon, "Gaussian 2", "l");
            legPID->AddEntry(doubleGausFit, "Double Gaussian", "l");
            legPID->Draw("same");
        }

        TF1 *fitFucnGausTempKaon = new TF1("fitFucnGausTempKaon", "gaus", -6, 6);
        fitFucnGausTempKaon->SetParameter(0, doubleGausFit->GetParameter(0));
        fitFucnGausTempKaon->SetParameter(1, doubleGausFit->GetParameter(1));
        fitFucnGausTempKaon->SetParameter(2, doubleGausFit->GetParameter(2));
        fitFucnGausTempKaon->SetLineColor(kBlue);
        fitFucnGausTempKaon->SetLineStyle(2);
        float purityKaon = fitFucnGausTempKaon->Integral(-2.0, 2.0) * 100 / doubleGausFit->Integral(-2.0, 2.0);
        cout << "Purity Kaon for pT bin " << ip << " (" << pT_bins[ip] << "-" << pT_bins[ip + 1] << " GeV/c) : " << purityKaon << " %" << endl;
        purityKaonVec.push_back(purityKaon);
        float purityKaonSelected = fitFucnGausTempKaon->Integral(kaonPIDrange[ip][0], kaonPIDrange[ip][1]) * 100 / doubleGausFit->Integral(kaonPIDrange[ip][0], kaonPIDrange[ip][1]);
        purityKaonVecSelected.push_back(purityKaonSelected);

        if (purityKaon < 10)
        {
            cout << "Mean of Kaon Gauss " << gausFitKaon->GetParameter(1) << " , mean of pion Gauss " << gausFitPion->GetParameter(1) << endl;
        }

        TF1 *fitFucnGausTempPion = new TF1("fitFucnGausTempPion", "gaus", -6, 6);
        fitFucnGausTempPion->SetParameter(0, doubleGausFit->GetParameter(3));
        fitFucnGausTempPion->SetParameter(1, doubleGausFit->GetParameter(4));
        fitFucnGausTempPion->SetParameter(2, doubleGausFit->GetParameter(5));
        fitFucnGausTempPion->SetLineColor(kGreen + 2);
        fitFucnGausTempPion->SetLineStyle(2);
        fitFucnGausTempPion->Draw("same");
        fitFucnGausTempKaon->Draw("same");
        lat.DrawLatexNDC(0.2, 0.71, Form("Mean1: %.1f, #sigma_{1}: %.1f", doubleGausFit->GetParameter(1), doubleGausFit->GetParameter(2)));
        lat.DrawLatexNDC(0.2, 0.62, Form("Mean2: %.1f, #sigma_{2}: %.1f", doubleGausFit->GetParameter(4), doubleGausFit->GetParameter(5)));

        mainPadNeg->cd(ip + 1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.02);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.07);
        hKaonTPC1Dneg->GetXaxis()->SetTitle("TPC N_{#sigma} Kaon");
        hKaonTPC1Dneg->GetYaxis()->SetTitle("Counts");
        hKaonTPC1Dneg->GetXaxis()->SetRangeUser(-5.5, 5.5);
        hKaonTPC1Dneg->SetMaximum(hKaonTPC1Dneg->GetMaximum() * 10);
        hKaonTPC1Dneg->SetLineColor(kRed + 1);
        hKaonTPC1Dneg->SetLineWidth(2);
        hKaonTPC1Dneg->Draw("HIST");

        if (ip < nPtBins_forCombinedTPCTOF)
        {
            mainPadTPCTOF->cd(ip + 1);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.02);
            gPad->SetBottomMargin(0.13);
            gPad->SetTopMargin(0.07);

            float ptLow = pTbins_forCombinedTPCTOF[ip];
            float ptHigh = pTbins_forCombinedTPCTOF[ip + 1];

            int binLowCombined = hnSigmaTPCTOFKaonBefore->GetZaxis()->FindBin(ptLow + 0.001);
            int binHighCombined = hnSigmaTPCTOFKaonBefore->GetZaxis()->FindBin(ptHigh - 0.001);

            hnSigmaTPCTOFKaonBefore->GetZaxis()->SetRange(binLowCombined, binHighCombined);
            TH2D *hnSigmaTPCTOFKaon = (TH2D *)hnSigmaTPCTOFKaonBefore->Project3D("xy");
            hnSigmaTPCTOFKaon->SetName(Form("hnSigmaTPCTOFKaon_%.1f_%.1f", ptLow, ptHigh));
            SetHistoQA(hnSigmaTPCTOFKaon);
            hnSigmaTPCTOFKaon->SetTitle("");
            hnSigmaTPCTOFKaon->GetXaxis()->SetTitle("TPC N_{#sigma} Ka");
            hnSigmaTPCTOFKaon->GetYaxis()->SetTitle("TOF N_{#sigma} Ka");
            // hnSigmaTPCTOFKaon->GetXaxis()->SetRangeUser(-5, 5);
            // hnSigmaTPCTOFKaon->GetYaxis()->SetRangeUser(-5, 5);
            hnSigmaTPCTOFKaon->SetStats(0);
            hnSigmaTPCTOFKaon->Draw("colz");

            lat2.SetTextColor(kWhite);
            lat2.DrawLatexNDC(0.3, 0.85, Form("p_{T} = %.2f - %.2f GeV/c", ptLow, ptHigh));
            TLine *line0x = new TLine(0, -5, 0, 5);
            line0x->SetLineColor(kGreen + 2);
            line0x->SetLineStyle(2);
            line0x->Draw("same");
            TLine *line0y = new TLine(-5, 0, 5, 0);
            line0y->SetLineColor(kGreen + 2);
            line0y->SetLineStyle(2);
            line0y->Draw("same");
        }
    }
    cPIDTPCKaPos->cd();
    TPad *titlePad6 = new TPad("titlePad6", "", 0.0, 0.94, 1.0, 1.0);
    titlePad6->SetFillStyle(0);
    titlePad6->SetFrameFillStyle(0);
    titlePad6->SetBorderMode(0);
    titlePad6->SetMargin(0, 0, 0, 0);
    titlePad6->Draw();
    titlePad6->cd();
    TLatex latCanvasTitle6;
    latCanvasTitle6.SetNDC();
    latCanvasTitle6.SetTextFont(42);
    latCanvasTitle6.SetTextSize(0.45);
    latCanvasTitle6.SetTextAlign(22);
    latCanvasTitle6.DrawLatex(0.5, 0.5, "n#sigma_{TPC} Ka");
    cPIDTPCKaPos->Modified();
    cPIDTPCKaPos->Update();
    cPIDTPCKaPos->SaveAs("doublePhi/PIDPlots/TPC_Kaon_pos.png");

    cPIDTPCKaNeg->cd();
    TPad *titlePad6Neg = new TPad("titlePad6Neg", "", 0.0, 0.94, 1.0, 1.0);
    titlePad6Neg->SetFillStyle(0);
    titlePad6Neg->SetFrameFillStyle(0);
    titlePad6Neg->SetBorderMode(0);
    titlePad6Neg->SetMargin(0, 0, 0, 0);
    titlePad6Neg->Draw();
    titlePad6Neg->cd();
    TLatex latCanvasTitle6Neg;
    latCanvasTitle6Neg.SetNDC();
    latCanvasTitle6Neg.SetTextFont(42);
    latCanvasTitle6Neg.SetTextSize(0.45);
    latCanvasTitle6Neg.SetTextAlign(22);
    latCanvasTitle6Neg.DrawLatex(0.5, 0.5, "n#sigma_{TPC} Ka");
    cPIDTPCKaNeg->Modified();
    cPIDTPCKaNeg->Update();
    cPIDTPCKaNeg->SaveAs("doublePhi/PIDPlots/TPC_Kaon_neg.png");

    cPIDTPCTOFKa->cd();
    TPad *titlePad7 = new TPad("titlePad7", "", 0.0, 0.94, 1.0, 1.0);
    titlePad7->SetFillStyle(0);
    titlePad7->SetFrameFillStyle(0);
    titlePad7->SetBorderMode(0);
    titlePad7->SetMargin(0, 0, 0, 0);
    titlePad7->Draw();
    titlePad7->cd();
    TLatex latCanvasTitle7;
    latCanvasTitle7.SetNDC();
    latCanvasTitle7.SetTextFont(42);
    latCanvasTitle7.SetTextSize(0.45);
    latCanvasTitle7.SetTextAlign(22);
    latCanvasTitle7.DrawLatex(0.5, 0.5, "n#sigma_{TPC+TOF} Ka");
    cPIDTPCTOFKa->Modified();
    cPIDTPCTOFKa->Update();
    cPIDTPCTOFKa->SaveAs("doublePhi/PIDPlots/TOF_Kaon.png");

    TCanvas *cPurityKaon = new TCanvas("cPurityKaon", "PurityKaon", 720, 720);
    SetCanvasStyle(cPurityKaon, 0.13, 0.03, 0.03, 0.13);
    TGraph *gPurityKaon = new TGraph(total_pT_bins, pT_bins, &purityKaonVec[0]);
    TGraph *gPurityKaonSelected = new TGraph(total_pT_bins, pT_bins, &purityKaonVecSelected[0]);
    SetGraphStyle(gPurityKaon);
    SetGraphStyle(gPurityKaonSelected);
    gPurityKaon->SetTitle(0);
    gPurityKaon->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gPurityKaon->GetYaxis()->SetTitle("Kaon Purity (%)");
    gPurityKaon->GetYaxis()->SetRangeUser(30, 105);
    gPurityKaon->SetMarkerStyle(20);
    gPurityKaon->Draw("APL");
    gPurityKaonSelected->SetMarkerStyle(21);
    gPurityKaonSelected->SetMarkerColor(kRed);
    gPurityKaonSelected->SetLineColor(kRed);
    gPurityKaonSelected->Draw("PLsame");
    TLegend *legPurityKaon = new TLegend(0.4, 0.2, 0.7, 0.4);
    legPurityKaon->SetBorderSize(0);
    legPurityKaon->SetTextFont(42);
    legPurityKaon->SetTextSize(0.04);
    legPurityKaon->SetHeader("Kaon TPC PID");
    legPurityKaon->AddEntry(gPurityKaon, "|N_{#sigma_{TPC}}| < #pm2", "p");
    legPurityKaon->AddEntry(gPurityKaonSelected, "#it{p_{T}} dependent selection", "p");
    legPurityKaon->Draw("same");
    TLine *lineat90Kaon = new TLine(0, 90, 11.0, 90);
    lineat90Kaon->SetLineColor(kGray + 2);
    lineat90Kaon->SetLineStyle(2);
    lineat90Kaon->Draw("same");
    cPurityKaon->SaveAs("doublePhi/PIDPlots/KaonPurity.png");
}