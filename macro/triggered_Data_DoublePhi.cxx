#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"
using namespace std;

void triggered_Data_DoublePhi()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetGridColor(15);

    // *******************LHC25ah*********************
    string runNumbers[] = {"569945", "569947", "569978", "569980", "569981", "569982", "570011", "570012", "570025", "570036", "570049", "570050"};

    string path = "/home/sawan/TriggeredDataFiles/";
    int totalFiles = sizeof(runNumbers) / sizeof(runNumbers[0]);

    TH1F *hselectivity = new TH1F("hselectivity", "Selectivity;Run Number;Selectivity", totalFiles, 1, totalFiles + 1);

    TH1F *hPhiMassFit = new TH1F("hPhiMassFit", "", totalFiles, 1, totalFiles + 1);
    TH1F *hPhiYieldFit = new TH1F("hPhiYieldFit", "", totalFiles, 1, totalFiles + 1);
    TH1F *hPhiMassResolution = new TH1F("hPhiMassResolution", "", totalFiles, 1, totalFiles + 1);
    TH1F *hChi2byNDF = new TH1F("hChi2byNDF", "", totalFiles, 1, totalFiles + 1);

    //******Fitting for Phi*********************

    // TCanvas *c2 = new TCanvas("c2", "Phi Mass", 1080, 720);
    // SetCanvasStyle(c2, 0.13, 0.06, 0.05, 0.13);
    // c2->Divide(4, 3);
    TCanvas *c3 = new TCanvas("c3", "Double Phi Mass", 1080, 720);
    SetCanvasStyle(c3, 0.13, 0.06, 0.05, 0.13);
    c3->Divide(4, 3);
    // TCanvas *c4 = new TCanvas("c4", "Nsigma TPC Kaon", 1080, 720);
    // SetCanvasStyle(c4, 0.13, 0.14, 0.05, 0.13);
    // c4->Divide(4, 3);
    // TCanvas *c5 = new TCanvas("c5", "Nsigma TOF Kaon", 1080, 720);
    // SetCanvasStyle(c5, 0.13, 0.14, 0.05, 0.13);
    // c5->Divide(4, 3);
    TCanvas *c6 = new TCanvas("c6", "Nsigma TPC Kaon", 1440, 1080);
    SetCanvasStyle(c6, 0.13, 0.14, 0.05, 0.13);
    TCanvas *c7 = new TCanvas("c7", "Nsigma TOF Kaon", 1440, 1080);
    SetCanvasStyle(c7, 0.13, 0.14, 0.05, 0.13);

    c6->cd();
    TPad *mainPad6 = new TPad("mainPad6", "", 0.0, 0.0, 1.0, 0.94);
    mainPad6->SetMargin(0.0, 0.0, 0.0, 0.0);
    mainPad6->Draw();
    mainPad6->cd();
    mainPad6->Divide(4, 3, 0.001, 0.001);

    c7->cd();
    TPad *mainPad7 = new TPad("mainPad7", "", 0.0, 0.0, 1.0, 0.94);
    mainPad7->SetMargin(0.0, 0.0, 0.0, 0.0);
    mainPad7->Draw();
    mainPad7->cd();
    mainPad7->Divide(4, 3, 0.001, 0.001);
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.08);

    for (int ifiles = 0; ifiles < totalFiles; ifiles++)
    {
        TFile *file = new TFile((path + "AnalysisResults_" + runNumbers[ifiles] + ".root").c_str(), "READ");
        if (file->IsZombie())
        {
            std::cerr << "Error: Could not open file triggered_data/AnalysisResults_fullrun.root\n";
            return;
        }

        TH1F *hevent_processed = (TH1F *)file->Get("lf-doublephi-filter/hProcessedEvents");
        TH2F *hphiMass = (TH2F *)file->Get("lf-doublephi-filter/hInvMassPhi");
        TH2F *hdoublePhiMass = (TH2F *)file->Get("lf-doublephi-filter/hInvMassDoublePhi");
        TH2F *hNsigmaTPCKaon = (TH2F *)file->Get("lf-doublephi-filter/hNsigmaPtkaonTPC");
        TH2F *hNsigmaTOFKaon = (TH2F *)file->Get("lf-doublephi-filter/hNsigmaPtkaonTOF");

        if (hevent_processed == nullptr || hphiMass == nullptr || hdoublePhiMass == nullptr || hNsigmaTPCKaon == nullptr || hNsigmaTOFKaon == nullptr)
        {
            std::cerr << "Error: Could not find one or more histograms in the file\n";
            return;
        }

        hselectivity->SetBinContent(ifiles + 1, hevent_processed->GetBinContent(3) / hevent_processed->GetBinContent(1));
        cout << "Selectivity for run " << runNumbers[ifiles] << " is " << hevent_processed->GetBinContent(3) / hevent_processed->GetBinContent(1) << endl;
        hselectivity->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        // /*

        // c2->cd(ifiles + 1);
        // int lowbinpT = hphiMass->GetYaxis()->FindBin(0.0 + 0.001);
        // int highbinpT = hphiMass->GetYaxis()->FindBin(10.0 - 0.001);
        // TH1D *hphiMassProj = hphiMass->ProjectionX(Form("hphiMassProj_%d", ifiles), lowbinpT, highbinpT);
        // TF1 *fitFcn = new TF1(Form("fitfunc_%d", ifiles), voigtpol2, 1.00, 1.4, 7);      // sig+bkg fit function
        // TF1 *fitFcnBkg = new TF1(Form("fitfunc1_%d", ifiles), polynomial2, 1.0, 1.4, 3); // only residualbkg
        // TF1 *fitFcnSig = new TF1(Form("fitFcnSig_%d", ifiles), voigt, 1.0, 1.4, 4);      // only signal
        // SetHistoQA(hphiMassProj);
        // fitFcn->SetParNames("Yield", "Mass", "Gwidth", "Lwidth", "p0", "p1", "p2");
        // hphiMassProj->GetYaxis()->SetTitleOffset(1.3);
        // hphiMassProj->GetYaxis()->SetMaxDigits(3);
        // hphiMassProj->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
        // hphiMassProj->GetXaxis()->SetNdivisions(505);
        // hphiMassProj->GetYaxis()->SetTitle("Counts");
        // hphiMassProj->SetMinimum(0);
        // hphiMassProj->Draw("pe");

        // // for voigtian distribution
        // fitFcn->SetParameter(0, 10000);     // yield
        // fitFcn->SetParLimits(0, 1000, 1e5); // yield
        // fitFcn->SetParameter(1, 1.019);     // mass peak
        // // if (ifiles == 4)
        // //     fitFcn->SetParLimits(1, 1.0064, 1.0316);
        // // else
        // fitFcn->SetParLimits(1, 1.006, 1.0316);

        // // fitFcn->SetParameter(3, 0.0042);   //lorentzian width
        // fitFcn->FixParameter(3, 0.0042);         // Lwidth
        // fitFcn->SetParameter(2, 0.0012);         //  Gwidth
        // fitFcn->SetParLimits(2, 0.001, 0.004);   // Gwidth. Double_t voigtpol2(Double_t *x, Double_t *par)
        // TFitResultPtr r;
        // r = hphiMassProj->Fit(fitFcn, "REBMS", "", 1.001, 1.039);

        // fitFcnBkg->SetParameters(fitFcn->GetParameter(4), fitFcn->GetParameter(5), fitFcn->GetParameter(6));
        // fitFcnSig->SetParameters(fitFcn->GetParameter(0), fitFcn->GetParameter(1), fitFcn->GetParameter(2), fitFcn->GetParameter(3));
        // fitFcnBkg->SetLineColor(kBlue);
        // fitFcnSig->SetLineColor(kGreen + 2);
        // fitFcnBkg->SetLineStyle(2);
        // fitFcnSig->SetLineStyle(2);
        // fitFcnSig->SetNpx(10000);
        // fitFcnBkg->Draw("same");
        // fitFcnSig->Draw("same");
        // hPhiMassFit->SetBinContent(ifiles + 1, fitFcn->GetParameter(1));
        // hPhiMassFit->SetBinError(ifiles + 1, fitFcn->GetParError(1));
        // hPhiMassFit->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        // double yieldIntegral = fitFcnSig->Integral(1.019 - 2 * 0.0042, 1.019 + 2 * 0.0042) / hevent_processed->GetBinContent(1);
        // TMatrixDSym cov = r->GetCovarianceMatrix();
        // TMatrixDSym cov1;
        // TMatrixDSym cov2;
        // cov.GetSub(0, 2, 0, 2, cov1);
        // cov.GetSub(3, 6, 3, 6, cov2);
        // Double_t *b = cov1.GetMatrixArray();
        // Double_t *a = cov2.GetMatrixArray();
        // Double_t *para = fitFcn->GetParameters();
        // double yieldIntegralError = fitFcnSig->IntegralError(1.001, 1.039, &para[0], b);

        // hPhiYieldFit->SetBinContent(ifiles + 1, yieldIntegral);
        // hPhiYieldFit->SetBinError(ifiles + 1, yieldIntegralError / hevent_processed->GetBinContent(1));
        // hPhiYieldFit->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        // hPhiMassResolution->SetBinContent(ifiles + 1, fitFcn->GetParameter(2));
        // hPhiMassResolution->SetBinError(ifiles + 1, fitFcn->GetParError(2));
        // hPhiMassResolution->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        // hChi2byNDF->SetBinContent(ifiles + 1, fitFcn->GetChisquare() / fitFcn->GetNDF());
        // hChi2byNDF->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        c3->cd(ifiles + 1);
        int lowbinpT2 = hdoublePhiMass->GetYaxis()->FindBin(2.0 + 0.001);
        int highbinpT2 = hdoublePhiMass->GetYaxis()->FindBin(10.0 - 0.001);
        TH1D *hdoublePhiMassProj = hdoublePhiMass->ProjectionX(Form("hdoublePhiMassProj_%d", ifiles), lowbinpT2, highbinpT2);
        SetHistoQA(hdoublePhiMassProj);
        hdoublePhiMassProj->GetYaxis()->SetTitleOffset(1.3);
        hdoublePhiMassProj->GetYaxis()->SetMaxDigits(3);
        hdoublePhiMassProj->GetXaxis()->SetTitle("#it{M}_{#Phi#Phi} (GeV/#it{c}^{2})");
        hdoublePhiMassProj->GetXaxis()->SetNdivisions(505);
        hdoublePhiMassProj->Rebin(10);
        hdoublePhiMassProj->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", (hdoublePhiMassProj->GetXaxis()->GetBinWidth(1) * 1000)));
        hdoublePhiMassProj->GetXaxis()->SetRangeUser(2.6, 2.9);
        hdoublePhiMassProj->GetXaxis()->SetNdivisions(508);
        hdoublePhiMassProj->Draw("pe");

        // c4->cd(ifiles + 1);
        // // SetHistoQA(hNsigmaTPCKaon);
        // hNsigmaTPCKaon->SetTitle(0);
        // hNsigmaTPCKaon->GetXaxis()->SetTitle("#it{n}_{#sigma} (TPC)");
        // hNsigmaTPCKaon->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        // hNsigmaTPCKaon->GetYaxis()->SetTitleOffset(1.3);
        // hNsigmaTPCKaon->GetYaxis()->SetMaxDigits(3);
        // hNsigmaTPCKaon->GetXaxis()->SetNdivisions(505);
        // hNsigmaTPCKaon->GetXaxis()->SetRangeUser(-5.0, 5.0);
        // hNsigmaTPCKaon->GetZaxis()->SetMaxDigits(3);
        // hNsigmaTPCKaon->SetStats(0);
        // hNsigmaTPCKaon->Draw("colz");

        // c5->cd(ifiles + 1);
        // // SetHistoQA(hNsigmaTOFKaon);
        // hNsigmaTOFKaon->SetTitle(0);
        // hNsigmaTOFKaon->GetXaxis()->SetTitle("#it{n}_{#sigma} (TOF)");
        // hNsigmaTOFKaon->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        // hNsigmaTOFKaon->GetYaxis()->SetTitleOffset(1.3);
        // hNsigmaTOFKaon->GetYaxis()->SetMaxDigits(3);
        // hNsigmaTOFKaon->GetXaxis()->SetNdivisions(505);
        // hNsigmaTOFKaon->GetXaxis()->SetRangeUser(-5.0, 5.0);
        // hNsigmaTOFKaon->GetZaxis()->SetMaxDigits(3);
        // hNsigmaTOFKaon->SetStats(0);
        // hNsigmaTOFKaon->Draw("colz");

        mainPad6->cd(ifiles + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.02);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.07);
        int lowbinpt_nsigma = hNsigmaTPCKaon->GetYaxis()->FindBin(0.5 + 0.001);
        int highbinpt_nsigma = hNsigmaTPCKaon->GetYaxis()->FindBin(5.0 - 0.001);
        TH1F *hNsigmaTPC = (TH1F *)hNsigmaTPCKaon->ProjectionX(Form("hNsigmaTPC_%d", ifiles), lowbinpt_nsigma, highbinpt_nsigma);
        SetHistoQA(hNsigmaTPC);
        hNsigmaTPC->SetTitle(Form("Run %s", runNumbers[ifiles].c_str()));
        hNsigmaTPC->GetXaxis()->SetTitle("#it{n}_{#sigma} (TPC)");
        hNsigmaTPC->GetYaxis()->SetTitle("Counts");
        hNsigmaTPC->Draw();
        // Draw a vertical line at x=0
        TLine *line = new TLine(0, 0, 0, hNsigmaTPC->GetMaximum());
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw("same");
        lat.DrawLatex(0.63, 0.7, Form("#mu: %.3f", hNsigmaTPC->GetMean()));
        lat.DrawLatex(0.63, 0.6, Form("#sigma: %.3f", hNsigmaTPC->GetRMS()));

        mainPad7->cd(ifiles + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.02);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.07);
        TH1F *hNsigmaTOF = (TH1F *)hNsigmaTOFKaon->ProjectionX(Form("hNsigmaTOF_%d", ifiles), lowbinpt_nsigma, highbinpt_nsigma);
        SetHistoQA(hNsigmaTOF);
        hNsigmaTOF->SetTitle(Form("Run %s", runNumbers[ifiles].c_str()));
        hNsigmaTOF->GetXaxis()->SetTitle("#it{n}_{#sigma} (TOF)");
        hNsigmaTOF->GetYaxis()->SetTitle("Counts");
        hNsigmaTOF->Draw();
        TLine *line2 = new TLine(0, 0, 0, hNsigmaTOF->GetMaximum());
        line2->SetLineColor(kRed);
        line2->SetLineStyle(2);
        line2->Draw("same");
        lat.DrawLatex(0.65, 0.7, Form("#mu: %.3f", hNsigmaTOF->GetMean()));
        lat.DrawLatex(0.65, 0.6, Form("#sigma: %.3f", hNsigmaTOF->GetRMS()));

        // */
    }

    c6->cd();
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
    latCanvasTitle6.DrawLatex(0.5, 0.5, "nSigma TPC Kaon");
    c6->Modified();
    c6->Update();
    c6->SaveAs("doublePhi/TriggerPlots/TPC_Kaon.png");

    c7->cd();
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
    latCanvasTitle7.DrawLatex(0.5, 0.5, "nSigma TOF Kaon");
    c7->Modified();
    c7->Update();
    c7->SaveAs("doublePhi/TriggerPlots/TOF_Kaon.png");

    TCanvas *cSelectivityDoublePhi = new TCanvas("cSelectivityDoublePhi", "Selectivity", 720, 720);
    gPad->SetGrid(1, 1);
    SetCanvasStyle(cSelectivityDoublePhi, 0.13, 0.03, 0.05, 0.14);
    SetHistoQA(hselectivity);
    hselectivity->SetMarkerSize(1.4);
    hselectivity->GetYaxis()->SetTitle("#Phi#Phi Selectivity");
    hselectivity->GetYaxis()->SetTitleOffset(1.3);
    hselectivity->GetXaxis()->SetTitleOffset(1.4);
    hselectivity->GetYaxis()->SetMaxDigits(3);
    hselectivity->SetMaximum(hselectivity->GetMaximum() * 1.5);
    hselectivity->SetMinimum(hselectivity->GetMinimum() * 0.5);
    hselectivity->Draw("p");
    cSelectivityDoublePhi->SaveAs("doublePhi/TriggerPlots/selectivity.png");

    // // /*
    // TCanvas *cPhiMass = new TCanvas("cPhiMass", "Phi Mass Peak", 720, 720);
    // gPad->SetGrid(1, 1);
    // SetCanvasStyle(cPhiMass, 0.19, 0.07, 0.05, 0.16);
    // SetHistoQA(hPhiMassFit);
    // hPhiMassFit->SetMarkerSize(1.4);
    // hPhiMassFit->GetYaxis()->SetTitle("#Phi Mass Peak (GeV/#it{c}^{2})");
    // hPhiMassFit->GetXaxis()->SetTitle("Run Number");
    // hPhiMassFit->GetYaxis()->SetTitleOffset(2.0);
    // hPhiMassFit->GetXaxis()->SetTitleOffset(1.5);
    // hPhiMassFit->GetYaxis()->SetMaxDigits(3);
    // hPhiMassFit->GetYaxis()->SetRangeUser(1.018, 1.0205);
    // hPhiMassFit->Draw("pe");
    // // TLine *linePhiPDGMass = new TLine(0, 1.019461, totalFiles + 1, 1.019461);
    // // linePhiPDGMass->SetLineColor(kRed);
    // // linePhiPDGMass->SetLineStyle(2);
    // // linePhiPDGMass->Draw("same");
    // TLegend *legend2 = new TLegend(0.5, 0.80, 0.85, 0.95);
    // // legend2->AddEntry((TObject *)0, "LHC25ah", "");
    // legend2->AddEntry(hPhiMassFit, "LHC25ah", "p");
    // // legend2->AddEntry(hPhiMassFit_skimmed, "LHC25ah_skimmed", "p");
    // // legend2->AddEntry(linePhiPDGMass, "#Phi PDG Mass", "l");
    // legend2->SetTextFont(42);
    // legend2->SetBorderSize(0);
    // legend2->SetFillStyle(0);
    // legend2->SetTextSize(0.035);
    // legend2->Draw();
    // // cPhiMass->SaveAs("doublePhi/phi_mass_peak.png");

    // TCanvas *cPhiMassResolution = new TCanvas("cPhiMassResolution", "Phi Mass Resolution", 720, 720);
    // gPad->SetGrid(1, 1);
    // SetCanvasStyle(cPhiMassResolution, 0.19, 0.07, 0.05, 0.16);
    // SetHistoQA(hPhiMassResolution);
    // hPhiMassResolution->SetMarkerSize(1.4);
    // hPhiMassResolution->GetYaxis()->SetTitle("#Phi Mass Resolution (GeV/#it{c}^{2})");
    // hPhiMassResolution->GetXaxis()->SetTitle("Run Number");
    // hPhiMassResolution->GetYaxis()->SetTitleOffset(2.0);
    // hPhiMassResolution->GetXaxis()->SetTitleOffset(1.5);
    // hPhiMassResolution->GetYaxis()->SetMaxDigits(3);
    // hPhiMassResolution->GetYaxis()->SetRangeUser(1e-3, 2e-3);
    // hPhiMassResolution->Draw("pe");
    // legend2->Draw();
    // // cPhiMassResolution->SaveAs("doublePhi/phi_mass_resolution.png");

    // TCanvas *cPhiYield = new TCanvas("cPhiYield", "Phi Yield", 720, 720);
    // gPad->SetGrid(1, 1);
    // SetCanvasStyle(cPhiYield, 0.19, 0.07, 0.05, 0.16);
    // SetHistoQA(hPhiYieldFit);
    // hPhiYieldFit->SetMarkerSize(1.4);
    // hPhiYieldFit->GetYaxis()->SetTitle("#Phi normalized raw yield");
    // hPhiYieldFit->GetXaxis()->SetTitle("Run Number");
    // hPhiYieldFit->GetYaxis()->SetTitleOffset(2.0);
    // hPhiYieldFit->GetXaxis()->SetTitleOffset(1.5);
    // hPhiYieldFit->GetYaxis()->SetMaxDigits(3);
    // // hPhiYieldFit->GetYaxis()->SetRangeUser(2000, 6000);
    // // hPhiYieldFit->GetYaxis()->SetRangeUser(2e-6, 3.2e-6);
    // hPhiYieldFit->Draw("pe");
    // legend2->Draw();
    // // cPhiYield->SaveAs("doublePhi/phi_yield.png");

    // TCanvas *cChi2byNDF = new TCanvas("cChi2byNDF", "Chi2 by NDF", 720, 720);
    // gPad->SetGrid(1, 1);
    // SetCanvasStyle(cChi2byNDF, 0.19, 0.07, 0.05, 0.16);
    // SetHistoQA(hChi2byNDF);
    // hChi2byNDF->SetMarkerSize(1.4);
    // hChi2byNDF->GetYaxis()->SetTitle("#chi^{2}/NDF");
    // hChi2byNDF->GetXaxis()->SetTitle("Run Number");
    // hChi2byNDF->GetYaxis()->SetTitleOffset(2.0);
    // hChi2byNDF->GetXaxis()->SetTitleOffset(1.5);
    // hChi2byNDF->GetYaxis()->SetMaxDigits(3);
    // // hChi2byNDF->GetYaxis()->SetRangeUser(0.5, 1.5);
    // hChi2byNDF->Draw("p");
    // // cChi2byNDF->SaveAs("doublePhi/chi2_by_ndf.png");
    // // */
}