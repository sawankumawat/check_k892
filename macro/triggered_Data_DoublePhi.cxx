#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"
using namespace std;

void triggered_Data_DoublePhi()
{
    gStyle->SetOptStat(0);
    gStyle->SetGridColor(15);

    // *******************LHC25ah*********************
    string runNumbers[] = {"564704", "564746", "564716", "564747", "564717", "564762", "564734", "564763", "564735", "564775"};
    string path = "/home/sawan/check_k892/data/doublePhi/TriggeredData2/LHC25ah/";
    string path_skimmed = "/home/sawan/check_k892/data/doublePhi/TriggeredData2/LHC25ah_Skimmed/";
    int totalFiles = sizeof(runNumbers) / sizeof(runNumbers[0]);

    TH1F *hTriggerEfficiencyPhi = new TH1F("hTriggerEfficiencyPhi", "Trigger Efficiency Phi;Run Number;Trigger Efficiency (#Phi)", totalFiles, 1, totalFiles + 1);
    TH1F *hTriggerEfficiencyDoublePhi = new TH1F("hTriggerEfficiencyDoublePhi", "Trigger Efficiency Double Phi;Run Number;Trigger Efficiency (Double #Phi)", totalFiles, 1, totalFiles + 1);

    TH1F *hTriggerEfficiencyPhi_skimmed = new TH1F("hTriggerEfficiencyPhi_skimmed", "Trigger Efficiency Phi (skimmed);Run Number;Trigger Efficiency (#Phi)", totalFiles, 1, totalFiles + 1);
    TH1F *hTriggerEfficiencyDoublePhi_skimmed = new TH1F("hTriggerEfficiencyDoublePhi_skimmed", "Trigger Efficiency Double Phi (skimmed);Run Number;Trigger Efficiency (Double #Phi)", totalFiles, 1, totalFiles + 1);

    TH1F *hselectivityDoublePhi = new TH1F("hselectivityDoublePhi", "Selectivity;Run Number;Selectivity", totalFiles, 1, totalFiles + 1); // Ratio of double phi triggered events in skimmed vs unskimmed data
    TH1F *hselectivityPhi = new TH1F("hselectivityPhi", "Selectivity;Run Number;Selectivity", totalFiles, 1, totalFiles + 1);             // Ratio of phi triggered events in skimmed vs unskimmed data

    TH1F *hPhiMassFit = new TH1F("hPhiMassFit", "", totalFiles, 1, totalFiles + 1);
    TH1F *hPhiMassFit_skimmed = new TH1F("hPhiMassFit_skimmed", "", totalFiles, 1, totalFiles + 1);
    TH1F *hPhiYieldFit = new TH1F("hPhiYieldFit", "", totalFiles, 1, totalFiles + 1);
    TH1F *hPhiYieldFit_skimmed = new TH1F("hPhiYieldFit_skimmed", "", totalFiles, 1, totalFiles + 1);
    TH1F *hPhiMassResolution = new TH1F("hPhiMassResolution", "", totalFiles, 1, totalFiles + 1);
    TH1F *hPhiMassResolution_skimmed = new TH1F("hPhiMassResolution_skimmed", "", totalFiles, 1, totalFiles + 1);
    TH1F *hChi2byNDF = new TH1F("hChi2byNDF", "", totalFiles, 1, totalFiles + 1);

    //******Fitting for Phi*********************
    TF1 *fitFcn = new TF1("fitfunc", voigtpol2, 1.0, 1.4, 7);       // sig+bkg fit function
    TF1 *fitFcnBkg = new TF1("fitfunc1", polynomial2, 1.0, 1.4, 3); // only residualbkg
    TF1 *fitFcnSig = new TF1("fitFcnSig", voigt, 1.0, 1.4, 4);      // only signal

    for (int ifiles = 0; ifiles < totalFiles; ifiles++)
    // for (int ifiles = 0; ifiles < 2; ifiles++)
    {
        TFile *file = new TFile((path + "AnalysisResults_" + runNumbers[ifiles] + ".root").c_str(), "READ");
        TFile *file_skimmed = new TFile((path_skimmed + "AnalysisResults_" + runNumbers[ifiles] + ".root").c_str(), "READ");
        if (file->IsZombie() || file_skimmed->IsZombie())
        {
            std::cerr << "Error: Could not open file triggered_data/AnalysisResults_fullrun.root\n";
            return;
        }

        TH1F *hevent_processed = (TH1F *)file->Get("lf-doublephi-filter/hProcessedEvents");
        TH2F *hphiMass = (TH2F *)file->Get("lf-doublephi-filter/hInvMassPhi");
        TH2F *hdoublePhiMass = (TH2F *)file->Get("lf-doublephi-filter/hInvMassDoublePhi");
        TH2F *hNsigmaTPCKaon = (TH2F *)file->Get("lf-doublephi-filter/hNsigmaPtkaonTPC");
        TH2F *hNsigmaTOFKaon = (TH2F *)file->Get("lf-doublephi-filter/hNsigmaPtkaonTOF");

        TH1F *hevent_processed_skimmed = (TH1F *)file_skimmed->Get("lf-doublephi-filter/hProcessedEvents");
        TH2F *hphiMass_skimmed = (TH2F *)file_skimmed->Get("lf-doublephi-filter/hInvMassPhi");

        if (hevent_processed == nullptr || hphiMass == nullptr || hdoublePhiMass == nullptr || hNsigmaTPCKaon == nullptr || hNsigmaTOFKaon == nullptr)
        {
            std::cerr << "Error: Could not find one or more histograms in the file\n";
            return;
        }

        hTriggerEfficiencyPhi->SetBinContent(ifiles + 1, hevent_processed->GetBinContent(2) / hevent_processed->GetBinContent(1));
        hTriggerEfficiencyDoublePhi->SetBinContent(ifiles + 1, hevent_processed->GetBinContent(3) / hevent_processed->GetBinContent(1));
        hTriggerEfficiencyPhi->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());
        hTriggerEfficiencyDoublePhi->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        hTriggerEfficiencyPhi_skimmed->SetBinContent(ifiles + 1, hevent_processed_skimmed->GetBinContent(2) / hevent_processed_skimmed->GetBinContent(1));
        hTriggerEfficiencyDoublePhi_skimmed->SetBinContent(ifiles + 1, hevent_processed_skimmed->GetBinContent(3) / hevent_processed_skimmed->GetBinContent(1));

        hselectivityDoublePhi->SetBinContent(ifiles + 1, hevent_processed_skimmed->GetBinContent(3) / hevent_processed->GetBinContent(3));
        hselectivityPhi->SetBinContent(ifiles + 1, hevent_processed_skimmed->GetBinContent(2) / hevent_processed->GetBinContent(2));

        // TCanvas *c1 = new TCanvas("c1", "Processed Events", 1080, 720);
        // c1->SetGrid(1, 1);
        // SetCanvasStyle(c1, 0.11, 0.06, 0.05, 0.10);
        // gPad->SetLogy();
        // SetHistoQA(hevent_processed);
        // hevent_processed->GetYaxis()->SetTitleOffset(1.1);
        // hevent_processed->Draw();
        // c1->SaveAs(Form("doublePhi/allPlots/processed_events_%s.png", runNumbers[ifiles].c_str()));

        // /*

        TCanvas *c2 = new TCanvas("c2", "Phi Mass", 720, 720);
        SetCanvasStyle(c2, 0.13, 0.06, 0.05, 0.13);
        int lowbinpT = hphiMass->GetYaxis()->FindBin(0.0 + 0.001);
        int highbinpT = hphiMass->GetYaxis()->FindBin(10.0 - 0.001);
        TH1D *hphiMassProj = hphiMass->ProjectionX("hphiMassProj", lowbinpT, highbinpT);
        TH1D *hphiMassProj_skimmed = hphiMass_skimmed->ProjectionX("hphiMassProj_skimmed", lowbinpT, highbinpT);
        SetHistoQA(hphiMassProj_skimmed);
        SetHistoQA(hphiMassProj);
        hphiMassProj->GetYaxis()->SetTitleOffset(1.3);
        hphiMassProj->GetYaxis()->SetMaxDigits(3);
        hphiMassProj->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
        hphiMassProj->GetXaxis()->SetNdivisions(505);
        hphiMassProj->GetYaxis()->SetTitle("Counts");
        hphiMassProj->SetMinimum(0);
        hphiMassProj->Draw("pe");

        // for voigtian distribution
        fitFcn->SetParameter(0, 5000);   // yield
        fitFcn->SetParLimits(0, 0, 1e5); // yield
        fitFcn->SetParameter(1, 1.019);  // mass peak
        if (ifiles == 4)
            fitFcn->SetParLimits(1, 1.0064, 1.0316);
        else
            fitFcn->SetParLimits(1, 1.006, 1.0316);

        // fitFcn->SetParameter(3, 0.0042);   //lorentzian width
        fitFcn->FixParameter(3, 0.0042);       // Lwidth
        fitFcn->SetParameter(2, 0.0012);       //  Gwidth
        fitFcn->SetParLimits(2, 0.001, 0.009); // Gwidth. Double_t voigtpol2(Double_t *x, Double_t *par)
        TFitResultPtr r = hphiMassProj->Fit("fitfunc", "REMS", "", 1.001, 1.039);
        fitFcnBkg->SetParameters(fitFcn->GetParameter(4), fitFcn->GetParameter(5), fitFcn->GetParameter(6));
        fitFcnSig->SetParameters(fitFcn->GetParameter(0), fitFcn->GetParameter(1), fitFcn->GetParameter(2), fitFcn->GetParameter(3));
        fitFcnBkg->SetLineColor(kBlue);
        fitFcnSig->SetLineColor(kGreen + 2);
        fitFcnBkg->SetLineStyle(2);
        fitFcnSig->SetLineStyle(2);
        fitFcnSig->SetNpx(10000);
        fitFcnBkg->Draw("same");
        fitFcnSig->Draw("same");
        hPhiMassFit->SetBinContent(ifiles + 1, fitFcn->GetParameter(1));
        hPhiMassFit->SetBinError(ifiles + 1, fitFcn->GetParError(1));
        hPhiMassFit->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        double yieldIntegral = fitFcnSig->Integral(1.019 - 2 * 0.0042, 1.019 + 2 * 0.0042) / hevent_processed->GetBinContent(1);
        TMatrixDSym cov = r->GetCovarianceMatrix();
        TMatrixDSym cov1;
        TMatrixDSym cov2;
        cov.GetSub(0, 2, 0, 2, cov1);
        cov.GetSub(3, 6, 3, 6, cov2);
        Double_t *b = cov1.GetMatrixArray();
        Double_t *a = cov2.GetMatrixArray();
        Double_t *para = fitFcn->GetParameters();
        double yieldIntegralError = fitFcnSig->IntegralError(1.001, 1.039, &para[0], b);

        hPhiYieldFit->SetBinContent(ifiles + 1, yieldIntegral);
        hPhiYieldFit->SetBinError(ifiles + 1, yieldIntegralError / hevent_processed->GetBinContent(1));
        hPhiYieldFit->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        hPhiMassResolution->SetBinContent(ifiles + 1, fitFcn->GetParameter(2));
        hPhiMassResolution->SetBinError(ifiles + 1, fitFcn->GetParError(2));
        hPhiMassResolution->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        cout << "Unskimmed data yield is " << fitFcnSig->Integral(1.019 - 2 * 0.0042, 1.019 + 2 * 0.0042) << endl;
        cout << "Unskimmed data triggered events " << hevent_processed->GetBinContent(2) << endl;
        cout << "Normalized yield " << yieldIntegral << endl;

        hChi2byNDF->SetBinContent(ifiles + 1, fitFcn->GetChisquare() / fitFcn->GetNDF());
        hChi2byNDF->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());
        c2->SaveAs(Form("doublePhi/allPlots/phi_mass_%s.png", runNumbers[ifiles].c_str()));

        TCanvas *c22 = new TCanvas("c22", "Phi Mass Skimmed", 720, 720);
        SetCanvasStyle(c22, 0.13, 0.06, 0.05, 0.13);
        SetHistoQA(hphiMassProj_skimmed);
        hphiMassProj_skimmed->GetYaxis()->SetTitleOffset(1.3);
        hphiMassProj_skimmed->GetYaxis()->SetMaxDigits(3);
        hphiMassProj_skimmed->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
        hphiMassProj_skimmed->GetXaxis()->SetNdivisions(505);
        hphiMassProj_skimmed->GetYaxis()->SetTitle("Counts");
        hphiMassProj_skimmed->SetMinimum(0);
        hphiMassProj_skimmed->Draw("pe");

        // for voigtian distribution
        fitFcn->SetParameter(0, 1000);   // yield
        fitFcn->SetParLimits(0, 0, 1e5); // yield
        fitFcn->SetParameter(1, 1.019);  // mass peak
        if (ifiles == 4)
            fitFcn->SetParLimits(1, 1.0064, 1.0316);
        else
            fitFcn->SetParLimits(1, 1.005, 1.0316);

        fitFcn->SetParameter(3, 0.0042);       // lorentzian width
        fitFcn->FixParameter(3, 0.0042);       // Lwidth
        fitFcn->SetParameter(2, 0.0012);       //  Gwidth
        fitFcn->SetParLimits(2, 0.001, 0.015); // Gwidth. Double_t voigtpol2(Double_t *x, Double_t *par)

        TFitResultPtr r1 = hphiMassProj_skimmed->Fit("fitfunc", "REMS", "", 1.001, 1.039);
        fitFcnBkg->SetParameters(fitFcn->GetParameter(4), fitFcn->GetParameter(5), fitFcn->GetParameter(6));
        fitFcnSig->SetParameters(fitFcn->GetParameter(0), fitFcn->GetParameter(1), fitFcn->GetParameter(2), fitFcn->GetParameter(3));
        fitFcnBkg->SetLineColor(kBlue);
        fitFcnSig->SetLineColor(kGreen + 2);
        fitFcnBkg->SetLineStyle(2);
        fitFcnSig->SetLineStyle(2);
        fitFcnSig->SetNpx(10000);
        fitFcnBkg->Draw("same");
        fitFcnSig->Draw("same");
        hPhiMassFit_skimmed->SetBinContent(ifiles + 1, fitFcn->GetParameter(1));
        hPhiMassFit_skimmed->SetBinError(ifiles + 1, fitFcn->GetParError(1));
        hPhiMassFit_skimmed->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        double yieldIntegral_skimmed = fitFcnSig->Integral(1.019 - 2 * 0.0042, 1.019 + 2 * 0.0042) / hevent_processed_skimmed->GetBinContent(1);
        TMatrixDSym cov_skimmed = r1->GetCovarianceMatrix();
        TMatrixDSym cov1_skimmed;
        TMatrixDSym cov2_skimmed;
        cov_skimmed.GetSub(0, 2, 0, 2, cov1_skimmed);
        cov_skimmed.GetSub(3, 6, 3, 6, cov2_skimmed);
        Double_t *b_skimmed = cov1_skimmed.GetMatrixArray();
        Double_t *a_skimmed = cov2_skimmed.GetMatrixArray();
        Double_t *para_skimmed = fitFcn->GetParameters();
        double yieldIntegralError_skimmed = fitFcnSig->IntegralError(1.001, 1.039, &para_skimmed[0], b_skimmed);

        cout << "skimmed data yield is " << fitFcnSig->Integral(1.019 - 2 * 0.0042, 1.019 + 2 * 0.0042) << endl;
        cout << "skimmed data triggered events " << hevent_processed_skimmed->GetBinContent(2) << endl;
        cout << "Normalized yield skimmed " << yieldIntegral_skimmed << endl;

        hPhiYieldFit_skimmed->SetBinContent(ifiles + 1, yieldIntegral_skimmed);
        hPhiYieldFit_skimmed->SetBinError(ifiles + 1, yieldIntegralError / (hevent_processed_skimmed->GetBinContent(1) * 50)); // FIXME: temporary fix of 50 as the skimmed data has 50 times less events
        // hPhiYieldFit_skimmed->SetBinError(ifiles + 1, 0);
        hPhiYieldFit_skimmed->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());

        hPhiMassResolution_skimmed->SetBinContent(ifiles + 1, fitFcn->GetParameter(2));
        hPhiMassResolution_skimmed->SetBinError(ifiles + 1, fitFcn->GetParError(2));
        hPhiMassResolution_skimmed->GetXaxis()->SetBinLabel(ifiles + 1, runNumbers[ifiles].c_str());
        c22->SaveAs(Form("doublePhi/allPlots/phi_mass_skimmed_%s.png", runNumbers[ifiles].c_str()));

        // */

        // TCanvas *c3 = new TCanvas("c3", "Double Phi Mass", 720, 720);
        // SetCanvasStyle(c3, 0.13, 0.06, 0.05, 0.13);
        // int lowbinpT2 = hdoublePhiMass->GetYaxis()->FindBin(2.0 + 0.001);
        // int highbinpT2 = hdoublePhiMass->GetYaxis()->FindBin(10.0 - 0.001);
        // TH1D *hdoublePhiMassProj = hdoublePhiMass->ProjectionX("hdoublePhiMassProj", lowbinpT2, highbinpT2);
        // SetHistoQA(hdoublePhiMassProj);
        // hdoublePhiMassProj->GetYaxis()->SetTitleOffset(1.3);
        // hdoublePhiMassProj->GetYaxis()->SetMaxDigits(3);
        // hdoublePhiMassProj->GetXaxis()->SetTitle("#it{M}_{#Phi#Phi} (GeV/#it{c}^{2})");
        // hdoublePhiMassProj->GetXaxis()->SetNdivisions(505);
        // hdoublePhiMassProj->Rebin(10);
        // hdoublePhiMassProj->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", (hdoublePhiMassProj->GetXaxis()->GetBinWidth(1) * 1000)));
        // hdoublePhiMassProj->GetXaxis()->SetRangeUser(2.6, 2.9);
        // hdoublePhiMassProj->GetXaxis()->SetNdivisions(508);
        // hdoublePhiMassProj->Draw("pe");
        // c3->SaveAs(Form("doublePhi/allPlots/double_phi_mass_pt0to10_%s.png", runNumbers[ifiles].c_str()));

        // TCanvas *c4 = new TCanvas("c4", "Nsigma TPC Kaon", 720, 720);
        // SetCanvasStyle(c4, 0.13, 0.14, 0.05, 0.13);
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
        // c4->SaveAs(Form("doublePhi/allPlots/nsigma_tpc_kaon_%s.png", runNumbers[ifiles].c_str()));

        // TCanvas *c5 = new TCanvas("c5", "Nsigma TOF Kaon", 720, 720);
        // SetCanvasStyle(c5, 0.13, 0.14, 0.05, 0.13);
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
        // c5->SaveAs(Form("doublePhi/allPlots/nsigma_tof_kaon_%s.png", runNumbers[ifiles].c_str()));

        // TCanvas *c6 = new TCanvas("c6", "Nsigma TPC", 720, 720);
        // SetCanvasStyle(c6, 0.13, 0.14, 0.05, 0.13);
        // int lowbinpt_nsigma = hNsigmaTPCKaon->GetYaxis()->FindBin(0.5 + 0.001);
        // int highbinpt_nsigma = hNsigmaTPCKaon->GetYaxis()->FindBin(5.0 - 0.001);
        // TH1F *hNsigmaTPC = (TH1F *)hNsigmaTPCKaon->ProjectionX("hNsigmaTPC", lowbinpt_nsigma, highbinpt_nsigma);
        // hNsigmaTPC->Draw();
        // // Draw a vertical line at x=0
        // TLine *line = new TLine(0, 0, 0, hNsigmaTPC->GetMaximum());
        // line->SetLineColor(kRed);
        // line->SetLineStyle(2);
        // line->Draw("same");
        // c6->SaveAs(Form("doublePhi/allPlots/nsigma_tpc1D_%s.png", runNumbers[ifiles].c_str()));

        // TCanvas *c7 = new TCanvas("c7", "Nsigma TOF", 720, 720);
        // SetCanvasStyle(c7, 0.13, 0.14, 0.05, 0.13);
        // TH1F *hNsigmaTOF = (TH1F *)hNsigmaTOFKaon->ProjectionX("hNsigmaTOF", lowbinpt_nsigma, highbinpt_nsigma);
        // hNsigmaTOF->Draw();
        // TLine *line2 = new TLine(0, 0, 0, hNsigmaTOF->GetMaximum());
        // line2->SetLineColor(kRed);
        // line2->SetLineStyle(2);
        // line2->Draw("same");
        // c7->SaveAs(Form("doublePhi/allPlots/nsigma_tof1D_%s.png", runNumbers[ifiles].c_str()));
    }

    TCanvas *cSelectivityDoublePhi = new TCanvas("cSelectivityDoublePhi", "Selectivity", 720, 720);
    gPad->SetGrid(1, 1);
    SetCanvasStyle(cSelectivityDoublePhi, 0.19, 0.07, 0.05, 0.16);
    SetHistoQA(hselectivityDoublePhi);
    hselectivityDoublePhi->SetMarkerSize(1.4);
    hselectivityDoublePhi->GetYaxis()->SetTitle("#Phi#Phi triggered events (skimmed / unskimmed)");
    hselectivityDoublePhi->GetYaxis()->SetTitleOffset(2.0);
    hselectivityDoublePhi->GetXaxis()->SetTitleOffset(1.5);
    hselectivityDoublePhi->GetYaxis()->SetMaxDigits(3);
    // hselectivityDoublePhi->GetYaxis()->SetRangeUser(0.22, 0.32);
    hselectivityDoublePhi->GetYaxis()->SetRangeUser(0.7, 1.1);
    hselectivityDoublePhi->Draw("p");
    cSelectivityDoublePhi->SaveAs("doublePhi/selectivityDoublePhi.png");

    TCanvas *cSelectivityPhi = new TCanvas("cSelectivityPhi", "Selectivity", 720, 720);
    gPad->SetGrid(1, 1);
    SetCanvasStyle(cSelectivityPhi, 0.19, 0.07, 0.05, 0.16);
    SetHistoQA(hselectivityPhi);
    hselectivityPhi->SetMarkerSize(1.4);
    hselectivityPhi->GetYaxis()->SetTitle("#Phi triggered events (skimmed / unskimmed)");
    hselectivityPhi->GetYaxis()->SetTitleOffset(2.0);
    hselectivityPhi->GetXaxis()->SetTitleOffset(1.5);
    hselectivityPhi->GetYaxis()->SetMaxDigits(3);
    // hselectivityPhi->GetYaxis()->SetRangeUser(0.55, 0.75);
    hselectivityPhi->GetYaxis()->SetRangeUser(0.0, 0.2);
    hselectivityPhi->Draw("p");
    cSelectivityPhi->SaveAs("doublePhi/selectivityPhi.png");

    TCanvas *cEfficiencyPhi = new TCanvas("cEfficiencyPhi", "Trigger Efficiency Phi", 720, 720);
    // gPad->SetLogy();
    gPad->SetGrid(1, 1);
    SetCanvasStyle(cEfficiencyPhi, 0.19, 0.07, 0.05, 0.16);
    SetHistoQA(hTriggerEfficiencyPhi);
    hTriggerEfficiencyPhi->SetMarkerSize(1.4);
    hTriggerEfficiencyPhi->GetYaxis()->SetTitle("#Phi triggered events / total events");
    hTriggerEfficiencyPhi->GetYaxis()->SetTitleOffset(2.0);
    hTriggerEfficiencyPhi->GetXaxis()->SetTitleOffset(1.5);
    hTriggerEfficiencyPhi->GetYaxis()->SetMaxDigits(3);
    hTriggerEfficiencyPhi->GetYaxis()->SetRangeUser(0.14e-3, 0.18e-3);
    // hTriggerEfficiencyPhi->GetYaxis()->SetRangeUser(1e-4, 2e-3);
    // hTriggerEfficiencyPhi->LabelsOption("v");
    hTriggerEfficiencyPhi->Draw("p");
    SetHistoQA(hTriggerEfficiencyPhi_skimmed);
    hTriggerEfficiencyPhi_skimmed->SetMarkerColor(kRed);
    hTriggerEfficiencyPhi_skimmed->SetMarkerStyle(21);
    hTriggerEfficiencyPhi_skimmed->SetMarkerSize(1.4);
    hTriggerEfficiencyPhi_skimmed->Draw("p same");
    TLegend *legend = new TLegend(0.5, 0.80, 0.85, 0.95);
    // legend->AddEntry((TObject *)0, "LHC25ah", "");
    legend->AddEntry(hTriggerEfficiencyPhi, "LHC25ah", "p");
    // legend->AddEntry(hTriggerEfficiencyPhi_skimmed, "LHC25ah_skimmed", "p");
    legend->SetTextFont(42);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->Draw();
    cEfficiencyPhi->SaveAs("doublePhi/trigger_efficiency_phi.png");

    TCanvas *cEfficiencyDoublePhi = new TCanvas("cEfficiencyDoublePhi", "Trigger Efficiency Double Phi", 720, 720);
    // gPad->SetLogy();
    gPad->SetGrid(1, 1);
    SetCanvasStyle(cEfficiencyDoublePhi, 0.19, 0.07, 0.05, 0.16);
    SetHistoQA(hTriggerEfficiencyDoublePhi);
    hTriggerEfficiencyDoublePhi->SetMarkerSize(1.4);
    hTriggerEfficiencyDoublePhi->GetYaxis()->SetTitle("#Phi#Phi triggered events / total Events");
    hTriggerEfficiencyDoublePhi->GetYaxis()->SetTitleOffset(2.0);
    hTriggerEfficiencyDoublePhi->GetXaxis()->SetTitleOffset(1.5);
    hTriggerEfficiencyDoublePhi->GetYaxis()->SetMaxDigits(3);
    hTriggerEfficiencyDoublePhi->GetYaxis()->SetRangeUser(14.62e-6, 16.89e-6);
    // hTriggerEfficiencyDoublePhi->GetYaxis()->SetRangeUser(1e-6, 1e-2);
    hTriggerEfficiencyDoublePhi->Draw("p");
    SetHistoQA(hTriggerEfficiencyDoublePhi_skimmed);
    hTriggerEfficiencyDoublePhi_skimmed->SetMarkerColor(kRed);
    hTriggerEfficiencyDoublePhi_skimmed->SetMarkerStyle(21);
    hTriggerEfficiencyDoublePhi_skimmed->SetMarkerSize(1.4);
    hTriggerEfficiencyDoublePhi_skimmed->Draw("p same");
    legend->Draw();
    cEfficiencyDoublePhi->SaveAs("doublePhi/trigger_efficiency_double_phi.png");

    TCanvas *cPhiMass = new TCanvas("cPhiMass", "Phi Mass Peak", 720, 720);
    gPad->SetGrid(1, 1);
    SetCanvasStyle(cPhiMass, 0.19, 0.07, 0.05, 0.16);
    SetHistoQA(hPhiMassFit);
    SetHistoQA(hPhiMassFit_skimmed);
    hPhiMassFit->SetMarkerSize(1.4);
    hPhiMassFit->GetYaxis()->SetTitle("#Phi Mass Peak (GeV/#it{c}^{2})");
    hPhiMassFit->GetXaxis()->SetTitle("Run Number");
    hPhiMassFit->GetYaxis()->SetTitleOffset(2.0);
    hPhiMassFit->GetXaxis()->SetTitleOffset(1.5);
    hPhiMassFit->GetYaxis()->SetMaxDigits(3);
    hPhiMassFit->GetYaxis()->SetRangeUser(1.018, 1.0205);
    hPhiMassFit->Draw("pe");
    hPhiMassFit_skimmed->SetMarkerColor(kBlue);
    hPhiMassFit_skimmed->SetMarkerStyle(25);
    hPhiMassFit_skimmed->SetMarkerSize(1.4);
    hPhiMassFit_skimmed->Draw("pe same");
    // TLine *linePhiPDGMass = new TLine(0, 1.019461, totalFiles + 1, 1.019461);
    // linePhiPDGMass->SetLineColor(kRed);
    // linePhiPDGMass->SetLineStyle(2);
    // linePhiPDGMass->Draw("same");
    TLegend *legend2 = new TLegend(0.5, 0.80, 0.85, 0.95);
    // legend2->AddEntry((TObject *)0, "LHC25ah", "");
    legend2->AddEntry(hPhiMassFit, "LHC25ah", "p");
    legend2->AddEntry(hPhiMassFit_skimmed, "LHC25ah_skimmed", "p");
    // legend2->AddEntry(linePhiPDGMass, "#Phi PDG Mass", "l");
    legend2->SetTextFont(42);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(0.035);
    legend2->Draw();
    cPhiMass->SaveAs("doublePhi/phi_mass_peak.png");

    TCanvas *cPhiMassResolution = new TCanvas("cPhiMassResolution", "Phi Mass Resolution", 720, 720);
    gPad->SetGrid(1, 1);
    SetCanvasStyle(cPhiMassResolution, 0.19, 0.07, 0.05, 0.16);
    SetHistoQA(hPhiMassResolution);
    SetHistoQA(hPhiMassResolution_skimmed);
    hPhiMassResolution->SetMarkerSize(1.4);
    hPhiMassResolution->GetYaxis()->SetTitle("#Phi Mass Resolution (GeV/#it{c}^{2})");
    hPhiMassResolution->GetXaxis()->SetTitle("Run Number");
    hPhiMassResolution->GetYaxis()->SetTitleOffset(2.0);
    hPhiMassResolution->GetXaxis()->SetTitleOffset(1.5);
    hPhiMassResolution->GetYaxis()->SetMaxDigits(3);
    hPhiMassResolution->GetYaxis()->SetRangeUser(1e-3, 2e-3);
    hPhiMassResolution->Draw("pe");
    hPhiMassResolution_skimmed->SetMarkerColor(kBlue);
    hPhiMassResolution_skimmed->SetMarkerStyle(25);
    hPhiMassResolution_skimmed->SetMarkerSize(1.4);
    hPhiMassResolution_skimmed->Draw("pe same");
    legend2->Draw();
    cPhiMassResolution->SaveAs("doublePhi/phi_mass_resolution.png");

    TCanvas *cPhiYield = new TCanvas("cPhiYield", "Phi Yield", 720, 720);
    gPad->SetGrid(1, 1);
    SetCanvasStyle(cPhiYield, 0.19, 0.07, 0.05, 0.16);
    SetHistoQA(hPhiYieldFit);
    SetHistoQA(hPhiYieldFit_skimmed);
    hPhiYieldFit->SetMarkerSize(1.4);
    hPhiYieldFit->GetYaxis()->SetTitle("#Phi normalized raw yield");
    hPhiYieldFit->GetXaxis()->SetTitle("Run Number");
    hPhiYieldFit->GetYaxis()->SetTitleOffset(2.0);
    hPhiYieldFit->GetXaxis()->SetTitleOffset(1.5);
    hPhiYieldFit->GetYaxis()->SetMaxDigits(3);
    // hPhiYieldFit->GetYaxis()->SetRangeUser(2000, 6000);
    hPhiYieldFit->GetYaxis()->SetRangeUser(2e-6, 3.2e-6);
    hPhiYieldFit->Draw("pe");
    hPhiYieldFit_skimmed->SetMarkerColor(kBlue);
    hPhiYieldFit_skimmed->SetMarkerStyle(25);
    hPhiYieldFit_skimmed->SetMarkerSize(1.4);
    hPhiYieldFit_skimmed->Draw("pe same");
    legend2->Draw();
    cPhiYield->SaveAs("doublePhi/phi_yield.png");

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
    // cChi2byNDF->SaveAs("doublePhi/chi2_by_ndf.png");
}