#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"
using namespace std;
Double_t pol2(double *x, double *par);
Double_t pol3(double *x, double *par);
Double_t BW_pol3(double *x, double *par);
Double_t BW_pol2(double *x, double *par);
Double_t relBW(double *x, double *par);

void doublephi()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    bool MixedEventExist = false;
    bool isQA = false;
    bool isSavePlots = true;
    bool isPhiQA = false;
    // TString dataFilePath = "../data/doublePhi/";
    TString dataFilePath = "/home/sawan/alice/practice/";
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.SetTextFont(42);

    // TString fileName = "merge"; // merged 23 and 24 dataset output for PID 0
    // TString fileName = "508501"; // LHC23 data (508501)
    // TString fileName = "508495"; // LHC24 data
    // TString fileName = "LHC25ac_allWagons"; // Triggered data (LHC25ac)
    // TString fileName = "Analysis541300PID4New"; // Triggered data (LHC25a5)
    // TString fileName = "Analysis549273PID4New"; // Triggered data (LHC25ah)
    TString fileName = "Analysis_ac_ah_MergedPID4ME"; // Triggered data (LHC25ah+ac)

    TString subWagon = "";
    // TString subWagon = "_LoosePID";
    // TString subWagon = "_DeepAngle";
    // TString subWagon = "_StrategyPID1";
    // TString subWagon = "_DeltaRDaughter2";
    // TString subWagon = "_DeltaRDaughter4";
    // TString subWagon = "_id37219"; //18, 19 (Sourav bhaiya train outputs)

    TString outputPath = "../output/doublePhi/" + fileName + "/";
    gSystem->Exec("mkdir -p " + outputPath);
    TString outputPathQA = outputPath + "QA/";
    gSystem->Exec("mkdir -p " + outputPathQA);
    TString outputPhiFits = outputPath + "PhiFits/";
    gSystem->Exec("mkdir -p " + outputPhiFits);

    TFile *file = new TFile(dataFilePath + fileName + ".root");
    TFile *file2 = new TFile("../output/doublePhi/Analysis541300PID6/PhiPhi.root"); // Very loose cut, for efficiency calculation
    if (file->IsZombie() || file2->IsZombie())
    {
        std::cerr << "Error: Could not open file " << "\n";
        return;
    }
    TFile *outputFile = new TFile(outputPath + "PhiPhi" + subWagon + ".root", "recreate");

    THnSparseF *hunlike = (THnSparseF *)file->Get("doublephimeson" + subWagon + "/SEMassUnlike"); // The axes are invariant mass, deltaRKaons(minimum), pT, deltaRPhi , MassDifference, PtCorrelation
    THnSparseF *hmixed;
    if (MixedEventExist)
        hmixed = (THnSparseF *)file->Get("doublephimeson" + subWagon + "/MEMassUnlike");

    if (hunlike == nullptr)
    {
        std::cerr << "Error: Could not find histogram 'doublephimeson" + subWagon + "/SEMassUnlike' in file\n";
        return;
    }
    TH1D *hDeltaRKaon = hunlike->Projection(1, "E");
    TH1D *hDeltaRPhiPhi = hunlike->Projection(3, "E");
    TH1D *hDeltaMassDifference = hunlike->Projection(4, "E");
    TH1D *hPtCorrelation = hunlike->Projection(5, "E");
    hDeltaRPhiPhi->Write("hDeltaRPhiPhi");
    hDeltaMassDifference->Write("hDeltaMassDifference");

    int deltaRKaonLow = hunlike->GetAxis(1)->FindBin(0.02 + 0.0001);
    int deltaRKaonHigh = hunlike->GetAxis(1)->FindBin(10.0 - 0.0001);

    int lowbinpT = hunlike->GetAxis(2)->FindBin(2.0 + 0.01);
    int highbinpT = hunlike->GetAxis(2)->FindBin(100.0 - 0.01);

    int deltaRPhiLow = hunlike->GetAxis(3)->FindBin(0.2 + 0.0001);
    int deltaRPhiHigh = hunlike->GetAxis(3)->FindBin(10.0 - 0.0001);

    int deltaMLow = hunlike->GetAxis(4)->FindBin(0.0 + 0.00001);
    int deltaMHigh = hunlike->GetAxis(4)->FindBin(0.02 - 0.00001); // 0.01 is good cut, but losing 60% statistics

    int ptCorrelationLow = hunlike->GetAxis(5)->FindBin(0.0 + 0.0001);
    int ptCorrelationHigh = hunlike->GetAxis(5)->FindBin(2.5 - 0.0001);

    hunlike->GetAxis(1)->SetRange(deltaRKaonLow, deltaRKaonHigh);
    hunlike->GetAxis(2)->SetRange(lowbinpT, highbinpT);
    hunlike->GetAxis(3)->SetRange(deltaRPhiLow, deltaRPhiHigh);
    hunlike->GetAxis(4)->SetRange(deltaMLow, deltaMHigh);
    hunlike->GetAxis(5)->SetRange(ptCorrelationLow, ptCorrelationHigh);

    if (MixedEventExist)
    {
        hmixed->GetAxis(1)->SetRange(deltaRKaonLow, deltaRKaonHigh);
        hmixed->GetAxis(2)->SetRange(lowbinpT, highbinpT);
        hmixed->GetAxis(3)->SetRange(deltaRPhiLow, deltaRPhiHigh);
        hmixed->GetAxis(4)->SetRange(deltaMLow, deltaMHigh);
        hmixed->GetAxis(5)->SetRange(ptCorrelationLow, ptCorrelationHigh);
    }

    TH1D *hmass = hunlike->Projection(0, "E");
    TH1D *hmassClone = (TH1D *)hmass->Clone();
    if (hmass == nullptr)
    {
        std::cerr << "Error: Could not create projection histogram\n";
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "#Phi#Phi invariant mass", 720, 720);
    SetCanvasStyle(c1, 0.17, 0.03, 0.05, 0.15);
    SetHistoQA(hmass);
    hmass->GetXaxis()->SetTitle("#it{M}_{#Phi#Phi} (GeV/#it{c}^{2})");
    hmass->SetMarkerStyle(20);
    hmass->GetYaxis()->SetMaxDigits(3);
    hmass->SetMarkerSize(1.0);
    hmass->GetYaxis()->SetTitleOffset(1.8);
    cout << "The bin width is " << hmass->GetBinWidth(1) * 1000 << " MeV/c^2" << endl;
    hmass->Rebin(3);
    hmass->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hmass->GetBinWidth(1) * 1000));
    hmass->Write("rawInvMass_PhiPhi");
    // hmass->GetXaxis()->SetRangeUser(2.5, 2.9);
    hmass->GetXaxis()->SetRangeUser(2.654, 2.774);
    hmass->Draw();

    // // TF1 *BWpol = new TF1("BWpol", BW_pol2, 2.68, 2.73, 6);
    // TF1 *BWpol = new TF1("BWpol", BW_pol3, 2.68, 2.73, 7); // 2.685 - 2.727 (fitting is fine)
    // if (BWpol->GetNpar() == 4)
    //     BWpol->SetParNames("Yield", "Mass", "Width", "p0", "p1", "p2", "p3");
    // else
    //     BWpol->SetParNames("Yield", "Mass", "Width", "p0", "p1", "p2");

    // BWpol->SetParameter(0, 10);
    // BWpol->SetParameter(1, 2.71);
    // BWpol->SetParameter(2, 0.016);
    // BWpol->SetParameter(3, 1);
    // BWpol->SetParameter(4, 0);
    // BWpol->SetParameter(5, 0);
    // BWpol->SetParameter(6, 0);
    // BWpol->SetParLimits(0, 0, 1e2);
    // BWpol->SetParLimits(1, 2.7, 2.72);
    // BWpol->SetParLimits(2, 0.002, 0.03);
    // //Bkg parameters
    // BWpol->SetParameter(3, 1.87330e+04);
    // BWpol->SetParameter(4, 7.50253e+05);
    // BWpol->SetParameter(5, -5.52942e+05);
    // if (BWpol->GetNpar() == 7)
    //     BWpol->SetParameter(6, 1.01426e+05);
    // hmass->Fit(BWpol, "RI");
    // TF1 *BWfunc = new TF1("BWfunc", BW, BWpol->GetXmin(), BWpol->GetXmax(), 3);
    // BWfunc->SetParameter(0, BWpol->GetParameter(0));
    // BWfunc->SetParameter(1, BWpol->GetParameter(1));
    // BWfunc->SetParameter(2, BWpol->GetParameter(2));
    // BWfunc->SetLineColor(kRed);
    // BWfunc->SetLineStyle(2);
    // BWfunc->Draw("same");
    // // TF1 *bkgfunc = new TF1("bkgfunc", pol2, BWpol->GetXmin(), BWpol->GetXmax(), 3);
    // TF1 *bkgfunc = new TF1("bkgfunc", pol3, BWpol->GetXmin(), BWpol->GetXmax(), 4);
    // bkgfunc->SetParameter(0, BWpol->GetParameter(3));
    // bkgfunc->SetParameter(1, BWpol->GetParameter(4));
    // bkgfunc->SetParameter(2, BWpol->GetParameter(5));
    // if (bkgfunc->GetNpar() == 4)
    //     bkgfunc->SetParameter(3, BWpol->GetParameter(6));
    // bkgfunc->SetLineColor(kGreen + 3);
    // bkgfunc->SetLineStyle(2);
    // bkgfunc->Draw("same");

    if (isSavePlots)
        c1->SaveAs(outputPath + "rawInvMass_PhiPhi.png");

    TLegend *leg = new TLegend(0.25, 0.18, 0.5, 0.3);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry(hmass, "Same event #Phi#Phi", "lpe");
    leg->Draw();

    TCanvas *cDeltaRKaon = new TCanvas("cDeltaRKaon", "#DeltaR Kaon vs p_{T}", 720, 720);
    SetCanvasStyle(cDeltaRKaon, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hDeltaRKaon);
    hDeltaRKaon->GetXaxis()->SetTitle("#DeltaR Kaon (cm)");
    hDeltaRKaon->SetMarkerStyle(20);
    hDeltaRKaon->GetYaxis()->SetMaxDigits(3);
    hDeltaRKaon->SetMarkerSize(1.0);
    hDeltaRKaon->Draw("HIST");
    if (isSavePlots)
        cDeltaRKaon->SaveAs(outputPath + "deltaR_kaon.png");

    TCanvas *cDeltaRPhi = new TCanvas("cDeltaRPhi", "#DeltaR #phi vs p_{T}", 720, 720);
    SetCanvasStyle(cDeltaRPhi, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hDeltaRPhiPhi);
    hDeltaRPhiPhi->SetTitle(0);
    hDeltaRPhiPhi->GetXaxis()->SetTitle("#DeltaR #phi (cm)");
    hDeltaRPhiPhi->SetMarkerStyle(20);
    hDeltaRPhiPhi->GetYaxis()->SetMaxDigits(3);
    hDeltaRPhiPhi->SetMarkerSize(1.0);
    hDeltaRPhiPhi->Draw("HIST");
    if (isSavePlots)
        cDeltaRPhi->SaveAs(outputPath + "deltaR_phi.png");

    TCanvas *cDeltaMassDifference = new TCanvas("cDeltaMassDifference", "#DeltaM vs p_{T}", 720, 720);
    SetCanvasStyle(cDeltaMassDifference, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hDeltaMassDifference);
    hDeltaMassDifference->GetXaxis()->SetTitle("#DeltaM (GeV/#  it{c}^{2})");
    hDeltaMassDifference->SetMarkerStyle(20);
    hDeltaMassDifference->GetYaxis()->SetMaxDigits(3);
    hDeltaMassDifference->SetMarkerSize(1.0);
    hDeltaMassDifference->Draw("HIST");
    if (isSavePlots)
        cDeltaMassDifference->SaveAs(outputPath + "deltaM.png");

    TCanvas *cPtCorrelation = new TCanvas("cPtCorrelation", "Pt Correlation", 720, 720);
    SetCanvasStyle(cPtCorrelation, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hPtCorrelation);
    hPtCorrelation->GetXaxis()->SetTitle("Pt Correlation (GeV/#it{c})");
    hPtCorrelation->SetMarkerStyle(20);
    hPtCorrelation->GetYaxis()->SetMaxDigits(3);
    hPtCorrelation->SetMarkerSize(1.0);
    hPtCorrelation->GetXaxis()->SetRangeUser(0.0, 10.0);
    int binMax = hPtCorrelation->GetMaximumBin();
    double maxXvalue = hPtCorrelation->GetXaxis()->GetBinCenter(binMax);
    hPtCorrelation->Draw("HIST");
    lat.DrawLatex(0.7, 0.82, Form("Peak %.1f", maxXvalue));
    lat.DrawLatex(0.7, 0.75, Form("Mean %.1f", hPtCorrelation->GetMean()));
    if (isSavePlots)
        cPtCorrelation->SaveAs(outputPath + "pt_correlation.png");

    TH1D *hmassmixed;
    if (MixedEventExist)
    {
        TCanvas *cSigBkg = new TCanvas("cSigBkg", "#Phi#Phi invariant mass with mixed event background", 720, 720);
        SetCanvasStyle(cSigBkg, 0.17, 0.03, 0.05, 0.15);
        SetHistoQA(hmassClone);
        hmassClone->GetXaxis()->SetTitle("#it{M}_{#Phi#Phi} (GeV/#it{c}^{2})");
        hmassClone->GetYaxis()->SetTitleOffset(1.8);
        hmassClone->SetMarkerStyle(20);
        hmassClone->GetYaxis()->SetMaxDigits(3);
        hmassClone->SetMarkerSize(1.0);
        // hmassClone->SetMaximum(1.2 * hmassClone->GetMaximum());
        hmassClone->Draw("pe");
        hmassmixed = hmixed->Projection(0, "E");
        SetHistoQA(hmassmixed);
        int normlow = hmassmixed->GetXaxis()->FindBin(2.8 + 0.01);
        int normhigh = hmassmixed->GetXaxis()->FindBin(2.9 - 0.01);
        auto signal_integral = hmassClone->Integral(normlow, normhigh);
        auto bkg_integral = hmassmixed->Integral(normlow, normhigh);
        auto normfactor = signal_integral / bkg_integral;
        hmassmixed->Scale(normfactor);
        // cout<<"Signal integral: "<<signal_integral<<endl;
        // cout<<"Bkg integral: "<<bkg_integral<<endl;
        // cout<<"Normalization factor for mixed event background: "<<normfactor<<endl;
        hmassClone->Rebin(3);
        hmassmixed->Rebin(3);
        hmassClone->GetXaxis()->SetRangeUser(2.5, 2.9);
        hmassmixed->GetXaxis()->SetRangeUser(2.5, 2.9);
        hmassClone->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hmassClone->GetBinWidth(1) * 1000));

        hmassmixed->SetLineColor(kRed);
        hmassmixed->SetMarkerColor(kRed);
        hmassmixed->Draw("pe same");

        TLegend *leg2 = new TLegend(0.25, 0.18, 0.5, 0.3);
        leg2->SetFillStyle(0);
        leg2->SetBorderSize(0);
        leg2->SetTextFont(42);
        leg2->SetTextSize(0.035);
        leg2->AddEntry(hmassClone, "Same event #Phi#Phi", "lpe");
        leg2->AddEntry(hmassmixed, "Mixed event #Phi#Phi", "lpe");
        leg2->Draw();

        if (isSavePlots)
            cSigBkg->SaveAs(outputPath + "rawInvMass_PhiPhi_withME.png");

        TCanvas *c2 = new TCanvas("c2", "#Phi#Phi invariant mass after bkg subtraction", 720, 720);
        SetCanvasStyle(c2, 0.15, 0.03, 0.05, 0.15);
        TH1D *hmassSubtracted = (TH1D *)hmassClone->Clone();
        hmassSubtracted->Add(hmassmixed, -1);
        SetHistoQA(hmassSubtracted);
        hmassSubtracted->GetXaxis()->SetTitle("#it{M}_{#Phi#Phi} (GeV/#it{c}^{2})");
        hmassSubtracted->SetMarkerStyle(20);
        hmassSubtracted->GetYaxis()->SetMaxDigits(3);
        hmassSubtracted->SetMarkerSize(1.0);
        hmassSubtracted->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hmassSubtracted->GetBinWidth(1) * 1000));
        hmassSubtracted->Write("bkgsubInvMass_PhiPhi");
        // hmassSubtracted->GetXaxis()->SetRangeUser(2.654, 2.774);
        hmassSubtracted->GetXaxis()->SetRangeUser(2.5, 2.9);
        hmassSubtracted->Draw();
        if (isSavePlots)
            c2->SaveAs(outputPath + "bkgsubInvMass_PhiPhi.png");
    }

    //===============Phi meson QA plots==================
    TH2F *hPhiMassvsPt = (TH2F *)file->Get("doublephimeson" + subWagon + "/hPhiMass2");
    if (hPhiMassvsPt == nullptr)
    {
        std::cerr << "Error: Could not find the phi invariant mass histogram in file\n";
        return;
    }

    if (isPhiQA)
    {

        // float pTbinsPhi[] = {0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 4.0, 5.0, 7.0, 10.0};
        float pTbinsPhi[] = {0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0};

        // float pTbinsPhi[] = {0.0, 100.0};
        int nPtBinsPhi = sizeof(pTbinsPhi) / sizeof(pTbinsPhi[0]) - 1;
        TH1F *hPhiMassFit = new TH1F("hPhiMassFit", "hPhiMassFit", nPtBinsPhi, pTbinsPhi);
        TH1F *hPhiMassResolutionFit = new TH1F("hPhiMassResolutionFit", "hPhiMassResolutionFit", nPtBinsPhi, pTbinsPhi);
        TH3F *hPhiPhiMassCorrelation = (TH3F *)file->Get("doublephimeson" + subWagon + "/hPhiMass3");
        TH2F *hDeltaRKaonPvspt = (TH2F *)file->Get("doublephimeson" + subWagon + "/hDeltaRkaonplusvspt");
        TH2F *hDeltaRKaonMvspt = (TH2F *)file->Get("doublephimeson" + subWagon + "/hDeltaRkaonminusvspt");
        if (hPhiPhiMassCorrelation == nullptr || hDeltaRKaonPvspt == nullptr || hDeltaRKaonMvspt == nullptr)
        {
            std::cerr << "Error: Could not find the correlation plots in file\n";
            return;
        }
        TH1F *hPhiYieldFit = new TH1F("hPhiYieldFit", "hPhiYieldFit", nPtBinsPhi, pTbinsPhi);
        TH1F *hPhiYieldFitLoosePID = (TH1F *)file2->Get("FittedPhiYield");
        if (hPhiYieldFitLoosePID == nullptr)
        {
            std::cerr << "Error: Could not find the FittedPhiYield histogram in file2 (loose PID cuts)\n";
            return;
        }
        TH1F *hEfficiency = new TH1F("hEfficiency", "hEfficiency", nPtBinsPhi, pTbinsPhi);
        TH1F *hChi2byNDF = new TH1F("hChi2byNDF", "hChi2byNDF", nPtBinsPhi, pTbinsPhi);
        TH1F *hPurity = new TH1F("hPurity", "hPurity", nPtBinsPhi, pTbinsPhi);
        TCanvas *cPhiSignal = new TCanvas("", "", 1440, 720);
        SetCanvasStyle(cPhiSignal, 0.13, 0.06, 0.05, 0.13);
        cPhiSignal->Divide(4, 4);
        TCanvas *cPhiMassCorrelation = new TCanvas("", "cPhiMassCorrelation", 1440, 720);
        SetCanvasStyle(cPhiMassCorrelation, 0.13, 0.06, 0.05, 0.13);
        cPhiMassCorrelation->Divide(4, 4);
        TCanvas *cDeltaRKaonP = new TCanvas("", "cDeltaRKaonP", 1440, 720);
        SetCanvasStyle(cDeltaRKaonP, 0.13, 0.06, 0.05, 0.13);
        cDeltaRKaonP->Divide(3, 4);
        TCanvas *cDeltaRPhi = new TCanvas("", "cDeltaRPhi", 1440, 720);
        SetCanvasStyle(cDeltaRPhi, 0.13, 0.06, 0.05, 0.13);
        cDeltaRPhi->Divide(3, 4);

        for (int ibinsPhi = 0; ibinsPhi < nPtBinsPhi; ibinsPhi++)
        {
            cPhiSignal->cd(ibinsPhi + 1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);

            int pTbinLow = hPhiMassvsPt->GetYaxis()->FindBin(pTbinsPhi[ibinsPhi] + 0.001);
            int pTbinHigh = hPhiMassvsPt->GetYaxis()->FindBin(pTbinsPhi[ibinsPhi + 1] - 0.001);
            hunlike->GetAxis(2)->SetRange(hunlike->GetAxis(2)->FindBin(pTbinsPhi[ibinsPhi] + 0.001), hunlike->GetAxis(2)->FindBin(pTbinsPhi[ibinsPhi + 1] - 0.001)); // Reset pT axis range
            TH1D *hDeltaRPhivsPt = hunlike->Projection(3, "E");

            TH1D *hphiMassProj = hPhiMassvsPt->ProjectionX(Form("hphiMassProj_%.1f_%.1f", pTbinsPhi[ibinsPhi], pTbinsPhi[ibinsPhi + 1]), pTbinLow, pTbinHigh);
            // hPhiPhiMassCorrelation->GetZaxis()->SetRange(pTbinLow, pTbinHigh);
            // TH2F *hPhiPhiMassCorrProj = (TH2F *)hPhiPhiMassCorrelation->Project3D("xy");
            // SetHistoQA2D(hPhiPhiMassCorrProj);
            TH1D *hDeltaRKaonPvsptProj = hDeltaRKaonPvspt->ProjectionX(Form("hDeltaRKaonPvsptProj_%.1f_%.1f", pTbinsPhi[ibinsPhi], pTbinsPhi[ibinsPhi + 1]), pTbinLow, pTbinHigh);
            TH1D *hDeltaRKaonMvsptProj = hDeltaRKaonMvspt->ProjectionX(Form("hDeltaRKaonMvsptProj_%.1f_%.1f", pTbinsPhi[ibinsPhi], pTbinsPhi[ibinsPhi + 1]), pTbinLow, pTbinHigh);
            SetHistoQA(hDeltaRKaonPvsptProj);
            SetHistoQA(hDeltaRKaonMvsptProj);

            SetHistoQA(hphiMassProj);
            hphiMassProj->GetYaxis()->SetTitleOffset(1.3);
            hphiMassProj->GetYaxis()->SetMaxDigits(3);
            hphiMassProj->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
            hphiMassProj->GetXaxis()->SetNdivisions(505);
            // hphiMassProj->Rebin(4);
            hphiMassProj->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", (hphiMassProj->GetXaxis()->GetBinWidth(1) * 1000)));
            hphiMassProj->SetMaximum(1.2 * hphiMassProj->GetMaximum());
            hphiMassProj->SetMinimum(0);
            hphiMassProj->Draw("pe");
            lat.DrawLatex(0.65, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", pTbinsPhi[ibinsPhi], pTbinsPhi[ibinsPhi + 1]));

            float fitRangeLow = 1.001;
            float fitRangeHigh = 1.039;
            // if (ibinsPhi == 0 || ibinsPhi == nPtBinsPhi - 1)
            if (ibinsPhi == 0)
            {
                fitRangeLow = 1.005;
                fitRangeHigh = 1.033;
            }
            else
            {
                fitRangeLow = 1.001;
                fitRangeHigh = 1.039;
            }

            //******Fitting for Phi*********************
            TF1 *fitFcn = new TF1("fitfunc", voigtpol2, fitRangeLow, fitRangeHigh, 7);       // sig+bkg fit function
            TF1 *fitFcnBkg = new TF1("fitfunc1", polynomial2, fitRangeLow, fitRangeHigh, 3); // only residualbkg
            TF1 *fitFcnSig = new TF1("fitFcnSig", voigt, fitRangeLow, fitRangeHigh, 4);      // only signal

            // for voigtian distribution
            fitFcn->SetParameter(0, 5000); // yield
            // if (ibinsPhi == 3 || ibinsPhi == nPtBinsPhi - 1)
            if (ibinsPhi == 3)
                fitFcn->SetParLimits(0, 0, 1e4); // yield
            else
                fitFcn->SetParLimits(0, 0, 1e6);    // yield
            fitFcn->SetParameter(1, 1.019);         // mass peak
            fitFcn->SetParLimits(1, 1.006, 1.0316); // mass peak
            fitFcn->SetParameter(2, 0.0012);        //  Gaussian width (Detector resolution)
            fitFcn->SetParLimits(2, 0.001, 0.009);  // Gaussian width.
            // fitFcn->SetParameter(3, 0.0042);   //lorentzian width (Resonance width)
            fitFcn->FixParameter(3, 0.0042); // lorentzian width

            TFitResultPtr r = hphiMassProj->Fit("fitfunc", "REBMS");
            fitFcnBkg->SetParameters(fitFcn->GetParameter(4), fitFcn->GetParameter(5), fitFcn->GetParameter(6));
            fitFcnSig->SetParameters(fitFcn->GetParameter(0), fitFcn->GetParameter(1), fitFcn->GetParameter(2), fitFcn->GetParameter(3));
            fitFcnBkg->SetLineColor(kBlue);
            fitFcnSig->SetLineColor(kGreen + 2);
            fitFcnBkg->SetLineStyle(2);
            fitFcnSig->SetLineStyle(2);
            // fitFcnSig->SetNpx(10000);
            fitFcnBkg->Draw("same");
            fitFcnSig->Draw("same");

            // TCanvas *cIndividualPlots = new TCanvas("", "", 720, 720);
            // SetCanvasStyle(cIndividualPlots, 0.14, 0.03, 0.05, 0.14);
            // SetHistoQA(hphiMassProj);
            // hphiMassProj->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
            // hphiMassProj->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", (hphiMassProj->GetXaxis()->GetBinWidth(1) * 1000)));
            // hphiMassProj->SetMaximum(1.2 * hphiMassProj->GetMaximum());
            // hphiMassProj->SetMinimum(0);
            // hphiMassProj->GetXaxis()->SetRangeUser(0.99, 1.05);
            // hphiMassProj->Draw("pe");
            // fitFcn->Draw("same");
            // fitFcnBkg->Draw("same");
            // fitFcnSig->Draw("same");

            // TLegend *lp = DrawLegend(0.18, 0.75, 0.6, 0.92);
            // lp->SetTextSize(0.03);
            // lp->SetTextFont(42);
            // lp->SetFillStyle(0);
            // lp->AddEntry(hphiMassProj, "K^{+}K^{-} invariant mass", "pe");
            // lp->AddEntry(fitFcn, "Combined fit", "l");
            // lp->AddEntry(fitFcnSig, "Signal (Voigt)", "l");
            // lp->AddEntry(fitFcnBkg, "Residual Background (Pol2)", "l");
            // lp->Draw("same");
            // lat.SetTextSize(0.08);
            // lat.DrawLatex(0.20, 0.71, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", pTbinsPhi[ibinsPhi], pTbinsPhi[ibinsPhi + 1]));
            // if (isSavePlots)
            //     cIndividualPlots->SaveAs(outputPhiFits + Form("phi_mass_fit_%.1f_%.1f.png", pTbinsPhi[ibinsPhi], pTbinsPhi[ibinsPhi + 1]));

            double mean = fitFcn->GetParameter(1);
            // double mean = 1.019; // PDG value
            double widthPDG = 0.0042;
            float nSigmaForPurity = 2.0;
            double areabkg = fitFcnBkg->Integral(mean - nSigmaForPurity * widthPDG, mean + nSigmaForPurity * widthPDG) / (hphiMassProj->GetXaxis()->GetBinWidth(1));
            // double areabkg = fitFcnBkg->Integral(1.001, 1.039) / (hphiMassProj->GetXaxis()->GetBinWidth(1));
            double areasigbkg = hphiMassProj->Integral(hphiMassProj->FindBin(mean - nSigmaForPurity * widthPDG), hphiMassProj->FindBin(mean + nSigmaForPurity * widthPDG));
            // double areasigbkg = hphiMassProj->Integral(hphiMassProj->FindBin(1.001), hphiMassProj->FindBin(1.039));
            double areasig = areasigbkg - areabkg;
            double purity = areasig * 100 / areasigbkg;
            hPurity->SetBinContent(ibinsPhi + 1, purity);
            hPurity->SetBinError(ibinsPhi + 1, 0);

            hPhiMassFit->SetBinContent(ibinsPhi + 1, fitFcn->GetParameter(1));
            hPhiMassFit->SetBinError(ibinsPhi + 1, fitFcn->GetParError(1));

            TMatrixDSym cov = r->GetCovarianceMatrix();
            TMatrixDSym cov1;
            TMatrixDSym cov2;
            cov.GetSub(0, 3, 0, 3, cov1);
            cov.GetSub(3, 6, 3, 6, cov2);
            Double_t *b = cov1.GetMatrixArray();
            Double_t *a = cov2.GetMatrixArray();
            Double_t *para = fitFcn->GetParameters();
            float ptBinWidth = pTbinsPhi[ibinsPhi + 1] - pTbinsPhi[ibinsPhi];
            double yieldIntegral = fitFcnSig->Integral(mean - 3 * widthPDG, mean + 3 * widthPDG) / ptBinWidth;
            double yieldIntegralError = fitFcnSig->IntegralError(mean - 3 * widthPDG, mean + 3 * widthPDG, &para[0], b) / ptBinWidth;

            cout << "pT bin " << ibinsPhi << ", yield " << yieldIntegral << ", error " << yieldIntegralError << endl;

            hPhiYieldFit->SetBinContent(ibinsPhi + 1, yieldIntegral);
            hPhiYieldFit->SetBinError(ibinsPhi + 1, yieldIntegralError);
            if (hPhiYieldFitLoosePID->GetNbinsX() == nPtBinsPhi)
            {
                double yieldIntegral2 = hPhiYieldFitLoosePID->GetBinContent(ibinsPhi + 1);
                double yieldIntegralError2 = hPhiYieldFitLoosePID->GetBinError(ibinsPhi + 1);
                double ratio = yieldIntegral / yieldIntegral2;
                double ratioError = ratio * sqrt(pow(yieldIntegralError / yieldIntegral2, 2) + pow(yieldIntegralError2 * yieldIntegral / (yieldIntegral2 * yieldIntegral2), 2));
                hEfficiency->SetBinContent(ibinsPhi + 1, ratio * 100);
                hEfficiency->SetBinError(ibinsPhi + 1, ratioError * 100);
            }

            hPhiMassResolutionFit->SetBinContent(ibinsPhi + 1, fitFcn->GetParameter(2));
            hPhiMassResolutionFit->SetBinError(ibinsPhi + 1, fitFcn->GetParError(2));

            hChi2byNDF->SetBinContent(ibinsPhi + 1, fitFcn->GetChisquare() / fitFcn->GetNDF());

            hphiMassProj->Write();

            // cPhiMassCorrelation->cd(ibinsPhi + 1);
            // gPad->SetLeftMargin(0.15);
            // gPad->SetBottomMargin(0.15);
            // gPad->SetRightMargin(0.05);
            // gPad->SetTopMargin(0.05);
            // hPhiPhiMassCorrProj->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
            // hPhiPhiMassCorrProj->GetYaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
            // hPhiPhiMassCorrProj->Draw("COLZ");
            // TLine *lineHorizontalPDG = new TLine(1.019461, 1.0, 1.019461, 1.04);
            // lineHorizontalPDG->SetLineColor(kRed);
            // lineHorizontalPDG->SetLineStyle(2);
            // TLine *lineVerticalPDG = new TLine(1.0, 1.019461, 1.04, 1.019461);
            // lineVerticalPDG->SetLineColor(kRed);
            // lineVerticalPDG->SetLineStyle(2);
            // lineHorizontalPDG->Draw("same");
            // lineVerticalPDG->Draw("same");
            // lat.DrawLatex(0.6, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", pTbinsPhi[ibinsPhi], pTbinsPhi[ibinsPhi + 1]));

            cDeltaRKaonP->cd(ibinsPhi + 1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            hDeltaRKaonPvsptProj->GetXaxis()->SetRangeUser(0, 2);
            hDeltaRKaonPvsptProj->GetXaxis()->SetTitle("#DeltaR K^{+} (cm)");
            hDeltaRKaonPvsptProj->Draw("HIST");
            lat.DrawLatex(0.6, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", pTbinsPhi[ibinsPhi], pTbinsPhi[ibinsPhi + 1]));

            cDeltaRPhi->cd(ibinsPhi + 1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            hDeltaRPhivsPt->SetTitle(0);
            hDeltaRPhivsPt->GetXaxis()->SetTitle("#DeltaR #phi (cm)");
            hDeltaRPhivsPt->Draw("HIST");
            lat.SetTextFont(22);
            lat.DrawLatex(0.6, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", pTbinsPhi[ibinsPhi], pTbinsPhi[ibinsPhi + 1]));
        }
        cPhiSignal->Write("AllPhiMassFits");
        if (isSavePlots)
        {
            cPhiSignal->SaveAs(outputPath + "/phi_mass_allPtBins.png");
            cDeltaRKaonP->SaveAs(outputPath + "/deltaR_kaon_plus_vs_pt.png");
            // cPhiMassCorrelation->SaveAs(outputPath + "/phi_phi_mass_correlation.png");
            cDeltaRPhi->SaveAs(outputPath + "/deltaR_phi_vs_pt.png");
        }

        TCanvas *cPhiPhiMassCorrelation = new TCanvas("cPhiPhiMassCorrelation", "cPhiPhiMassCorrelation", 720, 720);
        SetCanvasStyle(cPhiPhiMassCorrelation, 0.19, 0.15, 0.05, 0.13);
        double deltaMLow = hPhiPhiMassCorrelation->GetZaxis()->FindBin(0.0);
        double deltaMHigh = hPhiPhiMassCorrelation->GetZaxis()->FindBin(0.015); // 0.01 is good cut (60% statistics are lost)
        cout << "Entries before the cut " << hPhiPhiMassCorrelation->Integral() << endl;

        hPhiPhiMassCorrelation->GetZaxis()->SetRange(deltaMLow, deltaMHigh);
        cout << "Entries after the cut " << hPhiPhiMassCorrelation->Integral() << endl;

        TH2F *hPhiPhiMassCorrProj = (TH2F *)hPhiPhiMassCorrelation->Project3D("xy");
        SetHistoQA2D(hPhiPhiMassCorrProj);
        hPhiPhiMassCorrProj->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
        hPhiPhiMassCorrProj->GetYaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
        hPhiPhiMassCorrProj->GetYaxis()->SetTitleOffset(1.9);
        hPhiPhiMassCorrProj->GetXaxis()->SetNdivisions(505);
        hPhiPhiMassCorrProj->Draw("COLZ");
        TLine *lineHorizontalPDG = new TLine(1.019461, 1.0, 1.019461, 1.04);
        lineHorizontalPDG->SetLineColor(kRed);
        lineHorizontalPDG->SetLineStyle(2);
        TLine *lineVerticalPDG = new TLine(1.0, 1.019461, 1.04, 1.019461);
        lineVerticalPDG->SetLineColor(kRed);
        lineVerticalPDG->SetLineStyle(2);
        lineHorizontalPDG->Draw("same");
        lineVerticalPDG->Draw("same");
        if (isSavePlots)
            cPhiPhiMassCorrelation->SaveAs(outputPath + "/phi_phi_mass_correlation.png");

        TCanvas *cPurity = new TCanvas("cPurity", "Phi Purity", 720, 720);
        SetCanvasStyle(cPurity, 0.17, 0.03, 0.05, 0.15);
        SetHistoQA(hPurity);
        hPurity->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hPurity->GetYaxis()->SetTitle("Purity (%)");
        hPurity->SetMarkerStyle(20);
        hPurity->SetMarkerSize(1);
        hPurity->SetLineWidth(2);
        hPurity->GetYaxis()->SetRangeUser(0, 100);
        hPurity->Draw("HIST");
        hPurity->Write("PhiPurity");
        if (isSavePlots)
            cPurity->SaveAs(outputPath + "/phi_purity.png");

        TCanvas *cEfficiency = new TCanvas("cEfficiency", "Phi Efficiency", 720, 720);
        SetCanvasStyle(cEfficiency, 0.17, 0.03, 0.05, 0.15);
        SetHistoQA(hEfficiency);
        hEfficiency->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hEfficiency->GetYaxis()->SetTitle("Efficiency (PID selection / PID loose) (%)");
        hEfficiency->SetMarkerStyle(20);
        hEfficiency->SetMarkerSize(1);
        hEfficiency->SetLineWidth(2);
        hEfficiency->GetYaxis()->SetRangeUser(0, 100);
        hEfficiency->Draw("HIST");
        hEfficiency->Write("PhiEfficiency");
        if (isSavePlots)
            cEfficiency->SaveAs(outputPath + "/phi_efficiency.png");

        // TCanvas *cFitMass = new TCanvas("cFitMass", "Phi Mass Fit Parameters", 720, 720);
        // SetCanvasStyle(cFitMass, 0.17, 0.03, 0.05, 0.15);
        // SetHistoQA(hPhiMassFit);
        // hPhiMassFit->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        // hPhiMassFit->GetYaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
        // hPhiMassFit->SetMarkerStyle(20);
        // hPhiMassFit->SetMarkerSize(1);
        // hPhiMassFit->SetLineWidth(2);
        // hPhiMassFit->Write("FittedPhiMass");
        // hPhiMassFit->Draw("pe");
        // TLine *linePDG = new TLine(pTbinsPhi[0], 1.019461, pTbinsPhi[nPtBinsPhi], 1.019461);
        // linePDG->SetLineColor(kRed);
        // linePDG->SetLineStyle(2);
        // linePDG->SetLineWidth(2);
        // linePDG->Draw("same");
        // TLegend *legMass = new TLegend(0.25, 0.18, 0.5, 0.3);
        // legMass->SetFillStyle(0);
        // legMass->SetBorderSize(0);
        // legMass->SetTextFont(42);
        // legMass->SetTextSize(0.035);
        // legMass->AddEntry(hPhiMassFit, "Fitted #phi mass", "lpe");
        // legMass->AddEntry(linePDG, "#phi PDG mass", "l");
        // legMass->Draw();
        // if (isSavePlots)
        //     cFitMass->SaveAs(outputPath + "/phi_mass_fit.png");

        // TCanvas *cFitMassResolution = new TCanvas("cFitMassResolution", "Phi Mass Resolution Fit Parameters", 720, 720);
        // SetCanvasStyle(cFitMassResolution, 0.17, 0.03, 0.05, 0.15);
        // SetHistoQA(hPhiMassResolutionFit);
        // hPhiMassResolutionFit->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        // hPhiMassResolutionFit->GetYaxis()->SetTitle("#sigma_{#phi} (GeV/#it{c}^{2})");
        // hPhiMassResolutionFit->SetMarkerStyle(20);
        // hPhiMassResolutionFit->SetMarkerSize(1);
        // hPhiMassResolutionFit->SetLineWidth(2);
        // hPhiMassResolutionFit->Draw("pe");
        // hPhiMassResolutionFit->Write("FittedPhiMassResolution");
        // if (isSavePlots)
        //     cFitMassResolution->SaveAs(outputPath + "/phi_mass_resolution_fit.png");

        TCanvas *cFitYield = new TCanvas("cFitYield", "Phi Yield Fit Parameters", 720, 720);
        SetCanvasStyle(cFitYield, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hPhiYieldFit);
        // gPad->SetLogy();
        hPhiYieldFit->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hPhiYieldFit->GetYaxis()->SetTitle("Raw yield");
        hPhiYieldFit->SetMarkerStyle(20);
        hPhiYieldFit->SetMarkerSize(1);
        hPhiYieldFit->SetLineWidth(2);
        hPhiYieldFit->Draw("HIST");
        hPhiYieldFit->Write("FittedPhiYield");
        if (isSavePlots)
            cFitYield->SaveAs(outputPath + "/phi_yield_fit.png");

        // TCanvas *cFitChi2Ndf = new TCanvas("cFitChi2Ndf", "Phi Chi2/NDF Fit Parameters", 720, 720);
        // SetCanvasStyle(cFitChi2Ndf, 0.17, 0.03, 0.05, 0.15);
        // SetHistoQA(hChi2byNDF);
        // hChi2byNDF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        // hChi2byNDF->GetYaxis()->SetTitle("#chi^{2}/NDF");
        // hChi2byNDF->SetMarkerStyle(20);
        // hChi2byNDF->SetMarkerSize(1);
        // hChi2byNDF->SetLineWidth(2);
        // hChi2byNDF->Draw("HIST");
        // hChi2byNDF->Write("FittedPhiChi2NDF");
        // if (isSavePlots)
        //     cFitChi2Ndf->SaveAs(outputPath + "/phi_chi2NDF_fit.png");
    } // end of isPhiQA

    //==========================================================================//
    //=============================== QA Plots =================================//
    //==========================================================================//

    if (isQA)
    {
        // TCanvas *c2DPhiMassvsPt = new TCanvas("c2DPhiMassvsPt", "Phi Mass vs Pt", 720, 720);
        // SetCanvasStyle(c2DPhiMassvsPt, 0.12, 0.15, 0.05, 0.12);
        // SetHistoQA(hPhiMassvsPt);
        // hPhiMassvsPt->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
        // hPhiMassvsPt->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        // hPhiMassvsPt->GetYaxis()->SetTitleOffset(1.3);
        // hPhiMassvsPt->GetYaxis()->SetMaxDigits(3);
        // hPhiMassvsPt->GetXaxis()->SetNdivisions(505);
        // hPhiMassvsPt->GetZaxis()->SetMaxDigits(3);
        // hPhiMassvsPt->Draw("colz");
        // if (isSavePlots)
        //     c2DPhiMassvsPt->SaveAs(outputPathQA + "2DPhiMassvsPt.png");

        // TCanvas *cDeltaR = new TCanvas("cDeltaR", "#DeltaR of #Phi#Phi pairs", 720, 720);
        // SetCanvasStyle(cDeltaR, 0.15, 0.03, 0.05, 0.15);
        // SetHistoQA(hDeltaRPhiPhi);
        // hDeltaRPhiPhi->GetXaxis()->SetTitle("#DeltaR = #sqrt{(#Delta#eta)^{2}_{#Phi} + (#Delta#varphi)^{2}_{#Phi}}");
        // hDeltaRPhiPhi->GetYaxis()->SetTitle("Counts");
        // hDeltaRPhiPhi->Draw("HIST");
        // if (isSavePlots)
        //     cDeltaR->SaveAs(outputPathQA + "DeltaR_PhiPhi.png");

        // // TCanvas *cNumPhi = new TCanvas("cNumPhi", "Number of #Phi mesons in an event", 720, 720);
        // // SetCanvasStyle(cNumPhi, 0.15, 0.05, 0.05, 0.15);
        // // SetHistoQA(hNumPhi);
        // // gPad->SetLogy();
        // // hNumPhi->GetXaxis()->SetTitle("Number of #Phi mesons in an event");
        // // hNumPhi->GetYaxis()->SetTitle("Counts");
        // // hNumPhi->GetXaxis()->SetRangeUser(1, 10.0);
        // // hNumPhi->Draw("HIST");
        // // if (isSavePlots)
        // //     cNumPhi->SaveAs(outputPathQA + "Num_Phi_per_event.png");

        // TCanvas *cDeltaMassDifference = new TCanvas("cDeltaMassDifference", "Mass Difference of #Phi from PDG", 720, 720);
        // SetCanvasStyle(cDeltaMassDifference, 0.15, 0.03, 0.05, 0.15);
        // SetHistoQA(hDeltaMassDifference);
        // hDeltaMassDifference->GetXaxis()->SetTitle("#DeltaM = #sqrt{(#Phi_{1} - #Phi_{PDG})^{2} + (#Phi_{2} - #Phi_{PDG})^{2}} (GeV/#it{c}^{2})");
        // hDeltaMassDifference->GetYaxis()->SetTitle("Counts");
        // hDeltaMassDifference->Draw("HIST");
        // if (isSavePlots)
        //     cDeltaMassDifference->SaveAs(outputPathQA + "DeltaMassDifference_PhiPhi.png");

        // TH1F *hDeltaRKaonPlus = (TH1F *)file->Get("doublephimeson" + subWagon + "/hDeltaRkaonplus");
        // TH1F *hDeltaRKaonMinus = (TH1F *)file->Get("doublephimeson" + subWagon + "/hDeltaRkaonminus");

        // if (hDeltaRKaonPlus == nullptr || hDeltaRKaonMinus == nullptr)
        // {
        //     std::cerr << "Warning: Could not find histogram 'doublephimeson" + subWagon + "/hDeltaRkaonplus' or 'doublephimeson" + subWagon + "/hDeltaRkaonminus' in file\n";
        // }

        // TCanvas *cDeltaRKaon = new TCanvas("cDeltaRKaon", "#DeltaR of Kaons", 720, 720);
        // SetCanvasStyle(cDeltaRKaon, 0.15, 0.03, 0.05, 0.15);
        // SetHistoQA(hDeltaRKaonPlus);
        // SetHistoQA(hDeltaRKaonMinus);
        // hDeltaRKaonPlus->GetXaxis()->SetTitle("#DeltaR = #sqrt{(#Delta#eta)^{2}_{K} + (#Delta#varphi)^{2}_{K}}");
        // hDeltaRKaonPlus->GetYaxis()->SetTitle("Counts");
        // hDeltaRKaonPlus->SetMarkerStyle(20);
        // hDeltaRKaonPlus->SetMarkerSize(1.0);
        // hDeltaRKaonPlus->SetLineColor(kBlue);
        // hDeltaRKaonPlus->SetMarkerColor(kBlue);
        // hDeltaRKaonPlus->Rebin(2);
        // hDeltaRKaonMinus->Rebin(2);
        // hDeltaRKaonPlus->Draw();
        // hDeltaRKaonMinus->SetMarkerStyle(21);
        // hDeltaRKaonMinus->SetMarkerSize(1.0);
        // hDeltaRKaonMinus->SetLineColor(kRed);
        // hDeltaRKaonMinus->SetMarkerColor(kRed);
        // hDeltaRKaonMinus->Draw("same");
        // TLegend *legDeltaR = new TLegend(0.55, 0.7, 0.85, 0.85);
        // legDeltaR->SetFillStyle(0);
        // legDeltaR->SetBorderSize(0);
        // legDeltaR->SetTextFont(42);
        // legDeltaR->SetTextSize(0.035);
        // legDeltaR->AddEntry(hDeltaRKaonPlus, "K^{+}", "lpe");
        // legDeltaR->AddEntry(hDeltaRKaonMinus, "K^{-}", "lpe");
        // legDeltaR->Draw();
        // if (isSavePlots)
        //     cDeltaRKaon->SaveAs(outputPathQA + "DeltaR_Kaons.png");

        float pTbins[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0};
        int nPtBins = sizeof(pTbins) / sizeof(float) - 1;
        TLatex lat;
        lat.SetNDC();
        // lat.SetTextFont(42);
        lat.SetTextSize(0.08);

        TH2D *hnSigmaTPCKaonPlus = (TH2D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTPCKaonPlusBefore");
        if (hnSigmaTPCKaonPlus == nullptr)
        {
            std::cerr << "Error: Could not find histogram 'doublephimeson" + subWagon + "/hnsigmaTPCKaonPlus' in file\n";
            return;
        }
        TH2D *hnSigmaTPCKaonMinus = (TH2D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTPCKaonMinusBefore");
        if (hnSigmaTPCKaonMinus == nullptr)
        {
            std::cerr << "Error: Could not find histogram 'doublephimeson" + subWagon + "/hnsigmaTPCKaonMinus' in file\n";
            return;
        }
        TH3D *hnSigmaTPCTOFKaon = (TH3D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTPCTOFKaonBefore");
        if (hnSigmaTPCTOFKaon == nullptr)
        {
            std::cerr << "Error: Could not find histogram 'doublephimeson" + subWagon + "/hnsigmaTPCTOFKaon' in file\n";
            return;
        }
        TH2D *hnSigmaTOFKaonPlus = (TH2D *)file->Get("doublephimeson" + subWagon + "/hnsigmaTOFKaonPlusBefore");
        if (hnSigmaTOFKaonPlus == nullptr)
        {
            std::cerr << "Error: Could not find histogram 'doublephimeson" + subWagon + "/hnsigmaTOFKaonPlus' in file\n";
        }

        // TCanvas *cnSigmaTPCKaonPlus = new TCanvas("cnSigmaTPCKaonPlus", "nSigmaTPC Kaon Plus", 720, 720);
        // SetCanvasStyle(cnSigmaTPCKaonPlus, 0.12, 0.15, 0.05, 0.13);
        // SetHistoQA(hnSigmaTPCKaonPlus);
        // hnSigmaTPCKaonPlus->GetXaxis()->SetTitle("#it{n}_{#sigma} (TPC)");
        // hnSigmaTPCKaonPlus->GetYaxis()->SetTitle("#it{p} (GeV/#it{c})");
        // hnSigmaTPCKaonPlus->Draw("COLZ");
        // if (isSavePlots)
        //     cnSigmaTPCKaonPlus->SaveAs(outputPathQA + "2DSigmaTPC_KaonPlusvsPt.png");

        // TCanvas *cnSigmaTPCKaonMinus = new TCanvas("cnSigmaTPCKaonMinus", "nSigmaTPC Kaon Minus", 720, 720);
        // SetCanvasStyle(cnSigmaTPCKaonMinus, 0.12, 0.15, 0.05, 0.13);
        // SetHistoQA(hnSigmaTPCKaonMinus);
        // hnSigmaTPCKaonMinus->GetXaxis()->SetTitle("#it{n}_{#sigma} (TPC)");
        // hnSigmaTPCKaonMinus->GetYaxis()->SetTitle("#it{p} (GeV/#it{c})");
        // hnSigmaTPCKaonMinus->Draw("COLZ");
        // if (isSavePlots)
        //     cnSigmaTPCKaonMinus->SaveAs(outputPathQA + "2DSigmaTPC_KaonMinusvsPt.png");

        TCanvas *cnSigmaTPCKaonPlusvsPt = new TCanvas("cnSigmaTPCKaonPlusvsPt", "nSigmaTPC Kaon Plus vs Pt", 1440, 720);
        SetCanvasStyle(cnSigmaTPCKaonPlusvsPt, 0.11, 0.06, 0.05, 0.10);
        cnSigmaTPCKaonPlusvsPt->Divide(4, 5);
        for (int i = 0; i < nPtBins; i++)
        {
            cnSigmaTPCKaonPlusvsPt->cd(i + 1);
            gPad->SetBottomMargin(0.13);
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.03);
            gPad->SetTopMargin(0.01);
            gPad->SetGridx();
            float ptLow = pTbins[i];
            float ptHigh = pTbins[i + 1];
            int binLow = hnSigmaTPCKaonPlus->GetYaxis()->FindBin(ptLow + 0.001);
            int binHigh = hnSigmaTPCKaonPlus->GetYaxis()->FindBin(ptHigh - 0.001);
            TH1D *hKaonTPC1D = hnSigmaTPCKaonPlus->ProjectionX(Form("hKaonTPC1D_%.1f_%.1f", ptLow, ptHigh), binLow, binHigh);
            SetHistoQA(hKaonTPC1D);
            // hKaonTPC1D->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));
            hKaonTPC1D->GetXaxis()->SetTitle("#it{n}_{#sigma} (TPC)");
            hKaonTPC1D->GetYaxis()->SetTitle("Counts");
            hKaonTPC1D->SetMarkerStyle(20);
            hKaonTPC1D->SetMarkerSize(0.6);
            if (i > nPtBins - 5)
                hKaonTPC1D->Rebin(5);
            else
                hKaonTPC1D->Rebin(2);
            hKaonTPC1D->SetMaximum(1.5 * hKaonTPC1D->GetMaximum());
            hKaonTPC1D->Write(Form("hKaonTPC1DPos_%.1f_%.1f", ptLow, ptHigh));
            hKaonTPC1D->Draw("pe");
            lat.DrawLatex(0.5, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));
            TLine *line0 = new TLine(0, 0, 0, hKaonTPC1D->GetMaximum());
            line0->SetLineColor(kRed);
            line0->SetLineStyle(2);
            line0->Draw("same");
        }
        if (isSavePlots)
            cnSigmaTPCKaonPlusvsPt->SaveAs(outputPathQA + "nSigmaTPC_KaonPlus_vs_Pt.png");

        if (hnSigmaTOFKaonPlus != nullptr)
        {
            TCanvas *cnSigmaTOFKaonPlusvsPt = new TCanvas("cnSigmaTOFKaonPlusvsPt", "nSigmaTOF Kaon Plus vs Pt", 1080, 720);
            SetCanvasStyle(cnSigmaTOFKaonPlusvsPt, 0.11, 0.06, 0.05, 0.10);
            cnSigmaTOFKaonPlusvsPt->Divide(4, 4);
            for (int i = 0; i < nPtBins; i++)
            {
                gPad->SetBottomMargin(0.12);
                gPad->SetLeftMargin(0.12);
                gPad->SetRightMargin(0.03);
                gPad->SetTopMargin(0.01);
                gPad->SetGridx();
                cnSigmaTOFKaonPlusvsPt->cd(i + 1);
                float ptLow = pTbins[i];
                float ptHigh = pTbins[i + 1];
                int binLow = hnSigmaTOFKaonPlus->GetYaxis()->FindBin(ptLow + 0.001);
                int binHigh = hnSigmaTOFKaonPlus->GetYaxis()->FindBin(ptHigh - 0.001);
                TH1D *hKaonTOF1D = hnSigmaTOFKaonPlus->ProjectionX(Form("hKaonTOF1D_%.1f_%.1f", ptLow, ptHigh), binLow, binHigh);
                SetHistoQA(hKaonTOF1D);
                // hKaonTOF1D->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));
                hKaonTOF1D->GetXaxis()->SetTitle("#it{n}_{#sigma} (TOF)");
                hKaonTOF1D->GetYaxis()->SetTitle("Counts");
                hKaonTOF1D->SetMarkerStyle(20);
                hKaonTOF1D->SetMarkerSize(0.6);
                if (i > nPtBins - 5)
                    hKaonTOF1D->Rebin(4);
                hKaonTOF1D->SetMaximum(1.2 * hKaonTOF1D->GetMaximum());
                // hKaonTOF1D->GetXaxis()->SetRangeUser(-1.5,1.5);
                hKaonTOF1D->Draw("pe");
                lat.DrawLatex(0.5, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));
                TLine *line0 = new TLine(0, 0, 0, hKaonTOF1D->GetMaximum());
                line0->SetLineColor(kRed);
                line0->SetLineStyle(2);
                line0->Draw("same");
            }
            if (isSavePlots)
                cnSigmaTOFKaonPlusvsPt->SaveAs(outputPathQA + "nSigmaTOF_KaonPlus_vs_Pt.png");
        }

        // TCanvas *cnSigmaTPCKaonMinusvsPt = new TCanvas("cnSigmaTPCKaonMinusvsPt", "nSigmaTPC Kaon Minus vs Pt", 1440, 720);
        // SetCanvasStyle(cnSigmaTPCKaonMinusvsPt, 0.11, 0.06, 0.05, 0.10);
        // cnSigmaTPCKaonMinusvsPt->Divide(4, 4);
        // for (int i = 0; i < nPtBins; i++)
        // {
        //     gPad->SetBottomMargin(0.12);
        //     gPad->SetLeftMargin(0.12);
        //     gPad->SetRightMargin(0.03);
        //     gPad->SetTopMargin(0.01);
        //     gPad->SetGridx();
        //     cnSigmaTPCKaonMinusvsPt->cd(i + 1);
        //     float ptLow = pTbins[i];
        //     float ptHigh = pTbins[i + 1];
        //     int binLow = hnSigmaTPCKaonMinus->GetYaxis()->FindBin(ptLow + 0.001);
        //     int binHigh = hnSigmaTPCKaonMinus->GetYaxis()->FindBin(ptHigh - 0.001);
        //     TH1D *hKaonTPC1D = hnSigmaTPCKaonMinus->ProjectionX(Form("hKaonTPC1D_%.1f_%.1f", ptLow, ptHigh), binLow, binHigh);
        //     SetHistoQA(hKaonTPC1D);
        //     // hKaonTPC1D->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));
        //     hKaonTPC1D->GetXaxis()->SetTitle("#it{n}_{#sigma} (TPC)");
        //     hKaonTPC1D->GetYaxis()->SetTitle("Counts");
        //     hKaonTPC1D->SetMarkerStyle(20);
        //     hKaonTPC1D->SetMarkerSize(0.6);
        //     if (i > nPtBins - 5)
        //         hKaonTPC1D->Rebin(4);
        //     hKaonTPC1D->SetMaximum(1.2 * hKaonTPC1D->GetMaximum());
        //     // hKaonTPC1D->GetXaxis()->SetRangeUser(-1.5,1.5);
        //     hKaonTPC1D->Draw("pe");
        //     lat.DrawLatex(0.5, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));
        //     TLine *line0 = new TLine(0, 0, 0, hKaonTPC1D->GetMaximum());
        //     line0->SetLineColor(kRed);
        //     line0->SetLineStyle(2);
        //     line0->Draw("same");
        // }
        // if (isSavePlots)
        //     cnSigmaTPCKaonMinusvsPt->SaveAs(outputPathQA + "nSigmaTPC_KaonMinus_vs_Pt.png");

        // TCanvas *cnSigmaTPCTOFKaon = new TCanvas("cnSigmaTPCTOFKaon", "nSigmaTPC and TOF Kaon", 720, 720);
        // SetCanvasStyle(cnSigmaTPCTOFKaon, 0.15, 0.03, 0.05, 0.15);
        // SetHistoQA(hnSigmaTPCTOFKaon);
        // hnSigmaTPCTOFKaon->GetXaxis()->SetTitle("#it{n}_{#sigma} (TPC)");
        // hnSigmaTPCTOFKaon->GetYaxis()->SetTitle("#it{n}_{#sigma} (TOF)");
        // hnSigmaTPCTOFKaon->GetZaxis()->SetTitle("#it{p} (GeV/#it{c})");
        // hnSigmaTPCTOFKaon->Draw("BOX2Z");
        // if(isSavePlots)
        // cnSigmaTPCTOFKaon->SaveAs(outputPathQA + "3DSigmaTPCTOF_Kaon_vs_P.png");

        TCanvas *cnSigmaTPCTOFKaonvsPt = new TCanvas("cnSigmaTPCTOFKaonvsPt", "nSigmaTPC and TOF Kaon vs Pt", 1080, 720);
        SetCanvasStyle(cnSigmaTPCTOFKaonvsPt, 0.11, 0.06, 0.05, 0.10);
        cnSigmaTPCTOFKaonvsPt->Divide(2, 3);
        float pTbins_forCombinedTPCTOF[] = {0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0};
        int nPtBins_forCombinedTPCTOF = sizeof(pTbins_forCombinedTPCTOF) / sizeof(float) - 1;
        for (int i = 0; i < nPtBins_forCombinedTPCTOF; i++)
        {
            gPad->SetBottomMargin(0.12);
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.03);
            gPad->SetTopMargin(0.01);
            gPad->SetGrid(1, 1);
            cnSigmaTPCTOFKaonvsPt->cd(i + 1);
            float ptLow = pTbins_forCombinedTPCTOF[i];
            float ptHigh = pTbins_forCombinedTPCTOF[i + 1];
            int binLow = hnSigmaTPCTOFKaon->GetZaxis()->FindBin(ptLow + 0.001);
            int binHigh = hnSigmaTPCTOFKaon->GetZaxis()->FindBin(ptHigh - 0.001);
            hnSigmaTPCTOFKaon->GetZaxis()->SetRange(binLow, binHigh);
            TH2D *hKaonTPCtofpT = (TH2D *)hnSigmaTPCTOFKaon->Project3D("xy");
            hKaonTPCtofpT->SetName(Form("hKaonTPCtofpT_%.1f_%.1f", ptLow, ptHigh));
            SetHistoQA(hKaonTPCtofpT);
            // hKaonTPCtofpT->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));
            hKaonTPCtofpT->GetXaxis()->SetTitle("#it{n}_{#sigma} (TPC)");
            hKaonTPCtofpT->GetYaxis()->SetTitle("#it{n}_{#sigma} (TOF)");
            hKaonTPCtofpT->GetYaxis()->SetTitleOffset(1.3);
            hKaonTPCtofpT->Write(Form("hKaonTPCtofpT_%.1f_%.1f", ptLow, ptHigh));
            hKaonTPCtofpT->Draw("COLZ");
            lat.DrawLatex(0.5, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));
            TLine *line0 = new TLine(0, -3, 0, 3);
            line0->SetLineColor(kRed);
            line0->SetLineStyle(2);
            line0->Draw("same");
            TLine *line1 = new TLine(-3, 0, 3, 0);
            line1->SetLineColor(kRed);
            line1->SetLineStyle(2);
            line1->Draw("same");
        }
        if (isSavePlots)
            cnSigmaTPCTOFKaonvsPt->SaveAs(outputPathQA + "nSigmaTPC_TOF_Kaon_vs_Pt.png");
    }
}

Double_t pol2(double *x, double *par)
{
    double poly2 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
    return poly2;
}

Double_t pol3(double *x, double *par)
{
    double poly3 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
    return poly3;
}

Double_t relBW(double *x, double *par)
{
    double yield = par[0];
    double mass = par[1];
    double width = par[2];

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));
    return fit;
}

Double_t BW_pol2(double *x, double *par)
{
    return relBW(x, par) + pol2(x, &par[3]);
}

Double_t BW_pol3(double *x, double *par)
{
    return (relBW(x, par) + pol3(x, &par[3]));
}

//   TIter nextkey(file->GetListOfKeys());
//     TKey *key;
//     // print all keys in the file
//     while ((key = (TKey *)nextkey()))
//     {
//         std::cout << key->GetName() << "\n";
//     }

//     // now get keys inside the first key
//     TKey *key1 = (TKey *)file->GetListOfKeys()->At(0);
//     TDirectory *dir = (TDirectory *)key1->ReadObj();
//     TIter nextkey1(dir->GetListOfKeys());
//     TKey *key2;
//     // print all keys in the directory
//     while ((key2 = (TKey *)nextkey1()))
//     {
//         std::cout << key2->GetName() << "\n";
//     }