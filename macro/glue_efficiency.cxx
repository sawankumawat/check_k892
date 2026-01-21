#include <iostream>
using namespace std;
#include "src/style.h"

void mcQAplots(TFile *f, string path);

void glue_efficiency()
{
    gStyle->SetOptStat(1110);
    // gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // TFile *f = new TFile("../mc/LHC24l1/590796.root", "read");
    TFile *f = new TFile("/home/sawan/Downloads/AnalysisResults.root", "read");
    // TFile *f = new TFile("/home/sawan/alice/practice/AnalysisResults.root", "read");

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string histpathf01710 = "higher-mass-resonances_id45289/hMChists";
    string histpathf21525 = "higher-mass-resonances_f21525_id45289/hMChists";
    // string histpatha21320 = "higher-mass-resonances_a21320/hMChists";
    // string histpathf21270 = "higher-mass-resonances_f21270/hMChists";
    // mcQAplots(f, histpathf01710);

    THnSparseD *GenpTf01710 = (THnSparseD *)f->Get(Form("%s/Genf17102", histpathf01710.c_str()));     // axis: multiplicity, pt, helicity angle
    THnSparseD *recptf01710 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpathf01710.c_str()));  // axis: multiplicity, pt, mass, helicity angle
    THnSparseD *GenpTf21525 = (THnSparseD *)f->Get(Form("%s/Genf17102", histpathf21525.c_str()));     // axis: multiplicity, pt, helicity angle
    THnSparseD *recpt1f21525 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpathf21525.c_str())); // axis: multiplicity, pt, mass, helicity angle
    // THnSparseD *GenpTa21320 = (THnSparseD *)f->Get(Form("%s/Genf17102", histpatha21320.c_str()));     // axis: multiplicity, pt, helicity angle
    // THnSparseD *recpta21320 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpatha21320.c_str()));  // axis: multiplicity, pt, mass, helicity angle
    // THnSparseD *GenpTf21270 = (THnSparseD *)f->Get(Form("%s/Genf17102", histpathf21270.c_str()));     // axis: multiplicity, pt, helicity angle
    // THnSparseD *recptf21270 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpathf21270.c_str()));  // axis: multiplicity, pt, mass, helicity angle

    if (GenpTf01710 == nullptr || recptf01710 == nullptr || GenpTf21525 == nullptr || recpt1f21525 == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    float pTbins[] = {1, 2, 3, 5, 7, 10, 15};
    int sizePtBins = sizeof(pTbins) / sizeof(pTbins[0]);
    float cosThetaBins[] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};
    int sizeCosThetaBins = sizeof(cosThetaBins) / sizeof(cosThetaBins[0]);

    TH1D *hgenpt1710 = GenpTf01710->Projection(1);            // project on pt axis
    TH1D *hgencosThetaStar1710 = GenpTf01710->Projection(2);  // project on helicity angle axis
    TH1D *hrecpt1710 = recptf01710->Projection(1);            // project on pt axis
    TH1D *hreccosThetaStar1710 = recptf01710->Projection(3);  // project on helicity angle axis
    TH1D *hgenpt1525 = GenpTf21525->Projection(1);            // project on pt axis
    TH1D *hreccosThetaStar1525 = recpt1f21525->Projection(3); // project on helicity angle axis
    TH1D *hrecpt1525 = recpt1f21525->Projection(1);           // project on pt axis
    TH1D *hgencosThetaStar1525 = GenpTf21525->Projection(2);  // project on helicity angle axis
    // TH1D *hgenpt1320 = GenpTa21320->Projection(1);            // project on pt axis
    // TH1D *hreccosThetaStar1320 = recpta21320->Projection(3);  // project on helicity angle axis
    // TH1D *hrecpt1320 = recpta21320->Projection(1);           // project on pt axis
    // TH1D *hgenpt1270 = GenpTf21270->Projection(1);            // project on pt axis
    // TH1D *hreccosThetaStar1270 = recptf21270->Projection(3);  // project on helicity angle axis
    // TH1D *hrecpt1270 = recptf21270->Projection(1);           // project on pt axis
    TH1D *heff1710 = new TH1D("heff1710", "Efficiency", sizePtBins - 1, pTbins);
    TH1D *heff1525 = new TH1D("heff1525", "Efficiency for f21525", sizePtBins - 1, pTbins);
    // TH1D *heff1320 = new TH1D("heff1320", "Efficiency for a21320", sizePtBins - 1, pTbins);
    // TH1D *heff1270 = new TH1D("heff1270", "Efficiency for f21270", sizePtBins - 1, pTbins);
    // TH1D *heffCosTheta1710 = new TH1D("heffCosTheta1710", "Efficiency CosThetaStar", sizeCosThetaBins - 1, cosThetaBins);
    // TH1D *heffCosTheta1525 = new TH1D("heffCosTheta1525", "Efficiency CosThetaStar for f21525", sizeCosThetaBins - 1, cosThetaBins);
    // TH1D *heffCosTheta1320 = new TH1D("heffCosTheta1320", "Efficiency CosThetaStar for a21320", sizeCosThetaBins - 1, cosThetaBins);
    // TH1D *heffCosTheta1270 = new TH1D("heffCosTheta1270", "Efficiency CosThetaStar for f21270", sizeCosThetaBins - 1, cosThetaBins);
    TH1D *heffpTfiner1710 = (TH1D *)hrecpt1710->Clone();
    heffpTfiner1710->Divide(hgenpt1710);
    TH1D *heffpTfiner1525 = (TH1D *)hrecpt1525->Clone();
    heffpTfiner1525->Divide(hgenpt1525);

    TH1D *hEventLoss1710 = new TH1D("hEventLoss1710", "Event Loss for f0(1710)", sizePtBins - 1, pTbins);
    TH1D *hSignalLoss1710 = new TH1D("hSignalLoss1710", "Signal Loss for f0(1710)", sizePtBins - 1, pTbins);

    TH1F *fracUncert1710 = new TH1F("fracUncert1710", "Fractional Uncertainty f0(1710)", sizePtBins - 1, pTbins);
    TH1F *fracUncert1525 = new TH1F("fracUncert1525", "Fractional Uncertainty f2(1525)", sizePtBins - 1, pTbins);

    // Event loss histograms
    TH1F *hAllGenColl = (TH1F *)f->Get(Form("%s/MCcorrections/MultiplicityGen", histpathf01710.c_str()));
    TH1F *hAllGenColl1Rec = (TH1F *)f->Get(Form("%s/MCcorrections/MultiplicityRec", histpathf01710.c_str()));

    if (hAllGenColl == nullptr || hAllGenColl1Rec == nullptr)
    {
        cout << "Error reading event loss histograms" << endl;
        cout << "Path is " << Form("%s/MCcorrections/MultiplicityRec", histpathf01710.c_str()) << endl;
        return;
    }

    // Signal loss histograms
    TH2F *hAllGenKstar = (TH2F *)f->Get(Form("%s/MCcorrections/hSignalLossDenominator", histpathf01710.c_str()));
    TH2F *hAllGenKstar1Rec = (TH2F *)f->Get(Form("%s/MCcorrections/hSignalLossNumerator", histpathf01710.c_str()));

    for (int i = 0; i < sizePtBins - 1; i++)
    {
        // get bin content accroding to pT bins and error according to bayesian method
        int lowptbin = hgenpt1710->GetXaxis()->FindBin(pTbins[i] + 0.01);
        int highptbin = hgenpt1710->GetXaxis()->FindBin(pTbins[i + 1] - 0.01);
        double genYield1710 = hgenpt1710->Integral(lowptbin, highptbin);
        double recYield1710 = hrecpt1710->Integral(lowptbin, highptbin);
        double recYieldError1710 = TMath::Sqrt(abs(((recYield1710 + 1) / (genYield1710 + 2)) * ((recYield1710 + 2) / (genYield1710 + 3) - (recYield1710 + 1) / (genYield1710 + 2))));
        if (genYield1710 > 0)
        {
            heff1710->SetBinContent(i + 1, recYield1710 / genYield1710);
            heff1710->SetBinError(i + 1, recYieldError1710);
            fracUncert1710->SetBinContent(i + 1, recYieldError1710 / (recYield1710 / genYield1710));
        }

        double genYield1525 = hgenpt1525->Integral(lowptbin, highptbin);
        double recYield1525 = hrecpt1525->Integral(lowptbin, highptbin);
        double recYieldError1525 = TMath::Sqrt(abs(((recYield1525 + 1) / (genYield1525 + 2)) * ((recYield1525 + 2) / (genYield1525 + 3) - (recYield1525 + 1) / (genYield1525 + 2))));
        if (genYield1525 > 0)
        {
            heff1525->SetBinContent(i + 1, recYield1525 / genYield1525);
            heff1525->SetBinError(i + 1, recYieldError1525);
            fracUncert1525->SetBinContent(i + 1, recYieldError1525 / (recYield1525 / genYield1525));
        }
        int multlbinlow = hAllGenKstar->GetXaxis()->FindBin(0.01);
        int multlbinhigh = hAllGenKstar->GetXaxis()->FindBin(100.0 - 0.01);
        
        // Event loss
        double eventLossNum1710 = hAllGenColl1Rec->Integral(multlbinlow, multlbinhigh);
        double eventLossDen1710 = hAllGenColl->Integral(multlbinlow, multlbinhigh);
        double eventLoss1710 = eventLossNum1710 / eventLossDen1710;
        hEventLoss1710->SetBinContent(i + 1, eventLoss1710);

        // Signal loss
        TH1D *hSignalLossDenPt1710 = hAllGenKstar->ProjectionX(Form("SignalLossDenPt1710_%d", i), multlbinlow, multlbinhigh);
        TH1D *hSignalLossNumPt1710 = hAllGenKstar1Rec->ProjectionX(Form("SignalLossNumPt1710_%d", i), multlbinlow, multlbinhigh);
        double signalLossDen1710 = hSignalLossDenPt1710->Integral(hSignalLossDenPt1710->GetXaxis()->FindBin(pTbins[i] + 0.01), hSignalLossDenPt1710->GetXaxis()->FindBin(pTbins[i + 1] - 0.01));
        double signalLossNum1710 = hSignalLossNumPt1710->Integral(hSignalLossNumPt1710->GetXaxis()->FindBin(pTbins[i] + 0.01), hSignalLossNumPt1710->GetXaxis()->FindBin(pTbins[i + 1] - 0.01));
        double signalLoss1710 = signalLossNum1710 / signalLossDen1710;
        hSignalLoss1710->SetBinContent(i + 1, signalLoss1710);

        cout << "pT bin " << i << " event loss " << eventLoss1710 << " signal loss " << signalLoss1710 << endl;

        // double genYield1320 = hgenpt1320->Integral(lowptbin, highptbin);
        // double recYield1320 = hrecpt1320->Integral(lowptbin, highptbin);
        // double recYieldError1320 = TMath::Sqrt(abs(((recYield1320 + 1) / (genYield1320 + 2)) * ((recYield1320 + 2) / (genYield1320 + 3) - (recYield1320 + 1) / (genYield1320 + 2))));
        // if (genYield1320 > 0)
        // {
        //     heff1320->SetBinContent(i + 1, recYield1320 / genYield1320);
        //     heff1320->SetBinError(i + 1, recYieldError1320);
        // }
        // double genYield1270 = hgenpt1270->Integral(lowptbin, highptbin);
        // double recYield1270 = hrecpt1270->Integral(lowptbin, highptbin);
        // double recYieldError1270 = TMath::Sqrt(abs(((recYield1270 + 1) / (genYield1270 + 2)) * ((recYield1270 + 2) / (genYield1270 + 3) - (recYield1270 + 1) / (genYield1270 + 2))));
        // if (genYield1270 > 0)
        // {
        //     heff1270->SetBinContent(i + 1, recYield1270 / genYield1270);
        //     heff1270->SetBinError(i + 1, recYieldError1270);
        // }
    }

    // for (int i = 0; i < hreccosThetaStar1710->GetNbinsX(); i++)
    // {
    //     // get bin content accroding to cosTheta bins and error according to bayesian method
    //     int lowcosbin = hgencosThetaStar1710->GetXaxis()->FindBin(cosThetaBins[i] + 0.001);
    //     int highcosbin = hgencosThetaStar1710->GetXaxis()->FindBin(cosThetaBins[i + 1] - 0.001);
    //     double genYieldCos = hgencosThetaStar1710->Integral(lowcosbin, highcosbin);
    //     double recYieldCos = hreccosThetaStar1710->Integral(lowcosbin, highcosbin);
    //     double recYieldErrorCos = TMath::Sqrt(abs(((recYieldCos + 1) / (genYieldCos + 2)) * ((recYieldCos + 2) / (genYieldCos + 3) - (recYieldCos + 1) / (genYieldCos + 2))));
    //     if (genYieldCos > 0)
    //     {
    //         heffCosTheta1710->SetBinContent(i + 1, recYieldCos / genYieldCos);
    //         heffCosTheta1710->SetBinError(i + 1, recYieldErrorCos);
    //     }

    //     double genYieldCosf2 = hgencosThetaStar1525->Integral(lowcosbin, highcosbin);
    //     double recYieldCosf2 = hreccosThetaStar1525->Integral(lowcosbin, highcosbin);
    //     double recYieldErrorCosf2 = TMath::Sqrt(abs(((recYieldCosf2 + 1) / (genYieldCosf2 + 2)) * ((recYieldCosf2 + 2) / (genYieldCosf2 + 3) - (recYieldCosf2 + 1) / (genYieldCosf2 + 2))));
    //     if (genYieldCosf2 > 0)
    //     {
    //         heffCosTheta1525->SetBinContent(i + 1, recYieldCosf2 / genYieldCosf2);
    //         heffCosTheta1525->SetBinError(i + 1, recYieldErrorCosf2);
    //     }
    // }

    // TFile *feff = new TFile("injected_mc_plots/efficiency.root", "RECREATE");
    // heff1710->Write();
    // heffCosTheta1710->Write();

    TCanvas *cefficiency = new TCanvas("", "Efficiency", 720, 720);
    SetCanvasStyle(cefficiency, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(heff1710);
    heff1710->GetYaxis()->SetTitle("Efficiency");
    heff1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    heff1710->GetYaxis()->SetTitleOffset(1.5);
    heff1710->GetYaxis()->SetMaxDigits(3);
    // heff1710->SetMaximum(heff1710->GetMaximum() * 1.6);
    heff1710->SetMaximum(61e-3);
    heff1710->GetXaxis()->SetRangeUser(0, 20.5);
    heff1710->Draw("HIST");
    // TCanvas *cefficiencyf2 = new TCanvas("", "Efficiency for f21525", 720, 720);
    // SetCanvasStyle(cefficiencyf2, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(heff1525);
    // heff1525->GetYaxis()->SetTitle("Efficiency");
    // heff1525->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // heff1525->GetYaxis()->SetTitleOffset(1.5);
    // heff1525->GetYaxis()->SetMaxDigits(3);
    heff1525->SetLineColor(kRed);
    heff1525->SetMarkerStyle(21);
    heff1525->SetMarkerColor(kRed);
    heff1525->Draw("HIST same");
    // SetHistoQA(heff1320);
    // heff1320->SetLineColor(kGreen + 2);
    // heff1320->SetMarkerStyle(22);
    // heff1320->SetMarkerColor(kGreen + 2);
    // heff1320->Draw("HIST same");
    // SetHistoQA(heff1270);
    // heff1270->SetLineColor(kBlue);
    // heff1270->SetMarkerStyle(23);
    // heff1270->SetMarkerColor(kBlue);
    // heff1270->Draw("HIST same");
    TLegend *leg = new TLegend(0.45, 0.3, 0.65, 0.45);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->AddEntry(heff1710, "f_{0}(1710)", "p");
    leg->AddEntry(heff1525, "f_{2}(1525)", "p");
    // leg->AddEntry(heff1320, "a_{2}(1320)", "p");
    // leg->AddEntry(heff1270, "f_{0}(1270)", "p");
    leg->Draw();
    cefficiency->SaveAs("injected_mc_plots/efficiency.png");

    TCanvas *ceventSigLoss = new TCanvas("", "Event and Signal Loss", 720, 720);
    SetCanvasStyle(ceventSigLoss, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hEventLoss1710);
    hEventLoss1710->GetYaxis()->SetTitle("Event and signal loss");
    hEventLoss1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hEventLoss1710->GetYaxis()->SetTitleOffset(1.5);
    hEventLoss1710->GetYaxis()->SetMaxDigits(3);
    // hEventLoss1710->SetMaximum(1.8);
    // hEventLoss1710->GetXaxis()->SetRangeUser(0, 20.5);
    hEventLoss1710->Draw("HIST");
    SetHistoQA(hSignalLoss1710);
    hSignalLoss1710->SetLineColor(kRed);
    hSignalLoss1710->SetMarkerColor(kRed);
    hSignalLoss1710->Draw("p same");

    // TCanvas *cfracUncert = new TCanvas("", "Fractional Uncertainty", 720, 720);
    // SetCanvasStyle(cfracUncert, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(fracUncert1710);
    // fracUncert1710->GetYaxis()->SetTitle("Fractional Uncertainty");
    // fracUncert1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // fracUncert1710->GetYaxis()->SetTitleOffset(1.5);
    // fracUncert1710->GetYaxis()->SetMaxDigits(3);
    // fracUncert1710->SetMaximum(fracUncert1710->GetMaximum() * 1.2);
    // fracUncert1710->GetXaxis()->SetRangeUser(0, 20.5);
    // fracUncert1710->Draw();
    // SetHistoQA(fracUncert1525);
    // fracUncert1525->SetLineColor(2);
    // fracUncert1525->SetMarkerColor(2);
    // fracUncert1525->Draw("same");
    // leg->Draw();
    // cfracUncert->SaveAs("injected_mc_plots/fractionalUncertainty.png");

    // TCanvas *cefficiencyfiner = new TCanvas("", "Efficiency (Finer Bins)", 720, 720);
    // SetCanvasStyle(cefficiencyfiner, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(heffpTfiner1710);
    // heffpTfiner1710->GetYaxis()->SetTitle("Efficiency");
    // heffpTfiner1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // heffpTfiner1710->GetYaxis()->SetTitleOffset(1.5);
    // heffpTfiner1710->GetYaxis()->SetMaxDigits(3);
    // heffpTfiner1710->SetMaximum(heffpTfiner1710->GetMaximum() * 1.2);
    // heffpTfiner1710->GetXaxis()->SetRangeUser(0, 20.5);
    // heffpTfiner1710->Draw();
    // SetHistoQA(heffpTfiner1525);
    // heffpTfiner1525->SetLineColor(2);
    // heffpTfiner1525->Draw("same");
    // leg->Draw();
    // cefficiencyfiner->SaveAs("injected_mc_plots/efficiency_finer.png");

    // TCanvas *ceffCosTheta = new TCanvas("", "Efficiency CosThetaStar", 720, 720);
    // SetCanvasStyle(ceffCosTheta, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(heffCosTheta1710);
    // heffCosTheta1710->GetYaxis()->SetTitle("Efficiency");
    // heffCosTheta1710->GetXaxis()->SetTitle("Cos(#theta*)");
    // heffCosTheta1710->GetYaxis()->SetTitleOffset(1.5);
    // heffCosTheta1710->GetYaxis()->SetMaxDigits(3);
    // heffCosTheta1710->SetMaximum(heffCosTheta1710->GetMaximum() * 1.15);
    // heffCosTheta1710->Draw();
    // SetHistoQA(heffCosTheta1525);
    // heffCosTheta1525->SetLineColor(kRed);
    // heffCosTheta1525->SetMarkerColor(kRed);
    // heffCosTheta1525->Draw("same");
    // leg->Draw();
    // ceffCosTheta->SaveAs("injected_mc_plots/efficiencyCosThetaf0.png");

    // TCanvas *ceffCosThetaf2 = new TCanvas("", "Efficiency CosThetaStar for f21525", 720, 720);
    // SetCanvasStyle(ceffCosThetaf2, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(heffCosTheta1525);
    // heffCosTheta1525->GetYaxis()->SetTitle("Efficiency");
    // heffCosTheta1525->GetXaxis()->SetTitle("Cos(#theta*)");
    // heffCosTheta1525->GetYaxis()->SetTitleOffset(1.5);
    // heffCosTheta1525->GetYaxis()->SetMaxDigits(3);
    // heffCosTheta1525->Draw();
    // ceffCosThetaf2->SaveAs("injected_mc_plots/efficiencyCosThetaf21525.png");

    // TCanvas *cgenCosThetaStar = new TCanvas("", "Helicity Angle Generated", 720, 720);
    // SetCanvasStyle(cgenCosThetaStar, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(hgencosThetaStar1710);
    // hgencosThetaStar1710->GetYaxis()->SetTitle("Counts");
    // hgencosThetaStar1710->GetXaxis()->SetTitle("Generated cos(#theta*)");
    // hgencosThetaStar1710->GetYaxis()->SetTitleOffset(1.5);
    // hgencosThetaStar1710->GetYaxis()->SetMaxDigits(3);
    // hgencosThetaStar1710->GetXaxis()->SetRangeUser(-1, 1);
    // hgencosThetaStar1710->SetMinimum(0);
    // hgencosThetaStar1710->SetMaximum(1.2 * hgencosThetaStar1710->GetMaximum());
    // hgencosThetaStar1710->Draw();
    // SetHistoQA(hgencosThetaStar1525);
    // hgencosThetaStar1525->SetLineColor(kRed);
    // hgencosThetaStar1525->Draw("same");
    // leg->Clear();
    // leg->AddEntry(hgencosThetaStar1710, "f_{0}(1710)", "l");
    // leg->AddEntry(hgencosThetaStar1525, "f_{2}(1525)", "l");
    // leg->Draw();
    // cgenCosThetaStar->SaveAs("injected_mc_plots/genCosThetaStar.png");

    // TCanvas *cgenpt = new TCanvas("", "Gen p_{T}", 720, 720);
    // SetCanvasStyle(cgenpt, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(hgenpt1710);
    // hgenpt1710->GetYaxis()->SetTitle("Counts");
    // hgenpt1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // hgenpt1710->GetYaxis()->SetTitleOffset(1.5);
    // hgenpt1710->GetYaxis()->SetMaxDigits(3);
    // hgenpt1710->GetXaxis()->SetRangeUser(0, 20.5);
    // hgenpt1710->Draw();
    // SetHistoQA(hgenpt1525);
    // hgenpt1525->SetLineColor(2);
    // hgenpt1525->SetMarkerColor(2);
    // hgenpt1525->Draw("same");
    // leg->SetY1(0.6);
    // leg->SetY2(0.8);
    // leg->Draw();
    // cgenpt->SaveAs("injected_mc_plots/genptf01710.png");

    // TCanvas *crecpt = new TCanvas("", "Rec p_{T}", 720, 720);
    // SetCanvasStyle(crecpt, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(hrecpt1710);
    // hrecpt1710->GetYaxis()->SetTitle("Counts");
    // hrecpt1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // hrecpt1710->GetYaxis()->SetTitleOffset(1.5);
    // hrecpt1710->GetYaxis()->SetMaxDigits(3);
    // hrecpt1710->GetXaxis()->SetRangeUser(0, 20.5);
    // hrecpt1710->SetMaximum(1.2 * hrecpt1710->GetMaximum());
    // hrecpt1710->Draw();
    // SetHistoQA(hrecpt1525);
    // hrecpt1525->SetLineColor(2);
    // hrecpt1525->SetMarkerColor(2);
    // hrecpt1525->Draw("same");
    // SetHistoQA(hrecpt1320);
    // hrecpt1320->SetLineColor(kGreen + 2);
    // hrecpt1320->SetMarkerColor(kGreen + 2);
    // hrecpt1320->Draw("same");
    // SetHistoQA(hrecpt1270);
    // hrecpt1270->SetLineColor(kBlue);
    // hrecpt1270->SetMarkerColor(kBlue);
    // hrecpt1270->Draw("same");
    // leg->SetY1(0.4);
    // leg->SetY2(0.6);
    // leg->Draw();
    // // crecpt->SaveAs("injected_mc_plots/recpt.png");

    // TCanvas *crecCosThetaStar = new TCanvas("", "Helicity Angle Reconstructed", 720, 720);
    // SetCanvasStyle(crecCosThetaStar, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(hreccosThetaStar1710);
    // hreccosThetaStar1710->GetYaxis()->SetTitle("Counts");
    // hreccosThetaStar1710->GetXaxis()->SetTitle("Reconstructed cos(#theta*)");
    // hreccosThetaStar1710->GetYaxis()->SetTitleOffset(1.5);
    // hreccosThetaStar1710->GetYaxis()->SetMaxDigits(3);
    // hreccosThetaStar1710->GetXaxis()->SetRangeUser(-1, 1);
    // hreccosThetaStar1710->SetMinimum(0);
    // hreccosThetaStar1710->SetMaximum(1.2 * hreccosThetaStar1710->GetMaximum());
    // hreccosThetaStar1710->Draw();
    // SetHistoQA(hreccosThetaStar1525);
    // hreccosThetaStar1525->SetLineColor(2);
    // hreccosThetaStar1525->Draw("same");
    // leg->Draw();
    // crecCosThetaStar->SaveAs("injected_mc_plots/recCosThetaStar.png");
}

void mcQAplots(TFile *f, string path)
{
    THnSparseF *genRapidity2D = (THnSparseF *)f->Get(Form("%s/GenRapidity2", path.c_str()));
    THnSparseF *genEta2D = (THnSparseF *)f->Get(Form("%s/GenEta2", path.c_str()));
    TH1F *genRapidity = (TH1F *)genRapidity2D->Projection(0);
    genRapidity->Rebin(3);
    TH1F *genEta = (TH1F *)genEta2D->Projection(0);
    genEta->Rebin(3);
    TH1F *genPhi = (TH1F *)f->Get(Form("%s/GenPhi", path.c_str()));
    genPhi->Rebin(5);
    TH1F *recRapidity = (TH1F *)f->Get(Form("%s/RecRapidity2", path.c_str()));
    recRapidity->Rebin(3);
    TH1F *recEta = (TH1F *)f->Get(Form("%s/RecEta2", path.c_str()));
    recEta->Rebin(3);
    TH1F *recPhi = (TH1F *)f->Get(Form("%s/RecPhi", path.c_str()));
    recPhi->Rebin(5);

    if (genRapidity == nullptr || genEta == nullptr || genPhi == nullptr ||
        recRapidity == nullptr || recEta == nullptr || recPhi == nullptr)
    {
        cout << "Error reading histograms" << endl;
        return;
    }

    string histNames[] = {
        "genRapidity", "genEta", "genPhi",
        "recRapidity", "recEta", "recPhi"};
    string histXaxis[] = {
        "y", "#eta", "#Phi",
        "y", "#eta", "#Phi"};

    TH1F *histograms[] = {
        genRapidity, genEta, genPhi,
        recRapidity, recEta, recPhi};

    for (int i = 0; i < 6; i++)
    {
        TCanvas *c = new TCanvas("", histNames[i].c_str(), 720, 720);
        SetCanvasStyle(c, 0.15, 0.05, 0.05, 0.15);
        histograms[i]->SetName(histNames[i].c_str());
        SetHistoQA(histograms[i]);
        histograms[i]->GetXaxis()->SetTitle(histXaxis[i].c_str());
        histograms[i]->GetYaxis()->SetTitle("Counts");
        histograms[i]->GetYaxis()->SetTitleOffset(1.5);
        histograms[i]->GetYaxis()->SetMaxDigits(3);
        if (i == 1)
            histograms[i]->GetYaxis()->SetRangeUser(0, histograms[i]->GetMaximum() * 1.3);
        histograms[i]->Draw();
        c->SaveAs(Form("injected_mc_plots/%s.png", histNames[i].c_str()));
    }
}