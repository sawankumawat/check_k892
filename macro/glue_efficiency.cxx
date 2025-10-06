#include <iostream>
using namespace std;
#include "src/style.h"

void mcQAplots(TFile *f, string path);

void glue_efficiency()
{
    gStyle->SetOptStat(1110);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile *f = new TFile("../mc/LHC24l1/463655.root", "read");

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string histpath = "higher-mass-resonances/hMChists";
    string histpath2 = "higher-mass-resonances_f21525/hMChists";
    mcQAplots(f, histpath);

    THnSparseD *GenpT = (THnSparseD *)f->Get(Form("%s/Genf17102", histpath.c_str()));        // axis: multiplicity, pt, helicity angle
    THnSparseD *recpt1 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpath.c_str()));    // axis: multiplicity, pt, mass, helicity angle
    THnSparseD *GenpTf2 = (THnSparseD *)f->Get(Form("%s/Genf17102", histpath2.c_str()));     // axis: multiplicity, pt, helicity angle
    THnSparseD *recpt1f2 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpath2.c_str())); // axis: multiplicity, pt, mass, helicity angle

    if (GenpT == nullptr || recpt1 == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    float pTbins[] = {0, 1, 2, 3, 5, 7, 10, 15, 20};
    int sizePtBins = sizeof(pTbins) / sizeof(pTbins[0]);
    float cosThetaBins[] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};
    int sizeCosThetaBins = sizeof(cosThetaBins) / sizeof(cosThetaBins[0]);

    TH1D *hgenpt = GenpT->Projection(1);                // project on pt axis
    TH1D *hgencosThetaStar = GenpT->Projection(2);      // project on helicity angle axis
    TH1D *hrecpt = recpt1->Projection(1);               // project on pt axis
    TH1D *hreccosThetaStar = recpt1->Projection(3);     // project on helicity angle axis
    TH1D *hgenptf2 = GenpTf2->Projection(1);            // project on pt axis
    TH1D *hreccosThetaStarf2 = recpt1f2->Projection(3); // project on helicity angle axis
    TH1D *hrecptf2 = recpt1f2->Projection(1);           // project on pt axis
    TH1D *hgencosThetaStarf2 = GenpTf2->Projection(2);  // project on helicity angle axis
    TH1D *heff = new TH1D("heff", "Efficiency", sizePtBins - 1, pTbins);
    TH1D *hefff2 = new TH1D("hefff2", "Efficiency for f21525", sizePtBins - 1, pTbins);
    TH1D *heffCosTheta = new TH1D("heffCosTheta", "Efficiency CosThetaStar", sizeCosThetaBins - 1, cosThetaBins);
    TH1D *heffCosThetaf2 = new TH1D("heffCosThetaf2", "Efficiency CosThetaStar for f21525", sizeCosThetaBins - 1, cosThetaBins);
    TH1D *heffpTfinerf0 = (TH1D *)hrecpt->Clone();
    heffpTfinerf0->Divide(hgenpt);
    TH1D *heffpTfinerf2 = (TH1D *)hrecptf2->Clone();
    heffpTfinerf2->Divide(hgenptf2);

    TH1F *fracUncertf0 = new TH1F("fracUncertf0", "Fractional Uncertainty f0(1710)", sizePtBins - 1, pTbins);
    TH1F *fracUncertf2 = new TH1F("fracUncertf2", "Fractional Uncertainty f2(1525)", sizePtBins - 1, pTbins);

    for (int i = 0; i < hrecpt->GetNbinsX(); i++)
    {
        // get bin content accroding to pT bins and error according to bayesian method
        int lowptbin = hgenpt->GetXaxis()->FindBin(pTbins[i] + 0.001);
        int highptbin = hgenpt->GetXaxis()->FindBin(pTbins[i + 1] - 0.001);
        double genYield = hgenpt->Integral(lowptbin, highptbin);
        double recYield = hrecpt->Integral(lowptbin, highptbin);
        double recYieldError = TMath::Sqrt(abs(((recYield + 1) / (genYield + 2)) * ((recYield + 2) / (genYield + 3) - (recYield + 1) / (genYield + 2))));
        if (genYield > 0)
        {
            heff->SetBinContent(i + 1, recYield / genYield);
            heff->SetBinError(i + 1, recYieldError);
            fracUncertf0->SetBinContent(i + 1, recYieldError / (recYield / genYield));
        }

        double genYieldf2 = hgenptf2->Integral(lowptbin, highptbin);
        double recYieldf2 = hrecptf2->Integral(lowptbin, highptbin);
        double recYieldErrorf2 = TMath::Sqrt(abs(((recYieldf2 + 1) / (genYieldf2 + 2)) * ((recYieldf2 + 2) / (genYieldf2 + 3) - (recYieldf2 + 1) / (genYieldf2 + 2))));
        if (genYieldf2 > 0)
        {
            hefff2->SetBinContent(i + 1, recYieldf2 / genYieldf2);
            hefff2->SetBinError(i + 1, recYieldErrorf2);
            fracUncertf2->SetBinContent(i + 1, recYieldErrorf2 / (recYieldf2 / genYieldf2));
        }
    }

    for (int i = 0; i < hreccosThetaStar->GetNbinsX(); i++)
    {
        // get bin content accroding to cosTheta bins and error according to bayesian method
        int lowcosbin = hgencosThetaStar->GetXaxis()->FindBin(cosThetaBins[i] + 0.001);
        int highcosbin = hgencosThetaStar->GetXaxis()->FindBin(cosThetaBins[i + 1] - 0.001);
        double genYieldCos = hgencosThetaStar->Integral(lowcosbin, highcosbin);
        double recYieldCos = hreccosThetaStar->Integral(lowcosbin, highcosbin);
        double recYieldErrorCos = TMath::Sqrt(abs(((recYieldCos + 1) / (genYieldCos + 2)) * ((recYieldCos + 2) / (genYieldCos + 3) - (recYieldCos + 1) / (genYieldCos + 2))));
        if (genYieldCos > 0)
        {
            heffCosTheta->SetBinContent(i + 1, recYieldCos / genYieldCos);
            heffCosTheta->SetBinError(i + 1, recYieldErrorCos);
        }

        double genYieldCosf2 = hgencosThetaStarf2->Integral(lowcosbin, highcosbin);
        double recYieldCosf2 = hreccosThetaStarf2->Integral(lowcosbin, highcosbin);
        double recYieldErrorCosf2 = TMath::Sqrt(abs(((recYieldCosf2 + 1) / (genYieldCosf2 + 2)) * ((recYieldCosf2 + 2) / (genYieldCosf2 + 3) - (recYieldCosf2 + 1) / (genYieldCosf2 + 2))));
        if (genYieldCosf2 > 0)
        {
            heffCosThetaf2->SetBinContent(i + 1, recYieldCosf2 / genYieldCosf2);
            heffCosThetaf2->SetBinError(i + 1, recYieldErrorCosf2);
        }
    }

    TFile *feff = new TFile("injected_mc_plots/efficiency.root", "RECREATE");
    heff->Write();
    heffCosTheta->Write();

    TCanvas *cefficiency = new TCanvas("", "Efficiency", 720, 720);
    SetCanvasStyle(cefficiency, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(heff);
    heff->GetYaxis()->SetTitle("Efficiency");
    heff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    heff->GetYaxis()->SetTitleOffset(1.5);
    heff->GetYaxis()->SetMaxDigits(3);
    heff->SetMaximum(heff->GetMaximum() * 1.25);
    heff->GetXaxis()->SetRangeUser(0, 20.5);
    heff->Draw();

    // TCanvas *cefficiencyf2 = new TCanvas("", "Efficiency for f21525", 720, 720);
    // SetCanvasStyle(cefficiencyf2, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hefff2);
    // hefff2->GetYaxis()->SetTitle("Efficiency");
    // hefff2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // hefff2->GetYaxis()->SetTitleOffset(1.5);
    // hefff2->GetYaxis()->SetMaxDigits(3);
    hefff2->SetLineColor(kRed);
    hefff2->SetMarkerStyle(21);
    hefff2->SetMarkerColor(kRed);
    hefff2->Draw("same");
    TLegend *leg = new TLegend(0.4, 0.4, 0.6, 0.6);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->AddEntry(heff, "f_{0}(1710)", "p");
    leg->AddEntry(hefff2, "f_{2}(1525)", "p");
    leg->Draw();
    cefficiency->SaveAs("injected_mc_plots/efficiency.png");

    TCanvas *cfracUncert = new TCanvas("", "Fractional Uncertainty", 720, 720);
    SetCanvasStyle(cfracUncert, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(fracUncertf0);
    fracUncertf0->GetYaxis()->SetTitle("Fractional Uncertainty");
    fracUncertf0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fracUncertf0->GetYaxis()->SetTitleOffset(1.5);
    fracUncertf0->GetYaxis()->SetMaxDigits(3);
    fracUncertf0->SetMaximum(fracUncertf0->GetMaximum() * 1.2);
    fracUncertf0->GetXaxis()->SetRangeUser(0, 20.5);
    fracUncertf0->Draw();
    SetHistoQA(fracUncertf2);
    fracUncertf2->SetLineColor(2);
    fracUncertf2->SetMarkerColor(2);
    fracUncertf2->Draw("same");
    leg->Draw();
    cfracUncert->SaveAs("injected_mc_plots/fractionalUncertainty.png");

    TCanvas *cefficiencyfiner = new TCanvas("", "Efficiency (Finer Bins)", 720, 720);
    SetCanvasStyle(cefficiencyfiner, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(heffpTfinerf0);
    heffpTfinerf0->GetYaxis()->SetTitle("Efficiency");
    heffpTfinerf0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    heffpTfinerf0->GetYaxis()->SetTitleOffset(1.5);
    heffpTfinerf0->GetYaxis()->SetMaxDigits(3);
    heffpTfinerf0->SetMaximum(heffpTfinerf0->GetMaximum() * 1.2);
    heffpTfinerf0->GetXaxis()->SetRangeUser(0, 20.5);
    heffpTfinerf0->Draw();
    SetHistoQA(heffpTfinerf2);
    heffpTfinerf2->SetLineColor(2);
    heffpTfinerf2->Draw("same");
    leg->Draw();
    cefficiencyfiner->SaveAs("injected_mc_plots/efficiency_finer.png");

    TCanvas *ceffCosTheta = new TCanvas("", "Efficiency CosThetaStar", 720, 720);
    SetCanvasStyle(ceffCosTheta, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(heffCosTheta);
    heffCosTheta->GetYaxis()->SetTitle("Efficiency");
    heffCosTheta->GetXaxis()->SetTitle("Cos(#theta*)");
    heffCosTheta->GetYaxis()->SetTitleOffset(1.5);
    heffCosTheta->GetYaxis()->SetMaxDigits(3);
    heffCosTheta->SetMaximum(heffCosTheta->GetMaximum() * 1.15);
    heffCosTheta->Draw();
    SetHistoQA(heffCosThetaf2);
    heffCosThetaf2->SetLineColor(kRed);
    heffCosThetaf2->SetMarkerColor(kRed);
    heffCosThetaf2->Draw("same");
    leg->Draw();
    ceffCosTheta->SaveAs("injected_mc_plots/efficiencyCosThetaf0.png");

    TCanvas *ceffCosThetaf2 = new TCanvas("", "Efficiency CosThetaStar for f21525", 720, 720);
    SetCanvasStyle(ceffCosThetaf2, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(heffCosThetaf2);
    heffCosThetaf2->GetYaxis()->SetTitle("Efficiency");
    heffCosThetaf2->GetXaxis()->SetTitle("Cos(#theta*)");
    heffCosThetaf2->GetYaxis()->SetTitleOffset(1.5);
    heffCosThetaf2->GetYaxis()->SetMaxDigits(3);
    heffCosThetaf2->Draw();
    ceffCosThetaf2->SaveAs("injected_mc_plots/efficiencyCosThetaf21525.png");

    TCanvas *cgenCosThetaStar = new TCanvas("", "Helicity Angle Generated", 720, 720);
    SetCanvasStyle(cgenCosThetaStar, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hgencosThetaStar);
    hgencosThetaStar->GetYaxis()->SetTitle("Counts");
    hgencosThetaStar->GetXaxis()->SetTitle("Generated cos(#theta*)");
    hgencosThetaStar->GetYaxis()->SetTitleOffset(1.5);
    hgencosThetaStar->GetYaxis()->SetMaxDigits(3);
    hgencosThetaStar->GetXaxis()->SetRangeUser(-1, 1);
    hgencosThetaStar->SetMinimum(0);
    hgencosThetaStar->SetMaximum(1.2 * hgencosThetaStar->GetMaximum());
    hgencosThetaStar->Draw();
    SetHistoQA(hgencosThetaStarf2);
    hgencosThetaStarf2->SetLineColor(kRed);
    hgencosThetaStarf2->Draw("same");
    leg->Clear();
    leg->AddEntry(hgencosThetaStar, "f_{0}(1710)", "l");
    leg->AddEntry(hgencosThetaStarf2, "f_{2}(1525)", "l");
    leg->Draw();
    cgenCosThetaStar->SaveAs("injected_mc_plots/genCosThetaStar.png");

    TCanvas *cgenpt = new TCanvas("", "Gen p_{T}", 720, 720);
    SetCanvasStyle(cgenpt, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hgenpt);
    hgenpt->GetYaxis()->SetTitle("Counts");
    hgenpt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hgenpt->GetYaxis()->SetTitleOffset(1.5);
    hgenpt->GetYaxis()->SetMaxDigits(3);
    hgenpt->GetXaxis()->SetRangeUser(0, 20.5);
    hgenpt->Draw();
    SetHistoQA(hgenptf2);
    hgenptf2->SetLineColor(2);
    hgenptf2->SetMarkerColor(2);
    hgenptf2->Draw("same");
    leg->SetY1(0.6);
    leg->SetY2(0.8);
    leg->Draw();
    cgenpt->SaveAs("injected_mc_plots/genpt.png");

    TCanvas *crecpt = new TCanvas("", "Rec p_{T}", 720, 720);
    SetCanvasStyle(crecpt, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hrecpt);
    hrecpt->GetYaxis()->SetTitle("Counts");
    hrecpt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hrecpt->GetYaxis()->SetTitleOffset(1.5);
    hrecpt->GetYaxis()->SetMaxDigits(3);
    hrecpt->GetXaxis()->SetRangeUser(0, 20.5);
    hrecpt->SetMaximum(1.2 * hrecpt->GetMaximum());
    hrecpt->Draw();
    SetHistoQA(hrecptf2);
    hrecptf2->SetLineColor(2);
    hrecptf2->SetMarkerColor(2);
    hrecptf2->Draw("same");
    leg->SetY1(0.4);
    leg->SetY2(0.6);
    leg->Draw();
    crecpt->SaveAs("injected_mc_plots/recpt.png");

    TCanvas *crecCosThetaStar = new TCanvas("", "Helicity Angle Reconstructed", 720, 720);
    SetCanvasStyle(crecCosThetaStar, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hreccosThetaStar);
    hreccosThetaStar->GetYaxis()->SetTitle("Counts");
    hreccosThetaStar->GetXaxis()->SetTitle("Reconstructed cos(#theta*)");
    hreccosThetaStar->GetYaxis()->SetTitleOffset(1.5);
    hreccosThetaStar->GetYaxis()->SetMaxDigits(3);
    hreccosThetaStar->GetXaxis()->SetRangeUser(-1, 1);
    hreccosThetaStar->SetMinimum(0);
    hreccosThetaStar->SetMaximum(1.2 * hreccosThetaStar->GetMaximum());
    hreccosThetaStar->Draw();
    SetHistoQA(hreccosThetaStarf2);
    hreccosThetaStarf2->SetLineColor(2);
    hreccosThetaStarf2->Draw("same");
    leg->Draw();
    crecCosThetaStar->SaveAs("injected_mc_plots/recCosThetaStar.png");
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