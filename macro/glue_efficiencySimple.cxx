#include <iostream>
using namespace std;
#include "src/style.h"

void mcQAplots(TFile *f, string path);

void glue_efficiencySimple()
{
    // gStyle->SetOptStat(1110);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TString savepath = "/home/sawan/Downloads/";
    // TString savepath = "/home/sawan/alice/practice/MCPlots/Gap4/";
    // TString savepath = "/home/sawan/check_k892/mc/LHC24l1/";

    // TFile *f = new TFile("/home/sawan/check_k892/mc/LHC24l1/463655.root", "read");
    // TFile *f = new TFile("/home/sawan/alice/practice/AnalysisResults_Gap4.root", "read");
    TFile *f = new TFile("/home/sawan/Downloads/AnalysisResults.root", "read");

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string histpathf01710 = "higher-mass-resonances/hMChists";
    string histpathf21525 = "higher-mass-resonances_f21525/hMChists";
    // string histpatha21320 = "higher-mass-resonances_a21320/hMChists";
    // string histpathf21270 = "higher-mass-resonances_f21270/hMChists";
    // mcQAplots(f, histpathf01710);
    // mcQAplots(f, histpathf21525);

    THnSparseD *GenpTf01710 = (THnSparseD *)f->Get(Form("%s/Genf17102", histpathf01710.c_str()));
    THnSparseD *recptf01710 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpathf01710.c_str()));

    THnSparseD *GenpTf21525 = (THnSparseD *)f->Get(Form("%s/Genf17102", histpathf21525.c_str()));
    THnSparseD *recptf21525 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpathf21525.c_str()));

    if (GenpTf01710 == nullptr || recptf01710 == nullptr || GenpTf21525 == nullptr || recptf21525 == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    float pTbins[] = {1, 2, 3, 5, 7, 10, 15};
    int sizePtBins = sizeof(pTbins) / sizeof(pTbins[0]);

    TH1D *hgenpt1710 = GenpTf01710->Projection(1); // project on pt axis
    TH1D *hrecpt1710 = recptf01710->Projection(1); // project on pt axis
    TH1D *hgenpt1525 = GenpTf21525->Projection(1); // project on pt axis
    TH1D *hrecpt1525 = recptf21525->Projection(1); // project on pt axis
    TH1D *heff1710 = new TH1D("heff1710", "Efficiency f0(1710)", sizePtBins - 1, pTbins);
    TH1D *heff1525 = new TH1D("heff1525", "Efficiency f2(1525)", sizePtBins - 1, pTbins);
    TH1D *hGenpt1710ptBins = new TH1D("hGenpt1710ptBins", "Generated pT distribution f0(1710)", 30, 0, 15);
    TH1D *hRecpt1710ptBins = new TH1D("hRecpt1710ptBins", "Reconstructed pT distribution f0(1710)", 30, 0, 15);
    TH1D *hGenpt1525ptBins = new TH1D("hGenpt1525ptBins", "Generated pT distribution f2(1525)", 30, 0, 15);
    TH1D *hRecpt1525ptBins = new TH1D("hRecpt1525ptBins", "Reconstructed pT distribution f2(1525)", 30, 0, 15);

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
        }

        double genYield1525 = hgenpt1525->Integral(lowptbin, highptbin);
        double recYield1525 = hrecpt1525->Integral(lowptbin, highptbin);
        double recYieldError1525 = TMath::Sqrt(abs(((recYield1525 + 1) / (genYield1525 + 2)) * ((recYield1525 + 2) / (genYield1525 + 3) - (recYield1525 + 1) / (genYield1525 + 2))));
        if (genYield1525 > 0)
        {
            heff1525->SetBinContent(i + 1, recYield1525 / genYield1525);
            heff1525->SetBinError(i + 1, recYieldError1525);
        }
    }

    for (int inum = 0; inum < 30; inum++)
    {
        int lowptbin = hgenpt1710->GetXaxis()->FindBin(inum / 2 + 0.01);
        int highptbin = hgenpt1710->GetXaxis()->FindBin((inum + 1) / 2 - 0.01);
        double genYield1710 = hgenpt1710->Integral(lowptbin, highptbin);
        hGenpt1710ptBins->SetBinContent(inum + 1, genYield1710);
        double recYield1710 = hrecpt1710->Integral(lowptbin, highptbin);
        hRecpt1710ptBins->SetBinContent(inum + 1, recYield1710);
    }

    TCanvas *cpT = new TCanvas("cpT", "pT distributions", 720, 720);
    SetCanvasStyle(cpT, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hgenpt1710);
    SetHistoQA(hrecpt1710);
    hgenpt1710->SetLineColor(kBlue);
    hgenpt1710->SetMarkerColor(kBlue);
    hgenpt1710->Rebin(2);
    hgenpt1710->SetMinimum(0);
    hgenpt1710->SetMaximum(hgenpt1710->GetMaximum() * 1.3);
    hgenpt1710->GetXaxis()->SetRangeUser(0, 16);
    hgenpt1710->Draw("HIST");
    hrecpt1710->SetLineColor(kRed);
    hrecpt1710->SetMarkerColor(kRed);
    hrecpt1710->Scale(10);
    hrecpt1710->Rebin(2);
    hrecpt1710->Draw("HIST SAME");

    TLegend *legpt = new TLegend(0.6, 0.8, 0.9, 0.92);
    legpt->SetBorderSize(0);
    legpt->SetTextSize(0.035);
    legpt->SetFillStyle(0);
    legpt->AddEntry((TObject *)0, "f_{0}(1710)", "");
    legpt->AddEntry(hgenpt1710, "Gen #it{p}_{T}", "l");
    legpt->AddEntry(hrecpt1710, "Rec #it{p}_{T} (x10)", "l");
    legpt->Draw();
    // cpT->SaveAs(savepath + "pT_dist_f0.png");

    TCanvas *cpt2 = new TCanvas("cpt2", "pT distributions in bins", 720, 720);
    SetCanvasStyle(cpt2, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hgenpt1525);
    SetHistoQA(hrecpt1525);
    hgenpt1525->SetLineColor(kBlue);
    hgenpt1525->SetMarkerColor(kBlue);
    hgenpt1525->Rebin(2);
    hgenpt1525->SetMinimum(0);
    hgenpt1525->SetMaximum(hgenpt1525->GetMaximum() * 1.3);
    hgenpt1525->Draw("HIST");
    hrecpt1525->SetLineColor(kRed);
    hrecpt1525->SetMarkerColor(kRed);
    hrecpt1525->Scale(10);
    hrecpt1525->Rebin(2);
    hrecpt1525->Draw("HIST SAME");

    TLegend *legpt2 = new TLegend(0.6, 0.8, 0.9, 0.92);
    legpt2->SetBorderSize(0);
    legpt2->SetTextSize(0.035);
    legpt2->SetFillStyle(0);
    legpt2->AddEntry((TObject *)0, "f_{2}(1525)", "");
    legpt2->AddEntry(hgenpt1525, "Gen #it{p}_{T}", "l");
    legpt2->AddEntry(hrecpt1525, "Rec #it{p}_{T} (x10)", "l");
    legpt2->Draw();
    // cpt2->SaveAs(savepath + "pT_dist_f2.png");

    TFile *foutput = new TFile(savepath + "efficiency.root", "RECREATE");

    TCanvas *cefficiency = new TCanvas("", "Efficiency", 720, 720);
    SetCanvasStyle(cefficiency, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(heff1710);
    heff1710->GetYaxis()->SetTitle("Efficiency");
    heff1710->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    heff1710->GetYaxis()->SetTitleOffset(1.5);
    heff1710->GetYaxis()->SetMaxDigits(3);
    heff1710->SetMaximum(71e-3);
    heff1710->SetMarkerSize(1.5);
    // heff1710->GetXaxis()->SetRangeUser(0, 20.5);
    heff1710->Write("hefficiencyf0");
    heff1710->Draw("pe");
    SetHistoQA(heff1525);
    heff1525->SetLineColor(kRed);
    heff1525->SetMarkerColor(kRed);
    heff1525->SetMarkerSize(1.5);
    heff1525->Draw("SAME pe");
    heff1525->Write("hefficiencyf2");
    TLegend *leg = new TLegend(0.6, 0.8, 0.9, 0.92);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->SetFillStyle(0);
    leg->AddEntry(heff1710, "f_{0}(1710)", "p");
    leg->AddEntry(heff1525, "f_{2}(1525)", "p");
    leg->Draw();
    // cefficiency->SaveAs(savepath + "efficiency_comparison.png");
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
        // c->SaveAs(Form("injected_mc_plots/%s.png", histNames[i].c_str()));
    }
}
