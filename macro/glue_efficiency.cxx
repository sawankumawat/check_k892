#include <iostream>
using namespace std;
#include "src/style.h"

void mcQAplots(TFile *f, string path);

void glue_efficiency()
{
    gStyle->SetOptStat(1110);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile *f = new TFile("../mc/LHC24l1/454849.root", "read");

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string histpath = "higher-mass-resonances/hMChists";
    string histpath_noycut = "higher-mass-resonances_f01710_noycut/hMChists";
    mcQAplots(f, histpath);

    THnSparseD *GenpT = (THnSparseD *)f->Get(Form("%s/Genf17102", histpath.c_str()));      // axis: multiplicity, pt, helicity angle
    THnSparseD *recpt1 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpath.c_str())); // axis: multiplicity, pt, mass, helicity angle

    if (GenpT == nullptr || recpt1 == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    double pTbins[] = {0, 1, 2, 3, 5, 7, 9, 12, 15, 20};
    int size = sizeof(pTbins) / sizeof(pTbins[0]);

    TH1D *hgenpt = GenpT->Projection(1);  // project on pt axis
    TH1D *hgencosThetaStar = GenpT->Projection(2); // project on helicity angle axis
    TH1D *hrecpt = recpt1->Projection(1); // project on pt axis
    TH1D *hreccosThetaStar = recpt1->Projection(3); // project on helicity angle axis
    TH1D *heff = new TH1D("heff", "Efficiency", size - 1, pTbins);
    for (int i = 0; i < hrecpt->GetNbinsX(); i++)
    {
        // get bin content accroding to pT bins and error according to bayesian method
        int lowptbin = hgenpt->GetXaxis()->FindBin(pTbins[i] + 0.01);
        int highptbin = hgenpt->GetXaxis()->FindBin(pTbins[i + 1] - 0.01);
        double genYield = hgenpt->Integral(lowptbin, highptbin);
        double recYield = hrecpt->Integral(lowptbin, highptbin);
        double recYieldError = TMath::Sqrt(((recYield + 1) / (genYield + 2)) * ((recYield + 2) / (genYield + 3) - (recYield + 1) / (genYield + 2)));
        if (genYield > 0)
        {
            heff->SetBinContent(i + 1, recYield / genYield);
            heff->SetBinError(i + 1, recYieldError);
        }
    }

    TCanvas *cefficiency = new TCanvas("", "Efficiency", 720, 720);
    SetCanvasStyle(cefficiency, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(heff);
    heff->GetYaxis()->SetTitle("Efficiency");
    heff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    heff->GetYaxis()->SetTitleOffset(1.5);
    heff->GetYaxis()->SetMaxDigits(3);
    heff->Draw();
    cefficiency->SaveAs("injected_mc_plots/efficiency.png");

    TCanvas *cgenpt = new TCanvas("", "Gen p_{T}", 720, 720);
    SetCanvasStyle(cgenpt, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hgenpt);
    hgenpt->GetYaxis()->SetTitle("Counts");
    hgenpt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hgenpt->GetYaxis()->SetTitleOffset(1.5);
    hgenpt->GetYaxis()->SetMaxDigits(3);
    hgenpt->Draw();
    cgenpt->SaveAs("injected_mc_plots/genpt.png");

    TCanvas *cgenCosThetaStar = new TCanvas("", "Helicity Angle Generated", 720, 720);
    SetCanvasStyle(cgenCosThetaStar, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hgencosThetaStar);
    hgencosThetaStar->GetYaxis()->SetTitle("Counts");
    hgencosThetaStar->GetXaxis()->SetTitle("Cos(#theta*)");
    hgencosThetaStar->GetYaxis()->SetTitleOffset(1.5);
    hgencosThetaStar->GetYaxis()->SetMaxDigits(3);
    hgencosThetaStar->GetXaxis()->SetRangeUser(-1, 1);
    hgencosThetaStar->SetMinimum(0);
    hgencosThetaStar->SetMaximum(1.2 * hgencosThetaStar->GetMaximum());
    hgencosThetaStar->Draw();
    cgenCosThetaStar->SaveAs("injected_mc_plots/genCosThetaStar.png");

    TCanvas *crecpt = new TCanvas("", "Rec p_{T}", 720, 720);
    SetCanvasStyle(crecpt, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hrecpt);
    hrecpt->GetYaxis()->SetTitle("Counts");
    hrecpt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hrecpt->GetYaxis()->SetTitleOffset(1.5);
    hrecpt->GetYaxis()->SetMaxDigits(3);
    hrecpt->Draw();
    crecpt->SaveAs("injected_mc_plots/recpt.png");

    TCanvas *crecCosThetaStar = new TCanvas("", "Helicity Angle Reconstructed", 720, 720);
    SetCanvasStyle(crecCosThetaStar, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hreccosThetaStar);
    hreccosThetaStar->GetYaxis()->SetTitle("Counts");
    hreccosThetaStar->GetXaxis()->SetTitle("Cos(#theta*)");
    hreccosThetaStar->GetYaxis()->SetTitleOffset(1.5);
    hreccosThetaStar->GetYaxis()->SetMaxDigits(3);
    hreccosThetaStar->GetXaxis()->SetRangeUser(-1, 1);
    hreccosThetaStar->SetMinimum(0);
    hreccosThetaStar->SetMaximum(1.2 * hreccosThetaStar->GetMaximum());
    hreccosThetaStar->Draw();
    crecCosThetaStar->SaveAs("injected_mc_plots/recCosThetaStar.png");
}

void mcQAplots(TFile *f, string path)
{
    THnSparseF *genRapidity2D = (THnSparseF *)f->Get(Form("%s/GenRapidity", path.c_str()));
    THnSparseF *genEta2D = (THnSparseF *)f->Get(Form("%s/GenEta", path.c_str()));
    TH1F *genRapidity = (TH1F *)genRapidity2D->Projection(0);
    TH1F *genEta = (TH1F *)genEta2D->Projection(0);
    TH1F *genPhi = (TH1F *)f->Get(Form("%s/GenPhi", path.c_str()));
    TH1F *recRapidity = (TH1F *)f->Get(Form("%s/RecRapidity", path.c_str()));
    TH1F *recEta = (TH1F *)f->Get(Form("%s/RecEta", path.c_str()));
    TH1F *recPhi = (TH1F *)f->Get(Form("%s/RecPhi", path.c_str()));

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