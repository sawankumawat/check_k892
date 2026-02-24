#include <TGenPhaseSpace.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "src/style.h"
#include <Math/VectorUtil.h>
#include "spectra/YieldMean.C"

using namespace std;
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void compare_efficiency()
{
    gStyle->SetOptStat(0);
    TString path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra/";
    TFile *fEffToy = new TFile("efficiencyPtSmear.root", "read");
    TFile *fEffMC = new TFile(path + "spectra_Default2.root", "read");
    if (fEffToy->IsZombie() || fEffMC->IsZombie())
    {
        cout << "Error opening pt smearing results file" << endl;
        return;
    }

    TH1F *hEff_1710_toy = (TH1F *)fEffToy->Get("hEff_f0");
    TH1F *hEff_1525_toy = (TH1F *)fEffToy->Get("hEff_f2");
    TH1F *hEff_1710_MC = (TH1F *)fEffMC->Get("hefficiencyf0");
    TH1F *hEff_1525_MC = (TH1F *)fEffMC->Get("hefficiencyf2");
    if (hEff_1710_toy == nullptr || hEff_1525_toy == nullptr)
    {
        cout << "Error: Histograms not found in the file." << endl;
        return;
    }

    cout << "Number of bins in toy model efficiency " << hEff_1710_toy->GetNbinsX() << " ,number of bins in MC efficiency " << hEff_1710_MC->GetNbinsX() << endl;

    TCanvas *cEffResults = new TCanvas("cEffResults", "Efficiency vs p_{T} Comparison", 720, 720);
    SetCanvasStyle(cEffResults, 0.15, 0.05, 0.05, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cEffResults, pad1Size, pad2Size);
    cEffResults->cd(1);
    SetHistoQA(hEff_1710_toy);
    SetHistoQA(hEff_1525_toy);
    SetHistoQA(hEff_1710_MC);
    SetHistoQA(hEff_1525_MC);
    hEff_1710_toy->SetMarkerSize(1.3);
    hEff_1525_toy->SetMarkerSize(1.3);
    hEff_1710_MC->SetMarkerSize(1.3);
    hEff_1525_MC->SetMarkerSize(1.3);
    hEff_1710_MC->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hEff_1710_MC->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hEff_1710_MC->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hEff_1710_MC->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hEff_1710_MC->GetYaxis()->SetTitleOffset(1.3 * pad1Size);
    hEff_1710_MC->SetLineColor(kBlue);
    hEff_1710_MC->SetMarkerColor(kBlue);
    hEff_1710_MC->GetYaxis()->SetMaxDigits(3);
    hEff_1710_MC->SetMaximum(0.06);
    hEff_1710_MC->SetMinimum(-0.002);
    hEff_1710_MC->Draw();
    hEff_1525_MC->SetLineColor(kRed);
    hEff_1525_MC->SetMarkerColor(kRed);
    hEff_1525_MC->Draw("SAME");
    hEff_1710_toy->SetLineColor(kGreen + 2);
    hEff_1710_toy->SetMarkerColor(kGreen + 2);
    hEff_1710_toy->SetMarkerStyle(21);
    hEff_1710_toy->Draw("SAME");
    hEff_1525_toy->SetLineColor(kMagenta);
    hEff_1525_toy->SetMarkerColor(kMagenta);
    hEff_1525_toy->SetMarkerStyle(21);
    hEff_1525_toy->Draw("SAME");

    TLegend *leg2 = new TLegend(0.2, 0.75, 0.90, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.04);
    leg2->SetTextFont(42);
    leg2->SetNColumns(2);
    leg2->AddEntry(hEff_1525_MC, "f_{2}(1525) (MC)", "l");
    leg2->AddEntry(hEff_1710_MC, "f_{0}(1710) (MC)", "l");
    leg2->AddEntry(hEff_1525_toy, "f_{2}(1525) (Toy model)", "l");
    leg2->AddEntry(hEff_1710_toy, "f_{0}(1710) (Toy model)", "l");
    leg2->Draw();

    cEffResults->cd(2);
    TH1F *hEffRatio1710 = (TH1F *)hEff_1710_MC->Clone("hEffRatio1710");
    TH1F *hEffRatio1525 = (TH1F *)hEff_1525_MC->Clone("hEffRatio1525");
    for (int ibin = 1; ibin <= hEffRatio1710->GetNbinsX(); ibin++)
    {
        if (hEff_1710_toy->GetBinContent(ibin) != 0)
        {
            double binToy = hEff_1710_toy->GetBinContent(ibin);
            double binMC = hEff_1710_MC->GetBinContent(ibin);
            double ratioMCbyToy = binMC / binToy;
            hEffRatio1710->SetBinContent(ibin, ratioMCbyToy);
            hEffRatio1710->SetBinError(ibin, 0);
            cout << "pT bins " << ibin << " , Toy Eff " << binToy << " , MC Eff " << binMC << " , Ratio " << ratioMCbyToy << endl;
        }
        if (hEff_1525_toy->GetBinContent(ibin) == 0)
            continue;
        double binToy1525 = hEff_1525_toy->GetBinContent(ibin);
        double binMC1525 = hEff_1525_MC->GetBinContent(ibin);
        double ratioMCbyToy1525 = binMC1525 / binToy1525;
        // double ratioMCbyToy1525 = binToy1525 / binMC1525;
        hEffRatio1525->SetBinContent(ibin, ratioMCbyToy1525);
        hEffRatio1525->SetBinError(ibin, 0);
    }
    SetHistoQA(hEffRatio1710);
    hEffRatio1710->GetYaxis()->SetTitle("MC/Toy model");
    hEffRatio1710->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    hEffRatio1710->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hEffRatio1710->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hEffRatio1710->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hEffRatio1710->GetYaxis()->SetTitleOffset(1.6 * pad2Size);
    hEffRatio1710->SetLineColor(kBlue);
    hEffRatio1710->SetMarkerColor(kBlue);
    hEffRatio1710->GetYaxis()->SetNdivisions(505);
    hEffRatio1710->SetMaximum(0.53);
    hEffRatio1710->SetMinimum(0.31);
    hEffRatio1710->Draw("HIST");
    SetHistoQA(hEffRatio1525);
    hEffRatio1525->SetLineColor(kRed);
    hEffRatio1525->SetMarkerColor(kRed);
    hEffRatio1525->Draw("HIST SAME");
    TLine *lineRatio = new TLine(0, 1, 15, 1);
    lineRatio->SetLineStyle(1);
    lineRatio->SetLineColor(kBlack);
    lineRatio->Draw();
    cEffResults->SaveAs(path + "pTSmearing/f0_f2_Average_eff_vs_pT.png");
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.02);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad1->SetLeftMargin(0.12);
    pad2->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.13);
}