#include "../src/style.h"
#include "../spectra/YieldMean.C"

using namespace std;
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void compare_efficiency2()
{
    gStyle->SetOptStat(0);

    TString MCEffpath1 = "/home/sawan/alice/practice/MCPlots/Gap2/";
    TString MCEffpath2 = "/home/sawan/alice/practice/MCPlots/Gap2_Rapidity/";
    // TString MCEffpath1 = "/home/sawan/alice/practice/MCPlots/Gap2_Rapidity/";
    // TString MCEffpath2 = "/home/sawan/check_k892/mc/LHC24l1/";
    // TString MCEffpath3 = "/home/sawan/alice/practice/MCPlots/Gap2_Rapidity/";

    TFile *fEffDefault = new TFile(MCEffpath1 + "efficiency.root", "read");
    TFile *fEffVar1 = new TFile(MCEffpath2 + "efficiency.root", "read");
    // TFile *fEffVar2 = new TFile(MCEffpath3 + "efficiency.root", "read");

    if (fEffVar1->IsZombie() || fEffDefault->IsZombie())
    {
        cout << "Error opening pt smearing results file" << endl;
        return;
    }

    TH1F *hEff_f0_default = (TH1F *)fEffDefault->Get("hefficiencyf0");
    TH1F *hEff_f2_default = (TH1F *)fEffDefault->Get("hefficiencyf2");
    TH1F *hEff_f0_var1 = (TH1F *)fEffVar1->Get("hefficiencyf0");
    TH1F *hEff_f2_var1 = (TH1F *)fEffVar1->Get("hefficiencyf2");
    // TH1F *hEff_f0_var2 = (TH1F *)fEffVar2->Get("hefficiencyf0");
    // TH1F *hEff_f2_var2 = (TH1F *)fEffVar2->Get("hefficiencyf2");
    if (hEff_f0_var1 == nullptr || hEff_f2_var1 == nullptr || hEff_f0_default == nullptr || hEff_f2_default == nullptr)
    {
        cout << "Error: Histograms not found in the file." << endl;
        return;
    }

    TCanvas *cEffResults = new TCanvas("cEffResults", "Efficiency vs p_{T} Comparison", 720, 720);
    SetCanvasStyle(cEffResults, 0.15, 0.05, 0.05, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cEffResults, pad1Size, pad2Size);
    cEffResults->cd(1);
    SetHistoQA(hEff_f0_var1);
    SetHistoQA(hEff_f0_default);
    hEff_f0_default->SetMarkerSize(1.3);
    hEff_f0_default->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hEff_f0_default->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hEff_f0_default->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hEff_f0_default->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hEff_f0_default->GetYaxis()->SetTitleOffset(1.3 * pad1Size);
    hEff_f0_default->SetLineColor(kBlue);
    hEff_f0_default->SetMarkerColor(kBlue);
    hEff_f0_default->GetYaxis()->SetMaxDigits(3);
    hEff_f0_default->SetMaximum(0.07);
    hEff_f0_default->SetMinimum(-0.002);
    hEff_f0_default->Draw();
    hEff_f0_var1->SetMarkerSize(1.3);
    hEff_f2_var1->SetMarkerSize(1.3);
    hEff_f0_var1->SetLineColor(kGreen + 2);
    hEff_f0_var1->SetMarkerColor(kGreen + 2);
    hEff_f0_var1->SetMarkerStyle(21);
    hEff_f0_var1->Draw("SAME");

    TLegend *leg2 = new TLegend(0.15, 0.77, 0.4, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.035);
    leg2->SetTextFont(42);
    leg2->AddEntry((TObject *)0, "f_{0}(1710) AccxEff", "");
    leg2->AddEntry(hEff_f0_default, "Low #it{p}_{T} rise", "lp");
    leg2->AddEntry(hEff_f0_var1, "Flat #it{p}_{T}", "lp");
    // leg2->AddEntry(hEff_f0_default, "2023 (new)", "lp");
    // leg2->AddEntry(hEff_f0_var1, "2022 (old)", "lp");
    leg2->Draw();

    cEffResults->cd(2);
    TH1F *hEffRatio1710 = (TH1F *)hEff_f0_default->Clone("hEffRatio1710");
    TH1F *hEffRatio1525 = (TH1F *)hEff_f2_default->Clone("hEffRatio1525");
    for (int ibin = 1; ibin <= hEffRatio1710->GetNbinsX(); ibin++)
    {
        if (hEff_f0_var1->GetBinContent(ibin) != 0)
        {
            double binToy = hEff_f0_var1->GetBinContent(ibin);
            double binMC = hEff_f0_default->GetBinContent(ibin);
            double ratioMCbyToy = binMC / binToy;
            hEffRatio1710->SetBinContent(ibin, ratioMCbyToy);
            hEffRatio1710->SetBinError(ibin, 0);
            cout << "pT bins " << ibin << " , Toy Eff " << binToy << " , MC Eff " << binMC << " , Ratio " << ratioMCbyToy << endl;
        }
        if (hEff_f2_var1->GetBinContent(ibin) == 0)
            continue;
        double binToy1525 = hEff_f2_var1->GetBinContent(ibin);
        double binMC1525 = hEff_f2_default->GetBinContent(ibin);
        double ratioMCbyToy1525 = binMC1525 / binToy1525;
        // double ratioMCbyToy1525 = binToy1525 / binMC1525;
        hEffRatio1525->SetBinContent(ibin, ratioMCbyToy1525);
        hEffRatio1525->SetBinError(ibin, 0);
    }
    SetHistoQA(hEffRatio1710);
    hEffRatio1710->GetYaxis()->SetTitle("Ratio");
    hEffRatio1710->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    hEffRatio1710->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hEffRatio1710->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hEffRatio1710->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hEffRatio1710->GetYaxis()->SetTitleOffset(1.6 * pad2Size);
    hEffRatio1710->SetLineColor(kBlue);
    hEffRatio1710->SetMarkerColor(kBlue);
    hEffRatio1710->GetYaxis()->SetNdivisions(505);
    hEffRatio1710->SetMaximum(1.75);
    hEffRatio1710->SetMinimum(0.52);
    // hEffRatio1710->SetMaximum(0.53);
    // hEffRatio1710->SetMinimum(0.31);
    hEffRatio1710->Draw("pe");
    TLine *lineRatio = new TLine(1, 1, 15, 1);
    lineRatio->SetLineStyle(2);
    lineRatio->SetLineColor(kBlack);
    lineRatio->Draw();
    cEffResults->SaveAs(MCEffpath1 + "eff_compf0_gap2_woRap.png");
    // cEffResults->SaveAs(MCEffpath1 + "eff_compf0_gap2_4.png");
    // cEffResults->SaveAs(MCEffpath1 + "eff_compf0_2k22vs23.png");

    TCanvas *cEffResultsf2 = new TCanvas("cEffResultsf2", "Efficiency vs p_{T} Comparison for f_{2}(1525)", 720, 720);
    SetCanvasStyle(cEffResultsf2, 0.15, 0.05, 0.05, 0.15);
    canvas_style(cEffResultsf2, pad1Size, pad2Size);
    cEffResultsf2->cd(1);
    SetHistoQA(hEff_f2_var1);
    SetHistoQA(hEff_f2_default);
    hEff_f2_var1->SetMarkerSize(1.3);
    hEff_f2_default->SetMarkerSize(1.3);
    hEff_f2_default->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hEff_f2_default->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hEff_f2_default->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hEff_f2_default->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hEff_f2_default->GetYaxis()->SetTitleOffset(1.3 * pad1Size);
    hEff_f2_default->SetLineColor(kBlue);
    hEff_f2_default->SetMarkerColor(kBlue);
    hEff_f2_default->GetYaxis()->SetMaxDigits(3);
    hEff_f2_default->SetMaximum(0.07);
    hEff_f2_default->SetMinimum(-0.002);
    hEff_f2_default->Draw();
    hEff_f2_var1->SetLineColor(kMagenta);
    hEff_f2_var1->SetMarkerColor(kMagenta);
    hEff_f2_var1->SetMarkerStyle(21);
    hEff_f2_var1->Draw("SAME");

    TLegend *leg3 = new TLegend(0.15, 0.77, 0.4, 0.9);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.035);
    leg3->SetTextFont(42);
    leg3->AddEntry((TObject *)0, "f_{2}'(1525) AccxEff", "");
    leg3->AddEntry(hEff_f2_default, "Low #it{p}_{T} rise", "lp");
    leg3->AddEntry(hEff_f2_var1, "Flat #it{p}_{T}", "lp");
    // leg3->AddEntry(hEff_f2_default, "Gap2", "lp");
    // leg3->AddEntry(hEff_f2_var1, "Gap4", "lp");
    // leg3->AddEntry(hEff_f2_default, "2023 (new)", "lp");
    // leg3->AddEntry(hEff_f2_var1, "2022 (old)", "lp");
    leg3->Draw();

    cEffResultsf2->cd(2);
    SetHistoQA(hEffRatio1525);
    hEffRatio1525->GetYaxis()->SetTitle("Ratio");
    hEffRatio1525->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    hEffRatio1525->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hEffRatio1525->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hEffRatio1525->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hEffRatio1525->GetYaxis()->SetTitleOffset(1.6 * pad2Size);
    hEffRatio1525->SetLineColor(kBlue);
    hEffRatio1525->SetMarkerColor(kBlue);
    hEffRatio1525->GetYaxis()->SetNdivisions(505);
    hEffRatio1525->SetMaximum(1.75);
    hEffRatio1525->SetMinimum(0.52);
    hEffRatio1525->Draw("pe");
    lineRatio->Draw();
    cEffResultsf2->SaveAs(MCEffpath1 + "eff_compf2_gap2_woRap.png");
    // cEffResultsf2->SaveAs(MCEffpath1 + "eff_compf2_gap2_4.png");
    // cEffResultsf2->SaveAs(MCEffpath1 + "eff_compf2_2k22vs23.png");
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
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}