#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"

void compare_mass_yield_residual()
{
    gStyle->SetOptStat(0);
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/";
    TFile *fDefault = new TFile((path + "FitParam_Default2_resparams.root").c_str(), "READ");
    if (fDefault->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hMass1525 = (TH1F *)fDefault->Get("Mult_0_100/hMass_1525");
    TH1F *hMass1710 = (TH1F *)fDefault->Get("Mult_0_100/hMass_1710");
    TH1F *hMass1525Res = (TH1F *)fDefault->Get("Mult_0_100/hMass_res_1525");
    TH1F *hMass1710Res = (TH1F *)fDefault->Get("Mult_0_100/hMass_res_1710");

    TH1F *hRawYield1525 = (TH1F *)fDefault->Get("Mult_0_100/hYield_1525");
    TH1F *hRawYield1710 = (TH1F *)fDefault->Get("Mult_0_100/hYield_1710");
    TH1F *hYield1525Res = (TH1F *)fDefault->Get("Mult_0_100/hYield_res_1525");
    TH1F *hYield1710Res = (TH1F *)fDefault->Get("Mult_0_100/hYield_res_1710");
    if (hMass1525 == nullptr)
    {
        cout << "Error: One of the histograms not found!" << endl;
        return;
    }
    TCanvas *cMass1710 = new TCanvas("cMass1710", "Mass Comparison", 720, 720);
    SetCanvasStyle(cMass1710, 0.15, 0.01, 0.05, 0.13);
    SetHistoQA(hMass1710);
    hMass1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass1710->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass1710->GetYaxis()->SetTitleOffset(1.5);
    hMass1710->GetYaxis()->SetRangeUser(1.56, 1.89);
    hMass1710->SetLineColor(kBlue);
    hMass1710->SetMarkerColor(kBlue);
    hMass1710->SetMarkerStyle(20);
    hMass1710->SetMarkerSize(1.5);
    hMass1710->Draw("pe");
    SetHistoQA(hMass1710Res);
    hMass1710Res->SetLineColor(kRed);
    hMass1710Res->SetMarkerColor(kRed);
    hMass1710Res->SetMarkerStyle(25);
    hMass1710Res->SetMarkerSize(1.5);
    hMass1710Res->Draw("pe same");
    TLine *line1710Mass = new TLine(1, f1710Mass, 15, f1710Mass);
    line1710Mass->SetLineStyle(2);
    line1710Mass->SetLineColor(2);
    line1710Mass->Draw("same");
    TBox *band1710Mass = new TBox(1, f1710Mass - f1710MassErr, 15, f1710Mass + f1710MassErr);
    band1710Mass->SetFillStyle(3001);
    band1710Mass->SetFillColorAlpha(kRed, 0.2); // shaded
    band1710Mass->SetLineColor(kRed);
    band1710Mass->SetLineWidth(1);
    band1710Mass->Draw("same");
    TLegend *leg1710Mass = new TLegend(0.22, 0.67, 0.9, 0.93);
    leg1710Mass->SetBorderSize(0);
    leg1710Mass->SetFillStyle(0);
    leg1710Mass->SetTextSize(0.03);
    // leg1710Mass->AddEntry((TObject *)0, "ALICE", "");
    leg1710Mass->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg1710Mass->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    leg1710Mass->AddEntry(hMass1710, "f_{0}(1710) Mass (before residual bkg subtraction)", "p");
    leg1710Mass->AddEntry(hMass1710Res, "f_{0}(1710) Mass (after residual bkg subtraction)", "p");
    band1710Mass->SetLineWidth(0);
    leg1710Mass->AddEntry(line1710Mass, "PDG value", "l");
    leg1710Mass->Draw();
    cMass1710->SaveAs((path + "mult_0-100/Spectra/Massf0_compare_res.png").c_str());

    TCanvas *cMass1525 = new TCanvas("cMass1525", "Mass Comparison", 720, 720);
    SetCanvasStyle(cMass1525, 0.15, 0.03, 0.05, 0.13);
    SetHistoQA(hMass1525);
    hMass1525->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass1525->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass1525->GetYaxis()->SetTitleOffset(1.5);
    hMass1525->GetYaxis()->SetRangeUser(1.485, 1.575);
    hMass1525->SetMarkerStyle(20);
    hMass1525->SetMarkerSize(1.5);
    hMass1525->SetLineColor(kBlue);
    hMass1525->SetMarkerColor(kBlue);
    hMass1525->Draw("pe");
    SetHistoQA(hMass1525Res);
    hMass1525Res->SetLineColor(kRed);
    hMass1525Res->SetMarkerColor(kRed);
    hMass1525Res->SetMarkerStyle(25);
    hMass1525Res->SetMarkerSize(1.5);
    hMass1525Res->Draw("pe same");
        TBox *band1525Mass = new TBox(1, f1525Mass - f1525MassErr, 15, f1525Mass + f1525MassErr);
    band1525Mass->SetFillStyle(3001);
    band1525Mass->SetFillColorAlpha(kRed, 0.2); // shaded
    band1525Mass->SetLineColor(kRed);
    band1525Mass->SetLineWidth(1);
    band1525Mass->Draw("same");
    // leg1710Mass->Draw();
    TLine *line1525Mass = new TLine(1, f1525Mass, 15, f1525Mass);
    line1525Mass->SetLineStyle(2);
    line1525Mass->SetLineColor(2);
    line1525Mass->Draw();
    leg1710Mass->Draw();
        cMass1525->SaveAs((path + "mult_0-100/Spectra/Massf2_compare_res.png").c_str());

    TCanvas *cYield1710 = new TCanvas("cYield1710", "Yield Comparison", 720, 720);
    SetCanvasStyle(cYield1710, 0.15, 0.01, 0.05, 0.13);
    SetHistoQA(hRawYield1710);
    gPad->SetLogy();
    hRawYield1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRawYield1710->GetYaxis()->SetTitle("Raw Yield");
    hRawYield1710->GetYaxis()->SetTitleOffset(1.5);
    hRawYield1710->GetYaxis()->SetRangeUser(3e-9, 2e-4);
    hRawYield1710->SetLineColor(kBlue);
    hRawYield1710->SetMarkerColor(kBlue);
    hRawYield1710->SetMarkerStyle(20);
    hRawYield1710->SetMarkerSize(1.5);
    hRawYield1710->Draw("pe");
    SetHistoQA(hYield1710Res);
    hYield1710Res->SetLineColor(kRed);
    hYield1710Res->SetMarkerColor(kRed);
    hYield1710Res->SetMarkerStyle(25);
    hYield1710Res->SetMarkerSize(1.5);
    hYield1710Res->Draw("pe same");
    leg1710Mass->Draw();
    cYield1710->SaveAs((path + "mult_0-100/Spectra/Yieldf0_compare_res.png").c_str());

    TCanvas *cYield1525 = new TCanvas("cYield1525", "Yield Comparison", 720, 720);
    SetCanvasStyle(cYield1525, 0.15, 0.03, 0.05, 0.13);
    SetHistoQA(hRawYield1525);
    gPad->SetLogy();
    hRawYield1525->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRawYield1525->GetYaxis()->SetTitle("Raw Yield");
    hRawYield1525->GetYaxis()->SetTitleOffset(1.5);
    hRawYield1525->GetYaxis()->SetRangeUser(3e-9, 9e-4);
    hRawYield1525->SetLineColor(kBlue);
    hRawYield1525->SetMarkerColor(kBlue);
    hRawYield1525->SetMarkerStyle(20);
    hRawYield1525->SetMarkerSize(1.5);
    hRawYield1525->Draw("pe");
    SetHistoQA(hYield1525Res);
    hYield1525Res->SetLineColor(kRed);
    hYield1525Res->SetMarkerColor(kRed);
    hYield1525Res->SetMarkerStyle(25);
    hYield1525Res->SetMarkerSize(1.5);
    hYield1525Res->Draw("pe same");
    leg1710Mass->Draw();
    cYield1525->SaveAs((path + "mult_0-100/Spectra/Yieldf2_compare_res.png").c_str());
}