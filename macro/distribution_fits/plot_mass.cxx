#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void plot_mass()
{
    gStyle->SetOptStat(0);
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra";
    TFile *f = new TFile((path + "/spectra_Default2.root").c_str(), "READ");
    TFile *fexpol1 = new TFile((path + "/spectra_FitChKstar.root").c_str(), "READ");
    TFile *fhera = new TFile((path + "/spectra_FitExpoHERA.root").c_str(), "READ");
    // TFile *fVariation = new TFile((path + "/../../FitParamBeforeCombinatorialSubtraction.root").c_str(), "READ");
    TFile *fVariation = new TFile((path + "/../../FitParam_Mix.root").c_str(), "READ");
    TFile *fRawDefault = new TFile((path + "/../../FitParam_Default2.root").c_str(), "READ");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    // pt integrated values of mass of f0(1710)
    double f1710Mass_integrated = 1.7091;
    double f1710Mass_StatErr = 0.0045;
    double f1710Mass_SysErr = 0.0142;
    TH1F *hIntegratedMass = new TH1F("hIntegratedMass", "", 1, 0.003, 15);
    hIntegratedMass->SetBinContent(1, f1710Mass_integrated);
    hIntegratedMass->SetBinError(1, f1710Mass_StatErr);
    TH1F *hIntegratedMassSys = new TH1F("hIntegratedMassSys", "", 1, 0.003, 15);
    hIntegratedMassSys->SetBinContent(1, f1710Mass_integrated);
    hIntegratedMassSys->SetBinError(1, f1710Mass_SysErr);

    TH1F *hMass1710 = (TH1F *)f->Get("hMass_1710");
    TH1F *hMass1525 = (TH1F *)f->Get("hMass_1525");
    if (hMass1710 == nullptr || hMass1525 == nullptr)
    {
        cout << "Mass histograms not found" << endl;
        return;
    }
    TFile *fsys = new TFile((path + "/SystematicPlots/SystematicUncertainties.root").c_str(), "READ");
    TH1F *hMass1710Sys = (TH1F *)fsys->Get("TotalSys_Smooth_Mass1710");
    TH1F *hMass1525Sys = (TH1F *)fsys->Get("TotalSys_Smooth_Mass1525");
    if (hMass1710Sys == nullptr || hMass1525Sys == nullptr)
    {
        cout << "Systematic uncertainty histograms not found" << endl;
        return;
    }

    int binsSys = hMass1710Sys->GetNbinsX();
    int binsMass = hMass1710->GetNbinsX();
    // if (binsSys != binsMass)
    // {
    //     cout << "Error: Number of bins in mass histogram and systematic uncertainty histogram do not match" << endl;
    //     cout<<"Bins in systematics histogram: "<<binsSys<<", Bins in mass histogram: "<<binsMass<<endl;
    //     return;
    // }

    for (int ibin = 1; ibin <= binsSys; ibin++)
    {
        double sys1710 = hMass1710Sys->GetBinContent(ibin);
        double sys1525 = hMass1525Sys->GetBinContent(ibin);
        double mass1710 = hMass1710->GetBinContent(ibin + 1);
        double mass1525 = hMass1525->GetBinContent(ibin + 1);
        hMass1710Sys->SetBinContent(ibin, mass1710);
        hMass1525Sys->SetBinContent(ibin, mass1525);
        hMass1710Sys->SetBinError(ibin, sys1710 * mass1710);
        hMass1525Sys->SetBinError(ibin, sys1525 * mass1525);
        if (ibin == 1)
        {
            hMass1710->SetBinError(ibin + 1, hMass1710->GetBinError(ibin + 1) * 3);
        }
    }

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.025);

    /*
    TCanvas *cMass1710 = new TCanvas("cMass1710", "Mass vs pT", 720, 720);
    SetCanvasStyle(cMass1710, 0.15, 0.01, 0.05, 0.13);
    SetHistoQA(hMass1710);
    hMass1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass1710->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass1710->GetYaxis()->SetTitleOffset(1.5);
    hMass1710->GetYaxis()->SetRangeUser(1.56, 1.89);
    hMass1710->SetLineColor(kBlue);
    hMass1710->SetMarkerColor(kBlue);
    hMass1710->SetMarkerStyle(20);
    hMass1710->SetMarkerSize(1.7);
    hMass1710->Draw("pe");
    hMass1710Sys->SetLineColor(kBlue);
    hMass1710Sys->SetMarkerColor(kBlue);
    hMass1710Sys->SetFillStyle(0);
    // // hMass1710Sys->Draw("e2 same");
    // TH1F *hMass1710Expol1 = (TH1F *)fexpol1->Get("hMass_1710");
    // TH1F *hMass1710HERA = (TH1F *)fhera->Get("hMass_1710");
    // SetHistoQA(hMass1710Expol1);
    // SetHistoQA(hMass1710HERA);
    // hMass1710Expol1->SetLineColor(kGreen + 2);
    // hMass1710Expol1->SetMarkerColor(kGreen + 2);
    // hMass1710Expol1->SetMarkerStyle(21);
    // hMass1710Expol1->SetMarkerSize(1.7);
    // hMass1710Expol1->Draw("p same");
    // hMass1710HERA->SetLineColor(kRed);
    // hMass1710HERA->SetMarkerColor(kRed);
    // hMass1710HERA->SetMarkerStyle(22);
    // hMass1710HERA->SetMarkerSize(1.7);
    // hMass1710HERA->Draw("p same");
    TH1F *hMassBeforeCombSub = (TH1F *)fVariation->Get("Mult_0_100/hMass_1710");
    if (hMassBeforeCombSub == nullptr)
    {
        cout << "Mass histogram before combinatorial subtraction not found" << endl;
        return;
    }
    SetHistoQA(hMassBeforeCombSub);
    hMassBeforeCombSub->SetLineColor(kGreen + 2);
    hMassBeforeCombSub->SetMarkerColor(kGreen + 2);
    hMassBeforeCombSub->SetMarkerStyle(21);
    hMassBeforeCombSub->SetMarkerSize(1.7);
    hMassBeforeCombSub->Draw("p same");
    // TLine *line1710PDGMass = new TLine(1, f1710Mass, 15, f1710Mass);
    // line1710PDGMass->SetLineStyle(2);
    // line1710PDGMass->SetLineColor(2);
    // line1710PDGMass->Draw("same");
    // TBox *band1710Mass = new TBox(1, f1710Mass - f1710MassErr, 15, f1710Mass + f1710MassErr);
    // band1710Mass->SetFillStyle(3001);
    // band1710Mass->SetFillColorAlpha(kRed, 0.2); // shaded
    // band1710Mass->SetLineColor(kRed);
    // band1710Mass->SetLineWidth(1);
    // band1710Mass->Draw("same");
    // TLine *line1710HERA = new TLine(1, 1.708, 15, 1.708);
    // line1710HERA->SetLineStyle(2);
    // line1710HERA->SetLineColor(kCyan);
    // line1710HERA->Draw("same");
    // TBox *band1710HERA = new TBox(1, 1.708 - 0.014, 15, 1.708 + 0.014);
    // band1710HERA->SetFillStyle(3001);
    // band1710HERA->SetFillColorAlpha(kCyan, 0.2); // shaded
    // band1710HERA->SetLineColor(kCyan);
    // band1710HERA->SetLineWidth(1);
    // band1710HERA->Draw("same");
    // TLine *line1710HERA = new TLine(1, 1.692, 15, 1.692);
    // line1710HERA->SetLineStyle(2);
    // line1710HERA->SetLineColor(kCyan);
    // line1710HERA->Draw("same");
    // TBox *band1710HERA = new TBox(1, 1.692 - 0.006, 15, 1.692 + 0.006);
    // band1710HERA->SetFillStyle(3001);
    // band1710HERA->SetFillColorAlpha(kCyan, 0.2); // shaded
    // band1710HERA->SetLineColor(kCyan);
    // band1710HERA->SetLineWidth(1);
    // band1710HERA->Draw("same");
    // cout << "mass of f01710 in pdg is " << f1710Mass << " GeV/c2" << endl;
    // SetHistoQA(hIntegratedMass);
    // hIntegratedMass->SetMarkerStyle(21);
    // hIntegratedMass->SetMarkerColor(kGreen + 2);
    // hIntegratedMass->SetLineColor(kGreen + 2);
    // hIntegratedMass->SetMarkerSize(1.5);
    // hIntegratedMass->Draw("p same");
    // SetHistoQA(hIntegratedMassSys);
    // hIntegratedMassSys->SetLineColor(kGreen + 2);
    // hIntegratedMassSys->SetMarkerColor(kGreen + 2);
    // hIntegratedMassSys->SetFillStyle(0);
    // hIntegratedMassSys->Draw("e2 same");
    TLegend *leg1710Mass = new TLegend(0.45, 0.67, 0.9, 0.93);
    leg1710Mass->SetBorderSize(0);
    leg1710Mass->SetFillStyle(0);
    leg1710Mass->SetTextSize(0.03);
    // leg1710Mass->AddEntry((TObject *)0, "ALICE", "");
    // leg1710Mass->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    // leg1710Mass->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    leg1710Mass->SetHeader("f_{0}(1710) Mass");
    // leg1710Mass->AddEntry(hMass1710, "Modified Boltzmann", "pe");
    // leg1710Mass->AddEntry(hMass1710Expol1, "Expol1", "p");
    // leg1710Mass->AddEntry(hMass1710HERA, "Exponential HERA", "p");

    leg1710Mass->AddEntry(hMassBeforeCombSub, "Before combinatorial subtraction", "p");
    leg1710Mass->AddEntry(hMass1710, "After combinatorial subtraction", "p");

    // leg1710Mass->AddEntry(hIntegratedMass, "f_{0}(1710) (integrated)", "p");
    // leg1710Mass->AddEntry(hMass1710, "This analysis", "pe");
    // leg1710Mass->AddEntry(line1710PDGMass, "PDG value", "l");

    // band1710Mass->SetLineWidth(0);
    // leg1710Mass->AddEntry(line1710PDGMass, "PDG value", "l");
    // // leg1710Mass->AddEntry(line1710HERA, "PRL 101, 112003 (2008)", "l");
    // leg1710Mass->AddEntry(line1710HERA, "pt-integrated", "l");
    leg1710Mass->Draw();
    // lat.DrawLatex(0.5, 0.2, "Uncertainties: stat.(bars), sys.(boxes)");

    cMass1710->SaveAs((path + "/plots/Mass_f0_CompareBeforeCombBkg.png").c_str());
    // cMass1710->SaveAs((path + "/plots/Mass_f0_sys.png").c_str());
    // cMass1710->SaveAs((path + "/plots/Mass_f0_sys_compareIntegrated.png").c_str());
    // cMass1710->SaveAs((path + "/plots/Massf0_DifferentFunctions.png").c_str());

    TCanvas *cMass1525 = new TCanvas("cMass1525", "Mass vs pT", 720, 720);
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
    hMass1525Sys->SetLineColor(kBlue);
    hMass1525Sys->SetMarkerColor(kBlue);
    hMass1525Sys->SetFillStyle(0);
    // hMass1525Sys->Draw("e2 same");

    TH1F *hMass1525BeforeCombSub = (TH1F *)fVariation->Get("Mult_0_100/hMass_1525");
    if (hMass1525BeforeCombSub == nullptr)
    {
        cout << "Mass histogram before combinatorial subtraction not found" << endl;
        return;
    }
    SetHistoQA(hMass1525BeforeCombSub);
    hMass1525BeforeCombSub->SetLineColor(kGreen + 2);
    hMass1525BeforeCombSub->SetMarkerColor(kGreen + 2);
    hMass1525BeforeCombSub->SetMarkerStyle(21);
    hMass1525BeforeCombSub->SetMarkerSize(1.7);
    hMass1525BeforeCombSub->Draw("p same");

    TBox *band1525Mass = new TBox(1, f1525Mass - f1525MassErr, 15, f1525Mass + f1525MassErr);
    band1525Mass->SetFillStyle(3001);
    band1525Mass->SetFillColorAlpha(kRed, 0.2); // shaded
    band1525Mass->SetLineColor(kRed);
    band1525Mass->SetLineWidth(1);
    // band1525Mass->Draw("same");

    TLine *line1525Mass = new TLine(1, f1525Mass, 15, f1525Mass);
    line1525Mass->SetLineStyle(2);
    line1525Mass->SetLineColor(2);
    // line1525Mass->Draw();
    TLegend *leg1525Mass = new TLegend(0.45, 0.67, 0.9, 0.93);
    leg1525Mass->SetBorderSize(0);
    leg1525Mass->SetFillStyle(0);
    leg1525Mass->SetTextSize(0.03);
    // leg1525Mass->AddEntry((TObject *)0, "ALICE", "");
    // leg1525Mass->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    // leg1525Mass->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    // leg1525Mass->AddEntry(hMass1525, "f_{2}'(1525)", "p");
    // band1525Mass->SetLineWidth(0);
    // leg1525Mass->AddEntry(line1525Mass, "PDG value", "l");
    // lat.DrawLatex(0.5, 0.2, "Uncertainties: stat.(bars), sys.(boxes)");

    leg1525Mass->SetHeader("f_{2}'(1525) Mass");
    leg1525Mass->AddEntry(hMass1525BeforeCombSub, "Before combinatorial subtraction", "p");
    leg1525Mass->AddEntry(hMass1525, "After combinatorial subtraction", "p");
    leg1525Mass->Draw();
    cMass1525->SaveAs((path + "/plots/Mass_f2_CompareBeforeCombBkg.png").c_str());
    // cMass1525->SaveAs((path + "/plots/Mass_f2_sys.png").c_str());

    */

    // // Plot width and compare with pt-integrated values
    // TFile *fWidthFree = new TFile((path + "/../../FitParam_f0WidthFree2.root").c_str(), "READ");
    // if (fWidthFree->IsZombie())
    // {
    //     cout << "File not found " << (path + "/../../FitParam_f0WidthFree2.root").c_str() << endl;
    //     return;
    // }
    // TH1F *hWidth1710 = (TH1F *)fWidthFree->Get("Mult_0_100/hWidth_1710");
    // if (hWidth1710 == nullptr)
    // {
    //     cout << "Width histogram not found" << endl;
    //     return;
    // }
    // double f1710Width_integrated = 0.1456;
    // double f1710Width_StatErr = 0.0151;
    // double f1710Width_SysErr = 0.0419;
    // TH1F *hIntegratedWidth = new TH1F("hIntegratedWidth", "", 1, 0.003, 15);
    // hIntegratedWidth->SetBinContent(1, f1710Width_integrated);
    // hIntegratedWidth->SetBinError(1, f1710Width_StatErr);
    // TH1F *hIntegratedWidthSys = new TH1F("hIntegratedWidthSys", "", 1, 0.003, 15);
    // hIntegratedWidthSys->SetBinContent(1, f1710Width_integrated);
    // hIntegratedWidthSys->SetBinError(1, f1710Width_SysErr);

    // TCanvas *cWidth1710 = new TCanvas("cWidth1710", "Width vs pT", 720, 720);
    // SetCanvasStyle(cWidth1710, 0.15, 0.01, 0.05, 0.13);
    // TH1F *histDummy = new TH1F("histDummy", "", hMass1710->GetNbinsX(), hMass1710->GetXaxis()->GetXmin(), hMass1710->GetXaxis()->GetXmax());
    // SetHistoQA(histDummy);
    // histDummy->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // histDummy->GetYaxis()->SetTitle("Width (GeV/#it{c}^{2})");
    // histDummy->GetYaxis()->SetTitleOffset(1.5);
    // histDummy->GetYaxis()->SetRangeUser(0.0, 0.409);
    // histDummy->Draw();
    // SetHistoQA(hWidth1710);
    // // hWidth1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // // hWidth1710->GetYaxis()->SetTitle("Width (GeV/#it{c}^{2})");
    // // hWidth1710->GetYaxis()->SetTitleOffset(1.5);
    // // hWidth1710->GetYaxis()->SetRangeUser(0.0, 0.429);
    // hWidth1710->SetLineColor(kBlue);
    // hWidth1710->SetMarkerColor(kBlue);
    // hWidth1710->SetMarkerStyle(20);
    // hWidth1710->SetMarkerSize(1.5);
    // hWidth1710->GetXaxis()->SetRangeUser(0, 16);
    // hWidth1710->Draw("pe same");
    // TLine *line1710Width = new TLine(1, f1710Width, 15, f1710Width);
    // line1710Width->SetLineStyle(2);
    // line1710Width->SetLineColor(2);
    // line1710Width->Draw("same");
    // TBox *band1710Width = new TBox(1, f1710Width - f1710WidthErr, 15, f1710Width + f1710WidthErr);
    // band1710Width->SetFillStyle(3001);
    // band1710Width->SetFillColorAlpha(kRed, 0.2); // shaded
    // band1710Width->SetLineColor(kRed);
    // band1710Width->SetLineWidth(1);
    // band1710Width->Draw("same");
    // cout << "Width of f01710 in pdg is " << f1710Width << " GeV/c2" << endl;
    // SetHistoQA(hIntegratedWidth);
    // hIntegratedWidth->SetMarkerStyle(21);
    // hIntegratedWidth->SetMarkerColor(kGreen + 2);
    // hIntegratedWidth->SetLineColor(kGreen + 2);
    // hIntegratedWidth->SetMarkerSize(1.5);
    // hIntegratedWidth->Draw("p same");
    // SetHistoQA(hIntegratedWidthSys);
    // hIntegratedWidthSys->SetLineColor(kGreen + 2);
    // hIntegratedWidthSys->SetMarkerColor(kGreen + 2);
    // hIntegratedWidthSys->SetFillStyle(0);
    // hIntegratedWidthSys->Draw("e2 same");
    // TLegend *leg1710Width = new TLegend(0.45, 0.67, 0.9, 0.93);
    // leg1710Width->SetBorderSize(0);
    // leg1710Width->SetFillStyle(0);
    // leg1710Width->SetTextSize(0.035);
    // leg1710Width->AddEntry((TObject *)0, "ALICE", "");
    // leg1710Width->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    // leg1710Width->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    // leg1710Width->AddEntry(hWidth1710, "f_{0}(1710)", "p");
    // leg1710Width->AddEntry(hIntegratedWidth, "f_{0}(1710) (integrated)", "p");
    // // leg1710Width->AddEntry(hWidth1710, "This analysis", "pe");
    // // leg1710Width->AddEntry(line1710Width, "PDG value", "l");
    // band1710Width->SetLineWidth(0);
    // leg1710Width->AddEntry(line1710Width, "PDG value", "l");
    // lat.DrawLatex(0.5, 0.2, "Uncertainties: stat.(bars), sys.(boxes)");
    // leg1710Width->Draw();
    // cWidth1710->SaveAs((path + "/plots/Width_f0_sys_compareIntegrated.png").c_str());

    TH1F *hRawf2Default = (TH1F *)fRawDefault->Get("Mult_0_100/hYield_1525");
    TH1F *hRawf2Variation = (TH1F *)fVariation->Get("Mult_0_100/hYield_1525");
    if (hRawf2Default == nullptr || hRawf2Variation == nullptr)
    {
        cout << "Yield histograms not found" << endl;
        return;
    }
    TCanvas *cYield1525 = new TCanvas("cYield1525", "Yield vs pT", 720, 720);
    SetCanvasStyle(cYield1525, 0.15, 0.03, 0.05, 0.13);
    double pad1Size, pad2Size;
    canvas_style(cYield1525, pad1Size, pad2Size);
    cYield1525->cd(1);
    gPad->SetLogy();
    SetHistoQA(hRawf2Default);
    hRawf2Default->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRawf2Default->GetYaxis()->SetTitle("1/N_{evt} #times d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hRawf2Default->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hRawf2Default->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hRawf2Default->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hRawf2Default->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hRawf2Default->GetYaxis()->SetTitleOffset(1.6 * pad1Size);
    // hRawf2Default->GetYaxis()->SetRangeUser(1e-7, 1e-3);
    hRawf2Default->SetLineColor(kBlue);
    hRawf2Default->SetMarkerColor(kBlue);
    hRawf2Default->SetMarkerStyle(20);
    hRawf2Default->SetMarkerSize(1.5);
    hRawf2Default->Draw("pe");
    hRawf2Variation->SetLineColor(kGreen + 2);
    hRawf2Variation->SetMarkerColor(kGreen + 2);
    hRawf2Variation->SetMarkerStyle(21);
    hRawf2Variation->SetMarkerSize(1.5);
    hRawf2Variation->Draw("p same");
    TLegend *legYield1525 = new TLegend(0.45, 0.67, 0.9, 0.93);
    legYield1525->SetBorderSize(0);
    legYield1525->SetFillStyle(0);
    legYield1525->SetTextSize(0.04);
    legYield1525->AddEntry((TObject *)0, "f_{2}'(1525) Yield", "");
    legYield1525->AddEntry(hRawf2Default, "Rotational", "p");
    legYield1525->AddEntry(hRawf2Variation, "Mixed-event", "p");
    // legYield1525->AddEntry(hRawf2Default, "After combinatorial subtraction", "p");
    // legYield1525->AddEntry(hRawf2Variation, "Before combinatorial subtraction", "p");
    legYield1525->Draw();

    cYield1525->cd(2);
    TH1F *hRatiof2 = (TH1F *)hRawf2Default->Clone("hRatiof2");
    int totalBins = hRawf2Default->GetNbinsX();
    // for (int ibin = 4; ibin <= totalBins; ibin++)
    cout << "f2(1525) Yield Ratio (Rotational / Mixed-event):" << endl;
    for (int ibin = 3; ibin <= totalBins; ibin++)
    {
        double defaultYield = hRawf2Default->GetBinContent(ibin);
        // double variationYield = hRawf2Variation->GetBinContent(ibin - 3);
        double variationYield = hRawf2Variation->GetBinContent(ibin - 2);
        double ratio = (variationYield != 0) ? defaultYield / variationYield : 0;
        cout << "Bin " << ibin << ", Ratio = " << ratio << endl;
        hRatiof2->SetBinContent(ibin, ratio);
        double defaultError = hRawf2Default->GetBinError(ibin);
        // double variationError = hRawf2Variation->GetBinError(ibin - 3);
        double variationError = hRawf2Variation->GetBinError(ibin - 2);
        double ratioError = 0;
        if (defaultYield != 0)
        {
            ratioError = ratio * sqrt(pow(variationError / variationYield, 2) + pow(defaultError / defaultYield, 2));
        }
        hRatiof2->SetBinError(ibin, ratioError);
    }
    SetHistoQA(hRatiof2);
    hRatiof2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // hRatiof2->GetYaxis()->SetTitle("Variation / Default");
    hRatiof2->GetYaxis()->SetTitle("Rotational / ME");
    hRatiof2->GetYaxis()->SetTitleSize(0.03 / pad2Size);
    hRatiof2->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof2->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof2->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatiof2->GetYaxis()->SetTitleOffset(2.2 * pad2Size);
    hRatiof2->SetLineColor(kBlack);
    hRatiof2->SetMarkerColor(kBlack);
    hRatiof2->SetMarkerStyle(20);
    hRatiof2->SetMarkerSize(1.5);
    hRatiof2->GetYaxis()->SetRangeUser(0.75, 1.25);
    hRatiof2->GetYaxis()->SetNdivisions(505);
    hRatiof2->Draw("pe");
    TLine *lineOnef2 = new TLine(hRatiof2->GetXaxis()->GetXmin(), 1, hRatiof2->GetXaxis()->GetXmax(), 1);
    lineOnef2->SetLineStyle(2);
    lineOnef2->SetLineColor(kRed);
    lineOnef2->Draw("same");
    cYield1525->SaveAs((path + "/plots/Yield_f2_CompareMIX.png").c_str());
    // cYield1525->SaveAs((path + "/plots/Yield_f2_CompareBeforeCombBkg.png").c_str());

    TH1F *hRawf0Default = (TH1F *)fRawDefault->Get("Mult_0_100/hYield_1710");
    TH1F *hRawf0Variation = (TH1F *)fVariation->Get("Mult_0_100/hYield_1710");
    if (hRawf0Default == nullptr || hRawf0Variation == nullptr)
    {
        cout << "Yield histograms not found" << endl;
        return;
    }
    TCanvas *cYield1710 = new TCanvas("cYield1710", "Yield vs pT", 720, 720);
    SetCanvasStyle(cYield1710, 0.15, 0.03, 0.05, 0.13);
    canvas_style(cYield1710, pad1Size, pad2Size);
    cYield1710->cd(1);
    gPad->SetLogy();
    SetHistoQA(hRawf0Default);
    hRawf0Default->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRawf0Default->GetYaxis()->SetTitle("1/N_{evt} #times d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hRawf0Default->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hRawf0Default->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hRawf0Default->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hRawf0Default->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hRawf0Default->GetYaxis()->SetTitleOffset(1.5 * pad1Size);
    // hRawf0Default->GetYaxis()->SetRangeUser(1e-7, 1e-3);
    hRawf0Default->SetLineColor(kBlue);
    hRawf0Default->SetMarkerColor(kBlue);
    hRawf0Default->SetMarkerStyle(20);
    hRawf0Default->SetMarkerSize(1.5);
    hRawf0Default->Draw("pe");
    hRawf0Variation->SetLineColor(kGreen + 2);
    hRawf0Variation->SetMarkerColor(kGreen + 2);
    hRawf0Variation->SetMarkerStyle(21);
    hRawf0Variation->SetMarkerSize(1.5);
    hRawf0Variation->Draw("p same");
    TLegend *legYield1710 = new TLegend(0.45, 0.67, 0.9, 0.93);
    legYield1710->SetBorderSize(0);
    legYield1710->SetFillStyle(0);
    legYield1710->SetTextSize(0.04);
    legYield1710->AddEntry((TObject *)0, "f_{0}(1710) Yield", "");
    legYield1710->AddEntry(hRawf0Default, "Rotational", "p");
    legYield1710->AddEntry(hRawf0Variation, "Mixed-event", "p");
    // legYield1710->AddEntry(hRawf0Default, "After combinatorial subtraction", "p");
    // legYield1710->AddEntry(hRawf0Variation, "Before combinatorial subtraction", "p");
    legYield1710->Draw();

    cYield1710->cd(2);
    TH1F *hRatiof0 = (TH1F *)hRawf0Default->Clone("hRatiof0");
    totalBins = hRawf0Default->GetNbinsX();
    cout << "f0(1710) Yield Ratio (Rotational / Mixed-event):" << endl;
    for (int ibin = 3; ibin <= totalBins; ibin++)
    {
        double defaultYield = hRawf0Default->GetBinContent(ibin);
        double variationYield = hRawf0Variation->GetBinContent(ibin - 2);
        double ratio = (variationYield != 0) ? defaultYield / variationYield : 0;
        hRatiof0->SetBinContent(ibin, ratio);
        cout << "Bin " << ibin << ", Ratio = " << ratio << endl;
        double defaultError = hRawf0Default->GetBinError(ibin);
        double variationError = hRawf0Variation->GetBinError(ibin - 2);
        double ratioError = 0;
        if (defaultYield != 0)
        {
            ratioError = ratio * sqrt(pow(variationError / variationYield, 2) + pow(defaultError / defaultYield, 2));
        }
        hRatiof0->SetBinError(ibin, ratioError);
    }

    SetHistoQA(hRatiof0);
    hRatiof0->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // hRatiof0->GetYaxis()->SetTitle("Variation / Default");
    hRatiof0->GetYaxis()->SetTitle("Rotational / ME");
    hRatiof0->GetYaxis()->SetTitleSize(0.03 / pad2Size);
    hRatiof0->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof0->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof0->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatiof0->GetYaxis()->SetTitleOffset(2.2 * pad2Size);
    hRatiof0->SetLineColor(kBlack);
    hRatiof0->SetMarkerColor(kBlack);
    hRatiof0->SetMarkerStyle(20);
    hRatiof0->SetMarkerSize(1.5);
    hRatiof0->GetYaxis()->SetRangeUser(0.75, 1.25);
    hRatiof0->GetYaxis()->SetNdivisions(505);
    hRatiof0->Draw("pe");
    TLine *lineOnef0 = new TLine(hRatiof0->GetXaxis()->GetXmin(), 1, hRatiof0->GetXaxis()->GetXmax(), 1);
    lineOnef0->SetLineStyle(2);
    lineOnef0->SetLineColor(kRed);
    lineOnef0->Draw("same");
    cYield1710->SaveAs((path + "/plots/Yield_f0_CompareMIX.png").c_str());
    // cYield1710->SaveAs((path + "/plots/Yield_f0_CompareBeforeCombBkg.png").c_str());
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
    pad1->SetLeftMargin(0.14);
    pad2->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.13);
}