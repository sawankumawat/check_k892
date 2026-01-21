#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"

void plotReweightedMC()
{
    string savePath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra";
    TFile *fReweightf0 = new TFile("ReweightingFactorsf0.root", "read");
    TFile *fReweightf2 = new TFile("ReweightingFactorsf2.root", "read");
    if (fReweightf0->IsZombie() || fReweightf2->IsZombie())
    {
        cout << "Error opening reweighting files" << endl;
        return;
    }

    TH1F *hGenReweighted = (TH1F *)fReweightf0->Get("Genf17102_proj_1_i3");
    TH1F *hRecReweighted = (TH1F *)fReweightf0->Get("Recf1710_pt2_proj_1_i3");
    TH1F *hYieldReweighted = (TH1F *)fReweightf0->Get("hYield1710Corrected_i3");
    TH1F *hGenReweighted2 = (TH1F *)fReweightf2->Get("Genf17102_proj_1_i3");
    TH1F *hRecReweighted2 = (TH1F *)fReweightf2->Get("Recf1710_pt2_proj_1_i3");
    TH1F *hYieldReweighted2 = (TH1F *)fReweightf2->Get("hYield1525Corrected_i3");

    TH1F *hGenUnweighted = (TH1F *)fReweightf0->Get("Genf17102_proj_1_i0");
    TH1F *hRecUnweighted = (TH1F *)fReweightf0->Get("Recf1710_pt2_proj_1_i0");
    TH1F *hGenUnweighted2 = (TH1F *)fReweightf2->Get("Genf17102_proj_1_i0");
    TH1F *hRecUnweighted2 = (TH1F *)fReweightf2->Get("Recf1710_pt2_proj_1_i0");

    if (hGenReweighted == nullptr || hRecReweighted == nullptr || hYieldReweighted == nullptr)
    {
        cout << "Error reading reweighted histograms from file" << endl;
        return;
    }

    TCanvas *cReweighted = new TCanvas("cReweighted", "Reweighted Efficiency for f0(1710)", 720, 720);
    SetCanvasStyle(cReweighted, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hGenReweighted);
    gPad->SetLogy();
    hGenReweighted->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hGenReweighted->GetYaxis()->SetTitle("1/N_{ev} * d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hGenReweighted->GetYaxis()->SetTitleOffset(1.5);
    hGenReweighted->SetMaximum(hGenReweighted->GetMaximum() * 1.5);
    hGenReweighted->SetMarkerStyle(53);
    hGenReweighted->SetMarkerColor(kGreen);
    hGenReweighted->SetLineColor(kGreen);
    hGenReweighted->GetXaxis()->SetRangeUser(0.0, 15.5);
    hGenReweighted->SetMinimum(5e-9);
    hGenReweighted->SetMaximum(hGenReweighted->GetMaximum() * 500);
    hGenReweighted->SetMarkerSize(1.2);
    hGenReweighted->Draw("pe");
    SetHistoQA(hRecReweighted);
    hRecReweighted->SetMarkerStyle(53);
    hRecReweighted->SetMarkerColor(kRed);
    hRecReweighted->SetLineColor(kRed);
    hRecReweighted->SetMarkerSize(1.2);
    hRecReweighted->Draw("pe same");
    SetHistoQA(hYieldReweighted);
    hYieldReweighted->SetMarkerStyle(20);
    hYieldReweighted->SetMarkerColor(kBlack);
    hYieldReweighted->SetLineColor(kBlack);
    hYieldReweighted->SetMarkerSize(1.2);
    hYieldReweighted->Draw("pe same");
    double integralFactor = 20861874; // this is multiplicity in MC
    SetHistoQA(hGenUnweighted);
    hGenUnweighted->Scale(1.0 / integralFactor);
    hGenUnweighted->SetMarkerStyle(53);
    hGenUnweighted->SetMarkerColor(kBlue);
    hGenUnweighted->SetLineColor(kBlue);
    hGenUnweighted->SetMarkerSize(1.2);
    hGenUnweighted->Draw("pe same");
    SetHistoQA(hRecUnweighted);
    hRecUnweighted->Scale(1.0 / integralFactor);
    hRecUnweighted->SetMarkerStyle(53);
    hRecUnweighted->SetMarkerColor(kMagenta);
    hRecUnweighted->SetLineColor(kMagenta);
    hRecUnweighted->SetMarkerSize(1.2);
    hRecUnweighted->Draw("pe same");
    TLegend *legReweighted = new TLegend(0.35, 0.73, 0.9, 0.93);
    legReweighted->SetBorderSize(0);
    legReweighted->SetFillStyle(0);
    legReweighted->SetTextSize(0.03);
    legReweighted->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    legReweighted->AddEntry(hGenUnweighted, "Unweighted Generated distribution", "pe");
    legReweighted->AddEntry(hRecUnweighted, "Unweighted Reconstructed distribution", "pe");
    legReweighted->AddEntry(hYieldReweighted, "Data (f_{0}(1710))", "pe");
    legReweighted->AddEntry(hGenReweighted, "Final Generated distribution", "pe");
    legReweighted->AddEntry(hRecReweighted, "Final Reconstructed distribution", "pe");
    legReweighted->Draw();
    cReweighted->SaveAs((savePath + "/ReweightedEfficiencyf0.pdf").c_str());

    TCanvas *cReweightedf2 = new TCanvas("cReweightedf2", "Reweighted Efficiency for f2'(1525)", 720, 720);
    SetCanvasStyle(cReweightedf2, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hGenReweighted2);
    gPad->SetLogy();
    hGenReweighted2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hGenReweighted2->GetYaxis()->SetTitle("1/N_{ev} * d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hGenReweighted2->GetYaxis()->SetTitleOffset(1.5);
    hGenReweighted2->SetMaximum(hGenReweighted2->GetMaximum() * 1.5);
    hGenReweighted2->SetMarkerStyle(53);
    hGenReweighted2->SetMarkerColor(kGreen);
    hGenReweighted2->SetLineColor(kGreen);
    hGenReweighted2->GetXaxis()->SetRangeUser(0.0, 15.5);
    hGenReweighted2->SetMinimum(1e-9);
    hGenReweighted2->SetMaximum(hGenReweighted2->GetMaximum() * 500);
    hGenReweighted2->SetMarkerSize(1.2);
    hGenReweighted2->Draw("pe");
    SetHistoQA(hRecReweighted2);
    hRecReweighted2->SetMarkerStyle(53);
    hRecReweighted2->SetMarkerColor(kRed);
    hRecReweighted2->SetLineColor(kRed);
    hRecReweighted2->SetMarkerSize(1.2);
    hRecReweighted2->Draw("pe same");
    SetHistoQA(hYieldReweighted2);
    hYieldReweighted2->SetMarkerStyle(20);
    hYieldReweighted2->SetMarkerColor(kBlack);
    hYieldReweighted2->SetLineColor(kBlack);
    hYieldReweighted2->SetMarkerSize(1.2);
    hYieldReweighted2->Draw("pe same");
    SetHistoQA(hGenUnweighted2);
    hGenUnweighted2->Scale(1.0 / integralFactor);
    hGenUnweighted2->SetMarkerStyle(53);
    hGenUnweighted2->SetMarkerColor(kBlue);
    hGenUnweighted2->SetLineColor(kBlue);
    hGenUnweighted2->SetMarkerSize(1.2);
    hGenUnweighted2->Draw("pe same");
    SetHistoQA(hRecUnweighted2);
    hRecUnweighted2->Scale(1.0 / integralFactor);
    hRecUnweighted2->SetMarkerStyle(53);
    hRecUnweighted2->SetMarkerColor(kMagenta);
    hRecUnweighted2->SetLineColor(kMagenta);
    hRecUnweighted2->SetMarkerSize(1.2);
    hRecUnweighted2->Draw("pe same");
    TLegend *legReweighted2 = new TLegend(0.35, 0.73, 0.9, 0.93);
    legReweighted2->SetBorderSize(0);
    legReweighted2->SetFillStyle(0);
    legReweighted2->SetTextSize(0.03);
    legReweighted2->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    legReweighted2->AddEntry(hGenUnweighted2, "Unweighted Generated distribution", "pe");
    legReweighted2->AddEntry(hRecUnweighted2, "Unweighted Reconstructed distribution", "pe");
    legReweighted2->AddEntry(hYieldReweighted2, "Data (f_{2}'(1525))", "pe");
    legReweighted2->AddEntry(hGenReweighted2, "Final Generated distribution", "pe");
    legReweighted2->AddEntry(hRecReweighted2, "Final Reconstructed distribution", "pe");
    legReweighted2->Draw();
    cReweightedf2->SaveAs((savePath + "/ReweightedEfficiencyf2.pdf").c_str());

    // TCanvas *cReweightFactor = new TCanvas("cReweightFactor", "Reweighting Factor for f0(1710)", 720, 720);
    // SetCanvasStyle(cReweightFactor, 0.18, 0.03, 0.05, 0.14);
    // TH1F *hReweightFactor = (TH1F *)hGenReweighted->Clone("hYield1710Corrected_correction_i3");
    // SetHistoQA(hReweightFactor);
    // hReweightFactor->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // hReweightFactor->GetYaxis()->SetTitle("Reweight Factor");
    // hReweightFactor->GetYaxis()->SetTitleOffset(1.5);
    // hReweightFactor->SetMaximum(1.5);
    // hReweightFactor->SetMinimum(0.3);
    // hReweightFactor->Draw("HIST");

    // TFile *file2 = new TFile((savePath + "/spectra.root").c_str(), "read");
    // TH1F *hYield1710Corrected = (TH1F *)file2->Get("hYield1710Corrected");
    // if (hYield1710Corrected == nullptr)
    // {
    //     cout << "Error reading corrected yield histogram from file" << endl;
    //     return;
    // }
    // TCanvas *cYieldCompare = new TCanvas("cYieldCompare", "Corrected Yield Comparison for f0(1710)", 720, 720);
    // SetCanvasStyle(cYieldCompare, 0.18, 0.03, 0.05, 0.14);
    // SetHistoQA(hYield1710Corrected);
    // gPad->SetLogy();
    // hYield1710Corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // hYield1710Corrected->GetYaxis()->SetTitle("1/N_{ev} * d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    // hYield1710Corrected->GetYaxis()->SetTitleOffset(1.5);
    // hYield1710Corrected->SetMaximum(hYield1710Corrected->GetMaximum() * 1.5);
    // hYield1710Corrected->SetMarkerStyle(20);
    // hYield1710Corrected->SetMarkerColor(kRed);
    // hYield1710Corrected->SetLineColor(kRed);
    // hYield1710Corrected->GetXaxis()->SetRangeUser(0.0, 15.5);
    // hYield1710Corrected->SetMinimum(1e-6);
    // hYield1710Corrected->SetMaximum(hYield1710Corrected->GetMaximum() * 12);
    // hYield1710Corrected->SetMarkerSize(1.2);
    // hYield1710Corrected->Draw("pe");
    // hYieldReweighted->Draw("pe same");

    // TLegend *legYieldCompare = new TLegend(0.4, 0.73, 0.9, 0.93);
    // legYieldCompare->SetBorderSize(0);
    // legYieldCompare->SetFillStyle(0);
    // legYieldCompare->SetTextSize(0.035);
    // legYieldCompare->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    // legYieldCompare->AddEntry(hYield1710Corrected, "Corrected Yield (Standard)", "pe");
    // legYieldCompare->AddEntry(hYieldReweighted, "Corrected Yield (Reweighted)", "pe");
    // legYieldCompare->Draw();
    // cYieldCompare->SaveAs((savePath + "/CorrectedYieldComparisonf0.pdf").c_str());
}