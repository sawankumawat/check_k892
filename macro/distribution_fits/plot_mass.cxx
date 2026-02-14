#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"

void plot_mass()
{
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra";
    TFile *f = new TFile((path + "/spectra_Default2.root").c_str(), "READ");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
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
    hMass1710->SetMarkerSize(1.5);
    hMass1710->Draw("pe");
    hMass1710Sys->SetLineColor(kBlue);
    hMass1710Sys->SetMarkerColor(kBlue);
    hMass1710Sys->SetFillStyle(0);
    hMass1710Sys->Draw("e2 same");
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
    TLegend *leg1710Mass = new TLegend(0.45, 0.67, 0.9, 0.93);
    leg1710Mass->SetBorderSize(0);
    leg1710Mass->SetFillStyle(0);
    leg1710Mass->SetTextSize(0.035);
    leg1710Mass->AddEntry((TObject *)0, "ALICE", "");
    leg1710Mass->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg1710Mass->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    leg1710Mass->AddEntry(hMass1710, "f_{0}(1710)", "p");
    // leg1710Mass->AddEntry(hMass1710, "This analysis", "pe");
    // leg1710Mass->AddEntry(line1710Mass, "PDG value", "l");
    band1710Mass->SetLineWidth(0);
    leg1710Mass->AddEntry(line1710Mass, "PDG value", "l");
    leg1710Mass->Draw();
    cMass1710->SaveAs((path + "/plots/Mass_f0_sys.png").c_str());

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
    hMass1525Sys->Draw("e2 same");
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
    TLegend *leg1525Mass = new TLegend(0.45, 0.67, 0.9, 0.93);
    leg1525Mass->SetBorderSize(0);
    leg1525Mass->SetFillStyle(0);
    leg1525Mass->SetTextSize(0.035);
    leg1525Mass->AddEntry((TObject *)0, "ALICE", "");
    leg1525Mass->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg1525Mass->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    leg1525Mass->AddEntry(hMass1525, "f_{2}'(1525)", "p");
    band1525Mass->SetLineWidth(0);
    leg1525Mass->AddEntry(line1525Mass, "PDG value", "l");
    leg1525Mass->Draw();
    cMass1525->SaveAs((path + "/plots/Mass_f2_sys.png").c_str());
}