
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>

void compareFitParameters()
{
    gStyle->SetOptStat(0);

    TFile *fFixed = TFile::Open("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/yield_0_100_default.root", "READ");
    TFile *fFree = TFile::Open("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/yield_0_100.root", "READ");

    if (!fFixed || fFixed->IsZombie() || !fFree || fFree->IsZombie())
    {
        std::cerr << "Error: Cannot open input ROOT files!" << std::endl;
        return;
    }

    TDirectory *dirFixed = (TDirectory *)fFixed->Get("mult_0-100");
    TDirectory *dirFree = (TDirectory *)fFree->Get("mult_0-100");

    if (!dirFixed || !dirFree)
    {
        std::cerr << "Error: mult_0-100 directory not found!" << std::endl;
        return;
    }

    const char *histNames[3] = {"width", "sigma", "yield_integral"};
    const char *yTitles[3] = {"Width", "Detector resolution (#sigma)", "Raw yield"};

    TH1 *hFixed[3] = {nullptr, nullptr, nullptr};
    TH1 *hFree[3] = {nullptr, nullptr, nullptr};

    for (int i = 0; i < 3; i++)
    {

        TH1 *h1 = (TH1 *)dirFixed->Get(histNames[i]);
        TH1 *h2 = (TH1 *)dirFree->Get(histNames[i]);

        if (!h1 || !h2)
        {
            std::cerr << "Error: Histogram " << histNames[i] << " not found!" << std::endl;
            return;
        }

        hFixed[i] = (TH1 *)h1->Clone(Form("%s_fixed", histNames[i]));
        hFree[i] = (TH1 *)h2->Clone(Form("%s_free", histNames[i]));

        hFixed[i]->SetDirectory(nullptr);
        hFree[i]->SetDirectory(nullptr);

        hFixed[i]->SetLineColor(kRed + 1);
        hFixed[i]->SetMarkerColor(kRed + 1);
        hFixed[i]->SetMarkerStyle(20);
        hFixed[i]->SetLineWidth(2);

        hFree[i]->SetLineColor(kBlue + 1);
        hFree[i]->SetMarkerColor(kBlue + 1);
        hFree[i]->SetMarkerStyle(24);
        hFree[i]->SetLineWidth(2);

        hFixed[i]->GetYaxis()->SetTitle(yTitles[i]);
        hFixed[i]->GetXaxis()->SetTitle("pT (GeV/c)");

        if (i == 1)
        {
            // remove the bin errors from 2nd, 5, 6, and 7th bins
            hFree[i]->SetBinError(2, 0);
            hFree[i]->SetBinError(5, 0);
            hFree[i]->SetBinError(6, 0);
            hFree[i]->SetBinError(7, 0);
        }
    }

    TCanvas *c = new TCanvas("cCompare", "Width fixed vs Width free", 1400, 720);
    c->Divide(3, 1, 0.01, 0.01);

    // First two plots: Width and sigma
    for (int i = 0; i < 2; i++)
    {
        c->cd(i + 1);

        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.05);

        hFixed[i]->SetTitle("");
        if (i == 0)
        {
            hFixed[i]->GetYaxis()->SetRangeUser(28e-3, 95e-3);
        }
        else
        {
            hFixed[i]->GetYaxis()->SetRangeUser(-1e-3, 26e-3);
        }
        hFixed[i]->Draw("E1");
        hFree[i]->Draw("E1 SAME");

        if (i == 0)
        {
            TLegend *leg = new TLegend(0.18, 0.72, 0.88, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.06);
            leg->SetHeader("FT0M: 0-100%");
            leg->AddEntry(hFixed[i], "Width fixed", "lep");
            leg->AddEntry(hFree[i], "Width free", "lep");
            leg->Draw();
        }
    }

    // Third plot: Yield with ratio panel
    c->cd(3);

    TPad *padTop = new TPad("padTop", "Yield", 0.0, 0.30, 1.0, 1.0);
    TPad *padBottom = new TPad("padBottom", "Ratio", 0.0, 0.0, 1.0, 0.30);

    padTop->SetBottomMargin(0.02);
    padTop->SetLeftMargin(0.16);
    padTop->SetRightMargin(0.03);
    padTop->SetTopMargin(0.05);
    padTop->SetLogy(1);
    padTop->SetTicks(1, 1);

    padBottom->SetTopMargin(0.03);
    padBottom->SetBottomMargin(0.30);
    padBottom->SetLeftMargin(0.16);
    padBottom->SetRightMargin(0.03);
    padBottom->SetTicks(1, 1);

    padTop->Draw();
    padBottom->Draw();

    // Upper panel: Yield comparison
    padTop->cd();

    hFixed[2]->SetTitle("");
    hFixed[2]->GetXaxis()->SetLabelSize(0);
    hFixed[2]->GetXaxis()->SetTitleSize(0);

    hFixed[2]->Draw("E1");
    hFree[2]->Draw("E1 SAME");

    // Lower panel: Ratio
    padBottom->cd();

    TH1 *hRatio = (TH1 *)hFree[2]->Clone("hRatio");
    hRatio->SetDirectory(nullptr);

    if (hFree[2]->GetNbinsX() != hFixed[2]->GetNbinsX())
    {
        std::cerr << "Error: Yield histograms have different binning!" << std::endl;
        return;
    }

    hRatio->Divide(hFixed[2]);

    hRatio->SetTitle("");
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetLineWidth(2);

    hRatio->GetYaxis()->SetTitle("Free / Fixed");
    hRatio->GetXaxis()->SetTitle("pT (GeV/c)");

    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetTitleOffset(0.55);

    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleOffset(1.0);

    hRatio->SetMinimum(0.55);
    hRatio->SetMaximum(1.6);

    hRatio->Draw("E1");

    TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kRed + 1);
    line->Draw("SAME");

    c->Update();

    // c->SaveAs("WidthFixed_vs_WidthFree.pdf");
    c->SaveAs("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/WidthFixed_vs_WidthFree.pdf");

    fFixed->Close();
    fFree->Close();
}
