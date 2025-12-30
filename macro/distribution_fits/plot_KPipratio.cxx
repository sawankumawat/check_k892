#include <iostream>
#include <utility>
#include "../src/style.h"

using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.38);
    pad1->SetLeftMargin(0.14);
    pad2->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.04);
}

void plot_KPipratio()
{
    // TFile *file = new TFile("../spectra/HEPData-pp5TeV.root", "read");
    // if (file->IsZombie())
    // {
    //     cout << "Error opening file" << endl;
    //     return;
    // }
    // TGraphErrors *gKtoPiRatio = (TGraphErrors *)file->Get("Table 8/Graph1D_y1");
    // TGraphErrors *gptoPiRatio = (TGraphErrors *)file->Get("Table 10/Graph1D_y1");
    // if (gKtoPiRatio == nullptr || gptoPiRatio == nullptr)
    // {
    //     cout << "Error reading graph" << endl;
    //     return;
    // }

    // TCanvas *cRatio1 = new TCanvas("cRatio1", "K/pi ratio", 1080, 720);
    // SetCanvasStyle(cRatio1, 0.16, 0.03, 0.03, 0.16);
    // cRatio1->Divide(2, 1, 0, 0);
    // cRatio1->cd(1);
    // gPad->SetLeftMargin(0.16);
    // SetGrapherrorStyle(gKtoPiRatio);
    // gKtoPiRatio->SetMarkerStyle(22);
    // gKtoPiRatio->SetMarkerSize(1);
    // gKtoPiRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // gKtoPiRatio->GetYaxis()->SetTitle("Ratio");
    // gKtoPiRatio->Draw("AP");
    // TLatex lat;
    // lat.SetTextFont(42);
    // lat.SetTextSize(0.04);
    // lat.SetNDC();
    // lat.DrawLatex(0.2, 0.9, "pp #sqrt{#it{s}} = 5.02 TeV");
    // lat.DrawLatex(0.2, 0.85, "V0M Multiplicity: 0-100%");
    // lat.DrawLatex(0.2, 0.8, "K^{+} + K^{-} / #pi^{+} + #pi^{-}");
    // cRatio1->cd(2);
    // gPad->SetRightMargin(0.01);
    // SetGrapherrorStyle(gptoPiRatio);
    // gptoPiRatio->SetMarkerStyle(22);
    // gptoPiRatio->SetMarkerSize(1);
    // gptoPiRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // gptoPiRatio->SetMaximum(0.38);
    // gptoPiRatio->Draw("AP");
    // lat.SetTextSize(0.045);
    // lat.DrawLatex(0.2, 0.9, "pp #sqrt{#it{s}} = 5.02 TeV");
    // lat.DrawLatex(0.2, 0.85, "V0M Multiplicity: 0-100%");
    // lat.DrawLatex(0.2, 0.8, "p + #bar{p} / #pi^{+} + #pi^{-}");
    // cRatio1->SaveAs("KtoPi_ptRatio_pp5TeV.png");

    int colors[] = {kBlack, kViolet, kBlue, kRed + 1, kCyan + 2, 28, kAzure + 7, kGreen + 2, kOrange + 1, kPink + 7, kGray + 2};
    int markers[] = {21, 22, 23, 33, 34, 43, 45, 29, 39, 47};

    TFile *file = new TFile("../spectra/HEP_PiKP_pp13TeV.root", "read");
    if (file->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TGraphAsymmErrors *hSpectraPion[10];
    TGraphAsymmErrors *hSpectraKaon[10];
    TGraphAsymmErrors *hSpectraProton[10];
    for (int i = 0; i < 9; i++)
    {
        {
            hSpectraPion[i] = (TGraphAsymmErrors *)file->Get(Form("Table 1/Graph1D_y%d", 1 + i));
            hSpectraKaon[i] = (TGraphAsymmErrors *)file->Get(Form("Table 3/Graph1D_y%d", 1 + i));
            hSpectraProton[i] = (TGraphAsymmErrors *)file->Get(Form("Table 5/Graph1D_y%d", 1 + i));
        }

        if (hSpectraPion[i] == nullptr || hSpectraKaon[i] == nullptr || hSpectraProton[i] == nullptr)
        {
            cout << "Error reading spectra for mult class " << i << endl;
            return;
        }
    }
    hSpectraPion[9] = (TGraphAsymmErrors *)hSpectraPion[0]->Clone("hSpectraPion_minBias");
    TGraphAsymmErrors *hPionRatioToMinBias = (TGraphAsymmErrors *)file->Get("Table 7/Graph1D_y1");
    hSpectraKaon[9] = (TGraphAsymmErrors *)hSpectraKaon[0]->Clone("hSpectraKaon_minBias");
    TGraphAsymmErrors *hKaonRatioToMinBias = (TGraphAsymmErrors *)file->Get("Table 9/Graph1D_y1");
    hSpectraProton[9] = (TGraphAsymmErrors *)hSpectraProton[0]->Clone("hSpectraProton_minBias");
    TGraphAsymmErrors *hProtonRatioToMinBias = (TGraphAsymmErrors *)file->Get("Table 11/Graph1D_y1");

    if (hPionRatioToMinBias == nullptr || hKaonRatioToMinBias == nullptr || hProtonRatioToMinBias == nullptr)
    {
        cout << "Error reading ratio to min bias graphs" << endl;
        return;
    }

    for (int ibin = 0; ibin < hSpectraPion[9]->GetN(); ibin++)
    {
        double x, y, xminbias, yminbiasRatio;
        hSpectraPion[9]->GetPoint(ibin, x, y);
        hPionRatioToMinBias->GetPoint(ibin, xminbias, yminbiasRatio);
        if (x != xminbias)
        {
            cout << "Error: x values do not match for pion at bin " << ibin << endl;
            return;
        }
        double YieldMinBias = y / yminbiasRatio;
        hSpectraPion[9]->SetPoint(ibin, x, YieldMinBias);
        // hSpectraPion[9]->SetPointError(ibin, 0, 0);

        hSpectraKaon[9]->GetPoint(ibin, x, y);
        hKaonRatioToMinBias->GetPoint(ibin, xminbias, yminbiasRatio);
        if (x != xminbias)
        {
            cout << "Error: x values do not match for kaon at bin " << ibin << endl;
            return;
        }
        YieldMinBias = y / yminbiasRatio;
        hSpectraKaon[9]->SetPoint(ibin, x, YieldMinBias);
        // hSpectraKaon[9]->SetPointError(ibin, 0, 0);

        hSpectraProton[9]->GetPoint(ibin, x, y);
        hProtonRatioToMinBias->GetPoint(ibin, xminbias, yminbiasRatio);
        if (x != xminbias)
        {
            cout << "Error: x values do not match for proton at bin " << ibin << endl;
            return;
        }
        YieldMinBias = y / yminbiasRatio;
        hSpectraProton[9]->SetPoint(ibin, x, YieldMinBias);
        // hSpectraProton[9]->SetPointError(ibin, 0, 0);
    }

    TGraphErrors *hKtoPiRatio[10];
    TGraphErrors *hPtoPiRatio[10];
    for (int i = 0; i < 10; i++)
    {
        hKtoPiRatio[i] = new TGraphErrors();
        hPtoPiRatio[i] = new TGraphErrors();

        const int nKaonBins = hSpectraKaon[i]->GetN();
        const int nProtonBins = hSpectraProton[i]->GetN();

        for (int ibin = 0; ibin < nKaonBins; ibin++)
        {
            double xkaon, ykaon;
            hSpectraKaon[i]->GetPoint(ibin, xkaon, ykaon);
            double piony = hSpectraPion[i]->Eval(xkaon);
            hKtoPiRatio[i]->SetPoint(ibin, xkaon, ykaon / piony);
            hKtoPiRatio[i]->SetPointError(ibin, 0, 0);

            if (hKtoPiRatio[i] == nullptr)
            {
                cout << "Error creating K/pi ratio graph for mult class " << i << endl;
                return;
            }
        }

        for (int ibin = 0; ibin < nProtonBins; ibin++)
        {
            double xproton, yproton;
            hSpectraProton[i]->GetPoint(ibin, xproton, yproton);
            double piony = hSpectraPion[i]->Eval(xproton);
            hPtoPiRatio[i]->SetPoint(ibin, xproton, yproton / piony);
            hPtoPiRatio[i]->SetPointError(ibin, 0, 0);
        }
    }
    gStyle->SetPalette(kRainBow);
    TCanvas *cRatio = new TCanvas("cRatio2", "K/pi ratio", 720, 720);
    double pad1Size, pad2Size;
    canvas_style(cRatio, pad1Size, pad2Size);
    cRatio->cd(1);
    gPad->SetLeftMargin(0.16);
    // gPad->SetLogy();
    SetGrapherrorStyle(hKtoPiRatio[0]);
    hKtoPiRatio[0]->SetMarkerStyle(markers[0]);
    hKtoPiRatio[0]->SetMarkerSize(1);
    hKtoPiRatio[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hKtoPiRatio[0]->GetYaxis()->SetTitle("Spectra ratio");
    hKtoPiRatio[0]->SetMaximum(0.9);
    hKtoPiRatio[0]->GetYaxis()->SetTitleOffset(1.3);
    hKtoPiRatio[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hKtoPiRatio[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hKtoPiRatio[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hKtoPiRatio[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hKtoPiRatio[0]->Draw("AP");

    string multClassess[] = {"V0M 0-1%", "V0M 1-5%", "V0M 5-10%", "V0M 10-20%", "V0M 20-30%", "V0M 30-40%", "V0M 40-50%", "V0M 50-70%", "V0M 70-100%", "V0M 0-100%"};
    TLegend *leg = new TLegend(0.5, 0.75, 0.95, 0.95);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetNColumns(3);
    leg->AddEntry(hKtoPiRatio[0], multClassess[0].c_str(), "p");
    for (int i = 1; i < 10; i++)
    {
        SetGrapherrorStyle(hKtoPiRatio[i]);
        hKtoPiRatio[i]->SetMarkerStyle(markers[i]);
        hKtoPiRatio[i]->SetMarkerSize(1);
        hKtoPiRatio[i]->SetLineColor(colors[i]);
        hKtoPiRatio[i]->SetMarkerColor(colors[i]);
        // Older ROOT versions may not support PLC/PMC draw options; use plain markers
        hKtoPiRatio[i]->Draw("P same");
        leg->AddEntry(hKtoPiRatio[i], multClassess[i].c_str(), "p");
    }
    leg->Draw();
    TLatex lat;
    lat.SetTextFont(42);
    lat.SetTextSize(0.05);
    lat.SetNDC();
    lat.DrawLatex(0.2, 0.9, "pp #sqrt{#it{s}} = 13 TeV");
    lat.DrawLatex(0.2, 0.83, "(K^{+} + K^{-}) / (#pi^{+} + #pi^{-})");

    cRatio->cd(2);
    gPad->SetRightMargin(0.01);
    TGraphErrors *hKtoPiRatioToMinBias[10];
    for (int i = 0; i < 10; i++)
    {
        hKtoPiRatioToMinBias[i] = new TGraphErrors();
        // hKtoPiRatioToMinBias[i]->Divide(hKtoPiRatio[10]);
        const int nRatioBins = hKtoPiRatio[i]->GetN();
        for (int ibin = 0; ibin < nRatioBins; ibin++)
        {
            double x, y;
            hKtoPiRatio[i]->GetPoint(ibin, x, y);
            double yMinBias;
            hKtoPiRatio[9]->GetPoint(ibin, x, yMinBias);
            hKtoPiRatioToMinBias[i]->SetPoint(ibin, x, y / yMinBias);
            hKtoPiRatioToMinBias[i]->SetPointError(ibin, 0, 0);
        }
    }
    SetGrapherrorStyle(hKtoPiRatioToMinBias[0]);
    hKtoPiRatioToMinBias[0]->SetMarkerStyle(markers[0]);
    hKtoPiRatioToMinBias[0]->SetMarkerSize(1);
    hKtoPiRatioToMinBias[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hKtoPiRatioToMinBias[0]->GetYaxis()->SetTitle("Ratio to MB");
    hKtoPiRatioToMinBias[0]->SetMaximum(1.2);
    hKtoPiRatioToMinBias[0]->SetMinimum(0.7);
    hKtoPiRatioToMinBias[0]->GetYaxis()->SetTitleOffset(0.5);
    hKtoPiRatioToMinBias[0]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hKtoPiRatioToMinBias[0]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hKtoPiRatioToMinBias[0]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hKtoPiRatioToMinBias[0]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hKtoPiRatioToMinBias[0]->GetYaxis()->SetNdivisions(505);
    hKtoPiRatioToMinBias[0]->Draw("AP");
    for (int i = 1; i < 10; i++)
    {
        SetGrapherrorStyle(hKtoPiRatioToMinBias[i]);
        hKtoPiRatioToMinBias[i]->SetMarkerStyle(markers[i]);
        hKtoPiRatioToMinBias[i]->SetMarkerSize(1);
        hKtoPiRatioToMinBias[i]->SetLineColor(colors[i]);
        hKtoPiRatioToMinBias[i]->SetMarkerColor(colors[i]);
        hKtoPiRatioToMinBias[i]->Draw("P same");
    }
    TLine *line = new TLine(0, 1, 20, 1);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw();
    cRatio->SaveAs("KtoPi_ptRatio_pp13TeV.png");

    TCanvas *cRatio2 = new TCanvas("cRatio3", "p/pi ratio", 720, 720);
    canvas_style(cRatio2, pad1Size, pad2Size);
    cRatio2->cd(1);
    gPad->SetLeftMargin(0.16);
    // gPad->SetLogy();
    SetGrapherrorStyle(hPtoPiRatio[0]);
    hPtoPiRatio[0]->SetMarkerStyle(markers[0]);
    hPtoPiRatio[0]->SetMarkerSize(1);
    hPtoPiRatio[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPtoPiRatio[0]->GetYaxis()->SetTitle("Spectra ratio");
    hPtoPiRatio[0]->SetMaximum(0.4);
    hPtoPiRatio[0]->SetMinimum(0.001);
    hPtoPiRatio[0]->GetYaxis()->SetTitleOffset(1.3);
    hPtoPiRatio[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hPtoPiRatio[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hPtoPiRatio[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hPtoPiRatio[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hPtoPiRatio[0]->Draw("AP");
    leg->Clear();
    leg->AddEntry(hPtoPiRatio[0], multClassess[0].c_str(), "p");
    for (int i = 1; i < 10; i++)
    {
        SetGrapherrorStyle(hPtoPiRatio[i]);
        hPtoPiRatio[i]->SetMarkerStyle(markers[i]);
        hPtoPiRatio[i]->SetMarkerSize(1);
        hPtoPiRatio[i]->SetLineColor(colors[i]);
        hPtoPiRatio[i]->SetMarkerColor(colors[i]);
        // Older ROOT versions may not support PLC/PMC draw options; use plain markers
        hPtoPiRatio[i]->Draw("P same");
        leg->AddEntry(hPtoPiRatio[i], multClassess[i].c_str(), "p");
    }
    leg->Draw();
    lat.SetTextFont(42);
    lat.SetTextSize(0.05);
    lat.SetNDC();
    lat.DrawLatex(0.2, 0.9, "pp #sqrt{#it{s}} = 13 TeV");
    lat.DrawLatex(0.2, 0.83, "(p + #bar{p}) / (#pi^{+} + #pi^{-})");
    cRatio2->cd(2);
    gPad->SetRightMargin(0.01);
    TGraphErrors *hPtoPiRatioToMinBias[10];
    for (int i = 0; i < 10; i++)
    {
        hPtoPiRatioToMinBias[i] = new TGraphErrors();
        // hPtoPiRatioToMinBias[i]->Divide(hPtoPiRatio[10]);
        const int nRatioBins = hPtoPiRatio[i]->GetN();
        for (int ibin = 0; ibin < nRatioBins; ibin++)
        {
            double x, y;
            hPtoPiRatio[i]->GetPoint(ibin, x, y);
            double yMinBias;
            hPtoPiRatio[9]->GetPoint(ibin, x, yMinBias);
            hPtoPiRatioToMinBias[i]->SetPoint(ibin, x, y / yMinBias);
            hPtoPiRatioToMinBias[i]->SetPointError(ibin, 0, 0);
        }
    }
    SetGrapherrorStyle(hPtoPiRatioToMinBias[0]);
    hPtoPiRatioToMinBias[0]->SetMarkerStyle(markers[0]);
    hPtoPiRatioToMinBias[0]->SetMarkerSize(1);
    hPtoPiRatioToMinBias[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPtoPiRatioToMinBias[0]->GetYaxis()->SetTitle("Ratio to MB");
    hPtoPiRatioToMinBias[0]->GetYaxis()->SetTitleOffset(0.5);
    hPtoPiRatioToMinBias[0]->SetMaximum(1.35);
    hPtoPiRatioToMinBias[0]->SetMinimum(0.55);
    hPtoPiRatioToMinBias[0]->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hPtoPiRatioToMinBias[0]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hPtoPiRatioToMinBias[0]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hPtoPiRatioToMinBias[0]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hPtoPiRatioToMinBias[0]->GetYaxis()->SetNdivisions(505);
    hPtoPiRatioToMinBias[0]->Draw("AP");
    for (int i = 1; i < 10; i++)
    {
        SetGrapherrorStyle(hPtoPiRatioToMinBias[i]);
        hPtoPiRatioToMinBias[i]->SetMarkerStyle(markers[i]);
        hPtoPiRatioToMinBias[i]->SetMarkerSize(1);
        hPtoPiRatioToMinBias[i]->SetLineColor(colors[i]);
        hPtoPiRatioToMinBias[i]->SetMarkerColor(colors[i]);
        hPtoPiRatioToMinBias[i]->Draw("P same");
    }
    line->Draw();
    cRatio2->SaveAs("PtoPi_ptRatio_pp13TeV.png");
}