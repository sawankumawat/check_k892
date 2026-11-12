#include <iostream>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TString.h>

using namespace std;

TFile *OpenFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (f->IsZombie())
    {
        cout << "Error: File not found: " << path << endl;
        return nullptr;
    }
    return f;
}

void SetupPadStyles(TPad *pad1, TPad *pad2)
{
    pad1->SetPad(0, 0.3, 1, 1);
    pad2->SetPad(0, 0.0, 1, 0.3);

    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.001);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);

    pad2->SetTopMargin(0.001);
    pad2->SetBottomMargin(0.32);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}

void PlotSpeciesSpectra(TFile *fEPOS, const TString &speciesName, double brCorrection = 1.0)
{
    // 10 Multiplicity classes (Cent0 to Cent9) corresponding to FT0M multiplicity bins
    const int nMult = 10;
    const char *multLabels[nMult] = {"0-1%", "1-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-100%"};

    TCanvas *cGrid = new TCanvas(Form("cGrid_%s", speciesName.Data()), Form("%s EPOS Spectra", speciesName.Data()), 2000, 1000);
    cGrid->Divide(4, 3); // 4 columns x 3 rows grid layout

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);

    for (int i = 0; i < nMult; ++i)
    {
        cGrid->cd(i + 1);

        // Define sub-pads for main spectrum (top) and ratio (bottom)
        TPad *padMain = new TPad(Form("padMain_%d", i), "", 0, 0.3, 1, 1);
        TPad *padRatio = new TPad(Form("padRatio_%d", i), "", 0, 0, 1, 0.3);
        SetupPadStyles(padMain, padRatio);

        padMain->Draw();
        padRatio->Draw();

        // ---------------- Main Spectrum Panel ----------------
        padMain->cd();
        padMain->SetLogy(1);

        TString offHistName = Form("%s_Mult_Spectra_OFF/h%s_Pt_OFF_Cent%d", speciesName.Data(), speciesName.Data(), 9-i);
        TString onHistName = Form("%s_Mult_Spectra_ON_reg_res/h%s_Pt_ON_reg_res_Cent%d", speciesName.Data(), speciesName.Data(), 9-i);

        TH1F *hOFF = (TH1F *)fEPOS->Get(offHistName);
        TH1F *hON = (TH1F *)fEPOS->Get(onHistName);

        if (!hOFF || !hON)
        {
            cout << "Warning: Missing histogram for bin " << i << " (" << speciesName << ")" << endl;
            continue;
        }

        // Clone to protect base pointers
        TH1F *hOFF_Plot = (TH1F *)hOFF->Clone(Form("hOFF_Plot_%s_%d", speciesName.Data(), i));
        TH1F *hON_Plot = (TH1F *)hON->Clone(Form("hON_Plot_%s_%d", speciesName.Data(), i));
        // hON_Plot-
        // FIXME: For Kstar Is BR correction missing? For phi the overall spectra with ON is not much lower, so I do not think it is missing

        // Scale by branching ratio if specified (e.g. for Phi UrQMD ON)
        if (brCorrection != 1.0)
        {
            hON_Plot->Scale(1.0 / brCorrection);
        }

        // Formatting Main Spectra
        hOFF_Plot->SetTitle("");
        hOFF_Plot->SetMarkerStyle(20);
        hOFF_Plot->SetMarkerSize(0.5);
        hOFF_Plot->SetMarkerColor(kRed + 1);
        hOFF_Plot->SetLineColor(kRed + 1);

        hON_Plot->SetMarkerStyle(21);
        hON_Plot->SetMarkerSize(0.5);
        hON_Plot->SetMarkerColor(kBlue + 1);
        hON_Plot->SetLineColor(kBlue + 1);

        hOFF_Plot->GetYaxis()->SetTitle("d^{2}N / (d#it{p}_{T} d#it{y}) (GeV/#it{c})^{-1}");
        hOFF_Plot->GetYaxis()->SetTitleSize(0.05 / 0.7);
        hOFF_Plot->GetYaxis()->SetLabelSize(0.045 / 0.7);
        hOFF_Plot->GetYaxis()->SetTitleOffset(1.1 * 0.7);

        if (speciesName == "Phi")
        {
            hOFF_Plot->GetXaxis()->SetRangeUser(0, 6.0);
        }

        hOFF_Plot->Draw("pe");
        hON_Plot->Draw("pe same");

        TLegend *leg = new TLegend(0.45, 0.65, 0.90, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.045 / 0.7);
        leg->SetHeader(Form("%s: %s", speciesName.Data(), multLabels[i]));
        leg->AddEntry(hOFF_Plot, "EPOS UrQMD OFF", "pe");
        leg->AddEntry(hON_Plot, "EPOS UrQMD ON", "pe");
        leg->Draw();

        // ---------------- Ratio Panel ----------------
        padRatio->cd();
        TH1F *hRatio = (TH1F *)hON_Plot->Clone(Form("hRatio_%s_%d", speciesName.Data(), i));
        hRatio->Divide(hOFF_Plot);

        hRatio->SetTitle("");
        hRatio->SetMarkerStyle(20);
        hRatio->SetMarkerSize(0.5);
        hRatio->SetMarkerColor(kBlack);
        hRatio->SetLineColor(kBlack);
        hRatio->SetLineWidth(2);

        hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hRatio->GetXaxis()->SetTitleSize(0.05 / 0.3);
        hRatio->GetXaxis()->SetLabelSize(0.045 / 0.3);
        hRatio->GetXaxis()->SetTitleOffset(0.9);

        hRatio->GetYaxis()->SetTitle("ON / OFF");
        hRatio->GetYaxis()->SetTitleSize(0.045 / 0.3);
        hRatio->GetYaxis()->SetLabelSize(0.04 / 0.3);
        hRatio->GetYaxis()->SetTitleOffset(1.2 * 0.3);
        hRatio->GetYaxis()->SetNdivisions(505);
        hRatio->GetYaxis()->SetRangeUser(0.38, 1.3);

        if (speciesName == "Phi")
        {
            hRatio->GetXaxis()->SetRangeUser(0, 6.0);
        }

        hRatio->Draw("HIST");

        TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
        line->SetLineColor(kGray + 2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    // ---------------- Save Canvas ----------------
    cGrid->Update();

    TString outputDir = "/home/sawan/Storage/check_k892/macro/SpectraRatios/Plots/EPOSComparison";
    gSystem->mkdir(outputDir, kTRUE);

    TString pdfName = Form("%s/%s_SpectraComparison.pdf", outputDir.Data(), speciesName.Data());
    TString pngName = Form("%s/%s_SpectraComparison.png", outputDir.Data(), speciesName.Data());
    TString rootName = Form("%s/%s_SpectraComparison.root", outputDir.Data(), speciesName.Data());

    cGrid->SaveAs(pdfName);
    cGrid->SaveAs(pngName);
    cGrid->SaveAs(rootName);

    cout << "Saved plots for " << speciesName << " in " << outputDir << endl;
}

void SpectraEPOSComparison()
{
    gStyle->SetOptStat(0);

    // Open EPOS file referenced in your code
    TFile *fEPOS = OpenFile("/home/sawan/Storage/EPOS4_Resonance/yield_outfile_light_Bongfix3.root");
    if (!fEPOS)
        return;

    // 1. Plot Kstar for all 10 multiplicity classes (Grid 4x3)
    PlotSpeciesSpectra(fEPOS, "Kstar", 0.667);

    // 2. Plot Phi for all 10 multiplicity classes (Grid 4x3)
    // Applying 0.415 BR correction for UrQMD ON case as done in your code
    double phiBR = 0.415;
    PlotSpeciesSpectra(fEPOS, "Phi", phiBR);
}