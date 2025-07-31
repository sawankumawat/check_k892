#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
}
void compare_rawCorrecYield()
{
    bool isCorrectedYield = false;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // string path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/448490/kstarqa/hInvMass"; // 2022 data
    // string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/449695/kstarqa/hInvMass"; // 2023 data
    // string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/451003/kstarqa/hInvMass"; // 2024 data

    string path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459845/kstarqa/hInvMass"; // 2022 data
    string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459908/kstarqa/hInvMass"; // 2023 data
    string path3 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/460233/kstarqa/hInvMass"; // 2024 data
    TString outputPath = path3;

    TFile *fspectra1 = (isCorrectedYield) ? new TFile((path1 + "/corrected_spectra.root").c_str(), "read") : new TFile((path1 + "/yield.root").c_str(), "read");
    TFile *fspectra2 = (isCorrectedYield) ? new TFile((path2 + "/corrected_spectra.root").c_str(), "read") : new TFile((path2 + "/yield.root").c_str(), "read");
    TFile *fspectra3 = (isCorrectedYield) ? new TFile((path3 + "/corrected_spectra.root").c_str(), "read") : new TFile((path3 + "/yield.root").c_str(), "read");

    if (fspectra1->IsZombie() || fspectra2->IsZombie() || fspectra3->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    TH1F *hmult1[numofmultbins + 1];
    TH1F *hmult2[numofmultbins + 1];
    TH1F *hmult3[numofmultbins + 1];

    TH1F *hefficiency1[numofmultbins + 1];
    TH1F *hefficiency2[numofmultbins + 1];
    TH1F *hefficiency3[numofmultbins + 1];

    for (int imult = 0; imult < numofmultbins + 1; imult++)
    {
        double multlow = (imult == 0) ? 0 : mult_classes[imult - 1];
        double multhigh = (imult == 0) ? 100 : mult_classes[imult];

        hmult1[imult] = (isCorrectedYield) ? (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", multlow, multhigh)) : (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));
        hmult2[imult] = (isCorrectedYield) ? (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", multlow, multhigh)) : (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));
        hmult3[imult] = (isCorrectedYield) ? (TH1F *)fspectra3->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", multlow, multhigh)) : (TH1F *)fspectra3->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));

        if (isCorrectedYield)
        {
            hefficiency1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/heff", multlow, multhigh));
            hefficiency2[imult] = (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/heff", multlow, multhigh));
            hefficiency3[imult] = (TH1F *)fspectra3->Get(Form("mult_%.0f-%.0f/heff", multlow, multhigh));
        }

        if (hmult1[imult] == nullptr)
        {
            cout << "Histogram others not found" << endl;
            return;
        }

        TH1F *hratio1 = (TH1F *)hmult2[imult]->Clone(Form("ratio_mult_%.0f-%.0f", multlow, multhigh));
        TH1F *hratio2 = (TH1F *)hmult3[imult]->Clone(Form("ratio_mult_%.0f-%.0f", multlow, multhigh));
        hratio1->Divide(hmult1[imult]);
        hratio2->Divide(hmult1[imult]);

        TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
        SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
        double pad1Size, pad2Size;
        canvas_style(c1, pad1Size, pad2Size);
        c1->cd(1);
        gPad->SetLogy(1);
        SetHistoStyle(hmult1[imult], 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
        hmult1[imult]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        hmult1[imult]->SetMaximum(hmult1[imult]->GetMaximum() * 10);
        hmult1[imult]->SetMinimum(hmult1[imult]->GetMinimum() * 1);
        hmult1[imult]->GetYaxis()->SetTitleOffset(1.30);
        hmult1[imult]->GetXaxis()->SetTitleOffset(1.02);
        hmult1[imult]->SetMarkerStyle(20);
        hmult1[imult]->SetMarkerSize(1);
        hmult1[imult]->GetXaxis()->SetRangeUser(0, 10);
        hmult1[imult]->Draw("pe");
        hmult2[imult]->SetMarkerStyle(21);
        hmult2[imult]->SetMarkerSize(1);
        hmult2[imult]->SetMarkerColor(kBlue);
        hmult2[imult]->SetLineColor(kBlue);
        hmult2[imult]->SetLineWidth(2);
        hmult2[imult]->Draw("pe same");
        hmult3[imult]->SetMarkerStyle(22);
        hmult3[imult]->SetMarkerSize(1);
        hmult3[imult]->SetMarkerColor(kRed);
        hmult3[imult]->SetLineColor(kRed);
        hmult3[imult]->SetLineWidth(2);
        hmult3[imult]->Draw("pe same");

        TLegend *leg = new TLegend(0.46, 0.7, 0.9, 0.91);
        SetLegendStyle(leg);
        leg->SetHeader(Form("Multiplicity: %.0f-%.0f%%", multlow, multhigh));
        leg->AddEntry(hmult1[imult], "2022 data", "lpe");
        leg->AddEntry(hmult2[imult], "2023 data", "lpe");
        leg->AddEntry(hmult3[imult], "2024 data", "lpe");
        // leg->AddEntry(hmult1[imult], "Mixed-event", "lpe");
        // leg->AddEntry(hmult2[imult], "Rotated pairs", "lpe");
        // leg->AddEntry(hmult3[imult], "Like-sign pairs", "lpe");
        leg->SetTextSize(0.04);
        leg->SetTextSize(0.05);
        leg->Draw();

        c1->cd(2);
        TH1F *hdummy = (TH1F *)hmult1[0]->Clone();
        for (int i = 0; i < hdummy->GetNbinsX(); i++)
        {
            hdummy->SetBinContent(i + 1, 0);
            hdummy->SetBinError(i + 1, 0);
        }

        SetHistoQA(hratio1);
        hratio1->GetYaxis()->SetTitleSize(0.035 / pad2Size);
        hratio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        hratio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        hratio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        hratio1->SetMarkerStyle(20);
        hratio1->SetMarkerSize(1.0);
        hratio1->SetMarkerColor(kBlue);
        hratio1->SetLineColor(kBlue);
        // hratio1->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
        // hratio1->GetYaxis()->SetTitle("#frac{2023 data}{2022 data}");
        hratio1->GetYaxis()->SetTitle("Ratio to 2022");
        hratio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hratio1->GetXaxis()->CenterTitle(1);
        hratio1->GetYaxis()->SetTitleOffset(0.6);
        hratio1->GetXaxis()->SetTitleOffset(1.1);
        hratio1->GetYaxis()->SetNdivisions(505);
        hratio1->SetMaximum(hratio2->GetMaximum() * 1.1);
        hratio1->GetXaxis()->SetRangeUser(0, 10);
        hratio1->Draw("p");
        hratio2->SetMarkerStyle(21);
        hratio2->SetMarkerSize(1.0);
        hratio2->SetMarkerColor(kRed);
        hratio2->SetLineColor(kRed);
        hratio2->Draw("p same");

        TLine *line = new TLine(0, 1, 10, 1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(1);
        line->Draw();
        if (isCorrectedYield)
        {
            c1->SaveAs(outputPath + Form("/CorrectedYield_%.0f-%.0f.png", multlow, multhigh));
        }
        else
        {
            c1->SaveAs(outputPath + Form("/RawYield_%.0f-%.0f.png", multlow, multhigh));
        }

        if (isCorrectedYield)
        {
            TCanvas *c2 = new TCanvas("c2", "c2", 720, 720);
            SetCanvasStyle(c2, 0.15, 0.03, 0.03, 0.15);
            SetHistoQA(hefficiency1[imult]);
            SetHistoQA(hefficiency2[imult]);
            hefficiency1[imult]->SetMarkerStyle(20);
            hefficiency1[imult]->GetYaxis()->SetTitle("Acceptance x Efficiency");
            hefficiency1[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hefficiency1[imult]->Draw("pe");
            hefficiency2[imult]->SetMarkerStyle(21);
            hefficiency2[imult]->SetMarkerColor(kBlue);
            hefficiency2[imult]->SetLineColor(kBlue);
            hefficiency2[imult]->Draw("pe same");
            TLegend *leg2 = new TLegend(0.46, 0.4, 0.9, 0.61);
            SetLegendStyle(leg2);
            leg2->SetHeader(Form("Multiplicity: %.0f-%.0f%%", multlow, multhigh));
            leg2->AddEntry(hefficiency1[imult], "2022 dataset", "lpe");
            leg2->AddEntry(hefficiency2[imult], "2023 dataset", "lpe");
            leg2->SetTextSize(0.04);
            leg2->Draw();
            c2->SaveAs(outputPath + Form("/EfficiencyMult_%.0f-%.0f.png", multlow, multhigh));
        }
    }
}