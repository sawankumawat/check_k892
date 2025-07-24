#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

void compare_fitParams()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    string path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/448490/kstarqa/hInvMass"; // 2022 data
    string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/449695/kstarqa/hInvMass"; // 2023 data
    // string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/451003/kstarqa/hInvMass"; // 2024 data
    TString outputPath = path2 + "/compare_fitParams";
    gSystem->mkdir(outputPath, kTRUE);
    string outputtype = "png"; // pdf, eps

    TFile *fspectra1 = new TFile((path1 + "/yield.root").c_str(), "read");
    TFile *fspectra2 = new TFile((path2 + "/yield.root").c_str(), "read");

    if (fspectra1->IsZombie() || fspectra2->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    TH1F *rawSpectra1[numofmultbins + 1];
    TH1F *rawSpectra2[numofmultbins + 1];
    TH1F *chi2byNDF1[numofmultbins + 1];
    TH1F *chi2byNDF2[numofmultbins + 1];
    TH1F *significance1[numofmultbins + 1];
    TH1F *significance2[numofmultbins + 1];
    TH1F *mass1[numofmultbins + 1];
    TH1F *mass2[numofmultbins + 1];

    for (int imult = 0; imult < numofmultbins + 1; imult++)
    {
        double multlow = (imult == 0) ? 0 : mult_classes[imult - 1];
        double multhigh = (imult == 0) ? 100 : mult_classes[imult];

        rawSpectra1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));
        rawSpectra2[imult] = (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));
        chi2byNDF1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/chi2byNDF", multlow, multhigh));
        chi2byNDF2[imult] = (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/chi2byNDF", multlow, multhigh));
        significance1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/significance", multlow, multhigh));
        significance2[imult] = (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/significance", multlow, multhigh));
        mass1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/mass", multlow, multhigh));
        mass2[imult] = (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/mass", multlow, multhigh));

        if (rawSpectra1[imult] == nullptr || rawSpectra2[imult] == nullptr ||
            chi2byNDF1[imult] == nullptr || chi2byNDF2[imult] == nullptr ||
            significance1[imult] == nullptr || significance2[imult] == nullptr ||
            mass1[imult] == nullptr || mass2[imult] == nullptr)
        {
            cout << "Error: Histograms not found for multiplicity bin " << imult << endl;
            return;
        }

        TCanvas *csig = new TCanvas("", "", 720, 720);
        SetCanvasStyle(csig, 0.16, 0.05, 0.055, 0.15);
        SetHistoQA(significance1[imult]);
        SetHistoQA(significance2[imult]);
        significance1[imult]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        significance1[imult]->GetYaxis()->SetTitle("Significance");
        significance1[imult]->SetStats(0);
        significance1[imult]->Draw();
        significance2[imult]->SetLineColor(kBlue);
        significance2[imult]->Draw("same");

        TLegend *leg = new TLegend(0.65, 0.8, 0.9, 0.9);
        leg->SetTextFont(42);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(significance1[imult], "2022 dataset", "l");
        leg->AddEntry(significance2[imult], "2023 dataset", "l");
        leg->Draw();
        csig->SaveAs(outputPath + (Form("/significance_%.0f-%.0f.", multlow, multhigh) + outputtype).c_str());

        // // // chisquare_NDF vs pt
        TCanvas *cChiSquare = new TCanvas("cChiSquare", "Chi Square vs pT", 720, 720);
        SetCanvasStyle(cChiSquare, 0.16, 0.05, 0.055, 0.15);
        SetHistoQA(chi2byNDF1[imult]);
        SetHistoQA(chi2byNDF2[imult]);
        chi2byNDF1[imult]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        chi2byNDF1[imult]->GetYaxis()->SetTitle("#chi^{2}/NDF ");
        chi2byNDF1[imult]->SetMaximum(5.5);
        chi2byNDF1[imult]->SetMinimum(0);
        chi2byNDF1[imult]->SetStats(0);
        chi2byNDF1[imult]->Draw("p");
        chi2byNDF2[imult]->SetLineColor(kBlue);
        chi2byNDF2[imult]->SetMarkerColor(kBlue);
        chi2byNDF2[imult]->SetMarkerStyle(22);
        chi2byNDF2[imult]->Draw("p same");
        leg->Clear();
        leg->AddEntry(chi2byNDF1[imult], "2022 dataset", "p");
        leg->AddEntry(chi2byNDF2[imult], "2023 dataset", "p");
        leg->Draw();
        cChiSquare->SaveAs(outputPath + (Form("/chi_%.0f-%.0f.", multlow, multhigh) + outputtype).c_str());

        // // mass vs pt
        TCanvas *cmass = new TCanvas("cmass", "Mass vs pT", 720, 720);
        SetCanvasStyle(cmass, 0.17, 0.05, 0.055, 0.15);
        SetHistoQA(mass1[imult]);
        SetHistoQA(mass2[imult]);
        mass1[imult]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        mass1[imult]->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
        mass1[imult]->GetYaxis()->SetRangeUser(0.878, 0.909);
        mass1[imult]->GetYaxis()->SetTitleOffset(1.8);
        mass1[imult]->SetStats(0);
        mass1[imult]->Draw("pe");
        mass2[imult]->SetLineColor(kBlue);
        mass2[imult]->SetMarkerColor(kBlue);
        mass2[imult]->SetMarkerStyle(22);
        mass2[imult]->Draw("pe same");

        TLine *line = new TLine(mass1[imult]->GetXaxis()->GetXmin(), 0.895, mass1[imult]->GetXaxis()->GetXmax(), 0.895);
        line->SetLineStyle(2);
        line->SetLineColor(2);
        line->SetLineWidth(3);
        line->Draw();
        TLegend *massleg = new TLegend(0.2, 0.75, 0.5, 0.9);
        massleg->SetTextSize(0.035);
        massleg->SetFillStyle(0);
        massleg->SetTextFont(42);
        massleg->AddEntry(mass1[imult], "2022 dataset", "p");
        massleg->AddEntry(mass2[imult], "2023 dataset", "p");
        massleg->AddEntry(line, "PDG Mass", "l");
        massleg->Draw("l");
        cmass->SaveAs(outputPath + (Form("/mass_%.0f-%.0f.", multlow, multhigh) + outputtype).c_str());

        // // // // Yield vs pT (integral method)
        // TCanvas *csigYield = new TCanvas("csigYield", "Yield vs pT (integral method)", 720, 720);
        // SetCanvasStyle(csigYield, 0.16, 0.05, 0.055, 0.15);
        // SetHistoQA(rawSpectra1[imult]);
        // SetHistoQA(rawSpectra2[imult]);
        // rawSpectra1[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        // rawSpectra1[imult]->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
        // gPad->SetLogy(1);
        // rawSpectra1[imult]->SetStats(0);
        // rawSpectra1[imult]->Draw("pe");
        // rawSpectra2[imult]->SetLineColor(kBlue);
        // rawSpectra2[imult]->SetMarkerColor(kBlue);
        // rawSpectra2[imult]->SetMarkerStyle(22);
        // rawSpectra2[imult]->Draw("pe same");
        // leg->Draw();
        // csigYield->SaveAs(outputPath + (Form("/yield_integral_%.0f-%.0f.", multlow, multhigh) + outputtype).c_str());
    }
}