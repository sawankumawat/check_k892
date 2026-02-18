#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"

void compare_mass_yield_coherent()
{
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra/";
    TFile *fDefault = new TFile((path + "spectra_Default2.root").c_str(), "READ");
    TFile *fCoherent = new TFile((path + "spectra_coherent.root").c_str(), "READ");
    if (fDefault->IsZombie() || fCoherent->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hMass_f1525_Default = (TH1F *)fDefault->Get("hMass_1525");
    TH1F *hMass_f1710_Default = (TH1F *)fDefault->Get("hMass_1710");
    TH1F *hMass_f1525_Coherent = (TH1F *)fCoherent->Get("hMass1525");
    TH1F *hMass_f1710_Coherent = (TH1F *)fCoherent->Get("hMass1710");
    if (hMass_f1525_Default == nullptr || hMass_f1710_Default == nullptr || hMass_f1525_Coherent == nullptr || hMass_f1710_Coherent == nullptr)
    {
        cout << "Error: One of the histograms not found!" << endl;
        return;
    }
    TCanvas *cMassf0 = new TCanvas("cMassf0", "Mass Comparison", 720, 720);
    SetCanvasStyle(cMassf0, 0.15, 0.03, 0.05, 0.14);
    SetHistoQA(hMass_f1710_Default);
    SetHistoQA(hMass_f1710_Coherent);
    hMass_f1710_Default->SetLineColor(kRed);
    hMass_f1710_Default->SetMarkerColor(kRed);
    hMass_f1710_Default->SetMarkerSize(1.5);
    hMass_f1710_Default->SetBinError(2, hMass_f1710_Default->GetBinError(2) * 4);
    hMass_f1710_Coherent->SetBinError(2, hMass_f1710_Coherent->GetBinError(2) * 4);
    hMass_f1710_Coherent->SetLineColor(kBlue);
    hMass_f1710_Coherent->SetMarkerColor(kBlue);
    hMass_f1710_Coherent->SetMarkerSize(1.5);
    hMass_f1710_Default->Draw("pe");
    hMass_f1710_Coherent->Draw("pe same");
    TLine *line1710Mass = new TLine(1, f1710Mass, 15, f1710Mass);
    line1710Mass->SetLineStyle(2);
    line1710Mass->SetLineColor(2);
    line1710Mass->Draw("same");
    TBox *band1710Mass = new TBox(1, f1710Mass - f1710MassErr, 15, f1710Mass + f1710MassErr);
    band1710Mass->SetFillColor(kRed); // shaded
    band1710Mass->SetFillStyle(3001); // hatched
    band1710Mass->SetLineColor(kRed);
    band1710Mass->SetLineWidth(1);
    band1710Mass->Draw("same");
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->SetTextFont(42);
    legend->SetHeader("f_{0}(1710) mass");
    legend->AddEntry(hMass_f1710_Default, "Default fit", "p");
    legend->AddEntry(hMass_f1710_Coherent, "Coherent Sum", "p");
    legend->AddEntry(line1710Mass, "PDG mass", "l");
    legend->Draw();
    cMassf0->SaveAs((path + "MassComparison_f1710.png").c_str());

    TCanvas *cMassf2 = new TCanvas("cMassf2", "Mass Comparison", 720, 720);
    SetCanvasStyle(cMassf2, 0.15, 0.03, 0.05, 0.14);
    SetHistoQA(hMass_f1525_Default);
    SetHistoQA(hMass_f1525_Coherent);
    hMass_f1525_Default->SetLineColor(kRed);
    hMass_f1525_Default->SetMarkerColor(kRed);
    hMass_f1525_Default->SetMarkerSize(1.5);
    hMass_f1525_Coherent->SetLineColor(kBlue);
    hMass_f1525_Coherent->SetMarkerColor(kBlue);
    hMass_f1525_Coherent->SetMarkerSize(1.5);
    hMass_f1525_Default->Draw("pe");
    hMass_f1525_Coherent->Draw("pe same");
    TLine *line1525Mass = new TLine(1, f1525Mass, 15, f1525Mass);
    line1525Mass->SetLineStyle(2);
    line1525Mass->SetLineColor(2);
    line1525Mass->Draw("same");
    TBox *band1525Mass = new TBox(1, f1525Mass - f1525MassErr, 15, f1525Mass + f1525MassErr);
    band1525Mass->SetFillColor(kRed); // shaded
    band1525Mass->SetFillStyle(3001); // hatched
    band1525Mass->SetLineColor(kRed);
    band1525Mass->SetLineWidth(1);
    band1525Mass->Draw("same");
    TLegend *legend2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(0.035);
    legend2->SetTextFont(42);
    legend2->SetHeader("f_{2}'(1525) mass");
    legend2->AddEntry(hMass_f1525_Default, "Default fit", "p");
    legend2->AddEntry(hMass_f1525_Coherent, "Coherent Sum", "p");
    legend2->AddEntry(line1525Mass, "PDG mass", "l");
    legend2->Draw();
    cMassf2->SaveAs((path + "MassComparison_f1525.png").c_str());

    // Now compare the <pT> values for f0(1710) and f2(1525) between the two fits
    TGraphErrors *gMeanpTBothDefault = (TGraphErrors *)fDefault->Get("gMeanPt_f0f2");
    TGraphErrors *gMeanpTBothCoherent = (TGraphErrors *)fCoherent->Get("gMeanPt_f0f2");
    if (gMeanpTBothDefault == nullptr || gMeanpTBothCoherent == nullptr)
    {
        cout << "Error: gMeanPt_f0f2 not found in one of the files!" << endl;
        return;
    }
        double f1525MeanPt_Default = gMeanpTBothDefault->GetY()[0];
        double f1525MeanPt_Default_err = gMeanpTBothDefault->GetErrorY(0);
        double f1710MeanPt_Default = gMeanpTBothDefault->GetY()[1];
        double f1710MeanPt_Default_err = gMeanpTBothDefault->GetErrorY(1);
    double f1525MeanPt_Coherent = gMeanpTBothCoherent->GetY()[0];
    double f1525MeanPt_Coherent_err = gMeanpTBothCoherent->GetErrorY(0);
    double f1710MeanPt_Coherent = gMeanpTBothCoherent->GetY()[1];
    double f1710MeanPt_Coherent_err = gMeanpTBothCoherent->GetErrorY(1);

    TFile *flightFlavourHadrons = new TFile("../spectra/LightFlavourHadronsProduction.root", "read");
    if (flightFlavourHadrons->IsZombie())
    {
        cout << "Error opening light flavour hadrons production file" << endl;
        return;
    }
    // There are graphs of <pT> in the file in table 26 to 34 for different particles (pion, kaon, K0s, K*(892), phi, proton, Lambda, sigma, omega)
    //(π+/π−,K+/K−,KS0​,K∗(892),ϕ(1020),p/pˉ​,Λ/Λˉ,Σ+/Σ−,Ω−/Ωˉ+)

    int totalParticles = 9;
    string particles[9] = {"Pion", "Kaon", "K0s", "Kstar", "Phi", "Proton", "Lambda", "Xi", "Omega"};
    string particlesLatex[9] = {"#pi", "K", "K^{0}_{S}", "K^{*0}", "#phi", "p", "#Lambda", "#Xi^{-}", "#Omega"};
    double particleMass[9] = {0.13957, 0.49367, 0.49761, 0.89166, 1.01946, 0.93827, 1.11568, 1.3217, 1.67245}; // in GeV/c2
    int colors[9] = {kBlack, kBlue, kGreen + 2, kOrange + 7, kViolet + 7, kCyan + 2, kMagenta + 2, kGray + 2, kPink + 2};
    int markers[9] = {21, 22, 23, 33, 34, 43, 45, 29, 39};
    vector<vector<float>> dNdyvalues_13TeV = {
        {4.775, 0.001, 0.243, 1.0},  // 4.775 ± 0.001 ± 0.243
        {6.205, 0.004, 0.303, 1e-1}, // (6.205 ± 0.004 ± 0.303) × 10⁻¹
        {3.192, 0.004, 0.111, 1e-1}, // (3.192 ± 0.004 ± 0.111) × 10⁻¹
        {2.098, 0.016, 0.200, 1e-1}, // (2.098 ± 0.016 ± 0.200) × 10⁻¹
        {2.750, 0.002, 0.188, 1e-1}, // (2.750 ± 0.002 ± 0.188) × 10⁻¹
        {3.734, 0.040, 0.213, 1e-2}, // (3.734 ± 0.040 ± 0.213) × 10⁻²
        {1.807, 0.005, 0.102, 1e-1}, // (1.807 ± 0.005 ± 0.102) × 10⁻¹
        {1.980, 0.012, 0.082, 1e-2}, // (1.980 ± 0.012 ± 0.082) × 10⁻²
        {1.846, 0.046, 0.122, 1e-3}  // (1.846 ± 0.046 ± 0.122) × 10⁻³
    };

    double meanPtAt13TeV[9];
    double meanPtAt13TeV_err[9];
    TGraphErrors *gMeanPt[9];
    TGraphErrors *gMeanPtvsMassMesons = new TGraphErrors();
    TGraphErrors *gMeanPtvsMassBaryons = new TGraphErrors();
    for (int i = 0; i < totalParticles; i++)
    {
        gMeanPt[i] = (TGraphErrors *)flightFlavourHadrons->Get(Form("Table %d/Graph1D_y1", 26 + i));
        if (gMeanPt[i] == nullptr)
        {
            cout << "Error reading graph for " << particlesLatex[i] << endl;
            return;
        }
        gMeanPt[i]->SetMarkerColor(colors[i]);
        gMeanPt[i]->SetLineColor(colors[i]);
        gMeanPt[i]->SetMarkerStyle(markers[i]);
        gMeanPt[i]->SetMarkerSize(1);

        // Now from graph store the value of mean pT at 13 TeV
        int nPoints = gMeanPt[i]->GetN();
        double *x = gMeanPt[i]->GetX();
        double *y = gMeanPt[i]->GetY();
        for (int j = 0; j < nPoints; j++)
        {
            if (x[j] == 13.0)
            {
                meanPtAt13TeV[i] = y[j];
                meanPtAt13TeV_err[i] = gMeanPt[i]->GetErrorY(j);

                cout << "Mean pT of " << particlesLatex[i] << " at 13 TeV is " << meanPtAt13TeV[i] << " +- " << meanPtAt13TeV_err[i] << endl;
                break;
            }
        }
        if (i < 5)
        {
            gMeanPtvsMassMesons->SetPoint(i, particleMass[i], meanPtAt13TeV[i]);
            gMeanPtvsMassMesons->SetPointError(i, 0, meanPtAt13TeV_err[i]);
        }
        else
        {
            gMeanPtvsMassBaryons->SetPoint(i - 5, particleMass[i], meanPtAt13TeV[i]);
            gMeanPtvsMassBaryons->SetPointError(i - 5, 0, meanPtAt13TeV_err[i]);
        }
    }
    SetGrapherrorStyle(gMeanPtvsMassMesons);
    gMeanPtvsMassMesons->SetMarkerStyle(22);
    gMeanPtvsMassMesons->SetMarkerColor(kMagenta);
    gMeanPtvsMassMesons->SetLineColor(kMagenta);
    gMeanPtvsMassMesons->SetMarkerSize(1.5);
    SetGrapherrorStyle(gMeanPtvsMassBaryons);
    gMeanPtvsMassBaryons->SetMarkerStyle(34);
    gMeanPtvsMassBaryons->SetMarkerColor(kRed);
    gMeanPtvsMassBaryons->SetLineColor(kRed);
    gMeanPtvsMassBaryons->SetMarkerSize(1.5);

    TCanvas *cMeanPt = new TCanvas("cMeanPt", "Mean pT vs mass", 720, 720);
    SetCanvasStyle(cMeanPt, 0.14, 0.03, 0.05, 0.14);
    gMeanPtvsMassMesons->GetXaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    gMeanPtvsMassMesons->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/#it{c})");
    gMeanPtvsMassMesons->GetYaxis()->SetTitleOffset(1.3);
    gMeanPtvsMassMesons->SetMinimum(0.1);
    gMeanPtvsMassMesons->SetMaximum(3.55);
    gMeanPtvsMassMesons->GetXaxis()->SetLimits(0, 1.99);
    gMeanPtvsMassMesons->Draw("AP");
    gMeanPtvsMassBaryons->Draw("P SAME");

    TF1 *pol1_meson = new TF1("pol1_meson", "pol1", 0.1, 1.8);
    pol1_meson->SetLineColor(kYellow + 2);
    pol1_meson->SetLineStyle(2);
    gMeanPtvsMassMesons->Fit(pol1_meson, "R");

    TF1 *pol1_baryon = new TF1("pol1_baryon", "pol1", 0.1, 1.8);
    pol1_baryon->SetLineColor(kYellow + 2);
    pol1_baryon->SetLineStyle(2);
    gMeanPtvsMassBaryons->Fit(pol1_baryon, "R");

    // Draw the last marker (f2(1525)) and its error bar
    double f2_mass = 1.5173;
    double f2_meanpt = 1.6633;
    double f2_err = 0.043027;
    double f2_meanpt2 = f1525MeanPt_Coherent;
    double f2_err2 = f1525MeanPt_Coherent_err;
    int f2_marker = 20;        // choose a unique marker style for f2(1525)
    int f2_color = kGreen + 1; // choose a unique color for f2(1525)
    int f2_color2 = kCyan + 2;
    int f2_marker2 = 21;
    TMarker *marker_f2 = new TMarker(f2_mass, f2_meanpt, f2_marker);
    marker_f2->SetMarkerColor(f2_color);
    marker_f2->SetMarkerSize(1.5);
    marker_f2->Draw("SAME");
    TLine *errBar_f2 = new TLine(f2_mass, f2_meanpt - f2_err, f2_mass, f2_meanpt + f2_err);
    errBar_f2->SetLineColor(f2_color);
    errBar_f2->SetLineWidth(2);
    errBar_f2->Draw("SAME");
    double cap_f2 = 0.01;
    TLine *capLow_f2 = new TLine(f2_mass - cap_f2, f2_meanpt - f2_err, f2_mass + cap_f2, f2_meanpt - f2_err);
    TLine *capHigh_f2 = new TLine(f2_mass - cap_f2, f2_meanpt + f2_err, f2_mass + cap_f2, f2_meanpt + f2_err);
    capLow_f2->SetLineColor(f2_color);
    capHigh_f2->SetLineColor(f2_color);
    capLow_f2->SetLineWidth(2);
    capHigh_f2->SetLineWidth(2);
    capLow_f2->Draw("SAME");
    capHigh_f2->Draw("SAME");

    // TMarker *marker_f2_2 = new TMarker(f2_mass, f2_meanpt2, f2_marker2);
    // marker_f2_2->SetMarkerColor(f2_color2);
    // marker_f2_2->SetMarkerSize(1.5);
    // marker_f2_2->Draw("SAME");
    // TLine *errBar_f2_2 = new TLine(f2_mass, f2_meanpt2 - f2_err2, f2_mass, f2_meanpt2 + f2_err2);
    // errBar_f2_2->SetLineColor(f2_color2);
    // errBar_f2_2->SetLineWidth(2);
    // errBar_f2_2->Draw("SAME");
    // TLine *capLow_f2_2 = new TLine(f2_mass - cap_f2, f2_meanpt2 - f2_err2, f2_mass + cap_f2, f2_meanpt2 - f2_err2);
    // TLine *capHigh_f2_2 = new TLine(f2_mass - cap_f2, f2_meanpt2 + f2_err2, f2_mass + cap_f2, f2_meanpt2 + f2_err2);
    // capLow_f2_2->SetLineColor(f2_color2);
    // capHigh_f2_2->SetLineColor(f2_color2);
    // capLow_f2_2->SetLineWidth(2);
    // capHigh_f2_2->SetLineWidth(2);
    // capLow_f2_2->Draw("SAME");
    // capHigh_f2_2->Draw("SAME");

    // Draw the last marker (f0(1710)) and its error bar
    double f0_mass = 1.710;
    double f0_meanpt = 2.35239;
    double f0_err = 0.100896;
    double f0_meanpt2 = f1710MeanPt_Coherent;
    double f0_err2 = f1710MeanPt_Coherent_err;
    int f0_marker = 21;   // choose a unique marker style for f0(1710)
    int f0_marker2 = 20;
    int f0_color = kBlue; // choose a unique color for f0(1710)
    int f0_color2 = kBrown;
    TMarker *marker_f0 = new TMarker(f0_mass, f0_meanpt, f0_marker);
    marker_f0->SetMarkerColor(f0_color);
    marker_f0->SetMarkerSize(1.5);
    marker_f0->Draw("SAME");
    TLine *errBar_f0 = new TLine(f0_mass, f0_meanpt - f0_err, f0_mass, f0_meanpt + f0_err);
    errBar_f0->SetLineColor(f0_color);
    errBar_f0->SetLineWidth(2);
    errBar_f0->Draw("SAME");
    double cap_f0 = 0.01;
    TLine *capLow_f0 = new TLine(f0_mass - cap_f0, f0_meanpt - f0_err, f0_mass + cap_f0, f0_meanpt - f0_err);
    TLine *capHigh_f0 = new TLine(f0_mass - cap_f0, f0_meanpt + f0_err, f0_mass + cap_f0, f0_meanpt + f0_err);
    capLow_f0->SetLineColor(f0_color);
    capHigh_f0->SetLineColor(f0_color);
    capLow_f0->SetLineWidth(2);
    capHigh_f0->SetLineWidth(2);
    capLow_f0->Draw("SAME");
    capHigh_f0->Draw("SAME");

    TMarker *marker_f0_2 = new TMarker(f0_mass, f0_meanpt2, f0_marker2);
    marker_f0_2->SetMarkerColor(f0_color2);
    marker_f0_2->SetMarkerSize(1.5);
    marker_f0_2->Draw("SAME");
    TLine *errBar_f0_2 = new TLine(f0_mass, f0_meanpt2 - f0_err2, f0_mass, f0_meanpt2 + f0_err2);
    errBar_f0_2->SetLineColor(f0_color2);
    errBar_f0_2->SetLineWidth(2);
    errBar_f0_2->Draw("SAME");
    TLine *capLow_f0_2 = new TLine(f0_mass - cap_f0, f0_meanpt2 - f0_err2, f0_mass + cap_f0, f0_meanpt2 - f0_err2);
    TLine *capHigh_f0_2 = new TLine(f0_mass - cap_f0, f0_meanpt2 + f0_err2, f0_mass + cap_f0, f0_meanpt2 + f0_err2);
    capLow_f0_2->SetLineColor(f0_color2);
    capHigh_f0_2->SetLineColor(f0_color2);
    capLow_f0_2->SetLineWidth(2);
    capHigh_f0_2->SetLineWidth(2);
    capLow_f0_2->Draw("SAME");
    capHigh_f0_2->Draw("SAME");
    // Add particle names below each point
    TLatex latex;
    latex.SetTextAlign(22);
    latex.SetTextSize(0.035);
    for (int i = 0; i < totalParticles; i++)
    {
        double x = particleMass[i];
        double y = meanPtAt13TeV[i];
        // latex.SetTextColor(colors[i]);
        if (i == 2)
            latex.DrawLatex(x + 0.08, y + 0.18, particlesLatex[i].c_str());
        else if (i <= 4)
            latex.DrawLatex(x, y + 0.18, particlesLatex[i].c_str());
        else
            latex.DrawLatex(x, y - 0.16, particlesLatex[i].c_str());
    }
    // latex.SetTextColor(kRed); // match f2(1525) marker color
    latex.DrawLatex(1.5173, f1525MeanPt_Default + 0.22, "f'_{2}(1525)");
    // latex.SetTextColor(kBlue); // match f0(1710) marker color
    latex.DrawLatex(1.710, f1710MeanPt_Default + 0.29, "f_{0}(1710)");

    TLegend *legend4 = new TLegend(0.18, 0.65, 0.65, 0.92);
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->SetTextSize(0.035);
    legend4->AddEntry(gMeanPtvsMassMesons, "Mesons (13 TeV)", "p");
    legend4->AddEntry(gMeanPtvsMassBaryons, "Baryons (13 TeV)", "p");
    legend4->AddEntry(marker_f2, "f'_{2}(1525)", "p");
    legend4->AddEntry(marker_f0, "f_{0}(1710) (Default fit)", "p");
    legend4->AddEntry(marker_f0_2, "f_{0}(1710) (Coherent fit)", "p");
    legend4->AddEntry(pol1_meson, "Pol 1", "l");
    legend4->Draw();
    cMeanPt->SaveAs((path + "MeanPt_vs_Mass_compare.png").c_str());
}