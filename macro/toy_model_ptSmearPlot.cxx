#include <TGenPhaseSpace.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "src/style.h"
#include <Math/VectorUtil.h>
#include "spectra/YieldMean.C"

using namespace std;
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void toy_model_ptSmearPlot()
{
    gStyle->SetOptStat(0);
    TString path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra/";
    TFile *f1710 = new TFile(path + "pTSmearing/f1710_ptSmearing_results.root", "read");
    TFile *f1525 = new TFile(path + "pTSmearing/f1525_ptSmearing_results.root", "read");
    TFile *fEffReweighted = new TFile(path + "ReweightedSpectra.root", "read");
    if (f1710->IsZombie() || f1525->IsZombie() || fEffReweighted->IsZombie())
    {
        cout << "Error opening pt smearing results file" << endl;
        return;
    }

    TFile *fout = new TFile(path + "pTSmearing/ptSmearing_AllPlots.root", "recreate");
    TH1F *hrec_pT_1710 = (TH1F *)f1710->Get("recpt");
    TH1F *hrec_pT_1525 = (TH1F *)f1525->Get("recpt");
    TH1F *hEff_1710 = (TH1F *)f1710->Get("eff");
    TH1F *hEff_1525 = (TH1F *)f1525->Get("eff");
    TH1F *hK0sEff = (TH1F *)f1710->Get("K0s_eff");
    hK0sEff->Write("K0s_efficiency_MC");
    TH1F *hEffReweighted_1710 = (TH1F *)fEffReweighted->Get("Eff_f0Reweighted");
    TH1F *hEffReweighted_1525 = (TH1F *)fEffReweighted->Get("Eff_f2Reweighted");
    if (hrec_pT_1710 == nullptr || hrec_pT_1525 == nullptr || hEff_1710 == nullptr || hEff_1525 == nullptr)
    {
        cout << "Error: Histograms not found in the file." << endl;
        return;
    }
    cout << "Bins in toy MC " << hEff_1710->GetNbinsX() << " , bin in reweighted MC " << hEffReweighted_1710->GetNbinsX() << endl;

    TCanvas *cRecpT = new TCanvas("cRecpT", "Reconstructed p_{T} Distribution Comparison", 720, 720);
    SetCanvasStyle(cRecpT, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hrec_pT_1710);
    hrec_pT_1710->SetLineColor(kBlue);
    hrec_pT_1710->SetMarkerColor(kBlue);
    hrec_pT_1710->GetYaxis()->SetMaxDigits(3);
    hrec_pT_1710->Write("toyModel_pt_f01710");
    hrec_pT_1710->Draw();
    hrec_pT_1525->SetLineColor(kRed);
    hrec_pT_1525->SetMarkerColor(kRed);
    hrec_pT_1525->Write("toyModel_pt_f21525");
    hrec_pT_1525->Draw("SAME");

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->SetTextFont(42);
    leg->AddEntry(hrec_pT_1525, "f_{2}(1525)", "lep");
    leg->AddEntry(hrec_pT_1710, "f_{0}(1710)", "lep");
    leg->Draw();
    cRecpT->SaveAs(path + "pTSmearing/f0_f2_rec_pT_compare.png");

    TCanvas *cEffResults = new TCanvas("cEffResults", "Average Efficiency vs p_{T} bin Comparison", 720, 720);
    SetCanvasStyle(cEffResults, 0.15, 0.05, 0.05, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cEffResults, pad1Size, pad2Size);
    cEffResults->cd(1);
    SetHistoQA(hEff_1710);
    // gPad->SetLogy();
    hEff_1710->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hEff_1710->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hEff_1710->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hEff_1710->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hEff_1710->GetYaxis()->SetTitleOffset(1.3 * pad1Size);
    hEff_1710->SetLineColor(kBlue);
    hEff_1710->SetMarkerColor(kBlue);
    hEff_1710->GetYaxis()->SetMaxDigits(3);
    // hEff_1710->SetMaximum(0.29);
    hEff_1710->SetMaximum(0.06);
    hEff_1710->SetMinimum(0.001);
    hEff_1710->Write("toyModel_efficiency_f01710");
    hEff_1710->Draw();
    hEff_1525->SetLineColor(kRed);
    hEff_1525->SetMarkerColor(kRed);
    hEff_1525->Write("toyModel_efficiency_f21525");
    hEff_1525->Draw("SAME");
    hEffReweighted_1710->SetLineColor(kGreen + 2);
    hEffReweighted_1710->SetMarkerColor(kGreen + 2);
    hEffReweighted_1710->Write("MC_efficiency_f01710");
    hEffReweighted_1710->Draw("HIST SAME");
    hEffReweighted_1525->SetLineColor(kCyan);
    hEffReweighted_1525->SetMarkerColor(kCyan);
    hEffReweighted_1525->Write("MC_efficiency_f21525");
    hEffReweighted_1525->Draw("HIST SAME");
    TLegend *leg2 = new TLegend(0.55, 0.4, 0.88, 0.6);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.03);
    leg2->SetTextFont(42);
    leg2->AddEntry(hEff_1525, "f_{2}(1525) (Toy model)", "l");
    leg2->AddEntry(hEff_1710, "f_{0}(1710) (Toy model)", "l");
    leg2->AddEntry(hEffReweighted_1525, "f_{2}(1525) (MC)", "l");
    leg2->AddEntry(hEffReweighted_1710, "f_{0}(1710) (MC)", "l");
    leg2->Draw();

    cEffResults->cd(2);
    TH1F *hEffRatio1710 = (TH1F *)hEffReweighted_1710->Clone("hEffRatio1710");
    TH1F *hEffRatio1525 = (TH1F *)hEffReweighted_1525->Clone("hEffRatio1525");
    for (int ibin = 1; ibin <= hEffRatio1710->GetNbinsX(); ibin++)
    {
        if (hEff_1710->GetBinContent(ibin) != 0)
        {
            double binToy = hEff_1710->GetBinContent(ibin);
            double binMC = hEffReweighted_1710->GetBinContent(ibin);
            double ratioMCbyToy = binMC / binToy;
            hEffRatio1710->SetBinContent(ibin, ratioMCbyToy);
            hEffRatio1710->SetBinError(ibin, 0);
            cout << "pT bins " << ibin << " , Toy Eff " << binToy << " , MC Eff " << binMC << " , Ratio " << ratioMCbyToy << endl;
        }
        if (hEff_1525->GetBinContent(ibin) == 0)
            continue;
        double binToy1525 = hEff_1525->GetBinContent(ibin);
        double binMC1525 = hEffReweighted_1525->GetBinContent(ibin);
        double ratioMCbyToy1525 = binMC1525 / binToy1525;
        hEffRatio1525->SetBinContent(ibin, ratioMCbyToy1525);
        hEffRatio1525->SetBinError(ibin, 0);
    }
    SetHistoQA(hEffRatio1710);
    hEffRatio1710->GetYaxis()->SetTitle("MC / Toy model");
    hEffRatio1710->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    hEffRatio1710->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hEffRatio1710->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hEffRatio1710->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hEffRatio1710->GetYaxis()->SetTitleOffset(1.5 * pad2Size);
    hEffRatio1710->SetLineColor(kBlue);
    hEffRatio1710->SetMarkerColor(kBlue);
    hEffRatio1710->GetYaxis()->SetNdivisions(505);
    hEffRatio1710->SetMaximum(0.5);
    hEffRatio1710->SetMinimum(0.05);
    hEffRatio1710->Draw("HIST");
    SetHistoQA(hEffRatio1525);
    hEffRatio1525->SetLineColor(kRed);
    hEffRatio1525->SetMarkerColor(kRed);
    hEffRatio1525->SetMaximum(0.5);
    hEffRatio1525->SetMinimum(0.05);
    hEffRatio1525->Draw("HIST SAME");
    TLine *lineRatio = new TLine(0, 1, 15, 1);
    lineRatio->SetLineStyle(1);
    lineRatio->SetLineColor(kBlack);
    lineRatio->Draw();
    cEffResults->SaveAs(path + "pTSmearing/f0_f2_Average_eff_vs_pT.png");

    // Open the raw spectra and correct using toy model efficiency
    string path2 = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/";
    TFile *inputFile = new TFile((path2 + "FitParamdefault.root").c_str(), "READ");
    if (inputFile->IsZombie())
    {
        cerr << "Error opening file" << endl;
    }
    TH1F *hRawf0 = (TH1F *)inputFile->Get("Mult_0_100/hYield_1710");
    TH1F *hRawf2 = (TH1F *)inputFile->Get("Mult_0_100/hYield_1525");
    if (hRawf0 == nullptr || hRawf2 == nullptr)
    {
        cout << "Error reading raw spectra histograms" << endl;
        return;
    }
    TH1F *hYield1710Corrected = (TH1F *)hRawf0->Clone("hYield1710Corrected");
    TH1F *hYield1525Corrected = (TH1F *)hRawf2->Clone("hYield1525Corrected");

    int totalBins = hRawf0->GetNbinsX();
    for (int ibins = 0; ibins < totalBins; ibins++)
    {
        double recYieldf0 = hRawf0->GetBinContent(ibins + 1);
        double recYieldErrorf0 = hRawf0->GetBinError(ibins + 1);
        double efficiencyf0 = hEff_1710->GetBinContent(ibins + 1);
        double recYieldf2 = hRawf2->GetBinContent(ibins + 1);
        double recYieldErrorf2 = hRawf2->GetBinError(ibins + 1);
        double efficiencyf2 = hEff_1525->GetBinContent(ibins + 1);

        if (efficiencyf0 > 0)
        {
            hYield1710Corrected->SetBinContent(ibins + 1, recYieldf0 / efficiencyf0);
            hYield1710Corrected->SetBinError(ibins + 1, recYieldErrorf0);
        }
        if (efficiencyf2 > 0)
        {
            hYield1525Corrected->SetBinContent(ibins + 1, recYieldf2 / efficiencyf2);
            hYield1525Corrected->SetBinError(ibins + 1, recYieldErrorf2);
        }
    }

    TH1F *h1 = (TH1F *)hYield1525Corrected->Clone("h1");
    TH1F *h2 = (TH1F *)hYield1525Corrected->Clone("h2");
    TH1F *h3 = (TH1F *)hYield1710Corrected->Clone("h3");
    TH1F *h4 = (TH1F *)hYield1710Corrected->Clone("h4");

    for (int i = 1; i <= h2->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr = (0.1 * h2->GetBinContent(i));
        h2->SetBinError(i, systemerr);
    }
    /*************meanpT*****************byresonance*******************package*************************/
    Double_t min = 0.0;
    Double_t max = 15.0;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "REBMS0+";
    TString logfilename = "log.root";
    Double_t minfit = 1.0;
    Double_t maxfit = 10.0;

    TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
    fitFcn->SetParameter(0, 5.0);
    // fitFcn->SetParameter(1, 0.05);
    fitFcn->SetParameter(1, 0.5);
    fitFcn->FixParameter(2, 1.525);
    fitFcn->SetParameter(3, 0.35);
    fitFcn->SetParNames("n", "dn/dy", "mass", "T");

    TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    TF1 *fitFcn2 = new TF1("fitfunc2", FuncLavy, 0.0, 15.0, 4);
    fitFcn2->SetParameter(0, 10.0);
    // fitFcn2->SetParameter(1, 0.05);
    fitFcn2->SetParameter(1, 0.005);
    fitFcn2->FixParameter(2, 1.710);
    fitFcn2->SetParameter(3, 0.5);
    fitFcn2->SetParLimits(3, 0.1, 1.0);
    fitFcn2->SetParNames("n", "dn/dy", "mass", "T");

    TH1 *hout2 = YieldMean(h3, h4, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    TCanvas *cCorrectedf2Fit = new TCanvas("cCorrectedf2Fit", "Corrected #it{p}_{T} distribution with fit", 720, 720);
    SetCanvasStyle(cCorrectedf2Fit, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(h1);
    h1->SetMarkerSize(1.5);
    // h1->SetMaximum(h1->GetMaximum() * 12);
    h1->SetMaximum(1e-3);
    h1->SetMinimum(1e-7);
    h1->SetMarkerColor(kBlue);
    h1->SetLineColor(kBlue);
    h1->Draw("pe");
    fitFcn->SetLineWidth(2);
    fitFcn->SetLineStyle(2);
    fitFcn->SetLineColor(kRed);
    fitFcn->Draw("l same");
    TLegend *legend22 = new TLegend(0.55, 0.7, 0.88, 0.88);
    legend22->SetBorderSize(0);
    legend22->SetFillStyle(0);
    legend22->SetTextSize(0.035);
    legend22->SetTextFont(42);
    legend22->AddEntry(h1, "f_{2}(1525) (13.6 TeV)", "p");
    legend22->AddEntry(fitFcn, "Levy-Tsallis fit", "l");
    legend22->Draw();
    cCorrectedf2Fit->SaveAs(path + "pTSmearing/LevyFitf2_ptSmearing.png");

    TCanvas *cCorrectedf0Fit = new TCanvas("cCorrectedf0Fit", "Corrected #it{p}_{T} distribution with fit", 720, 720);
    SetCanvasStyle(cCorrectedf0Fit, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(h3);
    h3->SetMarkerSize(1.5);
    // h3->SetMaximum(h3->GetMaximum() * 12);
    h3->SetMaximum(1e-4);
    h3->SetMinimum(1e-7);
    h3->SetMarkerColor(kBlue);
    h3->SetLineColor(kBlue);
    h3->Draw("pe");
    fitFcn2->SetLineWidth(2);
    fitFcn2->SetLineStyle(2);
    fitFcn2->SetLineColor(kRed);
    fitFcn2->Draw("l same");
    TLegend *legend44 = new TLegend(0.55, 0.7, 0.88, 0.88);
    legend44->SetBorderSize(0);
    legend44->SetFillStyle(0);
    legend44->SetTextSize(0.035);
    legend44->SetTextFont(42);
    legend44->AddEntry(h3, "f_{0}(1710) (13.6 TeV)", "p");
    legend44->AddEntry(fitFcn2, "Levy-Tsallis fit", "l");
    legend44->Draw();
    cCorrectedf0Fit->SaveAs(path + "pTSmearing/LevyFitf0_ptSmearing.png");

    // Now plot the <pT> as a function of particle masses
    TFile *flightFlavourHadrons = new TFile("spectra/LightFlavourHadronsProduction.root", "read");
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
    gMeanPtvsMassMesons->SetMinimum(0.0);
    gMeanPtvsMassMesons->SetMaximum(3.75);
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
    double f2_meanpt = hout->GetBinContent(5);
    double f2_err = hout->GetBinContent(6);
    int f2_marker = 20;        // choose a unique marker style for f2(1525)
    int f2_color = kGreen + 1; // choose a unique color for f2(1525)
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

    // Draw the last marker (f0(1710)) and its error bar
    double f0_mass = 1.710;
    double f0_meanpt = hout2->GetBinContent(5);
    double f0_err = hout2->GetBinContent(6);
    int f0_marker = 21;   // choose a unique marker style for f0(1710)
    int f0_color = kBlue; // choose a unique color for f0(1710)
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
    latex.DrawLatex(1.5173, hout->GetBinContent(5) + 0.22, "f'_{2}(1525)");
    // latex.SetTextColor(kBlue); // match f0(1710) marker color
    latex.DrawLatex(1.710, hout2->GetBinContent(5) + 0.29, "f_{0}(1710)");

    TLegend *legend4 = new TLegend(0.18, 0.65, 0.65, 0.92);
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->SetTextSize(0.035);
    legend4->AddEntry(gMeanPtvsMassMesons, "Mesons (13 TeV)", "p");
    legend4->AddEntry(gMeanPtvsMassBaryons, "Baryons (13 TeV)", "p");
    legend4->AddEntry(marker_f2, "f'_{2}(1525) (13.6 TeV)", "p");
    legend4->AddEntry(marker_f0, "f_{0}(1710) (13.6 TeV)", "p");
    legend4->AddEntry(pol1_meson, "Pol 1", "l");
    legend4->Draw();
    cMeanPt->SaveAs(path + "/pTSmearing/MeanPt_vs_Mass.png");

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
    pad1->SetLeftMargin(0.12);
    pad2->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.13);
}