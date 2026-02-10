#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"
#include "../spectra/YieldMean.C"

using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

// Function to find the highest available index in the root file
int FindHighestIndex(TFile *file, const string &baseHistoName)
{
    int maxIndex = -1;
    for (int i = 10; i >= 0; i--) // Check from i10 down to i0
    {
        string histoName = baseHistoName + "i" + to_string(i);
        TObject *obj = file->Get(histoName.c_str());
        if (obj != nullptr)
        {
            maxIndex = i;
            cout << "Highest index found: " << maxIndex << endl;
            break; // Found the highest index
        }
    }
    return maxIndex;
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void plotReweightedSpectra_single()
{
    bool otherQAPlots = true;
    gStyle->SetOptStat(0);
    // gStyle->SetOptFit(1111);

    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/";

    string savePath = path + "/mult_0-100/Spectra";
    TFile *fReweightf0 = new TFile((path + "mult_0-100/Spectra/ReweighFacf0_Default2.root").c_str(), "read");
    TFile *fReweightf2 = new TFile((path + "mult_0-100/Spectra/ReweighFacf2_Default2.root").c_str(), "read");

    // TFile *fReweightf0 = new TFile(Form("%s/ReweighFacf0_%s.root", savePath.c_str(), CurrentVariation.c_str()), "read");
    // TFile *fReweightf2 = new TFile(Form("%s/ReweighFacf2_%s.root", savePath.c_str(), CurrentVariation.c_str()), "read");

    if (fReweightf0->IsZombie() || fReweightf2->IsZombie())
    {
        cout << "Error opening reweighting files" << endl;
        return;
    }
    TFile *fOutput = new TFile(Form("%s/ReweightedSpectra.root", savePath.c_str()), "RECREATE");

    // Find the highest available index for reweighted histograms
    int maxIndexReweightedf0 = FindHighestIndex(fReweightf0, "Genf17102_proj_1_");
    int maxIndexReweightedf2 = FindHighestIndex(fReweightf2, "Genf17102_proj_1_");
    if (maxIndexReweightedf0 == -1 || maxIndexReweightedf2 == -1)
    {
        cout << "Error: No reweighted histogram with pattern Genf17102_proj_1_i* found in file" << endl;
        return;
    }
    cout << "Using index i" << maxIndexReweightedf0 << " for reweighted histograms" << endl;

    string indexStr = "i" + to_string(maxIndexReweightedf0);
    string indexStr2 = "i" + to_string(maxIndexReweightedf2);
    TH1F *hGenReweighted = (TH1F *)fReweightf0->Get(Form("Genf17102_proj_1_%s", indexStr.c_str()));
    TH1F *hRecReweighted = (TH1F *)fReweightf0->Get(Form("Recf1710_pt2_proj_1_%s", indexStr.c_str()));
    TH1F *hYieldReweighted = (TH1F *)fReweightf0->Get(Form("hYield1710Corrected_%s", indexStr.c_str()));
    TH1F *hGenReweighted2 = (TH1F *)fReweightf2->Get(Form("Genf17102_proj_1_%s", indexStr2.c_str()));
    TH1F *hRecReweighted2 = (TH1F *)fReweightf2->Get(Form("Recf1710_pt2_proj_1_%s", indexStr2.c_str()));
    TH1F *hYieldReweighted2 = (TH1F *)fReweightf2->Get(Form("hYield1525Corrected_%s", indexStr2.c_str()));

    // // Ensure all histograms have positive values for log scale drawing
    // for (int i = 1; i <= hGenReweighted->GetNbinsX(); i++)
    // {
    //     if (hGenReweighted->GetBinContent(i) <= 0)
    //         cout << "Warning: Bin " << i << " in hGenReweighted has non-positive content with value " << hGenReweighted->GetBinContent(i) << endl;
    //     if (hRecReweighted->GetBinContent(i) <= 0)
    //         cout << "Warning: Bin " << i << " in hRecReweighted has non-positive content with value " << hRecReweighted->GetBinContent(i) << endl;
    // }
    // for (int i = 1; i <= hGenReweighted2->GetNbinsX(); i++)
    // {
    //     if (hGenReweighted2->GetBinContent(i) <= 0)
    //         cout << "Warning: Bin " << i << " in hGenReweighted2 has non-positive content with value " << hGenReweighted2->GetBinContent(i) << endl;
    //     if (hRecReweighted2->GetBinContent(i) <= 0)
    //         cout << "Warning: Bin " << i << " in hRecReweighted2 has non-positive content with value " << hRecReweighted2->GetBinContent(i) << endl;
    // }

    TH1F *hGenUnweighted = (TH1F *)fReweightf0->Get("Genf17102_proj_1_i0");
    TH1F *hRecUnweighted = (TH1F *)fReweightf0->Get("Recf1710_pt2_proj_1_i0");
    TH1F *hGenUnweighted2 = (TH1F *)fReweightf2->Get("Genf17102_proj_1_i0");
    TH1F *hRecUnweighted2 = (TH1F *)fReweightf2->Get("Recf1710_pt2_proj_1_i0");

    if (hGenReweighted == nullptr || hRecReweighted == nullptr || hYieldReweighted == nullptr)
    {
        cout << "Error reading reweighted histograms from file " << Form("%s/ReweightedSpectra.root", savePath.c_str()) << endl;
        return;
    }

    TCanvas *cReweighted = new TCanvas("cReweighted", "Reweighted Efficiency for f0(1710)", 720, 720);
    SetCanvasStyle(cReweighted, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hGenReweighted);

    // Clone unweighted histograms before scaling for plotting; use clones for efficiency counts
    TH1F *hGenUnweightedCounts = (TH1F *)hGenUnweighted->Clone("hGenUnweightedCounts");
    TH1F *hRecUnweightedCounts = (TH1F *)hRecUnweighted->Clone("hRecUnweightedCounts");
    TH1F *hGenUnweightedCounts2 = (TH1F *)hGenUnweighted2->Clone("hGenUnweightedCounts2");
    TH1F *hRecUnweightedCounts2 = (TH1F *)hRecUnweighted2->Clone("hRecUnweightedCounts2");
    if (hGenUnweightedCounts)
        hGenUnweightedCounts->Sumw2();
    if (hRecUnweightedCounts)
        hRecUnweightedCounts->Sumw2();
    if (hGenUnweightedCounts2)
        hGenUnweightedCounts2->Sumw2();
    if (hRecUnweightedCounts2)
        hRecUnweightedCounts2->Sumw2();
    gPad->SetLogy();
    hGenReweighted->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hGenReweighted->GetYaxis()->SetTitle("1/N_{evt} #times d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hGenReweighted->GetYaxis()->SetTitleOffset(1.7);
    hGenReweighted->SetMaximum(hGenReweighted->GetMaximum() * 1.5);
    hGenReweighted->SetMarkerStyle(53);
    hGenReweighted->SetMarkerColor(kGreen);
    hGenReweighted->SetLineColor(kGreen);
    hGenReweighted->GetXaxis()->SetRangeUser(0.0, 15.5);
    hGenReweighted->SetMinimum(5e-9);
    hGenReweighted->SetMaximum(hGenReweighted->GetMaximum() * 500);
    hGenReweighted->SetMarkerSize(1.5);
    hGenReweighted->Draw("pe");
    SetHistoQA(hRecReweighted);
    hRecReweighted->SetMarkerStyle(53);
    hRecReweighted->SetMarkerColor(kRed);
    hRecReweighted->SetLineColor(kRed);
    hRecReweighted->SetMarkerSize(1.5);
    hRecReweighted->Draw("pe same");
    SetHistoQA(hYieldReweighted);
    hYieldReweighted->SetMarkerStyle(20);
    hYieldReweighted->SetMarkerColor(kBlack);
    hYieldReweighted->SetLineColor(kBlack);
    hYieldReweighted->SetMarkerSize(1.5);
    hYieldReweighted->Draw("pe same");
    double integralFactor = 20861874; // this is multiplicity in MC
    // double integralFactor = 1.0; // this is multiplicity in MC
    SetHistoQA(hGenUnweighted);
    // hGenUnweighted->Scale(1.0 / integralFactor);
    hGenUnweighted->SetMarkerStyle(53);
    hGenUnweighted->SetMarkerColor(kBlue);
    hGenUnweighted->SetLineColor(kBlue);
    hGenUnweighted->SetMarkerSize(1.5);
    hGenUnweighted->Draw("pe same");
    SetHistoQA(hRecUnweighted);
    // hRecUnweighted->Scale(1.0 / integralFactor);
    hRecUnweighted->SetMarkerStyle(53);
    hRecUnweighted->SetMarkerColor(kMagenta);
    hRecUnweighted->SetLineColor(kMagenta);
    hRecUnweighted->SetMarkerSize(1.5);
    hRecUnweighted->Draw("pe same");
    TLegend *legReweighted = new TLegend(0.35, 0.73, 0.9, 0.93);
    legReweighted->SetBorderSize(0);
    legReweighted->SetFillStyle(0);
    legReweighted->SetTextSize(0.03);
    legReweighted->SetHeader("pp, #sqrt{#it{s}} = 13.6 TeV");
    legReweighted->AddEntry(hGenUnweighted, "Generated distribution", "p");
    legReweighted->AddEntry(hRecUnweighted, "Reconstructed distribution", "p");
    legReweighted->AddEntry(hYieldReweighted, "Data (f_{0}(1710))", "p");
    legReweighted->AddEntry(hGenReweighted, "Reweighted Generated distribution", "p");
    legReweighted->AddEntry(hRecReweighted, "Reweighted Reconstructed distribution", "p");
    legReweighted->Draw();
    cReweighted->SaveAs((savePath + "/plots/ReweightedEfficiencyf0.png").c_str());

    TCanvas *cReweightedf2 = new TCanvas("cReweightedf2", "Reweighted Efficiency for f2'(1525)", 720, 720);
    SetCanvasStyle(cReweightedf2, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hGenReweighted2);
    gPad->SetLogy();
    hGenReweighted2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hGenReweighted2->GetYaxis()->SetTitle("1/N_{evt} #times d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hGenReweighted2->GetYaxis()->SetTitleOffset(1.7);
    // hGenReweighted2->SetMaximum(hGenReweighted2->GetMaximum() * 1.5);
    hGenReweighted2->SetMarkerStyle(53);
    hGenReweighted2->SetMarkerColor(kGreen);
    hGenReweighted2->SetLineColor(kGreen);
    hGenReweighted2->GetXaxis()->SetRangeUser(0.0, 15.5);
    hGenReweighted2->SetMaximum(hGenReweighted2->GetMaximum() * 500);
    hGenReweighted2->SetMinimum(2e-7);
    hGenReweighted2->SetMarkerSize(1.5);
    hGenReweighted2->Draw("pe");
    SetHistoQA(hRecReweighted2);
    hRecReweighted2->SetMarkerStyle(53);
    hRecReweighted2->SetMarkerColor(kRed);
    hRecReweighted2->SetLineColor(kRed);
    hRecReweighted2->SetMarkerSize(1.5);
    hRecReweighted2->Draw("pe same");
    SetHistoQA(hYieldReweighted2);
    hYieldReweighted2->SetMarkerStyle(20);
    hYieldReweighted2->SetMarkerColor(kBlack);
    hYieldReweighted2->SetLineColor(kBlack);
    hYieldReweighted2->SetMarkerSize(1.5);
    hYieldReweighted2->Draw("pe same");
    SetHistoQA(hGenUnweighted2);
    // hGenUnweighted2->Scale(1.0 / integralFactor);
    hGenUnweighted2->SetMarkerStyle(53);
    hGenUnweighted2->SetMarkerColor(kBlue);
    hGenUnweighted2->SetLineColor(kBlue);
    hGenUnweighted2->SetMarkerSize(1.5);
    hGenUnweighted2->Draw("pe same");
    SetHistoQA(hRecUnweighted2);
    // hRecUnweighted2->Scale(1.0 / integralFactor);
    hRecUnweighted2->SetMarkerStyle(53);
    hRecUnweighted2->SetMarkerColor(kMagenta);
    hRecUnweighted2->SetLineColor(kMagenta);
    hRecUnweighted2->SetMarkerSize(1.5);
    hRecUnweighted2->Draw("pe same");
    TLegend *legReweighted2 = new TLegend(0.35, 0.73, 0.9, 0.93);
    legReweighted2->SetBorderSize(0);
    legReweighted2->SetFillStyle(0);
    legReweighted2->SetTextSize(0.03);
    legReweighted2->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    legReweighted2->AddEntry(hGenUnweighted2, "Generated distribution", "p");
    legReweighted2->AddEntry(hRecUnweighted2, "Reconstructed distribution", "p");
    legReweighted2->AddEntry(hYieldReweighted2, "Data (f_{2}'(1525))", "p");
    legReweighted2->AddEntry(hGenReweighted2, "Reweighted Generated distribution", "p");
    legReweighted2->AddEntry(hRecReweighted2, "Reweighted Reconstructed distribution", "p");
    legReweighted2->Draw();
    cReweightedf2->SaveAs((savePath + "/plots/ReweightedEfficiencyf2.png").c_str());
    hYieldReweighted->Write("f01710_Reweighted_Yield");
    hYieldReweighted2->Write("f21525_Reweighted_Yield");

    ////****************************Levy-Tsallis fit**********************************////
    TFile *fSys = new TFile((path + "mult_0-100/Spectra/SystematicPlots/SystematicUncertainties.root").c_str(), "read");
    if (fSys->IsZombie())
    {
        cout << "Error opening systematic uncertainty file" << endl;
        return;
    }
    TH1F *hYieldSysf0 = (TH1F *)fSys->Get("TotalSys_Smooth_Yield1710");
    TH1F *hYieldSysf2 = (TH1F *)fSys->Get("TotalSys_Smooth_Yield1525");
    if (hYieldSysf0 == nullptr || hYieldSysf2 == nullptr)
    {
        cout << "Error reading systematic uncertainty histograms from file " << Form("%s/SystematicUncertainties.root", (path + "mult_0-100/Spectra/").c_str()) << endl;
        return;
    }

    cout << "Total bins in the yield histogram is " << hYieldReweighted->GetNbinsX() << endl;

    // For f2'(1525)
    TH1F *hf21 = (TH1F *)hYieldReweighted2->Clone("hf21");
    TH1F *hf22 = (TH1F *)hYieldReweighted2->Clone("hf22");

    // For f0(1710)
    TH1F *hf01 = (TH1F *)hYieldReweighted->Clone("hf01");
    TH1F *hf02 = (TH1F *)hYieldReweighted->Clone("hf02");

    // Enable error tracking for histograms with manual error setting
    hf01->Sumw2();
    hf02->Sumw2();
    hf21->Sumw2();
    hf22->Sumw2();

    for (int i = 1; i <= hf21->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double sys1710 = hYieldSysf0->GetBinContent(i);
        double sys1525 = hYieldSysf2->GetBinContent(i);
        double yield1710 = hf01->GetBinContent(i + 1);
        double yield1525 = hf21->GetBinContent(i + 1);
        cout << "Bin " << i << ": Yield f0(1710) = " << yield1710 << " with sys error = " << sys1710 << " and Yield f2'(1525) = " << yield1525 << " with sys error = " << sys1525 << endl;
        hf02->SetBinContent(i + 1, yield1710);
        hf02->SetBinError(i + 1, sys1710 * yield1710);
        hf22->SetBinContent(i + 1, yield1525);
        hf22->SetBinError(i + 1, sys1525 * yield1525);
    }
    Double_t min = 0.0;
    Double_t max = 15.0;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "REBMS0+";
    // Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 1.0;
    Double_t maxfit = 10.0;

    TF1 *fitFcnf0 = new TF1("fitfunc2", FuncLavy, 0.0, 15.0, 4);
    fitFcnf0->SetParameter(0, 5.0);
    // fitFcnf0->SetParameter(1, 0.05);
    fitFcnf0->SetParameter(1, 0.5);
    fitFcnf0->FixParameter(2, 1.710);
    fitFcnf0->SetParameter(3, 0.35);
    fitFcnf0->SetParNames("n", "dn/dy", "mass", "T");

    TF1 *fitFcnf2 = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
    fitFcnf2->SetParameter(0, 5.0);
    // fitFcnf2->SetParameter(1, 0.05);
    fitFcnf2->SetParameter(1, 0.5);
    fitFcnf2->FixParameter(2, 1.525);
    fitFcnf2->SetParameter(3, 0.35);
    fitFcnf2->SetParNames("n", "dn/dy", "mass", "T");

    TH1 *hout = YieldMean(hf01, hf02, fitFcnf0, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    cout << "Yield dN/dy = " << hout->GetBinContent(1) << " +- " << hout->GetBinContent(2) << endl;
    cout << "Mean pT = " << hout->GetBinContent(5) << " +- " << hout->GetBinContent(6) << endl;

    TH1 *hout2 = YieldMean(hf21, hf22, fitFcnf2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    cout << "Yield dN/dy of f2'(1525) = " << hout2->GetBinContent(1) << " +- " << hout2->GetBinContent(2) << endl;
    cout << "Mean pT of f2'(1525) = " << hout2->GetBinContent(5) << " +- " << hout2->GetBinContent(6) << endl;

    TCanvas *cCorrectedf0Fit = new TCanvas("cCorrectedf0Fit", "Corrected #it{p}_{T} distribution with fit", 720, 720);
    SetCanvasStyle(cCorrectedf0Fit, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(hf01);
    hf01->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hf01->GetYaxis()->SetTitle("BR #times 1/N_{evt} #times d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    // hf01->GetYaxis()->SetTitle("1/N_{evt} #times d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hf01->GetYaxis()->SetTitleOffset(1.7);
    hf01->SetMarkerSize(1.5);
    // hf01->SetMaximum(hf01->GetMaximum() * 2);
    hf01->SetMaximum(2e-3);
    hf01->SetMinimum(9e-8);
    hf01->SetMarkerColor(kBlue);
    hf01->SetLineColor(kBlue);
    hf01->Draw("pe");
    hf02->SetMarkerColor(kBlue);
    hf02->SetLineColor(kBlue);
    hf02->SetFillStyle(0);
    hf02->Draw("e2 same");
    fitFcnf0->SetLineWidth(2);
    fitFcnf0->SetLineStyle(2);
    fitFcnf0->SetLineColor(kRed);
    fitFcnf0->Draw("l same");

    // Create legend with physics information
    TLegend *leg = new TLegend(0.53, 0.7, 0.9, 0.93);
    leg->AddEntry((TObject *)0, "ALICE", "");
    leg->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    leg->AddEntry(hf01, "f_{0}(1710) spectra", "pe");
    leg->AddEntry(fitFcnf0, "Levy-Tsallis fit", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->Draw();
    cCorrectedf0Fit->SaveAs((savePath + "/plots/LevyFitf0_reweighted.png").c_str());

    TCanvas *cCorrectedf2Fit = new TCanvas("cCorrectedf2Fit", "Corrected #it{p}_{T} distribution with fit", 720, 720);
    SetCanvasStyle(cCorrectedf2Fit, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(hf21);
    hf21->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hf21->GetYaxis()->SetTitle("1/N_{evt} #times d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hf21->GetYaxis()->SetTitleOffset(1.7);
    hf21->SetMarkerSize(1.5);
    // hf21->SetMaximum(hf21->GetMaximum() * 7);
    hf21->SetMaximum(4e-2);
    hf21->SetMinimum(2e-7);
    hf21->SetMarkerColor(kBlue);
    hf21->SetLineColor(kBlue);
    hf21->Draw("pe");
    hf22->SetMarkerColor(kBlue);
    hf22->SetLineColor(kBlue);
    hf22->SetFillStyle(0);
    hf22->Draw("e2 same");
    fitFcnf2->SetLineWidth(2);
    fitFcnf2->SetLineStyle(2);
    fitFcnf2->SetLineColor(kRed);
    fitFcnf2->Draw("l same");

    // Create legend with physics information
    TLegend *leg2 = new TLegend(0.53, 0.7, 0.9, 0.93);
    leg2->AddEntry((TObject *)0, "ALICE", "");
    leg2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg2->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    leg2->AddEntry(hf21, "f_{2}'(1525) spectra", "pe");
    leg2->AddEntry(fitFcnf2, "Levy-Tsallis fit", "l");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.035);
    leg2->Draw();
    cCorrectedf2Fit->SaveAs((savePath + "/plots/LevyFitf2_reweighted.png").c_str());

    // /*
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
    // vector<vector<float>> dNdyvalues_13TeV = {
    //     {4.775, 0.001, 0.243, 1.0},  // 4.775 ± 0.001 ± 0.243
    //     {6.205, 0.004, 0.303, 1e-1}, // (6.205 ± 0.004 ± 0.303) × 10⁻¹
    //     {3.192, 0.004, 0.111, 1e-1}, // (3.192 ± 0.004 ± 0.111) × 10⁻¹
    //     {2.098, 0.016, 0.200, 1e-1}, // (2.098 ± 0.016 ± 0.200) × 10⁻¹
    //     {2.750, 0.002, 0.188, 1e-1}, // (2.750 ± 0.002 ± 0.188) × 10⁻¹
    //     {3.734, 0.040, 0.213, 1e-2}, // (3.734 ± 0.040 ± 0.213) × 10⁻²
    //     {1.807, 0.005, 0.102, 1e-1}, // (1.807 ± 0.005 ± 0.102) × 10⁻¹
    //     {1.980, 0.012, 0.082, 1e-2}, // (1.980 ± 0.012 ± 0.082) × 10⁻²
    //     {1.846, 0.046, 0.122, 1e-3}  // (1.846 ± 0.046 ± 0.122) × 10⁻³
    // };

    double meanPtAt13TeV[9];
    double meanPtAt13TeV_err[9];
    TGraphErrors *gMeanPt[9];
    TGraphAsymmErrors *gMeanPtvsMassMesons = new TGraphAsymmErrors();
    TGraphAsymmErrors *gMeanPtvsMassBaryons = new TGraphAsymmErrors();
    for (int i = 0; i < totalParticles; i++)
    {
        gMeanPt[i] = (TGraphErrors *)flightFlavourHadrons->Get(Form("Table %d/Graph1D_y1", 26 + i));
        cout << "Table " << 26 + i << " for particle " << particles[i] << endl;
        cout << "<pT> " << gMeanPt[i]->GetY()[0] << " +- " << gMeanPt[i]->GetErrorY(0) << endl;
        if (gMeanPt[i] == nullptr)
        {
            cout << "Error reading graph for " << particlesLatex[i] << endl;
            return;
        }
        gMeanPt[i]->SetMarkerColor(colors[i]);
        gMeanPt[i]->SetLineColor(colors[i]);
        gMeanPt[i]->SetMarkerStyle(markers[i]);
        gMeanPt[i]->SetMarkerSize(1.7);

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
            gMeanPtvsMassMesons->SetPointError(i, 0.02, 0.02, meanPtAt13TeV_err[i], meanPtAt13TeV_err[i]);
        }
        else
        {
            gMeanPtvsMassBaryons->SetPoint(i - 5, particleMass[i], meanPtAt13TeV[i]);
            gMeanPtvsMassBaryons->SetPointError(i - 5, 0.02, 0.02, meanPtAt13TeV_err[i], meanPtAt13TeV_err[i]);
        }
    }
    SetGraphStyleCommon(gMeanPtvsMassMesons);
    gMeanPtvsMassMesons->SetMarkerStyle(22);
    gMeanPtvsMassMesons->SetMarkerColor(kMagenta);
    gMeanPtvsMassMesons->SetLineColor(kMagenta);
    gMeanPtvsMassMesons->SetLineWidth(2);
    gMeanPtvsMassMesons->SetMarkerSize(1.7);
    gMeanPtvsMassMesons->SetFillStyle(0);
    SetGraphStyleCommon(gMeanPtvsMassBaryons);
    gMeanPtvsMassBaryons->SetMarkerStyle(34);
    gMeanPtvsMassBaryons->SetMarkerColor(kRed);
    gMeanPtvsMassBaryons->SetLineColor(kRed);
    gMeanPtvsMassBaryons->SetLineWidth(2);
    gMeanPtvsMassBaryons->SetMarkerSize(1.7);
    gMeanPtvsMassBaryons->SetFillStyle(0);

    TCanvas *cMeanPt = new TCanvas("cMeanPt", "Mean pT vs mass", 720, 720);
    SetCanvasStyle(cMeanPt, 0.14, 0.03, 0.05, 0.14);
    gMeanPtvsMassMesons->GetXaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    gMeanPtvsMassMesons->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/#it{c})");
    gMeanPtvsMassMesons->GetYaxis()->SetTitleOffset(1.3);
    gMeanPtvsMassMesons->SetMinimum(0.12);
    gMeanPtvsMassMesons->SetMaximum(3.59);
    gMeanPtvsMassMesons->GetXaxis()->SetLimits(0, 1.99);
    gMeanPtvsMassMesons->Draw("A 2"); // draw boxes for systematic errors
    TGraphAsymmErrors *gMeanPtvsMassMesonsSys = (TGraphAsymmErrors *)gMeanPtvsMassMesons->Clone("gMeanPtvsMassMesonsSys");
    gMeanPtvsMassMesonsSys->SetLineWidth(0); // to avoid showing error bars
    gMeanPtvsMassMesonsSys->Draw("P SAME");
    gMeanPtvsMassBaryons->Draw("2 SAME");
    TGraphAsymmErrors *gMeanPtvsMassBaryonsSys = (TGraphAsymmErrors *)gMeanPtvsMassBaryons->Clone("gMeanPtvsMassBaryonsSys");
    gMeanPtvsMassBaryonsSys->SetLineWidth(0); // to avoid showing error bars
    gMeanPtvsMassBaryonsSys->Draw("P SAME");

    TF1 *pol1_meson = new TF1("pol1_meson", "pol1", 0.1, 1.89);
    pol1_meson->SetLineColor(kYellow + 2);
    pol1_meson->SetLineStyle(2);
    gMeanPtvsMassMesons->Fit(pol1_meson, "R");

    TF1 *pol1_baryon = new TF1("pol1_baryon", "pol1", 0.1, 1.89);
    pol1_baryon->SetLineColor(kYellow + 2);
    pol1_baryon->SetLineStyle(2);
    gMeanPtvsMassBaryons->Fit(pol1_baryon, "R");

    // Draw the last marker (f2(1525)) and its error bar
    double f2_mass = 1.5173;
    double f2_meanpt = hout2->GetBinContent(5);
    double f2_Staterr = hout2->GetBinContent(6);
    double f2_SysErrLow = hout2->GetBinContent(7);
    double f2_SysErrHigh = hout2->GetBinContent(8);
    // cout<<"Sys Low "<<f2_SysErrLow<<" Sys High "<<f2_SysErrHigh<<endl;

    int f2_marker = 20;        // choose a unique marker style for f2(1525)
    int f2_color = kGreen + 1; // choose a unique color for f2(1525)
    TGraphErrors *graph_f2 = new TGraphErrors(1);
    graph_f2->SetPoint(0, f2_mass, f2_meanpt);
    graph_f2->SetPointError(0, 0, f2_Staterr);
    graph_f2->SetMarkerStyle(f2_marker);
    graph_f2->SetMarkerColor(f2_color);
    graph_f2->SetMarkerSize(1.7);
    graph_f2->SetLineColor(f2_color);
    graph_f2->SetLineWidth(2);
    graph_f2->Draw("P SAME");

    // Draw systematic error as boxes using TGraphAsymmErrors
    TGraphAsymmErrors *graph_f2_sys = new TGraphAsymmErrors(1);
    graph_f2_sys->SetPoint(0, f2_mass, f2_meanpt);
    graph_f2_sys->SetPointError(0, 0.02, 0.02, f2_SysErrLow, f2_SysErrHigh);
    graph_f2_sys->SetLineColor(f2_color);
    graph_f2_sys->SetLineWidth(2);
    graph_f2_sys->SetFillStyle(0);
    graph_f2_sys->Draw("2 SAME");

    // Draw the last marker (f0(1710)) and its error bar
    double f0_mass = 1.710;
    double f0_meanpt = hout->GetBinContent(5);
    double f0_Staterr = hout->GetBinContent(6);
    double f0_SysErrLow = hout->GetBinContent(7);
    double f0_SysErrHigh = hout->GetBinContent(8);
    int f0_marker = 21;   // choose a unique marker style for f0(1710)
    int f0_color = kBlue; // choose a unique color for f0(1710)
    TGraphErrors *graph_f0 = new TGraphErrors(1);
    graph_f0->SetPoint(0, f0_mass, f0_meanpt);
    graph_f0->SetPointError(0, 0, f0_Staterr);
    graph_f0->SetMarkerStyle(f0_marker);
    graph_f0->SetMarkerColor(f0_color);
    graph_f0->SetMarkerSize(1.7);
    graph_f0->SetLineColor(f0_color);
    graph_f0->SetLineWidth(2);
    graph_f0->Draw("P SAME");

    TGraphAsymmErrors *graph_f0_sys = new TGraphAsymmErrors(1);
    graph_f0_sys->SetPoint(0, f0_mass, f0_meanpt);
    graph_f0_sys->SetPointError(0, 0.02, 0.02, f0_SysErrLow, f0_SysErrHigh);
    graph_f0_sys->SetLineColor(f0_color);
    graph_f0_sys->SetLineWidth(2);
    graph_f0_sys->SetFillStyle(0);
    graph_f0_sys->Draw("2 SAME");

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
    latex.DrawLatex(1.5173, hout2->GetBinContent(5) + 0.22, "f'_{2}(1525)");
    latex.DrawLatex(1.710, hout->GetBinContent(5) + 0.2, "f_{0}(1710)");

    TLegend *legend4 = new TLegend(0.17, 0.65, 0.63, 0.92);
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->SetTextSize(0.035);
    legend4->AddEntry(gMeanPtvsMassMesons, "Mesons (13 TeV)", "p");
    legend4->AddEntry(gMeanPtvsMassBaryons, "Baryons (13 TeV)", "p");
    legend4->AddEntry(graph_f2, "f'_{2}(1525) (13.6 TeV)", "p");
    legend4->AddEntry(graph_f0, "f_{0}(1710) (13.6 TeV)", "p");
    legend4->AddEntry(pol1_meson, "Pol 1", "l");
    legend4->Draw();

    TLegend *legend5 = new TLegend(0.52, 0.78, 0.87, 0.92);
    legend5->SetBorderSize(0);
    legend5->SetFillStyle(0);
    legend5->SetTextSize(0.035);
    legend5->AddEntry((TObject *)0, "ALICE", "");
    legend5->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legend5->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    legend5->Draw();
    // cMeanPt->SaveAs((savePath + "/plots/MeanPt_vs_Mass_reweighted.png").c_str());

    /*
    // Reading by from the root file
    float ptBins[] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0}; // 2022 dataset
    int nBins = sizeof(ptBins) / sizeof(ptBins[0]) - 1;

    TH1F *hEfficiencyReweighted = new TH1F("hEfficiencyReweighted", "Reweighted Efficiency vs pT", nBins, ptBins);
    TH1F *hEfficiencyUnweighted = new TH1F("hEfficiencyUnweighted", "Unweighted Efficiency vs pT", nBins, ptBins);
    TH1F *hEfficiencyReweighted2 = new TH1F("hEfficiencyReweighted2", "Reweighted Efficiency vs pT for f2(1525)", nBins, ptBins);
    TH1F *hEfficiencyUnweighted2 = new TH1F("hEfficiencyUnweighted2", "Unweighted Efficiency vs pT for f2(1525)", nBins, ptBins);

    for (int i = 0; i < nBins; i++)
    {
        float lowpt = ptBins[i] + 0.01;
        float highpt = ptBins[i + 1] - 0.01; //
        // For reweighted efficiency
        int lowptBin = hGenReweighted->GetXaxis()->FindBin(lowpt);
        int highptBin = hGenReweighted->GetXaxis()->FindBin(highpt);
        Double_t recErrW = 0, genErrW = 0;
        Double_t recYieldReweighted = hRecReweighted->IntegralAndError(lowptBin, highptBin, recErrW);
        Double_t genYieldReweighted = hGenReweighted->IntegralAndError(lowptBin, highptBin, genErrW);
        if (genYieldReweighted > 0 && recYieldReweighted > 0)
        {
            Double_t effRW = recYieldReweighted / genYieldReweighted;
            Double_t errRW = effRW * TMath::Sqrt(TMath::Power(recErrW / recYieldReweighted, 2) + TMath::Power(genErrW / genYieldReweighted, 2));
            hEfficiencyReweighted->SetBinContent(i + 1, effRW);
            hEfficiencyReweighted->SetBinError(i + 1, errRW);
        }

        // Unweighted efficiency from pre-scaled count clones; use binomial approximation
        Double_t recErrU = 0, genErrU = 0;
        Double_t recCounts = hRecUnweightedCounts ? hRecUnweightedCounts->IntegralAndError(lowptBin, highptBin, recErrU) : 0.0;
        Double_t genCounts = hGenUnweightedCounts ? hGenUnweightedCounts->IntegralAndError(lowptBin, highptBin, genErrU) : 0.0;
        if (genCounts > 0 && recCounts >= 0)
        {
            Double_t effU = recCounts / genCounts;
            // Binomial error with genCounts as trials; guard numerical range
            Double_t errU = TMath::Sqrt(TMath::Max(0.0, effU * (1.0 - effU) / genCounts));
            hEfficiencyUnweighted->SetBinContent(i + 1, effU);
            hEfficiencyUnweighted->SetBinError(i + 1, errU);
        }

        // Now calculate for f2'(1525) as well
        // For reweighted efficiency
        lowptBin = hGenReweighted2->GetXaxis()->FindBin(lowpt);
        highptBin = hGenReweighted2->GetXaxis()->FindBin(highpt);
        Double_t recErrW2 = 0, genErrW2 = 0;
        Double_t recYieldReweighted2 = hRecReweighted2->IntegralAndError(lowptBin, highptBin, recErrW2);
        Double_t genYieldReweighted2 = hGenReweighted2->IntegralAndError(lowptBin, highptBin, genErrW2);
        if (genYieldReweighted2 > 0 && recYieldReweighted2 > 0)
        {
            Double_t effRW2 = recYieldReweighted2 / genYieldReweighted2;
            Double_t errRW2 = effRW2 * TMath::Sqrt(TMath::Power(recErrW2 / recYieldReweighted2, 2) + TMath::Power(genErrW2 / genYieldReweighted2, 2));
            hEfficiencyReweighted2->SetBinContent(i + 1, effRW2);
            hEfficiencyReweighted2->SetBinError(i + 1, errRW2);
        }
        // Unweighted efficiency from pre-scaled count clones; use binomial approximation
        Double_t recErrU2 = 0, genErrU2 = 0;
        Double_t recCounts2 = hRecUnweightedCounts2 ? hRecUnweightedCounts2->IntegralAndError(lowptBin, highptBin, recErrU2) : 0.0;
        Double_t genCounts2 = hGenUnweightedCounts2 ? hGenUnweightedCounts2->IntegralAndError(lowptBin, highptBin, genErrU2) : 0.0;
        if (genCounts2 > 0 && recCounts2 >= 0)
        {
            Double_t effU2 = recCounts2 / genCounts2;
            // Binomial error with genCounts as trials; guard numerical range
            Double_t errU2 = TMath::Sqrt(TMath::Max(0.0, effU2 * (1.0 - effU2) / genCounts2));
            hEfficiencyUnweighted2->SetBinContent(i + 1, effU2);
            hEfficiencyUnweighted2->SetBinError(i + 1, errU2);
        }
    }

    fOutput->cd();
    TCanvas *cNewEfficiency = new TCanvas("cNewEfficiency", "New Efficiency Comparison for f0(1710)", 720, 720);
    SetCanvasStyle(cNewEfficiency, 0.18, 0.03, 0.05, 0.14);
    double pad1Size, pad2Size;
    canvas_style(cNewEfficiency, pad1Size, pad2Size);
    cNewEfficiency->cd(1);
    SetHistoQA(hEfficiencyReweighted);
    hEfficiencyReweighted->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hEfficiencyReweighted->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hEfficiencyReweighted->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hEfficiencyReweighted->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hEfficiencyReweighted->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hEfficiencyReweighted->GetYaxis()->SetTitle("Acceptance x Efficiency");
    hEfficiencyReweighted->GetYaxis()->SetTitleOffset(1.2);
    hEfficiencyReweighted->SetMaximum(hEfficiencyReweighted->GetMaximum() * 1.8);
    hEfficiencyReweighted->SetMarkerStyle(20);
    hEfficiencyReweighted->SetMarkerColor(kRed);
    hEfficiencyReweighted->SetLineColor(kRed);
    hEfficiencyReweighted->SetMarkerSize(1.5);
    hEfficiencyReweighted->Draw("pe");
    hEfficiencyReweighted->Write("Eff_f0Reweighted");
    SetHistoQA(hEfficiencyUnweighted);
    hEfficiencyUnweighted->SetMarkerStyle(20);
    hEfficiencyUnweighted->SetMarkerColor(kBlue);
    hEfficiencyUnweighted->SetLineColor(kBlue);
    hEfficiencyUnweighted->SetMarkerSize(1.5);
    hEfficiencyUnweighted->Draw("pe same");
    TLegend *legNewEfficiency = new TLegend(0.4, 0.67, 0.9, 0.93);
    legNewEfficiency->SetBorderSize(0);
    legNewEfficiency->SetFillStyle(0);
    legNewEfficiency->SetTextSize(0.05);
    legNewEfficiency->AddEntry((TObject *)nullptr, "pp #sqrt{#it{s}} = 13.6 TeV", "");
    legNewEfficiency->AddEntry((TObject *)nullptr, "f_{0}(1710)", "");
    legNewEfficiency->AddEntry(hEfficiencyUnweighted, "Unweighted Efficiency", "pe");
    legNewEfficiency->AddEntry(hEfficiencyReweighted, "Reweighted Efficiency", "pe");
    legNewEfficiency->Draw();

    cNewEfficiency->cd(2);
    TH1F *hEfficiencyRatio = (TH1F *)hEfficiencyReweighted->Clone("hEfficiencyRatio");
    hEfficiencyRatio->Divide(hEfficiencyUnweighted);
    SetHistoQA(hEfficiencyRatio);
    hEfficiencyRatio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hEfficiencyRatio->GetYaxis()->SetTitleSize(0.02 / pad2Size);
    hEfficiencyRatio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hEfficiencyRatio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hEfficiencyRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hEfficiencyRatio->GetYaxis()->SetTitle("Reweighted / Unweighted");
    hEfficiencyRatio->GetYaxis()->SetTitleOffset(0.9);
    hEfficiencyRatio->SetMaximum(1.3);
    hEfficiencyRatio->SetMinimum(0.7);
    hEfficiencyRatio->SetMarkerStyle(20);
    hEfficiencyRatio->SetMarkerColor(kBlack);
    hEfficiencyRatio->SetLineColor(kBlack);
    hEfficiencyRatio->SetMarkerSize(1.5);
    hEfficiencyRatio->GetYaxis()->SetNdivisions(505);
    hEfficiencyRatio->Draw("lp");
    TLine *lineUnity = new TLine(hEfficiencyRatio->GetXaxis()->GetXmin(), 1.0, hEfficiencyRatio->GetXaxis()->GetXmax(), 1.0);
    lineUnity->SetLineColor(kRed);
    lineUnity->SetLineStyle(2);
    lineUnity->Draw("same");

    cNewEfficiency->SaveAs((savePath + "/plots/EfficiencyComparisonf0.png").c_str());

    TCanvas *cNewEfficiency2 = new TCanvas("cNewEfficiency2", "New Efficiency Comparison for f2(1525)", 720, 720);
    SetCanvasStyle(cNewEfficiency2, 0.18, 0.03, 0.05, 0.14);
    canvas_style(cNewEfficiency2, pad1Size, pad2Size);
    cNewEfficiency2->cd(1);
    SetHistoQA(hEfficiencyReweighted2);
    hEfficiencyReweighted2->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hEfficiencyReweighted2->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hEfficiencyReweighted2->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hEfficiencyReweighted2->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hEfficiencyReweighted2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hEfficiencyReweighted2->GetYaxis()->SetTitle("Acceptance x Efficiency");
    hEfficiencyReweighted2->GetYaxis()->SetTitleOffset(1.2);
    hEfficiencyReweighted2->SetMaximum(hEfficiencyReweighted2->GetMaximum() * 1.8);
    hEfficiencyReweighted2->SetMarkerStyle(20);
    hEfficiencyReweighted2->SetMarkerColor(kRed);
    hEfficiencyReweighted2->SetLineColor(kRed);
    hEfficiencyReweighted2->SetMarkerSize(1.5);
    hEfficiencyReweighted2->Write("Eff_f2Reweighted");
    hEfficiencyReweighted2->Draw("pe");
    SetHistoQA(hEfficiencyUnweighted2);
    hEfficiencyUnweighted2->SetMarkerStyle(20);
    hEfficiencyUnweighted2->SetMarkerColor(kBlue);
    hEfficiencyUnweighted2->SetLineColor(kBlue);
    hEfficiencyUnweighted2->SetMarkerSize(1.5);
    hEfficiencyUnweighted2->Draw("pe same");

    TLegend *legNewEfficiency2 = new TLegend(0.4, 0.67, 0.9, 0.93);
    legNewEfficiency2->SetBorderSize(0);
    legNewEfficiency2->SetFillStyle(0);
    legNewEfficiency2->SetTextSize(0.05);
    legNewEfficiency2->AddEntry((TObject *)nullptr, "pp #sqrt{#it{s}} = 13.6 TeV", "");
    legNewEfficiency2->AddEntry((TObject *)nullptr, "f'_{2}(1525)", "");
    legNewEfficiency2->AddEntry(hEfficiencyUnweighted2, "Unweighted Efficiency", "pe");
    legNewEfficiency2->AddEntry(hEfficiencyReweighted2, "Reweighted Efficiency", "pe");
    legNewEfficiency2->Draw();

    cNewEfficiency2->cd(2);
    TH1F *hEfficiencyRatio2 = (TH1F *)hEfficiencyReweighted2->Clone("hEfficiencyRatio2");
    hEfficiencyRatio2->Divide(hEfficiencyUnweighted2);
    SetHistoQA(hEfficiencyRatio2);
    hEfficiencyRatio2->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hEfficiencyRatio2->GetYaxis()->SetTitleSize(0.02 / pad2Size);
    hEfficiencyRatio2->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hEfficiencyRatio2->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hEfficiencyRatio2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hEfficiencyRatio2->GetYaxis()->SetTitle("Reweighted / Unweighted");
    hEfficiencyRatio2->GetYaxis()->SetTitleOffset(0.9);
    hEfficiencyRatio2->SetMaximum(1.3);
    hEfficiencyRatio2->SetMinimum(0.7);
    hEfficiencyRatio2->SetMarkerStyle(20);
    hEfficiencyRatio2->SetMarkerColor(kBlack);
    hEfficiencyRatio2->SetLineColor(kBlack);
    hEfficiencyRatio2->SetMarkerSize(1.5);
    hEfficiencyRatio2->GetYaxis()->SetNdivisions(505);
    hEfficiencyRatio2->Draw("p");
    TLine *lineUnity2 = new TLine(hEfficiencyRatio2->GetXaxis()->GetXmin(), 1.0, hEfficiencyRatio2->GetXaxis()->GetXmax(), 1.0);
    lineUnity2->SetLineColor(kRed);
    lineUnity2->SetLineStyle(2);
    lineUnity2->Draw("same");
    cNewEfficiency2->SaveAs((savePath + "/plots/EfficiencyComparisonf2.png").c_str());
    if (!otherQAPlots)
    {
        cNewEfficiency->Close();
        cNewEfficiency2->Close();
    }

    TCanvas *cReweightFactor = new TCanvas("cReweightFactor", "Reweighting Factor for f0(1710)", 720, 720);
    SetCanvasStyle(cReweightFactor, 0.18, 0.03, 0.05, 0.14);
    TH1F *hReweightFactor = (TH1F *)hGenReweighted->Clone("hYield1710Corrected_correction_i3");
    SetHistoQA(hReweightFactor);
    hReweightFactor->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hReweightFactor->GetYaxis()->SetTitle("Reweight Factor");
    hReweightFactor->GetYaxis()->SetTitleOffset(1.5);
    hReweightFactor->SetMaximum(1.5);
    hReweightFactor->SetMinimum(0.3);
    hReweightFactor->Draw("HIST");
    if (!otherQAPlots)
        cReweightFactor->Close();

    TFile *file2 = new TFile((savePath + "/spectra_.root").c_str(), "read");
    TH1F *hYield1710Corrected = (TH1F *)file2->Get("hYield1710Corrected");
    if (hYield1710Corrected == nullptr)
    {
        cout << "Error reading corrected yield histogram from file" << endl;
        return;
    }
    TCanvas *cYieldCompare = new TCanvas("cYieldCompare", "Corrected Yield Comparison for f0(1710)", 720, 720);
    SetCanvasStyle(cYieldCompare, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hYield1710Corrected);
    gPad->SetLogy();
    hYield1710Corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hYield1710Corrected->GetYaxis()->SetTitle("1/N_{ev} * d^{2}N/(d#it{p}_{T}dy) (GeV/#it{c})^{-1}");
    hYield1710Corrected->GetYaxis()->SetTitleOffset(1.5);
    hYield1710Corrected->SetMaximum(hYield1710Corrected->GetMaximum() * 1.5);
    hYield1710Corrected->SetMarkerStyle(20);
    hYield1710Corrected->SetMarkerColor(kRed);
    hYield1710Corrected->SetLineColor(kRed);
    hYield1710Corrected->GetXaxis()->SetRangeUser(0.0, 15.5);
    hYield1710Corrected->SetMinimum(1e-6);
    hYield1710Corrected->SetMaximum(hYield1710Corrected->GetMaximum() * 12);
    hYield1710Corrected->SetMarkerSize(1.5);
    hYield1710Corrected->Draw("pe");
    hYieldReweighted->Draw("pe same");

    TLegend *legYieldCompare = new TLegend(0.4, 0.73, 0.9, 0.93);
    legYieldCompare->SetBorderSize(0);
    legYieldCompare->SetFillStyle(0);
    legYieldCompare->SetTextSize(0.035);
    legYieldCompare->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    legYieldCompare->AddEntry(hYield1710Corrected, "Corrected spectra", "pe");
    legYieldCompare->AddEntry(hYieldReweighted, "Reweighted corrected spectra", "pe");
    legYieldCompare->Draw();
    cYieldCompare->SaveAs((savePath + "/plots/CorrectedYieldComparisonf0.png").c_str());
    if (!otherQAPlots)
        cYieldCompare->Close();
    */

    cout << "Finished processing plotting the reweighted MC spectra" << endl;
}

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
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.14);
    pad2->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.04);
}