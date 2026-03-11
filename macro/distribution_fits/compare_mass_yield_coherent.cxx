#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"
#include "../spectra/YieldMean2.C"

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

void compare_mass_yield_coherent()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra/";
    TFile *fDefault = new TFile((path + "spectra_Default4.root").c_str(), "READ");
    TFile *fReweightedf0 = new TFile((path + "ReweighFacf0_Default4.root").c_str(), "READ");
    TFile *fReweightedf2 = new TFile((path + "ReweighFacf2_Default4.root").c_str(), "READ");
    TFile *fCoherent = new TFile((path + "spectra_coherent.root").c_str(), "READ");
    // TFile *fCoherent = new TFile((path + "spectra_CoherentSumPhi0.root").c_str(), "READ");

    // Find the highest available index for reweighted histograms
    int maxIndexReweightedf0 = FindHighestIndex(fReweightedf0, "Genf17102_proj_1_");
    int maxIndexReweightedf2 = FindHighestIndex(fReweightedf2, "Genf17102_proj_1_");
    if (maxIndexReweightedf0 == -1 || maxIndexReweightedf2 == -1)
    {
        cout << "Error: No reweighted histogram with pattern Genf17102_proj_1_i* found in file" << endl;
        return;
    }
    cout << "Using index i" << maxIndexReweightedf0 << " for reweighted histograms" << endl;

    string indexStr = "i" + to_string(maxIndexReweightedf0);
    string indexStr2 = "i" + to_string(maxIndexReweightedf2);
    TH1F *hYieldReweightedf0 = (TH1F *)fReweightedf0->Get(Form("hYield1710Corrected_%s", indexStr.c_str()));
    TH1F *hYieldReweightedf2 = (TH1F *)fReweightedf2->Get(Form("hYield1525Corrected_%s", indexStr2.c_str()));

    int totalBins = hYieldReweightedf0->GetNbinsX();
    // cout << "Total bins in reweighted histogram: " << totalBins << endl;
    // hYieldReweightedf0->SetBinContent(totalBins -1, hYieldReweightedf0->GetBinContent(totalBins -1) * 1.15);
    // hYieldReweightedf2->SetBinContent(totalBins -1, hYieldReweightedf2->GetBinContent(totalBins -1) * 1.40);

    if (fDefault->IsZombie() || fCoherent->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hMass_f1525_Default = (TH1F *)fDefault->Get("hMass_1525");
    TH1F *hMass_f1710_Default = (TH1F *)fDefault->Get("hMass_1710");
    // TH1F *hMass_f1525_Coherent = (TH1F *)fCoherent->Get("hMass_1525");
    // TH1F *hMass_f1710_Coherent = (TH1F *)fCoherent->Get("hMass_1710");
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
    // cMassf0->SaveAs((path + "MassCompare_f1710_CoherentPhi0.png").c_str());

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
    // cMassf2->SaveAs((path + "MassComparison_f1525.png").c_str());
    // cMassf2->SaveAs((path + "MassCompare_f1525_CoherentPhi0.png").c_str());

    // Compare the spectra for default and variation case for f2(1525) and f0(1710)
    TH1F *SpectraCoherent_f1525 = (TH1F *)fCoherent->Get("hYield1525Corrected");
    TH1F *SpectraCoherent_f1710 = (TH1F *)fCoherent->Get("hYield1710Corrected");

    // TFile *fSigEventLossOld = new TFile("Loss_phi_mult0-100.root", "READ");
    // TFile *fSigEventLossNew = new TFile("../SignalLossPhiINEL.root", "READ");

    // TH1F *hEvbySigLossOld = (TH1F *)fSigEventLossOld->Get("hSignalLoss");
    // TH1F *hEvbySigLossNew = (TH1F *)fSigEventLossNew->Get("hSignalLoss");

    // TH1F *hNewbyOld = (TH1F *)hEvbySigLossNew->Clone("hNewbyOld");
    // hNewbyOld->Divide(hEvbySigLossOld);

    // for (int i = 1; i <= hNewbyOld->GetNbinsX(); i++)
    // {
    //     double binContentf2 = SpectraCoherent_f1525->GetBinContent(i + 1);
    //     double binContentf0 = SpectraCoherent_f1710->GetBinContent(i + 1);
    //     double ratio = hNewbyOld->GetBinContent(i);
    //     double newBinContentf2 = binContentf2 * ratio;
    //     double newBinContentf0 = binContentf0 * ratio;
    //     SpectraCoherent_f1525->SetBinContent(i + 1, newBinContentf2);
    //     SpectraCoherent_f1710->SetBinContent(i + 1, newBinContentf0);
    // }

    TFile *fSys = new TFile((path + "SystematicPlots/SystematicUncertainties.root").c_str(), "read");
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

    TH1F *hf21 = (TH1F *)SpectraCoherent_f1525->Clone("hf21");
    TH1F *hf22 = (TH1F *)SpectraCoherent_f1525->Clone("hf22");
    TH1F *hf01 = (TH1F *)SpectraCoherent_f1710->Clone("hf01");
    TH1F *hf02 = (TH1F *)SpectraCoherent_f1710->Clone("hf02");
    hf21->Sumw2();
    hf22->Sumw2();
    hf01->Sumw2();
    hf02->Sumw2();

    for (int i = 1; i <= hf21->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double sys1710 = hYieldSysf0->GetBinContent(i);
        double sys1525 = hYieldSysf2->GetBinContent(i);
        double yield1710 = hf01->GetBinContent(i + 1);
        double yield1525 = hf21->GetBinContent(i + 1);

        hf02->SetBinContent(i + 1, yield1710);
        hf02->SetBinError(i + 1, sys1710 * yield1710);
        hf22->SetBinContent(i + 1, yield1525);
        hf22->SetBinError(i + 1, sys1525 * yield1525);
    }

    Double_t min = 0.0;
    Double_t max = 15.0;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    // Option_t *opt = "REBMS0+";
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 1.0;
    Double_t maxfit = 15.0;

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

    TH1 *houtf0 = YieldMean(hf01, hf02, fitFcnf0, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *houtf2 = YieldMean(hf21, hf22, fitFcnf2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    // Now compare the <pT> values for f0(1710) and f2(1525) between the two fits
    TGraphErrors *gMeanpTBothDefault = (TGraphErrors *)fDefault->Get("gMeanPt_f0f2");
    TGraphErrors *gMeanpTBothCoherent = (TGraphErrors *)fCoherent->Get("gMeanPt_f0f2");
    if (gMeanpTBothDefault == nullptr || gMeanpTBothCoherent == nullptr)
    {
        cout << "Error: gMeanPt_f0f2 not found in one of the files!" << endl;
        return;
    }
    // double f1525MeanPt_Default = gMeanpTBothDefault->GetY()[0];
    // double f1525MeanPt_Default_err = gMeanpTBothDefault->GetErrorY(0);
    // double f1710MeanPt_Default = gMeanpTBothDefault->GetY()[1];
    // double f1710MeanPt_Default_err = gMeanpTBothDefault->GetErrorY(1);

    // Taking values after the reweighting for f0 and f2
    double f1525MeanPt_Default = 1.6474;
    double f1525MeanPt_Default_err = 0.0446258;
    double f1710MeanPt_Default = 2.35281;
    double f1710MeanPt_Default_err = 0.0938778;

    // double f1525MeanPt_Coherent = gMeanpTBothCoherent->GetY()[0];
    // double f1525MeanPt_Coherent_err = gMeanpTBothCoherent->GetErrorY(0);
    // double f1710MeanPt_Coherent = gMeanpTBothCoherent->GetY()[1];
    // double f1710MeanPt_Coherent_err = gMeanpTBothCoherent->GetErrorY(1);

    double f1525MeanPt_Coherent = houtf2->GetBinContent(5);
    double f1525MeanPt_Coherent_err = houtf2->GetBinContent(6);
    double f1710MeanPt_Coherent = houtf0->GetBinContent(5);
    double f1710MeanPt_Coherent_err = houtf0->GetBinContent(6);

    cout << "f2(1525) <pT> Default: " << f1525MeanPt_Default << " +- " << f1525MeanPt_Default_err << endl;
    cout << "f2(1525) <pT> Coherent: " << f1525MeanPt_Coherent << " +- " << f1525MeanPt_Coherent_err << endl;
    cout << "f0(1710) <pT> Default: " << f1710MeanPt_Default << " +- " << f1710MeanPt_Default_err << endl;
    cout << "f0(1710) <pT> Coherent: " << f1710MeanPt_Coherent << " +- " << f1710MeanPt_Coherent_err << endl;
    cout << "Relative change in <pT> for f2(1525): " << (f1525MeanPt_Coherent - f1525MeanPt_Default) / f1525MeanPt_Default * 100 << " %" << endl;
    cout << "Relative change in <pT> for f0(1710): " << (f1710MeanPt_Coherent - f1710MeanPt_Default) / f1710MeanPt_Default * 100 << " %" << endl;

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
    // gMeanPtvsMassMesons->SetMaximum(4.59);
    gMeanPtvsMassMesons->SetMaximum(4.09);
    gMeanPtvsMassMesons->GetXaxis()->SetLimits(0, 1.99);
    gMeanPtvsMassMesons->Draw("A 2"); // draw boxes for systematic errors
    TGraphAsymmErrors *gMeanPtvsMassMesonsSys = (TGraphAsymmErrors *)gMeanPtvsMassMesons->Clone("gMeanPtvsMassMesonsSys");
    gMeanPtvsMassMesonsSys->SetLineWidth(0); // to avoid showing error bars
    gMeanPtvsMassMesonsSys->Draw("P SAME");
    gMeanPtvsMassBaryons->Draw("2 SAME");
    TGraphAsymmErrors *gMeanPtvsMassBaryonsSys = (TGraphAsymmErrors *)gMeanPtvsMassBaryons->Clone("gMeanPtvsMassBaryonsSys");
    gMeanPtvsMassBaryonsSys->SetLineWidth(0); // to avoid showing error bars
    gMeanPtvsMassBaryonsSys->Draw("P SAME");

    TF1 *pol1_meson = new TF1("pol1_meson", "pol1", 0.1, 1.8);
    pol1_meson->SetLineColor(kYellow + 2);
    pol1_meson->SetLineStyle(2);
    gMeanPtvsMassMesons->Fit(pol1_meson, "R");

    TF1 *pol1_baryon = new TF1("pol1_baryon", "pol1", 0.1, 1.8);
    pol1_baryon->SetLineColor(kYellow + 2);
    pol1_baryon->SetLineStyle(2);
    gMeanPtvsMassBaryons->Fit(pol1_baryon, "R");

    // Draw the f2(1525) point using graph + systematic box style
    double f2_mass = 1.5173;
    double f2_meanpt = f1525MeanPt_Default;
    double f2_err = f1525MeanPt_Default_err;
    double f2_meanpt2 = f1525MeanPt_Coherent;
    double f2_err2 = f1525MeanPt_Coherent_err;
    int f2_marker = 20;        // choose a unique marker style for f2(1525)
    int f2_color = kGreen + 1; // choose a unique color for f2(1525)
    int f2_color2 = kCyan + 2;
    int f2_marker2 = 21;
    TGraphErrors *graph_f2 = new TGraphErrors(1);
    graph_f2->SetPoint(0, f2_mass, f2_meanpt);
    graph_f2->SetPointError(0, 0, f2_err);
    graph_f2->SetMarkerStyle(f2_marker);
    graph_f2->SetMarkerColor(f2_color);
    graph_f2->SetMarkerSize(1.5);
    graph_f2->SetLineColor(f2_color);
    graph_f2->SetLineWidth(2);
    graph_f2->Draw("P SAME");

    TGraphAsymmErrors *graph_f2_sys = new TGraphAsymmErrors(1);
    graph_f2_sys->SetPoint(0, f2_mass, f2_meanpt);
    graph_f2_sys->SetPointError(0, 0.02, 0.02, f2_err, f2_err);
    graph_f2_sys->SetLineColor(f2_color);
    graph_f2_sys->SetLineWidth(2);
    graph_f2_sys->SetFillStyle(0);
    // graph_f2_sys->Draw("2 SAME");

    // Draw the f2(1525) coherent point
    TGraphErrors *graph_f2_2 = new TGraphErrors(1);
    graph_f2_2->SetPoint(0, f2_mass, f2_meanpt2);
    graph_f2_2->SetPointError(0, 0, f2_err2);
    graph_f2_2->SetMarkerStyle(25);
    graph_f2_2->SetMarkerColor(kCyan + 1);
    graph_f2_2->SetMarkerSize(1.5);
    graph_f2_2->SetLineColor(kCyan + 1);
    graph_f2_2->SetLineWidth(2);
    graph_f2_2->Draw("P SAME");

    // Draw the f0(1710) default/coherent points using graph + box style
    double f0_mass = 1.710;
    double f0_meanpt = f1710MeanPt_Default;
    double f0_err = f1710MeanPt_Default_err;
    double f0_meanpt2 = f1710MeanPt_Coherent;
    double f0_err2 = f1710MeanPt_Coherent_err;
    int f0_marker = 21; // choose a unique marker style for f0(1710)
    int f0_marker2 = 20;
    int f0_color = kBlue; // choose a unique color for f0(1710)
    int f0_color2 = kBrown;
    TGraphErrors *graph_f0 = new TGraphErrors(1);
    graph_f0->SetPoint(0, f0_mass, f0_meanpt);
    graph_f0->SetPointError(0, 0, f0_err);
    graph_f0->SetMarkerStyle(f0_marker);
    graph_f0->SetMarkerColor(f0_color);
    graph_f0->SetMarkerSize(1.5);
    graph_f0->SetLineColor(f0_color);
    graph_f0->SetLineWidth(2);
    graph_f0->Draw("P SAME");

    TGraphAsymmErrors *graph_f0_sys = new TGraphAsymmErrors(1);
    graph_f0_sys->SetPoint(0, f0_mass, f0_meanpt);
    graph_f0_sys->SetPointError(0, 0.02, 0.02, f0_err, f0_err);
    graph_f0_sys->SetLineColor(f0_color);
    graph_f0_sys->SetLineWidth(2);
    graph_f0_sys->SetFillStyle(0);
    // graph_f0_sys->Draw("2 SAME");

    TGraphErrors *graph_f0_2 = new TGraphErrors(1);
    graph_f0_2->SetPoint(0, f0_mass, f0_meanpt2);
    graph_f0_2->SetPointError(0, 0, f0_err2);
    graph_f0_2->SetMarkerStyle(f0_marker2);
    graph_f0_2->SetMarkerColor(f0_color2);
    graph_f0_2->SetMarkerSize(1.5);
    graph_f0_2->SetLineColor(f0_color2);
    graph_f0_2->SetLineWidth(2);
    graph_f0_2->Draw("P SAME");

    TGraphAsymmErrors *graph_f0_2_sys = new TGraphAsymmErrors(1);
    graph_f0_2_sys->SetPoint(0, f0_mass, f0_meanpt2);
    graph_f0_2_sys->SetPointError(0, 0.02, 0.02, f0_err2, f0_err2);
    graph_f0_2_sys->SetLineColor(f0_color2);
    graph_f0_2_sys->SetLineWidth(2);
    graph_f0_2_sys->SetFillStyle(0);
    // graph_f0_2_sys->Draw("2 SAME");
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

    TLegend *legend4 = new TLegend(0.18, 0.6, 0.65, 0.92);
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->SetTextSize(0.035);
    legend4->AddEntry(gMeanPtvsMassMesons, "Mesons (13 TeV)", "p");
    legend4->AddEntry(gMeanPtvsMassBaryons, "Baryons (13 TeV)", "p");
    legend4->AddEntry(graph_f2, "f'_{2}(1525) (Default)", "p");
    legend4->AddEntry(graph_f2_2, "f'_{2}(1525) (Coherent)", "p");
    legend4->AddEntry(graph_f0, "f_{0}(1710) (Default)", "p");
    legend4->AddEntry(graph_f0_2, "f_{0}(1710) (Coherent)", "p");
    legend4->AddEntry(pol1_meson, "Pol 1", "l");
    legend4->Draw();
    // cMeanPt->SaveAs((path + "plots/MeanPt_vs_Mass_CoherentCompare.png").c_str());
    // cMeanPt->SaveAs((path + "plots/MeanPt_vs_Mass_CoherentComparePhi0.png").c_str());

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
    hf01->SetMinimum(2e-8);
    hf01->SetMarkerColor(kBlue);
    hf01->SetLineColor(kBlue);
    hf01->SetStats(1);
    hf01->Draw("pe");
    hf02->SetMarkerColor(kBlue);
    hf02->SetLineColor(kBlue);
    hf02->SetFillStyle(0);
    hf02->SetMinimum(2e-8);
    hf02->Draw("e2 same");
    fitFcnf0->SetLineWidth(2);
    fitFcnf0->SetLineStyle(2);
    fitFcnf0->SetLineColor(kRed);
    fitFcnf0->Draw("l same");

    // Create legend with physics information
    TLegend *leg = new TLegend(0.19, 0.2, 0.4, 0.45);
    // leg->AddEntry((TObject *)0, "ALICE", "");
    leg->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    leg->AddEntry(hf01, "f_{0}(1710) spectra", "pe");
    leg->AddEntry(fitFcnf0, "Levy-Tsallis fit", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);
    leg->Draw();

    TPaveText *box = new TPaveText(0.60, 0.7, 0.98, 0.93, "NDC");
    box->SetFillColor(0); // white background
    box->SetBorderSize(1);
    box->SetTextAlign(12); // left align
    box->SetTextSize(0.03);
    box->SetTextFont(42);
    box->AddText(Form("#chi^{2} / ndf                      %.3f / %d", fitFcnf0->GetChisquare(), fitFcnf0->GetNDF()));
    box->AddText(Form("n                        %.3f #pm %.3f", fitFcnf0->GetParameter(0), fitFcnf0->GetParError(0)));
    box->AddText(Form("dN/dy        %.5f #pm %.5f", fitFcnf0->GetParameter(1), fitFcnf0->GetParError(1)));
    box->AddText(Form("T                        %.3f #pm %.3f", fitFcnf0->GetParameter(3), fitFcnf0->GetParError(3)));
    box->Draw();
    // cCorrectedf0Fit->SaveAs((path + "plots/Yieldf0_CoherentFit.png").c_str());
    // cCorrectedf0Fit->SaveAs((path + "plots/Yieldf0_CoherentFitPhi0.png").c_str());

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
    hf21->SetMinimum(2e-8);
    hf21->SetMarkerColor(kBlue);
    hf21->SetLineColor(kBlue);
    hf21->SetStats(1);
    hf21->Draw("pe");
    hf22->SetMarkerColor(kBlue);
    hf22->SetLineColor(kBlue);
    hf22->SetFillStyle(0);
    hf22->SetMinimum(2e-8);
    hf22->Draw("e2 same");
    fitFcnf2->SetLineWidth(2);
    fitFcnf2->SetLineStyle(2);
    fitFcnf2->SetLineColor(kRed);
    fitFcnf2->Draw("l same");

    // Create legend with physics information
    TLegend *leg2 = new TLegend(0.19, 0.2, 0.4, 0.45);
    // leg2->AddEntry((TObject *)0, "ALICE", "");
    leg2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg2->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    leg2->AddEntry(hf21, "f_{2}'(1525) spectra", "pe");
    leg2->AddEntry(fitFcnf2, "Levy-Tsallis fit", "l");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.032);
    leg2->Draw();

    TPaveText *box2 = new TPaveText(0.60, 0.7, 0.98, 0.93, "NDC");
    box2->SetFillColor(0); // white background
    box2->SetBorderSize(1);
    box2->SetTextAlign(12); // left align
    box2->SetTextSize(0.03);
    box2->SetTextFont(42);
    box2->AddText(Form("#chi^{2} / ndf                      %.3f / %d", fitFcnf2->GetChisquare(), fitFcnf2->GetNDF()));
    box2->AddText(Form("n                        %.3f #pm %.3f", fitFcnf2->GetParameter(0), fitFcnf2->GetParError(0)));
    box2->AddText(Form("dN/dy        %.5f #pm %.5f", fitFcnf2->GetParameter(1), fitFcnf2->GetParError(1)));
    box2->AddText(Form("T                        %.3f #pm %.3f", fitFcnf2->GetParameter(3), fitFcnf2->GetParError(3)));
    box2->Draw();
    // cCorrectedf2Fit->SaveAs((path + "plots/Yieldf2_CoherentFit.png").c_str());
    // cCorrectedf2Fit->SaveAs((path + "plots/Yieldf2_CoherentFitPhi0.png").c_str());

    /// Compare the spectra from default
    ////****************************Spectra compare************************************
    TCanvas *cSpectraf2 = new TCanvas("cSpectraf2", "Spectra comparison f2", 720, 720);
    SetCanvasStyle(cSpectraf2, 0.15, 0.05, 0.05, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cSpectraf2, pad1Size, pad2Size);
    cSpectraf2->cd(1);
    gPad->SetLogy();
    SetHistoQA(hYieldReweightedf2);
    SetHistoQA(hf21);
    hYieldReweightedf2->SetMarkerSize(1.3);
    hf21->SetMarkerSize(1.3);
    hYieldReweightedf2->SetLineColor(kBlue);
    hYieldReweightedf2->SetMarkerColor(kBlue);
    hYieldReweightedf2->SetMarkerStyle(20);
    hYieldReweightedf2->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldReweightedf2->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hYieldReweightedf2->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hYieldReweightedf2->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldReweightedf2->GetYaxis()->SetTitleOffset(1.55 * pad1Size);
    hYieldReweightedf2->SetMaximum(5e-2);
    // hYieldReweightedf2->SetMaximum(hYieldReweightedf2->GetMaximum() * 2);
    hYieldReweightedf2->SetMinimum(2e-7);
    hYieldReweightedf2->Draw();
    hf21->SetLineColor(kRed);
    hf21->SetMarkerColor(kRed);
    hf21->SetMarkerStyle(21);
    hf21->SetMinimum(2e-7);
    hf21->Draw("SAME");
    TLegend *legf2 = new TLegend(0.55, 0.7, 0.9, 0.9);
    legf2->SetBorderSize(0);
    legf2->SetFillStyle(0);
    legf2->SetTextSize(0.055);
    legf2->SetTextFont(42);
    legf2->AddEntry((TObject *)0, "f'_{2}(1525) spectra", "");
    legf2->AddEntry(hYieldReweightedf2, "Default", "p");
    legf2->AddEntry(hf21, "Coherent", "p");
    legf2->Draw();
    cSpectraf2->cd(2);
    TH1F *hRatiof2 = (TH1F *)hYieldReweightedf2->Clone("hRatiof2");
    hRatiof2->Divide(hf21);
    SetHistoQA(hRatiof2);
    hRatiof2->GetYaxis()->SetTitle("Default/Coherent");
    hRatiof2->GetYaxis()->SetTitleSize(0.033 / pad2Size);
    hRatiof2->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof2->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof2->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatiof2->GetYaxis()->SetTitleOffset(1.7 * pad2Size);
    hRatiof2->SetLineColor(kBlue);
    hRatiof2->SetMarkerColor(kBlue);
    hRatiof2->GetYaxis()->SetNdivisions(505);
    hRatiof2->SetMaximum(2.58);
    hRatiof2->SetMinimum(0.19);
    // hRatiof2->SetMaximum(1.58);
    // hRatiof2->SetMinimum(0.65);
    hRatiof2->Draw();
    TLine *lineRatiof2 = new TLine(0, 1, 15, 1);
    lineRatiof2->SetLineStyle(2);
    lineRatiof2->SetLineColor(kBlack);
    lineRatiof2->Draw();
    // cSpectraf2->SaveAs((path + "plots/SpectraComparison_f2_Coherent.png").c_str());
    // cSpectraf2->SaveAs((path + "plots/SpectraComparison_f2_CoherentPhi0.png").c_str());

    TCanvas *cSpectraf0 = new TCanvas("cSpectraf0", "Spectra comparison f0", 720, 720);
    SetCanvasStyle(cSpectraf0, 0.15, 0.05, 0.05, 0.15);
    canvas_style(cSpectraf0, pad1Size, pad2Size);
    cSpectraf0->cd(1);
    gPad->SetLogy();
    SetHistoQA(hYieldReweightedf0);
    SetHistoQA(hf01);
    hYieldReweightedf0->SetMarkerSize(1.3);
    hf01->SetMarkerSize(1.3);
    hYieldReweightedf0->SetLineColor(kBlue);
    hYieldReweightedf0->SetMarkerColor(kBlue);
    hYieldReweightedf0->SetMarkerStyle(20);
    hYieldReweightedf0->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldReweightedf0->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hYieldReweightedf0->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hYieldReweightedf0->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldReweightedf0->GetYaxis()->SetTitleOffset(1.55 * pad1Size);
    hYieldReweightedf0->SetMaximum(hYieldReweightedf0->GetMaximum() * 0.1);
    hYieldReweightedf0->SetMinimum(8e-8);
    hYieldReweightedf0->Draw();
    hf01->SetLineColor(kRed);
    hf01->SetMarkerColor(kRed);
    hf01->SetMarkerStyle(21);
    hf01->SetMinimum(2e-7);
    hf01->Draw("SAME");
    TLegend *legf0 = new TLegend(0.55, 0.7, 0.9, 0.9);
    legf0->SetBorderSize(0);
    legf0->SetFillStyle(0);
    legf0->SetTextSize(0.055);
    legf0->SetTextFont(42);
    legf0->AddEntry((TObject *)0, "f_{0}(1710) spectra", "");
    legf0->AddEntry(hYieldReweightedf0, "Default", "p");
    legf0->AddEntry(hf01, "Coherent", "p");
    legf0->Draw();
    cSpectraf0->cd(2);
    TH1F *hRatiof0 = (TH1F *)hYieldReweightedf0->Clone("hRatiof0");
    hRatiof0->Divide(hf01);
    SetHistoQA(hRatiof0);
    hRatiof0->GetYaxis()->SetTitle("Default/Coherent");
    hRatiof0->GetYaxis()->SetTitleSize(0.033 / pad2Size);
    hRatiof0->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof0->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof0->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatiof0->GetYaxis()->SetTitleOffset(1.7 * pad2Size);
    hRatiof0->SetLineColor(kBlue);
    hRatiof0->SetMarkerColor(kBlue);
    hRatiof0->GetYaxis()->SetNdivisions(505);
    hRatiof0->SetMaximum(1.58);
    hRatiof0->SetMinimum(0.19);
    // hRatiof0->SetMaximum(1.58);
    // hRatiof0->SetMinimum(0.65);
    hRatiof0->Draw();
    TLine *lineRatiof0 = new TLine(0, 1, 15, 1);
    lineRatiof0->SetLineStyle(2);
    lineRatiof0->SetLineColor(kBlack);
    lineRatiof0->Draw();
    // cSpectraf0->SaveAs((path + "plots/SpectraComparison_f0_Coherent.png").c_str());
    // cSpectraf0->SaveAs((path + "plots/SpectraComparison_f0_CoherentPhi0.png").c_str());
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
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
}