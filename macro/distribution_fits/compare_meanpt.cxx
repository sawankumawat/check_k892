#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"
#include "../spectra/YieldMean2.C"

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

void compare_meanpt()
{
    bool otherQAPlots = true;
    gStyle->SetOptStat(0);
    // gStyle->SetOptFit(1111);

    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/";

    string savePath = path + "/mult_0-100/Spectra";
    TFile *fReweightf0 = new TFile((path + "mult_0-100/Spectra/ReweighFacf0_Default2.root").c_str(), "read");
    TFile *fReweightf2 = new TFile((path + "mult_0-100/Spectra/ReweighFacf2_Default2.root").c_str(), "read");

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

    TH1F *hGenUnweighted = (TH1F *)fReweightf0->Get("Genf17102_proj_1_i0");
    TH1F *hRecUnweighted = (TH1F *)fReweightf0->Get("Recf1710_pt2_proj_1_i0");
    TH1F *hGenUnweighted2 = (TH1F *)fReweightf2->Get("Genf17102_proj_1_i0");
    TH1F *hRecUnweighted2 = (TH1F *)fReweightf2->Get("Recf1710_pt2_proj_1_i0");

    if (hGenReweighted == nullptr || hRecReweighted == nullptr || hYieldReweighted == nullptr)
    {
        cout << "Error reading reweighted histograms from file " << Form("%s/ReweightedSpectra.root", savePath.c_str()) << endl;
        return;
    }

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
    TH1F *hf23 = (TH1F *)hYieldReweighted2->Clone("hf23");
    TH1F *hf24 = (TH1F *)hYieldReweighted2->Clone("hf24");

    // For f0(1710)
    TH1F *hf01 = (TH1F *)hYieldReweighted->Clone("hf01");
    TH1F *hf02 = (TH1F *)hYieldReweighted->Clone("hf02");
    TH1F *hf03 = (TH1F *)hYieldReweighted->Clone("hf03");
    TH1F *hf04 = (TH1F *)hYieldReweighted->Clone("hf04");

    // Enable error tracking for histograms with manual error setting
    hf01->Sumw2();
    hf02->Sumw2();
    hf03->Sumw2();
    hf04->Sumw2();
    hf21->Sumw2();
    hf22->Sumw2();
    hf23->Sumw2();
    hf24->Sumw2();

    double relUncertLowpTExtrapolationf2 = 11.0711 / 100;
    double relUncertLowpTExtrapolationf0 = 7.24382 / 100;

    for (int i = 1; i <= hf21->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double sys1710 = hYieldSysf0->GetBinContent(i);
        double sys1525 = hYieldSysf2->GetBinContent(i);
        double yield1710 = hf01->GetBinContent(i + 1);
        double yield1525 = hf21->GetBinContent(i + 1);

        hf02->SetBinContent(i + 1, yield1710);
        hf02->SetBinError(i + 1, sys1710 * yield1710);
        // hf02->SetBinError(i + 1, totalRelUncertf0 * yield1710);
        hf22->SetBinContent(i + 1, yield1525);
        hf22->SetBinError(i + 1, sys1525 * yield1525);
        // hf22->SetBinError(i + 1, totalRelUncertf2 * yield1525);

        hf04->SetBinContent(i + 1, yield1710);
        hf04->SetBinError(i + 1, sys1710 * yield1710);
        hf24->SetBinContent(i + 1, yield1525);
        hf24->SetBinError(i + 1, sys1525 * yield1525);
    }
    Double_t min = 0.0;
    Double_t min2 = 1.0;
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

    TF1 *fitFcnf0_2 = new TF1("fitfunc2_2", FuncLavy, 0.0, 15.0, 4);
    fitFcnf0_2->SetParameter(0, 5.0);
    fitFcnf0_2->SetParameter(1, 0.5);
    fitFcnf0_2->FixParameter(2, 1.710);
    fitFcnf0_2->SetParameter(3, 0.35);
    fitFcnf0_2->SetParNames("n", "dn/dy", "mass", "T");

    TF1 *fitFcnf2 = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
    fitFcnf2->SetParameter(0, 5.0);
    // fitFcnf2->SetParameter(1, 0.05);
    fitFcnf2->SetParameter(1, 0.5);
    fitFcnf2->FixParameter(2, 1.525);
    fitFcnf2->SetParameter(3, 0.35);
    fitFcnf2->SetParNames("n", "dn/dy", "mass", "T");

    TF1 *fitFcnf2_2 = new TF1("fitfunc2", FuncLavy, 0.0, 15.0, 4);
    fitFcnf2_2->SetParameter(0, 5.0);
    fitFcnf2_2->SetParameter(1, 0.5);
    fitFcnf2_2->FixParameter(2, 1.525);
    fitFcnf2_2->SetParameter(3, 0.35);
    fitFcnf2_2->SetParNames("n", "dn/dy", "mass", "T");

    TH1 *houtf0 = YieldMean(hf01, hf02, fitFcnf0, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *houtf2 = YieldMean(hf21, hf22, fitFcnf2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    TH1 *houtf0_noextrapol = YieldMean(hf01, hf02, fitFcnf0_2, min2, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *houtf2_noextrapol = YieldMean(hf21, hf22, fitFcnf2_2, min2, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);


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
    gMeanPtvsMassMesons->SetMaximum(4.59);
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
    // Get the value of <pT> at the mass of f0(1710) from the fit function
    double meanPt_f0_fromFit_meson = pol1_meson->Eval(1.710);
    cout << "Mean pT of f0(1710) from fit is " << meanPt_f0_fromFit_meson << endl;

    TF1 *pol1_baryon = new TF1("pol1_baryon", "pol1", 0.1, 1.89);
    pol1_baryon->SetLineColor(kYellow + 2);
    pol1_baryon->SetLineStyle(2);
    gMeanPtvsMassBaryons->Fit(pol1_baryon, "R");
    double meanPt_f0_fromFit_baryon = pol1_baryon->Eval(1.710);
    cout << "Mean pT of f0(1710) from fit is " << meanPt_f0_fromFit_baryon << endl;

    // Draw the last marker (f2(1525)) and its error bar
    double f2_mass = 1.5173;
    double f2_meanpt = houtf2->GetBinContent(5);
    double f2_Staterr = houtf2->GetBinContent(6);
    double f2_SysErrLow = sqrt(houtf2->GetBinContent(7) * houtf2->GetBinContent(7) + relUncertLowpTExtrapolationf2 * relUncertLowpTExtrapolationf2);
    double f2_SysErrHigh = sqrt(houtf2->GetBinContent(8) * houtf2->GetBinContent(8) + relUncertLowpTExtrapolationf2 * relUncertLowpTExtrapolationf2);
    cout << "Sys Low f2(1525) " << houtf2->GetBinContent(7) << " Sys High " << houtf2->GetBinContent(8) << endl;
    cout << "Sys total low " << f2_SysErrLow << " Sys total high " << f2_SysErrHigh << endl;
    cout << "Mean pT of f2(1525) is " << f2_meanpt << " with stat error " << f2_Staterr << endl;

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

    // Draw for f0(1710) without the extrapolation
    double f2_meanpt_noextrapol = houtf2_noextrapol->GetBinContent(5);
    double f2_staterr_noextrapol = houtf2_noextrapol->GetBinContent(6);
    double f2_syserrLow_noextrapol = houtf2_noextrapol->GetBinContent(7);
    double f2_syserrHigh_noextrapol = houtf2_noextrapol->GetBinContent(8);
    TGraphErrors *graph_f2_noextrapol = new TGraphErrors(1);
    graph_f2_noextrapol->SetPoint(0, f2_mass, f2_meanpt_noextrapol);
    graph_f2_noextrapol->SetPointError(0, 0, f2_staterr_noextrapol);
    graph_f2_noextrapol->SetMarkerStyle(23);
    graph_f2_noextrapol->SetMarkerColor(kYellow + 3);
    graph_f2_noextrapol->SetMarkerSize(1.7);
    graph_f2_noextrapol->SetLineColor(kYellow + 3);
    graph_f2_noextrapol->SetLineWidth(2);
    graph_f2_noextrapol->Draw("P SAME");

    TGraphAsymmErrors *graph_f2_sys_noextrapol = new TGraphAsymmErrors(1);
    graph_f2_sys_noextrapol->SetPoint(0, f2_mass, f2_meanpt_noextrapol);
    graph_f2_sys_noextrapol->SetPointError(0, 0.02, 0.02, f2_syserrLow_noextrapol, f2_syserrHigh_noextrapol);
    graph_f2_sys_noextrapol->SetLineColor(kYellow + 3);
    graph_f2_sys_noextrapol->SetLineWidth(2);
    graph_f2_sys_noextrapol->SetFillStyle(0);
    graph_f2_sys_noextrapol->Draw("2 SAME");

    // Draw the last marker (f0(1710)) and its error bar
    double f0_mass = 1.710;
    double f0_meanpt = houtf0->GetBinContent(5);
    double f0_Staterr = houtf0->GetBinContent(6);
    double f0_SysErrLow = sqrt(houtf0->GetBinContent(7) * houtf0->GetBinContent(7) + relUncertLowpTExtrapolationf0 * relUncertLowpTExtrapolationf0);
    double f0_SysErrHigh = sqrt(houtf0->GetBinContent(8) * houtf0->GetBinContent(8) + relUncertLowpTExtrapolationf0 * relUncertLowpTExtrapolationf0);
    cout << "Sys Low f0(1710) " << houtf0->GetBinContent(7) << " Sys High " << houtf0->GetBinContent(8) << endl;
    cout << "Sys total low " << f0_SysErrLow << " Sys total high " << f0_SysErrHigh << endl;
    cout << "Mean pT of f0(1710) is " << f0_meanpt << " with stat error " << f0_Staterr << endl;

    double SigmaDeviationMeson = fabs(meanPt_f0_fromFit_meson - f0_meanpt) / sqrt(f0_SysErrLow * f0_SysErrLow + f0_SysErrHigh * f0_SysErrHigh + f0_Staterr * f0_Staterr);
    cout << "Sigma deviation of f0(1710) from mesonic fit is " << SigmaDeviationMeson << endl;
    // cout << "Difference in mean value "<< fabs(meanPt_f0_fromFit_meson - f0_meanpt) << " and total error " << sqrt(f0_SysErrLow * f0_SysErrLow + f0_SysErrHigh * f0_SysErrHigh + f0_Staterr * f0_Staterr) << endl;
    double SigmaDeviationBaryon = fabs(meanPt_f0_fromFit_baryon - f0_meanpt) / sqrt(f0_SysErrLow * f0_SysErrLow + f0_SysErrHigh * f0_SysErrHigh + f0_Staterr * f0_Staterr);
    cout << "Sigma deviation of f0(1710) from baryonic fit is " << SigmaDeviationBaryon << endl;
    // cout << "Difference in mean value "<< fabs(meanPt_f0_fromFit_baryon - f0_meanpt) << " and total error " << sqrt(f0_SysErrLow * f0_SysErrLow + f0_SysErrHigh * f0_SysErrHigh + f0_Staterr * f0_Staterr) << endl;

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

    // Draw for f0(1710) without the extrapolation
    double f0_meanpt_noextrapol = houtf0_noextrapol->GetBinContent(5);
    double f0_staterr_noextrapol = houtf0_noextrapol->GetBinContent(6);
    double f0_syserrLow_noextrapol = houtf0_noextrapol->GetBinContent(7);
    double f0_syserrHigh_noextrapol = houtf0_noextrapol->GetBinContent(8);
    TGraphErrors *graph_f0_noextrapol = new TGraphErrors(1);
    graph_f0_noextrapol->SetPoint(0, f0_mass, f0_meanpt_noextrapol);
    graph_f0_noextrapol->SetPointError(0, 0, f0_staterr_noextrapol);
    graph_f0_noextrapol->SetMarkerStyle(22);
    graph_f0_noextrapol->SetMarkerColor(kViolet);
    graph_f0_noextrapol->SetMarkerSize(1.7);
    graph_f0_noextrapol->SetLineColor(kViolet);
    graph_f0_noextrapol->SetLineWidth(2);
    graph_f0_noextrapol->Draw("P SAME");
    TGraphAsymmErrors *graph_f0_sys_noextrapol = new TGraphAsymmErrors(1);
    graph_f0_sys_noextrapol->SetPoint(0, f0_mass, f0_meanpt_noextrapol);
    graph_f0_sys_noextrapol->SetPointError(0, 0.02, 0.02, f0_syserrLow_noextrapol, f0_syserrHigh_noextrapol);
    graph_f0_sys_noextrapol->SetLineColor(kViolet);
    graph_f0_sys_noextrapol->SetLineWidth(2);
    graph_f0_sys_noextrapol->SetFillStyle(0);
    graph_f0_sys_noextrapol->Draw("2 SAME");

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
    latex.DrawLatex(1.5, houtf2->GetBinContent(5) + 0.8, "f'_{2}(1525)");
    latex.DrawLatex(1.710, houtf0->GetBinContent(5) + 0.8, "f_{0}(1710)");

    TLegend *legend4 = new TLegend(0.17, 0.6, 0.63, 0.92);
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->SetTextSize(0.033);
    legend4->AddEntry(gMeanPtvsMassMesons, "Mesons (13 TeV)", "p");
    legend4->AddEntry(gMeanPtvsMassBaryons, "Baryons (13 TeV)", "p");
    legend4->AddEntry(graph_f2, "f'_{2}(1525) (With extrapolation)", "p");
    legend4->AddEntry(graph_f2_noextrapol, "f'_{2}(1525) (Without extrapolation)", "p");
    legend4->AddEntry(graph_f0, "f_{0}(1710) (With extrapolation)", "p");
    legend4->AddEntry(graph_f0_noextrapol, "f_{0}(1710) (Without extrapolation)", "p");
    legend4->AddEntry(pol1_meson, "Pol 1", "l");
    legend4->Draw();

    TLegend *legend5 = new TLegend(0.52, 0.78, 0.87, 0.92);
    legend5->SetBorderSize(0);
    legend5->SetFillStyle(0);
    legend5->SetTextSize(0.03);
    legend5->AddEntry((TObject *)0, "ALICE", "");
    legend5->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legend5->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    // legend5->Draw();
    cMeanPt->SaveAs((savePath + "/plots/MeanPt_compare_extrapol.png").c_str());
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