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
    TFile *fReweightf0 = new TFile((path + "mult_0-100/Spectra/ReweighFacf0_Default4.root").c_str(), "read");
    TFile *fReweightf2 = new TFile((path + "mult_0-100/Spectra/ReweighFacf2_Default4.root").c_str(), "read");
    TFile *fSpectraMC = new TFile((path + "mult_0-100/Spectra/spectra_Default4.root").c_str(), "read"); // Without reweighting the MC

    TFile *fSpectraToy = new TFile((path + "mult_0-100/Spectra/spectra_ToyMC2.root").c_str(), "read");
    TFile *fSingleBWf2 = new TFile((path + "f2_fit_1rBW.root").c_str(), "read");

    TFile *fReweightf0_var = new TFile((path + "mult_0-100/Spectra/ReweighFacf0_DefaultWithoutSigLoss.root").c_str(), "read");
    TFile *fReweightf2_var = new TFile((path + "mult_0-100/Spectra/ReweighFacf2_DefaultWithoutSigLoss.root").c_str(), "read");

    if (fReweightf0->IsZombie() || fReweightf2->IsZombie() || fSpectraToy->IsZombie())
    {
        cout << "Error opening reweighting or toy spectra files" << endl;
        return;
    }

    TH1F *hYield1525 = (TH1F *)fSpectraMC->Get("hYield1525Corrected");
    TH1F *hYield1710 = (TH1F *)fSpectraMC->Get("hYield1710Corrected");

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
    TH1F *hYieldReweightedf0 = (TH1F *)fReweightf0->Get(Form("hYield1710Corrected_%s", indexStr.c_str()));
    TH1F *hYieldReweightedf2 = (TH1F *)fReweightf2->Get(Form("hYield1525Corrected_%s", indexStr2.c_str()));

    if (hYieldReweightedf0 == nullptr)
    {
        cout << "Error reading reweighted histograms from file " << Form("%s/ReweightedSpectra.root", savePath.c_str()) << endl;
        return;
    }

    int maxIndexReweightedf0_var = FindHighestIndex(fReweightf0_var, "Genf17102_proj_1_");
    int maxIndexReweightedf2_var = FindHighestIndex(fReweightf2_var, "Genf17102_proj_1_");

    // // For comparison of with/without signal loss corrections
    // TH1F *hYieldToyf0 = (TH1F *)fReweightf0_var->Get(Form("hYield1710Corrected_%s", ("i" + to_string(maxIndexReweightedf0_var)).c_str()));
    // TH1F *hYieldToyf2 = (TH1F *)fReweightf2_var->Get(Form("hYield1525Corrected_%s", ("i" + to_string(maxIndexReweightedf2_var)).c_str()));

    //// For toy model comparison
    // TH1F *hYieldToyf0 = (TH1F *)fSpectraToy->Get("hYield1710Corrected");
    // TH1F *hYieldToyf2 = (TH1F *)fSpectraToy->Get("hYield1525Corrected");

    //// For comparison with single BW fit results
    TH1F *hYieldToyf2 = (TH1F *)fSingleBWf2->Get("hCorrectedYield_1525");
    TH1F *hYieldToyf0 = (TH1F *)fSingleBWf2->Get("hCorrectedYield_1525"); // dummy histogram

    // if (hYieldToyf0 == nullptr || hYieldToyf2 == nullptr)
    // {
    //     cout << "Error reading toy spectra histograms from file " << Form("%s/spectra_ToyMC.root", savePath.c_str()) << endl;
    //     return;
    // }

    // TFile *fSigEventLossOld = new TFile("Loss_phi_mult0-100.root", "READ");
    // TFile *fSigEventLossNew = new TFile("../SignalLossPhiINEL.root", "READ");

    // TH1F *hEvbySigLossOld = (TH1F *)fSigEventLossOld->Get("hSignalLoss");
    // TH1F *hEvbySigLossNew = (TH1F *)fSigEventLossNew->Get("hSignalLoss");

    // TH1F *hNewbyOld = (TH1F *)hEvbySigLossNew->Clone("hNewbyOld");
    // hNewbyOld->Divide(hEvbySigLossOld);

    // for (int i = 1; i <= hNewbyOld->GetNbinsX(); i++)
    // {
    //     double binContentf2 = hYieldToyf2->GetBinContent(i + 1);
    //     double binContentf0 = hYieldToyf0->GetBinContent(i + 1);
    //     double ratio = hNewbyOld->GetBinContent(i);
    //     double newBinContentf2 = binContentf2 * ratio;
    //     double newBinContentf0 = binContentf0 * ratio;
    //     hYieldToyf2->SetBinContent(i + 1, newBinContentf2);
    //     hYieldToyf0->SetBinContent(i + 1, newBinContentf0);
    // }

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

    cout << "Total bins in the yield histogram is " << hYieldReweightedf0->GetNbinsX() << endl;
    cout << "Total bins in the toy yield histogram is " << hYieldToyf0->GetNbinsX() << endl;

    ////****************************Spectra compare************************************
    TCanvas *cSpectraf2 = new TCanvas("cSpectraf2", "Spectra comparison f2", 720, 720);
    SetCanvasStyle(cSpectraf2, 0.15, 0.05, 0.05, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cSpectraf2, pad1Size, pad2Size);
    cSpectraf2->cd(1);
    gPad->SetLogy();
    // TH1F *hYieldDefaultClonef2 = (TH1F *)hYield1525->Clone("hYieldDefaultClonef2"); // For toy model
    TH1F *hYieldDefaultClonef2 = (TH1F *)hYieldReweightedf2->Clone("hYieldDefaultClonef2"); // For with/without signal loss
    SetHistoQA(hYieldDefaultClonef2);
    SetHistoQA(hYieldToyf2);
    hYieldDefaultClonef2->SetMarkerSize(1.3);
    hYieldToyf2->SetMarkerSize(1.3);
    hYieldDefaultClonef2->SetLineColor(kBlue);
    hYieldDefaultClonef2->SetMarkerColor(kBlue);
    hYieldDefaultClonef2->SetMarkerStyle(20);
    hYieldDefaultClonef2->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldDefaultClonef2->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hYieldDefaultClonef2->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hYieldDefaultClonef2->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldDefaultClonef2->GetYaxis()->SetTitleOffset(1.55 * pad1Size);
    hYieldDefaultClonef2->SetMaximum(hYieldDefaultClonef2->GetMaximum() * 200);
    hYieldDefaultClonef2->SetMinimum(2e-7);
    hYieldDefaultClonef2->Draw();
    hYieldToyf2->SetLineColor(kRed);
    hYieldToyf2->SetMarkerColor(kRed);
    hYieldToyf2->SetMarkerStyle(21);
    hYieldToyf2->SetMinimum(2e-7);
    hYieldToyf2->Draw("SAME");
    TLegend *legf2 = new TLegend(0.55, 0.7, 0.9, 0.9);
    legf2->SetBorderSize(0);
    legf2->SetFillStyle(0);
    legf2->SetTextSize(0.055);
    legf2->SetTextFont(42);
    legf2->SetHeader("f_{2}(1525) spectra");

    // legf2->AddEntry(hYieldDefaultClonef2, "Without Signal Loss", "p");
    // legf2->AddEntry(hYieldToyf2, "With Signal Loss", "p");
    //  legf2->Draw();

    // legf2->AddEntry(hYieldDefaultClonef2, "f_{2}(1525) (MC)", "p");
    // legf2->AddEntry(hYieldToyf2, "f_{2}(1525) (Toy model)", "p");
    // legf2->Draw();

    legf2->AddEntry(hYieldDefaultClonef2, "Default", "p");
    legf2->AddEntry(hYieldToyf2, "Single BW Fit", "p");
    legf2->Draw();

    cSpectraf2->cd(2);
    TH1F *hRatiof2 = (TH1F *)hYieldDefaultClonef2->Clone("hRatiof2");
    hRatiof2->Divide(hYieldToyf2);
    SetHistoQA(hRatiof2);
    // hRatiof2->GetYaxis()->SetTitle("MC/Toy model");
    // hRatiof2->GetYaxis()->SetTitle("Default/NoSigLoss");
    hRatiof2->GetYaxis()->SetTitle("Default/1rBW Fit");
    hRatiof2->GetYaxis()->SetTitleSize(0.033 / pad2Size);
    hRatiof2->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof2->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof2->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatiof2->GetYaxis()->SetTitleOffset(1.4 * pad2Size);
    hRatiof2->SetLineColor(kBlue);
    hRatiof2->SetMarkerColor(kBlue);
    hRatiof2->GetYaxis()->SetNdivisions(505);
    hRatiof2->SetMaximum(1.79);
    hRatiof2->SetMinimum(0.35);
    // hRatiof2->SetMaximum(4.1);
    // hRatiof2->SetMinimum(0.25);
    hRatiof2->Draw();
    TLine *lineRatiof2 = new TLine(0, 1, 15, 1);
    lineRatiof2->SetLineStyle(2);
    lineRatiof2->SetLineColor(kBlack);
    lineRatiof2->Draw();

    TCanvas *cSpectraf0 = new TCanvas("cSpectraf0", "Spectra comparison f0", 720, 720);
    SetCanvasStyle(cSpectraf0, 0.15, 0.05, 0.05, 0.15);
    canvas_style(cSpectraf0, pad1Size, pad2Size);
    cSpectraf0->cd(1);
    gPad->SetLogy();
    // TH1F *hYieldf0DefaultClone = (TH1F *)hYield1710->Clone("hYieldf0DefaultClone"); // For toy model
    TH1F *hYieldf0DefaultClone = (TH1F *)hYieldReweightedf0->Clone("hYieldf0DefaultClone"); // For with/without signal loss
    SetHistoQA(hYieldf0DefaultClone);
    SetHistoQA(hYieldToyf0);
    hYieldf0DefaultClone->SetMarkerSize(1.3);
    hYieldToyf0->SetMarkerSize(1.3);
    hYieldf0DefaultClone->SetLineColor(kBlue);
    hYieldf0DefaultClone->SetMarkerColor(kBlue);
    hYieldf0DefaultClone->SetMarkerStyle(20);
    hYieldf0DefaultClone->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldf0DefaultClone->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hYieldf0DefaultClone->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hYieldf0DefaultClone->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hYieldf0DefaultClone->GetYaxis()->SetTitleOffset(1.55 * pad1Size);
    // hYieldf0DefaultClone->SetMaximum(hYieldf0DefaultClone->GetMaximum() * 200);
    hYieldf0DefaultClone->SetMaximum(hYieldf0DefaultClone->GetMaximum() * 0.02);
    hYieldf0DefaultClone->SetMinimum(8e-8);
    hYieldf0DefaultClone->Draw();
    hYieldToyf0->SetLineColor(kRed);
    hYieldToyf0->SetMarkerColor(kRed);
    hYieldToyf0->SetMarkerStyle(21);
    hYieldToyf0->SetMinimum(2e-7);
    hYieldToyf0->Draw("SAME");
    TLegend *legf0 = new TLegend(0.55, 0.7, 0.9, 0.9);
    legf0->SetBorderSize(0);
    legf0->SetFillStyle(0);
    legf0->SetTextSize(0.055);
    legf0->SetTextFont(42);
    legf0->SetHeader("f_{0}(1710) spectra");
    // legf0->AddEntry(hYieldf0DefaultClone, "f_{0}(1710) (MC)", "p");
    // legf0->AddEntry(hYieldToyf0, "f_{0}(1710) (Toy model)", "p");
    legf0->AddEntry(hYieldf0DefaultClone, "Without Signal Loss", "p");
    legf0->AddEntry(hYieldToyf0, "With Signal Loss", "p");
    legf0->Draw();

    cSpectraf0->cd(2);
    TH1F *hRatiof0 = (TH1F *)hYieldf0DefaultClone->Clone("hRatiof0");
    hRatiof0->Divide(hYieldToyf0);
    SetHistoQA(hRatiof0);
    // hRatiof0->GetYaxis()->SetTitle("MC/Toy model");
    hRatiof0->GetYaxis()->SetTitle("Default/NoSigLoss");
    hRatiof0->GetYaxis()->SetTitleSize(0.033 / pad2Size);
    hRatiof0->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof0->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hRatiof0->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatiof0->GetYaxis()->SetTitleOffset(1.4 * pad2Size);
    hRatiof0->SetLineColor(kBlue);
    hRatiof0->SetMarkerColor(kBlue);
    hRatiof0->GetYaxis()->SetNdivisions(505);
    hRatiof0->SetMaximum(1.79);
    hRatiof0->SetMinimum(0.35);
    // hRatiof0->SetMaximum(4.1);
    // hRatiof0->SetMinimum(0.25);
    hRatiof0->Draw();
    TLine *lineRatiof0 = new TLine(0, 1, 15, 1);
    lineRatiof0->SetLineStyle(2);
    lineRatiof0->SetLineColor(kBlack);
    lineRatiof0->Draw();

    cSpectraf2->SaveAs((savePath + "/plots/SpectraCompare1rBWFit_f2.png").c_str());

    // cSpectraf0->SaveAs((savePath + "/plots/SpectraCompareSigLoss_f0.png").c_str());
    // cSpectraf2->SaveAs((savePath + "/plots/SpectraCompareSigLoss_f2.png").c_str());

    // cSpectraf0->SaveAs((savePath + "/plots/SpectraCompareToy_f0.png").c_str());
    // cSpectraf2->SaveAs((savePath + "/plots/SpectraCompareToy_f2.png").c_str());

    ////****************************<pT> calcuation**********************************////
    // For f2'(1525)
    TH1F *hf21 = (TH1F *)hYieldReweightedf2->Clone("hf21");
    TH1F *hf22 = (TH1F *)hYieldReweightedf2->Clone("hf22");

    //*********For comparison with fit range variation*************
    // TH1F *hf23 = (TH1F *)hYieldReweightedf2->Clone("hf23");
    // TH1F *hf24 = (TH1F *)hYieldReweightedf2->Clone("hf24");

    // *********For comparison with toy model, or compare with/without signal loss*************
    TH1F *hf23 = (TH1F *)hYieldToyf2->Clone("hf23");
    TH1F *hf24 = (TH1F *)hYieldToyf2->Clone("hf24");

    // //*********For comparison with unweighted efficiency**************
    // TH1F *hf23 = (TH1F *)hYield1525->Clone("hf23");
    // TH1F *hf24 = (TH1F *)hYield1525->Clone("hf24");

    // For f0(1710)
    TH1F *hf01 = (TH1F *)hYieldReweightedf0->Clone("hf01");
    TH1F *hf02 = (TH1F *)hYieldReweightedf0->Clone("hf02");

    //*********For comparison with fit range variation*************
    // TH1F *hf03 = (TH1F *)hYieldReweightedf0->Clone("hf03");
    // TH1F *hf04 = (TH1F *)hYieldReweightedf0->Clone("hf04");

    // *********For comparison with toy model or with/without signal loss*************
    TH1F *hf03 = (TH1F *)hYieldToyf0->Clone("hf03");
    TH1F *hf04 = (TH1F *)hYieldToyf0->Clone("hf04");

    // //*********For comparison with unweighted efficiency**************
    // TH1F *hf03 = (TH1F *)hYield1710->Clone("hf03");
    // TH1F *hf04 = (TH1F *)hYield1710->Clone("hf04");

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

        double yieldToy1710 = hf03->GetBinContent(i + 1);
        double yieldToy1525 = hf23->GetBinContent(i + 1);

        hf02->SetBinContent(i + 1, yield1710);
        hf02->SetBinError(i + 1, sys1710 * yield1710);
        hf22->SetBinContent(i + 1, yield1525);
        hf22->SetBinError(i + 1, sys1525 * yield1525);

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
    // Option_t *opt = "REBMS0+";
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 1.0;
    // Double_t maxfit = 10.0;
    Double_t maxfit = 15.0;

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

    // TH1 *houtf0_noextrapol = YieldMean(hf01, hf02, fitFcnf0_2, min2, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // TH1 *houtf2_noextrapol = YieldMean(hf21, hf22, fitFcnf2_2, min2, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *houtf0_noextrapol = YieldMean(hf03, hf04, fitFcnf0_2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *houtf2_noextrapol = YieldMean(hf23, hf24, fitFcnf2_2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

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

                // cout << "Mean pT of " << particlesLatex[i] << " at 13 TeV is " << meanPtAt13TeV[i] << " +- " << meanPtAt13TeV_err[i] << endl;
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

    TF1 *pol1_meson = new TF1("pol1_meson", "pol1", 0.1, 1.89);
    pol1_meson->SetLineColor(kYellow + 2);
    pol1_meson->SetLineStyle(2);
    gMeanPtvsMassMesons->Fit(pol1_meson, "R");
    // Get the value of <pT> at the mass of f0(1710) from the fit function
    double meanPt_f0_fromFit_meson = pol1_meson->Eval(1.710);
    // cout << "Mean pT of f0(1710) from fit is " << meanPt_f0_fromFit_meson << endl;

    TF1 *pol1_baryon = new TF1("pol1_baryon", "pol1", 0.1, 1.89);
    pol1_baryon->SetLineColor(kYellow + 2);
    pol1_baryon->SetLineStyle(2);
    gMeanPtvsMassBaryons->Fit(pol1_baryon, "R");
    double meanPt_f0_fromFit_baryon = pol1_baryon->Eval(1.710);
    // cout << "Mean pT of f0(1710) from fit is " << meanPt_f0_fromFit_baryon << endl;

    // Draw the last marker (f2(1525)) and its error bar
    double f2_mass = 1.5173;
    double f2_meanpt = houtf2->GetBinContent(5);
    double f2_Staterr = houtf2->GetBinContent(6);
    double f2_SysErr = sqrt(houtf2->GetBinContent(7) * houtf2->GetBinContent(7) + relUncertLowpTExtrapolationf2 * relUncertLowpTExtrapolationf2);
    // double f2_SysErrLow = sqrt(houtf2->GetBinContent(8) * houtf2->GetBinContent(8) + relUncertLowpTExtrapolationf2 * relUncertLowpTExtrapolationf2);
    cout << "Mean pT of f2(1525) is " << f2_meanpt << " stat error " << f2_Staterr << " sys error " << f2_SysErr << endl;
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
    graph_f2_sys->SetPointError(0, 0.02, 0.02, f2_SysErr, f2_SysErr);
    graph_f2_sys->SetLineColor(f2_color);
    graph_f2_sys->SetLineWidth(2);
    graph_f2_sys->SetFillStyle(0);
    // graph_f2_sys->Draw("2 SAME");

    // Draw for <pT> of f0(1710) for variation case
    /// For toy model
    double f2_meanpt_noextrapol = houtf2_noextrapol->GetBinContent(5);
    double f2_staterr_noextrapol = houtf2_noextrapol->GetBinContent(6);
    double f2_syserr_noextrapol = houtf2_noextrapol->GetBinContent(7);
    // double f2_syserrHigh_noextrapol = houtf2_noextrapol->GetBinContent(8);

    // // Directly giving the values
    // double f2_meanpt_noextrapol = 1.636;
    // double f2_staterr_noextrapol = 0.043;
    // double f2_syserr_noextrapol = 0.062;
    // // double f2_syserrLow_noextrapol = 0.033;

    cout << "f2(1525) mean pT in variation is " << f2_meanpt_noextrapol << endl;
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
    graph_f2_sys_noextrapol->SetPointError(0, 0.02, 0.02, f2_syserr_noextrapol, f2_syserr_noextrapol);
    graph_f2_sys_noextrapol->SetLineColor(kYellow + 3);
    graph_f2_sys_noextrapol->SetLineWidth(2);
    graph_f2_sys_noextrapol->SetFillStyle(0);
    // graph_f2_sys_noextrapol->Draw("2 SAME");

    // Draw the last marker (f0(1710)) and its error bar
    double f0_mass = 1.710;
    double f0_meanpt = houtf0->GetBinContent(5);
    double f0_Staterr = houtf0->GetBinContent(6);
    double f0_SysErr = sqrt(houtf0->GetBinContent(7) * houtf0->GetBinContent(7) + relUncertLowpTExtrapolationf0 * relUncertLowpTExtrapolationf0);
    // double f0_SysErrLow = sqrt(houtf0->GetBinContent(8) * houtf0->GetBinContent(8) + relUncertLowpTExtrapolationf0 * relUncertLowpTExtrapolationf0);
    cout << "Mean pT of f0(1710) is " << f0_meanpt << " stat error " << f0_Staterr << " Sys " << f0_SysErr << endl;

    double SigmaDeviationMeson = fabs(meanPt_f0_fromFit_meson - f0_meanpt) / sqrt(f0_SysErr * f0_SysErr + f0_Staterr * f0_Staterr);
    double SigmaDeviationBaryon = fabs(meanPt_f0_fromFit_baryon - f0_meanpt) / sqrt(f0_SysErr * f0_SysErr + f0_Staterr * f0_Staterr);

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
    graph_f0_sys->SetPointError(0, 0.02, 0.02, f0_SysErr, f0_SysErr);
    graph_f0_sys->SetLineColor(f0_color);
    graph_f0_sys->SetLineWidth(2);
    graph_f0_sys->SetFillStyle(0);
    // graph_f0_sys->Draw("2 SAME");

    //// Draw for <pT> for f0(1710) for the variation case
    double f0_meanpt_noextrapol = houtf0_noextrapol->GetBinContent(5);
    double f0_staterr_noextrapol = houtf0_noextrapol->GetBinContent(6);
    double f0_syserr_noextrapol = houtf0_noextrapol->GetBinContent(7);
    // double f0_syserrLow_noextrapol = houtf0_noextrapol->GetBinContent(8);

    // double f0_meanpt_noextrapol = 2.336;
    // double f0_staterr_noextrapol = 0.105;
    // double f0_syserr_noextrapol = 0.169;
    // double f0_syserrLow_noextrapol = 0.104;
    cout << "f0(1710) mean pT in variation is " << f0_meanpt_noextrapol << endl;

    TGraphErrors *graph_f0_noextrapol = new TGraphErrors(1);
    graph_f0_noextrapol->SetPoint(0, f0_mass, f0_meanpt_noextrapol);
    graph_f0_noextrapol->SetPointError(0, 0, f0_staterr_noextrapol);
    graph_f0_noextrapol->SetMarkerStyle(23);
    graph_f0_noextrapol->SetMarkerColor(kCyan + 1);
    graph_f0_noextrapol->SetMarkerSize(1.7);
    graph_f0_noextrapol->SetLineColor(kCyan + 1);
    graph_f0_noextrapol->SetLineWidth(2);
    // graph_f0_noextrapol->Draw("P SAME"); // Commented for 1rBW f2 fit comparison
    TGraphAsymmErrors *graph_f0_sys_noextrapol = new TGraphAsymmErrors(1);
    graph_f0_sys_noextrapol->SetPoint(0, f0_mass, f0_meanpt_noextrapol);
    graph_f0_sys_noextrapol->SetPointError(0, 0.02, 0.02, f0_syserr_noextrapol, f0_syserr_noextrapol);
    graph_f0_sys_noextrapol->SetLineColor(kCyan + 1);
    graph_f0_sys_noextrapol->SetLineWidth(2);
    graph_f0_sys_noextrapol->SetFillStyle(0);
    // graph_f0_sys_noextrapol->Draw("2 SAME");

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
    latex.DrawLatex(1.5, houtf2->GetBinContent(5) + 0.4, "f'_{2}(1525)");
    latex.DrawLatex(1.710, houtf0->GetBinContent(5) + 0.4, "f_{0}(1710)");

    TLegend *legend4 = new TLegend(0.17, 0.6, 0.63, 0.92);
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->SetTextSize(0.033);
    legend4->AddEntry(gMeanPtvsMassMesons, "Mesons (13 TeV)", "p");
    legend4->AddEntry(gMeanPtvsMassBaryons, "Baryons (13 TeV)", "p");

    //// For comparison between weighted and unweighted efficiency
    // legend4->AddEntry(graph_f2_noextrapol, "f'_{2}(1525)", "p");
    // legend4->AddEntry(graph_f2, "f'_{2}(1525) (Reweighted)", "p");
    // legend4->AddEntry(graph_f0_noextrapol, "f_{0}(1710)", "p");
    // legend4->AddEntry(graph_f0, "f_{0}(1710) (Reweighted)", "p");

    // //// For comparison with toy model
    // legend4->AddEntry(graph_f2, "f'_{2}(1525) (MC)", "p");
    // legend4->AddEntry(graph_f2_noextrapol, "f'_{2}(1525) (Toy Model)", "p");
    // legend4->AddEntry(graph_f0, "f_{0}(1710) (MC)", "p");
    // legend4->AddEntry(graph_f0_noextrapol, "f_{0}(1710) (Toy Model)", "p");
    // legend4->Draw();

    //// For comparison with 1rBW fit to f2(1525)
    legend4->AddEntry(graph_f2, "f'_{2}(1525) (Default)", "p");
    legend4->AddEntry(graph_f2_noextrapol, "f'_{2}(1525) (1rBW Fit)", "p");
    legend4->AddEntry(graph_f0, "f_{0}(1710) (Default)", "p");
    legend4->Draw();

    //// For comparison with and without extrapolation
    // legend4->AddEntry(graph_f2, "f'_{2}(1525) (With extrapolation)", "p");
    // legend4->AddEntry(graph_f2_noextrapol, "f'_{2}(1525) (Without extrapolation)", "p");
    // legend4->AddEntry(graph_f0, "f_{0}(1710) (With extrapolation)", "p");
    // legend4->AddEntry(graph_f0_noextrapol, "f_{0}(1710) (Without extrapolation)", "p");

    // // For comparison with different fit ranges
    // legend4->AddEntry(graph_f2, "f'_{2}(1525) (Fit 0-15 GeV/c)", "p");
    // legend4->AddEntry(graph_f2_noextrapol, "f'_{2}(1525) (Fit 0-10 GeV/c)", "p");
    // legend4->AddEntry(graph_f0, "f_{0}(1710) (Fit 0-15 GeV/c)", "p");
    // legend4->AddEntry(graph_f0_noextrapol, "f_{0}(1710) (Fit 0-10 GeV/c)", "p");
    // legend4->AddEntry(pol1_meson, "Pol 1", "l");
    // legend4->Draw();

    TLegend *legend5 = new TLegend(0.55, 0.78, 0.90, 0.92);
    legend5->SetBorderSize(0);
    legend5->SetFillStyle(0);
    legend5->SetTextSize(0.03);
    legend5->AddEntry((TObject *)0, "ALICE", "");
    legend5->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legend5->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    legend5->Draw();

    cMeanPt->SaveAs((savePath + "/plots/MeanPt_compare_1rBWFit_f2.png").c_str());
    // cMeanPt->SaveAs((savePath + "/plots/MeanPt_compare_SigLoss.png").c_str());
    // cMeanPt->SaveAs((savePath + "/plots/MeanPt_compare_extrapol.png").c_str());
    // cMeanPt->SaveAs((savePath + "/plots/MeanPtVsMass_ReweightingCompare.png").c_str());
    // cMeanPt->SaveAs((savePath + "/plots/MeanPtVsMass_ToyModelCompare.png").c_str());
    // cMeanPt->SaveAs((savePath + "/plots/MeanPtVsMass_FitRangeCompare.png").c_str());

    // // Similarly plot the dN/dy/(2J+1) as a function of particle mass
    // vector<vector<float>> dNdyvalues_13TeV = {
    //     {4.775, 0.001, 0.243, 1.0},  // 4.775 ± 0.001 ± 0.243 (pion) (it is pi+ + pi-)
    //     {6.205, 0.004, 0.303, 1e-1}, // (6.205 ± 0.004 ± 0.303) × 10⁻¹ (kaon) (K+ + K-)
    //     {3.192, 0.004, 0.111, 1e-1}, // (3.192 ± 0.004 ± 0.111) × 10⁻¹ (K0s)
    //     {2.098, 0.016, 0.200, 1e-1}, // (2.098 ± 0.016 ± 0.200) × 10⁻¹ (Kstar) (it is K*0 + anti-K*0)
    //     {3.734, 0.040, 0.213, 1e-2}, // (3.734 ± 0.040 ± 0.213) × 10⁻² (Phi)
    //     {2.750, 0.002, 0.188, 1e-1}, // (2.750 ± 0.002 ± 0.188) × 10⁻¹ (proton) (it is p + p_bar)
    //     {1.807, 0.005, 0.102, 1e-1}, // (1.807 ± 0.005 ± 0.102) × 10⁻¹ (Lambda) (it is Lambda + anti-Lambda)
    //     {1.980, 0.012, 0.082, 1e-2}, // (1.980 ± 0.012 ± 0.082) × 10⁻² (Xi) (it is Xi- + anti-Xi+)
    //     {1.846, 0.046, 0.122, 1e-3}  // (1.846 ± 0.046 ± 0.122) × 10⁻³ (Omega) (it is Omega- + anti-Omega+)
    // };
    // int Jvalues[9] = {0, 0, 0, 1, 1, 1 / 2, 1 / 2, 1 / 2, 3 / 2};
    // int averageFactor[9] = {2, 2, 1, 2, 1, 2, 2, 2, 2};
    // TGraphAsymmErrors *gdNdyvsMassMesons = new TGraphAsymmErrors();
    // TGraphAsymmErrors *gdNdyvsMassBaryons = new TGraphAsymmErrors();
    // double dNdyAll[9];
    // for (int i = 0; i < totalParticles; i++)
    // {
    //     double dNdy = dNdyvalues_13TeV[i][0] * dNdyvalues_13TeV[i][3] / ((2 * Jvalues[i] + 1) * averageFactor[i]); // apply the scaling factor
    //     dNdyAll[i] = dNdy;
    //     // double dNdy_err_stat = dNdyvalues_13TeV[i][1] * dNdyvalues_13TeV[i][3] / (2 * Jvalues[i] + 1);
    //     double dNdy_err_sys = dNdyvalues_13TeV[i][2] * dNdyvalues_13TeV[i][3] / ((2 * Jvalues[i] + 1) * averageFactor[i]);
    //     cout << "Particle " << particlesLatex[i] << " dN/dy/(2J+1) = " << dNdy << " with sys error " << dNdy_err_sys << endl;
    //     if (i < 5)
    //     {
    //         gdNdyvsMassMesons->SetPoint(i, particleMass[i], dNdy);
    //         gdNdyvsMassMesons->SetPointError(i, 0.02, 0.02, dNdy_err_sys, dNdy_err_sys);
    //     }
    //     else
    //     {
    //         gdNdyvsMassBaryons->SetPoint(i - 5, particleMass[i], dNdy);
    //         gdNdyvsMassBaryons->SetPointError(i - 5, 0.02, 0.02, dNdy_err_sys, dNdy_err_sys);
    //     }
    // }
    // SetGraphStyleCommon(gdNdyvsMassMesons);
    // gdNdyvsMassMesons->SetMarkerStyle(22);
    // gdNdyvsMassMesons->SetMarkerColor(kMagenta);
    // gdNdyvsMassMesons->SetLineColor(kMagenta);
    // gdNdyvsMassMesons->SetLineWidth(2);
    // gdNdyvsMassMesons->SetMarkerSize(1.7);
    // gdNdyvsMassMesons->SetFillStyle(0);
    // SetGraphStyleCommon(gdNdyvsMassBaryons);
    // gdNdyvsMassBaryons->SetMarkerStyle(34);
    // gdNdyvsMassBaryons->SetMarkerColor(kRed);
    // gdNdyvsMassBaryons->SetLineColor(kRed);
    // gdNdyvsMassBaryons->SetLineWidth(2);
    // gdNdyvsMassBaryons->SetMarkerSize(1.7);
    // gdNdyvsMassBaryons->SetFillStyle(0);
    // TCanvas *cdNdy = new TCanvas("cdNdy", "dN/dy vs mass", 720, 720);
    // SetCanvasStyle(cdNdy, 0.14, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    // gdNdyvsMassMesons->GetXaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    // gdNdyvsMassMesons->GetYaxis()->SetTitle("dN/dy #times 1/(2J+1)");
    // gdNdyvsMassMesons->GetYaxis()->SetTitleOffset(1.3);
    // gdNdyvsMassMesons->SetMinimum(2e-5);
    // gdNdyvsMassMesons->SetMaximum(3e2);
    // gdNdyvsMassMesons->GetXaxis()->SetLimits(0, 1.99);
    // gdNdyvsMassMesons->Draw("A 2");
    // TGraphAsymmErrors *gdNdyvsMassMesonsSys = (TGraphAsymmErrors *)gdNdyvsMassMesons->Clone("gdNdyvsMassMesonsSys");
    // gdNdyvsMassMesonsSys->SetLineWidth(0);
    // gdNdyvsMassMesonsSys->Draw("P SAME");
    // // gdNdyvsMassBaryons->Draw("2 SAME");
    // TGraphAsymmErrors *gdNdyvsMassBaryonsSys = (TGraphAsymmErrors *)gdNdyvsMassBaryons->Clone("gdNdyvsMassBaryonsSys");
    // gdNdyvsMassBaryonsSys->SetLineWidth(0);
    // // gdNdyvsMassBaryonsSys->Draw("P SAME");
    // // Draw the f2(1525) marker and error bars
    // //pol1 fit for K=1 mesons i.e. 1 strange quarks
    // double K1_masses[3] = { 0.49367, 0.49761, 0.89166}; // Kaon, K0s, K*0
    // double K1_yield[3] = {dNdyAll[1], dNdyAll[2], dNdyAll[3]};
    // TGraph *graph_k1 = new TGraph(3, K1_masses, K1_yield);
    // TF1 *pol1_k1 = new TF1("pol1_k1", "expo", 0.3, 1.4);
    // pol1_k1->SetLineStyle(2);
    // pol1_k1->SetLineColor(kYellow + 2);
    // graph_k1->Fit(pol1_k1, "R");
    // // Draw the same line in the current canvas
    // pol1_k1->Draw("SAME");
    // double f2_dNdy = houtf2->GetBinContent(1) / (2 * 2 + 1); // J = 2 for f2(1525)
    // double f2_dNdy_statErr = houtf2->GetBinContent(2) / (2 * 2 + 1);
    // double f2_dNdy_sysErrLow = houtf2->GetBinContent(3) / (2 * 2 + 1);
    // double f2_dNdy_sysErrHigh = houtf2->GetBinContent(4) / (2 * 2 + 1);
    // cout << "dN/dy f2(1525) " << f2_dNdy << " with stat error " << f2_dNdy_statErr << " sys err low " << f2_dNdy_sysErrLow << " sys err high " << f2_dNdy_sysErrHigh << endl;
    // TGraphErrors *graph_f2_dNdy = new TGraphErrors(1);
    // graph_f2_dNdy->SetPoint(0, f2_mass, f2_dNdy);
    // graph_f2_dNdy->SetPointError(0, 0, f2_dNdy_statErr);
    // graph_f2_dNdy->SetMarkerStyle(f2_marker);
    // graph_f2_dNdy->SetMarkerColor(f2_color);
    // graph_f2_dNdy->SetMarkerSize(1.7);
    // graph_f2_dNdy->SetLineColor(f2_color);
    // graph_f2_dNdy->SetLineWidth(2);
    // graph_f2_dNdy->Draw("P SAME");
    // TGraphAsymmErrors *graph_f2_dNdy_sys = new TGraphAsymmErrors(1);
    // graph_f2_dNdy_sys->SetPoint(0, f2_mass, f2_dNdy);
    // graph_f2_dNdy_sys->SetPointError(0, 0.02, 0.02, f2_dNdy_sysErrLow, f2_dNdy_sysErrHigh);
    // graph_f2_dNdy_sys->SetLineColor(f2_color);
    // graph_f2_dNdy_sys->SetLineWidth(2);
    // graph_f2_dNdy_sys->SetFillStyle(0);
    // graph_f2_dNdy_sys->Draw("2 SAME");

    // // Draw the f2(1525) marker and error bars for the variation case
    // double f2_dNdy_noextrapol = houtf2_noextrapol->GetBinContent(1) / (2 * 2 + 1); // J = 2 for f2(1525)
    // double f2_dNdy_statErr_noextrapol = houtf2_noextrapol->GetBinContent(2) / (2 * 2 + 1);
    // double f2_dNdy_sysErrLow_noextrapol = houtf2_noextrapol->GetBinContent(3) / (2 * 2 + 1);
    // double f2_dNdy_sysErrHigh_noextrapol = houtf2_noextrapol->GetBinContent(4) / (2 * 2 + 1);
    // cout << "dN/dy f2(1525) in toy model " << f2_dNdy_noextrapol << " with stat error " << f2_dNdy_statErr_noextrapol << " sys err low " << f2_dNdy_sysErrLow_noextrapol << " sys err high " << f2_dNdy_sysErrHigh_noextrapol << endl;
    // TGraphErrors *graph_f2_dNdy_noextrapol = new TGraphErrors(1);
    // graph_f2_dNdy_noextrapol->SetPoint(0, f2_mass, f2_dNdy_noextrapol);
    // graph_f2_dNdy_noextrapol->SetPointError(0, 0, f2_dNdy_statErr_noextrapol);
    // graph_f2_dNdy_noextrapol->SetMarkerStyle(23);
    // graph_f2_dNdy_noextrapol->SetMarkerColor(kYellow + 3);
    // graph_f2_dNdy_noextrapol->SetMarkerSize(1.7);
    // graph_f2_dNdy_noextrapol->SetLineColor(kYellow + 3);
    // graph_f2_dNdy_noextrapol->SetLineWidth(2);
    // graph_f2_dNdy_noextrapol->Draw("P SAME");
    // TGraphAsymmErrors *graph_f2_dNdy_sys_noextrapol = new TGraphAsymmErrors(1);
    // graph_f2_dNdy_sys_noextrapol->SetPoint(0, f2_mass, f2_dNdy_noextrapol);
    // graph_f2_dNdy_sys_noextrapol->SetPointError(0, 0.02, 0.02, f2_dNdy_sysErrLow_noextrapol, f2_dNdy_sysErrHigh_noextrapol);
    // graph_f2_dNdy_sys_noextrapol->SetLineColor(kYellow + 3);
    // graph_f2_dNdy_sys_noextrapol->SetLineWidth(2);
    // graph_f2_dNdy_sys_noextrapol->SetFillStyle(0);
    // graph_f2_dNdy_sys_noextrapol->Draw("2 SAME");

    // // Draw the f0(1710) marker and error bars
    // double f0_dNdy = houtf0->GetBinContent(1) / (2 * 0 + 1); // J = 0 for f0(1710)
    // double f0_dNdy_statErr = houtf0->GetBinContent(2) / (2 * 0 + 1);
    // double f0_dNdy_sysErrLow = houtf0->GetBinContent(3) / (2 * 0 + 1);
    // double f0_dNdy_sysErrHigh = houtf0->GetBinContent(4) / (2 * 0 + 1);
    // // cout << "dN/dy f0(1710) " << f0_dNdy << " with stat error " << f0_dNdy_statErr << " sys err low " << f0_dNdy_sysErrLow << " sys err high " << f0_dNdy_sysErrHigh << endl;
    // TGraphErrors *graph_f0_dNdy = new TGraphErrors(1);
    // graph_f0_dNdy->SetPoint(0, f0_mass, f0_dNdy);
    // graph_f0_dNdy->SetPointError(0, 0, f0_dNdy_statErr);
    // graph_f0_dNdy->SetMarkerStyle(f0_marker);
    // graph_f0_dNdy->SetMarkerColor(f0_color);
    // graph_f0_dNdy->SetMarkerSize(1.7);
    // graph_f0_dNdy->SetLineColor(f0_color);
    // graph_f0_dNdy->SetLineWidth(2);
    // graph_f0_dNdy->Draw("P SAME");
    // TGraphAsymmErrors *graph_f0_dNdy_sys = new TGraphAsymmErrors(1);
    // graph_f0_dNdy_sys->SetPoint(0, f0_mass, f0_dNdy);
    // graph_f0_dNdy_sys->SetPointError(0, 0.02, 0.02, f0_dNdy_sysErrLow, f0_dNdy_sysErrHigh);
    // graph_f0_dNdy_sys->SetLineColor(f0_color);
    // graph_f0_dNdy_sys->SetLineWidth(2);
    // graph_f0_dNdy_sys->SetFillStyle(0);
    // graph_f0_dNdy_sys->Draw("2 SAME");

    // // Draw the f0(1710) marker and error bars for the variation case
    // double f0_dNdy_noextrapol = houtf0_noextrapol->GetBinContent(1) / (2 * 0 + 1); // J = 0 for f0(1710)
    // double f0_dNdy_statErr_noextrapol = houtf0_noextrapol->GetBinContent(2) / (2 * 0 + 1);
    // double f0_dNdy_sysErrLow_noextrapol = houtf0_noextrapol->GetBinContent(3) / (2 * 0 + 1);
    // double f0_dNdy_sysErrHigh_noextrapol = houtf0_noextrapol->GetBinContent(4) / (2 * 0 + 1);
    // cout << "dN/dy f0(1710) in toy model " << f0_dNdy_noextrapol << " with stat error " << f0_dNdy_statErr_noextrapol << " sys err low " << f0_dNdy_sysErrLow_noextrapol << " sys err high " << f0_dNdy_sysErrHigh_noextrapol << endl;
    // TGraphErrors *graph_f0_dNdy_noextrapol = new TGraphErrors(1);
    // graph_f0_dNdy_noextrapol->SetPoint(0, f0_mass, f0_dNdy_noextrapol);
    // graph_f0_dNdy_noextrapol->SetPointError(0, 0, f0_dNdy_statErr_noextrapol);
    // graph_f0_dNdy_noextrapol->SetMarkerStyle(23);
    // graph_f0_dNdy_noextrapol->SetMarkerColor(kCyan + 1);
    // graph_f0_dNdy_noextrapol->SetMarkerSize(1.7);
    // graph_f0_dNdy_noextrapol->SetLineColor(kCyan + 1);
    // graph_f0_dNdy_noextrapol->SetLineWidth(2);
    // graph_f0_dNdy_noextrapol->Draw("P SAME");
    // TGraphAsymmErrors *graph_f0_dNdy_sys_noextrapol = new TGraphAsymmErrors(1);
    // graph_f0_dNdy_sys_noextrapol->SetPoint(0, f0_mass, f0_dNdy_noextrapol);
    // graph_f0_dNdy_sys_noextrapol->SetPointError(0, 0.02, 0.02, f0_dNdy_sysErrLow_noextrapol, f0_dNdy_sysErrHigh_noextrapol);
    // graph_f0_dNdy_sys_noextrapol->SetLineColor(kCyan + 1);
    // graph_f0_dNdy_sys_noextrapol->SetLineWidth(2);
    // graph_f0_dNdy_sys_noextrapol->SetFillStyle(0);
    // graph_f0_dNdy_sys_noextrapol->Draw("2 SAME");

    // // Add particle names below each point
    // TLatex latex2;
    // latex2.SetTextAlign(22);
    // latex2.SetTextSize(0.035);
    // for (int i = 0; i < totalParticles; i++)
    // {
    //     double x = particleMass[i];
    //     double y = dNdyvalues_13TeV[i][0] * dNdyvalues_13TeV[i][3] / (2 * Jvalues[i] + 1); // apply the scaling factor
    //     // latex2.SetTextColor(colors[i]);
    //     if (i == 0)
    //         latex2.DrawLatex(x, y + 1.5, particlesLatex[i].c_str());
    //     else if (i == 1)
    //         latex2.DrawLatex(x, y + 0.4, particlesLatex[i].c_str());
    //     else if (i == 2)
    //         latex2.DrawLatex(x + 0.08, y + 0.4, particlesLatex[i].c_str());
    //     else if (i == 3)
    //         latex2.DrawLatex(x - 0.02, y + 0.01, particlesLatex[i].c_str());
    //     else if (i == 4)
    //         latex2.DrawLatex(x, y + 0.01, particlesLatex[i].c_str());
    //     // else if (i == 5)
    //     //     latex2.DrawLatex(x, y + 0.15, particlesLatex[i].c_str());
    //     // else if (i == 6)
    //     //     latex2.DrawLatex(x, y + 0.12, particlesLatex[i].c_str());
    //     // else if (i == 7)
    //     //     latex2.DrawLatex(x, y + 0.02, particlesLatex[i].c_str());
    //     // else if (i == 8)
    //     //     latex2.DrawLatex(x - 0.015, y + 0.001, particlesLatex[i].c_str());
    // }
    // latex2.DrawLatex(1.5173, houtf2->GetBinContent(1) - 0.0045, "f'_{2}");
    // latex2.DrawLatex(1.710 + 0.015, houtf0->GetBinContent(1) + 0.001, "f_{0}");

    // TLegend *legend6 = new TLegend(0.3, 0.8, 0.95, 0.92);
    // // legend6->SetBorderSize(0);
    // legend6->SetFillStyle(0);
    // legend6->SetNColumns(2);
    // legend6->SetTextSize(0.03);
    // legend6->AddEntry(gdNdyvsMassMesons, "Mesons (13 TeV)", "p");
    // legend6->AddEntry(gdNdyvsMassBaryons, "Baryons (13 TeV)", "p");
    // legend6->AddEntry(graph_f2_dNdy, "f'_{2}(1525) (MC)", "p");
    // legend6->AddEntry(graph_f2_dNdy_noextrapol, "f'_{2}(1525) (Toy Model)", "p");
    // legend6->AddEntry(graph_f0_dNdy, "f_{0}(1710) (MC)", "p");
    // legend6->AddEntry(graph_f0_dNdy_noextrapol, "f_{0}(1710) (Toy Model)", "p");
    // legend6->Draw();

    // TCanvas *cdNdy_k1 = new TCanvas("cdNdy_k1", "dN/dy vs mass with K=1 fit", 720, 720);
    // SetCanvasStyle(cdNdy_k1, 0.14, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    // graph_k1->GetXaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    // graph_k1->GetYaxis()->SetTitle("dN/dy #times 1/(2J+1)");
    // graph_k1->GetYaxis()->SetTitleOffset(1.3);
    // graph_k1->SetMarkerStyle(20);
    // graph_k1->SetMarkerSize(1.7);
    // graph_k1->Draw("AP");
    // pol1_k1->Draw("SAME");

    TCanvas *cCorrectedFitf0Toy = new TCanvas("cCorrectedFitf0Toy", "Corrected Fit f0(1710) Toy", 720, 720);
    SetCanvasStyle(cCorrectedFitf0Toy, 0.14, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    hf03->SetMarkerStyle(24);
    hf03->SetMarkerColor(kBlue);
    hf03->SetLineColor(kBlue);
    hf03->SetMarkerSize(1.6);
    hf03->SetMaximum(9e-4);
    hf03->SetMinimum(3e-8);
    hf03->Draw("E1");
    fitFcnf0_2->SetLineColor(kRed);
    fitFcnf0_2->SetLineWidth(2);
    fitFcnf0_2->Draw("SAME");
    TLegend *legend6 = new TLegend(0.19, 0.25, 0.4, 0.5);
    legend6->SetBorderSize(0);
    legend6->SetFillStyle(0);
    legend6->SetTextSize(0.032);
    // legend6->AddEntry((TObject *)0, "ALICE", "");
    legend6->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legend6->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    // legend6->AddEntry(hf03, "f_{0}(1710) (Toy Model)", "p");
    legend6->AddEntry(hf03, "f_{0}(1710) (1rBW Fit)", "p");
    legend6->AddEntry(fitFcnf0_2, "Levy-Tsallis fit", "l");
    legend6->Draw();

    TPaveText *box = new TPaveText(0.60, 0.7, 0.98, 0.93, "NDC");
    box->SetFillColor(0); // white background
    box->SetBorderSize(1);
    box->SetTextAlign(12); // left align
    box->SetTextSize(0.03);
    box->SetTextFont(42);
    box->AddText(Form("#chi^{2} / ndf                      %.3f / %d", fitFcnf0_2->GetChisquare(), fitFcnf0_2->GetNDF()));
    box->AddText(Form("n                        %.3f #pm %.3f", fitFcnf0_2->GetParameter(0), fitFcnf0_2->GetParError(0)));
    box->AddText(Form("dN/dy        %.5f #pm %.5f", fitFcnf0_2->GetParameter(1), fitFcnf0_2->GetParError(1)));
    box->AddText(Form("T                        %.3f #pm %.3f", fitFcnf0_2->GetParameter(3), fitFcnf0_2->GetParError(3)));
    box->Draw();
    // cCorrectedFitf0Toy->SaveAs((savePath + "/plots/CorrectedFit_f0Toy.png").c_str());

    TCanvas *cCorrectedFitf2Toy = new TCanvas("cCorrectedFitf2Toy", "Corrected Fit f2(1525) Toy", 720, 720);
    SetCanvasStyle(cCorrectedFitf2Toy, 0.14, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    hf23->SetMarkerStyle(24);
    hf23->SetMarkerColor(kBlue);
    hf23->SetLineColor(kBlue);
    hf23->SetMarkerSize(1.6);
    hf23->SetMaximum(9e-2);
    hf23->SetMinimum(7e-8);
    hf23->Draw("E1");
    fitFcnf2_2->SetLineColor(kRed);
    fitFcnf2_2->SetLineWidth(2);
    fitFcnf2_2->Draw("SAME");
    TLegend *legend7 = new TLegend(0.19, 0.25, 0.4, 0.5);
    legend7->SetBorderSize(0);
    legend7->SetFillStyle(0);
    legend7->SetTextSize(0.032);
    // legend7->AddEntry((TObject *)0, "ALICE", "");
    legend7->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legend7->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    legend7->AddEntry(hf23, "f'_{2}(1525) (1rBW Fit)", "p");
    // legend7->AddEntry(hf23, "f'_{2}(1525) (Toy Model)", "p");
    legend7->AddEntry(fitFcnf2_2, "Levy-Tsallis fit", "l");
    legend7->Draw();

    TPaveText *box2 = new TPaveText(0.60, 0.7, 0.98, 0.93, "NDC");
    box2->SetFillColor(0); // white background
    box2->SetBorderSize(1);
    box2->SetTextAlign(12); // left align
    box2->SetTextSize(0.03);
    box2->SetTextFont(42);
    box2->AddText(Form("#chi^{2} / ndf                      %.3f / %d", fitFcnf2_2->GetChisquare(), fitFcnf2_2->GetNDF()));
    box2->AddText(Form("n                        %.3f #pm %.3f", fitFcnf2_2->GetParameter(0), fitFcnf2_2->GetParError(0)));
    box2->AddText(Form("dN/dy        %.5f #pm %.5f", fitFcnf2_2->GetParameter(1), fitFcnf2_2->GetParError(1)));
    box2->AddText(Form("T                        %.3f #pm %.3f", fitFcnf2_2->GetParameter(3), fitFcnf2_2->GetParError(3)));
    box2->Draw();
    cCorrectedFitf2Toy->SaveAs((savePath + "/plots/CorrectedFit_f21rBWFit.png").c_str());
    // cCorrectedFitf2Toy->SaveAs((savePath + "/plots/CorrectedFit_f2Toy.png").c_str());

    /*
    TCanvas *cCorrectedFitf0Default = new TCanvas("cCorrectedFitf0Default", "Corrected Fit f0(1710) Toy", 720, 720);
    SetCanvasStyle(cCorrectedFitf0Default, 0.14, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    hf01->SetMarkerStyle(24);
    hf01->SetMarkerColor(kBlue);
    hf01->SetLineColor(kBlue);
    hf01->SetMarkerSize(1.6);
    hf01->SetMaximum(9e-4);
    hf01->SetMinimum(3e-8);
    hf01->Draw("E1");
    fitFcnf0->SetLineColor(kRed);
    fitFcnf0->SetLineWidth(2);
    fitFcnf0->Draw("SAME");
    TLegend *legend8 = new TLegend(0.19, 0.25, 0.4, 0.5);
    legend8->SetBorderSize(0);
    legend8->SetFillStyle(0);
    legend8->SetTextSize(0.032);
    // legend8->AddEntry((TObject *)0, "ALICE", "");
    legend8->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legend8->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    legend8->AddEntry(hf01, "f_{0}(1710) (Default)", "p");
    legend8->AddEntry(fitFcnf0, "Levy-Tsallis fit", "l");
    legend8->Draw();

    TPaveText *box3 = new TPaveText(0.60, 0.7, 0.98, 0.93, "NDC");
    box3->SetFillColor(0); // white background
    box3->SetBorderSize(1);
    box3->SetTextAlign(12); // left align
    box3->SetTextSize(0.03);
    box3->SetTextFont(42);
    box3->AddText(Form("#chi^{2} / ndf                      %.3f / %d", fitFcnf0->GetChisquare(), fitFcnf0->GetNDF()));
    box3->AddText(Form("n                        %.3f #pm %.3f", fitFcnf0->GetParameter(0), fitFcnf0->GetParError(0)));
    box3->AddText(Form("dN/dy        %.5f #pm %.5f", fitFcnf0->GetParameter(1), fitFcnf0->GetParError(1)));
    box3->AddText(Form("T                        %.3f #pm %.3f", fitFcnf0->GetParameter(3), fitFcnf0->GetParError(3)));
    box3->Draw();
    // cCorrectedFitf0Default->SaveAs((savePath + "/plots/CorrectedFit_f0Default.png").c_str());

    TCanvas *cCorrectedFitf2Default = new TCanvas("cCorrectedFitf2Default", "Corrected Fit f2(1525) Toy", 720, 720);
    SetCanvasStyle(cCorrectedFitf2Default, 0.14, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    hf21->SetMarkerStyle(24);
    hf21->SetMarkerColor(kBlue);
    hf21->SetLineColor(kBlue);
    hf21->SetMarkerSize(1.6);
    hf21->SetMaximum(9e-3);
    hf21->SetMinimum(7e-8);
    hf21->Draw("E1");
    fitFcnf2->SetLineColor(kRed);
    fitFcnf2->SetLineWidth(2);
    fitFcnf2->Draw("SAME");
    TLegend *legend9 = new TLegend(0.19, 0.25, 0.4, 0.5);
    legend9->SetBorderSize(0);
    legend9->SetFillStyle(0);
    legend9->SetTextSize(0.032);
    // legend9->AddEntry((TObject *)0, "ALICE", "");
    legend9->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legend9->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    legend9->AddEntry(hf21, "f'_{2}(1525) (Default)", "p");
    legend9->AddEntry(fitFcnf2, "Levy-Tsallis fit", "l");
    legend9->Draw();

    TPaveText *box4 = new TPaveText(0.60, 0.7, 0.98, 0.93, "NDC");
    box4->SetFillColor(0); // white background
    box4->SetBorderSize(1);
    box4->SetTextAlign(12); // left align
    box4->SetTextSize(0.03);
    box4->SetTextFont(42);
    box4->AddText(Form("#chi^{2} / ndf                      %.3f / %d", fitFcnf2->GetChisquare(), fitFcnf2->GetNDF()));
    box4->AddText(Form("n                        %.3f #pm %.3f", fitFcnf2->GetParameter(0), fitFcnf2->GetParError(0)));
    box4->AddText(Form("dN/dy        %.5f #pm %.5f", fitFcnf2->GetParameter(1), fitFcnf2->GetParError(1)));
    box4->AddText(Form("T                        %.3f #pm %.3f", fitFcnf2->GetParameter(3), fitFcnf2->GetParError(3)));
    box4->Draw();
    // cCorrectedFitf2Default->SaveAs((savePath + "/plots/CorrectedFit_f2Default.png").c_str());
    */
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