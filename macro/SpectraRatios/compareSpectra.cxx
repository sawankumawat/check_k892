#include <iostream>
#include <iomanip>
#include "../src/style.h"
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

TH1D *GetHisto(TFile *f, const string &name)
{
    TH1D *histo = (TH1D *)f->Get(name.c_str());

    if (!histo || histo == nullptr)
    {
        cout << "Error: histo " << name << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetHistoQA(histo);
    histo->SetTitle(0);
    return histo;
}

TGraphErrors *GetGraph(TFile *f, const string &name)
{
    TGraphErrors *graph = (TGraphErrors *)f->Get(name.c_str());

    if (!graph || graph == nullptr)
    {
        cout << "Error: graph " << name << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetGraphErrorStyle(graph);
    graph->SetTitle(0);
    return graph;
}

TH1D *RebinToMatch(const TH1D *hInput, const TH1D *hReference, const TString &newName = "")
{
    if (!hInput || !hReference)
        return nullptr;

    int nBins = hReference->GetNbinsX();

    std::vector<double> bins(nBins + 1);

    for (int i = 1; i <= nBins + 1; i++)
        bins[i - 1] = hReference->GetBinLowEdge(i);

    TString hName = newName;
    if (hName.IsNull())
        hName = Form("%s_rebinned", hInput->GetName());

    TH1D *hClone = (TH1D *)hInput->Clone(Form("%s_tmp", hInput->GetName()));

    return (TH1D *)hClone->Rebin(nBins, hName.Data(), bins.data());
}

void compareSpectra()
{
    TFile *fSpectras[11];
    int multClasses[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    TH1D *hSpectra[11];
    for (int i = 0; i < 11; i++)
    {
        double multLow, multHigh;
        if (i == 0)
        {
            multLow = 0;
            multHigh = 100;
        }
        else
        {
            multLow = multClasses[i - 1];
            multHigh = multClasses[i];
        }
        fSpectras[i] = OpenFile(Form("../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/corrected_spectra_%.0f_%.0f.root", multLow, multHigh));
        hSpectra[i] = GetHisto(fSpectras[i], Form("mult_%.0f-%.0f/corrected_spectra_Integral_final", multLow, multHigh));
    }

    //===============================================
    // ===========Pythia tested locally================
    //===============================================
    TFile *fPythiaMonash = new TFile("../../pythia/Pythia_MonashLocal.root", "read");
    TFile *fPythiaMonashNoCR = new TFile("../../pythia/Pythia_MonashWoCRLocal.root", "read");
    TFile *fPythiaShoving = new TFile("../../pythia/Pythia_ShovingLocal.root", "read");
    TFile *fPythiaMonashRescattering = new TFile("../../pythia/Pythia_MonashRescatteringLocal.root", "read");
    TFile *fPythiaRopes = new TFile("../../pythia/Pythia_RopesLocal.root", "read");
    if (fPythiaMonash->IsZombie() || fPythiaMonashNoCR->IsZombie() || fPythiaShoving->IsZombie() || fPythiaMonashRescattering->IsZombie() || fPythiaRopes->IsZombie())
    {
        cout << "Error: Pythia local files not found" << endl;
        return;
    }
    enum PythiaModel
    {
        kPythiaMonashLocal,
        kPythiaMonashNoCRLocal,
        kPythiaShovingLocal,
        kPythiaRopesLocal,
        kNPythiaModels
    };

    const char *modelLabelLocal[kNPythiaModels] = {
        "Pythia Monash",
        "Pythia Monash No CR",
        "Pythia Shoving",
        "Pythia Ropes"};
    vector<TFile *> fPythiaModels = {fPythiaMonash, fPythiaMonashNoCR, fPythiaShoving, fPythiaRopes};
    TGraphErrors *gPythiaYieldLocal[kNPythiaModels];
    double dnch_detaRun3[] = {21.78, 18.48, 15.76, 13.89, 12.50, 10.86, 9.09, 7.63, 5.87, 3.69};

    TH1D *hSpectraHighestMult[kNPythiaModels];
    TH1D *hSpectraMult20to30[kNPythiaModels];

    int colorsPythia[kNPythiaModels] = {kCyan + 1, kYellow + 1, kMagenta, kGreen + 2};

    for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    {
        gPythiaYieldLocal[imodel] = GetGraph(fPythiaModels[imodel], "gMeanYield_kstar");

        int totalPoints = gPythiaYieldLocal[imodel]->GetN();
        cout << "Model: " << modelLabelLocal[imodel] << ", Total points: " << totalPoints << endl;
        for (int ipoint = 0; ipoint < totalPoints; ipoint++)
        {
            double x, y;
            gPythiaYieldLocal[imodel]->GetPoint(ipoint, x, y);
            cout << "Point " << ipoint << ": Multiplicity " << x << ", Yield: " << y << endl;
        }
        cout << endl;
    }

    // Monash tune
    TH1D *hSpectraTemp1Monash = GetHisto(fPythiaMonash, "hkstarPt_3"); // 3
    hSpectraTemp1Monash->Scale(1.0 / hSpectraTemp1Monash->Integral());
    TH1D *hSpectraTemp2Monash = GetHisto(fPythiaMonash, "hkstarPt_4"); // 4
    hSpectraTemp2Monash->Scale(1.0 / hSpectraTemp2Monash->Integral());
    hSpectraHighestMult[kPythiaMonashLocal] = (TH1D *)hSpectraTemp1Monash->Clone("hSpectraHighestMultMonash");
    hSpectraHighestMult[kPythiaMonashLocal]->Add(hSpectraTemp2Monash);
    hSpectraHighestMult[kPythiaMonashLocal]->Scale(0.5);

    hSpectraMult20to30[kPythiaMonashLocal] = GetHisto(fPythiaMonash, "hkstarPt_10");
    hSpectraMult20to30[kPythiaMonashLocal]->Scale(2.0 / hSpectraMult20to30[kPythiaMonashLocal]->Integral());

    // Shoving
    TH1D *hSpectraTemp1Shoving = GetHisto(fPythiaShoving, "hkstarPt_0"); // 3
    hSpectraTemp1Shoving->Scale(1.0 / hSpectraTemp1Shoving->Integral());
    TH1D *hSpectraTemp2Shoving = GetHisto(fPythiaShoving, "hkstarPt_1"); // 4
    hSpectraTemp2Shoving->Scale(1.0 / hSpectraTemp2Shoving->Integral());
    hSpectraHighestMult[kPythiaShovingLocal] = (TH1D *)hSpectraTemp1Shoving->Clone("hSpectraHighestMultShoving");
    hSpectraHighestMult[kPythiaShovingLocal]->Add(hSpectraTemp2Shoving);
    hSpectraHighestMult[kPythiaShovingLocal]->Scale(2.0);

    hSpectraMult20to30[kPythiaShovingLocal] = GetHisto(fPythiaShoving, "hkstarPt_11");
    hSpectraMult20to30[kPythiaShovingLocal]->Scale(2.0 / hSpectraMult20to30[kPythiaShovingLocal]->Integral());

    // Ropes
    TH1D *hSpectraTemp1Ropes = GetHisto(fPythiaRopes, "hkstarPt_0"); // 3
    hSpectraTemp1Ropes->Scale(1.0 / hSpectraTemp1Ropes->Integral());
    TH1D *hSpectraTemp2Ropes = GetHisto(fPythiaRopes, "hkstarPt_1"); // 4
    hSpectraTemp2Ropes->Scale(1.0 / hSpectraTemp2Ropes->Integral());
    hSpectraHighestMult[kPythiaRopesLocal] = (TH1D *)hSpectraTemp1Ropes->Clone("hSpectraHighestMultRopes");
    hSpectraHighestMult[kPythiaRopesLocal]->Add(hSpectraTemp2Ropes);
    hSpectraHighestMult[kPythiaRopesLocal]->Scale(2.0);

    hSpectraMult20to30[kPythiaRopesLocal] = GetHisto(fPythiaRopes, "hkstarPt_9");
    hSpectraMult20to30[kPythiaRopesLocal]->Scale(2.0 / hSpectraMult20to30[kPythiaRopesLocal]->Integral());

    // Monash without CR
    TH1D *hSpectraTemp1MonashNoCR = GetHisto(fPythiaMonashNoCR, "hkstarPt_0"); // 3
    hSpectraTemp1MonashNoCR->Scale(1.0 / hSpectraTemp1MonashNoCR->Integral());
    TH1D *hSpectraTemp2MonashNoCR = GetHisto(fPythiaMonashNoCR, "hkstarPt_1"); // 4
    hSpectraTemp2MonashNoCR->Scale(1.0 / hSpectraTemp2MonashNoCR->Integral());
    hSpectraHighestMult[kPythiaMonashNoCRLocal] = (TH1D *)hSpectraTemp1MonashNoCR->Clone("hSpectraHighestMultMonashNoCR");
    hSpectraHighestMult[kPythiaMonashNoCRLocal]->Add(hSpectraTemp2MonashNoCR);
    hSpectraHighestMult[kPythiaMonashNoCRLocal]->Scale(2.0);

    hSpectraMult20to30[kPythiaMonashNoCRLocal] = GetHisto(fPythiaMonashNoCR, "hkstarPt_11");
    hSpectraMult20to30[kPythiaMonashNoCRLocal]->Scale(2.0 / hSpectraMult20to30[kPythiaMonashNoCRLocal]->Integral());

    TCanvas *cSpectraHighestMult = new TCanvas("cSpectraHighestMult", "cSpectraHighestMult", 720, 720);
    SetCanvasStyle(cSpectraHighestMult, 0.15, 0.03, 0.05, 0.15);
    gPad->SetLogy();
    hSpectra[1]->SetMaximum(5 * hSpectra[1]->GetMaximum());
    hSpectra[1]->SetMinimum(9e-7);
    hSpectra[1]->Draw("pe"); // Mult 0-1, dNch/deta = 21.78
    // hSpectraHighestMult[kPythiaMonashLocal]->SetLineColor(kRed + 1);
    // hSpectraHighestMult[kPythiaMonashLocal]->Draw("HIST same");
    // hSpectraHighestMult[kPythiaShovingLocal]->SetLineColor(kBlue + 1);
    // hSpectraHighestMult[kPythiaShovingLocal]->Draw("HIST same");
    // hSpectraHighestMult[kPythiaRopesLocal]->SetLineColor(kGreen + 1);
    // hSpectraHighestMult[kPythiaRopesLocal]->SetLineStyle(2);
    // hSpectraHighestMult[kPythiaRopesLocal]->Draw("HIST same");
    // hSpectraHighestMult[kPythiaMonashNoCRLocal]->SetLineColor(kMagenta + 1);
    // hSpectraHighestMult[kPythiaMonashNoCRLocal]->SetLineStyle(2);
    // hSpectraHighestMult[kPythiaMonashNoCRLocal]->Draw("HIST same");
    TH1D *hSpectraHighestMultRebin = RebinToMatch(hSpectraHighestMult[kPythiaMonashLocal], hSpectra[1], "hSpectraHighestMultRebin");
    hSpectraHighestMultRebin->SetLineColor(kRed + 1);
    hSpectraHighestMultRebin->Draw("HIST same");

    // TLegend *legSpectraHighestMult = new TLegend(0.6, 0.5, 0.92, 0.92);
    // SetLegendStyle(legSpectraHighestMult);
    // legSpectraHighestMult->SetTextSize(0.03);
    // legSpectraHighestMult->SetHeader("Multiplicity 0-1%");
    // legSpectraHighestMult->AddEntry(hSpectra[1], "Data", "lpe");
    // legSpectraHighestMult->AddEntry(hSpectraHighestMult[kPythiaMonashLocal], "Pythia Monash", "l");
    // legSpectraHighestMult->AddEntry(hSpectraHighestMult[kPythiaMonashNoCRLocal], "Pythia Monash No CR", "l");
    // legSpectraHighestMult->AddEntry(hSpectraHighestMult[kPythiaShovingLocal], "Pythia Shoving", "l");
    // legSpectraHighestMult->AddEntry(hSpectraHighestMult[kPythiaRopesLocal], "Pythia Ropes", "l");
    // legSpectraHighestMult->Draw();
    // cSpectraHighestMult->SaveAs("Plots/spectra_highestMult_comparison.png");

    // TCanvas *cSpectra2 = new TCanvas("cSpectra2", "cSpectra2", 720, 720);
    // SetCanvasStyle(cSpectra2, 0.15, 0.03, 0.05, 0.15);
    // gPad->SetLogy();
    // hSpectra[6]->SetMaximum(5 * hSpectra[6]->GetMaximum());
    // hSpectra[6]->SetMinimum(9e-7);
    // hSpectra[6]->Draw("pe"); // Mult 20-30, dNch/deta = 10.86

    // hSpectraMult20to30[kPythiaMonashLocal]->SetLineColor(kRed + 1);
    // hSpectraMult20to30[kPythiaMonashLocal]->Draw("HIST same");
    // hSpectraMult20to30[kPythiaShovingLocal]->SetLineColor(kBlue + 1);
    // hSpectraMult20to30[kPythiaShovingLocal]->Draw("HIST same");
    // hSpectraMult20to30[kPythiaRopesLocal]->SetLineColor(kGreen + 1);
    // hSpectraMult20to30[kPythiaRopesLocal]->SetLineStyle(2);
    // hSpectraMult20to30[kPythiaRopesLocal]->Draw("HIST same");
    // hSpectraMult20to30[kPythiaMonashNoCRLocal]->SetLineColor(kMagenta + 1);
    // hSpectraMult20to30[kPythiaMonashNoCRLocal]->SetLineStyle(2);
    // hSpectraMult20to30[kPythiaMonashNoCRLocal]->Draw("HIST same");

    // TLegend *legClone1 = (TLegend *)legSpectraHighestMult->Clone("legClone1");
    // legClone1->SetHeader("Multiplicity 20-30%");
    // legClone1->Draw();
    // cSpectra2->SaveAs("Plots/spectra_mult20to30_comparison.png");

    // //======================================================
    // //    ===========EPOS local model files===========
    // //======================================================
    // TFile *fEPOS = OpenFile("ModelRootFiles/EPOS_finalQA_INELgt0Correct.root");
    // TH1D *hEPOS_Yield = GetHisto(fEPOS, "IST9_ITY80/hPtMB_kstar");
    // TCanvas *cEPOS = new TCanvas("cEPOS", "cEPOS", 720, 720);
    // SetCanvasStyle(cEPOS, 0.15, 0.03, 0.05, 0.15);
    // gPad->SetLogy();
    // hSpectra[0]->Draw("PE");

    // hEPOS_Yield->Scale(1.0 / hEPOS_Yield->Integral());
    // hEPOS_Yield->SetLineColor(kRed + 1);

    // TH1D *hEPOS_rebinned = RebinToMatch(hEPOS_Yield, hSpectra[0], "hEPOS_rebinned");
    // hEPOS_rebinned->Draw("HIST SAME");
}
