#include <iostream>
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

    TString hName = newName;
    if (hName.IsNull())
        hName = Form("%s_rebinned", hInput->GetName());

    TH1D *hRebinned = (TH1D *)hReference->Clone(hName);
    hRebinned->Reset();
    // hRebinned->Sumw2();

    for (int i = 1; i <= hReference->GetNbinsX(); ++i)
    {
        double low = hReference->GetBinLowEdge(i);
        double high = hReference->GetBinLowEdge(i + 1);

        int binLow = hInput->GetXaxis()->FindBin(low + 1e-8);
        int binHigh = hInput->GetXaxis()->FindBin(high - 1e-8);

        double err = 0.0;
        double integral = hInput->IntegralAndError(binLow, binHigh, err);

        double width = high - low;
        // cout << "Bin " << i << ", pT Low " << low << " High " << high << ", width " << width << endl;

        hRebinned->SetBinContent(i, integral / width);
        hRebinned->SetBinError(i, err / width);
    }

    return hRebinned;
}


void checkEPOS()
{
    TFile *fEPOS = OpenFile("EPOS_finalQA_CorrectpTCutPiKp.root");
    TGraph *gEPOSdNdy = GetGraph(fEPOS, "IST9_ITY80/kstar_vs_mult");
    const int totalCent = gEPOSdNdy->GetN();
    // cout << "Total centrality points in EPOS: " << totalCent << endl;
    TH1D *hEPOSYields[totalCent];

    for (int i = 0; i < totalCent; i++)
    {
        double mult, yield;
        gEPOSdNdy->GetPoint(i, mult, yield);
        hEPOSYields[i] = GetHisto(fEPOS, Form("IST9_ITY80/hPtMB_kstar_IST9_ITY80_Cent%d", i));
        hEPOSYields[i]->Scale(0.1); // Temporary multiplying by original bin width of 0.1 GeV/c.
        hEPOSYields[i]->Scale(0.5); // Average K*892 and anti-K*892
        cout << "Multiplicity: " << mult << ", EPOS Yield: " << yield / 2.0 << ", EPOS Yield from histogram: " << hEPOSYields[i]->Integral() << endl;
    }
    cout << endl;

    TFile *fData = OpenFile("../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/Results.root");
    TGraphErrors *gdNdyData = GetGraph(fData, "gMeanYieldRun3_sys");
    const int totalPoints = gdNdyData->GetN();
    for (int i = 0; i < totalPoints; i++)
    {
        double mult, yield;
        gdNdyData->GetPoint(i, mult, yield);
        cout << "Multiplicity: " << mult << ", Data Yield: " << yield << endl;
    }

    TFile *fSpectras[11];
    int multClasses[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    TH1D *hSpectra[11];
    TH1D *hSpectraSys[11];
    string sysPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/SystematicsPlots/SysUncert.root";
    TFile *fSystematics = OpenFile(sysPath);
    TH1D *hRelUncert = GetHisto(fSystematics, "hTotalSysSmoothed_0_100");

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
        hSpectraSys[i] = (TH1D *)hSpectra[i]->Clone(Form("hSpectraSys_%.0f-%.0f", multLow, multHigh));
        for (int bin = 1; bin <= hSpectraSys[i]->GetNbinsX(); ++bin)
        {
            double content = hSpectraSys[i]->GetBinContent(bin);
            double relUncert = hRelUncert->GetBinContent(bin);
            double sysError = content * relUncert;
            hSpectraSys[i]->SetBinError(bin, sysError);
        }
    }

    TCanvas *cSpectra = new TCanvas("cSpectra", "cSpectra", 720, 720);
    SetCanvasStyle(cSpectra, 0.15, 0.03, 0.05, 0.15);
    gPad->SetLogy();
    int whichMult = 10; // 0: Min Bias, 1: 0-1%, 2: 1-5%, 3: 5-10%, 4: 10-15%, 5: 15-20%, 6: 20-30%, 7: 30-40%, 8: 40-50%, 9: 50-70%, 10: 70-100%
    hSpectraSys[whichMult]->SetLineColor(kBlack);
    hSpectraSys[whichMult]->SetLineWidth(3);
    hSpectraSys[whichMult]->SetMarkerStyle(20);
    hSpectraSys[whichMult]->SetMarkerSize(1.2);
    hSpectraSys[whichMult]->Draw("PE");
    TH1D *hEPOSYieldsRebinned = RebinToMatch(hEPOSYields[19], hSpectra[whichMult], "hEPOSYieldsRebinned");
    hEPOSYieldsRebinned->SetLineColor(kRed + 1);
    hEPOSYieldsRebinned->SetLineStyle(2);
    hEPOSYieldsRebinned->SetLineWidth(3);
    hEPOSYieldsRebinned->Draw("HIST SAME");
    cout << "Yield from rebinned EPOS histogram: " << hEPOSYieldsRebinned->Integral("width") << endl;
    cout << "Yield from original EPOS histogram: " << hEPOSYields[19]->Integral() << endl;
}