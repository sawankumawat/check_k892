#include <iostream>
#include <iomanip>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/initializations.h"

void modelPrediction2()
{
    //_id53739 for Pythia Monash with rescattering
    // TFile *fmodel = new TFile("ModelRootFiles/Pythia_Shoving.root", "read");
    TFile *fmodel = new TFile("ModelRootFiles/EPOS.root", "read");
    if (fmodel->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }
    vector<string> particleTypes = {"Kstar", "Phi", "Pion", "Kaon", "Proton"};
    int selectParticle = 0;

    TH2D *hNch05VsFT0M = (TH2D *)fmodel->Get("mc-particle-prediction/multiplicity/vsETA05/FT0AC");
    TH2D *hParticlevsFT0M = (TH2D *)fmodel->Get(Form("mc-particle-prediction/prediction/pt/FT0AC/%s", particleTypes[selectParticle].c_str()));
    if (hNch05VsFT0M == nullptr || hParticlevsFT0M == nullptr)
    {
        cout << "Error: histograms not found" << endl;
        return;
    }

    TH1D *hNFT0M = hNch05VsFT0M->ProjectionX("hNFT0M", -1, -1);
    // Divide into 10 regions of equal number of events
    Double_t totalEntries = hNFT0M->Integral();
    Double_t entriesPerBin = totalEntries / 12.0;
    cout << "Entries per bin: " << entriesPerBin << endl;
    Double_t cumulativeEntries = 0;
    int counter = 0;
    vector<int> binEdges;
    binEdges.push_back(1); // Start with the first bin edge
    vector<double> xValuesFT0Mcorrespondingdnch_detaRun3;
    int PrevBin = 1;

    for (int i = 1; i <= hNFT0M->GetNbinsX(); ++i)
    {
        cumulativeEntries += hNFT0M->GetBinContent(i);
        if (cumulativeEntries >= entriesPerBin)
        {
            cout << " xValue = " << hNch05VsFT0M->GetXaxis()->GetBinCenter(i) << " entries = " << cumulativeEntries << endl;
            TH1D *hNch05 = hNch05VsFT0M->ProjectionY(Form("hNch05_%d", i), PrevBin, i);
            cout << "<dNch/d#eta> " << hNch05->GetMean() << ", total Entries " << hNch05->Integral() << endl;
            if (hNch05->GetMean() < 3.6)
                continue;

            binEdges.push_back(i);
            xValuesFT0Mcorrespondingdnch_detaRun3.push_back(hNch05->GetMean());

            PrevBin = i + 1; // advance to next bin to make ranges non-overlapping
            if (PrevBin > hNFT0M->GetNbinsX())
            {
                // nothing more to do
                break;
            }
            counter++;
            cumulativeEntries = 0; // reset for next bin
        }
    }
    cout << "Total bins with non-zero entries: " << counter << endl;

    TCanvas *cNch05 = new TCanvas("cNch05", "cNch05", 720, 720);
    hNch05VsFT0M->Draw("colz");

    int nGraphPoints = (int)xValuesFT0Mcorrespondingdnch_detaRun3.size();
    TGraphErrors *gMeanpTvsNch = new TGraphErrors(nGraphPoints);
    TGraphErrors *gMeanYield = new TGraphErrors(nGraphPoints);

    // Now get the pt distributions for the given bin ranges and calculate mean pT and yields
    for (int idx = 0; idx < nGraphPoints; ++idx)
    {
        int binStart = binEdges[idx];
        int binEnd = (idx + 1 < (int)binEdges.size()) ? binEdges[idx + 1] : hNFT0M->GetNbinsX();
        TH1D *hPt = hParticlevsFT0M->ProjectionX(Form("hPt_%d", idx + 1), binStart, binEnd);
        hPt->SetDirectory(nullptr);

        const double xVal = xValuesFT0Mcorrespondingdnch_detaRun3[idx];
        const double meanPt = (hPt->GetEntries() > 0) ? hPt->GetMean() : 0.0;
        const double rms = (hPt->GetEntries() > 1) ? hPt->GetRMS() : 0.0;
        gMeanpTvsNch->SetPoint(idx, xVal, meanPt);
        gMeanpTvsNch->SetPointError(idx, 0, (hPt->GetEntries() > 0) ? rms / sqrt(hPt->GetEntries()) : 0);

        double entries = hNFT0M->Integral(binStart, binEnd);
        double meanYield = (entries > 0) ? hPt->Integral() / entries : 0.0; // Normalize by number of events
        gMeanYield->SetPoint(idx, xVal, meanYield);
        gMeanYield->SetPointError(idx, 0, (entries > 0) ? sqrt(meanYield * (1 + meanYield) / entries) : 0);

        delete hPt;
    }

    TFile *fData = new TFile("levy_fit.root", "read");
    // TFile *fData = new TFile("PiKp_Run3_Results/Sawan/Pr_results.root", "read");
    if (fData->IsZombie())
    {
        cout << "Error: Data file not found" << endl;
        return;
    }

    TGraphErrors *gMeanpTRun3Data = (TGraphErrors *)fData->Get("gMeanpTRun3");
    TGraphErrors *gMeanYieldRun3Data = (TGraphErrors *)fData->Get("gMeanYieldRun3");

    // TGraphErrors *gMeanpTRun2Data = (TGraphErrors *)fData->Get("gMeanpTRun2");
    // TGraphErrors *gMeanYieldRun2Data = (TGraphErrors *)fData->Get("gMeanYieldRun2");

    if (gMeanpTRun3Data == nullptr)
    {
        cout << "Error: Data <p_T> graphs not found" << endl;
        return;
    }
    gStyle->SetCanvasPreferGL(kTRUE);

    TCanvas *cMeanpTvsNch = new TCanvas("cMeanpTvsNch", "cMeanpTvsNch", 720, 720);
    SetCanvasStyle(cMeanpTvsNch, 0.15, 0.03, 0.03, 0.15);
    SetGraphStyle(gMeanpTvsNch);
    gMeanpTvsNch->SetMarkerStyle(20);
    // gMeanpTvsNch->SetMarkerSize(1.2);
    // gMeanpTvsNch->SetFillStyle(3002);
    // gMeanpTvsNch->SetMarkerColor(kGreen + 2);
    gMeanpTvsNch->SetFillColorAlpha(kGreen + 2, 0.6);
    gMeanpTvsNch->SetLineColor(kGreen + 2);
    gMeanpTvsNch->SetLineWidth(2);
    // gMeanpTvsNch->Draw("AP");
    gMeanpTRun3Data->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/c)");
    gMeanpTRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanpTRun3Data->SetMarkerStyle(24);
    gMeanpTRun3Data->SetMarkerSize(1.2);
    gMeanpTRun3Data->SetMarkerColor(kRed);
    gMeanpTRun3Data->SetLineColor(kRed);
    gMeanpTRun3Data->GetYaxis()->SetRangeUser(0.31, 1.79);
    gMeanpTRun3Data->GetXaxis()->SetLimits(0, 28.9);
    gMeanpTRun3Data->Draw("AP");
    gMeanpTvsNch->Draw("E3 same");
    gMeanpTvsNch->Draw("l same");

    // gMeanpTRun2Data->SetMarkerStyle(25);
    // gMeanpTRun2Data->SetMarkerSize(1.2);
    // gMeanpTRun2Data->SetMarkerColor(kBlue);
    // gMeanpTRun2Data->SetLineColor(kBlue);
    // gMeanpTRun2Data->Draw("P same");

    TCanvas *cMeanYieldvsNch = new TCanvas("cMeanYieldvsNch", "cMeanYieldvsNch", 720, 720);
    SetCanvasStyle(cMeanYieldvsNch, 0.15, 0.03, 0.03, 0.15);
    SetGraphStyle(gMeanYield);
    gMeanYield->SetMarkerStyle(20);
    // gMeanYield->SetFillStyle(3002);
    // gMeanYield->SetMarkerSize(1.2);
    // gMeanYield->SetMarkerColor(kGreen + 2);
    gMeanYield->SetFillColorAlpha(kGreen + 2, 0.6);
    gMeanYield->SetLineColor(kGreen + 2);
    gMeanYield->SetLineWidth(2);
    // gMeanYield->Draw("AP");
    gMeanYieldRun3Data->GetYaxis()->SetTitle("dN/dy");
    gMeanYieldRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanYieldRun3Data->SetMarkerStyle(24);
    gMeanYieldRun3Data->SetMarkerSize(1.2);
    gMeanYieldRun3Data->SetMarkerColor(kRed);
    gMeanYieldRun3Data->SetLineColor(kRed);
    gMeanYieldRun3Data->GetYaxis()->SetRangeUser(0, 0.69);
    gMeanYieldRun3Data->GetXaxis()->SetLimits(0, 28.9);
    gMeanYieldRun3Data->Draw("AP");
    gMeanYield->Draw("E3 same");
    gMeanYield->Draw("l same");

    // gMeanYieldRun2Data->SetMarkerStyle(25);
    // gMeanYieldRun2Data->SetMarkerSize(1.2);
    // gMeanYieldRun2Data->SetMarkerColor(kBlue);
    // gMeanYieldRun2Data->SetLineColor(kBlue);
    // gMeanYieldRun2Data->Draw("P same");

    cout << "\nParticle type: " << particleTypes[selectParticle] << endl;
    cout << "Data file used: " << fData->GetName() << endl;
}