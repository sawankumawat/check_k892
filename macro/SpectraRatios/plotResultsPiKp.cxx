#include <iostream>
#include <iomanip>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/initializations.h"

void plotResultsPiKp()
{
    vector<string> particleType = {"Pi", "Ka", "Pr"};
    int selectParticle = 0;
    TFile *fresultsPos = new TFile(Form("PiKp_Run3_Results/Yeild_meanpT/Yields_Pos_%s_AllMult.root", particleType[selectParticle].c_str()), "read");
    TFile *fresultsNeg = new TFile(Form("PiKp_Run3_Results/Yeild_meanpT/Yields_Neg_%s_AllMult.root", particleType[selectParticle].c_str()), "read");
    if (fresultsPos->IsZombie() || fresultsNeg->IsZombie())
    {
        cout << "Error: Results files not found" << endl;
        return;
    }
    double dnch_detaRun3[] = {21.78, 18.48, 15.76, 13.89, 12.50, 10.86, 9.09, 7.63, 5.87, 3.69};
    double dnch_detaRun3_err[] = {0.38, 0.25, 0.22, 0.19, 0.17, 0.15, 0.13, 0.11, 0.09, 0.06};
    float multClasses[] = {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    int nBinsMult = sizeof(dnch_detaRun3) / sizeof(dnch_detaRun3[0]);
    TGraphErrors *gMeanpTRun3Data = new TGraphErrors(nBinsMult);
    TGraphErrors *gMeanYieldRun3Data = new TGraphErrors(nBinsMult);

    for (int ibins = 1; ibins <= nBinsMult; ++ibins)
    {
        TH1F *hResultsPos = (TH1F *)fresultsPos->Get(Form("IntegratedWithLevi_%.2f_%.2f", multClasses[ibins - 1], multClasses[ibins]));
        TH1F *hResultsNeg = (TH1F *)fresultsNeg->Get(Form("IntegratedWithLevi_%.2f_%.2f", multClasses[ibins - 1], multClasses[ibins]));
        double yield = (hResultsPos->GetBinContent(1) + hResultsNeg->GetBinContent(1)) / 2.0;
        double yieldErr = sqrt(pow(hResultsPos->GetBinContent(3), 2) + pow(hResultsNeg->GetBinContent(3), 2)) / 2.0;
        double meanpT = (hResultsPos->GetBinContent(5) + hResultsNeg->GetBinContent(5)) / 2.0;
        double meanpTErr = sqrt(pow(hResultsPos->GetBinContent(7), 2) + pow(hResultsNeg->GetBinContent(7), 2)) / 2.0;
        gMeanYieldRun3Data->SetPoint(ibins - 1, dnch_detaRun3[ibins - 1], yield);
        gMeanYieldRun3Data->SetPointError(ibins - 1, dnch_detaRun3_err[ibins - 1], yieldErr);
        gMeanpTRun3Data->SetPoint(ibins - 1, dnch_detaRun3[ibins - 1], meanpT);
        gMeanpTRun3Data->SetPointError(ibins - 1, dnch_detaRun3_err[ibins - 1], meanpTErr);
    }

    TCanvas *cMeanpTvsNch = new TCanvas("cMeanpTvsNch", "cMeanpTvsNch", 720, 720);
    gMeanpTRun3Data->SetMarkerStyle(20);
    gMeanpTRun3Data->SetTitle(0);
    gMeanpTRun3Data->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/c)");
    gMeanpTRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanpTRun3Data->Draw("AP");

    TCanvas *cMeanYieldvsNch = new TCanvas("cMeanYieldvsNch", "cMeanYieldvsNch", 720, 720);
    gMeanYieldRun3Data->SetMarkerStyle(20);
    gMeanYieldRun3Data->SetTitle(0);
    gMeanYieldRun3Data->GetYaxis()->SetTitle("dN/dy");
    gMeanYieldRun3Data->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    gMeanYieldRun3Data->Draw("AP");

    TFile *fOutput = new TFile(Form("PiKp_Run3_Results/Sawan/%s_results.root", particleType[selectParticle].c_str()), "recreate");
    gMeanpTRun3Data->Write("gMeanpTRun3");
    gMeanYieldRun3Data->Write("gMeanYieldRun3");

    TH1D *hSpectraPos[nBinsMult + 1];
    TH1D *hSpectraNeg[nBinsMult + 1];
    TH1D *hSpectraSum[nBinsMult + 1];
    for (int imult = 0; imult < nBinsMult + 1; ++imult)
    {
        int multlow, multhigh;
        if (imult == 0)
        {
            multlow = 0;
            multhigh = 100;
        }
        else
        {
            multlow = (int)multClasses[imult - 1];
            multhigh = (int)multClasses[imult];
        }
        TFile *fSpectraPos = new TFile(Form("PiKp_Run3_Results/Pos%s_LHC24_pass1_INELgt0_cent%d_%d.root", particleType[selectParticle].c_str(), multlow, multhigh), "read");
        TFile *fSpectraNeg = new TFile(Form("PiKp_Run3_Results/Neg%s_LHC24_pass1_INELgt0_cent%d_%d.root", particleType[selectParticle].c_str(), multlow, multhigh), "read");
        if (fSpectraPos->IsZombie() || fSpectraNeg->IsZombie())
        {
            cout << "Error: Spectra files not found for multiplicity class " << multlow << "-" << multhigh << endl;
            return;
        }

        hSpectraPos[imult] = (TH1D *)fSpectraPos->Get("cor_sys");
        hSpectraNeg[imult] = (TH1D *)fSpectraNeg->Get("cor_sys");
        hSpectraSum[imult] = (TH1D *)hSpectraPos[imult]->Clone(Form("hSpectraSum%d", imult));
        hSpectraSum[imult]->SetDirectory(nullptr);
        hSpectraSum[imult]->Add(hSpectraNeg[imult]);
        fOutput->cd();
        hSpectraSum[imult]->Write(Form("hCorrectedSpectra_%d_%d", multlow, multhigh));
    }

    fOutput->Close();
}