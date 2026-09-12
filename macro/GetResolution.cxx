#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/initializations.h"
#include <TSystem.h>
#include <TString.h>
#include <TStopwatch.h>

double voigtian(double *x, double *par)
{
    double amplitude = par[0];
    double mean = par[1];
    double sigma = par[2];
    double gamma = par[3];

    // Voigt profile using ROOT's TMath::Voigt
    return amplitude * TMath::Voigt(x[0] - mean, sigma, gamma);
}

void GetResolution()
{
    gStyle->SetOptFit(1111);
    bool isINEL = true;
    // TString MCpath = "../mc/LHC24f3c/750013.root"; // INEL>0, latest train with systematics
    TString MCpath = "../mc/LHC24f3c/755334.root"; // INEL, latest train with systematics
    TFile *fileMC = new TFile(MCpath, "READ");
    THnSparseF *hMCSignal = (THnSparseF *)fileMC->Get("kstarqa/hInvMass/h3KstarMassRec");
    // TString savePath = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ResolutionFromMC"; // INEL>0
    TString savePath = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/756343/kstarqa/hInvMass/ROTATED/ResolutionFromMC"; // INEL
    // Axes: pT, multiplicity, mass
    vector<float> Resolution;

    // float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    float mult_classes[] = {0.0};
    int totalMultClasses = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    cout << "Total pT bins " << Npt << ", Total multiplicity classes " << totalMultClasses << endl;

    // for (int imult = 0; imult < totalMultClasses + 1; imult++)
    for (int imult = 0; imult < 1; imult++)
    {
        int multlow, multhigh;
        if (imult == 0)
        {
            multlow = 0;
            multhigh = (isINEL) ? 120 : 100;
        }
        else
        {
            multlow = mult_classes[imult - 1];
            multhigh = mult_classes[imult];
        }

        int multLowBin = hMCSignal->GetAxis(1)->FindBin(multlow + 1e-3);
        int multHighBin = hMCSignal->GetAxis(1)->FindBin(multhigh - 1e-3);
        hMCSignal->GetAxis(1)->SetRange(multLowBin, multHighBin);

        for (int ipt = 0; ipt < Npt; ipt++)
        {
            float lowpt = pT_bins[ipt];
            float highpt = pT_bins[ipt + 1];

            cout << "Low pT: " << lowpt << ", High pT: " << highpt << ", Mult low: " << multlow << ", Mult high: " << multhigh << endl;

            int ptLowBin = hMCSignal->GetAxis(0)->FindBin(lowpt + 1e-3);
            int ptHighBin = hMCSignal->GetAxis(0)->FindBin(highpt - 1e-3);
            hMCSignal->GetAxis(0)->SetRange(ptLowBin, ptHighBin);

            TH1D *hMassProjection = hMCSignal->Projection(2, "E");
            hMassProjection->Rebin(2);

            float lowFitRange = 0.75;
            float highFitRange = 1.1;

            if (ipt == 0)
            {
                lowFitRange = 0.84;
                highFitRange = 0.95;
            }
            else if (ipt == 1)
            {
                lowFitRange = 0.83;
                highFitRange = 0.96;
            }
            else if (ipt == 2)
            {
                lowFitRange = 0.82;
                highFitRange = 0.96;
            }
            else if (ipt == 5)
            {
                lowFitRange = 0.79;
                highFitRange = 0.99;
            }
            else if (ipt == 18 || ipt == 19 || ipt == 22)
            {
                lowFitRange = 0.75;
                highFitRange = 1.05;
            }
            else
            {
                lowFitRange = 0.8;
                highFitRange = 1.00;
            }

            // Fit with voigtian
            TF1 *voigtFit = new TF1(Form("voigtFit_mult%d_pt%d", imult, ipt), voigtian, lowFitRange, highFitRange, 4);
            voigtFit->SetParameters(1, 0.8956, 0.05, 0.01);
            voigtFit->SetParNames("Amplitude", "Mean", "Sigma", "Gamma");
            voigtFit->FixParameter(3, 0.0471); // Fixing Gamma to known value from PDG
            hMassProjection->Fit(voigtFit, "R");
            hMassProjection->GetXaxis()->SetRangeUser(0.7, 1.1);

            TCanvas *cMassFit = new TCanvas(Form("cMassFit_mult%d_pt%d", imult, ipt), Form("Mass Fit for mult %d, pt %d", imult, ipt), 720, 720);
            SetHistoQA(hMassProjection);
            hMassProjection->Draw();
            voigtFit->Draw("SAME");
            cMassFit->SaveAs(Form("%s/mult_%d-%d/MassFit_pt_%.1f_%.1f.png", savePath.Data(), multlow, multhigh, lowpt, highpt));

            Resolution.push_back(voigtFit->GetParameter(2)); // Sigma is the resolution
        }
    }

    for (size_t i = 0; i < Resolution.size(); i++)
    {
        cout << "Resolution for bin " << i << ": " << Resolution[i] << endl;
    }
}