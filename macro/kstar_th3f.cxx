
#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "TArrow.h"
#include "src/style.h"
#include "src/fitfunc.h"

using namespace std;

void kstar_th3f()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    TFile *file = new TFile("/home/sawan/check_k892/data/kstar/LHC22o_pass7/Reflection/AnalysisResults_subhadeep_midbkg.root");
    if (file->IsZombie())
    {
        cerr << "File not found " << endl;
        return;
    }
    //// In the plot, xaxis: pT, yaxis: centrality, zaxis: invariant mass
    TString savepath = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/Reflection/";

    TString foldername = "kstar892-light-ion_id48735";
    TString misidentifiedResonances[] = {"hRhoToKpi", "hOmegaToKpi", "hPhiToKpi", "hKstarSelf", "hEtaToKpi", "hEtaPrimeToKpi"};
    TH3F *h3Rho = (TH3F *)file->Get(Form("%s/MC/Reflections/%s", foldername.Data(), misidentifiedResonances[0].Data()));
    TH3F *h3Omega = (TH3F *)file->Get(Form("%s/MC/Reflections/%s", foldername.Data(), misidentifiedResonances[1].Data()));
    TH3F *h3Phi = (TH3F *)file->Get(Form("%s/MC/Reflections/%s", foldername.Data(), misidentifiedResonances[2].Data()));
    TH3F *h3KstarSelf = (TH3F *)file->Get(Form("%s/MC/Reflections/%s", foldername.Data(), misidentifiedResonances[3].Data()));
    TH3F *h3Eta = (TH3F *)file->Get(Form("%s/MC/Reflections/%s", foldername.Data(), misidentifiedResonances[4].Data()));
    TH3F *h3EtaPrime = (TH3F *)file->Get(Form("%s/MC/Reflections/%s", foldername.Data(), misidentifiedResonances[5].Data()));
    if (h3Rho == nullptr || h3Omega == nullptr || h3Phi == nullptr || h3KstarSelf == nullptr || h3Eta == nullptr || h3EtaPrime == nullptr)
    {
        cerr << "One or more histograms not found!!!!!!!!!!!!" << endl;
        return;
    }
    double pT_bins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 10.0};
    // double pT_bins[] = {7.0, 10.0};
    int Npt = sizeof(pT_bins) / sizeof(pT_bins[0]) - 1;
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.04);
    vector<vector<float>> fitParameters;

    for (int ibin = 0; ibin < Npt; ibin++)
    {
        int lowpTbin = h3Rho->GetXaxis()->FindBin(pT_bins[ibin] + 1e-5);
        int highpTbin = h3Rho->GetXaxis()->FindBin(pT_bins[ibin + 1] - 1e-5);
        TH1D *h1Rho = h3Rho->ProjectionZ(Form("h1Rho_pT_%d-%d", int(pT_bins[ibin]), int(pT_bins[ibin + 1])), lowpTbin, highpTbin, -1, -1);
        TH1D *h1Omega = h3Omega->ProjectionZ(Form("h1Omega_pT_%d-%d", int(pT_bins[ibin]), int(pT_bins[ibin + 1])), lowpTbin, highpTbin, -1, -1);
        TH1D *h1Phi = h3Phi->ProjectionZ(Form("h1Phi_pT_%d-%d", int(pT_bins[ibin]), int(pT_bins[ibin + 1])), lowpTbin, highpTbin, -1, -1);
        TH1D *h1KstarSelf = h3KstarSelf->ProjectionZ(Form("h1KstarSelf_pT_%d-%d", int(pT_bins[ibin]), int(pT_bins[ibin + 1])), lowpTbin, highpTbin, -1, -1);
        TH1D *h1Eta = h3Eta->ProjectionZ(Form("h1Eta_pT_%d-%d", int(pT_bins[ibin]), int(pT_bins[ibin + 1])), lowpTbin, highpTbin, -1, -1);
        TH1D *h1EtaPrime = h3EtaPrime->ProjectionZ(Form("h1EtaPrime_pT_%d-%d", int(pT_bins[ibin]), int(pT_bins[ibin + 1])), lowpTbin, highpTbin, -1, -1);

        // TH1D *hSumAll = (TH1D *)h1Rho->Clone(Form("hSumAll_pT_%d-%d", int(pT_bins[ibin]), int(pT_bins[ibin + 1])));
        // hSumAll->Add(h1Omega);
        TH1D *hSumAll = (TH1D *)h1Omega->Clone(Form("hSumAll_pT_%d-%d", int(pT_bins[ibin]), int(pT_bins[ibin + 1])));
        hSumAll->Add(h1Phi);
        hSumAll->Add(h1KstarSelf);
        hSumAll->Add(h1Eta);
        hSumAll->Add(h1EtaPrime);

        TCanvas *cAddAll = new TCanvas("", "", 720, 720);
        SetCanvasStyle(cAddAll, 0.15, 0.06, 0.06, 0.13);
        SetHistoQA(hSumAll);
        hSumAll->SetMarkerStyle(20);
        // hSumAll->Rebin(4);
        hSumAll->Draw("pe");
        lat.DrawLatex(0.55, 0.96, Form("#it{p}_{T} = %.1f - %.1f GeV/c", pT_bins[ibin], pT_bins[ibin + 1]));

        float fitRangeLow = 0.64;
        float fitRangeHigh = 0.98;

        if (ibin == 0 || ibin == 1)
        {
            fitRangeLow = 0.65;
            fitRangeHigh = 0.92;
        }

        TF1 *CBRight = new TF1("CBRight", CrystalBallRight, fitRangeLow, fitRangeHigh, 5);
        CBRight->SetParameters(200000, 0.78, 0.07, 0.9, 0.5); // norm, mean, sigma, alpha, n
        CBRight->SetParLimits(1, 0.7, 0.9);
        CBRight->SetParLimits(2, 0.05, 0.09);
        CBRight->SetParLimits(3, 0.5, 1.2);
        hSumAll->Fit(CBRight, "REBMS");
        fitParameters.push_back({(float)CBRight->GetParameter(1), (float)CBRight->GetParameter(2), (float)CBRight->GetParameter(3), (float)CBRight->GetParameter(4)});
        cAddAll->SaveAs(Form("%sSumAll_pT_%d.png", savepath.Data(), ibin));
    }

    // Print results in a neat table
    cout << left << setw(20) << "pT range (GeV/c)" << setw(12) << "mean" << setw(12) << "sigma" << setw(12) << "alpha" << setw(12) << "n" << endl;
    cout << string(68, '-') << endl;
    for (int ibin = 0; ibin < Npt; ibin++)
    {
        std::ostringstream rng;
        rng << fixed << setprecision(2) << pT_bins[ibin] << "-" << pT_bins[ibin + 1];
        cout << left << setw(20) << rng.str()
             << right << setw(12) << setprecision(4) << fitParameters[ibin][0]
             << setw(12) << fitParameters[ibin][1]
             << setw(12) << fitParameters[ibin][2]
             << setw(12) << fitParameters[ibin][3] << endl;
    }
}

// Output table format:
// pT range (GeV/c)    mean        sigma       alpha       n
// --------------------------------------------------------------------
// 0.00-0.40                 0.7875     0.06119         1.2       146.8
// 0.40-0.80                 0.7527     0.07456      0.8809       139.5
// 0.80-1.20                 0.7678     0.08454      0.5499       127.8
// 1.20-1.60                 0.7684     0.08353      0.5506        1.23
// 1.60-2.00                 0.7686     0.08075      0.6135      0.6261
// 2.00-2.50                  0.769     0.07841      0.7388      0.3242
// 2.50-3.00                 0.7698     0.07658      0.8244      0.2302
// 3.00-3.50                 0.7709     0.07632      0.8935      0.1739
// 3.50-4.00                 0.7706     0.07539      0.9161      0.1594
// 4.00-4.50                 0.7715     0.07474      0.9677      0.1275
// 4.50-5.00                 0.7712     0.07387      0.9634      0.1449
// 5.00-6.00                 0.7712     0.07425      0.9524      0.1551
// 6.00-7.00                 0.7715     0.07425      0.9648      0.1463
// 7.00-10.00                0.7721     0.07309       1.031      0.1128
