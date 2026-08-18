#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

void toyModelMassShiftPhiData()
{
    // =========================================================
    // 1. Experimental Input Data (12 pT Bins)
    // =========================================================
    const std::vector<double> ptMin = {0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    const std::vector<double> ptMax = {0.6, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0};
    const std::vector<double> m2025 = {1.01868, 1.01875, 1.01893, 1.01900, 1.01898, 1.01907, 1.01915, 1.01921, 1.01927, 1.01931, 1.01936, 1.01931, 1.01929, 1.01925, 1.01912, 1.01915};
    const std::vector<double> m2026 = {1.01822, 1.01897, 1.01936, 1.01945, 1.01948, 1.01961, 1.01970, 1.01975, 1.01981, 1.01986, 1.01990, 1.01992, 1.01991, 1.01995, 1.01994, 1.01975};

    const int nBins = ptMin.size();
    const double mK = 0.493677; // Kaon PDG Mass (GeV/c^2)
    const int nEventsPerBin = 10000;

    // Scan range for relative momentum scale delta = dp / p
    const double deltaMin = -0.1; // -10%
    const double deltaMax = +0.1; // +10%
    const int nDelta = 5000;

    TRandom3 rand(12345);

    std::vector<double> ptCenters(nBins);
    std::vector<double> bestDeltas(nBins);
    std::vector<double> bestDeltasPercent(nBins);

    std::cout << "===========================================================================\n";
    std::cout << std::setw(12) << "pT Range"
              << std::setw(12) << "M_2025"
              << std::setw(12) << "M_2026"
              << std::setw(15) << "Shift (GeV)"
              << std::setw(18) << "Momentum Shift (%"
                                  ")\n";
    std::cout << "===========================================================================\n";

    // =========================================================
    // 2. Loop Over pT Bins to Extract Single-Track Delta
    // =========================================================
    for (int bin = 0; bin < nBins; ++bin)
    {
        ptCenters[bin] = (ptMin[bin] + ptMax[bin]) / 2.0;
        double targetMass2026 = m2026[bin];
        double baseMass2025 = m2025[bin];

        std::vector<TLorentzVector> kaon1, kaon2;
        kaon1.reserve(nEventsPerBin);
        kaon2.reserve(nEventsPerBin);

        // Generate synthetic decays with 2025 baseline kinematics in current pT bin
        for (int i = 0; i < nEventsPerBin; ++i)
        {
            double pt = rand.Uniform(ptMin[bin], ptMax[bin]);
            double phiAngle = rand.Uniform(-TMath::Pi(), TMath::Pi());
            double eta = rand.Uniform(-0.8, 0.8);

            double px = pt * std::cos(phiAngle);
            double py = pt * std::sin(phiAngle);
            double pz = pt * std::sinh(eta);
            double E = std::sqrt(px * px + py * py + pz * pz + baseMass2025 * baseMass2025);

            TLorentzVector phi(px, py, pz, E);
            TGenPhaseSpace decay;
            double masses[2] = {mK, mK};

            if (decay.SetDecay(phi, 2, masses))
            {
                decay.Generate();
                kaon1.push_back(*decay.GetDecay(0));
                kaon2.push_back(*decay.GetDecay(1));
            }
        }

        // Scan delta to match 2026 peak position
        double bestDelta = 0.0;
        double minDiff = 1e9;

        for (int i = 0; i < nDelta; ++i)
        {
            double delta = deltaMin + i * (deltaMax - deltaMin) / (nDelta - 1);
            TH1D hMass("hMass", "", 200, 1.01, 1.03);

            for (size_t j = 0; j < kaon1.size(); ++j)
            {
                double p1Corr = kaon1[j].P() * (1.0 + delta);
                double p2Corr = kaon2[j].P() * (1.0 + delta);

                TLorentzVector k1Corr, k2Corr;
                k1Corr.SetPtEtaPhiM(p1Corr / std::cosh(kaon1[j].Eta()), kaon1[j].Eta(), kaon1[j].Phi(), mK);
                k2Corr.SetPtEtaPhiM(p2Corr / std::cosh(kaon2[j].Eta()), kaon2[j].Eta(), kaon2[j].Phi(), mK);

                hMass.Fill((k1Corr + k2Corr).M());
            }

            // double currentPeak = hMass.GetXaxis()->GetBinCenter(hMass.GetMaximumBin());
            double currentPeak = hMass.GetMean();
            double diff = std::abs(currentPeak - targetMass2026);

            if (diff < minDiff)
            {
                minDiff = diff;
                bestDelta = delta;
            }
        }

        bestDeltas[bin] = bestDelta;
        bestDeltasPercent[bin] = bestDelta * 100.0;

        std::cout << Form("[%4.1f - %4.1f]", ptMin[bin], ptMax[bin]) << std::setw(6) << " "
                  << std::setw(8) << std::fixed << std::setprecision(5) << baseMass2025
                  << std::setw(12) << targetMass2026
                  << std::setw(14) << (targetMass2026 - baseMass2025)
                  << std::setw(16) << std::setprecision(4) << bestDeltasPercent[bin] << " %\n";
    }
    std::cout << "===========================================================================\n";

    // =========================================================
    // 3. Plot single-track momentum correction values as a function of pT
    // =========================================================
    TCanvas *c1 = new TCanvas("c1", "Momentum Shift Calibration", 800, 600);
    c1->SetGrid();

    TGraph *gCorrection = new TGraph(nBins, ptCenters.data(), bestDeltasPercent.data());
    gCorrection->SetTitle("Single-Track Momentum Correction Factor (#delta);#phi p_{T} (GeV/c);#delta = #Delta p / p (%)");
    gCorrection->SetMarkerStyle(20);
    gCorrection->SetMarkerSize(1.2);
    gCorrection->SetMarkerColor(kBlue + 2);
    gCorrection->SetLineColor(kBlue + 2);
    gCorrection->SetLineWidth(2);
    gCorrection->Draw("APL");
}
