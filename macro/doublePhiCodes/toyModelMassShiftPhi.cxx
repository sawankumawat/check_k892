#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <iostream>
#include <cmath>

void toyModelMassShiftPhi()
{
    // Measured peak positions in Data
    const double M2025 = 1.01893;              // GeV/c2
    const double M2026 = 1.01974;              // GeV/c2
    const double targetShift = M2026 - M2025; // Target shift to reproduce

    // PDG masses
    const double mPhiPDG = 1.019461;
    const double mK = 0.493677;

    const double ptMin = 0.8;
    const double ptMax = 1.2;
    const int nEvents = 10000;

    const double deltaMin = -0.1; // -10%
    const double deltaMax = +0.1; // +10%
    const int nDelta = 5000;

    TRandom3 rand(12345);
    std::vector<TLorentzVector> kaon1, kaon2;
    kaon1.reserve(nEvents);
    kaon2.reserve(nEvents);

    // 1. Generate decays at TRUE PDG mass
    for (int i = 0; i < nEvents; ++i)
    {
        double pt = rand.Exp(1.0) + ptMin;
        if (pt > ptMax)
            continue;

        double phiAngle = rand.Uniform(-TMath::Pi(), TMath::Pi());
        double eta = rand.Uniform(-0.8, 0.8);

        double px = pt * std::cos(phiAngle);
        double py = pt * std::sin(phiAngle);
        double pz = pt * std::sinh(eta);
        double E = std::sqrt(px * px + py * py + pz * pz + M2025 * M2025);

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

    // 2. Scan momentum scale factor (delta)
    double bestDelta = 0;
    double minDiff = 1e9;

    for (int i = 0; i < nDelta; ++i)
    {
        double delta = deltaMin + i * (deltaMax - deltaMin) / (nDelta - 1);
        TH1D hMass("hMass", "", 200, 1.00, 1.04);

        for (size_t j = 0; j < kaon1.size(); ++j)
        {
            double p1Corr = kaon1[j].P() * (1.0 + delta);
            double p2Corr = kaon2[j].P() * (1.0 + delta);

            TLorentzVector k1Corr, k2Corr;
            k1Corr.SetPtEtaPhiM(p1Corr / std::cosh(kaon1[j].Eta()), kaon1[j].Eta(), kaon1[j].Phi(), mK);
            k2Corr.SetPtEtaPhiM(p2Corr / std::cosh(kaon2[j].Eta()), kaon2[j].Eta(), kaon2[j].Phi(), mK);

            hMass.Fill((k1Corr + k2Corr).M());
        }

        // double currentShift = hMass.GetXaxis()->GetBinCenter(hMass.GetMaximumBin()) - M2025;
        double currentShift = hMass.GetMean() - M2025;
        double diff = std::abs(currentShift - targetShift);

        if (diff < minDiff)
        {
            minDiff = diff;
            bestDelta = delta;
        }
    }

    std::cout << "=====================================\n";
    std::cout << "Optimal momentum scale shift (delta) = " << bestDelta * 100.0 << " %\n";
    std::cout << "=====================================\n";

    // 3. Fill Histograms Before and After Correction
    TH1D *hBefore = new TH1D("hBefore", "#phi Mass Peak Comparison;M_{K^{+}K^{-}} (GeV/c^{2});Counts", 200, 1.00, 1.04);
    TH1D *hAfter = new TH1D("hAfter", "#phi Mass Peak Comparison;M_{K^{+}K^{-}} (GeV/c^{2});Counts", 200, 1.00, 1.04);

    for (size_t j = 0; j < kaon1.size(); ++j)
    {
        // Original spectrum (Before correction)
        TLorentzVector k1 = kaon1[j];
        TLorentzVector k2 = kaon2[j];
        hBefore->Fill((k1 + k2).M());

        // Recalculate daughter momenta using optimal bestDelta (After correction)
        double p1Corr = k1.P() * (1.0 + bestDelta);
        double p2Corr = k2.P() * (1.0 + bestDelta);

        TLorentzVector k1Corr, k2Corr;
        k1Corr.SetPtEtaPhiM(p1Corr / std::cosh(k1.Eta()), k1.Eta(), k1.Phi(), mK);
        k2Corr.SetPtEtaPhiM(p2Corr / std::cosh(k2.Eta()), k2.Eta(), k2.Phi(), mK);

        hAfter->Fill((k1Corr + k2Corr).M());
    }

    // 4. Plot both distributions on the same Canvas
    TCanvas *c1 = new TCanvas("c1", "Phi Mass Shift Comparison", 800, 600);
    c1->SetGrid();

    // Line styles and colors
    hBefore->SetLineColor(kBlue + 1);
    hBefore->SetLineWidth(2);
    hBefore->SetStats(0);

    hAfter->SetLineColor(kRed + 1);
    hAfter->SetLineWidth(2);
    hAfter->SetLineStyle(2); // Dashed line for clarity
    hAfter->SetStats(0);

    // Dynamic vertical scale adjustment
    double maxY = std::max(hBefore->GetMaximum(), hAfter->GetMaximum());
    hBefore->SetMaximum(maxY * 1.2);

    hBefore->Draw("HIST");
    hAfter->Draw("HIST SAME");

    // Add Legend
    TLegend *legend = new TLegend(0.60, 0.72, 0.88, 0.88);
    legend->SetBorderSize(1);
    legend->AddEntry(hBefore, Form("Before Correction: Mean = %.4f", hBefore->GetMean()), "l");
    legend->AddEntry(hAfter, Form("After Correction: Mean = %.4f", hAfter->GetMean()), "l");
    legend->AddEntry((TObject *)0, Form("Optimal #Delta = %.4f %%", bestDelta * 100.0), "");
    legend->Draw();
}