#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "src/style.h"

void toyModelMassShiftPhiDataMatch()
{
    // =========================================================
    // 1. Experimental Input Data
    // =========================================================

    const std::vector<double> ptMin = {
        0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0,
        4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0};

    const std::vector<double> ptMax = {
        0.8, 1.2, 1.6, 2.0, 2.5, 3.0,
        4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};

    // 2026 Mass Values
    std::vector<double> m2026 = {
        1.01892, 1.01936, 1.01945, 1.01948, 1.01961,
        1.0197, 1.01975, 1.01981, 1.01986, 1.0199,
        1.01991, 1.01992, 1.01994, 1.01998};

    // 2025 Mass Values
    std::vector<double> m2025 = {
        1.01879, 1.01897, 1.01905, 1.01902, 1.01912,
        1.01919, 1.01925, 1.0193, 1.01933, 1.01939,
        1.01938, 1.01932, 1.01929, 1.01923};

    const int nBins = ptMin.size();

    const double mK = 0.493677;
    const int nEventsPerBin = 5000;

    // =========================================================
    // Scan range for epsilon = dp/p
    // =========================================================

    const double epsilonMin = -0.1;
    const double epsilonMax = +0.1;
    const int nEpsilon = 1000;

    TRandom3 rand(12345);

    std::vector<double> ptCenters(nBins);
    std::vector<double> bestEpsilon(nBins);
    std::vector<double> bestEpsilonPercent(nBins);

    // Results of closure toy
    std::vector<double> correctedMass(nBins);
    std::vector<double> closureDifference(nBins);

    std::cout
        << "===============================================================================================\n";

    std::cout
        << std::setw(12) << "pT Range"
        << std::setw(12) << "M_2025"
        << std::setw(12) << "M_2026"
        << std::setw(14) << "Shift (MeV)"
        << std::setw(18) << "epsilon (%)"
        << std::setw(18) << "M_2025 corr."
        << std::setw(18) << "Closure diff."
        << "\n";

    std::cout
        << "===============================================================================================\n";

    // =========================================================
    // 2. Loop over pT bins
    // =========================================================

    for (int bin = 0; bin < nBins; ++bin)
    {
        ptCenters[bin] = 0.5 * (ptMin[bin] + ptMax[bin]);

        const double targetMass2026 = m2026[bin];
        const double baseMass2025 = m2025[bin];

        std::vector<TLorentzVector> kaon1;
        std::vector<TLorentzVector> kaon2;

        kaon1.reserve(nEventsPerBin);
        kaon2.reserve(nEventsPerBin);

        // =====================================================
        // Generate synthetic 2025-like phi -> K+ K- decays
        // =====================================================

        for (int i = 0; i < nEventsPerBin; ++i)
        {
            const double pt = rand.Uniform(ptMin[bin], ptMax[bin]);
            const double phiAngle = rand.Uniform(-TMath::Pi(), TMath::Pi());
            const double eta = rand.Uniform(-0.8, 0.8);
            const double px = pt * std::cos(phiAngle);
            const double py = pt * std::sin(phiAngle);
            const double pz = pt * std::sinh(eta);
            const double E = std::sqrt(px * px + py * py + pz * pz + baseMass2025 * baseMass2025);

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

        // =====================================================
        // 3. Scan epsilon
        //
        // pK(corrected) = (1 + epsilon) pK(2025)
        // =====================================================

        double bestEps = 0.0;
        double minDiff = 1e9;

        for (int i = 0; i < nEpsilon; ++i)
        {
            const double epsilon = epsilonMin + i * (epsilonMax - epsilonMin) / (nEpsilon - 1);
            TH1D hMass(Form("hMass_bin%d_eps%d", bin, i), "", 200, 1.01, 1.03);

            for (size_t j = 0; j < kaon1.size(); ++j)
            {
                // =================================================
                // Apply epsilon to THREE-MOMENTUM of each kaon
                // =================================================

                const double p1Corr = kaon1[j].P() * (1.0 + epsilon);
                const double p2Corr = kaon2[j].P() * (1.0 + epsilon);

                // =================================================
                // Recalculate energy using kaon mass
                // =================================================

                TLorentzVector k1Corr;
                TLorentzVector k2Corr;

                k1Corr.SetPtEtaPhiM(p1Corr / std::cosh(kaon1[j].Eta()), kaon1[j].Eta(), kaon1[j].Phi(), mK);
                k2Corr.SetPtEtaPhiM(p2Corr / std::cosh(kaon2[j].Eta()), kaon2[j].Eta(), kaon2[j].Phi(), mK);

                // =================================================
                // Reconstruct phi invariant mass
                // =================================================

                const double mass = (k1Corr + k2Corr).M();

                hMass.Fill(mass);
            }

            // -----------------------------------------------------
            // Use histogram mean to determine the response
            // -----------------------------------------------------

            const double currentMass = hMass.GetMean();
            const double diff = std::abs(currentMass - targetMass2026);

            if (diff < minDiff)
            {
                minDiff = diff;
                bestEps = epsilon;
            }
        }

        bestEpsilon[bin] = bestEps;
        bestEpsilonPercent[bin] = 100.0 * bestEps;

        // =========================================================
        // 4. CLOSURE TOY
        //
        // Apply the obtained epsilon to the 2025 kaons again
        // and check whether the peak reaches the 2026 value.
        // =========================================================

        TH1D hCorrected(Form("hCorrected_bin%d", bin), Form("%.1f < p_{T,#phi} < %.1f GeV/c", ptMin[bin], ptMax[bin]), 200, 1.01, 1.03);

        TH1D hNominal(Form("hNominal_bin%d", bin), "", 200, 1.01, 1.03);

        for (size_t j = 0; j < kaon1.size(); ++j)
        {
            // =====================================================
            // Nominal 2025 kaons
            // =====================================================

            hNominal.Fill((kaon1[j] + kaon2[j]).M());

            // =====================================================
            // Corrected 2025 kaons
            // =====================================================

            const double p1Corr = kaon1[j].P() * (1.0 + bestEps);
            const double p2Corr = kaon2[j].P() * (1.0 + bestEps);

            TLorentzVector k1Corr;
            TLorentzVector k2Corr;

            k1Corr.SetPtEtaPhiM(p1Corr / std::cosh(kaon1[j].Eta()), kaon1[j].Eta(), kaon1[j].Phi(), mK);

            k2Corr.SetPtEtaPhiM(p2Corr / std::cosh(kaon2[j].Eta()), kaon2[j].Eta(), kaon2[j].Phi(), mK);
            hCorrected.Fill((k1Corr + k2Corr).M());
        }

        // =========================================================
        // Determine corrected peak position
        // =========================================================

        const double toyMass2025 = hNominal.GetMean();
        const double toyMassCorrected = hCorrected.GetMean();
        correctedMass[bin] = toyMassCorrected;
        closureDifference[bin] = toyMassCorrected - targetMass2026;

        // =========================================================
        // Print result
        // =========================================================

        std::cout
            << Form("[%4.1f - %4.1f]",
                    ptMin[bin], ptMax[bin])
            << std::setw(6) << " "
            << std::setw(8)
            << std::fixed
            << std::setprecision(5)
            << baseMass2025
            << std::setw(12)
            << targetMass2026
            << std::setw(14)
            << std::setprecision(3)
            << 1000.0 *
                   (targetMass2026 - baseMass2025)
            << std::setw(15)
            << std::setprecision(4)
            << bestEpsilonPercent[bin]
            << " %"
            << std::setw(15)
            << std::setprecision(5)
            << toyMassCorrected
            << std::setw(15)
            << std::setprecision(3)
            << 1000.0 *
                   closureDifference[bin]
            << " MeV"
            << "\n";
    }

    std::cout
        << "===============================================================================================\n";

    // =========================================================
    // 5. Plot epsilon correction
    // =========================================================

    TCanvas *c1 = new TCanvas("c1", "Momentum Scale Correction", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.03, 0.05, 0.15);
    c1->SetGrid();

    TGraph *gCorrection = new TGraph(nBins, ptCenters.data(), bestEpsilonPercent.data());

    gCorrection->SetTitle("Kaon Momentum Scale Correction;"
                          "#phi p_{T} (GeV/c);"
                          "#epsilon = #Delta p / p (%)");
    SetGraphStyle(gCorrection);
    gCorrection->SetMarkerStyle(20);
    gCorrection->SetMarkerSize(1.2);
    gCorrection->SetMarkerColor(kBlue + 2);
    gCorrection->SetLineColor(kBlue + 2);
    gCorrection->SetLineWidth(2);

    gCorrection->Draw("APL");

    // =========================================================
    // 6. Plot 2025 vs 2026 vs corrected 2025
    // =========================================================

    TCanvas *c2 = new TCanvas("c2", "Phi Mass Closure", 720, 720);
    SetCanvasStyle(c2, 0.17, 0.03, 0.09, 0.13);
    c2->SetGrid();

    TGraph *g2025 = new TGraph(nBins, ptCenters.data(), m2025.data());
    TGraph *g2026 = new TGraph(nBins, ptCenters.data(), m2026.data());
    TGraph *gCorrected = new TGraph(nBins, ptCenters.data(), correctedMass.data());

    g2025->SetTitle(
        "#phi mass closure;"
        "#phi p_{T} (GeV/c);"
        "M_{#phi} (GeV/c^{2})");

    SetGraphStyle(g2025);
    SetGraphStyle(g2026);
    SetGraphStyle(gCorrected);

    g2025->SetMarkerStyle(20);
    g2025->SetMarkerColor(kRed + 1);
    g2025->SetLineColor(kRed + 1);

    g2026->SetMarkerStyle(21);
    g2026->SetMarkerColor(kBlue + 1);
    g2026->SetLineColor(kBlue + 1);

    gCorrected->SetMarkerStyle(24);
    gCorrected->SetMarkerColor(kGreen + 2);
    gCorrected->SetLineColor(kGreen + 2);

    g2025->SetMinimum(1.0176);
    g2025->SetMaximum(1.0206);
    g2025->GetYaxis()->SetTitleOffset(1.9);
    g2025->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    g2025->Draw("APL");
    g2026->Draw("PL SAME");
    gCorrected->Draw("PL SAME");

    TLegend *legend = new TLegend(0.4, 0.20, 0.7, 0.4);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);

    legend->AddEntry(g2025, "2025", "lp");
    legend->AddEntry(g2026, "2026", "lp");
    legend->AddEntry(gCorrected, "2025 + #epsilon correction", "lp");
    legend->Draw();
    c2->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/PhiInvMass/PhiMassClosure.png");
}
