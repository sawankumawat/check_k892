#include <TGenPhaseSpace.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "src/style.h"
#include <Math/VectorUtil.h>
using namespace std;

void toy_model_ptSmear()
{
    TRandom3 randGen(0); // Random seed
    TFile *feffK0s = new TFile("/home/sawan/check_k892/injected_mc/Ks_mc.root", "read");
    if (feffK0s->IsZombie())
    {
        cout << "Error opening efficiency file" << endl;
        return;
    }
    TH2D *hgen = (TH2D *)feffK0s->Get("lf-v0qaanalysis/Generated_MCGenRecoColl_INELgt0_K0Short");
    TH2D *hrec = (TH2D *)feffK0s->Get("lf-v0postprocessing/hMassVsPtK0Short");
    if (hgen == nullptr || hrec == nullptr)
    {
        cout << "Histograms not found" << endl;
        return;
    }
    TH1D *h1gen = hgen->ProjectionX();
    TH1D *h1rec = hrec->ProjectionX();
    float finepTbins[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                          1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                          2, 2.2, 2.4, 2.6, 2.8,
                          3, 3.2, 3.4, 3.6, 3.8,
                          4, 4.2, 4.4, 4.6, 4.8,
                          5, 5.5,
                          6, 7,
                          8, 9,
                          10, 12,
                          15};
    // TH1F *hOriginalK0sEff = (TH1F *)h1rec->Clone();
    // hOriginalK0sEff->Divide(h1gen);
    // TFile *fK0sEffTemp = new TFile("K0sEff.root", "recreate");
    // h1gen->Write("K0s_generated_pT");
    // h1rec->Write("K0s_reconstructed_pT");
    // hOriginalK0sEff->Write("K0s_efficiency");
    // fK0sEffTemp->Close();

    int nFineBins = sizeof(finepTbins) / sizeof(finepTbins[0]) - 1;
    TH1F *hEffK0s = new TH1F("hEffK0s", "K0s Acceptance #times Efficiency vs p_{T}; #it{p}_{T} (GeV/#it{c}); Acceptance #times Efficiency", nFineBins, finepTbins);
    for (int i = 1; i <= nFineBins; i++)
    {
        double genCount = h1gen->Integral(h1gen->GetXaxis()->FindBin(finepTbins[i - 1] + 0.0001), h1gen->GetXaxis()->FindBin(finepTbins[i] - 0.0001));
        double recCount = h1rec->Integral(h1rec->GetXaxis()->FindBin(finepTbins[i - 1] + 0.0001), h1rec->GetXaxis()->FindBin(finepTbins[i] - 0.0001));
        double efficiency = (genCount > 0) ? (recCount / genCount) : 0;
        hEffK0s->SetBinContent(i, efficiency);
    }

    gStyle->SetOptStat(0);
    // Define masses (GeV/c^2)
    double m_mother = 1.710;       // Mass of f1525
    double m_daughter1 = 0.497611; // daughter 1 (K0s mass)
    double m_daughter2 = 0.497611; // daughter 2 (K0s mass)

    float ptBins[] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0};
    int nPtBins = sizeof(ptBins) / sizeof(ptBins[0]) - 1;

    // Histogram for decay product pT
    TH1F *h_pT = new TH1F("h_pT", "p_{T} distribution of mother; p_{T} (GeV/c); Events", nPtBins, ptBins);
    TH1F *hrec_pT = new TH1F("hrec_pT", "Reconstructed p_{T} distribution of mother; p_{T} (GeV/c); Events", nPtBins, ptBins);

    // Loop over multiple decay events
    int nEvents = 1e6;
    TLorentzVector lvmother;
    ROOT::Math::PxPyPzMVector fourVecMother, fourVecDau1, fourVecDau2;
    ROOT::Math::XYZVector threeVecDauCM, threeVecMother;
    vector<double> EffpT1; // 1-2 GeV/c
    vector<double> EffpT2; // 2-3 GeV/c
    vector<double> EffpT3; // 3-5 GeV/c
    vector<double> EffpT4; // 5-7 GeV/c
    vector<double> EffpT5; // 7-10 GeV/c
    vector<double> EffpT6; // 10-15 GeV/c

    for (int i = 0; i < nEvents; i++)
    {
        // Generate in pt, y and phi
        double pT = randGen.Uniform(0, 20);                      // Uniform pT between 0 and 20 GeV/c
        double phi = randGen.Uniform(-TMath::Pi(), TMath::Pi()); // Uniform phi between -pi and pi
        double eta = randGen.Uniform(-1.2, 1.2);

        lvmother.SetPtEtaPhiM(pT, eta, phi, m_mother);

        // Setup the decay
        TGenPhaseSpace event;
        double masses[2] = {m_daughter1, m_daughter2};

        if (event.SetDecay(lvmother, 2, masses))
        {
            event.Generate(); // Generate the decay event

            // Get daughter particles
            TLorentzVector *ks1 = event.GetDecay(0);
            TLorentzVector *ks2 = event.GetDecay(1);
            TLorentzVector lvmother2 = *ks1 + *ks2; // Reconstruct mother from daughters
            fourVecMother = ROOT::Math::PxPyPzMVector(lvmother2.Px(), lvmother2.Py(), lvmother2.Pz(), lvmother2.M());
            fourVecDau1 = ROOT::Math::PxPyPzMVector(ks1->Px(), ks1->Py(), ks1->Pz(), m_daughter1);
            fourVecDau2 = ROOT::Math::PxPyPzMVector(ks2->Px(), ks2->Py(), ks2->Pz(), m_daughter2);

            h_pT->Fill(pT); // mother pT
            double gen_rapidity = fourVecMother.Y();

            if (abs(gen_rapidity) < 0.5)
            {
                double pt1_gen = ks1->Pt();
                double pt2_gen = ks2->Pt();

                // Apply pT smearing based on K0s efficiency histogram
                double pt1_rec = randGen.Gaus(pt1_gen, pt1_gen * 0.1); // Taking 10 MeV/c (1%) resolution
                double pt2_rec = randGen.Gaus(pt2_gen, pt2_gen * 0.1); // Taking 10 MeV/c (1%) resolution

                // Protect against negative pT
                if (pt1_rec < 0 || pt2_rec < 0)
                    continue;

                hrec_pT->Fill(fourVecMother.Pt());

                // calculate the efficiency of both k0s based on their pT
                // Without pT smearing
                // double eff_ks1 = hEffK0s->GetBinContent(hEffK0s->GetXaxis()->FindBin(pt1_gen));
                // double eff_ks2 = hEffK0s->GetBinContent(hEffK0s->GetXaxis()->FindBin(pt2_gen));

                // using pT smearing
                double eff_ks1 = hEffK0s->GetBinContent(hEffK0s->GetXaxis()->FindBin(pt1_rec));
                double eff_ks2 = hEffK0s->GetBinContent(hEffK0s->GetXaxis()->FindBin(pt2_rec));
                double combined_eff = eff_ks1 * eff_ks2;

                if (fourVecMother.Pt() >= 1 && fourVecMother.Pt() < 2)
                {
                    EffpT1.push_back(combined_eff);
                }
                else if (fourVecMother.Pt() >= 2 && fourVecMother.Pt() < 3)
                {
                    EffpT2.push_back(combined_eff);
                }
                else if (fourVecMother.Pt() >= 3 && fourVecMother.Pt() < 5)
                {
                    EffpT3.push_back(combined_eff);
                }
                else if (fourVecMother.Pt() >= 5 && fourVecMother.Pt() < 7)
                {
                    EffpT4.push_back(combined_eff);
                }
                else if (fourVecMother.Pt() >= 7 && fourVecMother.Pt() < 10)
                {
                    EffpT5.push_back(combined_eff);
                }
                else if (fourVecMother.Pt() >= 10 && fourVecMother.Pt() < 15)
                {
                    EffpT6.push_back(combined_eff);
                }
            }
        }
        lvmother.Clear();
    }

    // Calculate the efficiency in each pT bin
    auto calculateAverageEfficiency = [](const vector<double> &efficiencies)
    {
        if (efficiencies.empty())
            return 0.0;
        double sum = 0.0;
        for (double eff : efficiencies)
        {
            sum += eff;
        }
        return sum / efficiencies.size();
    };

    string savepath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra/pTSmearing/";
    TFile *fout = new TFile((savepath + "f1710_effToy_smeared.root").c_str(), "recreate");

    double avgEff1 = calculateAverageEfficiency(EffpT1);
    double avgEff2 = calculateAverageEfficiency(EffpT2);
    double avgEff3 = calculateAverageEfficiency(EffpT3);
    double avgEff4 = calculateAverageEfficiency(EffpT4);
    double avgEff5 = calculateAverageEfficiency(EffpT5);
    double avgEff6 = calculateAverageEfficiency(EffpT6);

    TH1F *hEffResults = new TH1F("hEffResults", "Acceptance #times Efficiency of f_{0}(1710); #it{p}_{T} bin (GeV/#it{c}); Acceptance #times Efficiency", nPtBins, ptBins);
    hEffResults->SetBinContent(1, avgEff1);
    hEffResults->SetBinContent(2, avgEff2);
    hEffResults->SetBinContent(3, avgEff3);
    hEffResults->SetBinContent(4, avgEff4);
    hEffResults->SetBinContent(5, avgEff5);
    hEffResults->SetBinContent(6, avgEff6);

    // TCanvas *cGenpT = new TCanvas("cGenpT", "Generated p_{T} Distribution", 720, 720);
    // SetCanvasStyle(cGenpT, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(h_pT);
    // h_pT->GetYaxis()->SetMaxDigits(3);
    // h_pT->Draw();

    TCanvas *cRecpT = new TCanvas("cRecpT", "Reconstructed #it{p}_{T} Distribution", 720, 720);
    SetCanvasStyle(cRecpT, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hrec_pT);
    hrec_pT->GetYaxis()->SetMaxDigits(3);
    hrec_pT->Draw();
    hrec_pT->Write("recpt");
    cRecpT->SaveAs((savepath + "f1710_reconstructed_pT_distribution.png").c_str());

    TCanvas *cEffResults = new TCanvas("cEffResults", "Average Efficiency vs #it{p}_{T} bin", 720, 720);
    SetCanvasStyle(cEffResults, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hEffResults);
    hEffResults->GetYaxis()->SetMaxDigits(3);
    hEffResults->Draw();
    hEffResults->Write("eff");
    // hEffK0s->SetLineColor(kRed);
    // hEffK0s->SetMarkerColor(kRed);
    // hEffK0s->Draw("SAME");
    cEffResults->SaveAs((savepath + "f1710_Average_Efficiency_vs_pT_bin.png").c_str());

    TCanvas *cEffK0s = new TCanvas("cEffK0s", "K0s Efficiency vs pT", 720, 720);
    SetCanvasStyle(cEffK0s, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hEffK0s);
    hEffK0s->GetYaxis()->SetMaxDigits(3);
    hEffK0s->Draw("HIST");
    hEffK0s->Write("K0s_eff");
    hOriginalK0sEff->SetLineColor(kRed);
    hOriginalK0sEff->SetMarkerColor(kRed);
    hOriginalK0sEff->Draw("SAME");
    cEffK0s->SaveAs((savepath + "K0s_Efficiency_vs_pT.png").c_str());
}