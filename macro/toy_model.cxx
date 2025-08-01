#include <TGenPhaseSpace.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "src/style.h"

using namespace std;

void toy_model()
{
    TRandom3 randGen(0); // Random seed

    gStyle->SetOptStat(0);
    // Define masses (GeV/c^2)
    double m_mother = 1.713;    // Mass of f1710
    double m_daughter1 = 0.493; // daughter 1 (Ks mass)
    double m_daughter2 = 0.493; // daughter 2 (Ks mass)

    // Histogram for decay product pT
    TH1F *h_pT = new TH1F("h_pT", "p_{T} distribution of mother; p_{T} (GeV/c); Events", 150, 0, 30);
    TH1F *hrec_pT = new TH1F("hrec_pT", "Reconstructed p_{T} distribution of mother; p_{T} (GeV/c); Events", 150, 0, 30);
    TH1F *h_angdist = new TH1F("h_angdist", "Angular distribution of daughter particles; Angular separation (rad); Events", 400, 0, 4);
    TH1F *h_rapidity = new TH1F("h_rapidity", "Rapidity distribution of mother; Rapidity; Events", 100, -1.0, 1.0);
    TH2F *hptvsrap = new TH2F("hptvsrap", "p_{T} vs Rapidity; p_{T} (GeV/c); Rapidity", 150, 0, 30, 100, -3.0, 3.0);
    TH1F *hcosThetaStar = new TH1F("h_cosThetaStar", "Cosine of the angle between mother and daughter in CM frame; cos(#theta^{*}); Events", 20, -1, 1);
    TH1F *hpseudoRapidity = new TH1F("h_pseudoRapidity", "Pseudo-rapidity distribution of mother; #eta; Events", 100, -1.0, 1.0);

    // Loop over multiple decay events
    int nEvents = 1e5;
    TLorentzVector lvmother;
    ROOT::Math::PxPyPzMVector fourVecMother, fourVecDau1, fourVecDau2, fourVecDauCM;
    ROOT::Math::XYZVector threeVecDauCM, threeVecMother;
    for (int i = 0; i < nEvents; i++)
    {
        // Generate random transverse momentum pT
        double pT = randGen.Uniform(0, 20);               // Uniform pT between 0 and 20 GeV/c
        double phi = randGen.Uniform(0, 2 * TMath::Pi()); //  Uniform phi between 0 and 2pi
        double eta = randGen.Uniform(-0.8, 0.8);          // according to ALICE TPC acceptance
        // double eta = randGen.Uniform(-1.0, 1.0);

        // Define the 4-momentum of the mother particle
        lvmother.SetPtEtaPhiM(pT, eta, phi, m_mother);

        // Setup the decay
        TGenPhaseSpace event;
        double masses[2] = {m_daughter1, m_daughter2};

        if (event.SetDecay(lvmother, 2, masses))
        {
            event.Generate(); // Generate the decay event

            // Get daughter particles
            TLorentzVector *p1 = event.GetDecay(0);
            TLorentzVector *p2 = event.GetDecay(1);
            TLorentzVector lvmother2 = *p1 + *p2; // Reconstruct mother from daughters
            fourVecMother = ROOT::Math::PxPyPzMVector(lvmother2.Px(), lvmother2.Py(), lvmother2.Pz(), lvmother2.M());
            fourVecDau1 = ROOT::Math::PxPyPzMVector(p1->Px(), p1->Py(), p1->Pz(), m_daughter1);
            fourVecDau2 = ROOT::Math::PxPyPzMVector(p2->Px(), p2->Py(), p2->Pz(), m_daughter2);
            ROOT::Math::Boost boost{fourVecMother.BoostToCM()}; // Boost to center of mass frame
            fourVecDauCM = boost(fourVecDau1);                  // Boost daughter momentum to the center of mass frame

            auto cosThetaStar = fourVecMother.Vect().Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(fourVecMother.Vect().Mag2()));

            // Fill histogram with transverse momentum of mother
            h_pT->Fill(pT);
            double eta1 = p1->Eta();
            double phi1 = p1->Phi();
            double eta2 = p2->Eta();
            double phi2 = p2->Phi();
            // double angsep = sqrt(pow(eta1 - eta2, 2) + pow(phi1 - phi2, 2));
            // h_angdist->Fill(angsep);

            // if (abs(p1->Eta()) > 0.8 || abs(p2->Eta()) > 0.8)
            // {
            //     continue;
            // }

            // if (abs(p1->Y()) > 0.5 || abs(p2->Y()) > 0.5)
            // {
            //     continue;
            // }

            double gen_rapidity = fourVecMother.Y();
            hpseudoRapidity->Fill(fourVecMother.Eta());
            if (abs(gen_rapidity) < 0.5)
            {
                h_rapidity->Fill(fourVecMother.Y());
                hrec_pT->Fill(fourVecMother.Pt());
                hptvsrap->Fill(fourVecMother.Pt(), gen_rapidity);
                hcosThetaStar->Fill(cosThetaStar); // Fill histogram with cos(theta*)
            }
        }
        lvmother.Clear();
    }

    // Draw histogram
    TCanvas *c1 = new TCanvas("c1", "p_{T} without rapidity cut", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(h_pT);
    h_pT->GetYaxis()->SetMaxDigits(3);
    h_pT->Draw();
    c1->SaveAs("toy_model_plots/pT_distribution.png");
    // TCanvas *c2 = new TCanvas("c2", "Angular Separation", 720, 720);
    // SetCanvasStyle(c2, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(h_angdist);
    // h_angdist->GetYaxis()->SetMaxDigits(3);
    // h_angdist->Draw();
    // c2->SaveAs("angular_separation.png");

    TCanvas *c3 = new TCanvas("c3", "Reconstructed pT", 720, 720);
    SetCanvasStyle(c3, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hrec_pT);
    hrec_pT->GetYaxis()->SetMaxDigits(3);
    hrec_pT->Draw();
    c3->SaveAs("toy_model_plots/reconstructed_pT_distribution.png");

    TCanvas *c4 = new TCanvas("c4", "Rapidity Distribution", 720, 720);
    SetCanvasStyle(c4, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(h_rapidity);
    h_rapidity->GetYaxis()->SetMaxDigits(3);
    h_rapidity->Draw();
    c4->SaveAs("toy_model_plots/rapidity_distribution.png");

    TCanvas *c5 = new TCanvas("c5", "pT vs Rapidity", 720, 720);
    SetCanvasStyle(c5, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hptvsrap);
    hptvsrap->GetYaxis()->SetMaxDigits(3);
    hptvsrap->Draw("colz");
    c5->SaveAs("toy_model_plots/pT_vs_rapidity.png");

    TCanvas *c6 = new TCanvas("c6", "Cosine of the angle between mother and daughter in CM frame", 720, 720);
    SetCanvasStyle(c6, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hcosThetaStar);
    hcosThetaStar->GetYaxis()->SetMaxDigits(3);
    hcosThetaStar->GetYaxis()->SetTitle("Counts");
    hcosThetaStar->GetXaxis()->SetTitle("cos(#theta)");
    hcosThetaStar->GetXaxis()->SetRangeUser(-1, 1);
    hcosThetaStar->SetMinimum(0);
    hcosThetaStar->SetMaximum(1.2 * hcosThetaStar->GetMaximum());
    hcosThetaStar->Draw();
    c6->SaveAs("toy_model_plots/cosThetaStar_distribution.png");

    TCanvas *c7 = new TCanvas("c7", "Pseudo-rapidity Distribution", 720, 720);
    SetCanvasStyle(c7, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hpseudoRapidity);
    hpseudoRapidity->GetYaxis()->SetMaxDigits(3);
    hpseudoRapidity->GetYaxis()->SetTitle("Counts");
    hpseudoRapidity->GetXaxis()->SetTitle("#eta");
    hpseudoRapidity->GetXaxis()->SetRangeUser(-1, 1);
    hpseudoRapidity->SetMinimum(0);
    hpseudoRapidity->SetMaximum(1.2 * hpseudoRapidity->GetMaximum());
    hpseudoRapidity->Draw();
    c7->SaveAs("toy_model_plots/pseudoRapidity_distribution.png");
}
