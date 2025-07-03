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
    TH1F *h_rapdity = new TH1F("h_rapidity", "Rapidity distribution of mother; Rapidity; Events", 100, -3.0, 3.0);
    TH2F *hptvsrap = new TH2F("hptvsrap", "p_{T} vs Rapidity; p_{T} (GeV/c); Rapidity", 150, 0, 30, 100, -3.0, 3.0);
    TH1F *hcosThetaStar = new TH1F("h_cosThetaStar", "Cosine of the angle between mother and daughter in CM frame; cos(#theta^{*}); Events", 20, -1, 1);

    // Loop over multiple decay events
    int nEvents = 1e7;
    TLorentzVector lvmother;
    ROOT::Math::PxPyPzMVector fourVecMother, fourVecDau1, fourVecDau2, fourVecDauCM;
    ROOT::Math::XYZVector threeVecDauCM, threeVecMother;
    for (int i = 0; i < nEvents; i++)
    {
        // Generate random transverse momentum pT
        double pT = randGen.Uniform(0, 30);               // Uniform pT between 0 and 30 GeV/c
        double phi = randGen.Uniform(0, 2 * TMath::Pi()); //  Uniform phi between 0 and 2pi
        double eta = randGen.Uniform(-0.8, 0.8);          // according to ALICE TPC acceptance

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
            fourVecDauCM = boost(fourVecDau1); // Boost daughter momentum to the center of mass frame
            threeVecDauCM = fourVecDauCM.Vect(); // Get the 3-vector of daughter in the frame of mother
            threeVecMother = fourVecMother.Vect(); // Get the 3-vector of mother
            auto cosThetaStar = threeVecMother.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(threeVecMother.Mag2()));
            
            // Fill histogram with transverse momentum of mother
            h_pT->Fill(pT);
            double eta1 = p1->Eta();
            double phi1 = p1->Phi();
            double eta2 = p2->Eta();
            double phi2 = p2->Phi();
            // double angsep = sqrt(pow(eta1 - eta2, 2) + pow(phi1 - phi2, 2));
            // h_angdist->Fill(angsep);
            
            double rec_rapidity = lvmother2.Rapidity();
            h_rapdity->Fill(rec_rapidity);
            if (abs(rec_rapidity) < 0.5)
            {
                hrec_pT->Fill(lvmother2.Pt());
                hptvsrap->Fill(lvmother2.Pt(), rec_rapidity);
                hcosThetaStar->Fill(cosThetaStar); // Fill histogram with cos(theta*)
            }
        }
        lvmother.Clear();
    }

    // Draw histogram
    TCanvas *c1 = new TCanvas("c1", "Decay Simulation", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(h_pT);
    h_pT->GetYaxis()->SetMaxDigits(3);
    h_pT->Draw();
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
    TCanvas *c4 = new TCanvas("c4", "Rapidity Distribution", 720, 720);
    SetCanvasStyle(c4, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(h_rapdity);
    h_rapdity->GetYaxis()->SetMaxDigits(3);
    h_rapdity->Draw();
    TCanvas *c5 = new TCanvas("c5", "pT vs Rapidity", 720, 720);
    SetCanvasStyle(c5, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hptvsrap);
    hptvsrap->GetYaxis()->SetMaxDigits(3);
    hptvsrap->Draw("colz");

    TCanvas *c6 = new TCanvas("c6", "Cosine of the angle between mother and daughter in CM frame", 720, 720);
    SetCanvasStyle(c6, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hcosThetaStar);
    hcosThetaStar->GetYaxis()->SetMaxDigits(3);
    hcosThetaStar->GetYaxis()->SetTitle("Counts");
    hcosThetaStar->GetXaxis()->SetTitle("cos(#theta^{*})");
    hcosThetaStar->GetXaxis()->SetRangeUser(-1, 1);
    hcosThetaStar->Draw();
}
