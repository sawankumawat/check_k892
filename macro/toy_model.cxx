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

    // Define masses (GeV/c^2)
    double m_mother = 1.713;    // Mass of f1710
    double m_daughter1 = 0.493; //  daughter 1 (Ks mass)
    double m_daughter2 = 0.493; //  daughter 2 (Ks mass)

    // Histogram for decay product pT
    TH1F *h_pT = new TH1F("h_pT", "p_{T} distribution of mother; p_{T} (GeV/c); Events", 150, 0, 30);
    TH1F *h_angdist = new TH1F("h_angdist", "Angular distribution of daughter particles; Angular separation (rad); Events", 500, 0, 3 * TMath::Pi());

    // Loop over multiple decay events
    int nEvents = 1e6;
    TLorentzVector lvmother;
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

            // Fill histogram with transverse momentum of mother
            h_pT->Fill(pT);
            double eta1 = p1->Eta();
            double phi1 = p1->Phi();
            double eta2 = p2->Eta();
            double phi2 = p2->Phi();
            double angsep = sqrt(pow(eta1 - eta2, 2) + pow(phi1 - phi2, 2));
            h_angdist->Fill(angsep);
        }
        lvmother.Clear();
    }

    // // Draw histogram
    // TCanvas *c1 = new TCanvas("c1", "Decay Simulation", 720, 720);
    // SetCanvasStyle(c1, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(h_pT);
    // h_pT->GetYaxis()->SetMaxDigits(3);
    // h_pT->Draw();
    TCanvas *c2 = new TCanvas("c2", "Angular Separation", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(h_angdist);
    h_angdist->GetYaxis()->SetMaxDigits(3);
    h_angdist->Draw();
}
