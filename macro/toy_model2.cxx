#include <TGenPhaseSpace.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "src/style.h"

#include <Math/VectorUtil.h>

using namespace std;

void toy_model2()
{
    TRandom3 randGen(0); // Random seed

    gStyle->SetOptStat(0);
    // Define masses (GeV/c^2)
    double m_mother = 2.714;     // Mass of f1710
    double m_daughter1 = 1.0194; // daughter 1 (Ks mass)
    double m_daughter2 = 1.0194; // daughter 2 (Ks mass)
    double m_pion = 0.13957;     // charged pion mass
    double m_kaon = 0.49367;     // charged kaon mass

    // Histogram for decay product pT
    TH1F *h_pT = new TH1F("h_pT", "p_{T} distribution of mother; p_{T} (GeV/c); Events", 150, 0, 30);
    TH1F *hrec_pT = new TH1F("hrec_pT", "Reconstructed p_{T} distribution of mother; p_{T} (GeV/c); Events", 150, 0, 30);
    TH1F *hPhiPhiPtCorr = new TH1F("hPhiPhiPtCorr", "pt correlaiton", 100, 0, 10);
    TH1F *hDeltaR = new TH1F("hDeltaR", "Delta R between two daughters", 200, 0, 4);
    TH2F *hDeltaRvsPt = new TH2F("hDeltaRvsPt", "Delta R vs pT of mother; p_{T} (GeV/c); #DeltaR", 100, 0, 20, 200, 0, 4);

    // Loop over multiple decay events
    int nEvents = 1e7;
    TLorentzVector lvmother;
    ROOT::Math::PxPyPzMVector fourVecMother, fourVecDau1, fourVecDau2;
    ROOT::Math::XYZVector threeVecDauCM, threeVecMother;
    TF1 ptdistribution("ptSpectra",
                       "x * TMath::Exp(-TMath::Sqrt(0.13957*0.13957 + x*x) / 0.4)",
                       0.1, 20.0); 

    for (int i = 0; i < nEvents; i++)
    {
        // Generate random transverse momentum pT
        double pT = ptdistribution.GetRandom();           // pT from thermal-like distribution
        double phi = randGen.Uniform(0, 2 * TMath::Pi()); //  Uniform phi between -pi and pi
        double eta = randGen.Uniform(-0.8, 0.8);          // according to ALICE TPC acceptance
        // double eta = randGen.Uniform(-2, 2);

        // // Generate in pt, y and phi
        // double pT = randGen.Uniform(0, 20);                      // Uniform pT between 0 and 20 GeV/c
        // double phi = randGen.Uniform(-TMath::Pi(), TMath::Pi()); // Uniform phi between -pi and pi
        // double y = randGen.Uniform(-1.0, 1.0);                   // Uniform rapidity (the desired range)

        // // Compute from pT, y, phi, and mass
        // double mT = TMath::Sqrt(pT * pT + m_mother * m_mother);
        // double px = pT * TMath::Cos(phi);
        // double py = pT * TMath::Sin(phi);
        // double pz = pT * TMath::SinH(eta);
        // double E = std::hypot(std::hypot(pT, pz), m_mother);
        // double Eprime = pT * TMath::CosH(y);
        // double E = TMath::Sqrt(Eprime * Eprime + m_mother * m_mother);
        // double eta = 0.5 * TMath::Log((E + pz) / (E - pz));

        // Define the 4-momentum of the mother particle
        lvmother.SetPtEtaPhiM(pT, eta, phi, m_mother);
        // lvmother.SetPxPyPzE(px, py, pz, E);

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

            double pTcorr = fourVecDau1.Pt() / (fourVecMother.Pt() - fourVecDau1.Pt());
            hPhiPhiPtCorr->Fill(pTcorr);
            double deltaR = sqrt(pow(ks1->Eta() - ks2->Eta(), 2) + pow(TVector2::Phi_mpi_pi(ks1->Phi() - ks2->Phi()), 2));
            hDeltaR->Fill(deltaR);
            hDeltaRvsPt->Fill(fourVecMother.Pt(), deltaR);

            // Fill histogram with transverse momentum of mother
            h_pT->Fill(pT);
            double eta1 = ks1->Eta();
            double phi1 = ks1->Phi();
            double eta2 = ks2->Eta();
            double phi2 = ks2->Phi();

            double gen_rapidity = fourVecMother.Y();

            if (abs(gen_rapidity) < 0.5)
            {
                hrec_pT->Fill(fourVecMother.Pt());
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
    // c1->SaveAs("toy_model_plots/pT_distribution.png");

    TCanvas *c3 = new TCanvas("c3", "Reconstructed pT", 720, 720);
    SetCanvasStyle(c3, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hrec_pT);
    hrec_pT->GetYaxis()->SetMaxDigits(3);
    hrec_pT->Draw();
    // c3->SaveAs("toy_model_plots/reconstructed_pT_distribution.png");

    TCanvas *c4 = new TCanvas("c4", "Pt correlation", 720, 720);
    SetCanvasStyle(c4, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hPhiPhiPtCorr);
    hPhiPhiPtCorr->GetYaxis()->SetMaxDigits(3);
    hPhiPhiPtCorr->GetXaxis()->SetTitle("p_{T1} / (p_{T mother} - p_{T1})");
    hPhiPhiPtCorr->GetYaxis()->SetTitle("Counts");
    hPhiPhiPtCorr->Draw();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.SetTextFont(42);
    int maxbin = hPhiPhiPtCorr->GetMaximumBin();
    double maxXvalue = hPhiPhiPtCorr->GetXaxis()->GetBinCenter(maxbin);
    lat.DrawLatex(0.7, 0.82, Form("Peak %.1f", maxXvalue));
    lat.DrawLatex(0.7, 0.75, Form("Mean %.1f", hPhiPhiPtCorr->GetMean()));
    c4->SaveAs("toy_model_plots/phi_phi_mass_correlation.png");

    TCanvas *c5 = new TCanvas("c5", "Delta R between two daughters", 720, 720);
    SetCanvasStyle(c5, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hDeltaR);
    hDeltaR->GetYaxis()->SetMaxDigits(3);
    hDeltaR->GetXaxis()->SetTitle("#DeltaR = #sqrt{(#Delta#eta)^{2} + (#Delta#varphi)^{2}}");
    hDeltaR->GetYaxis()->SetTitle("Counts");
    hDeltaR->Draw();
    c5->SaveAs("toy_model_plots/deltaR_between_daughters.png");

    // TCanvas *c6 = new TCanvas("c6", "Delta R vs pT of mother", 1440, 720);
    // SetCanvasStyle(c6, 0.12, 0.15, 0.05, 0.12);
    // c6->Divide(4, 3);
    // float pTbins[] = {0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0};
    // for (int i = 0; i < sizeof(pTbins) / sizeof(pTbins[0]) - 1; i++)
    // {
    //     c6->cd(i + 1);
    //     TH1D *hDeltaRbin = hDeltaRvsPt->ProjectionY(Form("hDeltaR_%.1f_%.1f", pTbins[i], pTbins[i + 1]),
    //                                                 hDeltaRvsPt->GetXaxis()->FindBin(pTbins[i]),
    //                                                 hDeltaRvsPt->GetXaxis()->FindBin(pTbins[i + 1]) - 1);
    //     SetHistoQA(hDeltaRbin);
    //     hDeltaRbin->GetYaxis()->SetMaxDigits(3);
    //     hDeltaRbin->GetXaxis()->SetTitle("#DeltaR = #sqrt{(#Delta#eta)^{2} + (#Delta#varphi)^{2}}");
    //     hDeltaRbin->GetYaxis()->SetTitle("Counts");
    //     hDeltaRbin->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c", pTbins[i], pTbins[i + 1]));
    //     hDeltaRbin->Draw();
    // }
    // c6->SaveAs("toy_model_plots/deltaR_vs_pTofMother.png");
}