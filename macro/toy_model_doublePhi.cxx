#include <TGenPhaseSpace.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "src/style.h"

#include <Math/VectorUtil.h>

using namespace std;

static double deltaR(double phiA, double etaA, double phiB, double etaB)
{
    return sqrt(pow(etaA - etaB, 2) + pow(TVector2::Phi_mpi_pi(phiA - phiB), 2));
}

void toy_model_doublePhi()
{
    TRandom3 randGen(0); // Random seed

    gStyle->SetOptStat(0);
    double m_mother = 2.714;     // Mass of T4s2++(2714)
    double m_daughter1 = 1.0194; // daughter 1 (Phi mass)
    double m_daughter2 = 1.0194; // daughter 2 (Phi mass)
    double m_kaon = 0.49367;     // charged kaon mass
    double mother_width = 0.012; // GeV (12 MeV)
    double phi_width = 0.0045;    // GeV (4.5 MeV)

    TFile *fout = new TFile("toy_model_plots/doublePhi/toyModel_doublePhi.root", "recreate");

    // Histogram for decay product pT
    TH1F *h_pT = new TH1F("h_pT", "p_{T} distribution of mother; p_{T} (GeV/c); Events", 100, 0, 10);
    TH1F *hrec_pT = new TH1F("hrec_pT", "Reconstructed p_{T} distribution of mother; p_{T} (GeV/c); Events", 150, 0, 30);
    TH1F *hPhiPhiPtCorr = new TH1F("hPhiPhiPtCorr", "pt correlaiton", 100, 0, 10);
    TH1F *hPhiPhiCorrA = new TH1F("hPhiPhiCorrA", "PhiPhiCorrA", 200, 0, 1);
    TH1F *hPhiPhiCorrZ = new TH1F("hPhiPhiCorrZ", "PhiPhiCorrZ", 200, 0, 1);
    TH1F *hDeltaR = new TH1F("hDeltaR", "Delta R between two daughters", 200, 0, 4);
    TH2F *hDeltaRvsPt = new TH2F("hDeltaRvsPt", "Delta R vs pT of mother; p_{T} (GeV/c); #DeltaR", 100, 0, 20, 200, 0, 4);
    TH1F *hPtPt = new TH1F("hPtPt", "#it{p}_{T}(#phi_{1}) x #it{p}_{T}(#phi_{2}); #it{p}_{T}(#phi_{1}) x #it{p}_{T}(#phi_{2}); Counts", 500, 0, 5);
    TH1F *hDeltaR_Kaons = new TH1F("hDeltaR_Kaons", "Delta R between kaons from different phis", 200, 0, 4);
    TH1F *hDeltaRkaonplus = new TH1F("hDeltaRkaonplus", "Delta R between K+ from two phis", 200, 0, 4);
    TH1F *hDeltaRkaonminus = new TH1F("hDeltaRkaonminus", "Delta R between K- from two phis", 200, 0, 4);

    // helper deltaR function defined at file scope to avoid name collisions

    // Loop over multiple decay events
    int nEvents = 1e6;
    TLorentzVector lvmother;
    ROOT::Math::PxPyPzMVector fourVecMother, fourVecDau1, fourVecDau2;
    ROOT::Math::XYZVector threeVecDauCM, threeVecMother;
    TF1 ptdistribution("ptSpectra", "x * TMath::Exp(-TMath::Sqrt(0.13957*0.13957 + x*x) / 0.4)", 0.0, 30.0);

    for (int i = 0; i < nEvents; i++)
    {
        // Generate random transverse momentum pT
        double pT = ptdistribution.GetRandom();           // pT from thermal-like distribution
        double phi = randGen.Uniform(0, 2 * TMath::Pi()); //  Uniform phi between -pi and pi
        double eta = randGen.Uniform(-0.8, 0.8);          // according to ALICE TPC acceptance
        // double eta = randGen.Uniform(-2, 2);

        // sample mother mass with a Breit-Wigner (small width)
        double m_mother_sample = randGen.BreitWigner(m_mother, mother_width * 0.5);
        if (m_mother_sample <= 0) m_mother_sample = m_mother;
        lvmother.SetPtEtaPhiM(pT, eta, phi, m_mother_sample);
        // lvmother.SetPxPyPzE(px, py, pz, E);

        // Setup the decay with sampled daughter masses (phi width)
        TGenPhaseSpace event;
        double m_phi1 = randGen.BreitWigner(m_daughter1, phi_width * 0.5);
        double m_phi2 = randGen.BreitWigner(m_daughter2, phi_width * 0.5);
        // ensure physical masses and above KK threshold for subsequent decays
        const double kkThreshold = 2.0 * m_kaon + 1e-6;
        if (m_phi1 < kkThreshold) m_phi1 = kkThreshold;
        if (m_phi2 < kkThreshold) m_phi2 = kkThreshold;
        double masses[2] = {m_phi1, m_phi2};

        if (event.SetDecay(lvmother, 2, masses))
        {
            event.Generate(); // Generate the decay event

            // Get daughter particles
            TLorentzVector *phi1 = event.GetDecay(0);
            TLorentzVector *dphi2 = event.GetDecay(1);
            TLorentzVector lvmother2 = *phi1 + *dphi2; // Reconstruct mother from daughters
            fourVecMother = ROOT::Math::PxPyPzMVector(lvmother2.Px(), lvmother2.Py(), lvmother2.Pz(), lvmother2.M());
            fourVecDau1 = ROOT::Math::PxPyPzMVector(phi1->Px(), phi1->Py(), phi1->Pz(), m_daughter1);
            fourVecDau2 = ROOT::Math::PxPyPzMVector(dphi2->Px(), dphi2->Py(), dphi2->Pz(), m_daughter2);

            double pTcorr = fourVecDau1.Pt() / (fourVecMother.Pt() - fourVecDau1.Pt());
            hPhiPhiPtCorr->Fill(pTcorr);

            double pTcorrA = abs(fourVecDau1.Pt() - fourVecDau2.Pt()) / (fourVecDau1.Pt() + fourVecDau2.Pt());
            hPhiPhiCorrA->Fill(pTcorrA);

            double pTcorrZ = fourVecDau1.Pt() / (fourVecDau1.Pt() + fourVecDau2.Pt());
            hPhiPhiCorrZ->Fill(pTcorrZ);

            double drPhi = sqrt(pow(phi1->Eta() - dphi2->Eta(), 2) + pow(TVector2::Phi_mpi_pi(phi1->Phi() - dphi2->Phi()), 2));
            hDeltaR->Fill(drPhi);

            hDeltaRvsPt->Fill(fourVecMother.Pt(), drPhi);
            hPtPt->Fill(phi1->Pt() * dphi2->Pt());
            // hPtPt->Fill(phi1->Pt());

            // Fill histogram with transverse momentum of mother
            h_pT->Fill(fourVecDau1.Pt());

            // Decay each phi -> K+ K- and fill deltaR between kaons
            TGenPhaseSpace eventPhi1, eventPhi2;
            double kaonMasses[2] = {m_kaon, m_kaon};
            TLorentzVector phi1LV = *phi1;
            TLorentzVector phi2LV = *dphi2;
            bool okPhi1 = false, okPhi2 = false;
            TLorentzVector k11, k12, k21, k22;

            if (eventPhi1.SetDecay(phi1LV, 2, kaonMasses))
            {
                eventPhi1.Generate();
                TLorentzVector *tk1 = eventPhi1.GetDecay(0);
                TLorentzVector *tk2 = eventPhi1.GetDecay(1);
                k11 = *tk1;
                k12 = *tk2;
                // same-phi kaon deltaR can be filled later if needed
                okPhi1 = true;
            }

            if (eventPhi2.SetDecay(phi2LV, 2, kaonMasses))
            {
                eventPhi2.Generate();
                TLorentzVector *tk1 = eventPhi2.GetDecay(0);
                TLorentzVector *tk2 = eventPhi2.GetDecay(1);
                k21 = *tk1;
                k22 = *tk2;
                // same-phi kaon deltaR can be filled later if needed
                okPhi2 = true;
            }

            // If both phis decayed to kaons, fill deltaR between all inter-phi kaon combinations
            if (okPhi1 && okPhi2)
            {
                ROOT::Math::PtEtaPhiMVector kplus1_v(k11.Pt(), k11.Eta(), k11.Phi(), k11.M());
                ROOT::Math::PtEtaPhiMVector kminus1_v(k12.Pt(), k12.Eta(), k12.Phi(), k12.M());
                ROOT::Math::PtEtaPhiMVector kplus2_v(k21.Pt(), k21.Eta(), k21.Phi(), k21.M());
                ROOT::Math::PtEtaPhiMVector kminus2_v(k22.Pt(), k22.Eta(), k22.Phi(), k22.M());

                const auto minKaonDeltaR = [&](const ROOT::Math::PtEtaPhiMVector &kplus1,
                                               const ROOT::Math::PtEtaPhiMVector &kminus1,
                                               const ROOT::Math::PtEtaPhiMVector &kplus2,
                                               const ROOT::Math::PtEtaPhiMVector &kminus2)
                {
                    const double dRkplus = deltaR(kplus1.Phi(), kplus1.Eta(), kplus2.Phi(), kplus2.Eta());
                    const double dRkminus = deltaR(kminus1.Phi(), kminus1.Eta(), kminus2.Phi(), kminus2.Eta());

                    hDeltaRkaonplus->Fill(dRkplus);
                    hDeltaRkaonminus->Fill(dRkminus);

                    double minDR = dRkplus;
                    minDR = std::min(minDR, dRkminus);
                    minDR = std::min(minDR, deltaR(kplus1.Phi(), kplus1.Eta(), kminus1.Phi(), kminus1.Eta()));
                    minDR = std::min(minDR, deltaR(kplus1.Phi(), kplus1.Eta(), kminus2.Phi(), kminus2.Eta()));
                    minDR = std::min(minDR, deltaR(kplus2.Phi(), kplus2.Eta(), kminus1.Phi(), kminus1.Eta()));
                    minDR = std::min(minDR, deltaR(kplus2.Phi(), kplus2.Eta(), kminus2.Phi(), kminus2.Eta()));

                    return minDR;
                };

                double minDR = minKaonDeltaR(kplus1_v, kminus1_v, kplus2_v, kminus2_v);
                hDeltaR_Kaons->Fill(minDR);
            }

            double eta1_v = phi1->Eta();
            double phi1_v = phi1->Phi();
            double eta2_v = dphi2->Eta();
            double phi2_v = dphi2->Phi();

            double gen_rapidity = fourVecMother.Y();

            // if (abs(gen_rapidity) < 0.5)
            // {
            //     hrec_pT->Fill(fourVecMother.Pt());
            // }
        }
        lvmother.Clear();
    }

    // // Draw histogram
    // TCanvas *cpT = new TCanvas("cpT", "p_{T} distribution of Phi", 720, 720);
    // SetCanvasStyle(cpT, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(h_pT);
    // h_pT->GetYaxis()->SetMaxDigits(3);
    // h_pT->Draw();
    // // cpT->SaveAs("toy_model_plots/doublePhi/MotherpT_distribution.png");

    TCanvas *cPhiPhiPtCorr = new TCanvas("cPhiPhiPtCorr", "Pt correlation", 720, 720);
    SetCanvasStyle(cPhiPhiPtCorr, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hPhiPhiPtCorr);
    hPhiPhiPtCorr->GetYaxis()->SetMaxDigits(3);
    hPhiPhiPtCorr->GetXaxis()->SetTitle("p_{T1} / (p_{T mother} - p_{T1})");
    hPhiPhiPtCorr->GetYaxis()->SetTitle("Counts");
    hPhiPhiPtCorr->Draw();
    hPhiPhiPtCorr->Write("hPtCorr");
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.SetTextFont(42);
    int maxbin = hPhiPhiPtCorr->GetMaximumBin();
    double maxXvalue = hPhiPhiPtCorr->GetXaxis()->GetBinCenter(maxbin);
    lat.DrawLatex(0.7, 0.82, Form("Peak %.1f", maxXvalue));
    lat.DrawLatex(0.7, 0.75, Form("Mean %.1f", hPhiPhiPtCorr->GetMean()));
    cPhiPhiPtCorr->SaveAs("toy_model_plots/doublePhi/phi_phi_mass_correlation.png");

    TCanvas *cPhiPhiCorrA = new TCanvas("cPhiPhiCorrA", "Pt correlation A", 720, 720);
    SetCanvasStyle(cPhiPhiCorrA, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hPhiPhiCorrA);
    hPhiPhiCorrA->GetYaxis()->SetMaxDigits(3);
    hPhiPhiCorrA->GetXaxis()->SetTitle("|p_{T1} - p_{T2}| / (p_{T1} + p_{T2})");
    hPhiPhiCorrA->GetYaxis()->SetTitle("Counts");
    hPhiPhiCorrA->Draw();
    hPhiPhiCorrA->Write("hPtCorrA");
    cPhiPhiCorrA->SaveAs("toy_model_plots/doublePhi/phi_phi_correlation_A.png");

    TCanvas *cPhiPhiCorrZ = new TCanvas("cPhiPhiCorrZ", "Pt correlation Z", 720, 720);
    SetCanvasStyle(cPhiPhiCorrZ, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hPhiPhiCorrZ);
    hPhiPhiCorrZ->GetYaxis()->SetMaxDigits(3);
    hPhiPhiCorrZ->GetXaxis()->SetTitle("p_{T1} / (p_{T1} + p_{T2})");
    hPhiPhiCorrZ->GetYaxis()->SetTitle("Counts");
    hPhiPhiCorrZ->Draw();
    hPhiPhiCorrZ->Write("hPtCorrZ");
    cPhiPhiCorrZ->SaveAs("toy_model_plots/doublePhi/phi_phi_correlation_Z.png");

    TCanvas *cDeltaRPhiPhi = new TCanvas("cDeltaRPhiPhi", "Delta R between two daughters", 720, 720);
    SetCanvasStyle(cDeltaRPhiPhi, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hDeltaR);
    hDeltaR->GetYaxis()->SetMaxDigits(3);
    hDeltaR->GetXaxis()->SetTitle("#DeltaR = #sqrt{(#Delta#eta)^{2} + (#Delta#varphi)^{2}}");
    hDeltaR->GetYaxis()->SetTitle("Counts");
    hDeltaR->Draw();
    cDeltaRPhiPhi->SaveAs("toy_model_plots/doublePhi/deltaR_between_daughters.png");

    // TCanvas *cDeltaRPhiPhivsPt = new TCanvas("cDeltaRPhiPhivsPt", "Delta R vs pT of mother", 1440, 720);
    // SetCanvasStyle(cDeltaRPhiPhivsPt, 0.12, 0.15, 0.05, 0.12);
    // cDeltaRPhiPhivsPt->Divide(4, 3);
    // float pTbins[] = {0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0};
    // for (int i = 0; i < sizeof(pTbins) / sizeof(pTbins[0]) - 1; i++)
    // {
    //     cDeltaRPhiPhivsPt->cd(i + 1);
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
    // // cDeltaRPhiPhivsPt->SaveAs("toy_model_plots/doublePhi/deltaR_vs_pTofMother.png");

    TCanvas *cminDeltaRKaons = new TCanvas("cminDeltaRKaons", "Delta R between kaons", 720, 720);
    SetCanvasStyle(cminDeltaRKaons, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hDeltaR_Kaons);
    hDeltaR_Kaons->GetYaxis()->SetMaxDigits(3);
    hDeltaR_Kaons->GetXaxis()->SetTitle("Min. #DeltaR");
    hDeltaR_Kaons->GetYaxis()->SetTitle("Counts");
    hDeltaR_Kaons->Draw();
    cminDeltaRKaons->SaveAs("toy_model_plots/doublePhi/deltaR_between_kaons.png");
}