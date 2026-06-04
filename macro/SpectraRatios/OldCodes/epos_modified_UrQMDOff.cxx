#include <vector>
#include <iostream>
#include "TDatabasePDG.h"
using namespace std;

// =====================================================
// centrality bin from percentile
// =====================================================
int GetCentralityBin(double cent)
{
    if (cent < 1)
        return 0;
    else if (cent < 5)
        return 1;
    else if (cent < 10)
        return 2;
    else if (cent < 15)
        return 3;
    else if (cent < 20)
        return 4;
    else if (cent < 30)
        return 5;
    else if (cent < 40)
        return 6;
    else if (cent < 50)
        return 7;
    else if (cent < 70)
        return 8;
    else
        return 9;
}

// int GetCentralityBin(double cent)
// {
//     if (cent < 0.96)
//         return 0;
//     else if (cent < 4.80)
//         return 1;
//     else if (cent < 9.60)
//         return 2;
//     else if (cent < 14.42)
//         return 3;
//     else if (cent < 19.25)
//         return 4;
//     else if (cent < 28.92)
//         return 5;
//     else if (cent < 38.62)
//         return 6;
//     else if (cent < 48.38)
//         return 7;
//     else if (cent < 68.09)
//         return 8;
//     else
//         return 9;
// }

// =====================================================
// MAIN
// =====================================================
void epos_modified_UrQMDOff()
{
    TDatabasePDG *pdg = new TDatabasePDG();

    // =====================================================
    // READ ALL FILES
    // =====================================================
    TChain chain("teposevent0");
    chain.Add("/home/sawan/Storage/EPOS_localOutputs/mergedEPOS_UrQMD.root");

    Long64_t nEv = chain.GetEntries();
    cout << "Total events = " << nEv << endl;
    nEv = 50000; // for testing

    // =====================================================
    // branches
    // =====================================================
    Int_t np;
    Float_t px[200000], py[200000], pz[200000];
    Float_t e[200000];
    Int_t id[200000];
    Int_t ist[200000];
    Int_t ity[200000];

    chain.SetBranchAddress("np", &np);
    chain.SetBranchAddress("px", px);
    chain.SetBranchAddress("py", py);
    chain.SetBranchAddress("pz", pz);
    chain.SetBranchAddress("e", e);
    chain.SetBranchAddress("id", id);
    chain.SetBranchAddress("ist", ist);
    chain.SetBranchAddress("ity", ity);

    // =====================================================
    // FT0 multiplicity
    // =====================================================
    TH1F *hFT0 = new TH1F("hFT0", "FT0 multiplicity", 500, 0, 500);
    vector<int> eventMult(nEv);
    TH1F *hIST = new TH1F("hIST", "IST distribution", 20, -0.5, 19.5);
    TH1F *hITY = new TH1F("hITY", "ITY distribution", 100, -0.5, 99.5);

    cout << "First pass: computing FT0..." << endl;

    int eposid[] = {120, -120, 130, -130, 1120, -1120};
    int pdgid[] = {211, -211, 321, -321, 2212, -2212};

    // =====================================================
    // FIRST PASS → FT0 multiplicity
    // =====================================================
    for (Long64_t ev = 0; ev < nEv; ev++)
    {
        chain.GetEntry(ev);
        int multFT0 = 0;

        for (int i = 0; i < np; i++)
        {
            int eposID = TMath::Abs(id[i]);
            if (eposID != 120 && eposID != 130 && eposID != 1120)
                continue;

            // int PID;

            // for (int j = 0; j < 6; j++)
            // {
            //     if (eposID == TMath::Abs(eposid[j]))
            //     {
            //         PID = pdgid[j];
            //         break;
            //     }
            // }

            // int charge = pdg->GetParticle(PID)->Charge() / 3.;

            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
            if (pt < 0.1)
                continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));

            if ((eta > 3.5 && eta < 4.9) || (eta > -3.3 && eta < -2.1))
                if (ist[i] == 0) // for final state particles
                                 // if (charge != 0)
                    multFT0++;
        }

        hFT0->Fill(multFT0);
        eventMult[ev] = multFT0;
    }

    // =====================================================
    // multiplicity → percentile
    // =====================================================
    vector<double> multToCent(hFT0->GetNbinsX() + 1, 100);
    double total = hFT0->Integral();
    double cum = 0;

    for (int b = hFT0->GetNbinsX(); b >= 1; b--)
    {
        cum += hFT0->GetBinContent(b);
        multToCent[b] = 100.0 * cum / total;
    }

    // =====================================================
    // storage
    // =====================================================
    double sumDNdEta[10] = {0};
    int nEventsCent[10] = {0};

    double phiYield[10] = {0};
    double kstarYield[10] = {0};
    double kshortYield[10] = {0};
    double protonYield[10] = {0};
    double pionYield[10] = {0};
    double kaonYield[10] = {0};
    double XiStarYield[10] = {0};
    double ChargedKstarYield[10] = {0};

    double sumPtPhi[10] = {0};
    double sumPtKstar[10] = {0};
    double sumPtKshort[10] = {0};
    double sumPtProton[10] = {0};
    double sumPtPion[10] = {0};
    double sumPtKaon[10] = {0};
    double sumPtXiStar[10] = {0};
    double sumPtChargedKstar[10] = {0};

    double countPtPhi[10] = {0};
    double countPtKstar[10] = {0};
    double countPtKshort[10] = {0};
    double countPtProton[10] = {0};
    double countPtPion[10] = {0};
    double countPtKaon[10] = {0};
    double countPtXiStar[10] = {0};
    double countPtChargedKstar[10] = {0};

    TH1F *hPtMBPhi = new TH1F("hPtMBPhi", "#phi p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBKstar = new TH1F("hPtMBKstar", "K*^{0} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBKshort = new TH1F("hPtMBKshort", "K^{0}_{S} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBProton = new TH1F("hPtMBProton", "p/#bar{p} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBPion = new TH1F("hPtMBPion", "#pi^{#pm} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBKaon = new TH1F("hPtMBKaon", "K^{#pm} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBXiStar = new TH1F("hPtMBXiStar", "#Xi(1530) p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBChargedKstar = new TH1F("hPtMBChargedKstar", "K*^{#pm} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);

    cout << "Second pass: physics observables..." << endl;

    // =====================================================
    // SECOND PASS
    // =====================================================
    for (Long64_t ev = 0; ev < nEv; ev++)
    {
        chain.GetEntry(ev);

        int multFT0 = eventMult[ev];
        int bin = hFT0->FindBin(multFT0);
        double cent = multToCent[bin];
        int cbin = GetCentralityBin(cent);

        nEventsCent[cbin]++;

        int nchMid = 0;
        int atLeastOnePiKp_in_ModEta1 = 0;

        for (int i = 0; i < np; i++)
        {
            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));

            int eposID = TMath::Abs(id[i]);

            // ----------------------
            // INEL > 0 selection
            // ----------------------
            if (eposID == 120 || eposID == 130 || eposID == 1120)
            {
                if (fabs(eta) < 1.0 && ist[i] == 0)
                {
                    atLeastOnePiKp_in_ModEta1++;
                }
            }
        }

        for (int i = 0; i < np; i++)
        {
            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
            // if (pt < 0.1)
            //     continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));
            double y = 0.5 * log((e[i] + pz[i]) / (e[i] - pz[i]));

            int eposID = TMath::Abs(id[i]);

            // int PID;
            // int charge = 0;
            // if (eposID == 120 || eposID == 130 || eposID == 1120)
            // {
            //     for (int j = 0; j < 6; j++)
            //     {
            //         if (eposID == TMath::Abs(eposid[j]))
            //         {
            //             PID = pdgid[j];
            //             break;
            //         }
            //     }

            //     charge = pdg->GetParticle(PID)->Charge() / 3.;
            // }

            // ----------------------
            // charged multiplicity
            // ----------------------
            if (eposID == 120 || eposID == 130 || eposID == 1120)
                if (fabs(eta) < 0.5)
                    if (ist[i] == 0)
                        // if (charge != 0)
                        if (pt >= 0.1)
                            nchMid++;

            // =================================================
            // URQMD OFF → final hadrons
            // =================================================
            hIST->Fill(ist[i]);
            hITY->Fill(ity[i]);

            // if (fabs(y) > 0.5)
            //     continue;

            // ----------------------
            // Kshort
            // ----------------------
            if (TMath::Abs(id[i]) == 20 && ist[i] == 0 && (fabs(y) < 0.5))
            {
                // if (atLeastOnePiKp_in_ModEta1 > 0) // INEL > 0
                {
                    kshortYield[cbin]++;
                    sumPtKshort[cbin] += pt;
                    countPtKshort[cbin]++;
                    hPtMBKshort->Fill(pt);
                }
            }

            // -----------------------
            // proton, pion, kaon
            // -----------------------
            if (eposID == 120 && ist[i] == 0 && (fabs(y) < 0.5))
            {
                pionYield[cbin]++;
                sumPtPion[cbin] += pt;
                countPtPion[cbin]++;
                hPtMBPion->Fill(pt);
            }

            if (eposID == 130 && ist[i] == 0 && (fabs(y) < 0.5))
            {
                kaonYield[cbin]++;
                sumPtKaon[cbin] += pt;
                countPtKaon[cbin]++;
                hPtMBKaon->Fill(pt);
            }

            if (eposID == 1120 && ist[i] == 0 && (fabs(y) < 0.5))
            {
                protonYield[cbin]++;
                sumPtProton[cbin] += pt;
                countPtProton[cbin]++;
                hPtMBProton->Fill(pt);
            }

            if (ist[i] != 7) // ist 7 = UrQMD off, ist = 9 UrQMD on
                continue;

            // cout << "EPOS ID: " << eposID << ", IST: " << ist[i] << ", Rapidity: " << y << endl;

            if (atLeastOnePiKp_in_ModEta1 == 0) // INEL > 0 selection
                continue;

            // ----------------------
            // PHI
            // ----------------------
            if (TMath::Abs(id[i]) == 331 && (fabs(y) < 0.5))

            {
                phiYield[cbin]++;
                sumPtPhi[cbin] += pt;
                countPtPhi[cbin]++;
                hPtMBPhi->Fill(pt);
            }

            // ----------------------
            // KSTAR
            // ----------------------
            if (TMath::Abs(id[i]) == 231 && (fabs(y) < 0.5))
            {
                kstarYield[cbin]++;
                sumPtKstar[cbin] += pt;
                countPtKstar[cbin]++;
                hPtMBKstar->Fill(pt);
            }

            // ----------------------
            // Xi(1530)
            // ----------------------
            if (TMath::Abs(id[i]) == 1331 && (fabs(y) < 0.5))
            {
                XiStarYield[cbin]++;
                sumPtXiStar[cbin] += pt;
                countPtXiStar[cbin]++;
                hPtMBXiStar->Fill(pt);
            }

            // ----------------------
            // Charged Kstar
            // ----------------------
            if (TMath::Abs(id[i]) == 131 && (fabs(y) < 0.5))
            {
                ChargedKstarYield[cbin]++;
                sumPtChargedKstar[cbin] += pt;
                countPtChargedKstar[cbin]++;
                hPtMBChargedKstar->Fill(pt);
            }
        }

        sumDNdEta[cbin] += nchMid;
    }

    // =====================================================
    // branching ratios
    // =====================================================
    const double BR_phi = 0.492;
    const double BR_kstar = 0.666;
    const double BR_kshort = 1.0;
    const double BR_proton = 1.0;
    const double BR_pion = 1.0;
    const double BR_kaon = 1.0;
    const double BR_XiStar = 0.666;
    const double BR_ChargedKstar = 0.333;

    double meanDNdEta[10];
    double meanPtPhi[10] = {0};
    double meanPtKstar[10] = {0};
    double meanPtKshort[10] = {0};
    double meanPtProton[10] = {0};
    double meanPtPion[10] = {0};
    double meanPtKaon[10] = {0};
    double meanPtXiStar[10] = {0};
    double meanPtChargedKstar[10] = {0};

    for (int i = 0; i < 10; i++)
    {
        if (nEventsCent[i] > 0)
        {
            meanDNdEta[i] = sumDNdEta[i] / nEventsCent[i];

            phiYield[i] = phiYield[i] / (nEventsCent[i] * BR_phi);
            kstarYield[i] = kstarYield[i] / (nEventsCent[i] * BR_kstar * 2); // Average of K*0 and K*0bar
            kshortYield[i] = kshortYield[i] / (nEventsCent[i] * BR_kshort);
            protonYield[i] = protonYield[i] / (nEventsCent[i] * BR_proton);
            pionYield[i] = pionYield[i] / (nEventsCent[i] * BR_pion);
            kaonYield[i] = kaonYield[i] / (nEventsCent[i] * BR_kaon);
            XiStarYield[i] = XiStarYield[i] / (nEventsCent[i] * BR_XiStar);
            ChargedKstarYield[i] = ChargedKstarYield[i] / (nEventsCent[i] * BR_ChargedKstar);

            if (countPtPhi[i] > 0)
                meanPtPhi[i] = sumPtPhi[i] / countPtPhi[i];
            if (countPtKstar[i] > 0)
                meanPtKstar[i] = sumPtKstar[i] / countPtKstar[i];
            if (countPtKshort[i] > 0)
                meanPtKshort[i] = sumPtKshort[i] / countPtKshort[i];
            if (countPtProton[i] > 0)
                meanPtProton[i] = sumPtProton[i] / countPtProton[i];
            if (countPtPion[i] > 0)
                meanPtPion[i] = sumPtPion[i] / countPtPion[i];
            if (countPtKaon[i] > 0)
                meanPtKaon[i] = sumPtKaon[i] / countPtKaon[i];
            if (countPtXiStar[i] > 0)
                meanPtXiStar[i] = sumPtXiStar[i] / countPtXiStar[i];
            if (countPtChargedKstar[i] > 0)
                meanPtChargedKstar[i] = sumPtChargedKstar[i] / countPtChargedKstar[i];
        }
    }

    // =====================================================
    // print
    // =====================================================
    cout << "\n===== URQMD OFF RESULTS =====\n";
    for (int i = 0; i < 10; i++)
    {
        cout << "bin " << i
             << "  <dNch/deta>=" << meanDNdEta[i]
             << "  phi=" << phiYield[i]
             << "  kstar=" << kstarYield[i]
             << endl;
    }

    // =====================================================
    // graphs
    // =====================================================
    TGraph *gPhi = new TGraph(10, meanDNdEta, phiYield);
    TGraph *gKstar = new TGraph(10, meanDNdEta, kstarYield);
    TGraph *gKshort = new TGraph(10, meanDNdEta, kshortYield);
    TGraph *gProton = new TGraph(10, meanDNdEta, protonYield);
    TGraph *gPion = new TGraph(10, meanDNdEta, pionYield);
    TGraph *gKaon = new TGraph(10, meanDNdEta, kaonYield);
    TGraph *gXiStar = new TGraph(10, meanDNdEta, XiStarYield);
    TGraph *gChargedKstar = new TGraph(10, meanDNdEta, ChargedKstarYield);

    TGraph *gMeanPtPhi = new TGraph(10, meanDNdEta, meanPtPhi);
    TGraph *gMeanPtKstar = new TGraph(10, meanDNdEta, meanPtKstar);
    TGraph *gMeanPtKshort = new TGraph(10, meanDNdEta, meanPtKshort);
    TGraph *gMeanPtProton = new TGraph(10, meanDNdEta, meanPtProton);
    TGraph *gMeanPtPion = new TGraph(10, meanDNdEta, meanPtPion);
    TGraph *gMeanPtKaon = new TGraph(10, meanDNdEta, meanPtKaon);
    TGraph *gMeanPtXiStar = new TGraph(10, meanDNdEta, meanPtXiStar);
    TGraph *gMeanPtChargedKstar = new TGraph(10, meanDNdEta, meanPtChargedKstar);

    gPhi->SetName("phi_vs_mult_urqmdOFF");
    gKstar->SetName("kstar_vs_mult_urqmdOFF");
    gKshort->SetName("kshort_vs_mult_urqmdOFF");
    gProton->SetName("proton_vs_mult_urqmdOFF");
    gPion->SetName("pion_vs_mult_urqmdOFF");
    gKaon->SetName("kaon_vs_mult_urqmdOFF");
    gXiStar->SetName("xistar_vs_mult_urqmdOFF");
    gChargedKstar->SetName("chargedkstar_vs_mult_urqmdOFF");

    gMeanPtPhi->SetName("meanpt_phi_vs_mult_urqmdOFF");
    gMeanPtKstar->SetName("meanpt_kstar_vs_mult_urqmdOFF");
    gMeanPtKshort->SetName("meanpt_kshort_vs_mult_urqmdOFF");
    gMeanPtProton->SetName("meanpt_proton_vs_mult_urqmdOFF");
    gMeanPtPion->SetName("meanpt_pion_vs_mult_urqmdOFF");
    gMeanPtKaon->SetName("meanpt_kaon_vs_mult_urqmdOFF");
    gMeanPtXiStar->SetName("meanpt_xistar_vs_mult_urqmdOFF");
    gMeanPtChargedKstar->SetName("meanpt_chargedkstar_vs_mult_urqmdOFF");

    // =====================================================
    // normalize minimum-bias pT spectra
    // =====================================================
    const double totalEvents = static_cast<double>(nEv);
    hPtMBPhi->Scale(1.0 / (totalEvents * BR_phi * hPtMBPhi->GetBinWidth(1)));
    hPtMBKstar->Scale(1.0 / (totalEvents * BR_kstar * 2 * hPtMBKstar->GetBinWidth(1)));
    hPtMBKshort->Scale(1.0 / (totalEvents * BR_kshort * hPtMBKshort->GetBinWidth(1)));
    hPtMBProton->Scale(1.0 / (totalEvents * BR_proton * hPtMBProton->GetBinWidth(1)));
    hPtMBPion->Scale(1.0 / (totalEvents * BR_pion * hPtMBPion->GetBinWidth(1)));
    hPtMBKaon->Scale(1.0 / (totalEvents * BR_kaon * hPtMBKaon->GetBinWidth(1)));
    hPtMBXiStar->Scale(1.0 / (totalEvents * BR_XiStar * hPtMBXiStar->GetBinWidth(1)));
    hPtMBChargedKstar->Scale(1.0 / (totalEvents * BR_ChargedKstar * hPtMBChargedKstar->GetBinWidth(1)));

    // =====================================================
    // save
    // =====================================================
    TFile fout("EPOS_UrQMD_Off.root", "RECREATE");
    hFT0->Write();
    hIST->Write();
    hITY->Write();
    gPhi->Write();
    gKstar->Write();
    gKshort->Write();
    gProton->Write();
    gPion->Write();
    gKaon->Write();
    gXiStar->Write();
    gChargedKstar->Write();
    gMeanPtPhi->Write();
    gMeanPtKstar->Write();
    gMeanPtKshort->Write();
    gMeanPtProton->Write();
    gMeanPtPion->Write();
    gMeanPtKaon->Write();
    gMeanPtXiStar->Write();
    gMeanPtChargedKstar->Write();
    hPtMBPhi->Write();
    hPtMBKstar->Write();
    hPtMBKshort->Write();
    hPtMBProton->Write();
    hPtMBPion->Write();
    hPtMBKaon->Write();
    hPtMBXiStar->Write();
    hPtMBChargedKstar->Write();
    fout.Close();
    cout << "\nSaved EPOS_UrQMD_Off.root\n";
}
