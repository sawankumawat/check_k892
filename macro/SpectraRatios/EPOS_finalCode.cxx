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
void EPOS_finalCode()
{
    TDatabasePDG *pdg = new TDatabasePDG();

    // =====================================================
    // READ ALL FILES
    // =====================================================
    TChain chain("teposevent0");
    // chain.Add("/home/sawan/Storage/EPOS_localOutputs/mergedEPOS_UrQMD.root");
    chain.Add("/home/sawan/Storage/EPOS_localOutputs/merged.root");

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
    TH2F *hRapvsIST = new TH2F("hRapvsIST", "Rapidity vs IST;IST;y", 20, -0.5, 19.5, 100, -1, 1);
    TH2F *hIdvsIST = new TH2F("hIdvsIST", "EPOS ID vs IST;IST;EPOS ID", 20, -0.5, 19.5, 3000, -1500, 1500);

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

    double phiYieldUrQMDOff[10] = {0};
    double kstarYieldUrQMDOff[10] = {0};
    double XiStarYieldUrQMDOff[10] = {0};
    double ChargedKstarYieldUrQMDOff[10] = {0};

    double phiYieldNoRescattering[10] = {0};
    double kstarYieldNoRescattering[10] = {0};
    double XiStarYieldNoRescattering[10] = {0};
    double ChargedKstarYieldNoRescattering[10] = {0};

    double phiYieldRescattering[10] = {0};
    double kstarYieldRescattering[10] = {0};
    double XiStarYieldRescattering[10] = {0};
    double ChargedKstarYieldRescattering[10] = {0};

    double sumPtPhi[10] = {0};
    double sumPtKstar[10] = {0};
    double sumPtKshort[10] = {0};
    double sumPtProton[10] = {0};
    double sumPtPion[10] = {0};
    double sumPtKaon[10] = {0};
    double sumPtXiStar[10] = {0};
    double sumPtChargedKstar[10] = {0};

    double sumPtPhiUrQMDOff[10] = {0};
    double sumPtKstarUrQMDOff[10] = {0};
    double sumPtXiStarUrQMDOff[10] = {0};
    double sumPtChargedKstarUrQMDOff[10] = {0};

    double sumPtPhiNoRescattering[10] = {0};
    double sumPtKstarNoRescattering[10] = {0};
    double sumPtXiStarNoRescattering[10] = {0};
    double sumPtChargedKstarNoRescattering[10] = {0};

    double sumPtPhiRescattering[10] = {0};
    double sumPtKstarRescattering[10] = {0};
    double sumPtXiStarRescattering[10] = {0};
    double sumPtChargedKstarRescattering[10] = {0};

    double countPtPhi[10] = {0};
    double countPtKstar[10] = {0};
    double countPtKshort[10] = {0};
    double countPtProton[10] = {0};
    double countPtPion[10] = {0};
    double countPtKaon[10] = {0};
    double countPtXiStar[10] = {0};
    double countPtChargedKstar[10] = {0};

    double countPtPhiUrQMDOff[10] = {0};
    double countPtKstarUrQMDOff[10] = {0};
    double countPtXiStarUrQMDOff[10] = {0};
    double countPtChargedKstarUrQMDOff[10] = {0};

    double countPtPhiNoRescattering[10] = {0};
    double countPtKstarNoRescattering[10] = {0};
    double countPtXiStarNoRescattering[10] = {0};
    double countPtChargedKstarNoRescattering[10] = {0};

    double countPtPhiRescattering[10] = {0};
    double countPtKstarRescattering[10] = {0};
    double countPtXiStarRescattering[10] = {0};
    double countPtChargedKstarRescattering[10] = {0};

    TH1F *hPtMBPhi = new TH1F("hPtMBPhi", "#phi p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBKstar = new TH1F("hPtMBKstar", "K*^{0} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBKshort = new TH1F("hPtMBKshort", "K^{0}_{S} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBProton = new TH1F("hPtMBProton", "p/#bar{p} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBPion = new TH1F("hPtMBPion", "#pi^{#pm} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBKaon = new TH1F("hPtMBKaon", "K^{#pm} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBXiStar = new TH1F("hPtMBXiStar", "#Xi(1530) p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBChargedKstar = new TH1F("hPtMBChargedKstar", "K*^{#pm} p_{T} (minimum bias);p_{T} (GeV/c);Counts", 200, 0, 20);

    TH1F *hPtMBPhiUrQMDOff = new TH1F("hPtMBPhiUrQMDOff", "#phi p_{T} UrQMD off;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBKstarUrQMDOff = new TH1F("hPtMBKstarUrQMDOff", "K*^{0} p_{T} UrQMD off;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBXiStarUrQMDOff = new TH1F("hPtMBXiStarUrQMDOff", "#Xi(1530) p_{T} UrQMD off;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBChargedKstarUrQMDOff = new TH1F("hPtMBChargedKstarUrQMDOff", "K*^{#pm} p_{T} UrQMD off;p_{T} (GeV/c);Counts", 200, 0, 20);

    TH1F *hPtMBPhiNoRescattering = new TH1F("hPtMBPhiNoRescattering", "#phi p_{T} no-rescattering;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBKstarNoRescattering = new TH1F("hPtMBKstarNoRescattering", "K*^{0} p_{T} no-rescattering;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBXiStarNoRescattering = new TH1F("hPtMBXiStarNoRescattering", "#Xi(1530) p_{T} no-rescattering;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBChargedKstarNoRescattering = new TH1F("hPtMBChargedKstarNoRescattering", "K*^{#pm} p_{T} no-rescattering;p_{T} (GeV/c);Counts", 200, 0, 20);

    TH1F *hPtMBPhiRescattering = new TH1F("hPtMBPhiRescattering", "#phi p_{T} rescattering;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBKstarRescattering = new TH1F("hPtMBKstarRescattering", "K*^{0} p_{T} rescattering;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBXiStarRescattering = new TH1F("hPtMBXiStarRescattering", "#Xi(1530) p_{T} rescattering;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F *hPtMBChargedKstarRescattering = new TH1F("hPtMBChargedKstarRescattering", "K*^{#pm} p_{T} rescattering;p_{T} (GeV/c);Counts", 200, 0, 20);

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
            // URQMD ON → final hadrons
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

            hRapvsIST->Fill(ist[i], y);
            hIdvsIST->Fill(ist[i], eposID);
            if (atLeastOnePiKp_in_ModEta1 == 0) // INEL > 0 condition
                continue;

            if (ist[i] == 8) // ist 7 = UrQMD ON, ist = 9 UrQMD on
            {
                // ----------------------
                // PHI
                // ----------------------
                if (eposID == 331 || eposID == 33100)
                    if (fabs(y) < 0.5)
                    {
                        phiYieldUrQMDOff[cbin]++;
                        sumPtPhiUrQMDOff[cbin] += pt;
                        countPtPhiUrQMDOff[cbin]++;
                        hPtMBPhiUrQMDOff->Fill(pt);
                    }

                // ----------------------
                // KSTAR
                // ----------------------
                if (eposID == 231 || eposID == 23100)
                    if (fabs(y) < 0.5)
                    {
                        kstarYieldUrQMDOff[cbin]++;
                        sumPtKstarUrQMDOff[cbin] += pt;
                        countPtKstarUrQMDOff[cbin]++;
                        hPtMBKstarUrQMDOff->Fill(pt);
                    }

                // ----------------------
                // Xi(1530)
                // ----------------------
                if (eposID == 1331 || eposID == 133100)
                    if (fabs(y) < 0.5)
                    {
                        XiStarYieldUrQMDOff[cbin]++;
                        sumPtXiStarUrQMDOff[cbin] += pt;
                        countPtXiStarUrQMDOff[cbin]++;
                        hPtMBXiStarUrQMDOff->Fill(pt);
                    }

                // ----------------------
                // Charged Kstar
                // ----------------------
                if (eposID == 131 || eposID == 13100)
                    if (fabs(y) < 0.5)
                    {
                        ChargedKstarYieldUrQMDOff[cbin]++;
                        sumPtChargedKstarUrQMDOff[cbin] += pt;
                        countPtChargedKstarUrQMDOff[cbin]++;
                        hPtMBChargedKstarUrQMDOff->Fill(pt);
                    }
            }

            // if (ev < 500) // for debugging
            // {
            //     // if (ist[i] == 9 && ity[i] == 80)
            //     //     cout << "Event " << ev << ": No-rescattering, EPOS ID = " << id[i] << ", pt = " << pt << ", y = " << y << endl;
            //     // if (ist[i] == 9 && ity[i] == 81)
            //     //     cout << "Event " << ev << ": rescattering, EPOS ID = " << id[i] << ", pt = " << pt << ", y = " << y << endl;

            //     // if (ist[i] == 9)
            //     // cout << "Event " << ev << ": UrQMD ON, EPOS ID = " << id[i] << ", pt = " << pt << ", y = " << y << endl;

            //     if (ist[i] == 8)
            //     cout << "Event " << ev << ": UrQMD OFF, EPOS ID = " << id[i] << ", pt = " << pt << ", y = " << y << endl;
            // }

            if (ist[i] == 9) // ist 7 or 8 = UrQMD Off, ist = 9 UrQMD ON
            {
                // ----------------------
                // PHI
                // ----------------------
                if (eposID == 331 || eposID == 33100 || eposID == 33101)
                    if (fabs(y) < 0.5)
                    {
                        phiYield[cbin]++;
                        sumPtPhi[cbin] += pt;
                        countPtPhi[cbin]++;
                        hPtMBPhi->Fill(pt);
                    }

                // ----------------------
                // KSTAR
                // ----------------------
                if (eposID == 231 || eposID == 23100 || eposID == 23101)
                    if (fabs(y) < 0.5)
                    {
                        kstarYield[cbin]++;
                        sumPtKstar[cbin] += pt;
                        countPtKstar[cbin]++;
                        hPtMBKstar->Fill(pt);
                    }

                // ----------------------
                // Xi(1530)
                // ----------------------
                if (eposID == 1331 || eposID == 133100 || eposID == 133101)
                    if (fabs(y) < 0.5)
                    {
                        XiStarYield[cbin]++;
                        sumPtXiStar[cbin] += pt;
                        countPtXiStar[cbin]++;
                        hPtMBXiStar->Fill(pt);
                    }

                // ----------------------
                // Charged Kstar
                // ----------------------
                if (eposID == 131 || eposID == 13100 || eposID == 13101)
                    if (fabs(y) < 0.5)
                    {
                        ChargedKstarYield[cbin]++;
                        sumPtChargedKstar[cbin] += pt;
                        countPtChargedKstar[cbin]++;
                        hPtMBChargedKstar->Fill(pt);
                    }
            }

            if (ist[i] == 9 && ity[i] == 80) // no-rescattering
            {
                if (eposID == 331 || eposID == 33100 || eposID == 33101)
                    if (fabs(y) < 0.5)
                    {
                        phiYieldNoRescattering[cbin]++;
                        sumPtPhiNoRescattering[cbin] += pt;
                        countPtPhiNoRescattering[cbin]++;
                        hPtMBPhiNoRescattering->Fill(pt);
                    }

                if (eposID == 231 || eposID == 23100 || eposID == 23101)
                    if (fabs(y) < 0.5)
                    {
                        kstarYieldNoRescattering[cbin]++;
                        sumPtKstarNoRescattering[cbin] += pt;
                        countPtKstarNoRescattering[cbin]++;
                        hPtMBKstarNoRescattering->Fill(pt);
                    }

                if (eposID == 1331 || eposID == 133100 || eposID == 133101)
                    if (fabs(y) < 0.5)
                    {
                        XiStarYieldNoRescattering[cbin]++;
                        sumPtXiStarNoRescattering[cbin] += pt;
                        countPtXiStarNoRescattering[cbin]++;
                        hPtMBXiStarNoRescattering->Fill(pt);
                    }

                if (eposID == 131 || eposID == 13100 || eposID == 13101)
                    if (fabs(y) < 0.5)
                    {
                        ChargedKstarYieldNoRescattering[cbin]++;
                        sumPtChargedKstarNoRescattering[cbin] += pt;
                        countPtChargedKstarNoRescattering[cbin]++;
                        hPtMBChargedKstarNoRescattering->Fill(pt);
                    }
            }

            if (ist[i] == 9 && ity[i] == 81) // rescattering
            {
                if (eposID == 331 || eposID == 33100 || eposID == 33101)
                    if (fabs(y) < 0.5)
                    {
                        phiYieldRescattering[cbin]++;
                        sumPtPhiRescattering[cbin] += pt;
                        countPtPhiRescattering[cbin]++;
                        hPtMBPhiRescattering->Fill(pt);
                    }

                if (eposID == 231 || eposID == 23100 || eposID == 23101)
                    if (fabs(y) < 0.5)
                    {
                        kstarYieldRescattering[cbin]++;
                        sumPtKstarRescattering[cbin] += pt;
                        countPtKstarRescattering[cbin]++;
                        hPtMBKstarRescattering->Fill(pt);
                    }

                if (eposID == 1331 || eposID == 133100 || eposID == 133101)
                    if (fabs(y) < 0.5)
                    {
                        XiStarYieldRescattering[cbin]++;
                        sumPtXiStarRescattering[cbin] += pt;
                        countPtXiStarRescattering[cbin]++;
                        hPtMBXiStarRescattering->Fill(pt);
                    }

                if (eposID == 131 || eposID == 13100 || eposID == 13101)
                    if (fabs(y) < 0.5)
                    {
                        ChargedKstarYieldRescattering[cbin]++;
                        sumPtChargedKstarRescattering[cbin] += pt;
                        countPtChargedKstarRescattering[cbin]++;
                        hPtMBChargedKstarRescattering->Fill(pt);
                    }
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

    double meanDNdEta[10] = {0};
    double meanPtPhi[10] = {0};
    double meanPtKstar[10] = {0};
    double meanPtKshort[10] = {0};
    double meanPtProton[10] = {0};
    double meanPtPion[10] = {0};
    double meanPtKaon[10] = {0};
    double meanPtXiStar[10] = {0};
    double meanPtChargedKstar[10] = {0};

    double meanPtPhiUrQMDOff[10] = {0};
    double meanPtKstarUrQMDOff[10] = {0};
    double meanPtXiStarUrQMDOff[10] = {0};
    double meanPtChargedKstarUrQMDOff[10] = {0};

    double meanPtPhiNoRescattering[10] = {0};
    double meanPtKstarNoRescattering[10] = {0};
    double meanPtXiStarNoRescattering[10] = {0};
    double meanPtChargedKstarNoRescattering[10] = {0};

    double meanPtPhiRescattering[10] = {0};
    double meanPtKstarRescattering[10] = {0};
    double meanPtXiStarRescattering[10] = {0};
    double meanPtChargedKstarRescattering[10] = {0};

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

            phiYieldUrQMDOff[i] = phiYieldUrQMDOff[i] / (nEventsCent[i] * BR_phi);
            kstarYieldUrQMDOff[i] = kstarYieldUrQMDOff[i] / (nEventsCent[i] * BR_kstar * 2);
            XiStarYieldUrQMDOff[i] = XiStarYieldUrQMDOff[i] / (nEventsCent[i] * BR_XiStar);
            ChargedKstarYieldUrQMDOff[i] = ChargedKstarYieldUrQMDOff[i] / (nEventsCent[i] * BR_ChargedKstar);

            phiYieldNoRescattering[i] = phiYieldNoRescattering[i] / (nEventsCent[i] * BR_phi);
            kstarYieldNoRescattering[i] = kstarYieldNoRescattering[i] / (nEventsCent[i] * BR_kstar * 2);
            XiStarYieldNoRescattering[i] = XiStarYieldNoRescattering[i] / (nEventsCent[i] * BR_XiStar);
            ChargedKstarYieldNoRescattering[i] = ChargedKstarYieldNoRescattering[i] / (nEventsCent[i] * BR_ChargedKstar);

            phiYieldRescattering[i] = phiYieldRescattering[i] / (nEventsCent[i] * BR_phi);
            kstarYieldRescattering[i] = kstarYieldRescattering[i] / (nEventsCent[i] * BR_kstar * 2);
            XiStarYieldRescattering[i] = XiStarYieldRescattering[i] / (nEventsCent[i] * BR_XiStar);
            ChargedKstarYieldRescattering[i] = ChargedKstarYieldRescattering[i] / (nEventsCent[i] * BR_ChargedKstar);

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

            if (countPtPhiUrQMDOff[i] > 0)
                meanPtPhiUrQMDOff[i] = sumPtPhiUrQMDOff[i] / countPtPhiUrQMDOff[i];
            if (countPtKstarUrQMDOff[i] > 0)
                meanPtKstarUrQMDOff[i] = sumPtKstarUrQMDOff[i] / countPtKstarUrQMDOff[i];
            if (countPtXiStarUrQMDOff[i] > 0)
                meanPtXiStarUrQMDOff[i] = sumPtXiStarUrQMDOff[i] / countPtXiStarUrQMDOff[i];
            if (countPtChargedKstarUrQMDOff[i] > 0)
                meanPtChargedKstarUrQMDOff[i] = sumPtChargedKstarUrQMDOff[i] / countPtChargedKstarUrQMDOff[i];

            if (countPtPhiNoRescattering[i] > 0)
                meanPtPhiNoRescattering[i] = sumPtPhiNoRescattering[i] / countPtPhiNoRescattering[i];
            if (countPtKstarNoRescattering[i] > 0)
                meanPtKstarNoRescattering[i] = sumPtKstarNoRescattering[i] / countPtKstarNoRescattering[i];
            if (countPtXiStarNoRescattering[i] > 0)
                meanPtXiStarNoRescattering[i] = sumPtXiStarNoRescattering[i] / countPtXiStarNoRescattering[i];
            if (countPtChargedKstarNoRescattering[i] > 0)
                meanPtChargedKstarNoRescattering[i] = sumPtChargedKstarNoRescattering[i] / countPtChargedKstarNoRescattering[i];

            if (countPtPhiRescattering[i] > 0)
                meanPtPhiRescattering[i] = sumPtPhiRescattering[i] / countPtPhiRescattering[i];
            if (countPtKstarRescattering[i] > 0)
                meanPtKstarRescattering[i] = sumPtKstarRescattering[i] / countPtKstarRescattering[i];
            if (countPtXiStarRescattering[i] > 0)
                meanPtXiStarRescattering[i] = sumPtXiStarRescattering[i] / countPtXiStarRescattering[i];
            if (countPtChargedKstarRescattering[i] > 0)
                meanPtChargedKstarRescattering[i] = sumPtChargedKstarRescattering[i] / countPtChargedKstarRescattering[i];
        }
    }

    // =====================================================
    // print
    // =====================================================
    cout << "\n===== URQMD ON RESULTS =====\n";
    for (int i = 0; i < 10; i++)
    {
        cout << "bin " << i
             << "  <dNch/deta>=" << meanDNdEta[i]
             << "  phi=" << phiYield[i]
             << "  kstar=" << kstarYield[i]
             << endl;
    }

    cout << "\n===== URQMD OFF RESULTS =====\n";
    for (int i = 0; i < 10; i++)
    {
        cout << "bin " << i
             << "  <dNch/deta>=" << meanDNdEta[i]
             << "  phi=" << phiYieldUrQMDOff[i]
             << "  kstar=" << kstarYieldUrQMDOff[i]
             << endl;
    }

    cout << "\n===== UrQMD ON: No Rescattering RESULTS =====\n";
    for (int i = 0; i < 10; i++)
    {
        cout << "bin " << i
             << "  <dNch/deta>=" << meanDNdEta[i]
             << "  phi=" << phiYieldNoRescattering[i]
             << "  kstar=" << kstarYieldNoRescattering[i]
             << endl;
    }

    cout << "\n===== UrQMD ON: Rescattering RESULTS =====\n";
    for (int i = 0; i < 10; i++)
    {
        cout << "bin " << i
             << "  <dNch/deta>=" << meanDNdEta[i]
             << "  phi=" << phiYieldRescattering[i]
             << "  kstar=" << kstarYieldRescattering[i]
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

    TGraph *gPhiUrQMDOff = new TGraph(10, meanDNdEta, phiYieldUrQMDOff);
    TGraph *gKstarUrQMDOff = new TGraph(10, meanDNdEta, kstarYieldUrQMDOff);
    TGraph *gXiStarUrQMDOff = new TGraph(10, meanDNdEta, XiStarYieldUrQMDOff);
    TGraph *gChargedKstarUrQMDOff = new TGraph(10, meanDNdEta, ChargedKstarYieldUrQMDOff);

    TGraph *gPhiNoRescattering = new TGraph(10, meanDNdEta, phiYieldNoRescattering);
    TGraph *gKstarNoRescattering = new TGraph(10, meanDNdEta, kstarYieldNoRescattering);
    TGraph *gXiStarNoRescattering = new TGraph(10, meanDNdEta, XiStarYieldNoRescattering);
    TGraph *gChargedKstarNoRescattering = new TGraph(10, meanDNdEta, ChargedKstarYieldNoRescattering);

    TGraph *gPhiRescattering = new TGraph(10, meanDNdEta, phiYieldRescattering);
    TGraph *gKstarRescattering = new TGraph(10, meanDNdEta, kstarYieldRescattering);
    TGraph *gXiStarRescattering = new TGraph(10, meanDNdEta, XiStarYieldRescattering);
    TGraph *gChargedKstarRescattering = new TGraph(10, meanDNdEta, ChargedKstarYieldRescattering);

    TGraph *gMeanPtPhiUrQMDOff = new TGraph(10, meanDNdEta, meanPtPhiUrQMDOff);
    TGraph *gMeanPtKstarUrQMDOff = new TGraph(10, meanDNdEta, meanPtKstarUrQMDOff);
    TGraph *gMeanPtXiStarUrQMDOff = new TGraph(10, meanDNdEta, meanPtXiStarUrQMDOff);
    TGraph *gMeanPtChargedKstarUrQMDOff = new TGraph(10, meanDNdEta, meanPtChargedKstarUrQMDOff);

    TGraph *gMeanPtPhiNoRescattering = new TGraph(10, meanDNdEta, meanPtPhiNoRescattering);
    TGraph *gMeanPtKstarNoRescattering = new TGraph(10, meanDNdEta, meanPtKstarNoRescattering);
    TGraph *gMeanPtXiStarNoRescattering = new TGraph(10, meanDNdEta, meanPtXiStarNoRescattering);
    TGraph *gMeanPtChargedKstarNoRescattering = new TGraph(10, meanDNdEta, meanPtChargedKstarNoRescattering);

    TGraph *gMeanPtPhiRescattering = new TGraph(10, meanDNdEta, meanPtPhiRescattering);
    TGraph *gMeanPtKstarRescattering = new TGraph(10, meanDNdEta, meanPtKstarRescattering);
    TGraph *gMeanPtXiStarRescattering = new TGraph(10, meanDNdEta, meanPtXiStarRescattering);
    TGraph *gMeanPtChargedKstarRescattering = new TGraph(10, meanDNdEta, meanPtChargedKstarRescattering);

    gPhi->SetName("phi_vs_mult_urqmdON");
    gKstar->SetName("kstar_vs_mult_urqmdON");
    gKshort->SetName("kshort_vs_mult_urqmdON");
    gProton->SetName("proton_vs_mult_urqmdON");
    gPion->SetName("pion_vs_mult_urqmdON");
    gKaon->SetName("kaon_vs_mult_urqmdON");
    gXiStar->SetName("xistar_vs_mult_urqmdON");
    gChargedKstar->SetName("chargedkstar_vs_mult_urqmdON");

    gPhiUrQMDOff->SetName("phi_vs_mult_urqmdOFF");
    gKstarUrQMDOff->SetName("kstar_vs_mult_urqmdOFF");
    gXiStarUrQMDOff->SetName("xistar_vs_mult_urqmdOFF");
    gChargedKstarUrQMDOff->SetName("chargedkstar_vs_mult_urqmdOFF");

    gMeanPtPhi->SetName("meanpt_phi_vs_mult_urqmdON");
    gMeanPtKstar->SetName("meanpt_kstar_vs_mult_urqmdON");
    gMeanPtKshort->SetName("meanpt_kshort_vs_mult_urqmdON");
    gMeanPtProton->SetName("meanpt_proton_vs_mult_urqmdON");
    gMeanPtPion->SetName("meanpt_pion_vs_mult_urqmdON");
    gMeanPtKaon->SetName("meanpt_kaon_vs_mult_urqmdON");
    gMeanPtXiStar->SetName("meanpt_xistar_vs_mult_urqmdON");
    gMeanPtChargedKstar->SetName("meanpt_chargedkstar_vs_mult_urqmdON");

    gMeanPtPhiUrQMDOff->SetName("meanpt_phi_vs_mult_urqmdOFF");
    gMeanPtKstarUrQMDOff->SetName("meanpt_kstar_vs_mult_urqmdOFF");
    gMeanPtXiStarUrQMDOff->SetName("meanpt_xistar_vs_mult_urqmdOFF");
    gMeanPtChargedKstarUrQMDOff->SetName("meanpt_chargedkstar_vs_mult_urqmdOFF");

    gPhiNoRescattering->SetName("phi_vs_mult_urqmdON_noRescattering");
    gKstarNoRescattering->SetName("kstar_vs_mult_urqmdON_noRescattering");
    gXiStarNoRescattering->SetName("xistar_vs_mult_urqmdON_noRescattering");
    gChargedKstarNoRescattering->SetName("chargedkstar_vs_mult_urqmdON_noRescattering");

    gPhiRescattering->SetName("phi_vs_mult_urqmdON_rescattering");
    gKstarRescattering->SetName("kstar_vs_mult_urqmdON_rescattering");
    gXiStarRescattering->SetName("xistar_vs_mult_urqmdON_rescattering");
    gChargedKstarRescattering->SetName("chargedkstar_vs_mult_urqmdON_rescattering");

    gMeanPtPhiNoRescattering->SetName("meanpt_phi_vs_mult_urqmdON_noRescattering");
    gMeanPtKstarNoRescattering->SetName("meanpt_kstar_vs_mult_urqmdON_noRescattering");
    gMeanPtXiStarNoRescattering->SetName("meanpt_xistar_vs_mult_urqmdON_noRescattering");
    gMeanPtChargedKstarNoRescattering->SetName("meanpt_chargedkstar_vs_mult_urqmdON_noRescattering");

    gMeanPtPhiRescattering->SetName("meanpt_phi_vs_mult_urqmdON_rescattering");
    gMeanPtKstarRescattering->SetName("meanpt_kstar_vs_mult_urqmdON_rescattering");
    gMeanPtXiStarRescattering->SetName("meanpt_xistar_vs_mult_urqmdON_rescattering");
    gMeanPtChargedKstarRescattering->SetName("meanpt_chargedkstar_vs_mult_urqmdON_rescattering");

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

    hPtMBPhiNoRescattering->Scale(1.0 / (totalEvents * BR_phi * hPtMBPhiNoRescattering->GetBinWidth(1)));
    hPtMBKstarNoRescattering->Scale(1.0 / (totalEvents * BR_kstar * 2 * hPtMBKstarNoRescattering->GetBinWidth(1)));
    hPtMBXiStarNoRescattering->Scale(1.0 / (totalEvents * BR_XiStar * hPtMBXiStarNoRescattering->GetBinWidth(1)));
    hPtMBChargedKstarNoRescattering->Scale(1.0 / (totalEvents * BR_ChargedKstar * hPtMBChargedKstarNoRescattering->GetBinWidth(1)));

    hPtMBPhiRescattering->Scale(1.0 / (totalEvents * BR_phi * hPtMBPhiRescattering->GetBinWidth(1)));
    hPtMBKstarRescattering->Scale(1.0 / (totalEvents * BR_kstar * 2 * hPtMBKstarRescattering->GetBinWidth(1)));
    hPtMBXiStarRescattering->Scale(1.0 / (totalEvents * BR_XiStar * hPtMBXiStarRescattering->GetBinWidth(1)));
    hPtMBChargedKstarRescattering->Scale(1.0 / (totalEvents * BR_ChargedKstar * hPtMBChargedKstarRescattering->GetBinWidth(1)));

    // =====================================================
    // save
    // =====================================================
    TFile *fout = new TFile("EPOS_final.root", "RECREATE");
    hFT0->Write();
    hIST->Write();
    hITY->Write();
    hRapvsIST->Write();
    hIdvsIST->Write();
    TDirectory *dirUrQMDOn = (TDirectory *)fout->mkdir("UrQMD_ON");
    dirUrQMDOn->cd();
    // Mean Yield
    gPhi->Write();
    gKstar->Write();
    gKshort->Write();
    gProton->Write();
    gPion->Write();
    gKaon->Write();
    gXiStar->Write();
    gChargedKstar->Write();
    // Mean pT
    gMeanPtPhi->Write();
    gMeanPtKstar->Write();
    gMeanPtKshort->Write();
    gMeanPtProton->Write();
    gMeanPtPion->Write();
    gMeanPtKaon->Write();
    gMeanPtXiStar->Write();
    gMeanPtChargedKstar->Write();
    // Minimum-bias pT spectra
    hPtMBPhi->Write();
    hPtMBKstar->Write();
    hPtMBKshort->Write();
    hPtMBProton->Write();
    hPtMBPion->Write();
    hPtMBKaon->Write();
    hPtMBXiStar->Write();
    hPtMBChargedKstar->Write();

    fout->cd();
    TDirectory *dirUrQMDOff = (TDirectory *)fout->mkdir("UrQMD_OFF");
    dirUrQMDOff->cd();
    // Mean Yield
    gKstarUrQMDOff->Write();
    gPhiUrQMDOff->Write();
    gXiStarUrQMDOff->Write();
    gChargedKstarUrQMDOff->Write();

    // Mean pT
    hPtMBPhiUrQMDOff->Write();
    hPtMBKstarUrQMDOff->Write();
    hPtMBXiStarUrQMDOff->Write();
    hPtMBChargedKstarUrQMDOff->Write();

    // Minimum-bias pT spectra
    gMeanPtPhiUrQMDOff->Write();
    gMeanPtKstarUrQMDOff->Write();
    gMeanPtXiStarUrQMDOff->Write();
    gMeanPtChargedKstarUrQMDOff->Write();

    fout->cd();
    TDirectory *dirNoRescattering = (TDirectory *)fout->mkdir("UrQMD_ON_NoRescattering");
    dirNoRescattering->cd();

    // Mean Yield
    gPhiNoRescattering->Write();
    gKstarNoRescattering->Write();
    gXiStarNoRescattering->Write();
    gChargedKstarNoRescattering->Write();

    // Mean pT
    gMeanPtPhiNoRescattering->Write();
    gMeanPtKstarNoRescattering->Write();
    gMeanPtXiStarNoRescattering->Write();
    gMeanPtChargedKstarNoRescattering->Write();

    // Minimum-bias pT spectra
    hPtMBPhiNoRescattering->Write();
    hPtMBKstarNoRescattering->Write();
    hPtMBXiStarNoRescattering->Write();
    hPtMBChargedKstarNoRescattering->Write();

    fout->cd();
    TDirectory *dirRescattering = (TDirectory *)fout->mkdir("UrQMD_ON_Rescattering");
    dirRescattering->cd();

    // Mean Yield
    gPhiRescattering->Write();
    gKstarRescattering->Write();
    gXiStarRescattering->Write();
    gChargedKstarRescattering->Write();

    // Mean pT
    gMeanPtPhiRescattering->Write();
    gMeanPtKstarRescattering->Write();
    gMeanPtXiStarRescattering->Write();
    gMeanPtChargedKstarRescattering->Write();

    // Minimum-bias pT spectra
    hPtMBPhiRescattering->Write();
    hPtMBKstarRescattering->Write();
    hPtMBXiStarRescattering->Write();
    hPtMBChargedKstarRescattering->Write();

    fout->Close();
    cout << "\nSaved EPOS_final.root.\n";
}
