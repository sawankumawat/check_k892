#include <vector>
#include <iostream>
#include "TDatabasePDG.h"
using namespace std;

// =====================================================
// centrality bin from percentile
// =====================================================
// int GetCentralityBin(double cent)
// {
//     if (cent < 1)
//         return 0;
//     else if (cent < 5)
//         return 1;
//     else if (cent < 10)
//         return 2;
//     else if (cent < 15)
//         return 3;
//     else if (cent < 20)
//         return 4;
//     else if (cent < 30)
//         return 5;
//     else if (cent < 40)
//         return 6;
//     else if (cent < 50)
//         return 7;
//     else if (cent < 70)
//         return 8;
//     else
//         return 9;
// }

int GetCentralityBin(double cent)
{
    if (cent < 0.96)
        return 0;
    else if (cent < 4.80)
        return 1;
    else if (cent < 9.60)
        return 2;
    else if (cent < 14.42)
        return 3;
    else if (cent < 19.25)
        return 4;
    else if (cent < 28.92)
        return 5;
    else if (cent < 38.62)
        return 6;
    else if (cent < 48.38)
        return 7;
    else if (cent < 68.09)
        return 8;
    else
        return 9;
}

// =====================================================
// MAIN
// =====================================================
void epos_modified()
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

            int PID;

            for (int j = 0; j < 6; j++)
            {
                if (eposID == TMath::Abs(eposid[j]))
                {
                    PID = pdgid[j];
                    break;
                }
            }

            int charge = pdg->GetParticle(PID)->Charge() / 3.;

            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
            if (pt < 0.1)
                continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));

            if ((eta > 3.5 && eta < 4.9) || (eta > -3.3 && eta < -2.1))
                if (ist[i] == 0) // for final state particles
                    if (charge != 0)
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
                if (fabs(eta) < 1.0)
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
            if (pt < 0.1)
                continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));
            double y = 0.5 * log((e[i] + pz[i]) / (e[i] - pz[i]));

            int eposID = TMath::Abs(id[i]);

            int PID;
            int charge = 0;
            if (eposID == 120 || eposID == 130 || eposID == 1120)
            {
                for (int j = 0; j < 6; j++)
                {
                    if (eposID == TMath::Abs(eposid[j]))
                    {
                        PID = pdgid[j];
                        break;
                    }
                }

                charge = pdg->GetParticle(PID)->Charge() / 3.;
            }

            // ----------------------
            // charged multiplicity
            // ----------------------
            if (eposID == 120 || eposID == 130 || eposID == 1120)
                if (fabs(eta) < 0.5)
                    if (ist[i] == 0)
                        if (charge != 0)
                            nchMid++;

            // =================================================
            // URQMD OFF → final hadrons
            // =================================================
            hIST->Fill(ist[i]);
            hITY->Fill(ity[i]);
            if (ist[i] != 7) // ist 7 = UrQMD off, ist = 9 UrQMD on
                continue;

            // ----------------------
            // PHI
            // ----------------------
            if (TMath::Abs(id[i]) == 331)
                if (fabs(y) < 0.5)
                    if (atLeastOnePiKp_in_ModEta1 > 0) // INEL > 0
                        phiYield[cbin]++;

            // ----------------------
            // KSTAR
            // ----------------------
            if (TMath::Abs(id[i]) == 231)
                if (fabs(y) < 0.5)
                    if (atLeastOnePiKp_in_ModEta1 > 0) // INEL > 0
                        kstarYield[cbin]++;
        }

        sumDNdEta[cbin] += nchMid;
    }

    // =====================================================
    // branching ratios
    // =====================================================
    const double BR_phi = 0.492;
    const double BR_kstar = 0.666;

    double meanDNdEta[10];

    for (int i = 0; i < 10; i++)
    {
        if (nEventsCent[i] > 0)
        {
            meanDNdEta[i] = sumDNdEta[i] / nEventsCent[i];

            phiYield[i] = (phiYield[i] / nEventsCent[i]) / BR_phi;
            kstarYield[i] = (kstarYield[i] / nEventsCent[i]) / BR_kstar;
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
    TGraph *gK = new TGraph(10, meanDNdEta, kstarYield);

    gPhi->SetName("phi_vs_mult_urqmdOFF");
    gK->SetName("kstar_vs_mult_urqmdOFF");

    // =====================================================
    // save
    // =====================================================
    TFile fout("resonance_vs_mult_urqmdOFF.root", "RECREATE");
    hFT0->Write();
    hIST->Write();
    hITY->Write();
    gPhi->Write();
    gK->Write();
    fout.Close();

    cout << "\nSaved resonance_vs_mult_urqmdOFF.root\n";
}
