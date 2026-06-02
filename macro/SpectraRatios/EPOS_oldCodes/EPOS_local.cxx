#include <cmath>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TString.h>
#include <TTree.h>
#include <TSystem.h>
using namespace std;

// =====================================================
// centrality bin from percentile
// =====================================================
int GetCentralityBin(double cent)
{
    if (cent < 10)
        return 0;
    else if (cent < 20)
        return 1;
    else if (cent < 30)
        return 2;
    else if (cent < 40)
        return 3;
    else if (cent < 50)
        return 4;
    else if (cent < 60)
        return 5;
    else if (cent < 70)
        return 6;
    else if (cent < 80)
        return 7;
    else if (cent < 90)
        return 8;
    else
        return 9;
}

// =====================================================
// MAIN
// =====================================================
void EPOS_local()
{
    TFile *inputFile = new TFile("/home/sawan/Storage/EPOS_localOutputs/mergedEPOS_UrQMD.root", "READ");
    TTree *tree = (TTree *)inputFile->Get("teposevent0");
    if (!tree)
    {
        cout << "Error: could not find tree teposevent0 in the root file" << endl;
        return;
    }

    // Long64_t nEv = tree->GetEntries();
    // cout << "Total events = " << nEv << endl;
    // int nEv = 500; // for testing

    // =====================================================
    // branches
    // =====================================================
    Int_t np;
    Float_t px[200000], py[200000], pz[200000];
    Float_t e[200000];
    Int_t id[200000];
    Int_t ist[200000];

    tree->SetBranchAddress("np", &np);
    tree->SetBranchAddress("px", px);
    tree->SetBranchAddress("py", py);
    tree->SetBranchAddress("pz", pz);
    tree->SetBranchAddress("e", e);
    tree->SetBranchAddress("id", id);
    tree->SetBranchAddress("ist", ist);

    // =====================================================
    // FT0 multiplicity
    // =====================================================
    TH1F *hFT0 = new TH1F("hFT0", "FT0 multiplicity", 500, 0, 500);
    vector<int> eventMult(nEv);

    cout << "First pass: computing FT0..." << endl;

    // =====================================================
    // FIRST PASS → FT0 multiplicity
    // =====================================================
    for (Long64_t ev = 0; ev < nEv; ev++)
    {
        tree->GetEntry(ev);
        int multFT0 = 0;

        for (int i = 0; i < np; i++)
        {
            int eposID = TMath::Abs(id[i]);
            // bool isCharged =
            //     (eposID == 120) ||  // pi
            //     (eposID == 130) ||  // K
            //     (eposID == 1120) || // proton
            //     (eposID == 20) ||   // electron
            //     (eposID == 14) ||   // muon
            //     (eposID == 110) ||  // Sigma
            //     (eposID == 240) ||  // Xi
            //     (eposID == 1220);   // Omega

            // if (!isCharged)
            //     continue;

            if (eposID != 120 && eposID != 130 && eposID != 1120) // Pi, K and Pr
                continue;

            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));

            if ((eta > 3.3 && eta < 5.0) || (eta > -3.3 && eta < -2.1))
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

    const int nSpecies = 8;
    const char *speciesName[nSpecies] = {"pion", "kaon", "proton", "phi", "kstar", "k0s", "xistar", "chKstar"};
    const int speciesId[nSpecies] = {120, 130, 1120, 331, 231, 20, 1331, 131};
    const double speciesBR[nSpecies] = {1.0, 1.0, 1.0, 0.492, 0.666, 1.0, 0.666, 0.333};

    double speciesYield[nSpecies][10] = {{0}};
    double speciesMeanPt[nSpecies][10] = {{0}};
    TH1F *hSpeciesPt[nSpecies][10] = {{nullptr}};

    for (int s = 0; s < nSpecies; s++)
    {
        for (int i = 0; i < 10; i++)
        {
            hSpeciesPt[s][i] = new TH1F(Form("h%sPt_bin%d", speciesName[s], i),
                                        Form("%s p_{T} bin %d;p_{T} (GeV/c);Counts", speciesName[s], i),
                                        100, 0, 10);
        }
    }

    cout << "Second pass: physics observables..." << endl;

    // =====================================================
    // SECOND PASS
    // =====================================================
    for (Long64_t ev = 0; ev < nEv; ev++)
    {
        tree->GetEntry(ev);

        int multFT0 = eventMult[ev];
        int bin = hFT0->FindBin(multFT0);
        double cent = multToCent[bin];
        int cbin = GetCentralityBin(cent);

        nEventsCent[cbin]++;

        int nchMid = 0;

        for (int i = 0; i < np; i++)
        {
            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));
            // if ((e[i] - fabs(pz[i])) <= 1e-10)
            //     continue;
            double y = 0.5 * log((e[i] + pz[i]) / (e[i] - pz[i]));

            int eposID = TMath::Abs(id[i]);

            // ----------------------
            // charged multiplicity
            // ----------------------

            if (eposID == 120 || eposID == 130 || eposID == 1120)
                if (fabs(eta) < 0.5)
                    nchMid++;

            // if (ist[i] == 7)
            // {
            //     bool isCharged =
            //         (eposID == 120) ||  // pi
            //         (eposID == 130) ||  // K
            //         (eposID == 1120) || // proton
            //         (eposID == 20) ||   // electron
            //         (eposID == 14) ||   // muon
            //         (eposID == 110) ||  // Sigma
            //         (eposID == 240) ||  // Xi
            //         (eposID == 1220);   // Omega

            //     if (isCharged && fabs(eta) < 0.5)
            //         nchMid++;
            // }
            // =================================================
            // final charged particles
            // =================================================

            // =================================================
            // URQMD OFF → final hadrons
            // =================================================
            // // if (id[i] != 120 && id[i] != 130 && id[i] != 1120)
            // {
            //     cout << "Event " << ev << "  Track " << i
            //          << "  ID=" << id[i]
            //          << "  ist=" << ist[i]
            //          << "  pT=" << sqrt(px[i] * px[i] + py[i] * py[i])
            //          << "  y=" << y
            //          << endl;
            // }

            if (ist[i] != 7) // ist 7 = UrQMD off, ist = 9 UrQMD on
                continue;

            const double pT = sqrt(px[i] * px[i] + py[i] * py[i]);
            for (int s = 0; s < nSpecies; s++)
            {
                if (TMath::Abs(id[i]) == speciesId[s])
                {
                    if (fabs(y) < 0.5)
                    {
                        speciesYield[s][cbin]++;
                        hSpeciesPt[s][cbin]->Fill(pT);
                    }
                }
            }
        }

        sumDNdEta[cbin] += nchMid;
    }

    // =====================================================
    // branching ratios
    // =====================================================
    double meanDNdEta[10] = {0};

    for (int i = 0; i < 10; i++)
    {
        if (nEventsCent[i] > 0)
        {
            meanDNdEta[i] = sumDNdEta[i] / nEventsCent[i];

            for (int s = 0; s < nSpecies; s++)
            {
                speciesYield[s][i] = (speciesYield[s][i] / nEventsCent[i]) / speciesBR[s];
                speciesMeanPt[s][i] = hSpeciesPt[s][i]->GetEntries() > 0 ? hSpeciesPt[s][i]->GetMean() : 0;

                const double norm = nEventsCent[i] * speciesBR[s] * hSpeciesPt[s][i]->GetBinWidth(1);
                if (norm > 0)
                    hSpeciesPt[s][i]->Scale(1.0 / norm);
            }
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
             << endl;

        for (int s = 0; s < nSpecies; s++)
        {
            cout << "    " << speciesName[s]
                 << "=" << speciesYield[s][i]
                 << "  <pT>=" << speciesMeanPt[s][i];
            if (s < nSpecies - 1)
                cout << "  |";
        }
        cout << endl;
    }

    // =====================================================
    // graphs
    // =====================================================
    TGraph *gYield[nSpecies];
    TGraph *gMeanPt[nSpecies];
    for (int s = 0; s < nSpecies; s++)
    {
        gYield[s] = new TGraph(10, meanDNdEta, speciesYield[s]);
        gMeanPt[s] = new TGraph(10, meanDNdEta, speciesMeanPt[s]);
        gYield[s]->SetName(Form("%s_vs_mult", speciesName[s]));
        gMeanPt[s]->SetName(Form("%s_meanpt_vs_mult", speciesName[s]));
    }

    // =====================================================
    // save
    // =====================================================
    TFile fout("resonance_vs_mult_urqmdOFF.root", "RECREATE");
    hFT0->Write();
    for (int s = 0; s < nSpecies; s++)
    {
        gYield[s]->Write();
        gMeanPt[s]->Write();
        for (int i = 0; i < 10; i++)
        {
            hSpeciesPt[s][i]->Write();
        }
    }
    fout.Close();

    cout << "\nSaved resonance_vs_mult_urqmdOFF.root\n";
}
