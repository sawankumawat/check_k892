#include <vector>
#include <iostream>
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

// =====================================================
// MAIN
// =====================================================
void epos()
{
    // =====================================================
    // READ ALL FILES
    // =====================================================
    TChain chain("teposevent0");
    chain.Add("/home/sawan/Storage/EPOS_localOutputs/mergedEPOS_UrQMD.root");
    // chain.Add("/home/sawan/merged.root");

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

    chain.SetBranchAddress("np", &np);
    chain.SetBranchAddress("px", px);
    chain.SetBranchAddress("py", py);
    chain.SetBranchAddress("pz", pz);
    chain.SetBranchAddress("e", e);
    chain.SetBranchAddress("id", id);
    chain.SetBranchAddress("ist", ist);

    // =====================================================
    // FT0 multiplicity
    // =====================================================
    TH1F *hFT0 = new TH1F("hFT0", "FT0 multiplicity", 5000, 0, 5000);
    vector<int> eventMult(nEv);

    cout << "First pass: computing FT0..." << endl;

    // =====================================================
    // FIRST PASS → FT0 multiplicity
    // =====================================================
    for (Long64_t ev = 0; ev < nEv; ev++)
    {
        chain.GetEntry(ev);
        int multFT0 = 0;

        for (int i = 0; i < np; i++)
        {
            int pid = TMath::Abs(id[i]);
            if (pid != 120 && pid != 130 && pid != 1120)
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

        for (int i = 0; i < np; i++)
        {
            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));
            double y = 0.5 * log((e[i] + pz[i]) / (e[i] - pz[i]));

            int pid = TMath::Abs(id[i]);

            // ----------------------
            // charged multiplicity
            // ----------------------
            if (pid == 120 || pid == 130 || pid == 1120)
                if (fabs(eta) < 0.5)
                    nchMid++;

            // =================================================
            // URQMD OFF → final hadrons
            // =================================================
            if (ist[i] != 7) // ist 7 = UrQMD off, ist = 9 UrQMD on
                continue;

            // ----------------------
            // PHI
            // ----------------------
            if (TMath::Abs(id[i]) == 331)
                if (fabs(y) < 0.5)
                    phiYield[cbin]++;

            // ----------------------
            // KSTAR
            // ----------------------
            if (TMath::Abs(id[i]) == 231)
                if (fabs(y) < 0.5)
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
            kstarYield[i] = kstarYield[i] / (nEventsCent[i] * BR_kstar * 2); // Average of K*0 and K*0bar
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
    gPhi->Write();
    gK->Write();
    fout.Close();

    cout << "\nSaved resonance_vs_mult_urqmdOFF.root\n";
}
