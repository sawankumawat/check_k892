#include <vector>
#include <iostream>
#include <algorithm>
#include "TDatabasePDG.h"
using namespace std;

// =====================================================
// MAIN
// =====================================================
void epos_modified2()
{
    TDatabasePDG *pdg = new TDatabasePDG();

    // =====================================================
    // READ ALL FILES
    // =====================================================
    TChain chain("teposevent0");
    chain.Add("/home/sawan/Storage/EPOS_localOutputs/mergedEPOS_UrQMD.root");

    Long64_t nEv = chain.GetEntries();
    cout << "Total events = " << nEv << endl;
    nEv = 100000; // for testing

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
    TH1F *hFT0 = new TH1F("hFT0", "FT0 multiplicity", 500, 0, 500);
    TH1F *hFT0Eta0p5 = new TH1F("hFT0Eta0p5", "FT0 multiplicity with |eta|<0.5 cut", 500, 0, 500);
    TH2F *hFT0_vs_FT0Eta0p5 = new TH2F("hFT0_vs_FT0Eta0p5", "FT0 vs FT0 with |eta|<0.5 cut", 500, 0, 500, 500, 0, 500);

    vector<int> eventMult(nEv);
    vector<int> eventMultEta0p5(nEv);

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
        int multFT0Eta0p5 = 0;

        for (int i = 0; i < np; i++)
        {
            int eposID = TMath::Abs(id[i]);
            if (eposID != 120 && eposID != 130 && eposID != 1120)
                continue;

            if (ist[i] != 0) // for final state particles
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

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));

            if ((eta > 3.5 && eta < 4.9) || (eta > -3.3 && eta < -2.1))
                if (charge != 0)
                    multFT0++;

            if (abs(eta) < 0.5)
                if (charge != 0)
                    multFT0Eta0p5++;
        }

        hFT0->Fill(multFT0);
        hFT0Eta0p5->Fill(multFT0Eta0p5);
        hFT0_vs_FT0Eta0p5->Fill(multFT0, multFT0Eta0p5);

        eventMult[ev] = multFT0;
        eventMultEta0p5[ev] = multFT0Eta0p5;
    }

    // =====================================================
    // 10 percentile classes from the right edge of FT0
    // =====================================================
    vector<int> orderedEvents(nEv);
    for (Long64_t ev = 0; ev < nEv; ev++)
        orderedEvents[ev] = ev;

    sort(orderedEvents.begin(), orderedEvents.end(), [&](int a, int b)
         { return eventMult[a] > eventMult[b]; });

    const int nClasses = 10;
    double meanChMult[nClasses] = {0};
    double classX[nClasses] = {0};

    for (int c = 0; c < nClasses; c++)
    {
        Long64_t start = (nEv * c) / nClasses;
        Long64_t stop = (nEv * (c + 1)) / nClasses;

        if (start >= stop)
            continue;

        long long sumFT0Eta0p5 = 0;
        long long sumFT0 = 0;

        for (Long64_t idx = start; idx < stop; idx++)
        {
            int ev = orderedEvents[idx];
            sumFT0 += eventMult[ev];
            sumFT0Eta0p5 += eventMultEta0p5[ev];
        }

        const double nClassEvents = static_cast<double>(stop - start);
        meanChMult[c] = sumFT0Eta0p5 / nClassEvents;
        classX[c] = 5.0 + 10.0 * c;

        cout << "class " << c
             << " (" << 10 * c << "-" << 10 * (c + 1) << "% from highest FT0)"
             << "  <FT0>=" << (sumFT0 / nClassEvents)
             << "  <Nch(|eta|<0.5)>=" << meanChMult[c]
             << endl;
    }

    TGraph *gMeanChMult = new TGraph(nClasses, classX, meanChMult);
    gMeanChMult->SetName("mean_ch_mult_vs_ft0_percentile");
    gMeanChMult->SetTitle("Average charged multiplicity vs FT0 percentile;FT0 percentile class;#LTN_{ch}#GT_{|#eta|<0.5}");

    TH1F *hMeanChMult = new TH1F("hMeanChMult", "Average charged multiplicity vs FT0 percentile;FT0 percentile class;#LTN_{ch}#GT_{|#eta|<0.5}", 10, 0, 100);
    for (int c = 0; c < nClasses; c++)
        hMeanChMult->SetBinContent(c + 1, meanChMult[c]);

    // TFile fout("epos_modified2_mult.root", "RECREATE");
    // hFT0->Write();
    // hFT0Eta0p5->Write();
    // hFT0_vs_FT0Eta0p5->Write();
    // hMeanChMult->Write();
    // gMeanChMult->Write();
    // fout.Close();
}
