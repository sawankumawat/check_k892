#include <vector>
#include <iostream>
#include <map>
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMath.h"
using namespace std;

// Changes in the code
// 1. Modified the pT cut from 0.15 to 0.1 GeV
// 2. For Pi,K,p and K0s, instead of checking first 3 EPOS digit, we check full epos ID only
// 3. For Pi,K,p, we choose the rapidity region of |y|<0.3 and |eta|<0.8 just like data and then later correct for it during dN/dy calculation. For other particles, we keep the same rapidity cut of |y|<0.5

// In charged-particle pseudorapidity analysis, no pT cut is used.

// static const int NCENT = 10;
static const int NCENT = 21;
static const int NIST = 10;

struct Species
{
    string name;
    vector<int> eposIDs; // matching EPOS ids
    double BR;

    double yieldIST[NIST][NCENT] = {{0}};
    double yieldIST9_ITY80[NCENT] = {0};
    double yieldIST9_ITY81[NCENT] = {0};

    double sumPtIST[NIST][NCENT] = {{0}};
    double sumPtIST9_ITY80[NCENT] = {0};
    double sumPtIST9_ITY81[NCENT] = {0};

    double countPtIST[NIST][NCENT] = {{0}};
    double countPtIST9_ITY80[NCENT] = {0};
    double countPtIST9_ITY81[NCENT] = {0};

    double yieldNoIST[NCENT] = {0};
    double sumPtNoIST[NCENT] = {0};
    double countPtNoIST[NCENT] = {0};

    double yieldIST0or8[NCENT] = {0};
    double sumPtIST0or8[NCENT] = {0};
    double countPtIST0or8[NCENT] = {0};

    TH1F *hPtCentIST[NIST][NCENT] = {{nullptr}};
    TH1F *hPtCentIST9_ITY80[NCENT] = {nullptr};
    TH1F *hPtCentIST9_ITY81[NCENT] = {nullptr};
    TH1F *hPtCentNoIST[NCENT] = {nullptr};
    TH1F *hPtCentIST0or8[NCENT] = {nullptr};
};

int GetCentralityBin(double cent)
{
    if (cent < 0.5)
        return 0;
    else if (cent < 1)
        return 1;
    else if (cent < 2)
        return 2;
    else if (cent < 3)
        return 3;
    else if (cent < 4)
        return 4;
    else if (cent < 5)
        return 5;
    else if (cent < 8)
        return 6;
    else if (cent < 10)
        return 7;
    else if (cent < 12)
        return 8;
    else if (cent < 15)
        return 9;
    else if (cent < 18)
        return 10;
    else if (cent < 20)
        return 11;
    else if (cent < 25)
        return 12;
    else if (cent < 30)
        return 13;
    else if (cent < 35)
        return 14;
    else if (cent < 40)
        return 15;
    else if (cent < 45)
        return 16;
    else if (cent < 50)
        return 17;
    else if (cent < 60)
        return 18;
    else if (cent < 70)
        return 19;
    else
        return 20;
}

int findSpeciesIndex(const vector<Species> &species, int eposID)
{
    for (size_t i = 0; i < species.size(); ++i)
    {
        // If the species if pion,kaon, proton and Kshort, then do not multiply by 100
        if (eposID == 120 || eposID == 130 || eposID == 1120 || eposID == 20)
        {
            if (eposID == species[i].eposIDs[0])
                return (int)i;
        }
        else
        {
            for (int id : species[i].eposIDs)
            {
                if (eposID == id || eposID == id * 100)
                    return (int)i;
            }
        }
    }
    return -1;
}

void EPOS_finalCodeQA()
{
    TDatabasePDG *pdg = new TDatabasePDG();

    TChain chain("teposevent0");
    // chain.Add("/home/sawan/Storage/EPOS_localOutputs/mergedEPOS_UrQMD.root");
    chain.Add("/home/sawan/Storage/EPOS_localOutputs/mergedEPOS_UrQMD2.root");
    // chain.Add("/home/sawan/Storage/EPOS_localOutputs/merged_EPOSV0.root");

    Long64_t nEv = chain.GetEntries();
    cout << "Total events = " << nEv << endl;
    nEv = 1000; // for quick test, comment out for full processing

    Int_t np;
    vector<Float_t> px(200000), py(200000), pz(200000), e(200000);
    vector<Int_t> id(200000), ist(200000), ity(200000);

    chain.SetBranchAddress("np", &np);
    chain.SetBranchAddress("px", &px[0]);
    chain.SetBranchAddress("py", &py[0]);
    chain.SetBranchAddress("pz", &pz[0]);
    chain.SetBranchAddress("e", &e[0]);
    chain.SetBranchAddress("id", &id[0]);
    chain.SetBranchAddress("ist", &ist[0]);
    chain.SetBranchAddress("ity", &ity[0]);

    TH1F *hFT0 = new TH1F("hFT0", "FT0 multiplicity", 500, 0, 500);
    vector<int> eventMult(nEv);
    TH1F *hIST = new TH1F("hIST", "IST distribution", 20, -0.5, 19.5);
    TH1F *hITY = new TH1F("hITY", "ITY distribution", 100, -0.5, 99.5);
    TH2F *hRapvsIST = new TH2F("hRapvsIST", "Rapidity vs IST;IST;y", 20, -0.5, 19.5, 100, -1, 1);
    TH2F *hIdvsIST = new TH2F("hIdvsIST", "EPOS ID vs IST;IST;EPOS ID", 20, -0.5, 19.5, 3000, -1500, 1500);

    cout << "First pass: computing FT0..." << endl;

    // species list
    vector<Species> species;
    auto addSpecies = [&](const string &name, const vector<int> &ids, double BR)
    {
        Species s;
        s.name = name;
        s.eposIDs = ids;
        s.BR = BR;
        string title = name + " p_{T};p_{T} (GeV/c);Counts";
        string baseName = string("hPtMB_") + name;
        for (int istBin = 0; istBin < NIST; ++istBin)
        {
            for (int icent = 0; icent < NCENT; ++icent)
            {
                s.hPtCentIST[istBin][icent] = new TH1F(Form("%s_IST%d_Cent%d", baseName.c_str(), istBin, icent), title.c_str(), 200, 0, 20);

                s.hPtCentIST[istBin][icent]->SetDirectory(nullptr);
            }
        }

        for (int icent = 0; icent < NCENT; ++icent)
        {
            s.hPtCentIST9_ITY80[icent] = new TH1F(Form("%s_IST9_ITY80_Cent%d", baseName.c_str(), icent), title.c_str(), 200, 0, 20);

            s.hPtCentIST9_ITY81[icent] = new TH1F(Form("%s_IST9_ITY81_Cent%d", baseName.c_str(), icent), title.c_str(), 200, 0, 20);

            s.hPtCentNoIST[icent] = new TH1F(Form("%s_NoIST_Cent%d", baseName.c_str(), icent), title.c_str(), 200, 0, 20);

            s.hPtCentIST0or8[icent] = new TH1F(Form("%s_IST0or8_Cent%d", baseName.c_str(), icent), title.c_str(), 200, 0, 20);

            s.hPtCentIST9_ITY80[icent]->SetDirectory(nullptr);
            s.hPtCentIST9_ITY81[icent]->SetDirectory(nullptr);
            s.hPtCentNoIST[icent]->SetDirectory(nullptr);
            s.hPtCentIST0or8[icent]->SetDirectory(nullptr);
        }

        species.push_back(s);
    };

    // add species with EPOS id prefixes and BR
    addSpecies("phi", {331}, 0.492);
    addSpecies("kstar", {231}, 0.666);
    addSpecies("kshort", {20}, 1.0);
    addSpecies("proton", {1120}, 1.0);
    addSpecies("pion", {120}, 1.0);
    addSpecies("kaon", {130}, 1.0);
    addSpecies("xistar", {1331}, 0.666);
    addSpecies("chargedkstar", {131}, 0.333);

    // EPOS ids used for FT0 selection (final hadrons IDs)
    vector<int> fftIDs = {120, -120, 130, -130, 1120, -1120};

    // FIRST PASS
    for (Long64_t ev = 0; ev < nEv; ev++)
    {
        chain.GetEntry(ev);
        int multFT0 = 0;
        for (int i = 0; i < np; i++)
        {
            int eposID = TMath::Abs(id[i]);
            if (eposID != 120 && eposID != 130 && eposID != 1120) // checking Pi,K,P
                continue;

            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
            // if (pt < 0.1)
            //     continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));

            if ((eta > 3.5 && eta < 4.9) || (eta > -3.3 && eta < -2.1))
                if (ist[i] == 0)
                    multFT0++;
        }
        hFT0->Fill(multFT0);
        eventMult[ev] = multFT0;
    }

    // multiplicity -> percentile
    vector<double> multToCent(hFT0->GetNbinsX() + 1, 100);
    double total = hFT0->Integral();
    double cum = 0;
    for (int b = hFT0->GetNbinsX(); b >= 1; b--)
    {
        cum += hFT0->GetBinContent(b);
        multToCent[b] = 100.0 * cum / total;
    }

    // storage for centroids
    double sumDNdEta[NCENT] = {0};
    int nEventsCent[NCENT] = {0};

    cout << "Second pass: physics observables..." << endl;

    // SECOND PASS
    for (Long64_t ev = 0; ev < nEv; ev++)
    {
        chain.GetEntry(ev);
        int multFT0 = eventMult[ev];
        int bin = hFT0->FindBin(multFT0);
        double cent = multToCent[bin];
        int cbin = GetCentralityBin(cent);
        // nEventsCent[cbin]++;

        int nchMid = 0;
        int atLeastOnePiKp_in_ModEta1 = 0;

        for (int i = 0; i < np; i++)
        {
            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
            // if (pt < 0.1)
            //     continue;

            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));
            int eposID = TMath::Abs(id[i]);
            if (eposID == 120 || eposID == 130 || eposID == 1120)
            {
                if (fabs(eta) < 1.0 && ist[i] == 0)
                    atLeastOnePiKp_in_ModEta1++;
            }
        }

        if (atLeastOnePiKp_in_ModEta1 > 0)
        {
            nEventsCent[cbin]++;
        }

        for (int i = 0; i < np; i++)
        {
            double p = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
            if (p == fabs(pz[i]))
                continue;

            double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
            double eta = 0.5 * log((p + pz[i]) / (p - pz[i]));
            double y = 0.5 * log((e[i] + pz[i]) / (e[i] - pz[i]));
            int eposID = TMath::Abs(id[i]);

            // charged multiplicity
            // if ((eposID == 120 || eposID == 130 || eposID == 1120) && fabs(eta) < 0.5 && ist[i] == 0 && pt >= 0.1)
            if ((eposID == 120 || eposID == 130 || eposID == 1120) && fabs(eta) < 0.5 && ist[i] == 0)
                nchMid++;

            hIST->Fill(ist[i]);
            hITY->Fill(ity[i]);
            hRapvsIST->Fill(ist[i], y);
            hIdvsIST->Fill(ist[i], eposID);

            // if (pt < 0.1)
            //     continue;

            bool isPiKp = (eposID == 120 || eposID == 130 || eposID == 1120);

            bool acceptParticle = false;

            if (isPiKp)
            {
                // if (fabs(y) < 0.3 && fabs(eta) < 0.8)
                // if (fabs(y) < 0.5 && fabs(eta) < 0.8)
                if (fabs(y) < 0.5)
                    acceptParticle = true;
            }
            else
            {
                if (fabs(y) < 0.5)
                    acceptParticle = true;
            }

            // if (fabs(y) < 0.5 && atLeastOnePiKp_in_ModEta1 > 0)
            if (acceptParticle && atLeastOnePiKp_in_ModEta1 > 0)
            {
                int sidx = findSpeciesIndex(species, eposID);
                if (sidx >= 0)
                {

                    Species &s = species[sidx];
                    s.yieldNoIST[cbin] += 1;
                    s.sumPtNoIST[cbin] += pt;
                    s.countPtNoIST[cbin] += 1;

                    s.hPtCentNoIST[cbin]->Fill(pt);

                    if (ist[i] == 0 || ist[i] == 8)
                    {
                        s.yieldIST0or8[cbin] += 1;
                        s.sumPtIST0or8[cbin] += pt;
                        s.countPtIST0or8[cbin] += 1;
                        s.hPtCentIST0or8[cbin]->Fill(pt);
                    }

                    if (ist[i] >= 0 && ist[i] < NIST)
                    {
                        s.yieldIST[ist[i]][cbin] += 1;
                        s.sumPtIST[ist[i]][cbin] += pt;
                        s.countPtIST[ist[i]][cbin] += 1;
                        s.hPtCentIST[ist[i]][cbin]->Fill(pt);
                    }
                    // ist==9 and ity==80 -> IST9_ITY80
                    if (ist[i] == 9 && ity[i] == 80)
                    {
                        s.yieldIST9_ITY80[cbin] += 1;
                        s.sumPtIST9_ITY80[cbin] += pt;
                        s.countPtIST9_ITY80[cbin] += 1;
                        s.hPtCentIST9_ITY80[cbin]->Fill(pt);
                    }
                    // ist==9 and ity==81 -> IST9_ITY81
                    if (ist[i] == 9 && ity[i] == 81)
                    {
                        s.yieldIST9_ITY81[cbin] += 1;
                        s.sumPtIST9_ITY81[cbin] += pt;
                        s.countPtIST9_ITY81[cbin] += 1;
                        s.hPtCentIST9_ITY81[cbin]->Fill(pt);
                    }
                }
            }
        }
        if (atLeastOnePiKp_in_ModEta1 > 0)
        {
            sumDNdEta[cbin] += nchMid;
        }
    }

    // Normalize and compute meanPt
    // branching ratios were set per species

    // Normalize and compute meanPt

    vector<double> meanDNdEta(NCENT, 0.0);
    for (int i = 0; i < NCENT; ++i)
    {
        if (nEventsCent[i] > 0)
            meanDNdEta[i] = sumDNdEta[i] / nEventsCent[i];
    }

    //==========================================================
    // Normalize yields and calculate <pT>
    //==========================================================

    for (auto &s : species)
    {
        double deltay = 1.0;
        // if (s.name == "pion" || s.name == "kaon" || s.name == "proton")
        //     deltay = 0.6;

        for (int i = 0; i < NCENT; i++)
        {
            if (nEventsCent[i] == 0)
                continue;

            // ---------- normalize dN/dy ----------
            for (int istBin = 0; istBin < NIST; istBin++)
                s.yieldIST[istBin][i] /= (nEventsCent[i] * s.BR * deltay);

            s.yieldIST9_ITY80[i] /= (nEventsCent[i] * s.BR * deltay);
            s.yieldIST9_ITY81[i] /= (nEventsCent[i] * s.BR * deltay);
            s.yieldNoIST[i] /= (nEventsCent[i] * s.BR * deltay);
            s.yieldIST0or8[i] /= (nEventsCent[i] * s.BR * deltay);

            // ---------- calculate <pT> ----------
            for (int istBin = 0; istBin < NIST; istBin++)
            {
                if (s.countPtIST[istBin][i] > 0)
                    s.sumPtIST[istBin][i] /= s.countPtIST[istBin][i];
            }

            if (s.countPtIST9_ITY80[i] > 0)
                s.sumPtIST9_ITY80[i] /= s.countPtIST9_ITY80[i];

            if (s.countPtIST9_ITY81[i] > 0)
                s.sumPtIST9_ITY81[i] /= s.countPtIST9_ITY81[i];

            if (s.countPtNoIST[i] > 0)
                s.sumPtNoIST[i] /= s.countPtNoIST[i];

            if (s.countPtIST0or8[i] > 0)
                s.sumPtIST0or8[i] /= s.countPtIST0or8[i];
        }
    }

    //==========================================================
    // Normalize pT histograms
    //==========================================================

    for (auto &s : species)
    {
        double deltay = 1.0;
        // if (s.name=="pion" || s.name=="kaon" || s.name=="proton")
        //     deltay = 0.6;

        for (int icent = 0; icent < NCENT; ++icent)
        {
            if (nEventsCent[icent] == 0)
                continue;

            double norm = nEventsCent[icent] * s.BR * deltay;

            for (int istBin = 0; istBin < NIST; ++istBin)
            {
                s.hPtCentIST[istBin][icent]->Scale(
                    1.0 / (norm * s.hPtCentIST[istBin][icent]->GetBinWidth(1)));
            }

            s.hPtCentIST9_ITY80[icent]->Scale(
                1.0 / (norm * s.hPtCentIST9_ITY80[icent]->GetBinWidth(1)));

            s.hPtCentIST9_ITY81[icent]->Scale(
                1.0 / (norm * s.hPtCentIST9_ITY81[icent]->GetBinWidth(1)));

            s.hPtCentNoIST[icent]->Scale(
                1.0 / (norm * s.hPtCentNoIST[icent]->GetBinWidth(1)));

            s.hPtCentIST0or8[icent]->Scale(
                1.0 / (norm * s.hPtCentIST0or8[icent]->GetBinWidth(1)));
        }
    }

    // prepare output file and write histograms/graphs
    TFile *fout = new TFile("EPOS_finalQA_ptCut.root", "RECREATE");
    hFT0->Write();
    hIST->Write();
    hITY->Write();
    hRapvsIST->Write();
    hIdvsIST->Write();

    // create graphs per species and write them into directories named by IST/ITY
    TDirectory *dirIST[NIST];
    for (int istBin = 0; istBin < NIST; ++istBin)
    {
        dirIST[istBin] = (TDirectory *)fout->mkdir((string("IST") + to_string(istBin)).c_str());
    }
    TDirectory *dirIST9_ITY80 = (TDirectory *)fout->mkdir("IST9_ITY80");
    TDirectory *dirIST9_ITY81 = (TDirectory *)fout->mkdir("IST9_ITY81");
    TDirectory *dirNoIST = (TDirectory *)fout->mkdir("NoIST");
    TDirectory *dirIST0or8 = (TDirectory *)fout->mkdir("IST0or8");

    for (auto &s : species)
    {
        double x[NCENT];
        double yIST[NIST][NCENT];
        double y_no[NCENT];
        double y_res[NCENT];
        double yieldNoIST[NCENT];
        double meanNoIST[NCENT];
        double yieldIST0or8[NCENT];
        double meanIST0or8[NCENT];
        for (int i = 0; i < NCENT; i++)
        {
            x[i] = meanDNdEta[i];
            for (int istBin = 0; istBin < NIST; ++istBin)
                yIST[istBin][i] = s.yieldIST[istBin][i];
            y_no[i] = s.yieldIST9_ITY80[i];
            y_res[i] = s.yieldIST9_ITY81[i];
            yieldNoIST[i] = s.yieldNoIST[i];
            meanNoIST[i] = s.sumPtNoIST[i];
            yieldIST0or8[i] = s.yieldIST0or8[i];
            meanIST0or8[i] = s.sumPtIST0or8[i];
        }

        // graphs (use same base name; will be written into separate directories)
        TGraph *gNo = new TGraph(NCENT, x, y_no);
        TGraph *gRes = new TGraph(NCENT, x, y_res);
        string gbase = s.name + "_vs_mult";
        gNo->SetName(gbase.c_str());
        gRes->SetName(gbase.c_str());

        // Write MB histograms and graphs into respective IST/ITY-named directories
        for (int istBin = 0; istBin < NIST; ++istBin)
        {
            TGraph *gIST = new TGraph(NCENT, x, yIST[istBin]);
            gIST->SetName(gbase.c_str());
            fout->cd();
            dirIST[istBin]->cd();
            for (int icent = 0; icent < NCENT; ++icent)
            {
                s.hPtCentIST[istBin][icent]->Write();
            }
            gIST->Write();
            // also write mean-pT graph for this IST bin
            double meanPtArr[NCENT];
            for (int ii = 0; ii < NCENT; ++ii)
            {
                meanPtArr[ii] = s.sumPtIST[istBin][ii];
            }
            TGraph *gMeanPtIST = new TGraph(NCENT, x, meanPtArr);
            string meanName = string("meanpt_") + s.name + string("_vs_mult");
            gMeanPtIST->SetName(meanName.c_str());
            gMeanPtIST->Write();
            delete gMeanPtIST;
            delete gIST;
        }

        // IST9_ITY80: write MB histogram, yield graph and mean-pT graph
        fout->cd();
        dirIST9_ITY80->cd();
        for (int icent = 0; icent < NCENT; ++icent)
        {
            s.hPtCentIST9_ITY80[icent]->Write();
        }
        gNo->Write();
        double meanNo[NCENT];
        for (int ii = 0; ii < NCENT; ++ii)
        {
            meanNo[ii] = s.sumPtIST9_ITY80[ii];
        }
        TGraph *gMeanNo = new TGraph(NCENT, x, meanNo);
        string meanNameNo = string("meanpt_") + s.name + string("_vs_mult");
        gMeanNo->SetName(meanNameNo.c_str());
        gMeanNo->Write();
        delete gMeanNo;

        // IST9_ITY81: write MB histogram, yield graph and mean-pT graph
        fout->cd();
        dirIST9_ITY81->cd();
        for (int icent = 0; icent < NCENT; ++icent)
        {
            s.hPtCentIST9_ITY81[icent]->Write();
        }
        gRes->Write();
        double meanRes[NCENT];
        for (int ii = 0; ii < NCENT; ++ii)
        {
            meanRes[ii] = s.sumPtIST9_ITY81[ii];
        }
        TGraph *gMeanRes = new TGraph(NCENT, x, meanRes);
        string meanNameRes = string("meanpt_") + s.name + string("_vs_mult");
        gMeanRes->SetName(meanNameRes.c_str());
        gMeanRes->Write();
        delete gMeanRes;

        dirNoIST->cd();

        TGraph *gYieldNoIST = new TGraph(NCENT, x, yieldNoIST);
        gYieldNoIST->SetName((s.name + "_vs_mult").c_str());
        gYieldNoIST->Write();

        TGraph *gMeanPtNoIST = new TGraph(NCENT, x, meanNoIST);
        gMeanPtNoIST->SetName((string("meanpt_") + s.name + "_vs_mult").c_str());
        gMeanPtNoIST->Write();

        for (int icent = 0; icent < NCENT; ++icent)
        {
            s.hPtCentNoIST[icent]->Write();
        }

        delete gYieldNoIST;
        delete gMeanPtNoIST;

        dirIST0or8->cd();

        TGraph *gYieldIST0or8 = new TGraph(NCENT, x, yieldIST0or8);
        gYieldIST0or8->SetName((s.name + "_vs_mult").c_str());
        gYieldIST0or8->Write();

        TGraph *gMeanPtIST0or8 = new TGraph(NCENT, x, meanIST0or8);
        gMeanPtIST0or8->SetName((string("meanpt_") + s.name + "_vs_mult").c_str());
        gMeanPtIST0or8->Write();

        for (int icent = 0; icent < NCENT; ++icent)
        {
            s.hPtCentIST0or8[icent]->Write();
        }

        delete gYieldIST0or8;
        delete gMeanPtIST0or8;

        fout->cd();
    }

    fout->Close();
    cout << "\nSaved EPOS_finalQA_ptCut.root.\n";
}
