#include <iostream>
#include <fstream>
#include <cmath>
#include "src/style.h"

using namespace std;

double BreitWigner(double *x, double *par) // fitting function
{
    double fing1;
    double fing2;
    double fing3;
    double bk1;
    double bk;
    double bk3;
    fing1 = (0.5 * par[2] * par[1]) / ((TMath::Pi()) * ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]));
    return (fing1);
}

double Gaussian(double *x, double *par) // fitting function
{
    double fing1;
    double fing2;
    double fing3;
    double bk1;
    double bk;
    double bk3;
    fing1 = par[2] * (1 / sqrt(2 * TMath::Pi())) * (1 / par[1]) * exp(-0.5 * pow(((x[0] - par[0]) / par[1]), 2));
    return fing1;
}

void toy_model_ptSmear_prottay()
{

    // accepting histograms for K0s, kaon and pion eff---------------------------------------
    gStyle->SetOptStat(0);

    TFile *fInputFile = TFile::Open("K0sEff.root");

    TH1D *heff = (TH1D *)fInputFile->Get("K0s_efficiency");
    TH1D *hgenksh = (TH1D *)fInputFile->Get("K0s_generated_pT");
    TH1D *hrecksh = (TH1D *)fInputFile->Get("K0s_reconstructed_pT");
    //----------------------------------------------------------------------

    double NPT[7] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0}; // pt bins of f1 in data

    // histograms for efficiency caluclation of f1-------------------------
    // TProfile *GTE = new TProfile("GTE", "GTE", 6, NPT);
    TH1D *hf1gen = new TH1D("f1gen", "f2(1525) Gen", 6, NPT);
    TH1D *hf1rec = new TH1D("f1rec", "f2(1525) Rec", 6, NPT);
    // TH1D *hkshkaonmass = new TH1D("kshkaonmass", "kshkaonmass", 300, 0, 3);
    hf1gen->Sumw2();
    hf1rec->Sumw2();

    //////////////////////////////////////////////////////////

    // TLorentz Vectors**********************************
    TLorentzVector F1;
    TLorentzVector *KSh1;
    TLorentzVector *KSh2;

    TLorentzVector ksh1;
    TLorentzVector ksh2;
    TLorentzVector kshksh;

    TLorentzVector *piPlus1;
    TLorentzVector *piMinus1;
    TLorentzVector *piPlus2;
    TLorentzVector *piMinus2;

    TLorentzVector pionPlus1;
    TLorentzVector pionMinus1;
    TLorentzVector pionPlus2;
    TLorentzVector pionMinus2;

    ////////////////////////////////////////////////////////////

    // For decay*************************************
    TGenPhaseSpace event;
    Double_t weight;
    TGenPhaseSpace eventPi1;
    TGenPhaseSpace eventPi2;
    ////////////////////////////////////////////////////

    //***Breit Wigner and gaussian distributions generated*********************
    // for f2(1525)
    TF1 *fE = new TF1("fE", BreitWigner, 1.285, 1.765, 3);
    fE->FixParameter(0, 1.5173);
    fE->FixParameter(1, 0.0844);
    fE->FixParameter(2.0, 1); // amplitude

    // // for f0(1710)
    // TF1 *fE = new TF1("fE", BreitWigner, 1.283, 2.183, 3);
    // fE->FixParameter(0, 1.733);
    // fE->FixParameter(1, 0.150);
    // fE->FixParameter(2.0, 1); // amplitude

    // For K0s
    TF1 *f = new TF1("f", Gaussian, 0.48, 0.52, 3);
    f->FixParameter(0, 0.497);
    f->FixParameter(1, 0.005);
    f->FixParameter(2, 1.0); // amplitude

    //***************************************

    //*****variables***********************************
    double ptf1, etaf1, phif1, momentumf1, energyf1;
    TRandom *eventGenerator = new TRandom();
    TRandom3 *randomGenerator = new TRandom3();
    double massf1 = 0.0;
    double massksh1 = 0.0;
    double massksh2 = 0.0;
    double probacceptkshort1 = 0.0, probacceptkshort2 = 0.0;
    double kshpt1 = 0.0, kshpt2 = 0.0;
    int ptBinksh1 = 0, ptBinksh2 = 0;
    double effksh1 = 0.0, effksh2 = 0.0, finaleff = 0.0;
    double genkshcontent1 = 0.0, genkshcontent2 = 0.0;
    double reckshcontent1 = 0.0, reckshcontent2 = 0.0;

    // std::vector<double> kshKaonMassValues;
    // double kshKaonMass = 0.0;
    ///////////////////////////////////////////////////////

    for (int i = 0; i < 1000000; i++)
    {
        massf1 = fE->GetRandom(); // mass of f1 from breitwigner
        massksh1 = f->GetRandom();
        massksh2 = f->GetRandom();
        ptf1 = gRandom->Uniform(0.0, 15.0);  // pt
        etaf1 = gRandom->Uniform(-1.2, 1.2); // eta
        phif1 = gRandom->Uniform(0, 6.28);   // phi

        if (ptf1 < 0.15)
            continue;

        // momentumf1 = ptf1 / TMath::Sin(2.0 * TMath::ATan(exp(-2.0 * etaf1)));
        momentumf1 = ptf1 / TMath::Sin(2.0 * TMath::ATan(exp(etaf1)));
        energyf1 = sqrt((momentumf1 * momentumf1) + (massf1 * massf1));

        F1.SetPtEtaPhiE(ptf1, etaf1, phif1, energyf1); // f1 generated

        if (TMath::Abs(F1.Rapidity()) >= 0.5)
            continue;       // checking rapidity cut of f1
        hf1gen->Fill(ptf1); // counting no. of generated f1

        Double_t masses4[2] = {0.497, 0.497}; // masses of decaying particle
        event.SetDecay(F1, 2, masses4);
        weight = event.Generate();
        KSh1 = event.GetDecay(0);
        KSh2 = event.GetDecay(1);
        ksh1 = *KSh1; // kshort
        ksh2 = *KSh2; // kshort

        kshpt1 = ksh1.Pt();
        kshpt2 = ksh2.Pt();

        // kinematics cut for daughters and pairs*********************
        if (kshpt1 < 0.15 || kshpt2 < 0.15)
            continue;

        // We have not applied eta cut on K0s, but applied only rapidity cut on it. The eta cut is applied on K0s daughters
        // if (TMath::Abs(ksh1.Eta()) > 0.8)
        //     continue;
        // if (TMath::Abs(ksh2.Eta()) > 0.8)
        //     continue;

        if (TMath::Abs(ksh1.Rapidity()) > 0.5)
            continue;
        if (TMath::Abs(ksh2.Rapidity()) > 0.5)
            continue;

        // secondary decay: K0s -> pi+ pi-
        Double_t pionMasses[2] = {0.13957, 0.13957};
        if (!eventPi1.SetDecay(ksh1, 2, pionMasses))
            continue;
        if (!eventPi2.SetDecay(ksh2, 2, pionMasses))
            continue;

        eventPi1.Generate();
        eventPi2.Generate();

        piPlus1 = eventPi1.GetDecay(0);
        piMinus1 = eventPi1.GetDecay(1);
        piPlus2 = eventPi2.GetDecay(0);
        piMinus2 = eventPi2.GetDecay(1);

        pionPlus1 = *piPlus1;
        pionMinus1 = *piMinus1;
        pionPlus2 = *piPlus2;
        pionMinus2 = *piMinus2;

        if (TMath::Abs(pionPlus1.Eta()) > 0.8)
            continue;
        if (TMath::Abs(pionMinus1.Eta()) > 0.8)
            continue;
        if (TMath::Abs(pionPlus2.Eta()) > 0.8)
            continue;
        if (TMath::Abs(pionMinus2.Eta()) > 0.8)
            continue;

        //****************************************

        // calculation of efficiency of ksh, kaon and pions*********************

        TAxis *xAxisksh1 = hgenksh->GetXaxis();
        int xBinNumberksh1 = xAxisksh1->FindBin(kshpt1);
        // TAxis *yAxisksh1 = hgenksh->GetYaxis();
        // int yBinNumberksh1 = yAxisksh1->FindBin(ksh1.Eta());

        TAxis *xAxisksh2 = hgenksh->GetXaxis();
        int xBinNumberksh2 = xAxisksh2->FindBin(kshpt2);
        // TAxis *yAxisksh2 = hgenksh->GetYaxis();
        // int yBinNumberksh2 = yAxisksh2->FindBin(ksh2.Eta());

        // Commented since eta information is not available
        // genkshcontent1 = hgenksh->GetBinContent(xBinNumberksh1, yBinNumberksh1);
        // genkshcontent2 = hgenksh->GetBinContent(xBinNumberksh2, yBinNumberksh2);

        // reckshcontent1 = hrecksh->GetBinContent(xBinNumberksh1, yBinNumberksh1);
        // reckshcontent2 = hrecksh->GetBinContent(xBinNumberksh2, yBinNumberksh2);

        genkshcontent1 = hgenksh->GetBinContent(xBinNumberksh1);
        genkshcontent2 = hgenksh->GetBinContent(xBinNumberksh2);

        reckshcontent1 = hrecksh->GetBinContent(xBinNumberksh1);
        reckshcontent2 = hrecksh->GetBinContent(xBinNumberksh2);

        if (genkshcontent1 == 0 || genkshcontent2 == 0)
            continue;

        effksh1 = reckshcontent1 / genkshcontent1;
        effksh2 = reckshcontent2 / genkshcontent2;

        //*******************************************************************8

        //***accepting from random numbers****************************
        probacceptkshort1 = randomGenerator->Uniform(0.0, 1.0);
        probacceptkshort2 = randomGenerator->Uniform(0.0, 1.0);

        // cout<<"prbs are:"<<probacceptpi<<" "<<probacceptpineg<<endl;
        if (probacceptkshort1 > effksh1 || probacceptkshort2 > effksh2)
            continue;

        finaleff = effksh1 * effksh2;
        hf1rec->Fill(ptf1); // final number of reconstructed f1
        // GTE->Fill(ptf1, finaleff, 1);
    }

    // TCanvas *cGTE = new TCanvas("", "GTE", 720, 720);
    // GTE->Draw();
    TCanvas *cgen = new TCanvas("", "Generated f2", 720, 720);
    SetCanvasStyle(cgen, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hf1gen);
    hf1gen->Draw();
    TCanvas *cRec = new TCanvas("", "Reconstructed f2", 720, 720);
    SetCanvasStyle(cRec, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hf1rec);
    hf1rec->Draw();

    TCanvas *cEff = new TCanvas("", "Efficiency", 720, 720);
    SetCanvasStyle(cEff, 0.15, 0.05, 0.05, 0.15);
    TH1D *hEff = (TH1D *)hf1rec->Clone("hEff_f2");
    hEff->Divide(hf1gen);
    SetHistoQA(hEff);
    hEff->GetYaxis()->SetTitle("Acceptance x Efficiency");
    hEff->Draw();

    // TFile *f2 = new TFile("efficiencyPtSmear.root", "UPDATE");
    // hEff->Write();
    // f2->Close();
}



