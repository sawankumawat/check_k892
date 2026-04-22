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

Double_t voigt(Double_t *x, Double_t *par)
{
    return (par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]));
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
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void toy_model_ptSmear_prottayPiKp()
{

    // accepting histograms for K0s, kaon and pion eff---------------------------------------
    gStyle->SetOptStat(0);

    // TFile *fInputFile = TFile::Open("K0sEff.root");
    // TFile *fInputFile = TFile::Open("PiKpGenRec.root");
    TFile *fInputFile = TFile::Open("AnalysisMCpass1pass2.root");

    // TH1D *hgenksh = (TH1D *)fInputFile->Get("K0s_generated_pT");
    // TH1D *hrecksh = (TH1D *)fInputFile->Get("K0s_reconstructed_pT");

    // TH1D *hgenpip = (TH1D *)fInputFile->Get("piplus_generated");
    // TH1D *hgenpim = (TH1D *)fInputFile->Get("piminus_generated");
    // TH1D *hgenkp = (TH1D *)fInputFile->Get("Kaon_generated");
    // TH1D *hgenkm = (TH1D *)fInputFile->Get("Kaon_minus_generated");

    // TH1D *hrecpip = (TH1D *)fInputFile->Get("piplus_reconstructed");
    // TH1D *hrecpim = (TH1D *)fInputFile->Get("piminus_reconstructed");
    // TH1D *hreckp = (TH1D *)fInputFile->Get("Kaon_reconstructed");
    // TH1D *hreckm = (TH1D *)fInputFile->Get("Kaon_minus_reconstructed");

    // For other file
    string folderName = "RsnOut_f1_defaultfinal_2.0_3.0_4.0_0.03_15.000_0_0_1_1.0_0.06_0.0_0.000_0.0_0.00";
    TList *list = (TList *)fInputFile->Get(folderName.c_str());
    if (list == nullptr)
    {
        cout << "Error: list/folder not found: " << folderName << endl;
        cout << "Available top-level keys in file:" << endl;
        TIter nextkey(fInputFile->GetListOfKeys());
        TKey *key;
        while ((key = (TKey *)nextkey()))
            cout << "  " << key->GetName() << endl;
        return;
    }

    auto GetHistFromList = [&](const char *baseName) -> TH2F *
    {
        TObject *obj = list->FindObject(baseName);
        if (!obj)
            obj = list->FindObject(Form("%s/%s", folderName.c_str(), baseName));
        return dynamic_cast<TH2F *>(obj);
    };

    TH2F *h2Dgenpip = GetHistFromList("KStarPlusMinusppMC_Pionposgenerated");
    TH2F *h2Dgenpim = GetHistFromList("KStarPlusMinusppMC_Pionneggenerated");
    TH2F *h2Dgenkp = GetHistFromList("KStarPlusMinusppMC_Kposgenerated");
    TH2F *h2Dgenkm = GetHistFromList("KStarPlusMinusppMC_Kneggenerated");

    TH2F *h2Drecpip = GetHistFromList("KStarPlusMinusppMC_Pionposreconstructed");
    TH2F *h2Drecpim = GetHistFromList("KStarPlusMinusppMC_Pionnegreconstructed");
    TH2F *h2Dreckp = GetHistFromList("KStarPlusMinusppMC_Kposreconstructed");
    TH2F *h2Dreckm = GetHistFromList("KStarPlusMinusppMC_Knegreconstructed");

    if (h2Dgenpip == nullptr || h2Dgenpim == nullptr || h2Dgenkp == nullptr || h2Dgenkm == nullptr ||
        h2Drecpip == nullptr || h2Drecpim == nullptr || h2Dreckp == nullptr || h2Dreckm == nullptr)
    {
        cout << "Error: One of the histograms not found! Missing:" << endl;
        if (!h2Dgenpip)
            cout << "  KStarPlusMinusppMC_Pionposgenerated" << endl;
        if (!h2Dgenpim)
            cout << "  KStarPlusMinusppMC_Pionneggenerated" << endl;
        if (!h2Dgenkp)
            cout << "  KStarPlusMinusppMC_Kposgenerated" << endl;
        if (!h2Dgenkm)
            cout << "  KStarPlusMinusppMC_Kneggenerated" << endl;
        if (!h2Drecpip)
            cout << "  KStarPlusMinusppMC_Pionposreconstructed" << endl;
        if (!h2Drecpim)
            cout << "  KStarPlusMinusppMC_Pionnegreconstructed" << endl;
        if (!h2Dreckp)
            cout << "  KStarPlusMinusppMC_Kposreconstructed" << endl;
        if (!h2Dreckm)
            cout << "  KStarPlusMinusppMC_Knegreconstructed" << endl;

        cout << "Objects available inside list " << folderName << ":" << endl;
        TIter next(list);
        TObject *obj = nullptr;
        while ((obj = next()))
            cout << "  " << obj->GetName() << " (" << obj->ClassName() << ")" << endl;
        return;
    }

    TH1D *hgenpip = h2Dgenpip->ProjectionX("hgenpip");
    TH1D *hgenpim = h2Dgenpim->ProjectionX("hgenpim");
    TH1D *hgenkp = h2Dgenkp->ProjectionX("hgenkp");
    TH1D *hgenkm = h2Dgenkm->ProjectionX("hgenkm");

    TH1D *hrecpip = h2Drecpip->ProjectionX("hrecpip");
    TH1D *hrecpim = h2Drecpim->ProjectionX("hrecpim");
    TH1D *hreckp = h2Dreckp->ProjectionX("hreckp");
    TH1D *hreckm = h2Dreckm->ProjectionX("hreckm");
    //----------------------------------------------------------------------

    double NPT[] = {0.439, 0.547, 0.646, 0.735, 0.843, 0.942, 1.05, 1.14, 1.24, 1.34, 1.43, 1.54, 1.63, 1.75, 1.94, 2.09, 2.29, 2.49, 2.68, 2.9, 3.24, 3.77, 4.25, 4.77, 5.52, 7.03}; // For Phi 13 TeV
    // double NPT[] = {0.0874, 0.197, 0.328, 0.437, 0.546, 0.656, 0.743, 0.831, 0.94, 1.09, 1.27, 1.51, 1.7, 1.9, 2.21, 2.62, 2.99, 3.39, 3.78, 4.5, 5.51, 6.49, 7.52, 9.01, 11, 13.5}; // For K* 13 TeV
    // double NPT[7] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0}; // pt bins of f1 in data

    int sizeNPT = sizeof(NPT) / sizeof(NPT[0]);

    // histograms for efficiency caluclation of f1-------------------------
    // TProfile *GTE = new TProfile("GTE", "GTE", 6, NPT);
    TH1D *hf1gen = new TH1D("f1gen", "f2(1525) Gen", sizeNPT - 1, NPT);
    TH1D *hf1rec = new TH1D("f1rec", "f2(1525) Rec", sizeNPT - 1, NPT);
    // TH1D *hkshkaonmass = new TH1D("kshkaonmass", "kshkaonmass", 300, 0, 3);
    hf1gen->Sumw2();
    hf1rec->Sumw2();

    //////////////////////////////////////////////////////////

    // TLorentz Vectors**********************************
    TLorentzVector F1;
    TLorentzVector *KaonPos;
    TLorentzVector *KaonNeg;

    TLorentzVector kaonpos;
    TLorentzVector kaonneg;
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

    // //***Breit Wigner and gaussian distributions generated*********************
    // // for K*
    // TF1 *fE = new TF1("fE", BreitWigner, 0.65, 1.5, 3);
    // fE->FixParameter(0, 0.895);
    // fE->FixParameter(1, 0.047);
    // fE->FixParameter(2.0, 1); // amplitude

    // For phi meson
    TF1 *fE = new TF1("fE", voigt, 1.0, 2.0, 3);
    fE->FixParameter(0, 1.0);      // amplitude
    fE->FixParameter(1, 1.01946);  // mass
    fE->FixParameter(2, 0.0012);   // Gaussian width (detector resolution)
    fE->FixParameter(3, 0.004249); // Lorentzian width (resonance width)

    //***************************************

    //*****variables***********************************
    double ptf1, etaf1, phif1, momentumf1, energyf1;
    TRandom *eventGenerator = new TRandom();
    TRandom3 *randomGenerator = new TRandom3();
    double massf1 = 0.0;
    double probacceptkaonPos = 0.0, probacceptkaonNeg = 0.0;
    double probacceptpionPos = 0.0, probacceptpionNeg = 0.0;
    double kaonptpos = 0.0, kaonptneg = 0.0;
    double effkaonpos = 0.0, effkaonneg = 0.0, finaleff = 0.0;
    double effpionpos = 0.0, effpionneg = 0.0;
    double genkaoncontentpos = 0.0, genkaoncontentneg = 0.0;
    double genpioncontentpos = 0.0, genpioncontentneg = 0.0;
    double reckaoncontentpos = 0.0, reckaoncontentneg = 0.0;
    double reckpioncontentpos = 0.0, reckpioncontentneg = 0.0;

    // std::vector<double> kshKaonMassValues;
    // double kshKaonMass = 0.0;
    ///////////////////////////////////////////////////////

    for (int i = 0; i < 1000000; i++)
    {
        massf1 = fE->GetRandom();            // mass of f1 from breitwigner
        ptf1 = gRandom->Uniform(0.0, 10.0);  // pt
        etaf1 = gRandom->Uniform(-0.8, 0.8); // eta
        phif1 = gRandom->Uniform(0, 6.28);   // phi

        if (ptf1 < 0.1)
            continue;

        // momentumf1 = ptf1 / TMath::Sin(2.0 * TMath::ATan(exp(-2.0 * etaf1)));
        momentumf1 = ptf1 / TMath::Sin(2.0 * TMath::ATan(exp(etaf1)));
        energyf1 = sqrt((momentumf1 * momentumf1) + (massf1 * massf1));

        F1.SetPtEtaPhiE(ptf1, etaf1, phif1, energyf1); // f1 generated

        if (TMath::Abs(F1.Rapidity()) >= 0.5)
            continue;       // checking rapidity cut of f1
        hf1gen->Fill(ptf1); // counting no. of generated f1

        Double_t masses4[2] = {0.49367, 0.49367}; // KK (for phi)
        // Double_t masses4[2] = {0.49367, 0.13957}; // Kpi (for k892)
        event.SetDecay(F1, 2, masses4);
        weight = event.Generate();
        KaonPos = event.GetDecay(0);
        KaonNeg = event.GetDecay(1);
        kaonpos = *KaonPos;
        kaonneg = *KaonNeg;

        kaonptpos = kaonpos.Pt();
        kaonptneg = kaonneg.Pt();

        // kinematics cut for daughters and pairs*********************
        if (kaonptpos < 0.1 || kaonptneg < 0.1)
            continue;

        // We have not applied eta cut on K0s, but applied only rapidity cut on it. The eta cut is applied on K0s daughters
        if (TMath::Abs(kaonpos.Eta()) > 0.8)
            continue;
        if (TMath::Abs(kaonneg.Eta()) > 0.8)
            continue;

        //****************************************

        // calculation of efficiency of ksh, kaon and pions*********************

        // // For K*(892)
        // TAxis *xAxiskaonPos = hgenkp->GetXaxis();
        // int xBinNumberkaonPos = xAxiskaonPos->FindBin(kaonptpos);

        // TAxis *xAxiskaonNeg = hgenpim->GetXaxis();
        // int xBinNumberkaonNeg = xAxiskaonNeg->FindBin(kaonptneg);

        // TAxis *xAxispionPos = hgenpip->GetXaxis();
        // int xBinNumberpionPos = xAxispionPos->FindBin(kaonptpos);

        // TAxis *xAxispionNeg = hgenpim->GetXaxis();
        // int xBinNumberpionNeg = xAxispionNeg->FindBin(kaonptneg);

        // genkaoncontentpos = hgenkp->GetBinContent(xBinNumberkaonPos);
        // genkaoncontentneg = hgenkm->GetBinContent(xBinNumberkaonNeg);
        // genpioncontentpos = hgenpip->GetBinContent(xBinNumberpionPos);
        // genpioncontentneg = hgenpim->GetBinContent(xBinNumberpionNeg);

        // reckaoncontentpos = hreckp->GetBinContent(xBinNumberkaonPos);
        // reckaoncontentneg = hreckm->GetBinContent(xBinNumberkaonNeg);
        // reckpioncontentpos = hrecpip->GetBinContent(xBinNumberpionPos);
        // reckpioncontentneg = hrecpim->GetBinContent(xBinNumberpionNeg);

        // if (genkaoncontentpos == 0 || genkaoncontentneg == 0 || genpioncontentpos == 0 || genpioncontentneg == 0)
        //     continue;

        // effkaonpos = reckaoncontentpos / genkaoncontentpos;
        // effkaonneg = reckaoncontentneg / genkaoncontentneg;
        // effpionpos = reckpioncontentpos / genpioncontentpos;
        // effpionneg = reckpioncontentneg / genpioncontentneg;

        //// For phi meson
        TAxis *xAxiskaonPos = hgenkp->GetXaxis();
        int xBinNumberkaonPos = xAxiskaonPos->FindBin(kaonptpos);

        TAxis *xAxiskaonNeg = hgenkm->GetXaxis();
        int xBinNumberkaonNeg = xAxiskaonNeg->FindBin(kaonptneg);

        genkaoncontentpos = hgenkp->GetBinContent(xBinNumberkaonPos);
        genkaoncontentneg = hgenkm->GetBinContent(xBinNumberkaonNeg);

        reckaoncontentpos = hreckp->GetBinContent(xBinNumberkaonPos);
        reckaoncontentneg = hreckm->GetBinContent(xBinNumberkaonNeg);

        if (genkaoncontentpos == 0 || genkaoncontentneg == 0)
            continue;

        effkaonpos = reckaoncontentpos / genkaoncontentpos;
        effkaonneg = reckaoncontentneg / genkaoncontentneg;

        //*******************************************************************8

        //***accepting from random numbers****************************
        probacceptkaonPos = randomGenerator->Uniform(0.0, 1.0);
        probacceptkaonNeg = randomGenerator->Uniform(0.0, 1.0);
        if (probacceptkaonPos > effkaonpos || probacceptkaonNeg > effkaonneg)
            continue;

        // // For pion only (for K*)
        // probacceptpionPos = randomGenerator->Uniform(0.0, 1.0);
        // probacceptpionNeg = randomGenerator->Uniform(0.0, 1.0);
        // if (probacceptpionPos > effpionpos || probacceptpionNeg > effpionneg)
        //     continue;

        finaleff = effkaonpos * effkaonneg; // For phi meson
        // finaleff = effkaonpos * effpionneg; // For K* meson
        hf1rec->Fill(ptf1); // final number of reconstructed f1
        // GTE->Fill(ptf1, finaleff, 1);
    }

    // // TCanvas *cGTE = new TCanvas("", "GTE", 720, 720);
    // // GTE->Draw();
    // TCanvas *cgen = new TCanvas("", "Generated f2", 720, 720);
    // SetCanvasStyle(cgen, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(hf1gen);
    // hf1gen->Draw();
    // TCanvas *cRec = new TCanvas("", "Reconstructed f2", 720, 720);
    // SetCanvasStyle(cRec, 0.15, 0.05, 0.05, 0.15);
    // SetHistoQA(hf1rec);
    // hf1rec->Draw();

    TCanvas *cEff = new TCanvas("", "Efficiency", 720, 720);
    SetCanvasStyle(cEff, 0.15, 0.05, 0.05, 0.15);
    double pad1Size, pad2Size;
    canvas_style(cEff, pad1Size, pad2Size);
    cEff->cd(1);
    TH1D *hEff = (TH1D *)hf1rec->Clone("hEff_f2");
    hEff->Divide(hf1gen);
    SetHistoQA(hEff);
    hEff->GetYaxis()->SetTitle("Acceptance x Efficiency");
    hEff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hEff->SetMaximum(0.63);
    hEff->SetMinimum(-0.03);
    hEff->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hEff->GetYaxis()->SetLabelSize(0.035 / pad1Size);
    hEff->GetXaxis()->SetLabelSize(0.035 / pad1Size);
    hEff->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hEff->GetYaxis()->SetTitleOffset(1.55 * pad1Size);
    hEff->Draw();

    double effArray[] = {0.00323, 0.0172, 0.0452, 0.0817, 0.12, 0.156, 0.185, 0.214, 0.238, 0.262, 0.283, 0.296, 0.309, 0.327, 0.346, 0.366, 0.381, 0.405, 0.414, 0.435, 0.454, 0.477, 0.488, 0.488, 0.495, 0.491}; // For phi meson
    double effArrayKstar[] = {0.03, 0.0356, 0.0441, 0.0516, 0.0628, 0.0797, 0.103, 0.132, 0.162, 0.208, 0.262, 0.31, 0.344, 0.379, 0.418, 0.455, 0.491, 0.514, 0.526, 0.541, 0.547, 0.553, 0.548, 0.548, 0.548};    // For K* meson
    TH1F *hEffMC = new TH1F("hEffMC", "Efficiency from MC", sizeNPT - 1, NPT);
    for (int i = 1; i <= sizeNPT - 1; i++)
    {
        hEffMC->SetBinContent(i, effArray[i - 1]);
        hEffMC->SetBinError(i, 0.00001);
    }
    SetHistoQA(hEffMC);
    hEffMC->SetLineColor(kRed);
    hEffMC->SetMarkerColor(kRed);
    hEffMC->SetMarkerStyle(20);
    hEffMC->Draw("pe same");

    TLegend *legend = new TLegend(0.17, 0.75, 0.5, 0.93);
    legend->AddEntry((TObject *)0, "Acc x Efficiency #Phi", "");
    // legend->AddEntry((TObject *)0, "Acc x Efficiency K*(892)", "");
    legend->AddEntry(hEffMC, "MC", "pe");
    legend->AddEntry(hEff, "Toy Model", "pe");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.05);
    legend->Draw();

    cEff->cd(2);

    TH1F *hRatio = (TH1F *)hEff->Clone("hRatio");
    hRatio->Divide(hEffMC);
    SetHistoQA(hRatio);
    hRatio->GetYaxis()->SetTitle("Toy / MC");
    hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    // hRatio->SetMaximum(0.58);
    // hRatio->SetMinimum(-0.05);
    hRatio->SetMaximum(1.38);
    hRatio->SetMinimum(0.25);
    hRatio->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    hRatio->GetYaxis()->SetLabelSize(0.035 / pad2Size);
    hRatio->GetXaxis()->SetLabelSize(0.035 / pad2Size);
    hRatio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatio->GetYaxis()->SetTitleOffset(1.4 * pad2Size);
    hRatio->Draw("pe");

    // TFile *f2 = new TFile("efficiencyPtSmear.root", "UPDATE");
    // hEff->Write();
    // f2->Close();
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.14);
    pad2->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
}
// Efficiency of Kstar from xy scan
// 0.0874	0.03
// 0.197	0.0356
// 0.328	0.0441
// 0.437	0.0516
// 0.546	0.0628
// 0.656	0.0797
// 0.743	0.103
// 0.831	0.132
// 0.94	    0.162
// 1.09	    0.208
// 1.27	    0.262
// 1.51	    0.31
// 1.7	    0.344
// 1.9	    0.379
// 2.21	    0.418
// 2.62	    0.455
// 2.99	    0.491
// 3.39	    0.514
// 3.78	    0.526
// 4.5	    0.541
// 5.51	    0.547
// 6.49	    0.553
// 7.52	    0.548
// 9.01	    0.548
// 11	    0.548
// 13.5	    0.548

//// Efficiency of Phi from xy scan
// 0.439	0.00323
// 0.547	0.0172
// 0.646	0.0452
// 0.735	0.0817
// 0.843	0.12
// 0.942	0.156
// 1.05	0.185
// 1.14	0.214
// 1.24	0.238
// 1.34	0.262
// 1.43	0.283
// 1.54	0.296
// 1.63	0.309
// 1.75	0.327
// 1.94	0.346
// 2.09	0.366
// 2.29	0.381
// 2.49	0.405
// 2.68	0.414
// 2.9	0.435
// 3.24	0.454
// 3.77	0.477
// 4.25	0.488
// 4.77	0.488
// 5.52	0.495
// 7.03	0.491
