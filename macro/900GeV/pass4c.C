

/////////////////////////////ALL IN ONE MACRO USING MIXED EVENT COMBINATORIAL BKG///////////////////////////////

#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "style.h"
#include "fitfunc.h"
using namespace std;

//*************************************************************************************

void pass4c()

{
    TString ResBkg = "LIKE"; // residual bkg to be used  //LIKE or MIX
    // const char* mixls[2] = {"LIKE", "MIX"};

    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);

    const Int_t Npt = 13; // no. of pt bins(max)

    //** For projection of required signals in different pt bins used as array elements**************************
    double Low_pt[Npt] =  {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0};
    double High_pt[Npt] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 8.0};

    double x[Npt] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.25, 2.75, 5.5};

    double xerr[Npt] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 2.5};

    //*****************************************************************************************************************

    //**Total no. of pt bins as an array and defining different histograms for extracting fitting parameters********************************************************************************************************************
    double NPT[Npt + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 8.0};

    TH1F *hYieldpar = new TH1F("hYieldpar", "hYieldpar", Npt, NPT);                   // for ptspectra from fitting parameter directly
    TH1F *hintegral_yield = new TH1F("hintegral_yield", "hintegral_yield", Npt, NPT); // pt spectra from function integration
    TH1F *hChiSquare = new TH1F("hChiSquare", "hChiSquare", Npt, NPT);                // for chisquare
    TH1F *hsgnfcance = new TH1F("hsgnfcance", "hsgnfcance", Npt, NPT);                // for significance of signal
    TH1F *hmass = new TH1F("hmass", "hmass", Npt, NPT);                               // for mass from fit
    TH1F *hwidth = new TH1F("hwidth", "hwidth", Npt, NPT);                            // for width from fit
    TH1F *hYbincount = new TH1F("hYbincount", "hYbincount", Npt, NPT);                // Yield calculation using bin counting method
    TH1F *hFrac_stat_error = new TH1F("hFrac_stat_error", "hFrac_stat_error", Npt, NPT);
    TH1F *herrormass = new TH1F("herrormass", "", Npt, 0.5, 16);   // for error band in mass
    TH1F *herrorwidth = new TH1F("herrorwidth", "", Npt, 0.5, 16); // for error band in width

    ////Setting bin content for mass and width error/////////////////////////////////////////////////
    for (int i = 1; i <= Npt; i++)
    {
        herrormass->SetBinContent(i, 0.895);
        herrormass->SetBinError(i, 0.0015);
        herrorwidth->SetBinContent(i, 0.047);
        herrorwidth->SetBinError(i, 0.005);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////

    //*************************************************************************************************************

    //**For normalisation **************************************************************************************

    double lownorm[Npt + 5] = {0.74, 0.74, 0.74, 0.74, 0.74, 0.74, 0.74, 1.05, 1.05, 1.05, 1.05, 1.05, 1.2};
    double highnorm[Npt + 5] = {0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 1.06, 1.06, 1.06, 1.06, 1.06, 1.3};

    //**********************************************************************************************************

    //**For fitting ranges in different pt bins******************************************************************************

    double lowfitrange[Npt + 1];
    double highfitrange[Npt + 1];

    // For ME free width (pol 3)
    double lowfitrangeme[Npt + 1] = {0.69, 0.66, 0.65, 0.68, 0.68, 0.70, 0.71, 0.73, 0.75, 0.74, 0.74, 0.74};

    double highfitrangeme[Npt + 1] = {1.15, 1.18, 1.20, 1.15, 1.12, 1.03, 1.05, 1.20, 1.2, 1.14, 1.2};

    // // For LS free width (normal and tight)
    double lowfitrangels[Npt] = {0.72, 0.73, 0.74, 0.755, 0.76, 0.78, 0.77, 0.77, 0.79, 0.79, 0.795, 0.795, 0.79};

    double highfitrangels[Npt] = {1.239, 1.2, 1.2, 1.18, 1.2, 1.2, 1.1, 1.1, 1.1, 1.12, 1.132, 1.18, 1.1};

        // For LS free width (PID cuts)
    // double lowfitrangels[Npt] = {0.72, 0.73, 0.75, 0.755, 0.76, 0.78, 0.775, 0.78, 0.79, 0.795, 0.795, 0.795, 0.79};

    // double highfitrangels[Npt] = {1.239, 1.2, 1.2, 1.18, 1.2, 1.2, 1.09, 1.1, 1.16, 1.2, 1.25, 1.25, 1.25};

    

    if (ResBkg.CompareTo("MIX") == 0)
    {
        for (int i = 0; i < Npt + 1; i++)
        {
            lowfitrange[i] = lowfitrangeme[i];
            highfitrange[i] = highfitrangeme[i];
        }
    }
    else if (ResBkg.CompareTo("LIKE") == 0)
    {
        for (int i = 0; i < Npt + 1; i++)
        {
            lowfitrange[i] = lowfitrangels[i];
            highfitrange[i] = highfitrangels[i];
        }
    }

   

    //************************************************************************************************************************

    //**For different centrality in different pt bins *********************************************************************

    double Low_Centr[Npt + 1] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double High_Centr[Npt + 1] = {100, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90};

    //**********************************************************************************************************************

    //**PDG Mass and Width of required signal for significance calculation ******************************************

    float masspdg = 0.895;  // in GeV/c^2
    float widthpdg = 0.047; // in 1 sigma GeV/c^2

    //**************************************************************************

    //**Initialisation values for yield of BW and coefficients of POL2 and POL3***************************************************************

    // Double_t Para2[Npt + 5] = {1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04};

    // 0-10

    // Double_t RBPar3L[Npt + 5] = {-100, -100, -100, -100, -100, -100, -20, -100, -100, -100, -20, -20, -20, -20};
    // Double_t RBPar3H[Npt + 5] = {100, 100, 100, 100, 100, 100, 200, 100, 100, 100, 200, 200, 200, 200};

    // Double_t RBPar4L[Npt + 5] = {-100, -100, -100, -100, -100, -100, -20, -100, -100, -100, -20, -20, -20, -20};
    // Double_t RBPar4H[Npt + 5] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100};

    // Double_t RBPar5L[Npt + 5] = {-200, -200, -200, -200, -200, -200, -50, -200, -200, -200, -50, -50, -50, -50};
    // Double_t RBPar5H[Npt + 5] = {200, 200, 200, 200, 200, 200, 30, 200, 200, 200, 30, 30, 30, 30};

    //****************************************************************************************************************************************

    //**Variables and arrays needed *******************************************************************************************

    Double_t integralsignalfunc[Npt];                         // for calculation of area under signal only (after fitting and extraction BW signal)
    Double_t ptcenter[Npt];                                   // stores mean pt for each pt bin
    Double_t ptbinwidth[Npt];                                 // stores pt binwidth
    Double_t interror[100];                                   // stores error while area calculation of signal for yield
    const int N = 1;                                          // for centrality loop below used
    Double_t dy = 1.0;                                        // for rapidity difference
    TFitResultPtr r;                                          // for fitting using TFitter
    TH1D *fHistTotal[Npt];                                    // for sig+bg
    TH1D *fHistBkg[Npt];                                      // for mixedbg
    TH1D *fHistbkgLS[Npt];                                    // for like sign
    TH1D *fHistlike[Npt];                                     // for resultant like bg
    TH1D *histPP[Npt];                                        // for like pp bg
    TH1D *histMM[Npt];                                        // for like mm bg
    TH1D *hfbkg;                                              // for normalised mixed bkg
    TH1D *hfbkg_like;                                         // for resultant like bkg
    char name[100];                                           // for giving name to canvases and histograms
    Double_t lowpt, highpt;                                   // for extracting pt from the above defined pt array
    int lbin, hbin;                                           // corresponding bin for a given pt
    int rebin = 4;                                            // for rebinning  //=1 for no rebin
    int lownormbin1, lownormbin2, highnormbin1, highnormbin2; // normalisation bins when using mixed background
    Double_t sigbkg_integral, bkg_integral, normfactor;       // integrals and normalisation factor(scaling factor) using above normalisation bins

    // details of input root file////////////////////
    int Nbins_file = 900;
    float lowrange_file = 0.6;
    float highrange_file = 1.5;
    float binwidth_file = (highrange_file - lowrange_file) * rebin / Nbins_file;
    ///////////////////////////////////////////////

    // fiiting parameters array///////////////////////
    Double_t Yield[Npt];
    Double_t Mass[Npt];
    Double_t Width[Npt];
    Double_t poly0[Npt];
    Double_t poly1[Npt];
    Double_t poly2[Npt];
    Double_t poly3[Npt];
    Double_t ErrorMass[Npt];
    Double_t ErrorWidth[Npt];
    Double_t ErrorYield[Npt];
    Double_t Chi2Ndf[Npt];
    /////////////////////////////////////////////////

    Double_t yieldcalc, yielderror; // calculation of raw yield from function integration
    Double_t yieldcalc1[Npt], yielderror1[Npt];
    Double_t Yield_value_par, Yield_error_par; // calculation of raw yield from fitting parameter
    Double_t BR = 0.66;                        // branching ratio

    Double_t significance_den, significance_num, ratio, ratio2; // calculate signal+bkg and only signal integral respectively and then significance(ratio)
    int bmin, bmax;
    Double_t hBCError_1, bkgvalue, Integral_BW_withsigma, fYield_BinCount, YieldIntegral_BW, Yfraction_cBW, sum_tail_correction, Total_Ybincounting, Tail_correction_plusm, Tail_correction_minusm, Error_2, Final_pro_error;
    Double_t nlow, nhigh;
    Double_t Yield_bincount_hist;

    // float al=0.65;
    // float bh=1.2;

    float al = 0.65;
    // float bh=1.05;
    float bh = 1.3;
    //***************************************************************************************************************************

    //**Canvas definitions and initialisations********************************************************************************

    // TCanvas *cgridls1 = new TCanvas("cgridls1", "cgridls1", 1200, 1200);
    // cgridls1->Divide(4, 4);
    // TCanvas *cgridls2 = new TCanvas("cgridls2", "cgridls2", 1200, 1200);
    // cgridls2->Divide(4, 4);
    // TCanvas *cgridme1 = new TCanvas("cgridme1", "cgridme1", 1200, 1200);
    // cgridme1->Divide(4, 4);
    // TCanvas *cgridme2 = new TCanvas("cgridme2", "cgridme2", 1200, 1200);
    // cgridme2->Divide(4, 4);

    // TCanvas *cbkgls1 = new TCanvas("cbkgls1", "cbkgls1", 1200, 1200);
    // cbkgls1->Divide(4, 4);
    // TCanvas *cbkgls2 = new TCanvas("cbkgls2", "cbkgls2", 1200, 1200);
    // cbkgls2->Divide(4, 4);
    // TCanvas *cbkgme1 = new TCanvas("cbkgme1", "cbkgme1", 1200, 1200);
    // cbkgme1->Divide(4, 4);
    // TCanvas *cbkgme2 = new TCanvas("cbkgme2", "cbkgme2", 1200, 1200);
    // cbkgme2->Divide(4, 4);

    TCanvas *cinv[Npt]; // for output canvases on screen containing fitted signal after subtraction

    for (Int_t ip = 0; ip < Npt; ip++)
    {
        TString cName = TString::Format("cinv_pt_%2.1f-%2.1f", Low_pt[ip], High_pt[ip]);
        cinv[ip] = new TCanvas(Form("cinv%d", ip), cName.Data(), 10, 10, 1200, 800);
        SetCanvasStyle2(cinv[ip], 0.15, 0.05, 0.08, 0.13);
    }

    TCanvas *cSigbkg[Npt]; // for output canvases on screen containing signal with bkg(after norm in case of mix)

    for (Int_t ip = 0; ip < Npt; ip++)
    {
        TString cNam = TString::Format("cSigbkg_pt_%2.1f-%2.1f", Low_pt[ip], High_pt[ip]);
        cSigbkg[ip] = new TCanvas(Form("cSigbkg%d", ip), cNam.Data(), 10, 10, 1200, 800);
        SetCanvasStyle2(cSigbkg[ip], 0.15, 0.05, 0.08, 0.13);
    }

    //**********************************************************************************************************************************

    //**Histogram arrays initialisation*************************************************************************************************

    for (Int_t m = 0; m < Npt; m++)
    {
        sprintf(name, "fHistNum%d", m);
        fHistTotal[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
        sprintf(name, "fHistbkg%d", m);
        fHistBkg[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
        sprintf(name, "histPP_pt%d", m);
        histPP[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
        sprintf(name, "histMM_pt%d", m);
        histMM[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
        sprintf(name, "histlike_pt%d", m);
        fHistlike[m] = new TH1D(name, "inv_mass", Nbins_file, lowrange_file, highrange_file);
    }

    //***************************************************************************************************************************************

    //**For accessing the histograms stored within root file*****************************************************************************

    TFile *fInputFile = new TFile("../data/LHC22c_pass4/AnalysisResults.root", "Read"); // Reading the File
    // TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(3);  // for lf-analysis
    TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(4);  // for lf-analysis-tight
    // TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(5);  // for pid cuts
    TList *fInputList = (TList *)fInputFile->Get(key->GetName());

    // Reading the file for yield of published pp 13 TeV results.

    // TFile *file2 = new TFile("../HEPData-ins1797443-v1-root.root", "READ");
    // TKey *key2 = (TKey *)file2->GetListOfKeys()->At(3);
    // TGraph *grsourav = (TGraph *)file2->Get("Table 4/Graph1D_y1");
    // grsourav->SetTitle(0);
    // grsourav->GetXaxis()->SetTitle(0);
    // grsourav->GetYaxis()->SetTitle(0);

    // const string kOutputName = "lf-k892analysis";
    const string kOutputName = "lf-k892analysis_tight";  // for tighter cuts
    // const string kOutputName = "lf-k892analysis_pid";  // for pic cuts
    const string kOutputName2 = "QAbefore";
    const string kOutputName3 = "QAafter";

    //**To calculate total number of events for which histograms were filled*************************************************************
    // TH1F* hEVent = (TH1F *) fInputList->FindObject("hAEventsVsMulti");
    // Double_t Event=hEVent->Integral(lc,hc);
    Double_t Event = 3.146236e7; // In the root file it is hCounterAcc
    //*************************************************************************************************************************

    cout << "*****************number of events********************:" << Event << endl;

    //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************
    TH3F *fHistNum = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassDS", kOutputName.c_str()));
    TH3F *fHistDen = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassME", kOutputName.c_str()));
    TH3F *fHistLS = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassLS", kOutputName.c_str()));

    // cout << " THE NUMBER OF BINS IN THE HISTOGRAM IS " << fHistNum->GetNbinsZ()<<endl;

    // // QA before

    // TH1F *histtrkpT_pi = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_pi", kOutputName.c_str(), kOutputName2.c_str()));
    // TH1F *histtrkpT_ka = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_ka", kOutputName.c_str(), kOutputName2.c_str()));
    // TH2F *histTofTpc = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_TPC_Map_pi_all", kOutputName.c_str(), kOutputName2.c_str()));
    // TH2F *histTofpi = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_pi_all", kOutputName.c_str(), kOutputName2.c_str()));
    // TH2F *histTpcpi = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigma_pi_all", kOutputName.c_str(), kOutputName2.c_str()));
    // TH2F *histTofka = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_ka_all", kOutputName.c_str(), kOutputName2.c_str()));
    // TH2F *histTpcka = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigmaka_all", kOutputName.c_str(), kOutputName2.c_str()));

    // // QA after

    // TH1F *histtrkpT_pi2 = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_pi", kOutputName.c_str(), kOutputName3.c_str()));
    // TH1F *histtrkpT_ka2 = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_ka", kOutputName.c_str(), kOutputName3.c_str()));
    // TH2F *histTofTpc2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_TPC_Map_pi_all", kOutputName.c_str(), kOutputName3.c_str()));
    // TH2F *histTofpi2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_pi_all", kOutputName.c_str(), kOutputName3.c_str()));
    // TH2F *histTpcpi2 = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigma_pi_all", kOutputName.c_str(), kOutputName3.c_str()));
    // TH2F *histTofka2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_ka_all", kOutputName.c_str(), kOutputName3.c_str()));
    // TH2F *histTpcka2 = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigmaka_all", kOutputName.c_str(), kOutputName3.c_str()));

    //**********************************************************************************************************************************

    gstyle(); // this is not gStyle, it is defined in the header file style.h

    //*********************************

    for (Int_t ip = 0; ip < Npt; ip++) // start pt bin loop
    // for (Int_t ip = 0; ip < 10; ip++)
    {
        lowpt = Low_pt[ip];
        highpt = High_pt[ip];

        fHistTotal[ip] = fHistNum->ProjectionZ(Form("hSig_%d", ip), -1, -1, fHistNum->GetYaxis()->FindBin(lowpt), fHistNum->GetYaxis()->FindBin(highpt), "E");
        fHistBkg[ip] = fHistDen->ProjectionZ(Form("hbkg_%d", ip), -1, -1, fHistDen->GetYaxis()->FindBin(lowpt), fHistDen->GetYaxis()->FindBin(highpt), "E");
        fHistbkgLS[ip] = fHistLS->ProjectionZ(Form("hbkgLS_%d", ip), -1, -1, fHistLS->GetYaxis()->FindBin(lowpt), fHistLS->GetYaxis()->FindBin(highpt), "E");

        // cout<<"The value of bins in z (i.e Energy bins) is "<<fHistNum->GetNbinsZ()<<endl;

        fHistTotal[ip]->Rebin(rebin);
        fHistBkg[ip]->Rebin(rebin);
        fHistbkgLS[ip]->Rebin(rebin);

        //**Cloning sig+bkg histogram for like sign or mixed event subtraction *********************************************************
        TH1D *hfsig = (TH1D *)fHistTotal[ip]->Clone();
        //*****************************************************************************************************************************

        //**If mixed bkg is used *******************************************************************************

        if (ResBkg.CompareTo("MIX") == 0)
        {

            Double_t init = lownorm[ip];
            Double_t fin = highnorm[ip];

            TAxis *axis1 = fHistTotal[ip]->GetXaxis();
            lownormbin1 = axis1->FindBin(init);
            highnormbin1 = axis1->FindBin(fin);
            sigbkg_integral = (fHistTotal[ip]->Integral(lownormbin1, highnormbin1));

            TAxis *axis2 = fHistBkg[ip]->GetXaxis();
            lownormbin2 = axis2->FindBin(init);
            highnormbin2 = axis2->FindBin(fin);
            bkg_integral = (fHistBkg[ip]->Integral(lownormbin2, highnormbin2));

            normfactor = sigbkg_integral / bkg_integral; // scaling factor for mixed bkg
            hfbkg = (TH1D *)fHistBkg[ip]->Clone();

            hfbkg->Scale(normfactor);

            hfsig->Add(hfbkg, -1);
        }
        else if (ResBkg.CompareTo("MIX") != 0)
        {
            hfbkg_like = (TH1D *)fHistbkgLS[ip]->Clone();
            hfsig->Add(hfbkg_like, -1);
        }

        //****pt bincenter and pt binwidth**********************************************************************

        ptcenter[ip] = (Low_pt[ip] + High_pt[ip]) / 2;
        ptbinwidth[ip] = High_pt[ip] - Low_pt[ip];

        //******************************************************************************************************

        TF1 *fitFcn = new TF1("fitfunc", BreitWignerpoly2, lowfitrange[ip], highfitrange[ip], 7); // sig+bkg fit function
        TF1 *fitFcn1 = new TF1("fitfunc1", polynomial3, lowfitrange[ip], highfitrange[ip], 4);    // only residualbkg
        TF1 *fitFcn2 = new TF1("fitFcn2", BW, lowfitrange[ip], highfitrange[ip], 3);              // only signal

        fitFcn->SetParLimits(0, 0.87, 0.95); // Mass
        fitFcn->SetParLimits(2, 0, 100000);  // Yield
        fitFcn->SetParameter(1, 0.047);      // width
        // fitFcn->SetParLimits(1, 0.035, 0.0590);      // width
        fitFcn->FixParameter(1, 0.047);      // width
        fitFcn->SetParameter(2, 1000); // yield

        // fitFcn->SetParameter(3, 100);
        // fitFcn->SetParameter(4, 100);
        // fitFcn->SetParameter(5, 100);
        // fitFcn->SetParameter(6, 100);

        // fitFcn->SetParLimits(3, -1.5, 1.5);
        // fitFcn->SetParLimits(4, 0, 1000);
        // fitFcn->SetParLimits(5, 0, 1000);
        // fitFcn->SetParLimits(6, -7, 7);

        /*
        fitFcn->SetParLimits(3,-200,200);
        fitFcn->SetParLimits(4,-200,200);
        fitFcn->SetParLimits(5,-200,200);
        */

        fitFcn->SetParNames("Mass","Width","Yield","A","B","C","D");
        // fitFcn->SetParNames("Mass", "Width", "Yield", "A", "B", "C");
        // fitFcn->SetParNames("Yield","Mass","Width","C","B","A","D");
        r = hfsig->Fit(fitFcn, "RS+"); // signal after bkg subtraction
        // gStyle->SetOptFit(0);

        //****************************************************************************************************************

        //**Extraction of fitting parameters******************************************************************************

        Double_t *par = fitFcn->GetParameters();

        Mass[ip] = fitFcn->GetParameter(0);
        Width[ip] = fitFcn->GetParameter(1);
        Yield[ip] = fitFcn->GetParameter(2);
        /*poly0[ip]=fitFcn->GetParameter(3);
        poly1[ip]=fitFcn->GetParameter(4);
        poly2[ip]=fitFcn->GetParameter(5);*/
        // poly2[ip] = fitFcn->GetParameter(3);
        // poly1[ip] = fitFcn->GetParameter(4);
        // poly0[ip] = fitFcn->GetParameter(5);
        // poly3[ip] = fitFcn->GetParameter(6);

        fitFcn2->SetParameters(&par[0]);
        fitFcn1->SetParameters(&par[3]);

        ErrorMass[ip] = fitFcn->GetParError(0);
        ErrorWidth[ip] = fitFcn->GetParError(1);
        ErrorYield[ip] = fitFcn->GetParError(2);
        Chi2Ndf[ip] = (fitFcn->GetChisquare()) / (fitFcn->GetNDF());

        //******************************************************************************************************************

        //**ERROR BIN COUNTING METHOD CALCULATION*****************************************************************************

        TF1 *fitFcn2_plusm = new TF1("fitFcn2_plusm", BW, lowfitrange[ip], highfitrange[ip], 3);
        TF1 *fitFcn2_minusm = new TF1("fitFcn2_minusm", BW, lowfitrange[ip], highfitrange[ip], 3);
        fitFcn2_plusm->FixParameter(0, Mass[ip] + ErrorMass[ip]);
        fitFcn2_plusm->FixParameter(1, 0.047);
        fitFcn2_plusm->FixParameter(2, Yield[ip]);

        fitFcn2_minusm->FixParameter(0, Mass[ip] - ErrorMass[ip]);
        fitFcn2_minusm->FixParameter(1, 0.047);
        fitFcn2_minusm->FixParameter(2, Yield[ip]);

        //*********************************************************************************************************************

        //**Calculation of significance and storing chi2 and sig in respective histograms*****************************************************

        bmin = hfsig->GetXaxis()->FindBin(masspdg - 2 * widthpdg);
        bmax = hfsig->GetXaxis()->FindBin(masspdg + 2 * widthpdg);

        significance_num = (fitFcn2->Integral(masspdg - 2 * widthpdg, masspdg + 2 * widthpdg)) / (binwidth_file);
        significance_den = TMath::Sqrt(fHistTotal[ip]->Integral(bmin, bmax));

        ratio = significance_num / significance_den; // significance of signal

        hsgnfcance->SetBinContent(ip + 1, ratio);
        hChiSquare->SetBinContent(ip + 1, Chi2Ndf[ip]); // storing both significance and chi2 in histogram

        //*****************************************************************************************************************************************

        //**Calculation of Yield using bin counting method and storing it in histogram***********************************************************

        Yield_bincount_hist = hfsig->IntegralAndError(bmin, bmax, hBCError_1);
        bkgvalue = fitFcn1->Integral(hfsig->GetBinLowEdge(bmin), hfsig->GetBinLowEdge(bmax + 1));
        Integral_BW_withsigma = fitFcn2->Integral(hfsig->GetBinLowEdge(bmin), hfsig->GetBinLowEdge(bmax + 1));
        fYield_BinCount = Yield_bincount_hist - (bkgvalue / binwidth_file);
        YieldIntegral_BW = fitFcn2->Integral(0.635, 5) / binwidth_file;
        Yfraction_cBW = (Integral_BW_withsigma / YieldIntegral_BW);

        sum_tail_correction = (fitFcn2->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + fitFcn2->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / binwidth_file;

        nlow = (fitFcn2->Integral(0.635, hfsig->GetBinLowEdge(bmin))) / binwidth_file;
        nhigh = (fitFcn2->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / binwidth_file;
        nlow = nlow / (Event * ptbinwidth[ip] * dy * BR * 2 * ptcenter[ip]);
        nhigh = nhigh / (Event * ptbinwidth[ip] * dy * BR * 2 * ptcenter[ip]);

        Total_Ybincounting = (sum_tail_correction + fYield_BinCount) / (Event * ptbinwidth[ip] * dy * BR * 2 * ptcenter[ip]);

        cout << "***************************************************************" << endl;
        cout << "****fraction of nlow for bin***********:"
             << " " << ip << " " << nlow / Total_Ybincounting << endl;
        cout << "****fraction of nhigh for bin***********:"
             << " " << ip << " " << nhigh / Total_Ybincounting << endl;
        cout << "***************************************************************" << endl;
        Tail_correction_plusm = (fitFcn2_plusm->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + (fitFcn2_plusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5))) / binwidth_file;
        Tail_correction_minusm = ((fitFcn2_minusm->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + fitFcn2_minusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / binwidth_file);
        Error_2 = sum_tail_correction - Tail_correction_plusm;
        Final_pro_error = TMath::Sqrt(Error_2 * Error_2 + hBCError_1 * hBCError_1) / (Event * ptbinwidth[ip] * dy * BR * 2 * ptcenter[ip]);

        ////Uncorrected Yield/////////////////////////////////////////////////////////////////////////////////

        hYbincount->SetBinContent(ip + 1, Total_Ybincounting);
        hYbincount->SetBinError(ip + 1, Final_pro_error);
        cout << "--------Total Value from bin counting----------" << (sum_tail_correction + fYield_BinCount) << endl;
        cout << "--------Value from bin counting----------" << Total_Ybincounting << endl;

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        // Fractional stat error///////////////////////////////////////////////////////////////////////////////

        hFrac_stat_error->SetBinContent(ip + 1, Final_pro_error / Total_Ybincounting);
        cout << "--------Frac error from bin counting----------" << (Final_pro_error / Total_Ybincounting) << endl;
        //////////////////////////////////////////////////////////////////////////////////////////////////////

        //****************************************************************************************************************************************

        //**Calculation for raw pt spectra using function integration and filling it in histogram*********************************************

        integralsignalfunc[ip] = (fitFcn2->Integral((masspdg - 5 * widthpdg), (masspdg + 5 * widthpdg)));
        // cout<<"The value of integral is "<<integralsignalfunc[ip]<<endl;
        TMatrixDSym cov = r->GetCovarianceMatrix();
        TMatrixDSym cov1;
        TMatrixDSym cov2;
        cov.GetSub(0, 2, 0, 2, cov1);
        cov.GetSub(3, 5, 3, 5, cov2);
        Double_t *b = cov1.GetMatrixArray();
        Double_t *a = cov2.GetMatrixArray();
        Double_t *para = fitFcn->GetParameters();
        interror[ip] = fitFcn2->IntegralError((masspdg - 5 * widthpdg), (masspdg + 5 * widthpdg), &para[0], b);

        yieldcalc = integralsignalfunc[ip] / (Event * ptbinwidth[ip] * dy * BR * /*2*/ binwidth_file /* * ptcenter[ip]*/); // raw yield calculation
        yieldcalc1[ip] = integralsignalfunc[ip] / (Event * ptbinwidth[ip] * dy * BR * /*2*/ binwidth_file);                // raw yield calculation
        // cout<<"THE VALUE OF YIELD FOR PT BIN FOR LOW PT"<<Low_pt[ip]<< " is "<<yieldcalc1[ip]<<endl;
        yielderror = interror[ip] / (Event * ptbinwidth[ip] * dy * BR * /*2*/ binwidth_file /* * ptcenter[ip] */); // raw yield error
        yielderror1[ip] = interror[ip] / (Event * ptbinwidth[ip] * dy * BR * /*2*/ binwidth_file);                 // raw yield error

        /*if (ip==10 || ip==11 || ip==12)
          {
          yieldcalc=Total_Ybincounting;
          yielderror=Final_pro_error;
          }
        */
        hintegral_yield->SetBinContent(ip + 1, yieldcalc);
        hintegral_yield->SetBinError(ip + 1, yielderror); // filling histogram including error

        cout << "Total Yield from function integration:" << integralsignalfunc[ip] << " error : " << interror[ip] << " " << ip << endl;
        cout << "Yield from function integration:" << yieldcalc << " error : " << interror[ip] << " " << ip << endl;
        cout << "--------Frac error from integral method----------" << (yielderror / yieldcalc) << endl;

        //*********************************************************************************************************************

        //**Calculation for raw pt spectra using fitting parameter and filling it in histogram*********************************************

        Yield_value_par = Yield[ip] / (Event * ptbinwidth[ip] * dy * BR * 2 * binwidth_file * ptcenter[ip]);
        Yield_error_par = ErrorYield[ip] / (Event * ptbinwidth[ip] * dy * BR * 2 * binwidth_file * ptcenter[ip]);

        hYieldpar->SetBinContent(ip + 1, Yield_value_par);
        hYieldpar->SetBinError(ip + 1, Yield_error_par);

        cout << "Total Yield from fitting:" << Yield[ip] << endl;
        cout << "Yield from fitting:" << Yield_value_par << endl;
        //***********************************************************************************************************************************

        //**Filling mass and width fitting parameter in histogram*******************************************************************************

        hmass->SetBinContent(ip + 1, Mass[ip]);
        hmass->SetBinError(ip + 1, ErrorMass[ip]);

        hwidth->SetBinContent(ip + 1, Width[ip]);
        hwidth->SetBinError(ip + 1, ErrorWidth[ip]);

        //*****************************************************************************************************************************

        //**Setting plot parameters style*************************************************************************************************

        // SetHistoStyle(hfsig,1,20,0.6,0.04,0.04,0.03,0.03,1.0,1.4);
        // SetHistoStyle(fHistTotal[ip],1,24,0.5,0.04,0.04,0.03,0.03,1.0,1.4);
        SetHistoStyle(hfsig, 1, 20, 1.0, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
        SetHistoStyle(fHistTotal[ip], 1, 24, 1.0, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

        hfsig->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
        hfsig->GetYaxis()->SetMaxDigits(2);
        hfsig->GetYaxis()->SetTitle("Counts");

        if (ResBkg.CompareTo("MIX") == 0)
        {

            SetHistoStyle(hfbkg, 2, 24, 0.5, 0.05, 0.05, 0.03, 0.03, 1.0, 1.4);
            SetHistoStyle(hfbkg, kRed, 24, 1.0, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

            hfbkg->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
            hfsig->GetYaxis()->SetMaxDigits(2);
            hfbkg->GetYaxis()->SetTitle("Counts");
        }

        else if (ResBkg.CompareTo("LIKE") == 0)
        {
            // SetHistoStyle(hfbkg_like,2,24,0.5,0.05,0.05,1.0,1.1);
            SetHistoStyle(hfbkg_like, 2, 24, 0.5, 0.05, 0.05, 0.03, 0.03, 1.0, 1.4);
            hfbkg_like->GetXaxis()->SetTitle("M_{K#pi} (GeV/#it{c}^{2})");
            hfsig->GetYaxis()->SetMaxDigits(2);
            hfbkg_like->GetYaxis()->SetTitle("Counts");
        }

        fitFcn1->SetLineColor(4);
        fitFcn1->SetLineStyle(2);
        fitFcn1->SetLineWidth(2);
        fitFcn2->SetLineColor(6);
        fitFcn2->SetLineStyle(2);
        fitFcn2->SetLineWidth(2);

        //*******************************************************************************************************************************

        //**Plot of histograms and graphs*********************************************************************************************

        cinv[ip]->cd();
        TLegend *pag = new TLegend(0.18, 0.76, 0.32, 0.9);
        hfsig->SetMinimum(0);
        hfsig->Draw("e");
        fitFcn->Draw("same");
        fitFcn1->Draw("same");
        fitFcn2->Draw("same");

        pag->AddEntry(fitFcn, "BW+pol3");
        pag->AddEntry(fitFcn1, "BW");
        pag->AddEntry(fitFcn2, "pol3");
        pag->Draw();
        t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
        t2->DrawLatex(0.65, 0.97,
                      Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
                           High_pt[ip]));
        hfsig->GetXaxis()->SetRangeUser(al, bh);
        cinv[ip]->SaveAs(Form("../output/LHC22c_pass4/hfitsig_pt%d.png", ip + 1));

        // if (ip < 16)
        // {
        //     cgridls1->cd(ip + 1);
        //     hfsig->Draw("e");
        //     fitFcn->Draw("same");
        //     fitFcn1->Draw("same");
        //     fitFcn2->Draw("same");
        //     t2->DrawLatex(0.65, 0.97,
        //                   Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
        //                        High_pt[ip]));
        // }
        // else if (ip >= 16)
        // {
        //     cgridls2->cd(ip - 15);
        //     hfsig->Draw("e");
        //     fitFcn->Draw("same");
        //     fitFcn1->Draw("same");
        //     fitFcn2->Draw("same");
        //     t2->DrawLatex(0.65, 0.97,
        //                   Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
        //                        High_pt[ip]));
        // }

        // cgridls1->SaveAs("gridme1.png");
        // cgridls2->SaveAs("gridme2.png");

        // signal extraction before bkg subtraction////////////////////////////////////////////////////////////////

        cSigbkg[ip]->cd();

        // For mixed event bkg///////////////////////////////////////////////////
        if (ResBkg.CompareTo("MIX") == 0)
        {
            hfbkg->GetXaxis()->SetRangeUser(al, bh);
            fHistTotal[ip]->Draw("e");
            fHistTotal[ip]->GetYaxis()->SetTitle("Counts");
            hfbkg->Draw("same");
            sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);

            TLegend *leg112 = new TLegend(0.680554, 0.8112735, 0.892902, 0.9138954, NULL, "brNDC");
            SetLegendStyle(leg112, 0.05, 1);
            leg112->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
            leg112->AddEntry(hfbkg, "Mixed-event bkg", "p");
            leg112->Draw();

            sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
            TLatex *ltx = new TLatex(0.35, 0.93, name);
            ltx->SetNDC();
            ltx->SetTextFont(22);
            ltx->SetTextSize(0.05);
            ltx->Draw();

            Double_t max = hfbkg->GetMaximum();
            Double_t min = fHistTotal[ip]->GetMaximum();
            /*if (max>min)
            hfbkg->GetYaxis()->SetRangeUser(0,max);
            else
            hfbkg->GetYaxis()->SetRangeUser(0,min);
            */
            cSigbkg[ip]->SaveAs(Form("../output/LHC22c_pass4/hsigbkg_pt%d.png", ip + 1));

            // if (ip < 16)
            // {
            //     cbkgme1->cd(ip + 1);
            //     fHistTotal[ip]->Draw("e");
            //     hfbkg->Draw("same");
            //     t2->DrawLatex(0.65, 0.97,
            //                   Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
            //                        High_pt[ip]));
            // }
            // else if (ip >= 16)
            // {
            //     cbkgme2->cd(ip - 15);
            //     fHistTotal[ip]->Draw("e");
            //     hfbkg->Draw("same");
            //     t2->DrawLatex(0.65, 0.97,
            //                   Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
            //                        High_pt[ip]));
            // }

            // cbkgme1->SaveAs("bkgme1.png");
            // cbkgme2->SaveAs("bkgme2.png");
        }

        ////////////////////////////////////////////////////////////////////////

        // For like sign bkg///////////////////////////////////////////////////
        else if (ResBkg.CompareTo("LIKE") == 0)
        {
            // hfbkg_like->GetXaxis()->SetRangeUser(al, bh);
            fHistTotal[ip]->Draw("e");
            fHistTotal[ip]->GetYaxis()->SetTitle("Counts");
            hfbkg_like->Draw("same");
            TLegend *leg11 = new TLegend(0.680554, 0.8112735, 0.892902, 0.9138954, NULL, "brNDC");
            SetLegendStyle(leg11, 0.04, 1);
            leg11->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
            leg11->AddEntry(hfbkg_like, "LikeSignPairs", "p");
            leg11->Draw();
            sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", Low_pt[ip], High_pt[ip]);
            TLatex *ltx = new TLatex(0.35, 0.93, name);
            ltx->SetNDC();
            ltx->SetTextFont(22);
            ltx->SetTextSize(0.05);
            ltx->Draw();

            cSigbkg[ip]->SaveAs(Form("../output/LHC22c_pass4/hsigbkg_pt%d.png", ip + 1));

            // if (ip < 16)
            // {
            //     cbkgls1->cd(ip + 1);
            //     fHistTotal[ip]->Draw("e");
            //     hfbkg_like->Draw("same");
            //     t2->DrawLatex(0.65, 0.97,
            //                   Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
            //                        High_pt[ip]));
            // }
            // else if (ip >= 16)
            // {
            //     cbkgls2->cd(ip - 15);
            //     fHistTotal[ip]->Draw("e");
            //     hfbkg_like->Draw("same");
            //     t2->DrawLatex(0.65, 0.97,
            //                   Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
            //                        High_pt[ip]));
            // }

            // cbkgls1->SaveAs("bkgls1.png");
            // cbkgls2->SaveAs("bkgls2.png");
        }

        ////////////////////////////////////////////////////////////////////////

    } // pt loop ends

    TCanvas *sig = new TCanvas ();
    SetCanvasStyle2(sig, 0.2, 0.05, 0.08, 0.2);
    hsgnfcance->Draw();
    sig->SaveAs("../output/LHC22c_pass4/significance.png");


    // TFile *filecmp = new TFile("filels.root", "RECREATE");

    TCanvas *cgr = new TCanvas("cgr", "graphs", 900, 800);
    SetCanvasStyle2(cgr, 0.2, 0.05, 0.08, 0.2);

    // chisquare_NDF vs pt
    TGraph *grchi = new TGraph(Npt, x, Chi2Ndf);
    grchi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    grchi->GetYaxis()->SetTitle("#chi^{2}/NDF ");
    SetGraphStyle(grchi, 1, 1);
    grchi->Draw("AP");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    cgr->SaveAs("../output/LHC22c_pass4/chindfvspt.png");
    // grchi->Write("chils");

    // mass vs pt
    TGraphErrors *grmass = new TGraphErrors(Npt, x, Mass, 0, ErrorMass);
    grmass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    grmass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    SetGrapherrorStyle(grmass, 1, 1);
    grmass->Draw("ap");
    grmass->Write("massls");
    TLegend *massleg = new TLegend(0.65, 0.45, 0.9, 0.65);
    massleg->SetFillColor(0);
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    TLine *line = new TLine(grmass->GetXaxis()->GetXmin(), 0.895, grmass->GetXaxis()->GetXmax(), 0.895);
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->SetLineWidth(3);
    line->Draw();
    massleg->AddEntry(line, "PDG Mass", "l");
    massleg->Draw("l");
    cgr->SaveAs("../output/LHC22c_pass4/mass_pt.png");

    // Width vs pT

    // TGraphErrors *grwidth = new TGraphErrors(Npt, x, Width, 0, ErrorWidth);
    // grwidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // grwidth->GetYaxis()->SetTitle("Width (GeV)");
    // SetGrapherrorStyle(grwidth, 1, 1);
    // grwidth->Draw("ap");
    // grwidth->Write("widthls");
    // TLegend *widthleg = new TLegend(0.45, 0.75, 0.7, 0.85);
    // widthleg->SetFillColor(0);
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    // TLine *line2 = new TLine(grmass->GetXaxis()->GetXmin(), 0.047, grmass->GetXaxis()->GetXmax(), 0.047);
    // line2->SetLineStyle(2);
    // line2->SetLineColor(2);
    // line2->SetLineWidth(3);
    // line2->Draw();
    // widthleg->AddEntry(line, "PDG Width", "l");
    // widthleg->Draw("l");
    // cgr->SaveAs("../output/LHC22c_pass4/width_pt.png");

    // yield vs pt

    TCanvas *cyield1 = new TCanvas("cyield1", "yeild", 1000, 800);
    gPad->SetLogy();
    SetCanvasStyle2(cyield1, 0.2, 0.05, 0.08, 0.2);
    TGraphErrors *gryield1 = new TGraphErrors(Npt, x, yieldcalc1, xerr, yielderror1);
    TLegend *legyield = new TLegend(0.5, 0.7, 0.9, 0.9);
    gryield1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gryield1->GetYaxis()->SetTitle("Yield (GeV/c)^{-1}");
    SetGrapherrorStyle(gryield1, 1, 1);
    legyield->AddEntry(gryield1, "pp 13.6 TeV ");

    gryield1->Draw("ap");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    legyield->SetTextSize(0.03);
    legyield->SetTextFont(2);
    legyield->Draw();
    cyield1->SaveAs("../output/LHC22c_pass4/yield1.png");



    // filecmp->Close();
}
