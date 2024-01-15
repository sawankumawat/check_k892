
#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/common.h"
using namespace std;

//*************************************************************************************

void kstar()

{
    const string kResBkg = "MIX";

    // const string kfoldername = "lf-k892analysis";
    // const string kfoldername = "lf-k892analysis_PID_TPC_30";
    const string kfoldername = "lf-k892analysis_PID_combined_20";
    // const string kfoldername = "lf-k892analysis_mix_10";

    const int kRebin = 3;
    const float txtsize = 0.045;

    const Int_t Npt = 21; // 12 for pbpb and 18 and 21 for pp

    //** For projection of required signals in different pt bins used as array elements**************************
    // for pp data low and high IR
    // double pT_bins[Npt + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 8.0};
    double pT_bins[Npt + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 8.0, 10, 12, 15};
    // for pbpb data LHC23zzh_pass1
    // double pT_bins[Npt + 1] = {0.4, 0.8, 1.2, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0};

    //*************************************************************************************************************

    //**For normalisation **************************************************************************************

    // dataset LHC23_pass1_lowB_highIR_sampling
    double lownorm[Npt + 15] = {1.1, 0.68, 0.68, 0.68, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
    double highnorm[Npt + 15] = {1.15, 0.73, 0.73, 0.73, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

    // dataset LHC23_pass1_lowB_lowIR
    // double lownorm[Npt + 15] = {0.68, 0.68, 0.68, 0.68, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
    // double highnorm[Npt + 15] = {0.73, 0.73, 0.73, 0.73, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

    // dataset LHC23zzh_pass1
    // double lownorm[Npt + 15] = {0.7, 0.7, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
    // double highnorm[Npt + 15] = {0.745, 0.745, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

    //**********************************************************************************************************

    //**For fitting ranges in different pt bins******************************************************************************

    //*************************MIXED EVENTS****************************************
    // // dataset LHC23_pass1_lowB_highIR_sampling
     double lowfitrangeme[Npt + 15]   =  {0.72, 0.67, 0.7, 0.8, 0.72, 0.72, 0.77, 0.74, 0.77, 0.78, 0.78, 0.7, 0.75, 0.7, 0.7, 0.72, 0.73, 0.7, 0.7, 0.73, 0.7, 0.7};
    double highfitrangeme[Npt + 15] = {1.1, 1.06, 1.1, 1.08, 1.06, 1.04, 1.06, 1.04, 1.04, 1.04, 1.02, 1.05, 1.1, 1.08, 1.1, 1.1, 1.1, 1.04, 1.06, 1.1, 1.1, 1.1, 1.1};

    // dataset LHC23_pass1_lowB_lowIR
    //   double lowfitrangeme[Npt + 15] = {0.75, 0.65, 0.78, 0.78, 0.72, 0.72, 0.73, 0.73, 0.75, 0.76, 0.76, 0.7, 0.75, 0.75, 0.73, 0.73, 0.73, 0.68, 0.8, 0.75, 0.75, 0.75};
    // double highfitrangeme[Npt + 15] = {1.1, 1.06, 1.1, 1.06, 1.04, 1.04, 1.04, 1.04, 1.04, 1.05, 1.05, 1.1, 1.1, 1.1, 1.1, 1.15, 1.12, 1.05, 1.0, 1.1, 1.1, 1.1, 1.1};

    // dataset LHC23zzh_pass1
    // double lowfitrangeme[Npt + 15] = {0.85, 0.85, 0.85, 0.8, 0.85, 0.85, 0.85, 0.82, 0.82, 0.73, 0.73, 0.7, 0.75, 0.75, 0.75, 0.73, 0.7, 0.67, 0.75, 0.75, 0.75, 0.75};
    // double highfitrangeme[Npt + 15] = {0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 1.05, 1.05, 1.1, 1.1, 1.1, 1.1, 1.15, 1.12, 1.13, 1.1, 1.1, 1.1, 1.1, 1.1};

    // ***************************LIKE SIGN************************************

    // For LS default (pol2)
    double lowfitrangels[Npt + 15] = {0.72, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.78, 0.75, 0.78, 0.74, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75};
    double highfitrangels[Npt + 15] = {1.1, 1.04, 1.1, 1.1, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1};

    double lowfitrange[Npt + 15];
    double highfitrange[Npt + 15];
    for (int i = 0; i < Npt; i++)
    {
        lowfitrange[i] = (kResBkg == "MIX") ? lowfitrangeme[i] : lowfitrangels[i];
        highfitrange[i] = (kResBkg == "MIX") ? highfitrangeme[i] : highfitrangels[i];
    }

    //************************************************************************************************************************

    //**For different centrality in different pt bins *********************************************************************

    double Low_Centr[Npt + 15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double High_Centr[Npt + 15] = {100, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90};

    //**PDG Mass and Width of required signal for significance calculation ******************************************

    float masspdg = 0.895;  // in GeV/c^2
    float widthpdg = 0.047; // in 1 sigma GeV/c^2

    //**************************************************************************

    // Histograms initialization

    TH1F *hYieldpar = new TH1F("hYieldpar", "hYieldpar", Npt, pT_bins);                   // for ptspectra from fitting parameter directly
    TH1F *hintegral_yield = new TH1F("hintegral_yield", "hintegral_yield", Npt, pT_bins); // pt spectra from function integration
    TH1F *hChiSquare = new TH1F("hChiSquare", "hChiSquare", Npt, pT_bins);                // for chisquare
    TH1F *hsignificance = new TH1F("hsignificance", "hsignificance", Npt, pT_bins);       // for significance of signal
    TH1F *hmass = new TH1F("hmass", "hmass", Npt, pT_bins);                               // for mass from fit
    TH1F *hwidth = new TH1F("hwidth", "hwidth", Npt, pT_bins);                            // for width from fit
    TH1F *hYbincount = new TH1F("hYbincount", "hYbincount", Npt, pT_bins);                // Yield calculation using bin counting method
    TH1F *hFrac_stat_error = new TH1F("hFrac_stat_error", "hFrac_stat_error", Npt, pT_bins);
    TH1F *herrormass = new TH1F("herrormass", "", Npt, 0.5, 16);   // for error band in mass
    TH1F *herrorwidth = new TH1F("herrorwidth", "", Npt, 0.5, 16); // for error band in width

    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.06);
    t2->SetTextFont(42);

    ////Setting bin content for mass and width error/////////////////////////////////////////////////
    for (int i = 1; i <= Npt; i++)
    {
        herrormass->SetBinContent(i, masspdg);
        herrormass->SetBinError(i, 0.0015);
        herrorwidth->SetBinContent(i, widthpdg);
        herrorwidth->SetBinError(i, 0.005);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////

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
    char name[100];                                           // for giving name to canvases and histograms
    Double_t lowpt, highpt;                                   // for extracting pt from the above defined pt array
    int lbin, hbin;                                           // corresponding bin for a given pt
    int lownormbin1, lownormbin2, highnormbin1, highnormbin2; // normalisation bins when using mixed background
    Double_t sigbkg_integral, bkg_integral, normfactor;       // integrals and normalisation factor(scaling factor) using above normalisation bins

    // details of input root file////////////////////
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
    float al = 0.65;
    float bh = 1.4;

    //***************************************************************************************************************************

    //**Canvas definitions and initialisations********************************************************************************

    TCanvas *cinv[Npt]; // for output canvases on screen containing fitted signal after subtraction

    for (Int_t ip = 0; ip < Npt; ip++)
    {
        TString cName = TString::Format("cinv_pt_%2.1f-%2.1f", pT_bins[ip], pT_bins[ip + 1]);
        cinv[ip] = new TCanvas(Form("cinv%d", ip), cName.Data(), 10, 10, 1200, 1100);
        SetCanvasStyle2(cinv[ip], 0.15, 0.05, 0.08, 0.13);
    }

    TCanvas *cSigbkg[Npt]; // for output canvases on screen containing signal with bkg(after norm in case of mix)

    for (Int_t ip = 0; ip < Npt; ip++)
    {
        TString cNam = TString::Format("cSigbkg_pt_%2.1f-%2.1f", pT_bins[ip], pT_bins[ip + 1]);
        cSigbkg[ip] = new TCanvas(Form("cSigbkg%d", ip), cNam.Data(), 10, 10, 1200, 1100);
        SetCanvasStyle2(cSigbkg[ip], 0.15, 0.05, 0.08, 0.13);
    }
    //**For accessing the histograms stored within root file******************************

    // TFile *fInputFile = new TFile("/home/sawan/check_k892/data/pbpb/23zzk_pass1_relval/AnalysisResults.root", "Read");
    TFile *fInputFile = new TFile(kDataFilename.c_str(), "Read");
    if (!fInputFile)
    {
        cerr << "File not found " << endl;
        return;
    }

    TH1F *hentries = (TH1F *)fInputFile->Get("event-selection-task/hColCounterAcc");
    double Event = hentries->GetEntries();
    cout << "*****************number of events********************:" << Event << endl;

    //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************

    TH3F *fHistNum = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassDS", kfoldername.c_str()));
    TH3F *fHistNum_anti = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassDSAnti", kfoldername.c_str()));
    fHistNum->Add(fHistNum_anti, 1);
    fHistNum->Scale(0.5);
    TH3F *fHistDen = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassME", kfoldername.c_str()));
    TH3F *fHistLS = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassLS", kfoldername.c_str()));
    TH3F *fHistLS_anti = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassLSAnti", kfoldername.c_str()));
    fHistLS->Add(fHistLS_anti, 1);
    fHistLS->Scale(0.5);

    // cout << " THE NUMBER OF BINS IN THE HISTOGRAM IS " << fHistNum->GetNbinsZ()<<endl;

    gstyle(); // this is not gStyle, it is defined in the header file style.h
    gStyle->SetOptFit(0);
    //***************************************************************************************************
    TCanvas *cgrid1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid_bkg1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid_bkg2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    cgrid1->Divide(kcanvasdivide[0], kcanvasdivide[1]);
    cgrid2->Divide(kcanvasdivide[0], kcanvasdivide[1]);
    cgrid_bkg1->Divide(kcanvasdivide[0], kcanvasdivide[1]);
    cgrid_bkg2->Divide(kcanvasdivide[0], kcanvasdivide[1]);

    for (Int_t ip = 0; ip < Npt; ip++) // start pt bin loop
    // for (Int_t ip = 10; ip < Npt; ip++) // start pt bin loop
    {
        lowpt = pT_bins[ip];
        highpt = pT_bins[ip + 1];

        fHistTotal[ip] = fHistNum->ProjectionZ(Form("hSig_%d", ip), -1, -1, fHistNum->GetYaxis()->FindBin(lowpt + 0.001), fHistNum->GetYaxis()->FindBin(highpt - 0.001), "E"); // ProjectionZ("title", xrangelow, xrangehigh, yrangelow, yrangehigh). range -1, -1 is used to take full range in the x axis. In the 3D histogram, x axis is centrality (V0M), y axis is pT, z axis is invariant mass
        fHistBkg[ip] = fHistDen->ProjectionZ(Form("hbkg_%d", ip), -1, -1, fHistDen->GetYaxis()->FindBin(lowpt + 0.001), fHistDen->GetYaxis()->FindBin(highpt - 0.001), "E");
        fHistbkgLS[ip] = fHistLS->ProjectionZ(Form("hbkgLS_%d", ip), -1, -1, fHistLS->GetYaxis()->FindBin(lowpt + 0.001), fHistLS->GetYaxis()->FindBin(highpt - 0.001), "E");

        auto binwidth_file = (fHistTotal[ip]->GetXaxis()->GetXmax() - fHistTotal[ip]->GetXaxis()->GetXmin()) * kRebin / fHistTotal[ip]->GetXaxis()->GetNbins();
        // cout << "no. of bins " << fHistTotal[ip]->GetXaxis()->GetNbins() << " xmin value " << fHistTotal[ip]->GetXaxis()->GetXmin() << " xmax value " << fHistTotal[ip]->GetXaxis()->GetXmax() << endl;

        //**Cloning sig+bkg histogram for like sign or mixed event subtraction *********************************************************
        TH1D *hfsig = (TH1D *)fHistTotal[ip]->Clone();
        //*****************************************************************************************************************************

        if (kResBkg == "MIX")
        {
            sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(lownorm[ip]), fHistTotal[ip]->GetXaxis()->FindBin(highnorm[ip])));
            bkg_integral = (fHistBkg[ip]->Integral(fHistBkg[ip]->GetXaxis()->FindBin(lownorm[ip]), fHistBkg[ip]->GetXaxis()->FindBin(highnorm[ip])));
            normfactor = sigbkg_integral / bkg_integral; // scaling factor for mixed bkg
            hfbkg = (TH1D *)fHistBkg[ip]->Clone();

            hfbkg->Scale(normfactor);
            hfbkg->Rebin(kRebin);
            hfsig->Rebin(kRebin);

            hfsig->Add(hfbkg, -1);
        }
        else
        {
            // sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(lownorm[ip]), fHistTotal[ip]->GetXaxis()->FindBin(highnorm[ip])));
            // bkg_integral = (fHistbkgLS[ip]->Integral(fHistbkgLS[ip]->GetXaxis()->FindBin(lownorm[ip]), fHistbkgLS[ip]->GetXaxis()->FindBin(highnorm[ip])));
            // normfactor = sigbkg_integral / bkg_integral; // scaling factor for mixed bkg
            // hfbkg = (TH1D *)fHistbkgLS[ip]->Clone();
            // hfbkg->Scale(normfactor);

            hfbkg = (TH1D *)fHistbkgLS[ip]->Clone();
            hfbkg->Rebin(kRebin);
            hfsig->Rebin(kRebin);
            hfsig->Add(hfbkg, -1);
        }

        fHistTotal[ip]->Rebin(kRebin);

        //****pt bincenter and pt binwidth**********************************************************************

        ptcenter[ip] = (pT_bins[ip] + pT_bins[ip + 1]) / 2;
        ptbinwidth[ip] = pT_bins[ip + 1] - pT_bins[ip];

        //*****************************************************************************************************

        TF1 *fitFcn = new TF1("fitfunc", BreitWignerpoly2, lowfitrange[ip], highfitrange[ip], 6); // sig+bkg fit function
        TF1 *fitFcn1 = new TF1("fitfunc1", polynomial2, lowfitrange[ip], highfitrange[ip], 3);    // only residualbkg
        TF1 *fitFcn2 = new TF1("fitFcn2", BW, lowfitrange[ip], highfitrange[ip], 3);              // only signal

        fitFcn->SetParLimits(0, 0.80, 0.98); // Mass
        fitFcn->SetParameter(2, 30000);      // yield
        fitFcn->SetParLimits(2, 0, 10e9);    // Yield
        fitFcn->SetParameter(1, 0.047);      // width
        // fitFcn->FixParameter(1, 0.047);      // width
        // fitFcn->SetParLimits(1, 0, 0.120); // width

        // fitFcn->SetParNames("Mass","Width","Yield","A","B","C","D");
        fitFcn->SetParNames("Mass", "Width", "Yield", "C", "B", "A");
        r = hfsig->Fit(fitFcn, "REBMS0+"); // signal after bkg subtraction

        //****************************************************************************************************************

        //**Extraction of fitting parameters******************************************************************************

        Double_t *par = fitFcn->GetParameters();

        Mass[ip] = fitFcn->GetParameter(0);
        Width[ip] = fitFcn->GetParameter(1);
        Yield[ip] = fitFcn->GetParameter(2);
        poly2[ip] = fitFcn->GetParameter(3);
        poly1[ip] = fitFcn->GetParameter(4);
        poly0[ip] = fitFcn->GetParameter(5);
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

        significance_den = TMath::Sqrt(fHistTotal[ip]->Integral(bmin, bmax));
        significance_num = (fitFcn2->Integral(masspdg - 2 * widthpdg, masspdg + 2 * widthpdg)) / (binwidth_file);

        ratio = significance_num / significance_den; // significance of signal

        hsignificance->SetBinContent(ip + 1, ratio);
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

        // cout << "***************************************************************" << endl;
        // cout << "****fraction of nlow for bin***********:"
        //      << " " << ip << " " << nlow / Total_Ybincounting << endl;
        // cout << "****fraction of nhigh for bin***********:"
        //      << " " << ip << " " << nhigh / Total_Ybincounting << endl;
        // cout << "***************************************************************" << endl;
        Tail_correction_plusm = (fitFcn2_plusm->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + (fitFcn2_plusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5))) / binwidth_file;
        Tail_correction_minusm = ((fitFcn2_minusm->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + fitFcn2_minusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / binwidth_file);
        Error_2 = sum_tail_correction - Tail_correction_plusm;
        Final_pro_error = TMath::Sqrt(Error_2 * Error_2 + hBCError_1 * hBCError_1) / (Event * ptbinwidth[ip] * dy * BR * 2 * ptcenter[ip]);

        ////Uncorrected Yield/////////////////////////////////////////////////////////////////////////////////

        hYbincount->SetBinContent(ip + 1, Total_Ybincounting);
        hYbincount->SetBinError(ip + 1, Final_pro_error);
        // cout << "--------Total Value from bin counting----------" << (sum_tail_correction + fYield_BinCount) << endl;
        // cout << "--------Value from bin counting----------" << Total_Ybincounting << endl;

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        // Fractional stat error///////////////////////////////////////////////////////////////////////////////

        hFrac_stat_error->SetBinContent(ip + 1, Final_pro_error / Total_Ybincounting);
        // cout << "--------Frac error from bin counting----------" << (Final_pro_error / Total_Ybincounting) << endl;
        //////////////////////////////////////////////////////////////////////////////////////////////////////

        //****************************************************************************************************************************************

        //**Calculation for raw pt spectra using function integration and filling it in histogram*********************************************

        integralsignalfunc[ip] = (fitFcn2->Integral((masspdg - 5 * widthpdg), (masspdg + 5 * widthpdg)));
        TMatrixDSym cov = r->GetCovarianceMatrix();
        TMatrixDSym cov1;
        TMatrixDSym cov2;
        cov.GetSub(0, 2, 0, 2, cov1);
        cov.GetSub(3, 5, 3, 5, cov2);
        Double_t *b = cov1.GetMatrixArray();
        Double_t *a = cov2.GetMatrixArray();
        Double_t *para = fitFcn->GetParameters();
        interror[ip] = fitFcn2->IntegralError((masspdg - 5 * widthpdg), (masspdg + 5 * widthpdg), &para[0], b);

        yieldcalc = integralsignalfunc[ip] / (Event * ptbinwidth[ip] * dy * BR * binwidth_file); // raw yield calculation
        yielderror = interror[ip] / (Event * ptbinwidth[ip] * dy * BR * binwidth_file);          // raw yield error

        hintegral_yield->SetBinContent(ip + 1, yieldcalc);
        hintegral_yield->SetBinError(ip + 1, yielderror); // filling histogram including error

        //**Filling mass and width fitting parameter in histogram*******************************************************************************

        hmass->SetBinContent(ip + 1, Mass[ip]);
        hmass->SetBinError(ip + 1, ErrorMass[ip]);

        hwidth->SetBinContent(ip + 1, Width[ip]);
        hwidth->SetBinError(ip + 1, ErrorWidth[ip]);

        //*****************************************************************************************************************************

        //**Setting plot parameters style*************************************************************************************************

        SetHistoStyle(hfsig, 1, 20, 1, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
        SetHistoStyle(fHistTotal[ip], 1, 8, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

        hfsig->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
        hfsig->GetYaxis()->SetMaxDigits(2);
        hfsig->GetYaxis()->SetTitle("Counts");

        SetHistoStyle(hfbkg, kRed, 24, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);
        hfbkg->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
        // hfsig->GetYaxis()->SetMaxDigits(2);
        hfbkg->GetYaxis()->SetTitle("Counts");

        fitFcn1->SetLineColor(4);
        fitFcn1->SetLineStyle(2);
        fitFcn1->SetLineWidth(4);
        fitFcn2->SetLineColor(6);
        fitFcn2->SetLineStyle(2);
        fitFcn2->SetLineWidth(4);
        fitFcn->SetLineWidth(4);

        //*******************************************************************************************************************************

        //**Plot of histograms and graphs*********************************************************************************************
        auto chibyndf = fitFcn->GetChisquare() / fitFcn->GetNDF();

        // inv mass histograms after the background subraction
        (multipanel_plots == 1) ? (ip < kupperpad * klowerpad) ? cgrid1->cd(ip + 1) : cgrid2->cd(ip + 1 - kupperpad * klowerpad) : cinv[ip]->cd();
        gPad->SetRightMargin(0.015);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        TLegend *pag = new TLegend(0.4, 0.7, 0.91, 0.9);
        TLegend *pag2 = new TLegend(0.2, 0.7, 0.45, 0.9);
        pag->SetBorderSize(0);
        pag->SetTextFont(42);
        pag->SetTextSize(0.04);
        pag->SetFillStyle(0);
         pag2->SetBorderSize(0);
        pag2->SetTextFont(42);
        pag2->SetTextSize(0.04);
        pag2->SetFillStyle(0);
        // if (ip == 0)
        // {
        //     hfsig->SetMaximum(hfsig->GetMaximum() * 1.8);
        // }
        // else if (ip == 1)
        //     hfsig->SetMaximum(hfsig->GetMaximum() * 1);
        // else if (ip == 2)
        //     hfsig->SetMaximum(hfsig->GetMaximum() * 0.9);

        // else
            hfsig->SetMaximum(hfsig->GetMaximum() * 1.5);
        hfsig->SetMarkerSize(0.5);
        fitFcn->SetLineWidth(2);
        fitFcn1->SetLineWidth(2);
        fitFcn2->SetLineWidth(2);
        hfsig->GetXaxis()->SetRangeUser(0.645, 1.145);
        hfsig->Draw("e");
        fitFcn->Draw("same");
        fitFcn1->Draw("same");
        fitFcn2->Draw("same");
        pag2->AddEntry(fitFcn, "BW+pol2");
        pag2->AddEntry(fitFcn1, "BW");
        pag2->AddEntry(fitFcn2, "pol2");
        pag->AddEntry((TObject *)0, Form("Mass: %.3f #pm %.3f", Mass[ip], ErrorMass[ip]), "");
        pag->AddEntry((TObject *)0, Form("Width: %.3f #pm %.3f", Width[ip], ErrorWidth[ip]), "");
        pag->AddEntry((TObject *)0, Form("Yield: %.1e #pm %.1e", yieldcalc * Event, yielderror * Event), "");
        pag->AddEntry((TObject *)0, Form("Probability: %f ", fitFcn->GetProb()), "");
        pag->AddEntry((TObject *)0, Form("#chi^{2}/NDF: %.2f ", chibyndf), "");
        pag->Draw();
        pag2->Draw();
        // t2->DrawLatex(0.26, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
        t2->DrawLatex(0.27, 0.95,
                      Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip],
                           pT_bins[ip + 1]));
        if (multipanel_plots == 0 && save_plots == 1)
            cinv[ip]->SaveAs(Form((kSignalOutput + "/hfitsig_pt%d.png").c_str(), ip + 1));
        // cinv[ip]->Close();

        // inv distribution before the background subtraction
        (multipanel_plots == 1) ? (ip < klowerpad * kupperpad) ? cgrid_bkg1->cd(ip + 1) : cgrid_bkg2->cd(ip + 1 - klowerpad * kupperpad) : cSigbkg[ip]->cd();
        TH1F *hbkg_nopeak = (TH1F *)hfbkg->Clone();
        hbkg_nopeak->SetLineColor(kRed);
        hbkg_nopeak->SetMarkerColor(kRed);
        hbkg_nopeak->SetFillColor(kRed);
        hbkg_nopeak->SetFillStyle(3001);
        for (int i = 0; i < hbkg_nopeak->GetNbinsX(); i++)
        {
            if (hbkg_nopeak->GetBinCenter(i + 1) < lownorm[ip] || hbkg_nopeak->GetBinCenter(i + 1) > highnorm[ip])
            {
                hbkg_nopeak->SetBinContent(i + 1, -999);
            }
        }
        gPad->SetRightMargin(0.01);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        fHistTotal[ip]->SetMarkerSize(1);
        fHistTotal[ip]->SetMaximum(fHistTotal[ip]->GetMaximum() * 1.3);
        fHistTotal[ip]->GetXaxis()->SetRangeUser(0.59, 1.31);
        fHistTotal[ip]->Draw("E");
        fHistTotal[ip]->GetYaxis()->SetTitle("Counts");
        TLegend *leg112 = new TLegend(0.60554, 0.7812735, 0.852902, 0.8938954, NULL, "brNDC");
        leg112->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
        SetLegendStyle(leg112);
        leg112->SetTextSize(txtsize);
        sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", pT_bins[ip], pT_bins[ip + 1]);
        TLatex *ltx = new TLatex(0.27, 0.95, name);
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.06);
        hfbkg->SetMarkerSize(1);
        hfbkg->Draw("E same");
        if (kResBkg == "MIX")
            hbkg_nopeak->Draw("BAR same");
        (kResBkg == "MIX") ? leg112->AddEntry(hfbkg, "Mixed-event bkg", "p") : leg112->AddEntry(hfbkg, "Like sign pairs", "p");
        leg112->Draw();
        ltx->Draw();
        if (multipanel_plots == 0 && save_plots == 1)
            cSigbkg[ip]->SaveAs(Form((kSignalOutput + "/hsigbkg_pt%d.png").c_str(), ip + 1));
        // cSigbkg[ip]->Close();

        ////////////////////////////////////////////////////////////////////////

    } // pt loop ends
    if (multipanel_plots == 1 && save_plots == 1)
    {
        cgrid1->SaveAs((kSignalOutput + "/grid1.png").c_str());
        cgrid2->SaveAs((kSignalOutput + "/grid2.png").c_str());
        cgrid_bkg1->SaveAs((kSignalOutput + "/grid_bkg1.png").c_str());
        cgrid_bkg2->SaveAs((kSignalOutput + "/grid_bkg2.png").c_str());
    }

    TFile *filecmp = new TFile((kSignalOutput + "/pp_lowIR.root").c_str(), "RECREATE");

    TCanvas *csig = new TCanvas("", "", 1200, 1000);
    SetCanvasStyle2(csig, 0.18, 0.05, 0.08, 0.15);
    // SetHistoStyle(hsignificance, 1, 20, 1, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
    // hsignificance->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // hsignificance->Draw();
    // hsignificance->Write("significance");
    // csig->SaveAs((kSignalOutput + "/significance.png").c_str());
    // csig->Clear();

    // // // chisquare_NDF vs pt

    // hChiSquare->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // hChiSquare->GetYaxis()->SetTitle("#chi^{2}/NDF ");
    // SetHistoQA(hChiSquare);
    // hChiSquare->Draw("Pe");
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    // csig->SaveAs((kSignalOutput + "/chi.png").c_str());
    // // hChiSquare->Write("chils");
    // csig->Clear();

    // // mass vs pt
    hmass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hmass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    // SetHistoQA(hmass);
    SetHistoStyle(hmass, 1, 20, 1, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
    hmass->SetMaximum(0.91);
    hmass->Draw("pe");
    hmass->Write("mass");
    TLegend *massleg = new TLegend(0.65, 0.2, 0.9, 0.3);
    SetLegendStyle(massleg);
    massleg->SetTextSize(txtsize);
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    TLine *line = new TLine(hmass->GetXaxis()->GetXmin(), 0.895, hmass->GetXaxis()->GetXmax(), 0.895);
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->SetLineWidth(3);
    line->Draw();
    massleg->AddEntry(line, "PDG Mass", "l");
    massleg->Draw("l");
    csig->SaveAs((kSignalOutput + "/mass.png").c_str());
    csig->Clear();

    // // // Width vs pT
    hwidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hwidth->GetYaxis()->SetTitle("Width (GeV)");
    // SetHistoQA(hwidth);
    SetHistoStyle(hwidth, 1, 20, 1, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
    hwidth->SetMaximum(hwidth->GetMaximum() * 2);
    hwidth->SetMinimum(0);
    hwidth->Draw("pe");
    hwidth->Write("width");
    TLegend *widthleg = new TLegend(0.2, 0.75, 0.4, 0.85);
    SetLegendStyle(widthleg);
    widthleg->SetTextSize(txtsize);
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    TLine *line2 = new TLine(hwidth->GetXaxis()->GetXmin(), 0.047, hwidth->GetXaxis()->GetXmax(), 0.047);
    line2->SetLineStyle(2);
    line2->SetLineColor(2);
    line2->SetLineWidth(3);
    line2->Draw();
    widthleg->AddEntry(line, "PDG Width", "l");
    widthleg->SetFillStyle(0);
    widthleg->Draw();
    csig->SaveAs((kSignalOutput + "/width_pt.png").c_str());
    csig->Clear();

    // // // Yield vs pT
    SetHistoStyle(hintegral_yield, 1, 20, 1, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
    hintegral_yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hintegral_yield->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    gPad->SetLogy();
    hintegral_yield->GetXaxis()->SetRangeUser(-0.1, 15.2);
    hintegral_yield->Draw("pe");
    hintegral_yield->Write("yield");
    TLegend *legyield = new TLegend(0.8, 0.8, 0.91, 0.9);
    SetLegendStyle(legyield);
    legyield->SetTextSize(txtsize);
    // legyield->AddEntry(hYieldpar, "pbpb 5.36 TeV");
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    legyield->Draw();
    csig->SaveAs((kSignalOutput + "/yield.png").c_str());
}
