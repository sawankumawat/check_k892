

/////////////////////////////ALL IN ONE MACRO USING MIXED EVENT COMBINATORIAL BKG///////////////////////////////

#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "src/style.h"
#include "src/fitfunc.h"
using namespace std;

//*************************************************************************************

void pass4m()

{
    TString ResBkg = "MIX"; // residual bkg to be used  //LIKE or MIX
    // const char* mixls[2] = {"LIKE", "MIX"};

    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);

    const Int_t Npt = 27; // no. of pt bins(max)

    //** For projection of required signals in different pt bins used as array elements**************************
    double Low_pt[Npt] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10, 12};
    double High_pt[Npt] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10, 12, 15};


    double x[Npt] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.1, 1.3, 1.5, 1.7, 1.9, 2.2, 2.6, 3, 3.4, 3.8, 4.5, 5.5, 6.5, 7.5, 9, 11, 13.5};


    double xerr[Npt] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 1, 1, 1.5};

    //*****************************************************************************************************************

    // //**Total no. of pt bins as an array and defining different histograms for extracting fitting parameters********************************************************************************************************************
    // double NPT[Npt+1] ={0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10, 12, 15};
    double NPT[Npt + 1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10, 12, 15};

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
    //PID2
    // double lownorm[Npt + 5]  = {0.74, 0.74, 0.74, 0.74, 0.74, 0.74, 0.74, 1.05, 1.05, 1.05, 1.1, 1.1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
    // double highnorm[Npt + 5] = {0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 1.06, 1.06, 1.06, 1.2, 1.2, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

    //PID_4, 3.5 & default & DCA xy tight
    double lownorm[Npt + 5]  = {0.74, 0.74, 0.74, 0.74, 0.74, 1.05, 1.05, 1.05, 1.05, 1.05, 1.1, 1.1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
    double highnorm[Npt + 5] = {0.75, 0.75, 0.75, 0.75, 0.75, 1.06, 1.06, 1.06, 1.06, 1.06, 1.2, 1.2, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

    // //DCA z tight and DCA all tight
    // double lownorm[Npt + 5]  = {1.05, 1.05, 1.05, 1.05, 0.74, 0.74, 1.05, 1.05, 1.05, 1.05, 1.1, 1.1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
    // double highnorm[Npt + 5] = {1.06, 1.06, 1.06, 1.06, 0.75, 0.75, 1.06, 1.06, 1.06, 1.06, 1.2, 1.2, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};



    //**********************************************************************************************************

    //**For fitting ranges in different pt bins******************************************************************************

    double lowfitrange[Npt + 1];
    double highfitrange[Npt + 1];

    //*****************MIXED EVENT***************************

    // For ME fixed width (pol 2)(PID 2)
    // double lowfitrangeme[Npt + 1] = {0.79, 0.79, 0.78, 0.79, 0.76, 0.75, 0.75, 0.75, 0.75, 0.74, 0.74, 0.74, 0.77, 0.79, 0.78, 0.78, 0.78, 0.76, 0.73, 0.74, 0.7, 0.72, 0.7, 0.7, 0.7, 0.7, 0.7};

    // double highfitrangeme[Npt + 1] = {1.06, 1.06, 1.04, 1.05, 1.08, 1.05, 1.05, 1.20, 1.2, 1.07, 1.1, 1.1, 1.1, 1.12, 1.1, 1.12, 1.1, 1.1, 1.14, 1.15, 1.15, 1.08, 1.18, 1.1, 1.08, 1.11, 1.08};

     // For ME fixed width (pol 2)(PID 4)
    // double lowfitrangeme[Npt + 1] = {0.79, 0.79, 0.78, 0.78, 0.78, 0.73, 0.75, 0.75, 0.75, 0.74, 0.74, 0.74, 0.77, 0.79, 0.78, 0.78, 0.78, 0.76, 0.73, 0.74, 0.7, 0.72, 0.7, 0.7, 0.7, 0.7, 0.7};

    // double highfitrangeme[Npt + 1] = {1.06, 1.0, 1.04, 1.02, 1.02, 1.05, 1.05, 1.05, 1.05, 1.08, 1.1, 1.1, 1.1, 1.12, 1.1, 1.12, 1.1, 1.1, 1.14, 1.15, 1.15, 1.08, 1.18, 1.1, 1.08, 1.11, 1.08};

    // For ME fixed width (pol 2)(PID 3.5)
    // double lowfitrangeme[Npt + 1] = {0.79, 0.79, 0.78, 0.78, 0.77, 0.73, 0.75, 0.75, 0.75, 0.77, 0.75, 0.74, 0.77, 0.79, 0.78, 0.78, 0.78, 0.76, 0.73, 0.74, 0.7, 0.72, 0.7, 0.7, 0.7, 0.7, 0.7};

    // double highfitrangeme[Npt + 1] = {1.06, 1.05, 1.04, 1.02, 1.05, 1.05, 1.05, 1.05, 1.05, 1.08, 1.05, 1.1, 1.1, 1.12, 1.1, 1.12, 1.1, 1.1, 1.14, 1.15, 1.15, 1.08, 1.18, 1.1, 1.08, 1.11, 1.08};


    // For ME fixed width (pol 2)(PID default)
    double lowfitrangeme[Npt + 1] = {0.79, 0.79, 0.78, 0.77, 0.77, 0.73, 0.75, 0.75, 0.75, 0.77, 0.75, 0.74, 0.77, 0.78, 0.78, 0.78, 0.78, 0.76, 0.73, 0.74, 0.7, 0.72, 0.7, 0.7, 0.7, 0.7, 0.7};

    double highfitrangeme[Npt + 1] = {1.06, 1.05, 1.04, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.08, 1.05, 1.1, 1.1, 1.1, 1.1, 1.12, 1.1, 1.1, 1.14, 1.15, 1.15, 1.08, 1.18, 1.1, 1.08, 1.11, 1.08};

    //  // For ME fixed width (pol 2)(DCA all tight)
    // double lowfitrangeme[Npt + 1] = {0.74, 0.74, 0.7, 0.7, 0.7, 0.76, 0.78, 0.73, 0.75, 0.77, 0.75, 0.74, 0.77, 0.78, 0.78, 0.78, 0.78, 0.76, 0.73, 0.74, 0.7, 0.72, 0.7, 0.7, 0.7, 0.7, 0.7};

    // double highfitrangeme[Npt + 1] = {1.0, 1.0, 1.0, 1.03, 1.09, 1.04, 1.07, 1.05, 1.05, 1.08, 1.05, 1.1, 1.1, 1.1, 1.1, 1.12, 1.1, 1.1, 1.14, 1.15, 1.15, 1.08, 1.18, 1.1, 1.08, 1.11, 1.08};

    // For ME fixed width (pol 2)(DCA xy tight)
    // double lowfitrangeme[Npt + 1] = {0.75, 0.78, 0.77, 0.74, 0.79, 0.77, 0.75, 0.73, 0.75, 0.77, 0.754, 0.74, 0.77, 0.78, 0.78, 0.78, 0.78, 0.76, 0.73, 0.74, 0.7, 0.72, 0.73, 0.7, 0.7, 0.7, 0.7};

    // double highfitrangeme[Npt + 1] = {1.01, 1.06, 1.04, 1.08, 1.09, 1.05, 1.07, 1.04, 1.05, 1.08, 1.05, 1.1, 1.1, 1.1, 1.1, 1.12, 1.1, 1.1, 1.14, 1.15, 1.15, 1.08, 1.18, 1.1, 1.08, 1.11, 1.08};

    // For ME fixed width (pol 2)(DCA z tight)
    // double lowfitrangeme[Npt + 1] = {0.75, 0.78, 0.77, 0.74, 0.79, 0.77, 0.75, 0.73, 0.75, 0.77, 0.754, 0.74, 0.77, 0.78, 0.78, 0.78, 0.78, 0.76, 0.73, 0.74, 0.7, 0.72, 0.73, 0.7, 0.7, 0.7, 0.7};

    // double highfitrangeme[Npt + 1] = {1.01, 1.06, 1.04, 1.08, 1.09, 1.05, 1.07, 1.04, 1.05, 1.08, 1.05, 1.1, 1.1, 1.1, 1.1, 1.12, 1.1, 1.1, 1.14, 1.15, 1.15, 1.08, 1.18, 1.1, 1.08, 1.11, 1.08};






    


    //***********************************LIKE SIGN*********************************


    // // For LS fixed width (pol2, Default)
    //    double lowfitrangels[Npt + 1] = {0.79, 0.79, 0.788, 0.77, 0.755, 0.74, 0.75, 0.745, 0.77, 0.77, 0.75, 0.76, 0.77, 0.754, 0.76, 0.78, 0.754, 0.75, 0.75, 0.76, 0.77, 0.765, 0.71, 0.72, 0.74, 0.74, 0.74};

    //     double highfitrangels[Npt + 1] = {1.11, 1.05, 1.04, 1.03, 1.04, 0.999, 1.049, 1.045, 1.05, 1.04, 1.06, 1.04, 1.04, 1.06, 1.05, 1.09, 1.12, 1.15, 1.12, 1.15, 1.12, 1.15, 1.18, 1.1, 1.1, 1.1, 1.1};

        // // For LS fixed width (pol3, default)
    // double lowfitrangels[Npt + 1] = {0.79, 0.79, 0.783, 0.77, 0.755, 0.74, 0.75, 0.745, 0.77, 0.77, 0.75, 0.78, 0.77, 0.754, 0.76, 0.78, 0.754, 0.75, 0.75, 0.76, 0.77, 0.765, 0.71, 0.72, 0.74, 0.74, 0.77};

    // double highfitrangels[Npt + 1] = {1.11, 1.05, 1.04, 1.04, 1.04, 1.05, 1.05, 1.045, 1.05, 1.07, 1.04, 1.04, 1.04, 1.06, 1.05, 1.09, 1.12, 1.15, 1.12, 1.15, 1.12, 1.15, 1.18, 1.1, 1.1, 1.1, 1.1};

      // // For LS fixed width (pol2, PID3.5)
    //    double lowfitrangels[Npt + 1] = {0.79, 0.79, 0.79, 0.77, 0.755, 0.74, 0.75, 0.745, 0.77, 0.77, 0.75, 0.76, 0.77, 0.754, 0.76, 0.78, 0.754, 0.75, 0.75, 0.76, 0.77, 0.765, 0.71, 0.72, 0.74, 0.74, 0.74};

    //     double highfitrangels[Npt + 1] = {1.11, 1.05, 1.038, 1.03, 1.01, 0.999, 1.049, 1.045, 1.05, 1.04, 1.06, 1.04, 1.04, 1.06, 1.05, 1.09, 1.12, 1.15, 1.12, 1.15, 1.12, 1.15, 1.18, 1.1, 1.1, 1.1, 1.1};

      // // For LS fixed width (pol2, PID4)
    //    double lowfitrangels[Npt + 1] = {0.79, 0.79, 0.79, 0.77, 0.755, 0.74, 0.75, 0.745, 0.77, 0.78, 0.75, 0.76, 0.77, 0.754, 0.76, 0.78, 0.754, 0.75, 0.75, 0.76, 0.77, 0.765, 0.71, 0.72, 0.74, 0.74, 0.74};

    //     double highfitrangels[Npt + 1] = {1.11, 1.05, 1.05, 1.04, 1.01, 0.999, 1.049, 1.045, 1.05, 1.04, 1.06, 1.04, 1.04, 1.06, 1.05, 1.09, 1.12, 1.15, 1.12, 1.15, 1.12, 1.15, 1.18, 1.1, 1.1, 1.25, 1.1};

     // // For LS fixed width (pol2, PID2)
    //    double lowfitrangels[Npt + 1] = {0.79, 0.79, 0.79, 0.77, 0.75, 0.74, 0.75, 0.745, 0.77, 0.76, 0.75, 0.76, 0.77, 0.754, 0.76, 0.78, 0.754, 0.75, 0.75, 0.76, 0.77, 0.765, 0.71, 0.72, 0.74, 0.74, 0.74};

    //     double highfitrangels[Npt + 1] = {1.06, 1.05, 1.05, 1.04, 1.01, 0.999, 1.049, 1.045, 1.05, 1.04, 1.06, 1.04, 1.04, 1.06, 1.05, 1.09, 1.12, 1.15, 1.12, 1.15, 1.12, 1.15, 1.18, 1.1, 1.1, 1.25, 1.1};

    //  // For LS fixed width (pol2, DCA z tight)
    //    double lowfitrangels[Npt + 1] = {0.78, 0.75, 0.75, 0.77, 0.77, 0.745, 0.75, 0.745, 0.77, 0.77, 0.75, 0.76, 0.77, 0.754, 0.76, 0.78, 0.754, 0.75, 0.75, 0.76, 0.77, 0.765, 0.71, 0.72, 0.74, 0.74, 0.74};

    //    double highfitrangels[Npt + 1] = {1.1, 1.01, 1.0, 1.03, 1.01, 1.02, 1.049, 1.045, 1.05, 1.04, 1.06, 1.04, 1.04, 1.06, 1.05, 1.09, 1.12, 1.15, 1.12, 1.15, 1.12, 1.15, 1.18, 1.1, 1.1, 1.1, 1.1};

    // // For LS fixed width (pol2, DCA xy tight)
    //    double lowfitrangels[Npt + 1] = {0.78, 0.77, 0.76, 0.755, 0.77, 0.745, 0.75, 0.745, 0.77, 0.77, 0.75, 0.76, 0.77, 0.754, 0.76, 0.78, 0.754, 0.75, 0.75, 0.76, 0.77, 0.765, 0.71, 0.72, 0.74, 0.74, 0.74};

    //    double highfitrangels[Npt + 1] = {1.04, 1.05, 1.01, 1.02, 1.04, 1.02, 1.049, 1.045, 1.05, 1.04, 1.06, 1.04, 1.04, 1.06, 1.05, 1.09, 1.12, 1.15, 1.12, 1.15, 1.12, 1.15, 1.18, 1.1, 1.1, 1.1, 1.1};

    // // // For LS fixed width (pol2, DCA all tight)
       double lowfitrangels[Npt + 1] = {0.76, 0.75, 0.74, 0.755, 0.77, 0.745, 0.75, 0.745, 0.77, 0.77, 0.75, 0.76, 0.77, 0.754, 0.76, 0.78, 0.754, 0.75, 0.75, 0.76, 0.77, 0.765, 0.71, 0.72, 0.74, 0.74, 0.74};

       double highfitrangels[Npt + 1] = {1.04, 1.02, 1.0, 1.02, 1.04, 1.02, 1.049, 1.045, 1.05, 1.04, 1.06, 1.04, 1.04, 1.06, 1.05, 1.09, 1.12, 1.15, 1.12, 1.15, 1.12, 1.15, 1.18, 1.1, 1.1, 1.1, 1.1};








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

    Double_t Para2[Npt + 5] = {1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04, 1.0e+04};

    // 0-10
    //for 3 sigma PID cuts (LS)
    // Double_t RBPar3L[Npt] = {-100 , -50, -30, -100, -100, -10, -10, -8, -40, -100, -8, -20, -50, -50, -50, -10, -10, -10, -50, -50, -50, -100, -50, -30, -50, -20, -20};
    // Double_t RBPar3H[Npt] = {100,    50,  30,  100,  100,  10,  10,  8,  40,  100,  8,  20,  50,  50, 50,  10,   10,  10,  50,  50,  50,  100,  50,  30,  50,  20,  20};

    //for 2.5 sigma PID cuts
    // Double_t RBPar3L[Npt] = {-10 , -50, -20, -10, -10, -100, -20, -20, -100, -10, -50000, -2000, -200, -10, -50, -10, -10, -100, -20, -50, -50, -100, -10, -50, -10, -20, -50};
    // Double_t RBPar3H[Npt] = {10,    50,  10,  20,  10,  100,  20,  20,  100,  10,  50000,  2000,  200,  10, 50,  10,   10,  100,  20,  50,  50,  100,  10,  50,  10,  20,  50};

    //for 2 sigma PID cuts
    //  Double_t RBPar3L[Npt] = {-10 , -10, -100, -10, -20, -10, -10, -150, -20, -10, -50000, -2000, -200, -100, -100, -10, -10, -30, -50, -100, -10, -10, -100, -10, -100, -20, -100};
    // Double_t RBPar3H[Npt] =  {10,    10,  10,   20,  20,  10,  10,  150,  20,  10,  50000,  2000,  200,  100, 100,    10,  10,  30,  50,  100,  10,  100,  10,  100,  50,  20,  100};

    //for DCA_all_tight cuts (ME)
    //  Double_t RBPar3L[Npt] = {-10 , -100, -100, -100, -20, -100, -100, -150, -20, -100, -500, -2000, -200, -100, -100, -100, -10, -1000, -50, -10, -10, -100, -100, -100, -100, -10, -100};
    // Double_t RBPar3H[Npt] =  {10,     10,   10,   20,  20,  100,  100,  150,  20,  100,  500,  2000,  200,  100, 100, 100,   10,   1000,  50,  10,  10,  100,  10,  100,  50,  10,  100};

    //for DCA_xy_tight cuts (ME)
    //  Double_t RBPar3L[Npt] = {-1000 , -100, -100, -100, -100, -100, -100, -150, -200, -100, -100, -2000, -200, -100, -20, -100, -20, -1000, -50, -10, -10, -100, -100, -100, -100, -51, -10};
    // Double_t RBPar3H[Npt] =  {1000,    100,  100,  100,  100,  100,  100,  150,  200,  100,  100,  2000,  200,  100, 20, 100,   20,    1000,  50,  10,  10,  100,  10,  100,  50,  10,  10};

    //for DCA_z_tight cuts (ME)
    //  Double_t RBPar3L[Npt] = {-10000 , -100, -100, -100, -10, -100, -100, -150, -200, -100, -100, -2000, -200, -100, -20, -100, -200, -200, -200, -10, -10, -200, -100, -100, -100, -51, -20};
    // Double_t RBPar3H[Npt] =  {10000,    100,  100,  100,  10,  100,  100,  150,  200,  100,  100,  2000,  200,  100, 20, 100,   200,    200,  200,  10,  10,  200,  10,  100,  50,  10,  20};

    //for DCA_all_tight (LS)
    //  Double_t RBPar3L[Npt] = {-10000 , -100, -100, -100, -10, -100, -100, -50, -10, -200, -10, -2000, -200, -100, -10, -100, -200, -200, -200, -100, -100, -200, -100, -100, -100, -51, -500};
    // Double_t RBPar3H[Npt] =  {10000,    100,  100,  100,  10,  100,  100,  50,  10,  200,  10,  2000,  200,  100, 10, 100,   200,   200,  200,  100,  100,  200,  10,  100,  50,  10,  500};


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
    int rebin = 2;                                            // for rebinning  //=1 for no rebin
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

    TCanvas *cgridls1 = new TCanvas("cgridls1", "cgridls1", 4000, 2300);
    cgridls1->Divide(4, 3);
    SetCanvasStyle2(cgridls1, 0.15, 0.05, 0.08, 0.13);
    TCanvas *cgridls2 = new TCanvas("cgridls2", "cgridls2", 4000, 2300);
    cgridls2->Divide(4, 3);
    SetCanvasStyle2(cgridls2, 0.15, 0.05, 0.08, 0.13);
    

    TCanvas *cgridme1 = new TCanvas("cgridme1", "cgridme1", 1200, 1200);
    cgridme1->Divide(4, 4);
    TCanvas *cgridme2 = new TCanvas("cgridme2", "cgridme2", 1200, 1200);
    cgridme2->Divide(4, 4);

    TCanvas *cbkgls1 = new TCanvas("cbkgls1", "cbkgls1", 2000, 1200);
    cbkgls1->Divide(4, 3);
    TCanvas *cbkgls2 = new TCanvas("cbkgls2", "cbkgls2", 2000, 1200);
    cbkgls2->Divide(4, 3);
    TCanvas *cbkgme1 = new TCanvas("cbkgme1", "cbkgme1", 1200, 1200);
    cbkgme1->Divide(4, 3);
    TCanvas *cbkgme2 = new TCanvas("cbkgme2", "cbkgme2", 1200, 1200);
    cbkgme2->Divide(4, 3);

    TCanvas *cinv[Npt]; // for output canvases on screen containing fitted signal after subtraction

    for (Int_t ip = 0; ip < Npt; ip++)
    {
        TString cName = TString::Format("cinv_pt_%2.1f-%2.1f", Low_pt[ip], High_pt[ip]);
        cinv[ip] = new TCanvas(Form("cinv%d", ip), cName.Data(), 10, 10, 1200, 1100);
        SetCanvasStyle2(cinv[ip], 0.15, 0.05, 0.08, 0.13);
    }

    TCanvas *cSigbkg[Npt]; // for output canvases on screen containing signal with bkg(after norm in case of mix)

    for (Int_t ip = 0; ip < Npt; ip++)
    {
        TString cNam = TString::Format("cSigbkg_pt_%2.1f-%2.1f", Low_pt[ip], High_pt[ip]);
        cSigbkg[ip] = new TCanvas(Form("cSigbkg%d", ip), cNam.Data(), 10, 10, 1200, 1100);
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

    // TFile *fInputFile = new TFile("../data/LHC22m_pass4_temp/PID_cuts/AnalysisResults.root", "Read"); 
    TFile *fInputFile = new TFile("../data/LHC22m_pass4_temp/PID/AnalysisResults.root", "Read"); 
    // TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(3);  // not there is PID for this dataset
    TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(4);
    // TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(5);
    // TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(6);
    // TKey *key = (TKey *)fInputFile->GetListOfKeys()->At(7);
    TList *fInputList = (TList *)fInputFile->Get(key->GetName());


    TFile *file2 = new TFile("../HEPData-ins1797443-v1-root.root", "READ");
    TKey *key2 = (TKey *)file2->GetListOfKeys()->At(3);
    TGraph *grsourav = (TGraph *)file2->Get("Table 4/Graph1D_y1");
    grsourav->SetTitle(0);
    grsourav->GetXaxis()->SetTitle(0);
    grsourav->GetYaxis()->SetTitle(0);

   // For DCA cuts
    // const string kOutputName = "lf-k892analysis_id4875";  //DCAz tight
    // const string kOutputName = "lf-k892analysis_id4876";  //DCAxy tight
    // const string kOutputName = "lf-k892analysis_id4877";  //DCAall tight

    //For PID cuts
    const string kOutputName = "lf-k892analysis_id3794";   //Default
    // const string kOutputName = "lf-k892analysis_id4878";   //PID3.5
    // const string kOutputName = "lf-k892analysis_id4879";   //PID 4
    // const string kOutputName = "lf-k892analysis_id4880";   //PID 2
    const string kOutputName2 = "QAbefore";
    const string kOutputName3 = "QAafter";

    //**To calculate total number of events for which histograms were filled*************************************************************
    // TH1F* hEVent = (TH1F *) fInputList->FindObject("hAEventsVsMulti");
    // Double_t Event=hEVent->Integral(lc,hc);
    Double_t Event = 1.702265e9;  // for PID cuts
    // Double_t Event = 1.708964e9;  // DCA cuts


    //*************************************************************************************************************************

    cout << "*****************number of events********************:" << Event << endl;

    //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************
    TH3F *fHistNum = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassDS", kOutputName.c_str()));
    // TH3F *fHistNum = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassDSAnti", kOutputName.c_str()));
    TH3F *fHistDen = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassME", kOutputName.c_str()));
    TH3F *fHistLS = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassLS", kOutputName.c_str()));
    // TH3F *fHistLS = (TH3F *)fInputFile->Get(Form("%s/h3k892invmassLSAnti", kOutputName.c_str()));

    // cout << " THE NUMBER OF BINS IN THE HISTOGRAM IS " << fHistNum->GetNbinsZ()<<endl;

    // // QA before

    TH1F *histtrkpT_pi = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_pi", kOutputName.c_str(), kOutputName2.c_str()));
    TH1F *histtrkpT_ka = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_ka", kOutputName.c_str(), kOutputName2.c_str()));
    TH1F *trDCAxy_ka = (TH1F *)fInputFile->Get(Form("%s/%s/trDCAxy_ka", kOutputName.c_str(), kOutputName2.c_str()));
    TH1F *trDCAz_ka = (TH1F *)fInputFile->Get(Form("%s/%s/trDCAz_ka", kOutputName.c_str(), kOutputName2.c_str()));
    TH1F *trDCAxy_pi = (TH1F *)fInputFile->Get(Form("%s/%s/trDCAxy_pi", kOutputName.c_str(), kOutputName2.c_str()));
    TH1F *trDCAz_pi = (TH1F *)fInputFile->Get(Form("%s/%s/trDCAz_pi", kOutputName.c_str(), kOutputName2.c_str()));
    TH2F *histTofTpc = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_TPC_Map_pi_all", kOutputName.c_str(), kOutputName2.c_str()));
    TH2F *toftpc_mappi_all = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_pi_all", kOutputName.c_str(), kOutputName2.c_str()));
    TH2F *histTpcpi = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigma_pi_all", kOutputName.c_str(), kOutputName2.c_str()));
    TH2F *histTofka = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_ka_all", kOutputName.c_str(), kOutputName2.c_str()));
    TH2F *histTpcka = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigmaka_all", kOutputName.c_str(), kOutputName2.c_str()));
     TH2F *toftpc_mapka_all = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_TPC_Mapka_all", kOutputName.c_str(), kOutputName2.c_str()));


    // // QA after

    TH1F *histtrkpT_pi2 = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_pi", kOutputName.c_str(), kOutputName3.c_str()));
    TH1F *histtrkpT_ka2 = (TH1F *)fInputFile->Get(Form("%s/%s/trkpT_ka", kOutputName.c_str(), kOutputName3.c_str()));
    TH1F *trkDCAxy_ka2 = (TH1F *)fInputFile->Get(Form("%s/%s/trDCAxy_ka", kOutputName.c_str(), kOutputName3.c_str()));
    TH1F *trDCAz_ka2 = (TH1F *)fInputFile->Get(Form("%s/%s/trDCAz_ka", kOutputName.c_str(), kOutputName3.c_str()));
    TH1F *trDCAxy_pi2 = (TH1F *)fInputFile->Get(Form("%s/%s/trDCAxy_pi", kOutputName.c_str(), kOutputName3.c_str()));
    TH1F *trDCAz_pi2 = (TH1F *)fInputFile->Get(Form("%s/%s/trDCAz_pi", kOutputName.c_str(), kOutputName3.c_str()));
    TH2F *histTofTpc2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_TPC_Map_pi_all", kOutputName.c_str(), kOutputName3.c_str()));
    TH2F *toftpc_mappi_all2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_pi_all", kOutputName.c_str(), kOutputName3.c_str()));
    TH2F *histTpcpi2 = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigma_pi_all", kOutputName.c_str(), kOutputName3.c_str()));
    TH2F *histTofka2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_Nsigma_ka_all", kOutputName.c_str(), kOutputName3.c_str()));
    TH2F *histTpcka2 = (TH2F *)fInputFile->Get(Form("%s/%s/TPC_Nsigmaka_all", kOutputName.c_str(), kOutputName3.c_str()));
    TH2F *toftpc_mapka_all2 = (TH2F *)fInputFile->Get(Form("%s/%s/TOF_TPC_Mapka_all", kOutputName.c_str(), kOutputName3.c_str()));


    //**********************************************************************************************************************************


    //**************************************************************************************************

      // ****Drawing histograms from the QAbefore and after folder

    // TCanvas *pq = new TCanvas("c", "Histograms", 900, 800);
    // SetCanvasStyle2(pq, 0.12, 0.12, 0.08, 0.13);
    // gStyle->SetOptStat(0);

    // // QAbefore
    // SetHistoQA(histtrkpT_pi);
    // pq->SetLogy();
    // histtrkpT_pi->Draw();
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDBtrkpT_pi.png");

    // SetHistoQA(histtrkpT_ka);
    // histtrkpT_ka->Draw();
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDBtrpt_ka.png");

    // // SetHistoQA(trDCAxy_ka);
    // // pq->SetLogy();
    // // trDCAxy_ka->Draw();
    // // pq->SaveAs("../output/LHC22m_pass4_temp/trDCAxy_ka.png");

    // // SetHistoQA(trDCAz_ka);
    // // pq->SetLogy();
    // // trDCAz_ka->Draw();
    // // pq->SaveAs("../output/LHC22m_pass4_temp/trDCAz_ka.png");

    // // SetHistoQA(trDCAxy_pi);
    // // pq->SetLogy();
    // // trDCAxy_pi->Draw();
    // // pq->SaveAs("../output/LHC22m_pass4_temp/trDCAxy_pi.png");

    // // SetHistoQA(trDCAz_pi);
    // // pq->SetLogy();
    // // trDCAz_pi->Draw();
    // // pq->SaveAs("../output/LHC22m_pass4_temp/trDCAz_pi.png");


    // SetHistoQA(histTofTpc);
    // pq->SetLogy(0);
    // pq->SetLogz();
    // histTofTpc->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDBtoftpc_pi_all.png");

    // SetHistoQA(toftpc_mappi_all);
    // toftpc_mappi_all->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/Btofpi.png");

    // SetHistoQA(histTpcpi);
    // histTpcpi->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDBtpcpi.png");

    // SetHistoQA(histTofka);
    // histTofka->Draw("colz");
    // histTofka->SetTitle("TOF NSigma for Kaon");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDBtofka.png");

    // SetHistoQA(histTpcka);
    // histTpcka->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDBtpcka.png");

    //  SetHistoQA(histTpcka);
    // histTpcka->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDBtpcka.png");

    // SetHistoQA(toftpc_mapka_all);
    // toftpc_mapka_all->SetTitle("TOF + TPC Combined PID for Kaon");
    // toftpc_mapka_all->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDBtoftpc_mapka_all.png");

    // // QA after
    // SetHistoQA(histtrkpT_pi2);
    // pq->SetLogy();
    // histtrkpT_pi2->Draw();
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDAtrkpT_pi.png");

    // SetHistoQA(histtrkpT_ka2);
    // histtrkpT_ka2->Draw();
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDAtrpt_ka.png");

    // SetHistoQA(histTofTpc2);
    // pq->SetLogy(0);
    // histTofTpc2->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDAtoftpc_mappi_all.png");

    // SetHistoQA(toftpc_mappi_all2);
    // toftpc_mappi_all2->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDAtofpi.png");

    // SetHistoQA(histTpcpi2);
    // histTpcpi2->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDAtpcpi.png");

    // SetHistoQA(histTofka2);
    // histTofka2->Draw("colz");
    // histTofka2->SetTitle("TOF NSigma for Kaon");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDAtofka.png");

    // SetHistoQA(histTpcka2);
    // histTpcka2->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDAtpcka.png");
    
    // SetHistoQA(toftpc_mapka_all2);
    // toftpc_mapka_all2->SetTitle("TOF + TPC Combined PID for Kaon");
    // toftpc_mapka_all2->Draw("colz");
    // pq->SaveAs("../output/LHC22m_pass4_temp/PIDAtoftpc_mapka_all.png");

    gstyle(); // this is not gStyle, it is defined in the header file style.h
    //***************************************************************************************************

    for (Int_t ip = 0; ip < Npt; ip++) // start pt bin loop
    // for (Int_t ip = 20; ip < Npt; ip++)
    {
        lowpt = Low_pt[ip];
        highpt = High_pt[ip];

        fHistTotal[ip] = fHistNum->ProjectionZ(Form("hSig_%d", ip), -1, -1, fHistNum->GetYaxis()->FindBin(lowpt), fHistNum->GetYaxis()->FindBin(highpt), "E");
        fHistBkg[ip] = fHistDen->ProjectionZ(Form("hbkg_%d", ip), -1, -1, fHistDen->GetYaxis()->FindBin(lowpt), fHistDen->GetYaxis()->FindBin(highpt), "E");
        fHistbkgLS[ip] = fHistLS->ProjectionZ(Form("hbkgLS_%d", ip), -1, -1, fHistLS->GetYaxis()->FindBin(lowpt), fHistLS->GetYaxis()->FindBin(highpt), "E");

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

        // TF1 *fitFcn = new TF1("fitfunc", BWExpo, lowfitrange[ip], highfitrange[ip], 7); // sig+bkg fit function
        // TF1 *fitFcn1 = new TF1("fitfunc1", Expo, lowfitrange[ip], highfitrange[ip], 4);    // only residualbkg
        // TF1 *fitFcn2 = new TF1("fitFcn2", BW, lowfitrange[ip], highfitrange[ip], 3);              // only signal
       
        TF1 *fitFcn = new TF1("fitfunc", BreitWignerpoly2, lowfitrange[ip], highfitrange[ip], 6); // sig+bkg fit function
        TF1 *fitFcn1 = new TF1("fitfunc1", polynomial2, lowfitrange[ip], highfitrange[ip], 3);    // only residualbkg
        TF1 *fitFcn2 = new TF1("fitFcn2", BW, lowfitrange[ip], highfitrange[ip], 3);              // only signal


        // TF1 *fitFcn = new TF1("fitfunc", voigtpol2, lowfitrange[ip], highfitrange[ip], 7); // sig+bkg fit function
        // TF1 *fitFcn1 = new TF1("fitfunc1", polynomial2, lowfitrange[ip], highfitrange[ip], 3);    // only residualbkg
        // TF1 *fitFcn2 = new TF1("fitFcn2", voigt, lowfitrange[ip], highfitrange[ip], 4);              // only signal


        fitFcn->SetParLimits(0, 0.83, 0.92); // Mass
        fitFcn->SetParLimits(2, 0, 1000000000);  // Yield
        fitFcn->SetParameter(1, 0.047);      // width
        // fitFcn->SetParLimits(1, 0.033, 0.061);  // width
        fitFcn->FixParameter(1, 0.047);      // width
        fitFcn->SetParameter(2, 1000); // yield

        //for voigitian function

        // fitFcn->SetParLimits(0, 1e1, 1.e8);  //yield
        // fitFcn->SetParameter(1, 0.895); //mass
        // fitFcn->SetParLimits(1, 0.83, 0.9); //mass
        // fitFcn->SetParameter(2, 0.004); //resolution
        // fitFcn->SetParameter(3, 0.0487);  //width
        // fitFcn->FixParameter(3, 0.0487);  


        // fitFcn->SetParLimits(3, RBPar3L[ip], RBPar3H[ip]);
        // fitFcn->SetParLimits(4, -10000, 10000);
        // fitFcn->SetParLimits(5, RBPar5L[0], RBPar5H[0]);
        // fitFcn->SetParLimits(6, RBPar6L[0], RBPar6H[0]);
        

        // fitFcn->SetParNames("Mass","Width","Yield","A","B","C","D");
        fitFcn->SetParNames("Mass", "Width", "Yield", "C", "B", "A");
        // fitFcn->SetParNames("Yield","Mass","resolution","Width","C","B","A");
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
        SetHistoStyle(hfsig, 1, 20, 1, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
        SetHistoStyle(fHistTotal[ip], 1, 8, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

        hfsig->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
        hfsig->GetYaxis()->SetMaxDigits(2);
        hfsig->GetYaxis()->SetTitle("Counts");

        if (ResBkg.CompareTo("MIX") == 0)
        {

            // SetHistoStyle(hfbkg, 2, 24, 0.5, 0.05, 0.05, 0.03, 0.03, 1.0, 1.4);
            SetHistoStyle(hfbkg, kRed, 24, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

            hfbkg->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
            hfsig->GetYaxis()->SetMaxDigits(2);
            hfbkg->GetYaxis()->SetTitle("Counts");
        }

        else if (ResBkg.CompareTo("LIKE") == 0)
        {
            // SetHistoStyle(hfbkg_like,2,24,0.5,0.05,0.05,1.0,1.1);
            SetHistoStyle(hfbkg_like, 2, 4, 1.5, 0.05, 0.05, 0.03, 0.03, 1.0, 1.4);
            hfbkg_like->GetXaxis()->SetTitle("M_{K#pi} (GeV/#it{c}^{2})");
            hfsig->GetYaxis()->SetMaxDigits(2);
            hfbkg_like->GetYaxis()->SetTitle("Counts");
        }

        fitFcn1->SetLineColor(4);
        fitFcn1->SetLineStyle(2);
        fitFcn1->SetLineWidth(4);
        fitFcn2->SetLineColor(6);
        fitFcn2->SetLineStyle(2);
        fitFcn2->SetLineWidth(4);
        fitFcn->SetLineWidth(4);

        //*******************************************************************************************************************************

        //**Plot of histograms and graphs*********************************************************************************************
        
        cinv[ip]->cd(); 
        TLegend *pag = new TLegend(0.16,0.8, 0.28, 0.9);
        hfsig->Draw("e");
        fitFcn->Draw("same");
        fitFcn1->Draw("same");
        fitFcn2->Draw("same");
        pag->AddEntry(fitFcn, "BW+pol2");
        pag->AddEntry(fitFcn1, "BW");
        pag->AddEntry(fitFcn2, "pol2");
        pag->Draw();
        t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
        t2->DrawLatex(0.65, 0.97,
                      Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
                           High_pt[ip]));
        hfsig->GetXaxis()->SetRangeUser(al, bh);
        cout<<"\n\n\n"<<endl;
        cinv[ip]->SaveAs(Form("../output/LHC22m_pass4_temp/PIDhfitsig_pt%d.png", ip + 1));
        cout<<"\n\n\n"<<endl;

        if (ip < 12)
        {
            cgridls1->cd(ip + 1);
            hfsig->Draw("e");
            fitFcn->Draw("same");
            fitFcn1->Draw("same");
            fitFcn2->Draw("same");
            t2->DrawLatex(0.65, 0.97,
                          Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
                               High_pt[ip]));
        }
        else if (ip >= 12 && ip <25)
        {
            cgridls2->cd(ip - 11);
            hfsig->Draw("e");
            fitFcn->Draw("same");
            fitFcn1->Draw("same");
            fitFcn2->Draw("same");
            t2->DrawLatex(0.65, 0.97,
                          Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
                               High_pt[ip]));
        }

      

        // // signal extraction before bkg subtraction////////////////////////////////////////////////////////////////

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

            // Double_t max = hfbkg->GetMaximum();
            // Double_t min = fHistTotal[ip]->GetMaximum();
            /*if (max>min)
            hfbkg->GetYaxis()->SetRangeUser(0,max);
            else
            hfbkg->GetYaxis()->SetRangeUser(0,min);
            */
            // cSigbkg[ip]->SaveAs(Form("../output/LHC22m_pass4_temp/PIDhsigbkg_pt%d.png", ip + 1));

            if (ip < 12)
            {
                cbkgme1->cd(ip + 1);
                fHistTotal[ip]->Draw("e");
                hfbkg->Draw("same");
                t2->DrawLatex(0.65, 0.97,
                              Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
                                   High_pt[ip]));
            }
            else if (ip >11 && ip < 24)
            {
                cbkgme2->cd(ip - 11);
                fHistTotal[ip]->Draw("e");
                hfbkg->Draw("same");
                t2->DrawLatex(0.65, 0.97,
                              Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
                                   High_pt[ip]));
            }

            
        }
        

        ////////////////////////////////////////////////////////////////////////

        // For like sign bkg///////////////////////////////////////////////////
        else if (ResBkg.CompareTo("LIKE") == 0)
        {
            hfbkg_like->GetXaxis()->SetRangeUser(al, bh);
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

            // cSigbkg[ip]->SaveAs(Form("../output/LHC22m_pass4_temp/PIDhsigbkg_pt%d.png", ip + 1));

            if (ip < 12)
            {
                cbkgls1->cd(ip + 1);
                fHistTotal[ip]->Draw("e");
                hfbkg_like->Draw("same");
                t2->DrawLatex(0.65, 0.97,
                              Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
                                   High_pt[ip]));
            }
            else if (ip >= 12 && ip < 25)
            {
                cbkgls2->cd(ip - 11);
                fHistTotal[ip]->Draw("e");
                hfbkg_like->Draw("same");
                t2->DrawLatex(0.65, 0.97,
                              Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", Low_pt[ip],
                                   High_pt[ip]));
            }

            
        }

        ////////////////////////////////////////////////////////////////////////

    } // pt loop ends

      cgridls1->SaveAs("../output/LHC22m_pass4_temp/DCA_all_tight_ME_pol2_sig1.jpg");
      cgridls2->SaveAs("../output/LHC22m_pass4_temp/DCA_all_tight_ME_pol2_sig2.jpg");
    //   cbkgls1->SaveAs("../output/LHC22m_pass4_temp/DCA_all_tight_LS_pol2_bkg1.png");
    //   cbkgls2->SaveAs("../output/LHC22m_pass4_temp/DCA_all_tight_LS_pol2_bkg2.png");
    //   cbkgme1->SaveAs("../output/LHC22m_pass4_temp/DCA_all_tight_ME_pol2_bkg1.png");
    //   cbkgme2->SaveAs("../output/LHC22m_pass4_temp/DCA_all_tight_ME_pol2_bkg2.png");

    TFile *filecmp = new TFile("../output/LHC22m_pass4_temp/DCA_all_tight_ME_pol2.root", "RECREATE");


    TCanvas *sig = new TCanvas ();
    SetCanvasStyle2(sig, 0.2, 0.05, 0.08, 0.2);
    hsgnfcance->Draw();
    hsgnfcance->Write("significance");
    sig->SaveAs("../output/LHC22m_pass4_temp/significance.png");

    TCanvas *cgr = new TCanvas("cgr", "graphs", 900, 800);
    SetCanvasStyle2(cgr, 0.2, 0.05, 0.08, 0.2);

    // // chisquare_NDF vs pt
    TGraph *grchi = new TGraph(Npt, x, Chi2Ndf);
    grchi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    grchi->GetYaxis()->SetTitle("#chi^{2}/NDF ");
    SetGraphStyle(grchi, 1, 1);
    grchi->Draw("AP");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    cgr->SaveAs("../output/LHC22m_pass4_temp/chindfvspt.png");
    grchi->Write("chils");

    // mass vs pt
    TGraphErrors *grmass = new TGraphErrors(Npt, x, Mass, 0, ErrorMass);
    grmass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    grmass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    SetGrapherrorStyle(grmass, 1, 1);
    grmass->Draw("ap");
    // grmass->Write("massls");
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
    cgr->SaveAs("../output/LHC22m_pass4_temp/mass_pt.png");

    // // Width vs pT

    TGraphErrors *grwidth = new TGraphErrors(Npt, x, Width, 0, ErrorWidth);
    grwidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    grwidth->GetYaxis()->SetTitle("Width (GeV)");
    SetGrapherrorStyle(grwidth, 1, 1);
    grwidth->Draw("ap");
    // grwidth->Write("widthls");
    TLegend *widthleg = new TLegend(0.45, 0.75, 0.7, 0.85);
    widthleg->SetFillColor(0);
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    TLine *line2 = new TLine(grmass->GetXaxis()->GetXmin(), 0.047, grmass->GetXaxis()->GetXmax(), 0.047);
    line2->SetLineStyle(2);
    line2->SetLineColor(2);
    line2->SetLineWidth(3);
    line2->Draw();
    widthleg->AddEntry(line, "PDG Width", "l");
    widthleg->Draw("l");
    // cgr->SaveAs("../output/LHC22m_pass4_temp/width_pt.png");

    // // yield vs pt

    TCanvas *cyield1 = new TCanvas("cyield1", "yeild", 1000, 800);
    gPad->SetLogy();
    SetCanvasStyle2(cyield1, 0.2, 0.05, 0.08, 0.2);
    TGraphErrors *gryield1 = new TGraphErrors(Npt, x, yieldcalc1, xerr, yielderror1);
    TLegend *legyield = new TLegend(0.5, 0.7, 0.9, 0.9);
    gryield1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gryield1->GetYaxis()->SetTitle("Yield (GeV/c)^{-1}");
    SetGrapherrorStyle(gryield1, 1, 1);
    grsourav->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    grsourav->GetYaxis()->SetTitle("(1/N_{ev})* d^{2}N/(dp_{T} dy) (Gev/c)^{-1}");
    grsourav->GetYaxis()->SetRangeUser(10e-8, 1);
    SetGraphStyle(grsourav, 2, 2);
    legyield->AddEntry(grsourav, "pp 13TeV (Published)");
    legyield->AddEntry(gryield1, "pp 13.6TeV (LHC22m_pass4_temp)");

     grsourav->Draw("ap");
    gryield1->Draw("psame");
    gryield1->Write("yieldls");
    // grsourav->Write("yieldpbl");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    legyield->SetTextSize(0.03);
    legyield->SetTextFont(2);
    legyield->Draw();
    cyield1->SaveAs("../output/LHC22m_pass4_temp/yield1.png");

    // // ratio of calculated by publsihed yields

    double valuey[Npt];
    double valueyerr[Npt];
    double ratioy[Npt];
    double ratioyerr[Npt];
    for (int i = 0; i < Npt; i++)
    {
        valuey[i] = grsourav->GetY()[i];
        // cout << "THE VALUE OF YIELD FROM SOURAV ANAYSIS  " << valuey[i] << endl;
        valueyerr[i] = grsourav->GetErrorY(i);
        ratioy[i] = yieldcalc1[i] / valuey[i];
        ratioyerr[i] = ratioy[i] * TMath::Sqrt((valueyerr[i] * valueyerr[i]) / (valuey[i] * valuey[i]) + (yielderror1[i] * yielderror1[i]) / (yieldcalc1[i] * yieldcalc1[i])); // error propagation
    }

    TCanvas *cratio = new TCanvas("yield", "yeild", 1000, 800);
    SetCanvasStyle2(cratio, 0.2, 0.05, 0.08, 0.2);

    TGraphErrors *grratio = new TGraphErrors(Npt, x, ratioy, 0, ratioyerr);
    grratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    grratio->GetYaxis()->SetTitle("Pseudo Efficiency (This Analysis/Published)");
    SetGrapherrorStyle(grratio, 1, 1);

    grratio->Draw("ap");
    // grratio->Write("ratiols");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    cratio->SaveAs("../output/LHC22m_pass4_temp/ratio_yield.png");

    // filecmp->Close();
}
