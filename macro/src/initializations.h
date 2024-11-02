#include "common.h"
#include "fitting_range.h"
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

//**Variables and arrays needed *******************************************************************************************

Double_t integralsignalfunc[Npt];                         // for calculation of area under signal only (after fitting and extraction BW signal)
Double_t ptcenter[Npt];                                   // stores mean pt for each pt bin
Double_t ptbinwidth[Npt];                                 // stores pt binwidth
Double_t interror[100];                                   // stores error while area calculation of signal for yield
const int N = 1;                                          // for centrality loop below used
Double_t dy = 1.0;                                        // for rapidity difference
TFitResultPtr r;                                          // for fitting using TFitter
TH1D *fHistTotal[Npt];                                    // for sig+bg
TH1D *fHistTotal_anti[Npt];                               // for sig+bg_anti
TH1D *fHistBkg[Npt];                                      // for mixedbg
TH1D *fHistbkgLS[Npt];                                    // for like sign
TH1D *fHistbkgLS_anti[Npt];                               // for like sign_anti
TH1D *fHistRotated1D[Npt];                                  // for rotated bkg
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

int bmin, bmax;
Double_t hBCError_1, bkgvalue, Integral_BW_withsigma, fYield_BinCount, YieldIntegral_BW, Yfraction_cBW, sum_tail_correction, Total_Ybincounting, Tail_correction_plusm, Tail_correction_minusm, Error_2, Final_pro_error;
Double_t nlow, nhigh;
Double_t Yield_bincount_hist;

//***************************************************************************************************************************

//**Canvas definitions and initialisations********************************************************************************

//***************************************************************************************************

