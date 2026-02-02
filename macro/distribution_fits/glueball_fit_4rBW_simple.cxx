#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

// Forward declarations
void glueball_fit_4rBW_simple();
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

int main()
{
    glueball_fit_4rBW_simple();
    return 0;
}

Double_t single_BW(double *x, double *par);
Double_t BWsumMassDepWidth(double *x, double *par);
Double_t BWsumMassDepWidth_exponential(double *x, double *par);
Double_t single_BW_mass_dep_spin0(double *x, double *par);
Double_t single_BW_mass_dep_spin2(double *x, double *par);
Double_t exponential_bkg_3(double *x, double *par); // 4 parameters
Double_t BWsum_expol_chkstar(double *x, double *par);
Double_t expol_chkstar(double *x, double *par);
Double_t BWsumMassDepWidth_simple_exponential(double *x, double *par);
Double_t simple_exponential(double *x, double *par);
Double_t BWsum_ConstWidth_ModBolt(double *x, double *par);
Double_t BWsum_constWidth(double *x, double *par);

void glueball_fit_4rBW_simple()
{
    // Note that the main source of systematic i.e. signal extraction is still remaining to be done
    std::vector<std::string> variations = {
        "", "_DCA0p1", "_TPCPID2", "_TPCPID5", "_TPCMinCls100", "_TPCMinCls60", "_DCAv0dau0p3", "_DCAv0dau1p0", "_Ks_selection2p5", "_Ks_selection5", "_cospa0p95", "_cospa0p992", "_decay_rad1p0", "_lambda_rej4", "_lambda_rej6", "_lifetime15", "_lifetime25"}; // All variations
    int totalVar = variations.size();

    for (int isysvars = 0; isysvars < 1; isysvars++)
    {

        //****************systematics train*******************************
        const string kvariation1 = variations[isysvars];
        //======== For fitting and normalization range variations ==========
        // const string kvariation1 = "_fitLow1p07";
        // const string kvariation1 = "_fitHigh2p17";
        // const string kvariation1 = "_fitHigh2p25"; // Low fit is 1.045 (without it same as default)
        // const string kvariation1 = "_normLeft"; //2.4-2.5 (fitlow 1.07 for 7-10 GeV/c)
        // const string kvariation1 = "_normRight"; // 2.6-2.7
        // const string kvariation1 = "_AllParametersFree";
        // const string kvariation1 = "_AllParametersFixed";
        //======== For fit function variation =============================
        // const string kvariation1 = "_FitChKstar";
        // const string kvariation1 = "_FitExpoHERA";
        // const string kvariation1 = "_ConstWidth";

        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/systematic2022_new/KsKs_Channel/higher-mass-resonances" + kvariation1; // for systematics studies (excluding signal extraction)

        string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/systematic2022_new/KsKs_Channel/higher-mass-resonances"; // for systematics studies (excluding signal extraction)

        //******************for default study with full train ************************
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances"; // 2022 dataset (old)
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/435448/KsKs_Channel/higher-mass-resonances"; //2022 dataset (new)
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/435450/KsKs_Channel/higher-mass-resonances"; // 2023 dataset
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/435449/KsKs_Channel/higher-mass-resonances"; // 2024 dataset
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/504802/KsKs_Channel/higher-mass-resonances_CS_Frame"; // 2023 dataset (derived data)
        string path2 = path;

        string savepath = path2 + "/fits/";

        gSystem->Exec(("mkdir -p " + savepath).c_str());

        TFile *f = new TFile((path + "/hglue_ROTATED_allPt.root").c_str(), "READ");
        // TFile *f = new TFile((path + "/hglue_ROTATED_allPt_normleft.root").c_str(), "READ");
        int colors[] = {kGreen + 4, 28, kMagenta, kBlue};
        double masses[] = {f1270Mass, a1320Mass, f1525Mass, f1710Mass};
        double widths[] = {f1270Width, a1320Width, f1525Width, f1710Width};
        string resonance_names[] = {"f_{2}(1270)", "a_{2}(1320)^{0}", "f'_{2}(1525)", "f_{0}(1710)"};
        string resonance_mass[] = {"1270", "1320", "1525", "1710"};

        if (f->IsZombie())
        {
            cout << "Error opening file" << endl;
            return;
        }
#define b_massdepWidth_modifiedBoltzmann
#define residual_subtracted
// #define doublepanelplot
#define multiPanelPlots
        // #define before_combinatorialPlots
        // #define singlePanelPlots

        TH1F *hmult = (TH1F *)f->Get("multiplicity_histogram");
        if (hmult == nullptr)
        {
            cout << "Multiplicity histogram not found" << endl;
            return;
        }
        int multlow, multhigh;
        double phi_mod, phi_mod2, total_events;

#ifdef multiPanelPlots
        TCanvas *cMultiPanelFit = new TCanvas("cMultiPanelFit", "Multi-Panel Fit Results", 1440, 720);
        cMultiPanelFit->Divide(3, 2);
        TCanvas *cMultiPanelResidual = new TCanvas("cMultiPanelResidual", "Multi-Panel Residuals", 1440, 720);
        cMultiPanelResidual->Divide(3, 2);
        TCanvas *cMultiPanelWithBkg = new TCanvas("cMultiPanelWithBkg", "Multi-Panel Fit with Bkg Results", 1440, 720);
        cMultiPanelWithBkg->Divide(3, 2);
#endif
        // (2022 data)
        vector<vector<pair<float, float>>> fitRanges = {
            {{1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}},  // 0-100%
            {{1.05, 2.20}, {1.09, 2.15}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}},  // 0-20%
            {{1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.03, 2.20}, {1.05, 2.20}},  // 20-50%
            {{1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}},  // 50-70%
            {{1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}}}; // 70-100%

        // // (2023 data)
        // vector<vector<pair<float, float>>> fitRanges = {
        //     {{1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}},  // 0-100%
        //     {{1.06, 2.20}, {1.09, 2.15}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}},  // 0-20%
        //     {{1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.03, 2.20}, {1.07, 2.20}},  // 20-50%
        //     {{1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}},  // 50-70%
        //     {{1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}, {1.05, 2.20}}}; // 70-100%

        TFile *fDefault = new TFile((savepath + "/FitParam" + kvariation1 + ".root").c_str(), "recreate");

        // for (int imult = 0; imult < Nmult + 1; imult++)
        for (int imult = 0; imult < 1; imult++)
        {
            if (imult == 0)
            {
                multlow = 0;
                multhigh = 100; // for all multiplicity
            }
            else
            {
                multlow = mult_classes[imult - 1];
                multhigh = mult_classes[imult];
            }
            // All the required histograms to be stored
            TH1F *hChi2NDF = new TH1F("hChi2NDF", "Chi2/NDF vs pT", Npt, pT_bins);
            TH1F *hSignificance[4];
            TH1F *hStatSignificance[4];
            TH1F *hMass[4];
            TH1F *hWidth[4];
            TH1F *hYield[4];

            for (int ires = 0; ires < 4; ires++)
            {
                hSignificance[ires] = new TH1F(Form("hSignificance%d", (int)masses[ires]), Form("Significance of %s vs pT", resonance_names[ires].c_str()), Npt, pT_bins);
                hStatSignificance[ires] = new TH1F(Form("hStatSignificance%d", (int)masses[ires]), Form("Statistical Significance of %s vs pT", resonance_names[ires].c_str()), Npt, pT_bins);
                hMass[ires] = new TH1F(Form("hMass%d", (int)masses[ires]), Form("Mass of %s vs pT", resonance_names[ires].c_str()), Npt, pT_bins);
                hWidth[ires] = new TH1F(Form("hWidth%d", (int)masses[ires]), Form("Width of %s vs pT", resonance_names[ires].c_str()), Npt, pT_bins);
                hYield[ires] = new TH1F(Form("hYield%d", (int)masses[ires]), Form("Raw Yield of %s vs pT", resonance_names[ires].c_str()), Npt, pT_bins);
            }

            // // Temporary for single multiplicity class checking
            // multlow = 0;
            // multhigh = 100;

            total_events = hmult->Integral(hmult->GetXaxis()->FindBin(multlow + 0.01), hmult->GetXaxis()->FindBin(multhigh - 0.01));

            string savepath_mult = savepath + Form("/mult_%d-%d", multlow, multhigh);
            gSystem->Exec(("mkdir -p " + savepath_mult).c_str());
            float maxRanges[] = {1.25, 1.4, 1.6, 1.8, 1.8, 2.0, 2.2};

            for (int ipt = 0; ipt < Npt; ipt++)
            // for (int ipt = 2; ipt < 3; ipt++)
            {
                float lowpT = pT_bins[ipt];
                float highpT = pT_bins[ipt + 1];

                // // Temporary for single bins checking
                // float lowpT = 10.0;
                // float highpT = 15.0;

                // ofstream file;
                // file.open((savepath_mult + Form("/fit_params_pT_%.1f-%.1f", lowpT, highpT) + kvariation1 + ".txt").c_str());

                // for (int irange = 0; irange < fitranges.size(); irange++)
                for (int irange = 0; irange < 1; irange++)
                {
                    // float fitlow = fitranges[irange][0];
                    // float fithigh = fitranges[irange][1];
                    float fitlow = fitRanges[imult][ipt].first;
                    float fithigh = fitRanges[imult][ipt].second;
                    cout << "fit low is " << fitlow << endl;
                    cout << "fit high is " << fithigh << endl;

                    cout << "Low fit range is " << fitlow << ", High fit range is " << fithigh << endl;

                    TH1F *hinvMass = (TH1F *)f->Get(Form("multiplicity_%d_%d/ksks_subtracted_invmass_pt_%.1f_%.1f", multlow, multhigh, lowpT, highpT));
                    // TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", lowpT, highpT));
                    TH1F *hraw = (TH1F *)f->Get(Form("multiplicity_%d_%d/ksks_invmass_pt_%.1f_%.1f", multlow, multhigh, lowpT, highpT));
                    // TH1F *hraw = (TH1F *)f->Get(Form("ksks_invmass_pt_%.1f_%.1f", lowpT, highpT));
                    TH1F *hbkg = (TH1F *)f->Get(Form("multiplicity_%d_%d/ksks_bkg_pt_%.1f_%.1f", multlow, multhigh, lowpT, highpT));
                    if (hinvMass == nullptr || hraw == nullptr || hbkg == nullptr)
                    {
                        cout << "Error opening histogram: " << Form("multiplicity_%d_%d/ksks_subtracted_invmass_pt_%.1f_%.1f", multlow, multhigh, lowpT, highpT) << endl;
                        return;
                    }

// //TCanvas *c_reducedFit = new TCanvas("c_reducedFit", "Reduced Fit", 720, 720);
// //SetCanvasStyle(c_reducedFit, 0.14, 0.03, 0.05, 0.13);
#ifdef singlePanelPlots
                    TCanvas *c = new TCanvas("", "", 720, 720);
                    SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
                    c->cd();
#endif
#ifdef multiPanelPlots
                    cMultiPanelFit->cd(ipt + 1);
                    gPad->SetTopMargin(0.06);
                    gPad->SetBottomMargin(0.14);
                    gPad->SetLeftMargin(0.13);
                    gPad->SetRightMargin(0.05);
#endif
                    // hinvMass->Rebin(2);

                    double binwidthfile = (hinvMass->GetXaxis()->GetXmax() - hinvMass->GetXaxis()->GetXmin()) / hinvMass->GetXaxis()->GetNbins();
                    cout << "Binwidth file is " << binwidthfile << endl;
                    cout << "bin width is " << binwidthfile << endl;
                    hinvMass->GetXaxis()->SetRangeUser(1.01, 2.20);
                    hinvMass->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
                    hinvMass->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidthfile * 1000));
                    TH1F *hinvMassClone = (TH1F *)hinvMass->Clone("hinvMassClone");
                    hinvMass->SetMaximum(maxRanges[ipt] * hinvMass->GetMaximum());
                    // hinvMass->SetMaximum(2.0 * hinvMass->GetMaximum());
#ifdef multiPanelPlots
                    hinvMass->GetYaxis()->SetTitleOffset(1.15);
                    hinvMass->SetMarkerSize(0.8);
#else
                    hinvMass->GetYaxis()->SetTitleOffset(1.35);
                    hinvMass->SetMarkerSize(1.3);
#endif
                    hinvMass->Draw("pe");
                    double resMaximumFactor = 2.2;
                    TH1F *hsubtracted = (TH1F *)hinvMass->Clone("hsubtracted");
                    TH1F *hsubtracted_res = (TH1F *)hinvMass->Clone("hsubtracted_res");
                    // gStyle->SetOptStat(1110);
                    gStyle->SetOptStat(0);
                    gStyle->SetOptFit(1111);
                    vector<tuple<float, int, float, float>> fit_parameters;

                    // // // //************************************************************************ */
                    // // // // **************** For BW sum with expol HERA ****************************

                    // ************************** fit with mass depndent width BW ***************************************

#ifdef b_massdepWidth_modifiedBoltzmann

                    //======================Default fit function==============================
                    TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_exponential, fitlow, fithigh, 16); // expol 3
                    TF1 *BEexpol_initial = new TF1("BEexpol_initial", BWsumMassDepWidth_exponential, fitlow, fithigh, 16);

                    // //======================Background function Charge K*======================
                    // TF1 *BEexpol = new TF1("BEexpol", BWsum_expol_chkstar, fitlow, fithigh, 16);
                    // TF1 *BEexpol_initial = new TF1("BEexpol_initial", BWsum_expol_chkstar, fitlow, fithigh, 16);

                    // //======================Background function expo HERA======================
                    // TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_simple_exponential, fitlow, fithigh, 15);
                    // TF1 *BEexpol_initial = new TF1("BEexpol_initial", BWsumMassDepWidth_simple_exponential, fitlow, fithigh, 15);

                    // //======================BW with constant width=============================
                    // TF1 *BEexpol = new TF1("BEexpol", BWsum_ConstWidth_ModBolt, fitlow, fithigh, 16);
                    // TF1 *BEexpol_initial = new TF1("BEexpol_initial", BWsum_ConstWidth_ModBolt, fitlow, fithigh, 16);

                    string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"}; // remove name "d" when using expo HERA
                    for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
                    {
                        BEexpol->SetParName(i, parnames[i].c_str());
                    }

                    double parameters[] = {1.1e4, f1270Mass, f1270Width, 5.8e3, a1320Mass, a1320Width, 1.6e4, f1525Mass, f1525Width, 3.1e3, f1710Mass, f1710Width};

                    if ((imult == 3 && ipt == 3) || (imult == 4 && ipt == 3) || (imult == 4 && ipt == 4))
                    {
                        parameters[0] = 3500;
                        parameters[3] = 2000;
                        parameters[6] = 7000;
                        parameters[9] = 2200;
                    }

                    int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

                    for (int i = 0; i < size_fitparams; i++)
                    {
                        BEexpol->SetParameter(i, parameters[i]);
                        BEexpol_initial->SetParameter(i, parameters[i]);
                    }

                    // //********systematic studies*************
                    double initial_param_bkg[] = {1.37518e5, 0.6, 6.071167, 1.04}; // default fit
                    // double initial_param_bkg[] = {9.5e6, -0.007, -2.4, -0.15}; // for expol charged K*
                    // double initial_param_bkg[] = {6.8e5, -0.2, 2.8}; // for Exponential HERA
                    // double initial_param_bkg[] = {1.4e5, -0.006, 6.071167, 1.0}; // const width BW

                    // Initial parameters for background
                    BEexpol_initial->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Free
                    BEexpol_initial->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Fix for medium train
                    BEexpol_initial->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Free
                    BEexpol_initial->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free
                    // BEexpol_initial->FixParameter(size_fitparams + 3, 1.1); // Fix for Constant width case

                    BEexpol_initial->FixParameter(2, f1270Width);
                    BEexpol_initial->FixParameter(5, a1320Width);
                    BEexpol_initial->FixParameter(8, f1525Width);

                    BEexpol_initial->FixParameter(1, f1270Mass);
                    BEexpol_initial->FixParameter(4, a1320Mass);
                    BEexpol_initial->FixParameter(7, f1525Mass);

                    BEexpol_initial->FixParameter(10, f1710Mass);
                    BEexpol_initial->FixParameter(11, f1710Width);
                    // TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "RELMS0");
                    hinvMass->Fit("BEexpol_initial", "REBMS0"); // fit with all fixed to get bkg params
                    // TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBS"); // comment while using toy mc and likelihood fits

                    //=============================================================//
                    // Again setting parameters for the next iteration in the fit.
                    //=============================================================//
                    for (int iparams = 0; iparams < 16; iparams++)
                    {
                        BEexpol->SetParameter(iparams, BEexpol_initial->GetParameter(iparams));
                    }

                    vector<vector<double>> par_limits = {{1, 1 * f1270Width}, {2, 3 * f1270WidthErr}, {4, 1 * a1320Width}, {5, 5 * a1320WidthErr}, {7, 1 * f1525Width}, {8, 3 * f1525WidthErr}, {10, 0.45 * f1710Width}, {11, 3 * f1710WidthErr}};

                    int limits_size = par_limits.size();
                    for (int i = 0; i < limits_size; i++)
                    {
                        int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                        BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
                    }

                    if (!(imult == 4 && ipt == 4) && !(imult == 1 && ipt == 0)) // second condition added for temp fix for 2023 data
                    {
                        BEexpol->SetParLimits(0, 0, 1e6);
                        BEexpol->SetParLimits(3, 0, 1e6);
                    }
                    // BEexpol->SetParLimits(6, 0, 1e6);
                    // BEexpol->SetParLimits(9, 0, 1e6);
                    // BEexpol->FixParameter(size_fitparams + 3, 1.0);

                    BEexpol->FixParameter(2, f1270Width);
                    BEexpol->FixParameter(5, a1320Width);
                    BEexpol->FixParameter(8, f1525Width);
                    // BEexpol->SetParameter(2, f1270Width);
                    // BEexpol->SetParameter(5, a1320Width);
                    // BEexpol->SetParameter(8, f1525Width);

                    // BEexpol->FixParameter(1, f1270Mass);
                    // BEexpol->FixParameter(4, a1320Mass);
                    // BEexpol->FixParameter(7, f1525Mass);

                    // BEexpol->FixParameter(10, f1710Mass);
                    BEexpol->FixParameter(11, f1710Width);

                    TFitResultPtr fitResultptr;
                    if (ipt == 0)
                        fitResultptr = hinvMass->Fit("BEexpol", "REBS"); // comment while using toy mc and likelihood fits
                    else
                        fitResultptr = hinvMass->Fit("BEexpol", "REBMS"); // comment while using toy mc and likelihood fits

                    double *obtained_parameters = BEexpol->GetParameters(); // comment while using toy mc and likelihood fits

                    //========================For default fit====================================
                    TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
                    TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); // Do not forget to change number of parameters in below for loop

                    // //=======================For Expol charged K*====================================
                    // TF1 *expol = new TF1("expol", expol_chkstar, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
                    // TF1 *expol_clone = new TF1("expol_clone", expol_chkstar, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);

                    // //=======================For Exponential (HERA)============================
                    // // since it has 3 parameters, changed the loop number also below
                    // TF1 *expol = new TF1("expol", simple_exponential, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                    // TF1 *expol_clone = new TF1("expol_clone", simple_exponential, BEexpol->GetXmin(), BEexpol->GetXmax(), 3); // Do not forget to change number of parameters in below for loop to 3

                    for (int i = 0; i < 4; i++)
                    {
                        expol->SetParameter(i, obtained_parameters[size_fitparams + i]);
                        expol_clone->SetParameter(i, obtained_parameters[size_fitparams + i]);
                    }
                    expol->SetLineColor(3);
                    expol->SetLineStyle(2);
                    expol_clone->SetLineColor(3);
                    expol_clone->SetLineStyle(2);
                    expol->Draw("same");

                    // // //=====================Default fit=======================================
                    // TF1 *onlyBW = new TF1("onlyBW", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
                    // TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);

                    // // =======================BW with constant width=============================
                    TF1 *onlyBW = new TF1("onlyBW", BWsum_constWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
                    TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum_constWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);

                    string parameter_names[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma"};
                    for (int i = 0; i < 12; i++)
                    {
                        onlyBW->SetParameter(i, obtained_parameters[i]);
                        onlyBW_clone->SetParameter(i, obtained_parameters[i]);
                        onlyBW_clone->SetParName(i, parameter_names[i].c_str());
                    }
                    for (int i = 0; i < limits_size; i++)
                    {
                        int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                        onlyBW->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
                        onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
                    }

                    onlyBW_clone->SetParLimits(0, 0.0, 1e6); // norm1270
                    onlyBW_clone->SetParLimits(3, 0.0, 1e6); // norm1320
                    // onlyBW_clone->SetParLimits(6, 0.0, 1e6);   // norm1525

                    onlyBW_clone->FixParameter(2, f1270Width);
                    onlyBW_clone->FixParameter(5, a1320Width);
                    onlyBW_clone->FixParameter(8, f1525Width);

                    // onlyBW_clone->FixParameter(1, f1270Mass);
                    // onlyBW_clone->FixParameter(4, a1320Mass);
                    // onlyBW_clone->FixParameter(7, f1525Mass);

                    // onlyBW_clone->FixParameter(10, f1710Mass);
                    onlyBW_clone->FixParameter(11, f1710Width);

                    onlyBW->SetLineColor(4);
                    onlyBW->SetLineStyle(2);
                    // onlyBW->Draw("same");

                    // // Now plot the indivial resonances
                    TF1 *singlefits[4];
                    for (int i = 0; i < 4; i++)
                    {
                        singlefits[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, 1.00, 3.0, 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, 1.00, 3.0, 3); // Default
                        // singlefits[i] = new TF1(Form("singlef%d", i), single_BW, 1.00, 3.0, 3); // Constant width BW
                        singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                        singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                        singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                        singlefits[i]->SetLineColor(colors[i]);
                        singlefits[i]->SetLineStyle(2);
                        if (i == 3)
                            singlefits[i]->SetLineWidth(3);
                        singlefits[i]->Draw("same");
                    }

                    TLegend *ltemp = new TLegend(0.25, 0.52, 0.55, 0.87);
                    ltemp->SetFillStyle(0);
                    ltemp->SetBorderSize(0);
                    ltemp->SetTextFont(42);
                    ltemp->SetTextSize(0.03);
                    ltemp->AddEntry((TObject *)0, "", "");
                    ltemp->AddEntry((TObject *)0, "", "");
                    ltemp->AddEntry(hinvMass, "Data (stat. uncert.)", "lpe");
                    ltemp->AddEntry(BEexpol, "4rBW + Residual BG", "l");
                    ltemp->AddEntry(expol, "Residual BG", "l");
                    ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
                    ltemp->AddEntry(singlefits[1], "a_{2}(1320)^{0}", "l");
                    ltemp->AddEntry(singlefits[2], "f'_{2}(1525)", "l");
                    ltemp->AddEntry(singlefits[3], "f_{0}(1710)", "l");
#ifdef multiPanelPlots
                    if (ipt == 0)
                        ltemp->Draw("same");
#else
                    ltemp->Draw("same");
#endif

                    TLatex lat1;
                    lat1.SetNDC();
#ifdef singlePanelPlots
                    lat1.SetTextSize(0.03);
#else
                    lat1.SetTextSize(0.05);
#endif
                    lat1.SetTextFont(42);
#ifdef singlePanelPlots
                    lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{#it{s}} = 13.6 TeV");
                    lat1.DrawLatex(0.255, 0.85, Form("FT0M: %d-%d%%), |y|<0.5", multlow, multhigh));
#endif
                    lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", lowpT, highpT));

                    for (int i = 0; i < 4; i++)
                    {
                        auto lowLimit = (masses[i] - 2 * widths[i] > 1.0) ? masses[i] - 2 * widths[i] : 1.0;
                        auto highLimit = (masses[i] + 2 * widths[i] < 3.0) ? masses[i] + 2 * widths[i] : 3.0;
                        double significance_num = singlefits[i]->Integral(lowLimit, highLimit) / binwidthfile;
                        int binlow = hraw->GetXaxis()->FindBin(lowLimit);
                        int binhigh = hraw->GetXaxis()->FindBin(highLimit);
                        double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
                        double significance = significance_num / significance_den;
                        double signal_counts = significance_num;
                        double background_counts = hraw->Integral(binlow, binhigh) - signal_counts;

                        cout << "numerator " << significance_num << " denominator " << significance_den << endl;
                        cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                        double amplitude = obtained_parameters[3 * i];
                        double amplitude_err = BEexpol->GetParError(3 * i);
                        double statSignificance = amplitude / amplitude_err;
                        cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
                        hSignificance[i]->SetBinContent(ipt + 1, significance);
                        hSignificance[i]->SetBinError(ipt + 1, 0);
                        hStatSignificance[i]->SetBinContent(ipt + 1, statSignificance);
                        hStatSignificance[i]->SetBinError(ipt + 1, 0);
                    }
#endif

                    // // // //******************************************************************************************
                    // // // //********************************* common for all fits ***************************************
                    gPad->Update();
                    TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
                    if (ptstats != nullptr)
                    {
                        // Customize precision for fit parameters display
                        TList *listOfLines = ptstats->GetListOfLines();
                        // You can iterate through and modify specific lines if needed

                        ptstats->SetX1NDC(0.6);
                        ptstats->SetX2NDC(0.99);
                        ptstats->SetY1NDC(0.4);
                        ptstats->SetY2NDC(0.92);
                        ptstats->Draw("same");
                    }
// c->SaveAs((savepath + Form("/rBWfit_pt_%.2f_%.2f_%s.png", lowpT, highpT, kvariation1.c_str())).c_str());
// c->SaveAs((savepath + "/rBWfit.png").c_str());
// c->SaveAs((savepath + Form("/rBWfit_fit_%.2f_%.2f.png", fitlow, fithigh)).c_str());
#ifdef singlePanelPlots
                    c->SaveAs((savepath_mult + Form("/rBWfit_pt_%.1f_%.1f.png", lowpT, highpT)).c_str());
#endif

                    double fitnorm1525 = BEexpol->GetParameter(6);
                    double fitnorm1525_err = BEexpol->GetParError(6);
                    double fitnorm1710 = BEexpol->GetParameter(9);
                    double fitnorm1710_err = BEexpol->GetParError(9);
                    double fitmass1525 = BEexpol->GetParameter(7);
                    double fitmass1525_err = BEexpol->GetParError(7);
                    double fitmass1710 = BEexpol->GetParameter(10);
                    double fitmass1710_err = BEexpol->GetParError(10);
                    double fitwidth1525 = BEexpol->GetParameter(8);
                    double fitwidth1525_err = BEexpol->GetParError(8);
                    double fitwidth1710 = BEexpol->GetParameter(11);
                    double fitwidth1710_err = BEexpol->GetParError(11);
                    double fitnorm1270 = BEexpol->GetParameter(0);
                    double fitnorm1270_err = BEexpol->GetParError(0);
                    double fitmass1270 = BEexpol->GetParameter(1);
                    double fitmass1270_err = BEexpol->GetParError(1);
                    double fitwidth1270 = BEexpol->GetParameter(2);
                    double fitwidth1270_err = BEexpol->GetParError(2);
                    double fitnorm1320 = BEexpol->GetParameter(3);
                    double fitnorm1320_err = BEexpol->GetParError(3);
                    double fitmass1320 = BEexpol->GetParameter(4);
                    double fitmass1320_err = BEexpol->GetParError(4);
                    double fitwidth1320 = BEexpol->GetParameter(5);
                    double fitwidth1320_err = BEexpol->GetParError(5);
                    double resdual_bkg_par1 = BEexpol->GetParameter(12);
                    double resdual_bkg_par2 = BEexpol->GetParameter(13);
                    double resdual_bkg_par3 = BEexpol->GetParameter(14);
                    double resdual_bkg_par4 = BEexpol->GetParameter(15);

                    if (fitnorm1270 < 0 || fitnorm1320 < 0 || fitnorm1525 < 0 || fitnorm1710 < 0)
                    {
                        cout << "Fit failed, amplitues are negative " << endl;
                        for (int i = 0; i < 10; i++)
                        {
                            cout << "Fit failed!!!!!!!!!!!!!!! " << endl;
                        }
                    }

#ifdef residual_subtracted
// Now subtract the residual background and plot
#ifdef singlePanelPlots
                    TCanvas *c2 = new TCanvas("", "", 720, 720);
                    SetCanvasStyle(c2, 0.14, 0.03, 0.05, 0.14);
#endif
#ifdef multiPanelPlots
                    cMultiPanelResidual->cd(ipt + 1);
                    gPad->SetTopMargin(0.06);
                    gPad->SetBottomMargin(0.14);
#endif
                    expol_clone->SetRange(0.99, 2.99);
                    hsubtracted->Add(expol_clone, -1);
                    hsubtracted->GetXaxis()->SetRangeUser(BEexpol->GetXmin(), BEexpol->GetXmax());
                    hsubtracted->SetMaximum(hsubtracted->GetMaximum() * resMaximumFactor);
                    hsubtracted->SetMinimum(-hsubtracted->GetMaximum() * 0.05);
                    hsubtracted->GetYaxis()->SetTitleOffset(1.0);
                    hsubtracted->Draw();
                    onlyBW_clone->SetNpx(1000);
                    hsubtracted->Fit("onlyBW_clone", "REBMS");
                    // onlyBW_clone->Draw("same");
                    double *obtained_parameters2 = onlyBW_clone->GetParameters();
                    TLine *line = new TLine(BEexpol->GetXmin() + 0.01, 0, BEexpol->GetXmax() - 0.01, 0);
                    line->SetLineColor(1);
                    line->SetLineStyle(4);
                    line->Draw("same");

                    // // Now plot the indivial resonances
                    TF1 *singlefits1[4];
                    for (int i = 0; i < 4; i++)
                    {
#ifdef b_massdepWidth_modifiedBoltzmann
                        singlefits1[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, 1.00, 2.5, 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, 1.00, 2.5, 3);
#else
                        singlefits1[i] = new TF1(Form("singlef%d", i), single_BW, onlyBW_clone->GetXmin(), onlyBW_clone->GetXmax(), 3);
#endif
                        singlefits1[i]->SetParameter(0, obtained_parameters2[3 * i]);
                        singlefits1[i]->SetParameter(1, obtained_parameters2[3 * i + 1]);
                        singlefits1[i]->SetParameter(2, obtained_parameters2[3 * i + 2]);
                        singlefits1[i]->SetLineColor(colors[i]);
                        singlefits1[i]->SetLineStyle(2);
                        singlefits1[i]->Draw("same");
                    }

                    TLegend *ltemp2 = new TLegend(0.20, 0.55, 0.42, 0.80);
                    ltemp2->SetFillStyle(0);
                    ltemp2->SetTextFont(42);
                    ltemp2->SetBorderSize(0);
                    ltemp2->SetTextSize(0.03);
                    // ltemp2->SetHeader("Residual bkg subtraction");
                    ltemp2->AddEntry(hsubtracted, "Data (stat. uncert.)", "lpe");
                    ltemp2->AddEntry(onlyBW_clone, "4rBW fit", "l");
                    ltemp2->AddEntry(singlefits1[0], "f_{2}(1270)", "l");
                    ltemp2->AddEntry(singlefits1[1], "a_{2}(1320)^{0}", "l");
                    ltemp2->AddEntry(singlefits1[2], "f'_{2}(1525)", "l");
                    ltemp2->AddEntry(singlefits1[3], "f_{0}(1710)", "l");
                    TLatex lat2;
                    lat2.SetNDC();
#ifdef multiPanelPlots
                    lat2.SetTextSize(0.05);
                    if (ipt == 0)
                        ltemp2->Draw("same");
#else
                    lat2.SetTextSize(0.03);
                    ltemp2->Draw("same");
#endif
                    lat2.SetTextFont(42);
#ifdef singlePanelPlots
                    lat2.DrawLatex(0.215, 0.90, "Residual BG subtracted");
                    lat2.DrawLatex(0.215, 0.86, Form("FT0M: %d-%d%%), |y|<0.5", multlow, multhigh));
#endif
                    lat2.DrawLatex(0.215, 0.82, Form("%.1f < p_{T} < %.1f GeV/c", lowpT, highpT));

// c2->SaveAs((savepath + "/rBWfit_residual_" + kvariation1 + ".png").c_str());
#ifdef singlePanelPlots
                    c2->SaveAs((savepath_mult + Form("/rBWfit_residual_pt_%.1f_%.1f.png", lowpT, highpT)).c_str());
#ifdef before_combinatorialPlots
                    TCanvas *cSigWithBkg = new TCanvas("", "", 720, 720);
                    SetCanvasStyle(cSigWithBkg, 0.14, 0.03, 0.05, 0.14);
                    cSigWithBkg->cd();
                    SetHistoQA(hraw);
                    SetHistoQA(hbkg);
                    hraw->GetXaxis()->SetRangeUser(1.01, 2.20);
                    hraw->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
                    hraw->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidthfile * 1000));
                    hraw->SetMaximum(hraw->GetMaximum() * 1.2);
                    hraw->SetMarkerSize(1.0);
                    hraw->Draw("pe");
                    hbkg->SetMarkerColor(kRed);
                    hbkg->SetLineColor(kRed);
                    hbkg->SetMarkerSize(1.0);
                    hbkg->Draw("pe,same");
                    TLegend *lSigWithBkg = new TLegend(0.55, 0.56, 0.85, 0.72);
                    lSigWithBkg->SetFillStyle(0);
                    lSigWithBkg->SetBorderSize(0);
                    lSigWithBkg->SetTextFont(42);
                    lSigWithBkg->SetTextSize(0.03);
                    lSigWithBkg->AddEntry(hraw, "Same event K^{0}_{s}K^{0}_{s} pair", "lpe");
                    lSigWithBkg->AddEntry(hbkg, "Same event rotational bkg", "lpe");
                    lSigWithBkg->Draw("same");

                    lat1.DrawLatex(0.58, 0.89, "ALICE");
                    lat1.DrawLatex(0.58, 0.835, "pp, #sqrt{#it{s}} = 13.6 TeV");
                    lat1.DrawLatex(0.58, 0.785, Form("FT0M: %d-%d%%), |y|<0.5", multlow, multhigh));
                    lat1.DrawLatex(0.58, 0.735, Form("%.1f < p_{T} < %.1f GeV/c", lowpT, highpT));
                    cSigWithBkg->SaveAs((savepath_mult + Form("/SignalWithBkg_pt_%.1f_%.1f.png", lowpT, highpT)).c_str());

                    TCanvas *cSignalWithoutFit = new TCanvas("", "", 720, 720);
                    SetCanvasStyle(cSignalWithoutFit, 0.14, 0.03, 0.05, 0.14);
                    cSignalWithoutFit->cd();
                    SetHistoQA(hinvMassClone);
                    hinvMassClone->GetXaxis()->SetRangeUser(1.05, 2.20);
                    hinvMassClone->Draw("pe");
                    lat1.DrawLatex(0.58, 0.89, "ALICE");
                    lat1.DrawLatex(0.58, 0.835, "pp, #sqrt{#it{s}} = 13.6 TeV");
                    lat1.DrawLatex(0.58, 0.785, Form("FT0M: %d-%d%%), |y|<0.5", multlow, multhigh));
                    lat1.DrawLatex(0.58, 0.735, Form("%.1f < p_{T} < %.1f GeV/c", lowpT, highpT));
                    cSignalWithoutFit->SaveAs((savepath_mult + Form("/SignalWithoutFit_pt_%.1f_%.1f.png", lowpT, highpT)).c_str());
#endif
#endif
#endif

                    float ptBinWidth = highpT - lowpT;
                    double fitrangelow = 1.001; // to avoid mass thresold (in denominator of BW function there is x^2 - m^2).
                    // double fitrangehigh = 2.499;
                    // double fitrangelow = BEexpol->GetXmin();
                    // double fitrangehigh = BEexpol->GetXmax();

                    TMatrixDSym cov = fitResultptr->GetCovarianceMatrix();
                    TMatrixDSym cov1;
                    TMatrixDSym cov2;
                    cov.GetSub(0, 11, 0, 11, cov1);
                    cov.GetSub(12, 15, 12, 15, cov2);
                    Double_t *b = cov1.GetMatrixArray();
                    Double_t *a = cov2.GetMatrixArray();
                    Double_t *para = onlyBW_clone->GetParameters();

                    float nsigma_yield = 5.0;
                    double fitRegion1270_low = f1270Mass - nsigma_yield * f1270Width;
                    double fitRegion1270_high = f1270Mass + nsigma_yield * f1270Width;
                    double fitRegion1320_low = a1320Mass - nsigma_yield * a1320Width;
                    double fitRegion1320_high = a1320Mass + nsigma_yield * a1320Width;
                    double fitRegion1525_low = f1525Mass - nsigma_yield * f1525Width;
                    double fitRegion1525_high = f1525Mass + nsigma_yield * f1525Width;
                    double fitRegion1710_low = f1710Mass - nsigma_yield * f1710Width;
                    double fitRegion1710_high = f1710Mass + nsigma_yield * f1710Width;
                    if (fitRegion1270_low < fitrangelow)
                        fitRegion1270_low = fitrangelow;
                    if (fitRegion1320_low < fitrangelow)
                        fitRegion1320_low = fitrangelow;
                    if (fitRegion1525_low < fitrangelow)
                        fitRegion1525_low = fitrangelow;
                    if (fitRegion1710_low < fitrangelow)
                        fitRegion1710_low = fitrangelow;
                    // if (fitRegion1270_high > fitrangehigh)
                    //     fitRegion1270_high = fitrangehigh;
                    // if (fitRegion1320_high > fitrangehigh)
                    //     fitRegion1320_high = fitrangehigh;
                    // if (fitRegion1525_high > fitrangehigh)
                    //     fitRegion1525_high = fitrangehigh;
                    // if (fitRegion1710_high > fitrangehigh)
                    //     fitRegion1710_high = fitrangehigh;

                    // // Yield calculation
                    double yield1270Integral = singlefits[0]->Integral(fitRegion1270_low, fitRegion1270_high) / (ptBinWidth * binwidthfile * total_events);
                    double yield1320Integral = singlefits[1]->Integral(fitRegion1320_low, fitRegion1320_high) / (ptBinWidth * binwidthfile * total_events);
                    double yield1525Integral = singlefits[2]->Integral(fitRegion1525_low, fitRegion1525_high) / (ptBinWidth * binwidthfile * total_events);
                    double yield1710Integral = singlefits[3]->Integral(fitRegion1710_low, fitRegion1710_high) / (ptBinWidth * binwidthfile * total_events);

                    // Create individual covariance matrices for each resonance (3x3 each)
                    TMatrixDSym cov1270(3), cov1320(3), cov1525(3), cov1710(3);
                    cov.GetSub(0, 2, 0, 2, cov1270);   // f1270: parameters 0-2
                    cov.GetSub(3, 5, 3, 5, cov1320);   // a1320: parameters 3-5
                    cov.GetSub(6, 8, 6, 8, cov1525);   // f1525: parameters 6-8
                    cov.GetSub(9, 11, 9, 11, cov1710); // f1710: parameters 9-11

                    Double_t *cov1270_array = cov1270.GetMatrixArray();
                    Double_t *cov1320_array = cov1320.GetMatrixArray();
                    Double_t *cov1525_array = cov1525.GetMatrixArray();
                    Double_t *cov1710_array = cov1710.GetMatrixArray();

                    Double_t *para1270 = singlefits[0]->GetParameters();
                    Double_t *para1320 = singlefits[1]->GetParameters();
                    Double_t *para1525 = singlefits[2]->GetParameters();
                    Double_t *para1710 = singlefits[3]->GetParameters();

                    double yield1270_err = singlefits[0]->IntegralError(fitRegion1270_low, fitRegion1270_high, para1270, cov1270_array) / (ptBinWidth * binwidthfile * total_events);
                    double yield1320_err = singlefits[1]->IntegralError(fitRegion1320_low, fitRegion1320_high, para1320, cov1320_array) / (ptBinWidth * binwidthfile * total_events);
                    double yield1525_err = singlefits[2]->IntegralError(fitRegion1525_low, fitRegion1525_high, para1525, cov1525_array) / (ptBinWidth * binwidthfile * total_events);
                    double yield1710_err = singlefits[3]->IntegralError(fitRegion1710_low, fitRegion1710_high, para1710, cov1710_array) / (ptBinWidth * binwidthfile * total_events);

                    cout << "Yield 1270: " << yield1270Integral << " +- " << yield1270_err << endl;
                    cout << "Yield 1320: " << yield1320Integral << " +- " << yield1320_err << endl;
                    cout << "Yield 1525: " << yield1525Integral << " +- " << yield1525_err << endl;
                    cout << "Yield 1710: " << yield1710Integral << " +- " << yield1710_err << endl;

                    hChi2NDF->SetBinContent(ipt + 1, BEexpol->GetChisquare() / BEexpol->GetNDF());
                    hChi2NDF->SetBinError(ipt + 1, 0);
                    for (int ires = 0; ires < 4; ires++)
                    {
                        double mass_temp = 0.0;
                        double massErr_temp = 0.0;
                        double width_temp = 0.0;
                        double widthErr_temp = 0.0;
                        double yield_temp = 0.0;
                        double yieldErr_temp = 0.0;

                        if (ires == 0)
                        {
                            mass_temp = fitmass1270;
                            massErr_temp = fitmass1270_err;
                            width_temp = fitwidth1270;
                            widthErr_temp = fitwidth1270_err;
                            yield_temp = yield1270Integral;
                            yieldErr_temp = yield1270_err;
                        }
                        else if (ires == 1)
                        {
                            mass_temp = fitmass1320;
                            massErr_temp = fitmass1320_err;
                            width_temp = fitwidth1320;
                            widthErr_temp = fitwidth1320_err;
                            yield_temp = yield1320Integral;
                            yieldErr_temp = yield1320_err;
                        }
                        else if (ires == 2)
                        {
                            mass_temp = fitmass1525;
                            massErr_temp = fitmass1525_err;
                            width_temp = fitwidth1525;
                            widthErr_temp = fitwidth1525_err;
                            yield_temp = yield1525Integral;
                            yieldErr_temp = yield1525_err;
                        }
                        else
                        {
                            mass_temp = fitmass1710;
                            massErr_temp = fitmass1710_err;
                            width_temp = fitwidth1710;
                            widthErr_temp = fitwidth1710_err;
                            yield_temp = yield1710Integral;
                            yieldErr_temp = yield1710_err;
                        }

                        hMass[ires]->SetBinContent(ipt + 1, mass_temp);
                        hMass[ires]->SetBinError(ipt + 1, massErr_temp);
                        hWidth[ires]->SetBinContent(ipt + 1, width_temp);
                        hWidth[ires]->SetBinError(ipt + 1, widthErr_temp);
                        hYield[ires]->SetBinContent(ipt + 1, yield_temp);
                        hYield[ires]->SetBinError(ipt + 1, yieldErr_temp);
                    }

                    // // // for .txt files for systematics
                    // file << "Norm range " << kNormRangepT[0][0] << " - " << kNormRangepT[0][1] << endl;
                    // file << "Fit range " << BEexpol->GetXmin() << " - " << BEexpol->GetXmax() << endl;
                    // file << "Significance " << significance << endl;
                    // file << "StatSignificance " << statSignificance << endl;
                    // file << "Fit parameters of f1710 " << endl;
                    // file << "Chi2NDF " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
                    // file << yield1710Integral << " ± " << yield1710_err << endl;
                    // file << yield1710Integral << " ± " << yield1710_err << endl;
                    // file << fitmass1710 << " ± " << fitmass1710_err << endl;
                    // file << fitwidth1710 << " ± " << fitwidth1710_err << endl;
                    // file << "for f1525" << endl;
                    // file << yield1525Integral << " ± " << yield1525_err << endl;
                    // file << yield1525Integral << " ± " << yield1525_err << endl;
                    // file << fitmass1525 << " ± " << fitmass1525_err << endl;
                    // file << fitwidth1525 << " ± " << fitwidth1525_err << endl;
                    // file << "for f1270" << endl;
                    // file << yield1270Integral << " ± " << yield1270_err << endl;
                    // file << yield1270Integral << " ± " << yield1270_err << endl;
                    // file << fitmass1270 << " ± " << fitmass1270_err << endl;
                    // file << fitwidth1270 << " ± " << fitwidth1270_err << endl;
                    // file << "for a1320" << endl;
                    // file << yield1320Integral << " ± " << yield1320_err << endl;
                    // file << yield1320Integral << " ± " << yield1320_err << endl;
                    // file << fitmass1320 << " ± " << fitmass1320_err << endl;
                    // file << fitwidth1320 << " ± " << fitwidth1320_err << endl;

#ifdef doublepanelplot
                    // gStyle->SetOptFit(0);
                    TCanvas *c1 = new TCanvas("", "", 720, 720);
                    SetCanvasStyle(c1, 0.14, 0.03, 0.05, 0.14);
                    double pad1Size, pad2Size;
                    canvas_style(c1, pad1Size, pad2Size);
                    c1->cd(1);
                    gPad->SetFrameLineWidth(2);
                    gPad->SetLeftMargin(0.14);
                    gPad->SetTickx(1);
                    SetHistoQA(hinvMass);
                    SetHistoQA(hsubtracted);
                    hinvMass->GetXaxis()->SetTitleSize(0.04 / pad1Size);
                    hinvMass->GetYaxis()->SetTitleSize(0.04 / pad1Size);
                    hinvMass->GetXaxis()->SetLabelSize(0.04 / pad1Size);
                    hinvMass->GetYaxis()->SetLabelSize(0.04 / pad1Size);
                    hinvMass->GetXaxis()->SetTitleOffset(1.02);
                    hinvMass->GetYaxis()->SetTitleOffset(0.87);
                    hinvMass->GetYaxis()->CenterTitle(0);
                    hinvMass->GetXaxis()->SetRangeUser(1.05, 2.15);
                    hsubtracted->GetXaxis()->SetRangeUser(BEexpol->GetXmin(), BEexpol->GetXmax());
                    hinvMass->SetMarkerSize(1.0);
                    hinvMass->SetLineColor(1);
                    hinvMass->SetMarkerColor(1);
                    hinvMass->SetStats(0);
                    // hinvMass->SetMinimum(-40000);
                    // hinvMass->SetMinimum(-10000);
                    hinvMass->SetMinimum(-4000);
                    // hinvMass->SetMaximum(0.9e6);
                    hinvMass->SetMaximum(hinvMass->GetMaximum() * 0.9);
                    hinvMass->GetYaxis()->SetMaxDigits(4);
                    hinvMass->GetYaxis()->SetNdivisions(506);
                    hinvMass->GetYaxis()->SetTitle(0);
                    hinvMass->Draw("pe");
                    expol->Draw("same");
                    onlyBW->SetLineWidth(0);
                    onlyBW->SetFillColor(4);
                    onlyBW->SetFillStyle(1001);
                    onlyBW->Draw("same");
                    onlyBW->SetNpx(1000);
                    gPad->Update();

                    TArrow *arrow = new TArrow(0.3300313, 0.4925806, 0.3300313, 0.4137097, 0.02, "|>");
                    arrow->SetFillColor(1);
                    arrow->SetFillStyle(1001);
                    arrow->SetLineWidth(2);
                    arrow->SetNDC();
                    arrow->Draw();
                    TLatex *tex = new TLatex(0.3347237, 0.5409677, "f_{2}(1270)/a_{2}^{0}(1320)");
                    tex->SetNDC();
                    tex->SetTextAlign(22);
                    tex->SetLineWidth(2);
                    tex->Draw();
                    arrow = new TArrow(0.4996435, 0.4056322, 0.4996435, 0.3264368, 0.02, "|>");
                    arrow->SetFillColor(1);
                    arrow->SetFillStyle(1001);
                    arrow->SetLineWidth(2);
                    arrow->SetNDC();
                    arrow->Draw();
                    tex = new TLatex(0.5, 0.4458621, "f^{,}_{2}(1525)");
                    tex->SetNDC();
                    tex->SetTextAlign(22);
                    tex->SetLineWidth(2);
                    tex->Draw();
                    arrow = new TArrow(0.6328134, 0.271954, 0.6328134, 0.1927586, 0.02, "|>");
                    arrow->SetFillColor(1);
                    arrow->SetFillStyle(1001);
                    arrow->SetLineWidth(2);
                    arrow->SetNDC();
                    arrow->Draw();
                    tex = new TLatex(0.6397772, 0.3093103, "f_{0}(1710)");
                    tex->SetNDC();
                    tex->SetTextAlign(22);
                    tex->SetLineWidth(2);
                    tex->Draw();

                    // TLegend *leg = new TLegend(0.65, 0.47, 0.99, 0.77);
                    TLegend *leg = new TLegend(0.65, 0.44, 0.99, 0.86);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetTextSize(0.055);
                    leg->SetBorderSize(0);
                    // leg->SetHeader("ALICE Performance");
                    leg->AddEntry(hinvMass, "Data (stat. uncert.)", "p");
                    leg->AddEntry(BEexpol, "4 BW + Res. Bkg", "l");
                    leg->AddEntry(expol, "Res. Bkg", "l");
                    leg->AddEntry(onlyBW, "Signal", "f");
                    leg->Draw("same");

                    TLatex lat5;
                    lat5.SetNDC();
                    lat5.SetTextSize(0.055);
                    lat5.SetTextFont(42);
                    // lat5.DrawLatex(0.32, 0.80, "ALICE Performance");
                    lat5.DrawLatex(0.32, 0.82, "pp, #sqrt{#it{s}} = 13.6 TeV");
                    lat5.DrawLatex(0.32, 0.73, "FT0M: 0-100%, |#it{y}| < 0.5");
                    lat5.DrawLatex(0.32, 0.65, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", lowpT, highpT));

                    // TLatex *text4 = new TLatex(0.65, 0.80, "ALICE, work in progress");
                    // text4->SetNDC();
                    // text4->SetTextSize(0.06);
                    // text4->SetTextFont(42);
                    // text4->Draw("same");

                    c1->cd(2);
                    gPad->SetFrameLineWidth(2);
                    gPad->SetLeftMargin(0.14);
                    gPad->SetTickx(1);
                    hsubtracted->GetYaxis()->SetTitleSize(0.04 / pad2Size);
                    hsubtracted->GetXaxis()->SetTitleSize(0.04 / pad2Size);
                    hsubtracted->GetXaxis()->SetLabelSize(0.04 / pad2Size);
                    hsubtracted->GetYaxis()->SetLabelSize(0.04 / pad2Size);
                    hsubtracted->GetYaxis()->SetNdivisions(505);
                    hsubtracted->GetXaxis()->SetRangeUser(1.05, 2.15);
                    hsubtracted->SetStats(0);
                    hsubtracted->SetMinimum(0);
                    // hsubtracted->SetMaximum(0.13e6);
                    hsubtracted->SetMaximum(hsubtracted->GetMaximum() * 0.62);
                    // hsubtracted->SetMaximum(hsubtracted->GetMaximum() * 0.67);
                    hsubtracted->SetMarkerSize(1.05);
                    // hsubtracted->GetYaxis()->SetMaxDigits(10);
                    hsubtracted->SetMarkerStyle(53);
                    hsubtracted->Draw("pe");

                    TLatex lat6;
                    lat6.SetNDC();
                    lat6.SetTextSize(0.06);
                    lat6.SetTextFont(42);
                    lat6.DrawLatex(0.59, 0.86, "Residual background subtraction");
                    // lat6.DrawLatex(0.65, 0.70, "pp #sqrt{#it{s}} = 13.6 TeV");
                    // lat6.DrawLatex(0.65, 0.60, "FT0M: 0-100%, |#it{y}|<0.5");
                    // lat6.DrawLatex(0.65, 0.50, Form("%.1f < #it{p}_{T} #leq %.1f GeV/c", lowpT, highpT));

                    // TLegend *leg = new TLegend(0.65, 0.50, 0.99, 0.83);
                    // leg->SetFillStyle(0);
                    // leg->SetTextFont(42);
                    // leg->SetTextSize(0.06);
                    // leg->SetBorderSize(0);
                    // leg->AddEntry(hinvMass, "Data (stat. uncert.)", "p");
                    // leg->AddEntry(BEexpol, "4 BW + Residual BG", "l");
                    // leg->AddEntry(onlyBW, "Signal", "f");
                    // leg->AddEntry(expol, "Residual BG", "l");
                    // leg->Draw("same");

                    TLatex *text5 = new TLatex(0.16, 0.89, "#times 10^{3}");
                    text5->SetNDC();
                    text5->SetTextSize(0.08);
                    text5->SetTextFont(42);
                    text5->Draw("same");

                    int linestyles[4] = {2, 2, 2, 7};
                    int linewidths[] = {2, 2, 2, 3};

                    for (int i = 0; i < 4; i++)
                    {
                        // singlefits1[i]->SetRange(masses[i] - 1 * widths[i], masses[i] + 1 * widths[i]);
                        singlefits1[i]->SetLineStyle(linestyles[i]);
                        singlefits1[i]->SetLineWidth(linewidths[i]);
                        singlefits1[i]->Draw("same");
                    }

                    TLegend *leg2 = new TLegend(0.63, 0.45, 0.99, 0.85);
                    leg2->SetFillStyle(0);
                    leg2->SetTextFont(42);
                    leg2->SetTextSize(0.055);
                    leg2->SetBorderSize(0);
                    leg2->AddEntry(singlefits1[0], "f_{2}(1270)", "l");
                    leg2->AddEntry(singlefits1[1], "a_{2}(1320)^{0}", "l");
                    leg2->AddEntry(singlefits1[2], "f'_{2}(1525)", "l");
                    leg2->AddEntry(singlefits1[3], "f_{0}(1710)", "l");
                    leg2->Draw("same");

                    c1->cd();
                    TPad *textPad1 = new TPad("textPad1", "", 0, 0, 1, 1);
                    textPad1->SetFillStyle(0); // Transparent
                    textPad1->SetFrameFillStyle(0);
                    textPad1->Draw();
                    textPad1->cd();

                    // Add left-side text
                    TLatex *textLeft1 = new TLatex(0.025, 0.4, Form("Counts / (%.0f MeV/#it{c}^{2})", binwidthfile * 1000));
                    textLeft1->SetTextAlign(12); // Left alignment
                    textLeft1->SetTextAngle(90); // Rotate the text by 90 degrees
                    textLeft1->SetTextSize(0.045);
                    textLeft1->SetTextFont(42);
                    textLeft1->Draw();

                    // c1->SaveAs("/home/sawan/Music/r4BWfit_doublepanel.png");
                    c1->SaveAs((savepath_mult + Form("/rBWfit_doublepanel_pt_%.1f_%.1f.png", lowpT, highpT)).c_str());

#endif

#ifdef multiPanelPlots
                    cMultiPanelWithBkg->cd(ipt + 1);
                    gPad->SetTopMargin(0.06);
                    gPad->SetBottomMargin(0.14);
                    gPad->SetLeftMargin(0.13);
                    gPad->SetRightMargin(0.05);
                    // SetHistoQA(hraw);
                    // SetHistoQA(hbkg);
                    hraw->GetYaxis()->SetTitleOffset(1.2);
                    hraw->Draw("pe");
                    hbkg->Draw("pe same");
                    if (ipt == 0)
                    {
                        TLegend *ltemp3 = new TLegend(0.5, 0.55, 0.92, 0.80);
                        ltemp3->SetFillStyle(0);
                        ltemp3->SetTextFont(42);
                        ltemp3->SetBorderSize(0);
                        ltemp3->SetTextSize(0.04);
                        // ltemp3->SetHeader("With background");
                        ltemp3->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV, FT0M: 0-100%");
                        ltemp3->AddEntry(hraw, "Same event K^{0}_{s}K^{0}_{s} pair", "p");
                        ltemp3->AddEntry(hbkg, "Same event rotational bkg", "p");
                        ltemp3->Draw("same");
                    }
                    lat2.DrawLatex(0.5, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", lowpT, highpT));
#endif

                } // =========End of fit range loop=============
                // file.close();

            } // ======End of pT loop============
            fDefault->cd();
            TDirectory *multDir = fDefault->mkdir(Form("Mult_%d_%d", multlow, multhigh), Form("Mult_%d_%d", multlow, multhigh));
            multDir->cd();
            hChi2NDF->Write();
            for (int ires = 0; ires < 4; ires++)
            {
                hMass[ires]->Write(Form("hMass_%s", resonance_mass[ires].c_str()));
                hWidth[ires]->Write(Form("hWidth_%s", resonance_mass[ires].c_str()));
                hYield[ires]->Write(Form("hYield_%s", resonance_mass[ires].c_str()));
                hSignificance[ires]->Write(Form("hSignificance_%s", resonance_mass[ires].c_str()));
                hStatSignificance[ires]->Write(Form("hStatSignificance_%s", resonance_mass[ires].c_str()));
            }
#ifdef multiPanelPlots
            cMultiPanelResidual->SaveAs((savepath_mult + Form("/rBWfit_residuals_multpanel_%s.png", kvariation1.c_str())).c_str());
            cMultiPanelFit->SaveAs((savepath_mult + Form("/rBWfit_fits_multpanel_%s.png", kvariation1.c_str())).c_str());
            cMultiPanelWithBkg->SaveAs((savepath_mult + Form("/rBWfit_withbkg_multpanel_%s.png", kvariation1.c_str())).c_str());
#endif
        } // end of multiplicity loop

    } // End of systematics loop
}
// end of main program

// *****************************************************************************************************
//************************************Fit functions************************************************* */

// We will define the single BW and the sum of 4 BWs
Double_t single_BW_hera(double *x, double *par)
{
    // normalization factor is missing and how to add it I am not sure
    double amplitude = par[0];
    double mass = par[1];
    double width = par[2];

    double den = (x[0] * x[0] - mass * mass) * (x[0] * x[0] - mass * mass) + mass * mass * width * width;
    double realnum = (mass * mass - x[0] * x[0]) * mass * TMath::Sqrt(width);
    double imagnum = mass * mass * width * TMath::Sqrt(width);

    double real3BW = realnum / den;
    double imag3BW = imagnum / den;
    double sig1 = amplitude * (real3BW * real3BW + imag3BW * imag3BW);

    return sig1;
}

Double_t single_BW_mass_dep_spin0(double *x, double *par)
{
    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    int spin = 0;
    double n1 = (2.0 * spin + 1.0) / 2.0;

    double yield = par[0];
    double mass = par[1];
    double width = par[2] * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));

    return fit;
}

Double_t single_BW_mass_dep_spin2(double *x, double *par)
{
    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    int spin = 2;
    double n1 = (2.0 * spin + 1.0) / 2.0;

    double yield = par[0];
    double mass = par[1];
    double width = par[2] * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));

    return fit;
}

Double_t BWsum_hera(double *x, double *par) // taken 4 resonances here
{
    // total 10 parameters, 2 for each resonance (total 4 resonance), 1 for normalization of spin2 particles, 1 for normalization of spin0 particle
    double mass1270 = par[0];
    double width1270 = par[1];
    double mass1320 = par[2];
    double width1320 = par[3];
    double mass1525 = par[4];
    double width1525 = par[5];
    double mass1710 = par[6];
    double width1710 = par[7];
    double a0 = par[8];
    double a1 = par[9];

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1320 = (x[0] * x[0] - mass1320 * mass1320) * (x[0] * x[0] - mass1320 * mass1320) + mass1320 * mass1320 * width1320 * width1320;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270);
    double realnum1320 = (mass1320 * mass1320 - x[0] * x[0]) * mass1320 * TMath::Sqrt(width1320);
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525);
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710);

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270);
    double imagnum1320 = mass1320 * mass1320 * width1320 * TMath::Sqrt(width1320);
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525);
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710);

    double real3BW = 5 * realnum1270 / den1270 - 3 * realnum1320 / den1320 + 2 * realnum1525 / den1525;

    double imag3BW = 5 * imagnum1270 / den1270 - 3 * imagnum1320 / den1320 + 2 * imagnum1525 / den1525;

    double sig1 = (real3BW * real3BW + imag3BW * imag3BW);

    double fit = a0 * sig1 + a1 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / (den1710 * den1710);

    return fit;
}

Double_t BWsum_hera_mass_dep(double *x, double *par) // taken 4 resonances here
{
    // total 10 parameters, 2 for each resonance (total 4 resonance), 1 for normalization of spin2 particles, 1 for normalization of spin0 particle
    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[0] * par[0] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[2] * par[2] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart4 = par[6] * par[6] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    double mass1270 = par[0];
    double width1270 = par[1] * (TMath::Power(par[0] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double mass1320 = par[2];
    double width1320 = par[3] * (TMath::Power(par[2] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double mass1525 = par[4];
    double width1525 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n1);
    double mass1710 = par[6];
    double width1710 = par[7] * (TMath::Power(par[6] / x[0], 1.0)) * TMath::Power((npart1) / (dpart4), n2);

    // double mass1270 = par[0];
    // double width1270 = par[1];
    // double mass1320 = par[2];
    // double width1320 = par[3];
    // double mass1525 = par[4];
    // double width1525 = par[5];
    // double mass1710 = par[6];
    // double width1710 = par[7];
    double a0 = par[8];
    double a1 = par[9];
    // double norm1 = par[10];
    // double norm2 = par[11];
    // double norm3 = par[12];
    double norm1 = 5;
    double norm2 = -3;
    double norm3 = 2;

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1320 = (x[0] * x[0] - mass1320 * mass1320) * (x[0] * x[0] - mass1320 * mass1320) + mass1320 * mass1320 * width1320 * width1320;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270);
    double realnum1320 = (mass1320 * mass1320 - x[0] * x[0]) * mass1320 * TMath::Sqrt(width1320);
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525);
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710);

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270);
    double imagnum1320 = mass1320 * mass1320 * width1320 * TMath::Sqrt(width1320);
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525);
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710);

    double real3BW = norm1 * realnum1270 / den1270 + norm2 * realnum1320 / den1320 + norm3 * realnum1525 / den1525;

    double imag3BW = norm1 * imagnum1270 / den1270 + norm2 * imagnum1320 / den1320 + norm3 * imagnum1525 / den1525;

    double sig1 = (real3BW * real3BW + imag3BW * imag3BW);

    double fit = a0 * sig1 + a1 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / (den1710 * den1710);

    return fit;
}

Double_t BWsum_hera_const(double *x, double *par) // taken 4 resonances here
{
    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
    double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    double yield1270 = par[0];
    double mass1270 = par[1];
    // double width1270 = par[2];
    double width1270 = par[2] * (TMath::Power(par[1] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double yield1320 = par[3];
    double mass1320 = par[4];
    // double width1320 = par[5];
    double width1320 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double yield1525 = par[6];
    double mass1525 = par[7];
    // double width1525 = par[8];
    double width1525 = par[8] * (TMath::Power(par[7] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n1);
    double yield1710 = par[9];
    double mass1710 = par[10];
    // double width1710 = par[11];
    double width1710 = par[11] * (TMath::Power(par[10] / x[0], 1.0)) * TMath::Power((npart1) / (dpart4), n2);

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1320 = (x[0] * x[0] - mass1320 * mass1320) * (x[0] * x[0] - mass1320 * mass1320) + mass1320 * mass1320 * width1320 * width1320;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270) / den1270;
    double realnum1320 = (mass1320 * mass1320 - x[0] * x[0]) * mass1320 * TMath::Sqrt(width1320) / den1320;
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525) / den1525;
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710) / den1710;

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270) / den1270;
    double imagnum1320 = mass1320 * mass1320 * width1320 * TMath::Sqrt(width1320) / den1320;
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525) / den1525;
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710) / den1710;

    double fit1270 = yield1270 * (realnum1270 * realnum1270 + imagnum1270 * imagnum1270);
    double fit1320 = yield1320 * (realnum1320 * realnum1320 + imagnum1320 * imagnum1320);
    double fit1525 = yield1525 * (realnum1525 * realnum1525 + imagnum1525 * imagnum1525);
    double fit1710 = yield1710 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710);
    double fit = fit1270 + fit1320 + fit1525 + fit1710;

    return fit;
}

Double_t coherent_sum(double *x, double *par) // taken 4 resonances here
{
    // Fit is |a1 *BW1 + a2*BW2 e^{i*phi1} + a3*BW3 e^{i*phi2}|^2 + |a4 *BW4|^2
    // Then the real and imaginary parts are separated
    // Total 14 parameters. 3 for each resonance and 2 for phases.

    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
    double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    // double temp = par[3];
    // double temp2 = par[6];

    double norm1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2] * (TMath::Power(par[1] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double norm1320 = par[3];
    // double norm1320 = par[0] * 3 / 5;
    double mass1320 = par[4];
    double width1320 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double norm1525 = par[6];
    // double norm1525 = par[0] * 2 / 5;
    double mass1525 = par[7];
    double width1525 = par[8] * (TMath::Power(par[7] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n1);
    double norm1710 = par[9];
    double mass1710 = par[10];
    double width1710 = par[11] * (TMath::Power(par[10] / x[0], 1.0)) * TMath::Power((npart1) / (dpart4), n2);

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1320 = (x[0] * x[0] - mass1320 * mass1320) * (x[0] * x[0] - mass1320 * mass1320) + mass1320 * mass1320 * width1320 * width1320;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270) / den1270;
    double realnum1320 = (mass1320 * mass1320 - x[0] * x[0]) * mass1320 * TMath::Sqrt(width1320) / den1320;
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525) / den1525;
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710) / den1710;

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270) / den1270;
    double imagnum1320 = mass1320 * mass1320 * width1320 * TMath::Sqrt(width1320) / den1320;
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525) / den1525;
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710) / den1710;

    double phase1 = par[12]; // this is angle phi in radians
    double phase2 = par[13];

    double real1 = realnum1270;
    double imag1 = imagnum1270;
    double real2 = realnum1320 * TMath::Cos(phase1) - imagnum1320 * TMath::Sin(phase1);
    double imag2 = realnum1320 * TMath::Sin(phase1) + imagnum1320 * TMath::Cos(phase1);
    double real3 = realnum1525 * TMath::Cos(phase2) - imagnum1525 * TMath::Sin(phase2);
    double imag3 = realnum1525 * TMath::Sin(phase2) + imagnum1525 * TMath::Cos(phase2);
    double real4 = realnum1710;
    double imag4 = imagnum1710;

    double real_sum = norm1270 * real1 + norm1320 * real2 + norm1525 * real3;
    double imag_sum = norm1270 * imag1 + norm1320 * imag2 + norm1525 * imag3;
    double fit = norm1710 * norm1710 * (real4 * real4 + imag4 * imag4) + (real_sum * real_sum + imag_sum * imag_sum);
    // double fit = norm1270 * (real1 * real1 + imag1 * imag1) + norm1320 * (real2 * real2 + imag2 * imag2) + norm1525 * (real3 * real3 + imag3 * imag3) + norm1710 * (real4 * real4 + imag4 * imag4); // independent sum (for cross check)

    return fit;
}

Double_t single_BW(double *x, double *par)
{
    double yield = par[0];
    double mass = par[1];
    double width = par[2];

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));
    return fit;
}

Double_t BWsum_constWidth(double *x, double *par)
{
    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2];
    double yield1320 = par[3];
    double mass1320 = par[4];
    double width1320 = par[5];
    double yield1525 = par[6];
    double mass1525 = par[7];
    double width1525 = par[8];
    double yield1710 = par[9];
    double mass1710 = par[10];
    double width1710 = par[11];

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1320 = yield1320 * mass1320 * width1320 * x[0] / (pow((x[0] * x[0] - mass1320 * mass1320), 2) + pow(mass1320 * width1320, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    double fit = fit1270 + fit1320 + fit1525 + fit1710;
    return fit;
}

Double_t BWsumMassDepWidth(double *x, double *par)
{
    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
    double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    // double phase_space = (x[0] / TMath::Sqrt(x[0] * x[0] + 15 * 15)) * (TMath::Exp(-TMath::Sqrt(x[0] * x[0] + 15 * 15) / 0.16)); // 160 MeV is the kinetic freeze-out temperature

    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2] * (TMath::Power(par[1] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double yield1320 = par[3];
    double mass1320 = par[4];
    double width1320 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double yield1525 = par[6];
    double mass1525 = par[7];
    double width1525 = par[8] * (TMath::Power(par[7] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n1);
    double yield1710 = par[9];
    double mass1710 = par[10];
    double width1710 = par[11] * (TMath::Power(par[10] / x[0], 1.0)) * TMath::Power((npart1) / (dpart4), n2);

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1320 = yield1320 * mass1320 * width1320 * x[0] / (pow((x[0] * x[0] - mass1320 * mass1320), 2) + pow(mass1320 * width1320, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    double fit = (fit1270 + fit1320 + fit1525 + fit1710);
    return fit;
}

// Now we will define the functions for the exponential background

Double_t simple_exponential(double *x, double *par) // 2 parameters
{
    return (par[0] * pow(x[0], par[1]) * TMath::Exp(-par[2] * x[0]));
}

Double_t exponential_bkg_1(double *x, double *par) // 3 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * (x[0] - 2.0 * 0.497)));
}

Double_t exponential_bkg_2(double *x, double *par) // 3 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[1])));
}

Double_t exponential_bkg_3(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[3])));
}

Double_t exponential_bkg_4(double *x, double *par) // 5 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * ((pow((x[0] - 2.0 * 0.497), par[3])) + pow((x[0] - 2.0 * 0.497), par[4]))));
}

Double_t exponential_bkg_5(double *x, double *par) // 3 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * x[0]));
}

Double_t exponential_bkg_6(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] + par[1] * x[0] + par[3] * x[0] * x[0]));
}

// Now we will define the functions for the Boltzmann background
Double_t Boltzmann_bkg_1(double *x, double *par) // 3 parameters
{
    double norm = par[0];
    double massks = 0.497;
    double n = par[1];
    double C = par[2];
    double y = norm * pow((x[0] - 2 * massks), n / 2) * pow(C, 3 / 2) * exp(-C * pow((x[0] - 2 * massks), n));
    return y;
}

Double_t Boltzmann_bkg_2(double *x, double *par) // 4 parameters
{
    double norm = par[0];
    double massks = 0.497;
    double n = par[1];
    double C = par[2];
    double n1 = par[3];
    double y = norm * pow((x[0] - 2 * massks), n / 2) * pow(C, 3 / 2) * exp(-C * pow((x[0] - 2 * massks), n1));
    return y;
}

// Now we will define the functions for the sum of the Breit-Wigner and the exponential background

Double_t single_BW_expol3(double *x, double *par)
{
    return (single_BW(x, par) + exponential_bkg_3(x, &par[3]));
}

Double_t single_BW_expol3_hera(double *x, double *par)
{
    return (single_BW_hera(x, par) + exponential_bkg_3(x, &par[3]));
}

Double_t BWsum_ConstWidth_ModBolt(double *x, double *par)
{
    return (BWsum_constWidth(x, par) + exponential_bkg_3(x, &par[12]));
}

Double_t BWsum_expol3_hera(double *x, double *par)
{
    return (BWsum_hera(x, par) + exponential_bkg_3(x, &par[12]));
}

// Now we will define the functions for the sum of the Breit-Wigner and the Boltzmann background

Double_t single_BW_boltzman_1(double *x, double *par)
{
    return (single_BW(x, par) + Boltzmann_bkg_1(x, &par[3]));
}

Double_t single_BW_boltzman_2(double *x, double *par)
{
    return (single_BW(x, par) + Boltzmann_bkg_2(x, &par[3]));
}

Double_t BWsum_boltzman_1(double *x, double *par)
{
    return (BWsumMassDepWidth(x, par) + Boltzmann_bkg_1(x, &par[12]));
}

Double_t BWsum_boltzman_2(double *x, double *par)
{
    return (BWsum_constWidth(x, par) + Boltzmann_bkg_2(x, &par[12]));
}
Double_t expol_chkstar(double *x, double *par)
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] * x[0] + x[0] * x[0] * par[3]));
}
Double_t BWsum_expol_chkstar(double *x, double *par)
{
    return (BWsumMassDepWidth(x, par) + expol_chkstar(x, &par[12]));
}

Double_t BWsumMassDepWidth_exponential(double *x, double *par)
{
    // return (BWsumMassDepWidth(x, par) + expol_chkstar(x, &par[12]));
    return (BWsumMassDepWidth(x, par) + exponential_bkg_3(x, &par[12]));
}

Double_t BWsumMassDepWidth_simple_exponential(double *x, double *par)
{
    return (BWsumMassDepWidth(x, par) + simple_exponential(x, &par[12]));
}

Double_t BWsum_modifiedBoltzmann_hera(double *x, double *par)
{
    return (BWsum_hera(x, par) + exponential_bkg_3(x, &par[10]));
}
Double_t BWsum_ModifiedBoltzmann_hera_mass_dep(double *x, double *par)
{
    return (BWsum_hera_mass_dep(x, par) + exponential_bkg_3(x, &par[10]));
}
Double_t BWsum_modifiedBoltzmann_hera_const(double *x, double *par)
{
    return (BWsum_hera_const(x, par) + exponential_bkg_3(x, &par[12]));
}
Double_t CoherentSum_modifiedBoltzmann(double *x, double *par)
{
    // return (coherent_sum(x, par) + expol_chkstar(x, &par[14]));
    return (coherent_sum(x, par) + exponential_bkg_3(x, &par[14]));
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1); // top pad
    TPad *pad2 = (TPad *)c->GetPad(2); // bottom pad
    pad2Size = 0.5;                    // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.5, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.5);
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.23);
    pad1->SetLeftMargin(0.125);
    pad2->SetLeftMargin(0.125);
    pad1->SetTopMargin(0.1);
    pad1->SetBottomMargin(0);
    pad2->SetTopMargin(0);
}