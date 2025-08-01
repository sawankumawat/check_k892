#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <TArrow.h>
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);
Double_t single_BW_hera(double *x, double *par);
Double_t single_BW(double *x, double *par);
Double_t BWsum_hera(double *x, double *par);
Double_t BWsum(double *x, double *par);
Double_t BWsumMassDepWidth(double *x, double *par);
Double_t BWsumMassDepWidth_exponential(double *x, double *par);
Double_t single_BW_mass_dep_spin0(double *x, double *par);
Double_t single_BW_mass_dep_spin2(double *x, double *par);
Double_t BWsumMassDepWidth_simple_exponential(double *x, double *par);
Double_t BWsum_modifiedBoltzmann_hera(double *x, double *par);
Double_t BWsum_ModifiedBoltzmann_hera_mass_dep(double *x, double *par);
Double_t BWsum_modifiedBoltzmann_hera_const(double *x, double *par);
Double_t CoherentSum_modifiedBoltzmann(double *x, double *par);

Double_t exponential_bkg_1(double *x, double *par); // 3 parameters
Double_t exponential_bkg_2(double *x, double *par); // 3 parameters
Double_t exponential_bkg_3(double *x, double *par); // 4 parameters
Double_t exponential_bkg_4(double *x, double *par); // 5 parameters
Double_t exponential_bkg_5(double *x, double *par); // 3 parameters
Double_t exponential_bkg_6(double *x, double *par); // 4 parameters

Double_t Boltzmann_bkg_1(double *x, double *par); // 3 parameters
Double_t Boltzmann_bkg_2(double *x, double *par); // 4 parameters
Double_t expol_chkstar(double *x, double *par);   // 4 parameters
Double_t BWsum_expol_chkstar(double *x, double *par);
Double_t simple_exponential(double *x, double *par); // 2 parameters
Double_t BWsum_hera_const(double *x, double *par);
Double_t BWsum_hera_mass_dep(double *x, double *par); // taken 4 resonances here
Double_t coherent_sum(double *x, double *par);        // taken 4 resonances here

Double_t single_BW_expol3(double *x, double *par);
Double_t single_BW_expol3_hera(double *x, double *par);
Double_t BWsum_expol3(double *x, double *par);
Double_t BWsum_expol3_hera(double *x, double *par);

Double_t single_BW_boltzman_1(double *x, double *par);
Double_t single_BW_boltzman_2(double *x, double *par);
Double_t BWsum_boltzman_1(double *x, double *par);
Double_t BWsum_boltzman_2(double *x, double *par);

void glueball_fit_4rBW()
{
    // systematic studies (signal extraction) ****************************
    // A. fit range: Default: 1.05-2.20 GeV/c^2, Variation1: 1.02-2.20 GeV/c^2, Variation2: 1.05-2.30 GeV/c^2, Variation3: 1.08-2.15 GeV/c^2, Variation4: 1.02-2.30 GeV/c^2
    // B. Norm range: Default: 2.50-2.60 GeV/c^2, Variation1: 2.40-2.50 GeV/c^2, Variation2: 2.60-2.70 GeV/c^2
    // C. Fit function: Default: 4rBW with mass dependent width + modified Boltzmann, Variation1: 4rBW with constant width + bkg, Variation2: 3rBW with mass dependen width + bkg, Variation3: 4rBW with mass dependent width + Expol1, Variation4: 4rBW with mass dependent width + Boltzmann
    // D. Fit paramters: Default: Width of spin-2 resonances fixed to PDG, Variation1: Width of spin-2 resonances free, Variation2: Both mass and width of spin-2 resonances fixed to PDG, Variation3: Width of f1710 fixed to PDG, Variation4: Mass of f1710 fixed to PDG
    // E. Combinatorial background: Default: Rotational, Variation1: Mixed (Not considered)

    // systematic studies (Track selection) ****************************
    // TrA. DCA track to PV: Deafult: 0.05 cm, Variation1: 0.04 cm, Variation2: 0.06 cm
    // TrB. TPC PID: Default: 3sigma, Variation1: 4sigma, Variation2: 5sigma
    // TrC. TPC crossed rows: Default 70, Variation1: 100, Variation2: 120
    // TrD. TPC crossed rows over findable clusters: Default: 0.8, Variation1: 0.9, Variation2: 1.0

    // systematic studies (Topological selection) ****************************
    // ToA. Cosine PA: Default: 0.97, Variation1: 0.95, Variation2: 0.99
    // ToB. Transeverse radius: Default: 0.5 cm, Variation1: 0.4 cm, Variation2: 0.6 cm
    // ToC. DCA b/w V0 daughters: Default: 0.5 cm, Variation1: 0.3 cm, Variation2: 1.0 cm
    // ToD. Lifetime: Default: 20 cm, Variation1: 15 cm, Variation2: 25 cm
    // ToE. Competing V0 rejection: Default: 5 MeV, Variation1: 4, Variation2: 6
    // ToF. Ks mass window: Default: 3sigma, Variation1: 4sigma, Variation2: 5sigma

    string sysallvar[] = {"default", "varTrA1", "varTrA2", "varTrB1", "varTrB2", "varTrC1", "varTrC2", "varTrD1", "varTrD2", "varToA1", "varToA2", "varToB1", "varToB2", "varToC1", "varToC2", "varToD1", "varToD2", "varToE1", "varToE2", "varToF1", "varToF2"}; // only for track, pid and topological variations

    string kvariation1names[] = {"_id24937", "_DCA0p04_id24940", "_DCA0p06_id24940", "_TPCPID4_id24937", "_TPCPID6_id24937", "_TPCcr100_id24937", "_TPCcr120_id24937", "_TPCcrfc0p9_id24940", "_TPCcrfc1p0_id24940", "_cospa0p95_id24938", "_cospa0p99_id24938", "_decay_rad0p4_id24938", "_decay_rad0p6_id24938", "_DCAv0dau0p3_id24938", "_DCAv0dau1_id24938", "_lifetime15_id24939", "_lifetime25_id24939", "_lambda_rej4_id24939", "_lambda_rej6_id24939", "_Ks_selection4_id24939", "_Ks_selection5_id24939"}; // only for track, pid and topological variations

    int sizeofsysallvar = sizeof(sysallvar) / sizeof(sysallvar[0]);
    int sizeofkvariation1names = sizeof(kvariation1names) / sizeof(kvariation1names[0]);
    if (sizeofsysallvar != sizeofkvariation1names)
    {
        cout << "Size of variations and names are not same" << endl;
        return;
    }

    // for (int isysvars = 0; isysvars < sizeofsysallvar; isysvars++)
    {

        // string kvariation1 = kvariation1names[isysvars];
        // string sysvar = sysallvar[isysvars];

        //****************systematics train*******************************
        // Defaults
        const string kvariation1 = "_id24937"; // first four are same
        // const string kvariation1 = "_id24938";
        // const string kvariation1 = "_id24939";
        // const string kvariation1 = "_id24940";

        // Track selections
        //  const string kvariation1 = "_DCA0p04_id24940"; //DCA track to PV
        //  const string kvariation1 = "_DCA0p06_id24940"; //DCA track to PV
        //  const string kvariation1 = "_TPCPID3_id24937"; //TPC PID
        //  const string kvariation1 = "b_massdepWidth_modifiedBoltzmann_TPCPID4_id24937"; //TPC PID
        //  const string kvariation1 = "_TPCPID6_id24937"; //TPC PID
        //  const string kvariation1 = "_TPCcr100_id24937"; //TPC crossed rows
        //  const string kvariation1 = "_TPCcr120_id24937"; //TPC crossed rows
        //  const string kvariation1 = "_TPCcrfc0p9_id24940"; //TPC crossed rows over findable clusters
        // const string kvariation1 = "_TPCcrfc1p0_id24940"; // TPC crossed rows over findable clusters

        // Topological selections
        //  const string kvariation1 = "_cospa0p95_id24938"; //Cosine PA
        //  const string kvariation1 = "_cospa0p99_id24938";
        //  const string kvariation1 = "_decay_rad0p4_id24938"; //Transverse radius
        //  const string kvariation1 = "_decay_rad0p6_id24938";
        //  const string kvariation1 = "_DCAv0dau0p3_id24938"; //DCA b/w V0 daughters
        //  const string kvariation1 = "_DCAv0dau1_id24938";
        //  const string kvariation1 = "_lifetime15_id24939"; //Lifetime
        //  const string kvariation1 = "_lifetime25_id24939";
        //  const string kvariation1 = "_lambda_rej4_id24939";
        //  const string kvariation1 = "_lambda_rej6_id24939";
        // const string kvariation1 = "_Ks_selection4_id24939";
        //  const string kvariation1 = "_Ks_selection5_id24939";

        // *******************variations for Ks cuts and angular separation************************
        ////Defaults
        // const string kvariation1 = "_id24794";
        // const string kvariation1 = "_id25081";
        ////Variations
        // const string kvariation1 = "_1Kscut_id24794";
        // const string kvariation1 = "_1p5Kscut_id24794";
        // const string kvariation1 = "_2Kscut_id24794";
        // const string kvariation1 = "_4Kscut_id24794";
        // const string kvariation1 = "_angsep_0p5_id25081";
        // const string kvariation1 = "_angsep_1_id25081";
        // const string kvariation1 = "_angsep_1p5_id25081";
        // const string kvariation1 = "_angsep_2_id25081";
        // const string kvariation1 = "_angsep_3_id25081";

        //*********for systematics and default study with full train ************************
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/351470/KsKs_Channel/higher-mass-resonances_PID3";
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial";
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances" + kvariation1;
        // string path2 = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937";
        string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances";
        string path2 = path;

        // //*********for temporary study with angular separation cuts************************
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/362701/KsKs_Channel/higher-mass-resonances" + kvariation1;
        // string path2 = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/362701/KsKs_Channel/higher-mass-resonances_id24794";

        string sysvar = "default"; // default

        ofstream file;
        // file.open((path2 + "/fits/4rBw_fits/fit_params_" + sysvar + ".csv").c_str());
        file.open((path2 + "/fits/4rBw_fits/fit_params_temp_" + sysvar + ".txt").c_str());

        string savepath = path2 + "/fits/4rBw_fits";

        gSystem->Exec(("mkdir -p " + savepath).c_str());

        // TFile *f = new TFile((path + "/hglue_ROTATED_norm_2.50_2.60_pt_3.00_30.00.root").c_str(), "READ"); // default
        // TFile *f = new TFile((path + "/hglue_ROTATED_norm_2.50_2.60_all_pT.root").c_str(), "READ"); //
        // TFile *f = new TFile((path + "/hglue_MIX_temp.root").c_str(), "READ");
        TFile *f = new TFile((path + "/hglue_ROTATED.root").c_str(), "READ");

        TFile *plots_4BW = new TFile("root_files/4rBW_plots_expol.root", "RECREATE");

        int colors[] = {kGreen + 4, 28, kMagenta, kBlue};
        double masses[] = {f1270Mass, a1320Mass, f1525Mass, f1710Mass};
        double widths[] = {f1270Width, a1320Width, f1525Width, f1710Width};
        string resonance_names[] = {"f_{2}(1270)", "a_{2}(1320)^{0}", "f'_{2}(1525)", "f_{0}(1710)"};
        double purity, significance, chi2ndf, statSignificance, signal_counts, background_counts;

        if (f->IsZombie())
        {
            cout << "Error opening file" << endl;
            return;
        }

// #define b_constantWidth_modified_Boltzmann // without mass dependent width
// #define b_massdepWidth_Standard_boltzman
// #define b_massdepWidth_modifiedBoltzmann
#define b_massdepWidth_expol2
        // #define b_massdepWidth_HERAexponential
        // #define b_modifiedBoltzmann_hera_const // for real + img part without interference
        // #define b_modifiedBoltzmann_hera // for real + img part with interference
        // #define b_modifiedBoltzmann_hera_mass_dep // for real + img part with interference and mass dependent width
        // #define b_coherentSum_modifiedBoltzmann // for coherent sum with phases

        // #define residual_subtracted
        // #define doublepanelplot

        TH1F *hmult = (TH1F *)f->Get("multiplicity_histogram");
        if (hmult == nullptr)
        {
            cout << "Multiplicity histogram not found" << endl;
            return;
        }
        int multlow = 0;
        int multhigh = 100;
        double total_events = hmult->Integral(hmult->GetXaxis()->FindBin(multlow), hmult->GetXaxis()->FindBin(multhigh));
        double phi_mod, phi_mod2;

        for (int ipt = 0; ipt < Npt; ipt++)
        {
            // TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", 3.0, 5.0));
            TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", pT_bins[ipt], pT_bins[ipt + 1]));
            TH1F *hraw = (TH1F *)f->Get(Form("ksks_invmass_pt_%.1f_%.1f", pT_bins[ipt], pT_bins[ipt + 1]));
            if (hinvMass == nullptr)
            {
                cout << "Error opening histogram" << endl;
                return;
            }

            // TCanvas *cverytemp = new TCanvas("", "", 720, 720);
            // SetCanvasStyle(cverytemp, 0.14, 0.03, 0.05, 0.14);
            // hraw->Draw("pe");

            if (hinvMass == nullptr)
            {
                cout << "Error opening histogram" << endl;
                return;
            }
            TCanvas *c = new TCanvas("", "", 720, 720);
            SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
            hinvMass->Rebin(2);
            double binwidthfile = (hinvMass->GetXaxis()->GetXmax() - hinvMass->GetXaxis()->GetXmin()) / hinvMass->GetXaxis()->GetNbins();
            // cout<<"x max is "<<hinvMass->GetXaxis()->GetXmax()<<endl;
            // cout<<"x min is "<<hinvMass->GetXaxis()->GetXmin()<<endl;
            // cout<<"n bins is "<<hinvMass->GetXaxis()->GetNbins()<<endl;
            cout << "bin width is " << binwidthfile << endl;
            hinvMass->GetXaxis()->SetRangeUser(1.00, 2.50);
            hinvMass->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
            hinvMass->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidthfile * 1000));
            hinvMass->SetMaximum(1.2 * hinvMass->GetMaximum());
            hinvMass->GetYaxis()->SetTitleOffset(1.35);
            hinvMass->SetMarkerSize(1.0);
            hinvMass->Draw("pe");
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

            TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_exponential, 1.05, 2.20, 16); // expol 3

            string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
            for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
            {
                BEexpol->SetParName(i, parnames[i].c_str());
            }

            // double parameters[] = {1500, f1270Mass, f1270Width, 1000, a1320Mass, a1320Width, 3700, f1525Mass, f1525Width, 1300, f1710Mass, f1710Width}; // default
            // double parameters[] = {2.31e4, f1270Mass, f1270Width, 1.36e4, a1320Mass, a1320Width, 4.73e4, f1525Mass, f1525Width, 1.45e4, f1710Mass, f1710Width}; // LHC23_pass4_thin
            double parameters[] = {3500, f1270Mass, f1270Width, 2000, a1320Mass, a1320Width, 7000, f1525Mass, f1525Width, 2200, f1710Mass, f1710Width}; // rebin twice
            // double parameters[] = {1400, f1270Mass, f1270Width, 1000, a1320Mass, a1320Width, 2500, f1525Mass, f1525Width, 800, f1710Mass, f1710Width}; // medium train
            // double parameters[] = {400, f1270Mass, f1270Width, 370, a1320Mass, a1320Width, 1200, f1525Mass, f1525Width, 450, f1710Mass, f1710Width}; // pt range 3-5 GeV/c, 5-8 GeV/c
            // double parameters[] = {700, f1270Mass, f1270Width, 706, a1320Mass, a1320Width, 2200, f1525Mass, f1525Width, 1000, f1710Mass, f1710Width}; // pt range 3-5 GeV/c, 5-8 GeV/c rebin 2
            // double parameters[] = {35, f1270Mass, f1270Width, 66, a1320Mass, a1320Width, 166, f1525Mass, f1525Width, 280, f1710Mass, f1710Width}; // pt range 8-15 GeV/c rebin 3
            // double parameters[] = {690, f1270Mass, f1270Width, 714, a1320Mass, a1320Width, 2300, f1525Mass, f1525Width, 500, f1710Mass, f1710Width}; // pt range 2-3 GeV/c rebin 2 (fit 1.05-2.25)
            // double parameters[] = {1470, f1270Mass, f1270Width, 714, a1320Mass, a1320Width, 1300, f1525Mass, f1525Width, 250, f1710Mass, f1710Width}; // pt range 1-2 GeV/c rebin 2 (fit 1.05-2.20)

            // double parameters[] = {2000, f1270Mass, f1270Width, 1714, a1320Mass, a1320Width, 5500, f1525Mass, f1525Width, 3000, f1710Mass, f1710Width}; // ME 3-30 GeV/c

            int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

            for (int i = 0; i < size_fitparams; i++)
            {
                BEexpol->SetParameter(i, parameters[i]);
            }
            vector<vector<double>> par_limits = {{1, 2 * f1270Width}, {4, 2 * a1320Width}, {7, 2 * f1525Width}, {10, 2 * f1710Width}, {11, 5 * f1710WidthErr}};

            int limits_size = par_limits.size();
            for (int i = 0; i < limits_size; i++)
            {
                int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            }

            // //********systematic studies*************
            // double initial_param_bkg[] = {2.37518e5, -0.102044, 3.071167, 1.354864}; // 0-30 GeV/c, 3sigma
            double initial_param_bkg[] = {7.37518e5, 0.0134, 3.071167, 1.04}; // rebin twice
            // double initial_param_bkg[] = {4.618e6, 0.00774, 2.82, 1.03}; // pass_4_thin
            // double initial_param_bkg[] = {2.57518e5, 0.0424, 3.171167, 0.964}; // 0-30 GeV/c, 3sigma
            // double initial_param_bkg[] = {4.8e4, -0.107, 3.51167, 1.4}; //pt range 3-5, 5-8 GeV/c
            // double initial_param_bkg[] = {7.0e4, -0.17, 3.867, 1.8}; // pt range 3-5, 5-8 GeV/c (rebin 2)
            // double initial_param_bkg[] = {3.7e5, 0.07, 3.7, 0.94}; //rot range 1 - 30
            // double initial_param_bkg[] = {1.1e5, -0.17, 4.00, 1.4}; //rot range 2 - 30
            // double initial_param_bkg[] = {7646, -0.06, 2.15, 1.4}; // pt range 8-15 GeV/c (rebin 3)
            // double initial_param_bkg[] = {8.2e4, -0.257, 5.715, 1.8}; // pt range 2-3 GeV/c (rebin 2)
            // double initial_param_bkg[] = {7.3e6, 0.57, 6.15, 0.4}; // pt range 1-2 GeV/c (rebin 2)
            // double initial_param_bkg[] = {2.3e5, -0.12, 3.15, 1.4}; // ME 3-30 GeV/c

            // Initial parameters for background
            BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Free
            BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Fix for medium train
            BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Free
            BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

            BEexpol->FixParameter(2, f1270Width);
            BEexpol->FixParameter(5, a1320Width);
            BEexpol->FixParameter(8, f1525Width);

            // BEexpol->FixParameter(1, f1270Mass);
            // BEexpol->FixParameter(4, a1320Mass);
            // BEexpol->FixParameter(7, f1525Mass);

            // BEexpol->FixParameter(10, f1710Mass);
            // BEexpol->FixParameter(11, f1710Width);

            TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
            chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
            cout << "chi2/ndf is " << chi2ndf << endl;
            string fitstatus = "Successfull";
            // cout<<"fit status code "<<fitResultptr->Status()<<endl;
            // if (fitResultptr->Status() != 4140)
            // {
            //     cout << "Fit failed or call limit !!!!!!!" << endl;
            //     fitstatus = "Failed";
            // }
            // fitResultptr->Print("V");

            double *obtained_parameters = BEexpol->GetParameters();
            TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
            TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
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

            TF1 *onlyBW = new TF1("onlyBW", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            string parameter_names[] = {"norm1270", "mass1270", "width1270", "norm1320", "mass1320", "width1320", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710"};
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

            onlyBW_clone->FixParameter(2, f1270Width);
            onlyBW_clone->FixParameter(5, a1320Width);
            onlyBW_clone->FixParameter(8, f1525Width);

            // onlyBW_clone->FixParameter(1, f1270Mass);
            // onlyBW_clone->FixParameter(4, a1320Mass);
            // onlyBW_clone->FixParameter(7, f1525Mass);

            // onlyBW_clone->FixParameter(11, f1710Width);
            // onlyBW_clone->FixParameter(10, f1710Mass);

            onlyBW->SetLineColor(4);
            onlyBW->SetLineStyle(2);
            // onlyBW->Draw("same");

            // // Now plot the indivial resonances
            TF1 *singlefits[4];
            for (int i = 0; i < 4; i++)
            {
                singlefits[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                singlefits[i]->SetLineColor(colors[i]);
                singlefits[i]->SetLineStyle(2);
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
            ltemp->Draw("same");

            TLatex lat1;
            lat1.SetNDC();
            lat1.SetTextSize(0.03);
            lat1.SetTextFont(42);
            lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{#it{s}} = 13.6 TeV");
            lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
            lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/#it{c}", pT_bins[ipt], pT_bins[ipt + 1]));

            for (int i = 0; i < 4; i++)
            {
                double significance_num = singlefits[i]->Integral(masses[i] - 2 * widths[i], masses[i] + 2 * widths[i]) / binwidthfile;
                int binlow = hraw->GetXaxis()->FindBin(masses[i] - 2 * widths[i]);
                int binhigh = hraw->GetXaxis()->FindBin(masses[i] + 2 * widths[i]);
                double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
                significance = significance_num / significance_den;
                signal_counts = significance_num;
                background_counts = hraw->Integral(binlow, binhigh) - signal_counts;

                cout << "numerator " << significance_num << " denominator " << significance_den << endl;
                cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                double amplitude = obtained_parameters[3 * i];
                double amplitude_err = BEexpol->GetParError(3 * i);
                statSignificance = amplitude / amplitude_err;
                cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
            }
#endif

            // // // //************************************************************************ */
            // // // // **************** For BW sum with exp + pol2 as used in Charged kstar **************************
#ifdef b_massdepWidth_expol2
            TF1 *BEexpol = new TF1("BEexpol", BWsum_expol_chkstar, 1.05, 2.20, 16); // expol 2
            string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
            for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
            {
                BEexpol->SetParName(i, parnames[i].c_str());
            }

            // double parameters[] = {2200, f1270Mass, f1270Width, 1800, a1320Mass, a1320Width, 3700, f1525Mass, f1525Width, 1500, f1710Mass, f1710Width}; // KsKs channel default
            // double parameters[] = {3500, f1270Mass, f1270Width, 2000, a1320Mass, a1320Width, 7500, f1525Mass, f1525Width, 2300, f1710Mass, f1710Width}; // KsKs channel default rebin
            // double parameters[] = {1000, f1270Mass, f1270Width, 1160, a1320Mass, a1320Width, 2500, f1525Mass, f1525Width, 850, f1710Mass, f1710Width}; // 3-5 GeV/c rebin twice

            double parameters[] = {3200, f1270Mass, f1270Width, 3000, a1320Mass, a1320Width, 7500, f1525Mass, f1525Width, 2200, f1710Mass, f1710Width}; // same bkg as f0 and f1

            int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

            for (int i = 0; i < size_fitparams; i++)
            {
                BEexpol->SetParameter(i, parameters[i]);
            }

            vector<vector<float>> par_limits = {{1, 2 * f1270Width}, {4, 2 * a1320Width}, {7, 2 * f1525Width}, {10, 2 * f1710Width}, {11, 10 * f1710WidthErr}};

            int limits_size = par_limits.size();
            for (int i = 0; i < limits_size; i++)
            {
                int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            }

            // double initial_param_bkg[] = {370.9, 0.0741322, 11.6184, -4.395450}; //default
            // double initial_param_bkg[] = {150.9, 0.039, 11.6184, -3.0450}; // default rebin
            // double initial_param_bkg[] = {250.9, 0.22, 11.24, -4.395450}; // 3-5 GeV/c rebin twice
            double initial_param_bkg[] = {9.5e6, -0.007, -2.4, -0.15}; // exact form as f0 and f1

            // Initial parameters for background
            BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Fix
            BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Free
            BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Fix
            BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

            BEexpol->FixParameter(2, f1270Width);
            BEexpol->FixParameter(5, a1320Width);
            BEexpol->FixParameter(8, f1525Width);

            // BEexpol->FixParameter(1, f1270Mass);
            // BEexpol->FixParameter(4, a1320Mass);
            // BEexpol->FixParameter(7, f1525Mass);

            // BEexpol->FixParameter(10, f1710Mass);
            // BEexpol->FixParameter(11, f1710Width);

            hinvMass->Fit("BEexpol", "REBMS");
            TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");

            string fitstatus = "Successfull";
            // if (fitResultptr->Status() != 4000)
            // {
            //     cout << "Fit failed or call limit" << endl;
            //     fitstatus = "Failed";
            // }
            // cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
            // fitResultptr->Print("V");
            chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();

            double *obtained_parameters = BEexpol->GetParameters();
            TF1 *expol = new TF1("expol", expol_chkstar, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
            TF1 *expol_clone = new TF1("expol_clone", expol_chkstar, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
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

            TF1 *onlyBW = new TF1("onlyBW", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            for (int i = 0; i < 12; i++)
            {
                onlyBW->SetParameter(i, obtained_parameters[i]);
                onlyBW_clone->SetParameter(i, obtained_parameters[i]);
                onlyBW_clone->SetParName(i, parnames[i].c_str());
            }

            onlyBW_clone->FixParameter(2, f1270Width);
            onlyBW_clone->FixParameter(5, a1320Width);
            onlyBW_clone->FixParameter(8, f1525Width);

            // onlyBW_clone->FixParameter(1, f1270Mass);
            // onlyBW_clone->FixParameter(4, a1320Mass);
            // onlyBW_clone->FixParameter(7, f1525Mass);

            // onlyBW_clone->FixParameter(11, f1710Width);
            // onlyBW_clone->FixParameter(10, f1710Mass);

            onlyBW->SetLineColor(4);
            onlyBW->SetLineStyle(2);
            // onlyBW->Draw("same");

            // // Now plot the indivial resonances
            TF1 *singlefits[4];
            for (int i = 0; i < 4; i++)
            {
                singlefits[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                singlefits[i]->SetLineColor(colors[i]);
                singlefits[i]->SetLineStyle(2);
                singlefits[i]->Draw("same");
            }

            TLegend *ltemp = new TLegend(0.25, 0.53, 0.55, 0.87);
            ltemp->SetFillStyle(0);
            ltemp->SetBorderSize(0);
            ltemp->SetTextFont(42);
            ltemp->SetTextSize(0.03);
            ltemp->AddEntry((TObject *)0, "", "");
            ltemp->AddEntry((TObject *)0, "", "");
            ltemp->AddEntry(hinvMass, "Data (stat. uncert.)", "lpe");
            ltemp->AddEntry(BEexpol, "4rBW + Residual BG", "l");
            // ltemp->AddEntry(onlyBW, "4rBW", "l");
            ltemp->AddEntry(expol, "Residual BG", "l");
            ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
            ltemp->AddEntry(singlefits[1], "a_{2}(1320)^{0}", "l");
            ltemp->AddEntry(singlefits[2], "f'_{2}(1525)", "l");
            ltemp->AddEntry(singlefits[3], "f_{0}(1710)", "l");
            ltemp->Draw("same");

            TLatex lat1;
            lat1.SetNDC();
            lat1.SetTextSize(0.03);
            lat1.SetTextFont(42);
            lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
            lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
            lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

            for (int i = 0; i < 4; i++)
            {
                double significance_num = singlefits[i]->Integral(masses[i] - 3 * widths[i], masses[i] + 3 * widths[i]) / binwidthfile;
                int binlow = hraw->GetXaxis()->FindBin(masses[i] - 3 * widths[i]);
                int binhigh = hraw->GetXaxis()->FindBin(masses[i] + 3 * widths[i]);
                double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
                significance = significance_num / significance_den;
                cout << "numerator " << significance_num << " denominator " << significance_den << endl;
                cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                double amplitude = obtained_parameters[3 * i];
                double amplitude_err = BEexpol->GetParError(3 * i);
                statSignificance = amplitude / amplitude_err;
                cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
            }

#endif

            // **********simple exponential HERA*********************

#ifdef b_massdepWidth_HERAexponential

            TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_simple_exponential, 1.05, 2.20, 15);

            string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b"};
            for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
            {
                BEexpol->SetParName(i, parnames[i].c_str());
            }

            double parameters[] = {3719, f1270Mass, f1270Width, 2000, a1320Mass, a1320Width, 7300, f1525Mass, f1525Width, 2227, f1710Mass, f1710Width}; // default

            int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

            for (int i = 0; i < size_fitparams; i++)
            {
                BEexpol->SetParameter(i, parameters[i]);
            }
            vector<vector<float>> par_limits = {{1, 1 * f1270Width}, {4, 1 * a1320Width}, {7, 1 * f1525Width}, {10, 1 * f1710Width}};
            int limits_size = par_limits.size();
            for (int i = 0; i < limits_size; i++)
            {
                int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            }

            double initial_param_bkg[] = {2.2e7, 0.7, 3.4}; // 0-30 GeV/c, 3sigma

            // Initial parameters for background
            BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]);
            BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]);
            BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]);

            BEexpol->FixParameter(2, f1270Width);
            BEexpol->FixParameter(5, a1320Width);
            BEexpol->FixParameter(8, f1525Width);

            // BEexpol->FixParameter(1, f1270Mass);
            // BEexpol->FixParameter(4, a1320Mass);
            // BEexpol->FixParameter(7, f1525Mass);

            // BEexpol->FixParameter(10, f1710Mass);
            // BEexpol->FixParameter(11, f1710Width);

            TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
            chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
            cout << "chi2/ndf is " << chi2ndf << endl;

            double *obtained_parameters = BEexpol->GetParameters();
            TF1 *expol = new TF1("expol", simple_exponential, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);             //
            TF1 *expol_clone = new TF1("expol_clone", simple_exponential, BEexpol->GetXmin(), BEexpol->GetXmax(), 3); //
            for (int i = 0; i < 3; i++)
            {
                expol->SetParameter(i, obtained_parameters[size_fitparams + i]);
                expol_clone->SetParameter(i, obtained_parameters[size_fitparams + i]);
            }
            expol->SetLineColor(3);
            expol->SetLineStyle(2);
            expol_clone->SetLineColor(3);
            expol_clone->SetLineStyle(2);
            expol->Draw("same");

            TF1 *onlyBW = new TF1("onlyBW", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            string parameter_names[] = {"norm1270", "mass1270", "width1270", "norm1320", "mass1320", "width1320", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710"};
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

            onlyBW_clone->FixParameter(2, f1270Width);
            onlyBW_clone->FixParameter(5, a1320Width);
            onlyBW_clone->FixParameter(8, f1525Width);

            // onlyBW_clone->FixParameter(1, f1270Mass);
            // onlyBW_clone->FixParameter(4, a1320Mass);
            // onlyBW_clone->FixParameter(7, f1525Mass);

            // onlyBW_clone->FixParameter(11, f1710Width);
            // onlyBW_clone->FixParameter(10, f1710Mass);

            onlyBW->SetLineColor(4);
            onlyBW->SetLineStyle(2);
            // onlyBW->Draw("same");

            // // Now plot the indivial resonances
            TF1 *singlefits[4];
            for (int i = 0; i < 4; i++)
            {
                singlefits[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                singlefits[i]->SetLineColor(colors[i]);
                singlefits[i]->SetLineStyle(2);
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
            ltemp->Draw("same");

            TLatex lat1;
            lat1.SetNDC();
            lat1.SetTextSize(0.03);
            lat1.SetTextFont(42);
            lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
            lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
            lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

            for (int i = 0; i < 4; i++)
            {
                double significance_num = singlefits[i]->Integral(masses[i] - 3 * widths[i], masses[i] + 3 * widths[i]) / binwidthfile;
                int binlow = hinvMass->GetXaxis()->FindBin(masses[i] - 3 * widths[i]);
                int binhigh = hinvMass->GetXaxis()->FindBin(masses[i] + 3 * widths[i]);
                double significance_den = TMath::Sqrt(hinvMass->Integral(binlow, binhigh));
                significance = significance_num / significance_den;
                cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                purity = significance_num * 100 / hinvMass->Integral(binlow, binhigh);
                cout << "Purity of " << resonance_names[i] << " is " << purity << endl;
            }
#endif

// // // Default fitting range is 1.02 to 2.20. Four types of fitting range variations: extend left (1.0), extend right (2.50), large range (1.0 to 2.50), small range (1.05 to 2.15)
#ifdef b_constantWidth_modified_Boltzmann

            TF1 *BEexpol = new TF1("BEexpol", BWsum_expol3, 1.05, 2.20, 16); // expol 3
            string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
            for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
            {
                BEexpol->SetParName(i, parnames[i].c_str());
            }

            // double parameters[] = {2263, f1270Mass, f1270Width, 1580, a1320Mass, a1320Width, 3800, f1525Mass, f1525Width, 1500, f1710Mass, f1710Width}; // default
            double parameters[] = {4500, f1270Mass, f1270Width, 3200, a1320Mass, a1320Width, 7500, f1525Mass, f1525Width, 3000, f1710Mass, f1710Width}; // default rebin 2
            // double parameters[] = {1400, f1270Mass, f1270Width, 1000, a1320Mass, a1320Width, 2200, f1525Mass, f1525Width, 750, f1710Mass, f1710Width}; // pt 3-5 GeV/c rebin twice

            int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

            for (int i = 0; i < size_fitparams; i++)
            {
                BEexpol->SetParameter(i, parameters[i]);
            }
            vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {4, 3 * a1320Width}, {7, 3 * f1525Width}, {10, 5 * f1710Width}, {11, 5 * f1710WidthErr}};
            int limits_size = par_limits.size();
            for (int i = 0; i < limits_size; i++)
            {
                int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            }

            // double initial_param_bkg[] = {3.073e5, -0.04378, 2.727, 1.011}; // rotational 0-30 GeV/c (KsKs channel)
            double initial_param_bkg[] = {5.43e5, -0.08, 2.627, 1.06}; // rotational 0-30 GeV/c rebin
            // double initial_param_bkg[] = {5.1e4, -0.593, 1.6, 3.781}; // 3-5 GeV/c rebin twice

            // Initial parameters for background
            BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5
            BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379
            BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569
            BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098

            BEexpol->FixParameter(2, f1270Width);
            BEexpol->FixParameter(5, a1320Width);
            BEexpol->FixParameter(8, f1525Width);

            // BEexpol->FixParameter(1, f1270Mass);
            // BEexpol->FixParameter(4, a1320Mass);
            // BEexpol->FixParameter(7, f1525Mass);

            // BEexpol->FixParameter(10, f1710Mass);
            // BEexpol->FixParameter(11, f1710Width);

            TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
            cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
            chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();

            double *obtained_parameters = BEexpol->GetParameters();
            TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
            TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
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

            TF1 *onlyBW = new TF1("onlyBW", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            string parameter_names[] = {"norm1270", "mass1270", "width1270", "norm1320", "mass1320", "width1320", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710"};
            for (int i = 0; i < 12; i++)
            {
                onlyBW->SetParameter(i, obtained_parameters[i]);
                onlyBW_clone->SetParameter(i, obtained_parameters[i]);
                onlyBW_clone->SetParName(i, parameter_names[i].c_str());
            }
            onlyBW->SetLineColor(4);
            onlyBW->SetLineStyle(2);
            // onlyBW->Draw("same");

            onlyBW_clone->FixParameter(2, f1270Width);
            onlyBW_clone->FixParameter(5, a1320Width);
            onlyBW_clone->FixParameter(8, f1525Width);

            // onlyBW_clone->FixParameter(11, f1710Width);
            // onlyBW_clone->FixParameter(1, f1270Mass);
            // onlyBW_clone->FixParameter(4, a1320Mass);
            // onlyBW_clone->FixParameter(7, f1525Mass);
            // onlyBW_clone->FixParameter(10, f1710Mass);

            // TLegend *ltemp = new TLegend(0.20, 0.67, 0.52, 0.92);
            // ltemp->SetFillStyle(0);
            // ltemp->SetTextFont(42);
            // ltemp->SetTextSize(0.035);
            // ltemp->AddEntry(hinvMass, "Data", "lpe");
            // ltemp->AddEntry(BEexpol, "4rBw + expol", "l");
            // ltemp->AddEntry(onlyBW, "4rBw", "l");
            // ltemp->AddEntry(expol, "expol", "l");
            // ltemp->Draw("same");

            // // Now plot the indivial resonances
            TF1 *singlefits[4];
            for (int i = 0; i < 4; i++)
            {
                singlefits[i] = new TF1(Form("singlef%d", i), single_BW, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                singlefits[i]->SetLineColor(colors[i]);
                singlefits[i]->SetLineStyle(2);
                singlefits[i]->Draw("same");
            }

            TLegend *ltemp = new TLegend(0.25, 0.55, 0.55, 0.87);
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
            ltemp->AddEntry(singlefits[2], "a_{2}^{0}(1320)", "l");
            ltemp->AddEntry(singlefits[2], "f'_{2}(1525)", "l");
            ltemp->AddEntry(singlefits[3], "f_{0}(1710)", "l");
            ltemp->Draw("same");

            TLatex lat1;
            lat1.SetNDC();
            lat1.SetTextSize(0.03);
            lat1.SetTextFont(42);
            lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
            lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
            // lat1.DrawLatex(0.255, )
            lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

            for (int i = 0; i < 4; i++)
            {
                double significance_num = singlefits[i]->Integral(masses[i] - 3 * widths[i], masses[i] + 3 * widths[i]) / binwidthfile;
                int binlow = hraw->GetXaxis()->FindBin(masses[i] - 3 * widths[i]);
                int binhigh = hraw->GetXaxis()->FindBin(masses[i] + 3 * widths[i]);
                double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
                significance = significance_num / significance_den;
                cout << "numerator " << significance_num << " denominator " << significance_den << endl;
                cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                double amplitude = obtained_parameters[3 * i];
                double amplitude_err = BEexpol->GetParError(3 * i);
                statSignificance = amplitude / amplitude_err;
                cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
            }
#endif

            // // // ************************************************************************************
            // // // **************** For BW sum with Boltzmann ****************************
#ifdef b_massdepWidth_Standard_boltzman

            TF1 *BEexpol = new TF1("BEexpol", BWsum_boltzman_1, 1.05, 2.20, 15); // expol 3

            string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "A", "n", "b", "c"};
            for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
            {
                BEexpol->SetParName(i, parnames[i].c_str());
            }

            // double parameters[] = {2000, f1270Mass, f1270Width, 1500, a1320Mass, a1320Width, 4000, f1525Mass, f1525Width, 1500, f1710Mass, f1710Width}; //default
            double parameters[] = {4800, f1270Mass, f1270Width, 2500, a1320Mass, a1320Width, 8000, f1525Mass, f1525Width, 3000, f1710Mass, f1710Width}; // default rebin
            // double parameters[] = {800, f1270Mass, f1270Width, 1300, a1320Mass, a1320Width, 2500, f1525Mass, f1525Width, 860, f1710Mass, f1710Width}; // 3-5 GeV/c rebin twice

            int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

            for (int i = 0; i < size_fitparams; i++)
            {
                BEexpol->SetParameter(i, parameters[i]);
            }
            // BEexpol->SetParLimits(9, 1400, 1600); // this is to constrain the width of f1710 (default)

            vector<vector<float>> par_limits = {{1, 2 * f1270Width}, {4, 2 * a1320Width}, {7, 2 * f1525Width}, {10, 5 * f1710Width}, {11, 10 * f1710WidthErr}};

            int limits_size = par_limits.size();
            for (int i = 0; i < limits_size; i++)
            {
                int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            }

            // double initial_param_bkg[] = {3.274e5, 0.7029, 4.391}; //default
            double initial_param_bkg[] = {6.44e5, 0.729, 4.391}; // default rebin
            // double initial_param_bkg[] = {1.274e5, 0.87029, 4.391};

            // Initial parameters for background
            BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5
            BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379
            BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569

            BEexpol->FixParameter(2, f1270Width);
            BEexpol->FixParameter(5, a1320Width);
            BEexpol->FixParameter(8, f1525Width);

            // BEexpol->FixParameter(1, f1270Mass);
            // BEexpol->FixParameter(4, a1320Mass);
            // BEexpol->FixParameter(7, f1525Mass);

            // BEexpol->FixParameter(10, f1710Mass);
            // BEexpol->FixParameter(11, f1710Width);

            TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
            cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
            chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
            string fitstatus = "Successfull";
            // cout<<"fit status code "<<fitResultptr->Status()<<endl;
            if (fitResultptr->Status() != 4140)
            {
                cout << "Fit failed or call limit !!!!!!!" << endl;
                fitstatus = "Failed";
            }

            double *obtained_parameters = BEexpol->GetParameters();
            TF1 *expol = new TF1("expol", Boltzmann_bkg_1, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
            TF1 *expol_clone = new TF1("expol_clone", Boltzmann_bkg_1, BEexpol->GetXmin(), BEexpol->GetXmax(), 3); //
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

            TF1 *onlyBW = new TF1("onlyBW", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            string parameter_names[] = {"norm1270", "mass1270", "width1270", "norm1320", "mass1320", "width1320", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710"};
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

            onlyBW_clone->FixParameter(2, f1270Width);
            onlyBW_clone->FixParameter(5, a1320Width);
            onlyBW_clone->FixParameter(8, f1525Width);

            onlyBW_clone->FixParameter(1, f1270Mass);
            onlyBW_clone->FixParameter(4, a1320Mass);
            onlyBW_clone->FixParameter(7, f1525Mass);

            onlyBW_clone->FixParameter(11, f1710Width);
            onlyBW_clone->FixParameter(10, f1710Mass);

            onlyBW->SetLineColor(4);
            onlyBW->SetLineStyle(2);
            // onlyBW->Draw("same");

            // // Now plot the indivial resonances
            TF1 *singlefits[4];
            for (int i = 0; i < 4; i++)
            {
                singlefits[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                singlefits[i]->SetLineColor(colors[i]);
                singlefits[i]->SetLineStyle(2);
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
            ltemp->Draw("same");

            TLatex lat1;
            lat1.SetNDC();
            lat1.SetTextSize(0.03);
            lat1.SetTextFont(42);
            lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
            lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
            lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

            for (int i = 0; i < 4; i++)
            {
                double significance_num = singlefits[i]->Integral(masses[i] - 3 * widths[i], masses[i] + 3 * widths[i]) / binwidthfile;
                int binlow = hraw->GetXaxis()->FindBin(masses[i] - 3 * widths[i]);
                int binhigh = hraw->GetXaxis()->FindBin(masses[i] + 3 * widths[i]);
                double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
                significance = significance_num / significance_den;
                cout << "numerator " << significance_num << " denominator " << significance_den << endl;
                cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                double amplitude = obtained_parameters[3 * i];
                double amplitude_err = BEexpol->GetParError(3 * i);
                statSignificance = amplitude / amplitude_err;
                cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
            }
#endif

            // // // *******************************************************************************************
            // // // //********************************* for real and img part sum with interference ***************************************

#ifdef b_modifiedBoltzmann_hera

            TF1 *BEexpol = new TF1("BEexpol", BWsum_modifiedBoltzmann_hera, 1.05, 2.20, 14); // expol 3
            string parnames[] = {"f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a0", "a1", "a", "b", "c", "d"};
            for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
            {
                BEexpol->SetParName(i, parnames[i].c_str());
            }

            double parameters[] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, f1710Mass, f1710Width, 702, 801};
            int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

            for (int i = 0; i < size_fitparams; i++)
            {
                BEexpol->SetParameter(i, parameters[i]);
            }
            vector<vector<float>> par_limits = {{0, 3 * f1270Width}, {2, 3 * a1320Width}, {4, 3 * f1525Width}, {6, 3 * f1710Width}, {7, 50 * f1710WidthErr}};
            // BEexpol->SetParLimits(8, 670, 730);
            // BEexpol->SetParLimits(9, 770, 830);

            int limits_size = par_limits.size();
            for (int i = 0; i < limits_size; i++)
            {
                int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            }

            // //********systematic studies*************
            double initial_param_bkg[] = {2.64e5, -0.078, 2.562, 1.116};

            // Initial parameters for background
            BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Free
            BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Fix for medium train
            BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Free
            BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

            BEexpol->FixParameter(1, f1270Width);
            BEexpol->FixParameter(3, a1320Width);
            // BEexpol->FixParameter(3, 0.1736);
            BEexpol->FixParameter(5, f1525Width);

            // BEexpol->FixParameter(0, f1270Mass);
            // BEexpol->FixParameter(2, a1320Mass);
            // BEexpol->FixParameter(4, f1525Mass);

            // BEexpol->FixParameter(6, f1710Mass);
            // BEexpol->FixParameter(7, f1710Width);

            TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
            chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
            cout << "chi2/ndf is " << chi2ndf << endl;

            double *obtained_parameters = BEexpol->GetParameters();
            TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
            TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
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

            TF1 *onlyBW = new TF1("onlyBW", BWsum_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 10);
            TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 10);
            string parameter_names[] = {"norm1270", "mass1270", "width1270", "norm1320", "mass1320", "width1320", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710"};
            for (int i = 0; i < 10; i++)
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

            onlyBW_clone->FixParameter(2, f1270Width);
            onlyBW_clone->FixParameter(5, a1320Width);
            onlyBW_clone->FixParameter(8, f1525Width);

            // onlyBW_clone->FixParameter(1, f1270Mass);
            // onlyBW_clone->FixParameter(4, a1320Mass);
            // onlyBW_clone->FixParameter(7, f1525Mass);

            // onlyBW_clone->FixParameter(11, f1710Width);
            // onlyBW_clone->FixParameter(10, f1710Mass);

            onlyBW->SetLineColor(4);
            onlyBW->SetLineStyle(2);
            // onlyBW->Draw("same");

            // // Now plot the indivial resonances
            TF1 *singlefits[4];
            for (int i = 0; i < 4; i++)
            {
                singlefits[i] = new TF1(Form("singlef%d", i), single_BW_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                singlefits[i]->SetLineColor(colors[i]);
                singlefits[i]->SetLineStyle(2);
                // singlefits[i]->Draw("same");
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
            // ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
            // ltemp->AddEntry(singlefits[1], "a_{2}(1320)^{0}", "l");
            // ltemp->AddEntry(singlefits[2], "f'_{2}(1525)", "l");
            // ltemp->AddEntry(singlefits[3], "f_{0}(1710)", "l");
            ltemp->Draw("same");

            TLatex lat1;
            lat1.SetNDC();
            lat1.SetTextSize(0.03);
            lat1.SetTextFont(42);
            lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{#it{s}} = 13.6 TeV");
            lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
            lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/#it{c}", pT_bins[ipt], pT_bins[ipt + 1]));

            for (int i = 0; i < 4; i++)
            {
                double significance_num = singlefits[i]->Integral(masses[i] - 2 * widths[i], masses[i] + 2 * widths[i]) / binwidthfile;
                int binlow = hraw->GetXaxis()->FindBin(masses[i] - 2 * widths[i]);
                int binhigh = hraw->GetXaxis()->FindBin(masses[i] + 2 * widths[i]);
                double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
                significance = significance_num / significance_den;
                signal_counts = significance_num;
                background_counts = hraw->Integral(binlow, binhigh) - signal_counts;

                cout << "numerator " << significance_num << " denominator " << significance_den << endl;
                cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                double amplitude = obtained_parameters[3 * i];
                double amplitude_err = BEexpol->GetParError(3 * i);
                statSignificance = amplitude / amplitude_err;
                cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
            }
#endif

            // // // //******************************************************************************************
            // // // //********************************* BW sum with mass dependent width ***************************************
#ifdef b_modifiedBoltzmann_hera_mass_dep
            TF1 *BEexpol = new TF1("BEexpol", BWsum_ModifiedBoltzmann_hera_mass_dep, 1.05, 2.20, 14); // expol 3
            string parnames[] = {"f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a0", "a1", "a", "b", "c", "d"};
            for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
            {
                BEexpol->SetParName(i, parnames[i].c_str());
            }

            double parameters[] = {f1270Mass, f1270Width, a1320Mass, a1320Width, f1525Mass, f1525Width, f1710Mass, f1710Width, 560, 920};
            int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

            for (int i = 0; i < size_fitparams; i++)
            {
                BEexpol->SetParameter(i, parameters[i]);
            }
            vector<vector<float>> par_limits = {{0, 3 * f1270Width}, {2, 3 * a1320Width}, {4, 3 * f1525Width}, {6, 3 * f1525Width}, {7, 10 * f1710WidthErr}};

            // int limits_size = par_limits.size();
            // for (int i = 0; i < limits_size; i++)
            // {
            //     int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            //     BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            // }

            double initial_param_bkg[] = {3.2e5, -0.02, 2.802, 1.08}; // rotational 0-30 GeV/c (KsKs channel)

            // Initial parameters for background
            BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Fix
            BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Free
            BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Fix
            BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

            BEexpol->FixParameter(1, f1270Width);
            BEexpol->FixParameter(3, a1320Width);
            BEexpol->FixParameter(5, f1525Width);

            // BEexpol->FixParameter(0, f1270Mass);
            // BEexpol->FixParameter(2, a1320Mass);
            // BEexpol->FixParameter(4, f1525Mass);

            // BEexpol->FixParameter(6, f1710Mass);
            // BEexpol->FixParameter(7, f1710Width);

            TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
            chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
            cout << "chi2/ndf is " << chi2ndf << endl;

            double *obtained_parameters = BEexpol->GetParameters();
            TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
            TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
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

            TF1 *onlyBW = new TF1("onlyBW", BWsum_hera_mass_dep, BEexpol->GetXmin(), BEexpol->GetXmax(), 10);
            TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum_hera_mass_dep, BEexpol->GetXmin(), BEexpol->GetXmax(), 10);
            string parameter_names[] = {"norm1270", "mass1270", "width1270", "norm1320", "mass1320", "width1320", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710"};
            for (int i = 0; i < 10; i++)
            {
                onlyBW->SetParameter(i, obtained_parameters[i]);
                onlyBW_clone->SetParameter(i, obtained_parameters[i]);
                onlyBW_clone->SetParName(i, parameter_names[i].c_str());
            }
            // for (int i = 0; i < limits_size; i++)
            // {
            //     int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            //     onlyBW->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            //     onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            // }

            onlyBW_clone->FixParameter(2, f1270Width);
            onlyBW_clone->FixParameter(5, a1320Width);
            onlyBW_clone->FixParameter(8, f1525Width);

            // onlyBW_clone->FixParameter(1, f1270Mass);
            // onlyBW_clone->FixParameter(4, a1320Mass);
            // onlyBW_clone->FixParameter(7, f1525Mass);

            // onlyBW_clone->FixParameter(11, f1710Width);
            // onlyBW_clone->FixParameter(10, f1710Mass);

            onlyBW->SetLineColor(4);
            onlyBW->SetLineStyle(2);
            onlyBW->Draw("same");

            // // Now plot the indivial resonances
            TF1 *singlefits[4];
            for (int i = 0; i < 4; i++)
            {
                singlefits[i] = new TF1(Form("singlef%d", i), single_BW_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                singlefits[i]->SetLineColor(colors[i]);
                singlefits[i]->SetLineStyle(2);
                // singlefits[i]->Draw("same");
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
            // ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
            // ltemp->AddEntry(singlefits[1], "a_{2}(1320)^{0}", "l");
            // ltemp->AddEntry(singlefits[2], "f'_{2}(1525)", "l");
            // ltemp->AddEntry(singlefits[3], "f_{0}(1710)", "l");
            ltemp->Draw("same");

            TLatex lat1;
            lat1.SetNDC();
            lat1.SetTextSize(0.03);
            lat1.SetTextFont(42);
            lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{#it{s}} = 13.6 TeV");
            lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
            lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/#it{c}", pT_bins[ipt], pT_bins[ipt + 1]));

            for (int i = 0; i < 4; i++)
            {
                double significance_num = singlefits[i]->Integral(masses[i] - 2 * widths[i], masses[i] + 2 * widths[i]) / binwidthfile;
                int binlow = hraw->GetXaxis()->FindBin(masses[i] - 2 * widths[i]);
                int binhigh = hraw->GetXaxis()->FindBin(masses[i] + 2 * widths[i]);
                double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
                significance = significance_num / significance_den;
                signal_counts = significance_num;
                background_counts = hraw->Integral(binlow, binhigh) - signal_counts;

                cout << "numerator " << significance_num << " denominator " << significance_den << endl;
                cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                double amplitude = obtained_parameters[3 * i];
                double amplitude_err = BEexpol->GetParError(3 * i);
                statSignificance = amplitude / amplitude_err;
                cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
            }
#endif

            // // // ****************************************************************
            // // // *************For BW sum using real and imaginary parts without interference ******************

#ifdef b_modifiedBoltzmann_hera_const

            TF1 *BEexpol = new TF1("BEexpol", BWsum_modifiedBoltzmann_hera_const, 1.05, 2.20, 16); // expol 3

            string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
            for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
            {
                BEexpol->SetParName(i, parnames[i].c_str());
            }

            double parameters[] = {6700, f1270Mass, f1270Width, 800, a1320Mass, a1320Width, 7800, f1525Mass, f1525Width, 3800, f1710Mass, f1710Width};

            int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

            for (int i = 0; i < size_fitparams; i++)
            {
                BEexpol->SetParameter(i, parameters[i]);
            }
            vector<vector<double>> par_limits = {{1, 1 * f1270Width}, {4, 1 * a1320Width}, {7, 1 * f1525Width}, {10, 2 * f1710Width}, {11, 5 * f1710WidthErr}};

            int limits_size = par_limits.size();
            for (int i = 0; i < limits_size; i++)
            {
                int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            }

            // //********systematic studies*************
            double initial_param_bkg[] = {6.8e5, -0.0234, 3.071167, 1.00}; // rebin twice

            // Initial parameters for background
            BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Free
            BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Fix for medium train
            BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Free
            BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

            BEexpol->FixParameter(2, f1270Width);
            BEexpol->FixParameter(5, a1320Width);
            BEexpol->FixParameter(8, f1525Width);

            // BEexpol->FixParameter(1, f1270Mass);
            // BEexpol->FixParameter(4, a1320Mass);
            // BEexpol->FixParameter(7, f1525Mass);

            // BEexpol->FixParameter(10, f1710Mass);
            // BEexpol->FixParameter(11, f1710Width);

            TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
            chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
            cout << "chi2/ndf is " << chi2ndf << endl;

            double *obtained_parameters = BEexpol->GetParameters();
            TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
            TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
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

            TF1 *onlyBW = new TF1("onlyBW", BWsum_hera_const, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum_hera_const, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
            string parameter_names[] = {"norm1270", "mass1270", "width1270", "norm1320", "mass1320", "width1320", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710"};
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

            onlyBW_clone->FixParameter(2, f1270Width);
            onlyBW_clone->FixParameter(5, a1320Width);
            onlyBW_clone->FixParameter(8, f1525Width);

            // onlyBW_clone->FixParameter(1, f1270Mass);
            // onlyBW_clone->FixParameter(4, a1320Mass);
            // onlyBW_clone->FixParameter(7, f1525Mass);

            // onlyBW_clone->FixParameter(11, f1710Width);
            // onlyBW_clone->FixParameter(10, f1710Mass);

            onlyBW->SetLineColor(4);
            onlyBW->SetLineStyle(2);
            // onlyBW->Draw("same");

            // // Now plot the indivial resonances
            TF1 *singlefits[4];
            for (int i = 0; i < 4; i++)
            {
                singlefits[i] = new TF1(Form("singlef%d", i), single_BW_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                singlefits[i]->SetLineColor(colors[i]);
                singlefits[i]->SetLineStyle(2);
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
            ltemp->Draw("same");

            TLatex lat1;
            lat1.SetNDC();
            lat1.SetTextSize(0.03);
            lat1.SetTextFont(42);
            lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{#it{s}} = 13.6 TeV");
            lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
            lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/#it{c}", pT_bins[ipt], pT_bins[ipt + 1]));

            for (int i = 0; i < 4; i++)
            {
                double significance_num = singlefits[i]->Integral(masses[i] - 2 * widths[i], masses[i] + 2 * widths[i]) / binwidthfile;
                int binlow = hraw->GetXaxis()->FindBin(masses[i] - 2 * widths[i]);
                int binhigh = hraw->GetXaxis()->FindBin(masses[i] + 2 * widths[i]);
                double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
                significance = significance_num / significance_den;
                signal_counts = significance_num;
                background_counts = hraw->Integral(binlow, binhigh) - signal_counts;

                cout << "numerator " << significance_num << " denominator " << significance_den << endl;
                cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                double amplitude = obtained_parameters[3 * i];
                double amplitude_err = BEexpol->GetParError(3 * i);
                statSignificance = amplitude / amplitude_err;
                cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
            }
#endif

#ifdef b_coherentSum_modifiedBoltzmann
            TF1 *BEexpol = new TF1("BEexpol", CoherentSum_modifiedBoltzmann, 1.05, 2.20, 18); // expol 3
            string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "#phi1", "#phi2", "a", "b", "c", "d"};
            for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
            {
                BEexpol->SetParName(i, parnames[i].c_str());
            }

            // double parameters[] = {45, f1270Mass, f1270Width, 30, a1320Mass, a1320Width, 70, f1525Mass, f1525Width, 35, f1710Mass, f1710Width, TMath::Pi(), 0.0};
            double parameters[] = {45, f1270Mass, f1270Width, 29, a1320Mass, a1320Width, 65, f1525Mass, f1525Width, 45, f1710Mass, f1710Width, TMath::Pi(), 0.0};

            int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

            for (int i = 0; i < size_fitparams; i++)
            {
                BEexpol->SetParameter(i, parameters[i]);
            }
            vector<vector<double>> par_limits = {{1, 1 * f1270Width}, {4, 1 * a1320Width}, {7, 1 * f1525Width}, {10, 2 * f1710Width}, {11, 1 * f1710WidthErr}};

            int limits_size = par_limits.size();
            for (int i = 0; i < limits_size; i++)
            {
                int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
                BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            }

            // //********systematic studies*************
            // double initial_param_bkg[] = {5.3e5, -0.036, 2.702, 1.09}; // for modified Boltzmann
            double initial_param_bkg[] = {1.14e7, 0.0012, -2.83, -0.03}; // for expol

            // BEexpol->FixParameter(12, TMath::Pi()); // Fix phase 1
            // BEexpol->FixParameter(13, 0.0);         // Fix phase 2
            // BEexpol->SetParLimits(12, -TMath::Pi(), TMath::Pi()); // Set limits for phase 1
            // BEexpol->SetParLimits(13, -TMath::Pi(), TMath::Pi()); // Set limits for phase 2

            // Initial parameters for background
            BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Free
            BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Fix for medium train
            BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Free
            BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

            BEexpol->FixParameter(2, f1270Width);
            BEexpol->FixParameter(5, a1320Width);
            BEexpol->FixParameter(8, f1525Width);

            // BEexpol->FixParameter(1, f1270Mass);
            // BEexpol->FixParameter(4, a1320Mass);
            // BEexpol->FixParameter(7, f1525Mass);

            // BEexpol->FixParameter(10, f1710Mass);
            // BEexpol->FixParameter(11, f1710Width);

            TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
            chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
            cout << "chi2/ndf is " << chi2ndf << endl;
            phi_mod = TVector2::Phi_0_2pi(BEexpol->GetParameter(12));
            phi_mod2 = TVector2::Phi_0_2pi(BEexpol->GetParameter(13));
            cout << "Phase 1 is " << phi_mod << " Phase 2 is " << phi_mod2 << endl;

            double *obtained_parameters = BEexpol->GetParameters();
            // TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
            // TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
            TF1 *expol = new TF1("expol", expol_chkstar, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
            TF1 *expol_clone = new TF1("expol_clone", expol_chkstar, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
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

            TF1 *onlyBW = new TF1("onlyBW", coherent_sum, BEexpol->GetXmin(), BEexpol->GetXmax(), 14);
            TF1 *onlyBW_clone = new TF1("onlyBW_clone", coherent_sum, BEexpol->GetXmin(), BEexpol->GetXmax(), 14);
            string parameter_names[] = {"norm1270", "mass1270", "width1270", "norm1320", "mass1320", "width1320", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710", "phase 1", "phase 2"};
            for (int i = 0; i < 14; i++)
            {
                onlyBW->SetParameter(i, obtained_parameters[i]);
                onlyBW_clone->SetParameter(i, obtained_parameters[i]);
                onlyBW_clone->SetParName(i, parameter_names[i].c_str());
            }
            // for (int i = 0; i < limits_size; i++)
            // {
            //     int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            //     onlyBW->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            //     onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            // }

            onlyBW_clone->FixParameter(2, f1270Width);
            onlyBW_clone->FixParameter(5, a1320Width);
            onlyBW_clone->FixParameter(8, f1525Width);

            // onlyBW_clone->FixParameter(1, f1270Mass);
            // onlyBW_clone->FixParameter(4, a1320Mass);
            // onlyBW_clone->FixParameter(7, f1525Mass);

            // onlyBW_clone->FixParameter(11, f1710Width);
            // onlyBW_clone->FixParameter(10, f1710Mass);

            onlyBW->SetLineColor(4);
            onlyBW->SetLineStyle(2);
            onlyBW->Draw("same");

            // // Now plot the indivial resonances
            TF1 *singlefits[4];
            for (int i = 0; i < 4; i++)
            {
                singlefits[i] = new TF1(Form("singlef%d", i), single_BW_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
                singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
                singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
                singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
                singlefits[i]->SetLineColor(colors[i]);
                singlefits[i]->SetLineStyle(2);
                // singlefits[i]->Draw("same");
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
            // ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
            // ltemp->AddEntry(singlefits[1], "a_{2}(1320)^{0}", "l");
            // ltemp->AddEntry(singlefits[2], "f'_{2}(1525)", "l");
            // ltemp->AddEntry(singlefits[3], "f_{0}(1710)", "l");
            ltemp->Draw("same");

            TLatex lat1;
            lat1.SetNDC();
            lat1.SetTextSize(0.03);
            lat1.SetTextFont(42);
            lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{#it{s}} = 13.6 TeV");
            lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
            lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/#it{c}", pT_bins[ipt], pT_bins[ipt + 1]));

            for (int i = 0; i < 4; i++)
            {
                double significance_num = singlefits[i]->Integral(masses[i] - 2 * widths[i], masses[i] + 2 * widths[i]) / binwidthfile;
                int binlow = hraw->GetXaxis()->FindBin(masses[i] - 2 * widths[i]);
                int binhigh = hraw->GetXaxis()->FindBin(masses[i] + 2 * widths[i]);
                double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
                significance = significance_num / significance_den;
                signal_counts = significance_num;
                background_counts = hraw->Integral(binlow, binhigh) - signal_counts;

                cout << "numerator " << significance_num << " denominator " << significance_den << endl;
                cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
                double amplitude = obtained_parameters[3 * i];
                double amplitude_err = BEexpol->GetParError(3 * i);
                statSignificance = amplitude / amplitude_err;
                cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
            }
#endif

            // // // //******************************************************************************************
            // // // //********************************* common for all fits ***************************************
            gPad->Update();
            TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
            ptstats->SetX1NDC(0.6);
            ptstats->SetX2NDC(0.99);
            ptstats->SetY1NDC(0.4);
            ptstats->SetY2NDC(0.92);
            ptstats->Draw("same");
            // c->SaveAs((savepath + Form("/rBWfit_pt_%.2f_%.2f_%s.png", pT_bins[ipt], pT_bins[ipt + 1], sysvar.c_str())).c_str());
            c->SaveAs((savepath + "/rBWfit.png").c_str());

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
            double fitrangelow = BEexpol->GetXmin();
            double fitrangehigh = BEexpol->GetXmax();
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

            // // Write data (values are separated by commas) for .csv file
            // file << "Norm range," << kNormRangepT[0][0] << " - " << kNormRangepT[0][1] << "\n";
            // file << "Fit range," << fitrangelow << " - " << fitrangehigh << "\n";
            // file << "Signal counts," << signal_counts << "\n";
            // file << "Background counts," << background_counts << "\n";
            // file << "S/B ratio," << (signal_counts / background_counts) << "\n";
            // file << std::fixed << std::setprecision(2);
            // file << "Significance," << significance << "\n";
            // file << "StatSignificance," << statSignificance << "\n";
            // file << "Chi2NDF," << (BEexpol->GetChisquare() / BEexpol->GetNDF()) << "\n";

            // file << std::fixed << std::setprecision(1);
            // // Fit parameters with ± sign
            // file << "Fit mass," << fitmass1710 * 1000 << " ± " << fitmass1710_err * 1000 << "\n";
            // file << "Fit width," << fitwidth1710 * 1000 << " ± " << fitwidth1710_err * 1000 << "\n";
            // file << "Fit norm," << fitnorm1710 << " ± " << fitnorm1710_err << "\n";

            // // // for .txt files for systematics
            // file << "Norm range " << kNormRangepT[0][0] << " - " << kNormRangepT[0][1] << endl;
            // file << "Fit range " << fitrangelow << " - " << fitrangehigh << endl;
            // file << "Significance " << significance << endl;
            // file << "StatSignificance " << statSignificance << endl;
            // file << "Fit parameters of f1710 " << endl;
            // file << "Chi2NDF " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
            // file << fitnorm1710 << " ± " << fitnorm1710_err << endl;
            // file << fitmass1710 << " ± " << fitmass1710_err << endl;
            // file << fitwidth1710 << " ± " << fitwidth1710_err << endl;

            // for table making
            file << std::fixed << std::setprecision(2);
            file << chi2ndf << endl;
            file << std::fixed << std::setprecision(0);
            file << fitmass1270 * 1000 << " ± " << fitmass1270_err * 1000 << endl;
            file << fitmass1320 * 1000 << " ± " << fitmass1320_err * 1000 << endl;
            file << fitmass1525 * 1000 << " ± " << fitmass1525_err * 1000 << endl;
            file << fitmass1710 * 1000 << " ± " << fitmass1710_err * 1000 << endl;
            file << fitwidth1710 * 1000 << " ± " << fitwidth1710_err * 1000 << endl;

#ifdef residual_subtracted
            // Now subtract the residual background and plot
            TCanvas *c2 = new TCanvas("", "", 720, 720);
            SetCanvasStyle(c2, 0.14, 0.03, 0.05, 0.14);
            expol_clone->SetRange(0.99, 2.99);
            hsubtracted->Add(expol_clone, -1);
            hsubtracted->GetXaxis()->SetRangeUser(BEexpol->GetXmin(), BEexpol->GetXmax());
            hsubtracted->SetMaximum(hsubtracted->GetMaximum() * 1.7);
            hsubtracted->Draw();
            hsubtracted->Write("3BW");
            // for (int i = 0; i < limits_size; i++)
            // {
            //     int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            //     onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            // }
            onlyBW_clone->SetNpx(1000);
            // onlyBW_clone->SetLineColor(1);
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
                singlefits1[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, onlyBW_clone->GetXmin(), onlyBW_clone->GetXmax(), 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, onlyBW_clone->GetXmin(), onlyBW_clone->GetXmax(), 3);
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

            TLegend *ltemp2 = new TLegend(0.20, 0.58, 0.42, 0.90);
            ltemp2->SetFillStyle(0);
            ltemp2->SetTextFont(42);
            ltemp2->SetBorderSize(0);
            ltemp2->SetTextSize(0.035);
            ltemp2->SetHeader("Residual bkg subtraction");
            ltemp2->AddEntry(hsubtracted, "Data (stat. uncert.)", "lpe");
            ltemp2->AddEntry(onlyBW_clone, "4rBW fit", "l");
            ltemp2->AddEntry(singlefits1[0], "f_{2}(1270)", "l");
            ltemp2->AddEntry(singlefits1[1], "a_{2}(1320)^{0}", "l");
            ltemp2->AddEntry(singlefits1[2], "f'_{2}(1525)", "l");
            ltemp2->AddEntry(singlefits1[3], "f_{0}(1710)", "l");
            ltemp2->Draw("same");
            TLatex lat2;
            lat2.SetNDC();
            lat2.SetTextSize(0.03);
            lat2.SetTextFont(42);
            // lat2.DrawLatex(0.20, 0.89, "KsKs unlike sign pairs");
            // lat2.DrawLatex(0.20, 0.85, "Residual background subtracted");
            c2->SaveAs((savepath + "/rBWfit_residual_" + sysvar + ".png").c_str());
            // c2->Close();
#endif

            TMatrixDSym cov = fitResultptr->GetCovarianceMatrix();
            TMatrixDSym cov1;
            TMatrixDSym cov2;
            cov.GetSub(0, 11, 0, 11, cov1);
            cov.GetSub(12, 15, 12, 15, cov2);
            Double_t *b = cov1.GetMatrixArray();
            Double_t *a = cov2.GetMatrixArray();
            Double_t *para = onlyBW_clone->GetParameters();

            float nsigma_yield = 3.0;

            cout << "Function integration for f0(1710) " << singlefits[3]->Integral(f1710Mass - nsigma_yield * f1710Width, f1710Mass + nsigma_yield * f1710Width)<<endl;
            cout << " Function integration for f2(1525) " << singlefits[2]->Integral(f1525Mass - nsigma_yield * f1525Width, f1525Mass + nsigma_yield * f1525Width)<<endl;

            // // Yield calculation
            // double yield1270 = singlefits[0]->Integral(f1270Mass - nsigma_yield * f1270Width, f1270Mass + nsigma_yield * f1270Width) / (binwidthfile * total_events);
            // double yield1320 = singlefits[1]->Integral(a1320Mass - nsigma_yield * a1320Width, a1320Mass + nsigma_yield * a1320Width) / (binwidthfile * total_events);
            double yield1525 = singlefits[2]->Integral(f1525Mass - nsigma_yield * f1525Width, f1525Mass + nsigma_yield * f1525Width) / (binwidthfile * total_events);
            double yield1710 = singlefits[3]->Integral(f1710Mass - nsigma_yield * f1710Width, f1710Mass + nsigma_yield * f1710Width) / (binwidthfile * total_events);

            // double yield1270 = singlefits[0]->Integral(f1270Mass - nsigma_yield * f1270Width, f1270Mass + nsigma_yield * f1270Width);
            // double yield1320 = singlefits[1]->Integral(a1320Mass - nsigma_yield * a1320Width, a1320Mass + nsigma_yield * a1320Width);
            // double yield1525 = singlefits[2]->Integral(f1525Mass - nsigma_yield * f1525Width, f1525Mass + nsigma_yield * f1525Width);
            // double yield1710 = singlefits[3]->Integral(f1710Mass - nsigma_yield * f1710Width, f1710Mass + nsigma_yield * f1710Width);

            // // Yield calculation after residual background subtraction
            // double yield1270_resSub = singlefits1[0]->Integral(f1270Mass - nsigma_yield * f1270Width, f1270Mass + nsigma_yield * f1270Width) / (binwidthfile * total_events);
            // double yield1320_resSub = singlefits1[1]->Integral(a1320Mass - nsigma_yield * a1320Width, a1320Mass + nsigma_yield * a1320Width) / (binwidthfile * total_events);
            // double yield1525_resSub = singlefits1[2]->Integral(f1525Mass - nsigma_yield * f1525Width, f1525Mass + nsigma_yield * f1525Width) / (binwidthfile * total_events);
            // double yield1710_resSub = singlefits1[3]->Integral(f1710Mass - nsigma_yield * f1710Width, f1710Mass + nsigma_yield * f1710Width) / (binwidthfile * total_events);

            // // Yield errors calculation
            // double yield1270_err = onlyBW_clone->IntegralError((f1270Mass - nsigma_yield * f1270Width), (f1270Mass + nsigma_yield * f1270Width), &para[0], b) / (binwidthfile * total_events);
            // double yield1320_err = onlyBW_clone->IntegralError((a1320Mass - nsigma_yield * a1320Width), (a1320Mass + nsigma_yield * a1320Width), &para[0], b) / (binwidthfile * total_events);
            double yield1525_err = onlyBW_clone->IntegralError((f1525Mass - nsigma_yield * f1525Width), (f1525Mass + nsigma_yield * f1525Width), &para[0], b) / (binwidthfile * total_events);
            double yield1710_err = onlyBW_clone->IntegralError((f1710Mass - nsigma_yield * f1710Width), (f1710Mass + nsigma_yield * f1710Width), &para[0], b) / (binwidthfile * total_events);

            cout << "Total event count is " << total_events << endl;
            // cout << "Yield of f1270 is " << yield1270 << " ± " << yield1270_err << endl;
            // cout << "Yield of a1320 is " << yield1320 << " ± " << yield1320_err << endl;
            cout << "Yield of f1525 is " << yield1525 << " ± " << yield1525_err << endl;
            cout << "Yield of f1710 is " << yield1710 << " ± " << yield1710_err << endl;
            // cout << "Total yield from function integration is " << yield1270 + yield1320 + yield1525 + yield1710 << endl;

            // cout << "Total yield after residual background subtraction " << yield1270_resSub + yield1320_resSub + yield1525_resSub + yield1710_resSub << endl;
            // cout << "Energy bin width is " << binwidthfile << endl;

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
            hinvMass->SetMinimum(-40000);
            // hinvMass->SetMinimum(-4000);
            // hinvMass->SetMaximum(0.9e6);
            hinvMass->SetMaximum(hinvMass->GetMaximum() * 1.0);
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
            leg->SetHeader("ALICE Performance");
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
            lat5.DrawLatex(0.32, 0.73, "FT0M 0-100%, |#it{y}| < 0.5");
            lat5.DrawLatex(0.32, 0.65, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", pT_bins[ipt], pT_bins[ipt + 1]));

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
            hsubtracted->SetMaximum(hsubtracted->GetMaximum() * 0.67);
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
            // lat6.DrawLatex(0.65, 0.60, "FT0M (0-100%), |#it{y}|<0.5");
            // lat6.DrawLatex(0.65, 0.50, Form("%.1f < #it{p}_{T} #leq %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

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
            c1->SaveAs((savepath + Form("/rBWfit_doublepanel_%s.eps", sysvar.c_str())).c_str());

#endif

            // **********************************************************************************************
            // *******************subtract the resonance peaks and fit the residual background*****************

            // // TCanvas *c3 = new TCanvas("", "", 720, 720);
            // // // Here we will subtract the resonances peaks and plot the residual background and the fit it
            // // SetCanvasStyle(c3, 0.14, 0.03, 0.05, 0.14);
            // // onlyBW->SetRange(f1270Mass - 3 * f1270Width, f1710Mass + 3 * f1710Width);
            // // hsubtracted_res->Add(onlyBW, -1);
            // // hsubtracted_res->Draw();
            // // hsubtracted_res->Fit("expol", "REBMSI");

            // // // Yield calculation
            // // double yield1270 = onlyBW_clone->Integral(f1270Mass - 2 * f1270Width, f1270Mass + 2 * f1270Width);
            // // double yield1320 = onlyBW_clone->Integral(a1320Mass - 2 * a1320Width, a1320Mass + 2 * a1320Width);
            // // double yield1525 = onlyBW_clone->Integral(f1525Mass - 2 * f1525Width, f1525Mass + 2 * f1525Width);
            // // double yield1710 = onlyBW_clone->Integral(f1710Mass - 2 * f1710Width, f1710Mass + 2 * f1710Width);

            // // double yield1270_err = onlyBW_clone->IntegralError((f1270Mass - 3 * f1270Width), (f1270Mass + 3 * f1270Width));
            // // double yield1320_err = onlyBW_clone->IntegralError((a1320Mass - 3 * a1320Width), (a1320Mass + 3 * a1320Width));
            // // double yield1525_err = onlyBW_clone->IntegralError((f1525Mass - 3 * f1525Width), (f1525Mass + 3 * f1525Width));
            // // double yield1710_err = onlyBW_clone->IntegralError((f1710Mass - 3 * f1710Width), (f1710Mass + 3 * f1710Width));

            // // cout << "Yield 1270: " << yield1270 << " +- " << yield1270_err << endl;
            // // cout << "Yield 1320: " << yield1320 << " +- " << yield1320_err << endl;
            // // cout << "Yield 1525: " << yield1525 << " +- " << yield1525_err << endl;
            // // cout << "Yield 1710: " << yield1710 << " +- " << yield1710_err << endl;
        }
        // plots_4BW->Close();
    }
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

    // cout << "x[0] is " << x[0] << endl;
    // cout<< "mass is "<<mass<<endl;
    // cout<<"num is "<<num<<endl;
    // cout<<"den is "<<den<<endl;

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

Double_t BWsum(double *x, double *par)
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

Double_t BWsum_expol3(double *x, double *par)
{
    return (BWsum(x, par) + exponential_bkg_3(x, &par[12]));
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
    return (BWsum(x, par) + Boltzmann_bkg_2(x, &par[12]));
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
    return (coherent_sum(x, par) + expol_chkstar(x, &par[14]));
    // return (coherent_sum(x, par) + exponential_bkg_3(x, &par[14]));
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