#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
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
Double_t BWsum_boltzmann_massdep(double *x, double *par);
Double_t single_BW_mass_dep_spin0(double *x, double *par);
Double_t single_BW_mass_dep_spin2(double *x, double *par);

Double_t exponential_bkg_1(double *x, double *par); // 3 parameters
Double_t exponential_bkg_2(double *x, double *par); // 3 parameters
Double_t exponential_bkg_3(double *x, double *par); // 4 parameters
Double_t exponential_bkg_4(double *x, double *par); // 5 parameters
Double_t exponential_bkg_5(double *x, double *par); // 3 parameters
Double_t exponential_bkg_6(double *x, double *par); // 4 parameters
Double_t expol_chkstar(double *x, double *par);     // 4 parameters

Double_t Boltzmann_bkg_1(double *x, double *par); // 3 parameters
Double_t Boltzmann_bkg_2(double *x, double *par); // 4 parameters

Double_t single_BW_expol3(double *x, double *par);
Double_t single_BW_expol3_hera(double *x, double *par);
Double_t BWsum_expol3(double *x, double *par);
Double_t BWsum_expol3_hera(double *x, double *par);
Double_t BWsum_expol_chkstar(double *x, double *par);
Double_t BWsumMassDepWidth_expol3(double *x, double *par);

Double_t single_BW_boltzman_1(double *x, double *par);
Double_t single_BW_boltzman_2(double *x, double *par);
Double_t BWsum_boltzman_1(double *x, double *par);
Double_t BWsum_boltzman_2(double *x, double *par);

void glueball_fit_3rBW()
{

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
    //  const string kvariation1 = "_TPCPID4_id24937"; //TPC PID
    //  const string kvariation1 = "_TPCPID6_id24937"; //TPC PID
    //  const string kvariation1 = "_TPCcr100_id24937"; //TPC crossed rows
    //  const string kvariation1 = "_TPCcr120_id24937"; //TPC crossed rows
    //  const string kvariation1 = "_TPCcrfc0p9_id24940"; //TPC crossed rows over findable clusters
    //  const string kvariation1 = "_TPCcrfc1p0_id24940"; //TPC crossed rows over findable clusters

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
    //  const string kvariation1 = "_Ks_selection4_id24939";
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
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances" + kvariation1;
    string path2 = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937";

    // //*********for temporary study with angular separation cuts************************
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/362701/KsKs_Channel/higher-mass-resonances" + kvariation1;
    // string path2 = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/362701/KsKs_Channel/higher-mass-resonances_id24794";

    string sysvar = "varC2"; // default

    ofstream file;
    file.open((path2 + "/fits/4rBw_fits/fit_params_" + sysvar + ".txt").c_str());

    string savepath = path2 + "/fits/4rBw_fits";

    gSystem->Exec(("mkdir -p " + savepath).c_str());

    TFile *f = new TFile((path + "/hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ"); // default
    // TFile *f = new TFile((path + "/hglue_MIX_norm_2.50_2.60_pt_3.00_5.00.root").c_str(), "READ"); // default
    // TFile *f = new TFile((path + "/hglue_ROTATED_norm_2.50_2.60_all_pT.root").c_str(), "READ"); //

    int colors[] = {4, 6, 28, 46};
    double masses[] = {f1270Mass, f1525Mass, f1710Mass};
    double widths[] = {f1270Width, f1525Width, f1710Width};
    string resonance_names[] = {"f_{2}(1270)", "f'_{2}(1525)", "f_{0}(1710)"};
    double purity, significance, chi2ndf, chi2, ndf, statSignificance;
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);
    t2->SetTextFont(42);

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

// #define b_modfied_boltzmann
#define b_massdepwidth
    // #define b_expol
    // #define b_boltzman_pure
    // #define b_boltzmann_pure_massdepwidth
// #define residual_subtracted
    // #define doublepanelplot

    for (int ipt = 0; ipt < Npt; ipt++)
    {

        // TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", 0.0, 30.0));
        TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", pT_bins[ipt], pT_bins[ipt + 1]));
        TH1F *hraw = (TH1F *)f->Get(Form("ksks_invmass_pt_%.1f_%.1f", pT_bins[ipt], pT_bins[ipt + 1]));
        if (hinvMass == nullptr)
        {
            cout << "Error opening histogram" << endl;
            return;
        }

        if (hinvMass == nullptr)
        {
            cout << "Error opening histogram" << endl;
            return;
        }
        TCanvas *c = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
        hinvMass->Rebin(2);
        double binwidthfile = (hinvMass->GetXaxis()->GetXmax() - hinvMass->GetXaxis()->GetXmin()) / hinvMass->GetXaxis()->GetNbins();
        hinvMass->GetXaxis()->SetRangeUser(1.00, 2.50);
        hinvMass->GetXaxis()->SetTitle("M_{K^{0}_{s}K^{0}_{s}} (GeV/c^{2})");
        hinvMass->GetYaxis()->SetTitle(Form("Counts/%.2f GeV/c^{2}", binwidthfile));
        hinvMass->Draw();
        // t2->DrawLatex(0.29, 0.96, Form("%.1f < #it{p}_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));
        TH1F *hsubtracted = (TH1F *)hinvMass->Clone("hsubtracted");
        TH1F *hsubtracted_res = (TH1F *)hinvMass->Clone("hsubtracted_res");

        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1111);
        vector<tuple<float, int, float, float>> fit_parameters;

// // //************************************************************************ */
// // // **************** For BW sum with expol HERA ****************************

// // Default fitting range is 1.02 to 2.20. Four types of fitting range variations: extend left (1.0), extend right (2.50), large range (1.0 to 2.50), small range (1.05 to 2.15)
#ifdef b_massdepwidth
        // TF1 *BEexpol = new TF1("BEexpol", BWsum_expol3, 1.05, 2.20, 13); // expol 3
        TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_expol3, 1.05, 2.20, 13); // expol 3
        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        double parameters[] = {8384, f1270Mass, f1270Width, 7858, f1525Mass, f1525Width, 3218, f1710Mass, f1710Width};
        // double parameters[] = {2500, f1270Mass, f1270Width, 7000, f1525Mass, f1525Width, 2500, f1710Mass, f1710Width}; // for preview presentation
        // double parameters[] = {1600, f1270Mass, f1270Width, 2200, f1525Mass, f1525Width, 1000, f1710Mass, f1710Width}; // 3-5 GeV/c rebin 2
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {4, 3 * f1525Width}, {10, 3 * f1710Width}};
        // vector<vector<float>> par_limits = {{2, 10 * f1270WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        // for rotational bkg with pt range 0-30 GeV/c , // 3BW/3BWAmp + expol/exponential1, 3BWAmp + expol, 3BWAmp + exponential3
        double initial_param_bkg[] = {4.1562e5, -0.09379, 3.869, 1.30982}; // rotational 1-30 GeV/c (KsKs channel)
        // double initial_param_bkg[] = {1.062e5, -0.09, 3.69, 1.572}; // 3-5 GeV/c rebin twice

        BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // 206 //5.845e5
        BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  // 0.04316 //-0.07378
        BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // 11.48 //2.685
        BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // -3.149 //1.176

        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(4, f1525Mass);

        BEexpol->FixParameter(2, f1270Width);
        BEexpol->FixParameter(5, f1525Width);

        // BEexpol->FixParameter(7, f1710Mass);
        // BEexpol->FixParameter(8, f1710Width);

        // BEexpol->FixParameter(1, a1320Mass);
        // BEexpol->FixParameter(2, a1320Width);

        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        // status codes: 4000 successful, 4 call limit, 4910 failed
        string fitstatus = "Successfull";
        if (fitResultptr->Status() != 4000)
        {
            cout << "Fit failed or call limit" << endl;
            fitstatus = "Failed";
        }
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;

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

        // TF1 *onlyBW = new TF1("onlyBW", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        // TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW = new TF1("onlyBW", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        onlyBW_clone->SetParNames("f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma");
        for (int i = 0; i < 9; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        // onlyBW->Draw("same");
        onlyBW_clone->FixParameter(2, f1270Width);
        onlyBW_clone->FixParameter(5, f1525Width);

        // TLegend *ltemp = new TLegend(0.20, 0.67, 0.52, 0.92);
        // ltemp->SetFillStyle(0);
        // ltemp->SetTextFont(42);
        // ltemp->SetTextSize(0.035);
        // ltemp->AddEntry(hinvMass, "Data", "lpe");
        // ltemp->AddEntry(BEexpol, "3rBW + expol", "l");
        // ltemp->AddEntry(onlyBW, "3rBW", "l");
        // ltemp->AddEntry(expol, "expol", "l");
        // ltemp->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits[3];
        for (int i = 0; i < 3; i++)
        {
            singlefits[i] = (i < 2) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
            singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
            singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
            singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
            singlefits[i]->SetLineColor(colors[i]);
            singlefits[i]->SetLineStyle(2);
            singlefits[i]->SetLineWidth(2);
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
        ltemp->AddEntry(BEexpol, "3rBW + Residual BG", "l");
        ltemp->AddEntry(expol, "Residual BG", "l");
        ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
        ltemp->AddEntry(singlefits[1], "f'_{2}(1525)", "l");
        ltemp->AddEntry(singlefits[2], "f_{0}(1710)", "l");
        ltemp->Draw("same");

        TLatex lat1;
        lat1.SetNDC();
        lat1.SetTextSize(0.03);
        lat1.SetTextFont(42);
        lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
        lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
        lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

        for (int i = 0; i < 3; i++)
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
#ifdef b_boltzman_pure

        TF1 *BEexpol = new TF1("BEexpol", BWsum_boltzman_1, 1.05, 2.20, 12); // Boltzmann
        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        // double parameters[] = {100, f1270Mass, f1270Width, 100, f1525Mass, f1525Width, 50, f1710Mass, f1710Width};
        double parameters[] = {6000, f1270Mass, f1270Width, 8000, f1525Mass, f1525Width, 4000, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {2, 10 * f1270WidthErr}, {4, 5 * f1525Width}, {7, 20 * f1710WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        BEexpol->SetParameter(size_fitparams + 0, 7.0618e5); // expol 3  // 7.618e5  // 5.562e5
        BEexpol->SetParameter(size_fitparams + 1, 0.6756);   // expol 3  //  0.6456  // -0.09379
        BEexpol->SetParameter(size_fitparams + 2, 4.238);    // expol 3  //4.238  // 2.569

        // // BEexpol->FixParameter(0, 6998);
        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(2, f1270Width);
        // // BEexpol->FixParameter(3, 7128);
        // BEexpol->FixParameter(4, f1525Mass);
        // BEexpol->FixParameter(5, f1525Width);
        // // BEexpol->FixParameter(6, 4058);
        // BEexpol->FixParameter(7, f1710Mass);
        // BEexpol->FixParameter(8, f1710Width);

        // BEexpol->FixParameter(1, a1320Mass);
        // BEexpol->FixParameter(2, a1320Width);

        hinvMass->Fit("BEexpol", "REMBS");

        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;

        double *obtained_parameters = BEexpol->GetParameters();
        TF1 *expol = new TF1("expol", Boltzmann_bkg_1, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
        TF1 *expol_clone = new TF1("expol_clone", Boltzmann_bkg_1, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
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

        TF1 *onlyBW = new TF1("onlyBW", BWsum_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        onlyBW_clone->SetParNames("norm1270", "mass1270", "width1270", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710");
        for (int i = 0; i < 9; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        // onlyBW->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits[3];
        for (int i = 0; i < 3; i++)
        {
            singlefits[i] = new TF1(Form("singlef%d", i), single_BW, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
            singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
            singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
            singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
            singlefits[i]->SetLineColor(colors[i]);
            singlefits[i]->SetLineStyle(2);
            singlefits[i]->SetLineWidth(2);
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
        ltemp->AddEntry(BEexpol, "3rBW + Residual BG", "l");
        ltemp->AddEntry(expol, "Residual BG", "l");
        ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
        ltemp->AddEntry(singlefits[1], "f'_{2}(1525)", "l");
        ltemp->AddEntry(singlefits[2], "f_{0}(1710)", "l");
        ltemp->Draw("same");

        TLatex lat1;
        lat1.SetNDC();
        lat1.SetTextSize(0.03);
        lat1.SetTextFont(42);
        lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
        lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
        lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

        for (int i = 0; i < 3; i++)
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

        // // // ************************************************************************************
        // // // ***** For BW sum with Boltzmann with mass dependent width ****************************

#ifdef b_boltzmann_pure_massdepwidth

        TF1 *BEexpol = new TF1("BEexpol", BWsum_boltzmann_massdep, 1.05, 2.20, 12); // Boltzmann
        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        // double parameters[] = {100, f1270Mass, f1270Width, 100, f1525Mass, f1525Width, 50, f1710Mass, f1710Width};
        double parameters[] = {10000, f1270Mass, f1270Width, 8000, f1525Mass, f1525Width, 5000, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {2, 10 * f1270WidthErr}, {4, 5 * f1525Width}, {7, 20 * f1710WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        BEexpol->SetParameter(size_fitparams + 0, 7.0618e5); // expol 3  // 7.618e5  // 5.562e5
        BEexpol->SetParameter(size_fitparams + 1, 0.6756);   // expol 3  //  0.6456  // -0.09379
        BEexpol->SetParameter(size_fitparams + 2, 4.238);    // expol 3  //4.238  // 2.569

        // // BEexpol->FixParameter(0, 6998);
        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(2, f1270Width);
        // // BEexpol->FixParameter(3, 7128);
        // BEexpol->FixParameter(4, f1525Mass);
        // BEexpol->FixParameter(5, f1525Width);
        // // BEexpol->FixParameter(6, 4058);
        // BEexpol->FixParameter(7, f1710Mass);
        // BEexpol->FixParameter(8, f1710Width);

        // BEexpol->FixParameter(1, a1320Mass);
        // BEexpol->FixParameter(2, a1320Width);

        hinvMass->Fit("BEexpol", "REMBS");

        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;

        double *obtained_parameters = BEexpol->GetParameters();
        TF1 *expol = new TF1("expol", Boltzmann_bkg_1, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);             //
        TF1 *expol_clone = new TF1("expol_clone", Boltzmann_bkg_1, BEexpol->GetXmin(), BEexpol->GetXmax(), 4); //
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

        TF1 *onlyBW = new TF1("onlyBW", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        onlyBW_clone->SetParNames("norm1270", "mass1270", "width1270", "norm12525", "mass1525", "width1525", "norm1710", "mass1710", "width1710");
        for (int i = 0; i < 9; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        // onlyBW->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits[3];
        for (int i = 0; i < 3; i++)
        {
            singlefits[i] = (i < 2) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
            singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
            singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
            singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
            singlefits[i]->SetLineColor(colors[i]);
            singlefits[i]->SetLineStyle(2);
            singlefits[i]->SetLineWidth(2);
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
        ltemp->AddEntry(BEexpol, "3rBW + Residual BG", "l");
        ltemp->AddEntry(expol, "Residual BG", "l");
        ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
        ltemp->AddEntry(singlefits[1], "f'_{2}(1525)", "l");
        ltemp->AddEntry(singlefits[2], "f_{0}(1710)", "l");
        ltemp->Draw("same");

        TLatex lat1;
        lat1.SetNDC();
        lat1.SetTextSize(0.03);
        lat1.SetTextFont(42);
        lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
        lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
        lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

        for (int i = 0; i < 3; i++)
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
        // // // //************************************************************************ */
        // // // // **************** For BW sum with exp + pol2 as used in Charged kstar **************************
#ifdef b_modfied_boltzmann
        TF1 *BEexpol = new TF1("BEexpol", BWsum_expol3, 1.05, 2.20, 13); // expol 1
        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        double parameters[] = {4000, f1270Mass, f1270Width, 3758, f1525Mass, f1525Width, 1500, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }

        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {4, 5 * f1525Width}, {7, 20 * f1710WidthErr}};
        // vector<vector<float>> par_limits = {{1, 3 * a1320Mass}, {4, 5 * f1525Width}, {7, 20 * f1710WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        double initial_param_bkg[] = {3.562e5, -0.009379, 2.9, 1.0982}; // rotational 1-30 GeV/c (KsKs channel)

        // for rotational bkg with pt range 0-30 GeV/c (KsKs channel)
        BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Fix
        BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Free
        BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Fix
        BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

        BEexpol->FixParameter(2, f1270Width);
        BEexpol->FixParameter(5, f1525Width);

        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(4, f1525Mass);

        // BEexpol->FixParameter(7, f1710Mass);
        // BEexpol->FixParameter(8, f1710Width);

        // BEexpol->FixParameter(1, a1320Mass);
        // BEexpol->FixParameter(2, a1320Width);

        // hinvMass->Fit("BEexpol", "REBMS");
        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        // status codes: 4000 successful, 4 call limit, 4910 failed
        string fitstatus = "Successfull";
        // cout<<"fit status code "<<fitResultptr->Status()<<endl;
        if (fitResultptr->Status() != 4140)
        {
            cout << "Fit failed or call limit !!!!!!!" << endl;
            fitstatus = "Failed";
        }
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
        chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
        chi2 = BEexpol->GetChisquare();
        ndf = BEexpol->GetNDF();

        double *obtained_parameters = BEexpol->GetParameters();
        TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
        TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
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

        string par_names2[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma"};

        TF1 *onlyBW = new TF1("onlyBW", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        for (int i = 0; i < 9; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
            onlyBW->SetParName(i, par_names2[i].c_str());
            onlyBW_clone->SetParName(i, par_names2[i].c_str());
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        // onlyBW->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits[3];
        for (int i = 0; i < 3; i++)
        {
            singlefits[i] = new TF1(Form("singlef%d", i), single_BW, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
            singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
            singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
            singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
            singlefits[i]->SetLineColor(colors[i]);
            singlefits[i]->SetLineStyle(2);
            singlefits[i]->SetLineWidth(2);
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
        ltemp->AddEntry(BEexpol, "3rBW + Residual BG", "l");
        ltemp->AddEntry(expol, "Residual BG", "l");
        ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
        ltemp->AddEntry(singlefits[1], "f'_{2}(1525)", "l");
        ltemp->AddEntry(singlefits[2], "f_{0}(1710)", "l");
        ltemp->Draw("same");

        TLatex lat1;
        lat1.SetNDC();
        lat1.SetTextSize(0.03);
        lat1.SetTextFont(42);
        lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
        lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
        // lat1.DrawLatex(0.255, )
        lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

        for (int i = 0; i < 3; i++)
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

// expol function with 3 rBW without mass dependent fit
#ifdef b_expol
        TF1 *BEexpol = new TF1("BEexpol", BWsum_expol_chkstar, 1.05, 2.20, 13); // expol 3
        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        // double parameters[] = {8384, f1270Mass, f1270Width, 7858, f1525Mass, f1525Width, 3218, f1710Mass, f1710Width};
        double parameters[] = {6200, f1270Mass, f1270Width, 7700, f1525Mass, f1525Width, 2500, f1710Mass, f1710Width}; // for preview presentation
        // double parameters[] = {1700, f1270Mass, f1270Width, 2300, f1525Mass, f1525Width, 500, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        // vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {2, 50*f1270WidthErr}, {4, 5 * f1525Width}, {7, 10 * f1710WidthErr}};
        vector<vector<float>> par_limits = {{2, 10 * f1270WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        // for rotational bkg with pt range 0-30 GeV/c , // 3BW/3BWAmp + expol/exponential1, 3BWAmp + expol, 3BWAmp + exponential3
        double initial_param_bkg[] = {2600, -0.02379, 8.369, -2.982}; // rotational 1-30 GeV/c (KsKs channel)
        // double initial_param_bkg[] = {2.062e5, -0.05379, 2.769, 1.1982}; // rotational 1-30 GeV/c (KsKs channel) // other file

        BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // 206 //5.845e5
        BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  // 0.04316 //-0.07378
        BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // 11.48 //2.685
        BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // -3.149 //1.176
        // BEexpol->SetParLimits(size_fitparams + 3, 0.99, 1.01);

        // // for rotational bkg with pt range 1-30 GeV/c
        // BEexpol->FixParameter(size_fitparams + 0, 5.927e5);  // 5.562e5   // Fix
        // BEexpol->SetParameter(size_fitparams + 1, -0.05466); // -0.09379  //Free
        // BEexpol->FixParameter(size_fitparams + 2, 3.26);     // 2.569     // Fix
        // BEexpol->SetParameter(size_fitparams + 3, 0.9221);   // 1.098     // Free

        // for ME bkg with pt range 1-30 GeV/c
        // Till now the problem with ME data is that, it goes negative. Now the expol fit does not go to neagtive values. Need to fix this.
        // BEexpol->SetParameter(size_fitparams + 0, -0); // 5.562e5   // Fix
        // BEexpol->SetParameter(size_fitparams + 1, 0.4);      // -0.09379  //Free
        // BEexpol->SetParameter(size_fitparams + 2, 15.5);     // 2.569     // Fix
        // BEexpol->SetParameter(size_fitparams + 3, 2.8);      // 1.098     // Free

        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(2, f1270Width);

        // BEexpol->FixParameter(4, f1525Mass);
        // BEexpol->FixParameter(5, f1525Width);

        // BEexpol->FixParameter(1, a1320Mass);
        // BEexpol->FixParameter(2, a1320Width);

        // BEexpol->FixParameter(7, f1710Mass);
        // BEexpol->FixParameter(8, f1710Width);

        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        // status codes: 4000 successful, 4 call limit, 4910 failed
        string fitstatus = "Successfull";
        if (fitResultptr->Status() != 4000)
        {
            cout << "Fit failed or call limit" << endl;
            fitstatus = "Failed";
        }
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;

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

        // TF1 *onlyBW = new TF1("onlyBW", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        // TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW = new TF1("onlyBW", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        onlyBW_clone->SetParNames("f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma");
        for (int i = 0; i < 9; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        // onlyBW->Draw("same");

        // TLegend *ltemp = new TLegend(0.20, 0.67, 0.52, 0.92);
        // ltemp->SetFillStyle(0);
        // ltemp->SetTextFont(42);
        // ltemp->SetTextSize(0.035);
        // ltemp->AddEntry(hinvMass, "Data", "lpe");
        // ltemp->AddEntry(BEexpol, "3rBW + expol", "l");
        // ltemp->AddEntry(onlyBW, "3rBW", "l");
        // ltemp->AddEntry(expol, "expol", "l");
        // ltemp->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits[3];
        for (int i = 0; i < 3; i++)
        {
            singlefits[i] = new TF1(Form("singlef%d", i), single_BW, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
            singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
            singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
            singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
            singlefits[i]->SetLineColor(colors[i]);
            singlefits[i]->SetLineStyle(2);
            singlefits[i]->SetLineWidth(2);
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
        ltemp->AddEntry(BEexpol, "3rBW + Residual BG", "l");
        ltemp->AddEntry(expol, "Residual BG", "l");
        ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
        ltemp->AddEntry(singlefits[1], "f'_{2}(1525)", "l");
        ltemp->AddEntry(singlefits[2], "f_{0}(1710)", "l");
        ltemp->Draw("same");

        TLatex lat1;
        lat1.SetNDC();
        lat1.SetTextSize(0.03);
        lat1.SetTextFont(42);
        lat1.DrawLatex(0.255, 0.89, "pp, #sqrt{s} = 13.6 TeV");
        lat1.DrawLatex(0.255, 0.85, "FT0M (0-100%), |y|<0.5");
        lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

        for (int i = 0; i < 3; i++)
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

        // // //******************************************************************************************
        // // //********************************* common for all fits ***************************************
        // // //****************************common for all fits************************************
        // // //************************** printing fit parameters in the file *******************************
        gPad->Update();
        TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
        ptstats->SetX1NDC(0.6);
        ptstats->SetX2NDC(0.99);
        ptstats->SetY1NDC(0.4);
        ptstats->SetY2NDC(0.92);
        ptstats->Draw("same");
        // hinvMass->GetXaxis()->SetRangeUser(BEexpol->GetXmin(), BEexpol->GetXmax());
        // hinvMass->SetMaximum(1.5 * hinvMass->GetMaximum());
        c->SaveAs((savepath + Form("/rBWfit_pt_%.2f_%.2f_%s.png", pT_bins[ipt], pT_bins[ipt + 1], sysvar.c_str())).c_str());

        double chi2_ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
        double fitnorm1525 = BEexpol->GetParameter(3);
        double fitnorm1525_err = BEexpol->GetParError(3);
        double fitnorm1710 = BEexpol->GetParameter(6);
        double fitnorm1710_err = BEexpol->GetParError(6);
        double fitmass1525 = BEexpol->GetParameter(4);
        double fitmass1525_err = BEexpol->GetParError(4);
        double fitmass1710 = BEexpol->GetParameter(7);
        double fitmass1710_err = BEexpol->GetParError(7);
        double fitwidth1525 = BEexpol->GetParameter(5);
        double fitwidth1525_err = BEexpol->GetParError(5);
        double fitwidth1710 = BEexpol->GetParameter(8);
        double fitwidth1710_err = BEexpol->GetParError(8);
        double fitnorm1270 = BEexpol->GetParameter(0);
        double fitnorm1270_err = BEexpol->GetParError(0);
        double fitmass1270 = BEexpol->GetParameter(1);
        double fitmass1270_err = BEexpol->GetParError(1);
        double fitwidth1270 = BEexpol->GetParameter(2);
        double fitwidth1270_err = BEexpol->GetParError(2);
        double fitrangelow = BEexpol->GetXmin();
        double fitrangehigh = BEexpol->GetXmax();
        double expol1 = BEexpol->GetParameter(9);
        double expol2 = BEexpol->GetParameter(10);
        double expol3 = BEexpol->GetParameter(11);
        double expol4 = BEexpol->GetParameter(12);

        // // file << fitstatus << endl;
        // file << std::fixed << std::setprecision(2);
        // file << fitrangelow << " - " << fitrangehigh << endl;
        // file << std::fixed << std::setprecision(1);
        // file << chi2_ndf << endl;
        // file << std::fixed << std::setprecision(0);
        // file << fitnorm1270 << " ± " << fitnorm1270_err << endl;
        // file << std::fixed << std::setprecision(1);
        // file << fitmass1270 * 1000.0 << " ± " << fitmass1270_err * 1000.0 << endl;
        // file << fitwidth1270 * 1000.0 << " ± " << fitwidth1270_err * 1000.0 << endl;
        // file << std::fixed << std::setprecision(0);
        // file << fitnorm1525 << " ± " << fitnorm1525_err << endl;
        // file << std::fixed << std::setprecision(1);
        // file << fitmass1525 * 1000.0 << " ± " << fitmass1525_err * 1000.0 << endl;
        // file << fitwidth1525 * 1000.0 << " ± " << fitwidth1525_err * 1000.0 << endl;

        // // f0(1710)
        // file << endl;
        // file << std::fixed << std::setprecision(1);
        // file << "Chi2/ndf " << endl;
        // file << chi2_ndf << endl;
        // file << std::fixed << std::setprecision(0);
        // file << fitnorm1710 << " ± " << fitnorm1710_err << endl;
        // file << std::fixed << std::setprecision(1);
        // file << fitmass1710 * 1000.0 << " ± " << fitmass1710_err * 1000.0 << endl;
        // file << fitwidth1710 * 1000.0 << " ± " << fitwidth1710_err * 1000.0 << endl;
        // file << endl;

        file << "Norm range " << kNormRangepT[0][0] << " - " << kNormRangepT[0][1] << endl;
        file << "Fit range " << fitrangelow << " - " << fitrangehigh << endl;
        file << "Significance " << significance << endl;
        file << "StatSignificance " << statSignificance << endl;
        file << "Fit parameters of f1710 " << endl;
        file << "Chi2NDF " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
        file << fitnorm1710 << " ± " << fitnorm1710_err << endl;
        file << fitmass1710 << " ± " << fitmass1710_err << endl;
        file << fitwidth1710 << " ± " << fitwidth1710_err << endl;

        // #if defined(b_exponential) || defined(b_modfied_boltzmann)
        //         file << expol1 << ", " << expol2 << ", " << expol3 << ", " << expol4 << endl;
        // #endif
        // #ifdef b_boltzman_pure
        //         file << expol1 << ", " << expol2 << ", " << expol3 << endl;
        // #endif
        //         file << endl;
        //         file << endl;
        //         file << "for analysis note " << endl;
        //         file << "$" << fitmass1270 << " \\pm " << fitmass1270_err << "$ & $" << fitwidth1270 << " \\pm " << fitwidth1270_err << "$ & " << chi2 << "/" << ndf << " & " << (int)significance[0] << endl;
        //         file << "$" << fitmass1525 << " \\pm " << fitmass1525_err << "$ & $" << fitwidth1525 << " \\pm " << fitwidth1525_err << "$ & -- & " << (int)significance[1] << endl;
        //         file << "$" << fitmass1710 << " \\pm " << fitmass1710_err << "$ & $" << fitwidth1710 << " \\pm " << fitwidth1710_err << "$ & -- & " << (int)significance[2] << endl;
        //         file << endl;
        //         file << endl;
        //         file << std::fixed << std::setprecision(4);
        //         file << "fit masses" << endl;
        //         file << fitmass1270 << ", " << fitmass1525 << ", " << fitmass1710 << endl;
        //         file << "fit masses errors" << endl;
        //         file << fitmass1270_err << ", " << fitmass1525_err << ", " << fitmass1710_err << endl;
        //         file << "fit widths" << endl;
        //         file << fitwidth1270 << ", " << fitwidth1525 << ", " << fitwidth1710 << endl;
        //         file << "fit widths errors" << endl;
        //         file << fitwidth1270_err << ", " << fitwidth1525_err << ", " << fitwidth1710_err << endl;

#ifdef residual_subtracted
        // Now subtract the residual background and plot
        TCanvas *c2 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c2, 0.14, 0.03, 0.05, 0.14);
        expol_clone->SetRange(0.99, 2.99);
        hsubtracted->Add(expol_clone, -1);
        hsubtracted->GetXaxis()->SetRangeUser(BEexpol->GetXmin(), BEexpol->GetXmax());
        hsubtracted->SetMaximum(hsubtracted->GetMaximum() * 1.5);
        hsubtracted->Draw();
        TH1F *hsubtracted_clone = (TH1F *)hsubtracted->Clone("hsubtracted_clone");
        // hsubtracted_clone->Write("3BW");

        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }
        onlyBW_clone->SetNpx(1000);
        hsubtracted->Fit("onlyBW_clone", "REBMS");
        double *obtained_parameters2 = onlyBW_clone->GetParameters();
        TLine *line = new TLine(BEexpol->GetXmin() + 0.01, 0, BEexpol->GetXmax() - 0.01, 0);
        line->SetLineColor(1);
        line->SetLineStyle(4);
        line->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits1[3];
        for (int i = 0; i < 3; i++)
        {
#if defined(b_massdepwidth) || defined(b_boltzmann_pure_massdepwidth)
            singlefits1[i] = (i < 2) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, onlyBW_clone->GetXmin(), onlyBW_clone->GetXmax(), 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, onlyBW_clone->GetXmin(), onlyBW_clone->GetXmax(), 3);
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

        TLegend *ltemp2 = new TLegend(0.20, 0.64, 0.42, 0.90);
        ltemp2->SetFillStyle(0);
        ltemp2->SetTextFont(42);
        ltemp2->SetTextSize(0.035);
        ltemp2->SetBorderSize(0);
        ltemp2->SetHeader("Residual bkg subtraction");
        ltemp2->AddEntry(hsubtracted, "Data (stat. uncert.)", "lpe");
        ltemp2->AddEntry(onlyBW_clone, "3rBW", "l");
        ltemp2->AddEntry(singlefits[0], "f_{2}(1270)", "l");
        ltemp2->AddEntry(singlefits[1], "f'_{2}(1525)", "l");
        ltemp2->AddEntry(singlefits[2], "f_{0}(1710)", "l");
        ltemp2->Draw("same");
        c2->SaveAs((savepath + "/rBWfit_residual_" + sysvar + ".png").c_str());

#endif

#ifdef doublepanelplot
        // gStyle->SetOptFit(0);
        TCanvas *c1 = new TCanvas("", "", 720, 720);
        double pad1Size, pad2Size;
        canvas_style(c1, pad1Size, pad2Size);
        c1->cd(1);
        gPad->SetTickx(1);
        hinvMass->GetXaxis()->SetTitleSize(0.04 / pad1Size);
        hinvMass->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        hinvMass->GetXaxis()->SetLabelSize(0.04 / pad1Size);
        hinvMass->GetYaxis()->SetLabelSize(0.04 / pad1Size);
        hinvMass->GetXaxis()->SetTitleOffset(1.02);
        hinvMass->GetYaxis()->SetTitleOffset(0.8);
        hinvMass->GetYaxis()->CenterTitle(0);
        hsubtracted->GetXaxis()->SetRangeUser(BEexpol->GetXmin() + 0.01, BEexpol->GetXmax() - 0.01);

        // hinvMass->SetMarkerStyle(22);
        hinvMass->SetMarkerSize(0.8);
        hinvMass->SetLineColor(1);
        hinvMass->SetMarkerColor(1);
        hinvMass->SetStats(0);
        hinvMass->SetMinimum(-100);
        // hinvMass->SetMaximum(0.9e6);
        hinvMass->SetMaximum(hinvMass->GetMaximum() * 1.1);
        hinvMass->GetYaxis()->SetMaxDigits(4);
        hinvMass->GetYaxis()->SetNdivisions(505);
        hinvMass->Draw("pe");
        expol->Draw("same");
        onlyBW->SetLineWidth(0);
        onlyBW->SetFillColor(4);
        onlyBW->SetFillStyle(1001);
        onlyBW->Draw("same");
        onlyBW->SetNpx(1000);
        gPad->Update();

        TArrow *arrow = new TArrow(0.3259053, 0.7255747, 0.3259053, 0.5158046, 0.02, "|>");
        arrow->SetFillColor(1);
        arrow->SetFillStyle(1001);
        arrow->SetLineWidth(2);
        arrow->SetNDC();
        arrow->Draw();
        TLatex *tex = new TLatex(0.3356546, 0.7816092, "f_{2}(1270)");
        tex->SetNDC();
        tex->SetTextAlign(22);
        tex->SetLineWidth(2);
        tex->Draw();
        arrow = new TArrow(0.5, 0.6048851, 0.5, 0.3951149, 0.02, "|>");
        arrow->SetFillColor(1);
        arrow->SetFillStyle(1001);
        arrow->SetLineWidth(2);
        arrow->SetNDC();
        arrow->Draw();
        tex = new TLatex(0.4916435, 0.6522989, "f'_{2}(1525)");
        tex->SetNDC();
        tex->SetTextAlign(22);
        tex->SetLineWidth(2);
        tex->Draw();
        arrow = new TArrow(0.6420613, 0.4109195, 0.6420613, 0.2011494, 0.02, "|>");
        arrow->SetFillColor(1);
        arrow->SetFillStyle(1001);
        arrow->SetLineWidth(2);
        arrow->SetNDC();
        arrow->Draw();
        tex = new TLatex(0.6448468, 0.4511494, "f_{0}(1710)");
        tex->SetNDC();
        tex->SetTextAlign(22);
        tex->SetLineWidth(2);
        tex->Draw();

        TLegend *leg = new TLegend(0.65, 0.47, 0.99, 0.77);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.06);
        leg->SetBorderSize(0);
        leg->AddEntry(hinvMass, "pp #sqrt{s} = 13.6 TeV", "lpe");
        leg->AddEntry(BEexpol, "3 BW + Residual BG", "l");
        leg->AddEntry(onlyBW, "Signal", "f");
        leg->AddEntry(expol, "Residual BG", "l");
        leg->Draw("same");

        TLatex *text4 = new TLatex(0.65, 0.80, "ALICE work in progress");
        text4->SetNDC();
        text4->SetTextSize(0.06);
        text4->SetTextFont(42);
        text4->Draw("same");

        c1->cd(2);
        gPad->SetTickx(1);
        hsubtracted->GetYaxis()->SetTitleSize(0.04 / pad2Size);
        hsubtracted->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        hsubtracted->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        hsubtracted->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        hsubtracted->GetYaxis()->SetNdivisions(504);
        hsubtracted->GetYaxis()->SetMaxDigits(4);
        hsubtracted->GetXaxis()->SetRangeUser(1.05, 2.15);
        hsubtracted->SetStats(0);
        hsubtracted->SetMinimum(0);
        // hsubtracted->SetMaximum(0.13e6);
        hsubtracted->SetMaximum(0.145e6);
        hsubtracted->SetMarkerSize(0.8);
        // hsubtracted->GetYaxis()->SetMaxDigits(10);
        hsubtracted->Draw("pe");

        TLatex lat5;
        lat5.SetNDC();
        lat5.SetTextSize(0.06);
        lat5.SetTextFont(42);
        lat5.DrawLatex(0.59, 0.80, "Residual background subtraction");
        lat5.DrawLatex(0.69, 0.70, "pp #sqrt{#it{s}} = 13.6 TeV");
        lat5.DrawLatex(0.69, 0.60, "FT0M (0-100%), |#it{y}|<0.5");
        lat5.DrawLatex(0.69, 0.50, Form("%.1f < #it{p}_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

        for (int i = 0; i < 3; i++)
        {
            singlefits1[i]->Draw("same");
        }
        c1->Update();
        // c1->SaveAs("/home/sawan/Music/r3BWfit_doublepanel.pdf");
        c1->SaveAs((savepath + Form("/rBWfit_doublepanel_pt_%.2f_%.2f.png", pT_bins[ipt], pT_bins[ipt + 1])).c_str());

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
    // plots_3BW->Close();
}
// end of main program

// *****************************************************************************************************
//************************************Fit functions************************************************* */

// We will define the single BW and the sum of 3 BWs
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

    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2];
    double yield1525 = par[3];
    double mass1525 = par[4];
    double width1525 = par[5];
    double yield1710 = par[6];
    double mass1710 = par[7];
    double width1710 = par[8];

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270);
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525);
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710);

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270);
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525);
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710);

    double fit = yield1270 * (realnum1270 * realnum1270 + imagnum1270 * imagnum1270) / (den1270 * den1270) + yield1525 * (realnum1525 * realnum1525 + imagnum1525 * imagnum1525) / (den1525 * den1525) + yield1710 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / (den1710 * den1710);

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
    double yield1525 = par[3];
    double mass1525 = par[4];
    double width1525 = par[5];
    double yield1710 = par[6];
    double mass1710 = par[7];
    double width1710 = par[8];

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    double fit = fit1270 + fit1525 + fit1710;
    return fit;
}

Double_t BWsumMassDepWidth(double *x, double *par)
{
    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2] * (TMath::Power(par[1] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double yield1525 = par[3];
    double mass1525 = par[4];
    double width1525 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double yield1710 = par[6];
    double mass1710 = par[7];
    double width1710 = par[8] * (TMath::Power(par[7] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n2);

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    double fit = fit1270 + fit1525 + fit1710;
    return fit;
}

// Now we will define the functions for the exponential background

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

Double_t expol_chkstar(double *x, double *par)
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] + x[0] * par[3]));
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
    return (BWsum(x, par) + exponential_bkg_3(x, &par[9]));
}

Double_t BWsum_expol3_hera(double *x, double *par)
{
    return (BWsum_hera(x, par) + exponential_bkg_3(x, &par[9]));
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
    return (BWsum(x, par) + Boltzmann_bkg_1(x, &par[9]));
}

Double_t BWsum_boltzmann_massdep(double *x, double *par)
{
    return (BWsumMassDepWidth(x, par) + Boltzmann_bkg_1(x, &par[9]));
}

Double_t BWsum_boltzman_2(double *x, double *par)
{
    return (BWsum(x, par) + Boltzmann_bkg_2(x, &par[9]));
}

Double_t BWsum_expol_chkstar(double *x, double *par)
{
    return (BWsum(x, par) + expol_chkstar(x, &par[9]));
}

Double_t BWsumMassDepWidth_expol3(double *x, double *par)
{
    // return (BWsumMassDepWidth(x, par) + expol_chkstar(x, &par[9]));
    return (BWsumMassDepWidth(x, par) + exponential_bkg_3(x, &par[9]));
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
