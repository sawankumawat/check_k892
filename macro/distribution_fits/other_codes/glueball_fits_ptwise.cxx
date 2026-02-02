#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);
Double_t single_BW(double *x, double *par);
Double_t BWsumMassDepWidth(double *x, double *par);
Double_t BWsumMassDepWidth_exponential(double *x, double *par);
Double_t single_BW_mass_dep_spin0(double *x, double *par);
Double_t single_BW_mass_dep_spin2(double *x, double *par);

Double_t exponential_bkg_3(double *x, double *par); // 4 parameters

Double_t single_BW_expol3(double *x, double *par);
Double_t single_BW_expol3_hera(double *x, double *par);

void glueball_fit_4rBW_ptwise()
{
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances"; // 2022 dataset (old)
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/435448/KsKs_Channel/higher-mass-resonances"; //2022 dataset (new)
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/435450/KsKs_Channel/higher-mass-resonances"; // 2023 dataset/
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/435449/KsKs_Channel/higher-mass-resonances"; // 2024 dataset

    string savepath = path + "/fits/4rBw_fits/fit_validation/pt_dependent";

    gSystem->Exec(("mkdir -p " + savepath).c_str());

    // TFile *f = new TFile((path + "/hglue_ROTATED_norm_2.50_2.60_pt_3.00_30.00.root").c_str(), "READ"); /default
    TFile *f = new TFile((path + "/hglue_ROTATED.root").c_str(), "READ");
    // TFile *f = new TFile((path + "/hglue_MIX_temp.root").c_str(), "READ");
    // TFile *f = new TFile((path + "/hglue_ROTATED_pT_0_1.root").c_str(), "READ");
    // TFile *f = new TFile((path + "/hglue_ROTATED_pT_2.0_12.0.root").c_str(), "READ");
    // TFile *f = new TFile((path + "/hglue_ROTATED_allPt.root").c_str(), "READ");

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

    TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_exponential, fitlow, fithigh, 16); // expol 3
    TF1 *BEexpol_initial = new TF1("BEexpol_initial", BWsumMassDepWidth_exponential, fitlow, fithigh, 16);
    TF1 *BEexpol_reduced = new TF1("BEexpol_reduced", BWsumMassDepWidth_exponential, fitlow, fithigh, 16);
    string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
    for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
    {
        BEexpol->SetParName(i, parnames[i].c_str());
        BEexpol_reduced->SetParName(i, parnames[i].c_str());
    }

    // double parameters[] = {3500, f1270Mass, f1270Width, 2000, a1320Mass, a1320Width, 7000, f1525Mass, f1525Width, 2200, f1710Mass, f1710Width}; // rebin twice (2022 and 2024 dataset)
    double parameters[] = {2.1e4, f1270Mass, f1270Width, 1.8e4, a1320Mass, a1320Width, 4.3e4, f1525Mass, f1525Width, 1.35e4, f1710Mass, f1710Width}; // LHC23_pass4_thin
    // double parameters[] = {1500, f1270Mass, f1270Width, 1000, a1320Mass, a1320Width, 3700, f1525Mass, f1525Width, 1300, f1710Mass, f1710Width}; // default
    // double parameters[] = {400, f1270Mass, f1270Width, 370, a1320Mass, a1320Width, 1200, f1525Mass, f1525Width, 450, f1710Mass, f1710Width}; // pt range 3-5 GeV/c, 5-8 GeV/c
    // double parameters[] = {700, f1270Mass, f1270Width, 706, a1320Mass, a1320Width, 2200, f1525Mass, f1525Width, 1000, f1710Mass, f1710Width}; // pt range 3-5 GeV/c, 5-8 GeV/c rebin 2
    // double parameters[] = {35, f1270Mass, f1270Width, 66, a1320Mass, a1320Width, 166, f1525Mass, f1525Width, 280, f1710Mass, f1710Width}; // pt range 8-15 GeV/c rebin 3
    // double parameters[] = {690, f1270Mass, f1270Width, 714, a1320Mass, a1320Width, 2300, f1525Mass, f1525Width, 500, f1710Mass, f1710Width}; // pt range 2-3 GeV/c rebin 2 (fit 1.05-2.25)
    // double parameters[] = {1470, f1270Mass, f1270Width, 714, a1320Mass, a1320Width, 1300, f1525Mass, f1525Width, 250, f1710Mass, f1710Width}; // pt range 1-2 GeV/c rebin 2 (fit 1.05-2.20), keep fit options as REBS (remove M)
    // double parameters[] = {2000, f1270Mass, f1270Width, 1714, a1320Mass, a1320Width, 5500, f1525Mass, f1525Width, 3000, f1710Mass, f1710Width}; // ME 3-30 GeV/c

    int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

    for (int i = 0; i < size_fitparams; i++)
    {
        BEexpol->SetParameter(i, parameters[i]);
        BEexpol_reduced->SetParameter(i, parameters[i]);
        BEexpol_initial->SetParameter(i, parameters[i]);
    }

    // //********systematic studies*************
    // double initial_param_bkg[] = {7.37518e5, 0.0134, 3.071167, 1.04}; // rebin twice (2022 dataset)
    double initial_param_bkg[] = {3.6e6, -0.04, 2.82, 1.05}; // rebin twice (2023 dataset)
    // double initial_param_bkg[] = {2.6e6, 0.20, 3.67, 0.8}; // rebin twice (2024 dataset)
    // double initial_param_bkg[] = {2.37518e5, -0.102044, 3.071167, 1.354864}; // 0-30 GeV/c, 3sigma
    // double initial_param_bkg[] = {4.618e6, 0.00774, 2.82, 1.03}; // pass_4_thin (old)
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
    BEexpol_initial->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Free
    BEexpol_initial->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Fix for medium train
    BEexpol_initial->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Free
    BEexpol_initial->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

    BEexpol_initial->FixParameter(2, f1270Width);
    BEexpol_initial->FixParameter(5, a1320Width);
    BEexpol_initial->FixParameter(8, f1525Width);

    BEexpol_initial->FixParameter(1, f1270Mass);
    BEexpol_initial->FixParameter(4, a1320Mass);
    BEexpol_initial->FixParameter(7, f1525Mass);

    BEexpol_initial->FixParameter(10, f1710Mass);
    BEexpol_initial->FixParameter(11, f1710Width);
    // TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "RELMS0");
    hinvMass->Fit("BEexpol_initial", "REBS0"); // comment while using toy mc and likelihood fits
    // TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBS"); // comment while using toy mc and likelihood fits

    //=============================================================//
    // Again setting parameters for the next iteration in the fit.
    //=============================================================//
    for (int iparams = 0; iparams < 16; iparams++)
    {
        BEexpol->SetParameter(iparams, BEexpol_initial->GetParameter(iparams));
    }

    vector<vector<double>> par_limits = {{1, 2 * f1270Width}, {2, 3 * f1270WidthErr}, {4, 2 * a1320Width}, {5, 5 * a1320WidthErr}, {7, 5 * f1525Width}, {8, 5 * f1525WidthErr}, {10, 1 * f1710Width}, {11, 5 * f1710WidthErr}};
    int limits_size = par_limits.size();
    for (int i = 0; i < limits_size; i++)
    {
        int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
        BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        BEexpol_reduced->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
    }

    BEexpol->SetParLimits(0, 0, 1e6);
    BEexpol->SetParLimits(3, 0, 1e6);
    // BEexpol->SetParLimits(6, 0, 1e6);

    // BEexpol->FixParameter(2, f1270Width);
    // BEexpol->FixParameter(5, a1320Width);
    // BEexpol->FixParameter(8, f1525Width);
    BEexpol->SetParameter(2, f1270Width);
    BEexpol->SetParameter(5, a1320Width);
    BEexpol->SetParameter(8, f1525Width);

    BEexpol->SetParameter(1, f1270Mass);
    BEexpol->SetParameter(4, a1320Mass);
    BEexpol->SetParameter(7, f1525Mass);

    BEexpol->SetParameter(10, f1710Mass);
    BEexpol->SetParameter(11, f1710Width);

    TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS"); // comment while using toy mc and likelihood fits
    double *obtained_parameters = BEexpol->GetParameters();         // comment while using toy mc and likelihood fits

    for (int ipt = 0; ipt < Npt; ipt++)
    {
        float lowpT = pT_bins[ipt];
        float highpT = pT_bins[ipt + 1];
        // float lowpT = 0.0;
        // float highpT = 30.0;
        ofstream file;
        file.open((savepath + Form("/fit_params_pT_%.1f-%.1f", lowpT, highpT) + sysvar + ".txt").c_str());

        vector<vector<float>> fitranges = {
            {1.05, 2.20}, // default
            {1.05, 2.25},
            {1.05, 2.15},
            {1.08, 2.20},
            {1.02, 2.20},
            {1.02, 2.25},
            {1.08, 2.15}};

        // Vectors to store fit parameters for all fit ranges
        vector<string> fit_range_labels;
        vector<string> chi2ndf_values;
        vector<string> norm1525_values;
        vector<string> mass1525_values;
        vector<string> norm1710_values;
        vector<string> mass1710_values;
        vector<string> width1710_values;

        float fitlow = fitranges[irange][0];
        float fithigh = fitranges[irange][1];

        cout << "Low fit range is " << fitlow << ", High fit range is " << fithigh << endl;

        TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", lowpT, highpT));
        TH1F *hraw = (TH1F *)f->Get(Form("ksks_invmass_pt_%.1f_%.1f", lowpT, highpT));
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
        // TCanvas *c_reducedFit = new TCanvas("c_reducedFit", "Reduced Fit", 720, 720);
        // SetCanvasStyle(c_reducedFit, 0.14, 0.03, 0.05, 0.13);
        TCanvas *c = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
        c->cd();
        // hinvMass->Rebin(2);
        double binwidthfile = (hinvMass->GetXaxis()->GetXmax() - hinvMass->GetXaxis()->GetXmin()) / hinvMass->GetXaxis()->GetNbins();
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

        // onlyBW_clone->FixParameter(10, f1710Mass);
        onlyBW_clone->FixParameter(11, f1710Width);

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
        lat1.DrawLatex(0.255, 0.815, Form("%.1f < p_{T} < %.1f GeV/#it{c}", lowpT, highpT));

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
        hsubtracted->Fit("onlyBW_clone", "REBMS");
        double *obtained_parameters2 = onlyBW_clone->GetParameters();
        TLine *line = new TLine(BEexpol->GetXmin() + 0.01, 0, BEexpol->GetXmax() - 0.01, 0);
        line->SetLineColor(1);
        line->SetLineStyle(4);
        line->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits1[4];
        for (int i = 0; i < 4; i++)
        {
            singlefits1[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, 1.00, 2.5, 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, 1.00, 2.5, 3);

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
        c2->SaveAs((savepath + Form("/rBWfit_residual_pt_%.1f_%.1f.png", lowpT, highpT)).c_str());

        //======================Yield calculation from function integration============================

        float ptBinWidth = highpT - lowpT;
        double fitrangelow = 1.001;
        double fitrangehigh = 2.499;
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

        float nsigma_yield = 3.0;
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
        if (fitRegion1270_high > fitrangehigh)
            fitRegion1270_high = fitrangehigh;
        if (fitRegion1320_low < fitrangelow)
            fitRegion1320_low = fitrangelow;
        if (fitRegion1320_high > fitrangehigh)
            fitRegion1320_high = fitrangehigh;
        if (fitRegion1525_low < fitrangelow)
            fitRegion1525_low = fitrangelow;
        if (fitRegion1525_high > fitrangehigh)
            fitRegion1525_high = fitrangehigh;
        if (fitRegion1710_low < fitrangelow)
            fitRegion1710_low = fitrangelow;
        if (fitRegion1710_high > fitrangehigh)
            fitRegion1710_high = fitrangehigh;

        // // Yield calculation
        double yield1270 = singlefits[0]->Integral(fitRegion1270_low, fitRegion1270_high) / (ptBinWidth * binwidthfile * total_events);
        double yield1320 = singlefits[1]->Integral(fitRegion1320_low, fitRegion1320_high) / (ptBinWidth * binwidthfile * total_events);
        double yield1525 = singlefits[2]->Integral(fitRegion1525_low, fitRegion1525_high) / (ptBinWidth * binwidthfile * total_events);
        double yield1710 = singlefits[3]->Integral(fitRegion1710_low, fitRegion1710_high) / (ptBinWidth * binwidthfile * total_events);
        // double yieldbkg1 = expol_clone->Integral(2.20, 2.30) / (ptBinWidth * binwidthfile * total_events);
        // double yieldbkg2 = expol_clone->Integral(2.30, 2.40) / (ptBinWidth * binwidthfile * total_events);
        // double yieldbkg3 = expol_clone->Integral(2.40, 2.50) / (ptBinWidth * binwidthfile * total_events);

        cout << "Total events: " << total_events << endl;
        // cout << "Cos theta bin width: " << ptBinWidth << endl;
        // cout << "Energy bin width: " << binwidthfile << endl;
        // cout << "Yield region for 1710 : " << f1710Mass - nsigma_yield * f1710Width << " to " << f1710Mass + nsigma_yield * f1710Width << endl;
        cout << "Area under the curve 1710: " << singlefits[3]->Integral(fitRegion1710_low, fitRegion1710_high) << endl;
        cout << "Amplitude for 1710 is " << obtained_parameters[9] << endl;
        cout << "Yield for 1710 is " << yield1710 << endl;
        cout << endl;
        cout << "Area under the curve for f1525: " << singlefits[2]->Integral(fitRegion1525_low, fitRegion1525_high) << endl;
        cout << "Amplitude for f1525 is " << obtained_parameters[6] << endl;
        cout << "Yield for f1525 is " << yield1525 << endl;
        cout << endl;

        // cout << "Yield in bkg region 1: " << yieldbkg1 << endl;
        // cout << "Yield in bkg region 2: " << yieldbkg2 << endl;
        // cout << "Yield in bkg region 3: " << yieldbkg3 << endl;

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

        // cout << "Yield 1270: " << yield1270 << " +- " << yield1270_err << endl;
        // cout << "Yield 1320: " << yield1320 << " +- " << yield1320_err << endl;
        // cout << "Yield 1525: " << yield1525 << " +- " << yield1525_err << endl;
        // cout << "Yield 1710: " << yield1710 << " +- " << yield1710_err << endl;

        //======================Yield calculation bin counting============================
        // Note: The reason why the yield of f0(1710) from bin counting is slightly different than function integration is because the fitting of f0(1710) is not good, so the area under the fit curve is not reliable.
        double hBCError1525, hBCError1710, hBCError1270, hBCError1320;
        int bmin1525 = hinvMass->GetXaxis()->FindBin(fitRegion1525_low);
        int bmax1525 = hinvMass->GetXaxis()->FindBin(fitRegion1525_high);
        int bmin1710 = hinvMass->GetXaxis()->FindBin(fitRegion1710_low);
        int bmax1710 = hinvMass->GetXaxis()->FindBin(fitRegion1710_high);
        int bmin1270 = hinvMass->GetXaxis()->FindBin(fitRegion1270_low);
        int bmax1270 = hinvMass->GetXaxis()->FindBin(fitRegion1270_high);
        int bmin1320 = hinvMass->GetXaxis()->FindBin(fitRegion1320_low);
        int bmax1320 = hinvMass->GetXaxis()->FindBin(fitRegion1320_high);
        double Yield_bincount_hist1525 = hinvMass->IntegralAndError(bmin1525, bmax1525, hBCError1525);
        double Yield_bincount_hist1710 = hinvMass->IntegralAndError(bmin1710, bmax1710, hBCError1710);
        double Yield_bincount_hist1270 = hinvMass->IntegralAndError(bmin1270, bmax1270, hBCError1270);
        double Yield_bincount_hist1320 = hinvMass->IntegralAndError(bmin1320, bmax1320, hBCError1320);
        double bkgvalue1525 = expol->Integral(hinvMass->GetBinLowEdge(bmin1525), hinvMass->GetBinLowEdge(bmax1525 + 1));
        double bkgvalue1710 = expol->Integral(hinvMass->GetBinLowEdge(bmin1710), hinvMass->GetBinLowEdge(bmax1710 + 1));
        double bkgvalue1270 = expol->Integral(hinvMass->GetBinLowEdge(bmin1270), hinvMass->GetBinLowEdge(bmax1270 + 1));
        double bkgvalue1320 = expol->Integral(hinvMass->GetBinLowEdge(bmin1320), hinvMass->GetBinLowEdge(bmax1320 + 1));
        double fYield_BinCount1525 = Yield_bincount_hist1525 - (bkgvalue1525 / binwidthfile);
        double fYield_BinCount1710 = Yield_bincount_hist1710 - (bkgvalue1710 / binwidthfile);
        double fYield_BinCount1270 = Yield_bincount_hist1270 - (bkgvalue1270 / binwidthfile);
        double fYield_BinCount1320 = Yield_bincount_hist1320 - (bkgvalue1320 / binwidthfile);

        // Compute resonance cross-feed in the counting window using the fitted line shapes
        const double window1525_low = hinvMass->GetBinLowEdge(bmin1525);
        const double window1525_high = hinvMass->GetBinLowEdge(bmax1525 + 1);
        const double window1710_low = hinvMass->GetBinLowEdge(bmin1710);
        const double window1710_high = hinvMass->GetBinLowEdge(bmax1710 + 1);
        const double window1270_low = hinvMass->GetBinLowEdge(bmin1270);
        const double window1270_high = hinvMass->GetBinLowEdge(bmax1270 + 1);
        const double window1320_low = hinvMass->GetBinLowEdge(bmin1320);
        const double window1320_high = hinvMass->GetBinLowEdge(bmax1320 + 1);

        double overlap_other1525 = 0.0; // This overlap is calculated using function integration. So bin counting is not suitable for multiple resonances.
        overlap_other1525 += singlefits[0]->Integral(window1525_low, window1525_high) / binwidthfile;
        overlap_other1525 += singlefits[1]->Integral(window1525_low, window1525_high) / binwidthfile;
        overlap_other1525 += singlefits[3]->Integral(window1525_low, window1525_high) / binwidthfile;

        double overlap_other1710 = 0.0;
        overlap_other1710 += singlefits[0]->Integral(window1710_low, window1710_high) / binwidthfile;
        overlap_other1710 += singlefits[1]->Integral(window1710_low, window1710_high) / binwidthfile;
        overlap_other1710 += singlefits[2]->Integral(window1710_low, window1710_high) / binwidthfile;

        double overlap_other1270 = 0.0;
        overlap_other1270 += singlefits[0]->Integral(window1270_low, window1270_high) / binwidthfile;
        overlap_other1270 += singlefits[1]->Integral(window1270_low, window1270_high) / binwidthfile;
        overlap_other1270 += singlefits[3]->Integral(window1270_low, window1270_high) / binwidthfile;

        double overlap_other1320 = 0.0;
        overlap_other1320 += singlefits[0]->Integral(window1320_low, window1320_high) / binwidthfile;
        overlap_other1320 += singlefits[1]->Integral(window1320_low, window1320_high) / binwidthfile;
        overlap_other1320 += singlefits[3]->Integral(window1320_low, window1320_high) / binwidthfile;

        double fYield_BinCount1525_clean = std::max(0.0, fYield_BinCount1525 - overlap_other1525);
        double fYield_BinCount1710_clean = std::max(0.0, fYield_BinCount1710 - overlap_other1710);
        double fYield_BinCount1270_clean = std::max(0.0, fYield_BinCount1270 - overlap_other1270);
        double fYield_BinCount1320_clean = std::max(0.0, fYield_BinCount1320 - overlap_other1320);

        double sum_tail_correction1525 = (singlefits[2]->Integral(1.00, hinvMass->GetBinLowEdge(bmin1525)) + singlefits[2]->Integral(hinvMass->GetBinLowEdge(bmax1525 + 1), 2.5)) / binwidthfile;
        double sum_tail_correction1710 = (singlefits[3]->Integral(1.00, hinvMass->GetBinLowEdge(bmin1710)) + singlefits[3]->Integral(hinvMass->GetBinLowEdge(bmax1710 + 1), 2.5)) / binwidthfile;
        double sum_tail_correction1270 = (singlefits[0]->Integral(1.00, hinvMass->GetBinLowEdge(bmin1270)) + singlefits[0]->Integral(hinvMass->GetBinLowEdge(bmax1270 + 1), 2.5)) / binwidthfile;
        double sum_tail_correction1320 = (singlefits[1]->Integral(1.00, hinvMass->GetBinLowEdge(bmin1320)) + singlefits[1]->Integral(hinvMass->GetBinLowEdge(bmax1320 + 1), 2.5)) / binwidthfile;
        double Total_Ybincounting1525 = (sum_tail_correction1525 + fYield_BinCount1525_clean) / (total_events * ptBinWidth);
        double Total_Ybincounting1710 = (sum_tail_correction1710 + fYield_BinCount1710_clean) / (total_events * ptBinWidth);
        double Total_Ybincounting1270 = (sum_tail_correction1270 + fYield_BinCount1270_clean) / (total_events * ptBinWidth);
        double Total_Ybincounting1320 = (sum_tail_correction1320 + fYield_BinCount1320_clean) / (total_events * ptBinWidth);

        // Error calculation
        TF1 *fitFcn2_plusm1525 = new TF1("fitFcn2_plusm1525", single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
        TF1 *fitFcn2_minusm1525 = new TF1("fitFcn2_minusm1525", single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
        TF1 *fitFcn2_plusm1710 = new TF1("fitFcn2_plusm1710", single_BW_mass_dep_spin0, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
        TF1 *fitFcn2_minusm1710 = new TF1("fitFcn2_minusm1710", single_BW_mass_dep_spin0, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
        TF1 *fitFcn2_plusm1270 = new TF1("fitFcn2_plusm1270", single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
        TF1 *fitFcn2_minusm1270 = new TF1("fitFcn2_minusm1270", single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
        TF1 *fitFcn2_plusm1320 = new TF1("fitFcn2_plusm1320", single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
        TF1 *fitFcn2_minusm1320 = new TF1("fitFcn2_minusm1320", single_BW_mass_dep_spin2, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);

        fitFcn2_plusm1525->FixParameter(0, singlefits[2]->GetParameter(0));
        fitFcn2_plusm1525->FixParameter(1, singlefits[2]->GetParameter(1) + singlefits[2]->GetParError(1));
        fitFcn2_plusm1525->FixParameter(2, f1525Width);

        fitFcn2_minusm1525->FixParameter(0, singlefits[2]->GetParameter(0));
        fitFcn2_minusm1525->FixParameter(1, singlefits[2]->GetParameter(1) - singlefits[2]->GetParError(1));
        fitFcn2_minusm1525->FixParameter(2, f1525Width);

        fitFcn2_plusm1710->FixParameter(0, singlefits[3]->GetParameter(0));
        fitFcn2_plusm1710->FixParameter(1, singlefits[3]->GetParameter(1) + singlefits[3]->GetParError(1));
        fitFcn2_plusm1710->FixParameter(2, f1710Width);
        fitFcn2_minusm1710->FixParameter(0, singlefits[3]->GetParameter(0));
        fitFcn2_minusm1710->FixParameter(1, singlefits[3]->GetParameter(1) - singlefits[3]->GetParError(1));
        fitFcn2_minusm1710->FixParameter(2, f1710Width);

        fitFcn2_plusm1270->FixParameter(0, singlefits[0]->GetParameter(0));
        fitFcn2_plusm1270->FixParameter(1, singlefits[0]->GetParameter(1) + singlefits[0]->GetParError(1));
        fitFcn2_plusm1270->FixParameter(2, f1270Width);
        fitFcn2_minusm1270->FixParameter(0, singlefits[0]->GetParameter(0));
        fitFcn2_minusm1270->FixParameter(1, singlefits[0]->GetParameter(1) - singlefits[0]->GetParError(1));
        fitFcn2_minusm1270->FixParameter(2, f1270Width);

        fitFcn2_plusm1320->FixParameter(0, singlefits[1]->GetParameter(0));
        fitFcn2_plusm1320->FixParameter(1, singlefits[1]->GetParameter(1) + singlefits[1]->GetParError(1));
        fitFcn2_plusm1320->FixParameter(2, a1320Width);
        fitFcn2_minusm1320->FixParameter(0, singlefits[1]->GetParameter(0));
        fitFcn2_minusm1320->FixParameter(1, singlefits[1]->GetParameter(1) - singlefits[1]->GetParError(1));
        fitFcn2_minusm1320->FixParameter(2, a1320Width);

        double Tail_correction_plusm1525 = (fitFcn2_plusm1525->Integral(1.00, hinvMass->GetBinLowEdge(bmin1525)) + (fitFcn2_plusm1525->Integral(hinvMass->GetBinLowEdge(bmax1525 + 1), 2.5))) / binwidthfile;
        double Tail_correction_minusm1525 = ((fitFcn2_minusm1525->Integral(1.00, hinvMass->GetBinLowEdge(bmin1525)) + fitFcn2_minusm1525->Integral(hinvMass->GetBinLowEdge(bmax1525 + 1), 2.5))) / binwidthfile;
        double Tail_correction_plusm1710 = (fitFcn2_plusm1710->Integral(1.00, hinvMass->GetBinLowEdge(bmin1710)) + (fitFcn2_plusm1710->Integral(hinvMass->GetBinLowEdge(bmax1710 + 1), 2.5))) / binwidthfile;
        double Tail_correction_minusm1710 = ((fitFcn2_minusm1710->Integral(1.00, hinvMass->GetBinLowEdge(bmin1710)) + fitFcn2_minusm1710->Integral(hinvMass->GetBinLowEdge(bmax1710 + 1), 2.5))) / binwidthfile;
        double Tail_correction_plusm1270 = (fitFcn2_plusm1270->Integral(1.00, hinvMass->GetBinLowEdge(bmin1270)) + (fitFcn2_plusm1270->Integral(hinvMass->GetBinLowEdge(bmax1270 + 1), 2.5))) / binwidthfile;
        double Tail_correction_minusm1270 = ((fitFcn2_minusm1270->Integral(1.00, hinvMass->GetBinLowEdge(bmin1270)) + fitFcn2_minusm1270->Integral(hinvMass->GetBinLowEdge(bmax1270 + 1), 2.5))) / binwidthfile;
        double Tail_correction_plusm1320 = (fitFcn2_plusm1320->Integral(1.00, hinvMass->GetBinLowEdge(bmin1320)) + (fitFcn2_plusm1320->Integral(hinvMass->GetBinLowEdge(bmax1320 + 1), 2.5))) / binwidthfile;
        double Tail_correction_minusm1320 = ((fitFcn2_minusm1320->Integral(1.00, hinvMass->GetBinLowEdge(bmin1320)) + fitFcn2_minusm1320->Integral(hinvMass->GetBinLowEdge(bmax1320 + 1), 2.5))) / binwidthfile;
        double Error1525 = sum_tail_correction1525 - Tail_correction_plusm1525;
        double Error1710 = sum_tail_correction1710 - Tail_correction_plusm1710;
        double Error1270 = sum_tail_correction1270 - Tail_correction_plusm1270;
        double Error1320 = sum_tail_correction1320 - Tail_correction_plusm1320;
        double Final_pro_error1525 = TMath::Sqrt(Error1525 * Error1525 + hBCError1525 * hBCError1525) / (total_events * ptBinWidth);
        double Final_pro_error1710 = TMath::Sqrt(Error1710 * Error1710 + hBCError1710 * hBCError1710) / (total_events * ptBinWidth);
        double Final_pro_error1270 = TMath::Sqrt(Error1270 * Error1270 + hBCError1270 * hBCError1270) / (total_events * ptBinWidth);
        double Final_pro_error1320 = TMath::Sqrt(Error1320 * Error1320 + hBCError1320 * hBCError1320) / (total_events * ptBinWidth);

        cout << "f1525 yield from intergration is " << yield1525 << " +- " << yield1525_err << endl;
        cout << "f1525 yield from bincounting is " << Total_Ybincounting1525 << " +- " << Final_pro_error1525 << endl;
        cout << endl;
        cout << "f01710 yield from intergration is " << yield1710 << " +- " << yield1710_err << endl;
        cout << "f1710 yield from bincounting is " << Total_Ybincounting1710 << " +- " << Final_pro_error1710 << endl;
        cout << endl;
        cout << "f1270 yield from intergration is " << yield1270 << " +- " << yield1270_err << endl;
        cout << "f1270 yield from bincounting is " << Total_Ybincounting1270 << " +- " << Final_pro_error1270 << endl;
        cout << endl;
        cout << "a1320 yield from intergration is " << yield1320 << " +- " << yield1320_err << endl;
        cout << "a1320 yield from bincounting is " << Total_Ybincounting1320 << " +- " << Final_pro_error1320 << endl;
        cout << endl;
        // cout<<"Fitting range low is "<<BEexpol->GetXmin()<<" and high is "<<BEexpol->GetXmax()<<endl;

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
        // file << yield1710 << " ± " << yield1710_err << endl;
        // file << Total_Ybincounting1525 << " ± " << Final_pro_error1710 << endl;
        // file << fitmass1710 << " ± " << fitmass1710_err << endl;
        // file << fitwidth1710 << " ± " << fitwidth1710_err << endl;
        // file << "for f1525" << endl;
        // file << yield1525 << " ± " << yield1525_err << endl;
        // file << Total_Ybincounting1525 << " ± " << Final_pro_error1525 << endl;
        // file << fitmass1525 << " ± " << fitmass1525_err << endl;
        // file << fitwidth1525 << " ± " << fitwidth1525_err << endl;
        // file << "for f1270" << endl;
        // file << yield1270 << " ± " << yield1270_err << endl;
        // file << Total_Ybincounting1270 << " ± " << Final_pro_error1270 << endl;
        // file << fitmass1270 << " ± " << fitmass1270_err << endl;
        // file << fitwidth1270 << " ± " << fitwidth1270_err << endl;
        // file << "for a1320" << endl;
        // file << yield1320 << " ± " << yield1320_err << endl;
        // file << Total_Ybincounting1320 << " ± " << Final_pro_error1320 << endl;
        // file << fitmass1320 << " ± " << fitmass1320_err << endl;
        // file << fitwidth1320 << " ± " << fitwidth1320_err << endl;
        // // file << "Background" << endl;
        // // file << yieldbkg1 << " ± " << expol_clone->GetParError(0) << endl;
        // // file << yieldbkg2 << " ± " << expol_clone->GetParError(1) << endl;
        // // file << yieldbkg3 << " ± " << expol_clone->GetParError(2) << endl;
    }
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

Double_t single_BW(double *x, double *par)
{
    double yield = par[0];
    double mass = par[1];
    double width = par[2];

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));
    return fit;
}

Double_t exponential_bkg_3(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[3])));
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

Double_t BWsumMassDepWidth_exponential(double *x, double *par)
{
    // return (BWsumMassDepWidth(x, par) + expol_chkstar(x, &par[12]));
    return (BWsumMassDepWidth(x, par) + exponential_bkg_3(x, &par[12]));
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