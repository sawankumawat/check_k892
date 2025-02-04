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
Double_t BWsumMassDepWidth_expol1(double *x, double *par);

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
    ofstream file;
    file.open("fit_params.txt");

    string savepath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/fits/4rBw_fits";

    TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60_fullpt.root", "READ"); // full pT range

    // TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60_pt1.root", "READ"); // pT 1 to 30 GeV/c

    TFile *plots_4BW = new TFile("root_files/4rBW_plots_expol.root", "RECREATE");

    int colors[] = {4, 6, 28, 46};
    double masses[] = {f1270Mass, a1320Mass, f1525Mass, f1710Mass};
    double widths[] = {f1270Width, a1320Width, f1525Width, f1710Width};
    string resonance_names[] = {"f_{2}(1270)", "a_{2}(1320)^{0}", "f'_{2}(1525)", "f_{0}(1710)"};
    double purity, significance, chi2ndf;

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
// #define b_boltzman
// #define b_modified_Boltzmann
    #define b_massdepwidth
    // #define b_expol1

    // #define residual_subtracted
    // #define doublepanelplot

    for (int ipt = 0; ipt < Npt; ipt++)
    {
        // TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", 1.0, 30.0));
        TH1F *hinvMass = (TH1F *)f->Get(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", pT_bins[ipt], pT_bins[ipt + 1]));

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
        hinvMass->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
        hinvMass->GetYaxis()->SetTitle(Form("Counts/%.2f GeV/#it{c}^{2}", binwidthfile));
        hinvMass->Draw();
        TH1F *hsubtracted = (TH1F *)hinvMass->Clone("hsubtracted");
        TH1F *hsubtracted_res = (TH1F *)hinvMass->Clone("hsubtracted_res");
        // gStyle->SetOptStat(1110);
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1111);
        vector<tuple<float, int, float, float>> fit_parameters;

// // // //************************************************************************ */
// // // // **************** For BW sum with expol HERA ****************************

// // // Default fitting range is 1.02 to 2.20. Four types of fitting range variations: extend left (1.0), extend right (2.50), large range (1.0 to 2.50), small range (1.05 to 2.15)
#ifdef b_modified_Boltzmann

        TF1 *BEexpol = new TF1("BEexpol", BWsum_expol3, 1.05, 2.25, 16); // expol 3
        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        // double parameters[] = {8384, f1270Mass, f1270Width, 8000, a1320Mass, a1320Width, 7858, f1525Mass, f1525Width, 3218, f1710Mass, f1710Width};
        double parameters[] = {6000, f1270Mass, f1270Width, 6000, a1320Mass, a1320Width, 8000, f1525Mass, f1525Width, 6000, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {7, 5 * f1525Width}, {10, 20 * f1710WidthErr}};
        int limits_size = par_limits.size();
        // for (int i = 0; i < limits_size; i++)
        // {
        //     int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
        //     BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        // }

        double initial_param_bkg[] = {7.791e5, 0.004, 2.902, 0.9726}; // rotational 0-30 GeV/c (KsKs channel)

        // Initial parameters for background
        BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Fix
        BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Free
        BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Fix
        BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

        // // for rotational bkg with pt range 1-30 GeV/c
        // BEexpol->FixParameter(size_fitparams + 0, 5.927e5);  // 5.562e5   // Fix
        // BEexpol->SetParameter(size_fitparams + 1, -0.05466); // -0.09379  //Free
        // BEexpol->FixParameter(size_fitparams + 2, 3.26);     // 2.569     // Fix
        // BEexpol->SetParameter(size_fitparams + 3, 0.9221);   // 1.098     // Free

        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(4, a1320Mass);
        // BEexpol->FixParameter(7, f1525Mass);

        // BEexpol->FixParameter(2, f1270Width);
        // BEexpol->FixParameter(5, a1320Width);
        // BEexpol->FixParameter(8, f1525Width);

        // BEexpol->FixParameter(10, f1710Mass);
        // BEexpol->FixParameter(11, f1710Width);

        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
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
            int binlow = hinvMass->GetXaxis()->FindBin(masses[i] - 3 * widths[i]);
            int binhigh = hinvMass->GetXaxis()->FindBin(masses[i] + 3 * widths[i]);
            double significance_den = TMath::Sqrt(hinvMass->Integral(binlow, binhigh));
            significance = significance_num / significance_den;
            cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
            purity = significance_num * 100 / hinvMass->Integral(binlow, binhigh);
            cout << "Purity of " << resonance_names[i] << " is " << purity << endl;
        }
#endif

        // ************************** fit with mass depndent width BW ***************************************

#ifdef b_massdepwidth

        // TF1 *BEexpol = new TF1("BEexpol", BWsum_expol3, 1.05, 2.20, 16); // expol 3
        // TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_expol1, 1.05, 2.20, 16); // expol 3

        TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_exponential, 1.05, 2.25, 16); // expol 3

        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        // double parameters[] = {8384, f1270Mass, f1270Width, 8000, a1320Mass, a1320Width, 7858, f1525Mass, f1525Width, 3218, f1710Mass, f1710Width};
        double parameters[] = {6000, f1270Mass, f1270Width, 6000, a1320Mass, a1320Width, 8000, f1525Mass, f1525Width, 6000, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        // vector<vector<float>> par_limits = {{1, 5 * f1270Width}, {2, 10 * f1270WidthErr}, {7, 5 * f1525Width}, {10, 20 * f1710WidthErr}};
        vector<vector<float>> par_limits = {{2, 8 * f1270WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        double initial_param_bkg[] = {7.25e5, 0.0242, 3.656, 0.982}; // rotational 1-30 GeV/c (KsKs channel)
        // double initial_param_bkg[] = {7.791e5, 0.004, 2.902, 0.9726}; // rotational 0-30 GeV/c (KsKs channel)

        // Initial parameters for background
        BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Fix
        BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Free
        BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Fix
        BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

        BEexpol->FixParameter(1, f1270Mass);
        BEexpol->FixParameter(4, a1320Mass);
        BEexpol->FixParameter(7, f1525Mass);

        BEexpol->FixParameter(2, f1270Width);
        BEexpol->FixParameter(5, a1320Width);
        BEexpol->FixParameter(8, f1525Width);

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
        // onlyBW_clone->FixParameter(1, f1270Mass);
        // onlyBW_clone->FixParameter(2, f1270Width);
        // onlyBW_clone->FixParameter(4, a1320Mass);
        // onlyBW_clone->FixParameter(5, a1320Width);
        // onlyBW_clone->FixParameter(7, f1525Mass);
        // onlyBW_clone->FixParameter(8, f1525Width);
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            onlyBW->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        // onlyBW->Draw("same");

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

        // // // //************************************************************************ */
        // // // // **************** For BW sum with exp + pol1 as used in Charged kstar **************************
#ifdef b_expol1
        TF1 *BEexpol = new TF1("BEexpol", BWsum_expol_chkstar, 1.05, 2.20, 16); // expol 3
        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "A", "n", "b", "c"};
        // for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        // {
        //     BEexpol->SetParName(i, parnames[i].c_str());
        // }

        double parameters[] = {6000, f1270Mass, f1270Width, 6000, a1320Mass, a1320Width, 8000, f1525Mass, f1525Width, 4000, f1710Mass, f1710Width}; // KsKs channel
        // double parameters[] = {6000, f1270Mass, f1270Width, 5000, 1.38, 0.03, 8000, f1525Mass, f1525Width, 4000, f1710Mass, f1710Width}; // KK channel
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {7, 5 * f1525Width}, {20, 20 * f1710WidthErr}}; // KsKs channel
        // vector<vector<float>> par_limits = {{10, 15 * f1710WidthErr}}; // KK channel
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        double initial_param_bkg[] = {3.86556e+02, -5.10155e-02, 1.01604e+01, -2.73381e+00}; // rotational 0-30 GeV/c
        // double initial_param_bkg[] = {343.793, -0.393, 10.072, -1.642}; // like-sign 0-30 GeV/c (KK channel)

        // Initial parameters for background
        BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // Fix
        BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  //Free
        BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // Fix
        BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // Free

        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(2, f1270Width);
        // BEexpol->FixParameter(4, a1320Mass);
        // BEexpol->FixParameter(5, a1320Width);
        // BEexpol->FixParameter(7, f1525Mass);
        // BEexpol->FixParameter(8, f1525Width);
        // BEexpol->FixParameter(10, f1710Mass);
        // BEexpol->FixParameter(11, f1710Width);

        hinvMass->Fit("BEexpol", "REBMS");
        hinvMass->SetMaximum(hinvMass->GetMaximum() * 1.3);
        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        // status codes: 4000 successful, 4 call limit, 4910 failed
        string fitstatus = "Successfull";
        // if (fitResultptr->Status() != 4000)
        // {
        //     cout << "Fit failed or call limit" << endl;
        //     fitstatus = "Failed";
        // }
        // cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
        fitResultptr->Print("V");

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

        TF1 *onlyBW = new TF1("onlyBW", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
        for (int i = 0; i < 12; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParName(i, parnames[i].c_str());
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        onlyBW->Draw("same");

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
            singlefits[i]->SetLineWidth(2);
            // singlefits[i]->Draw("same");
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
        // // // **************** For BW sum with Boltzmann ****************************
#ifdef b_boltzman
        // // int iteration = 0;
        // // for (int ipar1 = 700000; ipar1 < 800000; ipar1 += 10000) // loop for expol parameter 1
        // // {
        // //     for (double ipar2 = 0.4; ipar2 <= 0.8; ipar2 += 0.02) // loop for expol parameter 2
        // //     {
        // //         for (double ipar3 = 3.8; ipar3 < 4.6; ipar3 += 0.1) // loop for expol parameter 3
        // //         {

        // Default fitting range is 1.02 to 2.20. Four types of fitting range variations: extend left (1.0), extend right (2.50), large range (1.0 to 2.50), small range (1.05 to 2.15)

        TF1 *BEexpol = new TF1("BEexpol", BWsum_boltzman_1, 1.03, 2.20, 15); // expol 3
        string parnames[] = {"norm1270", "mass1270", "width1270", "norm1525", "mass1525", "width1525", "norm1710", "mass1710", "width1710", "A", "n", "B"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        double parameters[] = {6000, f1270Mass, f1270Width, 6000, a1320Mass, a1320Width, 8000, f1525Mass, f1525Width, 4000, f1710Mass, f1710Width};
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }
        vector<vector<float>> par_limits = {{1, 3 * f1270Width}, {2, 10 * f1270WidthErr}, {7, 5 * f1525Width}, {10, 20 * f1710WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        BEexpol->SetParameter(size_fitparams + 0, 7.618e5); // expol 3  // 7.618e5  // 5.562e5
        BEexpol->SetParameter(size_fitparams + 1, 0.6456);  // expol 3  //  0.6456  // -0.09379
        BEexpol->SetParameter(size_fitparams + 2, 4.238);   // expol 3  //4.238  // 2.569

        // BEexpol->SetParameter(size_fitparams + 0, ipar1); // expol 3  // 710000  // 5.562e5
        // BEexpol->FixParameter(size_fitparams + 1, ipar2); // expol 3  // -0.03  // -0.09379
        // BEexpol->FixParameter(size_fitparams + 2, ipar3); // expol 3  // 2.78  // 2.569

        // BEexpol->FixParameter(0, 6998);
        // BEexpol->FixParameter(1, f1270Mass);
        // BEexpol->FixParameter(2, f1270Width);
        // BEexpol->FixParameter(3, 7128);
        // BEexpol->FixParameter(4, f1525Mass);
        // BEexpol->FixParameter(5, f1525Width);
        // BEexpol->FixParameter(6, 4058);
        // BEexpol->FixParameter(7, f1710Mass);
        // BEexpol->FixParameter(8, f1710Width);

        // BEexpol->FixParameter(1, a1320Mass);
        // BEexpol->FixParameter(2, a1320Width);

        TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "REBMS");
        cout << "chi2/ndf is " << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
        string fitstatus = "Successfull";
        // cout<<"fit status code "<<fitResultptr->Status()<<endl;
        if (fitResultptr->Status() != 4140)
        {
            cout << "Fit failed or call limit !!!!!!!" << endl;
            fitstatus = "Failed";
        }

        // //             chi2ndf = BEexpol->GetChisquare() / BEexpol->GetNDF();
        // //             fit_parameters.push_back(make_tuple(chi2ndf, ipar1, ipar2, ipar3));
        // //             iteration++;
        // //             cout << "Iteration: " << iteration << endl;
        // //         }
        // //     }
        // // }

        // // // sort in asceding order w.r.t to the third array i.e. chi2/NDF
        // // sort(fit_parameters.begin(), fit_parameters.end(),
        // //      [](const auto &a, const auto &b)
        // //      {
        // //          return get<0>(a) < get<0>(b);
        // //      });
        //
        // // for (int i = 0; i < 20; i++)
        // // {
        // //     float best_ipar1 = std::get<1>(fit_parameters[i]);
        // //     float best_ipar2 = std::get<2>(fit_parameters[i]);
        // //     float best_ipar3 = std::get<3>(fit_parameters[i]);
        // //     float best_chi2ndf = std::get<0>(fit_parameters[i]);
        // //     cout << "ipar1: " << best_ipar1 << ",  ipar2: " << best_ipar2 << ",  ipar3: " << best_ipar3 << ", chi2/NDF: " << best_chi2ndf << endl;
        // // }

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

        TF1 *onlyBW = new TF1("onlyBW", BWsum_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsum_hera, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);

        for (int i = 0; i < 12; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParName(i, parnames[i].c_str());
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        onlyBW->Draw("same");

        TLegend *ltemp = new TLegend(0.20, 0.67, 0.52, 0.92);
        ltemp->SetFillStyle(0);
        ltemp->SetTextFont(42);
        ltemp->SetTextSize(0.035);
        ltemp->AddEntry(hinvMass, "Data", "lpe");
        ltemp->AddEntry(BEexpol, "4rBw + Boltzmann", "l");
        ltemp->AddEntry(onlyBW, "4rBw", "l");
        ltemp->AddEntry(expol, "Boltzmann", "l");
        ltemp->Draw("same");
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
        c->SaveAs((savepath + Form("/rBWfit_pt_%.2f_%.2f.png", pT_bins[ipt], pT_bins[ipt + 1])).c_str());

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

        file << fitstatus << "  " << "4rBW fits , pT range " << pT_bins[ipt] << " - " << pT_bins[ipt + 1] << " GeV/c" << endl;
        file << "fit range, fit params f1270 (norm, mass, and width), a1320, f1525, f1710, significance (f1710), chi2/ndf" << endl;
        file << std::fixed << std::setprecision(2);
        file << "Fit range: " << fitrangelow << " - " << fitrangehigh << endl;
        file << std::fixed << std::setprecision(0);
        file << fitnorm1270 << " ± " << fitnorm1270_err << endl;
        file << std::fixed << std::setprecision(1);
        file << fitmass1270 * 1000.0 << " \\pm " << fitmass1270_err * 1000.0 << endl;
        file << fitwidth1270 * 1000.0 << " \\pm " << fitwidth1270_err * 1000.0 << endl;
        file << std::fixed << std::setprecision(0);
        file << fitnorm1320 << " \\pm " << fitnorm1320_err << endl;
        file << std::fixed << std::setprecision(1);
        file << fitmass1320 * 1000.0 << " \\pm " << fitmass1320_err * 1000.0 << endl;
        file << fitwidth1320 * 1000.0 << " \\pm " << fitwidth1320_err * 1000.0 << endl;
        file << std::fixed << std::setprecision(0);
        file << fitnorm1525 << " \\pm " << fitnorm1525_err << endl;
        file << std::fixed << std::setprecision(1);
        file << fitmass1525 * 1000.0 << " \\pm " << fitmass1525_err * 1000.0 << endl;
        file << fitwidth1525 * 1000.0 << " \\pm " << fitwidth1525_err * 1000.0 << endl;

        // f0(1710) resonance
        file << endl;
        file << std::fixed << std::setprecision(1);
        file << "Chi2/NDF" << endl;
        file << BEexpol->GetChisquare() / BEexpol->GetNDF() << endl;
        file << std::fixed << std::setprecision(0);
        file << fitnorm1710 << " ± " << fitnorm1710_err << endl;
        file << std::fixed << std::setprecision(1);
        file << fitmass1710 * 1000.0 << " ± " << fitmass1710_err * 1000.0 << endl;
        file << fitwidth1710 * 1000.0 << " ± " << fitwidth1710_err * 1000.0 << endl;

        // file << std::fixed << std::setprecision(1);
        // file << significance << endl;
        // file << chi2ndf << endl;
        // file << endl;
        // file << endl;
        // file << "residual bkg parameters" << endl;
        // file << std::fixed << std::setprecision(3);
        // file << resdual_bkg_par1 << ", " << resdual_bkg_par2 << ", " << resdual_bkg_par3 << ", " << resdual_bkg_par4 << endl;
        // file << endl;
        // file << endl;
        // file << std::fixed << std::setprecision(4);
        // file << "fit masses" << endl;
        // file << fitmass1270 << ", " << fitmass1320 << ", " << fitmass1525 << ", " << fitmass1710 << endl;
        // file << "fit masses errors" << endl;
        // file << fitmass1270_err << ", " << fitmass1320_err << ", " << fitmass1525_err << ", " << fitmass1710_err << endl;
        // file << "fit widths" << endl;
        // file << fitwidth1270 << ", " << fitwidth1320 << ", " << fitwidth1525 << ", " << fitwidth1710 << endl;
        // file << "fit widths errors" << endl;
        // file << fitwidth1270_err << ", " << fitwidth1320_err << ", " << fitwidth1525_err << ", " << fitwidth1710_err << endl;

#ifdef residual_subtracted
        // Now subtract the residual background and plot
        TCanvas *c2 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c2, 0.14, 0.03, 0.05, 0.14);
        expol_clone->SetRange(0.99, 2.99);
        hsubtracted->Add(expol_clone, -1);
        hsubtracted->GetXaxis()->SetRangeUser(1.0, 2.5);
        hsubtracted->SetMaximum(hsubtracted->GetMaximum() * 1.5);
        hsubtracted->Draw();
        hsubtracted->Write("3BW");
        // for (int i = 0; i < limits_size; i++)
        // {
        //     int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
        //     onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        // }
        hsubtracted->Fit("onlyBW_clone", "REBMS");
        double *obtained_parameters2 = onlyBW_clone->GetParameters();
        TLine *line = new TLine(0.99, 0, 2.5, 0);
        line->SetLineColor(28);
        line->SetLineStyle(2);
        line->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits1[4];
        for (int i = 0; i < 4; i++)
        {
            singlefits1[i] = new TF1(Form("singlef%d", i), single_BW, BEexpol->GetXmin(), BEexpol->GetXmax(), 3);
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
        c2->SaveAs((savepath + "/rBWfit_residual_allfreeParams.png").c_str());
        // c2->Close();
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
        hinvMass->GetXaxis()->SetRangeUser(1, 2.15);
        hinvMass->SetMarkerSize(0.8);
        hinvMass->SetLineColor(1);
        hinvMass->SetMarkerColor(1);
        hinvMass->SetStats(0);
        hinvMass->SetMinimum(-40000);
        hinvMass->SetMaximum(0.9e6);
        hinvMass->GetYaxis()->SetMaxDigits(4);
        hinvMass->GetYaxis()->SetNdivisions(506);
        hinvMass->Draw("pe");
        expol->Draw("same");
        onlyBW->SetLineWidth(0);
        onlyBW->SetFillColor(4);
        onlyBW->SetFillStyle(1001);
        onlyBW->Draw("same");
        onlyBW->SetNpx(1000);
        gPad->Update();

        TArrow *arrow = new TArrow(0.3537604, 0.6235632, 0.3537604, 0.4137931, 0.02, "|>");
        arrow->SetNDC();
        // arrow->SetLineColor(1); // Set arrow color
        // arrow->SetFillColor(1); // Set arrow fill color
        arrow->SetLineWidth(2); // Set arrow width
        arrow->Draw();

        TLatex *text1 = new TLatex(0.362117, 0.6724138, "f_{2}(1270)/a_{2}^{0}(1320)");
        text1->SetNDC();
        text1->SetTextSize(0.05);
        text1->SetTextAlign(22);
        text1->Draw("same");

        TArrow *arrow2 = new TArrow(0.5125348, 0.5360345, 0.5125348, 0.3262644, 0.02, "|>");
        arrow2->SetNDC();
        arrow2->SetLineWidth(2);
        arrow2->Draw();

        TLatex *text2 = new TLatex(0.5097493, 0.5977011, "f'_{2}(1525)");
        text2->SetNDC();
        text2->SetTextSize(0.05);
        text2->SetTextAlign(22);
        text2->Draw("same");

        TArrow *arrow3 = new TArrow(0.6601671, 0.3965517, 0.6601671, 0.1867816, 0.02, "|>");
        arrow3->SetNDC();
        arrow3->SetLineWidth(2);
        arrow3->Draw();

        TLatex *text3 = new TLatex(0.6587744, 0.4281609, "f_{0}(1710)");
        text3->SetNDC();
        text3->SetTextSize(0.05);
        text3->SetTextAlign(22);
        text3->Draw("same");

        TLegend *leg = new TLegend(0.65, 0.47, 0.99, 0.77);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.06);
        leg->SetBorderSize(0);
        leg->AddEntry(hinvMass, "Data (stat. uncert.)", "lpe");
        leg->AddEntry(BEexpol, "4 BW + Residual BG", "l");
        leg->AddEntry(onlyBW, "Signal", "f");
        leg->AddEntry(expol, "Residual BG", "l");
        leg->Draw("same");

        // TLatex *text4 = new TLatex(0.65, 0.80, "ALICE, work in progress");
        // text4->SetNDC();
        // text4->SetTextSize(0.06);
        // text4->SetTextFont(42);
        // text4->Draw("same");

        TLatex lat4;
        lat4.SetNDC();
        lat4.SetTextSize(0.06);
        lat4.SetTextFont(42);
        lat4.DrawLatex(0.65, 0.80, "ALICE, work in progress");

        c1->cd(2);
        gPad->SetTickx(1);
        hsubtracted->GetYaxis()->SetTitleSize(0.04 / pad2Size);
        hsubtracted->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        hsubtracted->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        hsubtracted->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        hsubtracted->GetYaxis()->SetNdivisions(504);
        hsubtracted->GetXaxis()->SetRangeUser(1.0, 2.15);
        hsubtracted->SetStats(0);
        hsubtracted->SetMinimum(0);
        hsubtracted->SetMaximum(0.13e6);
        hsubtracted->SetMarkerSize(0.8);
        // hsubtracted->GetYaxis()->SetMaxDigits(10);
        hsubtracted->Draw("pe");

        // TLatex *text5 = new TLatex(0.55, 0.80, "Residual background subtraction");
        // text5->SetNDC();
        // text5->SetTextSize(0.06);
        // text5->SetTextFont(42);
        // text5->Draw("same");

        TLatex lat5;
        lat5.SetNDC();
        lat5.SetTextSize(0.06);
        lat5.SetTextFont(42);
        lat5.DrawLatex(0.59, 0.80, "Residual background subtraction");
        lat5.DrawLatex(0.69, 0.70, "pp #sqrt{#it{s}} = 13.6 TeV");
        lat5.DrawLatex(0.69, 0.60, "FT0M (0-100%), |#it{y}|<0.5");
        lat5.DrawLatex(0.69, 0.50, Form("%.1f < #it{p}_{T} < %.1f GeV/c", pT_bins[ipt], pT_bins[ipt + 1]));

        for (int i = 0; i < 4; i++)
        {
            singlefits1[i]->Draw("same");
        }
        c1->Update();
        // c1->SaveAs("/home/sawan/Music/r4BWfit_doublepanel.png");
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
    // plots_4BW->Close();
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

Double_t BWsum_hera(double *x, double *par) // taken 4 resonances here
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

    double fit1 = yield1270 * (realnum1270 * realnum1270 + imagnum1270 * imagnum1270) / (den1270 * den1270);
    double fit2 = yield1320 * (realnum1320 * realnum1320 + imagnum1320 * imagnum1320) / (den1320 * den1320);
    double fit3 = yield1525 * (realnum1525 * realnum1525 + imagnum1525 * imagnum1525) / (den1525 * den1525);
    double fit4 = yield1710 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / (den1710 * den1710);
    double fit = fit1 + fit2 + fit3 + fit4;

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

    double fit = fit1270 + fit1320 + fit1525 + fit1710;
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
    return (BWsum(x, par) + Boltzmann_bkg_1(x, &par[12]));
}

Double_t BWsum_boltzman_2(double *x, double *par)
{
    return (BWsum(x, par) + Boltzmann_bkg_2(x, &par[12]));
}
Double_t expol_chkstar(double *x, double *par)
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] + x[0] * par[3]));
}
Double_t BWsum_expol_chkstar(double *x, double *par)
{
    return (BWsum(x, par) + expol_chkstar(x, &par[12]));
}

Double_t BWsumMassDepWidth_exponential(double *x, double *par)
{
    // return (BWsumMassDepWidth(x, par) + expol_chkstar(x, &par[12]));
    return (BWsumMassDepWidth(x, par) + exponential_bkg_3(x, &par[12]));
}

Double_t BWsumMassDepWidth_expol1(double *x, double *par)
{
    return (BWsumMassDepWidth(x, par) + expol_chkstar(x, &par[12]));
    // return (BWsumMassDepWidth(x, par) + exponential_bkg_3(x, &par[12]));
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
