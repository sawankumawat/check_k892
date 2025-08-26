#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

Double_t single_BW(double *x, double *par);
Double_t BWsum(double *x, double *par);
Double_t BWsumMassDepWidth(double *x, double *par);
Double_t BWsum_boltzmann_massdep(double *x, double *par);
Double_t single_BW_mass_dep_spin0(double *x, double *par);
Double_t single_BW_mass_dep_spin2(double *x, double *par);

Double_t exponential_bkg_3(double *x, double *par); // 4 parameters

Double_t single_BW_expol3(double *x, double *par);
Double_t BWsum_expol3(double *x, double *par);
Double_t BWsum_expol_chkstar(double *x, double *par);
Double_t BWsumMassDepWidth_expol3(double *x, double *par);

Double_t single_BW_boltzman_1(double *x, double *par);
Double_t single_BW_boltzman_2(double *x, double *par);

void fit_check()
{
    ofstream file;
    file.open("fit_params.txt");
    string savepath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/fits/3rBW_fits";

    TFile *f = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60_fullpt.root", "READ"); // full pT range

    int colors[] = {4, 6, 28, 46};
    double masses[] = {f1270Mass, f1525Mass, f1710Mass};
    double widths[] = {f1270Width, f1525Width, f1710Width};
    string resonance_names[] = {"f_{2}(1270)", "f'_{2}(1525)", "f_{0}(1710)"};
    double purity, significance[3], chi2ndf, chi2, ndf;
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);
    t2->SetTextFont(42);

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

#define b_massdepwidth
#define residual_subtracted

    for (int ipt = 0; ipt < Npt; ipt++)
    {

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
        TF1 *onlyBW_clone;

// // Default fitting range is 1.02 to 2.20. Four types of fitting range variations: extend left (1.0), extend right (2.50), large range (1.0 to 2.50), small range (1.05 to 2.15)
#ifdef b_massdepwidth
        // TF1 *BEexpol = new TF1("BEexpol", BWsum_expol3, 1.05, 2.20, 13); // expol 3
        TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_expol3, 1.05, 2.20, 13); // expol 3
        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        double parameters[] = {6200, f1270Mass, f1270Width, 7700, f1525Mass, f1525Width, 2500, f1710Mass, f1710Width}; // for preview presentation
        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }

        vector<vector<float>> par_limits = {{2, 10 * f1270WidthErr}};
        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        double initial_param_bkg[] = {5.562e5, -0.09379, 2.569, 1.0982}; // rotational 1-30 GeV/c (KsKs channel)

        BEexpol->SetParameter(size_fitparams + 0, initial_param_bkg[0]); // 5.562e5   // 206 //5.845e5
        BEexpol->SetParameter(size_fitparams + 1, initial_param_bkg[1]); // -0.09379  // 0.04316 //-0.07378
        BEexpol->SetParameter(size_fitparams + 2, initial_param_bkg[2]); // 2.569     // 11.48 //2.685
        BEexpol->SetParameter(size_fitparams + 3, initial_param_bkg[3]); // 1.098     // -3.149 //1.176

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

        TF1 *onlyBW = new TF1("onlyBW", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        onlyBW_clone = new TF1("onlyBW_clone", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 9);
        onlyBW_clone->SetParNames("f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma");
        for (int i = 0; i < 9; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
        }
        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);

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
            significance[i] = significance_num / significance_den;
            cout << "Significance of " << resonance_names[i] << " is " << significance[i] << endl;
            purity = significance_num * 100 / hinvMass->Integral(binlow, binhigh);
            cout << "Purity of " << resonance_names[i] << " is " << purity << endl;
        }
#endif

        // // //********************************* common for all fits ***************************************
        // // //****************************common for all fits************************************
        gPad->Update();
        TPaveStats *ptstats = (TPaveStats *)hinvMass->FindObject("stats");
        ptstats->SetX1NDC(0.6);
        ptstats->SetX2NDC(0.99);
        ptstats->SetY1NDC(0.4);
        ptstats->SetY2NDC(0.92);
        ptstats->Draw("same");
        c->SaveAs((savepath + Form("/rBWfit_pt_%.2f_%.2f.png", pT_bins[ipt], pT_bins[ipt + 1])).c_str());

#ifdef residual_subtracted
        // Now subtract the residual background and plot
        TCanvas *c2 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c2, 0.14, 0.03, 0.05, 0.14);
        expol_clone->SetRange(0.99, 2.99);
        hsubtracted->Add(expol_clone, -1);
        hsubtracted->GetXaxis()->SetRangeUser(1.0, 2.5);
        hsubtracted->SetMaximum(hsubtracted->GetMaximum() * 1.5);
        hsubtracted->Draw();
        TH1F *hsubtracted_clone = (TH1F *)hsubtracted->Clone("hsubtracted_clone");

        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }
        hsubtracted->Fit("onlyBW_clone", "REBMS");
        double *obtained_parameters2 = onlyBW_clone->GetParameters();
        TLine *line = new TLine(0.99, 0, 2.5, 0);
        line->SetLineColor(28);
        line->SetLineStyle(2);
        line->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits1[3];
        for (int i = 0; i < 3; i++)
        {
            singlefits1[i] = (i < 2) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, onlyBW_clone->GetXmin(), onlyBW_clone->GetXmax(), 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, onlyBW_clone->GetXmin(), onlyBW_clone->GetXmax(), 3);

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
        c2->SaveAs((savepath + "/rBWfit_residual.png").c_str());

#endif
    }
}
// end of main program

// *****************************************************************************************************
//************************************Fit functions************************************************* */

// We will define the single BW and the sum of 3 BWs

Double_t single_BW(double *x, double *par)
{
    double yield = par[0];
    double mass = par[1];
    double width = par[2];

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));
    return fit;
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

Double_t exponential_bkg_3(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[3])));
}

// Now we will define the functions for the sum of the Breit-Wigner and the exponential background

Double_t single_BW_expol3(double *x, double *par)
{
    return (single_BW(x, par) + exponential_bkg_3(x, &par[3]));
}

Double_t BWsum_expol3(double *x, double *par)
{
    return (BWsum(x, par) + exponential_bkg_3(x, &par[9]));
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
