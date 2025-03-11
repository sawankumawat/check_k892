#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>
#include <TArrow.h>
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);
// systematic studies (signal extraction) ****************************
// A. fit range: Default: 1.05-2.20 GeV/c^2, Variation1: 1.02-2.20 GeV/c^2, Variation2: 1.05-2.30 GeV/c^2, Variation3: 1.08-2.15 GeV/c^2, Variation4: 1.02-2.30 GeV/c^2
// B. Norm range: Default: 2.50-2.60 GeV/c^2, Variation1: 2.40-2.50 GeV/c^2, Variation2: 2.60-2.70 GeV/c^2
// C. Residual bkg Fit function: Default: 4rBW with mass dependent width + modified Boltzmann, Variation1: 4rBW with constant width + bkg, Variation2: 3rBW with mass dependen width + bkg, Variation3: 4rBW with mass dependent width + Expol1, Variation4: 4rBW with mass dependent width + Boltzmann
// D. Fit paramters: Default: Width of spin-2 resonances fixed to PDG, Variation1: Width of spin-2 resonances free, Variation2: Both mass and width of spin-2 resonances fixed to PDG, Variation3: Width of f1710 fixed to PDG, Variation4: Mass of f1710 fixed to PDG
// E. Combinatorial background: Default: Rotational, Variation1: Mixed (Not considered)

// systematic studies (Track selection) ****************************
// TrA. DCA track to PV: Deafult: 0.05 cm, Variation1: 0.04 cm, Variation2: 0.06 cm
// TrC. TPC crossed rows: Default 70, Variation1: 100, Variation2: 120
// TrD. TPC crossed rows over findable clusters: Default: 0.8, Variation1: 0.9, Variation2: 1.0

// systematic studies (Topological selection) ****************************
// ToA. Cosine PA: Default: 0.97, Variation1: 0.95, Variation2: 0.99
// ToB. Transeverse radius: Default: 0.5 cm, Variation1: 0.4 cm, Variation2: 0.6 cm
// ToC. DCA b/w V0 daughters: Default: 0.5 cm, Variation1: 0.3 cm, Variation2: 1.0 cm
// ToD. Lifetime: Default: 20 cm, Variation1: 15 cm, Variation2: 25 cm
// ToE. Competing V0 rejection: Default: 5 MeV, Variation1: 4, Variation2: 6
// ToF. Ks mass window: Default: 3sigma, Variation1: 4sigma, Variation2: 5sigma

// systematic studies (PID) ****************************
// TrB. TPC PID: Default: 3sigma, Variation1: 4sigma, Variation2: 5sigma

using namespace std;

void systematics()
{
    gStyle->SetOptStat(0);
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/351471/KsKs_Channel/higher-mass-resonances_3sigmaKs/fits/4rBw_fits/systematics_onlywidthfixed/";
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937/fits/4rBw_fits/backup/";
    ifstream file;

    // int colors[] = {kBlue, kRed, kGreen + 3, kOrange + 7, kMagenta + 2, kCyan + 2, kBlack};
    vector<int> colors = {kBlue, kRed, kGreen + 3, kOrange + 7, kMagenta + 2, kCyan + 2, kBlack};
    int markers[] = {20, 21, 22, 23, 24, 25, 26};

    ofstream sysfile;
    sysfile.open(path + "temp.txt");

    // string all_variations[] = {"Signal extraction", "PID", "Track selection", "Topological selection"};

    // // // // Mapping variations to sources
    // // *********************signal extraction*******************************
    // int variations[] = {1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
    // int variationsNumber[] = {1, 4, 2, 4, 4};                                                                                                                         // Number of variations for each source
    // vector<string> names = {"default", "varA1", "varA2", "varA3", "varA4", "varB1", "varB2", "varC1", "varC2", "varC3", "varC4", "varD1", "varD2", "varD3", "varD4"}; // signal extraction
    // vector<string> varNames = {"Fit range", "Normalization range", "Fit function", "Fit parameters"};                                                                 // signal extraction
    // vector<vector<string>> individual_sources{{"1.05-2.20 GeV/c^{2} (def)", "1.02-2.20 GeV/c^{2}", "1.05-2.30 GeV/c^{2}", "1.08-2.15 GeV/c^{2}", "1.02-2.30 GeV/c^{2}"}, {"2.50-2.60 GeV/c^{2} (def)", "2.40-2.50 GeV/c^{2}", "2.60-2.70 GeV/c^{2}"}, {"Fit 1 (def)", "Fit 2", "Fit 3", "Fit 4", "Fit 5"}, {"#Gamma fix to PDG (spin-2) (def)", "All parameter free", "#Gamma & M fixed (spin-2)", "#Gamma f_{0}(1710) fix to PDG", "M f_{0} fix to PDG"}};
    // string combinedName = "Signal Extraction";

    // // *****************PID********************************
    // int variations[] = {1, 2, 2};
    // int variationsNumber[] = {1, 2};
    // vector<string> names = {"default", "varTrB1", "varTrB2"}; // PID
    // vector<string> varNames = {"TPC PID"};                    // PID
    // vector<vector<string>> individual_sources{{"3#sigma (def)", "4#sigma", "5#sigma"}};
    // string combinedName = "PID";

    // // //*****************Track selection**************************
    // int variations[] = {1, 2, 2, 3, 3, 4, 4};
    // int variationsNumber[] = {1, 2, 2, 2};
    // vector<string> names = {"default", "varTrA1", "varTrA2", "varTrC1", "varTrC2", "varTrD1", "varTrD2"};         // Track selection
    // vector<string> varNames = {"DCA track to PV", "TPC crossed rows", "TPC crossed rows over findable clusters"}; // Track selection
    // vector<vector<string>> individual_sources{{"0.05 cm (def)", "0.04 cm", "0.06 cm"}, {"70 (def)", "100", "120"}, {"0.8 (def)", "0.9", "1.0"}};
    // string combinedName = "Track Selection";

    // // // ***************Topological selection***************
    // int variations[] = {1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};
    // int variationsNumber[] = {1, 2, 2, 2, 2, 2, 2};
    // vector<string> names = {"default", "varToA1", "varToA2", "varToB1", "varToB2", "varToC1", "varToC2", "varToD1", "varToD2", "varToE1", "varToE2", "varToF1", "varToF2"}; // Topological selection
    // vector<string> varNames = {"Cosine PA", "Transverse radius", "DCA b/w V0 daughters", "Lifetime", "Competing V0 rejection", "Ks mass window"};                           // Topological selection
    // vector<vector<string>> individual_sources{{"0.97 (def)", "0.95", "0.99"}, {"0.5 cm (def)", "0.4 cm", "0.6 cm"}, {"0.5 cm (def)", "0.3 cm", "1.0 cm"}, {"20 cm (def)", "15 cm", "25 cm"}, {"5 MeV (def)", "4", "6"}, {"3#sigma (def)", "4#sigma", "5#sigma"}}; // Topological selection
    // string combinedName = "Topological Selection";

    // // ********All sources******************
    int variations[] = {1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15};
    int variationsNumber[] = {1, 4, 2, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    vector<string> names = {"default", "varA1", "varA2", "varA3", "varA4", "varB1", "varB2", "varC1", "varC2", "varC3", "varC4", "varD1", "varD2", "varD3", "varD4", "varTrB1", "varTrB2", "varTrA1", "varTrA2", "varTrC1", "varTrC2", "varTrD1", "varTrD2", "varToA1", "varToA2", "varToB1", "varToB2", "varToC1", "varToC2", "varToD1", "varToD2", "varToE1", "varToE2", "varToF1", "varToF2"}; // All sources
    vector<string> varNames = {"Fit range", "Normalization range", "Fit function", "Fit parameters", "TPC PID", "DCA track to PV", "TPC crossed rows", "TPC RowsClusters", "Cosine PA", "Transverse Radius", "DCA V0 Daughters", "Lifetime", "Competing V0 Rejection", "Ks Mass Window"};
    vector<vector<string>> individual_sources{{"1.05-2.20 GeV/c^{2} (def)", "1.02-2.20 GeV/c^{2}", "1.05-2.30 GeV/c^{2}", "1.08-2.15 GeV/c^{2}", "1.02-2.30 GeV/c^{2}"}, {"2.50-2.60 GeV/c^{2} (def)", "2.40-2.50 GeV/c^{2}", "2.60-2.70 GeV/c^{2}"}, {"Fit 1 (def)", "Fit 2", "Fit 3", "Fit 4", "Fit 5"}, {"#Gamma fix to PDG (spin-2) (def)", "All parameter free", "#Gamma & M fixed (spin-2)", "#Gamma f_{0}(1710) fix to PDG", "M f_{0} fix to PDG"}, {"3#sigma (def)", "4#sigma", "5#sigma"}, {"0.05 cm (def)", "0.04 cm", "0.06 cm"}, {"70 (def)", "100", "120"}, {"0.8 (def)", "0.9", "1.0"}, {"0.97 (def)", "0.95", "0.99"}, {"0.5 cm (def)", "0.4 cm", "0.6 cm"}, {"0.5 cm (def)", "0.3 cm", "1.0 cm"}, {"20 cm (def)", "15 cm", "25 cm"}, {"5 MeV (def)", "4", "6"}, {"3#sigma (def)", "4#sigma", "5#sigma"}}; // All sources
    string combinedName = "All Sources";

    const int size = varNames.size();
    // Vectors to store fit parameters and their uncertainties
    // get the size from the variations_map

    vector<vector<double>> allparams_mass(size + 1), allparams_width(size + 1);         // To store values
    vector<vector<double>> allparams_mass_unc(size + 1), allparams_width_unc(size + 1); // To store uncertainties
    // total variations size taking from the variations_map

    // Reading files
    for (int ivars = 0; ivars < names.size(); ivars++)
    {
        file.open(path + Form("fit_params_%s.txt", names[ivars].c_str()));
        if (!file.is_open())
        {
            cerr << "Error opening file: " << names[ivars] << endl;
            return;
        }

        string line;
        vector<double> values, uncertainties;

        while (getline(file, line))
        {
            istringstream iss(line);
            double num, unc;
            string pm; // String to capture ±

            if (iss >> num >> pm >> unc) // Read three parts
            {
                values.push_back(num);
                uncertainties.push_back(unc);
            }
        }

        file.close();
        // std::setprecision(2);
        // Ensure at least two values were read
        if (values.size() >= 2)
        {
            cout << names[ivars] << ": "
                 << Form("%.2f", values[values.size() - 2] * 1000) << " ± " << Form("%.2f", uncertainties[values.size() - 2] * 1000) << " & "
                 << Form("%.2f", values[values.size() - 1] * 1000) << " ± " << Form("%.2f", uncertainties[values.size() - 1] * 1000) << endl;
            
            // Store the second-last value and its uncertainty
            int varIndex = variations[ivars] - 1;
            allparams_mass[varIndex].push_back(values[values.size() - 2]);
            allparams_width[varIndex].push_back(values[values.size() - 1]);

            allparams_mass_unc[varIndex].push_back(uncertainties[values.size() - 2]);
            allparams_width_unc[varIndex].push_back(uncertainties[values.size() - 1]);
        }
        else
        {
            cerr << "Not enough numeric values found in " << names[ivars] << "!" << endl;
        }
    }

    // ofstream txtfile;
    // txtfile.open("/home/sawan/Videos/fit_params.txt");

    // Compute relative uncertainties
    vector<double> relative_uncertainty_mass, relative_uncertainty_width;
    for (int i = 1; i < allparams_mass.size(); i++) // Skip default (index 0)
    {
        double sum_mass{0}, sum_width{0};
        if (allparams_mass[i].size() != variationsNumber[i])
        {
            cerr << "Mismatch in expected variations count for source " << i << endl;
            continue;
        }
        cout << "Source " << i << ", number of variations: " << variationsNumber[i] << endl;

        for (int j = 0; j < allparams_mass[i].size(); j++)
        {
            double diff_mass = (allparams_mass[i][j] - allparams_mass[0][0]) / allparams_mass[0][0];
            double diff_width = (allparams_width[i][j] - allparams_width[0][0]) / allparams_width[0][0];
            sum_mass += pow(diff_mass, 2);
            sum_width += pow(diff_width, 2);
            cout << "default value mass: " << allparams_mass[0][0] << " ± " << allparams_mass_unc[0][0]
                 << "  variation value: " << allparams_mass[i][j] << " ± " << allparams_mass_unc[i][j] << endl;
            cout << "default value width: " << allparams_width[0][0] << " ± " << allparams_width_unc[0][0]
                 << "  variation value: " << allparams_width[i][j] << " ± " << allparams_width_unc[i][j] << endl;
        }

        double relUnc_mass = sqrt(sum_mass / variationsNumber[i]);
        relative_uncertainty_mass.push_back(relUnc_mass);
        double relUnc_width = sqrt(sum_width / variationsNumber[i]);
        relative_uncertainty_width.push_back(relUnc_width);

        cout << "Relative Uncertainty (Source " << i << "): " << relUnc_mass << endl;
    }

    // Compute quadrature sum of uncertainties
    double total_uncertainty_mass{0}, total_uncertainty_width{0};
    for (double ru : relative_uncertainty_mass)
    {
        total_uncertainty_mass += pow(ru, 2);
    }
    total_uncertainty_mass = sqrt(total_uncertainty_mass);
    cout << "Total Quadrature Sum Uncertainty (Mass): " << total_uncertainty_mass << endl;

    for (double ru : relative_uncertainty_width)
    {
        total_uncertainty_width += pow(ru, 2);
    }
    total_uncertainty_width = sqrt(total_uncertainty_width);
    cout << "Total Quadrature Sum Uncertainty (Width): " << total_uncertainty_width << endl;

    sysfile << "Total Quadrature Sum Uncertainty (Mass): " << total_uncertainty_mass << endl;
    sysfile << "Total Quadrature Sum Uncertainty (Width): " << total_uncertainty_width << endl;
    sysfile.close();

    // plots***********************
    /*

    // // now lets plot the the variations
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c1, 0.18, 0.03, 0.05, 0.14);
    TH1F *hmass_var[size], *hwidth_var[size], *hmass_width[size];
    for (int i = 0; i < size; i++)
    {
        hmass_var[i] = new TH1F(Form("hmass_var%d", i), "", 1, 0, 30);
        hmass_width[i] = new TH1F(Form("hmass_width%d", i), "", 2, 0, 2);
        hmass_var[i]->SetBinContent(1, relative_uncertainty_mass[i]);
        hmass_width[i]->SetBinContent(1, relative_uncertainty_mass[i]);
        hmass_width[i]->SetBinContent(2, relative_uncertainty_width[i]);
        SetHistoQA(hmass_var[i]);
        SetHistoQA(hmass_width[i]);
        hmass_var[i]->SetLineColor(colors[i]);
        hmass_width[i]->SetLineColor(colors[i]);
        hmass_var[i]->SetLineWidth(4);
        hmass_var[i]->SetLineStyle(2);
        hmass_var[i]->GetYaxis()->SetRangeUser(0, 0.0151);
        hmass_var[i]->GetYaxis()->SetTitleOffset(1.8);
        hmass_var[i]->GetYaxis()->SetTitle("Fractional Uncertainty (Mass)");
        hmass_var[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hmass_var[i]->Draw("same");
    }
    TH1F *hmass_total = new TH1F("hmass_total", "", 1, 0, 30);
    hmass_total->SetBinContent(1, total_uncertainty_mass);
    SetHistoQA(hmass_total);
    hmass_total->SetLineColor(1);
    hmass_total->SetLineWidth(4);
    // hmass_total->SetLineStyle(2);
    hmass_total->Draw("same");

    TLegend *l1 = new TLegend(0.25, 0.60, 0.55, 0.91);
    l1->SetFillStyle(0);
    l1->SetBorderSize(0);
    l1->SetTextFont(42);
    l1->SetTextSize(0.04);
    l1->SetHeader(combinedName.c_str());
    for (int i = 0; i < size; i++)
    {
        l1->AddEntry(hmass_var[i], varNames[i].c_str(), "l");
    }
    l1->AddEntry(hmass_total, "Total", "l");
    l1->Draw("same");

    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.18, 0.03, 0.05, 0.14);
    for (int i = 0; i < size; i++)
    {
        hwidth_var[i] = new TH1F(Form("hwidth_var%d", i), "", 1, 0, 30);
        hwidth_var[i]->SetBinContent(1, relative_uncertainty_width[i]);
        SetHistoQA(hwidth_var[i]);
        hwidth_var[i]->SetLineColor(colors[i]);
        hwidth_var[i]->SetLineWidth(4);
        hwidth_var[i]->SetLineStyle(2);
        hwidth_var[i]->GetYaxis()->SetRangeUser(0, 0.351);
        hwidth_var[i]->GetYaxis()->SetTitleOffset(1.8);
        hwidth_var[i]->GetYaxis()->SetTitle("Fractional Uncertainty (Width)");
        hwidth_var[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hwidth_var[i]->Draw("same");
    }
    TH1F *hwidth_total = new TH1F("hwidth_total", "", 1, 0, 30);
    hwidth_total->SetBinContent(1, total_uncertainty_width);
    SetHistoQA(hwidth_total);
    hwidth_total->SetLineColor(1);
    hwidth_total->SetLineWidth(4);
    hwidth_total->GetYaxis()->SetRangeUser(0, 0.351);
    hwidth_total->Draw("same");

    TLegend *l2 = new TLegend(0.25, 0.60, 0.55, 0.91);
    l2->SetFillStyle(0);
    l2->SetBorderSize(0);
    l2->SetTextFont(42);
    l2->SetTextSize(0.04);
    l2->SetHeader(combinedName.c_str());
    for (int i = 0; i < size; i++)
    {
        l2->AddEntry(hwidth_var[i], varNames[i].c_str(), "l");
    }
    l2->AddEntry(hwidth_total, "Total", "l");
    l2->Draw("same");

    c1->SaveAs(("/home/sawan/Videos/systematic_mass_" + combinedName + ".png").c_str());
    c2->SaveAs(("/home/sawan/Videos/systematic_width_" + combinedName + ".png").c_str());

    TCanvas *c3 = new TCanvas("c", "", 720, 720);
    SetCanvasStyle(c3, 0.18, 0.03, 0.001, 0.08);
    c3->Divide(2, 1, 0, 0); // No gap between pads

    c3->cd(1);
    gPad->SetLeftMargin(0.27);
    gPad->SetTopMargin(0.05);
    for (int i = 0; i < size; i++)
    {
        hmass_var[i]->GetXaxis()->SetBinLabel(1, "Mass (GeV/#it{c}^{2})");
        hmass_var[i]->GetXaxis()->SetLabelSize(0.1);
        hmass_var[i]->GetYaxis()->SetLabelSize(0.06);
        hmass_var[i]->GetYaxis()->SetTitleSize(0.07);
        hmass_var[i]->GetYaxis()->SetTitleOffset(1.9);
        hmass_var[i]->GetYaxis()->SetTitle("Fractional Uncertainty");
        hmass_var[i]->GetXaxis()->SetTitle(0);
        hmass_var[i]->Draw("same");
    }
    hmass_total->Draw("same");

    // Right plot (Width)
    c3->cd(2);
    gPad->SetRightMargin(0.28);
    gPad->SetTopMargin(0.05);
    SetHistoQA(hwidth_total);
    hwidth_total->GetXaxis()->SetBinLabel(1, "Width (GeV/#it{c}^{2})");
    hwidth_total->GetXaxis()->SetLabelSize(0.12);
    hwidth_total->GetYaxis()->SetLabelSize(0.07);
    hwidth_total->GetYaxis()->SetTitleSize(0.09);
    hwidth_total->GetYaxis()->SetTitleOffset(1.6);
    hwidth_total->SetLineWidth(4);
    hwidth_total->GetYaxis()->SetTitle("Fractional Uncertainty");
    hwidth_total->Draw("Y+");

    for (int i = 0; i < size; i++)
    {
        hwidth_var[i]->GetXaxis()->SetBinLabel(1, "Width (GeV/#it{c}^{2})");
        hwidth_var[i]->GetXaxis()->SetLabelSize(0.1);
        hwidth_var[i]->GetYaxis()->SetLabelSize(0.06);
        hwidth_var[i]->GetYaxis()->SetTitleSize(0.07);
        hwidth_var[i]->GetXaxis()->SetTitle(0);
        hwidth_var[i]->SetLineStyle(4);
        hwidth_var[i]->Draw("same");
    }

    TPad *padfullarea = new TPad("padfullarea", "", 0.0, 0.0, 1.0, 1.0);
    padfullarea->SetFillStyle(0);
    c3->cd();          // Go back to the main canvas
    padfullarea->Draw();
    padfullarea->cd(); // Activate the full canvas pad

    TLegend *l3 = new TLegend(0.18, 0.60, 0.75, 0.91);
    l3->SetFillStyle(0);
    l3->SetBorderSize(0);
    l3->SetTextFont(42);
    l3->SetTextSize(0.03);
    l3->SetHeader(combinedName.c_str());

    for (int i = 0; i < size; i++)
    {
        l3->AddEntry(hmass_var[i], varNames[i].c_str(), "l");
    }
    l3->AddEntry(hmass_total, "Total", "l");
    l3->Draw(); // No need for "same"

    c3->SaveAs(("/home/sawan/Videos/systematic_mass_width_doublepanel_" + combinedName + ".png").c_str());

    // // make double panel plots for the ratios
    // TH1F *mass_default = new TH1F("mass_default", "", 1, 0, 30);
    // TH1F *width_default = new TH1F("width_default", "", 1, 0, 30);
    // mass_default->SetBinContent(1, allparams_mass[0][0]);
    // mass_default->SetBinError(1, allparams_mass_unc[0][0]);
    // width_default->SetBinContent(1, allparams_width[0][0]);
    // width_default->SetBinError(1, allparams_width_unc[0][0]);
    // SetHistoQA(mass_default);
    // SetHistoQA(width_default);

    // // for mass
    // for (int i = 0; i < size; i++)
    // {
    //     TCanvas *cmassvar = new TCanvas(Form("cmassvar%d", i), "", 720, 720);
    //     double pad1Size, pad2Size;
    //     canvas_style(cmassvar, pad1Size, pad2Size);
    //     cmassvar->cd(1);
    //     mass_default->GetYaxis()->SetRangeUser(1.661, 1.779);
    //     mass_default->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    //     mass_default->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    //     mass_default->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    //     mass_default->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    //     mass_default->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    //     mass_default->GetYaxis()->SetTitleOffset(1.2);
    //     mass_default->Draw("pe");
    //     TH1F *mass_var[variationsNumber[i + 1]];
    //     for (int j = 0; j < variationsNumber[i + 1]; j++)
    //     {
    //         mass_var[j] = new TH1F(Form("mass_var%d", i), "", 1, 0, 30);
    //         mass_var[j]->SetBinContent(1, allparams_mass[i + 1][j]);
    //         mass_var[j]->SetBinError(1, allparams_mass_unc[i + 1][j]);
    //         SetHistoQA(mass_var[j]);
    //         mass_var[j]->SetLineColor(colors[j]);
    //         mass_var[j]->SetMarkerStyle(markers[j]);
    //         mass_var[j]->SetMarkerColor(colors[j]);
    //         mass_var[j]->Draw("pesame");
    //     }

    //     TLegend *leg = new TLegend(0.20, 0.5, 0.55, 0.9);
    //     leg->SetFillStyle(0);
    //     leg->SetBorderSize(0);
    //     leg->SetTextFont(42);
    //     leg->SetTextSize(0.04 / pad1Size);
    //     leg->SetHeader(varNames[i].c_str());
    //     leg->AddEntry(mass_default, individual_sources[i][0].c_str(), "lpe");
    //     for (int j = 0; j < variationsNumber[i + 1]; j++)
    //     {
    //         leg->AddEntry(mass_var[j], individual_sources[i][j + 1].c_str(), "lpe");
    //     }
    //     leg->Draw("same");

    //     cmassvar->cd(2);
    //     for (int j = 0; j < variationsNumber[i + 1]; j++)
    //     {
    //         TH1F *mass_ratio = (TH1F *)mass_default->Clone();
    //         TH1F *mass_var = new TH1F(Form("mass_var%d", i), "", 1, 0, 30);
    //         mass_var->SetBinContent(1, allparams_mass[i + 1][j]);
    //         mass_var->SetBinError(1, allparams_mass_unc[i + 1][j]);
    //         mass_ratio->Divide(mass_var);
    //         mass_ratio->SetMarkerStyle(markers[j]);
    //         mass_ratio->SetMarkerColor(colors[j]);
    //         mass_ratio->SetLineColor(colors[j]);
    //         mass_ratio->GetYaxis()->SetRangeUser(0.989, 1.011);
    //         mass_ratio->GetYaxis()->SetNdivisions(505);
    //         mass_ratio->GetYaxis()->SetTitle("Ratio");
    //         mass_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    //         mass_ratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    //         mass_ratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    //         mass_ratio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    //         mass_ratio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    //         mass_ratio->GetYaxis()->SetTitleOffset(0.55);
    //         mass_ratio->Draw("pe same");
    //     }
    //     TH1F *ratio1 = (TH1F *)mass_default->Clone();
    //     ratio1->Divide(mass_default);
    //     ratio1->SetMarkerStyle(20);
    //     ratio1->SetMarkerColor(kBlack);
    //     ratio1->SetLineColor(kBlack);
    //     ratio1->Draw("same");
    //     TLine *line = new TLine(0, 1, 30, 1);
    //     line->SetLineColor(kBlack);
    //     line->SetLineStyle(2);
    //     line->Draw("same");

    //     cmassvar->SaveAs(Form("/home/sawan/Videos/mass_var_%s_%d.png", combinedName.c_str(), i));
    // }

    // // for width
    // for (int i = 0; i < size; i++)
    // {
    //     TCanvas *cwidthvar = new TCanvas(Form("cwidthvar%d", i), "", 720, 720);
    //     double pad1Size, pad2Size;
    //     canvas_style(cwidthvar, pad1Size, pad2Size);
    //     cwidthvar->cd(1);
    //     width_default->GetYaxis()->SetRangeUser(0.10, 0.301);
    //     width_default->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    //     width_default->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    //     width_default->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    //     width_default->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    //     width_default->GetYaxis()->SetTitle("Width (GeV/#it{c}^{2})");
    //     width_default->GetYaxis()->SetTitleOffset(1.2);
    //     width_default->Draw("pe");
    //     TH1F *width_var[variationsNumber[i + 1]];
    //     for (int j = 0; j < variationsNumber[i + 1]; j++)
    //     {
    //         width_var[j] = new TH1F(Form("width_var%d", i), "", 1, 0, 30);
    //         width_var[j]->SetBinContent(1, allparams_width[i + 1][j]);
    //         width_var[j]->SetBinError(1, allparams_width_unc[i + 1][j]);
    //         SetHistoQA(width_var[j]);
    //         width_var[j]->SetLineColor(colors[j]);
    //         width_var[j]->SetMarkerStyle(markers[j]);
    //         width_var[j]->SetMarkerColor(colors[j]);
    //         width_var[j]->Draw("pesame");
    //     }

    //     TLegend *leg = new TLegend(0.20, 0.5, 0.55, 0.9);
    //     leg->SetFillStyle(0);
    //     leg->SetBorderSize(0);
    //     leg->SetTextFont(42);
    //     leg->SetTextSize(0.04 / pad1Size);
    //     leg->SetHeader(varNames[i].c_str());
    //     leg->AddEntry(width_default, individual_sources[i][0].c_str(), "lpe");
    //     for (int j = 0; j < variationsNumber[i + 1]; j++)
    //     {
    //         leg->AddEntry(width_var[j], individual_sources[i][j + 1].c_str(), "lpe");
    //     }
    //     leg->Draw("same");

    //     cwidthvar->cd(2);
    //     for (int j = 0; j < variationsNumber[i + 1]; j++)
    //     {
    //         TH1F *width_ratio = (TH1F *)width_default->Clone();
    //         TH1F *width_var = new TH1F(Form("width_var%d", i), "", 1, 0, 30);
    //         width_var->SetBinContent(1, allparams_width[i + 1][j]);
    //         width_var->SetBinError(1, allparams_width_unc[i + 1][j]);
    //         width_ratio->Divide(width_var);
    //         width_ratio->SetMarkerStyle(markers[j]);
    //         width_ratio->SetMarkerColor(colors[j]);
    //         width_ratio->SetLineColor(colors[j]);
    //         width_ratio->GetYaxis()->SetRangeUser(0.59, 1.51);
    //         width_ratio->GetYaxis()->SetNdivisions(505);
    //         width_ratio->GetYaxis()->SetTitle("Ratio");
    //         width_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    //         width_ratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    //         width_ratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    //         width_ratio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    //         width_ratio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    //         width_ratio->GetYaxis()->SetTitleOffset(0.55);
    //         width_ratio->Draw("pe same");
    //     }
    //     TH1F *ratio1 = (TH1F *)width_default->Clone();
    //     ratio1->Divide(width_default);
    //     ratio1->SetMarkerStyle(20);
    //     ratio1->SetMarkerColor(kBlack);
    //     ratio1->SetLineColor(kBlack);
    //     ratio1->Draw("same");
    //     TLine *line = new TLine(0, 1, 30, 1);
    //     line->SetLineColor(kBlack);
    //     line->SetLineStyle(2);
    //     line->Draw("same");

    //     cwidthvar->SaveAs(Form("/home/sawan/Videos/width_var_%s_%d.png", combinedName.c_str(), i));
    // }

    */
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.14);
    pad2->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.04);
}

// vector<string> systematic_sources = {"Fit Range", "Norm Range", "Residual Bkg Fit", "Fit Parameters", "DCA Track to PV", "TPC PID", "TPC Crossed Rows", "TPC RowsClusters", "Cosine PA", "Transverse Radius", "DCA V0 Daughters", "Lifetime", "Competing V0 Rejection", "Ks Mass Window"};

// map<string, vector<string>> variations_map = {
//     {"Fit Range", {"varA1", "varA2", "varA3", "varA4"}},
//     {"Norm Range", {"varB1", "varB2"}},
//     {"Residual Bkg Fit", {"varC1", "varC2"}},
//     {"Fit Parameters", {"varD1", "varD2", "varD3", "varD4"}},
//     {"DCA Track to PV", {"varTrA1", "varTrA2"}},
//     {"TPC PID", {"varTrB1", "varTrB2"}},
//     {"TPC Crossed Rows", {"varTrC1", "varTrC2"}},
//     {"TPC Rows/Clusters", {"varTrD1", "varTrD2"}},
//     {"Cosine PA", {"varToA1", "varToA2"}},
//     {"Transverse Radius", {"varToB1", "varToB2"}},
//     {"DCA V0 Daughters", {"varToC1", "varToC2"}},
//     {"Lifetime", {"varToD1", "varToD2"}},
//     {"Competing V0 Rejection", {"varToE1", "varToE2"}},
//     {"Ks Mass Window", {"varToF1", "varToF2"}}};
