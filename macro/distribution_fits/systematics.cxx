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

// systematic studies (signal extraction) ****************************
// A. fit range: Default: 1.05-2.20 GeV/c^2, Variation1: 1.02-2.20 GeV/c^2, Variation2: 1.05-2.30 GeV/c^2, Variation3: 1.08-2.15 GeV/c^2, Variation4: 1.02-2.30 GeV/c^2
// B. Norm range: Default: 2.50-2.60 GeV/c^2, Variation1: 2.40-2.50 GeV/c^2, Variation2: 2.60-2.70 GeV/c^2
// C. Residual bkg Fit function: Default: 4rBW with mass dependent width + modified Boltzmann, 4rBW with mass dependent width + Expol1, 4rBW with mass dependent width + Boltzmann
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

using namespace std;

void systematics()
{
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/351471/KsKs_Channel/higher-mass-resonances_3sigmaKs/fits/4rBw_fits/systematics_onlywidthfixed/";
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937/fits/4rBw_fits/backup/";
    ifstream file;

    int colors[] = {4, 6, 28, 46};

    // Mapping variations to sources
    int variations[] = {1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5};
    int variationsNumber[] = {1, 4, 2, 2, 4}; // Number of variations for each source

    vector<string> names = {"default", "varA1", "varA2", "varA3", "varA4", "varB1", "varB2",
                            "varC1", "varC2", "varD1", "varD2", "varD3", "varD4"};

    // Vectors to store fit parameters and their uncertainties
    vector<vector<double>> allparams_mass(5), allparams_width(5);
    vector<vector<double>> allparams_mass_unc(5), allparams_width_unc(5); // To store uncertainties

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

        // Ensure at least two values were read
        if (values.size() >= 2)
        {
            cout << "Last two values from " << names[ivars] << ": "
                 << values[values.size() - 2] << " ± " << uncertainties[values.size() - 2] << " & "
                 << values[values.size() - 1] << " ± " << uncertainties[values.size() - 1] << endl;

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

    // now lets plot the the variations
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c1, 0.18, 0.03, 0.05, 0.14);
    TH1F *hmass_var[4], *hwidth_var[4];
    for (int i = 0; i < 4; i++)
    {
        hmass_var[i] = new TH1F(Form("hmass_var%d", i), "", 1, 1, 10);
        hmass_var[i]->SetBinContent(1, relative_uncertainty_mass[i]);
        SetHistoQA(hmass_var[i]);
        hmass_var[i]->SetLineColor(colors[i]);
        hmass_var[i]->SetLineWidth(3);
        hmass_var[i]->SetLineStyle(2);
        hmass_var[i]->GetYaxis()->SetRangeUser(0, 0.0161);
        hmass_var[i]->GetYaxis()->SetTitleOffset(1.8);
        hmass_var[i]->GetYaxis()->SetTitle("Relative Uncertainty (Mass)");
        hmass_var[i]->GetXaxis()->SetTitle("p_{T} GeV/c");
        hmass_var[i]->Draw("same");
    }
    TH1F *hmass_total = new TH1F("hmass_total", "", 1, 1, 10);
    hmass_total->SetBinContent(1, total_uncertainty_mass);
    SetHistoQA(hmass_total);
    hmass_total->SetLineColor(1);
    hmass_total->SetLineWidth(3);
    // hmass_total->SetLineStyle(2);
    hmass_total->Draw("same");

    TLegend *l1 = new TLegend(0.25, 0.72, 0.55, 0.91);
    l1->SetFillStyle(0);
    l1->SetBorderSize(0);
    l1->SetTextFont(42);
    l1->SetTextSize(0.04);
    l1->AddEntry(hmass_var[0], "Fit range", "l");
    l1->AddEntry(hmass_var[1], "Normalization range", "l");
    l1->AddEntry(hmass_var[2], "Fit function", "l");
    l1->AddEntry(hmass_var[3], "Fit parameters", "l");
    l1->AddEntry(hmass_total, "Total", "l");
    l1->Draw("same");

    c1->SaveAs("/home/sawan/Music/systematic_mass.png");

    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.18, 0.03, 0.05, 0.14);
    for (int i = 0; i < 4; i++)
    {
        hwidth_var[i] = new TH1F(Form("hwidth_var%d", i), "", 1, 1, 10);
        hwidth_var[i]->SetBinContent(1, relative_uncertainty_width[i]);
        SetHistoQA(hwidth_var[i]);
        hwidth_var[i]->SetLineColor(colors[i]);
        hwidth_var[i]->SetLineWidth(3);
        hwidth_var[i]->SetLineStyle(2);
        hwidth_var[i]->GetYaxis()->SetRangeUser(0, 0.301);
        hwidth_var[i]->GetYaxis()->SetTitleOffset(1.8);
        hwidth_var[i]->GetYaxis()->SetTitle("Relative Uncertainty (Width)");
        hwidth_var[i]->GetXaxis()->SetTitle("p_{T} GeV/c");
        hwidth_var[i]->Draw("same");
    }
    TH1F *hwidth_total = new TH1F("hwidth_total", "", 1, 1, 10);
    hwidth_total->SetBinContent(1, total_uncertainty_width);
    SetHistoQA(hwidth_total);
    hwidth_total->SetLineColor(1);
    hwidth_total->SetLineWidth(3);
    // hwidth_total->SetLineStyle(2);
    hwidth_total->Draw("same");

    TLegend *l2 = new TLegend(0.25, 0.72, 0.55, 0.91);
    l2->SetFillStyle(0);
    l2->SetBorderSize(0);
    l2->SetTextFont(42);
    l2->SetTextSize(0.04);
    l2->AddEntry(hwidth_var[0], "Fit range", "l");
    l2->AddEntry(hwidth_var[1], "Normalization range", "l");
    l2->AddEntry(hwidth_var[2], "Fit function", "l");
    l2->AddEntry(hwidth_var[3], "Fit parameters", "l");
    l2->AddEntry(hwidth_total, "Total", "l");
    l2->Draw("same");

    c2->SaveAs("/home/sawan/Music/systematic_width.png");
}
