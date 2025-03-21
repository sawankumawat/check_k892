#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TH1F.h"
#include "TCanvas.h"
#include "../src/style.h"

using namespace std;

void systematics_plot()
{
    gStyle->SetOptStat(0);
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937/fits/4rBw_fits/backup2/";
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937/fits/4rBw_fits/";
    std::vector<std::string> filenames = {"sys_signal.txt", "sys_PID.txt", "sys_TrkSel.txt", "sys_TopSel.txt"};
    std::vector<std::string> binLabels = {"Signal extraction", "PID", "Track Sel.", "Topological Sel."};
    
    std::vector<double> massUncertainties;
    std::vector<double> widthUncertainties;

    for (const auto &filename : filenames)
    {
        std::ifstream file(path + filename);
        if (!file.is_open())
        {
            std::cerr << "Error opening file " << filename << std::endl;
            return;
        }
        // each file has simple mass and width uncertainties in two rows
        double massUncertainty, widthUncertainty;
        file >> massUncertainty >> widthUncertainty;
        cout << "Mass uncertainty from " << filename << " is " << massUncertainty << endl;
        cout << "Width uncertainty from " << filename << " is " << widthUncertainty << endl;
        massUncertainties.push_back(massUncertainty);
        widthUncertainties.push_back(widthUncertainty);
        file.close();
    }

    // calculate the total uncertainty using quadrature sum in both mass and width
    double totalMassUncertainty{0}, totalWidthUncertainty{0};
    for (const auto &massUncertainty : massUncertainties)
    {
        totalMassUncertainty += pow(massUncertainty, 2);
    }
    totalMassUncertainty = sqrt(totalMassUncertainty);
    cout << "Total Quadrature Sum Uncertainty (Mass): " << totalMassUncertainty << endl;
    for (const auto &widthUncertainty : widthUncertainties)
    {
        totalWidthUncertainty += pow(widthUncertainty, 2);
    }
    totalWidthUncertainty = sqrt(totalWidthUncertainty);
    cout << "Total Quadrature Sum Uncertainty (Width): " << totalWidthUncertainty << endl;

    TH1F *hMass = new TH1F("hMass", "Mass Uncertainties", 4, 0.5, 4.5);
    TH1F *hMassTotal = new TH1F("hMassTotal", "Total Mass Uncertainty", 1, 0.5, 4.5);
    TH1F *hWidth = new TH1F("hWidth", "Width Uncertainties", 4, 0.5, 4.5);
    TH1F *hWidthTotal = new TH1F("hWidthTotal", "Total Width Uncertainty", 1, 0.5, 4.5);
    hMassTotal->SetBinContent(1, totalMassUncertainty);
    hWidthTotal->SetBinContent(1, totalWidthUncertainty);

    for (int i = 0; i < 4; ++i)
    {
        hMass->SetBinContent(i + 1, massUncertainties[i]);
        hWidth->SetBinContent(i + 1, widthUncertainties[i]);
        hMass->GetXaxis()->SetBinLabel(i + 1, binLabels[i].c_str());
        hWidth->GetXaxis()->SetBinLabel(i + 1, binLabels[i].c_str());
    }

    // rotate the labels by 45 degrees
    //  hMass->LabelsOption("v"); // rotates by 90 degrees
    //  hWidth->LabelsOption("v");
    hMass->LabelsOption("d", "X"); // rotates by 45 degrees
    hWidth->LabelsOption("d", "X");
    // make the x axis labels left aligned

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.SetTextFont(42);
   
    TCanvas *c1 = new TCanvas("c1", "Relative Uncertainty Mass", 720, 720);
    SetCanvasStyle(c1, 0.18, 0.12, 0.05, 0.14);
    SetHistoQA(hMass);
    hMass->GetYaxis()->SetRangeUser(0, 0.0155);
    hMass->GetYaxis()->SetTitle("Fractional Uncertainty");
    hMass->GetYaxis()->SetTitleOffset(1.9);
    hMass->Draw();
    hMassTotal->SetLineColor(4);
    hMassTotal->SetLineWidth(3);
    hMassTotal->Draw("same");
    lat.DrawLatex(0.25, 0.88, "ALICE");
    lat.DrawLatex(0.25, 0.82, "pp #sqrt{#it{s}} = 13.6 TeV");
    lat.DrawLatex(0.25, 0.76, "FT0M (0-100%), |#it{y}|<0.5");
    lat.DrawLatex(0.25, 0.70, "3.0 < #it{p}_{T} < 5.0 GeV/#it{c}");
    lat.DrawLatex(0.25, 0.64, "f_{0}(1710) Mass");
    c1->SaveAs("/home/sawan/Videos/fract_uncert_mass_highpt1.png");

    TCanvas *c2 = new TCanvas("c2", "Relative Uncertainty Width", 720, 720);
    SetCanvasStyle(c2, 0.18, 0.12, 0.05, 0.14);
    SetHistoQA(hWidth);
    hWidth->GetYaxis()->SetRangeUser(0, 0.5);
    hWidth->GetYaxis()->SetTitle("Fractional Uncertainty");
    hWidth->GetYaxis()->SetTitleOffset(1.6);
    hWidth->Draw();
    hWidthTotal->SetLineColor(4);
    hWidthTotal->SetLineWidth(3);
    hWidthTotal->Draw("same");
    lat.DrawLatex(0.25, 0.88, "ALICE");
    lat.DrawLatex(0.25, 0.82, "pp #sqrt{#it{s}} = 13.6 TeV");
    lat.DrawLatex(0.25, 0.76, "FT0M (0-100%), |#it{y}|<0.5");
    lat.DrawLatex(0.25, 0.70, "3.0 < #it{p}_{T} < 5.0 GeV/#it{c}");
    lat.DrawLatex(0.25, 0.64, "f_{0}(1710) Width");
    c2->SaveAs("/home/sawan/Videos/fract_uncert_width_highpt1.png");
}
