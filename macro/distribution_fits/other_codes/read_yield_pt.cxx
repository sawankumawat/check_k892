#include <iostream>
#include <utility>
#include "../src/style.h"
#include "../src/common_glue.h"
#include "TF1.h"
#include "TMath.h"
#include "../spectra/YieldMean.C"

using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void read_yield_pt()
{
    gStyle->SetOptStat(0);
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/";
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/WidthFree/mult_0-100/";
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/435450/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/";

    // Temporarily openning f1525 yield file for single BW fit
    ifstream yieldfile;
    yieldfile.open("/home/sawan/check_k892/macro/distribution_fits/f1525_yield_fixWidth3.txt");
    if (!yieldfile.is_open())
    {
        cout << "Error opening f1525 yield file" << endl;
        return;
    }
    vector<double> f1525yieldSingleBW, f1525yieldSingleBWErr;
    double yield1525singleBW, yield1525singleBWErr;
    string pm;

    while (yieldfile >> yield1525singleBW >> pm >> yield1525singleBWErr)
    {
        f1525yieldSingleBW.push_back(yield1525singleBW);
        f1525yieldSingleBWErr.push_back(yield1525singleBWErr);
        // cout<<"Read f1525 yield: "<< yield1525singleBW << " +/- " << yield1525singleBWErr << endl;
    }

    TString outputPath = path + "/Spectra";
    if (gSystem->mkdir(outputPath, kTRUE))
    {
        std::cout << "Folder " << outputPath << " created successfully." << std::endl;
    }

    float ptBins[] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0};             // 2022 dataset
    float ptBins2[] = {0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 16.0}; // 2022 dataset
    // float ptBins[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};             // 2023 dataset
    // float ptBins2[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0}; // 2023 dataset

    int nBins = sizeof(ptBins) / sizeof(ptBins[0]) - 1;
    int nBins2 = sizeof(ptBins2) / sizeof(ptBins2[0]) - 1;

    // efficiency file
    TFile *feff = new TFile("/home/sawan/check_k892/mc/LHC24l1/463655.root", "read");
    if (feff->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    string histpath = "higher-mass-resonances/hMChists";
    string histpathf2 = "higher-mass-resonances_f21525/hMChists";
    TFile *fOutput = new TFile(Form("%s/spectra.root", outputPath.Data()), "RECREATE");

    THnSparseD *GenpTf0 = (THnSparseD *)feff->Get(Form("%s/Genf17102", histpath.c_str())); // axis: multiplicity, pt, helicity angle
    THnSparseD *GenpTf2 = (THnSparseD *)feff->Get(Form("%s/Genf17102", histpathf2.c_str()));
    THnSparseD *recpt1f0 = (THnSparseD *)feff->Get(Form("%s/Recf1710_pt2", histpath.c_str())); // axis: multiplicity, pt, mass, helicity angle
    THnSparseD *recpt1f2 = (THnSparseD *)feff->Get(Form("%s/Recf1710_pt2", histpathf2.c_str()));

    if (GenpTf0 == nullptr || recpt1f0 == nullptr || GenpTf2 == nullptr || recpt1f2 == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    TH1F *hefficiencyf0 = new TH1F("hefficiencyf0", "Efficiency vs pT", nBins, ptBins);
    TH1F *hefficiencyf2 = new TH1F("hefficiencyf2", "Efficiency vs pT for f21525", nBins, ptBins);

    TH1D *hgenptf0 = GenpTf0->Projection(1);  // project on pt axis
    TH1D *hrecptf0 = recpt1f0->Projection(1); // project on pt axis
    TH1D *hgenptf2 = GenpTf2->Projection(1);  // project on pt axis
    TH1D *hrecptf2 = recpt1f2->Projection(1); // project on pt axis

    for (int i = 0; i < nBins; i++)
    {
        // get bin content accroding to cosTheta bins and error according to bayesian method
        int lowpt = hgenptf0->GetXaxis()->FindBin(ptBins[i] + 0.01);
        int highpt = hgenptf0->GetXaxis()->FindBin(ptBins[i + 1] - 0.01);
        double genYieldf0 = hgenptf0->Integral(lowpt, highpt);
        double recYieldf0 = hrecptf0->Integral(lowpt, highpt);
        double recYieldError = TMath::Sqrt(((recYieldf0 + 1) / (genYieldf0 + 2)) * ((recYieldf0 + 2) / (genYieldf0 + 3) - (recYieldf0 + 1) / (genYieldf0 + 2)));
        if (genYieldf0 > 0)
        {
            hefficiencyf0->SetBinContent(i + 1, recYieldf0 / genYieldf0);
            hefficiencyf0->SetBinError(i + 1, recYieldError);
        }

        double efficiencyf2 = hgenptf2->Integral(lowpt, highpt);
        double recYieldf2 = hrecptf2->Integral(lowpt, highpt);
        double recYieldErrorf2 = TMath::Sqrt(((recYieldf2 + 1) / (efficiencyf2 + 2)) * ((recYieldf2 + 2) / (efficiencyf2 + 3) - (recYieldf2 + 1) / (efficiencyf2 + 2)));
        if (efficiencyf2 > 0)
        {
            hefficiencyf2->SetBinContent(i + 1, recYieldf2 / efficiencyf2);
            hefficiencyf2->SetBinError(i + 1, recYieldErrorf2);
        }
    }

    TH1D *hYield1710Raw = new TH1D("hYield1710Raw", "Yield vs cosTheta", nBins2, ptBins2);
    TH1D *hYield1525Raw = new TH1D("hYield1525Raw", "Yield vs cosTheta", nBins2, ptBins2);
    TH1D *hYield1270Raw = new TH1D("hYield1270Raw", "Yield vs cosTheta", nBins2, ptBins2);
    TH1D *hYield1320Raw = new TH1D("hYield1320Raw", "Yield vs cosTheta", nBins2, ptBins2);

    TH1D *hYield1710Corrected = new TH1D("hYield1710Corrected", "Yield vs cosTheta", nBins2, ptBins2);
    TH1D *hYield1525Corrected = new TH1D("hYield1525Corrected", "Yield vs cosTheta", nBins2, ptBins2);
    TH1D *hYield1270Corrected = new TH1D("hYield1270Corrected", "Yield vs cosTheta", nBins2, ptBins2);
    TH1D *hYield1320Corrected = new TH1D("hYield1320Corrected", "Yield vs cosTheta", nBins2, ptBins2);

    TH1D *hMass1710 = new TH1D("hMass1710", "Mass vs cosTheta", nBins2, ptBins2);
    TH1D *hWidth1710 = new TH1D("hWidth1710", "Width vs cosTheta", nBins2, ptBins2);
    TH1D *hMass1525 = new TH1D("hMass1525", "Mass vs cosTheta", nBins2, ptBins2);
    TH1D *hMass1270 = new TH1D("hMass1270", "Mass vs cosTheta", nBins2, ptBins2);
    TH1D *hMass1320 = new TH1D("hMass1320", "Mass vs cosTheta", nBins2, ptBins2);

    TH1F *hSignificance1710 = new TH1F("hSignificance1710", "Significance vs pT", nBins2, ptBins2);
    TH1F *hSignificance1525 = new TH1F("hSignificance1525", "Significance vs pT", nBins2, ptBins2);

    for (int i = 0; i < nBins; i++)
    {
        if (nBins != hefficiencyf0->GetNbinsX())
        {
            cout << "number of costheta bins " << nBins << ", number of efficiency bins " << hefficiencyf0->GetNbinsX() << endl;
            cout << "Bins mismatch " << endl;
            return;
        }

        ifstream infile;
        infile.open(path + Form("fit_params_pT_%.1f-%.1fdefault.txt", ptBins[i], ptBins[i + 1]));
        // cout << "path is " << path + "fit_params_cosTheta_" + to_string(ptBins[i]) + "_" + to_string(ptBins[i + 1]) + ".txt" << endl;

        std::string line;

        // Common fit info
        double significance = 0, statSignificance = 0, chi2ndf = 0;

        // Parameters for f1710
        double yield1 = 0, yield1_err = 0;
        double yield1Bin = 0, yield1Bin_err = 0;
        double mass1 = 0, mass1_err = 0;
        double width1 = 0, width1_err = 0;
        double bkgYield1 = 0, bkgYield1_err = 0, bkgYield2 = 0, bkgYield2_err = 0, bkgYield3 = 0, bkgYield3_err = 0;
        double signficance1{0}, signficance1_err{0};

        // Parameters for f1525
        double yield2 = 0, yield2_err = 0;
        double yield2Bin = 0, yield2Bin_err = 0;
        double mass2 = 0, mass2_err = 0;
        double width2 = 0, width2_err = 0;
        double signficance2{0}, signficance2_err{0};

        // Parameters for f1270
        double yield3 = 0, yield3_err = 0;
        double yield3Bin = 0, yield3Bin_err = 0;
        double mass3 = 0, mass3_err = 0;
        double width3 = 0, width3_err = 0;

        // Parameters for a1320
        double yield4 = 0, yield4_err = 0;
        double yield4Bin = 0, yield4Bin_err = 0;
        double mass4 = 0, mass4_err = 0;
        double width4 = 0, width4_err = 0;

        int paramIndex = 0;
        bool readingF1710 = false;
        bool readingF1525 = false;
        bool readingF1270 = false;
        bool readingA1320 = false;
        bool readingBkg = false;

        if (infile.is_open())
        {
            while (std::getline(infile, line))
            {
                std::istringstream iss(line);

                if (line.find("StatSignificance") != std::string::npos)
                    iss.ignore(100, ' '), iss >> statSignificance;

                else if (line.find("Significance") != std::string::npos && line.find("Stat") == std::string::npos)
                    iss.ignore(100, ' '), iss >> significance;

                else if (line.find("Chi2NDF") != std::string::npos)
                    iss.ignore(100, ' '), iss >> chi2ndf;

                else if (line.find("f1710") != std::string::npos)
                {
                    readingF1710 = true;
                    readingF1525 = false;
                    readingF1270 = false;
                    readingA1320 = false;
                    readingBkg = false;
                    paramIndex = 0;
                }

                else if (line.find("f1525") != std::string::npos)
                {
                    readingF1525 = true;
                    readingF1710 = false;
                    readingF1270 = false;
                    readingA1320 = false;
                    readingBkg = false;
                    paramIndex = 0;
                }
                else if (line.find("f1270") != std::string::npos)
                {
                    readingF1270 = true;
                    readingA1320 = false;
                    readingF1710 = false;
                    readingF1525 = false;
                    readingBkg = false;
                    paramIndex = 0;
                }
                else if (line.find("a1320") != std::string::npos)
                {
                    readingA1320 = true;
                    readingF1270 = false;
                    readingF1710 = false;
                    readingF1525 = false;
                    readingBkg = false;
                    paramIndex = 0;
                }
                else if (line.find("Norm range") != std::string::npos || line.find("Fit range") != std::string::npos)
                {
                    // Skip norm and fit range lines
                    continue;
                }
                else if (line.find("Yield region") != std::string::npos)
                {
                    // Skip yield region lines
                    continue;
                }

                else if (line.find("±") != std::string::npos)
                {
                    double val = 0, err = 0;
                    string pm;
                    std::istringstream(line) >> val >> pm >> err;

                    if (readingF1710)
                    {
                        switch (paramIndex++)
                        {
                        case 0:
                            yield1 = val;
                            yield1_err = err;
                            break;
                        case 1:
                            yield1Bin = val;
                            yield1Bin_err = err;
                            break;
                        // case 1:
                        case 2:
                            mass1 = val;
                            mass1_err = err;
                            break;
                        // case 2:
                        case 3:
                            width1 = val;
                            width1_err = err;
                            break;
                        }
                    }
                    else if (readingF1525)
                    {
                        switch (paramIndex++)
                        {
                        case 0:
                            yield2 = val;
                            yield2_err = err;
                            break;
                        case 1:
                            yield2Bin = val;
                            yield2Bin_err = err;
                            break;
                        // case 1:
                        case 2:
                            mass2 = val;
                            mass2_err = err;
                            break;
                        // case 2:
                        case 3:
                            width2 = val;
                            width2_err = err;
                            break;
                        }
                    }
                    else if (readingF1270)
                    {
                        switch (paramIndex++)
                        {
                        case 0:
                            yield3 = val;
                            yield3_err = err;
                            break;
                        case 1:
                            yield3Bin = val;
                            yield3Bin_err = err;
                            break;
                        // case 1:
                        case 2:
                            mass3 = val;
                            mass3_err = err;
                            break;
                        // case 2:
                        case 3:
                            width3 = val;
                            width3_err = err;
                            break;
                        }
                    }
                    else if (readingA1320)
                    {
                        switch (paramIndex++)
                        {
                        case 0:
                            yield4 = val;
                            yield4_err = err;
                            break;
                        case 1:
                            yield4Bin = val;
                            yield4Bin_err = err;
                            break;
                        // case 1:
                        case 2:
                            mass4 = val;
                            mass4_err = err;
                            break;
                        case 3:
                            width4 = val;
                            width4_err = err;
                            break;
                        }
                    }
                    else if (readingBkg)
                    {
                        switch (paramIndex++)
                        {
                        case 0:
                            bkgYield1 = val;
                            bkgYield1_err = err;
                            break;
                        case 1:
                            bkgYield2 = val;
                            bkgYield2_err = err;
                            break;
                        case 2:
                            bkgYield3 = val;
                            bkgYield3_err = err;
                            break;
                        }
                    }
                }
            }
        }

        double eff1710 = hefficiencyf0->GetBinContent(i + 1);
        double eff1710_err = hefficiencyf0->GetBinError(i + 1);
        double eff1525 = hefficiencyf2->GetBinContent(i + 1);
        double eff1525_err = hefficiencyf2->GetBinError(i + 1);
        double oneUponTriggerEfficiency = 70.0 / 46.1; // inverse of trigger efficiency
        double BR_f0 = 0.1667 / 2;
        double BR_f2 = 0.438 / 2;

        double corrected_yield1710 = yield1 / (eff1710 * oneUponTriggerEfficiency * BR_f0);
        double corrected_yield1710_err = sqrt(pow(yield1_err / eff1710, 2) + pow(yield1 * eff1710_err / pow(eff1710, 2), 2)) / (oneUponTriggerEfficiency * BR_f0);
        double corrected_yield1710Bin = yield1Bin / (eff1710 * oneUponTriggerEfficiency * BR_f0);
        double corrected_yield1710Bin_err = sqrt(pow(yield1Bin_err / eff1710, 2) + pow(yield1Bin * eff1710_err / pow(eff1710, 2), 2)) / (oneUponTriggerEfficiency * BR_f0);

        // //================Modifying f21525 yield for single BW fit check=================//
        // yield2 = f1525yieldSingleBW[i];
        // yield2_err = f1525yieldSingleBWErr[i];
        // //===============================================================================//

        double corrected_yield1525 = yield2 / (eff1525 * oneUponTriggerEfficiency * BR_f2);
        double corrected_yield1525_err = sqrt(pow(yield2_err / eff1525, 2) + pow(yield2 * eff1525_err / pow(eff1525, 2), 2)) / (oneUponTriggerEfficiency * BR_f2);
        cout << "Bin " << i << ", Raw yield f1525: " << yield2 << ", efficiency: " << eff1525 << endl;
        double corrected_yield1525Bin = yield2Bin / (eff1525 * oneUponTriggerEfficiency * BR_f2);
        double corrected_yield1525Bin_err = sqrt(pow(yield2Bin_err / eff1525, 2) + pow(yield2Bin * eff1525_err / pow(eff1525, 2), 2)) / (oneUponTriggerEfficiency * BR_f2);

        double corrected_yield1270 = yield3 / eff1525;
        double corrected_yield1270_err = sqrt(pow(yield3_err / eff1525, 2) + pow(yield3 * eff1525_err / pow(eff1525, 2), 2));
        double corrected_yield1270Bin = yield3Bin / eff1525;
        double corrected_yield1270Bin_err = sqrt(pow(yield3Bin_err / eff1525, 2) + pow(yield3Bin * eff1525_err / pow(eff1525, 2), 2));

        double corrected_yield1320 = yield4 / eff1525;
        double corrected_yield1320_err = sqrt(pow(yield4_err / eff1525, 2) + pow(yield4 * eff1525_err / pow(eff1525, 2), 2));
        double corrected_yield1320Bin = yield4Bin / eff1525;
        double corrected_yield1320Bin_err = sqrt(pow(yield4Bin_err / eff1525, 2) + pow(yield4Bin * eff1525_err / pow(eff1525, 2), 2));

        hYield1710Raw->SetBinContent(i + 2, yield1);
        hYield1710Raw->SetBinError(i + 2, yield1_err);
        hYield1710Corrected->SetBinContent(i + 2, corrected_yield1710);
        hYield1710Corrected->SetBinError(i + 2, corrected_yield1710_err);
        hMass1710->SetBinContent(i + 2, mass1);
        hMass1710->SetBinError(i + 2, mass1_err);
        hWidth1710->SetBinContent(i + 2, width1);
        hWidth1710->SetBinError(i + 2, width1_err);

        hYield1525Raw->SetBinContent(i + 2, yield2);
        hYield1525Raw->SetBinError(i + 2, yield2_err);
        hYield1525Corrected->SetBinContent(i + 2, corrected_yield1525);
        hYield1525Corrected->SetBinError(i + 2, corrected_yield1525_err);
        // cout<<"Bin "<<i<<", corrected f1525 yield: "<<corrected_yield1525<<" +/- "<<corrected_yield1525_err<<endl;
        hMass1525->SetBinContent(i + 2, mass2);
        hMass1525->SetBinError(i + 2, mass2_err);

        hYield1270Raw->SetBinContent(i + 2, yield3);
        hYield1270Raw->SetBinError(i + 2, yield3_err);
        hYield1270Corrected->SetBinContent(i + 2, corrected_yield1270);
        hYield1270Corrected->SetBinError(i + 2, corrected_yield1270_err);
        hMass1270->SetBinContent(i + 2, mass3);
        hMass1270->SetBinError(i + 2, mass3_err);

        hYield1320Raw->SetBinContent(i + 2, yield4);
        hYield1320Raw->SetBinError(i + 2, yield4_err);
        hYield1320Corrected->SetBinContent(i + 2, corrected_yield1320);
        hYield1320Corrected->SetBinError(i + 2, corrected_yield1320_err);
        hMass1320->SetBinContent(i + 2, mass4);
        hMass1320->SetBinError(i + 2, mass4_err);
    }

    TCanvas *cMass1710 = new TCanvas("cMass1710", "Mass vs #it{p}_{T} for f_{0}(1710)", 720, 720);
    SetCanvasStyle(cMass1710, 0.18, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    SetHistoQA(hMass1710);
    hMass1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass1710->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass1710->GetYaxis()->SetTitleOffset(1.6);
    hMass1710->GetYaxis()->SetRangeUser(1.56, 1.89);
    hMass1710->SetMarkerStyle(21);
    hMass1710->Write();
    hMass1710->Draw("pe");
    // TLine *line1710Mass = new TLine(0, f1710Mass, 12, f1710Mass);
    // line1710Mass->SetLineStyle(2);
    // line1710Mass->SetLineColor(2);
    // line1710Mass->Draw();
    TBox *band1710Mass = new TBox(1, f1710Mass - f1710MassErr, 15, f1710Mass + f1710MassErr);
    band1710Mass->SetFillColorAlpha(kRed, 0.2); // shaded
    band1710Mass->SetLineColor(kRed);
    band1710Mass->SetLineWidth(1);
    band1710Mass->Draw("same");
    TLegend *leg1710Mass = new TLegend(0.55, 0.75, 0.85, 0.9);
    leg1710Mass->SetBorderSize(0);
    leg1710Mass->SetFillStyle(0);
    leg1710Mass->SetTextSize(0.035);
    leg1710Mass->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    leg1710Mass->AddEntry(hMass1710, "This analysis", "pe");
    // leg1710Mass->AddEntry(line1710Mass, "PDG value", "l");
    band1710Mass->SetLineWidth(0);
    leg1710Mass->AddEntry(band1710Mass, "World-average value", "f");
    leg1710Mass->Draw();
    cMass1710->SaveAs(outputPath + "/Mass1710.png");

    TCanvas *cWidth1710 = new TCanvas("cWidth1710", "Width vs #it{p}_{T} for f_{0}(1710)", 720, 720);
    SetCanvasStyle(cWidth1710, 0.18, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    SetHistoQA(hWidth1710);
    hWidth1710->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hWidth1710->GetYaxis()->SetTitle("Width (GeV/#it{c}^{2})");
    hWidth1710->GetYaxis()->SetTitleOffset(1.6);
    hWidth1710->GetYaxis()->SetRangeUser(0.05, 0.32);
    hWidth1710->SetMarkerStyle(21);
    hWidth1710->Write();
    hWidth1710->Draw("pe");
    TLine *line1710Width = new TLine(0, f1710Width, 12, f1710Width);
    line1710Width->SetLineStyle(2);
    line1710Width->SetLineColor(2);
    line1710Width->Draw();
    TLegend *leg1710Width = new TLegend(0.55, 0.75, 0.85, 0.9);
    leg1710Width->SetBorderSize(0);
    leg1710Width->SetFillStyle(0);
    leg1710Width->SetTextSize(0.035);
    leg1710Mass->Draw();
    cWidth1710->SaveAs(outputPath + "/Width1710.png");

    TCanvas *cMass1525 = new TCanvas("cMass1525", "Mass vs #it{p}_{T} for f_{2}(1525)", 720, 720);
    SetCanvasStyle(cMass1525, 0.18, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    SetHistoQA(hMass1525);
    hMass1525->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass1525->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass1525->GetYaxis()->SetTitleOffset(1.6);
    hMass1525->GetYaxis()->SetRangeUser(1.49, 1.55);
    hMass1525->SetMarkerStyle(20);
    hMass1525->Write();
    hMass1525->Draw("pe");
    // TLine *line1525Mass = new TLine(0, f1525Mass, 12, f1525Mass);
    // line1525Mass->SetLineStyle(2);
    // line1525Mass->SetLineColor(2);
    // line1525Mass->Draw();
    TBox *band1525Mass = new TBox(1, f1525Mass - f1525MassErr, 15, f1525Mass + f1525MassErr);
    band1525Mass->SetFillColorAlpha(kRed, 0.2); // shaded
    band1525Mass->SetLineColor(kRed);
    band1525Mass->SetLineWidth(1);
    band1525Mass->Draw("same");
    leg1710Mass->Draw();
    cMass1525->SaveAs(outputPath + "/Mass1525.png");

    TCanvas *cMass1270 = new TCanvas("cMass1270", "Mass vs #it{p}_{T} for f_{2}(1270)", 720, 720);
    SetCanvasStyle(cMass1270, 0.18, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    SetHistoQA(hMass1270);
    hMass1270->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass1270->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass1270->GetYaxis()->SetTitleOffset(1.6);
    hMass1270->GetYaxis()->SetRangeUser(1.1, 1.4);
    hMass1270->SetMarkerStyle(20);
    hMass1270->Write();
    hMass1270->Draw("pe");
    TLine *line1270Mass = new TLine(0, f1270Mass, 12, f1270Mass);
    line1270Mass->SetLineStyle(2);
    line1270Mass->SetLineColor(2);
    line1270Mass->Draw();
    leg1710Mass->Draw();
    cMass1270->SaveAs(outputPath + "/Mass1270.png");

    TCanvas *cMass1320 = new TCanvas("cMass1320", "Mass vs #it{p}_{T} for a_{2}(1320)", 720, 720);
    SetCanvasStyle(cMass1320, 0.18, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    SetHistoQA(hMass1320);
    hMass1320->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass1320->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass1320->GetYaxis()->SetTitleOffset(1.6);
    hMass1320->GetYaxis()->SetRangeUser(1.15, 1.45);
    hMass1320->SetMarkerStyle(20);
    hMass1320->Write();
    hMass1320->Draw("pe");
    TLine *line1320Mass = new TLine(0, a1320Mass, 12, a1320Mass);
    line1320Mass->SetLineStyle(2);
    line1320Mass->SetLineColor(2);
    line1320Mass->Draw();
    leg1710Mass->Draw();
    cMass1320->SaveAs(outputPath + "/Mass1320.png");

    TCanvas *cYieldCorrectedf1525 = new TCanvas("cYieldCorrectedf1525", "Yield vs #it{p}_{T} for f_{2}(1525)", 720, 720);
    SetCanvasStyle(cYieldCorrectedf1525, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(hYield1525Corrected);
    hYield1525Corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hYield1525Corrected->GetYaxis()->SetTitle("1/N_{ev} * dN/(d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hYield1525Corrected->GetYaxis()->SetTitleOffset(1.6);
    // hYield1525Corrected->GetYaxis()->SetRangeUser(0.0, 460);
    hYield1525Corrected->SetMaximum(hYield1525Corrected->GetMaximum() * 1.5);
    hYield1525Corrected->SetMarkerStyle(20);
    // hYield1525Corrected->SetMinimum(-2e-7);
    hYield1525Corrected->Write();
    hYield1525Corrected->Draw("pe");
    cYieldCorrectedf1525->SaveAs(outputPath + "/CorrectedYieldf2.png");

    TCanvas *cYieldCorrected1710 = new TCanvas("cYieldCorrected1710", "Yield vs #it{p}_{T} for f_{0}(1710)", 720, 720);
    SetCanvasStyle(cYieldCorrected1710, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(hYield1710Corrected);
    hYield1710Corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hYield1710Corrected->GetYaxis()->SetTitle("1/N_{ev} * dN/(d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hYield1710Corrected->GetYaxis()->SetTitleOffset(1.6);
    // hYield1710Corrected->GetYaxis()->SetRangeUser(0.0, 460);
    hYield1710Corrected->SetMaximum(hYield1710Corrected->GetMaximum() * 1.5);
    hYield1710Corrected->SetMarkerStyle(21);
    // hYield1710Corrected->SetMinimum(-2e-7);
    hYield1710Corrected->Write();
    hYield1710Corrected->Draw("pe");
    cYieldCorrected1710->SaveAs(outputPath + "/CorrectedYieldf0.png");

    TCanvas *cYieldRatio = new TCanvas("cYieldRatio", "Yield ratio vs #it{p}_{T} for f_{0}(1710)/f_{2}'(1525)", 720, 720);
    SetCanvasStyle(cYieldRatio, 0.18, 0.03, 0.05, 0.14);
    // gPad->SetLogy();
    TH1D *hYieldRatio = (TH1D *)hYield1710Corrected->Clone("hYieldRatio");
    hYieldRatio->Divide(hYield1525Corrected);
    SetHistoQA(hYieldRatio);
    hYieldRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hYieldRatio->GetYaxis()->SetTitle("f_{0}(1710)/f_{2}'(1525)");
    hYieldRatio->GetYaxis()->SetTitleOffset(1.6);
    // hYieldRatio->GetYaxis()->SetRangeUser(0.0, 460);
    hYieldRatio->SetMaximum(2.1);
    hYieldRatio->SetMarkerStyle(21);
    // hYieldRatio->SetMinimum(-2e-7);
    hYieldRatio->Write();
    hYieldRatio->Draw("pe");
    cYieldRatio->SaveAs(outputPath + "/YieldRatiof0f2.png");

    TCanvas *cEfficiencyf0f2 = new TCanvas("cEfficiencyf0f2", "Efficiency vs #it{p}_{T} for f_{0}(1710) and f_{2}(1525)", 720, 720);
    SetCanvasStyle(cEfficiencyf0f2, 0.18, 0.03, 0.05, 0.14);
    TH1F *hefficiencyClone = new TH1F("", "", nBins2, ptBins2);
    SetHistoQA(hefficiencyClone);
    hefficiencyClone->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hefficiencyClone->GetYaxis()->SetTitle("Acceptance x Efficiency");
    hefficiencyClone->GetYaxis()->SetTitleOffset(1.6);
    hefficiencyClone->SetMaximum(hefficiencyf0->GetMaximum() * 2.0);
    hefficiencyClone->Draw();
    // gPad->SetLogy();
    SetHistoQA(hefficiencyf0);
    hefficiencyf0->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hefficiencyf0->GetYaxis()->SetTitle("Acceptance x Efficiency");
    hefficiencyf0->GetYaxis()->SetTitleOffset(1.6);
    hefficiencyf0->SetMaximum(hefficiencyf0->GetMaximum() * 1.5);
    hefficiencyf0->SetMarkerStyle(21);
    hefficiencyf0->SetLineColor(kRed);
    hefficiencyf0->SetMarkerColor(kRed);
    // hefficiencyf0->SetMinimum(-2e-7);
    hefficiencyf0->Write();
    hefficiencyf0->Draw("pe same");
    SetHistoQA(hefficiencyf2);
    hefficiencyf2->SetMarkerStyle(20);
    hefficiencyf2->SetLineColor(kBlue);
    hefficiencyf2->SetMarkerColor(kBlue);
    hefficiencyf2->Write();
    hefficiencyf2->Draw("pe same");
    TLegend *legEff = new TLegend(0.55, 0.75, 0.85, 0.9);
    legEff->SetBorderSize(0);
    legEff->SetFillStyle(0);
    legEff->SetTextSize(0.035);
    legEff->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    legEff->AddEntry(hefficiencyf0, "f_{0}(1710)", "pe");
    legEff->AddEntry(hefficiencyf2, "f_{2}'(1525)", "pe");
    legEff->Draw();
    cEfficiencyf0f2->SaveAs(outputPath + "/Efficiencyf0f2.png");

    TCanvas *cRawYieldf2 = new TCanvas("cRawYieldf2", "Raw #it{p}_{T} distribution for f_{2}(1525)", 720, 720);
    SetCanvasStyle(cRawYieldf2, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(hYield1525Raw);
    hYield1525Raw->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hYield1525Raw->GetYaxis()->SetTitle("1/N_{ev} * dN/(d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hYield1525Raw->GetYaxis()->SetTitleOffset(1.6);
    hYield1525Raw->SetMaximum(hYield1525Raw->GetMaximum() * 1.5);
    // hYield1525Raw->SetMarkerStyle(21);
    // hYield1525Raw->SetMinimum(-2e-7);
    hYield1525Raw->Write();
    hYield1525Raw->Draw("pe");
    cRawYieldf2->SaveAs(outputPath + "/RawYieldf2.png");

    TCanvas *cRawYieldf0 = new TCanvas("cRawYieldf0", "Raw #it{p}_{T} distribution for f_{0}(1710)", 720, 720);
    SetCanvasStyle(cRawYieldf0, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(hYield1710Raw);
    hYield1710Raw->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hYield1710Raw->GetYaxis()->SetTitle("1/N_{ev} * dN/(d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hYield1710Raw->GetYaxis()->SetTitleOffset(1.6);
    hYield1710Raw->SetMaximum(hYield1710Raw->GetMaximum() * 1.5);
    // hYield1710Raw->SetMinimum(-2e-7);
    hYield1710Raw->Write();
    hYield1710Raw->Draw("pe same");
    cRawYieldf0->SaveAs(outputPath + "/RawYieldf0.png");

    TH1F *h1 = (TH1F *)hYield1525Corrected->Clone("h1");
    TH1F *h2 = (TH1F *)hYield1525Corrected->Clone("h2");

    for (int i = 1; i <= h2->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr = (0.1 * h2->GetBinContent(i));
        h2->SetBinError(i, systemerr);
    }
    /*************meanpT*****************byresonance*******************package*************************/
    Double_t min = 0.0;
    Double_t max = 15.0;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "REBMS0+";
    TString logfilename = "log.root";
    Double_t minfit = 1.0;
    Double_t maxfit = 10.0;

    TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
    fitFcn->SetParameter(0, 5.0);
    // fitFcn->SetParameter(1, 0.05);
    fitFcn->SetParameter(1, 0.5);
    fitFcn->FixParameter(2, 1.525);
    fitFcn->SetParameter(3, 0.35);
    fitFcn->SetParNames("n", "dn/dy", "mass", "T");

    TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    cout << "Yield dN/dy = " << hout->GetBinContent(1) << " +- " << hout->GetBinContent(2) << endl;
    cout << "Mean pT = " << hout->GetBinContent(5) << " +- " << hout->GetBinContent(6) << endl;

    TF1 *fitFcn2 = new TF1("fitfunc2", FuncLavy, 0.0, 15.0, 4);
    fitFcn2->SetParameter(0, 5.0);
    // fitFcn2->SetParameter(1, 0.05);
    fitFcn2->SetParameter(1, 0.5);
    fitFcn2->FixParameter(2, 1.710);
    fitFcn2->SetParameter(3, 0.35);
    fitFcn2->SetParNames("n", "dn/dy", "mass", "T");

    TH1 *hout2 = YieldMean(hYield1710Corrected, hYield1710Corrected, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    cout << "Yield dN/dy of f0(1710) = " << hout2->GetBinContent(1) << " +- " << hout2->GetBinContent(2) << endl;
    cout << "Mean pT of f0(1710) = " << hout2->GetBinContent(5) << " +- " << hout2->GetBinContent(6) << endl;

    TCanvas *cCorrectedf2Fit = new TCanvas("cCorrectedf2Fit", "Corrected #it{p}_{T} distribution with fit", 720, 720);
    SetCanvasStyle(cCorrectedf2Fit, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(h1);
    h1->SetMaximum(h1->GetMaximum() * 12);
    h1->SetMinimum(9e-9);
    h1->Draw("pe");
    fitFcn->SetLineWidth(2);
    fitFcn->SetLineStyle(2);
    fitFcn->Draw("l same");
    TLegend *legend2 = new TLegend(0.55, 0.75, 0.85, 0.9);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(0.035);
    legend2->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    legend2->AddEntry(h1, "f_{2}'(1525)", "pe");
    legend2->AddEntry(fitFcn, "Levy fit", "l");
    legend2->Draw();
    cCorrectedf2Fit->SaveAs(outputPath + "/LevyFitf2.png");

    TCanvas *cCorrectedf0Fit = new TCanvas("cCorrectedf0Fit", "Corrected #it{p}_{T} distribution with fit", 720, 720);
    SetCanvasStyle(cCorrectedf0Fit, 0.18, 0.03, 0.05, 0.14);
    gPad->SetLogy();
    SetHistoQA(hYield1710Corrected);
    hYield1710Corrected->SetMaximum(hYield1710Corrected->GetMaximum() * 12);
    hYield1710Corrected->SetMinimum(9e-8);
    hYield1710Corrected->Draw("pe");
    fitFcn2->SetLineWidth(2);
    fitFcn2->SetLineStyle(2);
    fitFcn2->Draw("l same");
    TLegend *legend3 = new TLegend(0.55, 0.75, 0.85, 0.9);
    legend3->SetBorderSize(0);
    legend3->SetFillStyle(0);
    legend3->SetTextSize(0.035);
    legend3->SetHeader("pp #sqrt{#it{s}} = 13.6 TeV");
    legend3->AddEntry(hYield1710Corrected, "f_{0}(1710)", "pe");
    legend3->AddEntry(fitFcn2, "Levy fit", "l");
    legend3->Draw();
    cCorrectedf0Fit->SaveAs(outputPath + "/LevyFitf0.png");

    TFile *flightFlavourHadrons = new TFile("../spectra/LightFlavourHadronsProduction.root", "read");
    if (flightFlavourHadrons->IsZombie())
    {
        cout << "Error opening light flavour hadrons production file" << endl;
        return;
    }
    // There are graphs of <pT> in the file in table 26 to 34 for different particles (pion, kaon, K0s, K*(892), phi, proton, Lambda, sigma, omega)
    //(π+/π−,K+/K−,KS0​,K∗(892),ϕ(1020),p/pˉ​,Λ/Λˉ,Σ+/Σ−,Ω−/Ωˉ+)

    int totalParticles = 9;
    string particles[9] = {"Pion", "Kaon", "K0s", "Kstar", "Phi", "Proton", "Lambda", "Xi", "Omega"};
    string particlesLatex[9] = {"#pi", "K", "K^{0}_{S}", "K^{*0}", "#phi", "p", "#Lambda", "#Xi^{-}", "#Omega"};
    double particleMass[9] = {0.13957, 0.49367, 0.49761, 0.89166, 1.01946, 0.93827, 1.11568, 1.3217, 1.67245}; // in GeV/c2
    int colors[9] = {kBlack, kBlue, kGreen + 2, kOrange + 7, kViolet + 7, kCyan + 2, kMagenta + 2, kGray + 2, kPink + 2};
    int markers[9] = {21, 22, 23, 33, 34, 43, 45, 29, 39};
    vector<vector<float>> dNdyvalues_13TeV = {
        {4.775, 0.001, 0.243, 1.0},  // 4.775 ± 0.001 ± 0.243
        {6.205, 0.004, 0.303, 1e-1}, // (6.205 ± 0.004 ± 0.303) × 10⁻¹
        {3.192, 0.004, 0.111, 1e-1}, // (3.192 ± 0.004 ± 0.111) × 10⁻¹
        {2.098, 0.016, 0.200, 1e-1}, // (2.098 ± 0.016 ± 0.200) × 10⁻¹
        {2.750, 0.002, 0.188, 1e-1}, // (2.750 ± 0.002 ± 0.188) × 10⁻¹
        {3.734, 0.040, 0.213, 1e-2}, // (3.734 ± 0.040 ± 0.213) × 10⁻²
        {1.807, 0.005, 0.102, 1e-1}, // (1.807 ± 0.005 ± 0.102) × 10⁻¹
        {1.980, 0.012, 0.082, 1e-2}, // (1.980 ± 0.012 ± 0.082) × 10⁻²
        {1.846, 0.046, 0.122, 1e-3}  // (1.846 ± 0.046 ± 0.122) × 10⁻³
    };

    double meanPtAt13TeV[9];
    double meanPtAt13TeV_err[9];
    TGraphErrors *gMeanPt[9];
    TGraphErrors *gMeanPtvsMassMesons = new TGraphErrors();
    TGraphErrors *gMeanPtvsMassBaryons = new TGraphErrors();
    for (int i = 0; i < totalParticles; i++)
    {
        gMeanPt[i] = (TGraphErrors *)flightFlavourHadrons->Get(Form("Table %d/Graph1D_y1", 26 + i));
        if (gMeanPt[i] == nullptr)
        {
            cout << "Error reading graph for " << particlesLatex[i] << endl;
            return;
        }
        gMeanPt[i]->SetMarkerColor(colors[i]);
        gMeanPt[i]->SetLineColor(colors[i]);
        gMeanPt[i]->SetMarkerStyle(markers[i]);
        gMeanPt[i]->SetMarkerSize(1);

        // Now from graph store the value of mean pT at 13 TeV
        int nPoints = gMeanPt[i]->GetN();
        double *x = gMeanPt[i]->GetX();
        double *y = gMeanPt[i]->GetY();
        for (int j = 0; j < nPoints; j++)
        {
            if (x[j] == 13.0)
            {
                meanPtAt13TeV[i] = y[j];
                meanPtAt13TeV_err[i] = gMeanPt[i]->GetErrorY(j);

                cout << "Mean pT of " << particlesLatex[i] << " at 13 TeV is " << meanPtAt13TeV[i] << " +- " << meanPtAt13TeV_err[i] << endl;
                break;
            }
        }
        if (i < 5)
        {
            gMeanPtvsMassMesons->SetPoint(i, particleMass[i], meanPtAt13TeV[i]);
            gMeanPtvsMassMesons->SetPointError(i, 0, meanPtAt13TeV_err[i]);
        }
        else
        {
            gMeanPtvsMassBaryons->SetPoint(i - 5, particleMass[i], meanPtAt13TeV[i]);
            gMeanPtvsMassBaryons->SetPointError(i - 5, 0, meanPtAt13TeV_err[i]);
        }
    }
    SetGrapherrorStyle(gMeanPtvsMassMesons);
    gMeanPtvsMassMesons->SetMarkerStyle(22);
    gMeanPtvsMassMesons->SetMarkerColor(kMagenta);
    gMeanPtvsMassMesons->SetLineColor(kMagenta);
    gMeanPtvsMassMesons->SetMarkerSize(1.5);
    SetGrapherrorStyle(gMeanPtvsMassBaryons);
    gMeanPtvsMassBaryons->SetMarkerStyle(34);
    gMeanPtvsMassBaryons->SetMarkerColor(kRed);
    gMeanPtvsMassBaryons->SetLineColor(kRed);
    gMeanPtvsMassBaryons->SetMarkerSize(1.5);

    TCanvas *cMeanPt = new TCanvas("cMeanPt", "Mean pT vs mass", 720, 720);
    SetCanvasStyle(cMeanPt, 0.14, 0.03, 0.05, 0.14);
    gMeanPtvsMassMesons->GetXaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    gMeanPtvsMassMesons->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/#it{c})");
    gMeanPtvsMassMesons->GetYaxis()->SetTitleOffset(1.3);
    gMeanPtvsMassMesons->SetMinimum(0.0);
    gMeanPtvsMassMesons->SetMaximum(3.15);
    gMeanPtvsMassMesons->GetXaxis()->SetLimits(0, 1.8);
    gMeanPtvsMassMesons->Draw("AP");
    gMeanPtvsMassBaryons->Draw("P SAME");

    TF1 *pol1_meson = new TF1("pol1_meson", "pol1", 0.1, 1.8);
    pol1_meson->SetLineColor(kYellow + 2);
    pol1_meson->SetLineStyle(2);
    gMeanPtvsMassMesons->Fit(pol1_meson, "R");

    TF1 *pol1_baryon = new TF1("pol1_baryon", "pol1", 0.1, 1.8);
    pol1_baryon->SetLineColor(kYellow + 2);
    pol1_baryon->SetLineStyle(2);
    gMeanPtvsMassBaryons->Fit(pol1_baryon, "R");

    // // Draw each point with its own marker style and color
    // for (int i = 0; i < totalParticles; i++)
    // {
    //     cMeanPt->cd();
    //     TMarker *marker = new TMarker(particleMass[i], meanPtAt13TeV[i], markers[i]);
    //     marker->SetMarkerColor(colors[i]);
    //     marker->SetMarkerSize(2.0); // increased size
    //     marker->Draw("SAME");

    //     // Draw error bar for this point in the same color
    //     double x = particleMass[i];
    //     double y = meanPtAt13TeV[i];
    //     double errY = meanPtAt13TeV_err[i];
    //     TLine *errBar = new TLine(x, y - errY, x, y + errY);
    //     errBar->SetLineColor(colors[i]);
    //     errBar->SetLineWidth(2);
    //     errBar->Draw("SAME");
    //     // Add horizontal lines at the ends (caps)
    //     double cap = 0.01; // adjust cap width as needed
    //     TLine *capLow = new TLine(x - cap, y - errY, x + cap, y - errY);
    //     TLine *capHigh = new TLine(x - cap, y + errY, x + cap, y + errY);
    //     capLow->SetLineColor(colors[i]);
    //     capHigh->SetLineColor(colors[i]);
    //     capLow->SetLineWidth(2);
    //     capHigh->SetLineWidth(2);
    //     capLow->Draw("SAME");
    //     capHigh->Draw("SAME");
    // }

    // Draw the last marker (f2(1525)) and its error bar
    double f2_mass = 1.5173;
    double f2_meanpt = hout->GetBinContent(5);
    double f2_err = hout->GetBinContent(6);
    int f2_marker = 20;  // choose a unique marker style for f2(1525)
    int f2_color = kGreen + 1; // choose a unique color for f2(1525)
    TMarker *marker_f2 = new TMarker(f2_mass, f2_meanpt, f2_marker);
    marker_f2->SetMarkerColor(f2_color);
    marker_f2->SetMarkerSize(1.5);
    marker_f2->Draw("SAME");
    TLine *errBar_f2 = new TLine(f2_mass, f2_meanpt - f2_err, f2_mass, f2_meanpt + f2_err);
    errBar_f2->SetLineColor(f2_color);
    errBar_f2->SetLineWidth(2);
    errBar_f2->Draw("SAME");
    double cap_f2 = 0.01;
    TLine *capLow_f2 = new TLine(f2_mass - cap_f2, f2_meanpt - f2_err, f2_mass + cap_f2, f2_meanpt - f2_err);
    TLine *capHigh_f2 = new TLine(f2_mass - cap_f2, f2_meanpt + f2_err, f2_mass + cap_f2, f2_meanpt + f2_err);
    capLow_f2->SetLineColor(f2_color);
    capHigh_f2->SetLineColor(f2_color);
    capLow_f2->SetLineWidth(2);
    capHigh_f2->SetLineWidth(2);
    capLow_f2->Draw("SAME");
    capHigh_f2->Draw("SAME");

    // Draw the last marker (f0(1710)) and its error bar
    double f0_mass = 1.710;
    double f0_meanpt = hout2->GetBinContent(5);
    double f0_err = hout2->GetBinContent(6);
    int f0_marker = 21;   // choose a unique marker style for f0(1710)
    int f0_color = kBlue; // choose a unique color for f0(1710)
    TMarker *marker_f0 = new TMarker(f0_mass, f0_meanpt, f0_marker);
    marker_f0->SetMarkerColor(f0_color);
    marker_f0->SetMarkerSize(1.5);
    marker_f0->Draw("SAME");
    TLine *errBar_f0 = new TLine(f0_mass, f0_meanpt - f0_err, f0_mass, f0_meanpt + f0_err);
    errBar_f0->SetLineColor(f0_color);
    errBar_f0->SetLineWidth(2);
    errBar_f0->Draw("SAME");
    double cap_f0 = 0.01;
    TLine *capLow_f0 = new TLine(f0_mass - cap_f0, f0_meanpt - f0_err, f0_mass + cap_f0, f0_meanpt - f0_err);
    TLine *capHigh_f0 = new TLine(f0_mass - cap_f0, f0_meanpt + f0_err, f0_mass + cap_f0, f0_meanpt + f0_err);
    capLow_f0->SetLineColor(f0_color);
    capHigh_f0->SetLineColor(f0_color);
    capLow_f0->SetLineWidth(2);
    capHigh_f0->SetLineWidth(2);
    capLow_f0->Draw("SAME");
    capHigh_f0->Draw("SAME");

    // Add particle names below each point
    TLatex latex;
    latex.SetTextAlign(22);
    latex.SetTextSize(0.035);
    for (int i = 0; i < totalParticles; i++)
    {
        double x = particleMass[i];
        double y = meanPtAt13TeV[i];
        // latex.SetTextColor(colors[i]);
        if (i == 2)
            latex.DrawLatex(x + 0.08, y + 0.18, particlesLatex[i].c_str());
        else if (i <= 4)
            latex.DrawLatex(x, y + 0.18, particlesLatex[i].c_str());
        else
            latex.DrawLatex(x, y - 0.16, particlesLatex[i].c_str());
    }
    // latex.SetTextColor(kRed); // match f2(1525) marker color
    latex.DrawLatex(1.5173, hout->GetBinContent(5) + 0.22, "f'_{2}(1525)");
    // latex.SetTextColor(kBlue); // match f0(1710) marker color
    latex.DrawLatex(1.710, hout2->GetBinContent(5) + 0.29, "f_{0}(1710)");

    TLegend *legend4 = new TLegend(0.18, 0.65, 0.65, 0.92);
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->SetTextSize(0.035);
    legend4->AddEntry(gMeanPtvsMassMesons, "Mesons (13 TeV)", "p");
    legend4->AddEntry(gMeanPtvsMassBaryons, "Baryons (13 TeV)", "p");
    legend4->AddEntry(marker_f2, "f'_{2}(1525) (13.6 TeV)", "p");
    legend4->AddEntry(marker_f0, "f_{0}(1710) (13.6 TeV)", "p");
    legend4->AddEntry(pol1_meson, "Pol 1", "l");
    legend4->Draw();
    cMeanPt->SaveAs(outputPath + "/MeanPt_vs_Mass.png");

    TGraphErrors *gMeanPtBothTemp = new TGraphErrors();
    gMeanPtBothTemp->SetPoint(0, 1.5173, hout->GetBinContent(5)); // PDG mass of f2(1525) is 1.5173GeV/c2
    gMeanPtBothTemp->SetPointError(0, 0, hout->GetBinContent(6));
    gMeanPtBothTemp->SetPoint(1, 1.710, hout2->GetBinContent(5)); // PDG mass of f0(1710) is 1.710GeV/c2
    gMeanPtBothTemp->SetPointError(1, 0, hout2->GetBinContent(6));
    fOutput->cd();
    gMeanPtBothTemp->Write("gMeanPt_f0f2");

    cout << "dN/dy of f0(1710) = " << hout2->GetBinContent(1) << " +- " << hout2->GetBinContent(2) << endl;
    cout << "dN/dy of f2(1525) = " << hout->GetBinContent(1) << " +- " << hout->GetBinContent(2) << endl;
    cout << "Ratio dN/dy f0(1710)/f2(1525) = " << hout2->GetBinContent(1) / hout->GetBinContent(1) << " +- " << TMath::Sqrt(TMath::Power(hout2->GetBinContent(2) / hout->GetBinContent(1), 2) + TMath::Power(hout->GetBinContent(2) * hout2->GetBinContent(1) / (hout->GetBinContent(1) * hout->GetBinContent(1)), 2)) << endl;
    cout << "Ratio from thermal model " << 0.201577 / 0.781416 << endl;
}