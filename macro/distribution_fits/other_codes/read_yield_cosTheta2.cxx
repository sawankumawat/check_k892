#include <iostream>
#include "../src/style.h"
using namespace std;

struct FitHistograms
{
    TH1D *hYield1710, *hYield1525, *hMass1710, *hMass1525, *hWidth1710, *hWidth1525, *hYieldBkg1, *hYieldBkg2, *hYieldBkg3;
};

FitHistograms ReadFitData(const string &path, const string &namePrefix)
{
    float cosThetaBins[] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    int nBins = sizeof(cosThetaBins) / sizeof(cosThetaBins[0]) - 1;

    // Declare histogram pointers
    TH1D *hYield1710 = new TH1D((namePrefix + "_yield1710").c_str(), "Yield vs cosTheta", nBins, cosThetaBins);
    TH1D *hYield1525 = new TH1D((namePrefix + "_yield1525").c_str(), "Yield vs cosTheta", nBins, cosThetaBins);
    TH1D *hMass1710 = new TH1D((namePrefix + "_mass1710").c_str(), "Mass vs cosTheta", nBins, cosThetaBins);
    TH1D *hMass1525 = new TH1D((namePrefix + "_mass1525").c_str(), "Mass vs cosTheta", nBins, cosThetaBins);
    TH1D *hWidth1710 = new TH1D((namePrefix + "_width1710").c_str(), "Width vs cosTheta", nBins, cosThetaBins);
    TH1D *hWidth1525 = new TH1D((namePrefix + "_width1525").c_str(), "Width vs cosTheta", nBins, cosThetaBins);
    TH1D *hYieldBkg1 = new TH1D((namePrefix + "_yieldBkg1").c_str(), "Background Yield 1 vs cosTheta", nBins, cosThetaBins);
    TH1D *hYieldBkg2 = new TH1D((namePrefix + "_yieldBkg2").c_str(), "Background Yield 2 vs cosTheta", nBins, cosThetaBins);
    TH1D *hYieldBkg3 = new TH1D((namePrefix + "_yieldBkg3").c_str(), "Background Yield 3 vs cosTheta", nBins, cosThetaBins);

    for (int i = 0; i < nBins; i++)
    {
        ifstream infile;
        infile.open(path + Form("fit_params_cosTheta_%.1f_%.1f.txt", cosThetaBins[i], cosThetaBins[i + 1]));
        // cout << "path is " << path + "fit_params_cosTheta_" + to_string(cosThetaBins[i]) + "_" + to_string(cosThetaBins[i + 1]) + ".txt" << endl;

        // if (!infile.is_open())
        // {
        //     cerr << "Error opening file!" << endl;
        //     return;
        // }

        std::string line;

        // Common fit info
        double significance = 0, statSignificance = 0, chi2ndf = 0;

        // Parameters for f1710
        double yield1 = 0, yield1_err = 0;
        double mass1 = 0, mass1_err = 0;
        double width1 = 0, width1_err = 0;
        double bkgYield1 = 0, bkgYield1_err = 0, bkgYield2 = 0, bkgYield2_err = 0, bkgYield3 = 0, bkgYield3_err = 0;

        // Parameters for f1525
        double yield2 = 0, yield2_err = 0;
        double mass2 = 0, mass2_err = 0;
        double width2 = 0, width2_err = 0;

        int paramIndex = 0;
        bool readingF1710 = false;
        bool readingF1525 = false;
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
                    readingBkg = false;
                    paramIndex = 0;
                }

                else if (line.find("f1525") != std::string::npos)
                {
                    readingF1525 = true;
                    readingF1710 = false;
                    readingBkg = false;
                    paramIndex = 0;
                }
                else if (line.find("Background") != std::string::npos)
                {
                    readingBkg = true;
                    readingF1710 = false;
                    readingF1525 = false;
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
                            mass1 = val;
                            mass1_err = err;
                            break;
                        case 2:
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
                            mass2 = val;
                            mass2_err = err;
                            break;
                        case 2:
                            width2 = val;
                            width2_err = err;
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

        hYield1710->SetBinContent(i + 1, yield1);
        hYield1710->SetBinError(i + 1, yield1_err);
        hMass1710->SetBinContent(i + 1, mass1);
        hMass1710->SetBinError(i + 1, mass1_err);
        hWidth1710->SetBinContent(i + 1, width1);
        hWidth1710->SetBinError(i + 1, width1_err);

        hYield1525->SetBinContent(i + 1, yield2);
        hYield1525->SetBinError(i + 1, yield2_err);
        hMass1525->SetBinContent(i + 1, mass2);
        hMass1525->SetBinError(i + 1, mass2_err);
        hWidth1525->SetBinContent(i + 1, width2);
        hWidth1525->SetBinError(i + 1, width2_err);

        hYieldBkg1->SetBinContent(i + 1, bkgYield1);
        hYieldBkg1->SetBinError(i + 1, bkgYield1 / 20);
        hYieldBkg2->SetBinContent(i + 1, bkgYield2);
        hYieldBkg2->SetBinError(i + 1, bkgYield2 / 20);
        hYieldBkg3->SetBinContent(i + 1, bkgYield3);
        hYieldBkg3->SetBinError(i + 1, bkgYield3 / 20);

        // cout << "f1525 yield is " << yield2 << " +- " << yield2_err << endl;
        // cout << "f1710 yield is " << yield1 << " +- " << yield1_err << endl;
        // cout<<"f1525 mass is "<<mass2<<" +- "<<mass2_err<<endl;
        // cout<<"f1525 width is "<<width2<<" +- "<<width2_err<<endl;
    }

    return {hYield1710, hYield1525, hMass1710, hMass1525, hWidth1710, hWidth1525, hYieldBkg1, hYieldBkg2, hYieldBkg3};
}

void read_yield_cosTheta2()
{
    gStyle->SetOptStat(0);

    string path1 = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937/fits/4rBw_fits/";
    // string path2 = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/3pTcut/MIX/";
    string path2 = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/";

    FitHistograms h1 = ReadFitData(path1, "data1");
    FitHistograms h2 = ReadFitData(path2, "data2");

    // TCanvas *c = new TCanvas("c", "Yield Comparison", 720, 720);
    // SetCanvasStyle(c, 0.18, 0.03, 0.05, 0.14);

    // SetHistoQA(h1.yield1710);
    // SetHistoQA(h2.yield1710);

    // h1.yield1710->SetMarkerStyle(20);
    // h1.yield1710->SetMarkerColor(kBlue + 1);
    // h1.yield1710->SetLineColor(kBlue + 1);
    // h1.yield1710->GetXaxis()->SetTitle("cos(#theta)");
    // h1.yield1710->GetYaxis()->SetTitle("Yield");
    // h1.yield1710->Draw("pe");

    // h2.yield1710->SetMarkerStyle(21);
    // h2.yield1710->SetMarkerColor(kRed);
    // h2.yield1710->SetLineColor(kRed);
    // h2.yield1710->Draw("pe same");

    // auto leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    // leg->AddEntry(h1.yield1710, "Dataset 1", "p");
    // leg->AddEntry(h2.yield1710, "Dataset 2", "p");
    // leg->Draw();

    TCanvas *cYield = new TCanvas("cYield", "Yield vs cosTheta", 720, 720);
    SetCanvasStyle(cYield, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(h1.hYield1710);
    SetHistoQA(h1.hYield1525);
    SetHistoQA(h2.hYield1710);
    SetHistoQA(h2.hYield1525);
    SetHistoQA(h1.hYieldBkg1);
    SetHistoQA(h1.hYieldBkg2);
    SetHistoQA(h1.hYieldBkg3);
    h1.hYield1710->GetXaxis()->SetTitle("cos(#theta)");
    h1.hYield1710->GetYaxis()->SetTitle("1/N_{ev} * d^{2}N/(dy dCos#theta)");
    h1.hYield1710->GetYaxis()->SetTitleOffset(1.6);
    // h1.hYield1710->GetYaxis()->SetRangeUser(0.0, 460);
    h1.hYield1710->SetMaximum(h1.hYield1710->GetMaximum() * 2.0);
    h1.hYield1710->SetMarkerStyle(21);
    h1.hYield1710->SetMinimum(-1e-7);
    h1.hYield1710->Draw("pe");
    h1.hYield1525->SetMarkerColor(kRed);
    h1.hYield1525->SetLineColor(kRed);
    h1.hYield1525->Draw("pe same");
    h2.hYield1710->SetMarkerStyle(22);
    h2.hYield1710->SetMarkerColor(kBlue + 1);
    h2.hYield1710->SetLineColor(kBlue + 1);
    h2.hYield1710->Draw("pe same");
    h2.hYield1525->SetMarkerStyle(23);
    h2.hYield1525->SetMarkerColor(kGreen + 2);
    h2.hYield1525->SetLineColor(kGreen + 2);
    h2.hYield1525->Draw("pe same");

    // h1.hYieldBkg1->SetMarkerColor(kBlack);
    // h1.hYieldBkg1->SetLineColor(kBlack);
    // h1.hYieldBkg1->GetXaxis()->SetTitle("cos(#theta)");
    // h1.hYieldBkg1->GetYaxis()->SetTitle("Background Yield (1/N_{ev} * d^{2}N/(dy dCos#theta))");
    // h1.hYieldBkg1->GetYaxis()->SetTitleOffset(1.6);
    // h1.hYieldBkg1->SetMaximum(h1.hYieldBkg1->GetMaximum() * 2.0);
    // h1.hYieldBkg1->Draw("pe same");
    // h1.hYieldBkg2->SetMarkerColor(kMagenta - 2);
    // h1.hYieldBkg2->SetLineColor(kMagenta - 2);
    // h1.hYieldBkg2->Draw("pe same");
    // h1.hYieldBkg3->SetMarkerColor(kGreen + 2);
    // h1.hYieldBkg3->SetLineColor(kGreen + 2);
    // h1.hYieldBkg3->Draw("pe same");

    TLegend *legend = new TLegend(0.77, 0.83, 0.9, 0.93);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.037);
    legend->SetTextFont(42);
    legend->AddEntry(h1.hYield1710, "f_{0}(1710)", "lpe");
    legend->AddEntry(h1.hYield1525, "f'_{2}(1525)", "lpe");
    legend->Draw();
    // cYield->SaveAs((path + "yield_vs_cosTheta.png").c_str());

    TCanvas *cMass = new TCanvas("cMass", "Mass vs cosTheta", 720, 720);
    SetCanvasStyle(cMass, 0.16, 0.03, 0.05, 0.14);
    SetHistoQA(h1.hMass1710);
    SetHistoQA(h1.hMass1525);
    SetHistoQA(h2.hMass1710);
    SetHistoQA(h2.hMass1525);
    h1.hMass1710->SetMarkerStyle(21);
    h1.hMass1710->GetXaxis()->SetTitle("cos(#theta)");
    h1.hMass1710->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    h1.hMass1710->GetYaxis()->SetTitleOffset(1.5);
    h1.hMass1710->GetYaxis()->SetRangeUser(1.42, 1.79);
    h1.hMass1710->Draw("pe");
    h1.hMass1525->SetMarkerColor(kRed);
    h1.hMass1525->SetLineColor(kRed);
    h1.hMass1525->Draw("pe same");
    h2.hMass1710->SetMarkerStyle(22);
    h2.hMass1710->SetMarkerColor(kBlue + 1);
    h2.hMass1710->SetLineColor(kBlue + 1);
    h2.hMass1710->Draw("pe same");
    h2.hMass1525->SetMarkerStyle(23);
    h2.hMass1525->SetMarkerColor(kGreen + 2);
    h2.hMass1525->SetLineColor(kGreen + 2);
    h2.hMass1525->Draw("pe same");
    legend->Draw();
}
