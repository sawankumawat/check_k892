#include <iostream>
#include "../src/style.h"
using namespace std;

void read_yield_cosTheta()
{
    gStyle->SetOptStat(0);
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937/fits/4rBw_fits/";
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/nopTcut/";

    float cosThetaBins[] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    // float cosThetaBins[] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
    int nBins = sizeof(cosThetaBins) / sizeof(cosThetaBins[0]) - 1;

    TH1D *hYield1710 = new TH1D("hYield1710", "Yield vs cosTheta", nBins, cosThetaBins);
    TH1D *hMass1710 = new TH1D("hMass1710", "Mass vs cosTheta", nBins, cosThetaBins);
    TH1D *hWidth1710 = new TH1D("hWidth1710", "Width vs cosTheta", nBins, cosThetaBins);
    TH1D *hYield1525 = new TH1D("hYield1525", "Yield vs cosTheta", nBins, cosThetaBins);
    TH1D *hMass1525 = new TH1D("hMass1525", "Mass vs cosTheta", nBins, cosThetaBins);
    TH1D *hWidth1525 = new TH1D("hWidth1525", "Width vs cosTheta", nBins, cosThetaBins);
    TH1D *hYieldBkg1 = new TH1D("hYieldBkg1", "Background Yield 1 vs cosTheta", nBins, cosThetaBins);
    TH1D *hYieldBkg2 = new TH1D("hYieldBkg2", "Background Yield 2 vs cosTheta", nBins, cosThetaBins);
    TH1D *hYieldBkg3 = new TH1D("hYieldBkg3", "Background Yield 3 vs cosTheta", nBins, cosThetaBins);
    double totalYield1710, totalYield1525;

    for (int i = 0; i < nBins; i++)
    {
        ifstream infile;
        infile.open(path + Form("fit_params_cosTheta_%.1f_%.1f.txt", cosThetaBins[i], cosThetaBins[i + 1]));
        // cout << "path is " << path + "fit_params_cosTheta_" + to_string(cosThetaBins[i]) + "_" + to_string(cosThetaBins[i + 1]) + ".txt" << endl;

        if (!infile.is_open())
        {
            cerr << "Error opening file!" << endl;
            return;
        }
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

        totalYield1710 += yield1;
        totalYield1525 += yield2;

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

        cout<<"f1525 yield is "<<yield2<<" +- "<<yield2_err<<endl;
        cout<<"f1710 yield is "<<yield1<<" +- "<<yield1_err<<endl;
        // cout<<"f1525 mass is "<<mass2<<" +- "<<mass2_err<<endl;
        // cout<<"f1525 width is "<<width2<<" +- "<<width2_err<<endl;
    }

    cout << "Integrated yield of 1710 is " << totalYield1710 << endl;
    cout << "Integrated yield of 1525 is " << totalYield1525 << endl;

    TCanvas *cYield = new TCanvas("cYield", "Yield vs cosTheta", 720, 720);
    SetCanvasStyle(cYield, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hYield1710);
    SetHistoQA(hYield1525);
    hYield1710->GetXaxis()->SetTitle("cos(#theta)");
    hYield1710->GetYaxis()->SetTitle("1/N_{ev} * d^{2}N/(dy dCos#theta)");
    hYield1710->GetYaxis()->SetTitleOffset(1.6);
    // hYield1710->GetYaxis()->SetRangeUser(0.0, 460);
    hYield1710->SetMaximum(hYield1710->GetMaximum() * 2.0);
    hYield1710->SetMarkerStyle(21);
    hYield1710->SetMinimum(-1e-7);
    hYield1710->Draw("pe");
    hYield1525->SetMarkerColor(kRed);
    hYield1525->SetLineColor(kRed);
    hYield1525->Draw("pe same");
    SetHistoQA(hYieldBkg1);
    SetHistoQA(hYieldBkg2);
    SetHistoQA(hYieldBkg3);
    hYieldBkg1->SetMarkerColor(kBlack);
    hYieldBkg1->SetLineColor(kBlack);
    hYieldBkg1->GetXaxis()->SetTitle("cos(#theta)");
    hYieldBkg1->GetYaxis()->SetTitle("Background Yield (1/N_{ev} * d^{2}N/(dy dCos#theta))");
    hYieldBkg1->GetYaxis()->SetTitleOffset(1.6);
    hYieldBkg1->SetMaximum(hYieldBkg1->GetMaximum() * 2.0);
    hYieldBkg1->Draw("pe same");
    hYieldBkg2->SetMarkerColor(kMagenta - 2);
    hYieldBkg2->SetLineColor(kMagenta - 2);
    hYieldBkg2->Draw("pe same");
    hYieldBkg3->SetMarkerColor(kGreen + 2);
    hYieldBkg3->SetLineColor(kGreen + 2);
    hYieldBkg3->Draw("pe same");

    TLegend *legend = new TLegend(0.77, 0.83, 0.9, 0.93);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.037);
    legend->SetTextFont(42);
    legend->AddEntry(hYield1710, "f_{0}(1710)", "lpe");
    legend->AddEntry(hYield1525, "f'_{2}(1525)", "lpe");
    legend->Draw();
    // cYield->SaveAs((path + "yield_vs_cosTheta.png").c_str());

    TCanvas *cMass = new TCanvas("cMass", "Mass vs cosTheta", 720, 720);
    SetCanvasStyle(cMass, 0.16, 0.03, 0.05, 0.14);
    SetHistoQA(hMass1710);
    SetHistoQA(hMass1525);
    hMass1710->SetMarkerStyle(21);
    hMass1710->GetXaxis()->SetTitle("cos(#theta)");
    hMass1710->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    hMass1710->GetYaxis()->SetTitleOffset(1.5);
    hMass1710->GetYaxis()->SetRangeUser(1.42, 1.79);
    hMass1710->Draw("pe");
    hMass1525->SetMarkerColor(kRed);
    hMass1525->SetLineColor(kRed);
    hMass1525->Draw("pe same");
    legend->Draw();
    // cMass->SaveAs((path + "mass_vs_cosTheta.png").c_str());

    TCanvas *cWidth = new TCanvas("cWidth", "Width vs cosTheta", 720, 720);
    SetCanvasStyle(cWidth, 0.16, 0.03, 0.05, 0.14);
    SetHistoQA(hWidth1710);
    SetHistoQA(hWidth1525);
    hWidth1710->SetMarkerStyle(21);
    hWidth1710->GetXaxis()->SetTitle("cos(#theta)");
    hWidth1710->GetYaxis()->SetTitle("Width (GeV)");
    hWidth1710->GetYaxis()->SetTitleOffset(1.5);
    hWidth1710->GetYaxis()->SetRangeUser(0.0, 0.3);
    hWidth1710->Draw("pe");
    hWidth1525->SetMarkerColor(kRed);
    hWidth1525->SetLineColor(kRed);
    // hWidth1525->Draw("pe same");
    // legend->Draw();
    // cWidth->SaveAs((path + "width_vs_cosTheta.png").c_str());

    // TCanvas *cYieldBkg = new TCanvas("cYieldBkg", "Background Yield vs cosTheta", 720, 720);
    // SetCanvasStyle(cYieldBkg, 0.16, 0.03, 0.05, 0.14);
    // SetHistoQA(hYieldBkg1);
    // SetHistoQA(hYieldBkg2);
    // SetHistoQA(hYieldBkg3);
    // hYieldBkg1->SetMarkerStyle(21);
    // hYieldBkg1->GetXaxis()->SetTitle("cos(#theta)");
    // hYieldBkg1->GetYaxis()->SetTitle("Background Yield (1/N_{ev} * d^{2}N/(dy dCos#theta))");
    // hYieldBkg1->GetYaxis()->SetTitleOffset(1.6);
    // hYieldBkg1->SetMaximum(hYieldBkg1->GetMaximum() * 2.0);
    // hYieldBkg1->Draw("pe");
    // hYieldBkg2->SetMarkerColor(kRed);
    // hYieldBkg2->SetLineColor(kRed);
    // hYieldBkg2->Draw("pe same");
    // hYieldBkg3->SetMarkerColor(kGreen + 2);
    // hYieldBkg3->SetLineColor(kGreen + 2);
    // hYieldBkg3->Draw("pe same");
}
