#include <iostream>
#include "../src/style.h"
using namespace std;

void read_yield_cosTheta()
{
    gStyle->SetOptStat(0);
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/3pTcut/MIX/";
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/nopTcut/modified_boltzmann/";

    float cosThetaBins[] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    // float cosThetaBins[] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
    int nBins = sizeof(cosThetaBins) / sizeof(cosThetaBins[0]) - 1;

    // efficiency file
    TFile *feff = new TFile("/home/sawan/check_k892/mc/LHC24l1/463655.root", "read");
    if (feff->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    string histpath = "higher-mass-resonances/hMChists";
    string histpathf2 = "higher-mass-resonances_f21525/hMChists";

    THnSparseD *GenpTf0 = (THnSparseD *)feff->Get(Form("%s/Genf17102", histpath.c_str())); // axis: multiplicity, pt, helicity angle
    THnSparseD *GenpTf2 = (THnSparseD *)feff->Get(Form("%s/Genf17102", histpathf2.c_str()));
    THnSparseD *recpt1f0 = (THnSparseD *)feff->Get(Form("%s/Recf1710_pt2", histpath.c_str())); // axis: multiplicity, pt, mass, helicity angle
    THnSparseD *recpt1f2 = (THnSparseD *)feff->Get(Form("%s/Recf1710_pt2", histpathf2.c_str()));

    if (GenpTf0 == nullptr || recpt1f0 == nullptr || GenpTf2 == nullptr || recpt1f2 == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    TH1D *hgencosThetaStarf0 = GenpTf0->Projection(2);  // project on helicity angle axis
    TH1D *hreccosThetaStarf0 = recpt1f0->Projection(3); // project on helicity angle axis
    TH1D *heffCosThetaf0 = new TH1D("heffCosThetaf0", "Efficiency CosThetaStar", nBins, cosThetaBins);
    TH1D *hgencosThetaStarf2 = GenpTf2->Projection(2);  // project on helicity angle axis
    TH1D *hreccosThetaStarf2 = recpt1f2->Projection(3); // project on helicity angle axis
    TH1D *heffCosThetaf2 = new TH1D("heffCosThetaf2", "Efficiency CosThetaStar for f21525", nBins, cosThetaBins);

    for (int i = 0; i < nBins; i++)
    {
        // get bin content accroding to cosTheta bins and error according to bayesian method
        int lowcosbin = hgencosThetaStarf0->GetXaxis()->FindBin(cosThetaBins[i] + 0.01);
        int highcosbin = hgencosThetaStarf0->GetXaxis()->FindBin(cosThetaBins[i + 1] - 0.01);
        double genYieldCosf0 = hgencosThetaStarf0->Integral(lowcosbin, highcosbin);
        double recYieldCosf0 = hreccosThetaStarf0->Integral(lowcosbin, highcosbin);
        double recYieldErrorCos = TMath::Sqrt(((recYieldCosf0 + 1) / (genYieldCosf0 + 2)) * ((recYieldCosf0 + 2) / (genYieldCosf0 + 3) - (recYieldCosf0 + 1) / (genYieldCosf0 + 2)));
        if (genYieldCosf0 > 0)
        {
            heffCosThetaf0->SetBinContent(i + 1, recYieldCosf0 / genYieldCosf0);
            heffCosThetaf0->SetBinError(i + 1, recYieldErrorCos);
        }

        double genYieldCosf2 = hgencosThetaStarf2->Integral(lowcosbin, highcosbin);
        double recYieldCosf2 = hreccosThetaStarf2->Integral(lowcosbin, highcosbin);
        double recYieldErrorCosf2 = TMath::Sqrt(((recYieldCosf2 + 1) / (genYieldCosf2 + 2)) * ((recYieldCosf2 + 2) / (genYieldCosf2 + 3) - (recYieldCosf2 + 1) / (genYieldCosf2 + 2)));
        if (genYieldCosf2 > 0)
        {
            heffCosThetaf2->SetBinContent(i + 1, recYieldCosf2 / genYieldCosf2);
            heffCosThetaf2->SetBinError(i + 1, recYieldErrorCosf2);
        }
    }

    TH1D *hYield1710Raw = new TH1D("hYield1710Raw", "Yield vs cosTheta", nBins, cosThetaBins);
    TH1D *hYield1525Raw = new TH1D("hYield1525Raw", "Yield vs cosTheta", nBins, cosThetaBins);
    TH1D *hYield1710Corrected = new TH1D("hYield1710Corrected", "Yield vs cosTheta", nBins, cosThetaBins);
    TH1D *hYield1525Corrected = new TH1D("hYield1525Corrected", "Yield vs cosTheta", nBins, cosThetaBins);
    TH1D *hMass1710 = new TH1D("hMass1710", "Mass vs cosTheta", nBins, cosThetaBins);
    TH1D *hWidth1710 = new TH1D("hWidth1710", "Width vs cosTheta", nBins, cosThetaBins);
    TH1D *hMass1525 = new TH1D("hMass1525", "Mass vs cosTheta", nBins, cosThetaBins);
    TH1D *hWidth1525 = new TH1D("hWidth1525", "Width vs cosTheta", nBins, cosThetaBins);
    TH1D *hYieldBkg1 = new TH1D("hYieldBkg1", "Background Yield 1 vs cosTheta", nBins, cosThetaBins);
    TH1D *hYieldBkg2 = new TH1D("hYieldBkg2", "Background Yield 2 vs cosTheta", nBins, cosThetaBins);
    TH1D *hYieldBkg3 = new TH1D("hYieldBkg3", "Background Yield 3 vs cosTheta", nBins, cosThetaBins);
    double totalYield1710, totalYield1525;

    for (int i = 0; i < nBins; i++)
    {
        if (nBins != heffCosThetaf0->GetNbinsX())
        {
            cout << "number of costheta bins " << nBins << ", number of efficiency bins " << heffCosThetaf0->GetNbinsX() << endl;
            cout << "Bins mismatch " << endl;
            return;
        }

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
        double eff1710 = heffCosThetaf0->GetBinContent(i + 1);
        double eff1525 = heffCosThetaf2->GetBinContent(i + 1);

        hYield1710Raw->SetBinContent(i + 1, yield1);
        hYield1710Raw->SetBinError(i + 1, yield1_err);
        hYield1710Corrected->SetBinContent(i + 1, yield1/eff1710);
        hYield1710Corrected->SetBinError(i + 1, 0);
        hMass1710->SetBinContent(i + 1, mass1);
        hMass1710->SetBinError(i + 1, mass1_err);
        hWidth1710->SetBinContent(i + 1, width1);
        hWidth1710->SetBinError(i + 1, width1_err);

        hYield1525Raw->SetBinContent(i + 1, yield2);
        hYield1525Raw->SetBinError(i + 1, yield2_err);
        hYield1525Corrected->SetBinContent(i + 1, yield2/eff1525);
        hYield1525Corrected->SetBinError(i + 1, 0);
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

        cout << "f1525 yield is " << yield2 << " +- " << yield2_err << endl;
        cout << "f1710 yield is " << yield1 << " +- " << yield1_err << endl;
        // cout<<"f1525 mass is "<<mass2<<" +- "<<mass2_err<<endl;
        // cout<<"f1525 width is "<<width2<<" +- "<<width2_err<<endl;
    }

    cout << "Integrated yield of 1710 is " << totalYield1710 << endl;
    cout << "Integrated yield of 1525 is " << totalYield1525 << endl;

    TCanvas *cYield = new TCanvas("cYield", "Yield vs cosTheta", 720, 720);
    SetCanvasStyle(cYield, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hYield1710Raw);
    SetHistoQA(hYield1525Raw);
    hYield1710Raw->GetXaxis()->SetTitle("cos(#theta)");
    hYield1710Raw->GetYaxis()->SetTitle("1/N_{ev} * d^{2}N/(dy dCos#theta)");
    hYield1710Raw->GetYaxis()->SetTitleOffset(1.6);
    // hYield1710Raw->GetYaxis()->SetRangeUser(0.0, 460);
    hYield1710Raw->SetMaximum(hYield1710Raw->GetMaximum() * 2.0);
    hYield1710Raw->SetMarkerStyle(21);
    hYield1710Raw->SetMinimum(-1e-7);
    hYield1710Raw->Draw("pe");
    hYield1525Raw->SetMarkerColor(kRed);
    hYield1525Raw->SetLineColor(kRed);
    hYield1525Raw->Draw("pe same");
    SetHistoQA(hYieldBkg1);
    SetHistoQA(hYieldBkg2);
    SetHistoQA(hYieldBkg3);
    hYieldBkg1->SetMarkerColor(kBlack);
    hYieldBkg1->SetLineColor(kBlack);
    hYieldBkg1->GetXaxis()->SetTitle("cos(#theta)");
    hYieldBkg1->GetYaxis()->SetTitle("Background Yield (1/N_{ev} * dN/(dCos#theta))");
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
    legend->AddEntry(hYield1710Raw, "f_{0}(1710)", "lpe");
    legend->AddEntry(hYield1525Raw, "f'_{2}(1525)", "lpe");
    legend->Draw();
    // cYield->SaveAs((path + "yield_vs_cosTheta.png").c_str());

    TCanvas *cYieldCorrected = new TCanvas("cYieldCorrected", "Yield vs cosTheta (corrected)", 720, 720);
    SetCanvasStyle(cYieldCorrected, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hYield1710Corrected);
    SetHistoQA(hYield1525Corrected);
    hYield1710Corrected->GetXaxis()->SetTitle("cos(#theta)");
    hYield1710Corrected->GetYaxis()->SetTitle("1/N_{ev} * d^{2}N/(dy dCos#theta)");
    hYield1710Corrected->GetYaxis()->SetTitleOffset(1.6);
    // hYield1710Corrected->GetYaxis()->SetRangeUser(0.0, 460);
    hYield1710Corrected->SetMaximum(hYield1710Corrected->GetMaximum() * 2.0);
    hYield1710Corrected->SetMarkerStyle(21);
    hYield1710Corrected->SetMinimum(-1e-7);
    hYield1710Corrected->Draw("pe");
    hYield1525Corrected->SetMarkerColor(kRed);
    hYield1525Corrected->SetLineColor(kRed);
    hYield1525Corrected->Draw("pe same");
    legend->Draw();

    // TCanvas *cWidth = new TCanvas("cWidth", "Width vs cosTheta", 720, 720);
    // SetCanvasStyle(cWidth, 0.16, 0.03, 0.05, 0.14);
    // SetHistoQA(hWidth1710);
    // SetHistoQA(hWidth1525);
    // hWidth1710->SetMarkerStyle(21);
    // hWidth1710->GetXaxis()->SetTitle("cos(#theta)");
    // hWidth1710->GetYaxis()->SetTitle("Width (GeV)");
    // hWidth1710->GetYaxis()->SetTitleOffset(1.5);
    // hWidth1710->GetYaxis()->SetRangeUser(0.0, 0.3);
    // hWidth1710->Draw("pe");
    // hWidth1525->SetMarkerColor(kRed);
    // hWidth1525->SetLineColor(kRed);
    // // hWidth1525->Draw("pe same");
    // // legend->Draw();
    // // cWidth->SaveAs((path + "width_vs_cosTheta.png").c_str());

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
