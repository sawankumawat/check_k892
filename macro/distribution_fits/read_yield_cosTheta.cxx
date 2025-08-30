#include <iostream>
#include "../src/style.h"
#include "TF1.h"
#include "TMath.h"
using namespace std;

// Define Legendre polynomial function manually for lower orders
Double_t LegendrePolynomial(int n, Double_t x) {
    switch(n) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5 * (3*x*x - 1);
        case 3: return 0.5 * (5*x*x*x - 3*x);
        case 4: return 0.125 * (35*x*x*x*x - 30*x*x + 3);
        case 5: return 0.125 * (63*x*x*x*x*x - 70*x*x*x + 15*x);
        case 6: return 0.0625 * (231*x*x*x*x*x*x - 315*x*x*x*x + 105*x*x - 5);
        default: return 0.0;
    }
}

// Define Legendre series function
Double_t LegendreSeries(Double_t *x, Double_t *par) {
    double cosTheta = x[0];
    int Jmax = (int)par[0];  // first parameter is Jmax (integer)
    double sum = 0.0;
    for (int J = 0; J <= Jmax; J++) {
        sum += par[J+1] * LegendrePolynomial(J, cosTheta);  // par[J+1] = a_J
    }
    return sum;
}

void read_yield_cosTheta()
{
    gStyle->SetOptStat(0);
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/3pTcut/MIX/";
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/nopTcut/modified_boltzmann/";
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/AngularDistributions/without_ptCut/";
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/AngularDistributions/pTCut1/";

    // float cosThetaBins[] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    // float cosThetaBins[] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
    float cosThetaBins[] = {-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6};
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
        double eff1710_err = heffCosThetaf0->GetBinError(i + 1);
        double eff1525 = heffCosThetaf2->GetBinContent(i + 1);
        double eff1525_err = heffCosThetaf2->GetBinError(i + 1);
        double corrected_yield1710 = yield1 / eff1710;
        double corrected_yield1710_err = sqrt(pow(yield1_err / eff1710, 2) + pow(yield1 * eff1710_err / pow(eff1710, 2), 2));
        double corrected_yield1525 = yield2 / eff1525;
        double corrected_yield1525_err = sqrt(pow(yield2_err / eff1525, 2) + pow(yield2 * eff1525_err / pow(eff1525, 2), 2));

        hYield1710Raw->SetBinContent(i + 1, yield1);
        hYield1710Raw->SetBinError(i + 1, yield1_err);
        hYield1710Corrected->SetBinContent(i + 1, corrected_yield1710);
        hYield1710Corrected->SetBinError(i + 1, corrected_yield1710_err);
        hMass1710->SetBinContent(i + 1, mass1);
        hMass1710->SetBinError(i + 1, mass1_err);
        hWidth1710->SetBinContent(i + 1, width1);
        hWidth1710->SetBinError(i + 1, width1_err);

        hYield1525Raw->SetBinContent(i + 1, yield2);
        hYield1525Raw->SetBinError(i + 1, yield2_err);
        hYield1525Corrected->SetBinContent(i + 1, corrected_yield1525);
        hYield1525Corrected->SetBinError(i + 1, corrected_yield1525_err);
        hMass1525->SetBinContent(i + 1, mass2);
        hMass1525->SetBinError(i + 1, mass2_err);
        hWidth1525->SetBinContent(i + 1, width2);
        hWidth1525->SetBinError(i + 1, width2_err);

        hYieldBkg1->SetBinContent(i + 1, bkgYield1);
        hYieldBkg1->SetBinError(i + 1, bkgYield1);
        hYieldBkg2->SetBinContent(i + 1, bkgYield2);
        hYieldBkg2->SetBinError(i + 1, bkgYield2);
        hYieldBkg3->SetBinContent(i + 1, bkgYield3);
        hYieldBkg3->SetBinError(i + 1, bkgYield3);

        cout << "f1525 yield is " << yield2 << " +- " << yield2_err << endl;
        cout << "f1710 yield is " << yield1 << " +- " << yield1_err << endl;
        // cout<<"f1525 mass is "<<mass2<<" +- "<<mass2_err<<endl;
        // cout<<"f1525 width is "<<width2<<" +- "<<width2_err<<endl;
    }

    cout << "Integrated yield of 1710 is " << totalYield1710 << endl;
    cout << "Integrated yield of 1525 is " << totalYield1525 << endl;

    TCanvas *cYield = new TCanvas("cYield", "Yield vs cosTheta", 720, 720);
    SetCanvasStyle(cYield, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hYield1525Raw);
    SetHistoQA(hYield1710Raw);
    hYield1525Raw->GetXaxis()->SetTitle("cos(#theta)");
    hYield1525Raw->GetYaxis()->SetTitle("1/N_{ev} * dN/(dCos#theta)");
    hYield1525Raw->GetYaxis()->SetTitleOffset(1.6);
    // hYield1525Raw->GetYaxis()->SetRangeUser(0.0, 460);
    hYield1525Raw->SetMaximum(hYield1525Raw->GetMaximum() * 1.5);
    hYield1525Raw->SetMarkerStyle(21);
    hYield1525Raw->SetMinimum(-2e-7);
    hYield1525Raw->SetMarkerColor(kBlue);
    hYield1525Raw->SetLineColor(kBlue);
    hYield1525Raw->Draw("pe");
    hYield1710Raw->SetMarkerStyle(21);
    hYield1710Raw->SetMarkerColor(kRed);
    hYield1710Raw->SetLineColor(kRed);
    hYield1710Raw->Draw("pe same");
    SetHistoQA(hYieldBkg1);
    SetHistoQA(hYieldBkg2);
    SetHistoQA(hYieldBkg3);
    hYieldBkg1->SetMarkerStyle(22);
    hYieldBkg1->SetMarkerColor(kBlack);
    hYieldBkg1->SetLineColor(kBlack);
    hYieldBkg1->Draw("pe same");
    hYieldBkg2->SetMarkerStyle(23);
    hYieldBkg2->SetMarkerColor(kMagenta - 2);
    hYieldBkg2->SetLineColor(kMagenta - 2);
    hYieldBkg2->Draw("pe same");
    hYieldBkg3->SetMarkerStyle(24);
    hYieldBkg3->SetMarkerColor(kGreen + 2);
    hYieldBkg3->SetLineColor(kGreen + 2);
    hYieldBkg3->Draw("pe same");

    TLegend *legend = new TLegend(0.25, 0.75, 0.9, 0.93);
    legend->SetNColumns(2);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);
    legend->AddEntry(hYield1525Raw, "f'_{2}(1525)", "p");
    legend->AddEntry(hYield1710Raw, "f_{0}(1710)", "p");
    legend->AddEntry(hYieldBkg1, "Bkg (2.2-2.3 GeV)", "p");
    legend->AddEntry(hYieldBkg2, "Bkg (2.3-2.4 GeV)", "p");
    legend->AddEntry(hYieldBkg3, "Bkg (2.4-2.5 GeV)", "p");
    legend->Draw();
    cYield->SaveAs((path + "RawYield_vs_cosTheta.png").c_str());

    TCanvas *cYieldCorrected = new TCanvas("cYieldCorrected", "Yield vs cosTheta (corrected)", 720, 720);
    SetCanvasStyle(cYieldCorrected, 0.18, 0.03, 0.05, 0.14);
    SetHistoQA(hYield1525Corrected);
    SetHistoQA(hYield1710Corrected);
    hYield1525Corrected->GetXaxis()->SetTitle("cos(#theta)");
    hYield1525Corrected->GetYaxis()->SetTitle("1/N_{ev} * dN/(dCos#theta)");
    hYield1525Corrected->GetYaxis()->SetTitleOffset(1.6);
    // hYield1525Corrected->GetYaxis()->SetRangeUser(0.0, 460);
    hYield1525Corrected->SetMaximum(hYield1525Corrected->GetMaximum() * 1.5);
    hYield1525Corrected->SetMarkerStyle(21);
    hYield1525Corrected->SetMinimum(-2e-7);
    hYield1525Corrected->SetLineColor(kBlue);
    hYield1525Corrected->SetMarkerColor(kBlue);
    hYield1525Corrected->Draw("pe");
    hYield1525Corrected->SetMarkerStyle(21);
    hYield1710Corrected->SetMarkerColor(kRed);
    hYield1710Corrected->SetLineColor(kRed);
    hYield1710Corrected->Draw("pe same");
    TLegend *legend2 = new TLegend(0.75, 0.75, 0.95, 0.93);
    legend2->SetFillStyle(0);
    legend2->SetBorderSize(0);
    legend2->SetTextFont(42);
    legend2->SetTextSize(0.034);
    legend2->AddEntry(hYield1525Corrected, "f'_{2}(1525)", "p");
    legend2->AddEntry(hYield1710Corrected, "f_{0}(1710)", "p");
    legend2->Draw();
    cYieldCorrected->SaveAs((path + "CorrectedYield_vs_cosTheta.png").c_str());

    // Legendre polynomial fitting for both histograms
    cout << "\n==== Legendre Polynomial Fitting ====" << endl;
    
    // Define maximum order of Legendre polynomials to fit
    int Jmax = 4; // You can adjust this value
    int nParams = Jmax + 2; // one for Jmax + coefficients a_J
    
    // Fit hYield1710Corrected
    cout << "\nFitting f0(1710) corrected yield with Legendre polynomials (Jmax = " << Jmax << "):" << endl;
    TF1 *fLeg1710 = new TF1("fLeg1710", LegendreSeries, -0.6, 0.6, nParams);
    fLeg1710->SetParName(0,"Jmax");
    fLeg1710->FixParameter(0, Jmax); // keep Jmax fixed
    
    // Initial guesses for coefficients
    for (int J=0; J<=Jmax; J++) {
        fLeg1710->SetParName(J+1, Form("a_%d", J));
        fLeg1710->SetParameter(J+1, 1.0);
    }
    
    // Fit the histogram
    hYield1710Corrected->Fit(fLeg1710,"R");
    
    // Print coefficients for f1710
    cout << "Legendre coefficients for f0(1710):" << endl;
    for (int J=0; J<=Jmax; J++) {
        double coeff = fLeg1710->GetParameter(J+1);
        double err = fLeg1710->GetParError(J+1);
        printf("a_%d = %f ± %f\n", J, coeff, err);
    }
    
    // Fit hYield1525Corrected
    cout << "\nFitting f2(1525) corrected yield with Legendre polynomials (Jmax = " << Jmax << "):" << endl;
    TF1 *fLeg1525 = new TF1("fLeg1525", LegendreSeries, -1, 1, nParams);
    fLeg1525->SetParName(0,"Jmax");
    fLeg1525->FixParameter(0, Jmax); // keep Jmax fixed
    
    // Initial guesses for coefficients
    for (int J=0; J<=Jmax; J++) {
        fLeg1525->SetParName(J+1, Form("a_%d", J));
        fLeg1525->SetParameter(J+1, 1.0);
    }
    
    // Fit the histogram
    hYield1525Corrected->Fit(fLeg1525,"R");
    
    // Print coefficients for f1525
    cout << "Legendre coefficients for f2(1525):" << endl;
    for (int J=0; J<=Jmax; J++) {
        double coeff = fLeg1525->GetParameter(J+1);
        double err = fLeg1525->GetParError(J+1);
        printf("a_%d = %f ± %f\n", J, coeff, err);
    }
    
    // Create a canvas to show the fits
    TCanvas *cLegendre = new TCanvas("cLegendre", "Legendre Polynomial Fits", 1200, 600);
    cLegendre->Divide(2,1);
    
    // Plot f1710 fit
    cLegendre->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.14);
    hYield1710Corrected->SetTitle("f_{0}(1710) Legendre Fit");
    hYield1710Corrected->GetXaxis()->SetTitle("cos(#theta)");
    hYield1710Corrected->GetYaxis()->SetTitle("Corrected Yield");
    hYield1710Corrected->SetMarkerStyle(21);
    hYield1710Corrected->SetMarkerColor(kRed);
    hYield1710Corrected->SetLineColor(kRed);
    hYield1710Corrected->Draw("PE");
    fLeg1710->SetLineColor(kBlue);
    fLeg1710->SetLineWidth(2);
    fLeg1710->Draw("same");
    
    TLegend *leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->AddEntry(hYield1710Corrected, "Data", "p");
    leg1->AddEntry(fLeg1710, "Legendre Fit", "l");
    leg1->Draw();
    
    // Plot f1525 fit
    cLegendre->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.14);
    hYield1525Corrected->SetTitle("f'_{2}(1525) Legendre Fit");
    hYield1525Corrected->GetXaxis()->SetTitle("cos(#theta)");
    hYield1525Corrected->GetYaxis()->SetTitle("Corrected Yield");
    hYield1525Corrected->SetMarkerStyle(21);
    hYield1525Corrected->SetMarkerColor(kBlue);
    hYield1525Corrected->SetLineColor(kBlue);
    hYield1525Corrected->Draw("PE");
    fLeg1525->SetLineColor(kRed);
    fLeg1525->SetLineWidth(2);
    fLeg1525->Draw("same");
    
    TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->AddEntry(hYield1525Corrected, "Data", "p");
    leg2->AddEntry(fLeg1525, "Legendre Fit", "l");
    leg2->Draw();
    
    cLegendre->SaveAs((path + "Legendre_Fits_cosTheta.png").c_str());
    
    // Print fit quality information
    cout << "\nFit Quality Information:" << endl;
    cout << "f0(1710) - Chi2/NDF: " << fLeg1710->GetChisquare() << "/" << fLeg1710->GetNDF() 
         << " = " << fLeg1710->GetChisquare()/fLeg1710->GetNDF() << endl;
    cout << "f2(1525) - Chi2/NDF: " << fLeg1525->GetChisquare() << "/" << fLeg1525->GetNDF() 
         << " = " << fLeg1525->GetChisquare()/fLeg1525->GetNDF() << endl;

    // Analysis for spin determination
    cout << "\n==== Spin (J) Analysis ====" << endl;
    
    // For f0(1710) - expected to be J=0 (scalar)
    cout << "\nf0(1710) Analysis (Expected J=0):" << endl;
    cout << "Dominant coefficient should be a_0 for J=0 particle" << endl;
    double max_coeff_1710 = 0;
    int dominant_J_1710 = 0;
    for (int J=0; J<=Jmax; J++) {
        double coeff = TMath::Abs(fLeg1710->GetParameter(J+1));
        if (coeff > max_coeff_1710) {
            max_coeff_1710 = coeff;
            dominant_J_1710 = J;
        }
    }
    cout << "Dominant term: a_" << dominant_J_1710 << " = " << fLeg1710->GetParameter(dominant_J_1710+1) << endl;
    
    // For f2(1525) - expected to be J=2 (tensor)
    cout << "\nf2(1525) Analysis (Expected J=2):" << endl;
    cout << "Dominant coefficient should be a_0 and a_2 for J=2 particle" << endl;
    double max_coeff_1525 = 0;
    int dominant_J_1525 = 0;
    for (int J=0; J<=Jmax; J++) {
        double coeff = TMath::Abs(fLeg1525->GetParameter(J+1));
        if (coeff > max_coeff_1525) {
            max_coeff_1525 = coeff;
            dominant_J_1525 = J;
        }
    }
    cout << "Dominant term: a_" << dominant_J_1525 << " = " << fLeg1525->GetParameter(dominant_J_1525+1) << endl;
    
    // Ratio analysis
    cout << "\nCoefficient Ratios for Angular Distribution Analysis:" << endl;
    cout << "f0(1710) a_2/a_0 ratio: " << fLeg1710->GetParameter(3)/fLeg1710->GetParameter(1) << endl;
    cout << "f2(1525) a_2/a_0 ratio: " << fLeg1525->GetParameter(3)/fLeg1525->GetParameter(1) << endl;
    
    // For a J=2 particle, the angular distribution should have significant a_0 and a_2 terms
    // For a J=0 particle, only a_0 term should be significant

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
