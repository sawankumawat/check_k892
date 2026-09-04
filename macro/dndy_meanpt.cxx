#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include "src/style.h"
#include "src/initializations.h"
#include "SystematicHelper/helper.cxx"

TFile *openTFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (f->IsZombie())
    {
        cout << "Error: File not found: " << path << endl;
        return nullptr;
    }
    return f;
}
TH1D *GetHisto(TFile *f, const string &name);

void dndy_meanpt()
{
    string basePath = "../output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/";
    TFile *fTotalSys = openTFile(basePath + "SystematicsPlots/SysUncert.root");
    TFile *fUncorrSys = openTFile(basePath + "SystematicsPlots/UnCorrSystematics.root");

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    int totalMultBins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    for (int imult = 0; imult < totalMultBins; imult++)
    {
        int multLow = mult_classes[imult];
        int multHigh = mult_classes[imult + 1];

        TH1D *hTotalSys = GetHisto(fTotalSys, Form("hTotalSysSmoothed_%d_%d", multLow, multHigh));
        TH1D *hUncorrSys = GetHisto(fUncorrSys, Form("hUncorrelatedUncertainty_mult_%d_%d", multLow, multHigh));
        TFile *fSpectra = openTFile(basePath + Form("ROTATED/corrected_spectra_%d_%d.root", multLow, multHigh));
        TH1D *hSpectra = GetHisto(fSpectra, Form("mult_%d-%d/corrected_spectra_Integral_final", multLow, multHigh));

        double YieldTimespT = 0.0;
        double TotalYield = 0.0;
        double dndy_uncorr = 0.0;
        double meanpt_unc2 = 0.0;
        double meanpt = 0.0;

        for (int ibin = 1; ibin <= hSpectra->GetNbinsX(); ibin++)
        {
            double yield = hSpectra->GetBinContent(ibin);
            double pt = hSpectra->GetBinCenter(ibin);
            double ptBinWidth = hSpectra->GetBinWidth(ibin);

            TotalYield += yield * ptBinWidth;
            YieldTimespT += pt * yield * ptBinWidth;
        }

        if (TotalYield > 0)
            meanpt = YieldTimespT / TotalYield;

        // propagate yield-uncorrelated uncertainty

        for (int i = 1; i <= hSpectra->GetNbinsX(); i++)
        {
            // Uncorrelated uncertainty
            double relUnc = hUncorrSys->GetBinContent(i);

            // corrected yield
            double yield = hSpectra->GetBinContent(i);
            double pt = hSpectra->GetBinCenter(i);
            double ptBinWidth = hSpectra->GetBinWidth(i);

            // absolute uncertainty on corrected yield
            double absUnc = relUnc * yield;

            // dN/dy uncertainty (delta_dNdy = delta_yield * ptBinWidth + yield * delta_ptBinWidth, but delta_ptBinWidth is 0, so delta_dNdy = YieldUncertainty * ptBinWidth)
            dndy_uncorr += TMath::Power(absUnc * ptBinWidth, 2);

            // <pT> uncertainty (Note that in <pt> error calculation in the formula A/B both A and B are correlated, so there is minus sign in the formula. It is calculated as error = (deltaA * B - A * deltaB)/B^2, where A = sum(pt * yield * ptBinWidth), B = sum(yield * ptBinWidth))
            double contrib = absUnc * ptBinWidth * (pt - meanpt) / TotalYield;
            meanpt_unc2 += contrib * contrib;
        }

        double dn_uncorr = TMath::Sqrt(dndy_uncorr);
        double meanpt_unc = TMath::Sqrt(meanpt_unc2);

        cout << "\nRESULTS FOR yield " << multLow << "-" << multHigh << endl;
        cout << "dN/dy value and uncorrelated uncertainty = " << TotalYield << " ± " << dn_uncorr << endl;
        cout << "<pT> value and uncorrelated uncertainty = " << meanpt << " ± " << meanpt_unc << endl;
    }
}

TH1D *GetHisto(TFile *f, const string &name)
{
    TH1D *histo = (TH1D *)f->Get(name.c_str());

    if (!histo || histo == nullptr)
    {
        cout << "Error: histo " << name << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetHistoQA(histo);
    histo->SetTitle(0);
    return histo;
}