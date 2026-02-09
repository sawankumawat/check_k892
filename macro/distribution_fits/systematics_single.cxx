#include "../src/style.h"
using namespace std;

void systematics_single()
{
    bool isSimple = true;
    float bins[] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0};

    // This code will calculate the relative uncertainty in mass and raw yield for f0 and f2 resonance and also check Barlow for a single variation only.
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent";
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/systematic2022_new/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent";

    if (isSimple)
    {

        // TFile *fDefault = new TFile((path + "/FitParam.root").c_str(), "READ");
        // TFile *fVariation = new TFile((path + "/FitParam_TPCPID2.root").c_str(), "READ");

        TFile *fDefault = new TFile((path + "/mult_0-100/Spectra/spectra.root").c_str(), "READ");
        TFile *fVariation = new TFile((path + "/mult_0-100/Spectra/spectra_TPCMinCls100.root").c_str(), "READ");

        if (fDefault->IsZombie() || fVariation->IsZombie())
        {
            cout << "Error opening file" << endl;
            return;
        }

        // TH1F *hMass_f1525_Default = (TH1F *)fDefault->Get("Mult_0_100/hMass_1525");
        // TH1F *hMass_f1710_Default = (TH1F *)fDefault->Get("Mult_0_100/hMass_1710");
        // TH1F *hMass_f1525_Variation = (TH1F *)fVariation->Get("Mult_0_100/hMass_1525");
        // TH1F *hMass_f1710_Variation = (TH1F *)fVariation->Get("Mult_0_100/hMass_1710");

        TH1F *hMass_f1525_Default = (TH1F *)fDefault->Get("hMass_1525");
        TH1F *hMass_f1710_Default = (TH1F *)fDefault->Get("hMass_1710");
        TH1F *hMass_f1525_Variation = (TH1F *)fVariation->Get("hMass_1525");
        TH1F *hMass_f1710_Variation = (TH1F *)fVariation->Get("hMass_1710");

        // TH1F *hYield_f1525_Default = (TH1F *)fDefault->Get("Mult_0_100/hYield_1525");
        // TH1F *hYield_f1710_Default = (TH1F *)fDefault->Get("Mult_0_100/hYield_1710");
        // TH1F *hYield_f1525_Variation = (TH1F *)fVariation->Get("Mult_0_100/hYield_1525");
        // TH1F *hYield_f1710_Variation = (TH1F *)fVariation->Get("Mult_0_100/hYield_1710");

        TH1F *hYield_f1525_Default = (TH1F *)fDefault->Get("hYield1525Corrected");
        TH1F *hYield_f1710_Default = (TH1F *)fDefault->Get("hYield1710Corrected");
        TH1F *hYield_f1525_Variation = (TH1F *)fVariation->Get("hYield1525Corrected");
        TH1F *hYield_f1710_Variation = (TH1F *)fVariation->Get("hYield1710Corrected");

        if (hMass_f1525_Default == nullptr || hMass_f1710_Default == nullptr || hMass_f1525_Variation == nullptr || hMass_f1710_Variation == nullptr || hYield_f1525_Default == nullptr || hYield_f1710_Default == nullptr || hYield_f1525_Variation == nullptr || hYield_f1710_Variation == nullptr)
        {
            cout << "Error: One of the histograms not found!" << endl;
            return;
        }

        // writing a lambda function to calculate relative uncertainty and return it as a histogram
        auto calculateRelativeUncertainty = [=](TH1F *hDefault, TH1F *hVariation, const string &name) -> TH1F *
        {
            int totalBins = hDefault->GetNbinsX();
            TH1F *hRelativeUncertainty = new TH1F((name + "_RelativeUncertainty").c_str(), (name + " Relative Uncertainty").c_str(), 6, bins);
            for (int i = 1; i <= totalBins; i++)
            {
                double defaultValue = hDefault->GetBinContent(i);
                double variationValue = hVariation->GetBinContent(i);
                double relativeUncertainty = 0.0;
                if (defaultValue != 0)
                {
                    relativeUncertainty = fabs(variationValue - defaultValue) / fabs(defaultValue);
                }
                hRelativeUncertainty->SetBinContent(i, relativeUncertainty);
                // if (i == 2)
                {
                    cout << name << " - Default: " << defaultValue << ", Variation: " << variationValue << ", Relative Uncertainty: " << relativeUncertainty << endl;
                }
            }
            return hRelativeUncertainty;
        };
        TH1F *hRelUncMass_f1525 = calculateRelativeUncertainty(hMass_f1525_Default, hMass_f1525_Variation, "f1525_Mass");
        TH1F *hRelUncMass_f1710 = calculateRelativeUncertainty(hMass_f1710_Default, hMass_f1710_Variation, "f1710_Mass");
        TH1F *hRelUncYield_f1525 = calculateRelativeUncertainty(hYield_f1525_Default, hYield_f1525_Variation, "f1525_Yield");
        TH1F *hRelUncYield_f1710 = calculateRelativeUncertainty(hYield_f1710_Default, hYield_f1710_Variation, "f1710_Yield");

        TCanvas *cAll = new TCanvas("cAll", "Systematic Uncertainties", 1200, 800);
        cAll->Divide(2, 2);
        SetCanvasStyle(cAll, 0.15, 0.03, 0.05, 0.14);
        TLatex lat;
        lat.SetNDC();
        lat.SetTextFont(42);
        lat.SetTextSize(0.06);
        lat.SetTextColor(kRed);
        cAll->cd(1);
        SetHistoQA(hRelUncMass_f1525);
        hRelUncMass_f1525->SetLineWidth(3);
        hRelUncMass_f1525->Draw("hist");
        lat.DrawLatex(0.2, 0.85, "Mass f'_{2}(1525)");
        cAll->cd(2);
        SetHistoQA(hRelUncMass_f1710);
        hRelUncMass_f1710->SetLineWidth(3);
        hRelUncMass_f1710->Draw("hist");
        lat.DrawLatex(0.2, 0.85, "Mass f_{0}(1710)");
        cAll->cd(3);
        SetHistoQA(hRelUncYield_f1525);
        hRelUncYield_f1525->SetLineWidth(3);
        hRelUncYield_f1525->Draw("hist");
        lat.DrawLatex(0.2, 0.85, "Raw Yield f'_{2}(1525)");
        cAll->cd(4);
        SetHistoQA(hRelUncYield_f1710);
        hRelUncYield_f1710->SetLineWidth(3);
        hRelUncYield_f1710->Draw("hist");
        lat.DrawLatex(0.2, 0.85, "Raw Yield f_{0}(1710)");
    }
    else
    {
        // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent";
        path = path + "/mult_0-100/Spectra/";
        TFile *fDefault = new TFile((path + "/spectra.root").c_str(), "READ");
        string variations[] = {"_TPCPID2", "_TPCPID5"};
        int totalVariations = sizeof(variations) / sizeof(variations[0]);
        TFile *fVariations[totalVariations];

        // TH1F *hMass_f1525_Default = (TH1F *)fDefault->Get("Mult_0_100/hMass_1525");
        // TH1F *hMass_f1710_Default = (TH1F *)fDefault->Get("Mult_0_100/hMass_1710");
        // TH1F *hYield_f1525_Default = (TH1F *)fDefault->Get("Mult_0_100/hYield_1525");
        // TH1F *hYield_f1710_Default = (TH1F *)fDefault->Get("Mult_0_100/hYield_1710");

        TH1F *hMass_f1525_Default = (TH1F *)fDefault->Get("hMass_1525");
        TH1F *hMass_f1710_Default = (TH1F *)fDefault->Get("hMass_1710");
        TH1F *hYield_f1525_Default = (TH1F *)fDefault->Get("hYield1525Corrected");
        TH1F *hYield_f1710_Default = (TH1F *)fDefault->Get("hYield1710Corrected");

        TH1F *hMass_f1525_Variation[totalVariations];
        TH1F *hMass_f1710_Variation[totalVariations];
        TH1F *hYield_f1525_Variation[totalVariations];
        TH1F *hYield_f1710_Variation[totalVariations];
        TH1F *hRelUncMass_f1525[totalVariations];
        TH1F *hRelUncMass_f1710[totalVariations];
        TH1F *hRelUncYield_f1525[totalVariations];
        TH1F *hRelUncYield_f1710[totalVariations];

        // writing a lambda function to calculate relative uncertainty and return it as a histogram
        auto calculateRelativeUncertainty = [=](TH1F *hDefault, TH1F *hVariation, const string &name) -> TH1F *
        {
            int totalBins = hDefault->GetNbinsX();
            TH1F *hRelativeUncertainty = new TH1F((name + "_RelativeUncertainty").c_str(), (name + " Relative Uncertainty").c_str(), 6, bins);
            for (int i = 1; i <= totalBins; i++)
            {
                double defaultValue = hDefault->GetBinContent(i);
                double variationValue = hVariation->GetBinContent(i);
                double relativeUncertainty = 0.0;
                if (defaultValue != 0)
                {
                    relativeUncertainty = fabs(variationValue - defaultValue) / fabs(defaultValue);
                }
                hRelativeUncertainty->SetBinContent(i, relativeUncertainty);
                if (i == 1)
                {
                    cout << name << " - Default: " << defaultValue << ", Variation: " << variationValue << ", Relative Uncertainty: " << relativeUncertainty << endl;
                }
            }
            return hRelativeUncertainty;
        };
        for (int i = 0; i < totalVariations; i++)
        {
            // fVariations[i] = new TFile((path + "/FitParam" + variations[i] + ".root").c_str(), "READ");
            fVariations[i] = new TFile((path + "/spectra" + variations[i] + ".root").c_str(), "READ");
            if (fVariations[i]->IsZombie())
            {
                cout << "Error opening file for variation: " << variations[i] << endl;
                return;
            }
            // hMass_f1525_Variation[i] = (TH1F *)fVariations[i]->Get("Mult_0_100/hMass_1525");
            // hMass_f1710_Variation[i] = (TH1F *)fVariations[i]->Get("Mult_0_100/hMass_1710");
            // hYield_f1525_Variation[i] = (TH1F *)fVariations[i]->Get("Mult_0_100/hYield_1525");
            // hYield_f1710_Variation[i] = (TH1F *)fVariations[i]->Get("Mult_0_100/hYield_1710");

            hMass_f1525_Variation[i] = (TH1F *)fVariations[i]->Get("hMass_1525");
            hMass_f1710_Variation[i] = (TH1F *)fVariations[i]->Get("hMass_1710");
            hYield_f1525_Variation[i] = (TH1F *)fVariations[i]->Get("hYield1525Corrected");
            hYield_f1710_Variation[i] = (TH1F *)fVariations[i]->Get("hYield1710Corrected");

            if (hMass_f1525_Variation[i] == nullptr || hMass_f1710_Variation[i] == nullptr || hYield_f1525_Variation[i] == nullptr || hYield_f1710_Variation[i] == nullptr)
            {
                cout << "Error: One of the histograms not found for variation: " << variations[i] << endl;
                return;
            }
            hRelUncMass_f1525[i] = calculateRelativeUncertainty(hMass_f1525_Default, hMass_f1525_Variation[i], "f1525_Mass" + variations[i]);
            hRelUncMass_f1710[i] = calculateRelativeUncertainty(hMass_f1710_Default, hMass_f1710_Variation[i], "f1710_Mass" + variations[i]);
            hRelUncYield_f1525[i] = calculateRelativeUncertainty(hYield_f1525_Default, hYield_f1525_Variation[i], "f1525_Yield" + variations[i]);
            hRelUncYield_f1710[i] = calculateRelativeUncertainty(hYield_f1710_Default, hYield_f1710_Variation[i], "f1710_Yield" + variations[i]);
        }

        // Lambda function for quadrature sum of all relative uncertainties
        auto quadratureSum = [=](TH1F **hRelUncs, const string &name) -> TH1F *
        {
            int totalBins = hRelUncs[0]->GetNbinsX();
            TH1F *hQuadratureSum = new TH1F((name + "_QuadratureSum").c_str(), (name + " Quadrature Sum").c_str(), 6, bins);
            for (int i = 1; i <= totalBins; i++)
            {
                double sumSquares = 0.0;
                for (int j = 0; j < totalVariations; j++)
                {
                    double relUnc = hRelUncs[j]->GetBinContent(i);
                    sumSquares += relUnc * relUnc;
                }
                sumSquares = sumSquares / totalVariations; // taking average of squares
                double quadratureUnc = sqrt(sumSquares);
                hQuadratureSum->SetBinContent(i, quadratureUnc);
            }
            return hQuadratureSum;
        };
        TH1F *hQuadratureSumMass_f1525 = quadratureSum(hRelUncMass_f1525, "Mass_f1525");
        TH1F *hQuadratureSumMass_f1710 = quadratureSum(hRelUncMass_f1710, "Mass_f1710");
        TH1F *hQuadratureSumYield_f1525 = quadratureSum(hRelUncYield_f1525, "Yield_f1525");
        TH1F *hQuadratureSumYield_f1710 = quadratureSum(hRelUncYield_f1710, "Yield_f1710");

        TCanvas *cAll = new TCanvas("cAll", "Systematic Uncertainties", 1200, 800);
        cAll->Divide(2, 2);
        SetCanvasStyle(cAll, 0.15, 0.03, 0.05, 0.14);
        TLatex lat;
        lat.SetNDC();
        lat.SetTextFont(42);
        lat.SetTextSize(0.06);
        lat.SetTextColor(kRed);
        cAll->cd(1);
        SetHistoQA(hQuadratureSumMass_f1525);
        hQuadratureSumMass_f1525->SetLineWidth(3);
        hQuadratureSumMass_f1525->Draw("hist");
        lat.DrawLatex(0.2, 0.85, "Mass f'_{2}(1525)");
        cAll->cd(2);
        SetHistoQA(hQuadratureSumMass_f1710);
        hQuadratureSumMass_f1710->SetLineWidth(3);
        hQuadratureSumMass_f1710->Draw("hist");
        lat.DrawLatex(0.2, 0.85, "Mass f_{0}(1710)");
        cAll->cd(3);
        SetHistoQA(hQuadratureSumYield_f1525);
        hQuadratureSumYield_f1525->SetLineWidth(3);
        hQuadratureSumYield_f1525->Draw("hist");
        lat.DrawLatex(0.2, 0.85, "Raw Yield f'_{2}(1525)");
        cAll->cd(4);
        SetHistoQA(hQuadratureSumYield_f1710);
        hQuadratureSumYield_f1710->SetLineWidth(3);
        hQuadratureSumYield_f1710->Draw("hist");
        lat.DrawLatex(0.2, 0.85, "Raw Yield f_{0}(1710)");
    }
}