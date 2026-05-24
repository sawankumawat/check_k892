#include "helper.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <TH1.h>
#include <TH1D.h>

HistogramOperations::HistogramOperations() {} // Default constructor implementation

TH1D *HistogramOperations::CalculateRatio(TH1D *defaulthist, TH1 *varhist)
{
        TH1D *hRatio = (TH1D *)varhist->Clone();
        hRatio->Divide(defaulthist);
        return hRatio;
}

TH1D *HistogramOperations::RelativeUncertainty(TH1D *defaulthist, std::vector<TH1 *> &variationHistograms)
{
        int vecsize = variationHistograms.size();
        TH1D *hrelunsq[vecsize];
        for (int i = 0; i < vecsize; i++)
        {
                hrelunsq[i] = (TH1D *)variationHistograms[i]->Clone();
                for (int j = 0; j < hrelunsq[i]->GetNbinsX(); j++)
                {
                        double def = defaulthist->GetBinContent(j + 1);
                        double var = variationHistograms[i]->GetBinContent(j + 1);
                        double relunsq = 0;
                        if (def != 0)
                        {
                                relunsq = pow(def - var, 2) / pow(def, 2);
                        }
                        else
                        {
                                relunsq = 0;
                        }
                        hrelunsq[i]->SetBinContent(j + 1, relunsq);
                }
                if (i != 0)
                {
                        hrelunsq[0]->Add(hrelunsq[i]);
                }
        }

        hrelunsq[0]->Scale(1.0 / vecsize);

        for (int i = 0; i < hrelunsq[0]->GetNbinsX(); i++)
        {
                hrelunsq[0]->SetBinContent(i + 1, sqrt(hrelunsq[0]->GetBinContent(i + 1)));
        }

        return hrelunsq[0];
}

TH1D *HistogramOperations::sigma(std::vector<TH1D *> &variationHistograms)
{
        int vecsize = variationHistograms.size();

        // --- 1. Handle empty vector edge case ---
        if (vecsize == 0)
                return nullptr;

        // --- 2. Handle single histogram case ---
        if (vecsize == 1)
        {
                return (TH1D *)variationHistograms[0]->Clone();
        }

        // --- 3. Existing logic for multiple histograms ---
        TH1D *hsigma[vecsize];
        for (int i = 0; i < vecsize; i++)
        {
                hsigma[i] = (TH1D *)variationHistograms[i]->Clone();
                for (int j = 0; j < hsigma[i]->GetNbinsX(); j++)
                {
                        double var = variationHistograms[i]->GetBinContent(j + 1);
                        hsigma[i]->SetBinContent(j + 1, var * var);
                }
                if (i != 0)
                {
                        hsigma[0]->Add(hsigma[i]);
                }
        }

        for (int i = 0; i < hsigma[0]->GetNbinsX(); i++)
        {
                hsigma[0]->SetBinContent(i + 1, sqrt(hsigma[0]->GetBinContent(i + 1)));
        }

        // Optional: clean up the cloned array pointers to prevent memory leaks
        for (int i = 1; i < vecsize; i++)
        {
                delete hsigma[i];
        }

        return hsigma[0];
}

TH1D *HistogramOperations::smooth(TH1D *hist1, int n = 2)
{
        TH1D *hsmooth = (TH1D *)hist1->Clone();
        for (int i = 0; i < n; i++)
        {
                for (int j = 1; j < hsmooth->GetNbinsX() - 1; j++)
                {
                        double bin1 = hsmooth->GetBinContent(j);
                        double bin2 = hsmooth->GetBinContent(j + 1);
                        double bin3 = hsmooth->GetBinContent(j + 2);
                        double avg = (bin1 + bin2 + bin3) / 3;
                        hsmooth->SetBinContent(j + 1, avg);
                }
        }
        return hsmooth;
}

TH1D *HistogramOperations::barlowcheck(TH1D *hdefault, TH1 *hvariation, bool &checkbar)
{
        bool barlow = false;
        bool barlowtemp[4];
        TH1D *hbarlow = new TH1D("", "", 50, -50, 50);
        for (int i = 0; i < hvariation->GetNbinsX(); i++)
        {
                double defvalue = hdefault->GetBinContent(i + 1);
                double varvalue = hvariation->GetBinContent(i + 1);
                double delta = (varvalue - defvalue);
                double defstaterror = hdefault->GetBinError(i + 1);
                double varstaterror = hvariation->GetBinError(i + 1);
                double statsigma = sqrt(std::abs(pow(defstaterror, 2) - pow(varstaterror, 2)));
                double n = delta / statsigma;
                // cout<<"\n\n n is "<<n<<endl;
                hbarlow->Fill(n);
        }
        double mean = hbarlow->GetMean();
        double rms = hbarlow->GetRMS();
        double countstotal = hbarlow->Integral();
        double countsn1 = hbarlow->Integral(hbarlow->FindBin(-1.0), hbarlow->FindBin(1.0)) / countstotal;
        double countsn2 = hbarlow->Integral(hbarlow->FindBin(-2.0), hbarlow->FindBin(2.0)) / countstotal;
        // cout << "Mean: " << mean << ", RMS: " << rms << ", |n|<1 fraction: " << countsn1 * 100 << "%, |n|<2 fraction: " << countsn2 * 100 << "%" << endl;
        int conditionmet = 0;
        if (abs(mean) < 0.1)
                conditionmet++;
        if (abs(rms) < 1.1)
                conditionmet++;
        if (countsn1 > 0.60)
                conditionmet++;
        if (countsn2 > 0.90)
                conditionmet++;

        if (conditionmet >= 3)
                barlow = true;
        else
                barlow = false;

        // cout << "Conditions met: " << conditionmet << endl;

        checkbar = barlow;
        hbarlow->GetXaxis()->SetRangeUser(-100, 100);

        return hbarlow;
}

// std::vector<TH1 *> HistogramOperations::CalculateRatio(TH1D *defaulthist, std::vector<TH1 *> &variationHistograms)
// {
//         int vecsize = variationHistograms.size();
//         TH1D *hRatio[vecsize];
//         for (int i = 0; i < vecsize; i++)
//         {
//                 hRatio[i] = (TH1D *)variationHistograms[i]->Clone();
//                 for (int j = 0; j < hRatio[i]->GetNbinsX(); j++)
//                 {
//                         double def = defaulthist->GetBinContent(j + 1);
//                         double var = variationHistograms[i]->GetBinContent(j + 1);
//                         double ratio = var / def;
//                         hRatio[i]->SetBinContent(j + 1, ratio);
//                 }
//         }

//         std::vector<TH1 *> hRatioVec;
//         for (int i = 0; i < vecsize; i++)
//         {
//                 hRatioVec.push_back(hRatio[i]);
//         }

//         return hRatioVec;
// }
