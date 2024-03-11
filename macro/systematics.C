#include <iostream>
#include <cmath>
#include "src/style.h"
#include "src/fitfunc.h"

using namespace std;
// using namespace TMath;

void systematics()
{

    TCanvas *c = new TCanvas("c", "c", 1200, 900);
    gstyle(); // this is defined in header file
    SetCanvasStyle2(c, 0.2, 0.05, 0.08, 0.2);
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);

    // Default
    TFile *file1 = new TFile("/home/sawan/Documents/default/default_ME_pol2.root", "READ");
    if (!file1 || file1->IsZombie())
    {
        std::cout << "Error: Failed to open the default ROOT file!" << std::endl;
        // return 1;
    }

    // Group1 DCA
    TFile *file2 = new TFile("/home/sawan/Documents/DCA_variations/DCAz_ME_pol2.root", "READ");
    TFile *file3 = new TFile("/home/sawan/Documents/DCA_variations/DCAxy_ME_pol2.root", "READ");
    TFile *file4 = new TFile("/home/sawan/Documents/DCA_variations/DCAall_ME_pol2.root", "READ");
    if (!file2 || file2->IsZombie())
    {
        std::cout << "Error: Failed to open the group 1 DCA files!" << std::endl;
        // return 1;
    }

    // Group2 PID
    TFile *file5 = new TFile("/home/sawan/Documents/PID_variations/PID2_ME_pol2.root", "READ");
    TFile *file6 = new TFile("/home/sawan/Documents/PID_variations/PID2.5_ME_pol2.root", "READ");
    TFile *file7 = new TFile("/home/sawan/Documents/PID_variations/PID3.5_ME_pol2.root", "READ");
    TFile *file8 = new TFile("/home/sawan/Documents/PID_variations/PID4_ME_pol2.root", "READ");

    if (!file5 || file5->IsZombie())
    {
        std::cout << "Error: Failed to open the group 2 PID file!" << std::endl;
        // return 1;
    }

    // Group3 signal extraction range variations
    //  fit range variations
    TFile *file9 = new TFile("/home/sawan/Documents/signal_extraction/Fit_range_variations/fit_var_m20MeV_ME_pol2.root", "READ");
    TFile *file10 = new TFile("/home/sawan/Documents/signal_extraction/Fit_range_variations/fit_var_mp20MeV_ME_pol2.root", "READ");
    TFile *file11 = new TFile("/home/sawan/Documents/signal_extraction/Fit_range_variations/fit_var_pm20MeV_ME_pol2.root", "READ");
    TFile *file12 = new TFile("/home/sawan/Documents/signal_extraction/Fit_range_variations/fit_var_p20MeV_ME_pol2.root", "READ");

    if (!file9 || file9->IsZombie())
    {
        std::cout << "Error: failed to open the group 3 fit range variation file" << endl;
        // return 1;
    }

    // normalization range variations
    TFile *file13 = new TFile("/home/sawan/Documents/signal_extraction/normalization_range_variations/fit_var_m20MeV_ME_pol2.root", "READ");
    TFile *file14 = new TFile("/home/sawan/Documents/signal_extraction/normalization_range_variations/fit_var_mp20MeV_ME_pol2.root", "READ");
    TFile *file15 = new TFile("/home/sawan/Documents/signal_extraction/normalization_range_variations/fit_var_mp30MeV_ME_pol2.root", "READ");

    if (!file13 || file13->IsZombie())
    {
        std::cout << "Error: failed to open the group 3 normalization range variation file" << endl;
        // return 1;
    }

    // Group4 track selections
    TFile *file17 = new TFile("/home/sawan/Documents/track_selections/PVTrk/PVtrk_ME_pol2.root", "READ");

    if (!file17 || file17->IsZombie())
    {
        std::cout << "Error: failed to open the group 4 normalization range variation file" << endl;
        // return 1;
    }

    // Group5 signal extraction method variation
    // free width
    TFile *file18 = new TFile("/home/sawan/Documents/combinatorial_background/free_width/free_width_ME_pol2.root", "READ");

    // like sign
    TFile *file19 = new TFile("/home/sawan/Documents/combinatorial_background/like_sign/default_LS_pol2.root", "READ");

    // pol3
    TFile *file16 = new TFile("/home/sawan/Documents/signal_extraction/pol3/default_ME_pol3.root", "READ");

    const int npt = 13;
    int nvar;
    double x[npt] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.25, 2.75, 4};
    double xerr[npt] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 1};
    double NPT[npt + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 5.0};

    string var_names[5] = {"DCA variations", "PID variations", "Signal Extraction range variations", "Track selections variations", "signal extraction method variation"};

    TH1F *hfrac_total = new TH1F("hfrac_total", "fractional uncertainity", npt, NPT);
    TH1F *hfrac[5];
    for (int i = 0; i < 5; i++)
    {
        hfrac[i] = new TH1F(Form("hfrac%d", i), Form("hfrac_%d", i), npt, NPT);
    }
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

    TGraphErrors *denominator = (TGraphErrors *)file1->Get("yieldls");
    // denominator->Draw("ap");

    // Group 1 DCA
    TGraphErrors *numerator2 = (TGraphErrors *)file2->Get("yieldls");
    TGraphErrors *numerator3 = (TGraphErrors *)file3->Get("yieldls");
    TGraphErrors *numerator4 = (TGraphErrors *)file4->Get("yieldls");

    // Group 2 PID
    TGraphErrors *numerator5 = (TGraphErrors *)file5->Get("yieldls");
    TGraphErrors *numerator6 = (TGraphErrors *)file6->Get("yieldls");
    TGraphErrors *numerator7 = (TGraphErrors *)file7->Get("yieldls");
    TGraphErrors *numerator8 = (TGraphErrors *)file8->Get("yieldls");

    // Group 3 signal extraction range variations
    // fitrange
    TGraphErrors *numerator9 = (TGraphErrors *)file9->Get("yieldls");
    TGraphErrors *numerator10 = (TGraphErrors *)file10->Get("yieldls");
    TGraphErrors *numerator11 = (TGraphErrors *)file11->Get("yieldls");
    TGraphErrors *numerator12 = (TGraphErrors *)file12->Get("yieldls");
    // normalization range
    TGraphErrors *numerator13 = (TGraphErrors *)file13->Get("yieldls");
    TGraphErrors *numerator14 = (TGraphErrors *)file14->Get("yieldls");
    TGraphErrors *numerator15 = (TGraphErrors *)file15->Get("yieldls");

    // Group 4 track selections
    TGraphErrors *numerator17 = (TGraphErrors *)file17->Get("yieldls");

    // Group5 signal extraction method variation
    // free width
    TGraphErrors *numerator18 = (TGraphErrors *)file18->Get("yieldls");
    // like sign
    TGraphErrors *numerator19 = (TGraphErrors *)file19->Get("yieldls");
    // pol3
    TGraphErrors *numerator16 = (TGraphErrors *)file16->Get("yieldls");

    TGraphErrors *numeratorGraphs1[3] = {numerator2, numerator3, numerator4};
    TGraphErrors *numeratorGraphs2[4] = {numerator5, numerator6, numerator7, numerator8};
    TGraphErrors *numeratorGraphs3[7] = {numerator9, numerator10, numerator11, numerator12, numerator13, numerator14, numerator15};
    TGraphErrors *numeratorGraphs4[1] = {numerator17};
    TGraphErrors *numeratorGraphs5[3] = {numerator18, numerator19, numerator16};
    // TGraphErrors *numeratorGraphs5[1] = {numerator18};
    //  TGraphErrors *grratio[nvar];

    int anvar[5] = {3, 4, 7, 1, 3}; // no. of variations in a group

    // gPad->SetLogy(0);
    // gStyle->SetOptTitle(0);

    double Ydefault[npt] = {0};
    double Ydefaulterr[npt] = {0};
    double Yvariation[npt] = {0};
    // double variance[npt] = {0.0};
    double relative_error[5][npt] = {0};

    // for (int set = 0; set < 5; set++)  // changed here
    for (int set = 0; set < 4; set++)
    {
        double variance[npt] = {0.0};
        nvar = anvar[set];

        for (int j = 0; j < nvar; j++)
        {
            for (int i = 0; i < npt; i++)
            {
                if (set == 0)
                {
                    Yvariation[i] = numeratorGraphs1[j]->GetY()[i];
                }
                if (set == 1)
                {
                    Yvariation[i] = numeratorGraphs2[j]->GetY()[i];
                }
                if (set == 2)
                {
                    Yvariation[i] = numeratorGraphs3[j]->GetY()[i];
                    // continue;
                }
                if (set == 3)
                {
                    Yvariation[i] = numeratorGraphs4[j]->GetY()[i];
                }
                if (set == 4)
                {
                    Yvariation[i] = numeratorGraphs5[j]->GetY()[i];
                }

                Ydefault[i] = denominator->GetY()[i];
                // cout<<"The default value is "<<Ydefault[i]<<endl;
                Ydefaulterr[i] = denominator->GetErrorY(i);

                variance[i] += ((Ydefault[i] - Yvariation[i]) / Ydefault[i]) * ((Ydefault[i] - Yvariation[i]) / Ydefault[i]);
            }
        }

        for (int ipt = 0; ipt < npt; ipt++)
        {
            relative_error[set][ipt] = TMath::Sqrt(variance[ipt] / nvar);
            hfrac[set]->SetBinContent(ipt + 1, relative_error[set][ipt]);
            // cout << "The relative uncertainity for set"<<set<<" is "<<relative_error[set][ipt] << endl;
        }
        // cout<<"\n\n\n";
    }

    double systematics_square[npt] = {0};
    double systematics[npt] = {0};
    double frac_uncertainity[npt] = {0};
    for (int ipt = 0; ipt < npt; ipt++)
    {
        for (int set = 0; set < 4; set++)
        { // changed here
            // if(set != 2)
            systematics_square[ipt] += relative_error[set][ipt] * relative_error[set][ipt];
        }

        systematics[ipt] = (TMath::Sqrt(systematics_square[ipt])) * Ydefault[ipt];
        frac_uncertainity[ipt] = TMath::Sqrt(systematics_square[ipt]);
        hfrac_total->SetBinContent(ipt + 1, frac_uncertainity[ipt]);
        cout << frac_uncertainity[ipt] << endl;
    }

    TGraphErrors *gsys = new TGraphErrors(npt, x, Ydefault, xerr, systematics);
    SetGrapherrorStyle(gsys, 1, 1);
    gsys->SetMarkerSize(1.1);
    TGraphErrors *gyield = new TGraphErrors(npt, x, Ydefault, xerr, Ydefaulterr);
    SetGrapherrorStyle(gyield, 2, 2);
    gyield->SetMarkerSize(1.1);
    gsys->SetFillStyle(0);
    gPad->SetLogy();
    gsys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gsys->GetYaxis()->SetTitle("(1/N_{ev})* d^{2}N/(dp_{T} dy #epsilon BR) (Gev/c)^{-1}");
    gsys->Draw("ap2");
    gsys->Draw("psame");
    gyield->Draw("psame");
    // gyield->Draw("ap");
    c->SaveAs("/home/sawan/Documents/systematics.png");

    // plotting the fractional uncertainities.

    gPad->SetLogy(0);
    hfrac_total->SetLineColor(1);
    hfrac_total->SetLineStyle(2);
    hfrac_total->SetLineWidth(2);
    hfrac_total->GetYaxis()->SetRangeUser(0, 0.3);
    hfrac_total->Draw();
    leg->AddEntry(hfrac_total, "total systematic uncertainity");
    for (int i = 0; i < 4; i++) // changed here
    {
        hfrac[i]->SetLineColor(i + 2);
        hfrac[i]->SetLineWidth(2);
        hfrac[i]->Draw("same");
        leg->AddEntry(hfrac[i], var_names[i].c_str());
    }
    leg->Draw();

    c->SaveAs("/home/sawan/Documents/fraction_uncertainity.png");
}