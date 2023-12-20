#include <iostream>
#include <cmath>
#include <string>
#include "src/style.h"
#include "src/fitfunc.h"

using namespace std;
// using namespace TMath;

void ratios()
{

    TCanvas *c = new TCanvas("c", "c", 1200, 900);
    // gstyle(); // this is defined in header file
    SetCanvasStyle2(c, 0.2, 0.05, 0.08, 0.2);
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);
    TLegend *leg = new TLegend(0.5, 0.6, 0.9, 0.9);

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

    // Group3 signal extraction
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

    if (!file15 || file15->IsZombie())
    {
        std::cout << "Error: failed to open the group 3 normalization range variation file" << endl;
        // return 1;
    }

    // pol3
    TFile *file16 = new TFile("/home/sawan/Documents/signal_extraction/pol3/default_ME_pol3.root", "READ");

    // Group4 track selections
    TFile *file17 = new TFile("/home/sawan/Documents/track_selections/PVTrk/PVtrk_ME_pol2.root", "READ");

    if (!file17 || file17->IsZombie())
    {
        std::cout << "Error: failed to open the group 4 normalization range variation file" << endl;
        // return 1;
    }

    // Group5 combinatorial background
    // free width
    TFile *file18 = new TFile("/home/sawan/Documents/combinatorial_background/free_width/free_width_ME_pol2.root", "READ");

    // like sign
    TFile *file19 = new TFile("/home/sawan/Documents/combinatorial_background/like_sign/default_LS_pol2.root", "READ");

    string file_names[18] = {"DCAz loose", "DCAxy loose", "DCA loose", "PID2", "PID2.5", "PID3.5", "PID4", "fitting range plus 20 MeV", "fitting range wider by 40 MeV", "fitting range narrow by 40 MeV", "fitting range minus 20 MeV", "norm range minus 20 MeV", "norm range wider by 20 MeV", "norm range wider by 30 MeV", "pol3 bkg fit", "PVtrk", "free width", "like sign"};

    const int Npt = 13;
    int nvar;
    double x[Npt] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.25, 2.75, 4};
    double xerr[Npt] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 1};

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

    // Group 3 signal extraction range variation
    // fitrange
    TGraphErrors *numerator9 = (TGraphErrors *)file9->Get("yieldls");
    TGraphErrors *numerator10 = (TGraphErrors *)file10->Get("yieldls");
    TGraphErrors *numerator11 = (TGraphErrors *)file11->Get("yieldls");
    TGraphErrors *numerator12 = (TGraphErrors *)file12->Get("yieldls");
    // normalization range
    TGraphErrors *numerator13 = (TGraphErrors *)file13->Get("yieldls");
    TGraphErrors *numerator14 = (TGraphErrors *)file14->Get("yieldls");
    TGraphErrors *numerator15 = (TGraphErrors *)file15->Get("yieldls");
    // pol3
    // numerator15->Draw("ap");
    TGraphErrors *numerator16 = (TGraphErrors *)file16->Get("yieldls");

    // Group 4 track selections
    TGraphErrors *numerator17 = (TGraphErrors *)file17->Get("yieldls");

    // Group5 signal extraction method variation
    TGraphErrors *numerator18 = (TGraphErrors *)file18->Get("yieldls");
    TGraphErrors *numerator19 = (TGraphErrors *)file19->Get("yieldls");

    TGraphErrors *array[18] = {numerator2, numerator3, numerator4, numerator5, numerator6, numerator7, numerator8, numerator9, numerator10, numerator11, numerator12, numerator13, numerator14, numerator15, numerator16, numerator17, numerator18, numerator19};

    for (int i = 0; i < 18; i++)
    {
        double ratio[Npt] = {0};
        double ratioerr[Npt] = {0};

        for (int ipt = 0; ipt < Npt; ipt++)
        {
            ratio[ipt] = array[i]->GetY()[ipt] / denominator->GetY()[ipt];
            // cout<<ratio[ipt]<<"   ";

            ratioerr[ipt] = ratio[ipt] * (TMath::Sqrt((array[i]->GetErrorY(ipt) / array[i]->GetY()[ipt]) * (array[i]->GetErrorY(ipt) / array[i]->GetY()[ipt]) + (denominator->GetErrorY(ipt) / denominator->GetY()[ipt]) * (denominator->GetErrorY(ipt) / denominator->GetY()[ipt])));
            // cout<<ratioerr[ipt]<< "   ";
        }
        // Create two pads on the canvas: one for each graph
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1);
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.5);

        // Divide the canvas into two pads vertically (one below the other)
        pad1->SetBottomMargin(0); // Remove the bottom margin for the first pad
        pad1->Draw();
        pad2->SetTopMargin(0); // Remove the top margin for the second pad
        pad2->Draw();

        // Set the first pad as active
        pad1->cd();
        pad1->SetLogy();
        denominator->SetLineColor(1);
        denominator->SetMarkerColor(1);
        denominator->SetTitle(0);
        array[i]->SetLineColor(2);
        array[i]->SetMarkerColor(2);
        array[i]->SetMarkerStyle(33);
        array[i]->SetTitle(0);
        denominator->GetYaxis()->CenterTitle(1);
        denominator->GetYaxis()->SetTitleSize(0.05);
        denominator->GetYaxis()->SetTitleOffset(0.9);
        denominator->GetYaxis()->SetTitle("Corrected p_{T} Spectra");
        leg->AddEntry(denominator, "default");
        leg->AddEntry(array[i], file_names[i].c_str());
        denominator->Draw("ap");
        array[i]->Draw("psame");
        leg->Draw();

        pad2->cd();

        TGraphErrors *gratio = new TGraphErrors(Npt, x, ratio, xerr, ratioerr);
        // SetGrapherrorStyle(gratio, 1, 1);
        gratio->SetTitle(0);
        // gratio->GetYaxis()->SetRangeUser(0.6, 1.4);
        pad2->SetGrid(1,1);
        gratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        gratio->GetYaxis()->CenterTitle(1);
        gratio->GetYaxis()->SetTitleSize(0.05);
        gratio->GetYaxis()->SetTitle("Ratio (Variation / Default)");
        gratio->SetMarkerStyle(8);
        gratio->SetMarkerSize(1.5);
        gratio->Draw("ap");
        c->SaveAs(Form("/home/sawan/Documents/ratio%d.png", i));
        c->Clear();
        leg->Clear();
    }
}