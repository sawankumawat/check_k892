#include <iostream>
#include <cmath>
#include "src/style.h"
#include "src/fitfunc.h"

using namespace std;

void compare()
{

    TCanvas *c = new TCanvas("c", "c", 1200, 900);
    gstyle(); // this is defined in header file
    SetCanvasStyle2(c, 0.2, 0.05, 0.08, 0.2);
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);

    TFile *file1 = new TFile("/home/sawan/Documents/default_ME_pol2.root", "READ");
    TFile *file2 = new TFile("/home/sawan/Documents/DCAz_ME_pol2.root", "READ");
    TFile *file3 = new TFile("/home/sawan/Documents/DCAxy_ME_pol2.root", "READ");
    TFile *file4 = new TFile("/home/sawan/Documents/DCAall_ME_pol2.root", "READ");

    // TFile *file1 = new TFile("/home/sawan/Documents/default_ME_pol2.root", "READ");
    // TFile *file2 = new TFile("/home/sawan/Documents/PID2_ME_pol2.root", "READ");
    // TFile *file3 = new TFile("/home/sawan/Documents/PID2.5_ME_pol2.root", "READ");
    // TFile *file4 = new TFile("/home/sawan/Documents/PID3.5_ME_pol2.root", "READ");
    // TFile *file5 = new TFile("/home/sawan/Documents/PID4_ME_pol2.root", "READ");


    // mass vs pT

    // TGraphErrors *massls = (TGraphErrors *)file1->Get("massls");
    // TGraphErrors *massmix = (TGraphErrors *)file2->Get("massmix");
    // SetGrapherrorStyle(massls, 1, 1);
    // SetGrapherrorStyle(massmix, 2, 2);
    // massls->Draw("ap");
    // massmix->Draw("psame");
    // TLegend *massleg = new TLegend(0.65, 0.45, 0.9, 0.65);
    // massleg->SetFillColor(0);
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    // TLine *line = new TLine(massls->GetXaxis()->GetXmin(), 0.895, massls->GetXaxis()->GetXmax(), 0.895);
    // SetLineStyle(line, 4);
    // massleg->AddEntry(line, "PDG Mass", "l");
    // massleg->AddEntry(massls, "Like-Sign", "ep");
    // massleg->AddEntry(massmix, "Mixed-Event", "ep");
    // massleg->Draw("l");
    // c->SaveAs("../mass.png");

    // // width vs pT

    // TGraphErrors *widthls = (TGraphErrors *)file1->Get("widthls");
    // TGraphErrors *widthmix = (TGraphErrors *)file2->Get("widthmix");
    // SetGrapherrorStyle(widthls, 1, 1);
    // SetGrapherrorStyle(widthmix, 2, 2);
    // widthls->Draw("ap");
    // widthmix->Draw("psame");
    // TLegend *widthleg = new TLegend(0.45, 0.70, 0.7, 0.85);
    // widthleg->SetFillColor(0);
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    // TLine *line2 = new TLine(widthls->GetXaxis()->GetXmin(), 0.047, widthls->GetXaxis()->GetXmax(), 0.047);
    // SetLineStyle(line2, 4);
    // widthleg->AddEntry(line, "PDG Width", "l");
    // widthleg->AddEntry(widthls, "Like-Sign", "ep");
    // widthleg->AddEntry(widthmix, "Mixed-Event", "ep");
    // widthleg->Draw("l");
    // c->SaveAs("../width.png");

    // // chisq/NPT vs pT

    // TGraphErrors *chip2 = (TGraphErrors *)file5->Get("chils");
    // TGraphErrors *chip3 = (TGraphErrors *)file6->Get("chils");
    // SetGrapherrorStyle(chip2, 1, 1);
    // SetGrapherrorStyle(chip3, 2, 2);
    // gStyle->SetOptTitle(1);
    // chip2->SetTitle("Pol2 vs Pol3 fitting for 3 TOF cut");
    // chip2->Draw("ap");
    // chip3->Draw("psame");
    // TLegend *ly = new TLegend(0.65, 0.66, 0.87, 0.85);
    // ly->AddEntry(chip2, "pol2", "p");
    // ly->AddEntry(chip3, "pol3", "p");
    // ly->Draw();
    // c->SaveAs("/home/sawan/Aplots/chi_tof3.png");

    // // yield vs pT

    // TGraphErrors *yield1 = (TGraphErrors *)file1->Get("yieldls");
    // TGraphErrors *yield2 = (TGraphErrors *)file2->Get("yieldls");
    // // TGraphErrors *yieldp2_3 = (TGraphErrors *)file5->Get("yieldls");
    // SetGrapherrorStyle(yield1, 1, 1);
    // SetGrapherrorStyle(yield2, 2, 2);
    // // SetGrapherrorStyle(yieldp2_3, 4, 4);
    // gPad->SetLogy();
    // gStyle->SetOptTitle(1);
    // yield1->SetTitle("Pol2 fitting");
    // yield1->SetMarkerSize(1.2);
    // yield2->SetMarkerSize(1.2);
    // // yieldp2_3->SetMarkerSize(1.1);
    // yield1->Draw("ap");
    // yield2->Draw("psame");
    // // yieldp2_3->Draw("psame");
    // TLegend *ly = new TLegend(0.65, 0.66, 0.87, 0.85);
    // ly->AddEntry(yield1, "PID 3.5", "ep");
    // ly->AddEntry(yield2, "PID default", "ep");
    // // ly->AddEntry(yieldp2_3, "3 kaon TOF PID cut", "ep");
    // ly->Draw();
    // c->SaveAs("/home/sawan/Documents/yield_ME_PID_3.5vsdefault.png");

    // // yield ratio vs pT

    // TGraphErrors *ratiols = (TGraphErrors *)file1->Get("ratiols");
    // TGraphErrors *ratiomix = (TGraphErrors *)file2->Get("ratiomix");
    // SetGrapherrorStyle(ratiols, 1, 1);
    // SetGrapherrorStyle(ratiomix, 2, 2);
    // gPad->SetLogy(0);
    // ratiols->Draw("ap");
    // ratiomix->Draw("psame");
    // TLegend *lra = new TLegend(0.50, 0.7, 0.72, 0.9);
    // lra->AddEntry(ratiols, "Like-Sign", "p");
    // lra->AddEntry(ratiomix, "Mixed-Event", "p");
    // lra->Draw();
    // c->SaveAs("../ratio.png");

    // // significance vs pT
    // TH1F *significance2_2 = (TH1F *)file1->Get("significance");
    // TH1F *significance2_25 = (TH1F *)file3->Get("significance");
    // TH1F *significance2_3 = (TH1F *)file5->Get("significance");
    // significance2_2->SetLineColor(1);
    // significance2_25->SetLineColor(2);
    // significance2_3->SetLineColor(4);
    //  significance2_2->SetLineWidth(2);
    // significance2_25->SetLineWidth(2);
    // significance2_3->SetLineWidth(2);
    // gStyle->SetOptTitle(1);
    // significance2_2->SetTitle("Pol2 fitting");
    // significance2_2->SetMarkerSize(1.1);
    // significance2_25->SetMarkerSize(1.1);
    // significance2_3->SetMarkerSize(1.1);
    // significance2_2->Draw();
    // significance2_25->Draw("same");
    // significance2_3->Draw("same");
    // TLegend *ly = new TLegend(0.65, 0.66, 0.87, 0.85);
    // ly->AddEntry(significance2_2, "2 kaon TOF PID cut", "l");
    // ly->AddEntry(significance2_25, "2.5 kaon TOF PID cut", "l");
    // ly->AddEntry(significance2_3, "3 kaon TOF PID cut", "l");
    // ly->Draw();
    // c->SaveAs("/home/sawan/Aplots/significance_p2.png");

    // ratio of yields
    //  TFile *file = new TFile("/home/sawan/Documents/ratioDCAz.root", "RECREATE");

    const int npt = 13;
    double x[npt] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.25, 2.75, 5.5};

    TGraphErrors *denominator = (TGraphErrors *)file1->Get("yieldls");
    TGraphErrors *numerator1 = (TGraphErrors *)file2->Get("yieldls");
    TGraphErrors *numerator2 = (TGraphErrors *)file3->Get("yieldls");
    TGraphErrors *numerator3 = (TGraphErrors *)file4->Get("yieldls");
    // TGraphErrors *numerator4 = (TGraphErrors *)file5->Get("yieldls");


    TGraphErrors *numeratorGraphs[3] = {numerator1, numerator2, numerator3};
    // TGraphErrors *numeratorGraphs[4] = {numerator1, numerator2, numerator3, numerator4};
     TGraphErrors *grratio[3];

    gPad->SetLogy(0);
    gStyle->SetOptTitle(0);

    double yieldDCA_all[npt];
    double yieldDCA_all_err[npt];
    double yieldDCA_z[npt];
    double yieldDCA_zerr[npt];
    double ratioy[npt];
    double ratioyerr[npt];
    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < npt; i++)
        {
            yieldDCA_all[i] = numeratorGraphs[j]->GetY()[i];
            yieldDCA_z[i] = denominator->GetY()[i];
            yieldDCA_all_err[i] = numeratorGraphs[j]->GetErrorY(i);
            yieldDCA_zerr[i] = denominator->GetErrorY(i);
            ratioy[i] = yieldDCA_all[i] / yieldDCA_z[i];
            ratioyerr[i] = ratioy[i] * TMath::Sqrt((yieldDCA_zerr[i] * yieldDCA_zerr[i]) / (yieldDCA_z[i] * yieldDCA_z[i]) + (yieldDCA_all_err[i] * yieldDCA_all_err[i]) / (yieldDCA_all[i] * yieldDCA_all[i])); // error propagation
        }

        grratio[j]= new TGraphErrors(npt, x, ratioy, 0, ratioyerr);
        grratio[j]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        grratio[j]->GetYaxis()->SetTitle("Yield Ratio");
        grratio[j]->GetYaxis()->SetRangeUser(0.75, 1.28);
        // grratio[j]->GetYaxis()->SetRangeUser(0.5, 1.5);

        if(j != 2)
        SetGrapherrorStyle(grratio[j], j +1, j +1);
        else
        SetGrapherrorStyle(grratio[j], j +4, j +4);


        if (j == 0)
            grratio[j]->Draw("ap");
        else
            grratio[j]->Draw("p same");
    }

    TLegend *legend = new TLegend(0.65, 0.65, 0.85, 0.85); 
    legend->AddEntry(grratio[0], "DCAz_tight / default", "ep");
    legend->AddEntry(grratio[1], "DCAxy_tight / default", "ep");
    legend->AddEntry(grratio[2], "DCAall_tight / default", "ep");
    // legend->AddEntry(grratio[3], "PID4 / default", "ep");
    legend->Draw();

    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    c->SaveAs("/home/sawan/Documents/ratio_yield_DCA.png");
}