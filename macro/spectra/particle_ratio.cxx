#include <iostream>
#include "style.h"
#include "fitfunc.h"
#include <cmath>
using namespace std;

void graphstyle(TGraphAsymmErrors *gr)
{
    gr->SetTitle(0);
    gr->SetMarkerStyle(24);
    gr->SetMarkerColor(1);
    gr->SetMarkerSize(1.8);
    gr->SetLineWidth(2);
    gr->SetLineColor(1);
    gr->GetXaxis()->SetTitle("dN/d#eta");
    gr->GetXaxis()->CenterTitle(true);
    // gr->GetXaxis()->SetNdivisions(506);
    // gr->GetYaxis()->SetNdivisions(505);
    gr->GetXaxis()->SetLabelOffset(0.015);
    gr->GetXaxis()->SetLabelFont(42);
    gr->GetXaxis()->SetTitleFont(42);
    gr->GetXaxis()->SetLabelSize(0.05);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetTickLength(0.04);
    gr->GetXaxis()->SetTitleOffset(1.25);
    gr->GetYaxis()->SetTitleOffset(1.45);
    gr->GetYaxis()->CenterTitle(true);
    gr->GetYaxis()->SetDecimals(false);
    // gr->GetYaxis()->SetNdivisions(310);
    gr->GetYaxis()->SetLabelOffset(0.015);
    gr->GetYaxis()->SetLabelFont(42);
    gr->GetYaxis()->SetLabelSize(0.05);
    gr->GetYaxis()->SetTickLength(0.04);
    gr->GetYaxis()->SetTitleSize(0.05);
    // gr->GetYaxis()->SetTitle("");
    gr->GetYaxis()->SetTitleFont(42);
}

void particle_ratio()
{
    const int point = 1;
    double dNdYSawan = 0.0645353;
    double dNdYSawansys = 0.00580818;
    double dNdYSawanstat = 6.97089e-05;
    double dNdEtaSawan = 2.94;
    double dndeta_error[] = {0.05, 0.11};

    double dndy_pubkaon = 0.1825;
    double dndy_pubkaon_sys = 0.015;
    double dndy_pubkaon_stat = 0.004;

    double dndy_pubpion = 1.489;
    double dndy_pubpion_sys = 0.074;
    double dndy_pubpion_stat = 0.004;

    double ratio_kaon = dNdYSawan / dndy_pubkaon;
    double ratio_kaon_syserr = sqrt(pow(dNdYSawansys / dNdYSawan, 2) + pow(dndy_pubkaon_sys / dndy_pubkaon, 2));
    double ratio_kaon_staterr = sqrt(pow(dNdYSawanstat / dNdYSawan, 2) + pow(dndy_pubkaon_stat / dndy_pubkaon, 2));

    double ratio_pion = dNdYSawan / dndy_pubpion;
    double ratio_pion_syserr = sqrt(pow(dNdYSawansys / dNdYSawan, 2) + pow(dndy_pubpion_sys / dndy_pubpion, 2));
    double ratio_pion_staterr = sqrt(pow(dNdYSawanstat / dNdYSawan, 2) + pow(dndy_pubpion_stat / dndy_pubpion, 2));

    Double_t pointArray[] = {point};
    Double_t dNdYSawanArray[] = {dNdYSawan};
    Double_t ratio_kaonArray[] = {ratio_kaon};
    Double_t dndeta_error0Array[] = {dndeta_error[0]};
    Double_t dndeta_error1Array[] = {dndeta_error[1]};
    Double_t ratio_kaon_staterrArray[] = {ratio_kaon_staterr};
    Double_t ratio_kaon_syserrArray[] = {ratio_kaon_syserr};

    Double_t ratio_pionArray[] = {ratio_pion};
    Double_t ratio_pion_staterrArray[] = {ratio_pion_staterr};
    Double_t ratio_pion_syserrArray[] = {ratio_pion_syserr};

    auto grkaon_stat = new TGraphAsymmErrors(1, dNdYSawanArray, ratio_kaonArray, dndeta_error0Array, dndeta_error1Array, ratio_kaon_staterrArray, ratio_kaon_staterrArray);
    auto grkaon_sys = new TGraphAsymmErrors(1, dNdYSawanArray, ratio_kaonArray, dndeta_error0Array, dndeta_error1Array, ratio_kaon_syserrArray, ratio_kaon_syserrArray);

    auto grppion_stat = new TGraphAsymmErrors(1, dNdYSawanArray, ratio_pionArray, dndeta_error0Array, dndeta_error1Array, ratio_pion_staterrArray, ratio_pion_staterrArray);
    auto grppion_sys = new TGraphAsymmErrors(1, dNdYSawanArray, ratio_pionArray, dndeta_error0Array, dndeta_error1Array, ratio_pion_syserrArray, ratio_pion_syserrArray);

    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 800);
    SetCanvasStyle(c1, 0.15, 0.02, 0.02, 0.15);
    graphstyle(grkaon_stat);
    grkaon_stat->SetMarkerStyle(20);
    grkaon_stat->SetMarkerSize(1.5);
    grkaon_stat->SetMarkerColor(kBlack);
    grkaon_stat->SetLineColor(kBlack);
    grkaon_stat->SetLineWidth(2);
    grkaon_stat->GetYaxis()->SetTitle("K^{*0}/K");
    grkaon_stat->Draw("AP");
    grkaon_sys->SetFillStyle(0);
    grkaon_sys->SetLineColor(kBlack);
    grkaon_sys->SetLineWidth(2);
    grkaon_sys->Draw("e2SAME");
}