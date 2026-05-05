#include <iostream>
#include <vector>
#include "../macro/src/style.h"

void plot_dNbydy_Meanpt()
{
    gStyle->SetOptStat(0);

    const int n = 10;

    // Run 3 data
    double mult_run3[n] = {28.0094, 22.7302, 18.8714, 16.1524, 13.9717,
                           11.4394, 8.95528, 7.08418, 4.85499, 3.32155};

    double pt_run3[n] = {1.23787, 1.19252, 1.14361, 1.09794, 1.05306,
                         0.994975, 0.925666, 0.868026, 0.792003, 0.689025};

    double dNdy_run3[n] = {1.05595, 0.856883, 0.709512, 0.606239, 0.52315, 0.426186, 0.330544, 0.258226, 0.171096, 0.111019};

    // Run 2 data
    double mult_run2[n] = {28.3475, 22.7926, 18.8654, 15.9551, 13.7599,
                           11.1458, 8.4095, 6.58427, 4.53274, 3.1222};

    double pt_run2[n] = {1.2325, 1.18804, 1.13915, 1.0885, 1.04306,
                         0.980486, 0.901382, 0.842945, 0.774465, 0.663632};

    double dNdy_run2[n] = {1.06911, 0.85828, 0.709977, 0.59908, 0.515298, 0.415111, 0.309613, 0.238955, 0.15853, 0.103035};

    // Create graphs
    TGraph *gRun3_pt = new TGraph(n, mult_run3, pt_run3);
    TGraph *gRun2_pt = new TGraph(n, mult_run2, pt_run2);
    TGraph *gRun3_dNdy = new TGraph(n, mult_run3, dNdy_run3);
    TGraph *gRun2_dNdy = new TGraph(n, mult_run2, dNdy_run2);

    // Style
    gRun3_pt->SetMarkerStyle(20);
    gRun3_pt->SetMarkerColor(kRed);
    gRun3_pt->SetLineColor(kRed);
    gRun3_pt->SetMarkerSize(1.3);

    gRun2_pt->SetMarkerStyle(21);
    gRun2_pt->SetMarkerColor(kBlue);
    gRun2_pt->SetLineColor(kBlue);
    gRun2_pt->SetMarkerSize(1.3);

    gRun3_dNdy->SetMarkerStyle(20);
    gRun3_dNdy->SetMarkerColor(kRed);
    gRun3_dNdy->SetLineColor(kRed);
    gRun3_dNdy->SetMarkerSize(1.3);

    gRun2_dNdy->SetMarkerStyle(21);
    gRun2_dNdy->SetMarkerColor(kBlue);
    gRun2_dNdy->SetLineColor(kBlue);
    gRun2_dNdy->SetMarkerSize(1.3);

    // Canvas
    TCanvas *c = new TCanvas("c", "pT vs multiplicity", 720, 720);
    SetCanvasStyle(c, 0.14, 0.02, 0.05, 0.14);

    gRun3_pt->SetTitle(";dN_{ch}/d#eta; <p_{T}> (GeV/c)");
    gRun3_pt->Draw("AP"); // Axis + Points
    gRun2_pt->Draw("P SAME");

    // Legend
    TLegend *leg = new TLegend(0.6, 0.2, 0.85, 0.35);
    leg->AddEntry(gRun3_pt, "Run 3", "p");
    leg->AddEntry(gRun2_pt, "Run 2", "p");
    leg->Draw();
    c->SaveAs("ptSpectra/pt_vs_mult.png");

    TCanvas *c2 = new TCanvas("c2", "dN/dy vs multiplicity", 720, 720);
    SetCanvasStyle(c2, 0.14, 0.02, 0.05, 0.14);
    gRun3_dNdy->SetTitle(";dN_{ch}/d#eta; dN/dy");
    gRun3_dNdy->Draw("AP");
    gRun2_dNdy->Draw("P SAME");
    leg->Draw();
    c2->SaveAs("ptSpectra/dNdy_vs_mult.png");
}