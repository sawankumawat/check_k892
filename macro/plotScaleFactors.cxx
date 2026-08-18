#include <iostream>
#include "src/style.h"

void plotScaleFactors()
{
    //===========================================================
    // 1. First scale factor: getKaonMomentumScale()
    //===========================================================
    const int nBins1 = 15;

    double edges1[nBins1 + 1] = {
        0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0,
        4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};

    double scale1[nBins1] = {
        -0.008062,
        0.004021,
        0.008302,
        0.009222,
        0.010662,
        0.011902,
        0.012382,
        0.012182,
        0.012422,
        0.012703,
        0.012462,
        0.014223,
        0.014503,
        0.016463,
        0.019064};

    //===========================================================
    // 2. Second scale factor: kPtEdges / kEpsilon
    //===========================================================
    const int nBins2 = 17;

    double edges2[nBins2 + 1] = {
        0.8, 1.0, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0,
        6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0};

    double scale2[nBins2] = {
        0.0082316288136,
        0.0084477004712,
        0.00919943176938,
        0.0108425754669,
        0.0117196805635,
        0.0123951594009,
        0.0123680092232,
        0.0124683175364,
        0.0126507831441,
        0.0124286969418,
        0.0123413824165,
        0.0134753484815,
        0.0150030728106,
        0.014280642284,
        0.0149542169164,
        0.0174764570162,
        0.0196359204232};

    //===========================================================
    // Calculate bin centers
    //===========================================================
    double x1[nBins1];
    double x2[nBins2];

    for (int i = 0; i < nBins1; i++)
        x1[i] = 0.5 * (edges1[i] + edges1[i + 1]);

    for (int i = 0; i < nBins2; i++)
        x2[i] = 0.5 * (edges2[i] + edges2[i + 1]);

    //===========================================================
    // Create graphs
    //===========================================================
    TGraph *g1 = new TGraph(nBins1, x1, scale1);
    TGraph *g2 = new TGraph(nBins2, x2, scale2);

    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(1.1);
    g1->SetLineWidth(2);

    g2->SetMarkerStyle(21);
    g2->SetMarkerSize(1.1);
    g2->SetLineWidth(2);

    //===========================================================
    // Canvas
    //===========================================================
    TCanvas *c = new TCanvas("c", "Kaon Momentum Scale Factors",
                             900, 700);
    SetCanvasStyle(c, 0.15, 0.03, 0.03, 0.15);

    g1->SetTitle("");
    g1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    g1->GetYaxis()->SetTitle("Momentum scale factor (#epsilon)");
    
    g1->GetXaxis()->SetLimits(0.4, 20.0);
    g1->GetYaxis()->SetRangeUser(-0.01, 0.022);
    
    SetGraphStyle(g1);
    SetGraphStyle(g2);
    g1->GetYaxis()->SetTitleOffset(1.5);
    
    g1->Draw("APL");
    g2->SetLineColor(kBlue);
    g2->SetMarkerColor(kBlue);
    g2->Draw("PL SAME");

    //===========================================================
    // Legend
    //===========================================================
    TLegend *leg = new TLegend(0.5, 0.3, 0.8, 0.5);
    leg->SetFillStyle(0);
    leg->AddEntry(g1, "Toy model (Sawan)", "lp");
    leg->AddEntry(g2, "Sourav Bhaiya code", "lp");

    leg->SetBorderSize(0);
    leg->Draw();

    c->SetGridx();
    c->SetGridy();

    c->SaveAs("/home/sawan/Documents/kaonMomentumScaleComparison.png");
}