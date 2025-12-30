#include "src/style.h"
void BoltzmannYield()
{
    // Constants
    double T = 0.15;     // Temperature in GeV
    double m_pi = 0.139; // Mass of pion in GeV/c^2
    double m_p = 0.938;  // Mass of proton in GeV/c^2
    double m_k = 0.494;  // Mass of kaon in GeV/c^2
    double m_f2 = 1.525;
    double m_f0 = 1.710;

    // pT range from 0.01 to 20 GeV/c
    int nPoints = 500;
    double pTmin = 0.01;
    double pTmax = 20;
    double step = (pTmax - pTmin) / (nPoints - 1);

    TGraph *gr = new TGraph(nPoints);
    TGraph *gr_kaon = new TGraph(nPoints);
    TGraph *gr_f0 = new TGraph(nPoints);

    for (int i = 0; i < nPoints; ++i)
    {
        double pT = pTmin + i * step;

        // Yield formula: dN/dpT ~ pT * exp(-sqrt(pT^2 + m^2)/T)
        double yield_pi = pT * TMath::Exp(-TMath::Sqrt(pT * pT + m_pi * m_pi) / T);
        double yield_p = pT * TMath::Exp(-TMath::Sqrt(pT * pT + m_p * m_p) / T);
        double yield_k = pT * TMath::Exp(-TMath::Sqrt(pT * pT + m_k * m_k) / T);
        double yield_f0 = pT * TMath::Exp(-TMath::Sqrt(pT * pT + m_f0 * m_f0) / T);
        double yield_f2 = pT * TMath::Exp(-TMath::Sqrt(pT * pT + m_f2 * m_f2) / T);

        double ratio = yield_p / yield_pi;
        double ratio_kaon = yield_k / yield_pi;
        double ratio_f0_f2 = yield_f0 / yield_f2;
        gr->SetPoint(i, pT, ratio);
        gr_kaon->SetPoint(i, pT, ratio_kaon);
        gr_f0->SetPoint(i, pT, ratio_f0_f2);
    }

    // Plotting
    TCanvas *c1 = new TCanvas("c1", "Yield Ratio", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.03, 0.03, 0.15);
    SetGraphStyle(gr);
    gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gr->GetYaxis()->SetTitle("Yield Ratio");
    gr->SetTitle(0);
    gr->SetMaximum(1.6);
    gr->Draw("AL");
    SetGraphStyle(gr_kaon);
    gr_kaon->SetLineColor(kRed);
    gr_kaon->Draw("L same");
    SetGraphStyle(gr_f0);
    gr_f0->SetLineColor(kGreen + 2);
    gr_f0->Draw("L same");

    TLegend *leg = new TLegend(0.7, 0.75, 0.85, 0.9);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->AddEntry(gr, "Proton/Pion", "l");
    leg->AddEntry(gr_kaon, "Kaon/Pion", "l");
    leg->AddEntry(gr_f0, "f_{0}(1710)/f_{2}(1525)", "l");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.DrawLatex(0.19, 0.86, "dN/dp_{T} ~ p_{T} exp(-#frac{#sqrt{p_{T}^{2} + m^{2}}}{T})");
    latex.DrawLatex(0.19, 0.78, "T = 150 MeV");

    c1->SaveAs("Boltzmann_Yield_Ratios.png");
}
