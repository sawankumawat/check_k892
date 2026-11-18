#include <iostream>
#include "../../src/style.h"
#include "../../src/fitfunc.h"

void Efficiency()
{
    gStyle->SetOptStat(0);
    TFile *file = new TFile("/home/sawan/Storage/check_k892/mc/PhiPhi/805651.root", "READ");
    if (file->IsZombie())
    {
        std::cerr << "Error: Could not open the file." << std::endl;
        return;
    }

    THnSparse *hPhiPhiGen = (THnSparse *)file->Get("phi-phi-m-c-injected/Generated/h2_InvMassPhiPhi");
    THnSparse *hPhiPhiRec = (THnSparse *)file->Get("phi-phi-m-c-injected/Reconstructed/AllPhiPairs/h2_PhiPhi_All");
    THnSparse *hPhiPhiRecOnlyInjected = (THnSparse *)file->Get("phi-phi-m-c-injected/Reconstructed/PhiFromT4s/h2_PhiPhi_T4s");

    if (!hPhiPhiGen || !hPhiPhiRec || !hPhiPhiRecOnlyInjected)
    {
        std::cerr << "Error: Could not retrieve the THnSparse objects." << std::endl;
        return;
    }

    TH1D *h1Gen = hPhiPhiGen->Projection(1);                         // Projecting onto the second axis (pt)
    TH1D *h1Rec = hPhiPhiRec->Projection(1);                         // Projecting onto the second axis (pt)
    TH1D *h1RecOnlyInjected = hPhiPhiRecOnlyInjected->Projection(1); // Projecting onto the second axis (pt)

    double pTbins[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0};
    int nPtBins = sizeof(pTbins) / sizeof(pTbins[0]) - 1;

    TH1D *hEfficiency = new TH1D("hEfficiency", "Reconstruction Efficiency vs pT;#it{p}_{T} (GeV/#it{c});Efficiency", nPtBins, pTbins);
    TH1D *hEfficiencyOnlyInjected = new TH1D("hEfficiencyOnlyInjected", "Reconstruction Efficiency (Only Injected) vs pT;#it{p}_{T} (GeV/#it{c});Efficiency", nPtBins, pTbins);

    for (int i = 0; i < nPtBins; ++i)
    {
        double genCount = h1Gen->Integral(h1Gen->GetXaxis()->FindBin(pTbins[i]), h1Gen->GetXaxis()->FindBin(pTbins[i + 1]));
        double recCount = h1Rec->Integral(h1Rec->GetXaxis()->FindBin(pTbins[i]), h1Rec->GetXaxis()->FindBin(pTbins[i + 1]));
        double recOnlyInjectedCount = h1RecOnlyInjected->Integral(h1RecOnlyInjected->GetXaxis()->FindBin(pTbins[i]), h1RecOnlyInjected->GetXaxis()->FindBin(pTbins[i + 1]));

        double efficiency = (genCount > 0) ? recCount / genCount : 0.0;
        double efficiencyOnlyInjected = (genCount > 0) ? recOnlyInjectedCount / genCount : 0.0;
        hEfficiency->SetBinContent(i + 1, efficiency);
        hEfficiencyOnlyInjected->SetBinContent(i + 1, efficiencyOnlyInjected);

        double efficiencyerr = sqrt(abs(((recCount + 1) / (genCount + 2)) * ((recCount + 2) / (genCount + 3) - (recCount + 1) / (genCount + 2))));
        hEfficiency->SetBinError(i + 1, efficiencyerr);

        double efficiencyOnlyInjectedErr = sqrt(abs(((recOnlyInjectedCount + 1) / (genCount + 2)) * ((recOnlyInjectedCount + 2) / (genCount + 3) - (recOnlyInjectedCount + 1) / (genCount + 2))));
        hEfficiencyOnlyInjected->SetBinError(i + 1, efficiencyOnlyInjectedErr);
    }

    TCanvas *cEfficiency = new TCanvas("cEfficiency", "Reconstruction Efficiency", 720, 720);
    SetCanvasStyle(cEfficiency, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hEfficiency);
    hEfficiency->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hEfficiency->GetYaxis()->SetTitle("Acc x Eff");
    hEfficiency->SetMinimum(0);
    hEfficiency->SetMaximum(0.14);
    hEfficiency->Draw("pe");
    SetHistoQA(hEfficiencyOnlyInjected);
    hEfficiencyOnlyInjected->SetLineColor(kRed);
    hEfficiencyOnlyInjected->SetMarkerColor(kRed);
    hEfficiencyOnlyInjected->SetMarkerStyle(21);
    hEfficiencyOnlyInjected->Draw("pe same");

    TLegend *legend = new TLegend(0.45, 0.4, 0.6, 0.6);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(hEfficiency, "All Reconstructed", "lep");
    legend->AddEntry(hEfficiencyOnlyInjected, "Only Injected Reconstructed", "lep");
    legend->Draw();
    cEfficiency->SaveAs("Reconstruction_Efficiency.png");
}