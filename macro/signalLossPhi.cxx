#include <iostream>
#include "src/style.h"

using namespace std;

void signalLossPhi()
{
    gStyle->SetOptStat(0);
    // TFile *fSigLoss = new TFile("/home/sawan/check_k892/mc/PhiSigLossINEL.root", "READ");
    TFile *fSigLoss = new TFile("/home/sawan/Downloads/SigLossLHC24l1.root", "READ");
    if (fSigLoss->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }

    double pT_bins[] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0};

    string path1 = "kstarqa_AllEvents";
    string path2 = "kstarqa_AllEventsVz";
    // string path2 = "kstarqa_VzSel8";

    THnSparseF *hSpraseGenAll = (THnSparseF *)fSigLoss->Get(Form("%s/hInvMass/hk892GenpT2", path1.c_str()));
    THnSparseF *hSpraseGenVzSel = (THnSparseF *)fSigLoss->Get(Form("%s/hInvMass/hk892GenpT2", path2.c_str()));
    if (hSpraseGenAll == nullptr || hSpraseGenVzSel == nullptr)
    {
        cout << "Error reading histograms" << endl;
        return;
    }
    TH1D *hGenAll = hSpraseGenAll->Projection(0, "E");
    TH1D *hGenVzSel = hSpraseGenVzSel->Projection(0, "E");

    TH1D *hSignalLoss = new TH1D("hSignalLoss", "Signal Loss vs pT", sizeof(pT_bins) / sizeof(pT_bins[0]) - 1, pT_bins);
    for (int i = 1; i <= hSignalLoss->GetNbinsX(); i++)
    {
        double genAll = hGenAll->GetBinContent(i);
        double genVzSel = hGenVzSel->GetBinContent(i);
        double loss = genAll / genVzSel; // Signal loss
        hSignalLoss->SetBinContent(i, loss);
        double efficiencyerr = sqrt(abs(((genVzSel + 1) / (genAll + 2)) * ((genVzSel + 2) / (genAll + 3) - (genVzSel + 1) / (genAll + 2))));
        hSignalLoss->SetBinError(i, efficiencyerr);
    }

    TCanvas *cSignalLoss = new TCanvas("cSignalLoss", "Signal Loss vs pT", 720, 720);
    SetCanvasStyle(cSignalLoss, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hSignalLoss);
    hSignalLoss->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSignalLoss->GetYaxis()->SetTitle("Signal Loss");
    hSignalLoss->GetYaxis()->SetTitleOffset(1.6);
    hSignalLoss->SetMarkerSize(1.5);
    hSignalLoss->SetMaximum(1.38);
    hSignalLoss->SetMinimum(0.65);
    hSignalLoss->Draw("pe");
    // cSignalLoss->SaveAs("SignalLoss_phi.png");

    // TFile *fOutput = new TFile("SignalLossPhiINEL.root", "RECREATE");
    // fOutput->cd();
    // hSignalLoss->Write("hSignalLoss");
    // fOutput->Close();
}