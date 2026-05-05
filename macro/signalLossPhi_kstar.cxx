#include <iostream>
#include "src/style.h"

using namespace std;

void signalLossPhi_kstar()
{
    gStyle->SetOptStat(0);

    string common_path = "/home/sawan/check_k892/mc/LHC24f3c/";
    //===========2024 data==============
    // string MCpath = "658376.root"; // (For Phi signal loss, without TOF_overrideFT0) (wrong value, since AllEvents has RCT, which other two does not have it)
    // string MCpath = "659253.root"; // (For Phi signal loss) (do not use)
    // string MCpath = "661190.root"; // (SL with RCT (id52865), wo RCT (id52866))
    // string MCpath = "668755.root"; // (SL, with new addition |Vz| and inelgt0 cuts)
    string MCpath = "669655.root"; // (AllEvents, AllEventsVz, AllEventsVz2, sel8, VzSel8) (Without extra Vz cut in Gen MC)

    // TFile *fSigLoss = new TFile("/home/sawan/check_k892/mc/PhiSigLossINEL.root", "READ");
    // TFile *fSigLoss = new TFile("/home/sawan/Downloads/SigLossLHC24l1.root", "READ");
    TFile *fSigLoss = new TFile((common_path + MCpath).c_str(), "READ");
    if (fSigLoss->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }

    // double pT_bins[] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0};
    double pT_bins[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0};

    // string path1 = "kstarqa_AllEvents";
    string path1 = "kstarqa_AllEventsVz2";
    // string path2 = "kstarqa_sel8";
    string path2 = "kstarqa_VzSel8";

    THnSparseF *hSpraseGenAll = (THnSparseF *)fSigLoss->Get(Form("%s/hInvMass/hk892GenpT", path1.c_str()));
    THnSparseF *hSpraseGenVzSel = (THnSparseF *)fSigLoss->Get(Form("%s/hInvMass/hk892GenpT", path2.c_str()));
    if (hSpraseGenAll == nullptr || hSpraseGenVzSel == nullptr)
    {
        cout << "Error reading histograms" << endl;
        return;
    }
    TH1D *hGenAll = hSpraseGenAll->Projection(0, "E");
    hGenAll->SetName("hGenAll_proj");
    TH1D *hGenVzSel = hSpraseGenVzSel->Projection(0, "E");
    hGenVzSel->SetName("hGenVzSel_proj");

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
    hSignalLoss->SetMaximum(1.385);
    hSignalLoss->SetMinimum(1.29);
    hSignalLoss->Draw("pe");
    // cSignalLoss->SaveAs("SignalLoss_phi.png");

    // TFile *fOutput = new TFile("SignalLossPhiINEL.root", "RECREATE");
    // fOutput->cd();
    // hSignalLoss->Write("hSignalLoss");
    // fOutput->Close();

    // // Cleanup
    // delete hGenAll;
    // delete hGenVzSel;
    // delete hSignalLoss;
    // fSigLoss->Close();
    // delete fSigLoss;
}