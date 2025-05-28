#include <iostream>
#include "../src/style.h"
using namespace std;

void triggered_data()
{
    TFile *file = new TFile("../../data/doublePhi/trigger_data/AnalysisResults_fullrun.root");
    if (file->IsZombie())
    {
        std::cerr << "Error: Could not open file triggered_data/AnalysisResults_fullrun.root\n";
        return;
    }

    TH1F *hevent_processed = (TH1F *)file->Get("lf-doublephi-filter/hProcessedEvents");
    TH2F *hphiMass = (TH2F *)file->Get("lf-doublephi-filter/hInvMassPhi");
    TH2F *hdoublePhiMass = (TH2F *)file->Get("lf-doublephi-filter/hInvMassDoublePhi");
    TH2F *hNsigmaTPCKaon = (TH2F *)file->Get("lf-doublephi-filter/hNsigmaPtkaonTPC");
    TH2F *hNsigmaTOFKaon = (TH2F *)file->Get("lf-doublephi-filter/hNsigmaPtkaonTOF");
    if (hevent_processed == nullptr || hphiMass == nullptr || hdoublePhiMass == nullptr || hNsigmaTPCKaon == nullptr || hNsigmaTOFKaon == nullptr)
    {
        std::cerr << "Error: Could not find one or more histograms in the file\n";
        return;
    }
    gStyle->SetOptStat(1110);
    gStyle->SetGridColor(15);
    TCanvas *c1 = new TCanvas("c1", "Processed Events", 1080, 720);
    c1->SetGrid(1, 1);
    SetCanvasStyle(c1, 0.11, 0.06, 0.05, 0.10);
    gPad->SetLogy();
    SetHistoQA(hevent_processed);
    hevent_processed->GetYaxis()->SetTitleOffset(1.1);
    hevent_processed->Draw();
    c1->SaveAs("trigger_data/processed_events.png");

    TCanvas *c2 = new TCanvas("c2", "Phi Mass", 720, 720);
    SetCanvasStyle(c2, 0.13, 0.06, 0.05, 0.13);
    int lowbinpT = hphiMass->GetYaxis()->FindBin(0.0 + 0.001);
    int highbinpT = hphiMass->GetYaxis()->FindBin(10.0 - 0.001);
    TH1D *hphiMassProj = hphiMass->ProjectionX("hphiMassProj", lowbinpT, highbinpT);
    SetHistoQA(hphiMassProj);
    hphiMassProj->GetYaxis()->SetTitleOffset(1.3);
    hphiMassProj->GetYaxis()->SetMaxDigits(3);
    hphiMassProj->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
    hphiMassProj->GetXaxis()->SetNdivisions(505);
    hphiMassProj->GetYaxis()->SetTitle("Counts");
    hphiMassProj->Draw("pe");
    c2->SaveAs("trigger_data/phi_mass.png");

    TCanvas *c3 = new TCanvas("c3", "Double Phi Mass", 720, 720);
    SetCanvasStyle(c3, 0.13, 0.06, 0.05, 0.13);
    int lowbinpT2 = hdoublePhiMass->GetYaxis()->FindBin(2.0 + 0.001);
    int highbinpT2 = hdoublePhiMass->GetYaxis()->FindBin(5.0 - 0.001);
    TH1D *hdoublePhiMassProj = hdoublePhiMass->ProjectionX("hdoublePhiMassProj", lowbinpT2, highbinpT2);
    SetHistoQA(hdoublePhiMassProj);
    hdoublePhiMassProj->GetYaxis()->SetTitleOffset(1.3);
    hdoublePhiMassProj->GetYaxis()->SetMaxDigits(3);
    hdoublePhiMassProj->GetXaxis()->SetTitle("#it{M}_{#Phi#Phi} (GeV/#it{c}^{2})");
    hdoublePhiMassProj->GetXaxis()->SetNdivisions(505);
    hdoublePhiMassProj->Rebin(4);
    hdoublePhiMassProj->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", (hdoublePhiMassProj->GetXaxis()->GetBinWidth(1) * 1000)));
    hdoublePhiMassProj->GetXaxis()->SetRangeUser(2.6, 2.9);
    hdoublePhiMassProj->GetXaxis()->SetNdivisions(508);
    hdoublePhiMassProj->Draw("pe");
    c3->SaveAs("trigger_data/double_phi_mass_pt2to5.png");

    TCanvas *c4 = new TCanvas("c4", "Nsigma TPC Kaon", 720, 720);
    SetCanvasStyle(c4, 0.13, 0.14, 0.05, 0.13);
    // SetHistoQA(hNsigmaTPCKaon);
    hNsigmaTPCKaon->SetTitle(0);
    hNsigmaTPCKaon->GetXaxis()->SetTitle("#it{n}_{#sigma} (TPC)");
    hNsigmaTPCKaon->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hNsigmaTPCKaon->GetYaxis()->SetTitleOffset(1.3);
    hNsigmaTPCKaon->GetYaxis()->SetMaxDigits(3);
    hNsigmaTPCKaon->GetXaxis()->SetNdivisions(505);
    hNsigmaTPCKaon->GetXaxis()->SetRangeUser(-5.0, 5.0);
    hNsigmaTPCKaon->GetZaxis()->SetMaxDigits(3);
    hNsigmaTPCKaon->SetStats(0);
    hNsigmaTPCKaon->Draw("colz");
    c4->SaveAs("trigger_data/nsigma_tpc_kaon.png");

    TCanvas *c5 = new TCanvas("c5", "Nsigma TOF Kaon", 720, 720);
    SetCanvasStyle(c5, 0.13, 0.14, 0.05, 0.13);
    // SetHistoQA(hNsigmaTOFKaon);
    hNsigmaTOFKaon->SetTitle(0);
    hNsigmaTOFKaon->GetXaxis()->SetTitle("#it{n}_{#sigma} (TOF)");
    hNsigmaTOFKaon->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hNsigmaTOFKaon->GetYaxis()->SetTitleOffset(1.3);
    hNsigmaTOFKaon->GetYaxis()->SetMaxDigits(3);
    hNsigmaTOFKaon->GetXaxis()->SetNdivisions(505);
    hNsigmaTOFKaon->GetXaxis()->SetRangeUser(-5.0, 5.0);
    hNsigmaTOFKaon->GetZaxis()->SetMaxDigits(3);
    hNsigmaTOFKaon->SetStats(0);
    hNsigmaTOFKaon->Draw("colz");
    c5->SaveAs("trigger_data/nsigma_tof_kaon.png");

    TCanvas *c6 = new TCanvas("c6", "Nsigma TPC", 720, 720);
    SetCanvasStyle(c6, 0.13, 0.14, 0.05, 0.13);
    int lowbinpt_nsigma = hNsigmaTPCKaon->GetYaxis()->FindBin(0.5 + 0.001);
    int highbinpt_nsigma = hNsigmaTPCKaon->GetYaxis()->FindBin(5.0 - 0.001);
    TH1F *hNsigmaTPC = (TH1F *)hNsigmaTPCKaon->ProjectionX("hNsigmaTPC", lowbinpt_nsigma, highbinpt_nsigma);
    hNsigmaTPC->Draw();
    //Draw a vertical line at x=0
    TLine *line = new TLine(0, 0, 0, hNsigmaTPC->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->Draw("same");
    c6->SaveAs("trigger_data/nsigma_tpc1D.png");

    TCanvas *c7 = new TCanvas("c7", "Nsigma TOF", 720, 720);
    SetCanvasStyle(c7, 0.13, 0.14, 0.05, 0.13);
    TH1F *hNsigmaTOF = (TH1F *)hNsigmaTOFKaon->ProjectionX("hNsigmaTOF", lowbinpt_nsigma, highbinpt_nsigma);
    hNsigmaTOF->Draw();
    TLine *line2 = new TLine(0, 0, 0, hNsigmaTOF->GetMaximum());
    line2->SetLineColor(kRed);
    line2->SetLineStyle(2);
    line2->Draw("same");
    c7->SaveAs("trigger_data/nsigma_tof1D.png");

}