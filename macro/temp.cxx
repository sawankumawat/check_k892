#include "src/style.h"
void temp()
{
    TFile *fInputFile = TFile::Open("PiKpGenRec.root");
    TH1D *hgenpip = (TH1D *)fInputFile->Get("piplus_generated");
    TH1D *hgenpim = (TH1D *)fInputFile->Get("piminus_generated");
    TH1D *hgenkp = (TH1D *)fInputFile->Get("Kaon_generated");
    TH1D *hgenkm = (TH1D *)fInputFile->Get("Kaon_minus_generated");
    TH1D *hrecpip = (TH1D *)fInputFile->Get("piplus_reconstructed");
    TH1D *hrecpim = (TH1D *)fInputFile->Get("piminus_reconstructed");
    TH1D *hreckp = (TH1D *)fInputFile->Get("Kaon_reconstructed");
    TH1D *hreckm = (TH1D *)fInputFile->Get("Kaon_minus_reconstructed");

    TH1D *hEffPiplus = (TH1D *)hrecpip->Clone("hEffPiplus");
    hEffPiplus->Divide(hgenpip);
    TH1D *hEffKaonPlus = (TH1D *)hreckp->Clone("hEffKaonPlus");
    hEffKaonPlus->Divide(hgenkp);
    TH1D *hEffKaonMinus = (TH1D *)hreckm->Clone("hEffKaonMinus");
    hEffKaonMinus->Divide(hgenkm);

    TCanvas *cEffPiplus = new TCanvas("cEffPiplus", "Efficiency of all", 720, 720);
    SetCanvasStyle(cEffPiplus, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hEffPiplus);
    hEffPiplus->GetYaxis()->SetTitle("Efficiency");
    hEffPiplus->Draw();
    hEffKaonPlus->SetLineColor(kRed);
    hEffKaonPlus->SetMarkerColor(kRed);
    hEffKaonPlus->Draw("same");
    hEffKaonMinus->SetLineColor(kBlue);
    hEffKaonMinus->SetMarkerColor(kBlue);
    hEffKaonMinus->Draw("same");

    TLegend *leg = new TLegend(0.3, 0.25, 0.6, 0.5);
    leg->AddEntry(hEffPiplus, "Pion+", "pe");
    leg->AddEntry(hEffKaonPlus, "Kaon+", "pe");
    leg->AddEntry(hEffKaonMinus, "Kaon-", "pe");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->Draw();
}