#include <iostream>
#include "src/style.h"

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}

void compare_efficiency()
{
    // TString path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/664559/kstarqa/hInvMass"; // with TOF shift
    // TString path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/663738/kstarqa/hInvMass"; // without TOF shift
    // TString path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/664559/kstarqa/hInvMass"; // with OverrideFT0
    // TString path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/664559/kstarqa/hInvMass"; // without overrideFT0
    // TString path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa/hInvMass";
    // TString path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa_MIDptDep2/hInvMass";
    // TString path3 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa_MIDptDep2_small/hInvMass";
    // TString path4 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa_MIDptDep2_verySmall/hInvMass";
    // TString path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa_TOF3_withoutSquareCut/hInvMass";
    // TString path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/672297/kstarqa_MIDptDep2_TOF3/hInvMass";
    // TString path3 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/672297/kstarqa_MIDptDep2_small_TOF3/hInvMass";
    // TString path4 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/672297/kstarqa_MIDptDep2_0p3_TOF3/hInvMass";

    // TString path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa_TOF3_withoutSquareCut/hInvMass";
    // TString path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/675391/kstarqa_MIDNew_TOF3/hInvMass";
    // TString path3 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/675391/kstarqa_SquarePID_TOF3/hInvMass";
    // TString path4 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/675391/kstarqa_SquarePID_TOF3/hInvMass"; // dummy

    TString path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa/hInvMass";
    TString path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa_TOF3_withoutSquareCut/hInvMass";
    TString path3 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa_TOF3_withoutSquareCut/hInvMass"; // dummy
    TString path4 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/668605/kstarqa_TOF3_withoutSquareCut/hInvMass"; // dummy

    TFile *fileShift = new TFile(path1 + "/corrected_spectra.root", "read");
    TFile *fileNoShift = new TFile(path2 + "/corrected_spectra.root", "read");
    TFile *fileSmall = new TFile(path3 + "/corrected_spectra.root", "read");
    TFile *fileVerySmall = new TFile(path4 + "/corrected_spectra.root", "read");
    if (fileShift->IsZombie() || fileNoShift->IsZombie() || fileSmall->IsZombie() || fileVerySmall->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }
    TH1F *hEffDefault = (TH1F *)fileShift->Get("mult_0-100/heff");
    TH1F *hEffVar1 = (TH1F *)fileNoShift->Get("mult_0-100/heff");
    TH1F *hEffVar2 = (TH1F *)fileSmall->Get("mult_0-100/heff");
    TH1F *hEffVar3 = (TH1F *)fileVerySmall->Get("mult_0-100/heff");
    if (hEffDefault == nullptr || hEffVar1 == nullptr || hEffVar2 == nullptr || hEffVar3 == nullptr)
    {
        cout << "Error: histograms not found" << endl;
        return;
    }

    TCanvas *cEff = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cEff, 0.16, 0.06, 0.01, 0.14);
    double pad1Size, pad2Size;
    canvas_style(cEff, pad1Size, pad2Size);
    cEff->cd(1);
    SetHistoQA(hEffDefault);
    hEffDefault->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hEffDefault->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hEffDefault->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hEffDefault->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hEffDefault->GetYaxis()->SetTitleOffset(1.7 * pad1Size);
    hEffDefault->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hEffDefault->GetYaxis()->SetTitle("Efficiency");
    hEffDefault->SetMaximum(hEffDefault->GetMaximum() * 1.4);
    hEffDefault->Draw("ep");
    hEffVar1->SetMarkerStyle(21);
    hEffVar1->SetMarkerColor(kBlue);
    hEffVar1->SetLineColor(kBlue);
    hEffVar1->Draw("ep same");
    hEffVar2->SetMarkerStyle(22);
    hEffVar2->SetMarkerColor(kRed);
    hEffVar2->SetLineColor(kRed);
    // hEffVar2->Draw("ep same");
    hEffVar3->SetMarkerStyle(23);
    hEffVar3->SetMarkerColor(kGreen + 2);
    hEffVar3->SetLineColor(kGreen + 2);
    // hEffVar3->Draw("ep same");
    TLegend *legEff = new TLegend(0.20, 0.8, 0.5, 0.98);
    legEff->SetTextSize(0.04);
    legEff->SetFillStyle(0);
    legEff->SetBorderSize(0);
    // legEff->AddEntry(hEffDefault, "With overrideFT0", "p");
    // legEff->AddEntry(hEffVar1, "Without overrideFT0", "p");
    // legEff->AddEntry(hEffDefault, "With TOF shift", "p");
    // legEff->AddEntry(hEffVar1, "Without TOF shift", "p");

    // legEff->AddEntry(hEffVar1, "No MID", "p");
    // legEff->AddEntry(hEffDefault, "MID cut (1#sigma)", "p");
    // legEff->AddEntry(hEffVar2, "MID cut (0.5#sigma)", "p");
    // legEff->AddEntry(hEffVar3, "MID cut (0.1#sigma)", "p");

    // legEff->AddEntry(hEffVar1, "No MID", "p");
    // legEff->AddEntry(hEffDefault, "MID cut (1#sigma)", "p");
    // legEff->AddEntry(hEffVar2, "MID cut (0.5#sigma)", "p");
    // legEff->AddEntry(hEffVar3, "MID cut (0.3#sigma)", "p");

    // legEff->AddEntry(hEffDefault, "No MID", "p");
    // legEff->AddEntry(hEffVar1, "New MID", "p");
    // legEff->AddEntry(hEffVar2, "Square PID TOF3", "p");

    // legEff->AddEntry(hEffDefault, "No TOF shift", "p");
    // legEff->AddEntry(hEffVar1, "TOF shift", "p");

    legEff->AddEntry(hEffDefault, "TOF 2#sigma", "p");
    legEff->AddEntry(hEffVar1, "TOF 3#sigma", "p");

    legEff->Draw();

    cEff->cd(2);
    gPad->SetGridy(1);
    TH1F *hRatio1 = (TH1F *)hEffVar1->Clone("hRatio1");
    SetHistoQA(hRatio1);
    hRatio1->Divide(hEffDefault);
    hRatio1->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hRatio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hRatio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hRatio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio1->GetYaxis()->SetTitle("Ratio");
    // hRatio1->GetYaxis()->SetTitle("Shift / No Shift");
    hRatio1->GetYaxis()->SetTitleOffset(1.5 * pad2Size);
    hRatio1->SetMarkerStyle(20);
    hRatio1->SetMarkerSize(1.0);
    hRatio1->GetYaxis()->SetNdivisions(505);
    // hRatio1->GetYaxis()->SetRangeUser(0.7, 1.3);
    // hRatio1->GetYaxis()->SetRangeUser(0.7, 1.58);
    hRatio1->GetYaxis()->SetRangeUser(0.52, 1.38);
    hRatio1->Draw("ep");
    TH1F *hRatio2 = (TH1F *)hEffVar2->Clone("hRatio2");
    SetHistoQA(hRatio2);
    hRatio2->Divide(hEffDefault);
    hRatio2->SetMarkerStyle(21);
    hRatio2->SetMarkerColor(kRed);
    hRatio2->SetLineColor(kRed);
    // hRatio2->Draw("ep same");
    TH1F *hRatio3 = (TH1F *)hEffVar3->Clone("hRatio3");
    SetHistoQA(hRatio3);
    hRatio3->Divide(hEffDefault);
    hRatio3->SetMarkerStyle(22);
    hRatio3->SetMarkerColor(kGreen + 2);
    hRatio3->SetLineColor(kGreen + 2);
    // hRatio3->Draw("ep same");
    TLine *line = new TLine(0, 1, 10, 1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw();
    cEff->SaveAs(path1 + "/efficiency_comparison.png");
}