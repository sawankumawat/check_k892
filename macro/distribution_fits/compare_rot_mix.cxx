#include <iostream>
#include <vector>
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

void compare_rot_mix()
{
    gStyle->SetOptStat(0);

    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937";
    TFile *frot = new TFile((path + "/hglue_ROTATED_norm_2.50_2.60_pt_3.00_5.00.root").c_str(), "READ"); //
    TFile *fmix = new TFile((path + "/hglue_MIX_norm_2.50_2.60_pt_3.00_5.00.root").c_str(), "READ");     //
    if (frot->IsZombie() || fmix->IsZombie())
    {
        cout << "Files not found" << endl;
        return;
    }

    TH1F *hsignal = (TH1F *)frot->Get("ksks_invmass_pt_3.0_5.0");
    TH1F *hrot = (TH1F *)frot->Get("ksks_bkg_pt_3.0_5.0");
    TH1F *hmix = (TH1F *)fmix->Get("ksks_bkg_pt_3.0_5.0");
    TH1F *hrot_signal = (TH1F *)frot->Get("ksks_subtracted_invmass_pt_3.0_5.0");
    TH1F *hmix_signal = (TH1F *)fmix->Get("ksks_subtracted_invmass_pt_3.0_5.0");

    if (hsignal == nullptr || hrot == nullptr || hmix == nullptr)
    {
        cout << "Histogram not found" << endl;
        return;
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.15, 0.03, 0.05, 0.14);
    hsignal->GetXaxis()->SetRangeUser(1.00, 2.60);
    SetHistoQA(hsignal);
    hsignal->Draw("pe");

    SetHistoQA(hrot);
    hrot->GetXaxis()->SetRangeUser(1.00, 2.60);
    hrot->SetMarkerSize(1.0);
    hrot->SetMarkerColor(kBlue);
    hrot->SetLineColor(kBlue);
    hrot->SetMarkerStyle(23);
    hrot->Draw("same");

    SetHistoQA(hmix);
    hmix->GetXaxis()->SetRangeUser(1.00, 2.60);
    hmix->SetMarkerSize(0.8);
    hmix->SetMarkerColor(kRed);
    hmix->SetLineColor(kRed);
    hmix->SetMarkerStyle(21);
    hmix->Draw("same");

    TLegend *leg = new TLegend(0.55, 0.62, 0.92, 0.92);
    SetLegendStyle(leg);
    leg->SetHeader("3.0 < #it{p}_{T} < 5.0 GeV/c");
    leg->AddEntry(hsignal, "Raw K^{0}_{s}K^{0}_{s} invariant mass", "lpe");
    leg->AddEntry(hrot, "Rotational background", "lpe");
    leg->AddEntry(hmix, "Mixed-event background", "lpe");
    leg->Draw("same");

    c->SaveAs("/home/sawan/Music/compare_rot_mix.png");

    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.03, 0.05, 0.14);
    hrot_signal->GetXaxis()->SetRangeUser(1.00, 2.60);
    SetHistoQA(hrot_signal);
    hrot_signal->SetLineColor(kBlue);
    hrot_signal->SetMarkerColor(kBlue);
    hrot_signal->SetMarkerStyle(23);
    hrot_signal->SetMarkerSize(1.0);
    // hrot_signal->SetMinimum(-40000);
    hrot_signal->Draw("pe");

    SetHistoQA(hmix_signal);
    hmix_signal->GetXaxis()->SetRangeUser(1.00, 2.60);
    hmix_signal->SetMarkerSize(0.8);
    hmix_signal->SetMarkerColor(kRed);
    hmix_signal->SetLineColor(kRed);
    hmix_signal->SetMarkerStyle(21);
    hmix_signal->Draw("same");

    TLine *line = new TLine(1.00, 0, 2.60, 0);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw();

    TLegend *leg2 = new TLegend(0.30, 0.70, 0.92, 0.92);
    SetLegendStyle(leg2);
    leg2->SetHeader("3.0 < #it{p}_{T} < 5.0 GeV/c");
    leg2->AddEntry(hrot_signal, "Rotational background subtraction", "lpe");
    leg2->AddEntry(hmix_signal, "Mixed-event background subtraction", "lpe");
    leg2->Draw("same");

    c2->SaveAs("/home/sawan/Music/compare_rot_mix_subtracted.png");
}
