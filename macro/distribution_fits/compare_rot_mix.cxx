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
    TFile *frot = new TFile((path + "/hglue_ROTATED.root").c_str(), "READ"); //
    TFile *fmix = new TFile((path + "/hglue_MIX.root").c_str(), "READ");     //
    if (frot->IsZombie() || fmix->IsZombie())
    {
        cout << "Files not found" << endl;
        return;
    }

    TH1F *hsignal = (TH1F *)frot->Get("ksks_invmass_pt_0.0_30.0");
    TH1F *hrot = (TH1F *)frot->Get("ksks_bkg_pt_0.0_30.0");
    TH1F *hmix = (TH1F *)fmix->Get("ksks_bkg_pt_0.0_30.0");
    TH1F *hrot_signal = (TH1F *)frot->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    TH1F *hmix_signal = (TH1F *)fmix->Get("ksks_subtracted_invmass_pt_0.0_30.0");

    if (hsignal == nullptr || hrot == nullptr || hmix == nullptr)
    {
        cout << "Histogram not found" << endl;
        return;
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.15, 0.03, 0.05, 0.14);
    hsignal->GetXaxis()->SetRangeUser(1.00, 2.99);
    SetHistoQA(hsignal);
    hsignal->SetMaximum(1.1 * hsignal->GetMaximum());
    hsignal->Draw("pe");

    SetHistoQA(hrot);
    hrot->GetXaxis()->SetRangeUser(1.00, 2.99);
    hrot->SetMarkerSize(1.0);
    hrot->SetMarkerColor(kBlue);
    hrot->SetLineColor(kBlue);
    hrot->SetMarkerStyle(23);
    hrot->Draw("same");

    SetHistoQA(hmix);
    hmix->GetXaxis()->SetRangeUser(1.00, 2.99);
    hmix->SetMarkerSize(0.8);
    hmix->SetMarkerColor(kRed);
    hmix->SetLineColor(kRed);
    hmix->SetMarkerStyle(21);
    hmix->Draw("same");

    TH1F *hbkg_nopeak = (TH1F *)hmix->Clone();
    hbkg_nopeak->SetLineColor(40);
    hbkg_nopeak->SetMarkerColor(40);
    hbkg_nopeak->SetFillColor(40);
    // hbkg_nopeak->SetFillStyle(3001);
    for (int i = 0; i < hbkg_nopeak->GetNbinsX(); i++)
    {
        if (hbkg_nopeak->GetBinCenter(i + 1) < 2.5 || hbkg_nopeak->GetBinCenter(i + 1) > 2.6)
        {
            hbkg_nopeak->SetBinContent(i + 1, -999);
        }
    }
    hbkg_nopeak->Draw("BAR same");

    TLegend *legAlice = new TLegend(0.48, 0.79, 0.92, 0.92);
    SetLegendStyle(legAlice);
    legAlice->SetTextSize(0.035);
    legAlice->SetTextFont(42);
    legAlice->SetFillStyle(0);
    legAlice->SetBorderSize(0);
    legAlice->AddEntry((TObject *)0, "ALICE", "");
    legAlice->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    // legAlice->AddEntry((TObject *)0, "FT0M, 0-100%", "");
    // legAlice->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
    legAlice->AddEntry((TObject *)0, Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", 0.0, 30.0), "");
    legAlice->Draw("same");

    TLegend *leg = new TLegend(0.48, 0.57, 0.92, 0.77);
    SetLegendStyle(leg);
    leg->SetTextSize(0.035);
    leg->SetTextFont(42);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hsignal, "Raw K^{0}_{s}K^{0}_{s} invariant mass", "lp");
    leg->AddEntry(hrot, "Rotational background", "lp");
    leg->AddEntry(hmix, "Mixed-event background", "lp");
    leg->AddEntry(hbkg_nopeak, "Norm. region", "f");
    leg->Draw("same");

    c->SaveAs("/home/sawan/Music/compare_rot_mix.png");

    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.03, 0.05, 0.14);
    hrot_signal->GetXaxis()->SetRangeUser(1.00, 2.99);
    SetHistoQA(hrot_signal);
    hrot_signal->SetLineColor(kBlue);
    hrot_signal->SetMarkerColor(kBlue);
    hrot_signal->SetMarkerStyle(23);
    hrot_signal->SetMarkerSize(1.0);
    // hrot_signal->SetMinimum(-40000);
    hrot_signal->Draw("pe");

    SetHistoQA(hmix_signal);
    hmix_signal->GetXaxis()->SetRangeUser(1.00, 2.99);
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

    legAlice->Draw("same");

    TLegend *leg2 = new TLegend(0.25, 0.66, 0.92, 0.77);
    SetLegendStyle(leg2);
    leg2->SetTextSize(0.035);
    leg2->SetTextFont(42);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->AddEntry(hrot_signal, "Rotational background subtraction", "lp");
    leg2->AddEntry(hmix_signal, "Mixed-event background subtraction", "lp");
    leg2->Draw("same");

    c2->SaveAs("/home/sawan/Music/compare_rot_mix_subtracted.png");
}
