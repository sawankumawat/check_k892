#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/common_glue.h"
#include "../src/style.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1); // top pad
    TPad *pad2 = (TPad *)c->GetPad(2); // bottom pad
    pad2Size = 0.3;                    // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.33);
    pad2->SetLeftMargin(0.14);
    pad2->SetTopMargin(0.04);
    pad1->SetRightMargin(0.009);
    pad1->SetTopMargin(0.08);
    pad1->SetLeftMargin(0.14);
    pad1->SetBottomMargin(0.002);
}

void check()
{
    gStyle->SetOptStat(1110);
    TFile *frot = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60_fullpt.root", "READ");
    TFile *fmix = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_MIX_norm_2.50_2.60_fullpt.root", "READ");
    TFile *frot2 = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_ROTATED_norm_2.50_2.60_pt1.root", "READ");
    TFile *fmix2 = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/hglue_MIX_norm_2.50_2.60_pt1.root", "READ");
    if (frot->IsZombie() || fmix->IsZombie() || frot2->IsZombie() || fmix2->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hrot = (TH1F *)frot->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    TH1F *hmix = (TH1F *)fmix->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    TH1F *hrot2 = (TH1F *)frot2->Get("ksks_subtracted_invmass_pt_1.0_30.0");
    TH1F *hmix2 = (TH1F *)fmix2->Get("ksks_subtracted_invmass_pt_1.0_30.0");
    TH1F *hrot_unsubtracted = (TH1F *)frot->Get("ksks_bkg_pt_0.0_30.0");
    TH1F *hmix_unsubtracted = (TH1F *)fmix->Get("ksks_bkg_pt_0.0_30.0");
    TH1F *hrot_unsubtracted2 = (TH1F *)frot2->Get("ksks_subtracted_invmass_pt_1.0_30.0");
    TH1F *hmix_unsubtracted2 = (TH1F *)fmix2->Get("ksks_subtracted_invmass_pt_1.0_30.0");
    TH1F *hinvmass_raw = (TH1F *)frot->Get("ksks_invmass_pt_0.0_30.0");
    TH1F *hinvmass_raw2 = (TH1F *)frot2->Get("ksks_invmass_pt_1.0_30.0");
    TH1F *hrotbkg_wo_normalization = (TH1F *)frot->Get("bkg_without_normalization");
    TH1F *hmixbkg_wo_normalization = (TH1F *)fmix->Get("bkg_without_normalization");
    if (hrot == nullptr || hmix == nullptr || hrot_unsubtracted == nullptr || hmix_unsubtracted == nullptr || hinvmass_raw == nullptr || hrotbkg_wo_normalization == nullptr || hmixbkg_wo_normalization == nullptr || hrot_unsubtracted2 == nullptr || hmix_unsubtracted2 == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }

    TCanvas *crot = (TCanvas *)frot->Get("ksks_invmass_withbkg_pt_0.0_30.0");
    TCanvas *cmix = (TCanvas *)fmix->Get("ksks_invmass_withbkg_pt_0.0_30.0");
    if (crot == nullptr || cmix == nullptr)
    {
        cout << "Error opening canvas" << endl;
        return;
    }
    // Rename the canvases to avoid conflicts
    crot->SetName("crot_canvas");
    cmix->SetName("cmix_canvas");

    TCanvas *c = new TCanvas("", "comparison rotated vs ME", 720, 720);
    SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
    SetHistoQA(hrot);
    SetHistoQA(hmix);
    hrot->GetXaxis()->SetRangeUser(1, 3);
    hmix->GetXaxis()->SetRangeUser(1, 3);
    hrot->GetYaxis()->SetRangeUser(-0.1 * 1e6, 0.5 * 1e6);
    hrot->Draw();
    hmix->SetMarkerColor(2);
    hmix->SetMarkerStyle(22);
    hmix->SetLineColor(2);
    // hmix->Scale(2);
    hmix->Draw("same");

    TLegend *leg = new TLegend(0.20, 0.67, 0.52, 0.92);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry(hrot, "Rotated", "lpe");
    leg->AddEntry(hmix, "Mixed-event", "lpe");
    leg->Draw("same");

    cmix->cd();
    cmix->Draw();

    // Draw on crot canvas
    crot->cd();
    crot->Draw();

    string savepath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/check/";
    c->SaveAs((savepath + "compare_mix_rot_after_subtraction.png").c_str());
    crot->SaveAs((savepath + "rotated_bkg_withsignal.png").c_str());
    cmix->SaveAs((savepath + "mix_bkg_withsignal.png").c_str());

    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.14, 0.03, 0.05, 0.14);
    SetHistoQA(hinvmass_raw);
    SetHistoQA(hrot_unsubtracted);
    SetHistoQA(hmix_unsubtracted);
    hinvmass_raw->GetXaxis()->SetRangeUser(1, 3);
    hrot_unsubtracted->GetXaxis()->SetRangeUser(1, 3);
    hmix_unsubtracted->GetXaxis()->SetRangeUser(1, 3);
    hinvmass_raw->Draw();
    hrot_unsubtracted->SetMarkerColor(4);
    hrot_unsubtracted->SetMarkerStyle(21);
    hrot_unsubtracted->SetLineColor(4);
    hrot_unsubtracted->Draw("same");
    hmix_unsubtracted->SetMarkerColor(2);
    hmix_unsubtracted->SetMarkerStyle(22);
    hmix_unsubtracted->SetLineColor(2);
    // hmix_unsubtracted->Scale(2);
    hmix_unsubtracted->Draw("same");

    TLegend *leg2 = new TLegend(0.52, 0.52, 0.84, 0.83);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.035);
    leg2->SetBorderSize(0);
    leg2->AddEntry(hinvmass_raw, "Raw invariant distribution", "lpe");
    leg2->AddEntry(hrot_unsubtracted, "Rotated", "lpe");
    leg2->AddEntry(hmix_unsubtracted, "Mixed-event", "lpe");
    leg2->Draw("same");
    c2->SaveAs((savepath + "compare_mix_rot_before_subtraction.png").c_str());

    TCanvas *c3 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c3, 0.14, 0.03, 0.08, 0.14);
    double pad1Size, pad2Size;
    canvas_style(c3, pad1Size, pad2Size);
    c3->cd(1);
    SetHistoQA(hrotbkg_wo_normalization);
    SetHistoQA(hmixbkg_wo_normalization);
    hrotbkg_wo_normalization->GetXaxis()->SetRangeUser(1, 3);
    hmixbkg_wo_normalization->GetXaxis()->SetRangeUser(1, 3);
    hmixbkg_wo_normalization->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hmixbkg_wo_normalization->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hmixbkg_wo_normalization->GetYaxis()->SetTitle("Counts");
    hmixbkg_wo_normalization->GetYaxis()->SetTitleOffset(1.3);
    hmixbkg_wo_normalization->Draw();
    hmixbkg_wo_normalization->SetMarkerColor(2);
    hmixbkg_wo_normalization->SetMarkerStyle(22);
    hmixbkg_wo_normalization->SetLineColor(2);
    hrotbkg_wo_normalization->Draw("same");

    TLegend *leg3 = new TLegend(0.52, 0.52, 0.84, 0.83);
    leg3->SetFillStyle(0);
    leg3->SetTextFont(42);
    leg3->SetTextSize(0.035 / pad1Size);
    leg3->SetBorderSize(0);
    leg3->AddEntry(hrotbkg_wo_normalization, "Rotated", "lpe");
    leg3->AddEntry(hmixbkg_wo_normalization, "Mixed-event", "lpe");
    leg3->Draw("same");

    c3->cd(2);
    TH1F *hratio = (TH1F *)hrotbkg_wo_normalization->Clone("hratio");
    hratio->Divide(hmixbkg_wo_normalization);
    hratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hratio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hratio->GetYaxis()->SetTitle("Rotated/ME");
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hratio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hratio->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hratio->GetYaxis()->SetNdivisions(505);
    hratio->GetYaxis()->SetRangeUser(0.8, 1.25);
    hratio->Draw();
    TLine *line = new TLine(1, 1, 3, 1);
    line->SetLineColor(4);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    c3->SaveAs((savepath + "compare_mix_rot_bkg_wo_normalization.png").c_str());

    // //********************************************************************************** */
    // //read the MC file of 13 TeV and compare with the mixed event background to check the shape

    TFile *fmc = new TFile("/home/sawan/check_k892/mc/AnalysisResultsMC_KsKspp_13TeV_nopileupmother.root", "read");
    if (fmc->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    // // print the names of all lists in the root file
    // TIter next(fmc->GetListOfKeys());
    // TKey *key;
    // cout << "The folders in the root file are: \n";
    // while ((key = (TKey *)next()))
    // {
    //     cout << key->GetName() << endl;
    // }

    TList *list = (TList *)fmc->Get("RsnOut_output2023_3.0_3.0_0.0_0.00_12.000_0_1_1_1.0_1.00_0.0_0.000_0.0");
    if (list == nullptr)
    {
        cout << "Error opening list" << endl;
        return;
    }

    // // Iterate in the list to see all objects
    // TIter iterator2(list);
    // TObject *obj;

    // while ((obj = iterator2()))
    // {
    //     cout << "Object Name: " << obj->GetName()
    //          << ", Object Type: " << obj->ClassName() << endl;
    // }

    THnSparseF *hmc1525 = (THnSparseF *)list->FindObject("RsnTaskF0_f0_MotherMCf1525_MIX");
    if (hmc1525 == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }
    TH1F *hmc1525_proj = (TH1F *)hmc1525->Projection(0);
    if (hmc1525_proj == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }
    TCanvas *c4 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c4, 0.14, 0.03, 0.05, 0.14);
    hmc1525_proj->GetXaxis()->SetRangeUser(1, 3);
    SetHistoQA(hmc1525_proj);
    // hmc1525_proj->Scale(5.4); //pt 1-30 GeV/c
    hmc1525_proj->Scale(6.7); // pt 1-30 GeV/c
    hmc1525_proj->Draw();
    hmixbkg_wo_normalization->Draw("same");

    TLegend *leg4 = new TLegend(0.52, 0.6, 0.84, 0.89);
    leg4->SetFillStyle(0);
    leg4->SetTextFont(42);
    leg4->SetTextSize(0.035);
    leg4->SetBorderSize(0);
    leg4->AddEntry(hmc1525_proj, "ME pp 13 TeV MC", "lpe");
    leg4->AddEntry(hmixbkg_wo_normalization, "ME pp13.6 TeV data", "lpe");
    leg4->Draw("same");

    c4->SaveAs((savepath + "compare_mc_mix_bkg_wo_normalization.png").c_str());

    TCanvas *c5 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c5, 0.14, 0.03, 0.05, 0.14);
    hrot->Draw();
    hrot2->GetXaxis()->SetRangeUser(1, 3);
    SetHistoQA(hrot2);
    hrot2->SetMarkerColor(4);
    hrot2->SetLineColor(4);
    hrot2->SetMarkerStyle(21);
    hrot2->Draw("same");

    TLegend *leg5 = new TLegend(0.52, 0.52, 0.84, 0.83);
    leg5->SetFillStyle(0);
    leg5->SetTextFont(42);
    leg5->SetTextSize(0.035);
    leg5->SetBorderSize(0);
    leg5->AddEntry(hrot, "Rotated pt 0-30 GeV/c", "lpe");
    leg5->AddEntry(hrot2, "Rotated pt 1-30 GeV/c", "lpe");
    leg5->Draw("same");
    c5->SaveAs((savepath + "compare_rotated_pt1.png").c_str());

    // TCanvas *c6 = new TCanvas("", "", 720, 720);
    // SetCanvasStyle(c6, 0.14, 0.03, 0.05, 0.14);
}