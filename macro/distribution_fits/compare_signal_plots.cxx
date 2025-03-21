#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/common_glue.h"
#include "../src/style.h"
using namespace std;

void compare_signal_plots()
{
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances";
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/362701/KsKs_Channel/higher-mass-resonances";
    // string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/370825/KsKs_Channel/higher-mass-resonances";
    gStyle->SetOptStat(0);
    int colors[] = {2, 4, 6, 28, kGreen + 2, kBlue -7};
    int markers[] = {20, 21, 22, 23, 26, 27};

    // string path1 = path + "_2sigmaKs/";
    // string path2 = path + "_3sigmaKs/";
    // string path3 = path + "/";
    // string path4 = path + "_5SigmaKs/";

    // const string path1 = path + "_id25081/";
    // const string path2 = path + "_angsep_3_id25081/";
    // const string path3 = path + "_angsep_2_id25081/";
    // const string path4 = path + "_angsep_1p5_id25081/";
    // const string path5 = path + "_angsep_1_id25081/";
    // const string path6 = path + "_angsep_0p5_id25081/";

    // const string path1 = path + "_1Kscut_id24794/";
    // const string path2 = path + "_1p5Kscut_id24794/";
    // const string path3 = path + "_2Kscut_id24794/";
    // const string path4 = path + "_id24794/";
    // const string path5 = path + "_4Kscut_id24794/";

    // const string path1 = path + "_lambda_rej4_id24939/";
    // const string path2 = path + "_id24939/";
    // const string path3 = path + "_lambda_rej6_id24939/";

    // const string path1 = path + "_TPCPID3_id24937/";
    // const string path2 = path + "_TPCPID4_id24937/";
    // const string path3 = path + "_TPCPID6_id24937/";

    // const string path1 = path + "_no_lambda_cut/";
    // const string path2 = path + "_lambda_4/";
    // const string path3 = path + "/";
    // const string path4 = path + "_lambda_6/";
    // const string path5 = path + "_lambda_7/";

    const string path1 = path + "_id24937/";

    TFile *file1 = new TFile((path1 + "hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    TFile *file2 = new TFile((path1 + "hglue_ROTATED_norm_2.50_2.60_pt_0.50_30.00.root").c_str(), "READ");
    TFile *file3 = new TFile((path1 + "hglue_ROTATED_norm_2.50_2.60_pt_1.00_30.00.root").c_str(), "READ");
    TFile *file4 = new TFile((path1 + "hglue_ROTATED_norm_2.50_2.60_pt_2.00_30.00.root").c_str(), "READ");
    TFile *file5 = new TFile((path1 + "hglue_ROTATED_norm_2.50_2.60_pt_0.00_10.00.root").c_str(), "READ");

    // TFile *file1 = new TFile((path1 + "hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // TFile *file2 = new TFile((path2 + "hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // TFile *file3 = new TFile((path3 + "hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // TFile *file4 = new TFile((path4 + "hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // TFile *file5 = new TFile((path5 + "hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // // TFile *file6 = new TFile((path6 + "hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");

    // TFile *file1temp = new TFile((path1 + "backup/hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // TFile *file2temp = new TFile((path2 + "backup/hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // TFile *file3temp = new TFile((path3 + "backup/hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // TFile *file4temp = new TFile((path4 + "backup/hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // TFile *file5temp = new TFile((path5 + "backup/hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");
    // TFile *file6temp = new TFile((path6 + "backup/hglue_ROTATED_norm_2.50_2.60_pt_0.00_30.00.root").c_str(), "READ");

    if (file1->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    TH1F *h1 = (TH1F *)file1->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    TH1F *h2 = (TH1F *)file2->Get("ksks_subtracted_invmass_pt_0.5_30.0");
    TH1F *h3 = (TH1F *)file3->Get("ksks_subtracted_invmass_pt_1.0_30.0");
    TH1F *h4 = (TH1F *)file4->Get("ksks_subtracted_invmass_pt_2.0_30.0");
    TH1F *h5 = (TH1F *)file5->Get("ksks_subtracted_invmass_pt_0.0_10.0");


    // TH1F *h1 = (TH1F *)file1->Get("ksks_invmass_pt_0.0_30.0");
    // TH1F *h2 = (TH1F *)file2->Get("ksks_invmass_pt_0.0_30.0");
    // TH1F *h3 = (TH1F *)file3->Get("ksks_invmass_pt_0.0_30.0");
    // TH1F *h4 = (TH1F *)file4->Get("ksks_invmass_pt_0.0_30.0");
    // TH1F *h5 = (TH1F *)file5->Get("ksks_invmass_pt_0.0_30.0");
    // TH1F *h6 = (TH1F *)file6->Get("ksks_invmass_pt_0.0_30.0");

    // h1->Rebin(2);
    // h2->Rebin(2);
    // h3->Rebin(2);
    // h4->Rebin(2);
    // h5->Rebin(2);
    // // h6->Rebin(2);

    // TH1F *hbkg1 = (TH1F *)file1->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg2 = (TH1F *)file2->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg3 = (TH1F *)file3->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg4 = (TH1F *)file4->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg5 = (TH1F *)file5->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg6 = (TH1F *)file6->Get("ksks_bkg_pt_0.0_30.0");

    // hbkg1->Rebin(2);
    // hbkg2->Rebin(2);
    // hbkg3->Rebin(2);
    // hbkg4->Rebin(2);
    // hbkg5->Rebin(2);
    // hbkg6->Rebin(2);

    // TH1F *hbkg1 = (TH1F *)file1->Get("bkg_without_normalization");
    // TH1F *hbkg2 = (TH1F *)file2->Get("bkg_without_normalization");
    // TH1F *hbkg3 = (TH1F *)file3->Get("bkg_without_normalization");
    // TH1F *hbkg4 = (TH1F *)file4->Get("bkg_without_normalization");
    // TH1F *hbkg5 = (TH1F *)file5->Get("bkg_without_normalization");
    // TH1F *hbkg6 = (TH1F *)file6->Get("bkg_without_normalization");

    // TH1F *hbkg1temp = (TH1F *)file1temp->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg2temp = (TH1F *)file2temp->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg3temp = (TH1F *)file3temp->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg4temp = (TH1F *)file4temp->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg5temp = (TH1F *)file5temp->Get("ksks_bkg_pt_0.0_30.0");
    // TH1F *hbkg6temp = (TH1F *)file6temp->Get("ksks_bkg_pt_0.0_30.0");

    if (h1 == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.16, 0.03, 0.05, 0.14);
    h1->GetXaxis()->SetRangeUser(1.00, 2.60);
    SetHistoQA(h1);
    h1->SetMarkerSize(0.8);
    h1->SetMarkerColor(colors[0]);
    h1->SetLineColor(colors[0]);
    h1->SetMarkerStyle(markers[0]);
    h1->GetYaxis()->SetTitleOffset(1.6);
    h1->GetYaxis()->SetTitle(Form("Counts/%.2f GeV/c^{2}", h1->GetXaxis()->GetBinWidth(1)));
    h1->SetMinimum(-10000);
    h1->SetMaximum(h1->GetMaximum() * 1.2);
    h1->Draw();
    SetHistoQA(h2);
    h2->GetXaxis()->SetRangeUser(1.00, 2.60);
    h2->SetMarkerSize(0.8);
    h2->SetMarkerColor(colors[1]);
    h2->SetLineColor(colors[1]);
    h2->SetMarkerStyle(markers[1]);
    h2->Draw("same");
    SetHistoQA(h3);
    h3->GetXaxis()->SetRangeUser(1.00, 2.60);
    h3->SetMarkerSize(0.8);
    h3->SetMarkerColor(colors[2]);
    h3->SetLineColor(colors[2]);
    h3->SetMarkerStyle(markers[2]);
    h3->Draw("same");
    SetHistoQA(h4);
    h4->GetXaxis()->SetRangeUser(1.00, 2.60);
    h4->SetMarkerSize(0.8);
    h4->SetMarkerColor(colors[3]);
    h4->SetLineColor(colors[3]);
    h4->SetMarkerStyle(markers[3]);
    h4->Draw("same");
    SetHistoQA(h5);
    h5->GetXaxis()->SetRangeUser(1.00, 2.60);
    h5->SetMarkerColor(colors[4]);
    h5->SetLineColor(colors[4]);
    h5->SetMarkerStyle(markers[4]);
    h5->Draw("same");
    // SetHistoQA(h6);
    // h6->GetXaxis()->SetRangeUser(1.00, 2.60);
    // h6->SetMarkerColor(colors[5]);
    // h6->SetLineColor(colors[5]);
    // h6->SetMarkerStyle(markers[5]);
    // h6->Draw("same");

    // TLegend *leg = new TLegend(0.55, 0.62, 0.92, 0.92);
    TLegend *leg = new TLegend(0.4, 0.55, 0.92, 0.92);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    // leg->SetHeader("Normalization ranges");
    // leg->SetHeader("K^{0}_{s} selection mass window");
    // leg->SetHeader("Angular separation cut");
    // leg->AddEntry(h1, "No cut", "lpe");
    // leg->AddEntry(h2, "3.0", "lpe");
    // leg->AddEntry(h3, "2.0", "lpe");
    // leg->AddEntry(h4, "1.5", "lpe");
    // leg->AddEntry(h5, "1.0", "lpe");
    // leg->AddEntry(h6, "0.5", "lpe");

    // leg->AddEntry(h1, "1#sigma", "lpe");
    // leg->AddEntry(h2, "1.5#sigma", "lpe");
    // leg->AddEntry(h3, "2#sigma", "lpe");
    // leg->AddEntry(h4, "3#sigma", "lpe");
    // leg->AddEntry(h5, "4#sigma", "lpe");

    // leg->SetHeader("Lambda rejection");
    // leg->AddEntry(h1, "No cut", "lpe");
    // leg->AddEntry(h2, "4 MeV/c^{2}", "lpe");
    // leg->AddEntry(h3, "5 MeV/c^{2}", "lpe");
    // leg->AddEntry(h4, "6 MeV/c^{2}", "lpe");
    // leg->AddEntry(h5, "7 MeV/c^{2}", "lpe");

    // leg->SetHeader("TPC PID");
    // leg->AddEntry(h1, "3 #sigma", "lpe");
    // leg->AddEntry(h2, "4 #sigma", "lpe");
    // leg->AddEntry(h3, "5 #sigma", "lpe");

    leg->SetHeader("p_{T} ranges (statistics lost)");
    leg->AddEntry(h1, Form("0.0 - 30.0 GeV/c (%.1f %)", (h1->GetEntries() - h1->GetEntries()) * 100 / h1->GetEntries()), "lpe");
    leg->AddEntry(h2, Form("0.5 - 30.0 GeV/c (%.1f %)", (h1->GetEntries() - h2->GetEntries()) * 100 / h1->GetEntries()), "lpe");
    leg->AddEntry(h3, Form("1.0 - 30.0 GeV/c (%.1f %)", (h1->GetEntries() - h3->GetEntries()) * 100 / h1->GetEntries()), "lpe");
    leg->AddEntry(h4, Form("2.0 - 30.0 GeV/c (%.1f %)", (h1->GetEntries() - h4->GetEntries()) * 100 / h1->GetEntries()), "lpe");
    leg->AddEntry(h5, Form("0.0 - 10.0 GeV/c (%.1f %)", (h1->GetEntries() - h5->GetEntries()) * 100 / h1->GetEntries()), "lpe");
    leg->Draw("same");

    // cout << "entries in h1 " << h1->GetEntries() << endl;
    // cout << "entries in h2 " << h2->GetEntries() << endl;
    // cout << "entries in h3 " << h3->GetEntries() << endl;
    // cout << "entries in h4 " << h4->GetEntries() << endl;
    // cout << "entries in h5 " << h5->GetEntries() << endl;

    // // // c->SaveAs((path1 + "compare_norm_ranges.png").c_str());
    c->SaveAs("/home/sawan/Music/compare_rotpt.png");
    // c->SaveAs("/home/sawan/Music/compare_TPCPIDSignal.png");

    // TCanvas *c2 = new TCanvas("", "", 720, 720);
    // SetCanvasStyle(c2, 0.14, 0.03, 0.05, 0.14);
    // hbkg1->GetXaxis()->SetRangeUser(1.00, 2.60);
    // SetHistoQA(hbkg1);
    // hbkg1->SetMarkerSize(0.8);
    // hbkg1->SetMarkerColor(colors[0]);
    // hbkg1->SetLineColor(colors[0]);
    // hbkg1->SetMarkerStyle(markers[0]);
    // hbkg1->GetYaxis()->SetTitle(Form("Counts/%.2f GeV/c^{2}", hbkg1->GetXaxis()->GetBinWidth(1)));
    // hbkg1->SetMinimum(0);
    // hbkg1->SetMaximum(hbkg1->GetMaximum() * 1.4);
    // hbkg1->Draw();
    // SetHistoQA(hbkg2);
    // hbkg2->GetXaxis()->SetRangeUser(1.00, 2.60);
    // hbkg2->SetMarkerSize(0.8);
    // hbkg2->SetMarkerColor(colors[1]);
    // hbkg2->SetLineColor(colors[1]);
    // hbkg2->SetMarkerStyle(markers[1]);
    // hbkg2->Draw("same");
    // SetHistoQA(hbkg3);
    // hbkg3->GetXaxis()->SetRangeUser(1.00, 2.60);
    // hbkg3->SetMarkerSize(0.8);
    // hbkg3->SetMarkerColor(colors[2]);
    // hbkg3->SetLineColor(colors[2]);
    // hbkg3->SetMarkerStyle(markers[2]);
    // hbkg3->Draw("same");
    // SetHistoQA(hbkg4);
    // hbkg4->GetXaxis()->SetRangeUser(1.00, 2.60);
    // hbkg4->SetMarkerSize(0.8);
    // hbkg4->SetMarkerColor(colors[3]);
    // hbkg4->SetLineColor(colors[3]);
    // hbkg4->SetMarkerStyle(markers[3]);
    // hbkg4->Draw("same");
    // SetHistoQA(hbkg5);
    // hbkg5->GetXaxis()->SetRangeUser(1.00, 2.60);
    // hbkg5->SetMarkerColor(colors[4]);
    // hbkg5->SetLineColor(colors[4]);
    // hbkg5->SetMarkerStyle(markers[4]);
    // hbkg5->Draw("same");
    // // SetHistoQA(hbkg6);
    // // hbkg6->GetXaxis()->SetRangeUser(1.00, 2.60);
    // // hbkg6->SetMarkerColor(colors[5]);
    // // hbkg6->SetLineColor(colors[5]);
    // // hbkg6->SetMarkerStyle(markers[5]);
    // // hbkg6->Draw("same");

    // TLegend *leg2 = new TLegend(0.50, 0.62, 0.92, 0.92);
    // leg2->SetFillStyle(0);
    // leg2->SetTextFont(42);
    // leg2->SetTextSize(0.04);
    // leg2->SetBorderSize(0);
    // // leg2->SetHeader("Normalization ranges");
    // // leg2->SetHeader("K^{0}_{s} selection mass window");
    // // leg2->SetHeader("Angular separation cut");
    // // leg2->AddEntry(hbkg1, "No cut", "lpe");
    // // leg2->AddEntry(hbkg2, "3.0", "lpe");
    // // leg2->AddEntry(hbkg3, "2.0", "lpe");
    // // leg2->AddEntry(hbkg4, "1.5", "lpe");
    // // leg2->AddEntry(hbkg5, "1.0", "lpe");
    // // leg2->AddEntry(hbkg6, "0.5", "lpe");
    // // leg2->AddEntry(hbkg1, "1#sigma", "lpe");
    // // leg2->AddEntry(hbkg2, "1.5#sigma", "lpe");
    // // leg2->AddEntry(hbkg3, "2#sigma", "lpe");
    // // leg2->AddEntry(hbkg4, "3#sigma", "lpe");
    // // leg2->AddEntry(hbkg5, "4#sigma", "lpe");
    // leg2->SetHeader("Lambda rejection");
    // leg2->AddEntry(hbkg1, "No cut", "lpe");
    // leg2->AddEntry(hbkg2, "4 MeV/c^{2}", "lpe");
    // leg2->AddEntry(hbkg3, "5 MeV/c^{2}", "lpe");
    // leg2->AddEntry(hbkg4, "6 MeV/c^{2}", "lpe");
    // leg2->AddEntry(hbkg5, "7 MeV/c^{2}", "lpe");
    // leg2->Draw("same");

    // // // c2->SaveAs((path1 + "compare_norm_ranges_bkg.png").c_str());
    // c2->SaveAs("/home/sawan/Music/compare_LambdarejBkg.png");

    // TCanvas *call = new TCanvas("", "", 1440, 720);
    // SetCanvasStyle(call, 0.14, 0.03, 0.05, 0.14);
    // call->Divide(3, 2);
    // TLatex lat;
    // lat.SetNDC();
    // lat.SetTextSize(0.05);
    // lat.SetTextFont(42);
    // for (int i = 0; i < 6; i++)
    // {
    //     call->cd(i + 1);
    //     gPad->SetTickx(1);
    //     gPad->SetTicky(1);
    //     gPad->SetRightMargin(0.03);
    //     gPad->SetTopMargin(0.05);
    //     gPad->SetBottomMargin(0.15);
    //     gPad->SetLeftMargin(0.15);
    //     if (i == 0)
    //     {
    //         h1->SetMarkerColor(1);
    //         h1->SetLineColor(1);
    //         h1->SetMarkerStyle(20);
    //         hbkg1->SetMarkerColor(2);
    //         h1->Draw();
    //         hbkg1->Draw("same");
    //         hbkg1temp->SetMarkerColor(4);
    //         hbkg1temp->SetMarkerStyle(22);
    //         hbkg1temp->SetLineColor(4);
    //         hbkg1temp->Draw("same");

    //         TLegend *leg3 = new TLegend(0.22, 0.25, 0.52, 0.55);
    //         leg3->SetFillStyle(0);
    //         leg3->SetTextFont(42);
    //         leg3->SetTextSize(0.05);
    //         leg3->SetBorderSize(0);
    //         leg3->AddEntry(h1, "Signal", "lpe");
    //         leg3->AddEntry(hbkg1temp, "Rotational background", "lpe");
    //         leg3->AddEntry(hbkg1, "Rotational background (r cut)", "lpe");
    //         leg3->Draw("same");
    //         lat.DrawLatex(0.65, 0.80, "No angular cut");
    //     }
    //     if (i == 1)
    //     {
    //         h2->SetMinimum(0);
    //         h2->SetMarkerColor(1);
    //         h2->SetMarkerStyle(20);
    //         hbkg2->SetMarkerColor(2);
    //         h2->Draw();
    //         hbkg2->Draw("same");
    //         hbkg2temp->SetMarkerColor(4);
    //         hbkg2temp->SetMarkerStyle(22);
    //         hbkg2temp->Draw("same");
    //         lat.DrawLatex(0.7, 0.80, "3.0");
    //     }
    //     if (i == 2)
    //     {
    //         h3->SetMinimum(0);
    //         h3->SetMarkerColor(1);
    //         h3->SetMarkerStyle(20);
    //         hbkg3->SetMarkerColor(2);
    //         h3->Draw();
    //         hbkg3->Draw("same");
    //         hbkg3temp->SetMarkerColor(4);
    //         hbkg3temp->SetMarkerStyle(22);
    //         hbkg3temp->Draw("same");
    //         lat.DrawLatex(0.7, 0.80, "2.0");
    //     }
    //     if (i == 3)
    //     {
    //         h4->SetMinimum(0);
    //         h4->SetMarkerColor(1);
    //         h4->SetMarkerStyle(20);
    //         hbkg4->SetMarkerColor(2);
    //         h4->Draw();
    //         hbkg4->Draw("same");
    //         hbkg4temp->SetMarkerColor(4);
    //         hbkg4temp->SetMarkerStyle(22);
    //         hbkg4temp->Draw("same");
    //         lat.DrawLatex(0.7, 0.80, "1.5");
    //     }
    //     if (i == 4)
    //     {
    //         h5->SetMinimum(0);
    //         h5->SetMarkerColor(1);
    //         h5->SetMarkerStyle(20);
    //         hbkg5->SetMarkerColor(2);
    //         h5->Draw();
    //         hbkg5->Draw("same");
    //         hbkg5temp->SetMarkerColor(4);
    //         hbkg5temp->SetMarkerStyle(22);
    //         hbkg5temp->Draw("same");
    //         lat.DrawLatex(0.7, 0.80, "1.0");
    //     }
    //     if (i == 5)
    //     {
    //         h6->SetMinimum(0);
    //         h6->SetMarkerColor(1);
    //         h6->SetMarkerStyle(20);
    //         hbkg6->SetMarkerColor(2);
    //         h6->Draw();
    //         hbkg6->Draw("same");
    //         hbkg6temp->SetMarkerColor(4);
    //         hbkg6temp->SetMarkerStyle(22);
    //         hbkg6temp->Draw("same");
    //         lat.DrawLatex(0.7, 0.80, "0.5");
    //     }
    // }

    // call->SaveAs("/home/sawan/Music/compare_signal_bkg_angsep_all.png");
}