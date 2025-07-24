#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

void compare_QAplots()
{
    gStyle->SetOptStat(0);
    string datapath = "/home/sawan/check_k892/data/kstar/LHC22o_pass7/";
    string path1 = "448490"; // 2022 data
    string path2 = "449695"; // 2023 data
    // string path2 = "451003"; // 2024 data
    TString outputPath = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/" + path2 + "/kstarqa/hInvMass/QA";

    TFile *file1 = new TFile((datapath + path1 + ".root").c_str(), "read");
    TFile *file2 = new TFile((datapath + path2 + ".root").c_str(), "read");
    if (file1->IsZombie() || file2->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    TH1F *hNsigmaTOFKaon_neg1 = (TH1F *)file1->Get("kstarqa/hPID/Before/h1PID_TOF_neg_kaon");
    TH1F *hNsigmaTOFPion_neg1 = (TH1F *)file1->Get("kstarqa/hPID/Before/h1PID_TOF_neg_pion");
    TH1F *hNsigmaTOFKaon_pos1 = (TH1F *)file1->Get("kstarqa/hPID/Before/h1PID_TOF_pos_kaon");
    TH1F *hNsigmaTOFPion_pos1 = (TH1F *)file1->Get("kstarqa/hPID/Before/h1PID_TOF_pos_pion");
    TH1F *hNsigmaTPCKaon_neg1 = (TH1F *)file1->Get("kstarqa/hPID/Before/h1PID_TPC_neg_kaon");
    TH1F *hNsigmaTPCPion_neg1 = (TH1F *)file1->Get("kstarqa/hPID/Before/h1PID_TPC_neg_pion");
    TH1F *hNsigmaTPCKaon_pos1 = (TH1F *)file1->Get("kstarqa/hPID/Before/h1PID_TPC_pos_kaon");
    TH1F *hNsigmaTPCPion_pos1 = (TH1F *)file1->Get("kstarqa/hPID/Before/h1PID_TPC_pos_pion");

    TH1F *hNsigmaTOFKaon_neg2 = (TH1F *)file2->Get("kstarqa/hPID/Before/h1PID_TOF_neg_kaon");
    TH1F *hNsigmaTOFPion_neg2 = (TH1F *)file2->Get("kstarqa/hPID/Before/h1PID_TOF_neg_pion");
    TH1F *hNsigmaTOFKaon_pos2 = (TH1F *)file2->Get("kstarqa/hPID/Before/h1PID_TOF_pos_kaon");
    TH1F *hNsigmaTOFPion_pos2 = (TH1F *)file2->Get("kstarqa/hPID/Before/h1PID_TOF_pos_pion");
    TH1F *hNsigmaTPCKaon_neg2 = (TH1F *)file2->Get("kstarqa/hPID/Before/h1PID_TPC_neg_kaon");
    TH1F *hNsigmaTPCPion_neg2 = (TH1F *)file2->Get("kstarqa/hPID/Before/h1PID_TPC_neg_pion");
    TH1F *hNsigmaTPCKaon_pos2 = (TH1F *)file2->Get("kstarqa/hPID/Before/h1PID_TPC_pos_kaon");
    TH1F *hNsigmaTPCPion_pos2 = (TH1F *)file2->Get("kstarqa/hPID/Before/h1PID_TPC_pos_pion");

    hNsigmaTOFKaon_neg1->Scale(1.0 / hNsigmaTOFKaon_neg1->GetEntries());
    hNsigmaTOFPion_neg1->Scale(1.0 / hNsigmaTOFPion_neg1->GetEntries());
    hNsigmaTOFKaon_pos1->Scale(1.0 / hNsigmaTOFKaon_pos1->GetEntries());
    hNsigmaTOFPion_pos1->Scale(1.0 / hNsigmaTOFPion_pos1->GetEntries());
    hNsigmaTPCKaon_neg1->Scale(1.0 / hNsigmaTPCKaon_neg1->GetEntries());
    hNsigmaTPCPion_neg1->Scale(1.0 / hNsigmaTPCPion_neg1->GetEntries());
    hNsigmaTPCKaon_pos1->Scale(1.0 / hNsigmaTPCKaon_pos1->GetEntries());
    hNsigmaTPCPion_pos1->Scale(1.0 / hNsigmaTPCPion_pos1->GetEntries());

    hNsigmaTOFKaon_neg2->Scale(1.0 / hNsigmaTOFKaon_neg2->GetEntries());
    hNsigmaTOFPion_neg2->Scale(1.0 / hNsigmaTOFPion_neg2->GetEntries());
    hNsigmaTOFKaon_pos2->Scale(1.0 / hNsigmaTOFKaon_pos2->GetEntries());
    hNsigmaTOFPion_pos2->Scale(1.0 / hNsigmaTOFPion_pos2->GetEntries());
    hNsigmaTPCKaon_neg2->Scale(1.0 / hNsigmaTPCKaon_neg2->GetEntries());
    hNsigmaTPCPion_neg2->Scale(1.0 / hNsigmaTPCPion_neg2->GetEntries());
    hNsigmaTPCKaon_pos2->Scale(1.0 / hNsigmaTPCKaon_pos2->GetEntries());
    hNsigmaTPCPion_pos2->Scale(1.0 / hNsigmaTPCPion_pos2->GetEntries());

    if (hNsigmaTOFKaon_neg1 == nullptr)
    {
        cerr << "PID histograms not found in file 1!!!!!!!!!!!!" << endl;
        return;
    }
    if (hNsigmaTOFKaon_neg2 == nullptr)
    {
        cerr << "PID histograms not found in file 2!!!!!!!!!!!!" << endl;
        return;
    }

    TCanvas *cNsigmaTOFKaon_neg = new TCanvas("cNsigmaTOFKaon_neg", "Nsigma TOF Kaon Neg", 720, 720);
    SetCanvasStyle(cNsigmaTOFKaon_neg, 0.14, 0.05, 0.06, 0.14);
    SetHistoQA(hNsigmaTOFKaon_neg1);
    SetHistoQA(hNsigmaTOFKaon_neg2);
    hNsigmaTOFKaon_neg1->GetYaxis()->SetTitle("Counts");
    hNsigmaTOFKaon_neg1->GetXaxis()->SetTitle("n#sigma_{TOF} (K)");
    hNsigmaTOFKaon_neg1->GetXaxis()->SetRangeUser(-3, 3);
    hNsigmaTOFKaon_neg1->SetMaximum(1.3 * hNsigmaTOFKaon_neg1->GetMaximum());
    hNsigmaTOFKaon_neg1->Draw("HIST");
    hNsigmaTOFKaon_neg2->SetLineColor(kRed);
    hNsigmaTOFKaon_neg2->Draw("HIST same");
    TLegend *legend = new TLegend(0.2, 0.75, 0.5, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(hNsigmaTOFKaon_neg1, "2022 dataset", "l");
    legend->AddEntry(hNsigmaTOFKaon_neg2, "2023 dataset", "l");
    legend->Draw();
    TLine *lineverticalx0 = new TLine(0, hNsigmaTOFKaon_neg1->GetYaxis()->GetXmin(), 0, hNsigmaTOFKaon_neg1->GetYaxis()->GetXmax());
    lineverticalx0->SetLineStyle(2);
    lineverticalx0->SetLineColor(2);
    lineverticalx0->SetLineWidth(2);
    // lineverticalx0->Draw();
    cNsigmaTOFKaon_neg->SaveAs(outputPath + ("/compare_NsigmaTOFKaon_neg.png"));

    TCanvas *cNsigmaTOFPion_neg = new TCanvas("cNsigmaTOFPion_neg", "Nsigma TOF Pion Neg", 720, 720);
    SetCanvasStyle(cNsigmaTOFPion_neg, 0.14, 0.05, 0.06, 0.14);
    SetHistoQA(hNsigmaTOFPion_neg1);
    SetHistoQA(hNsigmaTOFPion_neg2);
    hNsigmaTOFPion_neg1->GetYaxis()->SetTitle("Counts");
    hNsigmaTOFPion_neg1->GetXaxis()->SetTitle("n#sigma_{TOF} (#pi)");
    hNsigmaTOFPion_neg1->GetXaxis()->SetRangeUser(-3, 3);
    hNsigmaTOFPion_neg1->SetMaximum(1.3 * hNsigmaTOFKaon_neg1->GetMaximum());
    hNsigmaTOFPion_neg1->Draw("HIST");
    hNsigmaTOFPion_neg2->SetLineColor(kRed);
    hNsigmaTOFPion_neg2->Draw("HIST same");
    legend->Draw();
    cNsigmaTOFPion_neg->SaveAs(outputPath + ("/compare_NsigmaTOFPion_neg.png"));

    TCanvas *cNsigmaTOFKaon_pos = new TCanvas("cNsigmaTOFKaon_pos", "Nsigma TOF Kaon Pos", 720, 720);
    SetCanvasStyle(cNsigmaTOFKaon_pos, 0.14, 0.05, 0.06, 0.14);
    SetHistoQA(hNsigmaTOFKaon_pos1);
    SetHistoQA(hNsigmaTOFKaon_pos2);
    hNsigmaTOFKaon_pos1->GetYaxis()->SetTitle("Counts");
    hNsigmaTOFKaon_pos1->GetXaxis()->SetTitle("n#sigma_{TOF} (K)");
    hNsigmaTOFKaon_pos1->GetXaxis()->SetRangeUser(-3, 3);
    hNsigmaTOFKaon_pos1->SetMaximum(1.3 * hNsigmaTOFKaon_pos1->GetMaximum());
    hNsigmaTOFKaon_pos1->Draw("HIST");
    hNsigmaTOFKaon_pos2->SetLineColor(kRed);
    hNsigmaTOFKaon_pos2->Draw("HIST same");
    legend->Draw();
    cNsigmaTOFKaon_pos->SaveAs(outputPath + ("/compare_NsigmaTOFKaon_pos.png"));

    TCanvas *cNsigmaTOFPion_pos = new TCanvas("cNsigmaTOFPion_pos", "Nsigma TOF Pion Pos", 720, 720);
    SetCanvasStyle(cNsigmaTOFPion_pos, 0.14, 0.05, 0.06, 0.14);
    SetHistoQA(hNsigmaTOFPion_pos1);
    SetHistoQA(hNsigmaTOFPion_pos2);
    hNsigmaTOFPion_pos1->GetYaxis()->SetTitle("Counts");
    hNsigmaTOFPion_pos1->GetXaxis()->SetTitle("n#sigma_{TOF} (#pi)");
    hNsigmaTOFPion_pos1->GetXaxis()->SetRangeUser(-3, 3);
    hNsigmaTOFPion_pos1->SetMaximum(1.3 * hNsigmaTOFPion_pos1->GetMaximum());
    hNsigmaTOFPion_pos1->Draw("HIST");
    hNsigmaTOFPion_pos2->SetLineColor(kRed);
    hNsigmaTOFPion_pos2->Draw("HIST same");
    legend->Draw();
    cNsigmaTOFPion_pos->SaveAs(outputPath + ("/compare_NsigmaTOFPion_pos.png"));

    TCanvas *cNsigmaTPCKaon_neg = new TCanvas("cNsigmaTPCKaon_neg", "Nsigma TPC Kaon Neg", 720, 720);
    SetCanvasStyle(cNsigmaTPCKaon_neg, 0.14, 0.05, 0.06, 0.14);
    SetHistoQA(hNsigmaTPCKaon_neg1);
    SetHistoQA(hNsigmaTPCKaon_neg2);
    hNsigmaTPCKaon_neg1->GetYaxis()->SetTitle("Counts");
    hNsigmaTPCKaon_neg1->GetXaxis()->SetTitle("n#sigma_{TPC} (K)");
    hNsigmaTPCKaon_neg1->GetXaxis()->SetRangeUser(-3, 3);
    hNsigmaTPCKaon_neg1->SetMaximum(1.3 * hNsigmaTPCKaon_neg1->GetMaximum());
    hNsigmaTPCKaon_neg1->Draw("HIST");
    hNsigmaTPCKaon_neg2->SetLineColor(kRed);
    hNsigmaTPCKaon_neg2->Draw("HIST same");
    legend->Draw();
    cNsigmaTPCKaon_neg->SaveAs(outputPath + ("/compare_NsigmaTPCKaon_neg.png"));

    TCanvas *cNsigmaTPCPion_neg = new TCanvas("cNsigmaTPCPion_neg", "Nsigma TPC Pion Neg", 720, 720);
    SetCanvasStyle(cNsigmaTPCPion_neg, 0.14, 0.05, 0.06, 0.14);
    SetHistoQA(hNsigmaTPCPion_neg1);
    SetHistoQA(hNsigmaTPCPion_neg2);
    hNsigmaTPCPion_neg1->GetYaxis()->SetTitle("Counts");
    hNsigmaTPCPion_neg1->GetXaxis()->SetTitle("n#sigma_{TPC} (#pi)");
    hNsigmaTPCPion_neg1->GetXaxis()->SetRangeUser(-3, 3);
    hNsigmaTPCPion_neg1->SetMaximum(1.3 * hNsigmaTPCPion_neg1->GetMaximum());
    hNsigmaTPCPion_neg1->Draw("HIST");
    hNsigmaTPCPion_neg2->SetLineColor(kRed);
    hNsigmaTPCPion_neg2->Draw("HIST same");
    legend->Draw();
    cNsigmaTPCPion_neg->SaveAs(outputPath + ("/compare_NsigmaTPCPion_neg.png"));

    TCanvas *cNsigmaTPCKaon_pos = new TCanvas("cNsigmaTPCKaon_pos", "Nsigma TPC Kaon Pos", 720, 720);
    SetCanvasStyle(cNsigmaTPCKaon_pos, 0.14, 0.05, 0.06, 0.14);
    SetHistoQA(hNsigmaTPCKaon_pos1);
    SetHistoQA(hNsigmaTPCKaon_pos2);
    hNsigmaTPCKaon_pos1->GetYaxis()->SetTitle("Counts");
    hNsigmaTPCKaon_pos1->GetXaxis()->SetTitle("n#sigma_{TPC} (K)");
    hNsigmaTPCKaon_pos1->GetXaxis()->SetRangeUser(-3, 3);
    hNsigmaTPCKaon_pos1->SetMaximum(1.3 * hNsigmaTPCKaon_pos1->GetMaximum());
    hNsigmaTPCKaon_pos1->Draw("HIST");
    hNsigmaTPCKaon_pos2->SetLineColor(kRed);
    hNsigmaTPCKaon_pos2->Draw("HIST same");
    legend->Draw();
    cNsigmaTPCKaon_pos->SaveAs(outputPath + ("/compare_NsigmaTPCKaon_pos.png"));

    TCanvas *cNsigmaTPCPion_pos = new TCanvas("cNsigmaTPCPion_pos", "Nsigma TPC Pion Pos", 720, 720);
    SetCanvasStyle(cNsigmaTPCPion_pos, 0.14, 0.05, 0.06, 0.14);
    SetHistoQA(hNsigmaTPCPion_pos1);
    SetHistoQA(hNsigmaTPCPion_pos2);
    hNsigmaTPCPion_pos1->GetYaxis()->SetTitle("Counts");
    hNsigmaTPCPion_pos1->GetXaxis()->SetTitle("n#sigma_{TPC} (#pi)");
    hNsigmaTPCPion_pos1->GetXaxis()->SetRangeUser(-3, 3);
    hNsigmaTPCPion_pos1->SetMaximum(1.3 * hNsigmaTPCPion_pos1->GetMaximum());
    hNsigmaTPCPion_pos1->Draw("HIST");
    hNsigmaTPCPion_pos2->SetLineColor(kRed);
    hNsigmaTPCPion_pos2->Draw("HIST same");
    legend->Draw();
    cNsigmaTPCPion_pos->SaveAs(outputPath + ("/compare_NsigmaTPCPion_pos.png"));
}