#include <iostream>
using namespace std;
#include "style.h"

void plots()
{
    gStyle->SetOptStat(1110);
    // gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // TFile *f = new TFile("/home/sawan/check_k892/mc/LHC24l1/316063_onlypass7.root", "read");
    // TFile *f = new TFile("/home/sawan/check_k892/mc/LHC24l1/337948.root", "read");
    // TFile *f = new TFile("/home/sawan/check_k892/mc/LHC24l1/340531.root", "read");
    TFile *f = new TFile("AnalysisResults_new2.root", "read");

    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string histpath = "higher-mass-resonances/hMChists";

    THnSparseD *GenpT = (THnSparseD *)f->Get(Form("%s/Genf1710", histpath.c_str()));
    THnSparseD *recpt1 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt1", histpath.c_str()));
    THnSparseD *recpt2 = (THnSparseD *)f->Get(Form("%s/Recf1710_pt2", histpath.c_str()));

    TH1F *recMass = (TH1F *)f->Get(Form("%s/Recf1710_mass", histpath.c_str()));
    TH1F *GenMass = (TH1F *)f->Get(Form("%s/Genf1710_mass", histpath.c_str()));
    TH1F *MC_mult = (TH1F *)f->Get(Form("%s/MC_mult", histpath.c_str()));
    TH1F *MC_mult_sel = (TH1F *)f->Get(Form("%s/MC_mult_after_event_sel", histpath.c_str()));
    TH1F *genPhi = (TH1F *)f->Get(Form("%s/GenPhi", histpath.c_str()));
    TH1F *recPhi = (TH1F *)f->Get(Form("%s/RecPhi", histpath.c_str()));
    TH1F *recRap = (TH1F *)f->Get(Form("%s/RecRapidity", histpath.c_str()));
    TH1F *recEta = (TH1F *)f->Get(Form("%s/RecEta", histpath.c_str()));
    THnSparseF *genRap = (THnSparseF *)f->Get(Form("%s/GenRapidity", histpath.c_str()));
    THnSparseF *genEta = (THnSparseF *)f->Get(Form("%s/GenEta", histpath.c_str()));

    // recpt1->Rebin(4);
    // recpt2->Rebin(4);
    // GenpT->Rebin(4);

    // cout << "GenpT: " << " Low edge: " << GenpT->GetBinLowEdge(1) << " High edge: " << GenpT->GetBinLowEdge(GenpT->GetNbinsX() + 1) << " Bin width: " << GenpT->GetBinWidth(1) << endl;
    // cout << "RecpT1: " << " Low edge: " << recpt1->GetBinLowEdge(1) << " High edge: " << recpt1->GetBinLowEdge(recpt1->GetNbinsX() + 1) << " Bin width: " << recpt1->GetBinWidth(1) << endl;

    if (GenpT == nullptr || GenMass == nullptr || MC_mult == nullptr || MC_mult_sel == nullptr || recpt1 == nullptr || recpt2 == nullptr || recMass == nullptr || genPhi == nullptr || recPhi == nullptr || recRap == nullptr || recEta == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    vector<string> histnames = {"Gen Mass f_{2}(1525)", "Multiplicity", "Multiplicity after event selection", "Rec f_{2}(1525) Mass", "Gen Phi", "Gen Rapidity", "Gen Eta", "Rec Phi", "Rec Rapidity", "Rec Eta"};
    vector<TH1F *> hists = {GenMass, MC_mult, MC_mult_sel, recMass, genPhi, recPhi, recRap, recEta};
    int counter = 0;
    vector<string> xaxis_names = {"Mass (GeV/c^{2})", "Multiplicity (%)", "Multiplicity (%)", "Mass (GeV/c^{2})", "#phi", "#phi", "Rapidity", "#eta"};

    for (auto hist : hists)
    {
        SetHistoQA(hist);
        hist->SetTitle(histnames[counter].c_str());
        hist->GetYaxis()->SetTitle("Counts");
        hist->GetXaxis()->SetTitle(xaxis_names[counter].c_str());
        TCanvas *c = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c, 0.15, 0.05, 0.05, 0.15);
        hist->SetMaximum(hist->GetMaximum() * 1.4);
        hist->Draw();

        std::string filename = Form("plots/%s.png", histnames[counter].c_str());
        // c->SaveAs(filename.c_str());
        counter++;
        c->Close();
    }

    int lowbin_mult = GenpT->GetAxis(0)->FindBin(0.0);
    int highbin_mult = GenpT->GetAxis(0)->FindBin(100.0);
    GenpT->GetAxis(0)->SetRange(lowbin_mult, highbin_mult);
    TH1D *GenpT_proj1 = GenpT->Projection(1);
    SetHistoQA(GenpT_proj1);

    TCanvas *c2 = new TCanvas("", "Generatee p_{T}", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.05, 0.05, 0.15);
    GenpT_proj1->GetYaxis()->SetTitle("Counts");
    GenpT_proj1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    GenpT_proj1->Draw();

    // lowbin_mult = GenpT->GetAxis(0)->FindBin(0.0);
    // highbin_mult = GenpT->GetAxis(0)->FindBin(80.0);
    // GenpT->GetAxis(0)->SetRange(lowbin_mult, highbin_mult);
    // TH1D *GenpT_proj2 = GenpT->Projection(1);
    // SetHistoQA(GenpT_proj2);
    // GenpT_proj2->SetLineColor(kRed);
    // // GenpT_proj2->Draw("same");

    // lowbin_mult = GenpT->GetAxis(0)->FindBin(20.0);
    // highbin_mult = GenpT->GetAxis(0)->FindBin(80.0);
    // GenpT->GetAxis(0)->SetRange(lowbin_mult, highbin_mult);
    // TH1D *GenpT_proj3 = GenpT->Projection(1);
    // SetHistoQA(GenpT_proj3);
    // GenpT_proj3->SetLineColor(kBlue);
    // // GenpT_proj3->Draw("same");

    // TLegend *ltemp = new TLegend(0.55, 0.65, 0.85, 0.90);
    // SetLegendStyle(ltemp);
    // ltemp->SetHeader("Multiplicity percentile");
    // ltemp->AddEntry(GenpT_proj1, "0-100%", "l");
    // ltemp->AddEntry(GenpT_proj2, "0-80%", "l");
    // ltemp->AddEntry(GenpT_proj3, "20-80%", "l");
    // // ltemp->Draw("same");

    // c2->SaveAs("plots/GenpT_comp.png");

    TCanvas *cmult = new TCanvas("", "Gen multiplicity", 720, 720);
    SetCanvasStyle(cmult, 0.15, 0.05, 0.05, 0.15);
    GenpT->GetAxis(0)->SetRange(lowbin_mult, highbin_mult);
    TH1D *genmult = GenpT->Projection(0);
    SetHistoQA(genmult);
    genmult->GetYaxis()->SetTitle("Counts");
    genmult->GetXaxis()->SetTitle("Multiplicity (%)");
    genmult->Draw();
    // cmult->SaveAs("plots/Gen_mult.png");

    // for reconstructed pt
    lowbin_mult = recpt1->GetAxis(0)->FindBin(0.0);
    highbin_mult = recpt1->GetAxis(0)->FindBin(100.0);
    recpt1->GetAxis(0)->SetRange(lowbin_mult, highbin_mult);
    TH1D *recpt1_proj1 = recpt1->Projection(1);
    SetHistoQA(recpt1_proj1);
    TCanvas *c3 = new TCanvas("", "Reconstructed p_{T}", 720, 720);
    SetCanvasStyle(c3, 0.15, 0.05, 0.05, 0.15);
    recpt1_proj1->GetYaxis()->SetTitle("Counts");
    recpt1_proj1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    recpt1_proj1->Draw();

    lowbin_mult = recpt1->GetAxis(0)->FindBin(0.0);
    highbin_mult = recpt1->GetAxis(0)->FindBin(80.0);
    recpt1->GetAxis(0)->SetRange(lowbin_mult, highbin_mult);
    TH1D *recpt1_proj2 = recpt1->Projection(1);
    SetHistoQA(recpt1_proj2);
    recpt1_proj2->SetLineColor(kRed);
    // recpt1_proj2->Draw("same");

    lowbin_mult = recpt1->GetAxis(0)->FindBin(20.0);
    highbin_mult = recpt1->GetAxis(0)->FindBin(80.0);
    recpt1->GetAxis(0)->SetRange(lowbin_mult, highbin_mult);
    TH1D *recpt1_proj3 = recpt1->Projection(1);
    SetHistoQA(recpt1_proj3);
    recpt1_proj3->SetLineColor(kBlue);
    // recpt1_proj3->Draw("same");

    TLegend *ltemp1 = new TLegend(0.55, 0.70, 0.85, 0.92);
    SetLegendStyle(ltemp1);
    ltemp1->SetHeader("Multiplicity percentile");
    ltemp1->AddEntry(recpt1_proj1, "0-100%", "l");
    ltemp1->AddEntry(recpt1_proj2, "0-80%", "l");
    ltemp1->AddEntry(recpt1_proj3, "20-80%", "l");
    // ltemp1->Draw("same");

    // c3->SaveAs("plots/RecpT_comp.png");

    TCanvas *c4 = new TCanvas("", "GenpT", 720, 720);
    SetCanvasStyle(c4, 0.15, 0.05, 0.05, 0.15);
    TH1F *efficiency = (TH1F *)recpt1_proj1->Clone();
    efficiency->Divide(GenpT_proj1);
    SetHistoQA(efficiency);
    efficiency->SetStats(0);
    efficiency->GetYaxis()->SetTitle("Acceptance #times Efficiency");
    efficiency->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    efficiency->Draw();
    // c4->SaveAs("plots/efficiency.png");

    // TCanvas *c5 = new TCanvas("", "Gen pT", 720, 720);
    // SetCanvasStyle(c5, 0.15, 0.05, 0.05, 0.15);
    // genRap->GetAxis(1)->SetRange(genRap->GetAxis(1)->FindBin(-1.0), genRap->GetAxis(1)->FindBin(1.0)); // pt axis
    // TH1D *genRap_proj0 = genRap->Projection(0);
    // SetHistoQA(genRap_proj0);
    // genRap_proj0->GetYaxis()->SetTitle("Counts");
    // genRap_proj0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // genRap_proj0->GetXaxis()->SetRangeUser(0, 10);
    // genRap_proj0->Draw();
    // c5->SaveAs("plots/GenpT.png");
    // cout<<"pt bin width: "<<genRap->GetAxis(1)->GetBinWidth(1)<<endl;

    // TCanvas *c6 = new TCanvas("", "Gen rapidity", 720, 720);
    // SetCanvasStyle(c6, 0.15, 0.05, 0.05, 0.15);
    // genRap->GetAxis(0)->SetRange(genRap->GetAxis(0)->FindBin(-0.1), genRap->GetAxis(0)->FindBin(0.1)); // rapidity axis
    // TH1D *genRap_proj1 = genRap->Projection(1);
    // SetHistoQA(genRap_proj1);
    // genRap_proj1->GetYaxis()->SetTitle("Counts");
    // genRap_proj1->GetXaxis()->SetTitle("Rapidity");
    // // genRap_proj1->GetXaxis()->SetRangeUser(0, 10);
    // genRap_proj1->Draw();
    // c6->SaveAs("plots/GenRapidity.png");

    TCanvas *c7 = new TCanvas("", "Gen pt vs Eta 2D", 720, 720);
    SetCanvasStyle(c7, 0.15, 0.05, 0.05, 0.15);
    TH2F *h2pteta = (TH2F *)genRap->Projection(1, 0);
    Set2Dstyle(h2pteta);
    h2pteta->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    h2pteta->GetYaxis()->SetTitle("y");
    h2pteta->Draw("colz");
}