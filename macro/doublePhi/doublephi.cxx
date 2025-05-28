#include <iostream>
#include "../src/style.h"
using namespace std;

void doublephi()
{
    // LHC24_pass1 data
    // TFile *file = new TFile("../../data/doublePhi/387839.root"); // PID 0
    // TFile *file = new TFile("../../data/doublePhi/388086.root"); // PID 1
    // TFile *file = new TFile("../../data/doublePhi/388535.root"); // PID 2
    // TFile *file = new TFile("../../data/doublePhi/388352.root"); // PID 3

    // LHC23_pass4 data
    // TFile *file = new TFile("../../data/doublePhi/389481.root"); // PID 0
    TFile *file = new TFile("../../data/doublePhi/merge.root"); // merged 23 and 24 dataset output for PID 0

    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << "\n";
        return;
    }

    TIter nextkey(file->GetListOfKeys());
    TKey *key;
    // print all keys in the file
    while ((key = (TKey *)nextkey()))
    {
        std::cout << key->GetName() << "\n";
    }

    // now get keys inside the first key
    TKey *key1 = (TKey *)file->GetListOfKeys()->At(0);
    TDirectory *dir = (TDirectory *)key1->ReadObj();
    TIter nextkey1(dir->GetListOfKeys());
    TKey *key2;
    // print all keys in the directory
    while ((key2 = (TKey *)nextkey1()))
    {
        std::cout << key2->GetName() << "\n";
    }

    THnSparseF *hunlike = (THnSparseF *)file->Get("doublephimeson/SEMassUnlike");
    THnSparseF *hmixed = (THnSparseF *)file->Get("doublephimeson/MEMassUnlike");
    if (hunlike == nullptr || hmixed == nullptr)
    {
        std::cerr << "Error: Could not find histogram 'doublephimeson/SEMassUnlike' or 'doublephimeson/MEMassUnlike' in file\n";
        return;
    }

    int lowbinpT = hunlike->GetAxis(1)->FindBin(2.0 + 0.001);
    int highbinpT = hunlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int deltaRlow = hunlike->GetAxis(4)->FindBin(0.0 + 0.001);
    int deltaRhigh = hunlike->GetAxis(4)->FindBin(2.0 - 0.001);

    hunlike->GetAxis(1)->SetRange(lowbinpT, highbinpT);
    hunlike->GetAxis(4)->SetRange(deltaRlow, deltaRhigh);
    hmixed->GetAxis(1)->SetRange(lowbinpT, highbinpT);
    hmixed->GetAxis(4)->SetRange(deltaRlow, deltaRhigh);

    TH1D *hmass = hunlike->Projection(0, "E");
    if (hmass == nullptr)
    {
        std::cerr << "Error: Could not create projection histogram\n";
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.17, 0.03, 0.05, 0.15);
    SetHistoQA(hmass);
    hmass->GetXaxis()->SetTitle("#it{M}_{#Phi#Phi} (GeV/#it{c}^{2})");
    hmass->SetMarkerStyle(20);
    hmass->GetYaxis()->SetMaxDigits(3);
    hmass->SetMarkerSize(1.0);

    TH1D *hmassmixed = hmixed->Projection(0, "E");
    SetHistoQA(hmassmixed);
    int normlow = hmassmixed->GetXaxis()->FindBin(2.9 + 0.001);
    int normhigh = hmassmixed->GetXaxis()->FindBin(3.0 - 0.001);
    auto signal_integral = hmass->Integral(normlow, normhigh);
    auto bkg_integral = hmassmixed->Integral(normlow, normhigh);
    auto normfactor = signal_integral / bkg_integral;
    hmassmixed->Scale(normfactor);

    hmassmixed->Rebin(4);
    hmass->Rebin(4);
    hmass->GetXaxis()->SetRangeUser(2.68, 2.84);
    hmass->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hmass->GetBinWidth(1) * 1000));
    hmass->GetYaxis()->SetTitleOffset(1.8);
    hmass->Draw();

    hmassmixed->SetLineColor(kRed);
    hmassmixed->SetMarkerColor(kRed);
    hmassmixed->Draw("same");

    TLegend *leg = new TLegend(0.25, 0.18, 0.5, 0.3);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry(hmass, "Same event #Phi#Phi pair", "lpe");
    leg->AddEntry(hmassmixed, "Mixed event #Phi#Phi pair", "lpe");
    leg->Draw();
    c1->SaveAs("raw_mass.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.03, 0.05, 0.15);
    TH1D *hmassclone = (TH1D *)hmass->Clone();
    hmassclone->Add(hmassmixed, -1);
    SetHistoQA(hmassclone);
    hmassclone->GetXaxis()->SetTitle("#it{M}_{#Phi#Phi} (GeV/#it{c}^{2})");
    hmassclone->SetMarkerStyle(20);
    hmassclone->GetYaxis()->SetMaxDigits(3);
    hmassclone->SetMarkerSize(1.0);
    hmassclone->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hmass->GetBinWidth(1) * 1000));
    hmassclone->GetXaxis()->SetRangeUser(2.68, 2.84);
    hmassclone->Draw();
    c2->SaveAs("raw_mass_subtracted.png");
}