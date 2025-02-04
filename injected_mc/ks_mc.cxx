#include <iostream>
#include "style.h"
using namespace std;

void ks_mc()
{
    gStyle->SetOptStat(0);
    TFile *f = new TFile("Ks_mc.root", "read");
    TH2D *hgen = (TH2D *)f->Get("lf-v0qaanalysis/Generated_MCGenRecoColl_INELgt0_K0Short");
    TH2D *hrec = (TH2D *)f->Get("lf-v0postprocessing/hMassVsPtK0Short");
    if (hgen == nullptr || hrec == nullptr)
    {
        cout << "Histograms not found" << endl;
        return;
    }
    TH1D *h1gen = hgen->ProjectionX();
    TH1D *h1rec = hrec->ProjectionX();

    // TCanvas *c1 = new TCanvas("", "", 720, 720);
    // SetCanvasStyle(c1, 0.15, 0.04, 0.06, 0.15);
    // SetHistoQA(h1gen);
    // h1gen->Draw();

    // TCanvas *c2 = new TCanvas("", "", 720, 720);
    // SetCanvasStyle(c2, 0.15, 0.04, 0.06, 0.15);
    // SetHistoQA(h1rec);
    // h1rec->Draw();

    TCanvas *c3 = new TCanvas("c3", "acc x eff for K^{0}_{s}", 720, 720);
    SetCanvasStyle(c3, 0.15, 0.04, 0.06, 0.15);
    TH1F *heff = (TH1F *)h1rec->Clone();
    heff->Divide(h1gen);
    SetHistoQA(heff);
    // heff->Rebin(2);
    heff->GetYaxis()->SetTitle(Form("(Acceptance #times Efficiency / %.1f GeV/#it{c})", heff->GetBinWidth(1)));
    heff->GetYaxis()->SetTitleOffset(1.5);
    heff->Draw();

    // Now lets calcualate acc x eff for the higher mass resonances
    TFile *f2 = new TFile("/home/sawan/check_k892/mc/LHC24l1/316063_onlypass7.root", "read");
    if (f2->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *GenpT = (TH1F *)f2->Get("higher-mass-resonances/hMChists/Genf1710");
    TH1F *recpt1 = (TH1F *)f2->Get("higher-mass-resonances/hMChists/Recf1710_pt1");
    // TH1F *recpt2 = (TH1F *)f2->Get("higher-mass-resonances/hMChists/Recf1710_pt2");

    TH1F *heff_reso = (TH1F *)recpt1->Clone();
    heff_reso->Divide(GenpT);
    TCanvas *c4 = new TCanvas("c4", "acc x eff for f_{0}(1710)", 720, 720);
    SetCanvasStyle(c4, 0.16, 0.04, 0.06, 0.15);
    SetHistoQA(heff_reso);
    // heff_reso->Rebin(2);
    heff_reso->GetYaxis()->SetTitle(Form("(Acceptance #times Efficiency / %.1f GeV/#it{c})", heff_reso->GetBinWidth(1)));
    heff_reso->GetYaxis()->SetTitleOffset(1.5);
    heff_reso->Draw();

    cout << "bin width for Ks: " << heff->GetBinWidth(1) << endl;
    cout << "bin width for f0(1710): " << heff_reso->GetBinWidth(1) << endl;

    TCanvas *ccompare = new TCanvas("ccompare", "Comparison of acc x eff for K^{0}_{s} and f_{0}(1710)", 720, 720);
    SetCanvasStyle(ccompare, 0.16, 0.04, 0.06, 0.15);
    heff->SetLineColor(2);
    heff_reso->SetLineColor(4);
    // lets multiply heff by itself
    for (int i = 1; i <= heff->GetNbinsX(); i++)
    {
        heff->SetBinContent(i, heff->GetBinContent(i) * heff->GetBinContent(i));
        // heff->SetBinContent(i, heff->GetBinContent(i));
    }
    heff->Draw();
    heff_reso->Draw("same");

    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    SetLegendStyle(leg);
    leg->SetHeader("Acc #times Eff");
    leg->AddEntry(heff, "K^{0}_{s} #times K^{0}_{s}", "l");
    leg->AddEntry(heff_reso, "f_{0}(1710)", "l");
    leg->Draw();

    ccompare->SaveAs("plots/acc_x_eff_compare.png");
}