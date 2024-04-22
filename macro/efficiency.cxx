#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

using namespace std;

void efficiency()
{
    // TFile *fileeff = new TFile("/home/sawan/check_k892/data/pp/kstar/199692.root", "READ");
    // TFile *fileeff = new TFile("/home/sawan/check_k892/200962.root", "READ");
    // TFile *fileeff = new TFile("/home/sawan/check_k892/phi_mc.root", "READ");
    TFile *fileeff = new TFile("/home/sawan/check_k892/mc/199692.root", "READ");
    // TFile *fileeff = new TFile("/home/sawan/check_k892/data/AnalysisResults.root", "READ");

    TFile *fileraw = new TFile((kSignalOutput + "/yield.root").c_str(), "READ");
    if (fileeff->IsZombie() || fileraw->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }
    // TH2F *h2gen = (TH2F *)fileeff->Get("lf-k892analysis/k892Gen");
    // TH2F *h2genanti = (TH2F*)fileeff->Get("lf-k892analysis/k892GenAnti");
    // TH2F *h2rec = (TH2F *)fileeff->Get("lf-k892analysis/k892Rec");
    // TH2F *h2recanti = (TH2F*)fileeff->Get("lf-k892analysis/k892RecAnti");
    // if(h2gen == nullptr || h2genanti == nullptr || h2rec == nullptr || h2recanti == nullptr) {
    //     cout << "Error reading histograms" << endl;
    //     return;
    // }

    // TH1D *h1gen = h2gen->ProjectionX("h1gen");
    // TH1D *h1genanti = h2genanti->ProjectionX("h1genanti");
    // TH1D *h1rec = h2rec->ProjectionX("h1rec");
    // TH1D *h1recanti = h2recanti->ProjectionX("h1recanti");

    TH1F *h1gen = (TH1F *)fileeff->Get("kstarqa_all_event_sel/histos/k892Gen");
    TH1F *h1rec = (TH1F *)fileeff->Get("kstarqa_all_event_sel/histos/h1recpt");
    //   TH1F *h1gen = (TH1F *)fileeff->Get("phianalysisrun3_id10133/h1PhiGen");
    // TH1F *h1rec = (TH1F *)fileeff->Get("phianalysisrun3_id10133/h1PhiRecsplit");
    // if (h1gen == nullptr || h1genanti == nullptr || h1rec == nullptr || h1recanti == nullptr)
    // {
    //     cout << "Error reading histograms" << endl;
    //     return;
    // }
    if(h1gen == nullptr || h1rec == nullptr) {
        cout << "Error reading efficiency histograms" << endl;
        return;
    }
    // h1gen->Add(h1genanti, 1);
    // h1rec->Add(h1recanti, 1);
    TH1F *hyieldraw = (TH1F *)fileraw->Get("hintegral_yield");
    if (hyieldraw == nullptr)
    {
        cout << "Error reading yield histogram" << endl;
        return;
    }
    TH1F *heff = (TH1F *)hyieldraw->Clone();
    TH1F *heff2 = (TH1F *)h1gen->Clone();
    cout<<"bins: "<<heff->GetNbinsX()<<endl;
    for (int i = 0; i < heff->GetNbinsX(); i++)
    // for (int i = 0; i < heff2->GetNbinsX(); i++)
    {
        lowpt = pT_bins[i];
        highpt = pT_bins[i + 1];
        cout<<"lowpt: "<<lowpt<<" highpt: "<<highpt<<endl;
        double efficiency = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt)) / h1gen->Integral(h1gen->GetXaxis()->FindBin(lowpt), h1gen->GetXaxis()->FindBin(highpt));
        // double efficiency2 = h1rec->Integral(i+1, i+2)/h1gen->Integral(i+1, i+2);
        cout << "Efficiency: " << efficiency << endl;
        // cout << "Efficiency2: " << efficiency2 << "\n\n"<< endl;
        heff->SetBinContent(i + 1, efficiency);
        // heff2->SetBinContent(i + 1, efficiency2);
        hyieldraw->SetBinContent(i + 1, hyieldraw->GetBinContent(i + 1) / efficiency);
    }
    TCanvas *c1 = new TCanvas();
    // heff->Scale(55.0/80.0);
    heff->Draw();
    // heff2->Draw();

    TCanvas *c2 = new TCanvas();
    gPad->SetLogy();
    hyieldraw->Draw();

    TFile *spectra = new TFile((kSignalOutput + "/spectra.root").c_str(), "RECREATE");
    heff->Write("heff");
    hyieldraw->Write("spectra");
}