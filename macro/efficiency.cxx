#include<iostream>
#include "src/style.h"
#include "src/initializations.h"

using namespace std;

void efficiency() {
    TFile *fileeff = new TFile("/home/sawan/check_k892/mc/192370.root", "READ");
    TFile *fileraw = new TFile((kSignalOutput + "/yield.root").c_str(), "READ");
    if(fileeff->IsZombie() || fileraw->IsZombie()) {
        cout << "Error opening files" << endl;
        return;
    }
    TH2F *h2gen = (TH2F*)fileeff->Get("lf-k892analysis/k892Gen");
    TH2F *h2genanti = (TH2F*)fileeff->Get("lf-k892analysis/k892GenAnti");
    TH2F *h2rec = (TH2F*)fileeff->Get("lf-k892analysis/k892Rec");
    TH2F *h2recanti = (TH2F*)fileeff->Get("lf-k892analysis/k892RecAnti");
    if(h2gen == nullptr || h2genanti == nullptr || h2rec == nullptr || h2recanti == nullptr) {
        cout << "Error reading histograms" << endl;
        return;
    }
    TH1D *h1gen = h2gen->ProjectionX("h1gen");
    TH1D *h1genanti = h2genanti->ProjectionX("h1genanti");
    TH1D *h1rec = h2rec->ProjectionX("h1rec");
    TH1D *h1recanti = h2recanti->ProjectionX("h1recanti");
    if(h1gen == nullptr || h1genanti == nullptr || h1rec == nullptr || h1recanti == nullptr) {
        cout << "Error reading histograms" << endl;
        return;
    }
    h1gen->Add(h1genanti, 1);
    h1rec->Add(h1recanti, 1);
    TH1F *hyieldraw = (TH1F*)fileraw->Get("hintegral_yield");
    if(hyieldraw == nullptr) {
        cout << "Error reading histogram" << endl;
        return;
    }
    TH1F *heff = (TH1F*)hyieldraw->Clone();
    
    for (int i = 0; i < heff->GetNbinsX(); i++)
    {
        lowpt = pT_bins[i];
        highpt = pT_bins[i + 1];
        double efficiency = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt)) / h1gen->Integral(h1gen->GetXaxis()->FindBin(lowpt), h1gen->GetXaxis()->FindBin(highpt));
        cout<<"Efficiency: "<<efficiency<<endl;
        heff->SetBinContent(i + 1, efficiency);
        hyieldraw->SetBinContent(i + 1, hyieldraw->GetBinContent(i + 1) / efficiency);
    }

    // heff->Draw();
    hyieldraw->Draw();

    
}