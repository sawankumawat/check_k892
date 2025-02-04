#include<iostream>
#include "../../macro/src/style.h"
using namespace std;

void check_mc() {
    TFile *f = new TFile("316648.root");
    TH2D *h2gen = (TH2D*)f->Get("kstarqa/hInvMass/hk892GenpT");
    TH2D *h2rec = (TH2D*)f->Get("kstarqa/hInvMass/h2KstarRecpt2");
    TH1F *hsplit = (TH1F*)f->Get("kstarqa/hInvMass/h1KSRecsplit");
    if(h2gen == nullptr || h2rec == nullptr || hsplit == nullptr) {
        cout << "Histograms not found" << endl;
        return;
    }

    TH1D *hgenpt = h2gen->ProjectionX();
    TH1D *hrecpt = h2rec->ProjectionX();

    TCanvas *c1 = new TCanvas ("", "", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.04, 0.06, 0.15);
    SetHistoQA(hgenpt);
    hgenpt->Draw();

    TCanvas *c2 = new TCanvas ("", "", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.04, 0.06, 0.15);
    SetHistoQA(hrecpt);
    hrecpt->Draw();

    cout<<"xlow in gen and rec and split: "<<hgenpt->GetXaxis()->GetBinLowEdge(1)<<" "<<hrecpt->GetXaxis()->GetBinLowEdge(1)<< " " <<hsplit->GetXaxis()->GetBinLowEdge(1)<<endl;
    cout<<"xup in gen and rec and split: "<<hgenpt->GetXaxis()->GetBinUpEdge(hgenpt->GetNbinsX())<<" "<<hrecpt->GetXaxis()->GetBinUpEdge(hrecpt->GetNbinsX())<<" "<<hsplit->GetXaxis()->GetBinUpEdge(hsplit->GetNbinsX())<<endl;
    cout<<"no of bins generated: "<<hgenpt->GetNbinsX()<<endl;
    cout<<"no of bins reconstructed: "<<hrecpt->GetNbinsX()<<endl;
    cout<<"no of bins split: "<<hsplit->GetNbinsX()<<endl;

    TCanvas *c3 = new TCanvas ("", "", 720, 720);
    SetCanvasStyle(c3, 0.15, 0.04, 0.06, 0.15);
    TH1F *heff = (TH1F *)hrecpt->Clone();
    heff->GetXaxis()->SetRangeUser(0, 10);
    hgenpt->GetXaxis()->SetRangeUser(0, 10);
    // heff->Divide(hgenpt);
    heff->Divide(hsplit);
    SetHistoQA(heff);
    heff->Rebin(5);
    heff->Draw();

}