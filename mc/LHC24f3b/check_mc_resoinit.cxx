#include <iostream>
#include "../../macro/src/style.h"
using namespace std;

void check_mc_resoinit()
{
    TFile *f = new TFile("279218.root", "read");
    TH2D *hrec = (TH2D *)f->Get("lf-k892analysis/k892Rec");
    THnSparseD *hgen = (THnSparseD *)f->Get("lf-k892analysis/k892Gen");
    if (hrec == nullptr || hgen == nullptr)
    {
        cout << "Histograms not found" << endl;
        return;
    }
    gStyle->SetOptStat(0);

    TH1D *h1rec = hrec->ProjectionX();
    int lowbinx = hgen->GetAxis(0)->FindBin(2 - 0.001);
    int highbinx = hgen->GetAxis(0)->FindBin(2 + 0.001);
    hgen->GetAxis(0)->SetRange(lowbinx, highbinx);

    TH1D *h1gen = (TH1D *)hgen->Projection(1);

    cout << "xlow in gen and rec: " << h1gen->GetXaxis()->GetBinLowEdge(1) << " " << h1rec->GetXaxis()->GetBinLowEdge(1) << endl;
    cout << "xup in gen and rec: " << h1gen->GetXaxis()->GetBinUpEdge(h1gen->GetNbinsX()) << " " << h1rec->GetXaxis()->GetBinUpEdge(h1rec->GetNbinsX()) << endl;
    cout << "no of bins generated: " << h1gen->GetNbinsX() << endl;
    cout << "no of bins reconstructed: " << h1rec->GetNbinsX() << endl;

    TCanvas *c1 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.04, 0.06, 0.15);
    SetHistoQA(h1rec);
    h1rec->GetYaxis()->SetTitle("Counts");
    h1rec->Draw();
    c1->SaveAs("plots/h1rec.png");

    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.04, 0.06, 0.15);
    SetHistoQA(h1gen);
    h1gen->GetYaxis()->SetTitle("Counts");
    h1gen->Draw();
    c2->SaveAs("plots/h1gen.png");

    TCanvas *c3 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c3, 0.15, 0.04, 0.06, 0.15);
    TH1F *heff = (TH1F *)h1rec->Clone();
    heff->Divide(h1gen);
    SetHistoQA(heff);
    heff->Rebin(2);
    heff->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    heff->Draw();
    c3->SaveAs("plots/heff.png");
}