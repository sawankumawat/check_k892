#include <iostream>
#include "../../macro/src/style.h"
using namespace std;

void ks_mc()
{
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

    TCanvas *c1 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.04, 0.06, 0.15);
    SetHistoQA(h1gen);
    h1gen->Draw();

    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.04, 0.06, 0.15);
    SetHistoQA(h1rec);
    h1rec->Draw();

    TCanvas *c3 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c3, 0.15, 0.04, 0.06, 0.15);
    TH1F *heff = (TH1F *)h1rec->Clone();
    heff->Divide(h1gen);
    SetHistoQA(heff);
    // heff->Rebin(5);
    heff->Draw();
}