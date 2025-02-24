#include <iostream>
using namespace std;
#include "style.h"

void analyze_lambda()
{
    TFile *f = new TFile("AnalysisResults_Lstar_MC.root", "read");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH3D *hgen3d = (TH3D *)f->Get("lf-lambda1520analysis_PID6/Result/MC/Genlambda1520pt");
    TH2D *hrec2d = (TH2D *)f->Get("lf-lambda1520analysis_PID6/Result/MC/lambda1520Reco");
    if (hgen3d == nullptr || hrec2d == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    TH1D *hgenpt = hgen3d->ProjectionY("hgenpt1", hgen3d->GetZaxis()->FindBin(0.0), hgen3d->GetZaxis()->FindBin(100.0), hgen3d->GetXaxis()->FindBin(2.0 - 0.1), hgen3d->GetXaxis()->FindBin(2.0 + 0.1), "E");

    TCanvas *cgen = new TCanvas("", "Generated p_{T}", 720, 720);
    SetCanvasStyle(cgen, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hgenpt);
    hgenpt->Draw();
    cgen->SaveAs("plots/GenpT_lambda.png");

    TH1D *hrecpt = hrec2d->ProjectionX();
    TCanvas *crec = new TCanvas("", "Reconstructed p_{T}", 720, 720);
    SetCanvasStyle(crec, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hrecpt);
    hrecpt->Draw("same");
    crec->SaveAs("plots/RecpT_lambda.png");

    cout<<"bin width gen : "<<hgenpt->GetBinWidth(1)<<endl;
    cout<<"bin width rec : "<<hrecpt->GetBinWidth(1)<<endl;

    TCanvas *ceff = new TCanvas("", "Efficiency", 720, 720);
    SetCanvasStyle(ceff, 0.15, 0.05, 0.05, 0.15);
    TH1D *heff = (TH1D *)hrecpt->Clone("heff");
    heff->Divide(hrecpt, hgenpt);
    SetHistoQA(heff);
    heff->Draw();
    ceff->SaveAs("plots/Efficiency_lambda.png");

}