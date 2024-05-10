#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

using namespace std;

void efficiency()
{
    TString outputfolder = kSignalOutput + "/efficiency";
    gSystem->mkdir(outputfolder, kTRUE);
    TFile *fileeff = new TFile(kMCFilename.c_str(), "READ");
    TFile *fileraw = new TFile((kSignalOutput + "/yield.root").c_str(), "READ");
    if (fileeff->IsZombie() || fileraw->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }

    TH1F *h1gen = (TH1F *)fileeff->Get("kstarqa/histos/k892Gen");
    // TH1F *h1rec = (TH1F *)fileeff->Get("kstarqa/histos/h3KstarRec");
    TH1F *h1rec = (TH1F *)fileeff->Get("kstarqa/histos/h1recpt");
   
    if (h1gen == nullptr || h1rec == nullptr)
    {
        cout << "Error reading efficiency histograms" << endl;
        return;
    }
    // h1gen->Add(h1genanti, 1);
    // h1rec->Add(h1recanti, 1);
    TH1F *hyieldraw = (TH1F *)fileraw->Get("yield_integral");
    // TH1F *hyieldraw = (TH1F *)fileraw->Get("yield_bincount");
    if (hyieldraw == nullptr)
    {
        cout << "Error reading yield histogram" << endl;
        return;
    }
    TH1F *heff = (TH1F *)hyieldraw->Clone();
    cout << "bins: " << heff->GetNbinsX() << endl;
    for (int i = 0; i < heff->GetNbinsX(); i++)
    {
        lowpt = pT_bins[i];
        highpt = pT_bins[i + 1];
        cout << "lowpt: " << lowpt << " highpt: " << highpt << endl;
        double nrec = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt));
        double ngen = h1gen->Integral(h1gen->GetXaxis()->FindBin(lowpt), h1gen->GetXaxis()->FindBin(highpt));
        double efficiency = nrec / ngen;
        double efficiencyerr = sqrt(((nrec + 1) / (ngen + 2)) * ((nrec + 2) / (ngen + 3) - (nrec + 1) / (ngen + 2)));
        cout << "Efficiency: " << efficiency << " +/- " << efficiencyerr << endl;
        heff->SetBinContent(i + 1, efficiency);
        // heff->SetBinError(i + 1, efficiencyerr);
        hyieldraw->SetBinContent(i + 1, hyieldraw->GetBinContent(i + 1) / efficiency);
    }
    TCanvas *c1 = new TCanvas();
    heff->Draw();
    c1->SaveAs(outputfolder + "/efficiency.png");

    TCanvas *c2 = new TCanvas();
    gPad->SetLogy();
    hyieldraw->Draw();
    c2->SaveAs(outputfolder + "/corrected_spectra.png");

    TFile *spectra = new TFile((outputfolder + "/corrected_spectra.root"), "RECREATE");
    heff->Write("heff");
    hyieldraw->Write("spectra");

    TCanvas *c3 = new TCanvas();
    TH1F *heffclone = (TH1F *)h1rec->Clone();
    heffclone->Divide(h1gen);
    heffclone->Draw();

    TCanvas *c4 = new TCanvas();
    h1gen->Draw();
    h1rec->SetLineColor(kRed);
    h1rec->Draw("same");
}