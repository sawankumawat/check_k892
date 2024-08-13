#include <iostream>
using namespace std;
#include "src/style.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "TMath.h"

void rootfit_mult()
{
    gStyle->SetOptStat(0);
    TFile *f = new TFile("../data/glueball/LHC22o_pass7_small/248760.root", "READ");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    TCanvas *cmultdist = new TCanvas("cmultdist", "cmultdist", 720, 720);
    gPad->SetLogy();
    string path_mult = "strangeness_tutorial/eventSelection/multdist_FT0M";
    TH1F *h_mult = (TH1F *)f->Get(path_mult.c_str());
    if (h_mult == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }
    SetCanvasStyle(cmultdist, 0.15, 0.03, 0.03, 0.14);
    SetHistoQA(h_mult);
    h_mult->GetXaxis()->SetTitle("Multiplicity");
    h_mult->GetYaxis()->SetTitle("Counts");
    h_mult->Draw("hist");

    // Convert TH1F to RooDataHist
    RooRealVar x("x", "Multiplicity", 0, 60000);
    RooDataHist data("data", "dataset with x", RooArgList(x), h_mult);

    // Define the parameters for the Negative Binomial Distribution (NBD)
    RooRealVar par0("par0", "Normalization", 1e7, 0, 1e10);
    RooRealVar par1("par1", "Number of Successes", 5, 0, 50);
    RooRealVar par2("par2", "Success Probability", 0.3, 0.0, 1.0);

    // Define the NBD using RooGenericPdf
    RooGenericPdf nbd("nbd", "par0*TMath::Gamma(x+par1)/(TMath::Gamma(x+1)*TMath::Gamma(par1))*TMath::Power(par2, par1)*TMath::Power(1-par2, x)",
                      RooArgSet(x, par0, par1, par2));

    // Fit the model to the data
    RooFitResult* result = nbd.fitTo(data, RooFit::Save());

    // Plot the results
    RooPlot* frame = x.frame();
    data.plotOn(frame);
    nbd.plotOn(frame);
    frame->Draw();

    // Optionally, print the fit results
    result->Print();
}
