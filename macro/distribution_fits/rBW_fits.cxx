#include <iostream>
#include "TDatabasePDG.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"

void rBW_fits()
{
    //constant parameters ****************************************
    const string kResBkg = "MIX";
    const string outputfolder_str = "../" + kSignalOutput + "/" + kchannel + "/" + kfoldername;
    gStyle->SetOptStat(1110);

    //pT loop ***************************************************
    for (Int_t ip = pt_start; ip < pt_end; ip++)
    {
        TFile *f = new TFile((outputfolder_str + "/hglue_KsKs_" + kResBkg + Form("_%f_%f.root", pT_bins[ip], pT_bins[ip + 1])).c_str(), "READ");

        if (f->IsZombie())
        {
            cout << "Error opening file" << endl;
            return;
        }

        TH1F *hinvMass = (TH1F *)f->Get("ksks_invmass");
        if(hinvMass == nullptr)
        {
            cout << "Error opening histogram" << endl;
            return;
        }

        TCanvas *c1 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c1, 0.12, 0.03, 0.05, 0.14);
        hinvMass->GetYaxis()->SetTitleOffset(1.1);
        hinvMass->GetXaxis()->SetTitleOffset(1.3);
        hinvMass->Draw();
    }
}