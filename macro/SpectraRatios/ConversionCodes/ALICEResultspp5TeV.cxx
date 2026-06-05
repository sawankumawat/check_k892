#include <iostream>
#include <vector>
#include "../../src/style.h"

using namespace std;

void ALICEResultspp5TeV()
{
    TFile *fALICE = new TFile("../HEP_data/pp5.02TeV_INELgt0.root", "read");

    if (!fALICE || fALICE->IsZombie())
    {
        cout << "Error: ALICE file not found" << endl;
        return;
    }

    vector<pair<int, TString>> tables = {
        {5, "Kstar_MeanYield"},
        {6, "Kstar_MeanpT"},
        {8, "Kstar_KaRatio"}};

    TGraphErrors *gKstar_MeanYield = (TGraphErrors *)fALICE->Get("Table 5/Graph1D_y1");
    TGraphErrors *gKstar_MeanpT = (TGraphErrors *)fALICE->Get("Table 6/Graph1D_y1");
    TGraphErrors *gKstar_KaRatio = (TGraphErrors *)fALICE->Get("Table 8/Graph1D_y1");

    if (!gKstar_MeanYield || !gKstar_MeanpT || !gKstar_KaRatio)
    {
        cout << "Error: One or more graphs not found in ALICE file" << endl;
        return;
    }

    TFile *fOutput = new TFile("pp5p02TeVALICE.root", "RECREATE");
    gKstar_MeanYield->Write("gKstar_MeanYield");
    gKstar_MeanpT->Write("gKstar_MeanpT");
    gKstar_KaRatio->Write("gKstar_KaRatio");

    fOutput->Close();
    fALICE->Close();

    cout << "Output written to pp5p02TeVALICE.root" << endl;
}