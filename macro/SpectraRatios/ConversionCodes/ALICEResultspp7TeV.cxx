#include <iostream>
#include <vector>
#include "../../src/style.h"

using namespace std;

struct DataSet
{
    int table = -1;
    TString name;

    TGraphErrors *g = nullptr;
    TH1D *hStat = nullptr;
    TH1D *hSys = nullptr;
    TH1D *hSysUncorr = nullptr;

    TGraphErrors *gStat = nullptr;
    TGraphErrors *gSys = nullptr;
    TGraphErrors *gSysUncorr = nullptr;
};

DataSet LoadTable(TFile *file, int table, const TString &name)
{
    DataSet d;

    d.table = table;
    d.name = name;

    d.g = (TGraphErrors *)file->Get(Form("Table %d/Graph1D_y1", table));
    d.hStat = (TH1D *)file->Get(Form("Table %d/Hist1D_y1_e1", table));
    d.hSys = (TH1D *)file->Get(Form("Table %d/Hist1D_y1_e2", table));
    d.hSysUncorr = (TH1D *)file->Get(Form("Table %d/Hist1D_y1_e3", table));

    // Only these are mandatory
    if (!d.g || !d.hStat || !d.hSys)
    {
        cout << "Error loading table " << table << " (" << name << ")" << endl;
        return d;
    }

    d.gStat = new TGraphErrors(d.g->GetN());
    d.gSys = new TGraphErrors(d.g->GetN());

    // Create only if e3 exists
    if (d.hSysUncorr)
        d.gSysUncorr = new TGraphErrors(d.g->GetN());

    return d;
}

void BuildErrorGraphs(DataSet &d)
{
    if (!d.g || !d.hStat || !d.hSys)
        return;

    for (int i = 0; i < d.g->GetN(); i++)
    {
        double x = 0.0;
        double y = 0.0;

        d.g->GetPoint(i, x, y);

        int bin = d.hStat->FindBin(x);

        double ex = d.g->GetErrorX(i);
        double eyStat = d.hStat->GetBinContent(bin);
        double eySys = d.hSys->GetBinContent(bin);

        d.gStat->SetPoint(i, x, y);
        d.gStat->SetPointError(i, ex, eyStat);

        d.gSys->SetPoint(i, x, y);
        d.gSys->SetPointError(i, ex, eySys);

        // Optional e3
        if (d.hSysUncorr && d.gSysUncorr)
        {
            double eySysUncorr = d.hSysUncorr->GetBinContent(bin);

            d.gSysUncorr->SetPoint(i, x, y);
            d.gSysUncorr->SetPointError(i, ex, eySysUncorr);
        }
    }
}

void WriteGraphs(TFile *outFile, const DataSet &d)
{
    if (!d.gStat || !d.gSys)
        return;

    outFile->cd();

    d.gStat->Write(Form("g%s_stat", d.name.Data()));
    d.gSys->Write(Form("g%s_sys", d.name.Data()));

    // Write only if available
    if (d.gSysUncorr)
        d.gSysUncorr->Write(Form("g%s_sysuncorr", d.name.Data()));
}

void ALICEResultspp7TeV()
{
    TFile *fALICE = new TFile("../HEP_data/HEPData_7TeV_KstarPhi_INELgt0.root", "read");

    if (!fALICE || fALICE->IsZombie())
    {
        cout << "Error: ALICE file not found" << endl;
        return;
    }

    vector<pair<int, TString>> tables = {
        {96, "Kstar_PiRatio"},
        {97, "Phi_PiRatio"}};

    vector<DataSet> data;
    data.reserve(tables.size());

    for (const auto &[table, name] : tables)
    {
        data.push_back(LoadTable(fALICE, table, name));
    }

    for (auto &d : data)
    {
        BuildErrorGraphs(d);
    }

    TFile *fOutput = new TFile("pp7TeVALICE.root", "RECREATE");

    for (const auto &d : data)
    {
        WriteGraphs(fOutput, d);
    }

    fOutput->Close();
    fALICE->Close();

    cout << "Output written to pp7TeVALICE.root" << endl;
}