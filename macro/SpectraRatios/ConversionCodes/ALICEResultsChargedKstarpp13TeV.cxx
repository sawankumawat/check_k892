#include <iostream>
#include <vector>
#include "../../src/style.h"

using namespace std;

struct DataSet
{

    TString table = "-1";
    TString name;

    TGraphErrors *g = nullptr;
    TH1D *hStat = nullptr;
    TH1D *hSys = nullptr;
    TH1D *hSysUncorr = nullptr;

    TGraphErrors *gStat = nullptr;
    TGraphErrors *gSys = nullptr;
    TGraphErrors *gSysUncorr = nullptr;
};

DataSet LoadTable(TFile *file, const TString &table, const TString &name, int GraphNumber = 1);
void BuildErrorGraphs(DataSet &d);
void WriteGraphs(TFile *outFile, const DataSet &d);

void ALICEResultsChargedKstarpp13TeV()
{
    TFile *fALICE = new TFile("../HEP_data/pp13TeV_ChargedKstar.root", "READ");
    if (!fALICE || fALICE->IsZombie())
    {
        cout << "Error: ALICE file not found" << endl;
        return;
    }

    vector<tuple<TString, int, TString>> tables = {
        {"19", 1, "ChKstar_MeanYield"},
        {"20", 1, "ChKstar_MeanpT"}};

    vector<DataSet> data;
    data.reserve(tables.size());

    for (const auto &[table, graphNumber, name] : tables)
    {
        data.push_back(LoadTable(fALICE, table, name, graphNumber));
    }

    for (auto &d : data)
    {
        BuildErrorGraphs(d);
    }

    TFile *fOutput = new TFile("pp13TeV_ChKstar.root", "RECREATE");

    for (const auto &d : data)
    {
        WriteGraphs(fOutput, d);
    }

    fOutput->Close();
    fALICE->Close();

    cout << "Output written to pp13TeVALICE.root" << endl;
}

DataSet LoadTable(TFile *file, const TString &table, const TString &name, int GraphNumber = 1)
{
    DataSet d;

    d.table = table;
    d.name = name;

    d.g = (TGraphErrors *)file->Get(Form("Table %s/Graph1D_y%d", table.Data(), GraphNumber));
    d.hStat = (TH1D *)file->Get(Form("Table %s/Hist1D_y%d_e1", table.Data(), GraphNumber));
    d.hSys = (TH1D *)file->Get(Form("Table %s/Hist1D_y%d_e2", table.Data(), GraphNumber));
    d.hSysUncorr = (TH1D *)file->Get(Form("Table %s/Hist1D_y%d_e3", table.Data(), GraphNumber));

    if (d.hSys == nullptr)
    {
        TH1D *hPlus = (TH1D *)file->Get(Form("Table %s/Hist1D_y%d_e2plus", table.Data(), GraphNumber));
        TH1D *hMinus = (TH1D *)file->Get(Form("Table %s/Hist1D_y%d_e2minus", table.Data(), GraphNumber));
        hPlus->Add(hMinus, -1.0);
        hPlus->Scale(0.5);
        d.hSys = hPlus;
    }

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
