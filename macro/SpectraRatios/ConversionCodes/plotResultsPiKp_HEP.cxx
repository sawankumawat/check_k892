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

void plotResultsPiKp_HEP()
{
    //========================Important information========================
    // In this code I have manually averaged the Pi,K,p yield by dividing by 2.0.
    //=========================================================================

    TFile *fALICE = new TFile("../HEP_data/HEPDataPiKp_Run3.root", "READ");
    if (!fALICE || fALICE->IsZombie())
    {
        cout << "Error: ALICE file not found" << endl;
        return;
    }

    vector<tuple<TString, int, TString>> tables = {
        {"17", 1, "Pion_MeanpT"},
        {"17", 2, "Kaon_MeanpT"},
        {"17", 3, "Proton_MeanpT"},

        {"20", 1, "Pion_MeanYield"},
        {"20", 2, "Kaon_MeanYield"},
        {"20", 3, "Proton_MeanYield"}};

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

    TFile *fOutput = new TFile("ppRun3_PiKpHEP.root", "RECREATE");

    for (const auto &d : data)
    {
        WriteGraphs(fOutput, d);
    }

    fOutput->Close();
    fALICE->Close();

    cout << "Output written to ppRun3_PiKpHEP.root" << endl;
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

        // double ex = d.g->GetErrorX(i);
        // double eyStat = d.hStat->GetBinContent(bin);
        // double eySys = d.hSys->GetBinContent(bin);

        double ex = d.g->GetErrorX(i);

        // Average the charged-particle yields
        if (d.name.Contains("MeanYield"))
        {
            y /= 2.0;
        }

        double eyStat = d.hStat->GetBinContent(bin);
        double eySys = d.hSys->GetBinContent(bin);

        // Errors also scale by the same factor
        if (d.name.Contains("MeanYield"))
        {
            eyStat /= 2.0;
            eySys /= 2.0;
        }

        d.gStat->SetPoint(i, x, y);
        d.gStat->SetPointError(i, ex, eyStat);

        d.gSys->SetPoint(i, x, y);
        d.gSys->SetPointError(i, ex, eySys);

        // // Optional e3
        // if (d.hSysUncorr && d.gSysUncorr)
        // {
        //     double eySysUncorr = d.hSysUncorr->GetBinContent(bin);

        //     d.gSysUncorr->SetPoint(i, x, y);
        //     d.gSysUncorr->SetPointError(i, ex, eySysUncorr);
        // }

        if (d.hSysUncorr && d.gSysUncorr)
        {
            double eySysUncorr = d.hSysUncorr->GetBinContent(bin);

            if (d.name.Contains("MeanYield"))
                eySysUncorr /= 2.0;

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
