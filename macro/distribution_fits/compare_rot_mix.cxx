#include <iostream>
#include <vector>
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

void DrawAndSaveComparison(TGraphErrors *graph_mix, TGraphErrors *graph_rot, double pdgValue, const string &outputfilename, const string &title)
{
    TCanvas *canvas = new TCanvas("", "", 720, 720);
    SetCanvasStyle(canvas, 0.15, 0.03, 0.05, 0.14);

    graph_mix->SetMarkerColor(kRed);
    graph_mix->SetLineColor(kRed);
    graph_mix->SetMarkerStyle(20);
    graph_mix->SetMarkerSize(1.5);
    if (title == "Width f1270")
    {
        graph_mix->SetMinimum(0.0);
        graph_mix->SetMaximum(0.4);
    }
    graph_mix->Draw("AP");

    graph_rot->SetMarkerColor(kBlue);
    graph_rot->SetLineColor(kBlue);
    graph_rot->SetMarkerStyle(21);
    graph_rot->SetMarkerSize(1.5);
    graph_rot->Draw("P SAME");

    TLine *linepdg = new TLine(0, pdgValue, 12, pdgValue);
    linepdg->SetLineStyle(2);
    linepdg->SetLineColor(kBlack);
    linepdg->SetLineWidth(2);
    linepdg->Draw();

    TLegend *legend = new TLegend(0.19, 0.79, 0.5, 0.93);
    SetLegendStyle(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(graph_mix, "MIX", "lpe");
    legend->AddEntry(graph_rot, "ROTATED", "lpe");
    legend->AddEntry(linepdg, ("PDG " + title).c_str(), "l");
    legend->Draw();

    canvas->SaveAs(outputfilename.c_str());
}

void compare_rot_mix()
{   
    gStyle->SetOptStat(0);
    string bkgmix = "MIX";
    string bkgrot = "ROTATED";
    const string outputfolder_str = "../" + kSignalOutput + "/" + kchannel + "/" + kfoldername;
    const string fits_folder_str = outputfolder_str + "/fits/";
    vector<string> param_names = {"massf1270", "widthf1270", "massf1525", "widthf1525", "massf1710", "widthf1710"};
    vector<double> pdgValues = {f1270Mass, f1270Width, f1525Mass, f1525Width, f1710Mass};
    vector<string> titles = {"Mass f1270", "Width f1270", "Mass f1525", "Width f1525", "Mass f1710"};

    TFile *fmix = new TFile((fits_folder_str + "fitparams_" + bkgmix + ".root").c_str(), "READ");
    TFile *frot = new TFile((fits_folder_str + "fitparams_" + bkgrot + ".root").c_str(), "READ");

    if (fmix->IsZombie() || frot->IsZombie())
    {
        cout << "Files not found" << endl;
        return;
    }

    vector<TGraphErrors *> graphs_mix;
    vector<TGraphErrors *> graphs_rot;

    for (const auto &name : param_names)
    {
        TGraphErrors *graph_mix = (TGraphErrors *)fmix->Get(name.c_str());
        TGraphErrors *graph_rot = (TGraphErrors *)frot->Get(name.c_str());

        if (graph_mix == nullptr || graph_rot == nullptr)
        {
            cout << "Graph " << name << " not found" << endl;
            return;
        }
        SetGrapherrorStyle(graph_mix);
        SetGrapherrorStyle(graph_rot);
        graphs_mix.push_back(graph_mix);
        graphs_rot.push_back(graph_rot);
    }

    for (size_t i = 0; i < param_names.size() - 1; i++)
    {
        string outputfilename = outputfolder_str + "/fits/compare_" + param_names[i] + ".png";
        DrawAndSaveComparison(graphs_mix[i], graphs_rot[i], pdgValues[i], outputfilename, titles[i]);
    }

    // //Compare invariant mass after ME subtraction in pass 6 and pass 7
    // TFile *fpass7 = new TFile("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/250337/KsKs_Channel/strangeness_tutorial/hglue_MIX.root", "READ");
    // TFile *fpass6 = new TFile("/home/sawan/check_k892/output/glueball/LHC220_pass6_small/230281/KsKs_Channel/strangeness_tutorial/hglue_MIX.root", "READ");
    // if(fpass7->IsZombie() || fpass6->IsZombie())
    // {
    //     cout << "Files not found" << endl;
    //     return;
    // }
    // TH1F *hinvpass7 = (TH1F *)fpass7->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    // TH1F *hinvpass6 = (TH1F *)fpass6->Get("ksks_subtracted_invmass_pt_0.0_30.0");
    // if(hinvpass7 == nullptr || hinvpass6 == nullptr)
    // {
    //     cout << "Histogram not found" << endl;
    //     return;
    // }
    // TCanvas *canvas = new TCanvas("", "", 720, 720);
    // SetCanvasStyle(canvas, 0.15, 0.03, 0.05, 0.14);
    // SetHistoQA(hinvpass7);
    // hinvpass7->SetMarkerColor(kRed);
    // hinvpass7->SetLineColor(kRed);
    // hinvpass7->SetMarkerStyle(20);
    // hinvpass7->GetYaxis()->SetRangeUser(-4000, 21000);
    // hinvpass7->Draw("EP");
    // SetHistoQA(hinvpass6);
    // hinvpass6->SetMarkerColor(kBlue);
    // hinvpass6->SetLineColor(kBlue);
    // hinvpass6->SetMarkerStyle(22);
    // hinvpass6->Draw("EP SAME");

    // TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    // SetLegendStyle(leg);
    // leg->AddEntry(hinvpass7, "Pass 7", "lpe");
    // leg->AddEntry(hinvpass6, "Pass 6", "lpe");
    // leg->Draw();
    // canvas->SaveAs("/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/250337/KsKs_Channel/strangeness_tutorial/compare_pass6_7.png");
    
}
