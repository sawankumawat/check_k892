#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

void plot_mass_width_ptwise_temp()
{
    float pt[] = {1, 2, 3, 4, 6, 10};
    float ptmean[] = {1.5, 2.5, 3.5, 5, 8};
    float binwidth[] = {0.5, 0.5, 0.5, 1, 2};
    const int noof_bins = (sizeof(pt) / sizeof(pt[0])) - 1;

    float mass1270[] = {1.194, 1.182, 1.247, 1.283, 1.287};
    float mass1270_err[] = {0.008, 0.010, 0.052, 0.009, 0.012};

    float mass1525[] = {1.508, 1.512, 1.513, 1.512, 1.517};
    float mass1525_err[] = {0.001, 0.001, 0.001, 0.001, 0.001};
    float width1525[] = {0.07253, 0.08429, 0.08777, 0.09489, 0.08765};
    float width1525_err[] = {0.00420, 0.00553, 0.00489, 0.00392, 0.0081};

    float mass1710[] = {1.778, 1.722, 1.683, 1.697, 1.693};
    float mass1710_err[] = {0.028, 0.013, 0.006, 0.004, 0.003};
    float width1710[] = {0.8477, 0.07214, 0.23, 0.1479, 0.1425};
    float width1710_err[] = {0.0898, 0.02846, 0.012, 0.0112, 0.0081};

    TGraphErrors *gr1270mass = new TGraphErrors(noof_bins, ptmean, mass1270, binwidth, mass1270_err);
    TGraphErrors *gr1525mass = new TGraphErrors(noof_bins, ptmean, mass1525, binwidth, mass1525_err);
    TGraphErrors *gr1525width = new TGraphErrors(noof_bins, ptmean, width1525, binwidth, width1525_err);
    TGraphErrors *gr1710mass = new TGraphErrors(noof_bins, ptmean, mass1710, binwidth, mass1710_err);
    TGraphErrors *gr1710width = new TGraphErrors(noof_bins, ptmean, width1710, binwidth, width1710_err);

    vector<TGraphErrors *> graphs = {gr1270mass, gr1525mass, gr1525width, gr1710mass, gr1710width};
    string canvasnames[] = {"mass1270", "mass1525", "width1525", "mass1710", "width1710"};
    float pdgvalues[] = {f1270Mass, f1525Mass, f1525Width, f1710Mass, f1710Width};
    float pdgvalueserr[] = {f1270MassErr, f1525MassErr, f1525WidthErr, f1710MassErr, f1710WidthErr};
    string yaxis_titles[] = {"Mass f1270 (GeV/c^{2})", "Mass f1525 (GeV/c^{2})", "Width f1525 (GeV/c^{2})", "Mass f1710 (GeV/c^{2})", "Width f1710 (GeV/c^{2})"};
    vector<vector<float>> plot_limits = {{1.1, 1.4}, {1.5, 1.53}, {0.06, 0.11}, {1.649, 1.849}, {0.0, 0.29}};

    for (int i = 0; i < graphs.size(); i++)
    {
        SetGrapherrorStyle(graphs[i]);
        graphs[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        graphs[i]->GetYaxis()->SetTitle(yaxis_titles[i].c_str());
        graphs[i]->GetXaxis()->SetLimits(0, 10);
        graphs[i]->GetYaxis()->SetMaxDigits(3);
        graphs[i]->GetYaxis()->SetRangeUser(plot_limits[i][0], plot_limits[i][1]);
        TCanvas *c = new TCanvas(canvasnames[i].c_str(), canvasnames[i].c_str(), 720, 720);
        SetCanvasStyle(c, 0.15, 0.03, 0.05, 0.14);
        graphs[i]->Draw("AP");

        TLine *line = new TLine(0, pdgvalues[i], 10, pdgvalues[i]);
        line->SetLineColor(2);
        // line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");
        // box of pdg error around the pdg line
        TBox *box = new TBox(0, pdgvalues[i] - pdgvalueserr[i] / 2, 10, pdgvalues[i] + pdgvalueserr[i] / 2);
        box->SetFillColor(3);
        box->SetFillStyle(3003);
        box->Draw("same");

        TLegend *l = new TLegend(0.70, 0.75, 0.9, 0.92);
        l->SetFillStyle(0);
        l->SetTextFont(42);
        l->SetBorderSize(0);
        l->SetTextSize(0.03);
        l->AddEntry(graphs[i], "Data", "lpe");
        l->AddEntry(line, "PDG value", "l");
        l->AddEntry(box, "PDG error", "f");
        l->Draw("same");

        string savepath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/260782/KsKs_Channel/strangeness_tutorial/fits/fit_func1_check2_code/ptfits";
        c->SaveAs((savepath + "/compare_"+ canvasnames[i] + ".png").c_str());
    }
}