#include <iostream>
#include "../macro/src/style.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void analyze_glue()
{
    gStyle->SetOptStat(0);
    TFile *f = new TFile("merged_output_glue.root");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    THnSparse *hdirect = (THnSparse *)f->Get("hKsKs_mass");
    THnSparse *hrot = (THnSparse *)f->Get("hKsKs_rotated");
    if (hdirect == nullptr || hrot == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }
    int lowbinpt = hdirect->GetAxis(1)->FindBin(0.0);
    int highbinpt = hdirect->GetAxis(1)->FindBin(30.0);
    hdirect->GetAxis(1)->SetRange(lowbinpt, highbinpt);
    TH1D *hdirectMass = hdirect->Projection(0, "E");
    hrot->GetAxis(1)->SetRange(lowbinpt, highbinpt);
    TH1D *hrotMass = hrot->Projection(0, "E");

    // Scale the rotational histogram
    double normalization_region[] = {2.2, 2.4};
    int lowbin = hrotMass->FindBin(normalization_region[0]);
    int highbin = hrotMass->FindBin(normalization_region[1]);
    double area_rot = hrotMass->Integral(lowbin, highbin);
    double area_direct = hdirectMass->Integral(lowbin, highbin);
    double scale_factor = area_direct / area_rot;
    hrotMass->Scale(scale_factor);

    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.14, 0.02, 0.06, 0.14);
    SetHistoQA(hdirectMass);
    SetHistoQA(hrotMass);
    hrotMass->SetLineColor(kRed);
    hrotMass->SetMarkerColor(kRed);
    hdirectMass->SetLineColor(kBlue);
    hdirectMass->SetMarkerColor(kBlue);
    hdirectMass->GetXaxis()->SetTitle("m_{K_{S}K_{S}} [GeV/c^{2}]");
    hdirectMass->GetYaxis()->SetTitle("Counts");
    hdirectMass->SetMarkerSize(0.5);
    hrotMass->SetMarkerSize(0.5);
    hdirectMass->Draw();
    hrotMass->Draw("same");
    TH1F *hnorm_region = (TH1F *)hdirectMass->Clone();
    for (int i = 0; i < hnorm_region->GetNbinsX(); i++)
    {
        if (hnorm_region->GetBinCenter(i) < normalization_region[0] || hnorm_region->GetBinCenter(i) > normalization_region[1])
            hnorm_region->SetBinContent(i, -999);
    }
    hnorm_region->SetFillColor(8);
    hnorm_region->SetFillStyle(3001);
    hnorm_region->SetLineWidth(0);
    hnorm_region->Draw("HIST same");

    TLegend *leg = new TLegend(0.5, 0.65, 0.9, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->AddEntry((TObject *)0, "Pythia simulation pp 13.6 TeV", "");
    leg->AddEntry(hdirectMass, "KsKs pairs", "l");
    leg->AddEntry(hrotMass, "Rotational background", "l");
    leg->AddEntry(hnorm_region, "Normalization region", "f");
    leg->Draw();

    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.14, 0.02, 0.06, 0.14);
    TH1F *hdirect_subtracted = (TH1F *)hdirectMass->Clone();
    hdirect_subtracted->Add(hrotMass, -1);
    SetHistoQA(hdirect_subtracted);
    hdirect_subtracted->GetXaxis()->SetTitle("m_{K_{S}K_{S}} [GeV/c^{2}]");
    hdirect_subtracted->GetYaxis()->SetTitle("Counts");
    hdirect_subtracted->Draw();
}