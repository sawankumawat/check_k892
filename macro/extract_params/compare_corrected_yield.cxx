#include <iostream>
#include "../src/style.h"
#include "../src/common.h"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);

    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.009);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.002);

    pad2->SetLeftMargin(0.1392758);
    pad2->SetRightMargin(0.009749304);
    pad2->SetTopMargin(0.1714559);
    pad2->SetBottomMargin(0.3304598);
}

void compare_corrected_yield()
{
    // published spectra
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TFile *fpub = new TFile("pp13TeV.root", "READ");
    TGraph *published_spectra = (TGraph *)fpub->Get("Table 4/Graph1D_y1");
    if (!published_spectra)
    {
        cout << "Graph not found" << endl;
        return;
    }
    int spectrano = 0;

    //*********Root files*************
    // TFile *fspectra = new TFile("/home/sawan/k892_postprocessing/output/LHC22o_pass6_small_INEL/common/spectra_1_LSBkg_0.root", "READ");
    TFile *fspectraqa = new TFile(("../" + koutputfolder + "/yield.root").c_str(), "READ"); // for raw yield
    // TFile *fefficiency = new TFile("/home/sawan/k892_postprocessing/output/LHC22o_pass6_small_INEL/common/efficiency_1.root", "READ"); // efficiency
    TFile *correc_spectra = new TFile(("../" + koutputfolder + "/efficiency/corrected_spectra.root").c_str(), "READ"); //corrected spectra
    if (fspectraqa->IsZombie() || correc_spectra->IsZombie())
    {
        cout << "File not found" << endl;
        return;
    }


    // TH1D *hrawbong = (TH1D *)fspectra->Get("lf-k892analysis/K892/0/hraw1Yields"); //raw yields from bong code
    // TH1D *hcorrectedbong = (TH1D *)fspectra->Get("lf-k892analysis/K892/0/hCorrectedYieldsWithoutLossFactors"); //corrected yields from bong code (without the correction factors)
    TH1D *hqa = (TH1D *)fspectraqa->Get("yield_integral");  //raw yields from this analysis
    // TH1D *heff = (TH1D *)fefficiency->Get("lf-k892analysis/hEfficiency_cen0"); //efficiency from Bong code
    TH1D *hcorrspec = (TH1D *)correc_spectra->Get("corrected_spectra"); //corrected spectra for this analysis

    if (hqa == nullptr || hcorrspec == nullptr)
    {
        cout << "Histogram not found" << endl;
        return;
    }


    // TH1F *hcorrectedspectra = (TH1F *)hqa->Clone();
    // hcorrectedspectra->Divide(heff);
    TH1F *hcorrectedspectra = (TH1F *)hcorrspec->Clone();

     TGraphErrors *gratio = new TGraphErrors();
    for (int i = 0; i < published_spectra->GetN(); i++)
    {
        double x, y;
        published_spectra->GetPoint(i, x, y);
        // cout << "The value of x of published is " << x << endl;
        double published = y;
        double thisanalysis = hcorrectedspectra->GetBinContent(i + 1);
        double published_error = published_spectra->GetErrorY(i);
        double thisanalysis_error = hcorrectedspectra->GetBinError(i + 1);
        // cout << "The value of x axis of this analysis is " << hcorrectedspectra->GetBinCenter(i + 1) << endl;
        double yieldratio = thisanalysis / published;
        cout<<"The value of yield ratio is "<<yieldratio<<endl;
        gratio->SetPoint(i, x, yieldratio);
        gratio->SetPointError(i, 0, (1.0 / published) * sqrt(thisanalysis_error * thisanalysis_error + (yieldratio * yieldratio) * (published_error * published_error)));
    }
  

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    c1->SetLogy(1);
    gPad->SetLogy(1);
    hcorrectedspectra->GetYaxis()->SetTitle(hcorrectedspectra->GetYaxis()->GetTitle());
    hcorrectedspectra->GetYaxis()->CenterTitle(1);
    hcorrectedspectra->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hcorrectedspectra->SetTitle(0);
    SetHistoStyle(hcorrectedspectra, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    hcorrectedspectra->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hcorrectedspectra->GetXaxis()->SetTitleOffset(1.02);
    hcorrectedspectra->GetYaxis()->SetTitleOffset(1.15);
    hcorrectedspectra->SetMinimum(1e-6);
    hcorrectedspectra->SetMarkerStyle(29);
    hcorrectedspectra->SetMarkerSize(1);
    hcorrectedspectra->SetMarkerColor(kRed);
    hcorrectedspectra->SetLineColor(kRed);
    hcorrectedspectra->GetXaxis()->SetRangeUser(0, 15);
    hcorrectedspectra->Draw("pe");
    SetgrStyle(published_spectra, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    published_spectra->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    published_spectra->SetMaximum(published_spectra->GetMaximum() * 2);
    published_spectra->GetYaxis()->SetTitleOffset(1.15);
    published_spectra->GetXaxis()->SetTitleOffset(1.02);
    published_spectra->SetMarkerStyle(22);
    published_spectra->SetMarkerSize(1);
    published_spectra->SetLineColor(1);
    published_spectra->SetMarkerColor(1);
    published_spectra->Draw("pe same");
    // hcorrectedbong->SetMarkerStyle(20);
    // hcorrectedbong->SetMarkerSize(1);
    // hcorrectedbong->SetMarkerColor(kBlue);
    // hcorrectedbong->SetLineColor(kBlue);
    // hcorrectedbong->Draw("pe same");
  

    TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.8);
    SetLegendStyle(leg);
    leg->AddEntry(published_spectra, " (published)", "lpe");
    leg->AddEntry(hcorrectedspectra, "(this analysis)", "lpe");
    leg->SetTextSize(0.05);
    leg->Draw();

    c1->cd(2);

    SetgrStyle(gratio, 1, 53, 1, 0.05, 0.05, 0.04 / pad2Size, 0.04 / pad2Size, 1.13, 1.8);
    gratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    gratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    gratio->SetMarkerStyle(23);
    gratio->SetMarkerSize(1.0);
    gratio->SetMarkerColor(kRed);
    gratio->SetLineColor(kRed);
    gratio->GetYaxis()->SetTitle("#frac{Corrected yield}{Published}");
    gratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gratio->GetXaxis()->CenterTitle(1);
    gratio->GetYaxis()->SetTitleOffset(0.45);
    gratio->GetYaxis()->SetNdivisions(505);
    gratio->GetXaxis()->SetRangeUser(0, 15);
    gratio->Draw("ap");

    // TCanvas *cpub = new TCanvas("cpub", "cpub", 720, 720);
    // gPad->SetLogy(1);
    // published_spectra->SetMarkerStyle(20);
    // published_spectra->SetMinimum(1e-6);
    // published_spectra->SetMaximum(1);
    // published_spectra->Draw("ape");
}