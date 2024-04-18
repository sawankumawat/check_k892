
#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/initializations.h"
using namespace std;

float parameter0(float mass, float width)
{
    double gamma = TMath::Sqrt(mass * mass * (mass * mass + width * width));
    double norm = 2.8284 * mass * width * gamma / (3.14 * TMath::Sqrt(mass * mass + gamma));
    return norm;
}

void glueball_KsKs_channel()

{
    const string kResBkg = "MIX";

    // Folder name inside the Analysis.root file *****************************************
    const string kfoldername = "strangeness_tutorial/hglueball";
    // const string kfoldername = "strangeness_tutorial_tight_cuts/hglueball";
    string folderpath = "/home/sawan/check_k892/data/pbpb/glueball/";
    // string filepath = "LHC23zzh_pass2";
    string filepath = "LHC23zs_pass3_QC1";
    // string filepath = "LHC23zzf_pass2_QC";
    // string filepath = "LHC23zzg_pass_QC";

    const int kRebin = 2;
    const float txtsize = 0.045;
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(110);

    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.06);
    t2->SetTextFont(42);

    cgrid1->Close();
    cgrid2->Close();
    cgrid_bkg1->Close();
    cgrid_bkg2->Close();

    // Input file

    TFile *fInputFile = new TFile((folderpath + filepath + "/AnalysisResults.root").c_str(), "Read");
    // TFile *fInputFile = new TFile(kDataFilename.c_str(), "Read");
    if (!fInputFile)
    {
        cerr << "File not found " << endl;
        return;
    }

    TH1F *hentries = (TH1F *)fInputFile->Get("event-selection-task/hColCounterAcc");
    double Event = hentries->GetEntries();
    cout << "*****************number of events********************:" << Event << endl;

    //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************

    TH3F *fHistNum = (TH3F *)fInputFile->Get(Form("%s/h3glueInvMassDS", kfoldername.c_str()));
    TH3F *fHistDen = (TH3F *)fInputFile->Get(Form("%s/h3glueInvMassME", kfoldername.c_str()));

    // cout << " THE NUMBER OF BINS IN THE HISTOGRAM IS " << fHistNum->GetNbinsZ()<<endl;

    // gstyle(); // this is not gStyle, it is defined in the header file style.h
    TH1D *fHistTotal;
    TH1D *fHistBkg;

    int ptbinlow = fHistNum->GetYaxis()->FindBin(0.2 + 0.001);
    int ptbinhigh = fHistNum->GetYaxis()->FindBin(10 - 0.001);
    int multbinlow = fHistNum->GetXaxis()->FindBin(2 + 0.001);
    int multbinhigh = fHistNum->GetXaxis()->FindBin(100 - 0.001);

    fHistTotal = fHistNum->ProjectionZ("hSig", -1, -1, ptbinlow, ptbinhigh, "E"); // multiplicity, pT, inv mass
    fHistBkg = fHistDen->ProjectionZ("hbkg", -1, -1, ptbinlow, ptbinhigh, "E");
    //  fHistTotal = fHistNum->ProjectionZ("hSig", -1, -1, -1, -1, "E"); // multiplicity, pT, inv mass
    // fHistBkg = fHistDen->ProjectionZ("hbkg", -1, -1, -1, -1, "E");

    auto binwidth_file = (fHistTotal->GetXaxis()->GetXmax() - fHistTotal->GetXaxis()->GetXmin()) * kRebin / fHistTotal->GetXaxis()->GetNbins();
    cout<<"*********The bin width is:  "<<binwidth_file<<"*********"<<endl;

    //**Cloning sig+bkg histogram for like sign or mixed event subtraction *********************************************************
    TH1D *hfsig = (TH1D *)fHistTotal->Clone();
    //*****************************************************************************************************************************
    float normalisationlow = 2.2;
    float normalisationhigh = 2.3;
    if (kResBkg == "MIX")
    {
        sigbkg_integral = (fHistTotal->Integral(fHistTotal->GetXaxis()->FindBin(normalisationlow), fHistTotal->GetXaxis()->FindBin(normalisationhigh)));
        bkg_integral = (fHistBkg->Integral(fHistBkg->GetXaxis()->FindBin(normalisationlow), fHistBkg->GetXaxis()->FindBin(normalisationhigh)));
        normfactor = sigbkg_integral / bkg_integral; // scaling factor for mixed bkg
        cout << "\n\n normalization factor " << 1. / normfactor << "\n\n";
        hfbkg = (TH1D *)fHistBkg->Clone();

        hfbkg->Scale(normfactor);
        hfbkg->Rebin(kRebin);
        hfsig->Rebin(kRebin);

        hfsig->Add(hfbkg, -1);
    }

    fHistTotal->Rebin(kRebin);

    //*****************************************************************************************************
    TCanvas *c1 = new TCanvas("", "", 1400, 1000);
    SetCanvasStyle2(c1, 0.15, 0.03, 0.01, 0.15);
    SetHistoQA(hfsig);
    hfsig->SetTitle(0);
    hfsig->SetMarkerStyle(8);
    hfsig->SetMarkerColor(kBlack);
    hfsig->SetLineColor(kBlack);
    hfsig->GetXaxis()->SetTitle("m_{K_{s}K_{s}} (GeV/c^{2})");
    hfsig->GetYaxis()->SetTitle("Counts");
    hfsig->Draw("e");
    gPad->Update();
    TPaveStats *ps = (TPaveStats *)hfsig->FindObject("stats");
    if (ps)
    {
        ps->SetTextSize(0.04);
        ps->SetTextFont(42);
        ps->SetX1NDC(0.6);
        ps->SetX2NDC(0.95);
        ps->SetY1NDC(0.35);
        ps->SetY2NDC(0.95);
    }
    gPad->Modified(); // Necessary to update the canvas with the new text size
    gPad->Update();

    // TF1 *fitBW = new TF1("fitBW", gluefit2bW, 1, 2.3, 9);
    TF1 *fitBW = new TF1("fitBW", gluefit3bW, 1, 2.4, 12);
    fitBW->SetParNames("Norm1", "mass1", "width1", "Norm2", "mass2", "width2", "norm3", "mass3", "width3", "p1", "p2");
    fitBW->SetParName(11, "p3");
    float f1270norm = parameter0(f1270Mass, f1270Width);
    float a1320norm = parameter0(a1320Mass, a1320Width);
    float f1525norm = parameter0(f1525Mass, f1525Width);
    float f1710norm = parameter0(f1710Mass, f1710Width);

    fitBW->SetParameter(0, 1.5 * f1270norm);
    fitBW->SetParameter(1, f1270Mass);
    fitBW->SetParameter(2, f1270Width);
    fitBW->SetParameter(3, 1.2 * f1525norm);
    fitBW->SetParameter(4, f1525Mass);
    fitBW->SetParameter(5, f1525Width);
    fitBW->SetParameter(6, 1.2 * f1710norm);
    fitBW->SetParameter(7, f1710Mass);
    fitBW->SetParameter(8, f1710Width);

    // fitBW->FixParameter(0, f1270norm);
    // fitBW->FixParameter(1, f1270Mass);
    fitBW->FixParameter(2, f1270Width);
    // fitBW->FixParameter(3, f1525norm);
    // fitBW->FixParameter(4, f1525Mass);
    fitBW->FixParameter(5, f1525Width);
    // fitBW->SetParLimits(4, f1525Mass - f1525Width * 3, f1525Mass + f1525Width * 3);
    fitBW->FixParameter(7, f1710Mass);
    // fitBW->FixParameter(8, f1710Width);
    // fitBW->SetParLimits(7, f1710Mass - f1710Width * 3, f1710Mass + f1525Width * 3);


    hfsig->Fit("fitBW", "REBMS");
    double *par = fitBW->GetParameters();
    TF1 *Bw1 = new TF1("Bw1", RelativisticBW, 1, 2.3, 3);
    TF1 *Bw2 = new TF1("Bw2", RelativisticBW, 1, 2.3, 3);
    TF1 *Bw3 = new TF1("Bw3", RelativisticBW, 1, 2.3, 3);
    TF1 *expo = new TF1("expo", exponential_bkg, 1, 2.3, 3);
    // fitBW->Draw("same");
    Bw1->SetParameters(&par[0]);
    Bw2->SetParameters(&par[3]);
    Bw3->SetParameters(&par[6]);
    expo->SetParameters(&par[9]);
    Bw1->SetLineColor(28);
    Bw2->SetLineColor(6);
    Bw3->SetLineColor(7);
    expo->SetLineColor(4);
    Bw1->Draw("same");
    Bw2->Draw("same");
    Bw3->Draw("same");
    expo->Draw("same");

    TLegend *lfit = new TLegend(0.3, 0.65, 0.55, 0.94);
    lfit->SetFillColor(0);
    // lfit->SetBorderSize(0);
    lfit->SetFillStyle(0);
    lfit->SetTextFont(42);
    lfit->SetTextSize(0.04);
    lfit->AddEntry(hfsig, "Data", "lpe");
    lfit->AddEntry(fitBW, "3rBW+expol", "l");
    lfit->AddEntry(Bw1, "rBW(a_{2}(1320))", "l");
    lfit->AddEntry(Bw2, "rBW(f_{2}(1525))", "l");
    lfit->AddEntry(Bw3, "rBW(f_{0}(1710))", "l");
    lfit->AddEntry(expo, "Expol", "l");
    lfit->Draw();

    c1->SaveAs(Form("%s/hglueball.png", kSignalOutput.c_str()));

    TCanvas *c2 = new TCanvas("", "", 1400, 1000);
    SetCanvasStyle2(c2, 0.15, 0.03, 0.01, 0.15);
    SetHistoQA(fHistTotal);
    SetHistoQA(hfbkg);

    TH1F *hbkg_nopeak = (TH1F *)hfbkg->Clone();
    hbkg_nopeak->SetLineColor(kRed);
    hbkg_nopeak->SetMarkerColor(kRed);
    hbkg_nopeak->SetFillColor(kRed);
    hbkg_nopeak->SetFillStyle(3001);
    for (int i = 0; i < hbkg_nopeak->GetNbinsX(); i++)
    {
        if (hbkg_nopeak->GetBinCenter(i + 1) < normalisationlow || hbkg_nopeak->GetBinCenter(i + 1) > normalisationhigh)
        {
            hbkg_nopeak->SetBinContent(i + 1, -999);
        }
    }

    fHistTotal->SetMarkerStyle(8);
    fHistTotal->SetMarkerColor(kBlack);
    hfbkg->SetMarkerStyle(8);
    hfbkg->SetMarkerColor(kRed);
    fHistTotal->Draw("E");
    fHistTotal->GetYaxis()->SetTitle("Counts");
    hfbkg->Draw("E same");
    hbkg_nopeak->Draw("BAR same");

    c2->SaveAs(Form("%s/hglueball_bkg.png", kSignalOutput.c_str()));

    ////////////////////////////////////////////////////////////////////////

    // } // pt loop ends
}
