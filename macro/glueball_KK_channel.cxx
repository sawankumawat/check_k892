
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

void glueball_KK_channel()

{
    const string kResBkg = "LIKE"; // "MIX" or "LIKE" or "ROTATED"
    const bool canvasgrids = 0;

    // Folder name inside the Analysis.root file *****************************************
    const string kfoldername = "phianalysisrun3";

    const int kRebin = 2;
    const float txtsize = 0.045;
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(110);

    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.06);
    t2->SetTextFont(42);
    if (canvasgrids == 1)
    {
        TCanvas *cgrid1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
        TCanvas *cgrid2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
        TCanvas *cgrid_bkg1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
        TCanvas *cgrid_bkg2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    }

    // Input file

    TFile *fInputFile = new TFile("/home/sawan/check_k892/data/pp/glueball/LHC22o_pass6/kk_channel/195807/AnalysisResults.root", "Read");
    if (fInputFile == nullptr)
    {
        cerr << "File not found " << endl;
        return;
    }

    TH1F *hentries = (TH1F *)fInputFile->Get("event-selection-task/hColCounterAcc");
    double Event = hentries->GetEntries();
    cout << "*****************number of events********************:" << Event << endl;

    //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************

    TH3F *fHistNum = (TH3F *)fInputFile->Get(Form("%s/h3PhiInvMassUnlikeSign", kfoldername.c_str()));
    TH3F *fHistDen = (TH3F *)fInputFile->Get(Form("%s/h3PhiInvMassMixed", kfoldername.c_str()));
    TH3F *fHistLS_pp = (TH3F *)fInputFile->Get(Form("%s/h3PhiInvMassLikeSignPP", kfoldername.c_str()));
    TH3F *fHistLS_mm = (TH3F *)fInputFile->Get(Form("%s/h3PhiInvMassLikeSignMM", kfoldername.c_str()));
    TH3F *fHist_rotated = (TH3F *)fInputFile->Get(Form("%s/h3PhiInvMassRotation", kfoldername.c_str()));
    if (fHistNum == nullptr || fHistDen == nullptr || fHistLS_pp == nullptr || fHistLS_mm == nullptr || fHist_rotated == nullptr)
    {
        cout << "One of the histograms is NULL" << endl;
        return;
    }

    // cout << " THE NUMBER OF BINS IN THE HISTOGRAM IS " << fHistNum->GetNbinsZ()<<endl;

    // gstyle(); // this is not gStyle, it is defined in the header file style.h
    TH1D *fHistTotal;
    TH1D *fHistBkg;

    int ptbinlow = fHistNum->GetYaxis()->FindBin(2 + 0.001);
    int ptbinhigh = fHistNum->GetYaxis()->FindBin(10 - 0.001);
    int multbinlow = fHistNum->GetXaxis()->FindBin(2 + 0.001);
    int multbinhigh = fHistNum->GetXaxis()->FindBin(100 - 0.001);

    fHistTotal = fHistNum->ProjectionZ("hSig", -1, -1, ptbinlow, ptbinhigh, "E"); // multiplicity, pT, inv mass
    TH1D *fHistBkg_ME = fHistDen->ProjectionZ("hbkg", -1, -1, ptbinlow, ptbinhigh, "E");
    TH1D *fHistBkg_pp = fHistLS_pp->ProjectionZ("hbkg_pp", -1, -1, ptbinlow, ptbinhigh, "E");
    TH1D *fHistBkg_mm = fHistLS_mm->ProjectionZ("hbkg_mm", -1, -1, ptbinlow, ptbinhigh, "E");
    TH1D *fHistBkg_rotated = fHist_rotated->ProjectionZ("hbkg_rotated", -1, -1, ptbinlow, ptbinhigh, "E");
    if (kResBkg == "MIX")
    {
        fHistBkg = fHistBkg_ME;
    }
    else if (kResBkg == "LIKE")
    {
        fHistBkg = fHistBkg_pp;
        for (int ibkg = 0; ibkg < fHistBkg_pp->GetNbinsX(); ibkg++)
        {
            fHistBkg->SetBinContent(ibkg + 1, 2 * sqrt(fHistBkg_pp->GetBinContent(ibkg + 1) * fHistBkg_mm->GetBinContent(ibkg + 1)));
        }
    }
    else
    {
        fHistBkg = fHistBkg_rotated;
    }

    auto binwidth_file = (fHistTotal->GetXaxis()->GetXmax() - fHistTotal->GetXaxis()->GetXmin()) * kRebin / fHistTotal->GetXaxis()->GetNbins();
    cout << "*********The bin width is:  " << binwidth_file << "*********" << endl;

    //**Cloning sig+bkg histogram for like sign or mixed event subtraction *********************************************************
    TH1D *hfsig = (TH1D *)fHistTotal->Clone();
    //*****************************************************************************************************************************
    float normalisationlow = 2.2;
    float normalisationhigh = 2.5;
    if (kResBkg == "MIX" || kResBkg == "ROTATED")
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
    else
    {
        hfbkg = (TH1D *)fHistBkg->Clone();
        hfbkg->Rebin(kRebin);
        hfsig->Rebin(kRebin);
        hfsig->Add(hfbkg, -1);
    }

    fHistTotal->Rebin(kRebin);

    //************************************Plotting from here****************************************************************
    TFile *fresults = new TFile(Form("/home/sawan/check_k892/output/pp/glueball/LHC22o_pass6/results_%s.root", kResBkg.c_str()), "recreate");
    TCanvas *c1 = new TCanvas("", "", 720, 720);
    SetCanvasStyle2(c1, 0.13, 0.03, 0.05, 0.13);
    SetHistoQA(hfsig);
    hfsig->SetTitle(0);
    hfsig->SetMarkerStyle(8);
    hfsig->SetMarkerColor(kBlack);
    hfsig->SetLineColor(kBlack);
    hfsig->GetXaxis()->SetTitle("m_{K^{#pm}K^{#mp}} (GeV/c^{2})");
    hfsig->GetYaxis()->SetTitle("Counts");
    hfsig->GetYaxis()->SetTitleOffset(1.3);
    hfsig->Draw("e");
    hfsig->Write("h_after_bkg_subraction");
    gPad->Update();
    // TPaveStats *ps = (TPaveStats *)hfsig->FindObject("stats");
    // if (ps)
    // {
    //     ps->SetTextSize(0.04);
    //     ps->SetTextFont(42);
    //     ps->SetX1NDC(0.6);
    //     ps->SetX2NDC(0.95);
    //     ps->SetY1NDC(0.35);
    //     ps->SetY2NDC(0.95);
    // }
    // gPad->Modified(); // Necessary to update the canvas with the new text size
    // gPad->Update();

    TF1 *fitBW = new TF1("fitBW", gluefit2bW, 1, 2.3, 9);
    // TF1 *fitBW = new TF1("fitBW", gluefit3bW, 1, 2.4, 12);
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
    // fitBW->SetParameter(6, 1.2 * f1710norm);
    // fitBW->SetParameter(7, f1710Mass);
    // fitBW->SetParameter(8, f1710Width);

    // fitBW->FixParameter(0, f1270norm);
    // fitBW->FixParameter(1, f1270Mass);
    fitBW->FixParameter(2, f1270Width);
    // fitBW->FixParameter(3, f1525norm);
    // fitBW->FixParameter(4, f1525Mass);
    // fitBW->FixParameter(5, f1525Width);
    // fitBW->SetParLimits(4, f1525Mass - f1525Width * 3, f1525Mass + f1525Width * 3);
    // fitBW->FixParameter(7, f1710Mass);
    // fitBW->FixParameter(8, f1710Width);
    // fitBW->SetParLimits(7, f1710Mass - f1710Width * 3, f1710Mass + f1525Width * 3);

    // hfsig->Fit("fitBW", "REBMS");
    // double *par = fitBW->GetParameters();
    // TF1 *Bw1 = new TF1("Bw1", RelativisticBW, 1, 2.3, 3);
    // TF1 *Bw2 = new TF1("Bw2", RelativisticBW, 1, 2.3, 3);
    // TF1 *Bw3 = new TF1("Bw3", RelativisticBW, 1, 2.3, 3);
    // TF1 *expo = new TF1("expo", exponential_bkg, 1, 2.3, 3);
    // // fitBW->Draw("same");
    // Bw1->SetParameters(&par[0]);
    // Bw2->SetParameters(&par[3]);
    // Bw3->SetParameters(&par[6]);
    // expo->SetParameters(&par[9]);
    // Bw1->SetLineColor(28);
    // Bw2->SetLineColor(6);
    // Bw3->SetLineColor(7);
    // expo->SetLineColor(4);
    // Bw1->Draw("same");
    // Bw2->Draw("same");
    // // Bw3->Draw("same");
    // expo->Draw("same");

    // TLegend *lfit = new TLegend(0.3, 0.65, 0.55, 0.94);
    // lfit->SetFillColor(0);
    // // lfit->SetBorderSize(0);
    // lfit->SetFillStyle(0);
    // lfit->SetTextFont(42);
    // lfit->SetTextSize(0.04);
    // lfit->AddEntry(hfsig, "Data", "lpe");
    // lfit->AddEntry(fitBW, "3rBW+expol", "l");
    // lfit->AddEntry(Bw1, "rBW(a_{2}(1320))", "l");
    // lfit->AddEntry(Bw2, "rBW(f_{2}(1525))", "l");
    // lfit->AddEntry(Bw3, "rBW(f_{0}(1710))", "l");
    // lfit->AddEntry(expo, "Expol", "l");
    // lfit->Draw();

    c1->SaveAs("/home/sawan/check_k892/output/pp/glueball/LHC22o_pass6/glueball_inv_sub.pdf");

    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle2(c2, 0.15, 0.03, 0.05, 0.13);
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
    fHistTotal->GetXaxis()->SetTitle("m_{K^{#pm}K^{#mp}} (GeV/c^{2})");
    fHistTotal->GetYaxis()->SetTitleOffset(1.4);
    hfbkg->Draw("E same");
    if (kResBkg == "LIKE")
    {
        SetHistoQA(fHistBkg_pp);
        SetHistoQA(fHistBkg_mm);
        fHistBkg_pp->SetMarkerStyle(20);
        fHistBkg_pp->SetMarkerColor(kBlue);
        fHistBkg_pp->Draw("E same");
        fHistBkg_mm->SetMarkerStyle(20);
        fHistBkg_mm->SetMarkerColor(kGreen);
        fHistBkg_mm->Draw("E same");
    }
    if (kResBkg == "MIX" || kResBkg == "ROTATED")
    {
        hbkg_nopeak->Draw("BAR same");
    }
    TLegend *l = new TLegend(0.55, 0.4, 0.86, 0.7);
    SetLegendStyle(l);
    l->AddEntry(fHistTotal, "Signal+bkg", "lpe");
    if (kResBkg == "MIX")
    {
        l->AddEntry(hfbkg, "Mixed event bkg", "lpe");
    }
    else if (kResBkg == "LIKE")
    {
        l->AddEntry(hfbkg, "2 #sqrt{K^{+}K^{+} #times K^{-}K^{-}}", "lpe");
        l->AddEntry(fHistBkg_pp, "K^{+}K^{+}", "lpe");
        l->AddEntry(fHistBkg_mm, "K^{-}K^{-}", "lpe");
    }
    else
    {
        l->AddEntry(hfbkg, "Rotated bkg", "lpe");
    }
    l->Draw();
    c2->Write("before_bkg_subtraction");

    c2->SaveAs("/home/sawan/check_k892/output/pp/glueball/LHC22o_pass6/glueball_inv.pdf");

    ////////////////////////////////////////////////////////////////////////

    // } // pt loop ends
}
