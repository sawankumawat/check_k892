
#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/common_glue.h"
#include "src/fitting_range_glue.h"

using namespace std;

float parameter0(float mass, float width)
{
    double gamma = TMath::Sqrt(mass * mass * (mass * mass + width * width));
    double norm = 2.8284 * mass * width * gamma / (3.14 * TMath::Sqrt(mass * mass + gamma));
    return norm;
}

void glueball_KsKs_channel()

{
    // change here ***********************************************************
    const string kResBkg = "MIX";
    // const string kResBkg = "ROTATED";
    // change here ***********************************************************

    TString outputfolder = kSignalOutput + "/" + kfoldername;
    // Create the folder using TSystem::mkdir()
    if (gSystem->mkdir(outputfolder, kTRUE))
    {
        std::cout << "Folder " << outputfolder << " created successfully." << std::endl;
    }
    else
    {
        std::cerr << "Creating folder " << outputfolder << std::endl;
    }
    // Folder name inside the Analysis.root file *****************************************

    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(110);

    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.06);
    t2->SetTextFont(42);

    // Input file
    TFile *fInputFile = new TFile(kDataFilename.c_str(), "Read");
    if (fInputFile->IsZombie())
    {
        cerr << "File not found " << endl;
        return;
    }

    TH1F *hentries = (TH1F *)fInputFile->Get("event-selection-task/hColCounterAcc");
    double Event = hentries->GetEntries();
    cout << "*******number of events from the event selection histogram is *******:" << Event << endl;
    TH1F *hmult = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/eventSelection/hmultiplicity").c_str());
    if (hmult == nullptr)
    {
        cout << "Multiplicity histogram not found" << endl;
        return;
    }
    double realevents = hmult->Integral(hmult->GetXaxis()->FindBin(0.0), hmult->GetXaxis()->FindBin(100.0));
    cout << "*******number of events from the multiplicity histogram is *******:" << realevents << endl;

    //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************

    THnSparseF *fHistNum = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassDS", kfoldername.c_str()));
    THnSparseF *fHistDen = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassME", kfoldername.c_str()));
    THnSparseF *fHistRot = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassRot", kfoldername.c_str()));
    cout << " The number of entries in histograms: \n"
         << "same event: " << fHistNum->GetEntries() << "\n"
         << "mixed event: " << fHistDen->GetEntries() << "\n"
         << "rotated bkg/2: " << fHistRot->GetEntries()/2 << endl;

    if (fHistNum == nullptr || fHistDen == nullptr)
    {
        cout << "Invariant mass histograms not found" << endl;
        return;
    }

    TH1D *fHistTotal[Npt];
    TH1D *fHistBkg[Npt];
    TH1D *fHistRotated[Npt];

    for (Int_t ip = pt_start; ip < pt_end; ip++) // start pt bin loop
    {

        int lowpt = pT_bins[ip];
        int highpt = pT_bins[ip + 1];
        int multlow = 0;
        int multhigh = 100;
        int lbin = fHistNum->GetAxis(1)->FindBin(lowpt + 1e-5);
        int hbin = fHistNum->GetAxis(1)->FindBin(highpt - 1e-5);

        fHistNum->GetAxis(1)->SetRange(lbin, hbin);
        fHistDen->GetAxis(1)->SetRange(lbin, hbin);
        fHistRot->GetAxis(1)->SetRange(lbin, hbin);

        int lbinmult = fHistNum->GetAxis(0)->FindBin(multlow + 1e-5);
        int hbinmult = fHistNum->GetAxis(0)->FindBin(multhigh - 1e-5);

        fHistNum->GetAxis(0)->SetRange(lbinmult, hbinmult);
        fHistDen->GetAxis(0)->SetRange(lbinmult, hbinmult);
        fHistRot->GetAxis(0)->SetRange(lbinmult, hbinmult);

        fHistTotal[ip] = fHistNum->Projection(2, "E");
        fHistBkg[ip] = fHistDen->Projection(2, "E");
        fHistRotated[ip] = fHistRot->Projection(2, "E");
        fHistTotal[ip]->SetName(Form("fHistTotal_%d", ip));
        fHistBkg[ip]->SetName(Form("fHistBkg_%d", ip));
        fHistRotated[ip]->SetName(Form("fHistRotated_%d", ip));

        auto energylow = fHistTotal[ip]->GetXaxis()->GetXmin();
        auto energyhigh = fHistTotal[ip]->GetXaxis()->GetXmax();
        cout << "energy low value is " << energylow << endl;
        cout << "energy high value is " << energyhigh << endl;

        auto binwidth_file = (fHistTotal[ip]->GetXaxis()->GetXmax() - fHistTotal[ip]->GetXaxis()->GetXmin()) * kRebin[ip] / fHistTotal[ip]->GetXaxis()->GetNbins();
        cout << "*********The bin width is:  " << binwidth_file << "*********" << endl;

        //**Cloning sig+bkg histogram for like sign or mixed event subtraction *********************************************************
        TH1D *hfsig = (TH1D *)fHistTotal[ip]->Clone();
        TH1D *hfbkg;

        //*****************************************************************************************************************************
        float normalisationlow = 2.0;
        float normalisationhigh = 2.1;
        if (kResBkg == "MIX")
        {
            auto sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(normalisationlow), fHistTotal[ip]->GetXaxis()->FindBin(normalisationhigh)));
            auto bkg_integral = (fHistBkg[ip]->Integral(fHistBkg[ip]->GetXaxis()->FindBin(normalisationlow), fHistBkg[ip]->GetXaxis()->FindBin(normalisationhigh)));
            auto normfactor = sigbkg_integral / bkg_integral; // scaling factor for mixed bkg
            cout << "\n\n normalization factor " << 1. / normfactor << "\n\n";
            hfbkg = (TH1D *)fHistBkg[ip]->Clone();

            hfbkg->Scale(normfactor);
            hfbkg->Rebin(kRebin[ip]);
            hfsig->Rebin(kRebin[ip]);

            hfsig->Add(hfbkg, -1);
        }
        else if (kResBkg == "ROTATED")
        {
            hfbkg = (TH1D *)fHistRotated[ip]->Clone();
            hfbkg->Scale(0.5);
            hfbkg->Rebin(kRebin[ip]);
            hfsig->Rebin(kRebin[ip]);
            hfsig->Add(hfbkg, -1);
        }

        fHistTotal[ip]->Rebin(kRebin[ip]);

        //*****************************************************************************************************
        TCanvas *c1 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c1, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hfsig);
        hfsig->SetTitle(0);
        hfsig->SetMarkerStyle(8);
        hfsig->GetYaxis()->SetMaxDigits(3);
        hfsig->GetYaxis()->SetTitleOffset(1.4);
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
        TF1 *fitBW = new TF1("fitBW", gluefit3bW, 1.1, 2.3, 12);
        // fitBW->SetParNames("Norm1", "mass1", "width1", "Norm2", "mass2", "width2", "norm3", "mass3", "width3", "p1", "p2");
        // fitBW->SetParName(11, "p3");
        float f1270norm = parameter0(f1270Mass, f1270Width);
        float a1320norm = parameter0(a1320Mass, a1320Width);
        float f1525norm = parameter0(f1525Mass, f1525Width);
        float f1710norm = parameter0(f1710Mass, f1710Width);

        fitBW->SetParameter(0, f1270norm);
        fitBW->SetParameter(1, f1270Mass);
        fitBW->SetParameter(2, f1270Width);
        fitBW->SetParameter(3, f1525norm);
        fitBW->SetParameter(4, f1525Mass);
        fitBW->SetParameter(5, f1525Width);
        fitBW->SetParameter(6, f1710norm);
        fitBW->SetParameter(7, f1710Mass);
        fitBW->SetParameter(8, f1710Width);

        fitBW->SetParLimits(1, f1270Mass - 1 * f1270Width, f1270Mass + 1 * f1270Width);
        fitBW->SetParLimits(4, f1525Mass - 2 * f1525Width, f1525Mass + 2 * f1525Width);
        fitBW->SetParLimits(7, f1710Mass - 1 * f1710Width, f1710Mass + 1 * f1710Width);
        fitBW->FixParameter(1, f1270Mass);

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

        c1->SaveAs((koutputfolder + "/hglueball_bkg." + koutputtype).c_str());

        TCanvas *c2 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c2, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(fHistTotal[ip]);
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

        fHistTotal[ip]->SetMarkerStyle(8);
        fHistTotal[ip]->SetMarkerColor(kBlack);
        hfbkg->SetMarkerStyle(8);
        hfbkg->SetMarkerColor(kRed);
        fHistTotal[ip]->GetYaxis()->SetMaxDigits(3);
        fHistTotal[ip]->GetYaxis()->SetTitleOffset(1.4);
        fHistTotal[ip]->Draw("E");
        fHistTotal[ip]->GetYaxis()->SetTitle("Counts");
        hfbkg->Draw("E same");
        if (kResBkg == "MIX")
            hbkg_nopeak->Draw("BAR same");

        c2->SaveAs((koutputfolder + "/hglueball_invmass." + koutputtype).c_str());
    } // pt bin loop end here

    ////////////////////////////////////////////////////////////////////////
}
