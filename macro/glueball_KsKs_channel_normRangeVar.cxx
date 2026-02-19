#include <iostream>
#include <cmath>
#include <TKey.h>
#include <TClass.h>
#include <TDirectory.h>
#include "TArrow.h"
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/common_glue.h"
#include "src/fitting_range_glue.h"

using namespace std;

// void printDirectoryContents(TDirectory *dir, int indent = 0);
float parameter0(float mass, float width);

void glueball_KsKs_channel()
{
    // change here ***********************************************************
    const string kResBkg = "MIX";
    // const string kResBkg = "ROTATED";
    const bool makeQAplots = false;
    const bool calculate_inv_mass = true;
    const bool save_invmass_distributions = true;
    const bool save_multiPanel_plots = false;
    // change here ***********************************************************

    std::vector<std::string> variations = {
        "", "_DCA0p1", "_TPCPID2", "_TPCPID5", "_TPCMinCls100", "_TPCMinCls60", "_DCAv0dau0p3", "_DCAv0dau1p0", "_Ks_selection2p5", "_Ks_selection5", "_cospa0p95", "_cospa0p992", "_decay_rad1p0", "_lambda_rej4", "_lambda_rej6", "_lifetime15", "_lifetime25"}; // All variations

    int totalVar = variations.size();

    for (int ivar = 0; ivar < 1; ivar++)
    {

        string foldername_final = kfoldername + variations[ivar]; // default

        TString outputfolder = kSignalOutput + "/" + kchannel + "/" + foldername_final;
        TString outputQAfolder = kSignalOutput + "/" + kchannel + "/" + foldername_final + "/QA";
        const string outputfolder_str = kSignalOutput + "/" + kchannel + "/" + foldername_final;
        const string outputQAfolder_str = kSignalOutput + "/" + kchannel + "/" + foldername_final + "/QA";
        // Create the folder using TSystem::mkdir()
        if (gSystem->mkdir(outputfolder, kTRUE))
        {
            std::cout << "Folder " << outputfolder << " created successfully." << std::endl;
        }
        if (gSystem->mkdir(outputQAfolder, kTRUE))
        {
            std::cout << "Folder " << outputQAfolder << " created successfully." << std::endl;
        }
        // Folder name inside the Analysis.root file *****************************************
        if (!save_invmass_distributions)
            gStyle->SetOptFit(1111);
        // gStyle->SetOptStat(1110);
        gStyle->SetOptStat(0);

        t2->SetNDC(); // to self adjust the text so that it remains in the box
        t2->SetTextSize(0.045);
        t2->SetTextFont(42);

        // Input file
        TFile *fInputFile = new TFile(kDataFilename.c_str(), "Read");
        if (fInputFile->IsZombie())
        {
            cerr << "File not found " << endl;
            return;
        }

        // showing all folders in the root file using keys
        TIter next(fInputFile->GetListOfKeys());
        TKey *key;
        cout << "The folders in the root file are: \n";
        while ((key = (TKey *)next()))
        {
            cout << key->GetName() << endl;
        }
        // showing all the folders in the root file as well as their contents
        // printDirectoryContents(fInputFile);
        // if (save_invmass_distributions)
        // {
        //     // TFile *fileInvDistPair = new TFile((outputfolder_str + "/hglue_" + kResBkg + "cosTheta" + ".root").c_str(), "RECREATE");
        // }
        // TFile *fileInvDistPair = new TFile((outputfolder_str + "/hglue_" + kResBkg + Form("_pT_%.1f_%.1f.root", pT_bins[0], pT_bins[1])).c_str(), "RECREATE");
        TFile *fileInvDistPair = new TFile((outputfolder_str + "/hglue_" + kResBkg + "_allPt.root").c_str(), "RECREATE");
        // TFile *fileInvDistPair = new TFile((outputfolder_str + "/hglue_" + kResBkg + "_allPt_normright.root").c_str(), "RECREATE");

        TH1F *hmult = (TH1F *)fInputFile->Get((foldername_final + "/eventSelection/hmultiplicity").c_str());
        hmult->Write("multiplicity_histogram");
        if (hmult == nullptr)
        {
            cout << "Multiplicity histogram not found" << endl;
            return;
        }
        int multlow = 0;
        int multhigh = 100;
        double realevents = hmult->Integral(hmult->GetXaxis()->FindBin(multlow), hmult->GetXaxis()->FindBin(multhigh));
        cout << "*******number of events from the multiplicity histogram is *******:" << realevents << endl;

        // const string tempfoldername = "higher-mass-resonances_id25081"; // for common rotational background

        if (calculate_inv_mass)
        {
            // TH1F *hentries = (TH1F *)fInputFile->Get("event-selection-task/hColCounterAcc");
            // double Event = hentries->GetEntries();
            // cout << "*******number of events from the event selection histogram is *******:" << Event << endl;

            //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************

            THnSparseF *fHistNum = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassDS", foldername_final.c_str()));
            THnSparseF *fHistME = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassME", foldername_final.c_str()));
            THnSparseF *fHistRot = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassRot", foldername_final.c_str()));

            if (fHistNum == nullptr || fHistME == nullptr || fHistRot == nullptr)
            {
                cout << "Invariant mass histograms not found" << endl;
                return;
            }
            cout << " The number of entries in histograms: \n"
                 << "same event: " << fHistNum->GetEntries() << "\n"
                 << "mixed event: " << fHistME->GetEntries() << "\n"
                 << "rotated bkg/2: " << fHistRot->GetEntries() / 2 << endl;

            TH1D *fHistTotal[Npt];
            TH1D *fHistBkg[Npt];
            TH1D *fHistRotated[Npt];

            TCanvas *cbkgall1;
            TCanvas *cbkgall2;
            TCanvas *csigall1;
            TCanvas *csigall2;

            if (save_multiPanel_plots)
            {
                cbkgall1 = new TCanvas("", "all_bins", 1440, 720);
                cbkgall2 = new TCanvas("", "all_bins", 1440, 720);
                csigall1 = new TCanvas("", "all_bins", 1440, 720);
                csigall2 = new TCanvas("", "all_bins", 1440, 720);
                SetCanvasStyle(cbkgall1, 0.15, 0.03, 0.05, 0.15);
                SetCanvasStyle(cbkgall2, 0.15, 0.03, 0.05, 0.15);
                SetCanvasStyle(csigall1, 0.15, 0.03, 0.05, 0.15);
                SetCanvasStyle(csigall2, 0.15, 0.03, 0.05, 0.15);
                cbkgall1->Divide(3, 2);
                cbkgall2->Divide(2, 2);
                csigall1->Divide(3, 2);
                csigall2->Divide(2, 2);
            }

            TH1D *hbkg_temp[Npt];
            TH1D *hbkg_nopeak_temp[Npt];
            TH1D *hsig_temp[Npt];

            // TFile *fileInvDistPair;
            // if (Npt == 1)
            // {
            //     fileInvDistPair = new TFile((outputfolder_str + "/hglue_" + kResBkg + Form("_norm_%.2f_%.2f_pt_%.2f_%.2f", kNormRangepT[0][0], kNormRangepT[0][1], pT_bins[0], pT_bins[1]) + ".root").c_str(), "RECREATE");
            // }
            // else
            // {
            //     fileInvDistPair = new TFile((outputfolder_str + "/hglue_" + kResBkg + Form("_norm_%.2f_%.2f_all_pT", kNormRangepT[0][0], kNormRangepT[0][1]) + ".root").c_str(), "RECREATE");
            // }

            /*

            // double cosThetaBins[] = {-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0}; // cos(theta) bins
            double cosThetaBins[] = {-1, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 1.0}; // cos(theta) bins
            // double cosThetaBins[] = {-0.6, 0.6};// cos(theta) bins
            double numberOfCosThetaBins = sizeof(cosThetaBins) / sizeof(cosThetaBins[0]) - 1; // number of cos(theta) bins

            for (int ip = 0; ip < numberOfCosThetaBins; ip++)
            {
                float lowTheta = cosThetaBins[ip];
                float highTheta = cosThetaBins[ip + 1];
                float lowpT = pT_bins[0];
                float highpT = pT_bins[1];
                cout << "low theta value is " << lowTheta << " high theta value is " << highTheta << endl;

                int lbin = fHistNum->GetAxis(3)->FindBin(lowTheta + 1e-3);
                int hbin = fHistNum->GetAxis(3)->FindBin(highTheta - 1e-3);
                int lbinpT = fHistNum->GetAxis(1)->FindBin(lowpT + 1e-3);
                int hbinpT = fHistNum->GetAxis(1)->FindBin(highpT - 1e-3);

                cout << "the value of lpt is " << fHistNum->GetAxis(1)->GetBinLowEdge(lbinpT) << endl;
                cout << "the value of hpt is " << fHistNum->GetAxis(1)->GetBinUpEdge(hbinpT) << endl;

                fHistNum->GetAxis(3)->SetRange(lbin, hbin);
                fHistME->GetAxis(3)->SetRange(lbin, hbin);
                fHistRot->GetAxis(3)->SetRange(lbin, hbin);

                fHistNum->GetAxis(1)->SetRange(lbinpT, hbinpT);
                fHistME->GetAxis(1)->SetRange(lbinpT, hbinpT);
                fHistRot->GetAxis(1)->SetRange(lbinpT, hbinpT);

                int lbinmult = fHistNum->GetAxis(0)->FindBin(multlow + 1e-5);
                int hbinmult = fHistNum->GetAxis(0)->FindBin(multhigh - 1e-5);

                fHistNum->GetAxis(0)->SetRange(lbinmult, hbinmult);
                fHistME->GetAxis(0)->SetRange(lbinmult, hbinmult);
                fHistRot->GetAxis(0)->SetRange(lbinmult, hbinmult);

                fHistTotal[ip] = fHistNum->Projection(2, "E");
                fHistBkg[ip] = fHistME->Projection(2, "E");
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

                if (kResBkg == "MIX" || kResBkg == "ROTATED")
                {
                    auto sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(normRangeLow), fHistTotal[ip]->GetXaxis()->FindBin(normRangeHigh)));

                    auto bkg_integral = (fHistBkg[ip]->Integral(fHistBkg[ip]->GetXaxis()->FindBin(normRangeLow), fHistBkg[ip]->GetXaxis()->FindBin(normRangeHigh)));

                    auto bkg_integral_rotated = (fHistRotated[ip]->Integral(fHistRotated[ip]->GetXaxis()->FindBin(normRangeLow), fHistRotated[ip]->GetXaxis()->FindBin(normRangeHigh)));

                    auto normfactor = sigbkg_integral / bkg_integral;                 // scaling factor for mixed bkg
                    auto normfactor_rotated = sigbkg_integral / bkg_integral_rotated; // scaling factor for rotated bkg
                    cout << "\n\n normalization factor " << 1. / normfactor << "\n\n";
                    if (kResBkg == "MIX")
                    {
                        hfbkg = (TH1D *)fHistBkg[ip]->Clone();
                        hfbkg->Write("bkg_without_normalization");
                        hfbkg->Scale(normfactor);
                    }
                    else
                    {
                        hfbkg = (TH1D *)fHistRotated[ip]->Clone();
                        hfbkg->Write("bkg_without_normalization");
                        hfbkg->Scale(normfactor_rotated);
                    }

                    hfbkg->Rebin(kRebin[ip]);
                    hfsig->Rebin(kRebin[ip]);

                    hfsig->Add(hfbkg, -1);
                }
                // else if (kResBkg == "ROTATED")
                // {
                //     hfbkg = (TH1D *)fHistRotated[ip]->Clone();
                //     hfbkg->Scale(0.5);
                //     hfbkg->Rebin(kRebin[ip]);
                //     hfsig->Rebin(kRebin[ip]);
                //     hfsig->Add(hfbkg, -1);
                // }

                fHistTotal[ip]->Rebin(kRebin[ip]);

                //*****************************************************************************************************
                TCanvas *c1 = new TCanvas("", "", 720, 720);
                SetCanvasStyle(c1, 0.15, 0.015, 0.05, 0.155);
                SetHistoQA(hfsig);
                hfsig->SetTitle(0);
                hfsig->SetMarkerStyle(20);
                hfsig->SetMarkerSize(1.1);
                hfsig->GetYaxis()->SetMaxDigits(3);
                hfsig->GetYaxis()->SetTitleOffset(1.5);
                hfsig->SetMarkerColor(kBlack);
                hfsig->SetLineColor(kBlack);
                hfsig->GetXaxis()->SetTitle("M_{K^{0}_{s}K^{0}_{s}} (GeV/c^{2})");
                hfsig->GetYaxis()->SetTitle(Form("Counts/%.0f MeV/c^{2}", binwidth_file * 1000));
                hfsig->GetXaxis()->SetRangeUser(1.00, 2.45);
                hfsig->Draw("e");
                TLine *linesig = new TLine(1.0, 0, 2.50, 0);
                linesig->SetLineColor(kRed);
                linesig->SetLineStyle(2);
                linesig->SetLineWidth(2);
                // linesig->Draw("same");
                // t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
                hfsig->Write(Form("ksks_subtracted_invmass_theta_%.1f_%.1f", lowTheta, highTheta));
                // gPad->Update();
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
                TLegend *lp2 = DrawLegend(0.55, 0.63, 0.85, 0.91);
                lp2->SetTextSize(0.035);
                lp2->SetTextFont(42);
                lp2->SetFillStyle(0);
                lp2->AddEntry((TObject *)0, "ALICE Performance", "");
                lp2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
                lp2->AddEntry((TObject *)0, "FT0M, 0-100%", "");
                lp2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
                lp2->AddEntry((TObject *)0, Form("%.1f < cos(#theta) < %.1f", lowTheta, highTheta), "");
                lp2->Draw("same");

                if (save_invmass_distributions)
                {
                    c1->SaveAs((outputfolder_str + "/hglueball_signal_" + kResBkg + Form("cosTheta_%.1f_%.1f_norm_%.1f_%.1f.", lowTheta, highTheta, normRangeLow, normRangeHigh) + koutputtype).c_str());
                }
                if (save_multiPanel_plots)
                {
                    if (ip < 6)
                        csigall1->cd(ip + 1);
                    else
                        csigall2->cd(ip - 5);

                    gPad->SetTopMargin(0.05);
                    gPad->SetRightMargin(0.03);
                    gPad->SetLeftMargin(0.12);
                    gPad->SetBottomMargin(0.12);
                    TH1F *hfsig_clone = (TH1F *)hfsig->Clone();
                    hfsig_clone->GetXaxis()->SetRangeUser(1.00, 2.40);
                    hfsig_clone->SetMarkerSize(0.8);
                    hfsig_clone->GetYaxis()->SetTitleOffset(1.0);
                    hfsig_clone->Draw("e");
                    // lp2->Draw("same");
                    t2->SetTextSize(0.05);
                    if (ip == 0)
                        t2->DrawLatex(0.27, 0.83, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[0], pT_bins[1]));
                    t2->DrawLatex(0.27, 0.73, Form("#bf{%.1f < cos#Theta < %.1f}", lowTheta, highTheta));
                }

                TCanvas *c2 = new TCanvas("", "", 720, 720);
                SetCanvasStyle(c2, 0.15, 0.01, 0.05, 0.135);
                SetHistoQA(fHistTotal[ip]);
                SetHistoQA(hfbkg);

                TH1F *hbkg_nopeak = (TH1F *)hfbkg->Clone();
                hbkg_nopeak->SetLineColor(kBlue - 7);
                hbkg_nopeak->SetMarkerColor(kBlue - 7);
                hbkg_nopeak->SetFillColor(kBlue - 7);
                // hbkg_nopeak->SetFillStyle(3001);
                for (int i = 0; i < hbkg_nopeak->GetNbinsX(); i++)
                {
                    if (hbkg_nopeak->GetBinCenter(i + 1) < normRangeLow || hbkg_nopeak->GetBinCenter(i + 1) > normRangeHigh)
                    {
                        hbkg_nopeak->SetBinContent(i + 1, -999);
                    }
                }

                fHistTotal[ip]->SetMarkerStyle(20);
                fHistTotal[ip]->SetMarkerColor(kBlack);
                fHistTotal[ip]->SetMarkerSize(1.1);
                hfbkg->SetMarkerStyle(20);
                hfbkg->SetMarkerSize(1.1);
                hfbkg->SetMarkerColor(kRed);
                hfbkg->SetLineColor(kRed);
                fHistTotal[ip]->GetYaxis()->SetMaxDigits(3);
                fHistTotal[ip]->GetYaxis()->SetTitleOffset(1.5);
                fHistTotal[ip]->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidth_file * 1000));
                fHistTotal[ip]->SetMaximum(1.4 * fHistTotal[ip]->GetMaximum());
                fHistTotal[ip]->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
                fHistTotal[ip]->Draw("E");
                fHistTotal[ip]->Write(Form("ksks_invmass_cosTheta_%.1f_%.1f", lowTheta, highTheta));
                hfbkg->Write(Form("ksks_bkg_cosTheta_%.1f_%.1f", lowTheta, highTheta));
                // if (save_invmass_distributions)
                // {
                //     c2->SaveAs((outputfolder_str + "/hglueball_invmass_only_." + Form("cosTheta_%.1f_%.1f_.", lowTheta, highTheta) + koutputtype).c_str());
                // }
                hfbkg->Draw("E same");
                if (kResBkg == "MIX" || kResBkg == "ROTATED")
                    hbkg_nopeak->Draw("BAR same");

                TLegend *leg = new TLegend(0.25, 0.2454598, 0.5445682, 0.3908046);
                leg->SetFillStyle(0);
                leg->SetBorderSize(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.035);
                leg->AddEntry(fHistTotal[ip], "Same-event pairs", "p");
                string bkgname = (kResBkg == "MIX") ? "Mixed event" : "Same-event rotated paris";
                leg->AddEntry(hfbkg, bkgname.c_str(), "p");
                // if (kResBkg == "MIX")
                hbkg_nopeak->SetLineWidth(0);
                leg->AddEntry(hbkg_nopeak, "Normalization region", "f");
                leg->Draw();
                lp2->Draw("same");

                // t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
                if (save_invmass_distributions)
                {
                    c2->SaveAs((outputfolder_str + "/hglueball_invmass_" + kResBkg + Form("cosTheta_%.1f_%.1f_norm_%.1f_%.1f.", lowTheta, highTheta, normRangeLow, normRangeHigh) + koutputtype).c_str());
                }
                c2->Write(Form("ksks_invmass_withbkg_cosTheta_%.1f_%.1f", lowTheta, highTheta));

                if (save_multiPanel_plots)
                {
                    if (ip < 6)
                        cbkgall1->cd(ip + 1);
                    else
                        cbkgall2->cd(ip - 5);

                    fHistTotal[ip]->Draw("E");
                    hfbkg->Draw("E same");
                    if (kResBkg == "MIX" || kResBkg == "ROTATED")
                        hbkg_nopeak->Draw("BAR same");
                    leg->Draw();
                    t2->DrawLatex(0.27, 0.8, Form("#bf{%.1f < cos#Theta < %.1f}", lowTheta, highTheta));

                    hbkg_temp[ip] = (TH1D *)hfbkg->Clone();
                    hbkg_nopeak_temp[ip] = (TH1D *)hbkg_nopeak->Clone();
                    hsig_temp[ip] = (TH1D *)hfsig->Clone();
                }
            } // cos theta loop ended

            if (save_multiPanel_plots)
            {
                csigall1->SaveAs((outputfolder_str + "/MultBins_glueballSignal_cosTheta1" + kResBkg + Form("_pT_%.0f-%.0f", pT_bins[0], pT_bins[1]) + ".png").c_str());
                csigall2->SaveAs((outputfolder_str + "/MultBins_glueballSignal_cosTheta2" + kResBkg + Form("_pT_%.0f-%.0f", pT_bins[0], pT_bins[1]) + ".png").c_str());
                cbkgall1->SaveAs((outputfolder_str + "/MultBins_glueballBackground_cosTheta1" + kResBkg + Form("_pT_%.0f-%.0f", pT_bins[0], pT_bins[1]) + ".png").c_str());
                cbkgall2->SaveAs((outputfolder_str + "/MultBins_glueballBackground_cosTheta2" + kResBkg + Form("_pT_%.0f-%.0f", pT_bins[0], pT_bins[1]) + ".png").c_str());
            }

            */

            // float pt_binsTemp[] = {0.0, 1.0, 2.0, 3.0};
            TCanvas *c1divide = new TCanvas("", "all_bins", 1080, 720);
            SetCanvasStyle(c1divide, 0.15, 0.03, 0.05, 0.15);
            c1divide->Divide(2, 2);
            TCanvas *c2divide = new TCanvas("", "all_bins", 1080, 720);
            SetCanvasStyle(c2divide, 0.15, 0.03, 0.05, 0.15);
            c2divide->Divide(2, 2);

            // /*

            // for (int imult = 0; imult < Nmult + 1; imult++)
            for (int imult = 0; imult < 1; imult++)
            {

                if (imult == 0)
                {
                    multlow = 0;
                    multhigh = 100; // for all multiplicity
                }
                else
                {
                    multlow = mult_classes[imult - 1];
                    multhigh = mult_classes[imult];
                }
                fileInvDistPair->cd();
                TDirectory *dir = fileInvDistPair->mkdir(Form("multiplicity_%.0f_%.0f", (float)multlow, (float)multhigh));
                dir->cd();

                for (Int_t ip = pt_start; ip < 1; ip++) // start pt bin loop
                {
                    float normRangesLow[] = {1.0, 1.7, 2.15, 2.5};
                    float normRangesHigh[] = {1.05, 1.8, 2.25, 2.6};
                    for (int inorm = 0; inorm < 4; inorm++)
                    {
                        float lowpt = pT_bins[ip];
                        float highpt = pT_bins[ip + 1];
                        // float lowpt = pt_binsTemp[ip];
                        // float highpt = 30.0;
                        float normRangeLow = normRangesLow[inorm];
                        float normRangeHigh = normRangesHigh[inorm];

                        cout << "low pt value is " << lowpt << " high pt value is " << highpt << endl;
                        int lbin = fHistNum->GetAxis(1)->FindBin(lowpt + 1e-3);
                        int hbin = fHistNum->GetAxis(1)->FindBin(highpt - 1e-3);

                        fHistNum->GetAxis(1)->SetRange(lbin, hbin);
                        fHistME->GetAxis(1)->SetRange(lbin, hbin);
                        fHistRot->GetAxis(1)->SetRange(lbin, hbin);

                        int lbinmult = fHistNum->GetAxis(0)->FindBin(multlow + 1e-3);
                        int hbinmult = fHistNum->GetAxis(0)->FindBin(multhigh - 1e-3);

                        fHistNum->GetAxis(0)->SetRange(lbinmult, hbinmult);
                        fHistME->GetAxis(0)->SetRange(lbinmult, hbinmult);
                        fHistRot->GetAxis(0)->SetRange(lbinmult, hbinmult);

                        fHistTotal[ip] = fHistNum->Projection(2, "E");
                        fHistBkg[ip] = fHistME->Projection(2, "E");
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

                        if (kResBkg == "MIX" || kResBkg == "ROTATED")
                        {
                            auto sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(normRangeLow), fHistTotal[ip]->GetXaxis()->FindBin(normRangeHigh)));
                            auto bkg_integral = (fHistBkg[ip]->Integral(fHistBkg[ip]->GetXaxis()->FindBin(normRangeLow), fHistBkg[ip]->GetXaxis()->FindBin(normRangeHigh)));
                            auto bkg_integral_rotated = (fHistRotated[ip]->Integral(fHistRotated[ip]->GetXaxis()->FindBin(normRangeLow), fHistRotated[ip]->GetXaxis()->FindBin(normRangeHigh)));
                            auto normfactor = sigbkg_integral / bkg_integral;                 // scaling factor for mixed bkg
                            auto normfactor_rotated = sigbkg_integral / bkg_integral_rotated; // scaling factor for rotated bkg
                            cout << "\n\n normalization factor " << 1. / normfactor << "\n\n";
                            if (kResBkg == "MIX")
                            {
                                hfbkg = (TH1D *)fHistBkg[ip]->Clone();
                                hfbkg->Write("bkg_without_normalization");
                                hfbkg->Scale(normfactor);
                            }
                            else
                            {
                                hfbkg = (TH1D *)fHistRotated[ip]->Clone();
                                hfbkg->Write("bkg_without_normalization");
                                hfbkg->Scale(normfactor_rotated);
                            }

                            hfbkg->Rebin(kRebin[ip]);
                            hfsig->Rebin(kRebin[ip]);

                            hfsig->Add(hfbkg, -1);
                        }
                        // else if (kResBkg == "ROTATED")
                        // {
                        //     hfbkg = (TH1D *)fHistRotated[ip]->Clone();
                        //     hfbkg->Scale(0.5);
                        //     hfbkg->Rebin(kRebin[ip]);
                        //     hfsig->Rebin(kRebin[ip]);
                        //     hfsig->Add(hfbkg, -1);
                        // }

                        fHistTotal[ip]->Rebin(kRebin[ip]);

                        //*****************************************************************************************************
                        // TCanvas *c1 = new TCanvas("", "", 720, 720);
                        // SetCanvasStyle(c1, 0.15, 0.015, 0.05, 0.155);
                        c1divide->cd(inorm + 1);
                        gPad->SetBottomMargin(0.15);
                        gPad->SetLeftMargin(0.15);
                        gPad->SetRightMargin(0.05);
                        gPad->SetTopMargin(0.05);
                        SetHistoQA(hfsig);
                        hfsig->SetTitle(0);
                        hfsig->SetMarkerStyle(20);
                        hfsig->SetMarkerSize(0.8);
                        hfsig->GetYaxis()->SetMaxDigits(3);
                        hfsig->GetYaxis()->SetTitleOffset(1.5);
                        hfsig->SetMarkerColor(kBlack);
                        hfsig->SetLineColor(kBlack);
                        hfsig->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
                        hfsig->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidth_file * 1000));
                        hfsig->GetXaxis()->SetRangeUser(1.00, 2.50);
                        hfsig->Draw("e");
                        TLine *linesig = new TLine(1.0, 0, 2.50, 0);
                        linesig->SetLineColor(kRed);
                        linesig->SetLineStyle(2);
                        linesig->SetLineWidth(2);
                        // linesig->Draw("same");
                        // t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
                        hfsig->Write(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", lowpt, highpt));
                        // gPad->Update();
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
                        // TLegend *lp2 = DrawLegend(0.55, 0.58, 0.85, 0.89);
                        TLegend *lp2 = DrawLegend(0.55, 0.65, 0.85, 0.92);
                        lp2->SetTextSize(0.04);
                        // lp2->SetTextSize(0.055);
                        lp2->SetTextFont(42);
                        lp2->SetFillStyle(0);
                        lp2->AddEntry((TObject *)0, "ALICE", "");
                        lp2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
                        lp2->AddEntry((TObject *)0, "FT0M, 0-100%", "");
                        lp2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
                        // lp2->AddEntry((TObject *)0, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", lowpt, highpt), "");
                        // if (inorm == 0)
                        //     lp2->Draw("same");

                        t2->DrawLatex(0.6, 0.88, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
                        t2->DrawLatex(0.5, 0.8, Form("#bf{Norm range: %.2f - %.2f GeV/#it{c}^{2}}", normRangeLow, normRangeHigh));
                        TLine *lineat0 = new TLine(1.0, 0, 2.50, 0);
                        lineat0->SetLineColor(kRed);
                        lineat0->SetLineStyle(2);
                        lineat0->SetLineWidth(2);
                        lineat0->Draw("same");

                        // if (save_invmass_distributions)
                        // {
                        //     c1->SaveAs((outputfolder_str + "/hglueball_signal_" + kResBkg + Form("pT_%.1f_%.1f_norm_%.2f_%.2f.", lowpt, highpt, normRangeLow, normRangeHigh) + koutputtype).c_str());
                        // }

                        // TCanvas *c2 = new TCanvas("", "", 720, 720);
                        // SetCanvasStyle(c2, 0.15, 0.01, 0.05, 0.135);
                        c2divide->cd(inorm + 1);
                        gPad->SetBottomMargin(0.15);
                        gPad->SetLeftMargin(0.15);
                        gPad->SetRightMargin(0.05);
                        gPad->SetTopMargin(0.05);
                        SetHistoQA(fHistTotal[ip]);
                        SetHistoQA(hfbkg);

                        TH1F *hbkg_nopeak = (TH1F *)hfbkg->Clone();
                        hbkg_nopeak->SetLineColor(kBlue - 7);
                        hbkg_nopeak->SetMarkerColor(kBlue - 7);
                        hbkg_nopeak->SetFillColor(kBlue - 7);
                        // hbkg_nopeak->SetFillStyle(3001);
                        for (int i = 0; i < hbkg_nopeak->GetNbinsX(); i++)
                        {
                            if (hbkg_nopeak->GetBinCenter(i + 1) < normRangeLow || hbkg_nopeak->GetBinCenter(i + 1) > normRangeHigh)
                            {
                                hbkg_nopeak->SetBinContent(i + 1, -999);
                            }
                        }

                        fHistTotal[ip]->SetMarkerStyle(20);
                        fHistTotal[ip]->SetMarkerColor(kBlack);
                        fHistTotal[ip]->SetMarkerSize(0.8);
                        hfbkg->SetMarkerStyle(20);
                        hfbkg->SetMarkerSize(0.8);
                        hfbkg->SetMarkerColor(kRed);
                        hfbkg->SetLineColor(kRed);
                        fHistTotal[ip]->GetYaxis()->SetMaxDigits(3);
                        fHistTotal[ip]->GetYaxis()->SetTitleOffset(1.5);
                        fHistTotal[ip]->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidth_file * 1000));
                        // fHistTotal[ip]->SetMaximum(1.2 * fHistTotal[ip]->GetMaximum());
                        fHistTotal[ip]->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
                        if (inorm == 0)
                            fHistTotal[ip]->SetMaximum(1.4 * fHistTotal[ip]->GetMaximum());
                        fHistTotal[ip]->Draw("E");
                        fHistTotal[ip]->Write(Form("ksks_invmass_pt_%.1f_%.1f", lowpt, highpt));
                        hfbkg->Write(Form("ksks_bkg_pt_%.1f_%.1f", lowpt, highpt));
                        // if (save_invmass_distributions)
                        // {
                        //     c2->SaveAs((outputfolder_str + "/hglueball_invmass_only_." + Form("pT_%.1f_%.1f_.", lowpt, highpt) + koutputtype).c_str());
                        // }
                        hfbkg->Draw("E same");
                        if (kResBkg == "MIX" || kResBkg == "ROTATED")
                            hbkg_nopeak->Draw("BAR same");

                        TLegend *leg = new TLegend(0.25, 0.2454598, 0.5445682, 0.3908046);
                        leg->SetFillStyle(0);
                        leg->SetBorderSize(0);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.035);
                        // leg->SetTextSize(0.055);
                        leg->AddEntry(fHistTotal[ip], "Same-event pairs", "p");
                        string bkgname = (kResBkg == "MIX") ? "Mixed event" : "Same-event rotated paris";
                        leg->AddEntry(hfbkg, bkgname.c_str(), "p");
                        // if (kResBkg == "MIX")
                        hbkg_nopeak->SetLineWidth(0);
                        leg->AddEntry(hbkg_nopeak, "Normalization region", "f");
                        if (inorm == 0)
                        {
                            leg->Draw();
                            // lp2->Draw("same");
                        }
                        // t2->DrawLatex(0.6, 0.6, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
                        // t2->DrawLatex(0.5, 0.6, Form("#bf{Norm range: %.2f - %.2f GeV/#it{c}^{2}}", normRangeLow, normRangeHigh));
                        t2->DrawLatex(0.6, 0.88, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
                        t2->DrawLatex(0.5, 0.8, Form("#bf{Norm range: %.2f - %.2f GeV/#it{c}^{2}}", normRangeLow, normRangeHigh));

                        // // t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
                        // if (save_invmass_distributions)
                        // {
                        //     c2->SaveAs((outputfolder_str + "/hglueball_invmass_" + kResBkg + Form("pT_%.1f_%.1f_norm_%.2f_%.2f.", lowpt, highpt, normRangeLow, normRangeHigh) + koutputtype).c_str());
                        // }
                        // c2->Write(Form("ksks_invmass_withbkg_pt_%.1f_%.1f", lowpt, highpt));

                        // cbkgall1->cd(ip + 1);
                        // fHistTotal[ip]->Draw("E");
                        // hfbkg->Draw("E same");
                        // if (kResBkg == "MIX" || kResBkg == "ROTATED")
                        //     hbkg_nopeak->Draw("BAR same");
                        // leg->Draw();
                        // t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
                        hbkg_temp[ip] = (TH1D *)hfbkg->Clone();
                        hbkg_nopeak_temp[ip] = (TH1D *)hbkg_nopeak->Clone();
                        hsig_temp[ip] = (TH1D *)hfsig->Clone();
                    }
                    c1divide->SaveAs((outputfolder_str + "/hglueball_signal_norm" + kResBkg + "." + koutputtype).c_str());
                    c2divide->SaveAs((outputfolder_str + "/hglueball_invmass_norm" + kResBkg + "." + koutputtype).c_str());
                } // pt bin loop end here
            } // end of multiplicity loop

            // */

            // TCanvas *cbkg = new TCanvas("", "", 1080, 720);
            // cbkg->Divide(2, 2);
            // SetCanvasStyle(cbkg, 0.18, 0.03, 0.03, 0.15);

            // for (Int_t ip = pt_start; ip < pt_end - 1; ip++) // start pt bin loop
            // {
            //     cbkg->cd(ip + 1);
            //     gPad->SetLeftMargin(0.15);
            //     gPad->SetRightMargin(0.03);
            //     gPad->SetTopMargin(0.05);
            //     gPad->SetBottomMargin(0.15);
            //     fHistTotal[ip]->SetMaximum(1.1 * fHistTotal[ip]->GetMaximum());
            //     fHistTotal[ip]->Draw("E");
            //     hbkg_temp[ip]->Draw("E same");
            //     if (kResBkg == "MIX" || kResBkg == "ROTATED")
            //         hbkg_nopeak_temp[ip]->Draw("BAR same");
            //     TLegend *leg = new TLegend(0.1851253, 0.2454598, 0.5445682, 0.3908046);
            //     leg->SetFillStyle(0);
            //     leg->SetBorderSize(0);
            //     leg->SetTextFont(42);
            //     leg->SetTextSize(0.055);
            //     leg->AddEntry(fHistTotal[ip], "Same event K^{0}_{s}K^{0}_{s} pair", "lpe");
            //     leg->AddEntry(hbkg_temp[ip], "Mixed event K^{0}_{s}K^{0}_{s} pair", "lpe");
            //     leg->AddEntry(hbkg_nopeak_temp[ip], "Norm. region", "f");
            //     if (ip == 0)
            //     {
            //         leg->Draw();
            //     }

            //     TLegend *lp2 = DrawLegend(0.6, 0.56, 0.92, 0.91);
            //     lp2->SetTextSize(0.055);
            //     lp2->SetTextFont(42);
            //     lp2->SetFillStyle(0);
            //     lp2->SetBorderSize(0);
            //     if (ip == 0)
            //     {
            //         lp2->AddEntry((TObject *)0, "ALICE", "");
            //         lp2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
            //         lp2->AddEntry((TObject *)0, "FT0M, 0-100%", "");
            //         lp2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
            //     }
            //     lp2->AddEntry((TObject *)0, Form("#it{p}_{T}: %.0f - %.0f GeV/#it{c}", lowpt, highpt), "");

            //     lp2->Draw("same");

            //     // t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
            // }
            // cbkg->SaveAs((outputfolder_str + "/hglueball_invmass_allbins." + kResBkg + ".pdf").c_str());

            // TCanvas *csignal = new TCanvas("", "", 1080, 720);
            // csignal->Divide(2, 2);
            // SetCanvasStyle(csignal, 0.18, 0.03, 0.03, 0.15);

            // for (Int_t ip = pt_start; ip < pt_end - 1; ip++) // start pt bin loop
            // {
            //     csignal->cd(ip + 1);
            //     gPad->SetLeftMargin(0.15);
            //     gPad->SetRightMargin(0.03);
            //     gPad->SetTopMargin(0.05);
            //     gPad->SetBottomMargin(0.15);
            //     hsig_temp[ip]->SetMaximum(1.1 * hsig_temp[ip]->GetMaximum());
            //     hsig_temp[ip]->Draw("E");
            //     // draw a line at counts 0
            //     TLine *line = new TLine(1.0, 0, 2.5, 0);
            //     line->SetLineColor(kRed);
            //     line->SetLineWidth(2);
            //     line->SetLineStyle(2);
            //     line->Draw("same");

            //     TLegend *leg = new TLegend(0.1851253, 0.5454598, 0.5445682, 0.6408046);
            //     leg->SetFillStyle(0);
            //     leg->SetBorderSize(0);
            //     leg->SetTextFont(42);
            //     leg->SetTextSize(0.055);
            //     leg->AddEntry(hsig_temp[ip], "Same event K^{0}_{s}K^{0}_{s} pair", "lpe");
            //     if (ip == 0)
            //     {
            //         leg->Draw();
            //     }

            //     TLegend *lp2 = DrawLegend(0.2, 0.56, 0.62, 0.91);
            //     lp2->SetTextSize(0.055);
            //     lp2->SetTextFont(42);
            //     lp2->SetFillStyle(0);
            //     lp2->SetBorderSize(0);
            //     lp2->AddEntry((TObject *)0, Form("#it{p}_{T}: %.0f - %.0f GeV/#it{c}", lowpt, highpt), "");
            //     lp2->Draw("same");
            // }
            // csignal->SaveAs((outputfolder_str + "/hglueball_invmass_allbins_signal." + kResBkg + ".pdf").c_str());
            // cbkgall1->SaveAs((outputfolder_str + "/hglueball_invmass_all." + kResBkg + ".png").c_str());
        }

        ////====================================================================================////
        ////**************************************QA PLOTS*********************************** **////
        ////====================================================================================////

        if (makeQAplots)
        {
            TFile *KsInvMass = new TFile((outputQAfolder_str + "/KsInvMass.root").c_str(), "RECREATE");
            TCanvas *c3 = new TCanvas("", "", 720, 720);
            SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);

            // Kshort pT and invariant mass distribution before the selections
            THnSparseF *hKshortPt = (THnSparseF *)fInputFile->Get((foldername_final + "/kzeroShort/hMassK0Shortbefore").c_str());
            if (hKshortPt == nullptr)
            {
                cout << "Kshort pT distribution not found" << endl;
                return;
            }
            TH1F *kshortpt_before = (TH1F *)hKshortPt->Projection(1, "E");
            TH1F *kshortmass_before = (TH1F *)hKshortPt->Projection(0, "E");
            SetHistoQA(kshortpt_before);
            SetHistoQA(kshortmass_before);
            kshortpt_before->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            kshortpt_before->GetYaxis()->SetTitle("Counts");
            kshortmass_before->GetXaxis()->SetTitle("M_{K^{0}_{s}} (GeV/c^{2})");
            kshortmass_before->GetYaxis()->SetTitle("Counts");
            kshortmass_before->GetXaxis()->SetRangeUser(0.2, 0.8);
            kshortpt_before->Draw("HIST");
            c3->SaveAs((outputQAfolder_str + "/kshort_pt_before." + koutputtype).c_str());
            c3->Clear();
            kshortmass_before->Draw("HIST");
            c3->SaveAs((outputQAfolder_str + "/kshort_mass_before." + koutputtype).c_str());

            // Kshort pT and invariant mass distribution after the selections
            THnSparseF *hKshortPt_after = (THnSparseF *)fInputFile->Get((foldername_final + "/kzeroShort/hMassK0ShortSelected").c_str());
            if (hKshortPt_after == nullptr)
            {
                cout << "Kshort pT distribution after the selections not found" << endl;
                return;
            }
            hKshortPt_after->Write("kshort_2dsparse");
            TH1F *kshortpt_after = (TH1F *)hKshortPt_after->Projection(1, "E");
            TH1F *kshortmass_after = (TH1F *)hKshortPt_after->Projection(0, "E");
            SetHistoQA(kshortpt_after);
            SetHistoQA(kshortmass_after);
            kshortpt_after->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            kshortpt_after->GetYaxis()->SetTitle("Counts");
            kshortmass_after->GetXaxis()->SetTitle("M_{K^{0}_{s}} (GeV/c^{2})");
            kshortmass_after->GetYaxis()->SetTitle("Counts");
            kshortmass_after->GetXaxis()->SetRangeUser(0.2, 0.8);
            kshortpt_after->Draw("HIST");
            kshortpt_after->Write("kspt");
            kshortmass_after->Write("ksmass");
            // KsInvMass->Close();
            c3->SaveAs((outputQAfolder_str + "/kshort_pt_after." + koutputtype).c_str());
            c3->Clear();
            kshortmass_after->Draw("HIST");
            c3->SaveAs((outputQAfolder_str + "/kshort_mass_after." + koutputtype).c_str());

            TF1 *fitKshort = new TF1("fitKshort", "gaus", 0.45, 0.55);
            kshortmass_after->Fit("fitKshort", "RM");
            kshortmass_after->GetXaxis()->SetRangeUser(0.4, 0.6);
            kshortmass_after->Draw("HIST");
            fitKshort->Draw("same");
            c3->SaveAs((outputQAfolder_str + "/kshort_mass_after_fit." + koutputtype).c_str());

            // Mulitplicity plot
            SetHistoQA(hmult);
            hmult->GetYaxis()->SetTitle("Counts");
            hmult->GetXaxis()->SetTitle("Multiplicity percentile");
            hmult->GetXaxis()->SetRangeUser(0, 105);
            hmult->SetStats(0);
            hmult->Draw();
            c3->SaveAs((outputQAfolder_str + "/hglueball_multiplicity_percentile." + koutputtype).c_str());

            // vtz distribution plot
            TH1F *hvtz = (TH1F *)fInputFile->Get((foldername_final + "/eventSelection/hVertexZRec").c_str());
            if (hvtz == nullptr)
            {
                cout << "Vertex Z distribution not found" << endl;
                return;
            }
            c3->Clear();
            SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            SetHistoQA(hvtz);
            hvtz->GetYaxis()->SetTitle("Counts");
            hvtz->GetXaxis()->SetTitle("Vertex Z (cm)");
            hvtz->SetStats(0);
            hvtz->Draw();
            c3->SaveAs((outputQAfolder_str + "/hglueball_vtz." + koutputtype).c_str());

            // // mass correlation plot between two Ks
            // TH2F *hmasscorr = (TH2F *)fInputFile->Get((foldername_final + "/kzeroShort/hMasscorrelationbefore").c_str());
            // if (hmasscorr == nullptr)
            // {
            //     cout << "Mass correlation plot not found" << endl;
            //     return;
            // }
            // hmasscorr->Write("ksks_mass_correlation");
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.14, 0.05, 0.15);
            // SetHistoQA(hmasscorr);
            // hmasscorr->GetYaxis()->SetTitle("M_{K^{0}_{s}} (GeV/c^{2})");
            // hmasscorr->GetXaxis()->SetTitle("M_{K^{0}_{s}} (GeV/c^{2})");
            // hmasscorr->GetXaxis()->SetTitleOffset(3.0);
            // hmasscorr->GetYaxis()->SetTitleOffset(3.0);
            // hmasscorr->GetXaxis()->SetTitleSize(0.04);
            // hmasscorr->GetYaxis()->SetTitleSize(0.04);
            // hmasscorr->GetZaxis()->SetTitle(Form("Counts/%.0f MeV/#it{c}^{2}", hmasscorr->GetXaxis()->GetBinWidth(1) * 1000));
            // hmasscorr->GetZaxis()->SetTitleSize(0.04);
            // hmasscorr->GetZaxis()->SetTitleOffset(2.0);
            // hmasscorr->GetXaxis()->SetRangeUser(0.475, 0.52);
            // hmasscorr->GetYaxis()->SetRangeUser(0.475, 0.52);
            // hmasscorr->GetXaxis()->SetMaxDigits(3);
            // hmasscorr->GetYaxis()->SetMaxDigits(3);
            // hmasscorr->GetXaxis()->SetNdivisions(505);
            // hmasscorr->GetYaxis()->SetNdivisions(505);
            // hmasscorr->Draw("surf1");
            // c3->SaveAs((outputQAfolder_str + "/ksks_mass_correlation." + koutputtype).c_str());

            // // kshort selection plots
            // // Armenteros alpha plot
            // TH1F *hArmenteros = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/halpha").c_str());
            // if (hArmenteros == nullptr)
            // {
            //     cout << "Armenteros alpha plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hArmenteros);
            // hArmenteros->GetYaxis()->SetTitle("Counts");
            // hArmenteros->GetXaxis()->SetTitle("#alpha");
            // hArmenteros->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_alpha." + koutputtype).c_str());

            // // DCA negative daughter to PV
            // TH1F *hDCAneg = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/hDCAnegtopv").c_str());
            // if (hDCAneg == nullptr)
            // {
            //     cout << "DCA negative daughter to PV plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hDCAneg);
            // hDCAneg->GetYaxis()->SetTitle("Counts");
            // hDCAneg->GetXaxis()->SetTitle("DCA neg. daughter to PV (cm)");
            // hDCAneg->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_DCAnegtopv." + koutputtype).c_str());

            // // DCA positive daughter to PV
            // TH1F *hDCApos = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/hDCApostopv").c_str());
            // if (hDCApos == nullptr)
            // {
            //     cout << "DCA positive daughter to PV plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hDCApos);
            // hDCApos->GetYaxis()->SetTitle("Counts");
            // hDCApos->GetXaxis()->SetTitle("DCA pos. daughter to PV (cm)");
            // hDCApos->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_DCApostopv." + koutputtype).c_str());

            // DCA daughters
            TH1F *hDCAdaughters = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/hDCAV0Daughters").c_str());
            if (hDCAdaughters == nullptr)
            {
                cout << "DCA daughters plot not found" << endl;
                return;
            }
            c3->Clear();
            SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            SetHistoQA(hDCAdaughters);
            hDCAdaughters->GetYaxis()->SetTitle("Counts");
            hDCAdaughters->GetXaxis()->SetTitle("DCA V0 daughters (cm)");
            hDCAdaughters->Draw();
            c3->SaveAs((outputQAfolder_str + "/kshort_DCAV0Daughters." + koutputtype).c_str());

            // Kshort lifetime
            TH1F *hKshortLifetime = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/hLT").c_str());
            if (hKshortLifetime == nullptr)
            {
                cout << "Kshort lifetime plot not found" << endl;
                return;
            }
            c3->Clear();
            SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            SetHistoQA(hKshortLifetime);
            hKshortLifetime->GetYaxis()->SetTitle("Counts");
            hKshortLifetime->GetXaxis()->SetTitle("K_{s}^{0} lifetime (cm)");
            hKshortLifetime->Draw();
            c3->SaveAs((outputQAfolder_str + "/kshort_lifetime." + koutputtype).c_str());

            // n sigma neg pion daugter before
            TH2F *hNSigmaNegPion_before = (TH2F *)fInputFile->Get((foldername_final + "/kzeroShort/hNSigmaNegPionK0s_before").c_str());
            if (hNSigmaNegPion_before == nullptr)
            {
                cout << "nSigma negative pion daughter plot before selection cuts not found" << endl;
                return;
            }
            c3->Clear();
            SetCanvasStyle(c3, 0.15, 0.12, 0.05, 0.15);
            TH1D *nsigma_projection_neg = hNSigmaNegPion_before->ProjectionX();
            SetHistoQA(nsigma_projection_neg);
            SetHistoQA(hNSigmaNegPion_before);
            hNSigmaNegPion_before->GetYaxis()->SetTitle("n#sigma_{#pi^{-}}");
            hNSigmaNegPion_before->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            nsigma_projection_neg->GetXaxis()->SetTitle("n#sigma_{#pi^{-}}");
            nsigma_projection_neg->GetYaxis()->SetTitle("Counts");
            // gPad->SetLogx();
            // gPad->SetLogz();
            gPad->SetLogy();
            hNSigmaNegPion_before->SetStats(0);
            // hNSigmaNegPion_before->Draw("colz");
            nsigma_projection_neg->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_nSigmaNegPion." + koutputtype).c_str());

            // n sigma pos pion daugter before
            TH2F *hNSigmaPosPion_before = (TH2F *)fInputFile->Get((foldername_final + "/kzeroShort/hNSigmaPosPionK0s_before").c_str());
            if (hNSigmaPosPion_before == nullptr)
            {
                cout << "nSigma positive pion daughter plot before selection cuts not found" << endl;
                return;
            }
            TH1D *nsigma_projection_pos = hNSigmaPosPion_before->ProjectionX();
            // c3->Clear();
            SetCanvasStyle(c3, 0.15, 0.12, 0.05, 0.15);
            SetHistoQA(nsigma_projection_pos);
            SetHistoQA(hNSigmaPosPion_before);
            hNSigmaPosPion_before->GetYaxis()->SetTitle("n#sigma_{#pi^{+}}");
            hNSigmaPosPion_before->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            nsigma_projection_pos->GetXaxis()->SetTitle("n#sigma_{#pi^{+}}");
            nsigma_projection_pos->GetYaxis()->SetTitle("Counts");
            hNSigmaPosPion_before->SetStats(0);
            nsigma_projection_neg->Draw("colz");
            nsigma_projection_pos->SetLineColor(kRed);
            nsigma_projection_pos->Draw("same");
            TLegend *leg3 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg3->SetFillStyle(0);
            leg3->SetBorderSize(0);
            leg3->SetTextFont(42);
            leg3->AddEntry(nsigma_projection_neg, "#pi^{-}", "l");
            leg3->AddEntry(nsigma_projection_pos, "#pi^{+}", "l");
            leg3->Draw();
            c3->SaveAs((outputQAfolder_str + "/kshort_nSigma_compare." + koutputtype).c_str());

            // // n sigma neg pion daugter after
            // TH2F *hNSigmaNegPion_after = (TH2F *)fInputFile->Get((foldername_final + "/kzeroShort/hNSigmaNegPionK0s_after").c_str());
            // if (hNSigmaNegPion_after == nullptr)
            // {
            //     cout << "nSigma negative pion daughter plot after selection cuts not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.12, 0.05, 0.15);
            // SetHistoQA(hNSigmaNegPion_after);
            // hNSigmaNegPion_after->GetYaxis()->SetTitle("n#sigma_{#pi^{-}}");
            // hNSigmaNegPion_after->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            // hNSigmaNegPion_after->SetStats(0);
            // hNSigmaNegPion_after->Draw("colz");
            // c3->SaveAs((outputQAfolder_str + "/kshort_nSigmaNegPion_after." + koutputtype).c_str());

            // // n sigma pos pion daugter after
            // TH2F *hNSigmaPosPion_after = (TH2F *)fInputFile->Get((foldername_final + "/kzeroShort/hNSigmaPosPionK0s_after").c_str());
            // if (hNSigmaPosPion_after == nullptr)
            // {
            //     cout << "nSigma positive pion daughter plot after selection cuts not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.12, 0.05, 0.15);
            // SetHistoQA(hNSigmaPosPion_after);
            // hNSigmaPosPion_after->GetYaxis()->SetTitle("n#sigma_{#pi^{+}}");
            // hNSigmaPosPion_after->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            // hNSigmaPosPion_after->SetStats(0);
            // hNSigmaPosPion_after->Draw("colz");
            // c3->SaveAs((outputQAfolder_str + "/kshort_nSigmaPosPion_after." + koutputtype).c_str());

            gPad->SetLogx(0);
            gPad->SetLogz(0);
            gPad->SetLogy(0);

            // // psi pair angle plot
            // TH1F *hPsiPair = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/hpsipair").c_str());
            // if (hPsiPair == nullptr)
            // {
            //     cout << "Psi pair angle plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hPsiPair);
            // hPsiPair->GetYaxis()->SetTitle("Counts");
            // hPsiPair->GetXaxis()->SetTitle("#Psi_{pair} (rad)");
            // hPsiPair->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_psiPair." + koutputtype).c_str());

            // v0 cos PA
            TH1F *hV0CosPA = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/hV0CosPA").c_str());
            if (hV0CosPA == nullptr)
            {
                cout << "V0 cos PA plot not found" << endl;
                return;
            }
            c3->Clear();
            SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            SetHistoQA(hV0CosPA);
            hV0CosPA->GetYaxis()->SetTitle("Counts");
            hV0CosPA->GetXaxis()->SetTitle("V0 cos PA (rad)");
            hV0CosPA->Draw();
            c3->SaveAs((outputQAfolder_str + "/kshort_v0CosPA." + koutputtype).c_str());

            // // v0 radius
            // TH1F *hV0Radius = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/hv0radius").c_str());
            // if (hV0Radius == nullptr)
            // {
            //     cout << "V0 radius plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hV0Radius);
            // hV0Radius->GetYaxis()->SetTitle("Counts");
            // hV0Radius->GetXaxis()->SetTitle("V0 radius (cm)");
            // hV0Radius->GetXaxis()->SetRangeUser(0, 60);
            // hV0Radius->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_v0Radius." + koutputtype).c_str());

            // // negative daughter eta
            // TH1F *hNegDaughterEta = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/negative_eta").c_str());
            // if (hNegDaughterEta == nullptr)
            // {
            //     cout << "Negative daughter eta plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hNegDaughterEta);
            // hNegDaughterEta->GetYaxis()->SetTitle("Counts");
            // hNegDaughterEta->GetXaxis()->SetTitle("Neg. daughter #eta");
            // hNegDaughterEta->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_negDaughterEta." + koutputtype).c_str());

            // // positive daughter eta
            // TH1F *hPosDaughterEta = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/positive_eta").c_str());
            // if (hPosDaughterEta == nullptr)
            // {
            //     cout << "Positive daughter eta plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hPosDaughterEta);
            // hPosDaughterEta->GetYaxis()->SetTitle("Counts");
            // hPosDaughterEta->GetXaxis()->SetTitle("Pos. daughter #eta");
            // hPosDaughterEta->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_posDaughterEta." + koutputtype).c_str());

            // // negative daughter phi
            // TH1F *hNegDaughterPhi = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/negative_phi").c_str());
            // if (hNegDaughterPhi == nullptr)
            // {
            //     cout << "Negative daughter phi plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hNegDaughterPhi);
            // hNegDaughterPhi->GetYaxis()->SetTitle("Counts");
            // hNegDaughterPhi->GetXaxis()->SetTitle("Neg. daughter #phi");
            // hNegDaughterPhi->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_negDaughterPhi." + koutputtype).c_str());

            // // positive daughter phi
            // TH1F *hPosDaughterPhi = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/positive_phi").c_str());
            // if (hPosDaughterPhi == nullptr)
            // {
            //     cout << "Positive daughter phi plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hPosDaughterPhi);
            // hPosDaughterPhi->GetYaxis()->SetTitle("Counts");
            // hPosDaughterPhi->GetXaxis()->SetTitle("Pos. daughter #phi");
            // hPosDaughterPhi->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_posDaughterPhi." + koutputtype).c_str());

            // // negative daughter pT
            // gPad->SetLogy();
            // TH1F *hNegDaughterPt = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/negative_pt").c_str());
            // if (hNegDaughterPt == nullptr)
            // {
            //     cout << "Negative daughter pT plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hNegDaughterPt);
            // hNegDaughterPt->GetYaxis()->SetTitle("Counts");
            // hNegDaughterPt->GetXaxis()->SetTitle("Neg. daughter p_{T} (GeV/c)");
            // hNegDaughterPt->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_negDaughterPt." + koutputtype).c_str());

            // // positive daughter pT
            // TH1F *hPosDaughterPt = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/positive_pt").c_str());
            // if (hPosDaughterPt == nullptr)
            // {
            //     cout << "Positive daughter pT plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hPosDaughterPt);
            // hPosDaughterPt->GetYaxis()->SetTitle("Counts");
            // hPosDaughterPt->GetXaxis()->SetTitle("Pos. daughter p_{T} (GeV/c)");
            // hPosDaughterPt->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_posDaughterPt." + koutputtype).c_str());
            // gPad->SetLogy(0);

            // // negative daughter rapidity
            // TH1F *hNegDaughterRapidity = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/negative_y").c_str());
            // if (hNegDaughterRapidity == nullptr)
            // {
            //     cout << "Negative daughter rapidity plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hNegDaughterRapidity);
            // hNegDaughterRapidity->GetYaxis()->SetTitle("Counts");
            // hNegDaughterRapidity->GetXaxis()->SetTitle("Neg. daughter y");
            // hNegDaughterRapidity->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_negDaughterRapidity." + koutputtype).c_str());

            // // positive daughter rapidity
            // TH1F *hPosDaughterRapidity = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/positive_y").c_str());
            // if (hPosDaughterRapidity == nullptr)
            // {
            //     cout << "Positive daughter rapidity plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hPosDaughterRapidity);
            // hPosDaughterRapidity->GetYaxis()->SetTitle("Counts");
            // hPosDaughterRapidity->GetXaxis()->SetTitle("Pos. daughter y");
            // hPosDaughterRapidity->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_posDaughterRapidity." + koutputtype).c_str());

            // // Lambda mass
            // TH1F *hLambdaMass = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/Mass_lambda").c_str());
            // if (hLambdaMass == nullptr)
            // {
            //     cout << "Lambda mass plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hLambdaMass);
            // hLambdaMass->GetYaxis()->SetTitle("Counts");
            // hLambdaMass->GetXaxis()->SetTitle("m_{#pi^{-}p} (GeV/c^{2})");
            // hLambdaMass->GetXaxis()->SetRangeUser(1.0, 1.25);
            // hLambdaMass->Draw();
            // TLine *pdgmass_lambda = new TLine(1.115683, 0, 1.115683, hLambdaMass->GetMaximum());
            // pdgmass_lambda->SetLineColor(kRed);
            // pdgmass_lambda->SetLineStyle(2);
            // pdgmass_lambda->SetLineWidth(2);
            // pdgmass_lambda->Draw("same");
            // c3->SaveAs((outputQAfolder_str + "/lambda_mass." + koutputtype).c_str());

            // // Anti-Lambda mass
            // TH1F *hAntiLambdaMass = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/mass_AntiLambda").c_str());
            // if (hAntiLambdaMass == nullptr)
            // {
            //     cout << "Anti-Lambda mass plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hAntiLambdaMass);
            // hAntiLambdaMass->GetYaxis()->SetTitle("Counts");
            // // hAntiLambdaMass->GetXaxis()->SetTitle("m_{#bar{#Lambda}} (GeV/c^{2})");
            // hAntiLambdaMass->GetXaxis()->SetTitle("m_{#pi^{+}#bar{p}} (GeV/c^{2})");
            // hAntiLambdaMass->GetXaxis()->SetRangeUser(1.0, 1.25);
            // hAntiLambdaMass->Draw();
            // pdgmass_lambda->Draw("same");
            // c3->SaveAs((outputQAfolder_str + "/antilambda_mass." + koutputtype).c_str());

            // // rapidity updated plot
            // TH1F *hRapidity = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/rapidity").c_str());
            // if (hRapidity == nullptr)
            // {
            //     cout << "Rapidity plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hRapidity);
            // hRapidity->GetYaxis()->SetTitle("Counts");
            // hRapidity->GetXaxis()->SetTitle("K_{s}^{0} rapidity");
            // hRapidity->Draw();
            // c3->SaveAs((outputQAfolder_str + "/kshort_rapidity." + koutputtype).c_str());

            // // TPC dE/dx plot
            // // gPad->SetLogx();
            // // gPad->SetLogz();
            // TH2F *hTPCenergyloss = (TH2F *)fInputFile->Get((foldername_final + "/kzeroShort/dE_by_dx_TPC").c_str());
            // if (hTPCenergyloss == nullptr)
            // {
            //     cout << "TPC energy loss plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.14, 0.05, 0.15);
            // SetHistoQA(hTPCenergyloss);
            // hTPCenergyloss->GetYaxis()->SetTitle("TPC dE/dx (a.u.)");
            // hTPCenergyloss->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            // hTPCenergyloss->GetXaxis()->SetRangeUser(0.1, 50);
            // hTPCenergyloss->SetStats(0);
            // hTPCenergyloss->Draw("colz");
            // c3->SaveAs((outputQAfolder_str + "/kshort_TPCdEdx." + koutputtype).c_str());

            // Number of ks produced in an event plot
            // gPad->SetLogy();
            TH1F *hNofKshort = (TH1F *)fInputFile->Get((foldername_final + "/kzeroShort/NksProduced").c_str());
            if (hNofKshort == nullptr)
            {
                cout << "Number of Kshort produced in an event plot not found" << endl;
                return;
            }
            c3->Clear();
            gPad->SetLogy();
            SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            SetHistoQA(hNofKshort);
            // hNofKshort->Scale(1.0 / hNofKshort->Integral());
            hNofKshort->GetYaxis()->SetTitle("Counts");
            hNofKshort->GetXaxis()->SetTitle("Number of K^{0}_{s} produced");
            hNofKshort->GetYaxis()->SetMaxDigits(3);
            hNofKshort->SetMaximum(hNofKshort->GetMaximum() * 100);
            hNofKshort->Draw("HIST text");
            // cout << "Percentage of events in which > 1 Ks are produced " << hNofKshort->Integral(3, hNofKshort->GetNbinsX()) / hNofKshort->Integral(2, hNofKshort->GetNbinsX()) * 100 << endl;
            c3->SaveAs((outputQAfolder_str + "/NKs_produced." + koutputtype).c_str());

            // // multiplicity distribution FT0M
            // TH1F *hMult_FT0M = (TH1F *)fInputFile->Get((foldername_final + "/eventSelection/multdist_FT0M").c_str());
            // if (hMult_FT0M == nullptr)
            // {
            //     cout << "Multiplicity distribution FT0M plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hMult_FT0M);
            // hMult_FT0M->GetYaxis()->SetTitle("Counts");
            // hMult_FT0M->GetXaxis()->SetTitle("Multiplicity FT0M");
            // hMult_FT0M->GetXaxis()->SetRangeUser(-10000, 45000);
            // hMult_FT0M->GetXaxis()->SetNdivisions(505);
            // hMult_FT0M->Draw();
            // c3->SaveAs((outputQAfolder_str + "/mult_FT0M." + koutputtype).c_str());

            // // Mass correlation between lambda and Kshort
            // gPad->SetLogx(0);
            // gPad->SetLogz(0);
            // gPad->SetLogy(0);
            // TH2F *hMassCorr_ks_lambda_before = (TH2F *)fInputFile->Get((foldername_final + "/kzeroShort/mass_lambda_kshort_after9").c_str());
            // if (hMassCorr_ks_lambda_before == nullptr)
            // {
            //     cout << "Mass correlation plot between lambda and Kshort before the cut not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.16, 0.13, 0.05, 0.15);
            // SetHistoQA(hMassCorr_ks_lambda_before);
            // hMassCorr_ks_lambda_before->GetYaxis()->SetTitle("m_{#Lambda} (GeV/c^{2})");
            // hMassCorr_ks_lambda_before->GetXaxis()->SetTitle("M_{K^{0}_{s}} (GeV/c^{2})");
            // hMassCorr_ks_lambda_before->GetYaxis()->SetTitleOffset(1.6);
            // hMassCorr_ks_lambda_before->GetXaxis()->SetRangeUser(0.25, 0.78);
            // hMassCorr_ks_lambda_before->GetYaxis()->SetRangeUser(1.05, 1.5);
            // hMassCorr_ks_lambda_before->GetZaxis()->SetMaxDigits(3);
            // hMassCorr_ks_lambda_before->Draw("colz");
            // c3->SaveAs((outputQAfolder_str + "/kshort_mass_correlation_lambda_ks_before." + koutputtype).c_str());

            // // Mass correlation between lambda and Kshort after the cut
            // TH2F *hMassCorr_ks_lambda_after = (TH2F *)fInputFile->Get((foldername_final + "/kzeroShort/mass_lambda_kshort_after10").c_str());
            // if (hMassCorr_ks_lambda_after == nullptr)
            // {
            //     cout << "Mass correlation plot between lambda and Kshort after the cut not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // gPad->SetLogy(0);
            // gPad->SetLogz();
            // SetCanvasStyle(c3, 0.16, 0.13, 0.05, 0.15);
            // SetHistoQA(hMassCorr_ks_lambda_after);
            // hMassCorr_ks_lambda_after->GetYaxis()->SetTitle("m_{#Lambda} (GeV/c^{2})");
            // hMassCorr_ks_lambda_after->GetXaxis()->SetTitle("M_{K^{0}_{s}} (GeV/c^{2})");
            // hMassCorr_ks_lambda_after->GetYaxis()->SetTitleOffset(1.6);
            // hMassCorr_ks_lambda_after->GetXaxis()->SetRangeUser(0.25, 0.78);
            // hMassCorr_ks_lambda_after->GetYaxis()->SetRangeUser(1.05, 1.5);
            // hMassCorr_ks_lambda_after->GetZaxis()->SetMaxDigits(3);
            // hMassCorr_ks_lambda_after->Draw("colz");
            // c3->SaveAs((outputQAfolder_str + "/kshort_mass_correlation_lambda_ks_after." + koutputtype).c_str());

            // // Events check histogram
            // gPad->SetLogy(0);
            // string xlable_eventscheck[] = {"No cut", "No cut", "time frame", "sel8", "pileup", "goodzvertex", "itstpc_match", "additional evsel", "additional evsel"};
            // TH1I *hEventsCheck = (TH1I *)fInputFile->Get((foldername_final + "/hglueball/heventscheck").c_str());
            // if (hEventsCheck == nullptr)
            // {
            //     cout << "Events check plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hEventsCheck);
            // hEventsCheck->GetYaxis()->SetTitle("Counts");
            // for (int i = 0; i < 9; i++)
            // {
            //     hEventsCheck->GetXaxis()->SetBinLabel(i + 1, xlable_eventscheck[i].c_str());
            // }
            // hEventsCheck->SetStats(0);
            // hEventsCheck->SetMaximum(1.2 * hEventsCheck->GetMaximum());
            // hEventsCheck->Draw("HIST text");
            // c3->SaveAs((outputQAfolder_str + "/events_check." + koutputtype).c_str());

            // // Events check in v0 selection histogram
            // string xlable_v0selection[] = {"No cut", "DCAv0_to_PV", "rapidity Ks", "p_{T} min", "DCAv0_daughters", "Cos PA", "Tran_rad_min", "Tran_rad_max", "max lifetime", "armenteros", "competing"};
            // TH1I *hEventsCheckV0 = (TH1I *)fInputFile->Get((foldername_final + "/hglueball/htrackscheck_v0").c_str());
            // if (hEventsCheckV0 == nullptr)
            // {
            //     cout << "Events check in V0 selection plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hEventsCheckV0);
            // hEventsCheckV0->GetYaxis()->SetTitle("Counts");
            // hEventsCheckV0->GetXaxis()->SetRangeUser(0, 12);
            // for (int i = 0; i < 11; i++)
            // {
            //     hEventsCheckV0->GetXaxis()->SetBinLabel(i + 1, xlable_v0selection[i].c_str());
            // }
            // hEventsCheckV0->SetStats(0);
            // hEventsCheckV0->SetMaximum(1.2 * hEventsCheckV0->GetMaximum());
            // hEventsCheckV0->Draw("HIST text");
            // c3->SaveAs((outputQAfolder_str + "/events_check_v0." + koutputtype).c_str());

            // // Events check in v0 daughters selection histogram
            // string xlable_v0daughterselection[] = {"No cut", "has TPC", "TPC CR", "TPC CRFC", "charge_sign", "charge_sign", "eta", "TPC clusters", "PID"};
            // TH1I *hEventsCheckV0Daughters = (TH1I *)fInputFile->Get((foldername_final + "/hglueball/htrackscheck_v0_daughters").c_str());
            // if (hEventsCheckV0Daughters == nullptr)
            // {
            //     cout << "Events check in V0 daughters selection plot not found" << endl;
            //     return;
            // }
            // c3->Clear();
            // SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
            // SetHistoQA(hEventsCheckV0Daughters);
            // hEventsCheckV0Daughters->GetYaxis()->SetTitle("Counts");
            // hEventsCheckV0Daughters->GetXaxis()->SetRangeUser(0, 10);
            // for (int i = 0; i < 9; i++)
            // {
            //     hEventsCheckV0Daughters->GetXaxis()->SetBinLabel(i + 1, xlable_v0daughterselection[i].c_str());
            // }
            // hEventsCheckV0Daughters->SetStats(0);
            // hEventsCheckV0Daughters->SetMaximum(1.2 * hEventsCheckV0Daughters->GetMaximum());
            // hEventsCheckV0Daughters->Draw("HIST text");
            // c3->SaveAs((outputQAfolder_str + "/events_check_v0_daughters." + koutputtype).c_str());
        }
    }
}

float parameter0(float mass, float width)
{
    double gamma = TMath::Sqrt(mass * mass * (mass * mass + width * width));
    double norm = 2.8284 * mass * width * gamma / (3.14 * TMath::Sqrt(mass * mass + gamma));
    return norm;
}

// void printDirectoryContents(TDirectory *dir, int indent = 0)
// {
//     // Get a list of all keys in the directory
//     TIter next(dir->GetListOfKeys());
//     TKey *key;

//     // Iterate over all keys
//     while ((key = (TKey *)next()))
//     {
//         // Print the name and class of the object
//         for (int i = 0; i < indent; i++)
//         {
//             std::cout << "  ";
//         }
//         std::cout << key->GetName() << " (" << key->GetClassName() << ")" << std::endl;

//         // If the object is a directory, recursively print its contents
//         TClass *cl = gROOT->GetClass(key->GetClassName());
//         if (cl->InheritsFrom(TDirectory::Class()))
//         {
//             TDirectory *subdir = (TDirectory *)key->ReadObj();
//             printDirectoryContents(subdir, indent + 1);
//         }
//     }
// }
