#include <iostream>
#include <cmath>
#include <TArrow.h>
#include <TSystem.h>
#include <TString.h>
#include <TStopwatch.h>
#include "initializations.h"

using namespace std;
bool plot_all = false;

void buildTemplate()
{
    TStopwatch timer;
    timer.Start();

    //*************************** change here ***************************************
    // const string kResBkg = "MIX";
    const string kResBkg = "LIKE";
    // const string kResBkg = "ROTATED";
    TString outputtype = "pdf";
    const float txtsize = 0.045; // text size in the plots
    bool isINEL = false;

    const TString kvariation = "";
    // const TString kvariation = "_TPC1p5_combined2";
    // const TString kvariation = "_TPC2p5_combined3p5";
    // const TString kvariation = "_DCAvar1";
    // const TString kvariation = "_DCAvar2";
    // const TString kvariation = "_NoPVContributor";
    //********************************************************************************

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    // float mult_classes[] = {0.0};
    int nmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1; // number of multiplicity bins

    t2->SetNDC();
    t2->SetTextSize(0.06);
    t2->SetTextFont(42);

    //***************************************************************************************************

    // Input file
    TFile *fInputFile = new TFile("/home/sawan/Storage/check_k892/mc/LHC24f3c/750013.root", "Read");
    // TFile *fInputFile = new TFile("/home/sawan/Storage/check_k892/mc/LHC24f3c/746556.root", "Read");
    if (fInputFile->IsZombie())
    {
        cerr << "File not found " << endl;
        return;
    }

    TH1F *hcent = (TH1F *)fInputFile->Get(Form("kstarqa%s/eventSelection/hMultiplicity", kvariation.Data()));
    if (hcent == nullptr)
    {
        cerr << "Histogram not found" << endl;
        return;
    }

    const string recpath = Form("kstarqa%s/hInvMass/h3KstarMassRec", kvariation.Data());
    THnSparseF *hRec = (THnSparseF *)fInputFile->Get(recpath.c_str());

    //**Invariant mass histograms for sig+bkg and mixed event bg*****************************************
    auto fDirectory = (TDirectoryFile *)fInputFile->Get(Form("kstarqa%s/hInvMass", kvariation.Data()));
    THnSparseF *fHistLikeMM = (THnSparseF *)fDirectory->Get("h3KstarInvMasslikeSignMM");
    THnSparseF *fHistLikePP = (THnSparseF *)fDirectory->Get("h3KstarInvMasslikeSignPP");
    THnSparseF *fHistUnlike = (THnSparseF *)fDirectory->Get("h3KstarInvMassUnlikeSign");
    THnSparseF *fHistMix = (THnSparseF *)fDirectory->Get("h3KstarInvMassMixed");
    THnSparseF *fHistRotated = (THnSparseF *)fDirectory->Get("h3KstarInvMassRotated");

    if (fHistUnlike == nullptr || fHistMix == nullptr || fHistLikeMM == nullptr || fHistLikePP == nullptr || fHistRotated == nullptr)
    {
        cerr << "Invariant mass histograms not found!!!!!!!!!!!!" << endl;
        return;
    }

    TFile *outPutSigMinusTrue = new TFile(Form("template/%s/SignalMinusTrue%s.root", kResBkg.c_str(), kvariation.Data()), "RECREATE");

    for (int imult = 0; imult < nmultbins + 1; imult++)
    {
        int multlow, multhigh;

        if (imult == 0)
        {
            multlow = 0;
            multhigh = (isINEL) ? 120 : 100; // for all multiplicity
        }
        else
        {
            multlow = mult_classes[imult - 1];
            multhigh = mult_classes[imult];
        }
        TDirectory *dir = outPutSigMinusTrue->mkdir(Form("%d-%d", multlow, multhigh));

        //*************************Create folders********************************************
        TString centRange = Form("%d_%d", multlow, multhigh);
        // TString Cenoutputfolder = Form("template/%s/%d-%d%%", kResBkg.c_str(), multlow, multhigh);
        TString Cenoutputfolder = Form("template/%s/%s/%d-%d%%", kResBkg.c_str(), kvariation.Data(), multlow, multhigh);

        if (gSystem->mkdir(Cenoutputfolder, kTRUE))
        {
            std::cout << "Folder " << Cenoutputfolder << " created successfully." << std::endl;
        }

        int hcentbinlow = hcent->FindBin(multlow + 1e-5);
        int hcentbinhigh = hcent->FindBin(multhigh - 1e-5);
        double Event = hcent->Integral(hcentbinlow, hcentbinhigh);

        cout << "*****************number of events********************:" << Event << endl;

        // gstyle(); // defined in style.h
        gStyle->SetOptStat(1110);

        for (Int_t ip = pt_start; ip < pt_end; ip++) // start pt bin loop
        {
            lowpt = pT_bins[ip];
            highpt = pT_bins[ip + 1];

            int lbinptRec = hRec->GetAxis(0)->FindBin(lowpt + 1e-5);
            int hbinptRec = hRec->GetAxis(0)->FindBin(highpt - 1e-5);
            int lowbinCentRec = hRec->GetAxis(1)->FindBin(multlow + 1e-5);
            int highbinCentRec = hRec->GetAxis(1)->FindBin(multhigh - 1e-5);

            hRec->GetAxis(0)->SetRange(lbinptRec, hbinptRec);
            hRec->GetAxis(1)->SetRange(lowbinCentRec, highbinCentRec);

            TH1D *h1rec = hRec->Projection(2, "E");
            h1rec->SetName(Form("h1rec_pt_%d", ip));

            int lbincent = fHistUnlike->GetAxis(0)->FindBin(multlow + 1e-5);
            int hbincent = fHistUnlike->GetAxis(0)->FindBin(multhigh - 1e-5);
            int lbinpt = fHistUnlike->GetAxis(1)->FindBin(lowpt + 1e-5);
            int hbinpt = fHistUnlike->GetAxis(1)->FindBin(highpt - 1e-5);

            fHistUnlike->GetAxis(0)->SetRange(lbincent, hbincent);
            fHistUnlike->GetAxis(1)->SetRange(lbinpt, hbinpt);

            // Signal and Combinatorial histograms
            fHistTotal[ip] = fHistUnlike->Projection(2, "E");
            fHistTotal[ip]->SetName(Form("fHistTotal_%d", ip));

            if (kResBkg == "MIX")
            {
                lbincent = fHistMix->GetAxis(0)->FindBin(multlow + 1e-5);
                hbincent = fHistMix->GetAxis(0)->FindBin(multhigh - 1e-5);
                lbinpt = fHistMix->GetAxis(1)->FindBin(lowpt + 1e-5);
                hbinpt = fHistMix->GetAxis(1)->FindBin(highpt - 1e-5);

                fHistMix->GetAxis(0)->SetRange(lbincent, hbincent);
                fHistMix->GetAxis(1)->SetRange(lbinpt, hbinpt);

                fHistBkg[ip] = fHistMix->Projection(2, "E");
                fHistBkg[ip]->SetName(Form("fHistMix_%d", ip));
                fHistMix->SetName(Form("fHistMix_%d", ip));
            }
            else if (kResBkg == "ROTATED")
            {
                lbincent = fHistRotated->GetAxis(0)->FindBin(multlow + 1e-5);
                hbincent = fHistRotated->GetAxis(0)->FindBin(multhigh - 1e-5);
                lbinpt = fHistRotated->GetAxis(1)->FindBin(lowpt + 1e-5);
                hbinpt = fHistRotated->GetAxis(1)->FindBin(highpt - 1e-5);

                fHistRotated->GetAxis(0)->SetRange(lbincent, hbincent);
                fHistRotated->GetAxis(1)->SetRange(lbinpt, hbinpt);

                fHistRotated1D[ip] = fHistRotated->Projection(2, "E");
                fHistRotated1D[ip]->SetName(Form("fHistRotated1D_%d", ip));
                fHistRotated->SetName(Form("fHistRotated_%d", ip));
            }
            else if (kResBkg == "LIKE")
            {
                lbincent = fHistLikeMM->GetAxis(0)->FindBin(multlow + 1e-5);
                hbincent = fHistLikeMM->GetAxis(0)->FindBin(multhigh - 1e-5);
                lbinpt = fHistLikeMM->GetAxis(1)->FindBin(lowpt + 1e-5);
                hbinpt = fHistLikeMM->GetAxis(1)->FindBin(highpt - 1e-5);
                fHistLikeMM->SetName(Form("fHistLikeMM_%d", ip));

                fHistLikeMM->GetAxis(0)->SetRange(lbincent, hbincent);
                fHistLikeMM->GetAxis(1)->SetRange(lbinpt, hbinpt);

                fHistbkgMM[ip] = fHistLikeMM->Projection(2, "E");
                fHistbkgMM[ip]->SetName(Form("fHistbkgMM_%d", ip));

                lbincent = fHistLikePP->GetAxis(0)->FindBin(multlow + 1e-5);
                hbincent = fHistLikePP->GetAxis(0)->FindBin(multhigh - 1e-5);
                lbinpt = fHistLikePP->GetAxis(1)->FindBin(lowpt + 1e-5);
                hbinpt = fHistLikePP->GetAxis(1)->FindBin(highpt - 1e-5);
                fHistLikePP->SetName(Form("fHistLikePP_%d", ip));

                fHistLikePP->GetAxis(0)->SetRange(lbincent, hbincent);
                fHistLikePP->GetAxis(1)->SetRange(lbinpt, hbinpt);

                fHistbkgPP[ip] = fHistLikePP->Projection(2, "E");
                fHistbkgPP[ip]->SetName(Form("fHistbkgPP_%d", ip));

                auto tempLS = (TH1D *)fHistbkgMM[ip]->Clone("tempLS");
                tempLS->Multiply(fHistbkgPP[ip]);

                fHistbkgLS[ip] = (TH1D *)tempLS->Clone(Form("fHistbkgLS_%d", ip));
                fHistbkgLS[ip]->Reset(); // Clear contents

                for (int i = 1; i < tempLS->GetNbinsX() + 1; i++)
                {
                    double ppnn = tempLS->GetBinContent(i);
                    double ppnnerr = tempLS->GetBinError(i);
                    double err = ppnnerr / std::sqrt(ppnn);
                    if (ppnn == 0)
                    {
                        fHistbkgLS[ip]->SetBinContent(i, 0);
                        fHistbkgLS[ip]->SetBinError(i, 0);
                    }
                    else
                    {
                        fHistbkgLS[ip]->SetBinContent(i, 2 * std::sqrt(ppnn));
                        fHistbkgLS[ip]->SetBinError(i, err);
                    }
                }
                fHistLikeMM->SetName(Form("fHistLikeMM_%d", ip));
                fHistLikePP->SetName(Form("fHistLikePP_%d", ip));
            }

            //**Cloning sig+bkg histogram for like sign, mixed event, or rotated subtraction ******************
            TH1D *hfsig = (TH1D *)fHistTotal[ip]->Clone();
            auto binwidth_file = (fHistTotal[ip]->GetXaxis()->GetXmax() - fHistTotal[ip]->GetXaxis()->GetXmin()) * kRebin[ip] / fHistTotal[ip]->GetXaxis()->GetNbins();
            cout << "The value of binwidth_file is: " << binwidth_file << endl;

            if (kResBkg == "MIX")
            {
                TH1D *bkgclonetemp = (TH1D *)fHistBkg[ip]->Clone();
                sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                bkg_integral = (bkgclonetemp->Integral(bkgclonetemp->GetXaxis()->FindBin(kNormRangepT[ip][0]), bkgclonetemp->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                normfactor = sigbkg_integral / bkg_integral;

                hfbkg = (TH1D *)bkgclonetemp->Clone();
                hfbkg->Scale(normfactor);
                hfbkg->Rebin(kRebin[ip]);
                hfsig->Rebin(kRebin[ip]);
                hfsig->Add(hfbkg, -1);
            }
            else if (kResBkg == "LIKE")
            {
                hfbkg = (TH1D *)fHistbkgLS[ip]->Clone();
                hfbkg->Rebin(kRebin[ip]);
                hfsig->Rebin(kRebin[ip]);
                hfsig->Add(hfbkg, -1);
            }
            else if (kResBkg == "ROTATED")
            {
                TH1D *bkgclonetemp = (TH1D *)fHistRotated1D[ip]->Clone();
                sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                bkg_integral = (bkgclonetemp->Integral(bkgclonetemp->GetXaxis()->FindBin(kNormRangepT[ip][0]), bkgclonetemp->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                normfactor = sigbkg_integral / bkg_integral;

                hfbkg = (TH1D *)bkgclonetemp->Clone();
                hfbkg->Scale(normfactor);
                hfbkg->Rebin(kRebin[ip]);
                hfsig->Rebin(kRebin[ip]);
                hfsig->Add(hfbkg, -1);
            }

            fHistTotal[ip]->Rebin(kRebin[ip]);
            h1rec->Rebin(kRebin[ip]);

            //**Signal Minus True Generation *****************************************************************
            TH1D *hSigminusTrue = (TH1D *)hfsig->Clone("hSigminusTrue");
            hSigminusTrue->Add(h1rec, -1);
            // hSigminusTrue->Add(hTotalRef, -1); // Kept commented as in original code

            TCanvas *cSigminusTrue = new TCanvas(Form("cSigminusTrue_pt_%d", ip + 1), Form("cSigminusTrue_pt_%d", ip + 1), 720, 720);
            cSigminusTrue->SetTopMargin(0.08);
            cSigminusTrue->SetBottomMargin(0.13);
            cSigminusTrue->SetLeftMargin(0.15);
            cSigminusTrue->SetRightMargin(0.03);

            hSigminusTrue->SetTitle(Form("%.1f < #it{p}_{T} (GeV/#it{c}) < %.1f; M_{K#pi} (GeV/c^{2}); Counts", lowpt, highpt));
            hSigminusTrue->Draw("ep");

            cSigminusTrue->SaveAs(Form((Cenoutputfolder + "/hSigminusTrue_pt%d." + outputtype).Data(), ip + 1));

            outPutSigMinusTrue->cd();
            hSigminusTrue->SetName(Form("hSigminusTrue_pt_%.1f_%.1f", lowpt, highpt));
            dir->cd();
            hSigminusTrue->Write();

            delete cSigminusTrue; // Free memory at end of each loop iteration

        } // pt loop ends
    }

    outPutSigMinusTrue->Close();

    // Stop the stopwatch
    timer.Stop();

    // Get the elapsed time
    Double_t realTime = timer.RealTime(); // Wall clock time in seconds
    Double_t cpuTime = timer.CpuTime();   // CPU time used in seconds

    // Print the elapsed times
    std::cout << "Real time elapsed: " << realTime << " seconds" << std::endl;
    std::cout << "CPU time used: " << cpuTime << " seconds" << std::endl;
}