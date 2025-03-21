
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

void glueball_KK_channel()
{
    // change here ***********************************************************
    // const string kResBkg = "MIX";
    const string kResBkg = "LIKE";
    // const string kResBkg = "ROTATED";
    const bool makeQAplots = false;
    const bool calculate_inv_mass = true;
    const bool save_invmass_distributions = true;
    // change here ***********************************************************

    TString outputfolder = kSignalOutput + "/" + kchannel + "/" + kfoldername;
    TString outputQAfolder = kSignalOutput + "/" + kchannel + "/" + kfoldername + "/QA";
    const string outputfolder_str = kSignalOutput + "/" + kchannel + "/" + kfoldername;
    const string outputQAfolder_str = kSignalOutput + "/" + kchannel + "/" + kfoldername + "/QA";
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
    gStyle->SetOptStat(1110);

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

    TH1F *hentries = (TH1F *)fInputFile->Get("event-selection-task/hColCounterAcc");
    double Event = hentries->GetEntries();
    cout << "*******number of events from the event selection histogram is *******:" << Event << endl;

    //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************

    THnSparseF *fHistUnlike = (THnSparseF *)fInputFile->Get(Form("%s/h3PhiInvMassUnlikeSign", kfoldername.c_str()));
    THnSparseF *fHistMix = (THnSparseF *)fInputFile->Get(Form("%s/h3PhiInvMassMixed", kfoldername.c_str()));
    THnSparseF *fHistRot = (THnSparseF *)fInputFile->Get(Form("%s/h3PhiInvMassRotation", kfoldername.c_str()));
    THnSparseF *fHistLike_pp = (THnSparseF *)fInputFile->Get(Form("%s/h3PhiInvMassLikeSignPP", kfoldername.c_str()));
    THnSparseF *fHistLike_mm = (THnSparseF *)fInputFile->Get(Form("%s/h3PhiInvMassLikeSignMM", kfoldername.c_str()));

    if (fHistUnlike == nullptr || fHistMix == nullptr || fHistRot == nullptr || fHistLike_pp == nullptr || fHistLike_mm == nullptr)
    {
        cout << "Invariant mass histograms not found" << endl;
        return;
    }

    cout << " The number of entries in histograms: \n"
         << "same event: " << fHistUnlike->GetEntries() << "\n"
         << "mixed event: " << fHistMix->GetEntries() << "\n"
         << "rotated bkg: " << fHistRot->GetEntries() << "\n"
         << "like sign pp: " << fHistLike_pp->GetEntries() << "\n"
         << "like sign mm: " << fHistLike_mm->GetEntries() << endl;

    TH1D *fHistTotal[Npt];
    TH1D *fHistME[Npt];
    TH1D *fHistRotated[Npt];
    TH1D *fHistLikepp[Npt];
    TH1D *fHistLikemm[Npt];
    TH1D *fHistLike[Npt];

    int multlow = 0;
    int multhigh = 100;

    TH1F *hmult = (TH1F *)fInputFile->Get(Form("%s/hmutiplicity", kfoldername.c_str()));
    if (hmult == nullptr)
    {
        cout << "Multiplicity histogram not found" << endl;
        return;
    }
    double realevents = hmult->Integral(hmult->GetXaxis()->FindBin(multlow), hmult->GetXaxis()->FindBin(multhigh));
    cout << "*******number of events from the multiplicity histogram is *******:" << realevents << endl;

    if (calculate_inv_mass)
    {
        TFile *fileInvDistPair;
        fileInvDistPair = (Npt > 1) ? new TFile((outputfolder_str + "/hglue_" + kResBkg + Form("_norm_%.2f_%.2f", kNormRangepT[0][0], kNormRangepT[0][1]) + "_allpt" + ".root").c_str(), "RECREATE") : new TFile((outputfolder_str + "/hglue_" + kResBkg + Form("_norm_%.2f_%.2f", kNormRangepT[0][0], kNormRangepT[0][1]) + Form("_pt_%.1f_%.1f", pT_bins[0], pT_bins[1]) + ".root").c_str(), "RECREATE");

        for (Int_t ip = pt_start; ip < pt_end; ip++) // start pt bin loop
        {

            float lowpt = pT_bins[ip];
            float highpt = pT_bins[ip + 1];
            int lbin = fHistUnlike->GetAxis(1)->FindBin(lowpt + 1e-5);
            int hbin = fHistUnlike->GetAxis(1)->FindBin(highpt - 1e-5);

            fHistUnlike->GetAxis(1)->SetRange(lbin, hbin);
            fHistMix->GetAxis(1)->SetRange(lbin, hbin);
            fHistRot->GetAxis(1)->SetRange(lbin, hbin);
            fHistLike_pp->GetAxis(1)->SetRange(lbin, hbin);
            fHistLike_mm->GetAxis(1)->SetRange(lbin, hbin);

            int lbinmult = fHistUnlike->GetAxis(0)->FindBin(multlow + 1e-5);
            int hbinmult = fHistUnlike->GetAxis(0)->FindBin(multhigh - 1e-5);

            fHistUnlike->GetAxis(0)->SetRange(lbinmult, hbinmult);
            fHistMix->GetAxis(0)->SetRange(lbinmult, hbinmult);
            fHistRot->GetAxis(0)->SetRange(lbinmult, hbinmult);
            fHistLike_pp->GetAxis(0)->SetRange(lbinmult, hbinmult);
            fHistLike_mm->GetAxis(0)->SetRange(lbinmult, hbinmult);

            fHistTotal[ip] = fHistUnlike->Projection(2, "E");
            fHistME[ip] = fHistMix->Projection(2, "E");
            fHistRotated[ip] = fHistRot->Projection(2, "E");
            fHistLikepp[ip] = fHistLike_pp->Projection(2, "E");
            fHistLikemm[ip] = fHistLike_mm->Projection(2, "E");
            fHistTotal[ip]->SetName(Form("fHistTotal_%d", ip));
            fHistME[ip]->SetName(Form("fHistME_%d", ip));
            fHistRotated[ip]->SetName(Form("fHistRotated_%d", ip));
            fHistLikepp[ip]->SetName(Form("fHistLikepp_%d", ip));
            fHistLikemm[ip]->SetName(Form("fHistLikemm_%d", ip));
            fHistLike[ip] = (TH1D *)fHistLikepp[ip]->Clone();

            for (int ibin = 0; ibin < fHistLikepp[ip]->GetNbinsX(); ibin++)
            {
                fHistLike[ip]->SetBinContent(ibin + 1, 2 * TMath::Sqrt(fHistLikepp[ip]->GetBinContent(ibin + 1) * fHistLikemm[ip]->GetBinContent(ibin + 1))); // direct sum of like sign pairs
            }

            fHistLike[ip]->SetName(Form("fHistLike_%d", ip));
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
                auto sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                auto bkg_integral = (fHistME[ip]->Integral(fHistME[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistME[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                auto bkg_integral_rotated = (fHistRotated[ip]->Integral(fHistRotated[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistRotated[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));

                auto normfactor = sigbkg_integral / bkg_integral;                 // scaling factor for mixed bkg
                auto normfactor_rotated = sigbkg_integral / bkg_integral_rotated; // scaling factor for rotated bkg
                cout << "\n\n normalization factor " << 1. / normfactor << "\n\n";
                if (kResBkg == "MIX")
                {
                    hfbkg = (TH1D *)fHistME[ip]->Clone();
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
            else if (kResBkg == "LIKE")
            {
                hfbkg = (TH1D *)fHistLike[ip]->Clone();
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
            hfsig->SetMarkerSize(0.5);
            hfsig->GetYaxis()->SetMaxDigits(3);
            hfsig->GetYaxis()->SetTitleOffset(1.5);
            hfsig->SetMarkerColor(kBlack);
            hfsig->SetLineColor(kBlack);
            hfsig->GetXaxis()->SetTitle("m_{K^{#pm}K^{#mp}} (GeV/c^{2})");
            hfsig->GetYaxis()->SetTitle(Form("Counts/%.3f GeV/c^{2}", binwidth_file));
            hfsig->GetXaxis()->SetRangeUser(1.0, 2.8);
            hfsig->Write(Form("ksks_subtracted_invmass_pt_%.1f_%.1f", lowpt, highpt));
            hfsig->SetMaximum(2.9e6);
            hfsig->Draw("e");
            t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
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
            if (save_invmass_distributions)
            {
                c1->SaveAs((outputfolder_str + "/hglueball_signal_" + kResBkg + Form("pT_%.1f_%.1f_norm_%.2f_%.2f.", pT_bins[ip], pT_bins[ip + 1], kNormRangepT[ip][0], kNormRangepT[ip][1]) + koutputtype).c_str());
            }

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
                if (hbkg_nopeak->GetBinCenter(i + 1) < kNormRangepT[ip][0] || hbkg_nopeak->GetBinCenter(i + 1) > kNormRangepT[ip][1])
                {
                    hbkg_nopeak->SetBinContent(i + 1, -999);
                }
            }

            fHistTotal[ip]->SetMarkerStyle(8);
            fHistTotal[ip]->SetMarkerColor(kBlack);
            fHistTotal[ip]->SetMarkerSize(0.5);
            hfbkg->SetMarkerStyle(8);
            hfbkg->SetMarkerSize(0.5);
            hfbkg->SetMarkerColor(kRed);
            hfbkg->SetLineColor(kRed);
            fHistTotal[ip]->GetYaxis()->SetMaxDigits(3);
            fHistTotal[ip]->GetYaxis()->SetTitleOffset(1.5);
            fHistTotal[ip]->GetYaxis()->SetTitle(Form("Counts/%.3f GeV/c^{2}", binwidth_file));
            // fHistTotal[ip]->SetMaximum(1.2 * fHistTotal[ip]->GetMaximum());
            fHistTotal[ip]->Draw("E");
            fHistTotal[ip]->Write(Form("ksks_invmass_pt_%.1f_%.1f", lowpt, highpt));
            hfbkg->Write(Form("ksks_bkg_pt_%.1f_%.1f", lowpt, highpt));
            if (save_invmass_distributions)
            {
                c2->SaveAs((outputfolder_str + "/hglueball_invmass_only_." + Form("pT_%.1f_%.1f_.", pT_bins[ip], pT_bins[ip + 1]) + koutputtype).c_str());
            }
            hfbkg->Draw("E same");
            if (kResBkg == "MIX" || kResBkg == "ROTATED")
                hbkg_nopeak->Draw("BAR same");

            TLegend *leg = new TLegend(0.2451253, 0.2054598, 0.5445682, 0.3908046);
            leg->SetFillStyle(0);
            leg->SetBorderSize(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.04);
            leg->AddEntry(fHistTotal[ip], "Signal", "lpe");
            string bkgname = (kResBkg == "MIX") ? "Mixed event" : (kResBkg == "LIKE") ? "Like-sign bkg"
                                                                                      : "Rotated bkg";
            leg->AddEntry(hfbkg, bkgname.c_str(), "lpe");
            if (kResBkg == "MIX")
                leg->AddEntry(hbkg_nopeak, "Norm. region", "f");
            leg->Draw();
            t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
            if (save_invmass_distributions)
            {
                c2->SaveAs((outputfolder_str + "/hglueball_invmass_" + kResBkg + Form("pT_%.1f_%.1f_norm_%.2f_%.2f.", pT_bins[ip], pT_bins[ip + 1], kNormRangepT[ip][0], kNormRangepT[ip][1]) + koutputtype).c_str());
            }
            c2->Write(Form("ksks_invmass_withbkg_pt_%.1f_%.1f", lowpt, highpt));

            // cdivide->cd(ip + 1);
            // fHistTotal[ip]->Draw("E");
            // hfbkg->Draw("E same");
            // if (kResBkg == "MIX" || kResBkg == "ROTATED")
            //     hbkg_nopeak->Draw("BAR same");
            // leg->Draw();
            // t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
        } // pt bin loop end here
    }
    ////////////////////////////////////////////////////////////////////////
    // QA plots
    // multiplicity percentile plot
    if (!makeQAplots)
    {
        cout << "QA plots are not made; quitting \n";
        return;
    }
    TCanvas *c3 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hmult);
    hmult->GetXaxis()->SetTitle("Multiplicity percentile");
    hmult->GetYaxis()->SetTitle("Counts");
    hmult->GetXaxis()->SetRangeUser(0, 150);
    hmult->Draw();
    c3->SaveAs((outputQAfolder_str + "/hglueball_mult." + koutputtype).c_str());

    // multiplicity distribution
    gPad->SetLogy();
    TH1F *hmultdist_FT0M = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/multdist_FT0M").c_str());
    if (hmultdist_FT0M == nullptr)
    {
        cout << "Multiplicity distribution histogram not found" << endl;
        return;
    }
    c3->Clear();
    SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hmultdist_FT0M);
    hmultdist_FT0M->GetXaxis()->SetTitle("Multiplicity");
    hmultdist_FT0M->GetYaxis()->SetTitle("Events");
    hmultdist_FT0M->GetXaxis()->SetNdivisions(505);
    hmultdist_FT0M->Scale(1. / hmultdist_FT0M->Integral());
    // hmultdist_FT0M->SetMinimum(1e-8);
    hmultdist_FT0M->Draw();
    c3->SaveAs((outputQAfolder_str + "/hglueball_multdist_FT0M." + koutputtype).c_str());

    // vertex z position plot
    gPad->SetLogy(0);
    TH1F *hvertexz = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/hVtxZ").c_str());
    if (hvertexz == nullptr)
    {
        cout << "Vertex z position histogram not found" << endl;
        return;
    }
    c3->Clear();
    SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hvertexz);
    hvertexz->GetXaxis()->SetTitle("Vertex z position (cm)");
    hvertexz->GetYaxis()->SetTitle("Counts");
    hvertexz->GetXaxis()->SetRangeUser(-15, 15);
    hvertexz->Draw();
    c3->SaveAs((outputQAfolder_str + "/hglueball_vertexz." + koutputtype).c_str());

    // eta distribution
    TH1F *heta = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/hEta").c_str());
    if (heta == nullptr)
    {
        cout << "Eta distribution histogram not found" << endl;
        return;
    }
    c3->Clear();
    SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(heta);
    heta->GetXaxis()->SetTitle("#eta");
    heta->GetYaxis()->SetTitle("Counts");
    heta->Draw();
    c3->SaveAs((outputQAfolder_str + "/hglueball_eta." + koutputtype).c_str());

    // DCAxy distribution
    gPad->SetLogy();
    TH1F *hdcaxy = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/hDcaxy").c_str());
    if (hdcaxy == nullptr)
    {
        cout << "DCAxy distribution histogram not found" << endl;
        return;
    }
    c3->Clear();
    SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hdcaxy);
    hdcaxy->GetXaxis()->SetTitle("DCA_{xy} (cm)");
    hdcaxy->GetYaxis()->SetTitle("Counts");
    hdcaxy->GetXaxis()->SetRangeUser(-0.2, 0.2);
    hdcaxy->Draw();
    c3->SaveAs((outputQAfolder_str + "/hglueball_dca." + koutputtype).c_str());

    // DCAz distribution
    TH1F *hdcaz = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/hDcaz").c_str());
    if (hdcaz == nullptr)
    {
        cout << "DCAz distribution histogram not found" << endl;
        return;
    }
    c3->Clear();
    SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hdcaz);
    hdcaz->GetXaxis()->SetTitle("DCA_{z} (cm)");
    hdcaz->GetYaxis()->SetTitle("Counts");
    hdcaz->GetXaxis()->SetRangeUser(-0.2, 0.2);
    hdcaz->Draw();
    c3->SaveAs((outputQAfolder_str + "/hglueball_dcaz." + koutputtype).c_str());

    // nsigma TPC distribution before
    gPad->SetLogy(0);
    gPad->SetLogx(0);
    gPad->SetLogz(1);
    TH2F *hnsigmaTPC_before = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/hNsigmaKaonTPC_before").c_str());
    if (hnsigmaTPC_before == nullptr)
    {
        cout << "Nsigma TPC distribution histogram not found" << endl;
        return;
    }
    c3->Clear();
    SetCanvasStyle(c3, 0.15, 0.18, 0.05, 0.15);
    SetHistoQA2D(hnsigmaTPC_before);
    hnsigmaTPC_before->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hnsigmaTPC_before->GetYaxis()->SetTitle("n#sigma_{TPC}");
    hnsigmaTPC_before->Draw("colz");
    c3->SaveAs((outputQAfolder_str + "/hglueball_nsigmaTPC_before." + koutputtype).c_str());

    // // nsigma TPC distribution after
    // TH2F *hnsigmaTPC_after = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/hNsigmaKaonTPC_after").c_str());
    // if (hnsigmaTPC_after == nullptr)
    // {
    //     cout << "Nsigma TPC distribution histogram not found" << endl;
    //     return;
    // }
    // c3->Clear();
    // SetCanvasStyle(c3, 0.15, 0.18, 0.05, 0.15);
    // SetHistoQA2D(hnsigmaTPC_after);
    // hnsigmaTPC_after->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // hnsigmaTPC_after->GetYaxis()->SetTitle("n#sigma_{TPC}");
    // hnsigmaTPC_after->Draw("colz");
    // c3->SaveAs((outputQAfolder_str + "/hglueball_nsigmaTPC_after." + koutputtype).c_str());

    // nsigma TOF distribution before
    TH2F *hnsigmaTOF_before = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/hNsigmaKaonTOF_before").c_str());
    if (hnsigmaTOF_before == nullptr)
    {
        cout << "Nsigma TOF distribution histogram not found" << endl;
        return;
    }
    c3->Clear();
    SetCanvasStyle(c3, 0.15, 0.18, 0.05, 0.15);
    SetHistoQA2D(hnsigmaTOF_before);
    hnsigmaTOF_before->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hnsigmaTOF_before->GetYaxis()->SetTitle("n#sigma_{TOF}");
    hnsigmaTOF_before->Draw("colz");
    c3->SaveAs((outputQAfolder_str + "/hglueball_nsigmaTOF_before." + koutputtype).c_str());

    // // nsigma TOF distribution after
    // TH2F *hnsigmaTOF_after = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/hNsigmaKaonTOF_after").c_str());
    // if (hnsigmaTOF_after == nullptr)
    // {
    //     cout << "Nsigma TOF distribution histogram not found" << endl;
    //     return;
    // }
    // c3->Clear();
    // SetCanvasStyle(c3, 0.15, 0.18, 0.05, 0.15);
    // SetHistoQA2D(hnsigmaTOF_after);
    // hnsigmaTOF_after->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // hnsigmaTOF_after->GetYaxis()->SetTitle("n#sigma_{TOF}");
    // hnsigmaTOF_after->Draw("colz");
    // c3->SaveAs((outputQAfolder_str + "/hglueball_nsigmaTOF_after." + koutputtype).c_str());

    // nsigma TPC vs nsigma TOF before
    TH2F *hnsigmaTPCvsTOF_before = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/hNsigmaKaonTOF_TPC_before").c_str());
    if (hnsigmaTPCvsTOF_before == nullptr)
    {
        cout << "Nsigma TPC vs TOF distribution before selection histogram not found" << endl;
        return;
    }
    c3->Clear();
    SetCanvasStyle(c3, 0.15, 0.18, 0.05, 0.15);
    SetHistoQA2D(hnsigmaTPCvsTOF_before);
    hnsigmaTPCvsTOF_before->GetXaxis()->SetTitle("n#sigma_{TOF}");
    hnsigmaTPCvsTOF_before->GetYaxis()->SetTitle("n#sigma_{TPC}");
    hnsigmaTPCvsTOF_before->Draw("colz");
    c3->SaveAs((outputQAfolder_str + "/hglueball_nsigmaTPCvsTOF_before." + koutputtype).c_str());

    // // nsigma TPC vs nsigma TOF after
    // TH2F *hnsigmaTPCvsTOF_after = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/hNsigmaKaonTOF_TPC_after").c_str());
    // if (hnsigmaTPCvsTOF_after == nullptr)
    // {
    //     cout << "Nsigma TPC vs TOF distribution after selection histogram not found" << endl;
    //     return;
    // }
    // c3->Clear();
    // SetCanvasStyle(c3, 0.15, 0.18, 0.05, 0.15);
    // SetHistoQA2D(hnsigmaTPCvsTOF_after);
    // hnsigmaTPCvsTOF_after->GetXaxis()->SetTitle("n#sigma_{TOF}");
    // hnsigmaTPCvsTOF_after->GetYaxis()->SetTitle("n#sigma_{TPC}");
    // hnsigmaTPCvsTOF_after->Draw("colz");
    // c3->SaveAs((outputQAfolder_str + "/hglueball_nsigmaTPCvsTOF_after." + koutputtype).c_str());

    // End of code **********************************************************************************************
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
