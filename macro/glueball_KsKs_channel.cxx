
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
    const bool makeQAplots = false;
    const bool save_invmass_distributions = true;
    // change here ***********************************************************

    TString outputfolder = kSignalOutput + "/" + kchannel + "/" + kfoldername;
    TString outputQAfolder = kSignalOutput + "/" + kchannel + "/" + kfoldername + "/QA";
    const string outputfolder_str = kSignalOutput + "/" + kchannel + "/" + kfoldername;
    const string outputQAfolder_str = kSignalOutput + "/" + kchannel + "/" + kfoldername + "/QA";
    // Create the folder using TSystem::mkdir()
    if (gSystem->mkdir(outputfolder, kTRUE) || gSystem->mkdir(outputQAfolder, kTRUE))
    {
        std::cout << "Folder " << outputfolder << " created successfully." << std::endl;
        std::cout << "Folder " << outputQAfolder << " created successfully." << std::endl;
    }
    else
    {
        std::cout << "Creating folder " << outputfolder << std::endl;
        std::cout << "Creating folder " << outputQAfolder << std::endl;
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

    TH1F *hentries = (TH1F *)fInputFile->Get("event-selection-task/hColCounterAcc");
    double Event = hentries->GetEntries();
    cout << "*******number of events from the event selection histogram is *******:" << Event << endl;

    //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************

    THnSparseF *fHistNum = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassDS", kfoldername.c_str()));
    THnSparseF *fHistDen = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassME", kfoldername.c_str()));
    THnSparseF *fHistRot = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassRot", kfoldername.c_str()));

    if (fHistNum == nullptr || fHistDen == nullptr || fHistRot == nullptr)
    {
        cout << "Invariant mass histograms not found" << endl;
        return;
    }
    cout << " The number of entries in histograms: \n"
         << "same event: " << fHistNum->GetEntries() << "\n"
         << "mixed event: " << fHistDen->GetEntries() << "\n"
         << "rotated bkg/2: " << fHistRot->GetEntries() / 2 << endl;

    TH1D *fHistTotal[Npt];
    TH1D *fHistBkg[Npt];
    TH1D *fHistRotated[Npt];
    TH1F *hmult = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/eventSelection/hmultiplicity").c_str());
    if (hmult == nullptr)
    {
        cout << "Multiplicity histogram not found" << endl;
        return;
    }
    int multlow = 0;
    int multhigh = 100;
    double realevents = hmult->Integral(hmult->GetXaxis()->FindBin(multlow), hmult->GetXaxis()->FindBin(multhigh));
    cout << "*******number of events from the multiplicity histogram is *******:" << realevents << endl;

    for (Int_t ip = pt_start; ip < pt_end; ip++) // start pt bin loop
    {

        float lowpt = pT_bins[ip];
        float highpt = pT_bins[ip + 1];
        cout << "low pt value is " << lowpt << " high pt value is " << highpt << endl;
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

        if (kResBkg == "MIX")
        {
            auto sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
            auto bkg_integral = (fHistBkg[ip]->Integral(fHistBkg[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistBkg[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
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
        TFile *fileInvDistPair = new TFile((outputfolder_str + "/hglue_KsKs_" + kResBkg + Form("_%f_%f.root", pT_bins[ip], pT_bins[ip + 1])).c_str(), "RECREATE");
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
        hfsig->GetYaxis()->SetTitle(Form("Counts/%.3f GeV/c^{2}", binwidth_file));
        hfsig->GetXaxis()->SetRangeUser(1.1, 2.3);
        hfsig->Draw("e");
        hfsig->Write("ksks_invmass");
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
            c1->SaveAs((outputfolder_str + "/hglueball_signal_" + kResBkg + Form("_%d.", ip) + koutputtype).c_str());

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
        hfbkg->SetMarkerStyle(8);
        hfbkg->SetMarkerColor(kRed);
        hfbkg->SetLineColor(kRed);
        fHistTotal[ip]->GetYaxis()->SetMaxDigits(3);
        fHistTotal[ip]->GetYaxis()->SetTitleOffset(1.4);
        fHistTotal[ip]->Draw("E");
        fHistTotal[ip]->GetYaxis()->SetTitle(Form("Counts/%.3f GeV/c^{2}", binwidth_file));
        hfbkg->Draw("E same");
        if (kResBkg == "MIX")
            hbkg_nopeak->Draw("BAR same");

        TLegend *leg = new TLegend(0.2451253, 0.2054598, 0.5445682, 0.3908046);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);
        leg->AddEntry(fHistTotal[ip], "Signal", "lpe");
        string bkgname = (kResBkg == "MIX") ? "Mixed event" : "Rotated bkg";
        leg->AddEntry(hfbkg, bkgname.c_str(), "lpe");
        if (kResBkg == "MIX")
            leg->AddEntry(hbkg_nopeak, "Norm. region", "f");
        leg->Draw();
        t2->DrawLatex(0.27, 0.96, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", lowpt, highpt));
        if (save_invmass_distributions)
            c2->SaveAs((outputfolder_str + "/hglueball_invmass_" + kResBkg + Form("_%d.", ip) + koutputtype).c_str());
    } // pt bin loop end here
    ////////////////////////////////////////////////////////////////////////
    // QA plots here
    if (makeQAplots)
    {
        TFile *KsInvMass = new TFile((outputQAfolder_str + "/KsInvMass.root").c_str(), "RECREATE");
        TCanvas *c3 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);

        // Kshort pT and invariant mass distribution before the selections
        THnSparseF *hKshortPt = (THnSparseF *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hMassK0Shortbefore").c_str());
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
        kshortmass_before->GetXaxis()->SetTitle("m_{K_{s}} (GeV/c^{2})");
        kshortmass_before->GetYaxis()->SetTitle("Counts");
        kshortmass_before->GetXaxis()->SetRangeUser(0.2, 0.8);
        kshortpt_before->Draw("HIST");
        c3->SaveAs((outputQAfolder_str + "/kshort_pt_before." + koutputtype).c_str());
        c3->Clear();
        kshortmass_before->Draw("HIST");
        c3->SaveAs((outputQAfolder_str + "/kshort_mass_before." + koutputtype).c_str());

        // Kshort pT and invariant mass distribution after the selections
        THnSparseF *hKshortPt_after = (THnSparseF *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hMassK0ShortSelected").c_str());
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
        kshortmass_after->GetXaxis()->SetTitle("m_{K_{s}} (GeV/c^{2})");
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
        hmult->Draw();
        c3->SaveAs((outputQAfolder_str + "/hglueball_multiplicity." + koutputtype).c_str());

        // vtz distribution plot
        TH1F *hvtz = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/eventSelection/hVertexZRec").c_str());
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
        hvtz->Draw();
        c3->SaveAs((outputQAfolder_str + "/hglueball_vtz." + koutputtype).c_str());

        // mass correlation plot
        TH2F *hmasscorr = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/hglueball/hmasscorrelation").c_str());
        if (hmasscorr == nullptr)
        {
            cout << "Mass correlation plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.12, 0.05, 0.15);
        SetHistoQA(hmasscorr);
        hmasscorr->GetYaxis()->SetTitle("m_{K_{s}} (GeV/c^{2})");
        hmasscorr->GetXaxis()->SetTitle("m_{K_{s}} (GeV/c^{2})");
        hmasscorr->GetXaxis()->SetRangeUser(0.475, 0.52);
        hmasscorr->GetYaxis()->SetRangeUser(0.475, 0.52);
        hmasscorr->GetXaxis()->SetMaxDigits(3);
        hmasscorr->GetYaxis()->SetMaxDigits(3);
        hmasscorr->GetXaxis()->SetNdivisions(505);
        hmasscorr->GetYaxis()->SetNdivisions(505);
        hmasscorr->Draw("colz");
        c3->SaveAs((outputQAfolder_str + "/hglueball_masscorrelation." + koutputtype).c_str());

        // kshort selection plots
        // Armenteros alpha plot
        TH1F *hArmenteros = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/halpha").c_str());
        if (hArmenteros == nullptr)
        {
            cout << "Armenteros alpha plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hArmenteros);
        hArmenteros->GetYaxis()->SetTitle("Counts");
        hArmenteros->GetXaxis()->SetTitle("#alpha");
        hArmenteros->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_alpha." + koutputtype).c_str());

        // DCA negative daughter to PV
        TH1F *hDCAneg = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hDCAnegtopv").c_str());
        if (hDCAneg == nullptr)
        {
            cout << "DCA negative daughter to PV plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hDCAneg);
        hDCAneg->GetYaxis()->SetTitle("Counts");
        hDCAneg->GetXaxis()->SetTitle("DCA neg. daughter to PV (cm)");
        hDCAneg->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_DCAnegtopv." + koutputtype).c_str());

        // DCA positive daughter to PV
        TH1F *hDCApos = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hDCApostopv").c_str());
        if (hDCApos == nullptr)
        {
            cout << "DCA positive daughter to PV plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hDCApos);
        hDCApos->GetYaxis()->SetTitle("Counts");
        hDCApos->GetXaxis()->SetTitle("DCA pos. daughter to PV (cm)");
        hDCApos->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_DCApostopv." + koutputtype).c_str());

        // DCA daughters
        TH1F *hDCAdaughters = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hDCAV0Daughters").c_str());
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
        TH1F *hKshortLifetime = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hLT").c_str());
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
        TH2F *hNSigmaNegPion_before = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hNSigmaNegPionK0s_before").c_str());
        if (hNSigmaNegPion_before == nullptr)
        {
            cout << "n sigma neg pion daughter plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.12, 0.05, 0.15);
        SetHistoQA(hNSigmaNegPion_before);
        hNSigmaNegPion_before->GetYaxis()->SetTitle("n#sigma_{#pi^{-}}");
        hNSigmaNegPion_before->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hNSigmaNegPion_before->Draw("colz");
        c3->SaveAs((outputQAfolder_str + "/kshort_nSigmaNegPion." + koutputtype).c_str());

        // n sigma pos pion daugter before
        TH2F *hNSigmaPosPion_before = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hNSigmaPosPionK0s_before").c_str());
        if (hNSigmaPosPion_before == nullptr)
        {
            cout << "n sigma pos pion daughter plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.12, 0.05, 0.15);
        SetHistoQA(hNSigmaPosPion_before);
        hNSigmaPosPion_before->GetYaxis()->SetTitle("n#sigma_{#pi^{+}}");
        hNSigmaPosPion_before->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hNSigmaPosPion_before->Draw("colz");
        c3->SaveAs((outputQAfolder_str + "/kshort_nSigmaPosPion." + koutputtype).c_str());

        // n sigma neg pion daugter after
        TH2F *hNSigmaNegPion_after = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hNSigmaNegPionK0s_after").c_str());
        if (hNSigmaNegPion_after == nullptr)
        {
            cout << "n sigma neg pion daughter plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.12, 0.05, 0.15);
        SetHistoQA(hNSigmaNegPion_after);
        hNSigmaNegPion_after->GetYaxis()->SetTitle("n#sigma_{#pi^{-}}");
        hNSigmaNegPion_after->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hNSigmaNegPion_after->Draw("colz");
        c3->SaveAs((outputQAfolder_str + "/kshort_nSigmaNegPion_after." + koutputtype).c_str());

        // n sigma pos pion daugter after
        TH2F *hNSigmaPosPion_after = (TH2F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hNSigmaPosPionK0s_after").c_str());
        if (hNSigmaPosPion_after == nullptr)
        {
            cout << "n sigma pos pion daughter plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.12, 0.05, 0.15);
        SetHistoQA(hNSigmaPosPion_after);
        hNSigmaPosPion_after->GetYaxis()->SetTitle("n#sigma_{#pi^{+}}");
        hNSigmaPosPion_after->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hNSigmaPosPion_after->Draw("colz");
        c3->SaveAs((outputQAfolder_str + "/kshort_nSigmaPosPion_after." + koutputtype).c_str());

        // psi pair angle plot
        TH1F *hPsiPair = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hpsipair").c_str());
        if (hPsiPair == nullptr)
        {
            cout << "Psi pair angle plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hPsiPair);
        hPsiPair->GetYaxis()->SetTitle("Counts");
        hPsiPair->GetXaxis()->SetTitle("#Psi_{pair} (rad)");
        hPsiPair->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_psiPair." + koutputtype).c_str());

        // v0 cos PA
        TH1F *hV0CosPA = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hV0CosPA").c_str());
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

        // v0 radius
        TH1F *hV0Radius = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/hv0radius").c_str());
        if (hV0Radius == nullptr)
        {
            cout << "V0 radius plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hV0Radius);
        hV0Radius->GetYaxis()->SetTitle("Counts");
        hV0Radius->GetXaxis()->SetTitle("V0 radius (cm)");
        hV0Radius->GetXaxis()->SetRangeUser(0, 60);
        hV0Radius->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_v0Radius." + koutputtype).c_str());

        // negative daughter eta
        TH1F *hNegDaughterEta = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/negative_eta").c_str());
        if (hNegDaughterEta == nullptr)
        {
            cout << "Negative daughter eta plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hNegDaughterEta);
        hNegDaughterEta->GetYaxis()->SetTitle("Counts");
        hNegDaughterEta->GetXaxis()->SetTitle("Neg. daughter #eta");
        hNegDaughterEta->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_negDaughterEta." + koutputtype).c_str());

        // positive daughter eta
        TH1F *hPosDaughterEta = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/positive_eta").c_str());
        if (hPosDaughterEta == nullptr)
        {
            cout << "Positive daughter eta plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hPosDaughterEta);
        hPosDaughterEta->GetYaxis()->SetTitle("Counts");
        hPosDaughterEta->GetXaxis()->SetTitle("Pos. daughter #eta");
        hPosDaughterEta->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_posDaughterEta." + koutputtype).c_str());

        // negative daughter phi
        TH1F *hNegDaughterPhi = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/negative_phi").c_str());
        if (hNegDaughterPhi == nullptr)
        {
            cout << "Negative daughter phi plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hNegDaughterPhi);
        hNegDaughterPhi->GetYaxis()->SetTitle("Counts");
        hNegDaughterPhi->GetXaxis()->SetTitle("Neg. daughter #phi");
        hNegDaughterPhi->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_negDaughterPhi." + koutputtype).c_str());

        // positive daughter phi
        TH1F *hPosDaughterPhi = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/positive_phi").c_str());
        if (hPosDaughterPhi == nullptr)
        {
            cout << "Positive daughter phi plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hPosDaughterPhi);
        hPosDaughterPhi->GetYaxis()->SetTitle("Counts");
        hPosDaughterPhi->GetXaxis()->SetTitle("Pos. daughter #phi");
        hPosDaughterPhi->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_posDaughterPhi." + koutputtype).c_str());

        // negative daughter pT
        TH1F *hNegDaughterPt = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/negative_pt").c_str());
        if (hNegDaughterPt == nullptr)
        {
            cout << "Negative daughter pT plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hNegDaughterPt);
        hNegDaughterPt->GetYaxis()->SetTitle("Counts");
        hNegDaughterPt->GetXaxis()->SetTitle("Neg. daughter p_{T} (GeV/c)");
        hNegDaughterPt->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_negDaughterPt." + koutputtype).c_str());

        // positive daughter pT
        TH1F *hPosDaughterPt = (TH1F *)fInputFile->Get((kfoldername_temp + kvariation + "/kzeroShort/positive_pt").c_str());
        if (hPosDaughterPt == nullptr)
        {
            cout << "Positive daughter pT plot not found" << endl;
            return;
        }
        c3->Clear();
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hPosDaughterPt);
        hPosDaughterPt->GetYaxis()->SetTitle("Counts");
        hPosDaughterPt->GetXaxis()->SetTitle("Pos. daughter p_{T} (GeV/c)");
        hPosDaughterPt->Draw();
        c3->SaveAs((outputQAfolder_str + "/kshort_posDaughterPt." + koutputtype).c_str());
    }
}
