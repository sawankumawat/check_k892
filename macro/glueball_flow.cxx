#include <iostream>
#include <cmath>
#include <TKey.h>
#include <TClass.h>
#include <TDirectory.h>
#include "TArrow.h"
#include "src/style.h"

void glueball_flow()
{
    gStyle->SetOptStat(1110);
    bool makeQAplots = false;
    TString dataPath = "../data/glueball/LHC22o_pass7_small";
    TString fileName = "699196";
    TString outputfolder = "../output/glueball/Flow/" + fileName;

    if (gSystem->mkdir(outputfolder, kTRUE))
    {
        std::cout << "Folder " << outputfolder << " created successfully." << std::endl;
    }

    TFile *fInputFile = new TFile((dataPath + "/" + fileName + ".root"), "Read");
    if (fInputFile->IsZombie())
    {
        cerr << "File not found " << endl;
        return;
    }

    TString rootFolderName = "higher-mass-resonances";

    //========Multiplicity histogram========
    TH1F *hMult = (TH1F *)fInputFile->Get((rootFolderName + "/eventSelection/hmultiplicity"));
    if (hMult == nullptr)
    {
        cout << "Multiplicity histogram not found" << endl;
        return;
    }
    //=========Multiplicity and pT range for analysis========
    int multLow = 0;
    int multHigh = 30;

    int ptLow = 2;
    int ptHigh = 30;

    int realEvents = hMult->Integral(hMult->GetXaxis()->FindBin(multLow), hMult->GetXaxis()->FindBin(multHigh));

    //========Invariant mass histograms for signal and background========
    THnSparseF *fHistNum = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassEPDS", rootFolderName.Data()));
    THnSparseF *fHistRot = (THnSparseF *)fInputFile->Get(Form("%s/hglueball/h3glueInvMassEPRot", rootFolderName.Data()));
    if (fHistNum == nullptr || fHistRot == nullptr)
    {
        cout << "Invariant mass histograms not found" << endl;
        return;
    }

    //========Find the bin range for the selected multiplicity and pT range========
    int lbinmult = fHistNum->GetAxis(0)->FindBin(multLow + 1e-3);
    int hbinmult = fHistNum->GetAxis(0)->FindBin(multHigh - 1e-3);

    int lbinpt = fHistNum->GetAxis(1)->FindBin(ptLow + 1e-3);
    int hbinpt = fHistNum->GetAxis(1)->FindBin(ptHigh - 1e-3);

    fHistNum->GetAxis(0)->SetRange(lbinmult, hbinmult);
    fHistRot->GetAxis(0)->SetRange(lbinmult, hbinmult);

    fHistNum->GetAxis(1)->SetRange(lbinpt, hbinpt);
    fHistRot->GetAxis(1)->SetRange(lbinpt, hbinpt);

    const int TotalAngBins = fHistNum->GetAxis(3)->GetNbins();
    cout << "Total number of angular bins is " << TotalAngBins << endl;
    TH1F *hInvMass[TotalAngBins];
    TH1F *hInvMassRot[TotalAngBins];
    TH1F *hInvMassSubtracted[TotalAngBins];

    //========Rebinning factor and normalization range========
    int rebinFactor = 1;
    double energyBinWidth = fHistNum->GetAxis(2)->GetBinWidth(2) * rebinFactor;
    cout << "Energy bin width is " << energyBinWidth << " GeV/c^2" << endl;
    double normLow = 2.3;
    double normHigh = 2.6;

    //===========Output root file==================
    TFile *fOutput = new TFile((outputfolder + "/glueball_flow.root"), "RECREATE");
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.SetTextFont(42);

    //======Projection of invariant mass for different v2 bins========
    for (int ibins = 0; ibins < TotalAngBins; ibins++)
    {
        int lbinang = fHistNum->GetAxis(3)->FindBin(fHistNum->GetAxis(3)->GetBinLowEdge(ibins) + 1e-3);
        int hbinang = fHistNum->GetAxis(3)->FindBin(fHistNum->GetAxis(3)->GetBinUpEdge(ibins) - 1e-3);

        fHistNum->GetAxis(3)->SetRange(lbinang, hbinang);
        fHistRot->GetAxis(3)->SetRange(lbinang, hbinang);

        hInvMass[ibins] = (TH1F *)fHistNum->Projection(2, "E");
        hInvMassRot[ibins] = (TH1F *)fHistRot->Projection(2, "E");
        hInvMass[ibins]->SetName(Form("hInvMass_%d", ibins));
        hInvMassRot[ibins]->SetName(Form("hInvMassRot_%d", ibins));

        int normBinLow = hInvMass[ibins]->GetXaxis()->FindBin(normLow + 1e-3);
        int normBinHigh = hInvMass[ibins]->GetXaxis()->FindBin(normHigh - 1e-3);
        double sigIntegral = hInvMass[ibins]->Integral(normBinLow, normBinHigh);
        double bkgIntegral = hInvMassRot[ibins]->Integral(normBinLow, normBinHigh);
        double scaleFactor = sigIntegral / bkgIntegral;
        TH1F *hRotScaled = (TH1F *)hInvMassRot[ibins]->Clone(Form("hRotScaled_%d", ibins));
        hRotScaled->Scale(scaleFactor);

        hInvMassSubtracted[ibins] = (TH1F *)hInvMass[ibins]->Clone(Form("hInvMassSubtracted_%d", ibins));
        hInvMassSubtracted[ibins]->Add(hRotScaled, -1);
        hInvMassSubtracted[ibins]->Rebin(rebinFactor);

        hInvMass[ibins]->Write();
        hInvMassRot[ibins]->Write();
        SetHistoQA(hInvMassSubtracted[ibins]);
        hInvMassSubtracted[ibins]->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", energyBinWidth * 1000));
        hInvMassSubtracted[ibins]->SetMarkerSize(0.5);
        hInvMassSubtracted[ibins]->Write();

        TCanvas *cShowNormRange = new TCanvas(Form("cShowNormRange_%d", ibins), "", 720, 720);
        SetCanvasStyle(cShowNormRange, 0.15, 0.03, 0.05, 0.15);
        TH1F *hbkg_nopeak = (TH1F *)hRotScaled->Clone(Form("hbkg_nopeak_%d", ibins));
        hbkg_nopeak->SetLineColor(kBlue - 7);
        hbkg_nopeak->SetMarkerColor(kBlue - 7);
        hbkg_nopeak->SetFillColor(kBlue - 7);
        for (int i = 0; i < hbkg_nopeak->GetNbinsX(); i++)
        {
            if (hbkg_nopeak->GetBinCenter(i + 1) < normLow || hbkg_nopeak->GetBinCenter(i + 1) > normHigh)
            {
                hbkg_nopeak->SetBinContent(i + 1, -999);
            }
        }

        SetHistoQA(hInvMass[ibins]);
        SetHistoQA(hRotScaled);
        hInvMass[ibins]->SetMarkerStyle(20);
        hInvMass[ibins]->SetMarkerColor(kBlack);
        hInvMass[ibins]->SetMarkerSize(0.8);
        hRotScaled->SetMarkerStyle(20);
        hRotScaled->SetMarkerSize(0.8);
        hRotScaled->SetMarkerColor(kRed);
        hRotScaled->SetLineColor(kRed);
        hInvMass[ibins]->GetYaxis()->SetMaxDigits(3);
        hInvMass[ibins]->GetYaxis()->SetTitleOffset(1.5);
        hInvMass[ibins]->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", energyBinWidth * 1000));
        hInvMass[ibins]->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
        hInvMass[ibins]->Draw("E");
        hRotScaled->Draw("E same");
        hbkg_nopeak->Draw("BAR same");

        TLegend *leg = new TLegend(0.25, 0.2454598, 0.5445682, 0.3908046);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->AddEntry(hInvMass[ibins], "Same-event pairs", "p");
        string bkgname = "Rotated pairs";
        leg->AddEntry(hRotScaled, bkgname.c_str(), "p");
        hbkg_nopeak->SetLineWidth(0);
        leg->AddEntry(hbkg_nopeak, "Normalization region", "f");
        leg->Draw();
        lat.DrawLatex(0.4, 0.9, Form("EP bin %d", ibins + 1));
        lat.DrawLatex(0.4, 0.85, Form("Multiplicity: %d-%d%%", multLow, multHigh));
        cShowNormRange->Write(Form("cShowNormRange_%d", ibins));
        cShowNormRange->SaveAs((outputfolder + Form("/cShowNormRange_%d.png", ibins)).Data());

        TCanvas *cSignal = new TCanvas(Form("cSignal_%d", ibins), "", 720, 720);
        SetCanvasStyle(cSignal, 0.15, 0.03, 0.05, 0.15);
        TH1F *hSignal = (TH1F *)hInvMassSubtracted[ibins]->Clone(Form("hSignal_%d", ibins));
        int rebinForSignal = 4;
        SetHistoQA(hSignal);
        hSignal->Rebin(rebinForSignal);
        hSignal->GetYaxis()->SetMaxDigits(3);
        hSignal->GetYaxis()->SetTitleOffset(1.5);
        hSignal->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", energyBinWidth * rebinForSignal * 1000));
        hSignal->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
        hSignal->GetXaxis()->SetRangeUser(1.0, 2.4);
        hSignal->Draw("E");
        lat.DrawLatex(0.4, 0.9, Form("EP bin %d", ibins + 1));
        lat.DrawLatex(0.4, 0.85, Form("Multiplicity: %d-%d%%", multLow, multHigh));
        lat.DrawLatex(0.4, 0.8, Form("p_{T} range: %d-%d GeV/c", ptLow, ptHigh));
        cSignal->Write(Form("cSignal_%d", ibins));
        cSignal->SaveAs((outputfolder + Form("/cSignal_%d.png", ibins)).Data());
    }

    //==========Some QA plots==========================
    if (makeQAplots)
    {
        if (gSystem->mkdir((outputfolder + "/QAplots").Data(), kTRUE))
        {
            std::cout << "Folder " << outputfolder + "/QAplots" << " created successfully." << std::endl;
        }
        string outputQAfolder_str = (outputfolder + "/QAplots").Data();
        string foldername_final = rootFolderName.Data();
        string koutputtype = "png";

        TCanvas *c3 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c3, 0.15, 0.03, 0.05, 0.15);

        TH2F *EpDet = (TH2F *)fInputFile->Get(rootFolderName + "/hglueball/EpDet");
        TH2F *EpRefA = (TH2F *)fInputFile->Get(rootFolderName + "/hglueball/EpRefA");
        TH2F *EPRefB = (TH2F *)fInputFile->Get(rootFolderName + "/hglueball/EpRefB");
        TH2F *EpResDetRefA = (TH2F *)fInputFile->Get(rootFolderName + "/hglueball/EpResDetRefA");
        TH2F *EpResDetRefB = (TH2F *)fInputFile->Get(rootFolderName + "/hglueball/EpResDetRefB");
        TH2F *EpResRefARefB = (TH2F *)fInputFile->Get(rootFolderName + "/hglueball/EpResRefARefB");
        if (EpDet == nullptr || EpRefA == nullptr || EPRefB == nullptr || EpResDetRefA == nullptr || EpResDetRefB == nullptr || EpResRefARefB == nullptr)
        {
            cout << "Event plane histograms not found" << endl;
            return;
        }

        c3->cd();

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
        SetHistoQA(hMult);
        hMult->GetYaxis()->SetTitle("Counts");
        hMult->GetXaxis()->SetTitle("Multiplicity percentile");
        hMult->GetXaxis()->SetRangeUser(0, 105);
        hMult->SetStats(0);
        hMult->Draw();
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
        gPad->SetLogy();
        hNSigmaNegPion_before->SetStats(0);
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

        gPad->SetLogx(0);
        gPad->SetLogz(0);
        gPad->SetLogy(0);

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

        // Number of ks produced in an event plot
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
        hNofKshort->GetYaxis()->SetTitle("Counts");
        hNofKshort->GetXaxis()->SetTitle("Number of K^{0}_{s} produced");
        hNofKshort->GetYaxis()->SetMaxDigits(3);
        hNofKshort->SetMaximum(hNofKshort->GetMaximum() * 100);
        hNofKshort->Draw("HIST text");
        c3->SaveAs((outputQAfolder_str + "/NKs_produced." + koutputtype).c_str());
    }
}