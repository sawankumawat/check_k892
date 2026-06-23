#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/initializations.h"
#include "spectra/YieldMean.C"

using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}
// // int colors[] = {kBlack, kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 1, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kBlue, kGray + 2};
// int colors[] = {kBlack, kBlue, kRed, kGreen + 2, kOrange + 7, kMagenta, kCyan + 2, kViolet, kOrange, kAzure + 1, kSpring + 2, kGray + 2};

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
    pad1->SetRightMargin(0.02);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.31);
    pad1->SetLeftMargin(0.12);
    pad2->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.04);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}

void plot_spectra()
{
    bool plotOnlyRaw = false;
    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);
    TString outputType = "png"; // pdf, png
    double fitRangeMax = 20.0;

    int colors[12];
    int nPaletteColors = TColor::GetNumberOfColors();

    for (int i = 0; i < 12; ++i)
    {
        int index = i * (nPaletteColors - 1) / 11;
        colors[i] = TColor::GetColorPalette(index);
    }

    //==============================Pt-dependent PID=======================
    // string path = "../output/kstar/LHC22o_pass7/586976/kstarqa_NoRCT/hInvMass"; // 2023 data
    // string path = "../output/kstar/LHC22o_pass7/586385/kstarqa/hInvMass"; // 2024 data
    TFile *fSysUncert = new TFile("../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/SystematicsPlots/SysUncert.root", "READ");
    if (fSysUncert->IsZombie())
    {
        cout << "Systematic uncertainty file not found" << endl;
        return;
    }
    TH1D *hTotalSysSmoothed = (TH1D *)fSysUncert->Get("hTotalSysSmoothed_0_100"); // Temporary assigning same to all multiplicity classes
    if (hTotalSysSmoothed == nullptr)
    {
        cout << "Histogram hTotalSysSmoothed_0_100 not found in the systematic uncertainty file" << endl;
        return;
    }

    // for (int ivar = 0; ivar < nSysVars; ivar++)
    {
        //================================After SQM=======================
        string path = "../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass"; // 2024 data
        // string path = "../output/kstar/LHC22o_pass7/682963/kstarqa_NoPVContributor/hInvMass"; // 2024 data
        // path = path + "/" + sysVars[ivar];
        TString pathLevyFits = path + "/LevyFits";
        if (gSystem->mkdir(pathLevyFits, kTRUE))
        {
            std::cout << "Folder " << pathLevyFits << " created successfully." << std::endl;
        }

        TFile *fspectra = (plotOnlyRaw) ? new TFile((path + "/yield.root").c_str(), "read") : new TFile((path + "/corrected_spectra.root").c_str(), "read");

        if (fspectra->IsZombie())
        {
            cout << "File not found" << endl;
            return;
        }
        int markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 32, 47};
        float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

        const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;
        TH1F *hmult[numofmultbins + 1];
        TH1F *hmultClone[numofmultbins + 1];

        hmult[0] = (plotOnlyRaw) ? (TH1F *)fspectra->Get("mult_0-100/yield_integral") : (TH1F *)fspectra->Get("mult_0-100/corrected_spectra_Integral_final");
        hmultClone[0] = (TH1F *)hmult[0]->Clone("hmultClone0");
        if (hmult[0] == nullptr)
        {
            cout << "Histogram 1 not found" << endl;
            return;
        }

        for (int i = 1; i < numofmultbins + 1; i++)
        {
            // hmult[i] = (TH1F *)fspectra->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", mult_classes[i - 1], mult_classes[i]));
            hmult[i] = (plotOnlyRaw) ? (TH1F *)fspectra->Get(Form("mult_%.0f-%.0f/yield_integral", mult_classes[i - 1], mult_classes[i])) : (TH1F *)fspectra->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral_final", mult_classes[i - 1], mult_classes[i]));

            if (hmult[i] == nullptr)
            {
                cout << "Histogram not found in path: " << Form("mult_%.0f-%.0f/corrected_spectra_Integral_final", mult_classes[i - 1], mult_classes[i]) << endl;
                return;
            }
            // hmult[i]->Scale(inelNormFactorRun2[i - 1]);
            hmultClone[i] = (TH1F *)hmult[i]->Clone(Form("hmultClone%d", i));
        }

        TH1F *hratios[numofmultbins];
        for (int i = 1; i < numofmultbins + 1; i++)
        {
            hratios[i] = (TH1F *)hmult[i]->Clone(Form("hratio%d", i));
            hratios[i]->Divide(hmult[0]);
        }

        TCanvas *c = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c, 0.25, 0.03, 0.03, 0.15);
        double pad1Size, pad2Size;
        canvas_style(c, pad1Size, pad2Size);
        c->cd(1);
        c->SetLogy();
        gPad->SetLogy();
        SetHistoQA(hmult[0]);
        hmult[1]->GetYaxis()->SetTitle(hmult[1]->GetYaxis()->GetTitle());
        hmult[1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        // int multiplicationFactors[] = {9, 8, 7, 5, 4, 3, 2, 1, 0}; // 6 is skipped because the multiplicity bin 10-15 and 15-20 are merged (in run 2)

        for (int i = 1; i < numofmultbins + 1; i++)
        {
            SetHistoQA(hmult[i]);
            hmult[i]->SetMarkerSize(1.2);
            hmult[i]->Scale(pow(2, numofmultbins - i));
            hmult[i]->GetXaxis()->SetTitleSize(0.045);
            hmult[i]->GetYaxis()->SetTitleSize(0.045);
            hmult[i]->GetYaxis()->SetTitleOffset(1.3);
            hmult[i]->SetMaximum(hmult[1]->GetMaximum() * 25);
            hmult[i]->SetMinimum(3e-8);
            hmult[i]->SetMarkerStyle(markers[i]);
            hmult[i]->SetLineColor(colors[i]);
            hmult[i]->SetMarkerColor(colors[i]);
            hmult[i]->Draw("pe same");
        }
        hmult[0]->SetMarkerStyle(markers[numofmultbins]);
        hmult[0]->SetMarkerSize(1.2);

        // TLegend *leg = new TLegend(0.37, 0.75, 0.95, 0.97); // for 4 columns legend
        TLegend *leg = new TLegend(0.47, 0.7, 0.98, 0.97);
        SetLegendStyle(leg);
        leg->SetNColumns(3);
        // leg->AddEntry(hmult[0], "0-100%", "lpe");
        for (int i = 1; i < numofmultbins + 1; i++)
        {
            leg->AddEntry(hmult[i], Form("%.0f-%.0f%%(#times2^{%d})", mult_classes[i - 1], mult_classes[i], numofmultbins - i), "lpe");
        }

        if (!plotOnlyRaw)
        {
            for (int imult = 0; imult < numofmultbins + 1; imult++)
            // for (int imult = 9; imult < 10; imult++)
            {
                // TH1F *h1 = (TH1F *)hmultClone[imult]->Clone("h1");
                // TH1F *h2 = (TH1F *)hmultClone[imult]->Clone("h2");
                TH1F *h1 = (TH1F *)hmult[imult]->Clone("h1");
                TH1F *h2 = (TH1F *)hmult[imult]->Clone("h2");

                for (int i = 1; i <= h2->GetNbinsX(); i++) // putting small systematic error by hand
                {
                    double systemerr = (hTotalSysSmoothed->GetBinContent(i) * h2->GetBinContent(i));
                    h2->SetBinError(i, systemerr);
                }
                /*************meanpT*****************byresonance*******************package*************************/
                Double_t min = 0.0;
                Double_t max = fitRangeMax;
                Double_t loprecision = 0.01;
                Double_t hiprecision = 0.5;
                Option_t *opt = "RI0+";
                TString logfilename = pathLevyFits + "/log_fit.root";
                Double_t minfit = 0.0;
                Double_t maxfit = fitRangeMax;

                // TF1 *fitFcn = new TF1(Form("fitfunc_%d_%d", imult, ivar), FuncLavy, 0.0, fitRangeMax, 4);
                TF1 *fitFcn = new TF1(Form("fitfunc_%d", imult), FuncLavy, 0.0, fitRangeMax, 4);
                fitFcn->SetParameter(0, 7.0);
                // fitFcn->SetParameter(1, 0.05);
                fitFcn->SetParameter(1, 0.5);
                fitFcn->FixParameter(2, 0.895);
                fitFcn->SetParameter(3, 0.35);
                fitFcn->SetParNames("n", "dn/dy", "mass", "T");
                fitFcn->SetLineColor(colors[imult]);

                TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
                c->cd(1);
                // fitFcn->SetLineColor(color);
                fitFcn->SetLineColor(colors[imult]);
                fitFcn->SetLineWidth(2);
                fitFcn->SetLineStyle(2);
                if (imult != 0)
                    fitFcn->Draw("l same");

                h2->SetFillStyle(0);
                h2->SetLineWidth(1);
                h2->SetLineColor(colors[imult]);
                h2->SetMarkerColor(colors[imult]);
                h2->SetMarkerSize(0);
                if (imult != 0)
                    h2->Draw("e2 same");

                if (imult == 0)
                    leg->AddEntry(fitFcn, "L#acute{e}vy-Tsallis", "l");
            }
        }
        c->cd(1);

        leg->SetTextSize(0.03);
        leg->Draw();

        c->cd(2);
        gPad->SetLogy(1);
        hratios[1]->SetMinimum(0.32);
        hratios[1]->GetYaxis()->SetMoreLogLabels();
        hratios[1]->GetYaxis()->SetNoExponent();
        hratios[1]->SetMaximum(8);
        hratios[1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hratios[1]->GetYaxis()->SetTitle("Ratio to INEL>0");
        hratios[1]->GetYaxis()->SetNdivisions(535);
        // hratios[0]->GetXaxis()->SetRangeUser(0, 25);

        for (int i = 1; i < numofmultbins; i++)
        {
            SetHistoQA(hratios[i]);
            hratios[i]->GetXaxis()->SetTitleSize(0.045 * pad1Size / pad2Size);
            hratios[i]->GetYaxis()->SetTitleSize(0.045 * pad1Size / pad2Size);
            hratios[i]->GetXaxis()->SetLabelSize(0.04 * pad1Size / pad2Size);
            hratios[i]->GetYaxis()->SetLabelSize(0.04 * pad1Size / pad2Size);
            hratios[i]->GetYaxis()->SetTitleOffset(0.55);
            hratios[i]->SetMarkerStyle(markers[i]);
            hratios[i]->SetLineColor(colors[i]);
            hratios[i]->SetMarkerColor(colors[i]);
            hratios[i]->Draw("pe same");
        }
        TLine *line = new TLine(0, 1, 20, 1);
        line->SetLineStyle(2);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->Draw();

        TString outputfolder = path;
        TString spectraPath = (plotOnlyRaw) ? outputfolder + "/spectra_raw." + outputType : outputfolder + "/spectra." + outputType;
        c->SaveAs(spectraPath.Data());

        if (!plotOnlyRaw)
        {
            TCanvas *ctemp = new TCanvas("ctemp", "ctemp", 720, 720);
            SetCanvasStyle(ctemp, 0.15, 0.03, 0.03, 0.15);
            gPad->SetLogy();

            double meanpT[numofmultbins], yield[numofmultbins];
            double meanpT_errStat[numofmultbins], yield_errStat[numofmultbins];
            double meanpT_errSys[numofmultbins], yield_errSys[numofmultbins];
            double nLevy[numofmultbins], nLevy_errStat[numofmultbins], nLevy_errSys[numofmultbins];
            double TLevy[numofmultbins], TLevy_errStat[numofmultbins], TLevy_errSys[numofmultbins];

            // Now lets calculate the mean pT and yield as a function of dN_charge/deta
            for (int imult = 1; imult < numofmultbins + 1; imult++)
            {
                TH1F *h1 = (TH1F *)hmultClone[imult]->Clone("h1");
                TH1F *h2 = (TH1F *)hmultClone[imult]->Clone("h2");

                for (int i = 1; i <= h2->GetNbinsX(); i++) // putting small systematic error by hand
                {
                    double systemerr = (hTotalSysSmoothed->GetBinContent(i) * h2->GetBinContent(i));
                    h2->SetBinError(i, systemerr);
                }
                /*************meanpT*****************byresonance*******************package*************************/
                Double_t min = 0.0;
                Double_t max = fitRangeMax;
                Double_t loprecision = 0.01;
                Double_t hiprecision = 0.5;
                Option_t *opt = "RI+";
                TString logfilename = pathLevyFits + "/log_mean.root";
                Double_t minfit = 0.0;
                Double_t maxfit = fitRangeMax;
                // Double_t maxfit=8.0;

                // TF1 *fitFcn = new TF1(Form("fitfunc_%d_%d", imult, ivar), FuncLavy, 0.0, fitRangeMax, 4);
                TF1 *fitFcn = new TF1(Form("fitfunc_%d", imult), FuncLavy, 0.0, fitRangeMax, 4);
                fitFcn->SetParameter(0, 5.0);
                fitFcn->SetParameter(1, 0.5);
                // fitFcn->SetParameter(1, 50);
                fitFcn->FixParameter(2, 0.895);
                fitFcn->SetParameter(3, 0.35);
                fitFcn->SetParNames("n", "dn/dy", "mass", "T");
                fitFcn->SetLineColor(colors[imult]);
                fitFcn->SetLineStyle(2);

                TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

                // Taking from Levy fit
                meanpT[imult - 1] = hout->GetBinContent(5);
                meanpT_errStat[imult - 1] = hout->GetBinContent(6);
                meanpT_errSys[imult - 1] = hout->GetBinContent(7);
                yield[imult - 1] = hout->GetBinContent(1);
                yield_errStat[imult - 1] = hout->GetBinContent(2);
                yield_errSys[imult - 1] = hout->GetBinContent(3);
                nLevy[imult - 1] = fitFcn->GetParameter(0);
                nLevy_errStat[imult - 1] = fitFcn->GetParError(0);
                TLevy[imult - 1] = fitFcn->GetParameter(3);
                TLevy_errStat[imult - 1] = fitFcn->GetParError(3);

                // //// Taking from histogram directly
                // meanpT[imult - 1] = hmultClone[imult]->GetMean();
                // meanpT_errStat[imult - 1] = hmultClone[imult]->GetMeanError();
                // yield[imult - 1] = hout->GetBinContent(1);
                // yield_errStat[imult - 1] = hout->GetBinContent(2);

                gStyle->SetOptFit(1111);
                gStyle->SetOptStat(1110);
                TCanvas *clevy = new TCanvas(Form("clevy_%d", imult), Form("clevy_%d", imult), 720, 720);
                SetCanvasStyle(clevy, 0.18, 0.03, 0.03, 0.15);
                gPad->SetLogy();
                h2->SetLineColor(colors[imult]);
                h2->SetMarkerColor(colors[imult]);
                h2->SetStats(1);
                h2->Draw();
                TString levyFitPath = pathLevyFits + Form("/levy_fit_mult_%.0f-%.0f.%s", mult_classes[imult - 1], mult_classes[imult], outputType.Data());
                clevy->SaveAs(levyFitPath.Data());
                delete clevy;

                cout << "Multiplicity class " << mult_classes[imult - 1] << " - " << mult_classes[imult] << endl;
                cout << "dN/dy: " << yield[imult - 1] << " +/- " << yield_errStat[imult - 1] << endl;
                cout << "<pT>: " << meanpT[imult - 1] << " +/- " << meanpT_errStat[imult - 1] << endl;
                cout << "\n\n";

                ctemp->cd();
                hmultClone[imult]->SetMarkerStyle(markers[imult]);
                hmultClone[imult]->SetMarkerSize(1.2);
                hmultClone[imult]->GetYaxis()->SetTitleOffset(1.3);
                hmultClone[imult]->SetMaximum(hmultClone[1]->GetMaximum() * 15);
                hmultClone[imult]->SetMinimum(2e-7);
                hmultClone[imult]->SetLineColor(colors[imult]);
                hmultClone[imult]->SetMarkerColor(colors[imult]);
                hmultClone[imult]->Draw("pe same");
                fitFcn->SetLineColor(colors[imult]);
                fitFcn->SetLineWidth(2);
                fitFcn->SetLineStyle(2);
                fitFcn->Draw("l same");
            }
            TString tempPath = outputfolder + "/spectra_temp." + outputType;
            ctemp->SaveAs(tempPath.Data());

            // double dnch_detaRun2[] = {26.18, 20.16, 16.4, 13.14, 10.3, 8.24, 6.62, 4.77, 2.76};
            double dnch_detaRun3[] = {21.78, 18.48, 15.76, 13.89, 12.50, 10.86, 9.09, 7.63, 5.87, 3.69};
            double dnch_detaRun3_err[] = {0.38, 0.25, 0.22, 0.19, 0.17, 0.15, 0.13, 0.11, 0.09, 0.06}; // (paper link: https://alice-publications.web.cern.ch/system/files/draft/10934/2025-03-03-dndeta_pp136_draft_250303.pdf)

            TGraphErrors *gMeanYieldRun3 = new TGraphErrors(numofmultbins, dnch_detaRun3, yield, dnch_detaRun3_err, yield_errStat);
            TGraphAsymmErrors *gMeanYieldRun3_sys = new TGraphAsymmErrors(numofmultbins, dnch_detaRun3, yield, dnch_detaRun3_err, dnch_detaRun3_err, yield_errSys, yield_errSys);
            TGraphErrors *gMeanpTRun3 = new TGraphErrors(numofmultbins, dnch_detaRun3, meanpT, dnch_detaRun3_err, meanpT_errStat);
            TGraphAsymmErrors *gMeanpTRun3_sys = new TGraphAsymmErrors(numofmultbins, dnch_detaRun3, meanpT, dnch_detaRun3_err, dnch_detaRun3_err, meanpT_errSys, meanpT_errSys);
            TGraphErrors *gNLevyRun3 = new TGraphErrors(numofmultbins, dnch_detaRun3, nLevy, dnch_detaRun3_err, nLevy_errStat);
            TGraphErrors *gTLevyRun3 = new TGraphErrors(numofmultbins, dnch_detaRun3, TLevy, dnch_detaRun3_err, TLevy_errStat);

            TFile *fRun2 = new TFile("spectra/pp13TeV_INELgt0.root", "read");
            TFile *fpp5020MeV = new TFile("spectra/pp5.02TeV_INELgt0.root", "read");
            if (fRun2->IsZombie() || fpp5020MeV->IsZombie())
            {
                cout << "Run 2 file not found" << endl;
                return;
            }
            TGraphErrors *gMeanYieldRun2 = (TGraphErrors *)fRun2->Get("Table 41/Graph1D_y1");
            TGraphErrors *gMeanpTRun2 = (TGraphErrors *)fRun2->Get("Table 39/Graph1D_y1");

            TH1D *hMeanYieldRun2_stat = (TH1D *)fRun2->Get("Table 41/Hist1D_y1_e1"); // stat error
            TH1D *hMeanYieldRun2_sys = (TH1D *)fRun2->Get("Table 41/Hist1D_y1_e2");  // sys error
            // TH1D *hMeanYieldRun2_sysUncorr = (TH1D *)fRun2->Get("Table 41/Hist1D_y1_e3"); // uncorrelated sys error
            TH1D *hMeanpTRun2_stat = (TH1D *)fRun2->Get("Table 39/Hist1D_y1_e1"); // stat error
            TH1D *hMeanpTRun2_sys = (TH1D *)fRun2->Get("Table 39/Hist1D_y1_e2");  // sys error
            // TH1D *hMeanpTRun2_sysUncorr = (TH1D *)fRun2->Get("Table 39/Hist1D_y1_e3"); // uncorrelated sys error

            TGraphErrors *gMeanYieldRun2_stat = new TGraphErrors(gMeanYieldRun2->GetN());
            TGraphErrors *gMeanpTRun2_stat = new TGraphErrors(gMeanpTRun2->GetN());
            TGraphErrors *gMeanYieldRun2_sys = new TGraphErrors(gMeanYieldRun2->GetN());
            TGraphErrors *gMeanpTRun2_sys = new TGraphErrors(gMeanpTRun2->GetN());

            for (int i = 0; i < gMeanYieldRun2->GetN(); ++i)
            {
                double x = 0.0, y = 0.0;
                gMeanYieldRun2->GetPoint(i, x, y);
                int bin = hMeanYieldRun2_stat->FindBin(x);
                double ex = gMeanYieldRun2->GetErrorX(i);
                double ey_stat = hMeanYieldRun2_stat->GetBinContent(bin);
                double ey_sys = hMeanYieldRun2_sys->GetBinContent(bin);
                // double ey_sysUncorr = hMeanYieldRun2_sysUncorr->GetBinContent(bin);

                gMeanYieldRun2_stat->SetPoint(i, x, y);
                gMeanYieldRun2_stat->SetPointError(i, ex, ey_stat);

                gMeanYieldRun2_sys->SetPoint(i, x, y);
                gMeanYieldRun2_sys->SetPointError(i, ex, ey_sys);

                // Similarly for mean pT
                gMeanpTRun2->GetPoint(i, x, y);
                double pt_stat = hMeanpTRun2_stat->GetBinContent(bin);
                double pt_sys = hMeanpTRun2_sys->GetBinContent(bin);
                // double pt_sysUncorr = hMeanpTRun2_sysUncorr->GetBinContent(bin);
                gMeanpTRun2_stat->SetPoint(i, x, y);
                gMeanpTRun2_stat->SetPointError(i, ex, pt_stat);

                gMeanpTRun2_sys->SetPoint(i, x, y);
                gMeanpTRun2_sys->SetPointError(i, ex, pt_sys);
            }

            TGraphErrors *gMeanYieldRun2_5020MeV = (TGraphErrors *)fpp5020MeV->Get("Table 5/Graph1D_y1");
            TGraphErrors *gMeanpTRun2_5020MeV = (TGraphErrors *)fpp5020MeV->Get("Table 6/Graph1D_y1");

            if (gMeanYieldRun2 == nullptr || gMeanpTRun2 == nullptr || gMeanYieldRun2_5020MeV == nullptr || gMeanpTRun2_5020MeV == nullptr)
            {
                cout << "Run 2 graph not found" << endl;
                return;
            }

            TFile *fLevyFit = new TFile(outputfolder + "/Results.root", "recreate");

            TCanvas *cMeanYield = new TCanvas("cMeanYield", "cMeanYield", 720, 720);
            SetCanvasStyle(cMeanYield, 0.15, 0.03, 0.03, 0.15);
            SetGraphStyle(gMeanYieldRun2_stat);
            gMeanYieldRun2_stat->SetMarkerStyle(20);
            gMeanYieldRun2_stat->SetMarkerSize(1.2);
            gMeanYieldRun2_stat->SetMarkerColor(kRed);
            gMeanYieldRun2_stat->SetLineColor(kRed);
            gMeanYieldRun2_stat->GetXaxis()->SetTitle("<#it{dN}_{ch}/d#eta>_{|#eta|< 0.5}");
            gMeanYieldRun2_stat->GetYaxis()->SetTitle("dN/dy");
            gMeanYieldRun2_stat->GetYaxis()->SetRangeUser(0.0, 0.89);
            gMeanYieldRun2_stat->SetTitle(0);
            gMeanYieldRun2_stat->Draw("AP");
            gMeanYieldRun2_stat->Write("gMeanYieldRun2_stat");
            gMeanYieldRun2_sys->SetFillStyle(0);
            gMeanYieldRun2_sys->SetLineWidth(2);
            gMeanYieldRun2_sys->SetLineColor(kRed);
            gMeanYieldRun2_sys->Draw("5 same");
            gMeanYieldRun2_sys->Write("gMeanYieldRun2_sys");

            SetGraphErrorStyle(gMeanYieldRun3);
            gMeanYieldRun3->SetMarkerStyle(21);
            gMeanYieldRun3->SetMarkerSize(1.2);
            gMeanYieldRun3->SetMarkerColor(kBlue);
            gMeanYieldRun3->SetLineColor(kBlue);
            gMeanYieldRun3->Draw("pe same");

            gMeanYieldRun3_sys->SetMarkerColor(kBlue);
            gMeanYieldRun3_sys->SetLineColor(kBlue);
            gMeanYieldRun3_sys->SetFillStyle(0);
            gMeanYieldRun3_sys->SetLineWidth(2);
            gMeanYieldRun3_sys->Draw("5 same");
            gMeanYieldRun3_sys->Write("gMeanYieldRun3_sys");
            gMeanYieldRun3->Write("gMeanYieldRun3_stat");

            gMeanYieldRun2_5020MeV->SetLineColor(kGreen + 2);
            gMeanYieldRun2_5020MeV->SetMarkerColor(kGreen + 2);
            gMeanYieldRun2_5020MeV->SetMarkerStyle(22);
            // gMeanYieldRun2_5020MeV->Draw("P same");

            TLegend *legMeanYield = new TLegend(0.2, 0.80, 0.45, 0.90);
            legMeanYield->SetTextSize(0.035);
            legMeanYield->SetBorderSize(0);
            legMeanYield->SetFillStyle(0);
            // legMeanYield->AddEntry(gMeanYieldRun2_5020MeV, "Run 2 (5.02 TeV)", "p");
            legMeanYield->SetHeader("pp collisions");
            legMeanYield->AddEntry(gMeanYieldRun3, "Run 3 (13.6 TeV)", "p");
            legMeanYield->AddEntry(gMeanYieldRun2_stat, "Run 2 (13 TeV)", "p");
            legMeanYield->Draw();

            TString yieldPath = outputfolder + "/mean_yield." + outputType;
            cMeanYield->SaveAs(yieldPath.Data());

            TCanvas *cMeanpT = new TCanvas("cMeanpT", "cMeanpT", 720, 720);
            SetCanvasStyle(cMeanpT, 0.15, 0.03, 0.03, 0.15);
            SetGraphStyle(gMeanpTRun2_stat);
            gMeanpTRun2_stat->SetMarkerStyle(20);
            gMeanpTRun2_stat->SetMarkerSize(1.2);
            gMeanpTRun2_stat->SetMarkerColor(kRed);
            gMeanpTRun2_stat->SetLineColor(kRed);
            gMeanpTRun2_stat->GetXaxis()->SetTitle("<#it{dN}_{ch}/d#eta>_{|#eta|< 0.5}");
            gMeanpTRun2_stat->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/c)");
            gMeanpTRun2_stat->GetYaxis()->SetRangeUser(0.25, 2.09);
            gMeanpTRun2_stat->SetTitle("");
            gMeanpTRun2_stat->Draw("ape");
            gMeanpTRun2_stat->Write("gMeanpTRun2_sys");
            gMeanpTRun2_sys->SetFillStyle(0);
            gMeanpTRun2_sys->SetLineWidth(2);
            gMeanpTRun2_sys->SetLineColor(kRed);
            gMeanpTRun2_sys->Draw("5 same");
            gMeanpTRun2_sys->Write("gMeanpTRun2_sys");

            SetGraphErrorStyle(gMeanpTRun3);
            gMeanpTRun3->SetMarkerStyle(21);
            gMeanpTRun3->SetMarkerSize(1.2);
            gMeanpTRun3->SetMarkerColor(kBlue);
            gMeanpTRun3->SetLineColor(kBlue);
            gMeanpTRun3->Draw("pe same");

            gMeanpTRun3_sys->SetLineColor(kBlue);
            gMeanpTRun3_sys->SetFillStyle(0);
            gMeanpTRun3_sys->SetLineWidth(2);
            gMeanpTRun3_sys->Draw("5 same");
            gMeanpTRun3_sys->Write("gMeanpTRun3_sys");

            gMeanpTRun2_5020MeV->SetLineColor(kGreen + 2);
            gMeanpTRun2_5020MeV->SetMarkerColor(kGreen + 2);
            gMeanpTRun2_5020MeV->SetMarkerStyle(22);
            // gMeanpTRun2_5020MeV->Draw("P same");
            legMeanYield->Draw();
            gMeanpTRun3->Write("gMeanpTRun3_stat");
            TString ptPath = outputfolder + "/mean_pT." + outputType;
            cMeanpT->SaveAs(ptPath.Data());

            TCanvas *cNLevy = new TCanvas("cNLevy", "cNLevy", 720, 720);
            SetCanvasStyle(cNLevy, 0.15, 0.03, 0.03, 0.15);
            SetGraphErrorStyle(gNLevyRun3);
            gNLevyRun3->SetMarkerStyle(21);
            gNLevyRun3->SetMarkerSize(1.2);
            gNLevyRun3->SetMarkerColor(kBlue);
            gNLevyRun3->SetLineColor(kBlue);
            gNLevyRun3->GetXaxis()->SetTitle("<#it{dN}_{ch}/d#eta>_{|#eta|< 0.5}");
            gNLevyRun3->GetYaxis()->SetTitle("n");
            gNLevyRun3->GetYaxis()->SetRangeUser(6.2, 9.5);
            gNLevyRun3->SetTitle("");
            gNLevyRun3->Draw("ape");
            gNLevyRun3->Write("gNLevyRun3_stat");
            TString nLevyPath = outputfolder + "/nLevy." + outputType;
            cNLevy->SaveAs(nLevyPath.Data());

            TCanvas *cTLevy = new TCanvas("cTLevy", "cTLevy", 720, 720);
            SetCanvasStyle(cTLevy, 0.15, 0.03, 0.03, 0.15);
            SetGraphErrorStyle(gTLevyRun3);
            gTLevyRun3->SetMarkerStyle(21);
            gTLevyRun3->SetMarkerSize(1.2);
            gTLevyRun3->SetMarkerColor(kBlue);
            gTLevyRun3->SetLineColor(kBlue);
            gTLevyRun3->GetXaxis()->SetTitle("<#it{dN}_{ch}/d#eta>_{|#eta|< 0.5}");
            gTLevyRun3->GetYaxis()->SetTitle("T (GeV)");
            gTLevyRun3->GetYaxis()->SetRangeUser(0.13, 0.48);
            gTLevyRun3->SetTitle("");
            gTLevyRun3->Draw("ape");
            gTLevyRun3->Write("gTLevyRun3_stat");
            TString TLevyPath = outputfolder + "/TLevy." + outputType;
            cTLevy->SaveAs(TLevyPath.Data());
        }
    }
} // End of the code

// IR study
//  string path = "../output/kstar/LHC22o_pass7/IR_study/463114/kstarqa/hInvMass"; // 1-2 MHz
//  string path = "../output/kstar/LHC22o_pass7/IR_study/535069/kstarqa/hInvMass"; // 14 kHz
//  string path = "../output/kstar/LHC22o_pass7/IR_study/535545/kstarqa/hInvMass"; // 70 kHz
//  string path = "../output/kstar/LHC22o_pass7/IR_study/535645/kstarqa/hInvMass"; // 135 kHz
//  string path = "../output/kstar/LHC22o_pass7/IR_study/LHC23z/kstarqa/hInvMass"; // 450 kHz
//  string path = "../output/kstar/LHC22o_pass7/IR_study/LHC23ls/kstarqa/hInvMass"; // 650 kHz
// string path = "../output/kstar/LHC22o_pass7/IR_study/466180/kstarqa_id33593/hInvMass"; // 2024 data (500 kHz)
// string path = "../output/kstar/LHC22o_pass7/IR_study/459908/kstarqa/hInvMass"; // 2023 data (135 kHz)
// string path = "../output/kstar/LHC22o_pass7/IR_study/459845/kstarqa/hInvMass"; // 2022 data (500 kHz)

// Checks on the data
// string path = "../output/kstar/LHC22o_pass7/checks/473185/kstarqa_MID_id33593/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

////*************************PID Variations for Kaon (without MID, multcentTable)**************************
// string path = "../output/kstar/LHC22o_pass7/999999/kstarqa/hInvMass"; // 2022 data
// string path = "../output/kstar/LHC22o_pass7/480448/kstarqa/hInvMass"; // 2023 data
// string path = "../output/kstar/LHC22o_pass7/480358/kstarqa/hInvMass"; // 2024 data

//*************************ItsTpcTracksCheck, betacutTOF******************************
// string path = "../output/kstar/LHC22o_pass7/481941/kstarqa_PIDKa1_itstpc/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0

// string sysVars[] = {"", "Norm1", "Norm2", "FitRange1", "FitRange2", "WidhtFree"};
// TString sysVars[] = {""};
// int nSysVars = sizeof(sysVars) / sizeof(sysVars[0]);

// double inelNormFactorRun2[] = {0.997814, 0.998632, 0.998465, 0.997509, 0.993852, 0.985782, 0.971972, 0.935197, 0.756786}; // this is event loss factor used in run 2

// ******************Correct placement of TPC crossed rows**************************
// string path = "../output/kstar/LHC22o_pass7/459845/kstarqa/hInvMass"; // 2022 data
// string path = "../output/kstar/LHC22o_pass7/459908/kstarqa_PIDKa2/hInvMass"; // 2023 data
// string path = "../output/kstar/LHC22o_pass7/460233/kstarqa_PIDKa2/hInvMass"; // 2024 data

//*********************PID Variations for Kaon (without MID)************************
// string path = "../output/kstar/LHC22o_pass7/480317/kstarqa/hInvMass"; // 2022 data
// string path = "../output/kstar/LHC22o_pass7/480447/kstarqa/hInvMass"; // 2023 data
// string path = "../output/kstar/LHC22o_pass7/480657/kstarqa/hInvMass"; // 2024 data
