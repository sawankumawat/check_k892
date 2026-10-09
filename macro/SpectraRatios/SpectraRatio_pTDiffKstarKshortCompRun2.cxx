#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../spectra/YieldMean.C"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

Double_t FuncLavy(Double_t *x, Double_t *par)
{
    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

TFile *OpenFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (f->IsZombie())
    {
        cout << "Error: File not found: " << path << endl;
        return nullptr;
    }
    return f;
}

TGraphErrors *ScaleGraph(const TGraphErrors *gr, double factor)
{
    auto *scaled = new TGraphErrors(gr->GetN());

    for (int i = 0; i < gr->GetN(); ++i) {
        double x, y;
        gr->GetPoint(i, x, y);

        scaled->SetPoint(i, x, y * factor);
        scaled->SetPointError(i, gr->GetErrorX(i), gr->GetErrorY(i) * factor);
    }

    return scaled;
}

void SpectraRatio_pTDiffKstarKshortCompRun2()
{
    gStyle->SetOptStat(0);
    TFile *fKstar = OpenFile("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/KstarSpectra.root");
    TFile *fKshort = OpenFile("K0s_Run3_Results/Spectra_K0Short.root");

    TFile *fKstarK0sRun2 = OpenFile("ConversionCodes/pp13TeV_KstarPhi.root");
    TGraphErrors *gKstarK0sRun2_HM_stat = (TGraphErrors *)fKstarK0sRun2->Get("gKstar_KshortRatioPtClass2_stat");
    TGraphErrors *gKstarK0sRun2_HM_sys = (TGraphErrors *)fKstarK0sRun2->Get("gKstar_KshortRatioPtClass2_sys");
    TGraphErrors *gKstarK0sRun2_LM_stat = (TGraphErrors *)fKstarK0sRun2->Get("gKstar_KshortRatioPtClass10_stat");
    TGraphErrors *gKstarK0sRun2_LM_sys = (TGraphErrors *)fKstarK0sRun2->Get("gKstar_KshortRatioPtClass10_sys");
    TGraphErrors *gDoubleRatio_HM_LM_Run2_stat = (TGraphErrors *)fKstarK0sRun2->Get("gKstar_KshortDoubleRatio_Class2by10_stat");
    TGraphErrors *gDoubleRatio_HM_LM_Run2_sys = (TGraphErrors *)fKstarK0sRun2->Get("gKstar_KshortDoubleRatio_Class2by10_sys");

    // Multiplicity classes
    vector<string> multClasses = {"0_1", "1_5", "5_10", "10_15", "15_20", "20_30", "30_40", "40_50", "50_70", "70_100"};
    const int nMult = multClasses.size();

    vector<TH1F *> hKstar_stat(nMult), hKstar_sys(nMult), hKstar_sysUncorr(nMult);
    vector<TH1F *> hKshort_stat(nMult), hKshort_sys(nMult), hKshort_sysUncorr(nMult);
    vector<TF1 *> fitFcn(nMult);

    vector<TH1F *> hKstarK0sRatio_stat(nMult);
    vector<TH1F *> hKstarK0sRatio_sys(nMult);
    vector<TH1F *> hKstarK0sRatio_Uncorr(nMult);

    Double_t min = 0.0;
    Double_t max = 10;
    Double_t loprecision = 0.001;
    Double_t hiprecision = 0.5;
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 10;

    // Load histograms and initialize fits/ratios for each multiplicity class
    for (int i = 0; i < nMult; ++i)
    {
        string mStr = multClasses[i];

        hKstar_stat[i] = (TH1F *)fKstar->Get(Form("SpectraStat_%s", mStr.c_str()));
        hKstar_sys[i] = (TH1F *)fKstar->Get(Form("SpectraSys_%s", mStr.c_str()));
        hKstar_sysUncorr[i] = (TH1F *)fKstar->Get(Form("SpectraUncorrSys_%s", mStr.c_str()));

        hKshort_stat[i] = (TH1F *)fKshort->Get(Form("hK0s_%s", mStr.c_str()));
        hKshort_sys[i] = (TH1F *)fKshort->Get(Form("hK0sSyst_%s", mStr.c_str()));
        hKshort_sysUncorr[i] = (TH1F *)fKshort->Get(Form("hK0sSyst_%s", mStr.c_str())); // Using total systematic if uncorr isn't separate

        fitFcn[i] = new TF1(Form("fitfunc_%s", mStr.c_str()), FuncLavy, 0.0, 10.0, 4);
        fitFcn[i]->SetParameter(0, 6.6);
        fitFcn[i]->SetParameter(1, 0.2);
        fitFcn[i]->SetParameter(2, 0.493677);
        fitFcn[i]->SetParameter(3, 0.18);
        fitFcn[i]->SetParNames("n", "dn/dy", "mass", "T");
        fitFcn[i]->SetLineColor(kRed + 1);
        fitFcn[i]->SetLineStyle(2);
        fitFcn[i]->SetLineWidth(2);

        TH1 *hout = YieldMean(hKshort_stat[i], hKshort_sys[i], fitFcn[i], min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
        fitFcn[i]->SetNpx(500);

        hKstarK0sRatio_stat[i] = (TH1F *)hKstar_stat[i]->Clone(Form("hKstarK0sRatio_stat_%s", mStr.c_str()));
        hKstarK0sRatio_sys[i] = (TH1F *)hKstar_sys[i]->Clone(Form("hKstarK0sRatio_sys_%s", mStr.c_str()));
        hKstarK0sRatio_Uncorr[i] = (TH1F *)hKstar_sysUncorr[i]->Clone(Form("hKstarK0sRatio_Uncorr_%s", mStr.c_str()));
    }

    auto getKaonErr = [&](TH1F *hErr, double xLow, double xHigh) -> double
    {
        const double eps = 1e-6;
        const int nBins = hErr->GetNbinsX();
        const double kaonMin = hErr->GetBinLowEdge(1);
        const double kaonMax = hErr->GetBinLowEdge(nBins) + hErr->GetBinWidth(nBins);

        if (xLow < kaonMin - eps || xHigh > kaonMax + eps)
        {
            return hErr->GetBinError(nBins);
        }

        for (int ibin = 1; ibin <= nBins; ++ibin)
        {
            double kaonLow = hErr->GetBinLowEdge(ibin);
            double kaonHigh = kaonLow + hErr->GetBinWidth(ibin);

            if (std::abs(xLow - kaonLow) < eps && std::abs(xHigh - kaonHigh) < eps)
            {
                return hErr->GetBinError(ibin);
            }
        }

        double sumErr = 0.0;
        int nOverlap = 0;

        for (int ibin = 1; ibin <= nBins; ++ibin)
        {
            double kaonLow = hErr->GetBinLowEdge(ibin);
            double kaonHigh = kaonLow + hErr->GetBinWidth(ibin);

            if (kaonHigh > xLow + eps && kaonLow < xHigh - eps)
            {
                sumErr += hErr->GetBinError(ibin);
                ++nOverlap;
            }
        }

        if (nOverlap > 0)
        {
            return sumErr / nOverlap;
        }

        return hErr->GetBinError(nBins);
    };

    auto getKaonYield = [&](TH1F *hKaon, TF1 *fFit, double xLow, double xHigh) -> double
    {
        const double eps = 1e-6;
        const int nBins = hKaon->GetNbinsX();

        for (int ibin = 1; ibin <= nBins; ++ibin)
        {
            double kaonLow = hKaon->GetBinLowEdge(ibin);
            double kaonHigh = kaonLow + hKaon->GetBinWidth(ibin);

            if (std::abs(xLow - kaonLow) < eps && std::abs(xHigh - kaonHigh) < eps)
            {
                return hKaon->GetBinContent(ibin);
            }
        }

        double xCenter = 0.5 * (xLow + xHigh);
        return fFit->Eval(xCenter);
    };

    int TotalBinsKstar = hKstar_stat[0]->GetNbinsX();

    // Compute ratios for all multiplicity classes
    for (int i = 0; i < nMult; ++i)
    {
        for (int ibin = 1; ibin <= TotalBinsKstar; ++ibin)
        {
            double kstarYield = hKstar_stat[i]->GetBinContent(ibin);
            double kstar_stat_err = hKstar_stat[i]->GetBinError(ibin);
            double kstar_sys_Uncorrerr = hKstar_sysUncorr[i]->GetBinError(ibin);
            double kstar_sys_err = hKstar_sys[i]->GetBinError(ibin);

            double xLow = hKstar_stat[i]->GetBinLowEdge(ibin);
            double xHigh = xLow + hKstar_stat[i]->GetBinWidth(ibin);

            double kaonYield = getKaonYield(hKshort_stat[i], fitFcn[i], xLow, xHigh);

            if (kstarYield <= 0 || kaonYield <= 0)
            {
                hKstarK0sRatio_stat[i]->SetBinContent(ibin, 0);
                hKstarK0sRatio_stat[i]->SetBinError(ibin, 0);

                hKstarK0sRatio_sys[i]->SetBinContent(ibin, 0);
                hKstarK0sRatio_sys[i]->SetBinError(ibin, 0);

                hKstarK0sRatio_Uncorr[i]->SetBinContent(ibin, 0);
                hKstarK0sRatio_Uncorr[i]->SetBinError(ibin, 0);
                continue;
            }

            double ratio = kstarYield / kaonYield;

            // 1. Statistical uncertainty
            double kaon_stat_err = getKaonErr(hKshort_stat[i], xLow, xHigh);
            double rel_kstar_stat = kstar_stat_err / kstarYield;
            double rel_kaon_stat = kaon_stat_err / kaonYield;
            double ratio_stat_err = ratio * std::sqrt(rel_kstar_stat * rel_kstar_stat + rel_kaon_stat * rel_kaon_stat);

            // 2. Systematic uncertainty
            double kaon_sys_err = getKaonErr(hKshort_sys[i], xLow, xHigh);
            double rel_kstar_sys = kstar_sys_err / kstarYield;
            double rel_kaon_sys = kaon_sys_err / kaonYield;
            double ratio_sys_err = ratio * std::sqrt(rel_kstar_sys * rel_kstar_sys + rel_kaon_sys * rel_kaon_sys);

            // 3. Uncorrelated systematic uncertainty
            double kaon_sysUncorr_err = getKaonErr(hKshort_sysUncorr[i], xLow, xHigh);
            double rel_kstar_sysUncorr = kstar_sys_Uncorrerr / kstarYield;
            double rel_kaon_sysUncorr = kaon_sysUncorr_err / kaonYield;
            double ratio_sysUncorr_err = ratio * std::sqrt(rel_kstar_sysUncorr * rel_kstar_sysUncorr + rel_kaon_sysUncorr * rel_kaon_sysUncorr);

            // Store results
            hKstarK0sRatio_stat[i]->SetBinContent(ibin, ratio);
            hKstarK0sRatio_stat[i]->SetBinError(ibin, ratio_stat_err);

            hKstarK0sRatio_sys[i]->SetBinContent(ibin, ratio);
            hKstarK0sRatio_sys[i]->SetBinError(ibin, ratio_sys_err);

            hKstarK0sRatio_Uncorr[i]->SetBinContent(ibin, ratio);
            hKstarK0sRatio_Uncorr[i]->SetBinError(ibin, ratio_sysUncorr_err);
        }
    }

    // Access high multiplicity (first index: 0) and low multiplicity (last index: nMult - 1)
    TH1F *hKstarK0sRatioHM_stat = hKstarK0sRatio_stat[1];
    TH1F *hKstarK0sRatioHM_sys = hKstarK0sRatio_sys[1];
    TH1F *hKstarK0sRatioHM_Uncorr = hKstarK0sRatio_Uncorr[1];

    
    TH1F *hKstarK0sRatioLM_stat = hKstarK0sRatio_stat[nMult - 1];
    TH1F *hKstarK0sRatioLM_sys = hKstarK0sRatio_sys[nMult - 1];
    TH1F *hKstarK0sRatioLM_Uncorr = hKstarK0sRatio_Uncorr[nMult - 1];
    
    TGraphErrors *gKstarK0sRun2_HM_stat_scaled = ScaleGraph(gKstarK0sRun2_HM_stat, 2.0);
    TGraphErrors *gKstarK0sRun2_HM_sys_scaled = ScaleGraph(gKstarK0sRun2_HM_sys, 2.0);

    TH1D *hKstarK0sRatioHM_stat_scaled = (TH1D *)hKstarK0sRatioHM_stat->Clone("hKstarK0sRatioHM_stat_scaled");
    TH1D *hKstarK0sRatioHM_sys_scaled = (TH1D *)hKstarK0sRatioHM_sys->Clone("hKstarK0sRatioHM_sys_scaled");
    hKstarK0sRatioHM_stat_scaled->Scale(2.0);
    hKstarK0sRatioHM_sys_scaled->Scale(2.0);


    // Plotting Canvas
    TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 720, 720);
    SetCanvasStyle(cRatio, 0.16, 0.06, 0.06, 0.13);
    double pad1Size, pad2Size;
    canvas_style(cRatio, pad1Size, pad2Size);
    cRatio->cd(1);

    SetHistoQA(hKstarK0sRatioHM_stat_scaled);
    hKstarK0sRatioHM_stat_scaled->GetYaxis()->SetTitle("K^{*0} / K^{0}_{S}");
    hKstarK0sRatioHM_stat_scaled->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hKstarK0sRatioHM_stat_scaled->SetMarkerStyle(20);
    hKstarK0sRatioHM_stat_scaled->SetMarkerSize(1.45);
    hKstarK0sRatioHM_stat_scaled->GetYaxis()->SetTitleOffset(1.60 * pad1Size);
    hKstarK0sRatioHM_stat_scaled->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat_scaled->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat_scaled->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat_scaled->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hKstarK0sRatioHM_stat_scaled->GetYaxis()->SetRangeUser(-0.08, 2.35);
    hKstarK0sRatioHM_stat_scaled->GetXaxis()->SetRangeUser(0.0, 10.0);
    hKstarK0sRatioHM_stat_scaled->SetMarkerColor(kRed + 1);
    hKstarK0sRatioHM_stat_scaled->SetLineColor(kRed + 1);
    hKstarK0sRatioHM_stat_scaled->Draw("pe");
    SetHistoQA(hKstarK0sRatioHM_sys_scaled);
    hKstarK0sRatioHM_sys_scaled->SetFillStyle(0);
    hKstarK0sRatioHM_sys_scaled->SetLineColor(kRed + 1);
    hKstarK0sRatioHM_sys_scaled->SetMarkerColor(kRed + 1);
    hKstarK0sRatioHM_sys_scaled->Draw("e2 same");
    hKstarK0sRatioHM_Uncorr->SetFillStyle(1001);
    hKstarK0sRatioHM_Uncorr->SetMarkerStyle(20);
    hKstarK0sRatioHM_Uncorr->SetMarkerColor(kRed + 1);
    hKstarK0sRatioHM_Uncorr->SetFillColorAlpha(kRed + 1, 0.35);
    // hKstarK0sRatioHM_Uncorr->Draw("e2 same");

    SetHistoQA(hKstarK0sRatioLM_stat);
    SetHistoQA(hKstarK0sRatioLM_sys);
    hKstarK0sRatioLM_stat->SetMarkerStyle(21);
    hKstarK0sRatioLM_stat->SetMarkerSize(1.45);
    hKstarK0sRatioLM_stat->SetMarkerColor(kBlue + 1);
    hKstarK0sRatioLM_stat->SetLineColor(kBlue + 1);
    hKstarK0sRatioLM_stat->GetXaxis()->SetRangeUser(0.0, 10.0);
    hKstarK0sRatioLM_stat->GetYaxis()->SetRangeUser(0.08, 1.25);
    hKstarK0sRatioLM_stat->Draw("pe same");
    hKstarK0sRatioLM_sys->SetFillStyle(0);
    hKstarK0sRatioLM_sys->SetLineColor(kBlue + 1);
    hKstarK0sRatioLM_sys->SetMarkerColor(kBlue + 1);
    hKstarK0sRatioLM_sys->Draw("e2 same");
    hKstarK0sRatioLM_Uncorr->SetFillStyle(1001);
    hKstarK0sRatioLM_Uncorr->SetMarkerStyle(21);
    hKstarK0sRatioLM_Uncorr->SetMarkerColor(kBlue + 1);
    hKstarK0sRatioLM_Uncorr->SetFillColorAlpha(kBlue + 1, 0.35);
    // hKstarK0sRatioLM_Uncorr->Draw("e2 same");

    SetGraphErrorStyle(gKstarK0sRun2_HM_stat_scaled);
    SetGraphErrorStyle(gKstarK0sRun2_HM_sys_scaled);
    SetGraphErrorStyle(gKstarK0sRun2_LM_stat);
    SetGraphErrorStyle(gKstarK0sRun2_LM_sys);
    gKstarK0sRun2_HM_stat_scaled->SetMarkerColor(kMagenta);
    gKstarK0sRun2_HM_stat_scaled->SetLineColor(kMagenta);
    gKstarK0sRun2_HM_sys_scaled->SetMarkerColor(kMagenta);
    gKstarK0sRun2_HM_sys_scaled->SetLineColor(kMagenta);
    gKstarK0sRun2_LM_stat->SetMarkerColor(kGreen + 1);
    gKstarK0sRun2_LM_stat->SetLineColor(kGreen + 1);
    gKstarK0sRun2_LM_sys->SetMarkerColor(kGreen + 1);
    gKstarK0sRun2_LM_sys->SetLineColor(kGreen + 1);
    gKstarK0sRun2_HM_stat_scaled->Draw("pe same");
    gKstarK0sRun2_HM_sys_scaled->SetFillStyle(0);
    gKstarK0sRun2_HM_sys_scaled->Draw("e2 same");
    gKstarK0sRun2_LM_stat->Draw("pe same");
    gKstarK0sRun2_LM_sys->SetFillStyle(0);
    gKstarK0sRun2_LM_sys->Draw("e2 same");

    TLegend *legendRatio = new TLegend(0.6, 0.75, 0.90, 0.95);
    legendRatio->SetBorderSize(0);
    legendRatio->SetFillStyle(0);
    legendRatio->SetTextFont(42);
    legendRatio->SetTextSize(0.04);
    legendRatio->AddEntry(hKstarK0sRatioHM_stat_scaled, "FT0M: (1-5%) x 2", "p");
    legendRatio->AddEntry(hKstarK0sRatioLM_stat, "FT0M: (70-100%)", "p");
    legendRatio->AddEntry(gKstarK0sRun2_HM_stat_scaled, "Run 2 HM (1-5%) x 2", "p");
    legendRatio->AddEntry(gKstarK0sRun2_LM_stat, "Run 2 LM (70-100%)", "p");
    legendRatio->Draw();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.05);
    lat.DrawLatex(0.22, 0.89, "ALICE");
    lat.DrawLatex(0.22, 0.83, "|y| < 0.5");
    lat.DrawLatex(0.22, 0.76, "pp #sqrt{#it{s}} = 13.6 TeV");

    cRatio->cd(2);
    TH1F *hDoubleRatioStat = (TH1F *)hKstarK0sRatioHM_stat->Clone("hDoubleRatioStat");
    hDoubleRatioStat->Divide(hKstarK0sRatioLM_stat);
    TH1F *hDoubleRatioSys = (TH1F *)hKstarK0sRatioHM_Uncorr->Clone("hDoubleRatioSys");
    hDoubleRatioSys->Divide(hKstarK0sRatioLM_Uncorr);
    SetHistoQA(hDoubleRatioStat);
    hDoubleRatioStat->GetYaxis()->SetTitle("#frac{High mult.}{Low mult.}");
    hDoubleRatioStat->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hDoubleRatioStat->SetMarkerStyle(20);
    hDoubleRatioStat->SetMarkerSize(1.45);
    hDoubleRatioStat->GetYaxis()->SetTitleOffset(1.90 * pad2Size);
    hDoubleRatioStat->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    hDoubleRatioStat->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hDoubleRatioStat->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hDoubleRatioStat->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hDoubleRatioStat->GetYaxis()->SetRangeUser(0.45, 1.75);
    hDoubleRatioStat->GetXaxis()->SetRangeUser(0.0, 10.0);
    hDoubleRatioStat->SetMarkerColor(kRed + 1);
    hDoubleRatioStat->SetLineColor(kRed + 1);
    hDoubleRatioStat->GetYaxis()->SetNdivisions(505);
    hDoubleRatioStat->Draw("pe");
    SetHistoQA(hDoubleRatioSys);
    hDoubleRatioSys->SetFillStyle(0);
    hDoubleRatioSys->SetLineColor(kRed + 1);
    hDoubleRatioSys->SetMarkerColor(kRed + 1);
    hDoubleRatioSys->Draw("e2 same");

    gDoubleRatio_HM_LM_Run2_stat->SetMarkerColor(kMagenta);
    gDoubleRatio_HM_LM_Run2_stat->SetLineColor(kMagenta);
    gDoubleRatio_HM_LM_Run2_sys->SetMarkerColor(kMagenta);
    gDoubleRatio_HM_LM_Run2_sys->SetLineColor(kMagenta);
    gDoubleRatio_HM_LM_Run2_stat->Draw("pe same");
    gDoubleRatio_HM_LM_Run2_sys->SetFillStyle(0);
    gDoubleRatio_HM_LM_Run2_sys->Draw("e2 same");

    TLine *line = new TLine(0.0, 1.0, 10.0, 1.0);
    line->SetLineColor(kGray + 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    cRatio->SaveAs("Plots/KstarK0sRatio_CompareRun2.png");
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}