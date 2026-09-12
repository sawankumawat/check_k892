#include <iostream>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "YieldMean.C"
using namespace std;

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void compare_yield_inelgt0_minBias()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    bool isSameBins = false;

    TFile *fSystematics = new TFile("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots/SysUncert.root", "READ");
    if (fSystematics->IsZombie())
    {
        cout << "Error: Systematics file not found" << endl;
        return;
    }

    TH1D *hRelUncertMultEstVars = (TH1D *)fSystematics->Get("hTotalSysSmoothed_0_100");
    if (hRelUncertMultEstVars == nullptr)
    {
        cout << "Error: Relative uncertainty histogram not found in systematics file" << endl;
        return;
    }

    string path1 = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED"; // 2024 data
    TString outputPath = path1 + "/spectra_compare";
    gSystem->mkdir(outputPath, kTRUE);

    TFile *fspectra1 = new TFile((path1 + "/corrected_spectra_0_100.root").c_str(), "read");

    if (fspectra1->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    TFile *fpub = new TFile("pp13TeV_INELgt0.root", "READ");
    if (fpub->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;
    cout << "number of multiplicity bins are " << numofmultbins << endl;

    TH1F *hmult1;
    TH1F *hmultClone1;
    TGraphErrors *gRun2_spectra[numofmultbins];
    TGraphErrors *gRun2_ratio[numofmultbins];
    TGraphErrors *gRun2_minBias[numofmultbins];

    hmult1 = (TH1F *)fspectra1->Get("mult_0-100/corrected_spectra_Integral_final");
    hmultClone1 = (TH1F *)hmult1->Clone("hmultClone0");

    if (hmult1 == nullptr)
    {
        cout << "Histogram 1 not found" << endl;
        return;
    }

    // It is seen that apart from 0-1 multiplicity class, all other give same min bias (0-100) yield.

    for (int imult = 0; imult < numofmultbins; imult++)
    {
        gRun2_spectra[imult] = (TGraphErrors *)fpub->Get(Form("Table %d/Graph1D_y1", imult + 1));
        gRun2_minBias[imult] = (TGraphErrors *)fpub->Get(Form("Table %d/Graph1D_y1", imult + 1));
        if (gRun2_spectra[imult] == nullptr)
        {
            cout << "Run2 yield graph not found for mult bin " << imult << endl;
            return;
        }
        gRun2_ratio[imult] = (TGraphErrors *)fpub->Get(Form("Table %d/Graph1D_y1", imult + 1 + 9)); // ratio of given multiplicity class to 0-100% class
        if (gRun2_ratio[imult] == nullptr)
        {
            cout << "Run2 ratio graph not found for mult bin " << imult << endl;
            return;
        }

        int numPoints = gRun2_minBias[imult]->GetN();
        for (int i = 0; i < numPoints; i++)
        {
            double x, ymult, yratio, xerror, yerror;
            gRun2_minBias[imult]->GetPoint(i, x, ymult);
            gRun2_ratio[imult]->GetPoint(i, x, yratio);
            xerror = gRun2_minBias[imult]->GetErrorX(i);
            yerror = gRun2_minBias[imult]->GetErrorY(i);
            // ymult /= inelNormRun2;

            // Now the yminBias is ratio of gRun2_minBias to gRun2_ratio
            double minBiasYield = ymult / yratio;
            double ymult_error = gRun2_minBias[imult]->GetErrorY(i);
            double yratio_error = gRun2_ratio[imult]->GetErrorY(i);
            double error = sqrt(pow(ymult_error / yratio, 2) + pow(ymult * yratio_error / (yratio * yratio), 2));
            gRun2_minBias[imult]->SetPoint(i, x, minBiasYield);
            gRun2_minBias[imult]->SetPointError(i, 0, error);
        }
    }

    TH1F *h1 = (TH1F *)hmultClone1->Clone("h1");
    TH1F *h2 = (TH1F *)hmultClone1->Clone("h2");

    for (int i = 1; i <= h1->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr1 = (hRelUncertMultEstVars->GetBinContent(i) * h2->GetBinContent(i));
        h2->SetBinError(i, systemerr1);
    }

    TF1 *fitFcn1 = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
    fitFcn1->SetParameter(0, 5.0);
    fitFcn1->SetParameter(1, 0.07);
    fitFcn1->FixParameter(2, 0.895);
    fitFcn1->SetParameter(3, 0.3);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");
    fitFcn1->SetLineColor(kBlack);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);

    /*************meanpT*****************byresonance*******************package*************************/
    Double_t min = 0;
    Double_t max = 10;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 10;
    // Double_t maxfit=8.0;

    TH1 *hout = YieldMean(h1, h2, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    TGraphErrors *gratio1 = new TGraphErrors();
    TGraphErrors *gratio1_sys = new TGraphErrors();
    int minBiasFromWhichGraph = 7;

    if (!isSameBins)
    {
        for (int i = 0; i < gRun2_minBias[minBiasFromWhichGraph]->GetN(); i++)
        {
            double x_run2, yield_run2, x_error, yieldError_run2;
            gRun2_minBias[minBiasFromWhichGraph]->GetPoint(i, x_run2, yield_run2);
            x_error = gRun2_minBias[minBiasFromWhichGraph]->GetErrorX(i);
            yieldError_run2 = gRun2_minBias[minBiasFromWhichGraph]->GetErrorY(i);

            if (i == 6 || i == 7)
                yield_run2 = yield_run2 * 0.92; // event loss for 5-10 and 10-15 multiplicity classes
            if (i == 8)
                yield_run2 = yield_run2 * 0.95; // event loss for 1-5 multiplicity class

            // 1. Evaluate Levy fit at Run 2 x-position
            double thisanalysis1 = fitFcn1->Eval(x_run2);

            // 2. Look up relative systematic uncertainty corresponding to x_run2
            int sysBin = hRelUncertMultEstVars->GetXaxis()->FindBin(x_run2);
            double relSysError = hRelUncertMultEstVars->GetBinContent(sysBin);
            double sysError_fit = thisanalysis1 * relSysError;

            // 3. Statistical Error Propagation for gratio1
            //    (Assuming yieldError_run2 contains stat error)
            double error1_stat = sqrt(pow(thisanalysis1 * yieldError_run2 / (yield_run2 * yield_run2), 2));

            // 4. Full Systematic Error Propagation for gratio1_sys
            double error1_sys = sqrt(pow(sysError_fit / yield_run2, 2) +
                                     pow(thisanalysis1 * yieldError_run2 / (yield_run2 * yield_run2), 2));

            // Set central points and error bands
            gratio1->SetPoint(i, x_run2, thisanalysis1 / yield_run2);
            gratio1->SetPointError(i, x_error, error1_stat);

            gratio1_sys->SetPoint(i, x_run2, thisanalysis1 / yield_run2);

            // Set x-width for the systematic error box (0.10 GeV/c default or bin width)
            double xBand = (x_error > 0) ? x_error : 0.10;
            gratio1_sys->SetPointError(i, xBand, error1_sys);

            cout << "Bin " << i << " at pT=" << x_run2
                 << " | Ratio = " << thisanalysis1 / yield_run2
                 << " | SysError = " << error1_sys << endl;
        }
    }
    else
    {
        if (gRun2_minBias[minBiasFromWhichGraph]->GetN() != hmult1->GetNbinsX())
        {
            cout << "Error: Number of points in Run2 minBias graph does not match number of bins in histogram." << endl;
            cout << "Number of points in graphs is " << gRun2_minBias[minBiasFromWhichGraph]->GetN() << endl;
            cout << "Number of bins in histogram is " << hmult1->GetNbinsX() << endl;
            return;
        }
        cout << "Number of points in Run2 minBias graph is " << gRun2_minBias[minBiasFromWhichGraph]->GetN() << endl;
        for (int i = 0; i < gRun2_minBias[minBiasFromWhichGraph]->GetN(); i++)
        {
            double x_run2, yield_run2, x_error, yieldError_run2;
            gRun2_minBias[minBiasFromWhichGraph]->GetPoint(i, x_run2, yield_run2);

            gRun2_minBias[minBiasFromWhichGraph]->SetPoint(i, x_run2, yield_run2);
            x_error = gRun2_minBias[minBiasFromWhichGraph]->GetErrorX(i);
            yieldError_run2 = gRun2_minBias[minBiasFromWhichGraph]->GetErrorY(i);
            double yieldError_run3 = h1->GetBinError(i + 1);
            double SysError_run3 = h2->GetBinError(i + 1);

            double yield_Run3 = hmultClone1->GetBinContent(i + 1);

            gratio1->SetPoint(i, x_run2, yield_Run3 / yield_run2);
            gratio1_sys->SetPoint(i, x_run2, yield_Run3 / yield_run2);

            double error1 = sqrt(pow(yieldError_run3 / yield_run2, 2) + pow(yield_Run3 * yieldError_run2 / (yield_run2 * yield_run2), 2));

            double error1_sys = sqrt(pow(SysError_run3 / yield_run2, 2) + pow(yield_Run3 * yieldError_run2 / (yield_run2 * yield_run2), 2));

            double sysError = hRelUncertMultEstVars->GetBinContent(i + 1); // systematic error from histogram
            gratio1->SetPointError(i, x_error, error1);

            double xBand = 0.5 * hmultClone1->GetXaxis()->GetBinWidth(i + 1);
            gratio1_sys->SetPointError(i, xBand, error1_sys);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    SetHistoStyle(hmultClone1, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    hmultClone1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hmultClone1->SetMaximum(hmultClone1->GetMaximum() * 3);
    hmultClone1->SetMinimum(hmultClone1->GetMinimum() * 0.9);
    hmultClone1->GetYaxis()->SetTitleOffset(1.30);
    hmultClone1->GetXaxis()->SetTitleOffset(1.02);
    hmultClone1->SetMarkerStyle(20);
    hmultClone1->SetMarkerSize(1);
    hmultClone1->GetXaxis()->SetRangeUser(0, 10);
    hmultClone1->SetLineColor(kBlue);
    hmultClone1->SetMarkerColor(kBlue);
    hmultClone1->Draw("pe");
    h2->SetMarkerColor(kBlue);
    h2->SetLineColor(kBlue);
    h2->SetFillStyle(0);
    // h2->SetFillColorAlpha(kBlue, 0.20);
    h2->SetLineWidth(1);
    h2->Draw("e2 same");

    fitFcn1->SetLineColor(kBlue);
    fitFcn1->SetLineStyle(2);
    fitFcn1->Draw("same");

    gPad->SetLogy(1);
    gRun2_minBias[minBiasFromWhichGraph]->SetMarkerStyle(22);
    gRun2_minBias[minBiasFromWhichGraph]->SetMarkerSize(1);
    gRun2_minBias[minBiasFromWhichGraph]->SetMarkerColor(kBlack);
    gRun2_minBias[minBiasFromWhichGraph]->SetLineColor(kBlack);
    gRun2_minBias[minBiasFromWhichGraph]->SetLineWidth(2);
    gRun2_minBias[minBiasFromWhichGraph]->Draw("pe same");

    TLegend *leg = new TLegend(0.233983, 0.066092, 0.614206, 0.265189);
    SetLegendStyle(leg);
    leg->SetTextSize(0.045);
    leg->AddEntry(hmultClone1, "pp 13.6 TeV", "p");
    leg->AddEntry(gRun2_minBias[minBiasFromWhichGraph], "pp 13 TeV (#it{PLB 807 (2020) 135501)}", "p");
    leg->AddEntry(fitFcn1, "Levy-Tsallis", "l");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.SetTextFont(42);
    lat.DrawLatex(0.65, 0.90, "ALICE");
    lat.DrawLatex(0.65, 0.83, "|y| < 0.5, INEL > 0");
    lat.DrawLatex(0.65, 0.76, "FT0M: 0-100%");
    lat.DrawLatex(0.65, 0.69, "K*(892)^{0}");


    c1->cd(2);
    SetGraphErrorStyle(gratio1);
    gratio1->GetYaxis()->SetTitleSize(0.035 / pad2Size);
    gratio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    gratio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    gratio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    gratio1->SetMarkerStyle(20);
    gratio1->SetMarkerSize(1.0);
    gratio1->SetMarkerColor(kBlue);
    gratio1->SetLineColor(kBlue);
    gratio1->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
    gratio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gratio1->GetXaxis()->CenterTitle(1);
    gratio1->GetYaxis()->SetTitleOffset(0.6);
    gratio1->GetXaxis()->SetTitleOffset(1.1);
    gratio1->GetYaxis()->SetNdivisions(506);
    gratio1->GetXaxis()->SetRangeUser(0, 10);
    gratio1->GetHistogram()->SetMaximum(1.45);
    gratio1->GetHistogram()->SetMinimum(0.65);
    gratio1->Draw("ap");
    gratio1_sys->SetMarkerStyle(20);
    gratio1_sys->SetMarkerSize(1.0);
    gratio1_sys->SetMarkerColor(kBlue);
    gratio1_sys->SetLineColor(kBlue);
    gratio1_sys->SetFillColor(kBlue);
    gratio1_sys->SetLineWidth(2);
    gratio1_sys->SetFillStyle(0);
    gratio1_sys->Draw("e2 same");

    TLine *line = new TLine(0, 1, 10, 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(1);
    line->Draw();

    // Draw a grey band with 20% uncertainty around the ratio of 1
    TBox *box = new TBox(0, 0.9, 10, 1.1);
    box->SetFillColor(kGray + 2);
    box->SetFillStyle(3003);
    box->Draw("same");
    c1->SaveAs(outputPath + "/YieldMinBiasRatio.pdf");
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