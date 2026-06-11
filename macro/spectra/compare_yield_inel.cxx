#include <iostream>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/initializations.h"
#include "YieldMean.C"

using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
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

void compare_yield_inel()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    bool isSameBins = true;
    bool removeNormFactorsRun2 = false;
    bool compareUsingINELgt0SigLoss = true; // Best comparison is using Event loss/signal loss with for INEL>0 method. (This is used in Bong-Hwi analysis note. so we can use it.)

    ////Normalization factors used in run 2 INEL study:
    double f_norm_run2 = 0.7448;
    double f_vtx_run2 = 0.931264; // (In Phi AN, f_vtx_run2 is called as f_vtx)

    double f_norm_run3 = 0.677; //(52.8 / 77.904) (https://indico.cern.ch/event/1649731/contributions/6943855/attachments/3222296/5743342/Gagliardi_APW_180226.pdf)

    string path1 = "../../output/kstar/LHC22o_pass7/697595/kstarqa/hInvMass";
    TString outputPath = path1 + "/spectra_compare";
    gSystem->mkdir(outputPath, kTRUE);

    //==============First lets calcuate signal loss for INEL=======================
    TFile *fSigLoss;
    string pathNum, pathDen;
    THnSparseF *hSpraseGenAll, *hSpraseGenVzSel;
    TH1D *hGenAll, *hGenVzSel, *hSignalLoss;
    double eventLoss{1};
    if (!compareUsingINELgt0SigLoss)
    {
        fSigLoss = new TFile("../../mc/LHC24f3c/697699.root", "READ"); // 2024 MC INEL (TOF_overrideFT0)

        pathNum = "kstarqa_AllEventsVz"; // vZ only
        pathDen = "kstarqa_Vz_sel8";     // vZ + Sel8

        // //No effect of RCT is there.
        // pathNum = "kstarqa_AllEventsVz_NoRCT"; // vZ only
        // pathDen = "kstarqa_Vz_sel8_NoRCT";     // vZ + Sel8

        // // Gives same result as above
        // pathNum = "kstarqa_AllEvents_NoRCT"; // All events (without Vz)
        // pathDen = "kstarqa_NoVz_NoRCT";      // Sel8 (without Vz)

        hSpraseGenAll = (THnSparseF *)fSigLoss->Get(Form("%s/hInvMass/hk892GenpT", pathNum.c_str()));
        hSpraseGenVzSel = (THnSparseF *)fSigLoss->Get(Form("%s/hInvMass/hk892GenpT", pathDen.c_str()));
        if (hSpraseGenAll == nullptr || hSpraseGenVzSel == nullptr)
        {
            cout << "Error reading histograms" << endl;
            return;
        }
        hGenAll = hSpraseGenAll->Projection(0, "E");
        hGenVzSel = hSpraseGenVzSel->Projection(0, "E");

        hSignalLoss = new TH1D("hSignalLoss", "Signal Loss vs pT", sizeof(pT_bins) / sizeof(pT_bins[0]) - 1, pT_bins);
        for (int i = 1; i <= hSignalLoss->GetNbinsX(); i++)
        {
            double genAll = hGenAll->GetBinContent(i);
            double genVzSel = hGenVzSel->GetBinContent(i);
            double loss = genAll / genVzSel; // Signal loss
            hSignalLoss->SetBinContent(i, loss);
            double efficiencyerr = sqrt(abs(((genVzSel + 1) / (genAll + 2)) * ((genVzSel + 2) / (genAll + 3) - (genVzSel + 1) / (genAll + 2))));
            hSignalLoss->SetBinError(i, efficiencyerr);
        }

        TH1F *hDen = (TH1F *)fSigLoss->Get(Form("%s/hInvMass/MCcorrections/MultiplicityGen", pathDen.c_str())); // Without sel8
        TH1F *hNum = (TH1F *)fSigLoss->Get(Form("%s/hInvMass/MCcorrections/MultiplicityRec", pathDen.c_str())); // With Sel8
        double eventLossNum = hNum->Integral();
        double eventLossDen = hDen->Integral();
        eventLoss = eventLossNum / eventLossDen;
        cout << "Event Loss for INEL (Vz + Sel8): " << eventLoss << endl;

        TCanvas *cSignalLoss = new TCanvas("cSignalLoss", "Signal Loss vs pT", 720, 720);
        SetCanvasStyle(cSignalLoss, 0.16, 0.06, 0.01, 0.14);
        SetHistoQA(hSignalLoss);
        hSignalLoss->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hSignalLoss->GetYaxis()->SetTitle("Signal Loss");
        hSignalLoss->GetYaxis()->SetTitleOffset(1.6);
        hSignalLoss->SetMarkerSize(1.5);
        hSignalLoss->SetMaximum(1.38);
        hSignalLoss->SetMinimum(1.28);
        hSignalLoss->Draw("pe");
        cSignalLoss->SaveAs(outputPath + "/SignalLoss.png");
        //=================Signal loss end=======================
    }

    // TFile *fspectra1 = new TFile((path1 + "/corrected_spectra_INEL.root").c_str(), "read");
    TFile *fspectra1 = new TFile((path1 + "/corrected_spectra_0_120.root").c_str(), "read");

    if (fspectra1->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    TFile *fpub = new TFile("pp13TeV_INEL.root", "READ");
    if (fpub->IsZombie())
    {
        cout << "Error: file not found" << endl;
        return;
    }

    TH1F *hmult1, *hmultClone1;

    if (compareUsingINELgt0SigLoss)
    {
        hmult1 = (TH1F *)fspectra1->Get("mult_0-120/corrected_spectra_Integral_final");
    }
    else
    {
        hmult1 = (TH1F *)fspectra1->Get("mult_0-120/corrected_spectra_Integral");
    }

    if (hmult1 == nullptr)
    {
        cout << "Error: histograms not found in the corrected spectra root file" << endl;
        return;
    }
    hmult1->Scale(2.0); // for INEL it is K* + K*_bar and not their avearge (but confirm it once)
    // hmult1->Scale(f_norm_run2 * f_vtx_run2); // Adding same correction factors as run2 for comparison

    if (!compareUsingINELgt0SigLoss)
    {
        hmult1->Scale(f_norm_run3);
        hmult1->Multiply(hSignalLoss);
        // hmult1->Scale(eventLoss);
    }

    hmultClone1 = (TH1F *)hmult1->Clone("hmultClone0");

    TGraphErrors *gRun2_minBias = (TGraphErrors *)fpub->Get("Table 4/Graph1D_y1");
    if (hmult1 == nullptr)
    {
        cout << "Histogram 1 not found" << endl;
        return;
    }

    TH1F *h1 = (TH1F *)hmultClone1->Clone("h1");
    TH1F *h2 = (TH1F *)hmultClone1->Clone("h2");

    for (int i = 1; i <= h1->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr1 = (0.1 * h1->GetBinContent(i));
        h1->SetBinError(i, systemerr1);
    }

    TF1 *fitFcn1 = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
    fitFcn1->SetParameter(0, 5.0);
    fitFcn1->SetParameter(1, 0.07);
    fitFcn1->FixParameter(2, 0.895);
    fitFcn1->SetParameter(3, 0.3);
    fitFcn1->SetParNames("n", "dn/dy", "mass", "T");
    // fitFcn1->SetLineColor(color);
    fitFcn1->SetLineColor(kBlack);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(2);

    /*************meanpT*****************byresonance*******************package*************************/
    Double_t min = 0;
    Double_t max = 15;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 15;
    // Double_t maxfit=8.0;

    TH1 *hout = YieldMean(h1, h1, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    TGraphErrors *gratio1 = new TGraphErrors();

    if (!isSameBins)
    {
        for (int i = 0; i < gRun2_minBias->GetN(); i++) // took 4th for less error bar
        {
            double x_run2, yield_run2, x_error, y_error_run2;
            gRun2_minBias->GetPoint(i, x_run2, yield_run2);
            x_error = gRun2_minBias->GetErrorX(i);
            y_error_run2 = gRun2_minBias->GetErrorY(i);

            if (removeNormFactorsRun2)
            {
                yield_run2 = yield_run2 / (f_norm_run2 * f_vtx_run2);
                y_error_run2 = y_error_run2 / (f_norm_run2 * f_vtx_run2);
                gRun2_minBias->SetPoint(i, x_run2, yield_run2);
                gRun2_minBias->SetPointError(i, x_error, y_error_run2);
            }

            double thisanalysis1 = fitFcn1->Eval(x_run2);
            gratio1->SetPoint(i, x_run2, thisanalysis1 / yield_run2);
            cout << "Ratio 1 is " << thisanalysis1 / yield_run2 << endl;
            double error1 = sqrt(pow(thisanalysis1 * y_error_run2 / (yield_run2 * yield_run2), 2));
            gratio1->SetPointError(i, x_error, error1);
        }
    }
    else
    {
        if (gRun2_minBias->GetN() != hmult1->GetNbinsX())
        {
            cout << "Error: Number of points in Run2 minBias graph does not match number of bins in histogram." << endl;
            cout << "Number of points in graphs is " << gRun2_minBias->GetN() << endl;
            cout << "Number of bins in histogram is " << hmult1->GetNbinsX() << endl;
            return;
        }
        for (int i = 0; i < gRun2_minBias->GetN(); i++)
        {
            double x_run2, yield_run2, x_error, y_error_run2;
            gRun2_minBias->GetPoint(i, x_run2, yield_run2);
            x_error = gRun2_minBias->GetErrorX(i);
            y_error_run2 = gRun2_minBias->GetErrorY(i);

            if (removeNormFactorsRun2)
            {
                yield_run2 = yield_run2 / (f_norm_run2 * f_vtx_run2);
                y_error_run2 = y_error_run2 / (f_norm_run2 * f_vtx_run2);
                gRun2_minBias->SetPoint(i, x_run2, yield_run2);
                gRun2_minBias->SetPointError(i, x_error, y_error_run2);
            }

            double binvalue = hmultClone1->GetBinContent(i + 1);
            gratio1->SetPoint(i, x_run2, binvalue / yield_run2);
            double error1 = sqrt(pow(binvalue * y_error_run2 / (yield_run2 * yield_run2), 2));
            gratio1->SetPointError(i, x_error, error1);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    SetHistoStyle(h1, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    // SetHistoStyle(h21, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    h1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    h1->SetMaximum(h1->GetMaximum() * 5);
    h1->SetMinimum(h1->GetMinimum() * 0.1);
    h1->GetYaxis()->SetTitleOffset(1.30);
    h1->GetXaxis()->SetTitleOffset(1.02);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1);
    h1->GetXaxis()->SetRangeUser(0, 15);
    h1->SetLineColor(kBlue);
    h1->SetMarkerColor(kBlue);
    h1->Draw("pe");

    fitFcn1->SetLineColor(kBlue);
    fitFcn1->SetLineStyle(2);
    fitFcn1->Draw("same");

    gPad->SetLogy(1);
    gRun2_minBias->SetMarkerStyle(22);
    gRun2_minBias->SetMarkerSize(1);
    gRun2_minBias->SetMarkerColor(kBlack);
    gRun2_minBias->SetLineColor(kBlack);
    gRun2_minBias->SetLineWidth(2);
    gRun2_minBias->Draw("pe same");

    TLegend *leg = new TLegend(0.4, 0.65, 0.9, 0.91);
    SetLegendStyle(leg);
    leg->AddEntry((TObject *)0, "INEL", "");
    leg->AddEntry(h1, "LHC24_pass1_minBias", "p");
    leg->AddEntry(gRun2_minBias, "pp 13 TeV (Published)", "p");
    leg->SetTextSize(0.05);
    leg->Draw();

    c1->cd(2);
    TH1F *hdummy = (TH1F *)h1->Clone();
    for (int i = 0; i < hdummy->GetNbinsX(); i++)
    {
        hdummy->SetBinContent(i + 1, 0);
        hdummy->SetBinError(i + 1, 0);
    }

    SetGraphErrorStyle(gratio1);
    // SetGraphErrorStyle(gratio2);
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
    gratio1->GetXaxis()->SetRangeUser(0, 15);
    // gratio1->GetHistogram()->SetMaximum(gratio1->GetHistogram()->GetMaximum() * 1.5);
    // gratio1->GetHistogram()->SetMinimum(gratio1->GetHistogram()->GetMinimum() * 0.5);
    gratio1->SetMinimum(0.2);
    gratio1->SetMaximum(1.95);
    gratio1->Draw("ap");

    TLine *line = new TLine(0, 1, 15, 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(1);
    line->Draw();
    TString saveName = (removeNormFactorsRun2) ? "YieldRatioINEL_woNormFactorsRun2.png" : "YieldRatioINEL.png";
    c1->SaveAs(outputPath + "/" + saveName);
}