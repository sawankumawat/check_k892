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

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void compare_yield_inel()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    bool isSameBins = false;

    ////Normalization factors used in run 2 INEL study:
    // double f_norm_run2 = 0.7448;
    // double f_vtx_run2 = 0.931264; // (In Phi AN, f_vtx_run2 is called as f_vtx)

    double f_norm_run3 = 0.677; //(52.8 / 77.904)

    string path1 = "../../output/kstar/LHC22o_pass7/788987/kstarqa/hInvMass/ROTATED";
    string pathSys = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/756343/kstarqa/hInvMass/SystematicsPlots";
    TString outputPath = path1 + "/spectra_compare";
    gSystem->mkdir(outputPath, kTRUE);

    TFile *fSys = new TFile((pathSys + "/SysUncert.root").c_str(), "READ");
    TH1D *hSys = (TH1D *)fSys->Get("hTotalSysSmoothed_0_120");

    //==============First lets calcuate signal loss for INEL=======================
    TFile *fSigLoss;
    // fSigLoss = new TFile("../../mc/LHC24f3c/755334.root", "READ"); // 2024 MC INEL (TOF_overrideFT0)
    // fSigLoss = new TFile("../../mc/LHC24f3c/790177.root", "READ"); // 2024 MC INEL lf-mc-centrality
    fSigLoss = new TFile("../../mc/LHC24f3c/801014.root", "READ"); // 2024 MC INEL lf-mc-centrality

    string pathNum, pathDen;
    THnSparseF *hSpraseGenAll, *hSpraseGenVzSel;
    TH1D *hGenAll, *hGenVzSel, *hSignalLoss;
    double eventLoss;

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

    hmult1 = (TH1F *)fspectra1->Get("mult_0-120/corrected_spectra_Integral_final");

    if (hmult1 == nullptr)
    {
        cout << "Error: histograms not found in the corrected spectra root file" << endl;
        return;
    }

    TH1F *hDen = (TH1F *)fSigLoss->Get("kstarqa/hInvMass/MCcorrections/MultiplicityGen");
    TH1F *hNum = (TH1F *)fSigLoss->Get("kstarqa/hInvMass/MCcorrections/MultiplicityRec");
    TH1F *hNumExcludingTVX = (TH1F *)fSigLoss->Get("kstarqa/hInvMass/MCcorrections/MultiplicityRecINEL");

    double eventLossNum = hNum->Integral();
    double eventLossDen = hDen->Integral();
    double eventLossNumExcludingTVX = hNumExcludingTVX->Integral();
    eventLoss = eventLossNum / eventLossDen;
    double eventLossExcludingTVX = eventLossNumExcludingTVX / eventLossDen;
    cout << "Event Loss for INEL: " << eventLoss << endl;
    cout << "Event Loss for INEL excluding TVX: " << eventLossExcludingTVX << endl;

    hmult1->Scale(2.0);         // for INEL it is K* + K*_bar and not their average
    hmult1->Scale(f_norm_run3); // visible cross section ratio
    hmult1->Scale(1.0 / eventLoss); // removing event loss

    // hmult1->Scale(eventLossExcludingTVX);          // Event loss without TriggerTVX

    // Data-based corrections
    // Selected events: 12478172244
    // TVX events: 9235908442
    double ratioSel_TVX = 9235908442.0 / 12478172244.0;
    hmult1->Scale(ratioSel_TVX); // selected events ==> TVX events
    // hmult1->Scale(1 / 0.95);     // Vertex-z correction

    // ==============================================================================
    //   Save corrected INEL spectrum with statistical and systematic uncertainties
    // ==============================================================================

    TH1F *hINEL_stat = (TH1F *)hmult1->Clone("hINEL_stat");
    hINEL_stat->SetTitle("INEL spectrum with statistical uncertainties");
    hINEL_stat->SetDirectory(0);

    // Systematic uncertainty histogram
    TH1F *hINEL_sys = (TH1F *)hmult1->Clone("hINEL_sys");
    hINEL_sys->SetTitle("INEL spectrum with systematic uncertainties");
    hINEL_sys->SetDirectory(0);

    if (hSys->GetNbinsX() != hINEL_sys->GetNbinsX())
    {
        cout << "Error: Systematic histogram and spectrum have different bin counts!" << endl;
        return;
    }

    for (int i = 1; i <= hINEL_sys->GetNbinsX(); i++)
    {
        double sysError = hSys->GetBinContent(i) * hINEL_sys->GetBinContent(i);
        hINEL_sys->SetBinError(i, sysError);
    }

    // Recreate ROOT file
    TFile *fINEL = new TFile((path1 + "/INEL_spectra.root").c_str(), "RECREATE");

    hINEL_stat->Write();
    hINEL_sys->Write();

    fINEL->Close();

    cout << "Saved INEL spectra with statistical and systematic uncertainties: " << path1 << "/INEL_spectra.root" << endl;

    hmultClone1 = (TH1F *)hmult1->Clone("hmultClone0");

    TGraphErrors *gRun2_minBias = (TGraphErrors *)fpub->Get("Table 4/Graph1D_y1");
    if (gRun2_minBias == nullptr)
    {
        cout << "Run2 minBias graph not found" << endl;
        return;
    }

    TH1F *h1 = (TH1F *)hmultClone1->Clone("h1");
    TH1F *h2 = (TH1F *)hmultClone1->Clone("h2");

    for (int i = 1; i <= h2->GetNbinsX(); i++) // putting small systematic error by hand
    {
        // double systemerr1 = (0.1 * h2->GetBinContent(i));
        double systemerr1 = hSys->GetBinContent(i) * h2->GetBinContent(i); // Using the smoothed systematic error from the SysUncert.root file
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
    Double_t max = 15;
    Double_t loprecision = 0.01;
    Double_t hiprecision = 0.1;
    Option_t *opt = "RI0+";
    TString logfilename = "log.root";
    Double_t minfit = 0;
    Double_t maxfit = 15;

    TH1 *hout = YieldMean(h1, h2, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    TGraphErrors *gratio1 = new TGraphErrors();
    TGraphErrors *gratio1_sys = new TGraphErrors();

    if (!isSameBins)
    {
        for (int i = 0; i < gRun2_minBias->GetN(); i++)
        {
            double x_run2, yield_run2, x_error, y_error_run2;
            gRun2_minBias->GetPoint(i, x_run2, yield_run2);
            x_error = gRun2_minBias->GetErrorX(i);
            y_error_run2 = gRun2_minBias->GetErrorY(i);

            double thisanalysis1 = fitFcn1->Eval(x_run2);

            int sysBin = h2->GetXaxis()->FindBin(x_run2);
            double sysError_fit = hSys->GetBinContent(sysBin) * thisanalysis1;

            double error1_stat = sqrt(pow(thisanalysis1 * y_error_run2 / (yield_run2 * yield_run2), 2));

            double error1_sys = sqrt(pow(sysError_fit / yield_run2, 2) + pow(thisanalysis1 * y_error_run2 / (yield_run2 * yield_run2), 2));

            gratio1->SetPoint(i, x_run2, thisanalysis1 / yield_run2);
            gratio1->SetPointError(i, x_error, error1_stat);

            gratio1_sys->SetPoint(i, x_run2, thisanalysis1 / yield_run2);

            double xBand = (x_error > 0) ? x_error : 0.10;
            gratio1_sys->SetPointError(i, xBand, error1_sys);

            // cout << "Ratio 1 is " << thisanalysis1 / yield_run2 << endl;
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

            double yieldError_run3 = h1->GetBinError(i + 1);
            double sysError_run3 = hSys->GetBinError(i + 1) * h1->GetBinContent(i + 1);
            double yield_Run3 = hmultClone1->GetBinContent(i + 1);

            gratio1->SetPoint(i, x_run2, yield_Run3 / yield_run2);
            gratio1_sys->SetPoint(i, x_run2, yield_Run3 / yield_run2);

            double error1 = sqrt(pow(yieldError_run3 / yield_run2, 2) + pow(yield_Run3 * y_error_run2 / (yield_run2 * yield_run2), 2));
            double error1_sys = sqrt(pow(sysError_run3 / yield_run2, 2) + pow(yield_Run3 * y_error_run2 / (yield_run2 * yield_run2), 2));

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
    hmultClone1->SetMaximum(hmultClone1->GetMaximum() * 5);
    hmultClone1->SetMinimum(hmultClone1->GetMinimum() * 0.15);
    hmultClone1->GetYaxis()->SetTitleOffset(1.30);
    hmultClone1->GetXaxis()->SetTitleOffset(1.02);
    hmultClone1->SetMarkerStyle(20);
    hmultClone1->SetMarkerSize(1);
    hmultClone1->GetXaxis()->SetRangeUser(0, 15);
    hmultClone1->SetLineColor(kBlue);
    hmultClone1->SetMarkerColor(kBlue);
    hmultClone1->Draw("pe");
    h2->SetMarkerColor(kBlue);
    h2->SetLineColor(kBlue);
    h2->SetFillStyle(0);
    h2->SetLineWidth(1);
    h2->Draw("e2 same");

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

    TLegend *leg = new TLegend(0.233983, 0.066092, 0.614206, 0.265189);
    SetLegendStyle(leg);
    leg->SetTextSize(0.045);
    leg->AddEntry(hmultClone1, "pp 13.6 TeV", "p");
    leg->AddEntry(gRun2_minBias, "pp 13 TeV (Eur. Phys. J. C 81:256)", "p");
    leg->AddEntry(fitFcn1, "Levy-Tsallis", "l");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.SetTextFont(42);
    lat.DrawLatex(0.65, 0.90, "ALICE");
    lat.DrawLatex(0.65, 0.83, "|y| < 0.5");
    lat.DrawLatex(0.65, 0.76, "INEL");
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
    gratio1->GetXaxis()->SetRangeUser(0, 15);
    gratio1->SetMinimum(0.2);
    gratio1->SetMaximum(1.95);
    gratio1->Draw("ap");

    gratio1_sys->SetMarkerStyle(20);
    gratio1_sys->SetMarkerSize(1.0);
    gratio1_sys->SetMarkerColor(kBlue);
    gratio1_sys->SetLineColor(kBlue);
    gratio1_sys->SetFillColor(kBlue);
    gratio1_sys->SetLineWidth(2);
    gratio1_sys->SetFillStyle(0);
    gratio1_sys->Draw("e2 same");

    TLine *line = new TLine(0, 1.01, 15, 1.01);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(1);
    line->Draw();

    // Draw a grey band with 20% uncertainty around the ratio of 1
    TBox *box = new TBox(0, 0.9, 15, 1.1);
    box->SetFillColor(kGray + 2);
    box->SetFillStyle(3003);
    box->Draw("same");

    TString saveName = "YieldRatioINEL_newMC.png";
    c1->SaveAs(outputPath + "/" + saveName);
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