#include <iostream>
#include "../src/style.h"
#include "../src/fitfunc.h"
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
}

void compare_yield_inelgt0()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    string path1 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459845/kstarqa/hInvMass"; // 2022 data
    string path2 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459908/kstarqa/hInvMass"; // 2023 data
    string path3 = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/460233/kstarqa/hInvMass"; // 2024 data
    TString outputPath = path3;

    TFile *fspectra1 = new TFile((path1 + "/corrected_spectra.root").c_str(), "read");
    TFile *fspectra2 = new TFile((path2 + "/corrected_spectra.root").c_str(), "read");
    TFile *fspectra3 = new TFile((path3 + "/corrected_spectra.root").c_str(), "read");

    if (fspectra1->IsZombie() || fspectra2->IsZombie() || fspectra3->IsZombie())
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

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;

    TH1F *hmult1[numofmultbins + 1];
    TH1F *hmult2[numofmultbins + 1];
    TH1F *hmult3[numofmultbins + 1];
    TH1F *hmultClone1[numofmultbins + 1];
    TH1F *hmultClone2[numofmultbins + 1];
    TH1F *hmultClone3[numofmultbins + 1];
    TGraphErrors *g_vom[numofmultbins];

    hmult1[0] = (TH1F *)fspectra1->Get("mult_0-100/corrected_spectra_Integral");
    hmult2[0] = (TH1F *)fspectra2->Get("mult_0-100/corrected_spectra_Integral");
    hmult3[0] = (TH1F *)fspectra3->Get("mult_0-100/corrected_spectra_Integral");
    hmultClone1[0] = (TH1F *)hmult1[0]->Clone("hmultClone0");
    hmultClone2[0] = (TH1F *)hmult2[0]->Clone("hmultClone0");
    hmultClone3[0] = (TH1F *)hmult3[0]->Clone("hmultClone0");

    if (hmult1[0] == nullptr)
    {
        cout << "Histogram 1 not found" << endl;
        return;
    }

    for (int imult = 1; imult < numofmultbins + 1; imult++)
    // for (int imult = 1; imult < 2; imult++)
    {
        hmult1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", mult_classes[imult - 1], mult_classes[imult]));
        hmult2[imult] = (TH1F *)fspectra2->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", mult_classes[imult - 1], mult_classes[imult]));
        hmult3[imult] = (TH1F *)fspectra3->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", mult_classes[imult - 1], mult_classes[imult]));
        hmultClone1[imult] = (TH1F *)hmult1[imult]->Clone(Form("hmultClone%d", imult));
        hmultClone2[imult] = (TH1F *)hmult2[imult]->Clone(Form("hmultClone%d", imult));
        hmultClone3[imult] = (TH1F *)hmult3[imult]->Clone(Form("hmultClone%d", imult));
        if (hmult1[imult] == nullptr)
        {
            cout << "Histogram others not found" << endl;
            return;
        }
        g_vom[imult - 1] = (TGraphErrors *)fpub->Get(Form("Table %d/Graph1D_y1", imult)); // vom1
        if (g_vom[imult - 1] == nullptr)
        {
            cout << "Run2 yield graph not found for mult bin " << imult << endl;
            return;
        }

        hmultClone1[imult]->Scale(0.5); // In run the average of K* and anti-K* is taken. so we have to scale it.
        hmultClone2[imult]->Scale(0.5);
        hmultClone3[imult]->Scale(0.5);

        TH1F *h1 = (TH1F *)hmultClone1[imult]->Clone("h1");
        TH1F *h2 = (TH1F *)hmultClone1[imult]->Clone("h2");

        TH1F *h21 = (TH1F *)hmultClone2[imult]->Clone("h21");
        TH1F *h22 = (TH1F *)hmultClone2[imult]->Clone("h22");

        TH1F *h31 = (TH1F *)hmultClone3[imult]->Clone("h31");
        TH1F *h32 = (TH1F *)hmultClone3[imult]->Clone("h32");

        for (int i = 1; i <= h1->GetNbinsX(); i++) // putting small systematic error by hand
        {
            double systemerr1 = (0.1 * h1->GetBinContent(i));
            double systemerr2 = (0.1 * h21->GetBinContent(i));
            double systemerr3 = (0.1 * h31->GetBinContent(i));
            h1->SetBinError(i, systemerr1);
            h21->SetBinError(i, systemerr2);
            h31->SetBinError(i, systemerr3);
        }

        int numColors = gStyle->GetNumberOfColors();
        int paletteIndex = (imult - 1) * numColors / numofmultbins;
        paletteIndex = std::min(paletteIndex, numColors - 1); // Ensure within bounds
        int color = gStyle->GetColorPalette(paletteIndex);

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

        TF1 *fitFcn2 = new TF1("fitfunc2", FuncLavy, 0.0, 15.0, 4);
        fitFcn2->SetParameter(0, 5.0);
        fitFcn2->SetParameter(1, 0.07);
        fitFcn2->FixParameter(2, 0.895);
        fitFcn2->SetParameter(3, 0.3);
        fitFcn2->SetParNames("n", "dn/dy", "mass", "T");
        fitFcn2->SetLineColor(kRed);
        fitFcn2->SetLineStyle(2);
        fitFcn2->SetLineWidth(2);

        TF1 *fitFcn3 = new TF1("fitfunc3", FuncLavy, 0.0, 15.0, 4);
        fitFcn3->SetParameter(0, 5.0);
        fitFcn3->SetParameter(1, 0.07);
        fitFcn3->FixParameter(2, 0.895);
        fitFcn3->SetParameter(3, 0.3);
        fitFcn3->SetParNames("n", "dn/dy", "mass", "T");
        fitFcn3->SetLineColor(kRed);
        fitFcn3->SetLineStyle(2);
        fitFcn3->SetLineWidth(2);

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
        TH1 *hout2 = YieldMean(h21, h21, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
        TH1 *hout3 = YieldMean(h31, h31, fitFcn3, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

        TGraphErrors *gratio1 = new TGraphErrors();
        TGraphErrors *gratio2 = new TGraphErrors();
        TGraphErrors *gratio3 = new TGraphErrors();
        for (int i = 0; i < g_vom[imult - 1]->GetN(); i++)
        {
            double x_run2, yield_run2, x_error, y_error_run2;
            g_vom[imult - 1]->GetPoint(i, x_run2, yield_run2);
            x_error = g_vom[imult - 1]->GetErrorX(i);
            y_error_run2 = g_vom[imult - 1]->GetErrorY(i);

            double thisanalysis1 = fitFcn1->Eval(x_run2);
            double thisanalysis2 = fitFcn2->Eval(x_run2);
            double thisanalysis3 = fitFcn3->Eval(x_run2);
            // cout << "run 2 x is " << x_run2 << " y run2 is " << yield_run2 << "run 3 " << thisanalysis1 << endl;
            gratio1->SetPoint(i, x_run2, thisanalysis1 / yield_run2);
            gratio2->SetPoint(i, x_run2, thisanalysis2 / yield_run2);
            gratio3->SetPoint(i, x_run2, thisanalysis3 / yield_run2);
            double error1 = sqrt(pow(thisanalysis1 * y_error_run2 / (yield_run2 * yield_run2), 2));
            double error2 = sqrt(pow(thisanalysis2 * y_error_run2 / (yield_run2 * yield_run2), 2));
            double error3 = sqrt(pow(thisanalysis3 * y_error_run2 / (yield_run2 * yield_run2), 2));
            gratio1->SetPointError(i, x_error, error1);
            gratio2->SetPointError(i, x_error, error2);
            gratio3->SetPointError(i, x_error, error3);

            // double error = sqrt(pow(binerror / yield_run2, 2) + pow(bincontent * y_error_run2 / (yield_run2 * yield_run2), 2));
        }

        TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
        SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
        double pad1Size, pad2Size;
        canvas_style(c1, pad1Size, pad2Size);
        c1->cd(1);
        SetHistoStyle(h1, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
        SetHistoStyle(h21, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
        SetHistoStyle(h31, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
        h1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        h1->SetMaximum(h1->GetMaximum() * 10);
        h1->SetMinimum(h1->GetMinimum() * 1);
        h1->GetYaxis()->SetTitleOffset(1.30);
        h1->GetXaxis()->SetTitleOffset(1.02);
        h1->SetMarkerStyle(20);
        h1->SetMarkerSize(1);
        h1->GetXaxis()->SetRangeUser(0, 10);
        h1->Draw("pe");
        h21->SetMarkerStyle(21);
        h21->SetMarkerSize(1);
        h21->SetMarkerColor(kRed);
        h21->SetLineColor(kRed);
        h21->Draw("pe same");
        h31->SetMarkerStyle(22);
        h31->SetMarkerSize(1);
        h31->SetMarkerColor(kGreen +2);
        h31->SetLineColor(kGreen +2);
        h31->Draw("pe same");
        fitFcn1->SetLineColor(kBlack);
        fitFcn1->SetLineStyle(2);
        fitFcn1->Draw("same");
        fitFcn2->SetLineColor(kRed);
        fitFcn2->SetLineStyle(2);
        fitFcn2->Draw("same");
        fitFcn2->SetLineColor(kGreen + 2);
        fitFcn2->SetLineStyle(2);
        fitFcn2->Draw("same");
        gPad->SetLogy(1);
        g_vom[imult - 1]->SetMarkerStyle(22);
        g_vom[imult - 1]->SetMarkerSize(1);
        g_vom[imult - 1]->SetMarkerColor(kBlue);
        g_vom[imult - 1]->SetLineColor(kBlue);
        g_vom[imult - 1]->SetLineWidth(2);
        g_vom[imult - 1]->Draw("pe same");

        TLegend *leg = new TLegend(0.4, 0.67, 0.9, 0.91);
        SetLegendStyle(leg);
        leg->SetHeader(Form("Multiplicity: %.0f-%.0f%%", mult_classes[imult - 1], mult_classes[imult]));
        leg->AddEntry(h1, "2022 data", "p");
        leg->AddEntry(h21, "2023 data", "p");
        leg->AddEntry(h31, "2024 data", "p");
        // leg->AddEntry(fitFcn, "Levy-Tsallis fit (pp 13.6 TeV)", "l");
        leg->AddEntry(g_vom[imult - 1], "pp 13 TeV (Published)", "p");
        leg->SetTextSize(0.05);
        leg->Draw();

        c1->cd(2);
        TH1F *hdummy = (TH1F *)h1->Clone();
        for (int i = 0; i < hdummy->GetNbinsX(); i++)
        {
            hdummy->SetBinContent(i + 1, 0);
            hdummy->SetBinError(i + 1, 0);
        }

        SetGrapherrorStyle(gratio1);
        SetGrapherrorStyle(gratio2);
        SetGrapherrorStyle(gratio3);
        gratio1->GetYaxis()->SetTitleSize(0.035 / pad2Size);
        gratio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        gratio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        gratio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        gratio1->SetMarkerStyle(20);
        gratio1->SetMarkerSize(1.0);
        gratio1->SetMarkerColor(1);
        gratio1->SetLineColor(1);
        gratio1->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
        gratio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        gratio1->GetXaxis()->CenterTitle(1);
        gratio1->GetYaxis()->SetTitleOffset(0.6);
        gratio1->GetXaxis()->SetTitleOffset(1.1);
        gratio1->GetYaxis()->SetNdivisions(506);
        gratio1->GetXaxis()->SetRangeUser(0, 10);
        gratio1->GetHistogram()->SetMaximum(gratio3->GetHistogram()->GetMaximum() * 1.12);
        gratio1->GetHistogram()->SetMinimum(gratio1->GetHistogram()->GetMinimum() * 0.9);
        // gratio1->SetMinimum(0.45);
        gratio1->Draw("ap");
        gratio2->SetMarkerStyle(21);
        gratio2->SetMarkerSize(1.0);
        gratio2->SetMarkerColor(2);
        gratio2->SetLineColor(2);
        gratio2->Draw("p same");
        gratio3->SetMarkerStyle(22);
        gratio3->SetMarkerSize(1.0);
        gratio3->SetMarkerColor(kGreen + 2);
        gratio3->SetLineColor(kGreen + 2);
        gratio3->Draw("p same");

        TLine *line = new TLine(0, 1, 10, 1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(1);
        line->Draw();
        c1->SaveAs(outputPath + Form("/YieldRatioMult_%.0f-%.0f.png", mult_classes[imult - 1], mult_classes[imult]));
    }
}