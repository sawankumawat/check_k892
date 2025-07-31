#include <iostream>
#include <cmath>
#include "TArrow.h"
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
}

void plot_spectra()
{
    bool plotOnlyRaw = false;
    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/448490/kstarqa/hInvMass"; // 2022 data
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/449695/kstarqa/hInvMass"; // 2023 data
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/451993/kstarqa/hInvMass"; // 2024 data

    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459845/kstarqa/hInvMass"; // 2022 data
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459908/kstarqa/hInvMass"; // 2023 data
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/460233/kstarqa/hInvMass"; // 2024 data

    //IR study
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/463114/kstarqa/hInvMass"; // 1-2 MHz
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535069/kstarqa/hInvMass"; // 14 kHz
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535545/kstarqa/hInvMass"; // 70 kHz
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535645/kstarqa/hInvMass"; // 135 kHz
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535999/kstarqa/hInvMass"; // 330 kHz
    string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/536106/kstarqa/hInvMass"; // 650 kHz

    // TFile *fspectra = new TFile((path + "/corrected_spectra.root").c_str(), "read");
    TFile *fspectra = (plotOnlyRaw) ? new TFile((path + "/yield.root").c_str(), "read") : new TFile((path + "/corrected_spectra.root").c_str(), "read");

    if (fspectra->IsZombie())
    {
        cout << "File not found" << endl;
        return;
    }
    int markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 32, 47};
    float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;
    TH1F *hmult[numofmultbins + 1];
    TH1F *hmultClone[numofmultbins + 1];
    hmult[0] = (plotOnlyRaw) ? (TH1F *)fspectra->Get("mult_0-100/yield_integral") : (TH1F *)fspectra->Get("mult_0-100/corrected_spectra_Integral");
    // hmult[0] = (TH1F *)fspectra->Get("mult_0-100/corrected_spectra_Integral");
    hmultClone[0] = (TH1F *)hmult[0]->Clone("hmultClone0");

    if (hmult[0] == nullptr)
    {
        cout << "Histogram 1 not found" << endl;
        return;
    }

    for (int i = 1; i < numofmultbins + 1; i++)
    {
        // hmult[i] = (TH1F *)fspectra->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", mult_classes[i - 1], mult_classes[i]));
        hmult[i] = (plotOnlyRaw) ? (TH1F *)fspectra->Get(Form("mult_%.0f-%.0f/yield_integral", mult_classes[i - 1], mult_classes[i])) : (TH1F *)fspectra->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", mult_classes[i - 1], mult_classes[i]));

        if (hmult[i] == nullptr)
        {
            cout << "Histogram others not found" << endl;
            return;
        }
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
        // hmult[i]->Scale(pow(2, numofmultbins + 1 - i));
        hmult[i]->Scale(pow(2, numofmultbins - i));
        // hmult[i]->Scale(pow(2, multiplicationFactors[i - 1]));

        hmult[i]->GetXaxis()->SetTitleSize(0.045);
        hmult[i]->GetYaxis()->SetTitleSize(0.045);
        hmult[i]->GetYaxis()->SetTitleOffset(1.3);
        hmult[i]->SetMaximum(hmult[1]->GetMaximum() * 15);
        hmult[i]->SetMinimum(2e-7);
        hmult[i]->SetMarkerStyle(markers[i - 1]);
        hmult[i]->Draw("pe same PLC PMC");
    }
    hmult[0]->SetMarkerStyle(markers[numofmultbins]);
    hmult[0]->SetMarkerSize(1.2);
    hmult[0]->Draw("pe same PLC PMC");

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
                double systemerr = (0.1 * h2->GetBinContent(i));
                h2->SetBinError(i, systemerr);
            }
            /*************meanpT*****************byresonance*******************package*************************/
            Double_t min = 0.0;
            Double_t max = 15.0;
            Double_t loprecision = 0.01;
            Double_t hiprecision = 0.1;
            Option_t *opt = "RI0+";
            TString logfilename = "log.root";
            Double_t minfit = 0;
            Double_t maxfit = 15;
            // Double_t maxfit=8.0;

            TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
            fitFcn->SetParameter(0, 5.0);
            // fitFcn->SetParameter(1, 0.05);
            fitFcn->SetParameter(1, 50);
            fitFcn->FixParameter(2, 0.895);
            fitFcn->SetParameter(3, 0.35);
            fitFcn->SetParNames("n", "dn/dy", "mass", "T");

            int numColors = gStyle->GetNumberOfColors();
            int paletteIndex = (imult - 1) * numColors / numofmultbins;
            paletteIndex = std::min(paletteIndex, numColors - 1); // Ensure within bounds
            int color = gStyle->GetColorPalette(paletteIndex);
            fitFcn->SetLineColor(color);

            TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
            c->cd(1);
            fitFcn->SetLineColor(color);
            fitFcn->SetLineWidth(2);
            fitFcn->SetLineStyle(2);
            fitFcn->Draw("l same");
        }
    }
    c->cd(1);

    TLegend *leg = new TLegend(0.37, 0.75, 0.95, 0.97);
    SetLegendStyle(leg);
    leg->SetNColumns(4);
    leg->AddEntry(hmult[0], "0-100%", "lpe");
    for (int i = 1; i < numofmultbins + 1; i++)
    {
        // leg->AddEntry(hmult[i], Form("%.0f-%.0f%%(#times2^{%d})", mult_classes[i - 1], mult_classes[i], numofmultbins + 1 - i), "lpe");
        leg->AddEntry(hmult[i], Form("%.0f-%.0f%%(#times2^{%d})", mult_classes[i - 1], mult_classes[i], numofmultbins - i), "lpe");
        // leg->AddEntry(hmult[i], Form("%.0f-%.0f%%(#times2^{%d})", mult_classes[i - 1], mult_classes[i], multiplicationFactors[i - 1]), "lpe");
    }
    leg->SetTextSize(0.03);
    leg->Draw();

    c->cd(2);
    gPad->SetLogy(1);
    hratios[1]->SetMinimum(0.1);
    hratios[1]->SetMaximum(10);
    hratios[1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hratios[1]->GetYaxis()->SetTitle("Ratio to 0-100%");
    hratios[1]->GetYaxis()->SetNdivisions(505);
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
        hratios[i]->Draw("pe same PLC PMC");
    }
    TLine *line = new TLine(0, 1, 20, 1);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(3);
    line->Draw();

    TString outputfolder = path;
    (plotOnlyRaw) ? c->SaveAs(outputfolder + "/spectra_raw.png") : c->SaveAs(outputfolder + "/spectra.png");

    if (!plotOnlyRaw)
    {

        TCanvas *ctemp = new TCanvas("ctemp", "ctemp", 720, 720);
        SetCanvasStyle(ctemp, 0.15, 0.03, 0.03, 0.15);
        gPad->SetLogy();

        double meanpT[numofmultbins], yield[numofmultbins];
        double meanpT_err[numofmultbins], yield_err[numofmultbins];

        // Now lets calculate the mean pT and yield as a function of dN_charge/deta
        for (int imult = 1; imult < numofmultbins + 1; imult++)
        {
            hmultClone[imult]->Scale(0.5); // In run the average of K* and anti-K* is taken. so we have to scale it.
            TH1F *h1 = (TH1F *)hmultClone[imult]->Clone("h1");
            TH1F *h2 = (TH1F *)hmultClone[imult]->Clone("h2");
            // TH1F *h1 = (TH1F *)hmult[imult]->Clone("h1");
            // TH1F *h2 = (TH1F *)hmult[imult]->Clone("h2");

            for (int i = 1; i <= h1->GetNbinsX(); i++) // putting small systematic error by hand
            {
                double systemerr = (0.1 * h1->GetBinContent(i));
                h1->SetBinError(i, systemerr);
            }
            /*************meanpT*****************byresonance*******************package*************************/
            Double_t min = 0.0;
            Double_t max = 15.0;
            Double_t loprecision = 0.01;
            Double_t hiprecision = 0.1;
            Option_t *opt = "RI0+";
            TString logfilename = "log.root";
            Double_t minfit = 0.0;
            Double_t maxfit = 16.0;
            // Double_t maxfit=8.0;

            TF1 *fitFcn = new TF1("fitfunc", FuncLavy, 0.0, 15.0, 4);
            fitFcn->SetParameter(0, 5.0);
            fitFcn->SetParameter(1, 0.05);
            // fitFcn->SetParameter(1, 50);
            fitFcn->FixParameter(2, 0.895);
            fitFcn->SetParameter(3, 0.35);
            fitFcn->SetParNames("n", "dn/dy", "mass", "T");

            int numColors = gStyle->GetNumberOfColors();
            int paletteIndex = (imult - 1) * numColors / numofmultbins;
            paletteIndex = std::min(paletteIndex, numColors - 1); // Ensure within bounds
            int color = gStyle->GetColorPalette(paletteIndex);
            fitFcn->SetLineColor(color);

            TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
            meanpT[imult - 1] = hout->GetBinContent(5);
            meanpT_err[imult - 1] = hout->GetBinContent(6);
            yield[imult - 1] = hout->GetBinContent(1);
            yield_err[imult - 1] = hout->GetBinContent(2);
            for (Int_t ip = 0; ip < 9; ip++)
            {
                cout << hout->GetBinContent(ip + 1) << endl;
            }

            ctemp->cd(0);
            hmultClone[imult]->SetMarkerStyle(markers[imult]);
            hmultClone[imult]->SetMarkerSize(1.2);
            hmultClone[imult]->GetYaxis()->SetTitleOffset(1.3);
            hmultClone[imult]->SetMaximum(hmultClone[1]->GetMaximum() * 15);
            hmultClone[imult]->SetMinimum(2e-7);
            hmultClone[imult]->SetLineColor(color);
            hmultClone[imult]->SetMarkerColor(color);
            hmultClone[imult]->Draw("pe same");
            fitFcn->SetLineColor(color);
            fitFcn->SetLineWidth(2);
            fitFcn->SetLineStyle(2);
            fitFcn->Draw("l same");
        }
        ctemp->SaveAs(outputfolder + "/spectra_temp.png");

        // double dnch_detaRun2[] = {26.18, 20.16, 16.4, 13.14, 10.3, 8.24, 6.62, 4.77, 2.76};
        double dnch_detaRun3[] = {21.78, 18.48, 15.76, 13.195, 10.86, 9.09, 7.63, 5.87, 3.69};
        TGraphErrors *gMeanYieldRun3 = new TGraphErrors(numofmultbins, dnch_detaRun3, yield, nullptr, yield_err);
        TGraphErrors *gMeanpTRun3 = new TGraphErrors(numofmultbins, dnch_detaRun3, meanpT, nullptr, meanpT_err);

        TFile *fRun2 = new TFile("spectra/pp13TeV_INELgt0.root", "read");
        if (fRun2->IsZombie())
        {
            cout << "Run 2 file not found" << endl;
            return;
        }
        TGraphErrors *gMeanYieldRun2 = (TGraphErrors *)fRun2->Get("Table 41/Graph1D_y1");
        TGraphErrors *gMeanpTRun2 = (TGraphErrors *)fRun2->Get("Table 39/Graph1D_y1");
        if (gMeanYieldRun2 == nullptr || gMeanpTRun2 == nullptr)
        {
            cout << "Run 2 graph not found" << endl;
            return;
        }

        TFile *fLevyFit = new TFile(outputfolder + "/levy_fit.root", "recreate");

        TCanvas *cMeanYield = new TCanvas("cMeanYield", "cMeanYield", 720, 720);
        SetCanvasStyle(cMeanYield, 0.15, 0.03, 0.03, 0.15);
        SetGrapherrorStyle(gMeanYieldRun2);
        gMeanYieldRun2->SetMarkerStyle(20);
        gMeanYieldRun2->SetMarkerSize(1.2);
        gMeanYieldRun2->SetMarkerColor(kRed);
        gMeanYieldRun2->SetLineColor(kRed);
        gMeanYieldRun2->GetXaxis()->SetTitle("dN_{ch}/d#eta");
        gMeanYieldRun2->GetYaxis()->SetTitle("dN/dy");
        gMeanYieldRun2->Draw("AP");
        gMeanYieldRun2->Write("gMeanYieldRun2");
        SetGrapherrorStyle(gMeanYieldRun3);
        gMeanYieldRun3->SetMarkerStyle(21);
        gMeanYieldRun3->SetMarkerSize(1.2);
        gMeanYieldRun3->SetMarkerColor(kBlue);
        gMeanYieldRun3->SetLineColor(kBlue);
        gMeanYieldRun3->Draw("P same");
        gMeanYieldRun3->Write("gMeanYieldRun3");
        TLegend *legMeanYield = new TLegend(0.2, 0.80, 0.45, 0.90);
        legMeanYield->SetTextSize(0.04);
        legMeanYield->SetBorderSize(0);
        legMeanYield->SetFillStyle(0);
        legMeanYield->AddEntry(gMeanYieldRun2, "Run 2", "p");
        legMeanYield->AddEntry(gMeanYieldRun3, "Run 3", "p");
        legMeanYield->Draw();
        cMeanYield->SaveAs(outputfolder + "/mean_yield_run2.png");

        TCanvas *cMeanpT = new TCanvas("cMeanpT", "cMeanpT", 720, 720);
        SetCanvasStyle(cMeanpT, 0.15, 0.03, 0.03, 0.15);
        SetGrapherrorStyle(gMeanpTRun2);
        gMeanpTRun2->SetMarkerStyle(20);
        gMeanpTRun2->SetMarkerSize(1.2);
        gMeanpTRun2->SetMarkerColor(kRed);
        gMeanpTRun2->SetLineColor(kRed);
        gMeanpTRun2->GetXaxis()->SetTitle("dN_{ch}/d#eta");
        gMeanpTRun2->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
        gMeanpTRun2->Draw("AP");
        gMeanpTRun2->Write("gMeanpTRun2");
        SetGrapherrorStyle(gMeanpTRun3);
        gMeanpTRun3->SetMarkerStyle(21);
        gMeanpTRun3->SetMarkerSize(1.2);
        gMeanpTRun3->SetMarkerColor(kBlue);
        gMeanpTRun3->SetLineColor(kBlue);
        gMeanpTRun3->Draw("P same");
        legMeanYield->Draw();
        gMeanpTRun3->Write("gMeanpTRun3");
        cMeanpT->SaveAs(outputfolder + "/mean_pT_run2.png");
    }
}