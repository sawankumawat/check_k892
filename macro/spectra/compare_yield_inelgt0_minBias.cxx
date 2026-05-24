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

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}

void compare_yield_inelgt0_minBias()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    bool isSameBins = true;
    // double inelNormRun2 = 0.892; event loss factor used in run 2
    // 663738 (OnlyTPC)
    // 664559 (TOFshift, TOFshiftMID)
    // 670168 (2023 data with two MC productions)
    // 668605 (2024 data: Base, MIDptDep2, MIDptDep2_small, MIDptDep2_verySmall, TOF3_withoutSquareCut)
    // 672297 (2024 data: MIDptDep2_TOF3, MIDptDep2_small_TOF3, MIDptDep2_0p3_TOF3)
    // 675391 (2024 data: MIDNew_TOF2, MIDNew_TOF3, SquarePID_TOF2, SquarePID_TOF3)

    string path1 = "../../output/kstar/LHC22o_pass7/672297/kstarqa_TOF3_withoutSquareCut/hInvMass"; // 2024 data
    string path2 = "../../output/kstar/LHC22o_pass7/672297/kstarqa_TOF3_withoutSquareCut/hInvMass"; // 2023 data
    TString outputPath = path2 + "/spectra_compare";
    gSystem->mkdir(outputPath, kTRUE);

    TFile *fspectra1 = new TFile((path1 + "/corrected_spectra.root").c_str(), "read");
    TFile *fspectra2 = new TFile((path2 + "/corrected_spectra.root").c_str(), "read");
    // TFile *fspectra3 = new TFile((path3 + "/corrected_spectra.root").c_str(), "read");

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

    TH1F *hmult1, *hmult2, *hmult3;
    TH1F *hmultClone1, *hmultClone2, *hmultClone3;
    TGraphErrors *gRun2_spectra[numofmultbins];
    TGraphErrors *gRun2_ratio[numofmultbins];
    TGraphErrors *gRun2_minBias[numofmultbins];

    hmult1 = (TH1F *)fspectra1->Get("mult_0-100/corrected_spectra_Integral_final");
    hmult2 = (TH1F *)fspectra2->Get("mult_0-100/corrected_spectra_Integral_final");
    // hmult1 = (TH1F *)fspectra1->Get("mult_0-100/corrected_spectra_Integral");
    // hmult2 = (TH1F *)fspectra2->Get("mult_0-100/corrected_spectra_Integral");
    // hmult3 = (TH1F *)fspectra3->Get("mult_0-100/corrected_spectra_Integral");
    hmultClone1 = (TH1F *)hmult1->Clone("hmultClone0");
    hmultClone2 = (TH1F *)hmult2->Clone("hmultClone0");
    // hmultClone3 = (TH1F *)hmult3->Clone("hmultClone0");

    if (hmult1 == nullptr)
    {
        cout << "Histogram 1 not found" << endl;
        return;
    }

    // TCanvas *crunMinBias = new TCanvas("", "", 720, 720);
    // TLegend *leg = new TLegend(0.4, 0.65, 0.9, 0.91);
    // leg->SetTextSize(0.035);
    // leg->SetBorderSize(0);
    // leg->SetFillStyle(0);
    // leg->SetNColumns(3);

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
        // TCanvas *crunMinBias = new TCanvas(Form("crunMinBias_%d", imult), Form("crunMinBias_%d", imult), 720, 720);
        // SetCanvasStyle(crunMinBias, 0.15, 0.03, 0.03, 0.15);
        // gPad->SetLogy();
        // SetGrapherrorStyle(gRun2_minBias[imult]);
        // gRun2_minBias[imult]->SetMarkerStyle(20 + imult);
        // gRun2_minBias[imult]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        // gRun2_minBias[imult]->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
        // if (imult == 0)
        // {
        //     gRun2_minBias[imult]->Draw("AP PLC PMC");
        // }
        // else
        // {
        //     gRun2_minBias[imult]->Draw("P PLC PMC same");
        // }
        // // crunMinBias->SaveAs(outputPath + Form("/run2_minBias_spectra_%d.png", imult));
        // leg->AddEntry(gRun2_minBias[imult], Form("%.0f-%.0f", mult_classes[imult], mult_classes[imult + 1]), "p");
    }
    // leg->Draw();
    // crunMinBias->SaveAs(outputPath + "/run2_minBias_spectra_%d.png");

    TH1F *h1 = (TH1F *)hmultClone1->Clone("h1");
    TH1F *h2 = (TH1F *)hmultClone1->Clone("h2");

    TH1F *h21 = (TH1F *)hmultClone2->Clone("h21");
    TH1F *h22 = (TH1F *)hmultClone2->Clone("h22");

    // TH1F *h31 = (TH1F *)hmultClone3->Clone("h31");
    // TH1F *h32 = (TH1F *)hmultClone3->Clone("h32");

    for (int i = 1; i <= h1->GetNbinsX(); i++) // putting small systematic error by hand
    {
        double systemerr1 = (0.1 * h1->GetBinContent(i));
        double systemerr2 = (0.1 * h21->GetBinContent(i));
        // double systemerr3 = (0.1 * h31->GetBinContent(i));
        h1->SetBinError(i, systemerr1);
        h21->SetBinError(i, systemerr2);
        // h31->SetBinError(i, systemerr3);
    }

    // int numColors = gStyle->GetNumberOfColors();
    // int paletteIndex = (imult - 1) * numColors / numofmultbins;
    // paletteIndex = std::min(paletteIndex, numColors - 1); // Ensure within bounds
    // int color = gStyle->GetColorPalette(paletteIndex);

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

    // TF1 *fitFcn3 = new TF1("fitfunc3", FuncLavy, 0.0, 15.0, 4);
    // fitFcn3->SetParameter(0, 5.0);
    // fitFcn3->SetParameter(1, 0.07);
    // fitFcn3->FixParameter(2, 0.895);
    // fitFcn3->SetParameter(3, 0.3);
    // fitFcn3->SetParNames("n", "dn/dy", "mass", "T");
    // fitFcn3->SetLineColor(kRed);
    // fitFcn3->SetLineStyle(2);
    // fitFcn3->SetLineWidth(2);

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

    TH1 *hout = YieldMean(h1, h1, fitFcn1, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    TH1 *hout2 = YieldMean(h21, h21, fitFcn2, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);
    // TH1 *hout3 = YieldMean(h31, h31, fitFcn3, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

    TGraphErrors *gratio1 = new TGraphErrors();
    TGraphErrors *gratio2 = new TGraphErrors();
    // TGraphErrors *gratio3 = new TGraphErrors();
    int minBiasFromWhichGraph = 7;

    if (!isSameBins)
    {
        for (int i = 0; i < gRun2_minBias[minBiasFromWhichGraph]->GetN(); i++) // took 4th for less error bar
        {
            double x_run2, yield_run2, x_error, y_error_run2;
            gRun2_minBias[minBiasFromWhichGraph]->GetPoint(i, x_run2, yield_run2);
            x_error = gRun2_minBias[minBiasFromWhichGraph]->GetErrorX(i);
            y_error_run2 = gRun2_minBias[minBiasFromWhichGraph]->GetErrorY(i);

            double thisanalysis1 = fitFcn1->Eval(x_run2);
            double thisanalysis2 = fitFcn2->Eval(x_run2);
            // double thisanalysis3 = fitFcn3->Eval(x_run2);
            gratio1->SetPoint(i, x_run2, thisanalysis1 / yield_run2);
            gratio2->SetPoint(i, x_run2, thisanalysis2 / yield_run2);
            cout << "Ratio 1 is " << thisanalysis1 / yield_run2 << endl;
            // cout << "Ratio 2 is " << thisanalysis2 / yield_run2 << endl;
            // gratio3->SetPoint(i, x_run2, thisanalysis3 / yield_run2);
            double error1 = sqrt(pow(thisanalysis1 * y_error_run2 / (yield_run2 * yield_run2), 2));
            double error2 = sqrt(pow(thisanalysis2 * y_error_run2 / (yield_run2 * yield_run2), 2));
            // double error3 = sqrt(pow(thisanalysis3 * y_error_run2 / (yield_run2 * yield_run2), 2));
            gratio1->SetPointError(i, x_error, error1);
            gratio2->SetPointError(i, x_error, error2);
            // gratio3->SetPointError(i, x_error, error3);
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
            double x_run2, yield_run2, x_error, y_error_run2;
            gRun2_minBias[minBiasFromWhichGraph]->GetPoint(i, x_run2, yield_run2);
            if (i == 4 || i == 5)
                yield_run2 = yield_run2 * 1.03; // to account for event loss in run 2 for 15-20 and 20-30 multiplicity classes
            if (i == 6 || i == 7)
                yield_run2 = yield_run2 * 0.94; // to account for event loss in run 2 for 5-10 and 10-15 multiplicity classes
            if (i == 12)
                yield_run2 = yield_run2 * 0.87; // to account for event loss in run 2 for 0-1 multiplicity class
            if (i == 9)
                yield_run2 = yield_run2 * 1.05; // to account for event loss in run 2 for 1-5 multiplicity class
            gRun2_minBias[minBiasFromWhichGraph]->SetPoint(i, x_run2, yield_run2);
            x_error = gRun2_minBias[minBiasFromWhichGraph]->GetErrorX(i);
            y_error_run2 = gRun2_minBias[minBiasFromWhichGraph]->GetErrorY(i);
            double binvalue = hmultClone1->GetBinContent(i + 1);
            double binvalue2 = hmultClone2->GetBinContent(i + 1);
            gratio1->SetPoint(i, x_run2, binvalue / yield_run2);
            gratio2->SetPoint(i, x_run2, binvalue2 / yield_run2);
            double error1 = sqrt(pow(binvalue * y_error_run2 / (yield_run2 * yield_run2), 2));
            double error2 = sqrt(pow(binvalue2 * y_error_run2 / (yield_run2 * yield_run2), 2));
            gratio1->SetPointError(i, x_error, error1);
            gratio2->SetPointError(i, x_error, error2);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    SetHistoStyle(h1, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    SetHistoStyle(h21, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    h1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    h1->SetMaximum(h1->GetMaximum() * 5);
    h1->SetMinimum(h1->GetMinimum() * 0.5);
    h1->GetYaxis()->SetTitleOffset(1.30);
    h1->GetXaxis()->SetTitleOffset(1.02);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1);
    h1->GetXaxis()->SetRangeUser(0, 10);
    h1->SetLineColor(kBlue);
    h1->SetMarkerColor(kBlue);
    h1->Draw("pe");
    h21->SetMarkerStyle(21);
    h21->SetMarkerSize(1);
    h21->SetMarkerColor(kRed);
    h21->SetLineColor(kRed);
    // h21->Draw("pe same");
    // h31->SetMarkerStyle(22);
    // h31->SetMarkerSize(1);
    // h31->SetMarkerColor(kGreen +2);
    // h31->SetLineColor(kGreen +2);
    // h31->Draw("pe same");
    fitFcn1->SetLineColor(kBlue);
    fitFcn1->SetLineStyle(2);
    fitFcn1->Draw("same");
    // fitFcn2->SetLineColor(kRed);
    // fitFcn2->SetLineStyle(2);
    // fitFcn2->Draw("same");
    // fitFcn3->SetLineColor(kGreen + 2);
    // fitFcn3->SetLineStyle(2);
    // fitFcn3->Draw("same");
    gPad->SetLogy(1);
    gRun2_minBias[minBiasFromWhichGraph]->SetMarkerStyle(22);
    gRun2_minBias[minBiasFromWhichGraph]->SetMarkerSize(1);
    gRun2_minBias[minBiasFromWhichGraph]->SetMarkerColor(kBlack);
    gRun2_minBias[minBiasFromWhichGraph]->SetLineColor(kBlack);
    gRun2_minBias[minBiasFromWhichGraph]->SetLineWidth(2);
    gRun2_minBias[minBiasFromWhichGraph]->Draw("pe same");

    TLegend *leg = new TLegend(0.4, 0.65, 0.9, 0.91);
    SetLegendStyle(leg);
    leg->SetHeader("INEL > 0");
    leg->AddEntry(h1, "2024 data", "p");
    // leg->AddEntry(h21, "2024 data", "p");
    // leg->AddEntry(h31, "2024 data", "p");
    leg->AddEntry(fitFcn1, "Levy-Tsallis", "l");
    leg->AddEntry(gRun2_minBias[minBiasFromWhichGraph], "pp 13 TeV (Published)", "p");
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
    // gratio1->GetHistogram()->SetMaximum(2.4);
    // gratio1->GetHistogram()->SetMinimum(0.3);
    gratio1->GetHistogram()->SetMaximum(1.45);
    gratio1->GetHistogram()->SetMinimum(0.65);
    // gratio1->GetHistogram()->SetMaximum(gratio1->GetHistogram()->GetMaximum() * 1.5);
    // gratio1->GetHistogram()->SetMinimum(gratio1->GetHistogram()->GetMinimum() * 0.5);
    // gratio1->SetMinimum(0.45);
    gratio1->Draw("ap");
    gratio2->SetMarkerStyle(21);
    gratio2->SetMarkerSize(1.0);
    gratio2->SetMarkerColor(kRed);
    gratio2->SetLineColor(kRed);
    // gratio2->Draw("p same");
    // gratio3->SetMarkerStyle(22);
    // gratio3->SetMarkerSize(1.0);
    // gratio3->SetMarkerColor(kGreen + 2);
    // gratio3->SetLineColor(kGreen + 2);
    // gratio3->Draw("p same");

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
    c1->SaveAs(outputPath + "/YieldMinBiasRatio.png");
}