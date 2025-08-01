#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);
void plot_spectra(vector<TString> paths, bool isCorrectedYield = false, vector<string> legendnames = {}, bool isSinglePanel = false);

std::vector<int> markerStyles = {20, 21, 22, 23, 29, 33, 34, 47, 48, 49};
std::vector<int> vibrantColors = {
    kBlack,      // #656364
    kBlue,       // #578dff
    kGreen + 3,  // #86c8dd
    kMagenta,    // #adad7d
    kRed,        // #c849a9
    kOrange + 3, // #c91116
    kSpring + 3, // #f5e002
    kCyan + 3    // #1845fb
};

void compare_multiple_spectra()
{
    bool isCorrectedYield = false; // set false for efficency and corrected yield plots.
    bool isSinglePanel = false;    // Set to true if you want to plot only the first panel
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    std::vector<string> runNumbers = {"535069", "535545", "535645", "535999", "536106", "463114"}; // Example run numbers
    std::vector<string> IR = {"14", "70", "135", "330", "650", "1500"};                            // interaction rates

    std::vector<TString> paths;
    for (const auto &run : runNumbers)
    {
        paths.emplace_back(Form("/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/%s/kstarqa/hInvMass", run.c_str()));
    }

    plot_spectra(paths, isCorrectedYield, IR, isSinglePanel);
}

void plot_spectra(vector<TString> paths, bool isCorrectedYield = false, vector<string> legendnames = {}, bool isSinglePanel = false)
{
    // float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    float mult_classes[] = {0};
    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;
    int totalfiles = paths.size();
    vector<TFile *> fspectra(totalfiles);
    vector<TH1F *> hmult[numofmultbins + 1];
    vector<TH1F *> hratio[numofmultbins + 1]; // Fixed: should be numofmultbins + 1
    for (int ifiles = 0; ifiles < totalfiles; ifiles++)
    {
        fspectra[ifiles] = (isCorrectedYield) ? new TFile((paths[ifiles] + "/corrected_spectra.root"), "read") : new TFile((paths[ifiles] + "/yield.root"), "read");
        if (fspectra[ifiles]->IsZombie())
        {
            cout << "File not found: " << paths[ifiles] << endl;
            return;
        }

        for (int imult = 0; imult < numofmultbins + 1; imult++)
        {
            double multlow = (imult == 0) ? 0 : mult_classes[imult - 1];
            double multhigh = (imult == 0) ? 100 : mult_classes[imult];
            hmult[imult].push_back((isCorrectedYield) ? (TH1F *)fspectra[ifiles]->Get(Form("mult_%.0f-%.0f/corrected_spectra_Integral", multlow, multhigh)) : (TH1F *)fspectra[ifiles]->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh)));
            // hmult[imult].push_back((isCorrectedYield) ? (TH1F *)fspectra[ifiles]->Get(Form("mult_%.0f-%.0f/heff", multlow, multhigh)) : (TH1F *)fspectra[ifiles]->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh)));
            if (hmult[imult][ifiles] == nullptr)
            {
                cout << "Histogram not found for mult bin " << imult << " in file " << paths[ifiles] << endl;
                return;
            }

            // Create ratio histogram (ratio to first file)
            hratio[imult].push_back((TH1F *)hmult[imult][ifiles]->Clone(Form("hratio_mult_%.0f-%.0f_file%d", multlow, multhigh, ifiles)));
            if (ifiles > 0)
            { // Only create meaningful ratios for files other than the first one
                hratio[imult][ifiles]->Divide(hmult[imult][0]);
            }
        }
    }

    for (int imult = 0; imult < numofmultbins + 1; imult++)
    {
        double multlow = (imult == 0) ? 0 : mult_classes[imult - 1];
        double multhigh = (imult == 0) ? 100 : mult_classes[imult];

        TCanvas *c1 = new TCanvas("", "", 720, 720);
        double pad1Size, pad2Size;
        if (!isSinglePanel)
        {
            SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
            canvas_style(c1, pad1Size, pad2Size);
            c1->cd(1);
        }
        else
        {
            SetCanvasStyle(c1, 0.15, 0.05, 0.03, 0.15);
            pad1Size = 4.0 / 4.5; // Full canvas size for single panel
            pad2Size = 4.0 / 4.5; // No second pad
        }
        TLegend *leg = new TLegend(0.6, 0.7, 0.93, 0.96);
        SetLegendStyle(leg);
        leg->SetNColumns(2);
        (isSinglePanel) ? leg->SetTextSize(0.035) : leg->SetTextSize(0.05);
        // leg->SetHeader(Form("Multiplicity: %.0f-%.0f%%", multlow, multhigh));
        leg->SetHeader("Interaction Rate (kHz)");

        gPad->SetLogy(1);
        for (int ifiles = 0; ifiles < totalfiles; ifiles++)
        {
            SetHistoQA(hmult[imult][ifiles]);
            hmult[imult][ifiles]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
            hmult[imult][ifiles]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
            hmult[imult][ifiles]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
            hmult[imult][ifiles]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
            hmult[imult][ifiles]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
            // hmult[imult][ifiles]->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
            hmult[imult][ifiles]->GetYaxis()->SetTitle("Acceptance x Efficiency");
            // hmult[imult][ifiles]->SetMinimum(hmult[imult][ifiles]->GetMinimum() * 1);
            if (!isSinglePanel)
            {
                hmult[imult][ifiles]->SetMaximum(hmult[imult][ifiles]->GetMaximum() * 10);
                hmult[imult][ifiles]->GetYaxis()->SetTitleOffset(1.30);
                hmult[imult][ifiles]->GetXaxis()->SetTitleOffset(1.18);
            }
            else
            {
                hmult[imult][ifiles]->SetMaximum(hmult[imult][ifiles]->GetMaximum() * 1.6);
                hmult[imult][ifiles]->GetYaxis()->SetTitleOffset(1.59);
                hmult[imult][ifiles]->GetXaxis()->SetTitleOffset(1.4);
            }
            hmult[imult][ifiles]->GetXaxis()->SetRangeUser(0, 10);
            hmult[imult][ifiles]->SetMarkerStyle(markerStyles[ifiles]);
            hmult[imult][ifiles]->SetMarkerColor(vibrantColors[ifiles]);
            hmult[imult][ifiles]->SetMarkerSize(1.5);
            hmult[imult][ifiles]->SetLineColor(vibrantColors[ifiles]);
            hmult[imult][ifiles]->SetMinimum(2e-6);
            hmult[imult][ifiles]->Draw("pe same");

            leg->AddEntry(hmult[imult][ifiles], Form("%s", legendnames[ifiles].c_str()), "p");
        }

        leg->Draw();
        TLatex *t2 = new TLatex();
        t2->SetNDC();
        (isSinglePanel) ? t2->SetTextSize(0.03) : t2->SetTextSize(0.045);
        t2->SetTextFont(42);
        t2->DrawLatex(0.22, 0.92, Form("%.0f-%.0f%%", multlow, multhigh));
        t2->DrawLatex(0.22, 0.86, Form("pp 13.6 TeV, 2023 pass4"));

        if (!isSinglePanel)
        {
            c1->cd(2);
            for (int ifiles = 1; ifiles < totalfiles; ifiles++) // Start from 1 to skip the reference file
            {
                SetHistoQA(hratio[imult][ifiles]);
                hratio[imult][ifiles]->GetYaxis()->SetTitleSize(0.03 / pad2Size);
                hratio[imult][ifiles]->GetXaxis()->SetTitleSize(0.04 / pad2Size);
                hratio[imult][ifiles]->GetYaxis()->SetLabelSize(0.04 / pad2Size);
                hratio[imult][ifiles]->GetXaxis()->SetLabelSize(0.04 / pad2Size);
                hratio[imult][ifiles]->SetMarkerStyle(markerStyles[ifiles]);
                hratio[imult][ifiles]->SetMarkerSize(1.5);
                hratio[imult][ifiles]->SetMarkerColor(vibrantColors[ifiles]);
                hratio[imult][ifiles]->SetLineColor(vibrantColors[ifiles]);
                hratio[imult][ifiles]->GetYaxis()->SetTitle(Form("Ratio to %s kHz", legendnames[0].c_str()));
                hratio[imult][ifiles]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
                hratio[imult][ifiles]->GetXaxis()->CenterTitle(1);
                hratio[imult][ifiles]->GetYaxis()->SetTitleOffset(0.75);
                hratio[imult][ifiles]->GetXaxis()->SetTitleOffset(1.1);
                hratio[imult][ifiles]->GetYaxis()->SetNdivisions(505);
                hratio[imult][ifiles]->GetXaxis()->SetRangeUser(0, 10);
                hratio[imult][ifiles]->SetMinimum(-0.09);
                hratio[imult][ifiles]->SetMaximum(2.16);
                hratio[imult][ifiles]->Draw("pe same");
            }
            TLine *line = new TLine(0, 1, 10, 1);
            line->SetLineStyle(2);
            line->SetLineColor(kBlack);
            line->SetLineWidth(2);
            line->Draw();
        }

        // c1->SaveAs(Form("%s/compareIR_Corrected_%.0f-%.0f.png", paths[0].Data(), multlow, multhigh));
    }
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