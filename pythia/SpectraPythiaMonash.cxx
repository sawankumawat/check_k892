#include <iostream>
#include <vector>
#include "../macro/src/style.h"
#include "../macro/spectra/YieldMean.C"
using namespace std;

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);
////==In the pythia code, I have used the Nch > 1. I think this may be the reason for higher dNch/deta values.

void SpectraPythiaMonash()
{
    // gStyle->SetOptStat(1110);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    int colors[] = {kRed, kOrange, kYellow, kGreen, kCyan, kBlue, kViolet, kMagenta, kPink - 1, kGray};
    string modelNames[] = {"MonashLocal", "MonashWoCRLocal", "ShovingLocal", "MonashRescatteringLocal", "MonashLocal2"};
    int chooseModel = 4; // 0: Monash with CR, 1: Monash without CR, 2: Shoving, 3: Monash with hadronic rescattering

    TFile *fpp13p6TeV = new TFile(Form("NewSimPythia8p3/%s.root", modelNames[chooseModel].c_str()), "READ");
    if (fpp13p6TeV->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string particles[] = {"Phi", "Kstar", "K0s_pt", "Proton", "Pion", "Kaon", "ChKstar", "XiStar"};
    string saveNames[] = {"phi", "kstar", "kshort", "proton", "pion", "kaon", "chargedkstar", "xistar"};

    int numParticles = sizeof(particles) / sizeof(particles[0]);
    // TFile *fOutput = new TFile("SpectraPythiaMonash.root", "RECREATE");
    TFile *fOutput = new TFile(Form("Pythia_%s.root", modelNames[chooseModel].c_str()), "RECREATE");

    TH2D *hmultFT0vshmultFT0Eta0p5 = (TH2D *)fpp13p6TeV->Get("hmultFT0vshmultFT0Eta0p5");
    TH1D *hmultpp136 = (TH1D *)hmultFT0vshmultFT0Eta0p5->ProjectionX("hmultpp136");

    if (hmultpp136 == nullptr || hmultFT0vshmultFT0Eta0p5 == nullptr)
    {
        cout << "Error reading hmultFT0vshmultFT0Eta0p5 histogram" << endl;
        return;
    }

    TGraphErrors *gMeanPt = new TGraphErrors();
    TGraphErrors *gdNdy = new TGraphErrors();

    double totalEvents = hmultpp136->Integral(1, hmultpp136->GetNbinsX());

    // ALICE-style percentiles
    const int NumPercent = 19;
    // double percentiles[NumPercent] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    double percentiles[NumPercent] = {0.8, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0};

    vector<int> regionEdges;
    double cumulative = 0;
    int idx = 0;

    // Integrate from HIGH multiplicity → LOW
    for (int i = hmultpp136->GetNbinsX(); i >= 1; i--)
    {
        cumulative += hmultpp136->GetBinContent(i);

        while (idx < NumPercent && cumulative >= (percentiles[idx] / 100) * totalEvents)
        {
            regionEdges.push_back(i);
            idx++;
        }
    }

    // Safety (in case of rounding)
    cout << "Size of regionEdges: " << regionEdges.size() << endl;
    while ((int)regionEdges.size() < NumPercent)
        regionEdges.push_back(1);

    // Create region histograms
    TH1D *hRegions[NumPercent];
    for (int i = 0; i < NumPercent; i++)
    {
        hRegions[i] = (TH1D *)hmultpp136->Clone(Form("hRegion_%d", i));
        hRegions[i]->Reset();
        hRegions[i]->SetFillColor(colors[i]);
        hRegions[i]->SetLineColor(kBlack);
        hRegions[i]->SetFillStyle(3001);
    }

    // Fill regions (high → low)
    int prevBin = hmultpp136->GetNbinsX();

    for (int r = 0; r < NumPercent; r++)
    {
        int lowBin = regionEdges[r];

        for (int b = lowBin; b <= prevBin; b++)
        {
            hRegions[r]->SetBinContent(b, hmultpp136->GetBinContent(b));
        }

        // cout <<"Total event in region " << percentiles[NumPercent - r -1] * 100 << "% - " << hmultpp136->Integral(lowBin, prevBin) << endl;

        prevBin = lowBin - 1;
    }

    // -------------------------------
    // Plot multiplicity distribution
    // -------------------------------
    TCanvas *cMult = new TCanvas("cMult", "Multiplicity", 720, 720);
    SetCanvasStyle(cMult, 0.14, 0.01, 0.01, 0.14);
    SetHistoQA(hmultpp136);
    gPad->SetLogy();
    gPad->SetLogx();
    hmultpp136->SetMaximum(hmultpp136->GetMaximum() * 20);
    hmultpp136->GetXaxis()->SetTitle("N_{ch}");
    hmultpp136->GetYaxis()->SetTitle("Events");
    hmultpp136->Draw();

    // const char *labels[NumPercent] = {"0-1%", "1-5%", "5-10%", "10-15%", "15-20%", "20-30%", "30-40%", "40-50%", "50-70%", "70-100%"};
    const char *labels[NumPercent] = {"0-0.8%", "0.8-1.5%", "1.5-2%", "2-3%", "3-5%", "5-7%", "7-10%", "10-13%", "13-16%", "16-20%", "20-25%", "25-30%", "30-35%", "35-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-100%"};

    TLegend *leg = new TLegend(0.3, 0.87, 0.95, 0.98);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->SetNColumns(5);

    for (int i = 0; i < NumPercent; i++)
    {
        hRegions[i]->Draw("HIST SAME");
        leg->AddEntry(hRegions[i], labels[i], "f");

        const double eventsInClass = hRegions[i]->Integral(1, hRegions[i]->GetNbinsX());
        cout << labels[i] << " events = " << eventsInClass
             << " (" << 100.0 * eventsInClass / totalEvents << "%)" << endl;
    }
    leg->Draw();

    cout << "\n===== <dN/deta>|eta|<0.5 =====" << endl;

    int highBin = hmultpp136->GetNbinsX();

    TH1D *hChMultEta0p5[NumPercent];

    for (int i = 0; i < NumPercent; i++)
    {
        int lowBin = regionEdges[i];

        hmultFT0vshmultFT0Eta0p5->GetXaxis()->SetRange(lowBin, highBin);

        // TH1D *proj = (TH1D *)hmultFT0vshmultFT0Eta0p5->ProjectionY(Form("proj_%d", i));
        hChMultEta0p5[i] = (TH1D *)hmultFT0vshmultFT0Eta0p5->ProjectionY(Form("hChMultEta0p5_%d", i));

        cout << labels[i] << " : " << hChMultEta0p5[i]->GetMean() << endl;

        highBin = lowBin - 1;
    }

    for (int ipar = 0; ipar < numParticles; ipar++)
    {
        cout << endl;
        TH2D *hmultFT0vskstar_pt = (TH2D *)fpp13p6TeV->Get(Form("hmultFT0vs%s", particles[ipar].c_str()));
        if (hmultFT0vskstar_pt == nullptr)
        {
            cout << "Error reading histogram for particle " << particles[ipar] << endl;
            continue;
        }

        if (hmultFT0vskstar_pt->GetNbinsX() != hmultpp136->GetNbinsX())
        {
            cout << "Warning: multiplicity-axis binning mismatch between hmultFT0vskstar_pt and hmultpp136" << endl;
        }

        // TCanvas *cChMultEta0p5 = new TCanvas("cChMultEta0p5", "cChMultEta0p5", 720, 720);
        // SetCanvasStyle(cChMultEta0p5, 0.14, 0.03, 0.01, 0.14);
        // gPad->SetLogy();
        // // gPad->SetLogx();
        // for (int i = 0; i < NumPercent; i++)
        // {
        //     SetHistoQA(hChMultEta0p5[i]);
        //     hChMultEta0p5[i]->GetXaxis()->SetRangeUser(0, 75);
        //     hChMultEta0p5[i]->SetMaximum(20e7);
        //     hChMultEta0p5[i]->SetLineColor(colors[i]);
        //     hChMultEta0p5[i]->SetFillColor(colors[i]);
        //     hChMultEta0p5[i]->SetFillStyle(3001);
        //     hChMultEta0p5[i]->GetXaxis()->SetTitle("dN/d#eta in |#eta|<0.5");
        //     hChMultEta0p5[i]->GetYaxis()->SetTitle("Entries");
        //     hChMultEta0p5[i]->Draw("same");
        // }
        // leg->Draw();

        // //Get the K* pT for the given multiplicity classes
        // TCanvas *cKstarPt = new TCanvas("cKstarPt", "Kstar pT vs Multiplicity", 720, 720);
        // SetCanvasStyle(cKstarPt, 0.14, 0.03, 0.01, 0.14);
        // gPad->SetLogy();

        highBin = hmultpp136->GetNbinsX();
        TH1D *hKstarPt[NumPercent];
        const double deltaY = 1.0; // Acceptance |y| < 0.5

        for (int i = 0; i < NumPercent; i++)
        {
            int lowBin = regionEdges[i];
            hmultFT0vskstar_pt->GetXaxis()->SetRange(lowBin, highBin);
            hKstarPt[i] = (TH1D *)hmultFT0vskstar_pt->ProjectionY(Form("hKstarPt_%d", i));

            const double classEvents = hRegions[i]->Integral(1, hRegions[i]->GetNbinsX());
            const double rawYield = hKstarPt[i]->Integral(1, hKstarPt[i]->GetNbinsX());
            const double dNdy = (classEvents > 0.0) ? rawYield / (classEvents * deltaY) : 0.0;

            gMeanPt->SetPoint(i, hChMultEta0p5[i]->GetMean(), hKstarPt[i]->GetMean());
            gMeanPt->SetPointError(i, 0.0, 0.0);
            gdNdy->SetPoint(i, hChMultEta0p5[i]->GetMean(), dNdy);
            gdNdy->SetPointError(i, 0.0, 0.0);

            fOutput->cd();

            // hKstarPt[i]->Write(Form("hKstarPt_%d", i));
            hKstarPt[i]->Write(Form("h%sPt_%d", saveNames[ipar].c_str(), i));

            cout << "dN_{ch}/d#eta " << hChMultEta0p5[i]->GetMean()
                 << " : <pT> = " << hKstarPt[i]->GetMean()
                 //  << ", raw yield = " << rawYield
                 //  << ", N_events(class) = " << classEvents
                 << ", dN/dy = " << dNdy << endl;
            highBin = lowBin - 1;
        }
        gMeanPt->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5>}");
        gMeanPt->GetYaxis()->SetTitle("<#it{p}_{T}>");
        gdNdy->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5>}");
        gdNdy->GetYaxis()->SetTitle("dN/dy");

        gMeanPt->Write(Form("gMeanpT_%s", saveNames[ipar].c_str()));
        gdNdy->Write(Form("gMeanYield_%s", saveNames[ipar].c_str()));
    }
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.05);
    pad2->SetRightMargin(0.05);
    pad2->SetBottomMargin(0.35);
    pad1->SetLeftMargin(0.175);
    pad2->SetLeftMargin(0.175);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.003);
    pad2->SetTopMargin(0.04);

    // Set ticks on individual pads
    pad1->SetTicks(1, 1);
    pad2->SetTicks(1, 1);
}