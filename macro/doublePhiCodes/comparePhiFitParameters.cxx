#include <iostream>
#include "src/style.h"

TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const std::string &name);

void comparePhiFitParameters()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/PhiInvMass/PeriodWise";

    string suffixes[] = {"25ac", "25ah", "25ai", "25am", "26ac", "26ad", "26ae", "26af", "26ag", "26ai"};
    const int totalFiles = sizeof(suffixes) / sizeof(suffixes[0]);
    int markerStyles[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
    int colors[] = {kRed, kBlue, kGreen + 2, kMagenta, kCyan + 2, kOrange + 7, kViolet + 1, kTeal + 3, kPink + 1, kAzure + 1};

    TCanvas *cMassVsPt = new TCanvas("cMassVsPt", "Phi Mass vs Pt", 1080, 720);
    SetCanvasStyle(cMassVsPt, 0.15, 0.13, 0.05, 0.15);
    TCanvas *cResolutionVsPt = new TCanvas("cResolutionVsPt", "Phi Resolution vs Pt", 1080, 720);
    SetCanvasStyle(cResolutionVsPt, 0.13, 0.13, 0.05, 0.15);
    TCanvas *cPurityVsPt = new TCanvas("cPurityVsPt", "Phi Purity vs Pt", 1080, 720);
    SetCanvasStyle(cPurityVsPt, 0.13, 0.13, 0.05, 0.15);

    TLegend *legend = new TLegend(0.9, 0.25, 1.03, 0.75);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);

    for (int ifiles = 0; ifiles < totalFiles; ifiles++)
    {
        string suffix = suffixes[ifiles];
        TFile *fInput = OpenFile("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/PhiInvMass/PeriodWise/PhiParams" + suffix + ".root");
        TGraphErrors *hPhiMass = GetHisto<TGraphErrors>(fInput, "gMassVsPt");
        TGraphErrors *hPhiResolution = GetHisto<TGraphErrors>(fInput, "gResolutionVsPt");
        TGraphErrors *hPhiPurity = GetHisto<TGraphErrors>(fInput, "gPurity");

        SetGraphErrorStyle(hPhiMass);
        SetGraphErrorStyle(hPhiResolution);
        SetGraphErrorStyle(hPhiPurity);
        hPhiMass->SetMarkerStyle(markerStyles[ifiles]);
        hPhiMass->SetMarkerColor(colors[ifiles]);
        hPhiMass->SetLineColor(colors[ifiles]);
        hPhiMass->GetYaxis()->SetRangeUser(1.0185, 1.0203);

        hPhiResolution->SetMarkerStyle(markerStyles[ifiles]);
        hPhiResolution->SetMarkerColor(colors[ifiles]);
        hPhiResolution->SetLineColor(colors[ifiles]);
        hPhiResolution->GetYaxis()->SetTitleOffset(1.3);

        hPhiPurity->GetYaxis()->SetTitleOffset(1.3);
        hPhiPurity->SetMarkerStyle(markerStyles[ifiles]);
        hPhiPurity->SetMarkerColor(colors[ifiles]);
        hPhiPurity->SetLineColor(colors[ifiles]);

        // Mass
        cMassVsPt->cd();
        hPhiMass->Draw(ifiles == 0 ? "AP" : "P SAME");

        // Resolution
        cResolutionVsPt->cd();
        hPhiResolution->Draw(ifiles == 0 ? "AP" : "P SAME");

        // Purity
        cPurityVsPt->cd();
        hPhiPurity->Draw(ifiles == 0 ? "AP" : "P SAME");
        legend->AddEntry(hPhiPurity, Form("%s", suffix.c_str()), "p");
    }
    cMassVsPt->cd();
    legend->Draw();
    cResolutionVsPt->cd();
    legend->Draw();
    cPurityVsPt->cd();
    legend->Draw();

    cMassVsPt->SaveAs(savepath + "/PhiMassVsPt_PeriodWise.png");
    cResolutionVsPt->SaveAs(savepath + "/PhiResolutionVsPt_PeriodWise.png");
    cPurityVsPt->SaveAs(savepath + "/PhiPurityVsPt_PeriodWise.png");
}

//==============End of the main code==================

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

template <typename T>
T *GetHisto(TFile *f, const std::string &name)
{
    T *histo = dynamic_cast<T *>(f->Get(name.c_str()));

    if (!histo)
    {
        std::cout << "Error: histo " << name
                  << " not found in file " << f->GetName() << std::endl;
        return nullptr;
    }

    return histo;
}
