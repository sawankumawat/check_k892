#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"

TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const std::string &name);

void phiphiCorrelation()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/TetraquarkFit";
    int rebinFactor = 12;

    ////=====New===========
    ////====2026 data========
    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResults26_latest.root");

    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassPhiPhiRefitted");

    TH1D *hPt = hUnlike->Projection(1, "E");
    TH1D *hDeltaM = hUnlike->Projection(2, "E");

    int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowDeltaM = hUnlike->GetAxis(2)->FindBin(0.01 + 0.00001);
    int highDeltaM = hUnlike->GetAxis(2)->FindBin(100.0 - 0.00001);

    int lowPhiMassCut = hUnlike->GetAxis(5)->FindBin(1.01 + 0.00001);
    int highPhiMassCut = hUnlike->GetAxis(5)->FindBin(1.03 - 0.00001);

    int lowWithoutPhiCut = hUnlike->GetAxis(5)->FindBin(0.0 + 0.00001);
    int highWithoutPhiCut = hUnlike->GetAxis(5)->FindBin(3.0 - 0.00001);

    hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    hUnlike->GetAxis(2)->SetRange(lowDeltaM, highDeltaM);
    hUnlike->GetAxis(5)->SetRange(lowWithoutPhiCut, lowPhiMassCut); // Phi1 < 1.01 GeV/c^2 (left corner region)
    hUnlike->GetAxis(6)->SetRange(lowWithoutPhiCut, lowPhiMassCut); // Phi2 < 1.01 GeV/c^2 (left corner region)

    TH1D *hInvMassBkgLeft = hUnlike->Projection(0, "E");
    SetHistoQA(hInvMassBkgLeft);

    hUnlike->GetAxis(5)->SetRange(lowWithoutPhiCut, lowPhiMassCut);   // Phi1 < 1.01 GeV/c^2 (top corner region)
    hUnlike->GetAxis(6)->SetRange(highPhiMassCut, highWithoutPhiCut); // Phi2 > 1.025 GeV/c^2 (top corner region)

    TH1D *hInvMassBkgTop = hUnlike->Projection(0, "E");
    SetHistoQA(hInvMassBkgTop);

    hUnlike->GetAxis(5)->SetRange(highPhiMassCut, highWithoutPhiCut); // Phi1 > 1.025 GeV/c^2 (right corner region)
    hUnlike->GetAxis(6)->SetRange(highPhiMassCut, highWithoutPhiCut); // Phi2 > 1.025 GeV/c^2 (right corner region)

    TH1D *hInvMassBkgRight = hUnlike->Projection(0, "E");
    SetHistoQA(hInvMassBkgRight);

    hUnlike->GetAxis(5)->SetRange(highPhiMassCut, highWithoutPhiCut); // Phi1 > 1.025 GeV/c^2 (bottom corner region)
    hUnlike->GetAxis(6)->SetRange(lowWithoutPhiCut, lowPhiMassCut);   // Phi2 < 1.01 GeV/c^2 (bottom corner region)

    TH1D *hInvMassBkgBottom = hUnlike->Projection(0, "E");
    SetHistoQA(hInvMassBkgBottom);

    // Add all 4 histograms in hInvMassBkg
    TH1D *hInvMassBkg = (TH1D *)hInvMassBkgLeft->Clone("hInvMassBkg");
    hInvMassBkg->Add(hInvMassBkgTop);
    hInvMassBkg->Add(hInvMassBkgRight);
    hInvMassBkg->Add(hInvMassBkgBottom);
    hInvMassBkg->Rebin(rebinFactor);

    TCanvas *cPhiPhiCorrelation = new TCanvas("cPhiPhiCorrelation", "Phi-Phi Correlation", 720, 720);
    SetCanvasStyle(cPhiPhiCorrelation, 0.15, 0.03, 0.05, 0.15);
    hUnlike->GetAxis(5)->SetRange(lowWithoutPhiCut, highWithoutPhiCut); 
    hUnlike->GetAxis(6)->SetRange(lowWithoutPhiCut, highWithoutPhiCut);   // Phi2 < 1.01 GeV/c^2 (bottom corner region)
    TH2D *hPhiPhiCorrelation = hUnlike->Projection(5, 6, "E");
    SetHistoQA(hPhiPhiCorrelation);
    hPhiPhiCorrelation->GetXaxis()->SetTitle("M_{#phi} (GeV/#it{c}^{2})");
    hPhiPhiCorrelation->GetYaxis()->SetTitle("M_{#phi} (GeV/#it{c}^{2})");
    hPhiPhiCorrelation->GetZaxis()->SetTitle("Counts");
    hPhiPhiCorrelation->Draw("colz");
    cPhiPhiCorrelation->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/TetraquarkFit/PhiPhiCorrelation.png");

    // Compare background with Non-SS template
    TCanvas *cBkgTemplate = new TCanvas("cBkgTemplate", "Background Template", 720, 720);
    SetCanvasStyle(cBkgTemplate, 0.15, 0.03, 0.05, 0.15);
    // hInvMassBkg->SetMinimum(1e3);
    // hInvMassBkg->SetMaximum(11.8e3);
    hInvMassBkg->SetMinimum(0.6e3);
    hInvMassBkg->SetMaximum(6.8e3);
    hInvMassBkg->SetMarkerColor(kRed);
    hInvMassBkg->SetLineColor(kRed);
    hInvMassBkg->GetXaxis()->SetRangeUser(2.41, 2.95);
    hInvMassBkg->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    // hInvMassBkg->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hInvMassBkg->GetBinWidth(1) * 1000));
    hInvMassBkg->GetYaxis()->SetTitle("Counts");
    hInvMassBkg->Draw("pe");

    TFile *fTemplate = OpenFile("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template/PhiPhiBkgTemplate_BW_ExtendedFitRange.root");
    TH1D *hNonSSTemplate = GetHisto<TH1D>(fTemplate, "h_N_BBOnly");
    SetHistoQA(hNonSSTemplate);
    hNonSSTemplate->SetMarkerColor(kBlue + 1);
    hNonSSTemplate->SetLineColor(kBlue + 1);
    hNonSSTemplate->SetLineWidth(2);
    hNonSSTemplate->Scale(hInvMassBkg->Integral(hInvMassBkg->GetXaxis()->FindBin(2.4), hInvMassBkg->GetXaxis()->FindBin(2.9)) / hNonSSTemplate->Integral(hNonSSTemplate->GetXaxis()->FindBin(2.4), hNonSSTemplate->GetXaxis()->FindBin(2.9)));
    hNonSSTemplate->Draw("pe same");

    TLegend *legendTemplate = new TLegend(0.5, 0.79, 0.9, 0.92);
    legendTemplate->SetBorderSize(0);
    legendTemplate->SetFillStyle(0);
    legendTemplate->SetTextFont(42);
    legendTemplate->SetTextSize(0.035);
    legendTemplate->AddEntry(hInvMassBkg, "K^{+}K^{-} background", "p");
    legendTemplate->AddEntry(hNonSSTemplate, "BB (K^{+}K^{-}K^{+}K^{-}) template", "p");
    legendTemplate->Draw();
    
    //Add a legend to show that K+K- background is corner region m1<1.01 , m2 <1.01, m1<1.01, m2>1.025, m1>1.025, m2>1.025, m1>1.025, m2<1.01
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(22);
    latex->SetTextSize(0.03);
    latex->DrawLatex(0.3, 0.45, "K^{+}K^{-} background: corner region:");
    latex->SetTextFont(42);
    latex->DrawLatex(0.3, 0.4, "m_{1} < 1.01 GeV/#it{c}^{2}, m_{2} < 1.01 GeV/#it{c}^{2}");
    latex->DrawLatex(0.3, 0.35, "m_{1} < 1.01 GeV/#it{c}^{2}, m_{2} > 1.03 GeV/#it{c}^{2}");
    latex->DrawLatex(0.3, 0.3, "m_{1} > 1.03 GeV/#it{c}^{2}, m_{2} > 1.03 GeV/#it{c}^{2}");
    latex->DrawLatex(0.3, 0.25, "m_{1} > 1.03 GeV/#it{c}^{2}, m_{2} < 1.01 GeV/#it{c}^{2}");
    latex->DrawLatex(0.3, 0.2, "#Delta M > 0.01 GeV/#it{c}^{2}, #it{p}_{T}^{#phi#phi} > 9 GeV/#it{c}");
    cBkgTemplate->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/TemplateDeltaMCompare/TemplateAwayRegionCompare.png");
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
