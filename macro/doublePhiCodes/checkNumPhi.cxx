#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"

TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const std::string &name);

void checkNumPhi()
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TString savepath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/TetraquarkFit";

    TFile *fInput26 = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResults26_latest.root");
    TFile *fInput25 = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/LHC25/AnalysisResults25_aiam.root"); // opti8, 25 ai+am

    TH1F *hNumPhi26 = GetHisto<TH1F>(fInput26, "doublephimeson/NPhiPerEvent");
    TH1F *hNumPhi25 = GetHisto<TH1F>(fInput25, "doublephimeson/NPhiPerEvent");

    // hNumPhi26->Scale(1.0 / hNumPhi26->Integral());
    // hNumPhi25->Scale(1.0 / hNumPhi25->Integral());

    TCanvas *cNumPhi = new TCanvas("cNumPhi", "Number of Phi Mesons per Event", 720, 720);
    SetCanvasStyle(cNumPhi, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hNumPhi26);
    SetHistoQA(hNumPhi25);
    hNumPhi26->GetXaxis()->SetTitle("Number of #phi mesons");
    hNumPhi26->GetYaxis()->SetTitle("Counts");
    hNumPhi25->SetLineColor(kRed);
    hNumPhi26->SetLineColor(kBlue);
    gPad->SetLogy();
    hNumPhi26->Draw("HIST");
    hNumPhi25->Draw("HIST SAME");

    TLegend *legend = new TLegend(0.45, 0.75, 0.9, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(hNumPhi26, "LHC26_skimmed", "l");
    legend->AddEntry(hNumPhi25, "LHC25(ai+am)_skimmed", "l");
    legend->Draw();
    cNumPhi->SaveAs((savepath + "/NumPhiPerEvent_Compare.png").Data());
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
