#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"

TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const std::string &name);

void doublePhiSimpleOpti6()
{
    ////=====Pair=========
    TFile *fInput = OpenFile("/home/sawan/alice/practice/OutputDoublePhi/Pair/processPairOpti6/AnalysisResults.root");

    THnSparseF *hUnlike = GetHisto<THnSparseF>(fInput, "doublephimeson/SEMassUnlike_VertexVars");

    // It has 8 axes:
    // 0: Invariant Mass
    // 1: Pt
    // 2: deltaM
    // 3: Phi1 Mass
    // 4: Phi2 Mass
    // 5: Decay length
    // 6: fit Chi2/NDF
    // 7: RMS DCA signal

    TH1D *hPt = hUnlike->Projection(1, "E");
    TH1D *hDeltaM = hUnlike->Projection(2, "E");
    TH1D *hPhi1Mass = hUnlike->Projection(3, "E");
    TH1D *hPhi2Mass = hUnlike->Projection(4, "E");
    TH1D *hDecayLength = hUnlike->Projection(5, "E");
    TH1D *hFitChi2NDF = hUnlike->Projection(6, "E");
    TH1D *hRMSDCA = hUnlike->Projection(7, "E");

    int lowpT = hUnlike->GetAxis(1)->FindBin(9.0 + 0.001);
    int highpT = hUnlike->GetAxis(1)->FindBin(100.0 - 0.001);

    int lowDeltaM = hUnlike->GetAxis(2)->FindBin(0.0 + 0.0001);
    int highDeltaM = hUnlike->GetAxis(2)->FindBin(0.005 - 0.0001);

    int lowPhi1Mass = hUnlike->GetAxis(3)->FindBin(1.005 + 0.0001);
    int highPhi1Mass = hUnlike->GetAxis(3)->FindBin(1.035 - 0.0001);

    int lowPhi2Mass = hUnlike->GetAxis(4)->FindBin(1.005 + 0.0001);
    int highPhi2Mass = hUnlike->GetAxis(4)->FindBin(1.035 - 0.0001);

    int lowDecayLength = hUnlike->GetAxis(5)->FindBin(0.0 + 0.0001);
    int highDecayLength = hUnlike->GetAxis(5)->FindBin(10.0 - 0.0001);

    int lowFitChi2NDF = hUnlike->GetAxis(6)->FindBin(0.0 + 0.0001);
    int highFitChi2NDF = hUnlike->GetAxis(6)->FindBin(10.0 - 0.0001);

    int lowRMSDCA = hUnlike->GetAxis(7)->FindBin(0.0 + 0.0001);
    int highRMSDCA = hUnlike->GetAxis(7)->FindBin(0.5 - 0.0001);

    hUnlike->GetAxis(1)->SetRange(lowpT, highpT);
    hUnlike->GetAxis(2)->SetRange(lowDeltaM, highDeltaM);
    // hUnlike->GetAxis(3)->SetRange(lowPhi1Mass, highPhi1Mass);
    // hUnlike->GetAxis(4)->SetRange(lowPhi2Mass, highPhiMass);
    // hUnlike->GetAxis(5)->SetRange(lowDecayLength, highDecayLength);
    // hUnlike->GetAxis(6)->SetRange(lowFitChi2NDF, highFitChi2NDF);
    // hUnlike->GetAxis(7)->SetRange(lowRMSDCA, highRMSDCA);

    TH1D *hInvMass = hUnlike->Projection(0, "E");
    SetHistoQA(hInvMass);
    hInvMass->Rebin(8);
    hInvMass->GetXaxis()->SetRangeUser(2.5, 2.9);
    hInvMass->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hInvMass->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/#it{c}^{2}", hInvMass->GetBinWidth(1) * 1000));

    TCanvas *cInvMass = new TCanvas("cInvMass", "Invariant Mass", 720, 720);
    SetCanvasStyle(cInvMass, 0.15, 0.03, 0.05, 0.15);
    hInvMass->Draw("pe");
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
