#include <iostream>
#include "src/style.h"

using namespace std;

void efficiency_phi()
{
    gStyle->SetOptStat(0);
    TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Efficiency_phi.root", "READ");
    // TFile *fileeff = new TFile("/home/sawan/Downloads/AnalysisResults1.root", "READ");

    if (fileeff->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }
    double pT_bins[] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0};
    // double pT_bins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 10.0};
    const int Npt = sizeof(pT_bins) / sizeof(pT_bins[0]) - 1;

    const string genpath = "kstarqa_NoRCT/hInvMass";
    const string recpath = "kstarqa_NoRCT/hInvMass";

    THnSparseF *hSpraseGen = (THnSparseF *)fileeff->Get(Form("%s/hk892GenpTCalib1", genpath.c_str()));
    THnSparseF *hSparseRec = (THnSparseF *)fileeff->Get(Form("%s/h2KstarRecptCalib2", recpath.c_str()));
    if (hSpraseGen == nullptr || hSparseRec == nullptr)
    {
        cout << "Error reading efficiency histograms" << endl;
        return;
    }

    ////============Event loss histograms (method 1)===================
    TH1F *hAllGenColl = (TH1F *)fileeff->Get(Form("%s/MCcorrections/MultiplicityGen", genpath.c_str()));
    TH1F *hAllGenColl1Rec = (TH1F *)fileeff->Get(Form("%s/MCcorrections/MultiplicityRec", recpath.c_str()));
    if (hAllGenColl == nullptr || hAllGenColl1Rec == nullptr)
    {
        cout << "Error reading event loss histograms" << endl;
        return;
    }

    // ////=============Event loss histograms (method 2)===================
    // TH1F *hAllGenColl = (TH1F *)fileeff->Get("kstarqa/eventSelection/eventsCheckGen");

    //// ============Signal loss histogram (method 1)===================
    TH2F *hAllGenKstar = (TH2F *)fileeff->Get(Form("%s/MCcorrections/hSignalLossDenominator", genpath.c_str()));
    TH2F *hAllGenKstar1Rec = (TH2F *)fileeff->Get(Form("%s/MCcorrections/hSignalLossNumerator", recpath.c_str()));

    // // //// ============Signal loss histogram (method 2)===================
    // TH2F *hAllGenKstar = (TH2F *)fileeff->Get(Form("%s/CorrFactors/h2dGenKstar", genpath.c_str()));
    // THnSparse *hAllGenKstar1Rec = (THnSparse *)fileeff->Get(Form("%s/hk892GenpT2", recpath.c_str()));

    if (hAllGenKstar == nullptr || hAllGenKstar1Rec == nullptr)
    {
        cout << "Error reading signal loss histograms" << endl;
        return;
    }
    int multlow{0}, multhigh{100};

    TH1D *h1gen = hSpraseGen->Projection(0, "E");
    TH1D *h1rec = hSparseRec->Projection(0, "E");
    TH1F *heff = new TH1F("Efficiency", "Efficiency", Npt, pT_bins);
    TH1F *hSignalLoss = new TH1F("SignalLoss", "Signal Loss", Npt, pT_bins);
    TH1F *heventloss = new TH1F("EventLoss", "Event Loss", Npt, pT_bins);
    TH1F *hRatioEvBySig;

    int lowbinMultGen = hSpraseGen->GetAxis(1)->FindBin(multlow + 1e-3);
    int highbinMultGen = hSpraseGen->GetAxis(1)->FindBin(multhigh - 1e-3);
    hSpraseGen->GetAxis(1)->SetRange(lowbinMultGen, highbinMultGen);

    int lowbinMultRec = hSparseRec->GetAxis(1)->FindBin(multlow + 1e-3);
    int highbinMultRec = hSparseRec->GetAxis(1)->FindBin(multhigh - 1e-3);
    hSparseRec->GetAxis(1)->SetRange(lowbinMultRec, highbinMultRec);

    //====================Event loss calculations (method 1)==================
    double eventLossNum = hAllGenColl1Rec->Integral(hAllGenColl1Rec->GetXaxis()->FindBin(multlow + 0.001), hAllGenColl1Rec->GetXaxis()->FindBin(multhigh - 0.001)); // without event splitting consideration
    double eventLossDen = hAllGenColl->Integral(hAllGenColl->GetXaxis()->FindBin(multlow + 0.001), hAllGenColl->GetXaxis()->FindBin(multhigh - 0.001));

    // //====================Event loss calculations (method 2)==================
    // double eventLossNum = hAllGenColl->GetBinContent(5);
    // double eventLossDen = hAllGenColl->GetBinContent(4);

    //===================Signal loss histograms (method 1)=======================
    TH1F *hSignalLossNumPt = (TH1F *)hAllGenKstar1Rec->ProjectionX("SignalLossNumPt", hAllGenKstar1Rec->GetYaxis()->FindBin(multlow + 0.01), hAllGenKstar1Rec->GetYaxis()->FindBin(multhigh - 0.01));
    TH1F *hSignalLossDenPt = (TH1F *)hAllGenKstar->ProjectionX("SignalLossDenPt", hAllGenKstar->GetYaxis()->FindBin(multlow + 0.01), hAllGenKstar->GetYaxis()->FindBin(multhigh - 0.01));

    // //===================Signal loss histograms (method 2)=======================
    // hAllGenKstar1Rec->GetAxis(1)->SetRange(hAllGenKstar1Rec->GetAxis(1)->FindBin(multlow + 0.01), hAllGenKstar1Rec->GetAxis(1)->FindBin(multhigh - 0.01));
    // TH1F *hSignalLossNumPt = (TH1F *)hAllGenKstar1Rec->Projection(0, "E");
    // TH1F *hSignalLossDenPt = (TH1F *)hAllGenKstar->ProjectionY("SignalLossDenPt", hAllGenKstar->GetXaxis()->FindBin(multlow + 0.01), hAllGenKstar->GetXaxis()->FindBin(multhigh - 0.01));

    if (hSignalLossNumPt == nullptr || hSignalLossDenPt == nullptr)
    {
        cout << "Error projecting signal loss histograms" << endl;
        return;
    }

    for (int i = 0; i < heff->GetNbinsX(); i++)
    {
        double lowpt = pT_bins[i];
        double highpt = pT_bins[i + 1];

        double nrec = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt));
        double ngen = h1gen->Integral(h1gen->GetXaxis()->FindBin(lowpt), h1gen->GetXaxis()->FindBin(highpt));
        double efficiency = nrec / ngen;
        double efficiencyerr = sqrt(abs(((nrec + 1) / (ngen + 2)) * ((nrec + 2) / (ngen + 3) - (nrec + 1) / (ngen + 2))));

        heff->SetBinContent(i + 1, efficiency);
        heff->SetBinError(i + 1, efficiencyerr);

        //================Event loss calculations==================
        heventloss->SetBinContent(i + 1, eventLossNum / eventLossDen);
        heventloss->SetBinError(i + 1, 0);

        //================Signal loss calculations==================
        double signalLossNum = hSignalLossNumPt->Integral(hSignalLossNumPt->GetXaxis()->FindBin(lowpt + 0.001), hSignalLossNumPt->GetXaxis()->FindBin(highpt - 0.001));
        double signalLossDen = hSignalLossDenPt->Integral(hSignalLossDenPt->GetXaxis()->FindBin(lowpt + 0.001), hSignalLossDenPt->GetXaxis()->FindBin(highpt - 0.001));
        double signalLoss = signalLossNum / signalLossDen;

        hSignalLoss->SetBinContent(i + 1, signalLoss);
        hSignalLoss->SetBinError(i + 1, 0);
    }
    hRatioEvBySig = (TH1F *)heventloss->Clone();
    hRatioEvBySig->Divide(hSignalLoss);

    TString savePath = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra/";

    // Plot other plots for all multiplicity bins
    TCanvas *cefficiency = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cefficiency, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(heff);
    heff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    heff->GetYaxis()->SetTitle("Acceptance x Efficiency");
    heff->GetYaxis()->SetTitleOffset(1.6);
    heff->SetMaximum(0.5);
    heff->SetMarkerSize(1.2);
    heff->Draw("pe");

    TCanvas *cEventLoss = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cEventLoss, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(heventloss);
    heventloss->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    heventloss->GetYaxis()->SetTitle("Event and signal Loss");
    heventloss->GetYaxis()->SetTitleOffset(1.6);
    heventloss->SetStats(0);
    heventloss->SetLineStyle(2);
    heventloss->GetYaxis()->SetRangeUser(0.25, 1.18);
    heventloss->SetLineWidth(3);
    heventloss->Draw("l");
    SetHistoQA(hSignalLoss);
    hSignalLoss->SetMarkerSize(1.5);
    hSignalLoss->SetMarkerColor(kRed);
    hSignalLoss->Draw("pe same");
    TLegend *legend = new TLegend(0.55, 0.72, 0.91, 0.96);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->AddEntry((TObject *)0, "ALICE", "");
    legend->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legend->AddEntry((TObject *)0, Form("FT0M: %d-%d%%", multlow, multhigh), "");
    legend->AddEntry(heventloss, "Event Loss", "l");
    legend->AddEntry(hSignalLoss, "Signal Loss", "p");
    legend->Draw();
    cEventLoss->SaveAs(savePath + "/EventAndSignalLoss_phi_mult0-100.png");

    TCanvas *cEventBySignalLoss = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cEventBySignalLoss, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hRatioEvBySig);
    hRatioEvBySig->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatioEvBySig->GetYaxis()->SetTitle("Event Loss / Signal Loss");
    hRatioEvBySig->GetYaxis()->SetTitleOffset(1.6);
    hRatioEvBySig->SetStats(0);
    hRatioEvBySig->SetMaximum(1.06);
    hRatioEvBySig->SetMinimum(0.6);
    hRatioEvBySig->SetMarkerSize(1.5);
    hRatioEvBySig->Draw("pe");
    TLegend *legend2 = new TLegend(0.4, 0.67, 0.91, 0.93);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(0.035);
    legend2->AddEntry((TObject *)0, "ALICE", "");
    legend2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legend2->AddEntry((TObject *)0, Form("FT0M: %d-%d%%", multlow, multhigh), "");
    legend2->AddEntry(hRatioEvBySig, "Event Loss / Signal Loss", "p");
    legend2->Draw();
    cEventBySignalLoss->SaveAs(savePath + "/EventBySignalLoss_phi_mult0-100.png");

    TFile *fLossPhi = new TFile(savePath + "/Loss_phi_mult0-100.root", "RECREATE");
    heff->Write("hEff");
    heventloss->Write("hEventLoss");
    hSignalLoss->Write("hSignalLoss");
    hRatioEvBySig->Write("hRatio");
    fLossPhi->Close();
}
