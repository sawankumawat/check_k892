#include <iostream>
#include <array>
#include <TFile.h>
#include "../src/style.h"

void phiMassShiftMC()
{
    gStyle->SetOptStat(0);

    // TString dataPath = "/home/sawan/Storage/check_k892/phiMC";
    TString dataPath = "/home/sawan/Downloads";
    TFile *fLHC25ac = new TFile(dataPath + "/AnalysisResults25_ac.root");
    TFile *fLHC25ai = new TFile(dataPath + "/AnalysisResults25_ai.root");
    if (!fLHC25ac || !fLHC25ai)
    {
        std::cerr << "Error: Could not open one or both ROOT files." << std::endl;
        return;
    }

    ////=========Calculated from toy model the mass shift for 2025 vs 2026 dataset=========
    constexpr int kNBins2 = 14;
    constexpr double kPtEdges2[kNBins2 + 1] = {0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    constexpr double kEpsilon[kNBins2] = {0.0024146539379, 0.00762912304405, 0.00817386742107, 0.00982766792836, 0.0106788436324, 0.0112954252224, 0.0115583787246, 0.0117416086153, 0.012422150094, 0.011968584565, 0.0125561781092, 0.0139372556405, 0.014995530318, 0.0165687496159};
    ////====================================================================================

    TH1D *hShiftKaonFromPhiToy_25_26 = new TH1D("hShiftKaonFromPhiToy_25_26", "hShiftKaonFromPhiToy_25_26", kNBins2, kPtEdges2);
    for (int ipt = 0; ipt < kNBins2; ++ipt)
    {
        hShiftKaonFromPhiToy_25_26->SetBinContent(ipt + 1, -kEpsilon[ipt]);
        hShiftKaonFromPhiToy_25_26->SetBinError(ipt + 1, 0.0);
    }

    constexpr int kNBins = 13;
    constexpr double kPtEdges[kNBins + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0};

    TH2D *h2D_dpx_vs_px_AllKaons = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpx_vs_px_AllKaons");
    TH2D *h2D_dpy_vs_py_AllKaons = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpy_vs_py_AllKaons");
    TH2D *h2D_dpz_vs_pz_AllKaons = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpz_vs_pz_AllKaons");
    TH2D *h2D_dpt_vs_pt_AllKaons = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpt_vs_pt_AllKaons");

    TH2D *h2D_dpx_vs_px_KaonFromPhi = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpx_vs_px_KaonFromPhi");
    TH2D *h2D_dpy_vs_py_KaonFromPhi = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpy_vs_py_KaonFromPhi");
    TH2D *h2D_dpz_vs_pz_KaonFromPhi = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpz_vs_pz_KaonFromPhi");
    TH2D *h2D_dpt_vs_pt_KaonFromPhi = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpt_vs_pt_KaonFromPhi");

    TH2D *h2D_dpx_vs_px_Phi = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpx_vs_px_Phi");
    TH2D *h2D_dpy_vs_py_Phi = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpy_vs_py_Phi");
    TH2D *h2D_dpz_vs_pz_Phi = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpz_vs_pz_Phi");
    TH2D *h2D_dpt_vs_pt_Phi = (TH2D *)fLHC25ac->Get("kstarqa/hShift/h2D_dpt_vs_pt_Phi");

    vector<array<double, 3>> pxpypzShift_AllKaons;
    vector<array<double, 3>> pxpypzShift_KaonFromPhi;
    vector<array<double, 3>> pxpypzShift_Phi;

    TH1D *hpxShift_AllKaons = new TH1D("hpxShift_AllKaons", "hpxShift_AllKaons", kNBins, kPtEdges);
    TH1D *hpyShift_AllKaons = new TH1D("hpyShift_AllKaons", "hpyShift_AllKaons", kNBins, kPtEdges);
    TH1D *hpzShift_AllKaons = new TH1D("hpzShift_AllKaons", "hpzShift_AllKaons", kNBins, kPtEdges);
    TH1D *hpxShift_KaonFromPhi = new TH1D("hpxShift_KaonFromPhi", "hpxShift_KaonFromPhi", kNBins, kPtEdges);
    TH1D *hpyShift_KaonFromPhi = new TH1D("hpyShift_KaonFromPhi", "hpyShift_KaonFromPhi", kNBins, kPtEdges);
    TH1D *hpzShift_KaonFromPhi = new TH1D("hpzShift_KaonFromPhi", "hpzShift_KaonFromPhi", kNBins, kPtEdges);
    TH1D *hpxShift_Phi = new TH1D("hpxShift_Phi", "hpxShift_Phi", kNBins, kPtEdges);
    TH1D *hpyShift_Phi = new TH1D("hpyShift_Phi", "hpyShift_Phi", kNBins, kPtEdges);
    TH1D *hpzShift_Phi = new TH1D("hpzShift_Phi", "hpzShift_Phi", kNBins, kPtEdges);

    for (int ipt = 0; ipt < kNBins; ++ipt)
    {
        double ptMin = kPtEdges[ipt] + 1e-6;
        double ptMax = kPtEdges[ipt + 1] - 1e-6;

        // Calculate the mean shift for All Kaons
        TH1D *h_dpx_AllKaons = h2D_dpx_vs_px_AllKaons->ProjectionX(Form("h_dpx_AllKaons_pt%d", ipt), h2D_dpx_vs_px_AllKaons->GetYaxis()->FindBin(ptMin), h2D_dpx_vs_px_AllKaons->GetYaxis()->FindBin(ptMax));
        TH1D *h_dpy_AllKaons = h2D_dpy_vs_py_AllKaons->ProjectionX(Form("h_dpy_AllKaons_pt%d", ipt), h2D_dpy_vs_py_AllKaons->GetYaxis()->FindBin(ptMin), h2D_dpy_vs_py_AllKaons->GetYaxis()->FindBin(ptMax));
        TH1D *h_dpz_AllKaons = h2D_dpz_vs_pz_AllKaons->ProjectionX(Form("h_dpz_AllKaons_pt%d", ipt), h2D_dpz_vs_pz_AllKaons->GetYaxis()->FindBin(ptMin), h2D_dpz_vs_pz_AllKaons->GetYaxis()->FindBin(ptMax));

        double mean_dpx_AllKaons = h_dpx_AllKaons->GetMean();
        double mean_dpx_AllKaons_error = h_dpx_AllKaons->GetMeanError();
        double mean_dpy_AllKaons = h_dpy_AllKaons->GetMean();
        double mean_dpy_AllKaons_error = h_dpy_AllKaons->GetMeanError();
        double mean_dpz_AllKaons = h_dpz_AllKaons->GetMean();
        double mean_dpz_AllKaons_error = h_dpz_AllKaons->GetMeanError();

        hpxShift_AllKaons->SetBinContent(ipt + 1, mean_dpx_AllKaons);
        hpyShift_AllKaons->SetBinContent(ipt + 1, mean_dpy_AllKaons);
        hpzShift_AllKaons->SetBinContent(ipt + 1, mean_dpz_AllKaons);
        hpxShift_AllKaons->SetBinError(ipt + 1, mean_dpx_AllKaons_error);
        hpyShift_AllKaons->SetBinError(ipt + 1, mean_dpy_AllKaons_error);
        hpzShift_AllKaons->SetBinError(ipt + 1, mean_dpz_AllKaons_error);

        pxpypzShift_AllKaons.push_back({mean_dpx_AllKaons, mean_dpy_AllKaons, mean_dpz_AllKaons});

        // Calculate the mean shift for Kaons from Phi
        TH1D *h_dpx_KaonFromPhi = h2D_dpx_vs_px_KaonFromPhi->ProjectionX(Form("h_dpx_KaonFromPhi_pt%d", ipt), h2D_dpx_vs_px_KaonFromPhi->GetYaxis()->FindBin(ptMin), h2D_dpx_vs_px_KaonFromPhi->GetYaxis()->FindBin(ptMax));
        TH1D *h_dpy_KaonFromPhi = h2D_dpy_vs_py_KaonFromPhi->ProjectionX(Form("h_dpy_KaonFromPhi_pt%d", ipt), h2D_dpy_vs_py_KaonFromPhi->GetYaxis()->FindBin(ptMin), h2D_dpy_vs_py_KaonFromPhi->GetYaxis()->FindBin(ptMax));
        TH1D *h_dpz_KaonFromPhi = h2D_dpz_vs_pz_KaonFromPhi->ProjectionX(Form("h_dpz_KaonFromPhi_pt%d", ipt), h2D_dpz_vs_pz_KaonFromPhi->GetYaxis()->FindBin(ptMin), h2D_dpz_vs_pz_KaonFromPhi->GetYaxis()->FindBin(ptMax));

        double mean_dpx_KaonFromPhi = h_dpx_KaonFromPhi->GetMean();
        double mean_dpy_KaonFromPhi = h_dpy_KaonFromPhi->GetMean();
        double mean_dpz_KaonFromPhi = h_dpz_KaonFromPhi->GetMean();
        double mean_dpx_KaonFromPhi_error = h_dpx_KaonFromPhi->GetMeanError();
        double mean_dpy_KaonFromPhi_error = h_dpy_KaonFromPhi->GetMeanError();
        double mean_dpz_KaonFromPhi_error = h_dpz_KaonFromPhi->GetMeanError();

        hpxShift_KaonFromPhi->SetBinContent(ipt + 1, mean_dpx_KaonFromPhi);
        hpyShift_KaonFromPhi->SetBinContent(ipt + 1, mean_dpy_KaonFromPhi);
        hpzShift_KaonFromPhi->SetBinContent(ipt + 1, mean_dpz_KaonFromPhi);
        hpxShift_KaonFromPhi->SetBinError(ipt + 1, mean_dpx_KaonFromPhi_error);
        hpyShift_KaonFromPhi->SetBinError(ipt + 1, mean_dpy_KaonFromPhi_error);
        hpzShift_KaonFromPhi->SetBinError(ipt + 1, mean_dpz_KaonFromPhi_error);

        pxpypzShift_KaonFromPhi.push_back({mean_dpx_KaonFromPhi, mean_dpy_KaonFromPhi, mean_dpz_KaonFromPhi});

        // Calculate the mean shift for Phi
        TH1D *h_dpx_Phi = h2D_dpx_vs_px_Phi->ProjectionX(Form("h_dpx_Phi_pt%d", ipt), h2D_dpx_vs_px_Phi->GetYaxis()->FindBin(ptMin), h2D_dpx_vs_px_Phi->GetYaxis()->FindBin(ptMax));
        TH1D *h_dpy_Phi = h2D_dpy_vs_py_Phi->ProjectionX(Form("h_dpy_Phi_pt%d", ipt), h2D_dpy_vs_py_Phi->GetYaxis()->FindBin(ptMin), h2D_dpy_vs_py_Phi->GetYaxis()->FindBin(ptMax));
        TH1D *h_dpz_Phi = h2D_dpz_vs_pz_Phi->ProjectionX(Form("h_dpz_Phi_pt%d", ipt), h2D_dpz_vs_pz_Phi->GetYaxis()->FindBin(ptMin), h2D_dpz_vs_pz_Phi->GetYaxis()->FindBin(ptMax));

        double mean_dpx_Phi = h_dpx_Phi->GetMean();
        double mean_dpy_Phi = h_dpy_Phi->GetMean();
        double mean_dpz_Phi = h_dpz_Phi->GetMean();
        double mean_dpx_Phi_error = h_dpx_Phi->GetMeanError();
        double mean_dpy_Phi_error = h_dpy_Phi->GetMeanError();
        double mean_dpz_Phi_error = h_dpz_Phi->GetMeanError();

        hpxShift_Phi->SetBinContent(ipt + 1, mean_dpx_Phi);
        hpyShift_Phi->SetBinContent(ipt + 1, mean_dpy_Phi);
        hpzShift_Phi->SetBinContent(ipt + 1, mean_dpz_Phi);
        hpxShift_Phi->SetBinError(ipt + 1, mean_dpx_Phi_error);
        hpyShift_Phi->SetBinError(ipt + 1, mean_dpy_Phi_error);
        hpzShift_Phi->SetBinError(ipt + 1, mean_dpz_Phi_error);

        pxpypzShift_Phi.push_back({mean_dpx_Phi, mean_dpy_Phi, mean_dpz_Phi});
    }

    TCanvas *cShiftAllKaons = new TCanvas("cShiftAllKaons", "cShiftAllKaons", 720, 720);
    SetCanvasStyle(cShiftAllKaons, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hpxShift_AllKaons);
    SetHistoQA(hpyShift_AllKaons);
    SetHistoQA(hpzShift_AllKaons);
    hpxShift_AllKaons->SetMarkerStyle(20);
    hpxShift_AllKaons->SetLineColor(kRed);
    hpxShift_AllKaons->SetMarkerColor(kRed);
    hpxShift_AllKaons->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hpxShift_AllKaons->GetYaxis()->SetTitle("<#Delta p_{x}> (GeV/c)");
    hpxShift_AllKaons->Draw("PE");
    hpyShift_AllKaons->SetMarkerStyle(21);
    hpyShift_AllKaons->SetLineColor(kBlue);
    hpyShift_AllKaons->SetMarkerColor(kBlue);
    hpyShift_AllKaons->Draw("PE SAME");
    hpzShift_AllKaons->SetMarkerStyle(22);
    hpzShift_AllKaons->SetLineColor(kGreen + 2);
    hpzShift_AllKaons->SetMarkerColor(kGreen + 2);
    hpzShift_AllKaons->Draw("PE SAME");

    TLegend *legShiftAllKaons = new TLegend(0.65, 0.78, 0.9, 0.9);
    legShiftAllKaons->SetBorderSize(0);
    legShiftAllKaons->SetFillStyle(0);
    legShiftAllKaons->AddEntry(hpxShift_AllKaons, "p_{x} shift", "lp");
    legShiftAllKaons->AddEntry(hpyShift_AllKaons, "p_{y} shift", "lp");
    legShiftAllKaons->AddEntry(hpzShift_AllKaons, "p_{z} shift", "lp");
    legShiftAllKaons->Draw();

    TCanvas *cShiftKaonFromPhi = new TCanvas("cShiftKaonFromPhi", "cShiftKaonFromPhi", 720, 720);
    SetCanvasStyle(cShiftKaonFromPhi, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hpxShift_KaonFromPhi);
    SetHistoQA(hpyShift_KaonFromPhi);
    SetHistoQA(hpzShift_KaonFromPhi);
    hpxShift_KaonFromPhi->SetMarkerStyle(20);
    hpxShift_KaonFromPhi->SetLineColor(kRed);
    hpxShift_KaonFromPhi->SetMarkerColor(kRed);
    hpxShift_KaonFromPhi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hpxShift_KaonFromPhi->GetYaxis()->SetTitle("<#Delta p_{x}> (GeV/c)");
    hpxShift_KaonFromPhi->Draw("PE");
    hpyShift_KaonFromPhi->SetMarkerStyle(21);
    hpyShift_KaonFromPhi->SetLineColor(kBlue);
    hpyShift_KaonFromPhi->SetMarkerColor(kBlue);
    hpyShift_KaonFromPhi->Draw("PE SAME");
    hpzShift_KaonFromPhi->SetMarkerStyle(22);
    hpzShift_KaonFromPhi->SetLineColor(kGreen + 2);
    hpzShift_KaonFromPhi->SetMarkerColor(kGreen + 2);
    hpzShift_KaonFromPhi->Draw("PE SAME");

    legShiftAllKaons->Draw();

    TCanvas *cShiftPhi = new TCanvas("cShiftPhi", "cShiftPhi", 720, 720);
    SetCanvasStyle(cShiftPhi, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hpxShift_Phi);
    SetHistoQA(hpyShift_Phi);
    SetHistoQA(hpzShift_Phi);
    hpxShift_Phi->SetMarkerStyle(20);
    hpxShift_Phi->SetLineColor(kRed);
    hpxShift_Phi->SetMarkerColor(kRed);
    hpxShift_Phi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hpxShift_Phi->GetYaxis()->SetTitle("<#Delta p_{x}> (GeV/c)");
    hpxShift_Phi->Draw("PE");
    hpyShift_Phi->SetMarkerStyle(21);
    hpyShift_Phi->SetLineColor(kBlue);
    hpyShift_Phi->SetMarkerColor(kBlue);
    hpyShift_Phi->Draw("PE SAME");
    hpzShift_Phi->SetMarkerStyle(22);
    hpzShift_Phi->SetLineColor(kGreen + 2);
    hpzShift_Phi->SetMarkerColor(kGreen + 2);
    hpzShift_Phi->Draw("PE SAME");

    legShiftAllKaons->Draw();
}