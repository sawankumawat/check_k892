#include <iostream>
#include <iomanip>
#include "../src/style.h"

void ParticleRatio()
{
    string KstarPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/";
    TFile *fKstar = new TFile((KstarPath + "Results.root").c_str(), "read");
    if (fKstar->IsZombie())
    {
        cout << "Error: Kstar file not found" << endl;
        return;
    }

    TFile *fPion = new TFile("PiKp_Run3_Results/Sawan/Pi_results.root", "read");
    TFile *fProton = new TFile("PiKp_Run3_Results/Sawan/Pr_results.root", "read");
    TFile *fKaon = new TFile("PiKp_Run3_Results/Sawan/Ka_results.root", "read");

    TFile *fChKstar = new TFile("ChargedKstar_Results/Sawan/ResultsChargedKstar.root", "read");
    TFile *fXiStar = new TFile("XiStar_Results/Sawan/ResultsXiStar.root", "read");
    TFile *fKshort = new TFile("K0s_Run3_Results/Sawan/ResultsK0s.root", "read");

    if (fPion->IsZombie() || fProton->IsZombie() || fKaon->IsZombie() || fChKstar->IsZombie() || fXiStar->IsZombie() || fKshort->IsZombie())
    {
        cout << "Error: Particle files not found" << endl;
        return;
    }

    TGraphErrors *gMPtKstar = (TGraphErrors *)fKstar->Get("gMeanpTRun3_sys");
    TGraphErrors *gMYieldKstar = (TGraphErrors *)fKstar->Get("gMeanYieldRun3_sys");
    if (gMPtKstar == nullptr || gMYieldKstar == nullptr)
    {
        cout << "Error: Kstar graphs not found" << endl;
        return;
    }
    TGraphErrors *gMPtPion = (TGraphErrors *)fPion->Get("gMeanpTRun3");
    TGraphErrors *gMPtProton = (TGraphErrors *)fProton->Get("gMeanpTRun3");
    TGraphErrors *gMPtKaon = (TGraphErrors *)fKaon->Get("gMeanpTRun3");
    TGraphErrors *gMPtChKstar = (TGraphErrors *)fChKstar->Get("gMeanpTRun3");
    TGraphErrors *gMPtXiStar = (TGraphErrors *)fXiStar->Get("gMeanpTRun3");
    TGraphErrors *gMPtKshort = (TGraphErrors *)fKshort->Get("gMeanpTRun3");

    TGraphErrors *gMYieldPion = (TGraphErrors *)fPion->Get("gMeanYieldRun3");
    TGraphErrors *gMYieldProton = (TGraphErrors *)fProton->Get("gMeanYieldRun3");
    TGraphErrors *gMYieldKaon = (TGraphErrors *)fKaon->Get("gMeanYieldRun3");
    TGraphErrors *gMYieldChKstar = (TGraphErrors *)fChKstar->Get("gMeanYieldRun3");
    TGraphErrors *gMYieldXiStar = (TGraphErrors *)fXiStar->Get("gMeanYieldRun3");
    TGraphErrors *gMYieldKshort = (TGraphErrors *)fKshort->Get("gMeanYieldRun3");

    if (gMPtPion == nullptr || gMPtProton == nullptr || gMPtKaon == nullptr || gMPtChKstar == nullptr || gMPtXiStar == nullptr || gMPtKshort == nullptr || gMYieldPion == nullptr || gMYieldProton == nullptr || gMYieldKaon == nullptr || gMYieldChKstar == nullptr || gMYieldXiStar == nullptr || gMYieldKshort == nullptr)
    {
        cout << "Error: Particle graphs not found" << endl;
        return;
    }

    //======================================================
    //    ===========EPOS local model files===========
    //======================================================
    TFile *fEPOS_UrQMDON = new TFile("ModelRootFiles/EPOS_final.root", "read");
    if (fEPOS_UrQMDON->IsZombie())
    {
        cout << "Error: EPOS file not found" << endl;
        return;
    }
    TGraphErrors *gMYieldKstarEPOS_UrQMDON = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON/kstar_vs_mult_urqmdON");
    TGraphErrors *hMYieldPhiEPOS_UrQMDON = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON/phi_vs_mult_urqmdON");
    TGraphErrors *gMYieldKshortEPOS_UrQMDON = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON/kshort_vs_mult_urqmdON");
    TGraphErrors *gMYieldProtonEPOS_UrQMDON = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON/proton_vs_mult_urqmdON");
    TGraphErrors *gMYieldPionEPOS_UrQMDON = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON/pion_vs_mult_urqmdON");
    TGraphErrors *gMYieldKaonEPOS_UrQMDON = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON/kaon_vs_mult_urqmdON");
    TGraphErrors *gMYieldXiStarEPOS_UrQMDON = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON/xistar_vs_mult_urqmdON");
    TGraphErrors *gMYieldChargedKstarEPOS_UrQMDON = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON/chargedkstar_vs_mult_urqmdON");

    TGraphErrors *gMYieldKstarEPOS_WithoutRescattering = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON_NoRescattering/kstar_vs_mult_urqmdON_noRescattering");
    TGraphErrors *gMYieldKstarEPOS_Rescattering = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON_Rescattering/kstar_vs_mult_urqmdON_rescattering");
    TGraphErrors *gMYieldKstarEPOS_NoUrQMD = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_OFF/kstar_vs_mult_urqmdOFF");

    TGraphErrors *gMeanPtKstarEPOS_UrQMDON = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON/meanpt_kstar_vs_mult_urqmdON");
    TGraphErrors *gMeanPtKstarEPOS_WithoutRescattering = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON_NoRescattering/meanpt_kstar_vs_mult_urqmdON_noRescattering");
    TGraphErrors *gMeanPtKstarEPOS_Rescattering = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_ON_Rescattering/meanpt_kstar_vs_mult_urqmdON_rescattering");
    TGraphErrors *gMeanPtKstarEPOS_NoUrQMD = (TGraphErrors *)fEPOS_UrQMDON->Get("UrQMD_OFF/meanpt_kstar_vs_mult_urqmdOFF");

    if (gMYieldKstarEPOS_UrQMDON == nullptr || hMYieldPhiEPOS_UrQMDON == nullptr || gMYieldKshortEPOS_UrQMDON == nullptr || gMYieldProtonEPOS_UrQMDON == nullptr || gMYieldPionEPOS_UrQMDON == nullptr || gMYieldKaonEPOS_UrQMDON == nullptr || gMYieldXiStarEPOS_UrQMDON == nullptr || gMYieldChargedKstarEPOS_UrQMDON == nullptr)
    {
        cout << "Error: EPOS graphs not found" << endl;
        return;
    }

    //======================================================
    //    ================Plots=====================
    //======================================================

    TCanvas *cdNdyKstar = new TCanvas("cdNdyKstar", "cdNdyKstar", 720, 720);
    SetCanvasStyle(cdNdyKstar, 0.15, 0.03, 0.03, 0.15);
    SetGraphErrorStyle(gMYieldKstar);
    SetGraphErrorStyle(gMYieldKstarEPOS_UrQMDON);
    SetGraphErrorStyle(gMYieldKstarEPOS_WithoutRescattering);
    SetGraphErrorStyle(gMYieldKstarEPOS_Rescattering);
    SetGraphErrorStyle(gMYieldKstarEPOS_NoUrQMD);
    gMYieldKstar->SetTitle(0);
    gMYieldKstar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMYieldKstar->GetYaxis()->SetTitle("dN/dy");
    gMYieldKstar->GetXaxis()->SetLimits(0, 27);
    gMYieldKstar->GetYaxis()->SetRangeUser(0, 0.8);
    gMYieldKstar->Draw("AP");
    gMYieldKstarEPOS_UrQMDON->Draw("l same");
    gMYieldKstarEPOS_WithoutRescattering->SetLineColor(kBlue);
    gMYieldKstarEPOS_WithoutRescattering->SetLineStyle(2);
    gMYieldKstarEPOS_WithoutRescattering->Draw("l same");
    gMYieldKstarEPOS_Rescattering->SetLineColor(kGreen + 2);
    gMYieldKstarEPOS_Rescattering->SetLineStyle(4);
    gMYieldKstarEPOS_Rescattering->Draw("l same");
    gMYieldKstarEPOS_NoUrQMD->SetLineColor(kMagenta);
    gMYieldKstarEPOS_NoUrQMD->SetLineStyle(3);
    gMYieldKstarEPOS_NoUrQMD->Draw("l same");

    TLegend *legend = new TLegend(0.2, 0.72, 0.5, 0.92);
    SetLegendStyle(legend);
    legend->SetTextSize(0.027);
    legend->SetHeader("#frac{K^{*0} + #bar{K}^{*0}}{2}", "C");
    legend->AddEntry(gMYieldKstar, "Data", "P");
    legend->AddEntry(gMYieldKstarEPOS_UrQMDON, "Status = 9", "L");
    legend->AddEntry(gMYieldKstarEPOS_WithoutRescattering, "Status = 9 && ity = 80", "L");
    legend->AddEntry(gMYieldKstarEPOS_Rescattering, "Status = 9 && ity = 81", "L");
    legend->AddEntry(gMYieldKstarEPOS_NoUrQMD, "Status = 7", "L");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    cdNdyKstar->SaveAs("Plots/MeanYield_Kstar_EPOS_UrQMDON.png");

    TCanvas *cMeanPtKstar = new TCanvas("cMeanPtKstar", "cMeanPtKstar", 720, 720);
    SetCanvasStyle(cMeanPtKstar, 0.15, 0.03, 0.03, 0.15);
    SetGraphErrorStyle(gMPtKstar);
    SetGraphErrorStyle(gMeanPtKstarEPOS_UrQMDON);
    SetGraphErrorStyle(gMeanPtKstarEPOS_WithoutRescattering);
    SetGraphErrorStyle(gMeanPtKstarEPOS_Rescattering);
    SetGraphErrorStyle(gMeanPtKstarEPOS_NoUrQMD);
    gMPtKstar->SetTitle(0);
    gMPtKstar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gMPtKstar->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    gMPtKstar->GetXaxis()->SetLimits(0, 27);
    gMPtKstar->GetYaxis()->SetRangeUser(0.0, 1.89);
    gMPtKstar->Draw("AP");
    gMeanPtKstarEPOS_UrQMDON->Draw("l same");
    gMeanPtKstarEPOS_WithoutRescattering->SetLineColor(kBlue);
    gMeanPtKstarEPOS_WithoutRescattering->SetLineStyle(2);
    gMeanPtKstarEPOS_WithoutRescattering->Draw("l same");
    gMeanPtKstarEPOS_Rescattering->SetLineColor(kGreen + 2);
    gMeanPtKstarEPOS_Rescattering->SetLineStyle(4);
    gMeanPtKstarEPOS_Rescattering->Draw("l same");
    gMeanPtKstarEPOS_NoUrQMD->SetLineColor(kMagenta);
    gMeanPtKstarEPOS_NoUrQMD->SetLineStyle(3);
    gMeanPtKstarEPOS_NoUrQMD->Draw("l same");
    legend->Draw();
    cMeanPtKstar->SaveAs("Plots/MeanPt_Kstar_EPOS_UrQMDON.png");

    //====================================================
    //   ================K*/K Ratio===================
    //====================================================

    TCanvas *cRatioKstarKaon = new TCanvas("cRatioKstarKaon", "cRatioKstarKaon", 720, 720);
    SetCanvasStyle(cRatioKstarKaon, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarKaon = new TGraphErrors(gMYieldKstar->GetN());
    gRatioKstarKaon->SetTitle(0);
    gRatioKstarKaon->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarKaon->GetYaxis()->SetTitle("dN/dy");
    gRatioKstarKaon->GetXaxis()->SetLimits(0, 27);
    if (gRatioKstarKaon->GetN() != gMYieldKaon->GetN())
    {
        cout << "Error: Kstar and Kaon graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarKaon->GetN(); ++i)
    {
        double xKaon, yKaon, xKstar, yKstar;
        gMYieldKaon->GetPoint(i, xKaon, yKaon);    // It is average of K+ and K- yields
        gMYieldKstar->GetPoint(i, xKstar, yKstar); // It is average of K*0 and K*0bar yields
        double yRatio = (yKaon != 0) ? yKstar / (yKaon) : 0;
        gRatioKstarKaon->SetPoint(i, xKstar, yRatio);
        double yKaonErr = gMYieldKaon->GetErrorY(i);
        double yKstarErr = gMYieldKstar->GetErrorY(i);
        double yRatioErr = (yKaon != 0) ? sqrt(pow(yKstarErr / yKaon, 2) + pow(yKstar * yKaonErr / (yKaon * yKaon), 2)) : 0;
        gRatioKstarKaon->SetPointError(i, 0, yRatioErr);
        // cout << "<dN_{ch}/d#eta> " << xKaon << ": Kstar error / value = " << (yKstarErr / yKstar) * 100 << "%, kaon error / value = " << (yKaonErr / yKaon) * 100 << "%, ratio error / value = " << (yRatioErr / yRatio) * 100 << "%" << endl;
    }
    SetGraphErrorStyle(gRatioKstarKaon);
    gRatioKstarKaon->GetYaxis()->SetRangeUser(0.2, 0.52);
    gRatioKstarKaon->GetXaxis()->SetLimits(0, 27);
    gRatioKstarKaon->Draw("APE");

    TGraphErrors *gRatioKstarKaonEPOS_UrQMDON = new TGraphErrors(gMYieldKstarEPOS_UrQMDON->GetN());
    TGraphErrors *gRatioKstarKaonEPOS_WithoutRescattering = new TGraphErrors(gMYieldKstarEPOS_WithoutRescattering->GetN());
    TGraphErrors *gRatioKstarKaonEPOS_Rescattering = new TGraphErrors(gMYieldKstarEPOS_Rescattering->GetN());
    TGraphErrors *gRatioKstarKaonEPOS_NoUrQMD = new TGraphErrors(gMYieldKstarEPOS_NoUrQMD->GetN());
    if (gMYieldKstarEPOS_UrQMDON->GetN() != gMYieldKaonEPOS_UrQMDON->GetN())
    {
        cout << "Error: Kstar and Kaon EPOS graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarKaonEPOS_UrQMDON->GetN(); ++i)
    {
        double xKaon, yKaon1, xKstar, yKstar1, yKstar2, yKstar3, yKstar4;
        gMYieldKaonEPOS_UrQMDON->GetPoint(i, xKaon, yKaon1);
        gMYieldKstarEPOS_UrQMDON->GetPoint(i, xKstar, yKstar1);
        double yRatio1 = (yKaon1 != 0) ? yKstar1 / (yKaon1) : 0;
        gRatioKstarKaonEPOS_UrQMDON->SetPoint(i, xKstar, yRatio1);
        gRatioKstarKaonEPOS_UrQMDON->SetPointError(i, 0, 0);
        // cout << "EPOS UrQMD ON - <dN_{ch}/d#eta> " << xKaon << ": Kstar yield = " << yKstar1 << ", Kaon yield = " << yKaon1 << endl;

        gMYieldKstarEPOS_WithoutRescattering->GetPoint(i, xKstar, yKstar2);
        double yRatio2 = (yKaon1 != 0) ? yKstar2 / (yKaon1) : 0;
        gRatioKstarKaonEPOS_WithoutRescattering->SetPoint(i, xKstar, yRatio2);
        gRatioKstarKaonEPOS_WithoutRescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_Rescattering->GetPoint(i, xKstar, yKstar3);
        double yRatio3 = (yKaon1 != 0) ? yKstar3 / (yKaon1) : 0;
        gRatioKstarKaonEPOS_Rescattering->SetPoint(i, xKstar, yRatio3);
        gRatioKstarKaonEPOS_Rescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_NoUrQMD->GetPoint(i, xKstar, yKstar4);
        double yRatio4 = (yKaon1 != 0) ? yKstar4 / (yKaon1) : 0;
        gRatioKstarKaonEPOS_NoUrQMD->SetPoint(i, xKstar, yRatio4);
        gRatioKstarKaonEPOS_NoUrQMD->SetPointError(i, 0, 0);
    }
    SetGraphErrorStyle(gRatioKstarKaonEPOS_UrQMDON);
    gRatioKstarKaonEPOS_UrQMDON->SetLineColor(kRed);
    gRatioKstarKaonEPOS_UrQMDON->SetLineStyle(2);
    gRatioKstarKaonEPOS_UrQMDON->Draw("l same");
    SetGraphErrorStyle(gRatioKstarKaonEPOS_WithoutRescattering);
    gRatioKstarKaonEPOS_WithoutRescattering->SetLineColor(kBlue);
    gRatioKstarKaonEPOS_WithoutRescattering->SetLineStyle(2);
    gRatioKstarKaonEPOS_WithoutRescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarKaonEPOS_Rescattering);
    gRatioKstarKaonEPOS_Rescattering->SetLineColor(kGreen + 2);
    gRatioKstarKaonEPOS_Rescattering->SetLineStyle(4);
    gRatioKstarKaonEPOS_Rescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarKaonEPOS_NoUrQMD);
    gRatioKstarKaonEPOS_NoUrQMD->SetLineColor(kMagenta);
    gRatioKstarKaonEPOS_NoUrQMD->SetLineStyle(5);
    gRatioKstarKaonEPOS_NoUrQMD->Draw("l same");
    TLegend *legendRatio = new TLegend(0.2, 0.72, 0.5, 0.92);
    SetLegendStyle(legendRatio);
    legendRatio->SetTextSize(0.027);
    legendRatio->AddEntry(gRatioKstarKaon, "#frac{K^{*0} + #bar{K}^{*0}}{2K}", "p");
    legendRatio->AddEntry(gRatioKstarKaonEPOS_UrQMDON, "EPOS UrQMD ON", "l");
    legendRatio->AddEntry(gRatioKstarKaonEPOS_WithoutRescattering, "EPOS NoRescattering", "l");
    legendRatio->AddEntry(gRatioKstarKaonEPOS_Rescattering, "EPOS Rescattering", "l");
    legendRatio->AddEntry(gRatioKstarKaonEPOS_NoUrQMD, "EPOS UrQMD OFF", "l");
    legendRatio->Draw();
    // cRatioKstarKaon->SaveAs("Plots/Ratio_KstarKaon_Run3.png");

    //====================================================
    //   ================K*/Pi Ratio===================
    //====================================================

    TCanvas *cRatioKstarPion = new TCanvas("cRatioKstarPion", "cRatioKstarPion", 720, 720);
    SetCanvasStyle(cRatioKstarPion, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarPion = new TGraphErrors(gMYieldKstar->GetN());
    gRatioKstarPion->SetTitle(0);
    gRatioKstarPion->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarPion->GetYaxis()->SetTitle("dN/dy");
    gRatioKstarPion->GetXaxis()->SetLimits(0, 27);
    if (gRatioKstarPion->GetN() != gMYieldPion->GetN())
    {
        cout << "Error: Kstar and Pion graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarPion->GetN(); ++i)
    {
        double xPion, yPion, xKstar, yKstar;
        gMYieldPion->GetPoint(i, xPion, yPion);    // It is average of Pi+ and Pi- yields
        gMYieldKstar->GetPoint(i, xKstar, yKstar); // It is average of K*0 and K*0bar yields
        double yRatio = (yPion != 0) ? yKstar / (yPion) : 0;
        gRatioKstarPion->SetPoint(i, xKstar, yRatio);
        double yPionErr = gMYieldPion->GetErrorY(i);
        double yKstarErr = gMYieldKstar->GetErrorY(i);
        double yRatioErr = (yPion != 0) ? sqrt(pow(yKstarErr / yPion, 2) + pow(yKstar * yPionErr / (yPion * yPion), 2)) : 0;
        gRatioKstarPion->SetPointError(i, 0, yRatioErr);
    }
    SetGraphErrorStyle(gRatioKstarPion);
    gRatioKstarPion->GetYaxis()->SetRangeUser(0.02, 0.12);
    gRatioKstarPion->GetXaxis()->SetLimits(0, 27);
    gRatioKstarPion->Draw("APE");

    TGraphErrors *gRatioKstarPionEPOS_UrQMDON = new TGraphErrors(gMYieldKstarEPOS_UrQMDON->GetN());
    TGraphErrors *gRatioKstarPionEPOS_WithoutRescattering = new TGraphErrors(gMYieldKstarEPOS_WithoutRescattering->GetN());
    TGraphErrors *gRatioKstarPionEPOS_Rescattering = new TGraphErrors(gMYieldKstarEPOS_Rescattering->GetN());
    TGraphErrors *gRatioKstarPionEPOS_NoUrQMD = new TGraphErrors(gMYieldKstarEPOS_NoUrQMD->GetN());
    if (gMYieldKstarEPOS_UrQMDON->GetN() != gMYieldPionEPOS_UrQMDON->GetN())
    {
        cout << "Error: Kstar and Pion EPOS graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarPionEPOS_UrQMDON->GetN(); ++i)
    {
        double xPion, yPion, xKstar, yKstar1, yKstar2, yKstar3, yKstar4;
        gMYieldPionEPOS_UrQMDON->GetPoint(i, xPion, yPion);
        gMYieldKstarEPOS_UrQMDON->GetPoint(i, xKstar, yKstar1);
        double yRatio = (yPion != 0) ? yKstar1 / (yPion) : 0;
        gRatioKstarPionEPOS_UrQMDON->SetPoint(i, xKstar, yRatio);
        gRatioKstarPionEPOS_UrQMDON->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_WithoutRescattering->GetPoint(i, xKstar, yKstar2);
        double yRatio2 = (yPion != 0) ? yKstar2 / (yPion) : 0;
        gRatioKstarPionEPOS_WithoutRescattering->SetPoint(i, xKstar, yRatio2);
        gRatioKstarPionEPOS_WithoutRescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_Rescattering->GetPoint(i, xKstar, yKstar3);
        double yRatio3 = (yPion != 0) ? yKstar3 / (yPion) : 0;
        gRatioKstarPionEPOS_Rescattering->SetPoint(i, xKstar, yRatio3);
        gRatioKstarPionEPOS_Rescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_NoUrQMD->GetPoint(i, xKstar, yKstar4);
        double yRatio4 = (yPion != 0) ? yKstar4 / (yPion) : 0;
        gRatioKstarPionEPOS_NoUrQMD->SetPoint(i, xKstar, yRatio4);
        gRatioKstarPionEPOS_NoUrQMD->SetPointError(i, 0, 0);
    }
    SetGraphErrorStyle(gRatioKstarPionEPOS_UrQMDON);
    gRatioKstarPionEPOS_UrQMDON->SetLineColor(kRed);
    gRatioKstarPionEPOS_UrQMDON->SetLineStyle(2);
    gRatioKstarPionEPOS_UrQMDON->Draw("l same");
    SetGraphErrorStyle(gRatioKstarPionEPOS_WithoutRescattering);
    gRatioKstarPionEPOS_WithoutRescattering->SetLineColor(kBlue);
    gRatioKstarPionEPOS_WithoutRescattering->SetLineStyle(2);
    gRatioKstarPionEPOS_WithoutRescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarPionEPOS_Rescattering);
    gRatioKstarPionEPOS_Rescattering->SetLineColor(kGreen + 2);
    gRatioKstarPionEPOS_Rescattering->SetLineStyle(4);
    gRatioKstarPionEPOS_Rescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarPionEPOS_NoUrQMD);
    gRatioKstarPionEPOS_NoUrQMD->SetLineColor(kMagenta);
    gRatioKstarPionEPOS_NoUrQMD->SetLineStyle(5);
    gRatioKstarPionEPOS_NoUrQMD->Draw("l same");
    TLegend *legendRatioPion = new TLegend(0.2, 0.72, 0.5, 0.92);
    SetLegendStyle(legendRatioPion);
    legendRatioPion->SetTextSize(0.027);
    legendRatioPion->AddEntry(gRatioKstarPion, "#frac{K^{*0} + #bar{K}^{*0}}{#pi^{+} + #pi^{-}}", "p");
    legendRatioPion->AddEntry(gRatioKstarPionEPOS_UrQMDON, "EPOS UrQMD ON", "l");
    legendRatioPion->AddEntry(gRatioKstarPionEPOS_WithoutRescattering, "EPOS NoRescattering", "l");
    legendRatioPion->AddEntry(gRatioKstarPionEPOS_Rescattering, "EPOS Rescattering", "l");
    legendRatioPion->AddEntry(gRatioKstarPionEPOS_NoUrQMD, "EPOS UrQMD OFF", "l");
    legendRatioPion->Draw();
    // cRatioKstarPion->SaveAs("Plots/Ratio_KstarPion_Run3.png");

    //====================================================
    //   ================K*0/K*+- Ratio===================
    //====================================================
    TCanvas *cRatioKstarChargedKstar = new TCanvas("cRatioKstarChargedKstar", "cRatioKstarChargedKstar", 720, 720);
    SetCanvasStyle(cRatioKstarChargedKstar, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarChargedKstar = new TGraphErrors(gMYieldKstar->GetN());
    gRatioKstarChargedKstar->SetTitle(0);
    gRatioKstarChargedKstar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarChargedKstar->GetYaxis()->SetTitle("dN/dy");
    gRatioKstarChargedKstar->GetXaxis()->SetLimits(0, 27);
    if (gRatioKstarChargedKstar->GetN() != gMYieldChKstar->GetN())
    {
        cout << "Error: Kstar and Charged Kstar graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarChargedKstar->GetN(); ++i)
    {
        double xChKstar, yChKstar, xKstar, yKstar;
        gMYieldChKstar->GetPoint(i, xChKstar, yChKstar); // It is average of K*+ and K*- yields
        gMYieldKstar->GetPoint(i, xKstar, yKstar);       // It is average of K*0 and K*0bar yields
        double yRatio = (yChKstar != 0) ? yKstar / (yChKstar) : 0;
        gRatioKstarChargedKstar->SetPoint(i, xKstar, yRatio);
        double yChKstarErr = gMYieldChKstar->GetErrorY(i);
        double yKstarErr = gMYieldKstar->GetErrorY(i);
        double yRatioErr = (yChKstar != 0) ? sqrt(pow(yKstarErr / yChKstar, 2) + pow(yKstar * yChKstarErr / (yChKstar * yChKstar), 2)) : 0;
        gRatioKstarChargedKstar->SetPointError(i, 0, yRatioErr);
    }
    SetGraphErrorStyle(gRatioKstarChargedKstar);
    gRatioKstarChargedKstar->GetYaxis()->SetRangeUser(0.12, 1.49);
    gRatioKstarChargedKstar->GetXaxis()->SetLimits(0, 27);
    gRatioKstarChargedKstar->Draw("APE");

    TGraphErrors *gRatioKstarChargedKstarEPOS_UrQMDON = new TGraphErrors(gMYieldKstarEPOS_UrQMDON->GetN());
    TGraphErrors *gRatioKstarChargedKstarEPOS_WithoutRescattering = new TGraphErrors(gMYieldKstarEPOS_WithoutRescattering->GetN());
    TGraphErrors *gRatioKstarChargedKstarEPOS_Rescattering = new TGraphErrors(gMYieldKstarEPOS_Rescattering->GetN());
    TGraphErrors *gRatioKstarChargedKstarEPOS_NoUrQMD = new TGraphErrors(gMYieldKstarEPOS_NoUrQMD->GetN());
    if (gMYieldKstarEPOS_UrQMDON->GetN() != gMYieldChargedKstarEPOS_UrQMDON->GetN())
    {
        cout << "Error: Kstar and Charged Kstar EPOS graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarChargedKstarEPOS_UrQMDON->GetN(); ++i)
    {
        double xChKstar, yChKstar, xKstar, yKstar1, yKstar2, yKstar3, yKstar4;
        gMYieldChargedKstarEPOS_UrQMDON->GetPoint(i, xChKstar, yChKstar);
        gMYieldKstarEPOS_UrQMDON->GetPoint(i, xKstar, yKstar1);
        double yRatio = (yChKstar != 0) ? yKstar1 / (yChKstar) : 0;
        gRatioKstarChargedKstarEPOS_UrQMDON->SetPoint(i, xKstar, yRatio);
        gRatioKstarChargedKstarEPOS_UrQMDON->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_WithoutRescattering->GetPoint(i, xKstar, yKstar2);
        double yRatio2 = (yChKstar != 0) ? yKstar2 / (yChKstar) : 0;
        gRatioKstarChargedKstarEPOS_WithoutRescattering->SetPoint(i, xKstar, yRatio2);
        gRatioKstarChargedKstarEPOS_WithoutRescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_Rescattering->GetPoint(i, xKstar, yKstar3);
        double yRatio3 = (yChKstar != 0) ? yKstar3 / (yChKstar) : 0;
        gRatioKstarChargedKstarEPOS_Rescattering->SetPoint(i, xKstar, yRatio3);
        gRatioKstarChargedKstarEPOS_Rescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_NoUrQMD->GetPoint(i, xKstar, yKstar4);
        double yRatio4 = (yChKstar != 0) ? yKstar4 / (yChKstar) : 0;
        gRatioKstarChargedKstarEPOS_NoUrQMD->SetPoint(i, xKstar, yRatio4);
        gRatioKstarChargedKstarEPOS_NoUrQMD->SetPointError(i, 0, 0);
    }
    SetGraphErrorStyle(gRatioKstarChargedKstarEPOS_UrQMDON);
    gRatioKstarChargedKstarEPOS_UrQMDON->SetLineColor(kRed);
    gRatioKstarChargedKstarEPOS_UrQMDON->SetLineStyle(2);
    gRatioKstarChargedKstarEPOS_UrQMDON->Draw("l same");
    SetGraphErrorStyle(gRatioKstarChargedKstarEPOS_WithoutRescattering);
    gRatioKstarChargedKstarEPOS_WithoutRescattering->SetLineColor(kBlue);
    gRatioKstarChargedKstarEPOS_WithoutRescattering->SetLineStyle(2);
    gRatioKstarChargedKstarEPOS_WithoutRescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarChargedKstarEPOS_Rescattering);
    gRatioKstarChargedKstarEPOS_Rescattering->SetLineColor(kGreen + 2);
    gRatioKstarChargedKstarEPOS_Rescattering->SetLineStyle(4);
    gRatioKstarChargedKstarEPOS_Rescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarChargedKstarEPOS_NoUrQMD);
    gRatioKstarChargedKstarEPOS_NoUrQMD->SetLineColor(kMagenta);
    gRatioKstarChargedKstarEPOS_NoUrQMD->SetLineStyle(4);
    gRatioKstarChargedKstarEPOS_NoUrQMD->Draw("l same");

    TLegend *legendRatioChargedKstar = new TLegend(0.2, 0.72, 0.5, 0.92);
    SetLegendStyle(legendRatioChargedKstar);
    legendRatioChargedKstar->SetTextSize(0.027);
    legendRatioChargedKstar->AddEntry(gRatioKstarChargedKstar, "#frac{K^{*0} + #bar{K}^{*0}}{K^{*+} + K^{*-}}", "p");
    legendRatioChargedKstar->AddEntry(gRatioKstarChargedKstarEPOS_UrQMDON, "EPOS UrQMD ON", "l");
    legendRatioChargedKstar->AddEntry(gRatioKstarChargedKstarEPOS_WithoutRescattering, "EPOS NoRescattering", "l");
    legendRatioChargedKstar->AddEntry(gRatioKstarChargedKstarEPOS_Rescattering, "EPOS Rescattering", "l");
    legendRatioChargedKstar->AddEntry(gRatioKstarChargedKstarEPOS_NoUrQMD, "EPOS UrQMD OFF", "l");
    legendRatioChargedKstar->Draw();
    // cRatioKstarChargedKstar->SaveAs("Plots/Ratio_KstarChargedKstar_Run3.png");

    //====================================================
    //   ================K*/Xi* Ratio===================
    //====================================================
    TCanvas *cRatioKstarXiStar = new TCanvas("cRatioKstarXiStar", "cRatioKstarXiStar", 720, 720);
    SetCanvasStyle(cRatioKstarXiStar, 0.15, 0.03, 0.03, 0.15);
    gPad->SetLogy();
    TGraphErrors *gRatioKstarXiStar = new TGraphErrors(gMYieldKstar->GetN());
    gRatioKstarXiStar->SetTitle(0);
    gRatioKstarXiStar->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarXiStar->GetYaxis()->SetTitle("dN/dy");
    gRatioKstarXiStar->GetXaxis()->SetLimits(0, 27);
    if (gRatioKstarXiStar->GetN() != gMYieldXiStar->GetN())
    {
        cout << "Error: Kstar and XiStar graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarXiStar->GetN(); ++i)
    {
        double xXiStar, yXiStar, xKstar, yKstar;
        gMYieldXiStar->GetPoint(i, xXiStar, yXiStar);
        gMYieldKstar->GetPoint(i, xKstar, yKstar);
        double yRatio = (yXiStar != 0) ? yKstar / (yXiStar) : 0;
        gRatioKstarXiStar->SetPoint(i, xKstar, yRatio);
        double yXiStarErr = gMYieldXiStar->GetErrorY(i);
        double yKstarErr = gMYieldKstar->GetErrorY(i);
        double yRatioErr = (yXiStar != 0) ? sqrt(pow(yKstarErr / yXiStar, 2) + pow(yKstar * yXiStarErr / (yXiStar * yXiStar), 2)) : 0;
        gRatioKstarXiStar->SetPointError(i, 0, yRatioErr);
    }
    SetGraphErrorStyle(gRatioKstarXiStar);
    gRatioKstarXiStar->GetYaxis()->SetRangeUser(0, 499);
    gRatioKstarXiStar->GetXaxis()->SetLimits(0, 27);
    gRatioKstarXiStar->Draw("APE");

    TGraphErrors *gRatioKstarXiStarEPOS_UrQMDON = new TGraphErrors(gMYieldKstarEPOS_UrQMDON->GetN());
    TGraphErrors *gRatioKstarXiStarEPOS_WithoutRescattering = new TGraphErrors(gMYieldKstarEPOS_WithoutRescattering->GetN());
    TGraphErrors *gRatioKstarXiStarEPOS_Rescattering = new TGraphErrors(gMYieldKstarEPOS_Rescattering->GetN());
    TGraphErrors *gRatioKstarXiStarEPOS_NoUrQMD = new TGraphErrors(gMYieldKstarEPOS_NoUrQMD->GetN());
    if (gMYieldKstarEPOS_UrQMDON->GetN() != gMYieldXiStarEPOS_UrQMDON->GetN())
    {
        cout << "Error: Kstar and XiStar EPOS graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarXiStarEPOS_UrQMDON->GetN(); ++i)
    {
        double xXiStar, yXiStar, xKstar, yKstar1, yKstar2, yKstar3, yKstar4;
        gMYieldXiStarEPOS_UrQMDON->GetPoint(i, xXiStar, yXiStar);
        gMYieldKstarEPOS_UrQMDON->GetPoint(i, xKstar, yKstar1);
        double yRatio = (yXiStar != 0) ? yKstar1 / (yXiStar) : 0;
        gRatioKstarXiStarEPOS_UrQMDON->SetPoint(i, xKstar, yRatio);
        gRatioKstarXiStarEPOS_UrQMDON->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_WithoutRescattering->GetPoint(i, xKstar, yKstar2);
        double yRatio2 = (yXiStar != 0) ? yKstar2 / (yXiStar) : 0;
        gRatioKstarXiStarEPOS_WithoutRescattering->SetPoint(i, xKstar, yRatio2);
        gRatioKstarXiStarEPOS_WithoutRescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_Rescattering->GetPoint(i, xKstar, yKstar3);
        double yRatio3 = (yXiStar != 0) ? yKstar3 / (yXiStar) : 0;
        gRatioKstarXiStarEPOS_Rescattering->SetPoint(i, xKstar, yRatio3);
        gRatioKstarXiStarEPOS_Rescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_NoUrQMD->GetPoint(i, xKstar, yKstar4);
        double yRatio4 = (yXiStar != 0) ? yKstar4 / (yXiStar) : 0;
        gRatioKstarXiStarEPOS_NoUrQMD->SetPoint(i, xKstar, yRatio4);
        gRatioKstarXiStarEPOS_NoUrQMD->SetPointError(i, 0, 0);
    }
    SetGraphErrorStyle(gRatioKstarXiStarEPOS_UrQMDON);
    gRatioKstarXiStarEPOS_UrQMDON->SetLineColor(kRed);
    gRatioKstarXiStarEPOS_UrQMDON->SetLineStyle(2);
    gRatioKstarXiStarEPOS_UrQMDON->Draw("l same");
    SetGraphErrorStyle(gRatioKstarXiStarEPOS_WithoutRescattering);
    gRatioKstarXiStarEPOS_WithoutRescattering->SetLineColor(kBlue);
    gRatioKstarXiStarEPOS_WithoutRescattering->SetLineStyle(2);
    gRatioKstarXiStarEPOS_WithoutRescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarXiStarEPOS_Rescattering);
    gRatioKstarXiStarEPOS_Rescattering->SetLineColor(kGreen + 2);
    gRatioKstarXiStarEPOS_Rescattering->SetLineStyle(4);
    gRatioKstarXiStarEPOS_Rescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarXiStarEPOS_NoUrQMD);
    gRatioKstarXiStarEPOS_NoUrQMD->SetLineColor(kMagenta);
    gRatioKstarXiStarEPOS_NoUrQMD->SetLineStyle(5);
    gRatioKstarXiStarEPOS_NoUrQMD->Draw("l same");
    TLegend *legendRatioXiStar = new TLegend(0.2, 0.72, 0.5, 0.92);
    SetLegendStyle(legendRatioXiStar);
    legendRatioXiStar->SetTextSize(0.027);
    legendRatioXiStar->AddEntry(gRatioKstarXiStar, "#frac{K^{*0} + #bar{K}^{*0}}{#Xi^{*0} + #bar{#Xi}^{*0}}", "p");
    legendRatioXiStar->AddEntry(gRatioKstarXiStarEPOS_UrQMDON, "EPOS UrQMD ON", "l");
    legendRatioXiStar->AddEntry(gRatioKstarXiStarEPOS_WithoutRescattering, "EPOS NoRescattering", "l");
    legendRatioXiStar->AddEntry(gRatioKstarXiStarEPOS_Rescattering, "EPOS Rescattering", "l");
    legendRatioXiStar->AddEntry(gRatioKstarXiStarEPOS_NoUrQMD, "EPOS UrQMD OFF", "l");
    legendRatioXiStar->Draw();
    // cRatioKstarXiStar->SaveAs("Plots/Ratio_KstarXiStar_Run3.png");

    //====================================================
    //   ================K*/Kshort Ratio===================
    //====================================================
    TCanvas *cRatioKstarKshort = new TCanvas("cRatioKstarKshort", "cRatioKstarKshort", 720, 720);
    SetCanvasStyle(cRatioKstarKshort, 0.15, 0.03, 0.03, 0.15);
    TGraphErrors *gRatioKstarKshort = new TGraphErrors(gMYieldKstar->GetN());
    gRatioKstarKshort->SetTitle(0);
    gRatioKstarKshort->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    gRatioKstarKshort->GetYaxis()->SetTitle("dN/dy");
    gRatioKstarKshort->GetXaxis()->SetLimits(0, 27);
    if (gRatioKstarKshort->GetN() != gMYieldKshort->GetN())
    {
        cout << "Error: Kstar and Kshort graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarKshort->GetN(); ++i)
    {
        double xKshort, yKshort, xKstar, yKstar;
        gMYieldKshort->GetPoint(i, xKshort, yKshort);
        gMYieldKstar->GetPoint(i, xKstar, yKstar);
        double yRatio = (yKshort != 0) ? yKstar / (yKshort) : 0;
        gRatioKstarKshort->SetPoint(i, xKstar, yRatio);
        double yKshortErr = gMYieldKshort->GetErrorY(i);
        double yKstarErr = gMYieldKstar->GetErrorY(i);
        double yRatioErr = (yKshort != 0) ? sqrt(pow(yKstarErr / yKshort, 2) + pow(yKstar * yKshortErr / (yKshort * yKshort), 2)) : 0;
        gRatioKstarKshort->SetPointError(i, 0, yRatioErr);
    }
    SetGraphErrorStyle(gRatioKstarKshort);
    gRatioKstarKshort->GetYaxis()->SetRangeUser(0.2, 0.52);
    gRatioKstarKshort->GetXaxis()->SetLimits(0, 27);
    gRatioKstarKshort->Draw("APE");

    TGraphErrors *gRatioKstarKshortEPOS_UrQMDON = new TGraphErrors(gMYieldKstarEPOS_UrQMDON->GetN());
    TGraphErrors *gRatioKstarKshortEPOS_WithoutRescattering = new TGraphErrors(gMYieldKstarEPOS_WithoutRescattering->GetN());
    TGraphErrors *gRatioKstarKshortEPOS_Rescattering = new TGraphErrors(gMYieldKstarEPOS_Rescattering->GetN());
    TGraphErrors *gRatioKstarKshortEPOS_NoUrQMD = new TGraphErrors(gMYieldKstarEPOS_NoUrQMD->GetN());
    if (gMYieldKstarEPOS_UrQMDON->GetN() != gMYieldKshortEPOS_UrQMDON->GetN())
    {
        cout << "Error: Kstar and Kshort EPOS graphs have different number of points" << endl;
        return;
    }
    for (int i = 0; i < gRatioKstarKshortEPOS_UrQMDON->GetN(); ++i)
    {
        double xKshort, yKshort, xKstar, yKstar1, yKstar2, yKstar3, yKstar4;
        gMYieldKshortEPOS_UrQMDON->GetPoint(i, xKshort, yKshort);
        gMYieldKstarEPOS_UrQMDON->GetPoint(i, xKstar, yKstar1);
        double yRatio = (yKshort != 0) ? yKstar1 / (yKshort) : 0;
        gRatioKstarKshortEPOS_UrQMDON->SetPoint(i, xKstar, yRatio);
        gRatioKstarKshortEPOS_UrQMDON->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_WithoutRescattering->GetPoint(i, xKstar, yKstar2);
        double yRatio2 = (yKshort != 0) ? yKstar2 / (yKshort) : 0;
        gRatioKstarKshortEPOS_WithoutRescattering->SetPoint(i, xKstar, yRatio2);
        gRatioKstarKshortEPOS_WithoutRescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_Rescattering->GetPoint(i, xKstar, yKstar3);
        double yRatio3 = (yKshort != 0) ? yKstar3 / (yKshort) : 0;
        gRatioKstarKshortEPOS_Rescattering->SetPoint(i, xKstar, yRatio3);
        gRatioKstarKshortEPOS_Rescattering->SetPointError(i, 0, 0);

        gMYieldKstarEPOS_NoUrQMD->GetPoint(i, xKstar, yKstar4);
        double yRatio4 = (yKshort != 0) ? yKstar4 / (yKshort) : 0;
        gRatioKstarKshortEPOS_NoUrQMD->SetPoint(i, xKstar, yRatio4);
        gRatioKstarKshortEPOS_NoUrQMD->SetPointError(i, 0, 0);
    }
    SetGraphErrorStyle(gRatioKstarKshortEPOS_UrQMDON);
    gRatioKstarKshortEPOS_UrQMDON->SetLineColor(kRed);
    gRatioKstarKshortEPOS_UrQMDON->SetLineStyle(2);
    gRatioKstarKshortEPOS_UrQMDON->Draw("l same");
    SetGraphErrorStyle(gRatioKstarKshortEPOS_WithoutRescattering);
    gRatioKstarKshortEPOS_WithoutRescattering->SetLineColor(kBlue);
    gRatioKstarKshortEPOS_WithoutRescattering->SetLineStyle(2);
    gRatioKstarKshortEPOS_WithoutRescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarKshortEPOS_Rescattering);
    gRatioKstarKshortEPOS_Rescattering->SetLineColor(kGreen + 2);
    gRatioKstarKshortEPOS_Rescattering->SetLineStyle(4);
    gRatioKstarKshortEPOS_Rescattering->Draw("l same");
    SetGraphErrorStyle(gRatioKstarKshortEPOS_NoUrQMD);
    gRatioKstarKshortEPOS_NoUrQMD->SetLineColor(kMagenta);
    gRatioKstarKshortEPOS_NoUrQMD->SetLineStyle(5);
    gRatioKstarKshortEPOS_NoUrQMD->Draw("l same");
    TLegend *legendRatioKshort = new TLegend(0.2, 0.72, 0.5, 0.92);
    SetLegendStyle(legendRatioKshort);
    legendRatioKshort->SetTextSize(0.027);
    legendRatioKshort->AddEntry(gRatioKstarKshort, "#frac{K^{*0} + #bar{K}^{*0}}{K_{S}^{0}}", "p");
    legendRatioKshort->AddEntry(gRatioKstarKshortEPOS_UrQMDON, "EPOS UrQMD ON", "l");
    legendRatioKshort->AddEntry(gRatioKstarKshortEPOS_WithoutRescattering, "EPOS NoRescattering", "l");
    legendRatioKshort->AddEntry(gRatioKstarKshortEPOS_Rescattering, "EPOS Rescattering", "l");
    legendRatioKshort->AddEntry(gRatioKstarKshortEPOS_NoUrQMD, "EPOS UrQMD OFF", "l");
    legendRatioKshort->Draw();
    // cRatioKstarKshort->SaveAs("Plots/Ratio_KstarKshort_Run3.png");
}