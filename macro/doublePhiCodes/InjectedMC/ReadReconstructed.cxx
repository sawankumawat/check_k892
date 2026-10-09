#include <iostream>
#include "../../src/style.h"
#include "../../src/fitfunc.h"

void ReadReconstructed()
{
    TFile *file = new TFile("/home/sawan/alice/practice/AnalysisResults.root", "READ");
    if (file->IsZombie())
    {
        std::cerr << "Error: Could not open the file." << std::endl;
        return;
    }

    gStyle->SetOptStat(1110);

    THnSparse *hT4sMassPtGen = (THnSparse *)file->Get("phi-phi-m-c/Generated/h2_InvMassPhiPhi");
    TH1D *hT4sMassPtGen_1D = hT4sMassPtGen->Projection(0);    // Projecting onto the first axis (mass)
    TH1D *hT4sMassPtGen_1D_Pt = hT4sMassPtGen->Projection(1); // Projecting onto the second axis (pt)
    TH1D *hT4sPhiGen = (TH1D *)file->Get("phi-phi-m-c/Generated/h1_PhiPhi_Phi");
    TH1D *hT4sEtaGen = (TH1D *)file->Get("phi-phi-m-c/Generated/h1_PhiPhi_Eta");

    TCanvas *cMassGen = new TCanvas("cMassGen", "Generated Mass Distribution", 720, 720);
    SetCanvasStyle(cMassGen, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hT4sMassPtGen_1D);
    gPad->SetLogy();
    hT4sMassPtGen_1D->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hT4sMassPtGen_1D->GetYaxis()->SetTitle("Counts");
    hT4sMassPtGen_1D->Draw();
    cMassGen->SaveAs("Generated_T4sMass.png");

    TCanvas *cPtGen = new TCanvas("cPtGen", "Generated Pt Distribution", 720, 720);
    SetCanvasStyle(cPtGen, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hT4sMassPtGen_1D_Pt);
    hT4sMassPtGen_1D_Pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hT4sMassPtGen_1D_Pt->GetYaxis()->SetTitle("Counts");
    hT4sMassPtGen_1D_Pt->SetMinimum(0);
    hT4sMassPtGen_1D_Pt->SetMaximum(hT4sMassPtGen_1D_Pt->GetMaximum() * 1.3); // Add some space above the highest bin
    hT4sMassPtGen_1D_Pt->Draw();
    cPtGen->SaveAs("Generated_T4sPt.png");

    TCanvas *cPhiGen = new TCanvas("cPhiGen", "Generated Phi Distribution", 720, 720);
    SetCanvasStyle(cPhiGen, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hT4sPhiGen);
    hT4sPhiGen->GetXaxis()->SetTitle("#it{#phi} (rad)");
    hT4sPhiGen->GetYaxis()->SetTitle("Counts");
    hT4sPhiGen->SetMinimum(0);
    hT4sPhiGen->SetMaximum(hT4sPhiGen->GetMaximum() * 1.3); // Add some space above the highest bin
    hT4sPhiGen->Draw();
    cPhiGen->SaveAs("Generated_T4sPhi.png");

    TCanvas *cEtaGen = new TCanvas("cEtaGen", "Generated Eta Distribution", 720, 720);
    SetCanvasStyle(cEtaGen, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hT4sEtaGen);
    hT4sEtaGen->GetXaxis()->SetTitle("#it{#eta}");
    hT4sEtaGen->GetYaxis()->SetTitle("Counts");
    hT4sEtaGen->SetMinimum(0);
    hT4sEtaGen->SetMaximum(hT4sEtaGen->GetMaximum() * 1.3); // Add some space above the highest bin
    hT4sEtaGen->Draw();
    cEtaGen->SaveAs("Generated_T4sEta.png");

    THnSparse *hPhiGenInjected = (THnSparse *)file->Get("phi-phi-m-c/Generated/h2_InvMassPhiInjected");
    THnSparse *hPhiGenDecayed = (THnSparse *)file->Get("phi-phi-m-c/Generated/h2_InvMassPhiFromT4s");
    THnSparse *hPhiGenMinBias = (THnSparse *)file->Get("phi-phi-m-c/Generated/h2_InvMassPhiMinBias");

    TH1D *hPhiGenInjectedMass = hPhiGenInjected->Projection(0);
    TH1D *hPhiGenDecayedMass = hPhiGenDecayed->Projection(0);
    TH1D *hPhiGenMinBiasMass = hPhiGenMinBias->Projection(0);
    TH1D *hPhiGenInjectedPt = hPhiGenInjected->Projection(1);
    TH1D *hPhiGenDecayedPt = hPhiGenDecayed->Projection(1);
    TH1D *hPhiGenMinBiasPt = hPhiGenMinBias->Projection(1);

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.05);

    // lambda function to set gpad margins
    auto setMargins = [](TPad *pad, double left, double right, double top, double bottom)
    {
        pad->SetLeftMargin(left);
        pad->SetRightMargin(right);
        pad->SetTopMargin(top);
        pad->SetBottomMargin(bottom);
    };

    TCanvas *cMassComparison = new TCanvas("cMassComparison", "Mass Comparison", 2000, 720);
    SetCanvasStyle(cMassComparison, 0.16, 0.06, 0.06, 0.13);
    cMassComparison->Divide(3, 1);
    SetHistoQA(hPhiGenInjectedMass);
    SetHistoQA(hPhiGenDecayedMass);
    SetHistoQA(hPhiGenMinBiasMass);
    cMassComparison->cd(1);
    gPad->SetLogy();
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    hPhiGenInjectedMass->GetXaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
    hPhiGenInjectedMass->GetYaxis()->SetTitle("Counts");
    hPhiGenInjectedMass->Draw();
    lat.DrawLatex(0.65, 0.7, "Injected #phi");
    cMassComparison->cd(2);
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    gPad->SetLogy();
    hPhiGenDecayedMass->GetXaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
    hPhiGenDecayedMass->GetYaxis()->SetTitle("Counts");
    hPhiGenDecayedMass->Draw();
    lat.DrawLatex(0.65, 0.7, "Decayed #phi");
    cMassComparison->cd(3);
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    gPad->SetLogy();
    hPhiGenMinBiasMass->GetXaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
    hPhiGenMinBiasMass->GetYaxis()->SetTitle("Counts");
    hPhiGenMinBiasMass->Draw();
    lat.DrawLatex(0.65, 0.7, "MinBias #phi");
    cMassComparison->SaveAs("Generated_PhiMassComparison.png");

    cMassComparison->SaveAs("Generated_PhiMassComparison.png");

    TCanvas *cPtComparison = new TCanvas("cPtComparison", "Pt Comparison", 2000, 720);
    SetCanvasStyle(cPtComparison, 0.16, 0.06, 0.01, 0.13);
    cPtComparison->Divide(3, 1);
    SetHistoQA(hPhiGenInjectedPt);
    SetHistoQA(hPhiGenDecayedPt);
    SetHistoQA(hPhiGenMinBiasPt);
    cPtComparison->cd(1);
    gPad->SetLogy();
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    hPhiGenInjectedPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiGenInjectedPt->GetYaxis()->SetTitle("Counts");
    // hPhiGenInjectedPt->SetMinimum(0);
    hPhiGenInjectedPt->Draw();
    lat.DrawLatex(0.65, 0.7, "Injected #phi");
    cPtComparison->cd(2);
    gPad->SetLogy();
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    hPhiGenDecayedPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiGenDecayedPt->GetYaxis()->SetTitle("Counts");
    // hPhiGenDecayedPt->SetMinimum(0);
    hPhiGenDecayedPt->Draw();
    lat.DrawLatex(0.65, 0.7, "Decayed #phi");
    cPtComparison->cd(3);
    gPad->SetLogy();
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    hPhiGenMinBiasPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiGenMinBiasPt->GetYaxis()->SetTitle("Counts");
    // hPhiGenMinBiasPt->SetMinimum(0);
    hPhiGenMinBiasPt->Draw();
    lat.DrawLatex(0.65, 0.7, "MinBias #phi");

    cPtComparison->SaveAs("Generated_PhiPtComparison.png");

    //=========================
    //    Reconstructed Histograms
    //=========================

    TH1D *hPhiRecMassInjected = (TH1D *)file->Get("phi-phi-m-c/Reconstructed/AllInjectedPhi/h1_Phi1_Mass");
    TH1D *hPhiRecpTInjected = (TH1D *)file->Get("phi-phi-m-c/Reconstructed/AllInjectedPhi/h1_Phi1_Pt");
    TH1D *hPhiRecMassDecayT4s = (TH1D *)file->Get("phi-phi-m-c/Reconstructed/PhiFromT4s/h1_Phi1_Mass");
    TH1D *hPhiRecpTDecayT4s = (TH1D *)file->Get("phi-phi-m-c/Reconstructed/PhiFromT4s/h1_Phi1_Pt");
    TH1D *hPhiRecMassMinBias = (TH1D *)file->Get("phi-phi-m-c/Reconstructed/MinimumBiasPhi/h1_Phi1_Mass");
    TH1D *hPhiRecpTMinBias = (TH1D *)file->Get("phi-phi-m-c/Reconstructed/MinimumBiasPhi/h1_Phi1_Pt");

    TCanvas *cMassRecComparison = new TCanvas("cMassRecComparison", "Reconstructed Mass Comparison", 2000, 720);
    SetCanvasStyle(cMassRecComparison, 0.16, 0.06, 0.06, 0.13);
    cMassRecComparison->Divide(3, 1);
    SetHistoQA(hPhiRecMassInjected);
    SetHistoQA(hPhiRecMassDecayT4s);
    SetHistoQA(hPhiRecMassMinBias);
    cMassRecComparison->cd(1);
    gPad->SetLogy();
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    hPhiRecMassInjected->GetXaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
    hPhiRecMassInjected->GetYaxis()->SetTitle("Counts");
    hPhiRecMassInjected->Draw();
    lat.DrawLatex(0.65, 0.7, "Injected #phi");
    cMassRecComparison->cd(2);
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    gPad->SetLogy();
    hPhiRecMassDecayT4s->GetXaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
    hPhiRecMassDecayT4s->GetYaxis()->SetTitle("Counts");
    hPhiRecMassDecayT4s->Draw();
    lat.DrawLatex(0.65, 0.7, "Decayed #phi");
    cMassRecComparison->cd(3);
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    gPad->SetLogy();
    hPhiRecMassMinBias->GetXaxis()->SetTitle("#it{M}_{#phi} (GeV/#it{c}^{2})");
    hPhiRecMassMinBias->GetYaxis()->SetTitle("Counts");
    hPhiRecMassMinBias->Draw();
    lat.DrawLatex(0.65, 0.7, "MinBias #phi");
    cMassRecComparison->SaveAs("Reconstructed_PhiMassComparison.png");

    TCanvas *cPtRecComparison = new TCanvas("cPtRecComparison", "Reconstructed Pt Comparison", 2000, 720);
    SetCanvasStyle(cPtRecComparison, 0.16, 0.06, 0.01, 0.13);
    cPtRecComparison->Divide(3, 1);
    SetHistoQA(hPhiRecpTInjected);
    SetHistoQA(hPhiRecpTDecayT4s);
    SetHistoQA(hPhiRecpTMinBias);
    cPtRecComparison->cd(1);
    gPad->SetLogy();
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    hPhiRecpTInjected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiRecpTInjected->GetYaxis()->SetTitle("Counts");
    // hPhiRecpTInjected->SetMinimum(0);
    hPhiRecpTInjected->Draw();
    lat.DrawLatex(0.65, 0.7, "Injected #phi");
    cPtRecComparison->cd(2);
    gPad->SetLogy();
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    hPhiRecpTDecayT4s->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiRecpTDecayT4s->GetYaxis()->SetTitle("Counts");
    // hPhiRecpTDecayT4s->SetMinimum(0);
    hPhiRecpTDecayT4s->Draw();
    lat.DrawLatex(0.65, 0.7, "Decayed #phi");
    cPtRecComparison->cd(3);
    gPad->SetLogy();
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    hPhiRecpTMinBias->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiRecpTMinBias->GetYaxis()->SetTitle("Counts");
    // hPhiRecpTMinBias->SetMinimum(0);
    hPhiRecpTMinBias->Draw();
    lat.DrawLatex(0.65, 0.7, "MinBias #phi");
    cPtRecComparison->SaveAs("Reconstructed_PhiPtComparison.png");

    THnSparse *hT4sMassRec_Injected = (THnSparse *)file->Get("phi-phi-m-c/Reconstructed/AllInjectedPhi/h2_PhiPhi_Injected");
    THnSparse *hT4sMassRec_Pure = (THnSparse *)file->Get("phi-phi-m-c/Reconstructed/PhiFromT4s/h2_PhiPhi_T4s");
    THnSparse *hT4sMassRec_MinBias = (THnSparse *)file->Get("phi-phi-m-c/Reconstructed/MinimumBiasPhi/h2_PhiPhi_MinBias");

    TH1D *hT4sMassRec_Injected_1D = hT4sMassRec_Injected->Projection(0);
    TH1D *hT4sMassRec_Pure_1D = hT4sMassRec_Pure->Projection(0);
    TH1D *hT4sMassRec_MinBias_1D = hT4sMassRec_MinBias->Projection(0);

    TH1D *hT4sPtRec_Injected_1D = hT4sMassRec_Injected->Projection(1);
    TH1D *hT4sPtRec_Pure_1D = hT4sMassRec_Pure->Projection(1);
    TH1D *hT4sPtRec_MinBias_1D = hT4sMassRec_MinBias->Projection(1);

    TCanvas *cMassRecComparison2 = new TCanvas("cMassRecComparison2", "Reconstructed Mass Comparison 2", 2000, 720);
    SetCanvasStyle(cMassRecComparison2, 0.16, 0.06, 0.06, 0.13);
    cMassRecComparison2->Divide(3, 1);
    SetHistoQA(hT4sMassRec_Injected_1D);
    SetHistoQA(hT4sMassRec_Pure_1D);
    SetHistoQA(hT4sMassRec_MinBias_1D);
    cMassRecComparison2->cd(1);
    gPad->SetLogy();
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    hT4sMassRec_Injected_1D->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hT4sMassRec_Injected_1D->GetYaxis()->SetTitle("Counts");
    hT4sMassRec_Injected_1D->SetMaximum(hT4sMassRec_Injected_1D->GetMaximum() * 13); // Add some space above the highest bin
    hT4sMassRec_Injected_1D->Draw();
    lat.DrawLatex(0.65, 0.7, "Injected #phi#phi");
    cMassRecComparison2->cd(2);
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    gPad->SetLogy();
    hT4sMassRec_Pure_1D->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hT4sMassRec_Pure_1D->GetYaxis()->SetTitle("Counts");
    hT4sMassRec_Pure_1D->Draw();
    lat.DrawLatex(0.65, 0.7, "Decayed #phi#phi");
    cMassRecComparison2->cd(3);
    setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    gPad->SetLogy();
    hT4sMassRec_MinBias_1D->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hT4sMassRec_MinBias_1D->GetYaxis()->SetTitle("Counts");
    hT4sMassRec_MinBias_1D->Draw();
    lat.DrawLatex(0.65, 0.7, "MinBias #phi#phi");
    cMassRecComparison2->SaveAs("Reconstructed_PhiPhiMassComparison2.png");

    // TCanvas *cPtRecComparison2 = new TCanvas("cPtRecComparison2", "Reconstructed Pt Comparison 2", 2000, 720);
    // SetCanvasStyle(cPtRecComparison2, 0.16, 0.06, 0.01, 0.13);
    // cPtRecComparison2->Divide(3, 1);
    // SetHistoQA(hT4sPtRec_Injected_1D);
    // SetHistoQA(hT4sPtRec_Pure_1D);
    // SetHistoQA(hT4sPtRec_MinBias_1D);
    // cPtRecComparison2->cd(1);
    // setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    // hT4sPtRec_Injected_1D->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // hT4sPtRec_Injected_1D->GetYaxis()->SetTitle("Counts");
    // hT4sPtRec_Injected_1D->SetMinimum(0);
    // hT4sPtRec_Injected_1D->SetMaximum(hT4sPtRec_Injected_1D->GetMaximum() * 1.5); // Add some space above the highest bin
    // hT4sPtRec_Injected_1D->Draw();
    // lat.DrawLatex(0.65, 0.7, "Injected #phi#phi");
    // cPtRecComparison2->cd(2);
    // setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    // hT4sPtRec_Pure_1D->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // hT4sPtRec_Pure_1D->GetYaxis()->SetTitle("Counts");
    // hT4sPtRec_Pure_1D->SetMinimum(0);
    // hT4sPtRec_Pure_1D->SetMaximum(hT4sPtRec_Pure_1D->GetMaximum() * 1.5); // Add some space above the highest bin
    // hT4sPtRec_Pure_1D->Draw();
    // lat.DrawLatex(0.65, 0.7, "Decayed #phi#phi");
    // cPtRecComparison2->cd(3);
    // setMargins(static_cast<TPad *>(gPad), 0.13, 0.03, 0.06, 0.12);
    // hT4sPtRec_MinBias_1D->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // hT4sPtRec_MinBias_1D->GetYaxis()->SetTitle("Counts");
    // hT4sPtRec_MinBias_1D->SetMinimum(0);
    // hT4sPtRec_MinBias_1D->Draw();
    // lat.DrawLatex(0.65, 0.7, "MinBias #phi#phi");
    // cPtRecComparison2->SaveAs("Reconstructed_PhiPhiPtComparison2.png");

    THnSparse *hPhiPhiRecAll = (THnSparse *)file->Get("phi-phi-m-c/Reconstructed/AllPhiPairs/h2_PhiPhi_All");
    TH1D *hPhiPhiRecAllMass = hPhiPhiRecAll->Projection(0);
    TH1D *hPhiPhiRecAllPt = hPhiPhiRecAll->Projection(1);

    TCanvas *cMassRecAll = new TCanvas("cMassRecAll", "Reconstructed Mass All", 720, 720);
    SetCanvasStyle(cMassRecAll, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hPhiPhiRecAllMass);
    gPad->SetLogy();
    hPhiPhiRecAllMass->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hPhiPhiRecAllMass->GetYaxis()->SetTitle("Counts");
    hPhiPhiRecAllMass->GetXaxis()->SetRangeUser(2.0, 3.5); // Set the x-axis range to focus on the region of interest
    hPhiPhiRecAllMass->Draw();
    cMassRecAll->SaveAs("Reconstructed_PhiPhiMassAll.png");

    TCanvas *cPtRecAll = new TCanvas("cPtRecAll", "Reconstructed Pt All", 720, 720);
    SetCanvasStyle(cPtRecAll, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hPhiPhiRecAllPt);
    hPhiPhiRecAllPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiPhiRecAllPt->GetYaxis()->SetTitle("Counts");
    hPhiPhiRecAllPt->SetMinimum(0);
    hPhiPhiRecAllPt->SetMaximum(hPhiPhiRecAllPt->GetMaximum() * 1.3); // Add some space above the highest bin
    hPhiPhiRecAllPt->Draw();
    cPtRecAll->SaveAs("Reconstructed_PhiPhiPtAll.png");

    //=========================
    // Reweighted distributions
    //=========================

    THnSparse *hPhiPhiInjectedGen = (THnSparse *)file->Get("phi-phi-m-c/Generated/h2_InvMassPhiPhiFromInjected");
    TH1D *hPhiPhiInjectedGenMass = hPhiPhiInjectedGen->Projection(0);
    TH1D *hPhiPhiInjectedGenPt = hPhiPhiInjectedGen->Projection(1);
    THnSparse *hPhiPtReweighted = (THnSparse *)file->Get("phi-phi-m-c/Generated/h2_InvMassPhiInjectedWeighted");
    TH1D *hPhiReweightedPt = hPhiPtReweighted->Projection(1);
    TH1D *hPhiPhiReweightedMass = (TH1D *)file->Get("phi-phi-m-c/Generated/h2_InvMassPhiPhiFromWeighted");
    TH1D *hPhiPhiReweightedPt = (TH1D *)file->Get("phi-phi-m-c/Generated/h2_InvMassPhiPhiPtFromWeighted");

    TCanvas *cMassInjectedGen = new TCanvas("cMassInjectedGen", "Generated Mass Injected", 720, 720);
    SetCanvasStyle(cMassInjectedGen, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hPhiPhiInjectedGenMass);
    gPad->SetLogy();
    hPhiPhiInjectedGenMass->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hPhiPhiInjectedGenMass->GetYaxis()->SetTitle("Counts");
    // hPhiPhiInjectedGenMass->GetXaxis()->SetRangeUser(2.0, 10.0); // Set the x-axis range to focus on the region of interest
    hPhiPhiInjectedGenMass->Draw();
    cMassInjectedGen->SaveAs("Generated_PhiPhiMassInjected.png");

    TCanvas *cPtInjectedGen = new TCanvas("cPtInjectedGen", "Generated Pt Injected", 720, 720);
    SetCanvasStyle(cPtInjectedGen, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hPhiPhiInjectedGenPt);
    gPad->SetLogy();
    hPhiPhiInjectedGenPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiPhiInjectedGenPt->GetYaxis()->SetTitle("Counts");
    // hPhiPhiInjectedGenPt->SetMinimum(0);
    hPhiPhiInjectedGenPt->SetMaximum(hPhiPhiInjectedGenPt->GetMaximum() * 1.3); // Add some space above the highest bin
    hPhiPhiInjectedGenPt->Draw();
    cPtInjectedGen->SaveAs("Generated_PhiPhiPtInjected.png");

    TCanvas *cPtReweighted = new TCanvas("cPtReweighted", "Reweighted Pt", 720, 720);
    SetCanvasStyle(cPtReweighted, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hPhiReweightedPt);
    gPad->SetLogy();
    hPhiReweightedPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiReweightedPt->GetYaxis()->SetTitle("Counts");
    // hPhiReweightedPt->SetMinimum(0);
    hPhiReweightedPt->SetMaximum(hPhiReweightedPt->GetMaximum() * 1.3); // Add some space above the highest bin
    hPhiReweightedPt->Draw();
    cPtReweighted->SaveAs("Generated_PhiPtReweighted.png");

    TCanvas *cPhiPhiReweightedMass = new TCanvas("cPhiPhiReweightedMass", "Reweighted Mass", 720, 720);
    SetCanvasStyle(cPhiPhiReweightedMass, 0.16, 0.06, 0.06, 0.13);
    hPhiPhiReweightedMass->Sumw2();
    hPhiPhiReweightedMass->Rebin(5);
    SetHistoQA(hPhiPhiReweightedMass);
    gPad->SetLogy();
    hPhiPhiReweightedMass->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hPhiPhiReweightedMass->GetYaxis()->SetTitle("Counts");
    hPhiPhiReweightedMass->GetXaxis()->SetRangeUser(2.0, 5.0); // Set the x-axis range to focus on the region of interest
    hPhiPhiReweightedMass->Draw();
    cPhiPhiReweightedMass->SaveAs("Generated_PhiPhiMassReweighted.png");

    TCanvas *cPhiPhiReweightedPt = new TCanvas("cPhiPhiReweightedPt", "Reweighted Pt", 720, 720);
    SetCanvasStyle(cPhiPhiReweightedPt, 0.16, 0.06, 0.06, 0.13);
    hPhiPhiReweightedPt->Sumw2();
    hPhiPhiReweightedPt->Rebin(5);
    SetHistoQA(hPhiPhiReweightedPt);
    gPad->SetLogy();
    hPhiPhiReweightedPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPhiPhiReweightedPt->GetYaxis()->SetTitle("Counts");
    // hPhiPhiReweightedPt->SetMinimum(0);
    hPhiPhiReweightedPt->SetMaximum(hPhiPhiReweightedPt->GetMaximum() * 1.3); // Add some space above the highest bin
    hPhiPhiReweightedPt->Draw();
    cPhiPhiReweightedPt->SaveAs("Generated_PhiPhiPtReweighted.png");

    //========================
    // Reweighted reconstructed distributions
    //========================

    //First lets draw hT4sMassRec_Injected_1D and hT4sPtRec_Injected_1D individually to see how it looks
    TCanvas *cPhiPhiMassRec_Injected = new TCanvas("cPhiPhiMassRec_Injected", "Reconstructed Mass Injected", 720, 720);
    SetCanvasStyle(cPhiPhiMassRec_Injected, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hT4sMassRec_Injected_1D);
    gPad->SetLogy();
    hT4sMassRec_Injected_1D->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    hT4sMassRec_Injected_1D->GetYaxis()->SetTitle("Counts");
    hT4sMassRec_Injected_1D->Rebin(5);
    hT4sMassRec_Injected_1D->Draw();
    cPhiPhiMassRec_Injected->SaveAs("Reconstructed_PhiPhiMassRec_Injected.png");

    TCanvas *cPhiPhiPtRec_Injected = new TCanvas("cPhiPhiPtRec_Injected", "Reconstructed Pt Injected", 720, 720);
    SetCanvasStyle(cPhiPhiPtRec_Injected, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(hT4sPtRec_Injected_1D);
    gPad->SetLogy();
    hT4sPtRec_Injected_1D->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hT4sPtRec_Injected_1D->GetYaxis()->SetTitle("Counts");
    hT4sPtRec_Injected_1D->Rebin(5);
    hT4sPtRec_Injected_1D->Draw();
    cPhiPhiPtRec_Injected->SaveAs("Reconstructed_PhiPhiPtRec_Injected.png");

    TH1D *h1_PhiPhiMassRecRecReweighted = (TH1D *)file->Get("phi-phi-m-c/Reconstructed/AllInjectedPhi/h1_PhiPhi_MassReweighted");
    TH1D *h1_PhiPhiPtRecRecReweighted = (TH1D *)file->Get("phi-phi-m-c/Reconstructed/AllInjectedPhi/h1_PhiPhi_PtReweighted");

    TCanvas *cPhiPhiMassRecReweighted = new TCanvas("cPhiPhiMassRecReweighted", "Reconstructed Mass Reweighted", 720, 720);
    SetCanvasStyle(cPhiPhiMassRecReweighted, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(h1_PhiPhiMassRecRecReweighted);
    gPad->SetLogy();
    h1_PhiPhiMassRecRecReweighted->Sumw2();
    h1_PhiPhiMassRecRecReweighted->GetXaxis()->SetTitle("#it{M}_{#phi#phi} (GeV/#it{c}^{2})");
    h1_PhiPhiMassRecRecReweighted->GetYaxis()->SetTitle("Counts");
    h1_PhiPhiMassRecRecReweighted->Rebin(20);
    h1_PhiPhiMassRecRecReweighted->GetXaxis()->SetRangeUser(2.0, 5.0); // Set the x-axis range to focus on the region of interest
    h1_PhiPhiMassRecRecReweighted->SetMinimum(1);
    h1_PhiPhiMassRecRecReweighted->Draw();
    cPhiPhiMassRecReweighted->SaveAs("Reconstructed_PhiPhiMassRecReweighted.png");

    TCanvas *cPhiPhiPtRecReweighted = new TCanvas("cPhiPhiPtRecReweighted", "Reconstructed Pt Reweighted", 720, 720);
    SetCanvasStyle(cPhiPhiPtRecReweighted, 0.16, 0.06, 0.06, 0.13);
    SetHistoQA(h1_PhiPhiPtRecRecReweighted);
    gPad->SetLogy();
    h1_PhiPhiPtRecRecReweighted->Sumw2();
    h1_PhiPhiPtRecRecReweighted->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h1_PhiPhiPtRecRecReweighted->GetYaxis()->SetTitle("Counts");
    h1_PhiPhiPtRecRecReweighted->Rebin(5);
    h1_PhiPhiPtRecRecReweighted->Draw();
    cPhiPhiPtRecReweighted->SaveAs("Reconstructed_PhiPhiPtRecReweighted.png");

    
}