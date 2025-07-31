#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

using namespace std;

void efficiency()
{
    // ****************Data files ********************
    //// This is with wrong placement of TPC crossed rows
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/448490/kstarqa/hInvMass"; // 2022 data
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/449695/kstarqa/hInvMass"; // 2023 data
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/451993/kstarqa/hInvMass"; // 2024 data

    // correct placement of TPC crossed rows
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459845/kstarqa/hInvMass"; // 2022 data
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459908/kstarqa/hInvMass"; // 2023 data
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/460233/kstarqa/hInvMass"; // 2024 data

    // IR study
    // string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/463114/kstarqa/hInvMass"; // 1-2 MHz
    //  string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535069/kstarqa/hInvMass"; // 14 kHz
    //  string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535545/kstarqa/hInvMass"; // 70 kHz
    //  string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535645/kstarqa/hInvMass"; // 135 kHz
    //  string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535999/kstarqa/hInvMass"; // 330 kHz
    string path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/536106/kstarqa/hInvMass"; // 650 kHz

    TString outputfolder = path + "/efficiency";
    gSystem->mkdir(outputfolder, kTRUE);

    //********MC files *****************
    // Wrong placement of TPC crossed rows
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/448748.root", "READ"); // 2022 MC (INEL)
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/449625.root", "READ"); // 2023 MC (INEL)
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/452215.root", "READ"); // 2024 MC  (INEL > 0)

    // Correct placement of TPC crossed rows
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/459909.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/460571.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/459912.root", "READ"); // 2024 MC

    // IR study
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/463114MC.root", "READ"); // 1-2 MHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/535069MC.root", "READ"); // 14 kHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/535545MC.root", "READ"); // 70 kHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/535645MC.root", "READ"); // 135 kHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/535999MC.root", "READ"); // 330 kHz
    TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/536106MC.root", "READ"); // 650 kHz

    TFile *fileraw = new TFile((path + "/yield.root").c_str(), "READ");

    if (fileeff->IsZombie() || fileraw->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }
    const string genpath = "kstarqa/hInvMass/hk892GenpT";
    const string recpath = "kstarqa/hInvMass/h2KstarRecpt2";
    // const string genpath = "kstarqa_PIDKa2/hInvMass/hk892GenpT";
    // const string recpath = "kstarqa_PIDKa2/hInvMass/h2KstarRecpt2";

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    int nmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1; // number of multiplicity bins

    THnSparseF *hSpraseGen = (THnSparseF *)fileeff->Get(genpath.c_str());
    THnSparseF *hSparseRec = (THnSparseF *)fileeff->Get(recpath.c_str());

    if (hSpraseGen == nullptr || hSparseRec == nullptr)
    {
        cout << "Error reading efficiency histograms" << endl;
        return;
    }
    TFile *spectra = new TFile((path + "/corrected_spectra.root").c_str(), "RECREATE");
    TH1F *hChi2byNDF[nmultbins + 1];
    TH1F *hMass[nmultbins + 1];
    TH1F *hSignificance[nmultbins + 1];
    TH1F *heff[nmultbins + 1];
    int markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 32, 47};

    for (int imult = 0; imult < nmultbins + 1; imult++)
    {

        int multlow, multhigh;
        if (imult == 0)
        {
            multlow = 0;
            multhigh = 100; // for all multiplicity
        }
        else
        {
            multlow = mult_classes[imult - 1];
            multhigh = mult_classes[imult];
        }
        hChi2byNDF[imult] = (TH1F *)fileraw->Get(Form("mult_%d-%d/chi2byNDF", multlow, multhigh));
        hMass[imult] = (TH1F *)fileraw->Get(Form("mult_%d-%d/mass", multlow, multhigh));
        hSignificance[imult] = (TH1F *)fileraw->Get(Form("mult_%d-%d/significance", multlow, multhigh));

        int lowbinMultGen = hSpraseGen->GetAxis(1)->FindBin(multlow + 1e-3);
        int highbinMultGen = hSpraseGen->GetAxis(1)->FindBin(multhigh - 1e-3);
        hSpraseGen->GetAxis(1)->SetRange(lowbinMultGen, highbinMultGen);

        int lowbinMultRec = hSparseRec->GetAxis(1)->FindBin(multlow + 1e-3);
        int highbinMultRec = hSparseRec->GetAxis(1)->FindBin(multhigh - 1e-3);
        hSparseRec->GetAxis(1)->SetRange(lowbinMultRec, highbinMultRec);

        TH1D *h1gen = hSpraseGen->Projection(0, "E");
        TH1D *h1rec = hSparseRec->Projection(0, "E");

        TH1F *hyieldBinCount = (TH1F *)fileraw->Get(Form("mult_%d-%d/yield_bincount", multlow, multhigh));
        TH1F *hyieldIntegral = (TH1F *)fileraw->Get(Form("mult_%d-%d/yield_integral", multlow, multhigh));
        if (hyieldBinCount == nullptr)
        {
            cout << "Error reading yield histogram" << endl;
            return;
        }
        // heff[imult] = (TH1F *)hyieldBinCount->Clone(); // Just taking bins from the yield histogram which will be set to efficiency value.
        heff[imult] = (TH1F *)hyieldIntegral->Clone(); // Just taking bins from the yield histogram which will be set to efficiency value.

        cout << "bins: " << heff[imult]->GetNbinsX() << endl;
        for (int i = 0; i < heff[imult]->GetNbinsX(); i++)
        {
            lowpt = pT_bins[i];
            highpt = pT_bins[i + 1];
            cout << "lowpt: " << lowpt << " highpt: " << highpt << endl;

            double nrec = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt));
            double ngen = h1gen->Integral(h1gen->GetXaxis()->FindBin(lowpt), h1gen->GetXaxis()->FindBin(highpt));
            double efficiency = nrec / ngen;
            double efficiencyerr = sqrt(abs(((nrec + 1) / (ngen + 2)) * ((nrec + 2) / (ngen + 3) - (nrec + 1) / (ngen + 2))));

            cout << "Efficiency: " << efficiency << " +/- " << efficiencyerr << endl;
            heff[imult]->SetBinContent(i + 1, efficiency);
            heff[imult]->SetBinError(i + 1, efficiencyerr);
            hyieldBinCount->SetBinContent(i + 1, hyieldBinCount->GetBinContent(i + 1) / efficiency);
            hyieldIntegral->SetBinContent(i + 1, hyieldIntegral->GetBinContent(i + 1) / efficiency);

            double errorinyieldBinCount = hyieldBinCount->GetBinError(i + 1);
            double rawyieldvalueBinCount = hyieldBinCount->GetBinContent(i + 1);
            hyieldBinCount->SetBinError(i + 1, sqrt(pow(errorinyieldBinCount / efficiency, 2) + pow(rawyieldvalueBinCount * efficiencyerr / (efficiency * efficiency), 2)));

            double errorinyieldIntegral = hyieldIntegral->GetBinError(i + 1);
            double rawyieldvalueIntegral = hyieldIntegral->GetBinContent(i + 1);
            hyieldIntegral->SetBinError(i + 1, sqrt(pow(errorinyieldIntegral / efficiency, 2) + pow(rawyieldvalueIntegral * efficiencyerr / (efficiency * efficiency), 2)));
        }

        // TCanvas *c1 = new TCanvas("", "", 720, 720);
        // SetCanvasStyle(c1, 0.16, 0.06, 0.01, 0.14);
        // SetHistoQA(heff[imult]);
        // heff[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        // heff[imult]->GetYaxis()->SetTitle("Acceptance x Efficiency");
        // heff[imult]->GetYaxis()->SetTitleOffset(1.6);
        // heff[imult]->Draw();
        // c1->SaveAs(outputfolder + Form("/efficiency_%d_%d.png", multlow, multhigh));

        TCanvas *c2 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c2, 0.18, 0.06, 0.01, 0.14);
        gPad->SetLogy();
        SetHistoQA(hyieldBinCount);
        hyieldBinCount->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hyieldBinCount->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
        hyieldBinCount->GetYaxis()->SetTitleOffset(1.8);
        hyieldBinCount->Draw("ep");
        // c2->SaveAs(outputfolder + Form("/corrected_spectra_BinCount_%d_%d.png", multlow, multhigh));

        TCanvas *c3 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c3, 0.18, 0.06, 0.01, 0.14);
        gPad->SetLogy();
        SetHistoQA(hyieldIntegral);
        hyieldIntegral->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hyieldIntegral->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
        hyieldIntegral->GetYaxis()->SetTitleOffset(1.8);
        hyieldIntegral->Draw("ep");
        c3->SaveAs(outputfolder + Form("/corrected_spectra_Integral_%d_%d.png", multlow, multhigh));

        TDirectory *dir = spectra->mkdir(Form("mult_%d-%d", multlow, multhigh));
        dir->cd();
        heff[imult]->Write("heff");
        hyieldBinCount->Write("corrected_spectra_BinCount");
        hyieldIntegral->Write("corrected_spectra_Integral");
        spectra->cd();
    }
    gStyle->SetPalette(kRainBow);

    // Plot other plots for all multiplicity bins
    TCanvas *cefficiency = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cefficiency, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(heff[1]);
    heff[1]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    heff[1]->GetYaxis()->SetTitle("Acceptance x Efficiency");
    heff[1]->GetYaxis()->SetTitleOffset(1.6);
    heff[1]->SetMaximum(0.65);
    heff[1]->Draw("pe");
    for (int imult = 1; imult < nmultbins + 1; imult++)
    {
        heff[imult]->SetMarkerStyle(markers[imult]);
        heff[imult]->SetMarkerSize(1.2);
        heff[imult]->Draw("pe same PLC PMC");
    }
    TF1 *pol1fit = new TF1("pol1fit", "pol1", 5.5, 12);
    pol1fit->FixParameter(1, 0.0); // Fixing the slope to 0
    pol1fit->SetLineColor(kBlack);
    pol1fit->SetLineWidth(2);
    pol1fit->SetLineStyle(2);
    heff[4]->Fit(pol1fit, "RQ0");
    cout << "Intercept value : " << pol1fit->GetParameter(0) << endl;
    TLine *line = new TLine(0, pol1fit->GetParameter(0), 15, pol1fit->GetParameter(0));
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw();
    TLegend *legall = new TLegend(0.20, 0.75, 0.92, 0.98);
    legall->SetTextSize(0.03);
    legall->SetNColumns(5);
    legall->SetFillStyle(0);
    legall->SetBorderSize(0);
    // legall->AddEntry(heff[0], "0-100%", "p");
    for (int imult = 1; imult < nmultbins + 1; imult++)
    {
        legall->AddEntry(heff[imult], Form("%.0f-%.0f%%", mult_classes[imult - 1], mult_classes[imult]), "p");
    }
    legall->Draw();
    cefficiency->SaveAs(outputfolder + "/efficiency_all_mult.png");

    TCanvas *cSignificance = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cSignificance, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hSignificance[1]);
    hSignificance[1]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSignificance[1]->GetYaxis()->SetTitle("Significance");
    hSignificance[1]->GetYaxis()->SetTitleOffset(1.6);
    hSignificance[1]->SetMaximum(750);
    hSignificance[1]->SetMinimum(-5);
    hSignificance[1]->Draw("p");
    for (int imult = 1; imult < nmultbins + 1; imult++)
    {
        hSignificance[imult]->SetMarkerStyle(markers[imult]);
        hSignificance[imult]->SetMarkerSize(1.2);
        hSignificance[imult]->Draw("p same PLC PMC");
    }
    legall->Draw();
    cSignificance->SaveAs(outputfolder + "/significance_all_mult.png");

    TCanvas *cChi2byNDF = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cChi2byNDF, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hChi2byNDF[1]);
    hChi2byNDF[1]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hChi2byNDF[1]->GetYaxis()->SetTitle("#chi^{2}/NDF");
    hChi2byNDF[1]->GetYaxis()->SetTitleOffset(1.6);
    hChi2byNDF[1]->SetMaximum(6.5);
    hChi2byNDF[1]->SetMinimum(0);
    hChi2byNDF[1]->SetStats(0);
    hChi2byNDF[1]->Draw("p");
    for (int imult = 1; imult < nmultbins + 1; imult++)
    {
        hChi2byNDF[imult]->SetMarkerStyle(markers[imult]);
        hChi2byNDF[imult]->SetMarkerSize(1.2);
        hChi2byNDF[imult]->Draw("p same PLC PMC");
    }
    legall->Draw();
    cChi2byNDF->SaveAs(outputfolder + "/chi2byNDF_all_mult.png");

    TCanvas *cMass = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cMass, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hMass[1]);
    hMass[1]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass[1]->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass[1]->GetYaxis()->SetTitleOffset(1.6);
    hMass[1]->GetYaxis()->SetRangeUser(0.878, 0.919);
    hMass[1]->SetStats(0);
    hMass[1]->Draw("pe");
    for (int imult = 1; imult < nmultbins + 1; imult++)
    {
        hMass[imult]->SetMarkerStyle(markers[imult]);
        hMass[imult]->SetMarkerSize(1.2);
        hMass[imult]->Draw("pe same PLC PMC");
    }
    TLine *linePDG = new TLine(0, masspdg, 20, masspdg);
    linePDG->SetLineStyle(2);
    linePDG->SetLineColor(2);
    linePDG->SetLineWidth(2);
    linePDG->Draw();
    legall->AddEntry(linePDG, "PDG Mass", "l");
    legall->Draw();
    cMass->SaveAs(outputfolder + "/mass_all_mult.png");
}