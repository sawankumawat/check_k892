#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

using namespace std;

void efficiency()
{
    bool makeQAplots = true;   // set to false if you do not want to make QA plots
    string outputtype = "png"; // pdf, eps
    // ****************Data files ********************
    //// This is with wrong placement of TPC crossed rows
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/448490/kstarqa/hInvMass"; // 2022 data
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/449695/kstarqa/hInvMass"; // 2023 data
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/451993/kstarqa/hInvMass"; // 2024 data

    // correct placement of TPC crossed rows
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459845/kstarqa/hInvMass"; // 2022 data
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459908/kstarqa/hInvMass"; // 2023 data
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/460233/kstarqa/hInvMass"; // 2024 data

    // IR study
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/463114/kstarqa/hInvMass"; // 1-2 MHz
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/536899/kstarqa/hInvMass"; // 1-1.3 MHz
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/537861/kstarqa/hInvMass"; // 2 MHz
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535069/kstarqa/hInvMass"; // 14 kHz
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535545/kstarqa/hInvMass"; // 70 kHz
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/535645/kstarqa/hInvMass"; // 135 kHz
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/LHC23z/kstarqa/hInvMass"; // 450 kHz
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/LHC23ls/kstarqa/hInvMass"; // 650 kHz
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/466180/kstarqa_id33593/hInvMass"; // 2024 data (500 kHz)
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/459908/kstarqa/hInvMass"; // 2023 data (135 kHz)

    //*******************Occupancy cut study********************************
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/Occupancy_effect/466154/kstarqa_OCC500_id33931/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/Occupancy_effect/463286/kstarqa/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/Occupancy_effect/463171/kstarqa/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

    //*****Checks in data: Cuts on Mothers particle, only TPC, PID Ka2, pT dependent DCAxy, RCT ***********
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/468837/kstarqa/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/468791/kstarqa/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/468697/kstarqa/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

    //*****Checks in data: MinClusters ITS > 5, PbPb cuts, Pt dependent PID*******************
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/472720/kstarqa/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/470913/kstarqa/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/471898/kstarqa/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

    //*******Checks in data: Fake tracks, MID, TPCChi2MinCut, TrackRapidityCut *****************
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/473153/kstarqa_MID_id33593/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
    //  string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/473237/kstarqa_MID_id33593/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
    string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/473185/kstarqa_MID_id33593/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

    TString outputfolder = data_path + "/efficiency";
    gSystem->mkdir(outputfolder, kTRUE);

    //********MC files *****************
    // Wrong placement of TPC crossed rows
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/448748.root", "READ"); // 2022 MC (INEL)
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/449625.root", "READ"); // 2023 MC (INEL)
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/452215.root", "READ"); // 2024 MC  (INEL > 0)

    // ******************TPC tune on data without NN ***********************************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/459909.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/460571.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/459912.root", "READ"); // 2024 MC

    // ***********************TPC tune on data with NN ***********************************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/TuneOnDataWithNN/465413.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/TuneOnDataWithNN/464570.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/TuneOnDataWithNN/465419.root", "READ"); // 2023 MC

    //*********MC without tune on data (with NN)****************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/NotTuneOnData/NN/462343.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/NotTuneOnData/NN/462503.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/NotTuneOnData/NN/462351.root", "READ"); // 2024 MC

    //*********MC without tune on data (without NN)****************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/NotTuneOnData/464270.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/NotTuneOnData/464356.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/NotTuneOnData/464278.root", "READ"); // 2024 MC

    // ****************IR study*************************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/463114MC.root", "READ"); // 1-2 MHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/536899MC.root", "READ"); // 1-1.3 MHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/537861MC.root", "READ"); // 2 MHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/535069MC.root", "READ"); // 14 kHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/535545MC.root", "READ"); // 70 kHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/535645MC.root", "READ"); // 135 kHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/LHC23zMC.root", "READ"); // 450 kHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/LHC23lsMC.root", "READ"); // 650 kHz
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/467403.root", "READ"); // 2024 data (500 kHz)
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/460571.root", "READ"); // 2024 data (135 kHz)

    // *****************Occupancy cut study*************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Occupancy_effect/466809.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Occupancy_effect/463585.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Occupancy_effect/464277.root", "READ"); // 2024 MC

    // *********************All Tuning check********************
    //************2022 MC**************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/467976.root", "READ"); // without tune and NN
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/467975.root", "READ"); // With NN
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/467977.root", "READ"); // With tune on data
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/467978.root", "READ"); // With NN and tune on data

    //************2023 MC**************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/468146.root", "READ"); // without tune and NN
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/468149.root", "READ"); // With NN
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/468147.root", "READ"); // With tune on data
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/468148.root", "READ"); // With NN and tune on data

    //*************2024 MC***************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/467981.root", "READ"); // without tune and NN
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/467980.root", "READ"); // With NN
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/467982.root", "READ"); // With tune on data
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/AllTuningCheck/467983.root", "READ"); // With NN and tune on data

    //********************MC closure test*************************/
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Occupancy_effect/466809.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Occupancy_effect/463585.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Occupancy_effect/464277.root", "READ"); // 2024 MC

    //**************Checks (Default, FakeTrack, MID, PIDKa2) and new MC closure test********************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/checks/472785.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/checks/472820.root", "READ"); // 2023 MC
    TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/checks/474287.root", "READ"); // 2024 MC
    //*******************************************************************************************************

    TFile *fileraw = new TFile((data_path + "/yield.root").c_str(), "READ");

    if (fileeff->IsZombie() || fileraw->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }
    const string genpath = "kstarqa_MID_id34109/hInvMass/hk892GenpT";
    const string recpath = "kstarqa_MID_id34109/hInvMass/h2KstarRecpt2";
    // const string genpath = "kstarqa_OCC500_id34026/hInvMass/hk892GenpT";
    // const string recpath = "kstarqa_OCC500_id34026/hInvMass/h2KstarRecpt2";
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
    TFile *spectra = new TFile((data_path + "/corrected_spectra.root").c_str(), "RECREATE");
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
    TF1 *pol1fit = new TF1("pol1fit", "pol1", 7.0, 14);
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

    if (makeQAplots)
    {
        string QaPath = genpath.substr(0, genpath.length() - 20);
        cout << "QaPath: " << QaPath << endl;
        TString output_QA_folder = outputfolder + "/QA";
        gSystem->mkdir(output_QA_folder, kTRUE);
        TH1F *hNsigmaTOFKaon_neg = (TH1F *)fileeff->Get(Form("%s/hPID/Before/h1PID_TOF_neg_kaon", QaPath.c_str()));
        TH1F *hNsigmaTOFPion_neg = (TH1F *)fileeff->Get(Form("%s/hPID/Before/h1PID_TOF_neg_pion", QaPath.c_str()));
        TH1F *hNsigmaTOFKaon_pos = (TH1F *)fileeff->Get(Form("%s/hPID/Before/h1PID_TOF_pos_kaon", QaPath.c_str()));
        TH1F *hNsigmaTOFPion_pos = (TH1F *)fileeff->Get(Form("%s/hPID/Before/h1PID_TOF_pos_pion", QaPath.c_str()));
        TH1F *hNsigmaTPCKaon_neg = (TH1F *)fileeff->Get(Form("%s/hPID/Before/h1PID_TPC_neg_kaon", QaPath.c_str()));
        TH1F *hNsigmaTPCPion_neg = (TH1F *)fileeff->Get(Form("%s/hPID/Before/h1PID_TPC_neg_pion", QaPath.c_str()));
        TH1F *hNsigmaTPCKaon_pos = (TH1F *)fileeff->Get(Form("%s/hPID/Before/h1PID_TPC_pos_kaon", QaPath.c_str()));
        TH1F *hNsigmaTPCPion_pos = (TH1F *)fileeff->Get(Form("%s/hPID/Before/h1PID_TPC_pos_pion", QaPath.c_str()));
        if (hNsigmaTOFKaon_neg == nullptr || hNsigmaTOFPion_neg == nullptr || hNsigmaTOFKaon_pos == nullptr || hNsigmaTOFPion_pos == nullptr || hNsigmaTPCKaon_neg == nullptr || hNsigmaTPCPion_neg == nullptr || hNsigmaTPCKaon_pos == nullptr || hNsigmaTPCPion_pos == nullptr)
        {
            cerr << "PID histograms before selection not found!!!!!!!!!!!!" << endl;
            return;
        }

        TH2F *hNsigmaTPCTOFKaon = (TH2F *)fileeff->Get(Form("%s/hPID/Before/hNsigma_TPC_TOF_Ka_before", QaPath.c_str()));
        TH2F *hNsigmaTPCTOFPion = (TH2F *)fileeff->Get(Form("%s/hPID/Before/hNsigma_TPC_TOF_Pi_before", QaPath.c_str()));
        TH2F *hNsigmaTPCKaon = (TH2F *)fileeff->Get(Form("%s/hPID/Before/hNsigmaTPC_Ka_before", QaPath.c_str()));
        TH2F *hNsigmaTPCPion = (TH2F *)fileeff->Get(Form("%s/hPID/Before/hNsigmaTPC_Pi_before", QaPath.c_str()));
        TH2F *hNsigmaTOFKaon = (TH2F *)fileeff->Get(Form("%s/hPID/Before/hNsigmaTOF_Ka_before", QaPath.c_str()));
        TH2F *hNsigmaTOFPion = (TH2F *)fileeff->Get(Form("%s/hPID/Before/hNsigmaTOF_Pi_before", QaPath.c_str()));
        if (hNsigmaTPCTOFKaon == nullptr || hNsigmaTPCTOFPion == nullptr || hNsigmaTPCKaon == nullptr || hNsigmaTPCPion == nullptr || hNsigmaTOFKaon == nullptr || hNsigmaTOFPion == nullptr)
        {
            cerr << "2D PID histograms before selection not found!!!!!!!!!!!!" << endl;
            return;
        }

        TCanvas *cNsigmaTPCTOFKaon = new TCanvas("cNsigmaTPCTOFKaon", "Nsigma TPC TOF Kaon", 720, 720);
        SetCanvasStyle(cNsigmaTPCTOFKaon, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTPCTOFKaon);
        hNsigmaTPCTOFKaon->GetXaxis()->SetTitle("n#sigma_{TPC}");
        hNsigmaTPCTOFKaon->GetYaxis()->SetTitle("n#sigma_{TOF}");
        hNsigmaTPCTOFKaon->Draw("colz");
        // cNsigmaTPCTOFKaon->SaveAs(output_QA_folder + ("/NsigmaTPCTOFKaon." + outputtype).c_str());

        TCanvas *cNsigmaTPCTOFPion = new TCanvas("cNsigmaTPCTOFPion", "Nsigma TPC TOF Pion", 720, 720);
        SetCanvasStyle(cNsigmaTPCTOFPion, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTPCTOFPion);
        hNsigmaTPCTOFPion->GetXaxis()->SetTitle("n#sigma_{TPC}");
        hNsigmaTPCTOFPion->GetYaxis()->SetTitle("n#sigma_{TOF}");
        hNsigmaTPCTOFPion->Draw("colz");
        // cNsigmaTPCTOFPion->SaveAs(output_QA_folder + ("/NsigmaTPCTOFPion." + outputtype).c_str());

        TCanvas *cNsigmaTPCKaon = new TCanvas("cNsigmaTPCKaon", "Nsigma TPC Kaon", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTPCKaon);
        hNsigmaTPCKaon->GetYaxis()->SetTitle("n#sigma_{TPC}");
        hNsigmaTPCKaon->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hNsigmaTPCKaon->Draw("colz");
        // cNsigmaTPCKaon->SaveAs(output_QA_folder + ("/NsigmaTPCKaon." + outputtype).c_str());

        TCanvas *cNsigmaTPCPion = new TCanvas("cNsigmaTPCPion", "Nsigma TPC Pion", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTPCPion);
        hNsigmaTPCPion->GetYaxis()->SetTitle("n#sigma_{TPC}");
        hNsigmaTPCPion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hNsigmaTPCPion->Draw("colz");
        // cNsigmaTPCPion->SaveAs(output_QA_folder + ("/NsigmaTPCPion." + outputtype).c_str());

        TCanvas *cNsigmaTOFKaon = new TCanvas("cNsigmaTOFKaon", "Nsigma TOF Kaon", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTOFKaon);
        hNsigmaTOFKaon->GetYaxis()->SetTitle("n#sigma_{TOF}");
        hNsigmaTOFKaon->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hNsigmaTOFKaon->Draw("colz");
        // cNsigmaTOFKaon->SaveAs(output_QA_folder + ("/NsigmaTOFKaon." + outputtype).c_str());

        TCanvas *cNsigmaTOFPion = new TCanvas("cNsigmaTOFPion", "Nsigma TOF Pion", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTOFPion);
        hNsigmaTOFPion->GetYaxis()->SetTitle("n#sigma_{TOF}");
        hNsigmaTOFPion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hNsigmaTOFPion->Draw("colz");
        // cNsigmaTOFPion->SaveAs(output_QA_folder + ("/NsigmaTOFPion." + outputtype).c_str());

        TCanvas *cNsigmaTOFKaon_neg = new TCanvas("cNsigmaTOFKaon_neg", "Nsigma TOF Kaon Neg", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon_neg, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTOFKaon_neg);
        hNsigmaTOFKaon_neg->GetYaxis()->SetTitle("Counts");
        hNsigmaTOFKaon_neg->GetXaxis()->SetTitle("n#sigma_{TOF} (K^{-})");
        hNsigmaTOFKaon_neg->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTOFKaon_neg->Draw();
        TLine *lineverticalx0 = new TLine(0, hNsigmaTOFKaon_neg->GetMinimum(), 0, hNsigmaTOFKaon_neg->GetMaximum());
        lineverticalx0->SetLineStyle(2);
        lineverticalx0->SetLineColor(2);
        lineverticalx0->SetLineWidth(3);
        lineverticalx0->Draw();
        cNsigmaTOFKaon_neg->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_neg." + outputtype).c_str());

        TCanvas *cNsigmaTOFPion_neg = new TCanvas("cNsigmaTOFPion_neg", "Nsigma TOF Pion Neg", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion_neg, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTOFPion_neg);
        hNsigmaTOFPion_neg->GetYaxis()->SetTitle("Counts");
        hNsigmaTOFPion_neg->GetXaxis()->SetTitle("n#sigma_{TOF} (#pi^{-})");
        hNsigmaTOFPion_neg->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTOFPion_neg->Draw();
        lineverticalx0->SetY1(hNsigmaTOFPion_neg->GetMinimum());
        lineverticalx0->SetY2(hNsigmaTOFPion_neg->GetMaximum());
        lineverticalx0->Draw();
        cNsigmaTOFPion_neg->SaveAs(output_QA_folder + ("/NsigmaTOFPion_neg." + outputtype).c_str());

        TCanvas *cNsigmaTOFKaon_pos = new TCanvas("cNsigmaTOFKaon_pos", "Nsigma TOF Kaon Pos", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon_pos, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTOFKaon_pos);
        hNsigmaTOFKaon_pos->GetYaxis()->SetTitle("Counts");
        hNsigmaTOFKaon_pos->GetXaxis()->SetTitle("n#sigma_{TOF} (K^{+})");
        hNsigmaTOFKaon_pos->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTOFKaon_pos->Draw();
        lineverticalx0->SetY1(hNsigmaTOFKaon_pos->GetMinimum());
        lineverticalx0->SetY2(hNsigmaTOFKaon_pos->GetMaximum());
        lineverticalx0->Draw();
        cNsigmaTOFKaon_pos->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_pos." + outputtype).c_str());

        TCanvas *cNsigmaTOFPion_pos = new TCanvas("cNsigmaTOFPion_pos", "Nsigma TOF Pion Pos", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion_pos, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTOFPion_pos);
        hNsigmaTOFPion_pos->GetYaxis()->SetTitle("Counts");
        hNsigmaTOFPion_pos->GetXaxis()->SetTitle("n#sigma_{TOF} (#pi^{+})");
        hNsigmaTOFPion_pos->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTOFPion_pos->Draw();
        lineverticalx0->SetY1(hNsigmaTOFPion_pos->GetMinimum());
        lineverticalx0->SetY2(hNsigmaTOFPion_pos->GetMaximum());
        lineverticalx0->Draw();
        cNsigmaTOFPion_pos->SaveAs(output_QA_folder + ("/NsigmaTOFPion_pos." + outputtype).c_str());

        TCanvas *cNsigmaTPCKaon_neg = new TCanvas("cNsigmaTPCKaon_neg", "Nsigma TPC Kaon Neg", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon_neg, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTPCKaon_neg);
        hNsigmaTPCKaon_neg->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCKaon_neg->GetXaxis()->SetTitle("n#sigma_{TPC} (K^{-})");
        hNsigmaTPCKaon_neg->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTPCKaon_neg->Draw();
        lineverticalx0->SetY1(hNsigmaTPCKaon_neg->GetMinimum());
        lineverticalx0->SetY2(hNsigmaTPCKaon_neg->GetMaximum());
        lineverticalx0->Draw();
        cNsigmaTPCKaon_neg->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_neg." + outputtype).c_str());

        TCanvas *cNsigmaTPCPion_neg = new TCanvas("cNsigmaTPCPion_neg", "Nsigma TPC Pion Neg", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion_neg, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTPCPion_neg);
        hNsigmaTPCPion_neg->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCPion_neg->GetXaxis()->SetTitle("n#sigma_{TPC} (#pi^{-})");
        hNsigmaTPCPion_neg->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTPCPion_neg->Draw();
        lineverticalx0->SetY1(hNsigmaTPCPion_neg->GetMinimum());
        lineverticalx0->SetY2(hNsigmaTPCPion_neg->GetMaximum());
        lineverticalx0->Draw();
        cNsigmaTPCPion_neg->SaveAs(output_QA_folder + ("/NsigmaTPCPion_neg." + outputtype).c_str());

        TCanvas *cNsigmaTPCKaon_pos = new TCanvas("cNsigmaTPCKaon_pos", "Nsigma TPC Kaon Pos", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon_pos, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTPCKaon_pos);
        hNsigmaTPCKaon_pos->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCKaon_pos->GetXaxis()->SetTitle("n#sigma_{TPC} (K^{+})");
        hNsigmaTPCKaon_pos->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTPCKaon_pos->Draw();
        lineverticalx0->SetY1(hNsigmaTPCKaon_pos->GetMinimum());
        lineverticalx0->SetY2(hNsigmaTPCKaon_pos->GetMaximum());
        lineverticalx0->Draw();
        cNsigmaTPCKaon_pos->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_pos." + outputtype).c_str());

        TCanvas *cNsigmaTPCPion_pos = new TCanvas("cNsigmaTPCPion_pos", "Nsigma TPC Pion Pos", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion_pos, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTPCPion_pos);
        hNsigmaTPCPion_pos->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCPion_pos->GetXaxis()->SetTitle("n#sigma_{TPC} (#pi^{+})");
        hNsigmaTPCPion_pos->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTPCPion_pos->Draw();
        lineverticalx0->SetY1(hNsigmaTPCPion_pos->GetMinimum());
        lineverticalx0->SetY2(hNsigmaTPCPion_pos->GetMaximum());
        lineverticalx0->Draw();
        cNsigmaTPCPion_pos->SaveAs(output_QA_folder + ("/NsigmaTPCPion_pos." + outputtype).c_str());
    }
}