#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

using namespace std;
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void efficiency()
{
    bool makePIDplots = false; // qa plots
    string outputtype = "png"; // pdf, eps

    int colors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

    // ****************Data files ********************

    // correct placement of TPC crossed rows
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459845/kstarqa/hInvMass"; // 2022 data
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/459908/kstarqa_PIDKa2/hInvMass"; // 2023 data
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/460233/kstarqa_PIDKa2/hInvMass"; // 2024 data

    //**********************************MID Cuts******************************************
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/477779/kstarqa_MID/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/478015/kstarqa_MID/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/477833/kstarqa/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

    //*************************PID Variations for Kaon (without MID)**************************
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/480317/kstarqa/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/480447/kstarqa/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/480657/kstarqa/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

    //================================QA checks==================
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/585940/kstarqa_VertexTOFMatched/hInvMass";

    //==============================Pt-dependent PID=======================
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/586976/kstarqa/hInvMass"; // 2023
    // string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/586385/kstarqa/hInvMass"; // 2024

    //================================Temporary=======================
    string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/589661/kstarqa_INEL/hInvMass";

    TString outputfolder = data_path + "/efficiency";
    gSystem->mkdir(outputfolder, kTRUE);

    //********MC files *****************

    // ***********************TPC tune on data with NN ***********************************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/TuneOnDataWithNN/465413.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/TuneOnDataWithNN/464570.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/TuneOnDataWithNN/465419.root", "READ"); // 2023 MC

    //********************MC closure test*************************/
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Occupancy_effect/466809.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Occupancy_effect/463585.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/Occupancy_effect/464277.root", "READ"); // 2024 MC

    //**************Checks (Default, FakeTrack, MID, PIDKa2) and new MC closure test********************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/checks/472785.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/checks/472820.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/checks/472788.root", "READ"); // 2024 MC

    //********************************MID Cuts******************************************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/478021.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/477593.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/478022.root", "READ"); // 2024 MC

    ////*************************PID Variations for Kaon (with MID)**************************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/480451.root", "READ"); // 2022 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/480565.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/480453.root", "READ"); // 2024 MC

    //**************************After Calibrated MC from Nicolo****************************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/483982.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/temp.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/484450.root", "READ"); // 2024 MC

    //*******************************************QA checks************************************************
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/585904.root", "READ"); // 2023 MC

    //================================Pt-dependent PID=======================
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/586966.root", "READ"); // 2023 MC
    // TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/586467.root", "READ"); // 2024 MC

    //================================Temporary=======================
    TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/589815.root", "READ"); // 2023 MC

    TFile *fileraw = new TFile((data_path + "/yield.root").c_str(), "READ");

    if (fileeff->IsZombie() || fileraw->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }
    const string genpath = "kstarqa_INEL/hInvMass";
    const string recpath = "kstarqa_INEL/hInvMass";
    // const string genpath = "kstarqa_OCC500_id34026/hInvMass/hk892GenpT";
    // const string recpath = "kstarqa_OCC500_id34026/hInvMass/h2KstarRecpt2";
    // const string genpath = "kstarqa_PIDKa2/hInvMass/hk892GenpT";
    // const string recpath = "kstarqa_PIDKa2/hInvMass/h2KstarRecpt2";

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    int nmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1; // number of multiplicity bins

    // THnSparseF *hSpraseGen = (THnSparseF *)fileeff->Get(Form("%s/hk892GenpT", genpath.c_str()));
    // THnSparseF *hSparseRec = (THnSparseF *)fileeff->Get(Form("%s/h2KstarRecpt2", recpath.c_str()));
    THnSparseF *hSpraseGen = (THnSparseF *)fileeff->Get(Form("%s/hk892GenpTCalib1", genpath.c_str()));
    THnSparseF *hSparseRec = (THnSparseF *)fileeff->Get(Form("%s/h2KstarRecptCalib2", recpath.c_str()));
    if (hSpraseGen == nullptr || hSparseRec == nullptr)
    {
        cout << "Error reading efficiency histograms" << endl;
        return;
    }

    // Event loss histograms
    // TH1F *hAllGenColl = (TH1F *)fileeff->Get(Form("%s/hAllGenCollisions", genpath.c_str()));
    // TH1F *hAllGenColl1Rec = (TH1F *)fileeff->Get(Form("%s/hAllGenCollisions1Rec", recpath.c_str()));
    TH1F *hAllGenColl = (TH1F *)fileeff->Get(Form("%s/MCcorrections/MultiplicityGen", genpath.c_str()));
    TH1F *hAllGenColl1Rec = (TH1F *)fileeff->Get(Form("%s/MCcorrections/MultiplicityRec", recpath.c_str())); // without event splitting correction
    TH1F *hAllGenColl1Rec_evsplit = (TH1F *)fileeff->Get(Form("%s/h1RecMult2", recpath.c_str()));            // with event splitting correction

    if (hAllGenColl == nullptr || hAllGenColl1Rec == nullptr)
    {
        cout << "Error reading event loss histograms" << endl;
        return;
    }

    // Signal loss histograms
    // TH2F *hAllGenKstar = (TH2F *)fileeff->Get(Form("%s/hAllKstarGenCollisisons", genpath.c_str()));
    // TH2F *hAllGenKstar1Rec = (TH2F *)fileeff->Get(Form("%s/hAllKstarGenCollisisons1Rec", recpath.c_str()));

    TH2F *hAllGenKstar = (TH2F *)fileeff->Get(Form("%s/MCcorrections/hSignalLossDenominator", genpath.c_str()));
    TH2F *hAllGenKstar1Rec = (TH2F *)fileeff->Get(Form("%s/MCcorrections/hSignalLossNumerator", recpath.c_str()));

    // // Get correlation of nchargeParticles vs multiplicity from MC
    // TH2F *hNchVsMultGen = (TH2F *)fileeff->Get(Form("%s/CorrFactors/hMultiplicityVsMultMC", genpath.c_str()));
    // if (hAllGenKstar == nullptr || hAllGenKstar1Rec == nullptr || hNchVsMultGen == nullptr)
    // {
    //     cout << "Error reading signal loss histograms" << endl;
    //     return;
    // }

    vector<double> nChParticlesFromMult;
    // int multLoopEnd = nmultbins + 1;
    int multLoopEnd = 1;
    int multlow, multhigh;

    TFile *spectra = new TFile((data_path + "/corrected_spectra.root").c_str(), "RECREATE");
    TH1F *hChi2byNDF[multLoopEnd];
    TH1F *hMass[multLoopEnd];
    TH1F *hWidth[multLoopEnd];
    TH1F *hSignificance[multLoopEnd];
    TH1F *heff[multLoopEnd];
    TH1F *hSignalLoss[multLoopEnd];
    int markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 32, 47};
    TH1F *heventloss[multLoopEnd];
    TH1D *h1gen;
    TH1D *h1rec;
    TH1F *hRatioEvBySig[multLoopEnd];

    for (int imult = 0; imult < multLoopEnd; imult++)
    {

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
        hWidth[imult] = (TH1F *)fileraw->Get(Form("mult_%d-%d/width", multlow, multhigh));
        hSignificance[imult] = (TH1F *)fileraw->Get(Form("mult_%d-%d/significance", multlow, multhigh));
        hSignalLoss[imult] = new TH1F(Form("signal_loss_%d", multhigh), "Signal Loss", Npt, pT_bins);

        int lowbinMultGen = hSpraseGen->GetAxis(1)->FindBin(multlow + 1e-3);
        int highbinMultGen = hSpraseGen->GetAxis(1)->FindBin(multhigh - 1e-3);
        hSpraseGen->GetAxis(1)->SetRange(lowbinMultGen, highbinMultGen);

        int lowbinMultRec = hSparseRec->GetAxis(1)->FindBin(multlow + 1e-3);
        int highbinMultRec = hSparseRec->GetAxis(1)->FindBin(multhigh - 1e-3);
        hSparseRec->GetAxis(1)->SetRange(lowbinMultRec, highbinMultRec);

        h1gen = hSpraseGen->Projection(0, "E");
        h1rec = hSparseRec->Projection(0, "E");

        TH1F *hyieldBinCount = (TH1F *)fileraw->Get(Form("mult_%d-%d/yield_bincount", multlow, multhigh));
        TH1F *hyieldIntegral = (TH1F *)fileraw->Get(Form("mult_%d-%d/yield_integral", multlow, multhigh));
        if (hyieldBinCount == nullptr)
        {
            cout << "Error reading yield histogram" << endl;
            return;
        }
        // heff[imult] = (TH1F *)hyieldBinCount->Clone(); // Just taking bins from the yield histogram which will be set to efficiency value.
        heff[imult] = (TH1F *)hyieldIntegral->Clone(); // Just taking bins from the yield histogram which will be set to efficiency value.

        // Event loss calculations
        heventloss[imult] = new TH1F(Form("hEventLoss_%d", imult), "Event Loss", Npt, pT_bins);

        double eventLossNum = hAllGenColl1Rec->Integral(hAllGenColl1Rec->GetXaxis()->FindBin(multlow + 0.001), hAllGenColl1Rec->GetXaxis()->FindBin(multhigh - 0.001)); // without event splitting consideration

        // double eventLossNum = hAllGenColl1Rec_evsplit->Integral(hAllGenColl1Rec_evsplit->GetXaxis()->FindBin(multlow + 0.001), hAllGenColl1Rec_evsplit->GetXaxis()->FindBin(multhigh - 0.001)); // with event splitting consideration (***Need to be fixed***)

        double eventLossDen = hAllGenColl->Integral(hAllGenColl->GetXaxis()->FindBin(multlow + 0.001), hAllGenColl->GetXaxis()->FindBin(multhigh - 0.001));

        // TH1F *hNchInMult = (TH1F *)hNchVsMultGen->ProjectionY(Form("hNchInMult_%d", imult), hNchVsMultGen->GetXaxis()->FindBin(multlow + 0.01), hNchVsMultGen->GetXaxis()->FindBin(multhigh - 0.01));
        // hAllGenColl1Rec->Multiply(hAllGenColl1Rec, hNchInMult);
        // hAllGenColl->Multiply(hAllGenColl, hNchInMult);
        // double eventLossNum = hAllGenColl1Rec->Integral();
        // double eventLossDen = hAllGenColl->Integral();
        // cout << "Multiplicity bin " << multlow << "-" << multhigh << ", Event Loss " << eventLossNum / eventLossDen << endl;

        // Signal loss calculations
        TH1F *hSignalLossNumPt = (TH1F *)hAllGenKstar1Rec->ProjectionX(Form("SignalLossNumPt_%d", imult), hAllGenKstar1Rec->GetYaxis()->FindBin(multlow + 0.01), hAllGenKstar1Rec->GetYaxis()->FindBin(multhigh - 0.01));
        TH1F *hSignalLossDenPt = (TH1F *)hAllGenKstar->ProjectionX(Form("SignalLossDenPt_%d", imult), hAllGenKstar->GetYaxis()->FindBin(multlow + 0.01), hAllGenKstar->GetYaxis()->FindBin(multhigh - 0.01));

        // TH1F *hSignalLossNumPt = (TH1F *)hAllGenKstar1Rec->ProjectionX();
        // TH1F *hSignalLossDenPt = (TH1F *)hAllGenKstar->ProjectionX();

        // cout << "bins: " << heff[imult]->GetNbinsX() << endl;
        for (int i = 0; i < heff[imult]->GetNbinsX(); i++)
        {
            lowpt = pT_bins[i];
            highpt = pT_bins[i + 1];
            // cout << "lowpt: " << lowpt << " highpt: " << highpt << endl;

            double nrec = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt));
            double ngen = h1gen->Integral(h1gen->GetXaxis()->FindBin(lowpt), h1gen->GetXaxis()->FindBin(highpt));
            double efficiency = nrec / ngen;
            double efficiencyerr = sqrt(abs(((nrec + 1) / (ngen + 2)) * ((nrec + 2) / (ngen + 3) - (nrec + 1) / (ngen + 2))));

            // cout << "Efficiency: " << efficiency << " +/- " << efficiencyerr << endl;
            heff[imult]->SetBinContent(i + 1, efficiency);
            heff[imult]->SetBinError(i + 1, efficiencyerr);
            hyieldBinCount->SetBinContent(i + 1, hyieldBinCount->GetBinContent(i + 1) / (efficiency * 2));
            hyieldIntegral->SetBinContent(i + 1, hyieldIntegral->GetBinContent(i + 1) / (efficiency * 2));
            // Note: the efficiency is multiplied by 2, because in run2 the average of K* + anit-K* was taken.

            double errorinyieldBinCount = hyieldBinCount->GetBinError(i + 1);
            double rawyieldvalueBinCount = hyieldBinCount->GetBinContent(i + 1);
            hyieldBinCount->SetBinError(i + 1, 0.5 * sqrt(pow(errorinyieldBinCount / efficiency, 2) + pow(rawyieldvalueBinCount * efficiencyerr / (efficiency * efficiency), 2))); // 0.5 for average of K*+ and K*-

            double errorinyieldIntegral = hyieldIntegral->GetBinError(i + 1);
            double rawyieldvalueIntegral = hyieldIntegral->GetBinContent(i + 1);
            hyieldIntegral->SetBinError(i + 1, 0.5 * sqrt(pow(errorinyieldIntegral / efficiency, 2) + pow(rawyieldvalueIntegral * efficiencyerr / (efficiency * efficiency), 2))); // 0.5 for average of K*+ and K*-
            // Event loss calculations
            // if (imult > 0)
            {
                heventloss[imult]->SetBinContent(i + 1, eventLossNum / eventLossDen);
                heventloss[imult]->SetBinError(i + 1, 0);
                // cout << "Event Loss (Multiplicity " << multlow << "-" << multhigh << "): " << heventloss[imult]->GetBinContent(imult) << endl;
            }

            // Signal loss calculations
            // if (imult > 0)
            {
                double signalLossNum = hSignalLossNumPt->Integral(hSignalLossNumPt->GetXaxis()->FindBin(lowpt + 0.001), hSignalLossNumPt->GetXaxis()->FindBin(highpt - 0.001));
                double signalLossDen = hSignalLossDenPt->Integral(hSignalLossDenPt->GetXaxis()->FindBin(lowpt + 0.001), hSignalLossDenPt->GetXaxis()->FindBin(highpt - 0.001));
                double signalLoss = signalLossNum / signalLossDen;
                cout << "signal loss value is " << signalLoss << endl;
                hSignalLoss[imult]->SetBinContent(i + 1, signalLoss);
                hSignalLoss[imult]->SetBinError(i + 1, 0);
            }
        }
        hRatioEvBySig[imult] = (TH1F *)heventloss[imult]->Clone(Form("hRatioEvBySig_%d", imult));
        hRatioEvBySig[imult]->Divide(hSignalLoss[imult]);

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

        TH1F *hlossCorrected_integral = (TH1F *)hyieldIntegral->Clone("hlossCorrected_integral");
        hlossCorrected_integral->Multiply(hRatioEvBySig[imult]);
        TH1F *hlossCorrected_bincount = (TH1F *)hyieldBinCount->Clone("hlossCorrected_bincount");
        hlossCorrected_bincount->Multiply(hRatioEvBySig[imult]);

        TDirectory *dir = spectra->mkdir(Form("mult_%d-%d", multlow, multhigh));
        dir->cd();
        heff[imult]->Write("heff");
        hyieldBinCount->Write("corrected_spectra_BinCount");
        hyieldIntegral->Write("corrected_spectra_Integral");
        hlossCorrected_integral->Write("corrected_spectra_Integral_final");
        hlossCorrected_bincount->Write("corrected_spectra_BinCount_final");
        spectra->cd();
    }
    gStyle->SetPalette(kRainBow);

    // TCanvas *cGeneratedMult = new TCanvas("", "", 720, 720);
    // SetCanvasStyle(cGeneratedMult, 0.16, 0.06, 0.01, 0.14);
    // SetHistoQA(h1gen);
    // h1gen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // h1gen->GetYaxis()->SetTitle("Generated Events");
    // h1gen->GetYaxis()->SetTitleOffset(1.6);
    // h1gen->Draw("ep");
    // cGeneratedMult->SaveAs(outputfolder + "/generated_pT.png");

    // TCanvas *cReconstructedMult = new TCanvas("", "", 720, 720);
    // SetCanvasStyle(cReconstructedMult, 0.16, 0.06, 0.01, 0.14);
    // SetHistoQA(h1rec);
    // h1rec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // h1rec->GetYaxis()->SetTitle("Reconstructed Events");
    // h1rec->GetYaxis()->SetTitleOffset(1.6);
    // h1rec->Draw("ep");
    // cReconstructedMult->SaveAs(outputfolder + "/reconstructed_pT.png");

    // //************Take ratio of min bias efficiency with any other efficiency file********************
    // TFile *fileEff2 = new TFile("/home/sawan/check_k892/output/kstar/LHC22o_pass7/480657/kstarqa_PIDKa1/hInvMass/efficiency/beforeCalibration/corrected_spectra.root"); // File for efficiency comparison
    // TH1F *hEffMinBias = (TH1F *)fileEff2->Get("mult_0-100/heff");
    // TH1F *hratio = (TH1F *)hEffMinBias->Clone("hratio");
    // hratio->Divide(heff[0]);
    // TCanvas *cRatioEff = new TCanvas("", "", 720, 720);
    // double pad1Size, pad2Size;
    // SetCanvasStyle(cRatioEff, 0.25, 0.03, 0.03, 0.15);
    // canvas_style(cRatioEff, pad1Size, pad2Size);
    // cRatioEff->cd(1);
    // SetHistoQA(heff[0]);
    // SetHistoQA(hEffMinBias);
    // heff[0]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    // heff[0]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    // heff[0]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    // heff[0]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    // heff[0]->SetMaximum(heff[0]->GetMaximum() * 1.2);
    // heff[0]->Draw("ep");
    // hEffMinBias->SetMarkerStyle(20);
    // hEffMinBias->SetMarkerColor(kRed);
    // hEffMinBias->SetLineColor(kRed);
    // hEffMinBias->Draw("ep same");
    // TLegend *legEffRatio = new TLegend(0.20, 0.75, 0.92, 0.98);
    // legEffRatio->SetTextSize(0.03);
    // legEffRatio->SetFillStyle(0);
    // legEffRatio->SetBorderSize(0);
    // legEffRatio->AddEntry(heff[0], "Efficiency 1", "p");
    // legEffRatio->AddEntry(hEffMinBias, "Efficiency 2", "p");
    // legEffRatio->Draw();

    // cRatioEff->cd(2);
    // hratio->GetYaxis()->SetTitleSize(0.026 / pad2Size);
    // hratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    // hratio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    // hratio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    // hratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // hratio->GetYaxis()->SetTitle("Ratio (Eff2 / Eff1)");
    // hratio->GetYaxis()->SetTitleOffset(0.75);
    // hratio->Draw();
    // TLine *line1 = new TLine(0, 1, 15, 1);
    // line1->SetLineColor(kBlack);
    // line1->SetLineStyle(2);
    // line1->Draw("same");
    // cRatioEff->SaveAs(outputfolder + "/efficiency_ratio.png");

    // Plot other plots for all multiplicity bins
    TCanvas *cefficiency = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cefficiency, 0.16, 0.06, 0.01, 0.14);
    for (int imult = 0; imult < multLoopEnd; imult++)
    {
        SetHistoQA(heff[imult]);
        heff[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        heff[imult]->GetYaxis()->SetTitle("Acceptance x Efficiency");
        heff[imult]->GetYaxis()->SetTitleOffset(1.6);
        heff[imult]->SetMaximum(0.65);
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
    for (int imult = 1; imult < multLoopEnd; imult++)
    {
        legall->AddEntry(heff[imult], Form("%.0f-%.0f%%", mult_classes[imult - 1], mult_classes[imult]), "p");
    }
    legall->Draw();
    cefficiency->SaveAs(outputfolder + "/efficiency_all_mult.png");

    // Generate colors
    int allColors[11];
    int numColors = gStyle->GetNumberOfColors();
    for (int icol = 0; icol < 11; icol++)
    {
        int paletteIndex = (icol)*numColors / 11;             // total 11 colors are assumed
        paletteIndex = std::min(paletteIndex, numColors - 1); // Ensure within bounds
        int color = gStyle->GetColorPalette(paletteIndex);
        allColors[icol] = color;
    }

    TCanvas *cEventLoss = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cEventLoss, 0.16, 0.06, 0.01, 0.14);
    legall->Clear();
    for (int imult = 0; imult < multLoopEnd; imult++)
    {
        if (imult == 0)
        {
            multlow = 0;
            multhigh = 100;
        }
        else
        {
            multlow = mult_classes[imult - 1];
            multhigh = mult_classes[imult];
        }
        SetHistoQA(heventloss[imult]);
        heventloss[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        heventloss[imult]->GetYaxis()->SetTitle("Event Loss");
        heventloss[imult]->GetYaxis()->SetTitleOffset(1.6);
        heventloss[imult]->SetStats(0);
        heventloss[imult]->SetMarkerStyle(markers[imult]);
        heventloss[imult]->SetMarkerSize(1.2);
        heventloss[imult]->SetLineColor(colors[imult]);
        heventloss[imult]->Draw("l same");
        legall->AddEntry(heventloss[imult], Form("%d-%d%%", multlow, multhigh), "l");
    }
    legall->Draw();
    // TLine *lineEventLoss = new TLine(0, 1.0, 100, 1.0);
    // lineEventLoss->SetLineStyle(2);
    // lineEventLoss->SetLineColor(2);
    // lineEventLoss->SetLineWidth(2);
    // lineEventLoss->Draw();
    cEventLoss->SaveAs(outputfolder + "/event_loss.png");

    TCanvas *cSignalLoss = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cSignalLoss, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hSignalLoss[0]);
    hSignalLoss[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSignalLoss[0]->GetYaxis()->SetTitle("Event and signal loss");
    hSignalLoss[0]->GetYaxis()->SetTitleOffset(1.6);
    hSignalLoss[0]->SetStats(0);
    hSignalLoss[0]->SetMaximum(1.4);
    hSignalLoss[0]->SetMinimum(0.3);
    hSignalLoss[0]->Draw("pe");
    legall->Clear();
    for (int imult = 0; imult < multLoopEnd; imult++)
    {
        if (imult == 0)
        {
            multlow = 0;
            multhigh = 100;
        }
        else
        {
            multlow = mult_classes[imult - 1];
            multhigh = mult_classes[imult];
        }
        hSignalLoss[imult]->SetMarkerStyle(markers[imult]);
        hSignalLoss[imult]->SetMarkerSize(1.2);
        hSignalLoss[imult]->SetLineColor(colors[imult]);
        hSignalLoss[imult]->SetMarkerColor(colors[imult]);
        hSignalLoss[imult]->Draw("pe same");
        heventloss[imult]->SetLineStyle(2);
        heventloss[imult]->SetLineWidth(2);
        heventloss[imult]->SetLineColor(colors[imult]);
        heventloss[imult]->Draw("l same");
        legall->AddEntry(hSignalLoss[imult], Form("%d-%d%%", multlow, multhigh), "p");
    }
    legall->AddEntry(heventloss[0], "Event Loss", "l");
    legall->Draw();
    cSignalLoss->SaveAs(outputfolder + "/signal_loss_all_mult.png");

    TCanvas *cEventBySignalLoss = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cEventBySignalLoss, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hRatioEvBySig[0]);
    hRatioEvBySig[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatioEvBySig[0]->GetYaxis()->SetTitle("Event Loss / Signal Loss");
    hRatioEvBySig[0]->GetYaxis()->SetTitleOffset(1.6);
    hRatioEvBySig[0]->SetStats(0);
    hRatioEvBySig[0]->SetMaximum(1.25);
    hRatioEvBySig[0]->SetMinimum(0.6);
    hRatioEvBySig[0]->Draw("pe");
    legall->Clear();
    for (int imult = 0; imult < multLoopEnd; imult++)
    {
        if (imult == 0)
        {
            multlow = 0;
            multhigh = 100;
        }
        else
        {
            multlow = mult_classes[imult - 1];
            multhigh = mult_classes[imult];
        }
        hRatioEvBySig[imult]->SetMarkerStyle(markers[imult]);
        hRatioEvBySig[imult]->SetMarkerSize(1.2);
        hRatioEvBySig[imult]->SetMarkerColor(colors[imult]);
        hRatioEvBySig[imult]->SetLineColor(colors[imult]);
        hRatioEvBySig[imult]->Draw("pe same");
        legall->AddEntry(hRatioEvBySig[imult], Form("%d-%d%%", multlow, multhigh), "p");
    }
    legall->Draw();
    cEventBySignalLoss->SaveAs(outputfolder + "/event_by_signal_loss.png");

    TCanvas *cSignificance = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cSignificance, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hSignificance[0]);
    hSignificance[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSignificance[0]->GetYaxis()->SetTitle("Significance");
    hSignificance[0]->GetYaxis()->SetTitleOffset(1.6);
    hSignificance[0]->SetMaximum(750);
    hSignificance[0]->SetMinimum(-5);
    hSignificance[0]->Draw("p");
    for (int imult = 0; imult < multLoopEnd; imult++)
    {
        hSignificance[imult]->SetMarkerStyle(markers[imult]);
        hSignificance[imult]->SetMarkerSize(1.2);
        hSignificance[imult]->Draw("p same PLC PMC");
    }
    legall->Draw();
    cSignificance->SaveAs(outputfolder + "/significance_all_mult.png");

    TCanvas *cChi2byNDF = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cChi2byNDF, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hChi2byNDF[0]);
    hChi2byNDF[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hChi2byNDF[0]->GetYaxis()->SetTitle("#chi^{2}/NDF");
    hChi2byNDF[0]->GetYaxis()->SetTitleOffset(1.6);
    hChi2byNDF[0]->SetMaximum(6.5);
    hChi2byNDF[0]->SetMinimum(0);
    hChi2byNDF[0]->SetStats(0);
    hChi2byNDF[0]->Draw("p");
    for (int imult = 0; imult < multLoopEnd; imult++)
    {
        hChi2byNDF[imult]->SetMarkerStyle(markers[imult]);
        hChi2byNDF[imult]->SetMarkerSize(1.2);
        hChi2byNDF[imult]->Draw("p same PLC PMC");
    }
    legall->Draw();
    cChi2byNDF->SaveAs(outputfolder + "/chi2byNDF_all_mult.png");

    TCanvas *cMass = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cMass, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hMass[0]);
    hMass[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass[0]->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass[0]->GetYaxis()->SetTitleOffset(1.6);
    hMass[0]->GetYaxis()->SetRangeUser(0.878, 0.919);
    hMass[0]->SetStats(0);
    hMass[0]->Draw("pe");
    for (int imult = 0; imult < multLoopEnd; imult++)
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

    TCanvas *cWidth = new TCanvas("", "", 720, 720);
    SetCanvasStyle(cWidth, 0.16, 0.06, 0.01, 0.14);
    SetHistoQA(hWidth[0]);
    hWidth[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hWidth[0]->GetYaxis()->SetTitle("Width (GeV/#it{c}^{2})");
    hWidth[0]->GetYaxis()->SetTitleOffset(1.6);
    hWidth[0]->GetYaxis()->SetRangeUser(0.047 - 0.04, 0.047 + 0.05);
    hWidth[0]->SetStats(0);
    hWidth[0]->Draw("pe");
    TLine *linePDGWidth = new TLine(0, widthpdg, 20, widthpdg);
    linePDGWidth->SetLineStyle(2);
    linePDGWidth->SetLineColor(2);
    linePDGWidth->SetLineWidth(2);
    linePDGWidth->Draw();
    legall->AddEntry(linePDGWidth, "PDG Width", "l");
    legall->Draw();
    for (int imult = 0; imult < multLoopEnd; imult++)
    {
        hWidth[imult]->SetMarkerStyle(markers[imult]);
        hWidth[imult]->SetMarkerSize(1.2);
        hWidth[imult]->Draw("pe same PLC PMC");
    }
    legall->Draw();
    cWidth->SaveAs(outputfolder + "/width_all_mult.png");

    if (makePIDplots)
    {
        string QaPath = genpath.substr(0, genpath.length() - 9);
        cout << "QaPath: " << QaPath << endl;
        TString output_QA_folder = data_path + "/QA/MC";
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

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.08, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);
}

//******************************************In Data**************************************************
//******************************************In Data**************************************************
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
// string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/IR_study/459845/kstarqa/hInvMass"; // 2022 data (500 kHz)

//*******************Occupancy cut study********************************
// string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/Occupancy_effect/466154/kstarqa_OCC500_id33931/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
// string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/Occupancy_effect/463286/kstarqa/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
// string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/Occupancy_effect/463171/kstarqa/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

////*************************PID Variations for Kaon (without MID, multcentTable)**************************
// string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/999999/kstarqa/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
// string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/480448/kstarqa/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
// string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/480358/kstarqa/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

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
// string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/checks/473185/kstarqa_MID_id33593/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

//*************************ItsTpcTracksCheck, |y|<0.5******************************
// string data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/481941/kstarqa_PIDKa1_itstpc/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0

//***********************************************In MC**************************************************
//***********************************************In MC**************************************************

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
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/460571.root", "READ"); // 2023 data (135 kHz)
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/IR_study/459909.root", "READ"); // 2022 data (500 kHz)

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

//*******************************MID Cuts (other reconstructed process function)********************** */
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/479098.root", "READ"); // 2022 MC
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/478491.root", "READ"); // 2023 MC
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/479101.root", "READ"); // 2024 MC

////*************************PID Variations for Kaon (without MID, multcentTable)**************************
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/999999.root", "READ"); // 2022 MC
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/480567.root", "READ"); // 2023 MC
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/480455.root", "READ"); // 2024 MC

////*************************ItsTpcTracksCheck, betacutTOF**********************************
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/482058.root", "READ"); // 2023 MC, No effect is seen

// ******************TPC tune on data without NN ***********************************
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/459909.root", "READ"); // 2022 MC
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/460571.root", "READ"); // 2023 MC
// TFile *fileeff = new TFile("/home/sawan/check_k892/mc/LHC24f3c/459912.root", "READ"); // 2024 MC
