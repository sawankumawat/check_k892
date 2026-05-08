#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

using namespace std;
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void efficiency()
{
    bool makePIDplots = true; // qa plots
    bool skipEfficiencyPlots = false;
    string outputtype = "png"; // pdf, eps
    bool isINEL = false;
    const string kResBkg = "MIX";
    // const string kResBkg = "LIKE";
    // const string kResBkg = "ROTATED";

    gStyle->SetOptFit(1111);
    gStyle->SetPalette(kRainBow);
    // gStyle->SetOptStat(1110);
    gStyle->SetOptStat(0);

    int colors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

    // ****************Data files ********************
    string common_data_path = "/home/sawan/check_k892/output/kstar/LHC22o_pass7/";
    string common_MC_path = "/home/sawan/check_k892/mc/LHC24f3c/";

    //================================After SQM=======================
    //============2023 data===========
    // string data_path = "660453/kstarqa_id35679/hInvMass"; // 2023 (pTDepPID, pTDepPID, pTDepPIDTOF)
    // string data_path = "670168/kstarqa/hInvMass"; // 2023 (Base, INEL)

    //===========Other checks==============
    // string data_path = "655628/kstarqa/hInvMass"; // 2024
    // string data_path = "661905/kstarqa/hInvMass"; // 2024 (Event_time dependency check in TOF)
    // string data_path = "658307/kstarqa/hInvMass"; // 2024 (LoossePID, pTDepPID, pTDepPIDTOF)
    // string data_path = "MC_closure/657468/kstarqa_MC_closure/hInvMass"; (MC closure)
    // string data_path = "668039/kstarqa/hInvMass"; // 2024 (With square PID: Base, MID, MIDptDep2)
    string data_path = "668605/kstarqa_TOF3_withoutSquareCut/hInvMass"; // 2024 (Base, MID, MID_small, MID_verySmall, MIDptDep, MIDptDep2, MIDptDep2_small, MIDptDep2_verySmall (With Square PID): TOF3, TOF3_withoutSquareCut)
    // string data_path = "672297/kstarqa_MIDptDep2_verySmall_TOF3/hInvMass"; // (Base (2sigma TOF), MIDptDep2_0p3_TOF3, MIDptDep2_small_TOF3, MIDptDep2_TOF3, MIDptDep2_verySmall_TOF3)
    // string data_path = "MC_closure/669655/kstarqa_MC_closure/hInvMass"; // (MC_closure)

    //===========2024 data==============
    // string data_path = "658306/kstarqa_INELgt0/hInvMass"; // 2024 (with DeepAngle, PV Contributor, INELgt0)
    // string data_path = "660943/kstarqa/hInvMass"; // 2024 (Base, INEL, MID, MIDptDep, pTDepPID, pTDepPIDTOF, LoosePID)
    // string data_path = "663738/kstarqa/hInvMass"; // 2024 (Base, OnlyTPC)
    // string data_path = "664559/kstarqa/hInvMass"; // 2024 (Base, INEL, TOFshift, TOFshiftMID)

    data_path = common_data_path + data_path + ((kResBkg == "MIX") ? "" : "/" + kResBkg);
    TString outputfolder = data_path + (isINEL ? "/efficiency/INEL" : "/efficiency");
    gSystem->mkdir(outputfolder, kTRUE);

    //********MC files *****************
    //================================After SQM=======================
    // string MCpath = "660483.root"; // temporary (pp reference 5.36 TeV)

    //==========2023 MC===========
    // string MCpath = "660557.root"; // 2023 MC (QA+PID checks, _id52344) (with TOF_OverrideFT0 wagon, not required for 2023 dataset)
    // string MCpath = "661782.root"; // 2023 MC (QA+PID checks, _NoOverrideFT0, with same runs as data (661782_selected.root, but with very low statistics))
    // string MCpath = "662467.root"; // 2023 MC (Without RCT)
    // string MCpath = "670292.root"; // 2023 MC (LHC26d1: Base, INEL)
    // string MCpath = "670198.root"; // 2023 MC (LHC23k4i: Base, INEL)

    //==========Other checks===========
    // string MCpath = "655737.root"; // 2024 MC (without TOF_OverrideFT0 wagon)
    // string MCpath = "657468.root"; // 2024 MC (with TOF_OverrideFT0 wagon)
    // string MCpath = "662038.root"; // 2024 MC (dependency check with CovMatrix and GlobalCounters)
    // string MCpath = "661190.root"; // 2024 MC (QA+PID, MID, MIDptDep)
    // string MCpath = "662971.root"; // 2024 MC (Base, MC_closure)
    // string MCpath = "663902.root"; // 2024 MC (Base, Without NN)
    // string MCpath = "669655.root"; // 2024 MC (MC_closure)
    // string MCpath = "670822.root"; // 2024 MC (Base, MID, MID_small, MID_verySmall, MIDptDep2, MIDptDep2_small, MIDptDep2_verySmall)
    string MCpath = "669641.root"; // 2024 MC (All square cut: Base, MID, MIDptDep2, TOF3, TOF3_withoutSquareCut)
    // string MCpath = "672420.root"; // 2024 MC (MIDptDep2_0p3_TOF3, MIDptDep2_small_TOF3, MIDptDep2_TOF3, MIDptDep2_verySmall_TOF3)
    // string MCpath = "673285.root"; // 2024 MC (TOF3_WithoutShift)

    //==========2024 MC===========
    // string MCpath = "658377.root"; // 2024 MC (QA+PID checks with INEL)
    // string MCpath = "663903.root"; // 2024 MC (Base, OnlyTPC)
    // string MCpath = "664785.root"; // 2024 MC (Base, TOFshift, TOFshiftMID)
    // string MCpath = "665348.root"; // 2024 MC (All are TOFshifted, INEL, MC_closure, corrections for INEL)
    // string MCpath = "665524.root"; // 2024 MC (ToF shift with PID also shifted (base), INEL)
    // string MCpath = "667875.root"; // 2024 MC (higher TOF shift for pions, Base, MID, MIDptDep, MIDptDep2)

    TFile *fileraw = (isINEL) ? new TFile((data_path + "/yield_INEL.root").c_str(), "READ") : new TFile((data_path + "/yield.root").c_str(), "READ"); // datafile
    MCpath = common_MC_path + MCpath;
    TFile *fileeff = new TFile(MCpath.c_str(), "READ"); // MC file

    if (fileeff->IsZombie() || fileraw->IsZombie())
    {
        cout << "Error opening files in the path " << data_path << " or " << MCpath << endl;
        return;
    }
    const string genpath = "kstarqa_TOF3_withoutSquareCut/hInvMass";
    const string recpath = "kstarqa_TOF3_withoutSquareCut/hInvMass";

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    int nmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1; // number of multiplicity bins

    if (!skipEfficiencyPlots)
    {

        //===============================Efficiency histograms=======================
        // THnSparseF *hSpraseGen = (THnSparseF *)fileeff->Get(Form("%s/hk892GenpT2", genpath.c_str()));
        // THnSparseF *hSparseRec = (THnSparseF *)fileeff->Get(Form("%s/h2KstarRecpt1", recpath.c_str()));
        THnSparseF *hSpraseGen = (THnSparseF *)fileeff->Get(Form("%s/hk892GenpTCalib1", genpath.c_str()));
        THnSparseF *hSparseRec = (THnSparseF *)fileeff->Get(Form("%s/h2KstarRecptCalib1", recpath.c_str()));
        if (hSpraseGen == nullptr || hSparseRec == nullptr)
        {
            cout << "Error reading efficiency histograms " << Form("%s/hk892GenpTCalib1", genpath.c_str()) << endl;
            return;
        }

        //==================================Event loss histograms===================
        // TH1F *hAllGenColl = (TH1F *)fileeff->Get(Form("%s/hAllGenCollisions", genpath.c_str()));
        // TH1F *hAllGenColl1Rec = (TH1F *)fileeff->Get(Form("%s/hAllGenCollisions1Rec", recpath.c_str()));
        TH1F *hAllGenColl = (TH1F *)fileeff->Get(Form("%s/MCcorrections/MultiplicityGen", genpath.c_str()));
        TH1F *hAllGenColl1Rec = (TH1F *)fileeff->Get(Form("%s/MCcorrections/MultiplicityRec", recpath.c_str())); // without event splitting correction

        if (hAllGenColl == nullptr || hAllGenColl1Rec == nullptr)
        {
            cout << "Error reading event loss histograms" << endl;
            return;
        }

        //==================================Signal loss histograms===================
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

        //// Event splitting (Gen. events with atleast 1 reconstruction / All reconstructed events)
        TH1F *hRecAll = (TH1F *)fileeff->Get(Form("%s/h1RecMult2", recpath.c_str()));  // Denominator (calibrated multiplicity)
        TH1F *hGen1Rec = (TH1F *)fileeff->Get(Form("%s/h1GenMult2", recpath.c_str())); // Numerator (calibrated multiplicity)
        // TH1F *hGen1Rec = (TH1F *)fileeff->Get(Form("%s/MCcorrections/MultiplicityRec", recpath.c_str())); // Numerator

        // TH1F *hRecAll = (TH1F *)fileeff->Get(Form("%s/h1RecMult", recpath.c_str())); // Denominator
        // TH1F *hGen1Rec = (TH1F *)fileeff->Get(Form("%s/hAllGenCollisions1Rec", recpath.c_str())); // Numerator

        vector<double> nChParticlesFromMult;
        int multLoopEnd = (isINEL) ? 1 : nmultbins + 1;
        // int multLoopEnd = 1;
        int multlow, multhigh;

        TFile *spectra = (isINEL) ? new TFile((data_path + "/corrected_spectra_INEL.root").c_str(), "RECREATE") : new TFile((data_path + "/corrected_spectra.root").c_str(), "RECREATE");
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
        TH1F *hEventSplit = new TH1F("hEventSplit", "Event Splitting", nmultbins, mult_classes);

        for (int imult = 0; imult < multLoopEnd; imult++)
        {
            if (imult == 0)
            {
                multlow = 0;
                multhigh = (isINEL) ? 120 : 100; // for all multiplicity
            }
            else
            {
                multlow = mult_classes[imult - 1];
                multhigh = mult_classes[imult];
            }
            cout << "Mult low is " << multlow << " and mult high is " << multhigh << endl;
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
                cout << "Error reading yield histogram, " << Form("mult_%d-%d/yield_bincount", multlow, multhigh) << " not found" << endl;
                return;
            }
            // heff[imult] = (TH1F *)hyieldBinCount->Clone(); // Just taking bins from the yield histogram which will be set to efficiency value.
            heff[imult] = (TH1F *)hyieldIntegral->Clone(); // Just taking bins from the yield histogram which will be set to efficiency value.
            if (hyieldIntegral->GetNbinsX() != Npt)
            {
                cerr << "The total number of bins in the yield histogram (data) and the fitting_range.h file do not match" << endl;
                return;
            }

            //=============Event loss calculations=====================
            heventloss[imult] = new TH1F(Form("hEventLoss_%d", imult), "Event Loss", Npt, pT_bins);

            double eventLossNum = hAllGenColl1Rec->Integral(hAllGenColl1Rec->GetXaxis()->FindBin(multlow + 0.001), hAllGenColl1Rec->GetXaxis()->FindBin(multhigh - 0.001));

            double eventLossDen = hAllGenColl->Integral(hAllGenColl->GetXaxis()->FindBin(multlow + 0.001), hAllGenColl->GetXaxis()->FindBin(multhigh - 0.001));

            // TH1F *hNchInMult = (TH1F *)hNchVsMultGen->ProjectionY(Form("hNchInMult_%d", imult), hNchVsMultGen->GetXaxis()->FindBin(multlow + 0.01), hNchVsMultGen->GetXaxis()->FindBin(multhigh - 0.01));
            // hAllGenColl1Rec->Multiply(hAllGenColl1Rec, hNchInMult);
            // hAllGenColl->Multiply(hAllGenColl, hNchInMult);
            // double eventLossNum = hAllGenColl1Rec->Integral();
            // double eventLossDen = hAllGenColl->Integral();
            // cout << "Multiplicity bin " << multlow << "-" << multhigh << ", Event Loss " << eventLossNum / eventLossDen << endl;

            //// Event splitting calculation
            double eventSplitDen = hRecAll->Integral(hRecAll->GetXaxis()->FindBin(multlow + 0.001), hRecAll->GetXaxis()->FindBin(multhigh - 0.001));
            double eventSplitNum = hGen1Rec->Integral(hGen1Rec->GetXaxis()->FindBin(multlow + 0.001), hGen1Rec->GetXaxis()->FindBin(multhigh - 0.001));
            double eventSplitting = eventSplitNum / eventSplitDen;
            if (imult > 0)
            {
                hEventSplit->SetBinContent(imult, eventSplitting);
                hEventSplit->SetBinError(imult, 0.00000001);
                cout << "Multiplicity bin " << multlow << "-" << multhigh << ", Event Splitting " << eventSplitting << endl;
            }
            else if (imult == 0)
            {
                for (int ibinstemp = 0; ibinstemp <= hEventSplit->GetNbinsX(); ibinstemp++)
                {
                    hEventSplit->SetBinContent(ibinstemp, eventSplitting);
                    hEventSplit->SetBinError(ibinstemp, 0.00000001);
                }
            }

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
                cout << "Multiplicity bin " << multlow << "-" << multhigh << ", pT bin " << lowpt << "-" << highpt << ", Efficiency " << efficiency << endl;
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
                    // cout << "signal loss value is " << signalLoss << endl;
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
            TH1F *hlossCorrected_bincount = (TH1F *)hyieldBinCount->Clone("hlossCorrected_bincount");

            hlossCorrected_integral->Multiply(hRatioEvBySig[imult]); // Event/Signal loss correction
            hlossCorrected_bincount->Multiply(hRatioEvBySig[imult]);
            hlossCorrected_integral->Scale(1.0 / eventSplitting); // Event splitting correction
            hlossCorrected_bincount->Scale(1.0 / eventSplitting);

            TDirectory *dir = spectra->mkdir(Form("mult_%d-%d", multlow, multhigh));
            dir->cd();
            heff[imult]->Write("heff");
            hyieldBinCount->Write("corrected_spectra_BinCount");
            hyieldIntegral->Write("corrected_spectra_Integral");
            hlossCorrected_integral->Write("corrected_spectra_Integral_final");
            hlossCorrected_bincount->Write("corrected_spectra_BinCount_final");
            spectra->cd();
        }

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

        //// Plot other plots for all multiplicity bins
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
        (multLoopEnd == 1) ? heff[0]->Fit(pol1fit, "RQ0") : heff[4]->Fit(pol1fit, "RQ0");
        cout << "Intercept value : " << pol1fit->GetParameter(0) << endl;
        TLine *line = new TLine(0, pol1fit->GetParameter(0), Npt, pol1fit->GetParameter(0));
        line->SetLineStyle(2);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->Draw();
        TLegend *legall = new TLegend(0.20, 0.8, 0.92, 0.92);
        legall->SetTextSize(0.03);
        legall->SetNColumns(5);
        legall->SetFillStyle(0);
        legall->SetBorderSize(0);
        legall->AddEntry(heff[0], "0-100%", "p");
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
            SetHistoQA(hSignalLoss[imult]);
            hSignalLoss[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hSignalLoss[imult]->GetYaxis()->SetTitle("Event and signal loss");
            hSignalLoss[imult]->GetYaxis()->SetTitleOffset(1.6);
            hSignalLoss[imult]->SetStats(0);
            hSignalLoss[imult]->SetMaximum(1.29);
            hSignalLoss[imult]->SetMinimum(0.3);
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
            SetHistoQA(hRatioEvBySig[imult]);
            hRatioEvBySig[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hRatioEvBySig[imult]->GetYaxis()->SetTitle("Event Loss / Signal Loss");
            hRatioEvBySig[imult]->GetYaxis()->SetTitleOffset(1.6);
            hRatioEvBySig[imult]->SetStats(0);
            hRatioEvBySig[imult]->SetMaximum(1.25);
            hRatioEvBySig[imult]->SetMinimum(0.6);
            hRatioEvBySig[imult]->SetMarkerStyle(markers[imult]);
            hRatioEvBySig[imult]->SetMarkerSize(1.2);
            hRatioEvBySig[imult]->SetMarkerColor(colors[imult]);
            hRatioEvBySig[imult]->SetLineColor(colors[imult]);
            hRatioEvBySig[imult]->Draw("pe same");
            legall->AddEntry(hRatioEvBySig[imult], Form("%d-%d%%", multlow, multhigh), "p");
        }
        legall->Draw();
        cEventBySignalLoss->SaveAs(outputfolder + "/event_by_signal_loss.png");

        TCanvas *cEventSplit = new TCanvas("", "", 720, 720);
        SetCanvasStyle(cEventSplit, 0.16, 0.06, 0.01, 0.14);
        SetHistoQA(hEventSplit);
        hEventSplit->GetXaxis()->SetTitle("Multiplicity (%)");
        hEventSplit->GetYaxis()->SetTitle("#epsilon_{Event_split}");
        hEventSplit->GetYaxis()->SetTitleOffset(1.7);
        hEventSplit->SetStats(0);
        hEventSplit->SetMarkerStyle(20);
        hEventSplit->SetMarkerSize(1.2);
        hEventSplit->SetMaximum(1.0037);
        hEventSplit->SetMinimum(0.993);
        hEventSplit->Draw("pe");
        if (multLoopEnd == 1)
        {
            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.04);
            lat.DrawLatex(0.2, 0.8, Form("Multiplicity: %.0f-%.0f%%", mult_classes[0], 100.0));
        }
        cEventSplit->SaveAs(outputfolder + "/event_splitting.png");

        TCanvas *cSignificance = new TCanvas("", "", 720, 720);
        SetCanvasStyle(cSignificance, 0.16, 0.06, 0.06, 0.14);
        // hSignificance[0]->Draw("p");
        legall->Clear();
        for (int imult = 0; imult < multLoopEnd; imult++)
        {
            SetHistoQA(hSignificance[imult]);
            hSignificance[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hSignificance[imult]->GetYaxis()->SetTitle("Significance");
            hSignificance[imult]->GetYaxis()->SetTitleOffset(1.6);
            hSignificance[imult]->SetMaximum(1250);
            hSignificance[imult]->SetMinimum(-5);
            hSignificance[imult]->SetMarkerStyle(markers[imult]);
            hSignificance[imult]->SetMarkerSize(1.2);
            hSignificance[imult]->Draw("p same PLC PMC");
            multlow = mult_classes[imult - 1];
            multhigh = mult_classes[imult];
            legall->AddEntry(hSignificance[imult], Form("%d-%d%%", multlow, multhigh), "p");
        }
        legall->Draw();
        cSignificance->SaveAs(outputfolder + "/significance_all_mult.png");

        TCanvas *cChi2byNDF = new TCanvas("", "", 720, 720);
        SetCanvasStyle(cChi2byNDF, 0.16, 0.06, 0.01, 0.14);
        // hChi2byNDF[0]->Draw("p");
        for (int imult = 0; imult < multLoopEnd; imult++)
        {
            SetHistoQA(hChi2byNDF[imult]);
            hChi2byNDF[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hChi2byNDF[imult]->GetYaxis()->SetTitle("#chi^{2}/NDF");
            hChi2byNDF[imult]->GetYaxis()->SetTitleOffset(1.6);
            hChi2byNDF[imult]->SetMaximum(6.5);
            hChi2byNDF[imult]->SetMinimum(0);
            hChi2byNDF[imult]->SetStats(0);
            hChi2byNDF[imult]->SetMarkerStyle(markers[imult]);
            hChi2byNDF[imult]->SetMarkerSize(1.2);
            hChi2byNDF[imult]->Draw("p same PLC PMC");
        }
        legall->Draw();
        cChi2byNDF->SaveAs(outputfolder + "/chi2byNDF_all_mult.png");

        TCanvas *cMass = new TCanvas("", "", 720, 720);
        SetCanvasStyle(cMass, 0.16, 0.06, 0.01, 0.14);
        for (int imult = 0; imult < multLoopEnd; imult++)
        {
            SetHistoQA(hMass[imult]);
            hMass[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hMass[imult]->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
            hMass[imult]->GetYaxis()->SetTitleOffset(1.6);
            hMass[imult]->GetYaxis()->SetRangeUser(0.878, 0.919);
            hMass[imult]->SetStats(0);
            hMass[imult]->SetMarkerStyle(markers[imult]);
            hMass[imult]->SetMarkerSize(1.2);
            hMass[imult]->Draw("pe same PLC PMC");
        }
        TLine *linePDG = new TLine(0, 0.895, 20, 0.895);
        linePDG->SetLineStyle(2);
        linePDG->SetLineColor(2);
        linePDG->SetLineWidth(2);
        linePDG->Draw();
        legall->AddEntry(linePDG, "PDG Mass", "l");
        legall->Draw();
        cMass->SaveAs(outputfolder + "/mass_all_mult.png");

        TCanvas *cWidth = new TCanvas("", "", 720, 720);
        SetCanvasStyle(cWidth, 0.16, 0.06, 0.01, 0.14);
        TLine *linePDGWidth = new TLine(0, widthpdg, 20, widthpdg);
        linePDGWidth->SetLineStyle(2);
        linePDGWidth->SetLineColor(2);
        linePDGWidth->SetLineWidth(2);
        linePDGWidth->Draw();
        legall->AddEntry(linePDGWidth, "PDG Width", "l");
        legall->Draw();
        for (int imult = 0; imult < multLoopEnd; imult++)
        {
            SetHistoQA(hWidth[imult]);
            hWidth[imult]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hWidth[imult]->GetYaxis()->SetTitle("Width (GeV/#it{c}^{2})");
            hWidth[imult]->GetYaxis()->SetTitleOffset(1.6);
            hWidth[imult]->GetYaxis()->SetRangeUser(0.047 - 0.04, 0.047 + 0.05);
            hWidth[imult]->SetStats(0);
            hWidth[imult]->SetMarkerStyle(markers[imult]);
            hWidth[imult]->SetMarkerSize(1.2);
            hWidth[imult]->Draw("pe same PLC PMC");
        }
        legall->Draw();
        cWidth->SaveAs(outputfolder + "/width_all_mult.png");
    }

    if (makePIDplots)
    {
        string QaPath = genpath.substr(0, genpath.length() - 9);
        cout << "QaPath: " << QaPath << endl;
        TString output_QA_folder = data_path + "/QA/MC";
        TString output_QA_folder_before = data_path + "/QA/MC/before";
        gSystem->mkdir(output_QA_folder, kTRUE);
        gSystem->mkdir(output_QA_folder_before, kTRUE);

        TH3F *hDCAxy3D = (TH3F *)fileeff->Get(Form("%s/eventSelection/hDcaxy_cent_pt", QaPath.c_str()));
        TH3F *hDCAz3D = (TH3F *)fileeff->Get(Form("%s/eventSelection/hDcaz_cent_pt", QaPath.c_str()));
        TH1F *hMult = (TH1F *)fileeff->Get(Form("%s/eventSelection/hMultiplicity", QaPath.c_str()));
        TH1F *hvz = (TH1F *)fileeff->Get(Form("%s/eventSelection/hVertexZRec", QaPath.c_str()));
        TH1F *hoccupancy = (TH1F *)fileeff->Get(Form("%s/eventSelection/hOccupancy", QaPath.c_str()));
        TH1F *hEventCut = (TH1F *)fileeff->Get(Form("%s/eventSelection/hEventCut", QaPath.c_str()));
        TH1F *htracksRecMC = (TH1F *)fileeff->Get(Form("%s/eventSelection/recMCparticles", QaPath.c_str()));
        TH1D *hDCAxy = hDCAxy3D->ProjectionX("hDCAxy", -1, -1, -1, -1);
        TH1D *hDCAz = hDCAz3D->ProjectionX("hDCAz", -1, -1, -1, -1);
        if (hDCAxy == nullptr || hDCAz == nullptr || hMult == nullptr || hvz == nullptr || hoccupancy == nullptr || hEventCut == nullptr || htracksRecMC == nullptr)
        {
            cerr << "Event selection histograms not found!!!!!!!!!!!!" << endl;
            return;
        }

        TH3F *hNsigmaTPCTOFKaon = (TH3F *)fileeff->Get(Form("%s/hPID/After/hNsigma_TPC_TOF_Ka_after", QaPath.c_str()));
        TH3F *hNsigmaTPCTOFPion = (TH3F *)fileeff->Get(Form("%s/hPID/After/hNsigma_TPC_TOF_Pi_after", QaPath.c_str()));
        TH3F *hNsigmaTPCKaon = (TH3F *)fileeff->Get(Form("%s/hPID/After/hTPCnsigKa_mult_pt", QaPath.c_str()));
        TH3F *hNsigmaTPCPion = (TH3F *)fileeff->Get(Form("%s/hPID/After/hTPCnsigPi_mult_pt", QaPath.c_str()));
        TH3F *hNsigmaTOFKaon = (TH3F *)fileeff->Get(Form("%s/hPID/After/hTOFnsigKa_mult_pt", QaPath.c_str()));
        TH3F *hNsigmaTOFPion = (TH3F *)fileeff->Get(Form("%s/hPID/After/hTOFnsigPi_mult_pt", QaPath.c_str()));
        if (hNsigmaTPCTOFKaon == nullptr || hNsigmaTPCTOFPion == nullptr || hNsigmaTPCKaon == nullptr || hNsigmaTPCPion == nullptr || hNsigmaTOFKaon == nullptr || hNsigmaTOFPion == nullptr)
        {
            cerr << "PID histograms after selection not found!!!!!!!!!!!!" << endl;
            return;
        }

        TH3F *hNsigmaTPCTOFKaon_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hNsigma_TPC_TOF_Ka_before", QaPath.c_str()));
        TH3F *hNsigmaTPCTOFPion_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hNsigma_TPC_TOF_Pi_before", QaPath.c_str()));
        TH3F *hNsigmaTPCKaon_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hTPCnsigKa_mult_pt", QaPath.c_str()));
        TH3F *hNsigmaTPCPion_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hTPCnsigPi_mult_pt", QaPath.c_str()));
        TH3F *hNsigmaTOFKaon_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hTOFnsigKa_mult_pt", QaPath.c_str()));
        TH3F *hNsigmaTOFPion_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hTOFnsigPi_mult_pt", QaPath.c_str()));
        if (hNsigmaTPCTOFKaon_before == nullptr || hNsigmaTPCTOFPion_before == nullptr || hNsigmaTPCKaon_before == nullptr || hNsigmaTPCPion_before == nullptr || hNsigmaTOFKaon_before == nullptr || hNsigmaTOFPion_before == nullptr)
        {
            cerr << "PID histograms before selection not found!!!!!!!!!!!!" << endl;
            return;
        }

        outputtype = "png";
        TCanvas *cMult = new TCanvas("cMult", "Multiplicity", 720, 720);
        SetCanvasStyle(cMult, 0.14, 0.03, 0.06, 0.14);
        SetHistoQA(hMult);
        // hMult->GetXaxis()->SetTitle("Multiplicity (%)");
        hMult->GetXaxis()->SetTitle("Multiplicity (%)");
        hMult->GetYaxis()->SetTitle("Counts");
        hMult->Draw();
        cMult->SaveAs(output_QA_folder + ("/Multiplicity." + outputtype).c_str());

        TCanvas *cvz = new TCanvas("cvz", "Vertex Z", 720, 720);
        SetCanvasStyle(cvz, 0.14, 0.03, 0.06, 0.14);
        SetHistoQA(hvz);
        hvz->GetXaxis()->SetTitle("Vertex Z (cm)");
        hvz->GetYaxis()->SetTitle("Counts");
        hvz->Draw();
        cvz->SaveAs(output_QA_folder + ("/VertexZ." + outputtype).c_str());

        TCanvas *cTracksData = new TCanvas("cTracksData", "Tracks Data", 1440, 720);
        SetCanvasStyle(cTracksData, 0.1, 0.06, 0.06, 0.16);
        SetHistoQA(htracksRecMC);
        // gPad->SetLogy(1);
        gPad->SetGrid(1, 0);
        gPad->SetLogy(1);
        htracksRecMC->SetTitle("Tracks Rec MC");
        htracksRecMC->GetYaxis()->SetTitle("Counts");
        htracksRecMC->GetYaxis()->SetTitleOffset(0.8);
        htracksRecMC->GetXaxis()->SetRangeUser(0, 17);
        // htracksRecMC->SetMaximum(htracksRecMC->GetMaximum() * 10);
        htracksRecMC->Draw();
        cTracksData->SaveAs(output_QA_folder + ("/TracksRecMC." + outputtype).c_str());

        TCanvas *cEventCut = new TCanvas("cEventCut", "Event Cut", 1080, 720);
        SetCanvasStyle(cEventCut, 0.1, 0.05, 0.06, 0.17);
        SetHistoQA(hEventCut);
        gPad->SetGrid(1, 0);
        // gPad->SetLogy(1);
        hEventCut->SetTitle("Event selections");
        hEventCut->GetYaxis()->SetTitle("Counts");
        hEventCut->GetYaxis()->SetTitleOffset(0.8);
        hEventCut->GetXaxis()->SetRangeUser(0, 11);
        hEventCut->Draw();
        cEventCut->SaveAs(output_QA_folder + ("/EventCut." + outputtype).c_str());

        //===============PID plots==========================
        // Plots to make. 2D plots TPCKa, TPCPi, TOFKa, TOFPi as a function of pT. Then 2D plot of TPC vs TOF of Ka and Pi at different pT ranges.

        double pTrangesForTPCandTOF[] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
        vector<vector<int>> multRangesForTPCandTOF = {{0, 1}, {1, 5}, {10, 20}, {30, 40}, {50, 70}, {70, 100}};
        TH2F *h2DTPCKaon = (TH2F *)hNsigmaTPCKaon->Project3D("xz");
        TH2F *h2DTPCPion = (TH2F *)hNsigmaTPCPion->Project3D("xz");
        TH2F *h2DTOFKaon = (TH2F *)hNsigmaTOFKaon->Project3D("xz");
        TH2F *h2DTOFPion = (TH2F *)hNsigmaTOFPion->Project3D("xz");
        TH2 *h2DNsigmaTPCTOFKaon[6];
        TH2 *h2DNsigmaTPCTOFPion[6];
        TH1F *h1DNsigmaTPCKaon[6];
        TH1F *h1DNsigmaTPCPion[6];
        TH1F *h1DNsigmaTOFKaon[6];
        TH1F *h1DNsigmaTOFPion[6];

        TH1F *h1DNsigmaTPCKaon_before[6];
        TH1F *h1DNsigmaTPCPion_before[6];
        TH1F *h1DNsigmaTOFKaon_before[6];
        TH1F *h1DNsigmaTOFPion_before[6];

        TH1F *h1DNsigmaTPCKaon_pt[6];
        TH1F *h1DNsigmaTPCPion_pt[6];
        TH1F *h1DNsigmaTOFKaon_pt[6];
        TH1F *h1DNsigmaTOFPion_pt[6];

        for (int i = 0; i < 6; i++)
        {
            int lowBinpT_z = hNsigmaTPCTOFKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i] + 0.001);
            int highBinpT_z = hNsigmaTPCTOFKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i + 1] - 0.001);

            hNsigmaTPCTOFKaon->GetZaxis()->SetRange(lowBinpT_z, highBinpT_z);
            hNsigmaTPCTOFPion->GetZaxis()->SetRange(lowBinpT_z, highBinpT_z);

            h2DNsigmaTPCTOFKaon[i] = (TH2 *)hNsigmaTPCTOFKaon->Project3D("yx"); // the axis are switched during projection, so we get TPC on y and TOF on x
            h2DNsigmaTPCTOFPion[i] = (TH2 *)hNsigmaTPCTOFPion->Project3D("yx");

            h2DNsigmaTPCTOFKaon[i]->SetName(Form("xy_Ka_%d", i));
            h2DNsigmaTPCTOFPion[i]->SetName(Form("xy_Pi_%d", i));

            // Now lets make 1D projections with limit on multiplicity ranges
            int lowBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multRangesForTPCandTOF[i][0] + 0.001);
            int highBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multRangesForTPCandTOF[i][1] - 0.001);

            h1DNsigmaTPCKaon[i] = (TH1F *)hNsigmaTPCKaon->ProjectionX(Form("h1D_TPC_Ka_%d", i), lowBinMult, highBinMult, -1, -1);
            h1DNsigmaTPCPion[i] = (TH1F *)hNsigmaTPCPion->ProjectionX(Form("h1D_TPC_Pi_%d", i), lowBinMult, highBinMult, -1, -1);
            h1DNsigmaTOFKaon[i] = (TH1F *)hNsigmaTOFKaon->ProjectionX(Form("h1D_TOF_Ka_%d", i), lowBinMult, highBinMult, -1, -1);
            h1DNsigmaTOFPion[i] = (TH1F *)hNsigmaTOFPion->ProjectionX(Form("h1D_TOF_Pi_%d", i), lowBinMult, highBinMult, -1, -1);

            h1DNsigmaTPCKaon_before[i] = (TH1F *)hNsigmaTPCKaon_before->ProjectionX(Form("h1D_TPC_Ka_before_%d", i), lowBinMult, highBinMult, -1, -1);
            h1DNsigmaTPCPion_before[i] = (TH1F *)hNsigmaTPCPion_before->ProjectionX(Form("h1D_TPC_Pi_before_%d", i), lowBinMult, highBinMult, -1, -1);
            h1DNsigmaTOFKaon_before[i] = (TH1F *)hNsigmaTOFKaon_before->ProjectionX(Form("h1D_TOF_Ka_before_%d", i), lowBinMult, highBinMult, -1, -1);
            h1DNsigmaTOFPion_before[i] = (TH1F *)hNsigmaTOFPion_before->ProjectionX(Form("h1D_TOF_Pi_before_%d", i), lowBinMult, highBinMult, -1, -1);

            int lowBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i] + 0.001);
            int highBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i + 1] - 0.001);

            h1DNsigmaTPCKaon_pt[i] = (TH1F *)hNsigmaTPCKaon->ProjectionX(Form("h1D_TPC_Ka_pt_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
            h1DNsigmaTPCPion_pt[i] = (TH1F *)hNsigmaTPCPion->ProjectionX(Form("h1D_TPC_Pi_pt_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
            h1DNsigmaTOFKaon_pt[i] = (TH1F *)hNsigmaTOFKaon->ProjectionX(Form("h1D_TOF_Ka_pt_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
            h1DNsigmaTOFPion_pt[i] = (TH1F *)hNsigmaTOFPion->ProjectionX(Form("h1D_TOF_Pi_pt_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        }

        TCanvas *cNsigmaTPCKaon = new TCanvas("cNsigmaTPCKaon", "Nsigma TPC Kaon", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(h2DTPCKaon);
        // gPad->SetLogz(1);
        h2DTPCKaon->GetYaxis()->SetTitle("n#sigma_{TPC}");
        h2DTPCKaon->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h2DTPCKaon->GetYaxis()->SetRangeUser(-3, 3);
        h2DTPCKaon->GetXaxis()->SetRangeUser(0, 5);
        h2DTPCKaon->Draw("colz");
        cNsigmaTPCKaon->SaveAs(output_QA_folder + ("/NsigmaTPCKaon2D." + outputtype).c_str());

        TCanvas *cNsigmaTPCPion = new TCanvas("cNsigmaTPCPion", "Nsigma TPC Pion", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(h2DTPCPion);
        // gPad->SetLogz(1);
        h2DTPCPion->GetYaxis()->SetTitle("n#sigma_{TPC}");
        h2DTPCPion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h2DTPCPion->GetYaxis()->SetRangeUser(-3, 3);
        h2DTPCPion->GetXaxis()->SetRangeUser(0, 5);
        h2DTPCPion->Draw("colz");
        cNsigmaTPCPion->SaveAs(output_QA_folder + ("/NsigmaTPCPion2D." + outputtype).c_str());

        TCanvas *cNsigmaTOFKaon = new TCanvas("cNsigmaTOFKaon", "Nsigma TOF Kaon", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(h2DTOFKaon);
        // gPad->SetLogz(1);
        h2DTOFKaon->GetYaxis()->SetTitle("n#sigma_{TOF}");
        h2DTOFKaon->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h2DTOFKaon->GetYaxis()->SetRangeUser(-3, 3);
        h2DTOFKaon->GetXaxis()->SetRangeUser(0, 5);
        h2DTOFKaon->Draw("colz");
        cNsigmaTOFKaon->SaveAs(output_QA_folder + ("/NsigmaTOFKaon2D." + outputtype).c_str());

        TCanvas *cNsigmaTOFPion = new TCanvas("cNsigmaTOFPion", "Nsigma TOF Pion", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(h2DTOFPion);
        // gPad->SetLogz(1);
        h2DTOFPion->GetYaxis()->SetTitle("n#sigma_{TOF}");
        h2DTOFPion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h2DTOFPion->GetYaxis()->SetRangeUser(-3, 3);
        h2DTOFPion->GetXaxis()->SetRangeUser(0, 5);
        h2DTOFPion->Draw("colz");
        cNsigmaTOFPion->SaveAs(output_QA_folder + ("/NsigmaTOFPion2D." + outputtype).c_str());

        TCanvas *cNsigmaTPCTOFKaon = new TCanvas("cNsigmaTPCTOFKaon", "Nsigma TPC vs TOF Kaon", 1080, 720);
        TCanvas *cNsigmaTPCTOFPion = new TCanvas("cNsigmaTPCTOFPion", "Nsigma TPC vs TOF Pion", 1080, 720);
        SetCanvasStyle(cNsigmaTPCTOFKaon, 0.14, 0.15, 0.06, 0.14);
        SetCanvasStyle(cNsigmaTPCTOFPion, 0.14, 0.15, 0.06, 0.14);
        cNsigmaTPCTOFKaon->Divide(3, 2);
        cNsigmaTPCTOFPion->Divide(3, 2);

        TLatex latPID;
        latPID.SetNDC();
        latPID.SetTextFont(42);
        latPID.SetTextSize(0.06);

        for (int i = 0; i < 6; i++)
        {
            cNsigmaTPCTOFKaon->cd(i + 1);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.15);
            gPad->SetBottomMargin(0.13);
            gPad->SetTopMargin(0.06);
            SetCanvasStyle(cNsigmaTPCTOFKaon, 0.14, 0.15, 0.06, 0.14);
            SetHistoQA2D(h2DNsigmaTPCTOFKaon[i]);
            h2DNsigmaTPCTOFKaon[i]->GetXaxis()->SetTitle("n#sigma_{TPC}");
            h2DNsigmaTPCTOFKaon[i]->GetYaxis()->SetTitle("n#sigma_{TOF}");
            h2DNsigmaTPCTOFKaon[i]->GetYaxis()->SetRangeUser(-3, 3);
            h2DNsigmaTPCTOFKaon[i]->GetXaxis()->SetRangeUser(-3, 3);
            h2DNsigmaTPCTOFKaon[i]->Draw("colz");
            latPID.DrawLatex(0.2, 0.85, Form("p_{T}: %.1f-%.1f GeV/c", pTrangesForTPCandTOF[i], pTrangesForTPCandTOF[i + 1]));
            // cout<<"pT range: "<<pTrangesForTPCandTOF[i]<<"-"<<pTrangesForTPCandTOF[i + 1]<<" GeV/c, mean TPC nSigma: "<<h2DNsigmaTPCTOFKaon[i]->GetMean(1)<<", mean TOF nSigma: "<<h2DNsigmaTPCTOFKaon[i]->GetMean(2)<<endl;

            cNsigmaTPCTOFPion->cd(i + 1);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.15);
            gPad->SetBottomMargin(0.13);
            gPad->SetTopMargin(0.06);
            SetCanvasStyle(cNsigmaTPCTOFPion, 0.14, 0.15, 0.06, 0.14);
            SetHistoQA2D(h2DNsigmaTPCTOFPion[i]);
            h2DNsigmaTPCTOFPion[i]->GetXaxis()->SetTitle("n#sigma_{TPC}");
            h2DNsigmaTPCTOFPion[i]->GetYaxis()->SetTitle("n#sigma_{TOF}");
            h2DNsigmaTPCTOFPion[i]->GetYaxis()->SetRangeUser(-3, 3);
            h2DNsigmaTPCTOFPion[i]->GetXaxis()->SetRangeUser(-3, 3);
            h2DNsigmaTPCTOFPion[i]->Draw("colz");
            latPID.DrawLatex(0.2, 0.85, Form("p_{T}: %.1f-%.1f GeV/c", pTrangesForTPCandTOF[i], pTrangesForTPCandTOF[i + 1]));
        }

        cNsigmaTPCTOFKaon->SaveAs(output_QA_folder + ("/NsigmaTPCTOFKaon_2D." + outputtype).c_str());
        cNsigmaTPCTOFPion->SaveAs(output_QA_folder + ("/NsigmaTPCTOFPion_2D." + outputtype).c_str());

        TCanvas *cNsigmaTPCKaon_1D = new TCanvas("cNsigmaTPCKaon_1D", "Nsigma TPC Kaon 1D", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon_1D, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTPCKaon_1D->SetGrid(1, 1);

        TCanvas *cNsigmaTPCKaon_1D_before = new TCanvas("cNsigmaTPCKaon_1D_before", "Nsigma TPC Kaon 1D", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon_1D_before, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTPCKaon_1D_before->SetGrid(1, 1);

        TCanvas *cNsigmaTPCKaon_1D_pt = new TCanvas("cNsigmaTPCKaon_1D_pt", "Nsigma TPC Kaon 1D", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon_1D_pt, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTPCKaon_1D_pt->SetGrid(1, 1);

        TLegend *leg1DPID = new TLegend(0.17, 0.82, 0.6, 0.92);
        leg1DPID->SetNColumns(3);
        leg1DPID->SetTextSize(0.028);
        leg1DPID->SetFillStyle(0);
        leg1DPID->SetBorderSize(0);
        leg1DPID->SetTextFont(42);

        TLegend *leg1DPID_pt = new TLegend(0.17, 0.82, 0.8, 0.92);
        leg1DPID_pt->SetNColumns(3);
        leg1DPID_pt->SetTextSize(0.028);
        leg1DPID_pt->SetFillStyle(0);
        leg1DPID_pt->SetBorderSize(0);
        leg1DPID_pt->SetTextFont(42);
        leg1DPID_pt->SetHeader("#it{p}_{T} ranges");

        for (int i = 0; i < 6; i++)
        {
            cNsigmaTPCKaon_1D->cd();
            SetHistoQA(h1DNsigmaTPCKaon[i]);
            h1DNsigmaTPCKaon[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCKaon[i]->GetXaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");
            h1DNsigmaTPCKaon[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCKaon[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            h1DNsigmaTPCKaon[i]->SetMaximum(h1DNsigmaTPCKaon[i]->GetMaximum() * 15);
            h1DNsigmaTPCKaon[i]->Draw("p same");
            leg1DPID->AddEntry(h1DNsigmaTPCKaon[i], Form("%d-%d%%", multRangesForTPCandTOF[i][0], multRangesForTPCandTOF[i][1]), "p");
            leg1DPID_pt->AddEntry(h1DNsigmaTPCKaon_pt[i], Form("%.1f-%.1f", pTrangesForTPCandTOF[i], pTrangesForTPCandTOF[i + 1]), "p");

            cNsigmaTPCKaon_1D_before->cd();
            SetHistoQA(h1DNsigmaTPCKaon_before[i]);
            h1DNsigmaTPCKaon_before[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCKaon_before[i]->GetXaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");
            h1DNsigmaTPCKaon_before[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCKaon_before[i]->GetXaxis()->SetRangeUser(-5, 5);
            h1DNsigmaTPCKaon_before[i]->SetMaximum(h1DNsigmaTPCKaon_before[i]->GetMaximum() * 15);
            h1DNsigmaTPCKaon_before[i]->Draw("p same");

            cNsigmaTPCKaon_1D_pt->cd();
            SetHistoQA(h1DNsigmaTPCKaon_pt[i]);
            h1DNsigmaTPCKaon_pt[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCKaon_pt[i]->GetXaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");
            h1DNsigmaTPCKaon_pt[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCKaon_pt[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            h1DNsigmaTPCKaon_pt[i]->SetMaximum(h1DNsigmaTPCKaon_pt[i]->GetMaximum() * 15);
            h1DNsigmaTPCKaon_pt[i]->Draw("p same");
        }
        cNsigmaTPCKaon_1D->cd();
        leg1DPID->Draw();
        TLine *lineTPCKaon = new TLine(0, 0, 0, h1DNsigmaTPCKaon[0]->GetMaximum());
        lineTPCKaon->SetLineStyle(2);
        lineTPCKaon->SetLineColor(2);
        lineTPCKaon->Draw();
        cNsigmaTPCKaon_1D->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_1D." + outputtype).c_str());

        cNsigmaTPCKaon_1D_before->cd();
        leg1DPID->Draw();
        lineTPCKaon->SetY2(h1DNsigmaTPCKaon_before[0]->GetMaximum());
        lineTPCKaon->Draw();
        cNsigmaTPCKaon_1D_before->SaveAs(output_QA_folder_before + ("/NsigmaTPCKaon_1D_before." + outputtype).c_str());

        cNsigmaTPCKaon_1D_pt->cd();
        leg1DPID_pt->Draw();
        lineTPCKaon->SetY2(h1DNsigmaTPCKaon_pt[0]->GetMaximum());
        lineTPCKaon->Draw();
        cNsigmaTPCKaon_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_1D_pt." + outputtype).c_str());

        TCanvas *cNsigmaTPCPion_1D = new TCanvas("cNsigmaTPCPion_1D", "Nsigma TPC Pion 1D", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion_1D, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTPCPion_1D->SetGrid(1, 1);

        TCanvas *cNsigmaTPCPion_1D_before = new TCanvas("cNsigmaTPCPion_1D_before", "Nsigma TPC Pion 1D", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion_1D_before, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTPCPion_1D_before->SetGrid(1, 1);

        TCanvas *cNsigmaTPCPion_1D_pt = new TCanvas("cNsigmaTPCPion_1D_pt", "Nsigma TPC Pion 1D", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion_1D_pt, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTPCPion_1D_pt->SetGrid(1, 1);

        for (int i = 0; i < 6; i++)
        {
            cNsigmaTPCPion_1D->cd();
            SetHistoQA(h1DNsigmaTPCPion[i]);
            h1DNsigmaTPCPion[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCPion[i]->GetXaxis()->SetTitle("n#sigma_{TPC} #pi^{#pm}");
            h1DNsigmaTPCPion[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCPion[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            h1DNsigmaTPCPion[i]->SetMaximum(h1DNsigmaTPCPion[i]->GetMaximum() * 15);
            h1DNsigmaTPCPion[i]->Draw("p same");

            cNsigmaTPCPion_1D_before->cd();
            SetHistoQA(h1DNsigmaTPCPion_before[i]);
            h1DNsigmaTPCPion_before[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCPion_before[i]->GetXaxis()->SetTitle("n#sigma_{TPC} #pi^{#pm}");
            h1DNsigmaTPCPion_before[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCPion_before[i]->GetXaxis()->SetRangeUser(-5, 5);
            h1DNsigmaTPCPion_before[i]->SetMaximum(h1DNsigmaTPCPion_before[i]->GetMaximum() * 15);
            h1DNsigmaTPCPion_before[i]->Draw("p same");

            cNsigmaTPCPion_1D_pt->cd();
            SetHistoQA(h1DNsigmaTPCPion_pt[i]);
            h1DNsigmaTPCPion_pt[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCPion_pt[i]->GetXaxis()->SetTitle("n#sigma_{TPC} #pi^{#pm}");
            h1DNsigmaTPCPion_pt[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCPion_pt[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            h1DNsigmaTPCPion_pt[i]->SetMaximum(h1DNsigmaTPCPion_pt[i]->GetMaximum() * 12);
            h1DNsigmaTPCPion_pt[i]->Draw("p same");
        }
        cNsigmaTPCPion_1D->cd();
        leg1DPID->Draw();
        TLine *lineTPCPion = new TLine(0, 0, 0, h1DNsigmaTPCPion[0]->GetMaximum());
        lineTPCPion->SetLineStyle(2);
        lineTPCPion->SetLineColor(2);
        lineTPCPion->Draw();
        cNsigmaTPCPion_1D->SaveAs(output_QA_folder + ("/NsigmaTPCPion_1D." + outputtype).c_str());

        cNsigmaTPCPion_1D_before->cd();
        leg1DPID->Draw();
        lineTPCPion->SetY2(h1DNsigmaTPCPion_before[0]->GetMaximum());
        lineTPCPion->Draw();
        cNsigmaTPCPion_1D_before->SaveAs(output_QA_folder_before + ("/NsigmaTPCPion_1D_before." + outputtype).c_str());

        cNsigmaTPCPion_1D_pt->cd();
        leg1DPID_pt->Draw();
        lineTPCPion->SetY2(h1DNsigmaTPCPion_pt[0]->GetMaximum());
        lineTPCPion->Draw();
        cNsigmaTPCPion_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTPCPion_1D_pt." + outputtype).c_str());

        TCanvas *cNsigmaTOFKaon_1D = new TCanvas("cNsigmaTOFKaon_1D", "Nsigma TOF Kaon 1D", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon_1D, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFKaon_1D->SetGrid(1, 1);

        TCanvas *cNsigmaTOFKaon_1D_before = new TCanvas("cNsigmaTOFKaon_1D_before", "Nsigma TOF Kaon 1D", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon_1D_before, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFKaon_1D_before->SetGrid(1, 1);

        TCanvas *cNsigmaTOFKaon_1D_pt = new TCanvas("cNsigmaTOFKaon_1D_pt", "Nsigma TOF Kaon 1D", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon_1D_pt, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFKaon_1D_pt->SetGrid(1, 1);

        TF1 *fitGaus = new TF1("fitGaus", "gaus", -2, 2);
        fitGaus->SetParameter(0, 1);
        fitGaus->SetParameter(1, -0.8);
        fitGaus->SetParameter(2, 1);
        for (int i = 0; i < 6; i++)
        {
            cNsigmaTOFKaon_1D->cd();
            SetHistoQA(h1DNsigmaTOFKaon[i]);
            h1DNsigmaTOFKaon[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFKaon[i]->GetXaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");
            h1DNsigmaTOFKaon[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFKaon[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            h1DNsigmaTOFKaon[i]->SetMaximum(h1DNsigmaTOFKaon[i]->GetMaximum() * 15);
            h1DNsigmaTOFKaon[i]->Draw("p same");

            cNsigmaTOFKaon_1D_before->cd();
            SetHistoQA(h1DNsigmaTOFKaon_before[i]);
            h1DNsigmaTOFKaon_before[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFKaon_before[i]->GetXaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");
            h1DNsigmaTOFKaon_before[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFKaon_before[i]->GetXaxis()->SetRangeUser(-5, 5);
            h1DNsigmaTOFKaon_before[i]->SetMaximum(h1DNsigmaTOFKaon_before[i]->GetMaximum() * 15);
            h1DNsigmaTOFKaon_before[i]->Draw("p same");

            cNsigmaTOFKaon_1D_pt->cd();
            SetHistoQA(h1DNsigmaTOFKaon_pt[i]);
            h1DNsigmaTOFKaon_pt[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFKaon_pt[i]->GetXaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");
            h1DNsigmaTOFKaon_pt[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFKaon_pt[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            h1DNsigmaTOFKaon_pt[i]->SetMaximum(h1DNsigmaTOFKaon_pt[i]->GetMaximum() * 100);
            h1DNsigmaTOFKaon_pt[i]->Draw("p same");
        }
        cNsigmaTOFKaon_1D->cd();
        // h1DNsigmaTOFKaon[3]->Fit(fitGaus, "REBMS");
        h1DNsigmaTOFKaon[3]->SetStats(1);
        leg1DPID->Draw();
        TLine *lineTOFKaon = new TLine(0, 0, 0, h1DNsigmaTOFKaon[0]->GetMaximum());
        lineTOFKaon->SetLineStyle(2);
        lineTOFKaon->SetLineColor(2);
        lineTOFKaon->Draw();
        cNsigmaTOFKaon_1D->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_1D." + outputtype).c_str());

        cNsigmaTOFKaon_1D_before->cd();
        leg1DPID->Draw();
        lineTOFKaon->SetY2(h1DNsigmaTOFKaon_before[0]->GetMaximum());
        lineTOFKaon->Draw();
        cNsigmaTOFKaon_1D_before->SaveAs(output_QA_folder_before + ("/NsigmaTOFKaon_1D_before." + outputtype).c_str());

        cNsigmaTOFKaon_1D_pt->cd();
        leg1DPID_pt->Draw();
        lineTOFKaon->SetY2(h1DNsigmaTOFKaon_pt[0]->GetMaximum());
        lineTOFKaon->Draw();
        cNsigmaTOFKaon_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_1D_pt." + outputtype).c_str());

        TCanvas *cNsigmaTOFPion_1D = new TCanvas("cNsigmaTOFPion_1D", "Nsigma TOF Pion 1D", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion_1D, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFPion_1D->SetGrid(1, 1);

        TCanvas *cNsigmaTOFPion_1D_before = new TCanvas("cNsigmaTOFPion_1D_before", "Nsigma TOF Pion 1D", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion_1D_before, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFPion_1D_before->SetGrid(1, 1);

        TCanvas *cNsigmaTOFPion_1D_pt = new TCanvas("cNsigmaTOFPion_1D_pt", "Nsigma TOF Pion 1D", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion_1D_pt, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFPion_1D_pt->SetGrid(1, 1);

        for (int i = 0; i < 6; i++)
        {
            cNsigmaTOFPion_1D->cd();
            SetHistoQA(h1DNsigmaTOFPion[i]);
            h1DNsigmaTOFPion[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFPion[i]->GetXaxis()->SetTitle("n#sigma_{TOF} #pi^{#pm}");
            h1DNsigmaTOFPion[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFPion[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            h1DNsigmaTOFPion[i]->SetMaximum(h1DNsigmaTOFPion[i]->GetMaximum() * 15);
            h1DNsigmaTOFPion[i]->Draw("p same");

            cNsigmaTOFPion_1D_before->cd();
            SetHistoQA(h1DNsigmaTOFPion_before[i]);
            h1DNsigmaTOFPion_before[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFPion_before[i]->GetXaxis()->SetTitle("n#sigma_{TOF} #pi^{#pm}");
            h1DNsigmaTOFPion_before[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFPion_before[i]->GetXaxis()->SetRangeUser(-5, 5);
            h1DNsigmaTOFPion_before[i]->SetMaximum(h1DNsigmaTOFPion_before[i]->GetMaximum() * 15);
            h1DNsigmaTOFPion_before[i]->Draw("p same");

            cNsigmaTOFPion_1D_pt->cd();
            SetHistoQA(h1DNsigmaTOFPion_pt[i]);
            h1DNsigmaTOFPion_pt[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFPion_pt[i]->GetXaxis()->SetTitle("n#sigma_{TOF} #pi^{#pm}");
            h1DNsigmaTOFPion_pt[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFPion_pt[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            h1DNsigmaTOFPion_pt[i]->SetMaximum(h1DNsigmaTOFPion_pt[i]->GetMaximum() * 25);
            h1DNsigmaTOFPion_pt[i]->Draw("p same");
        }
        cNsigmaTOFPion_1D->cd();
        // h1DNsigmaTOFPion[5]->Fit(fitGaus, "REBMS");
        leg1DPID->Draw();
        TLine *lineTOFPion = new TLine(0, 0, 0, h1DNsigmaTOFPion[0]->GetMaximum());
        lineTOFPion->SetLineStyle(2);
        lineTOFPion->SetLineColor(2);
        lineTOFPion->Draw();
        cNsigmaTOFPion_1D->SaveAs(output_QA_folder + ("/NsigmaTOFPion_1D." + outputtype).c_str());

        cNsigmaTOFPion_1D_before->cd();
        leg1DPID->Draw();
        lineTOFPion->SetY2(h1DNsigmaTOFPion_before[0]->GetMaximum());
        lineTOFPion->Draw();
        cNsigmaTOFPion_1D_before->SaveAs(output_QA_folder_before + ("/NsigmaTOFPion_1D_before." + outputtype).c_str());

        cNsigmaTOFPion_1D_pt->cd();
        leg1DPID_pt->Draw();
        lineTOFPion->SetY2(h1DNsigmaTOFPion_pt[0]->GetMaximum());
        lineTOFPion->Draw();
        cNsigmaTOFPion_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTOFPion_1D_pt." + outputtype).c_str());
    }

    cout << "\n================= End of the code =================" << endl;
    cout << "Data file used: " << data_path << endl;
    cout << "MC file used: " << MCpath << endl;
    cout << "Selection used: " << (isINEL ? "INEL" : "INEL > 0") << endl;
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
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
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

// correct placement of TPC crossed rows
// string data_path = "459845/kstarqa/hInvMass"; // 2022 data
// string data_path = "459908/kstarqa_PIDKa2/hInvMass"; // 2023 data
// string data_path = "460233/kstarqa_PIDKa2/hInvMass"; // 2024 data

//**********************************MID Cuts******************************************
// string data_path = "477779/kstarqa_MID/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
// string data_path = "478015/kstarqa_MID/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
// string data_path = "477833/kstarqa/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

//*************************PID Variations for Kaon (without MID)**************************
// string data_path = "480317/kstarqa/hInvMass"; // LHC22_pass7_medium dataset, INEL > 0
// string data_path = "480447/kstarqa/hInvMass"; // LHC23_pass4_thin_small dataset, INEL > 0
// string data_path = "480657/kstarqa/hInvMass"; // LHC24_pass1_minBias dataset, INEL > 0

//================================QA checks==================
// string data_path = "585940/kstarqa_VertexTOFMatched/hInvMass";

//==============================Pt-dependent PID=======================
// string data_path = "586976/kstarqa/hInvMass"; // 2023
// string data_path = "586385/kstarqa/hInvMass"; // 2024

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

// ***********************TPC tune on data with NN ***********************************
// string MCpath = "TuneOnDataWithNN/465413.root"; // 2022 MC
// string MCpath = "TuneOnDataWithNN/464570.root"; // 2023 MC
// string MCpath = "TuneOnDataWithNN/465419.root"; // 2023 MC

//********************MC closure test*************************/
// string MCpath = "Occupancy_effect/466809.root"; // 2022 MC
// string MCpath = "Occupancy_effect/463585.root"; // 2023 MC
// string MCpath = "Occupancy_effect/464277.root"; // 2024 MC

//**************Checks (Default, FakeTrack, MID, PIDKa2) and new MC closure test********************
// string MCpath = "checks/472785.root"; // 2022 MC
// string MCpath = "checks/472820.root"; // 2023 MC
// string MCpath = "checks/472788.root"; // 2024 MC

//********************************MID Cuts******************************************
// string MCpath = "478021.root"; // 2022 MC
// string MCpath = "477593.root"; // 2023 MC
// string MCpath = "478022.root"; // 2024 MC

////*************************PID Variations for Kaon (with MID)**************************
// string MCpath = "480451.root"; // 2022 MC
// string MCpath = "480565.root"; // 2023 MC
// string MCpath = "480453.root"; // 2024 MC

//**************************After Calibrated MC from Nicolo****************************
// string MCpath = "483982.root"; // 2023 MC
// string MCpath = "temp.root"; // 2023 MC
// string MCpath = "484450.root"; // 2024 MC

//*******************************************QA checks************************************************
// string MCpath = "585904.root"; // 2023 MC

//================================Pt-dependent PID=======================
// string MCpath = "586966.root"; // 2023 MC
// string MCpath = "586467.root"; // 2024 MC
