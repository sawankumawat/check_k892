
#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/initializations.h"
#include <TSystem.h>
#include <TString.h>
#include <TStopwatch.h>

using namespace std;

void kstar_sparse()
{
    TStopwatch timer;
    timer.Start();
    //*************************** change here ***************************************
    const string kResBkg = "MIX";
    // const string kResBkg = "LIKE";
    // const string kResBkg = "ROTATED";
    const string kbkg = "pol3";
    const string outputtype = "pdf"; // pdf, eps
    const bool save_bkg_plots = 1;   // save background plots
    const float txtsize = 0.045;     // text size in the plots
    bool makeQAplots = false;
    bool makeallpTplots = true; // make all pT plots
    bool calcInvMass = true;

    //********************************************************************************

    //*************************Create folders********************************************
    TString outputfolder = kSignalOutput + "/" + kfoldername;
    TString output_QA_folder = kSignalOutput + "/" + kfoldername + "/QA";
    // Create the folder using TSystem::mkdir()
    if (gSystem->mkdir(outputfolder, kTRUE))
    {
        std::cout << "Folder " << outputfolder << " created successfully." << std::endl;
    }
    if (gSystem->mkdir(output_QA_folder, kTRUE))
    {
        std::cout << "Folder " << output_QA_folder << " created successfully." << std::endl;
    }
    //***********************************************************************************

    TCanvas *cgrid1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid_bkg1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid_bkg2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);

    // some initializations ********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************

    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.06);
    t2->SetTextFont(42);

    TCanvas *cinv[Npt]; // for output canvases on screen containing fitted signal after subtraction

    for (Int_t ip = 0; ip < Npt; ip++)
    {
        TString cName = TString::Format("cinv_pt_%2.1f-%2.1f", pT_bins[ip], pT_bins[ip + 1]);
        cinv[ip] = new TCanvas(Form("cinv%d", ip), cName.Data(), 10, 10, 720, 720);
        SetCanvasStyle(cinv[ip], 0.15, 0.05, 0.08, 0.13);
    }

    TCanvas *cSigbkg[Npt]; // for output canvases on screen containing signal with bkg(after norm in case of mix)

    for (Int_t ip = 0; ip < Npt; ip++)
    {
        TString cNam = TString::Format("cSigbkg_pt_%2.1f-%2.1f", pT_bins[ip], pT_bins[ip + 1]);
        cSigbkg[ip] = new TCanvas(Form("cSigbkg%d", ip), cNam.Data(), 720, 720);
        SetCanvasStyle(cSigbkg[ip], 0.15, 0.06, 0.06, 0.13);
    }

    if (multipanel_plots)
    {
        cgrid1->Divide(kcanvasdivide[0], kcanvasdivide[1]);
        cgrid2->Divide(kcanvasdivide[0], kcanvasdivide[1]);
        cgrid_bkg1->Divide(kcanvasdivide[0], kcanvasdivide[1]);
        cgrid_bkg2->Divide(kcanvasdivide[0], kcanvasdivide[1]);
    }
    Double_t significance_den, significance_num, ratio, ratio2;

    //********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************

    // Input file
    TFile *fInputFile = new TFile(kDataFilename.c_str(), "Read");
    if (fInputFile->IsZombie())
    {
        cerr << "File not found " << endl;
        return;
    }

    // TH1F *hmult = (TH1F *)fInputFile->Get("kstarqa_id21631/eventSelection/hMultiplicity");
    // TH1F *hmult = (TH1F *)fInputFile->Get("kstarqa/eventSelection/hMultiplicity");
    string multpath = kfoldername.substr(0, kfoldername.length() - 9);
    TH1F *hmult = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hMultiplicity", multpath.c_str()));
    // cout << "Given path is " << kfoldername << endl;
    if (hmult == nullptr)
    {
        cerr << "Histogram not found" << endl;
        return;
    }
    double Event = hmult->GetEntries();
    cout << "*****************number of events********************:" << Event << endl;

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    // float mult_classes[] = {0.0, 100.0};
    int nmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1; // number of multiplicity bins
    int rebin_value;

    THnSparseF *fHistNum, *fHistDen, *fHistLSPP, *fHistLSMM, *fHistRotated;

    cout << "path for invariant mass histograms is " << Form("%s/h3KstarInvMassUnlikeSign", kfoldername.c_str()) << endl;

    fHistNum = (THnSparseF *)fInputFile->Get(Form("%s/h3KstarInvMassUnlikeSign", kfoldername.c_str()));
    if (kResBkg == "MIX")
    {
        fHistDen = (THnSparseF *)fInputFile->Get(Form("%s/h3KstarInvMassMixed", kfoldername.c_str()));
    }
    if (kResBkg == "LIKE")
    {
        fHistLSPP = (THnSparseF *)fInputFile->Get(Form("%s/h3KstarInvMasslikeSignPP", kfoldername.c_str()));
        fHistLSMM = (THnSparseF *)fInputFile->Get(Form("%s/h3KstarInvMasslikeSignMM", kfoldername.c_str()));
    }
    if (kResBkg == "ROTATED")
    {
        fHistRotated = (THnSparseF *)fInputFile->Get(Form("%s/h3KstarInvMassRotated", kfoldername.c_str()));
    }

    if (fHistNum == nullptr)
    {
        cerr << "Invariant mass histograms not found!!!!!!!!!!!!" << endl;
        return;
    }
    TFile *filecmp = new TFile((koutputfolder + "/yield.root").c_str(), "RECREATE");

    for (int imult = 0; imult < nmultbins + 1; imult++)
    // for (int imult = 0; imult < 1; imult++)
    {
        // basic checks
        if (kNormRangepT.size() < Npt || kFitRange.size() < Npt || kRebin.size() < Npt)
        {
            cerr << "Error: kNormRangepT, kFitRange, or kRebin arrays are not initialized for all pT bins." << endl;
            return;
        }
        //**************Invariant mass histograms for sig+bkg and mixed event bg******************
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

        double Event = hmult->Integral(hmult->GetXaxis()->FindBin(multlow + 1e-3), hmult->GetXaxis()->FindBin(multhigh - 1e-3));
        cout << "Event in mult bin " << imult << " is " << Event << endl;

        TString outputfolder_mult = kSignalOutput + "/" + kfoldername + Form("/mult_%d-%d", multlow, multhigh);
        if (gSystem->mkdir(outputfolder_mult, kTRUE))
        {
            std::cout << "Folder " << outputfolder_mult << " created successfully." << std::endl;
        }

        if (calcInvMass)
        {
            // gStyle->SetOptStat(1110);
            gStyle->SetOptStat(0);
            gStyle->SetOptFit(1);

            // Check if directory exists, if not create it
            TDirectory *dir = filecmp->GetDirectory(Form("mult_%d-%d", (int)multlow, (int)multhigh));
            if (!dir)
            {
                dir = filecmp->mkdir(Form("mult_%d-%d", (int)multlow, (int)multhigh));
            }
            filecmp->cd();
            dir->cd();
            TDirectory *subdir = dir->GetDirectory("SignalInAllPtBins");
            if (!subdir)
            {
                subdir = dir->mkdir("SignalInAllPtBins");
            }
            subdir->cd();

            for (Int_t ip = pt_start; ip < pt_end; ip++) // start pt bin loop
            {
                rebin_value = kRebin[ip][imult]; // rebinning value for the multiplicity bin
                // rebin_value = 2;                 // for medium dataset (temporarily set to 1)

                double lowfitrange = kFitRange[ip][0];
                double highfitrange = kFitRange[ip][1];

                if (imult == nmultbins)
                {
                    lowfitrange = kFitRange70to100[ip][0]; // for last mult bin, extend the range a bit
                    highfitrange = kFitRange70to100[ip][1];
                }
                if (imult == 1 && ip == 0)
                { // for multiplicity 0-1% and first pt bin
                    lowfitrange = 0.77;
                    highfitrange = 1.03;
                }

                lowpt = pT_bins[ip];
                highpt = pT_bins[ip + 1];
                int lbin = fHistNum->GetAxis(1)->FindBin(lowpt + 1e-3);
                int hbin = fHistNum->GetAxis(1)->FindBin(highpt - 1e-3);

                fHistNum->GetAxis(1)->SetRange(lbin, hbin);
                if (kResBkg == "MIX")
                {
                    fHistDen->GetAxis(1)->SetRange(lbin, hbin);
                }
                if (kResBkg == "LIKE")
                {
                    fHistLSPP->GetAxis(1)->SetRange(lbin, hbin);
                    fHistLSMM->GetAxis(1)->SetRange(lbin, hbin);
                }
                if (kResBkg == "ROTATED")
                {
                    fHistRotated->GetAxis(1)->SetRange(lbin, hbin);
                }

                int lbinmult = fHistNum->GetAxis(0)->FindBin(multlow + 1e-3);
                int hbinmult = fHistNum->GetAxis(0)->FindBin(multhigh - 1e-3);

                fHistNum->GetAxis(0)->SetRange(lbinmult, hbinmult);
                if (kResBkg == "MIX")
                {
                    fHistDen->GetAxis(0)->SetRange(lbinmult, hbinmult);
                }
                if (kResBkg == "LIKE")
                {
                    fHistLSPP->GetAxis(0)->SetRange(lbinmult, hbinmult);
                    fHistLSMM->GetAxis(0)->SetRange(lbinmult, hbinmult);
                }
                if (kResBkg == "ROTATED")
                {
                    fHistRotated->GetAxis(0)->SetRange(lbinmult, hbinmult);
                }

                fHistTotal[ip] = fHistNum->Projection(2, "E");
                if (kResBkg == "MIX")
                {
                    fHistBkg[ip] = fHistDen->Projection(2, "E");
                }
                if (kResBkg == "LIKE")
                {
                    fHistbkgLSPP[ip] = fHistLSPP->Projection(2, "E");
                    fHistbkgLSMM[ip] = fHistLSMM->Projection(2, "E");
                }
                if (kResBkg == "ROTATED")
                {
                    fHistRotated1D[ip] = fHistRotated->Projection(2, "E");
                }

                if (kResBkg == "LIKE")
                {
                    // Initialize fHistbkgLS[ip] by cloning one of the existing histograms
                    fHistbkgLS[ip] = (TH1D *)fHistbkgLSPP[ip]->Clone(Form("fHistbkgLS_%d_%d", imult, ip));
                    fHistbkgLS[ip]->Reset(); // Clear the content, keep the binning structure

                    for (int ibin = 0; ibin < fHistbkgLSPP[ip]->GetNbinsX(); ibin++)
                    {
                        double linkesignpp = fHistbkgLSPP[ip]->GetBinContent(ibin + 1);
                        double linkesignmm = fHistbkgLSMM[ip]->GetBinContent(ibin + 1);
                        fHistbkgLS[ip]->SetBinContent(ibin + 1, 2 * sqrt(linkesignpp * linkesignmm));
                        double binerrorpp = fHistbkgLSPP[ip]->GetBinError(ibin + 1);
                        double binerrormm = fHistbkgLSMM[ip]->GetBinError(ibin + 1);
                        fHistbkgLS[ip]->SetBinError(ibin + 1, sqrt(linkesignmm / linkesignpp) * binerrorpp + sqrt(linkesignpp / linkesignmm) * binerrormm);
                    }
                }
                fHistNum->SetName(Form("fHistNum_%d_%d", imult, ip));
                if (kResBkg == "MIX")
                {
                    fHistDen->SetName(Form("fHistDen_%d_%d", imult, ip));
                }
                if (kResBkg == "LIKE")
                {
                    fHistLSPP->SetName(Form("fHistLSPP_%d_%d", imult, ip));
                    fHistLSMM->SetName(Form("fHistLSMM_%d_%d", imult, ip));
                }
                if (kResBkg == "ROTATED")
                {
                    fHistRotated1D[ip]->SetName(Form("fHistRotated_%d_%d", imult, ip));
                }

                auto energylow = fHistTotal[ip]->GetXaxis()->GetXmin();
                auto energyhigh = fHistTotal[ip]->GetXaxis()->GetXmax();

                // cout<<"energy low value is "<<energylow<<endl;
                // cout<<"energy high value is "<<energyhigh<<endl;

                //**Cloning sig+bkg histogram for like sign or mixed event subtraction *********************************************************
                TH1D *hfsig = (TH1D *)fHistTotal[ip]->Clone();
                auto binwidth_file = (fHistTotal[ip]->GetXaxis()->GetXmax() - fHistTotal[ip]->GetXaxis()->GetXmin()) * rebin_value / fHistTotal[ip]->GetXaxis()->GetNbins();
                cout << "The value of binwidth_file is: " << binwidth_file << endl;
                //*****************************************************************************************************************************

                if (kResBkg == "MIX" || kResBkg == "ROTATED")
                {
                    TH1D *bkgclonetemp = (kResBkg == "MIX") ? (TH1D *)fHistBkg[ip]->Clone() : (TH1D *)fHistRotated1D[ip]->Clone();

                    sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                    bkg_integral = (bkgclonetemp->Integral(bkgclonetemp->GetXaxis()->FindBin(kNormRangepT[ip][0]), bkgclonetemp->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                    normfactor = sigbkg_integral / bkg_integral; // scaling factor for mixed bkg
                    cout << "\n\n normalization factor " << 1 / normfactor << "\n\n";
                    hfbkg = (TH1D *)bkgclonetemp->Clone();
                    hfbkg->Scale(normfactor);
                    hfbkg->Rebin(rebin_value);
                    hfsig->Rebin(rebin_value);
                    hfsig->Add(hfbkg, -1);
                }
                else if (kResBkg == "LIKE")
                {
                    // sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                    // bkg_integral = (fHistbkgLS[ip]->Integral(fHistbkgLS[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistbkgLS[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
                    // normfactor = sigbkg_integral / bkg_integral; // scaling factor for mixed bkg
                    // hfbkg = (TH1D *)fHistbkgLS[ip]->Clone();
                    // hfbkg->Scale(normfactor);

                    hfbkg = (TH1D *)fHistbkgLS[ip]->Clone();
                    hfbkg->Rebin(rebin_value);
                    hfsig->Rebin(rebin_value);
                    hfsig->Add(hfbkg, -1);
                }

                fHistTotal[ip]->Rebin(rebin_value);

                //**** pt binwidth************x*****************************
                ptbinwidth[ip] = pT_bins[ip + 1] - pT_bins[ip];
                // cout<<"the value of pt bin width is "<<ptbinwidth[ip]<<endl;

                //****************************************************************************************************

                TF1 *fitFcn, *fitFcn1;

                if (kbkg == "pol2")
                {
                    fitFcn = new TF1("fitfunc", BreitWignerpoly2, lowfitrange, highfitrange, 6);
                    fitFcn1 = new TF1("fitfunc1", polynomial2, lowfitrange, highfitrange, 3);
                }
                else if (kbkg == "pol3")
                {
                    fitFcn = new TF1("fitfunc", BreitWignerpoly3, lowfitrange, highfitrange, 7);
                    fitFcn1 = new TF1("fitfunc1", polynomial3, lowfitrange, highfitrange, 4);
                }
                else if (kbkg == "expol")
                {
                    fitFcn = new TF1("fitfunc", BWExpo, lowfitrange, highfitrange, 7);
                    fitFcn1 = new TF1("fitfunc1", Expo, lowfitrange, highfitrange, 4);
                }

                TF1 *fitFcn2 = new TF1("fitFcn2", BW, lowfitrange, highfitrange, 3); // only signal

                fitFcn->SetParameter(0, masspdg); // mass
                if (imult == 1 && ip < 3)
                {
                    fitFcn->SetParLimits(0, 0.885, 0.888); // Mass
                }
                else if (imult == 1 && ip > Npt - 4)
                {
                    fitFcn->SetParLimits(0, 0.895, 0.898); // Mass
                }
                else if (imult == nmultbins && ip > Npt - 4)
                {
                    fitFcn->SetParLimits(0, 0.895, 0.898); // Mass
                }
                else
                {
                    fitFcn->SetParLimits(0, 0.80, 0.98); // Mass
                }
                fitFcn->FixParameter(1, widthpdg); // width
                // fitFcn->SetParLimits(1, 0.03, 0.07); // width
                fitFcn->SetParameter(2, 1000);     // yield
                fitFcn->SetParLimits(2, 0.0, 1e8); // Yield

                // // //pol3 parameters
                // fitFcn->SetParameter(3, -1e6);
                // fitFcn->SetParameter(4, 1e6);
                // fitFcn->SetParameter(5, -1e6);
                // fitFcn->SetParameter(6, 1e6);

                fitFcn->SetParNames("Mass", "Width", "Yield", "A", "B", "C", "D");
                // Redirect standard output to /dev/null
                // int old_stdout = dup(1);
                // freopen("/dev/null", "w", stdout);

                r = hfsig->Fit(fitFcn, "REBMS0+"); // signal after bkg subtraction

                // Restore standard output
                // fflush(stdout);
                // dup2(old_stdout, 1);
                // close(old_stdout);

                //****************************************************************************************************************

                //**Extraction of fitting parameters******************************************************************************

                Double_t *par = fitFcn->GetParameters();

                Mass[ip] = fitFcn->GetParameter(0);
                Width[ip] = fitFcn->GetParameter(1);
                Yield[ip] = fitFcn->GetParameter(2);
                poly2[ip] = fitFcn->GetParameter(3);
                poly1[ip] = fitFcn->GetParameter(4);
                poly0[ip] = fitFcn->GetParameter(5);
                if (kbkg == "pol3" || kbkg == "expol")
                    poly3[ip] = fitFcn->GetParameter(6);

                fitFcn2->SetParameters(&par[0]);
                fitFcn1->SetParameters(&par[3]);

                ErrorMass[ip] = fitFcn->GetParError(0);
                ErrorWidth[ip] = fitFcn->GetParError(1);
                ErrorYield[ip] = fitFcn->GetParError(2);
                Chi2Ndf[ip] = (fitFcn->GetChisquare()) / (fitFcn->GetNDF());

                //******************************************************************************************************************

                //**ERROR BIN COUNTING METHOD CALCULATION*****************************************************************************

                TF1 *fitFcn2_plusm = new TF1("fitFcn2_plusm", BW, lowfitrange, highfitrange, 3);
                TF1 *fitFcn2_minusm = new TF1("fitFcn2_minusm", BW, lowfitrange, highfitrange, 3);
                fitFcn2_plusm->FixParameter(0, Mass[ip] + ErrorMass[ip]);
                fitFcn2_plusm->FixParameter(1, 0.047);
                fitFcn2_plusm->FixParameter(2, Yield[ip]);

                fitFcn2_minusm->FixParameter(0, Mass[ip] - ErrorMass[ip]);
                fitFcn2_minusm->FixParameter(1, 0.047);
                fitFcn2_minusm->FixParameter(2, Yield[ip]);

                //*********************************************************************************************************************

                //**Calculation of significance and storing chi2 and sig in respective histograms*****************************************************

                bmin = hfsig->GetXaxis()->FindBin(masspdg - 2 * widthpdg);
                bmax = hfsig->GetXaxis()->FindBin(masspdg + 2 * widthpdg);

                significance_num = (fitFcn2->Integral(masspdg - 2 * widthpdg, masspdg + 2 * widthpdg)) / (binwidth_file);
                significance_den = TMath::Sqrt(fHistTotal[ip]->Integral(bmin, bmax));

                ratio = significance_num / significance_den; // significance of signal

                hsignificance->SetBinContent(ip + 1, ratio);
                hChiSquare->SetBinContent(ip + 1, Chi2Ndf[ip]); // storing both significance and chi2 in histogram

                //*****************************************************************************************************************************************

                //**Calculation of Yield using bin counting method and storing it in histogram***********************************************************

                Yield_bincount_hist = hfsig->IntegralAndError(bmin, bmax, hBCError_1);
                bkgvalue = fitFcn1->Integral(hfsig->GetBinLowEdge(bmin), hfsig->GetBinLowEdge(bmax + 1));
                Integral_BW_withsigma = fitFcn2->Integral(hfsig->GetBinLowEdge(bmin), hfsig->GetBinLowEdge(bmax + 1));
                fYield_BinCount = Yield_bincount_hist - (bkgvalue / binwidth_file);
                YieldIntegral_BW = fitFcn2->Integral(energylow, energyhigh) / binwidth_file;
                Yfraction_cBW = (Integral_BW_withsigma / YieldIntegral_BW);

                sum_tail_correction = (fitFcn2->Integral(energylow, hfsig->GetBinLowEdge(bmin)) + fitFcn2->Integral(hfsig->GetBinLowEdge(bmax + 1), energyhigh)) / binwidth_file;

                nlow = (fitFcn2->Integral(energylow, hfsig->GetBinLowEdge(bmin))) / binwidth_file;
                nhigh = (fitFcn2->Integral(hfsig->GetBinLowEdge(bmax + 1), energyhigh)) / binwidth_file;
                nlow = nlow / (Event * ptbinwidth[ip] * dy * BR);
                nhigh = nhigh / (Event * ptbinwidth[ip] * dy * BR);

                Total_Ybincounting = (sum_tail_correction + fYield_BinCount) / (Event * ptbinwidth[ip] * dy * BR);

                // cout << "***************************************************************" << endl;
                // cout << "****fraction of nlow for bin***********:"
                //      << " " << ip << " " << nlow / Total_Ybincounting << endl;
                // cout << "****fraction of nhigh for bin***********:"
                //      << " " << ip << " " << nhigh / Total_Ybincounting << endl;
                // cout << "***************************************************************" << endl;
                Tail_correction_plusm = (fitFcn2_plusm->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + (fitFcn2_plusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5))) / binwidth_file;
                Tail_correction_minusm = ((fitFcn2_minusm->Integral(0.635, hfsig->GetBinLowEdge(bmin)) + fitFcn2_minusm->Integral(hfsig->GetBinLowEdge(bmax + 1), 5)) / binwidth_file);
                Error_2 = sum_tail_correction - Tail_correction_plusm;
                Final_pro_error = TMath::Sqrt(Error_2 * Error_2 + hBCError_1 * hBCError_1) / (Event * ptbinwidth[ip] * dy * BR);

                ////Uncorrected Yield/////////////////////////////////////////////////////////////////////////////////

                hYbincount->SetBinContent(ip + 1, Total_Ybincounting);
                hYbincount->SetBinError(ip + 1, Final_pro_error);
                // cout << "--------Total Value from bin counting----------" << (sum_tail_correction + fYield_BinCount) << endl;
                // cout << "--------Value from bin counting----------" << Total_Ybincounting << endl;

                //////////////////////////////////////////////////////////////////////////////////////////////////////

                // Fractional stat error///////////////////////////////////////////////////////////////////////////////

                hFrac_stat_error->SetBinContent(ip + 1, Final_pro_error / Total_Ybincounting);
                // cout << "--------Frac error from bin counting----------" << (Final_pro_error / Total_Ybincounting) << endl;
                //////////////////////////////////////////////////////////////////////////////////////////////////////

                //****************************************************************************************************************************************

                //**Calculation for raw pt spectra using function integration and filling it in histogram*********************************************

                integralsignalfunc[ip] = (fitFcn2->Integral((masspdg - 5 * widthpdg), (masspdg + 5 * widthpdg)));
                TMatrixDSym cov = r->GetCovarianceMatrix();
                TMatrixDSym cov1;
                TMatrixDSym cov2;
                cov.GetSub(0, 2, 0, 2, cov1);
                cov.GetSub(3, 6, 3, 6, cov2);
                Double_t *b = cov1.GetMatrixArray();
                Double_t *a = cov2.GetMatrixArray();
                Double_t *para = fitFcn->GetParameters();
                interror[ip] = fitFcn2->IntegralError((masspdg - 5 * widthpdg), (masspdg + 5 * widthpdg), &para[0], b);

                yieldcalc = integralsignalfunc[ip] / (Event * ptbinwidth[ip] * dy * BR * binwidth_file); // raw yield calculation
                yielderror = interror[ip] / (Event * ptbinwidth[ip] * dy * BR * binwidth_file);          // raw yield error

                hintegral_yield->SetBinContent(ip + 1, yieldcalc);
                hintegral_yield->SetBinError(ip + 1, yielderror); // filling histogram including error

                //**Filling mass and width fitting parameter in histogram*******************************************************************************

                hmass->SetBinContent(ip + 1, Mass[ip]);
                hmass->SetBinError(ip + 1, ErrorMass[ip]);

                hwidth->SetBinContent(ip + 1, Width[ip]);
                hwidth->SetBinError(ip + 1, ErrorWidth[ip]);

                //*****************************************************************************************************************************

                //**Setting plot parameters style*************************************************************************************************

                // SetHistoStyle(hfsig, 1, 20, 1, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
                // SetHistoStyle(fHistTotal[ip], 1, 8, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);
                SetHistoQA(hfsig);
                SetHistoQA(fHistTotal[ip]);

                hfsig->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
                hfsig->GetYaxis()->SetMaxDigits(2);
                hfsig->GetYaxis()->CenterTitle(1);
                hfsig->GetYaxis()->SetTitle(Form("Counts/%.0f MeV/#it{c}^{2}", binwidth_file * 1000));

                SetHistoQA(hfbkg);
                hfbkg->SetLineColor(kRed);
                hfbkg->SetMarkerColor(kRed);
                hfbkg->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
                // hfsig->GetYaxis()->SetMaxDigits(2);
                hfbkg->GetYaxis()->SetTitle(Form("Counts/%.0f MeV/#it{c}^{2}", binwidth_file * 1000));

                fitFcn1->SetLineColor(4);
                fitFcn1->SetLineStyle(2);
                fitFcn1->SetLineWidth(4);
                fitFcn2->SetLineColor(6);
                fitFcn2->SetLineStyle(2);
                fitFcn2->SetLineWidth(4);
                fitFcn->SetLineWidth(4);

                //*******************************************************************************************************************************

                //**Plot of histograms and graphs*********************************************************************************************
                auto chibyndf = fitFcn->GetChisquare() / fitFcn->GetNDF();

                // inv mass histograms after the background subraction
                (multipanel_plots == 1) ? (ip < kupperpad * klowerpad) ? cgrid1->cd(ip + 1) : cgrid2->cd(ip + 1 - kupperpad * klowerpad) : cinv[ip]->cd();
                gPad->SetRightMargin(0.015);
                gPad->SetLeftMargin(0.15);
                gPad->SetBottomMargin(0.15);
                // if (ip == 0)
                // {
                //     hfsig->SetMaximum(hfsig->GetMaximum() * 1.8);
                // }
                // else if (ip == 1)
                //     hfsig->SetMaximum(hfsig->GetMaximum() * 1);
                // else if (ip == 2)
                //     hfsig->SetMaximum(hfsig->GetMaximum() * 0.9);

                // else
                hfsig->SetMaximum(hfsig->GetMaximum() * 1.1);
                fitFcn->SetLineWidth(2);
                fitFcn1->SetLineWidth(2);
                fitFcn2->SetLineWidth(2);
                hfsig->GetXaxis()->SetRangeUser(0.7, 1.11);
                hfsig->SetMarkerSize(1.0);
                hfsig->Draw("e");
                fitFcn->Draw("same");
                fitFcn1->Draw("same");
                fitFcn2->Draw("same");

                TLegend *pag = new TLegend(0.14, 0.6, 0.35, 0.9);
                pag->SetBorderSize(0);
                pag->SetTextFont(42);
                pag->SetTextSize(0.035);
                pag->SetFillStyle(0);
                pag->AddEntry((TObject *)0, "ALICE", "");
                pag->AddEntry((TObject *)0, "pp, #sqrt{s} = 13.6 TeV", "");
                pag->AddEntry((TObject *)0, Form("FT0M, %d-%d%%", multlow, multhigh), "");
                pag->AddEntry((TObject *)0, "K*^{0}#rightarrow K#pi", "");
                pag->AddEntry((TObject *)0, Form("#it{p}_{T} = %.1f-%.1f GeV/c", pT_bins[ip], pT_bins[ip + 1]), "");
                pag->Draw();

                // TLegend *pag2 = new TLegend(0.2, 0.7, 0.45, 0.9);
                TLegend *pag2 = new TLegend(0.7, 0.65, 0.9, 0.9);
                pag2->SetBorderSize(0);
                pag2->SetTextFont(42);
                pag2->SetTextSize(0.035);
                pag2->SetFillStyle(0);
                pag2->AddEntry(fitFcn, "BW+pol3", "l");
                pag2->AddEntry(fitFcn1, "BW", "l");
                pag2->AddEntry(fitFcn2, "pol3", "l");
                // pag2->Draw();

                double fitprob = fitFcn->GetProb();
                // pag->AddEntry((TObject *)0, Form("Mass: %.3f #pm %.3f", Mass[ip], ErrorMass[ip]), "");
                // pag->AddEntry((TObject *)0, Form("Width: %.3f #pm %.3f", Width[ip], ErrorWidth[ip]), "");
                // pag->AddEntry((TObject *)0, Form("Yield: %.1e #pm %.1e", yieldcalc * Event, yielderror * Event), "");
                // pag->AddEntry((TObject *)0, Form("Probability: %f ", fitprob), "");
                // pag->AddEntry((TObject *)0, Form("#chi^{2}/NDF: %.2f ", chibyndf), "");
                // pag->Draw();
                // t2->DrawLatex(0.26, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
                // t2->DrawLatex(0.27, 0.95, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip], pT_bins[ip + 1]));
                if (multipanel_plots == 0 && save_plots == 1)
                    cinv[ip]->SaveAs(Form(outputfolder_mult + ("/hfitsig_pt%d." + outputtype).c_str(), ip + 1));
                if (multipanel_plots == 1)
                    cinv[ip]->Close();

                hfsig->Write(Form("hfsig_pt%d", ip + 1));

                // inv distribution before the background subtraction
                (multipanel_plots == 1) ? (ip < klowerpad * kupperpad) ? cgrid_bkg1->cd(ip + 1) : cgrid_bkg2->cd(ip + 1 - klowerpad * kupperpad) : cSigbkg[ip]->cd();
                TH1F *hbkg_nopeak = (TH1F *)hfbkg->Clone();
                hbkg_nopeak->SetLineColor(kRed);
                hbkg_nopeak->SetMarkerColor(kRed);
                hbkg_nopeak->SetFillColor(kRed);
                hbkg_nopeak->SetFillStyle(3001);
                for (int i = 0; i < hbkg_nopeak->GetNbinsX(); i++)
                {
                    if (hbkg_nopeak->GetBinCenter(i + 1) < kNormRangepT[ip][0] || hbkg_nopeak->GetBinCenter(i + 1) > kNormRangepT[ip][1])
                    {
                        hbkg_nopeak->SetBinContent(i + 1, -999);
                    }
                }
                gPad->SetRightMargin(0.04);
                gPad->SetLeftMargin(0.15);
                gPad->SetBottomMargin(0.15);
                fHistTotal[ip]->SetMaximum(fHistTotal[ip]->GetMaximum() * 1.15);
                fHistTotal[ip]->SetMarkerSize(1.0);
                hfbkg->SetMarkerSize(1.0);
                fHistTotal[ip]->GetXaxis()->SetRangeUser(0.7, 1.3);
                fHistTotal[ip]->Draw("pe");
                fHistTotal[ip]->GetYaxis()->SetTitle(Form("Counts/(%.0f MeV/#it{c}^{2})", binwidth_file * 1000));

                TLatex *ltx = new TLatex(0.27, 0.95, name);
                ltx->SetNDC();
                ltx->SetTextFont(22);
                ltx->SetTextSize(0.06);
                hfbkg->Draw("E same");
                fHistTotal[ip]->SetMarkerSize(1.0);
                hfbkg->SetMarkerSize(1.0);
                if (kResBkg == "MIX" || kResBkg == "ROTATED")
                    hbkg_nopeak->Draw("BAR same");
                // (kResBkg == "MIX") ? leg112->AddEntry(hfbkg, "Mixed-event bkg", "p") : leg112->AddEntry(hfbkg, "Like sign pairs", "p");
                ltx->Draw();
                pag->Draw();

                TLegend *pag3 = new TLegend(0.7, 0.65, 0.9, 0.9);
                pag3->SetBorderSize(0);
                pag3->SetTextFont(42);
                pag3->SetTextSize(0.035);
                pag3->SetFillStyle(0);
                pag3->AddEntry(fHistTotal[ip], "Signal + Bkg", "p");
                pag3->AddEntry(hfbkg, "Bkg", "p");
                pag3->AddEntry(hbkg_nopeak, "Norm. region", "f");
                pag3->Draw();

                if (multipanel_plots == 0 && save_plots == 1 && save_bkg_plots == 1)
                    cSigbkg[ip]->SaveAs(Form(outputfolder_mult + ("/hsigbkg_pt%d." + outputtype).c_str(), ip + 1));
                // cSigbkg[ip]->Close();

                // ////////////////////////////////////////////////////////////////////////

            } // pt loop ends
            if (multipanel_plots == 1 && save_plots == 1)
            {
                cgrid1->SaveAs(outputfolder_mult + ("/grid1." + outputtype).c_str());
                cgrid2->SaveAs(outputfolder_mult + ("/grid2." + outputtype).c_str());
                cgrid_bkg1->SaveAs(outputfolder_mult + ("/grid_bkg1." + outputtype).c_str());
                cgrid_bkg2->SaveAs(outputfolder_mult + ("/grid_bkg2." + outputtype).c_str());
            }
            if (multipanel_plots == 0)
            {
                cgrid1->Close();
                cgrid2->Close();
                cgrid_bkg1->Close();
                cgrid_bkg2->Close();
            }
            dir->cd();

            if (makeallpTplots)
            {
                TCanvas *csig = new TCanvas("", "", 720, 720);
                SetCanvasStyle(csig, 0.16, 0.05, 0.055, 0.15);
                SetHistoQA(hsignificance);
                hsignificance->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hsignificance->GetYaxis()->SetTitle("Significance");
                hsignificance->SetStats(0);
                hsignificance->Draw();
                hsignificance->Write("significance");
                csig->SaveAs(outputfolder_mult + ("/significance." + outputtype).c_str());
                csig->Clear();

                // // // chisquare_NDF vs pt
                hChiSquare->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hChiSquare->GetYaxis()->SetTitle("#chi^{2}/NDF ");
                SetHistoQA(hChiSquare);
                hChiSquare->SetMaximum(10.5);
                hChiSquare->SetMinimum(0);
                hChiSquare->SetStats(0);
                hChiSquare->Draw("p");
                // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
                csig->SaveAs(outputfolder_mult + ("/chi." + outputtype).c_str());
                hChiSquare->Write("chi2byNDF");
                csig->Clear();

                // // mass vs pt
                hmass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hmass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
                SetHistoQA(hmass);
                hmass->GetYaxis()->SetRangeUser(0.878, 0.909);
                hmass->SetStats(0);
                hmass->Draw("pe");
                hmass->Write("mass");
                TLegend *massleg = new TLegend(0.65, 0.2, 0.9, 0.3);
                SetLegendStyle(massleg);
                massleg->SetTextSize(txtsize);
                // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
                TLine *line = new TLine(hmass->GetXaxis()->GetXmin(), 0.895, hmass->GetXaxis()->GetXmax(), 0.895);
                line->SetLineStyle(2);
                line->SetLineColor(2);
                line->SetLineWidth(3);
                line->Draw();
                massleg->AddEntry(line, "PDG Mass", "l");
                massleg->Draw("l");
                csig->SaveAs(outputfolder_mult + ("/mass." + outputtype).c_str());
                csig->Clear();

                // // // Width vs pT
                hwidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hwidth->GetYaxis()->SetTitle("Width (GeV)");
                SetHistoQA(hwidth);
                hwidth->SetMaximum(hwidth->GetMaximum() * 2);
                hwidth->SetMinimum(0);
                hwidth->SetStats(0);
                hwidth->Draw("pe");
                hwidth->Write("width");
                TLegend *widthleg = new TLegend(0.2, 0.75, 0.4, 0.85);
                SetLegendStyle(widthleg);
                widthleg->SetTextSize(txtsize);
                // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
                TLine *line2 = new TLine(hwidth->GetXaxis()->GetXmin(), 0.047, hwidth->GetXaxis()->GetXmax(), 0.047);
                line2->SetLineStyle(2);
                line2->SetLineColor(2);
                line2->SetLineWidth(3);
                line2->Draw();
                widthleg->AddEntry(line, "PDG Width", "l");
                widthleg->SetFillStyle(0);
                widthleg->Draw();
                csig->SaveAs(outputfolder_mult + ("/width_pt." + outputtype).c_str());
                csig->Clear();

                // // // Yield vs pT (integral method)
                SetHistoQA(hintegral_yield);
                hintegral_yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                hintegral_yield->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
                gPad->SetLogy(1);
                // hintegral_yield->GetXaxis()->SetRangeUser(-0.1, 15.2);
                hintegral_yield->SetStats(0);
                hintegral_yield->Draw("pe");
                hintegral_yield->Write("yield_integral");
                // hintegral_yield->Write("yield");
                TLegend *legyield = new TLegend(0.8, 0.8, 0.91, 0.9);
                SetLegendStyle(legyield);
                legyield->SetTextSize(txtsize);
                // legyield->AddEntry(hYieldpar, "pbpb 5.36 TeV");
                // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
                legyield->Draw();
                csig->SaveAs(outputfolder_mult + ("/yield_integral." + outputtype).c_str());
                csig->Clear();

                // Yield vs pT (bin counting method)
                hYbincount->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                hYbincount->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
                SetHistoQA(hYbincount);
                hYbincount->SetStats(0);
                hYbincount->Draw("pe");
                hYbincount->Write("yield_bincount");
                csig->SaveAs(outputfolder_mult + ("/yield_bincount." + outputtype).c_str());
                csig->Close();

                hsignificance->Clear();
                hChiSquare->Clear();
                hmass->Clear();
                hwidth->Clear();
                hintegral_yield->Clear();
                hYbincount->Clear();
            }
            filecmp->cd();
        }
    }

    // Stop the stopwatch
    timer.Stop();

    // Get the elapsed time
    Double_t realTime = timer.RealTime(); // Wall clock time in seconds
    Double_t cpuTime = timer.CpuTime();   // CPU time used in seconds

    // Print the elapsed times
    std::cout << "Real time elapsed: " << realTime << " seconds" << std::endl;
    std::cout << "CPU time used: " << cpuTime << " seconds" << std::endl;

    if (makeQAplots)
    {
        gStyle->SetOptStat(0);
        TH1F *hDCAxy = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hDcaxy", multpath.c_str()));
        TH1F *hDCAz = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hDcaz", multpath.c_str()));
        TH1F *hMult = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hMultiplicity", multpath.c_str()));
        TH1F *hvz = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hVertexZRec", multpath.c_str()));
        TH1F *hoccupancy = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hOccupancy", multpath.c_str()));
        if (hDCAxy == nullptr || hDCAz == nullptr || hMult == nullptr || hvz == nullptr || hoccupancy == nullptr)
        {
            cerr << "Event selection histograms not found!!!!!!!!!!!!" << endl;
            return;
        }

        TH2F *hNsigmaTPCTOFKaon = (TH2F *)fInputFile->Get(Form("%s/hPID/After/hNsigma_TPC_TOF_Ka_after", multpath.c_str()));
        TH2F *hNsigmaTPCTOFPion = (TH2F *)fInputFile->Get(Form("%s/hPID/After/hNsigma_TPC_TOF_Pi_after", multpath.c_str()));
        TH2F *hNsigmaTPCKaon = (TH2F *)fInputFile->Get(Form("%s/hPID/After/hNsigmaKaonTPC_after", multpath.c_str()));
        TH2F *hNsigmaTPCPion = (TH2F *)fInputFile->Get(Form("%s/hPID/After/hNsigmaPionTPC_after", multpath.c_str()));
        TH2F *hNsigmaTOFKaon = (TH2F *)fInputFile->Get(Form("%s/hPID/After/hNsigmaKaonTOF_after", multpath.c_str()));
        TH2F *hNsigmaTOFPion = (TH2F *)fInputFile->Get(Form("%s/hPID/After/hNsigmaPionTOF_after", multpath.c_str()));
        if (hNsigmaTPCTOFKaon == nullptr || hNsigmaTPCTOFPion == nullptr || hNsigmaTPCKaon == nullptr || hNsigmaTPCPion == nullptr || hNsigmaTOFKaon == nullptr || hNsigmaTOFPion == nullptr)
        {
            cerr << "PID histograms after selection not found!!!!!!!!!!!!" << endl;
            return;
        }

        TH1F *hNsigmaTOFKaon_neg = (TH1F *)fInputFile->Get(Form("%s/hPID/Before/h1PID_TOF_neg_kaon", multpath.c_str()));
        TH1F *hNsigmaTOFPion_neg = (TH1F *)fInputFile->Get(Form("%s/hPID/Before/h1PID_TOF_neg_pion", multpath.c_str()));
        TH1F *hNsigmaTOFKaon_pos = (TH1F *)fInputFile->Get(Form("%s/hPID/Before/h1PID_TOF_pos_kaon", multpath.c_str()));
        TH1F *hNsigmaTOFPion_pos = (TH1F *)fInputFile->Get(Form("%s/hPID/Before/h1PID_TOF_pos_pion", multpath.c_str()));
        TH1F *hNsigmaTPCKaon_neg = (TH1F *)fInputFile->Get(Form("%s/hPID/Before/h1PID_TPC_neg_kaon", multpath.c_str()));
        TH1F *hNsigmaTPCPion_neg = (TH1F *)fInputFile->Get(Form("%s/hPID/Before/h1PID_TPC_neg_pion", multpath.c_str()));
        TH1F *hNsigmaTPCKaon_pos = (TH1F *)fInputFile->Get(Form("%s/hPID/Before/h1PID_TPC_pos_kaon", multpath.c_str()));
        TH1F *hNsigmaTPCPion_pos = (TH1F *)fInputFile->Get(Form("%s/hPID/Before/h1PID_TPC_pos_pion", multpath.c_str()));
        if (hNsigmaTOFKaon_neg == nullptr || hNsigmaTOFPion_neg == nullptr || hNsigmaTOFKaon_pos == nullptr || hNsigmaTOFPion_pos == nullptr || hNsigmaTPCKaon_neg == nullptr || hNsigmaTPCPion_neg == nullptr || hNsigmaTPCKaon_pos == nullptr || hNsigmaTPCPion_pos == nullptr)
        {
            cerr << "PID histograms before selection not found!!!!!!!!!!!!" << endl;
            return;
        }

        TCanvas *cDCAxy = new TCanvas("cDCAxy", "DCAxy", 720, 720);
        SetCanvasStyle(cDCAxy, 0.14, 0.03, 0.06, 0.14);
        SetHistoQA(hDCAxy);
        gPad->SetLogy();
        hDCAxy->GetXaxis()->SetTitle("DCA_{xy} (cm)");
        hDCAxy->GetYaxis()->SetTitle("Counts");
        hDCAxy->Draw();
        cDCAxy->SaveAs(output_QA_folder + ("/DCAxy." + outputtype).c_str());

        TCanvas *cDCAz = new TCanvas("cDCAz", "DCAz", 720, 720);
        SetCanvasStyle(cDCAz, 0.14, 0.03, 0.06, 0.14);
        gPad->SetLogy();
        SetHistoQA(hDCAz);
        hDCAz->GetXaxis()->SetTitle("DCA_{z} (cm)");
        hDCAz->GetYaxis()->SetTitle("Counts");
        hDCAz->Draw();
        cDCAz->SaveAs(output_QA_folder + ("/DCAz." + outputtype).c_str());

        gPad->SetLogy(0);
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

        TCanvas *cOccupancy = new TCanvas("cOccupancy", "Occupancy", 720, 720);
        SetCanvasStyle(cOccupancy, 0.14, 0.03, 0.06, 0.14);
        gPad->SetLogy();
        SetHistoQA(hoccupancy);
        hoccupancy->GetXaxis()->SetTitle("Occupancy");
        hoccupancy->GetYaxis()->SetTitle("Counts");
        hoccupancy->GetXaxis()->SetRangeUser(0, 2200);
        hoccupancy->GetXaxis()->SetNdivisions(505);
        hoccupancy->Draw();
        cOccupancy->SaveAs(output_QA_folder + ("/Occupancy." + outputtype).c_str());

        gPad->SetLogy(0);
        TCanvas *cNsigmaTPCTOFKaon = new TCanvas("cNsigmaTPCTOFKaon", "Nsigma TPC TOF Kaon", 720, 720);
        SetCanvasStyle(cNsigmaTPCTOFKaon, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTPCTOFKaon);
        hNsigmaTPCTOFKaon->GetXaxis()->SetTitle("n#sigma_{TPC}");
        hNsigmaTPCTOFKaon->GetYaxis()->SetTitle("n#sigma_{TOF}");
        hNsigmaTPCTOFKaon->Draw("colz");
        cNsigmaTPCTOFKaon->SaveAs(output_QA_folder + ("/NsigmaTPCTOFKaon." + outputtype).c_str());

        TCanvas *cNsigmaTPCTOFPion = new TCanvas("cNsigmaTPCTOFPion", "Nsigma TPC TOF Pion", 720, 720);
        SetCanvasStyle(cNsigmaTPCTOFPion, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTPCTOFPion);
        hNsigmaTPCTOFPion->GetXaxis()->SetTitle("n#sigma_{TPC}");
        hNsigmaTPCTOFPion->GetYaxis()->SetTitle("n#sigma_{TOF}");
        hNsigmaTPCTOFPion->Draw("colz");
        cNsigmaTPCTOFPion->SaveAs(output_QA_folder + ("/NsigmaTPCTOFPion." + outputtype).c_str());

        TCanvas *cNsigmaTPCKaon = new TCanvas("cNsigmaTPCKaon", "Nsigma TPC Kaon", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTPCKaon);
        hNsigmaTPCKaon->GetYaxis()->SetTitle("n#sigma_{TPC}");
        hNsigmaTPCKaon->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hNsigmaTPCKaon->Draw("colz");
        cNsigmaTPCKaon->SaveAs(output_QA_folder + ("/NsigmaTPCKaon." + outputtype).c_str());

        TCanvas *cNsigmaTPCPion = new TCanvas("cNsigmaTPCPion", "Nsigma TPC Pion", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTPCPion);
        hNsigmaTPCPion->GetYaxis()->SetTitle("n#sigma_{TPC}");
        hNsigmaTPCPion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hNsigmaTPCPion->Draw("colz");
        cNsigmaTPCPion->SaveAs(output_QA_folder + ("/NsigmaTPCPion." + outputtype).c_str());

        TCanvas *cNsigmaTOFKaon = new TCanvas("cNsigmaTOFKaon", "Nsigma TOF Kaon", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTOFKaon);
        hNsigmaTOFKaon->GetYaxis()->SetTitle("n#sigma_{TOF}");
        hNsigmaTOFKaon->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hNsigmaTOFKaon->Draw("colz");
        cNsigmaTOFKaon->SaveAs(output_QA_folder + ("/NsigmaTOFKaon." + outputtype).c_str());

        TCanvas *cNsigmaTOFPion = new TCanvas("cNsigmaTOFPion", "Nsigma TOF Pion", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(hNsigmaTOFPion);
        hNsigmaTOFPion->GetYaxis()->SetTitle("n#sigma_{TOF}");
        hNsigmaTOFPion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hNsigmaTOFPion->Draw("colz");
        cNsigmaTOFPion->SaveAs(output_QA_folder + ("/NsigmaTOFPion." + outputtype).c_str());

        TCanvas *cNsigmaTOFKaon_neg = new TCanvas("cNsigmaTOFKaon_neg", "Nsigma TOF Kaon Neg", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon_neg, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTOFKaon_neg);
        hNsigmaTOFKaon_neg->GetYaxis()->SetTitle("Counts");
        hNsigmaTOFKaon_neg->GetXaxis()->SetTitle("n#sigma_{TOF} (K^{-})");
        hNsigmaTOFKaon_neg->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTOFKaon_neg->Draw();
        TLine *lineverticalx0 = new TLine(0, hNsigmaTOFKaon_neg->GetYaxis()->GetXmin(), 0, hNsigmaTOFKaon_neg->GetYaxis()->GetXmax());
        lineverticalx0->SetLineStyle(2);
        lineverticalx0->SetLineColor(2);
        lineverticalx0->SetLineWidth(2);
        lineverticalx0->Draw();
        cNsigmaTOFKaon_neg->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_neg." + outputtype).c_str());

        TCanvas *cNsigmaTOFPion_neg = new TCanvas("cNsigmaTOFPion_neg", "Nsigma TOF Pion Neg", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion_neg, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTOFPion_neg);
        hNsigmaTOFPion_neg->GetYaxis()->SetTitle("Counts");
        hNsigmaTOFPion_neg->GetXaxis()->SetTitle("n#sigma_{TOF} (#pi^{-})");
        hNsigmaTOFPion_neg->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTOFPion_neg->Draw();
        cNsigmaTOFPion_neg->SaveAs(output_QA_folder + ("/NsigmaTOFPion_neg." + outputtype).c_str());

        TCanvas *cNsigmaTOFKaon_pos = new TCanvas("cNsigmaTOFKaon_pos", "Nsigma TOF Kaon Pos", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon_pos, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTOFKaon_pos);
        hNsigmaTOFKaon_pos->GetYaxis()->SetTitle("Counts");
        hNsigmaTOFKaon_pos->GetXaxis()->SetTitle("n#sigma_{TOF} (K^{+})");
        hNsigmaTOFKaon_pos->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTOFKaon_pos->Draw();
        cNsigmaTOFKaon_pos->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_pos." + outputtype).c_str());

        TCanvas *cNsigmaTOFPion_pos = new TCanvas("cNsigmaTOFPion_pos", "Nsigma TOF Pion Pos", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion_pos, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTOFPion_pos);
        hNsigmaTOFPion_pos->GetYaxis()->SetTitle("Counts");
        hNsigmaTOFPion_pos->GetXaxis()->SetTitle("n#sigma_{TOF} (#pi^{+})");
        hNsigmaTOFPion_pos->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTOFPion_pos->Draw();
        cNsigmaTOFPion_pos->SaveAs(output_QA_folder + ("/NsigmaTOFPion_pos." + outputtype).c_str());

        TCanvas *cNsigmaTPCKaon_neg = new TCanvas("cNsigmaTPCKaon_neg", "Nsigma TPC Kaon Neg", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon_neg, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTPCKaon_neg);
        hNsigmaTPCKaon_neg->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCKaon_neg->GetXaxis()->SetTitle("n#sigma_{TPC} (K^{-})");
        hNsigmaTPCKaon_neg->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTPCKaon_neg->Draw();
        cNsigmaTPCKaon_neg->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_neg." + outputtype).c_str());

        TCanvas *cNsigmaTPCPion_neg = new TCanvas("cNsigmaTPCPion_neg", "Nsigma TPC Pion Neg", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion_neg, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTPCPion_neg);
        hNsigmaTPCPion_neg->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCPion_neg->GetXaxis()->SetTitle("n#sigma_{TPC} (#pi^{-})");
        hNsigmaTPCPion_neg->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTPCPion_neg->Draw();
        cNsigmaTPCPion_neg->SaveAs(output_QA_folder + ("/NsigmaTPCPion_neg." + outputtype).c_str());

        TCanvas *cNsigmaTPCKaon_pos = new TCanvas("cNsigmaTPCKaon_pos", "Nsigma TPC Kaon Pos", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon_pos, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTPCKaon_pos);
        hNsigmaTPCKaon_pos->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCKaon_pos->GetXaxis()->SetTitle("n#sigma_{TPC} (K^{+})");
        hNsigmaTPCKaon_pos->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTPCKaon_pos->Draw();
        cNsigmaTPCKaon_pos->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_pos." + outputtype).c_str());

        TCanvas *cNsigmaTPCPion_pos = new TCanvas("cNsigmaTPCPion_pos", "Nsigma TPC Pion Pos", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion_pos, 0.14, 0.05, 0.06, 0.14);
        SetHistoQA(hNsigmaTPCPion_pos);
        hNsigmaTPCPion_pos->GetYaxis()->SetTitle("Counts");
        hNsigmaTPCPion_pos->GetXaxis()->SetTitle("n#sigma_{TPC} (#pi^{+})");
        hNsigmaTPCPion_pos->GetXaxis()->SetRangeUser(-3, 3);
        hNsigmaTPCPion_pos->Draw();
        cNsigmaTPCPion_pos->SaveAs(output_QA_folder + ("/NsigmaTPCPion_pos." + outputtype).c_str());
    }
}
