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

////Initial parameters from template fit
float mean_CB[] = {0.7875, 0.7527, 0.7678, 0.7684, 0.7686, 0.769, 0.7698, 0.7709, 0.7706, 0.7715, 0.7712, 0.7712, 0.7715, 0.7721};
float sigma_CB[] = {0.06119, 0.07456, 0.08454, 0.08353, 0.08075, 0.07841, 0.07658, 0.07632, 0.07539, 0.07474, 0.07387, 0.07425, 0.07425, 0.07309};
float alpha_CB[] = {1.2, 0.8809, 0.5499, 0.5506, 0.6135, 0.7388, 0.8244, 0.8935, 0.9161, 0.9677, 0.9634, 0.9524, 0.9648, 1.031};
float n_CB[] = {146.8, 139.5, 127.8, 1.23, 0.6261, 0.3242, 0.2302, 0.1739, 0.1594, 0.1275, 0.1449, 0.1551, 0.1463, 0.1128};

void kstar_sparse()
{
    TStopwatch timer;
    timer.Start();
    //*************************** change here ***************************************
    // TString sysVars[] = {"", "Norm1", "Norm2", "FitRange1", "FitRange2", "WidthFree"};
    TString sysVars[] = {""};
    int nSysVars = sizeof(sysVars) / sizeof(sysVars[0]);
    const string kResBkg = "MIX";
    // const string kResBkg = "LIKE";
    // const string kResBkg = "ROTATED";
    const string kbkg = "pol3";
    string outputtype = "pdf";     // pdf, eps
    const bool save_bkg_plots = 1; // save background plots
    const float txtsize = 0.045;   // text size in the plots
    bool makeQAplots = false;
    bool makeallpTplots = true; // make all pT plots
    bool calcInvMass = true;
    bool isINEL = false;

    int colors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 1, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

    TCanvas *cgrid1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid_bkg1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid2 = nullptr;
    TCanvas *cgrid_bkg2 = nullptr;
    // if (Npt > 16)
    if (Npt > 9)
    {
        cgrid2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
        cgrid_bkg2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    }

    // some initializations ********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************
    const string resBkgFolder = (kResBkg == "MIX") ? "" : "/" + kResBkg;
    const string kbkgFolder = (kbkg == "pol3") ? "" : "/" + kbkg;

    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.06);
    t2->SetTextFont(42);

    int sizeOfpTbins = sizeof(pT_bins) / sizeof(pT_bins[0]) - 1;
    if (sizeOfpTbins != Npt)
    {
        cerr << "Error: size of pT_bins array does not match Npt" << endl;
        return;
    }

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
        cgrid_bkg1->Divide(kcanvasdivide[0], kcanvasdivide[1]);
        // if (Npt > 16)
        if (Npt > 9)
        {
            cgrid2->Divide(kcanvasdivide[0], kcanvasdivide[1]);
            cgrid_bkg2->Divide(kcanvasdivide[0], kcanvasdivide[1]);
        }
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

    // float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    float mult_classes[] = {0.0};
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
    //********************************************************************************
    //*************************Create folders********************************************
    TString output_QA_folder = kSignalOutput + "/" + kfoldername + resBkgFolder + kbkgFolder + "/QA";
    if (gSystem->mkdir(output_QA_folder, kTRUE))
    {
        std::cout << "Folder " << output_QA_folder << " created successfully." << std::endl;
    }

    for (int ivar = 0; ivar < nSysVars; ivar++)
    // for (int ivar = 0; ivar < 1; ivar++)
    {
        if (nSysVars > 1 && (kResBkg != "MIX" || kbkg != "pol3"))
        {
            cerr << "Error: Systematic variations are only implemented for MIX background." << endl;
            return;
        }

        //***********************************************************************************
        if (calcInvMass)
        {
            TString outputRootDir = (ivar == 0)
                                        ? koutputfolder + resBkgFolder + kbkgFolder
                                        : koutputfolder + "/" + sysVars[ivar] + kbkgFolder;
            if (gSystem->mkdir(outputRootDir, kTRUE))
            {
                std::cout << "Folder " << outputRootDir << " created successfully." << std::endl;
            }

            //// Commented it out to store separate root files for each multiplicity bin (for systematics)
            // TFile *filecmp;
            // if (ivar == 0)
            // {
            //     filecmp = (isINEL) ? new TFile((koutputfolder + resBkgFolder + kbkgFolder + "/yield_INEL.root").c_str(), "RECREATE") : new TFile((koutputfolder + resBkgFolder + kbkgFolder + "/yield.root").c_str(), "RECREATE");
            // }
            // else
            // {
            //     filecmp = (isINEL) ? new TFile((koutputfolder + "/" + sysVars[ivar].Data() + kbkgFolder + "/yield_INEL.root").c_str(), "RECREATE") : new TFile((koutputfolder + "/" + sysVars[ivar].Data() + kbkgFolder + "/yield.root").c_str(), "RECREATE");
            // }
        }

        for (int imult = 0; imult < nmultbins + 1; imult++)
        // for (int imult = 0; imult < 1; imult++)
        {
            if (isINEL && imult != 0)
                break;

            // basic checks
            if (kNormRangepT.size() < Npt || kFitRange.size() < Npt || kRebin.size() < Npt)
            {
                cerr << "Error: kNormRangepT, kFitRange, or kRebin arrays are not initialized for all pT bins." << endl;
                return;
            }
            if (kNormRangepT.size() != kFitRange.size() || kNormRangepT.size() != kRebin.size())
            {
                cout << "!!!!!!!!!!!!!!!Warning !!!!!!!!!!" << endl;
                cout << "!!!!!!!!!!!!!!!Warning !!!!!!!!!!" << endl;
                cout << "Error: kNormRangepT, kFitRange, and kRebin arrays must have the same size." << endl;
                cout << "!!!!!!!!!!!!!!!Warning !!!!!!!!!!" << endl;
                cout << "!!!!!!!!!!!!!!!Warning !!!!!!!!!!" << endl;
                // return;
            }
            //**************Invariant mass histograms for sig+bkg and mixed event bg******************
            int multlow, multhigh;

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

            // TFile *filecmp = new TFile((koutputfolder + "/" + sysVars[ivar].Data() + kbkgFolder + Form("/yield_%d_%d.root", multlow, multhigh)).c_str(), "RECREATE");

            TFile *filecmp;
            if (ivar == 0)
            {
                filecmp = new TFile((koutputfolder + resBkgFolder + kbkgFolder + + Form("/yield_%d_%d.root", multlow, multhigh)).c_str(), "RECREATE");
            }
            else
            {
                filecmp = new TFile((koutputfolder + "/" + sysVars[ivar].Data() + kbkgFolder + + Form("/yield_%d_%d.root", multlow, multhigh)).c_str(), "RECREATE");
            }

            double Event = hmult->Integral(hmult->GetXaxis()->FindBin(multlow + 1e-5), hmult->GetXaxis()->FindBin(multhigh - 1e-5));
            cout << "Event in mult bin " << imult << " is " << Event << endl;
            TString outputfolder_mult;
            if (ivar == 0)
            {
                outputfolder_mult = kSignalOutput + "/" + kfoldername + resBkgFolder + kbkgFolder + Form("/mult_%d-%d", multlow, multhigh);
            }
            else
            {
                outputfolder_mult = kSignalOutput + "/" + kfoldername + "/" + sysVars[ivar] + kbkgFolder + Form("/mult_%d-%d", multlow, multhigh);
            }
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

                float prevWidth = widthpdg;

                for (Int_t ip = pt_start; ip < pt_end; ip++) // start pt bin loop
                {
                    rebin_value = kRebin[ip][imult]; // rebinning value for the multiplicity bin
                    // rebin_value = 2; // for medium dataset (temporarily set to 4)

                    double lowfitrange = kFitRange[ip][0];
                    double highfitrange = kFitRange[ip][1];
                    if (sysVars[ivar] == "FitRange1")
                    {
                        lowfitrange -= 0.02;
                        highfitrange -= 0.02;
                    }
                    else if (sysVars[ivar] == "FitRange2")
                    {
                        lowfitrange += 0.02;
                        highfitrange += 0.02;
                    }

                    //// For multiplicity estimator FV0A
                    // if (imult == 7 || imult == 8)
                    // {
                    //     if (ip == pt_end - 1)
                    //     {
                    //         lowfitrange = 0.76;
                    //         highfitrange = 1.04;
                    //     }
                    // }

                    lowpt = pT_bins[ip];
                    highpt = pT_bins[ip + 1];
                    int lbin = fHistNum->GetAxis(1)->FindBin(lowpt + 1e-5);
                    int hbin = fHistNum->GetAxis(1)->FindBin(highpt - 1e-5);

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

                    int lbinmult = fHistNum->GetAxis(0)->FindBin(multlow + 1e-5);
                    int hbinmult = fHistNum->GetAxis(0)->FindBin(multhigh - 1e-5);

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
                    float normRangeLow = kNormRangepT[ip][0];
                    float normRangeHigh = kNormRangepT[ip][1];
                    if (sysVars[ivar] == "Norm1")
                    {
                        normRangeLow = kNormRangepT_sysVar1[ip][0];
                        normRangeHigh = kNormRangepT_sysVar1[ip][1];
                    }
                    else if (sysVars[ivar] == "Norm2")
                    {
                        normRangeLow = kNormRangepT_sysVar2[ip][0];
                        normRangeHigh = kNormRangepT_sysVar2[ip][1];
                    }

                    if (kResBkg == "MIX" || kResBkg == "ROTATED")
                    {
                        TH1D *bkgclonetemp = (kResBkg == "MIX") ? (TH1D *)fHistBkg[ip]->Clone() : (TH1D *)fHistRotated1D[ip]->Clone();

                        sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(normRangeLow), fHistTotal[ip]->GetXaxis()->FindBin(normRangeHigh)));
                        bkg_integral = (bkgclonetemp->Integral(bkgclonetemp->GetXaxis()->FindBin(normRangeLow), bkgclonetemp->GetXaxis()->FindBin(normRangeHigh)));
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
                        // sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(normRangeLow), fHistTotal[ip]->GetXaxis()->FindBin(normRangeHigh)));
                        // bkg_integral = (fHistbkgLS[ip]->Integral(fHistbkgLS[ip]->GetXaxis()->FindBin(normRangeLow), fHistbkgLS[ip]->GetXaxis()->FindBin(normRangeHigh)));
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

                    TF1 *fitFcn, *fitFcn1, *fitFcn2;

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
                    else if (kbkg == "CBpol2")
                    {
                        fitFcn = new TF1("fitfunc", BWCBpol2, lowfitrange, highfitrange, 11);
                        fitFcn1 = new TF1("fitfunc1", CBRightpol2, lowfitrange, highfitrange, 8);
                    }
                    else if (kbkg == "CBpol3")
                    {
                        fitFcn = new TF1("fitfunc", BWCBpol3, lowfitrange, highfitrange, 12);
                        fitFcn1 = new TF1("fitfunc1", CBRightpol3, lowfitrange, highfitrange, 9);
                    }

                    fitFcn2 = new TF1("fitFcn2", BW, lowfitrange, highfitrange, 3); // only signal

                    fitFcn->SetParameter(0, masspdg); // mass
                    // fitFcn->SetParLimits(0, 0.7, 0.98); // mass limits
                    // fitFcn->SetParLimits(0, 0.8, 0.98); // mass limits

                    if (imult == 1 && ip < 3)
                    {
                        fitFcn->SetParLimits(0, 0.885, 0.888); // Mass
                    }
                    else if (imult == 1 && ip > Npt - 4)
                    {
                        fitFcn->SetParLimits(0, 0.893, 0.9); // Mass
                    }
                    else if (imult == nmultbins && ip > Npt - 4)
                    {
                        fitFcn->SetParLimits(0, 0.893, 0.9); // Mass
                    }
                    else if (imult == 1 && ip == 10)
                    {
                        fitFcn->SetParLimits(0, 0.89, 0.91); // Mass
                    }
                    else
                    {
                        fitFcn->SetParLimits(0, 0.80, 0.98); // Mass
                    }

                    fitFcn->SetParameter(1, widthpdg); // width
                    if (sysVars[ivar] != "WidthFree")
                        fitFcn->FixParameter(1, widthpdg); // width
                    else
                        fitFcn->SetParLimits(1, 0.044, 0.054); // width
                    // fitFcn->SetParLimits(1, prevWidth, prevWidth + 0.003); // For MC closure test in min bias

                    fitFcn->SetParameter(2, 10e4);     // yield
                    fitFcn->SetParLimits(2, 0.0, 1e8); // Yield2

                    // // //pol3 parameters
                    // fitFcn->SetParameter(3, -1e6);
                    // fitFcn->SetParameter(4, 1e6);
                    // fitFcn->SetParameter(5, -1e6);
                    // fitFcn->SetParameter(6, 1e6);

                    ////CBpol2 parameters
                    if (kbkg == "CBpol2" || kbkg == "CBpol3")
                    {
                        fitFcn->SetParameter(3, 1);
                        fitFcn->FixParameter(4, mean_CB[ip]);
                        fitFcn->FixParameter(5, sigma_CB[ip]);
                        fitFcn->FixParameter(6, alpha_CB[ip]);
                        fitFcn->FixParameter(7, n_CB[ip]);
                        fitFcn->SetParameter(8, 1e7);
                        fitFcn->SetParameter(9, 1e7);
                        fitFcn->SetParameter(10, 1e7);
                        if (kbkg == "CBpol3")
                            fitFcn->SetParameter(11, 1e7);
                    }

                    fitFcn->SetParNames("Mass", "Width", "Yield", "A", "B", "C", "D");
                    // Redirect standard output to /dev/null
                    // int old_stdout = dup(1);
                    // freopen("/dev/null", "w", stdout);

                    r = hfsig->Fit(fitFcn, "REBMS0+"); // signal after bkg subtraction

                    prevWidth = fitFcn->GetParameter(1);

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
                    fitFcn2_plusm->FixParameter(1, widthpdg);
                    fitFcn2_plusm->FixParameter(2, Yield[ip]);

                    fitFcn2_minusm->FixParameter(0, Mass[ip] - ErrorMass[ip]);
                    fitFcn2_minusm->FixParameter(1, widthpdg);
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

                    cout << "Yield from bin counting is " << Total_Ybincounting << " +/- " << Final_pro_error << endl;

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
                    cout << "Total yield from function integration is " << yieldcalc << " +/- " << yielderror << endl;

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
                    hfsig->GetYaxis()->SetMaxDigits(3);
                    hfsig->GetYaxis()->CenterTitle(1);
                    hfsig->GetYaxis()->SetTitleOffset(1.45);
                    hfsig->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidth_file * 1000));

                    SetHistoQA(hfbkg);
                    hfbkg->SetLineColor(kRed);
                    hfbkg->SetMarkerColor(kRed);
                    hfbkg->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
                    // hfsig->GetYaxis()->SetMaxDigits(3);
                    hfbkg->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidth_file * 1000));

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
                    hfsig->SetMaximum(hfsig->GetMaximum() * 1.3);
                    fitFcn->SetLineWidth(2);
                    fitFcn1->SetLineWidth(2);
                    fitFcn2->SetLineWidth(2);
                    hfsig->GetXaxis()->SetRangeUser(0.7, 1.11);
                    hfsig->SetMarkerSize(1.0);
                    if (ip != 0 && ip != 1 && ip != 2)
                        hfsig->SetMinimum(0);
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
                    if (!multipanel_plots)
                        pag->Draw();

                    TLegend *pag2 = new TLegend(0.2, 0.7, 0.45, 0.9);
                    // TLegend *pag2 = new TLegend(0.7, 0.65, 0.9, 0.9);
                    pag2->SetBorderSize(0);
                    pag2->SetTextFont(42);
                    pag2->SetTextSize(0.035);
                    pag2->SetFillStyle(0);
                    pag2->AddEntry(fitFcn, "BW+pol3", "l");
                    pag2->AddEntry(fitFcn1, "BW", "l");
                    pag2->AddEntry(fitFcn2, "pol3", "l");
                    pag2->Draw();
                    if (multipanel_plots)
                        t2->DrawLatex(0.27, 0.95, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip], pT_bins[ip + 1]));

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
                        if (hbkg_nopeak->GetBinCenter(i + 1) < normRangeLow || hbkg_nopeak->GetBinCenter(i + 1) > normRangeHigh)
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
                    fHistTotal[ip]->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidth_file * 1000));

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
                    if (!multipanel_plots)
                        pag->Draw();

                    TLegend *pag3;
                    pag3 = (!multipanel_plots) ? new TLegend(0.7, 0.65, 0.9, 0.9) : new TLegend(0.2, 0.7, 0.45, 0.9);
                    pag3->SetBorderSize(0);
                    pag3->SetTextFont(42);
                    pag3->SetTextSize(0.035);
                    pag3->SetFillStyle(0);
                    pag3->AddEntry(fHistTotal[ip], "Signal + Bkg", "p");
                    pag3->AddEntry(hfbkg, "Bkg", "p");
                    pag3->AddEntry(hbkg_nopeak, "Norm. region", "f");
                    pag3->Draw();

                    if (multipanel_plots)
                        t2->DrawLatex(0.27, 0.95, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip], pT_bins[ip + 1]));

                    if (multipanel_plots == 0 && save_plots == 1 && save_bkg_plots == 1)
                        cSigbkg[ip]->SaveAs(Form(outputfolder_mult + ("/hsigbkg_pt%d." + outputtype).c_str(), ip + 1));
                    // cSigbkg[ip]->Close();

                    // ////////////////////////////////////////////////////////////////////////

                } // pt loop ends
                if (multipanel_plots == 1 && save_plots == 1)
                {
                    cgrid1->SaveAs(outputfolder_mult + (Form("/grid1_mult%d_%d.", multlow, multhigh) + outputtype).c_str());
                    cgrid_bkg1->SaveAs(outputfolder_mult + (Form("/gridBkg1_mult%d_%d.", multlow, multhigh) + outputtype).c_str());

                    if (Npt >= klowerpad * kupperpad)
                    {
                        cgrid2->SaveAs(outputfolder_mult + (Form("/grid2_mult%d_%d.", multlow, multhigh) + outputtype).c_str());
                        cgrid_bkg2->SaveAs(outputfolder_mult + (Form("/gridBkg2_mult%d_%d.", multlow, multhigh) + outputtype).c_str());
                    }
                }
                if (multipanel_plots == 0)
                {
                    cgrid1->Close();
                    cgrid_bkg1->Close();

                    if (Npt > 9)
                    {
                        cgrid2->Close();
                        cgrid_bkg2->Close();
                    }
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
                    // csig->SaveAs(outputfolder_mult + ("/chi." + outputtype).c_str());
                    csig->SaveAs(outputfolder_mult + "/chi.png");
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
                    TLine *line = new TLine(hmass->GetXaxis()->GetXmin(), masspdg, hmass->GetXaxis()->GetXmax(), masspdg);
                    line->SetLineStyle(2);
                    line->SetLineColor(2);
                    line->SetLineWidth(3);
                    line->Draw();
                    massleg->AddEntry(line, "PDG Mass", "l");
                    massleg->Draw("l");
                    // csig->SaveAs(outputfolder_mult + ("/mass." + outputtype).c_str());
                    csig->SaveAs(outputfolder_mult + "/mass.png");
                    csig->Clear();

                    // // // Width vs pT
                    hwidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                    hwidth->GetYaxis()->SetTitle("Width (GeV)");
                    SetHistoQA(hwidth);
                    hwidth->GetYaxis()->SetRangeUser(0.04, 0.06);
                    hwidth->SetMinimum(0);
                    hwidth->SetStats(0);
                    hwidth->Draw("pe");
                    hwidth->Write("width");
                    TLegend *widthleg = new TLegend(0.2, 0.75, 0.4, 0.85);
                    SetLegendStyle(widthleg);
                    widthleg->SetTextSize(txtsize);
                    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
                    TLine *line2 = new TLine(hwidth->GetXaxis()->GetXmin(), widthpdg, hwidth->GetXaxis()->GetXmax(), widthpdg);
                    line2->SetLineStyle(2);
                    line2->SetLineColor(2);
                    line2->SetLineWidth(3);
                    line2->Draw();
                    widthleg->AddEntry(line2, "PDG Width", "l");
                    widthleg->SetFillStyle(0);
                    widthleg->Draw();
                    // csig->SaveAs(outputfolder_mult + ("/width_pt." + outputtype).c_str());
                    csig->SaveAs(outputfolder_mult + "/width_pt.png");
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
                    // csig->SaveAs(outputfolder_mult + ("/yield_integral." + outputtype).c_str());
                    csig->SaveAs(outputfolder_mult + "/yield_integral.png");
                    csig->Clear();

                    // Yield vs pT (bin counting method)
                    hYbincount->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                    hYbincount->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
                    SetHistoQA(hYbincount);
                    hYbincount->SetStats(0);
                    hYbincount->Draw("pe");
                    hYbincount->Write("yield_bincount");
                    // csig->SaveAs(outputfolder_mult + ("/yield_bincount." + outputtype).c_str());
                    csig->SaveAs(outputfolder_mult + "/yield_bincount.png");
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
        cout << "============= End of the code =============" << endl;
        cout << "Data file used: " << kDataFilename.c_str() << endl;
        cout << "Selection used: " << (isINEL ? "INEL" : "INEL > 0") << endl;
        cout << "Format of output plots: " << outputtype.c_str() << endl;
        cout << "Residual background function: " << kbkg.c_str() << endl;
        cout << "Number of pT bins: " << Npt << endl;
        (multipanel_plots) ? cout << "Canvas: " << klowerpad << "x" << kupperpad << " multi-panel" << endl : cout << "Canvas: Single panel plots" << endl;
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
        TH3F *hDCAxy3D = (TH3F *)fInputFile->Get(Form("%s/eventSelection/hDcaxy_cent_pt", multpath.c_str()));
        TH3F *hDCAz3D = (TH3F *)fInputFile->Get(Form("%s/eventSelection/hDcaz_cent_pt", multpath.c_str()));
        TH1F *hMult = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hMultiplicity", multpath.c_str()));
        TH1F *hvz = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hVertexZRec", multpath.c_str()));
        TH1F *hoccupancy = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hOccupancy", multpath.c_str()));
        TH1F *hEventCut = (TH1F *)fInputFile->Get(Form("%s/eventSelection/hEventCut", multpath.c_str()));
        TH1F *htracksData = (TH1F *)fInputFile->Get(Form("%s/eventSelection/tracksCheckData", multpath.c_str()));
        TH1D *hDCAxy = hDCAxy3D->ProjectionX("hDCAxy", -1, -1, -1, -1);
        TH1D *hDCAz = hDCAz3D->ProjectionX("hDCAz", -1, -1, -1, -1);
        if (hDCAxy == nullptr || hDCAz == nullptr || hMult == nullptr || hvz == nullptr || hoccupancy == nullptr || hEventCut == nullptr || htracksData == nullptr)
        {
            cerr << "Event selection histograms not found!!!!!!!!!!!!" << endl;
            return;
        }

        // For PID QA plots, another separate code is present. In this code, we need to replace TH2F with TH3F to plot PID.

        TH3F *hNsigmaTPCTOFKaon = (TH3F *)fInputFile->Get(Form("%s/hPID/After/hNsigma_TPC_TOF_Ka_after", multpath.c_str()));
        TH3F *hNsigmaTPCTOFPion = (TH3F *)fInputFile->Get(Form("%s/hPID/After/hNsigma_TPC_TOF_Pi_after", multpath.c_str()));
        TH3F *hNsigmaTPCKaon = (TH3F *)fInputFile->Get(Form("%s/hPID/After/hTPCnsigKa_mult_pt", multpath.c_str()));
        TH3F *hNsigmaTPCPion = (TH3F *)fInputFile->Get(Form("%s/hPID/After/hTPCnsigPi_mult_pt", multpath.c_str()));
        TH3F *hNsigmaTOFKaon = (TH3F *)fInputFile->Get(Form("%s/hPID/After/hTOFnsigKa_mult_pt", multpath.c_str()));
        TH3F *hNsigmaTOFPion = (TH3F *)fInputFile->Get(Form("%s/hPID/After/hTOFnsigPi_mult_pt", multpath.c_str()));
        if (hNsigmaTPCTOFKaon == nullptr || hNsigmaTPCTOFPion == nullptr || hNsigmaTPCKaon == nullptr || hNsigmaTPCPion == nullptr || hNsigmaTOFKaon == nullptr || hNsigmaTOFPion == nullptr)
        {
            cerr << "PID histograms after selection not found!!!!!!!!!!!!" << endl;
            return;
        }

        TH3F *hNsigmaTPCTOFKaon_before = (TH3F *)fInputFile->Get(Form("%s/hPID/Before/hNsigma_TPC_TOF_Ka_before", multpath.c_str()));
        TH3F *hNsigmaTPCTOFPion_before = (TH3F *)fInputFile->Get(Form("%s/hPID/Before/hNsigma_TPC_TOF_Pi_before", multpath.c_str()));
        TH3F *hNsigmaTPCKaon_before = (TH3F *)fInputFile->Get(Form("%s/hPID/Before/hTPCnsigKa_mult_pt", multpath.c_str()));
        TH3F *hNsigmaTPCPion_before = (TH3F *)fInputFile->Get(Form("%s/hPID/Before/hTPCnsigPi_mult_pt", multpath.c_str()));
        TH3F *hNsigmaTOFKaon_before = (TH3F *)fInputFile->Get(Form("%s/hPID/Before/hTOFnsigKa_mult_pt", multpath.c_str()));
        TH3F *hNsigmaTOFPion_before = (TH3F *)fInputFile->Get(Form("%s/hPID/Before/hTOFnsigPi_mult_pt", multpath.c_str()));
        if (hNsigmaTPCTOFKaon_before == nullptr || hNsigmaTPCTOFPion_before == nullptr || hNsigmaTPCKaon_before == nullptr || hNsigmaTPCPion_before == nullptr || hNsigmaTOFKaon_before == nullptr || hNsigmaTOFPion_before == nullptr)
        {
            cerr << "PID histograms before selection not found!!!!!!!!!!!!" << endl;
            return;
        }

        outputtype = "pdf";

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
        hMult->GetXaxis()->SetRangeUser(0, 109.5);
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
        gPad->SetLogy(1);
        SetHistoQA(hoccupancy);
        hoccupancy->GetXaxis()->SetTitle("Occupancy");
        hoccupancy->GetYaxis()->SetTitle("Counts");
        hoccupancy->GetXaxis()->SetRangeUser(0, 2200);
        hoccupancy->GetXaxis()->SetNdivisions(505);
        hoccupancy->Draw();
        cOccupancy->SaveAs(output_QA_folder + ("/Occupancy." + outputtype).c_str());

        TCanvas *cTracksData = new TCanvas("cTracksData", "Tracks Data", 1440, 720);
        SetCanvasStyle(cTracksData, 0.1, 0.03, 0.06, 0.16);
        SetHistoQA(htracksData);
        // gPad->SetLogy(1);
        gPad->SetGrid(1, 0);
        htracksData->SetTitle("Tracks Data");
        htracksData->GetYaxis()->SetTitle("Counts");
        htracksData->GetYaxis()->SetTitleOffset(0.8);
        htracksData->GetXaxis()->SetRangeUser(0, 8);
        // htracksData->SetMaximum(htracksData->GetMaximum() * 10);
        htracksData->Draw();
        cTracksData->SaveAs(output_QA_folder + ("/TracksData." + outputtype).c_str());

        TCanvas *cEventCut = new TCanvas("cEventCut", "Event Cut", 720, 720);
        SetCanvasStyle(cEventCut, 0.1, 0.05, 0.06, 0.17);
        SetHistoQA(hEventCut);
        gPad->SetGrid(1, 0);
        // gPad->SetLogy(1);
        hEventCut->GetXaxis()->SetBinLabel(4, "INEL > 0");
        hEventCut->SetBinContent(4, hEventCut->GetBinContent(9));
        hEventCut->SetTitle("Event selections (Data)");
        hEventCut->GetYaxis()->SetTitle("Counts");
        hEventCut->GetXaxis()->SetRangeUser(0, 4);
        hEventCut->SetMinimum(0);
        hEventCut->SetMaximum(hEventCut->GetBinContent(1) * 1.2);
        hEventCut->Draw();
        double firstBin = hEventCut->GetBinContent(1);

        TLatex latex;
        latex.SetTextSize(0.03);
        latex.SetTextAlign(22); // centered

        for (int i = 1; i <= 4; i++)
        {

            double content = hEventCut->GetBinContent(i);

            double percent = 0;
            if (firstBin > 0)
                percent = (content / firstBin) * 100.0;

            TString label = Form("%.1f%%", percent);

            double x = hEventCut->GetBinCenter(i);
            double y = content * 1.05; // slightly above bin

            latex.DrawLatex(x, y, label);
        }
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

            h2DNsigmaTPCTOFKaon[i] = (TH2 *)hNsigmaTPCTOFKaon->Project3D("xy");
            h2DNsigmaTPCTOFPion[i] = (TH2 *)hNsigmaTPCTOFPion->Project3D("xy");

            // ✅ Rename AFTER creation (otherwise all histograms will be same)
            h2DNsigmaTPCTOFKaon[i]->SetName(Form("xy_Ka_%d", i));
            h2DNsigmaTPCTOFPion[i]->SetName(Form("xy_Pi_%d", i));

            // Now lets make 1D projections with limit on multiplicity ranges
            int lowBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multRangesForTPCandTOF[i][0] + 0.001);
            int highBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multRangesForTPCandTOF[i][1] - 0.001);
            int lowBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i] + 0.001);
            int highBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i + 1] - 0.001);

            h1DNsigmaTPCKaon[i] = (TH1F *)hNsigmaTPCKaon->ProjectionX(Form("h1D_TPC_Ka_%d", i), lowBinMult, highBinMult, -1, -1);
            h1DNsigmaTPCPion[i] = (TH1F *)hNsigmaTPCPion->ProjectionX(Form("h1D_TPC_Pi_%d", i), lowBinMult, highBinMult, -1, -1);
            h1DNsigmaTOFKaon[i] = (TH1F *)hNsigmaTOFKaon->ProjectionX(Form("h1D_TOF_Ka_%d", i), lowBinMult, highBinMult, -1, -1);
            h1DNsigmaTOFPion[i] = (TH1F *)hNsigmaTOFPion->ProjectionX(Form("h1D_TOF_Pi_%d", i), lowBinMult, highBinMult, -1, -1);

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
            h2DNsigmaTPCTOFKaon[i]->GetYaxis()->SetTitle("n#sigma_{TOF}");
            h2DNsigmaTPCTOFKaon[i]->GetXaxis()->SetTitle("n#sigma_{TPC}");
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
            h2DNsigmaTPCTOFPion[i]->GetYaxis()->SetTitle("n#sigma_{TOF}");
            h2DNsigmaTPCTOFPion[i]->GetXaxis()->SetTitle("n#sigma_{TPC}");
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
        TLegend *leg1DPID = new TLegend(0.17, 0.82, 0.6, 0.92);
        leg1DPID->SetNColumns(3);
        leg1DPID->SetTextSize(0.028);
        leg1DPID->SetFillStyle(0);
        leg1DPID->SetBorderSize(0);
        leg1DPID->SetTextFont(42);

        for (int i = 0; i < 6; i++)
        {
            SetHistoQA(h1DNsigmaTPCKaon[i]);
            h1DNsigmaTPCKaon[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCKaon[i]->GetXaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");
            h1DNsigmaTPCKaon[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCKaon[i]->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h1DNsigmaTPCKaon[i]->SetMaximum(h1DNsigmaTPCKaon[i]->GetMaximum() * 4.7);
            h1DNsigmaTPCKaon[i]->Draw("p same");
            leg1DPID->AddEntry(h1DNsigmaTPCKaon[i], Form("%d-%d%%", multRangesForTPCandTOF[i][0], multRangesForTPCandTOF[i][1]), "p");
            // cout<<"pT range: "<<pTrangesForTPCandTOF[i]<<"-"<<pTrangesForTPCandTOF[i + 1]<<" GeV/c, mean TPC nSigma: "<<h1DNsigmaTPCKaon[i]->GetMean()<<endl;
        }
        leg1DPID->Draw();
        TLine *lineTPCKaon = new TLine(0, 0, 0, h1DNsigmaTPCKaon[0]->GetMaximum());
        lineTPCKaon->SetLineStyle(2);
        lineTPCKaon->SetLineColor(2);
        lineTPCKaon->Draw();
        cNsigmaTPCKaon_1D->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_1D." + outputtype).c_str());

        TCanvas *cNsigmaTPCPion_1D = new TCanvas("cNsigmaTPCPion_1D", "Nsigma TPC Pion 1D", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion_1D, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTPCPion_1D->SetGrid(1, 1);
        for (int i = 0; i < 6; i++)
        {
            SetHistoQA(h1DNsigmaTPCPion[i]);
            h1DNsigmaTPCPion[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCPion[i]->GetXaxis()->SetTitle("n#sigma_{TPC} #pi^{#pm}");
            h1DNsigmaTPCPion[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCPion[i]->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h1DNsigmaTPCPion[i]->SetMaximum(h1DNsigmaTPCPion[i]->GetMaximum() * 4.7);
            h1DNsigmaTPCPion[i]->Draw("p same");
        }
        leg1DPID->Draw();
        TLine *lineTPCPion = new TLine(0, 0, 0, h1DNsigmaTPCPion[0]->GetMaximum());
        lineTPCPion->SetLineStyle(2);
        lineTPCPion->SetLineColor(2);
        lineTPCPion->Draw();
        cNsigmaTPCPion_1D->SaveAs(output_QA_folder + ("/NsigmaTPCPion_1D." + outputtype).c_str());

        TCanvas *cNsigmaTOFKaon_1D = new TCanvas("cNsigmaTOFKaon_1D", "Nsigma TOF Kaon 1D", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon_1D, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFKaon_1D->SetGrid(1, 1);
        for (int i = 0; i < 6; i++)
        {
            SetHistoQA(h1DNsigmaTOFKaon[i]);
            h1DNsigmaTOFKaon[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFKaon[i]->GetXaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");
            h1DNsigmaTOFKaon[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFKaon[i]->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h1DNsigmaTOFKaon[i]->SetMaximum(h1DNsigmaTOFKaon[i]->GetMaximum() * 4.7);
            h1DNsigmaTOFKaon[i]->Draw("p same");
        }
        leg1DPID->Draw();
        TLine *lineTOFKaon = new TLine(0, 0, 0, h1DNsigmaTOFKaon[0]->GetMaximum());
        lineTOFKaon->SetLineStyle(2);
        lineTOFKaon->SetLineColor(2);
        lineTOFKaon->Draw();
        cNsigmaTOFKaon_1D->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_1D." + outputtype).c_str());

        TCanvas *cNsigmaTOFPion_1D = new TCanvas("cNsigmaTOFPion_1D", "Nsigma TOF Pion 1D", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion_1D, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFPion_1D->SetGrid(1, 1);
        for (int i = 0; i < 6; i++)
        {
            SetHistoQA(h1DNsigmaTOFPion[i]);
            h1DNsigmaTOFPion[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFPion[i]->GetXaxis()->SetTitle("n#sigma_{TOF} #pi^{#pm}");
            h1DNsigmaTOFPion[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFPion[i]->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h1DNsigmaTOFPion[i]->SetMaximum(h1DNsigmaTOFPion[i]->GetMaximum() * 4.7);
            h1DNsigmaTOFPion[i]->Draw("p same");
        }
        leg1DPID->Draw();
        TLine *lineTOFPion = new TLine(0, 0, 0, h1DNsigmaTOFPion[0]->GetMaximum());
        lineTOFPion->SetLineStyle(2);
        lineTOFPion->SetLineColor(2);
        lineTOFPion->Draw();
        cNsigmaTOFPion_1D->SaveAs(output_QA_folder + ("/NsigmaTOFPion_1D." + outputtype).c_str());

        TCanvas *cNsigmaTPCKaon_1D_pt = new TCanvas("cNsigmaTPCKaon_1D_pt", "Nsigma TPC Kaon 1D pT", 720, 720);
        SetCanvasStyle(cNsigmaTPCKaon_1D_pt, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTPCKaon_1D_pt->SetGrid(1, 1);

        TCanvas *cNsigmaTPCPion_1D_pt = new TCanvas("cNsigmaTPCPion_1D_pt", "Nsigma TPC Pion 1D pT", 720, 720);
        SetCanvasStyle(cNsigmaTPCPion_1D_pt, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTPCPion_1D_pt->SetGrid(1, 1);

        TCanvas *cNsigmaTOFKaon_1D_pt = new TCanvas("cNsigmaTOFKaon_1D_pt", "Nsigma TOF Kaon 1D pT", 720, 720);
        SetCanvasStyle(cNsigmaTOFKaon_1D_pt, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFKaon_1D_pt->SetGrid(1, 1);

        TCanvas *cNsigmaTOFPion_1D_pt = new TCanvas("cNsigmaTOFPion_1D_pt", "Nsigma TOF Pion 1D pT", 720, 720);
        SetCanvasStyle(cNsigmaTOFPion_1D_pt, 0.14, 0.05, 0.06, 0.14);
        cNsigmaTOFPion_1D_pt->SetGrid(1, 1);

        TLegend *leg1DPID_pt = new TLegend(0.17, 0.82, 0.8, 0.92);
        leg1DPID_pt->SetNColumns(3);
        leg1DPID_pt->SetTextSize(0.028);
        leg1DPID_pt->SetFillStyle(0);
        leg1DPID_pt->SetBorderSize(0);
        leg1DPID_pt->SetTextFont(42);
        leg1DPID_pt->SetHeader("#it{p}_{T} ranges");

        for (int i = 0; i < 6; i++)
        {
            cNsigmaTPCKaon_1D_pt->cd();
            SetHistoQA(h1DNsigmaTPCKaon_pt[i]);
            h1DNsigmaTPCKaon_pt[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCKaon_pt[i]->GetXaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");
            h1DNsigmaTPCKaon_pt[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCKaon_pt[i]->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h1DNsigmaTPCKaon_pt[i]->SetMaximum(h1DNsigmaTPCKaon_pt[i]->GetMaximum() * 25);
            h1DNsigmaTPCKaon_pt[i]->Draw("p same");
            leg1DPID_pt->AddEntry(h1DNsigmaTPCKaon_pt[i], Form("%.1f-%.1f", pTrangesForTPCandTOF[i], pTrangesForTPCandTOF[i + 1]), "p");

            cNsigmaTPCPion_1D_pt->cd();
            SetHistoQA(h1DNsigmaTPCPion_pt[i]);
            h1DNsigmaTPCPion_pt[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTPCPion_pt[i]->GetXaxis()->SetTitle("n#sigma_{TPC} #pi^{#pm}");
            h1DNsigmaTPCPion_pt[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTPCPion_pt[i]->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h1DNsigmaTPCPion_pt[i]->SetMaximum(h1DNsigmaTPCPion_pt[i]->GetMaximum() * 3);
            h1DNsigmaTPCPion_pt[i]->Draw("p same");

            cNsigmaTOFKaon_1D_pt->cd();
            SetHistoQA(h1DNsigmaTOFKaon_pt[i]);
            h1DNsigmaTOFKaon_pt[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFKaon_pt[i]->GetXaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");
            h1DNsigmaTOFKaon_pt[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFKaon_pt[i]->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h1DNsigmaTOFKaon_pt[i]->SetMaximum(h1DNsigmaTOFKaon_pt[i]->GetMaximum() * 70);
            h1DNsigmaTOFKaon_pt[i]->Draw("p same");

            cNsigmaTOFPion_1D_pt->cd();
            SetHistoQA(h1DNsigmaTOFPion_pt[i]);
            h1DNsigmaTOFPion_pt[i]->SetMarkerColor(colors[i]);
            h1DNsigmaTOFPion_pt[i]->GetXaxis()->SetTitle("n#sigma_{TOF} #pi^{#pm}");
            h1DNsigmaTOFPion_pt[i]->GetYaxis()->SetTitle("Counts");
            h1DNsigmaTOFPion_pt[i]->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h1DNsigmaTOFPion_pt[i]->SetMaximum(h1DNsigmaTOFPion_pt[i]->GetMaximum() * 13);
            h1DNsigmaTOFPion_pt[i]->Draw("p same");
        }

        cNsigmaTPCKaon_1D_pt->cd();
        leg1DPID_pt->Draw();
        TLine *lineTPCKaon_pt = new TLine(0, 0, 0, h1DNsigmaTPCKaon_pt[0]->GetMaximum());
        lineTPCKaon_pt->SetLineStyle(2);
        lineTPCKaon_pt->SetLineColor(2);
        lineTPCKaon_pt->Draw();
        cNsigmaTPCKaon_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_1D_pt." + outputtype).c_str());

        cNsigmaTPCPion_1D_pt->cd();
        leg1DPID_pt->Draw();
        TLine *lineTPCPion_pt = new TLine(0, 0, 0, h1DNsigmaTPCPion_pt[0]->GetMaximum());
        lineTPCPion_pt->SetLineStyle(2);
        lineTPCPion_pt->SetLineColor(2);
        lineTPCPion_pt->Draw();
        cNsigmaTPCPion_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTPCPion_1D_pt." + outputtype).c_str());

        cNsigmaTOFKaon_1D_pt->cd();
        leg1DPID_pt->Draw();
        TLine *lineTOFKaon_pt = new TLine(0, 0, 0, h1DNsigmaTOFKaon_pt[0]->GetMaximum());
        lineTOFKaon_pt->SetLineStyle(2);
        lineTOFKaon_pt->SetLineColor(2);
        lineTOFKaon_pt->Draw();
        cNsigmaTOFKaon_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_1D_pt." + outputtype).c_str());

        cNsigmaTOFPion_1D_pt->cd();
        leg1DPID_pt->Draw();
        TLine *lineTOFPion_pt = new TLine(0, 0, 0, h1DNsigmaTOFPion_pt[0]->GetMaximum());
        lineTOFPion_pt->SetLineStyle(2);
        lineTOFPion_pt->SetLineColor(2);
        lineTOFPion_pt->Draw();
        cNsigmaTOFPion_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTOFPion_1D_pt." + outputtype).c_str());
    }
}
