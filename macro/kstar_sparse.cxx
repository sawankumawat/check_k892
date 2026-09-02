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
    // TString sysVars[] = {"", "Norm1", "Norm2", "FitRange1", "FitRange2", "WidthFree"};
    TString sysVars[] = {""};
    int nSysVars = sizeof(sysVars) / sizeof(sysVars[0]);
    // const string kResBkg = "MIX";
    // const string kResBkg = "LIKE";
    string kResBkg = "ROTATED";

    string kbkg = "pol3";
    // string kbkg = "pol2";

    string outputtype = "pdf";     // pdf, eps
    const bool save_bkg_plots = 1; // save background plots
    const float txtsize = 0.045;   // text size in the plots
    bool makeallpTplots = true;    // make all pT plots
    bool calcInvMass = true;
    const bool multipanel_plots = 0;
    const bool save_plots = 1;
    bool isINEL = false;
    bool widthFixed = true; // width fixed to PDG value

    double ResolutionMCtrue[] = {0.00545008, 0.00565446, 0.0065543, 0.00658792, 0.00583034, 0.00517954, 0.00541337, 0.00556974, 0.00557882, 0.00564703, 0.00595414, 0.0061774, 0.00648009, 0.0066736, 0.00691093, 0.00727043, 0.00718425, 0.0074514, 0.00830118, 0.00842358, 0.00850154, 0.00891884, 0.0111535};
    // double ResolutionMCtrue[] = {0.009, 0.0057, 0.0065543, 0.0034, 0.0032, 0.0030, 0.0037, 0.0011, 0.00014, 0.00032, 0.0003, 0.0007, 0.0033, 0.0042, 0.0048, 0.0054, 0.0061, 0.0063, 0.0069, 0.0075, 0.0099, 0.011, 0.012};

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

    std::vector<TString> badFits;

    // for (Int_t ip = 0; ip < Npt; ip++)
    // {
    //     TString cName = TString::Format("cinv_pt_%2.1f-%2.1f", pT_bins[ip], pT_bins[ip + 1]);
    //     cinv[ip] = new TCanvas(Form("cinv%d", ip), cName.Data(), 10, 10, 720, 720);
    //     SetCanvasStyle(cinv[ip], 0.15, 0.05, 0.08, 0.13);
    // }

    // for (Int_t ip = 0; ip < Npt; ip++)
    // {
    //     TString cNam = TString::Format("cSigbkg_pt_%2.1f-%2.1f", pT_bins[ip], pT_bins[ip + 1]);
    //     cSigbkg[ip] = new TCanvas(Form("cSigbkg%d", ip), cNam.Data(), 720, 720);
    //     SetCanvasStyle(cSigbkg[ip], 0.15, 0.06, 0.06, 0.13);
    // }

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

    TFile *fTemplateFile = TFile::Open(Form("template/buildTemplate/template/%s/SignalMinusTrue%s.root", kResBkg.c_str(), kvariation.c_str()), "READ");
    if (!fTemplateFile || fTemplateFile->IsZombie())
    {
        cerr << "ERROR: SignalMinusTrue.root not found!" << endl;
        return;
    }
    cout << "Reflection template file opened successfully." << endl;

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
    // for (int ivar = 5; ivar < 6; ivar++)
    {
        if (nSysVars > 1 && (kResBkg != "ROTATED" || kbkg != "pol3"))
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

            if (kNormRangepT.size() < Npt || kRebin.size() < Npt)
            {
                cerr << "Error: kNormRangepT or kRebin arrays are not initialized for all pT bins." << endl;
                return;
            }

            if (kFitRange.size() <= imult || kFitRange[imult].size() < Npt)
            {
                cerr << "Error: kFitRange is not initialized for all multiplicity/pT bins." << endl;
                return;
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
                filecmp = new TFile((koutputfolder + resBkgFolder + kbkgFolder + Form("/yield_%d_%d.root", multlow, multhigh)).c_str(), "RECREATE");
            }
            else
            {
                filecmp = new TFile((koutputfolder + "/" + sysVars[ivar].Data() + kbkgFolder + +Form("/yield_%d_%d.root", multlow, multhigh)).c_str(), "RECREATE");
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

            std::vector<TCanvas *> c_fitsig, c_sigbkg;
            TString Cenoutputfolder = outputfolder_mult;
            TCanvas *cinv[Npt];    // for output canvases on screen containing fitted signal after subtraction
            TCanvas *cSigbkg[Npt]; // for output canvases on screen containing signal with bkg(after norm in case of mix)
            for (Int_t ip = 0; ip < Npt; ip++)
            {
                cinv[ip] = new TCanvas(Form("cinv%d", ip), "", 720, 720);
                cSigbkg[ip] = new TCanvas(Form("cSigbkg%d", ip), "", 720, 720);
            }

            if (calcInvMass)
            {
                // gStyle->SetOptStat(1110);
                gStyle->SetOptStat(0);
                gStyle->SetOptFit(0);

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

                    double lowfitrange = kFitRange[imult][ip][0];
                    double highfitrange = kFitRange[imult][ip][1];

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
                        normRangeLow = 1.15;
                        normRangeHigh = 1.20;
                    }
                    else if (sysVars[ivar] == "Norm2")
                    {
                        normRangeLow = 1.24;
                        normRangeHigh = 1.29;
                    }

                    if (kResBkg == "MIX" || kResBkg == "ROTATED")
                    {
                        TH1D *bkgclonetemp = (kResBkg == "MIX") ? (TH1D *)fHistBkg[ip]->Clone() : (TH1D *)fHistRotated1D[ip]->Clone();

                        sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(normRangeLow + 1e-5), fHistTotal[ip]->GetXaxis()->FindBin(normRangeHigh - 1e-5)));
                        bkg_integral = (bkgclonetemp->Integral(bkgclonetemp->GetXaxis()->FindBin(normRangeLow + 1e-5), bkgclonetemp->GetXaxis()->FindBin(normRangeHigh - 1e-5)));
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

                        hfbkg = (TH1D *)fHistbkgLS[ip]->Clone();
                        hfbkg->Rebin(rebin_value);
                        hfsig->Rebin(rebin_value);
                        hfsig->Add(hfbkg, -1);
                    }

                    fHistTotal[ip]->Rebin(rebin_value);

                    //**** pt binwidth************x*****************************
                    ptbinwidth[ip] = pT_bins[ip + 1] - pT_bins[ip];
                    // cout<<"the value of pt bin width is "<<ptbinwidth[ip]<<endl;

                    //******************************************************************************************//
                    //                               Fit function using template
                    //******************************************************************************************//

                    cout << "  pT bin " << ip << " (" << lowpt << " < pT < " << highpt << " GeV/c):" << endl;

                    // =====================================================================
                    // Load reflection template
                    // =====================================================================
                    TString templName = Form("%d-%d/hSigminusTrue_pt_%.1f_%.1f", multlow, multhigh, lowpt, highpt);
                    TH1D *hReflRaw = (TH1D *)fTemplateFile->Get(templName);
                    if (!hReflRaw)
                    {
                        cerr << "WARNING: Template '" << templName
                             << "' not found. Skipping pT bin " << ip << "." << endl;
                        continue;
                    }

                    TH1D *hReflection = (TH1D *)hReflRaw->Clone(Form("hReflection_ip%d", ip));

                    // if (kResBkg == "LIKE" && (imult == 1 || imult == 2 || imult == 3))
                    //     hReflection->Rebin(rebin_value * 2);
                    // else
                    hReflection->Rebin(rebin_value);

                    ////No difference is seen even if there is bin mismath error. So I have commented it out.
                    // TH1D *hReflection_truncated = new TH1D("hReflection_truncated", "Truncated Reflection Template", 300, 0.7, 1.3);
                    // for (int i = 1; i <= 300; ++i)
                    // {
                    //     hReflection_truncated->SetBinContent(i, hReflection->GetBinContent(i));
                    //     hReflection_truncated->SetBinError(i, hReflection->GetBinError(i));
                    // }
                    // hReflection_truncated->Rebin(rebin_value);

                    TCanvas *cRefl = new TCanvas(Form("cRefl_ip%d", ip), "Reflection template", 720, 720);
                    TH1D *hDatabyReflection = (TH1D *)hfsig->Clone(Form("hDatabyReflection_ip%d", ip));
                    hDatabyReflection->SetTitle(Form("%.1f < p_{T} (GeV/c) < %.1f; M_{K#pi} (GeV/c^{2}); Data / Reflection template", lowpt, highpt));

                    TH1D *hRefNorm = (TH1D *)hReflection->Clone(Form("hRefNorm_ip%d", ip));
                    // TH1D *hRefNorm = (TH1D *)hReflection_truncated->Clone(Form("hRefNorm_ip%d", ip));

                    cout << "Ref Norm = " << hDatabyReflection->Integral(hDatabyReflection->GetXaxis()->FindBin(0.7), hDatabyReflection->GetXaxis()->FindBin(0.8)) / (hRefNorm->Integral(hRefNorm->GetXaxis()->FindBin(0.7), hRefNorm->GetXaxis()->FindBin(0.8))) << endl;

                    hRefNorm->Scale(hDatabyReflection->Integral(hDatabyReflection->GetXaxis()->FindBin(0.7), hDatabyReflection->GetXaxis()->FindBin(0.8)) / (hRefNorm->Integral(hRefNorm->GetXaxis()->FindBin(0.7), hRefNorm->GetXaxis()->FindBin(0.8))));

                    // if (hDatabyReflection->GetNbinsX() != hReflection_truncated->GetNbinsX())
                    // {
                    //     cerr << "ERROR: Bin mismatch between data and reflection template for pT bin " << ip << endl;
                    //     cout << "Number of bins in data: " << hDatabyReflection->GetNbinsX() << ", Number of bins in template: " << hReflection_truncated->GetNbinsX() << endl;
                    // }

                    hDatabyReflection->Divide(hRefNorm);
                    hDatabyReflection->GetXaxis()->SetRangeUser(lowfitrange, highfitrange);
                    cRefl->cd();
                    hDatabyReflection->Draw();

                    TLatex latexref;
                    latexref.SetNDC();
                    latexref.SetTextSize(0.05);
                    latexref.SetTextAlign(22);
                    latexref.DrawLatex(0.5, 0.95, Form("%.1f < p_{T} (GeV/c) < %.1f", lowpt, highpt));

                    // auto c_clone_sigbyref = (TCanvas *)cRefl->Clone(Form("hsigbyref_pt_%d", ip + 1));
                    // c_sigbyref.push_back(c_clone_sigbyref);
                    // cRefl->SaveAs(Cenoutputfolder + Form("/hDatabyReflection_ip%d.%s", ip, outputtype.Data()));

                    for (int ib = 1; ib <= hReflection->GetNbinsX(); ib++)
                    {
                        if (hReflection->GetBinContent(ib) < 0)
                        {
                            hReflection->SetBinContent(ib, 0.0);
                            hReflection->SetBinError(ib, 0.0);
                        }
                    }

                    bool bTemplEmpty = (hReflection->Integral() <= 0);
                    if (bTemplEmpty)
                        cerr << "WARNING: Template empty in fit range for pT bin " << ip << endl;

                    // =====================================================================
                    // Superior Analytical Normalization Setup
                    // =====================================================================
                    const double fitLo = lowfitrange;
                    const double fitHi = highfitrange;
                    const double leftLo = fitLo;
                    const double leftHi = 0.82;
                    const double binWidthFit = hfsig->GetBinWidth(1);

                    int binLo = hReflection->GetXaxis()->FindBin(fitLo + 1e-6);
                    int binHi = hReflection->GetXaxis()->FindBin(fitHi - 1e-6);
                    double reflBinWidth = hReflection->GetBinWidth(1);
                    double reflNorm = std::max(hReflection->Integral(binLo, binHi) * reflBinWidth, 1e-12);

                    double total_in_fit = std::max(1.0, hfsig->Integral(
                                                            hfsig->GetXaxis()->FindBin(fitLo + 1e-5),
                                                            hfsig->GetXaxis()->FindBin(fitHi - 1e-5)));

                    // ==================================================================
                    // Define Universal Residual Background Evaluator Lambda
                    // ==================================================================
                    int nBkgPars = 1; // pol1 default
                    if (kbkg == "pol2")
                        nBkgPars = 2;
                    else if (kbkg == "pol3")
                        nBkgPars = 3;
                    else if (kbkg == "expol" || kbkg == "pol3Thresh")
                        nBkgPars = 4;

                    auto evalBkgDensity = [=](double x, const double *par) -> double
                    {
                        if (kbkg == "pol1")
                        {
                            double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
                            double chebRaw = 1.0 + par[0] * t;
                            return chebRaw / (fitHi - fitLo);
                        }
                        else if (kbkg == "pol2")
                        {
                            double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
                            double T1 = t, T2 = 2.0 * t * t - 1.0;
                            double chebRaw = 1.0 + par[0] * T1 + par[1] * T2;
                            double chebInt = ((fitHi - fitLo) / 2.0) * (2.0 - (2.0 / 3.0) * par[1]);
                            return chebRaw / std::max(chebInt, 1e-12);
                        }
                        else if (kbkg == "pol3")
                        {
                            double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
                            double T1 = t, T2 = 2.0 * t * t - 1.0, T3 = 4.0 * t * t * t - 3.0 * t;
                            double chebRaw = 1.0 + par[0] * T1 + par[1] * T2 + par[2] * T3;
                            double chebInt = ((fitHi - fitLo) / 2.0) * (2.0 - (2.0 / 3.0) * par[1]);
                            return chebRaw / std::max(chebInt, 1e-12);
                        }
                        else if (kbkg == "expol")
                        {
                            if (x <= 0)
                                return 0.0;
                            auto rawExpol = [](double val, const double *p)
                            {
                                if (val <= 0)
                                    return 0.0;
                                return std::pow(val, p[0]) * std::exp(-p[1] - p[2] * val - p[3] * val * val);
                            };
                            double rawVal = rawExpol(x, par);
                            int nSteps = 40;
                            double h = (fitHi - fitLo) / nSteps;
                            double sum = rawExpol(fitLo, par) + rawExpol(fitHi, par);
                            for (int j = 1; j < nSteps; j++)
                            {
                                double xj = fitLo + j * h;
                                sum += rawExpol(xj, par) * ((j % 2 == 0) ? 2.0 : 4.0);
                            }
                            double normInt = (h / 3.0) * sum;
                            return rawVal / std::max(normInt, 1e-12);
                        }
                        else if (kbkg == "pol3Thresh")
                        {
                            double m_thresh = massPi + massKa;
                            if (x < m_thresh)
                                return 0.0;
                            double dx = x - m_thresh;
                            double rawVal = par[0] + par[1] * dx + par[2] * dx * dx + par[3] * dx * dx * dx;

                            auto antideriv = [m_thresh](double val, const double *p)
                            {
                                if (val <= m_thresh)
                                    return 0.0;
                                double d = val - m_thresh;
                                return p[0] * d + 0.5 * p[1] * d * d + (1.0 / 3.0) * p[2] * d * d * d + 0.25 * p[3] * d * d * d * d;
                            };
                            double normInt = antideriv(fitHi, par) - antideriv(std::max(fitLo, m_thresh), par);
                            return rawVal / std::max(normInt, 1e-12);
                        }
                        double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
                        return (1.0 + par[0] * t) / (fitHi - fitLo);
                    };

                    // ==================================================================
                    // Left-sideband pre-fit
                    // ==================================================================
                    TF1 *fLeft = new TF1(Form("fLeft_ip%d", ip), [hReflection, reflNorm, fitLo, fitHi, binWidthFit](double *xx, double *pp) -> double
                                         {
                double x = xx[0];
                double Ncorr = pp[0], Nbkg = pp[1], c1 = pp[2];
                double corrDensity = hReflection->Interpolate(x) / reflNorm;
                double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
                double bkgDensity = (1.0 + c1 * t) / (fitHi - fitLo);
                return binWidthFit * (Ncorr * corrDensity + Nbkg * bkgDensity); }, leftLo, leftHi, 3);
                    fLeft->SetParameters(total_in_fit * 0.3, total_in_fit * 0.3, 0.0);
                    fLeft->SetParLimits(0, 0.0, 1e10);
                    fLeft->SetParLimits(1, 0.0, 1e10);
                    fLeft->SetParLimits(2, -1.0, 1.0);

                    hfsig->Fit(fLeft, "R Q N 0");
                    double N_temp_0 = fLeft->GetParameter(0);
                    cout << "  -> Left SB pre-fit: N_temp_0 = " << N_temp_0 << " +/- " << fLeft->GetParError(0) << endl;
                    delete fLeft;

                    // ==================================================================
                    // Full invariant-mass VOIGTIAN fit
                    // ==================================================================
                    double max_corr = total_in_fit;

                    TF1 *fTotal = new TF1(Form("fTotal_ip%d", ip), [hReflection, reflNorm, fitLo, fitHi, binWidthFit, evalBkgDensity](double *xx, double *pp) -> double
                                          {
                double x = xx[0];
                double Nsig = pp[0], mass = pp[1], width_val = pp[2], sigma = pp[3]; // Slot 3: Gaussian sigma
                double Ncorr = pp[4], Nbkg = pp[5];                   // Shifted up by 1

                double voigtRaw = TMath::Voigt(x - mass, sigma, width_val);

                // Fast Numerical Normalization dynamically valid for all pT (low and high sigma)
                int nSteps = 40; 
                double h_step = (fitHi - fitLo) / nSteps;
                double vInt = TMath::Voigt(fitLo - mass, sigma, width_val) + TMath::Voigt(fitHi - mass, sigma, width_val);
                for (int j = 1; j < nSteps; j++) {
                    double xj = fitLo + j * h_step;
                    vInt += TMath::Voigt(xj - mass, sigma, width_val) * ((j % 2 == 0) ? 2.0 : 4.0);
                }
                vInt *= h_step / 3.0;
                double sigDensity = voigtRaw / std::max(vInt, 1e-12);

                // 2. Unit-Normalized MC Correlated Template Density
                double corrDensity = 0.0;
                if (hReflection && hReflection->GetEntries() > 0) {
                    corrDensity = hReflection->Interpolate(x) / reflNorm;
                    if (corrDensity < 0) corrDensity = 0.0;
                }

                // 3. Normalized Residual Background Density (Starts at pp[6])
                double bkgDensity = evalBkgDensity(x, &pp[6]);

                return binWidthFit * (Nsig * sigDensity + Ncorr * corrDensity + Nbkg * bkgDensity); }, fitLo, fitHi, 6 + nBkgPars); // 6 base params + bkg shape params

                    fTotal->SetParNames("Nsig", "mass", "width", "sigma", "Ncorr", "Nbkg");
                    fTotal->SetParameter(0, total_in_fit * 0.40);
                    fTotal->SetParLimits(0, 0.0, total_in_fit);
                    fTotal->SetParameter(1, masspdg);
                    fTotal->SetParLimits(1, masspdg - 0.010, masspdg + 0.010);

                    if (sysVars[ivar] == "WidthFree")
                    {
                        widthFixed = false;
                    }
                    else
                    {
                        widthFixed = true;
                    }

                    // --- SLOT 2: WIDTH CONFIGURATION FOR SYSTEMATICS ---
                    fTotal->SetParameter(2, widthpdg);
                    if (widthFixed)
                        fTotal->FixParameter(2, widthpdg);
                    else
                        fTotal->SetParLimits(2, widthpdg - 0.005, widthpdg + 0.005);

                    if (ivar == 5 && ip == pt_end - 1 && imult == 10)
                        fTotal->SetParLimits(2, widthpdg - 0.010, widthpdg + 0.015);

                    // --- SLOT 3: GAUSSIAN SIGMA FLOATING ---
                    if (widthFixed)
                    {
                        fTotal->SetParameter(3, 0.002);
                        fTotal->SetParLimits(3, 0.00005, 0.030); // Allow up to 30 MeV for high-pT smearing
                    }
                    else
                    {
                        fTotal->FixParameter(3, ResolutionMCtrue[ip]); // Fix to MC resolution
                    }

                    fTotal->SetParameter(4, std::min(std::max(N_temp_0, 0.0), max_corr));
                    fTotal->SetParLimits(4, 0.0, 2 * N_temp_0);
                    fTotal->SetParameter(5, total_in_fit * 0.30);
                    fTotal->SetParLimits(5, 0.0, total_in_fit);

                    // Dynamically initialize background shape parameters starting at index 6
                    if (kbkg == "pol1")
                    {
                        fTotal->SetParName(6, "c1");
                        fTotal->SetParameter(6, -0.5);
                        fTotal->SetParLimits(6, -1.0, 1.0);
                    }
                    else if (kbkg == "pol2")
                    {
                        fTotal->SetParName(6, "c1");
                        fTotal->SetParName(7, "c2");
                        fTotal->SetParameter(6, -0.5);
                        fTotal->SetParLimits(6, -1.0, 1.0);
                        fTotal->SetParameter(7, 0.1);
                        fTotal->SetParLimits(7, -1.0, 1.0);
                    }
                    else if (kbkg == "pol3")
                    {
                        fTotal->SetParName(6, "c1");
                        fTotal->SetParName(7, "c2");
                        fTotal->SetParName(8, "c3");
                        fTotal->SetParameter(6, -0.5);
                        fTotal->SetParLimits(6, -1.0, 1.0);
                        fTotal->SetParameter(7, 0.1);
                        fTotal->SetParLimits(7, -1.0, 1.0);
                        fTotal->SetParameter(8, -0.05);
                        fTotal->SetParLimits(8, -1.0, 1.0);
                    }
                    else if (kbkg == "expol")
                    {
                        fTotal->SetParName(6, "p0");
                        fTotal->SetParName(7, "p1");
                        fTotal->SetParName(8, "p2");
                        fTotal->SetParName(9, "p3");
                        fTotal->SetParameter(6, 1.0);
                        fTotal->SetParLimits(6, -15.0, 15.0);
                        fTotal->SetParameter(7, 1.0);
                        fTotal->SetParLimits(7, -20.0, 20.0);
                        fTotal->SetParameter(8, 1.0);
                        fTotal->SetParLimits(8, -15.0, 15.0);
                        fTotal->SetParameter(9, 0.1);
                        fTotal->SetParLimits(9, -15.0, 15.0);
                    }
                    else if (kbkg == "pol3Thresh")
                    {
                        fTotal->SetParName(6, "p0");
                        fTotal->SetParName(7, "p1");
                        fTotal->SetParName(8, "p2");
                        fTotal->SetParName(9, "p3");
                        fTotal->SetParameter(6, 1.0);
                        fTotal->SetParLimits(6, 0.0, 1000.0);
                        fTotal->SetParameter(7, 1.0);
                        fTotal->SetParLimits(7, -50.0, 50.0);
                        fTotal->SetParameter(8, 0.0);
                        fTotal->SetParLimits(8, -50.0, 50.0);
                        fTotal->SetParameter(9, 0.0);
                        fTotal->SetParLimits(9, -50.0, 50.0);
                    }

                    if (bTemplEmpty)
                        fTotal->FixParameter(4, 0.0); // Shifted to Slot 4

                    // Execute fit with RQS options
                    TFitResultPtr fitResult = hfsig->Fit(fTotal, "R Q S+");
                    int fitStatus = fitResult;
                    int covStatus = fitResult->CovMatrixStatus();

                    // Superior Geometric Simplex Fallback Pipeline
                    if (fitStatus != 0 || covStatus < 2)
                    {
                        cout << "  Retrying: Simplex -> Migrad (Strategy 2)..." << endl;
                        ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Simplex");
                        hfsig->Fit(fTotal, "R Q N");
                        ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
                        ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
                        fitResult = hfsig->Fit(fTotal, "R Q S");
                        ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
                    }

                    fitStatus = fitResult;
                    covStatus = fitResult->CovMatrixStatus();
                    cout << "  Fit status: " << fitStatus << "  |  Cov status: " << covStatus << endl;
                    if (fitStatus != 0 || covStatus < 2)
                    {
                        cout << "  WARNING: Bad fit in pT bin " << ip << endl;
                        badFits.push_back(Form("Multiplicity: %d-%d %, pT Bin %d (%.1f - %.1f GeV/c)",
                                               multlow, multhigh, ip, lowpt, highpt));
                    }

                    // =====================================================================
                    // Extract parameters directly as physical integrated counts
                    // =====================================================================
                    double N_sig = fTotal->GetParameter(0);
                    double N_sig_err = fTotal->GetParError(0);
                    double mass_fit = fTotal->GetParameter(1);
                    double mass_err = fTotal->GetParError(1);
                    double width_fit = fTotal->GetParameter(2);
                    double width_err = fTotal->GetParError(2);
                    double sigma_fit = fTotal->GetParameter(3); // Gaussian Sigma
                    double sigma_err = fTotal->GetParError(3);
                    double N_temp_fit = fTotal->GetParameter(4); // Shifted to Slot 4
                    double N_res_fit = fTotal->GetParameter(5);  // Shifted to Slot 5

                    Mass[ip] = mass_fit;
                    ErrorMass[ip] = mass_err;
                    Width[ip] = width_fit;      // Now stores your fitted or varied width!
                    ErrorWidth[ip] = width_err; // Now stores your width error!
                    Yield[ip] = N_sig;

                    cout << "  Mass     = " << Mass[ip] << " +/- " << ErrorMass[ip] << " GeV/c^2" << endl;
                    cout << "  Width    = " << Width[ip] << " +/- " << ErrorWidth[ip] << " GeV/c^2" << endl;
                    cout << "  Res(Sig) = " << sigma_fit << " +/- " << sigma_err << " GeV/c^2" << endl;
                    cout << "  N_sig    = " << N_sig << " +/- " << N_sig_err << endl;
                    cout << "  N_temp   = " << N_temp_fit << endl;
                    cout << "  N_res    = " << N_res_fit << endl;

                    // =====================================================================
                    // Chi²/NDF
                    // =====================================================================
                    double chi2 = fTotal->GetChisquare();
                    int ndf = fTotal->GetNDF();
                    Chi2Ndf[ip] = (ndf > 0) ? chi2 / ndf : 0.0;
                    hChiSquare->SetBinContent(ip + 1, Chi2Ndf[ip]);
                    cout << "  chi2/NDF = " << Chi2Ndf[ip] << "  (" << ndf << " NDF)" << endl;

                    // =====================================================================
                    // Yield integration calculations (Voigtian Corrected)
                    // =====================================================================
                    auto templRawIntegral = [&](double lo, double hi) -> double
                    {
                        int bLo = hReflection->GetXaxis()->FindBin(lo + 1e-6);
                        int bHi = hReflection->GetXaxis()->FindBin(hi - 1e-6);
                        return hReflection->Integral(bLo, bHi) * reflBinWidth;
                    };

                    double m5lo = masspdg - 5 * width_fit, m5hi = masspdg + 5 * width_fit;
                    double m2lo = masspdg - 2 * width_fit, m2hi = masspdg + 2 * width_fit;

                    // Numerically integrate the exact fitted Voigtian shape
                    TF1 fTempVoigt("fTempVoigt", "[0]*TMath::Voigt(x - [1], [3], [2])", 0.0, 2.0);
                    fTempVoigt.SetParameters(1.0, mass_fit, width_fit, sigma_fit);

                    double voigtIntegralFitRange = fTempVoigt.Integral(fitLo, fitHi);
                    double f_5g = fTempVoigt.Integral(m5lo, m5hi) / std::max(voigtIntegralFitRange, 1e-12);
                    double f_2g = fTempVoigt.Integral(m2lo, m2hi) / std::max(voigtIntegralFitRange, 1e-12);

                    double N_sig_5g = N_sig * f_5g;
                    double N_sig_2g = N_sig * f_2g;

                    cout << "  Voigtian fractions: +/-5Γ=" << f_5g << "  +/-2Γ=" << f_2g << endl;

                    // =====================================================================
                    // Raw yield – integral method
                    // =====================================================================
                    yieldcalc = N_sig_5g / (Event * ptbinwidth[ip] * dy * BR);
                    yielderror = (N_sig_err * f_5g) / (Event * ptbinwidth[ip] * dy * BR);

                    hintegral_yield->SetBinContent(ip + 1, yieldcalc);
                    hintegral_yield->SetBinError(ip + 1, yielderror);
                    hmass->SetBinContent(ip + 1, Mass[ip]);
                    hmass->SetBinError(ip + 1, ErrorMass[ip]);
                    hwidth->SetBinContent(ip + 1, Width[ip]);
                    hwidth->SetBinError(ip + 1, ErrorWidth[ip]);

                    std::cout << "Yield (Functional integration): " << yieldcalc << " #pm " << yielderror << std::endl;

                    // =====================================================================
                    // Significance
                    // =====================================================================
                    bmin = hfsig->GetXaxis()->FindBin(masspdg - 2 * width_fit);
                    bmax = hfsig->GetXaxis()->FindBin(masspdg + 2 * width_fit);
                    significance_num = N_sig_2g;
                    significance_den = TMath::Sqrt(std::max(1.0, (double)fHistTotal[ip]->Integral(bmin, bmax)));
                    ratio = significance_num / significance_den;
                    hsignificance->SetBinContent(ip + 1, ratio);

                    // =====================================================================
                    // Raw yield – bin-counting method with Snapped Boundaries Fix
                    // =====================================================================
                    Yield_bincount_hist = hfsig->IntegralAndError(bmin, bmax, hBCError_1);

                    // --- SNAPPED BOUNDARIES FIX START ---
                    // IntegralAndError sums WHOLE bins, so the data integral's true range
                    // is [GetBinLowEdge(bmin), GetBinUpEdge(bmax)] -- not the exact continuous
                    // masspdg+/-2*width_fit window. We snap the continuous limits to match!
                    double m2lo_snap = hfsig->GetXaxis()->GetBinLowEdge(bmin);
                    double m2hi_snap = hfsig->GetXaxis()->GetBinUpEdge(bmax);

                    double corr_f_fit = templRawIntegral(fitLo, fitHi);
                    double corr_f_2g = templRawIntegral(m2lo_snap, m2hi_snap) / std::max(corr_f_fit, 1e-12);

                    TF1 fTempBkg("fTempBkg", [&](double *xx, double *pp)
                                 { return evalBkgDensity(xx[0], pp); }, fitLo, fitHi, nBkgPars);
                    fTempBkg.SetParameters(fTotal->GetParameters() + 6); // Shifted to Offset 6

                    double bkg_f_fit = fTempBkg.Integral(fitLo, fitHi);
                    double bkg_f_2g = fTempBkg.Integral(m2lo_snap, m2hi_snap) / std::max(bkg_f_fit, 1e-12);

                    // Calculate the snapped Voigtian fraction so tail correction is exact
                    double f_2g_snap = fTempVoigt.Integral(m2lo_snap, m2hi_snap) / std::max(voigtIntegralFitRange, 1e-12);
                    double N_sig_2g_snap = N_sig * f_2g_snap;

                    double N_temp_2g = N_temp_fit * corr_f_2g;
                    double N_res_2g = N_res_fit * bkg_f_2g;
                    fYield_BinCount = Yield_bincount_hist - N_temp_2g - N_res_2g;

                    double tail_correction = N_sig_5g - N_sig_2g_snap;
                    Total_Ybincounting = (fYield_BinCount + tail_correction) / (Event * ptbinwidth[ip] * dy * BR);
                    // --- SNAPPED BOUNDARIES FIX END ---

                    double Final_pro_error = hBCError_1 / (Event * ptbinwidth[ip] * dy * BR);

                    hYbincount->SetBinContent(ip + 1, Total_Ybincounting);
                    hYbincount->SetBinError(ip + 1, Final_pro_error);
                    hFrac_stat_error->SetBinContent(ip + 1, (Total_Ybincounting > 0) ? Final_pro_error / Total_Ybincounting : 0.0);

                    std::cout << "Yield (Bin Count): " << Total_Ybincounting << " #pm " << Final_pro_error << std::endl;
                    hsigma->SetBinContent(ip + 1, sigma_fit);
                    hsigma->SetBinError(ip + 1, sigma_err);

                    // =====================================================================
                    // PLOTTING – fitted invariant mass + Ratio Pad
                    // =====================================================================
                    if (multipanel_plots == 1)
                    {
                        if (ip < kupperpad * klowerpad)
                        {
                            cgrid1->cd(ip + 1);
                        }
                        else
                        {
                            cgrid2->cd(ip + 1 - kupperpad * klowerpad);
                        }
                    }
                    else
                    {
                        cinv[ip]->cd();
                        cinv[ip]->Clear();
                    }

                    // --- 1. Upper Pad: Main Fit ---
                    TPad *pad1 = new TPad(Form("pad1_%d", ip), Form("pad1_%d", ip), 0.0, 0.3, 1.0, 1.0);
                    pad1->SetBottomMargin(0.0);
                    pad1->SetLeftMargin(0.12);
                    pad1->SetRightMargin(0.035);
                    pad1->SetFillStyle(4000);
                    pad1->Draw();
                    pad1->cd();

                    if (hfsig->GetFunction(Form("fTotal_ip%d", ip)))
                    {
                        hfsig->GetFunction(Form("fTotal_ip%d", ip))->SetBit(TF1::kNotDraw);
                    }

                    hfsig->SetMarkerStyle(20);
                    hfsig->SetMarkerColor(kBlack);
                    hfsig->SetLineColor(kBlack);
                    hfsig->GetXaxis()->SetRangeUser(lowfitrange, highfitrange);
                    hfsig->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/c^{2}", binwidth_file * 1000));
                    hfsig->GetYaxis()->CenterTitle(1);
                    hfsig->GetYaxis()->SetMaxDigits(2);
                    hfsig->GetXaxis()->SetLabelSize(0);
                    hfsig->GetXaxis()->SetTitleSize(0);
                    hfsig->GetYaxis()->SetTitleSize(0.055);
                    hfsig->GetYaxis()->SetLabelSize(0.045);
                    hfsig->SetStats(0);

                    SetHistoQA(hfsig);
                    SetHistoQA(fHistTotal[ip]);
                    hfsig->SetMarkerSize(0.8);
                    fHistTotal[ip]->SetMarkerSize(0.8);

                    hfsig->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
                    hfsig->GetYaxis()->SetMaxDigits(3);
                    hfsig->GetYaxis()->CenterTitle(1);
                    hfsig->GetYaxis()->SetTitleOffset(1.1);
                    hfsig->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidth_file * 1000));

                    SetHistoQA(hfbkg);
                    hfbkg->SetLineColor(kRed);
                    hfbkg->SetMarkerColor(kRed);
                    hfbkg->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
                    // hfsig->GetYaxis()->SetMaxDigits(3);
                    hfbkg->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidth_file * 1000));

                    if (hfsig->GetMaximum() > 0)
                        hfsig->SetMinimum(-hfsig->GetMaximum() * 0.08);

                    hfsig->SetMaximum(hfsig->GetMaximum() * 1.3);
                    hfsig->Draw("E");

                    pad1->cd();
                    const int kSmoothNpx = 1000;
                    fTotal->SetLineColor(kRed);
                    fTotal->SetLineWidth(2); // fTotal->SetNpx(kSmoothNpx);
                    fTotal->Draw("SAME");

                    // Analytical Component Functions for Plotting
                    TF1 *fSig = new TF1(Form("fSig_ip%d", ip), [fitLo, fitHi, binWidthFit](double *xx, double *pp) -> double
                                        {
                double x = xx[0]; double Nsig_val = pp[0], mass_val = pp[1], width_val = pp[2], sigma_val = pp[3];
                double voigtRaw = TMath::Voigt(x - mass_val, sigma_val, width_val);
                int nSteps = 40; 
                double h_step = (fitHi - fitLo) / nSteps;
                double vInt = TMath::Voigt(fitLo - mass_val, sigma_val, width_val) + TMath::Voigt(fitHi - mass_val, sigma_val, width_val);
                for (int j = 1; j < nSteps; j++) {
                    double xj = fitLo + j * h_step;
                    vInt += TMath::Voigt(xj - mass_val, sigma_val, width_val) * ((j % 2 == 0) ? 2.0 : 4.0);
                }
                vInt *= h_step / 3.0;
                return binWidthFit * Nsig_val * voigtRaw / std::max(vInt, 1e-12); }, fitLo, fitHi, 4);
                    fSig->SetParameters(N_sig, mass_fit, width_fit, sigma_fit);
                    fSig->SetLineColor(kMagenta + 1);
                    fSig->SetLineStyle(kSolid);
                    fSig->SetLineWidth(2); // fSig->SetNpx(kSmoothNpx);
                    fSig->Draw("SAME");

                    TF1 *fCorr = new TF1(Form("fCorr_ip%d", ip), [hReflection, reflNorm, binWidthFit](double *xx, double *pp) -> double
                                         {
                double x = xx[0]; double Ncorr_val = pp[0];
                double corrDen = hReflection->Interpolate(x) / reflNorm;
                return binWidthFit * Ncorr_val * corrDen; }, fitLo, fitHi, 1);
                    fCorr->SetParameter(0, N_temp_fit);
                    fCorr->SetLineColor(kGreen + 2);
                    fCorr->SetLineStyle(kSolid);
                    fCorr->SetLineWidth(2); // fCorr->SetNpx(kSmoothNpx);
                    fCorr->Draw("SAME");

                    TF1 *fBkg = new TF1(Form("fBkg_ip%d", ip), [fitLo, fitHi, binWidthFit, evalBkgDensity](double *xx, double *pp) -> double
                                        {
                double x = xx[0];
                double Nbkg_val = pp[0];
                double bkgDen = evalBkgDensity(x, &pp[1]);
                return binWidthFit * Nbkg_val * bkgDen; }, fitLo, fitHi, 1 + nBkgPars);
                    fBkg->SetParameter(0, N_res_fit);
                    for (int ib = 0; ib < nBkgPars; ++ib)
                        fBkg->SetParameter(1 + ib, fTotal->GetParameter(6 + ib));
                    fBkg->SetLineColor(kBlue);
                    fBkg->SetLineStyle(kSolid);
                    fBkg->SetLineWidth(2); // fBkg->SetNpx(kSmoothNpx);
                    fBkg->Draw("SAME");

                    TF1 *fTotalBkg = new TF1(Form("fTotalBkg_ip%d", ip), [hReflection, reflNorm, fitLo, fitHi, binWidthFit, evalBkgDensity](double *xx, double *pp) -> double
                                             {
                double x = xx[0];
                double Ncorr_val = pp[0], Nbkg_val = pp[1];
                double corrDen = 0.0;
                if (hReflection && hReflection->GetEntries() > 0) {
                    corrDen = hReflection->Interpolate(x) / reflNorm;
                    if (corrDen < 0) corrDen = 0.0;
                }
                double bkgDen = evalBkgDensity(x, &pp[2]);
                return binWidthFit * (Ncorr_val * corrDen + Nbkg_val * bkgDen); }, fitLo, fitHi, 2 + nBkgPars);
                    fTotalBkg->SetParameter(0, N_temp_fit);
                    fTotalBkg->SetParameter(1, N_res_fit);
                    for (int ib = 0; ib < nBkgPars; ++ib)
                        fTotalBkg->SetParameter(2 + ib, fTotal->GetParameter(6 + ib));
                    fTotalBkg->SetLineColor(kOrange + 7);
                    fTotalBkg->SetLineStyle(kDashed);
                    fTotalBkg->SetLineWidth(2); // fTotalBkg->SetNpx(kSmoothNpx);
                    fTotalBkg->Draw("SAME");

                    auto mkLine = [](Color_t col, Style_t sty = kSolid, Width_t wid = 2) -> TLine *
                    {
                        auto *l = new TLine();
                        l->SetLineColor(col);
                        l->SetLineStyle(sty);
                        l->SetLineWidth(wid);
                        return l;
                    };

                    TLegend *legComp = new TLegend(0.13, 0.6, 0.4, 0.91);
                    legComp->SetBorderSize(0);
                    legComp->SetFillStyle(0);
                    legComp->SetTextFont(42);
                    legComp->SetTextSize(0.045);
                    legComp->AddEntry(hfsig, "Data", "ep");
                    legComp->AddEntry(mkLine(kRed), "Total Fit", "l");
                    legComp->AddEntry(mkLine(kMagenta + 1, kSolid), "Voigtian Sig", "l");
                    legComp->AddEntry(mkLine(kGreen + 2, kSolid), "MC Template", "l");

                    TString bkgLegendLabel = (kbkg == "expol") ? "Exp. Pol" : (kbkg == "pol3Thresh") ? "Pol3 Thresh"
                                                                                                     : Form("Cheby %s", kbkg.c_str());
                    // legComp->AddEntry(mkLine(kBlue, kSolid), Form("Residual Bkg (%s)", bkgLegendLabel.Data()), "l");
                    legComp->AddEntry(mkLine(kBlue, kSolid), bkgLegendLabel.Data(), "l");
                    legComp->AddEntry(mkLine(kOrange + 7, kDashed), "Total Bkg (Temp+Res)", "l");
                    legComp->Draw();

                    TLegend *legPars = new TLegend(0.49, 0.6, 0.93, 0.89);
                    legPars->SetBorderSize(0);
                    legPars->SetFillStyle(0);
                    legPars->SetTextFont(42);
                    legPars->SetTextSize(0.04);
                    legPars->AddEntry((TObject *)0, Form("Mass: %.3f #pm %.1e GeV/c^{2}", Mass[ip], ErrorMass[ip]), "");

                    if (widthFixed)
                    {
                        legPars->AddEntry((TObject *)0,
                                          Form("Width: %.1f MeV/c^{2} (fixed)", width_fit * 1000), "");
                    }
                    else
                    {
                        legPars->AddEntry((TObject *)0,
                                          Form("Width: %.1f #pm %.1f MeV/c^{2}",
                                               width_fit * 1000, width_err * 1000),
                                          "");
                    }
                    legPars->AddEntry((TObject *)0, Form("#sigma_{res}: %.4f #pm %.4f MeV/c^{2}", sigma_fit * 1000, sigma_err * 1000), "");
                    legPars->AddEntry((TObject *)0, Form("N_{sig}: %.0f #pm %.0f", N_sig, N_sig_err), "");
                    // legPars->AddEntry((TObject *)0, Form("N_{temp}: %.0f", N_temp_fit), "");
                    // legPars->AddEntry((TObject *)0, Form("N_{res}: %.0f", N_res_fit), "");
                    legPars->AddEntry((TObject *)0, Form("Fit Status %d", fitStatus), "");
                    legPars->AddEntry((TObject *)0, Form("Cov Status %d", covStatus), "");
                    legPars->AddEntry((TObject *)0, Form("#chi^{2}/NDF: %.2f", Chi2Ndf[ip]), "");
                    legPars->Draw();

                    t2->DrawLatex(0.25, 0.94, Form("#bf{%.1f < #it{p}_{T}(GeV/#it{c}) < %.1f}", pT_bins[ip], pT_bins[ip + 1]));

                    pad1->Update();
                    pad1->Modified();

                    // --- 2. Lower Pad: Data / Fit Ratio ---
                    if (multipanel_plots == 1)
                    {
                        if (ip < kupperpad * klowerpad)
                        {
                            cgrid1->cd(ip + 1);
                        }
                        else
                        {
                            cgrid2->cd(ip + 1 - kupperpad * klowerpad);
                        }
                    }
                    else
                    {
                        cinv[ip]->cd();
                    }
                    TPad *pad2 = new TPad(Form("pad2_%d", ip), Form("pad2_%d", ip), 0.0, 0.0, 1.0, 0.3);
                    pad2->SetTopMargin(0.0);
                    pad2->SetBottomMargin(0.35);
                    pad2->SetLeftMargin(0.12);
                    pad2->SetRightMargin(0.035);
                    pad2->SetFillStyle(4000);
                    pad2->Draw();
                    pad2->cd();

                    TGraphAsymmErrors *gRatio = new TGraphAsymmErrors();
                    int pt_idx = 0;
                    for (int ibin = 1; ibin <= hfsig->GetNbinsX(); ibin++)
                    {
                        double xp = hfsig->GetBinCenter(ibin);
                        if (xp < lowfitrange || xp > highfitrange)
                            continue;
                        double yp = hfsig->GetBinContent(ibin);
                        double yfit = fTotal->Eval(xp);
                        if (yfit > 0)
                        {
                            gRatio->SetPoint(pt_idx, xp, yp / yfit);
                            gRatio->SetPointError(pt_idx, 0, 0,
                                                  hfsig->GetBinErrorLow(ibin) / yfit,
                                                  hfsig->GetBinErrorUp(ibin) / yfit);
                            pt_idx++;
                        }
                    }

                    TH1D *hRatioFrame = new TH1D(Form("hRatioFrame_%d", ip), "", 100, lowfitrange, highfitrange);
                    hRatioFrame->GetXaxis()->SetTitle("M_{K#pi} (GeV/c^{2})");
                    hRatioFrame->GetYaxis()->SetTitle("Data / Fit");
                    hRatioFrame->SetStats(0);

                    hRatioFrame->GetXaxis()->SetTitleSize(0.13);
                    hRatioFrame->GetXaxis()->SetLabelSize(0.11);
                    hRatioFrame->GetXaxis()->SetTitleOffset(1.05);
                    hRatioFrame->GetYaxis()->SetTitleSize(0.12);
                    hRatioFrame->GetYaxis()->SetLabelSize(0.10);
                    hRatioFrame->GetYaxis()->SetTitleOffset(0.4);
                    hRatioFrame->GetYaxis()->SetNdivisions(505);

                    hRatioFrame->SetMinimum(0.85);
                    hRatioFrame->SetMaximum(1.15);
                    hRatioFrame->Draw("AXIS");

                    gRatio->SetMarkerStyle(20);
                    gRatio->SetMarkerSize(0.5);
                    gRatio->SetLineColor(kBlack);
                    gRatio->SetMarkerColor(kBlack);
                    gRatio->Draw("P SAME");

                    TLine *lineRatio = new TLine(lowfitrange, 1.0, highfitrange, 1.0);
                    lineRatio->SetLineColor(kRed);
                    lineRatio->SetLineStyle(2);
                    lineRatio->Draw("SAME");

                    // --- Finalize and Save ---
                    if (multipanel_plots == 1)
                    {
                        (ip < kupperpad * klowerpad) ? cgrid1->cd(ip + 1) : cgrid2->cd(ip + 1 - kupperpad * klowerpad);
                        cgrid1->Modified();
                        cgrid1->Update();
                        if (Npt > 9)
                        {
                            cgrid2->Modified();
                            cgrid2->Update();
                        }
                    }
                    else
                    {
                        cinv[ip]->cd();
                        cinv[ip]->Modified();
                        cinv[ip]->Update();
                        auto c_clone_fit = (TCanvas *)cinv[ip]->Clone(Form("hfitsig_pt_%d", ip + 1));
                        c_fitsig.push_back(c_clone_fit);
                        cinv[ip]->SaveAs(Form((Cenoutputfolder + "/hfitsig_pt%d." + outputtype).Data(), ip + 1));
                    }

                    // =====================================================================
                    // PLOTTING – signal + combinatorial background
                    // =====================================================================
                    if (multipanel_plots == 1)
                    {
                        (ip < klowerpad * kupperpad) ? cgrid_bkg1->cd(ip + 1) : cgrid_bkg2->cd(ip + 1 - klowerpad * kupperpad);
                    }
                    else
                    {
                        cSigbkg[ip]->cd();
                        cSigbkg[ip]->Clear();
                    }
                    TH1F *hbkg_nopeak = (TH1F *)hfbkg->Clone();
                    hbkg_nopeak->SetLineColor(kRed);
                    hbkg_nopeak->SetMarkerColor(kRed);
                    hbkg_nopeak->SetFillColor(kRed);
                    hbkg_nopeak->SetFillStyle(3001);
                    for (int ib = 0; ib < hbkg_nopeak->GetNbinsX(); ib++)
                    {
                        double bc = hbkg_nopeak->GetBinCenter(ib + 1);
                        if (bc < kNormRangepT[ip][0] || bc > kNormRangepT[ip][1])
                            hbkg_nopeak->SetBinContent(ib + 1, -999);
                    }
                    gPad->SetRightMargin(0.05);
                    gPad->SetLeftMargin(0.15);
                    gPad->SetTopMargin(0.09);
                    gPad->SetBottomMargin(0.12);

                    SetHistoStyle(fHistTotal[ip], 1, 8, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);
                    SetHistoStyle(hfbkg, kRed, 24, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

                    fHistTotal[ip]->SetMaximum(fHistTotal[ip]->GetMaximum() * 1.15);
                    fHistTotal[ip]->SetMarkerSize(0.5);
                    hfbkg->SetMarkerSize(0.5);
                    fHistTotal[ip]->Draw("E");
                    fHistTotal[ip]->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/c^{2}", binwidth_file * 1000));
                    fHistTotal[ip]->GetXaxis()->SetTitle("M_{K#pi} (GeV/c^{2})");
                    fHistTotal[ip]->GetXaxis()->SetLabelOffset(0.015);
                    fHistTotal[ip]->GetYaxis()->SetMaxDigits(3);
                    fHistTotal[ip]->SetStats(0);
                    hfbkg->Draw("E same");

                    TLegend *leg112 = new TLegend(0.18, 0.80, 0.50, 0.893, NULL, "brNDC");
                    leg112->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
                    SetLegendStyle(leg112);
                    leg112->SetTextSize(0.035);
                    if (kResBkg == "MIX")
                        leg112->AddEntry(hfbkg, "Mixed-event bkg", "p");
                    else if (kResBkg == "LIKE")
                        leg112->AddEntry(hfbkg, "Like-sign pairs", "p");
                    else if (kResBkg == "ROTATED")
                        leg112->AddEntry(hfbkg, "Rotated unlike-sign pairs", "p");
                    if (kResBkg == "MIX" || kResBkg == "ROTATED")
                    {
                        hbkg_nopeak->Draw("BAR same");
                        leg112->AddEntry(hbkg_nopeak, "Normalisation range", "f");
                    }
                    leg112->SetNColumns(2);
                    leg112->SetColumnSeparation(0.3);
                    leg112->Draw();

                    TLatex *ltx = new TLatex(0.27, 0.95,
                                             Form("%0.1f < #it{p}_{T}(GeV/#it{c}) < %0.1f", pT_bins[ip], pT_bins[ip + 1]));
                    ltx->SetNDC();
                    ltx->SetTextFont(22);
                    ltx->SetTextSize(0.06);
                    ltx->Draw();

                    if (multipanel_plots == 1)
                    {
                        (ip < klowerpad * kupperpad) ? cgrid_bkg1->cd(ip + 1) : cgrid_bkg2->cd(ip + 1 - klowerpad * kupperpad);
                        cgrid_bkg1->Modified();
                        cgrid_bkg1->Update();
                        if (Npt > 9)
                        {
                            cgrid_bkg2->Modified();
                            cgrid_bkg2->Update();
                        }
                    }
                    else
                    {
                        cSigbkg[ip]->Modified();
                        cSigbkg[ip]->Update();
                        gPad->Modified();
                        gPad->Update();
                        auto c_clone_sig = (TCanvas *)cSigbkg[ip]->Clone(Form("hsigbkg_pt_%d", ip + 1));
                        c_sigbkg.push_back(c_clone_sig);
                        cSigbkg[ip]->SaveAs(Form((Cenoutputfolder + "/hsigbkg_pt%d." + outputtype).Data(), ip + 1));
                        // cSigbkg[ip]->Close();
                    }

                    // delete fSig;
                    // delete fCorr;
                    // delete fBkg;
                    // delete fTotalBkg;
                    // delete fTotal;

                    // ////////////////////////////////////////////////////////////////////////
                } // pt loop ends

                //==============================================
                //              END OF PT LOOP
                //==============================================

                if (multipanel_plots == 1 && save_plots == 1)
                {
                    cgrid1->Modified();
                    cgrid1->Update();
                    cgrid_bkg1->Modified();
                    cgrid_bkg1->Update();
                    cgrid1->SaveAs(outputfolder_mult + (Form("/grid1_mult%d_%d.", multlow, multhigh) + outputtype).c_str());
                    cgrid_bkg1->SaveAs(outputfolder_mult + (Form("/gridBkg1_mult%d_%d.", multlow, multhigh) + outputtype).c_str());

                    if (Npt >= klowerpad * kupperpad)
                    {
                        cgrid2->Modified();
                        cgrid2->Update();
                        cgrid_bkg2->Modified();
                        cgrid_bkg2->Update();
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

                    ////// Sigma vs pT
                    hsigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                    hsigma->GetYaxis()->SetTitle("Resolution #sigma_{res} (GeV/c^{2})");
                    SetHistoQA(hsigma);
                    hsigma->SetMaximum(hsigma->GetMaximum() * 1.5);
                    hsigma->SetMinimum(0);
                    hsigma->Draw("pe");
                    hsigma->Write("sigma");
                    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
                    csig->SaveAs((Cenoutputfolder + "/sigma_pt." + outputtype).Data());
                    csig->Clear();

                    // // // Yield vs pT (integral method)
                    SetHistoQA(hintegral_yield);
                    hintegral_yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                    hintegral_yield->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
                    gPad->SetLogy(1);
                    hintegral_yield->GetYaxis()->SetTitleOffset(1.5);
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
                    TLatex latYield;
                    latYield.SetTextSize(0.05);
                    latYield.SetTextFont(42);
                    latYield.DrawLatexNDC(0.65, 0.83, "pp, INEL");
                    latYield.DrawLatexNDC(0.65, 0.71, "K*(892)^{0}");
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

                    // =====================================================================
                    // Yield Comparison with Ratio Plot
                    // =====================================================================
                    csig->Clear();

                    // --- 1. Upper Pad: Yield Spectra ---
                    TPad *padYield1 = new TPad("padYield1", "padYield1", 0.0, 0.3, 1.0, 1.0);
                    padYield1->SetBottomMargin(0.0);
                    padYield1->SetLeftMargin(0.18);
                    padYield1->SetLogy(1);
                    padYield1->Draw();
                    padYield1->cd();

                    hintegral_yield->SetMarkerColor(kRed);
                    hintegral_yield->SetLineColor(kRed);
                    hYbincount->SetMarkerColor(kBlue);
                    hYbincount->SetLineColor(kBlue);
                    hintegral_yield->SetStats(0);
                    hYbincount->SetStats(0);

                    hintegral_yield->GetXaxis()->SetLabelSize(0);
                    hintegral_yield->GetXaxis()->SetTitleSize(0);
                    hintegral_yield->GetYaxis()->SetTitleSize(0.05);
                    hintegral_yield->GetYaxis()->SetLabelSize(0.045);
                    hintegral_yield->GetYaxis()->SetTitleOffset(1.5);

                    hintegral_yield->Draw("pe");
                    hYbincount->Draw("pe same");

                    TLegend *legcomp = new TLegend(0.45, 0.75, 0.85, 0.88);
                    SetLegendStyle(legcomp);
                    legcomp->SetTextSize(0.045);
                    legcomp->AddEntry(hintegral_yield, "Integral yield", "ep");
                    legcomp->AddEntry(hYbincount, "Bin-count yield", "ep");
                    legcomp->Draw();

                    // --- 2. Lower Pad: Ratio (Bin-count / Integral) ---
                    csig->cd();
                    TPad *padYield2 = new TPad("padYield2", "padYield2", 0.0, 0.0, 1.0, 0.3);
                    padYield2->SetTopMargin(0.0);
                    padYield2->SetBottomMargin(0.35);
                    padYield2->SetLeftMargin(0.18);
                    padYield2->Draw();
                    padYield2->cd();

                    TH1D *hYieldRatio = (TH1D *)hYbincount->Clone("hYieldRatio");
                    hYieldRatio->Divide(hintegral_yield);
                    hYieldRatio->SetTitle("");
                    hYieldRatio->GetYaxis()->SetTitle("Bin-count / Integral");
                    hYieldRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");

                    hYieldRatio->GetXaxis()->SetTitleSize(0.13);
                    hYieldRatio->GetXaxis()->SetLabelSize(0.11);
                    hYieldRatio->GetXaxis()->SetTitleOffset(1.1);
                    hYieldRatio->GetYaxis()->SetTitleSize(0.09);
                    hYieldRatio->GetYaxis()->SetLabelSize(0.08);
                    hYieldRatio->GetYaxis()->SetTitleOffset(0.8);
                    hYieldRatio->GetYaxis()->SetNdivisions(505);

                    hYieldRatio->SetMinimum(0.85);
                    hYieldRatio->SetMaximum(1.15);
                    hYieldRatio->SetMarkerColor(kBlack);
                    hYieldRatio->SetLineColor(kBlack);
                    hYieldRatio->Draw("pe");

                    TLine *lineYR = new TLine(hYieldRatio->GetXaxis()->GetXmin(), 1.0,
                                              hYieldRatio->GetXaxis()->GetXmax(), 1.0);
                    lineYR->SetLineColor(kRed);
                    lineYR->SetLineStyle(2);
                    lineYR->Draw("same");

                    csig->SaveAs(outputfolder_mult + "/yield_compare.png");
                    csig->Close();

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

            // Delete at the end of the imult loop (outside the pT loop):
            for (Int_t ip = 0; ip < Npt; ip++)
            {
                delete cinv[ip];
                delete cSigbkg[ip];
            }
        } // End of multiplicity loop

        cout << "============= End of the code =============" << endl;
        cout << "Data file used: " << kDataFilename.c_str() << endl;
        cout << "Selection used: " << (isINEL ? "INEL" : "INEL > 0") << endl;
        cout << "Format of output plots: " << outputtype.c_str() << endl;
        cout << "Residual background function: " << kbkg.c_str() << endl;
        cout << "Number of pT bins: " << Npt << endl;
        (multipanel_plots) ? cout << "Canvas: " << klowerpad << "x" << kupperpad << " multi-panel" << endl : cout << "Canvas: Single panel plots" << endl;
    }

    // =====================================================================
    // BAD FITS SUMMARY
    // =====================================================================
    cout << "\n==========================================" << endl;
    cout << "           BAD FITS SUMMARY               " << endl;
    cout << "==========================================" << endl;
    if (badFits.empty())
    {
        cout << "All fits converged successfully!" << endl;
    }
    else
    {
        cout << "Total bad fits detected: " << badFits.size() << endl;
        for (const auto &bad : badFits)
        {
            cout << "  - " << bad << endl;
        }
    }
    cout << "==========================================\n"
         << endl;

    // Stop the stopwatch
    timer.Stop();

    // Get the elapsed time
    Double_t realTime = timer.RealTime(); // Wall clock time in seconds
    Double_t cpuTime = timer.CpuTime();   // CPU time used in seconds

    // Print the elapsed times
    std::cout << "Real time elapsed: " << realTime << " seconds" << std::endl;
    std::cout << "CPU time used: " << cpuTime << " seconds" << std::endl;
}
