
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
    const string outputtype = "png"; // pdf, eps
    const bool save_bkg_plots = 1;   // save background plots
    const float txtsize = 0.045;     // text size in the plots
    //********************************************************************************

    //*************************Create folders********************************************
    TString outputfolder = kSignalOutput + "/" + kfoldername;
    TString output_root_folder = kSignalOutput + "/" + kfoldername + "/rootfiles";
    // Create the folder using TSystem::mkdir()
    if (gSystem->mkdir(outputfolder, kTRUE))
    {
        std::cout << "Folder " << outputfolder << " created successfully." << std::endl;
    }
    if (gSystem->mkdir(output_root_folder, kTRUE))
    {
        std::cout << "Folder " << output_root_folder << " created successfully." << std::endl;
    }
    //***********************************************************************************

    TCanvas *cgrid1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid_bkg1 = new TCanvas("", "", kcanvaswidth, kcanvasheight);
    TCanvas *cgrid_bkg2 = new TCanvas("", "", kcanvaswidth, kcanvasheight);

    // some initializations ********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************

    double lowfitrange[Npt + 20];
    double highfitrange[Npt + 20];
    for (int i = 0; i < Npt; i++)
    {
        lowfitrange[i] = kFitRange[i][0];
        highfitrange[i] = kFitRange[i][1];
    }

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
        SetCanvasStyle(cSigbkg[ip], 0.15, 0.05, 0.08, 0.13);
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

    TH1F *hmult = (TH1F *)fInputFile->Get("event-selection-task/hColCounterAcc");
    if (hmult == nullptr)
    {
        cerr << "Histogram not found" << endl;
        return;
    }
    double Event = hmult->GetEntries();
    cout << "*****************number of events********************:" << Event << endl;

    //**Invariant mass histograms for sig+bkg and mixed event bg***********************************************************************

    THnSparseF *fHistNum = (THnSparseF *)fInputFile->Get(Form("%s/h3k892invmassDS", kfoldername.c_str()));
    THnSparseF *fHistDen = (THnSparseF *)fInputFile->Get(Form("%s/h3k892invmassME", kfoldername.c_str()));
    THnSparseF *fHistLS = (THnSparseF *)fInputFile->Get(Form("%s/h3k892invmassLS", kfoldername.c_str()));
    THnSparseF *fHistRotated = (THnSparseF *)fInputFile->Get(Form("%s/h3K892InvMassRotation", kfoldername.c_str()));

    if (fHistNum == nullptr || fHistDen == nullptr || fHistLS == nullptr || fHistRotated == nullptr)
    {
        cerr << "Invariant mass histograms not found!!!!!!!!!!!!" << endl;
        return;
    }

    // cout << " THE NUMBER OF BINS IN THE HISTOGRAM IS " << fHistNum->GetNbinsZ()<<endl;

    gstyle(); // this is not gStyle, it is defined in the header file style.h
    gStyle->SetOptStat(1110);
    // gStyle->SetOptFit(0);

    for (Int_t ip = pt_start; ip < pt_end; ip++) // start pt bin loop
    {
        lowpt = pT_bins[ip];
        highpt = pT_bins[ip + 1];
        int multlow = 0;
        int multhigh = 110;
        int lbin = fHistNum->GetAxis(1)->FindBin(lowpt + 1e-5);
        int hbin = fHistNum->GetAxis(1)->FindBin(highpt - 1e-5);

        fHistNum->GetAxis(1)->SetRange(lbin, hbin);
        fHistDen->GetAxis(1)->SetRange(lbin, hbin);
        fHistLS->GetAxis(1)->SetRange(lbin, hbin);
        fHistRotated->GetAxis(1)->SetRange(lbin, hbin);

        int lbinmult = fHistNum->GetAxis(0)->FindBin(multlow + 1e-5);
        int hbinmult = fHistNum->GetAxis(0)->FindBin(multhigh - 1e-5);
        fHistNum->GetAxis(0)->SetRange(multlow, multhigh);
        fHistDen->GetAxis(0)->SetRange(multlow, multhigh);
        fHistLS->GetAxis(0)->SetRange(multlow, multhigh);
        fHistRotated->GetAxis(0)->SetRange(multlow, multhigh);

        fHistTotal[ip] = fHistNum->Projection(2, "E");
        fHistBkg[ip] = fHistDen->Projection(2, "E");
        fHistbkgLS[ip] = fHistLS->Projection(2, "E");
        fHistRotated1D[ip] = fHistRotated->Projection(2, "E");
        fHistNum->SetName(Form("fHistNum_%d", ip));
        fHistDen->SetName(Form("fHistDen_%d", ip));
        fHistLS->SetName(Form("fHistLS_%d", ip));
        fHistRotated->SetName(Form("fHistRotated_%d", ip));

        auto energylow = fHistTotal[ip]->GetXaxis()->GetXmin();
        auto energyhigh = fHistTotal[ip]->GetXaxis()->GetXmax();

        // cout<<"energy low value is "<<energylow<<endl;
        // cout<<"energy high value is "<<energyhigh<<endl;

        //**Cloning sig+bkg histogram for like sign or mixed event subtraction *********************************************************
        TH1D *hfsig = (TH1D *)fHistTotal[ip]->Clone();
        auto binwidth_file = (fHistTotal[ip]->GetXaxis()->GetXmax() - fHistTotal[ip]->GetXaxis()->GetXmin()) * kRebin[ip] / fHistTotal[ip]->GetXaxis()->GetNbins();
        cout << "The value of binwidth_file is: " << binwidth_file << endl;
        //*****************************************************************************************************************************

        if (kResBkg == "MIX")
        {
            TH1D *bkgclonetemp = (TH1D *)fHistBkg[ip]->Clone();

            sigbkg_integral = (fHistTotal[ip]->Integral(fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]), fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1])));
            bkg_integral = (bkgclonetemp->Integral(bkgclonetemp->GetXaxis()->FindBin(kNormRangepT[ip][0]), bkgclonetemp->GetXaxis()->FindBin(kNormRangepT[ip][1])));
            normfactor = sigbkg_integral / bkg_integral; // scaling factor for mixed bkg
            cout << "\n\n normalization factor " << 1 / normfactor << "\n\n";
            hfbkg = (TH1D *)bkgclonetemp->Clone();
            hfbkg->Scale(normfactor);
            hfbkg->Rebin(kRebin[ip]);
            hfsig->Rebin(kRebin[ip]);
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
            hfbkg->Rebin(kRebin[ip]);
            hfsig->Rebin(kRebin[ip]);
            hfsig->Add(hfbkg, -1);
        }
        else if (kResBkg == "ROTATED")
        {
            hfbkg = (TH1D *)fHistRotated1D[ip]->Clone();
            hfbkg->Scale(0.5);
            hfbkg->Rebin(kRebin[ip]);
            hfsig->Rebin(kRebin[ip]);
            hfsig->Add(hfbkg, -1);
        }

        fHistTotal[ip]->Rebin(kRebin[ip]);

        //**** pt binwidth************x*****************************
        ptbinwidth[ip] = pT_bins[ip + 1] - pT_bins[ip];
        // cout<<"the value of pt bin width is "<<ptbinwidth[ip]<<endl;

        //****************************************************************************************************

        TF1 *fitFcn, *fitFcn1;

        if (kbkg == "pol2")
        {
            fitFcn = new TF1("fitfunc", BreitWignerpoly2, kFitRange[ip][0], kFitRange[ip][1], 6);
            fitFcn1 = new TF1("fitfunc1", polynomial2, kFitRange[ip][0], kFitRange[ip][1], 3);
        }
        else if (kbkg == "pol3")
        {
            fitFcn = new TF1("fitfunc", BreitWignerpoly3, kFitRange[ip][0], kFitRange[ip][1], 7);
            fitFcn1 = new TF1("fitfunc1", polynomial3, kFitRange[ip][0], kFitRange[ip][1], 4);
        }
        else if (kbkg == "expol")
        {
            fitFcn = new TF1("fitfunc", BWExpo, kFitRange[ip][0], kFitRange[ip][1], 7);
            fitFcn1 = new TF1("fitfunc1", Expo, kFitRange[ip][0], kFitRange[ip][1], 4);
        }

        TF1 *fitFcn2 = new TF1("fitFcn2", BW, kFitRange[ip][0], kFitRange[ip][1], 3); // only signal

        fitFcn->SetParLimits(0, 0.80, 0.98); // Mass
        fitFcn->SetParameter(0, 0.895);
        fitFcn->SetParLimits(2, 0, 10e9); // Yield
        fitFcn->FixParameter(1, 0.047);   // width

        fitFcn->SetParNames("Mass", "Width", "Yield", "A", "B", "C", "D");
        // Redirect standard output to /dev/null
        // int old_stdout = dup(1);
        // freopen("/dev/null", "w", stdout);

        r = hfsig->Fit(fitFcn, "REBMSQ+"); // signal after bkg subtraction

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

        TF1 *fitFcn2_plusm = new TF1("fitFcn2_plusm", BW, lowfitrange[ip], highfitrange[ip], 3);
        TF1 *fitFcn2_minusm = new TF1("fitFcn2_minusm", BW, lowfitrange[ip], highfitrange[ip], 3);
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

        SetHistoStyle(hfsig, 1, 20, 1, 0.05, 0.045, 0.045, 0.045, 1.13, 1.8);
        SetHistoStyle(fHistTotal[ip], 1, 8, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

        hfsig->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
        hfsig->GetYaxis()->SetMaxDigits(2);
        hfsig->GetYaxis()->CenterTitle(1);
        hfsig->GetYaxis()->SetTitle(Form("Counts/%.3f", binwidth_file));

        SetHistoStyle(hfbkg, kRed, 24, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);
        hfbkg->GetXaxis()->SetTitle("M_{K#pi} (Gev/#it{c}^{2})");
        // hfsig->GetYaxis()->SetMaxDigits(2);
        hfbkg->GetYaxis()->SetTitle(Form("Counts/%.3f", binwidth_file));

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
        TLegend *pag = new TLegend(0.4, 0.7, 0.91, 0.9);
        TLegend *pag2 = new TLegend(0.2, 0.7, 0.45, 0.9);
        pag->SetBorderSize(0);
        pag->SetTextFont(42);
        pag->SetTextSize(0.04);
        pag->SetFillStyle(0);
        pag2->SetBorderSize(0);
        pag2->SetTextFont(42);
        pag2->SetTextSize(0.04);
        pag2->SetFillStyle(0);
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
        hfsig->GetXaxis()->SetRangeUser(0.64, 1.35);
        hfsig->SetMarkerSize(0.5);
        hfsig->Draw("e");
        fitFcn->Draw("same");
        fitFcn1->Draw("same");
        fitFcn2->Draw("same");
        pag2->AddEntry(fitFcn, "BW+pol3");
        pag2->AddEntry(fitFcn1, "BW");
        pag2->AddEntry(fitFcn2, "pol3");
        double fitprob = fitFcn->GetProb();
        pag->AddEntry((TObject *)0, Form("Mass: %.3f #pm %.3f", Mass[ip], ErrorMass[ip]), "");
        pag->AddEntry((TObject *)0, Form("Width: %.3f #pm %.3f", Width[ip], ErrorWidth[ip]), "");
        pag->AddEntry((TObject *)0, Form("Yield: %.1e #pm %.1e", yieldcalc * Event, yielderror * Event), "");
        // pag->AddEntry((TObject *)0, Form("Probability: %f ", fitprob), "");
        pag->AddEntry((TObject *)0, Form("#chi^{2}/NDF: %.2f ", chibyndf), "");
        // pag->Draw();
        pag2->Draw();
        // t2->DrawLatex(0.26, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
        t2->DrawLatex(0.27, 0.95,
                      Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", pT_bins[ip],
                           pT_bins[ip + 1]));
        if (multipanel_plots == 0 && save_plots == 1)
            cinv[ip]->SaveAs(Form((koutputfolder + "/hfitsig_pt%d." + outputtype).c_str(), ip + 1));
        if (multipanel_plots == 1)
            cinv[ip]->Close();

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
        gPad->SetRightMargin(0.01);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        fHistTotal[ip]->SetMaximum(fHistTotal[ip]->GetMaximum() * 1.15);
        fHistTotal[ip]->SetMarkerSize(0.5);
        hfbkg->SetMarkerSize(0.5);
        // fHistTotal[ip]->GetXaxis()->SetRangeUser(0.6, 3.0);
        fHistTotal[ip]->Draw("E");
        fHistTotal[ip]->GetYaxis()->SetTitle(Form("Counts/%.3f", binwidth_file));
        TLegend *leg112 = new TLegend(0.60554, 0.7812735, 0.852902, 0.8938954, NULL, "brNDC");
        leg112->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
        SetLegendStyle(leg112);
        leg112->SetTextSize(txtsize);
        sprintf(name, "%0.2f<p_{T}(GeV/c)<%0.2f", pT_bins[ip], pT_bins[ip + 1]);
        TLatex *ltx = new TLatex(0.27, 0.95, name);
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.06);
        hfbkg->Draw("E same");
        fHistTotal[ip]->SetMarkerSize(0.3);
        hfbkg->SetMarkerSize(0.3);
        if (kResBkg == "MIX")
            hbkg_nopeak->Draw("BAR same");
        (kResBkg == "MIX") ? leg112->AddEntry(hfbkg, "Mixed-event bkg", "p") : leg112->AddEntry(hfbkg, "Like sign pairs", "p");
        // leg112->Draw();
        ltx->Draw();
        if (multipanel_plots == 0 && save_plots == 1 && save_bkg_plots == 1)
            cSigbkg[ip]->SaveAs(Form((koutputfolder + "/hsigbkg_pt%d." + outputtype).c_str(), ip + 1));
        cSigbkg[ip]->Close();

        // ////////////////////////////////////////////////////////////////////////

    } // pt loop ends
    if (multipanel_plots == 1 && save_plots == 1)
    {
        cgrid1->SaveAs((koutputfolder + "/grid1." + outputtype).c_str());
        cgrid2->SaveAs((koutputfolder + "/grid2." + outputtype).c_str());
        cgrid_bkg1->SaveAs((koutputfolder + "/grid_bkg1." + outputtype).c_str());
        cgrid_bkg2->SaveAs((koutputfolder + "/grid_bkg2." + outputtype).c_str());
    }
    if (multipanel_plots == 0)
    {
        cgrid1->Close();
        cgrid2->Close();
        cgrid_bkg1->Close();
        cgrid_bkg2->Close();
    }

    // TFile *filecmp = new TFile((koutputfolder + "/" + kDataset + ".root").c_str(), "RECREATE");

    TCanvas *csig = new TCanvas("", "", 720, 720);
    SetCanvasStyle(csig, 0.18, 0.05, 0.08, 0.15);
    SetHistoQA(hsignificance);
    hsignificance->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hsignificance->Draw();
    // hsignificance->Write("significance");
    csig->SaveAs((koutputfolder + "/significance." + outputtype).c_str());
    csig->Clear();

    // // // chisquare_NDF vs pt

    hChiSquare->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hChiSquare->GetYaxis()->SetTitle("#chi^{2}/NDF ");
    SetHistoQA(hChiSquare);
    hChiSquare->Draw("Pe");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    csig->SaveAs((koutputfolder + "/chi." + outputtype).c_str());
    // hChiSquare->Write("chils");
    csig->Clear();

    // // mass vs pt
    hmass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hmass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    SetHistoQA(hmass);
    hmass->SetMaximum(0.91);
    hmass->Draw("pe");
    // hmass->Write("mass");
    TLegend *massleg = new TLegend(0.65, 0.2, 0.9, 0.3);
    SetLegendStyle(massleg);
    massleg->SetTextSize(txtsize);
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    TLine *line = new TLine(hmass->GetXaxis()->GetXmin(), 0.895, hmass->GetXaxis()->GetXmax(), 0.895);
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->SetLineWidth(3);
    line->Draw();
    massleg->AddEntry(line, "PDG Mass", "l");
    massleg->Draw("l");
    csig->SaveAs((koutputfolder + "/mass." + outputtype).c_str());
    csig->Clear();

    // // // Width vs pT
    hwidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hwidth->GetYaxis()->SetTitle("Width (GeV)");
    SetHistoQA(hwidth);
    hwidth->SetMaximum(hwidth->GetMaximum() * 2);
    hwidth->SetMinimum(0);
    hwidth->Draw("pe");
    // hwidth->Write("width");
    TLegend *widthleg = new TLegend(0.2, 0.75, 0.4, 0.85);
    SetLegendStyle(widthleg);
    widthleg->SetTextSize(txtsize);
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    TLine *line2 = new TLine(hwidth->GetXaxis()->GetXmin(), 0.047, hwidth->GetXaxis()->GetXmax(), 0.047);
    line2->SetLineStyle(2);
    line2->SetLineColor(2);
    line2->SetLineWidth(3);
    line2->Draw();
    widthleg->AddEntry(line, "PDG Width", "l");
    widthleg->SetFillStyle(0);
    widthleg->Draw();
    csig->SaveAs((koutputfolder + "/width_pt." + outputtype).c_str());
    csig->Clear();

    // // // Yield vs pT (integral method)
    TFile *fyield = new TFile((koutputfolder + "/yield.root").c_str(), "RECREATE");
    SetHistoQA(hintegral_yield);
    hintegral_yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hintegral_yield->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    gPad->SetLogy(1);
    hintegral_yield->GetXaxis()->SetRangeUser(-0.1, 15.2);
    hintegral_yield->Draw("pe");
    hintegral_yield->Write("yield_integral");
    // hintegral_yield->Write("yield");
    TLegend *legyield = new TLegend(0.8, 0.8, 0.91, 0.9);
    SetLegendStyle(legyield);
    legyield->SetTextSize(txtsize);
    // legyield->AddEntry(hYieldpar, "pbpb 5.36 TeV");
    // t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    legyield->Draw();
    csig->SaveAs((koutputfolder + "/yield_integral." + outputtype).c_str());
    csig->Clear();

    // Yield vs pT (bin counting method)
    hYbincount->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hYbincount->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    SetHistoQA(hYbincount);
    hYbincount->Draw("pe");
    hYbincount->Write("yield_bincount");
    csig->SaveAs((koutputfolder + "/yield_bincount." + outputtype).c_str());
    csig->Close();

    // Stop the stopwatch
    timer.Stop();

    // Get the elapsed time
    Double_t realTime = timer.RealTime(); // Wall clock time in seconds
    Double_t cpuTime = timer.CpuTime();   // CPU time used in seconds

    // Print the elapsed times
    std::cout << "Real time elapsed: " << realTime << " seconds" << std::endl;
    std::cout << "CPU time used: " << cpuTime << " seconds" << std::endl;
}
