#include "src/common.h"
#include "src/logger.h"
#include "src/fitfunc.h"
#include "src/DrawingHelper.C"
TH1 *GetNorBkg(TH1 *hSig, TH1 *hBkg, int currentptbin, double normStart, double normEnd);
TF1 *CreateFitFunction(bool useBkgPol3);
TF1 *CreateBackgroundFitFunction(bool useBkgPol3);
TF1 *CreateSignalFitFunction();
TH1F *GenerateSignalHistogram(TH3F *hData, int iCent, int iPt);
TH1F *GenerateBackgroundHistogram(TH3F *hData, int iCent, int iPt);
void Signal(int chosenDir = 0, int chosenAnti = -1, int chosenCent = -1, string sys = "", int sysvar = 0)
{
    gErrorIgnoreLevel = kError; // Suppressing warning outputs
    LogCustom log(TLogLevel::INFO);
    bool RunSpecificFolder(false), RunSpecificAnti(false), RunSpecificCent(false), RunFitSys(false);
    if (chosenDir > 0)
        RunSpecificFolder = true;
    if (chosenAnti >= 0)
        RunSpecificAnti = true;
    if (chosenCent >= 0)
        RunSpecificCent = true;
    if (sys.length() > 0)
        RunFitSys = true;
    TString sysOption = sys;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    TFile inputFile(kDataFilename.data());
    // Output file naming and figure folder
    TFile *outputFile;
    TString outputFileName = kSignalOutputPart.data();
    TString figurePath = kFiguresFolder.data();
    if ((chosenDir == 0) && (chosenAnti < 0) && (chosenCent < 0) && !RunFitSys)
    {
        outputFile = (TFile *)new TFile(kSignalOutput.data(), "recreate");
    }
    else
    {
        outputFileName += Form("_%d", chosenDir);
        figurePath += Form("/sys%d", chosenDir);
        if (chosenAnti != -1)
        {
            outputFileName += Form("_%d", chosenAnti);
            figurePath += Form("_part%d", chosenAnti);
        }
        if (chosenCent != -1)
            outputFileName += Form("_%d", chosenCent);
        if (RunFitSys)
        {
            outputFileName += Form("_%s_%d", sys.c_str(), sysvar);
            figurePath += Form("_%s_%d", sys.c_str(), sysvar);
        }
        outputFileName += ".root";
        outputFile = (TFile *)new TFile(outputFileName, "recreate");
    }
    figurePath += "/";
    log.Get(TLogLevel::INFO) << "Output file: " << outputFileName;
    log.Get(TLogLevel::INFO) << "Figure path: " << figurePath;

    // Preparation
    std::vector<double> pt_points_e;
    for (int i = 0; i < (int)kpTbin.size() - 1; i++)
        pt_points_e.push_back(kpTbin.at(i + 1) - kpTbin[i]);
    // Arrays
    std::vector<Double_t> fitmeans;
    std::vector<Double_t> fitmeans_err;
    std::vector<Double_t> fitsigmas;
    std::vector<Double_t> fitsigmas_err;
    std::vector<Double_t> fitchi2s;
    std::vector<Double_t> significances;
    std::vector<Double_t> RawYields;
    std::vector<Double_t> RawYields_err;
    std::vector<Double_t> RawYields_BC;
    std::vector<Double_t> RawYields_BC_err;
    // Canvas
    TCanvas *cInvMasspT = GetCanvas("cInvMasspT");
    TCanvas *cInvMasspTFit = GetCanvas("cInvMasspTFit");
    TCanvas *cQA = GetCanvas("cQA");
    TCanvas *cBigSig = GetCanvas("cBigSig", 1920, 1920);
    TCanvas *cBigFit = GetCanvas("cBigFit", 1920, 1920);
    // Divide cBig by kpTbin size
    int dividing = 0;
    if (kNPtBins > 9)
        dividing = 4;
    if (kNPtBins > 16)
        dividing = 5;
    if (kNPtBins > 25)
        dividing = 6;
    cBigSig->Divide(dividing, dividing, 0.001, 0.001);
    cBigFit->Divide(dividing, dividing, 0.001, 0.001);
    // Text
    TLatex *t2 = new TLatex();
    t2->SetNDC();
    t2->SetTextSize(0.04);
    // Fitting functions (see header for details)
    bool useBkgPol3 = false;
    if (sysOption.Contains("pol3")) // bkg function systematics
    {
        log.Get(TLogLevel::INFO) << "Used Pol(3) background";
        useBkgPol3 = true;
    }
    TF1 *fitFcn = CreateFitFunction(useBkgPol3);
    TF1 *fitFcn_bkg = CreateBackgroundFitFunction(useBkgPol3);
    TF1 *fitFcn_sig = CreateSignalFitFunction();

    TH3F *hSig_temp = nullptr; // temporary histogram for signal
    TH3F *hBkg_temp = nullptr; // temporary histogram for background
    Double_t bcError;          // Error of bin counting
    bool useLSbkg = true; // changed here false to true
    bool noBkgsubtraction = false;
    if (sysOption.Contains("LSBkg")) // Like-sign background
    {
        log.Get(TLogLevel::INFO) << "Used Like-sign background";
        useLSbkg = true;
    }
    if (sysOption.Contains("BkgFit")) // Like-sign background
    {
        log.Get(TLogLevel::INFO) << "Used direct background fitting without subtraction";
        noBkgsubtraction = true;
    }
    // Systematic directory name
    const TString SystematicDirectory = Form("%s%s", kFilterListNames.data(), SystematicBins[chosenDir].c_str());
    for (auto list_key : *inputFile.GetListOfKeys())
    {
        /// Preliminary operation to read the list and create an output dir
        log.Get(TLogLevel::INFO) << "Checking " << list_key->GetName();
        if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos)
        {
            log.Get(TLogLevel::INFO) << "  Analysis folder not found, stop processing";
            continue;
        }
        log.Get(TLogLevel::INFO) << "  Analysis folder found";
        if (RunSpecificFolder && string(list_key->GetName()) != SystematicDirectory) // For multi-core run
        {
            log.Get(TLogLevel::INFO) << "  Not a target folder, skip";
            continue;
        }
        log.Get(TLogLevel::INFO) << "  Start processing " << list_key->GetName();

        TDirectory *baseDir = outputFile->mkdir(list_key->GetName());
        // cout<<"\n\n\nThe name of the key is "<<list_key->GetName()<<"\n\n\n";
        outputFile->cd(list_key->GetName());
        string partname = "k892";
        auto th3_Signal = (TH3F *)inputFile.Get(Form("%s/h3%sinvmassDS", list_key->GetName(), partname.data()));
        auto th3_Signal_Anti = (TH3F *)inputFile.Get(Form("%s/h3%sinvmassDSAnti", list_key->GetName(), partname.data()));
        auto th3_LS = (TH3F *)inputFile.Get(Form("%s/h3%sinvmassLS", list_key->GetName(), partname.data()));
        auto th3_LS_Anti = (TH3F *)inputFile.Get(Form("%s/h3%sinvmassLSAnti", list_key->GetName(), partname.data()));

        auto th3_ME = (TH3F *)inputFile.Get(Form("%s/h3%sinvmassME", list_key->GetName(), partname.data()));

        std::vector<TH3F *> hData = {th3_Signal, th3_Signal_Anti};
        std::vector<TH3F *> hBkg_LS = {th3_LS, th3_LS_Anti};
        TH3F *hBkg_ME = (TH3F *)th3_ME->Clone("hBkg_ME");

        for (int iAnti = 0; iAnti < kAntiBinSize; iAnti++)
        {
            if (RunSpecificAnti && (chosenAnti != iAnti))
                continue;
            auto isAntiSum = (kAntiBinSize < 2) ? true : false;
            TDirectory *partDir = baseDir->mkdir(kParticleType[iAnti].c_str());
            partDir->cd();
            for (int iCent = 0; iCent < kNMultiplicityBins + 1; iCent++)
            {
                double multiStart{kMultiplicityBins[iCent - 1]}, multiEnd{kMultiplicityBins[iCent]};
                if (RunSpecificCent && (chosenCent != iCent))
                    continue;
                log.Get(TLogLevel::INFO) << "Cent bin: " << iCent;
                if (!kIsINEL && (iCent < 1))
                {
                    multiStart = 0.0;
                    multiEnd = 100.0;
                    log.Get(TLogLevel::INFO) << "MB MODE";
                }
                if (kIsINEL)
                {
                    log.Get(TLogLevel::INFO) << "INEL";
                    multiStart = 0.0;
                    multiEnd = 110.0;
                    if (iCent > 0)
                    {
                        log.Get(TLogLevel::INFO) << "Skip (INEL)";
                        continue; // skip INEL>0
                    }
                }
                log.Get(TLogLevel::INFO) << "Multiplicity percentile: " << multiStart << " - " << multiEnd;
                hSig_temp = (TH3F *)hData[iAnti]->Clone("hSig_temp");
                hBkg_temp = (TH3F *)hBkg_LS[iAnti]->Clone("hBkg_temp");
                if (isAntiSum)
                {
                    if (iAnti > 0)
                        continue; // skip anti particle if requested
                    log.Get(TLogLevel::INFO) << "Normal + Anti";
                    hSig_temp->Add(hData[abs(1 - iAnti)]); // 0: 1, 1: 0
                    hBkg_temp->Add(hBkg_LS[abs(1 - iAnti)]);
                }
                if (!useLSbkg)
                    hBkg_temp = (TH3F *)hBkg_ME->Clone("hBkg_ME"); // use ME as bkg

                fitmeans.clear();
                fitmeans_err.clear();
                fitsigmas.clear();
                fitsigmas_err.clear();
                fitchi2s.clear();
                significances.clear();
                RawYields.clear();
                RawYields_err.clear();
                RawYields_BC.clear();
                RawYields_BC_err.clear();

                TDirectory *centDir = partDir->mkdir(Form("%d", iCent));
                centDir->cd();
                for (Int_t iPt = 0, nPt = kpTbin.size() - 1; iPt < nPt; ++iPt)
                {
                    log.Get(TLogLevel::INFO) << "Pt bin: " << iPt << " (" << kpTbin[iPt] << ", " << kpTbin[iPt + 1] << ")";
                    TH1F *hSig = GenerateSignalHistogram(hSig_temp, iCent, iPt);
                    TH1F *hBkg_orig = GenerateBackgroundHistogram(hBkg_temp, iCent, iPt);
                    TH1F *hBkg = (TH1F *)hBkg_orig->Clone(Form("hBkg_%d", iPt));  //hBkg is a pointer to the cloned histogram hBkg_%d. simply clone() can be used without the arguments, then the compiler automatically gives the cloned histogram some name to which the hBkg points.

                    // Format
                    hSig->GetXaxis()->SetTitle("m_{inv} (GeV/#it{c}^{2})");
                    auto binWidth = (hSig->GetXaxis()->GetXmax() - hSig->GetXaxis()->GetXmin()) * kRebin / hSig->GetXaxis()->GetNbins();
                    log.Get(TLogLevel::INFO) << "Bin width: " << binWidth;
                    hBkg->SetLineColor(2);

                    std::vector<double> NormRange = {kNormRangepT[iPt][0], kNormRangepT[iPt][1]};
                    if (sysOption.Contains("normrange")) // Normalisation range systematics
                    {
                        log.Get(TLogLevel::INFO) << "Normalisation range systematics";
                        if (sysOption.Contains("L"))
                        {
                            NormRange[0] += sysvar * binWidth;
                            std::cout<<"\n\nThe normalization range value is "<<NormRange[0]<<"\n\n";
                            log.Get(TLogLevel::INFO) << "  Lower bound: " << NormRange[0];
                        }
                        if (sysOption.Contains("R"))
                        {
                            NormRange[1] += sysvar * binWidth;
                            log.Get(TLogLevel::INFO) << "  Upper bound: " << NormRange[1];
                        }
                    }
                    log.Get(TLogLevel::INFO) << "Normalisation range: " << NormRange[0] << " - " << NormRange[1];
                    TH1F *hBkg_norm = (TH1F *)GetNorBkg(hSig, hBkg_orig, iPt, NormRange[0], NormRange[1]);

                    // rebin
                    hSig->Rebin(kRebin);
                    hBkg->Rebin(kRebin);
                    hBkg_norm->Rebin(kRebin);

                    // scale region
                    auto hBkg_norm_scale = (TH1F *)hBkg_norm->Clone(Form("hBkg_norm_scale_%d", iPt));
                    for (int iBin = 1; iBin <= hBkg_norm_scale->GetNbinsX(); iBin++)
                    {
                        if (hBkg_norm_scale->GetBinCenter(iBin) < kNormRangepT[iPt][0] || hBkg_norm_scale->GetBinCenter(iBin) > kNormRangepT[iPt][1])
                            hBkg_norm_scale->SetBinContent(iBin, -999);  // Value of bins other than normalization range is set to -999 such that only the normalization range is visible to be plotted as a bar graph.
                    }
                    hBkg_norm_scale->SetLineColor(kRed);
                    hBkg_norm_scale->SetMarkerColor(kRed);
                    hBkg_norm_scale->SetFillColor(kRed);
                    hBkg_norm_scale->SetFillStyle(3001);

                    hSig->Write(Form("hSignal_%i_%i", iCent, iPt));
                    hBkg->Write(Form("hBkg_%i_%i", iCent, iPt));
                    hBkg_norm->Write(Form("hBkg_norm_%i_%i", iCent, iPt));
                    hBkg_norm_scale->Write(Form("hBkg_norm_scale_%i_%i", iCent, iPt));

                    hSig->GetXaxis()->SetRangeUser(kDrawRange[0], kDrawRange[1]);
                    hSig->SetMaximum(hSig->GetMaximum() * 1.5);
                    cInvMasspT->cd();
                    hSig->Draw("E");
                    // hBkg->Draw("E same");
                    hBkg_norm->Draw("E same");
                    hBkg_norm_scale->Draw("BAR same");

                    auto legend = new TLegend(0.17, 0.7, 0.37, 0.9);
                    legend->SetFillStyle(0);
                    legend->SetBorderSize(0);
                    legend->SetHeader(Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", kpTbin[iPt], kpTbin[iPt + 1]));
                    legend->AddEntry(hSig, "data", "LE");
                    legend->AddEntry(hBkg, "Bkg", "LE");
                    // legend->AddEntry(temp, "Normalization Region", "F");
                    legend->Draw();
                    t2->DrawLatex(0.15, 0.95, "#bf{K(892)^{0} #rightarrow #pi + K}");
                    t2->DrawLatex(0.5, 0.95, Form("#bf{%.1f < T0M percentile < %.1f (%%)}", multiStart, multiEnd));
                    // t2->DrawLatex(0.75, 0.85, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", kpTbin[iPt], kpTbin[iPt + 1]));
                    cInvMasspT->Write(Form("cInvMasspT_%i_%i", iCent, iPt));
                    SaveCanvas(cInvMasspT, Form("cInvMasspT_cen%i_pT%i", iCent, iPt), Form("%s/%d/", figurePath.Data(), iCent), kOutputType.data());
                    // cInvMasspT->SaveAs(Form("%s/cInvMasspT_cen%i_%i.png", figurePath.Data(), iCent, iPt));
                    cBigSig->cd(iPt + 1);
                    hSig->Draw("E");
                    // hBkg->Draw("E same");
                    hBkg_norm->Draw("E same");
                    hBkg_norm_scale->Draw("BAR same");
                    legend->Draw();
                    t2->DrawLatex(0.15, 0.95, "#bf{K(892)^{0} #rightarrow #pi + K}");
                    t2->DrawLatex(0.5, 0.95, Form("#bf{%.1f < T0M percentile < %.1f (%%)}", multiStart, multiEnd));
                    // t2->DrawLatex(0.75, 0.85, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", kpTbin[iPt], kpTbin[iPt + 1]));

                    // Subtract background
                    TH1F *hSig_sub = (TH1F *)hSig->Clone(Form("hSig_sub_%i_%i", iCent, iPt));
                    if (!noBkgsubtraction)
                        hSig_sub->Add(hBkg_norm, -1);
                    hSig_sub->Write(Form("hSig_sub_%i_%i", iCent, iPt));

                    // Fit
                    std::vector<double> fitRange = {kFitRange[iPt][0], kFitRange[iPt][1]};
                    if (sysOption.Contains("fitrange")) // Fit range systematics
                    {
                        log.Get(TLogLevel::INFO) << "Fit range systematics";
                        if (sysOption.Contains("L"))
                        {
                            fitRange[0] += sysvar * binWidth;
                            log.Get(TLogLevel::INFO) << "  Lower bound: " << fitRange[0];
                        }
                        if (sysOption.Contains("R"))
                        {
                            fitRange[1] += sysvar * binWidth;
                            log.Get(TLogLevel::INFO) << "  Upper bound: " << fitRange[1];
                        }
                    }
                    if (useBkgPol3)
                        fitFcn->SetParLimits(3, kFitPol3ParLimit[iPt][0], kFitPol3ParLimit[iPt][1]);
                    log.Get(TLogLevel::INFO) << "Fit range: " << fitRange[0] << " - " << fitRange[1];

                    // Real fitting procedure
                    // auto r = hSig_sub->Fit(fitFcn, "QREBMS0+", "", fitRange[0], fitRange[1]); // original fit method (simple)

                    fitFcn->SetRange(fitRange[0], fitRange[1]);
                    fitFcn_sig->SetRange(fitRange[0], fitRange[1]);
                    fitFcn_bkg->SetRange(fitRange[0], fitRange[1]);

                    TH1 *hNoPeak = (TH1 *)hSig_sub->Clone(Form("hNoPeak_%i_%i", iCent, iPt));
                    for (int iBin = hSig_sub->GetXaxis()->FindBin(1.000001 * kPeakRange[0]); iBin <= hSig_sub->GetXaxis()->FindBin(0.999999 * kPeakRange[1]); iBin++)
                    {
                        hNoPeak->SetBinContent(iBin, 0);
                        hNoPeak->SetBinError(iBin, 0);
                    }
                    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // better minimizer!
                    auto fitFcn_bkg_tempfit = (TF1 *)fitFcn_bkg->Clone("fitFcn_bkg_tempfit");
                    fitFcn_bkg_tempfit->SetRange(fitRange[0], fitRange[1]);
                    hNoPeak->Fit(fitFcn_bkg_tempfit, "QREBMS0+", "", fitRange[0], fitRange[1]);
                    // TCanvas *c = new TCanvas();
                    // hNoPeak->Draw(); // to show the fit remove 0 from the fit type
                    // c->SaveAs(Form("Nopeakfit_%d.png", iPt));
                    hNoPeak->Write(Form("hNoPeak_%i_%i", iCent, iPt));
                    fitFcn_bkg_tempfit->Write(Form("fitFcn_bkg_tempfit%i_%i", iCent, iPt));
                    // fitFcn->SetParameter(3, fitFcn_bkg_tempfit->GetParameter(0));
                    // fitFcn->SetParameter(4, fitFcn_bkg_tempfit->GetParameter(1));
                    // fitFcn->SetParameter(5, fitFcn_bkg_tempfit->GetParameter(2));
                    // if (useBkgPol3)
                    //     fitFcn->SetParameter(6, fitFcn_bkg_tempfit->GetParameter(3));

                    auto r = hSig_sub->Fit(fitFcn, "REBMS0+", "", fitRange[0], fitRange[1]); // original fit method (simple)
                    // Get fit results
                    Double_t *par = fitFcn->GetParameters();
                    fitFcn_sig->SetParameters(&par[0]);
                    fitFcn_sig->SetParError(0, fitFcn->GetParError(0));
                    fitFcn_sig->SetParError(1, fitFcn->GetParError(1));
                    fitFcn_sig->SetParError(2, fitFcn->GetParError(2));
                    fitFcn_bkg->SetParameters(&par[3]);
                    fitFcn_bkg->SetParError(0, fitFcn->GetParError(3));
                    fitFcn_bkg->SetParError(1, fitFcn->GetParError(4));
                    fitFcn_bkg->SetParError(2, fitFcn->GetParError(5));
                    if (useBkgPol3)
                        fitFcn_bkg->SetParError(3, fitFcn->GetParError(6));

                    // cov
                    TMatrixDSym cov = r->GetCovarianceMatrix();
                    TMatrixDSym cov1;
                    TMatrixDSym cov2;
                    cov.GetSub(0, 2, 0, 2, cov1);
                    cov.GetSub(3, 5, 3, 5, cov2);
                    Double_t *b = cov1.GetMatrixArray();
                    Double_t *a = cov2.GetMatrixArray();

                    // fit results
                    auto Mass = fitFcn->GetParameter(0);
                    auto Width = fitFcn->GetParameter(1);
                    auto Yield = fitFcn->GetParameter(2);
                    auto poly2 = fitFcn->GetParameter(3);
                    auto poly1 = fitFcn->GetParameter(4);
                    auto poly0 = fitFcn->GetParameter(5);
                    if (!useBkgPol3)
                        log.Get(TLogLevel::INFO) << "Fit parameters: " << Mass << " " << Width << " " << Yield << " " << poly2 << " " << poly1 << " " << poly0;
                    else
                    {
                        auto poly3 = fitFcn->GetParameter(6);
                        log.Get(TLogLevel::INFO) << "Fit parameters: " << Mass << " " << Width << " " << Yield << " " << poly2 << " " << poly1 << " " << poly0 << " " << poly3;
                    }

                    auto ErrorMass = fitFcn->GetParError(0);
                    auto ErrorWidth = fitFcn->GetParError(1);
                    auto ErrorYield = fitFcn->GetParError(2);
                    auto Chi2Ndf = (fitFcn->GetChisquare()) / (fitFcn->GetNDF());

                    // data info
                    auto yStart = masspdg - 3 * widthpdg;
                    auto yEnd = masspdg + 3 * widthpdg;
                    // if (fitRange[0] > yStart)
                    //     yStart = fitRange[0];
                    // if (fitRange[1] < yEnd)
                    //     yEnd = fitRange[1];
                    auto bmin = hSig_sub->GetXaxis()->FindBin(yStart);
                    auto bmax = hSig_sub->GetXaxis()->FindBin(yEnd);
                    auto significance_den = TMath::Sqrt(hSig->Integral(bmin, bmax));
                    auto significance_num = (fitFcn_sig->Integral(yStart, yEnd)) / binWidth;
                    auto significance = significance_num / significance_den;

                    auto rawYields_fit = fitFcn_sig->Integral(yStart, yEnd) / binWidth;
                    auto rawYields_fit_err = fitFcn_sig->IntegralError(yStart, yEnd, &par[0], b) / binWidth;

                    // bin counting method
                    auto rawYields_bincount_temp = hSig_sub->IntegralAndError(bmin, bmax, bcError);
                    auto rawYields_bincount_bkg = fitFcn_bkg->Integral(yStart, yEnd) / binWidth;
                    auto rawYields_bincount_bkg_err = fitFcn_bkg->IntegralError(yStart, yEnd, &par[3], a) / binWidth;
                    auto rawYields_bincount = rawYields_bincount_temp - rawYields_bincount_bkg; // background subtraction
                    auto rawYields_bincount_err = TMath::Sqrt(bcError * bcError + rawYields_bincount_bkg_err * rawYields_bincount_bkg_err);

                    log.Get(TLogLevel::INFO) << "Bincoudnt yields: " << rawYields_bincount_temp << " fit bkg yields: " << rawYields_bincount_bkg << " fit yields: " << rawYields_fit;

                    // Fill results
                    fitmeans.push_back(Mass);
                    fitmeans_err.push_back(ErrorMass);
                    fitsigmas.push_back(Width);
                    fitsigmas_err.push_back(ErrorWidth);
                    fitchi2s.push_back(Chi2Ndf);
                    significances.push_back(significance);
                    RawYields.push_back(rawYields_fit);
                    RawYields_err.push_back(rawYields_fit_err);
                    RawYields_BC.push_back(rawYields_bincount);
                    RawYields_BC_err.push_back(rawYields_bincount_err);

                    // fit legend
                    auto legend_fit = new TLegend(0.5, 0.35, 0.9, 0.9);
                    legend_fit->SetTextFont(42);
                    legend_fit->SetTextSize(0.03);
                    legend_fit->SetFillStyle(0);
                    legend_fit->SetBorderSize(0);
                    legend_fit->SetHeader(Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", kpTbin[iPt], kpTbin[iPt + 1]));
                    legend_fit->AddEntry(hSig_sub, "data", "LE");
                    legend_fit->AddEntry(fitFcn, "fit", "L");
                    legend_fit->AddEntry(fitFcn_sig, "signal", "L");
                    legend_fit->AddEntry(fitFcn_bkg, "background", "L");
                    legend_fit->AddEntry((TObject *)0, Form("Mass = %.3f #pm %.3f GeV/#it{c}^{2}", Mass, ErrorMass), "");
                    legend_fit->AddEntry((TObject *)0, Form("Width = %.3f #pm %.3f GeV/#it{c}^{2}", Width, ErrorWidth), "");
                    legend_fit->AddEntry((TObject *)0, Form("Yield(BWfit) = %.0f #pm %.0f", rawYields_fit, rawYields_fit_err), "");
                    legend_fit->AddEntry((TObject *)0, Form("Yield(BC) = %.0f #pm %.0f", rawYields_bincount, rawYields_bincount_err), "");
                    legend_fit->AddEntry((TObject *)0, Form("Significance = %.2f", significance), "");
                    legend_fit->AddEntry((TObject *)0, Form("#chi^{2}/NDF = %.2f", Chi2Ndf), "");

                    // Draw
                    if (!noBkgsubtraction) {
                        hSig_sub->GetXaxis()->SetRangeUser(kDrawFitRange[0], kDrawFitRange[1]);
                        hSig_sub->SetMaximum(hSig_sub->GetMaximum() * 1.5);
                    }
                    else {
                        hSig_sub->GetXaxis()->SetRangeUser(fitRange[0], fitRange[1]);
                        hSig_sub->SetMaximum(hSig_sub->GetBinContent(hSig_sub->GetXaxis()->FindBin(masspdg)) * 1.1);
                    }
                    cInvMasspTFit->cd();
                    hSig_sub->Draw("E");
                    fitFcn->Draw("same");
                    fitFcn_sig->Draw("same");
                    fitFcn_bkg->Draw("same");
                    legend_fit->Draw();
                    t2->DrawLatex(0.15, 0.95, "#bf{K(892)^{0} #rightarrow #pi + K}");
                    t2->DrawLatex(0.5, 0.92, Form("#bf{%.1f < T0M percentile < %.1f (%%)}", multiStart, multiEnd));

                    cInvMasspTFit->Write(Form("cInvMasspTFit_%i_%i", iCent, iPt));
                    SaveCanvas(cInvMasspTFit, Form("cInvMasspTFit_cen%i_pT%i", iCent, iPt), Form("%s/%d/", figurePath.Data(), iCent), kOutputType.data());
                    // cInvMasspTFit->SaveAs(Form("%s/cInvMasspTFit_cen%i_%i.png", figurePath.Data(), iCent, iPt));

                    auto fitFcn_temp = fitFcn->Clone();
                    auto fitFcn_sig_temp = fitFcn_sig->Clone();
                    auto fitFcn_bkg_temp = fitFcn_bkg->Clone();
                    cBigFit->cd(iPt + 1);
                    hSig_sub->Draw("E");
                    fitFcn_temp->Draw("same");
                    fitFcn_sig_temp->Draw("same");
                    fitFcn_bkg_temp->Draw("same");
                    legend_fit->Draw();
                    t2->DrawLatex(0.15, 0.95, "#bf{K(892)^{0} #rightarrow #pi + K}");
                    t2->DrawLatex(0.5, 0.95, Form("#bf{%.1f < T0M percentile < %.1f (%%)}", multiStart, multiEnd));

                    // Fit
                } // end of loop over pT bins
                SaveCanvas(cBigSig, Form("cBigSig_cen%i", iCent), figurePath, kOutputType.data());
                SaveCanvas(cBigFit, Form("cBigFit_cen%i", iCent), figurePath, kOutputType.data());

                auto hRawYield = MakeHistfromArray("RawYields", RawYields, kpTbin, RawYields_err);
                hRawYield->SetMarkerStyle(20);
                hRawYield->SetMarkerColor(kBlack);
                hRawYield->SetLineColor(kBlack);
                hRawYield->Write("hRawYields");
                cQA->cd();
                hRawYield->Draw("E");
                cQA->Write("cRawYields");
                // cQA->SaveAs(Form("%s/cRawYields_cen%i.png", Form("%s/%d/", figurePath.Data(), iCent), iCent));
                SaveCanvas(cQA, Form("cRawYields_cen%i", iCent), Form("%s/%d/", figurePath.Data(), iCent), kOutputType.data());

                auto hRawYieldBC = MakeHistfromArray("RawYieldsBC", RawYields_BC, kpTbin, RawYields_BC_err);
                hRawYieldBC->SetMarkerStyle(20);
                hRawYieldBC->SetMarkerColor(kRed);
                hRawYieldBC->SetLineColor(kRed);
                hRawYieldBC->Write("hRawYieldsBC");
                cQA->cd();
                hRawYieldBC->Draw("E");
                cQA->Write("cRawYieldsBC");
                // cQA->SaveAs(Form("%s/cRawYieldsBC_cen%i.png", Form("%s/%d/", figurePath.Data(), iCent), iCent));
                SaveCanvas(cQA, Form("cRawYieldsBC_cen%i", iCent), Form("%s/%d/", figurePath.Data(), iCent), kOutputType.data());

                cQA->cd();
                hRawYield->SetMinimum(0);
                hRawYield->Draw("E");
                hRawYieldBC->Draw("E same");
                cQA->Write("cRawYieldsCompare");
                // cQA->SaveAs(Form("%s/cRawYieldsCompare_cen%i.png", Form("%s/%d/", figurePath.Data(), iCent), iCent));
                SaveCanvas(cQA, Form("cRawYieldsCompare_cen%i", iCent), Form("%s/%d/", figurePath.Data(), iCent), kOutputType.data());

                auto hFitMean = MakeHistfromArray("Fit means", fitmeans, kpTbin, fitmeans_err);
                hFitMean->Write("hFitMean");
                TF1 *pdgmass = new TF1("pdgmass", "[0]", 0, 20);
                pdgmass->SetParameter(0, masspdg);
                pdgmass->SetLineColor(kRed);
                pdgmass->SetLineStyle(2);
                cQA->cd();
                hFitMean->Draw("E");
                pdgmass->Draw("same");
                cQA->Write("cFitMean");
                // cQA->SaveAs(Form("%s/cFitMean_cen%i.png", Form("%s/%d/", figurePath.Data(), iCent), iCent));
                SaveCanvas(cQA, Form("cFitMean_cen%i", iCent), Form("%s/%d/", figurePath.Data(), iCent), kOutputType.data());

                auto hFitSigma = MakeHistfromArray("Fit sigmas", fitsigmas, kpTbin, fitsigmas_err);
                hFitSigma->Write("hFitSigma");
                TF1 *pdgwidth = new TF1("pdgwidth", "[0]", 0, 20);
                pdgwidth->SetParameter(0, widthpdg);
                pdgwidth->SetLineColor(kRed);
                pdgwidth->SetLineStyle(2);
                cQA->cd();
                hFitSigma->Draw("E");
                pdgwidth->Draw("same");
                cQA->Write("cFitSigma");
                // cQA->SaveAs(Form("%s/cFitSigma_cen%i.png", Form("%s/%d/", figurePath.Data(), iCent), iCent));
                SaveCanvas(cQA, Form("cFitSigma_cen%i", iCent), Form("%s/%d/", figurePath.Data(), iCent), kOutputType.data());

                auto hFitChi2 = MakeHistfromArray("Fit chi2s", fitchi2s, kpTbin);
                hFitChi2->Write("hFitChi2");
                cQA->cd();
                hFitChi2->Draw("E");
                cQA->Write("cFitChi2");
                // cQA->SaveAs(Form("%s/cFitChi2_cen%i.png", Form("%s/%d/", figurePath.Data(), iCent), iCent));
                SaveCanvas(cQA, Form("cFitChi2_cen%i", iCent), Form("%s/%d/", figurePath.Data(), iCent), kOutputType.data());

                auto hSignificance = MakeHistfromArray("Significance", significances, kpTbin);
                hSignificance->Write("hSignificance");
                cQA->cd();
                hSignificance->Draw("E");
                cQA->Write("cSignificance");
                // cQA->SaveAs(Form("%s/cSignificance_cen%i.png", Form("%s/%d/", figurePath.Data(), iCent), iCent));
                SaveCanvas(cQA, Form("cSignificance_cen%i", iCent), Form("%s/%d/", figurePath.Data(), iCent), kOutputType.data());

            } // end of loop over centrality bins
        }     // end of loop over anti particles
        break;
    } // end of loop over directories

    // th3_Signal->Write();
    outputFile->Close();
}

TH1 *GetNorBkg(TH1 *hSig, TH1 *hBkg, int currentptbin, double normStart, double normEnd)
{
    TH1 *hBkgTemp = (TH1 *)hBkg->Clone();
    Double_t normalization_data = 0.;
    Double_t normalization_mixed = 0.;

    normalization_data +=
        hSig->Integral(hSig->GetXaxis()->FindBin(normStart),
                       hSig->GetXaxis()->FindBin(normEnd));
    normalization_mixed +=
        hBkg->Integral(hBkg->GetXaxis()->FindBin(normStart),
                       hBkg->GetXaxis()->FindBin(normEnd));
    hBkgTemp->Scale(normalization_data / normalization_mixed);
    hBkgTemp->SetLineColor(kRed);
    hBkgTemp->SetMarkerColor(kRed);

    return hBkgTemp;
}

TF1 *CreateFitFunction(bool useBkgPol3)
{
    //Initially fixed range is used to define the fit function, but in main code it is changed for different pT using SetRange method.
    TF1 *fitFcn;
    if (useBkgPol3)
        fitFcn = new TF1("fitfunc", BreitWignerpoly3, kFitRange[0][0], kFitRange[0][1], 7); // sig+bkg fit function
    else
        fitFcn = new TF1("fitfunc", BreitWignerpoly2, kFitRange[0][0], kFitRange[0][1], 6); // sig+bkg fit function
    fitFcn->SetParLimits(0, 0.83, 0.92);                                                    // Mass
    fitFcn->SetParameter(0, masspdg);                                                       // Mass
    fitFcn->SetParameter(1, widthpdg);                                                      // width
    fitFcn->FixParameter(1, widthpdg);                                                      // width
    fitFcn->SetParameter(2, 30000);                                                         // yield
    fitFcn->SetParLimits(2, 1, 1e8);                                                 // Yield
    fitFcn->SetParNames("Mass", "Width", "Yield", "C", "B", "A");
    fitFcn->SetLineColor(kRed);
    fitFcn->SetLineStyle(1);
    fitFcn->SetLineWidth(2);

    return fitFcn;
}
TF1 *CreateBackgroundFitFunction(bool useBkgPol3)
{
    TF1 *fitFcn_bkg;
    if (useBkgPol3)
        fitFcn_bkg = new TF1("fitfunc_bkg", polynomial3, kFitRange[0][0], kFitRange[0][1], 4); // bkg fit function
    else
        fitFcn_bkg = new TF1("fitfunc_bkg", polynomial2, kFitRange[0][0], kFitRange[0][1], 3); // bkg fit function
    fitFcn_bkg->SetParNames("C", "B", "A");
    fitFcn_bkg->SetLineColor(kGreen);
    fitFcn_bkg->SetLineStyle(2);
    fitFcn_bkg->SetLineWidth(2);

    return fitFcn_bkg;
}
TF1 *CreateSignalFitFunction()
{
    TF1 *fitFcn_sig = new TF1("fitfunc_sig", BW, kFitRange[0][0], kFitRange[0][1], 3); // sig fit function
    fitFcn_sig->SetParNames("Mass", "Width", "Yield");
    fitFcn_sig->SetLineColor(kBlue);
    fitFcn_sig->SetLineStyle(2);
    fitFcn_sig->SetLineWidth(2);

    return fitFcn_sig;
}
TH1F *GenerateSignalHistogram(TH3F *hData, int iCent, int iPt)
{

    //Here the conditional (ternary) operator is used i.e. "? :" to create if else statement. The number 0.001 is used to avoid overflow and underflow bins
    TH1F *hSig = (iCent == 0) ? (TH1F *)hData->ProjectionZ( // MB case
                                    Form("hSig_%d", iPt), -1, -1, hData->GetYaxis()->FindBin(kpTbin[iPt] + 0.001),
                                    hData->GetYaxis()->FindBin(kpTbin[iPt + 1] - 0.001)) // from the 3D histogram hSig
                              : (TH1F *)hData->ProjectionZ(                              // Centrality dependent case
                                    Form("hSig_%d", iPt), hData->GetXaxis()->FindBin(kMultiplicityBins[iCent - 1] + 0.001),
                                    hData->GetXaxis()->FindBin(kMultiplicityBins[iCent] - 0.001), hData->GetYaxis()->FindBin(kpTbin[iPt] + 0.001),
                                    hData->GetYaxis()->FindBin(kpTbin[iPt + 1] - 0.001)); // from the 3D histogram hSig
    hSig->SetTitle(Form("hSig_%d_%d", iCent, iPt));
    hSig->SetName(Form("hSig_%d_%d", iCent, iPt));
    hSig->SetMarkerStyle(20);
    hSig->SetMarkerColor(kBlack);
    hSig->SetLineColor(kBlack);
    hSig->SetMarkerSize(1);
    hSig->SetLineWidth(2);
    hSig->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hSig->GetYaxis()->SetTitle("Counts");
    hSig->GetXaxis()->SetRangeUser(kDrawRange[0], kDrawRange[1]);

    return hSig;
}
TH1F *GenerateBackgroundHistogram(TH3F *hData, int iCent, int iPt)
{
    TH1F *hBkg = (iCent == 0) ? (TH1F *)hData->ProjectionZ( // MB case
                                    Form("hBkg_%d", iPt), -1, -1, hData->GetYaxis()->FindBin(kpTbin[iPt] + 0.001),
                                    hData->GetYaxis()->FindBin(kpTbin[iPt + 1] - 0.001)) // from the 3D histogram hBkg
                              : (TH1F *)hData->ProjectionZ(                              // Centrality dependent case
                                    Form("hBkg_%d", iPt), hData->GetXaxis()->FindBin(kMultiplicityBins[iCent - 1] + 0.001),
                                    hData->GetXaxis()->FindBin(kMultiplicityBins[iCent] - 0.001), hData->GetYaxis()->FindBin(kpTbin[iPt] + 0.001),
                                    hData->GetYaxis()->FindBin(kpTbin[iPt + 1] - 0.001)); // from the 3D histogram hBkg
    hBkg->SetTitle(Form("hBkg_%d_%d", iCent, iPt));
    hBkg->SetName(Form("hBkg_%d_%d", iCent, iPt));
    hBkg->SetMarkerStyle(20);
    hBkg->SetMarkerColor(kRed);
    hBkg->SetLineColor(kRed);
    hBkg->SetMarkerSize(1);
    hBkg->SetLineWidth(2);
    hBkg->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hBkg->GetYaxis()->SetTitle("Counts");
    hBkg->GetXaxis()->SetRangeUser(kDrawRange[0], kDrawRange[1]);

    return hBkg;
}