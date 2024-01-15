#include <AliPWGFunc.h>
#include <AliPWGHistoTools.h>

#include "AdditionalFunctions.h"
#include "src/YieldMeanNew.C"
#include "src/common.h"
#include "src/DrawingHelper.C"
#include "src/logger.h"

vector<TH1 *> GetSysSpectra(vector<vector<double>> multibin);
vector<TH1 *> GetSysSpectraNocor(vector<vector<double>> multibin);
vector<TH1 *> GetStatSpectra(vector<vector<double>> multibin);
TH1 *GetSpectrasys(double multi_start, double multi_end);
TH1 *GetSpectrasysNocor(double multi_start, double multi_end);
TH1 *GetSpectrastat(double multi_start, double multi_end);

TString multistring;
bool isINEL = false;
bool skipcustomfit = false;
enum
{
    kFitExpPt = 1,
    kFitLevi,
    fFitExpMt,
    kFitBoltzmann,
    kFitBlastWave,
    kFitBoseEinstein,
    kFitFermiDirac
};

std::vector<TString> functions = {"",
                                  "kFitExpPt",
                                  "kFitLevi",
                                  "fFitExpMt",
                                  "kFitBoltzmann",
                                  "kFitBlastWave",
                                  "kFitBoseEinstein",
                                  "kFitFermiDirac"};
vector<vector<vector<double>>> // { {first pT bin variation}, {last pTbin variation} } // for each function
    functionFitRangeVar = {
        {{}, {}},                            // null
        {{0.0}, {2, 2.5, 3, 5.0}}, // ExpPt, only front part.
        {{0.0}, {2, 2.5, 3, 5.0}}, // Levy, all region.
        {{0.0}, {2, 2.5, 3, 5.0}}, // ExpMt, only front part.
        {{0.0}, {2, 2.5, 3, 5.0}}, // Boltzmann, only front part.
        {{0.0}, {2, 2.5, 3, 5.0}}, // BGBW, all region.
        {{0.0}, {2, 2.5, 3, 5.0}}, // BoseEinstein, only front part.
        {{0.0}, {2, 2.5, 3, 5.0}}  // FermiDirac, only front part.
};
vector<vector<double>> centralvalue = {
    {0.0, 5.0}, // 0: MB
    {0.0, 5.0}, // 1: 0-10
    {0.0, 5.0}, // 2: 10-20
    {0.0, 5.0}, // 3: 20-30
    {0.0, 5.0}, // 4: 30-40
    {0.0, 5.0}, // 5: 40-50
    {0.0, 5.0}, // 6: 50-60
    {0.0, 5.0}, // 7: 60-70
    {0.0, 5.0}, // 8: 70-80
    {0.0, 5.0}, // 9: 80-90
    {0.0, 5.0}  // 10: 90-100
};
Int_t maxtrial = 10000;
TString fitoption = "0qiS"; // default "0q"
void YieldMeanK892(int chosenAnti = -1, int chosenCent = -1)
{
    LogCustom log(TLogLevel::INFO);
    if (kIsINEL)
    {
        centralvalue[0] = {0.0, 5.0};
    }
    int centralvaluebin = chosenCent;
    bool RunSpecificAnti(false), RunSpecificCent(false);
    if (chosenAnti >= 0)
        RunSpecificAnti = true;
    if (chosenCent >= 0)
        RunSpecificCent = true;

    TString outputFileName = kYieldMeanOutputPart.data();
    TString inputSysErrFileName = kSystematicsOutput.data();
    TString inputStatErrFileName = kSpectraOutputPart.data();
    if (!RunSpecificAnti && !RunSpecificCent)
    {
        outputFileName = kYieldMeanOutput.data();
        inputStatErrFileName = kSpectraOutput.data();
    }
    else
    {
        if (chosenAnti != -1)
        {
            outputFileName += Form("_%d", chosenAnti);
            inputStatErrFileName += Form("_%d", chosenAnti);
        }
        if (chosenCent != -1)
        {
            outputFileName += Form("_%d", chosenCent);
            inputStatErrFileName += "_1";
        }
        outputFileName += ".root";
        inputStatErrFileName += ".root";
    }
    TFile *outputFile = (TFile *)new TFile(outputFileName, "recreate");
    TFile SysErrFile(inputSysErrFileName);
    TFile StatErrFile(inputStatErrFileName);

    log.Get(TLogLevel::INFO) << "Input spectra file: " << inputStatErrFileName.Data();
    log.Get(TLogLevel::INFO) << "Input systematics file: " << inputSysErrFileName.Data();
    log.Get(TLogLevel::INFO) << "Output file: " << outputFileName.Data();

    TCanvas *cResults = new TCanvas("cResults", "cResults", 960, 720);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    cResults->SetTickx();
    cResults->SetTicky();
    cResults->SetTopMargin(0.05);
    cResults->SetLeftMargin(0.10);
    // cSigbkg->SetBottomMargin(0.01);
    cResults->SetRightMargin(0.01);
    cResults->SetFillStyle(0);

    AliPWGFunc *fm = new AliPWGFunc;
    fm->SetVarType(AliPWGFunc::VarType_t(AliPWGFunc::kdNdpt));
    TF1 *func = 0;
    int slopePar = 0;
    double fitchi2, fitNDF;
    TString nameinput;

    // Skipping some configuration
    // Multibin Customizing
    std::vector<vector<double>> functionFitCustom = {
        {0}, // null
        {0}, // ExpPt, only front part.
        {0}, // Levy, all region.
        {0}, // ExpMt, only front part.
        {0}, // Boltzmann, only front part.
        {0}, // BGBW, all region.
        {0}, // BoseEinstein, only front part.
        {0}  // FermiDirac, only front part.
    };

    for (int iAnti = 0; iAnti < kAntiBinSize; iAnti++)
    {
        if (RunSpecificAnti && (chosenAnti != iAnti))
            continue;
        auto isAntiSum = (kAntiBinSize < 2) ? true : false;
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
                multiEnd = 100.0;
                if (iCent > 0)
                {
                    log.Get(TLogLevel::INFO) << "Skip (INEL)";
                    continue; // skip INEL>0
                }
            }
            TH1 *hspectra_sys = (TH1 *)SysErrFile.Get(Form("hSystematic_%d_%d", iAnti, iCent));
            TH1 *hspectra_sys_nocor = (TH1 *)SysErrFile.Get(Form("hSystematic_NoCor_%d_%d", iAnti, iCent));
            TH1 *hspectra_sys_cor = (TH1 *)SysErrFile.Get(Form("hSystematic_Cor_%d_%d", iAnti, iCent));
            TH1 *hspectra_sys_multicor = (TH1 *)SysErrFile.Get(Form("hSystematic_MultiCor_%d_%d", iAnti, iCent));
            TH1 *hspectra_sys_multiuncor = (TH1 *)SysErrFile.Get(Form("hSystematic_MultiUnCor_%d_%d", iAnti, iCent));
            TH1 *hspectra_stat = (TH1 *)StatErrFile.Get(Form("%s/%s/%i/hCorrectedYields", kFilterListNames.data(), kParticleType[iAnti].c_str(), iCent));

            cResults->cd();
            cResults->Draw();
            cResults->SetLogy();

            std::vector<double> yield;
            std::vector<double> yield_state;
            std::vector<double> meanpt;
            std::vector<double> meanpt_state;
            std::vector<double> extraYields;

            std::vector<TF1 *> funcs;
            std::vector<TString> fitnames;

            int flinestyle = 1;
            TF1 *fout;
            for (int fitvar = 1; fitvar < functions.size(); fitvar++)
            {
                // if(fitvar > 2) continue;
                flinestyle = 1;
                for (auto const &fitstart : functionFitRangeVar[fitvar][0])
                { // fit bin start
                    for (auto const &fitend : functionFitRangeVar[fitvar][1])
                    { // fit bin end
                        nameinput = Form("%s %.2f - %.2f GeV/#it{c}", functions[fitvar].Data(), fitstart, fitend);
                        fitnames.push_back(nameinput);
                        if (fitvar == kFitLevi)
                        {
                            func = fm->GetLevi(masspdg, 0.4, 750, 3);
                            func->SetParLimits(1, 0.0001, 20000);
                            slopePar = 2;
                            // fitrange = {0.8,8.8};
                        }
                        if (fitvar == kFitExpPt)
                        {
                            func = fm->GetPTExp(0.2, 20);
                            slopePar = 1;
                            // fitrange = {0.8,1.6};
                        }
                        if (fitvar == fFitExpMt)
                        {
                            func = fm->GetMTExp(masspdg, 0.2, 20);
                            slopePar = 1;
                        }
                        if (fitvar == kFitBoltzmann)
                        {
                            func = fm->GetBoltzmann(masspdg, 0.2, 20);
                            slopePar = 1;
                            // fitrange = {0.8,1.6};
                        }
                        if (fitvar == kFitBlastWave)
                        {
                            func = fm->GetBGBW(masspdg, 0.6, 0.3, 1, 1e5); // beta, T, n, norm
                            func->SetParLimits(1, 0.1, 0.99);
                            func->SetParLimits(2, 0.01, 1);
                            func->SetParLimits(3, 0.01, 2);
                            slopePar = 2;
                            // fitrange = {0.8,8.8};
                        }
                        if (fitvar == kFitBoseEinstein)
                        {
                            func = fm->GetBoseEinstein(masspdg, 0.3, 20);
                            slopePar = 1;
                            // fitrange = {0.8,1.6};
                        }
                        if (fitvar == kFitFermiDirac)
                        {
                            func = fm->GetFermiDirac(masspdg, 0.3, 20);
                            slopePar = 1;
                            // fitrange = {0.8,1.6};
                        }
                        auto hfinal =
                            yieldmeannew::YieldMeanNew(hspectra_stat, hspectra_sys, hspectra_sys_nocor, hspectra_sys_cor, hspectra_sys_multicor, fout, func, 0.05, 10,
                                                       0.01, 0.1, false, Form("log_%d_%d.root", iAnti, iCent), "", fitoption.Data(), fitstart, fitend);
                        yield.push_back(hfinal->GetBinContent(1));
                        meanpt.push_back(hfinal->GetBinContent(5));
                        fout->SetLineWidth(1);
                        fout->SetLineStyle(flinestyle);
                        fout->SetRange(0, fitend);
                        outputFile->cd();
                        fout->Write(Form("%s_%s_%.2f-%.2f", multistring.Data(), functions[fitvar].Data(), fitstart, fitend));
                        fitchi2 = fout->GetChisquare();
                        fitNDF = fout->GetNDF();
                        log.Get(TLogLevel::INFO) << "Fit Result!!: " << fitchi2 / fitNDF << ", fitchi2: " << fitchi2 << ", fitNDF: " << fitNDF;
                        funcs.push_back(fout);
                        flinestyle++;
                    }
                }
            }
            int sysColorPallet = GetSerialColors(functions.size());
            // for memo, small
            TLatex *tm = new TLatex();
            tm->SetNDC();
            tm->SetTextSize(0.03);
            auto legendhead = new TLegend(.4, .65, .95, .9);
            legendhead->SetNColumns(2);
            legendhead->SetBorderSize(0);
            legendhead->SetFillStyle(0);
            double fitsyserrory;
            double fitsyserrorm;
            // QA
            if (!iCent)
                log.Get(TLogLevel::INFO) << "INEL or MB";
            else
                log.Get(TLogLevel::INFO) << "Mutli: " << kMultiplicityBins[iCent - 1] << " - " << kMultiplicityBins[iCent];
            legendhead->Clear();
            cResults->cd();
            hspectra_sys_nocor->GetXaxis()->SetRangeUser(0, 3.2);
            hspectra_sys_nocor->SetMaximum(5 *
                                           hspectra_sys_nocor->GetBinContent(2));
            hspectra_sys_nocor->SetMinimum(0.5 *
                                           hspectra_sys_nocor->GetBinContent(11));
            hspectra_sys_nocor->Draw("E");
            legendhead->AddEntry(hspectra_sys_nocor, "Data", "PLE");

            // Result?
            auto ymax = max_element(std::begin(yield), std::end(yield));
            auto ymin = min_element(std::begin(yield), std::end(yield));
            auto mmax = max_element(std::begin(meanpt), std::end(meanpt));
            auto mmin = min_element(std::begin(meanpt), std::end(meanpt));

            TH1F *systematicyields =
                new TH1F("systematic_yields", "", 100, *ymin * 0.85, *ymax * 1.15);
            systematicyields->GetXaxis()->SetTitle("<d#it{N}/d#it{y}>");
            systematicyields->GetYaxis()->SetTitle("Counts");
            TH1F *systematicmeanpTs =
                new TH1F("systematic_meanpTs", "", 100, *mmin * 0.85, *mmax * 1.15);
            systematicmeanpTs->GetXaxis()->SetTitle("<#it{p}_{T}>");
            systematicmeanpTs->GetYaxis()->SetTitle("Counts");

            double centralyield;
            double centralmeanpT;
            int loopcheck = 0;
            // loop, loop for the contents
            for (int fitvar = 1; fitvar < functions.size(); fitvar++)
            {
                log.Get(TLogLevel::INFO) << "fit funcion: " << functions[fitvar].Data();
                for (auto const &fitstart :
                     functionFitRangeVar[fitvar][0])
                { // fit bin start
                    for (auto const &fitend :
                         functionFitRangeVar[fitvar][1])
                    { // fit bin end
                        fitnames[loopcheck] =
                            Form("%s(dydn %.2f, mpt: %.2f)", fitnames[loopcheck].Data(),
                                 1e2 * yield[loopcheck], meanpt[loopcheck]);
                        // funcs[loopcheck]->SetLineColor(sysColorPallet + functions.size() -
                        //                                fitvar);
                        funcs[loopcheck]->SetLineColor(sysColorPallet + fitvar);
                        funcs[loopcheck]->Draw("same");
                        legendhead->AddEntry(funcs[loopcheck], fitnames[loopcheck], "L");

                        log.Get(TLogLevel::INFO) << "yields: " << yield[loopcheck] << " meanpT"
                             << meanpt[loopcheck]
                             << ", Chi^2: " << funcs[loopcheck]->GetChisquare()
                             << ", NDF: " << funcs[loopcheck]->GetNDF() << ", Chi^2/NDF: "
                             << funcs[loopcheck]->GetChisquare() /
                                    funcs[loopcheck]->GetNDF()
                            ;

                        if ((fitvar == kFitLevi) &&
                            (fitstart == centralvalue[centralvaluebin][0]) &&
                            (fitend == centralvalue[centralvaluebin][1]))
                        {
                            centralyield = yield[loopcheck];
                            centralmeanpT = meanpt[loopcheck];
                        }

                        if (fitvar != kFitBoseEinstein)
                            systematicyields->Fill(yield[loopcheck]);
                        if (fitvar != kFitBoseEinstein)
                            systematicmeanpTs->Fill(meanpt[loopcheck]);
                        loopcheck++;
                    }
                }
            }
            legendhead->Draw();
            SaveCanvas(cResults, Form("Spectrafit_%s_zoom_%d", kParticleType[iAnti].c_str(), iCent), Form("%s/LowpTExtrapolate/", kFiguresFolder.data()));
            outputFile->cd();
            cResults->Write(Form("Spectrafit_%s_zoom", multistring.Data()));
            systematicyields->Write(Form("Yieldsys_%s", multistring.Data()));
            systematicmeanpTs->Write(Form("MeanPtsys_%s", multistring.Data()));
            cResults->SetLogy(false);
            systematicyields->Draw();
            TLine *centralarrow1 = new TLine(centralyield, 0, centralyield, systematicyields->GetMaximum());
            centralarrow1->SetLineColor(kRed);
            centralarrow1->SetLineStyle(2);
            centralarrow1->Draw();
            SaveCanvas(cResults, Form("Yields_systematic_%s_zoom_%d", kParticleType[iAnti].c_str(), iCent),
                       Form("%s/LowpTExtrapolate/", kFiguresFolder.data()));
            log.Get(TLogLevel::INFO) << "mean: " << systematicyields->GetMean()
                 << ", stdev: " << systematicyields->GetRMS() << ", percentage: "
                 << systematicyields->GetRMS() / systematicyields->GetMean();
            systematicmeanpTs->Draw();
            TLine *centralarrow2 = new TLine(centralmeanpT, 0, centralmeanpT,
                                             systematicyields->GetMaximum());
            centralarrow2->SetLineColor(kRed);
            centralarrow2->SetLineStyle(2);
            centralarrow2->Draw();
            SaveCanvas(cResults,
                       Form("MeanPts_systematic_%s_%d",
                            kParticleType[iAnti].c_str(), iCent),
                       Form("%s/LowpTExtrapolate/", kFiguresFolder.data()));
            fitsyserrory = systematicyields->GetRMS();
            fitsyserrorm = systematicmeanpTs->GetRMS();

            // Systematic error
            func = fm->GetLevi(masspdg, 0.4, 750, 3);
            func->SetParLimits(1, 0.0001, 20000);
            slopePar = 2;
            auto hfinal =
                yieldmeannew::YieldMeanNew(hspectra_stat, hspectra_sys, hspectra_sys_nocor, hspectra_sys_cor, hspectra_sys_multicor, fout, func, 0.05, 10,
                                           0.01, 0.1, false, Form("log_%d_%d_%d.root", iAnti, iAnti, iCent), "", fitoption.Data(), centralvalue[centralvaluebin][0], centralvalue[centralvaluebin][1]);
            auto houty = new TH1D(Form("hYield_%s", multistring.Data()), "", 7, 0, 7);
            houty->SetBinContent(1, hfinal->GetBinContent(yieldmeannew::kYield));          // central value
            houty->SetBinContent(2, hfinal->GetBinContent(yieldmeannew::kYieldStat));      // stat.err.
            houty->SetBinContent(3, hfinal->GetBinContent(yieldmeannew::kYieldSysTot));    // tot syst. err.
            houty->SetBinContent(4, fitsyserrory);                                         // fit syst. err.
            houty->SetBinContent(5, hfinal->GetBinContent(yieldmeannew::kYieldSysHiCorr)); // extreme high
            houty->SetBinContent(6, hfinal->GetBinContent(yieldmeannew::kYieldSysLoCorr)); // extreme low
            houty->SetBinContent(7, hfinal->GetBinContent(yieldmeannew::kExtra));          // extrapolated yield
            outputFile->cd();
            houty->Write(Form("hYield_%s", multistring.Data()));
            log.Get(TLogLevel::INFO) << "dNdy: Central value: " << hfinal->GetBinContent(yieldmeannew::kYield) << ", stat.e: " << hfinal->GetBinContent(yieldmeannew::kYieldStat)
                 << ", fit sys err.: " << fitsyserrory << " total sys err.: " << std::sqrt(std::pow(fitsyserrory, 2) + std::pow(hfinal->GetBinContent(yieldmeannew::kYieldSysTot), 2))
                 << ", extreme high: " << hfinal->GetBinContent(yieldmeannew::kYieldSysHiCorr)
                 << ", extreme low: " << hfinal->GetBinContent(yieldmeannew::kYieldSysLoCorr);

            auto houtm = new TH1D(Form("hMeanPt_%s", multistring.Data()), "", 6, 0, 6);
            houtm->SetBinContent(1, hfinal->GetBinContent(yieldmeannew::kMean));            // central value
            houtm->SetBinContent(2, hfinal->GetBinContent(yieldmeannew::kMeanStat));        // stat.err.
            houtm->SetBinContent(3, hfinal->GetBinContent(yieldmeannew::kMeanSysTot));      // syst.err.
            houtm->SetBinContent(4, fitsyserrorm);                                          // fit syst. err.
            houtm->SetBinContent(5, hfinal->GetBinContent(yieldmeannew::kMeanSysHardCorr)); // extreme high
            houtm->SetBinContent(6, hfinal->GetBinContent(yieldmeannew::kMeanSysSoftCorr)); // extreme low
            outputFile->cd();
            houtm->Write(Form("hMeanPt_%s", multistring.Data()));

            log.Get(TLogLevel::INFO) << "Mean pT: Central value: " << hfinal->GetBinContent(yieldmeannew::kMean) << ", stat.e: " << hfinal->GetBinContent(yieldmeannew::kMeanStat)
                 << ", fit sys err.: " << fitsyserrorm << " total sys err.: " << std::sqrt(std::pow(fitsyserrorm, 2) + std::pow(hfinal->GetBinContent(yieldmeannew::kMeanSysTot), 2))
                 << ", extreme high: " << hfinal->GetBinContent(yieldmeannew::kMeanSysHardCorr)
                 << ", extreme low: " << hfinal->GetBinContent(yieldmeannew::kMeanSysSoftCorr);
        }
    }
}