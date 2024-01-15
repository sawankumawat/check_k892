#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>

#define DATASET_LHC22CDE_PASS4
// #define DATASET_LHC22M_PASS4
// #define DATASET_LHC22O_PASS4
// #define DATASET_LHC23ZX_PART
const std::string kRootExtension = ".root";

// Variables
// const bool kIsINEL = false;
const bool kIsINEL = true;
const std::string kBaseInputDir = "../";
#ifdef DATASET_LHC22CDE_PASS4
const std::string kBaseOutputDirTemp = "../output/pbpb";
#endif
#ifdef DATASET_LHC22M_PASS4
const std::string kBaseOutputDirTemp = "../output/LHC22m_pass4";
#endif
#ifdef DATASET_LHC22O_PASS4
const std::string kBaseOutputDirTemp = "../output/LHC22o_pass4";
#endif
#ifdef DATASET_LHC23ZX_PART
// const std::string kBaseOutputDirTemp = "../output/LHC22s_pass5";
const std::string kBaseOutputDirTemp = "../output/LHC23zx_part2";
#endif
const std::string kBaseOutputDir = (kIsINEL) ? kBaseOutputDirTemp + "_INEL/" : kBaseOutputDirTemp + "/";
const std::string kFiguresFolder = kBaseOutputDir + "figures/";

const std::string kOutputName = "lf-k892analysis";z

#ifdef DATASET_LHC22CDE_PASS4
// const std::string kDataFilename = kBaseInputDir + "data/LHC22cde_pass4/107888_89_90.root"; // without INELgt0
// const std::string kDataFilename = kBaseInputDir + "data/LHC22cde_pass4/109817_19_21.root"; // with INELgt0 TPCNN
// const std::string kDataFilename = kBaseInputDir + "data/LHC22cde_pass4/109818_20_22.root"; // with INELgt0 
// const std::string kDataFilename = kBaseInputDir + "data/LHC22cde_pass4/118014.root"; // with INELgt0 LHC22e
// const std::string kDataFilename = kBaseInputDir + "data/LHC22cde_pass4/118015_16.root"; // with INELgt0 LHC22cd
// const std::string kDataFilename = kBaseInputDir + "data/LHC22cde_pass4/118014_15_16.root"; // with INELgt0 LHC22cde
// const std::string kDataFilename = kBaseInputDir + "data/LHC22cde_pass4/128012_13_17.root"; // with right TPC-TOF Cut
const std::string kDataFilename = kBaseInputDir + "data/pbpb/AnalysisResults.root"; // pb-pb data


// const std::string kDataFilename = kBaseInputDir + "data/LHC22cde_pass4/LHC22e_partial.root"; // with INELgt0 LHC22e part

// const std::string kMCfilename = kBaseInputDir + "mc/LHC22h1b2/AnalysisResults_gen900.root"; // LHC22h1b2
// const std::string kMCfilename = kBaseInputDir + "mc/LHC22h1b11/AnalysisResults_gen900ideal.root"; // LHC22h1b11
// const std::string kMCfilename = kBaseInputDir + "mc/LHC22h1c1/117133.root"; // LHC22h1c1
// const std::string kMCfilename = kBaseInputDir + "mc/LHC22h1c1/evtCor.root"; // LHC22h1c1-LHC22e (evtCorrection)
// const std::string kMCfilename = kBaseInputDir + "mc/LHC22h1c1/LHC22cd_evtcor.root"; // LHC22h1c1-LHC22cd (evtCorrection)
// const std::string kMCfilename = kBaseInputDir + "mc/LHC22h1c1/LHC22e_partial.root"; // LHC22h1c1-LHC22cd (evtCorrection)
const std::string kMCfilename = kBaseInputDir + "mc/LHC22h1d1/128015.root"; // LHC22h1d1

#endif
#ifdef DATASET_LHC22M_PASS4
const std::string kDataFilename = kBaseInputDir + "data/LHC22m_pass4_small/107887.root"; // test
// const std::string kMCfilename = kBaseInputDir + "mc/LHC23f3/AnalysisResults_pp13.6TeV.root"; // test
// const std::string kMCfilename = kBaseInputDir + "mc/LHC23d1/AnalysisResults_LHC23d1j.root"; // test
const std::string kMCfilename = kBaseInputDir + "mc/LHC23d1k/117134.root"; // test
#endif
#ifdef DATASET_LHC22O_PASS4
// const std::string kDataFilename = kBaseInputDir + "data/LHC22o_pass4/2runs.root"; // test
const std::string kDataFilename = kBaseInputDir + "data/LHC22o_pass4/128014.root"; // test
// const std::string kMCfilename = kBaseInputDir + "mc/LHC23f4b2/117136.root"; // 
const std::string kMCfilename = kBaseInputDir + "mc/LHC23d1k/128016.root"; // test
#endif
#ifdef DATASET_LHC23ZX_PART
// const std::string kDataFilename = kBaseInputDir + "data/LHC23zx/AnalysisResults_K892.root"; // test
// const std::string kDataFilename = kBaseInputDir + "data/LHC22s_pass5/124251.root"; // test
// const std::string kDataFilename = kBaseInputDir + "data/LHC23zx/124625.root"; // test
const std::string kDataFilename = kBaseInputDir + "data/LHC23zx/AnalysisResults_K892_3.root"; // test
const std::string kMCfilename = kBaseInputDir + "mc/LHC23f4b2/117136.root"; //
#endif

const std::string kInitOutput = kBaseOutputDir + "common/init" + kRootExtension;
const std::string kSignalOutput = kBaseOutputDir + "common/signal" + kRootExtension;
const std::string kEfficiencyOutput = kBaseOutputDir + "common/efficiency" + kRootExtension;
const std::string kSpectraOutput = kBaseOutputDir + "common/spectra" + kRootExtension;
const std::string kSystematicsOutput = kBaseOutputDir + "common/systematics" + kRootExtension;
const std::string kYieldMeanOutput = kBaseOutputDir + "common/yieldmean" + kRootExtension;


const std::string kInitOutputPart = kInitOutput.substr(0, kInitOutput.size() - kRootExtension.size());
const std::string kSignalOutputPart = kSignalOutput.substr(0, kSignalOutput.size() - kRootExtension.size());
const std::string kEfficiencyOutputPart = kEfficiencyOutput.substr(0, kEfficiencyOutput.size() - kRootExtension.size());
const std::string kSpectraOutputPart = kSpectraOutput.substr(0, kSpectraOutput.size() - kRootExtension.size());

const string kYieldMeanOutputPart = kBaseOutputDir + "yieldmean";
const string kResultOutput = kBaseOutputDir + "results" + kRootExtension;
const string kResultOutputPart = kBaseOutputDir + "results";
const string kFinalOutput = kBaseOutputDir + "final" + kRootExtension;

// systematic details
const std::string kSystematicsDetailOutput = kBaseOutputDir + "systematic_details/systematics";

const std::string kFilterListResoTaskNames = "lf-reso2initializer";
const std::string kFilterListNames = "lf-k892analysis";
const std::string kFilterListNamesMC = "lf-k892analysis";

// Some constants
const float masspdg = 0.895;  // in GeV/c^2
const float widthpdg = 0.047; // in 1 sigma GeV/c^2
const float branchingratio = 0.66;
const float kMaterialBudgetError = 0.05; // FIXME: need to be updated

// Analysis
const bool kDoMCQA = false;
const std::vector<string> kParticleType = {"K892", "AntiK892"};

const int kAntiBinSize = 1; // if 2, analyse anti particle separately
const int kRebin = 4;
const float kINELnorm = 0.61;
#ifdef DATASET_LHC22CDE_PASS4
const std::vector<double> kMultiplicityBins = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
#endif
#ifdef DATASET_LHC22M_PASS4
const std::vector<double> kMultiplicityBins = {0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
#endif
#ifdef DATASET_LHC22O_PASS4
const std::vector<double> kMultiplicityBins = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
#endif
#ifdef DATASET_LHC23ZX_PART
const std::vector<double> kMultiplicityBins = {0, 100};
#endif

const std::vector<Double_t> kDrawRange = {0.65, 1.4};
const std::vector<Double_t> kDrawFitRange = {0.7, 1.2};
const std::vector<Double_t> kPeakRange = {0.84, 0.95};
// const std::vector<Double_t> kpTbin = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
//                                 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
//                                 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.5, 5.0, 6.0, 8.0};
#ifdef DATASET_LHC22CDE_PASS4
const std::vector<Double_t> kpTbin = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 5.0};
// const std::vector<Double_t> kpTbin = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0};
#endif
#ifdef DATASET_LHC22M_PASS4
const std::vector<Double_t> kpTbin = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};
#endif
#ifdef DATASET_LHC22O_PASS4
const std::vector<Double_t> kpTbin = {1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0};
#endif
#ifdef DATASET_LHC23ZX_PART
const std::vector<Double_t> kpTbin = {1.0, 10.0};
// const std::vector<Double_t> kpTbin = {1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0};
#endif
const double kLastpTbin = kpTbin[kpTbin.size() - 1];

const int kNPtBins = kpTbin.size();
const int kNMultiplicityBins = (kIsINEL) ? 1 : kMultiplicityBins.size() - 1;

// Normalization
//below the vector consists of vectors. to access it we need to use kNormRangepT[x][y]
#ifdef DATASET_LHC22CDE_PASS4
vector<vector<double>> kNormRangepT = { // 900 GeV
    {1.04, 1.07},  // 0.0-0.2
    {0.72, 0.75},  // 0.2-0.4
    {0.72, 0.75},  // 0.4-0.6
    {0.72, 0.75},  // 0.6-0.8
    {1.05, 1.08},  // 0.8-1.0
    {1.05, 1.08},  // 1.0-1.2
    {1.1, 1.2},  // 1.2-1.4
    {1.1, 1.2},  // 1.4-1.6
    {1.2, 1.3},  // 1.6-1.8
    {1.2, 1.3},  // 1.8-2.0
    {1.2, 1.3},  // 2.0-2.5
    {1.2, 1.3},  // 2.5-3.0
    {1.2, 1.3}   // 3.0-5.0
};
vector<vector<double>> kFitRange = { // 900 GeV
    {0.72, 1.1},  // 0.0-0.2
    {0.75, 1.04},  // 0.2-0.4
    {0.75, 1.1},  // 0.4-0.6
    {0.75, 1.1},  // 0.6-0.8
    {0.75, 1.04},  // 0.8-1.0
    {0.75, 1.04},  // 1.0-1.2
    {0.75, 1.04},  // 1.2-1.4
    {0.78, 1.04},  // 1.4-1.6
    {0.75, 1.04},  // 1.6-1.8
    {0.78, 1.1},  // 1.8-2.0
    {0.75, 1.1},  // 2.0-2.5
    {0.75, 1.1},  // 2.5-3.0
    {0.75, 1.1}   // 3.0-5.0
};
vector<vector<double>> kFitPol3ParLimit = { // 900 GeV
    {100, 200},  // 0.0-0.2
    {100, 200},  // 0.2-0.4
    {500, 1000},  // 0.4-0.6
    {100, 200},  // 0.6-0.8
    {100, 200},  // 0.8-1.0
    {100, 200},  // 1.0-1.2
    {100, 200},  // 1.2-1.4
    {100, 200},  // 1.4-1.6
    {100, 200},  // 1.6-1.8
    {50, 2000},  // 1.8-2.0
    {50, 2000},  // 2.0-2.5
    {100, 200},  // 2.5-3.0
    {100, 200}   // 3.0-5.0
};
#endif
#ifdef DATASET_LHC22M_PASS4
vector<vector<double>> kNormRangepT = { // 13.6 TeV
    {0.74, 0.75},  // 0.0-0.2
    {0.75, 0.78},  // 0.2-0.4
    {0.74, 0.75},  // 0.4-0.6
    {0.74, 0.75},  // 0.6-0.8
    {1.05, 1.06},  // 0.8-1.0
    {1.05, 1.06},  // 1.0-1.2
    {1.1, 1.2},  // 1.2-1.4
    {1.1, 1.2},  // 1.4-1.6
    {1.2, 1.3},  // 1.6-1.8
    {1.2, 1.3},  // 1.8-2.0
    {1.2, 1.3},  // 2.0-2.4
    {1.2, 1.3},  // 2.4-2.8
    {1.2, 1.3},  // 2.8-3.2
    {1.2, 1.3},  // 3.2-3.6
    {1.2, 1.3},  // 3.6-4.0
    {1.2, 1.3},  // 4.0-5.0
    {1.2, 1.3},  // 5.0-6.0
    {1.2, 1.3},  // 6.0-7.0
    {1.2, 1.3},  // 7.0-8.0
    {1.2, 1.3},  // 8.0-10.0
};
vector<vector<double>> kFitRange = { // 13.6 TeV
    {0.70, 1.15},  // 0.0-0.2
    {0.72, 1.05},  // 0.2-0.4
    {0.74, 1.15},  // 0.4-0.6
    {0.70, 1.15},  // 0.6-0.8
    {0.74, 1.15},  // 0.8-1.0
    {0.74, 1.12},  // 1.0-1.2
    {0.74, 1.12},  // 1.2-1.4
    {0.73, 1.12},  // 1.4-1.6
    {0.74, 1.12},  // 1.6-1.8
    {0.74, 1.12},  // 1.8-2.0
    {0.74, 1.12},  // 2.0-2.4
    {0.74, 1.12},  // 2.4-2.8
    {0.74, 1.12},  // 2.8-3.2
    {0.74, 1.12},  // 3.2-3.6
    {0.74, 1.12},  // 3.6-4.0
    {0.74, 1.12},  // 4.0-5.0
    {0.74, 1.12},  // 5.0-6.0
    {0.74, 1.12},  // 6.0-7.0
    {0.74, 1.12},  // 7.0-8.0
    {0.74, 1.12},  // 8.0-10.0
};
vector<vector<double>> kFitPol3ParLimit = { // 900 GeV
    {100, 200},  // 0.0-0.2
    {100, 200},  // 0.2-0.4
    {500, 1000},  // 0.4-0.6
    {100, 200},  // 0.6-0.8
    {100, 200},  // 0.8-1.0
    {100, 200},  // 1.0-1.2
    {100, 200},  // 1.2-1.4
    {100, 200},  // 1.4-1.6
    {100, 200},  // 1.6-1.8
    {50, 2000},  // 1.8-2.0
    {100, 200},  // 2.0-2.4
    {100, 200},  // 2.4-2.8
    {100, 200},  // 2.8-3.2
    {100, 200},  // 3.2-3.6
    {100, 200},  // 3.6-4.0
    {100, 200},  // 4.0-5.0
    {100, 200},  // 5.0-6.0
    {100, 200},  // 6.0-7.0
    {100, 200},  // 7.0-8.0
    {100, 200},  // 8.0-10.0
};
#endif
#ifdef DATASET_LHC23ZX_PART
vector<vector<double>> kNormRangepT = { // 5.6 TeV
    // {0.96, 1.0},  // 0.0-0.2
    // {0.75, 0.78},  // 0.2-0.4
    // {0.74, 0.75},  // 0.4-0.6
    // {0.95, 1.05},  // 0.6-0.8
    {0.7, 0.8},  // 0.8-1.0
    {1.05, 1.06},  // 0.8-1.0
    {1.05, 1.06},  // 1.0-1.2
    {1.1, 1.2},  // 1.2-1.4
    {1.1, 1.2},  // 1.4-1.6
    {1.2, 1.3},  // 1.6-1.8
    {1.2, 1.3},  // 1.8-2.0
    {1.2, 1.3},  // 2.0-2.4
    {1.2, 1.3},  // 2.4-2.8
    {1.2, 1.3},  // 2.8-3.2
    {1.2, 1.3},  // 3.2-3.6
    {1.2, 1.3},  // 3.6-4.0
    {1.2, 1.3},  // 4.0-5.0
    {1.2, 1.3},  // 5.0-6.0
    {1.2, 1.3},  // 6.0-7.0
    {1.2, 1.3},  // 7.0-8.0
    {1.2, 1.3},  // 8.0-10.0
    {1.2, 1.3},  // 10.0-15.0
};
vector<vector<double>> kFitRange = { // 13.6 TeV
    // {0.79, 1.05},  // 0.0-0.2
    // {0.79, 1.05},  // 0.2-0.4
    // {0.80, 1.05},  // 0.4-0.6
    // {0.75, 1.05},  // 0.6-0.8
    {0.75, 1.05},  // 0.8-1.0
    {0.78, 1.05},  // 1.0-1.2
    {0.78, 1.05},  // 1.2-1.4
    {0.74, 1.05},  // 1.4-1.6
    {0.79, 1.1},  // 1.6-1.8
    {0.79, 1.1},  // 1.8-2.0
    {0.78, 1.1},  // 2.0-2.4
    {0.78, 1.1},  // 2.4-2.8
    {0.78, 1.1},  // 2.8-3.2
    {0.78, 1.1},  // 3.2-3.6
    {0.78, 1.1},  // 3.6-4.0
    {0.78, 1.1},  // 4.0-5.0
    {0.78, 1.1},  // 5.0-6.0
    {0.82, 1.1},  // 6.0-7.0
    {0.82, 1.1},  // 7.0-8.0
    {0.82, 1.1},  // 8.0-10.0
    {0.82, 1.1},  // 10.0-15.0
};
vector<vector<double>> kFitPol3ParLimit = { // 13.6 TeV
    // {100, 200},  // 0.0-0.2
    // {100, 200},  // 0.2-0.4
    // {500, 1000},  // 0.4-0.6
    // {100, 200},  // 0.6-0.8
    {100, 200},  // 0.8-1.0
    {100, 200},  // 1.0-1.2
    {100, 200},  // 1.2-1.4
    {100, 200},  // 1.4-1.6
    {100, 200},  // 1.6-1.8
    {50, 2000},  // 1.8-2.0
    {100, 200},  // 2.0-2.4
    {100, 200},  // 2.4-2.8
    {100, 200},  // 2.8-3.2
    {100, 200},  // 3.2-3.6
    {100, 200},  // 3.6-4.0
    {100, 200},  // 4.0-5.0
    {100, 200},  // 5.0-6.0
    {100, 200},  // 6.0-7.0
    {100, 200},  // 7.0-8.0
    {100, 200},  // 8.0-10.0
    {100, 200},  // 10.0-15.0
};
#endif
#ifdef DATASET_LHC22M_PASS4
vector<vector<double>> kNormRangepT = { // 13.6 TeV
    {0.74, 0.75},  // 0.0-0.2
    {0.75, 0.78},  // 0.2-0.4
    {0.74, 0.75},  // 0.4-0.6
    {0.74, 0.75},  // 0.6-0.8
    {1.05, 1.06},  // 0.8-1.0
    {1.05, 1.06},  // 1.0-1.2
    {1.1, 1.2},  // 1.2-1.4
    {1.1, 1.2},  // 1.4-1.6
    {1.2, 1.3},  // 1.6-1.8
    {1.2, 1.3},  // 1.8-2.0
    {1.2, 1.3},  // 2.0-2.4
    {1.2, 1.3},  // 2.4-2.8
    {1.2, 1.3},  // 2.8-3.2
    {1.2, 1.3},  // 3.2-3.6
    {1.2, 1.3},  // 3.6-4.0
    {1.2, 1.3},  // 4.0-5.0
    {1.2, 1.3},  // 5.0-6.0
    {1.2, 1.3},  // 6.0-7.0
    {1.2, 1.3},  // 7.0-8.0
    {1.2, 1.3},  // 8.0-10.0
};
vector<vector<double>> kFitRange = { // 13.6 TeV
    {0.70, 1.15},  // 0.0-0.2
    {0.72, 1.05},  // 0.2-0.4
    {0.74, 1.15},  // 0.4-0.6
    {0.70, 1.15},  // 0.6-0.8
    {0.74, 1.15},  // 0.8-1.0
    {0.74, 1.12},  // 1.0-1.2
    {0.74, 1.12},  // 1.2-1.4
    {0.73, 1.12},  // 1.4-1.6
    {0.74, 1.12},  // 1.6-1.8
    {0.74, 1.12},  // 1.8-2.0
    {0.74, 1.12},  // 2.0-2.4
    {0.74, 1.12},  // 2.4-2.8
    {0.74, 1.12},  // 2.8-3.2
    {0.74, 1.12},  // 3.2-3.6
    {0.74, 1.12},  // 3.6-4.0
    {0.74, 1.12},  // 4.0-5.0
    {0.74, 1.12},  // 5.0-6.0
    {0.74, 1.12},  // 6.0-7.0
    {0.74, 1.12},  // 7.0-8.0
    {0.74, 1.12},  // 8.0-10.0
};
vector<vector<double>> kFitPol3ParLimit = { // 900 GeV
    {100, 200},  // 0.0-0.2
    {100, 200},  // 0.2-0.4
    {500, 1000},  // 0.4-0.6
    {100, 200},  // 0.6-0.8
    {100, 200},  // 0.8-1.0
    {100, 200},  // 1.0-1.2
    {100, 200},  // 1.2-1.4
    {100, 200},  // 1.4-1.6
    {100, 200},  // 1.6-1.8
    {50, 2000},  // 1.8-2.0
    {100, 200},  // 2.0-2.4
    {100, 200},  // 2.4-2.8
    {100, 200},  // 2.8-3.2
    {100, 200},  // 3.2-3.6
    {100, 200},  // 3.6-4.0
    {100, 200},  // 4.0-5.0
    {100, 200},  // 5.0-6.0
    {100, 200},  // 6.0-7.0
    {100, 200},  // 7.0-8.0
    {100, 200},  // 8.0-10.0
};
#endif

// Systematic configuration
const std::vector<string> SystematicBins = {
    "", // Dump
    "", "_PVTrk", "_DCAxy_loose", "_DCAz_loose", "_DCA_loose", "_PID25", "_PID35", "_PID20", "_PID40", "_Trk_loose"};
    // "_id5195", "_id5196", "_id5197", "_id5198", "_id5199", "_id5200", "_id5201", "_id5202", "_id5203", "_id5204"};
const std::vector<string> SystematicBinsName = {
    "", // Dump
    "", "_PVTrk", "_DCAxy_loose", "_DCAz_loose", "_DCA_loose", "_PID25", "_PID35", "_PID20", "_PID40", "_Trk_loose"};

enum {
    kECbegin = 0,
    kINEL = 1,
    kINEL10,
    kINELg0,
    kINELg010,
    kTrig,
    kINELg0Trig,
    kINELg010Trig,
    kECend,
};

// Plotter
const int kCanvasW = 720;
const int kCanvasH = 720;
const std::string kOutputType = "png";

// Common functions
bool IsFolderExist(TString path) {
    if (gSystem->Exec(Form("test -d %s", path.Data())) != 0) {
        return false;
    } else {
        return true;
    }
}
int MakeFolder(TString path) {
    if (IsFolderExist(path)) {
        return 0;
    } else {
        return gSystem->Exec(Form("mkdir -p %s", path.Data()));
    }
}

#endif
