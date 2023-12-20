const string kBaseInputDir = "../";
const string kBaseOutputDir = "../output/LHC22q/";
const string kFiguresFolder = kBaseOutputDir + "figures/";
const string kFiguresFolder2 = kBaseOutputDir + "figures3/";

const string kOutputName = "lf-k892analysis";

const string kDataFilename = kBaseInputDir + "data/LHC22m_pass4_temp/AnalysisResults.root";  // test
// const string kDataFilename = kBaseInputDir + "data/LHC22q_pass3/AnalysisResults.root"; // test
const string kSignalOutput = kBaseOutputDir + "common/signal.root";

// Analysis
const int kRebin = 5;
std::vector<Double_t> kDrawRange = {0.7, 1.4};
std::vector<Double_t> kDrawFitRange = {0.7, 1.2};
std::vector<Double_t> kpTbin = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                                1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
                                3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.5, 5.0, 6.0, 8.0};
// Plotter
const int kCanvasW = 1280;
const int kCanvasH = 720;