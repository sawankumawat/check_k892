const string kParticle = "kstar/";
// const string kParticle = "phi/";

// const string kcoll = "pbpb/";
const string kcoll = "pp/";

const bool multipanel_plots = 0;
const bool save_plots = 1;

// pbpb datasets
//  const string kDataset = "23zzk_pass1_relval";
//  const string kDataset = "LHC23zzh_cpass8";
//  const string kDataset = "pass1_golden_runs_QC_sampling";
//  const string kDataset = "LHC23zzh_pass1_small";
//  const string kDataset = "LHC23zzh_pass1";

// pp datasets
//  const string kDataset = "LHC23_pass1_lowB_highIR_sampling";
//  const string kDataset = "LHC23_pass1_lowB_lowIR";
//  const string kDataset = "LHC23_pass1_QC1_sampling";
//  const string kDataset = "LHC23_pass1_QC2_sampling";
//  const string kDataset = "LHC23zzs";
const string kDataset = "LHC23zf";

// const string kfoldername = (kParticle == "kstar/") ? "lf-k892analysis" : "phianalysisrun3_2_nsigma_cut";

const string kDataFilename = "../data/" + kcoll + kParticle + kDataset + "/AnalysisResults.root"; // data file
const string kSignalOutput = "../output/" + kcoll + kParticle + kDataset;                         // output directory

// // Analysis
// std::vector<Double_t> kpTbin = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
//                                 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
//                                 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.5, 5.0, 6.0, 8.0};

const int kcanvaswidth = 1400 * 1.1;
const int kcanvasheight = 1000 * 1.1;
const int klowerpad = 4;
const int kupperpad = 3;
const int kcanvasdivide[2] = {klowerpad, kupperpad};

const Int_t Npt = 21; // 12 for pbpb and 18 and 21 for pp
const int pt_start = 0;
const int pt_end = Npt;

float masspdg = 0.895;  // in GeV/c^2
float widthpdg = 0.047; // in 1 sigma GeV/c^2
