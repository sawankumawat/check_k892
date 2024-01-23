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
//  const string kDataset = "LHC23zzf_pass2_QC";  //latest

// pp datasets
//  const string kDataset = "LHC23_pass1_lowB_highIR_sampling";
//  const string kDataset = "LHC23_pass1_lowB_lowIR";
//  const string kDataset = "LHC23_pass1_QC1_sampling";
//  const string kDataset = "LHC23_pass1_QC2_sampling";
//  const string kDataset = "LHC23zzs"; // high IR (~650k)
// const string kDataset = "LHC23zb";  // high IR (~1000K)
// const string kDataset = "LHC23zf";  // low IR (~10K)
// const string kDataset = "LHC23zk";  // low IR (~10K)
// const string kDataset = "LHC23zm";  // IR (~50K)
// const string kDataset = "LHC23f";  // IR (~10K)
// const string kDataset = "LHC23r";  // IR (~330K)
// const string kDataset = "LHC23f_pass1"
// const string kDataset = "LHC23h";  // IR (~130 kHz)
const string kDataset = "LHC23t";  // IR (1.3 MHz)



const string kDataFilename = "../data/" + kcoll + kParticle + kDataset + "/AnalysisResults.root"; // data file
const string kSignalOutput = "../output/" + kcoll + kParticle + kDataset;                         // output directory


const int kcanvaswidth = 1400 * 1.1;
const int kcanvasheight = 1000 * 1.1;
const int klowerpad = 4;
const int kupperpad = 3;
const int kcanvasdivide[2] = {klowerpad, kupperpad};

float masspdg = 0.895;  // in GeV/c^2
float widthpdg = 0.047; // in 1 sigma GeV/c^2
