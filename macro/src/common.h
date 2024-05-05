const string kParticle = "kstar/";
// const string kParticle = "phi/";
// const string kParticle = "glueball/";

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
//  const string kDataset = "LHC23zzg_apass2";
//  const string kDataset = "LHC23zzf_pass2_QC";
//  const string kDataset = "LHC23zzh_pass2_small";

// pp datasets
//  const string kDataset = "LHC23_pass1_lowB_highIR_sampling";
//  const string kDataset = "LHC23_pass1_lowB_lowIR";
//  const string kDataset = "LHC23_pass1_QC1_sampling";
//  const string kDataset = "LHC23_pass1_QC2_sampling";
// const string kDataset = "LHC23zzs"; // high IR (~650k)
// const string kDataset = "LHC23zb";  // high IR (~1000K)
// const string kDataset = "LHC23zf";  // low IR (~10K)
// const string kDataset = "LHC23zk";  // low IR (~10K)
// const string kDataset = "LHC23zm";  // IR (~50K)
// const string kDataset = "LHC23f";  // IR (~10K)
// const string kDataset = "LHC23r";  // IR (~330K)
// const string kDataset = "LHC23f_pass1"
// const string kDataset = "LHC23h";  // IR (~130 kHz)
// const string kDataset = "LHC23t";  // IR (1.3 MHz)
//  const string kDataset = "LHC22o_apass4";
//  const string kDataset = "900GeV";
// const string kDataset = "LHC220_pass6_small/188648";
// const string kDataset = "LHC220_pass6_small/190048";
// const string kDataset = "LHC220_pass6_small/190185";
// const string kDataset = "LHC220_pass6_small/190325";
// const string kDataset = "LHC220_pass6_small/197586";
// const string kDataset = "LHC220_pass6_small/200200";
const string kDataset = "LHC220_pass6_small/201194";

const string kDataFilename = "../data/" + kcoll + kParticle + kDataset + "/AnalysisResults.root"; // data file
const string kSignalOutput = "../output/" + kcoll + kParticle + kDataset;                         // output directory
// const string kSignalOutput = "/home/sawan/check_k892/mc";                       // output directory for mc

const int klowerpad = 2;
const int kupperpad = 2;
const int kcanvaswidth = 600 * klowerpad;
const int kcanvasheight = 600 * kupperpad;
const int kcanvasdivide[2] = {klowerpad, kupperpad};

float masspdg = 0.895;  // in GeV/c^2
float widthpdg = 0.047; // in 1 sigma GeV/c^2

float f1270Mass = 1.275;
float f1270Width = 0.187;
float a1320Mass = 1.318;
float a1320Width = 0.105;
float f1525Mass = 1.518;
float f1525Width = 0.082;
float f1710Mass = 1.710;
float f1710Width = 0.150;
