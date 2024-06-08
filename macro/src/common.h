
// The variables that can be chaged are here ****************************************************
const string kParticle = "kstar/";
// const string kParticle = "glueball/";
const bool multipanel_plots = 0;
const bool save_plots = 1;
const string kfoldername_temp = "kstarqa";
const string kvariation = "_tpc2"; // change the variation here
////********************************************************************************************

// define datasets here
#define DATASET_LHC220_pass6_small
// #define DATASET_LHC220_pass4_small

const string kDataFilename_temp1 = "../data/" + kParticle;  // data file
const string kSignalOutput_temp = "../output/" + kParticle; // output folder

// define datasets name here
#ifdef DATASET_LHC220_pass6_small
const string kDataset_temp = "LHC220_pass6_small/";
// const string kMCDataset = "../mc/LHC24b1/";
const string kMCDataset = "../mc/LHC24b1b/";
#endif

// define datasets train output run number here
#ifdef DATASET_LHC220_pass6_small
// const string kDataFilename_temp2 = "208396.root"; // data file
// const string kDataFilename_temp2 = "210677.root"; // data file
// const string kDataFilename_temp2 = "211075.root"; // tpc cluster 120
const string kDataFilename_temp2 = "211557.root"; // global tracks w/o dca on
// const string kMCFilename_temp = "210563.root";    // only mc process (b1b)
// const string kMCFilename_temp = "211233.root";    // tpc cluster 120 + mass pi fix
const string kMCFilename_temp = "211346.root"; // split tracks
#endif

// final dataset name
const string kDataset = kDataFilename_temp1 + kDataset_temp;
const string kSignalOutput = kSignalOutput_temp + kDataset_temp + kDataFilename_temp2.substr(0, kDataFilename_temp2.rfind("."));
const string kDataFilename = kDataset + kDataFilename_temp2;
const string kMCFilename = kMCDataset + kMCFilename_temp;
const string kfoldername = kfoldername_temp + kvariation;
const string koutputfolder = kSignalOutput + "/" + kfoldername;

// Canvas dimensions
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

// old
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
// const string kDataset = "LHC220_pass6_small/201194";
