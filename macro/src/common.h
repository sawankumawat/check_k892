
// The variables that can be chaged are here ****************************************************
const string kParticle = "kstar/";
const bool multipanel_plots = 0;
const bool save_plots = 1;
// const string kfoldername_temp = "kstarqa_id21631/hInvMass";
const string kfoldername_temp = "kstarqa/hInvMass";
// const string kfoldername_temp = "lf-kstar892analysis";
// const string kfoldername_temp = "lf-k892analysis";

const string kvariation = ""; // change the variation here
// const string kvariation = "_Combined_cut"; // change the variation here
// const string kvariation = "_TOF"; // change the variation here
// const string kvariation = "_TOFHit"; // change the variation here
// const string kvariation = "_TPC20"; // change the variation here
// const string kvariation = "_no_addition_cuts_mother"; // change the variation here
// const string kvariation = "__alltof"; // change the variation here
////********************************************************************************************

// define datasets here
#define DATASET_LHC220_pass7

const string kDataFilename_temp1 = "../data/" + kParticle;  // data file
const string kSignalOutput_temp = "../output/" + kParticle; // output folder

#ifdef DATASET_LHC220_pass7
const string kDataset_temp = "LHC22o_pass7/IR_study/";
const string kMCDataset = "../mc/LHC24f3b/";
#endif

#ifdef DATASET_LHC220_pass7
// const string kDataFilename_temp2 = "269284.root"; // Marta code for pp 13.6 TeV
// const string kDataFilename_temp2 = "283597.root"; // DCA 0.1 cm
// const string kDataFilename_temp2 = "283598.root"; //  Neural network corrections, no combined PID cut, cuts on mother particle
// const string kDataFilename_temp2 = "447406.root"; //
// const string kDataFilename_temp2 = "448806.root"; // LHC22_pass7_small dataset, INEL > 0
// const string kDataFilename_temp2 = "448807.root"; // LHC22_pass7_small dataset, INEL

//***************TPC crossed rows placed in wrong place*****************************
// const string kDataFilename_temp2 = "448490.root"; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "448489.root"; // LHC22_pass7_medium dataset, INEL

// const string kDataFilename_temp2 = "449695.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "449694.root"; // LHC23_pass4_thin_small dataset, INEL

// const string kDataFilename_temp2 = "451993.root"; // 	LHC24_pass1_minBias dataset, INEL > 0
// const string kDataFilename_temp2 = "451992.root"; // 	LHC24_pass1_minBias dataset, INEL

// const string kDataFilename_temp2 = "451003.root"; // 	LHC24an_pass1_skimmed_small dataset, INEL > 0
// const string kDataFilename_temp2 = "451002.root"; // 	LHC24an_pass1_skimmed_small dataset, INEL

//*****************************Corrected TPC crossed rows********************************
// const string kDataFilename_temp2 = "459845.root"; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "459908.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "460233.root"; // 	LHC24_pass1_minBias dataset, INEL > 0

//*************************IR study********************************
// const string kDataFilename_temp2 = "463114.root"; // 1-2 MHz
// const string kDataFilename_temp2 = "535069.root"; // 14 kHz
// const string kDataFilename_temp2 = "535545.root"; // 70 kHz
const string kDataFilename_temp2 = "535645.root"; // 135 kHz
// const string kDataFilename_temp2 = "535999.root"; // 330 kHz
// const string kDataFilename_temp2 = "536106.root"; // 650 kHz

#endif

// final dataset name
const string kDataset = kDataFilename_temp1 + kDataset_temp;
const string kSignalOutput = kSignalOutput_temp + kDataset_temp + kDataFilename_temp2.substr(0, kDataFilename_temp2.rfind("."));
const string kDataFilename = kDataset + kDataFilename_temp2;
const string kfoldername = kfoldername_temp + kvariation;
const string koutputfolder = kSignalOutput + "/" + kfoldername;

// Canvas dimensions
const int klowerpad = 2;
const int kupperpad = 2;
const int kcanvaswidth = 720 * klowerpad;
const int kcanvasheight = 720 * kupperpad;
const int kcanvasdivide[2] = {klowerpad, kupperpad};

float masspdg = 0.895;  // in GeV/c^2
float widthpdg = 0.047; // in 1 sigma GeV/c^2

////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
//                               OLD DATASETS                                     //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////

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
