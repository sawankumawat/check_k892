#include <string>
#include <TLatex.h>

// #define KKchannel   // for Kaon Kaon channel
#define KsKschannel // for Kshort Kshort channel

const std::string kParticle = "glueball/";
const float txtsize = 0.045;
const std::string koutputtype = "png"; // pdf, eps
////********************************************************************************************
#ifdef KsKschannel
const std::string kchannel = "KsKs_Channel";
const std::string kfoldername_temp = "higher-mass-resonances";

// const std::string kvariation = ""; // No variation
// const std::string kvariation = "_lambda_7";
// const std::string kvariation = "_all_tight";
// const std::string kvariation = "_PID3";
// const std::string kvariation = "_3sigmaKs";
const std::string kvariation = "_CS_Frame";
#endif

TLatex *t2 = new TLatex();

// define datasets here
// #define DATASET_LHC220_pass6_small
#define DATASET_LHC22o_pass7_small
// #define DATASET_LHC220_pass4_small

const std::string kDataFilename_temp1 = "../data/" + kParticle;  // data file
const std::string kSignalOutput_temp = "../output/" + kParticle; // output folder

#ifdef DATASET_LHC22o_pass7_small
const std::string kDataset_temp = "LHC22o_pass7_small/";
#endif

#ifdef DATASET_LHC22o_pass7_small
// const std::string kDataFilename_temp2 = "433479.root"; // Full 2022 pass 7 (old)
// const std::string kDataFilename_temp2 = "435448.root"; // Full 2022 pass 7 (new)
// const std::string kDataFilename_temp2 = "435450.root"; // LHC23_pass4_thin
// const std::string kDataFilename_temp2 = "435449.root"; // LHC24_pass4_skimmed

//===========Strangeness derived data==========================
const std::string kDataFilename_temp2 = "504802.root"; // LHC23_pass4_thin

#endif

// final dataset name
const std::string kDataset = kDataFilename_temp1 + kDataset_temp;
const std::string kSignalOutput = kSignalOutput_temp + kDataset_temp + kDataFilename_temp2.substr(0, kDataFilename_temp2.rfind("."));
const std::string kDataFilename = kDataset + kDataFilename_temp2;
const std::string kfoldername = kfoldername_temp + kvariation;
const std::string koutputfolder = kSignalOutput + "/" + kfoldername;

// Canvas dimensions
const int klowerpad = 2;
const int kupperpad = 2;
const int kcanvaswidth = 600 * klowerpad;
const int kcanvasheight = 600 * kupperpad;
const int kcanvasdivide[2] = {klowerpad, kupperpad};

//******************Resonance PDG values********************
float f980Mass = 0.990;
float f980MassErr = 0.020;
float f980Width = 0.055;
float f1270Mass = 1.275;
float f1270MassErr = 8e-4;
float f1270Width = 0.187;
float f1270WidthErr = 0.0024;
float a1320Mass = 1.318;
float a1320MassErr = 6e-4;
float a1320Width = 0.105;
float a1320WidthErr = 0.0018;
float f1500Mass = 1.522;
float f1500Width = 0.108;
float f1525Mass = 1.517;
float f1525MassErr = 24e-4;
float f1525Width = 0.0844;
float f1525WidthErr = 0.0027;
float f1710Mass = 1.713;
float f1710MassErr = 8e-3;
float f1710Width = 0.150;
float f1710WidthErr = 0.011;
float f2020Mass = 1.982;
float f2020Width = 0.436;

//************************************OLD*************************************//
//************************************OLD*************************************//

// define datasets train output run number here
#ifdef DATASET_LHC220_pass6_small
// const std::string kDataFilename_temp2 = "215554.root"; // data file
// const std::string kDataFilename_temp2 = "221157.root"; // data file (with rotational background)
// const std::string kDataFilename_temp2 = "222487.root"; // data file (with larger pT range upto 30 GeV/c)
const std::string kDataFilename_temp2 = "230281.root"; // Medium data file (with larger pT range upto 30 GeV/c)
// const std::string kDataFilename_temp2 = "245896.root"; // Medium data file (Included  mulitiplicity distribution, dE/dx tpc plots)
// const std::string kDataFilename_temp2 = "248997.root"; // Small data file (Increased mult dist x axis range upto 70k, corrected the different entries in the pi- and pi+ daughters)
#endif

// #ifdef DATASET_LHC22o_pass7_small
// // const std::string kDataFilename_temp2 = "260782.root"; // Full dataset (corrected mass correlation plot)
// // const std::string kDataFilename_temp2 = "340845.root"; // Medium dataset (lambda rejection cut variation)
// // const std::string kDataFilename_temp2 = "341909.root"; // small dataset (for check only)
// // const std::string kDataFilename_temp2 = "341913.root"; // Full dataset (2 run numbers missing)
// // const std::string kDataFilename_temp2 = "356242.root"; // Medium dataset (with and without tighter cuts)
// // const std::string kDataFilename_temp2 = "351470.root"; // LHC23 pass4 thin (largest dataset)
// // const std::string kDataFilename_temp2 = "351471.root"; // Full pass 7 dataset
// // const std::string kDataFilename_temp2 = "359454.root"; // Medium dataset angular cut (wrong cut used)
// // const std::string kDataFilename_temp2 = "358932.root"; // Full train with systematics
// // const std::string kDataFilename_temp2 = "363021.root"; // small train (ks mass cuts and angular separation cuts)
// // const std::string kDataFilename_temp2 = "363594.root"; // medium train (ks mass cuts)
// // const std::string kDataFilename_temp2 = "362701.root"; // medium train (angular separation cuts)
// // const std::string kDataFilename_temp2 = "369624.root"; // full trian with Ks sigma = 3.7 MeV/c^2
// // const std::string kDataFilename_temp2 = "370825.root"; // LHC23_thin_small (competing cascade rejection cuts)

//****************systematics train*******************************
// const std::string kvariation = "_id24937"; // first four are same
// const std::string kvariation = "_id24938";
// const std::string kvariation = "_id24939";
// const std::string kvariation = "_id24940";

// const std::string kvariation = "_DCA0p04_id24940";
// const std::string kvariation = "_DCA0p06_id24940";
// const std::string kvariation = "_DCAv0dau0p3_id24938";
// const std::string kvariation = "_DCAv0dau1_id24938";
// const std::string kvariation = "_Ks_selection4_id24939";
// const std::string kvariation = "_Ks_selection5_id24939";
// const std::string kvariation = "_TPCPID3_id24937";
// const std::string kvariation = "_TPCPID4_id24937";
// const std::string kvariation = "_TPCPID6_id24937";
// const std::string kvariation = "_TPCcr100_id24937";
// const std::string kvariation = "_TPCcr120_id24937";
// const std::string kvariation = "_TPCcrfc0p9_id24940";
// const std::string kvariation = "_TPCcrfc1p0_id24940";
// const std::string kvariation = "_cospa0p95_id24938";
// const std::string kvariation = "_cospa0p99_id24938";
// const std::string kvariation = "_decay_rad0p4_id24938";
// const std::string kvariation = "_decay_rad0p6_id24938";
// const std::string kvariation = "_lambda_rej4_id24939";
// const std::string kvariation = "_lambda_rej6_id24939";
// const std::string kvariation = "_lifetime15_id24939";
// const std::string kvariation = "_lifetime25_id24939";

//************* angular and Ks cut medium train
// const std::string kvariation = "_id24794";
// const std::string kvariation = "_id25081";
// const std::string kvariation = "_1Kscut_id24794";
// const std::string kvariation = "_1p5Kscut_id24794";
// const std::string kvariation = "_2Kscut_id24794";
// const std::string kvariation = "_4Kscut_id24794";
// const std::string kvariation = "_angsep_0p5_id25081";
// const std::string kvariation = "_angsep_1_id25081";
// const std::string kvariation = "_angsep_1p5_id25081";
// const std::string kvariation = "_angsep_2_id25081";
// const std::string kvariation = "_angsep_3_id25081";

// ////********************************************************************************************
// #ifdef KKchannel
// const std::string kchannel = "KK_Channel";
// const std::string kfoldername_temp = "kaonkaonAnalysisRun3";
// const std::string kvariation = ""; // change the variation here
// // const std::string kvariation = "_Deep_angle"; // change the variation here
// #endif
