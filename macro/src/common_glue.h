// The variables that can be chaged are here ****************************************************
// #define KKchannel   // for Kaon Kaon channel
#define KsKschannel  //for Kshort Kshort channel

const string kParticle = "glueball/";
const float txtsize = 0.045;
const string koutputtype = "png"; // pdf, eps
////********************************************************************************************
#ifdef KsKschannel
const string kchannel = "KsKs_Channel";
// const string kfoldername_temp = "strangeness_tutorial";
const string kfoldername_temp = "higher-mass-resonances";

const string kvariation = ""; // change the variation here
// const string kvariation = "_all_tight"; // change the variation here
// const string kvariation = "_PID3"; // change the variation here
// const string kvariation = "_3sigmaKs"; // change the variation here
#endif
////********************************************************************************************
#ifdef KKchannel
const string kchannel = "KK_Channel";
const string kfoldername_temp = "kaonkaonAnalysisRun3";
const string kvariation = ""; // change the variation here
// const string kvariation = "_Deep_angle"; // change the variation here
#endif

////******some fixed variables******************************************************
TLatex *t2 = new TLatex();
////******some fixed variables******************************************************

// define datasets here
// #define DATASET_LHC220_pass6_small
#define DATASET_LHC22o_pass7_small
// #define DATASET_LHC220_pass4_small

const string kDataFilename_temp1 = "../data/" + kParticle;  // data file
const string kSignalOutput_temp = "../output/" + kParticle; // output folder

// define datasets name here
#ifdef DATASET_LHC220_pass6_small
const string kDataset_temp = "LHC220_pass6_small/";
// const string kMCDataset = "../mc/LHC24b1/";
const string kMCDataset = "../mc/LHC24b1b/";
#endif

#ifdef DATASET_LHC22o_pass7_small
const string kDataset_temp = "LHC22o_pass7_small/";
const string kMCDataset = "../mc/LHC24b1b/";
#endif

// define datasets train output run number here
#ifdef DATASET_LHC220_pass6_small
// const string kDataFilename_temp2 = "215554.root"; // data file
// const string kDataFilename_temp2 = "221157.root"; // data file (with rotational background)
// const string kDataFilename_temp2 = "222487.root"; // data file (with larger pT range upto 30 GeV/c)
const string kDataFilename_temp2 = "230281.root"; // Medium data file (with larger pT range upto 30 GeV/c)
// const string kDataFilename_temp2 = "245896.root"; // Medium data file (Included  mulitiplicity distribution, dE/dx tpc plots)
// const string kDataFilename_temp2 = "248997.root"; // Small data file (Increased mult dist x axis range upto 70k, corrected the different entries in the pi- and pi+ daughters)
const string kMCFilename_temp = "211346.root";    // MC file
#endif

#ifdef DATASET_LHC22o_pass7_small
// const string kDataFilename_temp2 = "260782.root"; // Full dataset (corrected mass correlation plot)
// const string kDataFilename_temp2 = "340845.root"; // Medium dataset (for check only)
// const string kDataFilename_temp2 = "341909.root"; // small dataset (for check only)
// const string kDataFilename_temp2 = "341913.root"; // Full dataset (2 run numbers missing)
// const string kDataFilename_temp2 = "356242.root"; // Medium dataset (with and without tighter cuts)
// const string kDataFilename_temp2 = "351470.root"; // LHC23 pass4 thin (largest dataset)
// const string kDataFilename_temp2 = "351471.root"; // Full pass 7 dataset
const string kDataFilename_temp2 = "359454.root"; // Medium dataset for angular cut (hard process removal)

const string kMCFilename_temp = "211346.root";    // MC file // temporary since injected MC is not available
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
