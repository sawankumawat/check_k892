// The variables that can be chaged are here ****************************************************
// #define KKchannel   // for Kaon Kaon channel
#define KsKschannel  //for Kshort Kshort channel

const string kParticle = "glueball/";
const float txtsize = 0.045;
const string koutputtype = "png"; // pdf, eps
////********************************************************************************************
#ifdef KsKschannel
const string kchannel = "KsKs_Channel";
const string kfoldername_temp = "strangeness_tutorial";
const string kvariation = ""; // change the variation here
// const string kvariation = "_full_ks_distribution"; // change the variation here
// const string kvariation = "_old_cuts"; // change the variation here
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
// const string kDataFilename_temp2 = "240457.root"; // data file (small)
// const string kDataFilename_temp2 = "240361.root"; // data file (medium)
// const string kDataFilename_temp2 = "244178.root"; // data file (medium)
// const string kDataFilename_temp2 = "247473.root"; // Medium data file medium (fixed the tpc energy loss and multiplicity distribution)
const string kDataFilename_temp2 = "250337.root"; // Medium data file medium (fixed rotation bkg and competing cascade cuts)
// const string kDataFilename_temp2 = "248760.root"; // Small data file medium (fixed different entries in the pi- and pi+ daughters, multiplicity distribution, number of Ks produced per event)
const string kMCFilename_temp = "211346.root";    // MC file
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
float f1500Mass = 1.522;
float f1500Width = 0.108;
float f1525Mass = 1.518;
float f1525Width = 0.082;
float f1710Mass = 1.710;
float f1710Width = 0.150;
