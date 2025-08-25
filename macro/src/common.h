
// The variables that can be chaged are here ****************************************************
const string kParticle = "kstar/";
const bool multipanel_plots = 0;
const bool save_plots = 1;
// const string kfoldername_temp = "kstarqa_id21631/hInvMass";
const string kfoldername_temp = "kstarqa/hInvMass";
// const string kfoldername_temp = "lf-kstar892analysis";
// const string kfoldername_temp = "lf-k892analysis";

// const string kvariation = ""; // change the variation here
// const string kvariation = "_CutsOnMotherParticle"; // change the variation here
// const string kvariation = "_PtDepDCAxy"; // change the variation here
// const string kvariation = "_RCT"; // change the variation here
// const string kvariation = "_onlyTPC"; // change the variation here
// const string kvariation = "_Min5ItsClusters"; // change the variation here
// const string kvariation = "_PbPbCuts"; // change the variation here
// const string kvariation = "_PtDepPID"; // change the variation here
// const string kvariation = "_FakeTracks_id33593"; // change the variation here
// const string kvariation = "_MID_id33593"; // change the variation here
// const string kvariation = "_TPCChi2Min_id34810"; // change the variation here
// const string kvariation = "_TrackRapidity0p3_id34810"; // change the variation here
// const string kvariation = "_MCclosure_id34109"; // change the variation here
// const string kvariation = "_PIDKa2"; // change the variation here
// const string kvariation = "_PIDKa1"; // change the variation here
const string kvariation = "_PIDKa1"; // change the variation here
////********************************************************************************************

// define datasets here
#define DATASET_LHC220_pass7

const string kDataFilename_temp1 = "../data/" + kParticle;  // data file
const string kSignalOutput_temp = "../output/" + kParticle; // output folder

#ifdef DATASET_LHC220_pass7
const string kDataset_temp = "LHC22o_pass7/";
// const string kDataset_temp = "LHC22o_pass7/checks/";
// const string kDataset_temp = "LHC22o_pass7/Occupancy_effect/";
// const string kDataset_temp = "LHC22o_pass7/IR_study/";
// const string kDataset_temp = "LHC22o_pass7/MC_closure/";
const string kMCDataset = "../mc/LHC24f3b/";
#endif

#ifdef DATASET_LHC220_pass7

//*****************************Corrected TPC crossed rows********************************
// const string kDataFilename_temp2 = "459845.root"; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "459908.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "460233.root"; // LHC24_pass1_minBias dataset, INEL > 0

//**********************************MID cuts***************************************
// const string kDataFilename_temp2 = "477779.root"; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "478015.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "477833.root"; // LHC24_pass1_minBias dataset, INEL > 0

//*************************PID Variations for Kaon (without MID)**************************
// const string kDataFilename_temp2 = "480317.root"; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "480447.root"; // LHC23_pass4_thin_small dataset, INEL > 0
const string kDataFilename_temp2 = "480657.root"; // LHC24_pass1_minBias dataset, INEL > 0

#endif

// final dataset name
const string kDataset = kDataFilename_temp1 + kDataset_temp;
const string kSignalOutput = kSignalOutput_temp + kDataset_temp + kDataFilename_temp2.substr(0, kDataFilename_temp2.rfind("."));
const string kDataFilename = kDataset + kDataFilename_temp2;
const string kfoldername = kfoldername_temp.substr(0, kfoldername_temp.length() - 9) + kvariation + kfoldername_temp.substr(kfoldername_temp.length() - 9, kfoldername_temp.length());
const string koutputfolder = kSignalOutput + "/" + kfoldername;

// Canvas dimensions
const int klowerpad = 4;
const int kupperpad = 4;
const int kcanvaswidth = 1440 * 2;
const int kcanvasheight = 720 * 2;
const int kcanvasdivide[2] = {klowerpad, kupperpad};

float masspdg = 0.896;   // in GeV/c^2
float widthpdg = 0.0473; // in 1 sigma GeV/c^2

////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
//                                   OLD Checks                                   //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////
//******************************MC Closure******************************/
// const string kDataFilename_temp2 = "474414.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "476106.root"; // LHC23_pass4_thin_small dataset, INEL > 0

//*****************************pp (reference) for OO, NeNe*************************************
// const string kDataFilename_temp2 = "477291.root"; // LHC23_pass4_thin_small dataset, INEL > 0

//*************************IR study********************************
// ********************2023 dataset******************************
// const string kDataFilename_temp2 = "463114.root"; // 1-2 MHz (combined 1-1.3 and 2 MHz)
// const string kDataFilename_temp2 = "536899.root"; // 1-1.3 MHz
// const string kDataFilename_temp2 = "537861.root"; // 2 MHz
// const string kDataFilename_temp2 = "535069.root"; // 14 kHz
// const string kDataFilename_temp2 = "535545.root"; // 70 kHz
// const string kDataFilename_temp2 = "535645.root"; // 135 kHz
// const string kDataFilename_temp2 = "535999.root"; // 330 kHz
// const string kDataFilename_temp2 = "LHC23z.root"; // 450 kHz
// const string kDataFilename_temp2 = "LHC23ls.root"; // 650 kHz
// ********************2024 dataset******************************

//******************* Occupancy cut study *********************************
// const string kDataFilename_temp2 = "466154.root"; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "696969.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "466180.root"; // LHC24_pass1_minBias dataset, INEL > 0

//*****Checks in data: Cuts on Mothers particle, only TPC, PID Ka2, pT dependent DCAxy, RCT ***********
// const string kDataFilename_temp2 = "468837.root"; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "468791.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "468697.root"; // LHC24_pass1_minBias dataset, INEL > 0

//*****Checks in data: MinClusters ITS > 5, PbPb cuts, Pt dependent PID*******************
// const string kDataFilename_temp2 = "472720.root"; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "470913.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "471898.root"; // LHC24_pass1_minBias dataset, INEL > 0

//*******Checks in data: Fake tracks, MID, TPCChi2MinCut, TrackRapidityCut *****************
// const string kDataFilename_temp2 = "473153.root"; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "473237.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "473185.root"; // LHC24_pass1_minBias dataset, INEL > 0

//*************************PID Variations for Kaon (without MID, multcentTable)**************************
// const string kDataFilename_temp2 = ""; // LHC22_pass7_medium dataset, INEL > 0
// const string kDataFilename_temp2 = "480448.root"; // LHC23_pass4_thin_small dataset, INEL > 0
// const string kDataFilename_temp2 = "480358.root"; // LHC24_pass1_minBias dataset, INEL > 0

//*************************ItsTpcTracksCheck, betacutTOF******************************
// const string kDataFilename_temp2 = "481941.root"; // LHC23_pass4_thin_small dataset, INEL > 0, No effect is seen