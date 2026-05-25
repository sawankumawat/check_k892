
// The variables that can be chaged are here ****************************************************
const string kParticle = "kstar/";
const bool multipanel_plots = 1;
const bool save_plots = 1;
// const string kfoldername_temp = "kstarqa_id21631/hInvMass";
const string kfoldername_temp = "kstarqa/hInvMass";
// const string kfoldername_temp = "lf-kstar892analysis";
// const string kfoldername_temp = "lf-k892analysis";

// const string kvariation = "";
// const string kvariation = "_MC_closure";
// const string kvariation = "_MC_closure_MID0p3";
// const string kvariation = "_MC_closure_MID";
// const string kvariation = "_MC_closure_NoITSROF";
// const string kvariation = "_MC_closure_PVContributor";
// const string kvariation = "_MC_closure_WithoutTOFShift";
// const string kvariation = "_MC_closure_OnlyTPC";
// const string kvariation = "_INEL";
// const string kvariation = "_INELgt0";
// const string kvariation = "_DeepAngle";
// const string kvariation = "_PVContributor";
// const string kvariation = "_LoosePID";
// const string kvariation = "_pTDepPID";
// const string kvariation = "_pTDepPIDTOF";
// const string kvariation = "_MID";
// const string kvariation = "_MID_small";
// const string kvariation = "_MID_verySmall";
// const string kvariation = "_MIDptDep";
// const string kvariation = "_MIDptDep2";
// const string kvariation = "_MIDptDep2_small";
// const string kvariation = "_MIDptDep2_verySmall";
// const string kvariation = "_TOF3";
// const string kvariation = "_TOF3_withoutSquareCut";
// const string kvariation = "_id35679";
// const string kvariation = "_OnlyTPC";
// const string kvariation = "_TOFshift";
// const string kvariation = "_MIDptDep2_0p3_TOF3";
// const string kvariation = "_MIDptDep2_small_TOF3";
// const string kvariation = "_MIDptDep2_TOF3";
// const string kvariation = "_MIDptDep2_verySmall_TOF3";
// const string kvariation = "_MIDNew_TOF2";
// const string kvariation = "_MIDNew_TOF3";
// const string kvariation = "_SquarePID_TOF2";
// const string kvariation = "_SquarePID_TOF3";

//// Systematic variations
// const string kvariation = "";
// const string kvariation = "_FT0C";
// const string kvariation = "_FV0A";
// const string kvariation = "_TPC1p5_combined2";
const string kvariation = "_TPC2p5_combined3p5";
// const string kvariation = "_DCAvar1";
// const string kvariation = "_DCAvar2";
// const string kvariation = "_NoPVContributor";
////********************************************************************************************

// define datasets here
#define DATASET_LHC220_pass7

const string kDataFilename_temp1 = "../data/" + kParticle;  // data file
const string kSignalOutput_temp = "../output/" + kParticle; // output folder

#ifdef DATASET_LHC220_pass7
const string kDataset_temp = "LHC22o_pass7/";
// const string kDataset_temp = "LHC22o_pass7/MC_closure/";
// const string kDataset_temp = "LHC22o_pass7/checks/";
// const string kDataset_temp = "LHC22o_pass7/Occupancy_effect/";
// const string kDataset_temp = "LHC22o_pass7/IR_study/";
#endif

#ifdef DATASET_LHC220_pass7


//*****************************After SQM***************************************
//==========2023 data===========
// const string kDataFilename_temp2 = "660453.root"; // (with All checks as above)
// const string kDataFilename_temp2 = "662039.root"; // (No RCT)
// const string kDataFilename_temp2 = "670168.root"; // (Base, INEL)

//===========MC closure===========
// const string kDataFilename_temp2 = "666966.root"; // pp 5.36 dataset  (MC for closure, mistakenly used TOFFT0)
// const string kDataFilename_temp2 = "657468.root"; // (MC for closure)
// const string kDataFilename_temp2 = "664785.root"; // (MC for closure)
// const string kDataFilename_temp2 = "665348.root"; // (MC for closure with TOF shift)
// const string kDataFilename_temp2 = "667890.root"; // (MC for closure with higher TOF shift)
// const string kDataFilename_temp2 = "669655.root"; // (MC_closure)
// const string kDataFilename_temp2 = "673285.root"; // (TOF3: MC_closure, MC_closure_INEL, MC_closure_MID0p3)
// const string kDataFilename_temp2 = "674418.root"; // (TOF3 with checks on Mother: MC_closure, MC_closure_INEL, MC_closure_MID0p3, MC_closure_MID, MC_closure_NoITSROF, MC_closure_PVContributor, MC_closure_WithoutTOFShift)
// const string kDataFilename_temp2 = "677471.root"; // (MC_closure, MC_closure_INEL, MC_closure_MID0p3, MC_closure_MID, MC_closure_NoITSROF, MC_closure_PVContributor, MC_closure_WithoutTOFShift, MC_closure_OnlyTPC)

//=========Other Checks============
// const string kDataFilename_temp2 = "655628.root"; // (Base, INEL)
// const string kDataFilename_temp2 = "658307.root"; // (LoosePID, pTDepPID, pTDepPIDTOF)
// const string kDataFilename_temp2 = "661905.root"; // (Event_Time dependency check in TOF)
// const string kDataFilename_temp2 = "658306.root"; // (DeepAngle, PV Contributor, INELgt0)
// const string kDataFilename_temp2 = "668039.root"; // (With square PID: Base, MID, MIDptDep2)
// const string kDataFilename_temp2 = "672297.root"; // (Base (2sigma TOF), MIDptDep2_0p3_TOF3, MIDptDep2_small_TOF3, MIDptDep2_TOF3, MIDptDep2_verySmall_TOF3)
// const string kDataFilename_temp2 = "675391.root"; // (MIDNew_TOF2, MIDNew_TOF3, SquarePID_TOF2, SquarePID_TOF3)
// const string kDataFilename_temp2 = "668605.root"; // (Base, MID, MID_small, MID_verySmall, MIDptDep, MIDptDep2, MIDptDep2_small, MIDptDep2_verySmall (With Square PID): TOF3, TOF3_withoutSquareCut)

//==========2024 data===========
// const string kDataFilename_temp2 = "663738.root"; // (Base, OnlyTPC)
// const string kDataFilename_temp2 = "660943.root"; // (Base, LoosePID, pTDepPID, pTDepPIDTOF, MIDptDep, MID)
// const string kDataFilename_temp2 = "664559.root"; // (Base, INEL, TOFshift, TOFshiftMID)
const string kDataFilename_temp2 = "679906.root"; // (Sys. train: Base (3sigma TOF), FT0C, FV0A, TPC1p5_combined2, TPC2p5_combined3p5)
// const string kDataFilename_temp2 = "682963.root"; // (Sys. train2: DCAvar1, DCAvar2, NoPVContributor)
#endif

// final dataset name
const string kDataset = kDataFilename_temp1 + kDataset_temp;
const string kSignalOutput = kSignalOutput_temp + kDataset_temp + kDataFilename_temp2.substr(0, kDataFilename_temp2.rfind("."));
const string kDataFilename = kDataset + kDataFilename_temp2;
const string kfoldername = kfoldername_temp.substr(0, kfoldername_temp.length() - 9) + kvariation + kfoldername_temp.substr(kfoldername_temp.length() - 9, kfoldername_temp.length());
const string koutputfolder = kSignalOutput + "/" + kfoldername;

// Canvas dimensions
const int klowerpad = 5;
const int kupperpad = 4;
// const int kcanvaswidth = 1440 * 2;
// const int kcanvasheight = 720 * 2;
const int kcanvaswidth = 1440;
const int kcanvasheight = 1080;
const int kcanvasdivide[2] = {klowerpad, kupperpad};

float masspdg = 0.896;   // in GeV/c^2
float widthpdg = 0.0487; // in 1 sigma GeV/c^2

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

//================Variation names===================
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
// const string kvariation = "_PIDKa1_NoMID"; // change the variation here
// const string kvariation = "_PIDKa2"; // change the variation here
// const string kvariation = "_BetaTOF0p5";      // change the variation here
// const string kvariation = "_GoodFT0vsPV";      // change the variation here
// const string kvariation = "_GoodITSLayersAll"; // change the variation here
// const string kvariation = "_ITSTPCRefit";       // change the variation here
// const string kvariation = "_VertexITSTPC";      // change the variation here
// const string kvariation = "_VertexTOFMatched";  // change the variation here
// const string kvariation = "_ptDepPID";  // change the variation here
// const string kvariation = "_NoRCT";  // change the variation here
// const string kvariation = "_hasITS";  // change the variation here

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
// const string kDataFilename_temp2 = "480657.root"; // LHC24_pass1_minBias dataset, INEL > 0

//****************************QA checks all ***************************************************
// const string kDataFilename_temp2 = "585940.root"; // LHC23_pass4_thin_small dataset, INEL > 0

//*****************************pT-dependent PID***************************************
// const string kDataFilename_temp2 = "586976.root"; //23 dataset
// const string kDataFilename_temp2 = "586385.root"; // 24 dataset
