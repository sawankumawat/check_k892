// Description: This file contains the fitting range for the different datasets and the different pT bins


#ifdef DATASET_LHC22cde
const Int_t Npt = 8;
const int pt_start = 0;
const int pt_end = 8;
#endif

#ifdef DATASET_LHC23_zzh_pass4
const Int_t Npt = 8;
const int pt_start = 0;
const int pt_end = 8;
#endif

// Projection of signal in different pT bins *************************************

#ifdef DATASET_LHC220_pass6_small
const Int_t Npt = 27;
double pT_bins[Npt + 1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0};
#endif

#ifdef DATASET_LHC220_pass6_small
const int pt_start = 0;
const int pt_end = Npt;
#endif
// for pp data low and high IR
// double pT_bins[Npt + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 8.0};

#ifdef DATASET_LHC23_zzh_pass4
double pT_bins[Npt + 1] = {0.6, 1.2, 1.8, 2.4, 3.0, 4.0, 6.0, 10.0};
#endif

#ifdef DATASET_LHC22cde
double pT_bins[Npt + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 5.0};
#endif

//****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************

// Normalization range *************************** Normalization range ******************************** Normalization range
//*************************************************pbpb datasets**************************************
// // dataset LHC23zzf_pass2_QC
// double lownorm[30]  = {1.1, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40};
// double highnorm[30] = {1.15, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45};

// // //  dataset LHC23zzh_pass2_small & LHC23zzg_apass2
// double lownorm[30]  = {0.7, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40};
// double highnorm[30] = {0.73, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45};

//******************************pp datasets ****************************************

// // dataset LHC23_pass1_lowB_highIR_sampling
// double lownorm[30] = {0.7, 0.68, 0.68, 0.68, 1.1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
// double highnorm[30] = {0.73, 0.73, 0.73, 0.73, 1.2, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

// dataset LHC23_pass1_lowB_lowIR
// double lownorm[30] = {1, 1, 0.68, 0.68, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
// double highnorm[30] = {1.05, 1.05, 0.73, 0.73, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

// LHC23zf
// double lownorm[30] = {1.05, 0.68, 0.68, 0.68, 1.1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
// double highnorm[30] = {1.1, 0.73, 0.73, 0.73, 1.2, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

// // LHC23f, LHC23h
// double lownorm[30] = {1.05, 0.68, 0.68, 0.68, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
// double highnorm[30] = {1.1, 0.73, 0.73, 0.73, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

// // LHC23r, LHC23zzs
// double lownorm[30] = {1.05, 0.68, 0.68, 0.68, 1.05, 1.05, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
// double highnorm[30] = {1.1, 0.73, 0.73, 0.73, 1.1, 1.1, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

// 900 GeV
// double lownorm[30] = {0.72, 0.72, 0.72, 0.72, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
// double highnorm[30] = {0.75, 0.75, 0.75, 0.75, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

// #ifdef DATASET_LHC220_pass6_small
// double lownorm[40] = {1.05, 1.05, 1.05, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20};
// double highnorm[40] = {1.10, 1.10, 1.10, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30};
// #endif

//****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************

// Fitting range ******************************* Fitting range ******************************* Fitting range**********

//*******************************************MIXED EVENTS****************************************

//****************************************** pbpb datasets ************************************************
// // dataset LHC23zzf_pass2_QC
// double lowfitrangeme[30]  = {0.86, 0.85, 0.81, 0.8, 0.78, 0.78, 0.76, 0.82, 0.82, 0.73, 0.73, 0.7, 0.75, 0.75, 0.75, 0.73, 0.7, 0.67, 0.75, 0.75, 0.75, 0.75};
// double highfitrangeme[30] = {0.95, 0.95, 0.997, 1.0, 1.05, 1.02, 1.02, 0.95, 0.95, 1.05, 1.05, 1.1, 1.1, 1.1, 1.1, 1.15, 1.12, 1.13, 1.1, 1.1, 1.1, 1.1, 1.1};

// //  dataset LHC23zzh_pass2_small
// double lowfitrangeme[30]  = {0.81, 0.85, 0.82, 0.8, 0.78, 0.78, 0.76, 0.82, 0.82, 0.73, 0.73, 0.7, 0.75, 0.75, 0.75, 0.73, 0.7, 0.67, 0.75, 0.75, 0.75, 0.75};
// double highfitrangeme[30] = {0.94, 0.95, 0.95, 1.0, 1.05, 1.02, 1.02, 0.95, 0.95, 1.05, 1.05, 1.1, 1.1, 1.1, 1.1, 1.15, 1.12, 1.13, 1.1, 1.1, 1.1, 1.1, 1.1};

// // dataset LHC23zzg_apass2
// double lowfitrangeme[30]  = {0.83, 0.86, 0.84, 0.8, 0.78, 0.76, 0.76, 0.82, 0.82, 0.73, 0.73, 0.7, 0.75, 0.75, 0.75, 0.73, 0.7, 0.67, 0.75, 0.75, 0.75, 0.75};
// double highfitrangeme[30] = {0.95, 0.95, 0.97, 1.0, 1.0, 1.02, 1.02, 0.95, 0.95, 1.05, 1.05, 1.1, 1.1, 1.1, 1.1, 1.15, 1.12, 1.13, 1.1, 1.1, 1.1, 1.1, 1.1};

//*****************************************************pp datasets ****************************************
// ****************************************************pp datasets ****************************************

// // // dataset LHC23_pass1_lowB_highIR_sampling
// //pol2
//  double lowfitrangeme[30]   =  {0.72, 0.68, 0.8, 0.8, 0.77, 0.785, 0.78, 0.74, 0.77, 0.78, 0.78, 0.7, 0.75, 0.7, 0.7, 0.72, 0.73, 0.7, 0.7, 0.73, 0.7, 0.7};
// double highfitrangeme[30]   =  {1.05, 1.05, 1.05, 1.05, 1.08, 1.025, 1.05, 1.04, 1.04, 1.04, 1.02, 1.05, 1.1, 1.08, 1.1, 1.1, 1.1, 1.04, 1.06, 1.1, 1.1, 1.1, 1.1};

// pol3
//   double lowfitrangeme[30]   =  {0.7, 0.68, 0.75, 0.75, 0.75, 0.75, 0.78, 0.74, 0.765, 0.775, 0.71, 0.7, 0.75, 0.76, 0.74, 0.73, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7};
//  double highfitrangeme[30]   =  {1.04, 1.05, 1.03, 1.05, 1.09, 1.02, 1.05, 1.04, 1.03, 1.03, 1.02, 1.05, 1.08, 1.05, 1.07, 1.1, 1.1, 1.04, 1.06, 1.1, 1.1, 1.1, 1.1};

// dataset LHC23_pass1_lowB_lowIR
//   double lowfitrangeme[30] = {0.75, 0.65, 0.78, 0.78, 0.72, 0.72, 0.73, 0.73, 0.75, 0.76, 0.76, 0.7, 0.75, 0.75, 0.73, 0.73, 0.73, 0.68, 0.8, 0.75, 0.75, 0.75};
// double highfitrangeme[30] = {1.1, 1.06, 1.1, 1.06, 1.04, 1.04, 1.04, 1.04, 1.04, 1.05, 1.05, 1.1, 1.1, 1.1, 1.1, 1.15, 1.12, 1.05, 1.0, 1.1, 1.1, 1.1, 1.1};

// pol3
//  double lowfitrangeme[30]   =  {0.7, 0.68, 0.75, 0.75, 0.75, 0.75, 0.78, 0.74, 0.765, 0.775, 0.71, 0.7, 0.75, 0.76, 0.74, 0.73, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7};
// double highfitrangeme[30]   =  {1.04, 1.05, 1.03, 1.05, 1.09, 1.02, 1.03, 1.04, 1.03, 1.03, 1.02, 1.05, 1.08, 1.05, 1.07, 1.1, 1.1, 1.04, 1.06, 1.1, 1.1, 1.1, 1.1};

// // dataset LHC23zf
// pol2
// double lowfitrangeme[30] =  {0.72, 0.68, 0.74, 0.76, 0.74, 0.74, 0.74, 0.74, 0.77, 0.78, 0.74, 0.70, 0.75, 0.70, 0.70, 0.72, 0.73, 0.70, 0.70, 0.73, 0.70, 0.70};
// double highfitrangeme[30] = {1.05, 1.05, 1.05, 1.05, 1.08, 1.05, 1.05, 1.04, 1.04, 1.04, 1.05, 1.05, 1.10, 1.08, 1.10, 1.10, 1.10, 1.04, 1.06, 1.10, 1.10, 1.10};

// // dataset LHC23zm
// // pol2
// double lowfitrangeme[30] =  {0.7, 0.68, 0.74, 0.78, 0.74, 0.74, 0.74, 0.74, 0.77, 0.78, 0.74, 0.70, 0.75, 0.70, 0.70, 0.72, 0.73, 0.71, 0.72, 0.73, 0.74, 0.70};
// double highfitrangeme[30] = {1.01, 1.05, 1.05, 1.03, 1.08, 1.05, 1.05, 1.04, 1.04, 1.04, 1.05, 1.05, 1.10, 1.08, 1.10, 1.10, 1.10, 1.03, 1.06, 1.10, 1.05, 1.10};

// // // dataset LHC23f
// // pol3
// double lowfitrangeme[30] =  {0.7, 0.68, 0.74, 0.76, 0.74, 0.76, 0.77, 0.775, 0.77, 0.78, 0.76, 0.70, 0.75, 0.70, 0.72, 0.74, 0.73, 0.71, 0.72, 0.73, 0.74, 0.70};
// double highfitrangeme[30] = {1.05, 1.05, 1.05, 1.03, 1.08, 1.02, 1.02, 1.02, 1.04, 1.04, 1.05, 1.05, 1.10, 1.08, 1.10, 1.10, 1.10, 1.03, 1.06, 1.10, 1.05, 1.10};

// // // dataset LHC23h
// // pol3
// double lowfitrangeme[30] =  {0.7, 0.68, 0.74, 0.72, 0.74, 0.76, 0.77, 0.775, 0.77, 0.78, 0.76, 0.70, 0.75, 0.70, 0.72, 0.72, 0.73, 0.71, 0.72, 0.73, 0.74, 0.70};
// double highfitrangeme[30] = {1.05, 1.05, 1.05, 1.06, 1.08, 1.02, 1.02, 1.02, 1.04, 1.04, 1.05, 1.05, 1.10, 1.08, 1.10, 1.10, 1.10, 1.03, 1.06, 1.10, 1.05, 1.10};

// // dataset LHC23r, LHC23zzs
// pol3
// double lowfitrangeme[30] =  {0.72, 0.72, 0.74, 0.73, 0.74, 0.76, 0.78, 0.78, 0.77, 0.78, 0.77, 0.70, 0.75, 0.70, 0.77, 0.75, 0.73, 0.71, 0.72, 0.73, 0.74, 0.70};
// double highfitrangeme[30] = {1.05, 1.07, 1.05, 1.06, 1.01, 1.01, 1.02, 1.02, 1.04, 1.04, 1.05, 1.05, 1.07, 1.07, 1.05, 1.10, 1.10, 1.03, 1.06, 1.10, 1.05, 1.10};

// // 900 GeV
// double lowfitrangeme[30] =  {0.72, 0.73, 0.75, 0.75, 0.75, 0.75, 0.75, 0.78, 0.75, 0.78, 0.75, 0.75, 0.75, 0.70, 0.77, 0.75, 0.73, 0.71, 0.72, 0.73, 0.74, 0.70};
// double highfitrangeme[30] = {1.1, 1.04, 1.06, 1.09, 1.04, 1.04, 1.04, 1.04, 1.04, 1.1, 1.1, 1.1, 1.1, 1.07, 1.05, 1.10, 1.10, 1.03, 1.06, 1.10, 1.05, 1.10};

#ifdef DATASET_LHC220_pass6_small
vector<vector<float>> kNormRangepT = {
    // run2 pT bins
    // 13.6 TeV
    {1.05, 1.10}, // 0.0-0.1
    {1.05, 1.10}, // 0.1-0.2
    {1.05, 1.10}, // 0.2-0.3
    {1.10, 1.20}, // 0.3-0.4
    {1.10, 1.20}, // 0.5-0.5
    {1.10, 1.20}, // 0.5-0.6
    {1.10, 1.20}, // 0.6-0.7
    {1.10, 1.20}, // 0.7-0.8
    {1.10, 1.20}, // 0.8-0.9
    {1.10, 1.20}, // 0.9-1.0
    {1.10, 1.20}, // 1.0-1.2
    {1.10, 1.20}, // 1.2-1.4
    {1.10, 1.20}, // 1.4-1.6
    {1.10, 1.20}, // 1.6-1.8
    {1.10, 1.20}, // 1.8-2.0
    {1.10, 1.20}, // 2.0-2.4
    // {1.10, 1.20}, // 2.2-2.4
    {1.10, 1.20}, // 2.4-2.8
    // {1.10, 1.20}, // 2.6-2.8
    {1.10, 1.20}, // 2.8-3.2
    {1.10, 1.20}, // 3.2-3.6
    {1.10, 1.20}, // 3.6-4.0
    {1.10, 1.20}, // 4.0-5.0
    {1.10, 1.20}, // 5.0-6.0
    {1.10, 1.20}, // 6.0-7.0
    {1.10, 1.20}, // 7.0-8.0
    {1.10, 1.20}, // 8.0-10.0
    {1.10, 1.20}, // 10.0-12.0
    {1.10, 1.20}, // 12.0-15.0
    {1.10, 1.20}, // 15.0-20.0
    {1.10, 1.20}, // 20.0-25.0
    {1.10, 1.20}, // 25.0-30.0
};
vector<vector<float>> kFitRange = {
    // run2 pT bins
    // ME background
    // 13.6 TeV
    {0.76, 1.03},  // 0.0-0.1
    {0.75, 1.04},  // 0.1-0.2
    {0.76, 1.04},  // 0.2-0.3
    {0.76, 1.05},  // 0.3-0.4
    {0.77, 0.99},  // 0.4-0.5
    {0.77, 0.99},  // 0.5-0.6
    {0.78, 0.99},  // 0.6-0.7
    {0.78, 0.99},  // 0.7-0.8
    {0.79, 0.99},  // 0.8-0.9
    {0.79, 0.99},  // 0.9-1.0
    {0.795, 0.99}, // 1.0-1.2
    {0.80, 0.99},  // 1.2-1.4
    {0.80, 0.99},  // 1.4-1.6
    {0.808, 0.99}, // 1.6-1.8
    {0.82, 0.99},  // 1.8-2.0
    {0.83, 0.99},  // 2.0-2.4
    // {0.825, 0.99}, // 2.2-2.4
    {0.83, 0.995}, // 2.4-2.8
    // {0.83, 0.99},  // 2.6-2.8
    {0.83, 0.99},  // 2.8-3.2
    {0.83, 0.99},  // 3.2-3.6
    {0.83, 0.99},  // 3.6-4.0
    {0.825, 0.97}, // 4.0-5.0
    {0.81, 0.99},  // 5.0-6.0
    {0.83, 0.99},  // 6.0-7.0
    {0.80, 0.99},  // 7.0-8.0
    {0.80, 0.99},  // 8.0-10.0
    {0.78, 1.00},  // 10.0-12.0
    {0.76, 1.04},  // 12.0-15.0
    {0.7, 1.08},   // 15.0-20.0
    {0.68, 1.1},   // 20.0-25.0
    {0.7, 1.08},   // 20.0-30.0
};

const std::vector<int> kRebin = {
    1, // 0.0-0.1
    1, // 0.1-0.2
    1, // 0.2-0.3
    1, // 0.3-0.4
    1, // 0.4-0.5
    1, // 0.5-0.6
    1, // 0.6-0.7
    1, // 0.7-0.8
    1, // 0.8-0.9
    1, // 0.9-1.0
    1, // 1.0-1.2
    1, // 1.2-1.4
    1, // 1.4-1.6
    1, // 1.6-1.8
    1, // 1.8-2.0
    1, // 2.0-2.4
    1, // 2.4-2.8
    1, // 2.8-3.2
    1, // 3.2-3.6
    1, // 3.6-4.0
    1, // 4.0-5.0
    1, // 5.0-6.0
    1, // 6.0-7.0
    1, // 7.0-8.0
    1, // 8.0-10.0
    1, // 10.0-12.0
    1, // 12.0-15.0
    2, // 15.0-20.0
    4, // 20.0-25.0
    6, // 25.0-30.0
};
#endif