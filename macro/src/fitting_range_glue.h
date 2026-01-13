#include <vector>

// Projection of signal in different pT bins *************************************

#ifdef DATASET_LHC22o_pass7_small
#ifdef KsKschannel
const Int_t Npt = 3;
// double pT_bins[Npt + 1] = {0.0, 30.0};
// double pT_bins[Npt + 1] = {2.0, 3.0, 5.0, 7.0, 10.0}; // 2022 data
// double pT_bins[Npt + 1] = {2.0, 3.0, 5.0, 7.0, 12.0}; // 2022 data
// double pT_bins[Npt + 1] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0}; // 2023 data
double pT_bins[Npt + 1] = {5.0, 7.0, 10.0, 15.0}; // 2023 data
// double pT_bins[Npt + 1] = {2.0, 12.0};
const int pt_start = 0;
const int pt_end = Npt;

const int Nmult = 4;
double mult_classes[Nmult + 1] = {0.0, 20.0, 50.0, 70.0, 100.0};
#endif
#ifdef KKchannel
const Int_t Npt = 1;
// const Int_t Npt = 8;
double pT_bins[Npt + 1] = {0.0, 30.0};
// float pT_bins[Npt + 1] = {0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 10.0, 20.0, 30.0};
const int pt_start = 0;
const int pt_end = Npt;
#endif
#endif

//******************************************************************************************************************************************************************************************************************************************************************************************************************************************//

#if defined(DATASET_LHC220_pass6_small) || defined(DATASET_LHC22o_pass7_small)
#ifdef KsKschannel
// const std::vector<vector<float>> kNormRangepT = {
//     // 13.6 TeV
//     // {2.20, 2.30}, // 0.0-30.0 GeV/c
//     // {2.80, 2.90}, // 0.0-1.0 GeV/c
//     {2.50, 2.60}, // 1.0-2.0 GeV/c
//     {2.20, 2.30}, // 2.0-3.0 GeV/c
//     {2.20, 2.30}, // 3.0-4.0 GeV/c
//     {2.20, 2.30}, // 4.0-6.0 GeV/c
//     {2.20, 2.30}, // 6.0-12.0 GeV/c
//     {2.20, 2.30}, // 12.0-20.0 GeV/c
//     {2.20, 2.30}, // 12.0-20.0 GeV/c

//     {2.20, 2.30}, // 12.0-20.0 GeV/c
//     {2.20, 2.30}, // 12.0-20.0 GeV/c
// };

const std::vector<std::vector<float>> kNormRangepT = {
    // for cos theta
    {2.50, 2.60}, //
    {2.50, 2.60}, //
    {2.50, 2.60}, //
    {2.50, 2.60}, //
    {2.50, 2.60}, //
    {2.50, 2.60}, //
    {2.50, 2.60}, //
    {2.50, 2.60}, //
    {2.50, 2.60}, //
    {2.50, 2.60}, //
    {2.50, 2.60}, //
};

const std::vector<int> kRebin = {
    // // for pT bins
    2, // 0.0-30.0 GeV/c
    2, // 0.0-1.0 GeV/c
    2, // 1.0-2.0 GeV/c
    2, // 2.0-3.0 GeV/c
    2, // 3.0-4.0 GeV/c
    2, // 4.0-6.0 GeV/c
    2, // 6.0-12.0 GeV/c
    2, // 12.0-20.0 GeV/c
    2, // 12.0-20.0 GeV/c
    2, // 12.0-20.0 GeV/c
    2, // 12.0-20.0 GeV/c

    // // // for cosTheta bins
    // 2, // -1.0 - -0.8
    // 2, // -0.8 - -0.6
    // 2, // -0.6 - -0.4
    // 2, // -0.4 - -0.2
    // 2, // -0.2 - 0.0
    // 2, // 0.0 - 0.2
    // 2, // 0.2 - 0.4
    // 2, // 0.4 - 0.6
    // 2, // 0.6 - 0.8
    // 2, // 0.8 - 1.0

};
#endif
#ifdef KKchannel
const std::vector<vector<double>> kNormRangepT = {
    // 13.6 TeV
    {2.50, 2.60}, // 0.0-30.0 GeV/c
    // {1.05, 1.10}, // 0.0-1.0 GeV/c
    // {2.00, 2.10}, // 1.0-2.0 GeV/c
    // {2.00, 2.10}, // 2.0-3.0 GeV/c
    // {2.00, 2.10}, // 3.0-4.0 GeV/c
    // {2.00, 2.10}, // 4.0-6.0 GeV/c
    // {2.00, 2.10}, // 6.0-10.0 GeV/c
    // {2.00, 2.10}, // 10.0-20.0 GeV/c
    // {2.00, 2.10}, // 20.0-30.0 GeV/c
};
const std::vector<vector<double>> kFitRange = {
    // 13.6 TeV
    {0.79, 1.02}, // 0.0-1.0 GeV/c
    {0.79, 1.02}, // 1.0-2.0 GeV/c
    {0.79, 1.02}, // 2.0-3.0 GeV/c
    {0.79, 1.02}, // 3.0-4.0 GeV/c
    {0.79, 1.02}, // 4.0-6.0 GeV/c
    {0.79, 1.02}, // 6.0-10.0 GeV/c
    {0.79, 1.02}, // 10.0-20.0 GeV/c
    {0.79, 1.02}, // 20.0-30.0 GeV/c

};
const std::vector<int> kRebin = {
    1, // 0.0-30.0 GeV/c
    // 1, // 0.0-15.0 GeV/c
    // 1, // 0.0-1.0 GeV/c
    // 1, // 1.0-2.0 GeV/c
    // 1, // 2.0-3.0 GeV/c
    // 1, // 3.0-4.0 GeV/c
    // 2, // 4.0-6.0 GeV/c
    // 3, // 6.0-10.0 GeV/c
    // 3, // 10.0-20.0 GeV/c
    // 3, // 20.0-30.0 GeV/c
};
#endif
#endif
