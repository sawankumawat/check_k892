

// Projection of signal in different pT bins *************************************

#ifdef DATASET_LHC220_pass6_small
const Int_t Npt = 1;
double pT_bins[Npt + 1] = {0.0, 15.0};
const int pt_start = 0;
const int pt_end = Npt;
#endif

//******************************************************************************************************************************************************************************************************************************************************************************************************************************************//

#ifdef DATASET_LHC220_pass6_small

const std::vector<vector<double>> kNormRangepT = {
    // 13.6 TeV
    {1.10, 1.15}, // 0.0-15.0 GeV/c
};
const std::vector<vector<double>> kFitRange = {
    // 13.6 TeV
    {0.79, 1.02}, // 0.0-15.0 GeV/c

};

const std::vector<int> kRebin = {
    1, // 0.0-15.0 GeV/c
};
#endif