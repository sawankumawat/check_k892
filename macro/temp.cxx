#include "src/style.h"
void temp()
{
    TFile *fSysUncert = new TFile("../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/SystematicsPlots/SysUncert.root", "READ");
    if (fSysUncert->IsZombie())
    {
        cout << "Systematic uncertainty file not found" << endl;
        return;
    }
    TH1D *hTotalSysSmoothed = (TH1D *)fSysUncert->Get("hTotalSysSmoothed_0_100"); // Temporary assigning same to all multiplicity classes
    if (hTotalSysSmoothed == nullptr)
    {
        cout << "Histogram hTotalSysSmoothed_0_100 not found in the systematic uncertainty file" << endl;
        return;
    }
    int totalBins = hTotalSysSmoothed->GetNbinsX();
    for (int i = 1; i <= totalBins; ++i)
    {
        double binCenter = hTotalSysSmoothed->GetBinCenter(i);
        double binLowEdge = hTotalSysSmoothed->GetBinLowEdge(i);
        double binHighEdge = hTotalSysSmoothed->GetBinLowEdge(i + 1);
        double binWidth = binHighEdge - binLowEdge;
        double binContent = hTotalSysSmoothed->GetBinContent(i);
        cout << "Bin " << i << ": pT bin = " << binLowEdge << " - " << binHighEdge << ", Sys Uncert = " << binContent << endl;
    }

    const int Npt = 34; // INEL with even higher pT bins
    double pT_bins[Npt + 1] = {0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0}; // INEL with even higher pT bins

    TH1D *hSysINELTemp = new TH1D("hSysINELTemp", "Systematic Uncertainty for INEL; pT (GeV/c); Systematic Uncertainty", Npt, pT_bins);
    double sysUncertValues[Npt] = {0.0815995, 0.0806168, 0.0796341, 0.0786513, 0.0776686, 0.0747205, 0.0698069, 0.0648933, 0.0599798, 0.0550662, 0.0521791, 0.0513185, 0.0504578, 0.0495972, 0.0496780, 0.0507001, 0.0520358, 0.0536851, 0.0539110, 0.0521150, 0.0423910, 0.0356319, 0.0370027, 0.0481308, 0.0622399, 0.0615530, 0.0677056, 0.0818712, 0.0883558, 0.0861680, 0.0926845, 0.0999785, 0.1090960, 0.1122135};
    for (int i = 0; i < Npt; ++i)
    {
        hSysINELTemp->SetBinContent(i + 1, sysUncertValues[i]);
    }
    TFile *fSysINELTemp = new TFile("../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/SystematicsPlots/SysUncert_INELTemp.root", "RECREATE");
    hSysINELTemp->Write();
    fSysINELTemp->Close();
}

// // pT bin : Estimated systematic uncertainty
// 0.00 - 0.02   : 0.0815995
// 0.02 - 0.04   : 0.0806168
// 0.04 - 0.06   : 0.0796341
// 0.06 - 0.08   : 0.0786513
// 0.08 - 0.10   : 0.0776686
// 0.10 - 0.20   : 0.0747205
// 0.20 - 0.30   : 0.0698069
// 0.30 - 0.40   : 0.0648933
// 0.40 - 0.50   : 0.0599798
// 0.50 - 0.60   : 0.0550662
// 0.60 - 0.70   : 0.0521791
// 0.70 - 0.80   : 0.0513185
// 0.80 - 0.90   : 0.0504578
// 0.90 - 1.00   : 0.0495972
// 1.00 - 1.20   : 0.0496780
// 1.20 - 1.40   : 0.0507001
// 1.40 - 1.60   : 0.0520358
// 1.60 - 1.80   : 0.0536851
// 1.80 - 2.00   : 0.0539110
// 2.00 - 2.40   : 0.0521150
// 2.40 - 2.80   : 0.0423910
// 2.80 - 3.20   : 0.0356319
// 3.20 - 3.60   : 0.0370027
// 3.60 - 4.00   : 0.0481308
// 4.00 - 5.00   : 0.0622399
// 5.00 - 6.00   : 0.0615530
// 6.00 - 7.00   : 0.0677056
// 7.00 - 8.00   : 0.0818712
// 8.00 - 10.00  : 0.0883558
// 10.00 - 12.00 : 0.0861680
// 12.00 - 15.00 : 0.0926845
// 15.00 - 20.00 : 0.0999785
// 20.00 - 25.00 : 0.1090960
// 25.00 - 30.00 : 0.1122135