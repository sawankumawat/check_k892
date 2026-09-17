using namespace std;

void saveKstarSpectra()
{
    TString SpectraPath = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED";
    TString UncertPath = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/SystematicsPlots";

    TFile *fTotalUncert = new TFile(UncertPath + "/SysUncert.root", "read");
    TFile *fUncorrUncert = new TFile(UncertPath + "/UnCorrSystematics.root", "read");
    if (fTotalUncert->IsZombie() || fUncorrUncert->IsZombie())
    {
        cerr << "Error: Systematics files not found!" << endl;
        return;
    }

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    int multBins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1; // number of multiplicity bins

    TFile *fSpectraSave = new TFile(SpectraPath + "/KstarSpectra.root", "RECREATE");

    for (int imult = 0; imult < multBins; imult++)
    {
        int multlow = mult_classes[imult];
        int multhigh = mult_classes[imult + 1];

        TFile *fSpectra = new TFile(SpectraPath + Form("/corrected_spectra_%d_%d.root", multlow, multhigh), "read");
        if (fSpectra->IsZombie())
        {
            cerr << "Error: Spectra file for multiplicity bin " << multlow << "-" << multhigh << " not found!" << endl;
            continue;
        }

        TH1F *hSpectra = (TH1F *)fSpectra->Get(Form("mult_%d-%d/corrected_spectra_Integral_final", multlow, multhigh));
        if (!hSpectra)
        {
            cerr << "Error: Histogram 'corrected_spectra_Integral_final' not found in file for multiplicity bin " << multlow << "-" << multhigh << endl;
            continue;
        }

        TH1F *hTotalUncert = (TH1F *)fTotalUncert->Get(Form("hTotalSysSmoothed_%d_%d", multlow, multhigh));
        TH1F *hUncorrUncert = (TH1F *)fUncorrUncert->Get(Form("hUncorrelatedUncertaintySmoothed_mult_%d_%d", multlow, multhigh));
        if (!hTotalUncert || !hUncorrUncert)
        {
            cerr << "Error: Uncertainty histograms not found for multiplicity bin " << multlow << "-" << multhigh << endl;
            continue;
        }

        TH1F *hSpectraTotalUncert = (TH1F *)hSpectra->Clone(Form("corrected_spectra_Integral_final_TotalUncert_%d_%d", multlow, multhigh));
        TH1F *hSpectraUncorrUncert = (TH1F *)hSpectra->Clone(Form("corrected_spectra_Integral_final_UncorrUncert_%d_%d", multlow, multhigh));

        for (int i = 1; i <= hSpectra->GetNbinsX(); i++)
        {
            double content = hSpectra->GetBinContent(i);
            double totalUncert = hTotalUncert->GetBinContent(i);
            double uncorrUncert = hUncorrUncert->GetBinContent(i);

            hSpectraTotalUncert->SetBinError(i, totalUncert * content);
            hSpectraUncorrUncert->SetBinError(i, uncorrUncert * content);
        }


        fSpectraSave->cd();
        hSpectra->Write(Form("SpectraStat_%d_%d", multlow, multhigh));
        hSpectraTotalUncert->Write(Form("SpectraSys_%d_%d", multlow, multhigh));
        hSpectraUncorrUncert->Write(Form("SpectraUncorrSys_%d_%d", multlow, multhigh));
    }
}