#include "spectra/YieldMean.C"
#include "spectra/ReweightEfficiency.C"
#include "src/style.h"
#include "src/initializations.h"

Double_t FuncLavy(Double_t *x, Double_t *par)
{

    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

void reweightFactor()
{
    string path = "../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass"; // 2024 data
    TFile *fspectra = new TFile((path + "/corrected_spectra.root").c_str(), "read");
    if (fspectra->IsZombie())
    {
        cout << "File not found" << endl;
        return;
    }

    TString reweightEffPath = path + "/ReweightEfficiency";
    if (gSystem->mkdir(reweightEffPath, kTRUE))
    {
        std::cout << "Folder " << reweightEffPath << " created successfully." << std::endl;
    }

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

    const int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;
    TH1D *hmult[numofmultbins + 1];
    TH1D *hmultClone[numofmultbins + 1];
    TH1D *hGen[numofmultbins + 1];
    TH1D *hRec[numofmultbins + 1];
    int multLow, multHigh;

    for (int i = 0; i < numofmultbins + 1; i++)
    {
        if (i == 0)
        {
            multLow = 0;
            multHigh = 100;
        }
        else
        {
            multLow = (int)mult_classes[i - 1];
            multHigh = (int)mult_classes[i];
        }

        hmult[i] = (TH1D *)fspectra->Get(Form("mult_%d-%d/corrected_spectra_Integral_final", multLow, multHigh));
        hGen[i] = (TH1D *)fspectra->Get(Form("mult_%d-%d/Generated", multLow, multHigh));
        hRec[i] = (TH1D *)fspectra->Get(Form("mult_%d-%d/Reconstructed", multLow, multHigh));

        if (hmult[i] == nullptr || hGen[i] == nullptr || hRec[i] == nullptr)
        {
            cout << "Histogram not found in path: " << Form("mult_%d-%d/corrected_spectra_Integral_final", multLow, multHigh) << endl;
            return;
        }
        hmultClone[i] = (TH1D *)hmult[i]->Clone(Form("hmultClone%d", i));

        TH1D *h1 = (TH1D *)hmultClone[i]->Clone("h1");
        TH1D *h2 = (TH1D *)hmultClone[i]->Clone("h2");

        for (int i = 1; i <= h1->GetNbinsX(); i++) // putting small systematic error by hand
        {
            double systemerr = (0.1 * h1->GetBinContent(i));
            h1->SetBinError(i, systemerr);
        }

        Double_t min = 0.0;
        Double_t max = 20.0;
        Double_t loprecision = 0.01;
        Double_t hiprecision = 0.5;
        Option_t *opt = "RI+";
        TString logfilename = reweightEffPath + "/log_mean.root";
        Double_t minfit = 0.0;
        Double_t maxfit = 20.0;

        TF1 *fitFcn = new TF1(Form("fitfunc_%d", i), FuncLavy, 0.0, 20.0, 4);
        fitFcn->SetParameter(0, 5.0);
        fitFcn->SetParameter(1, 0.5);
        fitFcn->FixParameter(2, 0.895);
        fitFcn->SetParameter(3, 0.35);
        fitFcn->SetParNames("n", "dn/dy", "mass", "T");

        TH1 *hout = YieldMean(h1, h2, fitFcn, min, max, loprecision, hiprecision, opt, logfilename, minfit, maxfit);

        TFile *fReweight = new TFile(Form("%s/ReweightFactor_mult_%d-%d.root", reweightEffPath.Data(), multLow, multHigh), "RECREATE");
        int someFactor = ReweightEfficiency(hmult[i], fitFcn, hGen[i], hRec[i], fReweight);
    }
}