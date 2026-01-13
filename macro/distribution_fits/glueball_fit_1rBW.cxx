#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

Double_t single_BW_mass_dep_spin2(double *x, double *par);
Double_t pol3func(double *x, double *par);
Double_t single_BW_pol3(double *x, double *par);
Double_t exponential_bkg_3(double *x, double *par);
Double_t single_BW_expol(double *x, double *par);
int colors[] = {kGreen + 4, 28, kMagenta, kBlue, kBlack};

void glueball_fit_1rBW()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/435450/KsKs_Channel/higher-mass-resonances"; // 2023 dataset
    TFile *file = new TFile((path + "/hglue_ROTATED_allPtMult2.root").c_str(), "READ");
    if (file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << "\n";
        return;
    }
    TH1F *hmult = (TH1F *)file->Get("multiplicity_histogram");
    if (hmult == nullptr)
    {
        cout << "Multiplicity histogram not found" << endl;
        return;
    }
    int multlow, multhigh;

    // // Temporary for single multiplicity class checking
    multlow = 0;
    multhigh = 100;

    int b1 = hmult->GetXaxis()->FindBin(multlow + 0.01);
    int b2 = hmult->GetXaxis()->FindBin(multhigh - 0.01);
    double total_events = hmult->Integral(b1, b2);

    double pos_sum = 0.0, neg_sum = 0.0;
    for (int b = b1; b <= b2; ++b)
    {
        double c = hmult->GetBinContent(b);
        if (c >= 0)
            pos_sum += c;
        else
            neg_sum += c;
    }
    // cout << "total entries (all bins): " << hmult->GetEntries()
    //      << ", selected sum: " << total_events
    //      << " (pos: " << pos_sum << ", neg: " << neg_sum << ")"
    //      << ", bins: [" << b1 << ", " << b2 << "]" << endl;

    // TCanvas *cMultDist = new TCanvas("", "Multiplicity Distribution", 720, 720);
    // SetCanvasStyle(cMultDist, 0.14, 0.03, 0.05, 0.13);
    // SetHistoQA(hmult);
    // hmult->Draw();
    TCanvas *cGlueMass = new TCanvas("", "Glueball Invariant Mass", 1880, 1080);
    SetCanvasStyle(cGlueMass, 0.14, 0.03, 0.05, 0.13);
    cGlueMass->Divide(3, 3);
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.05);
    lat.SetTextFont(42);
    ofstream outfile;
    outfile.open("f1525_yield.txt");

    for (int ipt = 0; ipt < Npt; ipt++)
    // for (int ipt = 0; ipt < 1; ipt++)
    {
        float lowpT = pT_bins[ipt];
        float highpT = pT_bins[ipt + 1];

        // // Temporary for single bins checking
        // float lowpT = 8.0;
        // float highpT = 10.0;

        double ptbinwidth = highpT - lowpT;

        TH1F *hinvMass = (TH1F *)file->Get(Form("multiplicity_%d_%d/ksks_subtracted_invmass_pt_%.1f_%.1f", multlow, multhigh, lowpT, highpT));
        if (hinvMass == nullptr)
        {
            cout << "Invariant mass histogram not found for pT bin " << lowpT << " - " << highpT << " GeV/c and multiplicity " << multlow << " - " << multhigh << endl;
            return;
        }
        double binwidthfile = (hinvMass->GetXaxis()->GetXmax() - hinvMass->GetXaxis()->GetXmin()) / hinvMass->GetXaxis()->GetNbins();
        cout << "bin width is " << binwidthfile << endl;
        hinvMass->GetXaxis()->SetRangeUser(1.3, 1.7);
        hinvMass->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
        hinvMass->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidthfile * 1000));

        cGlueMass->cd(ipt + 1);
        SetHistoQA(hinvMass);
        hinvMass->SetMinimum(0);
        hinvMass->SetMaximum(hinvMass->GetMaximum() * 1.5);
        hinvMass->Draw("pe");

        float fitlow = 1.35;
        float fithigh = 1.68;

        TF1 *BWpol3 = new TF1("BWpol3", single_BW_expol, fitlow, fithigh, 7);
        BWpol3->SetParNames("Yield", "Mass", "Width", "p0", "p1", "p2", "p3");
        BWpol3->SetParameter(0, 10000);     // yield
        BWpol3->SetParLimits(0, 0, 1e8);    // yield
        BWpol3->SetParameter(1, f1525Mass); // mass
        BWpol3->SetParLimits(1, f1525Mass - 1 * f1525Width, f1525Mass + 1 * f1525Width);
        BWpol3->FixParameter(2, f1525Width); // width
        // BWpol3->SetParLimits(2, f1525Width - 1 * f1525WidthErr, f1525Width + 1 * f1525WidthErr);
        // float bkgparameters[] = {7.6e5, -4.3e4, -2.6e5, 3.5e4}; // 2-3, 3-4 GeV/c (pol 3)
        // float bkgparameters[] = {1.37518e8, 0.6, 7.071167, 1.04}; //  (expol)
        float bkgparameters[] = {0.1, -2.6, 6.71167, 1.0}; //  (expol)
        BWpol3->SetParameter(3, bkgparameters[0]);
        BWpol3->SetParameter(4, bkgparameters[1]);
        BWpol3->SetParameter(5, bkgparameters[2]);
        BWpol3->SetParameter(6, bkgparameters[3]);
        TFitResultPtr fitResultptr = hinvMass->Fit(BWpol3, "REBS");

        TF1 *pol3 = new TF1("pol3", exponential_bkg_3, 1, 3, 4);
        pol3->SetParameter(0, BWpol3->GetParameter(3));
        pol3->SetParameter(1, BWpol3->GetParameter(4));
        pol3->SetParameter(2, BWpol3->GetParameter(5));
        pol3->SetParameter(3, BWpol3->GetParameter(6));
        pol3->SetLineColor(kBlue);
        pol3->SetLineStyle(2);
        pol3->SetRange(fitlow, fithigh);
        pol3->Draw("same");

        TF1 *BWonly = new TF1("BWonly", single_BW_mass_dep_spin2, 1, 3, 3);
        BWonly->SetParameter(0, BWpol3->GetParameter(0));
        BWonly->SetParameter(1, BWpol3->GetParameter(1));
        BWonly->SetParameter(2, BWpol3->GetParameter(2));
        BWonly->SetLineColor(kGreen + 3);
        BWonly->SetLineStyle(2);
        BWonly->SetRange(fitlow, fithigh);
        BWonly->Draw("same");
        lat.DrawLatex(0.12, 0.8, Form("#it{p}_{T}: %.1f - %.1f GeV/c", lowpT, highpT));

        double yieldIntegral = BWonly->Integral(f1525Mass - 3 * f1525Width, f1525Mass + 3 * f1525Width) / (binwidthfile * total_events * ptbinwidth);
        cout << "binwidthfile " << binwidthfile << ", total_events " << total_events << ", ptbinwidth " << ptbinwidth << endl;
        TMatrixDSym cov = fitResultptr->GetCovarianceMatrix();
        TMatrixDSym cov1525(3);
        cov.GetSub(0, 2, 0, 2, cov1525); // Extract covariance matrix for f2'(1525)
        Double_t *cov_array = cov1525.GetMatrixArray();
        Double_t *para = BWonly->GetParameters();
        double yield_err = BWonly->IntegralError(f1525Mass - 3 * f1525Width, f1525Mass + 3 * f1525Width, para, cov_array) / (ptbinwidth * binwidthfile * total_events);
        cout << "pt bin " << lowpT << " - " << highpT << " GeV/c: Yield = " << yieldIntegral << " +/- " << yield_err << endl;
        outfile << yieldIntegral << " ± " << yield_err << endl;
    }
    outfile.close();
    cGlueMass->SaveAs("f2_fit_1rBW.png");
}

Double_t single_BW_mass_dep_spin2(double *x, double *par)
{
    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    int spin = 2;
    double n1 = (2.0 * spin + 1.0) / 2.0;

    double yield = par[0];
    double mass = par[1];
    double width = par[2] * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));

    return fit;
}

Double_t pol3func(double *x, double *par)
{
    return (par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0]);
}

Double_t single_BW_pol3(double *x, double *par)
{
    return (single_BW_mass_dep_spin2(x, par) + pol3func(x, &par[3]));
}

Double_t exponential_bkg_3(double *x, double *par) // 4 parameters
{
    // return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[3]))); // expol
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] * x[0] + x[0] * x[0] * par[3])); // expol as used in charged kstar
}

Double_t single_BW_expol(double *x, double *par)
{
    return (single_BW_mass_dep_spin2(x, par) + exponential_bkg_3(x, &par[3]));
}