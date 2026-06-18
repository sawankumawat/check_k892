// #include "../src/common_glue.h"
// #include "../src/fitting_range_glue.h"
#include "../src/style.h"

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);
Double_t single_BW(double *x, double *par);
Double_t BWsumMassDepWidth(double *x, double *par);
Double_t BWsumMassDepWidth_exponential(double *x, double *par);
Double_t single_BW_mass_dep_spin0(double *x, double *par);
Double_t single_BW_mass_dep_spin2(double *x, double *par);
Double_t exponential_bkg_3(double *x, double *par);
Double_t BWsum_expol_chkstar(double *x, double *par);
Double_t expol_chkstar(double *x, double *par);
Double_t BWsumMassDepWidth_simple_exponential(double *x, double *par);
Double_t simple_exponential(double *x, double *par);
Double_t BWsum_ConstWidth_ModBolt(double *x, double *par);
Double_t BWsum_constWidth(double *x, double *par);
Double_t coherent_sum(double *x, double *par);
Double_t CoherentSum_modifiedBoltzmann(double *x, double *par);

int colors[] = {kGreen + 4, 28, kMagenta, kBlue};
double masses[] = {f1270Mass, a1320Mass, f1525Mass, f1710Mass};
double widths[] = {f1270Width, a1320Width, f1525Width, f1710Width};
string resonance_names[] = {"f_{2}(1270)", "a_{2}(1320)^{0}", "f'_{2}(1525)", "f_{0}(1710)"};
string resonance_mass[] = {"1270", "1320", "1525", "1710"};

void glueball_fit_flow()
{
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);
    TString inputFolder = "../../output/glueball/Flow/";
    TString file = "699196";
    TFile *fInputFile = new TFile((inputFolder + file + "/glueball_flow.root"), "Read");
    if (fInputFile->IsZombie())
    {
        cerr << "File not found " << endl;
        return;
    }
    TString outputFolder = inputFolder + file + "/Fit/";
    if (gSystem->mkdir(outputFolder, kTRUE))
    {
        std::cout << "Folder " << outputFolder << " created successfully." << std::endl;
    }

    TH1F *hMultiplicity = (TH1F *)fInputFile->Get("Multiplicity");
    if (hMultiplicity == nullptr)
    {
        cout << "Multiplicity histogram not found" << endl;
        return;
    }
    int totalEvents = hMultiplicity->Integral(hMultiplicity->GetXaxis()->FindBin(0.0 + 1e-6), hMultiplicity->GetXaxis()->FindBin(30.0 - 1e-6));
    cout << "Total events in multiplicity range 0-30%: " << totalEvents << endl;

    int totalFlowBins = 6;
    TH1F *hSignal[totalFlowBins];
    TH1F *hRaw[totalFlowBins];

    double f1525_3SigmaP = f1525Mass + 3 * f1525Width;
    double f1525_3SigmaM = f1525Mass - 3 * f1525Width;

    int RebinFactors[] = {3, 4, 5, 4, 4, 4};
    // float fitRangeLow[] = {1.25, 1.25, 1.25, 1.25, 1.25, 1.25};
    // float fitRangeHigh[] = {2.2, 2.2, 2.2, 2.2, 2.2, 2.2};

    double fitRangeLow[] = {f1525_3SigmaM, f1525_3SigmaM, f1525_3SigmaM, f1525_3SigmaM, f1525_3SigmaM, f1525_3SigmaM};
    double fitRangeHigh[] = {f1525_3SigmaP, f1525_3SigmaP, f1525_3SigmaP, f1525_3SigmaP, f1525_3SigmaP, f1525_3SigmaP};

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.SetTextFont(42);

    TH1F *hYield1525 = new TH1F("hYield1525", "Yield of f'_{2}(1525) vs Event Plane Bins", totalFlowBins, 0, totalFlowBins);
    double binwidthfile = 0.0;

    for (int ibins = 0; ibins < totalFlowBins; ibins++)
    {
        // Perform fitting for each event plane bin
        hSignal[ibins] = (TH1F *)fInputFile->Get(Form("hInvMassSubtracted_%d", ibins));
        hRaw[ibins] = (TH1F *)fInputFile->Get(Form("hInvMass_%d", ibins));
        if (hSignal[ibins] == nullptr || hRaw[ibins] == nullptr)
        {
            cout << "Invariant mass histogram for event plane bin " << ibins << " not found" << endl;
            return;
        }

        hSignal[ibins]->Rebin(RebinFactors[ibins]);
        double binwidthfile = hSignal[ibins]->GetXaxis()->GetBinWidth(1);

        // =========Making a canvas to show the signal and fit===========
        TCanvas *cFit = new TCanvas(Form("cFit_%d", ibins), "", 720, 720);
        SetCanvasStyle(cFit, 0.15, 0.03, 0.05, 0.15);
        SetHistoQA(hSignal[ibins]);
        hSignal[ibins]->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
        hSignal[ibins]->GetXaxis()->SetRangeUser(1.0, 2.5);
        hSignal[ibins]->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidthfile * 1000));
        hSignal[ibins]->SetMaximum(hSignal[ibins]->GetMaximum() * 1.4);
        hSignal[ibins]->Draw("E");
        lat.DrawLatex(0.35, 0.9, Form("EP bin %d", ibins + 1));
        lat.DrawLatex(0.35, 0.86, Form("Multiplicity: %d-%d%%", 0, 30));
        lat.DrawLatex(0.35, 0.83, Form("p_{T} range: %d-%d GeV/c", 2, 30));

        float fitlow = fitRangeLow[ibins];
        float fithigh = fitRangeHigh[ibins];

        //======================Default fit function==============================
        TF1 *BEexpol = new TF1("BEexpol", BWsumMassDepWidth_exponential, fitlow, fithigh, 16); // expol 3

        string parnames[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma", "a", "b", "c", "d"};
        for (int i = 0; i < sizeof(parnames) / sizeof(parnames[0]); i++)
        {
            BEexpol->SetParName(i, parnames[i].c_str());
        }

        double parameters[] = {0, f1270Mass, f1270Width, 0, a1320Mass, a1320Width, 450, f1525Mass, f1525Width, 0, f1710Mass, f1710Width};

        int size_fitparams = sizeof(parameters) / sizeof(parameters[0]);

        for (int i = 0; i < size_fitparams; i++)
        {
            BEexpol->SetParameter(i, parameters[i]);
        }

        // //********systematic studies*************
        double initial_param_bkg[] = {1e4, 0.5, 7.0, 1.8}; // default fit

        vector<vector<double>> par_limits = {{1, 1 * f1270Width}, {2, 5 * f1270WidthErr}, {4, 1 * a1320Width}, {5, 5 * a1320WidthErr}, {7, 1 * f1525Width}, {8, 1 * f1525WidthErr}, {10, 1 * f1710Width}, {11, 5 * f1710WidthErr}};

        int limits_size = par_limits.size();
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            BEexpol->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        // BEexpol->SetParLimits(0, 0, 1e6);
        // BEexpol->SetParLimits(3, 0, 1e6);
        BEexpol->SetParLimits(6, 0, 1e6);
        BEexpol->SetParLimits(9, 0, 1e6);

        BEexpol->FixParameter(0, 0); // Ignoring the f2(1270) and a2(1320) contribution for now.
        BEexpol->FixParameter(3, 0);

        BEexpol->FixParameter(1, f1270Mass);
        BEexpol->FixParameter(2, f1270Width);
        BEexpol->FixParameter(5, a1320Width);
        BEexpol->FixParameter(4, a1320Mass);

        BEexpol->SetParameter(7, f1525Mass);
        BEexpol->SetParameter(8, f1525Width);

        // BEexpol->SetParameter(10, f1710Mass);
        // BEexpol->SetParameter(11, f1710Width);

        BEexpol->FixParameter(9, 0); // Ignoring the f0(1710) contribution for now.
        BEexpol->FixParameter(10, f1710Mass);
        BEexpol->FixParameter(11, f1710Width);

        TFitResultPtr fitResultptr;
        fitResultptr = hSignal[ibins]->Fit("BEexpol", "REBMS");

        double *obtained_parameters = BEexpol->GetParameters();

        //========================For default fit====================================
        TF1 *expol = new TF1("expol", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);
        TF1 *expol_clone = new TF1("expol_clone", exponential_bkg_3, BEexpol->GetXmin(), BEexpol->GetXmax(), 4);

        for (int i = 0; i < 4; i++)
        {
            expol->SetParameter(i, obtained_parameters[size_fitparams + i]);
            expol_clone->SetParameter(i, obtained_parameters[size_fitparams + i]);
        }
        expol->SetLineColor(3);
        expol->SetLineStyle(2);
        expol_clone->SetLineColor(3);
        expol_clone->SetLineStyle(2);
        expol->Draw("same");

        // //=====================Default fit=======================================
        TF1 *onlyBW = new TF1("onlyBW", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);
        TF1 *onlyBW_clone = new TF1("onlyBW_clone", BWsumMassDepWidth, BEexpol->GetXmin(), BEexpol->GetXmax(), 12);

        string parameter_names[] = {"f_{2}(1270) Amp", "f_{2}(1270) Mass", "f_{2}(1270) #Gamma", "a_{2}(1320)^{0} Amp", "a_{2}(1320)^{0} Mass", "a_{2}(1320)^{0} #Gamma", "f'_{2}(1525) Amp", "f'_{2}(1525) Mass", "f'_{2}(1525) #Gamma", "f_{0}(1710) Amp", "f_{0}(1710) Mass", "f_{0}(1710) #Gamma"};
        for (int i = 0; i < 12; i++)
        {
            onlyBW->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParameter(i, obtained_parameters[i]);
            onlyBW_clone->SetParName(i, parameter_names[i].c_str());
        }
        for (int i = 0; i < limits_size; i++)
        {
            int param_index = static_cast<int>(par_limits[i][0]); // Cast the first element to int
            onlyBW->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
            onlyBW_clone->SetParLimits(par_limits[i][0], parameters[param_index] - par_limits[i][1], parameters[param_index] + par_limits[i][1]);
        }

        // onlyBW_clone->SetParLimits(0, 0.0, 1e6); // norm1270
        // onlyBW_clone->SetParLimits(3, 0.0, 1e6); // norm1320
        onlyBW_clone->SetParLimits(6, 0.0, 1e6); // norm1525
        onlyBW_clone->SetParLimits(9, 0.0, 1e6); // norm1710

        onlyBW_clone->FixParameter(0, 0); // Ignoring the f2(1270) and a2(1320) contribution for now.
        onlyBW_clone->FixParameter(3, 0);

        onlyBW_clone->FixParameter(1, f1270Mass);
        onlyBW_clone->FixParameter(2, f1270Width);
        onlyBW_clone->FixParameter(4, a1320Mass);
        onlyBW_clone->FixParameter(5, a1320Width);

        onlyBW_clone->SetParameter(7, f1525Mass);
        onlyBW_clone->SetParameter(8, f1525Width);

        onlyBW_clone->SetParameter(10, f1710Mass);
        onlyBW_clone->SetParameter(11, f1710Width);

        onlyBW->SetLineColor(4);
        onlyBW->SetLineStyle(2);
        // onlyBW->Draw("same");

        // // Now plot the indivial resonances
        TF1 *singlefits[4];
        for (int i = 2; i < 3; i++)
        {
            singlefits[i] = (i < 3) ? new TF1(Form("singlef%d", i), single_BW_mass_dep_spin2, 1.00, 3.0, 3) : new TF1(Form("singlef%d", i), single_BW_mass_dep_spin0, 1.00, 3.0, 3); // Default
            singlefits[i]->SetParameter(0, obtained_parameters[3 * i]);
            singlefits[i]->SetParameter(1, obtained_parameters[3 * i + 1]);
            singlefits[i]->SetParameter(2, obtained_parameters[3 * i + 2]);
            singlefits[i]->SetLineColor(colors[i]);
            singlefits[i]->SetLineStyle(2);
            // if (i == 3)
            //     singlefits[i]->SetLineWidth(3);
            singlefits[i]->Draw("same");
        }

        TLegend *ltemp = new TLegend(0.3, 0.65, 0.55, 0.9);
        ltemp->SetFillStyle(0);
        ltemp->SetBorderSize(0);
        ltemp->SetTextFont(42);
        ltemp->SetTextSize(0.028);
        ltemp->AddEntry((TObject *)0, "", "");
        ltemp->AddEntry((TObject *)0, "", "");
        ltemp->AddEntry(hSignal[ibins], "Data (stat. uncert.)", "lpe");
        ltemp->AddEntry(BEexpol, "4rBW + Residual BG", "l");
        ltemp->AddEntry(expol, "Residual BG", "l");
        // ltemp->AddEntry(singlefits[0], "f_{2}(1270)", "l");
        // ltemp->AddEntry(singlefits[1], "a_{2}(1320)^{0}", "l");
        ltemp->AddEntry(singlefits[2], "f'_{2}(1525)", "l");
        // ltemp->AddEntry(singlefits[3], "f_{0}(1710)", "l");
        ltemp->Draw("same");

        for (int i = 2; i < 3; i++)
        {
            auto lowLimit = (masses[i] - 2 * widths[i] > 1.0) ? masses[i] - 2 * widths[i] : 1.0;
            auto highLimit = (masses[i] + 2 * widths[i] < 3.0) ? masses[i] + 2 * widths[i] : 3.0;
            double significance_num = singlefits[i]->Integral(lowLimit, highLimit) / binwidthfile;
            int binlow = hRaw[ibins]->GetXaxis()->FindBin(lowLimit);
            int binhigh = hRaw[ibins]->GetXaxis()->FindBin(highLimit);
            double significance_den = TMath::Sqrt(hRaw[ibins]->Integral(binlow, binhigh));
            double significance = significance_num / significance_den;
            double signal_counts = significance_num;
            double background_counts = hRaw[ibins]->Integral(binlow, binhigh) - signal_counts;

            cout << "numerator " << significance_num << " denominator " << significance_den << endl;
            cout << "Significance of " << resonance_names[i] << " is " << significance << endl;
            double amplitude = obtained_parameters[3 * i];
            double amplitude_err = BEexpol->GetParError(3 * i);
            double statSignificance = amplitude / amplitude_err;
            cout << "Statistical significance of " << resonance_names[i] << " is " << statSignificance << endl;
        }
        // BEexpol->Draw("same");
        cFit->SaveAs((outputFolder + Form("cFit_%d.png", ibins)).Data());

        //=====Yield from function integration method=========
        double nsigma_yield = 3.0;
        double fitRegion1525_low = f1525Mass - nsigma_yield * f1525Width;
        double fitRegion1525_high = f1525Mass + nsigma_yield * f1525Width;
        double yield1525 = singlefits[2]->Integral(fitRegion1525_low, fitRegion1525_high) / (binwidthfile * totalEvents);

        // Create covariance matrices
        TMatrixDSym cov = fitResultptr->GetCovarianceMatrix();
        TMatrixDSym cov1270(3), cov1320(3), cov1525(3), cov1710(3);
        // cov.GetSub(0, 2, 0, 2, cov1270);   // f1270: parameters 0-2
        // cov.GetSub(3, 5, 3, 5, cov1320);   // a1320: parameters 3-5
        cov.GetSub(6, 8, 6, 8, cov1525); // f1525: parameters 6-8
        // cov.GetSub(9, 11, 9, 11, cov1710); // f1710: parameters 9-11

        // Covariance matrices for combined fits
        Double_t *cov1525_array = cov1525.GetMatrixArray();
        Double_t *para1525 = singlefits[2]->GetParameters();

        // Errors propagation for yields using the covariance matrix
        double yield1525_err = singlefits[2]->IntegralError(fitRegion1525_low, fitRegion1525_high, para1525, cov1525_array) / (binwidthfile * totalEvents);

        hYield1525->SetBinContent(ibins + 1, yield1525);
        hYield1525->SetBinError(ibins + 1, yield1525_err);
        cout << "Yield of f'_{2}(1525) in EP bin " << ibins + 1 << " is " << yield1525 << " ± " << yield1525_err << endl;
    }

    TCanvas *cYield = new TCanvas("cYield", "", 720, 720);
    SetCanvasStyle(cYield, 0.15, 0.03, 0.05, 0.15);
    SetHistoQA(hYield1525);
    hYield1525->GetXaxis()->SetTitle("Event Plane Bins");
    hYield1525->GetYaxis()->SetTitle("(1/N) dN/dy");
    hYield1525->Draw("E");
    cYield->SaveAs((outputFolder + "cYield.png").Data());
}

// end of main program

// *****************************************************************************************************
//************************************Fit functions************************************************* */

// We will define the single BW and the sum of 4 BWs
Double_t single_BW_hera(double *x, double *par)
{
    // normalization factor is missing and how to add it I am not sure
    double amplitude = par[0];
    double mass = par[1];
    double width = par[2];

    double den = (x[0] * x[0] - mass * mass) * (x[0] * x[0] - mass * mass) + mass * mass * width * width;
    double realnum = (mass * mass - x[0] * x[0]) * mass * TMath::Sqrt(width);
    double imagnum = mass * mass * width * TMath::Sqrt(width);

    double real3BW = realnum / den;
    double imag3BW = imagnum / den;
    double sig1 = amplitude * (real3BW * real3BW + imag3BW * imag3BW);

    return sig1;
}

Double_t single_BW_mass_dep_spin0(double *x, double *par)
{
    double num = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double den = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    int spin = 0;
    double n1 = (2.0 * spin + 1.0) / 2.0;

    double yield = par[0];
    double mass = par[1];
    double width = par[2] * (TMath::Power(mass / x[0], 1.0)) * TMath::Power((num) / (den), n1);

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));

    return fit;
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

Double_t BWsum_hera(double *x, double *par) // taken 4 resonances here
{
    // total 10 parameters, 2 for each resonance (total 4 resonance), 1 for normalization of spin2 particles, 1 for normalization of spin0 particle
    double mass1270 = par[0];
    double width1270 = par[1];
    double mass1320 = par[2];
    double width1320 = par[3];
    double mass1525 = par[4];
    double width1525 = par[5];
    double mass1710 = par[6];
    double width1710 = par[7];
    double a0 = par[8];
    double a1 = par[9];

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1320 = (x[0] * x[0] - mass1320 * mass1320) * (x[0] * x[0] - mass1320 * mass1320) + mass1320 * mass1320 * width1320 * width1320;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270);
    double realnum1320 = (mass1320 * mass1320 - x[0] * x[0]) * mass1320 * TMath::Sqrt(width1320);
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525);
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710);

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270);
    double imagnum1320 = mass1320 * mass1320 * width1320 * TMath::Sqrt(width1320);
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525);
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710);

    double real3BW = 5 * realnum1270 / den1270 - 3 * realnum1320 / den1320 + 2 * realnum1525 / den1525;

    double imag3BW = 5 * imagnum1270 / den1270 - 3 * imagnum1320 / den1320 + 2 * imagnum1525 / den1525;

    double sig1 = (real3BW * real3BW + imag3BW * imag3BW);

    double fit = a0 * sig1 + a1 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / (den1710 * den1710);

    return fit;
}

Double_t BWsum_hera_mass_dep(double *x, double *par) // taken 4 resonances here
{
    // total 10 parameters, 2 for each resonance (total 4 resonance), 1 for normalization of spin2 particles, 1 for normalization of spin0 particle
    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[0] * par[0] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[2] * par[2] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart4 = par[6] * par[6] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    double mass1270 = par[0];
    double width1270 = par[1] * (TMath::Power(par[0] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double mass1320 = par[2];
    double width1320 = par[3] * (TMath::Power(par[2] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double mass1525 = par[4];
    double width1525 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n1);
    double mass1710 = par[6];
    double width1710 = par[7] * (TMath::Power(par[6] / x[0], 1.0)) * TMath::Power((npart1) / (dpart4), n2);

    // double mass1270 = par[0];
    // double width1270 = par[1];
    // double mass1320 = par[2];
    // double width1320 = par[3];
    // double mass1525 = par[4];
    // double width1525 = par[5];
    // double mass1710 = par[6];
    // double width1710 = par[7];
    double a0 = par[8];
    double a1 = par[9];
    // double norm1 = par[10];
    // double norm2 = par[11];
    // double norm3 = par[12];
    double norm1 = 5;
    double norm2 = -3;
    double norm3 = 2;

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1320 = (x[0] * x[0] - mass1320 * mass1320) * (x[0] * x[0] - mass1320 * mass1320) + mass1320 * mass1320 * width1320 * width1320;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270);
    double realnum1320 = (mass1320 * mass1320 - x[0] * x[0]) * mass1320 * TMath::Sqrt(width1320);
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525);
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710);

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270);
    double imagnum1320 = mass1320 * mass1320 * width1320 * TMath::Sqrt(width1320);
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525);
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710);

    double real3BW = norm1 * realnum1270 / den1270 + norm2 * realnum1320 / den1320 + norm3 * realnum1525 / den1525;

    double imag3BW = norm1 * imagnum1270 / den1270 + norm2 * imagnum1320 / den1320 + norm3 * imagnum1525 / den1525;

    double sig1 = (real3BW * real3BW + imag3BW * imag3BW);

    double fit = a0 * sig1 + a1 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710) / (den1710 * den1710);

    return fit;
}

Double_t BWsum_hera_const(double *x, double *par) // taken 4 resonances here
{
    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
    double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    double yield1270 = par[0];
    double mass1270 = par[1];
    // double width1270 = par[2];
    double width1270 = par[2] * (TMath::Power(par[1] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double yield1320 = par[3];
    double mass1320 = par[4];
    // double width1320 = par[5];
    double width1320 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double yield1525 = par[6];
    double mass1525 = par[7];
    // double width1525 = par[8];
    double width1525 = par[8] * (TMath::Power(par[7] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n1);
    double yield1710 = par[9];
    double mass1710 = par[10];
    // double width1710 = par[11];
    double width1710 = par[11] * (TMath::Power(par[10] / x[0], 1.0)) * TMath::Power((npart1) / (dpart4), n2);

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1320 = (x[0] * x[0] - mass1320 * mass1320) * (x[0] * x[0] - mass1320 * mass1320) + mass1320 * mass1320 * width1320 * width1320;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270) / den1270;
    double realnum1320 = (mass1320 * mass1320 - x[0] * x[0]) * mass1320 * TMath::Sqrt(width1320) / den1320;
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525) / den1525;
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710) / den1710;

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270) / den1270;
    double imagnum1320 = mass1320 * mass1320 * width1320 * TMath::Sqrt(width1320) / den1320;
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525) / den1525;
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710) / den1710;

    double fit1270 = yield1270 * (realnum1270 * realnum1270 + imagnum1270 * imagnum1270);
    double fit1320 = yield1320 * (realnum1320 * realnum1320 + imagnum1320 * imagnum1320);
    double fit1525 = yield1525 * (realnum1525 * realnum1525 + imagnum1525 * imagnum1525);
    double fit1710 = yield1710 * (realnum1710 * realnum1710 + imagnum1710 * imagnum1710);
    double fit = fit1270 + fit1320 + fit1525 + fit1710;

    return fit;
}

Double_t coherent_sum(double *x, double *par) // taken 4 resonances here
{
    // Fit is |a1 *BW1 + a2*BW2 e^{i*phi1} + a3*BW3 e^{i*phi2}|^2 + |a4 *BW4|^2
    // Then the real and imaginary parts are separated
    // Total 14 parameters. 3 for each resonance and 2 for phases.

    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
    double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    // double temp = par[3];
    // double temp2 = par[6];

    double norm1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2] * (TMath::Power(par[1] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double norm1320 = par[3];
    // double norm1320 = par[0] * 3 / 5;
    double mass1320 = par[4];
    double width1320 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double norm1525 = par[6];
    // double norm1525 = par[0] * 2 / 5;
    double mass1525 = par[7];
    double width1525 = par[8] * (TMath::Power(par[7] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n1);
    double norm1710 = par[9];
    double mass1710 = par[10];
    double width1710 = par[11] * (TMath::Power(par[10] / x[0], 1.0)) * TMath::Power((npart1) / (dpart4), n2);

    double den1270 = (x[0] * x[0] - mass1270 * mass1270) * (x[0] * x[0] - mass1270 * mass1270) + mass1270 * mass1270 * width1270 * width1270;
    double den1320 = (x[0] * x[0] - mass1320 * mass1320) * (x[0] * x[0] - mass1320 * mass1320) + mass1320 * mass1320 * width1320 * width1320;
    double den1525 = (x[0] * x[0] - mass1525 * mass1525) * (x[0] * x[0] - mass1525 * mass1525) + mass1525 * mass1525 * width1525 * width1525;
    double den1710 = (x[0] * x[0] - mass1710 * mass1710) * (x[0] * x[0] - mass1710 * mass1710) + mass1710 * mass1710 * width1710 * width1710;

    double realnum1270 = (mass1270 * mass1270 - x[0] * x[0]) * mass1270 * TMath::Sqrt(width1270) / den1270;
    double realnum1320 = (mass1320 * mass1320 - x[0] * x[0]) * mass1320 * TMath::Sqrt(width1320) / den1320;
    double realnum1525 = (mass1525 * mass1525 - x[0] * x[0]) * mass1525 * TMath::Sqrt(width1525) / den1525;
    double realnum1710 = (mass1710 * mass1710 - x[0] * x[0]) * mass1710 * TMath::Sqrt(width1710) / den1710;

    double imagnum1270 = mass1270 * mass1270 * width1270 * TMath::Sqrt(width1270) / den1270;
    double imagnum1320 = mass1320 * mass1320 * width1320 * TMath::Sqrt(width1320) / den1320;
    double imagnum1525 = mass1525 * mass1525 * width1525 * TMath::Sqrt(width1525) / den1525;
    double imagnum1710 = mass1710 * mass1710 * width1710 * TMath::Sqrt(width1710) / den1710;

    double phase1 = par[12]; // this is angle phi in radians
    double phase2 = par[13];

    double real1 = realnum1270;
    double imag1 = imagnum1270;
    double real2 = realnum1320 * TMath::Cos(phase1) - imagnum1320 * TMath::Sin(phase1);
    double imag2 = realnum1320 * TMath::Sin(phase1) + imagnum1320 * TMath::Cos(phase1);
    double real3 = realnum1525 * TMath::Cos(phase2) - imagnum1525 * TMath::Sin(phase2);
    double imag3 = realnum1525 * TMath::Sin(phase2) + imagnum1525 * TMath::Cos(phase2);
    double real4 = realnum1710;
    double imag4 = imagnum1710;

    double real_sum = norm1270 * real1 + norm1320 * real2 + norm1525 * real3;
    double imag_sum = norm1270 * imag1 + norm1320 * imag2 + norm1525 * imag3;
    double fit = norm1710 * norm1710 * (real4 * real4 + imag4 * imag4) + (real_sum * real_sum + imag_sum * imag_sum);
    // double fit = norm1270 * (real1 * real1 + imag1 * imag1) + norm1320 * (real2 * real2 + imag2 * imag2) + norm1525 * (real3 * real3 + imag3 * imag3) + norm1710 * (real4 * real4 + imag4 * imag4); // independent sum (for cross check)

    return fit;
}

Double_t single_BW(double *x, double *par)
{
    double yield = par[0];
    double mass = par[1];
    double width = par[2];

    double fit = yield * mass * width * x[0] / (pow((x[0] * x[0] - mass * mass), 2) + pow(mass * width, 2));
    return fit;
}

Double_t BWsum_constWidth(double *x, double *par)
{
    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2];
    double yield1320 = par[3];
    double mass1320 = par[4];
    double width1320 = par[5];
    double yield1525 = par[6];
    double mass1525 = par[7];
    double width1525 = par[8];
    double yield1710 = par[9];
    double mass1710 = par[10];
    double width1710 = par[11];

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1320 = yield1320 * mass1320 * width1320 * x[0] / (pow((x[0] * x[0] - mass1320 * mass1320), 2) + pow(mass1320 * width1320, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    double fit = fit1270 + fit1320 + fit1525 + fit1710;
    return fit;
}

Double_t BWsumMassDepWidth(double *x, double *par)
{
    double npart1 = x[0] * x[0] - 4 * (0.4976 * 0.4976);
    double dpart1 = par[1] * par[1] - 4 * (0.4976 * 0.4976);
    double dpart2 = par[4] * par[4] - 4 * (0.4976 * 0.4976);
    double dpart3 = par[7] * par[7] - 4 * (0.4976 * 0.4976);
    double dpart4 = par[10] * par[10] - 4 * (0.4976 * 0.4976);

    Int_t j1 = 2;
    Int_t j2 = 0;
    double n1 = (2.0 * j1 + 1.0) / 2.0;
    double n2 = (2.0 * j2 + 1.0) / 2.0;

    // double phase_space = (x[0] / TMath::Sqrt(x[0] * x[0] + 15 * 15)) * (TMath::Exp(-TMath::Sqrt(x[0] * x[0] + 15 * 15) / 0.16)); // 160 MeV is the kinetic freeze-out temperature

    double yield1270 = par[0];
    double mass1270 = par[1];
    double width1270 = par[2] * (TMath::Power(par[1] / x[0], 1.0)) * TMath::Power((npart1) / (dpart1), n1);
    double yield1320 = par[3];
    double mass1320 = par[4];
    double width1320 = par[5] * (TMath::Power(par[4] / x[0], 1.0)) * TMath::Power((npart1) / (dpart2), n1);
    double yield1525 = par[6];
    double mass1525 = par[7];
    double width1525 = par[8] * (TMath::Power(par[7] / x[0], 1.0)) * TMath::Power((npart1) / (dpart3), n1);
    double yield1710 = par[9];
    double mass1710 = par[10];
    double width1710 = par[11] * (TMath::Power(par[10] / x[0], 1.0)) * TMath::Power((npart1) / (dpart4), n2);

    double fit1270 = yield1270 * mass1270 * width1270 * x[0] / (pow((x[0] * x[0] - mass1270 * mass1270), 2) + pow(mass1270 * width1270, 2));
    double fit1320 = yield1320 * mass1320 * width1320 * x[0] / (pow((x[0] * x[0] - mass1320 * mass1320), 2) + pow(mass1320 * width1320, 2));
    double fit1525 = yield1525 * mass1525 * width1525 * x[0] / (pow((x[0] * x[0] - mass1525 * mass1525), 2) + pow(mass1525 * width1525, 2));
    double fit1710 = yield1710 * mass1710 * width1710 * x[0] / (pow((x[0] * x[0] - mass1710 * mass1710), 2) + pow(mass1710 * width1710, 2));

    double fit = (fit1270 + fit1320 + fit1525 + fit1710);
    return fit;
}

// Now we will define the functions for the exponential background

Double_t simple_exponential(double *x, double *par) // 2 parameters
{
    return (par[0] * pow(x[0], par[1]) * TMath::Exp(-par[2] * x[0]));
}

Double_t exponential_bkg_1(double *x, double *par) // 3 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * (x[0] - 2.0 * 0.497)));
}

Double_t exponential_bkg_2(double *x, double *par) // 3 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[1])));
}

Double_t exponential_bkg_3(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * pow((x[0] - 2.0 * 0.497), par[3])));
}

Double_t exponential_bkg_4(double *x, double *par) // 5 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * ((pow((x[0] - 2.0 * 0.497), par[3])) + pow((x[0] - 2.0 * 0.497), par[4]))));
}

Double_t exponential_bkg_5(double *x, double *par) // 3 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(-par[2] * x[0]));
}

Double_t exponential_bkg_6(double *x, double *par) // 4 parameters
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] + par[1] * x[0] + par[3] * x[0] * x[0]));
}

// Now we will define the functions for the Boltzmann background
Double_t Boltzmann_bkg_1(double *x, double *par) // 3 parameters
{
    double norm = par[0];
    double massks = 0.497;
    double n = par[1];
    double C = par[2];
    double y = norm * pow((x[0] - 2 * massks), n / 2) * pow(C, 3 / 2) * exp(-C * pow((x[0] - 2 * massks), n));
    return y;
}

Double_t Boltzmann_bkg_2(double *x, double *par) // 4 parameters
{
    double norm = par[0];
    double massks = 0.497;
    double n = par[1];
    double C = par[2];
    double n1 = par[3];
    double y = norm * pow((x[0] - 2 * massks), n / 2) * pow(C, 3 / 2) * exp(-C * pow((x[0] - 2 * massks), n1));
    return y;
}

// Now we will define the functions for the sum of the Breit-Wigner and the exponential background

Double_t single_BW_expol3(double *x, double *par)
{
    return (single_BW(x, par) + exponential_bkg_3(x, &par[3]));
}

Double_t single_BW_expol3_hera(double *x, double *par)
{
    return (single_BW_hera(x, par) + exponential_bkg_3(x, &par[3]));
}

Double_t BWsum_ConstWidth_ModBolt(double *x, double *par)
{
    return (BWsum_constWidth(x, par) + exponential_bkg_3(x, &par[12]));
}

Double_t BWsum_expol3_hera(double *x, double *par)
{
    return (BWsum_hera(x, par) + exponential_bkg_3(x, &par[12]));
}

// Now we will define the functions for the sum of the Breit-Wigner and the Boltzmann background

Double_t single_BW_boltzman_1(double *x, double *par)
{
    return (single_BW(x, par) + Boltzmann_bkg_1(x, &par[3]));
}

Double_t single_BW_boltzman_2(double *x, double *par)
{
    return (single_BW(x, par) + Boltzmann_bkg_2(x, &par[3]));
}

Double_t BWsum_boltzman_1(double *x, double *par)
{
    return (BWsumMassDepWidth(x, par) + Boltzmann_bkg_1(x, &par[12]));
}

Double_t BWsum_boltzman_2(double *x, double *par)
{
    return (BWsum_constWidth(x, par) + Boltzmann_bkg_2(x, &par[12]));
}
Double_t expol_chkstar(double *x, double *par)
{
    return (par[0] * pow((x[0] - 2.0 * 0.497), par[1]) * TMath::Exp(par[2] * x[0] + x[0] * x[0] * par[3]));
}
Double_t BWsum_expol_chkstar(double *x, double *par)
{
    return (BWsumMassDepWidth(x, par) + expol_chkstar(x, &par[12]));
}

Double_t BWsumMassDepWidth_exponential(double *x, double *par)
{
    // return (BWsumMassDepWidth(x, par) + expol_chkstar(x, &par[12]));
    return (BWsumMassDepWidth(x, par) + exponential_bkg_3(x, &par[12]));
}

Double_t BWsumMassDepWidth_simple_exponential(double *x, double *par)
{
    return (BWsumMassDepWidth(x, par) + simple_exponential(x, &par[12]));
}

Double_t BWsum_modifiedBoltzmann_hera(double *x, double *par)
{
    return (BWsum_hera(x, par) + exponential_bkg_3(x, &par[10]));
}
Double_t BWsum_ModifiedBoltzmann_hera_mass_dep(double *x, double *par)
{
    return (BWsum_hera_mass_dep(x, par) + exponential_bkg_3(x, &par[10]));
}
Double_t BWsum_modifiedBoltzmann_hera_const(double *x, double *par)
{
    return (BWsum_hera_const(x, par) + exponential_bkg_3(x, &par[12]));
}
Double_t CoherentSum_modifiedBoltzmann(double *x, double *par)
{
    // return (coherent_sum(x, par) + expol_chkstar(x, &par[14]));
    return (coherent_sum(x, par) + exponential_bkg_3(x, &par[14]));
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1); // top pad
    TPad *pad2 = (TPad *)c->GetPad(2); // bottom pad
    pad2Size = 0.5;                    // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.5, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.5);
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.23);
    pad1->SetLeftMargin(0.125);
    pad2->SetLeftMargin(0.125);
    pad1->SetTopMargin(0.1);
    pad1->SetBottomMargin(0);
    pad2->SetTopMargin(0);
}