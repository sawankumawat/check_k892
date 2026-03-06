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
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void glueball_fit_1rBW()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances"; // 2023 dataset
    string savepath = path + "/fits/4rBw_fits/pt_dependent/";
    TFile *file = new TFile((path + "/hglue_ROTATED_allPtMult3.root").c_str(), "READ");
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

    // TCanvas *cMultDist = new TCanvas("", "Multiplicity Distribution", 720, 720);
    // SetCanvasStyle(cMultDist, 0.14, 0.03, 0.05, 0.13);
    // SetHistoQA(hmult);
    // hmult->Draw();

    TCanvas *cf2BWfit = new TCanvas("", "Glueball Invariant Mass", 1080, 720);
    SetCanvasStyle(cf2BWfit, 0.14, 0.03, 0.05, 0.13);
    cf2BWfit->Divide(3, 2);
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.05);
    lat.SetTextFont(42);

    // Efficiency and signal loss correction
    float ptBins[] = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0};             // 2022 dataset
    float ptBins2[] = {0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 15.5}; // 2022 dataset
    double Npt2 = sizeof(ptBins2) / sizeof(ptBins2[0]) - 1;

    TH1F *hMass = new TH1F("hMass", "Reconstructed f2(1525) mass", Npt2, ptBins2);
    TH1F *hWidth = new TH1F("hWidth", "Reconstructed f2(1525) width", Npt2, ptBins2);
    TH1F *hRawYield = new TH1F("hRawYield", "Reconstructed f2(1525) raw yield", Npt2, ptBins2);
    TH1F *hCorrectedYield = new TH1F("hCorrectedYield", "Reconstructed f2(1525) corrected yield", Npt2, ptBins2);
    TH1F *hSignificance = new TH1F("hSignificance", "Significance of f2(1525)", Npt2, ptBins2);

    TFile *fSigEventLoss = new TFile("../SignalLossPhiINEL.root", "READ");
    if (fSigEventLoss->IsZombie())
    {
        cerr << "Error opening event/signal loss file" << endl;
    }

    TH1F *hEvbySigLoss = (TH1F *)fSigEventLoss->Get("hSignalLoss");
    int binsLoss = hEvbySigLoss->GetNbinsX();
    if (binsLoss != (sizeof(ptBins) / sizeof(ptBins[0]) - 1))
    {
        cerr << "Mismatch in number of pT bins for event/signal loss" << endl;
        return;
    }

    int nBins = sizeof(ptBins) / sizeof(ptBins[0]) - 1;
    int nBins2 = sizeof(ptBins2) / sizeof(ptBins2[0]) - 1;

    // efficiency file
    TFile *feff = new TFile("/home/sawan/check_k892/mc/LHC24l1/463655.root", "read");
    if (feff->IsZombie())
    {
        cout << "Error opening MC file" << endl;
        return;
    }
    string histpathf2 = "higher-mass-resonances_f21525/hMChists";

    THnSparseD *GenpTf2 = (THnSparseD *)feff->Get(Form("%s/Genf17102", histpathf2.c_str()));
    THnSparseD *recpt1f2 = (THnSparseD *)feff->Get(Form("%s/Recf1710_pt2", histpathf2.c_str()));

    if (GenpTf2 == nullptr || recpt1f2 == nullptr)
    {
        cout << "Error reading MC histogram in path " << histpathf2 << endl;
        return;
    }

    TH1F *hefficiencyf2 = new TH1F("hefficiencyf2", "Efficiency vs pT for f21525", nBins, ptBins);

    TH1D *hgenptf2 = GenpTf2->Projection(1);  // project on pt axis
    TH1D *hrecptf2 = recpt1f2->Projection(1); // project on pt axis
    hgenptf2->GetXaxis()->SetRangeUser(0.0, 15.0);
    hrecptf2->GetXaxis()->SetRangeUser(0.0, 15.0);
    cout << "Number of bins in Generated " << hgenptf2->GetNbinsX() << endl;
    cout << "Number of bins in Reconstructed " << hrecptf2->GetNbinsX() << endl;
    cout << "Bin width of generated " << hgenptf2->GetXaxis()->GetBinWidth(1) << endl;

    for (int i = 0; i < nBins; i++)
    {
        // get bin content accroding to cosTheta bins and error according to bayesian method
        int lowpt = hgenptf2->GetXaxis()->FindBin(ptBins[i] + 0.01);
        int highpt = hgenptf2->GetXaxis()->FindBin(ptBins[i + 1] - 0.01);

        double efficiencyf2 = hgenptf2->Integral(lowpt, highpt);
        double recYieldf2 = hrecptf2->Integral(lowpt, highpt);
        double recYieldErrorf2 = TMath::Sqrt(((recYieldf2 + 1) / (efficiencyf2 + 2)) * ((recYieldf2 + 2) / (efficiencyf2 + 3) - (recYieldf2 + 1) / (efficiencyf2 + 2)));
        if (efficiencyf2 > 0)
        {
            hefficiencyf2->SetBinContent(i + 1, recYieldf2 / efficiencyf2);
            hefficiencyf2->SetBinError(i + 1, recYieldErrorf2);
        }
    }

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
        TH1F *hraw = (TH1F *)file->Get(Form("multiplicity_%d_%d/ksks_invmass_pt_%.1f_%.1f", multlow, multhigh, lowpT, highpT));
        double binwidthfile = (hinvMass->GetXaxis()->GetXmax() - hinvMass->GetXaxis()->GetXmin()) / hinvMass->GetXaxis()->GetNbins();
        cout << "bin width is " << binwidthfile << endl;
        hinvMass->GetXaxis()->SetRangeUser(1.3, 1.7);
        hinvMass->GetXaxis()->SetTitle("#it{M}_{K^{0}_{s}K^{0}_{s}} (GeV/#it{c}^{2})");
        hinvMass->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binwidthfile * 1000));

        cf2BWfit->cd(ipt + 1);
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

        // Significance calculation
        double lowRange = f1525Mass - 3 * f1525Width;
        double highRange = f1525Mass + 3 * f1525Width;
        double signal = BWonly->Integral(lowRange, highRange) / binwidthfile;
        int binlow = hraw->GetXaxis()->FindBin(lowRange);
        int binhigh = hraw->GetXaxis()->FindBin(highRange);
        double significance_den = TMath::Sqrt(hraw->Integral(binlow, binhigh));
        double significance = (significance_den > 0) ? signal / significance_den : 0;
        hSignificance->SetBinContent(ipt + 2, significance);

        double yieldIntegral = BWonly->Integral(f1525Mass - 3 * f1525Width, f1525Mass + 3 * f1525Width) / (binwidthfile * total_events * ptbinwidth);
        if (ipt == Npt - 1)
            yieldIntegral = yieldIntegral * 0.75;
        cout << "binwidthfile " << binwidthfile << ", total_events " << total_events << ", ptbinwidth " << ptbinwidth << endl;
        TMatrixDSym cov = fitResultptr->GetCovarianceMatrix();
        TMatrixDSym cov1525(3);
        cov.GetSub(0, 2, 0, 2, cov1525); // Extract covariance matrix for f2'(1525)
        Double_t *cov_array = cov1525.GetMatrixArray();
        Double_t *para = BWonly->GetParameters();
        double yield_err = BWonly->IntegralError(f1525Mass - 3 * f1525Width, f1525Mass + 3 * f1525Width, para, cov_array) / (ptbinwidth * binwidthfile * total_events);
        // cout << "pt bin " << lowpT << " - " << highpT << " GeV/c: Yield = " << yieldIntegral << " +/- " << yield_err << endl;

        // Corrections to the spectra
        double eff1525 = hefficiencyf2->GetBinContent(ipt + 1);
        double eff1525_err = hefficiencyf2->GetBinError(ipt + 1);
        double INELCrossSection = 77.904;
        double visibleCrossSection = 52.8;
        double oneUponTriggerEfficiency = INELCrossSection / visibleCrossSection;
        double BR_f2 = 0.222;
        double signalLossFactor = hEvbySigLoss->GetBinContent(ipt + 1);

        double correctedYield = yieldIntegral * signalLossFactor / (eff1525 * oneUponTriggerEfficiency * BR_f2);
        double corrected_yield1525_err = signalLossFactor * sqrt(pow(yield_err / eff1525, 2) + pow(yieldIntegral * eff1525_err / pow(eff1525, 2), 2)) / (oneUponTriggerEfficiency * BR_f2);

        hMass->SetBinContent(ipt + 2, BWpol3->GetParameter(1));
        hMass->SetBinError(ipt + 2, BWpol3->GetParError(1));
        hWidth->SetBinContent(ipt + 2, BWpol3->GetParameter(2));
        hWidth->SetBinError(ipt + 2, BWpol3->GetParError(2));

        hRawYield->SetBinContent(ipt + 2, yieldIntegral);
        hRawYield->SetBinError(ipt + 2, yield_err);

        hCorrectedYield->SetBinContent(ipt + 2, correctedYield);
        hCorrectedYield->SetBinError(ipt + 2, corrected_yield1525_err);
    } // end of pT loop

    TFile *OutputFile = new TFile((savepath + "f2_fit_1rBW.root").c_str(), "RECREATE");
    hMass->Write("hMass_1525");
    hWidth->Write("hWidth_1525");
    hRawYield->Write("hRawYield_1525");
    hCorrectedYield->Write("hCorrectedYield_1525");
    hSignificance->Write("hSignificance_1525");

    cf2BWfit->SaveAs((savepath + "mult_0-100/Spectra/plots/f2_fit_1rBW.png").c_str());

    TCanvas *cMass1525 = new TCanvas("cMass1525", "Mass vs pT", 720, 720);
    SetCanvasStyle(cMass1525, 0.15, 0.03, 0.05, 0.13);
    double pad1Size, pad2Size;
    canvas_style(cMass1525, pad1Size, pad2Size);
    cMass1525->cd(1);
    SetHistoQA(hMass);
    hMass->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMass->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    hMass->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hMass->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hMass->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hMass->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hMass->GetYaxis()->SetTitleOffset(1.3);
    hMass->GetYaxis()->SetRangeUser(1.485, 1.575);
    hMass->SetMarkerStyle(20);
    hMass->SetMarkerSize(1.5);
    hMass->SetLineColor(kBlue);
    hMass->SetMarkerColor(kBlue);
    hMass->Draw("pe");
    TFile *fDefault = new TFile(Form("%s/FitParam_Default2.root", savepath.c_str()), "READ");
    TH1F *hMassDefault = (TH1F *)fDefault->Get("Mult_0_100/hMass_1525");
    SetHistoQA(hMassDefault);
    hMassDefault->SetMarkerStyle(21);
    hMassDefault->SetMarkerSize(1.5);
    hMassDefault->SetLineColor(kRed);
    hMassDefault->SetMarkerColor(kRed);
    hMassDefault->Draw("pe same");

    TLegend *leg1525Mass = new TLegend(0.45, 0.65, 0.9, 0.93);
    leg1525Mass->SetBorderSize(0);
    leg1525Mass->SetFillStyle(0);
    leg1525Mass->SetTextSize(0.035);
    leg1525Mass->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    leg1525Mass->AddEntry((TObject *)0, "FT0M: 0-100%, |y|<0.5", "");
    leg1525Mass->AddEntry((TObject *)0, "f_{2}'(1525) Mass", "");
    leg1525Mass->AddEntry(hMassDefault, "Default Fit", "p");
    leg1525Mass->AddEntry(hMass, "1rBW Fit", "p");
    leg1525Mass->Draw();

    cMass1525->cd(2);
    TH1F *hRatio = (TH1F *)hMass->Clone("hRatio");
    int totalBins = hMass->GetNbinsX();
    for (int i = 1; i <= totalBins; i++)
    {
        double defaultVal = hMassDefault->GetBinContent(i);
        double defaultErr = hMassDefault->GetBinError(i);
        double newVal = hMass->GetBinContent(i + 1);
        double newErr = hMass->GetBinError(i + 1);
        if (defaultVal != 0)
        {
            double ratio = newVal / defaultVal;
            // cout << "Ratio is " << ratio << " for bin " << i << endl;
            double ratioErr = ratio * sqrt(pow(newErr / newVal, 2) + pow(defaultErr / defaultVal, 2));
            hRatio->SetBinContent(i + 1, ratio);
            hRatio->SetBinError(i + 1, ratioErr);
        }
        else
        {
            hRatio->SetBinContent(i + 1, 0);
            hRatio->SetBinError(i + 1, 0);
        }
    }
    SetHistoQA(hRatio);
    hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio->GetYaxis()->SetTitle("1rBW / Default");
    hRatio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hRatio->GetYaxis()->SetTitleSize(0.03 / pad2Size);
    hRatio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hRatio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hRatio->GetYaxis()->SetTitleOffset(0.8);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.5);
    hRatio->SetLineColor(kBlue);
    hRatio->SetMarkerColor(kBlue);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->SetMinimum(0.992);
    hRatio->SetMaximum(1.008);
    hRatio->Draw("pe");
    cMass1525->SaveAs((savepath + "mult_0-100/Spectra/plots/MassCompare_1rBW.png").c_str());
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

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.14);
    pad2->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.04);
}