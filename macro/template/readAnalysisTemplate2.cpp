#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <functional>
#include <TArrow.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <Math/MinimizerOptions.h>
#include "/home/subhadeep/softwares/alice_files/kstar_light_ion/kstar0/headers/style.h"
#include "initializations.h"

using namespace std;

bool plot_a4 = false;
bool plot_ppt = false;

void readAnalysisTemplate()
{
    TStopwatch timer;
    timer.Start();

    // Suppress standard ROOT terminal spam
    gErrorIgnoreLevel = kWarning;

    // Set global minimizer to Minuit2 / Migrad with Strategy 1
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);

    // =========================================================================
    // User configuration
    // =========================================================================
    const string kResBkg = "ROTATED"; // "MIX" | "LIKE" | "ROTATED"

    // Choose residual background model:
    // const string kbkg    = "pol1";
    // const string kbkg    = "pol2";
    const string kbkg = "pol3";
    // const string kbkg    = "expol";
    // const string kbkg    = "pol3Thresh";

    const bool widthFixed = true;
    // const bool widthFixed = false;

    // TString      outputtype = "pdf";
    TString outputtype = "png";
    const float txtsize = 0.045;

    // =========================================================================
    // Output folder
    // =========================================================================
    TString centRange = Form("%d_%d", lowcent, highcent);
    TString Cenoutputfolder = Form("output/%s/%s", kResBkg.c_str(), centRange.Data());
    if (gSystem->mkdir(Cenoutputfolder, kTRUE))
        cout << "Created output folder: " << Cenoutputfolder << endl;

    double lowfitrange[Npt + 20], highfitrange[Npt + 20];
    for (int i = 0; i < Npt; i++)
    {
        lowfitrange[i] = kFitRange[i][0];
        highfitrange[i] = kFitRange[i][1];
    }

    t2->SetNDC();
    t2->SetTextSize(0.06);
    t2->SetTextFont(42);

    // =========================================================================
    // Output canvases
    // =========================================================================
    TCanvas *cinv[Npt], *cSigbkg[Npt];
    for (int ip = 0; ip < Npt; ip++)
    {
        cinv[ip] = new TCanvas(Form("cinv%d", ip), TString::Format("cinv_pt_%2.1f-%2.1f", pT_bins[ip], pT_bins[ip + 1]).Data(), 10, 10, 720, 720);
        SetCanvasStyle(cinv[ip], 0.15, 0.05, 0.08, 0.13);

        cSigbkg[ip] = new TCanvas(Form("cSigbkg%d", ip), TString::Format("cSigbkg_pt_%2.1f-%2.1f", pT_bins[ip], pT_bins[ip + 1]).Data(), 720, 720);
        SetCanvasStyle(cSigbkg[ip], 0.15, 0.05, 0.08, 0.13);
    }

    Double_t significance_den, significance_num, ratio;

    // =========================================================================
    // Open input files
    // =========================================================================
    TFile *fInputFile = new TFile("../AnalysisResults728630.root", "READ");
    if (fInputFile->IsZombie())
    {
        cerr << "ERROR: Input file not found!" << endl;
        return;
    }

    TFile *fTemplateFile = TFile::Open(Form("../buildTemplate/template/%s/SignalMinusTrue.root", kResBkg.c_str()), "READ");
    if (!fTemplateFile || fTemplateFile->IsZombie())
    {
        cerr << "ERROR: SignalMinusTrue.root not found!" << endl;
        return;
    }
    cout << "Reflection template file opened successfully." << endl;

    // =========================================================================
    // Event counting
    // =========================================================================
    TH1F *hcent = (TH1F *)fInputFile->Get("kstar892-light-ion/eventSelection/hCentrality");
    if (!hcent)
    {
        cerr << "ERROR: hCentrality not found!" << endl;
        return;
    }
    int hcentbinlow = hcent->FindBin(lowcent + 1e-5);
    int hcentbinhigh = hcent->FindBin(highcent - 1e-5);
    double Event = hcent->Integral(hcentbinlow, hcentbinhigh);
    cout << "Number of events: " << Event << endl;

    // =========================================================================
    // Load invariant-mass 3D histograms
    // =========================================================================
    auto fDirectory = (TDirectoryFile *)fInputFile->Get("kstar892-light-ion/InvMass");
    TH3F *fHistLikeMM = (TH3F *)fDirectory->Get("h3KstarInvMasslikeSignMM");
    TH3F *fHistLikePP = (TH3F *)fDirectory->Get("h3KstarInvMasslikeSignPP");
    TH3F *fHistUnlike = (TH3F *)fDirectory->Get("h3KstarInvMassUnlikeSign");
    TH3F *fHistMix = (TH3F *)fDirectory->Get("h3KstarInvMassMixed");
    TH3F *fHistRotated = (TH3F *)fDirectory->Get("h3KstarInvMassRotated");
    if (!fHistUnlike || !fHistMix || !fHistLikeMM || !fHistLikePP || !fHistRotated)
    {
        cerr << "ERROR: Invariant-mass histograms not found!" << endl;
        return;
    }

    gstyle();
    gStyle->SetOptStat(0);

    std::vector<TCanvas *> c_sigbkg, c_fitsig, c_sigbyref;

    TH1D *hsigma = new TH1D("hsigma", "Resolution (#sigma_{res})", Npt, pT_bins);
    hsigma->Sumw2();

    // =========================================================================
    // =====================   pT   B I N   L O O P   =========================
    // =========================================================================
    for (Int_t ip = pt_start; ip < pt_end; ip++)
    {
        lowpt = pT_bins[ip];
        highpt = pT_bins[ip + 1];

        cout << "\n=============================================" << endl;
        cout << " pT bin [" << ip << "]: " << lowpt << " - " << highpt << " GeV/c" << endl;
        cout << "=============================================" << endl;

        // =====================================================================
        // Project 3D histogram -> 1D invariant mass
        // =====================================================================
        int lbincent = fHistUnlike->GetXaxis()->FindBin(lowcent + 1e-5);
        int hbincent = fHistUnlike->GetXaxis()->FindBin(highcent - 1e-5);
        int lbinpt = fHistUnlike->GetYaxis()->FindBin(lowpt + 1e-5);
        int hbinpt = fHistUnlike->GetYaxis()->FindBin(highpt - 1e-5);

        fHistTotal[ip] = fHistUnlike->ProjectionZ(
            Form("fHistUnlike_%d", ip), lbincent, hbincent, lbinpt, hbinpt);

        double energylow = fHistTotal[ip]->GetXaxis()->GetXmin();
        double energyhigh = fHistTotal[ip]->GetXaxis()->GetXmax();

        TH1D *hfsig = (TH1D *)fHistTotal[ip]->Clone();
        hfsig->Sumw2();
        double binwidth_file = (energyhigh - energylow) * kRebin[ip] / fHistTotal[ip]->GetNbinsX();
        cout << "  binwidth_file = " << binwidth_file * 1000 << " MeV/c^2" << endl;

        // =====================================================================
        // Combinatorial background subtraction
        // =====================================================================
        if (kResBkg == "MIX")
        {
            lbincent = fHistMix->GetXaxis()->FindBin(lowcent + 1e-5);
            hbincent = fHistMix->GetXaxis()->FindBin(highcent - 1e-5);
            lbinpt = fHistMix->GetYaxis()->FindBin(lowpt + 1e-5);
            hbinpt = fHistMix->GetYaxis()->FindBin(highpt - 1e-5);
            fHistBkg[ip] = fHistMix->ProjectionZ(
                Form("fHistMix_%d", ip), lbincent, hbincent, lbinpt, hbinpt);

            TH1D *bkgclone = (TH1D *)fHistBkg[ip]->Clone();
            sigbkg_integral = fHistTotal[ip]->Integral(
                fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]),
                fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1]));
            bkg_integral = bkgclone->Integral(
                bkgclone->GetXaxis()->FindBin(kNormRangepT[ip][0]),
                bkgclone->GetXaxis()->FindBin(kNormRangepT[ip][1]));
            normfactor = sigbkg_integral / bkg_integral;

            hfbkg = (TH1D *)bkgclone->Clone();
            hfbkg->Scale(normfactor);
            hfbkg->Rebin(kRebin[ip]);
            hfsig->Rebin(kRebin[ip]);
            hfsig->Add(hfbkg, -1);
            delete bkgclone;
        }
        else if (kResBkg == "LIKE")
        {
            lbincent = fHistLikeMM->GetXaxis()->FindBin(lowcent + 1e-5);
            hbincent = fHistLikeMM->GetXaxis()->FindBin(highcent - 1e-5);
            lbinpt = fHistLikeMM->GetYaxis()->FindBin(lowpt + 1e-5);
            hbinpt = fHistLikeMM->GetYaxis()->FindBin(highpt - 1e-5);
            fHistbkgMM[ip] = fHistLikeMM->ProjectionZ(Form("fHistLikeMM_%d", ip), lbincent, hbincent, lbinpt, hbinpt);

            lbincent = fHistLikePP->GetXaxis()->FindBin(lowcent + 1e-5);
            hbincent = fHistLikePP->GetXaxis()->FindBin(highcent - 1e-5);
            lbinpt = fHistLikePP->GetYaxis()->FindBin(lowpt + 1e-5);
            hbinpt = fHistLikePP->GetYaxis()->FindBin(highpt - 1e-5);
            fHistbkgPP[ip] = fHistLikePP->ProjectionZ(Form("fHistLikePP_%d", ip), lbincent, hbincent, lbinpt, hbinpt);

            auto tempLS = (TH1D *)fHistbkgMM[ip]->Clone("tempLS");
            tempLS->Multiply(fHistbkgPP[ip]);
            fHistbkgLS[ip] = (TH1D *)tempLS->Clone(Form("fHistbkgLS_%d", ip));
            fHistbkgLS[ip]->Reset();
            for (int ib = 1; ib <= tempLS->GetNbinsX(); ib++)
            {
                double ppnn = tempLS->GetBinContent(ib);
                double ppnnerr = tempLS->GetBinError(ib);
                if (ppnn <= 0)
                {
                    fHistbkgLS[ip]->SetBinContent(ib, 0);
                    fHistbkgLS[ip]->SetBinError(ib, 0);
                }
                else
                {
                    fHistbkgLS[ip]->SetBinContent(ib, 2 * std::sqrt(ppnn));
                    fHistbkgLS[ip]->SetBinError(ib, ppnnerr / std::sqrt(ppnn));
                }
            }
            delete tempLS;

            hfbkg = (TH1D *)fHistbkgLS[ip]->Clone();
            hfbkg->Rebin(kRebin[ip]);
            hfsig->Rebin(kRebin[ip]);
            hfsig->Add(hfbkg, -1);
        }
        else if (kResBkg == "ROTATED")
        {
            lbincent = fHistRotated->GetXaxis()->FindBin(lowcent + 1e-5);
            hbincent = fHistRotated->GetXaxis()->FindBin(highcent - 1e-5);
            lbinpt = fHistRotated->GetYaxis()->FindBin(lowpt + 1e-5);
            hbinpt = fHistRotated->GetYaxis()->FindBin(highpt - 1e-5);
            fHistRotated1D[ip] = fHistRotated->ProjectionZ(Form("fHistRotated_%d", ip), lbincent, hbincent, lbinpt, hbinpt);

            TH1D *bkgclone = (TH1D *)fHistRotated1D[ip]->Clone();
            sigbkg_integral = fHistTotal[ip]->Integral(
                fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][0]),
                fHistTotal[ip]->GetXaxis()->FindBin(kNormRangepT[ip][1]));
            bkg_integral = bkgclone->Integral(
                bkgclone->GetXaxis()->FindBin(kNormRangepT[ip][0]),
                bkgclone->GetXaxis()->FindBin(kNormRangepT[ip][1]));
            normfactor = sigbkg_integral / bkg_integral;

            hfbkg = (TH1D *)bkgclone->Clone();
            hfbkg->Scale(normfactor);
            hfbkg->Rebin(kRebin[ip]);
            hfsig->Rebin(kRebin[ip]);
            hfsig->Add(hfbkg, -1);
            delete bkgclone;
        }

        fHistTotal[ip]->Rebin(kRebin[ip]);
        ptbinwidth[ip] = pT_bins[ip + 1] - pT_bins[ip];

        // =====================================================================
        // Load reflection template
        // =====================================================================
        TString templName = Form("hSigminusTrue_pt_%.1f_%.1f", lowpt, highpt);
        TH1D *hReflRaw = (TH1D *)fTemplateFile->Get(templName);
        if (!hReflRaw)
        {
            cerr << "WARNING: Template '" << templName
                 << "' not found. Skipping pT bin " << ip << "." << endl;
            continue;
        }

        TH1D *hReflection = (TH1D *)hReflRaw->Clone(Form("hReflection_ip%d", ip));
        hReflection->Rebin(kRebin[ip]);

        TCanvas *cRefl = new TCanvas(Form("cRefl_ip%d", ip), "Reflection template", 720, 720);
        TH1D *hDatabyReflection = (TH1D *)hfsig->Clone(Form("hDatabyReflection_ip%d", ip));
        hDatabyReflection->SetTitle(Form("%.1f < p_{T} (GeV/c) < %.1f; M_{K#pi} (GeV/c^{2}); Data / Reflection template", lowpt, highpt));

        TH1D *hRefNorm = (TH1D *)hReflection->Clone(Form("hRefNorm_ip%d", ip));
        cout << "Ref Norm = " << hDatabyReflection->Integral(hDatabyReflection->GetXaxis()->FindBin(0.7), hDatabyReflection->GetXaxis()->FindBin(0.8)) / (hRefNorm->Integral(hRefNorm->GetXaxis()->FindBin(0.7), hRefNorm->GetXaxis()->FindBin(0.8))) << endl;

        hRefNorm->Scale(hDatabyReflection->Integral(hDatabyReflection->GetXaxis()->FindBin(0.7), hDatabyReflection->GetXaxis()->FindBin(0.8)) / (hRefNorm->Integral(hRefNorm->GetXaxis()->FindBin(0.7), hRefNorm->GetXaxis()->FindBin(0.8))));
        hDatabyReflection->Divide(hRefNorm);
        hDatabyReflection->GetXaxis()->SetRangeUser(kFitRange[ip][0], kFitRange[ip][1]);
        cRefl->cd();
        hDatabyReflection->Draw();

        TLatex latexref;
        latexref.SetNDC();
        latexref.SetTextSize(0.05);
        latexref.SetTextAlign(22);
        latexref.DrawLatex(0.5, 0.95, Form("%.1f < p_{T} (GeV/c) < %.1f", lowpt, highpt));

        auto c_clone_sigbyref = (TCanvas *)cRefl->Clone(Form("hsigbyref_pt_%d", ip + 1));
        c_sigbyref.push_back(c_clone_sigbyref);
        cRefl->SaveAs(Cenoutputfolder + Form("/hDatabyReflection_ip%d.%s", ip, outputtype.Data()));

        for (int ib = 1; ib <= hReflection->GetNbinsX(); ib++)
        {
            if (hReflection->GetBinContent(ib) < 0)
            {
                hReflection->SetBinContent(ib, 0.0);
                hReflection->SetBinError(ib, 0.0);
            }
        }

        bool bTemplEmpty = (hReflection->Integral() <= 0);
        if (bTemplEmpty)
            cerr << "WARNING: Template empty in fit range for pT bin " << ip << endl;

        // =====================================================================
        // Superior Analytical Normalization Setup
        // =====================================================================
        const double fitLo = kFitRange[ip][0];
        const double fitHi = kFitRange[ip][1];
        const double leftLo = fitLo;
        const double leftHi = 0.82;
        const double binWidthFit = hfsig->GetBinWidth(1);

        int binLo = hReflection->GetXaxis()->FindBin(fitLo + 1e-6);
        int binHi = hReflection->GetXaxis()->FindBin(fitHi - 1e-6);
        double reflBinWidth = hReflection->GetBinWidth(1);
        double reflNorm = std::max(hReflection->Integral(binLo, binHi) * reflBinWidth, 1e-12);

        double total_in_fit = std::max(1.0, hfsig->Integral(
                                                hfsig->GetXaxis()->FindBin(fitLo + 1e-5),
                                                hfsig->GetXaxis()->FindBin(fitHi - 1e-5)));

        // ==================================================================
        // Define Universal Residual Background Evaluator Lambda
        // ==================================================================
        int nBkgPars = 1; // pol1 default
        if (kbkg == "pol2")
            nBkgPars = 2;
        else if (kbkg == "pol3")
            nBkgPars = 3;
        else if (kbkg == "expol" || kbkg == "pol3Thresh")
            nBkgPars = 4;

        auto evalBkgDensity = [=](double x, const double *par) -> double
        {
            if (kbkg == "pol1")
            {
                double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
                double chebRaw = 1.0 + par[0] * t;
                return chebRaw / (fitHi - fitLo);
            }
            else if (kbkg == "pol2")
            {
                double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
                double T1 = t, T2 = 2.0 * t * t - 1.0;
                double chebRaw = 1.0 + par[0] * T1 + par[1] * T2;
                double chebInt = ((fitHi - fitLo) / 2.0) * (2.0 - (2.0 / 3.0) * par[1]);
                return chebRaw / std::max(chebInt, 1e-12);
            }
            else if (kbkg == "pol3")
            {
                double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
                double T1 = t, T2 = 2.0 * t * t - 1.0, T3 = 4.0 * t * t * t - 3.0 * t;
                double chebRaw = 1.0 + par[0] * T1 + par[1] * T2 + par[2] * T3;
                double chebInt = ((fitHi - fitLo) / 2.0) * (2.0 - (2.0 / 3.0) * par[1]);
                return chebRaw / std::max(chebInt, 1e-12);
            }
            else if (kbkg == "expol")
            {
                if (x <= 0)
                    return 0.0;
                auto rawExpol = [](double val, const double *p)
                {
                    if (val <= 0)
                        return 0.0;
                    return std::pow(val, p[0]) * std::exp(-p[1] - p[2] * val - p[3] * val * val);
                };
                double rawVal = rawExpol(x, par);
                int nSteps = 40;
                double h = (fitHi - fitLo) / nSteps;
                double sum = rawExpol(fitLo, par) + rawExpol(fitHi, par);
                for (int j = 1; j < nSteps; j++)
                {
                    double xj = fitLo + j * h;
                    sum += rawExpol(xj, par) * ((j % 2 == 0) ? 2.0 : 4.0);
                }
                double normInt = (h / 3.0) * sum;
                return rawVal / std::max(normInt, 1e-12);
            }
            else if (kbkg == "pol3Thresh")
            {
                double m_thresh = massPi + massKa;
                if (x < m_thresh)
                    return 0.0;
                double dx = x - m_thresh;
                double rawVal = par[0] + par[1] * dx + par[2] * dx * dx + par[3] * dx * dx * dx;

                auto antideriv = [m_thresh](double val, const double *p)
                {
                    if (val <= m_thresh)
                        return 0.0;
                    double d = val - m_thresh;
                    return p[0] * d + 0.5 * p[1] * d * d + (1.0 / 3.0) * p[2] * d * d * d + 0.25 * p[3] * d * d * d * d;
                };
                double normInt = antideriv(fitHi, par) - antideriv(std::max(fitLo, m_thresh), par);
                return rawVal / std::max(normInt, 1e-12);
            }
            double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
            return (1.0 + par[0] * t) / (fitHi - fitLo);
        };

        // ==================================================================
        // Left-sideband pre-fit
        // ==================================================================
        TF1 *fLeft = new TF1(Form("fLeft_ip%d", ip), [hReflection, reflNorm, fitLo, fitHi, binWidthFit](double *xx, double *pp) -> double
                             {
                double x = xx[0];
                double Ncorr = pp[0], Nbkg = pp[1], c1 = pp[2];
                double corrDensity = hReflection->Interpolate(x) / reflNorm;
                double t = 2.0 * (x - fitLo) / (fitHi - fitLo) - 1.0;
                double bkgDensity = (1.0 + c1 * t) / (fitHi - fitLo);
                return binWidthFit * (Ncorr * corrDensity + Nbkg * bkgDensity); }, leftLo, leftHi, 3);
        fLeft->SetParameters(total_in_fit * 0.3, total_in_fit * 0.3, 0.0);
        fLeft->SetParLimits(0, 0.0, 1e10);
        fLeft->SetParLimits(1, 0.0, 1e10);
        fLeft->SetParLimits(2, -1.0, 1.0);

        hfsig->Fit(fLeft, "R Q N 0");
        double N_temp_0 = fLeft->GetParameter(0);
        cout << "  -> Left SB pre-fit: N_temp_0 = " << N_temp_0 << " +/- " << fLeft->GetParError(0) << endl;
        delete fLeft;

        // ==================================================================
        // Full invariant-mass VOIGTIAN fit
        // ==================================================================
        double max_corr = total_in_fit;

        TF1 *fTotal = new TF1(Form("fTotal_ip%d", ip), [hReflection, reflNorm, fitLo, fitHi, binWidthFit, evalBkgDensity](double *xx, double *pp) -> double
                              {
                double x = xx[0];
                double Nsig = pp[0], mass = pp[1], width_val = pp[2], sigma = pp[3]; // Slot 3: Gaussian sigma
                double Ncorr = pp[4], Nbkg = pp[5];                   // Shifted up by 1

                double voigtRaw = TMath::Voigt(x - mass, sigma, width_val);

                // Fast Numerical Normalization dynamically valid for all pT (low and high sigma)
                int nSteps = 40; 
                double h_step = (fitHi - fitLo) / nSteps;
                double vInt = TMath::Voigt(fitLo - mass, sigma, width_val) + TMath::Voigt(fitHi - mass, sigma, width_val);
                for (int j = 1; j < nSteps; j++) {
                    double xj = fitLo + j * h_step;
                    vInt += TMath::Voigt(xj - mass, sigma, width_val) * ((j % 2 == 0) ? 2.0 : 4.0);
                }
                vInt *= h_step / 3.0;
                double sigDensity = voigtRaw / std::max(vInt, 1e-12);

                // 2. Unit-Normalized MC Correlated Template Density
                double corrDensity = 0.0;
                if (hReflection && hReflection->GetEntries() > 0) {
                    corrDensity = hReflection->Interpolate(x) / reflNorm;
                    if (corrDensity < 0) corrDensity = 0.0;
                }

                // 3. Normalized Residual Background Density (Starts at pp[6])
                double bkgDensity = evalBkgDensity(x, &pp[6]);

                return binWidthFit * (Nsig * sigDensity + Ncorr * corrDensity + Nbkg * bkgDensity); }, fitLo, fitHi, 6 + nBkgPars); // 6 base params + bkg shape params

        fTotal->SetParNames("Nsig", "mass", "width", "sigma", "Ncorr", "Nbkg");
        fTotal->SetParameter(0, total_in_fit * 0.40);
        fTotal->SetParLimits(0, 0.0, total_in_fit);
        fTotal->SetParameter(1, masspdg);
        fTotal->SetParLimits(1, masspdg - 0.010, masspdg + 0.010);

        // --- SLOT 2: WIDTH CONFIGURATION FOR SYSTEMATICS ---
        fTotal->SetParameter(2, widthpdg);
        if (widthFixed)
            fTotal->FixParameter(2, widthpdg);
        else
            fTotal->SetParLimits(2, widthpdg - 0.005, widthpdg + 0.005);
        // For Systematics: comment out FixParameter above and use:
        // fTotal->SetParLimits(2, 0.030, 0.070); // Let it float
        // OR: fTotal->FixParameter(2, widthpdg * 1.10); // Fix at +10% variation
        // ---------------------------------------------------

        // --- SLOT 3: GAUSSIAN SIGMA FLOATING ---
        fTotal->SetParameter(3, 0.002);
        fTotal->SetParLimits(3, 0.0001, 0.030); // Allow up to 30 MeV for high-pT smearing

        fTotal->SetParameter(4, std::min(std::max(N_temp_0, 0.0), max_corr));
        fTotal->SetParLimits(4, 0.0, 2 * N_temp_0);
        fTotal->SetParameter(5, total_in_fit * 0.30);
        fTotal->SetParLimits(5, 0.0, total_in_fit);

        // Dynamically initialize background shape parameters starting at index 6
        if (kbkg == "pol1")
        {
            fTotal->SetParName(6, "c1");
            fTotal->SetParameter(6, -0.5);
            fTotal->SetParLimits(6, -1.0, 1.0);
        }
        else if (kbkg == "pol2")
        {
            fTotal->SetParName(6, "c1");
            fTotal->SetParName(7, "c2");
            fTotal->SetParameter(6, -0.5);
            fTotal->SetParLimits(6, -1.0, 1.0);
            fTotal->SetParameter(7, 0.1);
            fTotal->SetParLimits(7, -1.0, 1.0);
        }
        else if (kbkg == "pol3")
        {
            fTotal->SetParName(6, "c1");
            fTotal->SetParName(7, "c2");
            fTotal->SetParName(8, "c3");
            fTotal->SetParameter(6, -0.5);
            fTotal->SetParLimits(6, -1.0, 1.0);
            fTotal->SetParameter(7, 0.1);
            fTotal->SetParLimits(7, -1.0, 1.0);
            fTotal->SetParameter(8, -0.05);
            fTotal->SetParLimits(8, -1.0, 1.0);
        }
        else if (kbkg == "expol")
        {
            fTotal->SetParName(6, "p0");
            fTotal->SetParName(7, "p1");
            fTotal->SetParName(8, "p2");
            fTotal->SetParName(9, "p3");
            fTotal->SetParameter(6, 1.0);
            fTotal->SetParLimits(6, -15.0, 15.0);
            fTotal->SetParameter(7, 1.0);
            fTotal->SetParLimits(7, -20.0, 20.0);
            fTotal->SetParameter(8, 1.0);
            fTotal->SetParLimits(8, -15.0, 15.0);
            fTotal->SetParameter(9, 0.1);
            fTotal->SetParLimits(9, -15.0, 15.0);
        }
        else if (kbkg == "pol3Thresh")
        {
            fTotal->SetParName(6, "p0");
            fTotal->SetParName(7, "p1");
            fTotal->SetParName(8, "p2");
            fTotal->SetParName(9, "p3");
            fTotal->SetParameter(6, 1.0);
            fTotal->SetParLimits(6, 0.0, 1000.0);
            fTotal->SetParameter(7, 1.0);
            fTotal->SetParLimits(7, -50.0, 50.0);
            fTotal->SetParameter(8, 0.0);
            fTotal->SetParLimits(8, -50.0, 50.0);
            fTotal->SetParameter(9, 0.0);
            fTotal->SetParLimits(9, -50.0, 50.0);
        }

        if (bTemplEmpty)
            fTotal->FixParameter(4, 0.0); // Shifted to Slot 4

        // Execute fit with RQS options
        TFitResultPtr fitResult = hfsig->Fit(fTotal, "R Q S+");
        int fitStatus = fitResult;
        int covStatus = fitResult->CovMatrixStatus();

        // Superior Geometric Simplex Fallback Pipeline
        if (fitStatus != 0 || covStatus < 2)
        {
            cout << "  Retrying: Simplex -> Migrad (Strategy 2)..." << endl;
            ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Simplex");
            hfsig->Fit(fTotal, "R Q N");
            ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
            ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
            fitResult = hfsig->Fit(fTotal, "R Q S");
            ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
        }

        fitStatus = fitResult;
        covStatus = fitResult->CovMatrixStatus();
        cout << "  Fit status: " << fitStatus << "  |  Cov status: " << covStatus << endl;
        if (fitStatus != 0 || covStatus < 2)
            cout << "  WARNING: Bad fit in pT bin " << ip << endl;

        // =====================================================================
        // Extract parameters directly as physical integrated counts
        // =====================================================================
        double N_sig = fTotal->GetParameter(0);
        double N_sig_err = fTotal->GetParError(0);
        double mass_fit = fTotal->GetParameter(1);
        double mass_err = fTotal->GetParError(1);
        double width_fit = fTotal->GetParameter(2);
        double width_err = fTotal->GetParError(2);
        double sigma_fit = fTotal->GetParameter(3); // Gaussian Sigma
        double sigma_err = fTotal->GetParError(3);
        double N_temp_fit = fTotal->GetParameter(4); // Shifted to Slot 4
        double N_res_fit = fTotal->GetParameter(5);  // Shifted to Slot 5

        Mass[ip] = mass_fit;
        ErrorMass[ip] = mass_err;
        Width[ip] = width_fit;      // Now stores your fitted or varied width!
        ErrorWidth[ip] = width_err; // Now stores your width error!
        Yield[ip] = N_sig;

        cout << "  Mass     = " << Mass[ip] << " +/- " << ErrorMass[ip] << " GeV/c^2" << endl;
        cout << "  Width    = " << Width[ip] << " +/- " << ErrorWidth[ip] << " GeV/c^2" << endl;
        cout << "  Res(Sig) = " << sigma_fit << " +/- " << sigma_err << " GeV/c^2" << endl;
        cout << "  N_sig    = " << N_sig << " +/- " << N_sig_err << endl;
        cout << "  N_temp   = " << N_temp_fit << endl;
        cout << "  N_res    = " << N_res_fit << endl;

        // =====================================================================
        // Chi²/NDF
        // =====================================================================
        double chi2 = fTotal->GetChisquare();
        int ndf = fTotal->GetNDF();
        Chi2Ndf[ip] = (ndf > 0) ? chi2 / ndf : 0.0;
        hChiSquare->SetBinContent(ip + 1, Chi2Ndf[ip]);
        cout << "  chi2/NDF = " << Chi2Ndf[ip] << "  (" << ndf << " NDF)" << endl;

        // =====================================================================
        // Yield integration calculations (Voigtian Corrected)
        // =====================================================================
        auto templRawIntegral = [&](double lo, double hi) -> double
        {
            int bLo = hReflection->GetXaxis()->FindBin(lo + 1e-6);
            int bHi = hReflection->GetXaxis()->FindBin(hi - 1e-6);
            return hReflection->Integral(bLo, bHi) * reflBinWidth;
        };

        double m5lo = masspdg - 5 * width_fit, m5hi = masspdg + 5 * width_fit;
        double m2lo = masspdg - 2 * width_fit, m2hi = masspdg + 2 * width_fit;

        // Numerically integrate the exact fitted Voigtian shape
        TF1 fTempVoigt("fTempVoigt", "[0]*TMath::Voigt(x - [1], [3], [2])", 0.0, 2.0);
        fTempVoigt.SetParameters(1.0, mass_fit, width_fit, sigma_fit);

        double voigtIntegralFitRange = fTempVoigt.Integral(fitLo, fitHi);
        double f_5g = fTempVoigt.Integral(m5lo, m5hi) / std::max(voigtIntegralFitRange, 1e-12);
        double f_2g = fTempVoigt.Integral(m2lo, m2hi) / std::max(voigtIntegralFitRange, 1e-12);

        double N_sig_5g = N_sig * f_5g;
        double N_sig_2g = N_sig * f_2g;

        cout << "  Voigtian fractions: +/-5Γ=" << f_5g << "  +/-2Γ=" << f_2g << endl;

        // =====================================================================
        // Raw yield – integral method
        // =====================================================================
        yieldcalc = N_sig_5g / (Event * ptbinwidth[ip] * dy * BR);
        yielderror = (N_sig_err * f_5g) / (Event * ptbinwidth[ip] * dy * BR);

        hintegral_yield->SetBinContent(ip + 1, yieldcalc);
        hintegral_yield->SetBinError(ip + 1, yielderror);
        hmass->SetBinContent(ip + 1, Mass[ip]);
        hmass->SetBinError(ip + 1, ErrorMass[ip]);
        hwidth->SetBinContent(ip + 1, Width[ip]);
        hwidth->SetBinError(ip + 1, ErrorWidth[ip]);

        std::cout << "Yield (Functional integration): " << yieldcalc << " #pm " << yielderror << std::endl;

        // =====================================================================
        // Significance
        // =====================================================================
        bmin = hfsig->GetXaxis()->FindBin(masspdg - 2 * width_fit);
        bmax = hfsig->GetXaxis()->FindBin(masspdg + 2 * width_fit);
        significance_num = N_sig_2g;
        significance_den = TMath::Sqrt(std::max(1.0, (double)fHistTotal[ip]->Integral(bmin, bmax)));
        ratio = significance_num / significance_den;
        hsignificance->SetBinContent(ip + 1, ratio);

        // =====================================================================
        // Raw yield – bin-counting method with Snapped Boundaries Fix
        // =====================================================================
        Yield_bincount_hist = hfsig->IntegralAndError(bmin, bmax, hBCError_1);

        // --- SNAPPED BOUNDARIES FIX START ---
        // IntegralAndError sums WHOLE bins, so the data integral's true range
        // is [GetBinLowEdge(bmin), GetBinUpEdge(bmax)] -- not the exact continuous
        // masspdg+/-2*width_fit window. We snap the continuous limits to match!
        double m2lo_snap = hfsig->GetXaxis()->GetBinLowEdge(bmin);
        double m2hi_snap = hfsig->GetXaxis()->GetBinUpEdge(bmax);

        double corr_f_fit = templRawIntegral(fitLo, fitHi);
        double corr_f_2g = templRawIntegral(m2lo_snap, m2hi_snap) / std::max(corr_f_fit, 1e-12);

        TF1 fTempBkg("fTempBkg", [&](double *xx, double *pp)
                     { return evalBkgDensity(xx[0], pp); }, fitLo, fitHi, nBkgPars);
        fTempBkg.SetParameters(fTotal->GetParameters() + 6); // Shifted to Offset 6

        double bkg_f_fit = fTempBkg.Integral(fitLo, fitHi);
        double bkg_f_2g = fTempBkg.Integral(m2lo_snap, m2hi_snap) / std::max(bkg_f_fit, 1e-12);

        // Calculate the snapped Voigtian fraction so tail correction is exact
        double f_2g_snap = fTempVoigt.Integral(m2lo_snap, m2hi_snap) / std::max(voigtIntegralFitRange, 1e-12);
        double N_sig_2g_snap = N_sig * f_2g_snap;

        double N_temp_2g = N_temp_fit * corr_f_2g;
        double N_res_2g = N_res_fit * bkg_f_2g;
        fYield_BinCount = Yield_bincount_hist - N_temp_2g - N_res_2g;

        double tail_correction = N_sig_5g - N_sig_2g_snap;
        Total_Ybincounting = (fYield_BinCount + tail_correction) / (Event * ptbinwidth[ip] * dy * BR);
        // --- SNAPPED BOUNDARIES FIX END ---

        double Final_pro_error = hBCError_1 / (Event * ptbinwidth[ip] * dy * BR);

        hYbincount->SetBinContent(ip + 1, Total_Ybincounting);
        hYbincount->SetBinError(ip + 1, Final_pro_error);
        hFrac_stat_error->SetBinContent(ip + 1, (Total_Ybincounting > 0) ? Final_pro_error / Total_Ybincounting : 0.0);

        std::cout << "Yield (Bin Count): " << Total_Ybincounting << " #pm " << Final_pro_error << std::endl;

        hsigma->SetBinContent(ip + 1, sigma_fit);
        hsigma->SetBinError(ip + 1, sigma_err);

        // =====================================================================
        // PLOTTING – fitted invariant mass + Ratio Pad
        // =====================================================================
        cinv[ip]->cd();
        cinv[ip]->Clear();

        // --- 1. Upper Pad: Main Fit ---
        TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.0);
        pad1->SetLeftMargin(0.12);
        pad1->SetRightMargin(0.035);
        pad1->Draw();
        pad1->cd();

        if (hfsig->GetFunction(Form("fTotal_ip%d", ip)))
        {
            hfsig->GetFunction(Form("fTotal_ip%d", ip))->SetBit(TF1::kNotDraw);
        }

        hfsig->SetMarkerStyle(20);
        hfsig->SetMarkerSize(0.5);
        hfsig->SetMarkerColor(kBlack);
        hfsig->SetLineColor(kBlack);
        hfsig->GetXaxis()->SetRangeUser(kFitRange[ip][0], kFitRange[ip][1]);
        hfsig->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/c^{2}", binwidth_file * 1000));
        hfsig->GetYaxis()->CenterTitle(1);
        hfsig->GetYaxis()->SetMaxDigits(2);
        hfsig->GetXaxis()->SetLabelSize(0);
        hfsig->GetXaxis()->SetTitleSize(0);
        hfsig->GetYaxis()->SetTitleSize(0.055);
        hfsig->GetYaxis()->SetLabelSize(0.045);
        hfsig->GetYaxis()->SetTitleOffset(1.05);
        hfsig->SetStats(0);

        if (hfsig->GetMaximum() > 0)
            hfsig->SetMinimum(-hfsig->GetMaximum() * 0.08);

        hfsig->SetMaximum(hfsig->GetMaximum() * 1.3);
        hfsig->Draw("E");

        const int kSmoothNpx = 1000;
        fTotal->SetLineColor(kRed);
        fTotal->SetLineWidth(2); // fTotal->SetNpx(kSmoothNpx);
        fTotal->Draw("SAME");

        // Analytical Component Functions for Plotting
        TF1 *fSig = new TF1(Form("fSig_ip%d", ip), [fitLo, fitHi, binWidthFit](double *xx, double *pp) -> double
                            {
                double x = xx[0]; double Nsig_val = pp[0], mass_val = pp[1], width_val = pp[2], sigma_val = pp[3];
                double voigtRaw = TMath::Voigt(x - mass_val, sigma_val, width_val);
                int nSteps = 40; 
                double h_step = (fitHi - fitLo) / nSteps;
                double vInt = TMath::Voigt(fitLo - mass_val, sigma_val, width_val) + TMath::Voigt(fitHi - mass_val, sigma_val, width_val);
                for (int j = 1; j < nSteps; j++) {
                    double xj = fitLo + j * h_step;
                    vInt += TMath::Voigt(xj - mass_val, sigma_val, width_val) * ((j % 2 == 0) ? 2.0 : 4.0);
                }
                vInt *= h_step / 3.0;
                return binWidthFit * Nsig_val * voigtRaw / std::max(vInt, 1e-12); }, fitLo, fitHi, 4);
        fSig->SetParameters(N_sig, mass_fit, width_fit, sigma_fit);
        fSig->SetLineColor(kMagenta + 1);
        fSig->SetLineStyle(kSolid);
        fSig->SetLineWidth(2); // fSig->SetNpx(kSmoothNpx);
        fSig->Draw("SAME");

        TF1 *fCorr = new TF1(Form("fCorr_ip%d", ip), [hReflection, reflNorm, binWidthFit](double *xx, double *pp) -> double
                             {
                double x = xx[0]; double Ncorr_val = pp[0];
                double corrDen = hReflection->Interpolate(x) / reflNorm;
                return binWidthFit * Ncorr_val * corrDen; }, fitLo, fitHi, 1);
        fCorr->SetParameter(0, N_temp_fit);
        fCorr->SetLineColor(kGreen + 2);
        fCorr->SetLineStyle(kSolid);
        fCorr->SetLineWidth(2); // fCorr->SetNpx(kSmoothNpx);
        fCorr->Draw("SAME");

        TF1 *fBkg = new TF1(Form("fBkg_ip%d", ip), [fitLo, fitHi, binWidthFit, evalBkgDensity](double *xx, double *pp) -> double
                            {
                double x = xx[0];
                double Nbkg_val = pp[0];
                double bkgDen = evalBkgDensity(x, &pp[1]);
                return binWidthFit * Nbkg_val * bkgDen; }, fitLo, fitHi, 1 + nBkgPars);
        fBkg->SetParameter(0, N_res_fit);
        for (int ib = 0; ib < nBkgPars; ++ib)
            fBkg->SetParameter(1 + ib, fTotal->GetParameter(6 + ib));
        fBkg->SetLineColor(kBlue);
        fBkg->SetLineStyle(kSolid);
        fBkg->SetLineWidth(2); // fBkg->SetNpx(kSmoothNpx);
        fBkg->Draw("SAME");

        TF1 *fTotalBkg = new TF1(Form("fTotalBkg_ip%d", ip), [hReflection, reflNorm, fitLo, fitHi, binWidthFit, evalBkgDensity](double *xx, double *pp) -> double
                                 {
                double x = xx[0];
                double Ncorr_val = pp[0], Nbkg_val = pp[1];
                double corrDen = 0.0;
                if (hReflection && hReflection->GetEntries() > 0) {
                    corrDen = hReflection->Interpolate(x) / reflNorm;
                    if (corrDen < 0) corrDen = 0.0;
                }
                double bkgDen = evalBkgDensity(x, &pp[2]);
                return binWidthFit * (Ncorr_val * corrDen + Nbkg_val * bkgDen); }, fitLo, fitHi, 2 + nBkgPars);
        fTotalBkg->SetParameter(0, N_temp_fit);
        fTotalBkg->SetParameter(1, N_res_fit);
        for (int ib = 0; ib < nBkgPars; ++ib)
            fTotalBkg->SetParameter(2 + ib, fTotal->GetParameter(6 + ib));
        fTotalBkg->SetLineColor(kOrange + 7);
        fTotalBkg->SetLineStyle(kDashed);
        fTotalBkg->SetLineWidth(2); // fTotalBkg->SetNpx(kSmoothNpx);
        fTotalBkg->Draw("SAME");

        auto mkLine = [](Color_t col, Style_t sty = kSolid, Width_t wid = 2) -> TLine *
        {
            auto *l = new TLine();
            l->SetLineColor(col);
            l->SetLineStyle(sty);
            l->SetLineWidth(wid);
            return l;
        };

        TLegend *legComp = new TLegend(0.13, 0.6, 0.4, 0.91);
        legComp->SetBorderSize(0);
        legComp->SetFillStyle(0);
        legComp->SetTextFont(42);
        legComp->SetTextSize(0.045);
        legComp->AddEntry(hfsig, "Data", "ep");
        legComp->AddEntry(mkLine(kRed), "Total Fit", "l");
        legComp->AddEntry(mkLine(kMagenta + 1, kSolid), "Voigtian Sig", "l");
        legComp->AddEntry(mkLine(kGreen + 2, kSolid), "MC Template", "l");

        TString bkgLegendLabel = (kbkg == "expol") ? "Exp. Pol" : (kbkg == "pol3Thresh") ? "Pol3 Thresh"
                                                                                         : Form("Cheby %s", kbkg.c_str());
        // legComp->AddEntry(mkLine(kBlue, kSolid), Form("Residual Bkg (%s)", bkgLegendLabel.Data()), "l");
        legComp->AddEntry(mkLine(kBlue, kSolid), bkgLegendLabel.Data(), "l");
        legComp->AddEntry(mkLine(kOrange + 7, kDashed), "Total Bkg (Temp+Res)", "l");
        legComp->Draw();

        TLegend *legPars = new TLegend(0.49, 0.6, 0.93, 0.89);
        legPars->SetBorderSize(0);
        legPars->SetFillStyle(0);
        legPars->SetTextFont(42);
        legPars->SetTextSize(0.04);
        legPars->AddEntry((TObject *)0, Form("Mass: %.3f #pm %.1e GeV/c^{2}", Mass[ip], ErrorMass[ip]), "");

        if (widthFixed)
        {
            legPars->AddEntry((TObject *)0,
                              Form("Width: %.1f MeV/c^{2} (fixed)", width_fit * 1000), "");
        }
        else
        {
            legPars->AddEntry((TObject *)0,
                              Form("Width: %.1f #pm %.1f MeV/c^{2}",
                                   width_fit * 1000, width_err * 1000),
                              "");
        }
        legPars->AddEntry((TObject *)0, Form("#sigma_{res}: %.1f #pm %.1f MeV/c^{2}", sigma_fit * 1000, sigma_err * 1000), "");
        legPars->AddEntry((TObject *)0, Form("N_{sig}: %.0f #pm %.0f", N_sig, N_sig_err), "");
        legPars->AddEntry((TObject *)0, Form("N_{temp}: %.0f", N_temp_fit), "");
        legPars->AddEntry((TObject *)0, Form("N_{res}: %.0f", N_res_fit), "");
        legPars->AddEntry((TObject *)0, Form("#chi^{2}/NDF: %.2f", Chi2Ndf[ip]), "");
        legPars->Draw();

        t2->DrawLatex(0.25, 0.94, Form("#bf{%.1f < #it{p}_{T}(GeV/#it{c}) < %.1f}", pT_bins[ip], pT_bins[ip + 1]));

        // --- 2. Lower Pad: Data / Fit Ratio ---
        cinv[ip]->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);
        pad2->SetTopMargin(0.0);
        pad2->SetBottomMargin(0.35);
        pad2->SetLeftMargin(0.12);
        pad2->SetRightMargin(0.035);
        pad2->Draw();
        pad2->cd();

        TGraphAsymmErrors *gRatio = new TGraphAsymmErrors();
        int pt_idx = 0;
        for (int ibin = 1; ibin <= hfsig->GetNbinsX(); ibin++)
        {
            double xp = hfsig->GetBinCenter(ibin);
            if (xp < kFitRange[ip][0] || xp > kFitRange[ip][1])
                continue;
            double yp = hfsig->GetBinContent(ibin);
            double yfit = fTotal->Eval(xp);
            if (yfit > 0)
            {
                gRatio->SetPoint(pt_idx, xp, yp / yfit);
                gRatio->SetPointError(pt_idx, 0, 0,
                                      hfsig->GetBinErrorLow(ibin) / yfit,
                                      hfsig->GetBinErrorUp(ibin) / yfit);
                pt_idx++;
            }
        }

        TH1D *hRatioFrame = new TH1D(Form("hRatioFrame_%d", ip), "", 100, kFitRange[ip][0], kFitRange[ip][1]);
        hRatioFrame->GetXaxis()->SetTitle("M_{K#pi} (GeV/c^{2})");
        hRatioFrame->GetYaxis()->SetTitle("Data / Fit");
        hRatioFrame->SetStats(0);

        hRatioFrame->GetXaxis()->SetTitleSize(0.13);
        hRatioFrame->GetXaxis()->SetLabelSize(0.11);
        hRatioFrame->GetXaxis()->SetTitleOffset(1.05);
        hRatioFrame->GetYaxis()->SetTitleSize(0.12);
        hRatioFrame->GetYaxis()->SetLabelSize(0.10);
        hRatioFrame->GetYaxis()->SetTitleOffset(0.4);
        hRatioFrame->GetYaxis()->SetNdivisions(505);

        hRatioFrame->SetMinimum(0.85);
        hRatioFrame->SetMaximum(1.15);
        hRatioFrame->Draw("AXIS");

        gRatio->SetMarkerStyle(20);
        gRatio->SetMarkerSize(0.5);
        gRatio->SetLineColor(kBlack);
        gRatio->SetMarkerColor(kBlack);
        gRatio->Draw("P SAME");

        TLine *lineRatio = new TLine(kFitRange[ip][0], 1.0, kFitRange[ip][1], 1.0);
        lineRatio->SetLineColor(kRed);
        lineRatio->SetLineStyle(2);
        lineRatio->Draw("SAME");

        // --- Finalize and Save ---
        cinv[ip]->cd();
        auto c_clone_fit = (TCanvas *)cinv[ip]->Clone(Form("hfitsig_pt_%d", ip + 1));
        c_fitsig.push_back(c_clone_fit);
        cinv[ip]->SaveAs(Form((Cenoutputfolder + "/hfitsig_pt%d." + outputtype).Data(), ip + 1));

        // =====================================================================
        // PLOTTING – signal + combinatorial background
        // =====================================================================
        cSigbkg[ip]->cd();
        TH1F *hbkg_nopeak = (TH1F *)hfbkg->Clone();
        hbkg_nopeak->SetLineColor(kRed);
        hbkg_nopeak->SetMarkerColor(kRed);
        hbkg_nopeak->SetFillColor(kRed);
        hbkg_nopeak->SetFillStyle(3001);
        for (int ib = 0; ib < hbkg_nopeak->GetNbinsX(); ib++)
        {
            double bc = hbkg_nopeak->GetBinCenter(ib + 1);
            if (bc < kNormRangepT[ip][0] || bc > kNormRangepT[ip][1])
                hbkg_nopeak->SetBinContent(ib + 1, -999);
        }
        gPad->SetRightMargin(0.05);
        gPad->SetLeftMargin(0.15);
        gPad->SetTopMargin(0.09);
        gPad->SetBottomMargin(0.12);

        SetHistoStyle(fHistTotal[ip], 1, 8, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);
        SetHistoStyle(hfbkg, kRed, 24, 1.5, 0.05, 0.05, 0.05, 0.05, 1.13, 1.4);

        fHistTotal[ip]->SetMaximum(fHistTotal[ip]->GetMaximum() * 1.15);
        fHistTotal[ip]->SetMarkerSize(0.5);
        hfbkg->SetMarkerSize(0.5);
        fHistTotal[ip]->Draw("E");
        fHistTotal[ip]->GetYaxis()->SetTitle(Form("Counts/%.1f MeV/c^{2}", binwidth_file * 1000));
        fHistTotal[ip]->GetXaxis()->SetTitle("M_{K#pi} (GeV/c^{2})");
        fHistTotal[ip]->GetXaxis()->SetLabelOffset(0.015);
        fHistTotal[ip]->GetYaxis()->SetMaxDigits(3);
        fHistTotal[ip]->SetStats(0);
        hfbkg->Draw("E same");

        TLegend *leg112 = new TLegend(0.18, 0.80, 0.50, 0.893, NULL, "brNDC");
        leg112->AddEntry(fHistTotal[ip], "Sig+bkg", "p");
        SetLegendStyle(leg112);
        leg112->SetTextSize(0.035);
        if (kResBkg == "MIX")
            leg112->AddEntry(hfbkg, "Mixed-event bkg", "p");
        else if (kResBkg == "LIKE")
            leg112->AddEntry(hfbkg, "Like-sign pairs", "p");
        else if (kResBkg == "ROTATED")
            leg112->AddEntry(hfbkg, "Rotated unlike-sign pairs", "p");
        if (kResBkg == "MIX" || kResBkg == "ROTATED")
        {
            hbkg_nopeak->Draw("BAR same");
            leg112->AddEntry(hbkg_nopeak, "Normalisation range", "f");
        }
        leg112->SetNColumns(2);
        leg112->SetColumnSeparation(0.3);
        leg112->Draw();

        TLatex *ltx = new TLatex(0.27, 0.95,
                                 Form("%0.1f < #it{p}_{T}(GeV/#it{c}) < %0.1f", pT_bins[ip], pT_bins[ip + 1]));
        ltx->SetNDC();
        ltx->SetTextFont(22);
        ltx->SetTextSize(0.06);
        ltx->Draw();

        auto c_clone_sig = (TCanvas *)cSigbkg[ip]->Clone(Form("hsigbkg_pt_%d", ip + 1));
        c_sigbkg.push_back(c_clone_sig);
        cSigbkg[ip]->SaveAs(Form((Cenoutputfolder + "/hsigbkg_pt%d." + outputtype).Data(), ip + 1));
        cSigbkg[ip]->Close();

        delete fSig;
        delete fCorr;
        delete fBkg;
        delete fTotalBkg;
        delete fTotal;

    } // =================== end pT loop ======================================

    // =========================================================================
    // Optional multi-panel summary canvases
    // =========================================================================
    if (plot_a4)
    {
        int A4w = 1000, A4h = 1414;
        int col = static_cast<int>(ceil(Npt / 7.0));
        auto *cAllSig = new TCanvas("cAll_sigbkg", "All sig+bkg", A4w, A4h);
        cAllSig->Divide(col, 7);
        for (int i = 0; i < (int)c_sigbkg.size(); i++)
        {
            cAllSig->cd(i + 1);
            c_sigbkg[i]->DrawClonePad();
        }
        cAllSig->SaveAs((Cenoutputfolder + "/hsigbkg_pt_group_all." + outputtype).Data());

        auto *cAllFit = new TCanvas("cAll_fitsig", "All fit", A4w, A4h);
        cAllFit->Divide(col, 7);
        for (int i = 0; i < (int)c_fitsig.size(); i++)
        {
            cAllFit->cd(i + 1);
            c_fitsig[i]->DrawClonePad();
        }
        cAllFit->SaveAs((Cenoutputfolder + "/hfitsig_pt_group_all." + outputtype).Data());

        auto *csigbyrefl = new TCanvas("cAll_sigbyrefl", "All sigByrefl", A4w, A4h);
        csigbyrefl->Divide(col, 7);
        for (int i = 0; i < (int)c_sigbyref.size(); i++)
        {
            csigbyrefl->cd(i + 1);
            c_sigbyref[i]->DrawClonePad();
        }
        csigbyrefl->SaveAs((Cenoutputfolder + "/hsigbyref_all." + outputtype).Data());
    }

    if (plot_ppt)
    {
        int col = static_cast<int>(ceil(Npt / 3.0));
        auto *cAllSig = new TCanvas("cAll_sigbkg_ppt", "All sig+bkg PPT", 1800, 1200);
        cAllSig->Divide(col, 3, 0.001, 0.001);
        for (int i = 0; i < (int)c_sigbkg.size(); i++)
        {
            cAllSig->cd(i + 1);
            c_sigbkg[i]->DrawClonePad();
        }
        cAllSig->SaveAs((Cenoutputfolder + "/hsigbkg_pt_group_all." + outputtype).Data());

        auto *cAllFit = new TCanvas("cAll_fitsig_ppt", "All fit PPT", 1800, 1200);
        cAllFit->Divide(col, 3, 0.001, 0.001);
        for (int i = 0; i < (int)c_fitsig.size(); i++)
        {
            cAllFit->cd(i + 1);
            c_fitsig[i]->DrawClonePad();
        }
        cAllFit->SaveAs((Cenoutputfolder + "/hfitsig_pt_group_all." + outputtype).Data());
    }

    // =========================================================================
    // Final result plots
    // =========================================================================
    TFile *fresults = new TFile((Cenoutputfolder + "/results.root").Data(), "RECREATE");
    TCanvas *csig = new TCanvas("csig", "", 720, 720);
    SetCanvasStyle(csig, 0.18, 0.05, 0.08, 0.15);

    SetHistoQA(hsignificance);
    hsignificance->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hsignificance->Draw();
    hsignificance->Write("significance");
    csig->SaveAs((Cenoutputfolder + "/significance." + outputtype).Data());
    csig->Clear();

    hChiSquare->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hChiSquare->GetYaxis()->SetTitle("#chi^{2}/NDF");
    SetHistoQA(hChiSquare);
    hChiSquare->Draw("P");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    hChiSquare->Write("chi_ndf");
    csig->SaveAs((Cenoutputfolder + "/chi." + outputtype).Data());
    csig->Clear();

    hmass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hmass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    SetHistoQA(hmass);
    hmass->SetMaximum(0.91);
    hmass->Draw("pe");
    hmass->Write("mass");
    TLine *lmass = new TLine(hmass->GetXaxis()->GetXmin(), masspdg,
                             hmass->GetXaxis()->GetXmax(), masspdg);
    lmass->SetLineStyle(2);
    lmass->SetLineColor(kRed);
    lmass->SetLineWidth(3);
    lmass->Draw();
    TLegend *massleg = new TLegend(0.65, 0.2, 0.9, 0.3);
    SetLegendStyle(massleg);
    massleg->SetTextSize(txtsize);
    massleg->AddEntry(lmass, "PDG Mass", "l");
    massleg->Draw("l");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    csig->SaveAs((Cenoutputfolder + "/mass." + outputtype).Data());
    csig->Clear();

    hwidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hwidth->GetYaxis()->SetTitle("Width (GeV/c^{2})");
    SetHistoQA(hwidth);
    hwidth->SetMaximum(hwidth->GetMaximum() * 2);
    hwidth->SetMinimum(0);
    hwidth->Draw("pe");
    hwidth->Write("width");
    TLine *lwidth = new TLine(hwidth->GetXaxis()->GetXmin(), widthpdg,
                              hwidth->GetXaxis()->GetXmax(), widthpdg);
    lwidth->SetLineStyle(2);
    lwidth->SetLineColor(kRed);
    lwidth->SetLineWidth(3);
    lwidth->Draw();
    TLegend *wleg = new TLegend(0.2, 0.75, 0.4, 0.85);
    SetLegendStyle(wleg);
    wleg->SetFillStyle(0);
    wleg->SetTextSize(txtsize);
    wleg->AddEntry(lwidth, "PDG Width", "l");
    wleg->Draw();
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    csig->SaveAs((Cenoutputfolder + "/width_pt." + outputtype).Data());
    csig->Clear();

    hsigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hsigma->GetYaxis()->SetTitle("Resolution #sigma_{res} (GeV/c^{2})");
    SetHistoQA(hsigma);
    hsigma->SetMaximum(hsigma->GetMaximum() * 1.5);
    hsigma->SetMinimum(0);
    hsigma->Draw("pe");
    hsigma->Write("sigma");
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    csig->SaveAs((Cenoutputfolder + "/sigma_pt." + outputtype).Data());
    csig->Clear();

    SetHistoQA(hintegral_yield);
    hintegral_yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hintegral_yield->GetYaxis()->SetTitle(
        "1/#it{N}_{Ev} d^{2}#it{N}/(d#it{y} d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    gPad->SetLogy(1);
    hintegral_yield->Scale(0.5);
    hintegral_yield->Draw("pe");
    hintegral_yield->Write("yield_integral");
    csig->SaveAs((Cenoutputfolder + "/yield_integral." + outputtype).Data());
    csig->Clear();

    hYbincount->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hYbincount->GetYaxis()->SetTitle(
        "1/#it{N}_{Ev} d^{2}#it{N}/(d#it{y} d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    SetHistoQA(hYbincount);
    hYbincount->Scale(0.5);
    hYbincount->Draw("pe");
    hYbincount->Write("yield_bincount");
    csig->SaveAs((Cenoutputfolder + "/yield_bincount." + outputtype).Data());
    csig->Clear();

    // =====================================================================
    // Yield Comparison with Ratio Plot
    // =====================================================================
    csig->Clear();

    // --- 1. Upper Pad: Yield Spectra ---
    TPad *padYield1 = new TPad("padYield1", "padYield1", 0.0, 0.3, 1.0, 1.0);
    padYield1->SetBottomMargin(0.0);
    padYield1->SetLeftMargin(0.18);
    padYield1->SetLogy(1);
    padYield1->Draw();
    padYield1->cd();

    hintegral_yield->SetMarkerColor(kRed);
    hintegral_yield->SetLineColor(kRed);
    hYbincount->SetMarkerColor(kBlue);
    hYbincount->SetLineColor(kBlue);
    hintegral_yield->SetStats(0);
    hYbincount->SetStats(0);

    hintegral_yield->GetXaxis()->SetLabelSize(0);
    hintegral_yield->GetXaxis()->SetTitleSize(0);
    hintegral_yield->GetYaxis()->SetTitleSize(0.05);
    hintegral_yield->GetYaxis()->SetLabelSize(0.045);
    hintegral_yield->GetYaxis()->SetTitleOffset(1.5);

    hintegral_yield->Draw("pe");
    hYbincount->Draw("pe same");

    TLegend *legcomp = new TLegend(0.45, 0.75, 0.85, 0.88);
    SetLegendStyle(legcomp);
    legcomp->SetTextSize(0.045);
    legcomp->AddEntry(hintegral_yield, "Integral yield", "ep");
    legcomp->AddEntry(hYbincount, "Bin-count yield", "ep");
    legcomp->Draw();

    // --- 2. Lower Pad: Ratio (Bin-count / Integral) ---
    csig->cd();
    TPad *padYield2 = new TPad("padYield2", "padYield2", 0.0, 0.0, 1.0, 0.3);
    padYield2->SetTopMargin(0.0);
    padYield2->SetBottomMargin(0.35);
    padYield2->SetLeftMargin(0.18);
    padYield2->Draw();
    padYield2->cd();

    TH1D *hYieldRatio = (TH1D *)hYbincount->Clone("hYieldRatio");
    hYieldRatio->Divide(hintegral_yield);
    hYieldRatio->SetTitle("");
    hYieldRatio->GetYaxis()->SetTitle("Bin-count / Integral");
    hYieldRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");

    hYieldRatio->GetXaxis()->SetTitleSize(0.13);
    hYieldRatio->GetXaxis()->SetLabelSize(0.11);
    hYieldRatio->GetXaxis()->SetTitleOffset(1.1);
    hYieldRatio->GetYaxis()->SetTitleSize(0.09);
    hYieldRatio->GetYaxis()->SetLabelSize(0.08);
    hYieldRatio->GetYaxis()->SetTitleOffset(0.8);
    hYieldRatio->GetYaxis()->SetNdivisions(505);

    hYieldRatio->SetMinimum(0.85);
    hYieldRatio->SetMaximum(1.15);
    hYieldRatio->SetMarkerColor(kBlack);
    hYieldRatio->SetLineColor(kBlack);
    hYieldRatio->Draw("pe");

    TLine *lineYR = new TLine(hYieldRatio->GetXaxis()->GetXmin(), 1.0,
                              hYieldRatio->GetXaxis()->GetXmax(), 1.0);
    lineYR->SetLineColor(kRed);
    lineYR->SetLineStyle(2);
    lineYR->Draw("same");

    csig->SaveAs((Cenoutputfolder + "/yield_compare." + outputtype).Data());
    csig->Close();

    fresults->Close();
    fTemplateFile->Close();

    timer.Stop();
    cout << "\n===== Timing =====" << endl;
    cout << "Real time : " << timer.RealTime() << " s" << endl;
    cout << "CPU time  : " << timer.CpuTime() << " s" << endl;
}