// DirectXRobustnessStudy_BW_expol3.C
//
// Standalone DIRECT X -> phi phi robustness study using:
//     Breit-Wigner + smooth background
//
// No template fitting is used.
//
// Keep this file in the same directory as:
//     DoublePhi_analysis_Voigtian.C
//
// Example:
// root -l -b -q 'DirectXRobustnessStudy_BW_expol3.C+("AnalysisResults.root","doublephimeson/SEMassDoublePhi",6,100,0.004,"XStudyBW_pt6")'
// root -l -b -q 'DirectXRobustnessStudy_BW_expol3.C+("AnalysisResults.root","doublephimeson/SEMassDoublePhi",9,100,0.004,"XStudyBW_pt9")'

#include "DoublePhi_analysis_Voigtian.C"

#include <iomanip>

// ============================================================================
// BW model used ONLY by this robustness study.
// p[0] = total signal yield N_X
// p[1] = full Breit-Wigner width Gamma_X
// p[2] = pole mass M_X
// p[3..6] = background parameters
// ============================================================================

namespace
{

    double gBWStudyBinWidth = 0.01;
    double gBWStudyFitMin = 2.5;
    double gBWStudyFitMax = 2.9;
    int gBWStudyBackgroundModel = kChebyshev3;

    double BWStudySignal(double *x, double *p)
    {
        const double nX = std::max(0.0, p[0]);
        const double gamma = std::max(1.0e-12, p[1]);
        const double mean = p[2];

        // TMath::BreitWigner uses the full width Gamma.
        // Multiplication by bin width makes N_X a yield-like normalization.
        return gBWStudyBinWidth * nX *
               TMath::BreitWigner(x[0], mean, gamma);
    }

    double BWStudyBackground(double *x, double *p)
    {
        return BackgroundValue(
            gBWStudyBackgroundModel,
            x[0],
            gBWStudyFitMin,
            gBWStudyFitMax,
            p);
    }

    double BWStudyTotal(double *x, double *p)
    {
        return BWStudySignal(x, p) + BWStudyBackground(x, p + 3);
    }

    // Fraction of a normalized non-relativistic Breit-Wigner inside [xmin,xmax].
    double BWIntegralFraction(double xmin,
                              double xmax,
                              double mean,
                              double gamma)
    {
        if (!(xmax > xmin) || !(gamma > 0.0))
            return 0.0;

        const double a = 2.0 * (xmin - mean) / gamma;
        const double b = 2.0 * (xmax - mean) / gamma;

        return (std::atan(b) - std::atan(a)) / TMath::Pi();
    }

    struct DirectFitSnapshotBW
    {
        bool ok = false;
        int status = -999;
        int covStatus = -999;

        double chi2 = 0.0;
        int ndf = 0;

        double nSig = 0.0;
        double nSigErr = 0.0;

        double gamma = 0.0;
        double gammaErr = 0.0;

        double mean = 0.0;
        double meanErr = 0.0;

        double corrNGamma = 0.0;
        double maxAbsCorrGammaBkg = 0.0;
        double maxAbsCorrNBkg = 0.0;

        std::vector<double> pars;
        std::vector<double> errs;
    };

    TF1 *MakeDirectBWFitFunction(TH1D *hMass,
                                 const char *name,
                                 double fitMin,
                                 double fitMax,
                                 double meanMin,
                                 double meanMax,
                                 double gammaMin,
                                 double gammaMax,
                                 int directBkgModel,
                                 const std::vector<double> *seed = nullptr,
                                 double gammaStart = -1.0,
                                 double meanStart = -1.0)
    {
        if (!hMass)
            return nullptr;

        gBWStudyBinWidth = hMass->GetXaxis()->GetBinWidth(1);
        gBWStudyFitMin = fitMin;
        gBWStudyFitMax = fitMax;
        gBWStudyBackgroundModel = directBkgModel;

        const int firstBin = hMass->FindBin(fitMin + 1.e-9);
        const int lastBin = hMass->FindBin(fitMax - 1.e-9);

        const double nTot = hMass->Integral(firstBin, lastBin);
        const int nBins = std::max(1, lastBin - firstBin + 1);
        const double average = std::max(1.0, nTot / nBins);

        const double meanInit =
            (meanStart > meanMin && meanStart < meanMax)
                ? meanStart
                : 0.5 * (meanMin + meanMax);

        const double gammaInit =
            (gammaStart > gammaMin && gammaStart < gammaMax)
                ? gammaStart
                : 0.5 * (gammaMin + gammaMax);

        const double fracInFit = std::max(
            1.0e-3,
            BWIntegralFraction(fitMin, fitMax, meanInit, gammaInit));

        const double signalYieldInit = 0.03 * nTot / fracInFit;

        auto *f = new TF1(
            name,
            BWStudyTotal,
            fitMin,
            fitMax,
            7);

        f->SetNpx(2000);

        f->SetParameter(0, signalYieldInit);
        f->SetParameter(1, gammaInit);
        f->SetParameter(2, meanInit);

        f->SetParName(0, "N_{X}");
        f->SetParName(1, "#Gamma_{X}");
        f->SetParName(2, "M_{X}");

        f->SetParLimits(
            0,
            0.0,
            std::max(10.0, 10.0 * nTot));

        f->SetParLimits(
            1,
            gammaMin,
            gammaMax);

        f->SetParLimits(
            2,
            meanMin,
            meanMax);

        // Background now starts at parameter 3 because there is no sigma_res.
        ConfigureBackgroundParameters(
            f,
            3,
            directBkgModel,
            average,
            HistogramMaximum(hMass));

        if (seed && seed->size() == 7)
        {
            for (int ip = 0; ip < 7; ++ip)
                f->SetParameter(ip, (*seed)[ip]);
        }

        return f;
    }

    DirectFitSnapshotBW FitDirectBWOneStart(
        TH1D *hMass,
        double fitMin,
        double fitMax,
        double meanMin,
        double meanMax,
        double gammaMin,
        double gammaMax,
        int directBkgModel,
        double fixedGamma = -1.0,
        const std::vector<double> *seed = nullptr,
        double gammaStart = -1.0,
        double meanStart = -1.0)
    {
        DirectFitSnapshotBW out;

        if (!hMass || hMass->Integral() <= 0.0)
            return out;

        static int counter = 0;

        TF1 *f = MakeDirectBWFitFunction(
            hMass,
            Form("fDirectBWStudy_%d", counter++),
            fitMin,
            fitMax,
            meanMin,
            meanMax,
            gammaMin,
            gammaMax,
            directBkgModel,
            seed,
            gammaStart,
            meanStart);

        if (!f)
            return out;

        if (fixedGamma > 0.0)
            f->FixParameter(1, fixedGamma);

        TVirtualFitter::SetDefaultFitter("Minuit");
        TVirtualFitter::SetMaxIterations(100000);
        ROOT::Math::MinimizerOptions::SetDefaultMinimizer(
            "Minuit",
            "Migrad");

        // First fit gives a settled starting point.
        hMass->Fit(f, "RQN0");

        // Second fit stores full result/covariance.
        TFitResultPtr r = hMass->Fit(f, "RQSN0");

        out.status =
            r.Get() ? r->Status() : static_cast<int>(r);

        out.covStatus =
            r.Get() ? r->CovMatrixStatus() : -999;

        out.ok =
            r.Get() &&
            out.status == 0;

        out.chi2 = f->GetChisquare();
        out.ndf = f->GetNDF();

        out.nSig = f->GetParameter(0);
        out.nSigErr = f->GetParError(0);

        out.gamma = f->GetParameter(1);
        out.gammaErr = f->GetParError(1);

        out.mean = f->GetParameter(2);
        out.meanErr = f->GetParError(2);

        out.pars.resize(7);
        out.errs.resize(7);

        for (int ip = 0; ip < 7; ++ip)
        {
            out.pars[ip] = f->GetParameter(ip);
            out.errs[ip] = f->GetParError(ip);
        }

        if (r.Get() &&
            r->CovMatrixStatus() > 0 &&
            fixedGamma <= 0.0)
        {

            const TMatrixDSym corr =
                r->GetCorrelationMatrix();

            out.corrNGamma =
                corr(0, 1);

            // Background parameters are 3..6.
            for (int ip = 3; ip < 7; ++ip)
            {

                out.maxAbsCorrGammaBkg =
                    std::max(
                        out.maxAbsCorrGammaBkg,
                        std::abs(corr(1, ip)));

                out.maxAbsCorrNBkg =
                    std::max(
                        out.maxAbsCorrNBkg,
                        std::abs(corr(0, ip)));
            }
        }

        delete f;
        return out;
    }

    // ============================================================================
    // Robust free fit wrapper.
    //
    // For the unconstrained nominal fit we try several physically reasonable
    // starting values in M_X and Gamma_X and keep the converged solution with the
    // smallest chi2.  This prevents the result from depending on the single
    // midpoint initialization.
    //
    // For a profile point (fixedGamma > 0), one refit is performed from the
    // nominal-fit seed because Gamma_X is fixed by construction.
    // ============================================================================
    DirectFitSnapshotBW FitDirectBWNoDraw(
        TH1D *hMass,
        double fitMin,
        double fitMax,
        double meanMin,
        double meanMax,
        double gammaMin,
        double gammaMax,
        int directBkgModel,
        double fixedGamma = -1.0,
        const std::vector<double> *seed = nullptr)
    {
        // Profile point: Gamma is fixed, so start from the supplied nominal seed.
        if (fixedGamma > 0.0 || seed)
        {
            return FitDirectBWOneStart(
                hMass,
                fitMin,
                fitMax,
                meanMin,
                meanMax,
                gammaMin,
                gammaMax,
                directBkgModel,
                fixedGamma,
                seed,
                fixedGamma > 0.0 ? fixedGamma : -1.0,
                -1.0);
        }

        DirectFitSnapshotBW best;
        bool haveBest = false;

        std::vector<double> gammaStarts = {
            0.015, 0.025, 0.030, 0.040, 0.045,
            0.5 * (gammaMin + gammaMax)};

        std::vector<double> meanStarts = {
            0.5 * (meanMin + meanMax),
            2.67, 2.69, 2.71};

        for (double g0 : gammaStarts)
        {
            if (!(g0 > gammaMin && g0 < gammaMax))
                continue;

            for (double m0 : meanStarts)
            {
                if (!(m0 > meanMin && m0 < meanMax))
                    continue;

                DirectFitSnapshotBW r = FitDirectBWOneStart(
                    hMass,
                    fitMin,
                    fitMax,
                    meanMin,
                    meanMax,
                    gammaMin,
                    gammaMax,
                    directBkgModel,
                    -1.0,
                    nullptr,
                    g0,
                    m0);

                if (!r.ok || !std::isfinite(r.chi2))
                    continue;

                if (!haveBest || r.chi2 < best.chi2)
                {
                    best = r;
                    haveBest = true;
                }
            }
        }

        if (!haveBest)
        {
            // Last-resort midpoint fit so the caller still gets diagnostic output.
            best = FitDirectBWOneStart(
                hMass,
                fitMin,
                fitMax,
                meanMin,
                meanMax,
                gammaMin,
                gammaMax,
                directBkgModel,
                -1.0,
                nullptr,
                0.5 * (gammaMin + gammaMax),
                0.5 * (meanMin + meanMax));
        }

        // Explicit boundary warning: covariance/correlations are not trustworthy
        // when Gamma_X is sitting on a hard fit limit.
        const double tol = 1.0e-4 * (gammaMax - gammaMin);
        if (best.ok &&
            (std::abs(best.gamma - gammaMin) < tol ||
             std::abs(best.gamma - gammaMax) < tol))
        {
            std::cout
                << "WARNING: best-fit Gamma_X = "
                << 1000.0 * best.gamma
                << " MeV/c2 is on a fit boundary ["
                << 1000.0 * gammaMin << ", "
                << 1000.0 * gammaMax
                << "] MeV/c2. Do NOT interpret the Hessian error or correlations "
                << "as a normal interior solution. Inspect the profile chi2 scan.\n";
        }

        return best;
    }

    // ============================================================================
    // PROFILE chi2(Gamma_X)
    //
    // For each fixed Gamma_X, refit:
    //   N_X, M_X, and all background parameters.
    // ============================================================================
    void DrawDirectBWGammaProfile(
        TH1D *hMass,
        double fitMin,
        double fitMax,
        double meanMin,
        double meanMax,
        double gammaMin,
        double gammaMax,
        int directBkgModel,
        const char *outPrefix,
        int nScan = 80)
    {
        if (!hMass || nScan < 5)
            return;

        const DirectFitSnapshotBW nominal =
            FitDirectBWNoDraw(
                hMass,
                fitMin,
                fitMax,
                meanMin,
                meanMax,
                gammaMin,
                gammaMax,
                directBkgModel);

        if (!nominal.ok)
        {
            std::cerr
                << "ERROR: nominal direct BW fit failed; "
                << "Gamma profile skipped.\n";
            return;
        }

        std::vector<double> scanGamma;
        std::vector<double> scanChi2;
        std::vector<double> scanNX;
        std::vector<double> scanMX;
        std::vector<int> scanStatus;
        std::vector<int> scanCovStatus;

        double minProfileChi2 = nominal.chi2;

        for (int i = 0; i < nScan; ++i)
        {

            const double gamma =
                gammaMin +
                (gammaMax - gammaMin) *
                    static_cast<double>(i) / (nScan - 1);

            const DirectFitSnapshotBW p =
                FitDirectBWNoDraw(
                    hMass,
                    fitMin,
                    fitMax,
                    meanMin,
                    meanMax,
                    gammaMin,
                    gammaMax,
                    directBkgModel,
                    gamma,
                    &nominal.pars);

            scanGamma.push_back(gamma);

            scanChi2.push_back(
                p.ok
                    ? p.chi2
                    : std::numeric_limits<double>::quiet_NaN());

            scanNX.push_back(p.nSig);
            scanMX.push_back(p.mean);

            scanStatus.push_back(p.status);
            scanCovStatus.push_back(p.covStatus);

            if (p.ok &&
                std::isfinite(p.chi2))
            {

                minProfileChi2 =
                    std::min(
                        minProfileChi2,
                        p.chi2);
            }
        }

        if (minProfileChi2 + 1.e-3 < nominal.chi2)
        {

            std::cout
                << "WARNING: fixed-Gamma scan found chi2 below "
                << "the nominal free fit by "
                << nominal.chi2 - minProfileChi2
                << ". Check minimizer convergence.\n";
        }

        auto *graph = new TGraph();

        graph->SetName(
            "gDirectBWGammaProfile");

        std::ofstream csv(
            Form(
                "%s_directBWGammaProfile_%s.csv",
                outPrefix,
                BackgroundTag(
                    directBkgModel)
                    .Data()));

        csv
            << "gamma_GeV,"
            << "chi2,"
            << "delta_chi2,"
            << "status,"
            << "cov_status,"
            << "NX,"
            << "MX\n";

        int point = 0;

        for (int i = 0; i < nScan; ++i)
        {

            if (!std::isfinite(scanChi2[i]))
            {

                csv
                    << scanGamma[i]
                    << ",nan,nan,"
                    << scanStatus[i] << ","
                    << scanCovStatus[i] << ","
                    << scanNX[i] << ","
                    << scanMX[i] << "\n";

                continue;
            }

            const double deltaChi2 =
                scanChi2[i] - minProfileChi2;

            graph->SetPoint(
                point++,
                1000.0 * scanGamma[i],
                deltaChi2);

            csv
                << scanGamma[i] << ","
                << scanChi2[i] << ","
                << deltaChi2 << ","
                << scanStatus[i] << ","
                << scanCovStatus[i] << ","
                << scanNX[i] << ","
                << scanMX[i] << "\n";
        }

        csv.close();

        auto *c = new TCanvas(
            "cDirectBWGammaProfile",
            "Direct BW + exp(pol3) profile scan of Gamma_X",
            900,
            700);

        c->SetLeftMargin(0.14);
        c->SetBottomMargin(0.13);
        c->SetGridy();

        graph->SetTitle(
            Form(
                "Direct Breit-Wigner profile: %s;"
                "#Gamma_{X} (MeV/#it{c}^{2});"
                "#Delta#chi^{2}",
                BackgroundLabel(
                    directBkgModel)));

        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.8);
        graph->SetLineWidth(2);

        graph->Draw("ALP");

        const double x1 =
            1000.0 * gammaMin;

        const double x2 =
            1000.0 * gammaMax;

        auto *line68 =
            new TLine(
                x1,
                1.0,
                x2,
                1.0);

        line68->SetLineStyle(2);
        line68->Draw("same");

        auto *line95 =
            new TLine(
                x1,
                3.84,
                x2,
                3.84);

        line95->SetLineStyle(3);
        line95->Draw("same");

        TLatex lat;

        lat.SetNDC();
        lat.SetTextSize(0.035);

        lat.DrawLatex(
            0.18,
            0.88,
            Form(
                "free fit: #Gamma_{X}=%.1f #pm %.1f MeV/#it{c}^{2}",
                1000.0 * nominal.gamma,
                1000.0 * nominal.gammaErr));

        lat.DrawLatex(
            0.18,
            0.83,
            Form(
                "#chi^{2}/ndf = %.1f/%d",
                nominal.chi2,
                nominal.ndf));

        lat.DrawLatex(
            0.18,
            0.78,
            "#Delta#chi^{2}=1: approximate 68% interval");

        lat.DrawLatex(
            0.18,
            0.73,
            "#Delta#chi^{2}=3.84: approximate 95% interval");

        c->SaveAs(
            Form(
                "%s_directBWGammaProfile_%s.png",
                outPrefix,
                BackgroundTag(
                    directBkgModel)
                    .Data()));

        c->SaveAs(
            Form(
                "%s_directBWGammaProfile_%s.pdf",
                outPrefix,
                BackgroundTag(
                    directBkgModel)
                    .Data()));

        delete line68;
        delete line95;
        delete c;
        delete graph;
    }

    // ============================================================================
    // BACKGROUND-MODEL STABILITY
    //
    // Same DIRECT Breit-Wigner signal model with each smooth background.
    // ============================================================================
    void DrawDirectBWBackgroundStability(
        TH1D *hMass,
        double fitMin,
        double fitMax,
        double meanMin,
        double meanMax,
        double gammaMin,
        double gammaMax,
        const char *outPrefix)
    {
        if (!hMass)
            return;

        constexpr int nModels = 8;

        auto *hMassStability =
            new TH1D(
                "hDirectBWMassStability",
                ";background model;"
                "M_{X} (GeV/#it{c}^{2})",
                nModels,
                0.5,
                nModels + 0.5);

        auto *hGammaStability =
            new TH1D(
                "hDirectBWGammaStability",
                ";background model;"
                "#Gamma_{X} (MeV/#it{c}^{2})",
                nModels,
                0.5,
                nModels + 0.5);

        auto *hYieldStability =
            new TH1D(
                "hDirectBWYieldStability",
                ";background model;"
                "N_{X}",
                nModels,
                0.5,
                nModels + 0.5);

        std::ofstream csv(
            Form(
                "%s_directBWBackgroundStability.csv",
                outPrefix));

        csv
            << "model_id,"
            << "model,"
            << "status,"
            << "cov_status,"
            << "chi2,"
            << "ndf,"
            << "chi2_ndf,"
            << "NX,"
            << "NX_err,"
            << "Gamma_GeV,"
            << "Gamma_err_GeV,"
            << "MX,"
            << "MX_err,"
            << "corr_NX_Gamma,"
            << "max_abs_corr_Gamma_bkg,"
            << "max_abs_corr_NX_bkg\n";

        for (int model = 0;
             model < nModels;
             ++model)
        {

            const DirectFitSnapshotBW r =
                FitDirectBWNoDraw(
                    hMass,
                    fitMin,
                    fitMax,
                    meanMin,
                    meanMax,
                    gammaMin,
                    gammaMax,
                    model);

            for (TH1D *h :
                 {hMassStability,
                  hGammaStability,
                  hYieldStability})
            {

                h->GetXaxis()->SetBinLabel(
                    model + 1,
                    BackgroundLabel(model));
            }

            if (r.ok)
            {

                hMassStability->SetBinContent(
                    model + 1,
                    r.mean);

                hMassStability->SetBinError(
                    model + 1,
                    r.meanErr);

                hGammaStability->SetBinContent(
                    model + 1,
                    1000.0 * r.gamma);

                hGammaStability->SetBinError(
                    model + 1,
                    1000.0 * r.gammaErr);

                hYieldStability->SetBinContent(
                    model + 1,
                    r.nSig);

                hYieldStability->SetBinError(
                    model + 1,
                    r.nSigErr);
            }

            csv
                << model << ","
                << BackgroundLabel(model) << ","
                << r.status << ","
                << r.covStatus << ","
                << r.chi2 << ","
                << r.ndf << ","
                << (r.ndf > 0
                        ? r.chi2 / r.ndf
                        : 0.0)
                << ","
                << r.nSig << ","
                << r.nSigErr << ","
                << r.gamma << ","
                << r.gammaErr << ","
                << r.mean << ","
                << r.meanErr << ","
                << r.corrNGamma << ","
                << r.maxAbsCorrGammaBkg << ","
                << r.maxAbsCorrNBkg
                << "\n";

            std::cout
                << "DIRECT BW STABILITY: "
                << BackgroundLabel(model)
                << "  status=" << r.status
                << "  cov=" << r.covStatus
                << "  chi2/ndf="
                << (r.ndf > 0
                        ? r.chi2 / r.ndf
                        : 0.0)
                << "  M=" << r.mean
                << "  Gamma="
                << 1000.0 * r.gamma
                << " +/- "
                << 1000.0 * r.gammaErr
                << " MeV"
                << "  NX="
                << r.nSig
                << " +/- "
                << r.nSigErr
                << "  rho(N,Gamma)="
                << r.corrNGamma
                << "  max|rho(Gamma,bkg)|="
                << r.maxAbsCorrGammaBkg
                << "  max|rho(N,bkg)|="
                << r.maxAbsCorrNBkg
                << "\n";
        }

        csv.close();

        auto *c = new TCanvas(
            "cDirectBWBackgroundStability",
            "Direct BW background stability",
            1500,
            1200);

        c->Divide(
            1,
            3,
            0.001,
            0.001);

        c->cd(1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.20);
        gPad->SetGridy();
        hMassStability->SetMarkerStyle(20);
        hMassStability->Draw("E1");

        c->cd(2);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.20);
        gPad->SetGridy();
        hGammaStability->SetMarkerStyle(20);
        hGammaStability->Draw("E1");

        c->cd(3);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.20);
        gPad->SetGridy();
        hYieldStability->SetMarkerStyle(20);
        hYieldStability->Draw("E1");

        c->SaveAs(
            Form(
                "%s_directBWBackgroundStability.png",
                outPrefix));

        c->SaveAs(
            Form(
                "%s_directBWBackgroundStability.pdf",
                outPrefix));

        delete c;
        delete hMassStability;
        delete hGammaStability;
        delete hYieldStability;
    }

} // namespace

// ============================================================================
// Public standalone entry point
// ============================================================================
void DirectXRobustnessStudy_BW_expol3(
    const char *inputFile = "../AnalysisResults_pid2003.root",
    const char *sparseName = "doublephimeson/SEMassDoublePhi",
    double ptMin = 9.0,
    double ptMax = 100.0,
    double deltaMMax = 0.005,
    const char *outPrefix = "DirectXStudyBW",
    int nominalDirectBkgModel = 3, // 3 = exp(pol3)
    double phiMassMin = 1.00,
    double phiMassMax = 1.04,
    double rapidityMin = 0.0,
    double rapidityMax = 0.8,
    double mPairMin = 2.5,
    double mPairMax = 2.9,
    double pairBinWidth = 0.010,
    double xMassMin = 2.63,
    double xMassMax = 2.75,
    double gammaMin = 0.001,
    double gammaMax = 0.20)
{
    // Force the direct background model to exp(pol3).
    nominalDirectBkgModel = 3;

    std::unique_ptr<TFile> f(
        TFile::Open(
            inputFile,
            "READ"));

    if (!f ||
        f->IsZombie())
    {

        std::cerr
            << "ERROR: cannot open "
            << inputFile
            << "\n";

        return;
    }

    THnSparseF *hSparse =
        dynamic_cast<THnSparseF *>(
            f->Get(
                sparseName));

    if (!hSparse)
    {

        std::cerr
            << "ERROR: cannot find THnSparseF "
            << sparseName
            << " in "
            << inputFile
            << "\n";

        return;
    }

    if (!(ptMax > ptMin) ||
        !(phiMassMax > phiMassMin) ||
        !(rapidityMax > rapidityMin) ||
        !(mPairMax > mPairMin) ||
        !(pairBinWidth > 0.0))
    {

        std::cerr
            << "ERROR: invalid analysis ranges.\n";

        return;
    }

    gRapidityMin = rapidityMin;
    gRapidityMax = rapidityMax;

    const int nPairBins =
        ExactNBinsFromWidth(
            mPairMin,
            mPairMax,
            pairBinWidth,
            "pairBinWidth");

    if (nPairBins <= 0)
        return;

    // DIRECT spectrum only.
    // No template construction and no subtraction.
    std::unique_ptr<TH1D> hMass(
        ProjectPairMassNative(
            hSparse,
            Form(
                "hDirectBWMass_pt%.1f_%.1f",
                ptMin,
                ptMax),
            ptMin,
            ptMax,
            phiMassMin,
            phiMassMax,
            0.0,
            deltaMMax,
            mPairMin,
            mPairMax,
            nPairBins));

    if (!hMass ||
        hMass->Integral() <= 0.0)
    {

        std::cerr
            << "ERROR: direct M(phi phi) histogram is empty.\n";

        return;
    }

    std::cout
        << "\n============================================================\n"
        << " DIRECT X ROBUSTNESS STUDY -- BREIT-WIGNER + exp(pol3)\n"
        << " NO TEMPLATE FITTING\n"
        << "============================================================\n"
        << "input       : "
        << inputFile << "\n"
        << "sparse      : "
        << sparseName << "\n"
        << "pT range    : ["
        << ptMin << ", "
        << ptMax
        << "] GeV/c\n"
        << "rapidity    : ["
        << rapidityMin << ", "
        << rapidityMax
        << ")\n"
        << "phi masses  : ["
        << phiMassMin << ", "
        << phiMassMax
        << "] GeV/c2\n"
        << "DeltaM      : < "
        << deltaMMax
        << " GeV/c2\n"
        << "Mphiphi fit : ["
        << mPairMin << ", "
        << mPairMax
        << "] GeV/c2\n"
        << "Gamma range : ["
        << 1000.0 * gammaMin
        << ", "
        << 1000.0 * gammaMax
        << "] MeV/c2\n"
        << "nominal bkg : "
        << BackgroundLabel(
               nominalDirectBkgModel)
        << "\n"
        << "============================================================\n";

    // Save exact direct spectrum.
    {
        TFile fout(
            Form(
                "%s_directBWStudy.root",
                outPrefix),
            "RECREATE");

        hMass->Write(
            "hMphiphi_direct_BW");

        fout.Close();
    }

    // Nominal fit printout.
    const DirectFitSnapshotBW nominal =
        FitDirectBWNoDraw(
            hMass.get(),
            mPairMin,
            mPairMax,
            xMassMin,
            xMassMax,
            gammaMin,
            gammaMax,
            nominalDirectBkgModel);

    std::cout
        << "\nNominal direct BW fit:\n"
        << "status        = "
        << nominal.status << "\n"
        << "cov status    = "
        << nominal.covStatus << "\n"
        << "M_X           = "
        << nominal.mean
        << " +/- "
        << nominal.meanErr
        << " GeV/c2\n"
        << "Gamma_X       = "
        << 1000.0 * nominal.gamma
        << " +/- "
        << 1000.0 * nominal.gammaErr
        << " MeV/c2\n"
        << "N_X           = "
        << nominal.nSig
        << " +/- "
        << nominal.nSigErr
        << "\n"
        << "chi2/ndf      = "
        << nominal.chi2
        << "/"
        << nominal.ndf
        << "\n"
        << "rho(N,Gamma)  = "
        << nominal.corrNGamma
        << "\n"
        << "max |rho(Gamma,bkg)| = "
        << nominal.maxAbsCorrGammaBkg
        << "\n"
        << "max |rho(N,bkg)|     = "
        << nominal.maxAbsCorrNBkg
        << "\n";

    if (std::abs(nominal.gamma - gammaMax) <
        1.0e-4 * (gammaMax - gammaMin))
    {
        std::cout
            << "NOTE: Gamma_X is at the UPPER fit limit. "
            << "The nominal width/error/correlations are boundary-limited; "
            << "use the profile plot before interpreting them.\n";
    }

    // Profile Gamma_X.
    DrawDirectBWGammaProfile(
        hMass.get(),
        mPairMin,
        mPairMax,
        xMassMin,
        xMassMax,
        gammaMin,
        gammaMax,
        nominalDirectBkgModel,
        outPrefix,
        80);

    // Background model stability.
    DrawDirectBWBackgroundStability(
        hMass.get(),
        mPairMin,
        mPairMax,
        xMassMin,
        xMassMax,
        gammaMin,
        gammaMax,
        outPrefix);

    std::cout
        << "\nStudy finished.\n"
        << "Main outputs:\n"
        << "  "
        << outPrefix
        << "_directBWGammaProfile_"
        << BackgroundTag(
               nominalDirectBkgModel)
               .Data()
        << ".pdf\n"
        << "  "
        << outPrefix
        << "_directBWBackgroundStability.pdf\n"
        << "  "
        << outPrefix
        << "_directBWBackgroundStability.csv\n"
        << "  "
        << outPrefix
        << "_directBWStudy.root\n";
}
