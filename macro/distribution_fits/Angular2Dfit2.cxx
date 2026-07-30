#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <complex>
#include <iomanip>

#include "TRandom3.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"
#include "../src/style.h"

struct Wave
{
    int l, m;
    bool operator<(const Wave &o) const
    {
        if (l != o.l)
            return l < o.l;
        return m < o.m;
    }
};

const std::vector<Wave> WAVES = {
    {0, 0}, {2, 0}, {2, 1}, {2, 2}, {4, 0}, {4, 1}, {4, 2}};

const std::map<Wave, double> NLM = {
    {{0, 0}, 1.0}, {{2, 0}, 1.0}, {{2, 1}, 2.0}, {{2, 2}, 2.0}, {{4, 0}, 1.0}, {{4, 1}, 2.0}, {{4, 2}, 2.0}};

double ReY(int l, int m, double ct, double phi)
{
    double st2 = TMath::Max(0.0, 1.0 - ct * ct);
    double st = TMath::Sqrt(st2);
    double pi = TMath::Pi();

    if (l == 0 && m == 0)
        return 1.0 / (2.0 * TMath::Sqrt(pi));
    if (l == 2)
    {
        if (m == 0)
            return 0.25 * TMath::Sqrt(5.0 / pi) * (3.0 * ct * ct - 1.0);
        if (m == 1)
            return 0.5 * TMath::Sqrt(15.0 / pi) * st * ct * TMath::Cos(phi);
        if (m == 2)
            return 0.25 * TMath::Sqrt(15.0 / pi) * st2 * TMath::Cos(2.0 * phi);
    }
    if (l == 4)
    {
        if (m == 0)
            return (3.0 / 16.0) * TMath::Sqrt(1.0 / pi) * (35.0 * std::pow(ct, 4) - 30.0 * ct * ct + 3.0);
        if (m == 1)
            return 0.75 * TMath::Sqrt(5.0 / (2.0 * pi)) * st * (7.0 * std::pow(ct, 3) - 3.0 * ct) * TMath::Cos(phi);
        if (m == 2)
            return (3.0 / 8.0) * TMath::Sqrt(5.0 / pi) * st2 * (7.0 * ct * ct - 1.0) * TMath::Cos(2.0 * phi);
        if (m == 3)
            return 0.75 * TMath::Sqrt(35.0 / (2.0 * pi)) * std::pow(st, 3) * ct * TMath::Cos(3.0 * phi);
        if (m == 4)
            return (3.0 / 16.0) * TMath::Sqrt(35.0 / pi) * std::pow(st, 4) * TMath::Cos(4.0 * phi);
    }
    throw std::runtime_error("Unsupported (l,m) combination");
}

double AB(std::complex<double> A, std::complex<double> B)
{
    return std::real(A * std::conj(B));
}

double I0(double ct, double phi, const std::map<Wave, double> &tvals)
{
    double out = 0.0;
    for (const auto &w : WAVES)
    {
        out += NLM.at(w) * tvals.at(w) * ReY(w.l, w.m, ct, phi);
    }
    return out;
}

double eff(double ct, double phi)
{
    return TMath::Max(0.02, 0.6 + 0.3 * TMath::Cos(phi) - 0.4 * ct * ct);
}

void sample(TRandom3 &rng, std::function<double(double, double)> weight_fn,
            int n_target, std::vector<double> &out_ct, std::vector<double> &out_phi, int batch = 200000)
{
    double wmax = 0.0;
    for (int i = 0; i < 200000; ++i)
    {
        double w = weight_fn(rng.Uniform(-1.0, 1.0), rng.Uniform(-TMath::Pi(), TMath::Pi()));
        if (w > wmax)
            wmax = w;
    }
    wmax *= 1.2;

    int n_have = 0;
    while (n_have < n_target)
    {
        for (int i = 0; i < batch; ++i)
        {
            double c = rng.Uniform(-1.0, 1.0);
            double p = rng.Uniform(-TMath::Pi(), TMath::Pi());
            if (rng.Uniform(0.0, wmax) < weight_fn(c, p))
            {
                out_ct.push_back(c);
                out_phi.push_back(p);
                n_have++;
                if (n_have == n_target)
                    break;
            }
        }
    }
}

struct ObjectiveData
{
    std::vector<double> Ncounts;
    std::vector<std::vector<double>> ReY_grid;
    std::vector<double> eps_lm;
};

static ObjectiveData gFitData;

double negLogL(const double *pars)
{
    double nbar = 0.0;
    for (size_t i = 0; i < WAVES.size(); ++i)
    {
        nbar += NLM.at(WAVES[i]) * gFitData.eps_lm[i] * pars[i];
    }

    double logL_sum = 0.0;
    size_t n_bins = gFitData.Ncounts.size();

    for (size_t bin = 0; bin < n_bins; ++bin)
    {
        double counts = gFitData.Ncounts[bin];
        if (counts <= 0)
            continue;

        double model = 0.0;
        for (size_t i = 0; i < WAVES.size(); ++i)
        {
            model += NLM.at(WAVES[i]) * pars[i] * gFitData.ReY_grid[i][bin];
        }

        if (model <= 0.0)
            return 1e12;

        logL_sum += counts * TMath::Log(model);
    }

    return -logL_sum + nbar;
}

void Angular2Dfit2()
{
    TRandom3 rng(1);

    // 1) True Physics Amplitudes
    std::complex<double> S0(5.0, 0.0);
    std::complex<double> D0 = 3.0 * std::exp(std::complex<double>(0, 0.5));
    std::complex<double> D1m = 2.0 * std::exp(std::complex<double>(0, 1.0));
    std::complex<double> D1p = 1.5 * std::exp(std::complex<double>(0, -0.3));

    double sqrt4pi = TMath::Sqrt(4.0 * TMath::Pi());
    std::map<Wave, double> t_true;

    t_true[{0, 0}] = (std::norm(S0) + std::norm(D0) + std::norm(D1m) + std::norm(D1p)) / sqrt4pi;
    t_true[{2, 0}] = (1.0 / 35.0) * (TMath::Sqrt(5.0) * (10.0 * std::norm(D0) + 5.0 * std::norm(D1m) + 5.0 * std::norm(D1p)) + 70.0 * AB(D0, S0)) / sqrt4pi;
    t_true[{2, 1}] = (TMath::Sqrt(2.0) / 35.0) * (35.0 * AB(S0, D1m) + 5.0 * TMath::Sqrt(5.0) * AB(D0, D1m)) / sqrt4pi;
    t_true[{2, 2}] = (TMath::Sqrt(15.0) / (35.0 * TMath::Sqrt(2.0))) * (5.0 * std::norm(D1m) - 5.0 * std::norm(D1p)) / sqrt4pi;
    t_true[{4, 0}] = (1.0 / 7.0) * (6.0 * std::norm(D0) - 4.0 * std::norm(D1m) - 4.0 * std::norm(D1p)) / sqrt4pi;
    t_true[{4, 1}] = (2.0 * TMath::Sqrt(15.0) / 7.0) * AB(D0, D1m) / sqrt4pi;
    t_true[{4, 2}] = (TMath::Sqrt(10.0) / 7.0) * (std::norm(D1m) - std::norm(D1p)) / sqrt4pi;

    // 2) Generate Toy Data
    int n_data = 20000;
    std::vector<double> data_ct, data_phi;
    sample(rng, [&](double c, double p)
           { return TMath::Max(0.0, I0(c, p, t_true)) * eff(c, p); }, n_data, data_ct, data_phi);

    // 3) Generate MC events
    int n_gen = 500000;
    std::vector<double> mc_ct, mc_phi;
    sample(rng, [&](double c, double p)
           { return eff(c, p); }, n_gen, mc_ct, mc_phi);

    double total_eff = 0.0;
    int n_eff_samples = 200000;
    for (int i = 0; i < n_eff_samples; ++i)
    {
        total_eff += eff(rng.Uniform(-1.0, 1.0), rng.Uniform(-TMath::Pi(), TMath::Pi()));
    }
    double mean_eff = total_eff / n_eff_samples;
    double Ngen_real = n_gen / mean_eff;

    gFitData.eps_lm.resize(WAVES.size(), 0.0);
    for (size_t i = 0; i < WAVES.size(); ++i)
    {
        double sum_rey = 0.0;
        for (size_t k = 0; k < mc_ct.size(); ++k)
        {
            sum_rey += ReY(WAVES[i].l, WAVES[i].m, mc_ct[k], mc_phi[k]);
        }
        gFitData.eps_lm[i] = (4.0 * TMath::Pi() / Ngen_real) * sum_rey;
    }

    // 4) Fill 2D Histogram Grid
    int nbin_ct = 20, nbin_phi = 20;
    TH2D *h_data = new TH2D("h_data", "Toy Data (2D);cos(#theta);#phi", nbin_ct, -1.0, 1.0, nbin_phi, -TMath::Pi(), TMath::Pi());
    for (size_t i = 0; i < data_ct.size(); ++i)
    {
        h_data->Fill(data_ct[i], data_phi[i]);
    }

    gFitData.ReY_grid.resize(WAVES.size());

    for (int ix = 1; ix <= nbin_ct; ++ix)
    {
        double ct_center = h_data->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= nbin_phi; ++iy)
        {
            double phi_center = h_data->GetYaxis()->GetBinCenter(iy);

            gFitData.Ncounts.push_back(h_data->GetBinContent(ix, iy));

            for (size_t i = 0; i < WAVES.size(); ++i)
            {
                gFitData.ReY_grid[i].push_back(ReY(WAVES[i].l, WAVES[i].m, ct_center, phi_center));
            }
        }
    }

    // 5) Minimize Likelihood via Migrad
    ROOT::Fit::Fitter fitter;
    ROOT::Math::Functor fcn(&negLogL, WAVES.size());

    std::vector<double> x0(WAVES.size());
    for (size_t i = 0; i < WAVES.size(); ++i)
    {
        x0[i] = t_true[WAVES[i]] * 0.5;
    }

    fitter.SetFCN(fcn, x0.data());

    for (size_t i = 0; i < WAVES.size(); ++i)
    {
        TString name = Form("t_%d%d", WAVES[i].l, WAVES[i].m);
        fitter.Config().ParSettings(i).SetName(name.Data());
        fitter.Config().ParSettings(i).SetValue(x0[i]);
        fitter.Config().ParSettings(i).SetStepSize(0.05);
    }

    fitter.Config().SetMinimizer("Minuit2", "Migrad");
    fitter.FitFCN();

    const ROOT::Fit::FitResult &result = fitter.Result();
    std::cout << "\nFit converged: " << (result.IsValid() ? "True" : "False") << "\n";

    std::cout << std::setw(10) << "wave"
              << std::setw(12) << "true"
              << std::setw(12) << "fitted"
              << std::setw(12) << "err" << "\n";
    std::cout << "--------------------------------------------------\n";

    std::map<Wave, double> t_fit;
    for (size_t i = 0; i < WAVES.size(); ++i)
    {
        t_fit[WAVES[i]] = result.Value(i);
        TString label = Form("(%d,%d)", WAVES[i].l, WAVES[i].m);
        std::cout << std::setw(10) << label
                  << std::setw(12) << std::fixed << std::setprecision(2) << t_true[WAVES[i]]
                  << std::setw(12) << result.Value(i)
                  << std::setw(12) << result.Error(i) << "\n";
    }

    // 6) Plotting
    double bin_area = (2.0 / nbin_ct) * (2.0 * TMath::Pi() / nbin_phi);
    TH2D *h_fit = new TH2D("h_fit", "Fitted Model Expected Counts;cos(#theta);#phi", nbin_ct, -1.0, 1.0, nbin_phi, -TMath::Pi(), TMath::Pi());

    for (int ix = 1; ix <= nbin_ct; ++ix)
    {
        double c = h_fit->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= nbin_phi; ++iy)
        {
            double p = h_fit->GetYaxis()->GetBinCenter(iy);
            double val = I0(c, p, t_fit) * eff(c, p) * bin_area;
            h_fit->SetBinContent(ix, iy, val);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Angular Fit Verification (with 3D Lego Overlay)", 1440, 1080);
    c1->Divide(2, 2);

    // Panel 1: 2D Data Heatmap
    c1->cd(1);
    SetHistoQA(h_data);
    h_data->SetStats(0);
    h_data->Draw("COLZ");

    // Panel 2: 2D Model Heatmap
    c1->cd(2);
    SetHistoQA(h_fit);
    h_fit->SetStats(0);
    h_fit->Draw("COLZ");

    // Panel 3: 1D Projection along cos(theta)
    c1->cd(3);
    TH1D *projX_data = h_data->ProjectionX("projX_data");
    TH1D *projX_fit = h_fit->ProjectionX("projX_fit");
    SetHistoQA(projX_data);
    SetHistoQA(projX_fit);

    projX_data->SetMarkerStyle(20);
    projX_data->SetMarkerSize(0.8);
    projX_data->SetLineColor(kBlack);
    projX_data->SetTitle("Projection along cos(#theta);cos(#theta);Events");
    projX_data->SetStats(0);
    projX_data->Draw("E1");

    projX_fit->SetLineColor(kRed);
    projX_fit->SetLineWidth(2);
    projX_fit->Draw("HIST SAME");

    TLegend *leg1 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg1->AddEntry(projX_data, "Data", "pe");
    leg1->AddEntry(projX_fit, "Fit Model", "l");
    leg1->Draw();

    // Panel 4: 1D Projection along phi
    c1->cd(4);
    TH1D *projY_data = h_data->ProjectionY("projY_data");
    TH1D *projY_fit = h_fit->ProjectionY("projY_fit");

    SetHistoQA(projY_data);
    SetHistoQA(projY_fit);

    projY_data->SetMarkerStyle(20);
    projY_data->SetMarkerSize(0.8);
    projY_data->SetLineColor(kBlack);
    projY_data->SetTitle("Projection along #phi;#phi;Events");
    projY_data->SetStats(0);
    projY_data->Draw("E1");

    projY_fit->SetLineColor(kRed);
    projY_fit->SetLineWidth(2);
    projY_fit->Draw("HIST SAME");

    TLegend *leg2 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg2->AddEntry(projY_data, "Data", "pe");
    leg2->AddEntry(projY_fit, "Fit Model", "l");
    leg2->Draw();

    c1->Update();

    TCanvas *c2 = new TCanvas("c2", "Angular Fit Verification (3D Lego Overlay)", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.03, 0.03, 0.15);
    TH2D *h_lego_data = (TH2D *)h_data->Clone("h_lego_data");
    SetHistoQA2D(h_lego_data);
    h_lego_data->SetTitle("3D Lego (Data) vs Fit Surface;cos(#theta);#phi;Events");
    h_lego_data->SetLineColor(kBlue + 1);
    h_lego_data->SetFillColor(kBlue - 9);
    h_lego_data->SetStats(0);
    h_lego_data->Draw("LEGO2");

    TH2D *h_lego_fit = (TH2D *)h_fit->Clone("h_lego_fit");
    SetHistoQA2D(h_lego_fit);
    h_lego_fit->SetLineColor(kRed);
    h_lego_fit->SetFillColorAlpha(kRed, 0.35); // Semi-transparent fitted surface
    h_lego_fit->Draw("SURF SAME");

    TLegend *leg_lego = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg_lego->AddEntry(h_lego_data, "Data Binned Counts", "f");
    leg_lego->AddEntry(h_lego_fit, "Fitted Model Surface", "f");
    leg_lego->Draw();
}