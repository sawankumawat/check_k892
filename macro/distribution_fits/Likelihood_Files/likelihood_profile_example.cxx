#include <iostream>
#include <vector>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TRandom3.h>

using namespace std;

void likelihood_profile_example()
{
    cout << "\n================================================" << endl;
    cout << "     Likelihood Profile Plotting Example" << endl;
    cout << "================================================\n"
         << endl;

    // Create example data: Gaussian signal + exponential background
    TH1F *h_data = new TH1F("h_data", "Example Data;Mass (GeV/c^{2});Counts", 100, 0, 10);

    TRandom3 rand(42);

    // Fill with background (exponential)
    for (int i = 0; i < 8000; i++)
    {
        double x = rand.Exp(2.0) + 1.0;
        if (x < 10.0)
            h_data->Fill(x);
    }

    // Fill with signal (Gaussian at 5.0 GeV with width 0.3 GeV)
    for (int i = 0; i < 1000; i++)
    {
        double x = rand.Gaus(5.0, 0.3);
        if (x > 0 && x < 10.0)
            h_data->Fill(x);
    }

    cout << "Generated data with " << h_data->GetEntries() << " entries" << endl;
    cout << "Signal: Gaussian peak at 5.0 GeV" << endl;
    cout << "Background: Exponential distribution" << endl
         << endl;

    // Define fitting function: signal + background
    TF1 *f_total = new TF1("f_total", "[0]*exp(-[1]*x) + [2]*exp(-0.5*((x-[3])/[4])^2)", 1, 9);
    f_total->SetParameters(1000, 0.5, 200, 5.0, 0.3);
    f_total->SetParNames("Bkg_Norm", "Bkg_Slope", "Signal_Norm", "Signal_Mean", "Signal_Width");

    cout << "Fitting full model..." << endl;
    TFitResultPtr fit_result = h_data->Fit("f_total", "RELBMSQ");

    if (fit_result->IsValid())
    {
        cout << "Fit successful!" << endl;
        cout << "Best-fit parameters:" << endl;
        for (int i = 0; i < f_total->GetNpar(); i++)
        {
            cout << "  " << f_total->GetParName(i) << " = "
                 << f_total->GetParameter(i) << " ± " << f_total->GetParError(i) << endl;
        }
        cout << "Chi2/NDF = " << f_total->GetChisquare() << "/" << f_total->GetNDF()
             << " = " << f_total->GetChisquare() / f_total->GetNDF() << endl;
        cout << "-2 log L = " << fit_result->MinFcnValue() << endl
             << endl;
    }

    // Function to create likelihood profile for a given parameter
    auto createProfile = [&](int param_index, const string &param_name, int n_points = 30) -> TGraph *
    {
        cout << "Creating likelihood profile for " << param_name << "..." << endl;

        double best_value = f_total->GetParameter(param_index);
        double param_error = f_total->GetParError(param_index);
        double nll_min = fit_result->MinFcnValue();

        // Define scan range (±3σ around best fit)
        double scan_range = 3.0;
        double param_min = best_value - scan_range * param_error;
        double param_max = best_value + scan_range * param_error;

        vector<double> param_values;
        vector<double> delta_nll_values;

        for (int i = 0; i < n_points; i++)
        {
            double test_value = param_min + i * (param_max - param_min) / (n_points - 1);

            // Fix parameter to test value
            f_total->SetParameter(param_index, test_value);
            f_total->FixParameter(param_index, test_value);

            // Refit with parameter fixed
            TFitResultPtr temp_fit = h_data->Fit("f_total", "RELBMSQ");
            double nll_test = temp_fit->MinFcnValue();

            param_values.push_back(test_value);
            delta_nll_values.push_back(nll_test - nll_min);

            if (i % 5 == 0)
            {
                cout << "  Point " << i + 1 << "/" << n_points
                     << ": " << param_name << " = " << test_value
                     << ", Δ(-2 log L) = " << nll_test - nll_min << endl;
            }
        }

        // Release parameter
        f_total->ReleaseParameter(param_index);
        f_total->SetParameter(param_index, best_value);

        // Create and return TGraph
        TGraph *profile = new TGraph(n_points, &param_values[0], &delta_nll_values[0]);
        profile->SetName(Form("profile_%s", param_name.c_str()));
        profile->SetTitle(Form("Likelihood Profile: %s;%s;#Delta(-2 log L)",
                               param_name.c_str(), param_name.c_str()));
        profile->SetMarkerStyle(20);
        profile->SetMarkerColor(kBlue);
        profile->SetLineColor(kBlue);
        profile->SetLineWidth(2);

        return profile;
    };

    cout << "\n================================================" << endl;
    cout << "Creating likelihood profiles for key parameters..." << endl;
    cout << "================================================" << endl;

    // Create profiles for signal parameters
    TGraph *profile_mean = createProfile(3, "Signal_Mean");
    TGraph *profile_width = createProfile(4, "Signal_Width");
    TGraph *profile_norm = createProfile(2, "Signal_Norm");

    // Create canvas with multiple pads
    TCanvas *c_profiles = new TCanvas("c_profiles", "Likelihood Profiles", 1200, 400);
    c_profiles->Divide(3, 1);

    // Plot signal mean profile
    c_profiles->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    profile_mean->Draw("ALP");

    // Add confidence level lines
    double x_min = profile_mean->GetXaxis()->GetXmin();
    double x_max = profile_mean->GetXaxis()->GetXmax();

    TLine *line1 = new TLine(x_min, 1.0, x_max, 1.0); // 1σ
    line1->SetLineColor(kGreen);
    line1->SetLineStyle(2);
    line1->SetLineWidth(2);
    line1->Draw("same");

    TLine *line2 = new TLine(x_min, 4.0, x_max, 4.0); // 2σ
    line2->SetLineColor(kOrange);
    line2->SetLineStyle(2);
    line2->SetLineWidth(2);
    line2->Draw("same");

    TLine *line3 = new TLine(x_min, 9.0, x_max, 9.0); // 3σ
    line3->SetLineColor(kRed);
    line3->SetLineStyle(2);
    line3->SetLineWidth(2);
    line3->Draw("same");

    profile_mean->GetYaxis()->SetRangeUser(0, 15);

    // Add legend
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(profile_mean, "Profile", "LP");
    leg->AddEntry(line1, "1#sigma (68% CL)", "L");
    leg->AddEntry(line2, "2#sigma (95% CL)", "L");
    leg->AddEntry(line3, "3#sigma (99.7% CL)", "L");
    leg->Draw();

    // Plot signal width profile
    c_profiles->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    profile_width->Draw("ALP");

    x_min = profile_width->GetXaxis()->GetXmin();
    x_max = profile_width->GetXaxis()->GetXmax();

    TLine *line1w = new TLine(x_min, 1.0, x_max, 1.0);
    line1w->SetLineColor(kGreen);
    line1w->SetLineStyle(2);
    line1w->SetLineWidth(2);
    line1w->Draw("same");

    TLine *line2w = new TLine(x_min, 4.0, x_max, 4.0);
    line2w->SetLineColor(kOrange);
    line2w->SetLineStyle(2);
    line2w->SetLineWidth(2);
    line2w->Draw("same");

    TLine *line3w = new TLine(x_min, 9.0, x_max, 9.0);
    line3w->SetLineColor(kRed);
    line3w->SetLineStyle(2);
    line3w->SetLineWidth(2);
    line3w->Draw("same");

    profile_width->GetYaxis()->SetRangeUser(0, 15);

    // Plot signal normalization profile
    c_profiles->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    profile_norm->Draw("ALP");

    x_min = profile_norm->GetXaxis()->GetXmin();
    x_max = profile_norm->GetXaxis()->GetXmax();

    TLine *line1n = new TLine(x_min, 1.0, x_max, 1.0);
    line1n->SetLineColor(kGreen);
    line1n->SetLineStyle(2);
    line1n->SetLineWidth(2);
    line1n->Draw("same");

    TLine *line2n = new TLine(x_min, 4.0, x_max, 4.0);
    line2n->SetLineColor(kOrange);
    line2n->SetLineStyle(2);
    line2n->SetLineWidth(2);
    line2n->Draw("same");

    TLine *line3n = new TLine(x_min, 9.0, x_max, 9.0);
    line3n->SetLineColor(kRed);
    line3n->SetLineStyle(2);
    line3n->SetLineWidth(2);
    line3n->Draw("same");

    profile_norm->GetYaxis()->SetRangeUser(0, 15);

    c_profiles->SaveAs("/home/sawan/check_k892/macro/distribution_fits/likelihood_profiles_example.png");
    cout << "\nLikelihood profiles saved as: likelihood_profiles_example.png" << endl;

    // Create the data fit plot
    TCanvas *c_fit = new TCanvas("c_fit", "Data and Fit", 800, 600);

    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.8);
    h_data->Draw("PE");

    f_total->SetLineColor(kRed);
    f_total->SetLineWidth(2);
    f_total->Draw("same");

    // Draw individual components
    TF1 *f_signal = new TF1("f_signal", "[0]*exp(-0.5*((x-[1])/[2])^2)", 1, 9);
    f_signal->SetParameters(f_total->GetParameter(2), f_total->GetParameter(3), f_total->GetParameter(4));
    f_signal->SetLineColor(kBlue);
    f_signal->SetLineStyle(2);
    f_signal->SetLineWidth(2);
    f_signal->Draw("same");

    TF1 *f_bkg = new TF1("f_bkg", "[0]*exp(-[1]*x)", 1, 9);
    f_bkg->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1));
    f_bkg->SetLineColor(kGreen);
    f_bkg->SetLineStyle(2);
    f_bkg->SetLineWidth(2);
    f_bkg->Draw("same");

    TLegend *leg_fit = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg_fit->SetFillStyle(0);
    leg_fit->SetBorderSize(0);
    leg_fit->AddEntry(h_data, "Data", "PE");
    leg_fit->AddEntry(f_total, "Total Fit", "L");
    leg_fit->AddEntry(f_signal, "Signal", "L");
    leg_fit->AddEntry(f_bkg, "Background", "L");
    leg_fit->Draw();

    // Add fit info
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.DrawLatex(0.15, 0.85, Form("Signal Mean = %.3f ± %.3f GeV", f_total->GetParameter(3), f_total->GetParError(3)));
    lat.DrawLatex(0.15, 0.80, Form("Signal Width = %.3f ± %.3f GeV", f_total->GetParameter(4), f_total->GetParError(4)));
    lat.DrawLatex(0.15, 0.75, Form("#chi^{2}/NDF = %.1f/%d = %.2f", f_total->GetChisquare(), f_total->GetNDF(), f_total->GetChisquare() / f_total->GetNDF()));

    c_fit->SaveAs("/home/sawan/check_k892/macro/distribution_fits/likelihood_fit_example.png");
    cout << "Fit plot saved as: likelihood_fit_example.png" << endl;

    cout << "\n================================================" << endl;
    cout << "How to interpret likelihood profiles:" << endl;
    cout << "================================================" << endl;
    cout << "• Δ(-2 log L) = 1.0 → 68% confidence interval (1σ)" << endl;
    cout << "• Δ(-2 log L) = 4.0 → 95% confidence interval (2σ)" << endl;
    cout << "• Δ(-2 log L) = 9.0 → 99.7% confidence interval (3σ)" << endl;
    cout << "\nThe width of the parabola indicates parameter uncertainty:" << endl;
    cout << "• Narrow profile → well-constrained parameter" << endl;
    cout << "• Wide profile → poorly-constrained parameter" << endl;
    cout << "• Asymmetric profile → non-Gaussian uncertainties" << endl;
    cout << "\nFor your glueball analysis:" << endl;
    cout << "• Mass profiles show how well resonance masses are determined" << endl;
    cout << "• Amplitude profiles show significance of each resonance" << endl;
    cout << "• Width profiles show how well natural widths are constrained" << endl;
    cout << "================================================" << endl;

    // // Clean up
    // delete h_data;
    // delete f_total;
    // delete f_signal;
    // delete f_bkg;
    // delete profile_mean;
    // delete profile_width;
    // delete profile_norm;
}

// Function to demonstrate profile interpretation
void explain_likelihood_profiles()
{
    cout << "\n================================================" << endl;
    cout << "     Likelihood Profile Interpretation Guide" << endl;
    cout << "================================================" << endl;

    cout << "\n1. WHAT IS A LIKELIHOOD PROFILE?" << endl;
    cout << "   - Shows how -2 log L changes as you vary one parameter" << endl;
    cout << "   - All other parameters reoptimized at each point" << endl;
    cout << "   - Minimum occurs at the best-fit value" << endl;

    cout << "\n2. CONFIDENCE INTERVALS:" << endl;
    cout << "   - Δ(-2 log L) = 1.0 → 68% CL (1σ)" << endl;
    cout << "   - Δ(-2 log L) = 4.0 → 95% CL (2σ)" << endl;
    cout << "   - Δ(-2 log L) = 9.0 → 99.7% CL (3σ)" << endl;

    cout << "\n3. SHAPE INTERPRETATION:" << endl;
    cout << "   • Parabolic → Gaussian uncertainties (asymptotic limit)" << endl;
    cout << "   • Asymmetric → Non-Gaussian errors, quote +/- separately" << endl;
    cout << "   • Flat bottom → Parameter not well constrained" << endl;
    cout << "   • Multiple minima → Degeneracies or local minima" << endl;

    cout << "\n4. FOR RESONANCE ANALYSIS:" << endl;
    cout << "   • Mass profiles: How well is the resonance position known?" << endl;
    cout << "   • Width profiles: How well is the natural width determined?" << endl;
    cout << "   • Amplitude profiles: How significant is the resonance?" << endl;

    cout << "\n5. PRACTICAL TIPS:" << endl;
    cout << "   • If profile doesn't reach minimum, expand scan range" << endl;
    cout << "   • If too noisy, increase number of scan points" << endl;
    cout << "   • If asymmetric, report asymmetric uncertainties" << endl;
    cout << "   • Compare with MINOS errors from fit result" << endl;

    cout << "================================================" << endl;
}
