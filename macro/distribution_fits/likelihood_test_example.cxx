#include <iostream>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TMath.h>

using namespace std;

// Example: Likelihood test to determine significance of a peak
void likelihood_test_example()
{
    cout << "\n========================================" << endl;
    cout << "   Likelihood Test Example (Δ(-2 log L))" << endl;
    cout << "========================================\n" << endl;

    // Create example data with a Gaussian peak + background
    TH1F *h_data = new TH1F("h_data", "Example Data;Mass (GeV/c^{2});Counts", 100, 0, 10);
    
    // Fill histogram with simulated data: signal + background
    TRandom3 rand(12345);
    // Background: exponential
    for (int i = 0; i < 10000; i++) {
        double x = rand.Exp(2.0) + 1.0; // exponential background
        if (x < 10.0) h_data->Fill(x);
    }
    // Signal: Gaussian peak at mass = 5 GeV
    for (int i = 0; i < 500; i++) {
        double x = rand.Gaus(5.0, 0.2); // Gaussian signal
        if (x > 0 && x < 10.0) h_data->Fill(x);
    }

    cout << "Step 1: Generated example data with signal + background" << endl;
    cout << "        - Background: exponential distribution" << endl;
    cout << "        - Signal: Gaussian peak at 5 GeV with width 0.2 GeV" << endl;
    cout << "        - Total events: " << h_data->GetEntries() << endl << endl;

    // Define fitting functions
    // Full model: Gaussian signal + exponential background
    TF1 *f_full = new TF1("f_full", "[0]*exp(-[1]*x) + [2]*exp(-0.5*((x-[3])/[4])^2)", 1, 9);
    f_full->SetParameters(1000, 0.5, 100, 5.0, 0.2); // Initial parameters
    f_full->SetParNames("Bkg_Norm", "Bkg_Slope", "Signal_Norm", "Signal_Mean", "Signal_Width");

    // Reduced model: only exponential background (no signal)
    TF1 *f_reduced = new TF1("f_reduced", "[0]*exp(-[1]*x)", 1, 9);
    f_reduced->SetParameters(1000, 0.5); // Initial parameters
    f_reduced->SetParNames("Bkg_Norm", "Bkg_Slope");

    cout << "Step 2: Defined fitting functions" << endl;
    cout << "        - Full model: signal + background (5 parameters)" << endl;
    cout << "        - Reduced model: background only (2 parameters)" << endl;
    cout << "        - Difference: 3 parameters (signal amplitude, mean, width)" << endl << endl;

    // Fit full model (signal + background)
    cout << "Step 3: Fitting full model (signal + background)..." << endl;
    TFitResultPtr fit_full = h_data->Fit("f_full", "RQLS");
    double logL_full = fit_full->MinFcnValue();
    int npar_full = f_full->GetNpar();
    cout << "        - -2 log L (full) = " << logL_full << endl;
    cout << "        - Number of parameters = " << npar_full << endl;
    cout << "        - Chi2/NDF = " << f_full->GetChisquare() / f_full->GetNDF() << endl << endl;

    // Fit reduced model (background only)
    cout << "Step 4: Fitting reduced model (background only)..." << endl;
    TFitResultPtr fit_reduced = h_data->Fit("f_reduced", "RQ+LS");
    double logL_reduced = fit_reduced->MinFcnValue();
    int npar_reduced = f_reduced->GetNpar();
    cout << "        - -2 log L (reduced) = " << logL_reduced << endl;
    cout << "        - Number of parameters = " << npar_reduced << endl;
    cout << "        - Chi2/NDF = " << f_reduced->GetChisquare() / f_reduced->GetNDF() << endl << endl;

    // Calculate likelihood test statistic
    double delta_2logL = logL_reduced - logL_full;
    int delta_npar = npar_full - npar_reduced;
    
    cout << "========================================" << endl;
    cout << "   LIKELIHOOD TEST RESULTS" << endl;
    cout << "========================================" << endl;
    cout << "Δ(-2 log L) = (-2 log L)_reduced - (-2 log L)_full" << endl;
    cout << "            = " << logL_reduced << " - " << logL_full << endl;
    cout << "            = " << delta_2logL << endl << endl;
    
    cout << "Degrees of freedom = " << delta_npar << " (difference in parameters)" << endl << endl;

    // For nested models, Δ(-2 log L) follows χ² distribution
    // Convert to significance in σ
    double significance_sigma = sqrt(delta_2logL);
    
    cout << "Statistical interpretation:" << endl;
    cout << "- Δ(-2 log L) follows χ² distribution with " << delta_npar << " DOF" << endl;
    cout << "- Approximate significance ≈ √[Δ(-2 log L)] = " << significance_sigma << " σ" << endl << endl;

    // Critical values for χ² with different DOF
    cout << "Critical values for χ² test:" << endl;
    if (delta_npar == 1) {
        cout << "- 1σ (68.3% CL): Δ(-2 log L) > 1.00" << endl;
        cout << "- 2σ (95.4% CL): Δ(-2 log L) > 3.84" << endl;
        cout << "- 3σ (99.7% CL): Δ(-2 log L) > 6.63" << endl;
        cout << "- 5σ (99.9999% CL): Δ(-2 log L) > 25.0" << endl;
    } else if (delta_npar == 3) {
        cout << "- 1σ (68.3% CL): Δ(-2 log L) > 3.53" << endl;
        cout << "- 2σ (95.4% CL): Δ(-2 log L) > 7.81" << endl;
        cout << "- 3σ (99.7% CL): Δ(-2 log L) > 11.34" << endl;
        cout << "- 5σ (99.9999% CL): Δ(-2 log L) > 28.74" << endl;
    }
    cout << endl;

    // Final conclusion
    cout << "CONCLUSION:" << endl;
    if (delta_2logL > 25.0 && delta_npar == 3) {
        cout << "★ DISCOVERY: Signal is HIGHLY SIGNIFICANT (>5σ)" << endl;
    } else if (delta_2logL > 11.34 && delta_npar == 3) {
        cout << "★ STRONG EVIDENCE: Signal is significant (>3σ)" << endl;
    } else if (delta_2logL > 7.81 && delta_npar == 3) {
        cout << "★ EVIDENCE: Signal is significant (>2σ)" << endl;
    } else if (delta_2logL > 3.53 && delta_npar == 3) {
        cout << "★ WEAK EVIDENCE: Signal significance (~1σ)" << endl;
    } else {
        cout << "✗ NO SIGNIFICANT SIGNAL: Background-only model preferred" << endl;
    }
    cout << "========================================" << endl << endl;

    // Create comparison plot
    TCanvas *c = new TCanvas("c", "Likelihood Test Example", 800, 600);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.8);
    h_data->Draw("PE");
    
    f_full->SetLineColor(kRed);
    f_full->SetLineWidth(2);
    f_full->Draw("same");
    
    f_reduced->SetLineColor(kBlue);
    f_reduced->SetLineWidth(2);
    f_reduced->SetLineStyle(2);
    f_reduced->Draw("same");

    // Add legend
    auto legend = new TLegend(0.6, 0.7, 0.85, 0.85);
    legend->AddEntry(h_data, "Data", "PE");
    legend->AddEntry(f_full, "Signal + Background", "L");
    legend->AddEntry(f_reduced, "Background only", "L");
    legend->Draw();

    // Add text with results
    auto text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.03);
    text->DrawLatex(0.15, 0.85, Form("#Delta(-2 log L) = %.1f", delta_2logL));
    text->DrawLatex(0.15, 0.80, Form("Significance #approx %.1f #sigma", significance_sigma));

    c->SaveAs("/home/sawan/check_k892/macro/distribution_fits/likelihood_test_example.png");
    
    cout << "Example plot saved as: likelihood_test_example.png" << endl;
    
    // Clean up
    delete h_data;
    delete f_full;
    delete f_reduced;
    delete c;
}

// Function to demonstrate likelihood test for your specific case
void likelihood_test_resonance_example()
{
    cout << "\n========================================" << endl;
    cout << "   Resonance Likelihood Test Guide" << endl;
    cout << "========================================" << endl;
    cout << "For your glueball analysis, here's how to implement Δ(-2 log L):" << endl << endl;

    cout << "1. FULL MODEL (4 resonances + background):" << endl;
    cout << "   - 4 Breit-Wigner resonances: f₂(1270), a₂(1320), f'₂(1525), f₀(1710)" << endl;
    cout << "   - Background: Modified Boltzmann or exponential" << endl;
    cout << "   - Fit and get: -2 log L_full" << endl << endl;

    cout << "2. REDUCED MODEL (3 resonances + background):" << endl;
    cout << "   - Remove one resonance (e.g., fix amplitude to 0)" << endl;
    cout << "   - Keep same background model" << endl;
    cout << "   - Fit and get: -2 log L_reduced" << endl << endl;

    cout << "3. CALCULATE TEST STATISTIC:" << endl;
    cout << "   Δ(-2 log L) = (-2 log L)_reduced - (-2 log L)_full" << endl << endl;

    cout << "4. INTERPRET RESULTS:" << endl;
    cout << "   - If removing resonance increases -2 log L significantly," << endl;
    cout << "     then the resonance is statistically significant" << endl;
    cout << "   - For 1 parameter difference (amplitude only):" << endl;
    cout << "     • Δ(-2 log L) > 3.84 → 2σ evidence" << endl;
    cout << "     • Δ(-2 log L) > 6.63 → 3σ evidence" << endl;
    cout << "     • Δ(-2 log L) > 25.0 → 5σ discovery" << endl << endl;

    cout << "5. CODE IMPLEMENTATION:" << endl;
    cout << "   // Full model fit" << endl;
    cout << "   TFitResultPtr fit_full = hist->Fit(\"full_function\", \"RQLS\");" << endl;
    cout << "   double logL_full = fit_full->MinFcnValue();" << endl << endl;
    
    cout << "   // Reduced model (fix resonance amplitude to 0)" << endl;
    cout << "   full_function->FixParameter(resonance_amplitude_index, 0.0);" << endl;
    cout << "   TFitResultPtr fit_reduced = hist->Fit(\"full_function\", \"RQLS\");" << endl;
    cout << "   double logL_reduced = fit_reduced->MinFcnValue();" << endl << endl;
    
    cout << "   // Calculate test statistic" << endl;
    cout << "   double delta_2logL = logL_reduced - logL_full;" << endl;
    cout << "   double significance = sqrt(delta_2logL);" << endl;
    cout << "========================================" << endl;
}
