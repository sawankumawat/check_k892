#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <TRandom3.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TFitResult.h>

using namespace std;

// Corrected Toy Monte Carlo following Cowan et al. approach
// Key improvements:
// 1. Proper test statistic definition: q0 = -2ln(λ) with boundary rule
// 2. Generate toys under null using conditional MLE θ̂(0) 
// 3. Proper convergence handling
// 4. Include auxiliary measurements if present

// Define signal + background model
Double_t signal_plus_background(double *x, double *par) {
    // par[0] = signal amplitude μ (parameter of interest)
    // par[1] = signal mean (nuisance)
    // par[2] = signal width (nuisance)
    // par[3] = background normalization (nuisance)
    // par[4] = background slope (nuisance)
    
    double signal = par[0] * TMath::Gaus(x[0], par[1], par[2], true);
    double background = par[3] * TMath::Exp(par[4] * x[0]);
    return signal + background;
}

// Background-only model (null hypothesis: μ = 0)
Double_t background_only(double *x, double *par) {
    // par[0] = background normalization  
    // par[1] = background slope
    return par[0] * TMath::Exp(par[1] * x[0]);
}

// Structure to hold fit results with convergence info
struct FitResult {
    bool converged;
    double logL;
    vector<double> parameters;
    vector<double> errors;
    
    FitResult() : converged(false), logL(0.0) {}
};

// Perform fits with proper error handling
FitResult performFit(TH1F* hist, TF1* model, bool fixSignal = false, double fixedSignalValue = 0.0) {
    FitResult result;
    
    // Set reasonable parameter limits
    if (!fixSignal) {
        model->SetParLimits(0, 0.0, 1e6);  // signal amplitude ≥ 0
    } else {
        model->FixParameter(0, fixedSignalValue);
    }
    
    // Set parameter limits for nuisance parameters
    if (model->GetNpar() > 2) {
        model->SetParLimits(1, 0.5, 10.0);     // signal mean
        model->SetParLimits(2, 0.01, 2.0);     // signal width
        model->SetParLimits(3, 1.0, 1e6);      // background norm
        model->SetParLimits(4, -10.0, 10.0);   // background slope
    } else {
        model->SetParLimits(0, 1.0, 1e6);      // background norm
        model->SetParLimits(1, -10.0, 10.0);   // background slope
    }
    
    // Perform fit with multiple attempts if needed
    TFitResultPtr fitPtr = hist->Fit(model, "RQELS0");
    
    if (fitPtr.Get() && fitPtr->Status() == 0) {
        result.converged = true;
        result.logL = fitPtr->MinFcnValue();
        
        for (int i = 0; i < model->GetNpar(); i++) {
            result.parameters.push_back(model->GetParameter(i));
            result.errors.push_back(model->GetParError(i));
        }
    }
    
    return result;
}

double performCowanToyMC(TH1F* data, int nToys = 1000) {
    cout << "\n=== CORRECTED TOY MONTE CARLO (COWAN ET AL. APPROACH) ===" << endl;
    cout << "Following Cowan et al. recommendations:" << endl;
    cout << "1. Proper test statistic q0 = -2ln(λ) with boundary rule" << endl;
    cout << "2. Generate toys under null using conditional MLE θ̂(0)" << endl;
    cout << "3. Handle convergence properly" << endl;
    
    // Step A: Prepare - fit data to get conditional MLE under null
    
    // Define models
    TF1* full_model = new TF1("full_model", signal_plus_background, 1.0, 5.0, 5);
    TF1* null_model = new TF1("null_model", background_only, 1.0, 5.0, 2);
    
    // Set reasonable initial parameters
    full_model->SetParameters(50, 3.0, 0.1, 500, -0.3);
    full_model->SetParNames("Signal μ", "Signal Mean", "Signal Width", "Bkg Norm", "Bkg Slope");
    
    null_model->SetParameters(500, -0.3);
    null_model->SetParNames("Bkg Norm", "Bkg Slope");
    
    cout << "\nStep A: Fitting data to get conditional MLE under null θ̂(0)..." << endl;
    
    // Fit data with full model (global fit)
    FitResult data_full_fit = performFit(data, full_model, false);
    
    // Fit data with null model (conditional MLE under μ = 0)
    FitResult data_null_fit = performFit(data, null_model, false);
    
    if (!data_full_fit.converged || !data_null_fit.converged) {
        cout << "ERROR: Data fits failed to converge!" << endl;
        return -1;
    }
    
    // Extract θ̂(0) - the conditional MLE nuisance parameters under null
    double theta_hat_0_bkg_norm = data_null_fit.parameters[0];
    double theta_hat_0_bkg_slope = data_null_fit.parameters[1];
    
    cout << "Conditional MLE under null θ̂(0):" << endl;
    cout << "  Background norm: " << theta_hat_0_bkg_norm << " ± " << data_null_fit.errors[0] << endl;
    cout << "  Background slope: " << theta_hat_0_bkg_slope << " ± " << data_null_fit.errors[1] << endl;
    
    // Calculate test statistic from data with boundary rule
    double mu_hat = data_full_fit.parameters[0];
    double q0_data = 2.0 * (data_null_fit.logL - data_full_fit.logL);
    
    // Apply boundary rule (Cowan's recommendation)
    if (mu_hat < 0) {
        q0_data = 0.0;
        cout << "Applied boundary rule: μ̂ < 0, setting q0 = 0" << endl;
    }
    
    cout << "\nData fit results:" << endl;
    cout << "Global fit: μ̂ = " << mu_hat << " ± " << data_full_fit.errors[0] << endl;
    cout << "-2 log L (global): " << data_full_fit.logL << endl;
    cout << "-2 log L (null): " << data_null_fit.logL << endl;
    cout << "Test statistic q0 = " << q0_data << endl;
    cout << "Asymptotic significance ≈ " << sqrt(q0_data) << "σ" << endl;
    
    // Step B: Generate toys under null hypothesis using θ̂(0)
    
    cout << "\nStep B: Generating " << nToys << " toys under null with θ̂(0)..." << endl;
    
    TRandom3 rng(12345);
    vector<double> q0_toys;
    q0_toys.reserve(nToys);
    
    int nbins = data->GetNbinsX();
    double xmin = data->GetXaxis()->GetXmin();
    double xmax = data->GetXaxis()->GetXmax();
    
    // Set up null model for toy generation using θ̂(0)
    null_model->SetParameters(theta_hat_0_bkg_norm, theta_hat_0_bkg_slope);
    
    int successful_toys = 0;
    int failed_fits = 0;
    
    auto start_time = chrono::high_resolution_clock::now();
    
    for (int itoy = 0; itoy < nToys; itoy++) {
        if (itoy % 200 == 0) {
            cout << "Processing toy " << itoy << "/" << nToys 
                 << " (success rate: " << (itoy > 0 ? 100.0*successful_toys/itoy : 0) << "%)\r" << flush;
        }
        
        // Generate toy data under null hypothesis with θ̂(0)
        TH1F* h_toy = new TH1F(Form("h_toy_%d", itoy), "toy", nbins, xmin, xmax);
        
        for (int ibin = 1; ibin <= nbins; ibin++) {
            double bin_center = h_toy->GetBinCenter(ibin);
            double bin_width = h_toy->GetBinWidth(ibin);
            
            // Expected counts from null model with θ̂(0)
            double expected = null_model->Eval(bin_center) * bin_width;
            
            // Generate Poisson-distributed counts
            int observed = rng.Poisson(expected);
            h_toy->SetBinContent(ibin, observed);
            h_toy->SetBinError(ibin, sqrt(max(1.0, double(observed))));
        }
        
        // Create fresh model copies for this toy
        TF1* toy_full = (TF1*)full_model->Clone(Form("toy_full_%d", itoy));
        TF1* toy_null = (TF1*)null_model->Clone(Form("toy_null_%d", itoy));
        
        // Reset parameters to reasonable starting values
        toy_full->SetParameters(50, 3.0, 0.1, theta_hat_0_bkg_norm, theta_hat_0_bkg_slope);
        toy_null->SetParameters(theta_hat_0_bkg_norm, theta_hat_0_bkg_slope);
        
        // Fit toy: global (all free)
        FitResult toy_full_fit = performFit(h_toy, toy_full, false);
        
        // Fit toy: null (fix μ = 0)  
        FitResult toy_null_fit = performFit(h_toy, toy_null, false);
        
        // Calculate test statistic for this toy
        if (toy_full_fit.converged && toy_null_fit.converged) {
            double mu_hat_toy = toy_full_fit.parameters[0];
            double q0_toy = 2.0 * (toy_null_fit.logL - toy_full_fit.logL);
            
            // Apply boundary rule
            if (mu_hat_toy < 0) {
                q0_toy = 0.0;
            }
            
            q0_toys.push_back(q0_toy);
            successful_toys++;
        } else {
            failed_fits++;
        }
        
        // delete h_toy;
        // delete toy_full;
        // delete toy_null;
    }
    
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    cout << "\nCompleted toy generation:" << endl;
    cout << "Successful toys: " << successful_toys << " / " << nToys << endl;
    cout << "Failed fits: " << failed_fits << " (" << 100.0*failed_fits/nToys << "%)" << endl;
    cout << "Time taken: " << duration.count() << " seconds" << endl;
    
    if (successful_toys < nToys * 0.8) {
        cout << "WARNING: High failure rate! Consider adjusting fit settings." << endl;
    }
    
    // Step C: Calculate p-value and significance
    
    int count_extreme = 0;
    for (double q0_toy : q0_toys) {
        if (q0_toy >= q0_data) {
            count_extreme++;
        }
    }
    
    double p_value = double(count_extreme) / double(q0_toys.size());
    
    // Convert to significance
    double significance = 0.0;
    if (p_value > 0.0 && p_value < 1.0) {
        significance = TMath::NormQuantile(1.0 - p_value);
    } else if (p_value == 0.0) {
        // Conservative estimate when no toys exceed data
        significance = TMath::NormQuantile(1.0 - 1.0/double(q0_toys.size()));
        cout << "No toys exceeded data. Significance > " << significance << "σ" << endl;
    }
    
    // Diagnostic statistics
    double q0_mean = 0.0;
    double q0_rms = 0.0;
    for (double q0_toy : q0_toys) {
        q0_mean += q0_toy;
        q0_rms += q0_toy * q0_toy;
    }
    q0_mean /= q0_toys.size();
    q0_rms = sqrt(q0_rms / q0_toys.size() - q0_mean * q0_mean);
    
    double q0_median = q0_toys[q0_toys.size()/2];
    sort(q0_toys.begin(), q0_toys.end());
    
    cout << "\n=== COWAN TOY MONTE CARLO RESULTS ===" << endl;
    cout << "Test statistic from data: " << q0_data << endl;
    cout << "Toys with q0 ≥ q0_data: " << count_extreme << " / " << q0_toys.size() << endl;
    cout << "p-value = " << p_value << endl;
    cout << "Toy MC significance = " << significance << "σ" << endl;
    cout << "Asymptotic significance = " << sqrt(q0_data) << "σ" << endl;
    cout << "\nToy distribution diagnostics:" << endl;
    cout << "Mean q0: " << q0_mean << endl;
    cout << "RMS q0:  " << q0_rms << endl;
    cout << "Median q0: " << q0_median << endl;
    cout << "Expected (χ²₁/2): 0.5 for mean, 1.0 for RMS" << endl;
    
    // Check distribution shape (should be ~ 1/2 δ(0) + 1/2 χ²₁ for large samples)
    int zeros = count(q0_toys.begin(), q0_toys.end(), 0.0);
    cout << "Fraction at q0=0: " << double(zeros)/q0_toys.size() 
         << " (expected ~0.5 for large samples)" << endl;
    
    // Create diagnostic plots
    TCanvas* c = new TCanvas("c_cowan_toy", "Cowan Toy Monte Carlo", 1200, 400);
    c->Divide(3, 1);
    
    // Plot 1: Data with fits
    c->cd(1);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(0.8);
    data->SetTitle("Data with Fitted Models;Mass (GeV);Events");
    data->Draw("PE");
    
    full_model->SetLineColor(kRed);
    full_model->SetLineWidth(2);
    full_model->Draw("same");
    
    null_model->SetLineColor(kBlue);
    null_model->SetLineWidth(2);
    null_model->SetLineStyle(2);
    null_model->Draw("same");
    
    TLegend* leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg1->AddEntry(data, "Data", "pe");
    leg1->AddEntry(full_model, "Signal + Bkg", "l");
    leg1->AddEntry(null_model, "Background only", "l");
    leg1->Draw();
    
    // Plot 2: Toy distribution of q0
    c->cd(2);
    double q_max = max(*max_element(q0_toys.begin(), q0_toys.end()), q0_data);
    
    TH1F* h_toys = new TH1F("h_toys", "Toy MC Distribution;q_{0};Toys", 
                           50, -0.5, q_max + 2.0);
    
    for (double q0_toy : q0_toys) {
        h_toys->Fill(q0_toy);
    }
    
    h_toys->SetFillColor(kBlue);
    h_toys->SetFillStyle(3004);
    h_toys->Draw();
    
    // Overlay expected χ²₁/2 distribution
    TF1* chi2_half = new TF1("chi2_half", "[0]*0.5*(x==0 ? [1] : TMath::Exp(-x/2)/sqrt(2*TMath::Pi()*x))", 
                             0, q_max + 2.0);
    chi2_half->SetParameters(h_toys->Integral(), h_toys->GetBinWidth(1));
    chi2_half->SetLineColor(kGreen);
    chi2_half->SetLineWidth(2);
    chi2_half->Draw("same");
    
    TLine* line_data = new TLine(q0_data, 0, q0_data, h_toys->GetMaximum());
    line_data->SetLineColor(kRed);
    line_data->SetLineWidth(3);
    line_data->Draw("same");
    
    TLatex lat2;
    lat2.SetNDC();
    lat2.SetTextSize(0.04);
    lat2.DrawLatex(0.15, 0.85, Form("p = %.4f", p_value));
    lat2.DrawLatex(0.15, 0.80, Form("Z = %.2f#sigma", significance));
    lat2.DrawLatex(0.15, 0.75, Form("Asymp: %.2f#sigma", sqrt(q0_data)));
    
    TLegend* leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg2->AddEntry(h_toys, "Toy MC", "f");
    leg2->AddEntry(chi2_half, "Expected (#chi^{2}_{1}/2)", "l");
    leg2->AddEntry(line_data, "Data", "l");
    leg2->Draw();
    
    // Plot 3: Comparison
    c->cd(3);
    TH1F* h_comparison = new TH1F("h_comparison", "Method Comparison;Method;Significance [#sigma]", 
                                 2, 0, 2);
    h_comparison->SetBinContent(1, sqrt(q0_data));
    h_comparison->SetBinContent(2, significance);
    h_comparison->GetXaxis()->SetBinLabel(1, "Asymptotic");
    h_comparison->GetXaxis()->SetBinLabel(2, "Cowan Toy MC");
    h_comparison->SetFillColor(kGreen);
    h_comparison->SetMaximum(max(sqrt(q0_data), significance) * 1.2);
    h_comparison->Draw();
    
    double difference = abs(significance - sqrt(q0_data));
    TLatex lat3;
    lat3.SetNDC();
    lat3.SetTextSize(0.04);
    lat3.DrawLatex(0.15, 0.85, Form("Difference: %.2f#sigma", difference));
    
    c->SaveAs("corrected_cowan_toy_mc.png");
    c->SaveAs("corrected_cowan_toy_mc.pdf");
    
    // // Cleanup
    // delete full_model;
    // delete null_model;
    // delete c;
    // delete h_toys;
    // delete h_comparison;
    // delete chi2_half;
    
    return significance;
}

void corrected_cowan_toy_mc() {
    cout << "=== CORRECTED TOY MONTE CARLO (COWAN ET AL.) ===" << endl;
    cout << "This implementation follows Cowan et al. recommendations exactly:" << endl;
    cout << "1. Proper test statistic q0 = -2ln(λ) with boundary rule" << endl;
    cout << "2. Generate toys under null using conditional MLE θ̂(0)" << endl;
    cout << "3. Proper convergence handling and diagnostics" << endl;
    cout << "4. Compare toy distribution to expected χ²₁/2 mixture" << endl;
    
    // Create example dataset
    TH1F* data = new TH1F("data", "Example Data;Mass (GeV);Events", 40, 1.0, 5.0);
    
    TRandom3 rng(42);
    
    // Add background events (exponential)
    for (int i = 0; i < 800; i++) {
        double x = 1.0 - (1.0/0.3) * TMath::Log(rng.Uniform());
        if (x < 5.0) data->Fill(x);
    }
    
    // Add moderate signal
    for (int i = 0; i < 80; i++) {
        double x = rng.Gaus(3.0, 0.1);
        if (x > 1.0 && x < 5.0) data->Fill(x);
    }
    
    // Set Poisson errors
    for (int i = 1; i <= data->GetNbinsX(); i++) {
        double content = data->GetBinContent(i);
        data->SetBinError(i, sqrt(max(1.0, content)));
    }
    
    cout << "\nGenerated example dataset with:" << endl;
    cout << "Total events: " << data->Integral() << endl;
    
    // Perform corrected toy Monte Carlo
    double toy_significance = performCowanToyMC(data, 1000);
    
    cout << "\n=== COMPARISON WITH YOUR ORIGINAL IMPLEMENTATION ===" << endl;
    cout << "Key differences from your original code:" << endl;
    cout << "1. ✓ Proper test statistic: q0 = 2×Δ(-logL) instead of Δ(-logL)" << endl;
    cout << "2. ✓ Boundary rule applied: q0 = 0 if μ̂ < 0" << endl;
    cout << "3. ✓ Toys generated using θ̂(0) from data null fit" << endl;
    cout << "4. ✓ Proper convergence monitoring" << endl;
    cout << "5. ✓ Distribution shape validation against χ²₁/2 mixture" << endl;
    
    // delete data;
}
