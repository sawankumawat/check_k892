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

// Simple example of toy Monte Carlo significance testing
// This demonstrates the robust method when asymptotic approximations may fail

// Define simple signal + background model
Double_t signal_plus_background(double *x, double *par) {
    // par[0] = signal amplitude
    // par[1] = signal mean 
    // par[2] = signal width
    // par[3] = background normalization
    // par[4] = background slope
    
    double signal = par[0] * TMath::Gaus(x[0], par[1], par[2]);
    double background = par[3] * TMath::Exp(par[4] * x[0]);
    return signal + background;
}

// Background-only model (null hypothesis)
Double_t background_only(double *x, double *par) {
    // par[0] = background normalization  
    // par[1] = background slope
    return par[0] * TMath::Exp(par[1] * x[0]);
}

double performToyMCTest(TH1F* data, int nToys = 1000) {
    cout << "\n=== TOY MONTE CARLO SIGNIFICANCE TEST EXAMPLE ===" << endl;
    cout << "Demonstrating robust method for low counts / non-Gaussian likelihood" << endl;
    
    // Define models
    TF1* full_model = new TF1("full_model", signal_plus_background, 1.0, 5.0, 5);
    TF1* null_model = new TF1("null_model", background_only, 1.0, 5.0, 2);
    
    // Set initial parameters for full model
    full_model->SetParameters(100, 3.0, 0.3, 1000, -0.5); // signal + background
    null_model->SetParameters(1000, -0.5); // background only
    
    // Fit data with both models
    cout << "Fitting data with full model (signal + background)..." << endl;
    TFitResultPtr full_fit = data->Fit(full_model, "RQELS");
    
    cout << "Fitting data with null model (background only)..." << endl;
    TFitResultPtr null_fit = data->Fit(null_model, "RQEL0S");
    
    // Calculate test statistic from data
    double logL_full_data = full_fit->MinFcnValue();
    double logL_null_data = null_fit->MinFcnValue();
    double q0_data = logL_null_data - logL_full_data;
    
    cout << "\nData fit results:" << endl;
    cout << "-2 log L (full model): " << logL_full_data << endl;
    cout << "-2 log L (null model): " << logL_null_data << endl;
    cout << "Test statistic q0 = " << q0_data << endl;
    cout << "Asymptotic significance ≈ " << sqrt(q0_data) << "σ" << endl;
    
    // Generate toy datasets under null hypothesis
    cout << "\nGenerating " << nToys << " toy datasets under null hypothesis..." << endl;
    
    TRandom3 rng(12345); // Fixed seed for reproducibility
    vector<double> q0_toys;
    q0_toys.reserve(nToys);
    
    int nbins = data->GetNbinsX();
    double xmin = data->GetXaxis()->GetXmin();
    double xmax = data->GetXaxis()->GetXmax();
    
    auto start_time = chrono::high_resolution_clock::now();
    
    for (int itoy = 0; itoy < nToys; itoy++) {
        if (itoy % 100 == 0) {
            cout << "Processing toy " << itoy << "/" << nToys << "\r" << flush;
        }
        
        // Create toy histogram
        TH1F* h_toy = new TH1F(Form("h_toy_%d", itoy), "toy", nbins, xmin, xmax);
        
        // Fill toy histogram by sampling from null model
        for (int ibin = 1; ibin <= nbins; ibin++) {
            double bin_center = h_toy->GetBinCenter(ibin);
            double bin_width = h_toy->GetBinWidth(ibin);
            
            // Expected counts from null model
            double expected = null_model->Eval(bin_center) * bin_width;
            
            // Generate Poisson-distributed counts
            int observed = rng.Poisson(expected);
            h_toy->SetBinContent(ibin, observed);
            h_toy->SetBinError(ibin, sqrt(max(1.0, double(observed)))); // Avoid zero errors
        }
        
        // Fit toy dataset with both models
        TF1* toy_full = (TF1*)full_model->Clone(Form("toy_full_%d", itoy));
        TF1* toy_null = (TF1*)null_model->Clone(Form("toy_null_%d", itoy));
        
        // Reset parameters to initial values
        toy_full->SetParameters(100, 3.0, 0.3, 1000, -0.5);
        toy_null->SetParameters(1000, -0.5);
        
        TFitResultPtr toy_null_fit = h_toy->Fit(toy_null, "RQEL0S");
        TFitResultPtr toy_full_fit = h_toy->Fit(toy_full, "RQEL0S");
        
        // Calculate test statistic for this toy
        if (toy_null_fit.Get() && toy_full_fit.Get() && 
            toy_null_fit->Status() == 0 && toy_full_fit->Status() == 0) {
            double q0_toy = toy_null_fit->MinFcnValue() - toy_full_fit->MinFcnValue();
            q0_toys.push_back(q0_toy);
        }
        
        delete h_toy;
        delete toy_full;
        delete toy_null;
    }
    
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    cout << "\nCompleted " << q0_toys.size() << " successful toy experiments in " 
         << duration.count() << " seconds" << endl;
    
    // Calculate p-value
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
        significance = TMath::NormQuantile(1.0 - 1.0/double(q0_toys.size()));
        cout << "No toys exceeded data. Significance > " << significance << "σ" << endl;
    }
    
    cout << "\n=== TOY MONTE CARLO RESULTS ===" << endl;
    cout << "Test statistic from data: " << q0_data << endl;
    cout << "Toys with q0 ≥ q0_data: " << count_extreme << " / " << q0_toys.size() << endl;
    cout << "p-value = " << p_value << endl;
    cout << "Toy MC significance = " << significance << "σ" << endl;
    
    // Create diagnostic plots
    TCanvas* c = new TCanvas("c_toy_example", "Toy Monte Carlo Example", 1200, 400);
    c->Divide(3, 1);
    
    // Plot 1: Data with fits
    c->cd(1);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(0.8);
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
    
    // Plot 2: Toy distribution
    c->cd(2);
    double q_min = *min_element(q0_toys.begin(), q0_toys.end());
    double q_max = max(*max_element(q0_toys.begin(), q0_toys.end()), q0_data);
    
    TH1F* h_toys = new TH1F("h_toys", "Toy MC Distribution;q_{0};Toys", 
                           50, q_min - 0.5, q_max + 1.0);
    
    for (double q0_toy : q0_toys) {
        h_toys->Fill(q0_toy);
    }
    
    h_toys->SetFillColor(kBlue);
    h_toys->SetFillStyle(3004);
    h_toys->Draw();
    
    TLine* line_data = new TLine(q0_data, 0, q0_data, h_toys->GetMaximum());
    line_data->SetLineColor(kRed);
    line_data->SetLineWidth(3);
    line_data->Draw("same");
    
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.15, 0.85, Form("p = %.4f", p_value));
    lat.DrawLatex(0.15, 0.80, Form("Z = %.2f#sigma", significance));
    
    // Plot 3: Cumulative distribution
    c->cd(3);
    vector<double> q0_sorted = q0_toys;
    sort(q0_sorted.begin(), q0_sorted.end());
    
    TH1F* h_cumul = new TH1F("h_cumul", "Cumulative Distribution;q_{0};P(Q_{0} #geq q_{0})", 
                            100, q_min, q_max + 1.0);
    
    for (int i = 0; i < h_cumul->GetNbinsX(); i++) {
        double q_test = h_cumul->GetBinCenter(i+1);
        int count = 0;
        for (double q_toy : q0_toys) {
            if (q_toy >= q_test) count++;
        }
        double p_test = double(count) / double(q0_toys.size());
        h_cumul->SetBinContent(i+1, p_test);
    }
    
    h_cumul->SetLineColor(kBlue);
    h_cumul->SetLineWidth(2);
    h_cumul->Draw();
    
    TLine* line_data2 = new TLine(q0_data, 0, q0_data, 1);
    line_data2->SetLineColor(kRed);
    line_data2->SetLineWidth(2);
    line_data2->Draw("same");
    
    TLine* line_pval = new TLine(q_min, p_value, q0_data, p_value);
    line_pval->SetLineColor(kRed);
    line_pval->SetLineWidth(2);
    line_pval->SetLineStyle(2);
    line_pval->Draw("same");
    
    c->SaveAs("toy_mc_significance_example.png");
    c->SaveAs("toy_mc_significance_example.pdf");
    
    // Cleanup
    delete full_model;
    delete null_model;
    delete c;
    delete h_toys;
    delete h_cumul;
    
    return significance;
}

void toy_mc_significance_example() {
    cout << "=== TOY MONTE CARLO SIGNIFICANCE TESTING EXAMPLE ===" << endl;
    cout << "This example demonstrates robust significance testing when" << endl;
    cout << "asymptotic approximations may fail due to:" << endl;
    cout << "- Low event counts" << endl;
    cout << "- Non-Gaussian likelihood shapes" << endl;
    cout << "- Complex parameter spaces" << endl;
    
    // Create example dataset with low counts
    TH1F* data = new TH1F("data", "Example Data;Mass (GeV);Events", 40, 1.0, 5.0);
    
    // Generate synthetic data with signal + background
    TRandom3 rng(42);
    
    // Add background events
    for (int i = 0; i < 500; i++) {
        double x = rng.Exp(2.0) + 1.0; // Exponential background
        if (x < 5.0) data->Fill(x);
    }
    
    // Add small signal (low counts scenario)
    for (int i = 0; i < 50; i++) {
        double x = rng.Gaus(3.0, 0.3); // Gaussian signal
        if (x > 1.0 && x < 5.0) data->Fill(x);
    }
    
    // Set Poisson errors
    for (int i = 1; i <= data->GetNbinsX(); i++) {
        double content = data->GetBinContent(i);
        data->SetBinError(i, sqrt(max(1.0, content)));
    }
    
    cout << "\nGenerated example dataset with:" << endl;
    cout << "Total events: " << data->Integral() << endl;
    cout << "This simulates a low-count scenario where asymptotic" << endl;
    cout << "approximations may not be reliable." << endl;
    
    // Perform toy Monte Carlo test
    double toy_significance = performToyMCTest(data, 1000);
    
    cout << "\n=== SUMMARY ===" << endl;
    cout << "Toy Monte Carlo provides robust p-values when:" << endl;
    cout << "1. Event counts are low (few events per bin)" << endl;
    cout << "2. Likelihood is non-Gaussian or irregular" << endl;
    cout << "3. You want to be conservative about claims" << endl;
    cout << "\nThe method works by:" << endl;
    cout << "1. Generating many toy datasets under null hypothesis (μ=0)" << endl;
    cout << "2. Computing test statistic q0 for each toy" << endl;
    cout << "3. p-value = fraction of toys with q0_toy ≥ q0_data" << endl;
    cout << "4. Converting p-value to significance via normal quantile" << endl;
    
    delete data;
}
