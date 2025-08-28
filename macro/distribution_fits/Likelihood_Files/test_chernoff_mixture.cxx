#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>

using namespace std;

// Function to generate Chernoff mixture samples
vector<double> generateChernoffMixture(int nSamples, int seed = 12345) {
    TRandom3 rng(seed);
    vector<double> samples;
    samples.reserve(nSamples);
    
    for (int i = 0; i < nSamples; i++) {
        double u = rng.Uniform(0.0, 1.0);
        if (u < 0.5) {
            // 50% probability: δ(0) - sample is exactly 0
            samples.push_back(0.0);
        } else {
            // 50% probability: χ²₁ distribution
            // Generate chi-squared(1) using transformation from normal
            double normal_sample = rng.Gaus(0.0, 1.0);
            samples.push_back(normal_sample * normal_sample);
        }
    }
    
    return samples;
}

void test_chernoff_mixture() {
    cout << "Testing Chernoff Mixture Implementation" << endl;
    cout << "=======================================" << endl;
    
    // Generate samples
    int nSamples = 100000;
    vector<double> chernoff_samples = generateChernoffMixture(nSamples);
    
    // Calculate empirical statistics
    double mean = 0.0;
    for (double sample : chernoff_samples) mean += sample;
    mean /= nSamples;
    
    double variance = 0.0;
    for (double sample : chernoff_samples) {
        variance += (sample - mean) * (sample - mean);
    }
    variance /= nSamples;
    double std_dev = sqrt(variance);
    
    // Count zeros
    int zero_count = 0;
    for (double sample : chernoff_samples) {
        if (sample == 0.0) zero_count++;
    }
    double zero_fraction = double(zero_count) / double(nSamples);
    
    // Theoretical values
    double theory_mean = 0.5;
    double theory_variance = 1.25;
    double theory_std = sqrt(theory_variance);
    double theory_zero_fraction = 0.5;
    
    cout << "Results for " << nSamples << " samples:" << endl;
    cout << "Mean:           " << mean << " (theory: " << theory_mean << ")" << endl;
    cout << "Std deviation:  " << std_dev << " (theory: " << theory_std << ")" << endl;
    cout << "Zero fraction:  " << zero_fraction << " (theory: " << theory_zero_fraction << ")" << endl;
    
    // Test p-value calculation
    vector<double> test_values = {0.0, 1.0, 2.71, 3.84, 6.63, 10.83};  // Common significance thresholds
    cout << "\nP-value calculations:" << endl;
    cout << "q0_test\tChernoff p-value\tPure χ²₁ p-value\tSignificance (σ)" << endl;
    
    for (double q0_test : test_values) {
        // Chernoff mixture p-value
        double chernoff_p = 0.0;
        if (q0_test > 0) {
            double chi2_tail = TMath::Prob(q0_test, 1);  // P(χ²₁ ≥ q0_test)
            chernoff_p = 0.5 * chi2_tail;
        } else {
            chernoff_p = 1.0;
        }
        
        // Pure chi-squared p-value for comparison
        double pure_chi2_p = (q0_test > 0) ? TMath::Prob(q0_test, 1) : 1.0;
        
        // Convert to significance
        double significance = (chernoff_p > 0) ? TMath::NormQuantile(1.0 - chernoff_p) : 0.0;
        
        cout << q0_test << "\t" << chernoff_p << "\t\t" << pure_chi2_p << "\t\t" << significance << endl;
    }
    
    // Create histogram and plot
    TCanvas *c = new TCanvas("c_test", "Chernoff Mixture Test", 800, 600);
    
    TH1F *h_chernoff = new TH1F("h_chernoff", "Chernoff Mixture Distribution;q_{0};Normalized Frequency", 
                                100, 0, 10);
    
    for (double sample : chernoff_samples) {
        if (sample <= 10.0) h_chernoff->Fill(sample);
    }
    
    // Normalize
    h_chernoff->Scale(1.0 / h_chernoff->Integral() / h_chernoff->GetBinWidth(1));
    h_chernoff->SetFillColor(kBlue);
    h_chernoff->SetFillStyle(3004);
    h_chernoff->Draw();
    
    // Overlay pure χ²₁ for comparison
    TF1 *chi2_1dof = new TF1("chi2_1dof", "0.5*exp(-0.5*x)/sqrt(2*TMath::Pi()*x)", 0.01, 10);
    chi2_1dof->SetLineColor(kRed);
    chi2_1dof->SetLineWidth(2);
    chi2_1dof->SetLineStyle(2);
    chi2_1dof->Draw("same");
    
    // Add legend
    TLegend *leg = new TLegend(0.6, 0.7, 0.89, 0.89);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h_chernoff, "Chernoff mixture", "f");
    leg->AddEntry(chi2_1dof, "#chi^{2}(1)", "l");
    leg->Draw();
    
    // Add statistics text
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.15, 0.85, Form("Mean: %.3f (theory: %.3f)", mean, theory_mean));
    lat.DrawLatex(0.15, 0.80, Form("Std: %.3f (theory: %.3f)", std_dev, theory_std));
    lat.DrawLatex(0.15, 0.75, Form("Zero fraction: %.3f", zero_fraction));
    
    c->SaveAs("chernoff_mixture_test.png");
    
    cout << "\nTest completed. Plot saved as 'chernoff_mixture_test.png'" << endl;
}
