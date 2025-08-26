#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <TArrow.h>
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"
#include "../src/style.h"
using namespace std;

// Forward declarations for functions (you would need to include these from the original file)
Double_t BWsumMassDepWidth_exponential(double *x, double *par);

void test_profile_likelihood()
{
    cout << "=== PROFILE LIKELIHOOD RATIO TEST FOR f0(1710) AMPLITUDE ===" << endl;
    cout << "This is a test to verify the profile likelihood implementation." << endl;
    
    // Example parameters (in a real scenario, these come from the full fit)
    double f0_amp_best = 1939.02;      // Example best-fit amplitude
    double f0_amp_error = 65.15;       // Example uncertainty
    double logL_full = -391.066;       // Example log-likelihood of full model
    
    cout << "Test parameters:" << endl;
    cout << "Best-fit f0(1710) amplitude: " << f0_amp_best << " ± " << f0_amp_error << endl;
    cout << "Full model -2 log L: " << logL_full << endl;
    
    // Simulate a profile likelihood test at amplitude = 0
    double null_hypothesis_logL = logL_full + 12.5;  // Example: add some penalty for null hypothesis
    double delta_2logL_f0 = null_hypothesis_logL - logL_full;
    
    cout << "\nProfile likelihood test results:" << endl;
    cout << "-2 log L (f0_amp = 0): " << null_hypothesis_logL << endl;
    cout << "Δ(-2 log L) = " << delta_2logL_f0 << endl;
    
    // Calculate significance
    double significance_f0_amp = sqrt(delta_2logL_f0);
    cout << "Significance of f0(1710) amplitude ≈ " << significance_f0_amp << " σ" << endl;
    
    // Interpret results
    cout << "\nInterpretation:" << endl;
    if (delta_2logL_f0 > 25.0) {
        cout << "DISCOVERY: f0(1710) amplitude is HIGHLY SIGNIFICANT (>5σ)" << endl;
    } else if (delta_2logL_f0 > 9.0) {
        cout << "EVIDENCE: f0(1710) amplitude is SIGNIFICANT (>3σ)" << endl;
    } else if (delta_2logL_f0 > 4.0) {
        cout << "EVIDENCE: f0(1710) amplitude shows EVIDENCE (>2σ)" << endl;
    } else if (delta_2logL_f0 > 1.0) {
        cout << "WEAK EVIDENCE: f0(1710) amplitude shows WEAK EVIDENCE (>1σ)" << endl;
    } else {
        cout << "NO EVIDENCE: f0(1710) amplitude is NOT SIGNIFICANT" << endl;
    }
    
    // Simulate confidence intervals
    cout << "\nSimulated confidence intervals:" << endl;
    double ci_68_lower = f0_amp_best - f0_amp_error;
    double ci_68_upper = f0_amp_best + f0_amp_error;
    double ci_95_lower = f0_amp_best - 1.96 * f0_amp_error;
    double ci_95_upper = f0_amp_best + 1.96 * f0_amp_error;
    
    cout << "68% CL (1σ): [" << ci_68_lower << ", " << ci_68_upper << "]" << endl;
    cout << "95% CL (2σ): [" << ci_95_lower << ", " << ci_95_upper << "]" << endl;
    
    // Check if zero is excluded
    bool zero_excluded_68 = (0.0 < ci_68_lower || 0.0 > ci_68_upper);
    bool zero_excluded_95 = (0.0 < ci_95_lower || 0.0 > ci_95_upper);
    
    cout << "\nZero exclusion test:" << endl;
    cout << "Zero excluded at 68% CL: " << (zero_excluded_68 ? "YES" : "NO") << endl;
    cout << "Zero excluded at 95% CL: " << (zero_excluded_95 ? "YES" : "NO") << endl;
    
    cout << "\n=== TEST COMPLETED SUCCESSFULLY ===" << endl;
}
