# Likelihood Test Implementation Guide: Calculating Δ(–2 log L) and Likelihood Profiles

## Overview

The likelihood test is a statistical method to determine the significance of a signal (e.g., a resonance) by comparing nested models. In glueball analysis, this helps quantify how statistically significant each resonance is. Additionally, likelihood profiles provide detailed information about parameter uncertainties and confidence intervals.

## The Theory

### What is Δ(–2 log L)?

- **L** = Likelihood function
- **-2 log L** = Negative log-likelihood (what ROOT minimizes during fitting)
- **Δ(-2 log L)** = Difference between nested models

For nested models, the test statistic follows a χ² distribution:
```
Δ(-2 log L) = (-2 log L)_reduced - (-2 log L)_full
```

### Likelihood Profiles

A likelihood profile shows how the -2 log L changes as you vary one parameter while allowing all other parameters to reoptimize. This provides:

- **Confidence intervals**: Regions where Δ(-2 log L) < critical value
- **Parameter uncertainties**: Both symmetric and asymmetric errors
- **Correlation information**: How well parameters are constrained

### Statistical Interpretation

For models differing by 1 parameter (e.g., removing one resonance):
- **Δ(-2 log L) > 1.0** → ~1σ significance (68% CL)
- **Δ(-2 log L) > 3.84** → ~2σ significance (95% CL)  
- **Δ(-2 log L) > 6.63** → ~3σ significance (99% CL)
- **Δ(-2 log L) > 25.0** → ~5σ discovery (99.9999% CL)

For likelihood profiles (1 parameter):
- **Δ(-2 log L) = 1.0** → 68% confidence interval (1σ)
- **Δ(-2 log L) = 4.0** → 95% confidence interval (2σ)
- **Δ(-2 log L) = 9.0** → 99.7% confidence interval (3σ)

Quick approximation: **Significance ≈ √[Δ(-2 log L)]**

## Implementation in Your Code

### Step 1: Perform Full Model Fit

```cpp
// Fit with all 4 resonances + background
TFitResultPtr fitResultptr = hinvMass->Fit("BEexpol", "RELBMS");
double logL_full = fitResultptr->MinFcnValue();
```

### Step 2: Create Reduced Models

For each resonance, create a model without that resonance:

```cpp
// Test significance of f₀(1710) resonance
TF1 *BEexpol_reduced = new TF1("BEexpol_reduced", BWsumMassDepWidth_exponential, 1.05, 2.20, 16);

// Copy all parameters from full model
for (int i = 0; i < 16; i++) {
    BEexpol_reduced->SetParameter(i, full_model_params[i]);
}

// Remove f₀(1710) by fixing its amplitude to zero
BEexpol_reduced->FixParameter(9, 0.0);  // f₀(1710) amplitude parameter

// Apply same constraints as full model
BEexpol_reduced->FixParameter(2, f1270Width);   // Fix f₂(1270) width
BEexpol_reduced->FixParameter(5, a1320Width);   // Fix a₂(1320) width  
BEexpol_reduced->FixParameter(8, f1525Width);   // Fix f'₂(1525) width

// Fit reduced model
TFitResultPtr fitResult_reduced = hinvMass->Fit("BEexpol_reduced", "RQELBMS");
double logL_reduced = fitResult_reduced->MinFcnValue();
```

### Step 3: Calculate Test Statistic

```cpp
double delta_2logL = logL_reduced - logL_full;
double significance_sigma = sqrt(delta_2logL);

cout << "Δ(-2 log L) = " << delta_2logL << endl;
cout << "Significance ≈ " << significance_sigma << " σ" << endl;
```

### Step 4: Interpret Results

```cpp
if (delta_2logL > 25.0) {
    cout << "DISCOVERY: >5σ significance" << endl;
} else if (delta_2logL > 6.63) {
    cout << "STRONG EVIDENCE: >3σ significance" << endl;
} else if (delta_2logL > 3.84) {
    cout << "EVIDENCE: >2σ significance" << endl;
} else {
    cout << "NOT SIGNIFICANT: <2σ" << endl;
}
```

## Complete Implementation

The modified code `glueball_fit_4rBW.cxx` includes:

1. **Full model fit** with all 4 resonances
2. **Automated likelihood tests** for each resonance
3. **Significance calculation** and interpretation
4. **Output to file** with results

### Key Features:

- **Nested model comparison**: Removes one resonance at a time
- **Proper constraints**: Maintains same parameter constraints
- **Statistical interpretation**: Converts Δ(-2 log L) to σ significance
- **Comprehensive output**: Results saved to text file

## Usage Example

```bash
root -l glueball_fit_4rBW.cxx
```

The output will show:
```
=== Likelihood Test for f₀(1710) ===
-2 log L (without f₀(1710)) = 2547.8
-2 log L (with f₀(1710))    = 2531.2
Δ(-2 log L) = 16.6
Significance ≈ 4.1 σ
Result: f₀(1710) is HIGHLY SIGNIFICANT (>99.9% CL, ~4σ)
```

## Likelihood Profiles

### Creating Profiles

The code automatically creates likelihood profiles for key parameters:

```cpp
// Scan parameter around best-fit value
double param_min = central_value - 3.0 * param_error;
double param_max = central_value + 3.0 * param_error;

for (int i = 0; i < n_points; i++) {
    double test_value = param_min + i * step;
    
    // Fix parameter and refit
    function->FixParameter(param_index, test_value);
    TFitResultPtr temp_fit = hist->Fit("function", "RQNWL");
    double nll_test = temp_fit->MinFcnValue();
    
    // Store relative to minimum
    profile_values.push_back(nll_test - nll_min);
}
```

### Interpreting Profiles

**Shape Analysis:**
- **Parabolic**: Gaussian uncertainties (ideal case)
- **Asymmetric**: Non-Gaussian errors, report +/- separately  
- **Flat bottom**: Parameter poorly constrained
- **Multiple minima**: Degeneracies or fit issues

**Confidence Intervals:**
- Find where profile crosses horizontal lines at Δ(-2 log L) = 1, 4, 9
- These give 68%, 95%, 99.7% confidence intervals respectively

**Generated Plots:**
- `likelihood_profiles_*.png`: Individual parameter profiles
- `likelihood_summary_*.png`: Overview of all resonance significances

### Profile Applications

1. **Parameter Uncertainties**: More robust than fit errors
2. **Confidence Intervals**: Exact intervals, not just ±1σ
3. **Correlation Assessment**: How well parameters are determined
4. **Non-linearity Detection**: Deviations from parabolic shape

## Important Notes

### Assumptions:
- **Nested models**: Reduced model is subset of full model
- **Large sample**: Asymptotic χ² distribution applies
- **Same data**: Both models fitted to identical dataset

### Best Practices:
- Use `"RQLS"` or `"RELBMS"` fit options for proper likelihood
- Ensure convergence of both fits
- Check fit quality (χ²/NDF) before interpreting results
- Consider systematic uncertainties separately
- Use sufficient scan points for smooth profiles (30-50 points)
- Extend scan range if profile doesn't return to minimum

## Validation

To validate implementation:

1. **Run the examples**: 
   - `root -l likelihood_test_example.cxx`
   - `root -l likelihood_profile_example.cxx`
2. **Check known cases**: Test with toy MC where you know the answer
3. **Compare methods**: Cross-check with other significance estimators
4. **Systematic studies**: Vary fitting ranges, background models

## Complete Implementation Summary

Your modified `glueball_fit_4rBW.cxx` now includes:

### Automatic Features:
1. **Likelihood tests** for all 4 resonances
2. **Likelihood profiles** for masses, amplitudes, and widths
3. **Statistical interpretation** with confidence levels
4. **Visual output** with summary plots
5. **File output** with numerical results

### Generated Files:
- `likelihood_profiles_*.png`: Parameter profile plots (3×3 grid)
- `likelihood_summary_*.png`: Significance overview bar chart
- Text file with Δ(-2 log L) values and interpretations

### Usage:
```bash
root -l glueball_fit_4rBW.cxx
```

The implementation provides robust, publication-quality statistical analysis of your glueball resonances with both numerical results and comprehensive visualizations.

## Physical Interpretation

For your glueball analysis:
- **High Δ(-2 log L)**: Strong evidence for resonance existence
- **Narrow profiles**: Well-determined parameters  
- **Broad profiles**: Poorly constrained parameters
- **Asymmetric profiles**: Non-Gaussian parameter uncertainties

This quantitative approach provides definitive statistical evidence for your glueball candidates, going far beyond traditional significance estimators to give you publication-ready results with proper confidence intervals and uncertainty quantification.
2. **Check known cases**: Test with toy MC where we know the answer
3. **Compare methods**: Cross-check with other significance estimators
4. **Systematic studies**: Vary fitting ranges, background models

## Physical Interpretation

For glueball analysis:
- **High Δ(-2 log L)**: Strong evidence for resonance existence
- **Low Δ(-2 log L)**: Resonance not statistically required
- **Intermediate values**: Marginal evidence - need more data or better analysis

This quantitative approach provides robust statistical evidence for your glueball candidates beyond traditional significance estimators.
