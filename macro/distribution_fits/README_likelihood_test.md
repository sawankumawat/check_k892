# Likelihood Test Implementation Guide: Calculating Δ(–2 log L)

## Overview

The likelihood test is a statistical method to determine the significance of a signal (e.g., a resonance) by comparing nested models. In glueball analysis, this helps quantify how statistically significant each resonance is.

## The Theory

### What is Δ(–2 log L)?

- **L** = Likelihood function
- **-2 log L** = Negative log-likelihood (what ROOT minimizes during fitting)
- **Δ(-2 log L)** = Difference between nested models

For nested models, the test statistic follows a χ² distribution:
```
Δ(-2 log L) = (-2 log L)_reduced - (-2 log L)_full
```

### Statistical Interpretation

For models differing by 1 parameter (e.g., removing one resonance):
- **Δ(-2 log L) > 1.0** → ~1σ significance (68% CL)
- **Δ(-2 log L) > 3.84** → ~2σ significance (95% CL)  
- **Δ(-2 log L) > 6.63** → ~3σ significance (99% CL)
- **Δ(-2 log L) > 25.0** → ~5σ discovery (99.9999% CL)

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

## Validation

To validate implementation:

1. **Run the example**: `root -l likelihood_test_example.cxx`
2. **Check known cases**: Test with toy MC where we know the answer
3. **Compare methods**: Cross-check with other significance estimators
4. **Systematic studies**: Vary fitting ranges, background models

## Physical Interpretation

For glueball analysis:
- **High Δ(-2 log L)**: Strong evidence for resonance existence
- **Low Δ(-2 log L)**: Resonance not statistically required
- **Intermediate values**: Marginal evidence - need more data or better analysis

This quantitative approach provides robust statistical evidence for your glueball candidates beyond traditional significance estimators.
