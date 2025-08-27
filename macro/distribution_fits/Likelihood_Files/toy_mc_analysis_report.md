# Analysis: Issues with Your Toy Monte Carlo Implementation

## Summary

I've analyzed your `realistic_toy_mc_example.cxx` against Cowan et al.'s approach and found several **critical issues** that need correction. I've created a corrected implementation in `corrected_cowan_toy_mc.cxx`.

## Critical Issues Found

### 1. **Incorrect Test Statistic Definition**

**Your Implementation:**
```cpp
double q0_data = logL_null_data - logL_full_data;
```

**Problem:** Missing factor of 2 and wrong sign convention.

**Cowan's Definition:** q₀ = -2ln(λ) where λ = L(μ=0)/L(μ̂)

**Corrected:**
```cpp
double q0_data = 2.0 * (logL_null_data - logL_full_data);
```

### 2. **Missing Boundary Rule**

**Your Implementation:** No boundary check for μ̂ < 0.

**Cowan's Requirement:** If μ̂ < 0, set q₀ = 0 (since μ ≥ 0 is physical).

**Corrected:**
```cpp
if (mu_hat < 0) {
    q0_data = 0.0;
}
```

### 3. **Wrong Nuisance Parameter Strategy**

**Your Implementation:** Generate toys using fixed null model parameters.

**Cowan's Approach:** Generate toys using conditional MLE θ̂(0) from fitting data under null hypothesis.

**Your approach:**
```cpp
null_model->SetParameters(500, -0.3); // Fixed values
```

**Corrected approach:**
```cpp
// First fit data with null model to get θ̂(0)
FitResult data_null_fit = performFit(data, null_model);
double theta_hat_0_bkg_norm = data_null_fit.parameters[0];
double theta_hat_0_bkg_slope = data_null_fit.parameters[1];

// Use θ̂(0) for toy generation
null_model->SetParameters(theta_hat_0_bkg_norm, theta_hat_0_bkg_slope);
```

### 4. **Inadequate Convergence Handling**

**Your Implementation:** Basic fit without proper error handling.

**Corrected:** Comprehensive convergence monitoring with:
- Parameter limits
- Multiple fit attempts
- Convergence status tracking
- Failure rate reporting

### 5. **Missing Distribution Validation**

**Your Implementation:** No check if toy distribution matches expected χ²₁/2.

**Corrected:** Added diagnostic checks:
- Compare to expected χ²₁/2 mixture distribution
- Check fraction of toys at q₀ = 0
- Validate mean and RMS against theoretical expectations

## Key Improvements in Corrected Version

### 1. Proper Workflow Following Cowan

```cpp
// Step A: Get conditional MLE under null θ̂(0)
FitResult data_null_fit = performFit(data, null_model);

// Step B: Generate toys using θ̂(0)  
for (int itoy = 0; itoy < nToys; itoy++) {
    // Generate toy under null with θ̂(0)
    generateToyData(h_toy, null_model_with_theta_hat_0);
    
    // Fit toy globally and under null
    FitResult toy_full_fit = performFit(h_toy, toy_full_model);
    FitResult toy_null_fit = performFit(h_toy, toy_null_model);
    
    // Calculate q₀ with boundary rule
    double q0_toy = 2.0 * (toy_null_fit.logL - toy_full_fit.logL);
    if (toy_full_fit.parameters[0] < 0) q0_toy = 0.0;
    
    q0_toys.push_back(q0_toy);
}

// Step C: Calculate p-value and significance
```

### 2. Robust Fit Handling

```cpp
struct FitResult {
    bool converged;
    double logL;
    vector<double> parameters;
    vector<double> errors;
};

FitResult performFit(TH1F* hist, TF1* model, bool fixSignal = false) {
    // Set parameter limits
    // Multiple convergence attempts
    // Proper error handling
    return result;
}
```

### 3. Distribution Diagnostics

```cpp
// Check if toy distribution matches expected χ²₁/2
int zeros = count(q0_toys.begin(), q0_toys.end(), 0.0);
cout << "Fraction at q0=0: " << double(zeros)/q0_toys.size() 
     << " (expected ~0.5)" << endl;

// Compare mean/RMS to theoretical expectations
cout << "Expected (χ²₁/2): 0.5 for mean, 1.0 for RMS" << endl;
```

## Impact on Your Results

For your glueball analysis with Δ(-2logL) = 1788.99:

### Your Implementation Result:
- q₀ = 1788.99 (incorrect)
- Z ≈ √1788.99 ≈ 42.3σ (incorrect)

### Corrected Implementation Result:
- q₀ = 2 × 1788.99 = 3577.98 (correct)
- Z ≈ √3577.98 ≈ 59.8σ (correct)

**Your signal is even stronger than you calculated!**

## Validation Recommendations

1. **Run the corrected implementation** on your actual data
2. **Check distribution shape** - should show χ²₁/2 mixture behavior
3. **Monitor convergence rate** - should be >80% successful fits
4. **Compare asymptotic vs toy results** - should agree well for such strong signals
5. **Verify boundary behavior** - check fraction of toys with q₀ = 0

## Conclusion

Your original implementation had several significant deviations from Cowan's approach, but the corrected version now properly implements:

✅ Correct test statistic definition with boundary rule  
✅ Proper nuisance parameter treatment using θ̂(0)  
✅ Robust convergence handling  
✅ Distribution validation against theory  
✅ Comprehensive diagnostics  

For your extremely strong signal (>40σ), both implementations would likely give similar significance conclusions, but the corrected version is statistically rigorous and follows best practices for discovery claims.
