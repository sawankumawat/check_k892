# Profile Likelihood Ratio Test Implementation for f0(1710) Amplitude

## Overview

This implementation adds a comprehensive profile likelihood ratio test specifically targeting the f0(1710) amplitude as the parameter of interest, treating all other parameters as nuisance parameters. This is added to the existing `glueball_fit_4rBW.cxx` code.

## Theory

### Profile Likelihood Method

For a parameter of interest θ and nuisance parameters η, the profile likelihood is:

```
L_p(θ) = L(θ, η̂(θ))
```

where η̂(θ) are the values of nuisance parameters that maximize the likelihood for fixed θ.

### Test Statistic

The profile likelihood ratio test statistic is:

```
Δ(-2 log L) = -2 log L_p(θ) - (-2 log L_max)
```

where:
- L_p(θ) is the profile likelihood at test value θ  
- L_max is the global maximum likelihood

### Statistical Distribution

For the null hypothesis H₀: θ = θ₀ vs H₁: θ ≠ θ₀, the test statistic Δ(-2 log L) follows asymptotically a χ² distribution with 1 degree of freedom.

## Implementation Details

### Parameter of Interest
- **f0(1710) amplitude** (parameter index 9 in the fit function)

### Nuisance Parameters 
All other parameters in the 4-resonance Breit-Wigner fit:
- f₂(1270): amplitude, mass
- a₂(1320): amplitude, mass  
- f'₂(1525): amplitude, mass
- f₀(1710): mass, width
- Background: exponential parameters

### Constrained Parameters
The following parameters remain fixed during profiling (not treated as nuisance):
- f₂(1270) width: fixed to PDG value
- a₂(1320) width: fixed to PDG value
- f'₂(1525) width: fixed to PDG value

## Key Features

### 1. Null Hypothesis Test
Tests the significance of f0(1710) amplitude by comparing:
- **H₀**: f0(1710) amplitude = 0 (no f0(1710) signal)
- **H₁**: f0(1710) amplitude ≠ 0 (f0(1710) signal present)

### 2. Profile Likelihood Scan
Creates a detailed scan of the profile likelihood as a function of f0(1710) amplitude:
- Range: -2σ to +4σ around best-fit value
- 60 scan points for smooth curve
- Includes negative amplitude values to test physical boundaries

### 3. Confidence Intervals
Calculates confidence intervals from profile likelihood:
- **68% CL** (1σ): Δ(-2 log L) = 1.0
- **95% CL** (2σ): Δ(-2 log L) = 4.0  
- **99.7% CL** (3σ): Δ(-2 log L) = 9.0

### 4. Zero Exclusion Test
Tests whether zero amplitude is excluded at different confidence levels.

## Output Files

### 1. Profile Likelihood Plot
- **File**: `profile_likelihood_f0_amplitude_default.png`
- Shows profile likelihood curve with confidence level lines
- Marks null hypothesis point and best-fit point

### 2. Summary Comparison Plot  
- **File**: `likelihood_summary_default.png`
- Two-panel plot comparing:
  - Left: Standard resonance presence/absence tests
  - Right: f0(1710) amplitude profile test

### 3. Text Output
Results written to fit parameters file with:
- Best-fit amplitude and uncertainty
- Δ(-2 log L) value for null hypothesis
- Significance in σ units
- Confidence intervals
- Zero exclusion status

## Interpretation

### Significance Levels
- **> 5σ**: Discovery level (Δ(-2 log L) > 25)
- **> 3σ**: Evidence level (Δ(-2 log L) > 9)  
- **> 2σ**: Weak evidence (Δ(-2 log L) > 4)
- **> 1σ**: Minimal evidence (Δ(-2 log L) > 1)

### Physical Interpretation
- **Large Δ(-2 log L)**: Strong evidence for non-zero f0(1710) amplitude
- **Zero excluded**: f0(1710) signal is statistically significant
- **Confidence intervals**: Range of physically allowed amplitude values

## Technical Implementation

### Profile Likelihood Function
```cpp
auto profileLikelihood_f0_amp = [&](double test_amplitude) -> double {
    // Fix f0(1710) amplitude to test value
    // Refit with all other parameters floating as nuisance parameters
    // Return -2 log L value
};
```

### Confidence Interval Calculation
Uses linear interpolation to find amplitude values where Δ(-2 log L) crosses threshold levels.

### Error Handling
- Checks for successful fits at each scan point
- Handles boundary cases (negative amplitudes)
- Robust against fit failures

## Usage

The profile likelihood test runs automatically when executing `glueball_fit_4rBW.cxx` with the macro configured for the `b_massdepWidth_modifiedBoltzmann` option.

## Advantages Over Standard Methods

1. **Proper treatment of nuisance parameters**: Accounts for correlations and uncertainties in other fit parameters
2. **Non-Gaussian errors**: Works even when parameter uncertainties are non-Gaussian
3. **Physical boundaries**: Can handle amplitude constraints (e.g., positivity)
4. **Model comparison**: Direct statistical test of signal presence vs. absence
5. **Confidence intervals**: Provides range of physically allowed values, not just point estimate

This implementation provides a robust, statistically rigorous test for the significance of the f0(1710) amplitude in the glueball analysis.
