# Toy Monte Carlo Significance Testing

## When Asymptotic Approximations May Fail

The standard likelihood ratio test assumes that the test statistic q₀ = -2Δlog L follows a χ² distribution (Wilks' theorem). However, this asymptotic approximation can fail when:

1. **Low event counts** - Few events per bin or low signal-to-background ratio
2. **Non-Gaussian likelihood** - Irregular or highly skewed likelihood shapes  
3. **Complex parameter spaces** - Multiple local minima or boundary effects
4. **Small samples** - Insufficient data for asymptotic regime

## Toy Monte Carlo Solution

When asymptotic approximations are questionable, toy Monte Carlo provides a robust alternative:

### Method:
1. **Generate toy datasets** under null hypothesis (μ = 0)
2. **Compute test statistic** q₀ for each toy dataset
3. **Calculate p-value** as fraction of toys with q₀ᵗᵒʸ ≥ q₀ᵈᵃᵗᵃ
4. **Convert to significance** using inverse normal CDF: Z = Φ⁻¹(1 - p)

### Implementation:

```cpp
// Example usage in your glueball fitting code:
double toy_significance = calculateToyMCSignificance(
    hinvMass,          // data histogram
    BEexpol_null,      // null model (signal amplitude = 0)  
    BEexpol_full,      // full model (signal + background)
    fitResult,         // best fit result
    1000,              // number of toy experiments
    false              // verbose output
);
```

### Key Features:

- **Robust**: Valid for any sample size or likelihood shape
- **Model-independent**: No assumptions about test statistic distribution
- **Conservative**: Typically gives more reliable p-values for edge cases
- **Slower**: Computationally intensive compared to asymptotic method

### Choosing Number of Toys:

- **100-500 toys**: Quick check, ~10% precision on p-value
- **1000 toys**: Standard analysis, ~3% precision  
- **5000+ toys**: Final results, ~1% precision
- **Rule of thumb**: Need ≥ 10/p toys for reliable p-value (e.g., 1000 toys for p ~ 0.01)

### Interpreting Results:

```
Asymptotic approximation: 2.3σ
Toy Monte Carlo:          1.8σ  
Difference:               0.5σ

→ Large difference suggests asymptotic may be unreliable
→ Use toy MC result for robust inference
```

### Example Files:

1. **`glueball_fit_4rBW.cxx`** - Main analysis with toy MC integrated
2. **`toy_mc_significance_example.cxx`** - Standalone demonstration
3. **`toy_model.cxx`** - Event generation utilities

### Performance Tips:

- Use parallel processing for large toy samples
- Cache fitted parameters to speed up toy generation
- Monitor fit convergence in toys (exclude failed fits)
- Use consistent random seeds for reproducibility

### When to Use:

**Use Toy MC when:**
- Signal significance < 3σ (borderline cases)
- Low statistics (< 100 signal events)
- Complex models with many parameters
- You need robust, conservative results

**Asymptotic OK when:**
- High statistics (> 500 events)
- Simple, well-behaved models
- Signal significance > 5σ
- Quick preliminary estimates

## References:

- Cousins, R.D. et al. "Evaluation of three methods for calculating statistical significance when incorporating a systematic uncertainty into a test of the background-only hypothesis for a Poisson process" Nucl.Instrum.Meth. A 595 (2008) 480
- Cranmer, K. et al. "HistFactory: A tool for creating statistical models for use with RooFit and RooStats" CERN-OPEN-2012-016
- Cowan, G. et al. "Asymptotic formulae for likelihood-based tests of new physics" Eur.Phys.J. C71 (2011) 1554
