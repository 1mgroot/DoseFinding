# Parameter Optimization and PoC Calibration Guide

This comprehensive guide covers both parameter optimization and Probability of Correct Selection (PoC) calibration for Bayesian dose-finding trials.

## Table of Contents

1. [Overview](#overview)
2. [Parameter Optimization](#parameter-optimization)
3. [PoC Calibration](#poc-calibration)
4. [Integration Workflow](#integration-workflow)
5. [Best Practices](#best-practices)
6. [Troubleshooting](#troubleshooting)

---

## Overview

### What is "Good Performance"?

A well-configured Bayesian dose-finding trial should:

1. **Complete trials**: Avoid excessive early termination
2. **Select correct doses**: Identify truly optimal doses
3. **Maximize utility**: Achieve high utility scores
4. **Control Type I error**: Limit false positive selections (~10%)
5. **Maintain power**: Achieve ≥50% true optimal selection in signal scenarios
6. **Allocate efficiently**: Balance participants across dose levels

### Two-Stage Optimization Process

1. **Parameter Optimization**: Find optimal threshold and credibility parameters (φ and c values)
2. **PoC Calibration**: Calibrate the PoC threshold (C_poc) using null/flat scenarios

---

## Parameter Optimization

### Key Parameters

#### Threshold Parameters (φ values)

- **φ_T (Toxicity Threshold)**: Upper limit for acceptable toxicity rate
  - Range: 0.25 - 0.45
  - Typical: 0.30
  - Lower = more conservative safety

- **φ_E (Efficacy Threshold)**: Lower limit for acceptable efficacy rate
  - Range: 0.05 - 0.25
  - Typical: 0.15
  - Higher = more stringent efficacy requirements

- **φ_I (Immune Response Threshold)**: Lower limit for acceptable immune response rate
  - Range: 0.15 - 0.35
  - Typical: 0.20
  - Higher = more stringent immune requirements

#### Credibility Parameters (c values)

- **c_T (Toxicity Credibility)**: Required probability that toxicity < φ_T
  - Range: 0.7 - 0.95
  - Typical: 0.8-0.9
  - Higher = more conservative

- **c_E (Efficacy Credibility)**: Required probability that efficacy > φ_E
  - Range: 0.6 - 0.9
  - Typical: 0.7-0.8
  - Higher = more stringent

- **c_I (Immune Credibility)**: Required probability that immune response > φ_I
  - Range: 0.6 - 0.9
  - Typical: 0.7-0.8
  - Higher = more stringent

#### Utility Functions

Define the value of different outcome combinations:

- **Conservative**: Lower rewards, more sensitive to toxicity
- **Balanced**: Moderate risk-reward trade-off (recommended starting point)
- **Aggressive**: Higher rewards, more toxicity tolerance

### Usage

#### Quick Optimization (2 hours)

Test 20 parameter combinations with 3 simulations each:

```r
# In run_optimization.R or your script
source("src/optimization/parameter_optimization.R")

# Run quick optimization
quick_results <- quick_optimization()

# Analyze results
analyze_results(quick_results)
```

#### Comprehensive Optimization (4-6 hours)

Test 50 parameter combinations with 5 simulations each:

```r
# Run comprehensive optimization
comprehensive_results <- comprehensive_optimization()

# Analyze results
analyze_results(comprehensive_results)
```

#### Test Specific Parameters

```r
# Test specific parameter combination
test_results <- test_specific_params(
  phi_T = 0.30,        # Toxicity threshold
  phi_E = 0.15,        # Efficacy threshold
  phi_I = 0.20,        # Immune response threshold
  c_T = 0.8,           # Toxicity credibility
  c_E = 0.7,           # Efficacy credibility
  c_I = 0.7,           # Immune credibility
  utility_type = "balanced"
)
```

### Evaluation Metrics

1. **Completion Rate**: Proportion of trials that complete without early termination
   - Target: > 0.8

2. **Correct Selection Rate**: Proportion selecting the true optimal dose
   - Target: > 0.7

3. **Average Final Utility**: Mean utility score across simulations
   - Target: > 60 (depends on utility table)

4. **Allocation Efficiency**: Balance of participant distribution
   - Measured by allocation variance

5. **Sample Size**: Average number of participants used
   - Lower is better (given other metrics are good)

### Parameter Tuning Strategies

#### Conservative Strategy
- Lower toxicity threshold (φ_T = 0.25-0.30)
- Higher credibility requirements (c_T, c_E, c_I = 0.8-0.95)
- Conservative utility function
- **Use when**: Safety is paramount, regulatory requirements are strict

#### Balanced Strategy (Recommended)
- Moderate thresholds (φ_T = 0.30-0.35, φ_E = 0.10-0.15)
- Moderate credibility (c_T, c_E, c_I = 0.7-0.8)
- Balanced utility function
- **Use when**: Standard dose-finding with typical risk-benefit trade-offs

#### Aggressive Strategy
- Higher toxicity threshold (φ_T = 0.35-0.45)
- Lower credibility requirements (c_T, c_E, c_I = 0.6-0.7)
- Aggressive utility function
- **Use when**: Severe disease, limited alternatives, benefit outweighs risk

### Visualization Outputs

Parameter optimization generates:

1. **Completion Rate vs Parameters**: Shows parameter impact on trial completion
2. **Selection Rate vs Parameters**: Shows parameter impact on correct dose selection
3. **Utility vs Parameters**: Shows parameter impact on final utility
4. **Sensitivity Analysis**: Shows importance of each parameter
5. **Top Combinations**: Best 10 parameter sets ranked by utility

---

## PoC Calibration

### Overview

PoC calibration controls Type I error rates by testing how often the trial incorrectly validates a dose selection in null/flat scenarios (where no dose is truly optimal).

### Core Methodology

#### Null/Flat Scenario Construction

For a trial with J dose levels:

1. **Immune Response**: P_I = (null_p_I, null_p_I, ..., null_p_I)
   - All doses have identical immune response
   - Set slightly above threshold to compensate for PAVA bias
   - Example: null_p_I = 0.25 when φ_I = 0.20

2. **Marginal Efficacy**: P_E = (null_p_E, null_p_E, ..., null_p_E)
   - All doses have identical efficacy
   - Set slightly above threshold to compensate for PAVA bias
   - Example: null_p_E = 0.30 when φ_E = 0.25

3. **Toxicity**: P_T = (null_p_T, null_p_T, ..., null_p_T)
   - All doses have identical safe toxicity
   - Set well below threshold
   - Example: null_p_T = 0.05 when φ_T = 0.30

#### Why Set Above Threshold?

**PAVA Bias**: The Pool Adjacent Violators Algorithm (PAVA) enforces monotonicity by averaging adjacent values. In flat scenarios, this creates a systematic downward bias, especially with small samples. Setting null values above thresholds compensates for this bias.

### Calibration Process

#### Step 1: Set Up Calibration Parameters

```r
# Define calibration parameters
calibration_params <- list(
  null_p_I = 0.25,      # True immune response (above threshold)
  null_p_E = 0.30,      # True efficacy (above threshold)
  null_p_T = 0.05,      # True toxicity (below threshold)
  n_simulations = 1000  # Number of simulations per C_poc value
)
```

#### Step 2: Create Null/Flat Scenario

```r
# Source calibration functions
source("src/optimization/poc_calibration_new.R")

# Create null scenario
null_scenario <- create_null_flat_scenario(
  n_doses = 5,
  phi_I = calibration_params$null_p_I,  # Note: parameter name is phi_I but value is null_p_I
  phi_E = calibration_params$null_p_E,
  tox_upper = 0.30,
  tox_flat = calibration_params$null_p_T
)
```

#### Step 3: Run Calibration

```r
# Test multiple C_poc thresholds
calibration_results <- calibrate_c_poc(
  null_scenario = null_scenario,
  c_poc_candidates = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
  n_simulations = calibration_params$n_simulations,
  base_config = trial_config
)
```

#### Step 4: Generate Detailed Report

```r
# Generate comprehensive calibration report
generate_calibration_report(
  calibration_results = calibration_results,
  null_scenario = null_scenario,
  base_config = trial_config,
  file_path = "results/calibration_report.txt"
)
```

#### Step 5: Use Calibrated C_poc

```r
# Update trial config with calibrated value
trial_config$c_poc <- calibration_results$optimal_c_poc
trial_config$c_poc_calibrated <- TRUE
trial_config$calibration_date <- Sys.Date()
```

### Interpretation of Results

#### Calibration Curve

The calibration curve shows C_poc vs PoC detection rate:

- **X-axis**: C_poc threshold (0.5 to 0.95)
- **Y-axis**: PoC detection rate (0 to 1)
- **Target**: Find C_poc where detection rate ≈ 10%
- **Trend**: Higher C_poc → Lower detection rate

#### Performance Indicators

**Good Calibration:**
- PoC detection rate ≈ 10% in null scenario
- True optimal selection ≥ 50% in signal scenarios
- Early termination rate < 80%
- Smooth calibration curve

**Poor Calibration:**
- PoC detection rate >> 10% (too liberal) or << 10% (too conservative)
- True optimal selection < 50%
- Early termination rate > 80%
- Erratic calibration curve

#### Detailed Report Contents

The calibration report includes:

1. **Configuration Summary**: Trial design parameters, thresholds, null scenario values
2. **Calibration Results**: Table of C_poc values with detection and termination rates
3. **Early Termination Analysis**: 
   - Statistics by C_poc value
   - Stage distribution of terminations
   - Sample size at termination
   - 3 detailed example cases showing why termination occurred
4. **Recommendations**: Optimal C_poc value and expected trial characteristics

### Parameter Selection Guidelines

#### Null Scenario Parameters

| Parameter | Range | Typical | Purpose |
|-----------|-------|---------|---------|
| null_p_I | 0.20-0.30 | 0.25 | True immune response in null scenario |
| null_p_E | 0.25-0.35 | 0.30 | True efficacy in null scenario |
| null_p_T | 0.03-0.08 | 0.05 | True toxicity in null scenario |

**Key Principle**: Set null_p_I and null_p_E slightly above their corresponding thresholds (φ_I, φ_E) to compensate for PAVA bias.

#### C_poc Candidates

- **Range**: 0.5 - 0.95
- **Recommended**: c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95)
- **Spacing**: Closer spacing near expected optimal value
- **Coverage**: Wide range ensures you capture the full relationship

---

## Integration Workflow

### Complete Optimization Pipeline

Follow this sequence for best results:

#### Phase 1: Parameter Optimization (Day 1)

```r
# 1. Run quick optimization to identify promising regions
source("src/optimization/parameter_optimization.R")
quick_results <- quick_optimization()

# 2. Analyze results and identify top performers
analyze_results(quick_results)

# 3. Test specific promising combinations
best_params <- test_specific_params(
  phi_T = 0.30,
  phi_E = 0.15,
  phi_I = 0.20,
  c_T = 0.8,
  c_E = 0.7,
  c_I = 0.7,
  utility_type = "balanced"
)

# 4. Update trial_config with optimized parameters
trial_config$phi_T <- 0.30
trial_config$phi_E <- 0.15
trial_config$phi_I <- 0.20
trial_config$c_T <- 0.8
trial_config$c_E <- 0.7
trial_config$c_I <- 0.7
```

#### Phase 2: PoC Calibration (Day 2)

```r
# 1. Source config with optimized parameters
source("src/core/config.R")
source("src/optimization/poc_calibration_new.R")

# 2. Set up null scenario
null_scenario <- create_null_flat_scenario(
  n_doses = 5,
  phi_I = 0.25,  # Above threshold φ_I = 0.20
  phi_E = 0.30,  # Above threshold φ_E = 0.25
  tox_upper = 0.30,
  tox_flat = 0.05
)

# 3. Run calibration
calibration_results <- calibrate_c_poc(
  null_scenario = null_scenario,
  c_poc_candidates = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
  n_simulations = 1000,
  base_config = trial_config
)

# 4. Generate and review report
generate_calibration_report(
  calibration_results = calibration_results,
  null_scenario = null_scenario,
  base_config = trial_config,
  file_path = "results/calibration_report.txt"
)

# 5. Update config with calibrated C_poc
trial_config$c_poc <- calibration_results$optimal_c_poc
trial_config$c_poc_calibrated <- TRUE
trial_config$calibration_date <- Sys.Date()
```

#### Phase 3: Validation (Day 3)

```r
# Run final validation with optimized and calibrated parameters
source("notebooks/simulation_notebook.qmd")

# Test on multiple scenarios:
# - Signal scenarios (dose 3 or 5 optimal)
# - Null scenarios (no dose optimal)
# - Edge cases (all toxic, all ineffective)
```

### Using Notebooks

#### PoC Calibration Notebook

```r
# Open in RStudio
# notebooks/poc_calibration_notebook.qmd

# The notebook provides:
# - Interactive parameter setting
# - Calibration execution
# - Automatic visualization
# - Report generation
```

#### Simulation Notebook

```r
# Open in RStudio
# notebooks/simulation_notebook.qmd

# The notebook provides:
# - Trial simulation with optimized parameters
# - Multi-scenario analysis
# - Performance visualization
# - Summary statistics
```

---

## Best Practices

### Simulation Settings

- **Minimum Simulations**: 1000 per C_poc value for calibration, 100+ for parameter testing
- **Recommended Simulations**: 2000+ for final calibration
- **Random Seeds**: Use reproducible but varied seeds (`10000 * index + sim`)
- **Verbose Logging**: Disable for large runs (`verbose_logging = FALSE`)

### Parameter Validation

- **Range Checks**: Ensure probabilities in [0,1]
- **Consistency Checks**: Verify null_p values relative to φ thresholds
- **Sensitivity Analysis**: Test parameter robustness
- **Cross-validation**: Test on held-out scenarios

### Result Documentation

Essential documentation:

1. **Parameter Choices**: Record rationale for all parameter selections
2. **Optimization Results**: Save all optimization outputs
3. **Calibration Results**: Save calibration curves and reports
4. **Validation Results**: Record performance on test scenarios
5. **Version Information**: Track code versions used

### File Management

```
results/
├── optimization/
│   ├── quick_optimization_results.RData
│   ├── comprehensive_results.RData
│   └── optimization_plots/
├── calibration/
│   ├── calibration_results.RData
│   ├── calibration_report.txt
│   ├── poc_calibration_curve.png
│   └── combined_performance_curves.png
└── validation/
    ├── signal_scenarios.RData
    ├── null_scenarios.RData
    └── validation_plots/
```

---

## Troubleshooting

### Parameter Optimization Issues

#### Trials Terminate Too Early (>50% termination)

**Symptoms**: Most trials end before final stage

**Causes**:
- Thresholds too strict (φ values)
- Credibility requirements too high (c values)
- Sample size too small

**Solutions**:
- Increase φ_I and φ_E (make less stringent)
- Decrease φ_T (allow more toxicity)
- Lower c_T, c_E, c_I values
- Increase cohort_size

#### Wrong Dose Selected Frequently (<50% correct)

**Symptoms**: Rarely selects true optimal dose

**Causes**:
- Utility function misaligned
- Insufficient sample size
- Thresholds don't match scenario characteristics

**Solutions**:
- Adjust utility table to better reward optimal outcomes
- Increase n_stages or cohort_size
- Re-examine φ and c values
- Check if scenarios are realistic

#### Low Utility Scores

**Symptoms**: Average utility < 50

**Causes**:
- Utility function too conservative
- Parameters too restrictive
- Poor dose allocation

**Solutions**:
- Use balanced or aggressive utility function
- Relax parameters slightly
- Check adaptive randomization is working

### PoC Calibration Issues

#### High Early Termination in Null Scenario (>80%)

**Symptoms**: Most null scenario trials terminate early

**Causes**:
- null_p values too close to thresholds
- PAVA bias not adequately compensated
- c_T, c_E, c_I too strict

**Solutions**:
- Increase null_p_I and null_p_E further above thresholds
  - Try null_p_I = φ_I + 0.10
  - Try null_p_E = φ_E + 0.10
- Lower c_I, c_E, c_T values
- Increase sample size (cohort_size)

#### PoC Detection Rate Too High (>20%)

**Symptoms**: False positive rate exceeds target

**Causes**:
- C_poc threshold too low
- Early termination removing difficult cases
- Insufficient simulations (high variance)

**Solutions**:
- Increase C_poc threshold
- Check early termination rates
- Increase n_simulations to 2000+

#### PoC Detection Rate Too Low (<5%)

**Symptoms**: Almost never validates selections

**Causes**:
- C_poc threshold too high
- null_p values set too low
- Sample size insufficient

**Solutions**:
- Decrease C_poc threshold
- Increase null_p values
- Increase sample size

#### Erratic Calibration Curve

**Symptoms**: No smooth trend between C_poc and detection rate

**Causes**:
- Insufficient simulations
- Parameter instability
- Implementation errors

**Solutions**:
- Increase n_simulations to 2000+
- Check parameter consistency
- Verify implementation against TRIAL_DESIGN.md
- Review verbose logs for patterns

### General Issues

#### Long Runtime

**Symptoms**: Optimization/calibration takes too long

**Solutions**:
- Start with quick_optimization (fewer combinations)
- Reduce n_simulations for initial exploration
- Use verbose_logging = FALSE
- Consider parallel processing (future enhancement)

#### Memory Issues

**Symptoms**: R runs out of memory

**Solutions**:
- Process results in batches
- Clear large objects after use: `rm(object); gc()`
- Reduce posterior sample size if possible
- Monitor memory: `pryr::mem_used()`

#### Inconsistent Results

**Symptoms**: Different runs give different rankings

**Causes**:
- Insufficient simulations
- High variance scenarios
- Seed management issues

**Solutions**:
- Increase n_simulations
- Use fixed seeds for reproducibility
- Check scenario definitions
- Average over multiple runs

---

## Summary

### Recommended Workflow

1. **Day 1 - Parameter Optimization**:
   - Run quick_optimization() (2 hours)
   - Identify top 3-5 parameter combinations
   - Test specific combinations with more simulations

2. **Day 2 - PoC Calibration**:
   - Set up null/flat scenario
   - Run calibrate_c_poc() with 1000+ simulations (3-4 hours)
   - Generate and review calibration report
   - Select optimal C_poc value

3. **Day 3 - Validation**:
   - Test on signal scenarios
   - Verify Type I error control (~10%)
   - Verify power (≥50% true optimal selection)
   - Document final configuration

### Key Success Metrics

- **Completion Rate**: > 70% in realistic scenarios
- **Correct Selection**: > 60% in signal scenarios
- **Type I Error**: ≈ 10% in null scenarios
- **Early Termination**: < 80% across scenarios
- **Utility**: Maximized subject to constraints

### Final Checklist

Before finalizing parameters:

- [ ] Optimization completed with multiple parameter sets tested
- [ ] PoC calibration completed with ≥1000 simulations
- [ ] Calibration report reviewed and understood
- [ ] Type I error rate controlled at ~10%
- [ ] Power validated at ≥50% for signal scenarios
- [ ] All parameters documented with rationale
- [ ] Results saved and version controlled
- [ ] Trial config updated with calibrated values
- [ ] Validation runs completed successfully
- [ ] Performance visualizations reviewed

---

## References

- **Implementation**: `src/optimization/parameter_optimization.R`, `src/optimization/poc_calibration_new.R`
- **Examples**: `notebooks/simulation_notebook.qmd`, `notebooks/poc_calibration_notebook.qmd`
- **Design Specification**: `docs/TRIAL_DESIGN.md`
- **Naming Conventions**: `NAMING_CONVENTION.md`
- **Updates**: `CALIBRATION_UPDATE_SUMMARY.md`

