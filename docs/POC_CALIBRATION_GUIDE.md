# PoC Calibration Guide

This guide provides comprehensive instructions for calibrating the Probability of Correct Selection (PoC) threshold using the null/flat scenario methodology.

## Overview

The PoC calibration system uses null/flat scenarios to control Type I error rates while maintaining power in signal scenarios. This approach ensures robust statistical performance across different trial conditions.

## Methodology

### Core Principle

The calibration methodology is based on the principle that in a null/flat scenario (where no dose is truly optimal), the PoC detection rate should be controlled at approximately 10%. This provides conservative Type I error control while maintaining adequate power for signal detection.

### Mathematical Framework

#### Null/Flat Scenario Construction

For a trial with J dose levels, the null scenario is constructed as follows:

1. **Immune Response**: P_I = (φ_I, φ_I, ..., φ_I)
   - All doses have identical immune response probability φ_I

2. **Marginal Efficacy**: P_E = (φ_E, φ_E, ..., φ_E)
   - All doses have identical marginal efficacy probability φ_E
   - Calculated using total probability formula:
     P_E(j) = P(E|I=0,j) · (1-P_I(j)) + P(E|I=1,j) · P_I(j)

3. **Toxicity**: P_T = (τ_flat, τ_flat, ..., τ_flat)
   - All doses have identical safe toxicity rate τ_flat

#### Calibration Process

1. **Parameter Selection**: Choose φ_I, φ_E, and τ_flat values
2. **C_poc Testing**: Test multiple C_poc thresholds (typically 0.5-0.95)
3. **Simulation**: Run 1000+ simulations for each threshold
4. **Analysis**: Find threshold achieving ~10% PoC detection rate
5. **Validation**: Test calibrated threshold on signal scenarios

## Implementation

### Step-by-Step Calibration

#### Step 1: Set Up Calibration Parameters

```r
# Define calibration parameters
calibration_params <- list(
  phi_I = 0.2,          # Immune response threshold
  phi_E = 0.25,         # Marginal efficacy threshold
  tox_upper = 0.30,     # Toxicity upper bound
  tox_flat = 0.05,      # Safe toxicity rate
  n_simulations = 1000  # Number of simulations
)
```

#### Step 2: Create Null/Flat Scenario

```r
# Create null scenario
null_scenario <- create_null_flat_scenario(
  n_doses = 5,
  phi_I = calibration_params$phi_I,
  phi_E = calibration_params$phi_E,
  tox_upper = calibration_params$tox_upper,
  tox_flat = calibration_params$tox_flat
)
```

#### Step 3: Run Calibration

```r
# Calibrate C_poc
calibration_results <- calibrate_c_poc(
  null_scenario = null_scenario,
  c_poc_candidates = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
  n_simulations = calibration_params$n_simulations
)
```

#### Step 4: Review and Use Results

After calibration is complete, use the optimal C_poc value in your trial simulations by updating your trial configuration:

```r
# Update trial config with calibrated C_poc
trial_config$c_poc <- calibration_results$optimal_c_poc
trial_config$c_poc_calibrated <- TRUE
trial_config$calibration_date <- Sys.Date()
```

## Parameter Selection Guidelines

### Choosing φ_I (Immune Response Threshold)

- **Range**: 0.15 - 0.35
- **Typical Value**: 0.20
- **Considerations**: 
  - Lower values make immune response criteria more stringent
  - Higher values allow more doses to be admissible
  - Should reflect clinical expectations for immune response

### Choosing φ_E (Marginal Efficacy Threshold)

- **Range**: 0.15 - 0.30
- **Typical Value**: 0.25
- **Considerations**:
  - Lower values make efficacy criteria more stringent
  - Higher values require stronger efficacy signals
  - Should reflect minimum clinically meaningful efficacy

### Choosing τ_flat (Safe Toxicity Rate)

- **Range**: 0.03 - 0.08
- **Typical Value**: 0.05
- **Considerations**:
  - Should be well below toxicity threshold
  - Represents acceptable toxicity in null scenario
  - Lower values provide more conservative safety profile

### Choosing C_poc Candidates

- **Range**: 0.5 - 0.95
- **Typical Values**: c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95)
- **Considerations**:
  - Wider range provides better calibration curve
  - Higher values provide more stringent PoC validation
  - Should include values around expected optimal threshold

## Interpretation of Results

### Calibration Curve

The calibration curve shows the relationship between C_poc threshold and PoC detection rate:

- **Target**: ~10% PoC detection rate in null scenario
- **Interpretation**: Higher C_poc → Lower PoC detection rate
- **Selection**: Choose C_poc closest to 10% target

### Validation Metrics

#### True Optimal Selection Rate

- **Target**: ≥50% in signal scenarios
- **Interpretation**: Higher rates indicate better power
- **Acceptable Range**: 50% - 80%

#### Early Termination Rate

- **Target**: <80% across scenarios
- **Interpretation**: High rates may indicate overly strict thresholds
- **Investigation**: Check CT/CE/CI thresholds if rates are too high

### Performance Indicators

#### Good Calibration

- PoC detection rate ≈ 10% in null scenario
- True optimal selection rate ≥ 50% in signal scenarios
- Early termination rate < 80% across scenarios
- Smooth calibration curve with clear trend

#### Poor Calibration

- PoC detection rate >> 10% or << 10% in null scenario
- True optimal selection rate < 50% in signal scenarios
- Early termination rate > 80% in most scenarios
- Flat or erratic calibration curve

## Troubleshooting

### Common Issues

#### High Early Termination Rates (>80%)

**Symptoms**: Most trials terminate early without selecting optimal dose

**Causes**:
- CT/CE/CI thresholds too strict
- φ_T/φ_E/φ_I thresholds too restrictive
- Insufficient signal strength in scenarios

**Solutions**:
- Increase c_T, c_E, c_I values
- Decrease φ_T, φ_E, φ_I values
- Check scenario parameter settings

#### Low True Optimal Selection Rates (<50%)

**Symptoms**: Signal scenarios fail to select correct optimal dose

**Causes**:
- C_poc threshold too high
- Insufficient sample size
- Weak signal scenarios

**Solutions**:
- Lower C_poc threshold
- Increase n_simulations
- Strengthen signal scenarios

#### Erratic Calibration Curve

**Symptoms**: No clear relationship between C_poc and detection rate

**Causes**:
- Insufficient simulations
- Parameter instability
- Implementation errors

**Solutions**:
- Increase n_simulations to 2000+
- Check parameter consistency
- Verify implementation correctness

### Diagnostic Tools

#### Stage Verification

```r
# Verify stage counts and allocation
stage_verification <- verify_stage_counts(config, n_simulations = 100)
```

#### Allocation Analysis

```r
# Check allocation patterns
allocation_verification <- verify_allocation_distribution(config, n_simulations = 100)
```

#### Early Termination Diagnosis

```r
# Diagnose early termination issues
termination_diagnosis <- diagnose_early_termination(config, n_simulations = 100)
```

## Best Practices

### Simulation Settings

- **Minimum Simulations**: 1000 per C_poc value
- **Recommended Simulations**: 2000+ for robust results
- **Random Seeds**: Use fixed seeds for reproducibility
- **Parallel Processing**: Consider for large simulation runs

### Parameter Validation

- **Range Checks**: Ensure all probabilities are in [0,1]
- **Consistency Checks**: Verify marginal probability calculations
- **Sensitivity Analysis**: Test parameter robustness

### Result Documentation

- **Save All Results**: Store calibration outputs for future reference
- **Document Parameters**: Record all parameter choices and rationale
- **Version Control**: Track parameter changes and their effects

## Example Workflows

### Basic Calibration

```r
# Simple calibration workflow
source("src/optimization/poc_calibration_new.R")

# Run calibration
results <- run_poc_calibration(
  phi_I = 0.2,
  phi_E = 0.25,
  tox_upper = 0.30,
  tox_flat = 0.05,
  n_simulations = 1000
)

# View results
print(results$calibration_results$optimal_c_poc)
```

### Advanced Calibration

For advanced calibration scenarios, you can customize the parameters in the notebook or create your own calibration scripts based on the provided examples in `notebooks/poc_calibration_notebook.qmd`.

## Output Files

The calibration process generates several output files:

### Calibration Results

- `calibration_results.RData`: Complete calibration results
- `calibration_summary.csv`: Summary of calibrated parameters
- `poc_calibration_curve.png`: Calibration curve visualization

### Validation Results

- `validation_results.RData`: Validation results on signal scenarios
- `validation_summary.csv`: Summary of validation metrics
- `performance_summary.csv`: Comprehensive performance metrics

### Diagnostic Results

- `diagnostic_summary.csv`: Summary of diagnostic checks
- `stage_verification.RData`: Stage count verification results
- `allocation_verification.RData`: Allocation pattern verification
- `termination_diagnosis.RData`: Early termination analysis

### Visualizations

- `performance_comparison.png`: Performance across scenarios
- `dose_selection_frequencies.png`: Dose selection patterns
- `average_participants_per_dose.png`: Allocation effectiveness
- `combined_performance_curves.png`: Performance trade-offs

## Integration with Trial Design

### Configuration Updates

After calibration, update your trial configuration:

```r
# Update config with calibrated C_poc
trial_config$c_poc <- calibrated_c_poc
trial_config$c_poc_calibrated <- TRUE
trial_config$calibration_date <- Sys.Date()
```

### Documentation

Document the calibration process:

1. **Parameter Choices**: Record rationale for φ_I, φ_E, τ_flat
2. **Calibration Results**: Document optimal C_poc and achieved rates
3. **Validation Results**: Record performance on signal scenarios
4. **Sensitivity Analysis**: Document parameter robustness

### Reproducibility

Ensure reproducibility by:

1. **Saving Seeds**: Record random seeds used
2. **Version Control**: Track code and parameter versions
3. **Documentation**: Maintain detailed calibration logs
4. **Backup Results**: Store all calibration outputs

## Conclusion

The PoC calibration methodology provides a robust framework for controlling Type I error rates while maintaining power in Bayesian dose-finding trials. By using null/flat scenarios and systematic threshold testing, this approach ensures reliable statistical performance across diverse trial conditions.

For additional support or questions about the calibration process, refer to the implementation files in `src/optimization/poc_calibration_new.R` or the interactive notebook at `notebooks/poc_calibration_notebook.qmd`.

