# How to Run the Bayesian Dose-Finding Trial Simulation

This guide explains the different ways to run simulations, optimize parameters, and calibrate the trial design.

## Quick Start

For first-time users, we recommend starting with the **Interactive Simulation Notebook** to understand how the trial works, then proceeding to **Parameter Optimization** and **PoC Calibration** for production use.

---

## Running Methods

### 1. Interactive Simulation Notebook (Recommended for Exploration)

The simulation notebook provides an interactive environment to configure, run, and visualize trial simulations.

**Location**: `notebooks/simulation_notebook.qmd`

**Steps**:

1. **Open the notebook** in RStudio:
   ```r
   # In RStudio, navigate to:
   # File > Open File > notebooks/simulation_notebook.qmd
   ```

2. **Set working directory** to project root:
   ```r
   setwd("/path/to/DoseFinding")
   ```

3. **Configure trial parameters** in the Configuration section:
   ```r
   # Modify trial_config list
   trial_config <- list(
     dose_levels = c(1, 2, 3, 4, 5),
     n_stages = 3,
     cohort_size = 15,
     phi_T = 0.30,  # Toxicity threshold
     phi_E = 0.25,  # Efficacy threshold
     phi_I = 0.20,  # Immune response threshold
     c_T = 0.8,     # Toxicity credibility
     c_E = 0.7,     # Efficacy credibility
     c_I = 0.7,     # Immune response credibility
     c_poc = 0.9,   # PoC threshold (should be calibrated)
     enable_early_termination = TRUE,
     log_early_termination = TRUE
   )
   ```

4. **Define scenario probabilities**:
   ```r
   # Immune response probabilities by dose
   p_YI <- c(0.2, 0.4, 0.6, 0.7, 0.8)
   
   # Toxicity and efficacy conditional probabilities
   # Rows: doses, Columns: [I=0, I=1]
   p_YT_given_I <- matrix(c(
     0.05, 0.10,  # Dose 1
     0.10, 0.15,  # Dose 2
     0.15, 0.25,  # Dose 3
     0.25, 0.35,  # Dose 4
     0.35, 0.50   # Dose 5
   ), ncol = 2, byrow = TRUE)
   ```

5. **Run code chunks sequentially** to:
   - Source required files
   - Run single or multiple trial simulations
   - Generate visualizations
   - View summary statistics

**Outputs**:
- Trial allocation plots
- Dose-response curves
- Selection frequency tables
- Summary statistics

**Use case**: Learning, scenario testing, single trial analysis

---

### 2. PoC Calibration Notebook (Recommended for Calibration)

The PoC calibration notebook guides you through calibrating the C_poc threshold using null/flat scenarios.

**Location**: `notebooks/poc_calibration_notebook.qmd`

**Steps**:

1. **Open the notebook** in RStudio:
   ```r
   # File > Open File > notebooks/poc_calibration_notebook.qmd
   ```

2. **Set working directory** to project root:
   ```r
   setwd("/path/to/DoseFinding")
   ```

3. **Configure calibration parameters**:
   ```r
   calibration_params <- list(
     null_p_I = 0.25,      # True immune response in null scenario
     null_p_E = 0.30,      # True efficacy in null scenario
     null_p_T = 0.05,      # True toxicity in null scenario
     n_simulations = 1000  # Number of simulations per C_poc value
   )
   ```

4. **Run calibration workflow**:
   - Create null/flat scenario
   - Test multiple C_poc candidates
   - Generate calibration curves
   - Review detailed report
   - Select optimal C_poc value

5. **Update trial config** with calibrated value:
   ```r
   trial_config$c_poc <- calibration_results$optimal_c_poc
   trial_config$c_poc_calibrated <- TRUE
   trial_config$calibration_date <- Sys.Date()
   ```

**Outputs**:
- Calibration curve (C_poc vs detection rate)
- Early termination analysis
- Detailed calibration report (text file)
- Combined performance curves

**Use case**: Calibrating C_poc threshold for Type I error control (~10%)

---

### 3. Direct R Script Execution

Run simulations programmatically for automation or batch processing.

**Steps**:

1. **Source configuration**:
   ```r
   source("src/core/config.R")
   ```

2. **Source main simulation**:
   ```r
   source("src/core/main.R")
   ```

3. **Run trial simulation**:
   ```r
   result <- run_trial_simulation(
     trial_config = trial_config,
     p_YI = p_YI,
     p_YT_given_I = p_YT_given_I,
     p_YE_given_I = p_YE_given_I,
     rho0 = 1.5,
     rho1 = 2.0,
     seed = 123  # Optional: for reproducibility
   )
   ```

4. **Access results**:
   ```r
   result$selected_dose
   result$final_utility
   result$data
   result$early_termination
   result$poc_results
   ```

**Use case**: Batch simulations, automation, integration with other tools

---

### 4. Parameter Optimization

Systematically search for optimal trial parameters (φ and c values).

**Location**: `src/optimization/parameter_optimization.R`

**Quick Optimization** (2 hours):

```r
source("src/optimization/parameter_optimization.R")

# Test 20 parameter combinations, 3 simulations each
results <- quick_optimization()

# Analyze results
analyze_results(results)
```

**Comprehensive Optimization** (4-6 hours):

```r
# Test 50 parameter combinations, 5 simulations each
results <- comprehensive_optimization()

# Analyze results
analyze_results(results)
```

**Test Specific Parameters**:

```r
results <- test_specific_params(
  phi_T = 0.30,
  phi_E = 0.15,
  phi_I = 0.20,
  c_T = 0.8,
  c_E = 0.7,
  c_I = 0.7,
  utility_type = "balanced"
)
```

**Outputs**:
- Completion rate vs parameters
- Correct selection rate vs parameters
- Utility vs parameters
- Sensitivity analysis plots
- Top 10 parameter combinations

**Use case**: Finding optimal trial design parameters

---

### 5. PoC Calibration Script

Calibrate C_poc threshold programmatically.

**Location**: `src/optimization/poc_calibration_new.R`

**Steps**:

```r
source("src/core/config.R")
source("src/optimization/poc_calibration_new.R")

# 1. Create null/flat scenario
null_scenario <- create_null_flat_scenario(
  n_doses = 5,
  phi_I = 0.25,      # Above threshold
  phi_E = 0.30,      # Above threshold
  tox_upper = 0.30,
  tox_flat = 0.05
)

# 2. Run calibration
calibration_results <- calibrate_c_poc(
  null_scenario = null_scenario,
  c_poc_candidates = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
  n_simulations = 1000,
  base_config = trial_config
)

# 3. Generate report
generate_calibration_report(
  calibration_results = calibration_results,
  null_scenario = null_scenario,
  base_config = trial_config,
  file_path = "results/calibration/calibration_report.txt"
)

# 4. Use optimal C_poc
optimal_c_poc <- calibration_results$optimal_c_poc
```

**Outputs**:
- Calibration results (RData)
- Calibration curve plot
- Detailed text report

**Use case**: Batch calibration, reproducible workflows

---

## Complete Workflow

For production use, follow this complete workflow:

### Phase 1: Initial Setup

1. Clone repository and set working directory
2. Install required R packages (see below)
3. Review `docs/TRIAL_DESIGN.md` for methodology

### Phase 2: Parameter Optimization (Day 1)

1. Run `quick_optimization()` to identify promising parameter regions
2. Analyze results and select top 3-5 combinations
3. Test specific combinations with `test_specific_params()`
4. Update `src/core/config.R` with optimized parameters

### Phase 3: PoC Calibration (Day 2)

1. Open `notebooks/poc_calibration_notebook.qmd`
2. Set null scenario parameters (null_p_I, null_p_E, null_p_T)
3. Run calibration with 1000+ simulations
4. Review calibration report and curves
5. Update config with calibrated C_poc value

### Phase 4: Validation (Day 3)

1. Open `notebooks/simulation_notebook.qmd`
2. Test on multiple scenarios:
   - Signal scenarios (dose 3, 5 optimal)
   - Null scenarios (no dose optimal)
   - Edge cases (all toxic, all ineffective)
3. Verify performance metrics:
   - Type I error ≈ 10% in null scenarios
   - Power ≥ 50% in signal scenarios
   - Early termination < 80%
4. Document final configuration

---

## Required R Packages

```r
# Install required packages
install.packages(c(
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "ggplot2",      # Visualization
  "gridExtra",    # Multiple plots
  "purrr",        # Functional programming
  "isotone",      # PAVA algorithm
  "Iso",          # Isotonic regression
  "copula"        # Copula functions (for data simulation)
))
```

---

## File Structure

### Core Simulation
- `src/core/config.R` - Trial configuration and parameters
- `src/core/main.R` - Main simulation function
- `src/core/simulate_data.R` - Data generation with Gumbel copula
- `src/core/model_utils.R` - Bayesian posterior calculations

### Decision Making
- `src/decision/dose_decision.R` - Admissibility, utility, PoC, early termination

### Optimization
- `src/optimization/parameter_optimization.R` - Parameter search and testing
- `src/optimization/poc_calibration_new.R` - PoC calibration system
- `src/optimization/run_optimization.R` - Optimization runner script

### Utilities
- `src/utils/helpers.R` - Plotting and helper functions
- `src/utils/plotting_extensions.R` - Extended plotting capabilities

### Notebooks
- `notebooks/simulation_notebook.qmd` - Interactive simulation
- `notebooks/poc_calibration_notebook.qmd` - Interactive calibration

### Tests
- `tests/test_main.R` - Main function tests
- `tests/test_dose_decision.R` - Decision logic tests
- `tests/test_workflow_order.R` - Workflow verification
- `tests/test_early_termination_poc.R` - Early termination and PoC tests

---

## Output Files

### Simulation Outputs
```
results/
├── plots/
│   ├── trial_summary_*.png       # Single trial visualizations
│   ├── dose_response_curves.png  # Dose-response relationships
│   └── obd_selection_*.png       # Selection frequency plots
└── Rplots.pdf                     # Additional plots
```

### Calibration Outputs
```
results/calibration/
├── calibration_results.RData      # Complete calibration data
├── calibration_report.txt         # Detailed text report
├── poc_calibration_curve.png      # Calibration curve
└── combined_performance_curves.png # Performance trade-offs
```

### Optimization Outputs
```
results/optimization/
├── quick_optimization_results.RData
├── comprehensive_results.RData
└── optimization_plots/            # Various analysis plots
```

---

## Common Parameters

### Trial Design Parameters

| Parameter | Description | Typical Range | Default |
|-----------|-------------|---------------|---------|
| `dose_levels` | Vector of dose levels | 3-7 doses | c(1,2,3) |
| `n_stages` | Number of trial stages | 2-4 | 3 |
| `cohort_size` | Patients per stage | 10-20 | 15 |

### Threshold Parameters (φ)

| Parameter | Description | Typical Range | Default |
|-----------|-------------|---------------|---------|
| `phi_T` | Toxicity threshold | 0.25-0.45 | 0.30 |
| `phi_E` | Efficacy threshold | 0.10-0.30 | 0.25 |
| `phi_I` | Immune threshold | 0.15-0.35 | 0.20 |

### Credibility Parameters (c)

| Parameter | Description | Typical Range | Default |
|-----------|-------------|---------------|---------|
| `c_T` | Toxicity credibility | 0.7-0.95 | 0.8 |
| `c_E` | Efficacy credibility | 0.6-0.9 | 0.7 |
| `c_I` | Immune credibility | 0.6-0.9 | 0.7 |

### PoC Parameters

| Parameter | Description | Typical Range | Default |
|-----------|-------------|---------------|---------|
| `c_poc` | PoC threshold | 0.5-0.95 | 0.9 (should calibrate) |
| `delta_poc` | PoC comparison threshold | 0.7-0.9 | 0.8 |

### Control Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `enable_early_termination` | Enable early stopping | TRUE |
| `log_early_termination` | Detailed termination logs | TRUE |
| `verbose_logging` | Print detailed progress | TRUE |

---

## Tips and Best Practices

### For Beginners
1. Start with `notebooks/simulation_notebook.qmd`
2. Run single trial with default parameters
3. Modify one parameter at a time
4. Observe how changes affect outcomes

### For Parameter Tuning
1. Use `quick_optimization()` first (2 hours)
2. Focus on promising regions
3. Test specific combinations
4. Document parameter rationale

### For Calibration
1. Set null_p values above thresholds to compensate for PAVA bias
2. Use ≥1000 simulations per C_poc value
3. Review detailed calibration report
4. Target ~10% Type I error rate

### For Large Simulations
1. Set `verbose_logging = FALSE` to reduce output
2. Save intermediate results frequently
3. Use reproducible seeds: `seed = 10000 * i + j`
4. Monitor memory usage

### For Reproducibility
1. Always record parameter settings
2. Use fixed seeds for critical runs
3. Save trial_config with results
4. Version control configuration files

---

## Troubleshooting

### Common Issues

**Issue**: "Cannot find file"
- **Solution**: Ensure working directory is set to project root
- **Check**: `getwd()` should show `.../DoseFinding`
- **Fix**: `setwd("/path/to/DoseFinding")`

**Issue**: "Package not found"
- **Solution**: Install missing packages
- **Check**: `installed.packages()`
- **Fix**: `install.packages("package_name")`

**Issue**: "Trials terminate early too often"
- **Solution**: Parameters too strict
- **Fix**: Lower c values, increase φ_I and φ_E, decrease φ_T

**Issue**: "PoC detection rate too high (>20%)"
- **Solution**: C_poc threshold too low
- **Fix**: Increase C_poc value

**Issue**: "Calibration takes too long"
- **Solution**: Reduce simulations for initial testing
- **Fix**: Start with n_simulations = 100, then increase to 1000+

---

## Further Reading

- **Methodology**: `docs/TRIAL_DESIGN.md` - Complete mathematical framework
- **Optimization Guide**: `docs/OPTIMIZATION_AND_CALIBRATION_GUIDE.md` - Detailed optimization instructions
- **Code Structure**: `docs/CODE_MAP.md` - Code organization overview
- **Naming Conventions**: `NAMING_CONVENTION.md` - Parameter naming guide
- **Recent Updates**: `CALIBRATION_UPDATE_SUMMARY.md` - Latest changes

---

## Support

For questions or issues:
1. Check existing documentation in `docs/`
2. Review test files in `tests/` for examples
3. Examine notebook code for usage patterns
4. Refer to `TRIAL_DESIGN.md` for methodology clarification

## Version Information

- **Last Updated**: October 2025
- **Major Features**: 
  - Multi-stage adaptive Bayesian design
  - Utility-based dose selection
  - Early termination with detailed logging
  - PoC calibration system with Type I error control
  - Parameter optimization framework
  - Interactive notebooks for simulation and calibration
