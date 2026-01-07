# Bayesian Adaptive Dose-Finding Trial Simulation

A comprehensive R implementation of Bayesian adaptive dose-finding clinical trials featuring multi-stage design, utility-based decision making, isotonic regression constraints, Probability of Correct Selection (PoC) calibration, and systematic parameter optimization.

## Overview

This package simulates adaptive Bayesian dose-finding trials with three correlated endpoints (immune response, toxicity, efficacy) using:
- **Gumbel copula** for data generation with realistic correlation structure
- **PAVA/BIVISO** isotonic regression to enforce monotone dose-response constraints
- **Utility-based adaptive randomization** to allocate patients to promising doses
- **Admissibility screening** to exclude unsafe/ineffective doses
- **PoC validation** for final dose selection with Type I error control
- **Early termination** when no acceptable doses remain

**Default configuration**: 5 doses, 5 stages, 15 patients per stage (75 total)

## Quick Start

### Prerequisites

Install required R packages:
```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "purrr",
  "isotone", "Iso", "gridExtra"
))
```

### Option 1: Interactive Simulation (Recommended for First-Time Users)

1. **Open RStudio** and set working directory to project root:
   ```r
   setwd("/path/to/DoseFinding")
   ```

2. **Open notebook**: `notebooks/simulation_notebook.qmd`

3. **Run chunks sequentially** to:
   - Configure trial parameters (dose levels, stages, thresholds)
   - Define scenario probabilities (immune response, toxicity, efficacy)
   - Execute single trial simulation
   - Visualize results (posterior distributions, allocation patterns, dose-response curves)

**When to use**: Learning the system, testing single scenarios, exploratory analysis

### Option 2: PoC Calibration (Required for Production Use)

1. **Open notebook**: `notebooks/poc_calibration_notebook.qmd`

2. **Configure null scenario** (all doses equally ineffective):
   ```r
   null_p_I <- 0.25   # Immune response probability
   null_p_E <- 0.30   # Efficacy probability
   null_p_T <- 0.05   # Toxicity probability
   ```

3. **Run calibration** across c_poc candidates (0.5, 0.6, ..., 0.95)

4. **Review outputs**:
   - Calibration curve (c_poc vs detection rate)
   - Detailed report with early termination analysis
   - Recommended c_poc value achieving ~10% Type I error

5. **Update configuration** with calibrated c_poc

**When to use**: Before running production simulations, to control false positive rate

### Option 3: Direct R Script Execution

```r
# Set working directory
setwd("/path/to/DoseFinding")

# Source configuration and main script
source("src/core/config.R")
source("src/core/main.R")

# Define scenario (dose 3 optimal)
p_YI <- c(0.10, 0.30, 0.50, 0.60, 0.70)

p_YT_given_I <- matrix(c(
  0.05, 0.10, 0.12, 0.18, 0.25,  # I=0
  0.08, 0.12, 0.15, 0.25, 0.35   # I=1
), nrow=5, ncol=2)

p_YE_given_I <- matrix(c(
  0.10, 0.20, 0.35, 0.45, 0.50,  # I=0
  0.30, 0.50, 0.70, 0.80, 0.75   # I=1
), nrow=5, ncol=2)

# Run simulation
result <- run_trial_simulation(
  trial_config = trial_config,
  p_YI = p_YI,
  p_YT_given_I = p_YT_given_I,
  p_YE_given_I = p_YE_given_I,
  rho0 = 1.5,
  rho1 = 2.0,
  seed = 123
)

# View results
result$final_od              # Selected optimal dose
result$final_utility         # Expected utility
result$poc_validated         # PoC threshold met?
result$terminated_early      # Early termination flag
```

**When to use**: Batch simulations, automation, integration with other tools

## Key Features

### 1. **Multi-Stage Bayesian Adaptive Design**
- Interim analyses after each stage
- Adaptive allocation based on accumulating data
- Early termination when admissible set becomes empty
- Stage-specific random seeds for reproducibility

### 2. **Three Correlated Endpoints**
- **Immune response (I)**: Binary indicator of immunogenic activity
- **Toxicity (T)**: Binary indicator of dose-limiting toxicity
- **Efficacy (E)**: Binary indicator of clinical efficacy
- Generated via Gumbel copula with dose-specific correlations (ρ₀, ρ₁)

### 3. **Bayesian Posterior Updates with Isotonic Constraints**
- Beta conjugate priors for computational efficiency
- **PAVA** (Pool Adjacent Violators Algorithm) for univariate monotonicity (immune response)
- **BIVISO** (Bivariate Isotonic Regression) for conditional probabilities (T|I, E|I)
- Marginal probabilities computed via mixing over immune response distribution

### 4. **Admissibility Screening**
- Doses must satisfy three criteria:
  - **Toxicity**: P(π_T ≤ φ_T | data) ≥ c_T
  - **Efficacy**: P(π_E ≥ φ_E | data) ≥ c_E
  - **Immune response**: P(π_I ≥ φ_I | data) ≥ c_I
- Empty admissible set triggers early termination

### 5. **Utility-Based Decision Making**
- 3D utility table: U[E, T, I] defines value for each outcome combination
- Expected utility calculated using posterior probabilities
- Adaptive allocation proportional to utility scores
- Final dose selection maximizes expected utility among admissible doses

### 6. **Probability of Correct Selection (PoC) Validation**
- **Method**: Posterior sample-based pairwise comparisons
- **Formula**: PoC = min{P(π_best > δ × π_competitor | data)} across all competing doses
- **Gating**: If PoC < c_poc, no optimal dose is selected (returns NA)
- **Calibration**: c_poc threshold determined via null scenario simulations to achieve ~10% Type I error

### 7. **PoC Calibration System**
- **Null/flat scenarios**: All doses have identical response probabilities
- **Procedure**: Test c_poc candidates (0.5, 0.6, ..., 0.95) with 1000+ simulations each
- **Target**: ~10% false positive rate (PoC detection in null scenario)
- **Outputs**: Calibration curves, detailed reports, optimal c_poc recommendation
- **Implementation**: `src/optimization/poc_calibration_new.R`

### 8. **Parameter Optimization Framework**
- **Search space**: φ_T, φ_E, φ_I (thresholds), c_T, c_E, c_I (credibility cutoffs)
- **Metrics**: Completion rate, correct selection rate, mean utility, allocation efficiency
- **Methods**: Grid search, random sampling, specific parameter testing
- **Outputs**: Performance plots, sensitivity analysis, top parameter combinations
- **Implementation**: `src/optimization/parameter_optimization.R`

### 9. **Interactive Notebooks**
- **Simulation notebook**: Single trial visualization and analysis
- **Calibration notebook**: PoC threshold calibration workflow
- Both built with Quarto for reproducible research

### 10. **Comprehensive Testing**
- Unit tests for decision logic, posterior calculations, utility functions
- Integration tests for complete workflow
- Workflow order verification (early termination timing, allocation logic)
- Edge case handling (empty admissible sets, extreme probabilities)

## 🚀 Quick Start

### Option 1: Interactive Notebook (Recommended)
```r
# Open notebooks/simulation_notebook.qmd
# Contains complete interactive examples and calibration framework
```

### Option 2: Direct Script Execution
```r
# Load core functions
source("src/core/config.R")
source("src/core/main.R")

# Run trial simulation
results <- run_trial_simulation(
  trial_config = trial_config,
  p_YI = p_YI,
  p_YT_given_I = p_YT_given_I,
  p_YE_given_I = p_YE_given_I,
  rho0 = rho0,
  rho1 = rho1
)
```

### Option 3: Calibration Demo
```r
# Run comprehensive calibration demo
source("examples/comprehensive_calibration_demo.R")
```

## 📁 Project Structure

```
DoseFinding/
├── src/                              # Source code
│   ├── core/                         # Core simulation engine
│   │   ├── main.R                    # Master workflow (run_trial_simulation)
│   │   ├── config.R                  # Trial configuration (5-dose, 5-stage default)
│   │   ├── simulate_data.R           # Gumbel copula data generation
│   │   └── model_utils.R             # Bayesian posteriors + PAVA/BIVISO
│   ├── decision/                     # Decision-making logic
<<<<<<< HEAD
│   │   └── dose_decision.R           # Admissibility, utility, PoC, early termination
│   ├── optimization/                 # Calibration & optimization tools
│   │   ├── poc_calibration_new.R     # PoC calibration system
│   │   ├── parameter_optimization.R  # Parameter search framework
│   │   └── run_optimization.R        # Optimization convenience wrappers
│   └── utils/                        # Visualization & helpers
│       ├── helpers.R                 # Core plotting functions
│       └── plotting_extensions.R     # Publication-ready plots
├── notebooks/                        # Interactive Quarto notebooks
│   ├── simulation_notebook.qmd       # Single trial simulation
│   ├── poc_calibration_notebook.qmd  # PoC calibration workflow
│   └── design_walkthrough.qmd        # Design documentation walkthrough
├── tests/                            # Test suite
│   ├── test_main.R                   # Integration tests
│   ├── test_dose_decision.R          # Decision logic tests
│   ├── test_early_termination_poc.R  # Early termination & PoC tests
│   └── test_workflow_order.R         # Workflow sequence verification
├── docs/                             # Documentation
│   ├── STATUS_AS_BUILT.md            # Current implementation status
│   ├── STAT_METHODS_AS_BUILT.md      # Statistical methods (bilingual)
│   ├── CODE_MAP.md                   # File structure guide
│   ├── HOW_TO_RUN.md                 # Detailed usage instructions
│   ├── Design1.tex                   # Design specification (LaTeX)
│   └── Design2.tex                   # Design specification (LaTeX)
├── examples/                         # Example scripts
│   └── plotting_demo.R               # Plotting examples
└── DoseFinding.Rproj                 # RStudio project file
```

## Documentation

Comprehensive documentation is available in the `docs/` directory:

- **[HOW_TO_RUN.md](docs/HOW_TO_RUN.md)** - Complete usage guide with examples for all workflows
- **[STATUS_AS_BUILT.md](docs/STATUS_AS_BUILT.md)** - Current implementation status, API contracts, known limitations
- **[STAT_METHODS_AS_BUILT.md](docs/STAT_METHODS_AS_BUILT.md)** - Statistical methods with code evidence (Chinese/English)
- **[CODE_MAP.md](docs/CODE_MAP.md)** - File structure and organization
- **[Design1.tex](docs/Design1.tex), [Design2.tex](docs/Design2.tex)** - Mathematical design specifications
=======
│   │   └── dose_decision.R           # Dose selection algorithms
│   ├── optimization/                 # Parameter optimization
│   │   ├── poc_calibration.R         # PoC calibration framework
│   │   └── early_termination_calibration.R  # Early termination calibration
│   └── utils/                        # Utility functions
│       ├── helpers.R                 # Helper functions and plotting
│       ├── plotting_extensions.R     # Plotting extensions
│       └── calibration_plots.R       # Calibration visualization
├── examples/                         # Example scripts
│   ├── comprehensive_calibration_demo.R  # Comprehensive calibration demo
│   ├── poc_calibration_demo.R        # PoC calibration demo
│   ├── flat_scenario_demo.R          # Flat scenario demo
│   └── bayesian_poc_demo.R           # Bayesian PoC demo
├── tests/                            # Test files
│   ├── test_comprehensive_calibration.R  # Comprehensive calibration tests
│   ├── test_poc_calibration.R        # PoC calibration tests
│   └── test_*.R                      # Other test files
├── docs/                             # Documentation
│   ├── PROJECT_OVERVIEW.md           # Project overview
│   ├── TRIAL_DESIGN.md               # Trial design specification
│   ├── NEXT_STEP_PLAN.md             # Implementation plan
│   └── CALIBRATION_IMPLEMENTATION_SUMMARY.md  # Calibration implementation summary
├── notebooks/                        # Interactive notebooks
│   └── simulation_notebook.qmd       # Interactive simulation notebook
└── results/                          # Generated outputs
    └── plots/                        # Generated plots
```

## ✨ Key Features

### 🎯 Trial Simulation
- **Multi-stage Bayesian adaptive design** with interim analyses
- **Utility-based dose selection** with customizable utility functions
- **Early termination criteria** for safety and efficacy
- **PoC validation** probability of correct selection validation
- **Adaptive randomization** based on utility scores

### 🔧 Calibration System
- **PoC calibration** target: 10% detection rate
- **Early termination calibration** target: 80% termination rate
- **Performance visualization** calibration curves and confidence intervals
- **Parameter optimization** systematic parameter tuning

### 📊 Visualization
- **Dose-response curves** toxicity, efficacy, and utility
- **Posterior distribution plots** modern styling
- **Calibration curves** threshold vs performance relationships
- **Allocation analysis** participant distribution

## 📚 Documentation

### Quick Start
- **QUICK_START.md** - 5-minute quick start guide
- **PROJECT_OVERVIEW.md** - Complete project overview and usage

### Detailed Documentation
- **TRIAL_DESIGN.md** - Complete trial design specification
- **NEXT_STEP_PLAN.md** - Implementation plan and status
- **CALIBRATION_IMPLEMENTATION_SUMMARY.md** - Calibration system implementation summary

## 🛠️ Requirements
>>>>>>> origin/main

## Complete Workflow (Production Use)

For production simulations, follow this recommended workflow:

### Phase 1: Parameter Optimization (Day 1)
```r
source("src/optimization/run_optimization.R")

# Quick exploration (2 hours)
results <- quick_optimization()
analyze_results(results)

# Select promising parameter regions
# Test specific combinations
results <- test_specific_params(
  phi_T = 0.30, phi_E = 0.15, phi_I = 0.20,
  c_T = 0.8, c_E = 0.7, c_I = 0.7
)
```

### Phase 2: PoC Calibration (Day 2)
```r
# Open notebooks/poc_calibration_notebook.qmd
# 1. Configure null scenario (all doses identical)
# 2. Run calibration across c_poc candidates
# 3. Review calibration report
# 4. Select c_poc achieving ~10% Type I error
# 5. Update trial_config with calibrated c_poc
```

### Phase 3: Validation (Day 3)
```r
# Open notebooks/simulation_notebook.qmd
# Test calibrated parameters on multiple scenarios:
# - Signal scenarios (dose 3, 5 optimal)
# - Null scenarios (no dose optimal)
# - Edge cases (all toxic, all ineffective)
# 
# Verify:
# - Type I error ≈ 10% in null scenarios
# - Power ≥ 50% in signal scenarios
# - Reasonable early termination rate (<80%)
```

## Key Configuration Parameters

### Trial Design
- `dose_levels`: Vector of dose labels (default: c(1,2,3,4,5))
- `n_stages`: Number of trial stages (default: 5)
- `cohort_size`: Patients per stage (default: 15)

### Admissibility Thresholds (φ)
- `phi_T`: Maximum acceptable toxicity probability (default: 0.35)
- `phi_E`: Minimum required efficacy probability (default: 0.10)
- `phi_I`: Minimum required immune response probability (default: 0.20)

### Credibility Cutoffs (c)
- `c_T`: Posterior probability threshold for toxicity safety (default: 0.5)
- `c_E`: Posterior probability threshold for efficacy (default: 0.5)
- `c_I`: Posterior probability threshold for immune response (default: 0.5)

### PoC Parameters
- `c_poc`: PoC threshold for final selection (default: 0.9, **should be calibrated**)
- `delta_poc`: Threshold for pairwise comparisons (default: 0.8)

### Control Flags
- `enable_early_termination`: Enable early stopping (default: TRUE)
- `log_early_termination`: Detailed termination logging (default: TRUE)
- `verbose_logging`: Print progress during simulation (default: TRUE, set FALSE for batch runs)

## Dependencies

- **R** (>= 4.0)
- **dplyr** - Data manipulation
- **tidyr** - Data tidying
- **ggplot2** - Visualization
- **gridExtra** - Multi-panel plots
- **purrr** - Functional programming
- **isotone** - PAVA algorithm
- **Iso** - Bivariate isotonic regression

## Current Implementation Status

### ✅ Fully Implemented
- Multi-stage Bayesian adaptive design (5 doses, 5 stages)
- Three correlated endpoints via Gumbel copula
- PAVA/BIVISO isotonic regression
- Admissibility screening with posterior probabilities
- Utility-based adaptive allocation
- Early termination when admissible set empty
- PoC validation using posterior samples (not normal approximation)
- PoC calibration system with null/flat scenarios
- Parameter optimization framework
- Interactive notebooks for simulation and calibration
- Comprehensive test suite
- Stage-specific seed management for reproducibility

### ⚠️ Known Limitations
- No control arm or γ_j logic (marked out-of-scope)
- No multiple comparison adjustment in admissibility screening
- No parallel computation support
- Calibration requires manual validation in signal scenarios

### 🔄 Recent Updates
- PoC changed from normal approximation to posterior sample-based method
- PoC can now gate selection (return NA if threshold not met)
- Configuration alignment: default now matches simulation notebook (5-dose)
- Added verbose_logging control for batch calibration runs
- Comprehensive calibration system with detailed reporting

## Testing

Run the test suite from project root:

```r
# Run all tests
source("tests/test_main.R")
source("tests/test_dose_decision.R")
source("tests/test_early_termination_poc.R")
source("tests/test_workflow_order.R")

# Or using testthat
testthat::test_dir("tests")
```

## Common Issues & Troubleshooting

**Issue**: "Cannot find file"
- **Solution**: Ensure working directory is set to project root: `setwd("/path/to/DoseFinding")`

**Issue**: Trials terminate early too often
- **Solution**: Parameters too strict. Lower c values, increase φ_E and φ_I, decrease φ_T

**Issue**: PoC detection rate too high in null scenarios (>20%)
- **Solution**: c_poc threshold too low. Increase c_poc value via calibration

**Issue**: Calibration takes too long
- **Solution**: Reduce n_simulations to 100 for testing, then increase to 1000+ for final calibration

## Citation

If you use this software in your research, please cite:

```
[Citation information to be added]
```

## 📋 Example Scripts

| Script | Function | Runtime |
|--------|----------|---------|
| `examples/plotting_demo.R` | Basic simulation and visualization | 1 minute |
| `examples/poc_calibration_demo.R` | PoC calibration demo | 5 minutes |
| `examples/comprehensive_calibration_demo.R` | Complete calibration system | 10 minutes |
| `notebooks/simulation_notebook.qmd` | Interactive notebook | Variable |

## 🧪 Testing

```r
# Run all tests
source("tests/test_comprehensive_calibration.R")

# Run specific tests
source("tests/test_poc_calibration.R")
source("tests/test_early_termination_poc.R")
```

## 📈 Project Status

### ✅ Completion Status (100%)
- ✅ Complete trial simulation workflow
- ✅ Bayesian posterior probability calculation
- ✅ Adaptive randomization algorithm
- ✅ Early termination mechanism
- ✅ PoC validation system
- ✅ Comprehensive calibration framework
- ✅ Performance visualization tools
- ✅ Complete integration testing

## 🤝 Support

For questions, please check:
1. **Quick Start**: `QUICK_START.md`
2. **Project Overview**: `docs/PROJECT_OVERVIEW.md`
3. **Example Scripts**: `examples/` directory
4. **Interactive Notebook**: `notebooks/` directory
5. **Test Files**: `tests/` directory

## 📄 License

This project is for research and educational purposes.

## Contact & Support

For questions or issues:
1. Check documentation in `docs/`
2. Review test files in `tests/` for usage examples
3. Examine notebooks for interactive workflows
4. Refer to `docs/STATUS_AS_BUILT.md` for implementation details

---

**Last Updated**: January 2026  
**Version**: 1.0 (Code-first documentation, aligned with implementation)
