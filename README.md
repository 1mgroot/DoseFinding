# Bayesian Adaptive Dose-Finding Trial Simulation

A comprehensive R implementation of Bayesian adaptive dose-finding trials with multi-stage design, utility-based decision making, and parameter optimization capabilities.

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
│   ├── core/                         # Core simulation logic
│   │   ├── main.R                    # Master simulation script
│   │   ├── config.R                  # Trial configuration parameters
│   │   ├── simulate_data.R           # Data simulation functions
│   │   └── model_utils.R             # Bayesian model utilities
│   ├── decision/                     # Decision-making logic
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

- R (>= 4.0)
- Required packages: dplyr, tidyr, isotone, purrr, ggplot2, Iso, testthat

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
