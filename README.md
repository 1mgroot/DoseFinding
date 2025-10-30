# Bayesian Adaptive Dose-Finding Trial Simulation

A comprehensive R implementation of Bayesian adaptive dose-finding trials with multi-stage design, utility-based decision making, and parameter optimization capabilities.

## Project Structure

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
│   │   ├── parameter_optimization.R  # Optimization algorithms
│   │   ├── run_optimization.R       # Optimization runner
│   │   └── poc_calibration_new.R    # PoC calibration functions
│   └── utils/                        # Utility functions
│       ├── helpers.R                 # Helper functions and plotting
│       └── plotting_extensions.R    # Extended plotting functions
├── tests/                            # Test files
│   ├── test_main.R                   # Main script tests
│   ├── test_dose_decision.R          # Decision logic tests
│   ├── test_early_termination_poc.R  # Early termination tests
│   └── test_workflow_order.R         # Workflow tests
├── docs/                             # Documentation
│   ├── DESIGN_NOTES.md               # Design implementation notes
│   ├── TRIAL_DESIGN.md               # Trial design specification
│   ├── TRIAL_DESIGN.html             # HTML version of design doc
│   ├── OPTIMIZATION_AND_CALIBRATION_GUIDE.md # Complete optimization & calibration guide
│   ├── UTILITY_CALCULATION_GUIDE.md  # Utility calculation guide
│   ├── CODE_MAP.md                   # Code structure overview
│   └── HOW_TO_RUN.md                 # Usage instructions
├── notebooks/                        # Interactive notebooks
│   ├── simulation_notebook.qmd      # Interactive simulation notebook
│   └── poc_calibration_notebook.qmd # PoC calibration notebook
└── results/                          # Generated outputs
    ├── plots/                        # Generated plots
    └── Rplots.pdf                    # Additional plots
```

## Quick Start

### Option 1: Interactive Notebook (Recommended)
1. Open `notebooks/simulation_notebook.qmd` in RStudio
2. Modify configuration parameters as needed
3. Run the notebook chunks sequentially

### Option 2: Direct Script Execution
1. Open `src/core/main.R` in RStudio
2. Source the script: `source("src/core/main.R")`

### Option 3: PoC Calibration Workflow
1. Open `notebooks/poc_calibration_notebook.qmd` in RStudio
2. Run the calibration workflow to find optimal C_poc threshold
3. Use calibrated parameters in your trial simulations

## Key Features

- **Multi-stage Bayesian adaptive design** with interim analyses
- **Utility-based dose selection** with customizable utility functions
- **Early termination criteria** for safety and efficacy
- **PoC calibration system** using null/flat scenarios for Type I error control
- **Parameter optimization** for trial design tuning
- **Comprehensive testing suite** with unit tests
- **Interactive simulation notebooks** for easy exploration
- **Extensive documentation** with design specifications

## Documentation

See the `docs/` directory for comprehensive documentation:
- `HOW_TO_RUN.md` - Detailed usage instructions and quick start guide
- `OPTIMIZATION_AND_CALIBRATION_GUIDE.md` - Complete parameter optimization and PoC calibration guide
- `TRIAL_DESIGN.md` - Complete trial design specification and mathematical framework
- `DESIGN_NOTES.md` - Implementation details and development notes
- `UTILITY_CALCULATION_GUIDE.md` - Utility function guide and examples
- `CODE_MAP.md` - Code structure overview

## Requirements

- R (>= 4.0)
- Required packages: dplyr, tidyr, isotone, purrr, ggplot2, Iso, testthat

## License

This project is for research and educational purposes.
