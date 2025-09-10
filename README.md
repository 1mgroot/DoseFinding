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
│   │   └── run_optimization.R       # Optimization runner
│   └── utils/                        # Utility functions
│       └── helpers.R                 # Helper functions and plotting
├── tests/                            # Test files
│   ├── test_main.R                   # Main script tests
│   ├── test_dose_decision.R          # Decision logic tests
│   ├── test_early_termination_poc.R  # Early termination tests
│   └── test_workflow_order.R         # Workflow tests
├── docs/                             # Documentation
│   ├── DESIGN_NOTES.md               # Design implementation notes
│   ├── TRIAL_DESIGN.md               # Trial design specification
│   ├── TRIAL_DESIGN.html             # HTML version of design doc
│   ├── PARAMETER_OPTIMIZATION_GUIDE.md # Optimization guide
│   ├── UTILITY_CALCULATION_GUIDE.md  # Utility calculation guide
│   ├── CODE_MAP.md                   # Code structure overview
│   └── HOW_TO_RUN.md                 # Usage instructions
├── notebooks/                        # Interactive notebooks
│   └── simulation_notebook.qmd      # Interactive simulation notebook
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

## Key Features

- **Multi-stage Bayesian adaptive design** with interim analyses
- **Utility-based dose selection** with customizable utility functions
- **Early termination criteria** for safety and efficacy
- **Parameter optimization** for trial design tuning
- **Comprehensive testing suite** with unit tests
- **Interactive simulation notebook** for easy exploration
- **Extensive documentation** with design specifications

## Documentation

See the `docs/` directory for comprehensive documentation:
- `TRIAL_DESIGN.md` - Complete trial design specification
- `DESIGN_NOTES.md` - Implementation details and notes
- `HOW_TO_RUN.md` - Detailed usage instructions
- `PARAMETER_OPTIMIZATION_GUIDE.md` - Optimization guide
- `UTILITY_CALCULATION_GUIDE.md` - Utility function guide

## Requirements

- R (>= 4.0)
- Required packages: dplyr, tidyr, isotone, purrr, ggplot2, Iso, testthat

## License

This project is for research and educational purposes.
