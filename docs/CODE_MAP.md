# Code Map

This document outlines the structure and purpose of each file in the DoseFinding project.

## Source Code (`src/`)

### Core Logic (`src/core/`)
-   `main.R`: Master script containing the `run_trial_simulation` function that executes the entire multi-stage trial simulation.
-   `config.R`: Configuration parameters for the trial, including dose levels, number of stages, cohort size, and the utility table.
-   `simulate_data.R`: Contains the `simulate_data_gumbel` function for generating patient data for each stage.
-   `model_utils.R`: Bayesian model utilities including posterior distributions, isotonic regression, and marginal probabilities.

### Decision Logic (`src/decision/`)
-   `dose_decision.R`: Functions for decision-making at each stage, including admissible set identification, expected utility calculation, and final optimal dose selection.

### Optimization (`src/optimization/`)
-   `parameter_optimization.R`: Systematic parameter optimization algorithms for trial design tuning.
-   `run_optimization.R`: Interface for running parameter optimization with different configurations.

### Utilities (`src/utils/`)
-   `helpers.R`: Helper functions, most notably `plot_posterior_summary` for visualizing results.

## Interactive Simulation (`notebooks/`)
-   `simulation_notebook.qmd`: Interactive Quarto notebook for easy configuration and execution of trial simulations.

## Testing (`tests/`)
-   `test_main.R`: Tests for the main simulation script.
-   `test_dose_decision.R`: Unit tests for dose decision logic.
-   `test_early_termination_poc.R`: Tests for early termination functionality.
-   `test_workflow_order.R`: Tests for workflow execution order.

## Documentation (`docs/`)
-   `DESIGN_NOTES.md`: Implementation details and design notes.
-   `TRIAL_DESIGN.md`: Complete trial design specification.
-   `PARAMETER_OPTIMIZATION_GUIDE.md`: Guide for parameter optimization.
-   `UTILITY_CALCULATION_GUIDE.md`: Guide for utility function calculations.
-   `HOW_TO_RUN.md`: Usage instructions and examples.

## Results (`results/`)
-   `plots/`: Generated visualization outputs.
-   `Rplots.pdf`: Additional plot outputs.