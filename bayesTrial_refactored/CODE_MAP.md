# Code Map

This document outlines the structure and purpose of each file in the `bayesTrial_refactored` directory.

## Core Scripts

-   `main.R`: This is the master script. It has been refactored to contain a primary function, `run_trial_simulation`, which executes the entire multi-stage trial simulation. It can be sourced by other scripts (like the notebook) or run directly for a test simulation.
-   `config.R`: Contains all the configuration parameters for the trial, including dose levels, number of stages, cohort size, and the utility table.
-   `simulate_data.R`: Contains the `simulate_data_gumbel` function, which generates patient data for each stage of the trial based on the predefined probabilities.
-   `model_utils.R`: A collection of utility functions for the Bayesian model, including functions to calculate posterior distributions (`simulate_beta_posterior`), apply isotonic regression (`apply_pava_on_samples`, `apply_biviso_on_matrix`), and compute marginal probabilities.
-   `dose_decision.R`: Contains all functions related to the decision-making process at the end of each stage, such as identifying the admissible set (`get_admissible_set`), calculating expected utility (`get_expected_utility`), and selecting the final optimal dose (`select_final_od`).
-   `helpers.R`: Contains helper functions, most notably `plot_posterior_summary` for visualizing the results.

## Interactive Simulation & Reporting

-   `simulation_notebook.qmd`: An interactive Quarto notebook that allows for easy configuration and execution of the trial simulation. It calls the `run_trial_simulation` function from `main.R` and provides a clean, user-friendly way to visualize the results.

## Testing

-   `test_dose_decision.R`: Contains unit tests for the dose decision logic.
-   `test_main.R`: A simple test script to ensure that the `main.R` script runs without errors.