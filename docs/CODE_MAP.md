# Code Map

This document outlines the structure and purpose of each file in the DoseFinding project.

## Source Code (`src/`)

### Core Logic (`src/core/`)
-   **`main.R`**: Master script containing the `run_trial_simulation()` function that executes the entire multi-stage trial simulation workflow.
    -   Orchestrates stage loops, data generation, posterior updates, admissible set screening, early termination checks, adaptive randomization, and final selection
    -   Returns comprehensive results including allocation history, posterior summaries, and final OD
    -   Evidence: L15-L202
-   **`config.R`**: Configuration parameters for the trial (default: 5 doses, 5 stages, cohort=15)
    -   Defines admissibility thresholds (phi_T, phi_E, phi_I) and credibility cutoffs (c_T, c_E, c_I)
    -   PoC parameters (c_poc, delta_poc) and early termination settings
    -   Utility table (3D array: E × T × I)
    -   Evidence: L8-L64
-   **`simulate_data.R`**: Data generation using Gumbel copula for correlated toxicity/efficacy endpoints
    -   `Gumbel()`: Calculates four-cell probabilities from marginal probabilities and correlation parameter
    -   `simulate_data_gumbel()`: Generates patient-level outcomes (Y_I, Y_T, Y_E) with optional seed
    -   Evidence: L1-L61
-   **`model_utils.R`**: Bayesian model utilities including posterior distributions, isotonic regression, and marginal calculations
    -   `simulate_beta_posterior()`: Beta conjugate posterior sampling (n_sims=1000)
    -   `apply_pava_on_samples()`: Univariate isotonic regression (immune response)
    -   `apply_biviso_on_matrix()`: Bivariate isotonic regression (toxicity, efficacy | immune)
    -   `compute_marginal_probability()`: Marginal probabilities via mixing
    -   Evidence: L1-L128

### Decision Logic (`src/decision/`)
-   **`dose_decision.R`**: Complete decision-making pipeline for admissibility, utility, early termination, and PoC validation
    -   `get_expected_utility()`: Two-layer expected utility calculation (I=0, I=1 scenarios)
    -   `get_admissible_set()`: Filters doses based on posterior probability thresholds
    -   `adaptive_randomization()`: Utility-proportional allocation over admissible set
    -   `check_early_termination()`: Triggers when admissible set empty
    -   `calculate_poc_probability()`: **Posterior sample-based** pairwise comparisons (not normal approximation)
    -   `select_final_od_with_poc()`: Final selection with PoC gating (can return NA)
    -   Evidence: L1-L425

### Optimization (`src/optimization/`)
-   **`poc_calibration_new.R`**: PoC calibration system using null/flat scenarios
    -   `create_null_flat_scenario()`: Constructs scenarios where all doses identical (P_I=φ_I, P_E=φ_E for all doses)
    -   `calibrate_c_poc()`: Tests multiple c_poc candidates with 1000+ simulations each
    -   `plot_calibration_curve()`: Visualizes c_poc vs PoC detection rate
    -   `generate_calibration_report()`: Detailed text report with early termination analysis
    -   Target: ~10% Type I error rate (PoC detection in null scenario)
    -   Evidence: L1-L591
-   **`parameter_optimization.R`**: Systematic parameter search framework
    -   `create_parameter_grids()`: Defines search space for phi_T, phi_E, phi_I, c_T, c_E, c_I, utility variants
    -   `run_parameter_optimization()`: Evaluates parameter combinations across multiple simulations
    -   `create_optimization_plots()`: Visualizes parameter sensitivity and performance metrics
    -   `find_best_parameters()`: Identifies optimal settings by completion rate, selection accuracy, utility
    -   Evidence: L1-L383
-   **`run_optimization.R`**: Convenience wrappers for parameter optimization
    -   `quick_optimization()`: Fast exploration (20 combinations, 3 sims each)
    -   `comprehensive_optimization()`: Thorough search (50 combinations, 5 sims each)
    -   `test_specific_params()`: Test individual parameter sets
    -   Evidence: L1-L173

### Utilities (`src/utils/`)
-   **`helpers.R`**: Core visualization and helper functions
    -   `plot_posterior_summary()`: Posterior means with credible intervals
    -   `compute_rn()`: r/n summaries for Beta posteriors
    -   Additional plotting utilities
-   **`plotting_extensions.R`**: Publication-ready plotting functions
    -   `plot_dose_response_curves()`: Toxicity, efficacy, utility vs dose
    -   `plot_method_comparison_bars()`: Compare methods across scenarios
    -   `plot_multi_scenario_curves()`: Multi-panel scenario comparisons

## Interactive Notebooks (`notebooks/`)
-   **`simulation_notebook.qmd`**: Interactive single trial simulation with visualization
    -   Configures 5-dose trial (aligned with default config)
    -   Runs `run_trial_simulation()` with specified parameters
    -   Generates publication-ready plots (posteriors, allocation, dose-response curves)
    -   Evidence: L1-L369
-   **`poc_calibration_notebook.qmd`**: Interactive PoC calibration workflow
    -   Creates null/flat scenarios using `create_null_flat_scenario()`
    -   Runs `calibrate_c_poc()` across c_poc candidates (0.5-0.95)
    -   Generates calibration curves and detailed reports
    -   Outputs optimal c_poc for ~10% Type I error rate
    -   Evidence: L1-L260

## Testing (`tests/`)
-   **`test_main.R`**: Integration tests for complete trial simulation
    -   Structure/field validation
    -   Allocation probability sums
-   **`test_dose_decision.R`**: Unit tests for decision logic
    -   Expected utility calculations
    -   Admissible set structure
    -   Adaptive allocation normalization
-   **`test_early_termination_poc.R`**: Tests for early termination and PoC scenarios
    -   Empty admissible set handling
    -   PoC validation logic
-   **`test_workflow_order.R`**: Workflow execution order verification
    -   Stage 1 equal allocation
    -   Adaptive allocation in later stages
    -   Early termination timing

## Documentation (`docs/`)
-   **`STATUS_AS_BUILT.md`**: Current implementation status (code-first, English)
    -   What the repo does today (capabilities & limitations)
    -   Canonical run paths (notebooks, scripts, APIs)
    -   Workflow contract (ordering & invariants)
    -   Configuration reality (default vs calibration configs)
    -   Calibration and optimization systems
-   **`STAT_METHODS_AS_BUILT.md`**: Statistical methods as implemented (code-first, Chinese/English)
    -   Data generation (Gumbel copula)
    -   Posterior & isotonic constraints (PAVA/BIVISO)
    -   Decision logic (utility, admissibility, PoC)
    -   Calibration & optimization methodology
-   **`CODE_MAP.md`**: This file - file structure and organization
-   **`HOW_TO_RUN.md`**: Usage instructions and examples
-   **`Design1.tex`, `Design2.tex`**: LaTeX design documents

## Results (`results/`)
-   **`plots/`**: Generated visualization outputs from notebooks
-   **`notebook_calibration/`**: Calibration reports and outputs

## Configuration Files
-   **`.cursorrules`**: Cursor AI coding guidelines and project context
-   **`DoseFinding.Rproj`**: RStudio project file