# Bayesian Adaptive Dose-Finding Trial Simulation

This R project simulates Bayesian adaptive dose-finding trials with three binary endpoints:

- immune response (`Y_I`)
- toxicity (`Y_T`)
- efficacy (`Y_E`)

The current design uses a 5-dose, 5-stage trial by default, with 15 patients per stage. The simulator updates Bayesian posterior estimates after each stage, screens out unacceptable doses, adaptively allocates future patients by expected utility, can stop early when no dose remains acceptable, and applies a final Probability of Correct Selection (PoC) gate before selecting an optimal dose.

## Start Here

For a full project re-orientation, read:

- [docs/CHATGPT_PROJECT_CONTEXT.md](docs/CHATGPT_PROJECT_CONTEXT.md): the canonical context packet for project Q&A
- [docs/PROJECT_READALOUD_EXPLANATION.md](docs/PROJECT_READALOUD_EXPLANATION.md): plain-language Chinese read-aloud explanation
- [docs/PROJECT_READALOUD_MOBILE.html](docs/PROJECT_READALOUD_MOBILE.html): phone-friendly read-aloud version
- [docs/HOW_TO_RUN.md](docs/HOW_TO_RUN.md): running simulations, calibration, optimization, and notebooks
- [docs/STAT_METHODS_AS_BUILT.md](docs/STAT_METHODS_AS_BUILT.md): statistical methods as implemented
- [docs/CODE_MAP.md](docs/CODE_MAP.md): file-by-file code map

The original design notes are kept in [docs/Design1.tex](docs/Design1.tex) and [docs/Design2.tex](docs/Design2.tex).

## Install Packages

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "purrr",
  "isotone", "Iso", "gridExtra", "testthat"
))
```

## Run One Trial

```r
setwd("/path/to/DoseFinding")

source("src/core/config.R")
source("src/core/main.R")

result <- run_trial_simulation(
  trial_config = trial_config,
  p_YI = p_YI,
  p_YT_given_I = p_YT_given_I,
  p_YE_given_I = p_YE_given_I,
  rho0 = rho0,
  rho1 = rho1,
  seed = 123
)

result$final_od
result$poc_validated
result$terminated_early
```

## Calibrate PoC

`src/optimization/poc_calibration.R` is the canonical PoC calibration implementation.

```r
source("src/core/config.R")
source("src/optimization/poc_calibration.R")

null_scenario <- create_null_flat_scenario(
  n_doses = 5,
  phi_I = 0.20,
  phi_E = 0.25,
  tox_flat = 0.05
)

calibration_results <- calibrate_c_poc(
  null_scenario = null_scenario,
  c_poc_candidates = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
  n_simulations = 1000,
  base_config = trial_config,
  target_rate = 0.10
)

calibration_results$optimal_c_poc
```

## Project Structure

```text
DoseFinding/
├── src/
│   ├── core/
│   │   ├── config.R
│   │   ├── main.R
│   │   ├── model_utils.R
│   │   └── simulate_data.R
│   ├── decision/
│   │   └── dose_decision.R
│   ├── optimization/
│   │   ├── poc_calibration.R
│   │   ├── early_termination_calibration.R
│   │   ├── parameter_optimization.R
│   │   └── run_optimization.R
│   └── utils/
│       ├── calibration_plots.R
│       ├── helpers.R
│       └── plotting_extensions.R
├── notebooks/
│   ├── simulation_notebook.qmd
│   ├── poc_calibration_notebook.qmd
│   ├── threshold_calibration_notebook.qmd
│   └── design_walkthrough.qmd
├── docs/
├── examples/
└── tests/
```

Generated notebooks, plots, calibration reports, RStudio state, and result files are intentionally ignored by git.

## Tests

```r
testthat::test_dir("tests")
```

Some calibration tests run small simulations and may take longer than ordinary unit tests.
