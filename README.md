# Bayesian Adaptive Dose-Finding Trial Simulation

This R project simulates Bayesian adaptive dose-finding trials with three binary endpoints:

- immune response (`Y_I`)
- toxicity (`Y_T`)
- efficacy (`Y_E`)

The current design uses a 5-dose, 5-stage trial by default, with 15 patients per stage. The simulator updates Bayesian posterior estimates after each stage, screens out unacceptable doses, adaptively allocates future patients by expected utility, can stop early when no dose remains acceptable, and applies a final Probability of Correct Selection (PoC) gate before selecting an optimal dose.

## Current Status

The project has been cleaned up and is runnable from the repository root.

- Canonical documentation starts at [docs/README.md](docs/README.md) and [docs/CHATGPT_PROJECT_CONTEXT.md](docs/CHATGPT_PROJECT_CONTEXT.md).
- Canonical PoC calibration code is [src/optimization/poc_calibration.R](src/optimization/poc_calibration.R).
- Generated notebooks, plots, calibration output, and RStudio state are ignored by git.
- Latest full local verification: `testthat::test_dir("tests")` passed with `FAIL 0 | WARN 0 | SKIP 0 | PASS 399`.

The main remaining work is statistical design confirmation, not basic code cleanup.

Before adding major new features, confirm with the project lead:

1. Should the current no-control-arm implementation be treated as the official design, or should the control arm and `gamma_j` allocation idea from `Design2.tex` be implemented?
2. Should PoC calibration target the design/PDF goal of familywise Type I error `0.05`, or the current code's null/flat PoC detection target `0.10`?
3. What should final PoC compare: immune response, marginal efficacy, expected utility, or another endpoint?
4. Which parameters should be production-calibrated: only `c_poc`, or also `c_T`, `c_E`, and `c_I`?
5. Which operating characteristics should appear in the final report?

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

Important: the example below uses the current code's `target_rate = 0.10` convention for null/flat PoC detection. The original design notes/PDF mention familywise Type I error `0.05`, so the final production target should be confirmed before reporting final calibration results.

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

Latest full local run: `FAIL 0 | WARN 0 | SKIP 0 | PASS 399`.
