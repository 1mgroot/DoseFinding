# Bayesian Adaptive Dose-Finding Trial Simulation

This repository contains an R simulation framework for Bayesian adaptive
dose-finding trials with three binary endpoints:

- immune response (`Y_I`)
- toxicity (`Y_T`)
- efficacy (`Y_E`)

The default design uses 5 doses, 5 stages, and 15 patients per stage. During a
simulated trial, the program learns from accumulated patient outcomes, updates
Bayesian posterior estimates, removes unacceptable doses, adaptively allocates
future patients toward better doses, stops early when no dose remains acceptable,
and applies a final Probability of Correct Selection (PoC) gate before selecting
an optimal dose (OD).

In plain terms: this project repeatedly simulates a clinical trial so we can
study which dose looks best after balancing safety, efficacy, immune response,
and strength of evidence.

## Current Status

The project is cleaned up and runnable from the repository root.

- Latest full local verification: `testthat::test_dir("tests")` passed with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 399`.
- Canonical project context is in
  [docs/CHATGPT_PROJECT_CONTEXT.md](docs/CHATGPT_PROJECT_CONTEXT.md).
- Documentation index is in [docs/README.md](docs/README.md).
- Main simulation entry point is `run_trial_simulation()` in
  [src/core/main.R](src/core/main.R).
- Canonical PoC calibration implementation is
  [src/optimization/poc_calibration.R](src/optimization/poc_calibration.R).
- Generated notebooks, plots, calibration outputs, result files, and RStudio
  state are ignored by git.

The main remaining work is statistical design confirmation and production-scale
operating characteristic simulation, not basic code repair.

## Start Here

For re-orientation, read these in order:

1. [docs/CHATGPT_PROJECT_CONTEXT.md](docs/CHATGPT_PROJECT_CONTEXT.md) - full
   project context packet, suitable to give another ChatGPT session.
2. [docs/PROJECT_READALOUD_EXPLANATION.md](docs/PROJECT_READALOUD_EXPLANATION.md)
   - plain-language Chinese explanation.
3. [docs/HOW_TO_RUN.md](docs/HOW_TO_RUN.md) - detailed run instructions.
4. [docs/STAT_METHODS_AS_BUILT.md](docs/STAT_METHODS_AS_BUILT.md) - statistical
   methods as implemented in code.
5. [docs/CODE_MAP.md](docs/CODE_MAP.md) - file-by-file code map.

The original design drafts are:

- [docs/Design1.tex](docs/Design1.tex): model layer, including immune response,
  toxicity, efficacy, and their probability structure.
- [docs/Design2.tex](docs/Design2.tex): decision layer, including admissible
  doses, utility, adaptive allocation, early stopping, and final OD selection.

Important: the current code implements the main adaptive dose-finding workflow,
but the design drafts and PDFs still contain unresolved choices. In particular,
confirm whether the final design should keep the current no-control-arm
implementation or add the control-arm / `gamma_j` allocation idea from
`Design2.tex`.

## Install R Packages

From R or RStudio:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "purrr",
  "isotone",
  "Iso",
  "gridExtra",
  "testthat"
))
```

## Verify the Repo

Set the working directory to the repository root, then run:

```r
testthat::test_dir("tests")
```

A recent full local run passed with:

```text
FAIL 0 | WARN 0 | SKIP 0 | PASS 399
```

For a quick smoke test:

```r
source("examples/simple_usage_example.R")
```

This runs one trial and writes example plots under `results/simple_example/`.
The `results/` directory is intentionally ignored by git.

## Core Workflow

The project workflow has three levels: one simulated trial, calibration, and
production operating characteristics.

### 1. Run One Simulated Trial

```r
setwd("/path/to/DoseFinding")

source("src/core/config.R")
source("src/core/main.R")

quiet_config <- trial_config
quiet_config$verbose_logging <- FALSE
quiet_config$log_early_termination <- FALSE

result <- run_trial_simulation(
  trial_config = quiet_config,
  p_YI = p_YI,
  p_YT_given_I = p_YT_given_I,
  p_YE_given_I = p_YE_given_I,
  rho0 = rho0,
  rho1 = rho1,
  seed = 123
)

result$final_od
result$terminated_early
isTRUE(result$poc_validated)
nrow(result$all_data)
```

Interpretation:

- `final_od`: selected optimal dose; `NA` means no final dose was selected.
- `terminated_early`: `TRUE` if the trial stopped because no dose remained
  acceptable.
- `poc_validated`: `TRUE` only when the final PoC evidence gate passed.
- `all_data`: simulated patient-level data.
- `all_alloc_probs`: allocation probabilities used at each stage.
- `posterior_summaries`: posterior estimates used for decision-making.

Use `isTRUE(result$poc_validated)` because early-stopped trials may not contain
a final PoC decision.

### 2. Understand the Trial Loop

One call to `run_trial_simulation()` does this:

1. Start with equal allocation across all doses.
2. Simulate one cohort of patient outcomes.
3. Update Bayesian posterior estimates for immune response, toxicity, and
   efficacy.
4. Apply monotonic smoothing: PAVA for immune response, BIVISO for conditional
   toxicity and efficacy.
5. Convert conditional toxicity and efficacy into marginal probabilities.
6. Build the admissible dose set using `phi_*` thresholds and `c_*` posterior
   credibility cutoffs.
7. Stop early if no dose is admissible.
8. Compute expected utility for each admissible dose.
9. Allocate the next cohort with higher probability to higher-utility doses.
10. At the final stage, select the best admissible dose only if the PoC gate
    passes.

### 3. Calibrate PoC

PoC calibration asks: under a null or flat scenario where no dose should clearly
win, how often does the design still claim a successful final dose?

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
calibration_results$achieved_rate
```

Important: current code examples use `target_rate = 0.10` for null/flat PoC
detection. The design drafts and source PDFs also discuss familywise Type I error
`0.05`. Confirm the final target before reporting production calibration
results.

### 4. Calibrate Early Termination

Early termination calibration asks: under unfavorable scenarios, how strict
should the admissibility thresholds be so bad trials stop at the desired rate?

```r
source("src/core/config.R")
source("src/optimization/early_termination_calibration.R")

early_results <- run_quick_early_termination_calibration(
  target_rate = 0.80,
  n_simulations = 100,
  verbose = FALSE
)

early_results$optimal_c_T
early_results$optimal_c_E
early_results$optimal_rate
```

Use larger simulation counts for final reporting. The quick version is for
checking that the workflow runs.

### 5. Run Production Operating Characteristics

After the design choices are confirmed, run larger simulations across scenario
families:

- null / flat scenarios
- unfavorable scenarios
- favorable signal scenarios
- boundary scenarios near the thresholds
- scenarios with different true optimal doses

Record at least:

- early termination rate
- PoC pass rate
- final OD selection frequencies
- correct selection rate when a true optimal dose is defined
- average sample size
- average allocation by dose
- mean final utility

## How the Program Judges a Dose

A dose must pass three layers.

First, it must be admissible:

```text
P(Tox < phi_T | data) > c_T
P(Eff > phi_E | data) > c_E
P(Imm > phi_I | data) > c_I
```

Second, among admissible doses, it receives an expected utility score from
`utility_table`. The score rewards efficacy and immune response, and penalizes
toxicity.

Third, at the end of a completed trial, the best utility dose must pass the PoC
gate. If the PoC probability is below `c_poc`, the code returns `NA` instead of
claiming a selected OD.

## Main Parameters

Simulation truth parameters:

- `p_YI`: true immune response probability by dose.
- `p_YT_given_I`: true toxicity probability by dose and immune status.
- `p_YE_given_I`: true efficacy probability by dose and immune status.
- `rho0`, `rho1`: toxicity-efficacy dependence parameters for `I = 0` and
  `I = 1`.

Clinical thresholds:

- `phi_T`: maximum acceptable toxicity.
- `phi_E`: minimum acceptable efficacy.
- `phi_I`: minimum acceptable immune response.

Posterior credibility cutoffs:

- `c_T`: required posterior confidence that toxicity is acceptable.
- `c_E`: required posterior confidence that efficacy is acceptable.
- `c_I`: required posterior confidence that immune response is acceptable.

Final evidence gate:

- `delta_poc`: pairwise comparison margin used by PoC.
- `c_poc`: required PoC probability for final selection.

Practical rule of thumb: `p_*` values define the simulated world; `phi_*` values
define clinical acceptability; `c_*` values define how much posterior confidence
is required before the trial acts.

## Examples and Notebooks

Runnable examples:

- [examples/simple_usage_example.R](examples/simple_usage_example.R): one
  complete simulation with basic plots.
- [examples/poc_calibration_demo.R](examples/poc_calibration_demo.R): small PoC
  calibration demo.
- [examples/comprehensive_calibration_demo.R](examples/comprehensive_calibration_demo.R):
  combined PoC and early termination calibration demo.

Interactive notebooks:

- [notebooks/simulation_notebook.qmd](notebooks/simulation_notebook.qmd):
  explore one or more simulations.
- [notebooks/poc_calibration_notebook.qmd](notebooks/poc_calibration_notebook.qmd):
  calibrate `c_poc`.
- [notebooks/threshold_calibration_notebook.qmd](notebooks/threshold_calibration_notebook.qmd):
  explore threshold calibration.
- [notebooks/design_walkthrough.qmd](notebooks/design_walkthrough.qmd):
  walkthrough of the design logic.

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
├── notebooks/
├── docs/
├── examples/
└── tests/
```

## Open Questions for the Project Lead

Before final production runs, confirm:

1. Should the current no-control-arm implementation be the official design, or
   should the control arm and `gamma_j` allocation idea from `Design2.tex` be
   implemented?
2. Should final PoC calibration target familywise Type I error `0.05`, or the
   current code convention of null/flat detection rate `0.10`?
3. What should final PoC compare: immune response, marginal efficacy, expected
   utility, or another endpoint?
4. Which parameters should be production-calibrated: only `c_poc`, or also
   `c_T`, `c_E`, and `c_I`?
5. Which operating characteristics are required in the final report?

## Troubleshooting

If R cannot find files, check that `getwd()` points to the repository root.

If examples print too much output, set:

```r
trial_config$verbose_logging <- FALSE
trial_config$log_early_termination <- FALSE
```

If a trial terminates early too often, the admissibility rules may be too strict
for the scenario. Review `phi_T`, `phi_E`, `phi_I`, `c_T`, `c_E`, and `c_I`.

If PoC passes too often in null/flat scenarios, increase `c_poc` or revisit the
PoC target definition.
