# Code Map

This document maps the current DoseFinding implementation. It intentionally
avoids hard-coded line numbers so it remains useful as functions move.

## Standard Workflow

Routine work is notebook-first. The dependency order is:

1. `notebooks/threshold_calibration_notebook.qmd`
2. `notebooks/poc_calibration_notebook.qmd`
3. `notebooks/simulation_notebook.qmd`
4. `notebooks/scenario_comparison_notebook.qmd`

The PoC notebook reads the threshold calibration RDS by default. The simulation
and scenario-comparison notebooks read both calibration RDS files by default.
Fallback values in their User Settings chunks are used only when saved
calibration results are unavailable.

`notebooks/design_walkthrough.qmd` is an explanatory notebook. It uses
standalone demonstration defaults and is not a production calibration report.

## Source Code

### Core (`src/core/`)

- `main.R`
  - `run_trial_simulation()` orchestrates stage allocation, data generation,
    posterior updates, admissibility, early termination, adaptive allocation,
    and final PoC-gated selection.
  - Returns final admissible and PoC-eligible sets, pairwise PoC probabilities,
    candidate utilities, allocation history, and posterior summaries.
- `config.R`
  - Standalone backend defaults, utility table, scenario inputs, and trial
    settings.
  - Notebook production runs may replace its credibility cutoffs with saved
    calibration results.
- `simulate_data.R`
  - Generates patient-level immune, toxicity, and efficacy outcomes.
  - `rho0 = rho1 = 0` is the active design default and gives conditional T/E
    independence given dose and immune response.
- `model_utils.R`
  - Beta posterior sampling, PAVA, BIVISO, and conditional-to-marginal
    probability calculations.

### Decisions (`src/decision/`)

- `dose_decision.R`
  - Computes utility per posterior draw and averages it as `E[U(pi)]`.
  - Builds the admissible set using the toxicity, efficacy, and immune criteria
    jointly.
  - Allocates by posterior probability of being optimal among admissible doses,
    with equal credit for ties within a posterior draw.
  - Builds the Design2 immune-response PoC-eligible set and selects the
    highest-utility eligible dose.

### Calibration (`src/optimization/`)

- `threshold_calibration.R`
  - Separately calibrates `c_T`, `c_I`, and `c_E` under endpoint-specific
    unfavorable scenarios.
  - The production target is an 80%-85% final admissible-set missing rate.
  - Inactive endpoint cutoffs remain fixed at baseline values.
- `poc_calibration.R`
  - Calibrates `c_poc` under a null/flat scenario.
  - Supports common random numbers, progress reporting, summary-only production
    runs, and readable history logs.

### Utilities (`src/utils/`)

- `helpers.R`: posterior summaries and core plotting helpers.
- `plotting_extensions.R`: dose-response and multi-scenario plots.

## Notebooks

- `threshold_calibration_notebook.qmd`
  - Saves `results/threshold_calibration/threshold_calibration_results.rds`.
- `poc_calibration_notebook.qmd`
  - Reads the threshold RDS and saves
    `results/notebook_calibration/poc_calibration_results.rds`.
- `simulation_notebook.qmd`
  - Runs 5 simulations in quick mode or 2,000 in production mode.
  - Writes aggregate metrics under `results/simulation/`.
- `scenario_comparison_notebook.qmd`
  - Runs multiple truth scenarios and writes comparison tables under
    `results/scenario_comparison/`.
- `design_walkthrough.qmd`
  - Explains the relationship between Design1, Design2, and the implementation.

Notebook source supports HTML rendering only. Generated HTML, figures, caches,
RDS files, CSV results, and plots are intentionally ignored by Git.

## Tests

- `test_main.R`: integration behavior and returned traceability fields.
- `test_dose_decision.R`: draw-average utility, admissibility, posterior
  optimality allocation, and final selection.
- `test_threshold_calibration.R`: scenario construction, target metric, and
  candidate selection.
- `test_poc_calibration.R`: null calibration and reporting.
- `test_notebook_workflow.R`: notebook settings, dependency reuse, output
  paths, and repository hygiene.
- `test_sample_size_invariants.R`: enrollment invariants.
- `test_workflow_order.R`: stage order and early termination placement.
- Additional focused tests cover flat scenarios, Bayesian PoC, and early
  termination.

## Documentation and Audit Evidence

- `README.md`: current status, calibrated values, and project entry point.
- `docs/HOW_TO_RUN.md`: supported notebook workflow.
- `docs/Design1.tex`, `docs/Design2.tex`: original design drafts.
- `_audit/decision_review/f004_f006_f007_repair_summary.md`: latest repair and
  verification summary on the current branch.
- `_audit/decision_review/poc_and_simulation_run_summary.md`: dated production
  calibration and simulation record; later rendering notes are marked as
  superseded.

The original design-alignment audit is preserved on commit `b57e150` in the
`codex/check-code-design-alignment` branch rather than duplicated in the current
working tree.
