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

## Standard User Workflow

Users should interact with this project through the notebooks. You do not need
to edit files in `src/` or call backend functions directly for routine use.

1. Open `DoseFinding.Rproj` in RStudio.
2. Open one of the workflow notebooks in `notebooks/`.
3. Edit only the **User Settings** chunk near the top of the notebook.
4. Click **Run All** or **Render**.
5. Review the tables, plots, and generated files under `results/`.

Each workflow notebook starts in `quick_mode <- TRUE` for a fast smoke test. Set
`quick_mode <- FALSE` before running production-scale calibration or reporting
final results.

## Notebook Decision Guide

- [notebooks/simulation_notebook.qmd](notebooks/simulation_notebook.qmd):
  run one adaptive trial and inspect allocation, posterior summaries, and final
  OD selection.
- [notebooks/poc_calibration_notebook.qmd](notebooks/poc_calibration_notebook.qmd):
  calibrate `c_poc` under a null/flat scenario.
- [notebooks/threshold_calibration_notebook.qmd](notebooks/threshold_calibration_notebook.qmd):
  calibrate `c_T`, `c_E`, and `c_I` for early termination behavior.
- [notebooks/design_walkthrough.qmd](notebooks/design_walkthrough.qmd):
  explanatory walkthrough of how the design maps to the implementation; not
  required for routine runs.

Recommended order for a new analysis:

1. Run the simulation notebook in quick mode.
2. Run PoC calibration in quick mode, then production mode.
3. Run threshold calibration if early stopping needs tuning.
4. Return to the simulation notebook with calibrated values.

## Repository Guide

- Documentation index: [docs/README.md](docs/README.md)
- Run guide: [docs/HOW_TO_RUN.md](docs/HOW_TO_RUN.md)
- Statistical methods: [docs/STAT_METHODS_AS_BUILT.md](docs/STAT_METHODS_AS_BUILT.md)
- Code map: [docs/CODE_MAP.md](docs/CODE_MAP.md)
- Backend simulation code: [src/core/main.R](src/core/main.R)
- Backend calibration code: [src/optimization/poc_calibration.R](src/optimization/poc_calibration.R)

Generated notebooks, plots, calibration outputs, result files, and RStudio state
are ignored by git.

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

From the repository root:

```bash
Rscript -e 'testthat::test_dir("tests")'
```

For notebook smoke testing, open each notebook in RStudio with
`quick_mode <- TRUE` and click **Render**.

## What the Notebooks Do

The simulation notebook runs this workflow automatically:

1. Equal allocation across all doses in stage 1.
2. Simulate patient outcomes for immune response, toxicity, and efficacy.
3. Update Bayesian posterior estimates.
4. Apply monotonic smoothing with PAVA/BIVISO.
5. Build the admissible dose set using `phi_*` thresholds and `c_*` posterior
   credibility cutoffs.
6. Stop early if no dose is admissible.
7. Allocate later stages toward higher-utility admissible doses.
8. At the final stage, select the best admissible dose only if the PoC gate
   passes.

The calibration notebooks use the same backend engine but expose only the
settings users usually need: dose levels, stage count, cohort size, scenario
probabilities, thresholds, candidate grids, target rates, and simulation count.

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

## Implementation Scope

- The current R code implements a no-control-arm adaptive simulation workflow.
- `Design2.tex` includes a control-arm / `gamma_j` allocation concept, but that
  extension is not implemented in the current code.
- The notebooks use `target_rate = 0.10` for null/flat PoC detection calibration
  by default. Protocol-specific applications should set the target required by
  the study design.

## Project Structure

```text
DoseFinding/
├── src/
│   ├── core/
│   ├── decision/
│   ├── optimization/
│   └── utils/
├── notebooks/
├── docs/
├── examples/
└── tests/
```

## Advanced Developer Use

The functions in `src/` remain available for automated scripts, tests, and
development work. Routine users should prefer the notebooks; developers can use
[docs/HOW_TO_RUN.md](docs/HOW_TO_RUN.md) and the files in `tests/` as reference
for programmatic usage.
