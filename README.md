# Bayesian Adaptive Dose-Finding Trial Simulation

This repository contains an R simulation framework for Bayesian adaptive
dose-finding trials with three binary endpoints:

- immune response (`Y_I`)
- toxicity (`Y_T`)
- efficacy (`Y_E`)

The project simulates a multi-stage clinical trial, updates Bayesian posterior
estimates as patient outcomes accumulate, removes unacceptable doses, adaptively
allocates later patients toward better doses, stops early when no dose remains
acceptable, and applies a final Probability of Correct Selection (PoC) gate
before recommending an optimal dose (OD).

In plain terms: this project lets you repeatedly test a dose-finding design
before using it, so you can understand whether the design tends to pick a good
dose while controlling safety, efficacy, immune response, and false PoC claims.

## Current Status

The standard workflow is notebook-first. Routine users should not edit files in
`src/` or call backend R functions directly.

Current calibrated defaults:

- `c_T = 0.55`
- `c_E = 0.50`
- `c_I = 0.70`
- `c_poc = 0.995`
- `delta_poc = 0.8`

The PoC calibration notebook is set up for a focused null/flat validation run
targeting about `10%` null PoC detection, using common random numbers to reduce
Monte Carlo noise when comparing `c_poc` candidates.

## Quick Start

1. Open `DoseFinding.Rproj` in RStudio.
2. Open a notebook from `notebooks/`.
3. Edit only the **User Settings** chunk near the top.
4. Keep `quick_mode <- TRUE` for a fast smoke test.
5. Click **Run All** or **Render**.
6. Review the rendered notebook and generated files under `results/`.

Switch `quick_mode <- FALSE` only after the quick run works and you are ready
for production-scale simulation or calibration.

## Which Notebook Should I Use?

Use [notebooks/simulation_notebook.qmd](notebooks/simulation_notebook.qmd) to
run one adaptive trial and inspect allocation, posterior summaries, early
termination, PoC validation, and final OD selection.

Use [notebooks/poc_calibration_notebook.qmd](notebooks/poc_calibration_notebook.qmd)
to calibrate `c_poc` under a null/flat scenario. This notebook also appends a
readable calibration history under `results/notebook_calibration/` and includes
an optional parameter search over `c_T`, `c_E`, and `c_I`.

Use [notebooks/threshold_calibration_notebook.qmd](notebooks/threshold_calibration_notebook.qmd)
to tune `c_T`, `c_E`, and `c_I` for early termination behavior under unfavorable
scenarios.

Use [notebooks/design_walkthrough.qmd](notebooks/design_walkthrough.qmd) to
understand how the design documents map to the implementation. It is
explanatory and not required for routine runs.

Recommended order for a new analysis:

1. Run the simulation notebook in quick mode.
2. Run PoC calibration in quick mode.
3. Run PoC calibration in production mode.
4. Run threshold calibration if early stopping needs tuning.
5. Return to the simulation notebook with the calibrated values.

## Main Workflow

The simulation notebook runs this trial process automatically:

1. Allocate stage 1 patients equally across doses.
2. Simulate immune response, toxicity, and efficacy outcomes.
3. Update Bayesian posterior estimates.
4. Apply monotonic smoothing with PAVA/BIVISO where required by the design.
5. Build the admissible dose set using clinical thresholds (`phi_*`) and
   posterior credibility cutoffs (`c_*`).
6. Stop early if no dose remains admissible.
7. Allocate later stages toward higher-utility admissible doses.
8. At the final stage, recommend the best admissible dose only if the PoC gate
   passes.

The calibration notebooks use the same backend engine, but expose only the
settings users normally need: dose levels, stage count, cohort size, scenario
probabilities, thresholds, candidate grids, target rates, seeds, and simulation
counts.

## Key Parameters

Simulation truth parameters define the world being simulated:

- `p_YI`: true immune response probability by dose.
- `p_YT_given_I`: true toxicity probability by dose and immune status.
- `p_YE_given_I`: true efficacy probability by dose and immune status.
- `rho0`, `rho1`: toxicity-efficacy dependence parameters for `I = 0` and
  `I = 1`.

Clinical thresholds define what is acceptable:

- `phi_T`: maximum acceptable toxicity.
- `phi_E`: minimum acceptable efficacy.
- `phi_I`: minimum acceptable immune response.

Posterior credibility cutoffs define how much evidence is required:

- `c_T`: required posterior confidence that toxicity is acceptable.
- `c_E`: required posterior confidence that efficacy is acceptable.
- `c_I`: required posterior confidence that immune response is acceptable.

Final PoC settings control final selection:

- `delta_poc`: pairwise comparison margin used by PoC.
- `c_poc`: required PoC probability for final OD selection.

Rule of thumb: `p_*` values define the simulated world, `phi_*` values define
clinical acceptability, and `c_*` values define how much posterior confidence is
needed before the trial acts.

## Outputs

Generated notebooks, plots, calibration logs, reports, and result files are
ignored by git. Common output locations:

```text
results/
├── plots/
├── notebook_calibration/
└── threshold_calibration/
```

Regenerate these outputs from the notebooks when needed.

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

Run the full test suite from the repository root:

```bash
Rscript -e 'testthat::test_dir("tests")'
```

Run a whitespace check before committing:

```bash
git diff --check
```

For notebook smoke testing, open each notebook in RStudio with
`quick_mode <- TRUE` and click **Render**.

## Documentation Map

- [docs/HOW_TO_RUN.md](docs/HOW_TO_RUN.md): detailed notebook workflow and
  troubleshooting.
- [docs/CODE_MAP.md](docs/CODE_MAP.md): file-by-file code organization.
- [docs/Design1.tex](docs/Design1.tex): original model-layer design draft.
- [docs/Design2.tex](docs/Design2.tex): original decision-layer design draft.

If documentation and code disagree, treat the current code as the implementation
source of truth.

## Implementation Scope

- The current R code implements a no-control-arm adaptive simulation workflow.
- `Design2.tex` includes a control-arm / `gamma_j` allocation idea, but that
  extension is not implemented in the current code.
- The notebooks use `target_rate = 0.10` for null/flat PoC detection calibration
  by default. Protocol-specific applications should use the target required by
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

The functions in `src/` remain available for scripts, tests, and development
work. Routine users should prefer the notebooks; developers can use
[docs/HOW_TO_RUN.md](docs/HOW_TO_RUN.md) and the files in `tests/` as reference
for programmatic usage.
