# How to Run the DoseFinding Notebooks

The standard user workflow is notebook-only. Routine users should not edit
files in `src/` or call backend functions from the R console.

## Quick Start

1. Open `DoseFinding.Rproj` in RStudio.
2. Open a notebook from `notebooks/`.
3. Edit only the **User Settings** chunk near the top.
4. Leave `quick_mode <- TRUE` for a fast smoke test.
5. Click **Run All** or **Render**.
6. Review generated tables, plots, and files under `results/`.

Switch `quick_mode <- FALSE` before running production-scale calibration or
reporting final results.

## Notebook Decision Guide

### 1. Run One Trial

Use `notebooks/simulation_notebook.qmd`.

This notebook is for learning the design, checking one scenario, and inspecting:

- final OD selection
- early termination status
- PoC validation status
- posterior summaries
- allocation by dose and stage
- dose-response and allocation plots

Only edit the **User Settings** chunk. Common settings include dose levels,
stage count, cohort size, thresholds, posterior credibility cutoffs, PoC values,
scenario probabilities, and seed.

### 2. Calibrate PoC

Use `notebooks/poc_calibration_notebook.qmd`.

This notebook calibrates `c_poc` under a null/flat scenario. It generates:

- calibration curve
- candidate-level PoC detection rates
- early termination summaries
- detailed text report under `results/notebook_calibration/`

Quick mode uses a small number of simulations for smoke testing. Production mode
uses the larger candidate grid and simulation count defined in the **User
Settings** chunk.

### 3. Calibrate Early Termination Thresholds

Use `notebooks/threshold_calibration_notebook.qmd`.

This notebook calibrates:

- `c_T` for toxicity-driven early stopping
- `c_E` for efficacy-driven early stopping
- `c_I` for immune-response-driven early stopping

It reports recommended thresholds, early stop rates, validation summaries, and a
saved report under `results/threshold_calibration/`.

### 4. Understand the Design

Use `notebooks/design_walkthrough.qmd`.

This notebook is explanatory. It connects `Design1.tex` and `Design2.tex` to the
implemented simulation behavior. It is not required for routine analyses.

## Output Files

Generated outputs are intentionally ignored by git. Common locations:

```text
results/
├── plots/
├── notebook_calibration/
└── threshold_calibration/
```

Regenerate these files from the notebooks when needed.

## Common User Settings

Trial scale:

- `dose_levels`: dose labels used by the design.
- `n_stages`: number of trial stages.
- `cohort_size`: patients enrolled per stage.

Clinical thresholds:

- `phi_T`: maximum acceptable toxicity.
- `phi_E`: minimum acceptable efficacy.
- `phi_I`: minimum acceptable immune response.

Posterior credibility cutoffs:

- `c_T`: required confidence that toxicity is acceptable.
- `c_E`: required confidence that efficacy is acceptable.
- `c_I`: required confidence that immune response is acceptable.

PoC settings:

- `c_poc`: final evidence threshold.
- `delta_poc`: pairwise comparison margin.
- `target_rate`: null/flat PoC detection target.

Simulation truth:

- `p_YI`: true immune response probabilities by dose.
- `p_YT_given_I`: true toxicity probabilities by dose and immune status.
- `p_YE_given_I`: true efficacy probabilities by dose and immune status.
- `rho0`, `rho1`: toxicity-efficacy dependence parameters.

## Troubleshooting

If R cannot find files, open `DoseFinding.Rproj` first and rerun the notebook.
The notebooks also try to locate the project root automatically.

If calibration takes too long, keep `quick_mode <- TRUE` while checking setup.
Use production mode only when the notebook runs successfully in quick mode.

If a trial terminates early too often, review `phi_T`, `phi_E`, `phi_I`,
`c_T`, `c_E`, and `c_I` in the notebook's **User Settings** chunk.

If PoC passes too often in null/flat scenarios, increase `c_poc`, raise the
candidate grid, or revisit the protocol's PoC target definition.

## Advanced Developer Use

The backend functions in `src/` remain available for scripts, tests, and
development. For examples of direct function usage, inspect the files in
`examples/` and `tests/`. This is not the standard user path.

## Validation

Run the full test suite from the repository root:

```bash
Rscript -e 'testthat::test_dir("tests")'
```

Run a whitespace check before committing:

```bash
git diff --check
```
