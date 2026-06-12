# How to Run the DoseFinding Notebooks

The standard user workflow is notebook-only. Routine users should not edit
files in `src/` or call backend functions from the R console.

## Quick Start

1. Open `DoseFinding.Rproj` in RStudio.
2. Open a notebook from `notebooks/`.
3. Edit only the **User Settings** chunk near the top.
4. Set `quick_mode <- TRUE` for a fast smoke test, or `FALSE` for production.
5. Click **Run All** or **Render**.
6. Review generated tables, plots, and files under `results/`.

Switch `quick_mode <- FALSE` before running production-scale calibration or
reporting final results.

## Notebook Decision Guide

### 1. Run Trial Simulations

Use `notebooks/simulation_notebook.qmd`.

This notebook is for learning the design, checking one scenario, and inspecting
repeated trial simulation behavior. Production mode runs `2,000` independent
trial simulations by default. Quick mode runs `5` simulations for a fast smoke
test.

The notebook reports:

- final OD selection
- early termination status
- PoC validation status
- Monte Carlo selection rates
- posterior summaries
- allocation by dose and stage
- dose-response and allocation plots

Only edit the **User Settings** chunk. Common settings include dose levels,
stage count, cohort size, thresholds, posterior credibility cutoffs, PoC values,
scenario probabilities, and seed.

By default, the simulation notebook reads
`results/threshold_calibration/threshold_calibration_results.rds` and
`results/notebook_calibration/poc_calibration_results.rds` when they exist, then
uses those calibrated `c_T`, `c_I`, `c_E`, and `c_poc` values. If either file is
missing, the notebook falls back to the values in **User Settings**.

The Monte Carlo summary, mean posterior summaries, mean allocation probability
plot, mean participant-allocation plots, and
`results/simulation/simulation_metrics.csv` describe all simulation replicates.
Participant allocation plots are unconditional means across all planned
simulations, so stages after early termination count as `0` for that simulation.
The stage enrollment summary also reports how many trials reached each stage;
among trials that reach a stage, the conditional mean should match the cohort
size.

The allocation plots intentionally keep all doses in one graph. When several
doses have the same value, the notebook uses a small display-only horizontal
dodge so overlapping points are visible. For cumulative participant allocation,
the notebook fills dose-stage combinations with `0` participants before
calculating cumulative counts. This makes doses with no new patients in a stage
show as flat lines instead of disappearing or being connected across missing
stages.

### 2. Compare Multiple Scenarios

Use `notebooks/scenario_comparison_notebook.qmd`.

This notebook is for running the same design across multiple truth scenarios and
organizing the final operating characteristics into tables. Edit the `scenarios`
list in the **User Settings** chunk. Each scenario can define different:

- `p_YI`
- `p_YT_given_I`
- `p_YE_given_I`
- `rho0`
- `rho1`

The notebook writes:

- scenario truth table
- simulation-level metrics
- selection-rate table
- final scenario comparison summary

Common output location:

```text
results/scenario_comparison/
```

### 3. Calibrate Thresholds

Use `notebooks/threshold_calibration_notebook.qmd`.

This notebook calibrates `c_T`, `c_I`, and `c_E` separately before PoC
calibration. It generates:

- endpoint-specific unfavorable scenarios
- candidate-level final admissible set missing rates
- recommended `c_T`, `c_I`, and `c_E`
- readable calibration history under `results/threshold_calibration/`
- saved RDS and CSV summaries under `results/threshold_calibration/`

The default target is a final admissible set missing rate of 80%-90%. Set
`quick_mode <- TRUE` when you only want a fast smoke test.

### 4. Calibrate PoC

Use `notebooks/poc_calibration_notebook.qmd`.

This notebook calibrates `c_poc` under a null/flat scenario after `c_T`, `c_I`,
and `c_E` have been selected. It generates:

- calibration curve
- candidate-level PoC detection rates
- early termination summaries
- detailed text report under `results/notebook_calibration/`
- readable calibration history under `results/notebook_calibration/`

The notebook is prefilled for a focused validation run around the current
candidate region. It uses common random numbers so `c_poc` candidates are
compared on the same simulated trial streams.

By default, the notebook reads
`results/threshold_calibration/threshold_calibration_results.rds` and uses the
recommended `c_T`, `c_I`, and `c_E` from the threshold calibration notebook. If
that file is not available, it falls back to the values in the PoC notebook's
**User Settings** chunk.

If no tested `c_poc` controls the null PoC detection rate, first test a higher
or denser `c_poc` candidate grid. If that still fails, rerun the separate
threshold calibration workflow or revisit the protocol's PoC target definition;
do not tune `c_T`, `c_I`, or `c_E` inside the PoC notebook.

### 5. Understand the Design

Use `notebooks/design_walkthrough.qmd`.

This notebook is explanatory. It connects `Design1.tex` and `Design2.tex` to the
implemented simulation behavior. It is not required for routine analyses.

## Output Files

Generated outputs are intentionally ignored by git. Common locations:

```text
results/
├── plots/
├── simulation/
├── scenario_comparison/
├── notebook_calibration/
└── threshold_calibration/
```

Regenerate these files from the notebooks when needed.

## Common User Settings

Trial scale:

- `dose_levels`: dose labels used by the design.
- `n_stages`: number of trial stages.
- `n_simulations`: number of independent trial simulations. The simulation
  notebook defaults to `2,000` in production mode and `5` in quick mode.
- `cohort_size`: patients enrolled per stage.
- `scenarios`: scenario comparison list; each item contains one set of true
  probability inputs for `p_YI`, `p_YT_given_I`, and `p_YE_given_I`.

Clinical thresholds:

- `phi_T`: maximum acceptable toxicity.
- `phi_E`: minimum acceptable efficacy.
- `phi_I`: minimum acceptable immune response.

Posterior credibility cutoffs:

- `c_T`: required confidence that toxicity is acceptable.
- `c_E`: required confidence that efficacy is acceptable.
- `c_I`: required confidence that immune response is acceptable.
- `target_missing_range`: threshold calibration target for the final admissible
  set missing rate under endpoint-specific unfavorable scenarios.

PoC settings:

- `c_poc`: final evidence threshold.
- `delta_poc`: pairwise comparison margin.
- `target_rate`: null/flat PoC detection target.
- `append_history_log`: whether to append each run to a readable Markdown log.
- `history_log_path`: path to the keep-growing PoC calibration history file.
- `calibration_seed`: base seed for reproducible PoC calibration.
- `use_common_random_numbers`: compares `c_poc` candidates with the same
  simulation seeds so rankings are less noisy.
- `show_progress`: whether to print periodic progress and ETA while PoC
  calibration is running.
- `progress_interval_seconds`: approximate interval for extra progress messages.
- `use_threshold_calibration_results`: whether PoC calibration should read the
  saved threshold calibration RDS before calibrating `c_poc`.
- `threshold_calibration_results_path`: path to the saved threshold calibration
  RDS file.
- `calibration_results_path`: path where the PoC notebook saves its calibrated
  `c_poc` RDS result.

Simulation calibration reuse:

- `use_calibration_results`: whether the simulation notebook should read saved
  calibration results before running trial simulations.
- `threshold_calibration_results_path`: threshold calibration RDS used by the
  simulation notebook.
- `poc_calibration_results_path`: PoC calibration RDS used by the simulation
  notebook.

Current calibrated defaults:

- `c_T = 0.35`, `c_E = 0.60`, `c_I = 0.50`
- `c_poc = 0.90`, `delta_poc = 0.8`
- The focused PoC search is set up to target about `10%` null/flat PoC detection.

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
`c_T`, `c_E`, and `c_I` in the notebook's **User Settings** chunk, then rerun
the threshold calibration notebook before recalibrating `c_poc`.

If PoC passes too often in null/flat scenarios, increase `c_poc`, raise the
candidate grid, or revisit the protocol's PoC target definition.

If an allocation plot looks odd, first confirm whether the relevant dose simply
received `0` participants in that stage. Flat cumulative segments are expected
in that case. The plot is a display of the simulated allocation data; the small
horizontal offset is only for readability.

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

Render the simulation notebook from the command line:

```bash
cd notebooks
quarto render simulation_notebook.qmd --to html
```
