# PoC Calibration and Production Simulation Run Summary

> Historical production-run record. The numerical calibration and simulation
> results remain valid. The PDF/double-render and output-freshness notes near the
> end were superseded later on 2026-07-13 by commit `a3c8f5a`; current notebook
> source renders HTML only. See `f004_f006_f007_repair_summary.md` for the latest
> repair status.

Generated: 2026-07-13 01:11 EDT

Branch: `codex/fix-design-alignment-dr-f001-f004`

Latest audited implementation commit before this run:

- `89f1713 Update audit status after design-alignment fixes`

## Purpose

This file records the post-fix calibration and simulation evidence after applying the DR-F001 through DR-F004 design-alignment fixes:

- DR-F001: active simulator defaults use `rho0 = 0`, `rho1 = 0` for conditional independence of toxicity and efficacy given immune response and dose.
- DR-F002: utility decisions use posterior draw-average utility, `E[U(pi)]`.
- DR-F003: adaptive allocation uses posterior probability of being optimal among admissible doses.
- DR-F004: threshold calibration targets final admissible-set missing rate in the 80%-85% range.

## Inputs

Threshold calibration source:

- File: `results/threshold_calibration/threshold_calibration_results.rds`
- Selected `c_T = 0.50`
- Selected `c_I = 0.60`
- Selected `c_E = 0.75`

PoC calibration notebook settings:

- Source notebook: `notebooks/poc_calibration_notebook.qmd`
- `use_threshold_calibration_results = TRUE`
- Threshold RDS path: `results/threshold_calibration/threshold_calibration_results.rds`
- Candidate `C_PoC` grid: `0.80, 0.90, 0.95, 0.98, 0.99, 0.995`
- Null PoC target rate: `<= 0.10`
- Production simulation count: `2,000`

Simulation notebook settings:

- Source notebook: `notebooks/simulation_notebook.qmd`
- `quick_mode = FALSE`
- `n_simulations = 2,000`
- `n_stages = 5`
- `cohort_size = 15`
- `seed = 11118`
- `use_calibration_results = TRUE`
- Threshold RDS path: `results/threshold_calibration/threshold_calibration_results.rds`
- PoC RDS path: `results/notebook_calibration/poc_calibration_results.rds`
- `rho0 = 0`
- `rho1 = 0`

## Commands Run

```sh
quarto render notebooks/poc_calibration_notebook.qmd
quarto render notebooks/simulation_notebook.qmd
```

Both notebooks completed successfully.

## PoC Calibration Results

Output file:

- `results/notebook_calibration/poc_calibration_results.rds`

Selected value:

- `optimal_c_poc = 0.80`

Target and achieved rates:

- Target null PoC detection rate: `0.10`
- Achieved rate at selected `C_PoC`: `0.062`
- `control_achieved = TRUE`
- Number of null simulations per candidate: `2,000`

Candidate PoC detection rates:

| C_PoC | Null PoC detection rate |
|---:|---:|
| 0.800 | 0.0620 |
| 0.900 | 0.0240 |
| 0.950 | 0.0105 |
| 0.980 | 0.0025 |
| 0.990 | 0.0005 |
| 0.995 | 0.0000 |

Interpretation:

- The selected `C_PoC = 0.80` is the smallest candidate in the grid that controls the implemented null PoC detection rate below the 10% target.
- This run used the newly selected `c_T`, `c_I`, and `c_E` from threshold calibration.

## Production Simulation Results

Output files:

- `results/simulation/simulation_metrics.csv`
- `results/simulation/dose_selection_summary.csv`
- `results/simulation/overdose_selection_summary.csv`
- `results/simulation/allocation_probability_summary.csv`
- `results/simulation/participant_allocation_summary.csv`
- `results/simulation/final_participant_allocation_summary.csv`
- `results/simulation/stage_enrollment_summary.csv`
- `notebooks/simulation_notebook.html`
- `notebooks/simulation_notebook.pdf` (historical output; later removed when
  notebook rendering became HTML-only)

Overall Monte Carlo summary:

| Quantity | Value |
|---|---:|
| Number of simulations | 2,000 |
| Early termination rate | 0.5945 |
| PoC validation rate | 0.3125 |
| Mean participants | 40.89 |
| No OD rate | 0.6875 |
| Overdose selection rate | 0.0095 |

Dose selection summary:

| Dose | Selected n | Selection probability | True marginal toxicity | Overdose |
|---:|---:|---:|---:|:---|
| 1 | 0 | 0.0000 | 0.053 | FALSE |
| 2 | 4 | 0.0020 | 0.106 | FALSE |
| 3 | 343 | 0.1715 | 0.135 | FALSE |
| 4 | 259 | 0.1295 | 0.222 | FALSE |
| 5 | 19 | 0.0095 | 0.320 | TRUE |

Stage enrollment summary:

| Stage | Active trial rate | Unconditional mean participants | Conditional mean participants |
|---|---:|---:|---:|
| Stage 1 | 1.0000 | 15.0000 | 15 |
| Stage 2 | 0.4770 | 7.1550 | 15 |
| Stage 3 | 0.4260 | 6.3900 | 15 |
| Stage 4 | 0.4135 | 6.2025 | 15 |
| Stage 5 | 0.4095 | 6.1425 | 15 |

Allocation probability check:

- Stage 1 remained equal allocation across all five doses: `0.20` each.
- Stage 2 mean active-trial allocation probabilities were concentrated on doses 2 and 3:
  - Dose 1: `0.0257`
  - Dose 2: `0.6222`
  - Dose 3: `0.3508`
  - Dose 4: `0.0013`
  - Dose 5: `0.0000`

Interpretation:

- The simulation notebook successfully loaded the new calibrated threshold and PoC values rather than using the user-setting fallbacks.
- The final overdose selection probability was low in this scenario: `0.95%`.
- The no-OD rate remained high at `68.75%`, consistent with the conservative threshold and PoC filters in this configured scenario.

## Historical Runtime Warnings and Reproducibility Notes

The double-render and PDF-specific warnings below describe this recorded run,
not the current HTML-only notebook workflow.

Both PoC and simulation renders emitted this knitr warning:

```text
You changed the working directory to /Users/jz/Development/DoseFinding (probably via setwd()).
It will be restored to /Users/jz/Development/DoseFinding/notebooks.
```

The warning is caused by notebook setup chunks calling `setwd(project_root)`. It did not stop execution.

Because both notebooks declare multiple output formats, `quarto render` executed the computational chunks more than once: once for HTML and once for PDF. For this run, the simulation seeds are deterministic by simulation index, so repeated execution should produce the same Monte Carlo results. However, the double execution is still a reproducibility and runtime-efficiency risk for production workflows.

The PoC notebook also emitted PDF glyph warnings for symbols such as check marks and less-than-or-equal signs. These are report-rendering warnings, not statistical computation failures.

## Git and Artifact Notes

Tracked audit artifact created by this update:

- `_audit/decision_review/poc_and_simulation_run_summary.md`

The rendered notebooks and `results/` outputs are ignored by Git in the current repository configuration and were not staged as production source changes.
