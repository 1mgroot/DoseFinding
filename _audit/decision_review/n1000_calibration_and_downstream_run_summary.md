# N=1000 Threshold Calibration and Downstream Production Run

Generated: 2026-07-28 00:05 EDT

Branch: `codex/fix-design-alignment-dr-f001-f004`

Latest implementation commit before this run:

- `fa542cf Remove stale project documentation and outputs`

## Purpose

This record captures the production rerun requested after confirming that the
existing threshold calibration output used only 500 simulations per candidate.
The threshold production default was raised to 1,000 simulations per candidate,
then every dependent notebook was rerun in dependency order.

No threshold candidate reused another candidate's trial results. The threshold
run executed 31,000 full trial simulations: 31 candidates times 1,000
simulations per candidate.

## Commands and Dependency Order

```sh
quarto render threshold_calibration_notebook.qmd --to html
quarto render poc_calibration_notebook.qmd --to html
quarto render simulation_notebook.qmd --to html
quarto render scenario_comparison_notebook.qmd --to html
```

All commands were run from `notebooks/` and completed successfully. Each render
emitted only the known knitr warning that the setup chunk temporarily changes
the working directory to the project root.

## Threshold Calibration

Runtime:

- Started: 2026-07-27 21:35:53 EDT
- Finished: 2026-07-27 23:30:31 EDT
- Full trial simulations: 31,000
- Simulations per candidate: 1,000
- Target final admissible-set missing rate: 0.80-0.85

| Parameter | Selected value | Final missing rate |
|---|---:|---:|
| `c_T` | 0.45 | 0.802 |
| `c_I` | 0.60 | 0.810 |
| `c_E` | 0.75 | 0.821 |

All three selected candidates are within the target range. Relative to the
previous 500-simulation run, `c_T` changed from 0.50 to 0.45; `c_I` and `c_E`
were unchanged.

Primary output:

- `results/threshold_calibration/threshold_calibration_results.rds`

## PoC Calibration

Runtime:

- Started: 2026-07-27 23:30:56 EDT
- Finished: 2026-07-27 23:39:24 EDT
- Null simulations: 2,000
- Full trial simulations: 2,000
- Common random numbers: enabled

The notebook loaded `c_T = 0.45`, `c_I = 0.60`, and `c_E = 0.75` from the new
threshold RDS.

| `c_poc` | Null PoC detection rate |
|---:|---:|
| 0.800 | 0.0775 |
| 0.900 | 0.0295 |
| 0.950 | 0.0125 |
| 0.980 | 0.0025 |
| 0.990 | 0.0005 |
| 0.995 | 0.0000 |

The selected value remains `c_poc = 0.80`; its 0.0775 null detection rate
controls the prespecified 0.10 target.

Primary output:

- `results/notebook_calibration/poc_calibration_results.rds`

## Production Simulation

The simulation notebook loaded both new calibration RDS files and completed
2,000 trial simulations.

| Quantity | Value |
|---|---:|
| Early termination rate | 0.4990 |
| PoC validation rate | 0.3990 |
| Mean participants | 46.425 |
| No OD rate | 0.6010 |
| Overdose selection rate | 0.0180 |

| Dose | Selected n | Selection probability | True marginal toxicity | Overdose |
|---:|---:|---:|---:|:---|
| 1 | 0 | 0.0000 | 0.053 | FALSE |
| 2 | 0 | 0.0000 | 0.106 | FALSE |
| 3 | 403 | 0.2015 | 0.135 | FALSE |
| 4 | 359 | 0.1795 | 0.222 | FALSE |
| 5 | 36 | 0.0180 | 0.320 | TRUE |

Primary output:

- `results/simulation/simulation_metrics.csv`

## Scenario Comparison

The scenario notebook loaded both new calibration RDS files and completed
6,000 trial simulations: 2,000 for each of three scenarios.

| Scenario | Early stop | PoC validation | No OD | True-optimal selection |
|---|---:|---:|---:|---:|
| Flat null | 0.9405 | 0.0100 | 0.9900 | 0.0100 |
| Increasing immune and efficacy | 0.4990 | 0.3990 | 0.6010 | 0.1795 |
| Middle dose best | 0.6645 | 0.2535 | 0.7465 | 0.1925 |

Primary output:

- `results/scenario_comparison/scenario_metrics.csv`

## Artifact and Git Notes

The four rendered HTML notebooks and all files under `results/` are intentionally
ignored by Git. This tracked audit record preserves the run settings, dependency
order, and headline results. The source change raises the threshold production
default to 1,000 in both the notebook and backend default settings, with the
workflow regression test updated to prevent a return to 500.
