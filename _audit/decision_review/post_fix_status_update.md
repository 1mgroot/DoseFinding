# Post-Fix Status Update

> Historical checkpoint from before production PoC calibration and simulation.
> Its “Next Steps” section was completed later on 2026-07-13. For the resulting
> production run, see `poc_and_simulation_run_summary.md`; for the latest repair
> status, see `f004_f006_f007_repair_summary.md`.

Date: 2026-07-13

Branch: `codex/fix-design-alignment-dr-f001-f004`

## What Changed

The implementation branch now contains the four design-alignment fixes and follow-up documentation/audit commits:

- `6acf71b` - active `rho0/rho1` defaults set to `0`, matching conditional independence of toxicity and efficacy given immune response and dose.
- `ebe58f9` - decision utility changed from plug-in `U(E[pi])` to posterior draw-average `E[U(pi)]`.
- `538eede` - adaptive allocation changed from normalized utility to posterior probability of being optimal among admissible doses, with equal tie splitting per posterior draw.
- `ec5afd9` - threshold calibration target changed to final admissible-set missing rate of `80%-85%`.
- `197fcfd` - fixed stale threshold calibration documentation about inactive cutoffs.
- `49143cb` - added n=1 threshold calibration hand-check audit artifact.
- `d0e9489` - removed PoC parameters from the threshold calibration notebook user settings.
- `f6d4e49` - added production threshold calibration summary.

## Problem Solved

The code now better matches the clarified statistical design:

1. The active simulator no longer introduces non-design T/E dependence by default.
2. Utility-based decisions now use draw-level posterior uncertainty rather than plugging in posterior means.
3. Adaptive randomization now follows posterior probability of being optimal.
4. Threshold calibration now targets the intended `80%-85%` final admissible-set missing rate.
5. Threshold calibration now keeps inactive endpoint cutoffs fixed at baseline values and the documentation says so.
6. The threshold calibration notebook no longer exposes `c_poc` or `delta_poc`, which belong to the PoC calibration workflow.

## Production Threshold Calibration Result

Rendered notebook:

- `notebooks/threshold_calibration_notebook.qmd`

Ignored runtime outputs produced:

- `notebooks/threshold_calibration_notebook.html`
- `results/threshold_calibration/threshold_calibration_summary.csv`
- `results/threshold_calibration/threshold_calibration_results.rds`
- `results/threshold_calibration/threshold_calibration_history.md`

Production settings:

- `quick_mode = FALSE`
- `n_sim_per_candidate = 500`
- `rho0 = 0`
- `rho1 = 0`
- target final admissible-set missing rate: `80%-85%`

Selected threshold values:

| cutoff | selected value | final admissible-set missing rate | count | approx MC SE |
|---|---:|---:|---:|---:|
| `c_T` | `0.50` | `84.4%` | `422/500` | `1.62%` |
| `c_I` | `0.60` | `81.4%` | `407/500` | `1.74%` |
| `c_E` | `0.75` | `83.6%` | `418/500` | `1.66%` |

Interpretation:

- All three selected threshold cutoffs landed inside the intended `80%-85%` target interval.
- This is materially better than the n=1 and n=5 smoke checks, which were useful for arithmetic verification but not stable operating-characteristic estimation.
- This is still endpoint-specific threshold calibration only. Joint validation and `C_PoC` calibration remain separate steps.

## Verification Run

Full test suite after the threshold notebook cleanup and production threshold render:

```text
PASS 644 | FAIL 0 | WARN 0 | SKIP 0
```

## Historical Next Steps (Completed Later on 2026-07-13)

1. Push `codex/fix-design-alignment-dr-f001-f004` to GitHub.
2. Run production PoC calibration using the new threshold calibration RDS.
3. If PoC calibration succeeds, run the simulation notebook using the newly calibrated threshold and PoC inputs.
4. Add audit summaries for PoC calibration and simulation notebook outputs.
