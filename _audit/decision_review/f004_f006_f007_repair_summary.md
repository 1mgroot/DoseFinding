# F-004, F-006, and F-007 Repair Summary

Generated: 2026-07-13

Branch: `codex/fix-design-alignment-dr-f001-f004`

## Findings Addressed

### F-004: Final Admissible-Set Traceability

Status after repair: addressed.

`run_trial_simulation()` now returns traceability fields for final selection:

- `final_admissible_set`
- `poc_eligible_set`
- `poc_pairwise_probs`
- `final_candidate_utilities`

Early-terminated trials also return explicit empty final candidate sets, plus
`poc_validated = FALSE` and `poc_probability = 0`.

### F-006: Stale HTML and Output Reproducibility

Status after repair: addressed for the rendered notebooks checked in this run.

Changes made:

- Removed default PDF rendering from heavy notebooks so `quarto render` does not
  execute simulation/calibration chunks twice.
- Re-rendered HTML for:
  - `notebooks/design_walkthrough.qmd`
  - `notebooks/poc_calibration_notebook.qmd`
  - `notebooks/simulation_notebook.qmd`
  - `notebooks/scenario_comparison_notebook.qmd`
- Fixed `scenario_comparison_notebook.qmd` so generated CSV outputs are written
  under root `results/scenario_comparison/` rather than notebook-relative
  `notebooks/results/scenario_comparison/`.

Verification:

- No rendered notebook or QMD search hit remained for old active values
  `rho0 = 1.5` or `rho1 = 2`.
- No QMD search hit remained for `pdf: default`.
- Root scenario comparison CSV output was refreshed:
  - `results/scenario_comparison/scenario_metrics.csv` had `6,000` rows.
  - Each of the three scenarios had `2,000` simulations.

Remaining note:

- Quarto still emits the known knitr `setwd()` warning. It did not stop any
  render in this run.

### F-007: Conditional vs Marginal Probability Wording

Status after repair: addressed for the reviewed user-facing and code comments.

Changes made:

- Updated simulator comments so `p_YT_given_I` and `p_YE_given_I` are described
  as conditional endpoint probabilities by dose and immune-response stratum.
- Updated README and HOW_TO_RUN wording to distinguish conditional inputs from
  marginal toxicity/efficacy values.
- Updated PoC calibration notebook wording so `null_p_E` and `null_p_T` are
  target marginal rates used to construct conditional scenario inputs.
- Updated design walkthrough comments to avoid calling conditional endpoint
  inputs marginal probabilities.

## Runtime Checks

Commands completed successfully:

```sh
Rscript -e 'testthat::test_file("tests/test_main.R"); testthat::test_file("tests/test_early_termination_poc.R"); testthat::test_file("tests/test_notebook_workflow.R")'
Rscript -e 'testthat::test_dir("tests")'
quarto render notebooks/design_walkthrough.qmd
quarto render notebooks/poc_calibration_notebook.qmd
quarto render notebooks/simulation_notebook.qmd
quarto render notebooks/scenario_comparison_notebook.qmd
```

Final full test result:

```text
PASS 658 | FAIL 0 | WARN 0 | SKIP 0
```

Key production render checks:

- PoC calibration: `c_poc = 0.80`, achieved null PoC detection rate `0.062`.
- Simulation notebook: `2,000` simulations; early termination rate `0.5945`,
  PoC validation rate `0.3125`, overdose selection rate `0.0095`.
- Scenario comparison: `6,000` total simulations, `2,000` per scenario.
