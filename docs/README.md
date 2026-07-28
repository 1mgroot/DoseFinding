# DoseFinding Documentation

This directory contains the public project documentation for the DoseFinding
simulation framework.

## Recommended Reading Order

For a project-status review:

1. `../README.md`
2. `../_audit/decision_review/f004_f006_f007_repair_summary.md`
3. `../_audit/decision_review/poc_and_simulation_run_summary.md`
4. The rendered notebook HTML files in dependency order:
   threshold calibration, PoC calibration, simulation, then scenario comparison

For implementation work:

1. `CODE_MAP.md`
2. `HOW_TO_RUN.md`
3. `Design1.tex` and `Design2.tex`
4. `../notebooks/design_walkthrough.qmd`

## Document Map

- `HOW_TO_RUN.md`: how to run simulations, calibration, and
  notebooks.
- `CODE_MAP.md`: file-by-file code organization.
- `Design1.tex`: original model-layer design, focused on immune response,
  toxicity, efficacy, and their probability structure.
- `Design2.tex`: original decision-layer design, focused on utility,
  admissible dose sets, adaptive allocation, early stopping, and final OD
  selection.
- `../README.md`: current workflow, calibrated values, and project scope.
- `../_audit/decision_review/`: dated audit evidence and run records. Historical
  files in that directory may be superseded by later records and say so at the
  top when applicable.

Generated notebook HTML and numerical results are intentionally untracked. The
current output roots are `notebooks/*.html` and `results/`. Do not use the
obsolete `notebooks/results/` path or old notebook PDF renders.

## Source of Truth

If documentation and code disagree, treat the current code as the implementation
source of truth. `Design1.tex` and `Design2.tex` are design drafts and may
include ideas that are not implemented.
