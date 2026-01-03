# STATUS AS BUILT (Code-First, Today)

## 1) What this repo does today
Bayesian multi-stage dose-finding simulator with immune, toxicity, and efficacy endpoints. It simulates data via a Gumbel copula, enforces monotone dose-response with PAVA/BIVISO, screens doses via admissibility thresholds, adaptively randomizes on expected utility, and makes a final selection with a PoC check (normal approximation).  
Evidence: `src/core/main.R run_trial_simulation L15-L167`; `src/core/simulate_data.R simulate_data_gumbel L9-L61`; `src/core/model_utils.R apply_pava_on_samples L33-L48`, `apply_biviso_on_matrix L50-L99`; `src/decision/dose_decision.R get_admissible_set L101-L140`, `adaptive_randomization L142-L155`, `calculate_poc_probability L310-L365`.

**Key capabilities**
- Multi-stage simulation with stage-wise seeds and per-dose posterior updating (main.R L15-L167).
- Monotone smoothing via PAVA (univariate) and BIVISO (dose × immune) before decision-making (model_utils.R L33-L99).
- Admissible-set screening on toxicity/efficacy/immune posterior probabilities (dose_decision.R L101-L140).
- Utility-based adaptive allocation and final OD selection with PoC probability check (dose_decision.R L142-L155, L396-L462).

**Key limitations**
- PoC uses a normal approximation on a combined efficacy metric with fixed sd=0.1 and still selects the best-utility dose even if PoC fails (dose_decision.R L310-L365, L396-L436) — potential Type I drift.
- No control-arm/γ_j logic is implemented despite being mentioned in docs (not present in decision or config code).
- Default config is 3-dose/3-stage; notebooks and calibration use different 5-dose, relaxed-threshold configs (config.R L8-L21 vs notebooks below).

## 2) Canonical run paths
- **Engine entrypoint (authoritative):** `run_trial_simulation(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1, seed=NULL)` (main.R L15-L167). Source from project root or setwd to root, then call the function with a config list.
- **Notebooks:**  
  - `notebooks/simulation_notebook.qmd`: sources the engine, sets a 5-dose, 5-stage, relaxed-threshold config, runs `run_trial_simulation`, and generates plots (simulation_notebook.qmd L37-L172).  
  - `notebooks/poc_calibration_notebook.qmd`: sources engine + calibration, builds null/flat scenarios, runs `calibrate_c_poc`, using 5 doses and relaxed thresholds (poc_calibration_notebook.qmd L62-L200).  
- **Authoritative path for results:** `run_trial_simulation` in `src/core/main.R` is the source of truth; notebooks are wrappers that change configs and add plotting.

## 3) Workflow contract (ordering + invariants)
- Order inside `run_trial_simulation` (main.R):  
  1) Stage loop over `n_stages` (L25-L52).  
  2) Simulate data for the stage via `simulate_data_gumbel` with stage-specific seed (L54-L66).  
  3) Posterior updates + PAVA/BIVISO for immune, tox, eff; compute marginals (L70-L98).  
  4) Admissible set via probability thresholds (L100-L109).  
  5) Early termination check; returns immediately if admissible set empty (L110-L126).  
  6) Adaptive randomization (stages < final) to set next-stage allocation probs (L132-L143).  
  7) Record allocation probs (L145).  
  8) Final selection with PoC validation after loop (L148-L167).  
- Invariants relied on:  
  - Stage-1 allocation uniform over doses (main.R L22-L52).  
  - When admissible set non-empty, adaptive allocation sums to 1 over admissible doses (dose_decision.R L146-L155).  
  - PoC probability computed for admissible doses; PoC threshold may fail but best dose still selected (dose_decision.R L310-L365, L396-L436).  
  - Returned list always includes `terminated_early` logical flag and data/alloc/posterior structures (main.R L115-L167).

## 4) Public API contracts
- `run_trial_simulation(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1, seed=NULL)`  
  - Requires `trial_config` with `dose_levels`, `n_stages`, `cohort_size`, `phi_T`, `phi_E`, `phi_I`, `c_T`, `c_E`, `c_I`, `c_poc`, `delta_poc`, `enable_early_termination`, `log_early_termination`, `utility_table` (config.R L8-L60).  
  - Returns list: `final_od`, `final_utility`, `poc_validated`, `poc_probability`, `selection_reason`, `all_data` (per-patient outcomes), `all_alloc_probs` (per-stage probabilities), `posterior_summaries`, `terminated_early`, `termination_stage`, `termination_reason` (main.R L115-L167).  
- Decision helpers:  
  - `get_admissible_set(posterior_summaries, config, verbose=TRUE)`: filters doses by posterior prob thresholds (dose_decision.R L101-L140).  
  - `adaptive_randomization(admissible_set, posterior_summaries, config)`: utility-proportional probs over admissible set (dose_decision.R L142-L155).  
  - `calculate_poc_probability(admissible_set, posterior_summaries, config)`: normal-approx PoC vs best utility dose (dose_decision.R L310-L365).  
  - `select_final_od_with_poc(...)`: picks best-utility dose, flags PoC, but still returns best dose if PoC fails (dose_decision.R L396-L436).  
- Sourcing pattern: main.R sources other modules with project-root relative paths (main.R L8-L13); tests setwd("..") before sourcing.

## 5) Configuration reality
- Default (`src/core/config.R`): 3 doses, 3 stages, cohort 6, φ_T=0.3 / c_T=0.9, φ_E=0.2 / c_E=0.9, φ_I=0.2 / c_I=0.8, c_poc=0.9, delta_poc=0.8, early termination enabled, utility table as defined (config.R L8-L60).  
- Notebooks:  
  - `simulation_notebook.qmd` uses 5 doses, 5 stages, cohort 15, relaxed cutoffs (c_T=c_E=c_I=0.5), different φ_T=0.35, φ_E=0.1 (simulation_notebook.qmd L37-L98).  
  - `poc_calibration_notebook.qmd` uses 5 doses with relaxed thresholds (c_T=c_E=0.3, c_I=0.2) and null scenario params (poc_calibration_notebook.qmd L62-L118).  
- Optimization/calibration scripts also assume 5-dose settings (e.g., poc_calibration_new.R base_config L179-L220).  
- Comparability: results differ between default 3-dose config and notebook/calibration 5-dose configs; align configs before comparing outcomes.

## 6) Reproducibility and randomness
- Stage seeds: base `seed` is incremented per stage (`stage_seed <- seed + stage`) (main.R L54-L66).  
- Data generation: `simulate_data_gumbel` calls `set.seed(seed)` if not NULL (simulate_data.R L19-L22); uses rbinom/rmultinom for outcomes (simulate_data.R L37-L52).  
- Within-stage RNG is reset each stage when a base seed is provided; no global RNG control beyond supplied seed.

## 7) Tests status (current)
- Run via: `R -q -e 'testthat::test_dir("tests")'`.  
- Coverage (examples):  
  - Structure/fields and allocation sums for `run_trial_simulation` (tests/test_main.R L11-L47).  
  - Expected utility, admissible set structure, adaptive allocation sums, final OD selection (tests/test_dose_decision.R L11-L65).  
  - Early termination and PoC scenario outputs, allocation sums when continuing (tests/test_early_termination_poc.R L15-L73).  
  - Workflow allocation checks with relaxed thresholds; tolerates early termination (tests/test_workflow_order.R L15-L52).  
- Gaps: no direct assertions on PAVA/BIVISO correctness, PoC approximation accuracy, or performance.

## 8) Known divergences vs design/docs
- PoC formula/gating: Code uses normal-approx PoC on combined efficacy with fixed sd=0.1 and still selects best dose if PoC fails (dose_decision.R L310-L365, L396-L436); TRIAL_DESIGN expects Pr(Π_i < δ Π_ij | D_n) gating. **Status: CONFLICT/DRIFT**.  
- Control arm / γ_j: Not implemented in code (no γ_j calc, no control params in config). **Status: MISSING**.  
- Config alignment: Defaults are 3-dose stringent; notebooks/calibration use 5-dose relaxed configs, so reported results diverge unless synced. **Status: DRIFT**.

## Next safe steps (incremental)
- Decide and document a single authoritative config (3-dose vs 5-dose) and propagate to notebooks/calibration/tests.  
- Clarify PoC contract (approx vs target formula) and, if keeping approximation, document its limits; optionally gate selection on PoC flag.  
- Add minimal γ_j/control-arm placeholders or clearly mark unsupported in docs/config.  
- Add unit tests for allocation sum invariants on non-empty admissible sets and for PoC probability monotonicity under known scenarios.  
- Profile PAVA/BIVISO cost in calibration/optimization runs; consider reducing n_sims or batching if slow.  
- Add a seed-handling note to docs/tests to avoid accidental identical draws across stages when reusing base seeds.  
- Ensure tests run from project root (setwd guard already in tests) and keep relative sourcing consistent.
