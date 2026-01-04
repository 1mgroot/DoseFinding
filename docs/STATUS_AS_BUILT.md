# STATUS AS BUILT (Code-First, Today)

## 1) What this repo does today
Bayesian multi-stage dose-finding simulator with immune, toxicity, and efficacy endpoints. It simulates data via a Gumbel copula, enforces monotone dose-response with PAVA/BIVISO, screens doses via admissibility thresholds, adaptively randomizes on expected utility, and makes a final selection with a PoC check using posterior samples.  
Evidence: `src/core/main.R run_trial_simulation L15-L167`; `src/core/simulate_data.R simulate_data_gumbel L9-L61`; `src/core/model_utils.R apply_pava_on_samples L33-L48`, `apply_biviso_on_matrix L50-L99`; `src/decision/dose_decision.R get_admissible_set L101-L140`, `adaptive_randomization L142-L155`, `calculate_poc_probability L310-L343`.

**Key capabilities**
- Multi-stage simulation with stage-wise seeds and per-dose posterior updating (main.R L15-L167).
- Monotone smoothing via PAVA (univariate) and BIVISO (dose × immune) before decision-making (model_utils.R L33-L99).
- Admissible-set screening on toxicity/efficacy/immune posterior probabilities (dose_decision.R L101-L140).
- Utility-based adaptive allocation and final OD selection with PoC probability check (dose_decision.R L142-L155, L364-L425).
- **PoC calibration system** using null/flat scenarios to control Type I error rate at ~10% (src/optimization/poc_calibration_new.R L1-L591).
- **Parameter optimization framework** for systematic threshold and credibility parameter tuning (src/optimization/parameter_optimization.R L1-L383).
- Interactive notebooks for simulation and PoC calibration (notebooks/simulation_notebook.qmd, notebooks/poc_calibration_notebook.qmd).

**Key improvements from previous versions**
- PoC now uses posterior samples directly (not normal approximation) comparing best dose vs competitors (dose_decision.R L310-L343).
- PoC validation can gate final selection: if PoC threshold not met, no optimal dose is selected (dose_decision.R L393-L399).
- Verbose logging control added to reduce output during batch calibration runs (config parameter: verbose_logging).
- Comprehensive PoC calibration workflow with detailed reporting and diagnostic plots.

**Current limitations**
- No control-arm/γ_j logic is implemented despite being mentioned in some design docs (not present in decision or config code).
- Default config in config.R is 5-dose/5-stage; consistent with notebooks and calibration scripts.

## 2) Canonical run paths
- **Engine entrypoint (authoritative):** `run_trial_simulation(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1, seed=NULL)` (main.R L15-L167). Source from project root or setwd to root, then call the function with a config list.
- **Notebooks:**  
  - `notebooks/simulation_notebook.qmd`: sources the engine, sets a 5-dose, 5-stage config with relaxed thresholds (c_T=c_E=c_I=0.5), runs `run_trial_simulation`, and generates publication-ready plots (simulation_notebook.qmd L15-L369).  
  - `notebooks/poc_calibration_notebook.qmd`: sources engine + calibration, builds null/flat scenarios via `create_null_flat_scenario`, runs `calibrate_c_poc` to find optimal c_poc threshold for ~10% Type I error, generates calibration curves and detailed reports (poc_calibration_notebook.qmd L1-L260).  
- **Calibration system:**
  - `src/optimization/poc_calibration_new.R`: Core calibration functions including `create_null_flat_scenario` (L16-L48), `calibrate_c_poc` (L180-L322), `plot_calibration_curve` (L325-L363), and `generate_calibration_report` (L366-L590).
  - Outputs detailed reports with early termination analysis, posterior summaries, and recommendations.
- **Parameter optimization:**
  - `src/optimization/parameter_optimization.R`: Systematic parameter search framework for phi_T, phi_E, phi_I, c_T, c_E, c_I values.
  - `src/optimization/run_optimization.R`: Convenience wrapper with `quick_optimization()` and `comprehensive_optimization()` functions.
- **Authoritative path for results:** `run_trial_simulation` in `src/core/main.R` is the source of truth; notebooks and optimization scripts are wrappers that configure and analyze results.

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
- Default (`src/core/config.R`): 5 doses, 5 stages, cohort 15, φ_T=0.35 / c_T=0.5, φ_E=0.1 / c_E=0.5, φ_I=0.2 / c_I=0.5, c_poc=0.9, delta_poc=0.8, early termination enabled, utility table as defined (config.R L8-L64).  
- Notebooks:  
  - `simulation_notebook.qmd` uses 5 doses, 5 stages, cohort 15, relaxed cutoffs (c_T=c_E=c_I=0.5), φ_T=0.35, φ_E=0.1 (simulation_notebook.qmd L37-L98) — **matches default config**.  
  - `poc_calibration_notebook.qmd` uses 5 doses with stricter thresholds (c_T=c_E=0.3, c_I=0.2) for null scenario calibration to achieve ~10% Type I error (poc_calibration_notebook.qmd L62-L120).  
- Calibration scripts:
  - `poc_calibration_new.R` base_config uses 5 doses with calibration-specific thresholds (c_T=c_E=0.3, c_I=0.3) designed for null scenario testing (poc_calibration_new.R L194-L220).
  - Calibration targets ~10% PoC detection rate in null/flat scenarios to control Type I error.
- Parameter optimization:
  - `parameter_optimization.R` tests multiple parameter combinations across ranges: φ_T=[0.25-0.45], φ_E=[0.05-0.25], φ_I=[0.15-0.35], c_T/c_E/c_I=[0.6-0.95] (parameter_optimization.R L20-L38).
- Comparability: Default config and simulation notebook are now aligned at 5-dose with relaxed thresholds; calibration uses stricter thresholds appropriate for Type I error control.

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
- PoC formula/gating: Code now uses **posterior sample-based PoC** comparing best dose vs competitors with Pr(π_best > δ × π_competitor | D_n) (dose_decision.R L310-L343). When PoC threshold not met, final selection can return NA for optimal dose (dose_decision.R L393-L399). **Status: RESOLVED** — now aligned with Bayesian framework.
- PoC calibration: Comprehensive calibration system implemented to find optimal c_poc threshold achieving ~10% Type I error rate in null scenarios (poc_calibration_new.R). **Status: IMPLEMENTED**.
- Control arm / γ_j: Not implemented in code (no γ_j calc, no control params in config). **Status: MISSING/OUT-OF-SCOPE**.  
- Config alignment: Default config (5-dose, relaxed) now aligned with simulation notebook; calibration uses appropriately strict thresholds for Type I error control. **Status: RESOLVED**.

## 9) Calibration and Optimization Systems

### PoC Calibration System
- **Purpose**: Determine optimal c_poc threshold to achieve ~10% Type I error rate
- **Methodology**: 
  - Creates null/flat scenarios where all doses have identical response probabilities at admissibility thresholds
  - Runs 1000+ simulations per c_poc candidate value
  - Measures PoC detection rate (false positive rate)
  - Selects c_poc achieving closest to 10% detection rate
- **Key functions** (poc_calibration_new.R):
  - `create_null_flat_scenario()`: Constructs null scenario using total probability formula
  - `calibrate_c_poc()`: Main calibration workflow testing multiple c_poc values
  - `generate_calibration_report()`: Detailed report with early termination analysis
- **Outputs**: Calibration curves, optimal c_poc value, detailed text report with example cases
- **Evidence**: poc_calibration_notebook.qmd successfully calibrates c_poc across candidates

### Parameter Optimization System
- **Purpose**: Systematically search parameter space to maximize trial performance metrics
- **Parameters tested**: φ_T, φ_E, φ_I (thresholds), c_T, c_E, c_I (credibility), utility table variations
- **Metrics evaluated**: Completion rate, correct dose selection rate, mean utility, PoC success rate, allocation efficiency
- **Key functions** (parameter_optimization.R):
  - `generate_parameter_combinations()`: Creates parameter grid with intelligent sampling
  - `run_parameter_combination()`: Evaluates single parameter set across multiple simulations
  - `create_optimization_plots()`: Visualizes parameter sensitivity and performance
  - `find_best_parameters()`: Identifies optimal parameter set by specified metric
- **Convenience wrappers** (run_optimization.R): `quick_optimization()`, `comprehensive_optimization()`, `test_specific_params()`

## Next safe steps (incremental)
- **Calibration validation**: Test calibrated c_poc in signal scenarios to ensure adequate power (>50% true optimal selection)
- **Parameter optimization**: Run comprehensive optimization to identify best threshold/credibility combinations for target scenarios
- Add minimal γ_j/control-arm placeholders or clearly mark unsupported in docs/config if needed for future work
- Add unit tests for PoC calculation correctness (posterior sample-based comparisons)
- Profile PAVA/BIVISO performance in large-scale calibration runs; current n_sims=1000 performs well
- Document seed-handling strategy: stage_seed = base_seed + stage ensures reproducibility without RNG conflicts
- Consider adding sensitivity analysis tools to assess robustness of calibrated parameters across scenario variations
