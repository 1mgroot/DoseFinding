# Test notebook-first user workflow defaults.

library(testthat)

if (basename(getwd()) == "tests") {
  setwd("..")
}

source("src/core/config.R")

workflow_notebooks <- c(
  simulation = "notebooks/simulation_notebook.qmd",
  scenario_comparison = "notebooks/scenario_comparison_notebook.qmd",
  poc_calibration = "notebooks/poc_calibration_notebook.qmd",
  threshold_calibration = "notebooks/threshold_calibration_notebook.qmd"
)

test_that("deprecated handoff documents and stale output paths stay absent", {
  deprecated_paths <- c(
    ".cursorrules",
    "docs/CHATGPT_PROJECT_CONTEXT.md",
    "docs/PROJECT_READALOUD_EXPLANATION.md",
    "docs/PROJECT_READALOUD_EXPLANATION.html",
    "docs/PROJECT_READALOUD_MOBILE.html",
    "docs/PROJECT_READALOUD_EXPLANATION_files",
    "docs/STAT_METHODS_AS_BUILT.md",
    "notebooks/poc_calibration_notebook.pdf",
    "notebooks/simulation_notebook.pdf",
    "notebooks/results",
    "notebooks/threshold_calibration_notebook_cache"
  )

  expect_false(
    any(file.exists(deprecated_paths)),
    info = paste(
      "Deprecated documentation or stale output paths were recreated:",
      paste(deprecated_paths[file.exists(deprecated_paths)], collapse = ", ")
    )
  )

  notebook_text <- vapply(
    c(workflow_notebooks, design_walkthrough = "notebooks/design_walkthrough.qmd"),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1)
  )

  expect_false(any(grepl("(?m)^\\s*pdf\\s*:", notebook_text, perl = TRUE)))
  expect_false(any(grepl("notebooks/results", notebook_text, fixed = TRUE)))
})

extract_qmd_chunk <- function(path, label) {
  lines <- readLines(path, warn = FALSE)
  start <- grep(paste0("^```\\{r ", label, "\\}"), lines)
  if (length(start) != 1) {
    stop("Could not find exactly one chunk named ", label, " in ", path)
  }

  remaining <- lines[(start + 1):length(lines)]
  end <- grep("^```\\s*$", remaining)
  if (length(end) == 0) {
    stop("Could not find closing fence for chunk ", label, " in ", path)
  }

  paste(remaining[seq_len(end[[1]] - 1)], collapse = "\n")
}

evaluate_user_settings <- function(path, quick_mode_override = NULL) {
  chunk <- extract_qmd_chunk(path, "user_settings")
  if (!is.null(quick_mode_override)) {
    chunk <- sub(
      "^quick_mode <- (TRUE|FALSE)",
      paste0("quick_mode <- ", if (quick_mode_override) "TRUE" else "FALSE"),
      chunk
    )
  }

  env <- new.env(parent = baseenv())
  eval(parse(text = chunk), envir = env)
  as.list(env)
}

test_that("workflow notebooks keep backend calls outside the user settings chunk", {
  backend_call_patterns <- c(
    "source\\s*\\(",
    "run_trial_simulation\\s*\\(",
    "calibrate_c_poc\\s*\\(",
    "calibrate_separate_thresholds\\s*\\("
  )

  for (path in workflow_notebooks) {
    user_settings <- extract_qmd_chunk(path, "user_settings")
    for (pattern in backend_call_patterns) {
      expect_false(
        grepl(pattern, user_settings),
        info = paste("Unexpected backend call in", path, "matching", pattern)
      )
    }
  }
})

test_that("simulation notebook defaults use the current calibrated fallback values", {
  settings <- evaluate_user_settings(workflow_notebooks[["simulation"]])
  simulation_settings <- settings$simulation_settings

  expect_type(settings$quick_mode, "logical")
  expect_equal(simulation_settings$dose_levels, trial_config$dose_levels)
  expect_equal(simulation_settings$n_stages, if (settings$quick_mode) 3 else 5)
  expect_equal(simulation_settings$n_simulations, if (settings$quick_mode) 5 else 2000)
  expect_equal(simulation_settings$cohort_size, trial_config$cohort_size)
  expect_equal(simulation_settings$phi_T, trial_config$phi_T)
  expect_equal(simulation_settings$c_T, 0.35)
  expect_equal(simulation_settings$phi_E, trial_config$phi_E)
  expect_equal(simulation_settings$c_E, 0.55)
  expect_equal(simulation_settings$phi_I, trial_config$phi_I)
  expect_equal(simulation_settings$c_I, 0.50)
  expect_equal(simulation_settings$c_poc, 0.90)
  expect_equal(simulation_settings$delta_poc, trial_config$delta_poc)
  expect_true(simulation_settings$use_calibration_results)
  expect_equal(
    simulation_settings$threshold_calibration_results_path,
    "results/threshold_calibration/threshold_calibration_results.rds"
  )
  expect_equal(
    simulation_settings$poc_calibration_results_path,
    "results/notebook_calibration/poc_calibration_results.rds"
  )

  quick_settings <- evaluate_user_settings(
    workflow_notebooks[["simulation"]],
    quick_mode_override = TRUE
  )$simulation_settings

  expect_equal(quick_settings$n_stages, 3)
  expect_equal(quick_settings$n_simulations, 5)
})

test_that("simulation notebook can reuse saved threshold and PoC calibration results", {
  calibration_chunk <- extract_qmd_chunk(workflow_notebooks[["simulation"]], "calibrated_inputs")

  expect_true(grepl("recommended_thresholds", calibration_chunk, fixed = TRUE))
  expect_true(grepl("optimal_c_poc", calibration_chunk, fixed = TRUE))
  expect_true(grepl("Calibration values used in this simulation", calibration_chunk, fixed = TRUE))

  threshold_file <- tempfile(fileext = ".rds")
  poc_file <- tempfile(fileext = ".rds")
  saveRDS(
    list(recommended_thresholds = list(c_T = 0.61, c_I = 0.72, c_E = 0.53)),
    threshold_file
  )
  saveRDS(list(optimal_c_poc = 0.987), poc_file)

  env <- new.env(parent = baseenv())
  env$simulation_settings <- list(
    use_calibration_results = TRUE,
    threshold_calibration_results_path = threshold_file,
    poc_calibration_results_path = poc_file,
    c_T = 0.55,
    c_I = 0.70,
    c_E = 0.50,
    c_poc = 0.995
  )
  env$kable <- function(...) invisible(NULL)

  capture.output(eval(parse(text = calibration_chunk), envir = env))

  expect_equal(env$simulation_settings$c_T, 0.61)
  expect_equal(env$simulation_settings$c_I, 0.72)
  expect_equal(env$simulation_settings$c_E, 0.53)
  expect_equal(env$simulation_settings$c_poc, 0.987)
})

test_that("simulation notebook runs repeated trial simulations and aggregate plots", {
  simulation_chunk <- extract_qmd_chunk(workflow_notebooks[["simulation"]], "simulation")
  results_chunk <- extract_qmd_chunk(workflow_notebooks[["simulation"]], "results")

  expect_true(grepl("n_simulations <- simulation_settings$n_simulations", simulation_chunk, fixed = TRUE))
  expect_true(grepl("lapply(seq_len(n_simulations)", simulation_chunk, fixed = TRUE))
  expect_true(grepl("seed = simulation_settings$seed + sim_id - 1", simulation_chunk, fixed = TRUE))
  expect_false(grepl("simulation_results[[1]]", simulation_chunk, fixed = TRUE))
  expect_false(grepl("Example Trial", results_chunk, fixed = TRUE))
  expect_true(grepl("simulation_metrics <- data.frame", simulation_chunk, fixed = TRUE))
  expect_true(grepl("simulation_metrics$selected_overdose", simulation_chunk, fixed = TRUE))
  expect_true(grepl("Number of simulations:", results_chunk, fixed = TRUE))
  expect_true(grepl("Final selection summary across simulations", results_chunk, fixed = TRUE))
  expect_true(grepl("Probability of selecting each dose as ODB", results_chunk, fixed = TRUE))
  expect_true(grepl("Probability of Selecting Each Dose as ODB", results_chunk, fixed = TRUE))
  expect_true(grepl("Probability of selecting an overdose as ODB", results_chunk, fixed = TRUE))
  expect_true(grepl("Probability of Selecting an Overdose as ODB", results_chunk, fixed = TRUE))
  expect_true(grepl("results/simulation/simulation_metrics.csv", results_chunk, fixed = TRUE))
  expect_true(grepl("results/simulation/dose_selection_summary.csv", results_chunk, fixed = TRUE))
  expect_true(grepl("results/simulation/overdose_selection_summary.csv", results_chunk, fixed = TRUE))
  expect_true(grepl("allocation_probability_summary <- bind_rows", results_chunk, fixed = TRUE))
  expect_true(grepl("mean_prob = mean(Prob)", results_chunk, fixed = TRUE))
  expect_true(grepl("Mean Allocation Probabilities Across Simulations", results_chunk, fixed = TRUE))
  expect_true(grepl("results/simulation/allocation_probability_summary.csv", results_chunk, fixed = TRUE))
  expect_true(grepl("results/simulation/participant_allocation_summary.csv", results_chunk, fixed = TRUE))
  expect_true(grepl("results/simulation/final_participant_allocation_summary.csv", results_chunk, fixed = TRUE))
  expect_true(grepl("results/simulation/stage_enrollment_summary.csv", results_chunk, fixed = TRUE))
  expect_true(grepl("Unconditional Mean Participant Allocation by Dose Level and Stage", results_chunk, fixed = TRUE))
  expect_true(grepl("Final mean participant allocation by dose across simulations", results_chunk, fixed = TRUE))
  expect_true(grepl("Final Mean Participants by Dose", results_chunk, fixed = TRUE))
  expect_true(grepl("Mean Final Immune Response Posterior Across Simulations", results_chunk, fixed = TRUE))
  expect_true(grepl("conditional_mean_participants", results_chunk, fixed = TRUE))
  expect_true(grepl("trials reached this stage", results_chunk, fixed = TRUE))
})

test_that("simulation notebook keeps zero-allocation stages in cumulative plots", {
  results_chunk <- extract_qmd_chunk(workflow_notebooks[["simulation"]], "results")

  expect_true(grepl("participant_allocation_grid <- expand.grid", results_chunk, fixed = TRUE))
  expect_true(grepl("d = trial_config$dose_levels", results_chunk, fixed = TRUE))
  expect_true(grepl("stage = seq_len(trial_config$n_stages)", results_chunk, fixed = TRUE))
  expect_true(grepl("left_join(participant_counts", results_chunk, fixed = TRUE))
  expect_true(grepl("dplyr::coalesce(n_participants, 0L)", results_chunk, fixed = TRUE))
})

test_that("scenario comparison notebook supports multiple probability scenarios", {
  settings <- evaluate_user_settings(workflow_notebooks[["scenario_comparison"]])
  comparison_settings <- settings$comparison_settings
  scenarios <- settings$scenarios

  expect_type(settings$quick_mode, "logical")
  expect_equal(comparison_settings$dose_levels, trial_config$dose_levels)
  expect_equal(comparison_settings$n_stages, if (settings$quick_mode) 3 else 5)
  expect_equal(comparison_settings$n_simulations, if (settings$quick_mode) 5 else 2000)
  expect_equal(comparison_settings$cohort_size, trial_config$cohort_size)
  expect_equal(comparison_settings$seed, 11118)
  expect_true(comparison_settings$use_calibration_results)
  expect_gte(length(scenarios), 2)

  for (scenario in scenarios) {
    expect_true(all(c("name", "p_YI", "p_YT_given_I", "p_YE_given_I", "rho0", "rho1") %in% names(scenario)))
    expect_length(scenario$p_YI, length(trial_config$dose_levels))
    expect_equal(dim(scenario$p_YT_given_I), c(length(trial_config$dose_levels), 2))
    expect_equal(dim(scenario$p_YE_given_I), c(length(trial_config$dose_levels), 2))
  }

  quick_settings <- evaluate_user_settings(
    workflow_notebooks[["scenario_comparison"]],
    quick_mode_override = TRUE
  )$comparison_settings

  expect_equal(quick_settings$n_stages, 3)
  expect_equal(quick_settings$n_simulations, 5)
})

test_that("scenario comparison notebook saves final comparison tables", {
  scenario_check_chunk <- extract_qmd_chunk(workflow_notebooks[["scenario_comparison"]], "scenario_check")
  run_chunk <- extract_qmd_chunk(workflow_notebooks[["scenario_comparison"]], "run_scenarios")
  results_chunk <- extract_qmd_chunk(workflow_notebooks[["scenario_comparison"]], "results")

  expect_true(grepl("calculate_true_utilities", scenario_check_chunk, fixed = TRUE))
  expect_true(grepl("calculate_utility_from_true_probs", scenario_check_chunk, fixed = TRUE))
  expect_true(grepl("scenario_true_utility_table", scenario_check_chunk, fixed = TRUE))
  expect_true(grepl("true_optimal_dose", scenario_check_chunk, fixed = TRUE))
  expect_true(grepl("true_utility_Dose_", scenario_check_chunk, fixed = TRUE))
  expect_true(grepl("run_trial_simulation", run_chunk, fixed = TRUE))
  expect_true(grepl("scenario_seed_stride <- max(100000, n_simulations + 1)", run_chunk, fixed = TRUE))
  expect_true(grepl("scenario_metrics <- bind_rows", run_chunk, fixed = TRUE))
  expect_true(grepl("Scenario comparison summary with true optimal dose", results_chunk, fixed = TRUE))
  expect_true(grepl("scenario_summary_table", results_chunk, fixed = TRUE))
  expect_true(grepl("tidyr::complete", results_chunk, fixed = TRUE))
  expect_true(grepl("selection_label", results_chunk, fixed = TRUE))
  expect_true(grepl("tidyr::pivot_wider", results_chunk, fixed = TRUE))
  expect_true(grepl("true_optimal_selection_table", results_chunk, fixed = TRUE))
  expect_true(grepl("selection_rate_true_optimal", results_chunk, fixed = TRUE))
  expect_true(grepl("selection_rate_Dose_", results_chunk, fixed = TRUE))
  expect_true(grepl("selection_rate_No_OD", results_chunk, fixed = TRUE))
  expect_true(grepl("scenario_output_dir <- file.path(project_root_abs, \"results\", \"scenario_comparison\")", results_chunk, fixed = TRUE))
  expect_true(grepl("scenario_output_path <- function(filename) file.path(scenario_output_dir, filename)", results_chunk, fixed = TRUE))
  expect_true(grepl("scenario_output_path(\"scenario_truth_table.csv\")", results_chunk, fixed = TRUE))
  expect_true(grepl("scenario_output_path(\"scenario_true_utility_table.csv\")", results_chunk, fixed = TRUE))
  expect_true(grepl("scenario_output_path(\"scenario_metrics.csv\")", results_chunk, fixed = TRUE))
  expect_true(grepl("scenario_output_path(\"scenario_selection_table.csv\")", results_chunk, fixed = TRUE))
  expect_true(grepl("scenario_output_path(\"scenario_summary_table.csv\")", results_chunk, fixed = TRUE))
})

test_that("PoC calibration notebook separates quick smoke settings from production calibration", {
  settings <- evaluate_user_settings(workflow_notebooks[["poc_calibration"]])
  poc_settings <- settings$poc_settings

  expect_false(settings$quick_mode)
  expect_equal(poc_settings$c_T, 0.35)
  expect_equal(poc_settings$c_E, 0.55)
  expect_equal(poc_settings$c_I, 0.50)
  expect_equal(poc_settings$delta_poc, trial_config$delta_poc)
  expect_equal(poc_settings$target_rate, 0.10)
  expect_true(poc_settings$append_history_log)
  expect_true(poc_settings$use_common_random_numbers)
  expect_true(poc_settings$show_progress)
  expect_equal(poc_settings$progress_interval_seconds, 300)
  expect_true(poc_settings$use_threshold_calibration_results)
  expect_equal(
    poc_settings$threshold_calibration_results_path,
    "results/threshold_calibration/threshold_calibration_results.rds"
  )
  expect_equal(
    poc_settings$calibration_results_path,
    "results/notebook_calibration/poc_calibration_results.rds"
  )
  expect_equal(poc_settings$calibration_seed, 11118)
  expect_equal(poc_settings$c_poc_candidates, c(0.80, 0.90, 0.95, 0.98, 0.99, 0.995))
  expect_true(trial_config$c_poc %in% poc_settings$c_poc_candidates)
  expect_false(any(poc_settings$c_poc_candidates >= 1))
  expect_equal(poc_settings$n_simulations, 2000)

  quick_settings <- evaluate_user_settings(
    workflow_notebooks[["poc_calibration"]],
    quick_mode_override = TRUE
  )$poc_settings

  expect_equal(quick_settings$c_poc_candidates, c(0.8, 0.9, 0.95))
  expect_equal(quick_settings$n_simulations, 5)
  expect_true(quick_settings$show_progress)
  expect_equal(quick_settings$progress_interval_seconds, 60)
})

test_that("threshold calibration notebook settings include baseline cutoffs", {
  settings <- evaluate_user_settings(workflow_notebooks[["threshold_calibration"]])
  threshold_settings <- settings$threshold_settings

  expect_type(settings$quick_mode, "logical")
  expect_equal(threshold_settings$quick_mode, settings$quick_mode)
  expect_equal(threshold_settings$dose_levels, trial_config$dose_levels)
  expect_equal(threshold_settings$cohort_size, trial_config$cohort_size)
  expect_equal(threshold_settings$phi_T, trial_config$phi_T)
  expect_equal(threshold_settings$phi_E, trial_config$phi_E)
  expect_equal(threshold_settings$phi_I, trial_config$phi_I)
  expect_false("c_poc" %in% names(threshold_settings))
  expect_false("delta_poc" %in% names(threshold_settings))
  expect_equal(threshold_settings$rho0, trial_config$rho0)
  expect_equal(threshold_settings$rho1, trial_config$rho1)
  expect_equal(threshold_settings$c_T, 0.35)
  expect_equal(threshold_settings$c_I, 0.50)
  expect_equal(threshold_settings$c_E, 0.55)
  expect_false("c_T_start" %in% names(threshold_settings))
  expect_false("c_E_start" %in% names(threshold_settings))
  expect_false("c_I_start" %in% names(threshold_settings))
  expect_equal(threshold_settings$target_missing_range, c(0.80, 0.85))
  expect_equal(threshold_settings$high_tox_p_I, 0.95)
  expect_equal(threshold_settings$high_tox_marginal_p_T, c(0.30, 0.50))
  expect_equal(threshold_settings$high_tox_marginal_p_E, 0.95)
  expect_equal(threshold_settings$low_immune_p_I, c(0.10, 0.15))
  expect_equal(threshold_settings$low_immune_marginal_p_T, 0.01)
  expect_equal(threshold_settings$low_immune_marginal_p_E, 0.95)
  expect_equal(threshold_settings$low_eff_p_I, 0.95)
  expect_equal(threshold_settings$low_eff_marginal_p_T, 0.01)
  expect_equal(threshold_settings$low_eff_marginal_p_E, c(0.10, 0.20))
  expect_equal(threshold_settings$n_sim_per_candidate, if (settings$quick_mode) 5 else 1000)
  expect_true(threshold_settings$show_progress)
  expect_equal(threshold_settings$progress_interval_seconds, 300)
  expect_true(threshold_settings$append_history_log)
})

test_that("PoC calibration notebook can reuse saved threshold calibration results", {
  threshold_chunk <- extract_qmd_chunk(workflow_notebooks[["poc_calibration"]], "threshold_inputs")

  expect_true(grepl("recommended_thresholds", threshold_chunk, fixed = TRUE))
  expect_true(grepl("poc_settings[[threshold_name]]", threshold_chunk, fixed = TRUE))
  expect_true(grepl("Threshold values used for PoC calibration", threshold_chunk, fixed = TRUE))

  threshold_file <- tempfile(fileext = ".rds")
  saveRDS(
    list(recommended_thresholds = list(c_T = 0.61, c_I = 0.72, c_E = 0.53)),
    threshold_file
  )

  env <- new.env(parent = baseenv())
  env$poc_settings <- list(
    use_threshold_calibration_results = TRUE,
    threshold_calibration_results_path = threshold_file,
    c_T = 0.55,
    c_I = 0.70,
    c_E = 0.50
  )
  env$kable <- function(...) invisible(NULL)

  capture.output(eval(parse(text = threshold_chunk), envir = env))

  expect_equal(env$poc_settings$c_T, 0.61)
  expect_equal(env$poc_settings$c_I, 0.72)
  expect_equal(env$poc_settings$c_E, 0.53)
})

test_that("threshold calibration notebook displays parameter explanations", {
  settings <- evaluate_user_settings(workflow_notebooks[["threshold_calibration"]])
  guide_chunk <- extract_qmd_chunk(workflow_notebooks[["threshold_calibration"]], "parameter_guide")
  env <- list2env(settings, parent = baseenv())
  env$project_root_abs <- normalizePath(".", winslash = "/", mustWork = TRUE)

  guide_output <- capture.output(eval(parse(text = guide_chunk), envir = env))

  expect_true(grepl("User settings and parameter meanings", guide_chunk, fixed = TRUE))
  expect_true(grepl("Clinical thresholds", guide_chunk, fixed = TRUE))
  expect_true(grepl("phi_T", guide_chunk, fixed = TRUE))
  expect_true(grepl("phi_E", guide_chunk, fixed = TRUE))
  expect_true(grepl("phi_I", guide_chunk, fixed = TRUE))
  expect_true(grepl("Maximum acceptable marginal toxicity probability", guide_chunk, fixed = TRUE))
  expect_true(grepl("Minimum acceptable marginal efficacy probability", guide_chunk, fixed = TRUE))
  expect_true(grepl("Minimum acceptable immune response probability", guide_chunk, fixed = TRUE))
  expect_true(grepl("Baseline toxicity credibility cutoff", guide_chunk, fixed = TRUE))
  expect_true(grepl("inactive endpoint cutoffs held at fixed baseline values", guide_chunk, fixed = TRUE))
  expect_true(any(grepl("User settings and parameter meanings", guide_output, fixed = TRUE)))
  expect_equal(nrow(env$parameter_rows), 33)
  expect_false(any(c("c_T_start", "c_E_start", "c_I_start", "c_poc", "delta_poc") %in% env$parameter_rows$parameter))
})
