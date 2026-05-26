# Test notebook-first user workflow defaults.

library(testthat)

if (basename(getwd()) == "tests") {
  setwd("..")
}

source("src/core/config.R")

workflow_notebooks <- c(
  simulation = "notebooks/simulation_notebook.qmd",
  poc_calibration = "notebooks/poc_calibration_notebook.qmd",
  threshold_calibration = "notebooks/threshold_calibration_notebook.qmd"
)

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

evaluate_user_settings <- function(path) {
  env <- new.env(parent = baseenv())
  eval(parse(text = extract_qmd_chunk(path, "user_settings")), envir = env)
  as.list(env)
}

test_that("workflow notebooks keep backend calls outside the user settings chunk", {
  backend_call_patterns <- c(
    "source\\s*\\(",
    "run_trial_simulation\\s*\\(",
    "calibrate_c_poc\\s*\\(",
    "run_poc_parameter_search\\s*\\(",
    "run_quick_early_termination_calibration\\s*\\("
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

test_that("simulation notebook defaults match the calibrated backend config", {
  settings <- evaluate_user_settings(workflow_notebooks[["simulation"]])
  simulation_settings <- settings$simulation_settings

  expect_true(settings$quick_mode)
  expect_equal(simulation_settings$dose_levels, trial_config$dose_levels)
  expect_equal(simulation_settings$cohort_size, trial_config$cohort_size)
  expect_equal(simulation_settings$phi_T, trial_config$phi_T)
  expect_equal(simulation_settings$c_T, trial_config$c_T)
  expect_equal(simulation_settings$phi_E, trial_config$phi_E)
  expect_equal(simulation_settings$c_E, trial_config$c_E)
  expect_equal(simulation_settings$phi_I, trial_config$phi_I)
  expect_equal(simulation_settings$c_I, trial_config$c_I)
  expect_equal(simulation_settings$c_poc, trial_config$c_poc)
  expect_equal(simulation_settings$delta_poc, trial_config$delta_poc)
})

test_that("PoC calibration notebook defaults protect the current calibrated run", {
  settings <- evaluate_user_settings(workflow_notebooks[["poc_calibration"]])
  poc_settings <- settings$poc_settings

  expect_false(settings$quick_mode)
  expect_equal(poc_settings$c_T, trial_config$c_T)
  expect_equal(poc_settings$c_E, trial_config$c_E)
  expect_equal(poc_settings$c_I, trial_config$c_I)
  expect_equal(poc_settings$delta_poc, trial_config$delta_poc)
  expect_equal(poc_settings$target_rate, 0.10)
  expect_true(trial_config$c_poc %in% poc_settings$c_poc_candidates)
  expect_true(trial_config$c_poc %in% poc_settings$parameter_search_c_poc_candidates)
  expect_true(poc_settings$append_history_log)
  expect_true(poc_settings$run_parameter_search)
  expect_true(poc_settings$use_common_random_numbers)
  expect_equal(poc_settings$calibration_seed, 10000)
  expect_equal(poc_settings$parameter_search_seed, poc_settings$calibration_seed)
  expect_equal(poc_settings$n_simulations, 500)
  expect_equal(poc_settings$parameter_search_n_simulations, 500)
  expect_equal(
    poc_settings$parameter_search_grid,
    data.frame(
      c_T = c(0.55, 0.55, 0.55),
      c_E = c(0.45, 0.50, 0.55),
      c_I = c(0.70, 0.70, 0.70)
    )
  )
})

test_that("threshold calibration notebook includes calibrated defaults in its grids", {
  settings <- evaluate_user_settings(workflow_notebooks[["threshold_calibration"]])
  threshold_settings <- settings$threshold_settings

  expect_false(settings$quick_mode)
  expect_equal(threshold_settings$dose_levels, trial_config$dose_levels)
  expect_equal(threshold_settings$cohort_size, trial_config$cohort_size)
  expect_equal(threshold_settings$phi_T, trial_config$phi_T)
  expect_equal(threshold_settings$phi_E, trial_config$phi_E)
  expect_equal(threshold_settings$phi_I, trial_config$phi_I)
  expect_equal(threshold_settings$c_poc, trial_config$c_poc)
  expect_equal(threshold_settings$delta_poc, trial_config$delta_poc)
  expect_true(trial_config$c_T %in% threshold_settings$c_T_candidates)
  expect_true(trial_config$c_E %in% threshold_settings$c_E_candidates)
  expect_true(trial_config$c_I %in% threshold_settings$c_I_candidates)
  expect_equal(threshold_settings$n_sim_per_candidate, 500)
  expect_equal(threshold_settings$validation_n_sim, 500)
})
