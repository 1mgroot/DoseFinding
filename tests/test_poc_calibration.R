# Test PoC Calibration Framework
# This file tests the PoC calibration functions.

# Load required libraries
library(testthat)
library(dplyr)

if (basename(getwd()) == "tests") {
  setwd("..")
}

# Source the functions to test
source("src/core/simulate_data.R")
source("src/core/config.R")
source("src/core/model_utils.R")
source("src/decision/dose_decision.R")
source("src/core/main.R")
source("src/optimization/poc_calibration.R")

# Test 1: run_calibration_simulation function
test_that("run_calibration_simulation produces valid results", {
  # Test with flat null scenario
  config <- flat_scenario_config
  config$c_poc <- 0.9
  
  # Run a single simulation
  result <- run_calibration_simulation(config, "flat_null", 1, seed = 123)
  
  # Check that result is logical
  expect_true(is.logical(result))
  expect_true(result %in% c(TRUE, FALSE))
})

# Test 2: calibrate_c_poc function with small range
test_that("calibrate_c_poc produces valid calibration results", {
  # Test with a small range and few simulations for speed
  c_poc_candidates <- c(0.8, 0.9, 0.95)
  n_simulations <- 10  # Very small for testing
  null_scenario <- create_default_null_scenario(flat_scenario_config)
  
  results <- calibrate_c_poc(
    null_scenario = null_scenario,
    target_rate = 0.10,
    base_config = flat_scenario_config,
    n_simulations = n_simulations,
    c_poc_candidates = c_poc_candidates,
    verbose = FALSE
  )
  results_table <- poc_calibration_results_table(results)
  
  # Check structure
  expect_true(is.list(results))
  expect_true("calibration_results" %in% names(results))
  expect_true("optimal_c_poc" %in% names(results))
  expect_true("achieved_rate" %in% names(results))
  expect_true("target_rate" %in% names(results))
  expect_true("n_simulations" %in% names(results))
  
  # Check calibration results table helper
  expect_true(is.data.frame(results_table))
  expect_equal(nrow(results_table), length(c_poc_candidates))
  expect_true("c_poc" %in% names(results_table))
  expect_true("poc_detection_rate" %in% names(results_table))
  expect_true("poc_detection_rate_lower" %in% names(results_table))
  expect_true("poc_detection_rate_upper" %in% names(results_table))
  
  # Check that detection rates are valid (0-1)
  expect_true(all(results_table$poc_detection_rate >= 0))
  expect_true(all(results_table$poc_detection_rate <= 1))
  
  # Check that optimal_c_poc is in the range
  expect_true(results$optimal_c_poc %in% c_poc_candidates)
  
  # Check that target_rate is correct
  expect_equal(results$target_rate, 0.10)
  
  # Check that n_simulations is correct
  expect_equal(results$n_simulations, n_simulations)
})

test_that("calibrate_c_poc reuses simulation seeds across c_poc candidates", {
  original_runner <- get("run_single_calibration_simulation", envir = .GlobalEnv)
  on.exit(assign("run_single_calibration_simulation", original_runner, envir = .GlobalEnv), add = TRUE)

  calls <- data.frame(c_poc = numeric(), seed = numeric())
  assign(
    "run_single_calibration_simulation",
    function(config, scenario_params, seed = NULL) {
      calls <<- rbind(calls, data.frame(c_poc = config$c_poc, seed = seed))
      list(
        metrics = list(
          terminated_early = FALSE,
          termination_stage = NA_integer_,
          final_od = 1L,
          poc_validated = seed %% 2 == 0,
          poc_probability = 0.5
        ),
        debug_info = list(posterior_summaries = NULL),
        success = TRUE
      )
    },
    envir = .GlobalEnv
  )

  results <- calibrate_c_poc(
    null_scenario = create_default_null_scenario(flat_scenario_config),
    target_rate = 0.10,
    base_config = flat_scenario_config,
    n_simulations = 3,
    c_poc_candidates = c(0.90, 0.95),
    verbose = FALSE,
    store_simulation_results = FALSE,
    common_random_numbers = TRUE,
    calibration_seed = 500
  )

  expect_equal(calls$c_poc, c(rep(0.90, 3), rep(0.95, 3)))
  expect_equal(calls$seed, c(501, 502, 503, 501, 502, 503))
  expect_true(results$common_random_numbers)
  expect_equal(results$calibration_seed, 500)
  expect_true(all(vapply(
    results$calibration_results,
    function(candidate) isTRUE(candidate$common_random_numbers),
    logical(1)
  )))
})

# Test 3: validate_calibration function
test_that("validate_calibration produces valid validation results", {
  # Create mock calibration results
  mock_calibration_results <- list(
    optimal_c_poc = 0.85,
    target_rate = 0.10,
    optimal_rate = 0.12
  )
  
  # Run validation with very few simulations for testing
  validation_results <- validate_calibration(mock_calibration_results, n_validation_simulations = 5)
  
  # Check structure
  expect_true(is.list(validation_results))
  expect_true("optimal_c_poc" %in% names(validation_results))
  expect_true("validation_rate" %in% names(validation_results))
  expect_true("validation_ci" %in% names(validation_results))
  expect_true("target_rate" %in% names(validation_results))
  expect_true("n_validation_simulations" %in% names(validation_results))
  
  # Check that validation_rate is valid (0-1)
  expect_true(validation_results$validation_rate >= 0)
  expect_true(validation_results$validation_rate <= 1)
  
  # Check that confidence interval is valid
  expect_true(length(validation_results$validation_ci) == 2)
  expect_true(validation_results$validation_ci[1] <= validation_results$validation_ci[2])
  expect_true(validation_results$validation_ci[1] >= 0)
  expect_true(validation_results$validation_ci[2] <= 1)
  
  # Check that optimal_c_poc matches input
  expect_equal(validation_results$optimal_c_poc, mock_calibration_results$optimal_c_poc)
  
  # Check that target_rate matches input
  expect_equal(validation_results$target_rate, mock_calibration_results$target_rate)
})

# Test 4: run_quick_calibration function
test_that("run_quick_calibration produces valid results", {
  # Run quick calibration with very few simulations
  results <- run_quick_calibration(target_rate = 0.10, n_simulations = 5, verbose = FALSE)
  
  # Check structure (same as calibrate_c_poc)
  expect_true(is.list(results))
  expect_true("calibration_results" %in% names(results))
  expect_true("optimal_c_poc" %in% names(results))
  expect_true("achieved_rate" %in% names(results))
  expect_true("target_rate" %in% names(results))
  expect_true("n_simulations" %in% names(results))
  
  # Check that results are valid
  results_table <- poc_calibration_results_table(results)
  expect_true(is.data.frame(results_table))
  expect_true(nrow(results_table) > 0)
  expect_true(results$optimal_c_poc >= 0.7)
  expect_true(results$optimal_c_poc <= 0.95)
  expect_equal(results$target_rate, 0.10)
  expect_equal(results$n_simulations, 5)
})

# Test 5: save and load calibration results
test_that("save and load calibration results work correctly", {
  # Create mock calibration results
  mock_results <- list(
    calibration_results = data.frame(
      c_poc = c(0.8, 0.9),
      poc_detection_rate = c(0.05, 0.15),
      poc_detection_rate_lower = c(0.02, 0.08),
      poc_detection_rate_upper = c(0.12, 0.25),
      n_simulations = c(100, 100)
    ),
    optimal_c_poc = 0.85,
    optimal_rate = 0.10,
    target_rate = 0.10,
    n_simulations = 100
  )
  
  # Test file path
  test_file <- "results/test_calibration_results.RData"
  
  # Save results
  save_calibration_results(mock_results, test_file)
  
  # Check that file was created
  expect_true(file.exists(test_file))
  
  # Load results
  loaded_results <- load_calibration_results(test_file)
  
  # Check that loaded results match original
  expect_true(is.list(loaded_results))
  expect_equal(loaded_results$optimal_c_poc, mock_results$optimal_c_poc)
  expect_equal(loaded_results$optimal_rate, mock_results$optimal_rate)
  expect_equal(loaded_results$target_rate, mock_results$target_rate)
  expect_equal(loaded_results$n_simulations, mock_results$n_simulations)
  
  # Clean up test file
  if (file.exists(test_file)) {
    file.remove(test_file)
  }
})

# Test 6: Edge cases
test_that("calibration functions handle edge cases", {
  # Test with empty c_poc_range
  expect_error(calibrate_c_poc(
    null_scenario = create_default_null_scenario(flat_scenario_config),
    target_rate = 0.10,
    base_config = flat_scenario_config,
    n_simulations = 1,
    c_poc_candidates = numeric(0)
  ))
  
  # Test with invalid target_rate
  expect_error(calibrate_c_poc(
    null_scenario = create_default_null_scenario(flat_scenario_config),
    target_rate = -0.1,
    base_config = flat_scenario_config,
    n_simulations = 1,
    c_poc_candidates = c(0.8)
  ))
  
  # Test with invalid n_simulations
  expect_error(calibrate_c_poc(
    null_scenario = create_default_null_scenario(flat_scenario_config),
    target_rate = 0.10,
    base_config = flat_scenario_config,
    n_simulations = 0,
    c_poc_candidates = c(0.8)
  ))
})

test_that("PoC calibration history log appends readable run summaries", {
  test_file <- tempfile("poc-calibration-history-", fileext = ".md")
  on.exit(unlink(test_file), add = TRUE)

  mock_results <- list(
    calibration_results = data.frame(
      c_poc = c(0.90, 0.99),
      poc_detection_rate = c(0.22, 0.12),
      poc_detection_rate_lower = c(0.15, 0.07),
      poc_detection_rate_upper = c(0.30, 0.20),
      poc_se = c(0.04, 0.03),
      completion_rate = c(0.80, 0.75),
      poc_rate_among_completed = c(0.275, 0.16),
      n_simulations = c(100, 100)
    ),
    optimal_c_poc = 0.99,
    target_rate = 0.10,
    achieved_rate = 0.12,
    optimal_rate = 0.12,
    control_achieved = FALSE,
    c_poc_candidates = c(0.90, 0.99),
    n_simulations = 100,
    common_random_numbers = TRUE,
    calibration_seed = 123
  )

  null_scenario <- create_null_flat_scenario(n_doses = 2)
  base_config <- within(flat_scenario_config, {
    dose_levels <- c(1, 2)
  })

  append_poc_calibration_log(
    calibration_results = mock_results,
    null_scenario = null_scenario,
    base_config = base_config,
    file_path = test_file,
    run_label = "test run 1"
  )
  append_poc_calibration_log(
    calibration_results = mock_results,
    null_scenario = null_scenario,
    base_config = base_config,
    file_path = test_file,
    run_label = "test run 2"
  )

  log_lines <- readLines(test_file)
  expect_equal(sum(grepl("^## ", log_lines)), 2)
  expect_true(any(grepl("CONTROL NOT ACHIEVED", log_lines, fixed = TRUE)))
  expect_true(any(grepl("| c_poc | PoC detection |", log_lines, fixed = TRUE)))
  expect_true(any(grepl("test run 2", log_lines, fixed = TRUE)))
})

test_that("PoC parameter search validates grids and returns ranked summaries", {
  test_file <- tempfile("poc-parameter-search-", fileext = ".md")
  on.exit(unlink(test_file), add = TRUE)

  null_scenario <- create_default_null_scenario(flat_scenario_config)
  search_grid <- data.frame(
    c_T = c(0.20, 0.50),
    c_E = c(0.50, 0.50),
    c_I = c(0.35, 0.35)
  )

  search_results <- run_poc_parameter_search(
    null_scenario = null_scenario,
    base_config = flat_scenario_config,
    search_grid = search_grid,
    c_poc_candidates = c(0.95),
    n_simulations = 1,
    target_rate = 0.10,
    log_file = test_file,
    append_log = TRUE,
    verbose = FALSE
  )

  expect_true(file.exists(test_file))
  expect_true(is.list(search_results))
  expect_true(is.data.frame(search_results$summary))
  expect_equal(nrow(search_results$summary), nrow(search_grid))
  expect_equal(length(search_results$details), nrow(search_grid))
  expect_true(search_results$common_random_numbers)
  expect_equal(search_results$calibration_seed, 10000)
  expect_equal(search_results$details[[1]]$calibration_seed, 10000)
  expect_true(all(c(
    "search_id", "optimal_c_poc", "achieved_rate",
    "control_achieved", "completion_rate", "ranking_gap"
  ) %in% names(search_results$summary)))

  expect_error(
    run_poc_parameter_search(
      null_scenario = null_scenario,
      base_config = flat_scenario_config,
      search_grid = data.frame(cohort_size = 15),
      c_poc_candidates = c(0.95),
      n_simulations = 1,
      append_log = FALSE,
      verbose = FALSE
    ),
    "Invalid variables"
  )

  expect_error(
    run_poc_parameter_search(
      null_scenario = null_scenario,
      base_config = flat_scenario_config,
      search_grid = search_grid,
      c_poc_candidates = c(0.95),
      n_simulations = 1,
      append_log = FALSE,
      verbose = FALSE,
      common_random_numbers = "yes"
    ),
    "common_random_numbers"
  )
})
