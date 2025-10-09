# Test PoC Calibration Framework
# This file tests the PoC calibration functions implemented in Phase 3

# Load required libraries
library(testthat)
library(dplyr)

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
  c_poc_range <- c(0.8, 0.9, 0.95)
  n_simulations <- 10  # Very small for testing
  
  results <- calibrate_c_poc(
    target_rate = 0.10,
    flat_scenario_config = flat_scenario_config,
    n_simulations = n_simulations,
    c_poc_range = c_poc_range
  )
  
  # Check structure
  expect_true(is.list(results))
  expect_true("calibration_results" %in% names(results))
  expect_true("optimal_c_poc" %in% names(results))
  expect_true("optimal_rate" %in% names(results))
  expect_true("target_rate" %in% names(results))
  expect_true("n_simulations" %in% names(results))
  
  # Check calibration_results data frame
  expect_true(is.data.frame(results$calibration_results))
  expect_equal(nrow(results$calibration_results), length(c_poc_range))
  expect_true("c_poc" %in% names(results$calibration_results))
  expect_true("poc_detection_rate" %in% names(results$calibration_results))
  expect_true("poc_detection_rate_lower" %in% names(results$calibration_results))
  expect_true("poc_detection_rate_upper" %in% names(results$calibration_results))
  
  # Check that detection rates are valid (0-1)
  expect_true(all(results$calibration_results$poc_detection_rate >= 0))
  expect_true(all(results$calibration_results$poc_detection_rate <= 1))
  
  # Check that optimal_c_poc is in the range
  expect_true(results$optimal_c_poc %in% c_poc_range)
  
  # Check that target_rate is correct
  expect_equal(results$target_rate, 0.10)
  
  # Check that n_simulations is correct
  expect_equal(results$n_simulations, n_simulations)
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
  results <- run_quick_calibration(target_rate = 0.10, n_simulations = 5)
  
  # Check structure (same as calibrate_c_poc)
  expect_true(is.list(results))
  expect_true("calibration_results" %in% names(results))
  expect_true("optimal_c_poc" %in% names(results))
  expect_true("optimal_rate" %in% names(results))
  expect_true("target_rate" %in% names(results))
  expect_true("n_simulations" %in% names(results))
  
  # Check that results are valid
  expect_true(is.data.frame(results$calibration_results))
  expect_true(nrow(results$calibration_results) > 0)
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
    target_rate = 0.10,
    flat_scenario_config = flat_scenario_config,
    n_simulations = 1,
    c_poc_range = numeric(0)
  ))
  
  # Test with invalid target_rate
  expect_error(calibrate_c_poc(
    target_rate = -0.1,
    flat_scenario_config = flat_scenario_config,
    n_simulations = 1,
    c_poc_range = c(0.8)
  ))
  
  # Test with invalid n_simulations
  expect_error(calibrate_c_poc(
    target_rate = 0.10,
    flat_scenario_config = flat_scenario_config,
    n_simulations = 0,
    c_poc_range = c(0.8)
  ))
})

cat("All PoC calibration tests passed!\n")
