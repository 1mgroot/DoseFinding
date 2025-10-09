# Test Comprehensive Calibration Workflow
# This file tests the complete calibration framework integration

# Load required libraries
library(testthat)

# Source required functions
source("src/core/config.R")
source("src/optimization/poc_calibration.R")
source("src/optimization/early_termination_calibration.R")
source("src/utils/calibration_plots.R")

test_that("PoC calibration framework works correctly", {
  # Test PoC calibration with minimal simulations
  cat("Testing PoC calibration framework...\n")
  
  # Run quick calibration
  results <- run_quick_calibration(
    target_rate = 0.10,
    n_simulations = 10  # Very small for testing
  )
  
  # Check results structure
  expect_true(is.list(results))
  expect_true("calibration_results" %in% names(results))
  expect_true("optimal_c_poc" %in% names(results))
  expect_true("optimal_rate" %in% names(results))
  expect_true("target_rate" %in% names(results))
  
  # Check calibration results data frame
  expect_true(is.data.frame(results$calibration_results))
  expect_true(nrow(results$calibration_results) > 0)
  expect_true(all(c("c_poc", "poc_detection_rate", "poc_detection_rate_lower", 
                   "poc_detection_rate_upper", "n_simulations") %in% 
                 colnames(results$calibration_results)))
  
  # Check optimal values are reasonable
  expect_true(results$optimal_c_poc >= 0.5 && results$optimal_c_poc <= 1.0)
  expect_true(results$optimal_rate >= 0 && results$optimal_rate <= 1)
  expect_equal(results$target_rate, 0.10)
  
  cat("✓ PoC calibration framework test passed\n")
})

test_that("Early termination calibration framework works correctly", {
  # Test early termination calibration with minimal simulations
  cat("Testing early termination calibration framework...\n")
  
  # Run quick calibration
  results <- run_quick_early_termination_calibration(
    target_rate = 0.80,
    n_simulations = 10  # Very small for testing
  )
  
  # Check results structure
  expect_true(is.list(results))
  expect_true("calibration_results" %in% names(results))
  expect_true("optimal_threshold" %in% names(results))
  expect_true("optimal_rate" %in% names(results))
  expect_true("target_rate" %in% names(results))
  expect_true("threshold_type" %in% names(results))
  
  # Check calibration results data frame
  expect_true(is.data.frame(results$calibration_results))
  expect_true(nrow(results$calibration_results) > 0)
  expect_true(all(c("threshold", "termination_rate", "termination_rate_lower", 
                   "termination_rate_upper", "n_simulations") %in% 
                 colnames(results$calibration_results)))
  
  # Check optimal values are reasonable
  expect_true(results$optimal_threshold >= 0.7 && results$optimal_threshold <= 1.0)
  expect_true(results$optimal_rate >= 0 && results$optimal_rate <= 1)
  expect_equal(results$target_rate, 0.80)
  expect_equal(results$threshold_type, "c_T")
  
  cat("✓ Early termination calibration framework test passed\n")
})

test_that("Calibration visualization functions work correctly", {
  # Test calibration plotting functions
  cat("Testing calibration visualization functions...\n")
  
  # Create mock calibration results for testing
  mock_poc_results <- list(
    calibration_results = data.frame(
      c_poc = seq(0.7, 0.9, by = 0.05),
      poc_detection_rate = c(0.05, 0.08, 0.12, 0.15, 0.18),
      poc_detection_rate_lower = c(0.03, 0.06, 0.10, 0.13, 0.16),
      poc_detection_rate_upper = c(0.07, 0.10, 0.14, 0.17, 0.20),
      n_simulations = rep(100, 5)
    ),
    optimal_c_poc = 0.8,
    optimal_rate = 0.12,
    target_rate = 0.10
  )
  
  # Test PoC calibration curve plotting
  expect_silent({
    poc_plot <- plot_poc_calibration_curve(mock_poc_results, target_rate = 0.10)
  })
  expect_true(inherits(poc_plot, "ggplot"))
  
  # Test calibration summary plotting
  expect_silent({
    summary_plot <- plot_calibration_summary(mock_poc_results)
  })
  expect_true(inherits(summary_plot, "ggplot"))
  
  cat("✓ Calibration visualization functions test passed\n")
})

test_that("Calibration persistence functions work correctly", {
  # Test save/load functionality
  cat("Testing calibration persistence functions...\n")
  
  # Create mock results
  mock_results <- list(
    calibration_results = data.frame(
      c_poc = c(0.8, 0.9),
      poc_detection_rate = c(0.10, 0.15),
      poc_detection_rate_lower = c(0.08, 0.13),
      poc_detection_rate_upper = c(0.12, 0.17),
      n_simulations = c(100, 100)
    ),
    optimal_c_poc = 0.8,
    optimal_rate = 0.10,
    target_rate = 0.10
  )
  
  # Test save functionality
  test_file <- "test_calibration_results.RData"
  expect_silent(save_calibration_results(mock_results, test_file))
  expect_true(file.exists(test_file))
  
  # Test load functionality
  loaded_results <- load_calibration_results(test_file)
  expect_true(is.list(loaded_results))
  expect_equal(loaded_results$optimal_c_poc, mock_results$optimal_c_poc)
  expect_equal(loaded_results$target_rate, mock_results$target_rate)
  
  # Clean up
  file.remove(test_file)
  
  cat("✓ Calibration persistence functions test passed\n")
})

test_that("Unfavorable scenario generation works correctly", {
  # Test unfavorable scenario data generation
  cat("Testing unfavorable scenario generation...\n")
  
  # Test data generation
  expect_silent({
    data <- generate_unfavorable_scenario_data(
      config = flat_scenario_config,
      n_patients_per_dose = 5,
      seed = 123
    )
  })
  
  # Check data structure
  expect_true(is.data.frame(data))
  expect_true(nrow(data) > 0)
  expect_true(all(c("d", "Y_I", "Y_T", "Y_E") %in% colnames(data)))
  
  # Test probability matrix creation
  expect_silent({
    probs <- create_unfavorable_probability_matrices(n_doses = 3)
  })
  
  # Check probability matrices
  expect_true(is.list(probs))
  expect_true(all(c("p_YI", "p_YT_given_I", "p_YE_given_I") %in% names(probs)))
  expect_equal(length(probs$p_YI), 3)
  expect_equal(ncol(probs$p_YT_given_I), 2)
  expect_equal(ncol(probs$p_YE_given_I), 2)
  
  # Check that probabilities are unfavorable (high toxicity, low efficacy)
  expect_true(all(probs$p_YI <= 0.15))  # Low immune response
  expect_true(all(probs$p_YT_given_I >= 0.70))  # High toxicity
  expect_true(all(probs$p_YE_given_I <= 0.15))  # Low efficacy
  
  cat("✓ Unfavorable scenario generation test passed\n")
})

test_that("Calibration integration works with main simulation", {
  # Test that calibrated parameters work with main simulation
  cat("Testing calibration integration with main simulation...\n")
  
  # Create config with calibrated parameters (using mock values)
  calibrated_config <- flat_scenario_config
  calibrated_config$c_poc <- 0.85
  calibrated_config$c_T <- 0.90
  
  # Run a minimal simulation to test integration
  expect_silent({
    # This should not error out
    test_simulation <- run_calibration_simulation(
      config = calibrated_config,
      scenario_type = "flat_null",
      n_simulations = 1,
      seed = 123
    )
  })
  
  # Check that simulation returns a logical value
  expect_true(is.logical(test_simulation))
  
  cat("✓ Calibration integration test passed\n")
})

test_that("Calibration configuration is properly set up", {
  # Test that calibration configuration is properly defined
  cat("Testing calibration configuration...\n")
  
  # Check that calibration_config exists and has required elements
  expect_true(exists("calibration_config"))
  expect_true(is.list(calibration_config))
  
  # Check required configuration elements
  required_elements <- c("poc_target_rate", "poc_tolerance", "early_termination_target_rate", 
                        "early_termination_tolerance", "n_calibration_simulations", 
                        "n_validation_simulations", "c_poc_range", "c_T_range", 
                        "c_E_range", "c_I_range")
  
  for (element in required_elements) {
    expect_true(element %in% names(calibration_config), 
                info = paste("Missing configuration element:", element))
  }
  
  # Check that values are reasonable
  expect_true(calibration_config$poc_target_rate > 0 && calibration_config$poc_target_rate < 1)
  expect_true(calibration_config$early_termination_target_rate > 0 && 
              calibration_config$early_termination_target_rate < 1)
  expect_true(calibration_config$n_calibration_simulations > 0)
  expect_true(calibration_config$n_validation_simulations > 0)
  
  cat("✓ Calibration configuration test passed\n")
})

# Run all tests
cat("\n=== RUNNING COMPREHENSIVE CALIBRATION TESTS ===\n")
test_results <- test_dir("tests", filter = "test_comprehensive_calibration")
cat("=== COMPREHENSIVE CALIBRATION TESTS COMPLETED ===\n")
