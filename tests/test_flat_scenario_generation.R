# Test Flat Scenario Generation Functions
# This file tests the flat scenario generation functions implemented in Phase 1

# Load required libraries
library(testthat)
library(dplyr)

# Source the functions to test
source("src/core/simulate_data.R")
source("src/core/config.R")

# Test 1: calculate_conditional_efficacy_flat function
test_that("calculate_conditional_efficacy_flat produces valid probabilities", {
  # Test with meeting requirements: phi_E = 0.25, phi_I = 0.20
  result <- calculate_conditional_efficacy_flat(0.25, 0.20)
  
  # Check that result is a matrix
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
  
  # Check that probabilities are valid (0-1)
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
  
  # Check that marginal efficacy is approximately correct
  p_E_given_I0 <- result[1, 1]
  p_E_given_I1 <- result[1, 2]
  marginal_efficacy <- (1 - 0.20) * p_E_given_I0 + 0.20 * p_E_given_I1
  expect_true(abs(marginal_efficacy - 0.25) < 0.01)
})

# Test 2: create_flat_probability_matrices function
test_that("create_flat_probability_matrices creates correct structure", {
  n_doses <- 3
  phi_I_lower <- 0.20
  phi_E_lower <- 0.25
  toxicity_low <- 0.05
  
  result <- create_flat_probability_matrices(n_doses, phi_I_lower, phi_E_lower, toxicity_low)
  
  # Check structure
  expect_true(is.list(result))
  expect_true("p_YI" %in% names(result))
  expect_true("p_YT_given_I" %in% names(result))
  expect_true("p_YE_given_I" %in% names(result))
  
  # Check dimensions
  expect_equal(length(result$p_YI), n_doses)
  expect_equal(nrow(result$p_YT_given_I), n_doses)
  expect_equal(ncol(result$p_YT_given_I), 2)
  expect_equal(nrow(result$p_YE_given_I), n_doses)
  expect_equal(ncol(result$p_YE_given_I), 2)
  
  # Check that all doses have identical probabilities
  expect_true(all(result$p_YI == phi_I_lower))
  expect_true(all(result$p_YT_given_I == toxicity_low))
})

# Test 3: generate_flat_scenario_data function
test_that("generate_flat_scenario_data produces flat scenario data", {
  # Use flat_scenario_config from config.R
  data <- generate_flat_scenario_data(
    config = flat_scenario_config,
    phi_I_lower = 0.20,
    phi_E_lower = 0.25,
    toxicity_low = 0.05,
    n_patients_per_dose = 20,
    seed = 123
  )
  
  # Check data structure
  expect_true(is.data.frame(data))
  expect_true("id" %in% names(data))
  expect_true("d" %in% names(data))
  expect_true("Y_I" %in% names(data))
  expect_true("Y_T" %in% names(data))
  expect_true("Y_E" %in% names(data))
  
  # Check that we have data for all doses
  expect_equal(length(unique(data$d)), length(flat_scenario_config$dose_levels))
  
  # Check that each dose has the expected number of patients
  for (dose in flat_scenario_config$dose_levels) {
    dose_data <- data[data$d == dose, ]
    expect_equal(nrow(dose_data), 20)
  }
  
  # Check metadata
  expect_equal(attr(data, "scenario_type"), "flat_null")
  expect_equal(attr(data, "phi_I_lower"), 0.20)
  expect_equal(attr(data, "phi_E_lower"), 0.25)
  expect_equal(attr(data, "toxicity_low"), 0.05)
})

# Test 4: validate_flat_scenario function
test_that("validate_flat_scenario correctly validates flat scenarios", {
  # Generate flat scenario data
  data <- generate_flat_scenario_data(
    config = flat_scenario_config,
    phi_I_lower = 0.20,
    phi_E_lower = 0.25,
    toxicity_low = 0.05,
    n_patients_per_dose = 100,  # Larger sample for more stable estimates
    seed = 123
  )
  
  # Validate the scenario
  validation <- validate_flat_scenario(data, 0.20, 0.25, 0.05, tolerance = 0.1)
  
  # Check validation structure
  expect_true(is.list(validation))
  expect_true("success" %in% names(validation))
  expect_true("immune_response_rates" %in% names(validation))
  expect_true("efficacy_rates" %in% names(validation))
  expect_true("toxicity_rates" %in% names(validation))
  expect_true("details" %in% names(validation))
  
  # Check that rates are calculated for all doses
  n_doses <- length(flat_scenario_config$dose_levels)
  expect_equal(length(validation$immune_response_rates), n_doses)
  expect_equal(length(validation$efficacy_rates), n_doses)
  expect_equal(length(validation$toxicity_rates), n_doses)
  
  # With large sample size and tolerance=0.1, validation should succeed
  expect_true(validation$success)
})

# Test 5: Integration test with main simulation
test_that("flat scenario integrates with main trial simulation", {
  # Generate flat scenario data
  data <- generate_flat_scenario_data(
    config = flat_scenario_config,
    phi_I_lower = 0.20,
    phi_E_lower = 0.25,
    toxicity_low = 0.05,
    n_patients_per_dose = 10,
    seed = 123
  )
  
  # Check that data can be used with existing simulation functions
  # This is a basic integration test - we're not running the full trial simulation
  # but checking that the data structure is compatible
  
  expect_true(nrow(data) > 0)
  expect_true(all(data$Y_I %in% c(0, 1)))
  expect_true(all(data$Y_T %in% c(0, 1)))
  expect_true(all(data$Y_E %in% c(0, 1)))
  expect_true(all(data$d %in% flat_scenario_config$dose_levels))
})

# Test 6: Edge cases
test_that("flat scenario functions handle edge cases", {
  # Test with extreme values
  result <- calculate_conditional_efficacy_flat(0.01, 0.01)
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
  
  result <- calculate_conditional_efficacy_flat(0.99, 0.99)
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
  
  # Test with single dose
  single_dose_config <- flat_scenario_config
  single_dose_config$dose_levels <- c(1)
  
  data <- generate_flat_scenario_data(
    config = single_dose_config,
    phi_I_lower = 0.20,
    phi_E_lower = 0.25,
    toxicity_low = 0.05,
    n_patients_per_dose = 5,
    seed = 123
  )
  
  expect_equal(length(unique(data$d)), 1)
  expect_equal(nrow(data), 5)
})

# Test 7: Reproducibility
test_that("flat scenario generation is reproducible", {
  # Generate data with same seed
  data1 <- generate_flat_scenario_data(
    config = flat_scenario_config,
    phi_I_lower = 0.20,
    phi_E_lower = 0.25,
    toxicity_low = 0.05,
    n_patients_per_dose = 10,
    seed = 123
  )
  
  data2 <- generate_flat_scenario_data(
    config = flat_scenario_config,
    phi_I_lower = 0.20,
    phi_E_lower = 0.25,
    toxicity_low = 0.05,
    n_patients_per_dose = 10,
    seed = 123
  )
  
  # Data should be identical
  expect_equal(data1, data2)
})

cat("All flat scenario generation tests passed!\n")
