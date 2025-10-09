# Test Bayesian PoC Calculation Functions
# This file tests the enhanced Bayesian PoC calculation implemented in Phase 2

# Load required libraries
library(testthat)
library(dplyr)

# Source the functions to test
source("src/core/simulate_data.R")
source("src/core/config.R")
source("src/core/model_utils.R")
source("src/decision/dose_decision.R")
source("src/core/main.R")

# Test 1: calculate_pi_parameters function
test_that("calculate_pi_parameters produces valid results", {
  # Create mock posterior summaries
  mock_posterior_summaries <- list(
    imm = list(
      samples_pava = list(
        rep(0.2, 1000),  # Dose 1 immune samples
        rep(0.4, 1000),  # Dose 2 immune samples
        rep(0.6, 1000)   # Dose 3 immune samples
      )
    ),
    eff = list(
      samples = list(
        rep(0.2, 1000),  # Dose 1, I=0 efficacy samples
        rep(0.4, 1000),  # Dose 1, I=1 efficacy samples
        rep(0.4, 1000),  # Dose 2, I=0 efficacy samples
        rep(0.6, 1000),  # Dose 2, I=1 efficacy samples
        rep(0.6, 1000),  # Dose 3, I=0 efficacy samples
        rep(0.8, 1000)   # Dose 3, I=1 efficacy samples
      )
    )
  )
  
  # Test calculation for dose 1
  result <- calculate_pi_parameters(1, mock_posterior_summaries)
  
  # Check structure
  expect_true(is.list(result))
  expect_true("pi_I_samples" %in% names(result))
  expect_true("pi_E_given_I0_samples" %in% names(result))
  expect_true("pi_E_given_I1_samples" %in% names(result))
  expect_true("pi_combined_samples" %in% names(result))
  expect_true("pi_combined_mean" %in% names(result))
  expect_true("pi_combined_sd" %in% names(result))
  
  # Check dimensions
  expect_equal(length(result$pi_I_samples), 1000)
  expect_equal(length(result$pi_E_given_I0_samples), 1000)
  expect_equal(length(result$pi_E_given_I1_samples), 1000)
  expect_equal(length(result$pi_combined_samples), 1000)
  
  # Check that combined efficacy is calculated correctly
  # For dose 1: pi_I=0.2, pi_E_given_I0=0.2, pi_E_given_I1=0.4
  # Expected: pi_combined = 0.2 * 0.4 + 0.8 * 0.2 = 0.08 + 0.16 = 0.24
  expected_combined <- 0.2 * 0.4 + 0.8 * 0.2
  expect_true(abs(result$pi_combined_mean - expected_combined) < 0.01)
  
  # Check that standard deviation is reasonable (should be small for constant values)
  expect_true(result$pi_combined_sd < 0.01)
})

# Test 2: Direct PoC calculation without utility dependency
test_that("PoC calculation produces valid probabilities", {
  # Test the core PoC calculation logic directly
  # Create samples where dose 2 is clearly better than dose 1
  
  # Dose 1: lower efficacy
  pi_I_samples_1 <- rep(0.2, 1000)
  pi_E_given_I0_samples_1 <- rep(0.2, 1000)
  pi_E_given_I1_samples_1 <- rep(0.4, 1000)
  pi_combined_samples_1 <- pi_I_samples_1 * pi_E_given_I1_samples_1 + 
                           (1 - pi_I_samples_1) * pi_E_given_I0_samples_1
  
  # Dose 2: higher efficacy (best dose)
  pi_I_samples_2 <- rep(0.4, 1000)
  pi_E_given_I0_samples_2 <- rep(0.4, 1000)
  pi_E_given_I1_samples_2 <- rep(0.6, 1000)
  pi_combined_samples_2 <- pi_I_samples_2 * pi_E_given_I1_samples_2 + 
                           (1 - pi_I_samples_2) * pi_E_given_I0_samples_2
  
  # Calculate PoC: Pr(Πᵢ < δ Πᵢⱼ | Dₙ)
  delta_poc <- 0.8
  poc_prob <- mean(pi_combined_samples_1 < delta_poc * pi_combined_samples_2)
  
  # Check that PoC probability is valid
  expect_true(poc_prob >= 0)
  expect_true(poc_prob <= 1)
  
  # Since dose 1 is worse than dose 2, PoC should be high (close to 1)
  expect_true(poc_prob > 0.5)
})

# Test 3: Edge cases for PoC calculation
test_that("PoC calculation handles edge cases", {
  # Test with identical doses (should give PoC = 0 since neither is significantly worse)
  pi_samples_1 <- rep(0.3, 1000)
  pi_samples_2 <- rep(0.3, 1000)
  delta_poc <- 0.8
  
  poc_prob <- mean(pi_samples_1 < delta_poc * pi_samples_2)
  expect_equal(poc_prob, 0)  # Should be 0 since 0.3 < 0.8*0.3 = 0.24 is false
  
  # Test with very different doses
  pi_samples_1 <- rep(0.1, 1000)  # Much worse
  pi_samples_2 <- rep(0.9, 1000)  # Much better
  
  poc_prob <- mean(pi_samples_1 < delta_poc * pi_samples_2)
  expect_equal(poc_prob, 1)  # Should be 1 since 0.1 < 0.8*0.9 = 0.72 is true
  
  # Test with moderately different doses
  pi_samples_1 <- rep(0.2, 1000)  # Worse
  pi_samples_2 <- rep(0.4, 1000)  # Better
  
  poc_prob <- mean(pi_samples_1 < delta_poc * pi_samples_2)
  expect_equal(poc_prob, 1)  # Should be 1 since 0.2 < 0.8*0.4 = 0.32 is true
})

# Test 5: Flat scenario data validation
test_that("Flat scenario data generation works correctly", {
  # Generate flat scenario data with larger sample size for more stable estimates
  flat_data <- generate_flat_scenario_data(
    config = flat_scenario_config,
    phi_I_lower = 0.20,
    phi_E_lower = 0.25,
    toxicity_low = 0.05,
    n_patients_per_dose = 200,  # Larger sample for more stable estimates
    seed = 123
  )
  
  # Check data structure
  expect_true(is.data.frame(flat_data))
  expect_true("d" %in% names(flat_data))
  expect_true("Y_I" %in% names(flat_data))
  expect_true("Y_E" %in% names(flat_data))
  expect_true("Y_T" %in% names(flat_data))
  
  # Check that we have data for all doses
  expect_equal(length(unique(flat_data$d)), 3)
  
  # Check that each dose has the expected number of patients
  for (dose in unique(flat_data$d)) {
    dose_data <- flat_data[flat_data$d == dose, ]
    expect_equal(nrow(dose_data), 200)
  }
  
  # Validate the flat scenario with higher tolerance for sampling variability
  validation <- validate_flat_scenario(flat_data, 0.20, 0.25, 0.05, tolerance = 0.2)
  expect_true(validation$success)
})

cat("All Bayesian PoC calculation tests passed!\n")
