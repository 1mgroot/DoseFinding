# Test Bayesian PoC Calculation Functions
# This file tests the Bayesian PoC calculation used by final dose selection.

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

make_design2_poc_posterior <- function(
  immune_means = c(0.20, 0.40, 0.21),
  eff_means = c(0.20, 0.25, 0.30, 0.30, 0.90, 0.90),
  tox_means = c(0.10, 0.10, 0.10, 0.10, 0.01, 0.01),
  n_samples = 1000
) {
  list(
    imm = data.frame(
      pava_mean = immune_means,
      samples_pava = I(lapply(immune_means, rep, times = n_samples))
    ),
    eff = data.frame(
      pava_mean = eff_means,
      samples_pava = I(lapply(eff_means, rep, times = n_samples))
    ),
    tox = data.frame(
      pava_mean = tox_means,
      samples_pava = I(lapply(tox_means, rep, times = n_samples))
    )
  )
}

make_design2_poc_config <- function(c_poc = 0.9, delta_poc = 0.8) {
  config <- trial_config
  config$c_poc <- c_poc
  config$delta_poc <- delta_poc
  config$log_early_termination <- FALSE
  config$utility_table <- array(0, dim = c(2, 2, 2))
  config$utility_table[2, 1, 1] <- 100
  config$utility_table[2, 1, 2] <- 100
  config
}

test_that("Design2 PoC builds P_final from immune response against dose 1", {
  posterior_summaries <- make_design2_poc_posterior()
  config <- make_design2_poc_config(c_poc = 0.9, delta_poc = 0.8)

  result <- calculate_poc_probability(
    admissible_set = c(1, 2, 3),
    posterior_summaries = posterior_summaries,
    config = config
  )

  expect_equal(unname(result$pairwise_probs), c(0, 1, 0))
  expect_equal(unname(result$P_final), 2)
  expect_equal(result$poc_probability, 1)
  expect_true(check_poc_threshold(result, config))
  expect_equal(result$reference_dose, 1)
})

test_that("Design2 PoC does not automatically pass a single admissible dose", {
  posterior_summaries <- make_design2_poc_posterior()
  config <- make_design2_poc_config(c_poc = 0.9, delta_poc = 0.8)

  result <- calculate_poc_probability(
    admissible_set = 3,
    posterior_summaries = posterior_summaries,
    config = config
  )

  expect_equal(unname(result$pairwise_probs), 0)
  expect_equal(result$P_final, numeric(0))
  expect_equal(result$poc_probability, 0)
  expect_false(check_poc_threshold(result, config))
})

test_that("final OD is selected from P_final, not the whole admissible set", {
  posterior_summaries <- make_design2_poc_posterior()
  config <- make_design2_poc_config(c_poc = 0.9, delta_poc = 0.8)

  result <- select_final_od_with_poc(
    admissible_set = c(2, 3),
    posterior_summaries = posterior_summaries,
    config = config,
    verbose = FALSE
  )

  expect_true(result$poc_validated)
  expect_equal(unname(result$P_final), 2)
  expect_equal(result$optimal_dose, 2)
  expect_gt(result$utilities[2], result$utilities[1])
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
  
  # Check that we have data for all configured doses
  expect_equal(length(unique(flat_data$d)), length(flat_scenario_config$dose_levels))
  
  # Check that each dose has the expected number of patients
  for (dose in unique(flat_data$d)) {
    dose_data <- flat_data[flat_data$d == dose, ]
    expect_equal(nrow(dose_data), 200)
  }
  
  # Validate the flat scenario with higher tolerance for sampling variability
  validation <- validate_flat_scenario(flat_data, 0.20, 0.25, 0.05, tolerance = 0.2)
  expect_true(validation$success)
})
