library(testthat)

# Set working directory to project root for proper path resolution
if (basename(getwd()) == "tests") {
  setwd("..")
}

source("src/core/config.R")
source("src/utils/helpers.R")
source("src/core/simulate_data.R")
source("src/core/model_utils.R")
source("src/decision/dose_decision.R")
source("src/core/main.R")

test_that("early termination triggers when admissible set is empty", {
  test_config <- trial_config
  test_config$dose_levels <- c(1, 2, 3)  # Use 3 doses for this test
  test_config$phi_T <- 0.1
  test_config$phi_E <- 0.8
  test_config$phi_I <- 0.8
  test_config$c_T <- 0.9
  test_config$c_E <- 0.9
  test_config$c_I <- 0.9

  test_p_YI <- c(0.1, 0.1, 0.1)
  test_p_YT_given_I <- matrix(c(
    0.3, 0.4, 0.5,
    0.4, 0.5, 0.6
  ), ncol = 2, byrow = TRUE)
  test_p_YE_given_I <- matrix(c(
    0.1, 0.1, 0.1,
    0.2, 0.2, 0.2
  ), ncol = 2, byrow = TRUE)

  result <- run_trial_simulation(test_config, test_p_YI, test_p_YT_given_I, test_p_YE_given_I, rho0, rho1, seed = 99)
  expect_true(result$terminated_early)
  expect_true(is.na(result$final_od))
  expect_true(is.na(result$termination_stage) | result$termination_stage >= 1)
})

test_that("PoC scenario returns structured outputs and valid allocation probabilities", {
  test_config <- trial_config
  test_config$dose_levels <- c(1, 2, 3)  # Use 3 doses for this test
  test_config$c_poc <- 0.8
  test_config$delta_poc <- 0.7

  test_p_YI <- c(0.2, 0.4, 0.6)
  test_p_YT_given_I <- matrix(c(
    0.1, 0.2, 0.3,
    0.2, 0.3, 0.4
  ), ncol = 2, byrow = TRUE)
  test_p_YE_given_I <- matrix(c(
    0.3, 0.5, 0.7,
    0.5, 0.7, 0.9
  ), ncol = 2, byrow = TRUE)

  result <- run_trial_simulation(test_config, test_p_YI, test_p_YT_given_I, test_p_YE_given_I, rho0, rho1, seed = 101)
  expect_true(is.logical(result$terminated_early))

  if (!result$terminated_early) {
    summed <- aggregate(Prob ~ Stage, data = result$all_alloc_probs, sum)
    expect_true(all(abs(summed$Prob - 1) < 1e-6))
    expect_true(length(result$final_od) == 1)
  }
})

test_that("baseline configuration yields logical outputs", {
  result <- run_trial_simulation(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1, seed = 202)
  expect_true(is.logical(result$terminated_early))
  if (!result$terminated_early) {
    expect_true(result$final_od %in% trial_config$dose_levels)
    summed <- aggregate(Prob ~ Stage, data = result$all_alloc_probs, sum)
    expect_true(all(abs(summed$Prob - 1) < 1e-6))
  }
})
