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

test_that("workflow produces equal allocation in stage 1 and valid probabilities thereafter", {
  test_config <- trial_config
  test_config$n_stages <- 2
  test_config$cohort_size <- 6
  # Use relaxed thresholds to ensure trial continues
  test_config$c_T <- 0.5
  test_config$c_E <- 0.5
  test_config$c_I <- 0.5

  test_p_YI <- c(0.3, 0.5, 0.7)
  test_p_YT_given_I <- matrix(c(
    0.2, 0.3, 0.4,
    0.3, 0.4, 0.5
  ), ncol = 2, byrow = TRUE)
  test_p_YE_given_I <- matrix(c(
    0.4, 0.6, 0.8,
    0.6, 0.8, 0.9
  ), ncol = 2, byrow = TRUE)

  results <- run_trial_simulation(test_config, test_p_YI, test_p_YT_given_I, test_p_YE_given_I, rho0, rho1, seed = 77)

  expect_true(is.logical(results$terminated_early))
  expect_s3_class(results$all_alloc_probs, "data.frame")

  # If trial didn't terminate early, check allocation structure
  if (!results$terminated_early && nrow(results$all_alloc_probs) > 0) {
    # Stage 1 should be equal allocation across doses
    stage1_data <- results$all_alloc_probs[results$all_alloc_probs$Stage == 1, ]
    if (nrow(stage1_data) > 0) {
      stage1_probs <- stage1_data$Prob
      expect_length(stage1_probs, length(test_config$dose_levels))
      expect_true(all(abs(stage1_probs - 1 / length(test_config$dose_levels)) < 1e-8))
    }
    
    # Any recorded allocation rows should sum to 1 by stage
    summed <- aggregate(Prob ~ Stage, data = results$all_alloc_probs, sum)
    expect_true(all(abs(summed$Prob - 1) < 1e-6))
  } else if (results$terminated_early) {
    # Early termination: all_alloc_probs may be empty, which is acceptable
    expect_true(nrow(results$all_alloc_probs) == 0 || 
                all(c("Stage", "Dose", "Prob") %in% names(results$all_alloc_probs)))
  }
})
