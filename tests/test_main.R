library(testthat)

# Set working directory to project root for proper path resolution
if (basename(getwd()) == "tests") {
  setwd("..")
}

source("src/core/config.R")
source("src/core/main.R")

quiet_trial_config <- within(trial_config, {
  verbose_logging <- FALSE
  log_early_termination <- FALSE
})

test_that("trial stage seeds are reproducible without adjacent simulation overlap", {
  first_trial_stage_seeds <- generate_trial_stage_seeds(11118, 5)
  second_trial_stage_seeds <- generate_trial_stage_seeds(11119, 5)

  expect_equal(generate_trial_stage_seeds(11118, 5), first_trial_stage_seeds)
  expect_length(unique(first_trial_stage_seeds), 5)
  expect_length(unique(c(first_trial_stage_seeds, second_trial_stage_seeds)), 10)
  expect_null(formals(simulate_data_gumbel)$seed)
})

test_that("run_trial_simulation returns expected structure", {
  result <- run_trial_simulation(
    quiet_trial_config,
    p_YI,
    p_YT_given_I,
    p_YE_given_I,
    rho0,
    rho1,
    seed = 11118
  )
  
  expect_type(result, "list")
  expect_true(all(c("final_od", "all_data", "all_alloc_probs", "posterior_summaries", "terminated_early") %in% names(result)))
  expect_length(result$final_od, 1)
  expect_true(is.logical(result$terminated_early))
  expect_s3_class(result$all_data, "data.frame")
  expect_s3_class(result$all_alloc_probs, "data.frame")
})

test_that("allocation probabilities are well-formed when trial continues", {
  result <- run_trial_simulation(
    quiet_trial_config,
    p_YI,
    p_YT_given_I,
    p_YE_given_I,
    rho0,
    rho1,
    seed = 11118
  )
  
  if (!result$terminated_early) {
    summed <- aggregate(Prob ~ Stage, data = result$all_alloc_probs, sum)
    expect_true(all(abs(summed$Prob - 1) < 1e-6))
  } else {
    expect_true(result$terminated_early)
  }
})
