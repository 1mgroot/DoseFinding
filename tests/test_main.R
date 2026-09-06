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

test_that("rho zero matches conditional independence for T and E", {
  p_t <- 0.25
  p_e <- 0.50
  cells <- as.numeric(Gumbel(p_t, p_e, 0))
  independent_cells <- c(
    (1 - p_t) * (1 - p_e),
    (1 - p_t) * p_e,
    p_t * (1 - p_e),
    p_t * p_e
  )

  expect_equal(cells, independent_cells)
  expect_equal(rho0, 0)
  expect_equal(rho1, 0)
  expect_equal(formals(simulate_data_gumbel)$rho0, 0)
  expect_equal(formals(simulate_data_gumbel)$rho1, 0)
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
  expected_names <- c(
    "final_od",
    "final_od_index",
    "final_utility",
    "poc_validated",
    "poc_probability",
    "selection_reason",
    "final_admissible_set",
    "final_admissible_indices",
    "poc_eligible_set",
    "poc_eligible_indices",
    "poc_pairwise_probs",
    "final_candidate_utilities",
    "all_data",
    "all_alloc_probs",
    "posterior_summaries",
    "terminated_early"
  )
  expect_true(all(expected_names %in% names(result)))
  expect_length(result$final_od, 1)
  expect_true(is.logical(result$terminated_early))
  expect_true(is.logical(result$poc_validated))
  expect_s3_class(result$all_data, "data.frame")
  expect_s3_class(result$all_alloc_probs, "data.frame")
  expect_true(all(result$final_admissible_set %in% quiet_trial_config$dose_levels))
  expect_equal(
    result$final_admissible_set,
    quiet_trial_config$dose_levels[result$final_admissible_indices]
  )
  expect_true(all(result$poc_eligible_set %in% result$final_admissible_set))
  expect_equal(
    result$poc_eligible_set,
    quiet_trial_config$dose_levels[result$poc_eligible_indices]
  )
  expect_length(result$poc_pairwise_probs, length(result$final_admissible_set))
  expect_length(result$final_candidate_utilities, length(result$final_admissible_set))
})

test_that("early termination returns traceable empty final candidate sets", {
  test_config <- within(quiet_trial_config, {
    dose_levels <- c(1, 2, 3)
    phi_T <- 0.1
    phi_E <- 0.8
    phi_I <- 0.8
    c_T <- 0.9
    c_E <- 0.9
    c_I <- 0.9
  })

  test_p_YI <- c(0.1, 0.1, 0.1)
  test_p_YT_given_I <- matrix(c(
    0.3, 0.4, 0.5,
    0.4, 0.5, 0.6
  ), ncol = 2, byrow = TRUE)
  test_p_YE_given_I <- matrix(c(
    0.1, 0.1, 0.1,
    0.2, 0.2, 0.2
  ), ncol = 2, byrow = TRUE)

  result <- run_trial_simulation(
    test_config,
    test_p_YI,
    test_p_YT_given_I,
    test_p_YE_given_I,
    rho0,
    rho1,
    seed = 11118
  )

  expect_true(result$terminated_early)
  expect_true(is.na(result$final_od))
  expect_true(is.na(result$final_od_index))
  expect_length(result$final_admissible_set, 0)
  expect_length(result$final_admissible_indices, 0)
  expect_length(result$poc_eligible_set, 0)
  expect_length(result$poc_eligible_indices, 0)
  expect_length(result$poc_pairwise_probs, 0)
  expect_length(result$final_candidate_utilities, 0)
  expect_false(result$poc_validated)
  expect_equal(result$poc_probability, 0)
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
