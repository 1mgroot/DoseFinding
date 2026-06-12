# Test separate threshold calibration helpers.

library(testthat)

if (basename(getwd()) == "tests") {
  setwd("..")
}

source("src/optimization/threshold_calibration.R")

test_that("conditional probability matrix preserves requested marginal probabilities", {
  p_immune <- c(0.10, 0.20, 0.30)
  marginal <- c(0.35, 0.45, 0.55)
  conditional <- make_conditional_probability_matrix(
    marginal_probability = marginal,
    p_immune = p_immune,
    immune_effect = 0.10,
    name = "toxicity"
  )

  expect_equal(dim(conditional), c(3, 2))
  expect_true(all(conditional >= 0 & conditional <= 1))
  expect_equal(
    calculate_marginal_probability(p_immune, conditional),
    marginal,
    tolerance = 1e-12
  )
})

test_that("threshold calibration scenarios use endpoint-specific unfavorable settings", {
  settings <- default_separate_threshold_settings(quick_mode = TRUE)

  tox <- create_threshold_scenario("toxicity", settings)
  immune <- create_threshold_scenario("immune", settings)
  eff <- create_threshold_scenario("efficacy", settings)

  expect_equal(tox$p_YI, rep(0.30, 5))
  expect_equal(tox$marginal_p_T, seq(0.30, 0.50, length.out = 5))
  expect_equal(tox$marginal_p_E, rep(0.40, 5))

  expect_equal(immune$p_YI, seq(0.10, 0.15, length.out = 5))
  expect_equal(immune$marginal_p_T, rep(0.15, 5))
  expect_equal(immune$marginal_p_E, rep(0.40, 5))

  expect_equal(eff$p_YI, rep(0.30, 5))
  expect_equal(eff$marginal_p_T, rep(0.15, 5))
  expect_equal(eff$marginal_p_E, seq(0.10, 0.20, length.out = 5))
})

test_that("single threshold calibration returns a structured result", {
  settings <- default_separate_threshold_settings(quick_mode = TRUE)
  settings$n_sim_per_candidate <- 1
  settings$c_T_candidates <- 0.55

  fixed_params <- list(c_T = settings$c_T_start, c_E = settings$c_E_start, c_I = settings$c_I_start)
  result <- calibrate_single_threshold(
    param_name = "c_T",
    endpoint = "toxicity",
    scenario = create_threshold_scenario("toxicity", settings),
    candidates = settings$c_T_candidates,
    fixed_params = fixed_params,
    settings = settings
  )

  expect_equal(result$param_name, "c_T")
  expect_equal(result$optimal_value, 0.55)
  expect_equal(nrow(result$results), 1)
  expect_true("final_admissible_missing_rate" %in% names(result$results))
  expect_true("target_endpoint_missing_rate" %in% names(result$results))
  expect_true("selected" %in% names(result$results))
  expect_false("target_distance" %in% names(result$results))
})

test_that("threshold candidate selection follows c cutoff direction", {
  in_range_table <- data.frame(
    param_value = c(0.45, 0.55, 0.65),
    final_admissible_missing_rate = c(0.78, 0.82, 0.88)
  )
  selected <- select_threshold_candidate(in_range_table, c(0.80, 0.90))
  expect_equal(selected$selected_index, 2)
  expect_match(selected$status, "least strict")

  below_range_table <- data.frame(
    param_value = c(0.45, 0.55, 0.65),
    final_admissible_missing_rate = c(0.40, 0.55, 0.70)
  )
  selected <- select_threshold_candidate(below_range_table, c(0.80, 0.90))
  expect_equal(selected$selected_index, 3)
  expect_match(selected$status, "strictest")

  above_range_table <- data.frame(
    param_value = c(0.45, 0.55, 0.65),
    final_admissible_missing_rate = c(0.92, 0.96, 0.99)
  )
  selected <- select_threshold_candidate(above_range_table, c(0.80, 0.90))
  expect_equal(selected$selected_index, 1)
  expect_match(selected$status, "least strict")
})
