library(testthat)

# Set working directory to project root for proper path resolution
if (basename(getwd()) == "tests") {
  setwd("..")
}

source("src/core/config.R")
source("src/decision/dose_decision.R")

test_that("get_expected_utility returns a single numeric value", {
  posterior_summaries <- list(
    tox = data.frame(pava_mean = c(0.1, 0.2)),
    eff = data.frame(pava_mean = c(0.6, 0.7)),
    imm = data.frame(pava_mean = c(0.3))
  )
  score <- get_expected_utility(1, posterior_summaries, trial_config)
  expect_type(score, "double")
  expect_length(score, 1)
})

test_that("get_expected_utility averages utility over posterior draws", {
  p_i <- c(0.20, 0.80, 0.50, 0.65)
  p_t0 <- c(0.10, 0.30, 0.15, 0.40)
  p_t1 <- c(0.20, 0.50, 0.25, 0.55)
  p_e0 <- c(0.20, 0.80, 0.45, 0.70)
  p_e1 <- c(0.30, 0.90, 0.60, 0.85)
  posterior_summaries <- list(
    tox = data.frame(
      pava_mean = c(mean(p_t0), mean(p_t1)),
      samples_pava = I(list(p_t0, p_t1))
    ),
    eff = data.frame(
      pava_mean = c(mean(p_e0), mean(p_e1)),
      samples_pava = I(list(p_e0, p_e1))
    ),
    imm = data.frame(
      pava_mean = mean(p_i),
      samples_pava = I(list(p_i))
    )
  )

  draw_utilities <- vapply(seq_along(p_i), function(i) {
    calculate_utility_from_probabilities(
      p_i[[i]], p_t0[[i]], p_t1[[i]], p_e0[[i]], p_e1[[i]],
      trial_config$utility_table
    )
  }, numeric(1))
  plugin_utility <- calculate_utility_from_probabilities(
    mean(p_i), mean(p_t0), mean(p_t1), mean(p_e0), mean(p_e1),
    trial_config$utility_table
  )

  expect_equal(get_expected_utility_draws(1, posterior_summaries, trial_config), draw_utilities)
  expect_equal(get_expected_utility(1, posterior_summaries, trial_config), mean(draw_utilities))
  expect_gt(abs(get_expected_utility(1, posterior_summaries, trial_config) - plugin_utility), 0.01)
})

test_that("get_admissible_set returns numeric indices with expected structure", {
  # Create test config with 3 doses to match test data
  test_config <- trial_config
  test_config$dose_levels <- c(1, 2, 3)
  
  posterior_summaries <- list(
    tox_marginal = data.frame(
      marginal_prob = c(0.1, 0.15, 0.2),
      samples = I(list(rep(0.1, 20), rep(0.15, 20), rep(0.2, 20)))
    ),
    eff_marginal = data.frame(
      marginal_prob = c(0.4, 0.5, 0.6),
      samples = I(list(rep(0.5, 20), rep(0.55, 20), rep(0.6, 20)))
    ),
    imm = data.frame(
      pava_mean = c(0.4, 0.5, 0.6),
      samples_pava = I(list(rep(0.5, 20), rep(0.55, 20), rep(0.6, 20)))
    )
  )
  admissible_set <- get_admissible_set(posterior_summaries, test_config, verbose = FALSE)
  expect_true(is.numeric(admissible_set))
  expect_true(all(admissible_set %in% seq_along(test_config$dose_levels) | length(admissible_set) == 0))
})

test_that("posterior optimality splits ties within each posterior draw", {
  utility_draw_matrix <- matrix(
    c(
      5, 3,
      2, 4,
      7, 7,
      8, 1
    ),
    ncol = 2,
    byrow = TRUE
  )

  expect_equal(
    posterior_optimality_from_utility_draws(utility_draw_matrix),
    c(0.625, 0.375)
  )
})

test_that("adaptive_randomization uses posterior probability of being optimal", {
  utility_table <- array(0, dim = c(2, 2, 2))
  utility_table[2, 1, 1] <- 100
  test_config <- trial_config
  test_config$dose_levels <- c(1, 2)
  test_config$utility_table <- utility_table

  dose1_eff <- c(0.9, 0.1, 0.5, 0.8)
  dose2_eff <- c(0.2, 0.8, 0.5, 0.3)
  posterior_summaries <- list(
    tox = data.frame(
      pava_mean = c(0, 0, 0, 0),
      samples_pava = I(list(rep(0, 4), rep(0, 4), rep(0, 4), rep(0, 4)))
    ),
    eff = data.frame(
      pava_mean = c(mean(dose1_eff), mean(dose1_eff), mean(dose2_eff), mean(dose2_eff)),
      samples_pava = I(list(dose1_eff, dose1_eff, dose2_eff, dose2_eff))
    ),
    imm = data.frame(
      pava_mean = c(0, 0),
      samples_pava = I(list(rep(0, 4), rep(0, 4)))
    )
  )
  admissible_set <- c(1, 2)

  alloc_probs <- adaptive_randomization(admissible_set, posterior_summaries, test_config)
  expect_type(alloc_probs, "double")
  expect_equal(alloc_probs[admissible_set], c(0.625, 0.375))
  expect_equal(sum(alloc_probs[admissible_set]), 1)
  expect_true(all(alloc_probs[-admissible_set] == 0))
})

test_that("select_final_od returns a single dose index from admissible set", {
  posterior_summaries <- list(
    tox = data.frame(pava_mean = c(0.1, 0.2, 0.25, 0.3)),
    eff = data.frame(pava_mean = c(0.6, 0.7, 0.65, 0.7)),
    imm = data.frame(pava_mean = c(0.4, 0.45))
  )
  admissible_set <- c(1, 2)
  final_od <- select_final_od(admissible_set, posterior_summaries, trial_config)
  expect_length(final_od, 1)
  expect_true(final_od %in% admissible_set | is.na(final_od))
})
