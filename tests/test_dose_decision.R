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

test_that("adaptive_randomization returns probabilities that sum to 1 over admissible doses", {
  posterior_summaries <- list(
    tox = data.frame(pava_mean = c(0.1, 0.2, 0.25, 0.3)),
    eff = data.frame(pava_mean = c(0.6, 0.7, 0.65, 0.7)),
    imm = data.frame(pava_mean = c(0.4, 0.45))
  )
  admissible_set <- c(1, 2)  # assume first two doses admissible
  alloc_probs <- adaptive_randomization(admissible_set, posterior_summaries, trial_config)
  expect_type(alloc_probs, "double")
  expect_equal(sum(alloc_probs[admissible_set]), 1)
  expect_true(all(alloc_probs >= 0))
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
