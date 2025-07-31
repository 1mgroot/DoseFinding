library(testthat)

source("dose_decision.R")
source("config.R")

cat("\n\n==============\nRunning tests for dose_decision.R\n==============\n")

test_that("get_utility_score returns a single value", {
  score <- get_utility_score(0.1, 0.8, 0.9, trial_config)
  expect_length(score, 1)
})

test_that("get_admissible_set returns a vector of dose indices", {
  posterior_summaries <- list(
    tox = data.frame(pava_mean = c(0.1, 0.2, 0.4)),
    eff = data.frame(pava_mean = c(0.3, 0.5, 0.1)),
    imm = data.frame(pava_mean = c(0.3, 0.6, 0.2))
  )
  admissible_set <- get_admissible_set(posterior_summaries, trial_config)
  expect_is(admissible_set, "numeric")
})

test_that("adaptive_randomization returns a vector of probabilities", {
  posterior_summaries <- list(
    tox = data.frame(pava_mean = c(0.1, 0.2, 0.4)),
    eff = data.frame(pava_mean = c(0.3, 0.5, 0.1)),
    imm = data.frame(pava_mean = c(0.3, 0.6, 0.2))
  )
  admissible_set <- get_admissible_set(posterior_summaries, trial_config)
  alloc_probs <- adaptive_randomization(admissible_set, posterior_summaries, trial_config)
  expect_is(alloc_probs, "numeric")
  expect_equal(sum(alloc_probs), 1)
})

test_that("select_final_od returns a single dose index", {
  posterior_summaries <- list(
    tox = data.frame(pava_mean = c(0.1, 0.2, 0.4)),
    eff = data.frame(pava_mean = c(0.3, 0.5, 0.1)),
    imm = data.frame(pava_mean = c(0.3, 0.6, 0.2))
  )
  admissible_set <- get_admissible_set(posterior_summaries, trial_config)
  final_od <- select_final_od(admissible_set, posterior_summaries, trial_config)
  expect_length(final_od, 1)
})

cat("\n==============\nTests for dose_decision.R passed\n==============\n")
