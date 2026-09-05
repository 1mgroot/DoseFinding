library(testthat)
library(dplyr)

if (basename(getwd()) == "tests") {
  setwd("..")
}

source("src/core/simulate_data.R")

test_that("dose labels are kept separate from probability array indices", {
  data <- simulate_data_gumbel(
    n_per_dose_vector = c(2, 3),
    dose_levels = c(10, 20),
    p_YI = c(0.2, 0.8),
    p_YT_given_I = matrix(c(0.1, 0.2, 0.2, 0.3), nrow = 2, byrow = TRUE),
    p_YE_given_I = matrix(c(0.3, 0.4, 0.6, 0.7), nrow = 2, byrow = TRUE),
    seed = 11118,
    id_start = 41
  )

  expect_equal(data$dose_index, c(1, 1, 2, 2, 2))
  expect_equal(data$dose_label, c(10, 10, 20, 20, 20))
  expect_equal(data$d, data$dose_label)
  expect_equal(data$id, 41:45)
})

test_that("a multi-stage trial assigns globally unique patient IDs", {
  source("src/core/config.R")
  source("src/core/main.R")

  test_config <- within(trial_config, {
    n_stages <- 3
    cohort_size <- 5
    enable_early_termination <- FALSE
    verbose_logging <- FALSE
    log_early_termination <- FALSE
    c_T <- -1
    c_E <- -1
    c_I <- -1
  })

  result <- run_trial_simulation(
    test_config,
    p_YI,
    p_YT_given_I,
    p_YE_given_I,
    rho0,
    rho1,
    seed = 11118
  )

  expect_equal(result$all_data$id, seq_len(nrow(result$all_data)))
  expect_equal(length(unique(result$all_data$id)), nrow(result$all_data))
  expect_equal(sort(unique(result$all_data$stage)), seq_len(test_config$n_stages))
})
