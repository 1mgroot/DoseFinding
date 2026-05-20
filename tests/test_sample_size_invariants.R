# Test sample-size invariants for the group-sequential design.

library(testthat)
library(dplyr)

if (basename(getwd()) == "tests") {
  setwd("..")
}

source("src/utils/helpers.R")
source("src/core/simulate_data.R")
source("src/core/model_utils.R")
source("src/decision/dose_decision.R")
source("src/core/main.R")

test_that("sample size invariants hold with early termination enabled", {
  test_config <- list(
    dose_levels = c(1, 2, 3, 4, 5),
    n_stages = 5,
    cohort_size = 15,
    phi_T = 0.30,
    c_T = 0.5,
    phi_E = 0.25,
    c_E = 0.5,
    phi_I = 0.20,
    c_I = 0.5,
    delta_poc = 0.8,
    c_poc = 0.7,
    enable_early_termination = TRUE,
    log_early_termination = FALSE,
    verbose_logging = FALSE
  )

  utility_table <- array(0, dim = c(2, 2, 2))
  utility_table[1, 1, 1] <- 0
  utility_table[2, 1, 1] <- 80
  utility_table[1, 2, 1] <- 0
  utility_table[2, 2, 1] <- 30
  utility_table[1, 1, 2] <- 10
  utility_table[2, 1, 2] <- 100
  utility_table[1, 2, 2] <- 0
  utility_table[2, 2, 2] <- 40
  test_config$utility_table <- utility_table

  test_scenario <- list(
    p_YI = c(0.1, 0.2, 0.3, 0.4, 0.5),
    p_YT_given_I = matrix(c(
      0.15, 0.10, 0.10, 0.05, 0.05,
      0.20, 0.15, 0.10, 0.05, 0.05
    ), nrow = 5, ncol = 2),
    p_YE_given_I = matrix(c(
      0.20, 0.25, 0.30, 0.35, 0.40,
      0.30, 0.35, 0.40, 0.45, 0.50
    ), nrow = 5, ncol = 2),
    rho0 = 1.5,
    rho1 = 2.0
  )

  max_n <- test_config$n_stages * test_config$cohort_size

  for (i in seq_len(50)) {
    seed <- 50000 + i
    result <- run_trial_simulation(
      test_config,
      test_scenario$p_YI,
      test_scenario$p_YT_given_I,
      test_scenario$p_YE_given_I,
      test_scenario$rho0,
      test_scenario$rho1,
      seed = seed
    )

    n_enrolled <- nrow(result$all_data)
    expect_true(n_enrolled <= max_n)
    expect_equal(n_enrolled %% test_config$cohort_size, 0)

    if (result$terminated_early) {
      expect_true(result$termination_stage %in% seq_len(test_config$n_stages))
      expect_equal(
        n_enrolled,
        result$termination_stage * test_config$cohort_size
      )
    } else {
      expect_equal(n_enrolled, max_n)
    }
  }
})
