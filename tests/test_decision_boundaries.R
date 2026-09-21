library(testthat)
if (basename(getwd()) == "tests") setwd("..")
source("src/core/config.R")
source("src/core/main.R")
source("src/optimization/threshold_calibration.R")

test_that("strict posterior cutoffs reject decimal equality from calibration grids", {
  cutoff <- seq(.35, .75, by = .05)[3]
  expect_true(.45 > cutoff) # The original implementation's counterexample.
  expect_false(exceeds_probability_cutoff(.45, cutoff))
  expect_identical(exceeds_probability_cutoff(c(.449,.450,.451), cutoff), c(FALSE,FALSE,TRUE))
  expect_true(exceeds_probability_cutoff(.45 + 1e-10, cutoff))
  expect_false(exceeds_probability_cutoff(1, 1))
  expect_false(exceeds_probability_cutoff(0, 0))

  config <- trial_config
  config$dose_levels <- 1
  config$c_T <- cutoff
  config$c_E <- .5
  config$c_I <- .5
  posterior <- list(
    tox_marginal = data.frame(samples = I(list(c(rep(.1,450), rep(.5,550))))),
    eff_marginal = data.frame(samples = I(list(rep(.8,1000)))),
    imm = data.frame(samples_pava = I(list(rep(.8,1000)))))
  expect_length(get_admissible_set(posterior, config, FALSE), 0)
  table <- endpoint_admissibility_table(posterior, config)
  expect_false(table$tox_pass)
  expect_false(table$overall_pass)
})

test_that("PoC uses the same strict cutoff at an exact posterior probability boundary", {
  config <- trial_config
  config$dose_levels <- 1:2
  config$c_poc <- seq(.35,.75,.05)[3]
  config$delta_poc <- .8
  posterior <- list(imm = data.frame(samples_pava = I(list(rep(.2,1000),
    c(rep(.5,450),rep(.1,550))))))
  env <- new.env(parent = environment(calculate_poc_probability))
  env$get_expected_utility <- function(...) 10
  calculate <- calculate_poc_probability
  environment(calculate) <- env
  result <- calculate(2L, posterior, config)
  expect_equal(unname(result$pairwise_probs), .45)
  expect_length(result$P_final, 0)
})

test_that("empty final admissible set means completed without OD rather than early stop", {
  config <- within(trial_config, {
    n_stages <- 1
    c_T <- 1
    verbose_logging <- FALSE
    log_early_termination <- FALSE
  })
  result <- run_trial_simulation(config,p_YI,p_YT_given_I,p_YE_given_I,0,0,11118)
  expect_false(result$terminated_early)
  expect_true(is.na(result$termination_stage))
  expect_equal(nrow(result$all_data), config$cohort_size)
  expect_length(result$final_admissible_set, 0)
  expect_false(result$poc_validated)
  expect_true(is.na(result$final_od))
  config$n_stages <- 2
  early <- run_trial_simulation(config,p_YI,p_YT_given_I,p_YE_given_I,0,0,11118)
  expect_true(early$terminated_early)
  expect_equal(early$termination_stage,1)
  expect_equal(nrow(early$all_data), config$cohort_size)
})
