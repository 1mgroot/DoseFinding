library(testthat)
if (basename(getwd()) == "tests") setwd("..")
source("src/core/config.R")
source("src/core/main.R")

test_that("multinomial allocation is reproducible and preserves the first stage", {
  config <- within(trial_config, {
    verbose_logging <- FALSE
    log_early_termination <- FALSE
    c_T <- 0
    c_E <- 0
    c_I <- 0
    adaptive_allocation_method <- "multinomial"
  })
  run <- function(config, seed = 11118) {
    run_trial_simulation(config, p_YI, p_YT_given_I, p_YE_given_I, 0, 0, seed)
  }
  random <- run(config)
  expect_identical(random, run(config))
  expect_equal(as.integer(table(random$all_data$stage)), rep(15L, 5))
  legacy_config <- config
  legacy_config$adaptive_allocation_method <- "largest_remainder"
  legacy <- run(legacy_config)
  expect_identical(subset(random$all_data, stage == 1),
                   subset(legacy$all_data, stage == 1))
  expect_identical(subset(random$all_alloc_probs, Stage <= 2),
                   subset(legacy$all_alloc_probs, Stage <= 2))
  expect_false(identical(random$all_data, legacy$all_data))
  excluded <- subset(random$all_alloc_probs, Prob == 0)
  for (i in seq_len(nrow(excluded))) {
    expect_false(any(random$all_data$stage == excluded$Stage[i] &
                     random$all_data$d == excluded$Dose[i]))
  }
  legacy_config$adaptive_allocation_method <- NULL
  expect_identical(legacy, run(legacy_config))
  config$adaptive_allocation_method <- "typo"
  expect_error(run(config), "adaptive_allocation_method")
})

test_that("patient counts follow probabilities rather than rounding small arms away", {
  # Exercise the trial's actual allocation branch with fixed adaptive probabilities.
  config <- within(trial_config, {
    verbose_logging <- FALSE
    log_early_termination <- FALSE
    n_stages <- 2
    c_T <- 0
    c_E <- 0
    c_I <- 0
    adaptive_allocation_method <- "multinomial"
  })
  env <- new.env(parent = environment(run_trial_simulation))
  env$adaptive_randomization <- function(...) c(.01, .49, .50, 0, 0)
  run <- run_trial_simulation
  environment(run) <- env
  counts <- sapply(1:60, function(seed) {
    result <- run(config, p_YI, p_YT_given_I, p_YE_given_I, 0, 0, seed)
    tabulate(result$all_data$d[result$all_data$stage == 2], nbins = 5)
  })
  expect_true(all(colSums(counts) == 15L))
  expect_true(all(counts[4:5, ] == 0L))
  expect_gt(sum(counts[1, ]), 0) # Legacy rounding assigns zero every time.
  expect_gt(var(counts[2, ]), 0)
  # Broad, fixed-seed Monte Carlo bound, not an exact reproduction of the code.
  expect_lt(max(abs(rowSums(counts) / sum(counts) - c(.01,.49,.5,0,0))), .07)
})
