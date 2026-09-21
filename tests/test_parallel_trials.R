library(testthat)
if (basename(getwd()) == "tests") setwd("..")
source("src/core/config.R")
source("src/core/main.R")
source("src/optimization/threshold_calibration.R")
source("src/optimization/poc_calibration.R")

test_that("parallel execution preserves seeded trial results and candidate summaries", {
  skip_on_os("windows")
  previous <- Sys.getenv("DOSEFINDING_WORKERS", unset = NA_character_)
  on.exit(if (is.na(previous)) Sys.unsetenv("DOSEFINDING_WORKERS") else
            Sys.setenv(DOSEFINDING_WORKERS = previous))
  config <- within(trial_config, { verbose_logging <- FALSE; log_early_termination <- FALSE })
  run <- function(id) run_trial_simulation(config, p_YI, p_YT_given_I, p_YE_given_I, 0, 0, 11118+id)
  Sys.setenv(DOSEFINDING_WORKERS = 1)
  serial <- map_independent_trials(1:12, run)
  settings <- default_separate_threshold_settings(FALSE)
  settings$show_progress <- FALSE
  scenario <- create_threshold_scenario("toxicity", settings)
  calibrate <- function() calibrate_single_threshold("c_T", "toxicity", scenario,
    c(.35,.45), settings, n_simulations=6)
  serial_cal <- calibrate()
  poc <- function() calibrate_c_poc(create_null_flat_scenario(), base_config = config,
    n_simulations = 6, c_poc_candidates = c(.8,.9), verbose = FALSE, progress = FALSE,
    store_simulation_results = TRUE)
  serial_poc <- poc()
  Sys.setenv(DOSEFINDING_WORKERS = 3)
  expect_identical(map_independent_trials(1:12, run), serial)
  expect_identical(calibrate(), serial_cal)
  parallel_poc <- poc()
  expect_identical(parallel_poc$calibration_results, serial_poc$calibration_results)
  expect_identical(parallel_poc$optimal_c_poc, serial_poc$optimal_c_poc)
})
