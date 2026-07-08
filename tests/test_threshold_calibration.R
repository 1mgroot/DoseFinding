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

  result <- calibrate_single_threshold(
    param_name = "c_T",
    endpoint = "toxicity",
    scenario = create_threshold_scenario("toxicity", settings),
    candidates = settings$c_T_candidates,
    settings = settings
  )

  expect_equal(result$param_name, "c_T")
  expect_equal(result$optimal_value, 0.55)
  expect_equal(nrow(result$results), 1)
  expect_equal(result$results$c_E, 0)
  expect_equal(result$results$c_I, 0)
  expect_equal(result$selection_metric, "final_admissible_missing_rate")
  expect_true("final_admissible_missing_rate" %in% names(result$results))
  expect_true("target_endpoint_missing_rate" %in% names(result$results))
  expect_true("selected" %in% names(result$results))
  expect_false("target_distance" %in% names(result$results))
})

test_that("threshold calibration reruns full trials for each c candidate", {
  original_runner <- get("run_threshold_calibration_simulation", envir = .GlobalEnv)
  on.exit(assign("run_threshold_calibration_simulation", original_runner, envir = .GlobalEnv), add = TRUE)

  recorded_calls <- data.frame(
    endpoint = character(),
    c_T = numeric(),
    c_I = numeric(),
    c_E = numeric(),
    seed = numeric()
  )
  assign(
    "run_threshold_calibration_simulation",
    function(config, scenario, seed = NULL) {
      recorded_calls <<- rbind(recorded_calls, data.frame(
        endpoint = scenario$endpoint,
        c_T = config$c_T,
        c_I = config$c_I,
        c_E = config$c_E,
        seed = seed
      ))
      list(
        terminated_early = FALSE,
        termination_stage = NA_integer_,
        final_admissible_set = 1L,
        final_admissible_missing = FALSE,
        target_endpoint_missing = FALSE,
        mean_admissible_count = 1,
        total_participants = config$cohort_size * config$n_stages,
        admissibility = data.frame(),
        success = TRUE
      )
    },
    envir = .GlobalEnv
  )

  settings <- default_separate_threshold_settings(quick_mode = TRUE)
  settings$n_sim_per_candidate <- 2
  settings$calibration_seed <- 1000
  settings$show_progress <- FALSE
  all_params <- c("c_T", "c_I", "c_E")
  checks <- list(
    c_T = list(endpoint = "toxicity", candidates = c(0.45, 0.55)),
    c_I = list(endpoint = "immune", candidates = c(0.50, 0.70)),
    c_E = list(endpoint = "efficacy", candidates = c(0.35, 0.50))
  )

  for (param_name in names(checks)) {
    recorded_calls <- recorded_calls[0, , drop = FALSE]
    candidates <- checks[[param_name]]$candidates

    result <- calibrate_single_threshold(
      param_name = param_name,
      endpoint = checks[[param_name]]$endpoint,
      scenario = create_threshold_scenario(checks[[param_name]]$endpoint, settings),
      candidates = candidates,
      settings = settings
    )

    expect_equal(nrow(recorded_calls), length(candidates) * settings$n_sim_per_candidate)
    expect_equal(result$n_simulations, settings$n_sim_per_candidate)
    expect_equal(recorded_calls[[param_name]], rep(candidates, each = settings$n_sim_per_candidate))
    for (inactive_param in setdiff(all_params, param_name)) {
      expect_equal(recorded_calls[[inactive_param]], rep(0, nrow(recorded_calls)))
    }
    seed_stride <- threshold_seed_stride(settings$n_sim_per_candidate)
    expect_equal(
      recorded_calls$seed,
      settings$calibration_seed +
        (rep(seq_along(candidates), each = settings$n_sim_per_candidate) - 1L) * seed_stride +
        rep(seq_len(settings$n_sim_per_candidate), times = length(candidates))
    )
  }

  expect_equal(threshold_seed_stride(100001), 100002)
  expect_equal(threshold_candidate_seed(1000, 2, 1, 100001), 101003)
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

test_that("threshold calibration history records runtime", {
  expect_equal(format_duration_seconds(65), "1m 05s")

  settings <- default_separate_threshold_settings(quick_mode = TRUE)
  result <- list(
    settings = settings,
    calibrations = list(
      c_T = list(
        param_name = "c_T",
        endpoint = "toxicity",
        optimal_value = 0.45,
        achieved_missing_rate = 0.82,
        target_missing_range = c(0.80, 0.90),
        status = "within target range",
        n_simulations = 5
      )
    ),
    recommended_thresholds = list(c_T = 0.45, c_I = 0.50, c_E = 0.35),
    run_duration_seconds = 65
  )

  log_file <- tempfile(fileext = ".md")
  append_threshold_calibration_log(result, file_path = log_file, run_label = "test run")
  log_text <- paste(readLines(log_file, warn = FALSE), collapse = "\n")
  expect_match(log_text, "- Runtime: 1m 05s", fixed = TRUE)
})
