# Separate threshold calibration for c_T, c_I, and c_E.
#
# This module calibrates endpoint-specific admissibility credibility cutoffs
# before PoC calibration. Each cutoff is calibrated in a scenario where only the
# corresponding endpoint is deliberately unfavorable.

library(dplyr)
library(ggplot2)

.threshold_project_root_candidates <- c(".", "..", "../..")
.threshold_project_root_matches <- .threshold_project_root_candidates[
  file.exists(file.path(.threshold_project_root_candidates, "DoseFinding.Rproj"))
]
if (length(.threshold_project_root_matches) == 0) {
  stop("Could not find project root containing DoseFinding.Rproj.")
}
.threshold_project_root <- .threshold_project_root_matches[[1]]

source(file.path(.threshold_project_root, "src/utils/helpers.R"))
source(file.path(.threshold_project_root, "src/core/simulate_data.R"))
source(file.path(.threshold_project_root, "src/core/model_utils.R"))
source(file.path(.threshold_project_root, "src/decision/dose_decision.R"))
source(file.path(.threshold_project_root, "src/core/main.R"))

make_probability_profile <- function(value, n_doses, name = "probability") {
  if (!is.numeric(value)) {
    stop(name, " must be numeric.")
  }
  if (length(value) == 1) {
    profile <- rep(value, n_doses)
  } else if (length(value) == 2) {
    profile <- seq(value[[1]], value[[2]], length.out = n_doses)
  } else if (length(value) == n_doses) {
    profile <- value
  } else {
    stop(name, " must have length 1, 2, or n_doses.")
  }
  if (any(is.na(profile)) || any(profile < 0 | profile > 1)) {
    stop(name, " values must be between 0 and 1.")
  }
  as.numeric(profile)
}

make_conditional_probability_matrix <- function(
  marginal_probability,
  p_immune,
  immune_effect = 0,
  name = "event"
) {
  n_doses <- length(p_immune)
  marginal_probability <- make_probability_profile(
    marginal_probability,
    n_doses,
    paste0(name, "_marginal_probability")
  )
  immune_effect <- make_probability_profile(
    immune_effect,
    n_doses,
    paste0(name, "_immune_effect")
  )

  p_event_given_i0 <- marginal_probability - p_immune * immune_effect
  p_event_given_i1 <- marginal_probability + (1 - p_immune) * immune_effect
  conditional_matrix <- matrix(
    c(p_event_given_i0, p_event_given_i1),
    nrow = n_doses,
    ncol = 2
  )

  if (any(conditional_matrix < 0 | conditional_matrix > 1)) {
    stop(
      name,
      " conditional probabilities are outside [0, 1]. ",
      "Reduce immune_effect or change the marginal probabilities."
    )
  }
  conditional_matrix
}

calculate_marginal_probability <- function(p_immune, conditional_matrix) {
  conditional_matrix[, 1] * (1 - p_immune) + conditional_matrix[, 2] * p_immune
}

make_default_utility_table <- function() {
  utility_table <- array(0, dim = c(2, 2, 2), dimnames = list(
    E = c(0, 1),
    T = c(0, 1),
    I = c(0, 1)
  ))
  utility_table[1, 1, 1] <- 0
  utility_table[2, 1, 1] <- 80
  utility_table[1, 2, 1] <- 0
  utility_table[2, 2, 1] <- 30
  utility_table[1, 1, 2] <- 10
  utility_table[2, 1, 2] <- 100
  utility_table[1, 2, 2] <- 0
  utility_table[2, 2, 2] <- 40
  utility_table
}

default_separate_threshold_settings <- function(quick_mode = TRUE) {
  list(
    quick_mode = quick_mode,
    dose_levels = c(1, 2, 3, 4, 5),
    n_stages = 5,
    cohort_size = 15,
    phi_T = 0.30,
    phi_E = 0.25,
    phi_I = 0.20,
    c_T_start = 0.55,
    c_E_start = 0.50,
    c_I_start = 0.70,
    c_poc = 0.995,
    delta_poc = 0.8,
    rho0 = 1.5,
    rho1 = 2.0,
    target_missing_range = c(0.80, 0.90),
    c_T_candidates = if (quick_mode) c(0.45, 0.55, 0.65) else seq(0.35, 0.75, by = 0.05),
    c_I_candidates = if (quick_mode) c(0.50, 0.70, 0.90) else seq(0.45, 0.95, by = 0.05),
    c_E_candidates = if (quick_mode) c(0.35, 0.50, 0.65) else seq(0.30, 0.80, by = 0.05),
    n_sim_per_candidate = if (quick_mode) 5 else 500,
    calibration_seed = 11118,
    high_tox_p_I = 0.30,
    high_tox_marginal_p_T = c(0.35, 0.60),
    high_tox_marginal_p_E = 0.40,
    low_immune_p_I = c(0.10, 0.15),
    low_immune_marginal_p_T = 0.15,
    low_immune_marginal_p_E = 0.40,
    low_eff_p_I = 0.30,
    low_eff_marginal_p_T = 0.15,
    low_eff_marginal_p_E = c(0.10, 0.20),
    toxicity_immune_effect = 0,
    efficacy_immune_effect = 0,
    output_dir = "results/threshold_calibration"
  )
}

create_threshold_trial_config <- function(settings, c_T, c_E, c_I) {
  config <- list(
    dose_levels = settings$dose_levels,
    n_stages = settings$n_stages,
    cohort_size = settings$cohort_size,
    phi_T = settings$phi_T,
    phi_E = settings$phi_E,
    phi_I = settings$phi_I,
    c_T = c_T,
    c_E = c_E,
    c_I = c_I,
    c_poc = settings$c_poc,
    delta_poc = settings$delta_poc,
    enable_early_termination = TRUE,
    log_early_termination = FALSE,
    verbose_logging = FALSE,
    utility_table = make_default_utility_table()
  )
  config
}

create_threshold_scenario <- function(endpoint, settings) {
  n_doses <- length(settings$dose_levels)

  if (endpoint == "toxicity") {
    p_immune <- make_probability_profile(settings$high_tox_p_I, n_doses, "high_tox_p_I")
    marginal_tox <- make_probability_profile(
      settings$high_tox_marginal_p_T,
      n_doses,
      "high_tox_marginal_p_T"
    )
    marginal_eff <- make_probability_profile(
      settings$high_tox_marginal_p_E,
      n_doses,
      "high_tox_marginal_p_E"
    )
    description <- "High toxicity: marginal P(T) above phi_T; immune and efficacy acceptable"
  } else if (endpoint == "immune") {
    p_immune <- make_probability_profile(settings$low_immune_p_I, n_doses, "low_immune_p_I")
    marginal_tox <- make_probability_profile(
      settings$low_immune_marginal_p_T,
      n_doses,
      "low_immune_marginal_p_T"
    )
    marginal_eff <- make_probability_profile(
      settings$low_immune_marginal_p_E,
      n_doses,
      "low_immune_marginal_p_E"
    )
    description <- "Low immune: P(I) below phi_I; toxicity and efficacy acceptable"
  } else if (endpoint == "efficacy") {
    p_immune <- make_probability_profile(settings$low_eff_p_I, n_doses, "low_eff_p_I")
    marginal_tox <- make_probability_profile(
      settings$low_eff_marginal_p_T,
      n_doses,
      "low_eff_marginal_p_T"
    )
    marginal_eff <- make_probability_profile(
      settings$low_eff_marginal_p_E,
      n_doses,
      "low_eff_marginal_p_E"
    )
    description <- "Low efficacy: marginal P(E) below phi_E; toxicity and immune acceptable"
  } else {
    stop("endpoint must be one of: toxicity, immune, efficacy.")
  }

  p_YT_given_I <- make_conditional_probability_matrix(
    marginal_probability = marginal_tox,
    p_immune = p_immune,
    immune_effect = settings$toxicity_immune_effect,
    name = "toxicity"
  )
  p_YE_given_I <- make_conditional_probability_matrix(
    marginal_probability = marginal_eff,
    p_immune = p_immune,
    immune_effect = settings$efficacy_immune_effect,
    name = "efficacy"
  )

  list(
    endpoint = endpoint,
    p_YI = p_immune,
    p_YT_given_I = p_YT_given_I,
    p_YE_given_I = p_YE_given_I,
    marginal_p_T = calculate_marginal_probability(p_immune, p_YT_given_I),
    marginal_p_E = calculate_marginal_probability(p_immune, p_YE_given_I),
    rho0 = settings$rho0,
    rho1 = settings$rho1,
    scenario_type = paste0("threshold_", endpoint),
    description = description
  )
}

endpoint_admissibility_table <- function(posterior_summaries, config) {
  rows <- lapply(seq_along(config$dose_levels), function(i) {
    tox_prob_safe <- mean(posterior_summaries$tox_marginal$samples[[i]] < config$phi_T)
    eff_prob_good <- mean(posterior_summaries$eff_marginal$samples[[i]] > config$phi_E)
    imm_prob_good <- mean(posterior_summaries$imm$samples_pava[[i]] > config$phi_I)
    data.frame(
      dose = config$dose_levels[[i]],
      dose_index = i,
      tox_prob_safe = tox_prob_safe,
      eff_prob_good = eff_prob_good,
      imm_prob_good = imm_prob_good,
      tox_pass = tox_prob_safe > config$c_T,
      eff_pass = eff_prob_good > config$c_E,
      imm_pass = imm_prob_good > config$c_I,
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  table$overall_pass <- table$tox_pass & table$eff_pass & table$imm_pass
  table
}

run_threshold_calibration_simulation <- function(config, scenario, seed = NULL) {
  results <- run_trial_simulation(
    trial_config = config,
    p_YI = scenario$p_YI,
    p_YT_given_I = scenario$p_YT_given_I,
    p_YE_given_I = scenario$p_YE_given_I,
    rho0 = scenario$rho0,
    rho1 = scenario$rho1,
    seed = seed
  )
  admissibility <- endpoint_admissibility_table(results$posterior_summaries, config)
  final_admissible_set <- admissibility$dose_index[admissibility$overall_pass]
  target_pass_column <- switch(
    scenario$endpoint,
    toxicity = "tox_pass",
    immune = "imm_pass",
    efficacy = "eff_pass"
  )

  list(
    terminated_early = results$terminated_early,
    termination_stage = results$termination_stage,
    final_admissible_set = final_admissible_set,
    final_admissible_missing = length(final_admissible_set) == 0,
    target_endpoint_missing = !any(admissibility[[target_pass_column]]),
    mean_admissible_count = length(final_admissible_set),
    total_participants = nrow(results$all_data),
    admissibility = admissibility,
    success = TRUE
  )
}

summarise_threshold_runs <- function(simulation_results) {
  final_missing <- vapply(
    simulation_results,
    function(x) isTRUE(x$final_admissible_missing),
    logical(1)
  )
  endpoint_missing <- vapply(
    simulation_results,
    function(x) isTRUE(x$target_endpoint_missing),
    logical(1)
  )
  early_stop <- vapply(
    simulation_results,
    function(x) isTRUE(x$terminated_early),
    logical(1)
  )
  admissible_count <- vapply(
    simulation_results,
    function(x) x$mean_admissible_count,
    numeric(1)
  )
  n_sim <- length(simulation_results)
  missing_rate <- mean(final_missing)
  missing_se <- sqrt(missing_rate * (1 - missing_rate) / n_sim)

  list(
    final_admissible_missing_rate = missing_rate,
    missing_rate_se = missing_se,
    missing_rate_ci_lower = max(0, missing_rate - 1.96 * missing_se),
    missing_rate_ci_upper = min(1, missing_rate + 1.96 * missing_se),
    target_endpoint_missing_rate = mean(endpoint_missing),
    early_stop_rate = mean(early_stop),
    mean_admissible_count = mean(admissible_count),
    n_simulations = n_sim
  )
}

calibrate_single_threshold <- function(
  param_name,
  endpoint,
  scenario,
  candidates,
  fixed_params,
  settings,
  n_simulations = settings$n_sim_per_candidate,
  target_missing_range = settings$target_missing_range,
  base_seed = settings$calibration_seed
) {
  if (!param_name %in% c("c_T", "c_I", "c_E")) {
    stop("param_name must be one of c_T, c_I, c_E.")
  }
  if (length(candidates) == 0) {
    stop("candidates must not be empty.")
  }

  candidate_results <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    params <- fixed_params
    params[[param_name]] <- candidates[[i]]
    config <- create_threshold_trial_config(
      settings,
      c_T = params$c_T,
      c_E = params$c_E,
      c_I = params$c_I
    )
    simulation_results <- lapply(seq_len(n_simulations), function(sim_index) {
      run_threshold_calibration_simulation(
        config = config,
        scenario = scenario,
        seed = base_seed + i * 100000 + sim_index
      )
    })
    summary <- summarise_threshold_runs(simulation_results)
    candidate_results[[i]] <- c(
      list(
        param_name = param_name,
        endpoint = endpoint,
        param_value = candidates[[i]],
        c_T = params$c_T,
        c_E = params$c_E,
        c_I = params$c_I
      ),
      summary
    )
  }

  result_table <- bind_rows(lapply(candidate_results, as.data.frame))
  target_center <- mean(target_missing_range)
  in_range <- which(
    result_table$final_admissible_missing_rate >= target_missing_range[[1]] &
      result_table$final_admissible_missing_rate <= target_missing_range[[2]]
  )
  if (length(in_range) > 0) {
    selected_index <- in_range[
      which.min(abs(result_table$final_admissible_missing_rate[in_range] - target_center))
    ]
    status <- "within target range"
  } else {
    selected_index <- which.min(abs(result_table$final_admissible_missing_rate - target_center))
    status <- "closest to target range"
  }

  list(
    param_name = param_name,
    endpoint = endpoint,
    scenario = scenario,
    candidates = candidates,
    results = result_table,
    optimal_value = result_table$param_value[[selected_index]],
    achieved_missing_rate = result_table$final_admissible_missing_rate[[selected_index]],
    target_missing_range = target_missing_range,
    status = status,
    selected_row = selected_index,
    n_simulations = n_simulations
  )
}

calibrate_separate_thresholds <- function(settings = default_separate_threshold_settings()) {
  scenarios <- list(
    c_T = create_threshold_scenario("toxicity", settings),
    c_I = create_threshold_scenario("immune", settings),
    c_E = create_threshold_scenario("efficacy", settings)
  )

  current_params <- list(
    c_T = settings$c_T_start,
    c_E = settings$c_E_start,
    c_I = settings$c_I_start
  )

  c_T_result <- calibrate_single_threshold(
    param_name = "c_T",
    endpoint = "toxicity",
    scenario = scenarios$c_T,
    candidates = settings$c_T_candidates,
    fixed_params = current_params,
    settings = settings
  )
  current_params$c_T <- c_T_result$optimal_value

  c_I_result <- calibrate_single_threshold(
    param_name = "c_I",
    endpoint = "immune",
    scenario = scenarios$c_I,
    candidates = settings$c_I_candidates,
    fixed_params = current_params,
    settings = settings,
    base_seed = settings$calibration_seed + 1000000
  )
  current_params$c_I <- c_I_result$optimal_value

  c_E_result <- calibrate_single_threshold(
    param_name = "c_E",
    endpoint = "efficacy",
    scenario = scenarios$c_E,
    candidates = settings$c_E_candidates,
    fixed_params = current_params,
    settings = settings,
    base_seed = settings$calibration_seed + 2000000
  )
  current_params$c_E <- c_E_result$optimal_value

  list(
    settings = settings,
    scenarios = scenarios,
    calibrations = list(c_T = c_T_result, c_I = c_I_result, c_E = c_E_result),
    recommended_thresholds = current_params
  )
}

threshold_calibration_summary_table <- function(calibration_results) {
  bind_rows(lapply(calibration_results$calibrations, function(result) {
    data.frame(
      parameter = result$param_name,
      endpoint = result$endpoint,
      selected_value = result$optimal_value,
      achieved_missing_rate = result$achieved_missing_rate,
      target_low = result$target_missing_range[[1]],
      target_high = result$target_missing_range[[2]],
      status = result$status,
      n_simulations = result$n_simulations,
      stringsAsFactors = FALSE
    )
  }))
}

append_threshold_calibration_log <- function(
  calibration_results,
  file_path = file.path(calibration_results$settings$output_dir, "threshold_calibration_history.md"),
  run_label = NULL
) {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  summary_table <- threshold_calibration_summary_table(calibration_results)
  header_lines <- character(0)
  if (!file.exists(file_path) || file.info(file_path)$size == 0) {
    header_lines <- c(
      "# Threshold Calibration History",
      "",
      "This log records separate c_T, c_I, and c_E calibration runs before PoC calibration.",
      ""
    )
  }
  title_suffix <- if (!is.null(run_label) && nchar(run_label) > 0) {
    paste0(" - ", run_label)
  } else {
    ""
  }
  lines <- c(
    header_lines,
    paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), title_suffix),
    "",
    paste0(
      "- Target final admissible set missing rate: ",
      calibration_results$settings$target_missing_range[[1]] * 100,
      "% to ",
      calibration_results$settings$target_missing_range[[2]] * 100,
      "%"
    ),
    paste0("- Simulations per candidate: ", calibration_results$settings$n_sim_per_candidate),
    paste0("- Recommended `c_T`: ", calibration_results$recommended_thresholds$c_T),
    paste0("- Recommended `c_I`: ", calibration_results$recommended_thresholds$c_I),
    paste0("- Recommended `c_E`: ", calibration_results$recommended_thresholds$c_E),
    "",
    "| parameter | endpoint | selected value | missing rate | target | status |",
    "|---|---|---:|---:|---:|---|"
  )
  for (i in seq_len(nrow(summary_table))) {
    row <- summary_table[i, , drop = FALSE]
    lines <- c(
      lines,
      paste0(
        "| ", row$parameter,
        " | ", row$endpoint,
        " | ", round(row$selected_value, 3),
        " | ", round(row$achieved_missing_rate * 100, 1), "%",
        " | ", round(row$target_low * 100), "%-",
        round(row$target_high * 100), "%",
        " | ", row$status,
        " |"
      )
    )
  }
  lines <- c(lines, "")
  write(lines, file = file_path, append = TRUE)
  invisible(file_path)
}

plot_threshold_calibration <- function(calibration_result) {
  df <- calibration_result$results
  ggplot(df, aes(x = param_value, y = final_admissible_missing_rate)) +
    geom_rect(
      aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = calibration_result$target_missing_range[[1]],
        ymax = calibration_result$target_missing_range[[2]]
      ),
      fill = "lightgreen",
      alpha = 0.25
    ) +
    geom_ribbon(
      aes(ymin = missing_rate_ci_lower, ymax = missing_rate_ci_upper),
      fill = "#2E86AB",
      alpha = 0.18
    ) +
    geom_line(color = "#2E86AB", linewidth = 1) +
    geom_point(color = "#2E86AB", size = 2.5) +
    geom_vline(
      xintercept = calibration_result$optimal_value,
      linetype = "dashed",
      color = "red"
    ) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(
      title = paste0(calibration_result$param_name, " Calibration"),
      subtitle = paste0(
        "Endpoint: ",
        calibration_result$endpoint,
        "; selected ",
        calibration_result$param_name,
        " = ",
        calibration_result$optimal_value
      ),
      x = calibration_result$param_name,
      y = "Final admissible set missing rate"
    ) +
    theme_bw(base_size = 13)
}
