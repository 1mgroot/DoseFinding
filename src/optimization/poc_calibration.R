# PoC Calibration System
# This script implements the null/flat scenario calibration methodology.

library(dplyr)
library(ggplot2)

# Source required functions relative to the project root.
# NOTE: config.R should be sourced only once at the top level (e.g., notebook or main script).
project_root_candidates <- c(".", "..", "../..")
project_root_matches <- project_root_candidates[
  file.exists(file.path(project_root_candidates, "DoseFinding.Rproj"))
]
if (length(project_root_matches) == 0) {
  stop("Could not find project root containing DoseFinding.Rproj.")
}
project_root <- project_root_matches[[1]]
source(file.path(project_root, "src/utils/helpers.R"))
source(file.path(project_root, "src/core/simulate_data.R"))
source(file.path(project_root, "src/core/model_utils.R"))
source(file.path(project_root, "src/decision/dose_decision.R"))
source(file.path(project_root, "src/core/main.R"))

# Create null/flat scenario parameters
create_null_flat_scenario <- function(
  n_doses = 5,
  phi_I = 0.2,          # Immune response threshold (flat across all doses)
  phi_E = 0.25,         # Marginal efficacy threshold (flat across all doses)
  tox_upper = 0.30,     # Toxicity upper bound
  tox_flat = 0.05       # Flat toxicity rate (safe level)
) {
  # Set P_I = (φ_I, φ_I, ..., φ_I) - all doses have same immune response
  p_YI <- rep(phi_I, n_doses)

  # Calculate P_E using total probability formula:
  # P_E(j) = P(E|I=0,j) * (1-P_I(j)) + P(E|I=1,j) * P_I(j)
  # To achieve marginal P_E = φ_E at all doses, set both conditional probabilities to φ_E
  # This is the simplest solution: 0.25 = 0.25 * (1-0.2) + 0.25 * 0.2
  p_YE_given_I <- matrix(rep(phi_E, n_doses * 2), nrow = n_doses, ncol = 2)

  # Set P_T = (tox_flat, tox_flat, ..., tox_flat) - all doses have same safe toxicity
  p_YT_given_I <- matrix(rep(tox_flat, n_doses * 2), nrow = n_doses, ncol = 2)

  # Correlation parameters (can be kept constant)
  rho0 <- 1.5
  rho1 <- 2

  return(list(
    p_YI = p_YI,
    p_YT_given_I = p_YT_given_I,
    p_YE_given_I = p_YE_given_I,
    rho0 = rho0,
    rho1 = rho1,
    scenario_type = "null_flat",
    description = paste0("Null/Flat: φ_I=", phi_I, ", φ_E=", phi_E, ", tox=", tox_flat)
  ))
}

# Run single simulation and extract metrics
run_single_calibration_simulation <- function(config, scenario_params, seed = NULL) {
  # Generate simulation-specific seed
  # Use seed directly as base seed for run_trial_simulation
  # run_trial_simulation will handle stage-specific seeds internally

  tryCatch({
    results <- run_trial_simulation(
      within(config, {
        verbose_logging <- FALSE
        log_early_termination <- FALSE
      }),
      scenario_params$p_YI,
      scenario_params$p_YT_given_I,
      scenario_params$p_YE_given_I,
      scenario_params$rho0,
      scenario_params$rho1,
      seed = seed
    )

    true_optimal_dose <- if (!is.null(scenario_params$true_optimal_dose)) {
      scenario_params$true_optimal_dose
    } else {
      NA_integer_
    }

    # Extract key metrics
    metrics <- list(
      terminated_early = results$terminated_early,
      termination_stage = ifelse(results$terminated_early, results$termination_stage, NA),
      final_od = ifelse(results$terminated_early, NA, results$final_od),
      poc_validated = ifelse(results$terminated_early, FALSE, results$poc_validated),
      poc_probability = ifelse(results$terminated_early, 0, results$poc_probability),
      total_participants = nrow(results$all_data),
      true_optimal_selected = if (!results$terminated_early && !is.na(true_optimal_dose)) {
        results$final_od == true_optimal_dose
      } else {
        NA
      }
    )

    # Allocation summary by dose
    allocation_summary <- results$all_data %>%
      group_by(d) %>%
      summarise(n_participants = n(), .groups = 'drop') %>%
      arrange(d)

    # Stage-by-stage allocation
    stage_allocation <- results$all_data %>%
      group_by(d, stage) %>%
      summarise(n_participants = n(), .groups = 'drop') %>%
      arrange(d, stage)

    return(list(
      metrics = metrics,
      allocation_summary = allocation_summary,
      stage_allocation = stage_allocation,
      success = TRUE,
      debug_info = list(
        terminated_early = results$terminated_early,
        termination_stage = ifelse(results$terminated_early, results$termination_stage, NA),
        termination_reason = ifelse(results$terminated_early, results$termination_reason, NA),
        posterior_summaries = results$posterior_summaries,
        all_data = results$all_data
      )
    ))

  }, error = function(e) {
    cat("Error in simulation:", e$message, "\n")
    return(list(
      metrics = list(
        terminated_early = TRUE,
        termination_stage = 1,
        final_od = NA,
        poc_validated = FALSE,
        poc_probability = 0,
        total_participants = 0,
        true_optimal_selected = FALSE
      ),
      allocation_summary = data.frame(d = 1:5, n_participants = 0),
      stage_allocation = data.frame(d = integer(0), stage = integer(0), n_participants = integer(0)),
      success = FALSE,
      debug_info = list(
        terminated_early = TRUE,
        termination_stage = 1,
        termination_reason = "Error in simulation",
        posterior_summaries = NULL,
        all_data = data.frame()
      )
    ))
  })
}

# Detailed debug logging for early termination cases
log_early_termination_context <- function(debug_info, config) {
  cat("\n--- Early Termination Debug ---\n")
  cat("Stage:", debug_info$termination_stage, "\n")
  if (!is.null(debug_info$termination_reason)) {
    cat("Reason:", debug_info$termination_reason, "\n")
  }

  # Allocation summaries
  if (!is.null(debug_info$all_data) && nrow(debug_info$all_data) > 0) {
    cat("Allocation by dose (cumulative):\n")
    alloc_dose <- as.data.frame(table(debug_info$all_data$d))
    if (nrow(alloc_dose) > 0) {
      for (r in seq_len(nrow(alloc_dose))) {
        cat("  Dose", alloc_dose[r, 1], ":", alloc_dose[r, 2], "\n")
      }
    }
    last_stage <- max(debug_info$all_data$stage, na.rm = TRUE)
    cat("Last stage:", last_stage, "allocation by dose:\n")
    alloc_last <- debug_info$all_data[debug_info$all_data$stage == last_stage, , drop = FALSE]
    alloc_last_tab <- as.data.frame(table(alloc_last$d))
    if (nrow(alloc_last_tab) > 0) {
      for (r in seq_len(nrow(alloc_last_tab))) {
        cat("  Dose", alloc_last_tab[r, 1], ":", alloc_last_tab[r, 2], "\n")
      }
    }
  } else {
    cat("No participant data collected.\n")
  }

  # Posterior marginal means
  ps <- debug_info$posterior_summaries
  if (!is.null(ps)) {
    cat("Posterior marginal means (at termination stage):\n")
    cat("  Tox:", paste(round(ps$tox_marginal$marginal_prob, 3), collapse = ", "), "\n")
    cat("  Eff:", paste(round(ps$eff_marginal$marginal_prob, 3), collapse = ", "), "\n")
    cat("  Imm:", paste(round(ps$imm$pava_mean, 3), collapse = ", "), "\n")

    cat("\nAdmissibility re-check at termination:\n")
    adm <- get_admissible_set(ps, config, verbose = TRUE)
    cat("Computed admissible set:", adm, "\n")
  }
  cat("--- End Early Termination Debug ---\n")
}

# Main calibration function for C_poc
poc_format_duration <- function(seconds) {
  if (is.null(seconds) || length(seconds) == 0 || is.na(seconds) || is.infinite(seconds)) {
    return("unknown")
  }
  seconds <- max(0, as.numeric(seconds))
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  secs <- round(seconds %% 60)
  if (hours > 0) {
    return(sprintf("%dh %02dm %02ds", hours, minutes, secs))
  }
  if (minutes > 0) {
    return(sprintf("%dm %02ds", minutes, secs))
  }
  sprintf("%ds", secs)
}

poc_progress_log <- function(..., enabled = TRUE) {
  if (!isTRUE(enabled)) {
    return(invisible(NULL))
  }
  cat(
    paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""), "\n"),
    file = stderr()
  )
  flush.console()
  invisible(NULL)
}

calibrate_c_poc <- function(
  null_scenario,
  c_poc_candidates = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
  n_simulations = 1000,
  base_config = NULL,
  target_rate = 0.10,
  debug_early_termination = FALSE,  # Changed default to FALSE to reduce console output
  max_debug_cases_per_candidate = 2,  # Reduced from 3 to 2
  verbose = TRUE,
  progress = FALSE,
  progress_interval_seconds = 300,
  progress_prefix = NULL,
  store_simulation_results = TRUE,
  common_random_numbers = TRUE,
  calibration_seed = 10000
) {
  if (missing(null_scenario) || is.null(null_scenario)) {
    stop("null_scenario is required. Use create_null_flat_scenario() to create one.")
  }
  if (length(c_poc_candidates) == 0) {
    stop("c_poc_candidates must contain at least one value.")
  }
  if (any(c_poc_candidates < 0 | c_poc_candidates > 1)) {
    stop("c_poc_candidates must be between 0 and 1.")
  }
  if (!is.numeric(n_simulations) || length(n_simulations) != 1 || n_simulations < 1) {
    stop("n_simulations must be a positive integer.")
  }
  if (!is.numeric(target_rate) || length(target_rate) != 1 || target_rate <= 0 || target_rate >= 1) {
    stop("target_rate must be between 0 and 1.")
  }
  if (!is.numeric(progress_interval_seconds) || length(progress_interval_seconds) != 1 || progress_interval_seconds < 1) {
    stop("progress_interval_seconds must be a positive number of seconds.")
  }
  if (!is.logical(store_simulation_results) || length(store_simulation_results) != 1) {
    stop("store_simulation_results must be TRUE or FALSE.")
  }
  if (!is.logical(common_random_numbers) || length(common_random_numbers) != 1) {
    stop("common_random_numbers must be TRUE or FALSE.")
  }
  if (!is.null(calibration_seed) && (
    !is.numeric(calibration_seed) ||
      length(calibration_seed) != 1 ||
      is.na(calibration_seed) ||
      calibration_seed < 0
  )) {
    stop("calibration_seed must be NULL or a non-negative number.")
  }
  if (!is.null(calibration_seed)) {
    calibration_seed <- as.numeric(calibration_seed)
  }
  if (isTRUE(common_random_numbers) && is.null(calibration_seed)) {
    stop("calibration_seed cannot be NULL when common_random_numbers is TRUE.")
  }
  if (!isTRUE(verbose)) {
    quiet_connection <- textConnection("calibration_output", "w", local = TRUE)
    sink(quiet_connection)
    on.exit({
      sink()
      close(quiet_connection)
    }, add = TRUE)
  }

  cat("Starting C_poc calibration...\n")
  cat("Testing", length(c_poc_candidates), "C_poc values\n")
  cat("Simulations per value:", n_simulations, "\n")
  cat(
    "Seed strategy:",
    if (isTRUE(common_random_numbers)) {
      "common random numbers across C_poc candidates"
    } else {
      "independent seeds across C_poc candidates"
    },
    "\n"
  )
  cat("Calibration seed:", if (is.null(calibration_seed)) "NULL" else calibration_seed, "\n")
  calibration_started_at <- Sys.time()
  next_progress_at <- calibration_started_at + progress_interval_seconds
  total_calibration_simulations <- length(c_poc_candidates) * n_simulations
  progress_label <- if (!is.null(progress_prefix) && nchar(progress_prefix) > 0) {
    paste0(progress_prefix, " ")
  } else {
    ""
  }

  # Default configuration if not provided
  if (is.null(base_config)) {
    base_config <- list(
      dose_levels = c(1, 2, 3, 4, 5),
      n_stages = 5,
      cohort_size = 15,
      phi_T = 0.30,
      c_T = 0.3,  # Relaxed for calibration (null scenario at boundary)
      phi_E = 0.20,
      c_E = 0.3,  # Relaxed for calibration (null scenario at boundary)
      phi_I = 0.20,
      c_I = 0.3,  # Relaxed for calibration (null scenario at boundary)
      delta_poc = 0.8,
      enable_early_termination = TRUE,
      log_early_termination = FALSE,
      verbose_logging = FALSE
    )

    # Add utility table
    utility_table <- array(0, dim = c(2, 2, 2))
    utility_table[1, 1, 1] <- 0   # E=0, T=0, I=0
    utility_table[2, 1, 1] <- 80  # E=1, T=0, I=0
    utility_table[1, 2, 1] <- 0   # E=0, T=1, I=0
    utility_table[2, 2, 1] <- 30  # E=1, T=1, I=0
    utility_table[1, 1, 2] <- 10  # E=0, T=0, I=1
    utility_table[2, 1, 2] <- 100 # E=1, T=0, I=1
    utility_table[1, 2, 2] <- 0   # E=0, T=1, I=1
    utility_table[2, 2, 2] <- 40  # E=1, T=1, I=1
    base_config$utility_table <- utility_table
  }

  # Print configuration parameters
  cat("\n=== Configuration Parameters Used ===\n")
  cat("Trial Design Parameters:\n")
  cat("  dose_levels:", paste(base_config$dose_levels, collapse = ", "), "\n")
  cat("  n_stages:", base_config$n_stages, "\n")
  cat("  cohort_size:", base_config$cohort_size, "\n")
  cat("\nToxicity Parameters:\n")
  cat("  phi_T:", base_config$phi_T, "\n")
  cat("  c_T:", base_config$c_T, "\n")
  cat("\nEfficacy Parameters:\n")
  cat("  phi_E:", base_config$phi_E, "\n")
  cat("  c_E:", base_config$c_E, "\n")
  cat("\nImmune Response Parameters:\n")
  cat("  phi_I:", base_config$phi_I, "\n")
  cat("  c_I:", base_config$c_I, "\n")
  cat("\nPoC Parameters:\n")
  cat("  delta_poc:", base_config$delta_poc, "\n")
  cat("  enable_early_termination:", base_config$enable_early_termination, "\n")
  cat("  log_early_termination:", base_config$log_early_termination, "\n")
  cat("\n=== Null Scenario Parameters ===\n")
  cat("  p_YI:", paste(null_scenario$p_YI, collapse = ", "), "\n")
  cat("  p_YT_given_I (rows=dose, cols=I=[0,1]):\n")
  for (i in seq_len(nrow(null_scenario$p_YT_given_I))) {
    cat("    Dose", i, ":", paste(null_scenario$p_YT_given_I[i,], collapse = ", "), "\n")
  }
  cat("  p_YE_given_I (rows=dose, cols=I=[0,1]):\n")
  for (i in seq_len(nrow(null_scenario$p_YE_given_I))) {
    cat("    Dose", i, ":", paste(null_scenario$p_YE_given_I[i,], collapse = ", "), "\n")
  }
  cat("  rho0:", null_scenario$rho0, "\n")
  cat("  rho1:", null_scenario$rho1, "\n")
  cat("\n")

  calibration_results <- list()
  candidate_seed_stride <- 100000

  simulation_seed <- function(candidate_index, simulation_index) {
    if (is.null(calibration_seed)) {
      return(NULL)
    }
    candidate_offset <- if (isTRUE(common_random_numbers)) {
      0
    } else {
      (candidate_index - 1) * candidate_seed_stride
    }
    calibration_seed + candidate_offset + simulation_index
  }

  for (i in seq_along(c_poc_candidates)) {
    c_poc <- c_poc_candidates[i]
    cat("=== Testing C_poc =", c_poc, "===\n")

    # Update config with current C_poc
    config <- base_config
    config$c_poc <- c_poc
    config$verbose_logging <- FALSE
    config$log_early_termination <- FALSE

    # Run simulations
    simulation_results <- if (isTRUE(store_simulation_results)) list() else NULL
    poc_detection_count <- 0
    early_termination_count <- 0
    n_completed <- 0
    debug_count <- 0

    for (sim in 1:n_simulations) {
      result <- run_single_calibration_simulation(
        config,
        null_scenario,
        seed = simulation_seed(i, sim)
      )
      if (isTRUE(store_simulation_results)) {
        simulation_results[[sim]] <- result
      }

      # Count PoC detection: ONLY when trial completes AND P_final is non-empty
      # Early terminated trials MUST have poc_validated = FALSE
      if (result$metrics$terminated_early) {
        early_termination_count <- early_termination_count + 1
      } else {
        n_completed <- n_completed + 1
      }
      if (!result$metrics$terminated_early && result$metrics$poc_validated) {
        poc_detection_count <- poc_detection_count + 1
      }

      # Debug output for first 3 reps of first c_poc candidate (reduced from 10 to avoid console overflow)
      if (i == 1 && sim <= 3) {
        cat("\n[DEBUG] c_poc =", c_poc, ", sim =", sim, "\n")
        cat("  Early terminated:", result$metrics$terminated_early, "\n")
        if (!result$metrics$terminated_early) {
          cat("  Final OD:", result$metrics$final_od, "\n")
          cat("  PoC validated:", result$metrics$poc_validated, "\n")
          cat("  PoC probability:", round(result$metrics$poc_probability, 3), "\n")
          # Show A_final, pairwise_probs, P_final from debug_info
          ps <- result$debug_info$posterior_summaries
          if (!is.null(ps)) {
            adm <- get_admissible_set(ps, config, verbose = FALSE)
            cat("  A_final:", adm, "\n")
            # We need to extract P_final info - this requires running the calculation again
            # For now, just show if PoC was validated (which means P_final was non-empty)
            cat("  P_final non-empty:", result$metrics$poc_validated, "\n")
          }
        } else {
          cat("  Termination stage:", result$metrics$termination_stage, "\n")
        }
      } else if (result$metrics$terminated_early && debug_early_termination && debug_count < max_debug_cases_per_candidate) {
        cat("\n[DEBUG] Early termination example (C_poc =", c_poc, ", sim =", sim, ")\n")
        log_early_termination_context(result$debug_info, config)
        debug_count <- debug_count + 1
      }

      now <- Sys.time()
      if (isTRUE(progress) && now >= next_progress_at) {
        completed_calibration_simulations <- (i - 1) * n_simulations + sim
        elapsed_seconds <- as.numeric(difftime(now, calibration_started_at, units = "secs"))
        eta_seconds <- if (completed_calibration_simulations > 0) {
          elapsed_seconds / completed_calibration_simulations *
            (total_calibration_simulations - completed_calibration_simulations)
        } else {
          NA_real_
        }
        poc_progress_log(
          progress_label,
          "still running: c_poc ", c_poc, " (", i, "/", length(c_poc_candidates), "), ",
          "simulation ", sim, "/", n_simulations, "; ",
          completed_calibration_simulations, "/", total_calibration_simulations,
          " calibration simulations done; elapsed ", poc_format_duration(elapsed_seconds),
          "; ETA ", poc_format_duration(eta_seconds),
          enabled = TRUE
        )
        next_progress_at <- now + progress_interval_seconds
      }
    }

    # Calculate metrics
    # PoC detection rate across ALL simulations (including early terminated = FALSE)
    poc_detection_rate <- poc_detection_count / n_simulations
    early_termination_rate <- early_termination_count / n_simulations
    completion_rate <- 1 - early_termination_rate

    # Among completed trials, what fraction had PoC validated?
    poc_rate_among_completed <- if (n_completed > 0) poc_detection_count / n_completed else 0

    # Monte Carlo standard error for PoC detection rate
    # SE = sqrt(p * (1-p) / N)
    poc_se <- sqrt(poc_detection_rate * (1 - poc_detection_rate) / n_simulations)
    poc_ci_lower <- max(0, poc_detection_rate - 1.96 * poc_se)
    poc_ci_upper <- min(1, poc_detection_rate + 1.96 * poc_se)

    calibration_results[[i]] <- list(
      c_poc = c_poc,
      poc_detection_rate = poc_detection_rate,
      poc_se = poc_se,
      poc_ci_lower = poc_ci_lower,
      poc_ci_upper = poc_ci_upper,
      early_termination_rate = early_termination_rate,
      completion_rate = completion_rate,
      poc_rate_among_completed = poc_rate_among_completed,
      n_simulations = n_simulations,
      n_completed = n_completed,
      early_termination_count = early_termination_count,
      simulation_results_stored = isTRUE(store_simulation_results),
      common_random_numbers = common_random_numbers,
      calibration_seed = calibration_seed,
      simulation_results = if (isTRUE(store_simulation_results)) simulation_results else list()
    )

    cat("  PoC detection rate:", round(poc_detection_rate, 3),
        "(SE:", round(poc_se, 4), ", 95% CI: [",
        round(poc_ci_lower, 3), ",", round(poc_ci_upper, 3), "])\n")
    cat("  Early termination rate:", round(early_termination_rate, 3), "\n")
    cat("  Completion rate:", round(completion_rate, 3), "\n")
    cat("  PoC rate among completed trials:", round(poc_rate_among_completed, 3), "\n")

    if (isTRUE(progress)) {
      completed_calibration_simulations <- i * n_simulations
      elapsed_seconds <- as.numeric(difftime(Sys.time(), calibration_started_at, units = "secs"))
      eta_seconds <- elapsed_seconds / completed_calibration_simulations *
        (total_calibration_simulations - completed_calibration_simulations)
      poc_progress_log(
        progress_label,
        "finished c_poc ", c_poc, " (", i, "/", length(c_poc_candidates), "); ",
        "PoC rate ", round(poc_detection_rate, 3),
        ", completion ", round(completion_rate, 3),
        "; ", completed_calibration_simulations, "/", total_calibration_simulations,
        " calibration simulations done; elapsed ", poc_format_duration(elapsed_seconds),
        "; ETA ", poc_format_duration(eta_seconds),
        enabled = TRUE
      )
    }

    # Sanity check: PoC rate should not exceed completion rate
    if (poc_detection_rate > completion_rate + 1e-6) {
      cat("  [WARNING] PoC rate (", round(poc_detection_rate, 3),
          ") exceeds completion rate (", round(completion_rate, 3), ") - this should not happen!\n")
    }
  }

  # Find optimal C_poc using proper Type I error control criterion
  # RULE: Choose smallest C_poc such that PoC detection rate ≤ target
  # If none satisfy, choose largest C_poc and report control not achieved
  poc_rates <- sapply(calibration_results, function(x) x$poc_detection_rate)
  completion_rates <- sapply(calibration_results, function(x) x$completion_rate)
  # Find all C_poc values that achieve Type I error control
  controlled_indices <- which(poc_rates <= target_rate)

  if (length(controlled_indices) > 0) {
    # Choose smallest C_poc that achieves control (most conservative)
    optimal_idx <- controlled_indices[1]
    optimal_c_poc <- c_poc_candidates[optimal_idx]
    control_achieved <- TRUE
  } else {
    # No C_poc achieves control → choose largest (most stringent)
    optimal_idx <- length(c_poc_candidates)
    optimal_c_poc <- c_poc_candidates[optimal_idx]
    control_achieved <- FALSE
  }

  cat("\n=== CALIBRATION RESULTS SUMMARY ===\n")
  cat("Target Type I error (PoC detection rate): ≤", target_rate * 100, "%\n")

  if (control_achieved) {
    cat("✓ TYPE I ERROR CONTROL ACHIEVED\n")
    cat("Optimal C_poc:", optimal_c_poc,
        "(smallest C_poc achieving control)\n")
  } else {
    cat("✗ TYPE I ERROR CONTROL NOT ACHIEVED\n")
    cat("Selected C_poc:", optimal_c_poc,
        "(largest tested, but still exceeds target)\n")
    cat("⚠ WARNING: All tested C_poc values exceed target Type I error!\n")
    cat("⚠ Consider testing higher C_poc values (e.g., 0.96, 0.97, 0.98)\n")
  }

  cat("Achieved PoC detection rate:", round(poc_rates[optimal_idx], 3),
      "(SE:", round(calibration_results[[optimal_idx]]$poc_se, 4), ")\n")
  cat("Achieved completion rate:", round(completion_rates[optimal_idx], 3), "\n")

  # Sanity checks and detailed results
  cat("\n=== TYPE I ERROR CONTROL TABLE ===\n")
  cat(sprintf("%-10s | %-20s | %-25s | %-15s\n",
              "C_poc", "PoC Rate (95% CI)", "Control Achieved?", "Status"))
  cat("--------------------------------------------------------------------------------\n")

  for (i in seq_along(calibration_results)) {
    res <- calibration_results[[i]]
    ci_str <- sprintf("[%.3f, %.3f]", res$poc_ci_lower, res$poc_ci_upper)
    control_ok <- res$poc_detection_rate <= target_rate
    control_str <- if (control_ok) "✓ Yes" else "✗ No"

    status_str <- ""
    if (i == optimal_idx) {
      status_str <- if (control_achieved) "← OPTIMAL" else "← SELECTED*"
    }

    cat(sprintf("%-10.4f | %-7.3f %-12s | %-25s | %-15s\n",
                res$c_poc,
                res$poc_detection_rate,
                ci_str,
                control_str,
                status_str))
  }
  cat("--------------------------------------------------------------------------------\n")
  if (!control_achieved) {
    cat("* Selected despite exceeding target (no C_poc achieved control)\n")
  }

  cat("\n=== SANITY CHECKS ===\n")
  cat("1. PoC rate <= completion rate for all c_poc values?\n")
  for (i in seq_along(calibration_results)) {
    res <- calibration_results[[i]]
    check_passed <- res$poc_detection_rate <= res$completion_rate + 1e-6
    status <- if (check_passed) "✓ PASS" else "✗ FAIL"
    cat("   c_poc =", res$c_poc, ":", status,
        "(PoC:", round(res$poc_detection_rate, 3),
        ", Completion:", round(res$completion_rate, 3), ")\n")
  }

  cat("\n2. Monotonicity check (higher c_poc should generally decrease PoC rate):\n")
  if (length(calibration_results) > 1) {
    for (i in 2:length(calibration_results)) {
      prev_rate <- calibration_results[[i-1]]$poc_detection_rate
      curr_rate <- calibration_results[[i]]$poc_detection_rate
      change <- curr_rate - prev_rate
      direction <- if (change <= 0) "✓ Non-increasing" else "⚠ Increased"
      cat("   c_poc:", calibration_results[[i-1]]$c_poc, "→", calibration_results[[i]]$c_poc,
          ":", direction, "(Δ =", round(change, 4), ")\n")
    }
  } else {
    cat("   Skipped: only one c_poc value tested.\n")
  }
  cat("\n=== END SANITY CHECKS ===\n")

  return(list(
    calibration_results = calibration_results,
    optimal_c_poc = optimal_c_poc,
    target_rate = target_rate,
    achieved_rate = poc_rates[optimal_idx],
    optimal_rate = poc_rates[optimal_idx],
    control_achieved = control_achieved,
    c_poc_candidates = c_poc_candidates,
    common_random_numbers = common_random_numbers,
    calibration_seed = calibration_seed,
    poc_detection_rates = poc_rates,
    n_simulations = n_simulations
  ))
}

poc_calibration_results_table <- function(calibration_results) {
  # Convert either supported PoC calibration result shape into a data frame.
  if (is.data.frame(calibration_results)) {
    return(calibration_results)
  }

  if (is.list(calibration_results) && "calibration_results" %in% names(calibration_results)) {
    raw_results <- calibration_results$calibration_results
  } else {
    raw_results <- calibration_results
  }

  if (is.data.frame(raw_results)) {
    return(raw_results)
  }

  if (!is.list(raw_results) || length(raw_results) == 0) {
    stop("calibration_results must contain a data frame or a non-empty list of candidate results.")
  }

  value_or <- function(x, name, default = NA_real_) {
    if (!is.null(x[[name]])) x[[name]] else default
  }

  rows <- lapply(raw_results, function(res) {
    data.frame(
      c_poc = value_or(res, "c_poc"),
      poc_detection_rate = value_or(res, "poc_detection_rate"),
      poc_detection_rate_lower = value_or(res, "poc_ci_lower"),
      poc_detection_rate_upper = value_or(res, "poc_ci_upper"),
      poc_se = value_or(res, "poc_se"),
      completion_rate = value_or(res, "completion_rate"),
      poc_rate_among_completed = value_or(res, "poc_rate_among_completed"),
      n_simulations = value_or(
        res,
        "n_simulations",
        length(value_or(res, "simulation_results", list()))
      ),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

create_default_null_scenario <- function(base_config = NULL) {
  n_doses <- if (!is.null(base_config$dose_levels)) length(base_config$dose_levels) else 5
  phi_I <- if (!is.null(base_config$phi_I)) base_config$phi_I else 0.20
  phi_E <- if (!is.null(base_config$phi_E)) base_config$phi_E else 0.25
  tox_flat <- if (!is.null(base_config$toxicity_low)) base_config$toxicity_low else 0.05

  create_null_flat_scenario(
    n_doses = n_doses,
    phi_I = phi_I,
    phi_E = phi_E,
    tox_flat = tox_flat
  )
}

run_calibration_simulation <- function(config, scenario_type = "flat_null", n_simulations = 1, seed = 123) {
  # Backward-compatible single-replicate helper built on the canonical simulator.
  if (scenario_type != "flat_null") {
    stop("Unsupported scenario type: ", scenario_type)
  }
  if (n_simulations != 1) {
    warning("run_calibration_simulation returns one logical result; use calibrate_c_poc() for repeated simulations.")
  }

  result <- run_single_calibration_simulation(
    config = within(config, {
      verbose_logging <- FALSE
      log_early_termination <- FALSE
    }),
    scenario_params = create_default_null_scenario(config),
    seed = seed
  )

  !result$metrics$terminated_early && isTRUE(result$metrics$poc_validated)
}

run_quick_calibration <- function(
  target_rate = 0.10,
  n_simulations = 100,
  c_poc_candidates = c(0.7, 0.8, 0.9, 0.95),
  base_config = NULL,
  verbose = TRUE
) {
  # Convenience wrapper for smoke tests and demos.
  calibrate_c_poc(
    null_scenario = create_default_null_scenario(base_config),
    c_poc_candidates = c_poc_candidates,
    n_simulations = n_simulations,
    base_config = base_config,
    target_rate = target_rate,
    debug_early_termination = FALSE,
    verbose = verbose
  )
}

validate_calibration <- function(
  calibration_results,
  n_validation_simulations = 1000,
  null_scenario = NULL,
  base_config = NULL
) {
  if (is.null(null_scenario)) {
    null_scenario <- create_default_null_scenario(base_config)
  }
  if (is.null(base_config)) {
    base_config <- list(
      dose_levels = seq_along(null_scenario$p_YI),
      n_stages = 5,
      cohort_size = 15,
      phi_T = 0.30,
      c_T = 0.3,
      phi_E = 0.20,
      c_E = 0.3,
      phi_I = 0.20,
      c_I = 0.3,
      delta_poc = 0.8,
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
    base_config$utility_table <- utility_table
  }

  optimal_c_poc <- calibration_results$optimal_c_poc
  if (is.null(optimal_c_poc)) {
    stop("calibration_results must include optimal_c_poc.")
  }

  config <- base_config
  config$c_poc <- optimal_c_poc

  detections <- vapply(seq_len(n_validation_simulations), function(i) {
    result <- run_single_calibration_simulation(config, null_scenario, seed = 900000 + i)
    !result$metrics$terminated_early && isTRUE(result$metrics$poc_validated)
  }, logical(1))

  validation_rate <- mean(detections)
  validation_ci <- binom.test(sum(detections), n_validation_simulations)$conf.int

  list(
    optimal_c_poc = optimal_c_poc,
    validation_rate = validation_rate,
    validation_ci = validation_ci,
    target_rate = calibration_results$target_rate,
    n_validation_simulations = n_validation_simulations
  )
}

save_calibration_results <- function(calibration_results, file_path = "results/poc_calibration_results.RData") {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  save(calibration_results, file = file_path)
  invisible(file_path)
}

load_calibration_results <- function(file_path = "results/poc_calibration_results.RData") {
  if (!file.exists(file_path)) {
    return(NULL)
  }
  load(file_path)
  calibration_results
}

poc_format_percent <- function(x, digits = 1) {
  if (is.null(x) || length(x) == 0) {
    return("NA")
  }
  if (length(x) > 1) {
    return(paste(vapply(x, poc_format_percent, character(1), digits = digits), collapse = ", "))
  }
  if (is.na(x)) {
    return("NA")
  }
  sprintf(paste0("%.", digits, "f%%"), x * 100)
}

poc_format_value <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0) {
    return("NA")
  }
  if (length(x) > 1) {
    return(paste(vapply(x, poc_format_value, character(1), digits = digits), collapse = ", "))
  }
  if (is.logical(x)) {
    return(if (isTRUE(x)) "TRUE" else "FALSE")
  }
  if (is.numeric(x)) {
    if (is.na(x)) {
      return("NA")
    }
    return(format(round(x, digits), nsmall = min(digits, 3), trim = TRUE))
  }
  as.character(x)
}

poc_config_markdown_lines <- function(base_config, null_scenario) {
  config_fields <- c(
    "dose_levels", "n_stages", "cohort_size",
    "phi_T", "c_T", "phi_E", "c_E", "phi_I", "c_I",
    "delta_poc", "enable_early_termination"
  )

  config_lines <- c("Trial/config parameters:")
  for (field in config_fields) {
    if (!is.null(base_config[[field]])) {
      config_lines <- c(
        config_lines,
        paste0("- `", field, "`: ", poc_format_value(base_config[[field]]))
      )
    }
  }

  scenario_lines <- c(
    "Null scenario:",
    paste0("- description: ", if (!is.null(null_scenario$description)) null_scenario$description else "NA"),
    paste0("- `p_YI`: ", poc_format_value(null_scenario$p_YI)),
    paste0("- `p_YT_given_I[, 1]`: ", poc_format_value(null_scenario$p_YT_given_I[, 1])),
    paste0("- `p_YT_given_I[, 2]`: ", poc_format_value(null_scenario$p_YT_given_I[, 2])),
    paste0("- `p_YE_given_I[, 1]`: ", poc_format_value(null_scenario$p_YE_given_I[, 1])),
    paste0("- `p_YE_given_I[, 2]`: ", poc_format_value(null_scenario$p_YE_given_I[, 2])),
    paste0("- `rho0`: ", poc_format_value(null_scenario$rho0)),
    paste0("- `rho1`: ", poc_format_value(null_scenario$rho1))
  )

  c(config_lines, "", scenario_lines)
}

append_poc_calibration_log <- function(
  calibration_results,
  null_scenario,
  base_config,
  file_path = "results/notebook_calibration/poc_calibration_history.md",
  run_label = NULL,
  notes = NULL
) {
  if (missing(calibration_results) || is.null(calibration_results)) {
    stop("calibration_results is required.")
  }
  if (missing(null_scenario) || is.null(null_scenario)) {
    stop("null_scenario is required.")
  }
  if (missing(base_config) || is.null(base_config)) {
    stop("base_config is required.")
  }
  if (!is.character(file_path) || length(file_path) != 1 || nchar(file_path) == 0) {
    stop("file_path must be a non-empty string.")
  }

  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)

  result_table <- poc_calibration_results_table(calibration_results)
  for (column in c(
    "poc_detection_rate_lower", "poc_detection_rate_upper",
    "completion_rate", "poc_rate_among_completed"
  )) {
    if (!column %in% names(result_table)) {
      result_table[[column]] <- NA_real_
    }
  }
  result_table$control_achieved <- result_table$poc_detection_rate <= calibration_results$target_rate
  result_table$is_selected <- result_table$c_poc == calibration_results$optimal_c_poc

  status <- if (isTRUE(calibration_results$control_achieved)) {
    "CONTROL ACHIEVED"
  } else {
    "CONTROL NOT ACHIEVED"
  }

  run_title <- if (!is.null(run_label) && nchar(run_label) > 0) {
    paste0(" - ", run_label)
  } else {
    ""
  }

  header_lines <- character(0)
  if (!file.exists(file_path) || file.info(file_path)$size == 0) {
    header_lines <- c(
      "# PoC Calibration History",
      "",
      "This file is appended by PoC calibration runs. Each entry records the tested parameters and outcome so future runs can avoid repeating unhelpful grids.",
      ""
    )
  }

  summary_lines <- c(
    paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), run_title),
    "",
    paste0("- Status: ", status),
    paste0("- Target Type I error: ", poc_format_percent(calibration_results$target_rate)),
    paste0("- Selected `c_poc`: ", poc_format_value(calibration_results$optimal_c_poc, digits = 4)),
    paste0("- Achieved PoC detection rate: ", poc_format_percent(calibration_results$achieved_rate)),
    paste0("- Simulations per candidate: ", poc_format_value(calibration_results$n_simulations)),
    paste0(
      "- Seed strategy: ",
      if (isTRUE(calibration_results$common_random_numbers)) {
        "common random numbers across `c_poc` candidates"
      } else {
        "independent seeds across `c_poc` candidates"
      },
      "; calibration seed = ",
      if (is.null(calibration_results$calibration_seed)) "NULL" else poc_format_value(calibration_results$calibration_seed, digits = 0)
    )
  )

  if (!is.null(notes) && nchar(notes) > 0) {
    summary_lines <- c(summary_lines, paste0("- Notes: ", notes))
  }

  table_lines <- c(
    "",
    "| c_poc | PoC detection | 95% CI | Completion | PoC among completed | Meets target? | Selected? |",
    "|---:|---:|---:|---:|---:|:---:|:---:|"
  )

  for (i in seq_len(nrow(result_table))) {
    row <- result_table[i, , drop = FALSE]
    ci <- if (!is.na(row$poc_detection_rate_lower) && !is.na(row$poc_detection_rate_upper)) {
      paste0(
        poc_format_percent(row$poc_detection_rate_lower),
        " to ",
        poc_format_percent(row$poc_detection_rate_upper)
      )
    } else {
      "NA"
    }
    table_lines <- c(
      table_lines,
      paste(
        "",
        poc_format_value(row$c_poc, digits = 4),
        poc_format_percent(row$poc_detection_rate),
        ci,
        poc_format_percent(row$completion_rate),
        poc_format_percent(row$poc_rate_among_completed),
        if (isTRUE(row$control_achieved)) "yes" else "no",
        if (isTRUE(row$is_selected)) "yes" else "no",
        "",
        sep = " | "
      )
    )
  }

  recommendation_lines <- character(0)
  if (!isTRUE(calibration_results$control_achieved)) {
    recommendation_lines <- c(
      "",
      "Parameter note: no tested `c_poc` controlled the null PoC detection rate. Next runs should test higher `c_poc` values and/or tune `c_T`, `c_E`, and `c_I` while keeping the clinical scenario definition fixed."
    )
  }

  log_lines <- c(
    header_lines,
    summary_lines,
    "",
    poc_config_markdown_lines(base_config, null_scenario),
    table_lines,
    recommendation_lines,
    ""
  )

  write(log_lines, file = file_path, append = TRUE)
  invisible(file_path)
}

default_poc_parameter_search_grid <- function() {
  expand.grid(
    c_T = c(0.20, 0.35, 0.50, 0.65),
    c_E = c(0.50, 0.65, 0.80),
    c_I = c(0.35, 0.50, 0.65),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

append_poc_parameter_search_log <- function(
  search_results,
  file_path = "results/notebook_calibration/poc_parameter_search_history.md",
  run_label = NULL,
  notes = NULL
) {
  if (missing(search_results) || is.null(search_results$summary)) {
    stop("search_results must be returned by run_poc_parameter_search().")
  }
  if (!is.character(file_path) || length(file_path) != 1 || nchar(file_path) == 0) {
    stop("file_path must be a non-empty string.")
  }

  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)

  summary_table <- search_results$summary
  run_title <- if (!is.null(run_label) && nchar(run_label) > 0) {
    paste0(" - ", run_label)
  } else {
    ""
  }

  header_lines <- character(0)
  if (!file.exists(file_path) || file.info(file_path)$size == 0) {
    header_lines <- c(
      "# PoC Parameter Search History",
      "",
      "This file is appended by batch PoC parameter searches. It records the tested cutoff grids and ranked outcomes.",
      ""
    )
  }

  best <- summary_table[1, , drop = FALSE]
  summary_lines <- c(
    paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), run_title),
    "",
    paste0("- Tested variables: ", paste(search_results$tested_variables, collapse = ", ")),
    paste0("- Fixed variables: ", paste(search_results$fixed_parameters, collapse = ", ")),
    paste0(
      "- Candidate `c_poc` values: ",
      paste(vapply(search_results$c_poc_candidates, poc_format_value, character(1), digits = 4), collapse = ", ")
    ),
    paste0("- Simulations per `c_poc`: ", search_results$n_simulations),
    paste0("- Target Type I error: ", poc_format_percent(search_results$target_rate)),
    paste0(
      "- Seed strategy: ",
      if (isTRUE(search_results$common_random_numbers)) {
        "common random numbers within each parameter row"
      } else {
        "independent seeds across `c_poc` candidates"
      },
      "; base calibration seed = ",
      if (is.null(search_results$calibration_seed)) "NULL" else poc_format_value(search_results$calibration_seed, digits = 0)
    ),
    paste0(
      "- Best ranked setting: row ", best$search_id,
      ", selected `c_poc` = ", poc_format_value(best$optimal_c_poc, digits = 4),
      ", achieved rate = ", poc_format_percent(best$achieved_rate),
      ", control achieved = ", if (isTRUE(best$control_achieved)) "yes" else "no"
    )
  )

  if (!is.null(notes) && nchar(notes) > 0) {
    summary_lines <- c(summary_lines, paste0("- Notes: ", notes))
  }

  parameter_columns <- intersect(
    c("c_T", "c_E", "c_I", "delta_poc"),
    names(summary_table)
  )
  table_header <- paste(
    c("rank", "search_id", parameter_columns, "selected c_poc", "PoC detection", "completion", "control"),
    collapse = " | "
  )
  table_rule <- paste(rep("---", length(strsplit(table_header, " \\| ")[[1]])), collapse = " | ")
  table_lines <- c("", paste0("| ", table_header, " |"), paste0("| ", table_rule, " |"))

  for (i in seq_len(nrow(summary_table))) {
    row <- summary_table[i, , drop = FALSE]
    param_values <- vapply(parameter_columns, function(col) poc_format_value(row[[col]]), character(1))
    table_lines <- c(
      table_lines,
      paste0(
        "| ",
        paste(
          c(
            i,
            row$search_id,
            param_values,
            poc_format_value(row$optimal_c_poc, digits = 4),
            poc_format_percent(row$achieved_rate),
            poc_format_percent(row$completion_rate),
            if (isTRUE(row$control_achieved)) "yes" else "no"
          ),
          collapse = " | "
        ),
        " |"
      )
    )
  }

  write(c(header_lines, summary_lines, table_lines, ""), file = file_path, append = TRUE)
  invisible(file_path)
}

run_poc_parameter_search <- function(
  null_scenario,
  base_config,
  search_grid = default_poc_parameter_search_grid(),
  c_poc_candidates = c(0.90, 0.95, 0.97, 0.98, 0.99, 0.995),
  n_simulations = 100,
  target_rate = 0.10,
  log_file = "results/notebook_calibration/poc_parameter_search_history.md",
  append_log = TRUE,
  run_label = NULL,
  verbose = TRUE,
  progress = verbose,
  progress_interval_seconds = 300,
  common_random_numbers = TRUE,
  calibration_seed = 10000
) {
  if (missing(null_scenario) || is.null(null_scenario)) {
    stop("null_scenario is required.")
  }
  if (missing(base_config) || is.null(base_config)) {
    stop("base_config is required.")
  }
  if (is.list(search_grid) && !is.data.frame(search_grid)) {
    search_grid <- expand.grid(search_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(search_grid) || nrow(search_grid) == 0) {
    stop("search_grid must be a non-empty data frame or list.")
  }

  allowed_variables <- c("c_T", "c_E", "c_I", "delta_poc")
  invalid_variables <- setdiff(names(search_grid), allowed_variables)
  if (length(invalid_variables) > 0) {
    stop(
      "search_grid can only tune these variables: ",
      paste(allowed_variables, collapse = ", "),
      ". Invalid variables: ",
      paste(invalid_variables, collapse = ", ")
    )
  }
  if (!is.numeric(n_simulations) || length(n_simulations) != 1 || n_simulations < 1) {
    stop("n_simulations must be a positive integer.")
  }
  if (!is.numeric(progress_interval_seconds) || length(progress_interval_seconds) != 1 || progress_interval_seconds < 1) {
    stop("progress_interval_seconds must be a positive number of seconds.")
  }
  if (!is.logical(common_random_numbers) || length(common_random_numbers) != 1) {
    stop("common_random_numbers must be TRUE or FALSE.")
  }
  if (!is.null(calibration_seed) && (
    !is.numeric(calibration_seed) ||
      length(calibration_seed) != 1 ||
      is.na(calibration_seed) ||
      calibration_seed < 0
  )) {
    stop("calibration_seed must be NULL or a non-negative number.")
  }
  if (!is.null(calibration_seed)) {
    calibration_seed <- as.numeric(calibration_seed)
  }
  if (isTRUE(common_random_numbers) && is.null(calibration_seed)) {
    stop("calibration_seed cannot be NULL when common_random_numbers is TRUE.")
  }

  fixed_parameters <- setdiff(
    c(
      "dose_levels", "n_stages", "cohort_size",
      "phi_T", "phi_E", "phi_I",
      "c_T", "c_E", "c_I", "delta_poc",
      "utility_table", "null_scenario", "rho0", "rho1"
    ),
    names(search_grid)
  )

  details <- vector("list", nrow(search_grid))
  summary_rows <- vector("list", nrow(search_grid))
  search_started_at <- Sys.time()
  simulations_per_row <- length(c_poc_candidates) * n_simulations
  total_search_simulations <- nrow(search_grid) * simulations_per_row
  row_seed_stride <- 10000000

  if (isTRUE(progress)) {
    poc_progress_log(
      "PoC parameter search workload: ",
      nrow(search_grid), " parameter rows × ",
      length(c_poc_candidates), " c_poc candidates × ",
      n_simulations, " simulations = ",
      total_search_simulations, " trial simulations. ",
      if (isTRUE(common_random_numbers)) {
        "Using common random numbers within each parameter row. "
      } else {
        "Using independent seeds across C_poc candidates. "
      },
      "Progress will be printed at row/candidate boundaries and about every ",
      poc_format_duration(progress_interval_seconds), " during long nodes.",
      enabled = TRUE
    )
  }

  for (i in seq_len(nrow(search_grid))) {
    row_values <- as.list(search_grid[i, , drop = FALSE])
    config <- base_config
    for (name in names(row_values)) {
      config[[name]] <- row_values[[name]]
    }

    if (isTRUE(verbose)) {
      cat("\n=== PoC parameter search row", i, "of", nrow(search_grid), "===\n")
      cat(
        paste(
          paste0(names(row_values), "=", unlist(row_values, use.names = FALSE)),
          collapse = ", "
        ),
        "\n"
      )
    }
    if (isTRUE(progress)) {
      elapsed_seconds <- as.numeric(difftime(Sys.time(), search_started_at, units = "secs"))
      completed_search_simulations <- (i - 1) * simulations_per_row
      eta_seconds <- if (completed_search_simulations > 0) {
        elapsed_seconds / completed_search_simulations *
          (total_search_simulations - completed_search_simulations)
      } else {
        NA_real_
      }
      poc_progress_log(
        "Starting search row ", i, "/", nrow(search_grid), " (",
        paste(paste0(names(row_values), "=", unlist(row_values, use.names = FALSE)), collapse = ", "),
        "); row workload ", simulations_per_row, " trial simulations; overall ",
        completed_search_simulations, "/", total_search_simulations,
        " done; elapsed ", poc_format_duration(elapsed_seconds),
        "; ETA ", poc_format_duration(eta_seconds),
        enabled = TRUE
      )
    }

    row_calibration_seed <- if (is.null(calibration_seed)) {
      NULL
    } else {
      calibration_seed + (i - 1) * row_seed_stride
    }

    calibration <- calibrate_c_poc(
      null_scenario = null_scenario,
      c_poc_candidates = c_poc_candidates,
      n_simulations = n_simulations,
      base_config = config,
      target_rate = target_rate,
      debug_early_termination = FALSE,
      verbose = FALSE,
      progress = progress,
      progress_interval_seconds = progress_interval_seconds,
      progress_prefix = paste0("[search row ", i, "/", nrow(search_grid), "]"),
      store_simulation_results = FALSE,
      common_random_numbers = common_random_numbers,
      calibration_seed = row_calibration_seed
    )

    result_table <- poc_calibration_results_table(calibration)
    selected_idx <- which(result_table$c_poc == calibration$optimal_c_poc)[1]
    selected_row <- result_table[selected_idx, , drop = FALSE]

    summary_row <- data.frame(
      search_id = i,
      optimal_c_poc = calibration$optimal_c_poc,
      achieved_rate = calibration$achieved_rate,
      target_rate = calibration$target_rate,
      control_achieved = calibration$control_achieved,
      completion_rate = selected_row$completion_rate,
      poc_rate_among_completed = selected_row$poc_rate_among_completed,
      n_simulations = n_simulations,
      stringsAsFactors = FALSE
    )
    summary_row <- cbind(search_grid[i, , drop = FALSE], summary_row)
    summary_row$ranking_gap <- if (isTRUE(calibration$control_achieved)) {
      target_rate - calibration$achieved_rate
    } else {
      calibration$achieved_rate - target_rate
    }

    details[[i]] <- list(
      search_id = i,
      tested_parameters = row_values,
      config = config,
      calibration_results = calibration,
      candidate_table = result_table,
      calibration_seed = row_calibration_seed
    )
    summary_rows[[i]] <- summary_row

    if (isTRUE(progress)) {
      elapsed_seconds <- as.numeric(difftime(Sys.time(), search_started_at, units = "secs"))
      completed_search_simulations <- i * simulations_per_row
      eta_seconds <- elapsed_seconds / completed_search_simulations *
        (total_search_simulations - completed_search_simulations)
      poc_progress_log(
        "Finished search row ", i, "/", nrow(search_grid), "; selected c_poc ",
        calibration$optimal_c_poc, "; achieved rate ", round(calibration$achieved_rate, 3),
        "; completion ", round(selected_row$completion_rate, 3),
        "; control ", if (isTRUE(calibration$control_achieved)) "yes" else "no",
        "; overall ", completed_search_simulations, "/", total_search_simulations,
        " trial simulations done; elapsed ", poc_format_duration(elapsed_seconds),
        "; ETA ", poc_format_duration(eta_seconds),
        enabled = TRUE
      )
    }
  }

  summary_table <- do.call(rbind, summary_rows)
  summary_table <- summary_table[order(
    !summary_table$control_achieved,
    summary_table$ranking_gap,
    -summary_table$completion_rate
  ), , drop = FALSE]
  rownames(summary_table) <- NULL

  search_results <- list(
    summary = summary_table,
    details = details,
    tested_variables = names(search_grid),
    fixed_parameters = fixed_parameters,
    c_poc_candidates = c_poc_candidates,
    n_simulations = n_simulations,
    target_rate = target_rate,
    common_random_numbers = common_random_numbers,
    calibration_seed = calibration_seed,
    row_seed_stride = row_seed_stride
  )

  if (isTRUE(append_log)) {
    append_poc_parameter_search_log(
      search_results = search_results,
      file_path = log_file,
      run_label = run_label
    )
  }

  if (isTRUE(progress)) {
    poc_progress_log(
      "PoC parameter search complete; total elapsed ",
      poc_format_duration(as.numeric(difftime(Sys.time(), search_started_at, units = "secs"))),
      ". Best row: ", search_results$summary$search_id[1],
      ", selected c_poc ", search_results$summary$optimal_c_poc[1],
      ", achieved rate ", round(search_results$summary$achieved_rate[1], 3),
      ", control ", if (isTRUE(search_results$summary$control_achieved[1])) "yes" else "no",
      enabled = TRUE
    )
  }

  search_results
}

# Create calibration curve plot
plot_calibration_curve <- function(calibration_results, file_path = NULL) {
  # Extract data for plotting
  c_poc_values <- sapply(calibration_results$calibration_results, function(x) x$c_poc)
  poc_rates <- sapply(calibration_results$calibration_results, function(x) x$poc_detection_rate)
  poc_se <- sapply(calibration_results$calibration_results, function(x) x$poc_se)
  completion_rates <- sapply(calibration_results$calibration_results, function(x) x$completion_rate)

  plot_data <- data.frame(
    c_poc = c_poc_values,
    poc_detection_rate = poc_rates,
    poc_se = poc_se,
    completion_rate = completion_rates,
    poc_ci_lower = pmax(0, poc_rates - 1.96 * poc_se),
    poc_ci_upper = pmin(1, poc_rates + 1.96 * poc_se)
  )

  p <- ggplot(plot_data, aes(x = c_poc, y = poc_detection_rate)) +
    geom_ribbon(aes(ymin = poc_ci_lower, ymax = poc_ci_upper),
                fill = "#2E86AB", alpha = 0.2) +
    geom_line(aes(y = completion_rate, color = "Completion Rate"),
              linewidth = 1, linetype = "dotted") +
    geom_line(aes(color = "PoC Detection Rate"), linewidth = 1.2) +
    geom_point(aes(color = "PoC Detection Rate"), size = 3) +
    geom_hline(yintercept = calibration_results$target_rate, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = calibration_results$optimal_c_poc,
               linetype = "dashed", color = "blue", linewidth = 1) +
    scale_color_manual(values = c("PoC Detection Rate" = "#2E86AB",
                                  "Completion Rate" = "#7209B7")) +
    labs(
      title = "PoC Calibration Curve (Null/Flat Scenario)",
      subtitle = paste0(
        if (calibration_results$control_achieved) {
          paste0("✓ Optimal C_poc = ", calibration_results$optimal_c_poc,
                 " (smallest achieving ≤10% Type I error)")
        } else {
          paste0("⚠ Selected C_poc = ", calibration_results$optimal_c_poc,
                 " (control NOT achieved - consider higher values)")
        },
        "\nAchieved: ", round(calibration_results$achieved_rate * 100, 1),
        "% (SE: ", round(calibration_results$calibration_results[[which(
          sapply(calibration_results$calibration_results, function(x) x$c_poc) ==
            calibration_results$optimal_c_poc
        )]]$poc_se * 100, 1), "%)"
      ),
      x = "C_poc Threshold",
      y = "Rate",
      color = "",
      caption = "Red line: Target 10% Type I error rate | Blue line: Selected C_poc\nShaded area: 95% CI (Monte Carlo SE) | Dotted line: Trial completion rate"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    scale_y_continuous(labels = scales::percent_format(),
                      limits = c(0, max(c(poc_rates + 2*poc_se, completion_rates)) * 1.1)) +
    scale_x_continuous(breaks = c_poc_values)

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 10, height = 7, dpi = 300)
    cat("Calibration curve saved to:", file_path, "\n")
  }

  return(p)
}

# Generate detailed calibration report
generate_calibration_report <- function(calibration_results, null_scenario, base_config, file_path = "calibration_report.txt") {
  # Create directory if it doesn't exist
  report_dir <- dirname(file_path)
  if (!dir.exists(report_dir)) {
    dir.create(report_dir, recursive = TRUE)
    cat("Created directory:", report_dir, "\n")
  }

  # Open file connection
  sink(file_path)

  cat("================================================================================\n")
  cat("                    POC CALIBRATION DETAILED REPORT                            \n")
  cat("================================================================================\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  # Section 1: Configuration Summary
  cat("================================================================================\n")
  cat("1. CONFIGURATION SUMMARY\n")
  cat("================================================================================\n\n")

  cat("Trial Design:\n")
  cat("  - Number of doses:", length(base_config$dose_levels), "\n")
  cat("  - Number of stages:", base_config$n_stages, "\n")
  cat("  - Cohort size:", base_config$cohort_size, "\n")
  cat("  - Total sample size (if completed):", base_config$n_stages * base_config$cohort_size, "\n\n")

  cat("Admissibility Thresholds:\n")
  cat("  - Toxicity (φ_T):", base_config$phi_T, "  (Probability threshold c_T:", base_config$c_T, ")\n")
  cat("  - Efficacy (φ_E):", base_config$phi_E, "  (Probability threshold c_E:", base_config$c_E, ")\n")
  cat("  - Immune Response (φ_I):", base_config$phi_I, "  (Probability threshold c_I:", base_config$c_I, ")\n\n")

  cat("PoC Parameters:\n")
  cat("  - delta_poc:", base_config$delta_poc, "\n")
  cat("  - Early termination enabled:", base_config$enable_early_termination, "\n\n")

  cat("Null Scenario Parameters:\n")
  cat("  - True immune response (p_YI):", paste(round(null_scenario$p_YI, 3), collapse=", "), "\n")
  cat("  - True toxicity (p_YT | I=0):", paste(round(null_scenario$p_YT_given_I[,1], 3), collapse=", "), "\n")
  cat("  - True toxicity (p_YT | I=1):", paste(round(null_scenario$p_YT_given_I[,2], 3), collapse=", "), "\n")
  cat("  - True efficacy (p_YE | I=0):", paste(round(null_scenario$p_YE_given_I[,1], 3), collapse=", "), "\n")
  cat("  - True efficacy (p_YE | I=1):", paste(round(null_scenario$p_YE_given_I[,2], 3), collapse=", "), "\n\n")

  # Section 2: Calibration Results Summary
  cat("================================================================================\n")
  cat("2. CALIBRATION RESULTS SUMMARY\n")
  cat("================================================================================\n\n")

  cat("Target Type I Error Rate: 10%\n")
  cat("Optimal C_poc:", calibration_results$optimal_c_poc, "\n")
  cat("Achieved PoC Detection Rate:", sprintf("%.1f%%", calibration_results$achieved_rate * 100), "\n\n")

  # Results table
  cat("Detailed Results by C_poc:\n")
  cat("--------------------------------------------------------------------------------------------------------\n")
  cat(sprintf("%-8s | %-22s | %-20s | %-20s | %-12s\n",
              "C_poc", "PoC Detection Rate", "Completion Rate", "PoC | Completed", "Status"))
  cat("--------------------------------------------------------------------------------------------------------\n")

  for (i in seq_along(calibration_results$calibration_results)) {
    res <- calibration_results$calibration_results[[i]]
    status <- if (res$c_poc == calibration_results$optimal_c_poc) "OPTIMAL" else ""

    # Format with MC standard error
    poc_str <- sprintf("%.1f%% ± %.1f%%",
                       res$poc_detection_rate * 100,
                       res$poc_se * 100)
    completion_str <- sprintf("%.1f%%", res$completion_rate * 100)
    poc_completed_str <- sprintf("%.1f%%", res$poc_rate_among_completed * 100)

    cat(sprintf("%-8.2f | %-22s | %-20s | %-20s | %-12s\n",
                res$c_poc,
                poc_str,
                completion_str,
                poc_completed_str,
                status))
  }
  cat("--------------------------------------------------------------------------------------------------------\n")
  cat("Note: PoC Detection Rate = Pr(PoC detected) across ALL trials (including early terminated)\n")
  cat("      PoC | Completed = Pr(PoC detected | trial completed without early termination)\n")
  cat("      Standard errors are Monte Carlo SE = sqrt(p*(1-p)/N)\n\n")

  # Section 3: Early Termination Analysis
  cat("================================================================================\n")
  cat("3. EARLY TERMINATION DETAILED ANALYSIS\n")
  cat("================================================================================\n\n")

  for (i in seq_along(calibration_results$calibration_results)) {
    res <- calibration_results$calibration_results[[i]]
    cat("--- C_poc =", res$c_poc, "---\n\n")

    n_sims <- length(res$simulation_results)
    n_early_term <- sum(sapply(res$simulation_results, function(x) x$metrics$terminated_early))
    n_completed <- n_sims - n_early_term

    cat("Overall Statistics:\n")
    cat("  - Total simulations:", n_sims, "\n")
    cat("  - Early terminations:", n_early_term, sprintf("(%.1f%%)", n_early_term/n_sims*100), "\n")
    cat("  - Completed trials:", n_completed, sprintf("(%.1f%%)", n_completed/n_sims*100), "\n")

    if (n_early_term > 0) {
      # Termination stage distribution
      term_stages <- sapply(res$simulation_results, function(x) {
        if (x$metrics$terminated_early) x$metrics$termination_stage else NA
      })
      term_stages <- term_stages[!is.na(term_stages)]

      cat("\n  Early Termination by Stage:\n")
      stage_table <- table(term_stages)
      for (stage in sort(unique(term_stages))) {
        count <- stage_table[as.character(stage)]
        cat(sprintf("    Stage %d: %d trials (%.1f%% of early terminations)\n",
                    stage, count, count/n_early_term*100))
      }

      # Sample size distribution for early terminations
      early_term_samples <- sapply(res$simulation_results, function(x) {
        if (x$metrics$terminated_early) x$metrics$total_participants else NA
      })
      early_term_samples <- early_term_samples[!is.na(early_term_samples)]

      cat("\n  Sample Size at Early Termination:\n")
      cat(sprintf("    Mean: %.1f patients (SD: %.1f)\n", mean(early_term_samples), sd(early_term_samples)))
      cat(sprintf("    Range: %d - %d patients\n", min(early_term_samples), max(early_term_samples)))

      # Analyze reasons for early termination
      cat("\n  Common Patterns in Early Termination:\n")

      # Get a few example cases with posterior summaries
      example_count <- 0
      max_examples <- 3

      for (sim_idx in seq_along(res$simulation_results)) {
        sim_res <- res$simulation_results[[sim_idx]]
        if (sim_res$metrics$terminated_early && example_count < max_examples) {
          example_count <- example_count + 1
          ps <- sim_res$debug_info$posterior_summaries

          if (!is.null(ps)) {
            cat(sprintf("\n    Example %d (Simulation %d, Stage %d):\n",
                        example_count, sim_idx, sim_res$metrics$termination_stage))

            # Posterior means
            cat("      Posterior Estimates (Mean):\n")
            cat("        Toxicity:  ", paste(sprintf("%.3f", ps$tox_marginal$marginal_prob), collapse=", "), "\n")
            cat("        Efficacy:  ", paste(sprintf("%.3f", ps$eff_marginal$marginal_prob), collapse=", "), "\n")
            cat("        Immune:    ", paste(sprintf("%.3f", ps$imm$pava_mean), collapse=", "), "\n")

            # Admissibility probabilities
            cat("      Admissibility Probabilities:\n")
            for (dose in seq_along(base_config$dose_levels)) {
              tox_prob <- mean(ps$tox_marginal$samples[[dose]] < base_config$phi_T)
              eff_prob <- mean(ps$eff_marginal$samples[[dose]] > base_config$phi_E)
              imm_prob <- mean(ps$imm$samples_pava[[dose]] > base_config$phi_I)

              admissible <- (tox_prob >= base_config$c_T) &&
                           (eff_prob >= base_config$c_E) &&
                           (imm_prob >= base_config$c_I)

              status_mark <- if (admissible) "✓" else "✗"
              cat(sprintf("        Dose %d %s: P(Tox<%.2f)=%.2f, P(Eff>%.2f)=%.2f, P(Imm>%.2f)=%.2f\n",
                          dose, status_mark,
                          base_config$phi_T, tox_prob,
                          base_config$phi_E, eff_prob,
                          base_config$phi_I, imm_prob))
            }
          }
        }
      }
    }

    # Completed trials analysis
    if (n_completed > 0) {
      cat("\n  Completed Trials Analysis:\n")

      completed_samples <- sapply(res$simulation_results, function(x) {
        if (!x$metrics$terminated_early) x$metrics$total_participants else NA
      })
      completed_samples <- completed_samples[!is.na(completed_samples)]

      cat(sprintf("    Mean sample size: %.1f patients (SD: %.1f)\n",
                  mean(completed_samples), sd(completed_samples)))

      poc_validated_count <- sum(sapply(res$simulation_results, function(x) {
        !x$metrics$terminated_early && x$metrics$poc_validated
      }))

      cat(sprintf("    PoC validated: %d trials (%.1f%% of completed trials)\n",
                  poc_validated_count, poc_validated_count/n_completed*100))
    }

    cat("\n")
  }

  # Section 4: Recommendations
  cat("================================================================================\n")
  cat("4. RECOMMENDATIONS\n")
  cat("================================================================================\n\n")

  optimal_result <- calibration_results$calibration_results[[which(
    sapply(calibration_results$calibration_results, function(x) x$c_poc) ==
      calibration_results$optimal_c_poc
  )]]

  cat("Based on the calibration results:\n\n")

  if (calibration_results$control_achieved) {
    cat("1. ✓ Type I Error Control Status: ACHIEVED\n")
    cat("   Selected C_poc =", calibration_results$optimal_c_poc, "\n")
    cat("   Selection criterion: Smallest C_poc achieving ≤10% Type I error\n")
    cat("   Achieved Type I error rate:", sprintf("%.1f%%", calibration_results$achieved_rate * 100),
        "(target: ≤10%)\n\n")
  } else {
    cat("1. ✗ Type I Error Control Status: NOT ACHIEVED\n")
    cat("   Selected C_poc =", calibration_results$optimal_c_poc, "(largest tested)\n")
    cat("   ⚠ WARNING: This still exceeds the 10% target!\n")
    cat("   Achieved Type I error rate:", sprintf("%.1f%%", calibration_results$achieved_rate * 100),
        "(target: ≤10%)\n")
    cat("   RECOMMENDATION: Test higher C_poc values (e.g., 0.96, 0.97, 0.98, 0.99)\n\n")
  }

  cat("2. Expected Trial Characteristics:\n")
  cat("   - Early termination rate:", sprintf("%.1f%%", optimal_result$early_termination_rate * 100), "\n")
  cat("   - Trial completion rate:", sprintf("%.1f%%", (1 - optimal_result$early_termination_rate) * 100), "\n\n")

  if (optimal_result$early_termination_rate > 0.7) {
    cat("⚠ WARNING: High early termination rate (>70%)\n")
    cat("   Consider:\n")
    cat("   - Relaxing admissibility thresholds (c_T, c_E, c_I)\n")
    cat("   - Adjusting null scenario parameters\n")
    cat("   - Increasing cohort size for better estimation\n\n")
  }

  cat("3. Next Steps:\n")
  cat("   - Validate calibration with alternative scenarios\n")
  cat("   - Test performance in signal scenarios\n")
  cat("   - Consider sensitivity analysis for key parameters\n\n")

  cat("================================================================================\n")
  cat("                           END OF REPORT                                       \n")
  cat("================================================================================\n")

  # Close file connection
  sink()

  cat("\nDetailed calibration report saved to:", file_path, "\n")

  return(invisible(file_path))
}
