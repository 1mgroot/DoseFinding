# PoC Calibration System - New Implementation
# This script implements the null/flat scenario calibration methodology per email specification

library(dplyr)
library(ggplot2)

# Source required functions
# NOTE: config.R should be sourced only once at the top level (e.g., notebook or main script)
source("src/utils/helpers.R")
source("src/core/simulate_data.R")
source("src/core/model_utils.R")
source("src/decision/dose_decision.R")
source("src/core/main.R")

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
      config, 
      scenario_params$p_YI, 
      scenario_params$p_YT_given_I, 
      scenario_params$p_YE_given_I, 
      scenario_params$rho0, 
      scenario_params$rho1,
      seed = seed
    )
    
    # Extract key metrics
    metrics <- list(
      terminated_early = results$terminated_early,
      termination_stage = ifelse(results$terminated_early, results$termination_stage, NA),
      final_od = ifelse(results$terminated_early, NA, results$final_od),
      poc_validated = ifelse(results$terminated_early, FALSE, results$poc_validated),
      poc_probability = ifelse(results$terminated_early, 0, results$poc_probability),
      total_participants = nrow(results$all_data),
      true_optimal_selected = ifelse(
        results$terminated_early, 
        FALSE, 
        results$final_od == scenario_params$true_optimal_dose
      )
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
calibrate_c_poc <- function(
  null_scenario,
  c_poc_candidates = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
  n_simulations = 1000,
  base_config = NULL,
  debug_early_termination = FALSE,  # Changed default to FALSE to reduce console output
  max_debug_cases_per_candidate = 2  # Reduced from 3 to 2
) {
  cat("Starting C_poc calibration...\n")
  cat("Testing", length(c_poc_candidates), "C_poc values\n")
  cat("Simulations per value:", n_simulations, "\n")
  
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
      log_early_termination = FALSE
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
  
  for (i in seq_along(c_poc_candidates)) {
    c_poc <- c_poc_candidates[i]
    cat("=== Testing C_poc =", c_poc, "===\n")
    
    # Update config with current C_poc
    config <- base_config
    config$c_poc <- c_poc
    
    # Run simulations
    simulation_results <- list()
    poc_detection_count <- 0
    debug_count <- 0
    
    # Use a large base seed (10000 * c_poc_index) to ensure different seeds for different c_poc values
    base_seed <- 10000 * i
    
    for (sim in 1:n_simulations) {
      # Each simulation gets a unique seed: base_seed + sim
      result <- run_single_calibration_simulation(config, null_scenario, seed = base_seed + sim)
      simulation_results[[sim]] <- result
      
      # Count PoC detection: ONLY when trial completes AND P_final is non-empty
      # Early terminated trials MUST have poc_validated = FALSE
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
    }
    
    # Calculate metrics
    # PoC detection rate across ALL simulations (including early terminated = FALSE)
    poc_detection_rate <- poc_detection_count / n_simulations
    early_termination_rate <- mean(sapply(simulation_results, function(x) x$metrics$terminated_early))
    completion_rate <- 1 - early_termination_rate
    
    # Among completed trials, what fraction had PoC validated?
    n_completed <- sum(!sapply(simulation_results, function(x) x$metrics$terminated_early))
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
      simulation_results = simulation_results
    )
    
    cat("  PoC detection rate:", round(poc_detection_rate, 3), 
        "(SE:", round(poc_se, 4), ", 95% CI: [", 
        round(poc_ci_lower, 3), ",", round(poc_ci_upper, 3), "])\n")
    cat("  Early termination rate:", round(early_termination_rate, 3), "\n")
    cat("  Completion rate:", round(completion_rate, 3), "\n")
    cat("  PoC rate among completed trials:", round(poc_rate_among_completed, 3), "\n")
    
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
  target_rate <- 0.10
  
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
    
    cat(sprintf("%-10.2f | %-7.3f %-12s | %-25s | %-15s\n", 
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
  for (i in 2:length(calibration_results)) {
    prev_rate <- calibration_results[[i-1]]$poc_detection_rate
    curr_rate <- calibration_results[[i]]$poc_detection_rate
    change <- curr_rate - prev_rate
    direction <- if (change <= 0) "✓ Non-increasing" else "⚠ Increased"
    cat("   c_poc:", calibration_results[[i-1]]$c_poc, "→", calibration_results[[i]]$c_poc,
        ":", direction, "(Δ =", round(change, 4), ")\n")
  }
  cat("\n=== END SANITY CHECKS ===\n")
  
  return(list(
    calibration_results = calibration_results,
    optimal_c_poc = optimal_c_poc,
    target_rate = target_rate,
    achieved_rate = poc_rates[optimal_idx],
    control_achieved = control_achieved,
    c_poc_candidates = c_poc_candidates,
    poc_detection_rates = poc_rates
  ))
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
              size = 1, linetype = "dotted") +
    geom_line(aes(color = "PoC Detection Rate"), size = 1.2) +
    geom_point(aes(color = "PoC Detection Rate"), size = 3) +
    geom_hline(yintercept = 0.10, linetype = "dashed", color = "red", size = 1) +
    geom_vline(xintercept = calibration_results$optimal_c_poc, 
               linetype = "dashed", color = "blue", size = 1) +
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
