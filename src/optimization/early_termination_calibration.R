# Early Termination Calibration Framework
# This file implements the calibration framework for early termination parameters

# Load required libraries
library(dplyr)
library(ggplot2)

# Source required functions
# Determine project root based on current working directory
project_root <- if (basename(getwd()) == "notebooks") ".." else "."

source(file.path(project_root, "src/core/simulate_data.R"))
source(file.path(project_root, "src/core/config.R"))
source(file.path(project_root, "src/core/model_utils.R"))
source(file.path(project_root, "src/decision/dose_decision.R"))
source(file.path(project_root, "src/core/main.R"))

run_early_termination_simulation <- function(config, scenario_type = "unfavorable", n_simulations = 1000, seed = 123) {
  # Run a single early termination simulation to determine if trial terminates early.
  #
  # Args:
  #   config: Trial configuration with current threshold values
  #   scenario_type: Type of scenario ("unfavorable", "flat_null", etc.)
  #   n_simulations: Number of simulations to run
  #   seed: Random seed for reproducibility
  #
  # Returns:
  #   logical: TRUE if trial terminated early in this simulation
  
  set.seed(seed)
  
  if (scenario_type == "unfavorable") {
    # Create unfavorable scenario where all doses are unsafe/inactive
    # This should lead to high early termination rates
    
    # Generate unfavorable scenario data (data not used directly, but function called for consistency)
    generate_unfavorable_scenario_data(
      config = config,
      n_patients_per_dose = config$cohort_size * config$n_stages,
      seed = seed
    )
    
    # Create probability matrices for the simulation
    unfavorable_probs <- create_unfavorable_probability_matrices(
      n_doses = length(config$dose_levels)
    )
    
    # Run trial simulation
    results <- run_trial_simulation(
      trial_config = config,
      p_YI = unfavorable_probs$p_YI,
      p_YT_given_I = unfavorable_probs$p_YT_given_I,
      p_YE_given_I = unfavorable_probs$p_YE_given_I,
      rho0 = config$rho0,
      rho1 = config$rho1
    )
    
    # Check if trial terminated early
    return(results$terminated_early)
    
  } else if (scenario_type == "flat_null") {
    # Use flat null scenario for early termination calibration
    # This should lead to moderate early termination rates
    
    # Generate flat scenario data
    data <- generate_flat_scenario_data(
      config = config,
      phi_I_lower = config$phi_I_lower,
      phi_E_lower = config$phi_E_lower,
      toxicity_low = config$toxicity_low,
      n_patients_per_dose = config$cohort_size * config$n_stages,
      seed = seed
    )
    
    # Create probability matrices for the simulation
    flat_probs <- create_flat_probability_matrices(
      n_doses = length(config$dose_levels),
      phi_I_lower = config$phi_I_lower,
      phi_E_lower = config$phi_E_lower,
      toxicity_low = config$toxicity_low
    )
    
    # Run trial simulation
    results <- run_trial_simulation(
      trial_config = config,
      p_YI = flat_probs$p_YI,
      p_YT_given_I = flat_probs$p_YT_given_I,
      p_YE_given_I = flat_probs$p_YE_given_I,
      rho0 = config$rho0,
      rho1 = config$rho1
    )
    
    # Check if trial terminated early
    return(results$terminated_early)
    
  } else {
    stop("Unsupported scenario type: ", scenario_type)
  }
}

generate_unfavorable_scenario_data <- function(config, n_patients_per_dose = 10, seed = 123) {
  # Generate unfavorable scenario data where all doses are unsafe/inactive
  #
  # Args:
  #   config: Trial configuration
  #   n_patients_per_dose: Number of patients per dose level
  #   seed: Random seed for reproducibility
  #
  # Returns:
  #   data.frame: Generated trial data
  
  set.seed(seed)
  
  n_doses <- length(config$dose_levels)
  
  # Create unfavorable probability matrices
  # High toxicity, low efficacy, low immune response for all doses
  p_YI_unfavorable <- rep(0.10, n_doses)  # Very low immune response
  
  # High toxicity for all doses, both I=0 and I=1
  p_YT_given_I_unfavorable <- matrix(rep(0.80, 2 * n_doses), ncol = 2, byrow = TRUE)
  
  # Low efficacy for all doses, both I=0 and I=1
  p_YE_given_I_unfavorable <- matrix(rep(0.10, 2 * n_doses), ncol = 2, byrow = TRUE)
  
  # Generate data using existing simulation function
  data <- simulate_data_gumbel(
    n_per_dose_vector = rep(n_patients_per_dose, n_doses),
    dose_levels = config$dose_levels,
    p_YI = p_YI_unfavorable,
    p_YT_given_I = p_YT_given_I_unfavorable,
    p_YE_given_I = p_YE_given_I_unfavorable,
    rho0 = config$rho0,
    rho1 = config$rho1
  )
  
  return(data)
}

create_unfavorable_probability_matrices <- function(n_doses) {
  # Create probability matrices for unfavorable scenario
  #
  # Args:
  #   n_doses: Number of dose levels
  #
  # Returns:
  #   list: Probability matrices for unfavorable scenario
  
  # Very low immune response for all doses
  p_YI <- rep(0.10, n_doses)
  
  # High toxicity for all doses, both I=0 and I=1
  p_YT_given_I <- matrix(rep(0.80, 2 * n_doses), ncol = 2, byrow = TRUE)
  
  # Low efficacy for all doses, both I=0 and I=1
  p_YE_given_I <- matrix(rep(0.10, 2 * n_doses), ncol = 2, byrow = TRUE)
  
  return(list(
    p_YI = p_YI,
    p_YT_given_I = p_YT_given_I,
    p_YE_given_I = p_YE_given_I
  ))
}

calibrate_early_termination <- function(target_rate = 0.80, 
                                       scenario_config, 
                                       n_simulations = 1000,
                                       scenario_type = "unfavorable",
                                       threshold_range = seq(0.7, 0.99, by = 0.01),
                                       threshold_type = "c_T") {
  # Calibrate early termination parameters to achieve target termination rate.
  #
  # Args:
  #   target_rate: Target early termination rate (default: 0.80)
  #   scenario_config: Configuration for unfavorable scenario
  #   n_simulations: Number of simulations per threshold value
  #   scenario_type: Type of scenario ("unfavorable", "flat_null")
  #   threshold_range: Range of threshold values to test
  #   threshold_type: Which threshold to calibrate ("c_T", "c_E", "c_I")
  #
  # Returns:
  #   list: Calibration results with optimal threshold and performance curves
  
  cat("Starting early termination calibration...\n")
  cat("Target termination rate:", target_rate, "\n")
  cat("Scenario type:", scenario_type, "\n")
  cat("Threshold type:", threshold_type, "\n")
  cat("Number of simulations per threshold:", n_simulations, "\n")
  cat("Threshold range:", min(threshold_range), "to", max(threshold_range), "\n\n")
  
  calibration_results <- data.frame(
    threshold = numeric(),
    termination_rate = numeric(),
    termination_rate_lower = numeric(),
    termination_rate_upper = numeric(),
    n_simulations = numeric(),
    stringsAsFactors = FALSE
  )
  
  total_iterations <- length(threshold_range)
  
  for (i in seq_along(threshold_range)) {
    threshold <- threshold_range[i]
    
    cat(sprintf("Progress: %d/%d - Testing %s = %.3f\n", i, total_iterations, threshold_type, threshold))
    
    # Update config with current threshold
    config <- scenario_config
    config[[threshold_type]] <- threshold
    
    # Run simulations
    termination_results <- replicate(n_simulations, {
      # Use different seeds for each simulation
      seed <- sample.int(1000000, 1)
      return(run_early_termination_simulation(config, scenario_type, 1, seed))
    })
    
    # Calculate termination rate and confidence interval
    termination_rate <- mean(termination_results)
    ci <- binom.test(sum(termination_results), n_simulations)$conf.int
    
    # Store results
    calibration_results <- rbind(calibration_results, data.frame(
      threshold = threshold,
      termination_rate = termination_rate,
      termination_rate_lower = ci[1],
      termination_rate_upper = ci[2],
      n_simulations = n_simulations,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("  Early termination rate = %.3f (95%% CI: [%.3f, %.3f])\n",
                termination_rate, ci[1], ci[2]))
    cat("\n")
  }
  
  # Find optimal threshold
  optimal_idx <- which.min(abs(calibration_results$termination_rate - target_rate))
  optimal_threshold <- calibration_results$threshold[optimal_idx]
  optimal_rate <- calibration_results$termination_rate[optimal_idx]
  
  cat("=== EARLY TERMINATION CALIBRATION COMPLETE ===\n")
  cat(sprintf("Optimal %s = %.3f\n", threshold_type, optimal_threshold))
  cat(sprintf("Achieved termination rate = %.3f (target: %.3f)\n", optimal_rate, target_rate))
  cat(sprintf("Difference from target = %.3f\n", abs(optimal_rate - target_rate)))
  cat("==============================================\n\n")
  
  return(list(
    calibration_results = calibration_results,
    optimal_threshold = optimal_threshold,
    optimal_rate = optimal_rate,
    target_rate = target_rate,
    threshold_type = threshold_type,
    scenario_type = scenario_type,
    n_simulations = n_simulations
  ))
}

validate_early_termination_calibration <- function(calibration_results, n_validation_simulations = 1000) {
  # Validate the calibrated early termination parameters with additional simulations.
  #
  # Args:
  #   calibration_results: Results from calibrate_early_termination
  #   n_validation_simulations: Number of validation simulations
  #
  # Returns:
  #   list: Validation results
  
  cat("Validating calibrated early termination parameters...\n")
  
  # Get the optimal config
  optimal_threshold <- calibration_results$optimal_threshold
  threshold_type <- calibration_results$threshold_type
  scenario_type <- calibration_results$scenario_type
  
  # Create config with optimal threshold
  config <- get("unfavorable_scenario_config", envir = .GlobalEnv)  # Use unfavorable scenario config as base
  config[[threshold_type]] <- optimal_threshold
  
  # Run validation simulations
  validation_terminations <- replicate(n_validation_simulations, {
    seed <- sample.int(1000000, 1)
    return(run_early_termination_simulation(config, scenario_type, 1, seed))
  })
  
  # Calculate validation statistics
  validation_rate <- mean(validation_terminations)
  validation_ci <- binom.test(sum(validation_terminations), n_validation_simulations)$conf.int
  
  cat(sprintf("Validation results:\n"))
  cat(sprintf("  %s = %.3f\n", threshold_type, optimal_threshold))
  cat(sprintf("  Termination rate = %.3f (95%% CI: [%.3f, %.3f])\n", 
              validation_rate, validation_ci[1], validation_ci[2]))
  cat(sprintf("  Target rate = %.3f\n", calibration_results$target_rate))
  cat(sprintf("  Difference = %.3f\n", abs(validation_rate - calibration_results$target_rate)))
  
  return(list(
    optimal_threshold = optimal_threshold,
    threshold_type = threshold_type,
    validation_rate = validation_rate,
    validation_ci = validation_ci,
    target_rate = calibration_results$target_rate,
    n_validation_simulations = n_validation_simulations
  ))
}

run_quick_early_termination_calibration <- function(target_rate = 0.80, n_simulations = 100) {
  # Run a quick early termination calibration with reduced simulations for testing.
  #
  # NOTE: Since achieving 80% early termination with only c_T calibration is mathematically
  # challenging in highly unfavorable scenarios, we calibrate BOTH c_T and c_E together
  # to find a feasible parameter combination.
  #
  # Args:
  #   target_rate: Target early termination rate
  #   n_simulations: Number of simulations per threshold value
  #
  # Returns:
  #   list: Quick calibration results
  
  cat("Running quick early termination calibration (optimizing c_T and c_E jointly)...\n")
  
  # For practical calibration, we test a grid of (c_T, c_E) combinations
  # to find which achieves closest to 80% termination rate
  c_T_range <- seq(0.60, 0.90, by = 0.10)
  c_E_range <- seq(0.60, 0.85, by = 0.10)
  
  calibration_results_grid <- expand.grid(c_T = c_T_range, c_E = c_E_range)
  calibration_results_grid$termination_rate <- NA
  
  for (i in seq_len(nrow(calibration_results_grid))) {
    c_T_val <- calibration_results_grid$c_T[i]
    c_E_val <- calibration_results_grid$c_E[i]
    
    cat(sprintf("Testing c_T=%.2f, c_E=%.2f\n", c_T_val, c_E_val))
    
    # Update config with both thresholds
    config <- get("unfavorable_scenario_config", envir = .GlobalEnv)
    config$c_T <- c_T_val
    config$c_E <- c_E_val
    
    # Run simulations
    termination_results <- replicate(n_simulations, {
      seed <- sample.int(1000000, 1)
      return(run_early_termination_simulation(config, "unfavorable", 1, seed))
    })
    
    termination_rate <- mean(termination_results)
    calibration_results_grid$termination_rate[i] <- termination_rate
    cat(sprintf("  → Termination rate = %.3f\n", termination_rate))
  }
  
  # Find the combination closest to target
  differences <- abs(calibration_results_grid$termination_rate - target_rate)
  best_idx <- which.min(differences)
  
  best_c_T <- calibration_results_grid$c_T[best_idx]
  best_c_E <- calibration_results_grid$c_E[best_idx]
  best_rate <- calibration_results_grid$termination_rate[best_idx]
  
  cat("\n=== EARLY TERMINATION CALIBRATION COMPLETE ===\n")
  cat(sprintf("Optimal c_T = %.2f, c_E = %.2f\n", best_c_T, best_c_E))
  cat(sprintf("Achieved termination rate = %.3f (target: %.3f)\n", best_rate, target_rate))
  cat(sprintf("Difference from target = %.3f\n", abs(best_rate - target_rate)))
  cat("==============================================\n\n")
  
  return(list(
    calibration_results = calibration_results_grid,
    optimal_c_T = best_c_T,
    optimal_c_E = best_c_E,
    optimal_threshold = best_c_T,  # Keep for backwards compatibility
    threshold_type = "c_T & c_E",
    optimal_rate = best_rate,
    target_rate = target_rate,
    scenario_type = "unfavorable",
    n_simulations = n_simulations
  ))
}

save_early_termination_results <- function(calibration_results, file_path = "results/early_termination_calibration_results.RData") {
  # Save early termination calibration results to file.
  #
  # Args:
  #   calibration_results: Results from calibrate_early_termination
  #   file_path: Path to save results
  
  # Create results directory if it doesn't exist
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  
  # Save results
  save(calibration_results, file = file_path)
  cat("Early termination calibration results saved to:", file_path, "\n")
}

load_early_termination_results <- function(file_path = "results/early_termination_calibration_results.RData") {
  # Load early termination calibration results from file.
  #
  # Args:
  #   file_path: Path to load results from
  #
  # Returns:
  #   list: Loaded calibration results
  
  if (file.exists(file_path)) {
    load(file_path)
    cat("Early termination calibration results loaded from:", file_path, "\n")
    return(get("calibration_results", envir = .GlobalEnv))
  } else {
    cat("No early termination calibration results found at:", file_path, "\n")
    return(NULL)
  }
}
