# PoC Calibration Framework
# This file implements the calibration framework for PoC (Probability of Correct Selection) parameters

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

run_calibration_simulation <- function(config, scenario_type = "flat_null", n_simulations = 1000, seed = 123) {
  # Run a single calibration simulation to determine if PoC is detected.
  #
  # Args:
  #   config: Trial configuration with current C_poc value
  #   scenario_type: Type of scenario ("flat_null", "unfavorable", etc.)
  #   n_simulations: Number of simulations to run
  #   seed: Random seed for reproducibility
  #
  # Returns:
  #   logical: TRUE if PoC was detected in this simulation
  
  set.seed(seed)
  
  if (scenario_type == "flat_null") {
    # Create a more lenient config for flat scenario calibration
    # The flat scenario should allow some doses to be admissible for PoC testing
    lenient_config <- config
    lenient_config$c_T <- 0.5  # More lenient safety threshold
    lenient_config$c_E <- 0.5  # More lenient efficacy threshold
    lenient_config$c_I <- 0.5  # More lenient immune response threshold
    
    # Generate flat scenario data
    data <- generate_flat_scenario_data(
      config = lenient_config,
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
    
    # Run trial simulation with lenient config
    results <- run_trial_simulation(
      trial_config = lenient_config,
      p_YI = flat_probs$p_YI,
      p_YT_given_I = flat_probs$p_YT_given_I,
      p_YE_given_I = flat_probs$p_YE_given_I,
      rho0 = config$rho0,
      rho1 = config$rho1
    )
    
    # Check if PoC was detected (trial completed with PoC validation)
    # If trial terminated early, PoC was not detected
    if (results$terminated_early) {
      return(FALSE)
    } else {
      return(results$poc_validated)
    }
    
  } else {
    stop("Unsupported scenario type: ", scenario_type)
  }
}

calibrate_c_poc <- function(target_rate = 0.10, 
                           flat_scenario_config, 
                           n_simulations = 1000,
                           c_poc_range = seq(0.5, 0.99, by = 0.01),
                           n_parallel = 1) {
  # Calibrate C_poc to achieve target PoC detection rate in flat scenarios.
  #
  # Args:
  #   target_rate: Target PoC detection rate (default: 0.10 for 10%)
  #   flat_scenario_config: Configuration for flat null scenario
  #   n_simulations: Number of simulations per C_poc value
  #   c_poc_range: Range of C_poc values to test
  #   n_parallel: Number of parallel simulations (future enhancement)
  #
  # Returns:
  #   list: Calibration results with optimal C_poc and performance curves
  
  cat("Starting PoC calibration...\n")
  cat("Target detection rate:", target_rate, "\n")
  cat("Number of simulations per C_poc:", n_simulations, "\n")
  cat("C_poc range:", min(c_poc_range), "to", max(c_poc_range), "\n\n")
  
  calibration_results <- data.frame(
    c_poc = numeric(),
    poc_detection_rate = numeric(),
    poc_detection_rate_lower = numeric(),
    poc_detection_rate_upper = numeric(),
    n_simulations = numeric(),
    stringsAsFactors = FALSE
  )
  
  total_iterations <- length(c_poc_range)
  
  for (i in seq_along(c_poc_range)) {
    c_poc <- c_poc_range[i]
    
    cat(sprintf("Progress: %d/%d - Testing C_poc = %.3f\n", i, total_iterations, c_poc))
    
    # Update config with current c_poc
    config <- flat_scenario_config
    config$c_poc <- c_poc
    
    # Run simulations
    poc_detections <- replicate(n_simulations, {
      # Use different seeds for each simulation
      seed <- sample.int(1000000, 1)
      return(run_calibration_simulation(config, "flat_null", 1, seed))
    })
    
    # Calculate detection rate and confidence interval
    detection_rate <- mean(poc_detections)
    ci <- binom.test(sum(poc_detections), n_simulations)$conf.int
    
    # Store results
    calibration_results <- rbind(calibration_results, data.frame(
      c_poc = c_poc,
      poc_detection_rate = detection_rate,
      poc_detection_rate_lower = ci[1],
      poc_detection_rate_upper = ci[2],
      n_simulations = n_simulations,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("  PoC detection rate = %.3f (95%% CI: [%.3f, %.3f])\n",
                detection_rate, ci[1], ci[2]))
    cat("\n")
  }
  
  # Find optimal c_poc
  optimal_idx <- which.min(abs(calibration_results$poc_detection_rate - target_rate))
  optimal_c_poc <- calibration_results$c_poc[optimal_idx]
  optimal_rate <- calibration_results$poc_detection_rate[optimal_idx]
  
  cat("=== CALIBRATION COMPLETE ===\n")
  cat(sprintf("Optimal C_poc = %.3f\n", optimal_c_poc))
  cat(sprintf("Achieved detection rate = %.3f (target: %.3f)\n", optimal_rate, target_rate))
  cat(sprintf("Difference from target = %.3f\n", abs(optimal_rate - target_rate)))
  cat("============================\n\n")
  
  return(list(
    calibration_results = calibration_results,
    optimal_c_poc = optimal_c_poc,
    optimal_rate = optimal_rate,
    target_rate = target_rate,
    n_simulations = n_simulations
  ))
}

validate_calibration <- function(calibration_results, n_validation_simulations = 1000) {
  # Validate the calibrated C_poc with additional simulations.
  #
  # Args:
  #   calibration_results: Results from calibrate_c_poc
  #   n_validation_simulations: Number of validation simulations
  #
  # Returns:
  #   list: Validation results
  
  cat("Validating calibrated C_poc...\n")
  
  # Get the optimal config
  optimal_c_poc <- calibration_results$optimal_c_poc
  config <- flat_scenario_config
  config$c_poc <- optimal_c_poc
  
  # Run validation simulations
  validation_detections <- replicate(n_validation_simulations, {
    seed <- sample.int(1000000, 1)
    return(run_calibration_simulation(config, "flat_null", 1, seed))
  })
  
  # Calculate validation statistics
  validation_rate <- mean(validation_detections)
  validation_ci <- binom.test(sum(validation_detections), n_validation_simulations)$conf.int
  
  cat(sprintf("Validation results:\n"))
  cat(sprintf("  C_poc = %.3f\n", optimal_c_poc))
  cat(sprintf("  Detection rate = %.3f (95%% CI: [%.3f, %.3f])\n", 
              validation_rate, validation_ci[1], validation_ci[2]))
  cat(sprintf("  Target rate = %.3f\n", calibration_results$target_rate))
  cat(sprintf("  Difference = %.3f\n", abs(validation_rate - calibration_results$target_rate)))
  
  return(list(
    optimal_c_poc = optimal_c_poc,
    validation_rate = validation_rate,
    validation_ci = validation_ci,
    target_rate = calibration_results$target_rate,
    n_validation_simulations = n_validation_simulations
  ))
}

run_quick_calibration <- function(target_rate = 0.10, n_simulations = 100) {
  # Run a quick calibration with reduced simulations for testing.
  #
  # Args:
  #   target_rate: Target PoC detection rate
  #   n_simulations: Number of simulations per C_poc value
  #
  # Returns:
  #   list: Quick calibration results
  
  cat("Running quick calibration (reduced simulations for testing)...\n")
  
  # Use a smaller range for quick calibration
  c_poc_range <- seq(0.7, 0.95, by = 0.05)
  
  results <- calibrate_c_poc(
    target_rate = target_rate,
    flat_scenario_config = flat_scenario_config,
    n_simulations = n_simulations,
    c_poc_range = c_poc_range
  )
  
  return(results)
}

save_calibration_results <- function(calibration_results, file_path = "results/poc_calibration_results.RData") {
  # Save calibration results to file.
  #
  # Args:
  #   calibration_results: Results from calibrate_c_poc
  #   file_path: Path to save results
  
  # Create results directory if it doesn't exist
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  
  # Save results
  save(calibration_results, file = file_path)
  cat("Calibration results saved to:", file_path, "\n")
}

load_calibration_results <- function(file_path = "results/poc_calibration_results.RData") {
  # Load calibration results from file.
  #
  # Args:
  #   file_path: Path to load results from
  #
  # Returns:
  #   list: Loaded calibration results
  
  if (file.exists(file_path)) {
    load(file_path)
    cat("Calibration results loaded from:", file_path, "\n")
    return(calibration_results)
  } else {
    cat("No calibration results found at:", file_path, "\n")
    return(NULL)
  }
}
