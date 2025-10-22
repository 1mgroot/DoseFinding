# Example: How to use the new logging control feature
# This example shows how to reduce logging during calibration

# Source the required files
source("src/core/config.R")
source("src/optimization/poc_calibration.R")

cat("=== Logging Control Example ===\n\n")

# Example 1: Run calibration with minimal logging (default behavior)
cat("Example 1: Running PoC calibration with minimal logging\n")
cat("This will run quietly without verbose output...\n\n")

# The calibration functions now automatically use minimal logging
# You can run calibration without seeing all the detailed workflow logs
quick_results <- run_quick_calibration(target_rate = 0.10, n_simulations = 10)

cat("Calibration completed with minimal logging!\n\n")

# Example 2: Run a single trial with verbose logging for debugging
cat("Example 2: Running a single trial with verbose logging for debugging\n")
cat("This will show detailed workflow information...\n\n")

# Create a verbose config for debugging
debug_config <- trial_config
debug_config$n_stages <- 2  # Short trial for example
debug_config$cohort_size <- 3

# Run with verbose logging
debug_results <- run_trial_simulation(debug_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1)

cat("\nDebug trial completed with verbose logging!\n\n")

# Example 3: Manual control of logging
cat("Example 3: Manual control of logging parameters\n")
cat("You can manually set verbose_logging = FALSE in any config:\n\n")

manual_config <- trial_config
manual_config$verbose_logging <- FALSE  # Disable verbose logging
manual_config$log_early_termination <- FALSE  # Also disable early termination logging

cat("Config with manual logging control:\n")
cat("verbose_logging:", manual_config$verbose_logging, "\n")
cat("log_early_termination:", manual_config$log_early_termination, "\n\n")

cat("=== Example Complete ===\n")
cat("The logging control feature allows you to:\n")
cat("1. Run calibrations quietly (automatic)\n")
cat("2. Debug individual trials with verbose output\n")
cat("3. Manually control logging levels as needed\n")
