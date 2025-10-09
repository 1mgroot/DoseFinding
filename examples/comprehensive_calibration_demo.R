# Comprehensive Calibration Demo
# This script demonstrates the complete calibration workflow for both PoC and early termination

# Load required libraries
library(ggplot2)
library(dplyr)

# Source required functions
source("src/core/config.R")
source("src/optimization/poc_calibration.R")
source("src/optimization/early_termination_calibration.R")
source("src/utils/calibration_plots.R")

# Set up output directory
output_dir <- "results/comprehensive_calibration"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== COMPREHENSIVE CALIBRATION DEMO ===\n")
cat("This demo will run both PoC and early termination calibration\n")
cat("Output directory:", output_dir, "\n\n")

# =============================================================================
# PHASE 1: PoC Calibration
# =============================================================================

cat("PHASE 1: PoC Calibration\n")
cat("========================\n")

# Run PoC calibration with reduced simulations for demo
cat("Running PoC calibration (reduced simulations for demo)...\n")
poc_calibration_results <- run_quick_calibration(
  target_rate = calibration_config$poc_target_rate,
  n_simulations = 200  # Reduced for demo
)

# Validate PoC calibration
cat("\nValidating PoC calibration...\n")
poc_validation_results <- validate_calibration(
  poc_calibration_results, 
  n_validation_simulations = 500
)

# Add validation results to calibration results
poc_calibration_results$validation_results <- poc_validation_results

# Save PoC calibration results
save_calibration_results(
  poc_calibration_results, 
  file.path(output_dir, "poc_calibration_results.RData")
)

cat("PoC calibration completed!\n")
cat("Optimal C_poc =", poc_calibration_results$optimal_c_poc, "\n")
cat("Achieved rate =", round(poc_calibration_results$optimal_rate, 3), "\n")
cat("Target rate =", poc_calibration_results$target_rate, "\n\n")

# =============================================================================
# PHASE 2: Early Termination Calibration
# =============================================================================

cat("PHASE 2: Early Termination Calibration\n")
cat("======================================\n")

# Run early termination calibration with reduced simulations for demo
cat("Running early termination calibration (reduced simulations for demo)...\n")
early_termination_results <- run_quick_early_termination_calibration(
  target_rate = calibration_config$early_termination_target_rate,
  n_simulations = 200  # Reduced for demo
)

# Validate early termination calibration
cat("\nValidating early termination calibration...\n")
early_termination_validation_results <- validate_early_termination_calibration(
  early_termination_results, 
  n_validation_simulations = 500
)

# Add validation results to calibration results
early_termination_results$validation_results <- early_termination_validation_results

# Save early termination calibration results
save_early_termination_results(
  early_termination_results, 
  file.path(output_dir, "early_termination_calibration_results.RData")
)

cat("Early termination calibration completed!\n")
cat("Optimal", early_termination_results$threshold_type, "=", early_termination_results$optimal_threshold, "\n")
cat("Achieved rate =", round(early_termination_results$optimal_rate, 3), "\n")
cat("Target rate =", early_termination_results$target_rate, "\n\n")

# =============================================================================
# PHASE 3: Visualization
# =============================================================================

cat("PHASE 3: Creating Calibration Visualizations\n")
cat("============================================\n")

# Create PoC calibration curve
cat("Creating PoC calibration curve...\n")
poc_plot <- plot_poc_calibration_curve(
  poc_calibration_results,
  target_rate = calibration_config$poc_target_rate,
  save_path = file.path(output_dir, "poc_calibration_curve.png")
)

# Create early termination calibration curve
cat("Creating early termination calibration curve...\n")
early_termination_plot <- plot_early_termination_curve(
  early_termination_results,
  target_rate = calibration_config$early_termination_target_rate,
  save_path = file.path(output_dir, "early_termination_calibration_curve.png")
)

# Create combined performance curves
cat("Creating combined performance curves...\n")
combined_calibration_data <- list(
  poc_calibration = poc_calibration_results,
  termination_calibration = early_termination_results
)

combined_plot <- plot_threshold_performance_curves(
  combined_calibration_data,
  save_path = file.path(output_dir, "combined_performance_curves.png")
)

# Create calibration summary plots
cat("Creating calibration summary plots...\n")
poc_summary_plot <- plot_calibration_summary(
  poc_calibration_results,
  save_path = file.path(output_dir, "poc_calibration_summary.png")
)

# Create comprehensive calibration report
cat("Creating comprehensive calibration report...\n")
calibration_report <- create_calibration_report(
  poc_calibration_results,
  output_dir = file.path(output_dir, "poc_report")
)

# =============================================================================
# PHASE 4: Integration Testing
# =============================================================================

cat("PHASE 4: Integration Testing\n")
cat("============================\n")

# Test calibrated parameters in a full trial simulation
cat("Testing calibrated parameters in full trial simulation...\n")

# Create config with calibrated parameters
calibrated_config <- flat_scenario_config
calibrated_config$c_poc <- poc_calibration_results$optimal_c_poc
calibrated_config$c_T <- early_termination_results$optimal_threshold

# Run a test simulation with calibrated parameters
cat("Running test simulation with calibrated parameters...\n")
test_results <- run_trial_simulation(
  trial_config = calibrated_config,
  p_YI = rep(calibrated_config$phi_I_lower, length(calibrated_config$dose_levels)),
  p_YT_given_I = matrix(rep(calibrated_config$toxicity_low, 2 * length(calibrated_config$dose_levels)), ncol = 2, byrow = TRUE),
  p_YE_given_I = calculate_conditional_efficacy_flat(calibrated_config$phi_E_lower, calibrated_config$phi_I_lower),
  rho0 = calibrated_config$rho0,
  rho1 = calibrated_config$rho1
)

# Save test results
save(test_results, file = file.path(output_dir, "calibrated_parameters_test.RData"))

cat("Integration test completed!\n")
cat("Test simulation results:\n")
cat("  Trial terminated early:", test_results$terminated_early, "\n")
cat("  PoC validated:", test_results$poc_validated, "\n")
cat("  Final dose selected:", test_results$final_dose, "\n")
cat("  Total participants:", nrow(test_results$all_data), "\n\n")

# =============================================================================
# PHASE 5: Summary Report
# =============================================================================

cat("PHASE 5: Summary Report\n")
cat("=======================\n")

# Create summary report
summary_report <- data.frame(
  Parameter = c("C_poc (PoC)", "C_T (Early Termination)", "PoC Target Rate", "PoC Achieved Rate", 
                "Early Termination Target", "Early Termination Achieved"),
  Value = c(
    round(poc_calibration_results$optimal_c_poc, 3),
    round(early_termination_results$optimal_threshold, 3),
    calibration_config$poc_target_rate,
    round(poc_calibration_results$optimal_rate, 3),
    calibration_config$early_termination_target_rate,
    round(early_termination_results$optimal_rate, 3)
  ),
  Status = c(
    "Calibrated",
    "Calibrated", 
    "Target",
    ifelse(abs(poc_calibration_results$optimal_rate - calibration_config$poc_target_rate) <= calibration_config$poc_tolerance, "✓ Achieved", "⚠ Needs Adjustment"),
    "Target",
    ifelse(abs(early_termination_results$optimal_rate - calibration_config$early_termination_target_rate) <= calibration_config$early_termination_tolerance, "✓ Achieved", "⚠ Needs Adjustment")
  )
)

# Save summary report
write.csv(summary_report, file.path(output_dir, "calibration_summary_report.csv"), row.names = FALSE)

# Print summary
cat("CALIBRATION SUMMARY:\n")
print(summary_report)
cat("\n")

# Save all results
save(poc_calibration_results, early_termination_results, summary_report, 
     file = file.path(output_dir, "comprehensive_calibration_results.RData"))

cat("=== COMPREHENSIVE CALIBRATION DEMO COMPLETED ===\n")
cat("All results saved to:", output_dir, "\n")
cat("Generated files:\n")
cat("  - poc_calibration_curve.png\n")
cat("  - early_termination_calibration_curve.png\n")
cat("  - combined_performance_curves.png\n")
cat("  - poc_calibration_summary.png\n")
cat("  - calibration_summary_report.csv\n")
cat("  - comprehensive_calibration_results.RData\n")
cat("  - poc_report/ (directory with detailed PoC report)\n")
cat("================================================\n")
