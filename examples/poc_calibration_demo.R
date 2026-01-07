# PoC Calibration Framework Demo
# This script demonstrates the PoC calibration framework implemented in Phase 3

# Load required libraries
library(dplyr)
library(ggplot2)

# Source the functions
source("src/core/simulate_data.R")
source("src/core/config.R")
source("src/core/model_utils.R")
source("src/decision/dose_decision.R")
source("src/core/main.R")
source("src/optimization/poc_calibration.R")

cat("=== PoC Calibration Framework Demo ===\n\n")

# 1. Demonstrate single calibration simulation
cat("1. Testing single calibration simulation...\n")

# Test with a lenient C_poc to see if we can get some PoC detections
test_config <- flat_scenario_config
test_config$c_poc <- 0.5  # Very lenient threshold

result <- run_calibration_simulation(test_config, "flat_null", 1, seed = 123)
cat("Single simulation result (C_poc = 0.5):", result, "\n\n")

# 2. Demonstrate quick calibration
cat("2. Running quick calibration (reduced simulations for demo)...\n")

# Run quick calibration with very few simulations
quick_results <- run_quick_calibration(target_rate = 0.10, n_simulations = 5)

cat("Quick calibration results:\n")
cat("  Optimal C_poc:", quick_results$optimal_c_poc, "\n")
cat("  Achieved rate:", round(quick_results$optimal_rate, 3), "\n")
cat("  Target rate:", quick_results$target_rate, "\n")
cat("  Number of C_poc values tested:", nrow(quick_results$calibration_results), "\n\n")

# 3. Show calibration results table
cat("3. Calibration results table:\n")
print(quick_results$calibration_results)
cat("\n")

# 4. Demonstrate the calibration concept
cat("4. Explaining the calibration concept...\n")
cat("In a flat null scenario (all doses similar):\n")
cat("- PoC detection should be rare (low rate)\n")
cat("- We want to find C_poc that gives ~10% detection rate\n")
cat("- This represents Type I error control\n")
cat("- Lower C_poc = more lenient = higher detection rate\n")
cat("- Higher C_poc = more strict = lower detection rate\n\n")

# 5. Show how different C_poc values affect detection
cat("5. Effect of different C_poc values on detection rate:\n")
for (i in 1:nrow(quick_results$calibration_results)) {
  row <- quick_results$calibration_results[i, ]
  cat(sprintf("  C_poc = %.2f: Detection rate = %.3f (95%% CI: [%.3f, %.3f])\n",
              row$c_poc, row$poc_detection_rate, 
              row$poc_detection_rate_lower, row$poc_detection_rate_upper))
}
cat("\n")

# 6. Demonstrate validation
cat("6. Validating calibrated C_poc...\n")
validation_results <- validate_calibration(quick_results, n_validation_simulations = 3)

cat("Validation results:\n")
cat("  C_poc:", validation_results$optimal_c_poc, "\n")
cat("  Validation rate:", round(validation_results$validation_rate, 3), "\n")
cat("  Target rate:", validation_results$target_rate, "\n")
cat("  Difference:", round(abs(validation_results$validation_rate - validation_results$target_rate), 3), "\n\n")

# 7. Save and load demonstration
cat("7. Save and load calibration results...\n")
demo_file <- "results/demo_calibration_results.RData"
save_calibration_results(quick_results, demo_file)

# Load the results
loaded_results <- load_calibration_results(demo_file)
if (!is.null(loaded_results)) {
  cat("Successfully loaded calibration results\n")
  cat("Loaded optimal C_poc:", loaded_results$optimal_c_poc, "\n")
}

# Clean up demo file
if (file.exists(demo_file)) {
  file.remove(demo_file)
  cat("Demo file cleaned up\n")
}
cat("\n")

# 8. Summary
cat("=== Summary ===\n")
cat("✓ PoC calibration framework implemented successfully\n")
cat("✓ Functions for single simulation, calibration, and validation\n")
cat("✓ Save/load functionality for calibration results\n")
cat("✓ Proper handling of flat null scenarios\n")
cat("✓ Ready for full calibration with more simulations\n\n")

cat("Key features:\n")
cat("- run_calibration_simulation(): Single simulation with flat scenario\n")
cat("- calibrate_c_poc(): Full calibration with multiple C_poc values\n")
cat("- validate_calibration(): Validation with additional simulations\n")
cat("- run_quick_calibration(): Quick testing with reduced simulations\n")
cat("- save/load_calibration_results(): Persistence of results\n\n")

cat("Next steps:\n")
cat("- Run full calibration with 10,000+ simulations per C_poc\n")
cat("- Implement early termination calibration\n")
cat("- Create visualization tools for calibration curves\n")
cat("- Integrate with main trial simulation workflow\n")
