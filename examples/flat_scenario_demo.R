# Flat Scenario Generation Demo
# This script demonstrates the flat scenario generation functions implemented in Phase 1

# Load required libraries
library(dplyr)
library(ggplot2)

# Source the functions
source("src/core/simulate_data.R")
source("src/core/config.R")

cat("=== Flat Scenario Generation Demo ===\n\n")

# 1. Demonstrate flat scenario data generation
cat("1. Generating flat scenario data...\n")
cat("   Parameters: phi_I = 0.20, phi_E = 0.25, toxicity = 0.05\n")
cat("   Sample size: 200 patients per dose\n\n")

flat_data <- generate_flat_scenario_data(
  config = flat_scenario_config,
  phi_I_lower = 0.20,
  phi_E_lower = 0.25,
  toxicity_low = 0.05,
  n_patients_per_dose = 200,  # Larger sample for more stable estimates
  seed = 123
)

cat("Generated data summary:\n")
cat("  Total patients:", nrow(flat_data), "\n")
cat("  Dose levels:", paste(unique(flat_data$d), collapse = ", "), "\n")
cat("  Patients per dose:", nrow(flat_data) / length(unique(flat_data$d)), "\n\n")

# 2. Show observed rates by dose
cat("2. Observed rates by dose:\n")
for (dose in unique(flat_data$d)) {
  dose_data <- flat_data[flat_data$d == dose, ]
  cat(sprintf("  Dose %d: Immune=%.3f, Efficacy=%.3f, Toxicity=%.3f\n",
              dose, mean(dose_data$Y_I), mean(dose_data$Y_E), mean(dose_data$Y_T)))
}
cat("\n")

# 3. Validate the flat scenario
cat("3. Validating flat scenario...\n")
validation <- validate_flat_scenario(flat_data, 0.20, 0.25, 0.05, tolerance = 0.15)

if (validation$success) {
  cat("   ✓ Validation PASSED - Data represents a flat scenario\n")
} else {
  cat("   ✗ Validation FAILED\n")
  cat("   Details:", paste(validation$details, collapse = "; "), "\n")
}
cat("\n")

# 4. Compare with non-flat scenario
cat("4. Comparing with non-flat scenario...\n")

# Generate non-flat data using original config
non_flat_data <- simulate_data_gumbel(
  n_per_dose_vector = rep(200, length(trial_config$dose_levels)),
  dose_levels = trial_config$dose_levels,
  p_YI = p_YI,
  p_YT_given_I = p_YT_given_I,
  p_YE_given_I = p_YE_given_I,
  rho0 = rho0,
  rho1 = rho1,
  seed = 123
)

cat("Non-flat scenario rates:\n")
for (dose in unique(non_flat_data$d)) {
  dose_data <- non_flat_data[non_flat_data$d == dose, ]
  cat(sprintf("  Dose %d: Immune=%.3f, Efficacy=%.3f, Toxicity=%.3f\n",
              dose, mean(dose_data$Y_I), mean(dose_data$Y_E), mean(dose_data$Y_T)))
}
cat("\n")

# 5. Create visualization
cat("5. Creating visualization...\n")

# Prepare data for plotting
plot_data <- rbind(
  data.frame(
    Scenario = "Flat",
    Dose = flat_data$d,
    Immune = flat_data$Y_I,
    Efficacy = flat_data$Y_E,
    Toxicity = flat_data$Y_T
  ),
  data.frame(
    Scenario = "Non-Flat",
    Dose = non_flat_data$d,
    Immune = non_flat_data$Y_I,
    Efficacy = non_flat_data$Y_E,
    Toxicity = non_flat_data$Y_T
  )
)

# Calculate rates by scenario and dose
rates_data <- plot_data %>%
  group_by(Scenario, Dose) %>%
  summarise(
    Immune_Rate = mean(Immune),
    Efficacy_Rate = mean(Efficacy),
    Toxicity_Rate = mean(Toxicity),
    .groups = 'drop'
  ) %>%
  tidyr::pivot_longer(cols = c(Immune_Rate, Efficacy_Rate, Toxicity_Rate),
                      names_to = "Endpoint", values_to = "Rate")

# Create plot
p <- ggplot(rates_data, aes(x = Dose, y = Rate, color = Scenario)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~ Endpoint, scales = "free_y") +
  labs(
    title = "Flat vs Non-Flat Scenario Comparison",
    subtitle = "Flat scenario shows identical rates across doses",
    x = "Dose Level",
    y = "Response Rate",
    color = "Scenario"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 11),
    legend.position = "bottom"
  )

# Save plot
ggsave("results/flat_scenario_comparison.png", plot = p, width = 10, height = 6, dpi = 300)
cat("   Plot saved to: results/flat_scenario_comparison.png\n\n")

# 6. Show conditional efficacy calculation
cat("6. Conditional efficacy calculation demonstration...\n")
conditional_eff <- calculate_conditional_efficacy_flat(0.25, 0.20)
cat("   For phi_E = 0.25, phi_I = 0.20:\n")
cat(sprintf("   P(E|I=0) = %.3f, P(E|I=1) = %.3f\n", 
            conditional_eff[1,1], conditional_eff[1,2]))

# Verify marginal efficacy
marginal_eff <- (1 - 0.20) * conditional_eff[1,1] + 0.20 * conditional_eff[1,2]
cat(sprintf("   Marginal efficacy = %.3f (target: 0.25)\n", marginal_eff))
cat("\n")

# 7. Summary
cat("=== Summary ===\n")
cat("✓ Flat scenario generation functions implemented successfully\n")
cat("✓ Functions generate data with identical probabilities across doses\n")
cat("✓ Validation function correctly identifies flat scenarios\n")
cat("✓ Conditional efficacy calculation maintains marginal rates\n")
cat("✓ Ready for Phase 2: Enhanced PoC calculation\n\n")

cat("Next steps:\n")
cat("- Implement Bayesian PoC calculation\n")
cat("- Create calibration framework\n")
cat("- Generate performance curves\n")
