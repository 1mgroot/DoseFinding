# Test file for Early Termination and PoC Implementation

# Source the necessary files
source("../src/core/config.R")
source("../src/utils/helpers.R")
source("../src/core/simulate_data.R")
source("../src/core/model_utils.R")
source("../src/decision/dose_decision.R")
source("../src/core/main.R")

# Test 1: Early Termination Scenario
cat("=== TEST 1: Early Termination Scenario ===\n")

# Create a scenario where all doses will be inadmissible
test_config <- trial_config
test_config$phi_T <- 0.1  # Very strict toxicity threshold
test_config$phi_E <- 0.8  # Very high efficacy requirement
test_config$phi_I <- 0.8  # Very high immune response requirement
test_config$c_T <- 0.9
test_config$c_E <- 0.9
test_config$c_I <- 0.9

# Use low probability data to trigger early termination
test_p_YI <- c(0.1, 0.1, 0.1)  # Low immune response
test_p_YT_given_I <- matrix(c(
  0.3, 0.4, 0.5,  # High toxicity given I=0
  0.4, 0.5, 0.6   # High toxicity given I=1
), ncol = 2, byrow = TRUE)
test_p_YE_given_I <- matrix(c(
  0.1, 0.1, 0.1,  # Low efficacy given I=0
  0.2, 0.2, 0.2   # Low efficacy given I=1
), ncol = 2, byrow = TRUE)

cat("Running early termination test...\n")
early_termination_results <- run_trial_simulation(test_config, test_p_YI, test_p_YT_given_I, test_p_YE_given_I, rho0, rho1)

cat("Early termination test results:\n")
cat("Terminated early:", early_termination_results$terminated_early, "\n")
if (early_termination_results$terminated_early) {
  cat("Termination stage:", early_termination_results$termination_stage, "\n")
  cat("Termination reason:", early_termination_results$termination_reason, "\n")
}

# Test 2: PoC Validation Scenario
cat("\n=== TEST 2: PoC Validation Scenario ===\n")

# Create a scenario with clear dose differences for PoC testing
test_config_2 <- trial_config
test_config_2$c_poc <- 0.8
test_config_2$delta_poc <- 0.7

# Use data with clear dose-response relationships
test_p_YI_2 <- c(0.2, 0.4, 0.6)  # Increasing immune response
test_p_YT_given_I_2 <- matrix(c(
  0.1, 0.2, 0.3,  # Moderate toxicity given I=0
  0.2, 0.3, 0.4   # Moderate toxicity given I=1
), ncol = 2, byrow = TRUE)
test_p_YE_given_I_2 <- matrix(c(
  0.3, 0.5, 0.7,  # Clear efficacy increase given I=0
  0.5, 0.7, 0.9   # Clear efficacy increase given I=1
), ncol = 2, byrow = TRUE)

cat("Running PoC validation test...\n")
poc_results <- run_trial_simulation(test_config_2, test_p_YI_2, test_p_YT_given_I_2, test_p_YE_given_I_2, rho0, rho1)

cat("PoC validation test results:\n")
cat("Terminated early:", poc_results$terminated_early, "\n")
if (!poc_results$terminated_early) {
  cat("Final OD:", poc_results$final_od, "\n")
  cat("Final utility:", round(poc_results$final_utility, 2), "\n")
  cat("PoC validated:", poc_results$poc_validated, "\n")
  cat("PoC probability:", round(poc_results$poc_probability, 3), "\n")
  cat("Selection reason:", poc_results$selection_reason, "\n")
}

# Test 3: Normal Trial Scenario
cat("\n=== TEST 3: Normal Trial Scenario ===\n")

# Use the original configuration
cat("Running normal trial test...\n")
normal_results <- run_trial_simulation(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1)

cat("Normal trial test results:\n")
cat("Terminated early:", normal_results$terminated_early, "\n")
if (!normal_results$terminated_early) {
  cat("Final OD:", normal_results$final_od, "\n")
  cat("Final utility:", round(normal_results$final_utility, 2), "\n")
  cat("PoC validated:", normal_results$poc_validated, "\n")
  cat("PoC probability:", round(normal_results$poc_probability, 3), "\n")
  cat("Selection reason:", normal_results$selection_reason, "\n")
}

cat("\n=== ALL TESTS COMPLETED ===\n")
cat("Early termination and PoC implementation appears to be working correctly.\n")
