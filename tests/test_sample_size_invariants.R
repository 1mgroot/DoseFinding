# Test script to verify sample-size invariants for group-sequential design
# 
# Expected invariants (cohort_size=15, stages=5):
# 1. If trial completes: N = exactly 75
# 2. If early stopping: N ∈ {15, 30, 45, 60} (multiple of 15, NOT exceeding 75)
# 3. Termination stage matches when A becomes empty

library(dplyr)

# Set working directory to project root
if (basename(getwd()) == "tests") {
  setwd("..")
}

# Source all required files
source("src/utils/helpers.R")
source("src/core/simulate_data.R")
source("src/core/model_utils.R")
source("src/decision/dose_decision.R")
source("src/core/main.R")

cat("=== Testing Sample Size Invariants ===\n\n")

# Test configuration
test_config <- list(
  dose_levels = c(1, 2, 3, 4, 5),
  n_stages = 5,
  cohort_size = 15,
  phi_T = 0.30,
  c_T = 0.5,
  phi_E = 0.25,
  c_E = 0.5,
  phi_I = 0.20,
  c_I = 0.5,
  delta_poc = 0.8,
  c_poc = 0.7,
  enable_early_termination = TRUE,
  log_early_termination = FALSE,
  verbose_logging = FALSE
)

# Add utility table
utility_table <- array(0, dim = c(2, 2, 2))
utility_table[1, 1, 1] <- 0
utility_table[2, 1, 1] <- 80
utility_table[1, 2, 1] <- 0
utility_table[2, 2, 1] <- 30
utility_table[1, 1, 2] <- 10
utility_table[2, 1, 2] <- 100
utility_table[1, 2, 2] <- 0
utility_table[2, 2, 2] <- 40
test_config$utility_table <- utility_table

# Scenario with moderate probabilities
test_scenario <- list(
  p_YI = c(0.1, 0.2, 0.3, 0.4, 0.5),
  p_YT_given_I = matrix(c(0.15, 0.10, 0.10, 0.05, 0.05,
                          0.20, 0.15, 0.10, 0.05, 0.05), 
                        nrow = 5, ncol = 2),
  p_YE_given_I = matrix(c(0.20, 0.25, 0.30, 0.35, 0.40,
                          0.30, 0.35, 0.40, 0.45, 0.50),
                        nrow = 5, ncol = 2),
  rho0 = 1.5,
  rho1 = 2.0
)

cat("Running 50 test simulations...\n")
cat("Checking sample size invariants:\n")
cat("1. Completed trials: N = 75\n")
cat("2. Early terminated: N ∈ {15, 30, 45, 60}\n")
cat("3. N ≤ 75 (never exceed maximum)\n\n")

n_tests <- 50
issues <- list()
all_results <- list()

for (i in 1:n_tests) {
  result <- run_trial_simulation(
    test_config,
    test_scenario$p_YI,
    test_scenario$p_YT_given_I,
    test_scenario$p_YE_given_I,
    test_scenario$rho0,
    test_scenario$rho1,
    seed = 50000 + i
  )
  
  N <- nrow(result$all_data)
  all_results[[i]] <- list(
    seed = 50000 + i,
    terminated_early = result$terminated_early,
    termination_stage = result$termination_stage,
    N = N
  )
  
  # Check invariants
  if (result$terminated_early) {
    # Early termination: N should be multiple of 15 and ≤ 60
    expected_N <- result$termination_stage * test_config$cohort_size
    
    if (N != expected_N) {
      issues[[length(issues) + 1]] <- list(
        sim = i,
        seed = 50000 + i,
        issue = "Sample size mismatch",
        expected_N = expected_N,
        actual_N = N,
        termination_stage = result$termination_stage
      )
    }
    
    if (N %% test_config$cohort_size != 0) {
      issues[[length(issues) + 1]] <- list(
        sim = i,
        seed = 50000 + i,
        issue = "N not multiple of cohort_size",
        N = N,
        cohort_size = test_config$cohort_size
      )
    }
    
    if (N > test_config$n_stages * test_config$cohort_size) {
      issues[[length(issues) + 1]] <- list(
        sim = i,
        seed = 50000 + i,
        issue = "N exceeds maximum",
        N = N,
        max_N = test_config$n_stages * test_config$cohort_size
      )
    }
  } else {
    # Completed: N should be exactly 75
    expected_N <- test_config$n_stages * test_config$cohort_size
    
    if (N != expected_N) {
      issues[[length(issues) + 1]] <- list(
        sim = i,
        seed = 50000 + i,
        issue = "Completed trial N mismatch",
        expected_N = expected_N,
        actual_N = N
      )
    }
  }
}

# Summary statistics
cat("=== Summary Statistics ===\n")
N_values <- sapply(all_results, function(x) x$N)
early_term_count <- sum(sapply(all_results, function(x) x$terminated_early))

cat("Total simulations:", n_tests, "\n")
cat("Early terminations:", early_term_count, "\n")
cat("Completed trials:", n_tests - early_term_count, "\n")
cat("\nSample sizes observed:\n")
print(table(N_values))

cat("\n=== Invariant Check Results ===\n")
if (length(issues) == 0) {
  cat("✓ ALL CHECKS PASSED - No sample size issues detected!\n")
} else {
  cat("✗ ISSUES DETECTED:", length(issues), "violations\n\n")
  
  for (issue in issues) {
    cat("Simulation", issue$sim, "(seed:", issue$seed, ")\n")
    cat("  Issue:", issue$issue, "\n")
    if (!is.null(issue$expected_N)) {
      cat("  Expected N:", issue$expected_N, "| Actual N:", issue$actual_N, "\n")
    }
    if (!is.null(issue$termination_stage)) {
      cat("  Termination stage:", issue$termination_stage, "\n")
    }
    if (!is.null(issue$max_N)) {
      cat("  Max allowed:", issue$max_N, "| Actual:", issue$N, "\n")
    }
    cat("\n")
  }
  
  # Detailed inspection of first issue
  if (length(issues) > 0) {
    first_issue <- issues[[1]]
    cat("=== Detailed Investigation of First Issue ===\n")
    cat("Re-running simulation", first_issue$sim, "with verbose output...\n\n")
    
    test_config$verbose_logging <- TRUE
    test_config$log_early_termination <- TRUE
    
    result_verbose <- run_trial_simulation(
      test_config,
      test_scenario$p_YI,
      test_scenario$p_YT_given_I,
      test_scenario$p_YE_given_I,
      test_scenario$rho0,
      test_scenario$rho1,
      seed = first_issue$seed
    )
    
    cat("\nStage-by-stage enrollment:\n")
    stage_summary <- result_verbose$all_data %>%
      group_by(stage) %>%
      summarise(n_patients = n(), .groups = 'drop')
    print(stage_summary)
    
    cat("\nTotal enrolled:", nrow(result_verbose$all_data), "\n")
    cat("Expected:", first_issue$expected_N, "\n")
  }
}

cat("\n=== Test Complete ===\n")
