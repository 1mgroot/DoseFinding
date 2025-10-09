get_expected_utility <- function(dose_idx, posterior_summaries, config) {
  # Get posterior probabilities for the given dose
  # Note: This assumes the posterior_summaries dataframes are ordered by dose and then by group
  pi_T_given_I0 <- posterior_summaries$tox$pava_mean[2 * dose_idx - 1]
  pi_T_given_I1 <- posterior_summaries$tox$pava_mean[2 * dose_idx]
  pi_E_given_I0 <- posterior_summaries$eff$pava_mean[2 * dose_idx - 1]
  pi_E_given_I1 <- posterior_summaries$eff$pava_mean[2 * dose_idx]
  pi_I <- posterior_summaries$imm$pava_mean[dose_idx]

  # Probabilities of T=0 and T=1 given I
  p_T_given_I0 <- c(1 - pi_T_given_I0, pi_T_given_I0)
  p_T_given_I1 <- c(1 - pi_T_given_I1, pi_T_given_I1)

  # Probabilities of E=0 and E=1 given I
  p_E_given_I0 <- c(1 - pi_E_given_I0, pi_E_given_I0)
  p_E_given_I1 <- c(1 - pi_E_given_I1, pi_E_given_I1)

  # Expected utility for I=0
  utility_I0 <- sum(config$utility_table[,,1] * (p_E_given_I0 %o% p_T_given_I0))

  # Expected utility for I=1
  utility_I1 <- sum(config$utility_table[,,2] * (p_E_given_I1 %o% p_T_given_I1))

  # Total expected utility
  total_utility <- (1 - pi_I) * utility_I0 + pi_I * utility_I1
  return(total_utility)
}

get_expected_utility_detailed <- function(dose_idx, posterior_summaries, config) {
  # Get posterior probabilities for the given dose
  pi_T_given_I0 <- posterior_summaries$tox$pava_mean[2 * dose_idx - 1]
  pi_T_given_I1 <- posterior_summaries$tox$pava_mean[2 * dose_idx]
  pi_E_given_I0 <- posterior_summaries$eff$pava_mean[2 * dose_idx - 1]
  pi_E_given_I1 <- posterior_summaries$eff$pava_mean[2 * dose_idx]
  pi_I <- posterior_summaries$imm$pava_mean[dose_idx]

  # Probabilities of T=0 and T=1 given I
  p_T_given_I0 <- c(1 - pi_T_given_I0, pi_T_given_I0)
  p_T_given_I1 <- c(1 - pi_T_given_I1, pi_T_given_I1)

  # Probabilities of E=0 and E=1 given I
  p_E_given_I0 <- c(1 - pi_E_given_I0, pi_E_given_I0)
  p_E_given_I1 <- c(1 - pi_E_given_I1, pi_E_given_I1)

  # Expected utility for I=0
  utility_I0 <- sum(config$utility_table[,,1] * (p_E_given_I0 %o% p_T_given_I0))

  # Expected utility for I=1
  utility_I1 <- sum(config$utility_table[,,2] * (p_E_given_I1 %o% p_T_given_I1))

  # Total expected utility
  total_utility <- (1 - pi_I) * utility_I0 + pi_I * utility_I1
  
  # Return detailed breakdown
  return(list(
    dose_idx = dose_idx,
    pi_I = pi_I,
    pi_T_given_I0 = pi_T_given_I0,
    pi_T_given_I1 = pi_T_given_I1,
    pi_E_given_I0 = pi_E_given_I0,
    pi_E_given_I1 = pi_E_given_I1,
    p_T_given_I0 = p_T_given_I0,
    p_T_given_I1 = p_T_given_I1,
    p_E_given_I0 = p_E_given_I0,
    p_E_given_I1 = p_E_given_I1,
    utility_I0 = utility_I0,
    utility_I1 = utility_I1,
    total_utility = total_utility
  ))
}

get_admissible_set <- function(posterior_summaries, config) {
  admissible_doses <- c()
  cat("\n--- Admissibility Check ---\n")
  
  # Log summary statistics for transparency
  cat("Summary: Toxicity marginal means:", round(posterior_summaries$tox_marginal$marginal_prob, 3), "\n")
  cat("Summary: Efficacy marginal means:", round(posterior_summaries$eff_marginal$marginal_prob, 3), "\n")
  cat("Summary: Immune response means:", round(posterior_summaries$imm$pava_mean, 3), "\n")
  
  for (i in 1:length(config$dose_levels)) {
    tox_samples <- posterior_summaries$tox_marginal$samples[[i]]
    eff_samples <- posterior_summaries$eff_marginal$samples[[i]]
    imm_samples <- posterior_summaries$imm$samples_pava[[i]]

    tox_prob_safe <- mean(tox_samples < config$phi_T)
    eff_prob_good <- mean(eff_samples > config$phi_E)
    imm_prob_good <- mean(imm_samples > config$phi_I)

    # Log detailed statistics for each dose
    cat(paste("Dose", i,":",
              "P(Tox < ", config$phi_T, ") = ", round(tox_prob_safe, 2),
              "(Threshold: ", config$c_T, ")",
              "P(Eff > ", config$phi_E, ") = ", round(eff_prob_good, 2),
              "(Threshold: ", config$c_E, ")",
              "P(Imm > ", config$phi_I, ") = ", round(imm_prob_good, 2),
              "(Threshold: ", config$c_I, ")\n"))

    if (tox_prob_safe > config$c_T && eff_prob_good > config$c_E && imm_prob_good > config$c_I) {
      admissible_doses <- c(admissible_doses, i)
    }
  }
  cat("--- End Admissibility Check ---\n")
  return(admissible_doses)
}

adaptive_randomization <- function(admissible_set, posterior_summaries, config) {
  n_doses <- length(config$dose_levels)
  alloc_probs <- numeric(n_doses)

  if (length(admissible_set) > 0) {
    utilities <- sapply(admissible_set, get_expected_utility, posterior_summaries, config)
    if (sum(utilities) > 0) {
      alloc_probs[admissible_set] <- utilities / sum(utilities)
    } else {
      alloc_probs[admissible_set] <- 1 / length(admissible_set) # Equal probability if all utilities are zero
    }
  }

  return(alloc_probs)
}

select_final_od <- function(admissible_set, posterior_summaries, config) {
  if (length(admissible_set) == 0) {
    return(NA)
  }

  utilities <- sapply(admissible_set, get_expected_utility, posterior_summaries, config)
  return(admissible_set[which.max(utilities)])
}

# New function to log utility calculations
log_utility_calculations <- function(posterior_summaries, config) {
  cat("\n--- Utility Score Calculations ---\n")
  
  # Display utility table for reference
  cat("Utility Table Reference:\n")
  cat("  I=0 (No Immune Response):\n")
  cat("    E=0, T=0:", config$utility_table[1,1,1], "  E=1, T=0:", config$utility_table[2,1,1], "\n")
  cat("    E=0, T=1:", config$utility_table[1,2,1], "  E=1, T=1:", config$utility_table[2,2,1], "\n")
  cat("  I=1 (Immune Response):\n")
  cat("    E=0, T=0:", config$utility_table[1,1,2], "  E=1, T=0:", config$utility_table[2,1,2], "\n")
  cat("    E=0, T=1:", config$utility_table[1,2,2], "  E=1, T=1:", config$utility_table[2,2,2], "\n")
  cat("\n")
  
  # Calculate utilities for all doses
  utility_summary <- data.frame(
    Dose = 1:length(config$dose_levels),
    Immune_Prob = numeric(length(config$dose_levels)),
    Tox_I0 = numeric(length(config$dose_levels)),
    Tox_I1 = numeric(length(config$dose_levels)),
    Eff_I0 = numeric(length(config$dose_levels)),
    Eff_I1 = numeric(length(config$dose_levels)),
    Utility_I0 = numeric(length(config$dose_levels)),
    Utility_I1 = numeric(length(config$dose_levels)),
    Total_Utility = numeric(length(config$dose_levels))
  )
  
  for (i in 1:length(config$dose_levels)) {
    utility_details <- get_expected_utility_detailed(i, posterior_summaries, config)
    
    # Store summary data
    utility_summary$Immune_Prob[i] <- utility_details$pi_I
    utility_summary$Tox_I0[i] <- utility_details$pi_T_given_I0
    utility_summary$Tox_I1[i] <- utility_details$pi_T_given_I1
    utility_summary$Eff_I0[i] <- utility_details$pi_E_given_I0
    utility_summary$Eff_I1[i] <- utility_details$pi_E_given_I1
    utility_summary$Utility_I0[i] <- utility_details$utility_I0
    utility_summary$Utility_I1[i] <- utility_details$utility_I1
    utility_summary$Total_Utility[i] <- utility_details$total_utility
    
    cat(paste("Dose", i, "Utility Calculation:\n"))
    cat(paste("  Immune response probability (π_I):", round(utility_details$pi_I, 3), "\n"))
    cat(paste("  Toxicity given I=0 (π_T|I=0):", round(utility_details$pi_T_given_I0, 3), "\n"))
    cat(paste("  Toxicity given I=1 (π_T|I=1):", round(utility_details$pi_T_given_I1, 3), "\n"))
    cat(paste("  Efficacy given I=0 (π_E|I=0):", round(utility_details$pi_E_given_I0, 3), "\n"))
    cat(paste("  Efficacy given I=1 (π_E|I=1):", round(utility_details$pi_E_given_I1, 3), "\n"))
    
    cat("  Probability distributions:\n")
    cat(paste("    P(T=0|I=0):", round(utility_details$p_T_given_I0[1], 3), 
              "P(T=1|I=0):", round(utility_details$p_T_given_I0[2], 3), "\n"))
    cat(paste("    P(T=0|I=1):", round(utility_details$p_T_given_I1[1], 3), 
              "P(T=1|I=1):", round(utility_details$p_T_given_I1[2], 3), "\n"))
    cat(paste("    P(E=0|I=0):", round(utility_details$p_E_given_I0[1], 3), 
              "P(E=1|I=0):", round(utility_details$p_E_given_I0[2], 3), "\n"))
    cat(paste("    P(E=0|I=1):", round(utility_details$p_E_given_I1[1], 3), 
              "P(E=1|I=1):", round(utility_details$p_E_given_I1[2], 3), "\n"))
    
    cat(paste("  Expected utility given I=0:", round(utility_details$utility_I0, 2), "\n"))
    cat(paste("  Expected utility given I=1:", round(utility_details$utility_I1, 2), "\n"))
    cat(paste("  Total expected utility:", round(utility_details$total_utility, 2), "\n"))
    cat("\n")
  }
  
  # Display utility summary table
  cat("Utility Summary Table:\n")
  cat("Dose | Immune | Tox(I=0) | Tox(I=1) | Eff(I=0) | Eff(I=1) | U(I=0) | U(I=1) | Total U\n")
  cat("-----|--------|----------|----------|----------|----------|--------|--------|--------\n")
  for (i in 1:nrow(utility_summary)) {
    cat(sprintf("%4d | %6.3f | %8.3f | %8.3f | %8.3f | %8.3f | %6.1f | %6.1f | %7.1f\n",
                utility_summary$Dose[i],
                utility_summary$Immune_Prob[i],
                utility_summary$Tox_I0[i],
                utility_summary$Tox_I1[i],
                utility_summary$Eff_I0[i],
                utility_summary$Eff_I1[i],
                utility_summary$Utility_I0[i],
                utility_summary$Utility_I1[i],
                utility_summary$Total_Utility[i]))
  }
  cat("\n")
  
  cat("--- End Utility Calculations ---\n")
}

# Early Termination Logic Functions

check_early_termination <- function(admissible_set, config) {
  # Check if trial should terminate early due to empty admissible set.
  #
  # Args:
  #   admissible_set: Vector of admissible dose indices
  #   config: Trial configuration
  #
  # Returns:
  #   logical: TRUE if trial should terminate early
  if (!config$enable_early_termination) {
    return(FALSE)
  }
  
  should_terminate <- length(admissible_set) == 0
  
  if (should_terminate && config$log_early_termination) {
    cat("\n--- EARLY TERMINATION TRIGGERED ---\n")
    cat("Reason: Admissible set is empty (no doses meet safety/efficacy criteria)\n")
    cat("Trial will terminate without selecting an Optimal Dose\n")
    cat("--- END EARLY TERMINATION ---\n\n")
  }
  
  return(should_terminate)
}

handle_trial_termination <- function(admissible_set, stage, config) {
  # Handle early trial termination.
  #
  # Args:
  #   admissible_set: Vector of admissible dose indices
  #   stage: Current trial stage
  #   config: Trial configuration
  #
  # Returns:
  #   list: Termination information
  termination_info <- list(
    terminated_early = TRUE,
    stage = stage,
    admissible_set = admissible_set,
    optimal_dose = NA,
    reason = "Empty admissible set"
  )
  
  if (config$log_early_termination) {
    cat("\n=== TRIAL TERMINATION SUMMARY ===\n")
    cat("Trial terminated early at stage:", stage, "\n")
    cat("Reason:", termination_info$reason, "\n")
    cat("No Optimal Dose selected\n")
    cat("=== END TRIAL TERMINATION ===\n\n")
  }
  
  return(termination_info)
}

# Probability of Correct Selection (PoC) Functions

calculate_pi_parameters <- function(dose_idx, posterior_summaries) {
  # Calculate Πᵢ parameters (combined efficacy measure) for a given dose.
  #
  # Args:
  #   dose_idx: Dose index
  #   posterior_summaries: Posterior probability summaries
  #
  # Returns:
  #   list: Πᵢ samples and summary statistics
  
  # Get posterior samples for this dose
  pi_I_samples <- posterior_summaries$imm$samples_pava[[dose_idx]]
  pi_E_given_I0_samples <- posterior_summaries$eff$samples[[2 * dose_idx - 1]]
  pi_E_given_I1_samples <- posterior_summaries$eff$samples[[2 * dose_idx]]
  
  # Calculate Πᵢ samples (combined efficacy measure)
  # Πᵢ = P(E|I=0) * P(I=0) + P(E|I=1) * P(I=1)
  pi_combined_samples <- pi_I_samples * pi_E_given_I1_samples + 
                        (1 - pi_I_samples) * pi_E_given_I0_samples
  
  return(list(
    pi_I_samples = pi_I_samples,
    pi_E_given_I0_samples = pi_E_given_I0_samples,
    pi_E_given_I1_samples = pi_E_given_I1_samples,
    pi_combined_samples = pi_combined_samples,
    pi_combined_mean = mean(pi_combined_samples),
    pi_combined_sd = sd(pi_combined_samples)
  ))
}

calculate_poc_probability <- function(admissible_set, posterior_summaries, config) {
  # Calculate Probability of Correct Selection (PoC) for admissible doses using proper Bayesian approach.
  #
  # Args:
  #   admissible_set: Vector of admissible dose indices
  #   posterior_summaries: Posterior probability summaries
  #   config: Trial configuration
  #
  # Returns:
  #   list: PoC probabilities for each admissible dose
  if (length(admissible_set) == 0) {
    return(list(poc_probabilities = numeric(0), max_poc = 0))
  }
  
  poc_probabilities <- numeric(length(admissible_set))
  
  # Calculate utilities to determine the best dose (reference)
  utilities <- sapply(admissible_set, get_expected_utility, posterior_summaries, config)
  best_dose_idx <- admissible_set[which.max(utilities)]
  
  # Calculate Πᵢⱼ parameters for the best dose (reference)
  best_dose_params <- calculate_pi_parameters(best_dose_idx, posterior_summaries)
  
  for (i in seq_along(admissible_set)) {
    dose_idx <- admissible_set[i]
    
    # Calculate Πᵢ parameters for this dose
    dose_params <- calculate_pi_parameters(dose_idx, posterior_summaries)
    
    # Calculate PoC probability: Pr(Πᵢ < δ Πᵢⱼ | Dₙ)
    # This represents the probability that this dose is significantly worse than the best dose
    # Using proper Bayesian calculation with posterior samples
    poc_prob <- mean(dose_params$pi_combined_samples < config$delta_poc * best_dose_params$pi_combined_samples)
    
    poc_probabilities[i] <- poc_prob
    
    # Log detailed PoC calculation if enabled
    if (config$log_early_termination) {
      cat(sprintf("PoC calculation for dose %d: Πᵢ=%.3f±%.3f, Πᵢⱼ=%.3f±%.3f, PoC=%.3f\n",
                  dose_idx, 
                  dose_params$pi_combined_mean, dose_params$pi_combined_sd,
                  best_dose_params$pi_combined_mean, best_dose_params$pi_combined_sd,
                  poc_prob))
    }
  }
  
  max_poc <- max(poc_probabilities, na.rm = TRUE)
  
  return(list(
    poc_probabilities = poc_probabilities,
    max_poc = max_poc,
    admissible_doses = admissible_set
  ))
}

check_poc_threshold <- function(poc_results, config) {
  # Check if PoC threshold is met for final dose selection.
  #
  # Args:
  #   poc_results: Results from calculate_poc_probability
  #   config: Trial configuration
  #
  # Returns:
  #   logical: TRUE if PoC threshold is met
  if (length(poc_results$poc_probabilities) == 0) {
    return(FALSE)
  }
  
  # Check if the best dose meets the PoC threshold
  # We want the probability of correct selection to be high
  poc_met <- poc_results$max_poc >= config$c_poc
  
  if (config$log_early_termination) {
    cat("\n--- PoC THRESHOLD CHECK ---\n")
    cat("Maximum PoC probability:", round(poc_results$max_poc, 3), "\n")
    cat("PoC threshold:", config$c_poc, "\n")
    cat("PoC threshold met:", poc_met, "\n")
    cat("--- END PoC CHECK ---\n\n")
  }
  
  return(poc_met)
}

# Enhanced Final Dose Selection with PoC
select_final_od_with_poc <- function(admissible_set, posterior_summaries, config) {
  # Select final Optimal Dose with PoC validation.
  #
  # Args:
  #   admissible_set: Vector of admissible dose indices
  #   posterior_summaries: Posterior probability summaries
  #   config: Trial configuration
  #
  # Returns:
  #   list: Selection results including PoC validation
  if (length(admissible_set) == 0) {
    return(list(
      optimal_dose = NA,
      poc_validated = FALSE,
      poc_probability = 0,
      reason = "No admissible doses"
    ))
  }
  
  # Calculate utilities for admissible doses
  utilities <- sapply(admissible_set, get_expected_utility, posterior_summaries, config)
  best_dose_idx <- admissible_set[which.max(utilities)]
  best_utility <- max(utilities)
  
  # Calculate PoC probabilities
  poc_results <- calculate_poc_probability(admissible_set, posterior_summaries, config)
  
  # Check PoC threshold
  poc_validated <- check_poc_threshold(poc_results, config)
  
  # Select optimal dose
  if (poc_validated) {
    optimal_dose <- best_dose_idx
    reason <- "PoC threshold met"
  } else {
    # If PoC not met, still select the best dose but note the issue
    optimal_dose <- best_dose_idx
    reason <- "PoC threshold not met, but selecting best available dose"
  }
  
  selection_result <- list(
    optimal_dose = optimal_dose,
    optimal_utility = best_utility,
    poc_validated = poc_validated,
    poc_probability = poc_results$max_poc,
    admissible_doses = admissible_set,
    utilities = utilities,
    reason = reason
  )
  
  # Log selection results
  if (config$log_early_termination) {
    cat("\n--- FINAL DOSE SELECTION WITH PoC ---\n")
    cat("Admissible doses:", admissible_set, "\n")
    cat("Utilities:", round(utilities, 2), "\n")
    cat("Selected dose:", optimal_dose, "\n")
    cat("Selected utility:", round(best_utility, 2), "\n")
    cat("PoC validated:", poc_validated, "\n")
    cat("Max PoC probability:", round(poc_results$max_poc, 3), "\n")
    cat("Selection reason:", reason, "\n")
    cat("--- END FINAL SELECTION ---\n\n")
  }
  
  return(selection_result)
}