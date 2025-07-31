library(dplyr)
library(tidyr)
library(isotone)
library(purrr)
library(ggplot2)
library(Iso)

source("/Users/jz/Development/DoseFinding/bayesTrial_refactored/config.R")
source("/Users/jz/Development/DoseFinding/bayesTrial_refactored/helpers.R")
source("/Users/jz/Development/DoseFinding/bayesTrial_refactored/simulate_data.R")
source("/Users/jz/Development/DoseFinding/bayesTrial_refactored/model_utils.R")
source("/Users/jz/Development/DoseFinding/bayesTrial_refactored/dose_decision.R")

run_trial_simulation <- function(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1) {
  all_data <- data.frame()
  all_alloc_probs <- data.frame()
  
  # Initial allocation is uniform
  alloc_probs <- rep(1/length(trial_config$dose_levels), length(trial_config$dose_levels))

  for (stage in 1:trial_config$n_stages) {
    cat(paste("
--- Stage", stage, "---
"))
    
    n_next_stage <- round(alloc_probs * trial_config$cohort_size)

    stage_data <- simulate_data_gumbel(
      n_per_dose_vector = n_next_stage,
      dose_levels = trial_config$dose_levels,
      p_YI = p_YI,
      p_YT_given_I = p_YT_given_I,
      p_YE_given_I = p_YE_given_I,
      rho0 = rho0,
      rho1 = rho1
    )
    stage_data$stage <- stage
    all_data <- rbind(all_data, stage_data)

    dose_grid <- expand.grid(d = trial_config$dose_levels)
    pi_I_stats <- compute_rn(all_data, outcome_col = "Y_I")
    pi_I_stats_complete <- left_join(dose_grid, pi_I_stats, by = "d") %>%
      replace_na(list(r = 0, n = 0))
    pi_I_post <- simulate_beta_posterior(pi_I_stats_complete)
    pi_I_post <- add_beta_variance(pi_I_post)
    pi_I_pava <- apply_pava_on_samples(pi_I_post)

    dose_group_grid <- expand.grid(d = trial_config$dose_levels, Y_I = 0:1)
    tox_stats <- compute_rn(all_data, outcome_col = "Y_T", group_col = "Y_I")
    tox_stats_complete <- left_join(dose_group_grid, tox_stats, by = c("d", "Y_I")) %>%
      replace_na(list(r = 0, n = 0))
    tox_post <- simulate_beta_posterior(tox_stats_complete)
    tox_post <- add_beta_variance(tox_post)
    tox_pava <- apply_biviso_on_matrix(tox_post)

    eff_stats <- compute_rn(all_data, outcome_col = "Y_E", group_col = "Y_I")
    eff_stats_complete <- left_join(dose_group_grid, eff_stats, by = c("d", "Y_I")) %>%
      replace_na(list(r = 0, n = 0))
    eff_post <- simulate_beta_posterior(eff_stats_complete)
    eff_post <- add_beta_variance(eff_post)
    eff_pava <- apply_biviso_on_matrix(eff_post)

    tox_marginal <- compute_marginal_probability(tox_pava, pi_I_pava)
    eff_marginal <- compute_marginal_probability(eff_pava, pi_I_pava)
    posterior_summaries <- list(
      tox = tox_pava, eff = eff_pava, imm = pi_I_pava,
      tox_marginal = tox_marginal, eff_marginal = eff_marginal
    )

    admissible_set <- get_admissible_set(posterior_summaries, trial_config)
    cat("Admissible set:", admissible_set, "
")

    log_utility_calculations(posterior_summaries, trial_config)

    if (stage < trial_config$n_stages) {
      alloc_probs <- adaptive_randomization(admissible_set, posterior_summaries, trial_config)
      cat("Allocation probabilities for next stage:", alloc_probs, "
")
    }
    
    all_alloc_probs <- rbind(all_alloc_probs, data.frame(Stage = stage, Dose = trial_config$dose_levels, Prob = alloc_probs))
  }

  final_od <- select_final_od(admissible_set, posterior_summaries, trial_config)

  return(list(
    final_od = final_od,
    all_data = all_data,
    all_alloc_probs = all_alloc_probs,
    posterior_summaries = posterior_summaries
  ))
}

# This part will only run when the script is sourced directly, for testing
if (sys.nframe() == 0) {
  results <- run_trial_simulation(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1)
  
  cat("
--- Final Results ---
")
  cat("Final OD:", results$final_od, "
")

  plot_posterior_summary(results$posterior_summaries$imm, title = "Immune Response vs Dose (PAVA Adjusted)", file_path = "results/immune_response_refactored.png")
  plot_posterior_summary(results$posterior_summaries$tox, title = "Toxicity Rate by Dose and Immune Status", group_col = "Y_I", file_path = "results/toxicity_refactored.png")
  plot_posterior_summary(results$posterior_summaries$eff, title = "Efficacy Rate by Dose and Immune Status", group_col = "Y_I", file_path = "results/efficacy_refactored.png")
}