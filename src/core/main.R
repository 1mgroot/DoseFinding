library(dplyr)
library(tidyr)
library(isotone)
library(purrr)
library(ggplot2)
library(Iso)

# Source files relative to the project root.
# NOTE: config.R should be sourced only once at the top level (e.g., notebook or main script).
project_root_candidates <- c(".", "..", "../..")
project_root_matches <- project_root_candidates[
  file.exists(file.path(project_root_candidates, "DoseFinding.Rproj"))
]
if (length(project_root_matches) == 0) {
  stop("Could not find project root containing DoseFinding.Rproj.")
}
project_root <- project_root_matches[[1]]
source(file.path(project_root, "src/core/simulate_data.R"))
source(file.path(project_root, "src/core/model_utils.R"))
source(file.path(project_root, "src/utils/helpers.R"))
source(file.path(project_root, "src/decision/dose_decision.R"))

validate_rng_seed <- function(seed, name = "seed") {
  if (is.null(seed)) {
    return(NULL)
  }
  max_seed <- .Machine$integer.max - 1L
  if (!is.numeric(seed) ||
      length(seed) != 1 ||
      is.na(seed) ||
      seed < 0 ||
      seed != floor(seed) ||
      seed > max_seed) {
    stop(name, " must be a non-negative integer no larger than ", max_seed, ".")
  }
  as.integer(seed)
}

generate_trial_stage_seeds <- function(seed, n_stages) {
  if (!is.numeric(n_stages) ||
      length(n_stages) != 1 ||
      is.na(n_stages) ||
      n_stages < 1 ||
      n_stages != floor(n_stages)) {
    stop("n_stages must be a positive integer.")
  }
  if (is.null(seed)) {
    return(rep(NA_integer_, n_stages))
  }

  set.seed(validate_rng_seed(seed))
  sample.int(.Machine$integer.max - 1L, n_stages)
}

dose_labels_from_indices <- function(dose_levels, dose_indices) {
  if (length(dose_indices) == 0) {
    return(dose_levels[integer(0)])
  }
  if (any(!is.na(dose_indices) &
          (dose_indices < 1 | dose_indices > length(dose_levels)))) {
    stop("dose_indices contains an index outside dose_levels.")
  }
  unname(dose_levels[dose_indices])
}

name_values_by_dose_label <- function(values, dose_labels) {
  if (length(values) != length(dose_labels)) {
    stop("values and dose_labels must have the same length.")
  }
  names(values) <- as.character(dose_labels)
  values
}

run_trial_simulation <- function(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1, seed = NULL) {
  all_data <- data.frame()
  all_alloc_probs <- data.frame()
  stage_seeds <- generate_trial_stage_seeds(seed, trial_config$n_stages)
  
  # Get verbose logging flag (default TRUE for backward compatibility)
  verbose <- if (is.null(trial_config$verbose_logging)) TRUE else trial_config$verbose_logging
  
  # Initial allocation is uniform
  alloc_probs <- rep(1/length(trial_config$dose_levels), length(trial_config$dose_levels))

  for (stage in 1:trial_config$n_stages) {
    if (verbose) {
      cat(paste("
--- Stage", stage, "---
"))
      if (stage == 1) {
        cat("Workflow: Step 1 - Equal randomization to all dose levels
")
      } else {
        cat("Workflow: Step 1 - Adaptive randomization (using probabilities from previous stage)
")
      }
    }
    
    # Record allocation probabilities at START of stage (before any early termination)
    all_alloc_probs <- rbind(all_alloc_probs, data.frame(Stage = stage, Dose = trial_config$dose_levels, Prob = alloc_probs))
    
    # Stage 1: Explicitly ensure equal allocation across all dose levels
    # Stages 2+: Use adaptive randomization probabilities
    if (stage == 1) {
      # Equal allocation: divide cohort_size evenly across all dose levels
      n_per_dose <- floor(trial_config$cohort_size / length(trial_config$dose_levels))
      remainder <- trial_config$cohort_size - (n_per_dose * length(trial_config$dose_levels))
      n_next_stage <- rep(n_per_dose, length(trial_config$dose_levels))
      # Distribute remainder across first few doses
      if (remainder > 0) {
        n_next_stage[1:remainder] <- n_next_stage[1:remainder] + 1
      }
    } else {
      # Adaptive allocation: ensure sum equals cohort_size exactly
      # Use floor + distribute remainder to avoid rounding errors
      n_exact <- alloc_probs * trial_config$cohort_size
      n_next_stage <- floor(n_exact)
      remainder <- trial_config$cohort_size - sum(n_next_stage)
      
      # Distribute remainder to doses with largest fractional parts
      if (remainder > 0) {
        fractional <- n_exact - n_next_stage
        top_indices <- order(fractional, decreasing = TRUE)[1:remainder]
        n_next_stage[top_indices] <- n_next_stage[top_indices] + 1
      }
    }
    
    # INVARIANT CHECK: Ensure cohort size sums correctly
    if (sum(n_next_stage) != trial_config$cohort_size) {
      stop(paste("Internal error: Stage", stage, "allocation sums to", sum(n_next_stage),
                 "but cohort_size is", trial_config$cohort_size))
    }

    # Generate stage-specific seed if base seed is provided.
    stage_seed <- if (!is.null(seed)) stage_seeds[[stage]] else NULL
    
    stage_data <- simulate_data_gumbel(
      n_per_dose_vector = n_next_stage,
      dose_levels = trial_config$dose_levels,
      p_YI = p_YI,
      p_YT_given_I = p_YT_given_I,
      p_YE_given_I = p_YE_given_I,
      rho0 = rho0,
      rho1 = rho1,
      seed = stage_seed,
      id_start = nrow(all_data) + 1L
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

    if (verbose) {
      cat("Workflow: Step 2 - Interim Analysis (update admissible set based on posterior probabilities)
")
    }
    admissible_set <- get_admissible_set(posterior_summaries, trial_config, verbose)
    if (verbose) {
      cat("Admissible set:", admissible_set, "
")
    }

    # Step 4: Early Termination Check (CORRECT PLACEMENT - after interim analysis, before adaptive randomization)
    # GROUP-SEQUENTIAL DESIGN: Check if A is empty after interim analysis of stage k
    # If A is empty, stop BEFORE enrolling stage k+1
    # Sample size at termination = stage * cohort_size (already enrolled through current stage)
    if (verbose) {
      cat("Workflow: Step 4 - Early Termination Check
")
    }
    if (check_early_termination(admissible_set, trial_config)) {
      termination_info <- handle_trial_termination(admissible_set, stage, trial_config)
      
      # INVARIANT VALIDATION: Check sample size at termination
      N_enrolled <- nrow(all_data)
      expected_N <- stage * trial_config$cohort_size
      if (N_enrolled != expected_N) {
        warning(paste("Sample size mismatch at early termination:",
                     "Expected", expected_N, "but enrolled", N_enrolled))
      }
      
      return(list(
        final_od = dose_labels_from_indices(trial_config$dose_levels, NA_integer_),
        final_od_index = NA_integer_,
        final_utility = NA_real_,
        poc_validated = FALSE,
        poc_probability = 0,
        selection_reason = termination_info$reason,
        final_admissible_set = dose_labels_from_indices(trial_config$dose_levels, admissible_set),
        final_admissible_indices = admissible_set,
        poc_eligible_set = trial_config$dose_levels[integer(0)],
        poc_eligible_indices = integer(0),
        poc_pairwise_probs = numeric(0),
        final_candidate_utilities = numeric(0),
        all_data = all_data,
        all_alloc_probs = all_alloc_probs,
        posterior_summaries = posterior_summaries,
        terminated_early = TRUE,
        termination_stage = stage,
        termination_reason = termination_info$reason
      ))
    }

    if (verbose) {
      log_utility_calculations(posterior_summaries, trial_config)
    }

    # Step 3: Adaptive Randomization (only if trial continues)
    if (stage < trial_config$n_stages) {
      if (verbose) {
        cat("Workflow: Step 3 - Adaptive Randomization (allocate by posterior probability of being optimal)
")
      }
      alloc_probs <- adaptive_randomization(admissible_set, posterior_summaries, trial_config)
      if (verbose) {
        cat("Allocation probabilities for next stage:", alloc_probs, "
")
      }
    }
  }

  # Step 5: Final Selection with PoC validation (CORRECT PLACEMENT - only at final stage)
  if (verbose) {
    cat("Workflow: Step 5 - Final Selection with PoC validation
")
  }
  final_selection <- select_final_od_with_poc(admissible_set, posterior_summaries, trial_config, verbose)

  # INVARIANT VALIDATION: Check sample size for completed trial
  N_enrolled <- nrow(all_data)
  expected_N <- trial_config$n_stages * trial_config$cohort_size
  if (N_enrolled != expected_N) {
    warning(paste("Sample size mismatch for completed trial:",
                 "Expected", expected_N, "but enrolled", N_enrolled))
  }

  final_od_index <- final_selection$optimal_dose
  final_admissible_indices <- final_selection$admissible_doses
  poc_eligible_indices <- final_selection$P_final
  final_admissible_labels <- dose_labels_from_indices(
    trial_config$dose_levels,
    final_admissible_indices
  )
  poc_eligible_labels <- dose_labels_from_indices(
    trial_config$dose_levels,
    poc_eligible_indices
  )

  return(list(
    final_od = dose_labels_from_indices(trial_config$dose_levels, final_od_index),
    final_od_index = final_od_index,
    final_utility = final_selection$optimal_utility,
    poc_validated = final_selection$poc_validated,
    poc_probability = final_selection$poc_probability,
    selection_reason = final_selection$reason,
    final_admissible_set = final_admissible_labels,
    final_admissible_indices = final_admissible_indices,
    poc_eligible_set = poc_eligible_labels,
    poc_eligible_indices = poc_eligible_indices,
    poc_pairwise_probs = name_values_by_dose_label(
      final_selection$pairwise_probs,
      final_admissible_labels
    ),
    final_candidate_utilities = name_values_by_dose_label(
      final_selection$utilities,
      final_admissible_labels
    ),
    all_data = all_data,
    all_alloc_probs = all_alloc_probs,
    posterior_summaries = posterior_summaries,
    terminated_early = FALSE,
    termination_stage = NA,
    termination_reason = NA
  ))
}

# This part will only run when the script is sourced directly, for testing
if (sys.nframe() == 0) {
  results <- run_trial_simulation(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1)
  
  cat("
--- Final Results ---
")
  
  if (results$terminated_early) {
    cat("Trial terminated early at stage:", results$termination_stage, "
")
    cat("Reason:", results$termination_reason, "
")
    cat("No Optimal Dose selected
")
  } else {
    cat("Final OD:", results$final_od, "
")
    cat("Final utility:", round(results$final_utility, 2), "
")
    cat("PoC validated:", results$poc_validated, "
")
    cat("PoC probability:", round(results$poc_probability, 3), "
")
    cat("Selection reason:", results$selection_reason, "
")
  }

  plot_posterior_summary(results$posterior_summaries$imm, title = "Immune Response vs Dose (PAVA Adjusted)", file_path = "results/immune_response_refactored.png", style = "modern")
  plot_posterior_summary(results$posterior_summaries$tox, title = "Toxicity Rate by Dose and Immune Status", group_col = "Y_I", file_path = "results/toxicity_refactored.png", style = "modern")
  plot_posterior_summary(results$posterior_summaries$eff, title = "Efficacy Rate by Dose and Immune Status", group_col = "Y_I", file_path = "results/efficacy_refactored.png", style = "modern")
}
