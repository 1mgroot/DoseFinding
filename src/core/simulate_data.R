Gumbel <- function(ttox, teff, c) {
  out <- matrix(0, nrow = 4, ncol = 1)
  out[1, ] <- (1 - ttox) * (1 - teff) + ttox * teff * (1 - ttox) * (1 - teff) * (exp(c) - 1) / (exp(c) + 1)
  out[2, ] <- (1 - ttox) * teff - ttox * teff * (1 - ttox) * (1 - teff) * (exp(c) - 1) / (exp(c) + 1)
  out[3, ] <- ttox * (1 - teff) - ttox * teff * (1 - ttox) * (1 - teff) * (exp(c) - 1) / (exp(c) + 1)
  out[4, ] <- ttox * teff + ttox * teff * (1 - ttox) * (1 - teff) * (exp(c) - 1) / (exp(c) + 1)
  return(out)
}
simulate_data_gumbel <- function(
    n_per_dose_vector = c(10, 10, 10),
    dose_levels = c(1, 2, 3),
    p_YI = c(0.2, 0.5, 0.8), # immune prob per dose
    p_YT_given_I, # marginal tox prob for I=0, I=1
    p_YE_given_I, # marginal eff prob for I=0, I=1
    rho0 = 1, # correlation under I=0
    rho1 = 1, # correlation under I=1
    seed = 123,
    debug = FALSE) {
  # Only set seed if provided (not NULL)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  J <- length(dose_levels)
  d_vec <- rep(dose_levels, times = n_per_dose_vector)
  n_total <- length(d_vec)
  # Generate Gumbel joint distributions for all dose levels
  pi0 <- Gumbel(p_YT_given_I[1], p_YE_given_I[1], rho0)
  pi1 <- Gumbel(p_YT_given_I[2], p_YE_given_I[2], rho1)
  pi0_mat <- matrix(0, nrow = 4, ncol = J)
  pi1_mat <- matrix(0, nrow = 4, ncol = J)
  for (j in seq_len(J)) {
    pi0_mat[, j] <- Gumbel(p_YT_given_I[j, 1], p_YE_given_I[j, 1], rho0)
    pi1_mat[, j] <- Gumbel(p_YT_given_I[j, 2], p_YE_given_I[j, 2], rho1)
  }
  records <- vector("list", n_total)
  for (i in seq_along(d_vec)) {
    d <- d_vec[i]
    I <- rbinom(1, 1, p_YI[d])
    if (I == 0) {
      res <- rmultinom(1, 1, pi0_mat[, d])
    } else {
      res <- rmultinom(1, 1, pi1_mat[, d])
    }
    Y_T <- res[3] + res[4]
    Y_E <- res[2] + res[4]
    records[[i]] <- data.frame(
      id = i,
      d = d,
      Y_I = I,
      Y_T = Y_T,
      Y_E = Y_E
    )
  }
  df <- bind_rows(records)
  if (debug) {
    assign("pi0", pi0_mat, envir = .GlobalEnv)
    assign("pi1", pi1_mat, envir = .GlobalEnv)
    assign("df_sim", df, envir = .GlobalEnv)
  }
  return(df)
}

# Flat Scenario Generation Functions for Calibration

calculate_conditional_efficacy_flat <- function(phi_E_lower, phi_I_lower) {
  # Calculate conditional efficacy probabilities to maintain marginal = phi_E_lower
  # Using total probability formula: P(E) = P(E|I=0)*P(I=0) + P(E|I=1)*P(I=1)
  #
  # Args:
  #   phi_E_lower: Target marginal efficacy rate for all doses
  #   phi_I_lower: Immune response rate for all doses
  #
  # Returns:
  #   matrix: Conditional efficacy probabilities [P(E|I=0), P(E|I=1)]
  
  # Assume P(E|I=1) = phi_E_lower + 0.05 (slightly higher for immune response)
  # Then solve for P(E|I=0) to maintain marginal = phi_E_lower
  
  p_E_given_I1 <- min(phi_E_lower + 0.05, 1.0)
  p_E_given_I0 <- (phi_E_lower - phi_I_lower * p_E_given_I1) / (1 - phi_I_lower)
  p_E_given_I0 <- max(0, min(p_E_given_I0, 1))  # Ensure valid probability
  
  return(matrix(c(p_E_given_I0, p_E_given_I1), nrow = 1))
}

create_flat_probability_matrices <- function(n_doses, phi_I_lower, phi_E_lower, toxicity_low) {
  # Create flat probability matrices for all doses
  #
  # Args:
  #   n_doses: Number of dose levels
  #   phi_I_lower: Immune response rate for all doses
  #   phi_E_lower: Marginal efficacy rate for all doses
  #   toxicity_low: Low toxicity rate for all doses
  #
  # Returns:
  #   list: Contains p_YI, p_YT_given_I, p_YE_given_I matrices
  
  # Immune response: flat across all doses
  p_YI_flat <- rep(phi_I_lower, n_doses)
  
  # Toxicity: low rate for all doses, both I=0 and I=1
  p_YT_given_I_flat <- matrix(rep(toxicity_low, 2 * n_doses), ncol = 2, byrow = TRUE)
  
  # Efficacy: use total probability formula to maintain marginal = phi_E_lower
  p_YE_given_I_flat <- calculate_conditional_efficacy_flat(phi_E_lower, phi_I_lower)
  # Replicate for all doses
  p_YE_given_I_flat <- matrix(rep(p_YE_given_I_flat, n_doses), ncol = 2, byrow = TRUE)
  
  return(list(
    p_YI = p_YI_flat,
    p_YT_given_I = p_YT_given_I_flat,
    p_YE_given_I = p_YE_given_I_flat
  ))
}

generate_flat_scenario_data <- function(config, phi_I_lower, phi_E_lower, toxicity_low = 0.05, 
                                       n_patients_per_dose = 10, seed = 123, debug = FALSE) {
  # Generate flat scenario data where all doses have identical probabilities at lower bounds
  #
  # Args:
  #   config: Trial configuration containing dose_levels, rho0, rho1
  #   phi_I_lower: Immune response rate for all doses (e.g., 0.20)
  #   phi_E_lower: Marginal efficacy rate for all doses (e.g., 0.25)
  #   toxicity_low: Low toxicity rate for all doses (e.g., 0.05)
  #   n_patients_per_dose: Number of patients per dose level
  #   seed: Random seed for reproducibility
  #   debug: Whether to store debug information
  #
  # Returns:
  #   data.frame: Simulated trial data with flat probabilities
  
  n_doses <- length(config$dose_levels)
  
  # Create flat probability matrices
  flat_probs <- create_flat_probability_matrices(n_doses, phi_I_lower, phi_E_lower, toxicity_low)
  
  # Generate data using existing simulation function
  data <- simulate_data_gumbel(
    n_per_dose_vector = rep(n_patients_per_dose, n_doses),
    dose_levels = config$dose_levels,
    p_YI = flat_probs$p_YI,
    p_YT_given_I = flat_probs$p_YT_given_I,
    p_YE_given_I = flat_probs$p_YE_given_I,
    rho0 = config$rho0,
    rho1 = config$rho1,
    seed = seed,
    debug = debug
  )
  
  # Add metadata about the flat scenario
  attr(data, "scenario_type") <- "flat_null"
  attr(data, "phi_I_lower") <- phi_I_lower
  attr(data, "phi_E_lower") <- phi_E_lower
  attr(data, "toxicity_low") <- toxicity_low
  
  return(data)
}

validate_flat_scenario <- function(data, phi_I_lower, phi_E_lower, toxicity_low, tolerance = 0.1) {
  # Validate that generated data represents a flat scenario
  #
  # Args:
  #   data: Generated trial data
  #   phi_I_lower: Expected immune response rate
  #   phi_E_lower: Expected marginal efficacy rate
  #   toxicity_low: Expected toxicity rate
  #   tolerance: Tolerance for validation (default 0.1 for sampling variability)
  #
  # Returns:
  #   list: Validation results with success flag and details
  
  validation_results <- list(
    success = TRUE,
    immune_response_rates = numeric(),
    efficacy_rates = numeric(),
    toxicity_rates = numeric(),
    marginal_efficacy_rates = numeric(),
    details = character()
  )
  
  # Calculate observed rates by dose
  dose_levels <- unique(data$d)
  
  for (dose in dose_levels) {
    dose_data <- data[data$d == dose, ]
    
    # Immune response rate
    immune_rate <- mean(dose_data$Y_I)
    validation_results$immune_response_rates <- c(validation_results$immune_response_rates, immune_rate)
    
    # Efficacy rate
    efficacy_rate <- mean(dose_data$Y_E)
    validation_results$efficacy_rates <- c(validation_results$efficacy_rates, efficacy_rate)
    
    # Toxicity rate
    toxicity_rate <- mean(dose_data$Y_T)
    validation_results$toxicity_rates <- c(validation_results$toxicity_rates, toxicity_rate)
    
    # Check if rates are within tolerance of expected values
    if (abs(immune_rate - phi_I_lower) > tolerance) {
      validation_results$success <- FALSE
      validation_results$details <- c(validation_results$details, 
                                   paste("Dose", dose, "immune rate", round(immune_rate, 3), 
                                         "differs from expected", phi_I_lower, "by", 
                                         round(abs(immune_rate - phi_I_lower), 3)))
    }
    
    if (abs(efficacy_rate - phi_E_lower) > tolerance) {
      validation_results$success <- FALSE
      validation_results$details <- c(validation_results$details, 
                                   paste("Dose", dose, "efficacy rate", round(efficacy_rate, 3), 
                                         "differs from expected", phi_E_lower, "by", 
                                         round(abs(efficacy_rate - phi_E_lower), 3)))
    }
    
    if (abs(toxicity_rate - toxicity_low) > tolerance) {
      validation_results$success <- FALSE
      validation_results$details <- c(validation_results$details, 
                                   paste("Dose", dose, "toxicity rate", round(toxicity_rate, 3), 
                                         "differs from expected", toxicity_low, "by", 
                                         round(abs(toxicity_rate - toxicity_low), 3)))
    }
  }
  
  # Check if rates are approximately flat across doses (within tolerance)
  # Use a more lenient check for flatness that accounts for sampling variability
  immune_range <- max(validation_results$immune_response_rates) - min(validation_results$immune_response_rates)
  efficacy_range <- max(validation_results$efficacy_rates) - min(validation_results$efficacy_rates)
  toxicity_range <- max(validation_results$toxicity_rates) - min(validation_results$toxicity_rates)
  
  if (immune_range > tolerance) {
    validation_results$success <- FALSE
    validation_results$details <- c(validation_results$details, 
                                   paste("Immune response rates vary by", round(immune_range, 3), 
                                         "across doses (tolerance:", tolerance, ")"))
  }
  
  if (efficacy_range > tolerance) {
    validation_results$success <- FALSE
    validation_results$details <- c(validation_results$details, 
                                   paste("Efficacy rates vary by", round(efficacy_range, 3), 
                                         "across doses (tolerance:", tolerance, ")"))
  }
  
  if (toxicity_range > tolerance) {
    validation_results$success <- FALSE
    validation_results$details <- c(validation_results$details, 
                                   paste("Toxicity rates vary by", round(toxicity_range, 3), 
                                         "across doses (tolerance:", tolerance, ")"))
  }
  
  return(validation_results)
}