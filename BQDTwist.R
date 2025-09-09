library(dplyr)
library(tidyr)
library(isotone)
library(purrr)

#' Simulate dataset for dose finding study
#' 
#' @param n_per_dose Number of subjects per dose level
#' @param dose_levels Vector of dose levels
#' @param p_YI Real immune response rate per dose level
#' @param p_YT_given_I Real toxicity response rate when I=0 or I=1
#' @param p_YE_given_I Real efficacy response rate when I=0 or I=1
#' @param seed Random seed for reproducibility
#' @param debug Whether to assign variables to global environment
#' @return Dataframe with simulated responses
simulate_data <- function(
  n_per_dose = 10,
  dose_levels = c(1, 2, 3),
  p_YI = c(0.2, 0.5, 0.8),        # Real immune response rate per dose level
  p_YT_given_I = c(0.2, 0.6),     # Real toxicity response rate when I=0, I=1
  p_YE_given_I = c(0.3, 0.7),     # Real efficacy response rate when I=0, I=1
  seed = 118,
  debug = FALSE
) {
  set.seed(seed)
  J <- length(dose_levels)
  n_total <- n_per_dose * J
  d <- rep(dose_levels, each = n_per_dose)
  
  # Generate Y_I based on p_YI[dose]
  p_YI_full <- setNames(p_YI, dose_levels)
  Y_I <- sapply(d, function(dj) rbinom(1, 1, p_YI_full[as.character(dj)]))
  
  # Generate Y_T, Y_E based on Y_I
  Y_T <- sapply(Y_I, function(I) rbinom(1, 1, p_YT_given_I[I + 1]))
  Y_E <- sapply(Y_I, function(I) rbinom(1, 1, p_YE_given_I[I + 1]))
  
  df <- data.frame(id = 1:n_total, d = d, Y_I = Y_I, Y_T = Y_T, Y_E = Y_E)
  
  if (debug) {
    assign("Y_I", Y_I, envir = .GlobalEnv)
    assign("Y_T", Y_T, envir = .GlobalEnv)
    assign("J", J, envir = .GlobalEnv)
  }
  return(df)
}

#' Compute response rates and counts
#' 
#' @param df Data frame with outcome data
#' @param outcome_col Column name for outcome variable
#' @param group_col Optional grouping column
#' @param dose_col Column name for dose variable
#' @return Data frame with response rates and counts by group
compute_rn <- function(df, outcome_col, group_col = NULL, dose_col = "d") {
  # Construct grouping columns
  group_vars <- if (!is.null(group_col) && nzchar(group_col)) {
    c(dose_col, group_col)
  } else {
    c(dose_col)
  }
  
  df %>%
    dplyr::group_by_at(group_vars) %>%
    dplyr::summarise(
      r = sum(.data[[outcome_col]], na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
}

#' Simulate beta posterior distribution
#' 
#' @param stats_df Data frame with response rates
#' @param alpha Prior alpha parameter
#' @param beta Prior beta parameter
#' @param n_sims Number of posterior samples
#' @return Data frame with posterior samples
simulate_beta_posterior <- function(stats_df, alpha = 1, beta = 1, n_sims = 1000) {
  stats_df |>
    dplyr::mutate(
      alpha_post = r + alpha,
      beta_post = n - r + beta,
      samples = purrr::map2(alpha_post, beta_post, ~ rbeta(n_sims, .x, .y))
    )
}

#' Add beta posterior variance and mean
#' 
#' @param posterior_df Data frame with posterior parameters
#' @return Data frame with added variance and mean
add_beta_variance <- function(posterior_df) {
  posterior_df |>
    dplyr::mutate(
      var_post = (alpha_post * beta_post) /
        ((alpha_post + beta_post)^2 * (alpha_post + beta_post + 1)),
      mean_post = alpha_post / (alpha_post + beta_post)
    )
}

#' Apply PAVA (Pool-Adjacent-Violators Algorithm) on posterior samples
#' 
#' @param posterior_df Data frame with posterior samples
#' @param dose_col Column name for dose variable
#' @param n_sims Number of posterior samples
#' @return Data frame with PAVA-adjusted samples and summaries
apply_pava_on_samples <- function(posterior_df, dose_col = "d", n_sims = 1000) {
  posterior_df <- posterior_df |> dplyr::arrange(.data[[dose_col]])
  
  # Construct posterior sample matrix (rows: samples, cols: dose levels)
  sample_matrix <- do.call(cbind, posterior_df$samples)  # n_sims x n_dose
  
  # Weight vector (one per dose)
  weights <- 1 / posterior_df$var_post
  dose_order <- posterior_df[[dose_col]]
  
  # Apply weighted PAVA to each sample path
  adjusted_matrix <- t(apply(sample_matrix, 1, function(row) {
    gpava(z = dose_order, y = row, weights = weights, solver = weighted.mean)$x
  }))
  
  samples_pava_list <- lapply(1:ncol(adjusted_matrix), function(j) adjusted_matrix[, j])
  posterior_df$samples_pava <- samples_pava_list
  
  # Add summary statistics (mean, confidence intervals)
  posterior_df$pava_mean <- colMeans(adjusted_matrix)
  posterior_df$pava_ci_lower <- apply(adjusted_matrix, 2, quantile, probs = 0.025)
  posterior_df$pava_ci_upper <- apply(adjusted_matrix, 2, quantile, probs = 0.975)
  
  return(posterior_df)
}

# Example usage
df <- simulate_data(
  n_per_dose = 10,
  dose_levels = c(1, 2, 3),
  p_YI = c(0.2, 0.3, 0.4),
  p_YT_given_I = c(0.1, 0.5),
  p_YE_given_I = c(0.3, 0.7),
  debug = TRUE
)

# Uncommented analysis steps
# imm_stats <- compute_rn(df, outcome_col = "Y_I")
# tox_stats <- compute_rn(df, group_col = "Y_I", outcome_col = "Y_T")
# eff_stats <- compute_rn(df, group_col = "Y_I", outcome_col = "Y_E")

# Analyze immune response
pi_I_stats <- compute_rn(df, outcome_col = "Y_I")
pi_I_post <- simulate_beta_posterior(pi_I_stats, alpha = 1, beta = 1, n_sims = 1000)
pi_I_post <- add_beta_variance(pi_I_post)
pi_I_pava <- apply_pava_on_samples(pi_I_post)

# Simulate enrollment times
# n, rate (number of enrollments per month)
x <- rexp(30, 3)
x <- cumsum(x)
# PFS -> time