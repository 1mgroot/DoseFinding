compute_rn <- function(df, outcome_col, group_col = NULL, dose_col = "d") {
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

simulate_beta_posterior <- function(stats_df, alpha = 1, beta = 1, n_sims = 1000) {
  stats_df |>
    dplyr::mutate(
      alpha_post = r + alpha,
      beta_post = n - r + beta,
      samples = purrr::map2(alpha_post, beta_post, ~ rbeta(n_sims, .x, .y))
    )
}

add_beta_variance <- function(posterior_df) {
  posterior_df |>
    dplyr::mutate(
      var_post = (alpha_post * beta_post) / ((alpha_post + beta_post)^2 * (alpha_post + beta_post + 1)),
      mean_post = alpha_post / (alpha_post + beta_post)
    )
}

apply_pava_on_samples <- function(posterior_df, dose_col = "d", n_sims = 1000) {
  posterior_df <- posterior_df |> dplyr::arrange(.data[[dose_col]])
  sample_matrix <- do.call(cbind, posterior_df$samples)
  weights <- 1 / posterior_df$var_post
  dose_order <- posterior_df[[dose_col]]
  adjusted_matrix <- t(apply(sample_matrix, 1, function(row) {
    #gpava(z = dose_order, y = row, weights = weights, solver = weighted.mean)$x
    Iso::pava(y = row, w = weights)
  }))
  samples_pava_list <- lapply(1:ncol(adjusted_matrix), function(j) adjusted_matrix[, j])
  posterior_df$samples_pava <- samples_pava_list
  posterior_df$pava_mean <- colMeans(adjusted_matrix)
  posterior_df$pava_ci_lower <- apply(adjusted_matrix, 2, quantile, probs = 0.025)
  posterior_df$pava_ci_upper <- apply(adjusted_matrix, 2, quantile, probs = 0.975)
  return(posterior_df)
}

apply_biviso_on_matrix <- function(posterior_df, dose_col = "d", group_col = "Y_I") {
  df <- posterior_df %>% arrange(.data[[dose_col]], .data[[group_col]])

  dose_levels <- sort(unique(df[[dose_col]]))
  group_levels <- sort(unique(df[[group_col]]))
  J <- length(dose_levels)
  G <- length(group_levels)

  mat <- matrix(df$mean_post, nrow = J, ncol = G)
  wts <- matrix(1 / df$var_post, nrow = J, ncol = G)

  result <- Iso::biviso(y = mat, w = wts)
  df$pava_mean <- as.vector(result)

  return(df)
}

compute_marginal_probability <- function(group_post_df, immune_post_df, dose_col = "d", group_col = "Y_I") {
  # Extract immune response probability per dose
  immune_probs <- immune_post_df %>% 
    dplyr::select(!!dose_col, pi_I = mean_post)
  
  # Merge immune probabilities into toxicity/efficacy posterior estimates
  df <- group_post_df %>% 
    dplyr::left_join(immune_probs, by = dose_col)
  
  # Compute marginal probability as weighted average:
  # If group_col == 1: weight = pi_I
  # If group_col == 0: weight = (1 - pi_I)
  df <- df %>% 
    dplyr::mutate(
      weighted_prob = dplyr::case_when(
        .data[[group_col]] == 1 ~ pi_I * mean_post,
        .data[[group_col]] == 0 ~ (1 - pi_I) * mean_post
      )
    )
  
  # Summarize weighted probabilities by dose
  marginal_summary <- df %>% 
    dplyr::group_by(!!sym(dose_col)) %>% 
    dplyr::summarise(marginal_prob = sum(weighted_prob), .groups = "drop")
  
  return(marginal_summary)
}