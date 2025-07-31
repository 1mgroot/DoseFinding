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

apply_biviso_on_matrix <- function(posterior_df, dose_col = "d", group_col = "Y_I", n_sims = 1000) {
  # Ensure all dose-group combinations exist
  dose_levels <- sort(unique(posterior_df[[dose_col]]))
  group_levels <- sort(unique(posterior_df[[group_col]]))
  
  # Create complete grid
  complete_grid <- expand.grid(d = dose_levels, Y_I = group_levels)
  
  # Merge with existing data, filling missing combinations with zeros
  df <- complete_grid %>%
    left_join(posterior_df, by = c("d", "Y_I")) %>%
    replace_na(list(r = 0, n = 0, alpha_post = 1, beta_post = 1, 
                   var_post = 0.25, mean_post = 0.5)) %>%
    arrange(.data[[dose_col]], .data[[group_col]])
  
  # Regenerate samples for missing combinations
  missing_rows <- is.na(df$samples)
  if (any(missing_rows)) {
    df$samples[missing_rows] <- lapply(which(missing_rows), function(i) {
      rbeta(n_sims, df$alpha_post[i], df$beta_post[i])
    })
  }

  J <- length(dose_levels)
  G <- length(group_levels)

  sample_matrix <- do.call(cbind, df$samples)
  wts <- matrix(1 / df$var_post, nrow = J, ncol = G)

  adjusted_matrix <- t(apply(sample_matrix, 1, function(row) {
    mat <- matrix(row, nrow = J, ncol = G)
    as.vector(Iso::biviso(y = mat, w = wts))
  }))

  samples_pava_list <- lapply(1:ncol(adjusted_matrix), function(j) adjusted_matrix[, j])
  df$samples_pava <- samples_pava_list
  df$pava_mean <- colMeans(adjusted_matrix)
  df$pava_ci_lower <- apply(adjusted_matrix, 2, quantile, probs = 0.025)
  df$pava_ci_upper <- apply(adjusted_matrix, 2, quantile, probs = 0.975)

  return(df)
}

compute_marginal_probability <- function(group_post_df, immune_post_df, dose_col = "d", group_col = "Y_I") {
  # Extract immune response probability per dose
  immune_samples <- do.call(cbind, immune_post_df$samples_pava)
  
  # Get samples for each group - need to ensure proper ordering
  group_df_ordered <- group_post_df %>% arrange(.data[[dose_col]], .data[[group_col]])
  group_samples_0 <- do.call(cbind, group_df_ordered[group_df_ordered[group_col] == 0, ]$samples_pava)
  group_samples_1 <- do.call(cbind, group_df_ordered[group_df_ordered[group_col] == 1, ]$samples_pava)

  # Validate dimensions
  if (ncol(immune_samples) != ncol(group_samples_0) || ncol(immune_samples) != ncol(group_samples_1)) {
    stop("Dimension mismatch in marginal probability calculation. Immune samples: ", 
         ncol(immune_samples), ", Group 0 samples: ", ncol(group_samples_0), 
         ", Group 1 samples: ", ncol(group_samples_1))
  }

  # Compute marginal probability samples
  marginal_samples <- (1 - immune_samples) * group_samples_0 + immune_samples * group_samples_1

  # Create a summary data frame
  marginal_summary <- data.frame(
    d = immune_post_df[[dose_col]],
    marginal_prob = colMeans(marginal_samples),
    samples = I(lapply(1:ncol(marginal_samples), function(i) marginal_samples[, i]))
  )
  
  return(marginal_summary)
}