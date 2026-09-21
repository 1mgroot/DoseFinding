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
  adjusted_matrix <- if (ncol(sample_matrix) == 1L) sample_matrix else t(apply(sample_matrix, 1, function(row) {
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
  group_levels <- 0:1
  if (anyDuplicated(posterior_df[c(dose_col, group_col)])) {
    stop("Duplicate dose/immune strata in posterior summaries.")
  }
  if (!all(posterior_df[[group_col]] %in% group_levels)) {
    stop("Immune strata must be coded 0 or 1.")
  }
  
  # Create complete grid
  complete_grid <- do.call(expand.grid, setNames(list(dose_levels, group_levels),
                                                 c(dose_col, group_col)))
  
  # Merge with existing data, filling missing combinations with zeros
  df <- complete_grid %>%
    left_join(posterior_df, by = c(dose_col, group_col)) %>%
    replace_na(list(r = 0, n = 0, alpha_post = 1, beta_post = 1)) %>%
    arrange(.data[[dose_col]], .data[[group_col]])
  moments <- add_beta_variance(df)
  df$var_post <- dplyr::coalesce(df$var_post, moments$var_post)
  df$mean_post <- dplyr::coalesce(df$mean_post, moments$mean_post)
  
  # Regenerate samples for missing combinations
  missing_rows <- lengths(df$samples) == 0 | vapply(df$samples, function(x) all(is.na(x)), FALSE)
  existing_lengths <- unique(lengths(df$samples[!missing_rows]))
  if (length(existing_lengths) > 1) stop("Posterior strata must have equal draw counts.")
  if (length(existing_lengths)) n_sims <- existing_lengths[[1]]
  if (any(missing_rows)) {
    df$samples[missing_rows] <- lapply(which(missing_rows), function(i) {
      rbeta(n_sims, df$alpha_post[i], df$beta_post[i])
    })
  }
  if (any(vapply(df$samples, function(x) any(!is.finite(x)), FALSE))) {
    stop("Posterior samples must be finite.")
  }

  J <- length(dose_levels)
  G <- length(group_levels)

  sample_matrix <- do.call(cbind, df$samples)
  
  # Create the weights matrix, ensuring it's filled by row to match the sample data
  wts_mat <- matrix(1 / df$var_post, nrow = J, ncol = G, byrow = TRUE)

  adjusted_matrix <- t(apply(sample_matrix, 1, function(row) {
    # Reshape the row into a matrix, filled by row
    mat <- matrix(row, nrow = J, ncol = G, byrow = TRUE)
    
    # Apply biviso
    result_mat <- if (J == 1L) {
      matrix(Iso::pava(as.numeric(mat), as.numeric(wts_mat)), nrow = 1L)
    } else {
      Iso::biviso(y = mat, w = wts_mat)
    }
    
    # Convert the result back to a vector, ensuring the order is preserved
    as.vector(t(result_mat))
  }))

  samples_pava_list <- lapply(1:ncol(adjusted_matrix), function(j) adjusted_matrix[, j])
  df$samples_pava <- samples_pava_list
  df$pava_mean <- colMeans(adjusted_matrix)
  df$pava_ci_lower <- apply(adjusted_matrix, 2, quantile, probs = 0.025)
  df$pava_ci_upper <- apply(adjusted_matrix, 2, quantile, probs = 0.975)

  return(df)
}

compute_marginal_probability <- function(group_post_df, immune_post_df, dose_col = "d", group_col = "Y_I") {
  if (anyDuplicated(immune_post_df[[dose_col]]) ||
      anyDuplicated(group_post_df[c(dose_col, group_col)])) {
    stop("Duplicate doses or dose/immune strata in marginalization inputs.")
  }
  if (!all(group_post_df[[group_col]] %in% 0:1)) stop("Immune strata must be coded 0 or 1.")
  if (length(unique(lengths(immune_post_df$samples_pava))) != 1L ||
      any(lengths(immune_post_df$samples_pava) == 0L)) {
    stop("Immune posterior draw counts must match and be nonzero for marginalization.")
  }
  # Extract immune response probability per dose
  immune_samples <- do.call(cbind, immune_post_df$samples_pava)
  
  # Get samples for each group - need to ensure proper ordering
  group_df_ordered <- group_post_df %>% arrange(.data[[dose_col]], .data[[group_col]])
  strata <- lapply(0:1, function(g) {
    x <- group_df_ordered[group_df_ordered[[group_col]] == g, , drop = FALSE]
    if (!setequal(x[[dose_col]], immune_post_df[[dose_col]])) {
      stop("Dose labels must match in both immune strata and the immune posterior.")
    }
    x <- x[match(immune_post_df[[dose_col]], x[[dose_col]]), , drop = FALSE]
    if (any(lengths(x$samples_pava) != nrow(immune_samples))) {
      stop("Posterior draw counts must match for marginalization.")
    }
    do.call(cbind, x$samples_pava)
  })
  group_samples_0 <- strata[[1]]
  group_samples_1 <- strata[[2]]

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
