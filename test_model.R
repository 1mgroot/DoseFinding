library(dplyr)
library(testthat)
library(tidyr)
library(isotone)
library(purrr)
library(ggplot2)
library(Iso)

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


# Gumbel
Gumbel <- function(ttox, teff, c) {
  out <- matrix(0, nrow = 4, ncol = 1)
  out[1, ] <- (1 - ttox) * (1 - teff) + ttox * teff * (1 - ttox) * (1 - teff) * (exp(c) - 1) / (exp(c) + 1)
  out[2, ] <- (1 - ttox) * teff - ttox * teff * (1 - ttox) * (1 - teff) * (exp(c) - 1) / (exp(c) + 1)
  out[3, ] <- ttox * (1 - teff) - ttox * teff * (1 - ttox) * (1 - teff) * (exp(c) - 1) / (exp(c) + 1)
  out[4, ] <- ttox * teff + ttox * teff * (1 - ttox) * (1 - teff) * (exp(c) - 1) / (exp(c) + 1)
  return(out)
}

# simulate_data_gumbel
simulate_data_gumbel <- function(
    n_per_dose = 10,
    dose_levels = c(1, 2, 3),
    p_YI = c(0.2, 0.5, 0.8), # immune prob per dose
    p_YT_given_I, # marginal tox prob for I=0, I=1
    p_YE_given_I, # marginal eff prob for I=0, I=1
    rho0 = 1, # correlation under I=0
    rho1 = 1, # correlation under I=1
    seed = 123,
    debug = FALSE) {
  set.seed(seed)
  J <- length(dose_levels)
  d_vec <- rep(dose_levels, each = n_per_dose)
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


# Visualization
plot_posterior_summary <- function(posterior_df,
                                   dose_col = "d",
                                   group_col = NULL,
                                   mean_col = "pava_mean",
                                   ci_lower_col = "pava_ci_lower",
                                   ci_upper_col = "pava_ci_upper",
                                   title = "Posterior Mean and 95% CI by Dose") {
  p <- ggplot(posterior_df, aes(x = .data[[dose_col]], y = .data[[mean_col]]))

  if (!is.null(group_col) && group_col %in% colnames(posterior_df)) {
    p <- p +
      aes(color = factor(.data[[group_col]]), group = .data[[group_col]]) +
      geom_point(size = 3) +
      geom_line()

    if (ci_lower_col %in% colnames(posterior_df) && ci_upper_col %in% colnames(posterior_df)) {
      p <- p +
        geom_errorbar(aes(ymin = .data[[ci_lower_col]], ymax = .data[[ci_upper_col]]), width = 0.2)
    }

    p <- p + scale_color_brewer(palette = "Set1", name = group_col)
  } else {
    p <- p +
      geom_point(size = 3, color = "blue") +
      geom_line(color = "blue")

    if (ci_lower_col %in% colnames(posterior_df) && ci_upper_col %in% colnames(posterior_df)) {
      p <- p +
        geom_errorbar(aes(ymin = .data[[ci_lower_col]], ymax = .data[[ci_upper_col]]),
                      width = 0.2, color = "blue")
    }
  }

  p +
    scale_x_continuous(breaks = unique(posterior_df[[dose_col]])) +
    labs(x = "Dose Level", y = "Posterior Mean (with 95% CI)", title = title) +
    theme_minimal(base_size = 14)
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

cat("🔍 Starting test suite...\n")

## ===== 1. Test compute_rn: basic aggregation of response and sample size =====
df_test <- data.frame(d = c(1, 1, 2), Y_T = c(1, 0, 1))
res_rn <- compute_rn(df_test, outcome_col = "Y_T")
expect_equal(res_rn$r, c(1, 1))
expect_equal(res_rn$n, c(2, 1))

cat("✅ compute_rn passed\n")

## ===== 2. Test simulate_beta_posterior + add_beta_variance: posterior mean correctness =====
beta_df <- simulate_beta_posterior(res_rn, alpha = 1, beta = 1, n_sims = 1000)
beta_df <- add_beta_variance(beta_df)
manual_mean <- with(res_rn, (r + 1) / (n + 2))
expect_true(all(abs(beta_df$mean_post - manual_mean) < 1e-6))

cat("✅ Beta posterior mean passed\n")

## ===== 3. Test apply_pava_on_samples: isotonic adjustment results in non-decreasing posterior means =====
posterior_df <- tibble::tibble(
  d = 1:3,
  alpha_post = c(2, 2, 2),
  beta_post = c(2, 2, 2),
  var_post = c(0.05, 0.05, 0.05),
  samples = list(
    rep(0.1, 1000),
    rep(0.6, 1000),
    rep(0.3, 1000)
  )
)
pava_res <- apply_pava_on_samples(posterior_df)
expect_true(all(diff(pava_res$pava_mean) >= -1e-8))

cat("✅ apply_pava_on_samples passed\n")

## ===== 4. Test simulate_data_gumbel: marginal and conditional probabilities are approximately correct =====
p_YT_given_I <- matrix(c(0.1, 0.3, 0.3, 0.5, 0.5, 0.7), ncol = 2, byrow = TRUE)
p_YE_given_I <- matrix(c(0.2, 0.4, 0.4, 0.6, 0.6, 0.8), ncol = 2, byrow = TRUE)

df <- simulate_data_gumbel(
  n_per_dose = 1000,
  dose_levels = c(1, 2, 3),
  p_YI = c(0.2, 0.4, 0.6),
  p_YT_given_I = p_YT_given_I,
  p_YE_given_I = p_YE_given_I,
  rho0 = 1.5,
  rho1 = 2,
  debug = FALSE
)

# Check immune rate matches specification
immune_rate <- df %>% group_by(d) %>% summarise(immune = mean(Y_I))
expect_true(all(abs(immune_rate$immune - c(0.2, 0.4, 0.6)) < 0.05))

# Check conditional toxicity matches input matrix
tox_rate <- df %>% group_by(d, Y_I) %>% summarise(tox = mean(Y_T), .groups = "drop")
diffs <- abs(tox_rate$tox - c(0.1, 0.3, 0.3, 0.5, 0.5, 0.7))
expect_true(all(diffs < 0.05))

cat("✅ simulate_data_gumbel passed\n")

## ===== 5. Test apply_biviso_on_matrix: matrix after adjustment should be isotonic (row- and column-wise) =====
biviso_df <- tibble::tibble(
  d = rep(1:3, each = 2),
  Y_I = rep(0:1, times = 3),
  mean_post = c(0.6, 0.2, 0.5, 0.3, 0.4, 0.1),
  var_post = rep(0.01, 6)
)
biviso_df$alpha_post <- 2
biviso_df$beta_post <- 2
biviso_df$samples <- list(rep(0.5, 1000))  # dummy

res <- apply_biviso_on_matrix(biviso_df)

# Reshape and verify monotonicity
mat <- matrix(res$pava_mean, nrow = 3, ncol = 2)
row_monotone <- all(apply(mat, 1, function(x) all(diff(x) >= -1e-6)))
col_monotone <- all(apply(mat, 2, function(x) all(diff(x) >= -1e-6)))
expect_true(row_monotone && col_monotone)

cat("✅ apply_biviso_on_matrix passed\n")

## ===== 6. Test compute_marginal_probability: weighted sum matches hand calculation =====
group_post_df <- tibble::tibble(
  d = c(1, 1),
  Y_I = c(0, 1),
  mean_post = c(0.1, 0.9)
)
immune_post_df <- tibble::tibble(
  d = 1,
  mean_post = 0.3  # P(Y_I = 1 | d = 1)
)
marginal_df <- compute_marginal_probability(group_post_df, immune_post_df)

# Manual check: 0.1 * 0.7 + 0.9 * 0.3 = 0.07 + 0.27 = 0.34
expect_true(abs(marginal_df$marginal_prob - 0.34) < 1e-6)

cat("✅ compute_marginal_probability passed\n")

cat("🎉 All tests passed successfully ✅\n")
