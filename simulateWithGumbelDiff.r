library(dplyr)
library(tidyr)
library(isotone)
library(purrr)
library(ggplot2)

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
        gpava(z = dose_order, y = row, weights = weights, solver = weighted.mean)$x
    }))
    samples_pava_list <- lapply(1:ncol(adjusted_matrix), function(j) adjusted_matrix[, j])
    posterior_df$samples_pava <- samples_pava_list
    posterior_df$pava_mean <- colMeans(adjusted_matrix)
    posterior_df$pava_ci_lower <- apply(adjusted_matrix, 2, quantile, probs = 0.025)
    posterior_df$pava_ci_upper <- apply(adjusted_matrix, 2, quantile, probs = 0.975)
    return(posterior_df)
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

# Gumbel +simulate_beta_posterior + PAVA
p_YT_given_I <- matrix(c(
    0.1, 0.3,
    0.3, 0.5,
    0.5, 0.7
), ncol = 2, byrow = TRUE)

p_YE_given_I <- matrix(c(
    0.2, 0.4,
    0.4, 0.6,
    0.6, 0.8
), ncol = 2, byrow = TRUE)

df <- simulate_data_gumbel(
    n_per_dose = 1000,
    dose_levels = c(1, 2, 3),
    p_YI = c(0.2, 0.4, 0.6),
    p_YT_given_I = p_YT_given_I,
    p_YE_given_I = p_YE_given_I,
    rho0 = 1.5,
    rho1 = 2,
    debug = TRUE
)

pi_I_stats <- compute_rn(df, outcome_col = "Y_I")
pi_I_post <- simulate_beta_posterior(pi_I_stats)
pi_I_post <- add_beta_variance(pi_I_post)
pi_I_pava <- apply_pava_on_samples(pi_I_post)

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
            geom_errorbar(aes(ymin = .data[[ci_lower_col]], ymax = .data[[ci_upper_col]]), width = 0.2) +
            geom_line() +
            scale_color_brewer(palette = "Set1", name = group_col)
    } else {
        p <- p +
            geom_point(size = 3, color = "blue") +
            geom_errorbar(aes(ymin = .data[[ci_lower_col]], ymax = .data[[ci_upper_col]]), width = 0.2, color = "blue") +
            geom_line(color = "blue")
    }
    p +
        scale_x_continuous(breaks = unique(posterior_df[[dose_col]])) +
        labs(x = "Dose Level", y = "Posterior Mean (with 95% CI)", title = title) +
        theme_minimal(base_size = 14)
}

plot_posterior_summary(pi_I_pava, title = "Immune Response vs Dose (PAVA Adjusted)")

tox_stats <- compute_rn(df, outcome_col = "Y_T", group_col = "Y_I")
tox_post <- simulate_beta_posterior(tox_stats)
tox_post <- add_beta_variance(tox_post)
tox_pava <- apply_pava_on_samples(tox_post)
plot_posterior_summary(tox_pava,
    title = "Toxicity Rate by Dose and Immune Status",
    group_col = "Y_I"
)

eff_stats <- compute_rn(df, outcome_col = "Y_E", group_col = "Y_I")
eff_post <- simulate_beta_posterior(eff_stats)
eff_post <- add_beta_variance(eff_post)
eff_pava <- apply_pava_on_samples(eff_post)
plot_posterior_summary(eff_pava,
    title = "Efficacy Rate by Dose and Immune Status",
    group_col = "Y_I"
)
