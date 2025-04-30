library(dplyr)

simulate_data <- function(
    n_per_dose = 10,
    dose_levels = c(1, 2, 3),
    p_YI = c(0.2, 0.5, 0.8),   # Real Immune response rate per dose level
    p_YT_given_I = c(0.2, 0.6), # Real Toxicity response rate，when I=0, I=1
    p_YE_given_I = c(0.3, 0.7), # Real Efficacy response rate，when I=0, I=1
    seed = 118,
    debug=FALSE
) {
  set.seed(seed)
  J <- length(dose_levels)
  n_total <- n_per_dose * J
  d <- rep(dose_levels, each = n_per_dose)

  # 根据 p_YI[dose] 生成 Y_I
  p_YI_full <- setNames(p_YI, dose_levels)
  Y_I <- sapply(d, function(dj) rbinom(1, 1, p_YI_full[as.character(dj)]))

  # 根据 Y_I 决定 Y_T, Y_E
  Y_T <- sapply(Y_I, function(I) rbinom(1, 1, p_YT_given_I[I + 1]))
  Y_E <- sapply(Y_I, function(I) rbinom(1, 1, p_YE_given_I[I + 1]))

  df = data.frame(id = 1:n_total, d = d, Y_I = Y_I, Y_T = Y_T, Y_E = Y_E)

  if (debug) {
    assign("Y_I", Y_I, envir = .GlobalEnv)
    assign("Y_T", Y_T, envir = .GlobalEnv)
    assign("J", J, envir = .GlobalEnv)
  }
  return(df)
}

compute_rn <- function(df, outcome_col, group_col = "Y_I", dose_col = "d") {
  # df: 数据框，必须包含 dose_col, group_col, outcome_col 三列
  # outcome_col: "Y_T" 或 "Y_E"
  df %>%
    group_by(!!sym(dose_col), !!sym(group_col)) %>%
    summarise(
      r = sum(!!sym(outcome_col)),
      n = n(),
      .groups = "drop"
    ) %>%
    rename(
      dose = !!sym(dose_col),
      group = !!sym(group_col)
    )
}

compute_beta_posterior <- function(stats_df, alpha = 1, beta = 1) {
  stats_df %>%
    mutate(
      alpha_post = r + alpha,
      beta_post = n - r + beta
    )
}

compute_marginal_probability <- function(post_df, pi_I_estimate) {
  post_df %>%
    mutate(
      post_mean = alpha_post / (alpha_post + beta_post)
    ) %>%
    select(dose, group, post_mean) %>%
    tidyr::pivot_wider(names_from = group, values_from = post_mean, names_prefix = "pi_T_") %>%
    mutate(
      pi_I = pi_I_estimate[as.character(dose)],
      pi_T_marg = pi_I * pi_T_1 + (1 - pi_I) * pi_T_0
    )
}



df <- simulate_data(
  n_per_dose = 10,
  dose_levels = c(1, 2, 3),
  p_YI = c(0.2, 0.3, 0.4),
  p_YT_given_I = c(0.1, 0.5),
  p_YE_given_I = c(0.3, 0.7),
  debug=TRUE
)
tox_stats <- compute_rn(df, outcome_col = "Y_T")
eff_stats <- compute_rn(df, outcome_col = "Y_E")

tox_post <- compute_beta_posterior(tox_stats, alpha = 1, beta = 1)



print(df)
print(tox_stats)
print(eff_stats)
print(tox_post)

