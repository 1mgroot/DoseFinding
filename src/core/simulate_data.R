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
  set.seed(seed)
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