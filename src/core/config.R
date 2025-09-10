library(dplyr)
library(tidyr)
library(isotone)
library(purrr)
library(ggplot2)
library(Iso)

trial_config <- list(
  dose_levels = c(1, 2, 3),
  n_stages = 3,
  cohort_size = 6,
  phi_T = 0.3, c_T = 0.9,
  phi_E = 0.2, c_E = 0.9,
  phi_I = 0.2, c_I = 0.8,
  # PoC parameters
  c_poc = 0.9,
  delta_poc = 0.8,  # Threshold for PoC comparison
  # Early termination parameters
  enable_early_termination = TRUE,
  log_early_termination = TRUE
)

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

p_YI = c(0.2, 0.4, 0.6)
rho0 = 1.5
rho1 = 2

# Utility table
# Rows: Efficacy (0, 1)
# Columns: Toxicity (0, 1)
# Slices: Immune Response (0, 1)
utility_table <- array(0, dim = c(2, 2, 2), dimnames = list(
  E = c(0, 1),
  T = c(0, 1),
  I = c(0, 1)
))

utility_table[1, 1, 1] <- 0   # E=0, T=0, I=0
utility_table[2, 1, 1] <- 80  # E=1, T=0, I=0
utility_table[1, 2, 1] <- 0   # E=0, T=1, I=0
utility_table[2, 2, 1] <- 30  # E=1, T=1, I=0

utility_table[1, 1, 2] <- 10  # E=0, T=0, I=1
utility_table[2, 1, 2] <- 100 # E=1, T=0, I=1
utility_table[1, 2, 2] <- 0   # E=0, T=1, I=1
utility_table[2, 2, 2] <- 40  # E=1, T=1, I=1

# Add utility table to trial_config
trial_config$utility_table <- utility_table