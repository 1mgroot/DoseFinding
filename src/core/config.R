library(dplyr)
library(tidyr)
library(isotone)
library(purrr)
library(ggplot2)
library(Iso)

# Trial configuration - aligned with simulation_notebook.qmd
trial_config <- list(
  dose_levels = c(1, 2, 3, 4, 5),
  n_stages = 5,
  cohort_size = 15,
  phi_T = 0.35, c_T = 0.5,   # Toxicity threshold and probability cutoff
  phi_E = 0.1,  c_E = 0.5,   # Efficacy threshold and probability cutoff
  phi_I = 0.20, c_I = 0.5,   # Immune response threshold and probability cutoff
  # PoC parameters
  c_poc = 0.9,
  delta_poc = 0.8,  # Threshold for PoC comparison
  # Early termination parameters
  enable_early_termination = TRUE,
  log_early_termination = TRUE
)

# Data simulation parameters (5-dose design with increasing dose-response)
p_YT_given_I <- matrix(c(
  # I=0 (No Immune Response)
  0.05, 0.10, 0.12, 0.18, 0.25,
  # I=1 (Immune Response)
  0.08, 0.12, 0.15, 0.25, 0.35
), nrow = 5, ncol = 2)

p_YE_given_I <- matrix(c(
  # I=0 (No Immune Response)
  0.10, 0.20, 0.35, 0.45, 0.50, 
  # I=1 (Immune Response)
  0.30, 0.50, 0.70, 0.80, 0.75  
), nrow = 5, ncol = 2)

p_YI <- c(0.10, 0.30, 0.50, 0.60, 0.70)  # Immune response probability per dose
rho0 <- 1.5  # Gumbel copula correlation under I=0
rho1 <- 2    # Gumbel copula correlation under I=1

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