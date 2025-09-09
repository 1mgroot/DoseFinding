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
  c_POC = 0.9
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
