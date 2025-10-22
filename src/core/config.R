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
  phi_E = 0.2, c_E = 0.75,  # Adjusted from 0.9 for calibration space
  phi_I = 0.1, c_I = 0.65,  # Adjusted from 0.7 for calibration space
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

# Add correlation parameters to trial_config
trial_config$rho0 <- rho0
trial_config$rho1 <- rho1

# Flat Scenario Parameters for Calibration
# These parameters define the null scenario for PoC calibration
flat_scenario_config <- list(
  # Scenario type
  scenario_type = "flat_null",
  
  # Lower bound parameters (from meeting requirements)
  # These define the null hypothesis: all doses are equally unfavorable
  # Using conditional probabilities that will result in appropriate marginal probabilities
  phi_I_lower = 0.25,  # Immune response rate for all doses (increased to improve admissibility)
  phi_E_lower = 0.20,  # Marginal efficacy rate for all doses
  toxicity_low = 0.02,  # Conditional toxicity rate for all doses (reduced to account for copula effect)
  
  # Trial configuration (inherited from main config)
  dose_levels = trial_config$dose_levels,
  n_stages = trial_config$n_stages,
  cohort_size = trial_config$cohort_size,
  
  # Thresholds (to be calibrated)
  phi_T = trial_config$phi_T,
  c_T = trial_config$c_T,
  phi_E = trial_config$phi_E,
  c_E = trial_config$c_E,
  phi_I = trial_config$phi_I,
  c_I = trial_config$c_I,
  
  # PoC parameters (to be calibrated)
  c_poc = trial_config$c_poc,
  delta_poc = trial_config$delta_poc,
  
  # Early termination parameters
  enable_early_termination = trial_config$enable_early_termination,
  log_early_termination = trial_config$log_early_termination,
  
  # Correlation parameters
  rho0 = trial_config$rho0,
  rho1 = trial_config$rho1,
  
  # Utility table
  utility_table = trial_config$utility_table
)

# Unfavorable Scenario Parameters for Early Termination Calibration
# These parameters define an unfavorable scenario for early termination calibration
unfavorable_scenario_config <- list(
  # Scenario type
  scenario_type = "unfavorable",
  
  # Unfavorable probabilities (more moderate to allow meaningful calibration)
  phi_I_lower = 0.15,  # Moderate immune response rate
  phi_E_lower = 0.25,  # Moderate efficacy rate (not completely inactive)
  toxicity_low = 0.15,  # Moderate toxicity rate (to allow some doses through)
  
  # Trial configuration (inherited from main config)
  dose_levels = trial_config$dose_levels,
  n_stages = trial_config$n_stages,
  cohort_size = trial_config$cohort_size,
  
  # Thresholds (to be calibrated)
  phi_T = trial_config$phi_T,
  c_T = trial_config$c_T,
  phi_E = trial_config$phi_E,
  c_E = trial_config$c_E,
  phi_I = trial_config$phi_I,
  c_I = trial_config$c_I,
  
  # PoC parameters (to be calibrated)
  c_poc = trial_config$c_poc,
  delta_poc = trial_config$delta_poc,
  
  # Early termination parameters
  enable_early_termination = trial_config$enable_early_termination,
  log_early_termination = trial_config$log_early_termination,
  
  # Correlation parameters
  rho0 = trial_config$rho0,
  rho1 = trial_config$rho1,
  
  # Utility table
  utility_table = trial_config$utility_table
)

# Calibration Parameters
calibration_config <- list(
  # PoC calibration targets
  poc_target_rate = 0.10,  # Target 10% PoC detection rate
  poc_tolerance = 0.02,     # ±2% tolerance
  
  # Early termination calibration targets
  early_termination_target_rate = 0.80,  # Target 80% early termination rate
  early_termination_tolerance = 0.05,    # ±5% tolerance
  
  # Simulation parameters
  n_calibration_simulations = 10000,  # Number of simulations for calibration
  n_validation_simulations = 5000,    # Number of simulations for validation
  
  # Parameter search ranges
  c_poc_range = seq(0.5, 0.99, by = 0.01),
  c_T_range = seq(0.7, 0.99, by = 0.01),
  c_E_range = seq(0.7, 0.99, by = 0.01),
  c_I_range = seq(0.6, 0.99, by = 0.01)
)