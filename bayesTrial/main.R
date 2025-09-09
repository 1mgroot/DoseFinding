source("config.R")
source("helpers.R")
source("simulate_data.R")
source("model_utils.R")
source("dose_decision.R")

# Simulate data
df <- simulate_data_gumbel(
  n_per_dose = 100,
  dose_levels = trial_config$dose_levels,
  p_YI = p_YI,
  p_YT_given_I = p_YT_given_I,
  p_YE_given_I = p_YE_given_I,
  rho0 = rho0,
  rho1 = rho1
)

# Analyze immune response
pi_I_stats <- compute_rn(df, outcome_col = "Y_I")
pi_I_post <- simulate_beta_posterior(pi_I_stats)
pi_I_post <- add_beta_variance(pi_I_post)
pi_I_pava <- apply_pava_on_samples(pi_I_post)

# Analyze toxicity
tox_stats <- compute_rn(df, outcome_col = "Y_T", group_col = "Y_I")
tox_post <- simulate_beta_posterior(tox_stats)
tox_post <- add_beta_variance(tox_post)
tox_pava <- apply_biviso_on_matrix(tox_post)

# Analyze efficacy
eff_stats <- compute_rn(df, outcome_col = "Y_E", group_col = "Y_I")
eff_post <- simulate_beta_posterior(eff_stats)
eff_post <- add_beta_variance(eff_post)
eff_pava <- apply_biviso_on_matrix(eff_post)

# Get admissible set
admissible_set <- get_admissible_set(NULL, trial_config)

# Get final OD
final_od <- select_final_od(admissible_set, NULL, trial_config)

# Print results
print("Final OD:")
print(final_od)

# Plot results
plot_posterior_summary(pi_I_pava, title = "Immune Response vs Dose (PAVA Adjusted)", file_path = "results/immune_response.png")
plot_posterior_summary(tox_pava, title = "Toxicity Rate by Dose and Immune Status", group_col = "Y_I", file_path = "results/toxicity.png")
plot_posterior_summary(eff_pava, title = "Efficacy Rate by Dose and Immune Status", group_col = "Y_I", file_path = "results/efficacy.png")
