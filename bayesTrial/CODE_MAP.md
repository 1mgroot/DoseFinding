## Code Origin Map

This document maps the code from the old R scripts to the new files in the `bayesTrial` project.

### `config.R`

- The `trial_config` list is new, based on the provided documentation.
- The `p_YT_given_I`, `p_YE_given_I`, `p_YI`, `rho0`, and `rho1` variables are from the bottom of the `addBiviso.R` script.

### `helpers.R`

- `plot_posterior_summary()`: This function is from `addBiviso.R` and `simulateWithGumbelDiff.r`.

### `simulate_data.R`

- `Gumbel()`: This function is from `addBiviso.R` and `SimulateWithGumbel.R`.
- `simulate_data_gumbel()`: This function is from `addBiviso.R` and `SimulateWithGumbel.R`.

### `model_utils.R`

- `compute_rn()`: This function is from `BQDTwist.R`, `addBiviso.R`, and `test_model.R`.
- `simulate_beta_posterior()`: This function is from `BQDTwist.R`, `addBiviso.R`, and `test_model.R`.
- `add_beta_variance()`: This function is from `BQDTwist.R`, `addBiviso.R`, and `test_model.R`.
- `apply_pava_on_samples()`: This function is from `BQDTwist.R`, `addBiviso.R`, and `test_model.R`.
- `apply_biviso_on_matrix()`: This function is from `addBiviso.R` and `test_model.R`.
- `compute_marginal_probability()`: This function is from `addBiviso.R` and `test_model.R`.

### `dose_decision.R`

- All functions in this file are new placeholders based on the provided documentation.

### `main.R`

- The code in this file is a combination of the analysis steps from the bottom of `addBiviso.R` and the new dose decision functions.

### `test_main.R`

- This is a new file to test the `main.R` script.
