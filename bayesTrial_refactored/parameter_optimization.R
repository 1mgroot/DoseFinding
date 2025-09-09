# Parameter Optimization System for Bayesian Dose-Finding Trial
# This script provides a systematic approach to optimize trial parameters

library(dplyr)
library(ggplot2)
library(gridExtra)
library(purrr)
library(tidyr)

# Source required functions
source("config.R")
source("helpers.R")
source("simulate_data.R")
source("model_utils.R")
source("dose_decision.R")
source("main.R")

# Define parameter grids for optimization
create_parameter_grids <- function() {
  # Threshold parameters (phi values)
  phi_grids <- list(
    # Toxicity threshold - conservative to very permissive
    phi_T = c(0.25, 0.30, 0.35, 0.40, 0.45),
    # Efficacy threshold - very low to moderate
    phi_E = c(0.05, 0.10, 0.15, 0.20, 0.25),
    # Immune response threshold - moderate to high
    phi_I = c(0.15, 0.20, 0.25, 0.30, 0.35)
  )
  
  # Credibility parameters (c values)
  c_grids <- list(
    # Toxicity credibility - conservative to moderate
    c_T = c(0.7, 0.8, 0.9, 0.95),
    # Efficacy credibility - moderate to high
    c_E = c(0.6, 0.7, 0.8, 0.9),
    # Immune response credibility - moderate to high
    c_I = c(0.6, 0.7, 0.8, 0.9)
  )
  
  # Utility table variations
  utility_variations <- list(
    # Conservative utility (lower rewards for toxicity)
    conservative = list(
      utility_table = array(c(
        0, 80, 0, 20,   # I=0: E=0,1 x T=0,1
        10, 100, 0, 30  # I=1: E=0,1 x T=0,1
      ), dim = c(2, 2, 2))
    ),
    # Balanced utility (current settings)
    balanced = list(
      utility_table = array(c(
        0, 80, 0, 30,   # I=0: E=0,1 x T=0,1
        10, 100, 0, 40  # I=1: E=0,1 x T=0,1
      ), dim = c(2, 2, 2))
    ),
    # Aggressive utility (higher rewards, less penalty for toxicity)
    aggressive = list(
      utility_table = array(c(
        0, 80, 0, 40,   # I=0: E=0,1 x T=0,1
        10, 100, 0, 50  # I=1: E=0,1 x T=0,1
      ), dim = c(2, 2, 2))
    )
  )
  
  return(list(
    phi_grids = phi_grids,
    c_grids = c_grids,
    utility_variations = utility_variations
  ))
}

# Generate parameter combinations
generate_parameter_combinations <- function(grids, n_combinations = 100) {
  # Create all possible combinations
  all_combinations <- expand.grid(
    phi_T = grids$phi_grids$phi_T,
    phi_E = grids$phi_grids$phi_E,
    phi_I = grids$phi_grids$phi_I,
    c_T = grids$c_grids$c_T,
    c_E = grids$c_grids$c_E,
    c_I = grids$c_grids$c_I,
    utility_type = names(grids$utility_variations)
  )
  
  # Sample if too many combinations
  if (nrow(all_combinations) > n_combinations) {
    set.seed(123)
    all_combinations <- all_combinations[sample(nrow(all_combinations), n_combinations), ]
  }
  
  return(all_combinations)
}

# Evaluation metrics for parameter performance
calculate_evaluation_metrics <- function(results, true_optimal_dose = 4) {
  metrics <- list()
  
  # 1. Trial completion rate
  metrics$completion_rate <- !results$terminated_early
  
  # 2. Optimal dose selection accuracy
  if (!results$terminated_early) {
    metrics$correct_selection <- results$final_od == true_optimal_dose
    metrics$selected_dose <- results$final_od
  } else {
    metrics$correct_selection <- FALSE
    metrics$selected_dose <- NA
  }
  
  # 3. Final utility score
  metrics$final_utility <- ifelse(!results$terminated_early, results$final_utility, 0)
  
  # 4. PoC validation success
  metrics$poc_success <- ifelse(!results$terminated_early, results$poc_validated, FALSE)
  
  # 5. Total participants used
  metrics$total_participants <- nrow(results$all_data)
  
  # 6. Allocation efficiency (how well participants were distributed)
  allocation_summary <- results$all_data %>%
    group_by(d) %>%
    summarise(n_participants = n(), .groups = 'drop')
  
  # Calculate allocation efficiency (prefer more balanced allocation)
  if (nrow(allocation_summary) > 0) {
    metrics$allocation_efficiency <- 1 - sd(allocation_summary$n_participants) / mean(allocation_summary$n_participants)
  } else {
    metrics$allocation_efficiency <- 0
  }
  
  # 7. Early termination stage (if applicable)
  metrics$termination_stage <- ifelse(results$terminated_early, results$termination_stage, NA)
  
  return(metrics)
}

# Run single parameter combination
run_parameter_combination <- function(params, data_params, n_simulations = 10) {
  # Create trial configuration
  config <- list(
    dose_levels = c(1, 2, 3, 4, 5),
    n_stages = 5,
    cohort_size = 15,
    phi_T = params$phi_T,
    c_T = params$c_T,
    phi_E = params$phi_E,
    c_E = params$c_E,
    phi_I = params$phi_I,
    c_I = params$c_I,
    c_poc = 0.9,
    delta_poc = 0.8,
    enable_early_termination = TRUE,
    log_early_termination = FALSE,  # Reduce output for batch runs
    utility_table = params$utility_table
  )
  
  # Run multiple simulations for robustness
  simulation_results <- list()
  metrics_summary <- list()
  
  for (i in 1:n_simulations) {
    tryCatch({
      results <- run_trial_simulation(config, data_params$p_YI, data_params$p_YT_given_I, 
                                    data_params$p_YE_given_I, data_params$rho0, data_params$rho1)
      metrics <- calculate_evaluation_metrics(results)
      
      simulation_results[[i]] <- results
      metrics_summary[[i]] <- metrics
    }, error = function(e) {
      cat("Error in simulation", i, ":", e$message, "\n")
      # Return default metrics for failed simulation
      metrics_summary[[i]] <<- list(
        completion_rate = FALSE,
        correct_selection = FALSE,
        selected_dose = NA,
        final_utility = 0,
        poc_success = FALSE,
        total_participants = 0,
        allocation_efficiency = 0,
        termination_stage = 1
      )
    })
  }
  
  # Aggregate metrics across simulations
  aggregated_metrics <- list(
    completion_rate = mean(sapply(metrics_summary, function(x) x$completion_rate)),
    correct_selection_rate = mean(sapply(metrics_summary, function(x) x$correct_selection)),
    mean_final_utility = mean(sapply(metrics_summary, function(x) x$final_utility)),
    poc_success_rate = mean(sapply(metrics_summary, function(x) x$poc_success)),
    mean_total_participants = mean(sapply(metrics_summary, function(x) x$total_participants)),
    mean_allocation_efficiency = mean(sapply(metrics_summary, function(x) x$allocation_efficiency)),
    mean_termination_stage = mean(sapply(metrics_summary, function(x) x$termination_stage), na.rm = TRUE)
  )
  
  return(list(
    params = params,
    aggregated_metrics = aggregated_metrics,
    individual_results = simulation_results,
    individual_metrics = metrics_summary
  ))
}

# Main optimization function
run_parameter_optimization <- function(n_combinations = 50, n_simulations = 5) {
  cat("Starting parameter optimization...\n")
  cat("Number of parameter combinations:", n_combinations, "\n")
  cat("Simulations per combination:", n_simulations, "\n")
  
  # Create parameter grids
  grids <- create_parameter_grids()
  
  # Generate parameter combinations
  param_combinations <- generate_parameter_combinations(grids, n_combinations)
  
  # Define data generation parameters (from simulation_notebook.qmd)
  data_params <- list(
    p_YI = c(0.10, 0.30, 0.50, 0.60, 0.70),
    p_YT_given_I = matrix(c(
      0.05, 0.10, 0.12, 0.18, 0.25,
      0.08, 0.12, 0.15, 0.25, 0.35
    ), nrow = 5, ncol = 2),
    p_YE_given_I = matrix(c(
      0.10, 0.20, 0.35, 0.45, 0.50,
      0.30, 0.50, 0.70, 0.80, 0.75
    ), nrow = 5, ncol = 2),
    rho0 = 1.5,
    rho1 = 2
  )
  
  # Add utility tables to parameter combinations
  param_combinations$utility_table <- lapply(param_combinations$utility_type, function(type) {
    grids$utility_variations[[type]]$utility_table
  })
  
  # Run optimization
  results <- list()
  for (i in 1:nrow(param_combinations)) {
    cat("Running combination", i, "of", nrow(param_combinations), "\n")
    
    params <- as.list(param_combinations[i, ])
    params$utility_table <- params$utility_table[[1]]  # Extract from list
    
    result <- run_parameter_combination(params, data_params, n_simulations)
    results[[i]] <- result
  }
  
  return(results)
}

# Visualization functions
create_optimization_plots <- function(optimization_results) {
  # Extract aggregated metrics
  metrics_df <- do.call(rbind, lapply(optimization_results, function(result) {
    c(result$params, result$aggregated_metrics)
  }))
  
  # Convert to data frame
  metrics_df <- as.data.frame(metrics_df)
  
  # Create plots
  plots <- list()
  
  # 1. Completion rate vs parameters
  p1 <- ggplot(metrics_df, aes(x = phi_T, y = completion_rate, color = factor(c_T))) +
    geom_point(size = 3) +
    facet_grid(phi_E ~ phi_I, labeller = label_both) +
    labs(title = "Trial Completion Rate vs Parameters",
         x = "Toxicity Threshold (φ_T)", y = "Completion Rate",
         color = "Toxicity Credibility (c_T)") +
    theme_minimal()
  
  # 2. Correct selection rate vs parameters
  p2 <- ggplot(metrics_df, aes(x = phi_T, y = correct_selection_rate, color = factor(c_T))) +
    geom_point(size = 3) +
    facet_grid(phi_E ~ phi_I, labeller = label_both) +
    labs(title = "Correct Dose Selection Rate vs Parameters",
         x = "Toxicity Threshold (φ_T)", y = "Correct Selection Rate",
         color = "Toxicity Credibility (c_T)") +
    theme_minimal()
  
  # 3. Mean utility vs parameters
  p3 <- ggplot(metrics_df, aes(x = phi_T, y = mean_final_utility, color = factor(c_T))) +
    geom_point(size = 3) +
    facet_grid(phi_E ~ phi_I, labeller = label_both) +
    labs(title = "Mean Final Utility vs Parameters",
         x = "Toxicity Threshold (φ_T)", y = "Mean Final Utility",
         color = "Toxicity Credibility (c_T)") +
    theme_minimal()
  
  # 4. Parameter importance analysis
  importance_df <- metrics_df %>%
    select(phi_T, phi_E, phi_I, c_T, c_E, c_I, completion_rate, correct_selection_rate, mean_final_utility) %>%
    gather(key = "parameter", value = "value", -completion_rate, -correct_selection_rate, -mean_final_utility)
  
  p4 <- ggplot(importance_df, aes(x = value, y = mean_final_utility, color = parameter)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE) +
    facet_wrap(~parameter, scales = "free_x") +
    labs(title = "Parameter Sensitivity Analysis",
         x = "Parameter Value", y = "Mean Final Utility") +
    theme_minimal()
  
  # 5. Best parameter combinations
  best_combinations <- metrics_df %>%
    arrange(desc(mean_final_utility)) %>%
    head(10)
  
  p5 <- ggplot(best_combinations, aes(x = reorder(paste("Combination", 1:nrow(best_combinations)), mean_final_utility), 
                                      y = mean_final_utility, fill = utility_type)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Top 10 Parameter Combinations by Utility",
         x = "Parameter Combination", y = "Mean Final Utility",
         fill = "Utility Type") +
    theme_minimal()
  
  return(list(
    completion_rate = p1,
    selection_rate = p2,
    utility = p3,
    sensitivity = p4,
    top_combinations = p5,
    metrics_df = metrics_df
  ))
}

# Function to find best parameters
find_best_parameters <- function(optimization_results, metric = "mean_final_utility") {
  # Extract aggregated metrics
  metrics_df <- do.call(rbind, lapply(optimization_results, function(result) {
    c(result$params, result$aggregated_metrics)
  }))
  
  # Find best combination based on specified metric
  best_idx <- which.max(metrics_df[[metric]])
  best_params <- metrics_df[best_idx, ]
  
  cat("Best parameters found:\n")
  cat("Metric used:", metric, "\n")
  cat("Toxicity threshold (φ_T):", best_params$phi_T, "\n")
  cat("Efficacy threshold (φ_E):", best_params$phi_E, "\n")
  cat("Immune response threshold (φ_I):", best_params$phi_I, "\n")
  cat("Toxicity credibility (c_T):", best_params$c_T, "\n")
  cat("Efficacy credibility (c_E):", best_params$c_E, "\n")
  cat("Immune response credibility (c_I):", best_params$c_I, "\n")
  cat("Utility type:", best_params$utility_type, "\n")
  cat("Mean final utility:", best_params$mean_final_utility, "\n")
  cat("Completion rate:", best_params$completion_rate, "\n")
  cat("Correct selection rate:", best_params$correct_selection_rate, "\n")
  
  return(best_params)
}

# Main execution function
main_optimization <- function() {
  # Run optimization
  results <- run_parameter_optimization(n_combinations = 30, n_simulations = 3)
  
  # Create visualizations
  plots <- create_optimization_plots(results)
  
  # Find best parameters
  best_params <- find_best_parameters(results)
  
  # Save results
  saveRDS(results, "optimization_results.rds")
  saveRDS(plots, "optimization_plots.rds")
  
  # Display plots
  print(plots$completion_rate)
  print(plots$selection_rate)
  print(plots$utility)
  print(plots$sensitivity)
  print(plots$top_combinations)
  
  return(list(
    results = results,
    plots = plots,
    best_params = best_params
  ))
}
