# Quick Parameter Optimization Runner
# This script provides a simple interface to run parameter optimization

# Load the optimization system
source("parameter_optimization.R")

# Function to run a quick optimization
quick_optimization <- function(n_combinations = 20, n_simulations = 3) {
  cat("=== Quick Parameter Optimization ===\n")
  cat("This will test", n_combinations, "parameter combinations\n")
  cat("with", n_simulations, "simulations each.\n")
  cat("Estimated time: ~", n_combinations * n_simulations * 2, "minutes\n\n")
  
  # Run optimization
  results <- run_parameter_optimization(n_combinations, n_simulations)
  
  # Create visualizations
  plots <- create_optimization_plots(results)
  
  # Find best parameters
  best_params <- find_best_parameters(results)
  
  # Save results
  saveRDS(results, "quick_optimization_results.rds")
  
  # Display summary
  cat("\n=== Optimization Summary ===\n")
  cat("Best mean final utility:", round(best_params$mean_final_utility, 2), "\n")
  cat("Best completion rate:", round(best_params$completion_rate, 3), "\n")
  cat("Best correct selection rate:", round(best_params$correct_selection_rate, 3), "\n")
  
  return(list(
    results = results,
    plots = plots,
    best_params = best_params
  ))
}

# Function to run comprehensive optimization
comprehensive_optimization <- function(n_combinations = 50, n_simulations = 5) {
  cat("=== Comprehensive Parameter Optimization ===\n")
  cat("This will test", n_combinations, "parameter combinations\n")
  cat("with", n_simulations, "simulations each.\n")
  cat("Estimated time: ~", n_combinations * n_simulations * 3, "minutes\n\n")
  
  # Run optimization
  results <- run_parameter_optimization(n_combinations, n_simulations)
  
  # Create visualizations
  plots <- create_optimization_plots(results)
  
  # Find best parameters
  best_params <- find_best_parameters(results)
  
  # Save results
  saveRDS(results, "comprehensive_optimization_results.rds")
  saveRDS(plots, "comprehensive_optimization_plots.rds")
  
  # Display summary
  cat("\n=== Optimization Summary ===\n")
  cat("Best mean final utility:", round(best_params$mean_final_utility, 2), "\n")
  cat("Best completion rate:", round(best_params$completion_rate, 3), "\n")
  cat("Best correct selection rate:", round(best_params$correct_selection_rate, 3), "\n")
  
  return(list(
    results = results,
    plots = plots,
    best_params = best_params
  ))
}

# Function to analyze existing results
analyze_results <- function(results_file = "optimization_results.rds") {
  if (file.exists(results_file)) {
    results <- readRDS(results_file)
    
    # Create visualizations
    plots <- create_optimization_plots(results)
    
    # Find best parameters
    best_params <- find_best_parameters(results)
    
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
  } else {
    cat("Results file not found:", results_file, "\n")
    cat("Please run optimization first.\n")
    return(NULL)
  }
}

# Function to test specific parameter combination
test_specific_params <- function(phi_T = 0.35, phi_E = 0.10, phi_I = 0.20,
                               c_T = 0.8, c_E = 0.7, c_I = 0.7,
                               utility_type = "balanced", n_simulations = 10) {
  cat("=== Testing Specific Parameters ===\n")
  cat("φ_T =", phi_T, "φ_E =", phi_E, "φ_I =", phi_I, "\n")
  cat("c_T =", c_T, "c_E =", c_E, "c_I =", c_I, "\n")
  cat("Utility type:", utility_type, "\n\n")
  
  # Create parameter grids
  grids <- create_parameter_grids()
  
  # Create params list
  params <- list(
    phi_T = phi_T,
    phi_E = phi_E,
    phi_I = phi_I,
    c_T = c_T,
    c_E = c_E,
    c_I = c_I,
    utility_type = utility_type,
    utility_table = grids$utility_variations[[utility_type]]$utility_table
  )
  
  # Define data parameters
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
  
  # Run parameter combination
  result <- run_parameter_combination(params, data_params, n_simulations)
  
  # Display results
  cat("=== Results ===\n")
  cat("Completion rate:", round(result$aggregated_metrics$completion_rate, 3), "\n")
  cat("Correct selection rate:", round(result$aggregated_metrics$correct_selection_rate, 3), "\n")
  cat("Mean final utility:", round(result$aggregated_metrics$mean_final_utility, 2), "\n")
  cat("PoC success rate:", round(result$aggregated_metrics$poc_success_rate, 3), "\n")
  cat("Mean total participants:", round(result$aggregated_metrics$mean_total_participants, 1), "\n")
  cat("Mean allocation efficiency:", round(result$aggregated_metrics$mean_allocation_efficiency, 3), "\n")
  
  return(result)
}

# Main execution - uncomment one of these lines to run:

# 1. Quick optimization (recommended for first run)
# quick_optimization()

# 2. Comprehensive optimization (more thorough but takes longer)
# comprehensive_optimization()

# 3. Test specific parameters
# test_specific_params(phi_T = 0.35, phi_E = 0.10, phi_I = 0.20, c_T = 0.8, c_E = 0.7, c_I = 0.7)

# 4. Analyze existing results
# analyze_results()

cat("Parameter optimization system loaded.\n")
cat("To run optimization, uncomment one of the function calls at the bottom of this file.\n")
cat("Recommended: Start with quick_optimization()\n")
