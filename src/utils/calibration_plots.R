# Calibration Visualization Functions
# This file implements visualization tools for calibration results

# Load required libraries
library(ggplot2)
library(dplyr)

plot_poc_calibration_curve <- function(calibration_results, target_rate = 0.10, save_path = NULL) {
  # Plot PoC calibration curve showing C_poc vs PoC detection rate
  #
  # Args:
  #   calibration_results: Results from calibrate_c_poc function
  #   target_rate: Target PoC detection rate (default: 0.10)
  #   save_path: Optional path to save the plot
  #
  # Returns:
  #   ggplot object: Calibration curve plot
  
  # Extract calibration data
  if (is.list(calibration_results) && "calibration_results" %in% names(calibration_results)) {
    data <- calibration_results$calibration_results
    optimal_c_poc <- calibration_results$optimal_c_poc
  } else {
    data <- calibration_results
    optimal_c_poc <- data$c_poc[which.min(abs(data$poc_detection_rate - target_rate))]
  }
  
  p <- ggplot(data, aes(x = .data$c_poc)) +
    geom_line(aes(y = .data$poc_detection_rate), color = "#009E73", size = 1.2) +
    geom_ribbon(aes(ymin = .data$poc_detection_rate_lower, ymax = .data$poc_detection_rate_upper), 
                alpha = 0.3, fill = "#009E73") +
    geom_hline(yintercept = target_rate, linetype = "dashed", color = "#E69F00", size = 1) +
    geom_vline(xintercept = optimal_c_poc, 
               linetype = "dashed", color = "#56B4E9", size = 1) +
    geom_point(aes(x = optimal_c_poc, y = data$poc_detection_rate[data$c_poc == optimal_c_poc]), 
               color = "#CC79A7", size = 4, shape = 16) +
    labs(
      title = "PoC Calibration Curve",
      subtitle = paste("Target detection rate:", target_rate, "| Optimal C_poc =", round(optimal_c_poc, 3)),
      x = "C_poc Threshold",
      y = "PoC Detection Rate",
      caption = paste("Based on", unique(data$n_simulations), "simulations per point")
    ) +
    scale_x_continuous(breaks = seq(0.5, 1.0, by = 0.1)) +
    scale_y_continuous(breaks = seq(0, 0.3, by = 0.05), limits = c(0, max(0.3, max(data$poc_detection_rate_upper, na.rm = TRUE)))) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      plot.caption = element_text(size = 10, hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", size = 0.5)
    )
  
  if (!is.null(save_path)) {
    # Create directory if it doesn't exist
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300)
    cat("PoC calibration curve saved to:", save_path, "\n")
  }
  
  return(p)
}

plot_early_termination_curve <- function(termination_results, target_rate = 0.80, save_path = NULL) {
  # Plot early termination calibration curve
  #
  # Args:
  #   termination_results: Results from early termination calibration
  #   target_rate: Target early termination rate (default: 0.80)
  #   save_path: Optional path to save the plot
  #
  # Returns:
  #   ggplot object: Early termination curve plot
  
  # Extract termination data
  if (is.list(termination_results) && "calibration_results" %in% names(termination_results)) {
    data <- termination_results$calibration_results
    optimal_threshold <- termination_results$optimal_threshold
  } else {
    data <- termination_results
    optimal_threshold <- data$threshold[which.min(abs(data$termination_rate - target_rate))]
  }
  
  has_threshold <- "threshold" %in% names(data)
  has_grid <- all(c("c_T", "c_E") %in% names(data))
  
  if (has_threshold) {
    p <- ggplot(data, aes(x = .data$threshold)) +
      geom_line(aes(y = .data$termination_rate), color = "#E69F00", size = 1.2) +
      geom_ribbon(aes(ymin = .data$termination_rate_lower, ymax = .data$termination_rate_upper), 
                  alpha = 0.3, fill = "#E69F00") +
      geom_hline(yintercept = target_rate, linetype = "dashed", color = "#009E73", size = 1) +
      geom_vline(xintercept = optimal_threshold, 
                 linetype = "dashed", color = "#56B4E9", size = 1) +
      geom_point(aes(x = optimal_threshold, y = data$termination_rate[data$threshold == optimal_threshold]), 
                 color = "#CC79A7", size = 4, shape = 16) +
      labs(
        title = "Early Termination Calibration Curve",
        subtitle = paste("Target termination rate:", target_rate, "| Optimal threshold =", round(optimal_threshold, 3)),
        x = "Termination Threshold",
        y = "Early Termination Rate",
        caption = paste("Based on", unique(data$n_simulations), "simulations per point")
      ) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
      theme_bw(base_size = 16) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.caption = element_text(size = 10, hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", size = 0.5)
      )
  } else if (has_grid) {
    optimal_c_T <- if (!is.null(termination_results$optimal_c_T)) termination_results$optimal_c_T else NA_real_
    optimal_c_E <- if (!is.null(termination_results$optimal_c_E)) termination_results$optimal_c_E else NA_real_
    threshold_type <- if (!is.null(termination_results$threshold_type)) termination_results$threshold_type else "c_T and c_E"
    
    caption_text <- NULL
    if ("n_simulations" %in% names(data)) {
      caption_text <- paste("Based on", unique(data$n_simulations), "simulations per point")
    }
    
    p <- ggplot(data, aes(x = .data$c_T, y = .data$c_E, fill = .data$termination_rate)) +
      geom_tile(color = "white") +
      geom_point(
        data = data.frame(c_T = optimal_c_T, c_E = optimal_c_E),
        aes(x = .data$c_T, y = .data$c_E),
        inherit.aes = FALSE,
        color = "#CC79A7",
        size = 4
      ) +
      labs(
        title = "Early Termination Calibration Heatmap",
        subtitle = paste("Target termination rate:", target_rate, "| Optimal", threshold_type),
        x = "c_T Threshold",
        y = "c_E Threshold",
        fill = "Termination Rate",
        caption = caption_text
      ) +
      scale_fill_gradient(low = "#F0F0F0", high = "#E69F00", limits = c(0, 1)) +
      theme_bw(base_size = 16) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.caption = element_text(size = 10, hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", size = 0.5)
      )
  } else {
    stop("Early termination calibration data must include 'threshold' or 'c_T' and 'c_E' columns.")
  }
  
  if (!is.null(save_path)) {
    # Create directory if it doesn't exist
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300)
    cat("Early termination curve saved to:", save_path, "\n")
  }
  
  return(p)
}

plot_threshold_performance_curves <- function(calibration_data, save_path = NULL) {
  # Plot combined threshold vs performance curves for multiple metrics
  #
  # Args:
  #   calibration_data: List containing both PoC and termination calibration results
  #   save_path: Optional path to save the plot
  #
  # Returns:
  #   ggplot object: Combined performance curves plot
  
  # Prepare data for combined plot
  combined_data <- data.frame()
  
  # Add PoC calibration data
  if ("poc_calibration" %in% names(calibration_data)) {
    poc_data <- calibration_data$poc_calibration$calibration_results
    poc_data$metric <- "PoC Detection Rate"
    poc_data$threshold <- poc_data$c_poc
    poc_data$rate <- poc_data$poc_detection_rate
    poc_data$rate_lower <- poc_data$poc_detection_rate_lower
    poc_data$rate_upper <- poc_data$poc_detection_rate_upper
    combined_data <- rbind(combined_data, poc_data[, c("threshold", "rate", "rate_lower", "rate_upper", "metric")])
  }
  
  # Add termination calibration data
  if ("termination_calibration" %in% names(calibration_data)) {
    term_data <- calibration_data$termination_calibration$calibration_results
    
    if ("threshold" %in% names(term_data)) {
      term_data$metric <- "Early Termination Rate"
      term_data$threshold <- term_data$threshold
      term_data$rate <- term_data$termination_rate
      term_data$rate_lower <- term_data$termination_rate_lower
      term_data$rate_upper <- term_data$termination_rate_upper
      combined_data <- rbind(combined_data, term_data[, c("threshold", "rate", "rate_lower", "rate_upper", "metric")])
    } else if (all(c("c_T", "c_E") %in% names(term_data))) {
      optimal_c_E <- calibration_data$termination_calibration$optimal_c_E
      if (is.null(optimal_c_E)) {
        optimal_c_E <- term_data$c_E[1]
      }
      
      term_slice <- term_data[term_data$c_E == optimal_c_E, , drop = FALSE]
      term_slice$metric <- "Early Termination Rate"
      term_slice$threshold <- term_slice$c_T
      term_slice$rate <- term_slice$termination_rate
      term_slice$rate_lower <- if ("termination_rate_lower" %in% names(term_slice)) {
        term_slice$termination_rate_lower
      } else {
        term_slice$termination_rate
      }
      term_slice$rate_upper <- if ("termination_rate_upper" %in% names(term_slice)) {
        term_slice$termination_rate_upper
      } else {
        term_slice$termination_rate
      }
      
      combined_data <- rbind(combined_data, term_slice[, c("threshold", "rate", "rate_lower", "rate_upper", "metric")])
    }
  }
  
  if (nrow(combined_data) == 0) {
    stop("No calibration data provided")
  }
  
  # Create combined plot
  p <- ggplot(combined_data, aes(x = .data$threshold, y = .data$rate, color = .data$metric)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = .data$rate_lower, ymax = .data$rate_upper, fill = .data$metric), 
                alpha = 0.3) +
    labs(
      title = "Threshold vs Performance Curves",
      subtitle = "Calibration Results for PoC Detection and Early Termination",
      x = "Threshold Parameter",
      y = "Performance Rate",
      color = "Metric",
      fill = "Metric"
    ) +
    scale_color_manual(values = c("PoC Detection Rate" = "#009E73", "Early Termination Rate" = "#E69F00")) +
    scale_fill_manual(values = c("PoC Detection Rate" = "#009E73", "Early Termination Rate" = "#E69F00")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", size = 0.5)
    )
  
  if (!is.null(save_path)) {
    # Create directory if it doesn't exist
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(save_path, plot = p, width = 12, height = 8, dpi = 300)
    cat("Combined performance curves saved to:", save_path, "\n")
  }
  
  return(p)
}

plot_calibration_summary <- function(calibration_results, save_path = NULL) {
  # Create a summary plot showing calibration results and validation
  #
  # Args:
  #   calibration_results: Results from calibration functions
  #   save_path: Optional path to save the plot
  #
  # Returns:
  #   ggplot object: Calibration summary plot
  
  # Extract key information
  if (is.list(calibration_results) && "calibration_results" %in% names(calibration_results)) {
    data <- calibration_results$calibration_results
    optimal_value <- calibration_results$optimal_c_poc
    target_rate <- calibration_results$target_rate
    achieved_rate <- calibration_results$optimal_rate
  } else {
    stop("Invalid calibration results format")
  }
  
  # Create summary data
  summary_data <- data.frame(
    Parameter = c("Optimal C_poc", "Target Rate", "Achieved Rate", "Difference"),
    Value = c(round(optimal_value, 3), 
              round(target_rate, 3), 
              round(achieved_rate, 3), 
              round(abs(achieved_rate - target_rate), 3)),
    Type = c("Threshold", "Target", "Achieved", "Error")
  )
  
  # Create bar plot
  p <- ggplot(summary_data, aes(x = .data$Parameter, y = .data$Value, fill = .data$Type)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_text(aes(label = Value), vjust = -0.3, size = 4, fontface = "bold") +
    labs(
      title = "Calibration Summary",
      subtitle = paste("PoC Calibration Results (", unique(data$n_simulations), "simulations per point)"),
      x = "Parameter",
      y = "Value",
      fill = "Type"
    ) +
    scale_fill_manual(values = c("Threshold" = "#009E73", "Target" = "#E69F00", 
                                "Achieved" = "#56B4E9", "Error" = "#CC79A7")) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  if (!is.null(save_path)) {
    # Create directory if it doesn't exist
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300)
    cat("Calibration summary saved to:", save_path, "\n")
  }
  
  return(p)
}

create_calibration_report <- function(calibration_results, output_dir = "results/calibration") {
  # Create a comprehensive calibration report with multiple plots
  #
  # Args:
  #   calibration_results: Results from calibration functions
  #   output_dir: Directory to save all plots
  #
  # Returns:
  #   list: All generated plots
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  plots <- list()
  
  # 1. Main calibration curve
  plots$calibration_curve <- plot_poc_calibration_curve(
    calibration_results, 
    save_path = file.path(output_dir, "poc_calibration_curve.png")
  )
  
  # 2. Calibration summary
  plots$summary <- plot_calibration_summary(
    calibration_results,
    save_path = file.path(output_dir, "calibration_summary.png")
  )
  
  # 3. Performance metrics table (if available)
  if ("validation_results" %in% names(calibration_results)) {
    # Create validation results table
    validation_data <- calibration_results$validation_results
    validation_table <- data.frame(
      Metric = c("C_poc", "Validation Rate", "Target Rate", "Difference", "95% CI Lower", "95% CI Upper"),
      Value = c(
        round(validation_data$optimal_c_poc, 3),
        round(validation_data$validation_rate, 3),
        round(validation_data$target_rate, 3),
        round(abs(validation_data$validation_rate - validation_data$target_rate), 3),
        round(validation_data$validation_ci[1], 3),
        round(validation_data$validation_ci[2], 3)
      )
    )
    
    # Save validation table
    write.csv(validation_table, file.path(output_dir, "validation_results.csv"), row.names = FALSE)
    cat("Validation results table saved to:", file.path(output_dir, "validation_results.csv"), "\n")
  }
  
  # 4. Save calibration data
  save(calibration_results, file = file.path(output_dir, "calibration_results.RData"))
  cat("Calibration results saved to:", file.path(output_dir, "calibration_results.RData"), "\n")
  
  cat("Calibration report created in:", output_dir, "\n")
  cat("Generated plots:\n")
  cat("  - poc_calibration_curve.png\n")
  cat("  - calibration_summary.png\n")
  cat("  - validation_results.csv (if available)\n")
  cat("  - calibration_results.RData\n")
  
  return(plots)
}
