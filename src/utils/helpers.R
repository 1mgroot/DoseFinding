# Enhanced plotting functions inspired by reference code
library(ggplot2)
library(dplyr)
plot_posterior_summary <- function(posterior_df,
                                   dose_col = "d",
                                   mean_col = "pava_mean", # Default to pava_mean
                                   ci_lower_col = "pava_ci_lower",
                                   ci_upper_col = "pava_ci_upper",
                                   group_col = NULL,
                                   title = "Posterior Mean and 95% CI by Dose",
                                   file_path = NULL,
                                   style = "modern") {
  
  if (style == "modern") {
    # Modern styling inspired by reference code
    p <- ggplot(posterior_df, aes(x = .data[[dose_col]], y = .data[[mean_col]]))

    if (!is.null(group_col) && group_col %in% colnames(posterior_df)) {
      p <- p +
        aes(color = factor(.data[[group_col]]), group = .data[[group_col]]) +
        geom_point(size = 3) +
        geom_line(linewidth = 1)

      if (ci_lower_col %in% colnames(posterior_df) && ci_upper_col %in% colnames(posterior_df)) {
        p <- p +
          geom_errorbar(aes(ymin = .data[[ci_lower_col]], ymax = .data[[ci_upper_col]]), 
                        width = 0.2, linewidth = 0.8)
      }

      p <- p + scale_color_manual(name = "", values = c("0" = "#E69F00", "1" = "#56B4E9"))
    } else {
      p <- p +
        geom_point(size = 3, color = "#009E73") +
        geom_line(color = "#009E73", linewidth = 1)

      if (ci_lower_col %in% colnames(posterior_df) && ci_upper_col %in% colnames(posterior_df)) {
        p <- p +
          geom_errorbar(aes(ymin = .data[[ci_lower_col]], ymax = .data[[ci_upper_col]]),
                        width = 0.2, color = "#009E73", linewidth = 0.8)
      }
    }

    p <- p +
      scale_x_continuous(breaks = unique(posterior_df[[dose_col]])) +
      labs(x = "Dose Level", y = "Probability", title = title) +
      theme_bw(base_size = 16) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.line = element_line(color = "black"),
            axis.line.y.right = element_line(color = "black"),
            axis.line.x.top = element_line(color = "black"))
    
  } else {
    # Original styling
    p <- ggplot(posterior_df, aes(x = .data[[dose_col]], y = .data[[mean_col]]))

    if (!is.null(group_col) && group_col %in% colnames(posterior_df)) {
      p <- p +
        aes(color = factor(.data[[group_col]]), group = .data[[group_col]]) +
        geom_point(size = 3) +
        geom_line()

      if (ci_lower_col %in% colnames(posterior_df) && ci_upper_col %in% colnames(posterior_df)) {
        p <- p +
          geom_errorbar(aes(ymin = .data[[ci_lower_col]], ymax = .data[[ci_upper_col]]), width = 0.2)
      }

      p <- p + scale_color_brewer(palette = "Set1", name = group_col)
    } else {
      p <- p +
        geom_point(size = 3, color = "blue") +
        geom_line(color = "blue")

      if (ci_lower_col %in% colnames(posterior_df) && ci_upper_col %in% colnames(posterior_df)) {
        p <- p +
          geom_errorbar(aes(ymin = .data[[ci_lower_col]], ymax = .data[[ci_upper_col]]),
                        width = 0.2, color = "blue")
      }
    }

    p <- p +
      scale_x_continuous(breaks = unique(posterior_df[[dose_col]])) +
      labs(x = "Dose Level", y = "Posterior Mean (with 95% CI)", title = title) +
      theme_minimal(base_size = 14)
  }

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 6, dpi = 400)
  }

  return(p)
}

# New function for creating dose-response curves similar to reference code
plot_dose_response_curves <- function(toxicity_data, efficacy_data, utility_data = NULL,
                                      title = "Dose-Response Curves", file_path = NULL) {
  
  # Prepare data for plotting
  dose_levels <- 1:length(toxicity_data)
  
  # Normalize utility data to 0-1 scale if provided (utilities are typically 0-100)
  utility_normalized <- if (!is.null(utility_data)) {
    # Check if utilities are on 0-100 scale and normalize if needed
    if (max(utility_data, na.rm = TRUE) > 1) {
      utility_data / 100
    } else {
      utility_data
    }
  } else {
    rep(NA, length(toxicity_data))
  }
  
  df_plot <- data.frame(
    dose = rep(dose_levels, 3),
    value = c(toxicity_data, efficacy_data, utility_normalized),
    group = rep(c("Toxicity", "Efficacy", "Utility"), each = length(toxicity_data))
  )
  
  # Remove utility if not provided or if all values are NA
  if (is.null(utility_data) || all(is.na(utility_normalized))) {
    df_plot <- df_plot[df_plot$group != "Utility", ]
  }
  
  df_plot$group <- factor(df_plot$group, levels = c("Toxicity", "Efficacy", "Utility"))
  
  # Determine y-axis limits based on data range
  max_prob_value <- max(c(toxicity_data, efficacy_data), na.rm = TRUE)
  max_util_value <- if (!is.null(utility_data) && !all(is.na(utility_normalized))) {
    max(utility_normalized, na.rm = TRUE)
  } else {
    0
  }
  y_max <- max(0.5, max(max_prob_value, max_util_value) * 1.1)  # Add 10% padding
  y_breaks <- if (y_max <= 0.5) {
    seq(0, 0.5, by = 0.1)
  } else if (y_max <= 1) {
    seq(0, 1, by = 0.2)
  } else {
    seq(0, ceiling(y_max), by = 0.5)
  }
  
  p <- ggplot(df_plot, aes(x = dose, y = value, color = group, shape = group, linetype = group)) +
    geom_line(linewidth = 1, na.rm = TRUE) +
    geom_point(size = 3, na.rm = TRUE) +
    geom_hline(yintercept = 0.25, color = "gray", linetype = "dashed", linewidth = 1) +
    labs(x = "Dose Level", y = if (is.null(utility_data) || all(is.na(utility_normalized))) "Probability" else "Probability / Normalized Utility", title = title) +
    scale_color_manual(name = "", values = c("Toxicity" = "#E69F00", "Efficacy" = "#56B4E9", "Utility" = "#009E73")) +
    scale_linetype_manual(name = "", values = c("Toxicity" = "dashed", "Efficacy" = "dashed", "Utility" = "solid")) +
    scale_shape_manual(name = "", values = c("Toxicity" = 17, "Efficacy" = 16, "Utility" = 15)) +
    scale_y_continuous(limits = c(0, y_max), breaks = y_breaks) +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(color = "black"),
          axis.line.y.right = element_line(color = "black"),
          axis.line.x.top = element_line(color = "black"))
  
  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 6, dpi = 400)
  }
  
  return(p)
}

# Function for creating method comparison bar charts
plot_method_comparison <- function(data, x_var, y_var, fill_var, 
                                  title = "Method Comparison", 
                                  y_label = "Value",
                                  file_path = NULL) {
  
  p <- ggplot(data, aes_string(x = x_var, y = y_var, fill = fill_var)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             color = "black", width = 0.8) +
    geom_text(aes_string(label = paste0("round(", y_var, ", 1)")),
              position = position_dodge(width = 0.8),
              vjust = -0.3, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, y = y_label) +
    scale_fill_manual(name = "", 
                      values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                      labels = c("Method 1", "Method 2", "Method 3", "Proposed")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top",
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 6, dpi = 400)
  }
  
  return(p)
}

# Function to create comprehensive trial summary plots
create_trial_summary_plots <- function(trial_results, file_prefix = "trial_summary") {
  
  plots <- list()
  
  # 1. Dose-response curves for true probabilities
  if (!is.null(trial_results$true_probabilities)) {
    true_probs <- trial_results$true_probabilities
    
    # Extract toxicity, efficacy, and utility data
    toxicity_data <- true_probs$toxicity
    efficacy_data <- true_probs$efficacy
    utility_data <- true_probs$utility
    
    plots$dose_response <- plot_dose_response_curves(
      toxicity_data, efficacy_data, utility_data,
      title = "True Dose-Response Curves",
      file_path = paste0("results/plots/", file_prefix, "_dose_response.png")
    )
  }
  
  # 2. Posterior summaries with modern styling
  if (!is.null(trial_results$posterior_summaries)) {
    plots$immune_response <- plot_posterior_summary(
      trial_results$posterior_summaries$imm,
      title = "Immune Response vs Dose (PAVA Adjusted)",
      file_path = paste0("results/plots/", file_prefix, "_immune_response.png"),
      style = "modern"
    )
    
    plots$toxicity <- plot_posterior_summary(
      trial_results$posterior_summaries$tox,
      title = "Toxicity Rate by Dose and Immune Status",
      group_col = "Y_I",
      file_path = paste0("results/plots/", file_prefix, "_toxicity.png"),
      style = "modern"
    )
    
    plots$efficacy <- plot_posterior_summary(
      trial_results$posterior_summaries$eff,
      title = "Efficacy Rate by Dose and Immune Status",
      group_col = "Y_I",
      file_path = paste0("results/plots/", file_prefix, "_efficacy.png"),
      style = "modern"
    )
  }
  
  # 3. Allocation over time
  if (!is.null(trial_results$all_alloc_probs)) {
    p_alloc_time <- ggplot(trial_results$all_alloc_probs, 
                          aes(x = Stage, y = Prob, color = factor(Dose))) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      labs(title = "Allocation Probabilities Over Time", 
           x = "Stage", y = "Allocation Probability", 
           color = "Dose Level") +
      scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
      theme_bw(base_size = 16) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.line = element_line(color = "black"))
    
    plots$allocation_time <- p_alloc_time
    
    ggsave(paste0("results/plots/", file_prefix, "_allocation_time.png"), 
           p_alloc_time, width = 10, height = 6, dpi = 400)
  }
  
  # 4. Participant allocation by dose and stage
  if (!is.null(trial_results$all_data)) {
    allocation_summary <- trial_results$all_data %>%
      group_by(d, stage) %>%
      summarise(n_participants = n(), .groups = 'drop') %>%
      mutate(d = factor(d), 
             stage = factor(stage, levels = 1:5, labels = paste("Stage", 1:5)))
    
    p_alloc <- ggplot(allocation_summary, aes(x = d, y = n_participants, fill = stage)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      labs(title = "Participant Allocation by Dose Level and Stage", 
           x = "Dose Level", y = "Number of Participants",
           subtitle = paste("Total participants:", sum(allocation_summary$n_participants))) +
      scale_fill_manual(name = "Stage", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
      theme_bw(base_size = 16) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            axis.line = element_line(color = "black"))
    
    plots$allocation_summary <- p_alloc
    
    ggsave(paste0("results/plots/", file_prefix, "_allocation_summary.png"), 
           p_alloc, width = 10, height = 6, dpi = 400)
  }
  
  return(plots)
}

# 多场景剂量-反应曲线图 (参考reference code)
plot_multi_scenario_curves <- function(scenarios_data, 
                                      title = "Dose-Response Curves Across Scenarios",
                                      file_path = NULL) {
  
  # scenarios_data should be a list with elements containing toxicity, efficacy, utility vectors
  n_scenarios <- length(scenarios_data)
  
  plot_list <- list()
  
  for (i in 1:n_scenarios) {
    scenario <- scenarios_data[[i]]
    
    # Prepare data for plotting
    dose_levels <- 1:length(scenario$toxicity)
    df_plot <- data.frame(
      dose = rep(dose_levels, 3),
      value = c(scenario$toxicity, scenario$efficacy, scenario$utility),
      group = rep(c("Toxicity", "Efficacy", "Utility"), each = length(dose_levels))
    )
    df_plot$group <- factor(df_plot$group, levels = c("Toxicity", "Efficacy", "Utility"))
    
    p <- ggplot(df_plot, aes(x = dose, y = value, color = group, shape = group, linetype = group)) +
      scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0.25, color = "gray", linetype = "dashed", linewidth = 1) +
      labs(x = "Dose Level", y = "Probability", title = paste("Scenario", i)) +
      scale_color_manual(name = "", values = c("Toxicity" = "#E69F00", "Efficacy" = "#56B4E9", "Utility" = "#009E73")) +
      scale_linetype_manual(name = "", values = c("Toxicity" = "dashed", "Efficacy" = "dashed", "Utility" = "solid")) +
      scale_shape_manual(name = "", values = c("Toxicity" = 17, "Efficacy" = 16, "Utility" = 15)) +
      theme_bw(base_size = 16) +
      theme(panel.grid = element_blank(),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            axis.line = element_line(color = "black"),
            axis.line.y.right = element_line(color = "black"),
            axis.line.x.top = element_line(color = "black"))
    
    # Add legend to first plot
    if (i == 1) {
      p <- p + theme(legend.position = "left")
    }
    
    # Add right y-axis label to middle plot
    if (i == ceiling(n_scenarios/2)) {
      p <- p + scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1), 
                                 position = "right", name = "Utility")
    }
    
    # Remove y-axis elements for non-edge plots
    if ((i - 1) %% 3 != 0) {
      p <- p + theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank())
    }
    
    plot_list[[i]] <- p
  }
  
  # Combine plots using patchwork
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = 3)
  } else {
    # Fallback: return list of plots
    combined_plot <- plot_list
  }
  
  if (!is.null(file_path)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      ggsave(file_path, combined_plot, width = 12, height = 3, dpi = 400)
    }
  }
  
  return(combined_plot)
}

# 方法对比柱状图 (参考reference code)
plot_method_comparison_bars <- function(data, x_var, y_var, fill_var,
                                      title = "Method Comparison",
                                      y_label = "Value (%)",
                                      limits = NULL,
                                      file_path = NULL) {
  
  p <- ggplot(data, aes_string(x = x_var, y = y_var, fill = fill_var)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             color = "black", width = 0.8) +
    geom_text(aes_string(label = paste0("round(", y_var, ", 1)")),
              position = position_dodge(width = 0.8),
              vjust = -0.3, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = limits) +
    labs(title = title, y = y_label) +
    scale_fill_manual(name = "", 
                      values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                      labels = c("Method 1", "Method 2", "Method 3", "Proposed")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top",
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 6.5, height = 2.5, dpi = 400)
  }
  
  return(p)
}

# 创建完整的试验评估图表集
create_comprehensive_evaluation_plots <- function(simulation_results, 
                                                 file_prefix = "evaluation") {
  
  plots <- list()
  
  # 1. 剂量选择率对比 (OBD Selection)
  if (!is.null(simulation_results$selection_rates)) {
    plots$obd_selection <- plot_method_comparison_bars(
      simulation_results$selection_rates,
      x_var = "scenario", y_var = "obd_rate", fill_var = "method",
      title = "OBD Selection Rate Comparison",
      y_label = "OBD Selection (%)",
      limits = c(0, 100),
      file_path = paste0("results/plots/", file_prefix, "_obd_selection.png")
    )
  }
  
  # 2. 安全性评估 (MTD Selection)
  if (!is.null(simulation_results$safety_rates)) {
    plots$mtd_selection <- plot_method_comparison_bars(
      simulation_results$safety_rates,
      x_var = "scenario", y_var = "mtd_rate", fill_var = "method",
      title = "MTD Selection Rate Comparison",
      y_label = "MTD Selection (%)",
      file_path = paste0("results/plots/", file_prefix, "_mtd_selection.png")
    )
  }
  
  # 3. 平均样本量对比
  if (!is.null(simulation_results$sample_sizes)) {
    plots$sample_sizes <- plot_method_comparison_bars(
      simulation_results$sample_sizes,
      x_var = "scenario", y_var = "avg_n", fill_var = "method",
      title = "Average Sample Size Comparison",
      y_label = "Average Sample Size",
      limits = c(0, 50),
      file_path = paste0("results/plots/", file_prefix, "_sample_sizes.png")
    )
  }
  
  # 4. 过量患者百分比
  if (!is.null(simulation_results$overdose_rates)) {
    plots$overdose_rates <- plot_method_comparison_bars(
      simulation_results$overdose_rates,
      x_var = "scenario", y_var = "overdose_pct", fill_var = "method",
      title = "Overdose Patient Percentage",
      y_label = "Overdose Pts (%)",
      limits = c(0, 25),
      file_path = paste0("results/plots/", file_prefix, "_overdose_rates.png")
    )
  }
  
  # 5. 试验持续时间
  if (!is.null(simulation_results$durations)) {
    plots$durations <- plot_method_comparison_bars(
      simulation_results$durations,
      x_var = "scenario", y_var = "duration", fill_var = "method",
      title = "Trial Duration Comparison",
      y_label = "Duration (months)",
      limits = c(0, 50),
      file_path = paste0("results/plots/", file_prefix, "_durations.png")
    )
  }
  
  # 6. 入组效率
  if (!is.null(simulation_results$efficiency)) {
    plots$efficiency <- plot_method_comparison_bars(
      simulation_results$efficiency,
      x_var = "scenario", y_var = "efficiency", fill_var = "method",
      title = "Enrollment Efficiency",
      y_label = "Relative Efficiency",
      file_path = paste0("results/plots/", file_prefix, "_efficiency.png")
    )
  }
  
  # 7. 多场景剂量-反应曲线
  if (!is.null(simulation_results$scenarios)) {
    plots$multi_scenarios <- plot_multi_scenario_curves(
      simulation_results$scenarios,
      title = "Dose-Response Curves Across Scenarios",
      file_path = paste0("results/plots/", file_prefix, "_multi_scenarios.png")
    )
  }
  
  return(plots)
}