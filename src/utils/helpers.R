plot_posterior_summary <- function(posterior_df,
                                   dose_col = "d",
                                   mean_col = "pava_mean", # Default to pava_mean
                                   ci_lower_col = "pava_ci_lower",
                                   ci_upper_col = "pava_ci_upper",
                                   group_col = NULL,
                                   title = "Posterior Mean and 95% CI by Dose",
                                   file_path = NULL) {
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

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 6)
  }

  p +
    scale_x_continuous(breaks = unique(posterior_df[[dose_col]])) +
    labs(x = "Dose Level", y = "Posterior Mean (with 95% CI)", title = title) +
    theme_minimal(base_size = 14)
}