library(dplyr)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

source("src/optimization/threshold_calibration.R")

audit_dir <- file.path(project_root, "_audit", "runtime", "ceti_n1_handcheck")
html_dir <- file.path(project_root, "_audit", "decision_review")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)

settings <- default_separate_threshold_settings(quick_mode = TRUE)
settings$n_sim_per_candidate <- 1
settings$show_progress <- FALSE
settings$output_dir <- audit_dir

started_at <- Sys.time()
calibration_results <- calibrate_separate_thresholds(settings)
finished_at <- Sys.time()
calibration_results$run_duration_seconds <- as.numeric(difftime(finished_at, started_at, units = "secs"))

summary_table <- threshold_calibration_summary_table(calibration_results)
candidate_table <- bind_rows(lapply(calibration_results$calibrations, function(result) {
  result$results %>%
    mutate(
      target_low = result$target_missing_range[[1]],
      target_high = result$target_missing_range[[2]],
      selection_status = result$status
    )
}))

scenario_table <- bind_rows(lapply(names(calibration_results$scenarios), function(param_name) {
  scenario <- calibration_results$scenarios[[param_name]]
  bind_rows(lapply(seq_along(settings$dose_levels), function(i) {
    p_i <- scenario$p_YI[[i]]
    p_t0 <- scenario$p_YT_given_I[i, 1]
    p_t1 <- scenario$p_YT_given_I[i, 2]
    p_e0 <- scenario$p_YE_given_I[i, 1]
    p_e1 <- scenario$p_YE_given_I[i, 2]
    data.frame(
      param_name = param_name,
      endpoint = scenario$endpoint,
      dose = settings$dose_levels[[i]],
      p_I = p_i,
      p_T_I0 = p_t0,
      p_T_I1 = p_t1,
      marginal_p_T_input = scenario$marginal_p_T[[i]],
      marginal_p_T_hand = (1 - p_i) * p_t0 + p_i * p_t1,
      p_E_I0 = p_e0,
      p_E_I1 = p_e1,
      marginal_p_E_input = scenario$marginal_p_E[[i]],
      marginal_p_E_hand = (1 - p_i) * p_e0 + p_i * p_e1,
      rho0 = scenario$rho0,
      rho1 = scenario$rho1,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  mutate(
    marginal_p_T_diff = marginal_p_T_hand - marginal_p_T_input,
    marginal_p_E_diff = marginal_p_E_hand - marginal_p_E_input
  )

gumbel_cells <- function(p_t, p_e, rho) {
  adj <- p_t * p_e * (1 - p_t) * (1 - p_e) * (exp(rho) - 1) / (exp(rho) + 1)
  c(
    `T=0,E=0` = (1 - p_t) * (1 - p_e) + adj,
    `T=0,E=1` = (1 - p_t) * p_e - adj,
    `T=1,E=0` = p_t * (1 - p_e) - adj,
    `T=1,E=1` = p_t * p_e + adj
  )
}

joint_check <- bind_rows(lapply(seq_len(nrow(scenario_table)), function(i) {
  row <- scenario_table[i, ]
  bind_rows(lapply(c(0, 1), function(immune_stratum) {
    p_t <- if (immune_stratum == 0) row$p_T_I0 else row$p_T_I1
    p_e <- if (immune_stratum == 0) row$p_E_I0 else row$p_E_I1
    cells <- gumbel_cells(p_t, p_e, 0)
    product_cells <- c(
      `T=0,E=0` = (1 - p_t) * (1 - p_e),
      `T=0,E=1` = (1 - p_t) * p_e,
      `T=1,E=0` = p_t * (1 - p_e),
      `T=1,E=1` = p_t * p_e
    )
    data.frame(
      param_name = row$param_name,
      endpoint = row$endpoint,
      dose = row$dose,
      immune_stratum = immune_stratum,
      p_T_given_I = p_t,
      p_E_given_I = p_e,
      cell = names(cells),
      rho0_cell_probability = as.numeric(cells),
      product_probability = as.numeric(product_cells),
      diff = as.numeric(cells - product_cells),
      stringsAsFactors = FALSE
    )
  }))
}))

joint_summary <- joint_check %>%
  group_by(param_name, endpoint, dose, immune_stratum) %>%
  summarise(
    max_abs_diff = max(abs(diff)),
    cell_sum = sum(rho0_cell_probability),
    .groups = "drop"
  )

target_range_distance_local <- function(values, target_range) {
  ifelse(
    values < target_range[[1]],
    target_range[[1]] - values,
    ifelse(values > target_range[[2]], values - target_range[[2]], 0)
  )
}

choose_expected_selected <- function(df, target_range, stricter_direction = "higher") {
  distance <- target_range_distance_local(df$final_admissible_missing_rate, target_range)
  in_range <- which(distance == 0)

  choose_least_strict <- function(indices) {
    values <- df$param_value[indices]
    if (stricter_direction == "higher") {
      indices[[which.min(values)]]
    } else {
      indices[[which.max(values)]]
    }
  }
  choose_most_strict <- function(indices) {
    values <- df$param_value[indices]
    if (stricter_direction == "higher") {
      indices[[which.max(values)]]
    } else {
      indices[[which.min(values)]]
    }
  }

  if (length(in_range) > 0) {
    selected_index <- choose_least_strict(in_range)
    status <- "within target range; least strict candidate"
  } else if (all(df$final_admissible_missing_rate < target_range[[1]])) {
    closest <- which(distance == min(distance))
    selected_index <- choose_most_strict(closest)
    status <- "below target range; strictest closest candidate"
  } else if (all(df$final_admissible_missing_rate > target_range[[2]])) {
    closest <- which(distance == min(distance))
    selected_index <- choose_least_strict(closest)
    status <- "above target range; least strict closest candidate"
  } else {
    closest <- which(distance == min(distance))
    selected_index <- choose_least_strict(closest)
    status <- "closest to target range; least strict tie-break"
  }

  list(selected_index = selected_index, status = status, distance = distance)
}

manual_selection_table <- bind_rows(lapply(split(candidate_table, candidate_table$param_name), function(df) {
  df <- df[order(df$param_value), , drop = FALSE]
  target_range <- c(df$target_low[[1]], df$target_high[[1]])
  expected <- choose_expected_selected(df, target_range)
  df$hand_missing_count <- df$final_admissible_missing_rate * df$n_simulations
  df$hand_target_endpoint_missing_count <- df$target_endpoint_missing_rate * df$n_simulations
  df$hand_distance_to_target <- expected$distance
  df$hand_selected <- seq_len(nrow(df)) == expected$selected_index
  df$hand_status <- expected$status
  df$selection_matches_hand <- df$selected == df$hand_selected
  df
})) %>%
  select(
    param_name,
    endpoint,
    param_value,
    c_T,
    c_I,
    c_E,
    n_simulations,
    hand_missing_count,
    final_admissible_missing_rate,
    hand_target_endpoint_missing_count,
    target_endpoint_missing_rate,
    target_blocks_overlap_rate,
    non_target_pair_empty_rate,
    early_stop_rate,
    target_low,
    target_high,
    hand_distance_to_target,
    selected,
    hand_selected,
    selection_matches_hand,
    selection_status,
    hand_status
  )

selected_comparison <- manual_selection_table %>%
  filter(selected | hand_selected) %>%
  select(
    param_name,
    endpoint,
    param_value,
    final_admissible_missing_rate,
    target_low,
    target_high,
    hand_distance_to_target,
    selected,
    hand_selected,
    selection_matches_hand,
    selection_status,
    hand_status
  )

scenario_expectation <- scenario_table %>%
  group_by(param_name, endpoint) %>%
  summarise(
    p_I_range = paste(range(p_I), collapse = " to "),
    marginal_p_T_range = paste(range(marginal_p_T_input), collapse = " to "),
    marginal_p_E_range = paste(range(marginal_p_E_input), collapse = " to "),
    max_abs_marginal_T_diff = max(abs(marginal_p_T_diff)),
    max_abs_marginal_E_diff = max(abs(marginal_p_E_diff)),
    design_expectation = dplyr::case_when(
      first(endpoint) == "toxicity" ~ "toxicity scenario: marginal P(T) is at/above phi_T while immune and efficacy are favorable",
      first(endpoint) == "immune" ~ "immune scenario: P(I) is below phi_I while toxicity and efficacy are favorable",
      first(endpoint) == "efficacy" ~ "efficacy scenario: marginal P(E) is below phi_E while toxicity and immune are favorable",
      TRUE ~ ""
    ),
    .groups = "drop"
  )

write.csv(summary_table, file.path(audit_dir, "summary_table.csv"), row.names = FALSE)
write.csv(candidate_table, file.path(audit_dir, "candidate_table_raw.csv"), row.names = FALSE)
write.csv(scenario_table, file.path(audit_dir, "scenario_probability_handcheck.csv"), row.names = FALSE)
write.csv(joint_check, file.path(audit_dir, "joint_cell_handcheck.csv"), row.names = FALSE)
write.csv(manual_selection_table, file.path(audit_dir, "manual_selection_handcheck.csv"), row.names = FALSE)
saveRDS(calibration_results, file.path(audit_dir, "calibration_results.rds"))

escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

table_html <- function(df, digits = 4) {
  df_out <- df
  df_out[] <- lapply(df_out, function(col) {
    if (is.numeric(col)) round(col, digits) else col
  })
  header <- paste0("<tr>", paste(sprintf("<th>%s</th>", escape_html(names(df_out))), collapse = ""), "</tr>")
  rows <- apply(df_out, 1, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", escape_html(row)), collapse = ""), "</tr>")
  })
  paste0("<table>", header, paste(rows, collapse = ""), "</table>")
}

formula_block <- function(lines) {
  paste0("<pre><code>", escape_html(paste(lines, collapse = "\n")), "</code></pre>")
}

html_path <- file.path(html_dir, "ceti_n1_handcheck.html")
html <- c(
  "<!doctype html>",
  "<html lang=\"zh-CN\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
  "<title>CETI n=1 Calibration Hand Check</title>",
  "<style>",
  "body{margin:0;background:#fbfaf7;color:#1f2933;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI','PingFang SC',sans-serif;line-height:1.58}main{max-width:1040px;margin:0 auto;padding:18px 14px 48px}section{background:#fff;border:1px solid #ddd6ca;border-radius:8px;padding:15px;margin:13px 0}h1{font-size:28px;line-height:1.12;margin:16px 0 8px}h2{font-size:20px;margin:0 0 9px}h3{font-size:16px;margin:12px 0 6px}.note{border-left:4px solid #0969da;background:#eef6ff;padding:10px 12px;border-radius:6px}.warn{border-left:4px solid #b54708;background:#fff7ed;padding:10px 12px;border-radius:6px}table{border-collapse:collapse;width:100%;font-size:12px;display:block;overflow-x:auto}th,td{border:1px solid #ddd;padding:5px 7px;text-align:left;white-space:nowrap}th{background:#f3f4f6}code{background:#f1f1f1;padding:1px 5px;border-radius:5px}pre{background:#111827;color:#f9fafb;padding:11px;border-radius:7px;overflow-x:auto;font-size:12px}.small{font-size:13px;color:#4b5563}",
  "</style></head><body><main>",
  "<h1>C_T / C_E / C_I n=1 Calibration Hand Check<br>从设计公式出发的小样本验算</h1>",
  "<p class=\"note\">This is an isolated audit artifact. It ran threshold calibration with <code>quick_mode=TRUE</code> and <code>n_sim_per_candidate=1</code>. It wrote only under <code>_audit/</code> and did not modify production outputs or design TeX files.</p>",
  "<section><h2>Run Settings / 运行设置</h2>",
  table_html(data.frame(
    setting = c("quick_mode", "n_sim_per_candidate", "calibration_seed", "rho0", "rho1", "target_missing_range", "c_T_candidates", "c_I_candidates", "c_E_candidates", "runtime_seconds"),
    value = c(
      settings$quick_mode,
      settings$n_sim_per_candidate,
      settings$calibration_seed,
      settings$rho0,
      settings$rho1,
      paste(settings$target_missing_range, collapse = ", "),
      paste(settings$c_T_candidates, collapse = ", "),
      paste(settings$c_I_candidates, collapse = ", "),
      paste(settings$c_E_candidates, collapse = ", "),
      round(calibration_results$run_duration_seconds, 3)
    )
  )),
  "</section>",
  "<section><h2>Design-side Probability Checks / 设计公式概率核对</h2>",
  "<p>For each endpoint-specific scenario, conditional probabilities are checked against the intended marginal probabilities.</p>",
  formula_block(c(
    "P(T | d) = (1 - P(I | d)) P(T | I=0,d) + P(I | d) P(T | I=1,d)",
    "P(E | d) = (1 - P(I | d)) P(E | I=0,d) + P(I | d) P(E | I=1,d)"
  )),
  table_html(scenario_expectation),
  "<h3>Dose-level marginal reconstruction</h3>",
  table_html(scenario_table %>% select(param_name, endpoint, dose, p_I, p_T_I0, p_T_I1, marginal_p_T_input, marginal_p_T_hand, marginal_p_T_diff, p_E_I0, p_E_I1, marginal_p_E_input, marginal_p_E_hand, marginal_p_E_diff)),
  "</section>",
  "<section><h2>rho=0 Joint-cell Check / 条件独立四格概率核对</h2>",
  "<p>With <code>rho=0</code>, the Gumbel adjustment term is zero, so every joint cell equals the product probability.</p>",
  formula_block(c(
    "adj = pT * pE * (1-pT) * (1-pE) * (exp(rho)-1)/(exp(rho)+1)",
    "rho = 0 => adj = 0",
    "P(T=t,E=e | I,d) = P(T=t | I,d) * P(E=e | I,d)"
  )),
  table_html(joint_summary),
  "</section>",
  "<section><h2>n=1 Selection Hand Calculation / n=1 候选选择手算</h2>",
  "<p>Because <code>n_simulations=1</code>, <code>final_admissible_missing_rate</code> can only be 0 or 1. With target range <code>[0.80, 0.85]</code>, rate 1 has distance 0.15 and rate 0 has distance 0.80. The manual selection below reproduces the implementation rule: choose in-range least strict if possible; otherwise choose closest to the target range with the coded tie-break.</p>",
  formula_block(c(
    "missing_count = final_admissible_missing_rate * n_simulations",
    "distance(rate, [0.80,0.85]) = 0 if 0.80 <= rate <= 0.85",
    "distance = 0.80 - rate if rate < 0.80",
    "distance = rate - 0.85 if rate > 0.85",
    "With n=1: distance(1) = 0.15, distance(0) = 0.80"
  )),
  "<h3>Selected candidates: implementation vs hand calculation</h3>",
  table_html(selected_comparison),
  "<h3>All candidate details</h3>",
  table_html(manual_selection_table),
  "</section>",
  "<section><h2>Calibration Summary / quick smoke result</h2>",
  table_html(summary_table),
  "</section>",
  "<section><h2>Files Written / 生成文件</h2>",
  "<ul>",
  "<li><code>_audit/decision_review/ceti_n1_handcheck.html</code></li>",
  "<li><code>_audit/runtime/ceti_n1_handcheck/summary_table.csv</code></li>",
  "<li><code>_audit/runtime/ceti_n1_handcheck/candidate_table_raw.csv</code></li>",
  "<li><code>_audit/runtime/ceti_n1_handcheck/scenario_probability_handcheck.csv</code></li>",
  "<li><code>_audit/runtime/ceti_n1_handcheck/joint_cell_handcheck.csv</code></li>",
  "<li><code>_audit/runtime/ceti_n1_handcheck/manual_selection_handcheck.csv</code></li>",
  "<li><code>_audit/runtime/ceti_n1_handcheck/calibration_results.rds</code></li>",
  "</ul>",
  "</section>",
  "</main></body></html>"
)
writeLines(html, html_path)

cat("CETI_N1_HANDCHECK_DONE\n")
cat("html_path=", html_path, "\n", sep = "")
cat("audit_dir=", audit_dir, "\n", sep = "")
print(summary_table)
