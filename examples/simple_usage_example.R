# 简单使用示例
# 这个脚本演示了如何使用贝叶斯剂量寻找试验仿真系统

# 加载必要的库
library(ggplot2)
library(dplyr)

# 设置工作目录
if (basename(getwd()) == "examples") {
  setwd("..")
}

# 加载核心函数
cat("正在加载核心函数...\n")
source("src/core/config.R")
source("src/core/main.R")
source("src/utils/helpers.R")

# 显示配置信息
cat("试验配置:\n")
cat("  剂量水平:", paste(trial_config$dose_levels, collapse = ", "), "\n")
cat("  试验阶段数:", trial_config$n_stages, "\n")
cat("  每阶段队列大小:", trial_config$cohort_size, "\n")
cat("  毒性阈值:", trial_config$phi_T, "\n")
cat("  疗效阈值:", trial_config$phi_E, "\n")
cat("  免疫反应阈值:", trial_config$phi_I, "\n")
cat("  PoC阈值:", trial_config$c_poc, "\n")
cat("  启用早期终止:", trial_config$enable_early_termination, "\n\n")

# 运行试验仿真
cat("正在运行试验仿真...\n")
results <- run_trial_simulation(
  trial_config = trial_config,
  p_YI = p_YI,
  p_YT_given_I = p_YT_given_I,
  p_YE_given_I = p_YE_given_I,
  rho0 = rho0,
  rho1 = rho1
)

# 显示结果
cat("=== 试验结果 ===\n")
cat("最终选择的剂量:", results$final_dose, "\n")
cat("是否早期终止:", results$terminated_early, "\n")
cat("PoC是否验证:", results$poc_validated, "\n")
cat("总参与者数:", nrow(results$all_data), "\n")
cat("可接受剂量集:", paste(results$admissible_set, collapse = ", "), "\n")

# 显示分配情况
if (!is.null(results$all_data)) {
  allocation_summary <- results$all_data %>%
    group_by(d) %>%
    summarise(n_participants = n(), .groups = 'drop')
  
  cat("\n=== 参与者分配 ===\n")
  for (i in 1:nrow(allocation_summary)) {
    cat("剂量", allocation_summary$d[i], ": ", allocation_summary$n_participants[i], " 参与者\n")
  }
}

# 显示后验概率摘要
if (!is.null(results$posterior_summaries)) {
  cat("\n=== 后验概率摘要 ===\n")
  
  # 免疫反应
  if (!is.null(results$posterior_summaries$imm)) {
    cat("免疫反应后验均值:\n")
    for (i in 1:nrow(results$posterior_summaries$imm)) {
      cat("  剂量", results$posterior_summaries$imm$d[i], ": ", 
          round(results$posterior_summaries$imm$pava_mean[i], 3), "\n")
    }
  }
  
  # 毒性
  if (!is.null(results$posterior_summaries$tox)) {
    cat("毒性后验均值:\n")
    for (i in 1:nrow(results$posterior_summaries$tox)) {
      cat("  剂量", results$posterior_summaries$tox$d[i], 
          " (I=", results$posterior_summaries$tox$Y_I[i], "): ", 
          round(results$posterior_summaries$tox$mean[i], 3), "\n")
    }
  }
  
  # 疗效
  if (!is.null(results$posterior_summaries$eff)) {
    cat("疗效后验均值:\n")
    for (i in 1:nrow(results$posterior_summaries$eff)) {
      cat("  剂量", results$posterior_summaries$eff$d[i], 
          " (I=", results$posterior_summaries$eff$Y_I[i], "): ", 
          round(results$posterior_summaries$eff$mean[i], 3), "\n")
    }
  }
}

# 创建简单的可视化
cat("\n正在创建可视化...\n")

# 创建结果目录
dir.create("results/simple_example", showWarnings = FALSE, recursive = TRUE)

# 1. 参与者分配图
if (!is.null(results$all_data)) {
  allocation_plot <- ggplot(results$all_data, aes(x = factor(d))) +
    geom_bar(fill = "steelblue", alpha = 0.7) +
    labs(title = "参与者分配 by 剂量水平",
         x = "剂量水平",
         y = "参与者数量") +
    theme_minimal()
  
  ggsave("results/simple_example/allocation_plot.png", allocation_plot, 
         width = 8, height = 6, dpi = 300)
  cat("  参与者分配图已保存: results/simple_example/allocation_plot.png\n")
}

# 2. 后验概率图
if (!is.null(results$posterior_summaries$imm)) {
  posterior_plot <- plot_posterior_summary(
    results$posterior_summaries$imm,
    title = "免疫反应后验概率",
    file_path = "results/simple_example/posterior_immune_response.png"
  )
  cat("  后验概率图已保存: results/simple_example/posterior_immune_response.png\n")
}

# 3. 试验总结
cat("\n=== 试验总结 ===\n")
cat("试验成功完成!\n")
cat("最终选择的剂量:", results$final_dose, "\n")
if (results$terminated_early) {
  cat("试验因安全考虑而早期终止\n")
} else {
  cat("试验按计划完成\n")
}
if (results$poc_validated) {
  cat("PoC验证通过\n")
} else {
  cat("PoC验证未通过\n")
}

cat("\n所有结果已保存到 results/simple_example/ 目录\n")
cat("示例完成!\n")
