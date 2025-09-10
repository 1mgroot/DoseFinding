# 扩展的绘图函数 - 基于reference code的图表类型
# 这些函数提供了reference code中展示的所有图表类型

library(ggplot2)
library(dplyr)

# 示例函数：创建模拟数据用于展示新图表
create_example_evaluation_data <- function() {
  
  # 创建示例场景数据
  scenarios <- list(
    list(
      toxicity = c(0.1, 0.18, 0.35, 0.40, 0.50),
      efficacy = c(0.35, 0.35, 0.37, 0.39, 0.39),
      utility = c(0.27, 0.23, 0.10, 0.13, 0.17)
    ),
    list(
      toxicity = c(0.05, 0.15, 0.25, 0.35, 0.50),
      efficacy = c(0.10, 0.35, 0.35, 0.38, 0.39),
      utility = c(0.07, 0.22, 0.22, 0.12, 0.06)
    ),
    list(
      toxicity = c(0.02, 0.06, 0.10, 0.20, 0.35),
      efficacy = c(0.05, 0.10, 0.35, 0.35, 0.40),
      utility = c(0.03, 0.07, 0.28, 0.22, 0.13)
    )
  )
  
  # 创建示例方法对比数据
  methods <- c("Current", "Proposed", "Reference")
  scenarios_names <- c("Scenario 1", "Scenario 2", "Scenario 3")
  
  # 剂量选择率数据
  selection_rates <- expand.grid(
    scenario = scenarios_names,
    method = methods,
    stringsAsFactors = FALSE
  )
  selection_rates$obd_rate <- c(45, 60, 55, 70, 65, 50, 80, 85, 75)
  
  # 安全性数据
  safety_rates <- expand.grid(
    scenario = scenarios_names,
    method = methods,
    stringsAsFactors = FALSE
  )
  safety_rates$mtd_rate <- c(30, 25, 35, 20, 15, 30, 15, 10, 20)
  
  # 样本量数据
  sample_sizes <- expand.grid(
    scenario = scenarios_names,
    method = methods,
    stringsAsFactors = FALSE
  )
  sample_sizes$avg_n <- c(25, 20, 30, 22, 18, 28, 18, 15, 25)
  
  # 过量患者数据
  overdose_rates <- expand.grid(
    scenario = scenarios_names,
    method = methods,
    stringsAsFactors = FALSE
  )
  overdose_rates$overdose_pct <- c(15, 10, 20, 12, 8, 18, 8, 5, 15)
  
  # 持续时间数据
  durations <- expand.grid(
    scenario = scenarios_names,
    method = methods,
    stringsAsFactors = FALSE
  )
  durations$duration <- c(12, 10, 15, 11, 9, 14, 9, 7, 12)
  
  # 效率数据
  efficiency <- expand.grid(
    scenario = scenarios_names,
    method = methods,
    stringsAsFactors = FALSE
  )
  efficiency$efficiency <- c(2.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.1, 2.1)
  
  return(list(
    scenarios = scenarios,
    selection_rates = selection_rates,
    safety_rates = safety_rates,
    sample_sizes = sample_sizes,
    overdose_rates = overdose_rates,
    durations = durations,
    efficiency = efficiency
  ))
}

# 快速演示函数
demo_new_plots <- function() {
  cat("=== 演示新的图表功能 ===\n")
  
  # 创建示例数据
  example_data <- create_example_evaluation_data()
  
  # 创建多场景剂量-反应曲线图
  cat("1. 创建多场景剂量-反应曲线图...\n")
  multi_scenario_plot <- plot_multi_scenario_curves(
    example_data$scenarios,
    title = "Dose-Response Curves Across Scenarios",
    file_path = "results/plots/demo_multi_scenarios.png"
  )
  
  # 创建方法对比图
  cat("2. 创建方法对比图...\n")
  obd_plot <- plot_method_comparison_bars(
    example_data$selection_rates,
    x_var = "scenario", y_var = "obd_rate", fill_var = "method",
    title = "OBD Selection Rate Comparison",
    y_label = "OBD Selection (%)",
    limits = c(0, 100),
    file_path = "results/plots/demo_obd_selection.png"
  )
  
  # 创建完整的评估图表集
  cat("3. 创建完整的评估图表集...\n")
  all_plots <- create_comprehensive_evaluation_plots(
    example_data,
    file_prefix = "demo"
  )
  
  cat("✅ 所有图表已创建完成！\n")
  cat("📁 图表保存在 results/plots/ 目录中\n")
  
  return(all_plots)
}

# 使用说明
print_plotting_guide <- function() {
  cat("=== Reference Code 图表类型分析 ===\n\n")
  
  cat("📊 Reference Code 中的图表类型：\n")
  cat("1. 多场景剂量-反应曲线图 (3×3网格布局)\n")
  cat("   - 展示9个不同场景的毒性、疗效和效用曲线\n")
  cat("   - 适合：参数敏感性分析、场景对比\n\n")
  
  cat("2. OBD选择率对比图\n")
  cat("   - 比较不同方法的最优生物剂量选择率\n")
  cat("   - 适合：方法性能评估、剂量选择准确性\n\n")
  
  cat("3. MTD选择率对比图\n")
  cat("   - 比较不同方法的最大耐受剂量选择率\n")
  cat("   - 适合：安全性评估、毒性控制\n\n")
  
  cat("4. 平均样本量对比图\n")
  cat("   - 比较不同方法的平均试验样本量\n")
  cat("   - 适合：试验效率评估、资源需求分析\n\n")
  
  cat("5. 过量患者百分比图\n")
  cat("   - 显示各方法中分配到过量剂量的患者比例\n")
  cat("   - 适合：安全性评估、患者保护\n\n")
  
  cat("6. 试验持续时间图\n")
  cat("   - 比较不同方法的试验持续时间\n")
  cat("   - 适合：时间效率评估、试验规划\n\n")
  
  cat("7. 入组效率图\n")
  cat("   - 计算入组效率 (样本量/持续时间)\n")
  cat("   - 适合：综合效率评估、试验优化\n\n")
  
  cat("🎯 推荐使用场景：\n")
  cat("- 参数优化结果展示\n")
  cat("- 方法对比分析\n")
  cat("- 临床试验报告\n")
  cat("- 学术论文图表\n")
  cat("- 监管提交材料\n\n")
  
  cat("💡 使用方法：\n")
  cat("1. 运行 demo_new_plots() 查看示例\n")
  cat("2. 使用 create_comprehensive_evaluation_plots() 创建完整图表集\n")
  cat("3. 根据需要调整数据和样式\n")
}
