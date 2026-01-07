# 绘图功能演示脚本
# 展示如何使用新的图表功能

# 加载必要的库和函数
source("src/utils/helpers.R")
source("src/utils/plotting_extensions.R")

# 显示使用指南
print_plotting_guide()

# 运行演示
cat("\n=== 开始演示 ===\n")
demo_plots <- demo_new_plots()

# 显示创建的图表
cat("\n=== 创建的图表 ===\n")
cat("1. 多场景剂量-反应曲线图: results/plots/demo_multi_scenarios.png\n")
cat("2. OBD选择率对比图: results/plots/demo_obd_selection.png\n")
cat("3. 完整评估图表集: results/plots/demo_*.png\n")

# 展示如何在你的项目中使用
cat("\n=== 在你的项目中使用 ===\n")
cat("# 1. 创建你的试验结果数据\n")
cat("my_results <- list(\n")
cat("  scenarios = your_scenario_data,\n")
cat("  selection_rates = your_selection_data,\n")
cat("  safety_rates = your_safety_data,\n")
cat("  # ... 其他数据\n")
cat(")\n\n")

cat("# 2. 创建所有图表\n")
cat("all_plots <- create_comprehensive_evaluation_plots(\n")
cat("  my_results,\n")
cat("  file_prefix = 'my_analysis'\n")
cat(")\n\n")

cat("# 3. 或者创建单个图表\n")
cat("single_plot <- plot_method_comparison_bars(\n")
cat("  data = my_data,\n")
cat("  x_var = 'scenario',\n")
cat("  y_var = 'success_rate',\n")
cat("  fill_var = 'method',\n")
cat("  title = 'My Analysis'\n")
cat(")\n")

