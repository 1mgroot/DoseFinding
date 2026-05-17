# 快速开始指南

## 项目简介
这是一个贝叶斯剂量寻找临床试验仿真项目，实现了自适应随机化和基于效用的剂量选择。

## 快速开始 (5分钟)

### 1. 环境准备
```r
# 安装依赖包
install.packages(c("dplyr", "ggplot2", "isotone", "purrr", "tidyr", "Iso"))
```

### 2. 运行基本仿真
```r
# 加载核心函数
source("src/core/config.R")
source("src/core/main.R")

# 运行试验仿真
results <- run_trial_simulation(
  trial_config = trial_config,
  p_YI = p_YI,
  p_YT_given_I = p_YT_given_I,
  p_YE_given_I = p_YE_given_I,
  rho0 = rho0,
  rho1 = rho1
)

# 查看结果
cat("最终选择的剂量:", results$final_dose, "\n")
cat("是否早期终止:", results$terminated_early, "\n")
cat("PoC是否验证:", results$poc_validated, "\n")
```

### 3. 运行校准演示
```r
# 运行综合校准演示 (约5-10分钟)
source("examples/comprehensive_calibration_demo.R")
```

### 4. 使用交互式笔记本
```r
# 打开仿真笔记本
# 文件: notebooks/simulation_notebook.qmd
# 包含完整的交互式示例
```

## 主要功能

### ✅ 试验仿真
- 多阶段自适应设计
- 贝叶斯后验计算
- 自适应随机化
- 早期终止机制
- PoC验证

### ✅ 校准系统
- PoC校准 (目标: 10%检测率)
- 早期终止校准 (目标: 80%终止率)
- 性能可视化
- 参数优化

### ✅ 可视化
- 剂量-反应曲线
- 后验分布图
- 校准曲线
- 分配分析

## 示例脚本

| 脚本 | 功能 | 运行时间 |
|------|------|----------|
| `examples/plotting_demo.R` | 基本仿真和可视化 | 1分钟 |
| `examples/poc_calibration_demo.R` | PoC校准演示 | 5分钟 |
| `examples/comprehensive_calibration_demo.R` | 完整校准系统 | 10分钟 |
| `notebooks/simulation_notebook.qmd` | 交互式笔记本 | 可变 |

## 配置参数

### 基本配置
```r
trial_config <- list(
  dose_levels = c(1, 2, 3),      # 剂量水平
  n_stages = 3,                  # 试验阶段数
  cohort_size = 6,               # 每阶段队列大小
  phi_T = 0.3, c_T = 0.9,       # 毒性阈值
  phi_E = 0.2, c_E = 0.9,       # 疗效阈值
  phi_I = 0.2, c_I = 0.8,       # 免疫反应阈值
  c_poc = 0.9,                  # PoC阈值
  enable_early_termination = TRUE  # 启用早期终止
)
```

## 常见问题

### Q: 如何修改试验参数？
A: 编辑 `src/core/config.R` 中的 `trial_config` 参数。

### Q: 如何运行校准？
A: 使用 `source("examples/comprehensive_calibration_demo.R")` 运行完整校准演示。

### Q: 如何查看结果？
A: 结果包含在返回的列表中，使用 `print(results)` 查看所有结果。

### Q: 如何生成图表？
A: 使用 `source("src/utils/helpers.R")` 加载绘图函数，或运行示例脚本。

## 详细文档

- **完整项目概览**: `docs/PROJECT_OVERVIEW.md`
- **试验设计规范**: `docs/TRIAL_DESIGN.md`
- **校准实施总结**: `docs/CALIBRATION_IMPLEMENTATION_SUMMARY.md`

## 支持

如有问题，请查看：
1. 示例脚本 (`examples/` 目录)
2. 测试文件 (`tests/` 目录)
3. 交互式笔记本 (`notebooks/` 目录)
4. 详细文档 (`docs/` 目录)
