# 贝叶斯剂量寻找试验仿真项目指南

## 📋 项目简介

这是一个基于贝叶斯方法的剂量寻找临床试验仿真项目，实现了多阶段自适应设计、效用驱动的剂量选择和参数校准系统。项目模拟包含毒性、疗效和免疫反应终点的临床试验。

### 🎯 核心功能
- **多阶段贝叶斯自适应设计** - 支持中期分析和自适应随机化
- **效用驱动的剂量选择** - 基于风险-效益权衡的决策
- **早期终止机制** - 安全性监控和试验效率优化
- **PoC验证系统** - 正确选择概率验证
- **参数校准框架** - PoC和早期终止参数优化

### 🎯 适用场景
- 临床试验设计和仿真
- 剂量寻找算法研究
- 自适应试验方法开发
- 统计方法验证

## 🚀 快速开始

### 环境配置
```r
# 安装依赖包
install.packages(c("dplyr", "ggplot2", "isotone", "purrr", "tidyr", "Iso"))
```

### 基本使用
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
cat("最终剂量:", results$final_od, "\n")
cat("是否早期终止:", results$terminated_early, "\n")
cat("PoC验证:", results$poc_validated, "\n")
```

## 🔧 核心功能说明

### 1. 多阶段工作流程
遵循 TRIAL_DESIGN.md 规范的五步工作流程：

1. **阶段1**: 所有剂量水平的等概率随机化
2. **中期分析**: 基于后验概率更新可接受剂量集
3. **早期终止检查**: 如果可接受剂量集为空则终止试验
4. **自适应随机化**: 基于效用分数分配患者（仅在试验继续时）
5. **最终选择**: 从可接受剂量集中选择最高效用的剂量 + PoC验证

### 2. 贝叶斯框架
- **后验分布**: 使用Beta先验计算
- **等张回归**: 应用PAVA确保单调性
- **边际概率**: 从联合分布计算

### 3. 效用计算
- **效用表**: 定义所有结果组合的值
- **期望效用**: 使用后验概率计算
- **总效用**: 结合免疫反应和非免疫反应场景

### 4. 校准系统
- **PoC校准**: 目标10%检测率
- **早期终止校准**: 目标80%终止率
- **性能可视化**: 校准曲线和置信区间

## 🏗️ 代码架构关系

### src/ 源代码模块

#### core/ - 核心仿真逻辑
- `config.R`: 所有配置参数、效用表和校准设置
- `main.R`: 主仿真函数 `run_trial_simulation()`，实现完整工作流程
- `simulate_data.R`: 使用Gumbel copula的数据生成函数
- `model_utils.R`: 贝叶斯模型工具和后验计算

#### decision/ - 决策算法
- `dose_decision.R`: 剂量决策逻辑、PoC计算、早期终止检查

#### optimization/ - 参数校准
- `poc_calibration.R`: PoC校准框架
- `early_termination_calibration.R`: 早期终止校准
- `run_optimization.R`: 校准运行脚本

#### utils/ - 工具函数
- `helpers.R`: 辅助函数和基础绘图
- `plotting_extensions.R`: 绘图扩展
- `calibration_plots.R`: 校准可视化

### notebooks/ 交互式笔记本

#### simulation_notebook.qmd - 主要仿真笔记本
- **与src代码的关系**: 直接调用 `src/core/main.R` 中的 `run_trial_simulation()` 函数
- **功能**: 提供交互式仿真环境，包含配置、仿真、结果可视化和校准演示
- **使用方式**: 
  ```r
  # 在notebook中直接运行
  results <- run_trial_simulation(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1)
  ```

### examples/ 示例脚本

| 脚本 | 功能 | 运行时间 | 调用关系 |
|------|------|----------|----------|
| `plotting_demo.R` | 基本仿真和可视化 | 1分钟 | 调用 `src/core/main.R` |
| `poc_calibration_demo.R` | PoC校准演示 | 5分钟 | 调用 `src/optimization/poc_calibration.R` |
| `comprehensive_calibration_demo.R` | 完整校准系统 | 10分钟 | 调用所有校准模块 |
| `flat_scenario_demo.R` | 平坦场景演示 | 3分钟 | 调用 `src/core/main.R` |
| `bayesian_poc_demo.R` | 贝叶斯PoC演示 | 2分钟 | 调用 `src/decision/dose_decision.R` |

### tests/ 测试文件
- `test_comprehensive_calibration.R`: 综合校准测试
- `test_poc_calibration.R`: PoC校准测试
- `test_early_termination_poc.R`: 早期终止和PoC测试
- `test_workflow_order.R`: 工作流程顺序验证

## 📊 使用场景

### 场景1: 运行单次试验仿真
```r
# 直接使用核心函数
source("src/core/config.R")
source("src/core/main.R")
results <- run_trial_simulation(trial_config, p_YI, p_YT_given_I, p_YE_given_I, rho0, rho1)
```

### 场景2: 使用交互式笔记本
```r
# 打开 notebooks/simulation_notebook.qmd
# 修改配置参数后运行代码块
# 获得交互式可视化和分析结果
```

### 场景3: 参数校准
```r
# 运行完整校准演示
source("examples/comprehensive_calibration_demo.R")

# 或单独运行PoC校准
source("src/optimization/poc_calibration.R")
poc_results <- run_quick_calibration(target_rate = 0.10, n_simulations = 100)
```

### 场景4: 批量仿真分析
```r
# 使用校准框架进行批量分析
source("src/optimization/run_optimization.R")
# 运行大规模仿真和参数优化
```

## ⚙️ 配置参数说明

### 试验参数配置
```r
trial_config <- list(
  dose_levels = c(1, 2, 3),      # 剂量水平
  n_stages = 3,                  # 试验阶段数
  cohort_size = 6,               # 每阶段队列大小
  phi_T = 0.3, c_T = 0.9,       # 毒性阈值
  phi_E = 0.2, c_E = 0.9,       # 疗效阈值
  phi_I = 0.2, c_I = 0.8,       # 免疫反应阈值
  c_poc = 0.9,                  # PoC阈值
  delta_poc = 0.8,              # PoC比较阈值
  enable_early_termination = TRUE,  # 启用早期终止
  log_early_termination = TRUE     # 启用详细日志
)
```

### 校准参数配置
```r
calibration_config <- list(
  poc_target_rate = 0.10,                    # PoC目标检测率
  poc_tolerance = 0.02,                      # PoC容差
  early_termination_target_rate = 0.80,      # 早期终止目标率
  early_termination_tolerance = 0.05,        # 早期终止容差
  n_calibration_simulations = 10000,         # 校准仿真次数
  n_validation_simulations = 5000            # 验证仿真次数
)
```

## 📈 结果解读

### 输出结果说明
```r
results <- run_trial_simulation(...)

# 主要结果字段
results$final_od              # 最终选择的剂量
results$final_utility         # 最终效用值
results$poc_validated        # PoC是否验证通过
results$poc_probability      # PoC概率值
results$terminated_early     # 是否早期终止
results$termination_stage    # 终止阶段
results$all_data             # 所有试验数据
results$posterior_summaries  # 后验摘要
results$all_alloc_probs      # 分配概率历史
```

### 可视化图表解读
- **剂量-反应曲线**: 显示毒性、疗效和效用的真实关系
- **后验摘要图**: 贝叶斯后验分布的可视化
- **分配分析**: 参与者分配的时间变化
- **校准曲线**: 参数阈值与性能的关系

### 关键指标含义
- **PoC概率**: 正确选择概率，>0.9表示高置信度
- **早期终止率**: 在不利场景中的终止比例
- **效用值**: 风险-效益综合评分
- **可接受剂量集**: 满足安全性、疗效和免疫反应标准的剂量

## ❓ 常见问题

### 安装问题
```r
# 确保所有依赖包已安装
install.packages(c("dplyr", "ggplot2", "isotone", "purrr", "tidyr", "Iso"))
```

### 运行问题
```r
# 确保工作目录正确
setwd("/path/to/DoseFinding")

# 检查路径设置
project_root <- "/Users/jz/Development/DoseFinding"
```

### 参数调整建议
- **增加仿真次数**: 提高校准精度，但增加计算时间
- **调整阈值**: 根据具体试验需求修改安全性和疗效标准
- **修改效用表**: 根据临床需求调整风险-效益权重

## 🔗 相关文档

- **详细设计规范**: `docs/TRIAL_DESIGN.md`
- **项目概览**: `docs/PROJECT_OVERVIEW.md`
- **快速开始**: `QUICK_START.md`
- **校准实施总结**: `docs/CALIBRATION_IMPLEMENTATION_SUMMARY.md`

## 📞 技术支持

如有问题，请查看：
1. 示例脚本 (`examples/` 目录)
2. 测试文件 (`tests/` 目录)
3. 交互式笔记本 (`notebooks/` 目录)
4. 详细文档 (`docs/` 目录)

---

**项目状态**: ✅ 100% 完成，所有功能已实现并经过测试验证



