# 贝叶斯自适应剂量寻找试验仿真

一个全面的R实现，用于贝叶斯自适应剂量寻找试验，具有多阶段设计、基于效用的决策制定和参数优化功能。

## 🚀 快速开始

### 方法1: 交互式笔记本 (推荐)
```r
# 打开 notebooks/simulation_notebook.qmd
# 包含完整的交互式示例和校准框架
```

### 方法2: 直接脚本执行
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
```

### 方法3: 校准演示
```r
# 运行综合校准演示
source("examples/comprehensive_calibration_demo.R")
```

## 📁 项目结构

```
DoseFinding/
├── src/                              # 源代码
│   ├── core/                         # 核心仿真逻辑
│   │   ├── main.R                    # 主仿真脚本
│   │   ├── config.R                  # 试验配置参数
│   │   ├── simulate_data.R           # 数据仿真函数
│   │   └── model_utils.R             # 贝叶斯模型工具
│   ├── decision/                     # 决策逻辑
│   │   └── dose_decision.R           # 剂量选择算法
│   ├── optimization/                 # 参数优化
│   │   ├── poc_calibration.R         # PoC校准框架
│   │   └── early_termination_calibration.R  # 早期终止校准
│   └── utils/                        # 工具函数
│       ├── helpers.R                 # 辅助函数和绘图
│       ├── plotting_extensions.R     # 绘图扩展
│       └── calibration_plots.R       # 校准可视化
├── examples/                         # 示例脚本
│   ├── comprehensive_calibration_demo.R  # 综合校准演示
│   ├── poc_calibration_demo.R        # PoC校准演示
│   ├── flat_scenario_demo.R          # 平坦场景演示
│   └── bayesian_poc_demo.R           # 贝叶斯PoC演示
├── tests/                            # 测试文件
│   ├── test_comprehensive_calibration.R  # 综合校准测试
│   ├── test_poc_calibration.R        # PoC校准测试
│   └── test_*.R                      # 其他测试文件
├── docs/                             # 文档
│   ├── PROJECT_OVERVIEW.md           # 项目概览
│   ├── TRIAL_DESIGN.md               # 试验设计规范
│   ├── NEXT_STEP_PLAN.md             # 实施计划
│   └── CALIBRATION_IMPLEMENTATION_SUMMARY.md  # 校准实施总结
├── notebooks/                        # 交互式笔记本
│   └── simulation_notebook.qmd       # 交互式仿真笔记本
└── results/                          # 生成输出
    └── plots/                        # 生成的图表
```

## ✨ 主要功能

### 🎯 试验仿真
- **多阶段贝叶斯自适应设计** 与中期分析
- **基于效用的剂量选择** 与可定制效用函数
- **早期终止标准** 用于安全性和疗效
- **PoC验证** 正确选择概率验证
- **自适应随机化** 基于效用分数

### 🔧 校准系统
- **PoC校准** 目标: 10%检测率
- **早期终止校准** 目标: 80%终止率
- **性能可视化** 校准曲线和置信区间
- **参数优化** 系统化参数调优

### 📊 可视化
- **剂量-反应曲线** 毒性、疗效和效用
- **后验分布图** 现代样式
- **校准曲线** 阈值vs性能关系
- **分配分析** 参与者分布

## 📚 文档

### 快速开始
- **QUICK_START.md** - 5分钟快速开始指南
- **PROJECT_OVERVIEW.md** - 完整项目概览和使用方法

### 详细文档
- **TRIAL_DESIGN.md** - 完整试验设计规范
- **NEXT_STEP_PLAN.md** - 实施计划和状态
- **CALIBRATION_IMPLEMENTATION_SUMMARY.md** - 校准系统实施总结

## 🛠️ 系统要求

- R (>= 4.0)
- 必需包: dplyr, tidyr, isotone, purrr, ggplot2, Iso, testthat

## 📋 示例脚本

| 脚本 | 功能 | 运行时间 |
|------|------|----------|
| `examples/plotting_demo.R` | 基本仿真和可视化 | 1分钟 |
| `examples/poc_calibration_demo.R` | PoC校准演示 | 5分钟 |
| `examples/comprehensive_calibration_demo.R` | 完整校准系统 | 10分钟 |
| `notebooks/simulation_notebook.qmd` | 交互式笔记本 | 可变 |

## 🧪 测试

```r
# 运行所有测试
source("tests/test_comprehensive_calibration.R")

# 运行特定测试
source("tests/test_poc_calibration.R")
source("tests/test_early_termination_poc.R")
```

## 📈 项目状态

### ✅ 完成状态 (100%)
- ✅ 完整试验仿真工作流程
- ✅ 贝叶斯后验概率计算
- ✅ 自适应随机化算法
- ✅ 早期终止机制
- ✅ PoC验证系统
- ✅ 综合校准框架
- ✅ 性能可视化工具
- ✅ 完整集成测试

## 🤝 支持

如有问题，请查看：
1. **快速开始**: `QUICK_START.md`
2. **项目概览**: `docs/PROJECT_OVERVIEW.md`
3. **示例脚本**: `examples/` 目录
4. **交互式笔记本**: `notebooks/` 目录
5. **测试文件**: `tests/` 目录

## 📄 许可证

本项目用于研究和教育目的。
