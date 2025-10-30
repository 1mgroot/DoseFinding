# PoC Calibration System Updates

## 更新时间
2025-10-29

## 主要更新

### 0. 命名约定改进 ✅
**问题**: `phi_I`, `phi_E`, `phi_T` 在两个地方使用，含义不同，容易混淆
**解决方案**: 
- 引入 `null_p_*` 命名约定表示零假设场景的真实值
  - `null_p_I` = 真实免疫反应率（用于生成数据）
  - `null_p_E` = 真实疗效率（用于生成数据）
  - `null_p_T` = 真实毒性率（用于生成数据）
- 保留 `phi_*` 表示admissibility阈值
  - `phi_I` = 免疫反应阈值（用于决策）
  - `phi_E` = 疗效阈值（用于决策）
  - `phi_T` = 毒性阈值（用于决策）

**参考**: 详见 `NAMING_CONVENTION.md`

### 1. 种子随机化修复 ✅
**问题**: 所有模拟都使用相同的seed (123)，导致结果完全相同
**解决方案**: 
- 在 `run_trial_simulation()` 添加seed参数
- 每个模拟使用唯一种子: `base_seed = 10000 * c_poc_index + sim`
- 保持向后兼容（seed = NULL 使用当前RNG状态）

**修改文件**:
- `src/core/main.R`
- `src/core/simulate_data.R`
- `src/optimization/poc_calibration_new.R`
- `src/decision/dose_decision.R` (bug修复)

### 2. 零假设场景参数调整 ✅
**问题**: phi_I = 0.20 恰好在阈值上，PAVA偏差导致100%早期终止
**解决方案**: 
- null_p_I: 0.20 → 0.25 (高于阈值留出margin)
- null_p_E: 0.25 → 0.30 (高于阈值留出margin)
- 同时改进命名：使用 `null_p_*` 表示零假设真实值，与 `phi_*` (阈值) 区分

**理由**: 
在flat scenario中，PAVA会引入向下偏差，特别是样本量小时。设置真实值高于阈值可以补偿这个偏差。

### 3. 详细报告生成功能 ✅
**新功能**: `generate_calibration_report()` 函数
**输出位置**: `notebooks/results/notebook_calibration/calibration_detailed_report.txt`

**报告内容**:
1. **配置摘要**
   - Trial design参数
   - Admissibility阈值
   - Null scenario真实参数

2. **校准结果摘要**
   - 所有C_poc值的结果表格
   - PoC检测率和早期终止率对比

3. **早期终止详细分析** (重点！)
   - 每个C_poc的统计信息
   - 按stage的终止分布
   - 终止时的样本量统计
   - **3个详细案例展示**:
     * 后验估计值（Toxicity, Efficacy, Immune）
     * 每个剂量的admissibility概率
     * ✓/✗ 标记哪些剂量通过检查
   - 完成试验的分析

4. **建议**
   - 推荐的C_poc值
   - 预期试验特征
   - 高早期终止率警告
   - 下一步行动建议

## 使用方法

### 运行校准（在notebook中）

```r
# 1. 运行校准
calibration_results <- calibrate_c_poc(
  null_scenario = null_scenario,
  c_poc_candidates = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
  n_simulations = 1000,
  base_config = trial_config
)

# 2. 生成详细报告
generate_calibration_report(
  calibration_results = calibration_results,
  null_scenario = null_scenario,
  base_config = trial_config,
  file_path = "notebooks/results/notebook_calibration/calibration_detailed_report.txt"
)
```

### 查看报告
报告以纯文本格式保存，易于阅读：
- 用任何文本编辑器打开
- 有清晰的章节分隔
- 包含表格和示例
- 中文注释友好

## 预期效果

### 修复前
- ❌ 所有模拟使用相同数据
- ❌ 100% 早期终止率
- ❌ Immune posterior: 0.15-0.18 (远低于0.20)
- ❌ P(Imm > 0.2) = 0.1-0.28 (全部失败)

### 修复后预期
- ✅ 每次模拟使用不同随机数据
- ✅ 早期终止率降低到合理范围 (期望 30-50%)
- ✅ Immune posterior 更接近真实值 0.25
- ✅ P(Imm > 0.2) 大部分能通过阈值 c_I = 0.3
- ✅ 能够完成足够多的试验进行校准

## 技术细节

### PAVA偏差问题
**PAVA (Pool Adjacent Violators Algorithm)** 强制单调性：
- 对于真正单调递增的数据有效（如toxicity, efficacy vs dose）
- 对于flat数据会引入系统性偏差
- 小样本时偏差更明显

**解决策略**:
1. **保留PAVA** (用户要求)
2. **调整零假设参数** 高于阈值以补偿偏差
3. **增大样本量** (如果需要，可以增加cohort_size)

### 种子策略
```
Simulation seed = 10000 * c_poc_index + simulation_number
Stage seed = simulation_seed + stage_number
```

这确保：
- 不同c_poc值的模拟独立
- 同一试验的不同stage有不同随机性
- 相同seed可重现结果

## 下一步

1. **运行更新后的notebook**
   ```bash
   # 在RStudio中打开并运行
   notebooks/poc_calibration_notebook.qmd
   ```

2. **检查报告**
   ```bash
   # 查看生成的详细报告
   cat notebooks/results/notebook_calibration/calibration_detailed_report.txt
   ```

3. **分析结果**
   - 早期终止率是否降到合理范围？
   - Immune posterior估计是否改善？
   - 是否有足够的完成试验进行校准？

4. **如果早期终止率还是太高**
   - 进一步增大 phi_I (0.25 → 0.28 或 0.30)
   - 降低 c_I (0.3 → 0.25 或 0.2)
   - 增大 cohort_size (15 → 20)

## 文件清单

### 修改的文件
- `src/core/main.R` - seed参数
- `src/core/simulate_data.R` - seed处理
- `src/optimization/poc_calibration_new.R` - seed传递 + 报告生成
- `src/decision/dose_decision.R` - NULL检查
- `notebooks/poc_calibration_notebook.qmd` - 参数调整 + 报告调用

### 新增功能
- `generate_calibration_report()` - 详细报告生成函数

### 输出文件
- `notebooks/results/notebook_calibration/calibration_detailed_report.txt` - 主要报告

## 问题排查

如果遇到问题：

1. **Seed不工作**: 检查 `simulate_data_gumbel` 是否收到seed参数
2. **报告生成失败**: 确保输出目录存在 `notebooks/results/notebook_calibration/`
3. **还是100%早期终止**: 进一步增大 phi_I 或降低 c_I
4. **报告太长**: 调整 `max_examples` 参数减少案例数量

