# Parameter Naming Convention

## 问题
之前使用 `phi_I`, `phi_E`, `phi_T` 表示两个不同的概念，容易混淆：
- 零假设场景的**真实概率值**（用于生成数据）
- Admissibility检查的**阈值**（用于决策）

## 解决方案：清晰的命名约定

### 1. `null_p_*` - 零假设场景真实值
用于**生成模拟数据**的真实概率：
- `null_p_I` = 真实免疫反应率 (例如: 0.25)
- `null_p_E` = 真实疗效率 (例如: 0.30)
- `null_p_T` = 真实毒性率 (例如: 0.05)

### 2. `phi_*` - Admissibility阈值
用于**决策制定**的阈值参数：
- `phi_I` = 免疫反应阈值 (例如: 0.20)
- `phi_E` = 疗效阈值 (例如: 0.25)
- `phi_T` = 毒性阈值 (例如: 0.30)

### 3. `c_*` - 概率阈值
用于判断的概率下限：
- `c_I` = P(Imm > phi_I) 必须超过此值 (例如: 0.3)
- `c_E` = P(Eff > phi_E) 必须超过此值 (例如: 0.3)
- `c_T` = P(Tox < phi_T) 必须超过此值 (例如: 0.3)

## 对比表格

| 参数类型 | 免疫反应 | 疗效 | 毒性 | 用途 |
|---------|---------|------|------|------|
| **零假设真实值** | null_p_I = 0.25 | null_p_E = 0.30 | null_p_T = 0.05 | 生成数据 |
| **Admissibility阈值** | phi_I = 0.20 | phi_E = 0.25 | phi_T = 0.30 | 决策判断 |
| **关系** | null_p_I > phi_I | null_p_E > phi_E | null_p_T < phi_T | 设计意图 |

## 为什么这样设计？

在零假设场景中，我们设置 `null_p_*` 略高于/低于 `phi_*` 阈值：

1. **免疫和疗效**: `null_p_I > phi_I`, `null_p_E > phi_E`
   - 留出margin补偿PAVA偏差
   - PAVA在flat场景中会产生向下偏差
   - 避免因统计偏差导致过多早期终止

2. **毒性**: `null_p_T < phi_T`
   - 零假设场景应该是安全的
   - 不希望因毒性而终止

## 代码示例

### 配置参数

```r
# 零假设场景参数（用于生成数据）
calibration_params <- list(
  null_p_I = 0.25,    # 真实免疫反应率
  null_p_E = 0.30,    # 真实疗效率
  null_p_T = 0.05     # 真实毒性率
)

# Trial配置（用于决策）
trial_config <- list(
  phi_I = 0.20,       # 免疫反应阈值
  c_I = 0.3,          # P(Imm > 0.20) 必须 > 0.3
  
  phi_E = 0.25,       # 疗效阈值
  c_E = 0.3,          # P(Eff > 0.25) 必须 > 0.3
  
  phi_T = 0.30,       # 毒性阈值
  c_T = 0.3           # P(Tox < 0.30) 必须 > 0.3
)
```

### 使用场景

```r
# 1. 创建零假设场景（使用 null_p_*）
null_scenario <- create_null_flat_scenario(
  phi_I = calibration_params$null_p_I,  # 注意：这里传入的是null_p_I
  phi_E = calibration_params$null_p_E,  # 虽然参数名还叫phi_I，但含义是真实值
  tox_flat = calibration_params$null_p_T
)

# 2. Admissibility检查（使用 phi_* 和 c_*）
for (dose in doses) {
  imm_prob <- mean(posterior_samples > trial_config$phi_I)  # 使用phi_I阈值
  is_admissible <- imm_prob >= trial_config$c_I            # 使用c_I概率阈值
}
```

## 影响的文件

### 已更新
- `notebooks/poc_calibration_notebook.qmd` - 使用新的命名约定

### 未更新（保持兼容）
- `src/optimization/poc_calibration_new.R` - `create_null_flat_scenario()` 函数
  - 参数名仍然是 `phi_I`, `phi_E`（为保持兼容性）
  - 但调用时传入 `null_p_I`, `null_p_E` 的值
  - 注释清楚说明这些是"真实值"而非"阈值"

## 关键原则

✅ **DO**:
- 用 `null_p_*` 表示零假设场景的真实概率
- 用 `phi_*` 表示admissibility检查的阈值
- 在注释中明确说明参数用途

❌ **DON'T**:
- 混用同一个变量名表示不同含义
- 在代码中使用模糊的参数名
- 假设读者能从上下文推断参数含义

## 理解检查

**问题**: 如果我想测试一个"真正的零假设"场景（所有剂量表现相同，在阈值边界上），应该设置什么参数？

**错误答案**: 
```r
null_p_I = 0.20  # 错误！会因PAVA偏差导致100%早期终止
```

**正确答案**: 
```r
null_p_I = 0.25  # 正确！留出margin补偿PAVA偏差
# 这样后验估计在PAVA偏差后仍能超过phi_I = 0.20
```

---

更新日期: 2025-10-29
