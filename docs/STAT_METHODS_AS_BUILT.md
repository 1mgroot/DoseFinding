# STAT METHODS AS BUILT (代码为准 / Code-first)

目标：从生物统计角度，基于当前代码描述已实现的模型、算法、公式和工作流（非设计意图）。所有关键点均给出代码证据：`path function Lstart-Lend`。

## 1. 总览 / Overview
本仓库实现了一个带免疫、毒性、疗效三终点的贝叶斯多阶段剂量探索模拟。流程：按阶段模拟数据（Gumbel copula），对各终点做后验估计并用单调约束（PAVA/BIVISO），按阈值筛选可行剂量集，基于期望效用自适应分配，最终结合PoC概率做最优剂量选择。  
证据：`src/core/main.R run_trial_simulation L15-L167`; `src/core/simulate_data.R simulate_data_gumbel L9-L61`; `src/core/model_utils.R apply_pava_on_samples L33-L48`, `apply_biviso_on_matrix L50-L99`; `src/decision/dose_decision.R get_admissible_set L101-L140`, `adaptive_randomization L142-L155`, `calculate_poc_probability L310-L343`, `select_final_od_with_poc L364-L425`.

## 2. 数据生成 / Data Generation
- Gumbel copula 生成 (四格结局：T/E/Imm)；按剂量和免疫状态给定边际毒性/疗效概率。  
  公式：对每剂量 j，I~Bern(p_YI[j])；条件分布用 Gumbel(ttox, teff, c) 计算四格概率，再用 rmultinom 抽取 (Y_T,Y_E)。  
  证据：`src/core/simulate_data.R simulate_data_gumbel L9-L61`, `Gumbel L1-L8`.
- 随机种子：如果提供 seed，则每阶段使用 seed+stage，阶段内先 set.seed，然后采样 rbinom/rmultinom。  
  证据：`main.R L54-L66`; `simulate_data_gumbel L19-L22`.

## 3. 后验与单调约束 / Posterior & Isotonic
- Beta 后验采样：r/n 汇总→alpha_post=r+1, beta_post=n-r+1→rbeta(n_sims)。  
  证据：`model_utils.R simulate_beta_posterior L16-L23`, `add_beta_variance L25-L31`.
- PAVA (单变量剂量序)：对免疫概率样本矩阵逐行做 Iso::pava，加权=1/var_post，输出 pava_mean/CI。  
  证据：`model_utils.R apply_pava_on_samples L33-L48`.
- BIVISO (剂量×免疫分组)：补全缺失格子，生成样本矩阵后对二维做 Iso::biviso，输出 pava_mean/CI。  
  证据：`model_utils.R apply_biviso_on_matrix L50-L99`.
- 边际毒性/疗效：混合免疫概率与条件概率样本，得到 marginal_prob 与样本。  
  证据：`model_utils.R compute_marginal_probability L101-L128`.

## 4. 决策逻辑 / Decision Logic
- 期望效用：对每剂量 i，用后验均值 π_T|I, π_E|I, π_I 计算两层期望：  
  Utility_I0 = Σ w * P(E|I=0)⊗P(T|I=0); Utility_I1 同理；Total = (1-π_I)*U_I0 + π_I*U_I1。  
  证据：`dose_decision.R get_expected_utility L1-L27`.
- 可行集 A：对每剂量计算 P(π_T < φ_T), P(π_E > φ_E), P(π_I > φ_I) ，均超过 c_T/c_E/c_I 才入集。  
  证据：`dose_decision.R get_admissible_set L101-L140`.
- 自适应分配：在可行集内按效用比例归一化；若效用全零则均分。  
  证据：`dose_decision.R adaptive_randomization L142-L155`.
- PoC 概率（当前实现 - 基于后验样本）：  
  - 最优剂量=可行集中效用最高的剂量 i*；  
  - 对每个竞争剂量 j，计算成对比较概率：Pr(π_i* > δ_poc × π_j | D_n)，使用边际疗效后验样本；  
  - PoC 概率 = min{成对比较概率}（与最差竞争者的比较）；若仅一个剂量则 PoC=1。  
  - 与阈值 c_poc 比较；若 PoC < c_poc，可能不选择最优剂量（返回 NA）。  
  证据：`dose_decision.R calculate_poc_probability L310-L343`, `check_poc_threshold L346-L361`.
- 最终选择（带 PoC 门控）：  
  - 若 PoC 通过阈值：选择效用最高剂量；  
  - 若 PoC 未通过：返回 optimal_dose=NA，reason="PoC threshold not met"。  
  证据：`dose_decision.R select_final_od_with_poc L364-L425`.

## 5. 工作流 / Workflow (run_trial_simulation)
1) 阶段循环，Stage1 等分分配；Stage2+ 用上一阶段 alloc_probs (main.R L22-L52).  
2) 模拟数据 (main.R L54-L66).  
3) 计算后验 + PAVA/BIVISO + 边际 (main.R L70-L98).  
4) 可行集判定 (main.R L100-L109).  
5) 早停检查，若 A 为空立即返回（all_alloc_probs 为空也可返回）(main.R L110-L126).  
6) 若未终止且非最后阶段，计算下一阶段分配概率 (main.R L132-L143).  
7) 记录分配概率 (main.R L145).  
8) 末尾做最终选择 + PoC (main.R L148-L167).  

## 6. 配置与参数 / Configuration Reality
- 默认 config：5 剂量，5 阶段，cohort=15，φ_T=0.35/c_T=0.5，φ_E=0.1/c_E=0.5，φ_I=0.2/c_I=0.5，c_poc=0.9，delta_poc=0.8，utility_table 如定义。  
  证据：`src/core/config.R L8-L64`.
- 笔记本 / 校准：  
  - `simulation_notebook.qmd`：5 剂量、5 阶段、cohort=15，c_T=c_E=c_I=0.5，φ_T=0.35，φ_E=0.1 (simulation_notebook.qmd L37-L98) — **与默认配置对齐**。  
  - `poc_calibration_notebook.qmd`：5 剂量，c_T=c_E=0.3，c_I=0.2，null 情景校准，阈值更严格以控制 Type I 错误率 (poc_calibration_notebook.qmd L62-L120)。  
  - `src/optimization/poc_calibration_new.R` base_config：5 剂量、校准专用阈值 (c_T=c_E=c_I=0.3) (poc_calibration_new.R L194-L220)。  
- 参数优化：`parameter_optimization.R` 测试多组参数组合：φ_T=[0.25-0.45], φ_E=[0.05-0.25], φ_I=[0.15-0.35], c_T/c_E/c_I=[0.6-0.95]。
- 比较性：默认配置与 simulation notebook 已对齐；校准脚本使用适合 Type I 错误控制的严格阈值。

## 7. 主要文件与函数角色 / File-Level Roles
- `src/core/main.R`: orchestrates full trial workflow; returns all artifacts (run_trial_simulation L15-L167).  
- `src/core/simulate_data.R`: Gumbel copula data generator with optional seed (simulate_data_gumbel L9-L61).  
- `src/core/model_utils.R`: posterior sampling, PAVA/BIVISO, marginals (apply_pava_on_samples L33-L48; apply_biviso_on_matrix L50-L99; compute_marginal_probability L101-L128).  
- `src/decision/dose_decision.R`: utility calc, admissible set, adaptive allocation, PoC, final selection (functions above).  
- `src/optimization/poc_calibration_new.R`: PoC calibration system (create_null_flat_scenario L16-L48; calibrate_c_poc L180-L322; generate_calibration_report L366-L590).
- `src/optimization/parameter_optimization.R`: parameter search framework (run_parameter_optimization L205-L249; create_optimization_plots L252-L326).
- `src/utils/helpers.R` / `plotting_extensions.R`: 结果可视化，非核心算法。  
- Notebooks:  
  - `simulation_notebook.qmd`: 交互式跑一次 trial，绘制后验与分配曲线，生成publication-ready图表。  
  - `poc_calibration_notebook.qmd`: 构造 null/flat 情景，调用 `calibrate_c_poc` 扫 c_poc 候选值，画校准曲线，生成详细报告。

## 8. 统计方法要点 / Key Statistical Details
- Copula: Gumbel 参数 c=rho0/rho1 控制毒性-疗效相关；四格概率显式公式 (simulate_data.R L1-L8).  
- 后验：Beta 共轭，Alpha=1, Beta=1；n_sims=1000 样本用于单调回归与阈值计算 (simulate_beta_posterior L16-L23).  
- 单调约束：PAVA 保持剂量递增；BIVISO 保持剂量和免疫分组的二维单调性。  
- 可行性判定：基于后验样本概率满足 φ/ c 阈值；无多重校正。  
- **PoC 方法（后验样本法）**：使用后验样本直接计算成对比较概率 Pr(π_best > δ × π_competitor | D_n)；PoC=min{成对概率}。基于边际疗效后验样本；完全贝叶斯，无正态近似。PoC 可门控选择：若 PoC < c_poc，则不选最优剂量。  
- **PoC 校准**：通过 null/flat 情景（所有剂量等概率在阈值处）校准 c_poc，目标：10% Type I 错误率（假阳性选择率）。校准方法：对多个 c_poc 候选值运行 1000+ 次模拟，选择最接近 10% PoC 检出率的 c_poc (poc_calibration_new.R L180-L322).  
- 效用：分层免疫的期望效用，utility_table 0–100 加权；自适应分配与最终选择都基于效用。

## 9. 随机性与可重复性 / Randomness
- base seed → 每阶段 seed+stage；阶段内 set.seed 再采样 (main.R L54-L66; simulate_data_gumbel L19-L22)。  
- 每阶段重置 RNG：若重复使用同一 base seed，会跨阶段生成可重复但非全局连续的流。  
- 未提供全局 RNG 控制或并行安全性保证。

## 10. 测试现状 / Tests
- 运行：`R -q -e 'testthat::test_dir("tests")'`.  
- 检查点：结构与字段、分配概率归一化、可行集/效用函数、早停/PoC 场景 (tests/test_main.R; tests/test_dose_decision.R; tests/test_early_termination_poc.R; tests/test_workflow_order.R)。  
- 缺口：未校验 PAVA/BIVISO 数值正确性、PoC 后验样本法准确性、性能优化。

## 11. 已知偏差 vs 文档
- PoC 公式/门控：代码**已更新**为后验样本法，计算 Pr(π_best > δ × π_competitor | D_n)；PoC 可门控选择（若不达标返回 NA）。**状态：RESOLVED**（dose_decision.R L310-L343, L364-L425）。  
- PoC 校准系统：**已实现** null/flat 情景校准框架，目标 10% Type I 错误率。**状态：IMPLEMENTED**（poc_calibration_new.R）。  
- Control arm / γ_j：未实现（代码无 γ_j/对照参数）。**状态：MISSING/OUT-OF-SCOPE**。  
- 配置一致性：默认 config 与 simulation notebook 现已对齐（5 剂量，放宽阈值）；校准用更严格阈值适配 null 情景。**状态：RESOLVED**。

## 12. 校准与优化系统 / Calibration & Optimization Systems

### PoC 校准方法论
- **Null 情景构造**：
  - P_I(j) = φ_I（所有剂量）  
  - P_E(j) = φ_E（所有剂量，用全概率公式保证边际等于阈值）  
  - P_T(j) = 常数安全值（如 0.05）  
  - 所有剂量无差异 → 真 null 情景
- **校准流程**：
  1. 对 c_poc 候选值（如 0.5, 0.6, ..., 0.95）各运行 1000 次模拟  
  2. 记录 PoC 检出率（完成试验且 PoC 验证通过的比例）  
  3. 选择最接近 10% 的 c_poc 作为最优值  
- **输出报告**：
  - 校准曲线图（c_poc vs PoC 检出率）  
  - 详细文本报告，包含早停分析、后验均值诊断、示例案例  
  - 推荐最优 c_poc 及预期试验特性（完成率、早停率）
- **验证步骤（建议）**：在信号情景中测试校准后的 c_poc，确保 ≥50% 真阳性选择率（power）。

### 参数优化框架
- **优化目标**：最大化试验性能指标（完成率、正确选择率、平均效用、分配效率）
- **参数空间**：
  - 阈值：φ_T ∈ [0.25, 0.45], φ_E ∈ [0.05, 0.25], φ_I ∈ [0.15, 0.35]  
  - 可信度：c_T, c_E, c_I ∈ [0.6, 0.95]  
  - 效用表变体：保守、平衡、激进（毒性惩罚不同）
- **评估指标**：
  - completion_rate（试验完成率）  
  - correct_selection_rate（真最优剂量选中率）  
  - mean_final_utility（平均最终效用）  
  - poc_success_rate（PoC 验证通过率）  
  - allocation_efficiency（分配均衡性，1 - CV）
- **工具函数**：
  - `generate_parameter_combinations()`: 生成/采样参数网格  
  - `run_parameter_combination()`: 对单组参数运行 N 次模拟并汇总指标  
  - `create_optimization_plots()`: 可视化参数敏感性与性能曲线  
  - `find_best_parameters()`: 按指定指标筛选最优参数组合

## 13. 下一步（安全增量）
- **验证校准**：用信号情景测试校准后 c_poc，确保 power ≥ 50%。  
- **参数优化**：运行全面优化搜索，为目标情景找到最佳阈值/可信度组合。  
- 若需要对照/γ_j，在 config 中显式标记未支持或添加占位符。  
- 补充 PoC 计算正确性的单元测试（成对比较、样本法）。  
- 性能分析：在大规模校准中 profiling PAVA/BIVISO；当前 n_sims=1000 表现良好。  
- 文档化 seed 策略：stage_seed = base_seed + stage，确保可重复性且无 RNG 冲突。  
- 敏感性分析工具：评估校准参数在情景变化下的稳健性。
