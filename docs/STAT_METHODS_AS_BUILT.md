# STAT METHODS AS BUILT (代码为准 / Code-first)

目标：从生物统计角度，基于当前代码描述已实现的模型、算法、公式和工作流（非设计意图）。所有关键点均给出代码证据：`path function Lstart-Lend`。

## 1. 总览 / Overview
本仓库实现了一个带免疫、毒性、疗效三终点的贝叶斯多阶段剂量探索模拟。流程：按阶段模拟数据（Gumbel copula），对各终点做后验估计并用单调约束（PAVA/BIVISO），按阈值筛选可行剂量集，基于期望效用自适应分配，最终结合PoC概率做最优剂量选择。  
证据：`src/core/main.R run_trial_simulation L15-L167`; `src/core/simulate_data.R simulate_data_gumbel L9-L61`; `src/core/model_utils.R apply_pava_on_samples L33-L48`, `apply_biviso_on_matrix L50-L99`; `src/decision/dose_decision.R get_admissible_set L101-L140`, `adaptive_randomization L142-L155`, `calculate_poc_probability L310-L365`, `select_final_od_with_poc L396-L436`.

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
- PoC 概率（当前实现）：  
  - 定义 combined efficacy π_combined = π_I*π_E|I=1 + (1-π_I)*π_E|I=0；  
  - 参考剂量=可行集中效用最高；计算 π_combined_ref；  
  - PoC ≈ 1 - Φ(π_ref*δ_poc - π_combined, sd=0.1)；max_poc 与阈值 c_poc 比较。  
  证据：`dose_decision.R calculate_poc_probability L310-L365`, `check_poc_threshold L367-L394`.
- 最终选择：总是选效用最高剂量；若 PoC 未过阈值仍返回该剂量但标记 reason。  
  证据：`dose_decision.R select_final_od_with_poc L396-L436`.

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
- 默认 config：3 剂量，3 阶段，cohort=6，φ_T=0.3/c_T=0.9，φ_E=0.2/c_E=0.9，φ_I=0.2/c_I=0.8，c_poc=0.9，delta_poc=0.8，utility_table 如定义。  
  证据：`src/core/config.R L8-L60`.
- 笔记本 / 校准：  
  - `simulation_notebook.qmd`：5 剂量、5 阶段、cohort=15，c_T=c_E=c_I=0.5，φ_T=0.35，φ_E=0.1 (simulation_notebook.qmd L37-L98)。  
  - `poc_calibration_notebook.qmd`：5 剂量，c_T=c_E=0.3，c_I=0.2，null 情景校准 (poc_calibration_notebook.qmd L62-L118)。  
  - `src/optimization/poc_calibration_new.R` base_config 亦为 5 剂量、放宽阈值 (poc_calibration_new.R L179-L220)。  
- 比较性：默认与 notebook/校准配置不一致，需对齐后再比较结果。

## 7. 主要文件与函数角色 / File-Level Roles
- `src/core/main.R`: orchestrates full trial workflow; returns all artifacts (run_trial_simulation L15-L167).  
- `src/core/simulate_data.R`: Gumbel copula data generator with optional seed (simulate_data_gumbel L9-L61).  
- `src/core/model_utils.R`: posterior sampling, PAVA/BIVISO, marginals (apply_pava_on_samples L33-L48; apply_biviso_on_matrix L50-L99; compute_marginal_probability L101-L128).  
- `src/decision/dose_decision.R`: utility calc, admissible set, adaptive allocation, PoC, final selection (functions above).  
- `src/utils/helpers.R` / `plotting_extensions.R`: 结果可视化，非核心算法。  
- Notebooks:  
  - `simulation_notebook.qmd`: 交互式跑一次 trial，绘制后验与分配曲线。  
  - `poc_calibration_notebook.qmd`: 构造 null/flat 情景，调用 `calibrate_c_poc` 扫 c_poc，画校准曲线。

## 8. 统计方法要点 / Key Statistical Details
- Copula: Gumbel 参数 c=rho0/rho1 控制毒性-疗效相关；四格概率显式公式 (simulate_data.R L1-L8).  
- 后验：Beta 共轭，Alpha=1, Beta=1；n_sims=1000 样本用于单调回归与阈值计算 (simulate_beta_posterior L16-L23).  
- 单调约束：PAVA 保持剂量递增；BIVISO 保持剂量和免疫分组的二维单调性。  
- 可行性判定：基于后验样本概率满足 φ/ c 阈值；无多重校正或错误率控制。  
- PoC 近似：正态近似差值，固定 sd=0.1，未用严格的 Π_i/Π_ij 贝叶斯比较；PoC 不阻止选择，只影响标记。  
- 效用：分层免疫的期望效用，utility_table 0–100 加权；自适应分配与最终选择都基于效用。

## 9. 随机性与可重复性 / Randomness
- base seed → 每阶段 seed+stage；阶段内 set.seed 再采样 (main.R L54-L66; simulate_data_gumbel L19-L22)。  
- 每阶段重置 RNG：若重复使用同一 base seed，会跨阶段生成可重复但非全局连续的流。  
- 未提供全局 RNG 控制或并行安全性保证。

## 10. 测试现状 / Tests
- 运行：`R -q -e 'testthat::test_dir("tests")'`.  
- 检查点：结构与字段、分配概率归一化、可行集/效用函数、早停/PoC 场景 (tests/test_main.R; tests/test_dose_decision.R; tests/test_early_termination_poc.R; tests/test_workflow_order.R)。  
- 缺口：未校验 PAVA/BIVISO 数值正确性、PoC 近似偏差、性能。

## 11. 已知偏差 vs 文档
- PoC 公式/门控：代码用正态近似并总是选效用最高剂量；TRIAL_DESIGN 要求 Pr(Π_i < δ Π_ij | D_n) 门控。**状态：DRIFT/CONFLICT**（dose_decision.R L310-L365, L396-L436）。  
- Control arm / γ_j：未实现（代码无 γ_j/对照参数）。**状态：MISSING**。  
- 配置不一致：默认 3 剂量 vs notebook/校准 5 剂量，阈值更宽松。**状态：DRIFT**。

## 12. 下一步（安全增量）
- 统一并声明“生产”配置（3剂量或5剂量），同步 notebooks/校准。  
- 明确 PoC 合约：保留近似需说明局限；若需要门控，调整选择逻辑。  
- 若将来要对照/γ_j，先在 config 中显式标记未支持。  
- 补充针对 PAVA/BIVISO 与 PoC 的数值/单调性单元测试。  
- 在文档/测试中记录 seed 策略，避免重复 seed 导致意外重现。
