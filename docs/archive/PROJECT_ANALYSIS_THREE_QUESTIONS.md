# 项目分析报告：从三个问题看 DoseFinding 项目

## 摘要

这个项目已经形成了一个贝叶斯自适应剂量探索试验仿真系统。当前代码实现了三终点建模、分阶段入组、后验更新、单调约束、可接受剂量筛选、早期终止、效用驱动随机化、最终 OD 选择以及 PoC 校准框架。

两份设计稿的分工很明确：`Design1.tex` 是统计模型层，说明免疫反应作为 mediator 时如何建模毒性和疗效；`Design2.tex` 是决策层，说明如何用效用、可接受集、组序贯流程和 PoC 规则完成剂量选择。补充提供的两个 PDF 与这两份 `.tex` 内容基本一致，可以作为正式排版版或汇报附件使用。

下一步不应再从零实现核心系统，而应进入“收口和验证”阶段：修复 README 冲突标记和测试路径问题，统一 PoC 定义，运行生产规模 operating characteristics 验证，并明确是否要实现对照臂和 `gamma_j` 分配逻辑。

## 1. 这个项目已经做了哪些事？

### 1.1 已实现的核心能力

项目当前已经实现了一个多阶段贝叶斯剂量寻找模拟器，核心入口是 `src/core/main.R` 中的 `run_trial_simulation()`。

已完成能力包括：

- **三终点数据生成**：模拟免疫反应 `Y_I`、毒性 `Y_T`、疗效 `Y_E`。免疫反应用 Bernoulli 生成，毒性和疗效在给定免疫状态下通过 Gumbel copula 生成相关二元结局。相关代码在 `src/core/simulate_data.R`。
- **贝叶斯后验更新**：用 Beta-Bernoulli 共轭模型，根据每个剂量、每个免疫分层的累计数据抽取后验样本。相关代码在 `src/core/model_utils.R`。
- **单调约束**：免疫反应使用 PAVA，毒性和疗效的条件概率使用 BIVISO，以保证剂量方向和免疫分层方向上的单调性。相关代码在 `src/core/model_utils.R`。
- **边际毒性和疗效计算**：用全概率公式把条件毒性/疗效与免疫概率混合，得到每个剂量的边际毒性和边际疗效。
- **可接受剂量集筛选**：每个剂量需要同时满足安全性、疗效、免疫反应三条后验概率门槛：`P(Tox < phi_T) > c_T`、`P(Eff > phi_E) > c_E`、`P(Imm > phi_I) > c_I`。相关代码在 `src/decision/dose_decision.R`。
- **早期终止**：若中期分析后可接受剂量集为空，试验立即终止，不选择 OD。
- **效用驱动的自适应随机化**：试验继续时，在可接受剂量集中按期望效用比例分配下一阶段患者。
- **最终 OD 选择和 PoC 门控**：最后阶段结束后，先选效用最高的剂量，再用 PoC 阈值判断是否真的选择 OD；若 PoC 不达标，则返回 `NA`。
- **PoC 校准系统**：通过 null/flat 场景扫描 `c_poc`，目标是在无真实差异场景下控制 PoC 检出率或 Type I error。主要实现包括 `src/optimization/poc_calibration_new.R`，同时还保留了较早的 `src/optimization/poc_calibration.R`。
- **早期终止和参数优化框架**：包含早停校准、阈值网格搜索、效用表变体、性能汇总和可视化。相关代码在 `src/optimization/early_termination_calibration.R`、`src/optimization/parameter_optimization.R`、`src/utils/calibration_plots.R`。
- **交互式和示例材料**：提供 Quarto notebooks、示例脚本、结果图表、方法说明和项目指南。

默认配置是 5 个剂量、5 个阶段、每阶段 15 人，总样本量最多 75。默认阈值在 `src/core/config.R` 中定义：`phi_T=0.35`、`phi_E=0.1`、`phi_I=0.20`、`c_T=c_E=c_I=0.5`、`c_poc=0.9`、`delta_poc=0.8`。

### 1.2 代码和文档资产

当前项目结构已经比较完整：

- `src/core/`：主模拟流程、配置、数据生成、贝叶斯后验和单调回归。
- `src/decision/`：效用、可接受集、自适应随机化、早停、PoC 和最终选择。
- `src/optimization/`：PoC 校准、早停校准、参数优化。
- `src/utils/`：结果图、校准图和辅助函数。
- `notebooks/`：单次仿真、PoC 校准、设计 walkthrough。
- `examples/`：基础使用、flat scenario、PoC、综合校准等演示。
- `tests/`：主流程、决策逻辑、PoC、早停、样本量不变量等测试。
- `docs/`：项目概览、代码地图、方法说明、运行指南、设计稿和状态文档。
- `results/`：已有的校准结果、图表和示例输出。

### 1.3 当前验证状态

我做了两类轻量检查：

1. 单次核心仿真可以运行。使用默认配置并关闭 verbose 后，`seed=123` 的一次运行结果为：`final_od=3`、未早停、`poc_validated=TRUE`、总样本量 75。
2. 直接运行 `R -q -e 'testthat::test_dir("tests")'` 时，测试总结果是 `PASS 28`、`FAIL 4`、`WARN 4`。失败原因不是仿真逻辑本身，而是部分测试文件在 `test_dir()` 改变工作目录后仍使用 `source("src/...")` 这种项目根目录相对路径，导致找不到文件。

另外，`README.md` 中当时还残留 Git 合并冲突标记。这属于必须先清理的项目卫生问题；后续文档整理已经清理。

## 2. `Design1.tex` 和 `Design2.tex` 设计了什么？

### 2.1 `Design1.tex`：统计模型设计

`Design1.tex` 设计的是模型层，也就是“数据如何产生、后验如何估计、约束如何施加”。

核心思想是把免疫反应 `Y_I` 作为 mediator：

```text
dose d -> Y_I
Y_I -> Y_T
Y_I -> Y_E
dose d -> Y_T
dose d -> Y_E
```

它把联合分布分解为：

```text
Pr(Y_I, Y_T, Y_E | d)
= Pr(Y_I | d) Pr(Y_T | Y_I, d) Pr(Y_E | Y_I, d)
```

具体设计包括：

- **免疫反应模型**：每个剂量 `j` 上，`Y_I | pi_Ij ~ Bernoulli(pi_Ij)`，`pi_Ij` 使用 Beta 先验，观测后得到 Beta 后验。
- **毒性条件模型**：对每个免疫状态 `I=0/1` 和每个剂量 `j`，建模 `pi_Tj^(I)`，也使用 Beta-Bernoulli 后验。
- **疗效条件模型**：同样对每个免疫状态和剂量建模 `pi_Ej^(I)`。
- **单调性约束**：
  - 免疫反应随剂量做 PAVA 单调调整。
  - 毒性和疗效在“剂量”和“免疫状态”两个方向上做 bivariate isotonic regression。
  - 约束形式是同一免疫层内随剂量递增，同时 `I=1` 层不低于 `I=0` 层。
- **边际化**：最终毒性和疗效用全概率公式混合：

```text
pi_Tj = pi_Ij * pi_Tj^(1) + (1 - pi_Ij) * pi_Tj^(0)
pi_Ej = pi_Ij * pi_Ej^(1) + (1 - pi_Ij) * pi_Ej^(0)
```

代码中的 `simulate_data.R`、`model_utils.R` 和 `main.R` 基本已经落实了这套设计。

### 2.2 `Design2.tex`：决策和试验流程设计

`Design2.tex` 设计的是决策层，也就是“什么时候继续、怎么分配、怎么选 OD”。

主要内容包括：

- **效用表**：根据 `Y_I`、`Y_E`、`Y_T` 的组合给分。例如疗效好、无毒性、免疫反应阳性时分数最高；毒性会降低分数；免疫阳性但无疗效也可能有少量分数。
- **期望效用**：对所有 `Y_I/Y_E/Y_T` 组合求加权和，得到每个剂量的 `U(d_j)`。
- **可接受集**：
  - 无对照臂时，剂量需满足边际安全性、免疫活性和疗效标准。
  - 有对照臂时，安全性和疗效要相对对照臂做比较。
- **组序贯设计**：
  1. 第一阶段对所有剂量等概率随机化。
  2. 从第二阶段起，根据中期数据更新可接受集 `A`。
  3. 若 `A` 为空，则早停且不选 OD。
  4. 若 `A` 非空，则在 `A` 中按后验/效用信息做自适应随机化。
  5. 试验完成后构造 PoC-eligible set，再从其中选择期望效用最高的 OD。
- **对照臂和 `gamma_j` 设想**：设计稿提到在有对照臂时，要控制对照臂分配比例，并用 `gamma_j = Pr(OD=d_j | D_n)` 辅助分配。
- **PoC 校准要求**：在 flat immune/toxicity/efficacy 场景下校准 `C_poc`，控制 familywise Type I error，并保证一定 power。

代码已实现无对照臂版本的效用、可接受集、早停、自适应随机化和最终选择。尚未实现对照臂、`gamma_j`、对照臂分配比例控制。

### 2.3 设计稿和代码的主要一致/偏差

一致部分：

- `Design1.tex` 的 mediator factorization、Beta 后验、PAVA/BIVISO、边际化都已经进入代码。
- `Design2.tex` 的无对照臂 admissible set、utility、group sequential、早停和自适应随机化都已经进入代码。
- PoC 校准的总体方向已经实现：通过 flat/null 场景扫描 `c_poc`，观察 PoC 检出率。

偏差或未定部分：

- **对照臂未实现**：`Design2.tex` 写了有对照臂版本，但当前配置和代码没有 `d0`、`pi_T0`、`pi_E0`、`gamma_j` 或对照分配规则。
- **PoC 定义需要统一**：`Design2.tex` 的公式更像是相对基准剂量或免疫概率的比较；当前代码则以效用最高剂量为 reference，用边际疗效样本计算 `Pr(pi_best > delta * pi_competitor | data)`，再取最小成对概率作为 PoC。
- **文档之间状态不完全一致**：一些文档写“100% 完成”，但 `NEXT_STEP_PLAN.md` 内部仍保留 pending 段落；`CALIBRATION_FIX_SUMMARY.md` 又指出早停校准的解释需要更谨慎。
- **测试启动方式不稳**：测试文件的路径假设不统一，导致标准 `testthat::test_dir("tests")` 入口失败。

### 2.4 两份 PDF 补充材料的作用

本次补充的两份 PDF 是：

- `/Users/jz/Downloads/Dose_Finding_Paper.pdf`
- `/Users/jz/Downloads/Dose_Finding_Paper Part 2.pdf`

两份 PDF 各 2 页，内容分别对应 `Design1.tex` 和 `Design2.tex` 的排版版本。它们没有引入一套新的算法，但能作为更正式的设计说明附件：

- `Dose_Finding_Paper.pdf` 强化了模型部分：mediator factorization、Beta-Bernoulli 后验、PAVA/BIVISO 单调约束、边际化公式。
- `Dose_Finding_Paper Part 2.pdf` 强化了决策部分：效用表、无对照/有对照 admissible set、组序贯流程、PoC-eligible set、对照臂分配原则。
- Part 2 PDF 明确写到 `C_PoC` 应在 flat immune response / toxicity / efficacy 曲线下校准，并控制 familywise Type I rate 为 `0.05`，同时希望 power 达到 `80%-90%` 或按实际情况降低。

这个 PDF 信息带来一个需要对齐的重点：当前代码和已有文档多处把 PoC 校准目标写成约 `10%` detection rate，而 PDF/设计稿写的是 familywise Type I rate `0.05`。这两个目标不一定等价，下一步需要由统计负责人确认最终采用哪个错误率定义和目标值。

## 3. 下一步要做什么？

### 3.1 第一优先级：项目卫生和可运行性收口

这些是最应该先做的短平快事项：

1. **清理 README 合并冲突**  
   `README.md` 当时仍含 Git 合并冲突标记。这会让项目看起来未合并完成，也可能误导使用者；后续整理已经清理。

2. **修复测试路径问题**  
   统一测试文件里的 source 路径。建议在测试 helper 中定义 project root，再用 `file.path(project_root, ...)` 引用源码。目标是让下面命令稳定通过：

   ```r
   testthat::test_dir("tests")
   ```

3. **统一校准入口**  
   当前同时存在 `poc_calibration.R` 和 `poc_calibration_new.R`。需要明确哪个是 canonical 版本，另一个要么归档、要么改名为 legacy，避免 notebooks、examples、docs 调用不同版本。

4. **同步状态文档**  
   把 `NEXT_STEP_PLAN.md`、`CALIBRATION_IMPLEMENTATION_SUMMARY.md`、`STATUS_AS_BUILT.md`、`README.md` 中互相矛盾的完成状态统一。建议以代码和最新验证结果为准。

### 3.2 第二优先级：统计定义和设计-代码对齐

1. **明确 PoC 的统计定义**  
   需要决定 PoC 到底比较什么：
   - 免疫反应概率？
   - 边际疗效概率？
   - 综合效用？
   - 相对基准剂量 `d1`/对照臂？
   - 还是相对当前效用最高剂量？

   这个决定会影响 `calculate_pi_parameters()`、`calculate_poc_probability()`、测试和论文/报告里的公式。

2. **确认 PoC 校准目标是 0.05 还是 0.10**  
   PDF/设计稿写的是 familywise Type I rate `0.05`；当前代码和项目文档常用的是 null/flat 场景下约 `10%` PoC detection rate。建议先明确统计目标，再决定 `calibration_config$poc_target_rate` 和报告口径。

3. **把 Design2 的 PoC 公式改成代码实际公式，或反过来改代码**  
   现在二者不是完全同一个规则。建议先和统计负责人确认，然后让设计稿、代码注释、测试断言保持一致。

4. **明确对照臂是否属于当前项目范围**  
   如果短期不实现，就在 README、config 和方法文档中显式写“当前仅支持无对照臂设计”。如果要实现，则需要补充 `d0`、对照臂后验、相对安全/疗效 admissibility、`gamma_j`、对照分配规则和相应测试。

### 3.3 第三优先级：生产规模 operating characteristics 验证

当前代码已经能跑，但“能跑”不等于“设计性能已经被证明”。下一步应该系统验证 operating characteristics。

建议至少输出以下场景：

- **Null/flat 场景**：所有剂量无真实差异，用于 Type I error / PoC false positive control。
- **Unfavorable 场景**：所有剂量临床上不可接受，用于检查早停和 futility control。
- **Favorable/signal 场景**：存在明确最优剂量，用于检查 power、正确选择率和平均样本量。
- **边界场景**：疗效、免疫、毒性接近阈值，用于检查敏感性。

建议报告以下指标：

- 早期终止率和终止阶段分布。
- 试验完成率。
- PoC 检出率。
- 最终 OD 选择率。
- 正确选择率。
- 平均最终效用。
- 各剂量平均入组人数。
- 平均样本量。
- Type I error 是否达到目标，如 PDF 中的 familywise `0.05`，或代码当前采用的 `0.10` PoC detection target。
- 有信号场景下 power 是否达到预设目标，如 80% 或统计负责人指定标准。

### 3.4 第四优先级：补足关键测试

建议新增或强化：

- PAVA 数值正确性测试：单调输入不变，非单调输入被正确 pool。
- BIVISO 数值正确性测试：剂量方向和免疫层方向同时满足约束。
- PoC 成对比较测试：构造人工 posterior samples，验证 PoC 概率和 `P_final`。
- 最终选择门控测试：PoC 通过选 OD，不通过返回 `NA`。
- seed 可重复性测试：同一 seed 完全复现，不同 seed 产生合理差异。
- notebooks/examples 的 smoke test：至少保证主要入口不因 source 路径失败。

### 3.5 推荐的下一步执行路线

**第 1 步：清理仓库**

- 修 README 冲突。
- 修测试路径。
- 确定 canonical PoC 校准脚本。
- 更新 README/状态文档。

**第 2 步：定稿统计规则**

- 确认 PoC 公式。
- 确认 PoC 校准目标：familywise Type I `0.05` 还是 PoC detection `0.10`。
- 确认是否支持对照臂。
- 让 `Design2.tex`、代码注释、测试和文档统一。

**第 3 步：跑小规模 OC 验证**

- 每个典型场景先跑 100-500 次仿真。
- 快速查看早停率、PoC 检出率、正确选择率是否符合直觉。
- 修正明显不合理的阈值和场景设置。

**第 4 步：跑生产规模 OC 验证**

- 每个关键场景跑 1000-10000 次仿真。
- 输出正式图表和 CSV/HTML 报告。
- 固化推荐参数，如 `c_poc`、`c_T/c_E/c_I`。

**第 5 步：准备交付**

- 完成一份最终统计方法报告。
- 完成一份用户运行指南。
- 确保 tests 全绿。
- 标记项目当前版本和已知不支持范围。

## 结论

这个项目的核心实现已经完成，重点不再是“有没有系统”，而是“系统的统计规则是否统一、验证是否充分、文档是否可信”。从现在开始，最有价值的工作是把设计稿、代码、测试和 operating characteristics 报告打通，形成一个可交付、可复现、可解释的剂量探索仿真工具。
