# DoseFinding 项目上下文包：给 ChatGPT 读取后回答项目问题

这份文档的目的：把 DoseFinding 项目的背景、目标、代码结构、统计设计、参数含义、当前完成状态、已知问题和下一步计划集中到一个文件里。把这份文档发给 ChatGPT 后，它应该能回答大多数和本项目有关的问题。

## 1. 项目一句话说明

DoseFinding 是一个用 R 写的贝叶斯自适应剂量寻找试验模拟项目。它模拟一个分阶段临床试验，在不同剂量下观察免疫反应、毒性和疗效三个二元终点，然后根据贝叶斯后验、单调约束、可接受剂量筛选、效用函数、自适应随机化和 PoC 门控来决定试验是否继续、下一阶段如何分配病人，以及最终是否选择一个最优剂量。

更口语化地说：这个项目想在电脑里反复模拟临床试验，帮助判断“哪个剂量在安全、有效、有免疫反应之间最平衡”。

## 2. 项目的核心目标

这个项目不是单纯找最大耐受剂量 MTD。它更接近寻找 Optimal Dose，简称 OD。

目标剂量不一定是最高剂量，也不一定只是毒性刚好可接受的剂量。它要综合考虑：

- 免疫反应是否足够；
- 疗效是否足够；
- 毒性是否可接受；
- 三者合起来的风险-收益效用是否最高；
- 最终证据是否足够强，能不能通过 PoC 检查。

最终目标是形成一个可信的剂量探索试验设计工具，用于评估不同试验参数在不同模拟场景下的 operating characteristics，例如早停率、正确选择率、PoC 通过率、平均样本量、各剂量入组人数等。

## 3. 三个核心终点

项目使用三个二元终点：

- `Y_I`：immune response，免疫反应。`1` 表示有免疫反应，`0` 表示没有。
- `Y_T`：toxicity，毒性。`1` 表示出现毒性，`0` 表示没有。
- `Y_E`：efficacy，疗效。`1` 表示有效，`0` 表示无效。

项目的一个关键设定是：免疫反应 `Y_I` 是 mediator，中介变量。剂量会影响免疫反应，免疫反应又会影响毒性和疗效。同时剂量也可能直接影响毒性和疗效。

概念图：

```text
dose d -> Y_I
Y_I -> Y_T
Y_I -> Y_E
dose d -> Y_T
dose d -> Y_E
```

## 4. 默认试验设置

默认配置在 `src/core/config.R`。

默认试验规模：

- `dose_levels = c(1, 2, 3, 4, 5)`：5 个剂量；
- `n_stages = 5`：5 个阶段；
- `cohort_size = 15`：每阶段 15 个病人；
- 最大样本量：`5 * 15 = 75`。

默认判断阈值：

- `phi_T = 0.35`：毒性上限；
- `phi_E = 0.1`：疗效下限；
- `phi_I = 0.20`：免疫反应下限；
- `c_T = 0.5`：毒性安全性的后验概率门槛；
- `c_E = 0.5`：疗效达标的后验概率门槛；
- `c_I = 0.5`：免疫反应达标的后验概率门槛；
- `c_poc = 0.9`：PoC 证据门槛；
- `delta_poc = 0.8`：PoC 成对比较时的比例参数。

默认模拟场景：

```r
p_YI <- c(0.10, 0.30, 0.50, 0.60, 0.70)
```

`p_YI` 是各剂量下真实免疫反应概率。

`p_YT_given_I` 是各剂量、各免疫状态下真实毒性概率。

`p_YE_given_I` 是各剂量、各免疫状态下真实疗效概率。

`rho0` 和 `rho1` 控制在 `I=0` 和 `I=1` 人群中，毒性和疗效之间的相关性。项目使用 Gumbel copula 生成相关二元结局。

## 5. 大流程：一次模拟试验如何运行

核心入口是 `src/core/main.R` 里的 `run_trial_simulation()`。

一次完整模拟大致如下：

1. 读取配置和真实场景参数。
2. 第一阶段把病人尽量平均分配到所有剂量。
3. 用 `simulate_data_gumbel()` 模拟这一阶段每个病人的 `Y_I`、`Y_T`、`Y_E`。
4. 累积当前所有阶段数据。
5. 用 Beta-Bernoulli 后验估计免疫、毒性、疗效概率。
6. 对免疫概率做 PAVA 单调调整。
7. 对毒性和疗效的条件概率做 BIVISO 二维单调调整。
8. 用全概率公式计算每个剂量的边际毒性和边际疗效。
9. 根据后验概率和阈值筛选可接受剂量集。
10. 如果可接受剂量集为空，则早期终止，不选择 OD。
11. 如果可接受剂量集不为空，计算每个可接受剂量的期望效用。
12. 非最后阶段时，按效用比例更新下一阶段的随机化分配概率。
13. 重复上述阶段直到试验结束或提前终止。
14. 如果试验完成，在最终可接受剂量集中选择效用最高剂量。
15. 用 PoC 检查最终证据是否足够强。
16. 若 PoC 通过，返回最终 OD；若 PoC 不通过，返回 `NA`，表示不选择 OD。

## 6. 程序如何判断一个剂量好不好

程序分三关判断。

### 第一关：是否进入可接受剂量集

一个剂量必须同时满足：

```text
P(Tox < phi_T | data) > c_T
P(Eff > phi_E | data) > c_E
P(Imm > phi_I | data) > c_I
```

毒性是越低越好，所以看是否低于 `phi_T`。疗效和免疫反应越高越好，所以看是否高于 `phi_E` 和 `phi_I`。

### 第二关：期望效用是否高

对可接受剂量，程序根据 `utility_table` 计算期望效用。效用表给不同结果组合打分，例如：

- 有疗效、没毒性、有免疫反应：最高分；
- 有疗效、没毒性、无免疫反应：高分；
- 有疗效但有毒性：扣分；
- 没疗效且有毒性：低分；
- 有免疫反应但暂无疗效：可能有少量分。

程序用每个剂量下不同结果组合的后验概率乘以分数，再求和，得到期望效用。

### 第三关：PoC 是否通过

即使某个剂量期望效用最高，程序还要判断证据是否足够强。当前代码大致做法是：以效用最高的剂量为 reference，用边际疗效后验样本与竞争剂量做成对比较，计算类似：

```text
Pr(pi_best > delta_poc * pi_competitor | data)
```

然后取最差竞争者比较结果作为 PoC 概率。如果 PoC 概率达到 `c_poc`，才选择 OD。否则返回 `NA`。

注意：PoC 的精确定义是当前项目最需要统一的地方之一。设计稿和代码的写法并不完全一致。

## 7. 主要参数作用和关系

### `p` 开头参数：模拟真实世界

`p_YI`、`p_YT_given_I`、`p_YE_given_I` 是模拟时设定的真实概率。它们不是试验已经知道的东西，而是电脑用来生成虚拟病人的“真实世界”。

- `p_YI`：真实免疫反应概率；
- `p_YT_given_I`：给定免疫状态后的真实毒性概率；
- `p_YE_given_I`：给定免疫状态后的真实疗效概率；
- `rho0`, `rho1`：毒性和疗效之间的相关结构参数。

### `phi` 参数：临床门槛

- `phi_T`：毒性上限，毒性低于它才算安全；
- `phi_E`：疗效下限，疗效高于它才算有效；
- `phi_I`：免疫反应下限，免疫反应高于它才算足够。

### `c` 参数：后验把握程度

- `c_T`：需要多有把握相信毒性安全；
- `c_E`：需要多有把握相信疗效达标；
- `c_I`：需要多有把握相信免疫反应达标。

`c` 越高，试验越保守，更容易排除剂量或早停。`c` 越低，试验越宽松。

### `utility_table`：风险收益价值观

`utility_table` 定义不同结果组合的分数，是项目用来权衡疗效、毒性和免疫反应的价值表。

### `c_poc` 和 `delta_poc`：最终证据门槛

- `delta_poc`：PoC 成对比较里的比例参数；
- `c_poc`：最终 PoC 概率必须超过的门槛。

整体关系：

```text
p 参数生成虚拟病人
    ↓
贝叶斯后验估计 pi_I, pi_T, pi_E
    ↓
phi 和 c 筛选可接受剂量集
    ↓
utility_table 给可接受剂量打综合分
    ↓
按效用自适应分配下一阶段病人
    ↓
c_poc 和 delta_poc 决定最后是否选择 OD
```

最容易混淆的例子：

```text
p_YI 是模拟世界里的真实免疫概率
phi_I 是试验判断免疫是否达标的标准
c_I 是我们要多有把握相信免疫超过 phi_I
```

毒性和疗效同理。

## 8. 代码结构

### `src/core/config.R`

配置文件。定义试验规模、阈值、PoC 参数、效用表、默认真实场景、flat/null 场景和 calibration 设置。

### `src/core/simulate_data.R`

数据生成。核心函数：

- `Gumbel()`：根据毒性、疗效边际概率和相关参数生成四格联合概率；
- `simulate_data_gumbel()`：按剂量和免疫状态生成病人级别的 `Y_I`、`Y_T`、`Y_E`；
- flat scenario 相关函数：用于 PoC 校准的 null/flat 场景生成。

### `src/core/model_utils.R`

贝叶斯后验和单调约束工具。核心功能：

- `compute_rn()`：计算成功数和样本量；
- `simulate_beta_posterior()`：Beta 后验抽样；
- `apply_pava_on_samples()`：对免疫概率做 PAVA；
- `apply_biviso_on_matrix()`：对毒性/疗效条件概率做 BIVISO；
- `compute_marginal_probability()`：把条件概率和免疫概率混合为边际概率。

### `src/decision/dose_decision.R`

决策逻辑。核心功能：

- `get_expected_utility()`：计算期望效用；
- `get_admissible_set()`：筛选可接受剂量集；
- `adaptive_randomization()`：按效用比例分配下一阶段病人；
- `check_early_termination()`：判断是否早停；
- `calculate_poc_probability()`：计算 PoC；
- `select_final_od_with_poc()`：最终选择 OD 或返回 `NA`。

### `src/core/main.R`

总流程。核心函数：

- `run_trial_simulation()`：串起数据生成、后验估计、可接受集、早停、自适应随机化、最终选择和 PoC。

### `src/optimization/`

校准和优化。

- `poc_calibration.R`：正式 PoC 校准系统；
- `early_termination_calibration.R`：早停校准；
- `parameter_optimization.R`：参数搜索；
- `run_optimization.R`：优化运行入口。

注意：PoC 校准现在只有一个正式入口：`src/optimization/poc_calibration.R`。

### `notebooks/`

交互式 Quarto 笔记本，例如单次仿真、PoC 校准和设计 walkthrough。

### `tests/`

测试文件，覆盖主流程、剂量决策、PoC、早停、样本量不变量等。

当前测试入口：`testthat::test_dir("tests")`。

## 9. Design1 和 Design2 的作用

### `docs/Design1.tex`

`Design1.tex` 是统计模型蓝图。

它回答：如何从数据估计免疫、毒性和疗效？

主要内容：

- mediator factorization：

```text
Pr(Y_I, Y_T, Y_E | d)
= Pr(Y_I | d) Pr(Y_T | Y_I, d) Pr(Y_E | Y_I, d)
```

- `Y_I` 使用 Beta-Bernoulli；
- `Y_T | Y_I` 使用分层 Beta-Bernoulli；
- `Y_E | Y_I` 使用分层 Beta-Bernoulli；
- 免疫概率用 PAVA；
- 毒性和疗效条件概率用 BIVISO；
- 用全概率公式得到边际毒性和边际疗效：

```text
pi_Tj = pi_Ij * pi_Tj^(1) + (1 - pi_Ij) * pi_Tj^(0)
pi_Ej = pi_Ij * pi_Ej^(1) + (1 - pi_Ij) * pi_Ej^(0)
```

代码对应：`simulate_data.R`、`model_utils.R`、`main.R`。

### `docs/Design2.tex`

`Design2.tex` 是试验决策蓝图。

它回答：估计完概率以后，试验怎么继续、怎么停、怎么分配、怎么选最终剂量？

主要内容：

- 效用表；
- 期望效用；
- 无对照臂的可接受剂量集；
- 有对照臂的可接受剂量集；
- 组序贯试验流程；
- 自适应随机化；
- 早期终止；
- PoC-eligible set；
- 最终 OD 选择；
- 对照臂和 `gamma_j` 的分配想法。

代码当前主要实现的是无对照臂版本。对照臂、`gamma_j` 和对照臂分配比例尚未真正实现。

## 10. simulate 和 Design1/Design2 的关系

一句话：

```text
Design1 = 模型层蓝图
Design2 = 决策层蓝图
simulate/main code = 按这两个蓝图跑出来的虚拟试验系统
```

`simulate_data.R` 负责模拟虚拟病人的结果。它更接近 Design1 的数据生成部分。

`model_utils.R` 负责后验估计和单调约束。它落实 Design1 的统计模型。

`dose_decision.R` 和 `main.R` 负责可接受集、早停、自适应分配、最终选择和 PoC。它们落实 Design2 的决策流程。

所以 simulate 不是和 Design1/Design2 并列的第三个设计。simulate 是把 Design1/Design2 里的想法变成可以反复运行的电脑试验。

## 11. 两份项目内部 PDF

用户提供了：

- `/Users/jz/Downloads/Dose_Finding_Paper.pdf`
- `/Users/jz/Downloads/Dose_Finding_Paper Part 2.pdf`

这两份 PDF 各 2 页，内容基本对应 `Design1.tex` 和 `Design2.tex` 的正式排版版。

重点补充信息：

- 第一份 PDF 对应模型层：mediator factorization、Beta-Bernoulli、PAVA/BIVISO、边际化。
- 第二份 PDF 对应决策层：效用、可接受集、组序贯、PoC、对照臂想法。
- 第二份 PDF 明确提到 `C_PoC` 应在 flat immune response / toxicity / efficacy 曲线下校准，并控制 familywise Type I rate 为 `0.05`，希望 power 达到 `80%-90%` 或按实际情况调整。

注意：当前代码和部分文档常使用 `0.10` PoC detection target，这和 PDF 中的 `0.05` familywise Type I rate 需要进一步确认和统一。

## 12. 两篇参考论文和本项目关系

项目开始时参考了两篇 BOIN 相关论文：

- `jrsssc_gBOIN.pdf`
- `BOIN_ETC.pdf`

### gBOIN

gBOIN 是 generalized Bayesian optimal interval design。它主要解决一期剂量探索，核心目标通常是找 MTD，即最大耐受剂量。

它的主要思想：

- model-assisted，既有模型依据，又尽量简单可执行；
- 不需要每次都重新拟合复杂模型；
- 用预先确定的升剂量/降剂量边界做决策；
- 可处理毒性等级、二元或连续终点；
- 最终可用单调回归选择 MTD。

本项目与 gBOIN 的相似点：

- 都是剂量探索；
- 都关心病人安全；
- 都是贝叶斯/自适应大方向；
- 都使用单调性思想。

本项目与 gBOIN 的不同点：

- gBOIN 主要找 MTD，本项目找 OD；
- gBOIN 主要围绕毒性，本项目同时看免疫、毒性、疗效；
- gBOIN 更像区间决策，本项目更像效用驱动决策；
- 本项目加入了免疫 mediator 和 PoC 门控。

### BOIN-ETC

BOIN-ETC 是 Bayesian optimal interval design considering efficacy and toxicity for dose combinations。它主要面向抗癌药联合治疗，寻找最优剂量组合 ODC。

它的主要思想：

- 同时考虑疗效和毒性；
- 不只是找最大耐受剂量；
- 使用 model-assisted 思路；
- 使用 waterfall approach 把二维剂量组合问题拆成一系列一维子问题；
- 最终用 utility 衡量风险-收益权衡。

本项目与 BOIN-ETC 的相似点：

- 都不满足于只找 MTD；
- 都关心疗效和毒性的风险-收益平衡；
- 都使用效用思想；
- 都带有 adaptive/model-assisted 的风格。

本项目与 BOIN-ETC 的不同点：

- BOIN-ETC 处理二维剂量组合，本项目默认是单轴 5 剂量；
- BOIN-ETC 核心终点是疗效和毒性，本项目多了免疫反应；
- BOIN-ETC 用 waterfall approach，本项目没有 waterfall；
- 本项目使用免疫 mediator、PAVA/BIVISO 和 PoC 门控。

简单总结：

```text
gBOIN 给了本项目安全性控制和 BOIN 思想背景。
BOIN-ETC 给了本项目疗效-毒性风险收益权衡背景。
本项目不是复刻任何一篇，而是面向免疫-毒性-疗效三终点场景的自定义模拟器。
```

## 13. 当前完成状态

已完成或基本完成：

- 核心单次试验模拟；
- 三终点数据生成；
- Beta 后验估计；
- PAVA/BIVISO 单调约束；
- 边际毒性/疗效计算；
- 可接受剂量筛选；
- 早期终止；
- 效用计算；
- 自适应随机化；
- 最终选择；
- PoC 门控；
- PoC 校准框架；
- 早停校准框架；
- 参数优化框架；
- notebooks、examples、tests 和多个文档。

轻量验证结果：

- 默认配置下，关闭 verbose 后，用 `seed=123` 跑一次 `run_trial_simulation()` 可以完成；
- 该次运行结果：`final_od = 3`，未早停，`poc_validated = TRUE`，样本量 75。

## 14. 当前主要问题

### 文档和代码入口已经清理

旧的重复文档已经从工作树移除，当前入口是 `README.md` 和 `docs/README.md`。PoC 校准入口也已统一到 `src/optimization/poc_calibration.R`，旧的 `_new` 版本不再保留。

### 测试入口

推荐运行：

```r
testthat::test_dir("tests")
```

测试目录现在包含 `tests/helper-project-root.R`，用于稳定项目根目录相对路径。

### PoC 定义和目标需要统一

当前代码的 PoC 实现、`Design2.tex` 的公式、PDF 中的校准目标，以及文档中的 `0.10` target 并不完全统一。

需要确认：

- PoC 比较的是免疫概率、边际疗效、效用，还是其他指标？
- reference 是 dose 1、control、效用最高剂量，还是其他？
- Type I error 目标是 `0.05` familywise rate，还是 `0.10` PoC detection rate？

### 对照臂未实现

Design2 中有对照臂版本和 `gamma_j` 分配想法。但当前代码主要实现无对照臂版本。

如果短期不实现，应在文档中明确写“当前版本只支持无对照臂设计”。如果要实现，需要新增 control arm 后验、相对 admissibility、`gamma_j` 和对照臂分配规则。

## 15. 下一步建议

优先级 1：项目卫生和可运行性

- 运行完整测试并记录失败点；
- 确认生成物是否需要重新渲染；
- 重新跑一次最小 operating-characteristics 验证。

优先级 2：统计规则定稿

- 统一 PoC 公式；
- 统一 PoC 校准目标；
- 决定是否支持对照臂；
- 让 Design2、代码注释、测试和文档一致。

优先级 3：运行 operating characteristics

至少设计并模拟：

- null/flat 场景；
- unfavorable 场景；
- favorable/signal 场景；
- 边界场景。

每类场景建议跑 1000 次或更多，报告：

- 早期终止率；
- 终止阶段分布；
- 试验完成率；
- PoC 通过率；
- 最终 OD 选择率；
- 正确选择率；
- 平均最终效用；
- 平均样本量；
- 各剂量平均入组人数；
- Type I error / power。

优先级 4：补充关键测试

- PAVA 数值正确性；
- BIVISO 数值正确性；
- PoC 成对比较；
- 最终选择门控；
- seed 可重复性；
- notebooks/examples smoke test。

## 16. 如果 ChatGPT 被问到常见问题，应这样回答

### 这个项目是干什么的？

回答：它是一个贝叶斯自适应剂量寻找试验模拟器，用来模拟分阶段临床试验，并综合免疫、毒性和疗效选择最优剂量。

### 它是不是 BOIN？

回答：不是直接复现 BOIN。它受 gBOIN 和 BOIN-ETC 启发，但加入了免疫 mediator、三终点建模、效用驱动分配和 PoC 门控。

### Design1 是什么？

回答：Design1 是模型层，定义免疫、毒性、疗效怎么建模，怎么做后验估计和单调约束。

### Design2 是什么？

回答：Design2 是决策层，定义效用、可接受集、组序贯流程、自适应随机化、早停、PoC 和最终 OD 选择。

### simulate 和 design 的关系是什么？

回答：Design1/2 是设计蓝图，simulate/main code 是把蓝图跑起来的虚拟试验系统。

### 程序怎么判断剂量好不好？

回答：先看是否进入可接受剂量集，再看期望效用是否最高，最后看 PoC 证据是否足够强。

### 当前做到哪里了？

回答：核心模拟系统已经实现并能跑通，文档入口和 PoC 校准入口已经清理。下一步主要是统一 PoC 统计目标，并运行生产规模 operating characteristics。

### 最大风险是什么？

回答：不是主流程不能跑，而是统计规则还需要最终定稿，特别是 PoC 校准目标、是否支持对照臂这两个问题。

## 17. 最短总结

这个项目是在模拟一种更聪明的临床剂量探索试验。它边做边学，看到某个剂量安全、有效、有免疫反应，就多分配病人；看到所有剂量都不行，就早停；最后只有在证据足够强时，才选择最优剂量。

它的主体已经搭好，但还需要整理说明书、修测试、统一 PoC 规则、跑大量模拟验证。完成这些后，它才会从“能运行的研究代码”变成“可信的试验设计工具”。
