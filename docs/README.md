# DoseFinding 文档入口

这个目录已经按用途重新整理。顶层只保留当前最值得读、最常用的文档；旧的阶段总结和重复解释没有删除，统一放进 `docs/archive/`。

## 先读哪个？

如果你想快速重新理解整个项目，读：

- `CHATGPT_PROJECT_CONTEXT.md`：最完整的项目上下文包。适合发给 ChatGPT，也适合自己用来恢复项目记忆。
- `PROJECT_READALOUD_EXPLANATION.md`：适合念出来听的中文解释稿，用高中生也能听懂的方式讲项目。
- `PROJECT_READALOUD_MOBILE.html`：同一份朗读稿的手机阅读版。

如果你想运行代码，读：

- `HOW_TO_RUN.md`：运行仿真、校准、优化和 notebook 的操作说明。

如果你想理解统计方法，读：

- `STAT_METHODS_AS_BUILT.md`：以当前代码为准的统计方法说明。
- `Design1.tex`：原始模型层设计，重点是免疫反应、毒性、疗效的建模关系。
- `Design2.tex`：原始决策层设计，重点是效用、可接受剂量集、自适应分配、早停和最终 OD 选择。

如果你想找代码位置，读：

- `CODE_MAP.md`：代码目录和主要文件职责。

## 当前推荐阅读顺序

1. `CHATGPT_PROJECT_CONTEXT.md`
2. `PROJECT_READALOUD_EXPLANATION.md` 或 `PROJECT_READALOUD_MOBILE.html`
3. `HOW_TO_RUN.md`
4. `STAT_METHODS_AS_BUILT.md`
5. `CODE_MAP.md`
6. `Design1.tex` 和 `Design2.tex`

## 归档区

`docs/archive/` 里是历史文档和阶段性总结。它们保留了开发过程中的信息，但有些内容已经和当前代码、当前判断不完全一致。以后除非需要追溯历史，不建议优先阅读。

归档文档包括：

- `archive/PROJECT_OVERVIEW.md`
- `archive/PROJECT_STATUS_FINAL.md`
- `archive/PROJECT_ANALYSIS_THREE_QUESTIONS.md`
- `archive/STATUS_AS_BUILT.md`
- `archive/NEXT_STEP_PLAN.md`
- `archive/CALIBRATION_FIX_SUMMARY.md`
- `archive/CALIBRATION_IMPLEMENTATION_SUMMARY.md`
- `archive/PROJECT_GUIDE.md`
- `archive/QUICK_START.md`
- `archive/generated/STAT_METHODS_AS_BUILT.html`

## 整理原则

- 顶层文档只保留“当前入口”和“仍有直接用途”的材料。
- 历史文档不删除，避免丢失以前的判断和开发记录。
- 如果文档之间冲突，优先以当前代码、`CHATGPT_PROJECT_CONTEXT.md` 和 `STAT_METHODS_AS_BUILT.md` 为准。
