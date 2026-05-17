# DoseFinding 文档入口

这个目录现在只保留当前仍然有用的文档。旧的阶段总结和重复文档已经从工作树删除；需要追溯时可以从 git 历史查看。

## 推荐阅读顺序

1. `CHATGPT_PROJECT_CONTEXT.md`
2. `PROJECT_READALOUD_EXPLANATION.md` 或 `PROJECT_READALOUD_MOBILE.html`
3. `HOW_TO_RUN.md`
4. `STAT_METHODS_AS_BUILT.md`
5. `CODE_MAP.md`
6. `Design1.tex` 和 `Design2.tex`

## 每个文档是干什么的？

- `CHATGPT_PROJECT_CONTEXT.md`：最完整的项目上下文包，适合发给 ChatGPT，让它回答项目相关问题。
- `PROJECT_READALOUD_EXPLANATION.md`：适合念出来听的中文解释稿，用高中生也能听懂的方式讲项目。
- `PROJECT_READALOUD_MOBILE.html`：同一份朗读稿的手机阅读版。
- `HOW_TO_RUN.md`：运行仿真、校准、优化和 notebook 的操作说明。
- `STAT_METHODS_AS_BUILT.md`：以当前代码为准的统计方法说明。
- `CODE_MAP.md`：代码目录和主要文件职责。
- `Design1.tex`：原始模型层设计，重点是免疫反应、毒性、疗效的建模关系。
- `Design2.tex`：原始决策层设计，重点是效用、可接受剂量集、自适应分配、早停和最终 OD 选择。

## 判断优先级

如果文档之间出现冲突，优先级是：

1. 当前代码
2. `CHATGPT_PROJECT_CONTEXT.md`
3. `STAT_METHODS_AS_BUILT.md`
4. `Design1.tex` / `Design2.tex`

旧文档已经不再作为当前事实来源保留在工作树中。
