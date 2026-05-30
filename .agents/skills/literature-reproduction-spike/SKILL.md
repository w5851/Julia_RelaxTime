---
name: literature-reproduction-spike
description: 在 Julia_RelaxTime 中执行隔离式文献复现 Spike。适用于指定论文图表、公式或数值结论的独立 tempN 沙箱复刻、文献事实表、未明示口径审计、复现证据链和二选一 verdict：与文献对齐 / 文献信息不足以复现。关键词：isolated reproduction spike, temp, 文献复现, 独立复现, 复刻, figure reproduction, formula audit, verdict
---

# literature-reproduction-spike

## Purpose

把“复现某篇文献的某个图、表、公式或数值结论”收束成一个可审计的隔离研究 Spike。

本 Skill 的目标不是把曲线调到像，而是在 `tempN/` 沙箱中独立建立证据链，并最终给出二选一结论：

- `与文献对齐`
- `文献信息不足以复现`

## Apply When

- 用户要求复现、复刻、审计或验证论文中的图、表、公式、阈值、相变线、截面、相移、热力学量或数值结论。
- 用户要求“新开 temp”“独立实现”“不要调用项目高层实现”“以文献为准”。
- 任务需要区分文献明示事实、直接推导事实和未明示实现口径。
- 任务结果需要可追溯到文献段落、公式、图注、脚本、CSV、图像和诊断输出。
- 目标是论文 Figure 或曲线时，可配合 `paper-figure-digitize` 做 figure crop、rough digitize、坐标校准或 paper-background overlay。

## Do Not Use

- 只是在主项目中实现已经确定的功能：优先 `julia-pro` 或相关实现 Skill。
- 只做综述、文献列表或方法比较：优先 `literature-review`、`deep-research` 或 `research-engineer`。
- 只做论文写作、图注或结果描述：优先 `paper-traceable-coauthor`。
- 用户明确要求直接修改主线代码，而不是先做隔离复现。

## Hard Rules

- 在根目录新建或使用独立 `tempN/` 沙箱；不要修改 `src/`、`tests/`、`docs/`、`config/` 等主项目文件，除非用户另行明确授权进入正式实现阶段。
- 主项目只读；最多单向参考依赖、通用工具、配置格式、基础 IO、绘图、数值积分或低层 helper。
- 不调用项目已有高层目标实现直接生成目标结果。
- 目标物理量、核心公式、求解流程和后处理逻辑必须按文献重新实现。
- 不把常见做法、经验猜测、作者代码可能做法或项目既有实现写成文献事实。
- 遇到不一致时先审计，不为了贴近文献图而隐式调参。
- `没做出来` 不是 `文献信息不足以复现`；只有完成关键未明示口径审计后才能下该 verdict。

## Sandbox Contract

每个复现沙箱至少包含：

- `README.md`：目标、文献对象、范围边界、非目标、验收标准。
- `literature_facts.md`：文献事实表。
- 主复现脚本和必要绘图脚本。
- `output/`：CSV、图、诊断输出、配置摘要和结论摘要。
- verdict 或 status 文档：记录当前结论、证据和残余 blocker。

推荐命名：

```text
tempN/<paper_or_topic>_independent_audit/
tempN/<paper_or_topic>_reproduction_spike/
```

## Standard Workflow

1. Lock target
   - 明确文献标题、年份、作者、目标 Figure/Table/Eq/数值结论。
   - 明确只复现哪一部分，以及不做哪些扩展。
   - 将目标定义写入 `README.md`。
   - 若目标来自论文图像或曲线，使用 `paper-figure-digitize` 提取 paper crop、rough CSV、overlay 和校准 metadata；这些产物必须写入当前 `tempN/` 沙箱。

2. Build literature fact table
   - 分类记录：`文献明确给出`、`由文献直接推得`、`文献未明示`。
   - 每条注明来源，例如 Eq、Fig caption、page、正文段落或必要上游引用。
   - 未明示项必须作为审计对象，不能默认吞掉。

3. Implement independent minimal loop
   - 在 `tempN/` 内实现最小闭环：背景态、极化函数、传播子、相移、热力学量、积分、根搜索、后处理等按任务需要裁剪。
   - 保留文献符号到代码变量的对应关系。
   - 输出中间诊断量，不只输出最终图。

4. Run first-pass reproduction
   - 先按最忠于文献明文的口径运行。
   - 产出图、CSV、关键中间量摘要和与文献目标的初步对照。
   - 若明显不一致，不直接调参，进入未明示口径审计。

5. Audit unspecified choices
   - 只参数化可能改变拓扑、阈值、量级或分支的关键未明示项。
   - 常见审计项：分支选择、正则化、`eta/iε`、相位 unwrap、Levinson 归一、初值、多解选择、积分范围、网格精度、单位换算、符号约定、是否使用上游引用定义。
   - 做最小必要变体扫描；每个变体记录配置、输出和是否改变结论。

6. Exclude false causes
   - 至少考虑：初值不佳、积分精度不足、网格太粗、单位换算错误、符号号差、相位/分支后处理不一致。
   - 必要时用 multistart、精度敏感性、中间复数实部/虚部、极点邻域和阈值位置诊断交叉验证。

7. Compare and decide
   - 尽量生成文献图裁剪、首轮复现和关键变体的并排对照；图像数字化、坐标映射或 paper-background overlay 可交给 `paper-figure-digitize`。
   - 给出严格二选一 verdict，并说明依据。

## Verdict Gate

只能在满足条件时使用 `与文献对齐`：

- 目标的关键定性特征对齐。
- 关键数值位置、阈值、峰、跳变或曲线形状在合理误差内对齐。
- 所需口径可由文献正文、公式、图注或必要上游引用唯一确定，或少量未明示项不改变结论。
- 没有依赖隐式调参或事后拟合。

只能在满足条件时使用 `文献信息不足以复现`：

- 已完成独立实现。
- 已审计关键未明示口径。
- 多个合理口径会导向显著不同结果。
- 文献文本无法唯一决定作者实际使用哪套口径。
- 因而无法从文献文本唯一且可靠地复现目标。

## Output Contract

最终汇报必须包含：

- verdict：`与文献对齐` 或 `文献信息不足以复现`
- 高密度证据摘要
- `tempN/` 路径
- 关键文件列表
- 已运行配置和变体
- 文献明确给出的信息
- 文献未明示且影响结果的信息
- 若 verdict 为 `文献信息不足以复现`：列出至少 3 条导致不可唯一复现的缺失信息
- 若 verdict 为 `与文献对齐`：列出对齐依据和误差范围

## Hand-off

- 若需要从论文 PDF 中提取曲线、裁剪图像、生成 digitized CSV 或把计算曲线叠加到论文截图背景上，配合使用 `paper-figure-digitize`；其输出应纳入本 Skill 的 `output/` 和 verdict 证据链。
- 若 verdict 通过且用户要求并入主项目，再进入正式实现任务；届时根据影响面选择 `julia-pro`、`transport-regression-keeper`、`baseline-regression-governance` 或 `api-doc-authoring`。
- 若复现产出需要论文表述或图注整理，交给 `paper-traceable-coauthor`。
- 若需要记录实验台账，可追加使用 `experiment-logbook-append`，但不要把沙箱过程史混入 API 或架构文档。
