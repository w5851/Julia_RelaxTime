---
name: paper-traceable-coauthor
description: 为 Julia_RelaxTime 论文写作生成可追溯草稿，只基于源码、结果文件、图表和明确提供的参考材料。适用于方法、实验设置、图注、结果描述、一致性审校。关键词：paper, methods, caption, evidence, traceable, manuscript
---

# paper-traceable-coauthor

## Purpose

把论文协作从“像论文”收束到“每个关键表述都能追溯到仓库证据”。

本 Skill 适合方法、实验设置、图注、结果说明和术语一致性检查，不负责凭空扩写研究结论。

## Apply When

- 需要写方法小节、实验设置、图表 caption、结果段落。
- 需要检查代码实现、结果文件、图表与论文表述是否一致。
- 需要把实验记录或结果表转成论文可编辑草稿。

## Do Not Use

- 任务核心是广泛文献调研：转 `literature-review` 或 `deep-research`。
- 任务核心是一般技术文档，而非论文/研究写作：转 `doc-coauthoring`。

## Hard Rules

- 不写没有证据支撑的结果或结论。
- 区分事实陈述与解释性表述。
- 单位、变量名、图号、章节引用要与证据一致。
- 证据不足时先列缺口，不强行成文。

## Evidence Priority

优先级从高到低：

1. 相关源码与公开入口
2. 结果表、CSV、数值输出
3. 图表与图注草稿
4. 仓库内现有研究笔记或手稿
5. 用户明确提供的外部参考

若高优先级证据与低优先级草稿冲突，以高优先级证据为准，并显式指出冲突。

## Standard Workflow

1. Bound the writing target
   - 明确是方法、设置、caption、结果、讨论还是一致性审校。
2. Gather traceable evidence
   - 只读取与目标段落直接相关的源码、结果、图表、草稿。
3. Check consistency first
   - 变量名、单位、流程顺序、参数来源、图号与标题。
4. Draft in two layers
   - 先给结构化提纲
   - 再给可编辑正文
5. Flag author-check items
   - 任何推断链过长、证据不足或需要作者判断的点都单列出来。

## Output Contract

输出默认分三段：

- 证据与缺口
- 提纲或结构草稿
- 正文草稿 + 需作者确认项

## Writing Defaults

- 默认使用中文。
- 保持科研写作口吻，不写营销式措辞。
- 结果描述优先“观察到什么”，解释其次；若解释证据不足，明确降级为假设或作者待确认。
