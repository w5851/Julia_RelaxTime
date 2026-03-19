---
name: citation-management
description: Citation operations skill for DOI/arXiv metadata normalization, BibTeX curation, duplicate cleanup, and reference consistency checks. Use when managing reference quality and bibliography integrity for papers and reports.
license: MIT
metadata:
  author: local-adapted-from-k-dense
  version: "1.0.0"
---

# Citation Management

## When to apply
- 需要把 DOI / arXiv / PMID 统一为规范引用条目。
- 需要清理 `.bib` 中重复、缺字段、格式不一致问题。
- 需要在论文提交前做引用一致性核查。
- 需要按期刊风格导出可用参考文献数据。

## Core workflow
1. 收集候选文献标识符（DOI、arXiv ID 等）。
2. 拉取并校验元数据（作者、题目、期刊、年份）。
3. 规范化 BibTeX 键与字段顺序。
4. 去重并处理版本冲突（预印本 vs 正式发表）。
5. 生成最终参考文献清单并抽样复核。

## Validation checklist
- 作者名拼写与顺序一致。
- 年份、卷期、页码或 article number 完整。
- DOI 可解析，arXiv 编号有效。
- 文内引用键与文末条目一一对应。

## Boundary
- 大规模文献检索与主题综合请配合 `literature-review`。
- 仅做写作润色时优先用 `writing-skills`。

## Hand-off
- 上游通常来自 `deep-research`、`literature-review` 或 `research-engineer` 产出的候选文献集。
- 本 skill 负责把参考文献数据清理到可引用状态，不负责替代主题综述或方法决策。
