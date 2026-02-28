```skill
---
name: literature-review
description: Structured multi-database literature review workflow for academic synthesis with reproducible search logs and thematic evidence mapping. Use when you need rigorous review methodology (scope, screening, synthesis), not general brainstorming.
license: MIT
metadata:
  author: local-adapted-from-k-dense
  version: "1.0.0"
---

# Literature Review

## When to apply
- 需要系统综述、范围综述、或研究现状盘点。
- 需要记录可复现的检索策略（数据库、时间窗、关键词）。
- 需要把文献按主题综合，而非逐篇罗列。
- 需要输出“结论 + 证据强度 + 局限性”。

## Core workflow
1. 定义研究问题与纳入/排除标准。
2. 设计多源检索策略（例如 arXiv / OpenAlex / Scholar）。
3. 去重并进行标题-摘要-全文分层筛选。
4. 提取关键字段（方法、数据、结论、限制）。
5. 按主题或方法路径做综合比较。
6. 输出可追踪结论与后续研究缺口。

## Quality checklist
- 明确检索日期、范围、关键词组合。
- 说明排除原因，避免选择偏差。
- 区分“证据事实”与“作者推断”。
- 给出冲突证据与不确定性来源。

## Boundary
- 通用跨域背景搜集优先用 `deep-research`。
- 方法选型并落地到工程决策优先用 `research-engineer`。
```