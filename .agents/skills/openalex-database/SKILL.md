```skill
---
name: openalex-database
description: OpenAlex-focused scholarly retrieval skill for work/author/institution queries, citation trend snapshots, and topic-level evidence discovery. Use when you need structured literature discovery and bibliometric signals from OpenAlex.
license: MIT
metadata:
  author: local-adapted-from-k-dense
  version: "1.0.0"
---

# OpenAlex Database

## When to apply
- 需要按主题快速检索高相关论文与引用网络。
- 需要按作者、机构、领域查看研究产出与趋势。
- 需要获取开放获取状态与基础文献计量信号。
- 需要为综述或选题提供结构化候选文献池。

## Core workflow
1. 设定查询主题与过滤条件（时间、领域、类型）。
2. 抓取 works / authors / institutions 结果。
3. 计算基础指标（被引、近年趋势、主题聚类）。
4. 输出候选阅读列表与优先级。
5. 将结果交给 `literature-review` 做系统综合。

## Output expectations
- Top 文献清单（含标识、年份、来源）。
- 主题分组与代表性文献。
- 近年趋势与潜在研究空白。
- 后续筛选建议（哪些文献需精读/剔除）。

## Boundary
- 该技能偏“检索与计量”，不替代深度方法评估。
- 方案落地与技术路线选择交给 `research-engineer`。
```