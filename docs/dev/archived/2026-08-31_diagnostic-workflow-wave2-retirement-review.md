---
title: Diagnostic workflow wave 2 retirement review
archived: true
original: docs/dev/active/2026-08-31_diagnostic-workflow-wave2-retirement-review.md
archived_date: 2026-08-31
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Diagnostic workflow wave 2 retirement review

创建日期：2026-08-31

状态：`implementation_pr_pending`

本任务承接 wave 1 workflow retirement（PR #296，merge SHA
`15e02b8e8f9b23629908b6ee164d1e858f7ac463`）。目标是让作者审核剩余纯手动 workflow 的
功能、历史证据、可退役理由和参数化合并边界。本任务先做 solver-free 静态审计，再执行
作者批准的精确 allowlist；不调用数值 solver，不删除历史 run/artifact。

## 固定边界

- wave 1 后 active workflow：43；带 `workflow_dispatch`：39；纯手动入口：25；
- 本轮审阅 17 个，其中 16 个退役、1 个保留为参数化入口；
- 不调用 PNJL equilibrium solver，不删除 run/artifact，不改变 reference 或数值结果；
- 不把结构相似度当作物理合同等价，不自动合并 workflow。

## 作者需要审核的三类

1. 8 个已闭环 Issue #130/CEP/legacy 专用入口：是否直接退役并保留 versioned definition；
2. 7 个 phase shadow/Maxwell local 入口：是否保留一个参数化入口，或整体退役；
3. 2 个 CEP oracle/pilot：是否仍有明确未来重放需求，否则可退役 active trigger。

逐文件功能和理由见：
`docs/analysis/governance/diagnostic_workflow_retirement_wave2_review_v1/manual_workflow_function_review.csv`。

## 作者决定（2026-08-31）

- 第一组 8 个 direct-retirement 候选全部退役；
- 4 个 phase-shadow workflow 全部退役；
- CEP deep-oracle 与 narrow-pilot 全部退役，未来若需复用须从 versioned definition 新建明确 scope 的 workflow；
- Maxwell-local 只保留 `pnjl-issue130-maxwell-cep-local-expansion.yml` 作为 target-list/failed-target-only 入口，pilot 和 endpoint-refinement 退役；
- 不存在明确的未来复用需求，不恢复任何已退役入口。

## Wave 2 实施结果

- [x] 16 个 YAML 已移动到
  `docs/analysis/governance/diagnostic_workflow_retirement_wave2_v1/definitions/`，原字节、hash、历史 run/artifact 和测试合同均保留；
- [x] `pnjl-phase-diagram.yml` 的 legacy audit 活动路径已改为 versioned historical definition；旧 evidence 中的原路径仅作历史字符串保留；
- [x] workflow-specific Python/Julia tests 已改读 versioned definitions；保留的 Maxwell-local expansion 未被改写；
- [x] 当前入口计数为 active 27、`workflow_dispatch` 23、纯手动 9；
- [x] evidence package 为
  `docs/analysis/governance/diagnostic_workflow_retirement_wave2_v1/`，manifest 明确 `solver_called=false`、不删除 run/artifact、不改变数值/reference；
- [ ] focused CI 通过并合并实施 PR 后，再将本任务标记为 accepted/closed；

## 验收条件

- [x] 17 个候选的功能、输入、后端脚本、证据和退役理由已记录；
- [x] 8 个 durable retain 入口已明确排除；
- [x] `solver_called=false`、`wave2_changes_executed=true` 已写入 manifest；
- [x] 作者按组决定 wave 2 的退役、参数化合并或保留；
- [x] 已按决定执行精确 allowlist，并生成 versioned evidence package；
- [ ] CI 通过、PR 合并并同步主分支台账。

## 归档时点收口附录（2026-08-31）

以上正文保留 wave2 审阅时点的原始状态。随后作者批准的精确 allowlist 已由 PR #298
合并到 `main@c9c76b12ba5a3199831246665721b54237068b34`：16 个 YAML 已移入
`docs/analysis/governance/diagnostic_workflow_retirement_wave2_v1/definitions/`，
保留一个 Maxwell-local target-list 入口；当前 workflow 计数为 active 27、
`workflow_dispatch` 23、纯手动 9。历史 run/artifact、原始 YAML 字节和 solver-free
provenance 均保留。台账现将本任务标为 `archived`；未来波次必须新建 active task，
不能把本归档任务重新作为执行入口。
