# Diagnostic workflow wave 2 retirement review

创建日期：2026-08-31

状态：`author_review`

本任务承接 wave 1 workflow retirement（PR #296，merge SHA
`15e02b8e8f9b23629908b6ee164d1e858f7ac463`）。目标是让作者审核剩余纯手动 workflow 的
功能、历史证据、可退役理由和参数化合并边界。本任务只做 solver-free 静态审计和文档
固化，不直接修改 `.github/workflows/`。

## 固定边界

- active workflow：43；带 `workflow_dispatch`：39；纯手动入口：25；
- 本轮审阅 17 个，明确保留 8 个；
- 不调用 PNJL equilibrium solver，不删除 run/artifact，不改变 reference 或数值结果；
- 不把结构相似度当作物理合同等价，不自动合并 workflow。

## 作者需要审核的三类

1. 8 个已闭环 Issue #130/CEP/legacy 专用入口：是否直接退役并保留 versioned definition；
2. 7 个 phase shadow/Maxwell local 入口：是否保留一个参数化入口，或整体退役；
3. 2 个 CEP oracle/pilot：是否仍有明确未来重放需求，否则可退役 active trigger。

逐文件功能和理由见：
`docs/analysis/governance/diagnostic_workflow_retirement_wave2_review_v1/manual_workflow_function_review.csv`。

## 验收条件

- [x] 17 个候选的功能、输入、后端脚本、证据和退役理由已记录；
- [x] 8 个 durable retain 入口已明确排除；
- [x] `solver_called=false`、`wave2_changes_executed=false` 已写入 manifest；
- [ ] 作者按组决定 wave 2 的退役、参数化合并或保留；
- [ ] 决策后另立精确 allowlist 的实施 PR。
