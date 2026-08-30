---
title: Issue #130：phase-reference runtime consumer compatibility audit 任务单
archived: true
original: docs/dev/active/2026-08-21_Issue130_phase-reference-runtime-consumer-audit.md
archived_date: 2026-08-23
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Issue #130：phase-reference runtime consumer compatibility audit 任务单

状态：accepted（PR253 已合并）。versioned candidate import PR #252 已合并，但 candidate
仍是 `runtime_consumption=false`；本任务完成 solver-free consumer/schema 审计，不切换
默认 reference，不启动 RS transport。后续适配合同由独立任务
`2026-08-23_Issue130_phase-reference-adapter-contract.md` 承接。

## 固定输入

- import merge SHA：`85172be24e5d5ea94239a76fcdd405fdf8e83ab3`。
- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。
- numerical run：`32354095831`；aggregate replay：`32451053476`。
- candidate：`data/reference/pnjl/issue130_phase_reference_v1/`。

## 目标与边界

- [x] candidate 三层表 finite、key unique、输入 hash 与 `runtime_consumption=false` 已核验。
- [x] 旧 reference 文件快照已记录；未修改旧文件。
- [x] gap transport、phase-guided、paper pipeline 和 legacy phase plot 的默认路径与 schema 已盘点。
- [x] verdict：`candidate_isolation_confirmed_requires_explicit_adapter`。
- [x] 结论：candidate 不会被现有 consumer 隐式选择，但不能直接喂给旧 consumer；必须另立显式 adapter/consumer contract。

非目标：不写 legacy reference，不在本任务中实现 adapter，不修改 Models/solver/Maxwell/transport，
不运行 PNJL、C0/C1/C2 或 RS production。adapter contract 和引用迁移不回写本审计包。

## 验收标准

- candidate integrity、旧 reference snapshot、source path checks 全部可追溯。
- `runtime_consumption=false`、`candidate_auto_selected=false`、`solver_called=false`。
- 若 schema 不兼容，verdict 必须为 `candidate_isolation_confirmed_requires_explicit_adapter`，
  并列出后续 adapter 设计边界，而不是自动转换或宣称可以 runtime 使用。
- focused Python test、docs consistency、active-docs、script governance、task ledger 和
  `git diff --check` 通过。

## 后续

作者已授权另立 adapter contract；在 adapter parity、显式 runtime switch 和审核完成前，
RS transport 继续由 `track:rs-transport` 阻塞，不因 import 或本 audit 自动解锁。
