---
name: codex-task-harness
description: 处理 Julia_RelaxTime 目标不明的任务续接、主线选择或跨主线协调，从指定线索及 task ledger 恢复当前任务、阻塞依赖和下一步。
---

# 跨主线恢复与分诊

从 `config/governance/task_tracks.toml` 恢复主线、当前任务和依赖；字段、状态、晋升条件和记录边界遵循[任务跟踪治理](../../../docs/dev/task_tracking_governance.md)。

## 恢复

1. 用户指定对话、文件或 issue 时先从该线索恢复目标；范围已经明确且不涉及跨主线协调时，检查相关 Git 状态后直接执行。
2. 仍有主线歧义或需要 ledger 协调时，运行 `julia --project=. scripts/dev/check_task_ledger.jl --preflight`；已知目标时加 `--track <id>`。预检摘要包含阻塞、下一步和证据位置；需要完整脏路径时加 `--full-paths`。
3. 只读选定 track、current item、阻塞依赖与必要证据，说明剩余工作、本轮范围和验证。ledger 校验失败仍可读取报告中的 Git 上下文，但不能把未经校验的任务状态当成有效决策。

## 分诊与写回

新工作按治理文档中的 `blocker`、`required_follow_up`、`independent` 或 `research` 分类。只有需要跨主线持续管理时才登记 ledger；普通工作由对话或已有 issue/PR 追踪。

对 ledger 管理的工作，仅回写已验证结果、确认的 blocker 或有效作者决策，并保留必要证据链接。分析型请求不改状态；主线切换、`accepted`/`promoted` 和历史记录保留遵循治理文档，不因新请求自动替换主线或授权生产。
