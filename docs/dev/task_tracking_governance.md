# 跨主线任务状态与分诊治理

`config/governance/task_tracks.toml` 保存主线、当前任务、依赖与恢复所需 evidence。明确范围的普通修复可直接由 issue/PR 追踪，不要求新增 track、item 或执行日志。

任务 DoD 留在任务单/issue，数值结论留在分析和产物证据中。当前主线、分支、SHA 和下一步只在 ledger 维护；本页不重复列动态摘要，测试也不固定这些事实。

## 状态语义

```text
inbox -> triaged -> ready -> active -> review -> accepted -> promoted -> archived
                         |         |          |
                         v         v          v
                      deferred  blocked    cancelled
```

- `accepted` 表示该项验收通过，不自动授权下一条生产链。
- `promoted` 表示适用的 promotion gate 已通过。
- `accepted -> archived` 适用于 `promotion_required=false`；需要晋升的工作保留 gate evidence。
- `blocked` 保留可解析依赖；`deferred` 记录父任务、原因、下一步与 backlog。
- 数值 verdict 标签与任务状态分别记录。

字段与转换由 `scripts/dev/check_task_ledger.jl` 校验；真实实现、验证或作者决策成立时更新。`current_sha` 锚定代码/证据来源，允许主线按事实变化，不表示 ledger 文件自身的提交。

## 跨主线分诊

| 分类 | 用途 |
| --- | --- |
| `blocker` | 影响当前推进，记录依赖和恢复路径 |
| `required_follow_up` | 挂到父任务，保留后续处理责任 |
| `independent` | 与当前数值主线分开处理 |
| `research` | 进入 backlog/明确的研究范围 |

只有需要管理跨主线关系时才登记这些分类；不把每个小修或新 PR 都转成 ledger 状态机任务。

## 恢复入口

请求含多个主线、下一项不明确或需要切换任务时，读取 ledger 和相关任务，并运行：

```powershell
julia --project=. scripts/dev/check_task_ledger.jl --preflight
```

简述当前任务、阻塞依赖、下一步和本轮范围即可。用户已给出明确任务文件或 issue 时直接执行，先检查其相关 Git 状态，不重复完整 harness。

分析型请求不修改 ledger；对 ledger 管理的任务，验证结果或 blocker 成立后同步事实。执行日志仅在用户或任务明确要求时追加。

## 历史与验证

终态 ID 和仍被依赖的 evidence 保持可解析。普通任务历史可通过 Git/归档文档恢复；只有外部保存或真实恢复需求才引入独立 hash archive 协议，不预设条目数量、行数或冲突次数门槛。

本地验证按 ledger/校验器的变化选择：

```powershell
julia --project=. scripts/dev/check_task_ledger.jl
julia --project=. tests/unit/config/test_task_ledger.jl
```

CI 当前为 advisory。schema、引用和状态转换测试保留，可变任务事实不进入 core 固定快照。