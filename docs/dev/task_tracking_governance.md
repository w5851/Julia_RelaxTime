# 跨主线任务状态与分诊治理

`config/governance/task_tracks.toml` 保存主线、当前任务、依赖与恢复所需 evidence。明确范围的普通工作使用对话计划或已有 issue/PR，不要求新增 track、item、任务文档或执行日志。

## 计划与持久记录

- 对话计划承载近期步骤；稳定设计、公式、验收与作者决策放入对应文档或产物证据。不要把同一份进度复制到计划、任务单、ledger 和日志。
- 现有 schema 对 `ready/active/blocked/review/accepted/promoted` 的 item 仍要求存在的 `task_file`。它可以引用已有设计或简短任务合同，不能用一次对话替代必需的文件引用。
- 登记了 ledger 的工作，在其中维护当前主线、代码/证据 SHA、依赖和下一步；只保存影响后续执行的状态与必要证据链接，不粘贴完整命令输出、实验流水或重复 DoD。
- 未冻结的任务说明维护当前范围、剩余工作、验收条件和直接证据，原位修订过期摘要。正式实验记录、作者决策和已冻结/哈希绑定文档保持原始可追溯性；需要提炼时另设简短入口并链接来源。
- 开发执行记录仅在用户或现有任务合同明确要求时，按指定文件的现有格式直接补充；不自动选择最新任务或创建日志。可按批次标识读取旧记录以查重、补充或关联勘误。日常实现和验证不自动触发追加；科研实验与结论留存遵循[科研证据治理](../analysis/governance/research_evidence.md)。

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

先从用户指定的对话、文件或 issue 恢复目标。范围明确且不涉及主线协调时，检查相关 Git 状态后直接执行。仍有多个候选主线、下一项不明确或需要切换任务时，再校验 ledger 和 Git 状态：

```powershell
julia --project=. scripts/dev/check_task_ledger.jl --preflight
```

已知目标 track 时可加 `--track <id>`。默认报告分支、HEAD、变更路径总数与前 10 项，以及选定 track/current item 的任务文件、阻塞、下一步和最多 5 条 evidence 入口。需要完整路径时加 `--full-paths`；完整 evidence 仍从 ledger 按需读取。

ledger 校验失败时仍报告 Git 上下文并返回失败，不输出可供续接的可信任务选择。根据有效报告只读取选定 track、current item、阻塞依赖和必要证据；不把整个 ledger、全部 active 文档或历史日志都送入上下文。

简述当前任务、阻塞依赖、下一步和本轮范围即可。若旧文档与现状冲突，以当前实现和有效决策为依据，明确剩余差异后再推进。

分析型请求不修改 ledger；对 ledger 管理的任务，验证结果或 blocker 成立后同步事实，不因新需求静默替换 primary track。

## 历史与验证

终态 ID 和仍被依赖的 evidence 保持可解析。普通任务历史可通过 Git/归档文档恢复；只有外部保存或真实恢复需求才引入独立 hash archive 协议，不预设条目数量、行数或冲突次数门槛。文件迁移先检查引用；压缩阅读入口不能破坏原始证据。

本地验证按 ledger/校验器的变化选择：

```powershell
julia --project=. scripts/dev/check_task_ledger.jl
julia --project=. tests/unit/config/test_task_ledger.jl
```

CI 当前为 advisory。schema、引用和状态转换测试保留，可变任务事实不进入 core 固定快照。
