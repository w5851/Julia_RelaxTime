---
name: codex-task-harness
description: 将宽泛或跨主线的 Julia_RelaxTime 需求收束为可执行 harness，并从 task ledger 恢复主线、依赖、下一步和四类新任务分诊。用于多个 active track、主线不明确、继续/下一步请求出现跨线歧义，或发现 blocker、required follow-up、independent、research 项；明确单一 active task 文件时交给 doc-implementation。
---

# Codex Task Harness

## Purpose

把“大而散”的请求压缩成可以立即执行的边界与验证框架，避免计划本身成为交付物。

跨主线状态的唯一机读来源是 `config/governance/task_tracks.toml`，校验入口是
`scripts/dev/check_task_ledger.jl`，状态职责和恢复输出格式见
`docs/dev/task_tracking_governance.md`。

## Hard rules

- 先写不允许改变的范围和契约，再写允许修改的内容。
- 只保留会影响执行的风险和假设。
- 给出最小必要验证，不以“以后补测试”代替验收。
- 只路由最少但足够的 skill；边界清楚后立即执行。
- 分析型请求不修改 task ledger；实现结果或已确认 blocker 才能更新它。
- `accepted` 不等于 `promoted`；作者审核、昂贵 Actions 和 promotion gate 不自动代办。

## Harness contract

- 目标：本轮必须完成的结果
- 允许修改：文件、模块、数据或文档范围
- 不要修改：显式非目标和外部影响边界
- 必须保持：稳定入口、求解语义、单位或数据契约
- 风险：最可能失败的 1–3 点
- 验证：测试层、治理检查、结果核对和文档义务
- 路由：本轮实际需要的其他 skill

## Repository routing

- 数值 baseline 制度：`baseline-regression-governance`
- transport/relaxtime 漂移：`transport-regression-keeper`
- API 文档：`api-doc-authoring`
- 前端数据契约：`frontend-contract-author`
- 论文证据写作：`paper-traceable-coauthor`
- 技术文档：`doc-coauthoring`
- 实验日志：`experiment-logbook-append`
- 执行台账：`doc-exec-log-append`
- Julia 实现：`julia-pro`

## Routing boundary

- 明确指定一个 active task 文件的“继续”“下一步”或“按计划推进”，直接路由
  `doc-implementation`，不重复触发本 harness。
- 存在多个 active track、primary 不明确、请求引入新 blocker/follow-up，或需要在多个待办中选下一项时，先执行本 harness，再把选定任务交给最小后续 skill。
- 新增 repair、PR 或 workflow 但没有明确主线时，先登记分诊类别和父任务，不自动创建 PR。

## Start protocol

1. 运行 `julia --project=. scripts/dev/check_task_ledger.jl`。
2. 读取 `primary_track`、对应 `current_task`、相关 active task 和必要 evidence。
3. 运行 `julia --project=. scripts/dev/check_task_ledger.jl --preflight`，记录 branch、HEAD 和 dirty paths；不覆盖用户改动。
4. 输出当前主线、当前任务、状态、阻塞依赖、下一步、本轮 scope、non-goals 和 deferred/research items。

## Triage contract

新发现任务只能归入一个分类：

- `blocker`：可暂时切换，但先保存原任务和恢复路径。
- `required_follow_up`：挂到父任务，不能因当前任务结束而消失。
- `independent`：使用独立 track、branch 和 PR。
- `research`：放入 backlog 或 `deferred`，不自动实现。

使用 `inbox -> triaged -> ready -> active -> blocked|review -> accepted -> promoted -> archived`
状态合同。`blocked` 必须有 `blocked_by`；`accepted` 必须有 evidence；`promoted` 必须有
promotion gate 和所需作者审核；`deferred` 必须有父任务、原因、下一步和 backlog 位置。

## End protocol

- 验证通过后才把确认的结果、evidence 或 blocker 写回 ledger。
- 若实现中途因 blocker 切换主线，先记录 blocker 和原任务，再切换；恢复时重新读取 primary/current task。
- 不把候选 verdict 写成状态，不把 feasibility 写成 convergence，不因新任务静默关闭旧任务。
- `doc-exec-log-append` 仍要求用户或活动任务显式授权；ledger 不复制执行日志和 DoD。

1. 用一句话锁定主目标。
2. 冻结文件、语义和外部影响边界。
3. 选择验证下限和文档义务。
4. 选择最小 skill 组合。
5. 执行并在范围或风险实质变化时重新锁定边界。
