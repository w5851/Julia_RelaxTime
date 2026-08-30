# 跨主线任务状态与分诊治理

## 目的

`config/governance/task_tracks.toml` 是跨主线状态的机读索引，用来恢复当前主线、当前任务、依赖、下一步和已确认 evidence。它不复制 active task 的 Definition of Done，不复制执行台账，也不替代数值产物、manifest 或 claim ledger。

本层解决的是“任务在哪里、为何暂停、下一步是什么”，不是“数值结果是否正确”。具体验收仍以 `docs/dev/active/`、`docs/dev/backlog/`、`docs/analysis/` 和正式产物 registry 为准。

ledger 顶层、track 和 item 字段是固定 schema；validator 会拒绝缺失或未知字段。`current_task`
只索引 item，不复制任务单内容；`current_sha` 锚定不可变代码、计算或 evidence SHA，不表示
ledger 自身提交；`current_branch`、`run_ids` 和 `evidence` 只保存恢复所需的引用与格式化事实。

## 状态语义

```text
inbox -> triaged -> ready -> active -> review -> accepted -> promoted -> archived
                         |         |          |
                         v         v          v
                      deferred  blocked    cancelled
```

- `accepted` 表示当前 DoD 或作者审核已通过；不等于允许进入下一条生产链。
- `promoted` 表示对应 promotion gate 已通过，并允许使用 `unlocks` 继续推进。
- `accepted -> archived` 只适用于 `promotion_required=false`；需要晋升的任务必须先留下 promotion/gate evidence。
- `archived` 表示主线与后续依赖均已收口。
- `blocked` 必须有可解析的 `blocked_by`；`deferred` 必须有父任务、原因、下一步和 backlog 文件。
- `full_hybrid_candidate`、`density_feasible_candidate`、`C2 limited feasibility` 等是 evidence/verdict 标签，不是 task ledger 状态。

状态转换合同和字段校验由 `scripts/dev/check_task_ledger.jl` 执行。首版没有自动写入器；状态只有在实现事实、验证结果或明确作者决策成立后才更新。

## 四类分诊

新发现任务必须归入且只能归入一个分类：

| 分类 | 处理规则 |
| --- | --- |
| `blocker` | 可以暂时切换当前工作，但先记录原任务、依赖和恢复路径。 |
| `required_follow_up` | 挂到父任务下，当前任务完成后仍必须处理，不能隐式消失。 |
| `independent` | 建议独立 track、branch 和 PR，不污染数值主线。 |
| `research` | 进入 backlog 或 `deferred`，不自动开始实现。 |

`current_task` 指向 `[[items]]` 的 item ID；每个 item 只能属于一个 `track_id`。文件路径、run、SHA、review 和 promotion 证据使用带前缀的 evidence 引用，例如 `file:...`、`run:...`、`review:...`、`promotion:...`。

## Harness 恢复协议

当请求涉及多个主线、主线不明确或出现新 blocker/follow-up 时，先运行：

```powershell
julia --project=. scripts/dev/check_task_ledger.jl
julia --project=. scripts/dev/check_task_ledger.jl --preflight
```

然后按顺序输出：

```text
Primary track:
Current task:
Status:
Blocked by:
Next action:
This-turn scope:
Non-goals:
Deferred/research items:
```

如果用户已经明确指定一个 active task 文件，继续请求仍由 `doc-implementation` 直接处理；harness 只负责确认是否存在跨主线冲突。`doc-exec-log-append` 仍要求用户或任务文档显式授权。

分析型请求不修改 ledger。实现型请求在验证结果或 blocker 已确认后更新 ledger，并在结束时报告恢复入口；作者审核、昂贵 Actions、reference promotion 和 failed Actions 自动分诊均不由本层代办。

## 当前三条主线

- `issue130-phase` 已完成 `full_hybrid_candidate` 作者审核、accepted reference promotion 和 legacy retirement；保留不可变 evidence，后续仅做治理归档。
- `rs-transport` 已完成 accepted-primary p104/p128 的 matched numerical convergence 并归档；旧 strict-era raw/figure 保持不变，accepted evidence 不替换正式结果。
- `plot-sop` 的公共 SOP 合同已登记为 `promoted`；后续逐图族迁移仍按独立 task/PR 推进，不由该 promotion 记录吸收历史图像或批量迁移工作。

## 验证与演进

本地检查：

```powershell
julia --project=. scripts/dev/check_task_ledger.jl
julia --project=. tests/unit/config/test_task_ledger.jl
```

首版 CI 为 advisory，仅在 ledger、harness、治理脚本、测试或 active/backlog 文档变化时运行。M4 才评估 PR 模板 `track_id`、长期未更新检查和 failed Actions 自动生成 blocker。
