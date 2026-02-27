---
name: doc-exec-log-append
description: 以 append-only 方式向执行台账追加记录，默认不回读历史正文，降低上下文污染。关键词：执行记录、台账追加、append-only、batch、不回读上下文。
---

# Skill: doc-exec-log-append

## Purpose

将“执行事实”以标准模板追加到执行台账末尾，且默认不读取台账历史内容，减少上下文污染并保持记录可追溯。

本 Skill 只处理“执行记录追加”，不负责开发任务定义、范围变更、DoD 设计。

默认协作行为：当 agent 在推进开发计划并产生可记录执行事实时，应自动触发本 Skill 追加台账；若目标台账不存在，先创建再追加。

---

## When to Use

当用户意图是以下之一时使用：

- “把这次执行追加到台账”
- “记录 Batch N* 执行结果”
- “写执行记录，不要改任务单”
- “按模板补一条执行日志”
- “继续推进开发任务”（隐式触发：有执行事实即追加）

关键词：

- 执行记录 / 执行台账
- append / 追加 / 批次记录
- Batch N1/N2/N3/N4
- 不回读上下文 / 降低上下文污染

---

## When NOT to Use

- 用户要修改“开发任务单”（目标、范围、DoD）时。
- 用户要归档文档（应使用 `doc-archive`）。
- 用户要求重写历史执行记录时（非追加场景）。

---

## Strong Constraints (Must Follow)

1. **Append-only**：只允许在目标执行台账文档末尾追加，不改写历史记录。
2. **No read-back by default**：追加时默认不读取执行台账历史正文。
3. **Self-contained entry**：每条记录必须自包含：目标、代码变更、验证命令、产物、结果、主线映射。
4. **Task/Log separation**：开发任务单与执行台账分离维护，禁止把任务定义写入执行台账。
5. **Auto log on execution progress**：当开发步骤产生可审计执行事实（命令、产物、结果）时，默认自动追加记录，不依赖用户显式要求“写台账”。

例外（仅以下情况允许读取台账）：

- 目标台账文件不存在，需要创建并写入初始头部模板；
- 用户明确要求“检查历史重复批次/冲突”；
- 文档格式损坏，需最小修复后再追加。

---

## Inputs

最小输入字段：

- `log_file`: 目标执行台账路径（通常在 `docs/dev/active/`）
- `batch_id`: 如 `Batch N2`
- `goal`: 本批次目标
- `code_changes`: 代码变更摘要
- `validation_commands`: 验证命令列表
- `artifacts`: 产物路径列表（如 `outputs/results/...`）
- `result`: 执行结果
- `mainline`: 主线映射（N1/N2/N3/N4）

可选字段：

- `timestamp`
- `notes`

---

## Standard Procedure

1. 定位任务文档与目标台账：
	- 优先使用显式输入 `log_file`；
	- 否则根据 `task_file` 自动推导同目录“执行台账”；
	- 若仍缺失，自动发现 `docs/dev/active` 最新“开发任务单”并推导。
2. 若台账不存在，创建标准台账骨架（仅初始化头部，不回读历史）。
3. 校验输入字段是否齐全（缺失则报错并要求补齐）。
4. 组装标准记录块（固定模板）。
5. 直接追加到文档末尾（append-only）。
6. 返回追加摘要（文件、Batch、产物数量、命令数量）。

---

## Output Template (Append Block)

直接复制并填充：

- [ ] 批次号：{{batch_id}}
	- 目标：{{goal}}
	- 代码变更：{{code_changes}}
	- 验证命令：
		- `{{cmd_1}}`
		- `{{cmd_2}}`
	- 产物：
		- `{{artifact_1}}`
		- `{{artifact_2}}`
	- 结果：{{result}}
	- 主线映射（N1/N2/N3/N4）：{{mainline}}

---

## Quality Checklist

- 记录是“追加”而非“覆盖”
- 必填字段完整
- 命令与产物可追溯
- 未混入任务目标/DoD 细节
- 仅在允许例外时读取台账历史

---

## Execution Binding (Implemented)

本仓库已落地脚本：`scripts/dev/append_exec_log.jl`，用于统一执行追加：

- 输入通过 CLI 参数传入；
- 脚本内部不读取台账历史正文；
- 仅做文件存在性/字段合法性校验后追加；
- 失败时给出结构化报错。
- 当未显式提供 `--log-file` 时，自动根据 `--task-file` 或 active 最新任务单推导台账路径；不存在则自动创建。

这将“no read-back”从约定升级为可执行约束，降低人工偏差。

---

## Example Invocation (Conceptual)

```powershell
julia --project=. scripts/dev/append_exec_log.jl \
	--log-file docs/dev/active/2026-02-26_多重派发重构与PNJL迁移下线执行台账.md \
	--batch "Batch N2" \
	--goal "移除旧入口并完成冒烟校验" \
	--code-change "删除 src/pnjl 兼容路径调用点" \
	--cmd "julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD" \
	--cmd "UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl" \
	--artifact "outputs/results/pnjl_prune_wave_snapshot_20260226_101500.txt" \
	--result "通过" \
	--mainline "N2"
```

---

## Collaboration Notes

- 若用户仅说“追加执行记录”，优先使用本 Skill。
- 若用户在推进开发计划但未明确提及台账，只要产生执行事实，也应默认调用本 Skill 追加。
- 若用户同时要求“更新任务单状态”，先追加执行台账，再由相关文档流程更新任务单。