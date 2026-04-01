# Wave-E New Agent Handoff (Temp)

> 用途：给零上下文 new agent 的启动交接。
> 状态：临时文档（Wave-E 启动后可并入任务单或归档）。

## 0. 当前仓库与分支状态

- 基线分支：`main` 已包含 Wave-D 合并结果（PR #45 已合并）。
- Wave-E 工作分支：`feat/pnjl-solver-decoupling-wave-e`
- 旧分支清理：`feat/pnjl-solver-decoupling-wave-d` 远端/本地均已删除。

## 1. Wave-E 执行主线（仅用以下文档）

1. `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-E任务单.md`
2. `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-e-design.md`
3. `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-e-implementation-plan.md`

参考台账：

- `docs/dev/active/2026-04-01_统一扫描入口兼容退场映射表_Wave-E.md`

参数契约参考：

- `docs/architecture/models_param_container_contract.md`

## 2. 已完成前置门（E0）

- PR #45（Wave-D）状态：`MERGED`
- URL：`https://github.com/w5851/Julia_RelaxTime/pull/45`
- merge commit：`ecf77f5529325f18953739324af6953a5411c838`

## 3. 下一步起点（建议）

从 E1 开始，严格 TDD：

1. 先补失败 integration 测试（统一 CLI 契约）
2. 再最小实现 `scripts/models/run_unified_scan.jl`
3. 跑本步定向测试 + smoke + governance
4. 回填任务单证据和勾选

## 4. Wave-E 关键约束（避免跑偏）

- compat 退场走阈值证据，不做无证据强删。
- 统一入口目标模型族：`PNJL/NJL/RPNJL/PNJLMagnetic/Rotation/GasLiquid`。
- `pnjl_aniso` 作为参数化模式（`profile/xi`），不新增独立 `model_kind`。
- 参数输入采用“统一容器 + 语义分层”：`control/background/model/numerics`。

## 5. 常用校验命令

- `julia --project=. scripts/dev/check_docs_consistency.jl`
- `julia --project=. scripts/dev/check_active_docs_governance.jl`
- `julia --project=. scripts/dev/check_unit_skip_policy.jl`

## 6. 交接完成判据（给 new agent）

- 能在零上下文下仅凭本文件 + 三件套定位 Wave-E 当前阶段。
- 能明确 E0 已完成，直接进入 E1。
- 能按任务单约束执行并持续回填证据。
