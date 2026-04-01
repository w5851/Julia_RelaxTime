# PNJL Solver Decoupling Wave-E Design

## Context

Wave-A/B/C/D 已完成契约冻结、治理收敛、扫描/工作流并轨与 hard-deprecate-first。

当前剩余问题不是“主路径是否存在”，而是“兼容层文件是否仍长期共存”。

Wave-E 的目标是完成兼容退场：

1. 仅保留单一路径（新治理层 + unified model-driven entrypoints）。
2. 去除 `scripts/pnjl` 的分散扫描脚本作为业务入口角色。
3. 把扫描参数化能力扩展到当前实现框架下全部模型族：`PNJL/NJL/RPNJL/PNJLMagnetic/Rotation/GasLiquid`，并明确 `pnjl_aniso` 作为参数化模式而非独立模型种类。

## Goals

1. 移除（或归档）已 hard-deprecated 的扫描/工作流 compat adapter 路径，避免双轨长期并存。
2. 建立单一“通用扫描脚本入口”（治理层主线），统一承载 T-mu / T-rho / phase-workflow 调度。
3. 在统一脚本与 `Models` 入口之间建立严格、可查询、可回归的模型能力矩阵。
4. 完成 PR #45 的收尾（review/checks/merge/文档状态一致性）并作为 Wave-E 启动前置门。

## Non-Goals

1. 不重写底层求解器数值内核。
2. 不在 Wave-E 引入新的物理模型类别（仅覆盖当前已实现模型）。
3. 不做与入口治理无关的大规模重构。

## Design Decisions

### E1: compat 退场采用“证据阈值 + 一次性切换”

只在满足以下阈值时执行删除/归档：

- PR #45 合并后主线 green；
- 新统一扫描脚本对旧入口具备命令级等价覆盖（代表点 + smoke + regression）；
- 文档与迁移映射完成 closed 状态回填。

### E2: 单一脚本入口采用“命令子域 + 模型能力矩阵”

新脚本建议采用统一入口（示例：`scripts/models/run_unified_scan.jl`），通过子命令区分任务域：

- `scan tmu`
- `scan trho`
- `workflow phase`

并以能力矩阵限制模型-模式组合，避免“参数合法但模型不支持”的隐性失败。

### E3: `pnjl_aniso` 定位为参数化 profile，不作为独立 model_kind

当前架构下各向异性属于 `xi`/配置维度，映射为 `PNJL` 或 `NJL` 等模型的参数化求解路径。
Wave-E 需在文档与 CLI 上显式这一点，防止把 `pnjl_aniso` 再次做成平行脚本分叉。

### E4: 新治理层文件最小集合原则

目标状态只保留：

- `src/models/entrypoints.jl`（统一入口）
- `src/models/scans/*`（统一扫描实现）
- `src/models/solver/*`（治理与求解契约）
- 单一脚本入口（CLI 薄适配）

其余 legacy 脚本转为删除或归档说明，不再作为生产入口。

### E5: 参数输入契约采用“统一容器 + 语义分层”

Wave-E 参数契约采用单一 NamedTuple 容器承载全部输入，但按语义分层：

- `control`（`mode/T/mu/rho/...`）
- `background`（`xi/eB/omega/profile/...`）
- `model`（`kind/profile`）
- `numerics`（网格与容差）

其中 `pnjl_aniso` 作为 `PNJL + profile/xi` 的参数化表达，不新增独立 `model_kind`。
该约束详见：`docs/architecture/models_param_container_contract.md`。

## Target Artifacts

1. Wave-E active 任务单（含 PR #45 收尾前置门）。
2. Wave-E 实施计划（TDD-first，E1->E2->E3->E4）。
3. 通用扫描脚本与能力矩阵文档。
4. compat 退场迁移表（最终 `removed/archived` 状态）。

## Verification Strategy

Wave-E 最小验证矩阵：

1. 统一脚本命令级定向测试（含模型矩阵校验）；
2. 旧脚本退场后 parity/regression 稳定性验证；
3. unit/integration/regression smoke；
4. docs governance + migration consistency；
5. PR #45 merge 后状态复核。

## Risks and Mitigations

1. Risk: 外部调用仍依赖旧脚本路径。
   - Mitigation: 先发布迁移指引与错误提示窗口，再执行 remove/archive。
2. Risk: 通用脚本过度泛化导致参数复杂度上升。
   - Mitigation: 子命令拆分 + 模型能力矩阵 + 显式参数校验。
3. Risk: 模型覆盖声明与真实实现不一致。
   - Mitigation: 用能力矩阵自动测试（model_kind × mode）做门禁。
