# Models Homomorphic Architecture Design (Post Plan-B)

## Context

在 Plan-B（dimension-agnostic + schema-driven）完成之前，模型层同构化容易反复返工；
在 Plan-B 完成之后，状态契约与结果契约稳定，适合进入结构治理。

因此本设计将“同构化”定义为 Plan-B 的后续工程：

- Gate A: Plan-B 收口（无固定维度主路径 + 旧架构移除）
- Gate B: 7 模型目录/API/测试骨架同构
- Gate C: Program Final（全模型回归 + docs/governance + 无 compat 残留）

## Goals

1. 在不强行同构物理内核的前提下，实现 7 模型的边界同构：接口、目录骨架、测试门禁。
2. 把“最小公共抽象”固化为唯一标准，减少模型扩展时的心智负担与散点改动。
3. 完成从“可维护”到“可预测”的代码组织升级（文件布局、命名、验收矩阵统一）。
4. 最终状态不保留旧兼容路径，避免双轨回潮。

## Non-Goals

1. 不统一模型内部数学表达或求解策略实现细节。
2. 不将 Rotation/GasLiquid 等强行改写为 PNJL 家族形态。
3. 不在本期引入新的物理模型类别。

## Implemented Model Inventory (Scope)

当前纳入同构化范围的 7 个已注册模型：

- `NJL`
- `NJL2`
- `PNJL`
- `PNJLMagnetic`
- `RPNJL`
- `Rotation`
- `GasLiquid`

## Design Principles

### P1. Interface Isomorphic, Kernel Heterogeneous

- 必须同构：`solve_gap / gap_residual / state_schema / result_contract`。
- 允许异构：模型内部方程、状态语义、数值策略。

### P2. Skeleton First

- 每个模型目录遵循统一骨架：
  - `core/`（物理内核）
  - `<ModelName>Model.jl`（适配器）
  - `workflows/`（模型工作流薄适配）
  - tests（unit/integration/regression 对齐模板）

### P3. Minimum Common Abstraction

- 仅抽象可跨 7 模型稳定复用的能力。
- 禁止为追求“看起来统一”而构建过宽基类与 flag-heavy 函数。

### P4. Drift Prevention by Gates

- Gate B 与 Gate C 采用硬门禁，防止执行过程中“只做一半同构”。

## Target End State

1. 7 模型均接入同一 schema/result 契约。
2. 7 模型目录骨架可一眼映射。
3. 测试模板统一：每模型至少具备
   - 接口契约 unit
   - 主链路 integration
   - 基线 regression
4. 文档只描述新主线；不再将 compat 作为常规路径。

## Dependency and Sequencing

- 本 spec 依赖 Plan-B spec 达成 Gate A 后执行。
- 执行顺序：A -> B -> C，禁止 A/B 并行硬推进。

## Acceptance Criteria

1. `Models.registered_model_kinds()` 中 7 模型全部通过同构验收矩阵。
2. 目录骨架差异仅出现在物理内核细节，不出现在接口层结构。
3. 不存在旧架构回退入口与兼容桥接残留。
4. 单元/集成/回归 + docs/governance 全绿。

## Risks and Mitigations

1. Risk: 过度抽象造成可读性下降。
   - Mitigation: 采用“禁止抽象清单”（见计划）并做 PR 评审门禁。
2. Risk: 同构化影响行为稳定。
   - Mitigation: 每模型保留固定回归点与接口一致性快测。
3. Risk: Gate A 未完全稳定就进入 Gate B。
   - Mitigation: 在计划中显式设置 A 的前置检查命令。
