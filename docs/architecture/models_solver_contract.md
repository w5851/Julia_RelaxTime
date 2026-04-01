# Models 求解接口契约（v0.1）

更新时间：2026-03-31

本文档定义 `Models` 统一入口下的求解接口规范。目标是：

- 在不重写现有数值内核的前提下，统一调用语义；
- 约束公开 API 的边界，降低调用方与具体后端的耦合；
- 为后续可选 SciML adapter 预留兼容空间。

## 1. 定位与范围

本规范适用于：

- `src/models/Models.jl`
- `src/models/entrypoints.jl`
- `src/models/solver/`
- `src/models/workflows/`

本规范不要求：

- 将当前求解内核全面迁移到 SciML 依赖栈；
- 替换现有 `NLsolve` 等后端实现。

## 2. 术语

- Problem-like：表示“问题定义”的类型或结构（例如 `FixedMu`、`FixedRho`）。
- Solve 入口：调用方面向任务闭环的统一执行函数（例如 `solve`、`solve_gap_and_transport`）。
- 后端：实际执行数值迭代/求根的实现（例如 `NLsolve.nlsolve`）。
- 结构化结果：面向调用方可直接消费的 `NamedTuple` 或等价稳定结构。

## 3. 强制规则（MUST）

### R1 统一入口规则

新增公开求解能力必须通过 `Models` 导出，并由 `src/models/entrypoints.jl` 或 `Models` 主表面统一暴露。

- 禁止新增“仅内部可见但被外部脚本直接依赖”的临时入口。
- 迁移期兼容入口允许保留，但不得作为新能力默认入口。

### R2 Problem-like 规则

每类约束求解任务必须有明确 Problem-like 契约（类型或结构），并与求解入口解耦。

- 推荐沿用 `ConstraintModes` 体系（如 `FixedMu`、`FixedRho`）。
- 新增模式时，必须在文档中声明输入域、单位与边界条件。

### R3 Solve 单入口语义规则

对外应保持统一的 solve 语义：

- 约束求解走 `solve(...)` / `solve_constraint(...)`；
- 流程编排走 `solve_*workflow(...)`（例如 `solve_gap_and_transport`）。

禁止把 workflow 内部辅助函数直接升级为公开主入口。

### R4 结果契约规则

公开求解入口必须返回结构化结果，并满足最小可诊断要求。

- 约束求解结果应包含收敛与残差信息（如 `converged`、`residual_norm`）。
- workflow 结果应包含复现实验元数据（如 profile、配置路径或等价信息）。
- 允许按场景扩展字段，但禁止破坏既有稳定字段语义。

### R5 后端解耦规则

调用方不得依赖具体后端实现细节。

- 允许内部使用 `NLsolve` 或其它后端；
- 后端特有参数不得直接外泄为跨模块硬依赖；
- 需要后端特化时，必须通过适配层或策略参数承接。

### R6 错误与告警规则

- 参数契约违规应尽早抛出 `ArgumentError`，并明确指出违规参数名；
- 数值失败应返回可诊断信息，必要时使用 `@warn` 给出可复现上下文；
- 禁止静默吞掉关键物理/数值异常。

### R7 单位与命名规则

输入输出命名必须保持单位语义显式：

- fm^-1 量：`*_fm` 或 `*_inv_fm`（按本仓库既有约定）；
- MeV 量：`*_MeV`；
- 高阶单位命名遵循现有模式（如 `sigma_fm4`、`coupling_inv_fm4`）。

### R8 测试与文档规则

新增公开入口必须同时满足：

- 至少一个 unit 测试覆盖主路径；
- 若跨模块链路有行为变化，补充 integration 或 regression 证据；
- 稳定公开入口需同步更新 `docs/api/`。

## 4. 当前实现映射

当前仓库已具备本规范基础：

- 统一入口：`src/models/Models.jl`、`src/models/entrypoints.jl`
- Problem-like：`src/models/solver/ConstraintModes.jl`
- Solve 路径：`src/models/solver/Solver.jl`
- Workflow 路径：`src/models/workflows/TransportWorkflow.jl`、`src/models/workflows/MesonMassWorkflow.jl`
- 迁移守卫：`scripts/dev/check_pnjl_migration_guard.jl`

## 5. 与 SciML 的关系

本规范采用 SciML-like 的“Problem + solve”语义，但不绑定 SciML 具体实现。

- 当前阶段：语义对齐（低风险）
- 后续阶段：可选 adapter 对接（按需求引入）

除非出现跨 ODE/PDE/Optimization/Sensitivity 的明确统一需求，不建议为“风格一致”而重写现有内核。

## 6. 变更流程

涉及本契约的变更，必须按以下顺序推进：

1. 先更新 ADR（`docs/decisions/`）。
2. 再更新本规范文档与 `docs/api/`。
3. 最后更新门禁脚本与测试。

## 7. 最小验收清单

- [ ] 新增入口已在 `Models` 导出并通过统一入口暴露。
- [ ] Problem-like 契约与单位声明完整。
- [ ] 返回结果具备收敛/诊断字段。
- [ ] 未新增后端实现耦合点。
- [ ] 已补充 unit 测试与必要 API 文档。

## 8. Wave-A 兼容层治理补充（2026-03-31）

- 主入口固定为 `Models.solve_constraint`。
- `solve_fixed*` 入口保留为兼容层，不再作为新调用首选。
- 兼容层映射与删除门槛由活动台账维护：
  - `docs/dev/active/2026-03-31_求解器兼容层迁移映射表.md`
