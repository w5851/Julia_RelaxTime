# PNJL求解器解耦框架兼容性评审与 Wave-A 实现计划

> 执行主线说明：本文档是当前唯一执行主线（任务勾选、进度、验收、证据回填）。
> 设计讨论与外部审阅依据保留在：
> `docs/dev/active/2026-03-31_PNJL求解器解耦框架约束与分层草案.md`。
> 架构决策与 superpowers 文档：
> `docs/decisions/0003-pnjl-solver-decoupling-governance-and-ad-implicit-contract.md`
> `docs/superpowers/specs/2026-03-31-pnjl-solver-decoupling-wave-a-design.md`
> `docs/superpowers/plans/2026-03-31-pnjl-solver-decoupling-wave-a-implementation-plan.md`

## 1. 背景与目标

- 对草案 `docs/dev/active/2026-03-31_PNJL求解器解耦框架约束与分层草案.md` 进行与当前 `src/models` 主线实现的兼容性核验。
- 重点确认“隐式微分 + 自动微分（AD/ID）”链路是否可在解耦框架下稳定延续。
- 输出可执行的 Wave-A 任务计划（含文件落点、验收标准、验证命令、回退路径）。

## 2. 范围与非范围

### 2.1 本期范围（Wave-A）

- [ ] 固化求解问题契约（ProblemSpec）与显式参数接口规范。
- [ ] 固化 AD/ID 契约（Dual-safe 条件函数、统一 forward/conditions 形态）。
- [ ] 固化候选治理接口（硬约束规则 + 选择器 + context）。
- [ ] 固化兼容层迁移映射（旧入口 -> 新入口 -> 删除门槛）。
- [ ] 增补最小验证矩阵（Jacobian smoke、维度一致性、选择稳定性、配置隔离）。

### 2.2 非范围（Wave-A 不做）

- [ ] 不进行一次性大重构或跨模块重写。
- [ ] 不替换当前底层求根后端（NLsolve 继续保留）。
- [ ] 不在 Wave-A 直接删除兼容层文件（只设删除门槛与迁移台账）。

## 3. 现状兼容性评审结论

### 3.1 总体结论

- [x] 草案方向与主线总体兼容：现有入口已具备统一化基础，能够承接 L0-L4 分层收敛。
- [x] AD/ID 主链可用：隐函数求解器与导数模块已经在主线运行，草案提出的“显式参数接口 + mode 可识别”方向正确。
- [ ] 仍存在结构性风险：双轨求解、全局可变配置、维度硬编码、候选治理 context 缺失，若不在 Wave-A 锁定会造成后续漂移。

### 3.2 证据映射（关键代码点）

- 统一入口：`src/models/entrypoints.jl`、`src/models/solver/Solver.jl`
- 约束模式体系：`src/models/solver/ConstraintModes.jl`
- 方程构造与 residual builder：`src/models/solver/Conditions.jl`
- 隐函数求解与导数入口：`src/models/implicit_gap.jl`、`src/models/solver/ImplicitSolver.jl`
- 导数主线：`src/models/derivatives/ThermoDerivatives.jl`、`src/models/derivatives/ConservedChargeSusceptibilities.jl`
- 约束与候选治理（另一条主链）：`src/models/constraint_solver.jl`
- 扫描调用收敛入口：`src/models/scans/TmuScan.jl`、`src/models/scans/TrhoScan.jl`

## 4. 兼容性风险清单（Wave-A 必须收口）

### R1 双轨求解职责重叠

- 现象：`constraint_solver.jl` 与 `solver/ImplicitSolver.jl` 同时承担约束求解、候选治理与择优逻辑。
- 风险：规则口径和失败语义分叉，回归难以判定“谁是基线”。
- Wave-A 要求：确定主从关系，冻结“唯一治理口径”接口。

### R2 显式参数接口不够硬

- 现象：大量通过 `build_residual!` 返回闭包，跨模式可复用但显式契约不足。
- 风险：AD/ID 接口难做一致性约束，后续 mode 扩展时容易反复改签名。
- Wave-A 要求：残差必须有“显式参数版本”，闭包仅做适配层。

### R3 全局可变配置导致可重入风险

- 现象：`IMPLICIT_CONFIG`、`CURRENT_XI/P_NUM/T_NUM` 使用全局 `Ref`。
- 风险：多扫描链路交错时出现状态污染，导致结果不可复现。
- Wave-A 要求：把配置收敛到 solver 实例或局部上下文。

### R4 维度与类型硬编码限制扩展

- 现象：多处写死 `SVector{5}` 和 8 维拼装。
- 风险：新增 mode 或新模型状态维度时需要跨文件改动，易出错。
- Wave-A 要求：在问题契约中显式 `x_dim/theta_dim`，并提供统一解包规则。

### R5 候选治理缺少 context 契约

- 现象：已有 `default_hard_constraint_rules` 与 `select_pressure_max_candidate`，但上下文（连续性优先、相变邻域、策略原因）未结构化。
- 风险：扫描/单点/回归的择优语义不一致。
- Wave-A 要求：规则与选择器接口显式接收 `context`。

## 5. Wave-A 契约设计（Spec）

## Chunk 1: 求解问题契约（ProblemSpec）

### 5.1 ProblemSpec 最小字段

```julia
ProblemSpec(
    mode,                 # ConstraintMode 或等价问题类型
    x_dim::Int,           # 状态维度
    theta_dim::Int,       # 参数维度
    residual!,            # 显式参数接口 residual!(F, x, theta, cfg, mode)
    forward_solve,        # forward_solve(theta, cfg, mode) -> (x, meta)
    conditions,           # conditions(theta, x, meta, cfg, mode) -> F
    unpack_solution,      # solution -> (x_state, mu_vec, ...)
    postprocess,          # -> omega, pressure, rho_norm, entropy, energy, masses
    hard_rules,           # Vector{HardRule}
    selector,             # CandidateSelector
)
```

### 5.2 约束

- `residual!` 必须保留显式参数版本；允许另外提供 closure adapter。
- `conditions` 路径必须 Dual-safe，不得做 `Float64` 强制转换。
- `x_dim/theta_dim` 与 `ConstraintMode` 一致，且由测试覆盖。

## Chunk 2: AD/ID 契约

### 5.3 统一接口

- `forward_solve(theta, cfg, mode) -> (x, meta)`
- `conditions(theta, x, meta, cfg, mode) -> F`

### 5.4 Dual-safe 清单

- 禁止在 `conditions` 对 `theta`/`x` 做 `Float64(...)`。
- 允许在 `forward_solve`（primal）做数值类型收敛。
- `gap_residual`/`model_pressure` 入口类型要保持可传播 Dual。

## Chunk 3: 候选治理契约

### 5.5 规则接口

```julia
hard_rule(candidate, params, context) -> (ok::Bool, reason::Symbol)
selector(candidates, params, context) -> selected_index::Int
```

### 5.6 context 最小字段

- `mode`
- `continuity_seed_available::Bool`
- `phase_hint::Symbol`
- `residual_norm_max::Float64`
- `prefer_continuity::Bool`

## Chunk 4: 兼容层迁移契约

### 5.7 迁移台账字段

- 旧入口（文件 + 函数）
- 新入口（文件 + 函数）
- 当前状态（active/deprecated/removable）
- 删除门槛（测试、文档、调用方清零）
- 目标波次（B/C/D）

## 6. Wave-A 执行任务单（Plan）

### Task A1：冻结 ProblemSpec 与显式参数残差契约

**目标**

- 把“模式 + 维度 + residual + postprocess + 治理策略”固化为统一结构，避免重复分支。

**预计改动文件**

- 新增：`src/models/solver/ProblemSpec.jl`
- 修改：`src/models/solver/Conditions.jl`
- 修改：`src/models/solver/ImplicitSolver.jl`
- 修改：`src/models/solver/Solver.jl`
- 修改：`src/models/Models.jl`（导出与 include 接线）

**可勾选步骤**

- [x] 定义 `ProblemSpec` 结构与构造器。
- [x] 为 `FixedMu/FixedRho/FixedAsymmetricRho/FixedEntropy/FixedSigma` 提供 spec builder。
- [x] 在 `Conditions` 增加显式参数 residual 入口。
- [x] 保留原 closure 构造函数并改为适配层。

**验收标准**

- [x] 各 mode 均可通过 spec 构建求解问题。
- [x] 不破坏现有 `Models.solve*` 对外语义。

### Task A2：固化 AD/ID 契约与 Dual-safe 边界

**目标**

- 统一 `forward_solve/conditions` 语义，显式标识可 Float64 路径与禁止路径。

**预计改动文件**

- 修改：`src/models/implicit_gap.jl`
- 修改：`src/models/solver/ImplicitSolver.jl`
- 修改：`src/models/derivatives/ThermoDerivatives.jl`
- 修改：`src/models/derivatives/ConservedChargeSusceptibilities.jl`
- 新增：`docs/architecture/models_ad_implicit_contract.md`

**可勾选步骤**

- [x] 抽取统一 `forward_solve` 适配。
- [x] 抽取统一 `conditions` 适配。
- [x] 补充 Dual-safe 清单文档与代码注释锚点。
- [x] 补充 smoke 测试覆盖 dual 传播。

**验收标准**

- [x] FixedMu 与 flavor-mu Jacobian smoke 通过。
- [x] `conditions` 路径无强制 `Float64` 回退。

### Task A3：配置状态治理（移除全局 Ref 热点）

**目标**

- 降低全局状态污染风险，保障并行/交错调用可复现。

**预计改动文件**

- 修改：`src/models/solver/ImplicitSolver.jl`
- 修改：`src/models/derivatives/ThermoDerivatives.jl`
- 可能修改：`src/models/derivatives/ConservedChargeSusceptibilities.jl`

**可勾选步骤**

- [x] 定义不可变 config 结构并绑定到 solver 实例。
- [x] 将全局配置 setter 改为实例化传参。
- [x] 为扫描路径补充“前后调用不串扰”测试。

**验收标准**

- [x] 连续两次不同配置调用结果互不污染。
- [x] 既有导数 API 兼容保留。

### Task A4：统一候选治理接口（规则 + 选择）

**目标**

- 把硬约束与择优从“函数散点”提升为可复用策略接口。

**预计改动文件**

- 修改：`src/models/constraint_solver.jl`
- 修改：`src/models/solver/ImplicitSolver.jl`
- 新增：`src/models/solver/CandidateGovernance.jl`
- 修改：`src/models/scans/TmuScan.jl`
- 修改：`src/models/scans/TrhoScan.jl`

**可勾选步骤**

- [x] 定义 `HardRule`/`Selector` 标准签名。
- [x] 接入 context（连续性、phase_hint、模式信息）。
- [x] 保留 `select_pressure_max_candidate` 作为默认策略实现。

**验收标准**

- [x] 单点与扫描路径采用同一策略接口。
- [x] 选择稳定性回归测试通过。

### Task A5：双轨收敛与兼容层治理台账

**目标**

- 明确 `solver/*` 主干地位，`constraint_solver.jl` 只保留兼容职责。

**预计改动文件**

- 新增：`docs/dev/active/2026-03-31_求解器兼容层迁移映射表.md`
- 修改：`src/models/constraint_solver.jl`
- 修改：`src/models/solver/Solver.jl`
- 修改：`docs/architecture/models_solver_contract.md`

**可勾选步骤**

- [x] 产出旧->新接口映射与删除门槛。
- [x] 对兼容入口加迁移注记（非破坏）。
- [x] 在文档中声明目标删除波次。

**验收标准**

- [x] 不新增兼容散点。
- [x] 调用方可追溯到主干入口。

### Task A6：最小验证矩阵落地

**目标**

- 用最小但有效的测试集防止契约漂移。

**预计改动文件**

- 新增：`tests/unit/models/test_problem_spec_contract.jl`
- 新增：`tests/integration/models/test_ad_implicit_contract_smoke.jl`
- 新增：`tests/regression/pnjl/test_candidate_governance_stability.jl`
- 新增：`tests/integration/models/test_solver_config_isolation_smoke.jl`

**可勾选步骤**

- [ ] 增加多 mode 维度一致性测试。
- [ ] 增加 FixedMu + flavor-mu Jacobian smoke。
- [ ] 增加候选治理稳定性回归。
- [ ] 增加配置隔离 smoke。

**验收标准**

- [ ] Wave-A 新增测试全部通过。
- [ ] 现有关键回归不退化。

## 7. 验证命令（建议）

- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_active_docs_governance.jl`

## 8. 里程碑

- [ ] M1：ProblemSpec + AD/ID 契约文档冻结（A1 + A2）。
- [ ] M2：配置治理 + 候选治理接口收敛（A3 + A4）。
- [ ] M3：兼容映射台账与删除门槛冻结（A5）。
- [ ] M4：最小验证矩阵通过并回填证据（A6）。

## 9. DoD（Wave-A）

- [ ] 已有统一 ProblemSpec 契约，且 mode 维度/参数维度可追踪。
- [ ] AD/ID 接口具备统一 forward/conditions 语义，Dual-safe 规则文档化。
- [ ] 候选治理规则与选择器接口统一，扫描与单点口径一致。
- [ ] 全局配置污染风险被治理（实例化配置或等价隔离）。
- [ ] 迁移映射表与删除门槛已固化。
- [ ] 最小验证矩阵全部通过并可复现。

## 10. 风险与回退

- [ ] 风险：契约重排触发现有调用方边缘兼容问题。
  - 回退：保留 closure adapter 与旧入口兼容壳，逐步迁移。
- [ ] 风险：配置治理改造影响现有导数缓存行为。
  - 回退：提供旧配置路径的临时适配器并加告警。
- [ ] 风险：候选治理统一后数值选择发生跳变。
  - 回退：保留 legacy selector 作为可切换策略并新增回归基线。

## 11. 执行记录

- [x] 2026-03-31：完成草案与 `src/models` 兼容性核验，确认“方向兼容、契约需加硬”。
- [x] 2026-03-31：形成 Wave-A 的 spec+plan（本文档），用于后续实现执行。
- [x] 2026-03-31：完成 A1（ProblemSpec 契约落地）最小闭环：新增 `ProblemSpec`/`build_problem_spec`，补齐 `explicit_residual/explicit_residual!`，并通过定向单测。
- [x] 2026-03-31：完成 A2（AD/ID 契约加固）最小闭环：新增 fixed-mu/flavor-mu 统一 adapter builder，补充 AD 契约 smoke 与架构文档，并通过定向验证。
- [x] 2026-03-31：完成 A3（配置隔离）最小闭环：新增 solver 配置隔离 smoke，`create_implicit_solver` 切换为实例级 config adapters，同时保留兼容 setter 路径。
- [x] 2026-03-31：完成 A4（候选治理接口统一）最小闭环：新增 `CandidateGovernance` 与 context 契约，`evaluate_hard_constraints`/`select_pressure_max_candidate` 支持参数化上下文，补齐单测与回归测试。
- [x] 2026-03-31：完成 A5（兼容层迁移映射）最小闭环：新增迁移台账文档、`solver_migration_map()` 以及 `COMPAT` 非破坏标注，验证统一入口与兼容入口可追溯且行为一致。
