# Solver 5-Doc Integration Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate the five active solver design docs into a single, executable, low-risk migration route without changing physics semantics.

**Architecture:** Keep a vector-first numerical kernel shared by solve and derivative paths, move NamedTuple to boundary APIs, and use a schema registry as the only named<->vector mapping source. Adopt spec-first contracts while preserving mode compatibility wrappers and explicit plugin boundaries for legacy/weighted fallback.

**Tech Stack:** Julia, NLsolve, ForwardDiff, ImplicitDifferentiation, existing Models/solver modules, unit+integration+regression tests.

---

## 输入文档评审结论（5 份）

- `docs/dev/active/2026-04-06_solver_target_architecture_blueprint.md`
  - 定位清晰：主链统一、策略中心、兼容层可拔插。
  - 待落地点：需要将 `primary_strategy` 数据结构明确到字段级。
- `docs/dev/active/2026-04-06_solver_constraintspec_interface_draft.md`
  - 主体正确：mode 降格为 spec 输入、`ConstraintDims(x_dim, theta_dim)` 合理。
  - 待统一点：第 5 节 residual 协议仍需明确“内核统一 `theta_vec`，语义层可用 `theta_named`”。
- `docs/dev/active/2026-04-06_solver_derivative_decoupling_interface_and_migration.md`
  - 解耦方向正确：导数层复用 residual 契约，不耦合 solver 迭代细节。
  - 待统一点：与三层文档保持同一 `residual_vec!(..., theta_vec, ...)` 口径。
- `docs/dev/active/2026-04-06_solver_three_layer_minimal_interfaces.md`
  - 可执行性最高：接口/约束/验收口径完整。
  - 建议作为主基线文档，其他文档对齐它。
- `docs/dev/active/2026-04-06_solver_project_sequence_diagram.md`
  - 评审表达好：可作为架构会与 PR 说明图。
  - 待补：异常分支时序（失败、fallback、治理拒绝、兼容插件介入）。

---

## Chunk 1: 文档口径收敛（先统一再开发）

### Task 1: 统一 residual 契约与输入分层描述

**Files:**
- Modify: `docs/dev/active/2026-04-06_solver_constraintspec_interface_draft.md`
- Modify: `docs/dev/active/2026-04-06_solver_derivative_decoupling_interface_and_migration.md`
- Modify: `docs/dev/active/2026-04-06_solver_three_layer_minimal_interfaces.md`

- [x] **Step 1: 写明唯一内核契约**
  - 统一为：`residual_vec!(F, x_vec, theta_vec, cfg)`。
- [x] **Step 2: 写明边界语义契约**
  - 统一为：外层可用 `theta_named`，进入内核前一次映射为 `theta_vec`。
- [x] **Step 3: 消除跨文档冲突描述**
  - 删除或改写 `residual!(..., theta_named, ...)` 作为内核默认表达。
- [x] **Step 4: 自查三文档一致性**
  - 检查术语：`solve_named/solve_vec/derive_named/derive_vec`。
- [ ] **Step 5: Commit**
  - `git commit -m "docs: align solver residual contracts across active design docs"`

### Task 2: 补齐异常分支时序图

**Files:**
- Modify: `docs/dev/active/2026-04-06_solver_project_sequence_diagram.md`

- [x] **Step 1: 增加异常分支图**
  - 包含主方法失败、fallback 尝试、治理拒绝、legacy 插件可选介入。
- [x] **Step 2: 对齐术语**
  - 使用 `primary_strategy`、`Candidate Governance`、`Legacy Adapter`。
- [ ] **Step 3: Commit**
  - `git commit -m "docs: add solver exception-path sequence diagram"`

---

## Chunk 2: 低风险代码骨架（不改物理公式）

### Task 3: 引入 schema 基础设施与通用适配器

**Files:**
- Create: `src/models/solver/SchemaAdapter.jl`
- Modify: `src/models/Models.jl`
- Test: `tests/unit/models/test_solver_schema_adapter.jl`

- [ ] **Step 1: 先写失败测试（映射往返）**
  - 覆盖 `named -> vec -> named` 一致性。
- [ ] **Step 2: 实现 `VarSchema` 与 `SchemaRegistry`**
  - 提供 `register_schema!`、`schema_for`、`validate_schema`。
- [ ] **Step 3: 实现通用转换函数**
  - `named_to_vec`、`vec_to_named`（泛型元素类型，不做 `Float64` 强转）。
- [ ] **Step 4: 跑单测并确认通过**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_schema_adapter.jl"; include("tests/unit/runtests.jl")'`
- [ ] **Step 5: Commit**
  - `git commit -m "feat: add schema registry and named-vector adapters for solver"`

### Task 4: 引入 `primary_strategy` 与最小策略对象

**Files:**
- Create: `src/models/solver/PrimaryStrategy.jl`
- Modify: `src/models/solver/GenericRootEngine.jl`
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/Models.jl`
- Test: `tests/unit/models/test_primary_strategy_contract.jl`

- [ ] **Step 1: 写失败测试（策略字段与默认行为）**
  - method/multi-seed/fallback 合并检查。
- [ ] **Step 2: 定义 `PrimaryStrategy` 数据结构**
  - 默认值对齐现有行为。
- [ ] **Step 3: 在 Solver 入口接入 `primary_strategy`**
  - 不改现有求解语义，仅做参数通道收敛。
- [ ] **Step 4: 跑单测并确认通过**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_primary_strategy_contract.jl,models/test_solver.jl"; include("tests/unit/runtests.jl")'`
- [ ] **Step 5: Commit**
  - `git commit -m "refactor: unify solver method multiseed fallback under primary_strategy"`

---

## Chunk 3: Spec-First 主链切入（兼容保留）

### Task 5: 增加 `solve_vec/solve_named` 双入口骨架

**Files:**
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/solver/ProblemSpec.jl`
- Modify: `src/models/Models.jl`
- Test: `tests/unit/models/test_solver_named_vec_parity.jl`

- [ ] **Step 1: 写失败测试（named 与 vec 输出一致）**
- [ ] **Step 2: 新增 `solve_vec`**
  - 内核入口，仅接受 `theta_vec`。
- [ ] **Step 3: 新增 `solve_named`**
  - 使用 schema 做单次边界转换。
- [ ] **Step 4: 保持旧 `solve(mode,...)` 兼容转发**
  - 旧路径内部转 spec/named/vec 新入口。
- [ ] **Step 5: 跑测试并提交**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_named_vec_parity.jl,models/test_problem_spec_contract.jl"; include("tests/unit/runtests.jl")'`
  - `git commit -m "feat: add named/vec dual solver entry with compatibility wrappers"`

### Task 6: `FixedMu` 向 ProblemSpec 主链对齐（A/B 护栏）

**Files:**
- Modify: `src/models/solver/ProblemSpec.jl`
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/constraint_solver.jl`
- Test: `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`
- Test: `tests/regression/models/test_dimension_agnostic_solver_regression.jl`

- [ ] **Step 1: 先加 A/B 开关与对比测试**
- [ ] **Step 2: 新增 `FixedMu` spec-first forward solve**
- [ ] **Step 3: 在开关下路由到新主链**
- [ ] **Step 4: 跑 integration + regression**
  - Run: `julia --project=. -e 'ENV["INTEGRATION_FILES"]="models/test_solver_auto_backend_semantic_parity.jl"; include("tests/integration/runtests.jl")'`
  - Run: `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_dimension_agnostic_solver_regression.jl"; include("tests/regression/runtests.jl")'`
- [ ] **Step 5: Commit**
  - `git commit -m "refactor: route fixedmu through problemspec chain behind parity guard"`

---

## Chunk 4: 导数层解耦落地

### Task 7: 增加 `derive_vec/derive_named`，复用同一 residual 内核

**Files:**
- Modify: `src/models/implicit_gap.jl`
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/Models.jl`
- Test: `tests/unit/models/test_implicit_gap.jl`
- Test: `tests/unit/models/test_implicit_gap_flavor_mu.jl`
- Test: `tests/integration/models/test_models_implicitdiff_flavor_mu_smoke.jl`

- [ ] **Step 1: 写失败测试（derive_named 与旧接口一致）**
- [ ] **Step 2: 实现 `derive_vec`（仅向量输入）**
- [ ] **Step 3: 实现 `derive_named`（边界转换一次）**
- [ ] **Step 4: 旧 `solve_with_derivatives` 改兼容转发**
- [ ] **Step 5: 跑 unit + integration 并提交**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_implicit_gap.jl,models/test_implicit_gap_flavor_mu.jl"; include("tests/unit/runtests.jl")'`
  - Run: `julia --project=. -e 'ENV["INTEGRATION_FILES"]="models/test_models_implicitdiff_flavor_mu_smoke.jl"; include("tests/integration/runtests.jl")'`
  - `git commit -m "refactor: decouple derivative engine with named/vec dual API"`

---

## Chunk 5: 插件边界与收口

### Task 8: 将 legacy/weighted fallback 明确为插件开关

**Files:**
- Modify: `src/models/solver/WeightedFallback.jl`
- Modify: `src/models/constraint_solver.jl`
- Modify: `src/models/solver/Solver.jl`
- Test: `tests/integration/pnjl/test_trho_scan_semantic_modes_smoke.jl`
- Test: `tests/integration/pnjl/test_trho_scan_solver_backend_models_smoke.jl`

- [ ] **Step 1: 写失败测试（默认主链不触发插件）**
- [ ] **Step 2: 增加显式开关与诊断标记**
- [ ] **Step 3: 保持开启开关时历史行为可用**
- [ ] **Step 4: 跑相关 integration 并提交**
  - Run: `julia --project=. -e 'ENV["INTEGRATION_FILES"]="pnjl/test_trho_scan_semantic_modes_smoke.jl,pnjl/test_trho_scan_solver_backend_models_smoke.jl"; include("tests/integration/runtests.jl")'`
  - `git commit -m "refactor: isolate fallback paths as explicit solver plugins"`

### Task 9: 全量收口验证

**Files:**
- Modify: `docs/api/models/solver/*` (as needed)
- Modify: `docs/dev/active/2026-04-06_*.md` (状态回写)

- [ ] **Step 1: 跑 smoke 回归闭环**
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] **Step 2: 跑治理脚本**
  - `julia --project=. scripts/dev/check_docs_consistency.jl`
  - `julia --project=. scripts/dev/check_legacy_solver_switch_leakage.jl`
- [ ] **Step 3: 更新文档状态与风险清单**
- [ ] **Step 4: Commit**
  - `git commit -m "docs: finalize integrated solver roadmap execution status"`

---

## 执行注意事项

- 遵循 DRY/YAGNI：只做接口收敛与边界清晰化，不重写物理公式。
- 每个 Task 结束即提交，避免大批量混改。
- 若出现性能回退，优先检查是否违反“零重复转换”准则。

## 完成判据（DoD）

- 5 份文档口径一致，术语与契约无冲突。
- `named` 仅在边界层，内核与 AD 共享 `vec` 契约。
- `primary_strategy` 成为唯一策略入口（含 multi-seed）。
- `FixedMu` 能通过 spec-first 主链（至少开关可切并通过 parity）。
- 导数接口完成解耦并保留兼容别名。
- fallback 成为显式插件开关，默认不污染主链。
