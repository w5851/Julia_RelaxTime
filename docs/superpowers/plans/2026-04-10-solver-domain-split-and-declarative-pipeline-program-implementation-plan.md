# Solver Domain Split and Declarative Pipeline Program Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 完成 `solver` 子域拆分与声明式研究流水线落地，并保持 `Models` 统一入口与可复现产物链路。

**Architecture:** 本计划按 Program 执行，拆成 Plan-A（`solver` 结构与边界硬化）和 Plan-B（声明式 Pipeline + CLI 迁移）。Plan-A 冻结 API 与 include 拓扑后，Plan-B 才允许接线，避免重构与编排改动互相干扰。Runner 统一负责 manifest 持久化，stage 仅负责业务产出。

**Tech Stack:** Julia 1.10+ (`src/models/solver/*`, `src/models/workflow/*`, `src/models/entrypoints.jl`, `src/models/Models.jl`, `scripts/pnjl/calculate_phase_structure.jl`), tests (`tests/unit`, `tests/integration`, `tests/regression`), docs (`docs/superpowers/specs`, `docs/api`).

---

## File Structure (Planned)

- Modify: `src/models/Models.jl`
- Create: `src/models/solver/api/SolverAPI.jl`
- Create: `src/models/solver/spec/ProblemSpec.jl`
- Create: `src/models/solver/spec/ConstraintComponents.jl`
- Create: `src/models/solver/spec/ConstraintModes.jl`
- Create: `src/models/solver/spec/Conditions.jl`
- Create: `src/models/solver/orchestrator/ProblemSpecOrchestrator.jl`
- Create: `src/models/solver/orchestrator/PrimaryStrategy.jl`
- Create: `src/models/solver/orchestrator/SeedStrategies.jl`
- Create: `src/models/solver/governance/CandidateGovernance.jl`
- Create: `src/models/solver/governance/WeightedFallback.jl`
- Create: `src/models/solver/runtime/GenericRootEngine.jl`
- Create: `src/models/solver/runtime/GapSolver.jl`
- Create: `src/models/solver/runtime/ConstraintSolver.jl`
- Create: `src/models/solver/runtime/ConstraintSolverCommon.jl`
- Create: `src/models/solver/runtime/ConstraintSolverFixedMu.jl`
- Create: `src/models/solver/runtime/ConstraintSolverFixedRho.jl`
- Create: `src/models/solver/runtime/ConstraintSolverFixedEntropy.jl`
- Create: `src/models/solver/runtime/ConstraintSolverFixedSigma.jl`
- Create: `src/models/solver/runtime/ConstraintSolverFixedAsymmetricRho.jl`
- Create: `src/models/solver/diagnostics/SolverDiagnostics.jl`
- Create: `src/models/solver/diagnostics/SolverDiagnosticsTypes.jl`
- Create: `src/models/solver/diagnostics/ThermoPostprocess.jl`
- Create: `src/models/solver/compat/ImplicitAdapters.jl`
- Create: `src/models/solver/compat/ImplicitGapLegacy.jl`
- Create: `src/models/solver/compat/SchemaAdapter.jl`
- Create: `src/models/solver/config/SolverRuntimeConfig.jl`
- Create: `src/models/solver/config/StateSchema.jl`
- Create: `src/models/workflow/PipelineTypes.jl`
- Create: `src/models/workflow/PipelineRunner.jl`
- Create: `src/models/workflow/StageCatalog.jl`
- Modify: `src/models/entrypoints.jl`
- Modify: `scripts/pnjl/calculate_phase_structure.jl`
- Create: `tests/unit/models/workflow/test_pipeline_types_contract.jl`
- Create: `tests/unit/models/workflow/test_pipeline_runner_behavior.jl`
- Create: `tests/integration/models/test_phase_pipeline_runner_smoke.jl`
- Create: `tests/regression/models/test_phase_pipeline_consistency.jl`
- Modify (if public docs changed): `docs/api/**/*.md`

## Chunk 1: Plan-A (Solver Domain Split + Boundary Hardening)

### Task A1: 先写失败测试锁定 include 拓扑与入口兼容

**Files:**
- Create: `tests/unit/models/solver/test_solver_include_topology.jl`
- Create: `tests/unit/models/solver/test_solver_api_facade_contract.jl`
- Create: `tests/unit/models/solver/test_solver_structure_only_behavior_lock.jl`

- [ ] **Step 1: 写 include 拓扑失败测试（旧路径不再被直接 include）**
- [ ] **Step 2: 写 API facade 失败测试（`solve/solve_multi/solve_constraint/solve_vec/solve_named` 均可从 `Models` 调用）**
- [ ] **Step 3: 写行为锁定失败测试（固定 smoke 输入关键字段快照）**
Contract:
- 锁定字段清单：`converged`, `solution` 长度, `residual_norm`, `omega`, `pressure`, `rho_norm`, `mu_vec`, `masses`
- 比较规则：
  - 布尔/维度类字段（`converged`, `solution` 长度）必须 exact match
  - 数值标量与向量字段使用 `rtol=1e-6`, `atol=1e-8`
  - `NaN` 字段按“同位置均为 `NaN`”判定通过
- 基线来源：使用迁移前主分支 smoke 样例输出（记录在同测试文件头部）
- [ ] **Step 4: 运行单测确认失败**
Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_include_topology.jl"); include("tests/unit/models/solver/test_solver_api_facade_contract.jl"); include("tests/unit/models/solver/test_solver_structure_only_behavior_lock.jl")'`
Expected: FAIL，报缺失新路径或导出契约未满足
- [ ] **Step 5: 提交测试基线**
Run: `git add tests/unit/models/solver/test_solver_include_topology.jl tests/unit/models/solver/test_solver_api_facade_contract.jl tests/unit/models/solver/test_solver_structure_only_behavior_lock.jl && git commit -m "test(models/solver): lock topology facade and structure-only behavior before split"`

### Task A2: 迁移 solver 文件到子域目录并建立 re-export 壳

**Files:**
- Modify: `src/models/Models.jl`
- Move: `src/models/solver/Solver.jl` -> `src/models/solver/api/SolverAPI.jl`
- Move: `src/models/solver/ProblemSpec.jl` -> `src/models/solver/spec/ProblemSpec.jl`
- Move: `src/models/solver/ConstraintComponents.jl` -> `src/models/solver/spec/ConstraintComponents.jl`
- Move: `src/models/solver/ConstraintModes.jl` -> `src/models/solver/spec/ConstraintModes.jl`
- Move: `src/models/solver/Conditions.jl` -> `src/models/solver/spec/Conditions.jl`
- Move: `src/models/solver/ProblemSpecOrchestrator.jl` -> `src/models/solver/orchestrator/ProblemSpecOrchestrator.jl`
- Move: `src/models/solver/PrimaryStrategy.jl` -> `src/models/solver/orchestrator/PrimaryStrategy.jl`
- Move: `src/models/solver/SeedStrategies.jl` -> `src/models/solver/orchestrator/SeedStrategies.jl`
- Move: `src/models/solver/CandidateGovernance.jl` -> `src/models/solver/governance/CandidateGovernance.jl`
- Move: `src/models/solver/WeightedFallback.jl` -> `src/models/solver/governance/WeightedFallback.jl`
- Move: `src/models/solver/GenericRootEngine.jl` -> `src/models/solver/runtime/GenericRootEngine.jl`
- Move: `src/models/solver/GapSolver.jl` -> `src/models/solver/runtime/GapSolver.jl`
- Move: `src/models/solver/ConstraintSolver.jl` -> `src/models/solver/runtime/ConstraintSolver.jl`
- Move: `src/models/solver/ConstraintSolverCommon.jl` -> `src/models/solver/runtime/ConstraintSolverCommon.jl`
- Move: `src/models/solver/ConstraintSolverFixedMu.jl` -> `src/models/solver/runtime/ConstraintSolverFixedMu.jl`
- Move: `src/models/solver/ConstraintSolverFixedRho.jl` -> `src/models/solver/runtime/ConstraintSolverFixedRho.jl`
- Move: `src/models/solver/ConstraintSolverFixedEntropy.jl` -> `src/models/solver/runtime/ConstraintSolverFixedEntropy.jl`
- Move: `src/models/solver/ConstraintSolverFixedSigma.jl` -> `src/models/solver/runtime/ConstraintSolverFixedSigma.jl`
- Move: `src/models/solver/ConstraintSolverFixedAsymmetricRho.jl` -> `src/models/solver/runtime/ConstraintSolverFixedAsymmetricRho.jl`
- Move: `src/models/solver/SolverDiagnostics.jl` -> `src/models/solver/diagnostics/SolverDiagnostics.jl`
- Move: `src/models/solver/SolverDiagnosticsTypes.jl` -> `src/models/solver/diagnostics/SolverDiagnosticsTypes.jl`
- Move: `src/models/solver/ThermoPostprocess.jl` -> `src/models/solver/diagnostics/ThermoPostprocess.jl`
- Move: `src/models/solver/ImplicitAdapters.jl` -> `src/models/solver/compat/ImplicitAdapters.jl`
- Move: `src/models/solver/ImplicitGapLegacy.jl` -> `src/models/solver/compat/ImplicitGapLegacy.jl`
- Move: `src/models/solver/SchemaAdapter.jl` -> `src/models/solver/compat/SchemaAdapter.jl`
- Move: `src/models/solver/SolverRuntimeConfig.jl` -> `src/models/solver/config/SolverRuntimeConfig.jl`
- Move: `src/models/solver/StateSchema.jl` -> `src/models/solver/config/StateSchema.jl`

- [ ] **Step 1: 先迁移 `spec` 子域文件（`ProblemSpec/ConstraintComponents/ConstraintModes/Conditions`）**
- [ ] **Step 2: 再迁移 `orchestrator` 与 `governance` 子域文件（`ProblemSpecOrchestrator/PrimaryStrategy/SeedStrategies/CandidateGovernance/WeightedFallback`）**
- [ ] **Step 3: 迁移 `runtime` 子域文件（保持实现不变）**
- [ ] **Step 4: 迁移 `diagnostics` 子域文件（保持实现不变）**
- [ ] **Step 5: 迁移 `compat` 与 `config` 子域文件（保持实现不变）**
- [ ] **Step 6: 更新 `src/models/Models.jl` include 顺序到新目录**
- [ ] **Step 7: 在 `solver/api/SolverAPI.jl` 统一保留公开入口函数定义与转发**
- [ ] **Step 8: 校验 facade 通过标准（导出含 `solve/solve_multi/solve_constraint/solve_vec/solve_named`，且对外求解入口不再散落在非 `api/` 文件）**
Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_api_facade_contract.jl")'`
Expected: PASS
- [ ] **Step 9: 执行 facade 单一职责校验（`SolverAPI.jl` 仅保留上述 5 个入口及轻量转发；业务 helper 必须下沉到 `spec/orchestrator/runtime/governance/diagnostics/compat/config`）**
Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_api_facade_contract.jl")'`
Expected: PASS，且无“api 层承载业务逻辑”断言失败
- [ ] **Step 10: 在 `docs/dev/active/2026-04-10_solver_planA_handoff_contract.md` 写入 old->new 映射表并标注已迁移状态**
- [ ] **Step 11: 运行 A1 测试并确认转绿**
Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_include_topology.jl"); include("tests/unit/models/solver/test_solver_api_facade_contract.jl")'`
Expected: PASS
- [ ] **Step 12: 运行行为锁定测试（固定 smoke 输入关键字段不漂移）**
Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_structure_only_behavior_lock.jl")'`
Expected: PASS
- [ ] **Step 13: 跑 unit smoke/core 门禁**
Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")' && julia --project=. -e 'ENV["UNIT_PROFILE"]="core"; include("tests/unit/runtests.jl")'`
Expected: PASS
- [ ] **Step 14: 提交结构重排与映射表**
Run: `git add src/models/Models.jl src/models/solver tests/unit/models/solver docs/dev/active/2026-04-10_solver_planA_handoff_contract.md && git commit -m "refactor(models/solver): split solver into domain subfolders with facade compatibility"`

### Task A3: 添加边界违规防回归测试并清理跨层直连

**Files:**
- Create: `tests/unit/models/solver/test_solver_boundary_rules.jl`
- Modify: `src/models/solver/**/*.jl` (仅依赖关系调整)

- [ ] **Step 1: 写失败测试规则-R1（允许 `api -> orchestrator/spec`，禁止反向依赖 `orchestrator/spec -> api`）**
- [ ] **Step 2: 写失败测试规则-R2（禁止 `runtime -> api`）**
- [ ] **Step 3: 写失败测试规则-R3（禁止 `governance -> ConstraintSolver*`）**
- [ ] **Step 4: 写失败测试规则-R4（`compat` 仅允许被 `api/orchestrator` 引用；`spec/governance/runtime/diagnostics/config` 任一引用都应失败）**
- [ ] **Step 5: 写失败测试规则-R5（仅 `orchestrator` 允许同时调用 `runtime` 与 `governance`）**
- [ ] **Step 6: 运行边界测试确认失败**
Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_boundary_rules.jl")'`
Expected: FAIL，提示非法依赖
- [ ] **Step 7: 修复 R1 违规依赖（反向依赖 `orchestrator/spec -> api`）**
- [ ] **Step 8: 修复 R2 违规依赖（`runtime -> api`）**
- [ ] **Step 9: 修复 R3 违规依赖（`governance -> ConstraintSolver*`）**
- [ ] **Step 10: 修复 R4 违规依赖（`compat` 被非 `api/orchestrator` 引用）**
- [ ] **Step 11: 修复 R5 违规依赖（非 `orchestrator` 同时桥接 `runtime+governance`）**
- [ ] **Step 12: 重跑边界测试确认通过**
Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_boundary_rules.jl")'`
Expected: PASS
- [ ] **Step 13: 跑 unit+integration core 门禁**
Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="core"; include("tests/unit/runtests.jl")' && julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'`
Expected: PASS
- [ ] **Step 14: 提交边界硬化**
Run: `git add tests/unit/models/solver/test_solver_boundary_rules.jl src/models/solver && git commit -m "refactor(models/solver): enforce domain boundary dependency rules"`

### Task A4: 冻结 Plan-A 交接契约

**Files:**
- Modify: `docs/dev/active/2026-04-10_solver_planA_handoff_contract.md`
- Modify (if needed): `docs/api/**/*.md`

- [ ] **Step 1: 若交接文档不存在则先创建，最小模板包含冻结 API 签名、include 拓扑、边界测试清单、old->new 映射表**
- [ ] **Step 2: 写明 Plan-A 冻结项（facade 签名、include 拓扑、边界测试清单）**
- [ ] **Step 3: 若 `Models` 稳定入口变化，同步 `docs/api`**
- [ ] **Step 4: 跑文档治理检查**
Run: `julia --project=. scripts/dev/check_docs_consistency.jl && julia --project=. scripts/dev/check_active_docs_governance.jl`
Expected: PASS
- [ ] **Step 5: 提交 Plan-A 冻结契约文档**
Run: `git add docs/dev/active/2026-04-10_solver_planA_handoff_contract.md docs/api && git commit -m "docs(dev): record Plan-A handoff contract for pipeline migration"`

## Chunk 2: Plan-B (Declarative Pipeline + CLI Migration)

### Task B1: 先写失败测试锁定 Pipeline 类型契约

**Files:**
- Create: `tests/unit/models/workflow/test_pipeline_types_contract.jl`
- Create: `src/models/workflow/PipelineTypes.jl`

- [ ] **Step 1: 写失败测试覆盖 `PipelineSpec/PipelineStage/PipelineContext/StageResult/PipelineArtifact` 字段与类型**
- [ ] **Step 2: 写失败测试覆盖 Symbol<->String 规范化规则**
- [ ] **Step 3: 运行测试确认失败**
Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_pipeline_types_contract.jl")'`
Expected: FAIL，提示缺失类型或字段
- [ ] **Step 4: 实现最小类型定义与规范化函数**
- [ ] **Step 5: 重跑测试确认通过**
Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_pipeline_types_contract.jl")'`
Expected: PASS
- [ ] **Step 6: 提交类型契约**
Run: `git add src/models/workflow/PipelineTypes.jl tests/unit/models/workflow/test_pipeline_types_contract.jl && git commit -m "feat(models/workflow): add pipeline core contract types and normalization"`

### Task B2: 实现 PipelineRunner（依赖校验、fail-fast、manifest 持久化）

**Files:**
- Create: `src/models/workflow/PipelineRunner.jl`
- Create: `tests/unit/models/workflow/test_pipeline_runner_behavior.jl`
- Create: `tests/unit/models/solver/test_solver_structure_only_behavior_lock.jl`

- [ ] **Step 1: 写失败测试覆盖依赖校验（缺依赖、循环依赖、重复 provide、重复 stage id）**
- [ ] **Step 2: 写失败测试覆盖失败路径（failed+skipped stage、error_kind/error_msg、required_outputs 仅成功路径校验）**
- [ ] **Step 3: 写失败测试覆盖 Runner 统一落盘 manifest（成功/失败都写）**
- [ ] **Step 4: 运行测试确认失败**
Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_pipeline_runner_behavior.jl")'`
Expected: FAIL
- [ ] **Step 5: 实现 Runner 最小功能并让测试转绿**
- [ ] **Step 6: 重跑测试确认通过**
Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_pipeline_runner_behavior.jl")'`
Expected: PASS
- [ ] **Step 7: 提交 Runner**
Run: `git add src/models/workflow/PipelineRunner.jl tests/unit/models/workflow/test_pipeline_runner_behavior.jl && git commit -m "feat(models/workflow): add fail-fast pipeline runner with manifest persistence"`

### Task B3: 实现 StageCatalog 并接入 `Models.run_phase_pipeline`

**Files:**
- Create: `src/models/workflow/StageCatalog.jl`
- Modify: `src/models/entrypoints.jl`
- Modify: `src/models/Models.jl`
- Create: `tests/integration/models/test_phase_pipeline_runner_smoke.jl`

- [ ] **Step 1: 写失败集成测试覆盖 7 个标准 stage 串行执行**
- [ ] **Step 2: 写失败测试覆盖 `run_phase_pipeline` 薄封装转调 Runner**
- [ ] **Step 3: 运行集成测试确认失败**
Run: `julia --project=. -e 'include("tests/integration/models/test_phase_pipeline_runner_smoke.jl")'`
Expected: FAIL
- [ ] **Step 4: 实现 StageCatalog 与 entrypoint 接线**
- [ ] **Step 5: 重跑集成测试确认通过**
Run: `julia --project=. -e 'include("tests/integration/models/test_phase_pipeline_runner_smoke.jl")'`
Expected: PASS
- [ ] **Step 6: 提交阶段编排接线**
Run: `git add src/models/workflow src/models/entrypoints.jl src/models/Models.jl tests/integration/models/test_phase_pipeline_runner_smoke.jl && git commit -m "refactor(models): route phase pipeline through declarative runner"`

### Task B4: 迁移 CLI 到 PipelineSpec 并补齐回归一致性测试

**Files:**
- Modify: `scripts/pnjl/calculate_phase_structure.jl`
- Create: `tests/regression/models/test_phase_pipeline_consistency.jl`

- [ ] **Step 1: 写失败回归测试覆盖旧/新 pipeline 关键字段一致性（rtol=1e-6, atol=1e-8）**
- [ ] **Step 2: 写失败测试覆盖 manifest_v1 必填字段、UTC 时间格式、SHA-256 hash 字段存在**
- [ ] **Step 3: 运行回归测试确认失败**
Run: `julia --project=. -e 'include("tests/regression/models/test_phase_pipeline_consistency.jl")'`
Expected: FAIL
- [ ] **Step 4: 改造 CLI：参数解析 -> 构造 `PipelineSpec` -> Runner 执行**
- [ ] **Step 5: 重跑回归测试确认通过**
Run: `julia --project=. -e 'include("tests/regression/models/test_phase_pipeline_consistency.jl")'`
Expected: PASS
- [ ] **Step 6: 提交 CLI 迁移与一致性回归**
Run: `git add scripts/pnjl/calculate_phase_structure.jl tests/regression/models/test_phase_pipeline_consistency.jl && git commit -m "refactor(scripts/pnjl): migrate phase CLI to PipelineSpec orchestration"`

## Chunk 3: Program Verification and Closure

### Task V1: 执行全量门禁矩阵并归档证据

**Files:**
- Modify: `docs/dev/active/2026-04-10_solver_planA_handoff_contract.md`
- Modify: `docs/api/**/*.md` (if needed)

- [ ] **Step 1: 跑 unit smoke/core**
Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")' && julia --project=. -e 'ENV["UNIT_PROFILE"]="core"; include("tests/unit/runtests.jl")'`
Expected: PASS
- [ ] **Step 2: 跑 integration smoke/core**
Run: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")' && julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'`
Expected: PASS
- [ ] **Step 3: 跑 regression smoke/core**
Run: `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")' && julia --project=. -e 'ENV["REGRESSION_PROFILE"]="core"; include("tests/regression/runtests.jl")'`
Expected: PASS
- [ ] **Step 4: 跑脚本 smoke（PipelineSpec 路径）**
Run: `julia --project=. scripts/pnjl/calculate_phase_structure.jl --preset=smoke`
Expected: PASS，生成 manifest 与 artifact 元数据
- [ ] **Step 5: 跑治理检查**
Run: `julia --project=. scripts/dev/check_unit_skip_policy.jl && julia --project=. scripts/dev/check_docs_consistency.jl && julia --project=. scripts/dev/check_active_docs_governance.jl && julia --project=. scripts/dev/check_pnjl_migration_guard.jl`
Expected: PASS
- [ ] **Step 6: 提交收尾**
Run: `git add docs/dev/active/2026-04-10_solver_planA_handoff_contract.md docs/api tests src/models scripts/pnjl && git commit -m "refactor(models/workflow): close solver split and declarative pipeline program gates"`
