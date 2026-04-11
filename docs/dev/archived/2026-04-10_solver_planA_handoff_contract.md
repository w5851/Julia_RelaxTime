---
title: Solver Plan-A 交接契约冻结（A4）
archived: true
original: docs/dev/active/2026-04-10_solver_planA_handoff_contract.md
archived_date: 2026-04-11
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Solver Plan-A 交接契约冻结（A4）

## 0) 任务状态（2026-04-11）

- [x] Plan-A A1：拓扑/Facade/行为锁定测试建立并转绿
- [x] Plan-A A2：solver 子域拆分与 include 拓扑迁移完成
- [x] Plan-A A3：R1-R5 边界规则测试与依赖清理完成
- [x] Plan-A A4：交接契约与 API 文档同步完成
- [x] Plan-B B1-B4：声明式 Pipeline、Runner、CLI 迁移与一致性回归完成
- [x] V1：关键门禁与治理检查完成（含 unit/integration/regression 关键验证与 docs/migration guard）

## 1) 冻结 API 签名（Models 边界）

以下签名以 `src/models/solver/api/SolverAPI.jl` + `src/models/Models.jl` 当前导出为准，作为 Plan-A 交接冻结口径：

```julia
# solve_constraint
solve_constraint(model::AbstractQCDModel, mode::FixedMu, T_fm::Real;
    problem_spec::Union{Nothing, ProblemSpec}=nothing,
    μ_fm::Real,
    kwargs...)

solve_constraint(model::AbstractQCDModel, mode::FixedRho, T_fm::Real;
    problem_spec::Union{Nothing, ProblemSpec}=nothing,
    kwargs...)

solve_constraint(model::AbstractQCDModel, mode::FixedEntropy, T_fm::Real; kwargs...)
solve_constraint(model::AbstractQCDModel, mode::FixedSigma, T_fm::Real; kwargs...)
solve_constraint(model::AbstractQCDModel, mode::FixedAsymmetricRho, T_fm::Real; kwargs...)

# solve
solve(model::AbstractPNJLModel, mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
solve(model::AbstractPNJLModel,
    mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma},
    T_fm::Real; kwargs...)
solve(mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
solve(mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, T_fm::Real; kwargs...)

# solve_multi
solve_multi(model::AbstractPNJLModel, mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
solve_multi(model::AbstractPNJLModel, mode::Union{FixedRho, FixedAsymmetricRho}, T_fm::Real; kwargs...)
solve_multi(mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
solve_multi(mode::Union{FixedRho, FixedAsymmetricRho}, T_fm::Real; kwargs...)

# solve_vec
solve_vec(model::AbstractPNJLModel, mode::FixedMu, theta_vec::AbstractVector{<:Real}; kwargs...)
solve_vec(model::AbstractPNJLModel,
    mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma},
    theta_vec::AbstractVector{<:Real}; kwargs...)

# solve_named
solve_named(model::AbstractPNJLModel, mode::FixedMu, theta_named::NamedTuple; kwargs...)
solve_named(model::AbstractPNJLModel,
    mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma},
    theta_named::NamedTuple; kwargs...)
```

补充冻结语义：

- `solve_constraint` 统一走 ProblemSpec 主链；`use_problem_spec` / `allow_legacy_path` / `warn_on_legacy_path` / `fixedmu_use_problem_spec` / `legacy_fallback_plugin` 已移除，传入即 `ArgumentError`。
- `solve_named(FixedMu, ...)` 要求 `:T_fm` + `:μ_fm`；非 `FixedMu` 要求 `:T_fm`。
- `solve_vec(FixedMu, ...)` 约定 `theta_vec` 长度为 2；非 `FixedMu` 约定长度为 1。

冻结锚点（A4）：`HEAD=954adadfd7edb7e0f82bba0c49a435a710f15585`。

## 2) A2 后 include 拓扑摘要

`Models` 入口下，solver include 主干（A2 收口后）为：

1. runtime: `GapSolver.jl` + `ConstraintSolver.jl`
2. contracts/spec: `ConstraintModes.jl` + `ConstraintComponents.jl` + `ProblemSpec.jl`
3. orchestrator: `PrimaryStrategy.jl` + `ProblemSpecOrchestrator.jl` + `SeedStrategies.jl`
4. diagnostics/governance/config: `SolverDiagnostics.jl` + `SolverDiagnosticsTypes.jl` + `SolverRuntimeConfig.jl` + `StateSchema.jl` + `CandidateGovernance.jl`
5. residual/root fallback: `Conditions.jl` + `GenericRootEngine.jl` + `WeightedFallback.jl`
6. API boundary: `SolverAPI.jl`（对外 `solve*` / `solve_constraint` 边界）

说明：A2 后 `ConstraintSolver` 保持“薄聚合 + mode/runtime 实现”定位，`ProblemSpecOrchestrator` 负责主链编排，`SolverAPI` 负责 `Models` 边界合同。

## 3) R1-R5 边界测试检查清单

来源风险定义：`docs/dev/archived/2026-03-31_PNJL求解器解耦框架兼容性评审与Wave-A实现计划.md`。

- [x] R1 双轨职责重叠：已验证 `solve_constraint`/`ProblemSpec` 为唯一主链，旧参数旁路移除。证据：`3a5758a`（Phase 3 contract hardening）、`954adad`（A4 契约冻结）。
- [x] R2 显式参数接口：已验证 residual 组装经 `build_conditions`/`build_residual!` 与契约入口对齐。证据：`766b31d`（boundary rules coverage）、`tests/unit/models/test_solver.jl`。
- [x] R3 全局可变配置：已验证并发路径不依赖旧全局可变状态。证据：`3a5758a`（runtime option/contract 收敛）、`tests/integration/models/test_solver_auto_backend_semantic_parity.jl`。
- [x] R4 维度与类型硬编码：已验证 `solve_vec`/`solve_named` 长度与键名守卫。证据：`3a5758a`（Result/contract 固化）、`tests/unit/models/test_problem_spec_contract.jl`。
- [x] R5 候选治理 context：已验证候选选择走治理契约（hard rules + selector + diagnostics 公共视图）。证据：`032ab81`、`173e8fd`、`tests/unit/models/test_candidate_governance_contract.jl`。

建议执行证据（交接最低集）：

- `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver.jl;models/test_problem_spec_contract.jl;models/test_candidate_governance_contract.jl"; include("tests/unit/runtests.jl")'`
- `julia --project=. -e 'ENV["INTEGRATION_FILES"]="models/test_solver_auto_backend_semantic_parity.jl;models/test_phase_solver_diagnostic_adapter_smoke.jl"; include("tests/integration/runtests.jl")'`
- `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_solver_attempt_engine_convergence_regression.jl;models/test_solver_diagnostic_exception_regression.jl"; include("tests/regression/runtests.jl")'`

## 4) old -> new 映射摘要

完整映射与历史波次细节见附录：

- `docs/dev/archived/2026-03-31_求解器兼容层迁移映射表.md`

本页仅保留交接摘要：

- `solve_fixedmu_constraint` -> `solve_constraint(model, FixedMu(), T_fm; μ_fm=...)`
- `solve_fixedrho_constraint` -> `solve_constraint(model, FixedRho(...), T_fm; ...)`
- `solve_fixedasymrho_constraint` -> `solve_constraint(model, FixedAsymmetricRho(...), T_fm; ...)`
- `solve_fixedentropy_constraint` -> `solve_constraint(model, FixedEntropy(...), T_fm; ...)`
- `solve_fixedsigma_constraint` -> `solve_constraint(model, FixedSigma(...), T_fm; ...)`

交接口径：旧 fixed-* 兼容入口不再作为公开主路径，`Models` 侧统一入口为 `solve*` + `solve_constraint` + `ProblemSpec` 契约。
