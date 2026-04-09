# Solver 稳定契约（Result / Diagnostic / Error）

本页定义 Phase 3 固化后的 solver 对外稳定契约，供 `phase` / `scans` / `workflows` 消费。

## 面向用户入口（稳定消费）

- `Models.SolverResult`
- `Models.coerce_solver_result(mode, raw_result; xi_override=nothing)`
- `Models.solver_result_view(result)`
- `Models.solver_contract_version(result)`
- `Models.solver_result_is_success(result; residual_norm_max=1e-6, require_converged=true)`

- `Models.SolverDiagnosticSummary`
- `Models.coerce_solver_diagnostic_summary(diag)`
- `Models.coerce_solver_diagnostic_public_view(diag)`
- `Models.to_public_namedtuple(summary_or_full)`
- `Models.solver_diagnostic_version(summary)`

## Result Contract（`SolverResult`）

当前稳定版本：`contract_version = :v1`（常量 `Models.SOLVER_CONTRACT_VERSION_V1`）。

`SolverResult` 的稳定字段集合（常量 `Models.SOLVER_RESULT_REQUIRED_FIELDS`）为：

- `contract_version::Symbol`
- `mode::ConstraintMode`
- `converged::Bool`
- `solution::Vector{Float64}`
- `x_state::AbstractVector{<:Real}`
- `mu_vec::AbstractVector{<:Real}`
- `omega::Float64`
- `pressure::Float64`
- `rho_norm::Float64`
- `entropy::Float64`
- `energy::Float64`
- `masses::AbstractVector{<:Real}`
- `iterations::Int`
- `residual_norm::Float64`
- `xi::Float64`

### 字段语义与单位

- `x_state`: 平衡态平均场变量向量（无量纲）。
- `mu_vec`: 三味化学势（fm^-1）。
- `omega`: 巨热势密度（fm^-4）。
- `pressure`: 压强（fm^-4）。
- `rho_norm`: 归一化重子密度（`ρ / ρ0`，无量纲）。
- `entropy`: 熵密度（fm^-3）。
- `energy`: 能量密度（fm^-4）。
- `masses`: 有效质量（fm^-1）。
- `residual_norm`: 求解残差范数（无量纲）。
- `xi`: 各向异性参数（无量纲）。

### 兼容策略

- 旧 raw NamedTuple 没有 `contract_version` 时，`coerce_solver_result` 自动补为 `:v1`。
- 不支持未知版本；传入 `contract_version != :v1` 会抛 `ArgumentError`。

## Batch-D 收口（2026-04-09）

### 已移除的过渡字段/参数

- `fixedmu_problem_spec_active`（结果兼容标记）已移除。
- `fixedmu_use_problem_spec`（true/false 开关）已移除。
- `legacy_fallback_plugin`（ProblemSpec 兜底插件开关）已移除。
- `allow_legacy_fallback`（mode solver 过渡参数）已从对外主链与 mode solver 入口移除。

### 迁移窗口期

- 迁移窗口：2026-04-09 至 2026-04-30。
- 在窗口期内，调用上述已移除参数会收到明确 `ArgumentError`，并提示统一使用 ProblemSpec 主链与稳定契约 API。

### 门禁上线日期

- `scripts/dev/check_solver_contract_leakage.jl` 已于 2026-04-09 上线。
- 门禁策略：在 `phase/scans/workflows` 禁止读取 solver 私有字段与已移除兼容标记。

### 替代方式

- 替代 `fixedmu_problem_spec_active`：直接以“调用成功 + 稳定契约字段存在性”判定（`SolverResult.contract_version` / 诊断公共字段）。
- 替代 `fixedmu_use_problem_spec=*`：直接调用 `solve_constraint(...)`，ProblemSpec 为唯一主链。
- 诊断消费统一走 `coerce_solver_diagnostic_public_view` 或 `to_public_namedtuple`。

示例：

```julia
raw = Models.solve_constraint(model, Models.FixedMu(), T_fm;
    μ_fm=μ_fm,
    seed_guess=seed,
    p_num=8,
    t_num=4,
)

res = Models.coerce_solver_result(Models.FixedMu(), raw)
view = Models.solver_result_view(res)

diag_public = haskey(raw, :diagnostic) ?
    Models.coerce_solver_diagnostic_public_view(raw.diagnostic) :
    nothing
```

## Diagnostic Contract（`SolverDiagnosticSummary`）

当前稳定版本：`diagnostic_version = :v1`（常量 `Models.SOLVER_DIAGNOSTIC_VERSION_V1`）。

### 公共稳定字段（允许上层依赖）

常量 `Models.SOLVER_DIAGNOSTIC_PUBLIC_FIELDS`：

- `diagnostic_version`
- `attempt_origin`
- `seed_source`
- `hard_constraint_ok`
- `failed_constraints`
- `error_kind`
- `error_msg`
- `selection_reason`
- `endpoint_cause`
- `continuity_distance`

### 内部调试字段（禁止上层硬依赖）

常量 `Models.SOLVER_DIAGNOSTIC_INTERNAL_FIELDS`：

- `selection_reason_source`

说明：`selection_reason_source` 仍可在 full 诊断或调试中出现，但运行时消费层应使用 `to_public_namedtuple` 或 `coerce_solver_diagnostic_public_view` 获得公共视图。

### summary / full 语义

- `diagnostic_level=:summary`：返回单个 summary（不含 `candidates`）。
- `diagnostic_level=:full`：返回 summary + `candidates` 列表；每个 candidate 可含内部调试字段。

## Error Taxonomy（错误分类）

`error_kind` 的稳定类别：

- `:none`：无错误。
- `:constraint_error`：参数/约束校验失败（`ArgumentError`、`DomainError` 等）。
- `:solver_error`：数值求解执行失败（一般异常）。
- `:interrupt`：求解被中断（`InterruptException`）。

`error_msg` 用于人类可读诊断，不建议上层做关键词分支判断；应优先按 `error_kind` 分类。
