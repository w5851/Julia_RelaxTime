---
title: Solver 与导数能力解耦：接口草案与迁移清单
archived: true
original: docs/dev/active/2026-04-06_solver_derivative_decoupling_interface_and_migration.md
archived_date: 2026-04-06
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Solver 与导数能力解耦：接口草案与迁移清单

更新日期：2026-04-06

> 状态说明（2026-04-06 收敛口径）：
> - 本文 `solve_point/derive_point` 为目标态命名。
> - 当前仓库稳定入口为 `solve_vec/solve_named` 与 `derive_vec/derive_named`。
> - 短期不强制新增 `solve_point/derive_point`，避免 API 面重复；若后续统一术语，可先做轻量别名。

## 1. 目标

- 将导数相关能力从 solver 主链中解耦。
- solver 仅负责点求解 `x*(theta)`。
- 导数层仅负责在 `F(x, theta)=0` 上求 `dx/dtheta`。
- 两层共享统一问题契约，不共享实现耦合（当前实现承载名为 `ProblemSpec`）。

## 2. 分层契约

### 2.1 求解层输出

```julia
struct SolvePointResult
    x_star::Vector{Float64}
    converged::Bool
    residual_norm::Float64
    diagnostics
end
```

### 2.2 导数层输出

```julia
struct DerivativeResult
    x_star::Vector{Float64}
    dx_dtheta::Matrix{Float64}
    method::Symbol
    diagnostics
end
```

## 3. 推荐 API

当前实现映射（落地名）：

```julia
solve_vec(model, mode, theta_vec; kwargs...)
solve_named(model, mode, theta_named; kwargs...)

derive_vec(model, theta_vec; kwargs...)
derive_named(model, theta_named; kwargs...)
```

目标态命名（保留为演进方向）：

```julia
solve_point(model, spec::ConstraintSpec, theta; policy=...) -> SolvePointResult

derive_point(model, spec::ConstraintSpec, theta;
             x_star=nothing,
             method::Symbol=:implicit_ad,
             linear_solver=:direct,
             verify_residual::Bool=true) -> DerivativeResult
```

输入约定建议：

- 对外层使用命名参数 `theta_named::NamedTuple`。
- 导数与求解执行层通过适配器将命名输入映射到 `theta` 向量。

语义约束：

- `derive_point` 不依赖 `solve_multi` 的内部候选流程细节。
- 若 `x_star` 未提供，导数层可先调用 `solve_point` 获取平衡点。
- 导数求解只依赖 `residual!` 与线性化对象，不依赖 solver 迭代历史。

## 4. 数学接口最小集合

导数层最小输入：

- `residual_vec!(F, x_vec, theta_vec, cfg)`
- `x_star`
- `theta_vec`（可由边界层 `theta_named` 一次映射得到）
- `x_dim/theta_dim`

隐函数一阶导形式：

- `Jx = dF/dx at (x_star, theta)`
- `Jtheta = dF/dtheta at (x_star, theta)`
- 求解 `Jx * dx_dtheta = -Jtheta`

该路径不需要耦合 NLsolve 内部迭代轨迹。

## 5. 当前仓库映射建议

- `src/models/implicit_gap.jl`：作为 DerivativeEngine 主实现。
- `src/models/solver/Solver.jl`：导数相关入口仅保留兼容转发，逐步下沉。
- `src/models/Models.jl`：当前统一导出 `solve_vec/solve_named` 与 `derive_vec/derive_named`，并保留旧导数 API 别名。

## 6. 函数签名级迁移清单

### 6.1 新增（面向目标态）

状态：暂缓（当前不作为必须落地项；优先保持现有命名收敛）。

文件：`src/models/solver/Solver.jl`

- 新增：
  - `solve_point(model::AbstractQCDModel, spec::ConstraintSpec, theta; kwargs...)`
  - `solve_point(spec::ConstraintSpec, theta; kwargs...)`

文件：`src/models/implicit_gap.jl`

- 新增：
  - `derive_point(model::AbstractQCDModel, spec::ConstraintSpec, theta; x_star=nothing, method::Symbol=:implicit_ad, linear_solver=:direct, verify_residual::Bool=true, kwargs...)`

文件：`src/models/Models.jl`

- 新增导出：
  - `solve_point`
  - `derive_point`

### 6.2 调整（兼容迁移）

文件：`src/models/solver/Solver.jl`

- 现有导数入口改为兼容壳：
  - `solve_with_derivatives(...)` -> 内部转发当前导数主入口（现状为 `solve_pnjl_with_derivatives` / `derive_vec` 路径）

文件：`src/models/implicit_gap.jl`

- 现有 `solve_pnjl_with_derivatives(...)` 标注为兼容别名：
  - 保留签名，内部转发当前导数主入口；未来若引入 `derive_point` 再统一切换。

### 6.3 逐步退役（后续阶段）

- 当上层 workflow 全部迁移后，逐步移除 solver 层导数特例分支与重复包装。

## 7. 验证与准入

最小准入标准：

- `solve_named/solve_vec` 与现有 `solve` 在同输入下 `x_star` 与 `residual_norm` 一致（在容差内）。
- `derive_named/derive_vec` 与现有导数接口在基准点的一阶导一致（在容差内）。
- 兼容 API（`solve_with_derivatives`、`solve_pnjl_with_derivatives`）行为不回归。

建议验证维度：

- `FixedMu` 基准点
- 至少一个非 `FixedMu` 约束点（后续扩展到 spec-first 后）

## 8. 风险与护栏

- 风险：导数层在近奇异 `Jx` 附近稳定性下降。
- 护栏：
  - `verify_residual=true` 默认开启；
  - diagnostics 记录 `cond(Jx)` 或线性求解状态；
  - 提供 `method=:finite_diff` 作为保底策略。

## 9. 结论

- 导数能力与 solver 解耦在当前仓库是可行且推荐的。
- 推荐先以“兼容转发 + 新接口并存”上线，再逐步收口旧导数包装。
