# TD 状态导数与 retired implicit compat

本页覆盖 solver 主题中最偏进阶的一组接口：TD-first 状态导数入口、residual problem builder，以及已退休的 legacy implicit solver compat wrapper 边界。

## 主要入口

当前推荐的公开导数入口包括：

- `solve_with_derivatives`
- `solve_pnjl_with_derivatives`
- `solve_pnjl_with_flavor_mu_derivatives`

以下 legacy factory 已从 `Models` 公共 export 面移除，并降级为 retired compat wrapper。它们仍可 qualified 调用，但只会抛出迁移错误，不再构造 `ImplicitFunction`：

- `Models.create_implicit_gap_solver`
- `Models.create_flavor_mu_implicit_gap_solver`
- `Models.create_pnjl_implicit_solver`

## 它们解决什么问题

普通 `solve_gap` / `solve` 负责给出数值解；而推荐导数入口负责进一步回答：

- 状态解 `x(T, μ)` 对参数 `θ` 的导数是多少
- 如何在不把通用嵌套 Dual 路线压到生产热路径上的前提下，使用 TaylorDiff series 获得 `dx/dθ`

因此它们更适合：

- 导数链路
- 涨落量构造
- residual adapter 审计（通过 `build_*_problem` 的 `forward_solve` / `conditions`）

## `solve_*with_derivatives`

这组入口是当前推荐的 PNJL 状态导数包装：

- `solve_pnjl_with_derivatives`
- `solve_pnjl_with_flavor_mu_derivatives`
- `solve_with_derivatives`

它们适合“不想自己管理低层求解对象，但需要直接拿到 `x` 和导数”的场景。

当前 PNJL 便捷导数包装默认使用 TaylorDiff series Newton；`derivative_backend=:auto` 与 `:taylordiff` 等价。旧 `derivative_backend=:forwarddiff` fallback 已下线，会抛出迁移错误。

## Minimal implicit problem/builder 契约

`ImplicitProblem`、`ImplicitSolverConfig`、`build_implicit_solver` 以及
`build_pnjl_fixedmu_problem` / `build_pnjl_flavor_mu_problem` / `build_njl_problem`
仍属于公开的低层构建契约，用来显式拼装 primal solve 与 conditions。

这组接口的定位是维护者/诊断用的 residual adapter surface；它们不是
ThermoDerivatives 或 PNJL 状态导数包装的生产 fallback。需要迁移旧调用时，
优先改到 `solve_pnjl_with_derivatives` / `solve_pnjl_with_flavor_mu_derivatives`。

注意：`build_implicit_solver` 作为兼容导出保留，但已退休；调用会抛出迁移错误。需要 adapter 审计时，直接调用 `problem.forward_solve(θ)` 和 `problem.conditions(θ, x, meta)`。

## `create_implicit_gap_solver`

这是已退休的 legacy 隐函数求解器工厂。它不再围绕模型对象、`gap_residual` 与 `solve_gap` 搭建 `ImplicitFunction`；qualified 调用会抛出迁移错误。

迁移路径：

- NJL/NJL2 residual adapter 审计：`build_njl_problem(model; ...)`
- PNJL fixed-μ residual adapter 审计：`build_pnjl_fixedmu_problem(model; ...)`
- PNJL 状态导数生产入口：`solve_pnjl_with_derivatives(...)`

## `create_pnjl_implicit_solver`

这是面向 PNJL 5 维状态的 retired compat wrapper。旧默认参数曾是：

- `θ = [T, μ]`
- `x = [φu, φd, φs, Φ, Φbar]`

如果你的目标就是 PNJL 平衡态对 `T`、`μ` 的生产导数，应使用 `solve_pnjl_with_derivatives`；如果只是检查 residual adapter，使用 `build_pnjl_fixedmu_problem`。

## `create_flavor_mu_implicit_gap_solver`

这个 retired wrapper 对应旧 flavor 化学势参数化：

- `θ = [T, μ_u, μ_d, μ_s]`

conserved-charge susceptibility 已切到 TaylorDiff / MixedTaylorJet，不再通过该 factory 做低阶对照。需要 flavor-mu 状态导数时使用 `solve_pnjl_with_flavor_mu_derivatives`；需要 residual adapter 审计时使用 `build_pnjl_flavor_mu_problem`。

## 与涨落/导数治理的关系

旧 `FluctuationADPath` 页里最有价值的内容，是对隐式求导边界的整理：

- 近奇异雅可比
- 多分支/相变切换点
- 非光滑算子
- 分母退化

这些边界不是旧架构专属问题，因此本主题直接吸收其语义：当前导数链路以 TD-first wrapper 和 residual problem builders 为边界；真正进入涨落和组合导数时仍需要状态码、奇异点和边界防护。
