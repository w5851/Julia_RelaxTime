# 隐式求解器工厂与导数接口

本页覆盖 solver 主题中最偏进阶的一组接口：TD-first 状态导数入口，以及 legacy implicit solver compat factory 的边界。

## 主要入口

当前推荐的公开导数入口包括：

- `solve_with_derivatives`
- `solve_pnjl_with_derivatives`
- `solve_pnjl_with_flavor_mu_derivatives`

以下 factory 已降级为 compat-only，不再从 `Models` 公共 export 面推荐使用，但仍可 qualified 调用：

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
- legacy 隐式微分验证（通过 compat factory）

## `solve_*with_derivatives`

这组入口是当前推荐的 PNJL 状态导数包装：

- `solve_pnjl_with_derivatives`
- `solve_pnjl_with_flavor_mu_derivatives`
- `solve_with_derivatives`

它们适合“不想自己管理低层求解对象，但需要直接拿到 `x` 和导数”的场景。

当前 PNJL 便捷导数包装默认使用 TaylorDiff series Newton；`derivative_backend=:auto` 与 `:taylordiff` 等价。旧 `derivative_backend=:forwarddiff` ImplicitFunction fallback 已下线，会抛出迁移错误。

## Minimal implicit problem/builder 契约

`ImplicitProblem`、`ImplicitSolverConfig`、`build_implicit_solver` 以及
`build_pnjl_fixedmu_problem` / `build_pnjl_flavor_mu_problem` / `build_njl_problem`
仍属于公开的低层构建契约，用来显式拼装 primal solve 与 conditions。

这组接口的定位是维护者/诊断用的 builder surface；它们不是
ThermoDerivatives 或 PNJL 状态导数包装的生产 fallback。需要迁移旧调用时，
优先改到 `solve_pnjl_with_derivatives` / `solve_pnjl_with_flavor_mu_derivatives`；
只有在做 legacy reference 或 adapter 诊断时才直接使用这些 builder。

## `create_implicit_gap_solver`

这是最通用的 legacy 隐函数求解器工厂。它围绕模型对象、`gap_residual` 与 `solve_gap` 搭建 `ImplicitFunction`。

可理解为：

- `forward_solve_impl` 负责 primal 数值解
- `conditions_impl` 负责对状态条件函数做隐式微分

## `create_pnjl_implicit_solver`

这是面向 PNJL 5 维状态的 legacy compat factory，默认参数是：

- `θ = [T, μ]`
- `x = [φu, φd, φs, Φ, Φbar]`

如果你的目标就是 PNJL 平衡态对 `T`、`μ` 的生产导数，应使用 `solve_pnjl_with_derivatives`；该 factory 只保留给显式 legacy reference。

## `create_flavor_mu_implicit_gap_solver`

这个版本把参数扩展到 flavor 化学势：

- `θ = [T, μ_u, μ_d, μ_s]`

它的意义在于为 legacy flavor-resolved reference 或 conserved-charge 低阶对照路径提供更细的参数接口；新生产路线优先使用 TaylorDiff / mixedjet。

## 与涨落/导数治理的关系

旧 `FluctuationADPath` 页里最有价值的内容，是对隐式求导边界的整理：

- 近奇异雅可比
- 多分支/相变切换点
- 非光滑算子
- 分母退化

这些边界不是旧架构专属问题，因此本主题直接吸收其语义：implicit solver factory 可作为导数链路基础，但真正进入涨落和组合导数时仍需要状态码、奇异点和边界防护。
