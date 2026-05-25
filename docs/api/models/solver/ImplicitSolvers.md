# 隐式求解器工厂与导数接口

本页覆盖 solver 主题中最偏进阶的一组公开表面：implicit solver factory 与隐函数导数接口。

## 主要入口

当前与隐式求导直接相关的公开入口包括：

- `create_implicit_gap_solver`
- `create_flavor_mu_implicit_gap_solver`
- `create_pnjl_implicit_solver`
- `solve_with_derivatives`
- `solve_pnjl_with_derivatives`
- `solve_pnjl_with_flavor_mu_derivatives`

## 它们解决什么问题

普通 `solve_gap` / `solve` 负责给出数值解；而这组接口负责进一步回答：

- 状态解 `x(T, μ)` 对参数 `θ` 的导数是多少
- 如何在不显式对 NLsolve 迭代过程求导的前提下，通过隐函数定理获得 `dx/dθ`

因此它们更适合：

- 导数链路
- 涨落量构造
- 隐式微分验证

## `create_implicit_gap_solver`

这是最通用的隐函数求解器工厂。它围绕模型对象、`gap_residual` 与 `solve_gap` 搭建 `ImplicitFunction`。

可理解为：

- `forward_solve_impl` 负责 primal 数值解
- `conditions_impl` 负责对状态条件函数做隐式微分

## `create_pnjl_implicit_solver`

这是面向 PNJL 5 维状态的更直接工厂，默认参数是：

- `θ = [T, μ]`
- `x = [φu, φd, φs, Φ, Φbar]`

如果你的目标就是 PNJL 平衡态对 `T`、`μ` 的导数，这个入口通常比手工组装更合适。

## `create_flavor_mu_implicit_gap_solver`

这个版本把参数扩展到 flavor 化学势：

- `θ = [T, μ_u, μ_d, μ_s]`

它的意义在于为 flavor-resolved 导数或 conserved-charge 路径提供更细的参数接口。

## `solve_*with_derivatives`

这组入口是在 solver factory 基础上的更高层包装：

- `solve_pnjl_with_derivatives`
- `solve_pnjl_with_flavor_mu_derivatives`
- `solve_with_derivatives`

它们适合“不想自己管理 `ImplicitFunction` 对象，但需要直接拿到 `x` 和导数”的场景。

当前 PNJL 便捷导数包装默认已经切到 TaylorDiff series Newton；`derivative_backend=:forwarddiff` 仍可显式回到旧 `ImplicitFunction + ImplicitDifferentiation` reference/fallback。`create_implicit_gap_solver` / `create_flavor_mu_implicit_gap_solver` 作为兼容工厂仍保留，主要供旧调用点和 fallback 使用。

## 与涨落/导数治理的关系

旧 `FluctuationADPath` 页里最有价值的内容，是对隐式求导边界的整理：

- 近奇异雅可比
- 多分支/相变切换点
- 非光滑算子
- 分母退化

这些边界不是旧架构专属问题，因此本主题直接吸收其语义：implicit solver factory 可作为导数链路基础，但真正进入涨落和组合导数时仍需要状态码、奇异点和边界防护。
