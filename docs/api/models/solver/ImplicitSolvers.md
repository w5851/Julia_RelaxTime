# TD 状态导数与 residual builders

本页覆盖 solver 主题中最偏进阶的一组接口：TD-first 状态导数入口，以及 residual problem builder 的审计边界。

## 主要入口

当前推荐的公开导数入口包括：

- `solve_with_derivatives`
- `solve_pnjl_with_derivatives`
- `solve_pnjl_with_flavor_mu_derivatives`

当前保留的 residual adapter 入口包括：

- `ImplicitProblem`
- `build_njl_problem`
- `build_pnjl_fixedmu_problem`
- `build_pnjl_flavor_mu_problem`

旧 implicit factory 与 retired builder 已移除，不再作为文档入口。

## 它们解决什么问题

普通 `solve_gap` / `solve` 负责给出数值解；推荐导数入口进一步回答：

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

当前 PNJL 便捷导数包装默认使用 TaylorDiff series Newton；`derivative_backend=:auto` 与 `:taylordiff` 等价。旧 `derivative_backend=:forwarddiff` fallback 已下线，会抛出迁移错误。该导数包装只接受 `thermo_backend=:models`，`solver_backend` 只接受 `:models` 或 `:auto`；旧 `:legacy` 后端不会被静默映射，会直接抛出 `ArgumentError`。

## Minimal residual problem 契约

`ImplicitProblem` 以及 `build_pnjl_fixedmu_problem` / `build_pnjl_flavor_mu_problem` / `build_njl_problem` 属于公开的低层 residual adapter 契约，用来显式拼装 primal solve 与 conditions。

这组接口的定位是维护者/诊断用的 adapter surface；它们不是 ThermoDerivatives 或 PNJL 状态导数包装的生产 fallback。需要迁移旧调用时，优先改到 `solve_pnjl_with_derivatives` / `solve_pnjl_with_flavor_mu_derivatives`。

需要 adapter 审计时，直接调用：

```julia
problem = Models.build_pnjl_fixedmu_problem(model; p_num=24, t_num=6)
x, meta = problem.forward_solve([T_fm, mu_fm])
residual = problem.conditions([T_fm, mu_fm], x, meta)
```

## 迁移路径

- NJL/NJL2 residual adapter 审计：`build_njl_problem(model; ...)`
- PNJL fixed-μ residual adapter 审计：`build_pnjl_fixedmu_problem(model; ...)`
- PNJL flavor-μ residual adapter 审计：`build_pnjl_flavor_mu_problem(model; ...)`
- PNJL 状态导数生产入口：`solve_pnjl_with_derivatives(...)`
- PNJL flavor-μ 状态导数生产入口：`solve_pnjl_with_flavor_mu_derivatives(...)`

## 与涨落/导数治理的关系

旧 `FluctuationADPath` 页里最有价值的内容，是对隐式求导边界的整理：

- 近奇异雅可比
- 多分支/相变切换点
- 非光滑算子
- 分母退化

这些边界不是旧架构专属问题，因此本主题直接吸收其语义：当前导数链路以 TD-first wrapper 和 residual problem builders 为边界；真正进入涨落和组合导数时仍需要状态码、奇异点和边界防护。
