# 守恒荷 Susceptibility 接口

本页集中说明守恒荷 susceptibility 主题的主接口。实现位于 `src/models/derivatives/ConservedChargeSusceptibilities.jl`。

## 主要导出

- `conserved_charge_susceptibility`
- `chi_BQS`
- `chi_B`, `chi1_B`, `chi2_B`, `chi3_B`, `chi4_B`
- `chi_Q`, `chi1_Q`, `chi2_Q`, `chi3_Q`, `chi4_Q`
- `chi_S`, `chi1_S`, `chi2_S`, `chi3_S`, `chi4_S`
- `chi11_BQ`, `chi11_BS`, `chi11_QS`

## `conserved_charge_susceptibility`

这是当前最底层但仍公开的统一入口：

- 输入：`T_fm, muB_fm, muQ_fm, muS_fm`
- 关键字：`orders=(i, j, k)`

当前支持范围：

- `B/Q/S` 单轴纯方向；`derivative_backend=:auto` / `:taylordiff` 时可用 TaylorDiff 走更高阶单方向导数
- mixed BQS 组合；`derivative_backend=:auto` / `:taylordiff` / `:mixedjet` 时使用项目内私有 multivariate Taylor jet，可表达例如 `(2,1,0)`、`(1,1,1)`、`(2,1,1)` 这类高阶 mixed 导数
- 显式 `derivative_backend=:forwarddiff` 只保留 mixed 总二阶 reference，包括 `chi11_BQ`、`chi11_BS`、`chi11_QS`

如果超出当前支持范围，接口会抛出参数错误，而不是静默退化。

`derivative_backend` 支持：

- `:auto`：默认策略。纯 B/Q/S 单方向走单变量 TaylorDiff fast path；mixed BQS 组合走内部 multivariate Taylor jet；`ForwardDiff` 仅作为显式 reference/fallback。
- `:forwarddiff`：保留旧的 ForwardDiff + 隐函数路径，适合作为低阶 reference/fallback。
- `:taylordiff`：纯 B/Q/S 单方向使用单变量 TaylorDiff fast path；mixed BQS 组合使用内部 multivariate Taylor jet。
- `:mixedjet`：显式请求 mixed jet backend。功能上也接受单方向 orders，但单方向会委派回单变量 TaylorDiff fast path，不走通用 `D=1` jet。

## `chi_BQS`

这是面向用户的统一别名入口，适合在脚本、扫描和文档里直接表达：

- `chi_BQS(...; orders=(0, 2, 0))`
- `chi_BQS(...; orders=(1, 1, 0))`
- `chi_BQS(...; orders=(2, 1, 0), derivative_backend=:mixedjet)`

当你需要一个统一接口处理纯方向、mixed 二阶和高阶 mixed BQS 组合时，优先使用它。

## `chi_B` 与 `chi1_B` 到 `chi4_B`

这组接口是 baryon 方向首选入口：

- `chi_B(...; order=n, derivative_backend=:auto)`
- `chi1_B`
- `chi2_B`
- `chi3_B`
- `chi4_B`

其中 `chi1_B` 到 `chi4_B` 是针对常用阶数的便捷封装，适合多数调用场景。

对于高阶单方向 susceptibility，推荐保留默认 `:auto`，或显式写成：

```julia
chi4 = Models.chi_B(T_fm, muB_fm; order=4, derivative_backend=:taylordiff)
```

TaylorDiff 路径还接受生产调优关键字：

- `linear_solve=:auto`：默认按阶数选择线性解策略；当前低阶到 `order=10` 走 `:refactor_each_order`，`order >= 16` 走 `:factorized_each_order`
- `linear_solve=:refactor_each_order | :factorized_each_order | :factorized_batched`：显式选择策略
- `series_iterations=nothing`：默认使用逐阶 Taylor 系数递推；提供非负整数时会额外执行对应次数的 Newton refinement
- `series_residual_tol=1e-7`：检查 primal residual 与按系数尺度归一化后的 series residual，病态点不满足时抛错

这些关键字只影响 TaylorDiff backend；ForwardDiff backend 会沿用旧路径。

对 mixed BQS 高阶导数，这些关键字同样传入内部 mixed jet gap-series Newton。jet 阶数等于 `sum(orders)`，变量维度等于非零 B/Q/S 轴数量；单方向不会走 mixed jet，以保留已验证的单变量 TaylorDiff 性能路径。

## `chi_Q` / `chi_S` 家族

电荷与奇异荷方向接口与 `chi_B` 家族保持同一口径：

- `chi_Q(...; order=n)` 及 `chi1_Q` 到 `chi4_Q`
- `chi_S(...; order=n)` 及 `chi1_S` 到 `chi4_S`

这使得脚本可以在不改变调用风格的前提下切换不同 conserved charge 方向。

## `chi11_BQ` / `chi11_BS` / `chi11_QS`

这三个接口是 mixed 二阶 susceptibility 的首选快捷入口：

- `chi11_BQ`
- `chi11_BS`
- `chi11_QS`

它们本质上都是 `conserved_charge_susceptibility(...; orders=...)` 的主题级快捷别名。

## 使用建议

- 只需要单一 baryon 方向时，优先用 `chi2_B`、`chi3_B`、`chi4_B` 这类明确接口
- 需要在 `B/Q/S` 与 mixed 组合之间切换时，优先用 `chi_BQS`
- 如果你直接调用 `flavor_pressure_derivatives`，它现在也支持 `derivative_backend=:auto|:forwarddiff|:taylordiff|:mixedjet`；其中 `:auto` 与 `:mixedjet` 都会走 TaylorDiff 方向组合，`ForwardDiff` 仍可作为对照参考
- 若你真正关心的是 flavor 化学势梯度或 Hessian，本主题不是首选入口，应等待 `derived/derivatives/` 子主题
