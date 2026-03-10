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

- `B/Q/S` 单轴纯方向 `1..4` 阶
- 总二阶组合，包括 `chi11_BQ`、`chi11_BS`、`chi11_QS`

如果超出当前支持范围，接口会抛出参数错误，而不是静默退化。

## `chi_BQS`

这是面向用户的统一别名入口，适合在脚本、扫描和文档里直接表达：

- `chi_BQS(...; orders=(0, 2, 0))`
- `chi_BQS(...; orders=(1, 1, 0))`

当你需要一个统一接口处理纯二阶与 mixed 二阶时，优先使用它。

## `chi_B` 与 `chi1_B` 到 `chi4_B`

这组接口是 baryon 方向首选入口：

- `chi_B(...; order=n)`
- `chi1_B`
- `chi2_B`
- `chi3_B`
- `chi4_B`

其中 `chi1_B` 到 `chi4_B` 是针对常用阶数的便捷封装，适合多数调用场景。

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
- 若你真正关心的是 flavor 化学势梯度或 Hessian，本主题不是首选入口，应等待 `derived/derivatives/` 子主题