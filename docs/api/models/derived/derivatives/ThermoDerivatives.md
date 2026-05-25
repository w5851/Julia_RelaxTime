# 热力学导数接口

本页集中说明 derivatives 主题中与热力学量导数直接相关的公开接口。实现位于 `src/models/derivatives/ThermoDerivatives.jl`。

当前 PNJL 默认导数后端是 `TaylorDiff + explicit Taylor-series gap Newton`。`ForwardDiff + ImplicitDifferentiation` 仍通过 `derivative_backend=:forwarddiff` 保留为 reference/fallback；默认 `:auto` 对 PNJL 解析为 `:taylordiff`。

## 主要导出

- `thermo_derivatives`
- `bulk_derivative_coeffs`
- `dP_dT`
- `dP_dmu`

## `thermo_derivatives`

接口形式：

- `thermo_derivatives(T_fm, mu_fm; xi=0.0, p_num, t_num, model=nothing, derivative_backend=:auto, linear_solve=:auto, series_residual_tol=1e-7)`

它返回完整的一阶热力学读出，包括：

- `pressure`
- `energy`
- `rho`
- `rho_norm`
- `entropy`
- `dP_dT`, `dP_dmu`
- `dEpsilon_dT`, `dEpsilon_dmu`
- `dn_dT`, `dn_dmu`
- `dP_depsilon_n`, `dP_dn_epsilon`
- 以及质量相关读出 `masses`, `dM_dT`, `dM_dmu`

这是 derivatives 主题最常用的综合入口。

`mu_fm` 在本页接口中表示对称 quark 化学势 `μq`，不是 `μB`。因此 `dP_dmu` 是沿 `μu=μd=μs=μq` 的导数，`rho` 是 `dP_dmu / 3` 对应的 baryon density。

## `bulk_derivative_coeffs`

这是针对体粘滞公式常用组合导数的便捷入口：

- `dP_depsilon_n`
- `dP_dn_epsilon`
- `dM_dT`
- `dM_dmu`

当你只关心组合导数，而不需要完整 `thermo_derivatives` 结果时，优先使用它。

## `dP_dT` / `dP_dmu`

这两个接口提供总压强对温度或化学势的高阶导数：

- `dP_dT(T_fm, mu_fm; order=1, derivative_backend=:auto, ...)`
- `dP_dmu(T_fm, mu_fm; order=1, derivative_backend=:auto, ...)`

它们适合：

- 单独分析压强导数
- 控制导数阶数
- 与中心差分或其它诊断方法做一致性比较

对 PNJL 单方向高阶导数，默认后端一次构造单变量 Taylor series，避免嵌套 Dual 随阶数膨胀。显式 `:forwarddiff` 路径仍可用于低阶对照，但不建议作为高阶热路径。

## 使用建议

- 首次进入 derivatives 主题，优先从 `thermo_derivatives` 开始
- 只需组合导数时，用 `bulk_derivative_coeffs`
- 只需单一标量高阶导数时，用 `dP_dT` 或 `dP_dmu`

## 测试与回归事实源

当前与该主线直接相关的验证资产包括：

- `tests/unit/pnjl/test_thermo_derivatives.jl`
- `tests/regression/pnjl/test_pnjl_thermo_derivatives_regression.jl`
- `tests/baselines/pnjl/baseline_pnjl_thermo_derivatives_fixedpoints_v1.csv`
