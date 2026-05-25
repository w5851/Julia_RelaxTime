# Derivatives 主题总览

本页优先回答“如何在 `Models` 入口下使用热力学导数、质量导数与 bulk viscosity 组合导数”，而不是先展开底层实现细节。

## 何时使用该主题

当你需要以下任一能力时，应从本主题开始：

- 在给定 `(T_fm, mu_fm)` 下获取压力、能量、熵、数密度及其导数
- 计算夸克有效质量对 `T` / `mu` 的一阶或二阶导数
- 读取 `(∂P/∂ε)_n` 与 `(∂P/∂n)_ε` 这类 bulk viscosity 组合导数
- 计算 `bulk_viscosity_coefficients` 或体粘滞公式中的 `B` 项

## 首选公开入口

主题主入口来自 `src/models/Models.jl` 对 `src/models/derivatives/ThermoDerivatives.jl` 的聚合导出。

优先关注：

- `thermo_derivatives`
- `mass_derivatives`
- `bulk_derivative_coeffs`
- `dP_dT`
- `dP_dmu`
- `bulk_viscosity_coefficients`
- `compute_B_bracket`

完整导出基线见 [generated/Exports.md](generated/Exports.md)。

## 最短使用路径

```julia
using .Models

T_fm = 150.0 / Models.Constants_PNJL.ħc_MeV_fm
mu_fm = 320.0 / Models.Constants_PNJL.ħc_MeV_fm

md = Models.mass_derivatives(T_fm, mu_fm; order=1)
td = Models.thermo_derivatives(T_fm, mu_fm)
bulk = Models.bulk_derivative_coeffs(T_fm, mu_fm)
cp = Models.bulk_viscosity_coefficients(T_fm, mu_fm)

d1T = Models.dP_dT(T_fm, mu_fm)
d1mu = Models.dP_dmu(T_fm, mu_fm)
```

## 输入与单位口径

- `T_fm`、`mu_fm` 使用 `fm^-1`
- `order` 当前主要用于 `mass_derivatives` 与 `dP_dT` / `dP_dmu` 的阶数控制
- PNJL derivatives 默认 `derivative_backend=:auto` 会解析为 `:taylordiff`；ThermoDerivatives 主题已下线 `derivative_backend=:forwarddiff`
- 本主题默认采用 `Models` 聚合导出接口，而不是旧 `PNJL.*` 兼容路径

## 推荐使用顺序

多数业务场景建议按以下顺序进入：

1. 先用 `thermo_derivatives` 获取基础热力学量与一阶导数
2. 需要质量信息时再用 `mass_derivatives`
3. 只关心 bulk viscosity 常用组合量时，用 `bulk_derivative_coeffs`
4. 只有在完整体粘滞公式装配时，再进入 `bulk_viscosity_coefficients` 与 `compute_B_bracket`

## 非首页首选入口

以下接口属于下游公式接口，不应作为 derivatives 首页首选入口：

- `compute_B_bracket`
