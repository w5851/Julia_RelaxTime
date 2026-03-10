# 质量导数接口

本页集中说明 derivatives 主题中与质量导数直接相关的公开接口。实现位于 `src/models/derivatives/ThermoDerivatives.jl`。

## 主要导出

- `mass_derivatives`

## `mass_derivatives`

接口形式：

- `mass_derivatives(T_fm, mu_fm; order=1|2, xi=0.0, p_num, t_num, model=nothing)`

它返回夸克有效质量以及对应导数。当前支持：

- 一阶：`masses`, `dM_dT`, `dM_dmu`
- 二阶：额外返回 `d2M_dT2`, `d2M_dTdmu`, `d2M_dmu2`

## 结果结构

一阶调用时，主题上最应关注的是：

- `masses`
- `dM_dT`
- `dM_dmu`

这些结果会被下游 `thermo_derivatives` 与 `bulk_viscosity_coefficients` 间接消费，因此 `mass_derivatives` 既是独立读出接口，也是下游公式链的一部分。

## 使用建议

- 只需要质量及其导数时，直接使用 `mass_derivatives`
- 需要同步读取压力、能量、熵、数密度及其导数时，应优先转到 [ThermoDerivatives.md](ThermoDerivatives.md)
- 只需常见 bulk viscosity 组合导数时，不必直接依赖 `mass_derivatives` 细节，可改用 `bulk_derivative_coeffs`

## 测试与回归事实源

当前与该接口直接相关的验证资产包括：

- `tests/unit/pnjl/test_thermo_derivatives.jl`
- `tests/integration/models/test_models_derivatives_dual_smoke.jl`
- `tests/baselines/pnjl/baseline_pnjl_thermo_derivatives_fixedpoints_v1.csv`