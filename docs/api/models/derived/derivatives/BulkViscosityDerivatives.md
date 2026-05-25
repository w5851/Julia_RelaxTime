# Bulk Viscosity 组合导数接口

本页说明 derivatives 主题中与 bulk viscosity 公式装配直接相关的公开接口。

## 主要导出

- `bulk_viscosity_coefficients`
- `compute_B_bracket`

## `bulk_viscosity_coefficients`

接口形式：

- `bulk_viscosity_coefficients(T_fm, mu_fm; xi=0.0, p_num, t_num, model=nothing, derivative_backend=:auto)`

它返回体粘滞公式所需的一组派生量，包括：

- `v_n_sq`
- `dμB_dT_sigma`
- `masses`
- `dM_dT`
- `dM_dμB`
- `ds_dT`, `ds_dμB`
- `dn_dT`, `dn_dμB`
- `c_p`
- `s`, `n_B`

这组结果面向下游 transport 或 bulk viscosity 公式装配，而不是新用户的最短入门接口。

PNJL 默认后端同样是 TaylorDiff series gap：先从压强 Taylor series 解析得到 `ds/dT`、`ds/dμB`、`dn/dT`、`dn/dμB`，再装配 bulk viscosity 所需组合。旧 `ForwardDiff + ImplicitDifferentiation` 路线已不再作为本接口 fallback。

## `compute_B_bracket`

接口形式：

- `compute_B_bracket(p, M, μq, T, v_n_sq, dμB_dT_sigma, dM_dT, dM_dμB; is_antiquark=false)`

它负责计算体粘滞公式中的 `B` 项，因此应理解为：

- 一个下游公式接口
- 依赖 `bulk_viscosity_coefficients` 或等价上游结果

通常不建议单独从它开始，而是先获取完整系数再进入该接口。

## `c_p` 读取约定

`c_p` 应通过 `bulk_viscosity_coefficients(...).c_p` 读取，不再提供独立兼容别名。

## 使用建议

- 需要 bulk viscosity 完整系数时，优先用 `bulk_viscosity_coefficients`
- 需要公式中的单个 `B` 项时，再调用 `compute_B_bracket`
- 需要 `c_p` 时，直接读取 `bulk_viscosity_coefficients(...).c_p`

## 测试与集成事实源

当前与该主线直接相关的验证资产包括：

- `tests/integration/relaxtime/test_bulk_viscosity_derivatives.jl`
- `tests/unit/pnjl/test_thermo_derivatives.jl`
