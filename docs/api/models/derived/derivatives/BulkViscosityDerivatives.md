# Bulk Viscosity 组合导数接口

本页说明 derivatives 主题中与 bulk viscosity 公式装配直接相关的公开接口。

## 主要导出

- `bulk_viscosity_coefficients`
- `compute_B_bracket`
- `legacy_transport_c_p`

## `bulk_viscosity_coefficients`

接口形式：

- `bulk_viscosity_coefficients(T_fm, mu_fm; xi=0.0, p_num, t_num, model=nothing)`

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

## `compute_B_bracket`

接口形式：

- `compute_B_bracket(p, M, μq, T, v_n_sq, dμB_dT_sigma, dM_dT, dM_dμB; is_antiquark=false)`

它负责计算体粘滞公式中的 `B` 项，因此应理解为：

- 一个下游公式接口
- 依赖 `bulk_viscosity_coefficients` 或等价上游结果

通常不建议单独从它开始，而是先获取完整系数再进入该接口。

## `legacy_transport_c_p`

这是明确的 compatibility / legacy 导出：

- 它本质上只是读取 `bulk_viscosity_coefficients(...).c_p`
- 应保留在 generated 导出全集中
- 但不应作为 derivatives 主题首页首选入口

文档和导航上必须显式降级它的优先级，避免误导为长期主 API。

## 使用建议

- 需要 bulk viscosity 完整系数时，优先用 `bulk_viscosity_coefficients`
- 需要公式中的单个 `B` 项时，再调用 `compute_B_bracket`
- 仅在兼容旧调用方时，才考虑 `legacy_transport_c_p`

## 测试与集成事实源

当前与该主线直接相关的验证资产包括：

- `tests/integration/relaxtime/test_bulk_viscosity_derivatives.jl`
- `tests/unit/pnjl/test_thermo_derivatives.jl`