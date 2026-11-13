# Bulk Viscosity 组合导数接口

本页说明 derivatives 主题中与 bulk viscosity 公式装配直接相关的公开接口。

## 主要导出

- `bulk_viscosity_coefficients`
- `compute_B_bracket`

## `bulk_viscosity_coefficients`

接口形式：

- `bulk_viscosity_coefficients(T_fm, mu_fm; xi=0.0, p_num, t_num, model=nothing, derivative_backend=:auto, base_state=nothing, base_masses=nothing, base_mu_vec=nothing, thermo_quadrature_kwargs=(;))`

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
- `base_state_source`, `base_state_polished`, `branch_locked`
- `base_residual_norm`
- `primal_solve_count`, `jacobian_factorization_count`, `derivative_series_count`

这组结果面向下游 transport 或 bulk viscosity 公式装配，而不是新用户的最短入门接口。

PNJL 默认后端同样是 TaylorDiff series gap：先从压强 Taylor series 解析得到 `ds/dT`、`ds/dμB`、`dn/dT`、`dn/dμB`，再装配 bulk viscosity 所需组合。旧 `ForwardDiff + ImplicitDifferentiation` 路线已不再作为本接口 fallback。

未提供 `base_state` 时，接口保持独立调用兼容性并自行求一次零阶 gap 解。已经有主平衡态的 workflow 应同时传入 `base_state`、`base_masses` 和 `base_mu_vec`：该路径不会重新进行全局 gap 求解，而是锁定输入分支，最多执行不允许跨支的局部 Newton polish。三个二阶方向 `T`、`mu`、`T+mu` 复用同一个 Jacobian 与一次分解，并同时读出压力和质量 series。

`solve_transport_from_equilibrium` 已采用上述复用路径，并要求 bulk 返回质量与主 equilibrium 质量一致；不匹配时该点的 bulk 计算失败，不会静默拼接两个相分支。

`thermo_quadrature_kwargs::NamedTuple` 接受热力学积分策略与控制参数，格式见
[低温积分示例](Overview.md#低温积分与同支响应)。它同时作用于基态残差、局部 polish、
Jacobian、三个 Taylor 方向及末端压强。已有低温平衡态的调用者应传入生成该态所用的
积分设置。本接口内部构造复用 context，对外使用 `base_state` 传入平衡态。

输入温度要求 `T_fm > 0`。`v_n_sq` 和 `dμB_dT_sigma` 是导数的比值，
在密度与分母同时趋零时可能返回 `NaN`、`Inf` 或对微小误差敏感的有限值。
此时应检查原始 `s`、`n_B` 和二阶压强响应；积分容差收敛本身不能确定未定比值的极限。
`c_p` 在 `abs(n_B) <= sqrt(eps(Float64))` 时按现有约定返回 `NaN`。

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
